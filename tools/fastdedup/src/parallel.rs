use std::hash::{BuildHasher, Hash};
use std::io::Read;

use anyhow::{Context, Result, bail};
use crossbeam_channel::{Receiver, Sender, bounded};
use needletail::parser::LineEnding;
use needletail::{parse_fastx_file, parse_fastx_reader};
use rustc_hash::{FxBuildHasher, FxHashSet};

use crate::hasher::SequenceHasher;
use crate::utils::{base_read_id, compute_write_buffer_size, validate_format_compatibility};
use crate::writer::{
    ByteCounter, OutputFormat, StreamCompression, prepare_writer, write_fasta_record,
};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

use super::DedupConfig;

// ==========================================
// TYPES
// ==========================================

/// A thread-safe hash verifier using hash-space sharding to avoid lock contention.
pub(crate) struct ShardedHashVerifier<T: Hash + Eq> {
    shards: Vec<FxHashSet<T>>,
    num_shards: usize,
}

impl<T: Hash + Eq> ShardedHashVerifier<T> {
    pub fn new(capacity: usize, num_shards: usize) -> Result<Self> {
        let shard_capacity = (capacity / num_shards).max(1);
        let mut shards = Vec::with_capacity(num_shards);
        for _ in 0..num_shards {
            let mut shard = FxHashSet::default();
            shard
                .try_reserve(shard_capacity)
                .map_err(|_| anyhow::anyhow!("OOM: sharded verifier shard allocation failed (capacity: {shard_capacity})"))?;
            shards.push(shard);
        }
        Ok(Self { shards, num_shards })
    }

    /// Returns `true` if the hash is newly inserted (unique), `false` if duplicate.
    pub fn verify(&mut self, hash: T) -> bool {
        let idx = self.shard_index(&hash);
        self.shards[idx].insert(hash)
    }

    fn shard_index(&self, hash: &T) -> usize {
        (FxBuildHasher.hash_one(hash) as usize) % self.num_shards
    }
}

/// Owned single-end record data that can be sent across threads.
///
/// FD-9: `crlf` carries the record's own line ending. The serial path writes
/// through needletail, which re-emits the line ending it parsed, so a CRLF
/// library round-tripped as CRLF at `-j 1` and as LF at `-j 4` — an
/// untested output difference between the two engines.
#[derive(Clone, Debug)]
pub(crate) struct RawRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub crlf: bool,
}

/// The line ending a record was parsed with.
#[inline]
fn ending(crlf: bool) -> &'static [u8] {
    if crlf { b"\r\n" } else { b"\n" }
}

/// Owned paired-end record data.
#[derive(Clone, Debug)]
pub(crate) struct RawPair {
    pub r1: RawRecord,
    pub r2: RawRecord,
}

/// Hashed single-end record, sent from worker to writer.
pub(crate) struct HashedRecord<T> {
    pub hash: T,
    pub record: RawRecord,
}

/// Hashed pair, sent from worker to writer.
pub(crate) struct HashedPair<T> {
    pub hash: T,
    pub pair: RawPair,
}

const CHUNK_SIZE: usize = 4096;

// ==========================================
// CHANNEL-BASED PIPE for gzippy streaming
// ==========================================

/// Writer side: implements `Write` by shipping whole chunks through a crossbeam channel.
///
/// FD-8: the bytes are appended into one reusable buffer and the buffer is
/// MOVED to the compression thread when it fills. The previous form allocated a
/// `Vec` per `write` call (`buf.to_vec()`) and the reader copied it again, so
/// every output byte was memcpy'd twice plus one allocation per chunk.
struct ChanWriter {
    tx: Sender<Vec<u8>>,
    buf: Vec<u8>,
}

/// Chunk size handed to the compression thread. Large enough that the channel
/// never sees a fragment, small enough that one buffer is ~256 KiB.
const CHUNK_BYTES: usize = 256 * 1024;

impl ChanWriter {
    fn new(tx: Sender<Vec<u8>>) -> Self {
        Self {
            tx,
            buf: Vec::with_capacity(CHUNK_BYTES),
        }
    }

    /// Move the buffered bytes to the compression thread.
    fn ship(&mut self) -> std::io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let filled = std::mem::replace(&mut self.buf, Vec::with_capacity(CHUNK_BYTES));
        self.tx.send(filled).map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::BrokenPipe, "gzippy compression thread died")
        })
    }
}

impl Write for ChanWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.buf.extend_from_slice(buf);
        if self.buf.len() >= CHUNK_BYTES {
            self.ship()?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.ship()
    }
}

impl Drop for ChanWriter {
    /// Close the pipe with the tail. `BufWriter` above this flushes on drop, but
    /// a caller that only drops the writer must still not lose the last chunk.
    fn drop(&mut self) {
        let _ = self.ship();
    }
}

/// Reader side: implements `Read` by receiving chunks from a crossbeam channel.
struct ChanReader {
    rx: Receiver<Vec<u8>>,
    current: Vec<u8>,
    pos: usize,
}

impl Read for ChanReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.pos >= self.current.len() {
            match self.rx.recv() {
                Ok(chunk) => {
                    self.current = chunk;
                    self.pos = 0;
                }
                Err(_) => return Ok(0), // EOF — channel closed
            }
        }
        let avail = &self.current[self.pos..];
        let n = buf.len().min(avail.len());
        buf[..n].copy_from_slice(&avail[..n]);
        self.pos += n;
        Ok(n)
    }
}

// ==========================================
// HELPERS
// ==========================================

/// Open a FASTX input as a byte stream.
///
/// Compression is left to needletail, which sniffs the magic bytes and wraps
/// the reader in `MultiGzDecoder`. A `.gz` suffix is not needed to decompress,
/// and the hand-rolled `flate2::read::GzDecoder` that used to sit here stopped
/// at the END OF THE FIRST GZIP MEMBER: BGZF (what `bgzip`/`samtools fastq`
/// emit) is multi-member, so `-j 2+` read the first member only and exited 0
/// with a truncated report.
fn open_input_reader(path: &str) -> Result<Box<dyn Read + Send + 'static>> {
    if crate::stream::is_stream(path) {
        return Ok(Box::new(std::io::stdin()));
    }
    let file = File::open(path).with_context(|| format!("Failed to open input file: {path}"))?;
    Ok(Box::new(file))
}

/// Wait on an optional gzippy compression thread handle.
fn wait_for_gzip_thread(
    handle_opt: Option<std::thread::JoinHandle<Result<()>>>,
    label: &str,
) -> Result<()> {
    if let Some(handle) = handle_opt {
        match handle.join() {
            Ok(Ok(())) => Ok(()),
            Ok(Err(e)) => Err(e.context(format!("{label} failed"))),
            Err(_) => bail!("{label} thread panicked"),
        }
    } else {
        Ok(())
    }
}

/// Prepare an output writer, optionally using gzippy for parallel gzip compression
/// when `num_threads > 1`, output is `.gz`, and `compression_level >= 6` (levels
/// 6+ with threads>1 produce standard gzip; levels <6 with threads>1 produce
/// gzippy's proprietary "GZ" multi-block format which is not compatible with
/// standard tools).
///
/// Falls back to `prepare_writer` (flate2) otherwise.
///
/// Returns (writer, output_format, optional_gzippy_thread_handle).
///
/// The gzippy thread handles compressing the stream in the background.
/// Caller MUST drop the returned writer BEFORE joining the handle.
fn create_output_writer(
    path: &str,
    compression_level: u32,
    num_threads: usize,
    out_format: OutputFormat,
    out_compression: StreamCompression,
    force: bool,
) -> Result<(
    Box<dyn Write>,
    OutputFormat,
    Option<std::thread::JoinHandle<Result<()>>>,
)> {
    if crate::stream::is_stream(path) {
        // D-4: stdout cannot be seeked or reopened, so there is no resume — but
        // there is no reason to give up parallel compression: gzippy writes to
        // any `Write + Send` and never seeks, so the stream path takes the same
        // channel pipe as the file path (FD-2). The level >= 6 guard matches the
        // file path: below it, threads produce gzippy's non-standard multi-block
        // format.
        if out_compression == StreamCompression::Gz
            && num_threads > 1
            && compression_level >= 6
        {
            let (tx, rx) = bounded::<Vec<u8>>(16);
            let writer: Box<dyn Write> = Box::new(ChanWriter::new(tx));
            let level = compression_level.min(12) as u8;
            let handle = std::thread::spawn(move || -> Result<()> {
                let reader = ChanReader {
                    rx,
                    current: Vec::new(),
                    pos: 0,
                };
                let out = std::io::BufWriter::with_capacity(1 << 20, std::io::stdout());
                gzippy::compress_to_writer_with_threads(reader, out, level, num_threads)
                    .context("gzippy compression failed for the output stream")?;
                Ok(())
            });
            return Ok((writer, out_format, Some(handle)));
        }
        let writer = crate::writer::prepare_stream_writer(out_compression, compression_level);
        return Ok((writer, out_format, None));
    }
    let output_path = Path::new(path);
    let output_format = OutputFormat::from_extension(output_path);
    let is_gz = OutputFormat::is_gz(output_path);

    // An existing `.gz` output cannot be resumed (multi-member gzip), and the
    // gzippy path truncates unconditionally: fall through to `prepare_writer`
    // so it refuses with the serial path's message instead of overwriting.
    if is_gz && num_threads > 1 && compression_level >= 6 && (force || !output_path.exists()) {
        // Gzippy parallel compression path — channel-based streaming pipe
        // Levels 6+ with threads>1 use PipelinedGzEncoder → standard gzip
        let (tx, rx) = bounded::<Vec<u8>>(16);
        let writer: Box<dyn Write> = Box::new(ChanWriter::new(tx));

        let out_path = path.to_string();
        let level = compression_level.min(12) as u8;
        let handle = std::thread::spawn(move || -> Result<()> {
            let reader = ChanReader {
                rx,
                current: Vec::new(),
                pos: 0,
            };
            let output_file = File::create(&out_path)
                .with_context(|| format!("Failed to create output file: {out_path}"))?;
            gzippy::compress_to_writer_with_threads(reader, output_file, level, num_threads)
                .with_context(|| format!("gzippy compression failed for {out_path}"))?;
            Ok(())
        });

        Ok((writer, output_format, Some(handle)))
    } else {
        // Standard path: flate2 GzEncoder (standard gzip at any level/threads)
        let (writer, format) = prepare_writer(output_path, force, compression_level)?;
        Ok((writer, format, None))
    }
}

// ==========================================
// RESUME (append / truncate-to-last-valid)
// ==========================================

/// Cut `path` back to the last complete record, as the serial path does after
/// a crash left a torn tail behind.
fn truncate_to_valid_len(path: &Path, valid_bytes: u64, label: &str, verbose: bool) -> Result<()> {
    let current_size = fs::metadata(path)?.len();
    if valid_bytes < current_size {
        if verbose {
            println!("Truncating {label} from {current_size} to {valid_bytes} bytes.");
        }
        let file = OpenOptions::new().write(true).open(path)?;
        file.set_len(valid_bytes)?;
    }
    Ok(())
}

/// Resume an existing uncompressed output when `--force` is absent: preload its
/// hashes into the verifier (so records already on disk are never appended
/// twice) and truncate a torn trailing record back to the last valid byte.
///
/// Mirrors the serial path's preload/truncate sequence, which the threaded
/// writer used to skip by forcing truncation. Nothing to do for a stream
/// (always written fresh) or a `.gz` output (`prepare_writer` refuses those).
fn resume_existing_output<T: SequenceHasher>(
    path: &str,
    cfg: &DedupConfig,
    verifier: &mut ShardedHashVerifier<T>,
) -> Result<()> {
    let bound = Path::new(path);
    if cfg.force || crate::stream::is_stream(path) || !bound.exists() || OutputFormat::is_gz(bound)
    {
        return Ok(());
    }

    if cfg.verbose {
        println!("Preloading sequences from existing output...");
    }

    let mut reader = parse_fastx_file(path).context("Error opening preload file")?;
    let mut valid_bytes = ByteCounter(0);
    let mut preloads = 0;

    while let Some(record_result) = reader.next() {
        let record = match record_result {
            Ok(record) => record,
            Err(e) => {
                if cfg.verbose {
                    eprintln!("Incomplete sequence detected at end of file ({e}).");
                    eprintln!("Calculating fail-safe truncation point...");
                }
                break;
            }
        };
        verifier.verify(T::hash_sequence(&record.seq()));
        preloads += 1;
        let _ = record.write(&mut valid_bytes, None);
    }

    if preloads > 0 && cfg.verbose {
        println!("{preloads} sequences preloaded.");
    }

    truncate_to_valid_len(bound, valid_bytes.0, "file", cfg.verbose)
}

/// Paired-end twin of [`resume_existing_output`]: preloads the pair hashes and
/// truncates both outputs so R1 and R2 stay synchronized.
fn resume_existing_paired_outputs<T: SequenceHasher>(
    path_r1: &str,
    path_r2: &str,
    cfg: &DedupConfig,
    verifier: &mut ShardedHashVerifier<T>,
) -> Result<()> {
    let bound_r1 = Path::new(path_r1);
    let bound_r2 = Path::new(path_r2);
    if cfg.force
        || crate::stream::is_stream(path_r1)
        || crate::stream::is_stream(path_r2)
        || !bound_r1.exists()
        || !bound_r2.exists()
        || OutputFormat::is_gz(bound_r1)
        || OutputFormat::is_gz(bound_r2)
    {
        return Ok(());
    }

    if cfg.verbose {
        println!("Preloading and synchronizing pairs from existing outputs...");
    }

    let mut reader_r1 = parse_fastx_file(path_r1).context("Error opening R1 preload file")?;
    let mut reader_r2 = parse_fastx_file(path_r2).context("Error opening R2 preload file")?;

    let mut valid_bytes_r1 = ByteCounter(0);
    let mut valid_bytes_r2 = ByteCounter(0);
    let mut preloads = 0;

    while let (Some(r1_result), Some(r2_result)) = (reader_r1.next(), reader_r2.next()) {
        let (record_r1, record_r2) = match (r1_result, r2_result) {
            (Ok(r1), Ok(r2)) => (r1, r2),
            _ => {
                if cfg.verbose {
                    eprintln!("Incomplete sequence or desynchronization detected at end of file.");
                    eprintln!("Calculating fail-safe truncation points for R1 and R2...");
                }
                break;
            }
        };
        verifier.verify(T::hash_pair(&record_r1.seq(), &record_r2.seq()));
        preloads += 1;
        let _ = record_r1.write(&mut valid_bytes_r1, None);
        let _ = record_r2.write(&mut valid_bytes_r2, None);
    }

    if preloads > 0 && cfg.verbose {
        println!("{preloads} pairs preloaded and synchronized.");
    }

    truncate_to_valid_len(bound_r1, valid_bytes_r1.0, "R1", cfg.verbose)?;
    truncate_to_valid_len(bound_r2, valid_bytes_r2.0, "R2", cfg.verbose)?;

    Ok(())
}

// ==========================================
// SINGLE-END PARALLEL
// ==========================================

/// Executes single-end deduplication using a multi-threaded pipeline.
///
/// One reader thread streams records (decompressing `.gz` via flate2),
/// `num_workers` worker threads hash, and the calling thread writes output.
/// When the output is `.gz` and `compression_level >= 6`, output is
/// compressed in parallel via gzippy; otherwise flate2 single-threaded.
pub(crate) fn execute_deduplication_parallel<T: SequenceHasher + Send + 'static>(
    input_path: &str,
    output_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    // --dry-run writes nothing: the serial counting pass reaches the identical
    // verdict (same refusals, same duplicate count) without creating a file.
    if cfg.dry_run {
        return crate::processor::execute_deduplication::<T>(input_path, output_path, cfg);
    }

    let num_workers = resolve_thread_count(cfg.threads);
    // FD-3: `-j N` sizes the worker pool, not the thread count. A run is
    // 1 reader + (N-1) workers + this writer thread + (N-2) gzippy compression
    // threads, i.e. 2N-1 threads. Workers leave one core to the reader/commit
    // thread and compression takes the rest, so deflate never competes with the
    // table insert for every core; the budget is a split, not a cap.
    let worker_threads = num_workers.saturating_sub(1).max(1);

    // Resume an existing output before any thread starts, so the preload pass
    // cannot race the unbounded worker channel.
    let mut verifier = ShardedHashVerifier::<T>::new(cfg.estimated_capacity, worker_threads)?;
    resume_existing_output(output_path, cfg, &mut verifier)?;

    let (chunk_tx, chunk_rx): (Sender<Vec<RawRecord>>, Receiver<Vec<RawRecord>>) =
        bounded(num_workers * 2);
    let (hash_tx, hash_rx): (Sender<HashedRecord<T>>, Receiver<HashedRecord<T>>) =
        crossbeam_channel::unbounded();

    // The reader reports, off the first record, whether the input carries
    // quality scores (FASTQ) or not (FASTA). The writer validates the
    // input/output format pairing before writing any record, so the threaded
    // verdict is identical to the serial path's.
    let (fmt_tx, fmt_rx): (Sender<bool>, Receiver<bool>) = bounded(1);

    // --- READER THREAD ---
    let reader_handle = std::thread::spawn({
        let chunk_tx = chunk_tx.clone();
        let input_owned = input_path.to_string();
        move || -> Result<usize> {
            let reader = open_input_reader(&input_owned)?;
            let mut fastx =
                parse_fastx_reader(reader).map_err(|e| anyhow::anyhow!("Parse error: {e}"))?;

            let mut total = 0;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);
            let mut format_reported = false;

            while let Some(record) = fastx.next() {
                let rec = record.context("Invalid sequence")?;
                if !format_reported {
                    let has_quality = rec.qual().is_some_and(|q| !q.is_empty());
                    let _ = fmt_tx.send(has_quality);
                    format_reported = true;
                }
                let id = rec.id().to_vec();
                let seq = rec.seq().to_vec();
                let qual = rec.qual().map(|q| q.to_vec()).unwrap_or_default();
                let crlf = matches!(rec.line_ending(), LineEnding::Windows);

                chunk.push(RawRecord { id, seq, qual, crlf });
                total += 1;

                if chunk.len() >= CHUNK_SIZE {
                    if chunk_tx.send(chunk).is_err() {
                        return Ok(total);
                    }
                    chunk = Vec::with_capacity(CHUNK_SIZE);
                }
            }
            if !chunk.is_empty() {
                let _ = chunk_tx.send(chunk);
            }

            Ok(total)
        }
    });

    // --- WORKER THREADS ---
    let mut worker_handles = Vec::with_capacity(worker_threads);
    for _ in 0..worker_threads {
        let chunk_rx = chunk_rx.clone();
        let hash_tx = hash_tx.clone();
        worker_handles.push(std::thread::spawn(move || {
            while let Ok(chunk) = chunk_rx.recv() {
                for record in chunk {
                    let hash = T::hash_sequence(&record.seq);
                    if hash_tx.send(HashedRecord { hash, record }).is_err() {
                        return;
                    }
                }
            }
        }));
    }

    // Drop original senders so channels close when producers finish
    drop(chunk_tx);
    drop(hash_tx);

    // --- WRITER (this thread) — runs concurrently with reader & workers ---
    let output_path_str = output_path;
    let output_path = Path::new(output_path_str);
    let gzip_threads = if num_workers > 1 {
        num_workers.saturating_sub(2).max(1)
    } else {
        0
    };
    let (writer, output_format, gzip_handle) = create_output_writer(
        output_path_str,
        cfg.compression_level,
        gzip_threads,
        cfg.out_format,
        cfg.out_compression,
        cfg.force,
    )?;
    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

    // Match the serial path: refuse a FASTA -> FASTQ conversion before any
    // record is written, rather than emitting FASTQ with an empty quality line.
    // `recv` fails only on an empty input (reader closed without a record), in
    // which case there is nothing to validate and nothing to write.
    if let Ok(has_quality) = fmt_rx.recv() {
        validate_format_compatibility(has_quality, output_format)?;
    }

    let mut duplicates = 0;
    let mut processed = 0;

    if cfg.verbose {
        println!(
            "Multi-threaded mode: {} workers, {} shards, {} compression threads",
            worker_threads, worker_threads, gzip_threads
        );
    }

    // Receive and process results from workers
    while let Ok(hr) = hash_rx.recv() {
        processed += 1;
        if verifier.verify(hr.hash) {
            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buffered_writer, &hr.record.id, &hr.record.seq)?;
                }
                OutputFormat::Fastq => {
                    let nl = ending(hr.record.crlf);
                    buffered_writer.write_all(b"@")?;
                    buffered_writer.write_all(&hr.record.id)?;
                    buffered_writer.write_all(nl)?;
                    buffered_writer.write_all(&hr.record.seq)?;
                    buffered_writer.write_all(nl)?;
                    buffered_writer.write_all(b"+")?;
                    buffered_writer.write_all(nl)?;
                    buffered_writer.write_all(&hr.record.qual)?;
                    buffered_writer.write_all(nl)?;
                }
            }
        } else {
            duplicates += 1;
        }
    }

    // Verify all output is flushed
    if processed > 0 {
        buffered_writer.flush().context("Failed to flush output")?;
    }
    // Drop the buffered writer to close the channel before joining the gzip thread
    drop(buffered_writer);
    wait_for_gzip_thread(gzip_handle, "gzippy compression")?;

    // Join reader to check for errors and get total count
    let total_reads = match reader_handle.join() {
        Ok(Ok(count)) => count,
        Ok(Err(e)) => return Err(e),
        Err(e) => bail!("Reader thread panicked: {:?}", e),
    };

    // Worker handles must finish since hash_rx closed
    for h in worker_handles {
        h.join().unwrap();
    }

    if cfg.verbose && OutputFormat::is_gz(output_path) {
        println!(
            "Output file written (gzip compressed): {}",
            output_path.display()
        );
    }

    Ok((total_reads, duplicates))
}

// ==========================================
// PAIRED-END PARALLEL
// ==========================================

/// Executes paired-end deduplication using a multi-threaded pipeline.
///
/// R1 and R2 inputs are streamed through flate2 for decompression (single-threaded).
/// When output paths are `.gz` and `compression_level >= 6`, each output stream is
/// compressed in parallel via gzippy (each gets `max(1, num_workers / 2)` threads);
/// otherwise flate2 single-threaded.
pub(crate) fn execute_paired_deduplication_parallel<T: SequenceHasher + Send + 'static>(
    input_r1_path: &str,
    input_r2_path: &str,
    output_r1_path: &str,
    output_r2_path: &str,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    // --dry-run writes nothing: the serial counting pass reaches the identical
    // verdict (same refusals, same duplicate count) without creating a file.
    if cfg.dry_run {
        return crate::processor::execute_paired_deduplication::<T>(
            input_r1_path,
            input_r2_path,
            output_r1_path,
            output_r2_path,
            cfg,
        );
    }

    let num_workers = resolve_thread_count(cfg.threads);
    // FD-3, paired twin: 1 reader + (N-1) workers + this writer thread + two
    // pools of max(1, (N-2)/2) compression threads each — 2N-1 threads for even
    // N, the same split as the single-end path. `-j` budgets the workers.
    let worker_threads = num_workers.saturating_sub(1).max(1);

    // Resume existing outputs before any thread starts, so the preload pass
    // cannot race the unbounded worker channel.
    let mut verifier = ShardedHashVerifier::<T>::new(cfg.estimated_capacity, worker_threads)?;
    resume_existing_paired_outputs(output_r1_path, output_r2_path, cfg, &mut verifier)?;

    let (chunk_tx, chunk_rx): (Sender<Vec<RawPair>>, Receiver<Vec<RawPair>>) =
        bounded(num_workers * 2);
    let (hash_tx, hash_rx): (Sender<HashedPair<T>>, Receiver<HashedPair<T>>) =
        crossbeam_channel::unbounded();

    // Same as the single-end path: the reader reports whether R1 carries
    // quality scores, and the writer validates the format pairing before
    // writing, so the threaded verdict matches the serial one.
    let (fmt_tx, fmt_rx): (Sender<bool>, Receiver<bool>) = bounded(1);

    // --- READER THREAD ---
    let reader_handle = std::thread::spawn({
        let chunk_tx = chunk_tx.clone();
        let r1_owned = input_r1_path.to_string();
        let r2_owned = input_r2_path.to_string();
        move || -> Result<usize> {
            let reader1 = open_input_reader(&r1_owned)?;
            let reader2 = open_input_reader(&r2_owned)?;

            let mut reader_r1 =
                parse_fastx_reader(reader1).map_err(|e| anyhow::anyhow!("Parse error: {e}"))?;
            let mut reader_r2 =
                parse_fastx_reader(reader2).map_err(|e| anyhow::anyhow!("Parse error: {e}"))?;

            let mut total = 0;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);
            let mut format_reported = false;

            while let Some(r1) = reader_r1.next() {
                let r1 = r1.context("Invalid R1 sequence")?;
                let r2 = match reader_r2.next() {
                    Some(r) => r.context("Invalid R2 sequence")?,
                    None => {
                        bail!("Desynchronization: R1 has more records than R2");
                    }
                };

                // The serial PE path validates against R1's quality only; do the
                // same here so both paths reach the same verdict.
                if !format_reported {
                    let has_quality = r1.qual().is_some_and(|q| !q.is_empty());
                    let _ = fmt_tx.send(has_quality);
                    format_reported = true;
                }

                // Validate pair synchronization
                let base_r1 = base_read_id(r1.id()).to_vec();
                let base_r2 = base_read_id(r2.id()).to_vec();
                if base_r1 != base_r2 {
                    bail!(
                        "Desynchronization: R1 base ID {:?} != R2 base ID {:?}",
                        String::from_utf8_lossy(&base_r1),
                        String::from_utf8_lossy(&base_r2)
                    );
                }

                let pair = RawPair {
                    r1: RawRecord {
                        id: r1.id().to_vec(),
                        seq: r1.seq().to_vec(),
                        qual: r1.qual().map(|q| q.to_vec()).unwrap_or_default(),
                        crlf: matches!(r1.line_ending(), LineEnding::Windows),
                    },
                    r2: RawRecord {
                        id: r2.id().to_vec(),
                        seq: r2.seq().to_vec(),
                        qual: r2.qual().map(|q| q.to_vec()).unwrap_or_default(),
                        crlf: matches!(r2.line_ending(), LineEnding::Windows),
                    },
                };

                chunk.push(pair);
                total += 1;

                if chunk.len() >= CHUNK_SIZE {
                    if chunk_tx.send(chunk).is_err() {
                        return Ok(total);
                    }
                    chunk = Vec::with_capacity(CHUNK_SIZE);
                }
            }

            // Check for R2 having extra records
            if reader_r2.next().is_some() {
                bail!("Desynchronization: R2 has more records than R1");
            }

            if !chunk.is_empty() {
                let _ = chunk_tx.send(chunk);
            }

            Ok(total)
        }
    });

    // --- WORKER THREADS ---
    let mut worker_handles = Vec::with_capacity(worker_threads);
    for _ in 0..worker_threads {
        let chunk_rx = chunk_rx.clone();
        let hash_tx = hash_tx.clone();
        worker_handles.push(std::thread::spawn(move || {
            while let Ok(chunk) = chunk_rx.recv() {
                for pair in chunk {
                    let hash = T::hash_pair(&pair.r1.seq, &pair.r2.seq);
                    if hash_tx.send(HashedPair { hash, pair }).is_err() {
                        return;
                    }
                }
            }
        }));
    }

    drop(chunk_tx);
    drop(hash_tx);

    // --- WRITER (this thread) — runs concurrently with reader & workers ---
    // Two pools share what the workers left, one per mate.
    let gzip_threads = if num_workers > 1 {
        (num_workers.saturating_sub(2) / 2).max(1)
    } else {
        0
    };
    let (writer_r1, output_format, gzip_handle1) = create_output_writer(
        output_r1_path,
        cfg.compression_level,
        gzip_threads,
        cfg.out_format,
        cfg.out_compression,
        cfg.force,
    )?;
    let (writer_r2, _output_format_r2, gzip_handle2) = create_output_writer(
        output_r2_path,
        cfg.compression_level,
        gzip_threads,
        cfg.out_format,
        cfg.out_compression,
        cfg.force,
    )?;

    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buf_r1 = BufWriter::with_capacity(buf_size, writer_r1);
    let mut buf_r2 = BufWriter::with_capacity(buf_size, writer_r2);

    // Match the serial PE path: refuse a FASTA -> FASTQ conversion before any
    // record is written, rather than emitting FASTQ with an empty quality line.
    if let Ok(has_quality) = fmt_rx.recv() {
        validate_format_compatibility(has_quality, output_format)?;
    }

    let mut duplicates = 0;
    let mut processed = 0;

    if cfg.verbose {
        println!(
            "Multi-threaded PE mode: {} workers, {} shards, {} compression threads per mate",
            worker_threads, worker_threads, gzip_threads
        );
    }

    while let Ok(hp) = hash_rx.recv() {
        processed += 1;
        if verifier.verify(hp.hash) {
            match output_format {
                OutputFormat::Fasta => {
                    write_fasta_record(&mut buf_r1, &hp.pair.r1.id, &hp.pair.r1.seq)?;
                    write_fasta_record(&mut buf_r2, &hp.pair.r2.id, &hp.pair.r2.seq)?;
                }
                OutputFormat::Fastq => {
                    let nl1 = ending(hp.pair.r1.crlf);
                    let nl2 = ending(hp.pair.r2.crlf);
                    buf_r1.write_all(b"@")?;
                    buf_r1.write_all(&hp.pair.r1.id)?;
                    buf_r1.write_all(nl1)?;
                    buf_r1.write_all(&hp.pair.r1.seq)?;
                    buf_r1.write_all(nl1)?;
                    buf_r1.write_all(b"+")?;
                    buf_r1.write_all(nl1)?;
                    buf_r1.write_all(&hp.pair.r1.qual)?;
                    buf_r1.write_all(nl1)?;

                    buf_r2.write_all(b"@")?;
                    buf_r2.write_all(&hp.pair.r2.id)?;
                    buf_r2.write_all(nl2)?;
                    buf_r2.write_all(&hp.pair.r2.seq)?;
                    buf_r2.write_all(nl2)?;
                    buf_r2.write_all(b"+")?;
                    buf_r2.write_all(nl2)?;
                    buf_r2.write_all(&hp.pair.r2.qual)?;
                    buf_r2.write_all(nl2)?;
                }
            }
        } else {
            duplicates += 1;
        }
    }

    if processed > 0 {
        buf_r1.flush()?;
        buf_r2.flush()?;
    }
    // Drop buffered writers to close the channels before joining gzip threads
    drop(buf_r1);
    drop(buf_r2);
    wait_for_gzip_thread(gzip_handle1, "gzippy compression (R1)")?;
    wait_for_gzip_thread(gzip_handle2, "gzippy compression (R2)")?;

    // Join reader to check for errors and get total count
    let total_reads = match reader_handle.join() {
        Ok(Ok(count)) => count,
        Ok(Err(e)) => return Err(e),
        Err(e) => bail!("Reader thread panicked: {:?}", e),
    };

    for h in worker_handles {
        h.join().unwrap();
    }

    Ok((total_reads, duplicates))
}

// ==========================================
// HELPERS
// ==========================================

fn resolve_thread_count(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
    } else {
        requested
    }
}
