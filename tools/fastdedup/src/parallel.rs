use std::hash::{BuildHasher, Hash};
use std::io::Read;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result, bail};
use crossbeam_channel::{Receiver, Sender, bounded};
use needletail::parse_fastx_reader;
use rustc_hash::{FxBuildHasher, FxHashSet};

use crate::hasher::SequenceHasher;
use crate::utils::{base_read_id, compute_write_buffer_size, validate_format_compatibility};
use crate::writer::{OutputFormat, prepare_writer, write_fasta_record, StreamCompression};
use std::fs::File;
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
#[derive(Clone, Debug)]
pub(crate) struct RawRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
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

/// Writer side: implements `Write` by sending chunks through a crossbeam channel.
struct ChanWriter {
    tx: Sender<Vec<u8>>,
}

impl Write for ChanWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.tx.send(buf.to_vec()).map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::BrokenPipe, "gzippy compression thread died")
        })?;
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
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

/// Open a FASTX input as a byte stream, using flate2 GzDecoder for gzip
/// decompression when the path ends with `.gz` (streaming, no subprocess).
fn open_input_reader(path: &str) -> Result<Box<dyn Read + Send + 'static>> {
    if crate::stream::is_stream(path) {
        // needletail sniffs gzip from magic bytes on any reader, so a
        // stdin stream decodes exactly like a mislabelled file.
        return Ok(Box::new(std::io::stdin()));
    }
    if path.ends_with(".gz") {
        let file = File::open(path).with_context(|| format!("Failed to open input file: {path}"))?;
        Ok(Box::new(flate2::read::GzDecoder::new(file)))
    } else {
        let file = File::open(path).with_context(|| format!("Failed to open input file: {path}"))?;
        Ok(Box::new(file))
    }
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
) -> Result<(
    Box<dyn Write>,
    OutputFormat,
    Option<std::thread::JoinHandle<Result<()>>>,
)> {
    if crate::stream::is_stream(path) {
        // D-4: stdout cannot be seeked or reopened, so there is no gzippy
        // thread and no resume — serial flate2 gzip when asked for.
        let writer = crate::writer::prepare_stream_writer(out_compression, compression_level);
        return Ok((writer, out_format, None));
    }
    let output_path = Path::new(path);
    let output_format = OutputFormat::from_extension(output_path);
    let is_gz = OutputFormat::is_gz(output_path);

    if is_gz && num_threads > 1 && compression_level >= 6 {
        // Gzippy parallel compression path — channel-based streaming pipe
        // Levels 6+ with threads>1 use PipelinedGzEncoder → standard gzip
        let (tx, rx) = bounded::<Vec<u8>>(16);
        let writer: Box<dyn Write> = Box::new(ChanWriter { tx });

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
        let (writer, format) = prepare_writer(output_path, true, compression_level)?;
        Ok((writer, format, None))
    }
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
    let num_workers = resolve_thread_count(cfg.threads);

    let (chunk_tx, chunk_rx): (Sender<Vec<RawRecord>>, Receiver<Vec<RawRecord>>) =
        bounded(num_workers * 2);
    let (hash_tx, hash_rx): (Sender<HashedRecord<T>>, Receiver<HashedRecord<T>>) =
        crossbeam_channel::unbounded();

    let has_quality_flag = Arc::new(AtomicBool::new(false));

    // --- READER THREAD ---
    let reader_handle = std::thread::spawn({
        let chunk_tx = chunk_tx.clone();
        let input_owned = input_path.to_string();
        let has_quality_flag = Arc::clone(&has_quality_flag);
        move || -> Result<usize> {
            let reader = open_input_reader(&input_owned)?;
            let mut fastx =
                parse_fastx_reader(reader).map_err(|e| anyhow::anyhow!("Parse error: {e}"))?;

            let mut total = 0;
            let mut has_quality = false;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);

            while let Some(record) = fastx.next() {
                let rec = record.context("Invalid sequence")?;
                let id = rec.id().to_vec();
                let seq = rec.seq().to_vec();
                let qual = rec.qual().map(|q| q.to_vec()).unwrap_or_default();

                if !has_quality {
                    has_quality = rec.qual().is_some_and(|q| !q.is_empty());
                }

                chunk.push(RawRecord { id, seq, qual });
                total += 1;

                if chunk.len() >= CHUNK_SIZE {
                    if chunk_tx.send(chunk).is_err() {
                        has_quality_flag.store(has_quality, Ordering::Relaxed);
                        return Ok(total);
                    }
                    chunk = Vec::with_capacity(CHUNK_SIZE);
                }
            }
            if !chunk.is_empty() {
                let _ = chunk_tx.send(chunk);
            }
            has_quality_flag.store(has_quality, Ordering::Relaxed);

            Ok(total)
        }
    });

    // --- WORKER THREADS ---
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
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
    let gzip_threads = if num_workers > 1 { num_workers } else { 0 };
    let (writer, output_format, gzip_handle) = create_output_writer(
        output_path_str,
        cfg.compression_level,
        gzip_threads,
        cfg.out_format,
        cfg.out_compression,
    )?;
    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buffered_writer = BufWriter::with_capacity(buf_size, writer);

    let output_has_quality = output_format == OutputFormat::Fastq;
    validate_format_compatibility(output_has_quality, output_format)?;

    let mut verifier = ShardedHashVerifier::<T>::new(cfg.estimated_capacity, num_workers)?;
    let mut duplicates = 0;
    let mut processed = 0;

    if cfg.verbose {
        println!(
            "Multi-threaded mode: {} workers, {} shards",
            num_workers, num_workers
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
                    buffered_writer.write_all(b"@")?;
                    buffered_writer.write_all(&hr.record.id)?;
                    buffered_writer.write_all(b"\n")?;
                    buffered_writer.write_all(&hr.record.seq)?;
                    buffered_writer.write_all(b"\n+\n")?;
                    buffered_writer.write_all(&hr.record.qual)?;
                    buffered_writer.write_all(b"\n")?;
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
    let num_workers = resolve_thread_count(cfg.threads);

    let (chunk_tx, chunk_rx): (Sender<Vec<RawPair>>, Receiver<Vec<RawPair>>) =
        bounded(num_workers * 2);
    let (hash_tx, hash_rx): (Sender<HashedPair<T>>, Receiver<HashedPair<T>>) =
        crossbeam_channel::unbounded();

    let has_quality_flag = Arc::new(AtomicBool::new(false));

    // --- READER THREAD ---
    let reader_handle = std::thread::spawn({
        let chunk_tx = chunk_tx.clone();
        let r1_owned = input_r1_path.to_string();
        let r2_owned = input_r2_path.to_string();
        let has_quality_flag = Arc::clone(&has_quality_flag);
        move || -> Result<usize> {
            let reader1 = open_input_reader(&r1_owned)?;
            let reader2 = open_input_reader(&r2_owned)?;

            let mut reader_r1 =
                parse_fastx_reader(reader1).map_err(|e| anyhow::anyhow!("Parse error: {e}"))?;
            let mut reader_r2 =
                parse_fastx_reader(reader2).map_err(|e| anyhow::anyhow!("Parse error: {e}"))?;

            let mut total = 0;
            let mut has_quality = false;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);

            while let Some(r1) = reader_r1.next() {
                let r1 = r1.context("Invalid R1 sequence")?;
                let r2 = match reader_r2.next() {
                    Some(r) => r.context("Invalid R2 sequence")?,
                    None => {
                        has_quality_flag.store(has_quality, Ordering::Relaxed);
                        bail!("Desynchronization: R1 has more records than R2");
                    }
                };

                // Validate pair synchronization
                let base_r1 = base_read_id(r1.id()).to_vec();
                let base_r2 = base_read_id(r2.id()).to_vec();
                if base_r1 != base_r2 {
                    has_quality_flag.store(has_quality, Ordering::Relaxed);
                    bail!(
                        "Desynchronization: R1 base ID {:?} != R2 base ID {:?}",
                        String::from_utf8_lossy(&base_r1),
                        String::from_utf8_lossy(&base_r2)
                    );
                }

                if !has_quality {
                    has_quality = r1.qual().is_some_and(|q| !q.is_empty());
                }

                let pair = RawPair {
                    r1: RawRecord {
                        id: r1.id().to_vec(),
                        seq: r1.seq().to_vec(),
                        qual: r1.qual().map(|q| q.to_vec()).unwrap_or_default(),
                    },
                    r2: RawRecord {
                        id: r2.id().to_vec(),
                        seq: r2.seq().to_vec(),
                        qual: r2.qual().map(|q| q.to_vec()).unwrap_or_default(),
                    },
                };

                chunk.push(pair);
                total += 1;

                if chunk.len() >= CHUNK_SIZE {
                    if chunk_tx.send(chunk).is_err() {
                        has_quality_flag.store(has_quality, Ordering::Relaxed);
                        return Ok(total);
                    }
                    chunk = Vec::with_capacity(CHUNK_SIZE);
                }
            }

            // Check for R2 having extra records
            if reader_r2.next().is_some() {
                has_quality_flag.store(has_quality, Ordering::Relaxed);
                bail!("Desynchronization: R2 has more records than R1");
            }

            if !chunk.is_empty() {
                let _ = chunk_tx.send(chunk);
            }
            has_quality_flag.store(has_quality, Ordering::Relaxed);

            Ok(total)
        }
    });

    // --- WORKER THREADS ---
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
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
    let gzip_threads = if num_workers > 1 {
        (num_workers / 2).max(1)
    } else {
        0
    };
    let (writer_r1, output_format, gzip_handle1) = create_output_writer(
        output_r1_path,
        cfg.compression_level,
        gzip_threads,
        cfg.out_format,
        cfg.out_compression,
    )?;
    let (writer_r2, _output_format_r2, gzip_handle2) = create_output_writer(
        output_r2_path,
        cfg.compression_level,
        gzip_threads,
        cfg.out_format,
        cfg.out_compression,
    )?;

    let buf_size = compute_write_buffer_size(cfg.read_length);
    let mut buf_r1 = BufWriter::with_capacity(buf_size, writer_r1);
    let mut buf_r2 = BufWriter::with_capacity(buf_size, writer_r2);

    let output_has_quality = output_format == OutputFormat::Fastq;
    validate_format_compatibility(output_has_quality, output_format)?;

    let mut verifier = ShardedHashVerifier::<T>::new(cfg.estimated_capacity, num_workers)?;
    let mut duplicates = 0;
    let mut processed = 0;

    if cfg.verbose {
        println!(
            "Multi-threaded PE mode: {} workers, {} shards",
            num_workers, num_workers
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
                    buf_r1.write_all(b"@")?;
                    buf_r1.write_all(&hp.pair.r1.id)?;
                    buf_r1.write_all(b"\n")?;
                    buf_r1.write_all(&hp.pair.r1.seq)?;
                    buf_r1.write_all(b"\n+\n")?;
                    buf_r1.write_all(&hp.pair.r1.qual)?;
                    buf_r1.write_all(b"\n")?;

                    buf_r2.write_all(b"@")?;
                    buf_r2.write_all(&hp.pair.r2.id)?;
                    buf_r2.write_all(b"\n")?;
                    buf_r2.write_all(&hp.pair.r2.seq)?;
                    buf_r2.write_all(b"\n+\n")?;
                    buf_r2.write_all(&hp.pair.r2.qual)?;
                    buf_r2.write_all(b"\n")?;
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
