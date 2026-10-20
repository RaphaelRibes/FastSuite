//! Input and output plumbing.
//!
//! Input: `needletail::parse_fastx_file` handles plain/.gz/.bz2/.xz/.zst
//! transparently via its internal auto-detection, so we just hand it a
//! path.
//!
//! Output: we inspect the extension and wrap the file in the appropriate
//! encoder. Supported:
//!   * `.gz`            → BGZF (gzip-compatible), parallel, via gzp
//!   * `.zst` / `.zstd` → zstd
//!   * anything else    → plain
//!
//! Gzip output is BGZF, which is multi-member gzip with a block-size field.
//! That is ordinary `.gz`: `zcat`, `gzip -d`, `pigz` and htslib all read it,
//! needletail decodes input with `MultiGzDecoder`, and tabix can index it.
//! Note that a consumer using `flate2::read::GzDecoder` rather than
//! `MultiGzDecoder` would stop after the first member.
//!
//! Output buffers are 1 MiB. On HPC-class hardware the cost is trivial
//! and the throughput win on compressed output is significant: gzip and
//! zstd both benefit from larger write chunks because they batch
//! internal state flushes.

use anyhow::{Context, Result};
use gzp::deflate::Bgzf;
use gzp::par::compress::ParCompressBuilder;
use gzp::par::decompress::ParDecompressBuilder;
use gzp::Compression as GzCompression;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, SeekFrom, Write};
use std::path::Path;

const OUTPUT_BUFFER: usize = 1 << 20; // 1 MiB
const INPUT_BUFFER: usize = 1 << 20; // 1 MiB

/// A boxed `Write` the rest of the pipeline writes FASTQ/FASTA records to.
pub type OutputWriter = Box<dyn Write + Send>;

/// Record format, inferred per input record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordFormat {
    Fastq,
    Fasta,
}

/// Compression levels. Each format carries its own level; the unused
/// one is ignored on any given file.
#[derive(Debug, Clone, Copy)]
pub struct CompressionLevels {
    pub gz: u32,
    pub zst: i32,
    /// Worker threads for gzip output. Must be >= 1.
    ///
    /// Measured on a 10 M-pair paired-end run: single-threaded deflate was
    /// 73% of total wall clock with compressed output, and 7x the cost fastp
    /// paid for the same work. Every other stage was at parity or better, so
    /// this is the only knob that moves the headline number.
    ///
    /// Each open output stream gets its own pool, so a paired run with
    /// `gz_threads = n` spawns 2n compression threads on top of the Rayon
    /// pool. `Cli::finalize` sizes it against `--threads` accordingly.
    pub gz_threads: usize,
}

impl Default for CompressionLevels {
    fn default() -> Self {
        Self {
            gz: 6,
            zst: 3,
            gz_threads: 1,
        }
    }
}

/// Open an output path, inferring compression from the extension.
pub fn open_output(path: &Path, levels: CompressionLevels) -> Result<OutputWriter> {
    let file = File::create(path).with_context(|| format!("creating output {}", path.display()))?;
    let buf = BufWriter::with_capacity(OUTPUT_BUFFER, file);
    encode_output(buf, path, levels)
}

/// Like [`open_output`] but writes to `write_path` while inferring the
/// compression format from `format_path`'s extension. Used when the write
/// path is a `.part` temporary file whose suffix differs from the final name.
pub fn open_output_as(
    write_path: &Path,
    format_path: &Path,
    levels: CompressionLevels,
) -> Result<OutputWriter> {
    let file = File::create(write_path)
        .with_context(|| format!("creating output {}", write_path.display()))?;
    let buf = BufWriter::with_capacity(OUTPUT_BUFFER, file);
    encode_output(buf, format_path, levels)
}

fn encode_output(
    buf: BufWriter<File>,
    format_path: &Path,
    levels: CompressionLevels,
) -> Result<OutputWriter> {
    let ext = format_path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "gz" => {
            // Parallel deflate, BGZF framing, libdeflate backend.
            //
            // BGZF is a strict subset of gzip: zcat, gzip, pigz, needletail
            // and htslib all read it, and it is additionally indexable with
            // tabix. libdeflate only backs gzp's block formats, because it
            // is a non-streaming API, and libdeflate is where the speed is:
            // measured at -j 1, serial zlib-ng deflate cost 160.6 s against
            // fastp's libdeflate at 59.3 s for the same 10 M-pair output.
            //
            // Cost: BGZF blocks are capped at 64 KiB, so files run a few
            // percent larger than large-block gzip.
            //
            // ParCompress finishes the stream in its Drop impl, which matters
            // because the pipeline drops writers before the atomic rename. A
            // half-written member can therefore not be renamed into place.
            let enc = ParCompressBuilder::<Bgzf>::new()
                .num_threads(levels.gz_threads.max(1))
                .context("configuring parallel gzip")?
                .compression_level(GzCompression::new(levels.gz))
                .from_writer(buf);
            Ok(Box::new(enc))
        }
        "zst" | "zstd" => {
            let enc = zstd::stream::write::Encoder::new(buf, levels.zst)
                .context("initialising zstd encoder")?
                .auto_finish();
            Ok(Box::new(enc))
        }
        _ => Ok(Box::new(buf)),
    }
}

/// A boxed `Read` the pipeline feeds into needletail. Parallel BGZF
/// decompression when the input is BGZF-structured, plain file otherwise
/// (needletail auto-detects gzip/bz2/xz/zstd itself).
pub type InputReader = Box<dyn Read + Send>;

/// BGZF magic: gzip header with FEXTRA flag and a `BC` subfield.
fn looks_like_bgzf(buf: &[u8]) -> bool {
    buf.len() >= 18
        && buf[0] == 0x1f
        && buf[1] == 0x8b
        && buf[3] & 0x04 != 0 // FEXTRA
        && buf[12] == b'B'
        && buf[13] == b'C'
}

/// Open an input path for reading, using parallel BGZF decompression when
/// the file is BGZF-structured (the tool's own `.gz` outputs are), falling
/// back to a plain buffered file so needletail's auto-detection handles
/// ordinary gzip/bz2/xz/zstd.
///
/// `gz_threads` is the size of the decompression pool; 0/1 selects the
/// single-threaded path (gzp's builder treats 0 as serial via
/// `maybe_par_from_reader`).
pub fn open_input(path: &Path, gz_threads: usize) -> Result<InputReader> {
    let mut file = File::open(path).with_context(|| format!("opening input {}", path.display()))?;
    let mut head = [0u8; 18];
    let n = std::io::Read::read(&mut file, &mut head)
        .with_context(|| format!("reading header of {}", path.display()))?;
    let _ = std::io::Seek::seek(&mut file, SeekFrom::Start(0));
    let head = &head[..n];

    if looks_like_bgzf(head) && gz_threads > 1 {
        let dec = ParDecompressBuilder::<Bgzf>::new()
            .num_threads(gz_threads)
            .context("configuring parallel gzip decompression")?
            .from_reader(file);
        Ok(Box::new(dec))
    } else {
        Ok(Box::new(BufReader::with_capacity(INPUT_BUFFER, file)))
    }
}

/// Flush threshold for the record staging buffer. A record is appended in
/// memory and handed to the compressor in one call once the buffer passes
/// this size, instead of six short calls per record.
pub const STAGE_FLUSH: usize = 256 << 10;

/// Record writer with a staging buffer.
///
/// Every FASTQ/FASTA output is one of these. The point is the syscall-free
/// append path: six `write_all` calls per record became one per ~256 KiB, and
/// each `ParCompress::write` (which memcpys into a `BytesMut` and re-checks the
/// 65 280-byte BGZF block boundary) now runs once per buffer rather than once
/// per record plus its quality slice. It also hands the compressor whole
/// blocks, which is what it wants.
///
/// Quality bytes are written exactly as ingested: the Phred offset is folded
/// into the quality-gate LUT and the mean sum, never into the data (see
/// `algorithm::Params::prepare`), so this is a pure append and the output
/// bytes are the input bytes.
pub struct Out {
    w: OutputWriter,
    stage: Vec<u8>,
}

impl Out {
    pub fn new(w: OutputWriter) -> Self {
        // Slack for the record that crosses the threshold, so a spill never
        // has to grow the buffer.
        Self {
            w,
            stage: Vec::with_capacity(STAGE_FLUSH + (8 << 10)),
        }
    }

    /// Append one FASTQ record: `@id\nseq\n+\nqual\n`.
    pub fn fastq(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        debug_assert_eq!(seq.len(), qual.len());
        self.stage.push(b'@');
        self.stage.extend_from_slice(strip_sigil(id, b'@'));
        self.stage.push(b'\n');
        self.stage.extend_from_slice(seq);
        self.stage.extend_from_slice(b"\n+\n");
        self.stage.extend_from_slice(qual);
        self.stage.push(b'\n');
        if self.stage.len() >= STAGE_FLUSH {
            self.spill()?;
        }
        Ok(())
    }

    /// Append one FASTA record: `>id\nseq\n`. Single-line FASTA (no 60/80-char
    /// wrapping) — this matches needletail's output convention and is what
    /// every downstream tool accepts.
    pub fn fasta(&mut self, id: &[u8], seq: &[u8]) -> Result<()> {
        self.stage.push(b'>');
        self.stage.extend_from_slice(strip_sigil(id, b'>'));
        self.stage.push(b'\n');
        self.stage.extend_from_slice(seq);
        self.stage.push(b'\n');
        if self.stage.len() >= STAGE_FLUSH {
            self.spill()?;
        }
        Ok(())
    }

    /// Hand the staged bytes to the compressor.
    fn spill(&mut self) -> Result<()> {
        if !self.stage.is_empty() {
            self.w.write_all(&self.stage)?;
            self.stage.clear();
        }
        Ok(())
    }

    /// Push everything staged, then flush the compressor. Must run before the
    /// atomic rename that commits the file.
    pub fn flush(&mut self) -> Result<()> {
        self.spill()?;
        self.w.flush()?;
        Ok(())
    }
}

fn strip_sigil<'a>(s: &'a [u8], sigil: u8) -> &'a [u8] {
    if s.first() == Some(&sigil) {
        &s[1..]
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;
    use std::path::PathBuf;

    use std::sync::atomic::{AtomicU64, Ordering};

    static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn tmp_path(name: &str) -> PathBuf {
        let id = TMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut p = std::env::temp_dir();
        p.push(format!("polyx-test-{}-{}-{}", std::process::id(), id, name));
        p
    }

    #[test]
    fn strip_sigil_removes_at_and_gt() {
        assert_eq!(strip_sigil(b"@read1", b'@'), b"read1");
        assert_eq!(strip_sigil(b">contig", b'>'), b"contig");
        assert_eq!(strip_sigil(b"read1", b'@'), b"read1");
        assert_eq!(strip_sigil(b"", b'@'), b"");
    }

    #[test]
    fn write_fastq_plain_roundtrip() {
        let p = tmp_path("plain.fq");
        {
            let mut w = Out::new(open_output(&p, CompressionLevels::default()).unwrap());
            w.fastq(b"read1", b"ACGT", b"IIII").unwrap();
            w.fastq(b"@read2", b"NN", b"##").unwrap();
            w.flush().unwrap();
        }
        let mut data = String::new();
        File::open(&p).unwrap().read_to_string(&mut data).unwrap();
        assert_eq!(data, "@read1\nACGT\n+\nIIII\n@read2\nNN\n+\n##\n");
        std::fs::remove_file(&p).ok();
    }

    /// Quality bytes are written verbatim: the Phred offset lives in the
    /// quality-gate LUT and the mean sum, never in the data, so egress is an
    /// append and not a per-byte re-encode.
    #[test]
    fn write_fastq_quality_bytes_are_verbatim() {
        let p = tmp_path("rawq.fq");
        {
            let mut w = Out::new(open_output(&p, CompressionLevels::default()).unwrap());
            w.fastq(b"r1", b"ACGT", b"IIII").unwrap();
            // A byte below any sane offset must round-trip untouched: the
            // clamp happens once, on ingest, where the LUT is indexed.
            w.fastq(b"r2", b"ACGT", b"!I#\x7f").unwrap();
            w.flush().unwrap();
        }
        let mut data = Vec::new();
        File::open(&p).unwrap().read_to_end(&mut data).unwrap();
        assert_eq!(data, b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n!I#\x7f\n");
        std::fs::remove_file(&p).ok();
    }

    /// A `Write` that records every call, to prove the record path batches.
    #[derive(Clone, Default)]
    struct Counting(std::sync::Arc<std::sync::Mutex<(Vec<u8>, usize)>>);

    impl Write for Counting {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            let mut g = self.0.lock().unwrap();
            g.0.extend_from_slice(buf);
            g.1 += 1;
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// Six `write_all`s per record became one per staging flush. 20 000
    /// records are ~2.2 MB, so a 256 KiB stage must need a handful of calls —
    /// the old code emitted 120 000.
    #[test]
    fn fastq_records_are_staged_not_written_one_call_per_field() {
        let sink = Counting::default();
        let mut w = Out::new(Box::new(sink.clone()));
        let seq = vec![b'A'; 50];
        let qual = vec![b'I'; 50];
        for i in 0..20_000 {
            w.fastq(format!("r{i}").as_bytes(), &seq, &qual).unwrap();
        }
        w.flush().unwrap();
        let (bytes, calls) = { sink.0.lock().unwrap().clone() };
        let expected_flushes = bytes.len() / STAGE_FLUSH + 1;
        assert!(
            calls <= expected_flushes + 1,
            "{calls} write calls for {} bytes; staging should need ~{expected_flushes}",
            bytes.len()
        );
        let mut expect = Vec::new();
        for i in 0..20_000 {
            expect.extend_from_slice(format!("@r{i}\n").as_bytes());
            expect.extend_from_slice(&seq);
            expect.extend_from_slice(b"\n+\n");
            expect.extend_from_slice(&qual);
            expect.push(b'\n');
        }
        assert_eq!(bytes, expect, "staged bytes differ from the record stream");
    }

    #[test]
    fn fasta_records_are_staged() {
        let sink = Counting::default();
        let mut w = Out::new(Box::new(sink.clone()));
        let seq = vec![b'A'; 100];
        for i in 0..10_000 {
            w.fasta(format!("c{i}").as_bytes(), &seq).unwrap();
        }
        w.flush().unwrap();
        let (bytes, calls) = { sink.0.lock().unwrap().clone() };
        assert!(calls <= bytes.len() / STAGE_FLUSH + 2, "{calls} write calls");
        assert!(bytes.starts_with(b">c0\nAAAA"), "content corrupted");
        assert!(bytes.ends_with(b"\n"));
    }

    /// Exercises the parallel path specifically: with more than one worker
    /// and enough data to fill several blocks, the output is multi-member,
    /// which is the case a single-member decoder would silently truncate.
    #[test]
    fn write_fastq_gz_multithreaded_roundtrip() {
        let p = tmp_path("rt-par.fq.gz");
        let levels = CompressionLevels {
            gz: 6,
            zst: 3,
            gz_threads: 4,
        };
        let seq = vec![b'A'; 200];
        let qual = vec![b'I'; 200];
        {
            let mut w = Out::new(open_output(&p, levels).unwrap());
            for i in 0..20_000 {
                w.fastq(format!("read{i}").as_bytes(), &seq, &qual).unwrap();
            }
            w.flush().unwrap();
        }
        let f = File::open(&p).unwrap();
        let mut dec = MultiGzDecoder::new(f);
        let mut s = String::new();
        dec.read_to_string(&mut s).unwrap();
        assert_eq!(
            s.lines().count(),
            20_000 * 4,
            "records lost across gzip members"
        );
        std::fs::remove_file(&p).ok();
    }

    /// BGZF input produced by the tool's own parallel writer must round-trip
    /// through `open_input` (parallel decompression path) identically.
    #[test]
    fn open_input_bgzf_roundtrip() {
        let p = tmp_path("bgzf-in.fq.gz");
        let levels = CompressionLevels {
            gz: 6,
            zst: 3,
            gz_threads: 4,
        };
        let seq = vec![b'A'; 150];
        let qual = vec![b'I'; 150];
        {
            let mut w = Out::new(open_output(&p, levels).unwrap());
            for i in 0..50_000 {
                w.fastq(format!("read{i}").as_bytes(), &seq, &qual).unwrap();
            }
            w.flush().unwrap();
        }
        // Must be detected as BGZF and take the parallel path.
        let mut head = [0u8; 18];
        let mut f = File::open(&p).unwrap();
        std::io::Read::read(&mut f, &mut head).unwrap();
        assert!(looks_like_bgzf(&head), "test file should be BGZF");

        let mut reader = open_input(&p, 4).unwrap();
        let mut all = Vec::new();
        reader.read_to_end(&mut all).unwrap();
        let text = String::from_utf8(all).unwrap();
        assert_eq!(
            text.lines().count(),
            50_000 * 4,
            "parallel decode lost records"
        );
        assert!(
            text.starts_with("@read0\n"),
            "content corrupted: {:?}",
            &text[..30]
        );
        std::fs::remove_file(&p).ok();
    }

    /// Plain (single-member) gzip must still decode through `open_input` via
    /// the needletail fallback (we hand back a plain file reader and
    /// needletail's gzip auto-detection takes over).
    #[test]
    fn open_input_plain_gz_fallback() {
        let p = tmp_path("plain-in.fq.gz");
        {
            let f = File::create(&p).unwrap();
            let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
            enc.finish().unwrap();
        }
        // NOT BGZF: must NOT take the parallel path; needletail decodes.
        let mut head = [0u8; 18];
        let mut f = File::open(&p).unwrap();
        std::io::Read::read(&mut f, &mut head).unwrap();
        assert!(
            !looks_like_bgzf(&head),
            "plain gzip must not look like BGZF"
        );

        let reader = open_input(&p, 4).unwrap();
        let mut parser = needletail::parse_fastx_reader(reader).unwrap();
        let rec = parser.next().unwrap().unwrap();
        assert_eq!(rec.seq(), &b"ACGT"[..]);
        assert!(parser.next().is_none(), "exactly one record");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn write_fastq_gz_roundtrip() {
        let p = tmp_path("rt.fq.gz");
        {
            let mut w = Out::new(open_output(&p, CompressionLevels::default()).unwrap());
            w.fastq(b"read1", b"ACGT", b"IIII").unwrap();
            w.flush().unwrap();
        }
        let f = File::open(&p).unwrap();
        let mut dec = MultiGzDecoder::new(f);
        let mut s = String::new();
        dec.read_to_string(&mut s).unwrap();
        assert_eq!(s, "@read1\nACGT\n+\nIIII\n");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn write_fastq_zst_roundtrip() {
        let p = tmp_path("rt.fq.zst");
        {
            let mut w = Out::new(open_output(&p, CompressionLevels::default()).unwrap());
            w.fastq(b"read1", b"ACGT", b"IIII").unwrap();
            w.flush().unwrap();
        }
        let f = File::open(&p).unwrap();
        let mut dec = zstd::stream::read::Decoder::new(f).unwrap();
        let mut s = String::new();
        dec.read_to_string(&mut s).unwrap();
        assert_eq!(s, "@read1\nACGT\n+\nIIII\n");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn write_fasta_plain() {
        let p = tmp_path("plain.fa");
        {
            let mut w = Out::new(open_output(&p, CompressionLevels::default()).unwrap());
            w.fasta(b">contig1", b"ACGT").unwrap();
            w.fasta(b">contig2", b"NN").unwrap();
            w.flush().unwrap();
        }
        let mut data = String::new();
        File::open(&p).unwrap().read_to_string(&mut data).unwrap();
        assert_eq!(data, ">contig1\nACGT\n>contig2\nNN\n");
        std::fs::remove_file(&p).ok();
    }
}
