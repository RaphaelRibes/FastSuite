//! Input: one pass, BGZF-parallel when possible, with explicit budgets.
//!
//! Supported inputs: plain FASTQ and `.gz` FASTQ. A `.gz` whose header says
//! BGZF (gzip FEXTRA + `BC` subfield — what Illumina writes and what
//! fastprep/fastdedup produce) is inflated block-parallel across
//! `gz_threads`. Any other `.gz` is inflated serially, which the report
//! records per file so a slow case is visible rather than mysterious. Other
//! extensions are refused next to the tool rather than half-supported.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use gzp::deflate::Bgzf;
use gzp::par::decompress::ParDecompressBuilder;

use crate::accum::Accum;

const INPUT_BUFFER: usize = 1 << 20;

/// Per-file reading budgets. Both stop the source, not the analysis, and the
/// report declares the run partial with the exact reason.
#[derive(Clone, Copy, Default, Debug)]
pub struct Limits {
    pub max_reads: Option<u64>,
    pub max_bytes: Option<u64>,
}

#[derive(Debug)]
pub struct ScanResult {
    /// Compressed bytes consumed for this file (equals file size when unbudgeted).
    pub bytes_read: u64,
    /// `None` when the whole file was read.
    pub stopped: Option<&'static str>,
    pub decompress: &'static str,
}

/// Counts bytes pulled from the underlying file, shared with the decoder so
/// the budget can be enforced mid-stream.
#[derive(Clone)]
struct ByteCounter(Arc<AtomicU64>);

struct Counting<R> {
    inner: R,
    n: Arc<AtomicU64>,
}

impl<R: Read> Read for Counting<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let k = self.inner.read(buf)?;
        self.n.fetch_add(k as u64, Ordering::Relaxed);
        Ok(k)
    }
}

/// BGZF magic: gzip header with the FEXTRA flag and a `BC` subfield.
fn looks_like_bgzf(buf: &[u8]) -> bool {
    buf.len() >= 18 && buf[0] == 0x1f && buf[1] == 0x8b && buf[3] & 0x04 != 0 && buf[12] == b'B' && buf[13] == b'C'
}

/// Open a plain or `.gz` FASTQ, choosing the fast path when it exists.
fn open(path: &Path, gz_threads: usize) -> Result<(Box<dyn BufRead + Send>, ByteCounter, &'static str)> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut head = [0u8; 18];
    let n = file
        .read(&mut head)
        .with_context(|| format!("reading header of {}", path.display()))?;
    let head = &head[..n];
    file.seek(SeekFrom::Start(0))
        .with_context(|| format!("rewinding {}", path.display()))?;

    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    let is_bgzf = is_gz && looks_like_bgzf(head);

    let counter = ByteCounter(Arc::new(AtomicU64::new(0)));
    let counting = Counting {
        inner: file,
        n: Arc::clone(&counter.0),
    };

    if is_bgzf && gz_threads > 1 {
        let dec = ParDecompressBuilder::<Bgzf>::new()
            .num_threads(gz_threads)
            .context("configuring parallel BGZF decompression")?
            .from_reader(counting);
        Ok((
            Box::new(BufReader::with_capacity(INPUT_BUFFER, dec)),
            counter,
            "bgzf-parallel",
        ))
    } else if is_gz {
        Ok((
            Box::new(BufReader::with_capacity(INPUT_BUFFER, MultiGzDecoder::new(counting))),
            counter,
            if is_bgzf { "bgzf-serial" } else { "gzip-serial" },
        ))
    } else {
        Ok((
            Box::new(BufReader::with_capacity(INPUT_BUFFER, counting)),
            counter,
            if gz_threads <= 1 { "plain" } else { "plain" },
        ))
    }
}

fn read_line(rdr: &mut dyn BufRead, out: &mut Vec<u8>) -> io::Result<usize> {
    out.clear();
    let n = rdr.read_until(b'\n', out)?;
    while matches!(out.last(), Some(b'\n') | Some(b'\r')) {
        out.pop();
    }
    Ok(n)
}

/// Read one file into `acc`, one pass, honouring the budgets.
pub fn scan_file(path: &Path, gz_threads: usize, limits: Limits, acc: &mut Accum) -> Result<ScanResult> {
    let ext_ok = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("fq") || e.eq_ignore_ascii_case("fastq"))
        .unwrap_or(true);
    if !ext_ok {
        bail!(
            "{}: fastcheck reads plain FASTQ and .gz FASTQ only — decompress the file first",
            path.display()
        );
    }

    let (mut rdr, counter, mode) = open(path, gz_threads)?;
    let (mut line, mut seq, mut plus, mut qual) = (
        Vec::with_capacity(512),
        Vec::with_capacity(512),
        Vec::with_capacity(512),
        Vec::with_capacity(512),
    );
    let mut reads: u64 = 0;
    let mut stopped = None;

    loop {
        if let Some(m) = limits.max_reads {
            if reads >= m {
                stopped = Some("reads budget");
                break;
            }
        }
        if let Some(b) = limits.max_bytes {
            if counter.0.load(Ordering::Relaxed) > b {
                stopped = Some("bytes budget");
                break;
            }
        }

        let n = read_line(&mut rdr, &mut line)
            .with_context(|| format!("reading {}", path.display()))?;
        if n == 0 {
            break; // clean EOF
        }
        if line.is_empty() {
            continue; // tolerate blank separators between records
        }
        if line[0] == b'>' {
            bail!(
                "{}: FASTA input is not supported by fastcheck (no quality string to check)",
                path.display()
            );
        }
        if line[0] != b'@' {
            bail!(
                "{}: expected a FASTQ record header starting with '@', found {:?}",
                path.display(),
                String::from_utf8_lossy(&line[..line.len().min(40)])
            );
        }

        read_line(&mut rdr, &mut seq)?;
        read_line(&mut rdr, &mut plus)?;
        read_line(&mut rdr, &mut qual)?;
        if seq.is_empty() {
            bail!(
                "{}: empty sequence line for record {:?}",
                path.display(),
                String::from_utf8_lossy(&line)
            );
        }
        if plus.first() != Some(&b'+') {
            bail!(
                "{}: FASTQ record {:?} is missing its '+' separator line",
                path.display(),
                String::from_utf8_lossy(&line)
            );
        }
        if qual.len() != seq.len() {
            bail!(
                "{}: record {:?} has seq length {} but quality length {} — multi-line/wrapped \
                 FASTQ is not supported, and guessing would produce wrong statistics",
                path.display(),
                String::from_utf8_lossy(&line),
                seq.len(),
                qual.len()
            );
        }

        acc.observe(&seq, &qual);
        reads += 1;
    }

    Ok(ScanResult {
        bytes_read: counter.0.load(Ordering::Relaxed),
        stopped,
        decompress: mode,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::accum::NT;
    use std::io::Write;

    fn write(dir: &Path, name: &str, body: &[u8]) -> std::path::PathBuf {
        let p = dir.join(name);
        std::fs::write(&p, body).unwrap();
        p
    }

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("fastcheck_reader_{}_{}", std::process::id(), tag));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    const FASTQ: &[u8] = b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGT\n+\n!!!!IIII\n";

    #[test]
    fn scans_plain_fastq() {
        let d = tmpdir("plain");
        let p = write(&d, "a.fastq", FASTQ);
        let mut acc = Accum::new(64);
        let r = scan_file(&p, 1, Limits::default(), &mut acc).unwrap();
        assert_eq!(acc.reads, 2);
        assert_eq!(acc.bases, 16);
        assert_eq!(r.stopped, None);
        assert_eq!(r.decompress, "plain");
    }

    #[test]
    fn reads_budget_stops_and_is_reported() {
        let d = tmpdir("budget");
        let p = write(&d, "a.fastq", FASTQ);
        let mut acc = Accum::new(64);
        let r = scan_file(
            &p,
            1,
            Limits {
                max_reads: Some(1),
                max_bytes: None,
            },
            &mut acc,
        )
        .unwrap();
        assert_eq!(acc.reads, 1, "exactly the budget, no more");
        assert_eq!(r.stopped, Some("reads budget"));
    }

    #[test]
    fn fasta_content_is_refused_with_a_clear_message() {
        // A file that lies about its format: .fastq extension, FASTA body.
        let d = tmpdir("fasta");
        let p = write(&d, "a.fastq", b">r1\nACGT\n");
        let mut acc = Accum::new(64);
        let err = scan_file(&p, 1, Limits::default(), &mut acc).unwrap_err();
        assert!(err.to_string().contains("FASTA"), "{err}");
    }

    #[test]
    fn unsupported_extension_is_refused() {
        let d = tmpdir("ext");
        let p = write(&d, "a.bam", b"not a fastq");
        let mut acc = Accum::new(64);
        let err = scan_file(&p, 1, Limits::default(), &mut acc).unwrap_err();
        assert!(err.to_string().contains("plain FASTQ"), "{err}");
    }

    #[test]
    fn wrapped_fastq_is_refused_not_guessed() {
        let d = tmpdir("wrapped");
        let p = write(&d, "a.fastq", b"@r1\nACGT\nACGT\n+\nIIII\nIIII\n");
        let mut acc = Accum::new(64);
        let err = scan_file(&p, 1, Limits::default(), &mut acc).unwrap_err();
        assert!(err.to_string().contains("wrapped"), "{err}");
    }

    #[test]
    fn plain_gzip_and_bgzf_yield_identical_counters() {
        let d = tmpdir("gz");
        let plain = write(&d, "a.fastq", FASTQ);

        // ordinary gzip
        let gz = d.join("a.fastq.gz");
        let mut enc = flate2::write::GzEncoder::new(File::create(&gz).unwrap(), flate2::Compression::default());
        enc.write_all(FASTQ).unwrap();
        enc.finish().unwrap();

        // BGZF, as fastprep/fastdedup write it
        let bgzf = d.join("b.fastq.gz");
        let f = File::create(&bgzf).unwrap();
        let mut w = gzp::par::compress::ParCompressBuilder::<Bgzf>::new()
            .num_threads(2)
            .unwrap()
            .from_writer(f);
        w.write_all(FASTQ).unwrap();
        drop(w);

        let mut a = Accum::new(64);
        let mut b = Accum::new(64);
        let mut c = Accum::new(64);
        let rp = scan_file(&plain, 1, Limits::default(), &mut a).unwrap();
        let rg = scan_file(&gz, 4, Limits::default(), &mut b).unwrap();
        let rb = scan_file(&bgzf, 4, Limits::default(), &mut c).unwrap();
        assert_eq!(rp.decompress, "plain");
        assert_eq!(rg.decompress, "gzip-serial");
        assert_eq!(rb.decompress, "bgzf-parallel");
        assert_eq!(a.reads, b.reads);
        assert_eq!(a.reads, c.reads);
        for p in 0..8 {
            for code in 0..NT {
                assert_eq!(a.base_at(p, code), b.base_at(p, code), "base {p}/{code}");
                assert_eq!(a.base_at(p, code), c.base_at(p, code), "base {p}/{code}");
            }
            assert_eq!(a.qual_at(p), b.qual_at(p), "qual at {p}");
            assert_eq!(a.qual_at(p), c.qual_at(p), "qual at {p}");
        }
        assert_eq!(a.qual_byte, b.qual_byte);
        assert_eq!(a.qual_byte, c.qual_byte);
    }
}
