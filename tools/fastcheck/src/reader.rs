//! Input: one pass, BGZF-parallel when possible, with explicit budgets.
//!
//! Supported inputs: plain FASTQ and `.gz` FASTQ. A `.gz` whose header says
//! BGZF (gzip FEXTRA + `BC` subfield — what Illumina writes and what
//! fastprep/fastdedup produce) is inflated block-parallel across
//! `gz_threads`. Any other `.gz` is inflated serially, which the report
//! records per file so a slow case is visible rather than mysterious. Other
//! extensions are refused next to the tool rather than half-supported.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use gzp::deflate::Bgzf;
use gzp::par::decompress::ParDecompressBuilder;

use crate::accum::Accum;
use crate::cli::Modules;
use crate::gated::Work;
use crate::overrep;

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
    /// Measurement work done per read. All zero on the ungated path, which does
    /// not count it: the counters exist for the gated run, where "did gating
    /// actually skip the work?" has to be answerable from outside.
    pub work: Work,
}

/// Counts bytes pulled from the underlying file, shared with the decoder so
/// the budget can be enforced mid-stream.
#[derive(Clone)]
pub(crate) struct ByteCounter(Arc<AtomicU64>);

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
    buf.len() >= 18
        && buf[0] == 0x1f
        && buf[1] == 0x8b
        && buf[3] & 0x04 != 0
        && buf[12] == b'B'
        && buf[13] == b'C'
}

/// gzip magic. Compression is decided from bytes, never from the filename, so
/// that a stream (`-`) and a mislabelled file behave the same way.
fn looks_like_gzip(buf: &[u8]) -> bool {
    buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
}

/// Read the first bytes of a stream, which cannot be rewound: they are chained
/// back in front of the reader (D-4, docs/DECISIONS.md).
fn peek<R: Read>(rdr: &mut R, n: usize) -> io::Result<Vec<u8>> {
    let mut head = vec![0u8; n];
    let mut have = 0;
    while have < n {
        let k = rdr.read(&mut head[have..])?;
        if k == 0 {
            break;
        }
        have += k;
    }
    head.truncate(have);
    Ok(head)
}

/// Open a plain or `.gz` FASTQ, from a file or from stdin (`-`).
///
/// `pub(crate)` so the insert-size pass reads through the same decompression
/// path (BGZF block-parallel included) instead of duplicating it.
pub(crate) fn open(
    path: &Path,
    gz_threads: usize,
) -> Result<(Box<dyn BufRead + Send>, ByteCounter, &'static str)> {
    let from_stdin = path.as_os_str() == "-";
    let mut source: Box<dyn Read + Send> = if from_stdin {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
    };

    let head = peek(&mut source, 18)?;
    let is_gz = looks_like_gzip(&head);
    let is_bgzf = looks_like_bgzf(&head);

    let counter = ByteCounter(Arc::new(AtomicU64::new(0)));
    // The peeked bytes go back in front; the counter sits *under* them so the
    // byte budget counts bytes that actually came from the file, not the peek.
    let counting = Counting {
        inner: io::Cursor::new(head).chain(source),
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
            if from_stdin { "stdin-bgzf-parallel" } else { "bgzf-parallel" },
        ))
    } else if is_gz {
        Ok((
            Box::new(BufReader::with_capacity(
                INPUT_BUFFER,
                MultiGzDecoder::new(counting),
            )),
            counter,
            if is_bgzf {
                "bgzf-serial"
            } else if from_stdin {
                "stdin-gzip-serial"
            } else {
                "gzip-serial"
            },
        ))
    } else {
        Ok((
            Box::new(BufReader::with_capacity(INPUT_BUFFER, counting)),
            counter,
            if from_stdin { "stdin-plain" } else { "plain" },
        ))
    }
}

pub(crate) fn read_line(rdr: &mut dyn BufRead, out: &mut Vec<u8>) -> io::Result<usize> {
    out.clear();
    let n = rdr.read_until(b'\n', out)?;
    while matches!(out.last(), Some(b'\n') | Some(b'\r')) {
        out.pop();
    }
    Ok(n)
}

/// Read one file into `acc`, one pass, honouring the budgets, measuring the modules
/// `gate` enables.
///
/// `None` is the ungated path — every per-read measurement runs (`Accum::observe`),
/// which is what any run without a gating flag takes, so its numbers are exactly what
/// they were before `--only` existed. `Some(modules)` runs `Accum::observe_gated`: the
/// quality line is still read and length-checked (it is part of the record's framing,
/// and a short quality line is an error, not a thing to guess at), but it is never
/// decoded, no position row is grown for a module that is off, and the grid, the polyX
/// runs, the fingerprint table and the k-mer counter are only touched when their module
/// is on.
pub fn scan_file_gated(
    path: &Path,
    gz_threads: usize,
    limits: Limits,
    acc: &mut Accum,
    gate: Option<Modules>,
) -> Result<ScanResult> {
    let from_stdin = path.as_os_str() == "-";
    // A named file must have a plausible FASTQ extension, so a .bam/.bz2/.xz
    // gets a clear message instead of being parsed as text. A stream has no
    // name, so it is taken at face value.
    let ext_ok = from_stdin
        || path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| {
                e.eq_ignore_ascii_case("gz")
                    || e.eq_ignore_ascii_case("fq")
                    || e.eq_ignore_ascii_case("fastq")
            })
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
    let mut work = Work::default();

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

        // Every record line carries the same "reading <path>" context as the header: a
        // truncated .gz used to report a bare "unexpected end of file" whenever the
        // truncation fell on the sequence, '+' or quality line, naming no file at all.
        read_line(&mut rdr, &mut seq)
            .with_context(|| format!("reading {}", path.display()))?;
        read_line(&mut rdr, &mut plus)
            .with_context(|| format!("reading {}", path.display()))?;
        read_line(&mut rdr, &mut qual)
            .with_context(|| format!("reading {}", path.display()))?;
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

        // --overrep-mode table packs each sequence 2 bits per base, so a read with an N or
        // an IUPAC code cannot be represented. Dropping it (what `let _ = observe(..)` used
        // to do) under-counts the unique-sequence table exactly on the libraries the mode
        // exists to name, so the mode refuses instead, naming the record and the offending
        // base's position. The other over-representation modes do not pack and are unaffected.
        if acc.track_overrep_table {
            if let Some(pos) = overrep::first_unpackable(&seq) {
                bail!(
                    "{}: record {:?} has {:?} at sequence position {}, which --overrep-mode \
                     table cannot pack (2-bit ACGT only). Refusing rather than silently \
                     dropping it from the unique-sequence table; remove reads with N or IUPAC \
                     codes, or use --overrep-mode fingerprint or second-pass.",
                    path.display(),
                    String::from_utf8_lossy(&line),
                    seq[pos] as char,
                    pos
                );
            }
        }

        // The record is measured here: gated when a module subset was selected, the
        // ungated (unchanged) path otherwise.
        let w = match &gate {
            Some(m) => acc.observe_gated(&line, &seq, &qual, m),
            None => {
                acc.observe(&line, &seq, &qual);
                Work::default()
            }
        };
        work.add(&w);
        reads += 1;
    }

    Ok(ScanResult {
        bytes_read: counter.0.load(Ordering::Relaxed),
        stopped,
        decompress: mode,
        work,
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
        let r = scan_file_gated(&p, 1, Limits::default(), &mut acc, None).unwrap();
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
        let r = scan_file_gated(
            &p,
            1,
            Limits {
                max_reads: Some(1),
                max_bytes: None,
            },
            &mut acc,
            None,
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
        let err = scan_file_gated(&p, 1, Limits::default(), &mut acc, None).unwrap_err();
        assert!(err.to_string().contains("FASTA"), "{err}");
    }

    #[test]
    fn unsupported_extension_is_refused() {
        let d = tmpdir("ext");
        let p = write(&d, "a.bam", b"not a fastq");
        let mut acc = Accum::new(64);
        let err = scan_file_gated(&p, 1, Limits::default(), &mut acc, None).unwrap_err();
        assert!(err.to_string().contains("plain FASTQ"), "{err}");
    }

    #[test]
    fn wrapped_fastq_is_refused_not_guessed() {
        let d = tmpdir("wrapped");
        let p = write(&d, "a.fastq", b"@r1\nACGT\nACGT\n+\nIIII\nIIII\n");
        let mut acc = Accum::new(64);
        let err = scan_file_gated(&p, 1, Limits::default(), &mut acc, None).unwrap_err();
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
        let rp = scan_file_gated(&plain, 1, Limits::default(), &mut a, None).unwrap();
        let rg = scan_file_gated(&gz, 4, Limits::default(), &mut b, None).unwrap();
        let rb = scan_file_gated(&bgzf, 4, Limits::default(), &mut c, None).unwrap();
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

    /// The gated scan measures what its modules say and nothing else, and the ungated
    /// scan counts no work at all (it does not keep counters).
    #[test]
    fn gated_scan_measures_only_the_selected_modules() {
        let d = tmpdir("gated");
        let p = write(&d, "a.fastq", FASTQ);

        let mut acc = Accum::new(64);
        let m = Modules { dup: true, ..Modules::none() };
        let r = scan_file_gated(&p, 1, Limits::default(), &mut acc, Some(m)).unwrap();
        assert_eq!(acc.reads, 2, "reads are always counted");
        assert_eq!(acc.bases, 16);
        assert_eq!(r.work.fingerprints, 2, "dup fingerprints every read");
        assert_eq!(r.work.base_cells, 0, "no composition row was grown");
        assert_eq!(r.work.quality_bytes, 0, "no quality byte was decoded");
        assert_eq!(acc.base_pos.len(), 0, "bases off: no row allocated");
        assert_eq!(acc.qual_pos.len(), 0, "quality off: no row allocated");
        assert_eq!(acc.len_hist.len(), 0, "length off: no histogram");
        assert_eq!(acc.reads_with_n, 0, "the fixture has no N");

        let mut u = Accum::new(64);
        let ru = scan_file_gated(&p, 1, Limits::default(), &mut u, None).unwrap();
        assert_eq!(ru.work, Work::default(), "the ungated path keeps no counters");
        assert_eq!(u.base_pos.len(), 8, "and measures everything, as it always did");
        assert_eq!(u.qual_pos.len(), 8);
        assert_eq!(u.len_hist[8], 2);

        // The quality line is still part of the framing: a gated run cannot skip it.
        let bad = write(&d, "bad.fastq", b"@r1\nACGTACGT\n+\nIIII\n");
        let mut acc = Accum::new(64);
        let err = scan_file_gated(&bad, 1, Limits::default(), &mut acc, Some(m)).unwrap_err();
        assert!(err.to_string().contains("quality length"), "{err}");
    }
}

#[cfg(test)]
mod regression_tests {
    use super::*;
    use std::io::Write;

    fn gz(bytes: &[u8]) -> Vec<u8> {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    /// C11: every record-line read carries the "reading <path>" context, not only the
    /// header. A truncated .gz used to report a bare "unexpected end of file" whenever the
    /// cut fell inside a record, naming no file at all.
    ///
    /// The fixture makes WHICH line fails deterministic: the first gzip member ends after
    /// 0, 1, 2 and 3 lines of a record (so the header, sequence, '+' and quality read
    /// respectively is the next one), and a truncated second member makes that read fail.
    #[test]
    fn a_truncated_gzip_names_the_file_on_every_record_line() {
        let d = std::env::temp_dir().join(format!("fastcheck_reader_trunc_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let path = d.join("t.fastq.gz");

        let cut_tail = &gz(b"@r2\nACGTACGT\n+\nIIIIIIII\n")[..12];
        let prefixes: [&[u8]; 4] = [
            b"",                       // the header read fails
            b"@r1\n",                  // the sequence read fails
            b"@r1\nACGTACGT\n",        // the '+' read fails
            b"@r1\nACGTACGT\n+\n",    // the quality read fails
        ];
        for (i, prefix) in prefixes.iter().enumerate() {
            let mut data = gz(prefix);
            data.extend_from_slice(cut_tail);
            std::fs::write(&path, &data).unwrap();
            let mut acc = Accum::new(64);
            let err = scan_file_gated(&path, 1, Limits::default(), &mut acc, None)
                .expect_err("a truncated member must be an error, not a short read");
            let msg = format!("{err:#}");
            assert!(
                msg.contains("reading") && msg.contains("t.fastq.gz"),
                "case {i}: the error must name the file it was reading: {msg:?}"
            );
        }
        let _ = std::fs::remove_dir_all(&d);
    }
}
