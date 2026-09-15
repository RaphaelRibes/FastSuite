//! Processing pipeline: reads → chunks → parallel trim → ordered write.
//!
//! # Ownership
//!
//! needletail hands back records as views into its internal buffer, which
//! is invalidated on the next `next()` call. To fan out to a rayon pool
//! we need `'static` data, so we copy seq + qual. The id is also copied
//! because we need it to survive past the decision for any kept reads;
//! the fields are small enough (~typically 50–200 bytes per read) that
//! the alloc cost is overshadowed by parsing + compression on any real
//! workload.
//!
//! # Pair alignment
//!
//! `chunk1.into_par_iter().zip(chunk2.into_par_iter())` in rayon builds
//! a joined parallel iterator that splits its combined index range
//! across tasks. Each task processes a contiguous slice `[i..j)` and
//! R1[i] is always paired with R2[i]. There is no reordering within a
//! chunk and no drift between chunks, regardless of thread count.

use anyhow::{bail, Context, Result};
use crossbeam_channel::bounded;
use needletail::{parse_fastx_reader, parser::SequenceRecord, FastxReader};
use rayon::prelude::*;
use std::io::IsTerminal;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::adapter::AdapterHint;
use crate::algorithm::{trim, DropReason, Params, TrimOutput, TrimResult};
use crate::cli::Runtime;
use crate::dedup::{AnchorSpec, Fingerprint, Seen, Verifier};
use crate::io::{
    open_input, open_output, open_output_as, write_fasta, write_fastq, OutputWriter, RecordFormat,
};
use crate::metadata::{self, Mate};
use crate::stats::Stats;

// ---------------------------------------------------------------------------
// Record types.
// ---------------------------------------------------------------------------

/// Chunk-owned byte storage (§10.4 zero-copy chunks). All reads of one
/// chunk share three contiguous buffers (id, seq, qual); each
/// [`ScanInput`] is three ranges into them instead of three heap
/// allocations per read. At 300 M reads this removes ~900 M allocator
/// calls (~60 GB of allocation traffic) from the parse path, which is
/// serial and therefore on the critical path for every stage.
#[derive(Default)]
struct ChunkBuf {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl ChunkBuf {
    fn reserve(&mut self, reads: usize, avg_len: usize) {
        let cap = reads * avg_len;
        self.id.reserve(cap / 4);
        self.seq.reserve(cap);
        self.qual.reserve(cap);
    }

    fn clear(&mut self) {
        self.id.clear();
        self.seq.clear();
        self.qual.clear();
    }
}

/// Owned copy of one record, suitable for moving across thread boundaries.
///
/// With the chunk arena the three byte slices live in the chunk's
/// [`ChunkBuf`]; the record itself is three `Range`s and a format tag.
struct ScanInput {
    id: std::ops::Range<usize>,
    seq: std::ops::Range<usize>,
    qual: std::ops::Range<usize>,
    format: RecordFormat,
}

impl ScanInput {
    fn push_needletail(buf: &mut ChunkBuf, rec: &SequenceRecord, phred_offset: u8) -> Result<Self> {
        let id_start = buf.id.len();
        buf.id.extend_from_slice(rec.id());
        let seq_start = buf.seq.len();
        buf.seq.extend_from_slice(rec.seq().as_ref());
        let (qual, format) = match rec.qual() {
            Some(q) => {
                // Decode ASCII Phred. saturating_sub handles corrupted
                // bytes below the offset (would otherwise underflow).
                let q_start = buf.qual.len();
                buf.qual
                    .extend(q.iter().map(|&b| b.saturating_sub(phred_offset)));
                (q_start..buf.qual.len(), RecordFormat::Fastq)
            }
            None => {
                let q_start = buf.qual.len();
                buf.qual.resize(q_start + buf.seq.len() - seq_start, 0u8);
                (q_start..buf.qual.len(), RecordFormat::Fasta)
            }
        };
        Ok(Self {
            id: id_start..buf.id.len(),
            seq: seq_start..buf.seq.len(),
            qual,
            format,
        })
    }
}

struct Processed {
    input: ScanInput,
    output: TrimOutput,
}

impl Processed {
    fn emit_to(&self, buf: &ChunkBuf, w: &mut OutputWriter, phred_offset: u8) -> Result<()> {
        match &self.output.result {
            TrimResult::Drop { .. } => Ok(()),
            TrimResult::Emit { .. } => {
                let t_5 = self.output.result.t5();
                let t_3 = self.output.result.t3();
                let seq = &buf.seq[self.input.seq.clone()];
                debug_assert!(
                    t_5 + t_3 <= seq.len(),
                    "trim exceeds read length: t_5={} t_3={} L={}",
                    t_5,
                    t_3,
                    seq.len()
                );
                let seq = &seq[t_5..seq.len() - t_3];
                match self.input.format {
                    RecordFormat::Fastq => {
                        let qual = &buf.qual[self.input.qual.clone()];
                        let raw_qual = &qual[t_5..qual.len() - t_3];
                        write_fastq(
                            w,
                            &buf.id[self.input.id.clone()],
                            seq,
                            raw_qual,
                            phred_offset,
                        )
                    }
                    RecordFormat::Fasta => write_fasta(w, &buf.id[self.input.id.clone()], seq),
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Chunk filling.
// ---------------------------------------------------------------------------

fn fill_chunk(
    reader: &mut dyn FastxReader,
    buf: &mut ChunkBuf,
    chunk: &mut Vec<ScanInput>,
    target: usize,
    phred_offset: u8,
) -> Result<usize> {
    buf.clear();
    chunk.clear();
    while chunk.len() < target {
        match reader.next() {
            None => break,
            Some(Err(e)) => return Err(anyhow::Error::new(e)).context("parsing FASTX"),
            Some(Ok(rec)) => chunk.push(ScanInput::push_needletail(buf, &rec, phred_offset)?),
        }
    }
    Ok(chunk.len())
}

fn process_chunk(
    buf: &ChunkBuf,
    chunk: Vec<ScanInput>,
    params: &Params,
    parallel: bool,
) -> Vec<Processed> {
    let work = |input: ScanInput| -> Processed {
        let seq = &buf.seq[input.seq.clone()];
        let qual = &buf.qual[input.qual.clone()];
        let output = trim(seq, qual, params, AdapterHint::Detect);
        Processed { input, output }
    };
    if parallel {
        chunk.into_par_iter().map(work).collect()
    } else {
        chunk.into_iter().map(work).collect()
    }
}

// ---------------------------------------------------------------------------
// Output safety helpers (tmp + atomic rename).
// ---------------------------------------------------------------------------

/// Compression used for metadata files regardless of FASTQ output level.
/// Metadata is a small random-access binary format; heavy compression is
/// counterproductive and the Python analyzer decompresses it entirely anyway.
const METADATA_COMPRESSION: crate::io::CompressionLevels = crate::io::CompressionLevels {
    gz: 1,
    zst: 1,
    // Metadata files are small; a compression pool would cost more in thread
    // spawn than it saves in deflate.
    gz_threads: 1,
};

/// Derive a `.part` write-side path for `p`. On success the caller renames
/// this path to `p`; on error the guard removes it.
/// Fingerprint set, monomorphised over the two widths.
///
/// The fingerprint is computed here on the ordered thread rather than in the
/// parallel map. xxh3 over a ~90-base window is roughly 20 ns while the set
/// lookup at multi-gigabyte table sizes is a DRAM round trip at ~80 ns, so
/// moving the hash into the workers would save a quarter of an already small
/// serial cost and would require threading a fingerprint field through
/// `Processed` and `ProcessedPair`. Not worth it until a profile says so.
enum Dedup {
    Bits64(Verifier<u64>, AnchorSpec),
    Bits128(Verifier<u128>, AnchorSpec),
}

impl Dedup {
    fn new(cfg: &crate::cli::DedupConfig, spec: AnchorSpec) -> Self {
        match cfg.bits {
            128 => Dedup::Bits128(Verifier::with_capacity(cfg.capacity, cfg.ceiling), spec),
            _ => Dedup::Bits64(Verifier::with_capacity(cfg.capacity, cfg.ceiling), spec),
        }
    }

    /// `true` when the read (or pair) is a duplicate and should be dropped.
    fn is_duplicate(&mut self, r1: &[u8], r2: Option<&[u8]>) -> bool {
        match self {
            Dedup::Bits64(v, spec) => {
                let fp = match r2 {
                    Some(b) => u64::of_pair(spec.slice(r1), spec.slice(b)),
                    None => u64::of(spec.slice(r1)),
                };
                v.offer(fp) == Seen::Duplicate
            }
            Dedup::Bits128(v, spec) => {
                let fp = match r2 {
                    Some(b) => u128::of_pair(spec.slice(r1), spec.slice(b)),
                    None => u128::of(spec.slice(r1)),
                };
                v.offer(fp) == Seen::Duplicate
            }
        }
    }

    fn saturated_at(&self) -> Option<u64> {
        match self {
            Dedup::Bits64(v, _) => v.saturated_at(),
            Dedup::Bits128(v, _) => v.saturated_at(),
        }
    }
}

/// Replaces a trim result with a duplicate drop. The read was never scanned,
/// so every trim count is zero.
fn duplicate_result() -> TrimResult {
    TrimResult::Drop {
        reason: DropReason::Duplicate,
        t5_adapter: 0,
        t3_adapter: 0,
        t5_poly: 0,
        t3_poly: 0,
        runs: vec![],
    }
}

fn part_path(p: &Path) -> std::path::PathBuf {
    let mut s = p.as_os_str().to_owned();
    s.push(".part");
    std::path::PathBuf::from(s)
}

/// Removes its path on drop unless `commit()` has been called. Ensures
/// partial output files are cleaned up on any error path.
struct TmpFile {
    path: std::path::PathBuf,
    committed: bool,
}

impl TmpFile {
    fn new(path: std::path::PathBuf) -> Self {
        Self {
            path,
            committed: false,
        }
    }

    /// Mark as successfully committed; `Drop` will no longer delete the file.
    fn commit(mut self) {
        self.committed = true;
    }
}

impl Drop for TmpFile {
    fn drop(&mut self) {
        if !self.committed {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

// ---------------------------------------------------------------------------
// Phred offset auto-detection.
// ---------------------------------------------------------------------------

/// Peek at the first `n_reads` records (up to 1000) of a FASTQ file and
/// infer the Phred offset from the raw quality bytes.
///
/// Heuristic: any byte `< '@'` (64) cannot appear in Phred+64 (which
/// starts at `'@'`) → must be Phred+33. If all observed bytes are ≥ 64,
/// fall back to Phred+64. FASTA files (no quality) default to 33.
pub fn detect_phred_offset(path: &Path) -> Result<u8> {
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("auto-detecting Phred offset from {}", path.display()))?;
    let mut min_q = u8::MAX;
    let mut n = 0usize;
    while n < 1000 {
        match reader.next() {
            None => break,
            Some(Err(e)) => return Err(anyhow::Error::new(e)),
            Some(Ok(rec)) => {
                if let Some(q) = rec.qual() {
                    for &b in q {
                        if b < min_q {
                            min_q = b;
                        }
                    }
                }
                n += 1;
            }
        }
    }
    if n == 0 || min_q == u8::MAX {
        return Ok(33); // FASTA or empty file: default
    }
    Ok(if min_q < 64 { 33 } else { 64 })
}

// ---------------------------------------------------------------------------
// Single-end.
// ---------------------------------------------------------------------------

pub fn run_single(in1: &Path, out1: &Path, params: &Params, rt: &Runtime) -> Result<Stats> {
    // Phred offset: auto-detect or use the configured value.
    let phred_offset = if rt.phred_auto {
        detect_phred_offset(in1)?
    } else {
        rt.phred_offset
    };

    let mut reader = parse_fastx_reader(open_input(in1, rt.compression.gz_threads)?)
        .with_context(|| format!("opening input {}", in1.display()))?;

    // Open all outputs to `.part` paths; rename atomically on success.
    let tmp_out = part_path(out1);
    let guard_out = TmpFile::new(tmp_out.clone());
    let mut writer = open_output_as(&tmp_out, out1, rt.compression)?;

    let tmp_meta_path = rt.metadata.as_deref().map(part_path);
    let guard_meta = tmp_meta_path.as_ref().map(|p| TmpFile::new(p.clone()));
    let mut meta = match (rt.metadata.as_deref(), tmp_meta_path.as_deref()) {
        (Some(_), Some(tp)) => {
            let mut w = open_output(tp, METADATA_COMPRESSION)?;
            metadata::write_file_header(&mut *w).context("writing metadata header")?;
            Some(w)
        }
        _ => None,
    };

    let mut stats = Stats::default();
    let parallel = rt.threads > 1;

    let params = params.clone();

    // The anchor window is resolved once for the run from the expected read
    // length; deriving it per read would make it vary with the polyG tail,
    // which is exactly the divergence it exists to absorb.
    let mut dedup = rt.dedup.as_ref().map(|cfg| {
        let spec = AnchorSpec::new(&params, cfg.read_len, params.adapter.is_some());
        Dedup::new(cfg, spec)
    });

    let (raw_tx, raw_rx) = bounded::<(ChunkBuf, Vec<ScanInput>)>(2);
    let (res_tx, res_rx) = bounded::<(ChunkBuf, Vec<Processed>)>(2);

    let params_arc = Arc::new(params);
    let thread_pool = rt.thread_pool.clone();
    let compute = std::thread::spawn(move || {
        while let Ok((buf, chunk)) = raw_rx.recv() {
            let results = match thread_pool.as_ref() {
                Some(pool) => pool.install(|| process_chunk(&buf, chunk, &params_arc, parallel)),
                None => process_chunk(&buf, chunk, &params_arc, parallel),
            };
            if res_tx.send((buf, results)).is_err() {
                break;
            }
        }
    });

    let stderr_is_tty = std::io::stderr().is_terminal();
    let mut chunks_written: u64 = 0;
    let mut progress_shown = false;

    let io_result = (|| -> Result<()> {
        let mut buf = ChunkBuf::default();
        buf.reserve(rt.chunk, 151);
        let mut chunk = Vec::with_capacity(rt.chunk);
        fill_chunk(&mut *reader, &mut buf, &mut chunk, rt.chunk, phred_offset)?;

        while !chunk.is_empty() {
            raw_tx.send((buf, chunk)).ok();

            let mut next_buf = ChunkBuf::default();
            next_buf.reserve(rt.chunk, 151);
            let mut next = Vec::with_capacity(rt.chunk);
            fill_chunk(
                &mut *reader,
                &mut next_buf,
                &mut next,
                rt.chunk,
                phred_offset,
            )?;

            let (res_buf, results) = res_rx.recv().context("compute thread died")?;
            for p in &results {
                // Deduplication runs on the ordered thread, so the decision
                // sequence is the input order regardless of thread count.
                let seq = &res_buf.seq[p.input.seq.clone()];
                let dup = dedup.as_mut().is_some_and(|d| d.is_duplicate(seq, None));
                let result = if dup {
                    duplicate_result()
                } else {
                    p.output.result.clone()
                };
                stats.record(seq.len(), &result, p.output.adapter_name.as_deref());
                if !dup {
                    p.emit_to(&res_buf, &mut writer, phred_offset)?;
                }
                if let Some(mw) = meta.as_mut() {
                    metadata::write_record(
                        &mut **mw,
                        &res_buf.id[p.input.id.clone()],
                        Mate::R1,
                        &result,
                        seq.len(),
                    )
                    .context("writing metadata record")?;
                }
            }
            chunks_written += 1;
            if stderr_is_tty && chunks_written % 100 == 0 {
                eprint!("\r  {} reads processed", stats.reads_in);
                progress_shown = true;
            }
            buf = next_buf;
            chunk = next;
        }
        Ok(())
    })();

    // Drop raw_tx so the compute thread sees the channel closed and exits.
    drop(raw_tx);
    // Drain any results the compute thread produced after our last recv.
    for (res_buf, results) in res_rx {
        for p in &results {
            let seq = &res_buf.seq[p.input.seq.clone()];
            stats.record(
                seq.len(),
                &p.output.result,
                p.output.adapter_name.as_deref(),
            );
            if io_result.is_ok() {
                let _ = p.emit_to(&res_buf, &mut writer, phred_offset);
            }
        }
    }
    compute.join().ok();

    if let Some(at) = dedup.as_ref().and_then(|d| d.saturated_at()) {
        eprintln!(
            "fastprep: warning: deduplication saturated at read {at}; fingerprints \
             seen after that point were not recorded, so the reported \
             duplication rate is a lower bound"
        );
    }

    if progress_shown {
        eprint!("\x1b[2K\r"); // ANSI: clear line + return to start
    }

    // Flush and drop writers before rename (zstd auto_finish finalises the
    // frame on drop; rename must happen after the frame is closed).
    writer.flush().ok();
    drop(writer);
    if let Some(mut mw) = meta.take() {
        mw.flush().ok();
    }

    io_result?;

    // Atomic rename: committed output replaces the destination path.
    std::fs::rename(&tmp_out, out1)
        .with_context(|| format!("committing {} → {}", tmp_out.display(), out1.display()))?;
    guard_out.commit();
    if let (Some(tp), Some(fp), Some(g)) = (tmp_meta_path, rt.metadata.as_deref(), guard_meta) {
        std::fs::rename(&tp, fp)
            .with_context(|| format!("committing {} → {}", tp.display(), fp.display()))?;
        g.commit();
    }

    Ok(stats)
}

// ---------------------------------------------------------------------------
// Paired-end.
// ---------------------------------------------------------------------------

#[derive(Debug, serde::Serialize)]
pub struct PairedStats {
    pub r1: Stats,
    pub r2: Stats,
    pub pairs_out: u64,
    pub pairs_dropped: u64,
    pub singletons_r1: u64,
    pub singletons_r2: u64,
    /// Pairs where one mate survived but was discarded because the
    /// corresponding --singletons* path was not provided.
    pub discarded_singletons_r1: u64,
    pub discarded_singletons_r2: u64,
}

enum PairOutcome {
    BothKept,
    BothDropped,
    R1KeptR2Dropped,
    R1DroppedR2Kept,
}

struct ProcessedPair {
    p1: Processed,
    p2: Processed,
}

pub fn run_paired(
    in1: &Path,
    in2: &Path,
    out1: &Path,
    out2: &Path,
    singletons1: Option<&Path>,
    singletons2: Option<&Path>,
    params_r1: &Params,
    params_r2: &Params,
    rt: &Runtime,
) -> Result<PairedStats> {
    // Phred offset: auto-detect from R1 (R2 must match).
    let phred_offset = if rt.phred_auto {
        detect_phred_offset(in1)?
    } else {
        rt.phred_offset
    };

    let mut r1 = parse_fastx_reader(open_input(in1, rt.compression.gz_threads)?)
        .with_context(|| format!("opening input {}", in1.display()))?;
    let mut r2 = parse_fastx_reader(open_input(in2, rt.compression.gz_threads)?)
        .with_context(|| format!("opening input {}", in2.display()))?;

    // Open all outputs to `.part` paths; rename atomically on success.
    let tmp_w1 = part_path(out1);
    let tmp_w2 = part_path(out2);
    let guard_w1 = TmpFile::new(tmp_w1.clone());
    let guard_w2 = TmpFile::new(tmp_w2.clone());
    let mut w1 = open_output_as(&tmp_w1, out1, rt.compression)?;
    let mut w2 = open_output_as(&tmp_w2, out2, rt.compression)?;

    let tmp_ws1 = singletons1.map(part_path);
    let tmp_ws2 = singletons2.map(part_path);
    let guard_ws1 = tmp_ws1.as_ref().map(|p| TmpFile::new(p.clone()));
    let guard_ws2 = tmp_ws2.as_ref().map(|p| TmpFile::new(p.clone()));
    let mut ws1 = tmp_ws1
        .as_deref()
        .zip(singletons1)
        .map(|(tp, fp)| open_output_as(tp, fp, rt.compression))
        .transpose()?;
    let mut ws2 = tmp_ws2
        .as_deref()
        .zip(singletons2)
        .map(|(tp, fp)| open_output_as(tp, fp, rt.compression))
        .transpose()?;

    let tmp_meta_path = rt.metadata.as_deref().map(part_path);
    let guard_meta = tmp_meta_path.as_ref().map(|p| TmpFile::new(p.clone()));
    let mut meta = match (rt.metadata.as_deref(), tmp_meta_path.as_deref()) {
        (Some(_), Some(tp)) => {
            let mut w = open_output(tp, METADATA_COMPRESSION)?;
            metadata::write_file_header(&mut *w).context("writing metadata header")?;
            Some(w)
        }
        _ => None,
    };

    let warn_r1_once = AtomicBool::new(false);
    let warn_r2_once = AtomicBool::new(false);

    let mut stats = PairedStats {
        r1: Stats::default(),
        r2: Stats::default(),
        pairs_out: 0,
        pairs_dropped: 0,
        singletons_r1: 0,
        singletons_r2: 0,
        discarded_singletons_r1: 0,
        discarded_singletons_r2: 0,
    };

    let parallel = rt.threads > 1;
    let params_r1 = params_r1.clone();
    let params_r2 = params_r2.clone();

    let mut dedup = rt.dedup.as_ref().map(|cfg| {
        let spec = AnchorSpec::new(&params_r1, cfg.read_len, params_r1.adapter.is_some());
        Dedup::new(cfg, spec)
    });

    type ChunkPair = ((ChunkBuf, Vec<ScanInput>), (ChunkBuf, Vec<ScanInput>));
    let (raw_tx, raw_rx) = bounded::<ChunkPair>(2);
    let (res_tx, res_rx) = bounded::<((ChunkBuf, ChunkBuf), Vec<ProcessedPair>)>(2);

    let params_r1_arc = Arc::new(params_r1);
    let params_r2_arc = Arc::new(params_r2);
    let thread_pool = rt.thread_pool.clone();
    let compute = std::thread::spawn(move || {
        while let Ok(((b1, c1), (b2, c2))) = raw_rx.recv() {
            let results = match thread_pool.as_ref() {
                Some(pool) => pool.install(|| {
                    pair_chunks(&b1, &b2, c1, c2, &params_r1_arc, &params_r2_arc, parallel)
                }),
                None => pair_chunks(&b1, &b2, c1, c2, &params_r1_arc, &params_r2_arc, parallel),
            };
            if res_tx.send(((b1, b2), results)).is_err() {
                break;
            }
        }
    });

    let stderr_is_tty = std::io::stderr().is_terminal();
    let mut chunks_written: u64 = 0;
    let mut progress_shown = false;

    let io_result = (|| -> Result<()> {
        let mut buf1 = ChunkBuf::default();
        buf1.reserve(rt.chunk, 151);
        let mut buf2 = ChunkBuf::default();
        buf2.reserve(rt.chunk, 151);
        let mut chunk1 = Vec::with_capacity(rt.chunk);
        let mut chunk2 = Vec::with_capacity(rt.chunk);
        fill_chunk(&mut *r1, &mut buf1, &mut chunk1, rt.chunk, phred_offset)?;
        fill_chunk(&mut *r2, &mut buf2, &mut chunk2, rt.chunk, phred_offset)?;
        let (n1, n2) = (chunk1.len(), chunk2.len());
        if n1 != n2 && !(n1 == 0 && n2 == 0) {
            bail!(
                "paired input de-synced on first chunk: R1 had {} records, R2 had {}.",
                n1,
                n2
            );
        }

        while !chunk1.is_empty() {
            raw_tx.send(((buf1, chunk1), (buf2, chunk2))).ok();

            let mut next_buf1 = ChunkBuf::default();
            next_buf1.reserve(rt.chunk, 151);
            let mut next_buf2 = ChunkBuf::default();
            next_buf2.reserve(rt.chunk, 151);
            let mut next1 = Vec::with_capacity(rt.chunk);
            let mut next2 = Vec::with_capacity(rt.chunk);
            let n1 = fill_chunk(&mut *r1, &mut next_buf1, &mut next1, rt.chunk, phred_offset)?;
            let n2 = fill_chunk(&mut *r2, &mut next_buf2, &mut next2, rt.chunk, phred_offset)?;
            if n1 != n2 {
                bail!(
                    "paired input de-synced: R1 chunk had {} records, R2 had {}. \
                     Are R1 and R2 the same length?",
                    n1,
                    n2
                );
            }

            let ((res_b1, res_b2), pairs) = res_rx.recv().context("compute thread died")?;
            write_pairs(
                &res_b1,
                &res_b2,
                pairs,
                &mut stats,
                &mut w1,
                &mut w2,
                ws1.as_mut(),
                ws2.as_mut(),
                meta.as_mut(),
                phred_offset,
                &warn_r1_once,
                &warn_r2_once,
                dedup.as_mut(),
            )?;

            chunks_written += 1;
            if stderr_is_tty && chunks_written % 100 == 0 {
                eprint!(
                    "\r  {} reads processed",
                    stats.r1.reads_in + stats.r2.reads_in
                );
                progress_shown = true;
            }
            buf1 = next_buf1;
            buf2 = next_buf2;
            chunk1 = next1;
            chunk2 = next2;
        }
        Ok(())
    })();

    drop(raw_tx);
    for ((res_b1, res_b2), pairs) in res_rx {
        let _ = write_pairs(
            &res_b1,
            &res_b2,
            pairs,
            &mut stats,
            &mut w1,
            &mut w2,
            ws1.as_mut(),
            ws2.as_mut(),
            meta.as_mut(),
            phred_offset,
            &warn_r1_once,
            &warn_r2_once,
            dedup.as_mut(),
        );
    }
    compute.join().ok();

    if let Some(at) = dedup.as_ref().and_then(|d| d.saturated_at()) {
        eprintln!(
            "fastprep: warning: deduplication saturated at read {at}; fingerprints \
             seen after that point were not recorded, so the reported \
             duplication rate is a lower bound"
        );
    }

    if progress_shown {
        eprint!("\x1b[2K\r");
    }

    w1.flush().ok();
    drop(w1);
    w2.flush().ok();
    drop(w2);
    if let Some(mut w) = ws1.take() {
        w.flush().ok();
    }
    if let Some(mut w) = ws2.take() {
        w.flush().ok();
    }
    if let Some(mut mw) = meta.take() {
        mw.flush().ok();
    }

    io_result?;

    // Atomic renames.
    std::fs::rename(&tmp_w1, out1)
        .with_context(|| format!("committing {} → {}", tmp_w1.display(), out1.display()))?;
    guard_w1.commit();
    std::fs::rename(&tmp_w2, out2)
        .with_context(|| format!("committing {} → {}", tmp_w2.display(), out2.display()))?;
    guard_w2.commit();
    if let (Some(tp), Some(fp), Some(g)) = (tmp_ws1, singletons1, guard_ws1) {
        std::fs::rename(&tp, fp)
            .with_context(|| format!("committing {} → {}", tp.display(), fp.display()))?;
        g.commit();
    }
    if let (Some(tp), Some(fp), Some(g)) = (tmp_ws2, singletons2, guard_ws2) {
        std::fs::rename(&tp, fp)
            .with_context(|| format!("committing {} → {}", tp.display(), fp.display()))?;
        g.commit();
    }
    if let (Some(tp), Some(fp), Some(g)) = (tmp_meta_path, rt.metadata.as_deref(), guard_meta) {
        std::fs::rename(&tp, fp)
            .with_context(|| format!("committing {} → {}", tp.display(), fp.display()))?;
        g.commit();
    }

    Ok(stats)
}

/// Write one batch of processed pairs to all output streams.
/// Extracted to avoid duplicating the logic between the main loop and the drain.
#[allow(clippy::too_many_arguments)]
fn write_pairs(
    buf1: &ChunkBuf,
    buf2: &ChunkBuf,
    pairs: Vec<ProcessedPair>,
    stats: &mut PairedStats,
    w1: &mut OutputWriter,
    w2: &mut OutputWriter,
    mut ws1: Option<&mut OutputWriter>,
    mut ws2: Option<&mut OutputWriter>,
    mut meta: Option<&mut OutputWriter>,
    phred_offset: u8,
    warn_r1_once: &AtomicBool,
    warn_r2_once: &AtomicBool,
    dedup: Option<&mut Dedup>,
) -> Result<()> {
    let mut dedup = dedup;
    for pp in pairs {
        // Pair-level decision: a duplicate pair drops both mates and yields
        // no singleton. Runs on this thread, so the sequence of decisions is
        // input order whatever the thread count.
        let s1 = &buf1.seq[pp.p1.input.seq.clone()];
        let s2 = &buf2.seq[pp.p2.input.seq.clone()];
        let dup = dedup
            .as_deref_mut()
            .is_some_and(|d| d.is_duplicate(s1, Some(s2)));
        let (res1, res2) = if dup {
            (duplicate_result(), duplicate_result())
        } else {
            (pp.p1.output.result.clone(), pp.p2.output.result.clone())
        };
        stats
            .r1
            .record(s1.len(), &res1, pp.p1.output.adapter_name.as_deref());
        stats
            .r2
            .record(s2.len(), &res2, pp.p2.output.adapter_name.as_deref());
        if let Some(mw) = meta.as_deref_mut() {
            metadata::write_record(
                &mut **mw,
                &buf1.id[pp.p1.input.id.clone()],
                Mate::R1,
                &res1,
                s1.len(),
            )
            .context("writing metadata record (R1)")?;
            metadata::write_record(
                &mut **mw,
                &buf2.id[pp.p2.input.id.clone()],
                Mate::R2,
                &res2,
                s2.len(),
            )
            .context("writing metadata record (R2)")?;
        }
        match outcome(&res1, &res2) {
            PairOutcome::BothKept => {
                pp.p1.emit_to(buf1, w1, phred_offset)?;
                pp.p2.emit_to(buf2, w2, phred_offset)?;
                stats.pairs_out += 1;
            }
            PairOutcome::BothDropped => {
                stats.pairs_dropped += 1;
            }
            PairOutcome::R1KeptR2Dropped => {
                if let Some(w) = ws1.as_deref_mut() {
                    pp.p1.emit_to(buf1, w, phred_offset)?;
                    stats.singletons_r1 += 1;
                } else {
                    stats.discarded_singletons_r1 += 1;
                    if !warn_r1_once.swap(true, Ordering::Relaxed) {
                        eprintln!(
                            "fastprep: warning: R1 survived while R2 was dropped, but \
                             --singletons1 was not provided; discarding. This warning \
                             will not repeat."
                        );
                    }
                }
            }
            PairOutcome::R1DroppedR2Kept => {
                if let Some(w) = ws2.as_deref_mut() {
                    pp.p2.emit_to(buf2, w, phred_offset)?;
                    stats.singletons_r2 += 1;
                } else {
                    stats.discarded_singletons_r2 += 1;
                    if !warn_r2_once.swap(true, Ordering::Relaxed) {
                        eprintln!(
                            "fastprep: warning: R2 survived while R1 was dropped, but \
                             --singletons2 was not provided; discarding. This warning \
                             will not repeat."
                        );
                    }
                }
            }
        }
    }
    Ok(())
}

fn outcome(r1: &TrimResult, r2: &TrimResult) -> PairOutcome {
    let keep1 = matches!(r1, TrimResult::Emit { .. });
    let keep2 = matches!(r2, TrimResult::Emit { .. });

    match (keep1, keep2) {
        (true, true) => PairOutcome::BothKept,
        (false, false) => PairOutcome::BothDropped,
        (true, false) => PairOutcome::R1KeptR2Dropped,
        (false, true) => PairOutcome::R1DroppedR2Kept,
    }
}

fn pair_chunks(
    buf1: &ChunkBuf,
    buf2: &ChunkBuf,
    c1: Vec<ScanInput>,
    c2: Vec<ScanInput>,
    params_r1: &Params,
    params_r2: &Params,
    parallel: bool,
) -> Vec<ProcessedPair> {
    debug_assert_eq!(c1.len(), c2.len());
    // Mates are trimmed independently. There is no cross-mate inference:
    // PE overlap detection was removed, so the adapter comes from the
    // configured set and nowhere else, and `trim` resolves it itself.
    let work = |(i1, i2): (ScanInput, ScanInput)| -> ProcessedPair {
        let s1 = &buf1.seq[i1.seq.clone()];
        let q1 = &buf1.qual[i1.qual.clone()];
        let s2 = &buf2.seq[i2.seq.clone()];
        let q2 = &buf2.qual[i2.qual.clone()];
        let out1 = trim(s1, q1, params_r1, AdapterHint::Detect);
        let out2 = trim(s2, q2, params_r2, AdapterHint::Detect);
        ProcessedPair {
            p1: Processed {
                input: i1,
                output: out1,
            },
            p2: Processed {
                input: i2,
                output: out2,
            },
        }
    };
    if parallel {
        c1.into_par_iter()
            .zip(c2.into_par_iter())
            .map(work)
            .collect()
    } else {
        c1.into_iter().zip(c2.into_iter()).map(work).collect()
    }
}

// ===========================================================================
// Integration tests.
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::MinLen;
    use crate::io::CompressionLevels;
    use flate2::read::GzDecoder;
    use std::fs::File;
    use std::io::{Read, Write};
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn tmp_path(name: &str) -> PathBuf {
        let id = TMP_COUNTER.fetch_add(1, AtomicOrdering::Relaxed);
        let mut p = std::env::temp_dir();
        p.push(format!("polyx-it-{}-{}-{}", std::process::id(), id, name));
        p
    }

    fn write_plain(path: &Path, body: &str) {
        let mut f = File::create(path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
    }

    fn read_plain(path: &Path) -> String {
        let mut s = String::new();
        File::open(path).unwrap().read_to_string(&mut s).unwrap();
        s
    }

    fn read_gz(path: &Path) -> String {
        let mut s = String::new();
        GzDecoder::new(File::open(path).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        s
    }

    fn default_rt() -> Runtime {
        Runtime {
            threads: 1,
            chunk: 8,
            compression: CompressionLevels::default(),
            phred_offset: 33,
            phred_auto: false,
            metadata: None,
            json: None,
            adapter_source: crate::cli::AdapterSource::None,
            adapter_error_rate: 0.1,
            min_adapter_overlap_se: 4,
            dedup: None,
            thread_pool: None,
        }
    }

    fn default_params() -> Params {
        Params::default()
    }

    // ------- Single-end plain --------------------------------------------

    #[test]
    fn single_end_plain_roundtrip() {
        // Two reads: one clean, one with a long polyA prefix that should
        // be dropped under defaults (70-A read, ΔL = 28, T_5 ≥ 70 > ΔL).
        let inp = tmp_path("in.fq");
        let out = tmp_path("out.fq");
        write_plain(
            &inp,
            "@clean\n\
             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
             +\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
             @doomed\n\
             AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\
             +\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        let stats = run_single(&inp, &out, &default_params(), &default_rt()).unwrap();
        assert_eq!(stats.reads_in, 2);
        assert_eq!(stats.reads_kept, 1);
        assert_eq!(stats.reads_dropped, 1);
        let body = read_plain(&out);
        assert!(body.contains("@clean"));
        assert!(!body.contains("@doomed"));
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }

    // ------- Single-end with gzip output ---------------------------------

    #[test]
    fn single_end_gz_output() {
        let inp = tmp_path("in.fq");
        let out = tmp_path("out.fq.gz");
        write_plain(
            &inp,
            "@r1\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        let stats = run_single(&inp, &out, &default_params(), &default_rt()).unwrap();
        assert_eq!(stats.reads_kept, 1);
        let body = read_gz(&out);
        assert!(body.contains("@r1"));
        assert!(body.contains("ACGTACGT"));
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }

    // ------- Single-end with zstd output ---------------------------------

    #[test]
    fn single_end_zst_output() {
        let inp = tmp_path("in.fq");
        let out = tmp_path("out.fq.zst");
        write_plain(
            &inp,
            "@r1\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        let stats = run_single(&inp, &out, &default_params(), &default_rt()).unwrap();
        assert_eq!(stats.reads_kept, 1);
        // Decompress and verify.
        let f = File::open(&out).unwrap();
        let mut dec = zstd::stream::read::Decoder::new(f).unwrap();
        let mut body = String::new();
        dec.read_to_string(&mut body).unwrap();
        assert!(body.contains("@r1"));
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }

    // ------- FASTA round-trip --------------------------------------------

    #[test]
    fn fasta_input_fasta_output() {
        // FASTA input, no quality. Q0 synthesis means all W = 0, so
        // matches contribute (1 − μ) = 0.9 each and mismatches contribute
        // exactly 0. The scan is purely match-count-driven, so mercy
        // never fires on a run of any length — the run must hit the
        // end of the scan window for the trim to register.
        //
        // To make this test deterministic: an all-A read of length 40
        // with ΔL = 16 triggers the budget short-circuit (drop); a
        // clean read is kept with no trim.
        let inp = tmp_path("in.fa");
        let out = tmp_path("out.fa");
        write_plain(
            &inp,
            ">clean\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
             >doomed\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        );
        let stats = run_single(&inp, &out, &default_params(), &default_rt()).unwrap();
        assert_eq!(stats.reads_in, 2);
        assert_eq!(stats.reads_kept, 1, "clean read should survive");
        assert_eq!(stats.reads_dropped, 1, "all-A read should drop by budget");
        let body = read_plain(&out);
        assert!(body.contains(">clean"));
        assert!(!body.contains(">doomed"));
        // Output is FASTA, not FASTQ — no `+` or `@` lines.
        assert!(!body.contains("\n+\n"));
        assert!(!body.contains("@"));
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }

    // ------- Paired-end with singleton routing ---------------------------

    #[test]
    fn paired_end_singleton_routing() {
        let in1 = tmp_path("r1.fq");
        let in2 = tmp_path("r2.fq");
        let out1 = tmp_path("out_r1.fq");
        let out2 = tmp_path("out_r2.fq");
        let sing1 = tmp_path("sing_r1.fq");
        let sing2 = tmp_path("sing_r2.fq");

        // Pair 1: both clean → both out.
        // Pair 2: R1 doomed (all-A), R2 clean → sing2 gets R2, sing1 empty.
        // Pair 3: R1 clean, R2 doomed (all-G) → sing1 gets R1, sing2 empty.
        write_plain(
            &in1,
            "@p1\n\
             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
             @p2\n\
             AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
             @p3\n\
             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        write_plain(
            &in2,
            "@p1\n\
             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
             @p2\n\
             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
             @p3\n\
             GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        );

        let p = default_params();
        let stats = run_paired(
            &in1,
            &in2,
            &out1,
            &out2,
            Some(&sing1),
            Some(&sing2),
            &p,
            &p,
            &default_rt(),
        )
        .unwrap();

        assert_eq!(stats.pairs_out, 1, "only p1 fully survives");
        assert_eq!(stats.singletons_r1, 1, "p3 R1 goes to sing1");
        assert_eq!(stats.singletons_r2, 1, "p2 R2 goes to sing2");
        assert_eq!(stats.discarded_singletons_r1, 0);
        assert_eq!(stats.discarded_singletons_r2, 0);

        let body_out1 = read_plain(&out1);
        let body_out2 = read_plain(&out2);
        let body_sing1 = read_plain(&sing1);
        let body_sing2 = read_plain(&sing2);

        assert!(
            body_out1.contains("@p1") && !body_out1.contains("@p2") && !body_out1.contains("@p3")
        );
        assert!(
            body_out2.contains("@p1") && !body_out2.contains("@p2") && !body_out2.contains("@p3")
        );
        assert!(body_sing1.contains("@p3"));
        assert!(body_sing2.contains("@p2"));

        for p in [&in1, &in2, &out1, &out2, &sing1, &sing2] {
            std::fs::remove_file(p).ok();
        }
    }

    // ------- Paired-end singleton discard (no --singletons*) -------------

    #[test]
    fn paired_end_discards_singleton_when_path_omitted() {
        let in1 = tmp_path("r1.fq");
        let in2 = tmp_path("r2.fq");
        let out1 = tmp_path("out_r1.fq");
        let out2 = tmp_path("out_r2.fq");

        // One pair: R1 doomed, R2 clean → R2 should be discarded.
        write_plain(
            &in1,
            "@p1\n\
             AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        write_plain(
            &in2,
            "@p1\n\
             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        let p = default_params();
        let stats =
            run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &default_rt()).unwrap();
        assert_eq!(stats.pairs_out, 0);
        assert_eq!(stats.singletons_r2, 0);
        assert_eq!(stats.discarded_singletons_r2, 1);
        assert_eq!(stats.discarded_singletons_r1, 0);

        for p in [&in1, &in2, &out1, &out2] {
            std::fs::remove_file(p).ok();
        }
    }

    // ------- Parallel execution preserves order --------------------------

    #[test]
    fn parallel_output_matches_serial() {
        let inp = tmp_path("big.fq");
        let out_serial = tmp_path("out_serial.fq");
        let out_parallel = tmp_path("out_parallel.fq");

        let mut body = String::new();
        for i in 0..200 {
            body.push_str(&format!(
                "@r{}\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
                i
            ));
        }
        write_plain(&inp, &body);

        let p = default_params();
        let rt_s = Runtime {
            threads: 1,
            chunk: 16,
            compression: CompressionLevels::default(),
            phred_offset: 33,
            phred_auto: false,
            metadata: None,
            json: None,
            adapter_source: crate::cli::AdapterSource::None,
            adapter_error_rate: 0.1,
            min_adapter_overlap_se: 4,
            dedup: None,
            thread_pool: None,
        };
        // Use a scoped pool so this test doesn't race with or pollute the
        // global rayon pool that other tests may rely on.
        let pool = std::sync::Arc::new(
            rayon::ThreadPoolBuilder::new()
                .num_threads(4)
                .build()
                .unwrap(),
        );
        let rt_p = Runtime {
            threads: 4,
            chunk: 16,
            compression: CompressionLevels::default(),
            phred_offset: 33,
            phred_auto: false,
            metadata: None,
            json: None,
            adapter_source: crate::cli::AdapterSource::None,
            adapter_error_rate: 0.1,
            min_adapter_overlap_se: 4,
            dedup: None,
            thread_pool: Some(pool),
        };
        run_single(&inp, &out_serial, &p, &rt_s).unwrap();
        run_single(&inp, &out_parallel, &p, &rt_p).unwrap();

        let a = read_plain(&out_serial);
        let b = read_plain(&out_parallel);
        assert_eq!(a, b, "parallel output must match serial");

        for p in [&inp, &out_serial, &out_parallel] {
            std::fs::remove_file(p).ok();
        }
    }

    // ------- Fractional min_len ------------------------------------------

    #[test]
    fn fractional_min_len_respected() {
        let inp = tmp_path("in.fq");
        let out = tmp_path("out.fq");
        write_plain(
            &inp,
            "@clean\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n\
             @trimmable\nAAAAAAAAAAACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n",
        );
        let mut p = default_params();
        p.min_len = MinLen::Frac(1.0);
        let stats = run_single(&inp, &out, &p, &default_rt()).unwrap();
        // Clean read survives (T = 0); polyA read drops (any trim > 0).
        assert_eq!(stats.reads_kept, 1);
        assert_eq!(stats.reads_dropped, 1);
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }
}
