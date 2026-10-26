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

use crate::algorithm::{trim, DropReason, Params, TrimOutput, TrimResult};
use crate::cli::{PeOverlapMode, Runtime};
use crate::dedup::{AnchorSpec, Fingerprint, Seen, Verifier};
use crate::io::{
    open_input, open_output_as, Out, OutputWriter, RecordFormat,
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

// Counts freshly allocated chunk arenas. Test builds only: recycling is
// invisible in the output, so only a counter can show it happens. Thread-local
// because the pipeline allocates on its own (ordered) thread, and a global
// would be polluted by the other tests running in parallel with this one.
#[cfg(test)]
thread_local! {
    static ARENAS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Take a recycled arena, or allocate one.
///
/// Spelled as a `match` rather than `unwrap_or_default()` on purpose: `Default`
/// would bypass the test-only `ARENAS` counter that the recycling test asserts
/// on.
#[inline]
fn pop_or_new(free: &mut Vec<ChunkBuf>) -> ChunkBuf {
    match free.pop() {
        Some(buf) => buf,
        None => ChunkBuf::new(),
    }
}

impl ChunkBuf {
    /// A fresh arena. Only the free-list miss path should reach this: a chunk
    /// at 1024 reads and 151 bp is ~464 KB per mate, so allocating one per
    /// chunk (and freeing the previous one only after reserving the next) was
    /// ~4.5 GB of large-object churn in a 5 M-pair run.
    fn new() -> Self {
        #[cfg(test)]
        ARENAS.with(|c| c.set(c.get() + 1));
        Self::default()
    }

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
                // Quality stays in its INPUT encoding: it is copied through,
                // used through a LUT shifted by `phred_offset`, and written
                // back verbatim. The one transformation kept is the clamp —
                // a byte below the offset used to decode to 0 and re-encode to
                // `offset`, so it is stored as `offset` here and the output
                // bytes are identical to the decode/encode path's.
                //
                // `max` rather than a branch: the whole slice is a memcpy plus
                // one SIMD max, where decoding was a subtract per byte and
                // egress an add per byte (two full passes over 151 bytes).
                let q_start = buf.qual.len();
                if phred_offset == 0 {
                    buf.qual.extend_from_slice(q);
                } else {
                    buf.qual.extend(q.iter().map(|&b| b.max(phred_offset)));
                }
                (q_start..buf.qual.len(), RecordFormat::Fastq)
            }
            None => {
                // FASTA carries no qualities. Zeros (not ASCII) because
                // `gate_lut[0]` is Q0 either way, and the mean-q filter is
                // skipped for FASTA records (`fasta_params`).
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
    /// Consensus-corrected `(seq, qual)` in FULL-read coordinates, present
    /// only when `--pe-correct` rewrote at least one retained base. `emit_to`
    /// slices the same trim window out of these as out of the buffer, so trim
    /// coordinates never move and the no-edit path is byte-identical.
    corrected: Option<(Vec<u8>, Vec<u8>)>,
}

impl Processed {
    fn emit_to(&self, buf: &ChunkBuf, w: &mut Out) -> Result<()> {
        match &self.output.result {
            TrimResult::Drop { .. } => Ok(()),
            TrimResult::Emit { .. } => {
                let t_5 = self.output.result.t5();
                let t_3 = self.output.result.t3();
                // With `--pe-correct` edits the full-read corrected copy stands
                // in for the buffer slices; the trim window is taken the same
                // way either way, so flag-off output is byte-identical.
                let seq_all: &[u8] = match &self.corrected {
                    Some((c_seq, _)) => c_seq.as_slice(),
                    None => &buf.seq[self.input.seq.clone()],
                };
                debug_assert!(
                    t_5 + t_3 <= seq_all.len(),
                    "trim exceeds read length: t_5={} t_3={} L={}",
                    t_5,
                    t_3,
                    seq_all.len()
                );
                let seq = &seq_all[t_5..seq_all.len() - t_3];
                let id = &buf.id[self.input.id.clone()];
                match self.input.format {
                    RecordFormat::Fastq => {
                        let qual_all: &[u8] = match &self.corrected {
                            Some((_, c_qual)) => c_qual.as_slice(),
                            None => &buf.qual[self.input.qual.clone()],
                        };
                        let raw_qual = &qual_all[t_5..qual_all.len() - t_3];
                        w.fastq(id, seq, raw_qual)
                    }
                    RecordFormat::Fasta => w.fasta(id, seq),
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

/// `params` with the mean-quality filter removed when the record is FASTA.
///
/// FASTA carries no qualities, and ingestion fills them with 0 (`ScanInput::
/// push_needletail`), so a mean-quality filter computed over a FASTA record
/// always fails: with the shipped `--mean-q 25`, *every* record of a FASTA
/// input was dropped. The filter is undefined without quality evidence, so it
/// is skipped for FASTA records rather than fed an invented quality; the polyX
/// gates are unaffected (matches are quality-independent by construction, and
/// a mismatch on a Q0 base was already free).
fn fasta_params(params: &Params) -> Option<Params> {
    params.mean_q.map(|_| Params {
        mean_q: None,
        ..params.clone()
    })
}

/// The configured adapter names in set order, for the report.
///
/// Called once per run. Per-read work carries an adapter *index*
/// (`TrimOutput::adapter`), so these strings are materialised only when the
/// report is printed — `Stats` used to allocate one per adapter-trimmed read,
/// twice.
fn adapter_names(params: &Params) -> Vec<String> {
    params
        .adapter
        .as_ref()
        .map(|set| set.adapters.iter().map(|a| a.name.clone()).collect())
        .unwrap_or_default()
}

fn params_for<'a>(params: &'a Params, fasta: &'a Option<Params>, fmt: RecordFormat) -> &'a Params {
    match (fmt, fasta) {
        (RecordFormat::Fasta, Some(p)) => p,
        _ => params,
    }
}

fn process_chunk(
    buf: &ChunkBuf,
    chunk: Vec<ScanInput>,
    params: &Params,
    parallel: bool,
) -> Vec<Processed> {
    let fasta = fasta_params(params);
    let work = |input: ScanInput| -> Processed {
        let seq = &buf.seq[input.seq.clone()];
        let qual = &buf.qual[input.qual.clone()];
        let output = trim(seq, qual, params_for(params, &fasta, input.format));
        Processed {
            input,
            output,
            corrected: None,
        }
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

/// Derive a `.part` write-side path for `p`. On success the caller renames
/// this path to `p`; on error the guard removes it.
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

/// What the quality-byte census concluded.
///
/// The census cannot *classify* a modern file: Phred+33 qualities at Q≥40 use
/// every byte from `'@'` up, exactly like Phred+64 at Q≥0. It can only prove
/// +33 (one byte below `'@'`) or abstain. The count is kept because D-6 says
/// every assumption names its evidence, and the old code parsed 1000 records
/// and then threw the counters away.
#[derive(Debug, Clone)]
pub struct PhredCensus {
    pub offset: u8,
    pub records: usize,
    pub bases: usize,
    /// Lowest quality byte seen; `None` when the input carries no qualities.
    pub min_q: Option<u8>,
    /// Bytes above `'J'` (74) — reported as evidence, never as proof.
    pub high: usize,
    /// A byte below `'@'` (64) cannot occur under Phred+64, so +33 is proven.
    pub proven: bool,
}

impl PhredCensus {
    /// The D-6 assumption record for the offset this census settled on.
    pub fn assumption(&self) -> crate::stats::Assumption {
        let basis;
        let consequence;
        if let Some(min_q) = self.min_q {
            if self.proven {
                basis = format!(
                    "census of {} records ({} quality bytes): lowest byte {} < 64, \
                     which cannot occur under Phred+64, so +33 is proven",
                    self.records, self.bases, min_q
                );
                consequence = "none; the byte range proves the offset".to_string();
            } else {
                basis = format!(
                    "census of {} records ({} quality bytes): lowest byte {} >= 64 and {} \
                     bytes above 74; a high-byte census cannot distinguish +33 from +64, \
                     so +33 is assumed",
                    self.records, self.bases, min_q, self.high
                );
                consequence = "if the input is legacy Phred+64 (Illumina 1.3/1.5), every \
                               quality is read 31 points high and --mean-q can drop the \
                               whole file; pass --phred-offset 64"
                    .to_string();
            }
        } else {
            basis = format!(
                "no quality bytes in the first {} records (FASTA or empty input); \
                 +33 is the FASTQ default",
                self.records
            );
            consequence = "none; no qualities are decoded".to_string();
        }
        crate::stats::Assumption {
            what: "phred_offset",
            value: serde_json::json!(self.offset),
            basis,
            consequence,
        }
    }
}

/// Peek at the first 1000 records of a FASTQ file and infer what the Phred
/// offset can be proven to be.
///
/// Heuristic: any byte `< '@'` (64) cannot appear in Phred+64 (which
/// starts at `'@'`) → Phred+33 is proven. If all observed bytes are ≥ 64,
/// nothing is proven; the offset stays at the modern default of 33 and the
/// assumption record says so. FASTA files (no quality) default to 33.
pub fn detect_phred_offset(path: &Path) -> Result<PhredCensus> {
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("auto-detecting Phred offset from {}", path.display()))?;
    let mut min_q = u8::MAX;
    let mut n_bases = 0usize;
    let mut n_high = 0usize;
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
                        if b > 74 {
                            n_high += 1;
                        }
                        n_bases += 1;
                    }
                }
                n += 1;
            }
        }
    }
    if n_bases == 0 || min_q == u8::MAX {
        return Ok(PhredCensus {
            offset: 33,
            records: n,
            bases: 0,
            min_q: None,
            high: 0,
            proven: false,
        }); // FASTA or empty file: modern default
    }
    // A byte below 64 cannot occur under Phred+64, so its presence PROVES
    // Phred+33. Its absence proves nothing: a Phred+33 file whose lowest
    // quality is Q31 also sits entirely at or above 64, and modern chemistry
    // produces Q>41 routinely (Q43 = 'L' = 76), so a high-byte census is not
    // evidence for +64 either. The old code took the absence as evidence for
    // 64; a later revision took a preponderance of bytes above 74 as evidence,
    // and that still reads an all-Q43 Phred+33 library as +64, subtracts 31
    // from every quality and lets --mean-q delete the entire file (reproduced
    // 2026-10-20: 5000/5000 reads dropped, exit 0).
    //
    // So the auto path never asserts +64. It reports +33 — the modern standard
    // — and records the census as the basis of that assumption, naming the
    // flag that settles it. This is exactly fastcheck::accum::offset_verdict's
    // rule: a wrong guess here is silent, total and unrecoverable, so
    // abstaining plus a labelled assumption is the only honest default. Legacy
    // Phred+64 data (Illumina 1.3/1.5) is still supported, by asking for it.
    Ok(PhredCensus {
        offset: 33,
        records: n,
        bases: n_bases,
        min_q: Some(min_q),
        high: n_high,
        proven: min_q < 64,
    })
}

/// The assumption record for an offset the caller stated explicitly.
fn explicit_phred_assumption(offset: u8) -> crate::stats::Assumption {
    crate::stats::Assumption {
        what: "phred_offset",
        value: serde_json::json!(offset),
        basis: format!("explicit --phred-offset {offset}"),
        consequence: "none; the caller stated the offset".to_string(),
    }
}

// ---------------------------------------------------------------------------
// Single-end.
// ---------------------------------------------------------------------------

/// The message a total-drop refusal carries.
///
/// Shared by the pre-commit refusal here and the post-run guard in `main`, so
/// the two paths cannot drift apart.
pub fn total_drop_message(reads_in: u64, mode: &str, outputs: &str) -> String {
    format!(
        "fastprep: dropped all {reads_in} {mode} reads — nothing was written to {outputs}. \
         Refusing to report success on a total drop. Check the per-cause drop \
         table above; the usual causes are --min-len above the trimmed insert \
         length, --mean-q above the library's mean quality, and a Phred offset \
         read the wrong way (pass --phred-offset 64 for legacy data). Pass \
         --allow-total-drop to accept an empty output."
    )
}

/// Refuse a record format the destination path contradicts.
///
/// fastprep preserves the input's record format (FASTA in, FASTA out) and
/// infers only the *compression* from the destination extension. That is a
/// footgun when the name says otherwise: a FASTA input written to `out.fq`
/// produced a FASTA file with a FASTQ name, and the next tool in the chain
/// refused it with a quality-length error that pointed at the wrong step.
/// CLI-SPEC §4 states named paths use the extension, so honour it here.
fn check_output_format(out: &Path, fmt: RecordFormat) -> Result<()> {
    let ext = out
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let says_fastq = matches!(ext.as_str(), "fq" | "fastq");
    let says_fasta = matches!(ext.as_str(), "fa" | "fasta" | "fna" | "faa");
    match fmt {
        RecordFormat::Fasta if says_fastq => anyhow::bail!(
            "fastprep: {} is a FASTA input but the output {} says FASTQ. \
             Qualities cannot be invented, so name the output .fa/.fasta \
             (or drop the extension) and the run will keep the record format.",
            "the input",
            out.display()
        ),
        RecordFormat::Fastq if says_fasta => anyhow::bail!(
            "fastprep: the input is FASTQ but the output {} says FASTA. \
             Name the output .fq/.fastq (or drop the extension) to keep the \
             record format.",
            out.display()
        ),
        _ => Ok(()),
    }
}

/// Refuse a total drop **before** the atomic rename.
///
/// The refusal used to run in `main` after `run_single`/`run_paired` had
/// already committed, so a run that dropped every read still replaced the
/// destination with a 0-byte file — and then claimed "nothing was written to
/// the output file". With `--force` that destroyed a good previous output.
fn refuse_before_commit(
    reads_in: u64,
    written: u64,
    allow: bool,
    mode: &str,
    outputs: &str,
) -> Result<()> {
    if allow || reads_in == 0 || written > 0 {
        return Ok(());
    }
    anyhow::bail!("{}", total_drop_message(reads_in, mode, outputs))
}

pub fn run_single(in1: &Path, out1: &Path, params: &Params, rt: &Runtime) -> Result<Stats> {
    // Phred offset: a census of the input, or the configured value. Either way
    // the choice is recorded as an assumption: a +64 file read as +33 is a
    // silent 31-point error in every quality and D-6 forbids guessing quietly.
    let (phred_offset, assumption) = if rt.phred_auto {
        let census = detect_phred_offset(in1)?;
        (census.offset, census.assumption())
    } else {
        (rt.phred_offset, explicit_phred_assumption(rt.phred_offset))
    };

    let mut reader = parse_fastx_reader(open_input(in1, rt.compression.gz_threads)?)
        .with_context(|| format!("opening input {}", in1.display()))?;

    // Open all outputs to `.part` paths; rename atomically on success.
    // Under --dry-run no output is opened at all: every counter is still
    // computed, nothing is written (CLI-SPEC §1).
    let tmp_out = part_path(out1);
    let guard_out = TmpFile::new(tmp_out.clone());
    let mut writer = if rt.dry_run {
        None
    } else {
        Some(Out::new(open_output_as(&tmp_out, out1, rt.compression)?))
    };

    let tmp_meta_path = rt.metadata.as_deref().map(part_path);
    let guard_meta = tmp_meta_path.as_ref().map(|p| TmpFile::new(p.clone()));
    let mut meta = match (rt.metadata.as_deref(), tmp_meta_path.as_deref()) {
        (Some(fp), Some(tp)) if !rt.dry_run => {
            let mut w = open_output_as(tp, fp, METADATA_COMPRESSION)?;
            metadata::write_file_header(&mut *w).context("writing metadata header")?;
            Some(w)
        }
        _ => None,
    };

    let mut stats = Stats {
        phred_offset,
        assumptions: vec![assumption],
        ..Default::default()
    };
    // Names are resolved once for the run; per-read work carries an index.
    stats.adapter.set_names(adapter_names(params));
    let parallel = rt.threads > 1;

    // Quality is never decoded: the offset is folded into the gate LUT and the
    // mean-q sum, so the params carry it. Set before the first chunk.
    let mut params = params.clone();
    params.set_phred_offset(phred_offset);
    let params = params;

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
        // Recycled arenas: the buffer a chunk was computed in comes back on
        // `res_rx` and is reused for the next chunk instead of being freed
        // after the next one was allocated.
        let mut free: Vec<ChunkBuf> = Vec::new();
        let mut buf = ChunkBuf::new();
        buf.reserve(rt.chunk, 151);
        let mut chunk = Vec::with_capacity(rt.chunk);
        fill_chunk(&mut *reader, &mut buf, &mut chunk, rt.chunk, phred_offset)?;
        if let Some(first) = chunk.first() {
            check_output_format(out1, first.format)?;
        }

        while !chunk.is_empty() {
            raw_tx.send((buf, chunk)).ok();

            let mut next_buf = pop_or_new(&mut free);
            next_buf.reserve(rt.chunk, 151);
            let mut next = Vec::with_capacity(rt.chunk);
            fill_chunk(
                &mut *reader,
                &mut next_buf,
                &mut next,
                rt.chunk,
                phred_offset,
            )?;

            let (mut res_buf, results) = res_rx.recv().context("compute thread died")?;
            for p in &results {
                // Deduplication runs on the ordered thread, so the decision
                // sequence is the input order regardless of thread count.
                let seq = &res_buf.seq[p.input.seq.clone()];
                let dup = dedup.as_mut().is_some_and(|d| d.is_duplicate(seq, None));
                // Borrowed, never cloned: `TrimResult::Emit` carries a heap
                // `Vec` of interior runs and this is the per-record path on the
                // ordered thread. `duplicate_result()` allocates nothing.
                let dup_result;
                let result: &TrimResult = if dup {
                    dup_result = duplicate_result();
                    &dup_result
                } else {
                    &p.output.result
                };
                stats.record(seq.len(), result, p.output.adapter);
                if !dup {
                    stats.record_margin(p.output.margin_5, p.output.margin_3);
                    if let Some(w) = writer.as_mut() {
                        p.emit_to(&res_buf, w)?;
                    }
                }
                if let Some(mw) = meta.as_mut() {
                    metadata::write_record(
                        &mut **mw,
                        &res_buf.id[p.input.id.clone()],
                        Mate::R1,
                        result,
                        seq.len(),
                    )
                    .context("writing metadata record")?;
                }
            }
            chunks_written += 1;
            if !rt.quiet && stderr_is_tty && chunks_written.is_multiple_of(100) {
                eprint!("\r  {} reads processed", stats.reads_in);
                progress_shown = true;
            }
            // The arena just written from goes back on the free list, cleared,
            // so the next chunk reuses its capacity instead of re-mapping.
            res_buf.clear();
            free.push(res_buf);
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
                p.output.adapter,
            );
            stats.record_margin(p.output.margin_5, p.output.margin_3);
            if io_result.is_ok() {
                if let Some(w) = writer.as_mut() {
                    p.emit_to(&res_buf, w)
                        .context("writing output after the compute thread finished")?;
                }
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
    if let Some(mut w) = writer.take() {
        // A flush error means the last staged block never reached the file.
        // This used to be `.ok()`: the truncated file was then renamed into
        // place and the run exited 0 (reproduced with `ulimit -f 1`: 1024 of
        // 6440 bytes kept, empty stderr, exit 0).
        w.flush().context("flushing output")?;
    }
    if let Some(mut mw) = meta.take() {
        // D-6: a metadata file that cannot be flushed is truncated, not
        // complete. This used to be `.ok()`, which reported success over a
        // partial file; the staging buffer makes that loss strictly larger.
        mw.flush().context("flushing metadata")?;
    }

    io_result?;

    // Atomic rename: committed output replaces the destination path.
    // A dry run never got here with a file to commit.
    if !rt.dry_run {
        // Refuse a total drop here, while the destination is still untouched:
        // the temp file is removed by its guard on the error path.
        refuse_before_commit(
            stats.reads_in,
            stats.reads_kept,
            rt.allow_total_drop,
            "single-end",
            "the output file",
        )?;
        std::fs::rename(&tmp_out, out1)
            .with_context(|| format!("committing {} → {}", tmp_out.display(), out1.display()))?;
        guard_out.commit();
        if let (Some(tp), Some(fp), Some(g)) = (tmp_meta_path, rt.metadata.as_deref(), guard_meta) {
            std::fs::rename(&tp, fp)
                .with_context(|| format!("committing {} → {}", tp.display(), fp.display()))?;
            g.commit();
        }
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
    /// Retained bases rewritten by `--pe-correct` consensus correction, both
    /// mates summed. 0 whenever the flag is off (the default).
    pub bases_corrected: u64,
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
    /// Retained bases rewritten by `--pe-correct` on this pair (both mates).
    /// Carried from the work closure so `write_pairs` sums it without
    /// re-deriving anything from the emitted buffers.
    corrected_bases: usize,
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
    // Phred offset: a census of R1 (R2 must match), or the configured value.
    let (phred_offset, phred_assumption) = if rt.phred_auto {
        let census = detect_phred_offset(in1)?;
        (census.offset, census.assumption())
    } else {
        (rt.phred_offset, explicit_phred_assumption(rt.phred_offset))
    };

    let mut r1 = parse_fastx_reader(open_input(in1, rt.compression.gz_threads)?)
        .with_context(|| format!("opening input {}", in1.display()))?;
    let mut r2 = parse_fastx_reader(open_input(in2, rt.compression.gz_threads)?)
        .with_context(|| format!("opening input {}", in2.display()))?;

    // Open all outputs to `.part` paths; rename atomically on success.
    // Under --dry-run nothing is opened: every counter is still computed.
    let tmp_w1 = part_path(out1);
    let tmp_w2 = part_path(out2);
    let guard_w1 = TmpFile::new(tmp_w1.clone());
    let guard_w2 = TmpFile::new(tmp_w2.clone());
    let mut w1 = if rt.dry_run {
        None
    } else {
        Some(Out::new(open_output_as(&tmp_w1, out1, rt.compression)?))
    };
    let mut w2 = if rt.dry_run {
        None
    } else {
        Some(Out::new(open_output_as(&tmp_w2, out2, rt.compression)?))
    };

    let tmp_ws1 = singletons1.map(part_path);
    let tmp_ws2 = singletons2.map(part_path);
    let guard_ws1 = tmp_ws1.as_ref().map(|p| TmpFile::new(p.clone()));
    let guard_ws2 = tmp_ws2.as_ref().map(|p| TmpFile::new(p.clone()));
    let mut ws1 = if rt.dry_run {
        None
    } else {
        tmp_ws1
            .as_deref()
            .zip(singletons1)
            .map(|(tp, fp)| open_output_as(tp, fp, rt.compression).map(Out::new))
            .transpose()?
    };
    let mut ws2 = if rt.dry_run {
        None
    } else {
        tmp_ws2
            .as_deref()
            .zip(singletons2)
            .map(|(tp, fp)| open_output_as(tp, fp, rt.compression).map(Out::new))
            .transpose()?
    };

    let tmp_meta_path = rt.metadata.as_deref().map(part_path);
    let guard_meta = tmp_meta_path.as_ref().map(|p| TmpFile::new(p.clone()));
    let mut meta = match (rt.metadata.as_deref(), tmp_meta_path.as_deref()) {
        (Some(fp), Some(tp)) if !rt.dry_run => {
            let mut w = open_output_as(tp, fp, METADATA_COMPRESSION)?;
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
        bases_corrected: 0,
    };

    let parallel = rt.threads > 1;
    // The offset census runs once, from R1, and applies to both mates.
    stats.r1.phred_offset = phred_offset;
    stats.r2.phred_offset = phred_offset;
    stats.r1.assumptions.push(phred_assumption.clone());
    stats.r2.assumptions.push(phred_assumption);
    // Names are resolved once for the run; per-read work carries an index.
    let names = adapter_names(params_r1);
    stats.r1.adapter.set_names(names.clone());
    stats.r2.adapter.set_names(names);
    // Quality is never decoded: the offset is folded into the gate LUT and the
    // mean-q sum, so both mates' params carry it before the first chunk.
    let mut params_r1 = params_r1.clone();
    params_r1.set_phred_offset(phred_offset);
    let mut params_r2 = params_r2.clone();
    params_r2.set_phred_offset(phred_offset);

    // PE overlap: resolved once per chunk, per pair, inside the work closure.
    // The mode resolved to `Never` (the default with auto-detection or
    // --no-adapter) ⇒ the detector is never called and `pe_adapter` stays
    // `None`, which is exactly the flag-off code path (byte-identical output).
    let pe_anchor = crate::pe_overlap::ANCHOR;

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
    let pe_overlap = rt.pe_overlap;
    let adapter_evidence = rt.adapter_evidence;
    let pe_ambiguous_smallest = rt.pe_ambiguous_smallest;
    // (gate, min_ovl) for `--pe-correct`; `None` = flag off = untouched path.
    let pe_correct = rt.pe_correct.map(|g| (g, rt.pe_correct_min_overlap));
    let compute = std::thread::spawn(move || {
        while let Ok(((b1, c1), (b2, c2))) = raw_rx.recv() {
            let results = match thread_pool.as_ref() {
                Some(pool) => pool.install(|| {
                    pair_chunks(
                        &b1,
                        &b2,
                        c1,
                        c2,
                        &params_r1_arc,
                        &params_r2_arc,
                        parallel,
                        pe_overlap,
                        pe_anchor,
                        adapter_evidence,
                        pe_ambiguous_smallest,
                        pe_correct,
                    )
                }),
                None => pair_chunks(
                    &b1,
                    &b2,
                    c1,
                    c2,
                    &params_r1_arc,
                    &params_r2_arc,
                    parallel,
                    pe_overlap,
                    pe_anchor,
                    adapter_evidence,
                    pe_ambiguous_smallest,
                    pe_correct,
                ),
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
        // One free list serves both mates: an arena is an arena.
        let mut free: Vec<ChunkBuf> = Vec::new();
        let mut buf1 = ChunkBuf::new();
        buf1.reserve(rt.chunk, 151);
        let mut buf2 = ChunkBuf::new();
        buf2.reserve(rt.chunk, 151);
        let mut chunk1 = Vec::with_capacity(rt.chunk);
        let mut chunk2 = Vec::with_capacity(rt.chunk);
        fill_chunk(&mut *r1, &mut buf1, &mut chunk1, rt.chunk, phred_offset)?;
        fill_chunk(&mut *r2, &mut buf2, &mut chunk2, rt.chunk, phred_offset)?;
        if let Some(first) = chunk1.first() {
            check_output_format(out1, first.format)?;
        }
        if let Some(first) = chunk2.first() {
            check_output_format(out2, first.format)?;
        }
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

            let mut next_buf1 = pop_or_new(&mut free);
            next_buf1.reserve(rt.chunk, 151);
            let mut next_buf2 = pop_or_new(&mut free);
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

            let ((mut res_b1, mut res_b2), pairs) = res_rx.recv().context("compute thread died")?;
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
                &warn_r1_once,
                &warn_r2_once,
                dedup.as_mut(),
            )?;

            chunks_written += 1;
            if !rt.quiet && stderr_is_tty && chunks_written.is_multiple_of(100) {
                eprint!(
                    "\r  {} reads processed",
                    stats.r1.reads_in + stats.r2.reads_in
                );
                progress_shown = true;
            }
            res_b1.clear();
            free.push(res_b1);
            res_b2.clear();
            free.push(res_b2);
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

    if let Some(mut w) = w1.take() {
        w.flush().context("flushing R1 output")?;
    }
    if let Some(mut w) = w2.take() {
        w.flush().context("flushing R2 output")?;
    }
    if let Some(mut w) = ws1.take() {
        w.flush().context("flushing R1 singletons")?;
    }
    if let Some(mut w) = ws2.take() {
        w.flush().context("flushing R2 singletons")?;
    }
    if let Some(mut mw) = meta.take() {
        // D-6: a metadata file that cannot be flushed is truncated, not
        // complete. This used to be `.ok()`, which reported success over a
        // partial file; the staging buffer makes that loss strictly larger.
        mw.flush().context("flushing metadata")?;
    }

    io_result?;

    // Atomic renames. A dry run wrote nothing, so it commits nothing.
    if !rt.dry_run {
        // Refuse a total drop before any rename: pairs plus salvaged
        // singletons are what the output actually holds.
        refuse_before_commit(
            stats.r1.reads_in,
            stats.pairs_out + stats.singletons_r1 + stats.singletons_r2,
            rt.allow_total_drop,
            "paired-end",
            "the R1/R2 outputs or either singleton file",
        )?;
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
    w1: &mut Option<Out>,
    w2: &mut Option<Out>,
    mut ws1: Option<&mut Out>,
    mut ws2: Option<&mut Out>,
    mut meta: Option<&mut OutputWriter>,
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
        // Borrowed, never cloned: `TrimResult::Emit` carries a heap `Vec` of
        // interior runs, and this runs per record on the ordered thread. The
        // duplicate case needs a zeroed result, which allocates nothing
        // (`runs: vec![]`).
        let dup1;
        let dup2;
        let (res1, res2): (&TrimResult, &TrimResult) = if dup {
            dup1 = duplicate_result();
            dup2 = duplicate_result();
            (&dup1, &dup2)
        } else {
            (&pp.p1.output.result, &pp.p2.output.result)
        };
        stats
            .r1
            .record(s1.len(), res1, pp.p1.output.adapter);
        stats
            .r2
            .record(s2.len(), res2, pp.p2.output.adapter);
        if !dup {
            stats
                .r1
                .record_margin(pp.p1.output.margin_5, pp.p1.output.margin_3);
            stats
                .r2
                .record_margin(pp.p2.output.margin_5, pp.p2.output.margin_3);
            // Only emitted bases count: a duplicate pair's corrected copies
            // never reach the output.
            stats.bases_corrected += pp.corrected_bases as u64;
        }
        if let Some(mw) = meta.as_deref_mut() {
            metadata::write_record(
                &mut **mw,
                &buf1.id[pp.p1.input.id.clone()],
                Mate::R1,
                res1,
                s1.len(),
            )
            .context("writing metadata record (R1)")?;
            metadata::write_record(
                &mut **mw,
                &buf2.id[pp.p2.input.id.clone()],
                Mate::R2,
                res2,
                s2.len(),
            )
            .context("writing metadata record (R2)")?;
        }
        match outcome(res1, res2) {
            PairOutcome::BothKept => {
                if let (Some(w1), Some(w2)) = (w1.as_mut(), w2.as_mut()) {
                    pp.p1.emit_to(buf1, w1)?;
                    pp.p2.emit_to(buf2, w2)?;
                }
                stats.pairs_out += 1;
            }
            PairOutcome::BothDropped => {
                stats.pairs_dropped += 1;
            }
            PairOutcome::R1KeptR2Dropped => {
                if let Some(w) = ws1.as_deref_mut() {
                    pp.p1.emit_to(buf1, w)?;
                    stats.singletons_r1 += 1;
                } else {
                    stats.discarded_singletons_r1 += 1;
                    // The mate was kept by the trimmer and thrown away here;
                    // `reads_out` must count records, not trim decisions.
                    stats
                        .r1
                        .unkeep(s1.len().saturating_sub(res1.t5() + res1.t3()));
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
                    pp.p2.emit_to(buf2, w)?;
                    stats.singletons_r2 += 1;
                } else {
                    stats.discarded_singletons_r2 += 1;
                    stats
                        .r2
                        .unkeep(s2.len().saturating_sub(res2.t5() + res2.t3()));
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

#[allow(clippy::too_many_arguments)]
fn pair_chunks(
    buf1: &ChunkBuf,
    buf2: &ChunkBuf,
    c1: Vec<ScanInput>,
    c2: Vec<ScanInput>,
    params_r1: &Params,
    params_r2: &Params,
    parallel: bool,
    pe_overlap: PeOverlapMode,
    pe_anchor: usize,
    adapter_evidence: crate::adapter::AdapterEvidence,
    pe_ambiguous_smallest: bool,
    pe_correct: Option<(f64, usize)>,
) -> Vec<ProcessedPair> {
    debug_assert_eq!(c1.len(), c2.len());
    // Mates are trimmed independently, except for the one cross-mate inference
    // `--pe-overlap` opts into: when the pair's overlap pins the adapter
    // boundary AND the readthrough past it is real adapter sequence, that
    // boundary is handed to each mate's trim as an override. The adapter set is
    // passed to the detector precisely so a boundary without adapter behind it
    // is not applied -- an overlap proves where the fragment ends, not that
    // anything was ligated past it. With the mode off the detector is never
    // called, so `pe_adapter` stays `None` on every Params and the trim path
    // below is exactly the pre-flag one (byte-identical output).
    let base1 = Arc::new(params_r1.clone());
    let base2 = Arc::new(params_r2.clone());
    let work = |(i1, i2): (ScanInput, ScanInput)| -> ProcessedPair {
        let s1 = &buf1.seq[i1.seq.clone()];
        let q1 = &buf1.qual[i1.qual.clone()];
        let s2 = &buf2.seq[i2.seq.clone()];
        let q2 = &buf2.qual[i2.qual.clone()];

        // The detection result is per-mate: the boundary falls at position I
        // on each, and each mate's adapter is read 3′, so mate i gets
        // `(0, l_i - I)` — the same number only when the mates are the same
        // length. `Arc::make_mut` clones only when the mode actually fired.
        let (p1, p2) = if pe_overlap.runs() {
            match crate::pe_overlap::detect_pe_overlap_trim_with(
                s1, s2, pe_anchor, params_r1.adapter.as_deref(), adapter_evidence,
                pe_ambiguous_smallest,
            ) {
                Some((t1, t2)) => {
                    let mut a = Arc::clone(&base1);
                    let mut b = Arc::clone(&base2);
                    Arc::make_mut(&mut a).pe_adapter = Some(t1);
                    Arc::make_mut(&mut b).pe_adapter = Some(t2);
                    (a, b)
                }
                // Detection refused (no overlap, ambiguous, chance hit, or no
                // adapter set to check the boundary against): the mates fall
                // back to the ordinary per-mate behaviour.
                None => (Arc::clone(&base1), Arc::clone(&base2)),
            }
        } else {
            (Arc::clone(&base1), Arc::clone(&base2))
        };

        let out1 = trim(s1, q1, &p1);
        let out2 = trim(s2, q2, &p2);

        // Opt-in posterior PE correction (`--pe-correct`), run strictly AFTER
        // the trim decisions so trimming is byte-identical with the flag off.
        // Full raw reads go in — `infer_insert` and the kernel both map
        // fragment position fpos to s1[fpos] / s2[insert-1-fpos], i.e. 5'-anchored
        // — and only the doubly-retained overlap is editable: positions a mate
        // trimmed away are masked first, so a discarded base can never rewrite a
        // kept one. `corrected` stays `None` (no allocation retained) unless an
        // edit actually fired, and inference refusal leaves output untouched.
        let (corr1, corr2, corrected_bases) = match pe_correct {
            Some((gate, min_ovl))
                if i1.format == RecordFormat::Fastq
                    && i2.format == RecordFormat::Fastq
                    && matches!(out1.result, TrimResult::Emit { .. })
                    && matches!(out2.result, TrimResult::Emit { .. }) =>
            {
                match crate::correct::infer_insert(s1, s2, min_ovl) {
                    Some(ins) => {
                        let keep1 = out1.result.t5()..s1.len().saturating_sub(out1.result.t3());
                        let keep2 = out2.result.t5()..s2.len().saturating_sub(out2.result.t3());
                        let mut c1s = s1.to_vec();
                        let mut c1q = q1.to_vec();
                        let mut c2s = s2.to_vec();
                        let mut c2q = q2.to_vec();
                        let edits = crate::correct::correct_overlap_retained(
                            &mut c1s,
                            &mut c1q,
                            &mut c2s,
                            &mut c2q,
                            keep1,
                            keep2,
                            p1.phred_offset,
                            ins,
                            gate,
                            min_ovl,
                        );
                        if edits > 0 {
                            (Some((c1s, c1q)), Some((c2s, c2q)), edits)
                        } else {
                            (None, None, 0)
                        }
                    }
                    None => (None, None, 0),
                }
            }
            _ => (None, None, 0),
        };
        ProcessedPair {
            p1: Processed {
                input: i1,
                output: out1,
                corrected: corr1,
            },
            p2: Processed {
                input: i2,
                output: out2,
                corrected: corr2,
            },
            corrected_bases,
        }
    };
    if parallel {
        c1.into_par_iter()
            .zip(c2.into_par_iter())
            .map(work)
            .collect()
    } else {
        c1.into_iter().zip(c2).map(work).collect()
    }
}

// ===========================================================================
// Integration tests.
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::MinLen;
    use crate::cli::PeOverlapMode;
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
            dry_run: false,
            quiet: false,
            allow_total_drop: false,
            adapter_source: crate::cli::AdapterSource::None,
            pe_overlap: PeOverlapMode::Never,
            adapter_evidence: crate::adapter::AdapterEvidence::Floor,
            pe_ambiguous_smallest: false,
            pe_correct: None,
            pe_correct_min_overlap: 30,
            adapter_error_rate: 0.1,
            min_adapter_overlap_se: 4,
            dedup: None,
            thread_pool: None,
        }
    }

    fn default_params() -> Params {
        Params::default()
    }

    // ------- Chunk arenas are recycled ------------------------------------

    /// A fresh `ChunkBuf` per chunk was ~464 KB alloc + first-touch per mate per
    /// chunk (928 KB paired), and the old one was only freed *after* the new
    /// one was reserved. Recycling needs a handful for a whole run. The counter
    /// is a thread-local so parallel tests cannot pollute it.
    #[test]
    fn chunk_arenas_are_recycled_not_reallocated() {
        let inp = tmp_path("arena.fq");
        let out = tmp_path("arena-out.fq");
        let mut body = String::new();
        for i in 0..200 {
            body.push_str(&format!(
                "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        write_plain(&inp, &body);
        let rt = Runtime { chunk: 8, ..default_rt() };
        ARENAS.with(|c| c.set(0));
        run_single(&inp, &out, &default_params(), &rt).unwrap();
        let n = ARENAS.with(|c| c.get());
        // 25 chunks; recycling needs one in flight, one being filled, one in
        // the result channel.
        assert!(n <= 3, "{n} chunk arenas allocated for 25 chunks (expected ≤ 3)");
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }

    #[test]
    fn paired_chunk_arenas_are_recycled_too() {
        let in1 = tmp_path("arena-r1.fq");
        let in2 = tmp_path("arena-r2.fq");
        let out1 = tmp_path("arena-o1.fq");
        let out2 = tmp_path("arena-o2.fq");
        let mut body = String::new();
        for i in 0..200 {
            body.push_str(&format!(
                "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        write_plain(&in1, &body);
        write_plain(&in2, &body);
        let p = default_params();
        let rt = Runtime { chunk: 8, ..default_rt() };
        ARENAS.with(|c| c.set(0));
        run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &rt).unwrap();
        let n = ARENAS.with(|c| c.get());
        // Two mates per chunk → two arenas in flight is still a handful.
        assert!(n <= 6, "{n} chunk arenas allocated for 25 paired chunks (expected ≤ 6)");
        for f in [&in1, &in2, &out1, &out2] {
            std::fs::remove_file(f).ok();
        }
    }

    /// `--interior` still works after the per-record `TrimResult` clone was
    /// removed: the runs are read through the borrowed result on the ordered
    /// thread, and the interior mass reaches the counters.
    #[test]
    fn interior_runs_reach_the_counters_without_cloning() {
        let inp = tmp_path("interior.fq");
        let out = tmp_path("interior-out.fq");
        // 16 random | polyA(12) interior | 16 more, ending in T so no 3' run.
        let seq = format!("ACGTACGTACGTACGT{}{}", "A".repeat(12), "CGTACGTACGTACGT");
        write_plain(
            &inp,
            &format!("@r\n{}\n+\n{}\n", seq, "I".repeat(seq.len())),
        );
        let mut p = default_params();
        p.scan_interior = true;
        let rt = default_rt();
        let stats = run_single(&inp, &out, &p, &rt).unwrap();
        assert_eq!(stats.reads_kept, 1);
        assert!(
            stats.bases_interior_mass >= 11,
            "interior mass {} — the runs were lost with the clone",
            stats.bases_interior_mass
        );
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
    }

    // ------- Ingest keeps quality raw ------------------------------------

    /// Quality is not decoded on ingest: the chunk arena holds the input bytes
    /// (clamped up to the offset), the gate LUT is indexed by them, and egress
    /// writes them back. This is the property that makes the whole round trip
    /// disappear; without it the arena holds decoded values and egress has to
    /// re-encode.
    #[test]
    fn ingest_keeps_raw_quality_bytes_and_clamps_below_the_offset() {
        let fq = b"@r\nACGT\n+\n I#I\n".to_vec(); // quals: ' ', 'I', '#', 'I'
        let mut parser = needletail::parse_fastx_reader(std::io::Cursor::new(fq)).unwrap();
        let rec = parser.next().unwrap().unwrap();
        let mut buf = ChunkBuf::default();
        let si = ScanInput::push_needletail(&mut buf, &rec, 33).unwrap();
        assert_eq!(
            &buf.qual[si.qual.clone()],
            b"!I#I",
            "ingest must store the input bytes, clamped up to the offset"
        );
        // And a second record appends rather than overwriting.
        let si2 = ScanInput::push_needletail(&mut buf, &rec, 33).unwrap();
        assert_eq!(&buf.qual[si2.qual.clone()], b"!I#I");
        assert_eq!(si2.qual.start, 4);
    }

    /// End to end: a FASTQ round trip must return the quality bytes it was
    /// given, including a sub-offset byte, which is clamped exactly as the old
    /// decode/encode path clamped it.
    #[test]
    fn quality_round_trips_through_the_pipeline() {
        let inp = tmp_path("q-rt.fq");
        let out = tmp_path("q-rt-out.fq");
        write_plain(
            &inp,
            "@r1\nACGTACGTACGTACGTACGTACGT\n+\n I#I5?+ABCDEFGHIJKLMNOPQ\n",
        );
        let mut rt = default_rt();
        rt.phred_offset = 33;
        let stats = run_single(&inp, &out, &default_params(), &rt).unwrap();
        assert_eq!(stats.reads_kept, 1);
        let body = read_plain(&out);
        assert!(
            body.contains("\n!I#I5?+ABCDEFGHIJKLMNOPQ\n"),
            "quality bytes changed: {body:?}"
        );
        std::fs::remove_file(&inp).ok();
        std::fs::remove_file(&out).ok();
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
        // This fixture is a total drop by design (R1 doomed, R2 discarded):
        // opt out of the refusal, which `run_paired` now raises BEFORE the
        // rename, so the counters can be inspected.
        let mut rt = default_rt();
        rt.allow_total_drop = true;
        let stats = run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &rt).unwrap();
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
            dry_run: false,
            quiet: false,
            allow_total_drop: false,
            adapter_source: crate::cli::AdapterSource::None,
            pe_overlap: PeOverlapMode::Never,
            adapter_evidence: crate::adapter::AdapterEvidence::Floor,
            pe_ambiguous_smallest: false,
            pe_correct: None,
            pe_correct_min_overlap: 30,
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
            dry_run: false,
            quiet: false,
            allow_total_drop: false,
            adapter_source: crate::cli::AdapterSource::None,
            pe_overlap: PeOverlapMode::Never,
            adapter_evidence: crate::adapter::AdapterEvidence::Floor,
            pe_ambiguous_smallest: false,
            pe_correct: None,
            pe_correct_min_overlap: 30,
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

    // ------- Phred-offset auto-detection ---------------------------------

    #[test]
    fn auto_phred_high_quality_only_is_33_not_64() {
        // Regression: a Phred+33 file whose lowest quality is Q35 ('D' = 68) has
        // every byte at or above 64. The old rule `min < 64 ? 33 : 64` read that
        // as Phred+64, subtracting 31 from every quality — and a --mean-q filter
        // then dropped the entire file. The absence of a sub-64 byte is not
        // evidence for +64.
        let inp = tmp_path("hq.fq");
        write_plain(&inp, "@r\nACGT\n+\nDDDD\n@r\nACGT\n+\nDDDD\n");
        let census = detect_phred_offset(&inp).unwrap();
        assert_eq!(census.offset, 33);
        assert!(!census.proven, "'D' = 68 proves nothing");
        std::fs::remove_file(&inp).ok();
    }

    #[test]
    fn auto_phred_never_asserts_64() {
        // A byte above 74 used to be taken as evidence for Phred+64. It is not:
        // under +33 the same byte is Q>41, which current chemistry produces
        // routinely, and believing +64 subtracts 31 from every quality — with
        // --mean-q on, that deletes the whole file (5000/5000 reads, exit 0,
        // reproduced 2026-10-20). The census cannot distinguish the two
        // encodings, so the auto path reports the modern default and says so;
        // --phred-offset 64 is the way to ask for legacy data.
        let inp = tmp_path("p64.fq");
        write_plain(&inp, "@r\nACGT\n+\nhhhh\n@r\nACGT\n+\nhhhh\n"); // 'h' = 104
        let census = detect_phred_offset(&inp).unwrap();
        assert_eq!(census.offset, 33);
        assert!(!census.proven);
        let a = census.assumption();
        assert!(
            a.basis.contains("cannot distinguish"),
            "the basis must admit the census cannot decide: {}",
            a.basis
        );
        assert!(a.consequence.contains("31"));
        std::fs::remove_file(&inp).ok();
    }

    #[test]
    fn auto_phred_sub_64_byte_proves_33() {
        // A byte below 64 cannot occur under Phred+64, so its presence proves +33.
        let inp = tmp_path("mix.fq");
        write_plain(&inp, "@r\nACGT\n+\n#DDD\n@r\nACGT\n+\nDDDD\n"); // '#' = 35
        let census = detect_phred_offset(&inp).unwrap();
        assert_eq!(census.offset, 33);
        assert!(census.proven, "a sub-64 byte proves +33");
        std::fs::remove_file(&inp).ok();
    }
}
