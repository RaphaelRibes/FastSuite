//! Adapter detection and trimming.
//!
//! **Sequence matching** is the only mode: slide each adapter over the read,
//! accept the largest overlap that passes the mismatch-rate threshold, and
//! trim from the match position. O(L × A) per read; negligible for Illumina
//! read/adapter lengths. Mates are trimmed independently — there is no
//! cross-mate inference — so PE overlap detection (which once compared the
//! R1 suffix with the reverse-complement R2 prefix to infer the insert) has
//! been removed; the adapter sequence comes from the configured set and
//! nowhere else.
//!
//! Auto-detection samples the first 10 k reads, counts catalogue 12-mers at
//! every offset (not just the 3′ end — a read that runs through the adapter into
//! the primer has no adapter at its terminus), and uses the top candidates as
//! adapters for subsequent sequence matching.

use memchr::memmem;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};

// ---------------------------------------------------------------------------
// Core types.
// ---------------------------------------------------------------------------

/// One named adapter sequence.
///
/// `PartialEq`/`Eq` are needed by [`Detection`] and by the dedup checks in
/// [`catalogue`] and [`identify_from_reads`], which compare by sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    /// Human-readable name used in stats and (optionally) metadata.
    /// Preset adapters use fixed strings, FASTA-loaded adapters use the
    /// sequence ID, and auto-detected adapters take the catalogue name of
    /// whichever sequence their fingerprint identified.
    pub name: String,
    pub seq: Vec<u8>,
}

/// The full adapter configuration passed to [`match_adapters`].
#[derive(Debug, Clone)]
pub struct AdapterSet {
    pub adapters: Vec<Adapter>,
    /// Maximum fraction of mismatches allowed in a match.
    pub error_rate: f64,
    /// Minimum overlap for sequence matching.
    pub min_overlap_se: usize,
    /// Seed index, built once at construction. `match_adapters` uses it to
    /// skip the ~118-position full scan on reads that share no seed with any
    /// adapter (the pigeonhole argument above makes this exact, not a
    /// heuristic). The seed width is derived from the error rate, so this is
    /// not a fixed 6-mer table.
    pub seeds: Arc<SeedTable>,
}

impl AdapterSet {
    /// Build a set and its seed index. The index is only a rejection aid
    /// (exact, never a heuristic), so it can always be built.
    pub fn new(adapters: Vec<Adapter>, error_rate: f64, min_overlap_se: usize) -> Self {
        let seeds = Arc::new(SeedTable::build(&adapters, error_rate));
        AdapterSet {
            adapters,
            error_rate,
            min_overlap_se,
            seeds,
        }
    }
}

/// Seed index for fast rejection of adapter-free reads.
///
/// Soundness of the pigeonhole argument, stated exactly: an alignment of `k`
/// bases with `mm = ceil(k·ε)` mismatches has at most `mm + 1` exactly
/// matching runs, totalling `k − mm` bases, so the longest of them is at least
/// `floor(k / (mm + 1))`. When that is `>= s` — the seed width — the alignment
/// contains an adapter `s`-mer that appears verbatim in the read, so a seed
/// lookup at that position finds the alignment and everything below `s` bases
/// of run can only be found by the direct scan.
///
/// The seed width is therefore a derived quantity, not a constant: `s = 6`
/// requires `k >= 6·(mm+1)`, which has no solution once `ε >= 1/6`, so at the
/// shipped `--adapter-error-rate 0.2` a fixed 6-mer table declared every
/// adapter unsound (`min_sound_k = usize::MAX`) and every read took the
/// reference full scan — 414 exact compares at L = 151 while only 2.3 % of
/// reads carry an adapter. Sizing `s` from the bound, per adapter, revives the
/// filter where it was dead (81 compares on random sequence, measured) and
/// leaves it untouched where it already worked.
#[derive(Debug)]
pub struct SeedTable {
    /// `table[index]` = index into `pairs` of the first `(adapter, offset)`
    /// for this seed code; `-1` when empty. Flat and ~11 KB so the per-read
    /// lookup stays in cache, unlike the old `Vec<Vec<…>>` which
    /// pointer-chased on every position of every read (profile: ~49% of
    /// runtime after the first seed optimisation).
    ///
    /// One array holds every width: each width owns a disjoint code range
    /// ([`code_base`]), because a 4-mer code and a 5-mer code with the same
    /// numeric value are different seeds and must not share a bucket.
    table: [i16; TABLE_LEN],
    /// Flat `(tagged code, adapter index, offset within adapter)` list, grouped
    /// by code so a lookup iterates a contiguous run.
    pairs: Vec<(u16, u8, u8)>,
    /// Smallest overlap length for which the seed guarantee holds, per
    /// adapter (`usize::MAX` = never: the adapter is too short to host the
    /// bound at this error rate and falls back to a direct scan).
    min_sound_k: Vec<usize>,
    /// Seed width chosen per adapter, `0` when the adapter is not seeded.
    width: Vec<usize>,
    /// The distinct widths in use, ascending. The scan iterates exactly these
    /// per position — normally one or two, because the wide adapters and a
    /// short one need different widths at the same error rate.
    widths: Vec<usize>,
}

/// Candidate seed widths, widest first.
///
/// A wider seed is exponentially more selective — a 4-mer is shared by almost
/// every read, a 6-mer by almost none — so each adapter gets the widest width
/// its own length and the error rate allow, and only shrinks when the
/// arithmetic forces it. At the shipped `--adapter-error-rate 0.2` the 33-base
/// TruSeq adapters and the 19-base Nextera adapter both take a 4-mer: with the
/// floor tolerance the pigeonhole bound holds from k = 16 for either length.
/// With the old fixed 6 every one of them was unsound and took the reference
/// full scan.
const SEED_WIDTHS: [usize; 4] = [6, 5, 4, 3];

/// Code ranges: width `s` owns `[code_base(s), code_base(s) + 4^s)`.
const TABLE_LEN: usize = 5440;

/// Base index of width `s`'s code range. Width 6 keeps the low range so the
/// original 12-bit table layout is unchanged for existing error rates.
#[inline]
fn code_base(s: usize) -> usize {
    match s {
        6 => 0,
        5 => 4096,
        4 => 5120,
        3 => 5376,
        _ => unreachable!("seed width is one of 3..=6"),
    }
}

#[inline]
fn code_mask(s: usize) -> u32 {
    (1u32 << (2 * s)) - 1
}

/// Smallest `k` at which an alignment of `k` bases with `floor(k·ε)`
/// mismatches is guaranteed to contain an exact run of at least `s` bases, or
/// `usize::MAX` when no `k <= a_len` satisfies it.
///
/// The test is `k >= s * (mm + 1)` — the pigeonhole bound above, in integers —
/// not `ceil(k / (mm + 1)) >= s`: the latter is optimistic by one seed width,
/// because the mismatch-free chunk may be one of the short ones when `k` is not
/// a multiple of `mm + 1`. Using it would make the filter unsound, which the
/// conformance test catches as a changed trim, not as a slow path.
///
/// The predicate is **not monotone in `k`**: `floor(k·ε)` steps up at
/// `k = m/ε`, and each step raises the right-hand side by `s` while `k` rose
/// by one, so a `k` that satisfies it can be followed by one that does not.
/// Returning the first satisfying `k` therefore claimed coverage for later,
/// unsound overlaps — at the shipped `s = 4, ε = 0.2` the rounds-up form was
/// false again at 21, 22, 23, 26, 27 and 31 — and those positions were then
/// skipped by both passes. The threshold returned here is instead the first
/// `k` after the LAST failing one, so every `k` the seed pass claims really is
/// covered and everything below falls into the direct scan.
fn min_sound_k_for(a_len: usize, error_rate: f64, s: usize) -> usize {
    if a_len < s {
        return usize::MAX;
    }
    let mut last_unsound = None;
    for k in s..=a_len {
        if !seed_bound_holds(k, error_rate, s) {
            last_unsound = Some(k);
        }
    }
    match last_unsound {
        // Nothing above the last failure can satisfy the bound either.
        Some(k) if k == a_len => usize::MAX,
        Some(k) => (k + 1).max(s),
        None => s,
    }
}

/// The pigeonhole predicate for one overlap length: with
/// `mm = max_mismatches(k, ε)` mismatches the alignment still contains an
/// exact run of at least `s` bases.
#[inline]
fn seed_bound_holds(k: usize, error_rate: f64, s: usize) -> bool {
    k >= s * (max_mismatches(k, error_rate) + 1)
}

/// Mismatches tolerated in an overlap of `k` bases: `floor(k·ε)`.
///
/// Rounded DOWN, never up. `ceil` let `--adapter-error-rate 0.2` admit 4
/// mismatches in 16 bases (25%) and 3 in 14 (21.4%), which is above the rate
/// the flag documents and multiplied the false-trim rate; the floor is what
/// the flag says and what the seed bound is derived from.
#[inline]
fn max_mismatches(k: usize, error_rate: f64) -> usize {
    (k as f64 * error_rate) as usize
}

impl SeedTable {
    fn build(adapters: &[Adapter], error_rate: f64) -> Self {
        // Per adapter, widest sound width. `0` (no width) means the adapter
        // cannot host any bound at this error rate; it then takes the direct
        // scan for every position, exactly as before.
        let width: Vec<usize> = adapters
            .iter()
            .map(|a| {
                SEED_WIDTHS
                    .iter()
                    .copied()
                    .find(|&s| {
                        a.seq.len() >= s && min_sound_k_for(a.seq.len(), error_rate, s) != usize::MAX
                    })
                    .unwrap_or(0)
            })
            .collect();
        let mut widths: Vec<usize> = width.iter().copied().filter(|&s| s > 0).collect();
        widths.sort_unstable();
        widths.dedup();

        let mut pairs: Vec<(u16, u8, u8)> = Vec::new();
        for ((ai, a), &s) in adapters.iter().enumerate().zip(width.iter()) {
            if s == 0 {
                continue;
            }
            let base = code_base(s);
            let mask = code_mask(s);
            for off in 0..=a.seq.len() - s {
                let code = seed_code(&a.seq[off..off + s], mask);
                pairs.push(((base + code as usize) as u16, ai as u8, off as u8));
            }
        }
        // Group by code so a lookup is a contiguous run. This matters: two
        // adapters (e.g. TruSeq R1/R2) share a long prefix, so the same seed
        // appears in both at the same offset — without the sort, table[code]
        // would point into R1's block and the run would stop before R2's
        // entry, silently missing valid R2 alignments.
        pairs.sort_unstable_by_key(|p| p.0);
        let mut table = [-1i16; TABLE_LEN];
        for (i, &(code, _, _)) in pairs.iter().enumerate() {
            if i == 0 || pairs[i - 1].0 != code {
                table[code as usize] = i as i16;
            }
        }
        let min_sound_k = adapters
            .iter()
            .zip(width.iter())
            .map(|(a, &s)| {
                if s == 0 {
                    usize::MAX
                } else {
                    min_sound_k_for(a.seq.len(), error_rate, s)
                }
            })
            .collect();
        SeedTable {
            table,
            pairs,
            min_sound_k,
            width,
            widths,
        }
    }
}

/// 2-bit encoding of a seed into a `2 * bases.len()`-bit code (A=0, C=1, G=2,
/// T=3), masked to the table's code space.
#[inline]
fn seed_code(bases: &[u8], mask: u32) -> u32 {
    let mut code: u32 = 0;
    for &b in bases {
        let c = match b {
            b'A' | b'a' => 0u32,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            _ => 3, // T, or any non-ACGT base (wildcard: may spuriously hit, rejected later)
        };
        code = (code << 2) | c;
    }
    code & mask
}

/// Byte → 2-bit base code lookup (0=A, 1=C, 2=G, 3=T; anything else → 3).
/// Used by the rolling seed scan to avoid a 4-way branch per base; the table
/// stays in L1 and the compiler turns the loop body into a couple of loads
/// and shifts.
static BASE2: [u8; 256] = {
    let mut t = [3u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t
};

/// Index reserved for the PE-overlap pseudo-adapter in [`AdapterMatch`] and in
/// the per-adapter counters.
///
/// The overlap detector pins a boundary without naming an adapter, and the
/// report still has to attribute the trim, so it owns the top index. Real
/// adapters are `0..adapters.len()`, which no set can push past this.
pub const PE_OVERLAP_ADAPTER: u8 = u8::MAX;

/// Result of a successful adapter match on one read.
///
/// The match carries an INDEX into the configured set, not the adapter's name:
/// a name is a heap `String`, and this value crosses a thread boundary and is
/// then counted for every trimmed read (462 792 of them at 5 M pairs). The name
/// is resolved once, at report time, from [`crate::algorithm::Params::adapter`]
/// and [`PE_OVERLAP_ADAPTER`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterMatch {
    /// Index into `AdapterSet::adapters`, or [`PE_OVERLAP_ADAPTER`].
    pub adapter: u8,
    /// Bases to trim from the 5′ end of this read.
    pub t5: usize,
    /// Bases to trim from the 3′ end of this read.
    pub t3: usize,
}

// ---------------------------------------------------------------------------
// Built-in presets.
// ---------------------------------------------------------------------------

/// All TruSeq + Nextera sequences (default "illumina" preset).
pub fn preset_illumina() -> Vec<Adapter> {
    vec![
        Adapter { name: "TruSeq/R1".into(),  seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec() },
        Adapter { name: "TruSeq/R2".into(),  seq: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec() },
        Adapter { name: "Nextera".into(),     seq: b"CTGTCTCTTATACACATCT".to_vec() },
    ]
}

pub fn preset_nextera() -> Vec<Adapter> {
    vec![
        Adapter { name: "Nextera".into(), seq: b"CTGTCTCTTATACACATCT".to_vec() },
    ]
}

pub fn preset_small_rna() -> Vec<Adapter> {
    vec![
        Adapter { name: "SmallRNA".into(), seq: b"TGGAATTCTCGGGTGCCAAGG".to_vec() },
    ]
}

/// Every built-in sequence, deduplicated. Turns a detected k-mer
/// fingerprint into a full-length adapter.
pub fn catalogue() -> Vec<Adapter> {
    let mut all = preset_illumina();
    for a in preset_nextera().into_iter().chain(preset_small_rna()) {
        if !all.iter().any(|x| x.seq == a.seq) {
            all.push(a);
        }
    }
    all
}

// ---------------------------------------------------------------------------
// Auto-detection: the k-mer identifies, the catalogue supplies the sequence.
//
// Detection previously handed the raw k-mer straight to the trimmer, which
// was wrong twice over. It is 12 bases, and matching a 12-mer within 2
// mismatches fires on ~1.6% of clean reads by chance. Worse, it is the read's
// *terminal* k-mer, so it is the END of whatever adapter fragment is present
// rather than its start: at an insert of 126 on 150-base reads it is
// `adapter[12:24]`, and trimming from its match position removed 12 of the 24
// adapter bases and left the other 12 behind. It only lined up when the
// insert happened to be exactly L - k.
//
// That first fix matched the fingerprint anywhere *inside the catalogue*, but
// still only COUNTED the terminal k-mer — which fails on the ordinary case of a
// read that runs through the adapter and continues into the downstream primer.
// Measured on ERR10123282 (MiSeq 2x251, Nextera): 26.9% of reads carry the full
// 19-base adapter, but only 3.05% end inside it, and those are spread over ~19
// distinct terminal k-mers (one per insert length), so no single one reaches the
// 1% threshold and nothing is detected. The library with the most readthrough
// was the one detection could not see.
//
// So fingerprints are now counted at EVERY offset, restricted to 12-mers that
// occur in the catalogue. Restricting to the catalogue keeps the counter bounded
// by the catalogue rather than by the read count (16.7M possible 12-mers would
// otherwise need a large table) and admits no non-catalogue noise, which is what
// the frequency threshold is guarding against.
// ---------------------------------------------------------------------------

/// Length of the k-mer used as a fingerprint. It only has to name an
/// adapter unambiguously, not to trim with, so 12 is ample: a given 12-mer
/// occurs in random sequence about once per 1.7e7 bases.
const DETECT_K: usize = 12;

/// A fingerprint must appear in at least this fraction of sampled reads.
const DETECT_MIN_FREQ: f64 = 0.01;

/// Distinct fingerprints counted per read. An adapter normally contributes one
/// or two; the cap bounds the work on a pathological read (a tandem repeat of a
/// catalogue sequence) without needing an allocation per read.
const DETECT_MAX_HITS_PER_READ: usize = 8;

/// What auto-detection concluded.
///
/// A run commits every read to a decision taken from one sample, so which
/// adapter was chosen and on what evidence must be visible rather than
/// inferred from the trim counts afterwards.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Detection {
    pub adapters: Vec<Adapter>,
    /// `(fingerprint, reads containing it)` for each k-mer that hit the
    /// catalogue, most frequent first.
    pub evidence: Vec<(Vec<u8>, usize)>,
    pub sampled: usize,
}

/// Every [`DETECT_K`]-mer occurring anywhere in the catalogue.
fn catalogue_fingerprints(cat: &[Adapter]) -> std::collections::HashSet<[u8; DETECT_K]> {
    let mut set = std::collections::HashSet::new();
    for a in cat {
        if a.seq.len() < DETECT_K {
            continue;
        }
        for w in a.seq.windows(DETECT_K) {
            set.insert(w.try_into().unwrap());
        }
    }
    set
}

/// Identify adapters by fingerprinting sampled reads against [`catalogue`].
///
/// A read contributes each catalogue fingerprint it contains **once**, at any
/// offset, so the evidence count stays "reads containing this fingerprint" —
/// which is what the caller reports as a percentage of the sample.
///
/// Returns empty `adapters` when nothing is recognised; the caller decides
/// whether that is fatal. Never guesses a sequence.
pub fn identify_from_reads(reads: &[&[u8]]) -> Detection {
    let mut det = Detection { adapters: vec![], evidence: vec![], sampled: reads.len() };
    if reads.is_empty() {
        return det;
    }

    let cat = catalogue();
    let fingerprints = catalogue_fingerprints(&cat);
    if fingerprints.is_empty() {
        return det;
    }

    let mut counts: HashMap<[u8; DETECT_K], usize> = HashMap::new();
    for read in reads {
        if read.len() < DETECT_K {
            continue;
        }
        let mut hits = [[0u8; DETECT_K]; DETECT_MAX_HITS_PER_READ];
        let mut n = 0usize;
        for w in read.windows(DETECT_K) {
            let kmer: [u8; DETECT_K] = w.try_into().unwrap();
            if !fingerprints.contains(&kmer) || hits[..n].contains(&kmer) {
                continue;
            }
            if n == DETECT_MAX_HITS_PER_READ {
                break;
            }
            hits[n] = kmer;
            n += 1;
        }
        for kmer in &hits[..n] {
            *counts.entry(*kmer).or_insert(0) += 1;
        }
    }

    let min_count = ((reads.len() as f64) * DETECT_MIN_FREQ).ceil() as usize;
    let mut candidates: Vec<([u8; DETECT_K], usize)> =
        counts.into_iter().filter(|(_, c)| *c >= min_count).collect();
    candidates.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    for (kmer, count) in candidates {
        // The fingerprint sits at an offset inside the adapter that depends
        // on the insert size, so match anywhere in the sequence, not at its
        // start.
        let hits: Vec<&Adapter> =
            cat.iter().filter(|a| memmem::find(&a.seq, &kmer).is_some()).collect();
        if hits.is_empty() {
            continue;
        }
        det.evidence.push((kmer.to_vec(), count));
        for h in hits {
            if !det.adapters.iter().any(|x| x.seq == h.seq) {
                det.adapters.push(h.clone());
            }
        }
    }
    det
}

// ---------------------------------------------------------------------------
// FASTA loader.
// ---------------------------------------------------------------------------

/// Load adapter sequences from a FASTA file. Sequence IDs become adapter
/// names. Non-unique IDs are accepted (counts are merged in stats).
pub fn load_fasta(path: &Path) -> Result<Vec<Adapter>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening adapter FASTA {}", path.display()))?;
    let reader = std::io::BufReader::new(file);

    let mut adapters: Vec<Adapter> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.context("reading adapter FASTA")?;
        let line = line.trim();
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                if !current_seq.is_empty() {
                    adapters.push(Adapter { name, seq: current_seq.clone() });
                    current_seq.clear();
                }
            }
            current_name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else if current_name.is_some() {
            current_seq.extend(line.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(name) = current_name {
        if !current_seq.is_empty() {
            adapters.push(Adapter { name, seq: current_seq });
        }
    }

    if adapters.is_empty() {
        anyhow::bail!("adapter FASTA {} contained no sequences", path.display());
    }
    Ok(adapters)
}

// ---------------------------------------------------------------------------
// Sequence matching.
// ---------------------------------------------------------------------------

/// Find the leftmost adapter match in `seq`, or `None`.
///
/// The scan is seed-guided: the read is walked once, a rolling `seed_len`-base
/// code is looked up in the adapter seed table, and only positions where the
/// read actually contains an adapter seed are tested with the exact
/// [`mismatches_upto`] comparison. This is *exact*, not a heuristic:
///
/// * any valid alignment of `k` bases with `mm = floor(k·ε)` mismatches has at
///   most `mm + 1` exactly-matching runs totalling `k − mm` bases, so one of
///   them is at least `floor(k / (mm + 1))` long;
/// * when that is `>= seed_len` (i.e. `k >= seed_len·(mm+1)`, exactly the
///   `min_sound_k` bound) it contains an adapter seed appearing verbatim in
///   the read, so the seed hit at its position finds the alignment;
/// * positions whose overlap `k` falls below `min_sound_k` (short truncated
///   tails, or adapters too short for the bound at this error rate) are
///   scanned directly.
///
/// "Leftmost" = smallest `p`, tie-broken by first adapter in set order, which
/// is what the reference implementation produces (largest `t3 = l − p`).
pub fn match_adapters(seq: &[u8], set: &AdapterSet) -> Option<AdapterMatch> {
    let l = seq.len();
    let min_overlap = set.min_overlap_se;
    let max_p = l.saturating_sub(min_overlap);

    // Global best: (p, adapter index). Smaller p wins; on a tie the adapter
    // that appears first in the set wins (matches the reference semantics).
    let mut best: Option<(usize, usize)> = None;
    let mut consider = |p: usize, ai: usize| {
        if best.is_none_or(|(bp, bai)| p < bp || (p == bp && ai < bai)) {
            best = Some((p, ai));
        }
    };

    // ---- Seed-guided pass (sound for k >= min_sound_k) -------------------
    // Walk the read once, rolling a 12-bit 6-mer code and deriving each active
    // seed width from it by masking; every table hit yields candidate
    // alignment starts p = q + 1 - s - offset for the adapter at that offset.
    // Only positions that provably need checking are tested.
    if !set.seeds.widths.is_empty() {
        let mut code6: u32 = 0;
        for q in 0..l {
            let c = BASE2[seq[q] as usize] as u32;
            code6 = ((code6 << 2) | c) & 0xFFF;
            for &s in &set.seeds.widths {
                if q + 1 < s {
                    continue; // window not full yet
                }
                let idx = (code_base(s) + (code6 & code_mask(s)) as usize) as u16;
                let start = set.seeds.table[idx as usize];
                if start < 0 {
                    continue;
                }
                let win_start = q + 1 - s; // start of the seed ending at q
                for &(pc, ai, off) in &set.seeds.pairs[start as usize..] {
                    if pc != idx {
                        break; // contiguous run for this code ends
                    }
                    let a = set.adapters[ai as usize].seq.len();
                    let ms = set.seeds.min_sound_k[ai as usize];
                    // Soundness requires k >= ms; p must be in range.
                    if win_start < off as usize {
                        continue; // p would be negative
                    }
                    let p = win_start - off as usize;
                    let k = a.min(l - p);
                    if k < min_overlap || k < ms || p > max_p {
                        continue;
                    }
                    let max_mm = max_mismatches(k, set.error_rate);
                    if mismatches_upto(&seq[p..p + k], &set.adapters[ai as usize].seq[..k], max_mm)
                        <= max_mm
                    {
                        consider(p, ai as usize);
                    }
                }
            }
        }
    }

    // ---- Direct pass for the seed-unsound region -------------------------
    // Positions with k = min(a, l-p) < min_sound_k: short truncated tails
    // (p > l - ms) and, for short adapters or high error rates, whole
    // adapters (a < ms). Scanned exactly like the reference implementation.
    for (ai, adapter) in set.adapters.iter().enumerate() {
        let a = adapter.seq.len();
        if a == 0 {
            continue;
        }
        let ms = set.seeds.min_sound_k[ai];
        if set.seeds.width[ai] == 0 {
            // This adapter could not host the bound at this error rate (too
            // short, or eps too high): every position is scanned directly,
            // reference-style. `ms` is `usize::MAX` here, so the seed pass
            // contributed nothing for it.
            for p in 0..=max_p {
                let k = a.min(l - p);
                if k < min_overlap {
                    break;
                }
                let max_mm = max_mismatches(k, set.error_rate);
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                    break; // leftmost match for this adapter
                }
            }
        } else {
            // Tail region only: p > l - ms  ⟺  p >= l - ms + 1. For short
            // reads (l < ms) this covers every position, so the range is
            // non-empty whenever the loop body can run.
            let regime_start = if l >= a { l - a + 1 } else { 0 };
            let start = regime_start.max(l.saturating_add(1).saturating_sub(ms));
            if start <= max_p {
                for p in start..=max_p {
                    let k = a.min(l - p);
                    if k < min_overlap {
                        break;
                    }
                    let max_mm = max_mismatches(k, set.error_rate);
                    if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                        consider(p, ai);
                        break; // leftmost match for this adapter
                    }
                }
            }
        }
    }

    best.map(|(p, ai)| AdapterMatch {
        adapter: ai as u8,
        t5: 0,
        t3: l - p,
    })
}

/// Mismatches between two equal-length slices, stopping once `cap` is
/// exceeded.
///
/// Counted in 32-byte blocks: the block body is a branchless compare-and-sum
/// that vectorises, while the saturation test sits between blocks rather
/// than inside the hot loop. The previous per-byte form had an unpredictable
/// branch on every base and was the single most expensive thing in the
/// program.
///
/// The block body is written as an explicit byte loop rather than
/// `iter().zip().filter(|(x, y)| x != y).count()`: the iterator closure chain
/// defeats LLVM's vectoriser (measured ~21% of total runtime in the adapter
/// scan), while the plain indexed loop vectorises with the portable baseline
/// ISA (measured ~40% faster than the closure chain). No `target-cpu` flag is
/// set: shipped binaries must run on any x86-64 machine.
#[inline]
fn mismatches_upto(a: &[u8], b: &[u8], cap: usize) -> usize {
    debug_assert_eq!(a.len(), b.len());
    #[cfg(test)]
    MM_CALLS.with(|c| c.set(c.get() + 1));
    let mut mm = 0usize;
    for (ca, cb) in a.chunks(32).zip(b.chunks(32)) {
        for i in 0..ca.len() {
            // Branchless compare: a differing byte adds 1. The compiler
            // turns this into a compare + pcmpeqb/pmovmskb-style popcount
            // over the 32-byte block.
            mm += (ca[i] != cb[i]) as usize;
        }
        if mm > cap {
            return mm;
        }
    }
    mm
}

// Counts exact compares. Test builds only, so the shipped binary pays nothing
// for the instrument — without it the prefilter's value is invisible: a dead
// prefilter and a live one produce the same trims, just at 414 calls per read
// instead of a few dozen.
#[cfg(test)]
thread_local! {
    static MM_CALLS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Reset the exact-compare counter, run `f`, and return how many compares it
/// took. `mismatches_upto` is otherwise invisible to a test.
#[cfg(test)]
fn count_mismatch_calls<R>(f: impl FnOnce() -> R) -> u64 {
    MM_CALLS.with(|c| c.set(0));
    let _ = f();
    MM_CALLS.with(|c| c.get())
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_set(seq: &[u8]) -> AdapterSet {
        AdapterSet::new(
            vec![Adapter {
                name: "test".into(),
                seq: seq.to_vec(),
            }],
            0.1,
            4,
        )
    }

    #[test]
    fn match_exact_3prime() {
        // Read ends with 10 exact adapter bases.
        let adapter = b"AGATCGGAAGAGCACACGTC";
        let mut seq = b"ACGTACGTACGT".to_vec();
        seq.extend_from_slice(&adapter[..10]);
        let set = simple_set(adapter);
        let m = match_adapters(&seq, &set).expect("should match");
        assert_eq!(m.t3, 10);
        assert_eq!(m.t5, 0);
    }

    /// The seed threshold must be monotone-safe: EVERY overlap the seed pass
    /// claims (`k >= min_sound_k`) has to satisfy the pigeonhole predicate, or
    /// that position is skipped by the seed pass (it is above the threshold)
    /// AND by the direct pass (which only covers `k < min_sound_k`).
    ///
    /// The pre-fix form returned the FIRST satisfying `k`; at the shipped
    /// `s = 4`, `eps = 0.2` the rounds-up predicate was false again at 21, 22,
    /// 23, 26, 27 and 31, so a real alignment at those overlaps could be
    /// missed.
    #[test]
    fn seed_bound_is_monotone_safe() {
        for &eps in &[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5] {
            for s in 3..=6 {
                for &a_len in &[19usize, 21, 34, 60] {
                    let ms = min_sound_k_for(a_len, eps, s);
                    if ms == usize::MAX {
                        continue;
                    }
                    assert!(
                        ms >= s && ms <= a_len,
                        "eps={eps} s={s} a={a_len}: threshold {ms} out of range"
                    );
                    for k in ms..=a_len {
                        assert!(
                            seed_bound_holds(k, eps, s),
                            "eps={eps} s={s} a={a_len}: k={k} is claimed seed-covered \
                             but the pigeonhole predicate fails"
                        );
                    }
                    if ms > s {
                        assert!(
                            !seed_bound_holds(ms - 1, eps, s),
                            "eps={eps} s={s} a={a_len}: threshold {ms} is not minimal"
                        );
                    }
                }
            }
        }
    }

    /// The tolerance is the documented rate, floored — never rounded up.
    /// `ceil` admitted 4 mismatches in 16 bases (25%) and 3 in 14 (21.4%) at a
    /// documented 20%.
    #[test]
    fn tolerance_never_exceeds_the_documented_rate() {
        for &eps in &[0.05, 0.1, 0.2, 0.25, 0.3] {
            for k in 1..=100 {
                let mm = max_mismatches(k, eps);
                assert!(
                    mm as f64 <= k as f64 * eps + 1e-9,
                    "eps={eps} k={k}: {mm} mismatches exceeds the documented rate"
                );
            }
        }
        assert_eq!(max_mismatches(16, 0.2), 3, "ceil gave 4 = 25%");
        assert_eq!(max_mismatches(14, 0.2), 2, "ceil gave 3 = 21.4%");
        assert_eq!(max_mismatches(21, 0.2), 4);
        assert_eq!(max_mismatches(30, 0.2), 6);
    }

    /// The implementation as it stood before the pigeonhole filter. The
    /// filter is a pure rejection optimisation, so the two must agree on
    /// every input; anything else is a silent change in what gets trimmed.
    fn match_adapters_reference(seq: &[u8], set: &AdapterSet) -> Option<AdapterMatch> {
        let l = seq.len();
        let mut best: Option<(usize, AdapterMatch)> = None;
        for (idx, adapter) in set.adapters.iter().enumerate() {
            let a = adapter.seq.len();
            if a == 0 {
                continue;
            }
            let max_p = l.saturating_sub(set.min_overlap_se);
            for p in 0..=max_p {
                let k = a.min(l - p);
                if k < set.min_overlap_se {
                    break;
                }
                let max_mm = max_mismatches(k, set.error_rate);
                let mut mm = 0usize;
                let mut matched = true;
                for i in 0..k {
                    if seq[p + i] != adapter.seq[i] {
                        mm += 1;
                        if mm > max_mm {
                            matched = false;
                            break;
                        }
                    }
                }
                if matched {
                    let t3 = l - p;
                    if best.as_ref().is_none_or(|(prev, _)| t3 > *prev) {
                        best = Some((
                            t3,
                            AdapterMatch {
                                adapter: idx as u8,
                                t5: 0,
                                t3,
                            },
                        ));
                    }
                    break;
                }
            }
        }
        best.map(|(_, m)| m)
    }

    #[test]
    fn pigeonhole_filter_matches_reference() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let presets: Vec<Vec<u8>> = vec![
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
            b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(),
            b"CTGTCTCTTATACACATCT".to_vec(),
        ];
        let bases = b"ACGT";
        for case in 0..4000 {
            let l = 5 + (next() % 165) as usize;
            let mut seq: Vec<u8> = (0..l).map(|_| bases[(next() % 4) as usize]).collect();
            // Plant an adapter in three quarters of the cases: as a suffix,
            // internally with sequence after it, and internally with errors.
            match case % 4 {
                0 => {}
                style => {
                    let ad = &presets[(next() % 3) as usize];
                    let mut copy = ad.clone();
                    if style == 3 {
                        for _ in 0..(next() % 5) {
                            let i = (next() as usize) % copy.len();
                            copy[i] = bases[(next() % 4) as usize];
                        }
                    }
                    let p = if l > copy.len() { (next() as usize) % (l - copy.len() + 1) } else { 0 };
                    for (o, &b) in copy.iter().enumerate() {
                        if p + o < seq.len() {
                            seq[p + o] = b;
                        }
                    }
                    if style == 1 {
                        seq.truncate((p + copy.len()).min(l));
                    }
                }
            }
            for &eps in &[0.0f64, 0.05, 0.1, 0.15, 0.2, 0.3] {
                for &ov in &[3usize, 6, 10, 20] {
                    let set = AdapterSet::new(
                        presets.iter().enumerate().map(|(i, s)| Adapter {
                            name: format!("a{i}"),
                            seq: s.clone(),
                        }).collect(),
                        eps,
                        ov,
                    );
                    assert_eq!(
                        match_adapters(&seq, &set),
                        match_adapters_reference(&seq, &set),
                        "case {case} eps {eps} overlap {ov} seq {}",
                        String::from_utf8_lossy(&seq),
                    );
                }
            }
        }
    }

    #[test]
    fn match_adapter_fully_internal() {
        // Simulates a short-insert read: insert = 12 bases, then full adapter
        // (20bp), then 18bp of primer junk. The adapter is fully inside the
        // read — the old tail-only algorithm would miss it.
        let adapter = b"AGATCGGAAGAGCACACGTC"; // 20bp
        let mut seq = b"ACGTACGTACGT".to_vec(); // 12bp insert
        seq.extend_from_slice(adapter);          // 20bp adapter
        seq.extend_from_slice(b"NNNNNNNNNNNNNNNNNN"); // 18bp junk
        // Total: 50bp. Adapter starts at position 12 → t3 = 50 - 12 = 38.
        let set = simple_set(adapter);
        let m = match_adapters(&seq, &set).expect("should find internal adapter");
        assert_eq!(m.t3, 38, "should trim from adapter start, not just the tail");
        assert_eq!(m.t5, 0);
    }

    #[test]
    fn match_with_one_mismatch() {
        let adapter = b"AGATCGGAAGAG";
        // Last 8 bases with 1 mismatch. floor(8 * 0.2) = 1 at the shipped
        // rate, so this is a match. (Before the tolerance was floored, the
        // fixture ran at 0.1 and relied on `ceil(0.8) = 1` — 1/8 = 12.5%,
        // above the rate the flag documents.)
        let mut seq = b"ACGTACGTACGTACGT".to_vec();
        seq.extend_from_slice(b"AGATCGGX"); // X is mismatch
        let set = AdapterSet::new(
            vec![Adapter {
                name: "test".into(),
                seq: adapter.to_vec(),
            }],
            0.2,
            4,
        );
        let m = match_adapters(&seq, &set).expect("should match with 1 mm");
        assert_eq!(m.t3, 8);
    }

    #[test]
    fn no_match_below_min_overlap() {
        // Only 3 adapter bases at end, min_overlap = 4.
        let adapter = b"AGATCGGAAGAG";
        let mut seq = b"ACGTACGTACGTACGTACGT".to_vec();
        seq.extend_from_slice(&adapter[..3]);
        let set = simple_set(adapter);
        assert!(match_adapters(&seq, &set).is_none());
    }

    #[test]
    fn no_match_too_many_mismatches() {
        let adapter = b"AGATCGGAAGAG";
        // 8 bases, all wrong.
        let mut seq = b"ACGTACGTACGT".to_vec();
        seq.extend_from_slice(b"TTTTTTTT");
        let set = simple_set(adapter);
        assert!(match_adapters(&seq, &set).is_none());
    }

    /// Two adapters sharing a long prefix: a read carrying the *second*
    /// adapter must be matched against the second one even though every seed
    /// 6-mer also appears in the first. This is the seed-table grouping
    /// regression: if entries are grouped per-adapter instead of per-code,
    /// the contiguous-run lookup stops at the first adapter's entry and the
    /// second adapter's identical-offset entries are never seen.
    #[test]
    fn shared_prefix_second_adapter_still_matches() {
        let a1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(); // TruSeq R1
        let a2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(); // TruSeq R2
        let set = AdapterSet::new(
            vec![
                Adapter {
                    name: "a1".into(),
                    seq: a1.clone(),
                },
                Adapter {
                    name: "a2".into(),
                    seq: a2.clone(),
                },
            ],
            0.1,
            4,
        );
        // Read = 100 random bases + a2 (full-length, no errors).
        let mut seq = b"ACGT".repeat(25);
        seq.extend_from_slice(&a2);
        let m = match_adapters(&seq, &set).expect("should match a2");
        assert_eq!(m.adapter, 1);
        assert_eq!(m.t3, a2.len());
        // Also the reverse direction must work.
        let mut seq = b"ACGT".repeat(25);
        seq.extend_from_slice(&a1);
        let m = match_adapters(&seq, &set).expect("should match a1");
        assert_eq!(m.adapter, 0);
        assert_eq!(m.t3, a1.len());
    }















    /// The fingerprint is a substring at an insert-dependent offset, so
    /// identification must work from any part of the adapter, not just its
    /// start. This is the case the old code got wrong.
    #[test]
    fn identifies_truseq_from_any_internal_fingerprint() {
        let ad = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        for off in 0..=(ad.len() - 12) {
            let mut reads: Vec<Vec<u8>> = Vec::new();
            for _ in 0..1000 {
                let mut r = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
                r.extend_from_slice(&ad[off..off + 12]);
                reads.push(r);
            }
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let det = identify_from_reads(&refs);
            assert!(
                det.adapters.iter().any(|a| a.seq == ad),
                "offset {off} failed to identify TruSeq",
            );
            // and the installed sequence is the FULL adapter, not the k-mer
            assert!(det.adapters.iter().all(|a| a.seq.len() > 12));
        }
    }

    #[test]
    fn identifies_nothing_on_random_reads() {
        let mut state = 0x1234_5678_9ABC_DEF0u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let reads: Vec<Vec<u8>> = (0..2000)
            .map(|_| (0..60).map(|_| b"ACGT"[(next() % 4) as usize]).collect())
            .collect();
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let det = identify_from_reads(&refs);
        assert!(det.adapters.is_empty(), "identified {:?} on random reads", det.adapters);
        assert_eq!(det.sampled, 2000);
    }

    #[test]
    fn detection_reports_its_evidence() {
        let ad = b"CTGTCTCTTATACACATCT"; // Nextera
        let reads: Vec<Vec<u8>> = (0..500)
            .map(|_| {
                let mut r = b"ACGTACGTACGTACGTACGT".to_vec();
                r.extend_from_slice(&ad[ad.len() - 12..]);
                r
            })
            .collect();
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let det = identify_from_reads(&refs);
        assert_eq!(det.adapters.len(), 1);
        assert_eq!(det.adapters[0].name, "Nextera");
        assert_eq!(det.evidence.len(), 1);
        assert_eq!(det.evidence[0].1, 500);
        assert_eq!(det.sampled, 500);
    }

    /// The case measured on ERR10123282 (MiSeq 2x251, Nextera): the read runs
    /// THROUGH the adapter and continues into the downstream primer, so its
    /// terminal k-mer is primer sequence and the adapter is only visible further
    /// in. Insert lengths vary across reads, so the terminal k-mer differs per
    /// read and no single one reaches the frequency threshold.
    ///
    /// Counted at the 3' terminus only, this detects nothing — on a library where
    /// 26.9% of reads carried the full 19-base adapter and 3.05% ended inside it,
    /// spread over ~19 distinct terminal k-mers.
    #[test]
    fn detection_finds_an_adapter_followed_by_primer_sequence() {
        let ad = b"CTGTCTCTTATACACATCT"; // Nextera, 19 bases
        let primer = b"CCGAGCCCA"; // as observed after the adapter in the real reads
        let reads: Vec<Vec<u8>> = (0..500)
            .map(|i| {
                let insert = 200 + (i % 20); // 20 offsets, ~5% of reads each
                let mut r: Vec<u8> =
                    (0..insert).map(|j| b"ACGT"[(i + j) % 4]).collect();
                r.extend_from_slice(ad);
                r.extend_from_slice(primer);
                r
            })
            .collect();
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let det = identify_from_reads(&refs);
        assert!(
            det.adapters.iter().any(|a| a.name == "Nextera"),
            "adapter followed by primer sequence was not detected: {:?}",
            det.adapters
        );
        // and the evidence must be strong, not a marginal hit
        assert!(
            det.evidence.iter().any(|(_, c)| *c >= 100),
            "expected a strong fingerprint, got {:?}",
            det.evidence
        );
    }

    #[test]
    fn catalogue_is_deduplicated() {
        let c = catalogue();
        for (i, a) in c.iter().enumerate() {
            assert!(!c[i + 1..].iter().any(|b| b.seq == a.seq), "duplicate {}", a.name);
        }
        assert!(c.len() >= 4);
    }

    /// The seed width is sized from the pigeonhole bound instead of being
    /// fixed at 6. Which width each error rate selects, and the overlap at
    /// which that width becomes sound, is the whole contract: a wrong
    /// `min_sound_k` silently changes what gets trimmed.
    ///
    /// Derivation used (sound, and the one the code implements): an alignment
    /// of `k` bases with `mm = ceil(k·eps)` mismatches has at most `mm + 1`
    /// exact runs totalling `k − mm` bases, so the longest run is at least
    /// `floor(k / (mm + 1))` and a seed of width `s` is found iff
    /// `k >= s · (mm + 1)`. At the shipped `--adapter-error-rate 0.2` no
    /// `k` satisfies that for `s = 6` (nor `s = 5`), which is why the old
    /// fixed 6 made every adapter take the reference full scan.
    #[test]
    fn seed_width_and_sound_bound_follow_the_error_rate() {
        let presets: Vec<Vec<u8>> = vec![
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(), // 34
            b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(), // 34
            b"CTGTCTCTTATACACATCT".to_vec(),               // 19
        ];
        let set_at = |eps: f64| {
            AdapterSet::new(
                presets
                    .iter()
                    .enumerate()
                    .map(|(i, s)| Adapter {
                        name: format!("a{i}"),
                        seq: s.clone(),
                    })
                    .collect(),
                eps,
                14,
            )
        };
        // (eps, [width per adapter], [min_sound_k per adapter])
        //
        // `min_sound_k` is the first k after the LAST k where the pigeonhole
        // predicate fails, so it is 16 at the shipped rate rather than the old
        // first-match 20 — and Nextera takes a 4-mer too, because the floor
        // tolerance makes the bound hold from 16 for a 19-base adapter.
        let expect: [(f64, [usize; 3], [usize; 3]); 6] = [
            (0.0, [6, 6, 6], [6, 6, 6]),
            (0.05, [6, 6, 6], [6, 6, 6]),
            (0.1, [6, 6, 6], [12, 12, 12]),
            // 6 needs k >= 6*(mm+1): still satisfiable at eps 0.15 for the
            // 19-base Nextera adapter (ms = 18) and, from k = 30, for the
            // 33-base presentation adapters (ms = 30).
            (0.15, [6, 6, 6], [30, 30, 18]),
            // The shipped rate: no width above 4 is satisfiable at all.
            (0.2, [4, 4, 4], [16, 16, 16]),
            (0.3, [3, 3, 3], [21, 21, 18]),
        ];
        for (eps, widths, ms) in expect {
            let set = set_at(eps);
            assert_eq!(set.seeds.width, widths, "eps {eps}: seed widths");
            assert_eq!(set.seeds.min_sound_k, ms, "eps {eps}: sound k");
        }
    }

    /// The prefilter must actually be LOADED at the shipped error rate, not
    /// merely sound on paper: a clean 151-base read used to take 414 exact
    /// compares per read (3 presets × 138 positions). Sizing the seed from the
    /// bound must leave only the truncated 3′ tail there.
    #[test]
    fn prefilter_is_live_at_the_shipped_error_rate() {
        let presets: Vec<Vec<u8>> = vec![
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
            b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(),
            b"CTGTCTCTTATACACATCT".to_vec(),
        ];
        let set = AdapterSet::new(
            presets
                .iter()
                .enumerate()
                .map(|(i, s)| Adapter {
                    name: format!("a{i}"),
                    seq: s.clone(),
                })
                .collect(),
            0.2,
            14,
        );
        // A 151-base random read with no adapter in it: the worst case for the
        // scan, since every seed hit and every tail position must be checked.
        // (A degenerate `ACGT` repeat would over-count: its 4-mers recur every
        // four bases and hit the table far more often than real sequence.)
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let read: Vec<u8> = (0..151)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect();
        let calls = count_mismatch_calls(|| {
            let _ = match_adapters(&read, &set);
        });
        // Measured: 81 with per-adapter widths (4-mers for the 34-base
        // presets, 3-mers for 19-base Nextera), 174 with one global 4-mer
        // width, 414 with the old fixed 6-mer. It cannot go to a handful: at
        // eps 0.2 a 4-mer is the widest SOUND seed (a 5-mer needs
        // k >= 5*(mm+1), which has no solution above eps = 0.2), and a 4-mer
        // occurs in a random read often enough that ~80 candidates are
        // unavoidable. What matters is that the bound is now live instead of
        // dead — and that every candidate is a real comparison, not a scan.
        assert!(
            calls < 120,
            "prefilter is not doing its job at eps 0.2: {calls} exact compares \
             (was 414 with the fixed 6-mer seed)"
        );
    }

    #[test]
    fn load_fasta_parses_correctly() {
        use std::io::Write;
        let mut p = std::env::temp_dir();
        p.push(format!("polyx-adapter-test-{}.fa", std::process::id()));
        {
            let mut f = std::fs::File::create(&p).unwrap();
            write!(f, ">TruSeq/R1\nAGATCGGAAGAGCACACGT\n>TruSeq/R2\nAGATCGGAAGAGCGTCGTG\n").unwrap();
        }
        let adapters = load_fasta(&p).unwrap();
        std::fs::remove_file(&p).ok();
        assert_eq!(adapters.len(), 2);
        assert_eq!(adapters[0].name, "TruSeq/R1");
        assert_eq!(adapters[1].seq, b"AGATCGGAAGAGCGTCGTG");
    }
}
