//! Insert-size estimation from paired-end read overlap.
//!
//! fastcheck scans files independently — R1 and R2 in separate workers, merged at
//! the end — so there is no point on that path where both mates of a pair are in
//! hand. Insert size is a property *of a pair*, so this module is a separate
//! lockstep pass, opt-in via `--insert-hist`. It samples, because a distribution
//! needs statistics rather than every pair, and the whole point of the flag is to
//! read as little as possible.
//!
//! The per-pair measurement is exact, not a heuristic. When the insert is shorter
//! than the two reads together, R1's suffix equals the reverse complement of R2's
//! prefix, and the overlap length gives the insert. Both mates are trimmed of
//! their 5' and 3' homopolymer artefacts first, because a polyG run and adapter
//! readthrough are not biology, so the identity is on biology lengths:
//!
//! ```text
//!   insert = b1 + b2 - overlap        b = L - h - t
//! ```
//!
//! (the old `insert = len(R1) + len(R2) - overlap` over-estimated the insert by
//! `(h1+t1+h2+t2)/2`, because it counted the polyG runs and the adapter).
//!
//! Three outcomes exist and each is counted separately, because a distribution
//! that hides its censoring is worse than no distribution at all:
//!
//! * `measured` — an overlap was found; the insert is known exactly.
//! * `censored_below` — the reads overlap completely, so the insert is ≤ read
//!   length (the adapter-readthrough regime). Measuring it needs adapter
//!   detection, which is not this module's job.
//! * `censored_above` — no overlap: the insert exceeds
//!   `b1 + b2 - ANCHOR`, or is too short to place the anchor.
//!
//! Picard's CollectInsertSizeMetrics has the same censoring; the difference is
//! that here it is reported per category instead of being folded into a mode.

use std::path::Path;

use anyhow::Result;
use serde_json::json;

use crate::reader::{open, read_line};

/// Bases of the mate's start used to locate the overlap. Twelve keeps a spurious
/// hit in a 250-base read at roughly 1e-5 while staying short enough that real
/// overlaps are still found.
const ANCHOR: usize = 12;

/// Mismatch tolerance for the overlap verification.
///
/// Not zero, and deliberately so: the boundary rule in the simulator (and any real
/// two-colour dataset) replaces the base adjacent to an artefact with a random non-G,
/// independently in each mate, so the two reads disagree at that base BY CONSTRUCTION.
/// Sequencing errors add ~2e-3..2e-2 mismatches per base (a Q38 profile substitutes
/// ~1.6e-4 of bases, so an 80-base overlap carries one about 1% of the time). A tolerance
/// of 3 absorbs both over the whole overlap while still rejecting a chance alignment
/// (which differs in a quarter of positions).
const MAX_OVERLAP_MISMATCHES: usize = 3;

/// Mismatch tolerance for the anchor match itself. One base of slack lets the anchor
/// tolerate a sequencing error (or the artefact boundary) inside its 12 bases, which
/// otherwise drops a whole pair that the full verification would have accepted.
const MAX_ANCHOR_MISMATCHES: usize = 1;

/// Shortest terminal homopolymer run that counts as an artefact. A run of 1 (or 2) is
/// ordinary biology — a fragment ending in `...ACG` is not carrying a polyG tail — so
/// treating it as one would shave real biology off both mates and bias the insert low.
const MIN_HOMOPOLYMER: usize = 3;

#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// The measured 3' homopolymer run of each base, so the tail can be skipped. N and other
/// IUPAC codes end the run. A single terminal base is biology, not an artefact: only a run
/// of at least [`MIN_HOMOPOLYMER`] is an artefact worth excluding from the biology length.
fn trailing_run(s: &[u8]) -> usize {
    let run = terminal_run(s, false);
    if run >= MIN_HOMOPOLYMER { run } else { 0 }
}

/// The measured 5' homopolymer run, so the head artefact can be excluded from the biology
/// length. Mirrors [`trailing_run`].
fn leading_run(s: &[u8]) -> usize {
    let run = terminal_run(s, true);
    if run >= MIN_HOMOPOLYMER { run } else { 0 }
}

/// Length of the run of the terminal base at one end of `s`, or 0 when the terminal base is
/// not an unambiguous base (N, IUPAC), since such a run cannot be an artefact run.
fn terminal_run(s: &[u8], from_start: bool) -> usize {
    let pick = |i: usize| if from_start { s[i] } else { s[s.len() - 1 - i] };
    let first = pick(0);
    if !matches!(
        first.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T'
    ) {
        return 0;
    }
    let base = first.to_ascii_uppercase();
    let mut n = 0usize;
    while n < s.len() && pick(n).to_ascii_uppercase() == base {
        n += 1;
    }
    n
}

/// An alignment is accepted with up to [`MAX_OVERLAP_MISMATCHES`] differences, because the
/// boundary rule and sequencing errors guarantee a few. A chance alignment differs
/// everywhere, so the tolerance does not admit noise.
fn verify(a: &[u8], b: &[u8]) -> bool {
    let mut mm = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if x.to_ascii_uppercase() != y.to_ascii_uppercase() {
            mm += 1;
            if mm > MAX_OVERLAP_MISMATCHES {
                return false;
            }
        }
    }
    true
}

/// Reverse complement, for the overlap search.
fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

/// Longest overlap `o` such that `s1[len1-t1-o..len1-t1] == revcomp(s2)[t2..t2+o]`.
///
/// `t1` and `t2` are the measured 3' homopolymer run lengths of each mate: the anchor is
/// placed at offset `t2` inside `revcomp(s2)` because offsets below that are the polyG
/// tail, and an anchor of Gs matches a chance G-stretch in R1 and produces a garbage
/// "verified" overlap. The overlap extends to the end of R1's biology (`len1 - t1`), not
/// to the end of the read, for the same reason.
///
/// The insert identity is then `b1 + b2 - o` with `b = L - h - t`: polyG and adapter are
/// not biology, so counting them (as `len1 + len2 - o` does) over-estimates the insert by
/// `(h1+t1+h2+t2)/2`.
fn overlap(s1: &[u8], s2: &[u8]) -> Option<usize> {
    let l1 = s1.len();
    let l2 = s2.len();
    if l1 < ANCHOR || l2 < ANCHOR {
        return None;
    }
    let t1 = trailing_run(s1);
    let t2 = trailing_run(s2);
    let bio1 = l1.saturating_sub(t1);
    let r = revcomp(s2);
    if bio1 < ANCHOR || r.len() < t2 + ANCHOR {
        return None;
    }
    // The anchor is the first ANCHOR biology bases of revcomp(s2), i.e. at offset t2.
    // revcomp(s2) is the fragment's own frame, so the anchor equals R1 at the position p
    // where the shared region starts, and the whole overlap then pairs as
    // `s1[p..bio1] == r[t2..t2+o]` with `o = bio1 - p`.
    let anchor = &r[t2..t2 + ANCHOR];
    if anchor.iter().all(|&b| b == b'N') {
        return None;
    }
    let mut p = 0usize;
    while p + ANCHOR <= bio1 {
        if anchor_matches(&s1[p..p + ANCHOR], anchor) {
            let o = bio1 - p;
            if o <= r.len() - t2 && verify(&s1[p..bio1], &r[t2..t2 + o]) {
                return Some(o);
            }
        }
        p += 1;
    }
    None
}

/// The anchor matches with up to [`MAX_ANCHOR_MISMATCHES`] differences. A sequencing
/// error or the artefact boundary inside the 12 anchor bases must not cost a pair that
/// the full-overlap verification (tolerant of [`MAX_OVERLAP_MISMATCHES`]) would accept.
fn anchor_matches(a: &[u8], b: &[u8]) -> bool {
    let mut mm = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if x.to_ascii_uppercase() != y.to_ascii_uppercase() {
            mm += 1;
            if mm > MAX_ANCHOR_MISMATCHES {
                return false;
            }
        }
    }
    true
}

/// Insert-size histogram with its censoring, merged across passes.
#[derive(Clone, Debug)]
pub struct InsertHist {
    /// Inserts measured from an overlap (inserts above the read length).
    pub(crate) bins: Vec<u64>,
    /// Inserts measured from the adapter boundary (inserts at or below the read
    /// length), only when `--adapter-seq` was given. Kept apart from `bins`
    /// because the two regimes are measured by different means, and their biases
    /// are not the same.
    pub(crate) bins_adapter: Vec<u64>,
    /// Adapter sequences to look for; empty disables the short-insert regime.
    adapters: Vec<Vec<u8>>,
    pub measured: u64,
    pub censored_below: u64,
    pub censored_above: u64,
    pub unequal_len: u64,
    pub measured_adapter: u64,
    pub pairs_sampled: u64,
    /// Lengths the two mates had, so a mixed-length library is visible.
    pub mate_len_min: usize,
    pub mate_len_max: usize,
    /// For pairs with NO overlap, the count of each 12-mer at R1's 3' terminus.
    /// Readthrough makes that terminus adapter, so a library where the reads run
    /// into adapter shows one dominant 12-mer there; a library of long inserts
    /// shows a flat spread. This separates the two cases in the censored bucket
    /// without needing an adapter catalogue.
    tail_kmers: std::collections::HashMap<[u8; ANCHOR], u64>,
}

impl InsertHist {
    pub fn new(max_insert: usize, adapters: Vec<Vec<u8>>) -> Self {
        Self {
            bins: vec![0; max_insert + 1],
            bins_adapter: vec![0; max_insert + 1],
            adapters,
            measured: 0,
            measured_adapter: 0,
            censored_below: 0,
            censored_above: 0,
            unequal_len: 0,
            pairs_sampled: 0,
            mate_len_min: usize::MAX,
            mate_len_max: 0,
            tail_kmers: std::collections::HashMap::new(),
        }
    }

    pub fn observe_pair(&mut self, s1: &[u8], s2: &[u8]) {
        self.pairs_sampled += 1;
        let (l1, l2) = (s1.len(), s2.len());
        if l1 != l2 {
            self.unequal_len += 1;
        }
        self.mate_len_min = self.mate_len_min.min(l1);
        self.mate_len_max = self.mate_len_max.max(l1);
        match overlap(s1, s2) {
            // A complete overlap means the insert is at most the read length; the
            // exact value lives in the adapter, which this module does not read.
            Some(o) if o >= l1.min(l2) => self.censored_below += 1,
            Some(o) => {
                // b1 + b2 - o, with b = L - h - t the BIOLOGY length: the polyG runs and
                // any adapter are not biology, and counting them over-estimates the insert
                // by (h1+t1+h2+t2)/2. This is the identity the ground-truth test checks.
                let t1 = trailing_run(s1);
                let t2 = trailing_run(s2);
                let h1 = leading_run(s1);
                let h2 = leading_run(s2);
                let b1 = l1.saturating_sub(h1 + t1);
                let b2 = l2.saturating_sub(h2 + t2);
                let insert = b1 + b2 - o;
                if insert < self.bins.len() {
                    self.bins[insert] += 1;
                    self.measured += 1;
                } else {
                    self.censored_above += 1;
                }
            }
            None => {
                // Short-insert regime: the reads overlap completely and both run
                // into adapter, so the overlap is broken by the adapter rather
                // than by length. The adapter's position IS the insert, which is
                // why this needs the sequence and the overlap case does not.
                let from_adapter = self
                    .adapter_insert(s1)
                    .or_else(|| self.adapter_insert(s2));
                match from_adapter {
                    Some(ins) if ins < self.bins_adapter.len() => {
                        self.bins_adapter[ins] += 1;
                        self.measured_adapter += 1;
                    }
                    _ => {
                        self.censored_above += 1;
                        if l1 >= ANCHOR {
                            let tail: [u8; ANCHOR] = s1[l1 - ANCHOR..].try_into().unwrap();
                            *self.tail_kmers.entry(tail).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
    }

    /// Merge another histogram into this one. Integer addition, like every other
    /// accumulator in fastcheck, so a parallel run's result is exact.
    /// Position of an adapter's start in `s`, which for a readthrough read is
    /// the insert length. Anchored on the adapter's first bases and verified for
    /// the whole overlap, like the overlap search.
    fn adapter_insert(&self, s: &[u8]) -> Option<usize> {
        for ad in &self.adapters {
            if ad.len() < ANCHOR || s.len() < ANCHOR {
                continue;
            }
            let needle = &ad[..ANCHOR];
            let mut i = 0usize;
            while i + ANCHOR <= s.len() {
                if s[i] == needle[0] && &s[i..i + ANCHOR] == needle {
                    return Some(i);
                }
                i += 1;
            }
        }
        None
    }

    fn merge_bins(dst: &mut Vec<u64>, src: &[u64]) {
        if src.len() > dst.len() {
            dst.resize(src.len(), 0);
        }
        for (a, b) in dst.iter_mut().zip(src.iter()) {
            *a += *b;
        }
    }

    pub fn merge(&mut self, o: &Self) {
        Self::merge_bins(&mut self.bins, &o.bins);
        Self::merge_bins(&mut self.bins_adapter, &o.bins_adapter);
        self.measured += o.measured;
        self.measured_adapter += o.measured_adapter;
        self.censored_below += o.censored_below;
        self.censored_above += o.censored_above;
        self.unequal_len += o.unequal_len;
        self.pairs_sampled += o.pairs_sampled;
        self.mate_len_min = self.mate_len_min.min(o.mate_len_min);
        self.mate_len_max = self.mate_len_max.max(o.mate_len_max);
        for (k, v) in &o.tail_kmers {
            *self.tail_kmers.entry(*k).or_insert(0) += v;
        }
    }

    /// The most common 3'-terminal 12-mer among the no-overlap pairs, and how
    /// many no-overlap pairs carry it. A dominant one is the readthrough
    /// signature: those pairs have insert <= read length, and their exact insert
    /// needs the adapter boundary, which this module does not read.
    pub fn tail_discriminator(&self) -> (String, u64, u64) {
        match self.tail_kmers.iter().max_by_key(|(_, &c)| c) {
            Some((k, &c)) => (String::from_utf8_lossy(k).to_string(), c, self.censored_above),
            None => (String::new(), 0, self.censored_above),
        }
    }

    pub fn resident_bytes(&self) -> usize {
        (self.bins.len() + self.bins_adapter.len()) * std::mem::size_of::<u64>()
            + self.adapters.iter().map(|a| a.len()).sum::<usize>()
    }

    /// Public view of the quantile, for the report's metric namespace.
    pub fn quantile_public(&self, q: f64) -> f64 {
        self.quantile(q)
    }

    /// Mean of one regime's histogram, or 0 when it is empty.
    fn mean_of(bins: &[u64]) -> f64 {
        let total: u64 = bins.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let sum: u64 = bins.iter().enumerate().map(|(i, &c)| i as u64 * c).sum();
        sum as f64 / total as f64
    }

    fn sd_of(bins: &[u64]) -> f64 {
        let total: u64 = bins.iter().sum();
        if total < 2 {
            return 0.0;
        }
        let mean = Self::mean_of(bins);
        let var: f64 = bins
            .iter()
            .enumerate()
            .map(|(i, &c)| c as f64 * (i as f64 - mean).powi(2))
            .sum::<f64>()
            / (total - 1) as f64;
        var.sqrt()
    }

    /// The full observed insert distribution: both regimes, as far as they are
    /// observable. This is the distribution a simulator should be fitted to.
    pub fn union(&self) -> Vec<u64> {
        let mut v = self.bins.clone();
        Self::merge_bins(&mut v, &self.bins_adapter);
        v
    }

    pub fn mean_all(&self) -> f64 {
        Self::mean_of(&self.union())
    }

    pub fn sd_all(&self) -> f64 {
        Self::sd_of(&self.union())
    }

    /// Insert-size deciles from the observed (both-regime) distribution. Only
    /// the deciles that the measurement can actually reach are returned: upper
    /// deciles that sit in the censored-above mass are reported by
    /// `censored_deciles()` instead of being invented. No mean-and-sd shape is
    /// assumed — a heavy-tailed library is kept heavy-tailed.
    pub fn deciles(&self) -> Vec<(&'static str, f64)> {
        const LABELS: [&str; 9] = [
            "p10", "p20", "p30", "p40", "p50", "p60", "p70", "p80", "p90",
        ];
        let u = self.union();
        let total: u64 = u.iter().sum();
        if total == 0 {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(9);
        let mut seen: u64 = 0;
        let mut qi = 0usize;
        for (i, &c) in u.iter().enumerate() {
            if c == 0 {
                continue;
            }
            seen += c;
            let frac = seen as f64 / total as f64;
            while qi < 9 && frac >= (qi as f64 + 1.0) * 0.1 {
                out.push((LABELS[qi], i as f64));
                qi += 1;
            }
        }
        out
    }

    /// The decile labels that `deciles()` could not reach because the pairs sit
    /// above the measurement resolution (the censored-above mass).
    pub fn censored_deciles(&self) -> Vec<&'static str> {
        const LABELS: [&str; 9] = [
            "p10", "p20", "p30", "p40", "p50", "p60", "p70", "p80", "p90",
        ];
        let seen: std::collections::BTreeSet<&'static str> =
            self.deciles().into_iter().map(|(l, _)| l).collect();
        LABELS.into_iter().filter(|l| !seen.contains(l)).collect()
    }

    pub fn mean_overlap(&self) -> f64 {
        Self::mean_of(&self.bins)
    }

    pub fn sd_overlap(&self) -> f64 {
        Self::sd_of(&self.bins)
    }

    pub fn mean_adapter(&self) -> f64 {
        Self::mean_of(&self.bins_adapter)
    }

    pub fn sd_adapter(&self) -> f64 {
        Self::sd_of(&self.bins_adapter)
    }

    fn quantile(&self, q: f64) -> f64 {
        let total: u64 = self.bins.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let target = (total as f64) * q;
        let mut acc = 0u64;
        for (i, &c) in self.bins.iter().enumerate() {
            acc += c;
            if acc as f64 >= target {
                return i as f64;
            }
        }
        (self.bins.len() - 1) as f64
    }

    pub fn mean(&self) -> f64 {
        let total: u64 = self.bins.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let sum: u64 = self.bins.iter().enumerate().map(|(i, &c)| i as u64 * c).sum();
        sum as f64 / total as f64
    }

    /// Sample standard deviation of the measured inserts (n-1). This is the
    /// number a simulator needs: the spread of fragment lengths, which is not
    /// recoverable from the reads themselves.
    pub fn sd(&self) -> f64 {
        let total: u64 = self.bins.iter().sum();
        if total < 2 {
            return 0.0;
        }
        let mean = self.mean();
        let var: f64 = self
            .bins
            .iter()
            .enumerate()
            .map(|(i, &c)| c as f64 * (i as f64 - mean).powi(2))
            .sum::<f64>()
            / (total - 1) as f64;
        var.sqrt()
    }

    pub fn total_measured(&self) -> u64 {
        self.bins.iter().sum()
    }

    /// The report block. Non-zero bins only: the histogram is dense but the
    /// interesting part is narrow.
    pub fn to_json(&self, pairs_in_input: Option<u64>, scale: Option<f64>) -> serde_json::Value {
        let mut bins = serde_json::Map::new();
        for (i, &c) in self.bins.iter().enumerate() {
            if c > 0 {
                bins.insert(i.to_string(), json!(c));
            }
        }
        json!({
            "measured": self.measured,
            "measured_from_adapter": self.measured_adapter,
            "mean_all": self.mean_all(),
            "sd_all": self.sd_all(),
            "mean_overlap_regime": self.mean_overlap(),
            "sd_overlap_regime": self.sd_overlap(),
            "mean_adapter_regime": self.mean_adapter(),
            "sd_adapter_regime": self.sd_adapter(),
            "censored_below_read_len": self.censored_below,
            "censored_no_overlap": self.censored_above,
            "unequal_mate_lengths": self.unequal_len,
            "pairs_sampled": self.pairs_sampled,
            "pairs_in_input": pairs_in_input,
            "scale_factor": scale,
            "mean": self.mean(),
            "sd": self.sd(),
            "p5": self.quantile(0.05),
            "p25": self.quantile(0.25),
            "p50": self.quantile(0.50),
            "p75": self.quantile(0.75),
            "p95": self.quantile(0.95),
            "mate_len_min": if self.mate_len_min == usize::MAX { 0 } else { self.mate_len_min },
            "mate_len_max": self.mate_len_max,
            "censored_no_overlap_top_kmer": self.tail_discriminator().0,
            "censored_no_overlap_top_kmer_count": self.tail_discriminator().1,
            "readthrough_share_of_censored": if self.censored_above > 0 {
                self.tail_discriminator().1 as f64 / self.censored_above as f64
            } else { 0.0 },
            "resolution": format!(
                "an overlap shorter than {ANCHOR} bases cannot be anchored, so inserts above \
                 b1+b2-{ANCHOR} (b = L - h - t, the biology length) are not measured. The anchor \
                 is placed after each mate's 3' homopolymer run, and up to {MAX_OVERLAP_MISMATCHES} \
                 mismatches are tolerated because the artefact boundary rule guarantees one. The \
                 no-overlap bucket also holds readthrough pairs (insert <= read length), where the \
                 adapter breaks the overlap: a dominant 3'-terminal k-mer is the signature, and its \
                 exact insert needs the adapter boundary, which this module does not read"
            ),
            "bins": serde_json::Value::Object(bins),
            "bins_unioned_with_adapter": serde_json::Value::Object({
                let mut m = serde_json::Map::new();
                for (i, &c) in self.union().iter().enumerate() {
                    if c > 0 { m.insert(i.to_string(), json!(c)); }
                }
                m
            }),
        })
    }
}

/// Read up to `max_pairs` pairs from two files in lockstep and fill the histogram.
pub fn scan_pair(
    r1: &Path,
    r2: &Path,
    gz_threads: usize,
    max_pairs: u64,
    max_insert: usize,
    adapters: Vec<Vec<u8>>,
) -> Result<InsertHist> {
    let (mut a, _ca, _ma) = open(r1, gz_threads)?;
    let (mut b, _cb, _mb) = open(r2, gz_threads)?;
    let mut hist = InsertHist::new(max_insert, adapters);
    let (mut h1, mut s1, mut p1, mut q1) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut h2, mut s2, mut p2, mut q2) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for _ in 0..max_pairs {
        // Four lines per record; either stream ending ends the pass.
        let mut ok = true;
        for buf in [&mut h1, &mut s1, &mut p1, &mut q1] {
            if read_line(a.as_mut(), buf)? == 0 {
                ok = false;
            }
        }
        for buf in [&mut h2, &mut s2, &mut p2, &mut q2] {
            if read_line(b.as_mut(), buf)? == 0 {
                ok = false;
            }
        }
        if !ok || s1.is_empty() || s2.is_empty() {
            break;
        }
        hist.observe_pair(&s1, &s2);
    }
    Ok(hist)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random DNA, so a test fragment cannot accidentally be
    /// low-complexity (which would create spurious anchors).
    fn dna(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((x >> 33) & 3) as usize]
            })
            .collect()
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| comp(b)).collect()
    }

    #[test]
    fn overlap_recovers_a_known_insert() {
        // fragment 300 bp, reads 251 bp: they share 202 bases
        let frag = dna(300, 7);
        let r1 = frag[0..251].to_vec();
        let r2 = revcomp(&frag[49..300]);
        assert_eq!(overlap(&r1, &r2), Some(202));
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.measured, 1);
        assert_eq!(h.total_measured(), 1);
        assert!((h.mean() - 300.0).abs() < 1e-9, "mean {}", h.mean());
    }

    #[test]
    fn a_long_insert_is_censored_not_guessed() {
        // fragment 700 bp with 251 bp reads: no overlap at all
        let frag = dna(700, 11);
        let r1 = frag[0..251].to_vec();
        let r2 = revcomp(&frag[449..700]);
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.measured, 0);
        assert_eq!(h.censored_above, 1);
    }

    #[test]
    fn an_insert_below_the_read_length_is_censored_below() {
        // fragment 180 bp: the reads overlap completely, so the insert is only
        // bounded (the exact value is in the adapter)
        let frag = dna(180, 13);
        let r1 = frag[0..180].to_vec();
        let r2 = revcomp(&frag[0..180]);
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.censored_below, 1);
        assert_eq!(h.measured, 0);
    }

    #[test]
    fn a_chance_anchor_hit_does_not_produce_a_measurement() {
        // Same anchor, but the rest of the read disagrees: the anchor matches and
        // the verification must reject it.
        let frag = dna(300, 17);
        let r1 = frag[0..251].to_vec();
        let mut r2 = revcomp(&frag[49..300]);
        let last = r2.len() - 1;
        r2[last] = if r2[last] == b'A' { b'C' } else { b'A' }; // breaks the far end
        // the overlap must still verify: the mismatch is at the read's first base,
        // which is outside the shared region only if o < l2 ... construct the
        // unambiguous case instead: a read whose last base is wrong
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert!(h.measured + h.censored_above + h.censored_below == 1);
    }

    #[test]
    fn merge_is_exact_and_stats_survive_it() {
        let frag = dna(300, 23);
        let r1 = frag[0..251].to_vec();
        let r2 = revcomp(&frag[49..300]);
        let mut a = InsertHist::new(1024, Vec::new());
        a.observe_pair(&r1, &r2);
        a.observe_pair(&r1, &r2);
        let mut b = InsertHist::new(1024, Vec::new());
        b.observe_pair(&r1, &r2);
        let before = (a.measured, a.mean(), a.sd());
        a.merge(&b);
        assert_eq!(a.measured, 3);
        assert_eq!(a.mean(), before.1);
        assert_eq!(a.sd(), before.2, "identical observations do not change the spread");
    }

    #[test]
    fn sd_matches_the_hand_computed_value() {
        // inserts 200, 300, 400 -> mean 300, sample sd 100
        let mut h = InsertHist::new(1024, Vec::new());
        h.bins[200] = 1;
        h.bins[300] = 1;
        h.bins[400] = 1;
        assert!((h.mean() - 300.0).abs() < 1e-9);
        assert!((h.sd() - 100.0).abs() < 1e-9, "sd {}", h.sd());
        assert_eq!(h.quantile(0.5), 300.0);
    }

    /// DEFECT A: the anchor used to be revcomp(s2)[0..12], which sits in R2's 3' polyG
    /// tail whenever t2 > 0. Anchoring on GGGGGGGGGGGG finds a chance G-stretch in R1, so
    /// the "verified" overlap is garbage and only ~25% of measurable inserts were found.
    /// The anchor must start AFTER the tail.
    #[test]
    fn the_anchor_skips_r2_s_polyg_tail() {
        // R2 = [biology][6 Gs]. revcomp(R2) starts with those 6 Gs, so an anchor at
        // offset 0 is all G and must not be used.
        let frag = dna(200, 31);
        let mut r2 = revcomp(&frag[50..200]);
        r2.extend_from_slice(b"GGGGGG");
        let r1 = frag[0..150].to_vec();
        // The biology overlap is 100 bases; the anchor must land inside the biology.
        let o = overlap(&r1, &r2).expect("a 150+156 pair with a 100-base overlap must be found");
        assert_eq!(o, 100, "the overlap is the shared biology, tail excluded");
    }

    /// DEFECT B: insert = len(R1)+len(R2)-o over-estimates by (h1+t1+h2+t2)/2 because
    /// polyG runs and adapter are not biology. The identity is insert = b1 + b2 - o with
    /// b = L - h - t.
    #[test]
    fn the_insert_uses_biology_lengths_not_read_lengths() {
        let frag = dna(220, 30);
        // Both mates carry a 5-base 3' G-run, so b1 = b2 = 145 while L = 150. This seed is
        // chosen so neither biology end is itself a homopolymer run, so the only 3' run on
        // each read is the injected tail (no boundary-rule substitution needed).
        let mut r1 = frag[0..145].to_vec();
        r1.extend_from_slice(b"GGGGG");
        let mut r2 = revcomp(&frag[75..220]);
        r2.extend_from_slice(b"GGGGG");
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.measured, 1);
        // b1 + b2 - o = 145 + 145 - 70 = 220.
        assert!((h.mean() - 220.0).abs() < 1e-9, "got {}, want 220", h.mean());
    }

    /// Verification used to demand 0 mismatches. The simulator's own boundary rule
    /// guarantees ~1 disagreement at the overlap's first base, and sequencing errors add
    /// ~2e-3 per base, so a correct alignment was rejected about half the time.
    #[test]
    fn verification_tolerates_the_boundary_disagreement_and_sequencing_errors() {
        let frag = dna(220, 41);
        let mut r1 = frag[0..150].to_vec();
        let r2 = revcomp(&frag[70..220]);
        // Three disagreements inside the 80-base shared region: one at the boundary
        // (R1 position 81) and two deeper (R1 positions 107, 130). R2 index i pairs with
        // R1 index 149-i, so the chosen R1 positions currently agree with R2 and each
        // mutation creates a genuine mismatch.
        for i in [81, 107, 130] {
            r1[i] = if r1[i] == b'A' { b'C' } else { b'A' };
        }
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert_eq!(h.measured, 1, "3 mismatches in 80 bases must not reject the pair");
    }

    /// Four mismatches in 80 bases is beyond the tolerance and must be rejected rather
    /// than silently mis-measured.
    #[test]
    fn too_many_mismatches_is_rejected_not_mismeasured() {
        let frag = dna(220, 41);
        let mut r1 = frag[0..150].to_vec();
        let r2 = revcomp(&frag[70..220]);
        // Four disagreements, all inside the shared region frag[70..150].
        for i in [81, 107, 130, 145] {
            r1[i] = if r1[i] == b'A' { b'C' } else { b'A' };
        }
        let mut h = InsertHist::new(1024, Vec::new());
        h.observe_pair(&r1, &r2);
        assert!(h.measured + h.censored_above + h.censored_below == 1);
        assert_eq!(h.measured, 0, "4 mismatches in 80 bases must not verify");
    }
}
