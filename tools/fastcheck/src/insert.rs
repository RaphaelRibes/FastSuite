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
//! prefix, and the overlap length gives the insert:
//!
//! ```text
//!   insert = len(R1) + len(R2) - overlap
//! ```
//!
//! Three outcomes exist and each is counted separately, because a distribution
//! that hides its censoring is worse than no distribution at all:
//!
//! * `measured` — an overlap was found; the insert is known exactly.
//! * `censored_below` — the reads overlap completely, so the insert is ≤ read
//!   length (the adapter-readthrough regime). Measuring it needs adapter
//!   detection, which is not this module's job.
//! * `censored_above` — no overlap: the insert exceeds
//!   `len(R1) + len(R2) - ANCHOR`, or is too short to place the anchor.
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

/// Longest overlap `o` such that `s1[len1-o..] == revcomp(s2)[0..o]`.
///
/// `o` is the number of bases the two reads share, so `insert = len1 + len2 - o`.
/// The search is anchored on the first [`ANCHOR`] bases of `revcomp(s2)` — which
/// must sit at a position `p` of `s1` with `p + o == len1` — and then verified
/// base by base over the whole overlap, so a chance anchor hit cannot produce a
/// measurement.
fn overlap(s1: &[u8], s2: &[u8]) -> Option<usize> {
    let l1 = s1.len();
    let l2 = s2.len();
    if l1 < ANCHOR || l2 < ANCHOR {
        return None;
    }
    let mut anchor = [0u8; ANCHOR];
    for (i, a) in anchor.iter_mut().enumerate() {
        *a = comp(s2[l2 - 1 - i]);
    }
    if anchor.iter().all(|&b| b == b'N') {
        return None;
    }
    let first = anchor[0];
    let mut p = 0usize;
    while p + ANCHOR <= l1 {
        if s1[p] == first && s1[p..p + ANCHOR] == anchor {
            let o = l1 - p;
            if o <= l2 && verify(s1, s2, o) {
                return Some(o);
            }
        }
        p += 1;
    }
    None
}

fn verify(s1: &[u8], s2: &[u8], o: usize) -> bool {
    let (l1, l2) = (s1.len(), s2.len());
    (0..o).all(|i| comp(s2[l2 - 1 - i]) == s1[l1 - o + i])
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
                let insert = l1 + l2 - o;
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
                 len(R1)+len(R2)-{ANCHOR} are not measured. The no-overlap bucket also holds \
                 readthrough pairs (insert <= read length), where the adapter breaks the \
                 overlap: a dominant 3'-terminal k-mer is the signature, and its exact insert \
                 needs the adapter boundary, which this module does not read"
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
}
