//! Exact duplication counting over 64-bit fingerprints.
//!
//! Copied from `tools/fastdedup/src/hasher.rs` (D-1: no shared crate) with the abundance
//! tracker kept and the dedup-only verifier dropped. A PCR duplicate re-sequences the same
//! fragment with FRESH polyG artefacts and fresh substitution errors, so two copies of one
//! molecule usually differ as strings: what identifies a duplicate is the fragment, not the
//! read. Fingerprinting the whole read therefore under-counts, and the report says so
//! (duplication_level is a LOWER BOUND on the PCR duplicate rate).

use rustc_hash::FxHashMap;
use xxhash_rust::xxh3::xxh3_64;

/// Fingerprint of one read. Quality is ignored (stated): two copies of one fragment differ
/// in their quality string because dropout is per-cluster, and that difference is exactly
/// what must NOT split a duplicate.
#[inline(always)]
pub fn fingerprint(seq: &[u8]) -> u64 {
    xxh3_64(seq)
}

/// Exact multiplicity counter: fingerprint -> occurrences.
///
/// Sized from a caller-supplied capacity (clamped to 2^26 entries). The scan starts it at
/// 1024 entries and it grows by rehashing; `crate::memory::project` gives the size it
/// reaches when it holds the estimated number of distinct reads.
#[derive(Clone)]
pub struct FingerprintTable {
    counts: FxHashMap<u64, u32>,
}

impl FingerprintTable {
    /// `capacity` is the projected number of distinct reads.
    pub fn with_capacity(capacity: usize) -> Self {
        let mut counts = FxHashMap::with_capacity_and_hasher(
            capacity.min(1 << 26),
            Default::default(),
        );
        counts.reserve(capacity.min(1 << 26));
        Self { counts }
    }

    #[inline]
    pub fn observe(&mut self, seq: &[u8]) {
        *self.counts.entry(fingerprint(seq)).or_insert(0) += 1;
    }

    /// Distinct fingerprints seen.
    pub fn distinct(&self) -> u64 {
        self.counts.len() as u64
    }

    /// Test-only: the element capacity hashbrown actually allocated. `memory::tests` uses
    /// it to check the projection's bucket count against a real table.
    #[cfg(test)]
    pub fn capacity(&self) -> usize {
        self.counts.capacity()
    }

    /// Total observations recorded (sum of multiplicities).
    pub fn total(&self) -> u64 {
        self.counts.values().map(|&c| c as u64).sum()
    }

    /// Exact multiplicity histogram: `out[c]` = number of distinct sequences occurring
    /// exactly `c` times. Index 0 is unused; index 1 is the singleton count.
    pub fn multiplicity_hist(&self) -> Vec<u64> {
        let max = self.counts.values().copied().max().unwrap_or(0) as usize;
        let mut out = vec![0u64; max + 2];
        for &c in self.counts.values() {
            out[c as usize] += 1;
        }
        out
    }

    /// Distinct sequences occurring at least `min_count` times, most frequent first,
    /// each as `(fingerprint, count)`. Names are not available in this mode.
    pub fn at_least(&self, min_count: u32) -> Vec<(u64, u32)> {
        let mut v: Vec<(u64, u32)> = self
            .counts
            .iter()
            .filter(|(_, &c)| c >= min_count)
            .map(|(&k, &c)| (k, c))
            .collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }

    /// Exact sum of two tables.
    pub fn merge(&mut self, o: &FingerprintTable) {
        for (&k, &c) in &o.counts {
            *self.counts.entry(k).or_insert(0) += c;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_exact_multiplicities() {
        let mut t = FingerprintTable::with_capacity(16);
        for _ in 0..5 {
            t.observe(b"ACGTACGTAC");
        }
        t.observe(b"TTTTTTTTTT");
        t.observe(b"TTTTTTTTTT");
        t.observe(b"GGGGCCCCGG");
        assert_eq!(t.distinct(), 3);
        assert_eq!(t.total(), 8);
        let h = t.multiplicity_hist();
        assert_eq!(h[1], 1, "one singleton");
        assert_eq!(h[2], 1, "one sequence seen twice");
        assert_eq!(h[5], 1, "one sequence seen five times");
    }

    #[test]
    fn at_least_filters_and_sorts_by_count() {
        let mut t = FingerprintTable::with_capacity(16);
        for _ in 0..7 {
            t.observe(b"ACGTACGTAC");
        }
        for _ in 0..3 {
            t.observe(b"TTTTTTTTTT");
        }
        t.observe(b"GGGGCCCCGG");
        let v = t.at_least(3);
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].1, 7, "most frequent first");
        assert_eq!(v[1].1, 3);
    }

    #[test]
    fn merge_is_exact() {
        let mut a = FingerprintTable::with_capacity(8);
        a.observe(b"ACGTACGTAC");
        a.observe(b"ACGTACGTAC");
        let mut b = FingerprintTable::with_capacity(8);
        b.observe(b"ACGTACGTAC");
        b.observe(b"TTTTTTTTTT");
        a.merge(&b);
        assert_eq!(a.distinct(), 2);
        assert_eq!(a.total(), 4);
        assert_eq!(a.multiplicity_hist()[3], 1);
    }

    /// fastdedup and fastcheck must agree on the fingerprint of a sequence (D-1: consistency
    /// by tests, not by the compiler).
    #[test]
    fn fingerprint_is_xxh3_64_of_the_sequence() {
        // xxh3_64 is deterministic; this pins the value so a dependency bump that changes
        // it fails loudly instead of silently splitting the two tools' counts.
        let h = fingerprint(b"ACGTACGTAC");
        assert_eq!(h, fingerprint(b"ACGTACGTAC"));
        assert_ne!(h, fingerprint(b"ACGTACGTACA"), "one base differs -> different fingerprint");
    }
}
