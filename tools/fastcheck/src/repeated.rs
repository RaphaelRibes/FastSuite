//! Bounded over-represented-sequence and k-mer counters.
//!
//! Streaming with singleton eviction: when the number of distinct keys exceeds
//! a cap, every key with count == 1 is dropped. A singleton can never be
//! over-represented, and losing its first occurrence costs at most 1 of a count
//! that must reach a large threshold to be reported, so memory stays bounded on
//! 50 M-read files while the frequent sequences are counted exactly.
//!
//! SCOPE (as of the exactness work): this counter is used for the **k-mer spectrum only**.
//! Full-sequence duplication and over-representation moved to the exact fingerprint table
//! (`crate::fingerprint`) and the packed sequence table (`crate::overrep`), because a
//! fingerprint is 8 B where a sequence is 46 B and the exactness is affordable. The
//! approximation that remains is labelled in the report (`kmer_spectrum.approximate`).

use std::collections::HashMap;
use std::hash::Hash;

/// Default cap on distinct keys before singleton eviction kicks in.
pub const DEFAULT_CAP: usize = 1_000_000;

/// A bounded counter keyed by any `Hash + Eq + Clone` key (a full read
/// sequence as `Vec<u8>`, or a packed k-mer as `u64`).
#[derive(Clone)]
pub struct RepeatedCounter<K> {
    map: HashMap<K, u32>,
    cap: usize,
    /// Total inserts seen (for the fraction).
    pub total: u64,
    /// Map length at which the next singleton eviction fires. Keeps the
    /// eviction amortized O(1) per insert (see `insert`).
    evict_at: usize,
}

impl<K: Hash + Eq + Clone + Ord> RepeatedCounter<K> {
    pub fn new(cap: usize) -> Self {
        Self {
            map: HashMap::with_capacity(cap.min(1 << 20)),
            cap,
            total: 0,
            evict_at: cap + cap / 4,
        }
    }

    /// Record one occurrence of `key`, evicting singletons if over the cap.
    ///
    /// AMORTIZED EVICTION: running `retain` on every insert while over the
    /// cap is O(n) per call — on a saturated k-mer space (100k reads x ~140
    /// 12-mers against 4^12 = 16.7M keys) the map sits at the cap
    /// permanently and every insert pays a full 1M-entry rehash: ~10^13
    /// operations, an unbounded spin at ~200% CPU with flat RSS (measured on
    /// DRR509062: 40 min of CPU on a 100k-pair unit, no output). Retain only
    /// frees singletons, so re-checking on every insert buys nothing once the
    /// map is diverse. Evict lazily instead: when the map grows `cap/4` past
    /// the last eviction, drop to `cap` in one pass. Amortized O(1) per
    /// insert, same memory bound, no pathological case.
    pub fn insert(&mut self, key: K) {
        self.total += 1;
        let e = self.map.entry(key).or_insert(0);
        *e += 1;
        if self.map.len() > self.evict_at {
            self.map.retain(|_, c| *c > 1);
            self.evict_at = self.map.len().max(self.cap) + self.cap / 4;
        }
    }

    /// Exact sum of two counters, under the same bound: every key's count is added in one
    /// batch, then singletons are evicted if the map is over `evict_at` — the rule `insert`
    /// applies, applied once per merge rather than once per insert.
    pub fn merge(&mut self, o: &RepeatedCounter<K>) {
        self.total += o.total;
        for (k, &c) in &o.map {
            *self.map.entry(k.clone()).or_insert(0) += c;
        }
        if self.map.len() > self.evict_at {
            self.map.retain(|_, c| *c > 1);
            self.evict_at = self.map.len().max(self.cap) + self.cap / 4;
        }
    }

    /// Keys with `count >= min_count`, most frequent first.
    pub fn report(&self, min_count: u32) -> Vec<(K, u32)> {
        let mut v: Vec<(K, u32)> =
            self.map.iter().filter(|(_, &c)| c >= min_count).map(|(k, &c)| (k.clone(), c)).collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }
}

/// Pack a k-mer into a u64 (2 bits per base: A=00, C=01, G=10, T=11).
/// Returns `None` if the slice is empty or contains a non-ACGT base.
pub fn pack_kmer(seq: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &b in seq {
        let c = match b.to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        v = (v << 2) | c;
    }
    Some(v)
}

/// Unpack a packed k-mer back to a string (for display).
pub fn unpack_kmer(v: u64, k: usize) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::with_capacity(k);
    for i in (0..k).rev() {
        out.push(bases[((v >> (2 * i)) & 0b11) as usize]);
    }
    String::from_utf8(out).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_only_frequent_sequences() {
        let mut c = RepeatedCounter::new(100);
        let common = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let rare = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        for _ in 0..5 {
            c.insert(common.to_vec());
        }
        for _ in 0..2 {
            c.insert(rare.to_vec());
        }
        // 100 genuinely distinct singletons.
        for i in 0..100u8 {
            c.insert(vec![b'A', b'C', b'G', b'T', i]);
        }
        let rep = c.report(3);
        assert_eq!(rep.len(), 1, "only the 5x sequence clears the threshold");
        assert_eq!(rep[0].0, common);
        assert_eq!(rep[0].1, 5);
    }

    #[test]
    fn eviction_keeps_memory_bounded_and_preserves_frequent() {
        // A realistic cap: eviction fires only after ~100 distinct keys, by
        // which time the frequent sequence (every 10th read) has a count well
        // above 1, so it survives. A cap of 10 would evict it between its own
        // insertions — that is a pathological cap, not the shipped default.
        let mut c = RepeatedCounter::new(100);
        let common = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        for i in 0..1000u16 {
            c.insert(vec![b'A', b'C', b'G', b'T', (i >> 8) as u8, i as u8]);
            if i % 10 == 0 {
                c.insert(common.to_vec());
            }
        }
        assert!(c.len() <= 100, "the map must stay under the cap, got {}", c.len());
        let rep = c.report(50);
        assert!(rep.iter().any(|(k, _)| k == common), "the frequent sequence survives eviction");
    }

    #[test]
    fn kmer_pack_roundtrip() {
        let v = pack_kmer(b"ACGT").unwrap();
        assert_eq!(v, 0b00011011, "A=00 C=01 G=10 T=11 -> 0b00011011");
        assert_eq!(unpack_kmer(v, 4), "ACGT");
        assert_eq!(pack_kmer(b"ACGN"), None, "N cannot be packed");
    }

    /// The approximation, stated as a bound rather than hidden: a singleton evicted before
    /// its second occurrence costs exactly 1 count. This is why the counter is only used
    /// below the 0.1% reporting threshold, where 1 count cannot change a verdict.
    #[test]
    fn eviction_costs_at_most_one_count_per_evicted_singleton() {
        let mut c = RepeatedCounter::new(4);
        // 4 distinct singletons fill the cap; the eviction is amortized (fires at
        // cap + cap/4 = 5), so the 5th and 6th distinct keys are the two that push
        // the map past that threshold and trigger the retain(count > 1).
        for i in 0..4u8 {
            c.insert(vec![i]);
        }
        assert_eq!(c.len(), 4);
        c.insert(vec![9]); // 5th key: at the threshold, not over it -> no eviction yet
        assert_eq!(c.len(), 5);
        c.insert(vec![10]); // 6th key: over the threshold -> retain drops all five singletons
        // Every key was a singleton, so the eviction takes the new one too: the map is
        // empty, not "the just-inserted key survives". The bound is on the COST, not on
        // which key is kept — each evicted singleton had exactly one occurrence.
        assert_eq!(c.len(), 0, "every key was a singleton, so all are evicted");
        // The frequent key is counted exactly once it clears 1.
        let mut c = RepeatedCounter::new(4);
        for _ in 0..3 {
            c.insert(b"ACGT".to_vec());
        }
        for i in 0..4u8 {
            c.insert(vec![i]);
        }
        c.insert(vec![9]);
        let rep = c.report(2);
        assert_eq!(rep.len(), 1);
        assert_eq!(rep[0].1, 3, "the frequent key survives with its exact count");
    }
}
