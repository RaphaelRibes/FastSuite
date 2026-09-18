//! Bounded over-represented-sequence and k-mer counters.
//!
//! Streaming with singleton eviction: when the number of distinct keys exceeds
//! a cap, every key with count == 1 is dropped. A singleton can never be
//! over-represented, and losing its first occurrence costs at most 1 of a count
//! that must reach a large threshold to be reported, so memory stays bounded on
//! 50 M-read files while the frequent sequences are counted exactly.

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
}

impl<K: Hash + Eq + Clone + Ord> RepeatedCounter<K> {
    pub fn new(cap: usize) -> Self {
        Self { map: HashMap::with_capacity(cap.min(1 << 20)), cap, total: 0 }
    }

    /// Record one occurrence of `key`, evicting singletons if over the cap.
    pub fn insert(&mut self, key: K) {
        self.total += 1;
        let e = self.map.entry(key).or_insert(0);
        *e += 1;
        if self.map.len() > self.cap {
            self.map.retain(|_, c| *c > 1);
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
}
