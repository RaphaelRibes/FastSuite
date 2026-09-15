use rustc_hash::{FxHashMap, FxHashSet};
use std::hash::Hash;
use xxhash_rust::xxh3::{xxh3_64, xxh3_128};
use anyhow::Result;

/// A hash verifier that checks for sequence uniqueness using a contiguous memory buffer (HashSet).
pub struct HashVerifier<T>
where
    T: Hash + Eq,
{
    memory: FxHashSet<T>,
}

impl<T> HashVerifier<T>
where
    T: Hash + Eq,
{
    /// Initializes the structure with an estimated capacity to minimize reallocations.
    ///
    /// Uses `try_reserve` to detect OOM instead of panicking.
    pub fn new(estimated_capacity: usize) -> Result<Self> {
        let mut memory = FxHashSet::default();
        memory
            .try_reserve(estimated_capacity)
            .map_err(|_| anyhow::anyhow!("OOM: hash verifier allocation failed (capacity: {estimated_capacity})"))?;
        Ok(Self { memory })
    }

    /// Checks if the hash is already in memory.
    /// Returns `true` if the hash is newly inserted (unique sequence), `false` if it was already present (duplicate).
    pub fn verify(&mut self, hash: T) -> bool {
        self.memory.insert(hash)
    }
}

/// Tracks sequence abundances using a HashMap of hash → count.
/// Used by the `--abundance` mode to annotate output headers with `;size=N`.
pub struct AbundanceTracker<T: Hash + Eq> {
    counts: FxHashMap<T, u32>,
}

impl<T: Hash + Eq> AbundanceTracker<T> {
    /// Initializes the tracker with an estimated capacity to minimize reallocations.
    ///
    /// Uses `try_reserve` to detect OOM instead of panicking.
    pub fn new(capacity: usize) -> Result<Self> {
        let mut counts = FxHashMap::default();
        counts
            .try_reserve(capacity)
            .map_err(|_| anyhow::anyhow!("OOM: abundance tracker allocation failed (capacity: {capacity})"))?;
        Ok(Self { counts })
    }

    /// Records an occurrence of `hash` and returns `(is_first_occurrence, count_after_insert)`.
    pub fn track(&mut self, hash: T) -> (bool, u32) {
        let entry = self.counts.entry(hash).or_insert(0);
        *entry += 1;
        let count = *entry;
        (count == 1, count)
    }

    /// Returns the abundance count for `hash` without mutating the tracker.
    /// Use this in the second pass of abundance mode to retrieve the final count.
    #[inline]
    pub fn get_count(&self, hash: T) -> u32 {
        self.counts.get(&hash).copied().unwrap_or(0)
    }

    /// Returns the number of unique sequences tracked.
    pub fn len(&self) -> usize {
        self.counts.len()
    }
}

/// A trait for hashing single and paired sequences across multiple backends.
///
/// **Paired-end note**: `hash_pair` is order-dependent : `hash_pair(R1, R2) != hash_pair(R2, R1)`.
/// This is intentional: R1 and R2 carry distinct biological meaning.
/// If your library prep can swap mates, pre-sort them before hashing.
pub trait SequenceHasher: Hash + Eq + Copy + Send + Sync {
    /// Hashes a single sequence.
    fn hash_sequence(seq: &[u8]) -> Self;
    /// Hashes a pair of sequences (e.g., paired-end reads).
    /// Order matters: hash_pair(a, b) != hash_pair(b, a).
    fn hash_pair(seq1: &[u8], seq2: &[u8]) -> Self;
}

impl SequenceHasher for u64 {
    #[inline(always)]
    fn hash_sequence(seq: &[u8]) -> Self {
        xxh3_64(seq)
    }

    #[inline(always)]
    fn hash_pair(seq1: &[u8], seq2: &[u8]) -> Self {
        let h1 = xxh3_64(seq1);
        let h2 = xxh3_64(seq2);
        h1 ^ h2.rotate_left(32)
    }
}

impl SequenceHasher for u128 {
    #[inline(always)]
    fn hash_sequence(seq: &[u8]) -> Self {
        xxh3_128(seq)
    }

    #[inline(always)]
    fn hash_pair(seq1: &[u8], seq2: &[u8]) -> Self {
        let h1 = xxh3_128(seq1);
        let h2 = xxh3_128(seq2);
        h1 ^ h2.rotate_left(64)
    }
}

/// The type of hash size to use depending on estimated sequence count.
#[derive(Debug)]
pub enum HashType {
    XXH3_64,
    XXH3_128,
}

impl HashType {
    /// Returns the number of bits for the hash type (64 or 128).
    pub fn to_num(&self) -> usize {
        match self {
            HashType::XXH3_64 => 64,
            HashType::XXH3_128 => 128,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_sequence_deterministic() {
        let seq = b"ATCGATCG";
        assert_eq!(u64::hash_sequence(seq), u64::hash_sequence(seq));
        assert_eq!(u128::hash_sequence(seq), u128::hash_sequence(seq));
    }

    #[test]
    fn test_hash_pair_deterministic() {
        let r1 = b"ATCGATCG";
        let r2 = b"GCTAGCTA";
        assert_eq!(u64::hash_pair(r1, r2), u64::hash_pair(r1, r2));
        assert_eq!(u128::hash_pair(r1, r2), u128::hash_pair(r1, r2));
    }

    #[test]
    fn test_hash_pair_order_dependent() {
        let r1 = b"ATCGATCG";
        let r2 = b"GCTAGCTA";
        assert_ne!(u64::hash_pair(r1, r2), u64::hash_pair(r2, r1));
        assert_ne!(u128::hash_pair(r1, r2), u128::hash_pair(r2, r1));
    }

    #[test]
    fn test_hash_pair_identical_mates_still_valid() {
        let seq = b"ATCGATCG";
        // With the seeded approach, identical mates should still produce a meaningful hash
        let h = u64::hash_pair(seq, seq);
        assert_ne!(h, 0);
        let h128 = u128::hash_pair(seq, seq);
        assert_ne!(h128, 0);
    }

    #[test]
    fn test_verifier_detects_duplicate() {
        let mut v = HashVerifier::<u64>::new(16).unwrap();
        let h = u64::hash_sequence(b"ATCG");
        assert!(v.verify(h)); // first insert: unique
        assert!(!v.verify(h)); // second insert: duplicate
    }

    #[test]
    fn test_verifier_distinct_sequences() {
        let mut v = HashVerifier::<u64>::new(16).unwrap();
        assert!(v.verify(u64::hash_sequence(b"ATCG")));
        assert!(v.verify(u64::hash_sequence(b"GCTA")));
    }

    // --- AbundanceTracker tests ---

    #[test]
    fn test_abundance_tracker_tracks_unique() {
        let mut t = AbundanceTracker::<u64>::new(16).unwrap();
        let (is_first, count) = t.track(42);
        assert!(is_first);
        assert_eq!(count, 1);
        assert_eq!(t.len(), 1);
    }

    #[test]
    fn test_abundance_tracker_tracks_duplicate() {
        let mut t = AbundanceTracker::<u64>::new(16).unwrap();
        t.track(42);
        let (is_first, count) = t.track(42);
        assert!(!is_first);
        assert_eq!(count, 2);
        assert_eq!(t.len(), 1);
    }

    #[test]
    fn test_abundance_tracker_multiple_distinct() {
        let mut t = AbundanceTracker::<u64>::new(16).unwrap();
        assert_eq!(t.track(1).1, 1);
        assert_eq!(t.track(2).1, 1);
        assert_eq!(t.track(3).1, 1);
        assert_eq!(t.len(), 3);
    }

    #[test]
    fn test_abundance_tracker_varying_counts() {
        let mut t = AbundanceTracker::<u64>::new(16).unwrap();
        t.track(10);
        t.track(10);
        t.track(10);
        t.track(20);
        t.track(20);
        assert_eq!(t.track(10).1, 4);
        assert_eq!(t.track(20).1, 3);
        assert_eq!(t.track(30).1, 1);
        assert_eq!(t.len(), 3);
    }
}
