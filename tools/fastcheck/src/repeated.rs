//! Exact 12-mer counting.
//!
//! k is fixed at 12, so the key space is 4¹² = 16 777 216 and a dense `Vec<u16>` of
//! that size is 32 MiB — **independent of the number of reads**, which is the property
//! the tool is sold on. There is no cap, no eviction and no sampling, so the spectrum is
//! exact and `distinct` (the number of 12-mers seen at least once) is countable from the
//! non-zero cells; neither was knowable from the bounded map this replaces
//! (`RepeatedCounter<u64>`, capped at 1e6 keys with singleton eviction, which also made
//! the block depend on read order and therefore on `-j`).
//!
//! The one bound left is the **width**, not the count: a cell is a `u16`, so a 12-mer
//! that occurs more than 65 535 times keeps its exact count in a small
//! `FxHashMap<u32, u32>` overflow (`hot`). Only the hottest k-mers of a library —
//! adapters, polyA — ever reach it, and the count stays exact; see the mean count per
//! 12-mer of ~167 on a 20 M-read cell, which is three orders of magnitude below the
//! ceiling for an ordinary k-mer.
//!
//! The table is allocated on first use rather than in `Accum::new`, so a gated run that
//! does not measure the k-mer module (`--only dup`) never pays for it.

use rustc_hash::FxHashMap;

/// Bases per k-mer. Fixed: the dense table's size is its key space.
pub const KMER_K: usize = 12;
/// 4^KMER_K: the number of cells in the dense table.
pub const KMER_CELLS: usize = 1 << (2 * KMER_K);

/// Exact 12-mer counter: a dense table plus a width-overflow map.
#[derive(Clone, Default)]
pub struct KmerTable {
    /// `cells[key]` = count of `key`, saturating at `u16::MAX` (the true total is then
    /// in `hot`). Empty until the first insert: a run that measures no k-mer allocates
    /// nothing.
    cells: Vec<u16>,
    /// Keys that exceeded `u16::MAX`, with their exact count.
    hot: FxHashMap<u32, u32>,
}

impl KmerTable {
    pub fn new() -> Self {
        Self::default()
    }

    /// One occurrence of a packed 12-mer (`pack_kmer`, 24 bits).
    #[inline]
    pub fn add(&mut self, key: u32) {
        debug_assert!((key as usize) < KMER_CELLS, "a 12-mer packs into 24 bits");
        if self.cells.is_empty() {
            self.cells = vec![0u16; KMER_CELLS];
        }
        let v = &mut self.cells[key as usize];
        if *v == u16::MAX {
            let e = self.hot.entry(key).or_insert(u16::MAX as u32);
            *e += 1;
        } else {
            *v += 1;
        }
    }

    /// Exact count of one key (0 if never seen).
    pub fn count(&self, key: u32) -> u32 {
        if self.cells.is_empty() {
            return 0;
        }
        self.total_at(key)
    }

    #[inline]
    fn total_at(&self, key: u32) -> u32 {
        let c = self.cells[key as usize];
        if c == u16::MAX {
            self.hot.get(&key).copied().unwrap_or(u16::MAX as u32)
        } else {
            c as u32
        }
    }

    /// Number of distinct 12-mers seen at least once, read off the non-zero cells.
    pub fn distinct(&self) -> u64 {
        self.cells.iter().filter(|&&c| c != 0).count() as u64
    }

    /// Resident bytes of the table, for the run's memory projection.
    pub fn resident_bytes(&self) -> usize {
        self.cells.len() * std::mem::size_of::<u16>()
            + self.hot.len() * (std::mem::size_of::<(u32, u32)>() + 1)
    }

    /// Exact sum of two tables. Cells that stay inside the `u16` width are summed in
    /// place; a cell that would overflow keeps its true total in `hot`, so the merged
    /// spectrum is the spectrum of the concatenated input — which is what makes the
    /// per-file accumulators mergeable and the block independent of `-j`.
    pub fn merge(&mut self, o: &KmerTable) {
        if o.cells.is_empty() {
            return;
        }
        if self.cells.is_empty() {
            self.cells = vec![0u16; KMER_CELLS];
        }
        for (i, &ov) in o.cells.iter().enumerate() {
            if ov == 0 {
                continue;
            }
            let k = i as u32;
            let sum = self.total_at(k) as u64 + o.total_at(k) as u64;
            if sum < u16::MAX as u64 {
                self.cells[i] = sum as u16;
                self.hot.remove(&k);
            } else {
                self.cells[i] = u16::MAX;
                self.hot.insert(k, sum as u32);
            }
        }
    }

    /// Keys with `count >= min_count`, most frequent first, ties broken by key so the
    /// order is a function of the data and not of insertion order.
    pub fn report(&self, min_count: u32) -> Vec<(u64, u32)> {
        let mut v: Vec<(u64, u32)> = Vec::new();
        for (i, &c) in self.cells.iter().enumerate() {
            if c == 0 || (c < u16::MAX && (c as u32) < min_count) {
                continue;
            }
            let n = self.total_at(i as u32);
            if n >= min_count {
                v.push((i as u64, n));
            }
        }
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
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
    fn kmer_pack_roundtrip() {
        let v = pack_kmer(b"ACGT").unwrap();
        assert_eq!(v, 0b00011011, "A=00 C=01 G=10 T=11 -> 0b00011011");
        assert_eq!(unpack_kmer(v, 4), "ACGT");
        assert_eq!(pack_kmer(b"ACGN"), None, "N cannot be packed");
    }

    /// The width bound, exercised: past 65 535 occurrences the cell saturates and the
    /// exact total moves to the overflow map. Nothing is lost and nothing is capped.
    #[test]
    fn the_table_counts_exactly_past_the_u16_width() {
        let mut t = KmerTable::new();
        let key = pack_kmer(b"ACGTACGTACGT").unwrap() as u32;
        for _ in 0..70_000 {
            t.add(key);
        }
        assert_eq!(t.count(key), 70_000, "the count is exact above u16::MAX");
        assert_eq!(t.distinct(), 1);
        assert_eq!(t.report(1), vec![(key as u64, 70_000)]);
    }

    /// The cap the dense table replaces was 1e6 keys with singleton eviction, so a
    /// spectrum with more distinct 12-mers than that lost its singletons entirely and
    /// could not report a distinct count at all. 1.1 M distinct keys, every one of them
    /// a singleton: all of them are still there.
    #[test]
    fn distinct_is_exact_past_the_old_one_million_key_cap() {
        let mut t = KmerTable::new();
        for i in 0..1_100_000u32 {
            t.add(i);
        }
        assert_eq!(t.distinct(), 1_100_000);
        assert!(t.report(2).is_empty(), "every key is a singleton");
        assert_eq!(t.count(0), 1);
        assert_eq!(t.count(1_099_999), 1);
    }

    /// Merge is exact across the width boundary in both directions: a hot key whose two
    /// halves each sit below the ceiling, and a cold key that must stay a plain sum.
    #[test]
    fn merge_is_exact_across_the_width_boundary() {
        let hot = pack_kmer(b"AAAAAAAAAAAA").unwrap() as u32;
        let cold = pack_kmer(b"CCCCCCCCCCCC").unwrap() as u32;
        let mut a = KmerTable::new();
        let mut b = KmerTable::new();
        for _ in 0..40_000 {
            a.add(hot);
            b.add(hot);
        }
        a.add(cold);
        b.add(cold);
        b.add(cold);
        a.merge(&b);
        assert_eq!(a.count(hot), 80_000, "the summed hot count is exact");
        assert_eq!(a.count(cold), 3);
        assert_eq!(a.distinct(), 2);
        assert_eq!(a.report(1), vec![(hot as u64, 80_000), (cold as u64, 3)]);

        // An empty side changes nothing, including the allocation.
        let empty = KmerTable::new();
        let mut c = a.clone();
        c.merge(&empty);
        assert_eq!(c.report(1), a.report(1));
    }

    /// Report order is data-determined: count descending, then key ascending, so two
    /// runs over the same input in a different order (or with a different `-j`) name the
    /// same k-mers in the same order.
    #[test]
    fn report_orders_by_count_then_key() {
        let mut t = KmerTable::new();
        for k in [9u32, 3, 7] {
            for _ in 0..2 {
                t.add(k);
            }
        }
        t.add(5);
        assert_eq!(t.report(1), vec![(3, 2), (7, 2), (9, 2), (5, 1)]);
        assert_eq!(t.report(2), vec![(3, 2), (7, 2), (9, 2)]);
    }
}
