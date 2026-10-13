//! Adapter content per position (PLAN-QC module 11).
//!
//! For each adapter in the catalogue and each position `p`, the exact count of reads whose
//! adapter match STARTS at `p`. This is a different question from "how many reads carry an
//! adapter": the per-position curve shows WHERE the readthrough begins, which is what
//! distinguishes a short-insert library (a spike at one position) from adapter
//! contamination spread across the read.
//!
//! The matcher semantics are fastprep's, copied (D-1) — same catalogue, same
//! `--adapter-error-rate` / `--min-adapter-overlap` defaults, same leftmost-wins rule — so
//! the curve and the trimmer cannot disagree about what an adapter is.

use crate::adapter::{self, AdapterSet};

/// `per_adapter[ai][pos]` = reads whose match for adapter `ai` starts at `pos`.
#[derive(Clone)]
pub struct AdapterPos {
    pub names: Vec<String>,
    pub counts: Vec<Vec<u64>>,
    /// Reads carrying at least one match, per adapter.
    pub reads_with_match: Vec<u64>,
    /// Total match observations (a read can match two adapters).
    pub total_matches: u64,
}

impl AdapterPos {
    pub fn new(set: &AdapterSet) -> Self {
        Self {
            names: set.adapters.iter().map(|a| a.name.clone()).collect(),
            // A match can start anywhere in the read, so the per-position rows grow with
            // the positions actually observed (lazily, like the accumulator's rows).
            counts: set.adapters.iter().map(|_| Vec::new()).collect(),
            reads_with_match: vec![0; set.adapters.len()],
            total_matches: 0,
        }
    }

    #[inline]
    pub fn observe(&mut self, seq: &[u8], set: &AdapterSet) {
        for (ai, p) in adapter::all_matches(seq, set) {
            if self.counts[ai].len() <= p {
                self.counts[ai].resize(p + 1, 0);
            }
            self.counts[ai][p] += 1;
            self.reads_with_match[ai] += 1;
            self.total_matches += 1;
        }
    }

    /// Exact sum of two counters.
    pub fn merge(&mut self, o: &AdapterPos) {
        for (ai, row) in o.counts.iter().enumerate() {
            if row.len() > self.counts[ai].len() {
                self.counts[ai].resize(row.len(), 0);
            }
            for (d, s) in self.counts[ai].iter_mut().zip(row.iter()) {
                *d += s;
            }
            self.reads_with_match[ai] += o.reads_with_match[ai];
        }
        self.total_matches += o.total_matches;
    }

    /// The report block: one entry per adapter with its non-zero per-position counts.
    pub fn to_json(&self, error_rate: f64, min_overlap: usize) -> serde_json::Value {
        let adapters: Vec<serde_json::Value> = self
            .names
            .iter()
            .enumerate()
            .map(|(ai, name)| {
                let mut per_pos = serde_json::Map::new();
                for (p, &c) in self.counts[ai].iter().enumerate() {
                    if c > 0 {
                        per_pos.insert((p + 1).to_string(), serde_json::json!(c));
                    }
                }
                serde_json::json!({
                    "name": name,
                    "reads_with_match": self.reads_with_match[ai],
                    "per_position": serde_json::Value::Object(per_pos),
                })
            })
            .collect();
        serde_json::json!({
            "error_rate": error_rate,
            "min_overlap": min_overlap,
            "total_matches": self.total_matches,
            "adapters": adapters,
        })
    }

    /// Extra positions-TSV columns: `adapter_<name>` per adapter, 1-based positions.
    pub fn tsv_columns(&self) -> Vec<(String, Vec<u64>)> {
        self.names
            .iter()
            .enumerate()
            .map(|(ai, name)| (format!("adapter_{name}"), self.counts[ai].clone()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set() -> AdapterSet {
        AdapterSet::new(adapter::preset_nextera(), 0.2, 14)
    }

    #[test]
    fn counts_the_start_position_of_each_match() {
        let mut c = AdapterPos::new(&set());
        // 120 filler bases then the full 19-base Nextera adapter.
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        c.observe(&seq, &set());
        assert_eq!(c.reads_with_match[0], 1);
        assert_eq!(c.counts[0][120], 1, "the match starts at 0-based 120 = position 121");
        assert_eq!(c.counts[0][0], 0);
    }

    #[test]
    fn a_clean_read_contributes_nothing() {
        let mut c = AdapterPos::new(&set());
        let seq: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        c.observe(&seq, &set());
        assert_eq!(c.total_matches, 0);
        assert!(c.counts[0].iter().all(|&c| c == 0));
    }

    #[test]
    fn merge_is_exact() {
        let mut a = AdapterPos::new(&set());
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        a.observe(&seq, &set());
        let mut b = AdapterPos::new(&set());
        b.observe(&seq, &set());
        a.merge(&b);
        assert_eq!(a.counts[0][120], 2);
        assert_eq!(a.reads_with_match[0], 2);
        assert_eq!(a.total_matches, 2);
    }

    #[test]
    fn to_json_uses_one_based_positions_and_omits_zeros() {
        let mut c = AdapterPos::new(&set());
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        c.observe(&seq, &set());
        let j = c.to_json(0.2, 14);
        let per = j["adapters"][0]["per_position"].as_object().unwrap();
        assert_eq!(per.len(), 1, "only the non-zero position");
        assert_eq!(per["121"], 1);
        assert_eq!(j["adapters"][0]["name"], "Nextera");
    }

    /// A truncated tail shorter than --min-adapter-overlap is not a match, so it does not
    /// appear at its position (the flag, honoured exactly).
    #[test]
    fn a_truncated_match_below_the_overlap_floor_is_not_counted() {
        let mut c = AdapterPos::new(&set());
        let seq = vec![b'A'; 140];
        c.observe(&seq, &set());
        assert_eq!(c.total_matches, 0);
        // 16 bases of adapter at the 3' end, floor 14: it IS a match at 0-based 134.
        let mut seq2 = vec![b'A'; 134];
        seq2.extend_from_slice(b"CTGTCTCTTATACACA");
        c.observe(&seq2, &set());
        assert_eq!(c.total_matches, 1);
        assert_eq!(c.counts[0][134], 1);
    }
}
