//! Over-represented sequences: the four modes of PLAN-QC §5.9.
//!
//! Exact identity needs one of:
//!
//! | mode | cost | report field |
//! |---|---|---|
//! | `none` (default) | zero | duplication levels still exact (§9) |
//! | `fingerprint` | +`unique x 8 B / 0.8` RAM | exact multiplicity, `names_available: false` |
//! | `second-pass` | one extra read of the file | `passes: 2` |
//! | `table` | 2-bit packed, `ceil(len/4)+8` B each | `unique_sequences_retained: true` |
//!
//! The table mode packs sequences 2 bits per base, so a 150 bp read costs 46 B including
//! the length header; 3e8 unique reads project ~11 GB and are refused above `--max-mem`
//! rather than OOM-killed at hour three. Second-pass keeps the fingerprint table in RAM
//! and re-reads the file only to recover the text of the fingerprints that cleared the
//! threshold, so its memory is the same as fingerprint mode.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{bail, Result};
use serde_json::json;

use crate::accum::Accum;
use crate::fingerprint;
use crate::reader;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverrepMode {
    None,
    Fingerprint,
    SecondPass,
    Table,
}

impl OverrepMode {
    pub fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "none" => OverrepMode::None,
            "fingerprint" => OverrepMode::Fingerprint,
            "second-pass" => OverrepMode::SecondPass,
            "table" => OverrepMode::Table,
            other => bail!(
                "--overrep-mode must be none, fingerprint, second-pass or table (got {other:?})"
            ),
        })
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            OverrepMode::None => "none",
            OverrepMode::Fingerprint => "fingerprint",
            OverrepMode::SecondPass => "second-pass",
            OverrepMode::Table => "table",
        }
    }

    /// Whether sequence text is recoverable in this mode.
    pub fn names_available(&self) -> bool {
        matches!(self, OverrepMode::SecondPass | OverrepMode::Table)
    }
}

/// 2-bit packed sequence with its length. `N` and every non-ACGT code are not packable and
/// force the `table` mode to refuse rather than to lose information silently.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Packed {
    len: u16,
    words: Vec<u64>,
}

impl Packed {
    pub fn pack(seq: &[u8]) -> Option<Self> {
        if seq.len() > u16::MAX as usize {
            return None;
        }
        let mut words = Vec::with_capacity(seq.len().div_ceil(32));
        let mut acc = 0u64;
        let mut n = 0u32;
        for &b in seq {
            let code = match b.to_ascii_uppercase() {
                b'A' => 0u64,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None, // not packable
            };
            acc |= code << (2 * n);
            n += 1;
            if n == 32 {
                words.push(acc);
                acc = 0;
                n = 0;
            }
        }
        if n > 0 {
            words.push(acc);
        }
        Some(Packed { len: seq.len() as u16, words })
    }

    pub fn unpack(&self) -> String {
        const BASES: [u8; 4] = *b"ACGT";
        let mut out = Vec::with_capacity(self.len as usize);
        for i in 0..self.len as usize {
            let w = self.words[i / 32];
            out.push(BASES[((w >> (2 * (i % 32))) & 0b11) as usize]);
        }
        String::from_utf8(out).unwrap_or_default()
    }

    /// Bytes this entry occupies: 2 bits per base plus the length and the map's overhead.
    pub fn bytes(&self) -> usize {
        2 + self.words.len() * 8 + 8 // len + payload + HashMap slot estimate
    }
}

/// The sequence table used by `table` mode: packed sequence -> multiplicity.
#[derive(Clone, Default)]
pub struct SequenceTable {
    counts: HashMap<Packed, u32>,
}

impl SequenceTable {
    pub fn observe(&mut self, seq: &[u8]) -> bool {
        match Packed::pack(seq) {
            Some(p) => {
                *self.counts.entry(p).or_insert(0) += 1;
                true
            }
            None => false, // a non-ACGT base: the caller must refuse, not lose data
        }
    }

    pub fn distinct(&self) -> u64 {
        self.counts.len() as u64
    }

    pub fn resident_bytes(&self) -> usize {
        self.counts.iter().map(|(k, _)| k.bytes()).sum()
    }

    /// Sequences occurring at least `min_count` times, most frequent first, as
    /// `(sequence text, count)`.
    pub fn at_least(&self, min_count: u32) -> Vec<(String, u32)> {
        let mut v: Vec<(String, u32)> = self
            .counts
            .iter()
            .filter(|(_, &c)| c >= min_count)
            .map(|(k, &c)| (k.unpack(), c))
            .collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }

    /// Exact sum of two tables (per-file tables are merged into the run total).
    pub fn merge(&mut self, o: &SequenceTable) {
        for (k, &c) in &o.counts {
            *self.counts.entry(k.clone()).or_insert(0) += c;
        }
    }
}

/// One entry of the over-representation report.
#[derive(Debug, Clone)]
pub struct OverrepEntry {
    pub seq: Option<String>,
    pub fingerprint: Option<String>,
    pub count: u32,
}

/// Build the over-representation report for one mode.
///
/// `acc` carries the exact fingerprint table. `files` and `threads` are needed only by
/// `second-pass`. `min_count` is the resolved threshold (flag, or max(10, 0.1% of reads)).
pub fn report(
    mode: OverrepMode,
    acc: &Accum,
    files: &[std::path::PathBuf],
    threads: usize,
    min_count: u32,
) -> Result<serde_json::Value> {
    let reads = acc.reads.max(1) as f64;
    match mode {
        OverrepMode::None => Ok(json!({
            "mode": "none",
            "names_available": false,
            "unique_sequences_retained": false,
            "min_count": min_count,
            "sequences": [],
            "note": "duplication levels are exact in this mode; over-represented sequence \
                     identity needs --overrep-mode fingerprint, second-pass or table",
        })),
        OverrepMode::Fingerprint => {
            let top = acc.seq_counts.at_least(min_count);
            let seqs: Vec<serde_json::Value> = top
                .iter()
                .take(20)
                .map(|(fp, c)| {
                    json!({
                        "fingerprint": format!("{fp:016x}"),
                        "count": c,
                        "frac": (*c as f64 / reads * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            Ok(json!({
                "mode": "fingerprint",
                "names_available": false,
                "unique_sequences_retained": false,
                "min_count": min_count,
                "sequences": seqs,
                "note": "multiplicities are exact; the sequence text is not recoverable from \
                         a 64-bit fingerprint. Use --overrep-mode table or second-pass to \
                         name them.",
            }))
        }
        OverrepMode::Table => {
            // The table was filled during the scan (accum.rs), so this mode only reports.
            let top = acc.overrep_table.at_least(min_count);
            let seqs: Vec<serde_json::Value> = top
                .iter()
                .take(20)
                .map(|(s, c)| {
                    json!({
                        "seq": s,
                        "count": c,
                        "frac": (*c as f64 / reads * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            Ok(json!({
                "mode": "table",
                "names_available": true,
                "unique_sequences_retained": true,
                "min_count": min_count,
                "unique_sequences": acc.overrep_table.distinct(),
                "sequences": seqs,
            }))
        }
        OverrepMode::SecondPass => {
            let wanted: std::collections::HashSet<u64> = acc
                .seq_counts
                .at_least(min_count)
                .into_iter()
                .map(|(fp, _)| fp)
                .collect();
            if wanted.is_empty() {
                return Ok(json!({
                    "mode": "second-pass",
                    "names_available": true,
                    "unique_sequences_retained": false,
                    "passes": 2,
                    "min_count": min_count,
                    "sequences": [],
                    "note": "no fingerprint cleared the threshold, so the second read found \
                             nothing to name",
                }));
            }
            // Second read: recover the text of the fingerprints above the threshold.
            let mut texts: HashMap<u64, String> = HashMap::new();
            for path in files {
                second_pass(path, threads, &wanted, &mut texts)?;
            }
            let mut seqs: Vec<serde_json::Value> = acc
                .seq_counts
                .at_least(min_count)
                .iter()
                .take(20)
                .map(|(fp, c)| {
                    json!({
                        "seq": texts.get(fp).cloned().unwrap_or_default(),
                        "count": c,
                        "frac": (*c as f64 / reads * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            seqs.sort_by(|a, b| {
                b["count"].as_u64().unwrap_or(0).cmp(&a["count"].as_u64().unwrap_or(0))
            });
            Ok(json!({
                "mode": "second-pass",
                "names_available": true,
                "unique_sequences_retained": false,
                "passes": 2,
                "min_count": min_count,
                "sequences": seqs,
            }))
        }
    }
}

/// Re-read one file and record the text of every sequence whose fingerprint is wanted.
///
/// This is the ONLY place in fastcheck that reads a file twice, and it is opt-in: the
/// report carries `passes: 2` so a consumer never mistakes the cost for the default.
fn second_pass(
    path: &Path,
    gz_threads: usize,
    wanted: &std::collections::HashSet<u64>,
    texts: &mut HashMap<u64, String>,
) -> Result<()> {
    let (mut rdr, _counter, _mode) = reader::open(path, gz_threads)?;
    let (mut line, mut seq, mut plus, mut qual) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    loop {
        if reader::read_line(&mut rdr, &mut line)? == 0 {
            break;
        }
        if line.is_empty() {
            continue;
        }
        reader::read_line(&mut rdr, &mut seq)?;
        reader::read_line(&mut rdr, &mut plus)?;
        reader::read_line(&mut rdr, &mut qual)?;
        if seq.is_empty() {
            continue;
        }
        let fp = fingerprint::fingerprint(&seq);
        if wanted.contains(&fp) {
            texts.entry(fp).or_insert_with(|| String::from_utf8_lossy(&seq).to_string());
        }
    }
    Ok(())
}

/// Projected bytes for `table` mode at a given unique-sequence count and read length.
/// `ceil(len/4) + 8` per sequence, times the hashbrown load factor of 8/7.
pub fn table_projection(unique: usize, read_len: usize) -> u64 {
    let per = (read_len.div_ceil(4) + 8) as u64;
    let need = (unique as u64 * 8).div_ceil(7);
    need.next_power_of_two().saturating_mul(per + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mode_parsing_rejects_unknown_values() {
        assert_eq!(OverrepMode::parse("none").unwrap(), OverrepMode::None);
        assert_eq!(OverrepMode::parse("second-pass").unwrap(), OverrepMode::SecondPass);
        let e = OverrepMode::parse("sampled").unwrap_err().to_string();
        assert!(e.contains("second-pass"), "{e}");
    }

    #[test]
    fn packed_roundtrip_and_non_acgt_refusal() {
        // N is not packable: the caller must refuse, not lose data.
        assert!(Packed::pack(b"ACGTN").is_none());
        assert!(Packed::pack(b"ACGTACGTNN").is_none());
        let p = Packed::pack(b"ACGTACGTAC").unwrap();
        assert_eq!(p.unpack(), "ACGTACGTAC");
        // 10 bases -> 1 word (32 bases fit in a u64) + 2 B length + 8 B slot.
        assert_eq!(p.bytes(), 2 + 8 + 8);
        // A 33-base read needs two words.
        let q = Packed::pack(&b"A".repeat(33)).unwrap();
        assert_eq!(q.bytes(), 2 + 2 * 8 + 8);
    }

    #[test]
    fn sequence_table_counts_and_names() {
        let mut t = SequenceTable::default();
        for _ in 0..4 {
            t.observe(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        }
        t.observe(b"TTTTAAAACCCCGGGG");
        assert_eq!(t.distinct(), 2);
        let top = t.at_least(3);
        assert_eq!(top.len(), 1);
        assert_eq!(top[0].0, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        assert_eq!(top[0].1, 4);
    }

    /// PLAN-QC §5.9 quotes ~11 GB for 3e8 unique reads at 150 bp — that is the packed
    /// payload alone (⌈150/4⌉ = 38 B). The projection adds the per-entry slot (8 B) and
    /// hashbrown's 8/7 load factor, which is what actually gets allocated, so it is ~2x
    /// the quote: the quote is a payload figure, this is a resident figure. Both are
    /// pinned here so neither silently changes.
    #[test]
    fn table_projection_matches_the_plan_quote() {
        // Payload-only, as PLAN-QC quotes it: ⌈len/4⌉ B per unique read.
        let payload = (300_000_000u64 * 150u64.div_ceil(4)) as f64 / (1u64 << 30) as f64;
        assert!((payload - 11.0).abs() < 1.5, "payload {payload:.2} GB, PLAN-QC says ~11 GB");
        // Resident, as allocated: the slot, the length header and the load factor.
        let bytes = table_projection(300_000_000, 150);
        let gb = bytes as f64 / (1u64 << 30) as f64;
        assert!((gb - 23.0).abs() < 2.0, "projected {gb:.2} GB");
    }
}
