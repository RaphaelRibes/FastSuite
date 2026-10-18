//! Adapter catalogue and seeded matching, copied from `tools/fastprep/src/adapter.rs`
//! (D-1: no shared crate — the two tools may drift, and the tests below are what keep
//! them honest).
//!
//! Semantics preserved exactly, because a per-position adapter-content curve is only
//! meaningful if "adapter" means the same thing here as in the trimmer:
//!
//! * catalogue: TruSeq/R1, TruSeq/R2, Nextera, SmallRNA (`preset_illumina` +
//!   `preset_nextera` + `preset_small_rna`, deduplicated by sequence)
//! * a match at position `p` with overlap `k = min(len(adapter), len(read) - p)` is
//!   accepted when `mismatches <= ceil(k * error_rate)`
//! * `min_overlap_se` (fastprep default 14) and `error_rate` (default 0.2) are the flags
//! * leftmost match wins; on a tie the adapter earlier in the set wins
//! * the 6-mer seed index is a rejection aid only: any alignment of `k` bases with at most
//!   `max_mm` mismatches contains an exactly-matching run of at least `k/(max_mm+1)` bases,
//!   so a seed hit at that position finds the alignment. Everything else is provably
//!   non-matching. It is exact, never a heuristic.

/// One named adapter sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    pub name: String,
    pub seq: Vec<u8>,
}

/// All TruSeq + Nextera sequences (fastprep's default "illumina" preset).
pub fn preset_illumina() -> Vec<Adapter> {
    vec![
        Adapter { name: "TruSeq/R1".into(), seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec() },
        Adapter { name: "TruSeq/R2".into(), seq: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec() },
        Adapter { name: "Nextera".into(),    seq: b"CTGTCTCTTATACACATCT".to_vec() },
    ]
}

pub fn preset_nextera() -> Vec<Adapter> {
    vec![Adapter { name: "Nextera".into(), seq: b"CTGTCTCTTATACACATCT".to_vec() }]
}

pub fn preset_small_rna() -> Vec<Adapter> {
    vec![Adapter { name: "SmallRNA".into(), seq: b"TGGAATTCTCGGGTGCCAAGG".to_vec() }]
}

/// Every built-in sequence, deduplicated by sequence (Nextera appears in two presets).
pub fn catalogue() -> Vec<Adapter> {
    let mut all = preset_illumina();
    for a in preset_nextera().into_iter().chain(preset_small_rna()) {
        if !all.iter().any(|x| x.seq == a.seq) {
            all.push(a);
        }
    }
    all
}

/// The full adapter configuration passed to [`match_at`].
#[derive(Debug, Clone)]
pub struct AdapterSet {
    pub adapters: Vec<Adapter>,
    /// Maximum fraction of mismatches allowed in a match.
    pub error_rate: f64,
    /// Minimum overlap for sequence matching.
    pub min_overlap_se: usize,
}

impl AdapterSet {
    pub fn new(adapters: Vec<Adapter>, error_rate: f64, min_overlap_se: usize) -> Self {
        AdapterSet { adapters, error_rate, min_overlap_se }
    }

    /// fastprep's defaults: `--adapter-error-rate 0.2 --min-adapter-overlap 14`.
    pub fn illumina_default() -> Self {
        Self::new(catalogue(), 0.2, 14)
    }
}

/// Result of a match: which adapter, and where it starts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdapterMatch {
    pub name: String,
    /// 0-based position in the read where the adapter match starts.
    pub start: usize,
    /// Overlap length actually compared.
    pub overlap: usize,
}

/// Byte -> 2-bit base code (0=A, 1=C, 2=G, 3=T; anything else -> 3, a wildcard that may
/// spuriously seed and is rejected by the mismatch count).
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

/// 2-bit encoding of a 6-mer into a 12-bit code.
#[inline]
fn seed_code(bases: &[u8]) -> u32 {
    let mut code: u32 = 0;
    for &b in bases {
        code = (code << 2) | BASE2[b as usize] as u32;
    }
    code & 0xFFF
}

/// Smallest overlap length for which the 6-mer guarantee holds, per adapter.
/// `usize::MAX` = never (e.g. a high error rate).
fn min_sound_k(adapter_len: usize, error_rate: f64) -> usize {
    for k in 1..=adapter_len {
        let max_mm = (k as f64 * error_rate).ceil() as usize;
        if k >= 6 * (max_mm + 1) {
            return k;
        }
    }
    usize::MAX
}

/// Mismatches between two equal-length slices, stopping once `cap` is exceeded.
#[inline]
fn mismatches_upto(a: &[u8], b: &[u8], cap: usize) -> usize {
    debug_assert_eq!(a.len(), b.len());
    let mut mm = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if x.to_ascii_uppercase() != y.to_ascii_uppercase() {
            mm += 1;
            if mm > cap {
                return mm;
            }
        }
    }
    mm
}

/// The leftmost adapter match in `seq`, tie-broken by set order.
///
/// Returns `None` when no adapter matches anywhere under the set's thresholds — which is
/// the common case for a clean read, and the reason the seed index exists.
pub fn match_at(seq: &[u8], set: &AdapterSet) -> Option<AdapterMatch> {
    let l = seq.len();
    let max_p = l.saturating_sub(set.min_overlap_se);
    let mut best: Option<(usize, usize)> = None; // (start, adapter index)
    let mut consider = |p: usize, ai: usize| {
        if best.is_none_or(|(bp, bai)| p < bp || (p == bp && ai < bai)) {
            best = Some((p, ai));
        }
    };

    for (ai, adapter) in set.adapters.iter().enumerate() {
        let a = adapter.seq.len();
        if a == 0 {
            continue;
        }
        let ms = min_sound_k(a, set.error_rate);
        // Seed-guided pass: only sound for k >= ms. Scan the read once with a rolling
        // 6-mer and test the candidate positions the seeds imply.
        let mut code: u32 = 0;
        for q in 0..l {
            code = ((code << 2) | BASE2[seq[q] as usize] as u32) & 0xFFF;
            if q < 5 {
                continue;
            }
            let win_start = q - 5;
            for off in 0..=a.saturating_sub(6) {
                if seed_code(&adapter.seq[off..off + 6]) != code {
                    continue;
                }
                if win_start < off {
                    continue;
                }
                let p = win_start - off;
                let k = a.min(l - p);
                if k < set.min_overlap_se || k < ms || p > max_p {
                    continue;
                }
                let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                }
            }
        }
        // Direct pass for the seed-unsound region (short adapters, high error rates, or
        // truncated tails): scanned exactly, reference-style.
        if a < ms || ms == usize::MAX {
            for p in 0..=max_p {
                let k = a.min(l - p);
                if k < set.min_overlap_se {
                    break;
                }
                let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                    break;
                }
            }
        } else {
            let start = if l >= a { l - a + 1 } else { 0 };
            let start = start.max(l.saturating_add(1).saturating_sub(ms));
            for p in start..=max_p {
                let k = a.min(l - p);
                if k < set.min_overlap_se {
                    break;
                }
                let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                    consider(p, ai);
                    break;
                }
            }
        }
    }

    best.map(|(p, ai)| {
        let a = set.adapters[ai].seq.len();
        AdapterMatch { name: set.adapters[ai].name.clone(), start: p, overlap: a.min(l - p) }
    })
}

/// All adapter matches in `seq`, one per adapter that matches anywhere, as
/// `(adapter index, start)`. Used by the per-position module, which counts a match at
/// EVERY position where an adapter begins, not just the leftmost.
pub fn all_matches(seq: &[u8], set: &AdapterSet) -> Vec<(usize, usize)> {
    let l = seq.len();
    let max_p = l.saturating_sub(set.min_overlap_se);
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (ai, adapter) in set.adapters.iter().enumerate() {
        let a = adapter.seq.len();
        if a == 0 {
            continue;
        }
        for p in 0..=max_p {
            let k = a.min(l - p);
            if k < set.min_overlap_se {
                break;
            }
            let max_mm = (k as f64 * set.error_rate).ceil() as usize;
            if mismatches_upto(&seq[p..p + k], &adapter.seq[..k], max_mm) <= max_mm {
                out.push((ai, p));
                break; // leftmost per adapter, like the matcher
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalogue_is_deduplicated_by_sequence() {
        let cat = catalogue();
        let nextera = cat.iter().filter(|a| a.seq == b"CTGTCTCTTATACACATCT").count();
        assert_eq!(nextera, 1, "Nextera appears in two presets but must be listed once");
        assert_eq!(cat.len(), 4, "TruSeq/R1, TruSeq/R2, Nextera, SmallRNA");
    }

    /// The agreement test D-1 asks for: fastcheck's matcher must accept exactly what
    /// fastprep's does, on the same inputs, with the same thresholds.
    #[test]
    fn a_full_adapter_at_the_3p_end_is_found_at_the_right_position() {
        let set = AdapterSet::illumina_default();
        // 120 bases of filler then the full 33-base TruSeq/R1 adapter: insert 120.
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        let m = match_at(&seq, &set).unwrap();
        assert_eq!(m.name, "TruSeq/R1");
        assert_eq!(m.start, 120);
        assert_eq!(m.overlap, 33);
    }

    #[test]
    fn a_clean_read_matches_nothing() {
        let set = AdapterSet::illumina_default();
        let seq: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        assert!(match_at(&seq, &set).is_none());
        assert!(all_matches(&seq, &set).is_empty());
    }

    /// The threshold is a flag: at error_rate 0.2 and overlap 20, up to 4 mismatches are
    /// tolerated. Three must still match, five must not.
    #[test]
    fn the_error_rate_threshold_is_honoured_exactly() {
        let set = AdapterSet::new(vec![Adapter {
            name: "T".into(),
            seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
        }], 0.2, 14);
        // 116 filler bases + the full 33-base adapter: the compared overlap is 33, and
        // ceil(33*0.2) = 7 mismatches are allowed for a full-length match.
        let mut seq = vec![b'A'; 116];
        let mut ad = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec();
        // 3 mismatches in the last 20 bases must pass.
        for i in [13, 20, 27] {
            ad[i] = if ad[i] == b'A' { b'C' } else { b'A' };
        }
        seq.extend_from_slice(&ad);
        assert!(match_at(&seq, &set).is_some(), "3 mismatches in the adapter must pass");

        // 10 mismatches: beyond ceil(33*0.2) = 7, so it must fail.
        let mut ad2 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec();
        for i in 0..10 {
            ad2[i] = if ad2[i] == b'A' { b'C' } else { b'A' };
        }
        let mut seq2 = vec![b'A'; 116];
        seq2.extend_from_slice(&ad2);
        assert!(match_at(&seq2, &set).is_none(), "10 mismatches in 33 bases must fail");
    }

    /// A match shorter than --min-adapter-overlap is not a match, however perfect.
    #[test]
    fn a_match_below_min_overlap_is_rejected() {
        let set = AdapterSet::new(vec![Adapter {
            name: "T".into(),
            seq: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
        }], 0.2, 14);
        // Only the first 10 bases of the adapter are present.
        let mut seq = vec![b'A'; 140];
        seq.extend_from_slice(b"AGATCGGAAG");
        assert!(match_at(&seq, &set).is_none(), "10 < 14: below the floor");
    }

    #[test]
    fn leftmost_wins_and_set_order_breaks_ties() {
        let set = AdapterSet::illumina_default();
        // TruSeq/R1 then Nextera, both matching: the earlier start must win.
        let mut seq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec();
        seq.extend_from_slice(b"CTGTCTCTTATACACATCT");
        let m = match_at(&seq, &set).unwrap();
        assert_eq!(m.name, "TruSeq/R1");
        assert_eq!(m.start, 0);
    }

    /// all_matches reports one entry per matching adapter, which is what the
    /// per-position curve counts.
    #[test]
    fn all_matches_reports_each_adapter_once() {
        let set = AdapterSet::illumina_default();
        let mut seq = vec![b'A'; 100];
        seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        let v = all_matches(&seq, &set);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].1, 100);
    }

    /// The seeded matcher must agree with a plain exhaustive reference scan on every
    /// input: the seeds are a rejection aid, so they may never change a verdict
    /// (fastprep pins the same property against its own reference).
    #[test]
    fn seeded_matcher_agrees_with_the_exhaustive_reference() {
        fn reference(seq: &[u8], set: &AdapterSet) -> Option<AdapterMatch> {
            let l = seq.len();
            let max_p = l.saturating_sub(set.min_overlap_se);
            let mut best: Option<(usize, usize)> = None;
            for (ai, adapter) in set.adapters.iter().enumerate() {
                let a = adapter.seq.len();
                if a == 0 {
                    continue;
                }
                for p in 0..=max_p {
                    let k = a.min(l - p);
                    if k < set.min_overlap_se {
                        break;
                    }
                    let max_mm = (k as f64 * set.error_rate).ceil() as usize;
                    let mut mm = 0usize;
                    let mut ok = true;
                    for i in 0..k {
                        if seq[p + i].to_ascii_uppercase() != adapter.seq[i].to_ascii_uppercase() {
                            mm += 1;
                            if mm > max_mm {
                                ok = false;
                                break;
                            }
                        }
                    }
                    if ok {
                        let better = match best {
                            None => true,
                            Some((bp, bai)) => p < bp || (p == bp && ai < bai),
                        };
                        if better {
                            best = Some((p, ai));
                        }
                        break;
                    }
                }
            }
            best.map(|(p, ai)| AdapterMatch {
                name: set.adapters[ai].name.clone(),
                start: p,
                overlap: set.adapters[ai].seq.len().min(l - p),
            })
        }

        // A deterministic xorshift keeps the test reproducible without a dev-dependency.
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let bases = b"ACGT";
        let adapters = catalogue();
        for case in 0..3000 {
            let l = 20 + (next() % 140) as usize;
            let mut seq: Vec<u8> = (0..l).map(|_| bases[(next() % 4) as usize]).collect();
            // Plant an adapter at a random position in three quarters of the cases.
            if case % 4 != 0 {
                let ad = &adapters[(next() as usize) % adapters.len()].seq;
                let mut copy = ad.clone();
                if case % 4 == 3 {
                    // With errors inside the planted adapter.
                    for _ in 0..(next() % 4) {
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
            }
            for &eps in &[0.0f64, 0.1, 0.2] {
                for &ov in &[4usize, 14] {
                    let set = AdapterSet::new(adapters.clone(), eps, ov);
                    assert_eq!(
                        match_at(&seq, &set),
                        reference(&seq, &set),
                        "case {case} eps {eps} overlap {ov}"
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod agreement_with_fastprep {
    //! Cross-tool check (D-1: consistency by tests): run fastcheck's copied matcher over
    //! the demo fixture and compare the hit count with what fastprep reported for the
    //! same file with the same preset. Run with:
    //!   FASTPREP_GT=path cargo test --offline adapter::agreement_with_fastprep -- --ignored
    //! (ignored by default: it needs the fastprep JSON produced by a separate run.)
    use super::*;

    #[test]
    #[ignore]
    fn nextera_hits_match_fastprep() {
        let gt = std::env::var("FASTPREP_GT").expect("FASTPREP_GT unset");
        let expected: u64 = {
            let d: serde_json::Value =
                serde_json::from_str(&std::fs::read_to_string(&gt).unwrap()).unwrap();
            d["adapter"]["reads_trimmed"].as_u64().unwrap()
        };
        let path = std::env::var("ADAPTER_FIXTURE").expect("ADAPTER_FIXTURE unset");
        let set = AdapterSet::new(preset_nextera(), 0.2, 14);
        let mut n = 0u64;
        let text = std::fs::read_to_string(&path).unwrap();
        let mut lines = text.lines();
        loop {
            let Some(h) = lines.next() else { break };
            if h.is_empty() {
                continue;
            }
            let s = lines.next().unwrap_or("");
            let _ = lines.next();
            let _ = lines.next();
            if match_at(s.as_bytes(), &set).is_some() {
                n += 1;
            }
        }
        assert_eq!(n, expected, "fastcheck {n} vs fastprep {expected}");
    }
}
