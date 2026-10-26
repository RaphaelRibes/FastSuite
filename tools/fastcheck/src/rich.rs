//! The rich-report content: polyX, over-represented sequences, k-mer spectrum
//! and heatmap, as both `metrics` (numeric, flows through `fasthelper
//! concatenate report`) and tool-specific blocks (the full distributions/grids,
//! carried in the JSON and HTML).

use serde_json::{json, Map, Value};

use crate::accum::Accum;
use crate::heatmap;
use crate::repeated;

/// Round a fraction to 4 decimals (the report's convention).
fn r4(x: f64) -> f64 {
    (x * 10000.0).round() / 10000.0
}

/// polyX metrics + the per-base 3′/5′ run-length block.
pub fn polyx_metrics(acc: &Accum) -> (Map<String, Value>, Value) {
    let reads = acc.reads.max(1) as f64;
    let (mut p3max, mut p5max) = (0usize, 0usize);
    let (mut p3gt5, mut p5gt1) = (0u64, 0u64);
    for b in 0..4 {
        for (len, &c) in acc.polyx3[b].iter().enumerate() {
            if c > 0 {
                p3max = p3max.max(len);
                if len > 5 {
                    p3gt5 += c;
                }
            }
        }
        for (len, &c) in acc.polyx5[b].iter().enumerate() {
            if c > 0 {
                p5max = p5max.max(len);
                if len > 1 {
                    p5gt1 += c;
                }
            }
        }
    }
    let mut m = Map::new();
    m.insert("polyx3_max".into(), json!(p3max));
    m.insert("polyx5_max".into(), json!(p5max));
    m.insert("polyx3_frac_gt_5".into(), json!(r4(p3gt5 as f64 / reads)));
    m.insert("polyx5_frac_gt_1".into(), json!(r4(p5gt1 as f64 / reads)));

    let names = ["A", "C", "G", "T"];
    let mut b3 = Map::new();
    let mut b5 = Map::new();
    for (bi, name) in names.iter().enumerate() {
        b3.insert((*name).into(), crate::report::hist_json(&acc.polyx3[bi]));
        b5.insert((*name).into(), crate::report::hist_json(&acc.polyx5[bi]));
    }
    (m, json!({ "3p": b3, "5p": b5 }))
}

/// The 2-colour first-cycle G dropout: the position-1 G fraction and the
/// library-G-baseline minus it (positive = a real dropout).
pub fn head_dropout(acc: &Accum) -> (f64, f64) {
    let g1 = if let Some(row) = acc.base_pos.first() {
        let depth: u64 = row.iter().sum();
        if depth > 0 {
            row[2] as f64 / depth as f64
        } else {
            0.25
        }
    } else {
        0.25
    };
    let g: u64 = acc.base_pos.iter().map(|r| r[2]).sum();
    let baseline = g as f64 / acc.bases.max(1) as f64;
    (g1, baseline - g1)
}

/// Duplication levels: the exact multiplicity histogram, plus the duplication level and the
/// distinct/read ratio. Over-represented sequences are a separate concern handled by
/// `crate::overrep` (§5.9) — this function only reports what the exact fingerprint table
/// knows. The top list names no sequences: that requires one of the --overrep-mode table
/// modes, and `names_available` says so.
pub fn repeated_metrics(acc: &Accum, min_count: u32) -> (Map<String, Value>, Value) {
    let reads = acc.reads.max(1) as f64;
    let distinct = acc.seq_counts.distinct();
    let hist = acc.seq_counts.multiplicity_hist();
    // The threshold is the caller's: `--overrep-min-count` when given, otherwise the
    // resolved 0.1%-of-reads floor. Computing a second threshold here made this metric
    // say "1 repeated sequence" while the `overrep` block it summarises said [].
    let overrep = acc.seq_counts.at_least(min_count);

    let mut m = Map::new();
    m.insert("distinct_seqs".into(), json!(distinct));
    m.insert(
        "duplication_level".into(),
        json!(r4(1.0 - distinct as f64 / reads)),
    );
    m.insert(
        "total_deduplicated_pct".into(),
        json!(r4(100.0 * (reads - distinct as f64) / reads)),
    );
    m.insert("dup_multiplicity_max".into(), json!(hist.len().saturating_sub(2)));
    m.insert("dup_singletons".into(), json!(hist.get(1).copied().unwrap_or(0)));
    m.insert("repeated_seq_count".into(), json!(overrep.len()));
    // Integer sum first: `f64::sum` over an empty iterator is -0.0 (the additive
    // identity), which serialised as "repeated_seq_frac": -0.0 in the report.
    let repeated_reads: u64 = overrep.iter().map(|(_, c)| *c as u64).sum();
    m.insert("repeated_seq_frac".into(), json!(r4(repeated_reads as f64 / reads)));

    // `names_available: false` is not a caveat, it is the contract of this mode: the
    // multiplicities are exact, the sequences are not recoverable from a hash.
    let list: Vec<Value> = overrep
        .iter()
        .take(20)
        .map(|(fp, c)| {
            json!({
                "fingerprint": format!("{fp:016x}"),
                "count": c,
                "frac": r4(*c as f64 / reads),
            })
        })
        .collect();
    (m, json!({ "names_available": false, "top": list }))
}

/// 12-mer spectrum: metrics + the top list.
///
/// Exact: the counter is a dense 4^12 `u16` table (32 MiB, allocated once, independent
/// of the number of reads) with an exact overflow map for k-mers that occur more than
/// 65 535 times. There is no cap and no eviction, so `approximate` is `false`, `distinct`
/// is a count of non-zero cells rather than an estimate, and the block is the same
/// whatever `-j` is.
pub fn kmer_metrics(acc: &Accum) -> (Map<String, Value>, Value) {
    let reads = acc.reads.max(1) as f64;
    // The denominator is the number of k-mer WINDOWS counted, not the number of reads: a
    // read of L bases carries L-11 windows, so a 12-mer repeating inside every read can be
    // counted L-11 times per read. Dividing by reads made the share unbounded (2 000 reads
    // of 40 bp: 16 000 windows, kmer_top_frac 8.0).
    let windows = acc.kmer_windows.max(1) as f64;
    let min_count = (reads * 0.001).max(10.0) as u32;
    let kmers = acc.kmer_counts.report(min_count);
    let top = kmers.first().map(|(_, c)| *c).unwrap_or(0);
    let mut m = Map::new();
    m.insert("kmer_top_count".into(), json!(top));
    m.insert(
        "kmer_top_frac".into(),
        json!(r4((top as f64 / windows).min(1.0))),
    );
    let list: Vec<Value> = kmers
        .iter()
        .take(20)
        .map(|(k, c)| {
            json!({
                "kmer": repeated::unpack_kmer(*k, 12),
                "count": c,
                "frac": r4((*c as f64 / windows).min(1.0)),
            })
        })
        .collect();
    (
        m,
        json!({
            "approximate": false,
            "exactness": "dense 4^12 = 16 777 216-cell u16 table, one allocation at scan \
                          start, merged cell by cell; counts above 65 535 are carried exactly \
                          in an overflow map, so only the cell WIDTH is bounded and the count \
                          is not",
            "distinct": acc.kmer_counts.distinct(),
            "top": list,
        }),
    )
}

/// Heatmap: metrics (peak + its cell) + the grid block.
pub fn heatmap_metrics(acc: &Accum) -> (Map<String, Value>, Value) {
    let (mut peak, mut px, mut py) = (0u64, 0usize, 0usize);
    for y in 0..heatmap::SIZE {
        for x in 0..heatmap::SIZE {
            let c = acc.heatmap[y * heatmap::SIZE + x];
            if c > peak {
                peak = c;
                px = x;
                py = y;
            }
        }
    }
    let mut m = Map::new();
    m.insert("heatmap_peak".into(), json!(peak));
    m.insert("heatmap_peak_x".into(), json!(px));
    m.insert("heatmap_peak_y".into(), json!(py));
    (m, json!({ "size": heatmap::SIZE, "counts": acc.heatmap }))
}
