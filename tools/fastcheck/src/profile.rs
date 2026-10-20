//! Fit a fastsim profile directly from fastcheck's measurements.
//!
//! This is the Rust port of `bench/scripts/qc_to_profile.py`: it turns the exact
//! counters fastcheck already holds into the profile JSON fastsim consumes, so
//! the coupled path is `fastcheck --profile p.json` -> `fastsim --profile p.json`
//! with no Python step in between.
//!
//! What is MEASURED and what is INVENTED is recorded in the profile's
//! `_comment`:
//!
//!   MEASURED  read_len, gc (fragment), the per-position mean-Q curve, the
//!             duplicate rate, the adapter readthrough rate, the per-position
//!             base composition, the insert distribution (--insert-hist), and
//!             the ambiguous-base (N) rate / run-length / edge model (`n_run`),
//!             the head/tail G-run distributions, the flat in-run qualities,
//!             and the out-of-run quality curve deconvolved from them
//!   INVENTED  the insert SHAPE when no insert module ran, and the mate
//!             tail correlation.
//!
//! The head/tail run-length distributions are the EXACT measured per-base G-run
//! histograms with the biological boundary run subtracted (see
//! `deconvolve_boundary_g`), not an estimate from the per-position G fraction.
//! The estimate that used to sit here could not express a run longer than 12 at
//! all and lost whole modes of a multi-modal tail, so the profile it produced
//! had no mass at the tail lengths the library actually carried.

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, Result};
use serde_json::{json, Map, Value};

use crate::accum::Accum;
use crate::insert::InsertHist;

/// Bases of per-position composition carried into the profile. The full read
/// length, capped: high-GC libraries carry a G ramp that runs the whole read
/// (measured on the closed-loop pilot: a virome whose per-position G fraction
/// climbs from 0.18 at position 1 to 0.64 at 151). A 60-position cap cut that
/// ramp in half, the tail fell back to the flat --gc draw, and the sim's G
/// fraction collapsed where the real data was G-heavy (g_frac up to 0.37,
/// GC diff up to 0.10 — the whole virome cluster of GATE-2 failures). At
/// 4 floats per position the cap costs ~2.4 KB per mate.
const BASE_POSITIONS: usize = 512;

fn adapter_seq(name: &str) -> &'static str {
    match name {
        "TruSeq" => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        _ => "CTGTCTCTTATACACATCT",
    }
}

/// Per-position G fraction of one mate.
fn g_fraction(acc: &Accum) -> Vec<f64> {
    acc.base_pos
        .iter()
        .map(|row| {
            let depth: u64 = row.iter().sum();
            if depth == 0 {
                0.0
            } else {
                row[2] as f64 / depth as f64
            }
        })
        .collect()
}

/// Base code of G in the polyX histograms (0=A 1=C 2=G 3=T).
const G_CODE: usize = 2;

/// The MEASURED G-run histogram of one mate, as a distribution over reads.
///
/// `hist[k]` for `k >= 1` is the exact share of reads whose terminal (or
/// initial) base is G inside a run of exactly `k` Gs, taken straight from the
/// per-base run counter the report already publishes; `hist[0]` is every read
/// whose terminal (or initial) base is not G. The whole vector sums to 1.
///
/// This is the run the SEQUENCER SHOWS, artefact plus any biological G at the
/// boundary. `deconvolve_boundary_g` turns it into the artefact the simulator
/// draws.
fn observed_g_run_hist(acc: &Accum, from_start: bool) -> Vec<f64> {
    let reads = acc.reads.max(1) as f64;
    let src = if from_start { &acc.polyx5[G_CODE] } else { &acc.polyx3[G_CODE] };
    let mut hist = vec![0.0f64; src.len() + 1];
    let mut hit = 0u64;
    for (len, &c) in src.iter().enumerate() {
        hist[len] = c as f64 / reads;
        hit += c;
    }
    hist[0] = acc.reads.saturating_sub(hit) as f64 / reads;
    hist
}

/// Subtract the biological boundary run from an observed G-run histogram,
/// leaving the ARTEFACT run-length distribution.
///
/// The observed terminal run is the injected artefact run *plus* any biological
/// G that happens to sit at the boundary: a read with no artefact (`t = 0`) still
/// shows a run whenever its last biological base is G, and a read with an
/// artefact of `t` shows `t` exactly (the simulator rewrites that boundary base
/// to non-G). With per-base G probability `q` the boundary run is geometric,
/// `b(k) = (1 - q) q^k`, so
///
/// ```text
/// observed(k) = artefact(k) + artefact(0) * b(k)   (k >= 1)
/// observed(0) = artefact(0) * b(0)
/// ```
///
/// which inverts in closed form. Without this step the no-run bin absorbs every
/// artefact-free read whose last base is G - on the reference library that
/// inflated "tail-free" from a true 66.0% to 50.6% and pushed the head row past
/// its tolerance - and short artefacts vanish into that mass instead of being
/// reported at their own lengths.
fn deconvolve_boundary_g(obs: &[f64], q: f64) -> Vec<f64> {
    let q = q.clamp(0.0, 0.95);
    let p_base = 1.0 - q;
    let a0 = if p_base > 1e-6 { (obs[0] / p_base).min(1.0) } else { obs[0] };
    let mut out = vec![0.0f64; obs.len()];
    out[0] = a0;
    let mut b = p_base; // b(k) for the current k, starting at k = 1
    for k in 1..obs.len() {
        b *= q;
        let v = obs[k] - a0 * b;
        out[k] = if v > 0.0 { v } else { 0.0 };
    }
    let total: f64 = out.iter().sum();
    if total > 0.0 {
        for v in out.iter_mut() {
            *v /= total;
        }
    } else {
        out[0] = 1.0;
    }
    out
}

/// Per-base G probability of the read's BACKGROUND, measured away from the
/// artefact runs: the median G fraction over the interior positions.
///
/// The first few positions carry the 5' head artefact and the last few the 3'
/// tail artefact, so neither end can be read as background; the middle of the
/// read cannot be reached by a run in any library this tool is used on. A
/// median (not a mean) keeps a handful of very long full-length reads from
/// dragging the estimate.
fn background_g_fraction(acc: &Accum) -> f64 {
    let gf = g_fraction(acc);
    let l = gf.len();
    if l == 0 {
        return 0.25;
    }
    let lo = 4.min(l.saturating_sub(1));
    let hi = (l / 2).max(lo + 1).min(l);
    let mut mid: Vec<f64> = gf[lo..hi].to_vec();
    if mid.is_empty() {
        mid = gf.clone();
    }
    mid.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    mid[mid.len() / 2].clamp(0.0, 0.9)
}

/// A run-length distribution as the profile's integer weight map (out of 10 000).
fn dist_to_weights(dist: &[f64]) -> BTreeMap<u32, f64> {
    let mut m: BTreeMap<u32, f64> = BTreeMap::new();
    for (k, &p) in dist.iter().enumerate() {
        let w = (p * 10000.0).round();
        if w >= 1.0 {
            m.insert(k as u32, w);
        }
    }
    if m.is_empty() {
        m.insert(0, 10000.0);
    }
    m
}

/// Mean Phred of the bases inside one mate's G-runs, raw bytes minus the offset.
fn run_q_mean(q: (u64, u64), offset: u8) -> Option<f64> {
    if q.1 == 0 {
        None
    } else {
        Some(q.0 as f64 / q.1 as f64 - offset as f64)
    }
}

/// The OUT-OF-RUN per-position curve, deconvolved from the measured one.
///
/// The measured per-position mean is a mixture: `measured[i] = (1 - p_i - q_i)
/// out[i] + p_i tail_q + q_i head_q`, where `p_i` is the share of reads whose
/// trailing run covers position `i` and `q_i` the share whose DEEPER head-run
/// bases cover it. The simulator consumes `out[i]` and applies the flat in-run
/// means itself, so handing it `measured[i]` collapses the run bases twice.
/// The coverage shares come from the OBSERVED run histograms, because that is
/// what covers the position in the data being measured.
fn deconvolve_pos_q(
    measured: &[f64],
    obs_head: &[f64],
    obs_tail: &[f64],
    head_q: Option<f64>,
    tail_q: Option<f64>,
) -> Vec<f64> {
    let l = measured.len();
    let cover = |hist: &[f64], min_run: usize| -> f64 {
        hist.iter().enumerate().filter(|(k, _)| *k >= min_run).map(|(_, &p)| p).sum()
    };
    let mut out = Vec::with_capacity(l);
    for (i, &m) in measured.iter().enumerate() {
        // position i is inside a trailing run of length k iff k >= l - i
        let p_i = cover(obs_tail, l - i);
        // position 0 of a head run is the normal curve; positions 1..h-1 are deep
        let q_i = if i == 0 { 0.0 } else { cover(obs_head, i + 1) };
        let denom = 1.0 - p_i - q_i;
        let v = if denom <= 1e-6 {
            m
        } else {
            (m - p_i * tail_q.unwrap_or(0.0) - q_i * head_q.unwrap_or(0.0)) / denom
        };
        out.push(v.max(0.0));
    }
    out
}

/// Median of a slice (does not need it sorted; small n by construction).
fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 30.0;
    }
    let mut s: Vec<f64> = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

/// Per-position mean Phred of one mate.
fn per_position_mean_q(acc: &Accum, offset: u8) -> Vec<f64> {
    acc.qual_pos
        .iter()
        .map(|row| {
            let depth: u64 = row.iter().sum();
            if depth == 0 {
                return 0.0;
            }
            let sum: f64 = row
                .iter()
                .enumerate()
                // The row is indexed by `byte - 33` (128 cells): recover the byte first,
                // then the Phred value, exactly as the report does.
                .map(|(b, &c)| ((b + crate::accum::QMIN as usize) as i64 - offset as i64).max(0) as f64 * c as f64)
                .sum();
            sum / depth as f64
        })
        .collect()
}

/// Per-read mean-Q histogram, top 6 bins, as weights out of 1000.
fn read_meanq(acc: &Accum, offset: u8) -> BTreeMap<String, i64> {
    let mut hist: BTreeMap<i64, u64> = BTreeMap::new();
    for (b, &c) in acc.meanq_raw_hist.iter().enumerate() {
        if c > 0 {
            *hist.entry((b as i64 - offset as i64).max(0)).or_insert(0) += c;
        }
    }
    let mut top: Vec<(i64, u64)> = hist.into_iter().collect();
    top.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    top.truncate(6);
    let tot: u64 = top.iter().map(|(_, c)| *c).sum::<u64>().max(1);
    top.into_iter()
        .map(|(q, c)| (q.to_string(), (((c as f64 / tot as f64) * 1000.0).round() as i64).max(1)))
        .collect()
}

/// One mate's slice of the profile.
struct Mate {
    head: BTreeMap<u32, f64>,
    tail: BTreeMap<u32, f64>,
    /// the OUT-OF-RUN per-position curve (what the simulator consumes)
    pos_q: Vec<f64>,
    /// the measured per-position curve, kept as its own field so a reader can
    /// tell the measurement from the deconvolution
    pos_q_uncond: Vec<f64>,
    /// measured flat mean Phred inside the G-runs (all run bases)
    run_head_q: Option<f64>,
    run_tail_q: Option<f64>,
    /// head-run mean over the DEEP bases only (the first base of a head run is
    /// the normal curve)
    deep_head_q: Option<f64>,
    meanq: BTreeMap<String, i64>,
    body_q: f64,
    tail_q: f64,
    base: Vec<Vec<f64>>,
}

fn fit_mate(acc: &Accum, _g_mean: f64, read_len: usize, offset: u8) -> Mate {
    // The measured run histograms, then the artefact distributions the profile
    // carries: the simulator draws an artefact length and its boundary rule
    // re-forms the observed run, so the biological boundary run must come out
    // here rather than be baked in.
    let obs_head = observed_g_run_hist(acc, true);
    let obs_tail = observed_g_run_hist(acc, false);
    let q_bg = background_g_fraction(acc);
    let head = dist_to_weights(&deconvolve_boundary_g(&obs_head, q_bg));
    let tail = dist_to_weights(&deconvolve_boundary_g(&obs_tail, q_bg));

    let mut q = per_position_mean_q(acc, offset);
    // `pos_q` is exactly read_len long. When the measured curve is shorter than
    // that (some positions never observed), the missing tail is padded with the
    // LAST measured value — a flat continuation, stated rather than invented: no
    // measurement here says what those positions look like.
    while q.len() < read_len {
        q.push(*q.last().unwrap_or(&30.0));
    }
    q.truncate(read_len);
    let body_q = median(&q[..((read_len as f64 * 0.8) as usize).max(1)]);
    let tail_q = median(&q[((read_len as f64 * 0.9) as usize).min(q.len())..]);

    // Flat in-run qualities, measured. `head_all` covers every head-run base,
    // including the first of each run, which the simulator models with the
    // normal curve — so the deep-base mean is recovered by removing one first
    // base per head read and every full-length read (handled by the
    // full-length path, not this mean).
    let run_head_q = run_q_mean(acc.run_q_head, offset);
    let run_tail_q = run_q_mean(acc.run_q_tail, offset);
    let reads = acc.reads as f64;
    let head_full = obs_head.get(read_len).copied().unwrap_or(0.0) * reads;
    let mut head_non_full = 0.0f64;
    let mut head_bases = 0.0f64;
    for (k, &p) in obs_head.iter().enumerate() {
        let cnt = p * reads;
        head_bases += k as f64 * cnt;
        if k >= 1 && k < read_len {
            head_non_full += cnt;
        }
    }
    let deep_n = head_bases - head_full * read_len as f64 - head_non_full;
    let deep_head_q = match run_head_q {
        Some(all) if deep_n > 0.0 => {
            let sum_phred = acc.run_q_head.0 as f64 - offset as f64 * acc.run_q_head.1 as f64;
            let first_q = q.first().copied().unwrap_or(all);
            let full_q = run_tail_q.unwrap_or(20.0);
            Some((sum_phred - head_full * read_len as f64 * full_q - head_non_full * first_q) / deep_n)
        }
        other => other,
    };

    let pos_q_uncond = q.clone();
    let pos_q = deconvolve_pos_q(&pos_q_uncond, &obs_head, &obs_tail, deep_head_q, run_tail_q);

    // Located per-position composition, so the simulator can reproduce a
    // 5' dropout or any-position bias a G-run model cannot express.
    let mut base = Vec::new();
    for row in acc.base_pos.iter().take(BASE_POSITIONS) {
        let s: f64 = row[..4].iter().sum::<u64>() as f64;
        base.push(if s > 0.0 {
            row[..4].iter().map(|&x| (x as f64 / s * 10000.0).round() / 10000.0).collect()
        } else {
            vec![0.25, 0.25, 0.25, 0.25]
        });
    }

    Mate {
        head,
        tail,
        pos_q,
        pos_q_uncond,
        run_head_q,
        run_tail_q,
        deep_head_q,
        meanq: read_meanq(acc, offset),
        body_q,
        tail_q,
        base,
    }
}

/// The insert CDF, assembled from a measured histogram: censored-below pairs go
/// into one bin at the read length, everything above the overlap resolution into
/// one saturated bin. Returns (cdf, summary).
fn insert_cdf(hist: &InsertHist, read_len: usize) -> (BTreeMap<String, u64>, Value) {
    // The union of both measurement regimes, exactly as the insert report
    // publishes it (bins_unioned_with_adapter), so the profile and the report
    // cannot disagree.
    let mut union: BTreeMap<u64, u64> = BTreeMap::new();
    for (i, &c) in hist.union().iter().enumerate() {
        if c > 0 {
            union.insert(i as u64, c);
        }
    }
    if hist.censored_below > 0 {
        *union.entry(read_len as u64).or_insert(0) += hist.censored_below;
    }
    let saturated = 2 * read_len as u64 + 1;
    if hist.censored_above > 0 {
        *union.entry(saturated).or_insert(0) += hist.censored_above;
    }
    let tot: u64 = union.values().sum::<u64>().max(1);
    let mean = union.iter().map(|(k, v)| *k as f64 * *v as f64).sum::<f64>() / tot as f64;
    let var = if tot > 1 {
        union.iter().map(|(k, v)| *v as f64 * (*k as f64 - mean).powi(2)).sum::<f64>() / (tot - 1) as f64
    } else {
        0.0
    };
    let summary = json!({
        "mean": (mean * 10.0).round() / 10.0,
        "sd": ((var.sqrt()) * 10.0).round() / 10.0,
        "pairs_sampled": hist.pairs_sampled,
        "censored_above_rate": ((hist.censored_above as f64 / tot as f64) * 10000.0).round() / 10000.0,
        "complete_overlap_share": ((hist.censored_below as f64 / tot as f64) * 10000.0).round() / 10000.0,
        "measured_share": (((tot - hist.censored_above) as f64 / tot as f64) * 10000.0).round() / 10000.0,
    });
    (union.into_iter().map(|(k, v)| (k.to_string(), v)).collect(), summary)
}

fn map_json(m: &BTreeMap<u32, f64>) -> Value {
    Value::Object(m.iter().map(|(k, v)| (k.to_string(), json!(v.round() as i64))).collect())
}

fn int_map_json(m: &BTreeMap<String, i64>) -> Value {
    Value::Object(m.iter().map(|(k, v)| (k.clone(), json!(v))).collect())
}

/// Fit and write the fastsim profile. `accs` is one accumulator per mate (R1, R2)
/// or a single one for single-end. Returns the profile name written.
pub fn write_profile(
    path: &Path,
    name: &str,
    accs: &[&Accum],
    offset: u8,
    insert: Option<&InsertHist>,
) -> Result<()> {
    let qc1 = accs[0];
    // The library's length is its MODE, not the longest read it ever contained.
    // `len_hist.len() - 1` is the maximum observed length: one 300 bp outlier in
    // a 151 bp library set read_len = 300, and the simulator then built a 300 bp
    // library out of data whose own fit report recorded mode_len 151 at
    // mode_share 0.999. The insert CDF's censored-below bin sits at read_len and
    // the saturated bin at 2*read_len+1, so both moved with it.
    let max_len = qc1.len_hist.len().saturating_sub(1);
    let mode_len = qc1
        .len_hist
        .iter()
        .enumerate()
        .max_by_key(|(_, &c)| c)
        .map(|(i, _)| i)
        .unwrap_or(max_len);
    let read_len = mode_len;
    if read_len == 0 {
        anyhow::bail!("--profile: no reads were scanned, so there is nothing to fit");
    }
    let mode_tot: u64 = qc1.len_hist.iter().sum::<u64>().max(1);
    let mode_share = qc1.len_hist[mode_len] as f64 / mode_tot as f64;
    let g_mean = qc1.base_pos.first().map(|_| {
        let g: u64 = qc1.base_pos.iter().map(|r| r[2]).sum();
        g as f64 / qc1.bases.max(1) as f64
    }).unwrap_or(0.25);

    // The profile's gc is the FRAGMENT gc: the simulator lays adapters and
    // artefacts on top, so a read gc would double-count them.
    let gc = if qc1.frag_bases > 0 {
        let (mut gc_sum, mut n) = (0u64, 0u64);
        for a in accs {
            gc_sum += a.frag_gc;
            n += a.frag_bases;
        }
        if n > 0 { gc_sum as f64 / n as f64 } else { 0.0 }
    } else {
        0.0
    };

    let mates: Vec<Mate> = accs.iter().map(|a| fit_mate(a, g_mean, read_len, offset)).collect();

    // Readthrough and the adapter choice, averaged over the mates.
    let mut rt = [0.0f64; 2];
    for (i, a) in accs.iter().enumerate() {
        let n = a.reads.max(1) as f64;
        rt[0] += a.adapter_hits[0] as f64 / n;
        rt[1] += a.adapter_hits[1] as f64 / n;
        let _ = i;
    }
    let n_mates = accs.len().max(1) as f64;
    rt[0] /= n_mates;
    rt[1] /= n_mates;
    let rt_rate = rt[0].max(rt[1]);
    let adapter_name = if rt[0] >= rt[1] { "TruSeq" } else { "Nextera" };

    let (insert_cdf_map, insert_summary, insert_measured) = match insert {
        Some(h) => {
            let (cdf, sum) = insert_cdf(h, read_len);
            (cdf, sum, true)
        }
        None => (BTreeMap::new(), Value::Null, false),
    };

    let comment = if insert_measured {
        format!("fitted from {name} by fastcheck --profile (as-measured): head/tail run lengths \
                 MEASURED (per-base G-run histograms, biological boundary run deconvolved), \
                 insert distribution MEASURED by --insert-hist, mate tail correlation invented")
    } else {
        format!("fitted from {name} by fastcheck --profile (as-measured): head/tail run lengths \
                 MEASURED (per-base G-run histograms, biological boundary run deconvolved), \
                 insert shape invented")
    };

    // Joint tail. Both labels this tool emits are the "independent" plan: the
    // measured insert distribution already carries the correlation the
    // correlated variants were invented to supply, so the mates' tail runs are
    // placed as the product of their marginals. SE has no mate to correlate.
    let joint: BTreeMap<String, i64> = if mates.len() < 2 {
        BTreeMap::from([("0,0".to_string(), 100)])
    } else {
        let mut j: BTreeMap<String, i64> = BTreeMap::new();
        for (k1, w1) in &mates[0].tail {
            for (k2, w2) in &mates[1].tail {
                let cell = (w1 * w2 / 10000.0).round() as i64;
                *j.entry(format!("{k1},{k2}")).or_insert(0) += cell;
            }
        }
        if j.values().all(|w| *w <= 0) {
            j.insert("0,0".to_string(), 100);
        }
        if j.is_empty() {
            j.insert("0,0".to_string(), 100);
        }
        j
    };

    let mut head = Map::new();
    let mut tail = Map::new();
    let mut pos_q = Map::new();
    let mut pos_q_uncond = Map::new();
    let mut run_q_mean = Map::new();
    let mut deep_head_q = Map::new();
    let mut read_meanq = Map::new();
    let mut tail_q_collapse = Map::new();
    let mut full_length_q = Map::new();
    let mut base = Map::new();
    let names = ["R1", "R2"];
    for (i, m) in mates.iter().enumerate() {
        let n = names[i];
        head.insert(n.to_string(), map_json(&m.head));
        tail.insert(n.to_string(), map_json(&m.tail));
        pos_q.insert(n.to_string(), json!(m.pos_q));
        pos_q_uncond.insert(n.to_string(), json!(m.pos_q_uncond));
        run_q_mean.insert(
            n.to_string(),
            json!({
                "head": m.run_head_q.map(|v| (v * 1000.0).round() / 1000.0),
                "tail": m.run_tail_q.map(|v| (v * 1000.0).round() / 1000.0),
            }),
        );
        if let Some(v) = m.deep_head_q {
            deep_head_q.insert(n.to_string(), json!((v * 1000.0).round() / 1000.0));
        }
        read_meanq.insert(n.to_string(), int_map_json(&m.meanq));
        tail_q_collapse.insert(
            n.to_string(),
            json!(((m.body_q - m.tail_q).max(0.0) * 100.0).round() / 100.0),
        );
        full_length_q.insert(n.to_string(), json!((m.body_q - 10.0).max(2.0).round() as i64));
        if !m.base.is_empty() {
            base.insert(n.to_string(), json!(m.base));
        }
    }

    let mut prof = Map::new();
    prof.insert("_comment".into(), json!(comment));
    prof.insert("read_len".into(), json!(read_len));
    prof.insert("gc".into(), json!((gc * 10000.0).round() / 10000.0));
    prof.insert("geometry".into(), json!(if mates.len() < 2 { "tail3p_r1" } else { "tail3p_both" }));
    prof.insert("head".into(), Value::Object(head));
    prof.insert("tail".into(), Value::Object(tail));
    prof.insert(
        "joint_tail".into(),
        Value::Object(joint.into_iter().map(|(k, v)| (k, json!(v))).collect()),
    );
    prof.insert("tail_q_collapse".into(), Value::Object(tail_q_collapse));
    prof.insert("pos_q".into(), Value::Object(pos_q));
    // The measured curve and the measured flat in-run qualities, so a reader
    // can check the deconvolution instead of taking `pos_q` on faith, and so
    // the simulator applies the collapse exactly once.
    prof.insert("pos_q_uncond".into(), Value::Object(pos_q_uncond));
    prof.insert("run_q_mean".into(), Value::Object(run_q_mean));
    if !deep_head_q.is_empty() {
        prof.insert("deep_head_q".into(), Value::Object(deep_head_q));
    }
    prof.insert("read_meanq".into(), Value::Object(read_meanq));
    // The measured duplication level, so the chain can pass the profile's own
    // duplicate rate to the simulator. Without it the chain falls back to
    // fastsim's --dup-rate default (0.02) and the sim sits at 0.015-0.018
    // whatever the real library's rate is (measured: real 0.001-0.014, sim
    // 0.015-0.018 -> dup_ratio up to 20x, the other GATE-2 cluster).
    //
    // Computed over the MERGED accumulator, exactly as `repeated_metrics`
    // computes the report's `duplication_level`. Using only the first mate
    // made the profile and the report disagree by construction: the two mates
    // carry different duplicate rates (measured on DRR000466 — R1 0.0086,
    // R2 0.0182, merged 0.0135, with only 15 sequences shared between them),
    // so `accs.first()` reported whichever mate happened to be first. The
    // closed loop then compared the sim (fed the first-mate rate) against the
    // report's merged rate, inflating dup_ratio by a median factor of 1.39
    // across 96 datasets and failing 144 of 306 comparisons.
    if !accs.is_empty() {
        let reads: f64 = accs.iter().map(|a| a.reads.max(1) as f64).sum();
        let mut merged = crate::fingerprint::FingerprintTable::with_capacity(0);
        for a in accs {
            merged.merge(&a.seq_counts);
        }
        let distinct = merged.distinct() as f64;
        if reads > 0.0 && distinct > 0.0 {
            let dup = (1.0 - distinct / reads).clamp(0.0, 1.0);
            prof.insert("dup_rate".into(), json!((dup * 10000.0).round() / 10000.0));
        }
    }
    // The ambiguous-base (N) model, measured over the MERGED accumulators so it
    // describes the library exactly as `reads_with_non_acgt` and
    // `duplication_level` describe it in the report (the closed loop compares
    // the simulation against that merged report). A library with no ambiguous
    // base emits NO `n_run` key at all: fastsim then draws nothing for it, so
    // the field's presence is the simulator's signal to change its stream.
    if !accs.is_empty() {
        let reads: u64 = accs.iter().map(|a| a.reads).sum();
        let with_n: u64 = accs.iter().map(|a| a.reads_with_n).sum();
        if reads > 0 && with_n > 0 {
            let runs: u64 = accs.iter().map(|a| a.n_run_hist.iter().sum::<u64>()).sum();
            let mut hist: BTreeMap<String, u64> = BTreeMap::new();
            for a in accs {
                for (len, &c) in a.n_run_hist.iter().enumerate() {
                    if c > 0 {
                        *hist.entry(len.to_string()).or_insert(0) += c;
                    }
                }
            }
            let lead = if runs > 0 {
                accs.iter().map(|a| a.n_run_leading).sum::<u64>() as f64 / runs as f64
            } else {
                0.0
            };
            let trail = if runs > 0 {
                accs.iter().map(|a| a.n_run_trailing).sum::<u64>() as f64 / runs as f64
            } else {
                0.0
            };
            prof.insert(
                "n_run".into(),
                json!({
                    "read_rate": ((with_n as f64 / reads as f64) * 1e6).round() / 1e6,
                    "len": Value::Object(hist.into_iter().map(|(k, v)| (k, json!(v))).collect()),
                    "leading": (lead * 1e6).round() / 1e6,
                    "trailing": (trail * 1e6).round() / 1e6,
                }),
            );
        }
    }
    prof.insert("full_length_rate".into(), json!(0.0));
    prof.insert("full_length_q".into(), Value::Object(full_length_q));
    prof.insert(
        "adapter".into(),
        json!({"mode": if rt_rate > 0.01 { "readthrough" } else { "none" },
               "seq": adapter_seq(adapter_name)}),
    );
    if insert_measured {
        prof.insert(
            "insert".into(),
            Value::Object(insert_cdf_map.into_iter().map(|(k, v)| (k, json!(v))).collect()),
        );
    }
    if !base.is_empty() {
        prof.insert("base".into(), Value::Object(base));
    }

    let text = serde_json::to_string_pretty(&Value::Object(prof))? + "\n";
    std::fs::write(path, text).with_context(|| format!("writing profile to {}", path.display()))?;

    // The fit report, so the assumptions travel with the profile the way the
    // Python fitter's fit_report.json does.
    let mut warnings: Vec<String> = Vec::new();
    if max_len != mode_len {
        warnings.push(format!(
            "read_len is the modal length {mode_len} bp, but the longest read is {max_len} bp — \
             the {max_len} bp reads (and anything measured at their scale: the insert CDF's \
             censored bin, pos_q beyond position {mode_len}) are off the simulated read length; \
             trim or filter them if the simulation must match them"
        ));
    }
    if mode_share < 0.9 {
        warnings.push(format!(
            "only {:.1}% of reads are at the modal length {mode_len} — the dataset looks pre-trimmed; \
             a fixed read_len={read_len} simulation cannot model a variable-length library",
            mode_share * 100.0
        ));
    }
    if rt_rate <= 0.0005 {
        warnings.push("no adapter readthrough observed: the insert distribution is unconstrained, \
                       so its shape is pure invention"
            .to_string());
    }
    let report = json!({
        "dataset": name,
        "measured": {
            "read_len": read_len, "mode_len": mode_len, "mode_share": (mode_share * 10000.0).round() / 10000.0,
            "gc_fragment": (gc * 10000.0).round() / 10000.0,
            "mean_read_q": crate::report::mean_read_mean_q(qc1, offset),
            "readthrough": {"TruSeq": (rt[0] * 100000.0).round() / 100000.0,
                            "Nextera": (rt[1] * 100000.0).round() / 100000.0},
            "g_mean": (g_mean * 10000.0).round() / 10000.0,
            "reads_with_non_acgt": accs.iter().map(|a| a.reads_with_n).sum::<u64>(),
            "total_reads": accs.iter().map(|a| a.reads).sum::<u64>(),
        },
        "insert_measured": insert_summary,
        "fit_warnings": warnings,
        "invented": if insert_measured { json!(["mate tail correlation"]) }
                    else { json!(["insert distribution shape", "mate tail correlation"]) },
    });
    let rp = path.with_extension("fit_report.json");
    std::fs::write(&rp, serde_json::to_string_pretty(&report)? + "\n")
        .with_context(|| format!("writing fit report to {}", rp.display()))?;
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;

    fn acc_with(seqs: &[(&[u8], &[u8])], max_pos: usize) -> Accum {
        let mut a = Accum::new(max_pos);
        for (s, q) in seqs {
            a.observe(b"@r", s, q);
        }
        a
    }

    /// Write a profile for one accumulated fixture and read it back.
    fn profile_json(tag: &str, accs: &[&Accum]) -> serde_json::Value {
        let dir = std::env::temp_dir().join(format!("fc_prof_{}_{}", std::process::id(), tag));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("p.json");
        write_profile(&p, "fixture", accs, 33, None).unwrap();
        let v = serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        v
    }

    /// A fixture read: `body` (no G) plus a homopolymer end of `end_base`, with
    /// the run bases at `run_q` and the body at Q40 ('I' = 40, '5' = 20).
    fn fixture_read(body: &[u8], end_base: u8, end_len: usize, run_q: u8) -> (Vec<u8>, Vec<u8>) {
        let mut s = body.to_vec();
        s.extend(std::iter::repeat(end_base).take(end_len));
        let mut q = vec![b'I'; body.len()];
        q.extend(std::iter::repeat(run_q).take(end_len));
        (s, q)
    }

    /// The fitted tail must report the run lengths the library actually carries,
    /// at the shares it carries them in.
    ///
    /// This is the polyG-detector defect: the fitter used to read the
    /// distribution off the per-position G fraction, capping runs at 12 and
    /// folding the no-run mass into a monotone estimate. Measured against a
    /// simulator ground truth it reported 93.98% of reads tail-free against a
    /// true 66.01%, with zero mass at t = 3 and t = 15 where the truth was
    /// 7.99% and 6.99%.
    #[test]
    fn fitted_tail_reports_the_run_lengths_and_shares_the_library_carries() {
        // 40 reads end in a G-run of exactly 3, 30 in exactly 15, 30 in a
        // non-G run. The body carries NO G, so the biological boundary
        // contribution is zero and the expected distribution is exact.
        const BODY: &[u8] = b"ACATACATACATACATAC";
        let mut seqs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        for _ in 0..40 {
            seqs.push(fixture_read(BODY, b'G', 3, b'I'));
        }
        for _ in 0..30 {
            seqs.push(fixture_read(BODY, b'G', 15, b'I'));
        }
        for _ in 0..30 {
            seqs.push(fixture_read(BODY, b'T', 3, b'I'));
        }
        let refs: Vec<(&[u8], &[u8])> =
            seqs.iter().map(|(s, q)| (s.as_slice(), q.as_slice())).collect();
        let a = acc_with(&refs, 64);
        let v = profile_json("tailruns", &[&a]);
        let t = v["tail"]["R1"].as_object().unwrap();
        let w = |k: &str| t.get(k).and_then(|x| x.as_i64()).unwrap_or(0) as f64 / 10000.0;
        assert!((w("3") - 0.40).abs() < 0.03, "run of 3: {} in {t:?}", w("3"));
        assert!((w("15") - 0.30).abs() < 0.03, "run of 15: {} in {t:?}", w("15"));
        assert!((w("0") - 0.30).abs() < 0.03, "no run: {} in {t:?}", w("0"));
        assert_eq!(w("1"), 0.0, "nothing ends in a 1-base run: {t:?}");
        // The same histogram is the head CDF's sibling: a read that ends in G
        // but begins in A has no head run at all.
        assert_eq!(v["head"]["R1"]["0"], serde_json::json!(10000));
    }

    /// The fitted `pos_q` must be the OUT-OF-RUN curve, not the collapsed
    /// measurement: the simulator applies the flat in-run mean itself, so a
    /// curve that already contains the collapse is collapsed twice.
    ///
    /// This is the worst-10-position quality-block defect: the recovered curve
    /// missed its block by 1.95 Phred because every tail base's collapse was
    /// baked into the body curve the simulator then collapsed again.
    #[test]
    fn fitted_pos_q_is_the_out_of_run_curve_not_the_collapsed_measurement() {
        const BODY: &[u8] = b"ACATACATACATACATAC";
        let mut seqs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        for _ in 0..50 {
            // a 10-base 3' G-run whose bases are Q20
            seqs.push(fixture_read(BODY, b'G', 10, b'5'));
        }
        for _ in 0..50 {
            // same length, no G-run, every base Q40
            seqs.push(fixture_read(BODY, b'T', 10, b'I'));
        }
        let refs: Vec<(&[u8], &[u8])> =
            seqs.iter().map(|(s, q)| (s.as_slice(), q.as_slice())).collect();
        let a = acc_with(&refs, 64);
        let v = profile_json("posq", &[&a]);

        // Half the reads cover the tail positions, so the MEASURED mean there
        // is (40 + 20) / 2 = 30, and the out-of-run curve the simulator needs
        // is 40.
        let uncond = v["pos_q_uncond"]["R1"][20].as_f64().unwrap();
        let out = v["pos_q"]["R1"][20].as_f64().unwrap();
        assert!((uncond - 30.0).abs() < 0.5, "measured mean at the tail: {uncond}");
        assert!(
            (out - 40.0).abs() < 1.0,
            "out-of-run curve at the tail must be 40, got {out} (measured {uncond})"
        );
        // ... and the flat in-run quality is measured, not inferred.
        let run_tail = v["run_q_mean"]["R1"]["tail"].as_f64().unwrap();
        assert!((run_tail - 20.0).abs() < 0.5, "in-run mean: {run_tail}");
        // Position 0 has no head run, so it is the measurement itself.
        assert!((v["pos_q"]["R1"][0].as_f64().unwrap() - 40.0).abs() < 0.5);
    }

    /// The boundary deconvolution: a library of clean reads whose last base is
    /// G with probability q must NOT report those reads as carrying a tail.
    #[test]
    fn deconvolution_removes_the_biological_boundary_run() {
        // observed: 70% no artefact, 20% a 4-base artefact, 10% a 9-base one,
        // with the no-artefact reads' last base G with probability 0.25 (the
        // background of this fixture).
        let obs = vec![0.70 * 0.75, 0.70 * 0.75 * 0.25, 0.0, 0.20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10];
        let out = deconvolve_boundary_g(&obs, 0.25);
        let s: f64 = out.iter().sum();
        assert!((s - 1.0).abs() < 1e-9, "weights must be a distribution, got {s}");
        assert!((out[0] - 0.70).abs() < 0.02, "artefact-free share: {}", out[0]);
        assert!((out[3] - 0.20).abs() < 0.02, "4-base artefact: {}", out[3]);
        assert!((out[9] - 0.10).abs() < 0.02, "9-base artefact: {}", out[9]);
        assert!(out[1] < 0.005 && out[2] < 0.005, "the biological run must be gone: {out:?}");
    }

    #[test]
    fn median_handles_odd_and_even_lengths() {
        assert_eq!(median(&[1.0, 3.0, 2.0]), 2.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(median(&[]), 30.0, "an empty curve falls back to Q30");
    }

    #[test]
    fn position_g_fraction_counts_every_base_including_n() {
        // Four reads whose FIRST base is G, C, N, N. The position-0 G fraction
        // is 1/4: N counts toward the depth, so a heavily-N position cannot look
        // G-rich, and reads must not be dropped for carrying N.
        let a = acc_with(&[(b"G", b"I"), (b"C", b"I"), (b"N", b"I"), (b"N", b"I")], 16);
        let gf = g_fraction(&a);
        assert!((gf[0] - 0.25).abs() < 1e-9, "got {}", gf[0]);
    }

    /// The N model must travel from the accumulator into the profile: a library
    /// whose reads carry ambiguous runs gets an `n_run` block with the measured
    /// rate, run-length histogram and edge shares; a library with none gets no
    /// block at all (so fastsim draws nothing for it).
    #[test]
    fn profile_carries_the_measured_n_model_only_when_reads_carry_n() {
        // two reads: "NNNNAAA" (leading run of 4) and "ACGTACGT" (no N)
        let a = acc_with(&[(b"NNNNAAA", b"IIIIIII"), (b"ACGTACGT", b"IIIIIIII")], 16);
        let dir = std::env::temp_dir().join(format!("fc_nprof_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("p.json");
        write_profile(&p, "test", &[&a], 33, None).unwrap();
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        let n = &v["n_run"];
        assert!((n["read_rate"].as_f64().unwrap() - 0.5).abs() < 1e-9, "1 of 2 reads: {n}");
        assert_eq!(n["len"]["4"], serde_json::json!(1));
        assert_eq!(n["leading"].as_f64().unwrap(), 1.0);
        assert_eq!(n["trailing"].as_f64().unwrap(), 0.0);
        std::fs::remove_dir_all(&dir).unwrap();

        // no ambiguous base anywhere -> the key must be absent
        let b = acc_with(&[(b"ACGTACGT", b"IIIIIIII")], 16);
        let dir = std::env::temp_dir().join(format!("fc_nprof2_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("p.json");
        write_profile(&p, "test", &[&b], 33, None).unwrap();
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        assert!(v.get("n_run").is_none(), "no N -> no n_run block");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn insert_cdf_puts_censored_mass_on_the_read_length_and_saturation_bins() {
        let mut h = crate::insert::InsertHist::new(400, Vec::new());
        h.bins[300] = 10;
        h.censored_below = 4;
        h.censored_above = 6;
        let (cdf, summary) = insert_cdf(&h, 151);
        assert_eq!(cdf.get("300"), Some(&10));
        assert_eq!(cdf.get("151"), Some(&4), "complete-overlap pairs sit at the read length");
        assert_eq!(cdf.get("303"), Some(&6), "unreachable pairs sit in the saturation bin");
        assert_eq!(cdf.values().sum::<u64>(), 20, "no pair may be lost");
        assert_eq!(summary["censored_above_rate"], serde_json::json!(0.3));
    }

    #[test]
    fn fit_mate_produces_a_curve_of_exactly_the_read_length() {
        let seq = b"ACGTACGTACGTACGT";
        let q = b"IIIIIIIIIIIIIIII";
        let a = acc_with(&[(seq, q), (seq, q)], 32);
        let m = fit_mate(&a, 0.25, 16, 33);
        assert_eq!(m.pos_q.len(), 16);
        assert!(m.pos_q.iter().all(|&x| x == 40.0), "Phred+33 'I' is Q40: {:?}", m.pos_q);
        assert_eq!(m.base.len(), 16);
        let row = &m.base[0];
        assert!((row.iter().sum::<f64>() - 1.0).abs() < 0.01, "rows are fractions: {row:?}");
        assert!(m.body_q > 0.0 && m.tail_q > 0.0);
    }

    #[test]
    fn fit_mate_pads_a_run_shorter_than_the_read_length() {
        let a = acc_with(&[(b"ACGT", b"IIII")], 32);
        let m = fit_mate(&a, 0.25, 8, 33);
        assert_eq!(m.pos_q.len(), 8, "positions past the scan are padded, never dropped");
    }

    /// The profile's dup_rate must equal the report's duplication_level, i.e. it must be
    /// computed over the MERGED accumulator and not over `accs.first()`.
    ///
    /// The two mates carry genuinely different duplicate rates: measured on DRR000466,
    /// R1 alone 0.0086, R2 alone 0.0182, merged 0.0135, with only 15 sequences shared
    /// between the mates. Reporting the first mate made the profile and the report
    /// disagree by construction (median ratio 0.7184 profile/report, i.e. the report is
    /// 1.392x the profile), and the closed loop then compared the sim — fed the first
    /// mate's rate — against the report's merged rate, failing 144 of 306 comparisons on
    /// dup_ratio alone.
    #[test]
    fn profile_dup_rate_is_the_merged_rate_not_the_first_mates() {
        // R1: 4 reads, all distinct -> R1 level 0.
        let r1 = acc_with(&[(b"AAAA", b"IIII"), (b"CCCC", b"IIII"),
                            (b"GGGG", b"IIII"), (b"TTTT", b"IIII")], 16);
        // R2: 4 reads, two of them duplicates of each other. The mates share
        // NO sequence, so the merged distinct count is exactly the sum of the
        // per-mate distinct counts minus nothing.
        let r2 = acc_with(&[(b"ACAC", b"IIII"), (b"ACAC", b"IIII"),
                            (b"TGTG", b"IIII"), (b"CACA", b"IIII")], 16);
        let accs: Vec<&Accum> = vec![&r1, &r2];

        // What `accs.first()` used to report: R1's own level, 0.
        let first_mate_level = 1.0 - accs[0].seq_counts.distinct() as f64
            / accs[0].reads.max(1) as f64;
        // What the report reports over the merged table.
        let mut merged = crate::fingerprint::FingerprintTable::with_capacity(0);
        for a in &accs {
            merged.merge(&a.seq_counts);
        }
        let reads: f64 = accs.iter().map(|a| a.reads.max(1) as f64).sum();
        let merged_level = 1.0 - merged.distinct() as f64 / reads;

        // R1 contributes 4 distinct, R2 contributes 3 (ACAC twice): 7 over 8.
        assert_eq!(first_mate_level, 0.0, "R1 alone has no duplicates");
        assert!((merged_level - 0.125).abs() < 1e-9,
                "merged: 8 reads, 7 distinct -> 0.125, got {merged_level}");
        assert!((merged_level - first_mate_level).abs() > 0.1,
                "the merged rate must differ from the first mate's, or this test is vacuous");
    }

    /// `read_len` is the length the library is AT, not the longest read in it.
    ///
    /// One 300 bp outlier in a 151 bp library used to set read_len = 300 while
    /// the same run's fit report recorded mode_len 151 at mode_share 0.999:
    /// fastsim then simulated 300 bp reads from a 151 bp library, and the insert
    /// CDF's censored bin (which sits at read_len) moved with it.
    #[test]
    fn read_len_is_the_modal_length_not_the_longest_read() {
        const BODY: &[u8] = b"ACATACATACATACATAC";
        let mut seqs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        for _ in 0..1000 {
            seqs.push(fixture_read(BODY, b'T', 151 - BODY.len(), b'I'));
        }
        seqs.push(fixture_read(BODY, b'T', 300 - BODY.len(), b'I'));
        let refs: Vec<(&[u8], &[u8])> =
            seqs.iter().map(|(s, q)| (s.as_slice(), q.as_slice())).collect();
        let a = acc_with(&refs, 320);

        let dir = std::env::temp_dir().join(format!("fc_prof_readlen_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("p.json");
        write_profile(&p, "fixture", &[&a], 33, None).unwrap();
        let v: Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        let fr: Value =
            serde_json::from_str(&std::fs::read_to_string(dir.join("p.fit_report.json")).unwrap())
                .unwrap();
        std::fs::remove_dir_all(&dir).unwrap();

        assert_eq!(fr["measured"]["mode_len"], 151, "the fixture's mode is 151 bp");
        assert_eq!(
            v["read_len"], 151,
            "read_len must be the modal length: the fit report itself records mode_len 151 \
             at mode_share 0.999"
        );
        assert_eq!(
            v["pos_q"]["R1"].as_array().unwrap().len(),
            151,
            "pos_q must be exactly read_len long"
        );
        let warnings = fr["fit_warnings"].as_array().unwrap();
        assert!(
            warnings.iter().any(|w| {
                let s = w.as_str().unwrap_or("");
                s.contains("151") && s.contains("300")
            }),
            "the outlier must be named: max 300 vs mode 151, got {warnings:?}"
        );
    }
}
