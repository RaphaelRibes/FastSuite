//! Fit a fastsim profile directly from fastcheck's measurements.
//!
//! This is the Rust port of `bench/scripts/qc_to_profile.py`: it turns the exact
//! counters fastcheck already holds into the profile JSON fastsim consumes, so
//! the coupled path is `fastcheck --profile p.json` -> `fastsim --profile p.json`
//! with no Python step in between.
//!
//! What is MEASURED, ESTIMATED and INVENTED is unchanged from the Python fitter
//! and is recorded in the profile's `_comment`:
//!
//!   MEASURED  read_len, gc (fragment), the per-position mean-Q curve, the
//!             duplicate rate, the adapter readthrough rate, the per-position
//!             base composition, and the insert distribution (--insert-hist)
//!   ESTIMATED the head/tail G-run distributions, from the per-position G
//!             fraction: P(run >= k) ~= g_frac(k-1) - g_mean
//!   INVENTED  the insert SHAPE when no insert module ran, and the mate
//!             tail correlation.

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, Result};
use serde_json::{json, Map, Value};

use crate::accum::Accum;
use crate::insert::InsertHist;

/// Bases of per-position composition carried into the profile: enough to hold a
/// located end bias (a dropout is a few bases), at ~4 floats each.
const BASE_POSITIONS: usize = 60;

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

/// Head (from_start) or tail run-length distribution, as a weight map.
///
/// P(run >= k) ~= g_frac(k-1) - g_mean, so P(run == k) is the drop between
/// successive values. Two physical constraints are enforced: P(run >= k) is
/// non-increasing in k (a run of length k+1 implies one of k, so a spike is
/// impossible and is clamped away by a cumulative minimum), and runs are capped
/// at 12 with the excess folded back into "no run".
fn run_length_cdf_from_gfraction(gf: &[f64], g_mean: f64, from_start: bool) -> BTreeMap<u32, f64> {
    const MAX_RUN: usize = 12;
    let l = gf.len();
    let at = |i: usize| if from_start { gf[i] } else { gf[l - 1 - i] };

    let upto = l.min(MAX_RUN + 1);
    let raw: Vec<f64> = (0..upto).map(|k| (at(k) - g_mean).max(0.0)).collect();
    let mut monotone = Vec::with_capacity(raw.len());
    let mut running_min = f64::INFINITY;
    for v in &raw {
        running_min = running_min.min(*v);
        monotone.push(running_min);
    }

    let mut dist: BTreeMap<u32, f64> = BTreeMap::new();
    dist.insert(0, (1.0 - monotone[0]).max(0.0));
    for k in 1..monotone.len() {
        let v = monotone[k - 1] - monotone[k];
        if v > 0.0 {
            *dist.entry(k as u32).or_insert(0.0) += v;
        }
    }
    if l > MAX_RUN + 1 {
        *dist.entry(0).or_insert(0.0) += (at(MAX_RUN + 1) - g_mean).max(0.0);
    }

    let total: f64 = dist.values().sum();
    if total <= 0.0 {
        return BTreeMap::from([(0, 100.0)]);
    }
    dist.into_iter()
        .filter(|(_, v)| v / total > 0.0005)
        .map(|(k, v)| (k, v / total * 10000.0))
        .collect()
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
                .map(|(b, &c)| ((b as i64 - offset as i64).max(0) as f64) * c as f64)
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
    pos_q: Vec<i64>,
    meanq: BTreeMap<String, i64>,
    body_q: f64,
    tail_q: f64,
    base: Vec<Vec<f64>>,
}

fn fit_mate(acc: &Accum, g_mean: f64, read_len: usize, offset: u8) -> Mate {
    let gf = g_fraction(acc);
    let head = run_length_cdf_from_gfraction(&gf, g_mean, true);
    let tail = run_length_cdf_from_gfraction(&gf, g_mean, false);

    let mut q = per_position_mean_q(acc, offset);
    while q.len() < read_len {
        q.push(*q.last().unwrap_or(&30.0)); // pad if the run is shorter
    }
    let pos_q: Vec<i64> = q[..read_len].iter().map(|x| x.round() as i64).collect();
    let body_q = median(&q[..((read_len as f64 * 0.8) as usize).max(1)]);
    let tail_q = median(&q[((read_len as f64 * 0.9) as usize).min(q.len())..]);

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

    Mate { head, tail, pos_q, meanq: read_meanq(acc, offset), body_q, tail_q, base }
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
    let read_len = (qc1.len_hist.len().saturating_sub(1)) as usize;
    if read_len == 0 {
        anyhow::bail!("--profile: no reads were scanned, so there is nothing to fit");
    }
    let mode_len = qc1
        .len_hist
        .iter()
        .enumerate()
        .max_by_key(|(_, &c)| c)
        .map(|(i, _)| i)
        .unwrap_or(read_len);
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
        format!("fitted from {name} by fastcheck --profile (as-measured): head/tail estimated, \
                 insert distribution MEASURED by --insert-hist, mate tail correlation invented")
    } else {
        format!("fitted from {name} by fastcheck --profile (as-observed): head/tail estimated, \
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
    prof.insert("read_meanq".into(), Value::Object(read_meanq));
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
            a.observe(s, q);
        }
        a
    }

    #[test]
    fn run_length_cdf_is_monotone_shaped_and_normalised() {
        // A library whose 3' end is G-rich: g_frac rises above the mean at the
        // last positions, which is exactly what a tail run looks like.
        let mut gf = vec![0.25f64; 20];
        gf[19] = 0.60;
        gf[18] = 0.45;
        let d = run_length_cdf_from_gfraction(&gf, 0.25, false);
        let total: f64 = d.values().sum();
        assert!((total - 10000.0).abs() < 1.0, "weights must sum to 10000, got {total}");
        assert!(d.contains_key(&0), "a run of 0 is always a possible outcome");
        assert!(d[&0] < 10000.0, "some mass must sit on a real run");
        assert!(!d.contains_key(&13), "runs are capped at 12");
    }

    #[test]
    fn run_length_cdf_folds_the_beyond_cap_excess_into_no_run() {
        // A flat G-rich tail: every position is above the mean, so most of the
        // mass belongs to runs longer than the cap and must land on 0.
        let gf = vec![0.50f64; 20];
        let d = run_length_cdf_from_gfraction(&gf, 0.25, false);
        let total: f64 = d.values().sum();
        assert!((total - 10000.0).abs() < 1.0);
        assert!(d[&0] > 5000.0, "the capped excess belongs in the no-run bin");
    }

    #[test]
    fn cdf_with_no_gt_baseline_is_all_no_run() {
        let gf = vec![0.25f64; 20];
        let d = run_length_cdf_from_gfraction(&gf, 0.25, false);
        assert_eq!(d.len(), 1);
        assert_eq!(d[&0], 10000.0);
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
        assert!(m.pos_q.iter().all(|&x| x == 40), "Phred+33 'I' is Q40: {:?}", m.pos_q);
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
}
