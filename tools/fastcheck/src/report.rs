//! Output: `fastsuite.report.v1` JSON, per-position TSV, and the summary line.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use serde_json::json;

use crate::accum::{quantile, Accum, OffsetVerdict, QBINS};
use crate::reader::ScanResult;

pub struct FileReport {
    pub path: PathBuf,
    pub acc: Accum,
    pub scan: ScanResult,
}

/// Total bases that are neither A, C, G nor T.
fn bases_other(a: &Accum) -> u64 {
    let acgt: u64 = (0..a.acgt_hist.len()).map(|i| i as u64 * a.acgt_hist[i]).sum();
    a.bases - acgt
}

fn gc_fraction(a: &Accum) -> f64 {
    let acgt: u64 = (0..a.acgt_hist.len()).map(|i| i as u64 * a.acgt_hist[i]).sum();
    if acgt == 0 {
        return 0.0;
    }
    let gc: u64 = (0..a.gc_hist.len()).map(|i| i as u64 * a.gc_hist[i]).sum();
    gc as f64 / acgt as f64
}

fn mean_len(a: &Accum) -> f64 {
    if a.reads == 0 {
        0.0
    } else {
        a.bases as f64 / a.reads as f64
    }
}

/// Mean over reads of the per-read mean quality. Exact: the per-read bin is
/// the floor of the raw byte mean, and `floor(x) - k == floor(x - k)` for an
/// integer offset `k`.
fn mean_read_mean_q(a: &Accum, offset: u8) -> f64 {
    if a.reads == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for (b, &c) in a.meanq_raw_hist.iter().enumerate() {
        let q = b as i64 - offset as i64;
        sum += q.max(0) as f64 * c as f64;
    }
    sum / a.reads as f64
}

fn frac_reads_at_least(a: &Accum, offset: u8, min_q: i64) -> f64 {
    if a.reads == 0 {
        return 0.0;
    }
    let mut n = 0u64;
    for (b, &c) in a.meanq_raw_hist.iter().enumerate() {
        if b as i64 - offset as i64 >= min_q {
            n += c;
        }
    }
    n as f64 / a.reads as f64
}

fn hist_json(hist: &[u64]) -> serde_json::Value {
    let mut m = serde_json::Map::new();
    for (i, &c) in hist.iter().enumerate() {
        if c > 0 {
            m.insert(i.to_string(), json!(c));
        }
    }
    serde_json::Value::Object(m)
}

fn byte_range(census: &[u64; QBINS]) -> (u64, u64) {
    let min = (0..QBINS).find(|&b| census[b] > 0).unwrap_or(0) as u64;
    let max = (0..QBINS).rev().find(|&b| census[b] > 0).unwrap_or(0) as u64;
    (min, max)
}

/// Per-position quality summary from the exact histogram at that position.
fn pos_stats(row: &[u64; QBINS], offset: u8) -> (f64, f64, f64, f64) {
    let depth: u64 = row.iter().sum();
    if depth == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    // Shift into Phred space once, then reuse the shared nearest-rank quantile.
    let mut q = vec![0u64; QBINS];
    let mut sum = 0.0;
    for (b, &c) in row.iter().enumerate() {
        let phred = (b as i64 - offset as i64).max(0) as usize;
        let phred = phred.min(QBINS - 1);
        q[phred] += c;
        sum += phred as f64 * c as f64;
    }
    (
        sum / depth as f64,
        quantile(&q, 0.10),
        quantile(&q, 0.50),
        quantile(&q, 0.90),
    )
}

fn write_out(path: &Path, body: &str) -> Result<()> {
    // D-4: `-` means "the stream", for reports as well as FASTQ.
    if path.as_os_str() == "-" {
        use std::io::Write;
        let mut out = std::io::stdout();
        out.write_all(body.as_bytes())
            .and_then(|_| out.flush())
            .context("writing the report to stdout")?;
        return Ok(());
    }
    std::fs::write(path, body).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

pub fn write_report(
    path: &Path,
    files: &[FileReport],
    total: &Accum,
    verdict: &OffsetVerdict,
    elapsed_s: f64,
    partial: Option<&str>,
    positions_tsv: Option<&Path>,
    insert: Option<(&crate::insert::InsertHist, u64)>,
) -> Result<()> {
    let offset = verdict.offset;
    let (byte_min, byte_max) = byte_range(&total.qual_byte);
    let bytes_read: u64 = files.iter().map(|f| f.scan.bytes_read).sum();

    let mut notes = vec![
        format!("Phred offset: {} ({})", offset, verdict.reason),
        format!(
            "%-GC denominators exclude N: {} of {} bases are non-ACGT",
            bases_other(total),
            total.bases
        ),
    ];
    if verdict.malformed {
        notes.push(format!(
            "quality bytes below 33 observed (min {byte_min}): the file is not valid Phred+33"
        ));
    }
    if total.reads_over_max_len > 0 {
        notes.push(format!(
            "{} reads longer than --max-read-len are counted in the summary but excluded from \
             the per-position matrices",
            total.reads_over_max_len
        ));
    }
    let modes: std::collections::BTreeSet<&str> = files.iter().map(|f| f.scan.decompress).collect();
    notes.push(format!("decompression modes: {}", modes.into_iter().collect::<Vec<_>>().join(", ")));
    if let Some(p) = positions_tsv {
        notes.push(format!("per-position matrices written to {}", p.display()));
    }

    let per_file: Vec<serde_json::Value> = files
        .iter()
        .map(|f| {
            json!({
                "name": f.path.to_string_lossy(),
                "reads": f.acc.reads,
                "bases": f.acc.bases,
                "mean_len": mean_len(&f.acc),
                "gc_fraction": gc_fraction(&f.acc),
                "mean_read_mean_q": mean_read_mean_q(&f.acc, offset),
                "decompression": f.scan.decompress,
                "bytes_read": f.scan.bytes_read,
                "stopped": f.scan.stopped,
            })
        })
        .collect();

    let mut metrics = serde_json::Map::new();
    let mut put = |k: &str, v: serde_json::Value| {
        metrics.insert(k.to_string(), v);
    };
    put("files", json!(files.len()));
    put("reads", json!(total.reads));
    put("bases", json!(total.bases));
    put("bases_acgt", json!(total.bases - bases_other(total)));
    put("bases_non_acgt", json!(bases_other(total)));
    put("reads_with_non_acgt", json!(total.reads_with_n));
    put("reads_over_max_read_len", json!(total.reads_over_max_len));
    // min_len is the first non-empty histogram bin — NOT a low quantile: the
    // 1.0 quantile is the maximum, which is exactly the bug this line replaces.
    put(
        "min_len",
        json!(total.len_hist.iter().position(|&c| c > 0).unwrap_or(0)),
    );
    put("max_len", json!(total.len_hist.len().saturating_sub(1)));
    put("median_len", json!(quantile(&total.len_hist, 0.5)));
    put("mean_len", json!(mean_len(total)));
    put("gc_fraction", json!(gc_fraction(total)));
    put("mean_read_mean_q", json!(mean_read_mean_q(total, offset)));
    put("reads_q20_frac", json!(frac_reads_at_least(total, offset, 20)));
    put("reads_q30_frac", json!(frac_reads_at_least(total, offset, 30)));
    put("phred_offset", json!(offset));
    put("phred_offset_certain", json!(if verdict.certain { 1 } else { 0 }));
    put("quality_byte_min", json!(byte_min));
    put("quality_byte_max", json!(byte_max));
    put("partial", json!(if partial.is_some() { 1 } else { 0 }));
    put("bytes_read", json!(bytes_read));
    put("elapsed_s", json!(elapsed_s));
    // Only when --insert-hist ran: a consumer must never read an absent metric
    // as a measured zero.
    if let Some((h, _)) = &insert {
        put("insert_measured_pairs", json!(h.measured));
        put("insert_censored_below_read_len", json!(h.censored_below));
        put("insert_censored_no_overlap", json!(h.censored_above));
        put("insert_mean", json!(h.mean()));
        put("insert_sd", json!(h.sd()));
        put("insert_p50", json!(h.quantile_public(0.50)));
    }

    let mut composition = serde_json::Map::new();
    for (k, idx) in [("a", 0usize), ("c", 1), ("g", 2), ("t", 3), ("other", 4)] {
        let n: u64 = (0..total.max_pos).map(|p| total.base_at(p, idx)).sum();
        composition.insert(k.to_string(), json!(n));
    }

    // Per-read mean quality in Phred space, from the exact raw-byte bins.
    let mut meanq_json = serde_json::Map::new();
    for (b, &c) in total.meanq_raw_hist.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let q = (b as i64 - offset as i64).max(0);
        let e = meanq_json.entry(q.to_string()).or_insert_with(|| json!(0u64));
        *e = json!(e.as_u64().unwrap_or(0) + c);
    }

    // D-6: every assumption this run rests on, stated as data rather than as
    // prose, so a consumer can act on it.
    let mut assumptions: Vec<serde_json::Value> = Vec::new();
    if !verdict.certain {
        assumptions.push(json!({
            "what": "phred_offset",
            "value": offset,
            "basis": verdict.reason,
            "consequence": "quality-derived metrics use this offset; they are wrong by 31 if the file is Phred+64",
        }));
    }

    let doc = json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastcheck",
        "tool_version": env!("CARGO_PKG_VERSION"),
        "command": std::env::args().collect::<Vec<_>>().join(" "),
        "mode": "qc",
        "inputs": { "files": files.iter().map(|f| f.path.to_string_lossy().to_string()).collect::<Vec<_>>() },
        "elapsed_s": elapsed_s,
        "metrics": serde_json::Value::Object(metrics),
        "assumptions": assumptions,
        "per_file": per_file,
        "base_composition": serde_json::Value::Object(composition),
        "insert_hist": match &insert {
            Some((h, pairs)) => {
                // The sample is a fraction of the input, so the histogram is
                // rescaled to the input's pair count and the factor is stated.
                let scale = if h.pairs_sampled > 0 {
                    Some(*pairs as f64 / h.pairs_sampled as f64)
                } else { None };
                h.to_json(Some(*pairs), scale)
            }
            None => serde_json::Value::Null,
        },
        "read_length_hist": hist_json(&total.len_hist),
        "read_mean_q_hist": serde_json::Value::Object(meanq_json),
        "partial": {
            "is_partial": partial.is_some(),
            "reason": partial,
            "bytes_read": bytes_read,
        },
        "notes": notes,
    });

    let text = serde_json::to_string_pretty(&doc)? + "\n";
    write_out(path, &text)?;
    Ok(())
}

pub fn write_positions_tsv(path: &Path, total: &Accum, offset: u8) -> Result<()> {
    let n_pos = total.covered_positions();
    let mut out = String::with_capacity(64 * (n_pos + 1));
    out.push_str("pos\ta\tc\tg\tt\tn\tdepth\tmean_q\tq10\tq50\tq90\n");
    for pos in 0..n_pos {
        let row = total.qual_at(pos);
        let depth: u64 = row.iter().sum();
        let (mean, q10, q50, q90) = pos_stats(&row, offset);
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.0}\t{:.0}\t{:.0}\n",
            pos + 1,
            total.base_at(pos, 0),
            total.base_at(pos, 1),
            total.base_at(pos, 2),
            total.base_at(pos, 3),
            total.base_at(pos, 4),
            depth,
            mean,
            q10,
            q50,
            q90
        ));
    }
    write_out(path, &out)?;
    Ok(())
}
