//! Output: `fastsuite.report.v1` JSON, per-position TSV, and the summary line.

use std::path::{Path, PathBuf};
use crate::rich;

use anyhow::{Context, Result};
use serde_json::json;

use crate::accum::{quantile, Accum, OffsetVerdict, QBYTES, QBINS};
use crate::cli::Modules;
use crate::reader::ScanResult;

pub struct FileReport {
    pub path: PathBuf,
    pub acc: Accum,
    pub scan: ScanResult,
}

/// Total bases that are neither A, C, G nor T.
pub(crate) fn r4(x: f64) -> f64 {
    (x * 10000.0).round() / 10000.0
}

pub(crate) fn bases_other(a: &Accum) -> u64 {
    let acgt: u64 = (0..a.acgt_hist.len()).map(|i| i as u64 * a.acgt_hist[i]).sum();
    a.bases - acgt
}

/// The input paths, for the modules that re-read the files (overrep second-pass).
fn files_paths(files: &[FileReport]) -> Vec<PathBuf> {
    files.iter().map(|f| f.path.clone()).collect()
}

pub(crate) fn gc_fraction(a: &Accum) -> f64 {
    let acgt: u64 = (0..a.acgt_hist.len()).map(|i| i as u64 * a.acgt_hist[i]).sum();
    if acgt == 0 {
        return 0.0;
    }
    let gc: u64 = (0..a.gc_hist.len()).map(|i| i as u64 * a.gc_hist[i]).sum();
    gc as f64 / acgt as f64
}

pub(crate) fn mean_len(a: &Accum) -> f64 {
    if a.reads == 0 {
        0.0
    } else {
        a.bases as f64 / a.reads as f64
    }
}

/// Mean over reads of the per-read mean quality. Exact: the per-read bin is
/// the floor of the raw byte mean, and `floor(x) - k == floor(x - k)` for an
/// integer offset `k`.
pub(crate) fn mean_read_mean_q(a: &Accum, offset: u8) -> f64 {
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

pub(crate) fn hist_json(hist: &[u64]) -> serde_json::Value {
    let mut m = serde_json::Map::new();
    for (i, &c) in hist.iter().enumerate() {
        if c > 0 {
            m.insert(i.to_string(), json!(c));
        }
    }
    serde_json::Value::Object(m)
}

fn byte_range(census: &[u64; QBYTES]) -> (u64, u64) {
    let min = (0..QBYTES).find(|&b| census[b] > 0).unwrap_or(0) as u64;
    let max = (0..QBYTES).rev().find(|&b| census[b] > 0).unwrap_or(0) as u64;
    (min, max)
}

/// Per-position quality summary from the exact histogram at that position.
///
/// The row is indexed by `byte - 33` (128 cells), so the byte is recovered before the
/// Phred offset is applied. A valid FASTQ can only land in those cells, which is why
/// this is the same mean and the same q10/q50/q90 it always was.
pub(crate) fn pos_stats(row: &[u64; QBINS], offset: u8) -> (f64, f64, f64, f64) {
    let depth: u64 = row.iter().sum();
    if depth == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    // Shift into Phred space once, then reuse the shared nearest-rank quantile.
    let mut q = vec![0u64; QBINS];
    let mut sum = 0.0;
    for (b, &c) in row.iter().enumerate() {
        let phred = ((b + crate::accum::QMIN as usize) as i64 - offset as i64).max(0) as usize;
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

/// Write `body` to `path`, or to stdout when `path` is `-`.
///
/// The one writer every report output uses: `--json`, `--html`, `--profile` and
/// `--positions-tsv` all stream to the same place, so a `-` cannot mean a file named `-`
/// in one of them.
pub(crate) fn write_out(path: &Path, body: &str) -> Result<()> {
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

/// CPU-seconds consumed by this process so far, for `busy_cores`.
///
/// Copied from `tools/fastprep/src/stats.rs` (D-1: no shared crate). `getrusage` on a
/// single process is cheap and needs no extra dependency: `libc` is already here for
/// `memory.rs`. A failing call reports 0 rather than a fabricated ratio.
fn cpu_seconds() -> f64 {
    let mut ru = libc::rusage {
        ru_utime: libc::timeval { tv_sec: 0, tv_usec: 0 },
        ru_stime: libc::timeval { tv_sec: 0, tv_usec: 0 },
        ..unsafe { std::mem::zeroed() }
    };
    // SAFETY: getrusage only writes into the struct we hand it; a failure leaves the
    // zeroed value, so the metric degrades to 0 rather than lying.
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut ru) } != 0 {
        return 0.0;
    }
    let secs = |tv: libc::timeval| tv.tv_sec as f64 + tv.tv_usec as f64 / 1e6;
    secs(ru.ru_utime) + secs(ru.ru_stime)
}

/// The cores this run kept busy: CPU-seconds per wall-second, rounded to 3 decimals.
///
/// `metrics.threads` says how many threads the run was ASKED for; a report that carried
/// only that could claim 8 threads while one core did the work. A zero-length run reports
/// 0 rather than an infinity, so consumers that sum the field do not break.
pub(crate) fn busy_cores(elapsed_s: f64) -> f64 {
    if elapsed_s <= 0.0 {
        return 0.0;
    }
    let busy = cpu_seconds() / elapsed_s;
    (busy * 1000.0).round() / 1000.0
}

fn base_position_json(a: &Accum) -> serde_json::Value {
    let n = measured_positions(a);
    let mut rows = serde_json::Map::with_capacity(n);
    for pos in 0..n {
        rows.insert(
            (pos + 1).to_string(),
            serde_json::json!([
                a.base_at(pos, 0),   // A
                a.base_at(pos, 1),   // C
                a.base_at(pos, 2),   // G
                a.base_at(pos, 3),   // T
            ]),
        );
    }
    serde_json::Value::Object(rows)
}

/// Positions the per-position matrices cover. Measured from the quality rows (grown
/// to the longest read, and every read contributes to both planes), falling back to
/// the composition rows when the quality module did not run. `q_pos > 0` whenever a
/// read was measured, so the fallback only ever applies to a gated run.
fn measured_positions(a: &Accum) -> usize {
    let q = a.covered_positions();
    if q > 0 {
        q
    } else {
        a.base_pos.iter().filter(|r| r.iter().any(|&c| c > 0)).count()
    }
}

pub fn write_report(
    path: &Path,
    files: &[FileReport],
    total: &Accum,
    verdict: &OffsetVerdict,
    elapsed_s: f64,
    partial: Option<&str>,
    positions_tsv: Option<&Path>,
    insert: Option<(&crate::insert::InsertHist, Option<u64>)>,
    processed: Option<&Accum>,
    overrep: crate::overrep::OverrepMode,
    overrep_min_count: u32,
    // The resolved `-j` count: the overrep second pass re-reads the files and
    // builds its BGZF inflate pool from it, so the one mode that pays a second read does
    // not pay it with one thread because a literal said so.
    threads: usize,
    // CPU-seconds per wall-second over the scan: `threads` is the request, this is what
    // it bought (fastprep's `busy_cores`).
    busy_cores: f64,
    adapter_content: Option<(&crate::adapter_pos::AdapterPos, f64, usize)>,
    modules: Option<&Modules>,
    tile: Option<&crate::tile::TileQual>,
) -> Result<()> {
    let offset = verdict.offset;
    let (byte_min, byte_max) = byte_range(&total.qual_byte);
    let bytes_read: u64 = files.iter().map(|f| f.scan.bytes_read).sum();

    // Module gating. `None` means no gating flag was given and every field below is
    // emitted exactly as it was before `--only` existed; `Some(m)` emits only what the
    // enabled modules measured. An absent field means "not measured" — never "measured
    // as zero" — which is the rule the adapter_content block already follows. The one
    // always-reported number is `phred_offset`: with the quality module off it is the
    // stated default 33, and the report says so in `assumptions`.
    let on = |f: fn(&Modules) -> bool| modules.map_or(true, |m| f(m));
    let (
        bases_on,
        q_on,
        len_on,
        gc_on,
        hm_on,
        dup_on,
        polyx_on,
        kmer_on,
    ) = (
        on(|m| m.bases),
        on(|m| m.quality),
        on(|m| m.length),
        on(|m| m.gc),
        on(|m| m.heatmap),
        on(|m| m.dup),
        on(|m| m.polyx),
        on(|m| m.kmer),
    );

    let mut notes = vec![
        if q_on {
            format!("Phred offset: {} ({})", offset, verdict.reason)
        } else {
            format!(
                "Phred offset: {offset} (stated default: the quality module did not run, \
                 so no byte census decided it)"
            )
        },
    ];
    if gc_on {
        notes.push(format!(
            "%-GC denominators exclude N: {} of {} bases are non-ACGT",
            bases_other(total),
            total.bases
        ));
    }
    if verdict.malformed {
        notes.push(format!(
            "quality bytes below 33 observed (min {byte_min}): the file is not valid Phred+33"
        ));
    }
    if total.qual_out_of_range > 0 {
        notes.push(format!(
            "{} quality bytes outside 33..=126 are in the byte census but excluded from the \
             per-position quality rows, which are indexed by byte - 33",
            total.qual_out_of_range
        ));
    }
    if total.reads_over_max_len > 0 {
        notes.push(format!(
            "{} reads longer than --max-read-len are counted in the summary and in \
             base_composition (which covers every base), but excluded from the per-position \
             matrices",
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
    // The resolved `-j` count, so a report says how it was run: it is a run parameter,
    // not a measurement, and two runs of one file at different -j are otherwise identical.
    put("threads", json!(threads));
    // The `-j` above is the REQUEST; this is what it bought. A value near 1 says the run
    // was serial however many cores it was told to use.
    put("busy_cores", json!(busy_cores));
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

    // Rich-report metrics (polyX, head dropout, repeats, k-mers, heatmap).
    let (polyx_m, polyx_block) = rich::polyx_metrics(total);
    let (hg, hd) = rich::head_dropout(total);
    // The SAME effective threshold the `overrep` block uses: one flag, one answer.
    let (rep_m, rep_block) = rich::repeated_metrics(total, overrep_min_count);
    let (kmer_m, kmer_block) = rich::kmer_metrics(total);
    let (hm_m, hm_block) = rich::heatmap_metrics(total);
    for (k, v) in [polyx_m, rep_m, kmer_m, hm_m].into_iter().flatten() {
        put(&k, v);
    }
    put("head_g_frac", json!(r4(hg)));
    put("head_dropout", json!(r4(hd)));

    // The composition covers EVERY base. `base_totals` is summed over the whole read,
    // including the bases past --max-read-len that no per-position row can hold; summing
    // `base_at` over the covered window instead made the block disagree with metrics.bases
    // (10 reads of 8 bp at --max-read-len 4 reported 40 bases of composition against 80).
    let mut composition = serde_json::Map::new();
    for (k, idx) in [("a", 0usize), ("c", 1), ("g", 2), ("t", 3), ("other", 4)] {
        composition.insert(k.to_string(), json!(total.base_totals[idx]));
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
    if !verdict.certain && q_on {
        assumptions.push(json!({
            "what": "phred_offset",
            "value": offset,
            "basis": verdict.reason,
            "consequence": "quality-derived metrics use this offset; they are wrong by 31 if the file is Phred+64",
        }));
    }
    if !q_on {
        assumptions.push(json!({
            "what": "quality",
            "value": "not measured",
            "basis": "--no-quality (or --only without quality) skipped the quality accumulators; \
                      the bytes were read for framing but never decoded",
            "consequence": "no quality-derived metric is reported; the Phred offset is the stated \
                            default 33, not a measurement",
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
        "base_position": base_position_json(total),
        "polyx": polyx_block,
        "repeated_sequences": rep_block,
        "overrep": crate::overrep::report(overrep, total, &files_paths(files), threads.max(1), overrep_min_count)?,
        "kmer_spectrum": kmer_block,
        "heatmap": hm_block,
        "processed": match processed {
            Some(pa) => {
                let (pm, pb) = rich::polyx_metrics(pa);
                let (phg, phd) = rich::head_dropout(pa);
                let (prm, prb) = rich::repeated_metrics(pa, overrep_min_count);
                let (pkm, pkb) = rich::kmer_metrics(pa);
                let (phm, phb) = rich::heatmap_metrics(pa);
                let mut pm = pm;
                pm.insert("head_g_frac".into(), json!(r4(phg)));
                pm.insert("head_dropout".into(), json!(r4(phd)));
                for (k, v) in [prm, pkm, phm].into_iter().flatten() {
                    pm.insert(k, v);
                }
                json!({
                    "reads": pa.reads,
                    "bases": pa.bases,
                    "gc_fraction": gc_fraction(pa),
                    "mean_read_mean_q": mean_read_mean_q(pa, offset),
                    "metrics": pm,
                    "polyx": pb,
                    "repeated_sequences": prb,
                    "kmer_spectrum": pkb,
                    "heatmap": phb,
                })
            }
            None => serde_json::Value::Null,
        },
        "insert_hist": match &insert {
            Some((h, pairs)) => {
                // The sample is a fraction of the input, so the histogram is rescaled to
                // the input's pair count — when that count is KNOWN. A sample cap or a
                // budget leaves it null, and so does the scale factor: a factor derived
                // from a different pass's budget rescaled 50 sampled pairs by 0.08.
                let scale = match (pairs, h.pairs_sampled) {
                    (Some(p), sampled) if sampled > 0 => Some(*p as f64 / sampled as f64),
                    _ => None,
                };
                h.to_json(*pairs, scale)
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

    // The adapter-content block is added only when the module ran: an absent key means
    // "not measured", which is different from a measured null (D-6 discipline).
    let mut doc = doc;
    if let Some((ap, error_rate, min_overlap)) = adapter_content {
        doc["adapter_content"] = ap.to_json(error_rate, min_overlap);
    }
    // The per-tile block follows the same rule: --tiles off leaves the accumulator None,
    // so the key is absent rather than a null or a measured zero.
    if let Some(t) = tile {
        doc["per_tile"] = t.to_json(offset);
    }

    // Gating: a module that did not run leaves no key behind. The values for gated-out
    // blocks were computed from empty accumulators, so they are dropped rather than
    // published — a consumer must never read an absent block as a measured zero. The
    // `modules` array is itself only emitted when a gating flag was given: with every
    // module on (or no flag at all) there is nothing to declare, and the default report
    // stays byte for byte what it was.
    if let Some(m) = modules {
        doc["modules"] = json!(m.enabled_names());

        // (key, the module that owns it)
        let blocks: [(&str, bool); 9] = [
            ("base_composition", bases_on),
            ("base_position", bases_on),
            ("read_length_hist", len_on),
            ("read_mean_q_hist", q_on),
            ("polyx", polyx_on),
            ("repeated_sequences", dup_on),
            ("overrep", dup_on),
            ("kmer_spectrum", kmer_on),
            ("heatmap", hm_on),
        ];
        let metrics_keys: [(&str, bool); 30] = [
            ("bases_acgt", gc_on),
            ("bases_non_acgt", gc_on),
            ("gc_fraction", gc_on),
            ("min_len", len_on),
            ("max_len", len_on),
            ("median_len", len_on),
            ("mean_read_mean_q", q_on),
            ("reads_q20_frac", q_on),
            ("reads_q30_frac", q_on),
            ("phred_offset_certain", q_on),
            ("quality_byte_min", q_on),
            ("quality_byte_max", q_on),
            ("head_g_frac", bases_on),
            ("head_dropout", bases_on),
            ("polyx3_max", polyx_on),
            ("polyx5_max", polyx_on),
            ("polyx3_frac_gt_5", polyx_on),
            ("polyx5_frac_gt_1", polyx_on),
            ("distinct_seqs", dup_on),
            ("duplication_level", dup_on),
            ("total_deduplicated_pct", dup_on),
            ("dup_multiplicity_max", dup_on),
            ("dup_singletons", dup_on),
            ("repeated_seq_count", dup_on),
            ("repeated_seq_frac", dup_on),
            ("kmer_top_count", kmer_on),
            ("kmer_top_frac", kmer_on),
            ("heatmap_peak", hm_on),
            ("heatmap_peak_x", hm_on),
            ("heatmap_peak_y", hm_on),
        ];

        let obj = doc.as_object_mut().expect("the report is an object");
        for (k, live) in blocks {
            if !live {
                obj.remove(k);
            }
        }
        if let Some(metrics) = obj.get_mut("metrics").and_then(|v| v.as_object_mut()) {
            for (k, live) in metrics_keys {
                if !live {
                    metrics.remove(k);
                }
            }
        }
        // The processed side ran under the same gate and loses the same fields.
        for side in ["per_file", "processed"] {
            let Some(side) = obj.get_mut(side) else { continue };
            let objects: Vec<&mut serde_json::Map<String, serde_json::Value>> = match side {
                serde_json::Value::Array(a) => a.iter_mut().filter_map(|v| v.as_object_mut()).collect(),
                serde_json::Value::Object(o) => vec![o],
                _ => Vec::new(),
            };
            for o in objects {
                if !q_on {
                    o.remove("mean_read_mean_q");
                }
                if !gc_on {
                    o.remove("gc_fraction");
                }
                for (k, live) in blocks {
                    if !live {
                        o.remove(k);
                    }
                }
                if let Some(metrics) = o.get_mut("metrics").and_then(|v| v.as_object_mut()) {
                    for (k, live) in metrics_keys {
                        if !live {
                            metrics.remove(k);
                        }
                    }
                }
            }
        }
    }

    let text = serde_json::to_string_pretty(&doc)? + "\n";
    write_out(path, &text)?;
    Ok(())
}

pub fn write_positions_tsv(
    path: &Path,
    total: &Accum,
    offset: u8,
    adapter_pos: Option<&crate::adapter_pos::AdapterPos>,
    modules: Option<&Modules>,
) -> Result<()> {
    // Per-position columns travel with their module, the way the adapter columns
    // already do: with `--no-quality` the file carries the composition and no quality
    // column at all, rather than four columns of zeros that read like measurements.
    let on = |f: fn(&Modules) -> bool| modules.map_or(true, |m| f(m));
    let (bases_on, q_on, polyx_on) = (on(|m| m.bases), on(|m| m.quality), on(|m| m.polyx));
    let n_pos = measured_positions(total);
    let reads = total.reads.max(1) as f64;
    // Per-position polyX reach: the fraction of reads whose trailing (3') or
    // leading (5') homopolymer run reaches this position, summed over bases.
    let mut p3_reach = vec![0u64; n_pos];
    let mut p5_reach = vec![0u64; n_pos];
    for b in 0..4 {
        for (len, &c) in total.polyx3[b].iter().enumerate() {
            if c == 0 || len == 0 {
                continue;
            }
            for p in 0..len.min(n_pos) {
                p3_reach[n_pos - 1 - p] += c;
            }
        }
        for (len, &c) in total.polyx5[b].iter().enumerate() {
            if c == 0 || len == 0 {
                continue;
            }
            for p in 0..len.min(n_pos) {
                p5_reach[p] += c;
            }
        }
    }
    // Adapter columns, only with --adapter-content: one per adapter, 1-based positions,
    // 0 when past the end of that adapter's observed rows.
    let adapter_cols: Vec<(String, Vec<u64>)> =
        adapter_pos.map(|ap| ap.tsv_columns()).unwrap_or_default();
    let mut out = String::with_capacity(64 * (n_pos + 1));
    out.push_str("pos");
    if bases_on {
        out.push_str("\ta\tc\tg\tt\tn");
    }
    if bases_on || q_on {
        out.push_str("\tdepth");
    }
    if q_on {
        out.push_str("\tmean_q\tq10\tq50\tq90");
    }
    if bases_on {
        out.push_str("\tg_frac");
    }
    if polyx_on {
        out.push_str("\tpolyx3_frac\tpolyx5_frac");
    }
    for (name, _) in &adapter_cols {
        out.push('\t');
        out.push_str(name);
    }
    out.push('\n');
    for pos in 0..n_pos {
        let row = total.qual_at(pos);
        let depth: u64 = if q_on {
            row.iter().sum()
        } else {
            (0..crate::accum::NT).map(|c| total.base_at(pos, c)).sum()
        };
        out.push_str(&(pos + 1).to_string());
        if bases_on {
            out.push_str(&format!(
                "\t{}\t{}\t{}\t{}\t{}",
                total.base_at(pos, 0),
                total.base_at(pos, 1),
                total.base_at(pos, 2),
                total.base_at(pos, 3),
                total.base_at(pos, 4),
            ));
        }
        if bases_on || q_on {
            out.push_str(&format!("\t{depth}"));
        }
        if q_on {
            let (mean, q10, q50, q90) = pos_stats(&row, offset);
            out.push_str(&format!("\t{mean:.4}\t{q10:.0}\t{q50:.0}\t{q90:.0}"));
        }
        if bases_on {
            let g_frac = if depth > 0 {
                total.base_at(pos, 2) as f64 / depth as f64
            } else {
                0.0
            };
            out.push_str(&format!("\t{g_frac:.4}"));
        }
        if polyx_on {
            out.push_str(&format!(
                "\t{:.4}\t{:.4}",
                p3_reach[pos] as f64 / reads,
                p5_reach[pos] as f64 / reads,
            ));
        }
        for (_, counts) in &adapter_cols {
            out.push('\t');
            out.push_str(&counts.get(pos).copied().unwrap_or(0).to_string());
        }
        out.push('\n');
    }
    write_out(path, &out)?;
    Ok(())
}
