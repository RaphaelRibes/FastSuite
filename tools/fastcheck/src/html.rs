//! A self-contained HTML QC report, in the shape people expect from FastQC.
//!
//! Every number here is an exact count from the same accumulator the JSON report
//! uses — nothing is sampled, binned down or approximated. The charts are inline
//! SVG: no JavaScript, no external asset, no font to fetch, so the report is a
//! single file that survives being emailed, opened offline, or printed.
//!
//! What it deliberately does NOT do is hand out pass/warn/fail verdicts from
//! thresholds nobody can inspect. FastQC's traffic lights encode opinions; this
//! shows the measurement and states, in numbers, the few observations that are
//! worth acting on.

use std::fmt::Write as _;
use std::path::Path;

use anyhow::{Context, Result};

use crate::accum::{quantile, Accum};
use crate::insert::InsertHist;
use crate::report::{bases_other, gc_fraction, mean_len, mean_read_mean_q, pos_stats, FileReport};

const CW: f64 = 940.0; // chart width in SVG user units
const CH: f64 = 250.0; // chart height

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

struct Series<'a> {
    name: &'a str,
    colour: &'a str,
    points: Vec<(f64, f64)>,
}

fn axes(x0: f64, y0: f64, x1: f64, y1: f64, x_max: f64, y_max: f64, x_label: &str, y_label: &str) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        r#"<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y0}" class="ax"/><line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" class="ax"/>"#
    );
    // horizontal gridlines at 0, 50% and 100% of the y range, labelled
    for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = y0 - f * (y0 - y1);
        let _ = write!(
            s,
            r#"<line x1="{x0}" y1="{y:.1}" x2="{x1}" y2="{y:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick end">{}</text>"#,
            x0 - 6.0,
            y + 4.0,
            fmt_num(f * y_max)
        );
    }
    // x ticks: first, middle, last position
    for (v, anchor) in [(0.0, "start"), (x_max / 2.0, "middle"), (x_max, "end")] {
        let x = x0 + (v / x_max) * (x1 - x0);
        let _ = write!(s, r#"<text x="{x:.1}" y="{:.1}" class="tick" text-anchor="{anchor}">{}</text>"#, y0 + 16.0, fmt_num(v));
    }
    let cx_lab = (x0 + x1) / 2.0;
    let cy_lab = (y0 + y1) / 2.0;
    let y_bottom = CH - 4.0;
    let _ = write!(
        s,
        r#"<text x="{cx_lab:.1}" y="{y_bottom:.1}" class="axis-label" text-anchor="middle">{x_label}</text><text x="14" y="{cy_lab:.1}" class="axis-label" text-anchor="middle" transform="rotate(-90 14 {cy_lab:.1})">{y_label}</text>"#
    );
    s
}

fn fmt_num(v: f64) -> String {
    if v.abs() >= 1000.0 {
        format!("{:.0}", v)
    } else if v.abs() >= 10.0 {
        format!("{:.0}", v)
    } else if v.abs() >= 1.0 {
        format!("{:.1}", v)
    } else {
        format!("{:.2}", v)
    }
}

fn line_chart(series: &[Series], y_max: f64, x_label: &str, y_label: &str) -> String {
    let (pad_l, pad_r, pad_t, pad_b) = (66.0, 16.0, 14.0, 36.0);
    let (x0, y0) = (pad_l, CH - pad_b);
    let (x1, y1) = (CW - pad_r, pad_t);
    let x_max = series
        .iter()
        .filter_map(|s| s.points.last().map(|p| p.0))
        .fold(1.0_f64, f64::max);
    let sx = |x: f64| x0 + (x / x_max) * (x1 - x0);
    let sy = |y: f64| y0 - (y.min(y_max) / y_max) * (y0 - y1);

    let mut out = format!(r#"<svg viewBox="0 0 {CW} {CH}" class="chart" role="img">"#);
    out.push_str(&axes(x0, y0, x1, y1, x_max, y_max, x_label, y_label));
    for s in series {
        let pts = s
            .points
            .iter()
            .map(|(x, y)| format!("{:.1},{:.1}", sx(*x), sy(*y)))
            .collect::<Vec<_>>()
            .join(" ");
        let _ = write!(out, r#"<polyline class="ln" points="{pts}" stroke="{}"/>"#, s.colour);
    }
    // legend
    let mut lx = x0 + 8.0;
    for s in series {
        let _ = write!(
            out,
            r#"<line x1="{lx:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" class="ln" stroke="{}"/><text x="{:.1}" y="{:.1}" class="tick">{}</text>"#,
            pad_t + 8.0,
            lx + 14.0,
            pad_t + 8.0,
            s.colour,
            lx + 18.0,
            pad_t + 12.0,
            esc(s.name)
        );
        lx += 34.0 + 6.5 * s.name.chars().count() as f64;
    }
    out.push_str("</svg>");
    out
}

fn bar_chart(points: &[(f64, f64)], x_label: &str, y_label: &str, highlight: Option<(f64, &str)>) -> String {
    let (pad_l, pad_r, pad_t, pad_b) = (66.0, 16.0, 14.0, 36.0);
    let (x0, y0) = (pad_l, CH - pad_b);
    let (x1, y1) = (CW - pad_r, pad_t);
    let x_max = points.iter().map(|p| p.0).fold(1.0_f64, f64::max);
    let y_max = points.iter().map(|p| p.1).fold(1.0_f64, f64::max).max(1.0);
    let bw = if points.is_empty() { 1.0 } else { ((x1 - x0) / points.len() as f64 * 0.86).max(0.6) };

    let mut out = format!(r#"<svg viewBox="0 0 {CW} {CH}" class="chart" role="img">"#);
    out.push_str(&axes(x0, y0, x1, y1, x_max, y_max, x_label, y_label));
    if let Some((cut, label)) = highlight {
        let cx = x0 + (cut / x_max) * (x1 - x0);
        let _ = write!(
            out,
            r#"<rect x="{cx:.1}" y="{y1:.1}" width="{:.1}" height="{:.1}" class="vband"/><text x="{:.1}" y="{:.1}" class="tick">{}</text>"#,
            (x1 - cx).max(0.0),
            y0 - y1,
            cx + 4.0,
            y1 + 12.0,
            esc(label)
        );
    }
    for (x, y) in points {
        let cx = x0 + (x / x_max) * (x1 - x0);
        let h = (y / y_max) * (y0 - y1);
        let _ = write!(
            out,
            r#"<rect x="{:.1}" y="{:.1}" width="{bw:.1}" height="{:.1}" class="bar"/>"#,
            cx - bw / 2.0,
            y0 - h,
            h
        );
    }
    out.push_str("</svg>");
    out
}

/// No-JS tabs: hidden radio inputs + labels in a row, panels shown by
/// `:checked ~ #panel`. The per-tab rules are emitted inline, scoped by prefix.
fn tabs(prefix: &str, panels: &[(&str, String)]) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        r#"<style>.tabs-{prefix} input{{display:none}}.tabs-{prefix} .tabpanel{{display:none}}"#
    );
    for (i, _) in panels.iter().enumerate() {
        let _ = write!(s, r#".tabs-{prefix} #t{prefix}-{i}:checked ~ #p{prefix}-{i}{{display:block}}"#);
    }
    s.push_str("</style>");
    let _ = write!(s, r#"<div class="tabs tabs-{prefix}">"#);
    for (i, (label, _)) in panels.iter().enumerate() {
        let _ = write!(
            s,
            r#"<input type="radio" name="t{prefix}" id="t{prefix}-{i}"{}><label for="t{prefix}-{i}">{}</label>"#,
            if i == 0 { " checked" } else { "" },
            esc(label)
        );
    }
    for (_, content) in panels {
        let _ = write!(s, r#"<div class="tabpanel">{content}</div>"#);
    }
    s.push_str("</div>");
    s
}

/// The 256x256 heatmap as inline SVG, log-scaled so the peak does not wash out
/// the rest. Only non-empty cells are emitted.
fn heatmap_svg(acc: &Accum) -> String {
    let size = crate::heatmap::SIZE;
    let cell = 2.0;
    let w = size as f64 * cell;
    let peak = acc.heatmap.iter().copied().max().unwrap_or(1).max(1);
    let log_peak = (1.0 + peak as f64).ln();
    let mut s = format!(r#"<svg viewBox="0 0 {w} {w}" class="heatmap" role="img">"#);
    for y in 0..size {
        for x in 0..size {
            let c = acc.heatmap[y * size + x];
            if c == 0 {
                continue;
            }
            let t = (1.0 + c as f64).ln() / log_peak;
            let v = (t * 200.0) as i32;
            let _ = write!(
                s,
                r#"<rect x="{:.1}" y="{:.1}" width="{cell}" height="{cell}" fill="rgb({},{},{})"/>"#,
                x as f64 * cell,
                y as f64 * cell,
                240 - v,
                240 - v,
                255 - (v / 3)
            );
        }
    }
    s.push_str("</svg>");
    s
}

fn table(rows: &[(&str, String)]) -> String {
    let mut s = String::from("<table class=\"kv\">");
    for (k, v) in rows {
        let _ = write!(s, "<tr><th>{}</th><td>{v}</td></tr>", esc(k));
    }
    s.push_str("</table>");
    s
}

const STYLE: &str = r#"
:root { --ink:#1d2530; --muted:#6b7580; --rule:#dde4ea; --accent:#2f6f9f; }
* { box-sizing: border-box; }
body { margin:0; padding:28px 34px 60px; color:var(--ink); background:#fff;
       font:14px/1.5 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
h1 { font-size:23px; margin:0 0 2px; color:var(--accent); }
h2 { font-size:15px; margin:30px 0 8px; padding-bottom:4px; border-bottom:1px solid var(--rule); }
.sub { color:var(--muted); font-size:13px; margin-bottom:20px; }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:22px; }
.chart { width:100%; height:auto; background:#fbfcfd; border:1px solid var(--rule); border-radius:4px; }
.ax { stroke:#9aa4ad; stroke-width:1; }
.grid { stroke:#eef2f5; stroke-width:1; }
.tick { fill:var(--muted); font-size:11px; }
.tick.end { text-anchor:end; }
.axis-label { fill:var(--muted); font-size:11.5px; }
.ln { fill:none; stroke-width:1.8; }
.bar { fill:var(--accent); }
.vband { fill:#e8d6d6; opacity:.55; }
table.kv { border-collapse:collapse; width:100%; }
table.kv th { text-align:left; font-weight:600; color:var(--muted); padding:3px 12px 3px 0; white-space:nowrap; vertical-align:top; }
table.kv td { padding:3px 0; }
table.files { border-collapse:collapse; width:100%; font-size:13px; }
table.files th, table.files td { text-align:right; padding:4px 10px 4px 0; border-bottom:1px solid var(--rule); }
table.files th:first-child, table.files td:first-child { text-align:left; }
ul.obs { margin:6px 0 0 18px; padding:0; }
ul.obs li { margin:3px 0; }
.note { background:#f7f9fb; border-left:3px solid var(--accent); padding:9px 12px; font-size:13px; color:#333; }
code { font-family:ui-monospace, SFMono-Regular, Menlo, monospace; font-size:12.5px; }
footer { margin-top:34px; color:var(--muted); font-size:12px; border-top:1px solid var(--rule); padding-top:8px; }
.tabs { margin:8px 0 12px; }
.tabs input { display:none; }
.tabs label { display:inline-block; padding:5px 16px; margin-right:4px; border:1px solid var(--rule);
              border-bottom:none; border-radius:4px 4px 0 0; cursor:pointer; font-size:13px;
              color:var(--muted); background:#f7f9fb; }
.tabs input:checked + label { background:#fff; color:var(--accent); font-weight:600; border-color:var(--accent); }
.tabs .tabpanel { display:none; border:1px solid var(--rule); border-radius:0 4px 4px 4px; padding:12px; background:#fff; }
.heatmap { width:100%; height:auto; background:#fff; border:1px solid var(--rule); border-radius:4px; }
table.rep { border-collapse:collapse; width:100%; font-size:12.5px; }
table.rep th, table.rep td { text-align:right; padding:3px 10px 3px 0; border-bottom:1px solid var(--rule); }
table.rep th:first-child, table.rep td:first-child { text-align:left; }
table.rep code { font-size:11.5px; }
.callout { background:#f7f9fb; border-left:3px solid var(--accent); padding:9px 12px; font-size:13px; color:#333; margin:10px 0; }
"#;

/// Write the report. `insert` is present when `--insert-hist` ran.
pub fn write_html(
    path: &Path,
    files: &[FileReport],
    total: &Accum,
    offset: u8,
    elapsed_s: f64,
    insert: Option<(&InsertHist, u64)>,
) -> Result<()> {
    let n_pos = total.covered_positions();
    let gc = gc_fraction(total);
    let depth_at = |pos: usize| -> u64 { total.base_at(pos, 0) + total.base_at(pos, 1) + total.base_at(pos, 2) + total.base_at(pos, 3) + total.base_at(pos, 4) };

    // ---- per-base sequence content (FastQC's signature plot) -----------------
    let mut seq_series = Vec::new();
    for (code, name, colour) in [(0usize, "A", "#2e7d5b"), (1, "C", "#2f6f9f"), (2, "G", "#c98a1a"), (3, "T", "#b23a48")] {
        let mut pts = Vec::with_capacity(n_pos);
        for pos in 0..n_pos {
            let d = depth_at(pos);
            if d > 0 {
                pts.push(((pos + 1) as f64, 100.0 * total.base_at(pos, code) as f64 / d as f64));
            }
        }
        seq_series.push(Series { name, colour, points: pts });
    }
    let seq_chart = line_chart(&seq_series, 100.0, "position in read (bp)", "% of bases at this position");

    // ---- per-base quality, with the 10th and 90th percentiles ----------------
    let mut q_pts = Vec::with_capacity(n_pos);
    let mut q_lo = Vec::with_capacity(n_pos);
    let mut q_hi = Vec::with_capacity(n_pos);
    for pos in 0..n_pos {
        let row = total.qual_at(pos);
        if row.iter().sum::<u64>() == 0 {
            continue;
        }
        let (mean, q10, _q50, q90) = pos_stats(&row, offset);
        q_pts.push(((pos + 1) as f64, mean));
        q_lo.push(((pos + 1) as f64, q10));
        q_hi.push(((pos + 1) as f64, q90));
    }
    let q_chart = line_chart(
        &[Series { name: "p90", colour: "#a9c4d9", points: q_hi },
          Series { name: "mean", colour: "#2f6f9f", points: q_pts },
          Series { name: "p10", colour: "#a9c4d9", points: q_lo }],
        40.0, "position in read (bp)", "Phred quality",
    );

    // ---- per-base N content --------------------------------------------------
    let mut n_pts = Vec::with_capacity(n_pos);
    for pos in 0..n_pos {
        let d = depth_at(pos);
        if d > 0 {
            n_pts.push(((pos + 1) as f64, 100.0 * total.base_at(pos, 4) as f64 / d as f64));
        }
    }
    let n_chart = line_chart(&[Series { name: "N", colour: "#6b7580", points: n_pts }],
                             5.0, "position in read (bp)", "% N at this position");

    // ---- distributions: length, per-read mean quality, per-read GC -----------
    let len_bars: Vec<(f64, f64)> = total.len_hist.iter().enumerate()
        .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
    let len_chart = bar_chart(&len_bars, "read length (bp)", "reads", None);

    let mq_bars: Vec<(f64, f64)> = total.meanq_raw_hist.iter().enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| (i as f64 - offset as f64, c as f64)).collect();
    let mq_chart = bar_chart(&mq_bars, "mean Phred quality of a read", "reads", None);

    let gc_bars: Vec<(f64, f64)> = total.gc_hist.iter().enumerate()
        .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
    let gc_chart = bar_chart(&gc_bars, "G+C bases in the read", "reads", None);

    // ---- insert sizes, when they were measured ------------------------------
    let (insert_chart, insert_table) = match insert {
        Some((ih, sampled)) => {
            let mut merged: std::collections::BTreeMap<u64, u64> = std::collections::BTreeMap::new();
            for (k, &v) in ih.bins.iter().enumerate() {
                if v > 0 {
                    *merged.entry(k as u64).or_insert(0) += v;
                }
            }
            for (k, &v) in ih.bins_adapter.iter().enumerate() {
                if v > 0 {
                    *merged.entry(k as u64).or_insert(0) += v;
                }
            }
            let bars: Vec<(f64, f64)> = merged.into_iter().map(|(k, v)| (k as f64, v as f64)).collect();
            let chart = bar_chart(&bars, "inferred insert length (bp)", "pairs", None);
            let sampled = ih.pairs_sampled.max(sampled);
            let pct = |n: u64| format!("{n} ({:.1}%)", 100.0 * n as f64 / sampled.max(1) as f64);
            let t = table(&[
                ("pairs sampled", format!("{sampled}")),
                ("measured from overlap", pct(ih.measured)),
                ("measured from the adapter", pct(ih.measured_adapter)),
                ("no overlap (insert above the resolution)", pct(ih.censored_above)),
                ("complete overlap, no adapter found", pct(ih.censored_below)),
                ("pairs with unequal mate lengths", pct(ih.unequal_len)),
                ("mean / sd", format!("{:.1} / {:.1}", ih.mean_all(), ih.sd_all())),
                ("deciles", ih.deciles().iter()
                    .map(|(l, v)| format!("{l} {v:.0}")).collect::<Vec<_>>().join("  ")),
            ]);
            (chart, t)
        }
        None => (String::new(), String::new()),
    };

    // ---- polyX: 4 base tabs, each with the 3' and 5' run-length distributions --
    let base_names = ["A", "C", "G", "T"];
    let mut polyx_panels: Vec<(&str, String)> = Vec::new();
    for (bi, name) in base_names.iter().enumerate() {
        let b3: Vec<(f64, f64)> = total.polyx3[bi].iter().enumerate()
            .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
        let b5: Vec<(f64, f64)> = total.polyx5[bi].iter().enumerate()
            .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
        let c3 = bar_chart(&b3, "3' run length (bp)", "reads", None);
        let c5 = bar_chart(&b5, "5' run length (bp)", "reads", None);
        polyx_panels.push((name, format!(r#"<div class="grid2"><div>{c3}</div><div>{c5}</div></div>"#)));
    }
    let (hg, hd) = crate::rich::head_dropout(total);
    let polyx_section = format!(
        "{}{}",
        tabs("polyx", &polyx_panels),
        format!(
            r#"<div class="callout">First base is G in {:.2}% of reads against {:.2}% for the library
                (head dropout {:.2} points) — a positive dropout is the 2-colour first-cycle G
                signature. The 5' tab for G shows the run-length distribution of that head.</div>"#,
            hg * 100.0, (hg + hd) * 100.0, hd * 100.0
        )
    );

    // ---- repeated sequences: full sequences + 12-mer spectrum -----------------
    let reads = total.reads.max(1) as f64;
    let min_count = (reads * 0.001).max(10.0) as u32;
    let seqs = total.seq_counts.report(min_count);
    let mut rep_rows = String::new();
    for (s, c) in seqs.iter().take(20) {
        let _ = write!(
            rep_rows,
            "<tr><td><code>{}</code></td><td>{}</td><td>{:.3}%</td></tr>",
            esc(&String::from_utf8_lossy(s)), c, 100.0 * *c as f64 / reads
        );
    }
    let rep_table = if seqs.is_empty() {
        "<p class=\"note\">No sequence appears more than 0.1% of reads — the library is not dominated by a repeated sequence.</p>".to_string()
    } else {
        format!(
            "<table class=\"rep\"><tr><th>sequence</th><th>count</th><th>% of reads</th></tr>{rep_rows}</table>"
        )
    };
    let kmers = total.kmer_counts.report(min_count);
    let kmer_bars: Vec<(f64, f64)> = kmers.iter().take(20)
        .map(|(k, c)| (crate::repeated::unpack_kmer(*k, 12).len() as f64, *c as f64)).collect();
    let kmer_chart = bar_chart(&kmer_bars, "top 12-mer (by length of its name)", "count", None);
    let dup = 1.0 - total.distinct_seqs as f64 / reads;
    let repeats_section = format!(
        r#"<div class="grid2"><div><h3>Over-represented full sequences</h3>{rep_table}</div>
           <div><h3>Top 12-mers</h3>{kmer_chart}</div></div>
           <div class="callout">{} distinct sequences of {reads:.0} reads — duplication level {:.2}%.</div>"#,
        total.distinct_seqs, dup * 100.0
    );

    // ---- heatmap: full view (side-by-side arrives with --processed-*) ---------
    let heatmap_section = format!(
        r#"<div class="callout">Each read is one cell at (T+C, T+G). A tight cluster is a uniform
           composition; a spread or a second cluster is a sub-population (e.g. adapter-dimers or a
           contaminant). Peak at ({}, {}).</div>{}"#,
        crate::rich::heatmap_metrics(total).0.get("heatmap_peak_x").map(|v| v.to_string()).unwrap_or_default(),
        crate::rich::heatmap_metrics(total).0.get("heatmap_peak_y").map(|v| v.to_string()).unwrap_or_default(),
        heatmap_svg(total)
    );

    // ---- observations: numbers, not verdicts --------------------------------
    let mut obs: Vec<String> = Vec::new();
    if n_pos > 0 {
        let d0 = depth_at(0);
        if d0 > 0 {
            let g0 = 100.0 * total.base_at(0, 2) as f64 / d0 as f64;
            let g_all = 100.0 * (0..n_pos).map(|p| total.base_at(p, 2)).sum::<u64>() as f64
                / (0..n_pos).map(depth_at).sum::<u64>().max(1) as f64;
            obs.push(format!(
                "first base is G in {g0:.2}% of reads against {g_all:.2}% for the library as a whole — \
                 a first-cycle G dropout shows up as the first being far higher"));
        }
        if n_pos >= 2 {
            let d_last = depth_at(n_pos - 1);
            if d_last > 0 {
                let g_last = 100.0 * total.base_at(n_pos - 1, 2) as f64 / d_last as f64;
                let d_mid = depth_at(n_pos / 2);
                let g_mid = if d_mid > 0 { 100.0 * total.base_at(n_pos / 2, 2) as f64 / d_mid as f64 } else { 0.0 };
                obs.push(format!(
                    "last base is G in {g_last:.2}% against {g_mid:.2}% mid-read — a polyG tail is the usual cause \
                     when the end is far higher"));
            }
        }
        let max_n = (0..n_pos).map(|p| {
            let d = depth_at(p);
            if d > 0 { 100.0 * total.base_at(p, 4) as f64 / d as f64 } else { 0.0 }
        }).fold(0.0_f64, f64::max);
        obs.push(format!("highest per-position N content: {max_n:.3}%"));
    }
    obs.push(format!(
        "non-ACGT bases: {} of {} ({:.4}%); %GC denominators exclude them",
        bases_other(total), total.bases,
        100.0 * bases_other(total) as f64 / total.bases.max(1) as f64
    ));
    if total.reads_over_max_len > 0 {
        obs.push(format!(
            "{} reads exceed --max-read-len and are counted in the summary but excluded from the \
             per-position charts", total.reads_over_max_len));
    }
    obs.push("every count in this report is exact: no sampling, no read was skipped".to_string());

    // ---- per-file table ------------------------------------------------------
    let mut file_rows = String::new();
    for f in files {
        let _ = write!(
            file_rows,
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.3}</td><td>{:.2}</td><td>{}</td></tr>",
            esc(&f.path.file_name().unwrap_or_default().to_string_lossy()),
            f.acc.reads, f.acc.bases, gc_fraction(&f.acc), mean_read_mean_q(&f.acc, offset),
            esc(&f.scan.decompress)
        );
    }

    let mut h = String::with_capacity(64 * 1024);
    h.push_str("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">");
    let _ = write!(h, "<title>fastcheck report — {} file(s)</title><style>{STYLE}</style></head><body>",
                   files.len());
    let _ = write!(h, "<h1>fastcheck report</h1><div class=\"sub\">{} file(s) · {} reads · {} bases · \
                       {:.2} s · Phred+{offset} · exact counts, no sampling</div>",
                   files.len(), total.reads, total.bases, elapsed_s);

    h.push_str("<h2>Summary</h2><div class=\"grid2\"><div>");
    h.push_str(&table(&[
        ("reads", format!("{}", total.reads)),
        ("bases", format!("{}", total.bases)),
        ("%GC", format!("{:.2}%", gc * 100.0)),
        ("mean read quality", format!("{:.2}", mean_read_mean_q(total, offset))),
        ("read length", format!("{}–{} (median {}, mean {:.1})",
            total.len_hist.iter().position(|&c| c > 0).unwrap_or(0),
            total.len_hist.len().saturating_sub(1),
            quantile(&total.len_hist, 0.5), mean_len(total))),
        ("reads containing N", format!("{}", total.reads_with_n)),
    ]));
    h.push_str("</div><div>");
    h.push_str(&table(&[
        ("files", format!("{}", files.len())),
        ("elapsed", format!("{elapsed_s:.2} s")),
        ("decompression", esc(&files.iter().map(|f| f.scan.decompress).collect::<std::collections::BTreeSet<_>>()
            .into_iter().collect::<Vec<_>>().join(", "))),
        ("insert sizes", if insert.is_some() { "measured".into() } else { "not measured (use --insert-hist)".to_string() }),
    ]));
    h.push_str("</div></div>");

    h.push_str("<table class=\"files\"><tr><th>file</th><th>reads</th><th>bases</th><th>%GC</th><th>mean Q</th><th>decompression</th></tr>");
    h.push_str(&file_rows);
    h.push_str("</table>");

    h.push_str("<h2>Per-base sequence content</h2>");
    let _ = write!(h, "<div class=\"grid2\"><div>{seq_chart}</div><div>{q_chart}</div></div>");
    h.push_str("<div class=\"note\">Left: the base composition at each position of the read. Right: \
                quality by position, with the 10th and 90th percentiles shaded by the same colour — \
                a widening spread is quality falling off unevenly between reads, not just a lower mean.</div>");

    h.push_str("<h2>N content by position</h2>");
    let _ = write!(h, "{n_chart}");

    h.push_str("<h2>Distributions</h2>");
    let _ = write!(h, "<div class=\"grid2\"><div>{len_chart}</div><div>{mq_chart}</div></div>");
    let _ = write!(h, "<div class=\"grid2\"><div>{gc_chart}</div><div>{}</div></div>",
                   if insert.is_some() { insert_chart.clone() } else { String::new() });

    h.push_str("<h2>polyX run lengths</h2>");
    h.push_str(&polyx_section);

    h.push_str("<h2>Repeated sequences</h2>");
    h.push_str(&repeats_section);

    h.push_str("<h2>Heatmap</h2>");
    h.push_str(&heatmap_section);

    if insert.is_some() {
        h.push_str("<h2>Insert sizes</h2>");
        h.push_str(&insert_table);
        let _ = write!(h, "<div class=\"grid2\"><div>{insert_chart}</div><div></div></div>");
    }

    h.push_str("<h2>Observations</h2><ul class=\"obs\">");
    for o in &obs {
        let _ = write!(h, "<li>{o}</li>");
    }
    h.push_str("</ul>");
    let _ = write!(h, "<footer>fastcheck</footer>");
    h.push_str("</body></html>\n");

    let mut fh = std::fs::File::create(path).with_context(|| format!("create {}", path.display()))?;
    std::io::Write::write_all(&mut fh, h.as_bytes())?;
    Ok(())
}
