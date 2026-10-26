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

use anyhow::Result;

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

fn line_chart(series: &[Series], y_max: f64, x_max: f64, x_label: &str, y_label: &str) -> String {
    let (pad_l, pad_r, pad_t, pad_b) = (66.0, 16.0, 14.0, 36.0);
    let (x0, y0) = (pad_l, CH - pad_b);
    let (x1, y1) = (CW - pad_r, pad_t);
    let x_max = x_max.max(1.0);
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

/// A bar chart with a text label under each bar (for k-mers and other
/// categorical x-axes). Labels are rotated 45° and truncated to fit.
fn labeled_bar_chart(bars: &[(String, f64)], y_label: &str) -> String {
    let (pad_l, pad_r, pad_t, pad_b) = (66.0, 16.0, 14.0, 46.0);
    let (x0, y0) = (pad_l, CH - pad_b);
    let (x1, y1) = (CW - pad_r, pad_t);
    let n = bars.len().max(1);
    let y_max = bars.iter().map(|(_, v)| *v).fold(1.0_f64, f64::max).max(1.0);
    let bw = ((x1 - x0) / n as f64 * 0.7).max(1.0);
    let mut out = format!(r#"<svg viewBox="0 0 {CW} {CH}" class="chart" role="img">"#);
    out.push_str(&axes(x0, y0, x1, y1, n as f64, y_max, "", y_label));
    for (i, (label, v)) in bars.iter().enumerate() {
        let cx = x0 + (i as f64 + 0.5) * (x1 - x0) / n as f64;
        let h = (v / y_max) * (y0 - y1);
        let _ = write!(
            out,
            r#"<rect x="{:.1}" y="{:.1}" width="{bw:.1}" height="{:.1}" class="bar"/>"#,
            cx - bw / 2.0,
            y0 - h,
            h
        );
        // label under the bar, rotated and truncated
        let short: String = label.chars().take(8).collect();
        let _ = write!(
            out,
            r#"<text x="{cx:.1}" y="{:.1}" class="tick" text-anchor="end" transform="rotate(-45 {cx:.1} {:.1})">{}</text>"#,
            y0 + 14.0,
            y0 + 14.0,
            esc(&short)
        );
    }
    out.push_str("</svg>");
    out
}

/// Tabs: buttons in a row, panels shown by a small JS handler (the same
/// mechanism MultiQC uses). The first tab is active by default.
fn tabs(prefix: &str, panels: &[(&str, String)]) -> String {
    let mut s = String::new();
    let _ = write!(s, r#"<div class="tabs" id="tabs-{prefix}" role="tablist" aria-label="{prefix}" aria-orientation="horizontal">"#);
    for (i, (label, _)) in panels.iter().enumerate() {
        let sel = i == 0;
        let _ = write!(
            s,
            r#"<button type="button" role="tab" id="tab-{prefix}-{i}" class="tab{}" aria-selected="{}" aria-controls="p{prefix}-{i}" tabindex="{}">{}</button>"#,
            if sel { " active" } else { "" },
            if sel { "true" } else { "false" },
            if sel { 0 } else { -1 },
            esc(label)
        );
    }
    for (i, (_, content)) in panels.iter().enumerate() {
        let _ = write!(
            s,
            r#"<div role="tabpanel" id="p{prefix}-{i}" class="tabpanel{}" aria-labelledby="tab-{prefix}-{i}">{content}</div>"#,
            if i == 0 { " active" } else { "" }
        );
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
            // A viridis-like palette (dark purple -> teal -> yellow), so even
            // sparse cells are clearly visible against the white background.
            const PALETTE: [(u8, u8, u8); 8] = [
                (68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98),
                (253, 231, 37), (240, 180, 20), (220, 120, 10), (180, 60, 5),
            ];
            let idx = (t * (PALETTE.len() - 1) as f64).round() as usize;
            let (r, g, b) = PALETTE[idx.min(PALETTE.len() - 1)];
            let _ = write!(
                s,
                r#"<rect x="{:.1}" y="{:.1}" width="{cell}" height="{cell}" fill="rgb({r},{g},{b})"/>"#,
                x as f64 * cell,
                y as f64 * cell
            );
        }
    }
    s.push_str("</svg>");
    s
}

/// A titled card around a chart, so sections read as labelled panels rather
/// than bare SVGs floating in a grid.
fn chart_card(title: &str, chart: &str) -> String {
    format!(
        r#"<div class="card"><div class="card-title">{}</div>{chart}</div>"#,
        esc(title)
    )
}

/// A titled card holding a note instead of a chart.
fn chart_card_note(title: &str, note: &str) -> String {
    format!(
        r#"<div class="card"><div class="card-title">{}</div><div class="note">{}</div></div>"#,
        esc(title),
        esc(note)
    )
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
       font:15px/1.5 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
h1 { font-size:25px; margin:0 0 2px; color:var(--accent); }
h2 { font-size:17px; margin:30px 0 8px; padding-bottom:4px; border-bottom:1px solid var(--rule); }
.sub { color:var(--muted); font-size:14px; margin-bottom:20px; }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:22px; }
.chart { width:100%; height:auto; background:#fbfcfd; border:1px solid var(--rule); border-radius:4px; }
.ax { stroke:#9aa4ad; stroke-width:1; }
.grid { stroke:#eef2f5; stroke-width:1; }
.tick { fill:#4a5560; font-size:12px; }
.tick.end { text-anchor:end; }
.axis-label { fill:#4a5560; font-size:12.5px; }
.ln { fill:none; stroke-width:1.8; }
.bar { fill:var(--accent); }
.vband { fill:#e8d6d6; opacity:.55; }
table.kv { border-collapse:collapse; width:100%; }
table.kv th { text-align:left; font-weight:600; color:var(--muted); padding:3px 12px 3px 0; white-space:nowrap; vertical-align:top; }
table.kv td { padding:3px 0; }
table.files { border-collapse:collapse; width:100%; font-size:14px; }
table.files th, table.files td { text-align:right; padding:4px 10px 4px 0; border-bottom:1px solid var(--rule); }
table.files th:first-child, table.files td:first-child { text-align:left; }
ul.obs { margin:6px 0 0 18px; padding:0; }
ul.obs li { margin:3px 0; }
.note { background:#f7f9fb; border-left:3px solid var(--accent); padding:9px 12px; font-size:13px; color:#333; }
code { font-family:ui-monospace, SFMono-Regular, Menlo, monospace; font-size:12.5px; }
footer { margin-top:34px; color:var(--muted); font-size:12px; border-top:1px solid var(--rule); padding-top:8px; }
.tabs { margin:8px 0 12px; }
.tabs .tab { display:inline-block; padding:5px 16px; margin-right:4px; border:1px solid var(--rule);
             border-bottom:none; border-radius:4px 4px 0 0; cursor:pointer; font-size:13px;
             color:var(--muted); background:#f7f9fb; font-family:inherit; }
.tabs .tab.active { background:#fff; color:var(--accent); font-weight:600; border-color:var(--accent); }
.tabs .tabpanel { display:none; border:1px solid var(--rule); border-radius:0 4px 4px 4px; padding:12px; background:#fff; }
.tabs .tabpanel.active { display:block; }
.heatmap { width:100%; height:auto; background:#fff; border:1px solid var(--rule); border-radius:4px; }
table.rep { border-collapse:collapse; width:100%; font-size:12.5px; }
table.rep th, table.rep td { text-align:right; padding:3px 10px 3px 0; border-bottom:1px solid var(--rule); }
table.rep th:first-child, table.rep td:first-child { text-align:left; }
table.rep code { font-size:11.5px; }
.callout { background:#f7f9fb; border-left:3px solid var(--accent); padding:9px 12px; font-size:13px; color:#333; margin:10px 0; }
nav.nav { position:sticky; top:0; background:#fff; border-bottom:1px solid var(--rule); padding:8px 0; margin-bottom:18px; z-index:10; }
nav.nav a { display:inline-block; margin-right:14px; color:var(--accent); text-decoration:none; font-size:13px; font-weight:600; }
nav.nav a:hover { text-decoration:underline; }
.ba-label { font-size:11px; font-weight:700; color:var(--muted); text-transform:uppercase; letter-spacing:0.5px; margin-bottom:4px; }
h3 { font-size:13.5px; margin:16px 0 6px; color:var(--ink); }
.card { border:1px solid var(--rule); border-radius:6px; padding:10px 12px 4px; background:#fff; }
.card-title { font-size:12.5px; font-weight:600; color:var(--ink); margin-bottom:2px; }
.grid2 .card { height:100%; }
.skip-link { position:absolute; left:-9999px; top:0; z-index:100; background:#fff; color:var(--accent); padding:8px 12px; }
.skip-link:focus { left:8px; top:8px; outline:2px solid var(--accent); }
html { scroll-behavior:smooth; }
@media (prefers-reduced-motion: reduce){ html { scroll-behavior:auto; } }
details.module { border:1px solid var(--rule); border-radius:6px; margin:0 0 18px; background:#fff; }
details.module > summary { cursor:pointer; list-style:none; padding:10px 14px; user-select:none; }
details.module > summary::-webkit-details-marker { display:none; }
details.module > summary h2 { display:inline-block; margin:0; border:none; padding:0; scroll-margin-top:60px; }
details.module > summary:hover { background:#f7f9fb; }
details.module[open] > summary { border-bottom:1px solid var(--rule); }
details.module > .body { padding:14px; }
nav.nav a { padding:4px 2px; }
nav.nav a[aria-current="true"] { color:var(--ink); border-bottom:2px solid var(--accent); }
button:focus-visible, a:focus-visible, summary:focus-visible { outline:2px solid #1a5fb4; outline-offset:2px; }
table.files caption { text-align:left; font-weight:600; color:var(--muted); font-size:13px; margin-bottom:4px; }
table.files thead th { border-bottom:2px solid var(--rule); }
table.files tbody tr:hover { background:#fff8e1; }
.banner { border:1px solid #e0b400; background:#fff8e1; border-radius:6px; padding:10px 12px; margin:0 0 18px; font-size:14px; }
"#;

/// Write the report. `insert` is present when `--insert-hist` ran, and `partial`
/// carries the reason a budget stopped the scan (`None` on a complete run), the
/// same value the JSON report publishes as `partial.reason`. `overrep_min_count` is
/// the resolved over-representation threshold, so the repeated-sequence section uses
/// the same one the `overrep` block does.
pub fn write_html(
    path: &Path,
    files: &[FileReport],
    total: &Accum,
    offset: u8,
    elapsed_s: f64,
    insert: Option<&InsertHist>,
    partial: Option<&str>,
    processed: Option<(&[FileReport], &Accum)>,
    overrep_min_count: u32,
) -> Result<()> {
    let n_pos = total.covered_positions();
    let gc = gc_fraction(total);
    let depth_at = |pos: usize| -> u64 { total.base_at(pos, 0) + total.base_at(pos, 1) + total.base_at(pos, 2) + total.base_at(pos, 3) + total.base_at(pos, 4) };

    // Per-mate accumulators (R1, R2) and the processed side, when present.
    let r1 = files.first().map(|f| &f.acc);
    let r2 = files.get(1).map(|f| &f.acc);
    let (p_files, p_total) = match processed {
        Some((pf, pt)) => (Some(pf), Some(pt)),
        None => (None, None),
    };
    let pr1 = p_files.and_then(|pf| pf.first()).map(|f| &f.acc);
    let pr2 = p_files.and_then(|pf| pf.get(1)).map(|f| &f.acc);

    // Per-position series for one mate: (seq A/C/G/T, N, quality mean/lo/hi).
    fn pos_series(acc: &Accum, n_pos: usize, offset: u8) -> (Vec<Series<'static>>, Vec<(f64, f64)>, Vec<(f64, f64)>, Vec<(f64, f64)>, Vec<(f64, f64)>) {
        let depth_at = |pos: usize| -> u64 {
            acc.base_at(pos, 0) + acc.base_at(pos, 1) + acc.base_at(pos, 2) + acc.base_at(pos, 3) + acc.base_at(pos, 4)
        };
        let mut seq = Vec::new();
        for (code, name, colour) in [(0usize, "A", "#0072b2"), (1, "C", "#e69f00"), (2, "G", "#009e73"), (3, "T", "#d55e00")] {
            let mut pts = Vec::new();
            for pos in 0..n_pos {
                let d = depth_at(pos);
                if d > 0 {
                    pts.push(((pos + 1) as f64, 100.0 * acc.base_at(pos, code) as f64 / d as f64));
                }
            }
            seq.push(Series { name, colour, points: pts });
        }
        let mut n_pts = Vec::new();
        for pos in 0..n_pos {
            let d = depth_at(pos);
            if d > 0 {
                n_pts.push(((pos + 1) as f64, 100.0 * acc.base_at(pos, 4) as f64 / d as f64));
            }
        }
        let (mut q_pts, mut q_lo, mut q_hi) = (Vec::new(), Vec::new(), Vec::new());
        for pos in 0..n_pos {
            let row = acc.qual_at(pos);
            if row.iter().sum::<u64>() == 0 {
                continue;
            }
            let (mean, q10, _q50, q90) = pos_stats(&row, offset);
            q_pts.push(((pos + 1) as f64, mean));
            q_lo.push(((pos + 1) as f64, q10));
            q_hi.push(((pos + 1) as f64, q90));
        }
        (seq, n_pts, q_pts, q_lo, q_hi)
    }

    // Render a per-position plot with R1 / R2 / both tabs (SE: no tabs).
    fn mate_plot(
        prefix: &str,
        r1: Option<&Accum>,
        r2: Option<&Accum>,
        both: &Accum,
        n_pos: usize,
        offset: u8,
        y_max: f64,
        x_label: &str,
        y_label: &str,
        kind: u8, // 0 = seq, 1 = N, 2 = quality
    ) -> String {
        let chart = |acc: &Accum| -> String {
            let (seq, n_pts, q_pts, q_lo, q_hi) = pos_series(acc, n_pos, offset);
            match kind {
                0 => line_chart(&seq, y_max, n_pos as f64, x_label, y_label),
                1 => {
                    let n_max = n_pts.iter().map(|p| p.1).fold(0.0_f64, f64::max);
                    let scale = (n_max * 2.0).max(5.0);
                    line_chart(&[Series { name: "N", colour: "#6b7580", points: n_pts }], scale, n_pos as f64, x_label, y_label)
                }
                _ => line_chart(
                    &[Series { name: "p90", colour: "#a9c4d9", points: q_hi },
                      Series { name: "mean", colour: "#2f6f9f", points: q_pts },
                      Series { name: "p10", colour: "#a9c4d9", points: q_lo }],
                    y_max, n_pos as f64, x_label, y_label,
                ),
            }
        };
        match (r1, r2) {
            (Some(a), Some(b)) => tabs(
                prefix,
                &[("R1", chart(a)), ("R2", chart(b)), ("both", chart(both))],
            ),
            (Some(a), None) => chart(a),
            _ => chart(both),
        }
    }

    // ---- per-position plots (R1 / R2 / both), nested in before/after ---------
    let has_ba = p_total.is_some();
    let pn = p_total.map(|pt| pt.covered_positions());
    // mate_plot returns R1/R2/both tabs; when processed data is present we wrap
    // that whole tab group in a before/after toggle so each graph appears once.
    let seq_chart = if has_ba {
        let (pt, pn) = (p_total.unwrap(), pn.unwrap());
        tabs("seqba", &[
            ("before", mate_plot("seqb", r1, r2, total, n_pos, offset, 100.0, "position in read (bp)", "% of bases at this position", 0)),
            ("after", mate_plot("seqa", pr1, pr2, pt, pn, offset, 100.0, "position in read (bp)", "% of bases at this position", 0)),
        ])
    } else {
        mate_plot("seq", r1, r2, total, n_pos, offset, 100.0, "position in read (bp)", "% of bases at this position", 0)
    };
    let n_chart = if has_ba {
        let (pt, pn) = (p_total.unwrap(), pn.unwrap());
        tabs("nba", &[
            ("before", mate_plot("nb", r1, r2, total, n_pos, offset, 5.0, "position in read (bp)", "% N at this position", 1)),
            ("after", mate_plot("na", pr1, pr2, pt, pn, offset, 5.0, "position in read (bp)", "% N at this position", 1)),
        ])
    } else {
        mate_plot("n", r1, r2, total, n_pos, offset, 5.0, "position in read (bp)", "% N at this position", 1)
    };
    let q_chart = if has_ba {
        let (pt, pn) = (p_total.unwrap(), pn.unwrap());
        tabs("qba", &[
            ("before", mate_plot("qb", r1, r2, total, n_pos, offset, 40.0, "position in read (bp)", "Phred quality", 2)),
            ("after", mate_plot("qa", pr1, pr2, pt, pn, offset, 40.0, "position in read (bp)", "Phred quality", 2)),
        ])
    } else {
        mate_plot("q", r1, r2, total, n_pos, offset, 40.0, "position in read (bp)", "Phred quality", 2)
    };

    // ---- distributions (read length / mean Q / GC), before/after ------------
    fn dist_for(acc: &Accum, offset: u8, note_insert: String) -> String {
        let lens: Vec<(f64, f64)> = acc.len_hist.iter().enumerate()
            .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
        let mqs: Vec<(f64, f64)> = acc.meanq_raw_hist.iter().enumerate()
            .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64 - offset as f64, c as f64)).collect();
        let gcs: Vec<(f64, f64)> = acc.gc_hist.iter().enumerate()
            .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
        format!(
            r#"<div class="grid2"><div>{}</div><div>{}</div></div><div class="grid2"><div>{}</div><div>{}</div></div>"#,
            chart_card("Read length", &bar_chart(&lens, "read length (bp)", "reads", None)),
            chart_card("Mean read quality", &bar_chart(&mqs, "mean Phred quality of a read", "reads", None)),
            chart_card("G+C content", &bar_chart(&gcs, "G+C bases in the read", "reads", None)),
            chart_card_note("Insert size", &note_insert),
        )
    }
    let insert_note = if insert.is_some() {
        "measured — see the Insert sizes section below.".to_string()
    } else {
        "not measured — run with --insert-hist".to_string()
    };
    let dist_before = dist_for(total, offset, insert_note.clone());
    let dist_section = if let Some(pt) = p_total {
        tabs("distba", &[
            ("before", dist_before.into()),
            ("after", dist_for(pt, offset, insert_note).into()),
        ])
    } else {
        dist_before
    };

    // ---- insert sizes, when they were measured ------------------------------
    let (insert_chart, insert_table) = match insert {
        Some(ih) => {
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
            // The denominator is the sample the pass actually read, never the
            // --insert-pairs cap: dividing by the cap printed "pairs sampled
            // 20000" for a 2000-pair run, with every bucket percentage 10x too
            // small (2000 (10.0%)), so the outcomes did not partition the sample.
            let sampled = ih.pairs_sampled;
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

    // ---- polyX: 4 base tabs, each with 3'/5' run-length distributions --------
    fn polyx_for(acc: &Accum, inner_prefix: &str) -> String {
        let mut panels: Vec<(&str, String)> = Vec::new();
        for (bi, name) in ["A", "C", "G", "T"].iter().enumerate() {
            let b3: Vec<(f64, f64)> = acc.polyx3[bi].iter().enumerate()
                .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
            let b5: Vec<(f64, f64)> = acc.polyx5[bi].iter().enumerate()
                .filter(|(_, &c)| c > 0).map(|(i, &c)| (i as f64, c as f64)).collect();
            panels.push((name, format!(r#"<div class="grid2"><div>{}</div><div>{}</div></div>"#,
                bar_chart(&b3, "3' run length (bp)", "reads", None), bar_chart(&b5, "5' run length (bp)", "reads", None))));
        }
        tabs(inner_prefix, &panels)
    }
    let polyx_before = polyx_for(total, "polyxb");
    let (hg, hd) = crate::rich::head_dropout(total);
    let polyx_callout = format!(
        r#"<div class="callout">First base is G in {:.2}% of reads against {:.2}% for the library
            (head dropout {:.2} points) — a positive dropout is the 2-colour first-cycle G
            signature. The 5' tab for G shows the run-length distribution of that head.</div>"#,
        hg * 100.0, (hg + hd) * 100.0, hd * 100.0
    );
    let polyx_section = if let Some(pt) = p_total {
        tabs("polyxba", &[
            ("before", polyx_before),
            ("after", polyx_for(pt, "polyxa")),
        ]) + &polyx_callout
    } else {
        polyx_before + &polyx_callout
    };

    // ---- repeated sequences: full sequences + 12-mer spectrum ----------------
    fn repeats_for(acc: &Accum, overrep_min_count: u32) -> String {
        let reads = acc.reads.max(1) as f64;
        // One threshold for over-represented SEQUENCES: the flag the caller resolved.
        // The k-mer spectrum keeps its own 0.1%-of-reads floor, exactly as the JSON
        // `kmer_spectrum` block computes it.
        let seqs = acc.seq_counts.at_least(overrep_min_count);
        let mut rep_rows = String::new();
        for (fp, c) in seqs.iter().take(20) {
            let _ = write!(
                rep_rows,
                "<tr><td><code>{fp:016x}</code></td><td>{}</td><td>{:.3}%</td></tr>",
                c,
                100.0 * *c as f64 / reads
            );
        }
        let rep_table = if seqs.is_empty() {
            format!(
                "<p class=\"note\">No sequence appears at least {overrep_min_count} times (the \
                 over-representation threshold) — the library is not dominated by a repeated \
                 sequence.</p>"
            )
        } else {
            format!(
                "<table class=\"rep\"><tr><th>fingerprint</th><th>count</th><th>% of reads</th></tr>{rep_rows}</table>"
            )
        };
        let min_count = (reads * 0.001).max(10.0) as u32;
        let kmers = acc.kmer_counts.report(min_count);
        let kmer_bars: Vec<(String, f64)> = kmers.iter().take(20)
            .map(|(k, c)| (crate::repeated::unpack_kmer(*k, 12), *c as f64)).collect();
        let kmer_chart = if kmer_bars.is_empty() {
            "<p class=\"note\">No 12-mer appears more than 0.1% of reads.</p>".to_string()
        } else {
            labeled_bar_chart(&kmer_bars, "count")
        };
        let dup = 1.0 - acc.seq_counts.distinct() as f64 / reads;
        format!(
            r#"<div class="grid2"><div><h3>Over-represented full sequences</h3>{rep_table}</div>
               <div><h3>Top 12-mers</h3>{kmer_chart}</div></div>
               <div class="callout">{} distinct sequences of {reads:.0} reads — duplication level {:.2}%.</div>"#,
            acc.seq_counts.distinct(), dup * 100.0
        )
    }
    let repeats_before = repeats_for(total, overrep_min_count);
    let repeats_section = if let Some(pt) = p_total {
        tabs("repba", &[
            ("before", repeats_before),
            ("after", repeats_for(pt, overrep_min_count)),
        ])
    } else {
        repeats_before
    };

    // ---- heatmap -------------------------------------------------------------
    let (hm_m, _) = crate::rich::heatmap_metrics(total);
    let hm_callout = format!(
        r#"<div class="callout">Each read is one cell at (T+C, T+G). A tight cluster is a uniform
           composition; a spread or a second cluster is a sub-population (e.g. adapter-dimers or a
           contaminant). Peak at ({}, {}).</div>"#,
        hm_m.get("heatmap_peak_x").map(|v| v.to_string()).unwrap_or_default(),
        hm_m.get("heatmap_peak_y").map(|v| v.to_string()).unwrap_or_default(),
    );
    let heatmap_section = if let Some(pt) = p_total {
        tabs("hmba", &[
            ("before", (hm_callout.clone() + &heatmap_svg(total))),
            ("after", (hm_callout + &heatmap_svg(pt))),
        ])
    } else {
        hm_callout + &heatmap_svg(total)
    };

    // ---- before/after delta: a compact number table only (no repeated charts) --
    let compare_table = match p_total {
        Some(pt) => {
            let rows_b = total.reads.max(1) as f64;
            let rows_a = pt.reads.max(1) as f64;
            let row = |label: &str, b: f64, a: f64, unit: &str| {
                let d = a - b;
                let sign = if d > 0.0 { "+" } else { "" };
                format!(
                    "<tr><th>{label}</th><td>{:.3}{unit}</td><td>{:.3}{unit}</td><td>{sign}{:.3}{unit}</td></tr>",
                    b, a, d
                )
            };
            let mut rows = String::new();
            rows.push_str(&row("reads", rows_b, rows_a, ""));
            rows.push_str(&row("%GC", gc_fraction(total) * 100.0, gc_fraction(pt) * 100.0, "%"));
            rows.push_str(&row("mean read Q", mean_read_mean_q(total, offset), mean_read_mean_q(pt, offset), ""));
            let (_hg, hd) = crate::rich::head_dropout(total);
            let (_phg, phd) = crate::rich::head_dropout(pt);
            rows.push_str(&row("head dropout (G)", hd * 100.0, phd * 100.0, "%"));
            format!(
                r#"<table class="rep"><tr><th>metric</th><th>before</th><th>after</th><th>delta</th></tr>{rows}</table>"#
            )
        }
        _ => String::new(),
    };

    // ---- observations: numbers, not verdicts --------------------------------
    let mut obs: Vec<String> = Vec::new();
    if n_pos > 0 {
        let d0 = depth_at(0);
        if d0 > 0 {
            let g0 = 100.0 * total.base_at(0, 2) as f64 / d0 as f64;
            let g_all = 100.0 * (0..n_pos).map(|p| total.base_at(p, 2)).sum::<u64>() as f64
                / (0..n_pos).map(depth_at).sum::<u64>().max(1) as f64;
            obs.push(format!(
                "first base is G in {g0:.2}% of reads against {g_all:.2}% for the library as a whole —                  a first-cycle G dropout shows up as the first being far higher"));
        }
        if n_pos >= 2 {
            let d_last = depth_at(n_pos - 1);
            if d_last > 0 {
                let g_last = 100.0 * total.base_at(n_pos - 1, 2) as f64 / d_last as f64;
                let d_mid = depth_at(n_pos / 2);
                let g_mid = if d_mid > 0 { 100.0 * total.base_at(n_pos / 2, 2) as f64 / d_mid as f64 } else { 0.0 };
                obs.push(format!(
                    "last base is G in {g_last:.2}% against {g_mid:.2}% mid-read — a polyG tail is the usual cause                      when the end is far higher"));
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
            "{} reads exceed --max-read-len and are counted in the summary but excluded from the              per-position charts", total.reads_over_max_len));
    }
    // Only a run that read every record may claim exactness: a budget-stopped
    // scan covers a prefix of each file, and saying otherwise next to a JSON
    // report that says partial.is_partial = true is a contradiction.
    if partial.is_none() {
        obs.push("every count in this report is exact: no sampling, no read was skipped".to_string());
    }

    // ---- per-file table (raw files, then processed files when present) -------
    let mut file_rows = String::new();
    for f in files {
        let _ = write!(
            file_rows,
            "<tr><td><code>{}</code></td><td>{}</td><td>{}</td><td>{:.2}%</td><td>{:.2}</td><td>{}</td></tr>",
            esc(&f.path.file_name().unwrap_or_default().to_string_lossy()),
            f.acc.reads, f.acc.bases, gc_fraction(&f.acc) * 100.0, mean_read_mean_q(&f.acc, offset),
            esc(&f.scan.decompress)
        );
    }
    for f in p_files.unwrap_or_default() {
        let _ = write!(
            file_rows,
            "<tr><td><code>{} (processed)</code></td><td>{}</td><td>{}</td><td>{:.2}%</td><td>{:.2}</td><td>{}</td></tr>",
            esc(&f.path.file_name().unwrap_or_default().to_string_lossy()),
            f.acc.reads, f.acc.bases, gc_fraction(&f.acc) * 100.0, mean_read_mean_q(&f.acc, offset),
            esc(&f.scan.decompress)
        );
    }

    // ---- assemble the page ---------------------------------------------------
    let mut h = String::with_capacity(64 * 1024);
    h.push_str("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">");
    let _ = write!(h, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>fastcheck report — {} file(s)</title><style>{STYLE}</style></head><body>", files.len());
    h.push_str(r##"<a class="skip-link" href="#main">Skip to report content</a>"##);
    // A partial run says so before any number: the same fact the JSON publishes
    // as partial.is_partial. Silence here made a --reads 2 page look like a
    // complete run of a two-read file.
    if let Some(reason) = partial {
        let _ = write!(
            h,
            "<div class=\"banner\" role=\"status\"><strong>Partial run — the scan stopped at the {reason}.</strong> \
             Every count below covers only the reads that were read, not the whole file. \
             The JSON report carries the same state (<code>partial.is_partial = true</code>, reason {reason:?}).</div>"
        );
    }
    let _ = write!(h, "<h1>fastcheck report</h1><div class=\"sub\">{} file(s) · {} reads · {} bases · \
                       {:.2} s · Phred+{offset} · {}</div>",
                   files.len(), total.reads, total.bases, elapsed_s,
                   match partial {
                       Some(reason) => format!("partial run, stopped at the {reason}"),
                       None => "exact counts, no sampling".to_string(),
                   });

    // nav: sticky section links, scrollspy-highlighted
    h.push_str(r##"<nav class="nav" aria-label="Report sections"><a href="#stats">General stats</a><a href="#summary">Summary</a><a href="#seq">Per-base content</a><a href="#n">N content</a><a href="#q">Quality</a><a href="#dist">Distributions</a><a href="#polyx">polyX</a><a href="#repeats">Repeats</a><a href="#heatmap">Heatmap</a>"##);
    if processed.is_some() {
        h.push_str(r##"<a href="#compare">Before / After</a>"##);
    }
    h.push_str(r##"<a href="#obs">Observations</a></nav>"##);

    // MultiQC-style glance: General statistics table at the very top.
    h.push_str(r##"<section class="module" aria-labelledby="stats"><h2 id="stats">General statistics</h2>
                  <table class="files"><caption>Per-file metrics</caption><thead><tr><th scope="col">file</th><th scope="col">reads</th><th scope="col">bases</th><th scope="col">%GC</th><th scope="col">mean Q</th><th scope="col">decompression</th></tr></thead><tbody>"##);
    h.push_str(&file_rows);
    h.push_str("</tbody></table></section>");

    h.push_str("<main id=\"main\">");

    h.push_str(r##"<details class="module" open><summary><h2 id="summary">Summary</h2></summary><div class="body"><div class="grid2"><div>"##);
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
    h.push_str("</div></div></div></details>");

    if processed.is_some() {
        h.push_str(r##"<details class="module" open><summary><h2 id="compare">Before / After</h2></summary><div class="body">"##);
        h.push_str(&compare_table);
        h.push_str("<p class=\"note\">The per-graph before/after comparison is available as a tab in each module above.</p>");
        h.push_str("</div></details>");
    }

    h.push_str(r##"<details class="module" open><summary><h2 id="seq">Per-base sequence content</h2></summary><div class="body">"##);
    let _ = write!(h, "{seq_chart}");
    h.push_str("<div class=\"note\">Base composition at each position. Use the R1 / R2 / both tabs to view a single mate or the merged library.</div></div></details>");

    h.push_str(r##"<details class="module" open><summary><h2 id="n">N content by position</h2></summary><div class="body">"##);
    let _ = write!(h, "{n_chart}");
    h.push_str("</div></details>");

    h.push_str(r##"<details class="module" open><summary><h2 id="q">Per-base quality</h2></summary><div class="body">"##);
    let _ = write!(h, "{q_chart}");
    h.push_str("<div class=\"note\">Mean quality by position with the 10th and 90th percentiles — a widening spread is quality falling off unevenly between reads, not just a lower mean.</div></div></details>");

    h.push_str(r##"<details class="module" open><summary><h2 id="dist">Distributions</h2></summary><div class="body">"##);
    let _ = write!(h, "{dist_section}");
    h.push_str("</div></details>");

    if insert.is_some() {
        h.push_str(r##"<details class="module" open><summary><h2 id="insert">Insert sizes</h2></summary><div class="body">"##);
        h.push_str(&insert_table);
        let _ = write!(h, "<div class=\"grid2\"><div>{insert_chart}</div><div></div></div>");
        h.push_str("</div></details>");
    }

    h.push_str(r##"<details class="module" open><summary><h2 id="polyx">polyX run lengths</h2></summary><div class="body">"##);
    h.push_str(&polyx_section);
    h.push_str("</div></details>");

    h.push_str(r##"<details class="module" open><summary><h2 id="repeats">Repeated sequences</h2></summary><div class="body">"##);
    h.push_str(&repeats_section);
    h.push_str("</div></details>");

    h.push_str(r##"<details class="module" open><summary><h2 id="heatmap">Heatmap</h2></summary><div class="body">"##);
    h.push_str(&heatmap_section);
    h.push_str("</div></details>");

    h.push_str(r##"<details class="module" open><summary><h2 id="obs">Observations</h2></summary><div class="body"><ul class="obs">"##);
    for o in &obs {
        let _ = write!(h, "<li>{o}</li>");
    }
    h.push_str("</ul></div></details>");

    h.push_str("</main>");
    let _ = write!(h, "<footer>fastcheck</footer>");
    // The only script in the report: tab switching (ARIA + keyboard) and the
    // sticky-nav scrollspy. Everything else is static inline SVG, so the
    // report still opens offline.
    h.push_str(
        r##"<script>
// Tabs: delegated, index-based. The i-th tab of a container shows the i-th
// panel of the SAME container. Only direct children count (:scope > ...), so
// nested tablists (e.g. R1/R2/both inside a before/after tab) never bleed
// into an outer container's indices, and duplicate panel ids from a
// before/after reuse can never resolve to the wrong panel.
function activateTab(list, i){
  var tabs = Array.prototype.slice.call(list.querySelectorAll(':scope > [role="tab"]'));
  var panels = Array.prototype.slice.call(list.querySelectorAll(':scope > [role="tabpanel"]'));
  tabs.forEach(function(b, j){
    var on = (j === i);
    b.classList.toggle('active', on);
    b.setAttribute('aria-selected', on ? 'true' : 'false');
    b.tabIndex = on ? 0 : -1;
  });
  panels.forEach(function(p, j){ p.classList.toggle('active', j === i); });
}
document.addEventListener('click', function(e){
  var el = e.target.closest ? e.target.closest('[role="tab"]') : null;
  if (!el) return;
  var list = el.closest('[role="tablist"]');
  if (!list) return;
  activateTab(list, Array.prototype.indexOf.call(list.querySelectorAll('[role="tab"]'), el));
});
document.addEventListener('keydown', function(e){
  if (!/^(ArrowLeft|ArrowRight|Home|End)$/.test(e.key)) return;
  var el = e.target.closest ? e.target.closest('[role="tab"]') : null;
  if (!el || !el.parentElement) return;
  var list = el.parentElement;
  if (list.getAttribute('role') !== 'tablist') return;
  e.preventDefault();
  var tabs = Array.prototype.slice.call(list.querySelectorAll(':scope > [role="tab"]'));
  var i = tabs.indexOf(el);
  if (e.key === 'ArrowRight') i = (i + 1) % tabs.length;
  else if (e.key === 'ArrowLeft') i = (i - 1 + tabs.length) % tabs.length;
  else if (e.key === 'Home') i = 0;
  else if (e.key === 'End') i = tabs.length - 1;
  tabs[i].focus();
  activateTab(list, i);
});
// Scrollspy: mark the nav link of the section currently in view.
(function(){
  var links = {};
  Array.prototype.forEach.call(document.querySelectorAll('nav.nav a'), function(a){
    var id = (a.getAttribute('href') || '').slice(1);
    if (id) links[id] = a;
  });
  function set(spyId){
    Array.prototype.forEach.call(document.querySelectorAll('nav.nav a'), function(a){
      var on = (a.getAttribute('href') || '').slice(1) === spyId;
      if (on) a.setAttribute('aria-current', 'true');
      else a.removeAttribute('aria-current');
    });
  }
  var spy;
  if ('IntersectionObserver' in window) {
    spy = new IntersectionObserver(function(entries){
      entries.forEach(function(en){
        if (en.isIntersecting) set(en.target.id);
      });
    }, { rootMargin: '-40% 0px -55% 0px', threshold: 0 });
    Array.prototype.forEach.call(document.querySelectorAll('section.module h2[id], details.module h2[id]'), function(h2){
      spy.observe(h2);
    });
  }
})();
</script>"##,
    );
    h.push_str("</body></html>\n");


    // `-` means the stream (D-4), like every other report output: this used to write a
    // 46 KB file named `-` and leave stdout empty.
    crate::report::write_out(path, &h)
}
