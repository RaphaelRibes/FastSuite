//! Report model, adapters and writers.
//!
//! The suite writes one JSON per tool per sample. fasthelper normalises all of
//! them into a single metric table so a run can be summarised in one file.
//!
//! Two shapes are understood, deliberately:
//!   * `fastsuite.report.v1` — the canonical wrapper (schema/tool/metrics)
//!   * legacy raw JSONs — fastmaps' heatmap.json grid and pre-wrapper
//!     fastprep stats objects, so FastHelper works with binaries (and
//!     report files) that predate the wrapper.

use anyhow::{Context, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::cli::ReportArgs;

pub const SCHEMA_VERSION: &str = "fastsuite.report.v1";

pub const SCHEMA_DOC: &str = r#"fastsuite.report.v1 — per-tool report contract
=======================================================================
Every FastSuite tool writes exactly one JSON per run with --json <PATH>:

{
  "schema":       "fastsuite.report.v1",
  "tool":         "fastdedup",              // binary name
  "tool_version": "1.3.0",                  // CARGO_PKG_VERSION
  "command":      "fastdedup -1 a.fq ...",  // exact argv, space-joined
  "mode":         "paired",                 // optional, tool-defined
  "dry_run":      false,
  "inputs":       {"r1": "...", "r2": "..."},
  "outputs":      {"r1": "...", "r2": "..."},
  "elapsed_s":    12.5,                     // wall clock, seconds
  "metrics": {                              // flat numeric namespace
    "reads_in": 1000000,
    "reads_out": 998000,
    "reads_dropped": 2000
  }
}

Rules
-----
1. `metrics` values are numbers only (booleans as 0/1). Nested structures
   belong in tool-specific blocks outside `metrics`.
2. Metric names are snake_case and tool-prefixed only when ambiguous
   (`reads_in`, not `fastdedup_reads_in`).
3. Recommended shared keys, honoured by every tool that can produce them:
   reads_in, reads_out, reads_dropped, bases_in, bases_out, bases_trimmed,
   elapsed_s is top-level, none of the metric keys are mandatory.
4. Adding a metric is a MINOR change; renaming/removing one is a MAJOR change
   (bump the schema string, keep the old key one release as an alias).
5. FastHelper also accepts legacy raw JSON (pre-wrapper fastprep stats
   objects, fastmaps heatmap.json) and derives the keys above; every new
   tool emits v1.

Producers today
---------------
  fastdedup   --json  -> v1 (native)
  fastprep    --json  -> v1 (native)
  fastmaps    --json  -> heatmap.json in --output-dir (v1 wrapper planned)
  fastsim     -       -> no JSON yet (planned: --json with truth summary)
"#;

/// One tool's report for one sample, normalised.
#[derive(Debug, Clone)]
pub struct ToolReport {
    pub sample: String,
    pub tool: String,
    pub version: Option<String>,
    pub source: PathBuf,
    pub metrics: BTreeMap<String, f64>,
    pub notes: Vec<String>,
}

fn num(v: &serde_json::Value) -> Option<f64> {
    match v {
        serde_json::Value::Number(n) => n.as_f64(),
        serde_json::Value::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
        _ => None,
    }
}

fn get<'a>(o: &'a serde_json::Value, key: &str) -> Option<&'a serde_json::Value> {
    o.get(key)
}

fn get_num(o: &serde_json::Value, key: &str) -> f64 {
    o.get(key).and_then(num).unwrap_or(0.0)
}

/// Sample name from a path: file stem, unless that stem is a generic report
/// name — then the parent directory carries the sample identity.
pub fn sample_from_path(path: &Path) -> String {
    const GENERIC: &[&str] = &[
        "report", "reports", "stats", "json", "out", "output", "heatmap", "fastprep",
        "fastdedup", "fastmaps", "fastsim", "fasthelper", "dedup", "prep", "maps", "result",
        "results", "summary", "manifest", "positions", "position", "pos", "matrix",
    ];
    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
        .to_string();
    if GENERIC.contains(&stem.to_ascii_lowercase().as_str()) {
        // Walk up until a directory that is not itself a generic report name:
        // runs/demo/<sample>/maps/heatmap.json must yield "<sample>", not "maps".
        let mut cur = path.parent();
        while let Some(dir) = cur {
            match dir.file_name().and_then(|s| s.to_str()) {
                Some(name)
                    if !name.is_empty()
                        && name != "."
                        && !GENERIC.contains(&name.to_ascii_lowercase().as_str()) =>
                {
                    return name.to_string();
                }
                Some(_) => cur = dir.parent(),
                None => break,
            }
        }
    }
    stem
}

/// Normalise one JSON file. Returns Ok(None) when the file is not a report we
/// recognise (caller decides whether that is worth a warning).
pub fn adapt_file(path: &Path, sample_hint: Option<&str>) -> Result<Option<ToolReport>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;
    let v: serde_json::Value = match serde_json::from_str(&text) {
        Ok(v) => v,
        Err(e) => return Ok(Some(ToolReport {
            sample: sample_hint.map(str::to_string).unwrap_or_else(|| sample_from_path(path)),
            tool: "unparsable".into(),
            version: None,
            source: path.to_path_buf(),
            metrics: BTreeMap::new(),
            notes: vec![format!("not valid JSON: {e}")],
        })),
    };
    if !v.is_object() {
        return Ok(None);
    }
    // fasthelper's own aggregate (report.json) carries the same schema string
    // as a per-tool report but is a collection, not a measurement. Skip it by
    // its generator/reports shape, so re-concatenating a run directory never
    // turns fasthelper's report.json into a phantom sample.
    if v.get("generator").and_then(|g| g.as_str()) == Some("fasthelper")
        || v.get("reports").is_some_and(|r| r.is_array())
    {
        return Ok(None);
    }
    let sample = sample_hint
        .map(str::to_string)
        .unwrap_or_else(|| sample_from_path(path));

    // --- canonical v1 wrapper ---
    if v.get("schema").and_then(|s| s.as_str()) == Some(SCHEMA_VERSION) {
        let tool = v
            .get("tool")
            .and_then(|t| t.as_str())
            .unwrap_or("unknown")
            .to_string();
        let mut metrics = BTreeMap::new();
        if let Some(m) = v.get("metrics").and_then(|m| m.as_object()) {
            for (k, val) in m {
                if let Some(f) = num(val) {
                    metrics.insert(k.clone(), f);
                }
            }
        }
        // Derived convenience keys, so every sample has a comparable pair.
        let rin = metrics.get("reads_in").copied();
        let rout = metrics.get("reads_out").copied();
        if let (Some(rin), Some(rout)) = (rin, rout) {
            if rin > 0.0 {
                metrics.insert("read_retention".into(), rout / rin);
            }
        }
        if let Some(e) = v.get("elapsed_s").and_then(num) {
            metrics.insert("elapsed_s".into(), e);
        }
        return Ok(Some(ToolReport {
            sample,
            tool,
            version: v
                .get("tool_version")
                .and_then(|s| s.as_str())
                .map(str::to_string),
            source: path.to_path_buf(),
            metrics,
            notes: Vec::new(),
        }));
    }

    // --- fastmaps heatmap.json (grid of (T+C, T+G) counts) ---
    if v.get("counts").is_some() && v.get("size").is_some() {
        let size = get_num(&v, "size") as i64;
        let windows: [(&str, &serde_json::Value); 3] = [
            ("all", &v),
            ("r1", get(&v, "counts_r1").unwrap_or(&serde_json::Value::Null)),
            ("r2", get(&v, "counts_r2").unwrap_or(&serde_json::Value::Null)),
        ];
        let mut metrics = BTreeMap::new();
        metrics.insert("total_reads".into(), get_num(&v, "total_reads"));
        metrics.insert("grid_size".into(), size as f64);
        let c = (size / 50).max(2); // corner window: 2% of the read length
        for (label, obj) in windows {
            let counts = match obj.get("counts") {
                Some(serde_json::Value::Array(a)) => a,
                _ => continue,
            };
            if counts.len() < (size * size) as usize {
                continue;
            }
            let (mut poly_a, mut poly_t, mut poly_g, mut poly_c, mut off_diag) =
                (0.0f64, 0.0, 0.0, 0.0, 0.0);
            let mut total = 0.0f64;
            for y in 0..size {
                for x in 0..size {
                    let n = counts[(y * size + x) as usize].as_f64().unwrap_or(0.0);
                    if n == 0.0 {
                        continue;
                    }
                    total += n;
                    let low_x = x <= c;
                    let high_x = x >= size - 1 - c;
                    let low_y = y <= c;
                    let high_y = y >= size - 1 - c;
                    if low_x && low_y {
                        poly_a += n;
                    }
                    if high_x && high_y {
                        poly_t += n;
                    }
                    if low_x && high_y {
                        poly_g += n;
                    }
                    if high_x && low_y {
                        poly_c += n;
                    }
                    if (x - y).abs() > size / 10 {
                        off_diag += n;
                    }
                }
            }
            if total == 0.0 {
                continue;
            }
            let p = |n: f64| n / total;
            let pre = if label == "all" { String::new() } else { format!("{label}_") };
            metrics.insert(format!("{pre}polyA_frac"), p(poly_a));
            metrics.insert(format!("{pre}polyT_frac"), p(poly_t));
            metrics.insert(format!("{pre}polyG_frac"), p(poly_g));
            metrics.insert(format!("{pre}polyC_frac"), p(poly_c));
            metrics.insert(format!("{pre}off_diagonal_frac"), p(off_diag));
            if label == "all" {
                metrics.insert("total_reads".into(), total);
            }
        }
        return Ok(Some(ToolReport {
            sample,
            tool: "fastmaps".into(),
            version: v.get("version").and_then(num).map(|n| format!("grid v{n}")),
            source: path.to_path_buf(),
            metrics,
            notes: Vec::new(),
        }));
    }

    // --- fastprep legacy stats: {r1: {...}, r2: {...}, pairs_out, ...} ---
    // Only reached by report files written before fastprep emitted the v1
    // wrapper (PLAN.md 1.1). A v1 report is claimed by the branch above, so
    // this note never lands on a current fastprep --json file.
    if let Some(r1) = v.get("r1").filter(|x| x.is_object()) {
        let r2 = v.get("r2").filter(|x| x.is_object());
        let mut metrics = BTreeMap::new();
        let mut reads_in = 0.0;
        let mut reads_kept = 0.0;
        let mut reads_dropped = 0.0;
        let mut bases_in = 0.0;
        let mut bases_trimmed = 0.0;
        let mut adapter_trimmed = 0.0;
        let mut poly_dropped = 0.0;
        let mut hard = 0.0;
        let mut mean_q = 0.0;
        let mut adapter_dropped = 0.0;
        let mut duplicate_dropped = 0.0;
        for mate in [Some(r1), r2].into_iter().flatten() {
            reads_in += get_num(mate, "reads_in");
            reads_kept += get_num(mate, "reads_kept");
            reads_dropped += get_num(mate, "reads_dropped");
            bases_in += get_num(mate, "bases_in");
            bases_trimmed += get_num(mate, "bases_trimmed_5") + get_num(mate, "bases_trimmed_3");
            if let Some(a) = get(mate, "adapter") {
                adapter_trimmed += get_num(a, "reads_trimmed");
            }
            for cause in [
                "reads_dropped_poly5",
                "reads_dropped_poly3",
                "reads_dropped_poly_int",
            ] {
                poly_dropped += get_num(mate, cause);
            }
            hard += get_num(mate, "reads_dropped_hard");
            mean_q += get_num(mate, "reads_dropped_mean_q");
            adapter_dropped += get_num(mate, "reads_dropped_adapter");
            duplicate_dropped += get_num(mate, "reads_dropped_duplicate");
        }
        // Per-cause drop counts are summed across mates: this is a read-level
        // budget, and reporting only the last mate silently halves it.
        metrics.insert("reads_dropped_hard".into(), hard);
        metrics.insert("reads_dropped_mean_q".into(), mean_q);
        metrics.insert("reads_dropped_adapter".into(), adapter_dropped);
        metrics.insert("reads_dropped_duplicate".into(), duplicate_dropped);
        metrics.insert("reads_in".into(), reads_in);
        metrics.insert("reads_out".into(), reads_kept);
        metrics.insert("reads_dropped".into(), reads_dropped);
        metrics.insert("bases_in".into(), bases_in);
        metrics.insert("bases_trimmed".into(), bases_trimmed);
        metrics.insert("reads_adapter_trimmed".into(), adapter_trimmed);
        metrics.insert("reads_dropped_polyx".into(), poly_dropped);
        if reads_in > 0.0 {
            metrics.insert("read_retention".into(), reads_kept / reads_in);
            metrics.insert("bases_trimmed_frac".into(), bases_trimmed / bases_in.max(1.0));
        }
        if v.get("pairs_out").is_some() {
            metrics.insert("pairs_out".into(), get_num(&v, "pairs_out"));
            metrics.insert("pairs_dropped".into(), get_num(&v, "pairs_dropped"));
            metrics.insert("singletons_r1".into(), get_num(&v, "singletons_r1"));
            metrics.insert("singletons_r2".into(), get_num(&v, "singletons_r2"));
        }
        return Ok(Some(ToolReport {
            sample,
            tool: "fastprep".into(),
            version: None,
            source: path.to_path_buf(),
            metrics,
            notes: vec!["legacy raw stats (no schema wrapper)".into()],
        }));
    }

    // --- anything else: take top-level numbers if it smells like a report ---
    let obj = v.as_object().unwrap();
    let numeric: BTreeMap<String, f64> = obj
        .iter()
        .filter_map(|(k, val)| num(val).map(|n| (k.clone(), n)))
        .collect();
    if numeric.is_empty() {
        return Ok(None);
    }
    Ok(Some(ToolReport {
        sample,
        tool: obj
            .get("tool")
            .and_then(|t| t.as_str())
            .unwrap_or("unknown")
            .to_string(),
        version: obj
            .get("tool_version")
            .and_then(|t| t.as_str())
            .map(str::to_string),
        source: path.to_path_buf(),
        metrics: numeric,
        notes: vec!["unrecognised schema: top-level numbers only".into()],
    }))
}

fn walk_json(dir: &Path, out: &mut Vec<PathBuf>, depth: usize) {
    if depth > 6 {
        return;
    }
    let Ok(rd) = std::fs::read_dir(dir) else { return };
    for entry in rd.flatten() {
        let p = entry.path();
        let name = entry.file_name().to_string_lossy().to_string();
        if p.is_dir() {
            if matches!(
                name.as_str(),
                "target" | ".git" | ".pixi" | "node_modules" | "logs" | "tmp"
            ) {
                continue;
            }
            walk_json(&p, out, depth + 1);
        } else if name.ends_with(".json")
            && !matches!(name.as_str(), "pixi.lock" | "Cargo.lock" | "manifest.json")
            && !p.to_string_lossy().contains(".rustc_info")
        {
            out.push(p);
        }
    }
}

/// Collect and normalise every report reachable from the given inputs.
pub fn collect(inputs: &[PathBuf], sample: Option<&str>) -> Result<Vec<ToolReport>> {
    let mut files = Vec::new();
    for input in inputs {
        if input.is_dir() {
            walk_json(input, &mut files, 0);
        } else if input.is_file() {
            files.push(input.clone());
        } else {
            anyhow::bail!("input not found: {}", input.display());
        }
    }
    files.sort();
    files.dedup();

    let mut reports = Vec::new();
    for (i, f) in files.iter().enumerate() {
        let hint = if files.len() == 1 { sample } else { None };
        match adapt_file(f, hint) {
            Ok(Some(r)) => reports.push(r),
            Ok(None) => {}
            Err(e) => eprintln!("fasthelper: skipping {}: {e}", f.display()),
        }
        let _ = i;
    }
    Ok(reports)
}

/// Wide metric matrix: one row per (sample, tool), union of metric columns.
pub struct Matrix {
    pub columns: Vec<String>,
    pub rows: Vec<(String, String, BTreeMap<String, f64>)>,
}

pub fn build_matrix(reports: &[ToolReport]) -> Matrix {
    let mut columns: BTreeSet<String> = BTreeSet::new();
    let mut rows: Vec<(String, String, BTreeMap<String, f64>)> = Vec::new();
    for r in reports {
        columns.extend(r.metrics.keys().cloned());
        rows.push((r.sample.clone(), r.tool.clone(), r.metrics.clone()));
    }
    Matrix {
        columns: columns.into_iter().collect(),
        rows,
    }
}

fn fmt_val(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else if v.abs() >= 1000.0 {
        format!("{v:.0}")
    } else if v.abs() >= 1.0 {
        format!("{v:.3}")
    } else {
        format!("{v:.5}")
    }
}

pub fn matrix_tsv(m: &Matrix) -> String {
    let mut s = String::from("sample\ttool");
    for c in &m.columns {
        s.push('\t');
        s.push_str(c);
    }
    s.push('\n');
    for (sample, tool, metrics) in &m.rows {
        s.push_str(sample);
        s.push('\t');
        s.push_str(tool);
        for c in &m.columns {
            s.push('\t');
            if let Some(v) = metrics.get(c) {
                s.push_str(&fmt_val(*v));
            }
        }
        s.push('\n');
    }
    s
}

pub fn matrix_md(reports: &[ToolReport]) -> String {
    let by_tool = group_by_tool(reports);
    let mut s = String::new();
    for (tool, reports) in &by_tool {
        let mut columns: BTreeSet<&str> = BTreeSet::new();
        for r in reports {
            columns.extend(r.metrics.keys().map(|k| k.as_str()));
        }
        s.push_str(&format!("## {tool}\n\n| metric |"));
        for r in reports {
            s.push_str(&format!(" {} |", r.sample));
        }
        s.push_str("\n|---|");
        for _ in reports {
            s.push_str("---|");
        }
        s.push('\n');
        for c in &columns {
            s.push_str(&format!("| `{c}` |"));
            for r in reports {
                s.push_str(&format!(
                    " {} |",
                    r.metrics.get(*c).map(|v| fmt_val(*v)).unwrap_or_default()
                ));
            }
            s.push('\n');
        }
        s.push('\n');
    }
    s
}

fn group_by_tool(reports: &[ToolReport]) -> BTreeMap<String, Vec<&ToolReport>> {
    let mut m: BTreeMap<String, Vec<&ToolReport>> = BTreeMap::new();
    for r in reports {
        m.entry(r.tool.clone()).or_default().push(r);
    }
    m
}

// ---------------------------------------------------------------- HTML ----

fn bar(value: f64, max: f64, colour: &str) -> String {
    let w = if max > 0.0 {
        ((value / max) * 180.0).clamp(0.0, 180.0)
    } else {
        0.0
    };
    format!(
        "<svg width='200' height='14' style='vertical-align:middle'>\
         <rect x='0' y='2' width='{w:.1}' height='10' fill='{colour}' rx='2'/></svg>"
    )
}

/// Self-contained HTML report: no scripts, no external assets.
pub fn render_html(title: &str, reports: &[ToolReport]) -> String {
    let by_tool = group_by_tool(reports);
    let mut out = String::new();
    out.push_str(&format!(
        "<!doctype html><html><head><meta charset='utf-8'><title>{title}</title>\
<style>\
body{{font-family:system-ui,sans-serif;margin:2rem;color:#1c1e21;background:#fff}}\
h1{{font-size:1.4rem;margin:0 0 .2rem}}h2{{font-size:1.05rem;margin:1.6rem 0 .4rem;border-bottom:1px solid #ddd;padding-bottom:.2rem}}\
table{{border-collapse:collapse;font-size:.82rem;margin:.4rem 0}}\
th,td{{border:1px solid #e3e3e3;padding:.25rem .5rem;text-align:right}}\
th:first-child,td:first-child{{text-align:left}}\
th{{background:#f6f7f9;font-weight:600}}\
tr:nth-child(even) td{{background:#fbfbfc}}\
.meta{{color:#616770;font-size:.85rem;margin-bottom:1rem}}\
.note{{color:#8a6d3b;font-size:.8rem}}\
code{{background:#f2f3f5;padding:0 .25rem;border-radius:3px}}\
</style></head><body>"
    ));
    out.push_str(&format!("<h1>{title}</h1>"));
    out.push_str(&format!(
        "<div class='meta'>samples: {} &middot; reports: {} &middot; tools: {}</div>",
        reports
            .iter()
            .map(|r| r.sample.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        reports.len(),
        by_tool.keys().cloned().collect::<Vec<_>>().join(", ")
    ));

    for (tool, reps) in &by_tool {
        let versions: BTreeSet<String> = reps
            .iter()
            .filter_map(|r| r.version.clone())
            .collect();
        out.push_str(&format!(
            "<h2>{tool}{}</h2>",
            if versions.is_empty() {
                String::new()
            } else {
                format!(" <span class='meta'>{}</span>", versions.into_iter().collect::<Vec<_>>().join(", "))
            }
        ));
        let mut columns: BTreeSet<&str> = BTreeSet::new();
        for r in reps {
            columns.extend(r.metrics.keys().map(|k| k.as_str()));
        }
        out.push_str("<table><tr><th>metric</th>");
        for r in reps {
            out.push_str(&format!("<th>{}</th>", r.sample));
        }
        out.push_str("</tr>");
        for c in &columns {
            let numeric: Vec<Option<f64>> = reps.iter().map(|r| r.metrics.get(*c).copied()).collect();
            let max = numeric.iter().flatten().cloned().fold(0.0f64, f64::max);
            let is_rate = c.ends_with("_frac") || c.ends_with("rate") || c.ends_with("retention");
            out.push_str(&format!("<tr><td><code>{c}</code></td>"));
            for v in &numeric {
                match v {
                    Some(v) => {
                        if is_rate && max > 0.0 {
                            out.push_str(&format!(
                                "<td>{} {}</td>",
                                fmt_val(*v),
                                bar(*v, max, "#4c78a8")
                            ));
                        } else {
                            out.push_str(&format!("<td>{}</td>", fmt_val(*v)));
                        }
                    }
                    None => out.push_str("<td></td>"),
                }
            }
            out.push_str("</tr>");
        }
        out.push_str("</table>");
        for r in reps {
            for n in &r.notes {
                out.push_str(&format!(
                    "<div class='note'>{} — {}: {}</div>",
                    r.sample, r.tool, n
                ));
            }
        }
    }
    out.push_str("</body></html>");
    out
}

pub fn render_json(title: &str, reports: &[ToolReport]) -> String {
    let samples: Vec<serde_json::Value> = reports
        .iter()
        .map(|r| {
            serde_json::json!({
                "sample": r.sample,
                "tool": r.tool,
                "tool_version": r.version,
                "source": r.source.to_string_lossy(),
                "metrics": r.metrics,
                "notes": r.notes,
            })
        })
        .collect();
    let doc = serde_json::json!({
        "schema": SCHEMA_VERSION,
        "generator": "fasthelper",
        "generator_version": env!("CARGO_PKG_VERSION"),
        "title": title,
        "reports": samples,
    });
    serde_json::to_string_pretty(&doc).unwrap_or_else(|_| "{}".into()) + "\n"
}

/// Write the aggregate report in the requested formats. Returns written paths.
pub fn write_reports(
    title: &str,
    reports: &[ToolReport],
    outdir: &Path,
    formats: &[String],
) -> Result<Vec<PathBuf>> {
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("creating {}", outdir.display()))?;
    let mut written = Vec::new();
    for f in formats {
        let (name, body) = match f.trim() {
            "json" => ("report.json", render_json(title, reports)),
            "tsv" => ("report.tsv", matrix_tsv(&build_matrix(reports))),
            "md" => ("report.md", matrix_md(reports)),
            "html" => ("report.html", render_html(title, reports)),
            other => anyhow::bail!("unknown format '{other}' (json,tsv,md,html)"),
        };
        let path = outdir.join(name);
        std::fs::write(&path, body).with_context(|| format!("writing {}", path.display()))?;
        written.push(path);
    }
    Ok(written)
}

/// One `metric OP value` threshold, e.g. `polyG_frac>0.05`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Op {
    Gt,
    Lt,
    Ge,
    Le,
    Eq,
    Ne,
}

impl Op {
    fn as_str(&self) -> &'static str {
        match self {
            Op::Gt => ">",
            Op::Lt => "<",
            Op::Ge => ">=",
            Op::Le => "<=",
            Op::Eq => "==",
            Op::Ne => "!=",
        }
    }
}

/// One `metric OP value` threshold, e.g. `polyG_frac>0.05`.
#[derive(Debug, Clone)]
pub struct Threshold {
    pub metric: String,
    pub op: Op,
    pub value: f64,
}

impl Threshold {
    /// Parse `metric>0.05` / `reads<=100` / ... Ops: > < >= <= == !=.
    pub fn parse(expr: &str) -> Result<Self> {
        // Two-char operators first, so `>=` is never read as `>`.
        const OPS: &[(&str, Op)] = &[
            (">=", Op::Ge),
            ("<=", Op::Le),
            ("==", Op::Eq),
            ("!=", Op::Ne),
            (">", Op::Gt),
            ("<", Op::Lt),
        ];
        for (op, kind) in OPS {
            if let Some(i) = expr.find(op) {
                let metric = expr[..i].trim().to_string();
                let value: f64 = expr[i + op.len()..]
                    .trim()
                    .parse()
                    .with_context(|| format!("parsing threshold '{expr}': value is not a number"))?;
                if metric.is_empty() {
                    anyhow::bail!("threshold '{expr}' has no metric name");
                }
                return Ok(Threshold {
                    metric,
                    op: kind.clone(),
                    value,
                });
            }
        }
        anyhow::bail!("threshold '{expr}' is not METRIC OP VALUE (ops: > < >= <= == !=)");
    }

    pub fn tripped(&self, actual: f64) -> bool {
        match self.op {
            Op::Gt => actual > self.value,
            Op::Lt => actual < self.value,
            Op::Ge => actual >= self.value,
            Op::Le => actual <= self.value,
            Op::Eq => (actual - self.value).abs() < 1e-12,
            Op::Ne => (actual - self.value).abs() >= 1e-12,
        }
    }
}

/// Evaluate thresholds against EVERY report, not a flattened map.
///
/// Each tripped metric is returned labelled with the sample it came from, so a
/// violation in one sample of a multi-sample run is never masked by another
/// sample's value ("last wins" hid exactly that). Returns the tripped
/// descriptions, e.g. `a: duplicate_rate=0.5 > 0.01`.
pub fn evaluate_thresholds(reports: &[ToolReport], thresholds: &[Threshold]) -> Vec<String> {
    let mut out = Vec::new();
    for r in reports {
        for t in thresholds {
            if let Some(actual) = r.metrics.get(&t.metric) {
                if t.tripped(*actual) {
                    out.push(format!(
                        "{}: {}={} {} {}",
                        r.sample,
                        t.metric,
                        actual,
                        t.op.as_str(),
                        t.value
                    ));
                }
            }
        }
    }
    out
}

pub fn cmd_concatenate(a: &ReportArgs) -> Result<()> {
    let reports = collect(&a.inputs, a.sample.as_deref())?;
    if reports.is_empty() {
        anyhow::bail!(
            "no recognisable reports found in {}",
            a.inputs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    let formats: Vec<String> = a.format.split(',').map(|s| s.trim().to_string()).collect();
    let title = "FastSuite report";
    let written = write_reports(title, &reports, &a.outdir, &formats)?;
    if !a.fail_on.is_empty() {
        let thresholds: Vec<Threshold> = a
            .fail_on
            .iter()
            .map(|s| Threshold::parse(s))
            .collect::<Result<Vec<_>>>()?;
        let tripped = evaluate_thresholds(&reports, &thresholds);
        if !tripped.is_empty() {
            for t in &tripped {
                eprintln!("fasthelper: threshold tripped: {t}");
            }
            eprintln!(
                "fasthelper: {} threshold(s) tripped — report written, exiting 2",
                tripped.len()
            );
            std::process::exit(2);
        }
    }
    if a.stdout {
        print!("{}", matrix_tsv(&build_matrix(&reports)));
    }
    println!(
        "fasthelper: concatenated {} report(s) from {} sample(s) -> {}",
        reports.len(),
        reports
            .iter()
            .map(|r| r.sample.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        written
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
    Ok(())
}
