//! `fasthelper compare`: a before/after side-by-side report from two report
//! sets. Reuses the same discovery and normalisation as `concatenate report`
//! (`collect`/`adapt_file`), so it works on any `fastsuite.report.v1` JSON —
//! typically a fastcheck run on the raw reads and one on the processed reads.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::report::{collect, ToolReport};

/// Merge every metric across a report set into one map (last value wins on a
/// key collision, which is fine for a single-tool comparison).
fn merge_metrics(reports: &[ToolReport]) -> BTreeMap<String, f64> {
    let mut m = BTreeMap::new();
    for r in reports {
        for (k, v) in &r.metrics {
            m.insert(k.clone(), *v);
        }
    }
    m
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// Write `compare.tsv` and `compare.html` to `outdir`.
pub fn cmd_compare(before: &[PathBuf], after: &[PathBuf], outdir: &Path) -> Result<()> {
    let b = collect(before, None)?;
    let a = collect(after, None)?;
    if b.is_empty() || a.is_empty() {
        anyhow::bail!(
            "compare needs at least one report on each side (before: {}, after: {})",
            b.len(),
            a.len()
        );
    }
    let bm = merge_metrics(&b);
    let am = merge_metrics(&a);

    let mut keys: BTreeSet<String> = bm.keys().cloned().collect();
    keys.extend(am.keys().cloned());
    let mut rows: Vec<(String, f64, f64, f64)> = Vec::new();
    for k in &keys {
        if let (Some(bv), Some(av)) = (bm.get(k).copied(), am.get(k).copied()) {
            rows.push((k.clone(), bv, av, av - bv));
        }
    }
    rows.sort_by(|a, b| a.0.cmp(&b.0));

    // ---- TSV -----------------------------------------------------------------
    let mut tsv = String::from("metric\tbefore\tafter\tdelta\n");
    for (k, bv, av, d) in &rows {
        tsv.push_str(&format!("{k}\t{bv}\t{av}\t{d}\n"));
    }
    std::fs::write(outdir.join("compare.tsv"), tsv)
        .with_context(|| format!("writing {}", outdir.join("compare.tsv").display()))?;

    // ---- HTML ----------------------------------------------------------------
    let mut tr = String::new();
    for (k, bv, av, d) in &rows {
        let sign = if *d > 0.0 { "+" } else { "" };
        tr.push_str(&format!(
            "<tr><th>{}</th><td>{:.4}</td><td>{:.4}</td><td>{sign}{:.4}</td></tr>",
            esc(k), bv, av, d
        ));
    }
    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>fasthelper compare</title>
<style>
body {{ margin:0; padding:28px 34px 60px; color:#1d2530; background:#fff;
       font:14px/1.5 -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }}
h1 {{ font-size:23px; margin:0 0 2px; color:#2f6f9f; }}
.sub {{ color:#6b7580; font-size:13px; margin-bottom:20px; }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }}
th, td {{ text-align:right; padding:4px 12px 4px 0; border-bottom:1px solid #dde4ea; }}
th:first-child, td:first-child {{ text-align:left; }}
th {{ color:#6b7580; font-weight:600; }}
</style></head><body>
<h1>fasthelper compare</h1>
<div class="sub">{} report(s) before · {} report(s) after · {} metrics in common</div>
<table><tr><th>metric</th><th>before</th><th>after</th><th>delta</th></tr>{tr}</table>
</body></html>
"#,
        b.len(),
        a.len(),
        rows.len()
    );
    std::fs::write(outdir.join("compare.html"), html)
        .with_context(|| format!("writing {}", outdir.join("compare.html").display()))?;

    eprintln!(
        "fasthelper: compare wrote {} metrics to {} and {}",
        rows.len(),
        outdir.join("compare.tsv").display(),
        outdir.join("compare.html").display()
    );
    Ok(())
}
