//! The `fastsuite.report.v1` wrapper and the corner metrics it carries.

use std::collections::BTreeMap;

/// Fraction of reads in each heatmap corner, per mate window.
///
/// Reproduces the convention fasthelper historically derived from
/// `heatmap.json` (tools/fasthelper/src/report.rs): a (T+C, T+G) grid of
/// side `size`, corner window `c = max(size/50, 2)`:
///   low_x/low_y  → polyA corner, high_x/high_y → polyT,
///   low_x/high_y → polyG,      high_x/low_y  → polyC,
/// and off-diagonal mass where `|x - y| > size/10`.
///
/// `grids` is `[("all", &counts), ("r1", &counts_r1), ("r2", &counts_r2)]`
/// with flat row-major index `y*size + x`. Keys use the "" / "r1_" / "r2_"
/// prefixes, exactly as fasthelper's adapter did, so consumers see no
/// difference in metric names.
pub fn corner_metrics(size: u32, grids: &[(&str, &[u64])]) -> BTreeMap<String, f64> {
    let c = (size / 50).max(2) as usize;
    let diag_lim = size / 10;
    let mut out = BTreeMap::new();
    for &(label, counts) in grids {
        let (mut poly_a, mut poly_t, mut poly_g, mut poly_c, mut off_diag) =
            (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut total = 0.0f64;
        for y in 0..size as usize {
            for x in 0..size as usize {
                let n = counts[y * size as usize + x] as f64;
                if n == 0.0 {
                    continue;
                }
                total += n;
                let (low_x, high_x) = (x <= c, x >= size as usize - 1 - c);
                let (low_y, high_y) = (y <= c, y >= size as usize - 1 - c);
                if low_x && low_y { poly_a += n; }
                if high_x && high_y { poly_t += n; }
                if low_x && high_y { poly_g += n; }
                if high_x && low_y { poly_c += n; }
                if (x as i64 - y as i64).abs() as u32 > diag_lim { off_diag += n; }
            }
        }
        if total == 0.0 {
            continue;
        }
        let p = |n: f64| n / total;
        let pre = if label.is_empty() { String::new() } else { format!("{label}_") };
        out.insert(format!("{pre}polyA_frac"), p(poly_a));
        out.insert(format!("{pre}polyT_frac"), p(poly_t));
        out.insert(format!("{pre}polyG_frac"), p(poly_g));
        out.insert(format!("{pre}polyC_frac"), p(poly_c));
        out.insert(format!("{pre}off_diagonal_frac"), p(off_diag));
        if label.is_empty() {
            out.insert("total_reads".into(), total);
        }
    }
    out
}

/// Write the `fastsuite.report.v1` wrapper (docs/REPORT-SCHEMA.md) to `path`.
///
/// fastmaps' grid JSON (`heatmap.json`) keeps its own shape; this wrapper is
/// the machine-readable layer, carrying the corner metrics the tool now owns
/// (they used to be derived by fasthelper from the grid). The grid file's
/// location is recorded under `grid` so a consumer can find both.
pub fn write_report_v1(
    path: &std::path::Path,
    command: &str,
    elapsed_s: f64,
    files: &[String],
    metrics: &BTreeMap<String, f64>,
) -> std::io::Result<()> {
    use crate::encode::json_escape;

    let metrics_json = metrics
        .iter()
        .map(|(k, v)| format!("    {}: {}", json_escape(k), fmt_f64(*v)))
        .collect::<Vec<_>>()
        .join(",\n");
    let files_json = files
        .iter()
        .map(|f| json_escape(f))
        .collect::<Vec<_>>()
        .join(", ");

    let doc = format!(
        concat!(
            "{{\n",
            "  \"schema\": \"fastsuite.report.v1\",\n",
            "  \"tool\": \"fastmaps\",\n",
            "  \"tool_version\": \"{}\",\n",
            "  \"command\": {},\n",
            "  \"mode\": \"composition\",\n",
            "  \"inputs\": {{\"files\": [{files}]}},\n",
            "  \"elapsed_s\": {},\n",
            "  \"metrics\": {{\n{metrics}\n  }},\n",
            "  \"grid\": {{\"format\": \"heatmap.json\", \"see\": \"docs/PLAN-QC.md\"}}\n",
            "}}\n"
        ),
        env!("CARGO_PKG_VERSION"),
        json_escape(command),
        fmt_f64(elapsed_s),
        files = files_json,
        metrics = metrics_json,
    );

    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    use std::io::Write;
    f.write_all(doc.as_bytes())?;
    f.flush()
}

/// Compact float rendering: integers print without a `.0`, everything else
/// with 6 significant decimals. Keeps the hand-written JSON small.
fn fmt_f64(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        let s = format!("{x:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// size=11, c=2 fixture: 4 reads of length 10.
    /// R1 (file A): TTTTTTTTTT (x=10,y=10) and GGGGGGGGGG (x=0,y=10);
    /// R2 (file B): GGGGGGGGGG (x=0,y=10) and AAAAAAAAAA (x=0,y=0).
    /// Matches the verified fastmaps numbers in the convergence plan:
    /// all: polyA .25, polyT .25, polyG .5, polyC 0, off_diag .5, total 4;
    /// r1: polyT .5, polyG .5; r2: polyA .5, polyG .5.
    fn grids() -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        let size = 11;
        let mut all = vec![0u64; size * size];
        let mut r1 = vec![0u64; size * size];
        let mut r2 = vec![0u64; size * size];
        all[10 * size + 10] += 1; // T at (10,10), from file A's R1
        all[10 * size + 0] += 2;  // G at (0,10), once per file
        r1[10 * size + 10] += 1;
        r1[10 * size + 0] += 1;
        r2[10 * size + 0] += 1;
        all[0 * size + 0] += 1;   // A at (0,0), from file B's R2
        r2[0 * size + 0] += 1;
        (all, r1, r2)
    }

    #[test]
    fn corner_metrics_matches_fasthelper_fixture() {
        let (all, r1, r2) = grids();
        let m = corner_metrics(11, &[("", &all), ("r1", &r1), ("r2", &r2)]);
        let near = |k: &str, want: f64| {
            assert!((m[k] - want).abs() < 1e-9, "{k}: got {}, want {want}", m[k]);
        };
        near("polyA_frac", 0.25);
        near("polyT_frac", 0.25);
        near("polyG_frac", 0.5);
        near("polyC_frac", 0.0);
        near("off_diagonal_frac", 0.5);
        near("total_reads", 4.0);
        near("r1_polyT_frac", 0.5);
        near("r1_polyG_frac", 0.5);
        near("r1_polyA_frac", 0.0);
        near("r2_polyA_frac", 0.5);
        near("r2_polyT_frac", 0.0);
        near("r2_polyG_frac", 0.5);
    }

    #[test]
    fn corner_metrics_empty_grid_is_omitted() {
        // The two empty windows contribute no metrics at all; the populated
        // one still yields the full set. A banded grid keeps all corners and
        // the off-diagonal region empty — the opposite of the plan's original
        // "uniform grid" premise, which is corner-enriched by construction
        // (one pixel per corner vs. nine per interior pixel).
        let size = 11;
        let banded = || -> Vec<u64> {
            let mut v = vec![0u64; size * size];
            for y in 0..size {
                for x in 0..size {
                    if (x as i64 - y as i64).abs() <= 1 {
                        v[y * size + x] = 1;
                    }
                }
            }
            v
        };
        let all = vec![0u64; size * size];
        let r1 = vec![0u64; size * size];
        let r2 = banded();
        let m = corner_metrics(size as u32, &[("", &all), ("r1", &r1), ("r2", &r2)]);
        assert!(!m.contains_key("polyT_frac"), "empty window must be skipped");
        assert!(!m.contains_key("r1_polyT_frac"));
        assert!(!m.contains_key("total_reads"), "total_reads lives in the all window");
        assert!(m.contains_key("r2_polyT_frac"));
        // The two high-y/low-x corners (polyG, polyC) and the off-diagonal
        // band are empty for a banded grid. The polyA corner holds 7 banded
        // pixels (a 3x3 window with the two |x-y|=2 pixels absent) out of 31
        // total; polyT holds 7 likewise.
        assert!((m["r2_polyG_frac"] - 0.0).abs() < 1e-9);
        assert!((m["r2_polyC_frac"] - 0.0).abs() < 1e-9);
        assert!((m["r2_off_diagonal_frac"] - 0.0).abs() < 1e-9);
        assert!((m["r2_polyT_frac"] - 7.0 / 31.0).abs() < 1e-9);
        assert!((m["r2_polyA_frac"] - 7.0 / 31.0).abs() < 1e-9);
    }

    #[test]
    fn report_v1_written_and_parseable() {
        let dir = std::env::temp_dir().join(format!("maps_v1_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("r.json");
        let mut metrics = std::collections::BTreeMap::new();
        metrics.insert("total_reads".into(), 4.0);
        metrics.insert("polyT_frac".into(), 0.25);
        metrics.insert("r2_polyG_frac".into(), 0.5);
        super::write_report_v1(
            &path,
            "fastmaps a_R1.fq.gz a_R2.fastq.gz -j 2",
            0.42,
            &["a_R1.fq.gz".into(), "a_R2.fastq.gz".into()],
            &metrics,
        )
        .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("\"schema\": \"fastsuite.report.v1\""), "{text}");
        assert!(text.contains("\"tool\": \"fastmaps\""));
        assert!(text.contains("\"mode\": \"composition\""));
        assert!(text.contains("\"elapsed_s\": 0.42"));
        assert!(text.contains("\"polyT_frac\": 0.25"));
        assert!(text.contains("\"r2_polyG_frac\": 0.5"));
        assert!(text.contains("\"a_R1.fq.gz\""));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
