//! The grid JSON and the v1 report must not disagree about one run.
//!
//! `heatmap.json` described the plane (`counts`, `counts_r1`, `counts_r2`) but
//! carried no metrics block, so a consumer that recomputed the corner fractions
//! from the grid — fasthelper's documented `concatenate` adapter does exactly
//! that — saw `polyA_frac 1.000` where `--json` reported `0`: a read with no
//! ACGT base lands on `(0,0)` like an all-A read, and only `no_call_reads`
//! (which the report published and the grid did not) can take it back out.
//!
//! Both shapes now carry the same flat numeric namespace, and this pins it key
//! by key — plus the arithmetic that shows the published `no_call_reads` is
//! what reconciles the two files.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_parity_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn parse(path: &Path) -> Value {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{} must exist: {e}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("{} must parse: {e}\n{text}", path.display()))
}

/// One run of `seqs` (header, sequence) at `-j 1`; returns (grid, report).
fn run_case(tag: &str, seqs: &[(&str, String)]) -> (Value, Value) {
    let dir = tmpdir(tag);
    let fq = dir.join("m_R1.fastq");
    let mut body = String::new();
    for (h, s) in seqs {
        body.push_str(&format!("@{h}\n{s}\n+\n{}\n", "I".repeat(s.len())));
    }
    std::fs::write(&fq, body).unwrap();

    let out = dir.join("out");
    let report = dir.join("r.json");
    let o = Command::new(bin())
        .args(["-o", out.to_str().unwrap(), "-j", "1", "--force"])
        .arg("--json")
        .arg(&report)
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "{tag}: run failed\nstderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    (parse(&out.join("heatmap.json")), parse(&report))
}

/// The worker count is execution metadata, not a measurement of the data: the
/// grid is part of the `-j`-invariance contract (two runs of the same command at
/// different `-j` are byte-identical), so it is the one report metric the grid
/// does not carry.
const REPORT_ONLY: &[&str] = &["threads"];

fn assert_same_metrics(g: &Value, r: &Value) {
    let gm = g["metrics"]
        .as_object()
        .unwrap_or_else(|| panic!("heatmap.json must publish the metric namespace: {g}"));
    let rm = r["metrics"].as_object().expect("report metrics");

    // Every shared key agrees, value for value.
    for (k, want) in rm {
        if REPORT_ONLY.contains(&k.as_str()) {
            assert!(
                !gm.contains_key(k),
                "the grid must not carry the run's {k} (-j invariance)"
            );
            continue;
        }
        let got = gm
            .get(k)
            .unwrap_or_else(|| panic!("heatmap.json omits the metric {k}: {g}"));
        assert_eq!(got, want, "metric {k} differs between heatmap.json and the report");
    }
    // ...and the grid carries nothing the report does not, so a key cannot be
    // published in one shape only.
    for k in gm.keys() {
        assert!(
            rm.contains_key(k),
            "heatmap.json publishes {k}, which the report does not"
        );
    }
}

/// The reported case: one all-N read. The report says `polyA_frac 0`; a
/// counts-only reader of the grid says 1.000; `no_call_reads` is the difference.
#[test]
fn no_call_reads_reconciles_the_grid_with_the_report() {
    let (g, r) = run_case("nocall", &[("n0", "N".repeat(150))]);
    assert_same_metrics(&g, &r);

    let no_call = g["metrics"]["no_call_reads"]
        .as_f64()
        .unwrap_or_else(|| panic!("heatmap.json must carry no_call_reads: {g}"));
    assert_eq!(no_call, 1.0);
    assert_eq!(
        r["metrics"]["no_call_reads"].as_f64().unwrap(),
        1.0,
        "the report's no_call_reads is the same number"
    );

    // Recompute the polyA corner from the grid, the way a counts-only adapter
    // does, and show the correction is exactly what the grid now publishes.
    let size = g["size"].as_u64().unwrap() as usize;
    let c = (size / 50).max(2);
    let counts = g["counts"].as_array().unwrap();
    let corner_a: f64 = counts
        .iter()
        .enumerate()
        .filter(|(i, _)| (i % size) <= c && (i / size) <= c)
        .map(|(_, n)| n.as_f64().unwrap())
        .sum();
    let total = g["total_reads"].as_f64().unwrap();
    let naive = corner_a / total;
    assert!(
        (naive - 1.0).abs() < 1e-9,
        "a counts-only reader sees polyA_frac 1.000, got {naive}"
    );
    let denom = total - no_call;
    let corrected = if denom > 0.0 { (corner_a - no_call) / denom } else { 0.0 };
    assert!(
        (corrected - r["metrics"]["polyA_frac"].as_f64().unwrap()).abs() < 1e-9,
        "grid polyA {corner_a} - no_call {no_call} over {total} - {no_call} must equal the \
         report's polyA_frac"
    );
    assert_eq!(
        r["metrics"]["polyA_frac"].as_f64().unwrap(),
        0.0,
        "an all-N read is not polyA-rich"
    );
}

/// A measurable read beside the no-call one: both shapes say 1/2, and a
/// counts-only reader says 2/3.
#[test]
fn mixed_run_agrees_key_by_key() {
    let mid = "ACGT".repeat(37) + "AC";
    let (g, r) = run_case(
        "mixed",
        &[
            ("a0", "A".repeat(150)),
            ("m0", mid),
            ("n0", "N".repeat(150)),
        ],
    );
    assert_same_metrics(&g, &r);
    assert_eq!(g["metrics"]["no_call_reads"].as_f64().unwrap(), 1.0);
    assert_eq!(g["metrics"]["polyA_frac"].as_f64().unwrap(), 0.5);
    assert_eq!(r["metrics"]["polyA_frac"].as_f64().unwrap(), 0.5);
    assert_eq!(g["metrics"]["total_reads"].as_f64().unwrap(), 3.0);
    assert_eq!(g["metrics"]["r1_no_call_reads"].as_f64().unwrap(), 1.0);
}

/// The grid stays byte-identical across `-j`; only the report's `threads`
/// changes. This is why `threads` is report-only.
#[test]
fn grid_is_j_invariant_and_only_threads_differs() {
    let dir = tmpdir("jinv");
    // Two files, so `-j 2` really resolves to two workers (the pool is clamped
    // by the file count).
    let mut files = Vec::new();
    for mate in ["R1", "R2"] {
        let fq = dir.join(format!("j_{mate}.fastq"));
        let s = "A".repeat(60);
        let mut body = String::new();
        for i in 0..8 {
            body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(60)));
        }
        std::fs::write(&fq, body).unwrap();
        files.push(fq);
    }

    let mut grids = Vec::new();
    let mut reports = Vec::new();
    for j in ["1", "2"] {
        let out = dir.join(format!("out{j}"));
        let report = dir.join(format!("r{j}.json"));
        let o = Command::new(bin())
            .args(["-o", out.to_str().unwrap(), "-j", j, "--force"])
            .arg("--json")
            .arg(&report)
            .args(&files)
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "-j {j} failed: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        grids.push(std::fs::read(out.join("heatmap.json")).unwrap());
        reports.push(parse(&report));
    }
    assert_eq!(grids[0], grids[1], "heatmap.json must not depend on -j");
    assert_ne!(
        reports[0]["metrics"]["threads"], reports[1]["metrics"]["threads"],
        "the report does say how the run was executed"
    );
    // Everything else in the report is -j-invariant too.
    let mut a = reports[0]["metrics"].as_object().unwrap().clone();
    let mut b = reports[1]["metrics"].as_object().unwrap().clone();
    a.remove("threads");
    b.remove("threads");
    assert_eq!(a, b, "only `threads` may differ between -j 1 and -j 2");
}
