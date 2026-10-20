//! Review 2026-10-21, items 1/5/6 — `--json` is the *report*, never the grid.
//!
//! Before this, `--json r.json` wrote the grid to `r.json` and the
//! `fastsuite.report.v1` wrapper to `r.report.json`: a suite-uniform caller
//! that asked every tool for `--json report.json` got a megabyte grid at exit
//! 0, and `-o` never received its `heatmap.json` at all. The contract (README,
//! CLI-SPEC §1, and every other tool's `--help`) is that `--json PATH` writes
//! the v1 report at PATH and the grid keeps its default
//! `<output-dir>/heatmap.json`.
//!
//! Also pinned here:
//!   * item 5 — `-o fresh --json fresh/r.json` used to fail because the report
//!     parent was validated before the output dir existed;
//!   * item 6 — the clobber guard used to ignore a `--json` target outside
//!     `-o`, so a second run overwrote it without `--force`.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_json_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fixture(dir: &Path, name: &str, read_num: &str) -> PathBuf {
    let fq = dir.join(name);
    let s = "A".repeat(60);
    let mut body = String::new();
    for i in 0..3 {
        body.push_str(&format!("@{read_num}{i}\n{s}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin()).args(args).output().unwrap()
}

fn parse(path: &Path) -> Value {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{} must exist: {e}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("{} must parse: {e}\n{text}", path.display()))
}

/// `--json PATH` writes the report at PATH; the grid stays in `-o`, and no
/// `<stem>.report.json` is invented next to either file.
#[test]
fn json_writes_the_v1_report_and_keeps_the_grid_in_the_output_dir() {
    let dir = tmpdir("report");
    let fq = fixture(&dir, "a_R1.fastq", "r");
    let out = dir.join("out");
    let report = dir.join("r.json");

    let o = run(&[
        "-o",
        out.to_str().unwrap(),
        "--json",
        report.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        o.status.success(),
        "run failed\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );

    // PATH is the v1 report, not the grid.
    let r = parse(&report);
    assert_eq!(r["schema"], "fastsuite.report.v1", "r.json must be the report, got {r}");
    assert_eq!(r["tool"], "fastmaps");
    assert!(r["metrics"]["total_reads"].as_f64().unwrap() == 3.0);

    // The grid is still in -o, under its historical name.
    let g = parse(&out.join("heatmap.json"));
    assert_eq!(g["version"], 2, "the grid JSON keeps its own shape: {g}");
    assert_eq!(g["size"], 61);
    assert_eq!(g["counts"].as_array().unwrap().len(), 61 * 61);

    // The wrapper's grid block names the file a consumer will find.
    assert_eq!(r["grid"]["format"], "heatmap.json", "{r}");

    // No PATH-stem wrapper anywhere.
    for stale in ["r.report.json", "heatmap.report.json"] {
        assert!(
            !dir.join(stale).exists(),
            "{stale} must not be created: --json PATH is the report itself"
        );
        assert!(!out.join(stale).exists(), "{stale} must not be created in -o");
    }
}

/// A `--corners` run still reports the real grid file name, and the corners
/// extraction lands in `-o` as before.
#[test]
fn corners_run_names_the_real_grid_file() {
    let dir = tmpdir("corners");
    let fq = fixture(&dir, "a_R1.fastq", "c");
    let out = dir.join("out");
    let report = dir.join("report.json");

    let o = run(&[
        "-o",
        out.to_str().unwrap(),
        "--corners",
        "0,0:5,5",
        "--json",
        report.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let r = parse(&report);
    assert_eq!(r["schema"], "fastsuite.report.v1");
    assert_eq!(
        r["grid"]["format"], "heatmap.json",
        "the grid block must name the file that exists: {r}"
    );
    assert!(out.join("heatmap.json").is_file(), "grid must be in -o");
    assert!(
        out.join("corners_0_0_5_5.fq").is_file(),
        "the corners file lands in -o"
    );
    assert!(!dir.join("report.report.json").exists());
}

/// Item 5: the report parent may be the output directory the same run creates.
#[test]
fn json_inside_a_fresh_output_dir_is_accepted() {
    let dir = tmpdir("fresh");
    let fq = fixture(&dir, "a_R1.fastq", "f");
    let out = dir.join("fresh");
    let report = out.join("r.json");

    let o = run(&[
        "-o",
        out.to_str().unwrap(),
        "--json",
        report.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        o.status.success(),
        "-o fresh --json fresh/r.json must work\nstderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert_eq!(parse(&report)["schema"], "fastsuite.report.v1");
    assert!(out.join("heatmap.json").is_file());
}

/// Item 6: a report target outside `-o` is guarded like every other output.
#[test]
fn json_target_outside_the_output_dir_is_guarded() {
    let dir = tmpdir("outside");
    let fq = fixture(&dir, "a_R1.fastq", "g");
    let gb = dir.join("gB");
    std::fs::create_dir_all(&gb).unwrap();
    let report = gb.join("x.json");

    let first = run(&[
        "-o",
        dir.join("gA").to_str().unwrap(),
        "--json",
        report.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    assert!(
        first.status.success(),
        "run 1 failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(report.is_file());

    // Second run into a *different* -o, same report target, no --force.
    let sentinel: &[u8] = b"sentinel: must survive a refused run\n";
    std::fs::write(&report, sentinel).unwrap();
    let second = run(&[
        "-o",
        dir.join("gA2").to_str().unwrap(),
        "--json",
        report.to_str().unwrap(),
        fq.to_str().unwrap(),
    ]);
    let stderr = String::from_utf8_lossy(&second.stderr).to_string();
    assert!(
        !second.status.success(),
        "run 2 must refuse to clobber {}\nstderr: {stderr}",
        report.display()
    );
    assert!(
        stderr.contains("x.json"),
        "the refusal must name the report file; stderr:\n{stderr}"
    );
    assert_eq!(
        std::fs::read(&report).unwrap(),
        sentinel,
        "the report was overwritten despite the refusal"
    );

    // --force rewrites it.
    let third = run(&[
        "-o",
        dir.join("gA2").to_str().unwrap(),
        "--json",
        report.to_str().unwrap(),
        "--force",
        fq.to_str().unwrap(),
    ]);
    assert!(
        third.status.success(),
        "--force must overwrite: {}",
        String::from_utf8_lossy(&third.stderr)
    );
    let r = parse(&report);
    assert_eq!(r["schema"], "fastsuite.report.v1", "--force must rewrite the report");
    assert_ne!(std::fs::read(&report).unwrap(), sentinel);
}
