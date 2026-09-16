//! CLI-level tests for the toolbox surface.
//!
//! `concatenate` is the name; `report` stays as an alias so nothing that was
//! written against the earlier build breaks (docs/DECISIONS.md D-8).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasthelper() -> &'static str {
    env!("CARGO_BIN_EXE_fasthelper")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fasthelper_cli_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// One minimal fastsuite.report.v1 file, as fastcheck would write it.
fn write_report(dir: &Path, name: &str, reads: u64, gc: f64) -> PathBuf {
    let p = dir.join(name);
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastcheck",
        "tool_version": "0.1.0",
        "elapsed_s": 1.0,
        "inputs": {"files": [name]},
        "metrics": {"reads": reads, "bases": reads * 150, "gc_fraction": gc}
    });
    std::fs::write(&p, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    p
}

fn run(args: &[&str]) -> (String, String, i32) {
    let out = Command::new(fasthelper()).args(args).output().expect("spawn fasthelper");
    (
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

#[test]
fn concatenate_merges_reports_into_one_row_per_sample() {
    let d = tmpdir("concat");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "sampleA.json", 100, 0.41);
    write_report(&reports, "sampleB.json", 200, 0.52);
    let out = d.join("out");

    let (stdout, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "tsv,json",
        "--stdout",
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    assert!(stdout.contains("sampleA") && stdout.contains("sampleB"), "{stdout}");

    let tsv = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines.len(), 3, "header + one row per sample: {tsv}");
    assert!(tsv.contains("sampleA"), "{tsv}");
    assert!(tsv.contains("sampleB"), "{tsv}");
    assert!(stdout.contains("concatenated 2 report(s)"), "{stdout}");
}

#[test]
fn report_is_still_a_working_alias() {
    let d = tmpdir("alias");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "s.json", 10, 0.5);

    let a = d.join("via_concatenate");
    let b = d.join("via_report");
    let (_, e1, c1) = run(&["concatenate", reports.to_str().unwrap(), "-o", a.to_str().unwrap()]);
    let (_, e2, c2) = run(&["report", reports.to_str().unwrap(), "-o", b.to_str().unwrap()]);
    assert_eq!(c1, 0, "{e1}");
    assert_eq!(c2, 0, "the alias must still work: {e2}");

    let ta = std::fs::read_to_string(a.join("report.tsv")).unwrap();
    let tb = std::fs::read_to_string(b.join("report.tsv")).unwrap();
    assert_eq!(ta, tb, "both spellings must produce the same table");
}

#[test]
fn help_lists_the_toolbox() {
    let (stdout, _, code) = run(&["--help"]);
    assert_eq!(code, 0);
    for cmd in ["concatenate", "run", "doctor", "schema"] {
        assert!(stdout.contains(cmd), "{cmd} missing from:\n{stdout}");
    }
}
