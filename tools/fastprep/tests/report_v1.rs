//! `--json` writes the `fastsuite.report.v1` wrapper (PLAN.md §1.1).
//!
//! fastprep used to write the raw `Stats` object here. It now writes the
//! documented wrapper, with the old raw object preserved under `legacy` for
//! one release, so `fasthelper concatenate` reads it through the canonical
//! v1 path instead of the legacy adapter.

use std::path::PathBuf;
use std::process::Command;

use serde_json::Value;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fastprep_v1_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A FASTQ record; the quality string is derived from the sequence length, so
/// a fixture edit can never produce a seq/qual mismatch.
fn rec(id: &str, seq: &str) -> String {
    format!("@{id}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
}

/// Four reads: one clean, one with a 5' polyA head (dropped), one with a 3'
/// polyT tail (trimmed, kept), one clean. Every counter is non-trivial.
fn four_reads() -> String {
    [
        rec("clean_1", "ACGTACGTACGTACGTACGTACGTACGTACGT"),
        rec(
            "polyA_head",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTACGTACGT",
        ),
        rec(
            "polyT_tail",
            "ACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
        rec("clean_2", "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"),
    ]
    .join("")
}

fn run(args: &[&str]) -> (String, String, i32) {
    let out = Command::new(fastprep())
        .args(args)
        .output()
        .expect("spawn fastprep");
    (
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

/// Depth-first walk of `metrics`, collecting the path of every value that is
/// an object or an array. The schema forbids nested structures here.
fn nested_metrics(value: &Value) -> Vec<String> {
    fn walk(v: &Value, path: &str, out: &mut Vec<String>) {
        match v {
            Value::Object(o) => {
                if !path.is_empty() {
                    out.push(path.to_string());
                }
                for (k, child) in o {
                    walk(child, &format!("{path}.{k}"), out);
                }
            }
            Value::Array(a) => {
                if !path.is_empty() {
                    out.push(path.to_string());
                }
                for (i, child) in a.iter().enumerate() {
                    walk(child, &format!("{path}[{i}]"), out);
                }
            }
            _ => {}
        }
    }
    let mut out = Vec::new();
    walk(value, "", &mut out);
    out
}

/// Single-end `--json`: the wrapper, its flat metrics, and the preserved
/// raw-stats object.
#[test]
fn single_end_json_is_the_v1_wrapper_with_flat_metrics() {
    let d = tmpdir("single");
    let fq = d.join("in.fastq");
    std::fs::write(&fq, four_reads()).unwrap();
    let json = d.join("report.json");
    let out = d.join("out.fastq");

    let (_, stderr, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--json",
        json.to_str().unwrap(),
        "--min-len",
        "0",
        "--no-adapter",
    ]);
    assert_eq!(code, 0, "fastprep must succeed: {stderr}");

    let doc: Value = serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(doc["schema"], "fastsuite.report.v1");
    assert_eq!(doc["tool"], "fastprep");
    assert_eq!(
        doc["tool_version"],
        env!("CARGO_PKG_VERSION"),
        "tool_version must be this crate's version"
    );
    assert_eq!(doc["mode"], "single");
    assert_eq!(doc["dry_run"], false);
    assert!(doc["command"].as_str().is_some_and(|c| !c.is_empty()));
    assert!(doc["inputs"]["r1"].as_str().is_some(), "{doc}");
    assert!(doc["outputs"]["r1"].as_str().is_some(), "{doc}");
    assert!(doc["elapsed_s"].is_number(), "{doc}");

    // Four reads in, and metrics are a flat numeric namespace.
    assert_eq!(doc["metrics"]["reads_in"], 4);
    assert!(
        nested_metrics(&doc["metrics"]).is_empty(),
        "metrics must stay flat: {:?}",
        nested_metrics(&doc["metrics"])
    );

    // The raw Stats object survives under `legacy`, with the old shape:
    // reads_in lives directly on it, as it did before the wrapper.
    let legacy = &doc["legacy"];
    assert_eq!(legacy["reads_in"], 4, "legacy must hold the old raw shape");
    assert_eq!(legacy["reads_kept"], doc["metrics"]["reads_out"]);
    assert!(legacy["bases_in"].is_number(), "{legacy}");

    std::fs::remove_dir_all(&d).unwrap();
}

/// Paired `--json`: `mode` is paired, reads_in sums both mates, and the
/// per-mate metrics are prefixed where they would otherwise be ambiguous.
#[test]
fn paired_json_wrapper_sums_mates_and_prefixes_ambiguous_metrics() {
    let d = tmpdir("paired");
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    std::fs::write(&r1, four_reads()).unwrap();
    std::fs::write(&r2, four_reads()).unwrap();
    let json = d.join("report.json");
    let o1 = d.join("o1.fastq");
    let o2 = d.join("o2.fastq");

    let (_, stderr, code) = run(&[
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-o",
        o1.to_str().unwrap(),
        "-p",
        o2.to_str().unwrap(),
        "--json",
        json.to_str().unwrap(),
        "--min-len",
        "0",
        "--no-adapter",
    ]);
    assert_eq!(code, 0, "fastprep must succeed: {stderr}");

    let doc: Value = serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(doc["schema"], "fastsuite.report.v1");
    assert_eq!(doc["tool"], "fastprep");
    assert_eq!(doc["mode"], "paired");
    assert!(doc["inputs"]["r2"].as_str().is_some(), "{doc}");
    assert!(doc["outputs"]["r2"].as_str().is_some(), "{doc}");

    // 4 reads per mate, 8 in total.
    assert_eq!(doc["metrics"]["reads_in"], 8);
    assert_eq!(doc["metrics"]["r1_reads_in"], 4);
    assert_eq!(doc["metrics"]["r2_reads_in"], 4);
    assert!(
        doc["metrics"]["pairs_out"].is_number(),
        "pair counters are top-level metrics: {doc}"
    );
    assert!(
        nested_metrics(&doc["metrics"]).is_empty(),
        "metrics must stay flat: {:?}",
        nested_metrics(&doc["metrics"])
    );

    // The legacy object keeps the paired shape, per mate.
    let legacy = &doc["legacy"];
    assert_eq!(legacy["r1"]["reads_in"], 4);
    assert_eq!(legacy["r2"]["reads_in"], 4);
    assert_eq!(legacy["pairs_out"], doc["metrics"]["pairs_out"]);

    std::fs::remove_dir_all(&d).unwrap();
}

/// `--dry-run` writes nothing at all, including the report (CLI-SPEC §1).
#[test]
fn dry_run_still_writes_no_json() {
    let d = tmpdir("dryrun");
    let fq = d.join("in.fastq");
    std::fs::write(&fq, four_reads()).unwrap();
    let json = d.join("report.json");
    let out = d.join("out.fastq");

    let (_, stderr, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--json",
        json.to_str().unwrap(),
        "--dry-run",
        "--min-len",
        "0",
        "--no-adapter",
    ]);
    assert_eq!(code, 0, "--dry-run must succeed: {stderr}");
    assert!(
        !json.exists(),
        "--dry-run must not write the --json report"
    );
    assert!(!out.exists(), "--dry-run must not write the output");

    std::fs::remove_dir_all(&d).unwrap();
}

/// The wrapper is written verbatim: `fasthelper concatenate` must read it
/// through the v1 path, which means the metrics it derives match the report.
#[test]
fn wrapper_carries_the_derived_retention_keys() {
    let d = tmpdir("retention");
    let fq = d.join("in.fastq");
    std::fs::write(&fq, four_reads()).unwrap();
    let json = d.join("report.json");
    let out = d.join("out.fastq");

    let (_, stderr, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--json",
        json.to_str().unwrap(),
        "--min-len",
        "0",
        "--no-adapter",
    ]);
    assert_eq!(code, 0, "fastprep must succeed: {stderr}");

    let doc: Value = serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    let m = &doc["metrics"];
    // The recommended shared keys, all present and consistent.
    assert_eq!(m["reads_in"], 4.0);
    assert_eq!(m["reads_out"], m["reads_in"].as_f64().unwrap() - m["reads_dropped"].as_f64().unwrap());
    assert!(m["bases_in"].is_number(), "{m}");
    assert!(m["bases_trimmed"].is_number(), "{m}");

    std::fs::remove_dir_all(&d).unwrap();
}
