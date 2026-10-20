//! P0 — mixed read lengths silently collapse the corner fractions.
//!
//! `size` is `max_seq_len + 1` and the corner window is `size / 50`, so the
//! grid is always sized for the *longest* read. A read shorter than that can
//! never reach a corner pixel, and leaves the fraction's denominator: ten
//! 150 bp all-T reads report `polyT_frac = 1.0` on their own, and the same ten
//! reads plus one 200 bp read report `polyT_frac = 0.0` — with no warning and
//! nothing in either output to notice it from.
//!
//! This fix does not change the corner definition. It publishes the lengths the
//! window was derived from (`min_seq_len`, `max_seq_len`, `corner_window`, in
//! both the grid JSON and the v1 wrapper) and warns on stderr when the spread
//! exceeds the window, so a length-biased run can no longer pass unnoticed.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

/// `n` reads of `len` bp, all `base`.
fn reads(base: &str, len: usize, n: usize, tag: &str) -> String {
    let s = base.repeat(len);
    let mut body = String::new();
    for i in 0..n {
        body.push_str(&format!("@{tag}{i}\n{s}\n+\n{}\n", "I".repeat(len)));
    }
    body
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_span_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run the binary on `body` with the grid redirected to `<dir>/grid.json`, so
/// the v1 wrapper is written too. Returns the run's dir and stderr.
fn run(tag: &str, body: &str) -> (PathBuf, String) {
    let dir = tmpdir(tag);
    let fq = dir.join("reads_R1.fastq");
    std::fs::write(&fq, body).unwrap();
    let out = dir.join("out");
    let grid = dir.join("grid.json");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--json")
        .arg(&grid)
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "run {tag} failed\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    (dir, String::from_utf8_lossy(&o.stderr).to_string())
}

fn parse(path: &Path) -> Value {
    let text = std::fs::read_to_string(path).unwrap();
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{} must parse: {e}\n{text}", path.display()))
}

/// A uniform-length run: the window is not biased, so nothing to warn about —
/// but the lengths it was derived from must still be published.
#[test]
fn uniform_lengths_publish_lengths_and_do_not_warn() {
    let (dir, stderr) = run("uniform", &reads("T", 150, 10, "u"));

    assert!(
        !stderr.contains("length-biased"),
        "a uniform-length run must not warn; stderr:\n{stderr}"
    );

    let m = parse(&dir.join("grid.report.json"));
    let m = &m["metrics"];
    assert_eq!(m["min_seq_len"].as_f64().unwrap(), 150.0, "{m}");
    assert_eq!(m["max_seq_len"].as_f64().unwrap(), 150.0, "{m}");
    // size = 151 → 151 / 50 = 3.
    assert_eq!(m["corner_window"].as_f64().unwrap(), 3.0, "{m}");
    // The ten all-T reads do reach the polyT corner here.
    assert!(
        (m["polyT_frac"].as_f64().unwrap() - 1.0).abs() < 1e-9,
        "uniform run: polyT_frac must be 1.0, got {}",
        m["polyT_frac"]
    );

    let g = parse(&dir.join("grid.json"));
    assert_eq!(g["min_seq_len"], 150, "grid JSON carries min_seq_len: {g}");
    assert_eq!(g["max_seq_len"], 150, "grid JSON carries max_seq_len: {g}");
    assert_eq!(g["corner_window"], 3, "grid JSON carries corner_window: {g}");
}

/// One 200 bp read in a 150 bp run: the corner fractions now describe only the
/// longest read, which is a silent trap. The run must say so, and both outputs
/// must carry the lengths that make the trap visible after the fact.
#[test]
fn mixed_lengths_warn_and_publish_the_spread() {
    let body = reads("T", 150, 10, "t") + &reads("A", 200, 1, "a");
    let (dir, stderr) = run("mixed", &body);

    assert!(
        stderr.contains("length-biased"),
        "a mixed-length run must warn about the corner bias; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("150..200"),
        "the warning must name the length spread; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("spread 50"),
        "the warning must state the spread against the window; stderr:\n{stderr}"
    );

    let m = parse(&dir.join("grid.report.json"));
    let m = &m["metrics"];
    assert_eq!(m["min_seq_len"].as_f64().unwrap(), 150.0, "{m}");
    assert_eq!(m["max_seq_len"].as_f64().unwrap(), 200.0, "{m}");
    // size = 201 → 201 / 50 = 4.
    assert_eq!(m["corner_window"].as_f64().unwrap(), 4.0, "{m}");

    // The bias itself, unchanged by this fix and asserted so it stays visible:
    // none of the ten 150 bp reads can reach the 200 bp corner window, so the
    // polyT corner is empty and the lone 200 bp all-A read takes polyA.
    assert!(
        m["polyT_frac"].as_f64().unwrap().abs() < 1e-9,
        "mixed run: the 150 bp reads fall outside the corner window, got polyT_frac = {}",
        m["polyT_frac"]
    );
    // (The wrapper renders 6 decimals, so the tolerance is 1e-6, not 1e-9.)
    assert!(
        (m["polyA_frac"].as_f64().unwrap() - 1.0 / 11.0).abs() < 1e-6,
        "mixed run: polyA_frac is 1 of 11 reads, got {}",
        m["polyA_frac"]
    );

    let g = parse(&dir.join("grid.json"));
    assert_eq!(g["min_seq_len"], 150, "grid JSON carries min_seq_len: {g}");
    assert_eq!(g["max_seq_len"], 200, "grid JSON carries max_seq_len: {g}");
    assert_eq!(g["corner_window"], 4, "grid JSON carries corner_window: {g}");
}
