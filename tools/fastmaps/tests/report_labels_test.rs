//! M3 — a 3' polyX tail is invisible to the corner fractions; say so.
//! M4 — publish the preview quota's real granularity and measure the maximum.
//! Cross-cutting — `metrics.threads` carries the resolved worker count.
//!
//! M3: 200 reads of 131 random bases + 20 trailing G give `r2_polyG_frac 0.0`
//! with no warning at all. The corner fractions measure *whole-read*
//! composition, so a polyG tail on a random body lands mid-plane and never
//! reaches a corner. The corner definition is not changed here; the run
//! publishes a `notes` entry and a stderr warning instead.
//!
//! M4: `preview_reads_per_pixel` = 5 reads like "5 previews per pixel", but
//! the quota is per `(pixel, kind)`: R1 and R2 each get 5, so one pixel holds
//! 10 (15 with an Unknown file). The label is published as
//! `preview_reads_per_pixel_per_kind` and the run's measured maximum as
//! `max_previews_per_pixel`.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_lbl_{tag}_{}", std::process::id()));
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

struct Run {
    out: PathBuf,
    report: PathBuf,
    stderr: String,
    ok: bool,
}

/// One binary run in `dir`, with its outputs in `dir/<tag>-out` and its report
/// in `dir/<tag>-report.json`. `dir` also holds the fixtures, so it is the
/// caller's, not this helper's, to create.
fn run(dir: &Path, tag: &str, files: &[PathBuf], extra: &[&str]) -> Run {
    let out = dir.join(format!("{tag}-out"));
    let report = dir.join(format!("{tag}-report.json"));
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--json")
        .arg(&report)
        .args(extra)
        .args(files)
        .output()
        .unwrap();
    Run {
        out,
        report,
        stderr: String::from_utf8_lossy(&o.stderr).to_string(),
        ok: o.status.success(),
    }
}

/// `n` reads of 131 random bases + a `tail_len` homopolymer of `base`.
fn tailed_library(dir: &Path, name: &str, n: usize, base: char, tail_len: usize) -> PathBuf {
    let fq = dir.join(name);
    let mut body = String::with_capacity(n * 320);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in 0..n {
        let mut seq = String::with_capacity(151);
        for _ in 0..131 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            seq.push(match state % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            });
        }
        for _ in 0..tail_len {
            seq.push(base);
        }
        body.push_str(&format!("@t{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

/// The reported blind spot: a 3' polyG tail with a 0.0 polyG corner fraction.
#[test]
fn a_terminal_poly_g_tail_is_noted_and_warned() {
    let dir = tmpdir("tailg");
    let fq = tailed_library(&dir, "tail_R2.fastq", 200, 'G', 20);
    let r = run(&dir, "tailg", &[fq], &[]);
    assert!(r.ok, "stderr:\n{}", r.stderr);

    let v = parse(&r.report);
    let m = &v["metrics"];
    assert!(
        m["r2_polyG_frac"].as_f64().unwrap().abs() < 1e-9,
        "the fixture must keep the reported blind spot (r2_polyG_frac 0), got {}",
        m["r2_polyG_frac"]
    );

    let notes = v["notes"]
        .as_array()
        .unwrap_or_else(|| panic!("the report must carry a `notes` array: {v}"));
    let hit = notes
        .iter()
        .filter_map(|n| n.as_str())
        .find(|n| n.contains("polyG"))
        .unwrap_or_else(|| panic!("no polyG note in {notes:?}"));
    assert!(
        hit.contains("200") && hit.contains("3'"),
        "the note must name the reads and the 3' end: {hit}"
    );

    assert!(
        r.stderr.contains("polyG") && r.stderr.contains("3'"),
        "the run must warn on stderr too, not only in the report; stderr:\n{}",
        r.stderr
    );
    assert!(
        r.stderr.contains("fastmaps: warning:"),
        "the warning must be a warning; stderr:\n{}",
        r.stderr
    );
}

/// A library whose 3' ends are random has nothing to report.
#[test]
fn a_library_without_tails_produces_no_notes() {
    let dir = tmpdir("clean");
    let fq = tailed_library(&dir, "clean_R1.fastq", 200, 'A', 0);
    let r = run(&dir, "clean", &[fq], &[]);
    assert!(r.ok, "stderr:\n{}", r.stderr);
    let v = parse(&r.report);
    let notes = v["notes"].as_array().expect("notes array");
    assert!(notes.is_empty(), "unexpected notes: {notes:?}");
    assert!(
        !r.stderr.contains("3' poly"),
        "stderr:\n{}",
        r.stderr
    );
}

/// A tail shorter than the detector's window is not called: a 6 bp run at the
/// 3' end of 200 reads is noise, not contamination.
#[test]
fn a_short_terminal_run_is_not_called_a_tail() {
    let dir = tmpdir("short");
    let fq = tailed_library(&dir, "s_R1.fastq", 200, 'G', 6);
    let r = run(&dir, "short", &[fq], &[]);
    assert!(r.ok, "stderr:\n{}", r.stderr);
    let v = parse(&r.report);
    assert!(v["notes"].as_array().unwrap().is_empty(), "{v}");
}

// ── M4 ──────────────────────────────────────────────────────────────────────

fn all_a(dir: &Path, name: &str, n: usize) -> PathBuf {
    let fq = dir.join(name);
    let s = "A".repeat(60);
    let mut body = String::new();
    for i in 0..n {
        body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

/// R1 + R2 at one pixel: 5 per (pixel, kind) = 10 in the pixel.
#[test]
fn the_preview_quota_is_published_per_pixel_and_kind() {
    let dir = tmpdir("pv2");
    let r1 = all_a(&dir, "p_R1.fastq", 20);
    let r2 = all_a(&dir, "p_R2.fastq", 20);
    let r = run(&dir, "pv2", &[r1, r2], &["--preview-reads", "5"]);
    assert!(r.ok, "stderr:\n{}", r.stderr);

    let g = parse(&r.out.join("heatmap.json"));
    assert_eq!(g["preview_reads_per_pixel"], 5, "{g}");
    assert_eq!(
        g["preview_reads_per_pixel_per_kind"], 5,
        "the grid must publish the quota's real granularity: {g}"
    );
    assert_eq!(
        g["max_previews_per_pixel"], 10,
        "R1 and R2 each keep 5 at one pixel, so the run's maximum is 10: {g}"
    );

    let m = parse(&r.report);
    let m = &m["metrics"];
    assert_eq!(m["preview_reads_per_pixel"].as_f64().unwrap(), 5.0, "{m}");
    assert_eq!(
        m["preview_reads_per_pixel_per_kind"].as_f64().unwrap(),
        5.0,
        "{m}"
    );
    assert_eq!(
        m["max_previews_per_pixel"].as_f64().unwrap(),
        10.0,
        "{m}"
    );
}

/// A third, Unknown file gets its own quota too: 15 at one pixel.
#[test]
fn an_unknown_file_gets_its_own_preview_quota() {
    let dir = tmpdir("pv3");
    let files = [
        all_a(&dir, "q_R1.fastq", 20),
        all_a(&dir, "q_R2.fastq", 20),
        all_a(&dir, "q_extra.fastq", 20),
    ];
    let r = run(&dir, "pv3", &files, &["--preview-reads", "5"]);
    assert!(r.ok, "stderr:\n{}", r.stderr);
    let g = parse(&r.out.join("heatmap.json"));
    assert_eq!(g["max_previews_per_pixel"], 15, "{g}");
    assert_eq!(g["preview_reads_per_pixel"], 5, "{g}");
}

/// `--preview-reads 0` writes no previews at all: the measured maximum is 0.
#[test]
fn no_previews_means_a_zero_maximum() {
    let dir = tmpdir("pv0");
    let fq = all_a(&dir, "z_R1.fastq", 5);
    let r = run(&dir, "pv0", &[fq], &["--preview-reads", "0"]);
    assert!(r.ok, "stderr:\n{}", r.stderr);
    let g = parse(&r.out.join("heatmap.json"));
    assert_eq!(g["max_previews_per_pixel"], 0, "{g}");
}

/// The CLI says what the knob means.
#[test]
fn help_states_the_per_kind_semantic() {
    let o = Command::new(bin()).arg("--help").output().unwrap();
    let help = String::from_utf8_lossy(&o.stdout).to_string();
    assert!(
        help.contains("per (pixel, kind)"),
        "--help must state the preview quota's granularity:\n{help}"
    );
}

// ── cross-cutting ───────────────────────────────────────────────────────────

#[test]
fn metrics_carry_the_resolved_thread_count() {
    let dir = tmpdir("threads");
    let r1 = all_a(&dir, "th_R1.fastq", 3);
    let r2 = all_a(&dir, "th_R2.fastq", 3);

    let two = run(&dir, "threads2", &[r1.clone(), r2.clone()], &["-j", "2"]);
    assert!(two.ok, "stderr:\n{}", two.stderr);
    assert_eq!(
        parse(&two.report)["metrics"]["threads"].as_f64().unwrap(),
        2.0,
        "-j 2 must be published as threads=2"
    );

    let one = run(&dir, "threads1", &[r1, r2], &["-j", "1"]);
    assert!(one.ok, "stderr:\n{}", one.stderr);
    assert_eq!(
        parse(&one.report)["metrics"]["threads"].as_f64().unwrap(),
        1.0
    );
}
