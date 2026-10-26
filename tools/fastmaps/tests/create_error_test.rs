//! M8 — an unwritable output path is an error, not a panic.
//!
//! `File::create(...).unwrap()` at `heatmap.png`, `heatmap.svg` and
//! `heatmap_3d.html` (and the spiller's row writes) aborted with exit 101 and
//! a panic message after a partial run. `--force` lets the run reach those
//! calls with a directory sitting where the file should go.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_io_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fixture(dir: &Path) -> PathBuf {
    let fq = dir.join("a_R1.fastq");
    let s = "A".repeat(60);
    let mut body = String::new();
    for i in 0..3 {
        body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

fn assert_refused_not_panicked(tag: &str, victim: &str, extra: &[&str]) {
    let dir = tmpdir(tag);
    let fq = fixture(&dir);
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    // A directory where the output file belongs.
    std::fs::create_dir(out.join(victim)).unwrap();

    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--force")
        .args(extra)
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();

    assert!(
        !o.status.success(),
        "{tag}: writing {victim} over a directory must fail\nstdout:\n{}\nstderr:\n{stderr}",
        String::from_utf8_lossy(&o.stdout)
    );
    assert_ne!(
        o.status.code(),
        Some(101),
        "{tag}: must be a refusal, not a panic\nstderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "{tag}: no panic in the output; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains(victim),
        "{tag}: the error must name {victim}; stderr:\n{stderr}"
    );
}

#[test]
fn png_over_a_directory_is_an_error() {
    assert_refused_not_panicked("png", "heatmap.png", &[]);
}

#[test]
fn html_over_a_directory_is_an_error() {
    assert_refused_not_panicked("html", "heatmap_3d.html", &[]);
}

#[test]
fn svg_over_a_directory_is_an_error() {
    assert_refused_not_panicked("svg", "heatmap.svg", &["--svg"]);
}

#[test]
fn data_over_a_directory_is_an_error() {
    // format.rs already returns an io::Result here; this pins that it stays
    // that way after the other writers were converted.
    assert_refused_not_panicked("data", "heatmap.data", &[]);
}

#[test]
fn grid_json_over_a_directory_is_an_error() {
    assert_refused_not_panicked("grid", "heatmap.json", &[]);
}
