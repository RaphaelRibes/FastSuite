//! M10 — the numeric knobs check their conversions.
//!
//! `--buffer-mb` was multiplied by 1 MiB unchecked: 2^44 MB wraps to a zero
//! budget, so the run "succeeded" with a threshold of 0. `read_num` bumped with
//! `checked_add(...).expect(...)`, i.e. a panic (exit 101) on a file with 2^32
//! reads.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_num_{tag}_{}", std::process::id()));
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

/// 2^44 MB is 2^64 bytes: it cannot be a budget, and must not become 0.
#[test]
fn an_overflowing_buffer_budget_is_refused() {
    let dir = tmpdir("wrap");
    let fq = fixture(&dir);
    let o = Command::new(bin())
        .arg("-o")
        .arg(dir.join("out"))
        .arg("--buffer-mb")
        .arg("17592186044416")
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        !o.status.success(),
        "a --buffer-mb that overflows must be refused, not wrapped to 0\nstderr:\n{stderr}"
    );
    assert_ne!(o.status.code(), Some(101), "stderr:\n{stderr}");
    assert!(
        stderr.contains("--buffer-mb") && stderr.contains("too large"),
        "the refusal must name the knob and the problem; stderr:\n{stderr}"
    );
    assert!(
        !dir.join("out").join("heatmap.json").exists(),
        "nothing may be written by a refused run"
    );
}

/// A large but representable budget still works (the check must not turn into
/// a cap on lawful values).
#[test]
fn a_large_representable_budget_still_runs() {
    let dir = tmpdir("big");
    let fq = fixture(&dir);
    let o = Command::new(bin())
        .arg("-o")
        .arg(dir.join("out"))
        .arg("--buffer-mb")
        .arg("4096")
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
}
