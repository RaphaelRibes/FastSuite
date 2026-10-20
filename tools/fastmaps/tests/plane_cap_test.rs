//! Review 2026-10-21, item 4 — the output plane is O(max_read_len²).
//!
//! Three dense `u64` grids of `size²` cells (plus their JSON rendering) are
//! materialised at the end of every run, whatever `--buffer-mb` says: 3 reads
//! of 4 000 bp peaked at 976 MB with `--buffer-mb 1` *and* 8, 8 000 bp at
//! 3 788 MB. `--buffer-mb` sizes the read-index buffer and cannot bound the
//! plane, so the plane needs its own published ceiling and a refusal that says
//! so. Silent OOM is the one unacceptable outcome.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_plane_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn long_read_fixture(dir: &PathBuf, len: usize) -> PathBuf {
    let fq = dir.join("long_R1.fastq");
    let s = "ACGT".repeat(len / 4);
    let mut body = String::new();
    for i in 0..3 {
        body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(s.len())));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

/// A plane over the cap is refused with a message naming the read length, the
/// cap and the knob that does *not* move it — and nothing is written.
#[test]
fn oversize_plane_is_refused_with_a_message() {
    let dir = tmpdir("4000");
    let fq = long_read_fixture(&dir, 4000);
    let out = dir.join("out");

    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--buffer-mb")
        .arg("1")
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();

    assert!(
        !o.status.success(),
        "a 4000 bp plane must be refused, not attempted\nstderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "the refusal must be a message, not a panic:\n{stderr}"
    );
    assert!(
        stderr.contains("4000"),
        "the refusal must name the read length; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("plane"),
        "the refusal must say what is too big; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("--buffer-mb"),
        "the refusal must state that --buffer-mb cannot bound the plane; stderr:\n{stderr}"
    );
    // The documented maximum read length (see --help below).
    assert!(
        stderr.contains("2047"),
        "the refusal must name the usable read-length limit; stderr:\n{stderr}"
    );
    // Nothing half-written: the check precedes the first output file.
    assert!(
        !out.join("heatmap.json").exists() && !out.join("heatmap.data").exists(),
        "a refused run must not leave outputs behind"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The cap must be discoverable from the tool's own help, not only from an
/// error message.
#[test]
fn the_plane_cap_is_documented_in_help() {
    let o = Command::new(bin()).arg("--help").output().unwrap();
    assert!(o.status.success());
    let help = String::from_utf8_lossy(&o.stdout).to_string();
    assert!(
        help.contains("2047"),
        "--help must document the maximum read length the plane allows:\n{help}"
    );
    assert!(
        help.contains("plane"),
        "--help must name the plane limit:\n{help}"
    );
}

/// The same read length through the shipped defaults, no --buffer-mb flag:
/// still a refusal, so the limit is the plane's and not the buffer's.
#[test]
fn default_buffer_mb_refuses_the_same_plane() {
    let dir = tmpdir("4000_default");
    let fq = long_read_fixture(&dir, 4000);
    let o = Command::new(bin())
        .arg("-o")
        .arg(dir.join("out"))
        .arg(&fq)
        .output()
        .unwrap();
    assert!(!o.status.success(), "the plane cap must not depend on --buffer-mb");
    let _ = std::fs::remove_dir_all(&dir);
}

/// A short-read run — the tool's actual domain — is unaffected.
#[test]
fn short_reads_are_unaffected_by_the_cap() {
    let dir = tmpdir("150");
    let fq = long_read_fixture(&dir, 150 + 146); // 296 bp, 4-base multiple + tail
    let out = dir.join("out");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "a 296 bp run is inside the cap and must succeed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(out.join("heatmap.json").is_file());
    let _ = std::fs::remove_dir_all(&dir);
}
