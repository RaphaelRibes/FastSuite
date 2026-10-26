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

/// M6 — the refusal must fire on the read that needs the plane, not after the
/// whole input has been scanned and spilled. `[1/1]` is the progress line the
/// worker prints only when it finishes the file, so its absence proves the run
/// stopped mid-file; the temp dir must be gone too, not left full of rows.
#[test]
fn an_oversize_read_is_refused_before_the_pass_completes() {
    let dir = tmpdir("early");
    let fq = dir.join("mixed_R1.fastq");
    let mut body = String::new();
    for i in 0..650 {
        let s = "A".repeat(60);
        body.push_str(&format!("@s{i}\n{s}\n+\n{}\n", "I".repeat(60)));
    }
    // The oversize read is last: the whole file used to be read and spilled
    // before `write_format` looked at the plane.
    let s = "ACGT".repeat(750); // 3000 bp
    body.push_str(&format!("@long\n{s}\n+\n{}\n", "I".repeat(s.len())));
    std::fs::write(&fq, body).unwrap();

    let out = dir.join("out");
    let temp = dir.join("temp");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--temp-dir")
        .arg(&temp)
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(!o.status.success(), "stderr:\n{stderr}");
    assert!(
        stderr.contains("3000") && stderr.contains("plane") && stderr.contains("--buffer-mb"),
        "the refusal must still explain the plane cap; stderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("[1/1]"),
        "the run must stop when the 3000 bp read is seen, not after the file is done; stderr:\n{stderr}"
    );
    assert!(
        !out.join("heatmap.json").exists() && !out.join("heatmap.data").exists(),
        "a refused run must not write outputs"
    );
    // No spill rows left behind at the top of --temp-dir.
    if let Ok(rd) = std::fs::read_dir(&temp) {
        for e in rd {
            let name = e.unwrap().file_name().to_string_lossy().to_string();
            assert!(
                !name.starts_with("heatmap_row_") && !name.starts_with("preview_row_"),
                "a refused run left a spill row behind: {name}"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
