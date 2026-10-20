//! Review 2026-10-21, item 3 — a record whose qual length != seq length.
//!
//! `encode_preview_row` writes `varint(seq_len)` and then seq+qual, so the
//! decoder trusts `slen` for both. A short qual string made the decoder slice
//! past the row buffer (panic, exit 101); a long one read a garbage length out
//! of the payload and panicked on a second slice. Worse, the same file exited
//! 0 under `--preview-reads 0`, because no preview row was ever decoded.
//!
//! The record is malformed FASTQ: it must be refused at read time, with the
//! file and the read number, whatever the preview/data mode.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_rec_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// `seq_len` seq bases with `qual_len` quality bytes, one read.
fn fixture(dir: &Path, seq_len: usize, qual_len: usize) -> PathBuf {
    let fq = dir.join("bad_R1.fastq");
    let body = format!("@r1\n{}\n+\n{}\n", "A".repeat(seq_len), "I".repeat(qual_len));
    std::fs::write(&fq, body).unwrap();
    fq
}

/// Run and require: non-zero exit, a message naming the file and read 1, and
/// no panic (a panic is exit 101 and the decoder's fingerprint, not a refusal).
fn assert_refused(tag: &str, seq_len: usize, qual_len: usize, extra: &[&str]) {
    let dir = tmpdir(tag);
    let fq = fixture(&dir, seq_len, qual_len);
    let out = dir.join("out");
    let mut c = Command::new(bin());
    c.arg("-o").arg(&out).args(extra).arg(&fq);
    let o = c.output().unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    let stdout = String::from_utf8_lossy(&o.stdout).to_string();

    assert!(
        !o.status.success(),
        "{tag}: a qual/seq length mismatch must be refused (seq {seq_len}, qual {qual_len})\nstdout: {stdout}\nstderr: {stderr}"
    );
    assert_ne!(
        o.status.code(),
        Some(101),
        "{tag}: must be a refusal, not a panic\nstderr: {stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "{tag}: the decoder must never panic; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("bad_R1.fastq"),
        "{tag}: the error must name the file; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("read 1"),
        "{tag}: the error must name the read; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains(&seq_len.to_string()) && stderr.contains(&qual_len.to_string()),
        "{tag}: the error must state both lengths; stderr:\n{stderr}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn short_qual_is_refused_with_file_and_read_number() {
    assert_refused("short", 60, 10, &[]);
}

#[test]
fn long_qual_is_refused_with_file_and_read_number() {
    assert_refused("long", 60, 200, &[]);
}

/// The mismatch must not depend on the preview stream existing.
#[test]
fn short_qual_is_refused_with_previews_diabled() {
    assert_refused("short_p0", 60, 10, &["--preview-reads", "0"]);
}

/// Nor on the per-read index stream: `--data-mode off` removes the render
/// pass, not the FASTQ validation.
#[test]
fn short_qual_is_refused_with_the_index_stream_off() {
    assert_refused("short_off", 60, 10, &["--data-mode", "off"]);
}

#[test]
fn long_qual_is_refused_with_previews_diabled() {
    assert_refused("long_p0", 60, 200, &["--preview-reads", "0"]);
}
