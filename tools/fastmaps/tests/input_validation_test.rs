//! M7 — an incomplete final FASTQ record must be refused, like a truncated gz.
//! M11 — FASTA input has no qualities; refuse it instead of reporting on it.
//!
//! M7: the plain reader ended the loop on any short read, so 650 headers with
//! 649 complete records exited 0 and produced a full report for 649 reads —
//! while the same truncation inside a `.gz` member exits 1. A record that
//! starts and does not finish is malformed input, not EOF.
//!
//! M11: DECISIONS.md D-4 says fastmaps needs qualities and refuses FASTA.
//! It used to read `>name` / bases / `>name2` / bases as if it were FASTQ and
//! exit 0 with a full report.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_input_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn record(i: usize) -> String {
    let s = "A".repeat(60);
    format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(60))
}

fn stderr_of(o: &std::process::Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

fn run(dir: &Path, fq: &Path) -> std::process::Output {
    Command::new(bin())
        .arg("-o")
        .arg(dir.join("out"))
        .arg(fq)
        .output()
        .unwrap()
}

/// 649 complete records, then a header with no body.
#[test]
fn a_truncated_final_record_is_refused() {
    let dir = tmpdir("trunc");
    let fq = dir.join("cut_R1.fastq");
    let mut body = String::new();
    for i in 0..649 {
        body.push_str(&record(i));
    }
    body.push_str("@r649\nACGTACGT\n"); // sequence read, then EOF
    std::fs::write(&fq, body).unwrap();

    let o = run(&dir, &fq);
    let stderr = stderr_of(&o);
    assert!(
        !o.status.success(),
        "an incomplete final record must be refused, not silently dropped\nstdout:\n{}\nstderr:\n{stderr}",
        String::from_utf8_lossy(&o.stdout)
    );
    assert_ne!(
        o.status.code(),
        Some(101),
        "a refusal, not a panic; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("cut_R1.fastq"),
        "the error must name the file; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("650"),
        "the error must name the record (read 650); stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("truncated"),
        "the error must say what is wrong; stderr:\n{stderr}"
    );
    assert!(
        !dir.join("out").join("heatmap.json").exists(),
        "a refused run must not write a full report"
    );
}

/// The same file with the last record completed is fine, and so is a trailing
/// blank line (blank separators between records stay tolerated).
#[test]
fn complete_input_and_trailing_blank_lines_still_run() {
    let dir = tmpdir("complete");
    let fq = dir.join("ok_R1.fastq");
    let mut body = String::new();
    for i in 0..10 {
        body.push_str(&record(i));
    }
    body.push('\n');
    std::fs::write(&fq, body).unwrap();
    let o = run(&dir, &fq);
    assert!(
        o.status.success(),
        "stderr:\n{}",
        stderr_of(&o)
    );
    assert!(dir.join("out").join("heatmap.json").is_file());
}

/// A FASTA file: `>header` then bases, no `+`, no qualities.
#[test]
fn fasta_input_is_refused() {
    let dir = tmpdir("fasta");
    let fa = dir.join("ref.fasta");
    std::fs::write(
        &fa,
        ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAA\n",
    )
    .unwrap();
    let o = run(&dir, &fa);
    let stderr = stderr_of(&o);
    assert!(
        !o.status.success(),
        "FASTA input must be refused (D-4: fastmaps needs qualities)\nstdout:\n{}\nstderr:\n{stderr}",
        String::from_utf8_lossy(&o.stdout)
    );
    assert_ne!(o.status.code(), Some(101), "stderr:\n{stderr}");
    assert!(
        stderr.contains("FASTA"),
        "the refusal must say FASTA; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("ref.fasta"),
        "the refusal must name the file; stderr:\n{stderr}"
    );
    assert!(
        !dir.join("out").join("heatmap.json").exists(),
        "a refused run must not write a report"
    );
}

/// A `.gz`-named FASTA is caught by content too, not by extension.
#[test]
fn fasta_content_in_a_fastq_named_file_is_refused() {
    let dir = tmpdir("fasta_named");
    let fa = dir.join("looks_like_R1.fastq");
    std::fs::write(&fa, ">chr1\nACGTACGTACGT\n").unwrap();
    let o = run(&dir, &fa);
    assert!(!o.status.success(), "stderr:\n{}", stderr_of(&o));
    assert!(stderr_of(&o).contains("FASTA"), "{}", stderr_of(&o));
}
