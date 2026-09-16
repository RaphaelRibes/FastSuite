//! D-6 policy (docs/DECISIONS.md): when a detection returns nothing, take the
//! safe action, label it, and exit 0 — `--strict` restores the refusal.
//!
//! This is the difference between a good sample failing a 700-sample array job
//! and a good sample being reported honestly.

use std::path::PathBuf;
use std::process::Command;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_policy_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Adapter-free reads: no fragment of the built-in catalogue appears here, so
/// auto-detection legitimately finds nothing.
fn write_clean_fastq(dir: &std::path::Path) -> PathBuf {
    let p = dir.join("clean.fastq");
    let mut body = String::new();
    // Deterministic pseudo-random ACGT, 60 bp each, 200 reads.
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    for i in 0..200 {
        let mut seq = String::with_capacity(60);
        for _ in 0..60 {
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
        body.push_str(&format!("@read{i}\n{seq}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&p, body).unwrap();
    p
}

fn run(dir: &std::path::Path, args: &[&str]) -> (String, i32) {
    let out = Command::new(fastprep())
        .args(args)
        .current_dir(dir)
        .output()
        .expect("spawn fastprep");
    (
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

#[test]
fn no_adapter_detected_continues_with_a_labelled_assumption() {
    let d = tmpdir("continue");
    let fq = write_clean_fastq(&d);
    let out = d.join("out.fastq");
    let (stderr, code) = run(
        &d,
        &[
            "-1",
            fq.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--min-len",
            "0",
        ],
    );
    assert_eq!(code, 0, "a clean sample must not fail the run: {stderr}");
    assert!(
        stderr.contains("assumption") && stderr.contains("NO adapter trimming"),
        "the assumption must be stated, not silent: {stderr}"
    );
    assert!(out.exists(), "output is still produced");
}

#[test]
fn strict_turns_the_assumption_into_a_failure() {
    let d = tmpdir("strict");
    let fq = write_clean_fastq(&d);
    let out = d.join("out.fastq");
    let (stderr, code) = run(
        &d,
        &[
            "-1",
            fq.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--min-len",
            "0",
            "--strict",
        ],
    );
    assert_ne!(code, 0, "--strict must refuse: {stderr}");
    assert!(
        stderr.contains("--strict") && stderr.contains("--adapter-preset"),
        "the refusal must say how to fix it: {stderr}"
    );
    assert!(!out.exists(), "nothing is written when refusing");
}
