use assert_cmd::Command;
use assert_cmd::cargo;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;
use std::io::Write;

fn cmd() -> Command {
    cargo::cargo_bin_cmd!("fastdedup")
}

/// Two reads, second identical to the first: 2 in, 1 out.
const FASTQ: &[u8] = b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGT\n+\n!!!!!!!!\n";

#[test]
fn pre_count_with_stdin_is_refused() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .arg("--pre-count")
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("--pre-count"));
}

#[test]
fn both_inputs_as_stdin_is_refused() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-2")
        .arg("-")
        .arg("-o")
        .arg(temp_dir.child("out1.fastq").path())
        .arg("-p")
        .arg(temp_dir.child("out2.fastq").path())
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be stdin"));
}

#[test]
fn abundance_from_stdin_requires_temp_dir() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .arg("-1")
        .arg("-")
        .arg("--abundance")
        .arg("-o")
        .arg(temp_dir.child("out.fastq").path())
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("--temp-dir"));
}

#[test]
fn dry_run_from_stdin() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .arg("-v")
        .write_stdin(FASTQ)
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 2"))
        .stdout(predicate::str::contains("Duplicates removed: 50.00%"));
}

#[test]
fn empty_stdin_is_a_clean_refusal() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .write_stdin(b"")
        .assert()
        .failure()
        .stderr(predicate::str::contains("empty"));
}
