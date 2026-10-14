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

fn gz(body: &[u8]) -> Vec<u8> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(body).unwrap();
    enc.finish().unwrap()
}

// ================================================================
// Task 2: stream refusals
// ================================================================

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

// ================================================================
// Task 3: stdin reading in the single-end serial path
// ================================================================

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

// ================================================================
// Task 4: stream output — flags, writer, summary-to-stderr rule
// ================================================================

#[test]
fn stdin_to_stdout_plain() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "2 reads in, 1 unique out"
    );
    assert!(out.starts_with('@'), "stdout must be FASTQ only, no summary");
    let err = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(
        err.contains("fastdedup: 2 reads in, 1 out"),
        "summary belongs on stderr when stdout carries data: {err}"
    );
}

// ================================================================
// Task 5: gzip on stdin is decoded by magic bytes
// ================================================================

#[test]
fn gzipped_stdin_is_decoded_from_magic_bytes() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .write_stdin(gz(FASTQ))
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "the gzip stream must be transparently decoded"
    );
}

// ================================================================
// Task 6: --out-compression gz writes gzip to stdout
// ================================================================

#[test]
fn stdout_can_be_gzip() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("--out-compression")
        .arg("gz")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let bytes = assert.get_output().stdout.clone();
    assert_eq!(bytes[0], 0x1f, "gzip magic byte");
    assert_eq!(bytes[1], 0x8b, "gzip magic byte");
}

// ================================================================
// Task 7: --out-format fasta on stdout, and FASTA-in → FASTQ-out refusal
// ================================================================

#[test]
fn stdout_can_be_fasta() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("--out-format")
        .arg("fasta")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(out.starts_with('>'), "FASTA output must start with '>'");
    assert!(!out.contains('+'), "FASTA has no '+' separator line");
}

#[test]
fn fasta_stdin_to_fastq_stdout_is_refused() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .write_stdin(b">r1\nACGTACGT\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "FASTA → FASTQ conversion not supported",
        ));
}

// ================================================================
// Task 8: --json - writes the report to stdout, summary to stderr
// ================================================================

#[test]
fn json_report_to_stdout_moves_summary_to_stderr() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .arg("--json")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(
        out.contains("fastsuite.report.v1"),
        "stdout must carry the report: {out}"
    );
    let err = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(
        err.contains("fastdedup: 2 reads in"),
        "summary must move to stderr: {err}"
    );
}

#[test]
fn json_and_fastq_on_stdout_is_refused() {
    let _temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("--json")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("--json (-)"));
}

// ================================================================
// Task 9: paired-end serial path streams (stdin R1, stream outputs)
// ================================================================

#[test]
fn paired_stdin_r1_with_file_r2() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let r2 = temp_dir.child("R2.fastq");
    r2.write_str("@p1 2:N:0:ATGC\nTTTTTTTT\n+\nJJJJJJJJ\n").unwrap();
    let out2 = temp_dir.child("out2.fastq");

    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg("-")
        .arg("-p")
        .arg(out2.path())
        .write_stdin(b"@p1 1:N:0:ATGC\nACGTACGT\n+\nIIIIIIII\n")
        .assert()
        .success();

    let out1 = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out1.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "one unique pair survives"
    );
    let out2_body = fs::read_to_string(out2.path()).unwrap();
    assert_eq!(
        out2_body.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "R2 file must receive the same single pair"
    );
}

// ================================================================
// Task 10: abundance mode streams (with --temp-dir)
// ================================================================

#[test]
fn abundance_from_stdin_with_temp_dir_works() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let spill = temp_dir.child("spill");
    spill.create_dir_all().unwrap();

    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("--abundance")
        .arg("--temp-dir")
        .arg(spill.path())
        .arg("-o")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .success();

    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(
        out.contains(";size=2"),
        "the duplicated read must be annotated with its abundance: {out}"
    );
}

// ================================================================
// Task 11: multi-threaded pipeline streams (stdin in, stdout out)
// ================================================================

#[test]
fn threaded_stdin_to_stdout() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("-j")
        .arg("4")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "the threaded pipeline must read stdin and write stdout identically"
    );
}
