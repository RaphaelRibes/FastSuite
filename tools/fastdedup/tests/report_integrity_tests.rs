//! Report-integrity and paired-end-synchronization regressions.
//!
//! Every test in this file fails on the binary as of 2026-10-26:
//!
//!   D1  `--json` naming an INPUT path replaced that input with the report.
//!   D2  the serial paired-end engines dropped a one-record surplus mate and
//!       exited 0 (the parallel engine refused).
//!   D3  a resume reported only the records written by the resuming run, so a
//!       complete 20 000-record output with one repaired tail printed
//!       "1 out, 19999 duplicates (100.00%)".
//!   D4  `metrics.abundance` was a JSON boolean (REPORT-SCHEMA rule 1: numbers).
//!   D5  `--abundance` reported the requested `-j`, not the single thread it
//!       runs, and `FASTSUITE_THREADS` was ignored entirely.

use assert_cmd::Command;
use assert_cmd::cargo;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;

fn cmd() -> Command {
    cargo::cargo_bin_cmd!("fastdedup")
}

// ================================================================
// FIXTURES
// ================================================================

/// Deterministic 60 bp ACGT sequence for record `i` of allele `mate`.
/// Different `mate` values never share a sequence, so R1 and R2 pairs are
/// distinct while carrying the same base read ID.
fn sequence(i: usize, mate: usize) -> String {
    const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut state = (i as u64 + 1)
        .wrapping_mul(2_654_435_761)
        .wrapping_add((mate as u64).wrapping_mul(40_503));
    let mut seq = Vec::with_capacity(60);
    for _ in 0..60 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        seq.push(NUCLEOTIDES[((state >> 33) & 3) as usize]);
    }
    String::from_utf8(seq).unwrap()
}

/// `n` FASTQ records; `mate` 0 = single-end, 1 = R1, 2 = R2.
fn fastq(n: usize, mate: usize) -> String {
    let mut body = String::new();
    for i in 0..n {
        let suffix = if mate == 0 {
            String::new()
        } else {
            format!("/{mate}")
        };
        body.push_str(&format!(
            "@read{i}{suffix}\n{}\n+\n{}\n",
            sequence(i, mate),
            "I".repeat(60)
        ));
    }
    body
}

fn count_records(path: &std::path::Path) -> usize {
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter(|l| l.starts_with('@'))
        .count()
}

fn read_report(path: &std::path::Path) -> serde_json::Value {
    let text = fs::read_to_string(path).unwrap();
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("unreadable report ({e}): {text}"))
}

/// "N threads," out of the human summary line.
fn summary_threads(stdout: &[u8]) -> u64 {
    let text = String::from_utf8_lossy(stdout);
    let idx = text
        .find(" threads,")
        .unwrap_or_else(|| panic!("no thread count in summary: {text}"));
    text[..idx]
        .rsplit(", ")
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap()
}

// ================================================================
// D1 — `--json` must not overwrite an input
// ================================================================

/// `--json victim.fq -1 victim.fq` used to replace the reads with the report
/// and exit 0: silent loss of the only copy of the data.
#[test]
fn json_naming_the_input_is_refused() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("victim.fq").write_str(&fastq(5, 0)).unwrap();
    let before = fs::read(temp.child("victim.fq").path()).unwrap();

    for force in [false, true] {
        let mut c = cmd();
        c.current_dir(temp.path())
            .arg("-1")
            .arg("victim.fq")
            .arg("-o")
            .arg("out.fq")
            .arg("--json")
            .arg("victim.fq");
        if force {
            c.arg("--force");
        }
        c.assert()
            .failure()
            .stderr(predicate::str::contains("--json"))
            .stderr(predicate::str::contains("victim.fq"))
            .stderr(predicate::str::contains("same file"));
    }

    assert_eq!(
        fs::read(temp.child("victim.fq").path()).unwrap(),
        before,
        "the input must not be replaced by the report"
    );
}

/// Same file, different spelling: the guard resolves both paths through the
/// filesystem (`path_identity`), so `./victim.fq` is not a distinct target.
#[test]
fn json_naming_the_input_through_another_spelling_is_refused() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("victim.fq").write_str(&fastq(4, 0)).unwrap();
    let before = fs::read(temp.child("victim.fq").path()).unwrap();

    cmd()
        .current_dir(temp.path())
        .arg("-1")
        .arg("victim.fq")
        .arg("-o")
        .arg("out.fq")
        .arg("--json")
        .arg("./victim.fq")
        .arg("--force")
        .assert()
        .failure()
        .stderr(predicate::str::contains("--json"));

    assert_eq!(fs::read(temp.child("victim.fq").path()).unwrap(), before);
}

/// The R2 input is an input too.
#[test]
fn json_naming_the_r2_input_is_refused() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("R1.fq").write_str(&fastq(4, 1)).unwrap();
    temp.child("R2.fq").write_str(&fastq(4, 2)).unwrap();
    let before = fs::read(temp.child("R2.fq").path()).unwrap();

    cmd()
        .arg("-1")
        .arg(temp.child("R1.fq").path())
        .arg("-2")
        .arg(temp.child("R2.fq").path())
        .arg("-o")
        .arg(temp.child("o1.fq").path())
        .arg("-p")
        .arg(temp.child("o2.fq").path())
        .arg("--json")
        .arg(temp.child("R2.fq").path())
        .arg("--force")
        .assert()
        .failure()
        .stderr(
            predicate::str::contains("--json").and(predicate::str::contains("--input-r2 (-2)")),
        );

    assert_eq!(fs::read(temp.child("R2.fq").path()).unwrap(), before);
}

/// Control: a report on its own path stays accepted and leaves the input intact.
#[test]
fn json_on_its_own_path_is_accepted() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("in.fq").write_str(&fastq(5, 0)).unwrap();
    let before = fs::read(temp.child("in.fq").path()).unwrap();
    let report = temp.child("report.json");

    cmd()
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(temp.child("out.fq").path())
        .arg("--json")
        .arg(report.path())
        .arg("--force")
        .arg("-q")
        .assert()
        .success();

    assert_eq!(read_report(report.path())["tool"], "fastdedup");
    assert_eq!(fs::read(temp.child("in.fq").path()).unwrap(), before);
    assert_eq!(count_records(temp.child("out.fq").path()), 5);
}

// ================================================================
// D2 — a surplus mate record must be refused, never dropped
// ================================================================

/// Runs the desynchronized pair and returns the two output paths.
fn run_desync_pair(
    temp: &assert_fs::TempDir,
    r1_records: usize,
    r2_records: usize,
    extra: &[&str],
) -> assert_cmd::assert::Assert {
    temp.child("R1.fq")
        .write_str(&fastq(r1_records, 1))
        .unwrap();
    temp.child("R2.fq")
        .write_str(&fastq(r2_records, 2))
        .unwrap();

    let mut c = cmd();
    c.arg("-1")
        .arg(temp.child("R1.fq").path())
        .arg("-2")
        .arg(temp.child("R2.fq").path());
    for arg in extra {
        c.arg(arg);
    }
    c.arg("-o")
        .arg(temp.child("o1.fq").path())
        .arg("-p")
        .arg(temp.child("o2.fq").path())
        .assert()
}

/// Three R1 records, two R2 records: the third R1 record used to be consumed
/// and thrown away by the failing loop iteration, the post-loop check then saw
/// EOF, and the run exited 0 with "2 reads in, 2 out" for a 5-record input.
#[test]
fn pe_surplus_on_r1_is_refused_at_every_thread_count() {
    for threads in ["1", "2", "4"] {
        let temp = assert_fs::TempDir::new().unwrap();
        let assert = run_desync_pair(&temp, 3, 2, &["-j", threads, "--force"]);
        assert.failure().stderr(
            predicate::str::contains("R1 has more records than R2")
                .and(predicate::str::contains("5 records consumed")),
        );
    }
}

/// The mirror case, R2 one record longer.
#[test]
fn pe_surplus_on_r2_is_refused_at_every_thread_count() {
    for threads in ["1", "2", "4"] {
        let temp = assert_fs::TempDir::new().unwrap();
        let assert = run_desync_pair(&temp, 2, 3, &["-j", threads, "--force"]);
        assert.failure().stderr(
            predicate::str::contains("R2 has more records than R1")
                .and(predicate::str::contains("5 records consumed")),
        );
    }
}

/// The default (`-j` absent, 1 thread) is the path the bug lived in.
#[test]
fn pe_surplus_is_refused_on_the_default_thread_count() {
    let temp = assert_fs::TempDir::new().unwrap();
    run_desync_pair(&temp, 3, 2, &["--force"])
        .failure()
        .stderr(predicate::str::contains("R1 has more records than R2"));
}

/// `--dry-run` reported "2 reads in" for the same 3+2 input: a dry run may
/// write nothing, but it must not hide the unpaired record either.
#[test]
fn pe_surplus_is_refused_in_dry_run() {
    for (r1, r2, needle) in [
        (3usize, 2usize, "R1 has more records than R2"),
        (2, 3, "R2 has more records than R1"),
    ] {
        let temp = assert_fs::TempDir::new().unwrap();
        run_desync_pair(&temp, r1, r2, &["--dry-run"])
            .failure()
            .stderr(predicate::str::contains(needle));
    }
}

/// Both engines must reach the same verdict, with the same wording: the serial
/// path said "detected at the end" only after a second silent record, the
/// parallel path refused at the first.
#[test]
fn pe_serial_and_parallel_refuse_identically() {
    let temp = assert_fs::TempDir::new().unwrap();
    let serial = run_desync_pair(&temp, 3, 2, &["-j", "1", "--force"]);
    let serial_err = String::from_utf8_lossy(&serial.get_output().stderr).to_string();
    serial.failure();

    let temp2 = assert_fs::TempDir::new().unwrap();
    let parallel = run_desync_pair(&temp2, 3, 2, &["-j", "4", "--force"]);
    let parallel_err = String::from_utf8_lossy(&parallel.get_output().stderr).to_string();
    parallel.failure();

    let needle = "Desynchronization: R1 has more records than R2";
    assert!(
        serial_err.contains(needle),
        "serial: expected {needle:?}, got {serial_err:?}"
    );
    assert!(
        parallel_err.contains(needle),
        "parallel: expected {needle:?}, got {parallel_err:?}"
    );
    assert_eq!(
        serial_err, parallel_err,
        "the two engines must refuse with the same message"
    );
}

/// Control: equal lengths still succeed, so the refusal is not a blanket one.
#[test]
fn pe_equal_length_mates_still_succeed() {
    let temp = assert_fs::TempDir::new().unwrap();
    run_desync_pair(&temp, 3, 3, &["--force"]).success();
    assert_eq!(count_records(temp.child("o1.fq").path()), 3);
    assert_eq!(count_records(temp.child("o2.fq").path()), 3);
}

// ================================================================
// D3 — the resume report describes the whole output and the whole input
// ================================================================

/// A completed 200-record run, its output truncated mid-record, then resumed:
/// the resuming run repairs one record but the output holds all 200 and the
/// input holds no duplicate. The report used to announce "1 out, 199
/// duplicates (99.50%)".
#[test]
fn resume_report_describes_the_whole_output() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("in.fq").write_str(&fastq(200, 0)).unwrap();
    let out = temp.child("out.fq");

    cmd()
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(out.path())
        .arg("--force")
        .arg("-q")
        .assert()
        .success();
    assert_eq!(count_records(out.path()), 200);

    // Torn tail: the last record loses 57 of its bytes, as a crash would.
    let len = fs::metadata(out.path()).unwrap().len();
    fs::OpenOptions::new()
        .write(true)
        .open(out.path())
        .unwrap()
        .set_len(len - 57)
        .unwrap();

    let report = temp.child("resume.json");
    let assert = cmd()
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(out.path())
        .arg("-j")
        .arg("1")
        .arg("--json")
        .arg(report.path())
        .assert()
        .success();

    let stdout = String::from_utf8_lossy(&assert.get_output().stdout).to_string();
    assert!(
        stdout.contains("200 reads in, 200 out, 0 duplicates"),
        "summary must describe the repaired output, got: {stdout}"
    );

    let metrics = read_report(report.path())["metrics"].clone();
    assert_eq!(metrics["reads_in"].as_u64(), Some(200));
    assert_eq!(
        metrics["reads_out"].as_u64(),
        Some(200),
        "reads_out counts the whole output, including preloaded records"
    );
    assert_eq!(
        metrics["reads_duplicate"].as_u64(),
        Some(0),
        "a repaired tail is not a duplicate"
    );
    assert!(
        metrics["duplicate_rate"].as_f64().unwrap().abs() < 1e-12,
        "duplicate_rate must describe the whole input, got {}",
        metrics["duplicate_rate"]
    );
    assert_eq!(count_records(out.path()), 200);
}

// ================================================================
// D4 — REPORT-SCHEMA rule 1: metrics are numbers, never booleans
// ================================================================

#[test]
fn abundance_metric_is_a_number_not_a_boolean() {
    for (abundance, expected) in [(false, 0u64), (true, 1u64)] {
        let temp = assert_fs::TempDir::new().unwrap();
        temp.child("in.fq").write_str(&fastq(5, 0)).unwrap();
        let report = temp.child("m.json");

        let mut c = cmd();
        c.arg("-1")
            .arg(temp.child("in.fq").path())
            .arg("-o")
            .arg(temp.child("out.fq").path())
            .arg("--json")
            .arg(report.path())
            .arg("--force")
            .arg("-q");
        if abundance {
            c.arg("--abundance");
        }
        c.assert().success();

        let value = read_report(report.path())["metrics"]["abundance"].clone();
        assert!(
            value.is_number(),
            "--abundance={abundance}: metrics.abundance must be 0/1, got {value}"
        );
        assert_eq!(value.as_u64(), Some(expected));
    }
}

// ================================================================
// D5 — report the pool that runs; FASTSUITE_THREADS is a default only
// ================================================================

/// `--abundance` warns that it runs single-threaded, then reported `-j 4`.
#[test]
fn abundance_reports_the_single_thread_it_runs() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("in.fq").write_str(&fastq(5, 0)).unwrap();
    let report = temp.child("ab.json");

    let assert = cmd()
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(temp.child("out.fq").path())
        .arg("--abundance")
        .arg("-j")
        .arg("4")
        .arg("--json")
        .arg(report.path())
        .arg("--force")
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "--abundance does not support multi-threading",
        ));

    assert_eq!(
        summary_threads(&assert.get_output().stdout),
        1,
        "the summary must announce the thread that runs"
    );
    assert_eq!(
        read_report(report.path())["metrics"]["threads"].as_u64(),
        Some(1)
    );
}

#[test]
fn fastsuite_threads_is_the_default_when_the_flag_is_absent() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("in.fq").write_str(&fastq(5, 0)).unwrap();
    let report = temp.child("env.json");

    cmd()
        .env("FASTSUITE_THREADS", "3")
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(temp.child("out.fq").path())
        .arg("--json")
        .arg(report.path())
        .arg("--force")
        .arg("-q")
        .assert()
        .success();

    assert_eq!(
        read_report(report.path())["metrics"]["threads"].as_u64(),
        Some(3),
        "FASTSUITE_THREADS is the suite-wide default when -j is absent"
    );
}

/// An explicit `-j` always wins over the environment.
#[test]
fn explicit_threads_override_fastsuite_threads() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("in.fq").write_str(&fastq(5, 0)).unwrap();
    let report = temp.child("flag.json");

    cmd()
        .env("FASTSUITE_THREADS", "3")
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(temp.child("out.fq").path())
        .arg("-j")
        .arg("2")
        .arg("--json")
        .arg(report.path())
        .arg("--force")
        .arg("-q")
        .assert()
        .success();

    assert_eq!(
        read_report(report.path())["metrics"]["threads"].as_u64(),
        Some(2)
    );
}

/// An unset or unparseable default is not fatal: 1 thread, as documented.
#[test]
fn unparseable_fastsuite_threads_falls_back_to_one() {
    let temp = assert_fs::TempDir::new().unwrap();
    temp.child("in.fq").write_str(&fastq(5, 0)).unwrap();
    let report = temp.child("bogus.json");

    cmd()
        .env("FASTSUITE_THREADS", "bogus")
        .arg("-1")
        .arg(temp.child("in.fq").path())
        .arg("-o")
        .arg(temp.child("out.fq").path())
        .arg("--json")
        .arg(report.path())
        .arg("--force")
        .arg("-q")
        .assert()
        .success();

    assert_eq!(
        read_report(report.path())["metrics"]["threads"].as_u64(),
        Some(1)
    );
}
