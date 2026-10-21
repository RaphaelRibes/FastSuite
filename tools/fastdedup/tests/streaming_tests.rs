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

#[test]
fn the_documented_stream_refusals_report_their_own_reason() {
    // README §83 lists five refusals. Each must name its OWN cause: until the
    // validation order was fixed, `-1 - -2 -` reported "an output file is
    // required" because that generic bail ran first and shadowed the real one.
    let cases: [(&[&str], &str); 5] = [
        (&["-1", "-", "-2", "-"], "one stream cannot carry two files"),
        (&["-1", "IN", "-o", "-", "-p", "-"], "cannot both be stdout"),
        (&["-1", "IN", "--json", "-", "-o", "-"], "cannot both write to stdout"),
        (&["-1", "-", "--pre-count"], "impossible on a stream"),
        (&["-1", "-", "--abundance"], "needs a second pass over the input"),
    ];
    let dir = assert_fs::TempDir::new().unwrap();
    let input = dir.child("in.fq");
    std::fs::write(input.path(), FASTQ).unwrap();

    for (i, (args, reason)) in cases.iter().enumerate() {
        let full: Vec<String> = args
            .iter()
            .map(|a| if *a == "IN" { input.path().to_string_lossy().into_owned() } else { a.to_string() })
            .collect();
        let out = cmd().args(&full).assert().failure();
        let stderr = String::from_utf8(out.get_output().stderr.clone()).unwrap();
        assert!(stderr.contains(reason), "case {i}: expected the refusal for '{reason}', got: {stderr}");
        assert!(!stderr.contains("An output file is required"), "case {i}: the generic bail shadowed the specific refusal: {stderr}");
    }
}

// ================================================================
// Task 2.5: output ORDER is not part of the contract, the SET is
// ================================================================

/// The kept SET must be identical across thread counts; the emit ORDER need not
/// be. This asserts the half that IS contractual, and deliberately does not
/// assert the half that is not — a byte-for-byte comparison here would be
/// asserting something the tool never promised, and would fail spuriously the
/// first time worker interleaving changed.
///
/// Measured on a 200k-pair corpus at -j 1/2/4/8: one distinct sorted digest
/// across all four runs, three distinct start-to-finish orderings. The corpus
/// here is small, so it does not reproduce the order difference itself; the
/// point of the test is the set equality plus the documented edge.
#[test]
fn the_kept_set_is_stable_across_thread_counts() {
    let dir = assert_fs::TempDir::new().unwrap();
    let input = dir.child("in.fq");

    // 400 reads, every fifth a duplicate of its predecessor, so the dedup path
    // has real work and the keep set is a strict subset of the input.
    let mut body = String::new();
    for i in 0..400usize {
        let idx = if i % 5 == 4 { i - 1 } else { i };
        body.push_str(&format!("@r{idx}\nACGTACGTACGT{idx:04}\n+\nIIIIIIIIIIIIIIII\n"));
    }
    fs::write(input.path(), body.as_bytes()).unwrap();

    let mut kept_sets = Vec::new();
    for j in ["1", "2", "4", "8"] {
        let out = dir.child(format!("out{j}.fq"));
        cmd()
            .args(["-1", input.path().to_str().unwrap(),
                   "-o", out.path().to_str().unwrap(),
                   "-j", j, "--force", "-q"])
            .assert()
            .success();
        let text = fs::read_to_string(out.path()).unwrap();
        // Sorted records: the order-insensitive view.
        let mut records: Vec<String> = text
            .split('@')
            .filter(|s| !s.trim().is_empty())
            .map(|s| s.trim().replace('\n', "\u{1f}"))
            .collect();
        records.sort();
        kept_sets.push(records);
    }

    let first_len = kept_sets[0].len();
    // Non-vacuity: 400 reads with every fifth a duplicate of its predecessor
    // leaves 320. Asserting the exact number means this test fails if the
    // parser above silently produced nothing, and proves dedup actually ran.
    assert_eq!(
        first_len, 320,
        "expected 320 kept records (400 reads, 80 duplicates); got {first_len}. \
         A different number means either dedup did not run or the record \
         parser in this test is wrong"
    );
    for (i, set) in kept_sets.iter().enumerate() {
        assert_eq!(
            set.len(), first_len,
            "thread count {} changed the number of kept records: {} vs {}",
            ["1", "2", "4", "8"][i], set.len(), first_len
        );
        assert_eq!(
            set, &kept_sets[0],
            "thread count {} changed the KEPT SET, which is contractual; only \
             the order is allowed to vary",
            ["1", "2", "4", "8"][i]
        );
    }
}
