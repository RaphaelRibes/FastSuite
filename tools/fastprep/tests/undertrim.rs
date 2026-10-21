//! Undertrimming regressions: the reads whose adapter is hidden, and the
//! two-axis expectations that keep the fixes honest.
//!
//! These are integration tests (the binary's own output), because the
//! property under test is what the user gets: a read whose adapter is
//! present must be trimmed, and a read whose adapter is absent must not be.

use std::path::{Path, PathBuf};
use std::process::Command;

const TRUSEQ_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fastprep_undertrim_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_fastq(path: &Path, records: &[(String, String)]) {
    let mut body = String::new();
    for (id, seq) in records {
        body.push_str(&format!("@{id}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(path, body).unwrap();
}

/// Run fastprep and return (stdout, stderr, exit code).
fn run(args: &[&str]) -> (String, String, i32) {
    let out = Command::new(fastprep())
        .args(args)
        .output()
        .expect("spawn fastprep");
    (
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

/// The length of each sequence in a FASTQ file.
fn seq_lengths(path: &Path) -> Vec<usize> {
    let text = std::fs::read_to_string(path).expect("read output");
    let lines: Vec<&str> = text.lines().collect();
    let mut out = Vec::new();
    let mut i = 1;
    while i < lines.len() {
        out.push(lines[i].len());
        i += 4;
    }
    out
}

/// Deterministic pseudo-random ACGT.
fn dna(len: usize, seed: u64) -> String {
    let mut x = seed | 1;
    let mut s = String::with_capacity(len);
    for _ in 0..len {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        s.push(match x % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        });
    }
    s
}

// ---------------------------------------------------------------------------
// UT-2: a terminal no-signal G-run must not hide the adapter.
// ---------------------------------------------------------------------------

/// A 150-base read whose last 20 bases are the TruSeq adapter, with the last
/// `g` of them replaced by G (NovaSeq two-colour no-signal overprint).
fn read_with_g_run(g: usize) -> String {
    let mut seq = dna(130, 0x5EED);
    let mut tail: Vec<u8> = TRUSEQ_R1[..20].to_vec();
    for j in tail.len() - g..tail.len() {
        tail[j] = b'G';
    }
    seq.push_str(std::str::from_utf8(&tail).unwrap());
    seq
}

#[test]
fn a_terminal_g_run_does_not_hide_the_adapter() {
    // Without the retry, a 20-base adapter tail whose last 6 bases are G is
    // not found at all: the G's sit inside the 3'-anchored window and spend
    // the mismatch budget. With the retry, the read is cut exactly.
    let d = tmpdir("g_run");
    let input = d.join("in.fq");
    let records: Vec<(String, String)> = (0..20)
        .map(|k| (format!("g{k}"), read_with_g_run(6)))
        .collect();
    write_fastq(&input, &records);

    // The shipped behaviour: the adapter is NOT found, so the read keeps its
    // full length (the G's alone are removed by polyX, if at all).
    let out_off = d.join("off.fq");
    let (_, _, rc) = run(&[
        "-1", input.to_str().unwrap(),
        "-o", out_off.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
    ]);
    assert_eq!(rc, 0, "fastprep must succeed");
    let len_off = seq_lengths(&out_off)[0];
    assert_eq!(
        len_off, 150,
        "without the retry the adapter is hidden by the G-run (today's behaviour)"
    );

    // With the retry: the adapter AND the G-run are cut, exactly.
    let out_on = d.join("on.fq");
    let (_, _, rc) = run(&[
        "-1", input.to_str().unwrap(),
        "-o", out_on.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--adapter-g-run-retry", "6",
        "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
    ]);
    assert_eq!(rc, 0, "fastprep must succeed");
    let len_on = seq_lengths(&out_on)[0];
    assert_eq!(
        len_on, 130,
        "with the retry the adapter and the G-run are both removed"
    );
}

#[test]
fn a_short_g_run_does_not_trigger_the_retry() {
    // g=4 is inside the mismatch budget already: the adapter matches without
    // any retry, so the retry must not change the result.
    let d = tmpdir("g_short");
    let input = d.join("in.fq");
    let records: Vec<(String, String)> = (0..8)
        .map(|k| (format!("s{k}"), read_with_g_run(4)))
        .collect();
    write_fastq(&input, &records);

    let out_a = d.join("a.fq");
    let out_b = d.join("b.fq");
    for (out, extra) in [(&out_a, vec![]), (&out_b, vec!["--adapter-g-run-retry", "6"])] {
        let mut args = vec![
            "-1", input.to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--adapter-preset", "illumina",
            "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
        ];
        args.extend(extra);
        let (_, _, rc) = run(&args);
        assert_eq!(rc, 0, "fastprep must succeed");
    }
    assert_eq!(
        std::fs::read(&out_a).unwrap(),
        std::fs::read(&out_b).unwrap(),
        "a 4-base G-run is inside the mismatch budget: the retry changes nothing"
    );
}

#[test]
fn a_clean_read_with_a_terminal_g_run_is_not_trimmed() {
    // The overtrim axis: a read that carries NO adapter but ends in a G-run
    // must not be cut. The retry matches on the pre-G slice, and a clean
    // read has no adapter there.
    let d = tmpdir("g_clean");
    let input = d.join("in.fq");
    let records: Vec<(String, String)> = (0..20)
        .map(|k| {
            let mut seq = dna(140, 0xC0FFEE + k);
            seq.push_str("GGGGGGGG");
            (format!("c{k}"), seq)
        })
        .collect();
    write_fastq(&input, &records);

    let out_a = d.join("a.fq");
    let out_b = d.join("b.fq");
    for (out, extra) in [(&out_a, vec![]), (&out_b, vec!["--adapter-g-run-retry", "6"])] {
        let mut args = vec![
            "-1", input.to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--adapter-preset", "illumina",
            "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
        ];
        args.extend(extra);
        let (_, _, rc) = run(&args);
        assert_eq!(rc, 0, "fastprep must succeed");
    }
    assert_eq!(
        std::fs::read(&out_a).unwrap(),
        std::fs::read(&out_b).unwrap(),
        "a clean read with a terminal G-run must be unaffected by the retry"
    );
}
