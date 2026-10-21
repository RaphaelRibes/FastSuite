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

// ---------------------------------------------------------------------------
// UT-3: a short adapter tail is trimmed only with corroboration.
// ---------------------------------------------------------------------------

/// A FASTQ quality string of a constant Phred value.
fn qual_str(q: u8, n: usize) -> String {
    std::iter::repeat((q + 33) as char).take(n).collect()
}

#[test]
fn a_short_exact_prefix_is_trimmed_when_the_quality_cliff_corroborates() {
    // A 10-base adapter tail with FLAT quality is not trimmed (no evidence).
    // The same tail with a sharp quality cliff at the boundary IS trimmed.
    let d = tmpdir("st_corr");
    let input = d.join("in.fq");
    let body = dna(140, 0x51_51);
    let tail = std::str::from_utf8(&TRUSEQ_R1[..10]).unwrap().to_string();

    // flat quality: no cliff -> not trimmed
    let flat = format!("{body}{tail}");
    // cliff: 4 low-quality bases just before the boundary
    let mut cliff_seq = body.clone();
    cliff_seq.push_str(&tail);
    let mut cliff_q = qual_str(40, 140);
    cliff_q.push_str(&qual_str(40, 10)); // placeholder, replaced below

    write_fastq(&input, &[(String::from("flat"), flat)]);
    let out_flat = d.join("flat.fq");
    let mut args = vec![
        "-1", input.to_str().unwrap(),
        "-o", out_flat.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--adapter-short-tail", "6",
        "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
    ];
    let (_, _, rc) = run(&args);
    assert_eq!(rc, 0);
    assert_eq!(
        seq_lengths(&out_flat)[0],
        150,
        "a flat-quality short tail is not evidence"
    );

    // Now the corroborated case, with an explicit quality cliff.
    let d2 = tmpdir("st_corr2");
    let input2 = d2.join("in.fq");
    // qualities: 40 for the body, then 10 for the 10 bases before the tail,
    // then 40 for the tail itself
    let mut seq2 = body.clone();
    seq2.push_str(&tail);
    let mut q2 = qual_str(40, 140 - 10);
    q2.push_str(&qual_str(10, 10));
    q2.push_str(&qual_str(40, 10));
    let mut body2 = String::new();
    body2.push_str(&format!("@cliff\n{seq2}\n+\n{q2}\n"));
    std::fs::write(&input2, body2).unwrap();
    let out_cliff = d2.join("cliff.fq");
    args[1] = input2.to_str().unwrap();
    args[3] = out_cliff.to_str().unwrap();
    let (_, _, rc) = run(&args);
    assert_eq!(rc, 0);
    assert_eq!(
        seq_lengths(&out_cliff)[0],
        140,
        "a short exact prefix with a quality cliff at the boundary is trimmed"
    );
}

#[test]
fn a_short_prefix_without_corroboration_is_not_trimmed() {
    // The overtrim guard: a clean read whose last 8 bases happen to equal an
    // adapter prefix, with FLAT quality, must not be cut.
    let d = tmpdir("st_clean");
    let input = d.join("in.fq");
    let body = dna(142, 0x7E57);
    let tail = std::str::from_utf8(&TRUSEQ_R1[..8]).unwrap().to_string();
    write_fastq(&input, &[(String::from("clean"), format!("{body}{tail}"))]);
    let out = d.join("o.fq");
    let (_, _, rc) = run(&[
        "-1", input.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--adapter-short-tail", "6",
        "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
    ]);
    assert_eq!(rc, 0);
    assert_eq!(
        seq_lengths(&out)[0],
        150,
        "an uncorroborated short prefix must not be trimmed"
    );
}

#[test]
fn the_short_tail_path_is_off_by_default() {
    let d = tmpdir("st_default");
    let input = d.join("in.fq");
    let body = dna(140, 0x51_51);
    let tail = std::str::from_utf8(&TRUSEQ_R1[..10]).unwrap().to_string();
    let seq = format!("{body}{tail}");
    let mut q = qual_str(40, 140 - 10);
    q.push_str(&qual_str(10, 10));
    q.push_str(&qual_str(40, 10));
    let mut body2 = String::new();
    body2.push_str(&format!("@x\n{seq}\n+\n{q}\n"));
    std::fs::write(&input, body2).unwrap();

    let out_a = d.join("a.fq");
    let out_b = d.join("b.fq");
    for (out, extra) in [(&out_a, vec![]), (&out_b, vec!["--adapter-short-tail", "6"])] {
        let mut args = vec![
            "-1", input.to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--adapter-preset", "illumina",
            "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
        ];
        args.extend(extra);
        let (_, _, rc) = run(&args);
        assert_eq!(rc, 0);
    }
    assert_eq!(
        seq_lengths(&out_a)[0],
        150,
        "without the flag the short tail is not trimmed"
    );
    assert_eq!(
        seq_lengths(&out_b)[0],
        140,
        "with the flag the corroborated short tail is trimmed"
    );
}

// ---------------------------------------------------------------------------
// UT-3: the dataset-level enrichment signal.
// ---------------------------------------------------------------------------

#[test]
fn enrichment_corroborates_a_short_prefix_that_chance_would_not_explain() {
    // 200 reads, of which 100 end with the exact 8-base TruSeq prefix.
    // Chance predicts 200 * 4^-8 = 0.0078 reads; 100 observed is an
    // enrichment of ~12800x. A short read ending in that prefix is trimmed.
    let d = tmpdir("enr");
    let input = d.join("in.fq");
    let tail = std::str::from_utf8(&TRUSEQ_R1[..8]).unwrap().to_string();
    let mut records: Vec<(String, String)> = Vec::new();
    for k in 0..100 {
        records.push((format!("hit{k}"), format!("{}{}", dna(142, k), tail)));
    }
    for k in 0..100 {
        records.push((format!("miss{k}"), dna(150, 1000 + k)));
    }
    write_fastq(&input, &records);

    let out = d.join("o.fq");
    let (_, _, rc) = run(&[
        "-1", input.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--adapter-short-tail", "6",
        "--adapter-min-enrichment", "10",
        "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
    ]);
    assert_eq!(rc, 0);
    let lens = seq_lengths(&out);
    let hits = &lens[..100];
    let misses = &lens[100..];
    let trimmed = hits.iter().filter(|&&l| l < 150).count();
    let missed_trimmed = misses.iter().filter(|&&l| l < 150).count();
    assert!(
        trimmed > 50,
        "an enriched short prefix must be trimmed: {trimmed}/100"
    );
    assert_eq!(missed_trimmed, 0, "unrelated reads must be untouched");
}

#[test]
fn enrichment_is_disabled_when_the_factor_is_zero() {
    let d = tmpdir("enr0");
    let input = d.join("in.fq");
    let tail = std::str::from_utf8(&TRUSEQ_R1[..8]).unwrap().to_string();
    let records: Vec<(String, String)> = (0..100)
        .map(|k| (format!("hit{k}"), format!("{}{}", dna(142, k), tail)))
        .collect();
    write_fastq(&input, &records);
    let out = d.join("o.fq");
    let (_, _, rc) = run(&[
        "-1", input.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--adapter-short-tail", "6",
        // factor 0 must DISABLE the signal, not make every prefix pass.
        "--no-dedup", "--no-mean-q", "--min-len", "0", "-q",
    ]);
    assert_eq!(rc, 0);
    assert!(
        seq_lengths(&out).iter().all(|&l| l == 150),
        "with the enrichment signal off, nothing is trimmed"
    );
}
