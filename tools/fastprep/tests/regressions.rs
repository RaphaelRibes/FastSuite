//! Regression tests for three default-path data-loss defects found in the
//! 2026-10-20 review. Each test fails on the pre-fix code:
//!
//! 1. clean reads at or below ~30 bp were dropped as `polyX 5′` because the
//!    budget short-circuit counted STEPS, not matches (algorithm.rs), and the
//!    3′ twin dropped the rest as `polyX 3′`;
//! 2. every record of a FASTA input was dropped: FASTA is ingested with
//!    quality 0 and the default `--mean-q 25` averages over it (pipeline.rs);
//! 3. one Q43 base (a real Illumina score, ASCII 'L' = 76) was taken as proof
//!    of Phred+64, reinterpreting the run 31 points lower and dropping the
//!    whole file on the mean-quality filter (pipeline.rs).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_regressions_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> (String, i32) {
    let out = Command::new(fastprep())
        .args(args)
        .output()
        .expect("spawn fastprep");
    (
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

fn records(path: &Path) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).expect("read output");
    let mut out = Vec::new();
    let mut lines = text.lines();
    while let Some(h) = lines.next() {
        let seq = lines.next().unwrap_or_default().to_string();
        // FASTA has no '+'/quality pair; FASTQ does.
        if h.starts_with('@') {
            lines.next();
            lines.next();
        }
        out.push((h.to_string(), seq));
    }
    out
}

/// Clean reads at 24, 27, 30, 33 and 36 bp — no homopolymer run, so no trim is
/// due. At `--min-len 0.9L` the trim budget is 2–3 bases, which is exactly the
/// width that used to trip the budget short-circuit.
#[test]
fn clean_short_reads_survive_the_default_min_len() {
    let d = tmpdir("short");
    let fq = d.join("short.fastq");
    let seqs = [
        "ACTCAGTACCTGATCAGTACCTGAT",         // 24
        "ACTCAGTACCTGATCAGTACCTGATCG",       // 26
        "ACTCAGTACCTGATCAGTACCTGATCAGTAC",   // 30
        "ACTCAGTACCTGATCAGTACCTGATCAGTACCTA", // 33
        "ACTCAGTACCTGATCAGTACCTGATCAGTACCTAGCTA", // 36
    ];
    let mut body = String::new();
    for (i, s) in seqs.iter().enumerate() {
        body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(s.len())));
    }
    std::fs::write(&fq, body).unwrap();

    let out = d.join("out.fastq");
    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "--json",
        d.join("r.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");
    let kept = records(&out);
    assert_eq!(
        kept.len(),
        seqs.len(),
        "clean reads were dropped as polyX: {err}\nkept: {kept:?}"
    );
    for (h, s) in &kept {
        assert!(
            seqs.contains(&s.as_str()),
            "kept record {h} is not one of the inputs: {s}"
        );
    }
    std::fs::remove_dir_all(&d).ok();
}

/// FASTA input has no qualities; the mean-quality filter is undefined for it
/// and must not drop the record.
#[test]
fn fasta_records_survive_the_default_mean_q() {
    let d = tmpdir("fasta");
    let fa = d.join("in.fasta");
    let seq = "ACTCAGTACCTGATCAGTACCTGATCAGTACCTAGCTAACTCAGTACCTGATCAGTACCTG\
               ATCAGTACCTAGCTAACTCAGTACCTGATCAGTACCTGATCAGTACCTAGCTAACTCA";
    std::fs::write(&fa, format!(">fa1\n{seq}\n")).unwrap();

    let out = d.join("out.fasta");
    let (err, code) = run(&[
        "-1",
        fa.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "--json",
        d.join("r.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");
    let kept = records(&out);
    assert_eq!(kept.len(), 1, "the FASTA record was dropped: {err}");
    assert_eq!(
        kept[0].1.replace([' ', '\n'], ""),
        seq.replace([' ', '\n'], ""),
        "the record must be kept unchanged (no quality evidence, no trim)"
    );
    std::fs::remove_dir_all(&d).ok();
}

/// A single Q43 base in an otherwise Phred+33 file is an outlier, not proof of
/// Phred+64: detecting +64 would subtract 31 from every quality and the default
/// mean-quality filter would then drop the whole file.
#[test]
fn a_single_q43_base_does_not_flip_to_phred64() {
    let d = tmpdir("phred");
    let fq = d.join("phred.fastq");
    let seq = "ACTCAGTACCTGATCAGTACCTGATCAGTACCTAGCTAACTCAGTACCTGATCAGTACCTGATCA";
    let mut body = String::new();
    for i in 0..50 {
        let mut q = vec![b'I'; seq.len()];
        if i == 0 {
            q[10] = b'L'; // Q43 under Phred+33
        }
        body.push_str(&format!(
            "@r{i}\n{seq}\n+\n{}\n",
            String::from_utf8(q).unwrap()
        ));
    }
    std::fs::write(&fq, body).unwrap();

    let out = d.join("out.fastq");
    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        // The fixture repeats one sequence 50 times; dedup is on by default and
        // would keep one — this test is about the Phred offset, not duplicates.
        "--no-dedup",
        "--json",
        d.join("r.json").to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");
    let kept = records(&out);
    assert_eq!(
        kept.len(),
        50,
        "a single Q43 base flipped the run to Phred+64 and the mean-q filter \
         dropped everything: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}
