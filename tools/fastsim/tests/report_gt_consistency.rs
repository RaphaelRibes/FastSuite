//! GT-vs-report consistency: every count the run publishes must describe
//! exactly what the GT and the FASTQ contain.
//!
//! Pins two reproduced defects. `adapter_pad_bases` counted every `N` of the
//! CONSTRUCTED read, so a template `N` (a 30-N run in `--reference`) was
//! subtracted from the GT's `r1_n`/`r2_n` columns and added to the pad metric:
//! a reference full of N produced all-zero GT N columns, a pad count of 28k and
//! "ambiguous 0 reads / 0 N bases" on stdout. And `full_artefact_pairs` counted
//! pairs where EITHER mate is full, while the GT's `full_flag` is the AND of the
//! two per-mate flags — a 10% full R1 head published 10% full pairs against a GT
//! of zero. Both counts are now published, each named for what it measures.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastsim")
}

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fastsim_gt_{}_{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A minimal profile: the two head histograms, an adapter mode, and an optional
/// insert CDF (absent -> the Gaussian fallback).
fn profile(
    dir: &Path,
    name: &str,
    read_len: u32,
    heads: (serde_json::Value, serde_json::Value),
    adapter: &str,
    insert: Option<serde_json::Value>,
) -> PathBuf {
    let q: Vec<f64> = vec![30.0; read_len as usize];
    let mut v = serde_json::json!({
        "read_len": read_len, "gc": 0.5, "geometry": "tail3p_both",
        "head": {"R1": heads.0, "R2": heads.1},
        "tail": {"R1": {"0": 100}, "R2": {"0": 100}},
        "joint_tail": {"0,0": 100},
        "tail_q_collapse": {"R1": 0.0, "R2": 0.0},
        "pos_q": {"R1": q.clone(), "R2": q},
        "read_meanq": {"R1": {"30": 100}, "R2": {"30": 100}},
        "full_length_rate": 0.0,
        "full_length_q": {"R1": 20, "R2": 20},
        "adapter": {"seq": "CTGTCTCTTATACACATCT", "mode": adapter},
    });
    if let Some(i) = insert {
        v["insert"] = i;
    }
    let p = dir.join(format!("{name}.json"));
    std::fs::write(&p, serde_json::to_string_pretty(&v).unwrap()).unwrap();
    p
}

struct Run {
    out: PathBuf,
    stderr: String,
    report: serde_json::Value,
}

fn run(dir: &Path, tag: &str, profile: &Path, pairs: u64, extra: &[&str]) -> Run {
    let out = dir.join(tag);
    let report = dir.join(format!("{tag}.json"));
    let o = Command::new(bin())
        .args([
            "--pairs", &pairs.to_string(), "--seed", "7",
            "--profile", profile.to_str().unwrap(),
            "--out", out.to_str().unwrap(),
            "--json", report.to_str().unwrap(),
            "--no-gzip", "--force", "--drop-rate", "0",
        ])
        .args(extra)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "fastsim {tag} failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    Run {
        out,
        stderr: String::from_utf8_lossy(&o.stderr).into_owned(),
        report: serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap(),
    }
}

fn seqs(out: &Path, mate: &str) -> Vec<String> {
    std::fs::read_to_string(format!("{}_{mate}.fastq", out.display()))
        .unwrap()
        .lines()
        .skip(1)
        .step_by(4)
        .map(str::to_string)
        .collect()
}

fn n_bases(seqs: &[String]) -> u64 {
    seqs.iter().map(|s| s.bytes().filter(|&b| b == b'N').count() as u64).sum()
}

fn gt_rows(out: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(format!("{}.gt.tsv", out.display()))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

fn sum_col(rows: &[Vec<String>], i: usize) -> u64 {
    rows.iter().map(|r| r[i].parse::<u64>().unwrap()).sum()
}

/// A reference whose contigs carry a 30-`N` run every 200 bases, so an insert
/// sampled from it has a realistic chance of containing template `N`.
fn n_reference(dir: &Path) -> PathBuf {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let base = |state: &mut u64| {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        b"ACGT"[(*state % 4) as usize] as char
    };
    let mut seq = String::new();
    for _ in 0..25 {
        for _ in 0..200 {
            seq.push(base(&mut state));
        }
        seq.push_str(&"N".repeat(30));
    }
    let p = dir.join("refN.fa");
    std::fs::write(&p, format!(">c1\n{seq}\n")).unwrap();
    p
}

/// S1: the GT's `r1_n`/`r2_n` columns and the report's N counters must account
/// for exactly the `N` bases the FASTQ carries, split by source: a template `N`
/// (from `--reference`) is a no-call in the read and belongs in the GT columns,
/// while the no-adapter length filler is padding and belongs in
/// `adapter_pad_bases` — never in the GT.
#[test]
fn gt_n_columns_and_pad_metric_split_the_two_n_sources() {
    let dir = scratch("n");
    let fa = n_reference(&dir);
    let fa = fa.to_str().unwrap();
    let none = (
        serde_json::json!({"0": 100}),
        serde_json::json!({"0": 100}),
    );

    // (A) template N only: a 100 bp insert covers the whole 40 bp window, so the
    // no-adapter filler never runs.
    let p = profile(&dir, "template_only", 40, none.clone(), "none", None);
    let a = run(
        &dir,
        "a",
        &p,
        500,
        &["--reference", fa, "--insert-mean", "100", "--insert-sd", "10"],
    );
    let (a1, a2) = (n_bases(&seqs(&a.out, "R1")), n_bases(&seqs(&a.out, "R2")));
    let rows = gt_rows(&a.out);
    let (g1, g2) = (sum_col(&rows, 10), sum_col(&rows, 11));
    assert!(a1 > 0 && a2 > 0, "the reference must put N in the reads: {a1} {a2}");
    assert_eq!(g1, a1, "template N in R1 must be in the GT's r1_n column");
    assert_eq!(g2, a2, "template N in R2 must be in the GT's r2_n column");
    assert_eq!(
        a.report["metrics"]["adapter_pad_bases"], 0,
        "no filler ran, so nothing is padding: {}",
        a.report["metrics"]
    );
    assert_eq!(a.report["metrics"]["ambiguous_bases"], g1 + g2);
    assert!(
        a.stderr.contains(&format!("/ {} N bases", g1 + g2)),
        "the summary must report the template N, not 0: {}",
        a.stderr
    );

    // (B) pad only: a 20 bp insert in a 40 bp window pads 20 bases per mate, and
    // with no reference every emitted N is filler.
    let pad = profile(
        &dir,
        "pad_only",
        40,
        none,
        "none",
        Some(serde_json::json!({"20": 100})),
    );
    let b = run(&dir, "b", &pad, 500, &[]);
    let (b1, b2) = (n_bases(&seqs(&b.out, "R1")), n_bases(&seqs(&b.out, "R2")));
    assert_eq!((b1, b2), (500 * 20, 500 * 20), "pad is 20 bases per mate");
    let rows = gt_rows(&b.out);
    assert_eq!(sum_col(&rows, 10) + sum_col(&rows, 11), 0, "padding is not a no-call");
    assert_eq!(b.report["metrics"]["adapter_pad_bases"], b1 + b2);
    assert_eq!(b.report["metrics"]["ambiguous_bases"], 0);
    assert!(
        b.stderr.contains("/ 0 N bases"),
        "pure padding is not an ambiguous base: {}",
        b.stderr
    );

    // (C) both at once: emitted N = template + pad, and GT plus the pad metric
    // must add back up to them exactly.
    let c = run(&dir, "c", &pad, 500, &["--reference", fa]);
    let (c1, c2) = (n_bases(&seqs(&c.out, "R1")), n_bases(&seqs(&c.out, "R2")));
    let rows = gt_rows(&c.out);
    let (g1, g2) = (sum_col(&rows, 10), sum_col(&rows, 11));
    assert!(g1 > 0 && g2 > 0, "template N survived: {g1} {g2}");
    assert_eq!(c.report["metrics"]["adapter_pad_bases"], 500 * 20 * 2);
    assert_eq!(c1, g1 + 500 * 20, "emitted R1 N = template + pad");
    assert_eq!(c2, g2 + 500 * 20, "emitted R2 N = template + pad");
    assert_eq!(c.report["metrics"]["ambiguous_bases"], g1 + g2);

    std::fs::remove_dir_all(&dir).unwrap();
}

/// S3: `full_artefact_pairs` must agree with the GT's `full_flag` (BOTH mates
/// full), while the old either-mate count stays available, named for what it is.
#[test]
fn full_pair_counters_agree_with_the_gt() {
    let dir = scratch("full");
    // 10% of R1 reads are a 60 bp head (the whole read) -> R1 full, R2 never.
    let p = profile(
        &dir,
        "head_r1",
        60,
        (
            serde_json::json!({"0": 90, "60": 10}),
            serde_json::json!({"0": 100}),
        ),
        "none",
        Some(serde_json::json!({"120": 100})),
    );

    // A 20% duplicate rate exercises the replay path too: a duplicate of a full
    // pair carries the stored artefact, so its GT row is full as well.
    let r = run(&dir, "paired", &p, 2000, &["--dup-rate", "0.2"]);
    let rows = gt_rows(&r.out);
    let count = |f: &dyn Fn(&Vec<String>) -> bool| rows.iter().filter(|x| f(x)).count();
    let (f1, f2, flagged) = (
        count(&|x| x[12] == "1"),
        count(&|x| x[13] == "1"),
        count(&|x| x[8] == "1"),
    );
    let either = count(&|x| x[12] == "1" || x[13] == "1");
    let both = count(&|x| x[12] == "1" && x[13] == "1");

    assert_eq!(flagged, both, "the GT's full_flag is the AND of r1_full/r2_full");
    assert_eq!(f1, either, "every full R1 row is an 'either' pair");
    assert_eq!(f2, 0, "the fixture never makes R2 full");
    assert_eq!(both, 0, "no pair has both mates full");
    assert!(either > 100, "the fixture must produce a visible either-rate: {either}");

    assert_eq!(r.report["metrics"]["full_artefact_pairs"], both);
    assert_eq!(r.report["metrics"]["full_artefact_pairs_either"], either);
    assert_eq!(r.report["metrics"]["full_artefact_rate"], 0.0);
    assert!(
        r.stderr.contains("full-artefact pairs (GT full_flag)"),
        "the summary must name which count it prints: {}",
        r.stderr
    );

    // Single-end: the pair IS R1, so the GT's full_flag is r1_full and the
    // both-mate counter must follow it.
    let s = run(&dir, "single", &p, 1000, &["--single", "--dup-rate", "0"]);
    let rows = gt_rows(&s.out);
    let single_full = rows.iter().filter(|x| x[8] == "1").count();
    assert!(single_full > 50, "single-end fixture must produce full reads");
    assert_eq!(s.report["metrics"]["full_artefact_pairs"], single_full);
    assert_eq!(s.report["metrics"]["full_artefact_pairs_either"], single_full);

    std::fs::remove_dir_all(&dir).unwrap();
}
