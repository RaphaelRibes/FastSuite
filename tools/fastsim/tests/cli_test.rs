//! CLI-level integration tests: parse args, run the binary end to end,
//! check the three output files and their consistency.

use std::path::Path;
use std::process::Command;

/// Minimal valid profile (150 bp, no artefacts, adapter off) written to a
/// temp file. Returns the path.
fn write_test_profile(dir: &Path) -> std::path::PathBuf {
    let q150 = vec!["30"; 150].join(",");
    let json = format!(
        r#"{{
        "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
        "head": {{"R1": {{"0": 100, "1": 50}}, "R2": {{"0": 100, "1": 40, "150": 10}}}},
        "tail": {{"R1": {{"0": 60, "1": 30, "2": 10}}, "R2": {{"0": 60, "1": 30, "2": 10}}}},
        "joint_tail": {{"0,0": 36, "1,1": 12, "2,2": 2}},
        "tail_q_collapse": {{"R1": 10.0, "R2": 2.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.01,
        "full_length_q": {{"R1": 20, "R2": 20}},
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"}}
    }}"#,
        q = q150
    );
    let p = dir.join("test_profile.json");
    std::fs::write(&p, json).unwrap();
    p
}

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastsim")
}

#[test]
fn smoke_run_produces_three_files() {
    let dir = std::env::temp_dir().join(format!("fastsim_it_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("sim");

    let status = Command::new(bin())
        .args([
            "--pairs",
            "1000",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--seed",
            "42",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let r1 = dir.join("sim_R1.fastq.gz");
    let r2 = dir.join("sim_R2.fastq.gz");
    let gt = dir.join("sim.gt.tsv");
    assert!(r1.exists(), "R1 missing");
    assert!(r2.exists(), "R2 missing");
    assert!(gt.exists(), "GT missing");

    // GT: header + 1000 rows.
    let gt_text = std::fs::read_to_string(&gt).unwrap();
    assert_eq!(gt_text.lines().count(), 1001);
    let header = gt_text.lines().next().unwrap();
    assert_eq!(header, "pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\tr1_n\tr2_n");

    // Cleanup.
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn single_end_writes_one_fastq_and_single_read_gt() {
        let dir = std::env::temp_dir().join(format!("fastsim_se_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let profile = write_test_profile(&dir);
        let out = dir.join("se");
        let status = Command::new(env!("CARGO_BIN_EXE_fastsim"))
            .args([
                "--pairs", "500", "--profile", profile.to_str().unwrap(),
                "--out", out.to_str().unwrap(), "--single", "--no-gzip",
            ])
            .status()
            .unwrap();
        assert!(status.success());

        // Only R1 is written; no R2.
        assert!(dir.join("se_R1.fastq").exists(), "R1 missing");
        assert!(!dir.join("se_R2.fastq").exists(), "R2 must not exist in single-end");

        // GT uses the single-read schema and has 500 rows + header.
        let gt = std::fs::read_to_string(dir.join("se.gt.tsv")).unwrap();
        let lines: Vec<&str> = gt.lines().collect();
        assert_eq!(lines.len(), 501);
        assert_eq!(
            lines[0],
            "pair_id\tinsert\tr1_h\tr1_t\tr1_adapter\tfull_flag\tis_dup\tr1_n"
        );
        // every row has 8 fields
        assert!(lines[1..].iter().all(|l| l.split('\t').count() == 8));

        // Adapter readthrough is labelled: some rows carry a nonzero r1_adapter
        // (the test profile has readthrough mode and short inserts).
        let with_adapter = lines[1..]
            .iter()
            .filter(|l| l.split('\t').nth(4).unwrap().parse::<u32>().unwrap() > 0)
            .count();
        assert!(with_adapter > 0, "no adapter-labelled rows in single-end GT");

        std::fs::remove_dir_all(&dir).unwrap();
    }

#[test]
fn gt_and_fastq_agree_on_full_flag_rate() {
    let dir = std::env::temp_dir().join(format!("fastsim_it2_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("sim");

    let status = Command::new(bin())
        .args([
            "--pairs",
            "10000",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--seed",
            "7",
            "--no-gzip",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // Plain FASTQ: count records (4 lines each).
    let r1 = std::fs::read_to_string(dir.join("sim_R1.fastq")).unwrap();
    let r2 = std::fs::read_to_string(dir.join("sim_R2.fastq")).unwrap();
    assert_eq!(r1.lines().count(), 40_000, "R1 records");
    assert_eq!(r2.lines().count(), 40_000, "R2 records");

    // GT full flag: R2 fixture heads reach 150 (10% of reads) -> R2 full
    // rate ~10%; R1 heads max 1 -> R1 never full (per-mate full-length).
    let gt_text = std::fs::read_to_string(dir.join("sim.gt.tsv")).unwrap();
    let mut full_r1 = 0usize;
    let mut full_r2 = 0usize;
    for line in gt_text.lines().skip(1) {
        let f = line.split('\t').collect::<Vec<_>>();
        let (r1h, r1t, r2h, r2t) = (
            f[2].parse::<u32>().unwrap(),
            f[3].parse::<u32>().unwrap(),
            f[4].parse::<u32>().unwrap(),
            f[5].parse::<u32>().unwrap(),
        );
        if r1h + r1t >= 150 {
            full_r1 += 1;
        }
        if r2h + r2t >= 150 {
            full_r2 += 1;
        }
    }
    assert_eq!(full_r1, 0, "R1 must never be full (head max 1)");
    let rate = full_r2 as f64 / 10_000.0;
    // R2 head CDF: P(h=150) = 10/(100+40+10) = 6.7%
    assert!(
        (rate - 0.067).abs() < 0.02,
        "R2 full rate {rate} (expect ~6.7%)"
    );

    // Every R1/R2 record is exactly 150 bases.
    for line in r1.lines().skip(1).step_by(4) {
        assert_eq!(line.len(), 150);
    }
    for line in r2.lines().skip(1).step_by(4) {
        assert_eq!(line.len(), 150);
    }

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn zero_pairs_is_an_error() {
    let status = Command::new(bin())
        .args(["--pairs", "0", "--profile", "x.json", "--out", "x"])
        .output()
        .unwrap();
    assert!(!status.status.success());
    let err = String::from_utf8_lossy(&status.stderr);
    assert!(err.contains("--pairs must be >= 1"), "stderr: {err}");
}

#[test]
fn missing_profile_is_an_error() {
    let status = Command::new(bin())
        .args(["--pairs", "10", "--out", "x"])
        .output()
        .unwrap();
    assert!(!status.status.success());
    let err = String::from_utf8_lossy(&status.stderr);
    // clap reports a missing required argument with its own wording; the
    // contract is "non-zero exit + the flag named on stderr".
    assert!(err.contains("--profile"), "stderr: {err}");
}

#[test]
fn unknown_flag_is_an_error() {
    let status = Command::new(bin())
        .args([
            "--pairs",
            "10",
            "--profile",
            "p.json",
            "--out",
            "x",
            "--bogus",
        ])
        .output()
        .unwrap();
    assert!(!status.status.success());
    let err = String::from_utf8_lossy(&status.stderr);
    // Suite error style (clap): "unexpected argument" + usage, exit 2.
    assert!(err.contains("unexpected argument"), "stderr: {err}");
}

#[test]
fn deterministic_across_runs() {
    let dir = std::env::temp_dir().join(format!("fastsim_it3_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    for n in ["a", "b"] {
        let out = dir.join(format!("sim_{n}"));
        let status = Command::new(bin())
            .args([
                "--pairs",
                "500",
                "--profile",
                profile.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
                "--seed",
                "99",
                "--no-gzip",
            ])
            .status()
            .unwrap();
        assert!(status.success());
    }
    let a = std::fs::read(dir.join("sim_a_R1.fastq")).unwrap();
    let b = std::fs::read(dir.join("sim_b_R1.fastq")).unwrap();
    assert_eq!(a, b, "same seed must produce identical R1 bytes");

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn dup_rate_reproduces_and_dups_copy_a_stored_template() {
    let dir = std::env::temp_dir().join(format!("fastsim_it4_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("sim");

    let status = Command::new(bin())
        .args([
            "--pairs",
            "20000",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--seed",
            "7",
            "--no-gzip",
            "--dup-rate",
            "0.10",
            "--dup-pool",
            "4096",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // 1. Rate: GT is_dup column should be ~10%.
    let gt = std::fs::read_to_string(dir.join("sim.gt.tsv")).unwrap();
    let mut n = 0u64;
    let mut nd = 0u64;
    for line in gt.lines().skip(1) {
        n += 1;
        let dup: u8 = line.split('\t').nth(9).unwrap().parse().unwrap();
        nd += u64::from(dup);
    }
    let rate = nd as f64 / n as f64;
    assert!((rate - 0.10).abs() < 0.01, "dup rate {rate} != 0.10");

    // 2. Same-template property: a duplicate re-sequences a previously
    // emitted insert, so its read shares the insert core with an earlier
    // read. Check via 60-mers: any 60-mer of the dup's R1 should appear in
    // the set of 60-mers of earlier reads. A random 60-mer collides at
    // ~4^-60, and with Q30 a shared 100+ bp core yields ~40 60-mers, so a
    // false positive is impossible and a true duplicate almost always
    // matches (each 60-mer survives its bases error-free with ~94% prob).
    let r1 = std::fs::read_to_string(dir.join("sim_R1.fastq")).unwrap();
    let seqs: Vec<&str> = r1.lines().skip(1).step_by(4).collect();
    assert_eq!(seqs.len() as u64, n, "one seq per GT row");
    let mut seen: std::collections::HashSet<&[u8]> = std::collections::HashSet::new();
    let mut dup_copies = 0usize;
    let mut matched = 0usize;
    for (i, line) in gt.lines().skip(1).enumerate() {
        let is_dup = line.split('\t').nth(9) == Some("1");
        let s = seqs[i].as_bytes();
        if is_dup {
            dup_copies += 1;
            if s.windows(60).any(|k| seen.contains(k)) {
                matched += 1;
            }
        }
        for k in s.windows(60) {
            seen.insert(k);
        }
    }
    assert!(dup_copies > 0, "no duplicates drawn at 10%");
    assert!(
        matched * 100 >= dup_copies * 90,
        "only {matched}/{dup_copies} dups share a 60-mer with an earlier read"
    );

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn drop_rate_produces_short_fragments_and_reads_stay_150bp() {
    let dir = std::env::temp_dir().join(format!("fastsim_it5_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("sim");

    let status = Command::new(bin())
        .args([
            "--pairs",
            "20000",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--seed",
            "11",
            "--no-gzip",
            "--drop-rate",
            "0.05",
            "--drop-min-len",
            "100",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // 1. Reads are exactly 150 bp — every one of them (the sequencer emits
    // read_len bases regardless of fragment length).
    for r in ["sim_R1.fastq", "sim_R2.fastq"] {
        let text = std::fs::read_to_string(dir.join(r)).unwrap();
        let seqs: Vec<&str> = text.lines().skip(1).step_by(4).collect();
        assert_eq!(seqs.len(), 20000, "{r}: one record per pair");
        for (i, s) in seqs.iter().enumerate() {
            assert_eq!(s.len(), 150, "{r}: record {i} has len {}", s.len());
        }
    }

    // 2. Drop class: ~5% of GT rows have insert < 100 (forced short
    // fragments), which makes the keep window [h, min(I, L-t)) < min_len on
    // both mates.
    let gt = std::fs::read_to_string(dir.join("sim.gt.tsv")).unwrap();
    let mut n = 0u64;
    let mut short = 0u64;
    for line in gt.lines().skip(1) {
        let p: Vec<&str> = line.split('\t').collect();
        n += 1;
        let insert: u32 = p[1].parse().unwrap();
        if insert < 100 {
            short += 1;
        }
    }
    let rate = short as f64 / n as f64;
    assert!((rate - 0.05).abs() < 0.01, "drop rate {rate} != 0.05");

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn profile_dup_rate_beats_the_cli_flag() {
    // main.rs does `profile.dup_rate.unwrap_or(cfg.dup_rate)`, and clap cannot
    // tell an explicit --dup-rate from its 0.02 default, so once the profile
    // carries a measured dup_rate the CLI flag is NEVER consulted. This pins
    // that precedence: profile dup_rate 0.25 + `--dup-rate 0.05` must yield a
    // GT duplicate rate that follows the PROFILE (~0.25), not the flag. A
    // harness pinning --dup-rate for reproducibility would otherwise silently
    // get the measured rate.
    let dir = std::env::temp_dir().join(format!("fastsim_it_duprate_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let q150 = vec!["30"; 150].join(",");
    let json = format!(
        r#"{{
        "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
        "dup_rate": 0.25,
        "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "joint_tail": {{"0,0": 100}},
        "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.0,
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"}}
    }}"#,
        q = q150
    );
    let profile = dir.join("duprate_profile.json");
    std::fs::write(&profile, json).unwrap();
    let out = dir.join("sim");

    let status = Command::new(bin())
        .args([
            "--pairs", "20000", "--profile", profile.to_str().unwrap(),
            "--out", out.to_str().unwrap(), "--seed", "7", "--no-gzip",
            "--dup-rate", "0.05", "--dup-pool", "4096",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // GT is_dup column (index 9) must track the profile's 0.25, not the flag.
    let gt = std::fs::read_to_string(dir.join("sim.gt.tsv")).unwrap();
    let mut n = 0u64;
    let mut nd = 0u64;
    for line in gt.lines().skip(1) {
        n += 1;
        let dup: u8 = line.split('\t').nth(9).unwrap().parse().unwrap();
        nd += u64::from(dup);
    }
    let rate = nd as f64 / n as f64;
    assert!(
        (rate - 0.25).abs() < 0.02,
        "observed dup rate {rate} must follow the profile's 0.25, not --dup-rate 0.05"
    );

    std::fs::remove_dir_all(&dir).unwrap();
}
