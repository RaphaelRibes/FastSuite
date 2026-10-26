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

/// Profile with a configurable read length (no insert CDF, adapter off), for
/// the short-read closed-loop test.
fn write_short_profile(dir: &Path, read_len: u32) -> std::path::PathBuf {
    let q = vec!["30"; read_len as usize].join(",");
    let json = format!(
        r#"{{
        "read_len": {l}, "gc": 0.5, "geometry": "tail3p_both",
        "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "tail": {{"R1": {{"0": 70, "2": 30}}, "R2": {{"0": 70, "2": 30}}}},
        "joint_tail": {{"0,0": 49, "2,2": 9}},
        "tail_q_collapse": {{"R1": 2.0, "R2": 1.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.0,
        "full_length_q": {{"R1": 20, "R2": 20}},
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
    }}"#,
        l = read_len,
        q = q
    );
    let p = dir.join(format!("profile_{read_len}.json"));
    std::fs::write(&p, json).unwrap();
    p
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
    assert_eq!(header, "pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\tr1_n\tr2_n\tr1_full\tr2_full");

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
            "pair_id\tinsert\tr1_h\tr1_t\tr1_adapter\tfull_flag\tis_dup\tr1_n\tr1_full"
        );
        // every row has 9 fields
        assert!(lines[1..].iter().all(|l| l.split('\t').count() == 9));

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

#[test]
fn short_read_profile_runs_without_insert_min() {
    // The documented closed loop is `fastcheck --profile p.json` then
    // `fastsim --profile p.json` and nothing else. A fixed 150 bp insert floor
    // refused every library with read_len <= 150 (fastsim requires
    // insert-min <= read_len); the floor now defaults to min(150, read_len),
    // so 100 bp and 75 bp profiles go straight in.
    for l in [100u32, 75] {
        let dir = std::env::temp_dir()
            .join(format!("fastsim_it_short{l}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let profile = write_short_profile(&dir, l);
        let out = dir.join("sim");
        let o = Command::new(bin())
            .args([
                "--pairs", "500", "--profile", profile.to_str().unwrap(),
                "--out", out.to_str().unwrap(), "--no-gzip", "--drop-rate", "0",
            ])
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "read_len {l} must simulate with the default insert floor; stderr: {}",
            String::from_utf8_lossy(&o.stderr)
        );

        // Every read is exactly read_len bases...
        let r1 = std::fs::read_to_string(dir.join("sim_R1.fastq")).unwrap();
        for s in r1.lines().skip(1).step_by(4) {
            assert_eq!(s.len(), l as usize, "read_len {l}: read is not full length");
        }
        // ...and the floor actually applied is the profile's read_len
        // (min(150, l) == l), neither 0 nor 150.
        let gt = std::fs::read_to_string(dir.join("sim.gt.tsv")).unwrap();
        let min_insert = gt
            .lines()
            .skip(1)
            .map(|line| line.split('\t').nth(1).unwrap().parse::<u32>().unwrap())
            .min()
            .unwrap();
        assert!(min_insert >= l, "read_len {l}: shortest insert {min_insert}");

        std::fs::remove_dir_all(&dir).unwrap();
    }
}

#[test]
fn explicit_insert_min_still_wins() {
    // The profile-derived default must not swallow an explicit flag: a value
    // above the profile's read_len is still refused, and the message names the
    // value the user passed (proving the flag, not the default, was used).
    let dir = std::env::temp_dir().join(format!("fastsim_it_imin_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_short_profile(&dir, 100);
    let out = dir.join("sim");
    let o = Command::new(bin())
        .args([
            "--pairs", "10", "--profile", profile.to_str().unwrap(),
            "--out", out.to_str().unwrap(), "--insert-min", "200",
        ])
        .output()
        .unwrap();
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--insert-min (200)"), "stderr: {err}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn zero_dup_pool_is_a_clean_cli_error() {
    // --dup-pool 0 used to reach History::push and panic (index out of bounds
    // on the empty buffer, then a modulo by zero); it must be a parse-time
    // error instead.
    let o = Command::new(bin())
        .args([
            "--pairs", "10", "--profile", "p.json", "--out", "x",
            "--dup-pool", "0",
        ])
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert_ne!(o.status.code(), Some(101), "must not be a panic");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--dup-pool must be >= 1"), "stderr: {err}");
    assert!(!err.contains("panicked"), "stderr: {err}");
}

/// Profile with a configurable `gc` (no insert CDF, adapter off, no artefacts).
fn gc_profile(dir: &Path, gc: f64) -> std::path::PathBuf {
    let q = vec!["30"; 150].join(",");
    let json = format!(
        r#"{{
        "read_len": 150, "gc": {gc}, "geometry": "tail3p_both",
        "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "joint_tail": {{"0,0": 100}},
        "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.0,
        "full_length_q": {{"R1": 20, "R2": 20}},
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
    }}"#,
        gc = gc,
        q = q
    );
    let p = dir.join(format!("gc_{}.json", (gc * 1000.0) as u32));
    std::fs::write(&p, json).unwrap();
    p
}

/// Emitted GC of a plain FASTQ's R1 records.
fn emitted_gc(path: &Path) -> f64 {
    let s = std::fs::read_to_string(path).unwrap();
    let mut gc = 0usize;
    let mut total = 0usize;
    for (i, line) in s.lines().enumerate() {
        if i % 4 != 1 {
            continue;
        }
        for b in line.bytes() {
            total += 1;
            if b == b'G' || b == b'C' {
                gc += 1;
            }
        }
    }
    gc as f64 / total as f64
}

/// A profile's `gc` is the authority; `--gc` only overrides it. Before this
/// pin, both were silently replaced by the flag's 0.49 default (a gc=0.30
/// profile emitted 0.4992, so the `gc` field was inert).
#[test]
fn profile_gc_drives_the_generator_and_the_flag_overrides_it() {
    let dir = std::env::temp_dir().join(format!("fastsim_gc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = gc_profile(&dir, 0.30);

    let run = |out: &Path, extra: &[&str]| {
        let o = Command::new(bin())
            .args([
                "--pairs",
                "2000",
                "--profile",
                profile.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
                "--seed",
                "5",
                "--no-gzip",
                "--no-adapter",
                "--drop-rate",
                "0",
                "--quiet",
            ])
            .args(extra)
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "fastsim {extra:?} failed: {}",
            String::from_utf8_lossy(&o.stderr)
        );
    };

    let from_profile = dir.join("p");
    run(&from_profile, &[]);
    let observed = emitted_gc(&std::path::PathBuf::from(format!(
        "{}_R1.fastq",
        from_profile.display()
    )));
    assert!(
        (observed - 0.30).abs() < 0.02,
        "profile gc=0.30 emitted {observed}"
    );

    // The flag still wins when given.
    let from_flag = dir.join("f");
    run(&from_flag, &["--gc", "0.60"]);
    let observed = emitted_gc(&std::path::PathBuf::from(format!(
        "{}_R1.fastq",
        from_flag.display()
    )));
    assert!(
        (observed - 0.60).abs() < 0.02,
        "--gc 0.60 emitted {observed}"
    );

    std::fs::remove_dir_all(&dir).unwrap();
}

/// A profile with an empty adapter sequence must be a clean refusal or a
/// no-adapter run, never a modulo-by-zero panic (exit 101).
#[test]
fn empty_profile_adapter_sequence_does_not_panic() {
    let dir = std::env::temp_dir().join(format!("fastsim_epad_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let q = vec!["30"; 150].join(",");
    let json = format!(
        r#"{{
        "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
        "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "joint_tail": {{"0,0": 100}},
        "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.0,
        "full_length_q": {{"R1": 20, "R2": 20}},
        "adapter": {{"seq": "", "mode": "readthrough"}}
    }}"#,
        q = q
    );
    let profile = dir.join("empty_adapter.json");
    std::fs::write(&profile, json).unwrap();
    // --insert-mean 60 forces readthrough: the window is longer than the insert.
    let o = Command::new(bin())
        .args([
            "--pairs",
            "200",
            "--insert-mean",
            "60",
            "--insert-sd",
            "10",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            dir.join("e").to_str().unwrap(),
            "--no-gzip",
            "--quiet",
        ])
        .output()
        .unwrap();
    assert_ne!(o.status.code(), Some(101), "must not be a panic");
    assert!(o.status.success(), "empty adapter should mean 'no adapter'");
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `is_dup=1` must always name a payload that exists earlier in the same file:
/// with an empty duplicate pool the pair is NOT a duplicate.
#[test]
fn duplicate_flag_always_has_an_earlier_copy() {
    let dir = std::env::temp_dir().join(format!("fastsim_dupflag_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("d");
    let o = Command::new(bin())
        .args([
            "--pairs",
            "40",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--seed",
            "3",
            "--dup-rate",
            "1.0",
            "--no-gzip",
            "--quiet",
        ])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let r1 = std::fs::read_to_string(dir.join("d_R1.fastq")).unwrap();
    let seqs: Vec<&str> = r1.lines().skip(1).step_by(4).collect();
    let gt = std::fs::read_to_string(dir.join("d.gt.tsv")).unwrap();
    for (i, line) in gt.lines().skip(1).enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        if f[9] == "1" {
            assert!(
                seqs[..i].contains(&seqs[i]),
                "pair {i} is flagged duplicate with no earlier byte-identical record"
            );
        }
    }
    // Pair 0 can never be a duplicate (nothing precedes it).
    assert_eq!(gt.lines().nth(1).unwrap().split('\t').nth(9).unwrap(), "0");
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `--single` writes no R2, so the report must not name one and must count one
/// read per pair (it used to claim 2*pairs and a nonexistent R2 path).
#[test]
fn single_mode_report_matches_the_files() {
    let dir = std::env::temp_dir().join(format!("fastsim_serep_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("s");
    let o = Command::new(bin())
        .args([
            "--pairs",
            "100",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--single",
            "--no-gzip",
            "--json",
            dir.join("r.json").to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let rep: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("r.json")).unwrap()).unwrap();
    assert!(rep["outputs"]["r2"].is_null(), "single mode names an R2: {rep}");
    assert_eq!(rep["metrics"]["reads_in"], 100, "reads_in must be one read per pair");
    assert!(!dir.join("s_R2.fastq").exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// S6: a profile with `read_len` 0 used to be accepted and wrote structurally
/// invalid FASTQ (empty sequence and quality lines) with exit 0. A zero-length
/// read is not a library; it must be refused before any output is opened.
#[test]
fn zero_read_len_profile_is_refused() {
    let dir = std::env::temp_dir().join(format!("fastsim_rl0_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let json = r#"{
        "read_len": 0, "gc": 0.5, "geometry": "tail3p_both",
        "head": {"R1": {"0": 100}, "R2": {"0": 100}},
        "tail": {"R1": {"0": 100}, "R2": {"0": 100}},
        "joint_tail": {"0,0": 100},
        "tail_q_collapse": {"R1": 0.0, "R2": 0.0},
        "pos_q": {"R1": [], "R2": []},
        "read_meanq": {"R1": {"30": 100}, "R2": {"30": 100}},
        "full_length_rate": 0.0,
        "full_length_q": {"R1": 20, "R2": 20},
        "adapter": {"seq": "CTGTCTCTTATACACATCT", "mode": "none"}
    }"#;
    let profile = dir.join("zero.json");
    std::fs::write(&profile, json).unwrap();

    let out = dir.join("z");
    let o = Command::new(bin())
        .args([
            "--pairs", "10", "--profile", profile.to_str().unwrap(),
            "--out", out.to_str().unwrap(), "--no-gzip",
        ])
        .output()
        .unwrap();
    assert!(!o.status.success(), "read_len 0 must be refused");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("read_len"), "the error must name read_len: {err}");
    assert!(
        !dir.join("z_R1.fastq").exists(),
        "a refused profile must not leave output behind"
    );

    std::fs::remove_dir_all(&dir).unwrap();
}

/// Cross-cutting: `-o` is a prefix, so its parent directories are created (a
/// harness writing into a fresh per-run directory should not have to mkdir).
#[test]
fn out_prefix_creates_missing_parent_directories() {
    let dir = std::env::temp_dir().join(format!("fastsim_mkdir_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    let out = dir.join("a/b/c/sim");
    assert!(!dir.join("a").exists(), "the parent must really be absent");
    let o = Command::new(bin())
        .args([
            "--pairs", "50", "--profile", profile.to_str().unwrap(),
            "--out", out.to_str().unwrap(), "--no-gzip", "-q",
        ])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "nested --out prefix: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    for suffix in ["_R1.fastq", "_R2.fastq", ".gt.tsv"] {
        assert!(dir.join(format!("a/b/c/sim{suffix}")).exists(), "{suffix} missing");
    }

    std::fs::remove_dir_all(&dir).unwrap();
}
