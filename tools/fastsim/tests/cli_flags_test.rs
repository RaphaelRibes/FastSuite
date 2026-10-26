//! fastsim's CLI flag contract, in one file: the FastSuite common flags
//! (`--force` and the outputs it guards, `-q/--quiet`), the flags a measured
//! profile makes inert, the `FASTSUITE_THREADS` env default, and the
//! clap-generated help/version surface. The parse-error and smoke-level
//! cases live in `cli_test.rs`.

use std::path::Path;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastsim")
}

/// Minimal valid profile (150 bp, no artefacts, adapter on) written to a
/// temp file. Returns the path.
fn write_test_profile(dir: &Path) -> std::path::PathBuf {
    let q150 = vec!["30"; 150].join(",");
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
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"}}
    }}"#,
        q = q150
    );
    let p = dir.join("flags_profile.json");
    std::fs::write(&p, json).unwrap();
    p
}

#[test]
fn force_refuses_to_clobber_without_the_flag() {
    let dir = std::env::temp_dir().join(format!("fastsim_force_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("f");
    let args = ["--pairs", "100", "--profile", profile.to_str().unwrap(),
                "--out", out.to_str().unwrap()];

    let first = Command::new(bin()).args(args).output().unwrap();
    assert!(first.status.success(), "first run must succeed");

    let second = Command::new(bin()).args(args).output().unwrap();
    assert!(!second.status.success(), "clobbered an existing prefix without --force");
    let err = String::from_utf8_lossy(&second.stderr);
    assert!(err.contains("--force"), "stderr must name --force: {err}");

    let third = Command::new(bin()).args(args).arg("--force").output().unwrap();
    assert!(third.status.success(), "--force must overwrite");

    // The short form -f works identically.
    let out2 = dir.join("f2");
    let a = ["--pairs", "100", "--profile", profile.to_str().unwrap(),
             "--out", out2.to_str().unwrap()];
    assert!(Command::new(bin()).args(a).output().unwrap().status.success());
    let again = Command::new(bin()).args(a).arg("-f").output().unwrap();
    assert!(again.status.success(), "-f must overwrite");

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn json_report_is_guarded_by_force_like_every_other_output() {
    // --json is an output the run creates, so the clobber guard must probe it
    // too: it used to check only the FASTQ/GT suffixes and replaced an existing
    // report with no --force and exit 0.
    let dir = std::env::temp_dir().join(format!("fastsim_force_json_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("j");
    let report = dir.join("j-report.json");
    std::fs::write(&report, "SENTINEL\n").unwrap();
    let args = [
        "--pairs",
        "50",
        "--profile",
        profile.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--json",
        report.to_str().unwrap(),
        "--no-gzip",
    ];

    let refused = Command::new(bin()).args(args).output().unwrap();
    assert!(
        !refused.status.success(),
        "an existing --json report was clobbered without --force"
    );
    let err = String::from_utf8_lossy(&refused.stderr);
    assert!(err.contains("--force"), "stderr must name --force: {err}");
    assert!(
        err.contains(&report.display().to_string()),
        "the refusal must name the report it protects: {err}"
    );
    assert_eq!(
        std::fs::read_to_string(&report).unwrap(),
        "SENTINEL\n",
        "the sentinel report must be untouched"
    );
    assert!(
        !dir.join("j.gt.tsv").exists(),
        "the guard must fire before any output is written"
    );

    let forced = Command::new(bin()).args(args).arg("--force").output().unwrap();
    assert!(
        forced.status.success(),
        "--force must overwrite the report: {}",
        String::from_utf8_lossy(&forced.stderr)
    );
    let written: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&report).unwrap()).unwrap();
    assert_eq!(written["tool"], "fastsim", "the report was replaced");

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn quiet_suppresses_the_summary_line() {
    let dir = std::env::temp_dir().join(format!("fastsim_quiet_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    let out = dir.join("q");
    let res = Command::new(bin())
        .args(["--pairs", "100", "--profile", profile.to_str().unwrap(),
               "--out", out.to_str().unwrap(), "-q"])
        .output()
        .unwrap();
    assert!(res.status.success());
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(!stderr.contains("fastsim: 100 pairs"), "summary leaked under -q: {stderr}");
    // The files are still written — -q only gates the summary line.
    assert!(dir.join("q.gt.tsv").exists(), "outputs missing under -q");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn fastsuite_threads_env_defaults_the_thread_count() {
    let dir = std::env::temp_dir().join(format!("fastsim_env_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    // FASTSUITE_THREADS=4, no -j: the (ignored) threads warning must name 4.
    let out = dir.join("e1");
    let res = Command::new(bin())
        .args(["--pairs", "10", "--profile", profile.to_str().unwrap(),
               "--out", out.to_str().unwrap()])
        .env("FASTSUITE_THREADS", "4")
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("--threads 4"), "env default not applied: {err}");

    // An explicit -j 1 beats the env var: no warning at all.
    let out2 = dir.join("e2");
    let res2 = Command::new(bin())
        .args(["--pairs", "10", "--profile", profile.to_str().unwrap(),
               "--out", out2.to_str().unwrap(), "-j", "1"])
        .env("FASTSUITE_THREADS", "4")
        .output()
        .unwrap();
    let err2 = String::from_utf8_lossy(&res2.stderr);
    assert!(!err2.contains("--threads"), "explicit -j was overridden: {err2}");

    // 0 / garbage in the env var is ignored (default 1).
    for bad in ["0", "not-a-number"] {
        let out3 = dir.join(format!("e_{bad}"));
        let res3 = Command::new(bin())
            .args(["--pairs", "10", "--profile", profile.to_str().unwrap(),
                   "--out", out3.to_str().unwrap()])
            .env("FASTSUITE_THREADS", bad)
            .output()
            .unwrap();
        let err3 = String::from_utf8_lossy(&res3.stderr);
        assert!(!err3.contains("--threads"), "{bad} must be ignored: {err3}");
    }

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn version_and_help_exit_zero() {
    for flag in ["-V", "--version"] {
        let out = Command::new(bin()).arg(flag).output().unwrap();
        assert!(out.status.success(), "{flag} must exit 0");
        let s = String::from_utf8_lossy(&out.stdout);
        assert!(s.starts_with("fastsim 0."), "stdout: {s}");
    }
    for flag in ["-h", "--help"] {
        let out = Command::new(bin()).arg(flag).output().unwrap();
        assert!(out.status.success(), "{flag} must exit 0");
        let s = String::from_utf8_lossy(&out.stdout);
        assert!(s.contains("--pairs") && s.contains("--force") && s.contains("--quiet"),
                "help missing a canonical flag: {s}");
    }
    // The compression level is a first-class flag, not a hidden constant.
    let help = Command::new(bin()).arg("--help").output().unwrap();
    let s = String::from_utf8_lossy(&help.stdout);
    assert!(s.contains("--gz-level"), "help must document --gz-level: {s}");
}

/// FS-1: `--gz-level` defaults to 6 (today's bytes), is threaded to both FASTQ
/// writers, and is recorded in the `--json` report.
#[test]
fn gz_level_defaults_to_six_and_is_reported() {
    let dir = std::env::temp_dir().join(format!("fastsim_gzlvl_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    let run = |out: &Path, extra: &[&str]| -> serde_json::Value {
        let rep = dir.join(format!("{}.json", out.display()));
        let mut args = vec![
            "--pairs".to_string(),
            "200".to_string(),
            "--profile".to_string(),
            profile.to_str().unwrap().to_string(),
            "--out".to_string(),
            out.to_str().unwrap().to_string(),
            "--json".to_string(),
            rep.to_str().unwrap().to_string(),
            "--force".to_string(),
            "--quiet".to_string(),
        ];
        args.extend(extra.iter().map(|s| s.to_string()));
        let res = Command::new(bin()).args(&args).output().unwrap();
        assert!(res.status.success(), "run {extra:?}: {}", String::from_utf8_lossy(&res.stderr));
        serde_json::from_slice(&std::fs::read(&rep).unwrap()).unwrap()
    };

    let dflt = run(&dir.join("dflt"), &[]);
    assert_eq!(dflt["metrics"]["gz_level"], 6, "default level must be 6");

    let one = run(&dir.join("one"), &["--gz-level", "1"]);
    assert_eq!(one["metrics"]["gz_level"], 1, "--gz-level must reach the report");

    // The smaller level is the larger file: the documented trade is real.
    let size = |prefix: &Path| {
        std::fs::metadata(format!("{}_R1.fastq.gz", prefix.display()))
            .unwrap()
            .len()
    };
    assert!(
        size(&dir.join("one")) > size(&dir.join("dflt")),
        "level 1 must trade size for speed"
    );

    let plain = run(&dir.join("plain"), &["--no-gzip"]);
    assert_eq!(plain["metrics"]["gz_bytes"], 0, "no gzip -> no gzip bytes");

    std::fs::remove_dir_all(&dir).unwrap();
}

/// FS-1: a level outside 1..=9 is refused at parse time, not silently clamped.
#[test]
fn gz_level_out_of_range_is_rejected() {
    let dir = std::env::temp_dir().join(format!("fastsim_gzlvl_bad_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);
    for bad in ["0", "10"] {
        let out = dir.join(format!("bad{bad}"));
        let res = Command::new(bin())
            .args([
                "--pairs", "10", "--profile", profile.to_str().unwrap(),
                "--out", out.to_str().unwrap(), "--gz-level", bad,
            ])
            .output()
            .unwrap();
        assert!(!res.status.success(), "--gz-level {bad} must be rejected");
        let err = String::from_utf8_lossy(&res.stderr);
        assert!(err.contains("gz-level"), "error must name the flag: {err}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// FS-11: the instrument counters (gz_bytes, substitutions, mean_q, rng_draws)
/// are reported, and enabling the draw counter must not move a single byte.
#[test]
fn instrument_counters_are_reported_and_stream_neutral() {
    let dir = std::env::temp_dir().join(format!("fastsim_stats_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    let run = |out: &Path, extra: &[&str]| -> serde_json::Value {
        let rep = dir.join(format!("{}.json", out.display()));
        let mut args = vec![
            "--pairs".to_string(), "300".to_string(),
            "--profile".to_string(), profile.to_str().unwrap().to_string(),
            "--out".to_string(), out.to_str().unwrap().to_string(),
            "--json".to_string(), rep.to_str().unwrap().to_string(),
            "--force".to_string(), "--quiet".to_string(),
            "--no-gzip".to_string(),
        ];
        args.extend(extra.iter().map(|s| s.to_string()));
        let res = Command::new(bin()).args(&args).output().unwrap();
        assert!(res.status.success(), "run {extra:?}: {}", String::from_utf8_lossy(&res.stderr));
        serde_json::from_slice(&std::fs::read(&rep).unwrap()).unwrap()
    };

    let quiet = run(&dir.join("quiet"), &[]);
    let counted = run(&dir.join("counted"), &["--stats"]);

    // Counting draws must not consume randomness: the GT is byte-identical.
    let gt = |out: &Path| std::fs::read(format!("{}.gt.tsv", out.display())).unwrap();
    assert_eq!(gt(&dir.join("quiet")), gt(&dir.join("counted")), "--stats moved the stream");

    // Task 2.3b: and neither is the READ STREAM. The GT alone is a weaker check
    // than it looks -- it records the simulated truth, so a change confined to
    // the emitted FASTQ (a different quality draw, say) could leave it equal.
    // `TrackedRng` is meant to be a pure pass-through, so both mates must match
    // byte for byte too. Measured on 200k pairs at -j 8: R1, R2 and the GT are
    // all identical, and CPU is identical to three decimals (1.345 s both).
    let reads = |out: &Path, mate: &str| {
        std::fs::read(format!("{}_{mate}.fastq", out.display())).unwrap()
    };
    for mate in ["R1", "R2"] {
        assert_eq!(
            reads(&dir.join("quiet"), mate),
            reads(&dir.join("counted"), mate),
            "--stats changed the emitted {mate} reads, not just the ground truth"
        );
    }

    assert_eq!(quiet["metrics"]["rng_draws"], 0, "draw counting is off by default");
    assert!(
        counted["metrics"]["rng_draws"].as_u64().unwrap() > 0,
        "--stats must count draws"
    );
    // ~300 bases/pair x 300 pairs: substitutions are rare but mean_q is not.
    let mq = counted["metrics"]["mean_q"].as_f64().unwrap();
    assert!(mq > 0.0 && mq <= 41.0, "mean_q {mq} out of range");
    assert!(counted["metrics"].get("substitutions").is_some());

    // mean_q must be the mean of the quality bytes actually emitted.
    let r1 = std::fs::read_to_string(format!("{}_R1.fastq", dir.join("counted").display())).unwrap();
    let quals: Vec<u8> = r1
        .lines()
        .skip(3)
        .step_by(4)
        .flat_map(|line| line.bytes())
        .collect();
    let mean = quals.iter().map(|&b| (b - 33) as f64).sum::<f64>() / quals.len() as f64;
    assert!(
        (mean - mq).abs() < 0.05,
        "reported mean_q {mq} vs emitted {mean} (R1 only, so a close match)"
    );

    std::fs::remove_dir_all(&dir).unwrap();
}

// ---------------------------------------------------------------------------
// GC authority, --insert-min, thread resolution, and the help text
// ---------------------------------------------------------------------------

/// The clap help, whitespace-flattened so a phrase that happens to straddle a
/// wrap column still matches.
fn help_text() -> String {
    let out = Command::new(bin()).arg("--help").output().unwrap();
    assert!(out.status.success());
    String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// Emitted GC of a plain FASTQ's sequence lines.
fn emitted_gc(path: &Path) -> f64 {
    let s = std::fs::read_to_string(path).unwrap();
    let (mut gc, mut total) = (0usize, 0usize);
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

/// Profile carrying a per-position base plane — the shape every fitted profile
/// has (fastcheck `base_position`). `read_len` 60 and an insert of exactly
/// 120 bp means the two 60-position planes meet with no interior gap, which is
/// the regime a 150 bp read against a 300 bp insert is in.
fn plane_profile(dir: &Path, name: &str, gc: f64) -> std::path::PathBuf {
    let read_len = 60usize;
    let q = vec!["30"; read_len].join(",");
    // A uniform plane: 25% of every base at every position, so the emitted GC is
    // 0.5 no matter what --gc asks for.
    let row = vec![0.25f32, 0.25, 0.25, 0.25];
    let json = format!(
        r#"{{
        "read_len": {l}, "gc": {gc}, "geometry": "tail3p_both",
        "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "joint_tail": {{"0,0": 100}},
        "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.0,
        "full_length_q": {{"R1": 20, "R2": 20}},
        "insert": {{"120": 100}},
        "base": {plane},
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
    }}"#,
        l = read_len,
        gc = gc,
        q = q,
        plane = serde_json::json!({
            "R1": vec![row.clone(); read_len],
            "R2": vec![row; read_len],
        })
    );
    let p = dir.join(format!("{name}.json"));
    std::fs::write(&p, json).unwrap();
    p
}

/// A profile with a measured insert CDF (and no base plane), `read_len` 150.
fn insert_cdf_profile(dir: &Path, name: &str) -> std::path::PathBuf {
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
        "insert": {{"150": 100}},
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"}}
    }}"#,
        q = q
    );
    let p = dir.join(format!("{name}.json"));
    std::fs::write(&p, json).unwrap();
    p
}

/// A profile that carries BOTH measured authorities the CLI flags cannot
/// override: a `dup_rate` (0.0) and an insert CDF (150 bp only).
fn measured_profile(dir: &Path, name: &str) -> std::path::PathBuf {
    let q: Vec<f64> = vec![30.0; 150];
    let v = serde_json::json!({
        "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
        "dup_rate": 0.0,
        "head": {"R1": {"0": 100}, "R2": {"0": 100}},
        "tail": {"R1": {"0": 100}, "R2": {"0": 100}},
        "joint_tail": {"0,0": 100},
        "tail_q_collapse": {"R1": 0.0, "R2": 0.0},
        "pos_q": {"R1": q.clone(), "R2": q},
        "read_meanq": {"R1": {"30": 100}, "R2": {"30": 100}},
        "full_length_rate": 0.0,
        "full_length_q": {"R1": 20, "R2": 20},
        "insert": {"150": 100},
        "adapter": {"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"},
    });
    let p = dir.join(format!("{name}.json"));
    std::fs::write(&p, serde_json::to_string_pretty(&v).unwrap()).unwrap();
    p
}

/// S2: a profile carrying a per-position base plane governs the composition, so
/// `--gc` cannot change a single emitted base. The run must say so instead of
/// silently publishing the requested value, and the report must publish the
/// ACHIEVED gc (the profile's), never the requested one.
#[test]
fn gc_flag_is_inert_when_the_profile_carries_a_base_plane() {
    let dir = std::env::temp_dir().join(format!("fastsim_plane_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let profile = plane_profile(&dir, "plane", 0.31);

    let run = |tag: &str, gc: &str| -> (Vec<u8>, serde_json::Value, String) {
        let out = dir.join(tag);
        let rep = dir.join(format!("{tag}.json"));
        let o = Command::new(bin())
            .args([
                "--pairs", "300", "--profile", profile.to_str().unwrap(),
                "--out", out.to_str().unwrap(), "--json", rep.to_str().unwrap(),
                "--no-gzip", "--force", "--drop-rate", "0", "--gc", gc,
            ])
            .output()
            .unwrap();
        assert!(o.status.success(), "run {tag}: {}", String::from_utf8_lossy(&o.stderr));
        let r1 = std::fs::read(format!("{}_R1.fastq", out.display())).unwrap();
        let report = serde_json::from_slice(&std::fs::read(&rep).unwrap()).unwrap();
        (r1, report, String::from_utf8_lossy(&o.stderr).into_owned())
    };

    let (r1_low, low, err_low) = run("low", "0.20");
    let (r1_high, high, _) = run("high", "0.90");
    assert_eq!(
        r1_low, r1_high,
        "--gc moved the emitted bases although the profile carries a base plane"
    );

    // The achieved composition is the plane's 0.5, not 0.20 or 0.90.
    let observed = emitted_gc(&std::path::PathBuf::from(format!(
        "{}_R1.fastq",
        dir.join("low").display()
    )));
    assert!(
        (observed - 0.5).abs() < 0.01,
        "the plane must govern the composition: emitted GC {observed}"
    );

    // The report publishes what the run achieved, never what was requested.
    for (label, rep) in [("low", &low), ("high", &high)] {
        let gc = rep["metrics"]["gc"].as_f64().unwrap();
        assert!(
            (gc - 0.31).abs() < 1e-9,
            "{label}: report publishes {gc}, must be the profile's 0.31, not the flag"
        );
    }

    // And the user is told the flag did nothing.
    assert!(
        err_low.contains("--gc") && err_low.contains("no effect"),
        "a no-effect --gc must warn: {err_low}"
    );
}

/// S4: `--insert-min` is documented as ignored when the profile carries a
/// measured insert CDF, so an explicit value above the read length must warn and
/// be dropped — not refuse a run whose insert distribution the flag cannot
/// influence. Without a CDF the value is live, so the refusal stays.
#[test]
fn insert_min_above_read_len_is_ignored_only_when_an_insert_cdf_exists() {
    let dir = std::env::temp_dir().join(format!("fastsim_imin_cdf_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let cdf = insert_cdf_profile(&dir, "cdf");

    let o = Command::new(bin())
        .args([
            "--pairs", "50", "--profile", cdf.to_str().unwrap(),
            "--out", dir.join("c").to_str().unwrap(), "--no-gzip",
            "--insert-min", "200", "--drop-rate", "0",
        ])
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        o.status.success(),
        "--insert-min is inert with an insert CDF, so it must not refuse the run: {err}"
    );
    assert!(
        err.contains("--insert-min") && err.contains("ignored"),
        "the ignored flag must be named in a warning: {err}"
    );
    // The CDF really governed: every insert is the measured 150 bp.
    let gt = std::fs::read_to_string(dir.join("c.gt.tsv")).unwrap();
    assert!(
        gt.lines().skip(1).all(|l| l.split('\t').nth(1) == Some("150")),
        "the measured insert CDF must set the insert length"
    );

    // An explicit --read-len matching the profile must not resurrect the
    // refusal from a second, parse-time copy of the rule.
    let o = Command::new(bin())
        .args([
            "--pairs", "10", "--profile", cdf.to_str().unwrap(),
            "--out", dir.join("c2").to_str().unwrap(), "--no-gzip",
            "--read-len", "150", "--insert-min", "200", "--drop-rate", "0",
        ])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "--read-len must not turn the inert --insert-min into a refusal: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    // No CDF: the floor is live, so a floor above read_len is still a refusal.
    let plain = write_test_profile(&dir);
    let o = Command::new(bin())
        .args([
            "--pairs", "10", "--profile", plain.to_str().unwrap(),
            "--out", dir.join("p").to_str().unwrap(), "--insert-min", "200",
        ])
        .output()
        .unwrap();
    assert!(!o.status.success(), "a live insert floor above read_len must be refused");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--insert-min (200)"), "stderr: {err}");

    std::fs::remove_dir_all(&dir).unwrap();
}

/// S3: the flags that lose to the profile must say so, like `--gc` does.
///
/// `--dup-rate` loses to a measured `dup_rate` and
/// `--insert-mean`/`--insert-sd`/`--insert-min` lose to a measured insert CDF;
/// all four used to be dropped in silence (only `--insert-min` above the read
/// length warned), so a harness pinning them for reproducibility published a
/// run the profile had shaped instead. Precedence is unchanged — the profile
/// still wins — and the warning goes to stderr so `-q` keeps stdout clean.
#[test]
fn profile_authority_flags_warn_when_inert() {
    let dir = std::env::temp_dir().join(format!("fastsim_inert_flags_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let profile = measured_profile(&dir, "measured");
    let out = dir.join("w");
    let base = [
        "--pairs",
        "200",
        "--profile",
        profile.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--no-gzip",
        "--drop-rate",
        "0",
        "--dup-rate",
        "0.9",
        "--insert-mean",
        "30",
        "--insert-sd",
        "1",
        "--insert-min",
        "20",
    ];

    let o = Command::new(bin()).args(base).output().unwrap();
    assert!(
        o.status.success(),
        "an inert flag must not refuse the run: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(o.stdout.is_empty(), "warnings must not go to stdout");
    for flag in ["--dup-rate", "--insert-mean", "--insert-sd"] {
        assert!(
            err.contains(flag) && err.contains("no effect"),
            "{flag} loses to the profile and must warn: {err}"
        );
    }
    assert!(
        err.contains("--insert-min"),
        "--insert-min loses to the insert CDF too: {err}"
    );

    // The profile really did win: its 150 bp insert and its 0.0 dup rate.
    let gt = std::fs::read_to_string(dir.join("w.gt.tsv")).unwrap();
    assert!(
        gt.lines().skip(1).all(|l| l.split('\t').nth(1) == Some("150")),
        "--insert-mean 30 must not move a single insert"
    );
    assert!(
        gt.lines().skip(1).all(|l| l.split('\t').nth(9) == Some("0")),
        "--dup-rate 0.9 must not create a duplicate the profile's 0.0 forbids"
    );

    // -q silences the summary, never the warning, and writes nothing to stdout.
    let quiet_out = dir.join("q");
    let o = Command::new(bin())
        .args([
            "--pairs",
            "50",
            "--profile",
            profile.to_str().unwrap(),
            "--out",
            quiet_out.to_str().unwrap(),
            "--no-gzip",
            "--drop-rate",
            "0",
            "--dup-rate",
            "0.9",
            "-q",
        ])
        .output()
        .unwrap();
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(o.stdout.is_empty(), "stdout must stay clean under -q");
    assert!(
        err.contains("--dup-rate") && err.contains("no effect"),
        "the warning must survive -q: {err}"
    );
    assert!(
        !err.contains("fastsim: 50 pairs"),
        "the summary line must not leak under -q: {err}"
    );

    std::fs::remove_dir_all(&dir).unwrap();
}

/// S5: `FASTSUITE_THREADS` is the default for an ABSENT thread flag, in every
/// spelling clap accepts — `-j 4`, `-j4`, `-j=4`, `--threads 4`, `--threads=4`.
/// The detached spellings used to be the only ones it recognised, so an explicit
/// `-j4` was silently overridden by the environment.
#[test]
fn env_threads_never_overrides_an_explicit_flag() {
    let dir = std::env::temp_dir().join(format!("fastsim_env_spell_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    for (i, extra) in [
        vec!["-j", "1"],
        vec!["-j1"],
        vec!["-j=1"],
        vec!["--threads", "1"],
        vec!["--threads=1"],
    ]
    .iter()
    .enumerate()
    {
        let out = dir.join(format!("s{i}"));
        let o = Command::new(bin())
            .args([
                "--pairs", "10", "--profile", profile.to_str().unwrap(),
                "--out", out.to_str().unwrap(), "-q",
            ])
            .args(extra)
            .env("FASTSUITE_THREADS", "4")
            .output()
            .unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(
            !err.contains("--threads"),
            "{extra:?} was overridden by FASTSUITE_THREADS=4: {err}"
        );
    }

    // An explicit -j 4 is the value used (warning says 4), not the env's 2.
    let o = Command::new(bin())
        .args([
            "--pairs", "10", "--profile", profile.to_str().unwrap(),
            "--out", dir.join("four").to_str().unwrap(), "-j4",
        ])
        .env("FASTSUITE_THREADS", "2")
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--threads 4"), "explicit -j4 must win: {err}");

    std::fs::remove_dir_all(&dir).unwrap();
}

/// Cross-cutting: the resolved thread count (`-j` or the env default) is
/// published in the report, so a harness can record what was asked for.
#[test]
fn resolved_thread_count_reaches_the_report() {
    let dir = std::env::temp_dir().join(format!("fastsim_threads_rep_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let profile = write_test_profile(&dir);

    let run = |tag: &str, extra: &[&str], env: Option<&str>| -> serde_json::Value {
        let out = dir.join(tag);
        let rep = dir.join(format!("{tag}.json"));
        let mut cmd = Command::new(bin());
        cmd.args([
            "--pairs", "10", "--profile", profile.to_str().unwrap(),
            "--out", out.to_str().unwrap(), "--json", rep.to_str().unwrap(),
            "--force", "-q",
        ])
        .args(extra)
        .env_remove("FASTSUITE_THREADS");
        if let Some(e) = env {
            cmd.env("FASTSUITE_THREADS", e);
        }
        let o = cmd.output().unwrap();
        assert!(o.status.success(), "{tag}: {}", String::from_utf8_lossy(&o.stderr));
        serde_json::from_slice(&std::fs::read(&rep).unwrap()).unwrap()
    };

    assert_eq!(run("flag", &["-j", "3"], Some("8"))["metrics"]["threads"], 3);
    assert_eq!(run("env", &[], Some("5"))["metrics"]["threads"], 5);
    assert_eq!(run("default", &[], None)["metrics"]["threads"], 1);

    std::fs::remove_dir_all(&dir).unwrap();
}

/// S7: the help text must state the truth. There is no "censored" column (a
/// length above the measurable ceiling lands in ONE saturated histogram bin),
/// and `gc` is a required profile key, not a 0.49 fallback. The `--gc` help also
/// names the base-plane precedence, since that is what makes the flag inert.
#[test]
fn help_text_states_the_measured_truth() {
    let help = help_text();
    assert!(
        !help.contains("censored"),
        "no censored column exists: {help}"
    );
    assert!(
        help.contains("saturated"),
        "the insert ceiling must be described as a saturated bin: {help}"
    );
    assert!(
        !help.contains("0.49"),
        "gc is a required profile key, not a 0.49 fallback: {help}"
    );
    assert!(
        help.contains("base plane"),
        "the --gc help must say a base plane wins: {help}"
    );
}
