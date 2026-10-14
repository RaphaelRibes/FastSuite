//! fastsim's CLI flag contract, in one file: the FastSuite common flags
//! (`--force`, `-q/--quiet`), the `FASTSUITE_THREADS` env default, and the
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
}
