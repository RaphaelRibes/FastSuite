//! Parity tests for the FastSuite common flags --force and -q/--quiet.
//! (Merged into cli_flags_test.rs by Task 7; kept separate here so the
//! hand-parser behaviour is pinned before the clap migration.)

use std::path::Path;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastsim")
}

fn write_test_profile(dir: &Path) -> std::path::PathBuf {
    let q150 = vec!["30"; 150].join(",");
    let json = format!(
        r#"{{"read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
        "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
        "joint_tail": {{"0,0": 100}},
        "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
        "pos_q": {{"R1": [{q}], "R2": [{q}]}},
        "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
        "full_length_rate": 0.0,
        "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"}}}}"#,
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
