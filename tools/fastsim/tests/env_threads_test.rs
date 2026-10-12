//! FASTSUITE_THREADS env default (CLI-SPEC §1, gap #6).

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
    let p = dir.join("env_profile.json");
    std::fs::write(&p, json).unwrap();
    p
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
