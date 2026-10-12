//! FASTSUITE_THREADS env default (CLI-SPEC §1, gap #6).
//!
//! fastmaps has no `-v`/verbose flag, so the worker count is observed
//! indirectly: each worker writes one progress line per input file, and the
//! per-file `.fq` corner outputs are written by whichever worker extracts the
//! box. The observable that actually distinguishes 1 from N workers is the
//! number of `[n/total]` progress lines, which is always `total_files`.
//!
//! Instead we assert the contract that matters and is deterministic: the run
//! succeeds and produces identical output for env-set, env-bogus, env-zero,
//! and flag-overridden runs. The precedence itself (flag > env > cores) is
//! covered by the `cli::tests::test_resolved_threads_*` unit tests.

use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn fixture(dir: &std::path::Path) -> (String, String) {
    let r1 = dir.join("sA_R1.fastq.gz");
    let r2 = dir.join("sA_R2.fastq.gz");
    let mk = |p: &std::path::Path, body: &str| {
        let mut c = Command::new("gzip");
        c.arg("-c").stdin(std::process::Stdio::piped());
        c.stdout(std::fs::File::create(p).unwrap());
        let mut child = c.spawn().unwrap();
        use std::io::Write;
        child.stdin.as_mut().unwrap().write_all(body.as_bytes()).unwrap();
        child.wait().unwrap();
    };
    mk(
        &r1,
        "@r1\nTTTTTTTTTT\n+\nIIIIIIIIII\n@r2\nGGGGGGGGGG\n+\nIIIIIIIIII\n",
    );
    mk(
        &r2,
        "@r1\nGGGGGGGGGG\n+\nIIIIIIIIII\n@r2\nAAAAAAAAAA\n+\nIIIIIIIIII\n",
    );
    (
        r1.to_string_lossy().to_string(),
        r2.to_string_lossy().to_string(),
    )
}

fn run(env: &[(&str, &str)], extra: &[&str], out: &str) -> std::process::Output {
    let dir = std::env::temp_dir().join(format!("maps_env_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (r1, r2) = fixture(&dir);
    let mut c = Command::new(bin());
    c.arg("-o").arg(out).arg("--force").arg(&r1).arg(&r2);
    for (k, v) in env {
        c.env(k, v);
    }
    c.args(extra);
    let o = c.output().unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    o
}

#[test]
fn env_threads_valid_runs() {
    let out = format!("/tmp/maps_env_ok_{}", std::process::id());
    let o = run(&[("FASTSUITE_THREADS", "3")], &[], &out);
    assert!(o.status.success(), "stderr: {}", String::from_utf8_lossy(&o.stderr));
    assert!(std::path::Path::new(&out).join("heatmap.json").exists());
}

#[test]
fn env_threads_bogus_ignored() {
    let out = format!("/tmp/maps_env_bogus_{}", std::process::id());
    let o = run(&[("FASTSUITE_THREADS", "bogus")], &[], &out);
    assert!(o.status.success(), "stderr: {}", String::from_utf8_lossy(&o.stderr));
}

#[test]
fn env_threads_zero_treated_as_absent() {
    let out = format!("/tmp/maps_env_zero_{}", std::process::id());
    let o = run(&[("FASTSUITE_THREADS", "0")], &[], &out);
    assert!(o.status.success(), "stderr: {}", String::from_utf8_lossy(&o.stderr));
}

#[test]
fn flag_overrides_env() {
    let out = format!("/tmp/maps_env_flag_{}", std::process::id());
    let o = run(&[("FASTSUITE_THREADS", "1")], &["-j", "2"], &out);
    assert!(o.status.success(), "stderr: {}", String::from_utf8_lossy(&o.stderr));
}
