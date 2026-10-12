//! `-q/--quiet` (CLI-SPEC §1, rule 5): the progress and end-of-run summary
//! stream is suppressed; warnings and errors are NOT. The suite rule is that
//! `--quiet` suppresses the summary line only, and the fastsuite skill's
//! pitfall is the reason it matters here — a bench harness that appends `-q`
//! used to get `error: unexpected argument` from fastprep (the flag did not
//! exist), so every row came back empty instead of measured.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_quiet_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_fastq(dir: &Path) -> PathBuf {
    let p = dir.join("in.fastq");
    // Random-ish ACGT: no built-in adapter fragment appears, so adapter
    // auto-detection legitimately finds nothing and emits a warning.
    let mut body = String::new();
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    for i in 0..200 {
        let mut seq = String::with_capacity(60);
        for _ in 0..60 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            seq.push(match state % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            });
        }
        body.push_str(&format!("@read{i}\n{seq}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&p, body).unwrap();
    p
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

#[test]
fn quiet_suppresses_the_summary_stream_only() {
    let d = tmpdir("summary");
    let fq = write_fastq(&d);
    let out = d.join("out.fastq");

    // The loud run: everything the summary stream is made of.
    let (loud, code) = run(&["-1", fq.to_str().unwrap(), "-o", out.to_str().unwrap()]);
    assert_eq!(code, 0, "the loud run must succeed: {loud}");
    for expected in [
        "fastprep:",
        "[R1]",
        "reads in:",
        "kept:",
        "dropped:",
        "elapsed:",
        "please cite",
    ] {
        assert!(
            loud.contains(expected),
            "the loud run must print {expected:?}: {loud}"
        );
    }

    // -q: nothing of that stream survives.
    let (quiet, code) = {
        let mut a = vec!["-1", fq.to_str().unwrap(), "-o", out.to_str().unwrap()];
        a.extend(["--min-len", "0"]);
        a.push("-q");
        a.push("--force");
        run(&a)
    };
    assert_eq!(code, 0, "-q must not change the outcome: {quiet}");
    for gone in [
        "single-end complete",
        "[R1]",
        "reads in:",
        "kept:",
        "dropped:",
        "elapsed:",
        "please cite",
    ] {
        assert!(
            !quiet.contains(gone),
            "-q must suppress {gone:?}: {quiet}"
        );
    }
    // The adapter assumption is a warning, not summary: it stays, and it is
    // the only thing left. (This is exactly the D-6 behaviour -q must not
    // be able to hide.)
    assert!(
        quiet.contains("assumption"),
        "-q must keep the D-6 assumption warning: {quiet}"
    );
    assert_eq!(
        quiet.lines().count(),
        1,
        "-q must leave only the warning line, got: {quiet}"
    );

    // The output is still written: -q gates the stream, not the work.
    assert!(out.exists(), "-q must still write the output");

    // --quiet is the same flag. With --no-adapter there is no assumption to
    // warn about, so stderr really is empty.
    let (long, code) = {
        let mut a = vec!["-1", fq.to_str().unwrap(), "-o", out.to_str().unwrap()];
        a.extend(["--min-len", "0", "--quiet", "--force", "--no-adapter"]);
        run(&a)
    };
    assert_eq!(code, 0, "--quiet must succeed: {long}");
    assert!(
        long.trim().is_empty(),
        "--quiet must behave like -q: {long}"
    );

    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn quiet_keeps_warnings_and_errors() {
    let d = tmpdir("warnings");
    let fq = write_fastq(&d);
    let out = d.join("out.fastq");

    // Adapter auto-detection finds nothing here and labels the assumption.
    // That is a warning, so -q must not hide it.
    let (warn_run, code) = {
        let mut a = vec!["-1", fq.to_str().unwrap(), "-o", out.to_str().unwrap()];
        a.extend(["--min-len", "0", "-q"]);
        run(&a)
    };
    assert_eq!(code, 0, "the assumption path must still exit 0: {warn_run}");
    assert!(
        warn_run.contains("assumption"),
        "-q must keep warnings, got: {warn_run}"
    );

    // A genuine error survives -q too: --force still refuses nowhere, so use
    // a bad range, which is a usage error on the summary-less path.
    let (err_run, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-len",
        "0",
        "-q",
        "--phred-sat",
        "0",
    ]);
    assert_ne!(code, 0, "an invalid --phred-sat must still fail");
    assert!(
        err_run.contains("--phred-sat"),
        "-q must keep errors, got: {err_run}"
    );

    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn quiet_and_dry_run_together_still_report() {
    // --dry-run's contract is that the report is still emitted (it is the
    // whole point of a dry run), so it composes with -q the same way a real
    // run's summary does: suppressed by -q, present without it.
    let d = tmpdir("dryrun");
    let fq = write_fastq(&d);
    let out = d.join("dry.fastq");

    let (plain, code) = {
        let mut a = vec!["-1", fq.to_str().unwrap(), "-o", out.to_str().unwrap()];
        a.extend(["--min-len", "0", "--dry-run"]);
        run(&a)
    };
    assert_eq!(code, 0, "the dry run must succeed: {plain}");
    assert!(plain.contains("reads in:"), "dry run must report: {plain}");

    let (quiet, code) = {
        let mut a = vec!["-1", fq.to_str().unwrap(), "-o", out.to_str().unwrap()];
        a.extend(["--min-len", "0", "--dry-run", "-q", "--no-adapter"]);
        run(&a)
    };
    assert_eq!(code, 0, "--dry-run -q must succeed: {quiet}");
    assert!(
        quiet.trim().is_empty(),
        "--dry-run -q must suppress the summary too: {quiet}"
    );
    assert!(!out.exists(), "--dry-run must still write nothing");

    std::fs::remove_dir_all(&d).unwrap();
}
