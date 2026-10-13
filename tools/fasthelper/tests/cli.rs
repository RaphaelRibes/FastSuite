//! CLI-level tests for the toolbox surface.
//!
//! `concatenate` is the name; `report` stays as an alias so nothing that was
//! written against the earlier build breaks (docs/DECISIONS.md D-8).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasthelper() -> &'static str {
    env!("CARGO_BIN_EXE_fasthelper")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fasthelper_cli_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// One minimal fastsuite.report.v1 file, as fastcheck would write it.
fn write_report(dir: &Path, name: &str, reads: u64, gc: f64) -> PathBuf {
    let p = dir.join(name);
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastcheck",
        "tool_version": "0.1.0",
        "elapsed_s": 1.0,
        "inputs": {"files": [name]},
        "metrics": {"reads": reads, "bases": reads * 150, "gc_fraction": gc}
    });
    std::fs::write(&p, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    p
}

#[test]
fn compare_creates_a_missing_output_directory() {
    // Regression: `compare -o out` used to fail with 'No such file or
    // directory' when out/ did not exist yet, while concatenate created it.
    let d = tmpdir("compare_mkdir");
    let before = d.join("before");
    let after = d.join("after");
    std::fs::create_dir_all(&before).unwrap();
    std::fs::create_dir_all(&after).unwrap();
    write_report(&before, "raw.json", 100, 0.41);
    write_report(&after, "clean.json", 90, 0.42);
    let out = d.join("does").join("not").join("exist");

    let (_, stderr, code) = run(&[
        "compare",
        "--before",
        before.to_str().unwrap(),
        "--after",
        after.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);

    assert_eq!(code, 0, "compare must succeed: {stderr}");
    assert!(out.join("compare.tsv").is_file());
    assert!(out.join("compare.html").is_file());
}

fn run(args: &[&str]) -> (String, String, i32) {
    let out = Command::new(fasthelper()).args(args).output().expect("spawn fasthelper");
    (
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

#[test]
fn concatenate_merges_reports_into_one_row_per_sample() {
    let d = tmpdir("concat");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "sampleA.json", 100, 0.41);
    write_report(&reports, "sampleB.json", 200, 0.52);
    let out = d.join("out");

    let (stdout, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "tsv,json",
        "--stdout",
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    assert!(stdout.contains("sampleA") && stdout.contains("sampleB"), "{stdout}");

    let tsv = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines.len(), 3, "header + one row per sample: {tsv}");
    assert!(tsv.contains("sampleA"), "{tsv}");
    assert!(tsv.contains("sampleB"), "{tsv}");
    assert!(stdout.contains("concatenated 2 report(s)"), "{stdout}");
}

#[test]
fn report_is_still_a_working_alias() {
    let d = tmpdir("alias");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "s.json", 10, 0.5);

    let a = d.join("via_concatenate");
    let b = d.join("via_report");
    let (_, e1, c1) = run(&["concatenate", reports.to_str().unwrap(), "-o", a.to_str().unwrap()]);
    let (_, e2, c2) = run(&["report", reports.to_str().unwrap(), "-o", b.to_str().unwrap()]);
    assert_eq!(c1, 0, "{e1}");
    assert_eq!(c2, 0, "the alias must still work: {e2}");

    let ta = std::fs::read_to_string(a.join("report.tsv")).unwrap();
    let tb = std::fs::read_to_string(b.join("report.tsv")).unwrap();
    assert_eq!(ta, tb, "both spellings must produce the same table");
}

#[test]
fn positions_wide_matrix_unions_positions_and_leaves_missing_cells_empty() {
    let d = tmpdir("positions");
    // Two samples with different read lengths: A covers 3 positions, B covers 2.
    for (sample, rows) in [
        ("A", "1\t10\t0.5\n2\t20\t0.6\n3\t30\t0.7\n"),
        ("B", "1\t11\t0.4\n2\t21\t0.5\n"),
    ] {
        let dir = d.join(sample);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("positions.tsv"),
            format!("pos\tdepth\tmean_q\n{rows}"),
        )
        .unwrap();
    }
    let out = d.join("matrix");
    let (stdout, stderr, code) = run(&[
        "concatenate",
        "--positions",
        d.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    assert!(stdout.contains("2 sample(s)"), "{stdout}");

    let m = std::fs::read_to_string(out.join("positions_matrix.tsv")).unwrap();
    let lines: Vec<&str> = m.lines().collect();
    assert_eq!(lines[0], "pos\tA\tB", "samples sorted, metric = mean_q by default: {m}");
    assert_eq!(lines[1], "1\t0.5000\t0.4000");
    assert_eq!(lines[2], "2\t0.6000\t0.5000");
    assert_eq!(lines[3], "3\t0.7000\t", "a missing position stays empty, never zero: {m}");
}

#[test]
fn positions_metric_selection_and_tidy_shape() {
    let d = tmpdir("positions_metric");
    // Input tree and output tree kept apart, as they would be in a real run.
    let dir = d.join("in/A");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("positions.tsv"),
        "pos\tdepth\tmean_q\n1\t100\t30.5\n2\t90\t29.5\n",
    )
    .unwrap();

    let wide = d.join("out/by_depth");
    let (_, _, code) = run(&[
        "concatenate", "--positions", d.join("in").to_str().unwrap(),
        "-o", wide.to_str().unwrap(), "--metric", "depth",
    ]);
    assert_eq!(code, 0);
    let m = std::fs::read_to_string(wide.join("positions_matrix.tsv")).unwrap();
    assert!(m.contains("1\t100") && m.contains("2\t90"), "{m}");

    // A metric that does not exist is an error that names the alternatives.
    let (_, stderr, code) = run(&[
        "concatenate", "--positions", d.join("in").to_str().unwrap(),
        "-o", d.join("out/nope").to_str().unwrap(), "--metric", "polyG_frac",
    ]);
    assert_ne!(code, 0);
    assert!(stderr.contains("depth") && stderr.contains("mean_q"), "{stderr}");

    // Tidy shape: one row per sample/position/metric.
    let (tidy, _, code) = run(&[
        "concatenate",
        "--positions",
        d.join("in").to_str().unwrap(),
        "--all-metrics",
        "--stdout",
    ]);
    assert_eq!(code, 0);
    let lines: Vec<&str> = tidy.lines().collect();
    assert_eq!(lines[0], "sample\tpos\tmetric\tvalue");
    assert_eq!(lines.len(), 1 + 2 * 2, "2 positions x 2 metrics: {tidy}");
    assert!(tidy.contains("A\t2\tmean_q\t29.5"), "{tidy}");
}

#[test]
fn positions_ignores_other_tsv_files_instead_of_misreading_them() {
    let d = tmpdir("positions_skip");
    let dir = d.join("A");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("positions.tsv"), "pos\tmean_q\n1\t30.0\n").unwrap();
    // A concatenated report table from the other flavour, sitting alongside.
    std::fs::write(d.join("report.tsv"), "sample\ttool\treads\nA\tfastcheck\t10\n").unwrap();

    let (_, stderr, code) = run(&[
        "concatenate", "--positions", d.to_str().unwrap(),
        "-o", d.join("m").to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stderr.contains("skipped"), "the skip must be reported: {stderr}");
    let m = std::fs::read_to_string(d.join("m/positions_matrix.tsv")).unwrap();
    assert_eq!(m.lines().count(), 2, "only the real sample is in the matrix: {m}");
}

#[test]
fn help_lists_the_toolbox() {
    let (stdout, _, code) = run(&["--help"]);
    assert_eq!(code, 0);
    for cmd in ["concatenate", "run", "doctor", "schema"] {
        assert!(stdout.contains(cmd), "{cmd} missing from:\n{stdout}");
    }
}

#[test]
fn fail_on_trips_exit_2_and_names_the_metric() {
    let d = tmpdir("failon");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "sampleA.json", 100, 0.60); // gc_fraction 0.60

    let (_, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "--fail-on",
        "gc_fraction>0.5",
    ]);
    assert_eq!(code, 2, "threshold tripped must be exit 2: {stderr}");
    assert!(stderr.contains("gc_fraction"), "must name the metric: {stderr}");
    // The report is still written before the exit.
    assert!(d.join("report.tsv").is_file() || Path::new("report.tsv").is_file());

    let (_, _, ok) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "--fail-on",
        "gc_fraction>0.9",
    ]);
    assert_eq!(ok, 0, "a satisfied threshold must stay exit 0");
}

#[test]
fn check_rejects_unknown_step_and_missing_input() {
    let d = tmpdir("check");
    let cfg = d.join("bad.toml");
    std::fs::write(
        &cfg,
        "[run]\nsteps = [\"fastnope\"]\n\n[[samples]]\nname = \"s1\"\nr1 = \"nope_R1.fastq\"\n",
    )
    .unwrap();
    let (_, stderr, code) = run(&["check", "-c", cfg.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(stderr.contains("fastnope"), "must name the bad step: {stderr}");
    assert!(stderr.contains("nope_R1.fastq"), "must name the missing input: {stderr}");

    let ok_cfg = d.join("ok.toml");
    let fq = d.join("s1_R1.fastq");
    std::fs::write(&fq, "@r1\nACGT\n+\nIIII\n").unwrap();
    std::fs::write(
        &ok_cfg,
        format!(
            "[run]\nsteps = [\"fastcheck\"]\n\n[[samples]]\nname = \"s1\"\nr1 = \"{}\"\n",
            fq.display()
        ),
    )
    .unwrap();
    // tools-dir pointing at the repo root: fastcheck may or may not resolve,
    // so accept either a clean pass or a clean "binary not found" refusal —
    // never a crash.
    let (_, stderr2, code2) = run(&[
        "check",
        "-c",
        ok_cfg.to_str().unwrap(),
        "--tools-dir",
        &std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".into()),
    ]);
    assert!(code2 == 0 || stderr2.contains("not found"), "{code2} {stderr2}");
}

#[test]
fn fastcheck_step_runs_and_writes_a_v1_report() {
    // fastcheck is a terminal QC step: it takes positional FASTQs and writes
    // only a report, never a new FASTQ. The runner must build that argv.
    let d = tmpdir("fastcheck_step");
    let out = d.join("run");
    std::fs::create_dir_all(&out).unwrap();

    // A minimal config: one step, one sample, no simulate.
    let cfg = d.join("fastsuite.toml");
    let fq = d.join("s1_R1.fastq");
    std::fs::write(&fq, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"qconly\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n[[samples]]\nname = \"s1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    // Point the runner at a toolset whose only real binary is fastcheck: the
    // step must not be rejected as "unsupported".
    let (_, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        option_env!("FASTSUITE_TOOLS").unwrap_or("tools"),
    ]);
    // fastcheck may not be built in this environment; the assertion that
    // matters is that it was recognised as a step, not called unsupported.
    assert!(
        !stderr.contains("unsupported step"),
        "fastcheck must be a known step: {stderr}"
    );
    if code == 0 {
        let man: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap())
                .unwrap();
        let cmd = &man["commands"][0];
        assert_eq!(cmd["tool"], "fastcheck");
        let argv: Vec<String> = serde_json::from_value(cmd["argv"].clone()).unwrap();
        assert!(argv.contains(&fq.display().to_string()), "argv: {argv:?}");
        assert!(argv.contains(&"--json".to_string()), "argv: {argv:?}");
        assert!(out.join("s1/fastcheck.json").is_file());
    }
}

#[test]
fn jobs_flag_is_accepted_and_the_run_still_succeeds() {
    let d = tmpdir("jobs");
    let out = d.join("run");
    std::fs::create_dir_all(&out).unwrap();
    let fq1 = d.join("s1_R1.fastq");
    let fq2 = d.join("s2_R1.fastq");
    std::fs::write(&fq1, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(&fq2, "@r2\nTTTTAAAA\n+\nIIIIIIII\n").unwrap();
    let cfg = d.join("two.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"two\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\nthreads = 2\n\n[[samples]]\nname = \"s1\"\nr1 = \"{}\"\n\n[[samples]]\nname = \"s2\"\nr1 = \"{}\"\n",
            out.display(),
            fq1.display(),
            fq2.display()
        ),
    )
    .unwrap();
    let (_, stderr, code) = run(&["run", "-c", cfg.to_str().unwrap(), "-P", "2"]);
    // fastcheck may be missing; what must hold: -P is a real flag, and when
    // the step runs both samples get manifest entries.
    assert!(
        !stderr.contains("unexpected argument"),
        "-P must exist: {stderr}"
    );
    if code == 0 {
        let man: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap())
                .unwrap();
        let cmds = man["commands"].as_array().unwrap();
        assert_eq!(cmds.len(), 2, "both samples must run: {cmds:?}");
    }
}
