//! CLI-level tests for the toolbox surface.
//!
//! `concatenate` is the name; `report` stays as an alias so nothing that was
//! written against the earlier build breaks (docs/DECISIONS.md D-8).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasthelper() -> &'static str {
    env!("CARGO_BIN_EXE_fasthelper")
}

/// tools/ — the directory `--tools-dir` expects (holds <tool>/target/<profile>/<tool>).
fn tools_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

/// Locate a built suite binary, or fail LOUDLY. These tests exercise the real
/// runner against the real fastcheck, so a missing binary is a broken
/// prerequisite, not a reason to pass silently: the test would otherwise go
/// green without ever proving the step ran.
fn require_tool(tool: &str) -> PathBuf {
    let dir = tools_dir().join(tool).join("target");
    for profile in ["release", "debug"] {
        let p = dir.join(profile).join(tool);
        if p.is_file() {
            return p;
        }
    }
    panic!(
        "{tool} is not built under {} — build it (cd tools/{tool} && cargo build --release) \
         before running this test; it is a real gate, not an optional one",
        tools_dir().display()
    );
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
    // Real gate: requires fastcheck to be built (see require_tool) and asserts
    // the step actually ran and wrote its report — no silent pass if it didn't.
    require_tool("fastcheck");
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

    let (_, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "the fastcheck step must run clean: {stderr}");
    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap())
            .unwrap();
    let cmd = &man["commands"][0];
    assert_eq!(cmd["tool"], "fastcheck");
    let argv: Vec<String> = serde_json::from_value(cmd["argv"].clone()).unwrap();
    assert!(argv.contains(&fq.display().to_string()), "argv: {argv:?}");
    assert!(argv.contains(&"--json".to_string()), "argv: {argv:?}");
    assert!(
        out.join("s1/fastcheck.json").is_file(),
        "fastcheck must have written its report, manifest: {man}"
    );
}

#[test]
fn jobs_flag_is_accepted_and_the_run_still_succeeds() {
    require_tool("fastcheck");
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
    let (_, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "-P",
        "2",
        "--tools-dir",
        tools_dir().to_str().unwrap(),
    ]);
    // -P is a real flag AND the run must succeed: both samples get a manifest entry.
    assert_eq!(code, 0, "-P must be accepted and the run must succeed: {stderr}");
    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap())
            .unwrap();
    let cmds = man["commands"].as_array().unwrap();
    assert_eq!(cmds.len(), 2, "both samples must run: {cmds:?}");
}

#[test]
fn from_flag_forces_a_step_rerun_and_resume_skips_matches() {
    require_tool("fastcheck");
    let d = tmpdir("resume");
    let out = d.join("run");
    std::fs::create_dir_all(&out).unwrap();
    let fq = d.join("s1_R1.fastq");
    std::fs::write(&fq, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    let cfg = d.join("one.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"r\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n[[samples]]\nname = \"s1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();
    // First run populates the manifest; it must actually succeed.
    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "first run must succeed: {e1}");
    // Second run: identical config -> the step is resumed (skipped).
    let (stdout2, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c2, 0, "{e2}");
    assert!(stdout2.contains("resumed s1/fastcheck"), "{stdout2}");
    // --from fastcheck forces it to run again.
    let (stdout3, e3, c3) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td, "--from", "fastcheck"]);
    assert_eq!(c3, 0, "{e3}");
    assert!(!stdout3.contains("resumed s1/fastcheck"), "{stdout3}");
}

#[test]
fn install_symlinks_binaries_and_legacy_names() {
    let d = tmpdir("install");
    let prefix = d.join("prefix");
    std::fs::create_dir_all(&prefix).unwrap();
    let (stdout, stderr, code) = run(&[
        "install",
        "--prefix",
        prefix.to_str().unwrap(),
        "--tools-dir",
        &std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".into()),
    ]);
    if code != 0 {
        // No binaries resolved here (empty tools dir): must still say so, not crash.
        assert!(
            stderr.contains("no FastSuite binaries") || stderr.contains("not found"),
            "{stderr}"
        );
        return;
    }
    let bin = prefix.join("bin");
    assert!(bin.is_dir(), "{stdout} {stderr}");
    // Whatever was found, the directory must contain at least one symlink.
    let n = std::fs::read_dir(&bin).unwrap().count();
    assert!(n > 0, "nothing installed");
}

/// One fastsuite.report.v1 file carrying a single named metric.
fn write_metric_report(dir: &Path, name: &str, metric: &str, value: f64) -> PathBuf {
    let p = dir.join(name);
    let mut metrics = serde_json::Map::new();
    metrics.insert(metric.to_string(), serde_json::json!(value));
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastdedup",
        "tool_version": "1.3.0",
        "metrics": metrics,
    });
    std::fs::write(&p, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    p
}

#[test]
fn fail_on_gates_every_sample_not_just_the_last() {
    // Regression (B1): the gate used to flatten all samples' metrics into one
    // "last wins" map, so a violating sample was missed whenever another
    // sample came after it. Every sample must be judged, and named.
    let d = tmpdir("failon_multi");
    let a = d.join("A");
    let b = d.join("B");
    std::fs::create_dir_all(&a).unwrap();
    std::fs::create_dir_all(&b).unwrap();
    write_metric_report(&a, "fastdedup.json", "duplicate_rate", 0.5); // violates
    write_metric_report(&b, "fastdedup.json", "duplicate_rate", 0.0); // clean

    let out = d.join("out");
    // A alone trips (the baseline the old code got right).
    let (_, stderr_a, code_a) = run(&[
        "concatenate", a.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--fail-on", "duplicate_rate>0.01",
    ]);
    assert_eq!(code_a, 2, "A alone must trip: {stderr_a}");

    // A first, then B: the old "last wins" map kept B=0.0 and missed A.
    let (_, stderr_ab, code_ab) = run(&[
        "concatenate", a.to_str().unwrap(), b.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--fail-on", "duplicate_rate>0.01",
    ]);
    assert_eq!(code_ab, 2, "A B must still trip on A: {stderr_ab}");
    assert!(stderr_ab.contains('A'), "must name the tripping sample: {stderr_ab}");

    // B first, then A: same, both orders.
    let (_, stderr_ba, code_ba) = run(&[
        "concatenate", b.to_str().unwrap(), a.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--fail-on", "duplicate_rate>0.01",
    ]);
    assert_eq!(code_ba, 2, "B A must still trip on A: {stderr_ba}");
    assert!(stderr_ba.contains('A'), "must name the tripping sample: {stderr_ba}");
}

#[test]
fn concatenate_ignores_fasthelpers_own_aggregate_report() {
    // Regression (B3): a run directory holds each sample's report plus
    // fasthelper's own aggregate report.json, which shares the
    // fastsuite.report.v1 schema. Re-concatenating the directory used to
    // ingest that aggregate as a phantom third sample ("run"/"unknown").
    let d = tmpdir("own_aggregate");
    let run_dir = d.join("run");
    std::fs::create_dir_all(run_dir.join("a")).unwrap();
    std::fs::create_dir_all(run_dir.join("b")).unwrap();
    write_metric_report(&run_dir.join("a"), "fastdedup.json", "duplicate_rate", 0.5);
    write_metric_report(&run_dir.join("b"), "fastdedup.json", "duplicate_rate", 0.0);

    // Let fasthelper write its own aggregate into the run directory.
    let (_, e, c) = run(&[
        "concatenate", run_dir.join("a").to_str().unwrap(), run_dir.join("b").to_str().unwrap(),
        "-o", run_dir.to_str().unwrap(),
    ]);
    assert_eq!(c, 0, "{e}");
    assert!(run_dir.join("report.json").is_file());

    // Re-concatenating the whole directory must see exactly the two samples.
    let (stdout, stderr, code) = run(&["concatenate", run_dir.to_str().unwrap(), "-o", d.join("out").to_str().unwrap(), "--stdout"]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        stdout.contains("from 2 sample(s)"),
        "the aggregate must not add a phantom sample: {stdout}"
    );
    assert!(
        !stdout.contains("run\tunknown"),
        "fasthelper's own report.json must be skipped: {stdout}"
    );
}

#[test]
fn check_accepts_the_fastsim_step() {
    // Regression (B4): `check`'s step list omitted fastsim, which `run` can
    // execute, so check rejected a config run accepted.
    let d = tmpdir("check_fastsim");
    let fq = d.join("s1_R1.fastq");
    std::fs::write(&fq, "@r1\nACGT\n+\nIIII\n").unwrap();
    let cfg = d.join("fastsim.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"s\"\nsteps = [\"fastsim\"]\n\n[[samples]]\nname = \"s1\"\nr1 = \"{}\"\n",
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "check", "-c", cfg.to_str().unwrap(),
        "--tools-dir", &std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".into()),
    ]);
    assert!(
        !stderr.contains("unknown step"),
        "fastsim must be a known step to check: {stderr}"
    );
    // Either it validates cleanly or it reports the (possibly absent) binary —
    // never "unknown step".
    if code != 0 {
        assert!(stderr.contains("not found"), "unexpected failure: {stderr}");
    }
}

/// A stand-in fastdedup: echoes its `--version`, copies the md5 of its -1
/// input into its -o output and writes a v1 report. A stale (un-rerun) output
/// is therefore visible as an unchanged md5.
const FAKE_DEDUP: &str = r#"#!/usr/bin/env bash
if [ "$1" = "--version" ]; then echo "fastdedup 9.9.9"; exit 0; fi
in=""; out=""; json=""
while [ $# -gt 0 ]; do
  case "$1" in
    -1) in="$2"; shift 2 ;;
    -o) out="$2"; shift 2 ;;
    --json) json="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$(dirname "$out")" "$(dirname "$json")"
md5sum "$in" > "$out"
printf '{"schema":"fastsuite.report.v1","tool":"fastdedup","tool_version":"9.9.9","command":"fastdedup stub","metrics":{"duplicate_rate":0.0}}\n' > "$json"
"#;

#[test]
fn resume_reruns_when_an_input_changes_in_place() {
    // Regression (B2): resume keyed only on (argv, binary, live version) and
    // never looked at the inputs, so overwriting a sample's FASTQ in place let
    // the next run "resume" and keep a report computed from the old reads.
    let d = tmpdir("resume_stale");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    let bin = tools.join("fastdedup");
    std::fs::write(&bin, FAKE_DEDUP).unwrap();
    let mut perm = std::fs::metadata(&bin).unwrap().permissions();
    std::os::unix::fs::PermissionsExt::set_mode(&mut perm, 0o755);
    std::fs::set_permissions(&bin, perm).unwrap();

    let out = d.join("run");
    let fq = d.join("a.fq");
    std::fs::write(&fq, "@r0\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"v\"\noutdir = \"{}\"\nsteps = [\"fastdedup\"]\n\n[[samples]]\nname = \"a\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap()]);
    assert_eq!(c1, 0, "first run must succeed: {e1}");
    let out_fq = out.join("a/a_dedup_R1.fastq.gz");
    let after1 = std::fs::read_to_string(&out_fq).unwrap();

    // Overwrite the input at the SAME path with different, longer reads.
    std::fs::write(&fq, "@r0\nTTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIII\n").unwrap();

    let (stdout2, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap()]);
    assert_eq!(c2, 0, "second run must succeed: {e2}");
    assert!(
        !stdout2.contains("resumed a/fastdedup"),
        "a changed input must invalidate the resume, not reuse stale output: {stdout2}"
    );
    let after2 = std::fs::read_to_string(&out_fq).unwrap();
    assert_ne!(after1, after2, "the output must be regenerated from the new reads");

    // A third run with the input unchanged DOES resume (the fix must not
    // disable resume, only make it input-sensitive).
    let (stdout3, e3, c3) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap()]);
    assert_eq!(c3, 0, "third run must succeed: {e3}");
    assert!(
        stdout3.contains("resumed a/fastdedup"),
        "an unchanged input must still resume: {stdout3}"
    );
}

/// A two-step config re-run into the same outdir: cleanup deletes each
/// intermediate FASTQ but used to leave the step's own --json behind, and
/// nothing passed --force, so fastprep's clobber guard aborted the second run
/// with exit 1. `--from` takes the same path.
#[test]
fn two_step_chain_runs_twice_into_the_same_outdir() {
    require_tool("fastdedup");
    require_tool("fastprep");
    let d = tmpdir("chain_twice");
    let out = d.join("out");
    let fq = d.join("s1_R1.fastq");
    std::fs::write(
        &fq,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let cfg = d.join("chain.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"chain\"\noutdir = \"{}\"\nsteps = [\"fastdedup\", \"fastprep\"]\nthreads = 2\ngz_level = 1\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();

    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "first run must succeed: {e1}");
    let (_, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(
        c2, 0,
        "a second run into the same outdir must succeed — the step's stale --json \
         and outputs must not trip its clobber guard: {e2}"
    );
    assert!(
        out.join("S1/S1_clean_R1.fastq.gz").is_file(),
        "the chain's final output must be there after the rerun"
    );
    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    let cmds = man["commands"].as_array().unwrap();
    assert_eq!(cmds.len(), 2, "one entry per step: {man}");
    for c in cmds {
        assert_eq!(c["exit"], 0, "both steps must have run clean on the rerun: {c}");
    }

    // Forcing the rerun from the second step must not re-trip the guard either.
    let (_, e3, c3) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", td, "--from", "fastprep",
    ]);
    assert_eq!(c3, 0, "--from must rerun the step cleanly: {e3}");
}

/// A step that writes no FASTQ, sitting mid-chain: it used to hand an empty
/// path to the next step (and let the cleanup delete the FASTQ that step still
/// needed). The chain must keep carrying the last FASTQ output forward.
#[test]
fn mid_chain_non_fastq_step_keeps_the_last_fastq_as_input() {
    require_tool("fastdedup");
    require_tool("fastcheck");
    require_tool("fastprep");
    let d = tmpdir("mid_nonfastq");
    let out = d.join("out");
    let fq = d.join("s1_R1.fastq");
    std::fs::write(
        &fq,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let cfg = d.join("mid.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"mid\"\noutdir = \"{}\"\nsteps = [\"fastdedup\", \"fastcheck\", \"fastprep\"]\nthreads = 2\ngz_level = 1\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "the chain must survive a mid-chain report-only step: {stderr}");
    assert!(out.join("S1/fastcheck.json").is_file(), "fastcheck must have run: {stderr}");

    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    let prep = man["commands"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["tool"] == "fastprep")
        .expect("fastprep must have run");
    let argv: Vec<String> = serde_json::from_value(prep["argv"].clone()).unwrap();
    let expected = out.join("S1/S1_dedup_R1.fastq.gz").display().to_string();
    assert!(
        argv.iter().all(|a| !a.is_empty()),
        "no step may be handed an empty path: {argv:?}"
    );
    assert!(
        argv.contains(&expected),
        "fastprep must read the last FASTQ-producing step's output ({expected}): {argv:?}"
    );
}

/// A step whose binary does not resolve must fail the run, not be skipped:
/// the run used to continue on the previous input and exit 0 while `check`
/// rejected the identical config.
#[test]
fn missing_step_binary_fails_the_run_like_check() {
    let d = tmpdir("missing_bin");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    // Only fastprep resolves here, so the run gets past its "no binaries at
    // all" guard and reaches the per-step resolution this test is about.
    let fake = tools.join("fastprep");
    std::fs::write(&fake, "#!/bin/sh\nexit 0\n").unwrap();
    let mut perm = std::fs::metadata(&fake).unwrap().permissions();
    std::os::unix::fs::PermissionsExt::set_mode(&mut perm, 0o755);
    std::fs::set_permissions(&fake, perm).unwrap();

    let fq = d.join("s1_R1.fastq");
    std::fs::write(&fq, "@r1\nACGT\n+\nIIII\n").unwrap();
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"m\"\noutdir = \"{}\"\nsteps = [\"fastcheck\", \"fastprep\"]\n\n[[samples]]\nname = \"s1\"\nr1 = \"{}\"\n",
            d.join("out").display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, run_err, run_code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap(),
    ]);
    let (_, check_err, check_code) = run(&[
        "check", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap(),
    ]);
    assert_eq!(
        run_code == 0,
        check_code == 0,
        "run must agree with check about an unresolvable step binary: \
         run {run_code} ({run_err}) / check {check_code} ({check_err})"
    );
    // (Only reachable when fastcheck is absent from tools/ and $PATH.)
    if run_code != 0 {
        assert!(
            run_err.contains("fastcheck"),
            "the failure must name the step: {run_err}"
        );
    }
}

/// A rerun must merge the manifest with what is on disk, not overwrite it:
/// a step this run did not reach still has artifacts in the outdir, and its
/// provenance used to be destroyed by the next run.
#[test]
fn manifest_keeps_entries_from_a_previous_run() {
    require_tool("fastdedup");
    require_tool("fastprep");
    let d = tmpdir("manifest_merge");
    let out = d.join("out");
    let fq = d.join("s1_R1.fastq");
    std::fs::write(
        &fq,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let body = |steps: &str| {
        format!(
            "[run]\nname = \"merge\"\noutdir = \"{}\"\nsteps = [{}]\nthreads = 2\ngz_level = 1\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            steps,
            fq.display()
        )
    };
    let both = d.join("both.toml");
    std::fs::write(&both, body("\"fastdedup\", \"fastprep\"")).unwrap();
    let one = d.join("one.toml");
    std::fs::write(&one, body("\"fastdedup\"")).unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();

    let (_, e1, c1) = run(&["run", "-c", both.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "{e1}");
    let (_, e2, c2) = run(&["run", "-c", one.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c2, 0, "{e2}");

    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    let tools: Vec<&str> = man["commands"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|c| c["tool"].as_str())
        .collect();
    assert!(
        tools.contains(&"fastprep"),
        "the previous run's fastprep entry must survive a rerun that does not \
         reach it: {man}"
    );
    assert_eq!(tools.iter().filter(|t| **t == "fastdedup").count(), 1, "{man}");
}

// ---------------------------------------------------------------- FH-2 ----

#[test]
fn check_refuses_duplicate_sample_names() {
    // FH-2 regression: two [[samples]] entries with the same name were taken
    // verbatim, yet the work dir and the resume key derive from the name alone
    // (runner.rs) — both workers created, unlinked and wrote the same paths.
    // The validator meant to catch this class passed it with exit 0.
    let d = tmpdir("fh2_dup");
    let a = d.join("a_R1.fastq");
    let b = d.join("b_R1.fastq");
    std::fs::write(&a, "@r1\nACGT\n+\nIIII\n").unwrap();
    std::fs::write(&b, "@r1\nACGT\n+\nIIII\n").unwrap();
    let cfg = d.join("dup.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"dup\"\nsteps = []\n\n[[samples]]\nname = \"lib1\"\nr1 = \"{}\"\n\n\
             [[samples]]\nname = \"lib1\"\nr1 = \"{}\"\n",
            a.display(),
            b.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&["check", "-c", cfg.to_str().unwrap()]);
    assert_ne!(code, 0, "a duplicate sample name must fail: {stderr}");
    assert!(
        stderr.contains("duplicate sample name") && stderr.contains("lib1"),
        "the refusal must name the duplicate: {stderr}"
    );

    // The same config with distinct names still validates: the new check must
    // not reject anything that was legitimately valid.
    let ok = d.join("ok.toml");
    std::fs::write(
        &ok,
        format!(
            "[run]\nname = \"ok\"\nsteps = [\"fastcheck\"]\n\n[steps.fastcheck]\nenabled = false\n\n[[samples]]\nname = \"lib1\"\nr1 = \"{}\"\n\n\
             [[samples]]\nname = \"lib2\"\nr1 = \"{}\"\n",
            a.display(),
            b.display()
        ),
    )
    .unwrap();
    let (_, stderr, code) = run(&["check", "-c", ok.to_str().unwrap()]);
    assert_eq!(code, 0, "distinct names must stay valid: {stderr}");
}

#[test]
fn check_still_accepts_the_repo_example_config() {
    // The duplicate-name refusal must not change the verdict on the config the
    // repo ships: `examples/fastsuite.toml` (simulate-only, three steps).
    for tool in ["fastcheck", "fastdedup", "fastprep"] {
        require_tool(tool);
    }
    let example = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("examples/fastsuite.toml");
    let (stdout, stderr, code) = run(&[
        "check",
        "-c",
        example.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "examples/fastsuite.toml must validate: {stderr}");
    assert!(
        stdout.contains("OK") && !stderr.contains("duplicate"),
        "{stdout}{stderr}"
    );
}

// ---------------------------------------------------------------- FH-3 ----

#[test]
fn report_json_keeps_its_alphabetical_field_order() {
    // The writer streams a borrowed view instead of building a Value tree, so
    // the order is now a struct declaration order rather than a BTreeMap's.
    // Alphabetical order is part of the file's stability contract.
    let d = tmpdir("fh3_order");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "sampleA.json", 100, 0.41);
    let out = d.join("out");
    let (_, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "json",
    ]);
    assert_eq!(code, 0, "{stderr}");
    let text = std::fs::read_to_string(out.join("report.json")).unwrap();
    let mut last = 0usize;
    for key in ["\"generator\":", "\"generator_version\":", "\"reports\":", "\"schema\":", "\"title\":"] {
        let at = text.find(key).unwrap_or_else(|| panic!("{key} missing"));
        assert!(at >= last, "{key} is out of alphabetical order:\n{text}");
        last = at;
    }
    for key in ["\"metrics\":", "\"notes\":", "\"sample\":", "\"source\":", "\"tool\":", "\"tool_version\":"] {
        assert!(text.contains(key), "{key} missing");
    }
    assert!(text.ends_with("}\n"), "the document must end with a newline");
}

// ---------------------------------------------------------------- FH-6 ----

#[test]
fn positions_all_metrics_stdout_matches_the_written_file() {
    // The tidy table is streamed to a BufWriter now; the file and the stdout
    // path must still produce the same bytes.
    let d = tmpdir("fh6_stdout");
    let src = d.join("p1");
    std::fs::create_dir_all(&src).unwrap();
    let mut tsv = String::from("pos\tmean_q\tdepth\n");
    for pos in 1..=5 {
        tsv.push_str(&format!("{pos}\t{}\t{}\n", 30.0 + pos as f64 / 4.0, 100 + pos));
    }
    std::fs::write(src.join("positions.tsv"), tsv).unwrap();

    let out = d.join("out");
    let (_, stderr, code) = run(&[
        "concatenate",
        src.to_str().unwrap(),
        "--positions",
        "--all-metrics",
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{stderr}");
    let written = std::fs::read_to_string(out.join("positions_matrix.tsv")).unwrap();

    let (stdout, stderr, code) = run(&[
        "concatenate",
        src.to_str().unwrap(),
        "--positions",
        "--all-metrics",
        "--stdout",
    ]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(stdout, written, "--stdout must be byte-identical to the file");
}

// ---------------------------------------------------------------- FH-8 ----

#[test]
fn concatenate_stdout_matches_the_written_tsv() {
    // --stdout used to build the matrix a second time; it now reuses the one
    // the writer built. The bytes must not change either way.
    let d = tmpdir("fh8_stdout");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "sampleA.json", 100, 0.41);
    write_report(&reports, "sampleB.json", 200, 0.52);
    let out = d.join("out");
    let (_, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "tsv",
    ]);
    assert_eq!(code, 0, "{stderr}");
    let written = std::fs::read_to_string(out.join("report.tsv")).unwrap();

    let (stdout, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "--format",
        "tsv",
        "--stdout",
    ]);
    assert_eq!(code, 0, "{stderr}");
    // stdout carries the table followed by the status line; the table itself
    // must be exactly the file.
    let table_end = stdout
        .find("fasthelper: concatenated")
        .unwrap_or_else(|| panic!("status line missing: {stdout}"));
    assert_eq!(
        &stdout[..table_end],
        written,
        "--stdout must be byte-identical to report.tsv"
    );
}

// ---------------------------------------------------------------- FH-9 ----

/// A stand-in fastdedup that records every invocation in `LOG` before doing
/// anything else, so `--version` forks are countable.
const COUNTING_DEDUP: &str = r#"#!/usr/bin/env bash
echo "$@" >> "@LOG@"
if [ "$1" = "--version" ]; then echo "fastdedup 9.9.9"; exit 0; fi
exit 0
"#;

#[test]
fn tool_version_is_resolved_once_per_run_not_once_per_sample() {
    // FH-9 regression: `--version` was forked per sample per step (30 000
    // forks at 10 k samples), including on --dry-run.
    let d = tmpdir("fh9_version");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    let log = d.join("invocations.log");
    let bin = tools.join("fastdedup");
    std::fs::write(&bin, COUNTING_DEDUP.replace("@LOG@", log.to_str().unwrap())).unwrap();
    let mut perm = std::fs::metadata(&bin).unwrap().permissions();
    std::os::unix::fs::PermissionsExt::set_mode(&mut perm, 0o755);
    std::fs::set_permissions(&bin, perm).unwrap();

    let mut cfg = String::from(
        "[run]\nname = \"v\"\nsteps = [\"fastdedup\"]\nthreads = 1\ngz_level = 1\n",
    );
    for s in ["s1", "s2", "s3"] {
        let fq = d.join(format!("{s}_R1.fastq"));
        std::fs::write(&fq, "@r0\nACGT\n+\nIIII\n").unwrap();
        cfg.push_str(&format!(
            "\n[[samples]]\nname = \"{s}\"\nr1 = \"{}\"\n",
            fq.display()
        ));
    }
    let path = d.join("cfg.toml");
    std::fs::write(&path, cfg).unwrap();

    let (stdout, stderr, code) = run(&[
        "run",
        "-c",
        path.to_str().unwrap(),
        "--tools-dir",
        tools.to_str().unwrap(),
        "--dry-run",
    ]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        stdout.matches("[dry-run]").count(),
        3,
        "all three samples must still be planned: {stdout}"
    );

    let calls = std::fs::read_to_string(&log).unwrap_or_default();
    let versions = calls.lines().filter(|l| l.trim() == "--version").count();
    assert_eq!(
        versions, 1,
        "one `--version` per distinct tool, not one per sample: {calls}"
    );
    assert_eq!(
        calls.lines().count(),
        versions,
        "dry-run must not execute any step: {calls}"
    );
}

// ------------------------------------------- orchestration (P0-P5) -------- //

/// Write a shell script and make it executable: the stand-in suite binaries
/// these tests point `--tools-dir` at.
fn write_exe(path: &Path, body: &str) {
    std::fs::write(path, body).unwrap();
    let mut perm = std::fs::metadata(path).unwrap().permissions();
    std::os::unix::fs::PermissionsExt::set_mode(&mut perm, 0o755);
    std::fs::set_permissions(path, perm).unwrap();
}

/// One sample, one paired FASTQ, reusable across the orchestration tests.
fn sample_fastq(dir: &Path) -> PathBuf {
    let fq = dir.join("s1_R1.fastq");
    std::fs::write(
        &fq,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    fq
}

#[test]
fn aggregate_ignores_reports_left_by_a_previous_configuration() {
    // P0: the aggregate re-walked the whole outdir, so run B — a chain WITHOUT
    // fastmaps — inherited run A's fastmaps row. A paper's table could carry a
    // measurement from a different configuration.
    require_tool("fastdedup");
    require_tool("fastprep");
    require_tool("fastmaps");
    let d = tmpdir("stale_aggregate");
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = |name: &str, steps: &str| {
        format!(
            "[run]\nname = \"{name}\"\noutdir = \"{}\"\nsteps = [{steps}]\nthreads = 2\ngz_level = 1\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        )
    };
    let a = d.join("a.toml");
    std::fs::write(&a, cfg("A", "\"fastdedup\", \"fastprep\", \"fastmaps\"")).unwrap();
    let b = d.join("b.toml");
    std::fs::write(&b, cfg("B", "\"fastdedup\", \"fastprep\"")).unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();

    let (_, e1, c1) = run(&["run", "-c", a.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "run A must succeed: {e1}");
    let tsv_a = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    assert!(
        tsv_a.contains("fastmaps"),
        "run A's report must carry its fastmaps row: {tsv_a}"
    );

    let (_, e2, c2) = run(&["run", "-c", b.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c2, 0, "run B must succeed: {e2}");
    let tsv_b = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    assert!(
        tsv_b.contains("fastdedup") && tsv_b.contains("fastprep"),
        "run B's own steps must still be reported: {tsv_b}"
    );
    assert!(
        !tsv_b.contains("fastmaps"),
        "run B does not run fastmaps, so its report must not carry run A's fastmaps row: {tsv_b}"
    );
}

/// A stand-in fastdedup that writes a partial output and a report, then fails —
/// the shape that used to be adopted by the next run as a completed step.
const FAKE_FAILING_DEDUP: &str = r#"#!/usr/bin/env bash
if [ "$1" = "--version" ]; then echo "fastdedup 9.9.9"; exit 0; fi
out=""; json=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    --json) json="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$(dirname "$out")" "$(dirname "$json")"
printf 'partial\n' > "$out"
printf '{"schema":"fastsuite.report.v1","tool":"fastdedup","tool_version":"9.9.9","metrics":{"duplicate_rate":0.0}}\n' > "$json"
exit 1
"#;

#[test]
fn resume_does_not_adopt_a_step_that_failed() {
    // P1: the resume entry was recorded before the exit code was tested, and a
    // resumed step was re-recorded with "exit": 0 — so a tool that left a
    // partial output and exited 1 was "resumed" as success, and its partial
    // output was fed to the next step.
    let d = tmpdir("resume_failed");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    write_exe(&tools.join("fastdedup"), FAKE_FAILING_DEDUP);
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"f\"\noutdir = \"{}\"\nsteps = [\"fastdedup\"]\nthreads = 1\ngz_level = 1\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();
    let td = tools.to_str().unwrap();

    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_ne!(c1, 0, "the failing step must fail the run: {e1} {c1}");

    let (stdout2, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert!(!stdout2.contains("resumed"), "a failed step must never be resumed: {stdout2}");
    assert_ne!(c2, 0, "the rerun must still fail, not report success: {e2}");

    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    for c in man["commands"].as_array().unwrap() {
        assert_ne!(
            c["exit"], 0,
            "the manifest must keep the recorded failure, not overwrite it with 0: {man}"
        );
    }
}

#[test]
fn the_documented_workflow_can_be_run_twice() {
    // P2: `fasthelper run` had no --force and the simulate step refused the
    // sim.gt.tsv the first run left behind, so the documented workflow exited 1
    // on the second run into the same outdir.
    require_tool("fastdedup");
    require_tool("fastprep");
    require_tool("fastmaps");
    require_tool("fastsim");
    let profile = tools_dir().join("fastsim/profiles/test_small.json");
    assert!(profile.is_file(), "the shipped test profile must exist: {}", profile.display());
    let d = tmpdir("documented_twice");
    let out = d.join("out");
    let cfg = d.join("fastsuite.toml");
    // The shape of examples/fastsuite.toml: the canonical chain plus a
    // simulated control sample, with per-step escape hatches.
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"demo\"\noutdir = \"{}\"\nthreads = 2\ngz_level = 1\n\
             steps = [\"fastdedup\", \"fastprep\", \"fastmaps\"]\n\n\
             [simulate]\nenabled = true\nname = \"sim_E\"\npairs = 200\nprofile = \"{}\"\nseed = 42\n\n\
             [steps.fastprep]\nargs = [\"--adapter-preset\", \"nextera\", \"--r2-scan-reverse\", \"--min-len\", \"0.6L\", \"--min-run\", \"6\"]\n\n\
             [steps.fastdedup]\nargs = [\"--pre-count\"]\n",
            out.display(),
            profile.display()
        ),
    )
    .unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();

    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "first run must succeed: {e1}");
    assert!(out.join("sim_E/sim.gt.tsv").is_file(), "the control sample must exist");

    let (_, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(
        c2, 0,
        "the documented workflow must be re-runnable into the same outdir: {e2}"
    );

    let (_, e3, c3) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td, "--force"]);
    assert_eq!(c3, 0, "--force must overwrite the previous run's outputs: {e3}");
    assert!(
        out.join("sim_E/sim_E_clean_R1.fastq.gz").is_file(),
        "the chain's final output must be there after the forced rerun"
    );
}

#[test]
fn unknown_step_fails_before_any_work() {
    // P3: an unknown step name was skipped with a warning while `check`
    // rejected the same config, so a typo silently dropped a step. It must be
    // refused up front, before a sample directory exists or a tool is forked.
    require_tool("fastcheck");
    let d = tmpdir("unknown_step");
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"t\"\noutdir = \"{}\"\nsteps = [\"fastcheck\", \"fastnope\"]\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();

    let (_, stderr, code) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_ne!(code, 0, "an unknown step must fail the run: {stderr}");
    assert!(
        stderr.contains("unknown step") && stderr.contains("fastnope"),
        "the refusal must name the unknown step: {stderr}"
    );
    assert!(
        !out.join("S1").exists(),
        "no sample may be started before the config is validated: {} exists",
        out.join("S1").display()
    );
    assert!(
        !out.join("report.tsv").exists(),
        "a rejected config must not write a report"
    );
}

/// A stand-in fastcheck: positional FASTQs in, a v1 report out via --json.
const FAKE_CHECK: &str = r#"#!/usr/bin/env bash
if [ "$1" = "--version" ]; then echo "fastcheck 9.9.9"; exit 0; fi
json=""
while [ $# -gt 0 ]; do
  case "$1" in
    --json) json="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$(dirname "$json")"
printf '{"schema":"fastsuite.report.v1","tool":"fastcheck","tool_version":"9.9.9","command":"fastcheck stub","metrics":{"reads":2}}\n' > "$json"
exit 0
"#;

#[test]
fn qc_only_config_needs_only_fastcheck() {
    // P4: the precondition demanded fastdedup OR fastprep OR fastmaps, so a
    // QC-only chain bailed with "no FastSuite binaries found" although the one
    // tool it needed was right there.
    let d = tmpdir("qc_only");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    write_exe(&tools.join("fastcheck"), FAKE_CHECK);
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"qc\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "a QC-only chain needs only fastcheck: {stderr}");
    assert!(out.join("S1/fastcheck.json").is_file(), "fastcheck must have run: {stderr}");
}

#[test]
fn dry_run_leaves_the_real_manifest_alone() {
    // P5: --dry-run rewrote manifest.json with entries claiming exit 0 and
    // elapsed 0.0, destroying the previous real run's provenance.
    require_tool("fastdedup");
    require_tool("fastprep");
    let d = tmpdir("dryrun_manifest");
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"p5\"\noutdir = \"{}\"\nsteps = [\"fastdedup\", \"fastprep\"]\nthreads = 2\ngz_level = 1\n\n[[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();

    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "{e1}");
    let before = std::fs::read_to_string(out.join("manifest.json")).unwrap();
    let man: serde_json::Value = serde_json::from_str(&before).unwrap();
    for c in man["commands"].as_array().unwrap() {
        assert!(
            c["elapsed_s"].as_f64().unwrap() > 0.0,
            "the real run must have measured something: {c}"
        );
    }

    let (_, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td, "--dry-run"]);
    assert_eq!(c2, 0, "{e2}");
    let after = std::fs::read_to_string(out.join("manifest.json")).unwrap();
    assert_eq!(
        before, after,
        "--dry-run must not rewrite the previous run's manifest"
    );
}

#[test]
fn concatenate_says_when_it_is_handed_an_aggregate() {
    // `collect` recognises per-tool reports only, so an aggregate used to fail
    // with "no recognisable reports found" — which denies that the file is a
    // report at all. It is one, of a different kind, and the message must say
    // so and name the input that can be merged.
    let d = tmpdir("concat_aggregate");
    let agg = d.join("report.json");
    std::fs::write(
        &agg,
        r#"{"generator":"fasthelper","generator_version":"0.1.0","reports":[]}"#,
    )
    .unwrap();
    let out = d.join("out");

    let (stdout, stderr, code) = run(&[
        "concatenate",
        agg.to_str().unwrap(),
        "--outdir",
        out.to_str().unwrap(),
        "--format",
        "tsv",
    ]);
    assert_ne!(code, 0, "an aggregate cannot be concatenated: {stdout}{stderr}");
    assert!(
        stderr.contains("aggregate"),
        "the refusal must identify the input as an aggregate, not as unrecognisable: {stderr}"
    );
    assert!(
        stderr.contains("report.json"),
        "it must name the offending file: {stderr}"
    );
    assert!(
        !out.join("report.tsv").exists(),
        "a refused concatenate must write nothing"
    );
}

// --------------------------------------------- 2026-10-21 review round ---- //

/// A two-read FASTQ, long enough for fastmaps to build a grid from.
fn write_two_read_fastq(path: &Path) {
    std::fs::write(
        path,
        "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
}

#[test]
fn two_sample_qc_report_keeps_each_sample_name() {
    // The generic-stem list had no `fastcheck` entry, so
    // `<outdir>/<sample>/fastcheck.json` resolved its sample to the stem
    // "fastcheck": two samples gave two identical rows in report.tsv and two
    // columns headed "fastcheck" in report.md/html. The run knows which sample
    // each report belongs to and must prefer that over path re-derivation.
    require_tool("fastcheck");
    let d = tmpdir("qc_two_samples");
    let out = d.join("out");
    let alpha = d.join("alpha_R1.fastq");
    let beta = d.join("beta_R1.fastq");
    write_two_read_fastq(&alpha);
    write_two_read_fastq(&beta);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"qc\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [[samples]]\nname = \"alpha\"\nr1 = \"{}\"\n\n\
             [[samples]]\nname = \"beta\"\nr1 = \"{}\"\n",
            out.display(),
            alpha.display(),
            beta.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "-P", "2",
        "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "a two-sample QC chain must run: {stderr}");

    let tsv = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    let rows: Vec<&str> = tsv.lines().skip(1).collect();
    assert_eq!(rows.len(), 2, "one row per sample: {tsv}");
    assert!(rows.iter().any(|r| r.starts_with("alpha\t")), "alpha lost its identity: {tsv}");
    assert!(rows.iter().any(|r| r.starts_with("beta\t")), "beta lost its identity: {tsv}");
    assert!(
        !rows.iter().any(|r| r.starts_with("fastcheck\t")),
        "the sample fell back to the generic report stem: {tsv}"
    );

    let html = std::fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        html.contains("<th>alpha</th>") && html.contains("<th>beta</th>"),
        "both samples must head their own column: {html}"
    );
    let md = std::fs::read_to_string(out.join("report.md")).unwrap();
    assert!(md.contains(" alpha |") && md.contains(" beta |"), "{md}");
}

#[test]
fn full_chain_with_a_fastsim_step_runs_twice_and_aggregates_every_step() {
    // fastsim as a chain step: it used to generate into the [simulate] control's
    // own prefix (clobber guard, exit 1) and, with no sim block, its outputs were
    // discarded (out1 empty) so the next step silently reprocessed the original
    // reads. Its artefacts now live in a step-private directory and become the
    // next step's input.
    for tool in ["fastsim", "fastdedup", "fastprep", "fastmaps", "fastcheck"] {
        require_tool(tool);
    }
    let profile = tools_dir().join("fastsim/profiles/test_small.json");
    assert!(profile.is_file(), "the shipped test profile must exist: {}", profile.display());
    let d = tmpdir("fastsim_chain");
    let out = d.join("out");
    let cfg = d.join("fastsuite.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"simchain\"\noutdir = \"{}\"\nthreads = 2\ngz_level = 1\n\
             steps = [\"fastsim\", \"fastdedup\", \"fastprep\", \"fastmaps\", \"fastcheck\"]\n\n\
             [simulate]\nenabled = true\nname = \"sim_E\"\npairs = 200\nprofile = \"{}\"\nseed = 42\n",
            out.display(),
            profile.display()
        ),
    )
    .unwrap();
    let td = tools_dir();
    let td = td.to_str().unwrap();
    let tools = ["fastsim", "fastdedup", "fastprep", "fastmaps", "fastcheck"];

    let (_, e1, c1) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c1, 0, "the chain with a fastsim step must run: {e1}");

    // Every step's report reached the aggregate, one row each.
    let tsv = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    let rows: Vec<&str> = tsv.lines().skip(1).collect();
    assert_eq!(rows.len(), tools.len(), "one row per step: {tsv}");
    for tool in tools {
        assert!(
            rows.iter().any(|r| *r == format!("sim_E\t{tool}") || r.starts_with(&format!("sim_E\t{tool}\t"))),
            "the {tool} step's report is missing from the aggregate: {tsv}"
        );
    }
    // The step generated into its own directory and the next step consumed it.
    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    let cmds = man["commands"].as_array().unwrap();
    assert_eq!(cmds.len(), tools.len(), "one manifest entry per step: {man}");
    for c in cmds {
        assert_eq!(c["exit"], 0, "every step must have run clean: {c}");
    }
    let dedup = cmds.iter().find(|c| c["tool"] == "fastdedup").unwrap();
    let argv: Vec<String> = serde_json::from_value(dedup["argv"].clone()).unwrap();
    let private_r1 = out.join("sim_E/fastsim/sim_R1.fastq.gz").display().to_string();
    assert!(
        argv.contains(&private_r1),
        "fastdedup must read the FASTQ the fastsim step generated ({private_r1}): {argv:?}"
    );

    // A second run into the same outdir, no --force: no clobber guard may trip.
    let (_, e2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_eq!(c2, 0, "the chain must be re-runnable into the same outdir: {e2}");
    let json2: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("report.json")).unwrap()).unwrap();
    assert_eq!(
        json2["reports"].as_array().unwrap().len(),
        tools.len(),
        "the rerun must aggregate every step: {json2}"
    );

    let (_, e3, c3) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td, "--force"]);
    assert_eq!(c3, 0, "--force must rerun the whole chain: {e3}");
    let tsv3 = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    assert_eq!(
        tsv3.lines().skip(1).count(),
        tools.len(),
        "the forced rerun must still aggregate one row per step: {tsv3}"
    );
}

/// A stand-in tool that exits 0 without writing anything: the shape that used
/// to leave a header-only report and exit 0.
const FAKE_SILENT_CHECK: &str = "#!/bin/sh\nexit 0\n";

#[test]
fn a_step_that_writes_no_json_fails_the_run() {
    let d = tmpdir("silent_step");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    write_exe(&tools.join("fastcheck"), FAKE_SILENT_CHECK);
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"silent\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools.to_str().unwrap(),
    ]);
    assert_ne!(
        code, 0,
        "a step that exited 0 without writing its --json must fail the run: {stderr}"
    );
    assert!(
        stderr.contains("fastcheck"),
        "the failure must name the step's tool: {stderr}"
    );
    assert!(
        stderr.contains("fastcheck.json"),
        "the failure must name the output that is missing: {stderr}"
    );
}

#[test]
fn fastcheck_assumptions_reach_the_aggregate() {
    // D-6: fastcheck states its assumptions as data. The aggregate dropped
    // every v1 top-level block, so the --no-quality caveat never reached a
    // published table.
    require_tool("fastcheck");
    let d = tmpdir("assumptions");
    let out = d.join("out");
    let fq = d.join("S1_R1.fastq");
    write_two_read_fastq(&fq);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"assumptions\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [steps.fastcheck]\nargs = [\"--no-quality\"]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{stderr}");

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("report.json")).unwrap()).unwrap();
    let row = json["reports"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["tool"] == "fastcheck")
        .expect("the fastcheck report must be aggregated");
    let notes = row["notes"].as_array().unwrap();
    assert!(
        notes
            .iter()
            .any(|n| n.as_str().unwrap_or_default().contains("no quality-derived metric is reported")),
        "the --no-quality assumption must be carried, not stripped: {row}"
    );
    let html = std::fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        html.contains("no quality-derived metric is reported"),
        "the caveat must be visible in the human-readable aggregate: {html}"
    );
}

#[test]
fn fastmaps_length_scalars_reach_the_aggregate() {
    // The corner fractions are length-biased, and fastmaps publishes the three
    // numbers that make the bias visible (min/max read length, corner window).
    // The grid reader dropped them, so the caveat was unpublishable downstream.
    require_tool("fastmaps");
    let d = tmpdir("maps_scalars");
    let out = d.join("out");
    let fq = d.join("S1_R1.fastq");
    write_two_read_fastq(&fq);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"maps\"\noutdir = \"{}\"\nsteps = [\"fastmaps\"]\nthreads = 1\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "{stderr}");

    let tsv = std::fs::read_to_string(out.join("report.tsv")).unwrap();
    let mut lines = tsv.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let row: Vec<&str> = lines.next().expect("one row").split('\t').collect();
    for key in ["min_seq_len", "max_seq_len", "corner_window"] {
        let i = header
            .iter()
            .position(|h| *h == key)
            .unwrap_or_else(|| panic!("{key} must be a column of the aggregate: {tsv}"));
        assert!(!row[i].is_empty(), "{key} must carry a value: {tsv}");
    }
    // 16 bp reads -> grid 17 -> corner window max(17/50, 2) = 2.
    let cell = |key: &str| {
        let i = header.iter().position(|h| *h == key).unwrap();
        row[i]
    };
    // The writer formats numeric cells with three decimals; compare values.
    let num = |key: &str| cell(key).parse::<f64>().unwrap();
    assert_eq!(num("min_seq_len"), 16.0, "{tsv}");
    assert_eq!(num("max_seq_len"), 16.0, "{tsv}");
    assert_eq!(num("corner_window"), 2.0, "{tsv}");
}

#[test]
fn from_an_unknown_step_fails_before_any_work() {
    // --from nosuchstep silently resumed everything and exited 0.
    require_tool("fastcheck");
    let d = tmpdir("from_unknown");
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"f\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--from", "nosuchstep",
        "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "an unknown --from step must fail the run: {stderr}");
    assert!(stderr.contains("nosuchstep"), "the refusal must name the step: {stderr}");
    assert!(
        !out.join("S1").exists(),
        "no sample may be started before --from is validated: {} exists",
        out.join("S1").display()
    );
    assert!(!out.join("report.tsv").exists(), "a rejected run must not write a report");
}

#[test]
fn check_refuses_a_sample_named_like_the_simulated_control() {
    // The control is appended to the sample list under [simulate] name, but the
    // duplicate check ran before that append: a [[samples]] entry called sim_E
    // put two workers on one directory and the control's report overwrote it.
    let d = tmpdir("dup_sim");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"dup\"\nsteps = [\"fastdedup\"]\n\n\
             [simulate]\nenabled = true\nname = \"sim_E\"\npairs = 10\n\n\
             [[samples]]\nname = \"sim_E\"\nr1 = \"{}\"\n",
            fq.display()
        ),
    )
    .unwrap();

    let (_, check_err, check_code) = run(&["check", "-c", cfg.to_str().unwrap()]);
    assert_ne!(check_code, 0, "check must refuse it: {check_err}");
    assert!(
        check_err.contains("duplicate sample name") && check_err.contains("sim_E"),
        "the refusal must name the collision: {check_err}"
    );

    let (_, run_err, run_code) = run(&["run", "-c", cfg.to_str().unwrap()]);
    assert_ne!(run_code, 0, "run must refuse it too: {run_err}");
    assert!(
        run_err.contains("duplicate sample name") && run_err.contains("sim_E"),
        "run must agree with check: {run_err}"
    );
}

#[test]
fn threads_zero_is_forwarded_to_every_tool() {
    // -j 0 means "all cores" to every tool, but the per-sample split turned it
    // into -j 1 before the children saw it (while the manifest recorded 0), and
    // fastcheck/fastmaps/fastsim never received -j at all.
    for tool in ["fastsim", "fastdedup", "fastprep", "fastmaps", "fastcheck"] {
        require_tool(tool);
    }
    let profile = tools_dir().join("fastsim/profiles/test_small.json");
    let d = tmpdir("threads_zero");
    let out = d.join("out");
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"j0\"\noutdir = \"{}\"\nthreads = 4\ngz_level = 1\n\
             steps = [\"fastsim\", \"fastdedup\", \"fastprep\", \"fastmaps\", \"fastcheck\"]\n\n\
             [simulate]\nenabled = true\nname = \"sim_E\"\npairs = 50\nprofile = \"{}\"\nseed = 42\n",
            out.display(),
            profile.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "-j", "0",
        "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "-j 0 must be accepted and the chain must run: {stderr}");

    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(man["threads"], 0, "the manifest must record what was asked for: {man}");
    let cmds = man["commands"].as_array().unwrap();
    assert_eq!(cmds.len(), 5, "{man}");
    for c in cmds {
        let argv: Vec<String> = serde_json::from_value(c["argv"].clone()).unwrap();
        let j = argv.iter().position(|a| a == "-j");
        let got = j
            .and_then(|i| argv.get(i + 1))
            .unwrap_or_else(|| panic!("{} must receive -j: {argv:?}", c["tool"]));
        if c["tool"] == "fastprep" {
            // fastprep documents --threads >= 1 and refuses 0, so the sentinel
            // is resolved for it rather than forwarded verbatim.
            assert!(got.parse::<usize>().unwrap() >= 1, "fastprep got -j {got}");
        } else {
            assert_eq!(got, "0", "{} must receive -j 0, not a rewritten 1: {argv:?}", c["tool"]);
        }
    }
}

// ------------------------------------------------ 2026-10-26 review round --- //

#[test]
fn compare_refuses_more_than_one_sample_per_side() {
    // H1: compare folded every report on a side into ONE flat metric map
    // ("last file wins"), so a before/ holding two samples produced a single
    // row carrying one sample's numbers, no sample names, and exit 0. compare
    // is a single-sample before/after tool: it must refuse, naming the samples.
    let d = tmpdir("compare_multi");
    for (side, sample, reads) in [
        ("before", "alpha", 100u64),
        ("before", "beta", 900),
        ("after", "alpha", 90),
    ] {
        let dir = d.join(side).join(sample);
        std::fs::create_dir_all(&dir).unwrap();
        write_report(&dir, "fastcheck.json", reads, 0.4);
    }
    let out = d.join("out");

    let (_, stderr, code) = run(&[
        "compare",
        "--before", d.join("before").to_str().unwrap(),
        "--after", d.join("after").to_str().unwrap(),
        "-o", out.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "two samples on one side must be refused: {stderr}");
    assert!(
        stderr.contains("alpha") && stderr.contains("beta"),
        "the refusal must name the samples it found: {stderr}"
    );
    assert!(
        stderr.contains("before"),
        "the refusal must name the side that holds them: {stderr}"
    );
    assert!(
        stderr.to_lowercase().contains("single-sample") || stderr.contains("one sample"),
        "it must say compare is a single-sample before/after tool: {stderr}"
    );
    assert!(
        !out.join("compare.tsv").exists(),
        "a refused compare must not publish a table built from one sample's numbers"
    );
}

#[test]
fn fail_on_with_an_unknown_metric_is_refused_with_the_available_metrics() {
    // H3: `--fail-on no_such_metric>0` matched nothing and exited 0 with no
    // diagnostic, so a typo silently disabled the gate.
    let d = tmpdir("failon_unknown");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "a.json", 100, 0.5);

    let (_, stderr, code) = run(&[
        "concatenate", reports.to_str().unwrap(),
        "-o", d.join("out").to_str().unwrap(),
        "--fail-on", "no_such_metric>0",
    ]);
    assert_ne!(code, 0, "an unknown --fail-on metric must be an error: {stderr}");
    assert!(
        stderr.contains("no_such_metric"),
        "the error must name the unknown metric: {stderr}"
    );
    for available in ["gc_fraction", "reads"] {
        assert!(
            stderr.contains(available),
            "the error must list the metrics that do exist ({available}): {stderr}"
        );
    }
    // A metric that does exist still gates as before.
    let (_, ok_err, ok_code) = run(&[
        "concatenate", reports.to_str().unwrap(),
        "-o", d.join("out2").to_str().unwrap(),
        "--fail-on", "gc_fraction>0.9",
    ]);
    assert_eq!(ok_code, 0, "a satisfied threshold must stay exit 0: {ok_err}");
}

#[test]
fn positions_refuses_fail_on_and_format_instead_of_dropping_them() {
    // H3: both flags were silently dropped under --positions: `--fail-on
    // no_such>0` exited 0 and `--format` produced none of the requested files.
    let d = tmpdir("positions_flags");
    let dir = d.join("A");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("positions.tsv"), "pos\tmean_q\n1\t30.0\n").unwrap();

    let (_, stderr, code) = run(&[
        "concatenate", "--positions", d.to_str().unwrap(),
        "-o", d.join("m").to_str().unwrap(),
        "--fail-on", "mean_q>0",
    ]);
    assert_ne!(code, 0, "a gate that cannot run must be refused: {stderr}");
    assert!(
        stderr.contains("fail-on") && stderr.contains("positions"),
        "the refusal must name the flag and the mode: {stderr}"
    );

    let (_, stderr2, code2) = run(&[
        "concatenate", "--positions", d.to_str().unwrap(),
        "-o", d.join("m2").to_str().unwrap(),
        "--format", "json",
    ]);
    assert_ne!(code2, 0, "an ignored --format must be refused: {stderr2}");
    assert!(
        stderr2.contains("format") && stderr2.contains("positions"),
        "the refusal must name the flag and the mode: {stderr2}"
    );

    // Neither refusal may leave a half-written matrix behind.
    assert!(!d.join("m/positions_matrix.tsv").exists());
    assert!(!d.join("m2/positions_matrix.tsv").exists());
}

#[test]
fn html_and_md_escape_tool_versions_and_metric_keys() {
    // H11: sample names and the title were escaped, but the tool name, its
    // version and the metric keys were interpolated raw into report.html and
    // report.md, so a report could inject markup into the published aggregate.
    let d = tmpdir("escape_tool");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fast<x>&tool",
        "tool_version": "9.9.9<b>",
        "metrics": {"m<key": 1.5, "reads": 10},
    });
    std::fs::write(reports.join("r.json"), serde_json::to_string(&doc).unwrap()).unwrap();

    let out = d.join("out");
    let (_, stderr, code) = run(&[
        "concatenate", reports.to_str().unwrap(),
        "-o", out.to_str().unwrap(),
        "--format", "html,md,json",
    ]);
    assert_eq!(code, 0, "{stderr}");

    let html = std::fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        !html.contains("<h2>fast<x>"),
        "the tool name must not reach the HTML raw: {html}"
    );
    assert!(
        html.contains("fast&lt;x&gt;&amp;tool"),
        "the tool name must be escaped in the HTML: {html}"
    );
    assert!(
        html.contains("9.9.9&lt;b&gt;"),
        "the tool version must be escaped in the HTML: {html}"
    );
    assert!(
        !html.contains("<code>m<key</code>"),
        "the metric key must not reach the HTML raw: {html}"
    );
    assert!(
        html.contains("m&lt;key"),
        "the metric key must be escaped in the HTML: {html}"
    );

    let md = std::fs::read_to_string(out.join("report.md")).unwrap();
    assert!(
        md.contains("## fast&lt;x&gt;&amp;tool"),
        "the tool heading must be escaped in the markdown: {md}"
    );
    assert!(
        !md.contains("## fast<x>&tool"),
        "the markdown heading must not carry raw markup: {md}"
    );
    assert!(
        md.contains("`m&lt;key`"),
        "the metric key must be escaped in the markdown: {md}"
    );

    // The JSON aggregate is data, not markup: it keeps the true names.
    let json = std::fs::read_to_string(out.join("report.json")).unwrap();
    assert!(json.contains("fast<x>&tool"), "report.json must stay raw: {json}");
}

#[test]
fn aggregate_carries_the_resolved_thread_count_a_step_reported() {
    // Cross-cutting: fastdedup reports the thread count it actually resolved
    // (`-j 0` = all cores), which is not the number the manifest records. The
    // aggregate must carry it, so a published table can state the parallelism.
    let d = tmpdir("agg_threads");
    let with = d.join("with");
    let without = d.join("without");
    std::fs::create_dir_all(&with).unwrap();
    std::fs::create_dir_all(&without).unwrap();
    for (dir, metrics) in [
        (&with, serde_json::json!({"reads_in": 10, "reads_out": 10, "threads": 7})),
        (&without, serde_json::json!({"reads_in": 10})),
    ] {
        let doc = serde_json::json!({
            "schema": "fastsuite.report.v1",
            "tool": "fastdedup",
            "tool_version": "1.3.0",
            "metrics": metrics,
        });
        std::fs::write(dir.join("fastdedup.json"), serde_json::to_string(&doc).unwrap()).unwrap();
    }

    let out = d.join("out");
    let (_, stderr, code) = run(&["concatenate", with.to_str().unwrap(), "-o", out.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let agg: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("report.json")).unwrap()).unwrap();
    assert_eq!(agg["threads"], 7.0, "the reported thread count must be carried: {agg}");

    let out2 = d.join("out2");
    let (_, stderr, code) = run(&["concatenate", without.to_str().unwrap(), "-o", out2.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let agg2: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out2.join("report.json")).unwrap()).unwrap();
    assert!(
        agg2.get("threads").is_none(),
        "no step reported threads, so the aggregate must not invent one: {agg2}"
    );
}

#[test]
fn fasthelpers_own_aggregate_passes_check_report() {
    // H12: report.json was stamped `schema: fastsuite.report.v1` while not
    // being a v1 report (no tool, no metrics), so `check --report` rejected
    // fasthelper's own output with four problems.
    let d = tmpdir("check_aggregate");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "a.json", 100, 0.41);
    let out = d.join("out");
    let (_, stderr, code) = run(&["concatenate", reports.to_str().unwrap(), "-o", out.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");

    let agg = out.join("report.json");
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&agg).unwrap()).unwrap();
    assert_eq!(
        doc["schema"], "fastsuite.aggregate.v1",
        "an aggregate is not a fastsuite.report.v1 report: {doc}"
    );

    let (stdout, stderr, code) = run(&["check", "--report", agg.to_str().unwrap()]);
    assert_eq!(code, 0, "check must accept fasthelper's own aggregate: {stderr}");
    assert!(stdout.contains("aggregate"), "the verdict must say what it checked: {stdout}");
}

#[test]
fn check_report_refuses_a_broken_aggregate() {
    let d = tmpdir("check_aggregate_bad");
    // A row without a tool: exactly the shape that made an aggregate unusable.
    let bad = d.join("report.json");
    std::fs::write(
        &bad,
        r#"{"schema":"fastsuite.aggregate.v1","generator":"fasthelper",
            "generator_version":"0.1.0","title":"t",
            "reports":[{"sample":"A","metrics":{"reads":1}}]}"#,
    )
    .unwrap();
    let (_, stderr, code) = run(&["check", "--report", bad.to_str().unwrap()]);
    assert_ne!(code, 0, "an aggregate row without a tool must be refused: {stderr}");
    assert!(
        stderr.contains("reports[0]") && stderr.contains("tool"),
        "the problem must locate the row and name the missing key: {stderr}"
    );
}

/// A stand-in fastcheck that exits 0 and writes a CORRUPT `--json` report: the
/// shape that used to become a zero-metric "unparsable" row and be resumed
/// forever (the resume signature is size+mtime, and neither ever changed).
const FAKE_CORRUPT_CHECK: &str = r#"#!/usr/bin/env bash
if [ "$1" = "--version" ]; then echo "fastcheck 9.9.9"; exit 0; fi
json=""
while [ $# -gt 0 ]; do
  case "$1" in
    --json) json="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$(dirname "$json")"
printf 'this is not json\n' > "$json"
exit 0
"#;

#[test]
fn a_step_whose_report_does_not_parse_fails_and_is_not_resumed() {
    // H2: a step that exits 0 with a corrupt --json passed the existence-only
    // check, was published as a zero-metric "unparsable" row, and was resumed
    // on every later run because size+mtime never changed.
    let d = tmpdir("corrupt_report");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    write_exe(&tools.join("fastcheck"), FAKE_CORRUPT_CHECK);
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"corrupt\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();
    let td = tools.to_str().unwrap();

    let (_, stderr, code) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert_ne!(code, 0, "an unparsable report must fail the step: {stderr}");
    assert!(
        stderr.contains("S1/fastcheck.json"),
        "the failure must name the report path: {stderr}"
    );
    assert!(
        stderr.contains("not valid JSON") || stderr.contains("unparsable"),
        "the failure must say the report does not parse: {stderr}"
    );
    // It must not be published as a measurement of anything.
    let agg = std::fs::read_to_string(out.join("report.json")).unwrap_or_default();
    assert!(
        !agg.contains("unparsable"),
        "a step whose report does not parse must not reach the aggregate: {agg}"
    );

    // And it must not be adopted by the next run.
    let (stdout2, stderr2, c2) = run(&["run", "-c", cfg.to_str().unwrap(), "--tools-dir", td]);
    assert!(
        !stdout2.contains("resumed"),
        "an unparsable report must never be resumed: {stdout2}"
    );
    assert_ne!(c2, 0, "the rerun must still fail: {stderr2}");
}

#[test]
fn from_a_disabled_step_is_refused() {
    // H7: `--from fastcheck` on a config where fastcheck is disabled was a
    // silent no-op — the step was skipped before `from_reached` was set, so
    // every step resumed and nothing reran, with exit 0.
    require_tool("fastcheck");
    let d = tmpdir("from_disabled");
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"f\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [steps.fastcheck]\nenabled = false\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run", "-c", cfg.to_str().unwrap(), "--from", "fastcheck",
        "--tools-dir", tools_dir().to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "--from a disabled step must be refused: {stderr}");
    assert!(
        stderr.contains("fastcheck") && stderr.contains("disabled"),
        "the refusal must name the step and say it is disabled: {stderr}"
    );
    assert!(
        !out.join("S1").exists(),
        "no sample may be started before --from is validated"
    );
}

#[test]
fn the_default_tools_dir_is_the_tools_directory() {
    // H9: the default resolved to tools/fasthelper (one parent short), so a run
    // without --tools-dir bailed "looked in .../tools/fasthelper" although the
    // binaries sit in tools/<tool>/target/release.
    require_tool("fastcheck");
    let d = tmpdir("default_tools_dir");
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"dt\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let o = Command::new(fasthelper())
        .args(["run", "-c", cfg.to_str().unwrap()])
        .env_remove("FASTSUITE_TOOLS_DIR")
        .env("PATH", "")
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert_eq!(o.status.code(), Some(0), "the default must find tools/: {stderr}");

    let man: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    let binary = man["commands"][0]["binary"].as_str().unwrap();
    assert!(
        Path::new(binary).starts_with(tools_dir()),
        "the resolved binary must come from tools/, not tools/fasthelper: {binary}"
    );
}

#[test]
fn doctor_exits_non_zero_when_tools_are_missing() {
    // H10: doctor printed "missing: ..." and exited 0.
    let d = tmpdir("doctor_missing");
    let empty = d.join("tools");
    std::fs::create_dir_all(&empty).unwrap();

    let o = Command::new(fasthelper())
        .args(["doctor", "--tools-dir", empty.to_str().unwrap()])
        .env_remove("FASTSUITE_TOOLS_DIR")
        .env("PATH", "")
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert_ne!(
        o.status.code(),
        Some(0),
        "doctor must fail when every tool is missing: {stderr}"
    );
    assert!(
        stderr.contains("fastcheck") && stderr.contains("not found"),
        "the failure must name the missing tools: {stderr}"
    );

    // With the suite built, the same command is a clean pass.
    let o = Command::new(fasthelper())
        .args(["doctor", "--tools-dir", tools_dir().to_str().unwrap()])
        .env_remove("FASTSUITE_TOOLS_DIR")
        .env("PATH", "")
        .output()
        .unwrap();
    let out = String::from_utf8_lossy(&o.stdout).to_string();
    assert_eq!(
        o.status.code(),
        Some(0),
        "a complete install must pass: {out}{}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(!out.contains("NOT FOUND"), "{out}");
}

#[test]
fn install_says_what_it_could_not_link_and_exits_non_zero() {
    // H10: install exited 0 having linked 1 of 6, saying nothing about the rest.
    let d = tmpdir("install_missing");
    let prefix = d.join("prefix");
    std::fs::create_dir_all(&prefix).unwrap();
    let empty = d.join("tools");
    std::fs::create_dir_all(&empty).unwrap();

    let o = Command::new(fasthelper())
        .args([
            "install", "--prefix", prefix.to_str().unwrap(),
            "--tools-dir", empty.to_str().unwrap(),
        ])
        .env_remove("FASTSUITE_TOOLS_DIR")
        .env("PATH", "")
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert_ne!(
        o.status.code(),
        Some(0),
        "install must fail when it cannot link the suite: {stderr}"
    );
    for missing in ["fastprep", "fastmaps"] {
        assert!(
            stderr.contains(missing),
            "install must name what it could not link ({missing}): {stderr}"
        );
    }

    // A complete install still succeeds, and now really links all six + legacy.
    for tool in ["fastcheck", "fastdedup", "fastprep", "fastmaps", "fastsim"] {
        require_tool(tool);
    }
    let o = Command::new(fasthelper())
        .args([
            "install", "--prefix", prefix.to_str().unwrap(),
            "--tools-dir", tools_dir().to_str().unwrap(),
        ])
        .env_remove("FASTSUITE_TOOLS_DIR")
        .env("PATH", "")
        .output()
        .unwrap();
    let combined = format!(
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    assert_eq!(o.status.code(), Some(0), "a complete install must pass: {combined}");
    for name in [
        "fastcheck", "fastdedup", "fastprep", "fastmaps", "fastsim", "fasthelper",
        "rpolyx", "fdedup", "novaseq-sim",
    ] {
        assert!(
            prefix.join("bin").join(name).is_symlink(),
            "{name} must be linked into the prefix: {combined}"
        );
    }
}

#[test]
fn check_validates_the_simulate_block() {
    // H6: `check -c` returned OK for a config `run` rejects at minute zero —
    // the [simulate] fastsim binary and profile path were never validated.
    let d = tmpdir("check_simulate");
    let empty = d.join("tools");
    std::fs::create_dir_all(&empty).unwrap();
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        "[run]\nname = \"sim\"\nsteps = []\n\n\
         [simulate]\nenabled = true\npairs = 10\nprofile = \"nope_profile.json\"\n",
    )
    .unwrap();

    let o = Command::new(fasthelper())
        .args([
            "check", "-c", cfg.to_str().unwrap(),
            "--tools-dir", empty.to_str().unwrap(),
        ])
        .env_remove("FASTSUITE_TOOLS_DIR")
        .env("PATH", "")
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert_ne!(o.status.code(), Some(0), "check must refuse it: {stderr}");
    assert!(
        stderr.contains("nope_profile.json"),
        "the missing simulate profile must be named: {stderr}"
    );
    assert!(
        stderr.contains("fastsim"),
        "the missing fastsim binary must be named: {stderr}"
    );
}

// ------------------------------------ 2026-10-26 contract/gate round ---- //

/// A stand-in fastcheck that writes a v1 report carrying a NESTED metric: the
/// shape `check --report` refuses (rule 1: a flat numeric namespace).
const FAKE_NESTED_METRIC_CHECK: &str = r#"#!/usr/bin/env bash
if [ "$1" = "--version" ]; then echo "fastcheck 9.9.9"; exit 0; fi
json=""
while [ $# -gt 0 ]; do
  case "$1" in
    --json) json="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$(dirname "$json")"
printf '{"schema":"fastsuite.report.v1","tool":"fastcheck","tool_version":"9.9.9","command":"fastcheck stub","metrics":{"reads":2,"adapter_content":{"r1":12}}}\n' > "$json"
exit 0
"#;

#[test]
fn a_step_whose_report_breaks_the_contract_fails_like_check_report() {
    // The production path ingested what the validator rejects: a report with
    // `"adapter_content": {"r1": 12}` was accepted by `concatenate` with exit 0
    // and the metric silently absent. The runner must apply the same rules.
    let d = tmpdir("nested_metric_run");
    let tools = d.join("tools");
    std::fs::create_dir_all(&tools).unwrap();
    write_exe(&tools.join("fastcheck"), FAKE_NESTED_METRIC_CHECK);
    let out = d.join("out");
    let fq = sample_fastq(&d);
    let cfg = d.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"nested\"\noutdir = \"{}\"\nsteps = [\"fastcheck\"]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            out.display(),
            fq.display()
        ),
    )
    .unwrap();

    let (_, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        tools.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "an invalid report must fail the step: {stderr}");
    assert!(
        stderr.contains("S1/fastcheck.json"),
        "the failure must name the report: {stderr}"
    );
    assert!(
        stderr.contains("adapter_content") && stderr.contains("not a number"),
        "the failure must carry the validator's reason: {stderr}"
    );
    let agg = std::fs::read_to_string(out.join("report.json")).unwrap_or_default();
    assert!(
        !agg.contains("adapter_content"),
        "a rejected report must not reach the aggregate: {agg}"
    );
}

#[test]
fn concatenate_names_every_metric_the_tolerant_reader_drops() {
    // The tolerant adapter still ingests such a file (concatenate is not the
    // gate), but a metric it cannot carry must be NAMED, never silently absent.
    let d = tmpdir("nested_metric_note");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastprep",
        "tool_version": "0.4.1",
        "command": "fastprep stub",
        "metrics": {
            "reads_in": 10, "reads_out": 10, "reads_dropped": 0,
            "adapter_content": {"r1": 12}
        }
    });
    std::fs::write(
        reports.join("fastprep.json"),
        serde_json::to_string(&doc).unwrap(),
    )
    .unwrap();
    let out = d.join("out");

    let (_, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "the tolerant path still ingests it: {stderr}");
    let agg: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("report.json")).unwrap()).unwrap();
    let row = &agg["reports"][0];
    assert!(
        row["metrics"].get("adapter_content").is_none(),
        "a nested value cannot be a metric: {agg}"
    );
    let notes = row["notes"].to_string();
    assert!(
        notes.contains("adapter_content") && notes.contains("dropped"),
        "the note must name the dropped metric: {notes}"
    );
}

/// A config whose only step carries `args`, for the `check -c` flag tests.
fn write_args_config(dir: &Path, args: &str) -> PathBuf {
    let fq = sample_fastq(dir);
    let cfg = dir.join("cfg.toml");
    std::fs::write(
        &cfg,
        format!(
            "[run]\nname = \"flags\"\noutdir = \"{}\"\nsteps = [\"fastprep\"]\n\n\
             [steps.fastprep]\nargs = [{args}]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            dir.join("out").display(),
            fq.display()
        ),
    )
    .unwrap();
    cfg
}

#[test]
fn check_refuses_a_flag_that_is_only_a_prefix_of_a_real_one() {
    // `--min-le` is a substring of `--min-len`, and the old test was
    // `help.contains(flag)`, so a truncated flag validated and the tool then
    // died on an unknown argument at the first step.
    require_tool("fastprep");
    let d = tmpdir("check_prefix_flag");
    let td = tools_dir();
    let cfg = write_args_config(&d, "\"--min-le\", \"100\"");
    let (_, stderr, code) = run(&[
        "check",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        td.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "a prefix of a real flag must be refused: {stderr}");
    assert!(
        stderr.contains("does not accept '--min-le'"),
        "the refusal must name the truncated flag: {stderr}"
    );

    // The whole flag it was truncated from still validates.
    let ok = write_args_config(&d, "\"--min-len\", \"100\"");
    let (_, ok_err, ok_code) = run(&[
        "check",
        "-c",
        ok.to_str().unwrap(),
        "--tools-dir",
        td.to_str().unwrap(),
    ]);
    assert_eq!(ok_code, 0, "a real flag must stay valid: {ok_err}");
}

#[test]
fn check_validates_every_flag_not_just_the_first() {
    // `.find()` returned only the first `--token`, so a typo in the SECOND one
    // was never looked at: the gate blessed a config the tool refuses.
    require_tool("fastprep");
    let d = tmpdir("check_second_flag");
    let td = tools_dir();
    let cfg = write_args_config(&d, "\"--min-len\", \"100\", \"--bogus-flag\", \"1\"");
    let (_, stderr, code) = run(&[
        "check",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        td.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "the second unknown flag must be refused: {stderr}");
    assert!(
        stderr.contains("does not accept '--bogus-flag'"),
        "the refusal must name the second flag: {stderr}"
    );
    assert!(
        !stderr.contains("does not accept '--min-len'"),
        "the valid first flag must not be blamed: {stderr}"
    );
}

#[test]
fn check_refuses_gz_level_threads_and_empty_steps_the_tools_refuse() {
    let d = tmpdir("check_ranges");
    let fq = sample_fastq(&d);
    let body = |run_extra: &str, steps: &str| {
        format!(
            "[run]\nname = \"r\"\noutdir = \"{}\"\n{run_extra}steps = [{steps}]\n\n\
             [[samples]]\nname = \"S1\"\nr1 = \"{}\"\n",
            d.join("out").display(),
            fq.display()
        )
    };
    let td = tools_dir();
    let td = td.to_str().unwrap();

    // gz_level 42: every tool bails on --gz-level outside 1..=9, and `check`
    // used to bless it.
    let gz = d.join("gz.toml");
    std::fs::write(&gz, body("gz_level = 42\n", "\"fastcheck\"")).unwrap();
    let (_, stderr, code) = run(&[
        "check",
        "-c",
        gz.to_str().unwrap(),
        "--tools-dir",
        td,
    ]);
    assert_ne!(code, 0, "gz_level 42 must be refused: {stderr}");
    assert!(
        stderr.contains("gz_level") && stderr.contains("[1, 9]"),
        "the refusal must state the range: {stderr}"
    );
    // `run` refuses the same value, before any tool is forked.
    let (_, run_err, run_code) = run(&[
        "run",
        "-c",
        gz.to_str().unwrap(),
        "--tools-dir",
        td,
    ]);
    assert_ne!(run_code, 0, "run must refuse it too: {run_err}");
    assert!(run_err.contains("gz_level"), "{run_err}");

    // A thread count no machine has: a typo gate, refused at minute zero.
    let th = d.join("threads.toml");
    std::fs::write(&th, body("threads = 100000\n", "\"fastcheck\"")).unwrap();
    let (_, stderr, code) = run(&[
        "check",
        "-c",
        th.to_str().unwrap(),
        "--tools-dir",
        td,
    ]);
    assert_ne!(code, 0, "an absurd thread count must be refused: {stderr}");
    assert!(stderr.contains("threads"), "{stderr}");

    // No steps at all: the run produced an empty report and exited 0.
    let empty = d.join("empty.toml");
    std::fs::write(&empty, body("", "")).unwrap();
    let (_, stderr, code) = run(&[
        "check",
        "-c",
        empty.to_str().unwrap(),
        "--tools-dir",
        td,
    ]);
    assert_ne!(code, 0, "an empty step list must be refused: {stderr}");
    assert!(
        stderr.contains("steps is empty"),
        "the refusal must say there is nothing to run: {stderr}"
    );
}

#[test]
fn compare_refuses_more_than_one_tool_per_side() {
    // One level below the multi-sample refusal: two tools of one sample were
    // merged last-wins, so one tool's numbers were published and the other's
    // dropped, with exit 0 and no hint that a report had been discarded.
    let d = tmpdir("compare_multitool");
    let before = d.join("before").join("alpha");
    let after = d.join("after").join("alpha");
    std::fs::create_dir_all(&before).unwrap();
    std::fs::create_dir_all(&after).unwrap();
    write_report(&before, "fastcheck.json", 100, 0.41);
    write_metric_report(&before, "fastdedup.json", "reads_in", 100.0);
    write_report(&after, "fastcheck.json", 90, 0.42);
    let out = d.join("out");

    let (_, stderr, code) = run(&[
        "compare",
        "--before",
        d.join("before").to_str().unwrap(),
        "--after",
        d.join("after").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "two tools on one side must be refused: {stderr}");
    assert!(
        stderr.contains("fastcheck") && stderr.contains("fastdedup"),
        "the refusal must name the tools it found: {stderr}"
    );
    assert!(
        stderr.contains("before"),
        "the refusal must name the side that holds them: {stderr}"
    );
    assert!(
        stderr.to_lowercase().contains("single-tool"),
        "it must say compare is a single-tool before/after tool: {stderr}"
    );
    assert!(
        !out.join("compare.tsv").exists(),
        "a refused compare must not publish a table built from one tool's numbers"
    );
}

#[test]
fn stdout_is_printed_before_a_threshold_exits_2() {
    // `process::exit(2)` ran before the `--stdout` block, so the gate that
    // tripped also threw away the table the caller asked for: exit 2, 0 bytes.
    let d = tmpdir("stdout_gate");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "a.json", 100, 0.41);
    let out = d.join("out");

    let (stdout, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "json",
        "--stdout",
        "--fail-on",
        "gc_fraction>0.1",
    ]);
    assert_eq!(code, 2, "a tripped threshold must exit 2: {stderr}");
    assert!(
        !stdout.trim().is_empty(),
        "--stdout must not be lost when the threshold trips: {stdout:?}"
    );
    assert!(
        stdout.contains("gc_fraction"),
        "--stdout must carry the table: {stdout}"
    );
    assert!(stderr.contains("threshold tripped"), "{stderr}");
}

#[test]
fn an_unknown_format_is_refused_before_anything_is_written() {
    // `--format json,bogus` wrote report.json and only then bailed on "bogus":
    // a usage error left a half-published directory behind.
    let d = tmpdir("format_bogus");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "a.json", 10, 0.5);
    let out = d.join("out");

    let (_, stderr, code) = run(&[
        "concatenate",
        reports.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "json,bogus",
    ]);
    assert_ne!(code, 0, "an unknown format must be an error: {stderr}");
    assert!(stderr.contains("bogus"), "{stderr}");
    assert!(
        !out.join("report.json").exists(),
        "a usage error must not publish a partial aggregate"
    );
    assert!(!out.join("report.tsv").exists());
}
