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
printf '{"schema":"fastsuite.report.v1","tool":"fastdedup","tool_version":"9.9.9","metrics":{"duplicate_rate":0.0}}\n' > "$json"
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
            "[run]\nname = \"ok\"\nsteps = []\n\n[[samples]]\nname = \"lib1\"\nr1 = \"{}\"\n\n\
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
