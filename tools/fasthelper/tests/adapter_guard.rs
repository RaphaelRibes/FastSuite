//! End-to-end tests for the post-trim residual-adapter guard.
//!
//! They drive the real `fasthelper run` against fixtures built inline, using the
//! suite's own `fastcheck`/`fastprep` binaries when they resolve under
//! `tools/`. These are real gates: if a required binary is absent the test FAILS
//! loudly (see `require_tools`) rather than passing without exercising anything.
//! The decision logic itself is also covered exhaustively by the unit tests in
//! `src/adapter_guard.rs`.

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

/// The binary `run` would resolve for `tool` under `tools`, or None.
fn resolve(tool: &str) -> Option<PathBuf> {
    let dir = tools_dir().join(tool).join("target");
    for profile in ["release", "debug"] {
        let p = dir.join(profile).join(tool);
        if p.is_file() {
            return Some(p);
        }
    }
    None
}

/// Fail LOUDLY when a required binary is missing. These are real end-to-end
/// gates, so an absent prerequisite must be a broken test, not a silent pass
/// that would go green without ever driving `fasthelper run`.
fn require_tools(tools: &[&str]) {
    let missing: Vec<&str> = tools
        .iter()
        .copied()
        .filter(|t| resolve(t).is_none())
        .collect();
    if !missing.is_empty() {
        panic!(
            "adapter_guard test needs {} built under {} — build them \
             (cd tools/<tool> && cargo build --release) before running this test",
            missing.join(", "),
            tools_dir().display()
        );
    }
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fh_ag_{}_{}", std::process::id(), tag));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Three reads carrying the Nextera adapter at the 3' end and three clean reads,
/// all six distinct, non-homopolymeric (so polyX trimming does not touch them).
fn fixture_fastq(path: &Path) {
    let adapter = "CTGTCTCTTATACACATCT"; // Nextera, 19 bp
    let fillers = [
        "ACGT".repeat(25),                 // clean 0
        "GATTACA".repeat(14) + "CC",       // clean 1
        "TTGCAACGTA".repeat(10),           // clean 2
        "CATGCATGCA".repeat(10),           // carries adapter
        "GGCTTAAGCC".repeat(10),           // carries adapter
        "ATATGCGCTA".repeat(10),           // carries adapter
    ];
    let mut reads: Vec<String> = Vec::new();
    for f in &fillers[..3] {
        reads.push(f.clone());
    }
    for f in &fillers[3..] {
        reads.push(format!("{f}{adapter}"));
    }
    let mut text = String::new();
    for (i, seq) in reads.iter().enumerate() {
        text.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(path, text).unwrap();
}

fn write_config(path: &Path, outdir: &Path, raw: &Path, steps: &str, prep_args: &str) {
    let cfg = format!(
        "[run]\nname = \"ag\"\noutdir = \"{}\"\nsteps = [{steps}]\n\n\
         [steps.fastprep]\nargs = [{prep_args}]\n\n\
         [[samples]]\nname = \"s1\"\nr1 = \"{}\"\n",
        outdir.display(),
        raw.display()
    );
    std::fs::write(path, cfg).unwrap();
}

fn run(args: &[&str]) -> (String, String, i32) {
    let out = Command::new(fasthelper())
        .args(args)
        .output()
        .expect("spawn fasthelper");
    (
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

fn read_manifest(outdir: &Path) -> serde_json::Value {
    let text = std::fs::read_to_string(outdir.join("manifest.json")).unwrap();
    serde_json::from_str(&text).unwrap()
}

/// Scenario 1 — the trimmed output still carries adapter: the guard warns with
/// the actionable cause, and can be made fatal on demand.
#[test]
fn residual_adapter_after_trimming_warns_with_an_actionable_cause() {
    require_tools(&["fastcheck", "fastprep"]);
    let d = tmpdir("residual");
    let raw = d.join("s1_R1.fastq");
    fixture_fastq(&raw);
    let out = d.join("run");
    let cfg = d.join("fastsuite.toml");
    // --no-adapter leaves the adapter in place: exactly the under-trim failure.
    write_config(&cfg, &out, &raw, "\"fastprep\"", "\"--no-adapter\"");

    let (stdout, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
    ]);

    assert_eq!(code, 0, "a warning must not fail the run: {stderr}");
    assert!(
        stderr.contains("adapter survived the trim"),
        "the guard must warn: {stderr}"
    );
    // The actionable causes the task names, by flag.
    assert!(stderr.contains("likely cause"), "{stderr}");
    assert!(stderr.contains("--adapter-error-rate"), "{stderr}");
    assert!(stderr.contains("--min-adapter-overlap"), "{stderr}");
    assert!(stderr.contains("--adapter-preset"), "{stderr}");
    assert!(
        stderr.contains("--no-adapter"),
        "the cause must name disabling the adapter set: {stderr}"
    );
    // Both measurements were taken and recorded (under logs/, so the aggregator
    // does not mistake them for samples).
    assert!(out.join("logs/s1.adapter-guard.raw.json").is_file());
    assert!(out.join("logs/s1.adapter-guard.trim.json").is_file());
    let man = read_manifest(&out);
    assert_eq!(man["adapter_guard"]["enabled"], true);
    assert_eq!(man["adapter_guard"]["tripped_samples"][0], "s1");
    let frac = man["commands"][0]["adapter_guard"]["raw"]["fraction"].as_f64().unwrap();
    assert!((frac - 0.5).abs() < 1e-9, "raw fraction should be 3/6: {man}");
    let _ = stdout;

    // --fail-on-adapter-residual turns the same trip into exit 2.
    let out2 = d.join("run_fail");
    let (_, stderr2, code2) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "-o",
        out2.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
        "--fail-on-adapter-residual",
    ]);
    assert_eq!(code2, 2, "the fail-on style condition must exit 2: {stderr2}");
    assert!(stderr2.contains("exiting 2"), "{stderr2}");
}

/// Scenario 2 — a clean trimmed output: the guard stays silent.
#[test]
fn clean_trimmed_output_is_silent() {
    require_tools(&["fastcheck", "fastprep"]);
    let d = tmpdir("clean");
    let raw = d.join("s1_R1.fastq");
    fixture_fastq(&raw);
    let out = d.join("run");
    let cfg = d.join("fastsuite.toml");
    // An explicit adapter set: fastprep removes the Nextera adapter it was told about.
    write_config(&cfg, &out, &raw, "\"fastprep\"", "\"--adapter-preset\", \"nextera\"");

    let (_, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
    ]);

    assert_eq!(code, 0, "{stderr}");
    assert!(
        !stderr.contains("adapter survived the trim"),
        "a clean trim must be silent: {stderr}"
    );
    assert!(
        !stderr.contains("adapter guard:"),
        "no guard noise on a clean trim: {stderr}"
    );
    let man = read_manifest(&out);
    assert!(
        man["adapter_guard"]["tripped_samples"].as_array().unwrap().is_empty(),
        "{man}"
    );
    // The measurement still ran (adapter was present in the input) and shows the fall.
    let raw_frac = man["commands"][0]["adapter_guard"]["raw"]["fraction"].as_f64().unwrap();
    let trim_frac = man["commands"][0]["adapter_guard"]["trimmed"]["fraction"].as_f64().unwrap();
    assert!(raw_frac > 0.4, "raw had adapter: {man}");
    assert!(trim_frac == 0.0, "trimmed is clean: {man}");
}

/// Scenario 3 — a run with no trimming step: the guard never fires.
#[test]
fn a_run_without_a_trimming_step_never_fires_the_guard() {
    require_tools(&["fastcheck"]);
    let d = tmpdir("notrim");
    let raw = d.join("s1_R1.fastq");
    fixture_fastq(&raw);
    let out = d.join("run");
    let cfg = d.join("fastsuite.toml");
    // fastcheck only: no step trims, so there is nothing to guard.
    write_config(&cfg, &out, &raw, "\"fastcheck\"", "");

    let (stdout, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
    ]);

    assert_eq!(code, 0, "{stderr}");
    assert!(!stderr.contains("adapter guard"), "{stderr}");
    assert!(!stderr.contains("adapter survived"), "{stderr}");
    assert!(!stdout.contains("adapter guard"), "{stdout}");
    // No guard scratch files were produced.
    assert!(!out.join("logs/s1.adapter-guard.raw.json").is_file());
}

/// The guard is opt-out: --no-adapter-guard silences even a bad trim.
#[test]
fn the_guard_can_be_disabled() {
    require_tools(&["fastcheck", "fastprep"]);
    let d = tmpdir("disabled");
    let raw = d.join("s1_R1.fastq");
    fixture_fastq(&raw);
    let out = d.join("run");
    let cfg = d.join("fastsuite.toml");
    write_config(&cfg, &out, &raw, "\"fastprep\"", "\"--no-adapter\"");

    let (_, stderr, code) = run(&[
        "run",
        "-c",
        cfg.to_str().unwrap(),
        "--tools-dir",
        tools_dir().to_str().unwrap(),
        "--no-adapter-guard",
    ]);

    assert_eq!(code, 0, "{stderr}");
    assert!(!stderr.contains("adapter"), "opted out: {stderr}");
    let man = read_manifest(&out);
    assert_eq!(man["adapter_guard"]["enabled"], false);
}

// ---------------------------------------------------------------------------
// The standalone `fasthelper guard` verb (Task 1.7).
// ---------------------------------------------------------------------------

/// A FASTQ record with a fixed length.
fn record(id: &str, seq: &str) -> String {
    format!("@{id}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
}

fn write_fastq(path: &Path, body: &str) {
    std::fs::write(path, body).unwrap();
}

#[test]
fn guard_trips_on_a_read_that_is_all_adapter() {
    require_tools(&["fastcheck"]);
    let d = tmpdir("trip");
    let f = d.join("adapter.fq");
    // 40 reads of 150 bases, each entirely the TruSeq sequence repeated.
    let mut body = String::new();
    for k in 0..40 {
        body.push_str(&record(&format!("a{k}"), &(TRUSEQ_R1.repeat(5))[..150]));
    }
    write_fastq(&f, &body);

    let out = Command::new(fasthelper())
        .args([
            "guard",
            "-1",
            f.to_str().unwrap(),
            "--adapter-preset",
            "illumina",
            "--fail-on-adapter-residual",
        ])
        .output()
        .expect("spawn fasthelper");
    assert_eq!(
        out.status.code(),
        Some(2),
        "an all-adapter read must trip the guard (exit 2); stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("RESIDUAL TOO HIGH"),
        "the message must name the failure: {err}"
    );
}

#[test]
fn guard_passes_on_a_clean_read_and_reports_the_residual() {
    require_tools(&["fastcheck"]);
    let d = tmpdir("clean");
    let f = d.join("clean.fq");
    let mut body = String::new();
    for k in 0..400 {
        body.push_str(&record(&format!("c{k}"), &random_dna(150, k)));
    }
    write_fastq(&f, &body);

    let out = Command::new(fasthelper())
        .args([
            "guard",
            "-1",
            f.to_str().unwrap(),
            "--adapter-preset",
            "illumina",
            "--fail-on-adapter-residual",
        ])
        .output()
        .expect("spawn fasthelper");
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(0),
        "a clean read must pass; stderr: {err}"
    );
    assert!(
        err.contains("-- OK"),
        "the report must be one readable line: {err}"
    );
    assert!(
        err.contains("% of 400 reads still carry adapter"),
        "the report must name the measured residual and the read count: {err}"
    );
}

#[test]
fn guard_exits_2_only_when_fail_is_requested() {
    require_tools(&["fastcheck"]);
    let d = tmpdir("nofail");
    let f = d.join("adapter.fq");
    let mut body = String::new();
    for k in 0..40 {
        body.push_str(&record(&format!("a{k}"), &(TRUSEQ_R1.repeat(5))[..150]));
    }
    write_fastq(&f, &body);

    let out = Command::new(fasthelper())
        .args([
            "guard",
            "-1",
            f.to_str().unwrap(),
            "--adapter-preset",
            "illumina",
        ])
        .output()
        .expect("spawn fasthelper");
    assert_eq!(
        out.status.code(),
        Some(0),
        "without --fail-on-adapter-residual a trip is a warning, not a failure"
    );
}

/// Deterministic pseudo-random ACGT of length `n`.
fn random_dna(n: usize, seed: u64) -> String {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            match x % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        })
        .collect()
}

const TRUSEQ_R1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
