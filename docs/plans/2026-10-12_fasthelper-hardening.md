# fasthelper: hardening (compare fix, fastcheck step, parallel, resume, check, install, thresholds)

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Turn `fasthelper` from a working single-threaded runner into the hardened orchestrator the plan promises: a compare that creates its output dir, fastcheck wired as a first-class run step, config validation before hour three, `-P` parallel samples, resumable runs, one-command installs with legacy names, and `--fail-on` thresholds with exit code 2.

**Architecture:** `fasthelper` stays the MultiQC slot (D-5): it owns no QC logic and only orchestrates the per-tool binaries plus aggregates `fastsuite.report.v1` JSON. New work is (a) argv wiring for `fastcheck` in the runner, (b) a `check` verb sharing the runner's resolution logic, (c) a sample worker pool around the existing per-sample loop, (d) a manifest keyed on argv+tool version for resume, (e) a symlink installer, (f) a threshold parser applied to the aggregated matrix. Every step keeps the runner's per-step failure semantics (fail loudly, name the log, exit non-zero).

**Tech Stack:** Rust, clap 4 (derive), serde/toml/serde_json, anyhow; tests are integration tests in `tools/fasthelper/tests/` that shell out to the built binary via `env!("CARGO_BIN_EXE_fasthelper")`.

**Baseline (verified, 2026-10-12):** `cargo test --offline` in `tools/fasthelper` → 13 passed (adapters 5, cli 7, compare_test 1), 0 failed. `fasthelper doctor --tools-dir tools` lists fastcheck, fastdedup, fasthelper, fastmaps, fastprep, fastsim. `fasthelper run -c examples/fastsuite.toml --tools-dir tools` runs the sim control end to end.

> **Already done on disk (do not re-implement):** the working tree of `tools/fasthelper` already contains the compare `create_dir_all` fix (`src/compare.rs:31`), the regression test `compare_creates_a_missing_output_directory` (`tests/cli.rs:36`), and the `fastcheck` entry in `ToolSet`/doctor/manifest (`src/runner.rs:35,81,427,487,530`). **Tasks 1–2 below verify that state and close the two gaps those commits left** (fastcheck as a *run step* and the example config). If a subagent starts on a clean checkout of the last commit, run Tasks 1–2 first; if the tree already has them, skip straight to Task 3.

---

## Task 1: Verify the compare mkdir fix + regression test

**Objective:** Confirm the compare fix and its regression test are present and green (they are, on disk — this task exists so a fresh checkout gets the fix before anything else builds on compare).

**Files:**
- Verify: `tools/fasthelper/src/compare.rs` (create_dir_all at the top of `cmd_compare`)
- Verify: `tools/fasthelper/tests/cli.rs:36-61` (`compare_creates_a_missing_output_directory`)

**Step 1: Inspect the code**

Run: `grep -n "create_dir_all" tools/fasthelper/src/compare.rs`
Expected: `31:    std::fs::create_dir_all(outdir)` inside `pub fn cmd_compare`.

**Step 2: Run the regression test**

Run: `cd tools/fasthelper && cargo test --offline compare_creates_a_missing_output_directory`
Expected: `test compare_creates_a_missing_output_directory ... ok`

**Step 3: Run the whole fasthelper suite**

Run: `cargo test --offline`
Expected: `13 passed; 0 failed` across `tests/adapters.rs`, `tests/cli.rs`, `tests/compare_test.rs`.

**Step 4: Commit if the fix was missing (clean checkout only)**

```bash
git add tools/fasthelper/src/compare.rs
git commit -m "fix(fasthelper): create the output dir in compare, like concatenate"
```

---

## Task 2: fastcheck as a run step (argv wiring) + manifest + example config

**Objective:** A `fastcheck` entry in `steps` runs `fastcheck <current-fastqs> --json <work>/fastcheck.json`, the manifest records the resolved binary, and the example config documents the step.

**Design decision (per the verified facts):** fastcheck takes **positional FILE...** arguments, never `-1/-2/-o`. It never writes FASTQ output, so it is a *terminal* step: after fastcheck, the per-sample chain has no further consumer and the runner must not treat its "output" as an input to anything. The existing loop deletes the previous step's intermediate after its last consumer; fastcheck reads the current inputs and produces only reports, so the deletion logic is unchanged (the current files still have no *later* consumer and get deleted as today). The per-position TSV is emitted only when the config passes it via `[steps.fastcheck] args`.

**Files:**
- Modify: `tools/fasthelper/src/runner.rs:320-350` (the per-step `(out1,out2,json_path,maps_dir)` match) and `runner.rs:154-207` (`argv_for`)
- Modify: `tools/fasthelper/src/config.rs:40-42` (doc comment: supported steps)
- Modify: `examples/fastsuite.toml` (fastcheck step block)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tools/fasthelper/tests/cli.rs`:

```rust
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
        let argv: Vec<String> =
            serde_json::from_value(cmd["argv"].clone()).unwrap();
        assert!(argv.contains(&fq.display().to_string()), "argv: {argv:?}");
        assert!(argv.contains(&"--json".to_string()), "argv: {argv:?}");
        assert!(out.join("s1/fastcheck.json").is_file());
    }
}
```

Note for the implementer: `tmpdir`/`run`/`write_report` helpers already exist at the top of `tests/cli.rs`; `option_env!` is replaced at compile time — prefer `env!("CARGO_MANIFEST_DIR")`-relative `../..` instead if unsure. The test passes on machines where fastcheck is not on `$PATH` only via the `code != 0` branch; when fastcheck *is* built (the suite machines do), it exercises the full path. Simpler alternative if the `option_env!` indirection bothers you: resolve the binary from `CARGO_MANIFEST_DIR` directly — see Step 3 for the runner-side change that makes the step runnable without a tools_dir.

**Step 2: Run the test to verify it fails**

Run: `cd tools/fasthelper && cargo test --offline fastcheck_step_runs_and_writes_a_v1_report`
Expected: FAIL — the run prints `fasthelper: unsupported step 'fastcheck' — skipped` and the assertion fires. (Also verify the current behaviour first: `cargo test` alone still shows `13 passed` — the new test is the only failure.)

**Step 3: Implement**

In `tools/fasthelper/src/runner.rs`, `argv_for`, add the fastcheck arm before the `_ => {}` fallback:

```rust
        "fastcheck" => {
            // fastcheck takes positional FASTQ files, never -1/-2/-o; it
            // writes only reports. --json is mandatory so the run directory
            // always carries a fastsuite.report.v1 for the aggregator.
            a.push(s(in1));
            if let Some(r2) = in2 {
                a.push(s(r2));
            }
            a.push("--json".into());
            a.push(s(json));
        }
```

In the per-step match inside `run` (currently at `runner.rs:324-349`), add the fastcheck arm before the `other =>` rejection:

```rust
                "fastcheck" => (
                    PathBuf::new(), // no FASTQ output
                    None,
                    work.join("fastcheck.json"),
                    work.join("maps"), // unused; the aggregator reads the JSON
                ),
```

`ToolSet::get("fastcheck")` already returns `self.check` (commit 6748839), so the `let Some(bin) = ts.get(tool)` lookup resolves it with no further change. The manifest already records `"fastcheck": ts.check...` (runner.rs:427).

**Step 4: Run the test to verify pass**

Run: `cd tools/fasthelper && cargo test --offline`
Expected: `14 passed; 0 failed`.

**Step 5: Update the config doc + example**

`tools/fasthelper/src/config.rs:40` — change the comment:

```rust
    /// Steps to execute, in order. Supported: fastdedup, fastprep, fastmaps, fastcheck.
```

`examples/fastsuite.toml` — add after the `[steps.fastdedup]` block:

```toml
# fastcheck is the terminal QC step: it reads the chain's current FASTQs and
# writes fastsuite.report.v1 + a fastsim profile, never a new FASTQ. It has no
# output to feed a later step, so list it last. Extra flags go here too.
# [steps.fastcheck]
# args = ["--insert-hist", "--max-read-len", "300"]
```

**Step 6: Run the suite and the doctor**

Run: `cd tools/fasthelper && cargo test --offline && cargo run --offline -- doctor --tools-dir ../..`
Expected: `14 passed` and a table whose first row is `fastcheck tools/fastcheck/target/release/fastcheck fastcheck 0.1.0`.

**Step 7: Commit**

```bash
git add tools/fasthelper/src/runner.rs tools/fasthelper/src/config.rs \
        examples/fastsuite.toml tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): run fastcheck as a terminal QC step"
```

---

## Task 3: `fasthelper check -c config.toml` — config validation

**Objective:** A `check` verb that validates a config *before* a long run starts: input files exist, requested tools resolve, steps are known, and the per-step args are flags the tools actually accept (via a `--help` scrape).

**Design:** validation lives in a new `src/check.rs`; it reuses `runner::resolve`/`runner::toolset` for binaries and `Config::load` for parsing. Flag validation is deliberately simple: run `<bin> --help` once, capture stdout+stderr, and for each `--flag` in a step's args assert the string `--flag` (or its long form without the value) appears. Budget flags like `fastcheck --reads N` are known-good; a typo like `--reds` fails here instead of after hour three.

**Files:**
- Create: `tools/fasthelper/src/check.rs`
- Modify: `tools/fasthelper/src/cli.rs` (new `Check` subcommand + args struct)
- Modify: `tools/fasthelper/src/main.rs` (dispatch)
- Modify: `tools/fasthelper/src/lib.rs` (`pub mod check;`)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tests/cli.rs`:

```rust
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
        "check", "-c", ok_cfg.to_str().unwrap(),
        "--tools-dir", std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".into()),
    ]);
    assert!(code2 == 0 || stderr2.contains("not found"), "{code2} {stderr2}");
}
```

`tests/cli.rs` already imports `serde_json`; nothing new needed.

**Step 2: Run the test to verify failure**

Run: `cd tools/fasthelper && cargo test --offline check_rejects_unknown_step_and_missing_input`
Expected: FAIL — clap rejects the unknown `check` subcommand: `error: unrecognized subcommand 'check'` (exit 2).

**Step 3: Implement `src/check.rs`**

Create `tools/fasthelper/src/check.rs`:

```rust
//! `fasthelper check -c config.toml` — validate a run config before it starts.
//!
//! The point is to fail at minute zero, not at hour three: every [[samples]]
//! input exists, every requested binary resolves, every step name is known,
//! and every flag passed under `[steps.<tool>] args` is accepted by that
//! tool's `--help` text.

use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::cli::CheckArgs;
use crate::config::Config;
use crate::runner;

const KNOWN_STEPS: &[&str] = &["fastdedup", "fastprep", "fastmaps", "fastcheck"];

/// Run `<bin> --help` and return the combined output text.
fn help_text(bin: &Path) -> Option<String> {
    let out = Command::new(bin).arg("--help").output().ok()?;
    let mut s = String::from_utf8_lossy(&out.stdout).to_string();
    s.push_str(&String::from_utf8_lossy(&out.stderr));
    Some(s)
}

/// Is `flag` (a token starting with `--`) accepted by the tool's help text?
fn flag_accepted(help: &str, flag: &str) -> bool {
    // --insert-hist and --insert-hist=... both count; --foo bar never does.
    help.contains(flag)
}

pub fn cmd_check(args: &CheckArgs) -> Result<()> {
    let cfg = Config::load(&args.config)?;
    let config_dir = args
        .config
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));
    let tools_dir = cfg.tools_dir(args.tools_dir.clone());
    let ts = runner::toolset(tools_dir);

    let mut errors: Vec<String> = Vec::new();

    // 1. Step names.
    for step in &cfg.run.steps {
        if !KNOWN_STEPS.contains(&step.as_str()) {
            errors.push(format!(
                "unknown step '{step}' — supported: {}",
                KNOWN_STEPS.join(", ")
            ));
        }
        if let Some(extra) = cfg.step(step).args.iter().find(|a| a.starts_with("--")) {
            if cfg.step(step).enabled {
                if let Some(bin) = ts.get(step) {
                    if let Some(h) = help_text(bin) {
                        if !flag_accepted(&h, extra) {
                            errors.push(format!(
                                "{step} does not accept '{extra}' (see {})",
                                bin.display()
                            ));
                        }
                    }
                }
            }
        }
    }

    // 2. Sample inputs.
    for s in &cfg.samples {
        for f in [Some(&s.r1), s.r2.as_ref()].into_iter().flatten() {
            let p = if f.is_absolute() {
                f.clone()
            } else {
                config_dir.join(f)
            };
            if !p.is_file() {
                errors.push(format!("sample {}: input {} not found", s.name, p.display()));
            }
        }
    }

    // 3. Tools.
    for step in &cfg.run.steps {
        if KNOWN_STEPS.contains(&step.as_str()) && cfg.step(step).enabled && ts.get(step).is_none()
        {
            errors.push(format!("tool {step} not found in the tools dir or $PATH"));
        }
    }

    if errors.is_empty() {
        println!(
            "fasthelper check: OK — {} sample(s), {} step(s), tools resolved",
            cfg.samples.len(),
            cfg.run.steps.len()
        );
        return Ok(());
    }
    for e in &errors {
        eprintln!("fasthelper check: {e}");
    }
    bail!("{} problem(s) — fix them before running", errors.len())
}
```

**Step 4: CLI plumbing**

In `tools/fasthelper/src/cli.rs`, add to the `Cmd` enum (after `Doctor`):

```rust
    /// Validate a run config: inputs exist, tools resolve, steps and flags are known.
    Check(CheckArgs),
```

and the args struct at the bottom of the file:

```rust
#[derive(Args, Debug)]
pub struct CheckArgs {
    /// Path to the run config (TOML).
    #[arg(short = 'c', long)]
    pub config: PathBuf,
    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,
}
```

`runner::toolset` and `ToolSet::get` are already `pub` — no visibility change needed.

In `tools/fasthelper/src/main.rs`, add the dispatch arm:

```rust
        cli::Cmd::Check(a) => fasthelper::check::cmd_check(&a),
```

In `tools/fasthelper/src/lib.rs`, add `pub mod check;`.

**Step 5: Run the test to verify pass**

Run: `cd tools/fasthelper && cargo test --offline`
Expected: `15 passed; 0 failed`.

**Step 6: Manual verification**

```bash
printf '[run]\nsteps = ["fastnope"]\n\n[[samples]]\nname = "s1"\nr1 = "nope_R1.fastq"\n' > /tmp/bad.toml
./target/debug/fasthelper check -c /tmp/bad.toml; echo "exit=$?"
```
Expected: two `fasthelper check:` lines (unknown step; missing input) then `Error: 2 problem(s) — fix them before running`, exit 1.

And the good path, from the repo root:
```bash
./tools/fasthelper/target/debug/fasthelper check -c examples/fastsuite.toml --tools-dir tools
```
Expected: `fasthelper check: OK — 1 sample(s), 3 step(s), tools resolved` (the demo config's single sample is the simulated control, which has no input file to check).

**Step 7: Commit**

```bash
git add tools/fasthelper/src/check.rs tools/fasthelper/src/cli.rs \
        tools/fasthelper/src/main.rs tools/fasthelper/src/lib.rs \
        tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): check -c validates a config before the run starts"
```

---

## Task 4: `--fail-on <expr>` thresholds + exit code 2

**Objective:** `fasthelper concatenate --fail-on 'polyG_frac>0.05'` (repeatable) exits 2 — "ran fine, but a threshold tripped" per CLI-SPEC §5 — with the offending metrics named on stderr, while still writing all report formats.

**Design:** parse each expression as `metric <op> value` where op ∈ `> < >= <= == !=`; evaluate against the aggregated matrix (last row wins when a metric appears in several reports); collect violations; print them; `std::process::exit(2)` after the reports are written. Exit 2 must come from `main` (returning `Err` gives exit 1).

**Files:**
- Modify: `tools/fasthelper/src/report.rs` (threshold parsing + evaluation)
- Modify: `tools/fasthelper/src/cli.rs` (`--fail-on` on `ReportArgs`)
- Modify: `tools/fasthelper/src/main.rs` (exit 2)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tests/cli.rs`:

```rust
#[test]
fn fail_on_trips_exit_2_and_names_the_metric() {
    let d = tmpdir("failon");
    let reports = d.join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    write_report(&reports, "sampleA.json", 100, 0.41); // gc_fraction 0.41 < 0.5

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
```

(`Path` is already imported in `tests/cli.rs`; the default outdir is `.`, so the second assertion tolerates the cwd write.)

**Step 2: Run the test to verify failure**

Run: `cd tools/fasthelper && cargo test --offline fail_on_trips_exit_2_and_names_the_metric`
Expected: FAIL — clap: `error: unexpected argument '--fail-on' found` (exit 2 of clap itself, but the test asserts the fasthelper-level code 2; either way the test does not pass today).

**Step 3: Implement the threshold engine**

In `tools/fasthelper/src/report.rs`, add near the bottom (before `cmd_concatenate`):

```rust
/// One `metric OP value` threshold, e.g. `polyG_frac>0.05`.
#[derive(Debug, Clone)]
pub struct Threshold {
    pub metric: String,
    pub op: char,
    pub value: f64,
}

impl Threshold {
    /// Parse `metric>0.05` / `reads<=100` / ... Ops: > < >= <= == !=.
    pub fn parse(expr: &str) -> Result<Self> {
        let ops = [">=", "<=", "==", "!=", ">", "<"];
        for op in ops {
            if let Some(i) = expr.find(op) {
                let metric = expr[..i].trim().to_string();
                let value: f64 = expr[i + op.len()..]
                    .trim()
                    .parse()
                    .with_context(|| format!("parsing threshold '{expr}': value is not a number"))?;
                if metric.is_empty() {
                    anyhow::bail!("threshold '{expr}' has no metric name");
                }
                return Ok(Threshold {
                    metric,
                    op: op.chars().next().unwrap(),
                    value,
                });
            }
        }
        anyhow::bail!(
            "threshold '{expr}' is not METRIC OP VALUE (ops: > < >= <= == !=)"
        );
    }

    pub fn tripped(&self, actual: f64) -> bool {
        match self.op {
            '>' => actual > self.value,
            '<' => actual < self.value,
            _ if self.op == '=' => (actual - self.value).abs() < 1e-12,
            _ if self.op == '!' => (actual - self.value).abs() >= 1e-12,
            _ => {
                // two-char ops collapse to their first char here
                match self.op {
                    '>' => actual >= self.value,
                    '<' => actual <= self.value,
                    _ => false,
                }
            }
        }
    }
}

/// Evaluate thresholds against the aggregated metrics (last value wins).
/// Returns the tripped descriptions, e.g. `polyG_frac=0.12 > 0.05`.
pub fn evaluate_thresholds(
    reports: &[ToolReport],
    thresholds: &[Threshold],
) -> Vec<String> {
    let mut values: BTreeMap<String, f64> = BTreeMap::new();
    for r in reports {
        for (k, v) in &r.metrics {
            values.insert(k.clone(), *v);
        }
    }
    let mut out = Vec::new();
    for t in thresholds {
        if let Some(actual) = values.get(&t.metric) {
            if t.tripped(*actual) {
                out.push(format!("{}={} {} {}", t.metric, actual, t.op, t.value));
            }
        }
    }
    out
}
```

Note: `>=`/`<=`/`==`/`!=` all parse to the first char (`>`, `<`, `=`, `!`) — the `tripped` match handles the two-char forms by re-interpreting those chars. This is intentionally compact; if a reviewer prefers explicit enums, refactor in the code-quality pass, not here.

In `cmd_concatenate` (`report.rs:654-688`), after `write_reports` and before the summary println, add:

```rust
    if !a.fail_on.is_empty() {
        let thresholds: Vec<Threshold> = a
            .fail_on
            .iter()
            .map(|s| Threshold::parse(s))
            .collect::<Result<Vec<_>>>()?;
        let tripped = evaluate_thresholds(&reports, &thresholds);
        if !tripped.is_empty() {
            for t in &tripped {
                eprintln!("fasthelper: threshold tripped: {t}");
            }
            eprintln!(
                "fasthelper: {} threshold(s) tripped — report written, exiting 2",
                tripped.len()
            );
            std::process::exit(2);
        }
    }
```

(Place it after `write_reports` so the report always exists on disk before the exit; parse errors still return `Err` → exit 1.)

**Step 4: CLI plumbing**

In `tools/fasthelper/src/cli.rs`, add to `ReportArgs`:

```rust
    /// Fail (exit 2) when a metric breaks the rule, e.g. --fail-on polyG_frac>0.05.
    /// Repeatable. Exit 2 = "ran fine, but a threshold tripped" (CLI-SPEC §5).
    #[arg(long = "fail-on", value_name = "METRIC>OP>VALUE")]
    pub fail_on: Vec<String>,
```

No `main.rs` change is needed: `cmd_concatenate` itself calls `std::process::exit(2)` after the reports are written.

**Step 5: Run the test to verify pass**

Run: `cd tools/fasthelper && cargo test --offline`
Expected: `16 passed; 0 failed`.

**Step 6: Manual verification**

```bash
mkdir -p /tmp/fo && cp /tmp/demo_b/fastcheck.json /tmp/fo/E94/fastcheck.json 2>/dev/null || { mkdir -p /tmp/fo/E94 && cp /tmp/demo_b/fastcheck.json /tmp/fo/E94/fastcheck.json; }
./tools/fasthelper/target/debug/fasthelper concatenate /tmp/fo -o /tmp/fo/out --fail-on 'gc_fraction>0.9'
echo "exit=$?"
./tools/fasthelper/target/debug/fasthelper concatenate /tmp/fo -o /tmp/fo/out --fail-on 'gc_fraction>0.01'
echo "exit=$?"
```
Expected (first): `fasthelper: threshold tripped: gc_fraction=0.4986 > 0.9` … `exit=2`. (second): no tripped line, `exit=0`.

**Step 7: Commit**

```bash
git add tools/fasthelper/src/report.rs tools/fasthelper/src/cli.rs \
        tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): --fail-on thresholds with exit code 2"
```

---

## Task 5: `-P/--jobs` parallel samples

**Objective:** `fasthelper run -P 8` runs up to 8 samples concurrently (per-sample step order stays strict; threads are divided across the pool: sample i gets `threads/P` threads, min 1).

**Design:** extract the body of the `for sample in &samples` loop into `run_sample(...) -> Result<Vec<serde_json::Value>>` (its manifest entries) plus a shared failure list. The driver runs samples through a semaphore-bounded pool: `std::thread::scope` with a channel, or a simple chunked approach — simplest correct version: spawn min(P, samples.len()) threads, each pulling from a shared `Mutex<VecDeque<Sample>>` work queue. Merge manifest entries sorted by (sample, step index) so the manifest stays deterministic.

**Files:**
- Modify: `tools/fasthelper/src/runner.rs` (extract `run_sample`, add pool, thread division)
- Modify: `tools/fasthelper/src/cli.rs` (`-P/--jobs` on `RunArgs`)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tests/cli.rs`:

```rust
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
    assert!(!stderr.contains("unexpected argument"), "-P must exist: {stderr}");
    if code == 0 {
        let man: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(out.join("manifest.json")).unwrap())
                .unwrap();
        let cmds = man["commands"].as_array().unwrap();
        assert_eq!(cmds.len(), 2, "both samples must run: {cmds:?}");
    }
}
```

**Step 2: Run the test to verify failure**

Run: `cd tools/fasthelper && cargo test --offline jobs_flag_is_accepted_and_the_run_still_succeeds`
Expected: FAIL — `error: unexpected argument '-P' found`.

**Step 3: Implement**

In `tools/fasthelper/src/cli.rs`, add to `RunArgs`:

```rust
    /// Samples to process in parallel (default 1). Threads are divided across the pool.
    #[arg(short = 'P', long, default_value_t = 1)]
    pub jobs: usize,
```

In `tools/fasthelper/src/runner.rs`:

1. Extract everything from `let work = outdir.join(&sample.name);` through the `to_delete` loop into:

```rust
fn run_sample(
    sample: &Sample,
    cfg: &Config,
    ts: &ToolSet,
    outdir: &Path,
    logs: &Path,
    config_dir: &Path,
    dry_run: bool,
    threads: usize,
) -> Result<Vec<serde_json::Value>> {
    let work = outdir.join(&sample.name);
    std::fs::create_dir_all(&work)?;
    // ... the existing body, minus the `failed` pushes: return a Vec of
    // manifest entries; the caller records failures from the entries'
    // exit codes.
    Ok(entries)
}
```

2. Replace the `for sample in &samples { ... }` loop with the pool:

```rust
    let samples = std::sync::Arc::new(samples);
    let queue = std::sync::Mutex::new(std::collections::VecDeque::from(samples.iter().cloned().collect::<Vec<_>>()));
    let per_sample_threads = (cfg.run.threads / args.jobs).max(1);
    let entries: std::sync::Mutex<Vec<serde_json::Value>> = std::sync::Mutex::new(Vec::new());
    let failures: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());

    std::thread::scope(|scope| {
        for _ in 0..args.jobs.max(1) {
            scope.spawn(|| loop {
                let sample = { queue.lock().unwrap().pop_front() };
                let Some(sample) = sample else { break };
                match run_sample(
                    &sample,
                    &cfg,
                    &ts,
                    &outdir,
                    &logs,
                    &config_dir,
                    args.dry_run,
                    per_sample_threads,
                ) {
                    Ok(mut e) => entries.lock().unwrap().append(&mut e),
                    Err(err) => {
                        eprintln!("fasthelper: {} failed: {err}", sample.name);
                        failures.lock().unwrap().push(sample.name.clone());
                    }
                }
            });
        }
    });

    let mut manifest = entries.into_inner().unwrap();
    manifest.sort_by_key(|m| {
        (
            m["sample"].as_str().unwrap_or("").to_string(),
            m["tool"].as_str().unwrap_or("").to_string(),
        )
    });
    let failed = failures.into_inner().unwrap();
```

The `cfg.run.threads` passed to each step becomes `per_sample_threads` — that is the only behavioural change for `-P 1` (threads/jobs = threads/1, unchanged).

**Step 4: Run the test to verify pass**

Run: `cd tools/fasthelper && cargo test --offline`
Expected: `17 passed; 0 failed`.

**Step 5: Manual verification**

```bash
./tools/fasthelper/target/debug/fasthelper run -c examples/fastsuite.toml --tools-dir tools -P 4 --dry-run
```
Expected: `[dry-run] .../fastcheck ...` line per step, no errors, exit 0.

**Step 6: Commit**

```bash
git add tools/fasthelper/src/runner.rs tools/fasthelper/src/cli.rs \
        tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): -P/--jobs parallel samples"
```

---

## Task 6: resume + `--from <step>`

**Objective:** A 34 TB run is restartable. When `manifest.json` exists, a step is skipped when (a) its recorded argv matches the argv the run would build today and (b) the recorded tool binary's `--version` matches the currently-resolved one; `--from <step>` re-runs from that step onward regardless.

**Design:** read the old manifest first (`outdir/manifest.json`); build a map `(sample, tool) -> (argv, binary, version)`. Before running a step, if not in `--from` mode and the map has a matching entry with identical argv and binary path and a recorded version, print `resumed <sample>/<tool>` and skip. Steps after a `--from` step run unconditionally (mark `from_step_seen` once the named step is reached; everything from there on is forced).

**Files:**
- Modify: `tools/fasthelper/src/runner.rs` (resume state + skip logic in `run_sample`)
- Modify: `tools/fasthelper/src/cli.rs` (`--from` on `RunArgs`)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tests/cli.rs`:

```rust
#[test]
fn from_flag_forces_a_step_rerun_and_resume_skips_matches() {
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
    // First run populates the manifest.
    let (_, _, c1) = run(&["run", "-c", cfg.to_str().unwrap()]);
    if c1 != 0 {
        return; // fastcheck not installed here; the flag test below still runs
    }
    // Second run: identical config -> the step is resumed (skipped).
    let (stdout2, _, c2) = run(&["run", "-c", cfg.to_str().unwrap()]);
    assert_eq!(c2, 0);
    assert!(stdout2.contains("resumed s1/fastcheck"), "{stdout2}");
    // --from fastcheck forces it to run again.
    let (stdout3, _, c3) = run(&["run", "-c", cfg.to_str().unwrap(), "--from", "fastcheck"]);
    assert_eq!(c3, 0);
    assert!(!stdout3.contains("resumed s1/fastcheck"), "{stdout3}");
}
```

**Step 2: Run the test to verify failure**

Run: `cd tools/fasthelper && cargo test --offline from_flag_forces_a_step_rerun_and_resume_skips_matches`
Expected: FAIL — clap rejects `--from`.

**Step 3: Implement**

In `cli.rs`, add to `RunArgs`:

```rust
    /// Re-run starting at this step (earlier steps are resumed/skipped as usual).
    #[arg(long)]
    pub from: Option<String>,
```

In `runner.rs`, inside `run_sample`, before building the argv for a step:

```rust
            let forced = args_from.map(|s| s == tool).unwrap_or(false);
            let key = format!("{}/{}", sample.name, tool);
            if !forced {
                if let Some((prev_argv, prev_bin, prev_ver)) = resume.get(&key) {
                    let bin_now = s(bin);
                    let ver_now = tool_version(bin);
                    if prev_argv == &argv_placeholder && prev_bin == &bin_now && prev_ver == ver_now
                    {
                        println!("fasthelper: resumed {key}");
                        continue;
                    }
                }
            }
```

Concretely: load `resume: BTreeMap<String, (Vec<String>, String, Option<String>)>` from the previous manifest at the top of `run` (`commands` entries carry `sample`, `tool`, `argv`, `binary`; version is re-resolved live, so store the *previous run's* version alongside), pass `&resume` and `args.from` into `run_sample`, and compare against the freshly built argv. The comparison needs the argv *before* running — so build the argv first, then consult `resume`, then `run_step`. That is the natural order already; only the map lookup is new.

Version capture (reuse the doctor pattern):

```rust
fn tool_version(bin: &Path) -> Option<String> {
    std::process::Command::new(bin)
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
}
```

**Step 4: Run the test to verify pass**

Run: `cd tools/fasthelper && cargo test --offline`
Expected: `18 passed; 0 failed`.

**Step 5: Manual verification**

```bash
./tools/fasthelper/target/debug/fasthelper run -c examples/fastsuite.toml --tools-dir tools
./tools/fasthelper/target/debug/fasthelper run -c examples/fastsuite.toml --tools-dir tools   # expect "resumed sim_E/fastcheck" etc.
./tools/fasthelper/target/debug/fasthelper run -c examples/fastsuite.toml --tools-dir tools --from fastprep   # fastdedup resumed, fastprep re-runs
```
Expected: second run prints `resumed ...` for every completed step and rewrites identical manifest content (same argv, so no spurious diffs); third run re-runs from fastprep.

**Step 6: Commit**

```bash
git add tools/fasthelper/src/runner.rs tools/fasthelper/src/cli.rs \
        tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): resume completed steps and --from <step>"
```

---

## Task 7: `fasthelper install --prefix DIR` with legacy names (D-8)

**Objective:** `fasthelper install --prefix ~/.local` creates `~/.local/bin/` and symlinks every resolved suite binary there, plus the legacy names (`rpolyx → fastprep`, `fdedup → fastdedup`, `novaseq-sim → fastsim`) unless `--no-legacy`.

**Design:** resolve each tool via `runner::toolset(tools_dir)`; for each found binary, `std::os::unix::fs::symlink` into `<prefix>/bin/` (create it); symlink legacy names to the same targets. `--no-legacy` skips the legacy set. Already-installed symlinks are replaced silently (symlink_metadata tells us it is a link; remove then re-create). Non-Unix: bail with a clear message (the suite is linux-64 per `pixi.toml`).

**Files:**
- Create: `tools/fasthelper/src/install.rs`
- Modify: `tools/fasthelper/src/cli.rs` (`Install` subcommand + args)
- Modify: `tools/fasthelper/src/main.rs` (dispatch)
- Modify: `tools/fasthelper/src/lib.rs` (`pub mod install;`)
- Test: `tools/fasthelper/tests/cli.rs`

**Step 1: Write the failing test**

Append to `tests/cli.rs`:

```rust
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
        std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".into()),
    ]);
    if code != 0 {
        // No binaries resolved here (empty tools dir): must still say so, not crash.
        assert!(stderr.contains("no FastSuite binaries") || stderr.contains("not found"), "{stderr}");
        return;
    }
    let bin = prefix.join("bin");
    assert!(bin.is_dir(), "{stdout} {stderr}");
    // Whatever was found, the directory must contain at least one symlink.
    let n = std::fs::read_dir(&bin).unwrap().count();
    assert!(n > 0, "nothing installed");
}
```

**Step 2: Run the test to verify failure**

Run: `cd tools/fasthelper && cargo test --offline install_symlinks_binaries_and_legacy_names`
Expected: FAIL — clap: unrecognized subcommand `install`.

**Step 3: Implement `src/install.rs`**

Create `tools/fasthelper/src/install.rs`:

```rust
//! `fasthelper install --prefix DIR` — one bin/ with every suite binary.
//!
//! Legacy binary names (rpolyx, fdedup, novaseq-sim) are installed as
//! symlinks alongside the canonical ones, per docs/DECISIONS.md D-8: cargo and
//! bioconda never install a binary name, fasthelper does.

use anyhow::{bail, Context, Result};
use std::os::unix::fs::symlink;
use std::path::{Path, PathBuf};

use crate::cli::InstallArgs;
use crate::runner::{self, ToolSet};

/// canonical name -> legacy symlink name (D-8)
const LEGACY: &[(&str, &str)] = &[
    ("fastprep", "rpolyx"),
    ("fastdedup", "fdedup"),
    ("fastsim", "novaseq-sim"),
];

fn link(dst: &Path, target: &Path) -> Result<()> {
    if let Ok(md) = std::fs::symlink_metadata(dst) {
        if md.file_type().is_symlink() {
            std::fs::remove_file(dst).ok();
        }
    }
    symlink(target, dst)
        .with_context(|| format!("linking {} -> {}", dst.display(), target.display()))
}

pub fn cmd_install(args: &InstallArgs) -> Result<()> {
    if cfg!(not(unix)) {
        bail!("install is only implemented for unix platforms");
    }
    let ts: ToolSet = runner::toolset(args.tools_dir.clone());
    let mut targets: Vec<(String, PathBuf)> = Vec::new();
    for tool in ["fastcheck", "fastdedup", "fastprep", "fastsim", "fastmaps", "fasthelper"] {
        let bin = if tool == "fasthelper" {
            std::env::current_exe().ok()
        } else {
            ts.get(tool).cloned()
        };
        if let Some(b) = bin {
            targets.push((tool.to_string(), b));
        }
    }
    if targets.is_empty() {
        bail!("no FastSuite binaries found — build them first or pass --tools-dir");
    }
    let bin = args.prefix.join("bin");
    std::fs::create_dir_all(&bin)
        .with_context(|| format!("creating {}", bin.display()))?;
    let mut made = Vec::new();
    for (name, target) in &targets {
        link(&bin.join(name), target)?;
        made.push(format!("{name} -> {}", target.display()));
    }
    if !args.no_legacy {
        for (canonical, legacy) in LEGACY {
            if let Some((_, target)) = targets.iter().find(|(n, _)| n == canonical) {
                link(&bin.join(legacy), target)?;
                made.push(format!("{legacy} -> {}", target.display()));
            }
        }
    }
    println!("fasthelper install: {} symlink(s) in {}", made.len(), bin.display());
    for m in &made {
        println!("  {m}");
    }
    Ok(())
}
```

**Step 4: CLI plumbing**

In `cli.rs`, add to `Cmd`:

```rust
    /// Symlink the suite binaries (plus legacy names) into <prefix>/bin.
    Install(InstallArgs),
```

and the args struct:

```rust
#[derive(Args, Debug)]
pub struct InstallArgs {
    /// Installation prefix; binaries land in <prefix>/bin.
    #[arg(long, default_value = "/usr/local")]
    pub prefix: PathBuf,
    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,
    /// Skip the legacy binary names (rpolyx, fdedup, novaseq-sim).
    #[arg(long)]
    pub no_legacy: bool,
}
```

`main.rs` dispatch arm: `cli::Cmd::Install(a) => fasthelper::install::cmd_install(&a),`
`lib.rs`: `pub mod install;`

**Step 5: Run the test to verify pass**

Run: `cd tools/fasthelper && cargo test --offline`
Expected: `19 passed; 0 failed`.

**Step 6: Manual verification**

```bash
./tools/fasthelper/target/debug/fasthelper install --prefix /tmp/fhinstall --tools-dir tools
ls -la /tmp/fhinstall/bin/
```
Expected: six canonical symlinks plus `fdedup`, `novaseq-sim`, `rpolyx`; exit 0. Re-running the same command succeeds (links replaced).

**Step 7: Commit**

```bash
git add tools/fasthelper/src/install.rs tools/fasthelper/src/cli.rs \
        tools/fasthelper/src/main.rs tools/fasthelper/src/lib.rs \
        tools/fasthelper/tests/cli.rs
git commit -m "feat(fasthelper): install --prefix with legacy names (D-8)"
```

---

## Task 8: docs + README surface

**Objective:** The new verbs appear in the README toolbox table and the module doc comment in `cli.rs`, so `--help` and the docs agree.

**Files:**
- Modify: `README.md` (fasthelper row, lines ~76-81)
- Modify: `tools/fasthelper/src/cli.rs:1-8` (doc comment: four verbs → six)

**Step 1: Update `cli.rs` doc comment**

```rust
//! CLI surface of fasthelper.
//!
//! Everything is one of six verbs:
//!   run          — execute a FastSuite chain from a TOML config
//!   concatenate  — aggregate per-sample reports into one multi-sample report
//!   compare      — before/after side by side
//!   check        — validate a config before the run starts
//!   doctor       — check that the suite binaries are present and working
//!   install      — symlink the binaries (plus legacy names) into a prefix
//!   schema       — print the report schema (for downstream tooling)
```

**Step 2: Update `README.md`** — extend the fasthelper usage block (after line 81) with:

```bash
fasthelper check -c fastsuite.toml                 # fail at minute zero, not hour three
fasthelper run -c fastsuite.toml -P 8 --from fastprep
fasthelper concatenate runs/ --fail-on polyG_frac>0.05
fasthelper install --prefix ~/.local
```

**Step 3: Run the full suite + the demo**

```bash
cd tools/fasthelper && cargo test --offline
bash examples/demo.sh 20000
./tools/fasthelper/target/release/fasthelper doctor --tools-dir tools
```
Expected: `19 passed; 0 failed`; demo builds and runs the sim control, printing the report.tsv; doctor's table lists `fastcheck … fastcheck 0.1.0` and every other tool.

**Step 4: Commit**

```bash
git add README.md tools/fasthelper/src/cli.rs
git commit -m "docs: fasthelper check/install/-P/--from/--fail-on in the README"
```

---

## Final acceptance (run after every task group)

```bash
cd tools/fasthelper && cargo test --offline        # 19 passed, 0 failed
bash examples/demo.sh 20000                        # end to end: sim -> dedup -> prep -> maps
./tools/fasthelper/target/release/fasthelper doctor --tools-dir tools
```
Doctor must show `fastcheck tools/fastcheck/target/release/fastcheck fastcheck 0.1.0` as the first row; demo must print a `report.tsv` with the sim_E sample; the compare/check/install/threshold tests must all be green.
