# fastsim: CLI convergence (clap, common flags) + deployment script path audit

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Migrate fastsim's hand-rolled argument parser to clap (same flag names, suite-standard help/errors), add the missing common flags (`--force`, `-q/--quiet`) and the `FASTSUITE_THREADS` env default, and fix the 33 deployment scripts that still point at the pre-rename `novaseq-sim` checkout.

**Architecture:** clap derive replaces `src/cli.rs`'s hand parser; the internal `Config` struct stays (main.rs depends on its field names) and gains `force`/`quiet`. `--force` refuses to clobber an existing output prefix; `-q` gates the end-of-run summary line only (warnings always print, per CLI-SPEC rule 5). `FASTSUITE_THREADS` is read when `-j` is absent — fastsim stays single-threaded and the value is only reported-and-ignored (CLI-SPEC gap #6, decided behaviour: never pretend to parallelise). Deployment scripts get a `FASTSUITE_ROOT` env indirection with a `/scratch/users/$USER/FastSuite` fallback (per the fastsuite skill's container-build guidance); no logic changes.

**Tech Stack:** Rust (clap 4.5.60 derive — the version already pinned by fastdedup/fastcheck; offline build verified), bash (deployment scripts only).

---

## Context (verified 2026-10-12, before any edits)

* `cargo test --offline` in `tools/fastsim`: **54 tests pass** (45 unit + 9 integration).
* The hand parser lives in `tools/fastsim/src/cli.rs` (`pub fn parse(args: &[String]) -> Result<Config>`, lines 83–172); `main.rs:26-28` collects `env::args().skip(1)` and calls it.
* Every existing flag: `--pairs --profile -o|--out|--out1|--out-prefix -j|--threads --json --seed --read-len --insert-mean --insert-sd --insert-min --gc --reference --adapter-seq --no-adapter --no-gzip --single --dup-rate --dup-pool --drop-rate --drop-min-len -V|--version -h|--help`. Post-parse validation (lines 142–170) enforces: `pairs >= 1`, `insert_sd > 0`, `gc/dup_rate/drop_rate ∈ [0,1]`, `drop_min_len >= 1`, `insert_min <= read_len`.
* Old error style (verified): `Error: unknown argument "--bogus" (see --help)`, exit 1. Suite style (fastdedup, verified): `error: unexpected argument '--bogus' found` + usage, exit 2.
* **A prototype of this entire migration was built and verified in a scratch copy** (`~/.hermes/cache/scratch/fastsim_clap_proto`): clap in, all 19 legacy invocation forms parse (parity script below), all 54 tests still pass (after two assertion-wording updates described in Task 7), plus 4 new flag tests. The snippets in this plan are the verified prototype code.
* **Help before the migration** (hand-rolled, 37 lines, captured to `docs/plans/2026-10-12_fastsim-cli-convergence.md` appendix below) vs **help after** (clap-generated, 29 lines). The em-dash and the `REQUIRED:` block disappear; `-f/--force` and `-q/--quiet` appear; defaults are shown by clap.
* `--dry-run` is **deliberately not added**: fastsim is a generator (it always writes), CLI-SPEC §2 already notes `-1/-2` are meaningless for it, and the report JSON's `"dry_run": false` field is the honest record. Say so in the plan; do not add the flag.
* stdin/stdout (`-`) support is **not** added: fastsim writes FASTQ and reads no FASTQ input (D-4 is moot for this tool).
* **Deployment scripts: 33 of 39 files under `tools/fastsim/scripts/` + the 3 `.sbatch` files at `tools/fastsim/` carry `novaseq-sim` paths** (49 occurrences total: 25× `/home/raphael/hermes_scratch/novaseq-sim`, 8× `/scratch/users/ribesr/novaseq-sim`, plus `target/release/novaseq-sim` binary names and profile paths). The Io checkout lives at `/scratch/users/ribesr/FastSuite` (fastsuite skill, verified); `$HOME/hermes_scratch/FastSuite` is absent there. Local checkout: `/home/raphael/hermes_scratch/FastSuite`.
* One dangerous detail: the scripts `cd` into the old checkout and then run `git add -A && git commit` — after the path fix they must commit inside FastSuite, which is correct, but only if the fix lands on the FastSuite checkout (the task explicitly permits it).
* `FASTSUITE_THREADS` is currently honoured by **no** binary (docs/PLAN.md phase 1.4 lists it as the target); fastsim is first.

---

## Phase 1 — `FASTSUITE_THREADS` env default (works with the hand parser; ~15 min)

Do this first: it is pure gain, needs no clap, and the hand parser reads env vars identically.

### Task 1: `FASTSUITE_THREADS` defaults `--threads` when `-j` is absent

**Objective:** `FASTSUITE_THREADS=4 fastsim …` (no `-j`) reports "`--threads 4` … ignored"; explicit `-j` always wins.

**Files:**
- Modify: `tools/fastsim/src/cli.rs` — add env resolution at the top of `parse()`
- Create: `tools/fastsim/tests/env_threads_test.rs`

**Step 1: Write failing test**

```rust
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
```

**Step 2: Run test to verify failure**

Run: `cd tools/fastsim && cargo test --offline --test env_threads_test`
Expected: FAIL — stderr has no `--threads 4` (env var not read yet).

**Step 3: Minimal implementation**

In `tools/fastsim/src/cli.rs`, inside `parse()`, right after `let mut cfg = Config::default();` add:

```rust
    // FASTSUITE_THREADS is the suite-wide default when -j/--threads is absent
    // (CLI-SPEC §1, gap #6). fastsim is single-threaded so the value is only
    // reported-and-ignored, but the contract must hold for the suite.
    let threads_given = args.iter().any(|a| a == "-j" || a == "--threads");
    if !threads_given {
        if let Some(n) = std::env::var("FASTSUITE_THREADS")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
        {
            cfg.threads = n;
        }
    }
```

**Step 4: Run test to verify pass**

Run: `cd tools/fastsim && cargo test --offline --test env_threads_test`
Expected: 1 passed.

Run: `cd tools/fastsim && cargo test --offline` — 55 passed (54 + 1), 0 failed.

**Step 5: Commit**

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastsim/src/cli.rs tools/fastsim/tests/env_threads_test.rs
git commit -m "fastsim: honour FASTSUITE_THREADS as the default thread count (CLI-SPEC gap #6)"
```

---

## Phase 2 — `--force` and `-q/--quiet` in the hand parser (~20 min)

Do these BEFORE the clap migration so the parity gate is measured against a binary that already has the new flags.

### Task 2: `--force` refuses to clobber the output prefix

**Objective:** a second run on the same `-o` prefix fails with a message naming `--force`; `--force`/`-f` overwrites.

**Files:**
- Modify: `tools/fastsim/src/cli.rs` — `Config.force: bool` field + parse arms + validation
- Modify: `tools/fastsim/src/main.rs` — clobber check at the top of `run()`
- Create: `tools/fastsim/tests/common_flags_test.rs`

**Step 1: Write failing test** — create `tools/fastsim/tests/common_flags_test.rs`:

```rust
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
```

**Step 2: Run test to verify failure**

Run: `cd tools/fastsim && cargo test --offline --test common_flags_test`
Expected: FAIL — `force_refuses_to_clobber_without_the_flag`: second run succeeds (no flag exists); `quiet_suppresses_the_summary_line`: summary line present.

**Step 3: Minimal implementation**

In `tools/fastsim/src/cli.rs`:

1. Add the field to `Config` (after `pub single: bool,`):

```rust
    /// Overwrite an existing output prefix instead of refusing.
    pub force: bool,
    /// Suppress the end-of-run summary line (warnings/errors still print).
    pub quiet: bool,
```

2. In `Default for Config`, after `single: false,`:

```rust
            force: false,
            quiet: false,
```

3. In `parse()`'s match, after the `"--single"` arm:

```rust
            "--force" | "-f" => cfg.force = true,
            "-q" | "--quiet" => cfg.quiet = true,
```

In `tools/fastsim/src/main.rs`, at the top of `run()`, before `let t0 = std::time::Instant::now();`:

```rust
    // --force: a generator always writes, so the clobber check is the whole
    // safety net. Probe the three artefacts the run will create.
    if !cfg.force {
        for suffix in [".gt.tsv", "_R1.fastq.gz", "_R2.fastq.gz", "_R1.fastq", "_R2.fastq"] {
            let probe = format!("{}{}", cfg.out, suffix);
            if std::path::Path::new(&probe).exists() {
                anyhow::bail!("{probe} exists; pass --force to overwrite");
            }
        }
    }
```

And gate the summary (the existing `eprintln!("fastsim: {} pairs, ...")` block near the end of `run()`):

```rust
    if !cfg.quiet {
        eprintln!(
            "fastsim: {} pairs, read_len {}, ...", // unchanged body
        );
    }
```

**Step 4: Run test to verify pass**

Run: `cd tools/fastsim && cargo test --offline --test common_flags_test`
Expected: 2 passed.

Run: `cd tools/fastsim && cargo test --offline` — 57 passed, 0 failed.

**Step 5: Commit**

```bash
git add tools/fastsim/src/cli.rs tools/fastsim/src/main.rs tools/fastsim/tests/common_flags_test.rs
git commit -m "fastsim: add --force (refuse to clobber) and -q/--quiet (suppress summary line)"
```

---

## Phase 3 — clap migration (~1 h, the big one)

**Pre-flight:** capture the before-help so the diff is reviewable:

```bash
cd ~/hermes_scratch/FastSuite/tools/fastsim
./target/debug/fastsim --help 2>/dev/null > /tmp/help_before.txt   # 37 lines, hand-rolled
```

### Task 3: Add the clap dependency

**Objective:** clap 4.5.60 (derive) declared, the version already pinned by fastdedup and fastcheck so the offline build resolves.

**Files:**
- Modify: `tools/fastsim/Cargo.toml`

**Step 1:** Add under `[dependencies]`:

```toml
clap = { version = "4.5.60", features = ["derive"] }
```

**Step 2: Verify the build resolves offline**

Run: `cd tools/fastsim && cargo build --offline`
Expected: `Finished dev profile` — clap 4.5.60 is in `~/.cargo/registry/cache` (fastdedup's lock pins it; verified).

**Step 3: Commit**

```bash
git add tools/fastsim/Cargo.toml tools/fastsim/Cargo.lock
git commit -m "fastsim: depend on clap (derive) for the suite-standard CLI"
```

### Task 4: Replace `src/cli.rs` with the clap derive parser

**Objective:** every flag name kept; `-h/--help/-V/--version` clap-generated; errors in suite style; validation moved into `Config::validate()`.

**Files:**
- Modify: `tools/fastsim/src/cli.rs` (full rewrite — the prototype, verified)

**Step 1: Write the new parser** — replace the whole file with:

```rust
//! CLI argument parsing (clap, FastSuite house style).
//!
//! Usage:
//! ```text
//! fastsim --pairs N --profile profiles/mandril_E.json -o prefix \
//!   [--seed S] [--read-len 150] [--insert-mean 300 --insert-sd 50] \
//!   [--gc 0.49] [--reference ref.fa] [--adapter-seq SEQ] [--no-adapter] \
//!   [--no-gzip] [--insert-min 150] [--dup-rate 0.02] \
//!   [--drop-rate 0.05 --drop-min-len 100] [--json report.json]
//! ```
//!
//! `-o/--out1` is the canonical output flag (a prefix for a simulator, not a
//! file); `-j/--threads` is accepted and ignored (the simulator is
//! single-threaded) with a warning; `--json` writes the suite report so
//! fasthelper can aggregate the synthetic control like any other sample.
//!
//! `--dry-run` is deliberately absent: a generator always writes, so there is
//! no report-only mode (CLI-SPEC §2 notes -1/-2 are meaningless here for the
//! same reason).

use anyhow::{bail, Result};
use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[command(
    version,
    about = "NovaSeq-like paired-end read simulator with ground truth",
    long_about = None,
    arg_required_else_help = true
)]
pub struct Cli {
    /// Number of read pairs to simulate.
    #[arg(long)]
    pub pairs: u64,

    /// Profile JSON (artefact CDFs, quality curves, adapter).
    #[arg(long)]
    pub profile: String,

    /// Output prefix: writes PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz, PREFIX.gt.tsv.
    #[arg(
        short = 'o',
        long,
        visible_alias = "out1",
        visible_alias = "out-prefix",
        value_name = "PREFIX"
    )]
    pub out: String,

    /// RNG seed (default 42).
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Read length (default: taken from the profile; an explicit value must match it).
    #[arg(long)]
    pub read_len: Option<u32>,

    /// Insert length mean (default 300). Ignored when the profile carries a measured insert CDF.
    #[arg(long, default_value_t = 300.0)]
    pub insert_mean: f64,

    /// Insert length sd (default 50). Ignored when the profile carries a measured insert CDF.
    #[arg(long, default_value_t = 50.0)]
    pub insert_sd: f64,

    /// Insert length floor (default 150). Ignored when the profile carries a measured insert CDF.
    #[arg(long, default_value_t = 150)]
    pub insert_min: u32,

    /// GC fraction for random DNA (default 0.49).
    #[arg(long, default_value_t = 0.49)]
    pub gc: f64,

    /// Sample inserts from a reference FASTA instead of random DNA.
    #[arg(long)]
    pub reference: Option<String>,

    /// Adapter sequence for readthrough (default: from the profile).
    #[arg(long)]
    pub adapter_seq: Option<String>,

    /// Disable adapter readthrough.
    #[arg(long)]
    pub no_adapter: bool,

    /// Write uncompressed FASTQ (default: gzip).
    #[arg(long)]
    pub no_gzip: bool,

    /// PCR duplicate fraction (default 0.02; mandril measures 1-4%).
    /// A duplicate re-sequences a previously emitted fragment with fresh
    /// quality, errors and artefacts (per-cluster dropout).
    #[arg(long, default_value_t = 0.02)]
    pub dup_rate: f64,

    /// Ring-buffer capacity for the duplicate pool (default 131072).
    #[arg(long, default_value_t = 131072)]
    pub dup_pool: u32,

    /// Fraction of pairs whose fragment is too short to survive trimming
    /// (insert < --drop-min-len; default 0.05).
    #[arg(long, default_value_t = 0.05)]
    pub drop_rate: f64,

    /// Insert ceiling for the forced-drop class (default 100, the bench min_length).
    #[arg(long, default_value_t = 100)]
    pub drop_min_len: u32,

    /// Accepted for FastSuite CLI compatibility; the simulator is
    /// single-threaded, so any value is reported and ignored.
    /// Falls back to the FASTSUITE_THREADS environment variable.
    #[arg(short = 'j', long, default_value_t = 1)]
    pub threads: usize,

    /// Optional FastSuite report (schema v1) destination.
    #[arg(long, value_name = "PATH")]
    pub json: Option<String>,

    /// Single-end mode: one read per fragment (no mate 2, no insert-size
    /// geometry; adapter readthrough still lands at the read's 3' end).
    #[arg(long)]
    pub single: bool,

    /// Force overwrite of existing output files (the default is to refuse).
    #[arg(long, short = 'f')]
    pub force: bool,

    /// Suppress the end-of-run summary line (warnings and errors still print).
    #[arg(long, short = 'q')]
    pub quiet: bool,
}

/// The resolved simulation configuration (field names unchanged so main.rs
/// keeps working verbatim).
#[derive(Debug, Clone)]
pub struct Config {
    pub pairs: u64,
    pub profile: String,
    pub out: String,
    pub seed: u64,
    pub read_len: Option<u32>,
    pub insert_mean: f64,
    pub insert_sd: f64,
    pub insert_min: u32,
    pub gc: f64,
    pub reference: Option<String>,
    pub adapter_seq: Option<String>,
    pub gzip: bool,
    pub dup_rate: f64,
    pub dup_pool: u32,
    pub drop_rate: f64,
    pub drop_min_len: u32,
    pub threads: usize,
    pub json: Option<String>,
    pub single: bool,
    pub force: bool,
    pub quiet: bool,
}

impl From<Cli> for Config {
    fn from(cli: Cli) -> Self {
        Config {
            pairs: cli.pairs,
            profile: cli.profile,
            out: cli.out,
            seed: cli.seed,
            read_len: cli.read_len,
            insert_mean: cli.insert_mean,
            insert_sd: cli.insert_sd,
            insert_min: cli.insert_min,
            gc: cli.gc,
            reference: cli.reference,
            adapter_seq: if cli.no_adapter {
                Some(String::new()) // empty = disabled, as the hand parser did
            } else {
                cli.adapter_seq
            },
            gzip: !cli.no_gzip,
            dup_rate: cli.dup_rate,
            dup_pool: cli.dup_pool,
            drop_rate: cli.drop_rate,
            drop_min_len: cli.drop_min_len,
            threads: cli.threads,
            json: cli.json,
            single: cli.single,
            force: cli.force,
            quiet: cli.quiet,
        }
    }
}

impl Config {
    /// The post-parse validation the hand parser performed (identical
    /// messages, so downstream scripts/tests keep matching).
    pub fn validate(&self) -> Result<()> {
        if self.pairs == 0 {
            bail!("--pairs must be >= 1");
        }
        if self.insert_sd <= 0.0 {
            bail!("--insert-sd must be > 0");
        }
        if !(0.0..=1.0).contains(&self.gc) {
            bail!("--gc must be in [0,1]");
        }
        if !(0.0..=1.0).contains(&self.dup_rate) {
            bail!("--dup-rate must be in [0,1]");
        }
        if !(0.0..=1.0).contains(&self.drop_rate) {
            bail!("--drop-rate must be in [0,1]");
        }
        if self.drop_min_len == 0 {
            bail!("--drop-min-len must be >= 1");
        }
        if let Some(l) = self.read_len {
            if self.insert_min > l {
                bail!(
                    "--insert-min ({}) should be <= read-len ({})",
                    self.insert_min,
                    l
                );
            }
        }
        Ok(())
    }
}

/// Parse the process argv into a validated `Config`.
/// `-V/--version` and `-h/--help` print and exit 0, clap-style.
pub fn parse() -> Result<Config> {
    let mut cli = Cli::parse();
    // FASTSUITE_THREADS is the suite-wide default when -j is absent
    // (CLI-SPEC §1). clap's `env` attribute can't express "only when the
    // flag is absent AND the value is valid", so resolve it manually.
    let threads_given = std::env::args().any(|a| a == "-j" || a == "--threads");
    if !threads_given {
        if let Some(n) = std::env::var("FASTSUITE_THREADS")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
        {
            cli.threads = n;
        }
    }
    let cfg: Config = cli.into();
    cfg.validate()?;
    Ok(cfg)
}
```

Note: `--no-adapter` is folded into `Config.adapter_seq = Some("")` exactly as the hand parser did (`main.rs` treats an empty adapter as `AdapterMode::None`) — the simulation path is byte-identical.

**Step 2: Update `main.rs` to the new entry point**

In `tools/fastsim/src/main.rs`, replace lines 25–29:

```rust
fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cfg = cli::parse(&args)?;
    run(&cfg)
}
```

with:

```rust
fn main() -> Result<()> {
    let cfg = cli::parse()?;
    run(&cfg)
}
```

Everything else in `main.rs` is untouched (the `Config` field set is a superset of the old one). The `--threads` warning and the `--force` clobber check added in Tasks 1–2 live in `run()` and stay as they are.

**Step 3: Run the parity gate**

Run the full legacy-invocation parity script (see Appendix A) against the new binary:

```bash
cd tools/fastsim
/tmp/parity_check.sh ./target/debug/fastsim
```

Expected: **19/19 OK** — every pre-migration invocation form still parses (verified in the prototype).

**Step 4: Run the test suite and fix the two assertion wordings**

Run: `cd tools/fastsim && cargo test --offline`
Expected: 2 FAILED — `unknown_flag_is_an_error` (clap says "unexpected argument", not "unknown argument") and `missing_profile_is_an_error` (clap prints "the following required arguments were not provided:\n  --profile <PROFILE>", not "--profile is required"). Both are contract-level message changes the migration is allowed to make (CLI-SPEC: the long name is canonical; the error STYLE is what converges). Update the two assertions in `tools/fastsim/tests/cli_test.rs`:

```rust
#[test]
fn missing_profile_is_an_error() {
    let status = Command::new(bin())
        .args(["--pairs", "10", "--out", "x"])
        .output()
        .unwrap();
    assert!(!status.status.success());
    let err = String::from_utf8_lossy(&status.stderr);
    // clap reports a missing required argument with its own wording; the
    // contract is "non-zero exit + the flag named on stderr".
    assert!(err.contains("--profile"), "stderr: {err}");
}
```

and in `unknown_flag_is_an_error`:

```rust
    assert!(err.contains("unexpected argument"), "stderr: {err}");
```

**Step 5: Run tests to verify pass**

Run: `cd tools/fastsim && cargo test --offline`
Expected: **58 passed** (45 unit + 13 integration: 9 original + 1 env + 2 common-flags + the two reworded assertions unchanged in count), 0 failed.

**Step 6: Verify the error-style convergence by hand**

```bash
cd tools/fastsim
./target/debug/fastsim --bogus
```

Expected (clap style, matches fastdedup, exit 2):

```
error: unexpected argument '--bogus' found

Usage: fastsim [OPTIONS] --pairs <PAIRS> --profile <PROFILE> --out <PREFIX>

For more information, try '--help'.
```

```bash
./target/debug/fastsim -V      # → "fastsim 0.1.0", exit 0
./target/debug/fastsim -h      # → clap help, exit 0
./target/debug/fastsim         # → clap help (arg_required_else_help), exit 2
```

**Step 7: Commit**

```bash
git add tools/fastsim/src/cli.rs tools/fastsim/src/main.rs tools/fastsim/tests/cli_test.rs
git commit -m "fastsim: clap migration — every flag name kept, help/version/errors now match the suite (CLI-SPEC gap #8)"
```

### Task 5: Merge the flag tests into one file

**Objective:** one `tests/cli_flags_test.rs` covering `--force`, `-q`, `FASTSUITE_THREADS` and help/version (the prototype's verified file, 4 tests).

**Files:**
- Create: `tools/fastsim/tests/cli_flags_test.rs`
- Delete: `tools/fastsim/tests/common_flags_test.rs`, `tools/fastsim/tests/env_threads_test.rs`

**Step 1:** Write `tools/fastsim/tests/cli_flags_test.rs` — the union of the three test bodies above plus:

```rust
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
```

(with the three shared helpers `bin()`, `write_test_profile()`, and the three test fns from Phases 1–2, each `#[test]`).

**Step 2:** Delete the two now-redundant files:

```bash
git rm tools/fastsim/tests/common_flags_test.rs tools/fastsim/tests/env_threads_test.rs
```

**Step 3: Verify**

Run: `cd tools/fastsim && cargo test --offline`
Expected: **58 passed**, 0 failed (same set, one file).

**Step 4: Commit**

```bash
git add tools/fastsim/tests/cli_flags_test.rs
git commit -m "fastsim: consolidate the common-flag tests into cli_flags_test.rs"
```

---

## Phase 4 — Deployment script path audit + FASTSUITE_ROOT indirection (~45 min)

**Scope discipline (from the fastsuite skill):** update ONLY the path strings. No logic changes, no reordering, no new steps. The `FASTSUITE_ROOT` indirection goes at the top of each script as a two-liner; the `cd` targets and the remote `scp`/`ssh` paths then use it.

### Task 6: Audit — list every file and every path form

**Objective:** a written audit (this section doubles as it) so no script is missed and nothing outside `tools/fastsim/scripts/` + the three `.sbatch` files is touched.

**Files:** audit only — read `tools/fastsim/scripts/*.sh`, `tools/fastsim/scripts/*.py`, `tools/fastsim/*.sbatch`.

**Step 1: Run the audit**

```bash
cd ~/hermes_scratch/FastSuite
grep -rn "novaseq" tools/fastsim/scripts/ tools/fastsim/*.sbatch | grep -v Binary | sort
```

**Expected (verified 2026-10-12): 33 files, 49 occurrences:**

| file | occurrences | what it carries |
|---|---|---|
| `scripts/rebuild_deploy.sh` | 5 | local cd + 3 remote scp targets + sbatch path |
| `scripts/local_real_profile.sh` | 4 | local cd, scp source, run line, python heredoc path |
| `scripts/deploy_resubmit.sh` | 4 | local cd + scp targets + sbatch path |
| `scripts/run_diag_q.sh` | 3 | local cd + scp target + ssh python path |
| `validate.sbatch` | 3 | `#SBATCH --job-name=novaseq-validate`, `cd`, `BIN=bin/novaseq-sim` |
| `validate2.sbatch` | 2 | `--job-name`, `cd` |
| `scripts/local_run_sim.sh`, `local_closed_loop.sh`, `laptop_closed_loop.sh`, `check_profile.sh` | 2 each | cd + binary name / remote path |
| `scripts/build_profile_srun.sh`, `resubmit_build.sh` | 1 each | remote `cd /scratch/users/ribesr/novaseq-sim` inside ssh |
| `scripts/{test_build_local,run_validation,run_validation2,run_plot_heads_tails,run_dump_profile_vals,run_dump_only,run_diag_residual,run_check_offset_bias,resim_validate,resim_only,resim_nogzip,regzip,rebuild_all,pixi_install_analysis,inspect_local,commit_validation,commit_final,check_sim_integrity,check_gzip}.sh` | 1 each | `cd /home/raphael/hermes_scratch/novaseq-sim` (or `.../local_work/sim`) |
| `scripts/diag_joint.py`, `diag_q.py`, `local_check.py` | 1 each | hardcoded profile path inside Python |
| `build_profile.sbatch` | 1 | `cd /scratch/users/ribesr/novaseq-sim` |

The three binary-name occurrences (`./target/release/novaseq-sim` in `local_closed_loop.sh`, `local_real_profile.sh`, `local_run_sim.sh`, `BIN=bin/novaseq-sim` in `validate.sbatch`) matter most: the binary is now `fastsim`, so those runs would fail outright even if the `cd` were fixed.

**Step 2: Commit the audit note**

Append the table above to `docs/plans/2026-10-12_fastsim-cli-convergence.md` (this file) — no code commit needed for the audit itself; proceed to Task 7.

**Audit executed 2026-10-12 (final count):** the same grep returns **55 occurrences across 35 files** (the table above is the plan's estimate; the repo carries 25 `.sh` with the plain local `cd` rather than 24, because `check_profile.sh` and `laptop_closed_loop.sh` each carry that cd alongside a second novaseq path). The per-file shape matches the table exactly, so the substitution below covers every occurrence. `src/fastq.rs`'s `header_format_matches_novaseq_style` test name, `README.md`'s "(formerly novaseq-sim)" and `CITATION.cff`'s keywords are historical and stay untouched.

### Task 7: The `FASTSUITE_ROOT` two-liner

**Objective:** every script resolves its checkout root once, local-first with the cluster fallback, and every `novaseq-sim` path string becomes `$FASTSUITE_ROOT` (or the renamed `fastsim`).

**Files:** all 33 files from the audit.

**Step 1: The canonical header** — insert immediately after the shebang (and comment block) of every local script:

```bash
# FastSuite checkout root: local first, then the Io cluster scratch fallback
# (fastsuite skill: $HOME/hermes_scratch/FastSuite is absent on Io).
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
```

And for the three `.sbatch` files (cluster-side; no `$HOME` there):

```bash
FASTSUITE_ROOT="${FASTSUITE_ROOT:-/scratch/users/$USER/FastSuite}"
```

**Step 2: The substitutions (per file, mechanical — no logic touched):**

| old string | new string | files |
|---|---|---|
| `cd /home/raphael/hermes_scratch/novaseq-sim` | `cd "$FASTSUITE_ROOT"` | 24 `.sh` files |
| `cd /home/raphael/hermes_scratch/novaseq-sim/local_work/sim` | `cd "$FASTSUITE_ROOT/tools/fastsim/local_work/sim"` | `check_gzip.sh`, `regzip.sh` |
| `cd /scratch/users/ribesr/novaseq-sim` (inside ssh '…') | `cd "$FASTSUITE_ROOT"` with the var **exported into the ssh command string** (see step 3) | `check_profile.sh`, `build_profile_srun.sh`, `resubmit_build.sh`, `run_diag_q.sh` |
| `Io:/scratch/users/ribesr/novaseq-sim/...` (scp targets) | `Io:$FASTSUITE_ROOT/...` — the remote side resolves the env var **on the cluster**, which has the same two-liner in its shell profile… **not guaranteed**; use the literal cluster path instead: `Io:/scratch/users/$USER/FastSuite/...` | `rebuild_deploy.sh`, `deploy_resubmit.sh`, `run_diag_q.sh`, `local_real_profile.sh` |
| `./target/release/novaseq-sim` | `./target/release/fastsim` | `local_closed_loop.sh`, `local_real_profile.sh`, `local_run_sim.sh` |
| `BIN=bin/novaseq-sim` | `BIN="$FASTSUITE_ROOT/tools/fastsim/target/release/fastsim"` | `validate.sbatch` |
| `#SBATCH --job-name=novaseq-validate*` | `#SBATCH --job-name=fastsim-validate*` | `validate.sbatch`, `validate2.sbatch` |
| `profiles/mandril_E.json` paths under `/home/raphael/hermes_scratch/novaseq-sim` | `"$FASTSUITE_ROOT/tools/fastsim/profiles/mandril_E.json"` — **but** `mandril_E.json` is not in the FastSuite checkout (`profiles/` holds `synthetic_novaseq.json` and `test_small.json`; the mandril profile lives on the cluster). Keep the cluster literal: `/scratch/users/$USER/FastSuite/tools/fastsim/...` only if it exists there; otherwise leave the profile path alone and note it | `diag_joint.py`, `diag_q.py`, `local_check.py`, `local_real_profile.sh` |

**Python files:** the three `.py` scripts open absolute paths. Convert to env-with-fallback, no logic change:

```python
import os
ROOT = os.environ.get("FASTSUITE_ROOT", "/scratch/users/ribesr/FastSuite")
prof = json.load(open(f"{ROOT}/tools/fastsim/profiles/mandril_E.json"))
```

**Step 3: The ssh/scp wrinkle** — a remote `cd` inside `ssh Io '...'` runs on the cluster, where `FASTSUITE_ROOT` may not be set. Two honest options, per file:
- **Preferred:** export the fallback literal — `ssh Io "cd $FASTSUITE_ROOT && ..."` only works if the remote shell expands it; since single quotes prevent local expansion, switch those commands to double quotes **only where the command has no other `$` or backtick** (verified: `check_profile.sh`, `build_profile_srun.sh`, `resubmit_build.sh`, `run_diag_q.sh` qualify).
- **Where the command embeds nested quotes** (the `check_profile.sh` python heredoc), leave the literal `/scratch/users/ribesr/FastSuite` and add a comment that this path is cluster-side.

**Step 4: Verify nothing else changed**

```bash
cd ~/hermes_scratch/FastSuite
git diff --stat tools/fastsim/scripts tools/fastsim/*.sbatch
# Expected: 33 files, small line counts (1–3 lines each), no logic hunks.
git diff tools/fastsim/scripts/resim_only.sh   # eyeball one end to end
```

Expected diff for `resim_only.sh` (whole file):

```diff
 #!/bin/bash
 # Re-simulate 20M with the offset-centering fix, then validate.
 set -e
-cd /home/raphael/hermes_scratch/novaseq-sim
+FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
+[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
+cd "$FASTSUITE_ROOT"
 ./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out local_work/sim/sim --seed 7 --no-gzip
 echo "SIM_DONE"
 ls -la local_work/sim/
```

**Step 5: Verify no `novaseq` remains**

```bash
grep -rn "novaseq" tools/fastsim/scripts/ tools/fastsim/*.sbatch | grep -v Binary
```

Expected: **no output** (the only remaining `novaseq` mentions in the repo are the historical ones in `docs/RENAME-MAP.md` and `docs/DECISIONS.md`, which are correct as history — do not touch them).

**Step 6: Syntax-check every script**

```bash
for f in tools/fastsim/scripts/*.sh tools/fastsim/*.sbatch; do bash -n "$f" || echo "SYNTAX FAIL: $f"; done
```

Expected: no failures.

**Step 7: Commit**

```bash
git add tools/fastsim/scripts tools/fastsim/*.sbatch
git commit -m "fastsim: deployment scripts — novaseq-sim paths → FASTSUITE_ROOT with /scratch/users/\$USER/FastSuite fallback (path strings only)"
```

---

## Phase 5 — Final verification (~10 min)

### Task 8: Full suite + closed-loop smoke test

**Step 1: Unit + integration tests**

```bash
cd ~/hermes_scratch/FastSuite/tools/fastsim
cargo test --offline
```

Expected: **58 passed; 0 failed** (45 unit + 13 integration).

**Step 2: Help before/after comparison** — capture both:

```bash
cd ~/hermes_scratch/FastSuite/tools/fastsim
./target/debug/fastsim --help 2>/dev/null | tee /tmp/help_after.txt | head -5
diff /tmp/help_before.txt /tmp/help_after.txt | head -20
```

Expected: the hand-rolled 37-line help is replaced by the 29-line clap help (appendix below); every flag still present; `-f, --force` and `-q, --quiet` added; `[aliases: --out1, --out-prefix]` shown.

**Step 3: Closed-loop smoke test** (the exact command from the task):

```bash
cd ~/hermes_scratch/FastSuite
tools/fastcheck/target/release/fastcheck --insert-hist --insert-pairs 2000 --profile /tmp/smoke_p.json /tmp/smoke_R1.fastq /tmp/smoke_R2.fastq
tools/fastsim/target/debug/fastsim --profile /tmp/smoke_p.json --pairs 20000 -o /tmp/sim2
```

Expected (verified in the prototype):

```
fastcheck: wrote fastsim profile to /tmp/smoke_p.json
fastsim: 20000 pairs, read_len 150, full-artefact 0 (0.0000%), adapter-readthrough 1028 (5.1400%), duplicates 362 (1.8100%), must-drop 1028 (5.1400%)
ls /tmp/sim2*   # sim2.gt.tsv sim2_R1.fastq.gz sim2_R2.fastq.gz
```

**Step 4: Regression — the parity gate once more on the release build**

```bash
cd tools/fastsim && cargo build --release --offline
/tmp/parity_check.sh ./target/release/fastsim   # 19/19 OK
```

**Step 5: Commit the verification log**

```bash
git add -A tools/fastsim
git commit -m "fastsim: CLI convergence complete — clap, --force/-q, FASTSUITE_THREADS, 58 tests green"
```

**Step 6: Update the docs the plan closes**

- `docs/CLI-SPEC.md` §4 gap #8: mark "fastsim adopting clap" **done**, and add `--force`/`-q` to the fastsim row (they are now ✓, not gap).
- `docs/CLI-SPEC.md` §4 gap #6 (`FASTSUITE_THREADS`): note "fastsim done; fastprep/fastmaps still to do".
- `docs/PLAN.md` phase 1.5: replace the "optional/cosmetic" note with a `[done 2026-10-12 — see docs/plans/2026-10-12_fastsim-cli-convergence.md]` pointer.

```bash
git add docs/CLI-SPEC.md docs/PLAN.md
git commit -m "docs: close CLI-SPEC gaps #6/#8 for fastsim, point PLAN 1.5 at the plan"
```

---

## Appendix A — the parity gate script

Save as `/tmp/parity_check.sh`; run against any fastsim binary:

```bash
#!/bin/bash
# Every invocation form the pre-clap suite used must still parse identically.
BIN=$1
fail=0
check() {
  desc="$1"; shift
  if "$BIN" "$@" >/dev/null 2>&1; then echo "OK   $desc"; else echo "FAIL $desc (rc=$?)"; fail=1; fi
}
P=/tmp/smoke_p.json
check "canonical -o"            --pairs 10 --profile $P -o /tmp/par_a
check "legacy --out"            --pairs 10 --profile $P --out /tmp/par_b
check "alias --out1"            --pairs 10 --profile $P --out1 /tmp/par_c
check "alias --out-prefix"      --pairs 10 --profile $P --out-prefix /tmp/par_d
check "seed"                    --pairs 10 --profile $P -o /tmp/par_e --seed 7
check "read-len match"          --pairs 10 --profile $P -o /tmp/par_f --read-len 150
check "insert overrides"        --pairs 10 --profile $P -o /tmp/par_g --insert-mean 300 --insert-sd 50 --insert-min 150
check "gc"                      --pairs 10 --profile $P -o /tmp/par_h --gc 0.5
check "no-adapter"              --pairs 10 --profile $P -o /tmp/par_i --no-adapter
check "adapter-seq"             --pairs 10 --profile $P -o /tmp/par_j --adapter-seq ACGTACGT
check "no-gzip"                 --pairs 10 --profile $P -o /tmp/par_k --no-gzip
check "dup controls"            --pairs 100 --profile $P -o /tmp/par_l --dup-rate 0.1 --dup-pool 1024
check "drop controls"           --pairs 100 --profile $P -o /tmp/par_m --drop-rate 0.1 --drop-min-len 80
check "single mode"             --pairs 100 --profile $P -o /tmp/par_n --single --no-gzip
check "json report"             --pairs 100 --profile $P -o /tmp/par_o --json /tmp/par_o.json
check "threads accepted"        --pairs 100 --profile $P -o /tmp/par_p -j 4
check "threads long"            --pairs 100 --profile $P -o /tmp/par_q --threads 2
check "help/version exit 0"     --help
check "combined"                --pairs 500 --profile $P -o /tmp/par_r --seed 42 --no-gzip --dup-rate 0.02 --drop-rate 0.05 --drop-min-len 100 --json /tmp/par_r.json --single
exit $fail
```

(The profile `/tmp/smoke_p.json` comes from the Task 8 fastcheck run; before that, use `tools/fastsim/profiles/test_small.json` — both work.)

## Appendix B — help output, before and after

**Before** (hand-rolled, 37 lines, captured 2026-10-12 from `target/debug/fastsim --help`):

```text
fastsim — NovaSeq-like paired-end read simulator with ground truth

USAGE:
    fastsim --pairs N --profile FILE --out PREFIX [OPTIONS]

REQUIRED:
    --pairs N              number of read pairs to simulate
    --profile FILE         profile JSON (artefact CDFs, quality curves, adapter)
    --out PREFIX           output prefix (writes PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz,
                           PREFIX.gt.tsv)

OPTIONS:
    --seed S               RNG seed (default 42)
    --read-len L           read length (default 150)
    --insert-mean M        insert length mean (default 300)
    --insert-sd S          insert length sd (default 50)
    --insert-min M         insert length floor (default 150)
    --gc G                 GC fraction for random DNA (default 0.49)
    --reference FILE       sample inserts from a reference FASTA instead of random DNA
    --adapter-seq SEQ      adapter sequence for readthrough (default: Nextera
                           CTGTCTCTTATACACATCT from profile)
    --no-adapter           disable adapter readthrough
    --no-gzip              write uncompressed FASTQ
    --dup-rate R           PCR duplicate fraction (default 0.02; mandril
                           measures 1-4%). A duplicate re-sequences a
                           previously emitted fragment with fresh quality,
                           errors and artefacts (per-cluster dropout).
    --dup-pool N           ring-buffer capacity for the duplicate pool
                           (default 131072)
    --drop-rate R          fraction of pairs whose fragment is too short to
                           survive trimming (insert < --drop-min-len, so the
                           GT keep window is below min_length on both mates
                           and a correct trimmer must drop the pair)
                           (default 0.05)
    --drop-min-len L       insert ceiling for the forced-drop class
                           (default 100, matches the bench min_length)
    -h, --help             this help
```

**After** (clap-generated, 29 lines, captured from the verified prototype):

```text
NovaSeq-like paired-end read simulator with ground truth

Usage: fastsim [OPTIONS] --pairs <PAIRS> --profile <PROFILE> --out <PREFIX>

Options:
      --pairs <PAIRS>                Number of read pairs to simulate
      --profile <PROFILE>            Profile JSON (artefact CDFs, quality curves, adapter)
  -o, --out <PREFIX>                 Output prefix: writes PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz, PREFIX.gt.tsv [aliases: --out1, --out-prefix]
      --seed <SEED>                  RNG seed (default 42) [default: 42]
      --read-len <READ_LEN>          Read length (default: taken from the profile; an explicit value must match it)
      --insert-mean <INSERT_MEAN>    Insert length mean (default 300). Ignored when the profile carries a measured insert CDF [default: 300]
      --insert-sd <INSERT_SD>        Insert length sd (default 50). Ignored when the profile carries a measured insert CDF [default: 50]
      --insert-min <INSERT_MIN>      Insert length floor (default 150). Ignored when the profile carries a measured insert CDF [default: 150]
      --gc <GC>                      GC fraction for random DNA (default 0.49) [default: 0.49]
      --reference <REFERENCE>        Sample inserts from a reference FASTA instead of random DNA
      --adapter-seq <ADAPTER_SEQ>    Adapter sequence for readthrough (default: from the profile)
      --no-adapter                   Disable adapter readthrough
      --no-gzip                      Write uncompressed FASTQ (default: gzip)
      --dup-rate <DUP_RATE>          PCR duplicate fraction (default 0.02; mandril measures 1-4%). A duplicate re-sequences a previously emitted fragment with fresh quality, errors and artefacts (per-cluster dropout) [default: 0.02]
      --dup-pool <DUP_POOL>          Ring-buffer capacity for the duplicate pool (default 131072) [default: 131072]
      --drop-rate <DROP_RATE>        Fraction of pairs whose fragment is too short to survive trimming (insert < --drop-min-len; default 0.05) [default: 0.05]
      --drop-min-len <DROP_MIN_LEN>  Insert ceiling for the forced-drop class (default 100, the bench min_length) [default: 100]
  -j, --threads <THREADS>            Accepted for FastSuite CLI compatibility; the simulator is single-threaded, so any value is reported and ignored [default: 1]
      --json <PATH>                  Optional FastSuite report (schema v1) destination
      --single                       Single-end mode: one read per fragment (no mate 2, no insert-size geometry; adapter readthrough still lands at the read's 3' end)
  -f, --force                        Force overwrite of existing output files (the default is to refuse)
  -q, --quiet                        Suppress the end-of-run summary line (warnings and errors still print)
  -h, --help                         Print help
  -V, --version                      Print version
```

**Error style, before vs after:**

```text
# before (hand parser)
Error: unknown argument "--bogus" (see --help)          # exit 1

# after (clap, matches fastdedup/fastprep/fastcheck)
error: unexpected argument '--bogus' found              # exit 2
Usage: fastsim [OPTIONS] --pairs <PAIRS> --profile <PROFILE> --out <PREFIX>
For more information, try '--help'.
```

## Appendix C — what this plan deliberately does NOT do

* **No `--dry-run`.** fastsim is a generator; it always writes. The report JSON records `"dry_run": false` and that is the honest contract (CLI-SPEC §2 already excuses `-1/-2` for the same reason).
* **No stdin/stdout (`-`) support.** fastsim reads no FASTQ; D-4's zero-intermediate pipe is moot for a simulator.
* **No `-v/--verbose`, `--gz-level`, `--temp-dir`, `--strict`.** Not needed by a generator; adding them would be YAGNI.
* **No threading.** `-j` is accepted and ignored **with a warning** — decided behaviour (fastsuite skill: "Do not present fastsim timings as a threading result"). Nobody should "fix" this by faking parallelism.
* **No logic changes in the deployment scripts.** Path strings only (fastsuite skill: "Update them only when the checkout itself is renamed" — the rename already happened; this is that update).
* **The mandril profile path inside the Python diag scripts** stays cluster-literal (`/scratch/users/…`) where the local checkout does not contain `profiles/mandril_E.json`; noted per-file rather than silently rewritten.
