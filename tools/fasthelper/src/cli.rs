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

use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "fasthelper",
    version,
    about = "FastSuite helper: run the chain, aggregate the reports, check the install",
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Cmd,
}

#[derive(Subcommand, Debug)]
pub enum Cmd {
    /// Run a FastSuite chain (dedup -> prep -> maps) described by a TOML config.
    Run(RunArgs),
    /// Concatenate per-sample reports into one table (the MultiQC slot):
    /// discovery, normalisation, html/tsv/md/json.
    #[command(visible_alias = "report")]
    Concatenate(ReportArgs),
    /// Compare two report sets (e.g. raw vs processed) side by side.
    Compare(CompareArgs),
    /// Check that the FastSuite binaries are found and runnable.
    Doctor(DoctorArgs),
    /// Validate a run config: inputs exist, tools resolve, steps and flags are known.
    Check(CheckArgs),
    /// Symlink the suite binaries (plus legacy names) into <prefix>/bin.
    Install(InstallArgs),
    /// Print the report schema (version, fields, how each tool contributes).
    Schema,
    /// Measure adapter residual in a trimmed FASTQ and fail if it is too high.
    ///
    /// This is the standalone form of the guard that `run --fail-on-adapter-
    /// residual` applies to a step: it measures a trimmed file directly, so a
    /// user can check their own output without running the whole chain.
    Guard(GuardArgs),
}

#[derive(Args, Debug)]
pub struct GuardArgs {
    /// The trimmed reads to check (mate 1).
    #[arg(short = '1', long = "in1")]
    pub in1: PathBuf,

    /// Mate 2, when the data is paired.
    #[arg(short = '2', long = "in2")]
    pub in2: Option<PathBuf>,

    /// The adapter catalogue used for the measurement.
    #[arg(long = "adapter-preset", default_value = "illumina")]
    pub adapter_preset: String,

    /// Residual fraction above which the guard trips (0.01 = 1%).
    #[arg(long, default_value_t = 0.01)]
    pub threshold: f64,

    /// Relative fall required from the raw fraction. Ignored here: there is no
    /// raw input to compare against, so the guard trips on the residual alone.
    #[arg(long, default_value_t = 0.5)]
    pub min_reduction: f64,

    /// Exit 2 when the residual exceeds the threshold (the `--fail-on` convention).
    #[arg(long)]
    pub fail_on_adapter_residual: bool,

    /// Thread count for the underlying fastcheck pass.
    #[arg(short = 'j', long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(Args, Debug)]
pub struct RunArgs {
    /// Path to the run config (TOML).
    #[arg(short = 'c', long)]
    pub config: PathBuf,

    /// Override the config's outdir.
    #[arg(short = 'o', long)]
    pub outdir: Option<PathBuf>,

    /// Override the config's thread count.
    #[arg(short = 'j', long)]
    pub threads: Option<usize>,

    /// Samples to process in parallel (default 1). Threads are divided across the pool.
    #[arg(short = 'P', long, default_value_t = 1)]
    pub jobs: usize,

    /// Re-run starting at this step (earlier steps are resumed/skipped as usual).
    #[arg(long)]
    pub from: Option<String>,

    /// Print the commands without executing them.
    #[arg(long)]
    pub dry_run: bool,

    /// Re-run every step even when a previous run's outputs are still there,
    /// overwriting them. Forwarded as --force to the tools that accept it.
    #[arg(long)]
    pub force: bool,

    /// Keep intermediate files instead of deleting them after the last consumer.
    #[arg(long)]
    pub keep_intermediate: bool,

    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,

    /// Disable the post-trim residual-adapter guard. The guard is on by default;
    /// it only runs when the chain has a trimming step (fastprep).
    #[arg(long)]
    pub no_adapter_guard: bool,

    /// Residual adapter fraction (of trimmed reads) that trips the guard.
    /// Conservative default: a working trim leaves essentially none.
    #[arg(long, value_name = "FLOAT", default_value_t = 0.05)]
    pub adapter_residual_threshold: f64,

    /// Relative fall required from the raw adapter fraction (0.5 = at least
    /// halved). If the trimmed fraction is above this, the guard trips.
    #[arg(long, value_name = "FLOAT", default_value_t = 0.5)]
    pub adapter_min_reduction: f64,

    /// Exit 2 (the --fail-on convention) when the residual-adapter guard trips,
    /// instead of printing a warning.
    #[arg(long)]
    pub fail_on_adapter_residual: bool,
}

#[derive(Args, Debug)]
pub struct ReportArgs {
    /// Report files or directories to scan (directories are walked for *.json).
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output directory for report.{json,tsv,md,html}.
    #[arg(short = 'o', long, default_value = ".")]
    pub outdir: PathBuf,

    /// Comma-separated formats to write (json,tsv,md,html).
    #[arg(long, default_value = "json,tsv,md,html")]
    pub format: String,

    /// Sample name to use for a single file input (default: file stem).
    #[arg(long)]
    pub sample: Option<String>,

    /// Concatenate per-position TSVs (from `fastcheck --positions-tsv`) into
    /// one positions x samples matrix instead of concatenating JSON reports.
    #[arg(long)]
    pub positions: bool,

    /// Column to place in the wide matrix (any column of the input table).
    #[arg(long, default_value = "mean_q", value_name = "NAME")]
    pub metric: String,

    /// Emit the tidy shape (sample, pos, metric, value) with every metric,
    /// instead of a wide single-metric matrix.
    #[arg(long)]
    pub all_metrics: bool,

    /// Print the aggregated table to stdout.
    #[arg(long)]
    pub stdout: bool,

    /// Fail (exit 2) when a metric breaks the rule, e.g. --fail-on polyG_frac>0.05.
    /// Repeatable. Exit 2 = "ran fine, but a threshold tripped" (CLI-SPEC §5).
    #[arg(long = "fail-on", value_name = "METRIC>OP>VALUE")]
    pub fail_on: Vec<String>,
}

#[derive(Args, Debug)]
pub struct CompareArgs {
    /// Before (e.g. raw) report files or directories to scan for *.json.
    #[arg(long, required = true)]
    pub before: Vec<PathBuf>,

    /// After (e.g. processed) report files or directories to scan for *.json.
    #[arg(long, required = true)]
    pub after: Vec<PathBuf>,

    /// Output directory for compare.{tsv,html}.
    #[arg(short = 'o', long, default_value = ".")]
    pub outdir: PathBuf,
}

#[derive(Args, Debug)]
pub struct DoctorArgs {
    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,
}

#[derive(Args, Debug)]
pub struct CheckArgs {
    /// Path to the run config (TOML).
    #[arg(short = 'c', long)]
    pub config: PathBuf,
    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,
}

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
