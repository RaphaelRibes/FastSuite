//! CLI surface of fasthelper.
//!
//! Everything is one of four verbs:
//!   run     — execute a FastSuite chain from a TOML config
//!   report  — aggregate per-tool reports into one multi-sample report
//!   doctor  — check that the suite binaries are present and working
//!   schema  — print the report schema (for downstream tooling)

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
    /// Print the report schema (version, fields, how each tool contributes).
    Schema,
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

    /// Keep intermediate files instead of deleting them after the last consumer.
    #[arg(long)]
    pub keep_intermediate: bool,

    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,
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
