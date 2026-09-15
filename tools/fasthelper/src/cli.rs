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
    /// Aggregate per-tool reports (fastsuite.report.v1 or legacy JSON) into one report.
    Report(ReportArgs),
    /// Check that the FastSuite binaries are found and runnable.
    Doctor(DoctorArgs),
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

    /// Print the aggregated table to stdout.
    #[arg(long)]
    pub stdout: bool,
}

#[derive(Args, Debug)]
pub struct DoctorArgs {
    /// Directory holding the tool checkouts (default: this crate's parent, then $FASTSUITE_TOOLS_DIR).
    #[arg(long, env = "FASTSUITE_TOOLS_DIR")]
    pub tools_dir: Option<PathBuf>,
}
