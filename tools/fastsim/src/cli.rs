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
//! `-o/--out` is the canonical output flag (a prefix for a simulator, not a
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
    /// messages, so downstream scripts/tests keep matching). The
    /// --profile/--out presence checks are gone: clap enforces both as
    /// required arguments before this runs.
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
    if cfg.threads > 1 {
        eprintln!(
            "fastsim: warning: simulation is single-threaded; --threads {} accepted for \
             FastSuite CLI compatibility and ignored",
            cfg.threads
        );
    }
    cfg.validate()?;
    Ok(cfg)
}
