//! CLI argument parsing (hand-parsed, no clap — house style).
//!
//! Usage:
//! ```text
//! fastsim --pairs N --profile profiles/mandril_E.json -o prefix \
//!   [--seed S] [--read-len 150] [--insert-mean 300 --insert-sd 50] \
//!   [--gc 0.49] [--reference ref.fa] [--adapter-seq SEQ] [--no-adapter] \
//!   [--no-gzip] [--insert-min 150] [--dup-rate 0.02] \
//!   [--drop-rate 0.05 --drop-min-len 100] [--json report.json]
//!
//! FastSuite compatibility: `-o/--out1` is the canonical output flag,
//! `-j/--threads` is accepted (and ignored), `--json` writes the suite
//! report so fasthelper can aggregate the synthetic control like any
//! other sample.
//! ```

use anyhow::{bail, Result};

#[derive(Debug, Clone)]
pub struct Config {
    pub pairs: u64,
    pub profile: String,
    pub out: String,
    pub seed: u64,
    /// Read length. `None` (the default) takes it from the profile, so a
    /// profile can be fed straight in; an explicit value must match the profile.
    pub read_len: Option<u32>,
    pub insert_mean: f64,
    pub insert_sd: f64,
    pub insert_min: u32,
    pub gc: f64,
    pub reference: Option<String>,
    pub adapter_seq: Option<String>,
    pub gzip: bool,
    /// Fraction of pairs that re-sequence a previously emitted fragment
    /// (PCR duplicates). The mandril E-dataset measures 1-4%; the default
    /// sits in the middle of that range.
    pub dup_rate: f64,
    /// Ring-buffer capacity for the duplicate pool (fragments to copy).
    pub dup_pool: u32,
    /// Fraction of pairs whose true fragment is too short to survive
    /// trimming: insert < drop_min_len forces the GT keep window below
    /// min_length on both mates, so a correct trimmer MUST drop the pair.
    pub drop_rate: f64,
    /// Insert-length ceiling for the forced-drop class (bench min_length).
    pub drop_min_len: u32,
    /// Accepted for FastSuite CLI compatibility; the simulator is
    /// single-threaded, so any value > 1 is reported and ignored.
    pub threads: usize,
    /// Optional FastSuite report (schema v1) destination.
    pub json: Option<String>,
    /// Single-end mode: one read per fragment (no mate 2, no insert-size
    /// geometry; adapter readthrough still lands at the read's 3′ end).
    pub single: bool,
    /// Overwrite an existing output prefix instead of refusing.
    pub force: bool,
    /// Suppress the end-of-run summary line (warnings/errors still print).
    pub quiet: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            pairs: 0,
            profile: String::new(),
            out: String::new(),
            seed: 42,
            read_len: None,
            insert_mean: 300.0,
            insert_sd: 50.0,
            insert_min: 150,
            gc: 0.49,
            reference: None,
            adapter_seq: None,
            gzip: true,
            dup_rate: 0.02,
            dup_pool: 1 << 17,
            drop_rate: 0.05,
            drop_min_len: 100,
            threads: 1,
            json: None,
            single: false,
            force: false,
            quiet: false,
        }
    }
}

pub fn parse(args: &[String]) -> Result<Config> {
    let mut cfg = Config::default();
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
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        let next = |i: &mut usize| -> Result<String> {
            *i += 1;
            if *i >= args.len() {
                bail!("missing value after {arg}");
            }
            Ok(args[*i].clone())
        };
        match arg.as_str() {
            "--pairs" => cfg.pairs = next(&mut i)?.parse()?,
            "--profile" => cfg.profile = next(&mut i)?,
            // FastSuite canonical: -o/--out1 (a prefix for a simulator, not one file).
            "-o" | "--out" | "--out1" | "--out-prefix" => cfg.out = next(&mut i)?,
            "-j" | "--threads" => {
                cfg.threads = next(&mut i)?.parse()?;
            }
            "--json" => cfg.json = Some(next(&mut i)?),
            "--seed" => cfg.seed = next(&mut i)?.parse()?,
            "--read-len" => cfg.read_len = Some(next(&mut i)?.parse()?),
            "--insert-mean" => cfg.insert_mean = next(&mut i)?.parse()?,
            "--insert-sd" => cfg.insert_sd = next(&mut i)?.parse()?,
            "--insert-min" => cfg.insert_min = next(&mut i)?.parse()?,
            "--gc" => cfg.gc = next(&mut i)?.parse()?,
            "--reference" => cfg.reference = Some(next(&mut i)?),
            "--adapter-seq" => cfg.adapter_seq = Some(next(&mut i)?),
            "--no-adapter" => {
                cfg.adapter_seq = Some(String::new()); // empty = disabled
            }
            "--no-gzip" => cfg.gzip = false,
            "--single" => cfg.single = true,
            "--force" | "-f" => cfg.force = true,
            "-q" | "--quiet" => cfg.quiet = true,
            "--dup-rate" => cfg.dup_rate = next(&mut i)?.parse()?,
            "--dup-pool" => cfg.dup_pool = next(&mut i)?.parse()?,
            "--drop-rate" => cfg.drop_rate = next(&mut i)?.parse()?,
            "--drop-min-len" => cfg.drop_min_len = next(&mut i)?.parse()?,
            "-V" | "--version" => {
                println!("fastsim {}", env!("CARGO_PKG_VERSION"));
                std::process::exit(0);
            }
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => bail!("unknown argument {other:?} (see --help)"),
        }
        i += 1;
    }

    if cfg.threads > 1 {
        eprintln!(
            "fastsim: warning: simulation is single-threaded; --threads {} accepted for \
             FastSuite CLI compatibility and ignored",
            cfg.threads
        );
    }

    if cfg.pairs == 0 {
        bail!("--pairs must be >= 1");
    }
    if cfg.profile.is_empty() {
        bail!("--profile is required");
    }
    if cfg.out.is_empty() {
        bail!("--out is required");
    }
    if cfg.insert_sd <= 0.0 {
        bail!("--insert-sd must be > 0");
    }
    if !(0.0..=1.0).contains(&cfg.gc) {
        bail!("--gc must be in [0,1]");
    }
    if !(0.0..=1.0).contains(&cfg.dup_rate) {
        bail!("--dup-rate must be in [0,1]");
    }
    if !(0.0..=1.0).contains(&cfg.drop_rate) {
        bail!("--drop-rate must be in [0,1]");
    }
    if cfg.drop_min_len == 0 {
        bail!("--drop-min-len must be >= 1");
    }
    if let Some(l) = cfg.read_len {
        if cfg.insert_min > l {
            bail!("--insert-min ({}) should be <= read-len ({})", cfg.insert_min, l);
        }
    }
    Ok(cfg)
}

fn print_usage() {
    eprintln!(
        "fastsim — NovaSeq-like paired-end read simulator with ground truth

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
    -h, --help             this help"
    );
}
