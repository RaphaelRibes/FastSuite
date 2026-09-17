use clap::Parser;
use std::path::PathBuf;

/// Exact, one-pass FASTQ quality control.
///
/// Every number is an integer count read off the whole file: no sampling, no
/// hidden thresholds, no PASS/FAIL verdicts. Modules live in the report as
/// measurements; the decision is yours.
#[derive(Parser, Debug)]
#[command(name = "fastcheck", version, arg_required_else_help = true)]
pub struct Cli {
    /// Input FASTQ files (plain or .gz; BGZF is decoded in parallel).
    #[arg(required = true, value_name = "FILE")]
    pub files: Vec<PathBuf>,

    /// Write the machine-readable report (fastsuite.report.v1) to PATH.
    #[arg(long, value_name = "PATH")]
    pub json: Option<PathBuf>,

    /// Write the per-position matrices to PATH (TSV, 1-based positions).
    #[arg(long, value_name = "PATH")]
    pub positions_tsv: Option<PathBuf>,

    /// Worker threads (0 = all cores). Split across files, and across BGZF
    /// blocks within a file.
    #[arg(short = 'j', long, default_value_t = 0)]
    pub threads: usize,

    /// Stop after N records per file. The report declares the run partial.
    #[arg(long, value_name = "N")]
    pub reads: Option<u64>,

    /// Stop after N compressed bytes per file. The report declares the run
    /// partial.
    #[arg(long, value_name = "N")]
    pub bytes: Option<u64>,

    /// Quality encoding: 33, 64, or auto (auto = decide from the exact byte
    /// census, and say so when the census cannot prove it).
    #[arg(long, value_name = "33|64|auto", default_value = "auto")]
    pub phred_offset: String,

    /// Longest read included in the per-position matrices. Longer reads are
    /// still counted everywhere else and reported as reads_over_max_len.
    #[arg(long, default_value_t = 1000, value_name = "N")]
    pub max_read_len: usize,

    /// Estimate the insert-size distribution from paired-end overlap (requires
    /// exactly two input files). The per-pair measurement is exact; the pass
    /// samples, because a distribution needs statistics rather than every pair.
    /// Adds an `insert_hist` block to the report — this is the measurement a
    /// simulator needs, since fragment length is not recoverable from the reads.
    #[arg(long)]
    pub insert_hist: bool,

    /// Adapter sequence for --insert-hist, repeatable. With it, inserts at or
    /// below the read length are measured from the adapter boundary instead of
    /// being counted as unmeasurable — the regime where the reads run into
    /// adapter and the overlap is broken by it rather than by length.
    #[arg(long, value_name = "SEQ")]
    pub adapter_seq: Vec<String>,

    /// Pairs to sample for --insert-hist (default 20000).
    #[arg(long, default_value_t = 20000, value_name = "N")]
    pub insert_pairs: u64,

    /// Suppress the end-of-run summary line.
    #[arg(short = 'q', long)]
    pub quiet: bool,
}

impl Cli {
    pub fn requested_offset(&self) -> anyhow::Result<Option<u8>> {
        match self.phred_offset.as_str() {
            "auto" => Ok(None),
            "33" => Ok(Some(33)),
            "64" => Ok(Some(64)),
            other => anyhow::bail!("--phred-offset must be 33, 64 or auto (got {other:?})"),
        }
    }

    pub fn thread_count(&self) -> usize {
        if self.threads > 0 {
            self.threads
        } else {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        }
    }
}
