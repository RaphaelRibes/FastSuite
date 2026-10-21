use crate::hasher::HashType;
use clap::{Parser, ValueEnum};

const ASCII_ART: &str = r#"
 /$$$$$$$$ /$$$$$$$                  /$$
| $$_____/| $$__  $$                | $$
| $$      | $$  \ $$  /$$$$$$   /$$$$$$$ /$$   /$$  /$$$$$$
| $$$$$   | $$  | $$ /$$__  $$ /$$__  $$| $$  | $$ /$$__  $$
| $$__/   | $$  | $$| $$$$$$$$| $$  | $$| $$  | $$| $$  \ $$
| $$      | $$  | $$| $$_____/| $$  | $$| $$  | $$| $$  | $$
| $$      | $$$$$$$/|  $$$$$$$|  $$$$$$$|  $$$$$$/| $$$$$$$/
|__/      |_______/  \_______/ \_______/ \______/ | $$____/
                                                  | $$
                                                  | $$
                                                  |__/
"#;

/// Command Line Interface (CLI) configuration for FDedup.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "A fast and memory-efficient FASTX PCR deduplication tool (Supports Paired-End)",
    before_help = ASCII_ART,
    arg_required_else_help = true
)]
#[command(help_expected = true)]
pub struct Cli {
    /// Path to the input FASTX file (R1 or Single-End). `-` reads stdin (D-4).
    #[arg(required = true, short = '1', long, visible_alias = "in1")]
    pub input: String,

    /// Path to the R2 input FASTX file (Optional, enables Paired-End mode)
    /// --output-r2 is required unless --dry-run is set.
    #[arg(short = '2', long, visible_alias = "in2")]
    pub input_r2: Option<String>,

    /// Path to the output file (R1 or Single-End). `-` writes stdout (D-4).
    /// Supported formats: .fastq, .fq, .fasta, .fa, .fna (+ .gz for compression)
    /// Optional when --dry-run is set.
    #[arg(short = 'o', long, visible_alias = "out1")]
    pub output: Option<String>,

    /// Path to the R2 output file (Required if --input-r2 is provided without --dry-run)
    /// Must have the same format as --output
    #[arg(short = 'p', long, visible_alias = "out2")]
    pub output_r2: Option<String>,

    /// Output record format when the output is a stream (`-`).
    /// Named paths keep inferring the format from their extension.
    #[arg(long, value_enum, default_value_t = OutFormat::Fastq)]
    pub out_format: OutFormat,

    /// Compression when the output is a stream (`-`).
    /// Named paths keep inferring compression from their extension.
    #[arg(long, value_enum, default_value_t = OutCompression::None)]
    pub out_compression: OutCompression,

    /// Force overwrite of output files if they exist
    #[arg(long, short)]
    pub force: bool,

    /// Enable verbose logs
    #[arg(long, short)]
    pub verbose: bool,

    /// Calculate the duplication rate without creating output files
    #[arg(long, short = 's')]
    pub dry_run: bool,

    /// Threshold for automatic hashing size selection (ignored if --hash is set)
    #[arg(long, short = 't', default_value_t = 0.001)]
    pub threshold: f64,

    /// Manually specify the hashing size (64 or 128 bits)
    #[arg(long, short = 'H')]
    pub hash: Option<HashMode>,

    /// Set the GZIP compression level (1-9).
    /// FastSuite canonical name: --gz-level (--compression kept as an alias).
    /// Note: levels < 6 disable multi-threaded gzip output (uses standard
    /// single-threaded flate2 instead).
    ///
    /// Measured on 5 M pairs at -j 8 (2026-10-21), because the trade is not
    /// obvious: level 1 is 1.16x FASTER (17.5 s vs 20.3 s) and uses 25x less
    /// peak memory (152 MB vs 3780 MB), but writes 1.87x more bytes (729 MB
    /// vs 390 MB). Both outputs are standard gzip. Level 1 is the right choice
    /// for a memory- or time-constrained run; level 6 is the right choice when
    /// output size matters.
    #[arg(
        long = "gz-level",
        visible_alias = "compression",
        short = 'c',
        default_value_t = 6,
        value_parser = clap::value_parser!(u32).range(1..=9)
    )]
    pub compression: u32,

    /// Expected read length in base pairs, used to tune I/O buffer sizes and capacity estimates.
    /// Does not need to be exact. Defaults to 150bp (Illumina short reads).
    /// Examples: 150 for Illumina, 10000 for ONT, 100000 for ultra-long ONT.
    #[arg(long, short = 'P', default_value_t = 150, value_parser = clap::value_parser!(usize))]
    pub read_length: usize,

    /// Enables abundance tracking mode: count occurrences of each unique sequence
    /// and annotate the output header with ;size=N.
    /// Note: incompatible with --threads (always runs single-threaded).
    #[arg(long)]
    pub abundance: bool,

    /// Separator used for abundance annotation in the sequence header (default: ";size=").
    /// Only meaningful when --abundance is set.
    #[arg(long, default_value = ";size=")]
    pub abundance_sep: String,

    /// Number of threads for parallel processing (0 = auto-detect). Default: 1 (single-threaded).
    /// When > 1, uses a producer/consumer pipeline for parallel hashing.
    #[arg(long, short = 'j', default_value_t = 1)]
    pub threads: usize,

    /// Pre-count the exact number of sequences by doing a fast first pass over the input.
    /// This gives the hash table an exact capacity, avoiding both OOM from over-estimation
    /// and reallocation overhead from under-estimation. Cost: one extra pass through the data.
    #[arg(long)]
    pub pre_count: bool,

    /// Directory for the temporary file used in abundance tracking.
    /// When set, abundance mode stores unique reads in an uncompressed temp file
    /// during the first pass and reads them back in the second pass (avoids
    /// re-decompressing the input gzip). Default: re-reads the input gzip.
    #[arg(long)]
    pub temp_dir: Option<String>,

    /// Write a machine-readable FastSuite report to PATH (report schema v1).
    /// --dry-run included: the report always carries the duplication statistics.
    #[arg(long, value_name = "PATH")]
    pub json: Option<String>,

    /// Suppress the end-of-run summary line.
    #[arg(long, short = 'q')]
    pub quiet: bool,
}

/// The selected hash mode (64-bit or 128-bit)
#[derive(ValueEnum, Clone, Debug)]
pub enum HashMode {
    #[value(name = "64")]
    Bit64,
    #[value(name = "128")]
    Bit128,
}

impl From<HashMode> for HashType {
    fn from(mode: HashMode) -> Self {
        match mode {
            HashMode::Bit64 => HashType::XXH3_64,
            HashMode::Bit128 => HashType::XXH3_128,
        }
    }
}

/// Record format for stream output (`-o -`): D-4 gives a stream no
/// extension to infer from, so the choice is explicit. Default: fastq.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutFormat {
    Fastq,
    Fasta,
}

/// Compression for stream output (`-o -`). Default: none — plain bytes in
/// the pipe; ask for `gz` when the consumer wants compressed bytes.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutCompression {
    None,
    Gz,
}
