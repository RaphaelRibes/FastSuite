//! CLI argument parsing (clap, FastSuite house style).
//!
//! Usage:
//! ```text
//! fastsim --pairs N --profile profiles/mandril_E.json -o prefix \
//!   [--seed S] [--read-len 150] [--insert-mean 300 --insert-sd 50] \
//!   [--gc G] [--reference ref.fa] [--adapter-seq SEQ] [--no-adapter] \
//!   [--no-gzip] [--gz-level 6] [--insert-min M] [--dup-rate 0.02] \
//!   [--drop-rate 0.05 --drop-min-len 100] [--json report.json] [--stats]
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
use clap::parser::ValueSource;
use clap::{CommandFactory, FromArgMatches, Parser};

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

    /// Insert length mean (default 300). Ignored when the profile carries a
    /// measured insert CDF.
    ///
    /// Inserts are measurable only up to 2 x read-len - 12 bp (288 bp at the
    /// default 150 bp read length): above that ceiling the mates share no
    /// 12-base anchor, so those lengths all land in ONE saturated histogram bin.
    /// There is no per-read censoring flag — a saturated bin is a measurement
    /// limit, and every length in it produces the same read (no overlap, no
    /// adapter).
    #[arg(long, default_value_t = 300.0)]
    pub insert_mean: f64,

    /// Insert length sd (default 50). Ignored when the profile carries a measured insert CDF.
    #[arg(long, default_value_t = 50.0)]
    pub insert_sd: f64,

    /// Insert length floor (default: min(150, the profile's read_len), so a
    /// short-read profile works unmodified). Ignored when the profile carries
    /// a measured insert CDF.
    #[arg(long)]
    pub insert_min: Option<u32>,

    /// GC fraction for random DNA (default: the profile's measured `gc`, a
    /// required profile key — there is no built-in fallback). An explicit
    /// --gc overrides the profile for a profile with no base plane.
    ///
    /// Ignored when the profile carries a per-position base plane, as every
    /// fitted profile does: the plane governs the composition position by
    /// position, so --gc has no effect, a warning says so, and the report
    /// publishes the profile's gc.
    #[arg(long)]
    pub gc: Option<f64>,

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

    /// DEFLATE level for the gzip outputs, 1..=9 (default 6: today's bytes).
    ///
    /// The level changes ONLY the encoding — the decompressed FASTQ is
    /// identical at every level. Deflate is ~75% of fastsim's wall time at
    /// level 6, and this tool's own backend (flate2/miniz_oxide) writes
    /// 38.75 MB per 200 000 pairs at level 6 and 52.63 MB at level 1 (+35.8%);
    /// the review's system-`gzip` proxy measured +17.9%, so size a campaign
    /// with the larger figure. End to end: 200 000 pairs go from ~14.8 s to
    /// ~2.3 s at level 1 (~6x), ~12 MB/s of packed input at 6 versus ~200 at 1.
    /// 4 is the balanced choice; 6 stays the default so existing campaigns
    /// keep their exact bytes.
    #[arg(
        long,
        value_name = "LEVEL",
        default_value_t = fastsim::fastq::DEFAULT_GZ_LEVEL,
        value_parser = clap::value_parser!(u32)
            .range(fastsim::fastq::MIN_GZ_LEVEL as i64..=fastsim::fastq::MAX_GZ_LEVEL as i64)
    )]
    pub gz_level: u32,

    /// Report instrument counters in the --json report: RNG words drawn,
    /// substitutions, mean Q and gzip bytes (FS-11). Counting is a no-op when
    /// this is absent, and it never changes the emitted bytes.
    #[arg(long)]
    pub stats: bool,

    /// PCR duplicate fraction (default 0.02; mandril measures 1-4%). Only the
    /// fallback: a profile carrying a measured dup_rate overrides it. A
    /// duplicate replays a previously emitted pair byte-for-byte.
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
    /// FASTSUITE_THREADS is used only when NO thread flag is given, in any
    /// spelling (-j 4, -j4, --threads 4, --threads=4).
    #[arg(short = 'j', long, default_value_t = 1)]
    pub threads: usize,

    /// Optional FastSuite report (schema v1) destination.
    ///
    /// `full_artefact_pairs` counts pairs where BOTH mates are entirely
    /// artefact — the number of GT rows whose `full_flag` is 1 — while
    /// `full_artefact_pairs_either` counts pairs where EITHER mate is.
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
    pub insert_min: Option<u32>,
    pub gc: Option<f64>,
    pub reference: Option<String>,
    pub adapter_seq: Option<String>,
    pub gzip: bool,
    pub gz_level: u32,
    pub stats: bool,
    pub dup_rate: f64,
    pub dup_pool: u32,
    pub drop_rate: f64,
    pub drop_min_len: u32,
    pub threads: usize,
    pub json: Option<String>,
    pub single: bool,
    pub force: bool,
    pub quiet: bool,
    /// Did the command line actually carry `--dup-rate`?
    ///
    /// clap hands over the RESOLVED value, so an explicit `--dup-rate 0.05` and
    /// the 0.02 default are indistinguishable in `Config`. A profile that
    /// measures its own rate overrides both, and `main` must warn only about a
    /// value the user actually typed (a default is inert by definition, and
    /// warning about it would fire on every fitted profile). `parse` fills these
    /// from clap's value sources; a programmatic `Config` leaves them false.
    pub dup_rate_given: bool,
    /// Same question for `--insert-mean` / `--insert-sd`, which a profile's
    /// measured insert CDF overrides.
    pub insert_mean_given: bool,
    pub insert_sd_given: bool,
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
            gz_level: cli.gz_level,
            stats: cli.stats,
            dup_rate: cli.dup_rate,
            dup_pool: cli.dup_pool,
            drop_rate: cli.drop_rate,
            drop_min_len: cli.drop_min_len,
            threads: cli.threads,
            json: cli.json,
            single: cli.single,
            force: cli.force,
            quiet: cli.quiet,
            // Filled by `parse` from clap's value sources: `From<Cli>` sees only
            // resolved values, where an explicit flag and its default are the
            // same thing.
            dup_rate_given: false,
            insert_mean_given: false,
            insert_sd_given: false,
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
        if self.dup_pool == 0 {
            bail!("--dup-pool must be >= 1");
        }
        if self.insert_sd <= 0.0 {
            bail!("--insert-sd must be > 0");
        }
        if let Some(gc) = self.gc {
            if !(0.0..=1.0).contains(&gc) {
                bail!("--gc must be in [0,1]");
            }
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
        // clap rejects an out-of-range --gz-level; this catches a Config built
        // programmatically (the tests do) so the encoder can never be handed a
        // level the DEFLATE backend would clamp silently.
        if !(fastsim::fastq::MIN_GZ_LEVEL..=fastsim::fastq::MAX_GZ_LEVEL).contains(&self.gz_level) {
            bail!(
                "--gz-level must be in {}..={}",
                fastsim::fastq::MIN_GZ_LEVEL,
                fastsim::fastq::MAX_GZ_LEVEL
            );
        }
        // --insert-min is validated in `main`, against the PROFILE's read_len:
        // the flag is ignored when the profile carries a measured insert CDF, and
        // only the profile can say whether it does. Checking it here as well would
        // refuse a run the profile makes valid (`--read-len 150 --insert-min 200`
        // with an insert CDF).
        Ok(())
    }
}

/// Is a thread-count flag present, in any spelling clap accepts?
///
/// `-j 4`, `-j4`, `-j=4`, `--threads 4` and `--threads=4` are the same request.
/// Recognising only the detached spellings let FASTSUITE_THREADS silently
/// override an explicit `-j4` / `--threads=4`.
fn thread_flag(args: &[String]) -> bool {
    args.iter().skip(1).any(|a| {
        a == "--threads"
            || a.starts_with("--threads=")
            // A short flag or cluster: only `-j` takes a value, so any `j` in a
            // short group is the thread flag (`-j4`, `-qj4`).
            || (a.starts_with('-') && !a.starts_with("--") && a[1..].contains('j'))
    })
}

/// Parse the process argv into a validated `Config`.
/// `-V/--version` and `-h/--help` print and exit 0, clap-style.
pub fn parse() -> Result<Config> {
    let matches = Cli::command().get_matches();
    let mut cli = Cli::from_arg_matches(&matches)?;
    // FASTSUITE_THREADS is the suite-wide default when NO thread flag is given
    // (CLI-SPEC §1). clap's `env` attribute can't express "only when the flag
    // is absent AND the value is valid", so resolve it manually.
    let args: Vec<String> = std::env::args().collect();
    if !thread_flag(&args) {
        if let Some(n) = std::env::var("FASTSUITE_THREADS")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
        {
            cli.threads = n;
        }
    }
    let mut cfg: Config = cli.into();
    // Which of the profile-overridden flags did the user actually type? The
    // resolved value cannot say (see `dup_rate_given`), so ask the matches.
    let given = |name: &str| matches.value_source(name) == Some(ValueSource::CommandLine);
    cfg.dup_rate_given = given("dup_rate");
    cfg.insert_mean_given = given("insert_mean");
    cfg.insert_sd_given = given("insert_sd");
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
