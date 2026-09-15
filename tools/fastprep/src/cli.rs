//! Command-line interface.
//!
//! Parses arguments, validates their mutual consistency, and builds the
//! two things the rest of the program consumes:
//!   * [`Mode`] — single-end vs paired-end, with all input/output paths.
//!   * Two [`Params`] structs — one per mate.
//!   * [`Runtime`] — I/O and threading knobs.

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::adapter::{self, Adapter, AdapterSet};
use crate::algorithm::{Direction, MatchRule, MinLen, Params};

// ---------------------------------------------------------------------------
// Adapter source enum (mutually exclusive on the CLI).
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, ValueEnum)]
pub enum AdapterPreset {
    /// TruSeq R1/R2 + Nextera.
    Illumina,
    /// Nextera/Transposase only.
    Nextera,
    /// Illumina Small RNA.
    SmallRna,
}

// ---------------------------------------------------------------------------
// CLI struct.
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(
    name = "fastprep",
    version,
    about = "Asymmetric quality-weighted polyX + adapter trimming for FASTQ",
    long_about = None,
)]
pub struct Cli {
    // ------------------------------------------------------------------ I/O
    /// Input FASTQ for R1 (single-end: the only input).
    /// Compression inferred from extension (.gz, .zst, .bz2, .xz).
    #[arg(short = '1', long = "in1", value_name = "PATH", help_heading = "Input / Output")]
    pub in1: PathBuf,

    /// Input FASTQ for R2 (triggers paired-end mode).
    #[arg(short = '2', long = "in2", value_name = "PATH", help_heading = "Input / Output")]
    pub in2: Option<PathBuf>,

    /// Output FASTQ for R1 (single-end: the only output).
    #[arg(short = 'o', long = "out1", value_name = "PATH", help_heading = "Input / Output")]
    pub out1: PathBuf,

    /// Output FASTQ for R2 (required in paired-end mode).
    #[arg(short = 'p', long = "out2", value_name = "PATH", help_heading = "Input / Output")]
    pub out2: Option<PathBuf>,

    /// Paired-end only: surviving R1 when R2 was dropped.
    #[arg(long, value_name = "PATH", help_heading = "Input / Output")]
    pub singletons1: Option<PathBuf>,

    /// Paired-end only: surviving R2 when R1 was dropped.
    #[arg(long, value_name = "PATH", help_heading = "Input / Output")]
    pub singletons2: Option<PathBuf>,

    /// Write per-read trim metadata to PATH (binary PXMT format).
    /// One record per touched read (trimmed or dropped).
    #[arg(long, value_name = "PATH", help_heading = "Input / Output")]
    pub metadata: Option<PathBuf>,

    /// Write machine-readable JSON statistics to PATH (in addition to stderr output).
    #[arg(long, value_name = "PATH", help_heading = "Input / Output")]
    pub json: Option<PathBuf>,

    // -------------------------------------------------------- Adapter trimming
    /// Built-in adapter preset. Mutually exclusive with --adapters.
    /// Values: illumina (TruSeq + Nextera), nextera, small-rna.
    /// Omit to auto-detect against the built-in adapter catalogue.
    #[arg(
        long = "adapter-preset",
        value_name = "NAME",
        conflicts_with = "adapters",
        help_heading = "Adapter Trimming"
    )]
    pub adapter_preset: Option<AdapterPreset>,

    /// FASTA file of adapter sequences. Mutually exclusive with --adapter-preset.
    /// Omit to auto-detect against the built-in adapter catalogue.
    #[arg(
        long = "adapters",
        value_name = "PATH",
        conflicts_with = "adapter_preset",
        help_heading = "Adapter Trimming"
    )]
    pub adapters: Option<PathBuf>,

    /// Disable adapter trimming entirely.
    #[arg(long = "no-adapter", help_heading = "Adapter Trimming")]
    pub no_adapter: bool,


    /// Maximum mismatch rate for adapter sequence matching [0, 1].
    #[arg(long, default_value_t = 0.1, value_name = "FLOAT", help_heading = "Adapter Trimming")]
    pub adapter_error_rate: f64,

    /// Minimum adapter overlap for SE / PE-fallback sequence matching.
    /// 8 (was 4): a 4-base window at 10% error tolerates 1 mismatch on 4
    /// bases — any read ending in "AGAT" (adapter prefix) was falsely
    /// trimmed (~10% of clean reads measured on the benchmark QC set).
    #[arg(long, default_value_t = 8, value_name = "N", help_heading = "Adapter Trimming")]
    pub min_adapter_overlap: usize,


    // -------------------------------------------------------- PolyX trimming
    /// Minimum post-trim read length. Integer N (bases) or fraction FL
    /// (e.g. `0.6L`). Reads below this after trimming are dropped.
    #[arg(long, value_name = "N|FL", default_value = "0.6L", help_heading = "PolyX Trimming")]
    pub min_len: MinLen,

    /// Mercy budget Δ_m: cumulative gap tolerated before declaring end-of-run.
    #[arg(long, default_value_t = 3.0, help_heading = "PolyX Trimming")]
    pub mercy: f64,

    /// Phred saturation threshold t.
    #[arg(long, default_value_t = 20, value_name = "N", help_heading = "PolyX Trimming")]
    pub phred_sat: u8,

    /// Quality-gate shape x (1 = linear, <1 = lenient, >1 = strict).
    #[arg(long, default_value_t = 1.0, value_name = "FLOAT", help_heading = "PolyX Trimming")]
    pub gate_shape: f64,

    /// Match leniency μ (0..=1).
    #[arg(long, default_value_t = 0.1, value_name = "FLOAT", help_heading = "PolyX Trimming")]
    pub match_leniency: f64,

    /// Mismatch strength ν (0..=1).
    #[arg(long, default_value_t = 0.5, value_name = "FLOAT", help_heading = "PolyX Trimming")]
    pub mismatch_strength: f64,

    /// Minimum homopolymer length m_h.
    #[arg(long, default_value_t = 6, value_name = "N", help_heading = "PolyX Trimming")]
    pub min_run: usize,

    /// Trim the NovaSeq first-cycle 5′ G-dropout: a short (1..=N) leading
    /// G-run with normal quality is removed exactly, no quality gate.
    /// Must be < --min-run so the scan stays in charge of long runs.
    #[arg(long, value_name = "N", help_heading = "PolyX Trimming")]
    pub head_dropout: Option<usize>,

    /// Target letter for anchored mode (repeatable, max 5 distinct).
    /// Omit for generic polyX (match on identical neighbours).
    #[arg(long, value_name = "LETTER", help_heading = "PolyX Trimming")]
    pub target: Vec<char>,

    /// Disable the primary (5′) end-scan.
    #[arg(long, help_heading = "PolyX Trimming")]
    pub no_5: bool,

    /// Disable the secondary (3′) end-scan.
    #[arg(long, help_heading = "PolyX Trimming")]
    pub no_3: bool,

    /// Disable the interior scan.
    #[arg(long, help_heading = "PolyX Trimming")]
    pub no_interior: bool,

    /// Single-end: scan 3′→5′ instead of 5′→3′.
    #[arg(long, conflicts_with_all = ["r1_scan_reverse", "r2_scan_reverse"], help_heading = "PolyX Trimming")]
    pub scan_reverse: bool,

    /// Paired-end: scan R1 3′→5′ instead of 5′→3′.
    #[arg(long, help_heading = "PolyX Trimming")]
    pub r1_scan_reverse: bool,

    /// Paired-end: scan R2 3′→5′ instead of 5′→3′. Recommended for
    /// Nextera/NextSeq/NovaSeq libraries (polyG at 5′ end of R2).
    #[arg(long, help_heading = "PolyX Trimming")]
    pub r2_scan_reverse: bool,

    // ------------------------------------------------------- Quality filtering
    /// Minimum whole-read mean Phred quality (raw, not ASCII offset).
    #[arg(long, value_name = "FLOAT", help_heading = "Quality Filtering")]
    pub mean_q: Option<f64>,

    /// ASCII offset used to decode FASTQ quality scores.
    /// 33 for modern Illumina; 64 for legacy Illumina 1.3/1.5.
    /// Omit to accept the default of 33, or use --auto-phred to detect.
    #[arg(long, value_name = "N", help_heading = "Quality Filtering",
          conflicts_with = "auto_phred")]
    pub phred_offset: Option<u8>,

    /// Auto-detect Phred quality offset from the first 1000 reads of the
    /// input. Mutually exclusive with --phred-offset.
    #[arg(long, help_heading = "Quality Filtering")]
    pub auto_phred: bool,

    // -----------------------------------------------------------------Runtime
    /// Unconditionally trim N bases from the 5′ end (before any scan).
    #[arg(long, value_name = "N", help_heading = "Runtime")]
    pub trim_5: Option<usize>,

    /// Unconditionally trim N bases from the 3′ end (before any scan).
    #[arg(long, value_name = "N", help_heading = "Runtime")]
    pub trim_3: Option<usize>,

    /// Number of worker threads (1 = serial).
    #[arg(short = 'j', long, default_value_t = 1, help_heading = "Runtime")]
    pub threads: usize,

    /// Records per chunk submitted to the worker pool.
    #[arg(long, default_value_t = 1024, help_heading = "Runtime")]
    pub chunk: usize,

    /// Compression level for gzip output (1..=9). FastSuite canonical name
    /// --gz-level; --compression is accepted as an alias.
    #[arg(
        long = "gz-level",
        visible_alias = "compression",
        default_value_t = 6,
        value_name = "N",
        help_heading = "Runtime"
    )]
    pub gz_level: u32,

    /// Compression level for zstd output (1..=22).
    #[arg(long, default_value_t = 3, value_name = "N", help_heading = "Runtime")]
    pub zst_level: i32,

    /// Remove exact duplicate reads before any trimming.
    ///
    /// Off by default: this is the only stage whose memory is not bounded by
    /// the chunk size. The fingerprint table grows with the number of UNIQUE
    /// reads, roughly 4.6 GB at 3e8 reads with 64-bit fingerprints. The
    /// projected size is checked against the enforced memory limit before
    /// anything is allocated.
    #[arg(long, help_heading = "Deduplication")]
    pub dedup: bool,

    /// Fingerprint width in bits. Wider means fewer collisions and twice the
    /// memory. Default picks the narrowest width whose collision probability
    /// stays under --dedup-threshold for the estimated read count.
    #[arg(long, value_name = "BITS", value_parser = ["64", "128"],
          help_heading = "Deduplication")]
    pub dedup_bits: Option<String>,

    /// Collision probability the automatic width selection must stay under.
    #[arg(long, default_value_t = 0.001, value_name = "P",
          help_heading = "Deduplication")]
    pub dedup_threshold: f64,

    /// Refuse to start if the projected fingerprint table exceeds this many
    /// megabytes. Defaults to the enforced memory limit of the process.
    #[arg(long, value_name = "MB", help_heading = "Deduplication")]
    pub dedup_max_mem: Option<u64>,

    /// Worker threads for gzip output. Defaults to half of --threads,
    /// minimum 1.
    ///
    /// Gzip output is the dominant cost on compressed workloads. Each output
    /// stream gets its own pool, so a paired run uses twice this number on
    /// top of the Rayon pool; the default leaves room for both inside a
    /// --threads-sized allocation.
    #[arg(long, value_name = "N", help_heading = "Runtime")]
    pub gz_threads: Option<usize>,
}

// ---------------------------------------------------------------------------
// Mode and Runtime.
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Mode {
    Single {
        in1: PathBuf,
        out1: PathBuf,
        params: Params,
    },
    Paired {
        in1: PathBuf,
        in2: PathBuf,
        out1: PathBuf,
        out2: PathBuf,
        singletons1: Option<PathBuf>,
        singletons2: Option<PathBuf>,
        params_r1: Params,
        params_r2: Params,
    },
}

/// How adapter sequences are supplied at runtime.
///
/// There is no auto-detect variant: detection resolves to a `Fixed` set in
/// `Cli::finalize`, so by the time the pipeline runs the adapters are known.
#[derive(Debug, Clone)]
pub enum AdapterSource {
    /// Pre-loaded set, from a preset, a FASTA, or catalogue detection.
    Fixed(Arc<AdapterSet>),
    /// Adapter trimming disabled.
    None,
}

/// Resolved deduplication settings.
#[derive(Debug, Clone, Copy)]
pub struct DedupConfig {
    /// Fingerprint width in bits: 64 or 128.
    pub bits: u32,
    /// Pre-sizing hint for the fingerprint table.
    pub capacity: usize,
    /// Expected read length, sampled from the first record. The anchor window
    /// is sized from this once for the run: deriving it per read would make
    /// its width vary with the polyG tail, which is the divergence it exists
    /// to absorb.
    pub read_len: usize,
    /// Resident-usage bound above which the table stops accepting inserts.
    /// See the saturation argument in the specification: freezing is
    /// one-sided, so it can miss duplicates but never drop a read that
    /// should have been kept.
    pub ceiling: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct Runtime {
    pub threads: usize,
    pub chunk: usize,
    pub compression: crate::io::CompressionLevels,
    pub phred_offset: u8,
    /// When `true`, `run_single`/`run_paired` auto-detect the Phred offset from
    /// the first 1000 reads before processing. `phred_offset` is ignored.
    pub phred_auto: bool,
    pub metadata: Option<PathBuf>,
    /// If `Some`, machine-readable JSON statistics are written here after
    /// processing, in addition to the human-readable stderr output.
    pub json: Option<PathBuf>,
    pub adapter_source: AdapterSource,

    /// Deduplication settings; `None` when `--dedup` was not given.
    pub dedup: Option<DedupConfig>,
    pub adapter_error_rate: f64,
    pub min_adapter_overlap_se: usize,
    /// Scoped Rayon thread pool for parallel chunk processing.
    /// `None` → use the global pool (configured via `rayon::build_global`).
    /// Set in tests to isolate parallelism per test case.
    pub thread_pool: Option<std::sync::Arc<rayon::ThreadPool>>,
}

// ---------------------------------------------------------------------------
// Validation and build.
// ---------------------------------------------------------------------------

/// Rough record count from file size, for pre-sizing the fingerprint table.
///
/// Deliberately conservative on compressed input: assuming a 4x ratio
/// under-estimates the count, so the table is under-sized and grows rather
/// than over-allocating up front. Getting this wrong costs a rehash, not
/// correctness.
fn estimate_records(path: &Path, read_len: usize) -> Result<usize> {
    let bytes = std::fs::metadata(path)
        .with_context(|| format!("sizing {}", path.display()))?
        .len();
    let per_record = (read_len * 2 + 50).max(1) as u64;
    let compressed = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| matches!(e, "gz" | "zst" | "zstd" | "bz2" | "xz"));
    let raw = if compressed { bytes.saturating_mul(4) } else { bytes };
    Ok((raw / per_record) as usize)
}

/// Read up to `n` sequences from a FASTX file, for adapter fingerprinting.
///
/// Sampling here rather than inside the pipeline costs one extra parse of
/// ~10k records, which is a few megabytes and page-cached by the time the
/// real pass reaches them. In exchange the pipeline needs no peek buffer, no
/// drain path, and no `AutoDetect` state: by the time it runs, the adapter
/// set is already fixed.
fn peek_sequences(path: &Path, n: usize) -> Result<Vec<Vec<u8>>> {
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("opening {} for adapter detection", path.display()))?;
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        match reader.next() {
            Some(rec) => out.push(
                rec.with_context(|| format!("reading {}", path.display()))?
                    .seq()
                    .to_vec(),
            ),
            None => break,
        }
    }
    Ok(out)
}

impl Cli {
    pub fn finalize(self) -> Result<(Mode, Runtime)> {
        // Paired-end consistency.
        let paired = self.in2.is_some();
        if paired && self.out2.is_none() {
            bail!("paired-end input (-2/--in2) requires matching -p/--out2");
        }
        if !paired && self.out2.is_some() {
            bail!("-p/--out2 requires matching -2/--in2");
        }
        if !paired && (self.singletons1.is_some() || self.singletons2.is_some()
            || self.r1_scan_reverse || self.r2_scan_reverse)
        {
            bail!("--singletons*, --r1-scan-reverse and --r2-scan-reverse apply to paired-end mode only");
        }

        // Numeric ranges.
        if !(0.0..=1.0).contains(&self.match_leniency) { bail!("--match-leniency must be in [0, 1]"); }
        if !(0.0..=1.0).contains(&self.mismatch_strength) { bail!("--mismatch-strength must be in [0, 1]"); }
        if self.mercy < 0.0 { bail!("--mercy must be ≥ 0"); }
        if self.gate_shape <= 0.0 { bail!("--gate-shape must be > 0"); }
        if self.phred_sat == 0 { bail!("--phred-sat must be > 0"); }
        if self.phred_sat > 93 {
            bail!("--phred-sat must be ≤ 93 (practical Phred+33 max; use ≤ 62 for Phred+64)");
        }
        let phred_offset = match self.phred_offset {
            None => 33u8, // default
            Some(v) => {
                if !matches!(v, 33 | 64) { bail!("--phred-offset must be 33 or 64"); }
                v
            }
        };
        if self.min_run < 2 { bail!("--min-run must be ≥ 2"); }
        if let Some(hd) = self.head_dropout {
            if hd < 1 { bail!("--head-dropout must be ≥ 1"); }
            if hd >= self.min_run {
                bail!("--head-dropout must be < --min-run (the scan owns runs ≥ min-run)");
            }
        }
        if self.threads == 0 { bail!("--threads must be ≥ 1"); }
        if self.chunk == 0 { bail!("--chunk must be ≥ 1"); }
        if !(1..=9).contains(&self.gz_level) { bail!("--gz-level must be in [1, 9]"); }
        if !(1..=22).contains(&self.zst_level) { bail!("--zst-level must be in [1, 22]"); }
        if !(0.0..1.0).contains(&self.dedup_threshold) {
            bail!("--dedup-threshold must be in [0, 1)");
        }
        if let Some(n) = self.gz_threads { if n == 0 { bail!("--gz-threads must be ≥ 1"); } }
        if let Some(mq) = self.mean_q { if mq < 0.0 { bail!("--mean-q must be ≥ 0"); } }
        if !(0.0..=1.0).contains(&self.adapter_error_rate) { bail!("--adapter-error-rate must be in [0, 1]"); }
        if self.min_adapter_overlap == 0 { bail!("--min-adapter-overlap must be ≥ 1"); }

        if self.target.len() > 5 { bail!("--target may be given at most 5 times"); }
        let rule = build_match_rule(&self.target)?;

        // Build adapter source.
        let adapter_source = if self.no_adapter {
            AdapterSource::None
        } else if let Some(preset) = self.adapter_preset {
            let adapters: Vec<Adapter> = match preset {
                AdapterPreset::Illumina => adapter::preset_illumina(),
                AdapterPreset::Nextera  => adapter::preset_nextera(),
                AdapterPreset::SmallRna => adapter::preset_small_rna(),
            };
            AdapterSource::Fixed(Arc::new(AdapterSet::new(
                adapters,
                self.adapter_error_rate,
                self.min_adapter_overlap,
            )))
        } else if let Some(path) = self.adapters {
            let adapters = adapter::load_fasta(&path)?;
            AdapterSource::Fixed(Arc::new(AdapterSet::new(
                adapters,
                self.adapter_error_rate,
                self.min_adapter_overlap,
            )))
        } else {
            // Auto-detect. The terminal k-mer is a fingerprint against the
            // built-in catalogue; the catalogue supplies the full sequence
            // that actually gets trimmed. The k-mer is never used as an
            // adapter itself: it is 12 bases, which both over-trims (a 12-mer
            // within 2 mismatches hits ~1.6% of clean reads by chance) and
            // under-trims (it is the END of the adapter fragment present, so
            // trimming from its position leaves the earlier bases behind).
            const PEEK: usize = 10_000;
            let mut seqs = peek_sequences(&self.in1, PEEK)?;
            if let Some(p2) = self.in2.as_ref() {
                // Both mates: TruSeq R1 and R2 diverge after a 13-base shared
                // prefix, so a fingerprint from R1's divergent tail would
                // leave R2 untrimmed.
                seqs.extend(peek_sequences(p2, PEEK)?);
            }
            let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
            let det = adapter::identify_from_reads(&refs);

            if det.adapters.is_empty() {
                bail!(
                    "adapter auto-detection found nothing in {} sampled reads.\n\
                     No sequence was guessed, because trimming on a guess is \
                     worse than not trimming. Pass one of:\n  \
                     --adapter-preset <illumina|nextera|small-rna>\n  \
                     --adapters <FASTA>   (sequence IDs become adapter names)\n  \
                     --no-adapter         (disable adapter trimming)",
                    det.sampled
                );
            }

            // A run commits every read to this one decision, so say what it
            // was and on what evidence.
            eprintln!("fastprep: adapter auto-detection ({} reads sampled)", det.sampled);
            for (kmer, count) in &det.evidence {
                eprintln!(
                    "  fingerprint {} in {} reads ({:.2}%)",
                    String::from_utf8_lossy(kmer),
                    count,
                    100.0 * *count as f64 / det.sampled.max(1) as f64
                );
            }
            for a in &det.adapters {
                eprintln!("  using {} ({} bp)", a.name, a.seq.len());
            }

            AdapterSource::Fixed(Arc::new(AdapterSet::new(
                det.adapters,
                self.adapter_error_rate,
                self.min_adapter_overlap,
            )))
        };

        // Build base Params.
        let mut base = Params {
            min_len: self.min_len,
            mercy: self.mercy,
            phred_sat: self.phred_sat,
            gate_shape: self.gate_shape,
            match_leniency: self.match_leniency,
            mismatch_strength: self.mismatch_strength,
            min_run: self.min_run,
            head_dropout: self.head_dropout,
            rule,
            scan_5: !self.no_5,
            scan_3: !self.no_3,
            scan_interior: !self.no_interior,
            direction: Direction::Forward,
            mean_q: self.mean_q,
            trim_5: self.trim_5,
            trim_3: self.trim_3,
            adapter: match &adapter_source {
                AdapterSource::Fixed(arc) => Some(Arc::clone(arc)),
                _ => None,
            },
            extra_committed: 0,
            original_len: None,
            gate_lut: [0.0; 256],
        };
        base.prepare();

        // Deduplication. Resolved here so the run refuses before allocating
        // rather than after filling the table.
        let dedup = if !self.dedup {
            None
        } else {
            let read_len = peek_sequences(&self.in1, 1)?
                .first()
                .map_or(150, |s| s.len());
            let est = estimate_records(&self.in1, read_len)?;
            let capacity = est.clamp(1, crate::memory::MAX_PREALLOC_ENTRIES);
            let bits: u32 = match self.dedup_bits.as_deref() {
                Some("64") => 64,
                Some("128") => 128,
                _ => crate::dedup::select_width(est as u64, self.dedup_threshold),
            };
            let proj = crate::memory::project(capacity, (bits / 8) as usize);
            let limit = crate::memory::discover();

            // An explicit --dedup-max-mem overrides discovery, but only
            // downwards: it cannot authorise exceeding a cgroup limit.
            let cap_bytes = self
                .dedup_max_mem
                .map(|mb| mb.saturating_mul(1024 * 1024))
                .into_iter()
                .chain(limit.as_ref().map(|l| l.bytes))
                .min();

            match crate::memory::admit(proj.bytes, limit.as_ref()) {
                crate::memory::Admission::Refuse { projected, limit, source } => bail!(
                    "deduplication needs ~{:.1} GB ({} buckets x {} B, {}-bit fingerprints)\n\
                     but this process is limited to {:.1} GB ({})\n\n  \
                     --dedup-bits 64      narrower fingerprints, half the memory\n  \
                     --no-dedup           run without deduplication",
                    crate::memory::mib(projected) / 1024.0,
                    proj.buckets,
                    bits / 8 + 1,
                    bits,
                    crate::memory::mib(limit) / 1024.0,
                    source,
                ),
                crate::memory::Admission::ProceedWithWarning { headroom_bytes } => {
                    eprintln!(
                        "fastprep: warning: deduplication table projected at {:.0} MB, \
                         leaving only {:.0} MB headroom",
                        crate::memory::mib(proj.bytes),
                        crate::memory::mib(headroom_bytes),
                    );
                }
                crate::memory::Admission::Proceed => {}
            }

            eprintln!(
                "fastprep: deduplication on, {}-bit fingerprints, table pre-sized for \
                 {} reads ({} buckets, {:.0} MB)",
                bits, capacity, proj.buckets, crate::memory::mib(proj.bytes),
            );

            Some(DedupConfig { bits, capacity, read_len, ceiling: cap_bytes })
        };

        let mode = if paired {
            let mut params_r1 = base.clone();
            params_r1.direction = if self.r1_scan_reverse { Direction::Reverse } else { Direction::Forward };
            let mut params_r2 = base;
            params_r2.direction = if self.r2_scan_reverse { Direction::Reverse } else { Direction::Forward };
            Mode::Paired {
                in1: self.in1,
                in2: self.in2.unwrap(),
                out1: self.out1,
                out2: self.out2.unwrap(),
                singletons1: self.singletons1,
                singletons2: self.singletons2,
                params_r1,
                params_r2,
            }
        } else {
            let mut params = base;
            params.direction = if self.scan_reverse { Direction::Reverse } else { Direction::Forward };
            Mode::Single { in1: self.in1, out1: self.out1, params }
        };

        let runtime = Runtime {
            threads: self.threads,
            chunk: self.chunk,
            compression: crate::io::CompressionLevels {
                gz: self.gz_level,
                zst: self.zst_level,
                // Half of --threads, so that a paired run's two compression
                // pools plus the Rayon pool stay inside the allocation.
                gz_threads: self.gz_threads.unwrap_or((self.threads / 2).max(1)),
            },
            phred_offset,
            phred_auto: self.auto_phred,
            metadata: self.metadata,
            json: self.json,
            adapter_source,
            dedup,
            adapter_error_rate: self.adapter_error_rate,
            min_adapter_overlap_se: self.min_adapter_overlap,
            thread_pool: None,
        };

        Ok((mode, runtime))
    }
}

fn build_match_rule(letters: &[char]) -> Result<MatchRule> {
    if letters.is_empty() {
        return Ok(MatchRule::Homopolymer);
    }
    let mut bytes: Vec<u8> = Vec::with_capacity(letters.len());
    for &c in letters {
        let up = c.to_ascii_uppercase();
        if !matches!(up, 'A' | 'C' | 'G' | 'T' | 'N') {
            bail!("--target {:?} is not one of A/C/G/T/N", c);
        }
        let b = up as u8;
        if bytes.contains(&b) { bail!("--target {:?} given more than once", c); }
        bytes.push(b);
    }
    Ok(MatchRule::anchored(&bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_rule_empty_is_homopolymer() {
        assert!(matches!(build_match_rule(&[]).unwrap(), MatchRule::Homopolymer));
    }

    #[test]
    fn build_rule_accepts_acgt() {
        for c in ['A', 'C', 'G', 'T', 'N', 'a', 'g'] {
            assert!(build_match_rule(&[c]).is_ok(), "should accept {}", c);
        }
    }

    #[test]
    fn build_rule_rejects_invalid() {
        assert!(build_match_rule(&['X']).is_err());
        assert!(build_match_rule(&['U']).is_err());
    }

    #[test]
    fn build_rule_rejects_duplicates() {
        assert!(build_match_rule(&['A', 'A']).is_err());
        assert!(build_match_rule(&['A', 'a']).is_err());
    }
}  
