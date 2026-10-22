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
// PE-overlap mode (tri-state, clap ValueEnum).
// ---------------------------------------------------------------------------

/// When the paired-end overlap detector may pin the adapter boundary.
///
/// `auto` (the default) tracks *whether the user configured an adapter set*,
/// because the safety guard needs a set to check the boundary against: `auto`
/// runs only when the adapter source is EXPLICIT (`--adapter-preset` /
/// `--adapters`), and stays off with auto-detection or `--no-adapter`. A
/// detected set is not an explicit one — [`AdapterSource::Fixed`] carries the
/// distinction in `explicit` — so detection alone cannot switch the stage on;
/// the mode's measured value on the cohort is a null, and running it on every
/// pair for a set the user never named was the bug this field exists to stop.
/// `always` is the explicit opt-in to run regardless of the source. Either way
/// the detector refuses a boundary when no set is configured
/// (`detect_pe_overlap_trim(.., None)` is `None`), so `always` cannot trim
/// biology either — it only opts into running the detector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum PeOverlapMode {
    /// Never run the mode (equivalent to the flag being absent).
    Never,
    /// Run only when the adapter source is explicit (the default).
    Auto,
    /// Run whenever an overlap is detected, even with no adapter set.
    Always,
}

impl PeOverlapMode {
    /// Whether the mode may run at all, given the resolved adapter source.
    ///
    /// `Never` never runs; `Always` always runs (the caller has opted in);
    /// `Auto` runs only for an explicit source — the guard then has a set to
    /// check the boundary against, and the user asked for that set. A
    /// *detected* set has one too, but nobody asked to run the stage, so it
    /// does not count: measured against the cohort the mode is a null, and its
    /// per-pair cost (a `revcomp` allocation plus ~140 anchor probes) is not.
    /// This is the one place the tri-state is resolved, so the mode's reach
    /// cannot drift from the source that justifies it.
    pub fn runs_with(self, source: &AdapterSource) -> bool {
        match self {
            PeOverlapMode::Never => false,
            PeOverlapMode::Always => true,
            PeOverlapMode::Auto => matches!(source, AdapterSource::Fixed { explicit: true, .. }),
        }
    }

    /// Whether the mode runs, for callers that already know the source.
    pub fn runs(self) -> bool {
        !matches!(self, PeOverlapMode::Never)
    }
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
    /// Not required with --auto-params: a recommendation query has no run.
    // Option, not PathBuf: --auto-params is a query with no run, so in1/out1
    // are absent there. They are checked manually in finalize() AFTER the
    // auto-params refusal — a plain PathBuf is implicitly `required`, and
    // clap forbids `required` together with `required_unless_present`.
    #[arg(short = '1', long = "in1", value_name = "PATH", help_heading = "Input / Output")]
    pub in1: Option<PathBuf>,

    /// Input FASTQ for R2 (triggers paired-end mode).
    #[arg(short = '2', long = "in2", value_name = "PATH", help_heading = "Input / Output")]
    pub in2: Option<PathBuf>,

    /// Output FASTQ for R1 (single-end: the only output).
    /// Not required with --auto-params: a recommendation query has no run.
    #[arg(short = 'o', long = "out1", value_name = "PATH", help_heading = "Input / Output")]
    pub out1: Option<PathBuf>,

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

    /// Treat an assumption as fatal. By default, when a detection returns
    /// nothing (e.g. no adapter found in the sampled reads) fastprep takes the
    /// safe action, records the assumption and carries on; with --strict it
    /// refuses to run instead.
    #[arg(long, help_heading = "Runtime")]
    pub strict: bool,

    /// Overwrite an existing -o/-p (and --metadata/--json/--singletons*)
    /// instead of refusing. Without it fastprep exits 1 naming the file
    /// that is in the way, so a re-run cannot silently destroy a finished
    /// result.
    #[arg(long, short = 'f', help_heading = "Runtime")]
    pub force: bool,

    /// Compute every counter, read the input, and write nothing: no output
    /// FASTQ, no metadata, no --json report. The summary still goes to
    /// stderr, so a dry run is exactly what a real run would have reported.
    #[arg(long, help_heading = "Runtime")]
    pub dry_run: bool,

    /// Suppress the progress and end-of-run summary stream. Warnings and
    /// errors still print: -q never hides a reason a run went wrong.
    #[arg(long, short = 'q', help_heading = "Runtime")]
    pub quiet: bool,

    // --------------------------------------------------- Model recommendation
    /// Recommend the trimming config from a fastcheck QC report of the input
    /// (`fastcheck --json qc.json R1 R2`). Requires --print-args.
    ///
    /// REFUSES until a Phase D model file exists. The shipped default is
    /// NOT the swept winner on the current corpus: over 101 real ENA
    /// datasets it is the per-dataset winner on 1 of 101, and a per-dataset
    /// winner from the 55-config grid beats it on the other 100, so there
    /// is real signal to switch on but no model trained on a corpus large
    /// enough to beat the default on held-out studies (docs/AUTO-PARAMS.md).
    /// `fastprep --auto-params <qc.json> --print-args` — the flag parses,
    /// validates the QC path, and refuses with a message pointing at the
    /// doc, exactly as CLI-SPEC §1 rule 4 requires of a correct refusal.
    #[arg(long, value_name = "QC.json", help_heading = "Model Recommendation")]
    pub auto_params: Option<PathBuf>,

    /// With --auto-params: print the recommended argv instead of running.
    /// The only mode --auto-params supports today, and the mode it will
    /// always support: the recommendation is inspectable before it is used.
    #[arg(long, requires = "auto_params", help_heading = "Model Recommendation")]
    pub print_args: bool,


    /// Maximum mismatch rate for adapter sequence matching [0, 1].
    #[arg(long, default_value_t = 0.2, value_name = "FLOAT", help_heading = "Adapter Trimming")]
    pub adapter_error_rate: f64,

    /// Paired-end only: infer the adapter boundary from the mate pair's own
    /// overlap, and trim from there.
    ///
    /// Values:
    /// * `never`  — the mode never runs.
    /// * `auto`   — **the default.** Runs only when an explicit adapter set is
    ///   configured (`--adapter-preset` or `--adapters`). Off with
    ///   auto-detection or `--no-adapter`, because the safety guard needs a set
    ///   to check the boundary against.
    /// * `always` — runs whenever an overlap is detected, even with no adapter
    ///   set (the caller explicitly opts in to the risk).
    ///
    /// The bare spelling `--pe-overlap` (no value) is accepted as an alias for
    /// `always`, so existing invocations that asked for the mode keep working.
    ///
    /// When the insert is shorter than the read, R1's 3′ suffix and the
    /// reverse-complement prefix of R2 describe the same fragment, so the mate
    /// pair pins the insert size and hence where the adapter readthrough begins
    /// on each mate. The detector anchors on a 12-base match (the same value
    /// fastcheck's insert-size estimator uses), verifies the whole shared
    /// window base for base, and refuses anything it cannot pin uniquely: no
    /// overlap, no anchor, a chance anchor hit, or more than one admissible
    /// alignment all fall back to the ordinary per-mate behaviour, silently and
    /// exactly. A boundary is applied **only if** the readthrough past it is
    /// adapter sequence from the configured set (`readthrough_is_adapter`), so
    /// an overlap alone never trims biology; with no set configured the mode
    /// refuses a boundary outright, structurally.
    ///
    /// Measured as a null on a cohort built to be scorable (Cell B,
    /// docs/OVERTRIM-PAIRED.md): exact matches identical with and without the
    /// mode, overtrim +106, +0.000 pp against a pre-registered 0.5 pp threshold.
    /// The 102-dataset corpus cannot answer the question -- adapter readthrough
    /// is configured in 53 of its 102 profiles but materialises in only 11,395
    /// of 10.2M pairs (0.112%), because it needs an insert shorter than the read
    /// -- so the corpus serves as the guard cell, not the value test. See
    /// docs/plans/2026-10-19_004500-pe-overlap-with-without.md and
    /// bench/scripts/trim_report.py's FASTPREP_VARIANTS slot (`fp-pe`).
    #[arg(long, value_name = "MODE", value_enum, default_value_t = PeOverlapMode::Auto,
          num_args = 0..=1, default_missing_value = "always",
          help_heading = "Adapter Trimming")]
    pub pe_overlap: PeOverlapMode,

    /// Minimum adapter overlap for sequence matching.
    /// 14 at the shipped `--adapter-error-rate 0.2`: a 4-base window admits a
    /// mismatch on 4 bases, so any read ending in "AGAT" (adapter prefix) was
    /// falsely trimmed (~10% of clean reads measured on the benchmark QC set).
    #[arg(long, default_value_t = 14, value_name = "N", help_heading = "Adapter Trimming")]
    pub min_adapter_overlap: usize,

    /// How much evidence an adapter trim requires.
    ///
    /// `floor` (the default) is today's behaviour: a 3'-anchored match of at
    /// least `--min-adapter-overlap` bases. `boundary` also accepts a shorter
    /// tail when the mate geometry has already pinned the boundary — the tail
    /// must still BE the adapter prefix at that boundary, only the length
    /// floor is lifted.
    ///
    /// This exists because the PE-overlap detector pins the insert boundary
    /// exactly and then checks it with the general matcher, which applies the
    /// overlap floor to the tail: a 10-base tail is provably unmatchable
    /// there, so the boundary is discarded. Measured on the paired A/B
    /// (`bench/scripts/undertrim_ab.py --paired`), the PE path currently buys
    /// 0.06 pp of residual-adapter reduction because of that.
    #[arg(long, value_name = "POLICY", value_enum,
          default_value_t = crate::adapter::AdapterEvidence::Floor,
          help_heading = "Adapter Trimming")]
    pub adapter_evidence: crate::adapter::AdapterEvidence,

    /// Take the smallest insert when the mate overlap is ambiguous.
    ///
    /// The overlap detector refuses an alignment that verifies at more than one
    /// insert size (a repeated fragment): choosing would be a guess. With this
    /// flag it chooses anyway — the SMALLEST verifying insert, i.e. the largest
    /// readthrough — but only when the adapter-checked boundary guard confirms
    /// that adapter starts exactly behind it on the pinned mate. Any other
    /// ambiguous pair keeps today's refusal.
    ///
    /// The guard still needs an adapter set: with none, the boundary is refused
    /// outright whatever this flag says (DECISIONS.md D-14).
    ///
    /// The default is off (ambiguity refuses).
    #[arg(long, help_heading = "Adapter Trimming")]
    pub pe_ambiguous_smallest: bool,

    /// Retry the adapter match before a terminal no-signal G-run.
    ///
    /// NovaSeq's two-colour overprint replaces the last bases of a read with
    /// G. The adapter matcher is 3'-anchored, so those G's sit inside the
    /// alignment window and spend the mismatch budget — with a 20-base
    /// adapter tail whose last 6 bases are G, the match fails outright
    /// (0 of 2 000 exact cuts, against 2 000 of 2 000 for 4 or fewer).
    ///
    /// With this flag, a read whose 3' end carries a G-run of at least N
    /// bases is re-matched on the slice before the run, and the G-run is
    /// consumed with the adapter (the polyX stage removes it anyway). The
    /// full overlap floor still applies to the adapter bases themselves.
    ///
    /// The default is off (no retry).
    #[arg(long, value_name = "N", help_heading = "Adapter Trimming")]
    pub adapter_g_run_retry: Option<usize>,

    /// Retry the adapter match with quality- and N-aware mismatch scoring.
    ///
    /// The matcher charges an `N` a full mismatch, which spends the mismatch
    /// budget on a base the sequencer never called. With this flag, a read
    /// the ordinary matcher rejected is re-scored: an ambiguous base (`N`) or
    /// a base called at quality <= N costs HALF a mismatch instead of one.
    ///
    /// A window with fewer than 12 real (unambiguous) calls is refused
    /// outright, so a read of N's cannot match anything.
    ///
    /// The default is off (no retry).
    #[arg(long, value_name = "LOW_Q", help_heading = "Adapter Trimming")]
    pub adapter_scored_retry: Option<u8>,

    /// Trim corroborated short adapter tails (6-13 bases).
    ///
    /// Tails shorter than `--min-adapter-overlap` are refused by the ordinary
    /// matcher before any mismatch count is consulted, and most of the
    /// remaining undertrimming lives there. With this flag, an EXACT adapter
    /// prefix of at least MIN bases at the 3' end is trimmed when the quality
    /// drops sharply at the boundary: the bases just before it are at least
    /// `--adapter-q-drop` lower in mean quality than the bases before those.
    ///
    /// Exactness plus the quality cliff is deliberate. Lowering the overlap
    /// floor without corroboration is the `mo=6` arm: 1.29 % false trims on
    /// clean reads, a 180x overtrim regression.
    ///
    /// The default is off.
    #[arg(long, value_name = "MIN", help_heading = "Adapter Trimming")]
    pub adapter_short_tail: Option<usize>,

    /// The minimum quality drop required for a short-tail trim.
    ///
    /// The mean quality of the 4 bases before the candidate boundary must be
    /// at least this far below the mean quality of the 10 bases before
    /// those. Higher is more conservative.
    #[arg(long, value_name = "Q", default_value_t = 8,
          requires = "adapter_short_tail", help_heading = "Adapter Trimming")]
    pub adapter_q_drop: u8,

    /// The dataset-level prefix enrichment a short tail must clear.
    ///
    /// The expected number of sampled reads ending with a given 6-base
    /// adapter prefix is `n * 4^-6`; the observed count divided by that is
    /// the enrichment. A factor of 10 means the prefix occurs ten times
    /// more often than chance predicts.
    ///
    /// This is the corroboration signal that works on data whose qualities
    /// are flat, where the quality cliff cannot be seen.
    ///
    /// The default is off.
    #[arg(long, value_name = "FACTOR", default_value_t = 0,
          requires = "adapter_short_tail", help_heading = "Adapter Trimming")]
    pub adapter_min_enrichment: u32,

    /// Prefer the LONGEST corroborated short-tail candidate instead of the
    /// shortest.
    ///
    /// `--adapter-short-tail` enumerates nested EXACT adapter prefixes and
    /// trims the shortest corroborated one: the smaller claim about how much
    /// of the read is adapter. With this flag the longest corroborated
    /// candidate wins instead, so a tail whose whole length is adapter is
    /// removed whole rather than down to the first corroborated prefix.
    ///
    /// Which bias is right is a property of the library: shortest-wins is the
    /// conservative one where trimming biology costs more than leaving
    /// adapter, longest-wins where residual adapter costs more than the
    /// residual chance of an exact longer prefix. Neither mode accepts an
    /// uncorroborated candidate.
    ///
    /// The default is off (shortest wins).
    #[arg(long, requires = "adapter_short_tail", help_heading = "Adapter Trimming")]
    pub short_tail_longest: bool,


    // -------------------------------------------------------- PolyX trimming
    /// Minimum post-trim read length. Integer N (bases) or fraction FL of the
    /// INSERT length (the read after adapter readthrough is removed; the
    /// original read length when no adapter is trimmed). A read whose
    /// post-trim length is below this floor is dropped; a read that would need
    /// more trim than the floor allows is dropped rather than trimmed past it.
    /// Integer floors are absolute: `--min-len 100` drops a 40 bp read.
    #[arg(long, value_name = "N|FL", default_value = "0.9L", help_heading = "PolyX Trimming")]
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
    #[arg(long, default_value_t = 12, value_name = "N", help_heading = "PolyX Trimming")]
    pub min_run: usize,

    /// Trim the NovaSeq first-cycle 5′ G-dropout: a short (1..=N) leading
    /// G-run with normal quality is removed exactly, no quality gate.
    /// Must be < --min-run so the scan stays in charge of long runs.
    ///
    /// `4` IS the default. The rule is a prior, not a detector — a leading
    /// biological G is indistinguishable from a 1 bp dropout by base identity —
    /// so trimming it is wrong on roughly three times in four affected reads.
    /// Whether that makes the rule wrong depends on the PRICE of the two
    /// errors, not on the rate: with undertrim costed at r x overtrim, the rule
    /// pays whenever r > P(biology)/P(artefact) = 2.6, measured three ways
    /// (posterior 0.72/0.28; the 102-library table 21.94/8.33; the 101-library
    /// axis sweep 3.70). At 7:1 — the standing preference "prefer less data to
    /// wrong data" — the shipped 4 scores +39.1 of mates against +16.4 for 0,
    /// i.e. a 22.7-point margin; at 1:1 (total error) the ordering flips and 0
    /// is the better setting (7.44% against 24.64% total error). Both readings
    /// are true; the default encodes the 7:1 one. `--head-dropout 0` disables
    /// the rule for a pipeline that prices a lost base and a kept artefact
    /// equally. See docs/DECISIONS.md D-16, docs/TRIM-ACCURACY.md and
    /// bench/results/phase-d/head-dropout-axis.tsv.
    #[arg(long, value_name = "N", default_value = "4", help_heading = "PolyX Trimming")]
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

    /// Disable the interior scan (the shipped default: the interior scan is OFF).
    #[arg(long, default_value_t = true, num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set, help_heading = "PolyX Trimming")]
    pub no_interior: bool,

    /// Re-enable the interior scan (off by default in the shipped config).
    #[arg(long = "interior", overrides_with = "no_interior", help_heading = "PolyX Trimming")]
    pub interior: bool,

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

    /// Boundary-anchored cut: place the polyX cut at the homopolymer run
    /// boundary (the last matching base) instead of the mercy walk's end.
    /// The mercy walk absorbs ~2 flanking non-run bases past a clean run
    /// before Δ_m trips, so the whole-run cut overshoots the run by that much.
    /// Anchoring trims exactly the run extent, dropping the overshoot.
    ///
    /// **Default ON since 0.4.0.** Measured on 102 datasets at the 7:1 trim
    /// cost anchor: +0.0121 per mate, 0.32 pp of mates moved from overtrim to
    /// exact, no undertrim cost (`docs/plans/2026-10-21_fastprep-exactness.md`
    /// §7). `--boundary-anchored` is therefore accepted as a no-op for existing
    /// invocations; `--boundary-anchored=false` (or `--no-boundary-anchored`)
    /// restores the pre-0.4.0 cut.
    #[arg(long, num_args = 0..=1, default_missing_value = "true",
          value_name = "BOOL", help_heading = "PolyX Trimming")]
    pub boundary_anchored: Option<bool>,

    /// Restore the pre-0.4.0 whole-run cut, including the mercy walk's ~2 bp
    /// flank overshoot.
    #[arg(long, help_heading = "PolyX Trimming")]
    pub no_boundary_anchored: bool,

    /// After a qualifying polyX end cut, extend it by up to N bases toward the
    /// read interior while each extension base is the run's own base (case
    /// insensitive) AND is called below --phred-sat.
    ///
    /// The shipped mercy cut already overshoots a clean run by ~2 flanking
    /// bases, so on a read whose run is followed by more of the same base at
    /// low quality the cut stops short of the artefact evidence. This extends
    /// it, one base at a time, over nothing but same-base low-quality calls.
    ///
    /// Accepted values are 0, 1 and 2. `0` (the default) disables the
    /// extension, leaving every existing behaviour unchanged. An interior run
    /// is never touched, a cut of 0 stays 0, `N` is never the run base nor an
    /// extension base, and an extension that would push the emitted read under
    /// --min-len drops the read (DropReason::Length) rather than violating the
    /// floor.
    #[arg(long, value_name = "N", default_value_t = 0,
          value_parser = parse_overtrim_margin, help_heading = "PolyX Trimming")]
    pub overtrim_margin: usize,

    // ------------------------------------------------------- Quality filtering
    /// Minimum whole-read mean Phred quality (raw, not ASCII offset).
    #[arg(long, value_name = "FLOAT", default_value = "25", help_heading = "Quality Filtering")]
    pub mean_q: Option<f64>,

    /// Disable the mean-quality filter (the shipped default is 25).
    #[arg(long = "no-mean-q", help_heading = "Quality Filtering")]
    pub no_mean_q: bool,

    /// Succeed even when the run drops every read.
    ///
    /// By default a run that outputs zero reads exits non-zero: an empty
    /// output file plus exit 0 is the one combination a pipeline cannot
    /// distinguish from success, and every cause of a total drop (a Phred
    /// offset read the wrong way, a floor above the read length, a mean-Q gate
    /// on a degraded library) is a configuration error worth stopping for.
    #[arg(long = "allow-total-drop", help_heading = "Quality Filtering")]
    pub allow_total_drop: bool,

    /// ASCII offset used to decode FASTQ quality scores.
    /// 33 for modern Illumina; 64 for legacy Illumina 1.3/1.5.
    /// Omit to accept the default of 33, or use --auto-phred to census the
    /// input (see the assumption the report records).
    #[arg(long, value_name = "N", help_heading = "Quality Filtering",
          conflicts_with = "auto_phred")]
    pub phred_offset: Option<u8>,

    /// Census the first 1000 records of the input for a quality byte below
    /// `'@'` (64), which would prove Phred+33. Bytes at or above 64 prove
    /// nothing (a Phred+33 Q>=31 library looks identical), so the offset stays
    /// at 33 and the report records that as an assumption, with the census as
    /// its basis. Mutually exclusive with --phred-offset.
    #[arg(long, default_value_t = true, num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set, help_heading = "Quality Filtering")]
    pub auto_phred: bool,

    /// Disable the Phred-offset census (the shipped default is on).
    #[arg(long = "no-auto-phred", help_heading = "Quality Filtering")]
    pub no_auto_phred: bool,

    // -----------------------------------------------------------------Runtime
    /// Unconditionally trim N bases from the 5′ end (before any scan).
    #[arg(long, value_name = "N", help_heading = "Runtime")]
    pub trim_5: Option<usize>,

    /// Unconditionally trim N bases from the 3′ end (before any scan).
    #[arg(long, value_name = "N", help_heading = "Runtime")]
    pub trim_3: Option<usize>,

    /// Number of worker threads (1 = serial).
    ///
    /// Throughput saturates around 8 workers: measured on a 200k-pair run,
    /// `-j 16` cost 12.6% more CPU than `-j 8` (5.631 vs 5.004 cpu-s) for 0.2%
    /// less wall time (0.651 vs 0.649 s). Past the plateau, extra threads buy
    /// nothing and burn CPU, so `-j 8` is the sensible ceiling on a typical box.
    /// The report's `busy_cores` metric (CPU-seconds per wall-second) shows
    /// where a run actually landed; a value near 1 means it ran effectively
    /// serial whatever `-j` said.
    #[arg(short = 'j', long, default_value_t = 1, help_heading = "Runtime")]
    pub threads: usize,

    /// Records per chunk submitted to the worker pool.
    #[arg(long, default_value_t = 1024, help_heading = "Runtime")]
    pub chunk: usize,

    /// Compression level for gzip output (1..=9). FastSuite canonical name
    /// --gz-level; --compression is accepted as an alias.
    ///
    /// Measured on 5 M pairs at -j 8 (2026-10-21): level 1 is 2.07x faster
    /// than level 6 (3.81 s vs 7.90 s) for 8% more bytes (441 MB vs 407 MB),
    /// at the same peak RSS. Unlike fastdedup, fastprep's parallel BGZF writer
    /// has no level floor, so level 1 keeps full parallel compression.
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
    #[arg(long, default_value_t = true, num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set, help_heading = "Deduplication")]
    pub dedup: bool,

    /// Disable deduplication (the shipped default: on).
    #[arg(long = "no-dedup", overrides_with = "dedup", help_heading = "Deduplication")]
    pub no_dedup: bool,

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
/// `explicit` is what distinguishes the two ways a set can come into being,
/// because they are not interchangeable for decisions that need a *user's*
/// configuration rather than any set at all: a preset or a FASTA is something
/// the user asked for, a catalogue detection is the tool's own guess. Only
/// `--adapter-preset` / `--adapters` set it to `true`.
#[derive(Debug, Clone)]
pub enum AdapterSource {
    /// Pre-loaded set, from a preset, a FASTA, or catalogue detection.
    Fixed {
        set: Arc<AdapterSet>,
        /// `true` only when the user named the set on the command line. A
        /// detected set is `false`, which is what keeps `--pe-overlap auto`
        /// from running a stage the user never asked for (see
        /// [`PeOverlapMode::runs_with`]).
        explicit: bool,
    },
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
    /// `--dry-run`: compute every counter and write nothing at all. The
    /// summary still goes to stderr, so a dry run reports exactly what a
    /// real run would have.
    pub dry_run: bool,
    /// `-q/--quiet`: suppress the progress and end-of-run summary stream.
    /// Warnings and errors still print — quiet never hides a reason a run
    /// went wrong (CLI-SPEC §1 rule 5).
    pub quiet: bool,
    /// `--allow-total-drop`: exit 0 even when no read survives.
    pub allow_total_drop: bool,
    pub adapter_source: AdapterSource,
    /// `--pe-overlap`: the tri-state as resolved in `finalize()`. `Auto` is
    /// resolved against the adapter source there (explicit `Fixed` set ⇒
    /// [`PeOverlapMode::Auto`] stays; detected set or `--no-adapter` ⇒
    /// [`PeOverlapMode::Never`]), so the pipeline reads the mode directly and
    /// the mode's reach cannot drift from the source that justifies it.
    pub pe_overlap: PeOverlapMode,
    /// `--adapter-evidence`: what the PE-overlap guard accepts once the
    /// boundary is pinned. `Floor` is today's behaviour; `Boundary` also
    /// accepts a tail shorter than `min_adapter_overlap_se`, because the mate
    /// geometry already fixed the boundary.
    pub adapter_evidence: crate::adapter::AdapterEvidence,
    /// `--pe-ambiguous-smallest`: let the overlap detector resolve an ambiguous
    /// alignment to the smallest verifying insert, when the boundary guard
    /// confirms adapter behind it. `false` (the default) keeps the refusal.
    pub pe_ambiguous_smallest: bool,

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

/// How many reads are sampled for the dataset-level enrichment signal.
const SAMPLE: usize = 10_000;

/// Peek a sample of reads for the dataset-level enrichment signal.
///
/// Both mates are sampled, because TruSeq R1 and R2 diverge after a 13-base
/// shared prefix: a fingerprint from R1's divergent tail would leave R2
/// uncorroborated.
fn peek_for_enrichment(
    in1: &Path,
    in2: Option<&Path>,
    n: usize,
) -> Result<Vec<Vec<u8>>> {
    let mut seqs = peek_sequences(in1, n)?;
    if let Some(p2) = in2 {
        seqs.extend(peek_sequences(p2, n)?);
    }
    Ok(seqs)
}

impl Cli {
    /// Every path this run would create or replace, in the order they are
    /// committed. Used by the `--force` clobber check: refusing on the first
    /// one that exists is what makes the message name a real file.
    fn out_paths(&self) -> Vec<PathBuf> {
        let mut v = Vec::new();
        if let Some(p) = &self.out1 { v.push(p.clone()); }
        if let Some(p) = &self.out2 { v.push(p.clone()); }
        if let Some(p) = &self.singletons1 { v.push(p.clone()); }
        if let Some(p) = &self.singletons2 { v.push(p.clone()); }
        if let Some(p) = &self.metadata { v.push(p.clone()); }
        if let Some(p) = &self.json { v.push(p.clone()); }
        v
    }

    /// --auto-params: the honest refusal, checked before anything else runs.
    ///
    /// The model does not exist yet (docs/AUTO-PARAMS.md), so the flag cannot
    /// recommend anything. What it CAN do is fail correctly and early: refuse
    /// when the QC report is missing or unparseable, and refuse with a
    /// message that names the reason and points at the doc. CLI-SPEC §1 rule
    /// 4: a correct refusal is a non-zero exit with a message saying what to
    /// do. When a Phase D model file ships, this becomes the tree-walk over
    /// it — the checks above (path exists, parses, has `metrics`) are the same
    /// ones the walk needs, so they stay.
    fn validate_auto_params(&self) -> Result<()> {
        let Some(qc) = self.auto_params.as_ref() else {
            // --print-args without --auto-params: clap's `requires` already
            // catches this; the guard is for safety in case that changes.
            if self.print_args {
                bail!("--print-args requires --auto-params <QC.json>");
            }
            return Ok(());
        };
        if !self.print_args {
            bail!(
                "--auto-params requires --print-args: the recommendation is \\
                 inspectable before it is used, so there is no silent-apply mode"
            );
        }
        let text = std::fs::read_to_string(qc)
            .with_context(|| format!("reading QC report {}", qc.display()))?;
        let doc: serde_json::Value = serde_json::from_str(&text)
            .with_context(|| format!("{} is not valid JSON", qc.display()))?;
        let metrics = doc.get("metrics").and_then(|m| m.as_object());
        if metrics.is_none() {
            bail!(
                "{} has no `metrics` block — pass a fastcheck QC report \\
                 (`fastcheck --json qc.json R1 R2`), not an arbitrary JSON",
                qc.display()
            );
        }
        // The QC report parses and has the fields a model would read. The
        // refusal is about the MODEL, not the input, so say so precisely.
        bail!(
            "--auto-params: no model is trained yet, so there is no \
             recommendation to print.\n  \
             The shipped default is a good CENTRAL configuration, not a \
             validated universal one: on the 101 real ENA datasets the \
             per-dataset winner beats it on 100 of 101, though the median \
             gap is only 0.0049 and 66 of 101 sit inside a 0.01 margin. \
             There is nothing to switch on until a Phase D model trained \
             on a larger corpus exists and beats the default on held-out \
             studies.\n  \
             See docs/AUTO-PARAMS.md for what the flag would do and what \
             Phase D would change. Until then: run with the default, or \
             pass the ten knobs explicitly."
        );
    }

    pub fn finalize(self) -> Result<(Mode, Runtime)> {
        // The recommendation refusal runs first: --auto-params is a query,
        // not a run, and it must refuse before any output-path or adapter
        // logic can fail for an unrelated reason.
        self.validate_auto_params()?;

        // --auto-params refused above when no model exists, so reaching here
        // without it means a real run: in1/out1 are then mandatory, and the
        // error names them the way the rest of the suite does.
        let (in1, out1) = match (self.in1.clone(), self.out1.clone()) {
            (Some(i), Some(o)) => (i, o),
            (None, _) => bail!("missing required argument -1/--in1 (a run needs an input)"),
            (_, None) => bail!("missing required argument -o/--out1 (a run needs an output)"),
        };

        // Paired-end consistency.
        let paired = self.in2.is_some();
        if paired && self.out2.is_none() {
            bail!("paired-end input (-2/--in2) requires matching -p/--out2");
        }
        if !paired && self.out2.is_some() {
            bail!("-p/--out2 requires matching -2/--in2");
        }
        // The paired-only flags are refused in single-end mode. `--pe-overlap`
        // counts as asked-for when it is anything other than the default
        // `auto`: `auto` resolves to off on single-end anyway (there is no
        // mate), and refusing a bare invocation would be a spurious error. The
        // bare spelling maps to `always`, so it still refuses as it always did.
        if !paired && (self.singletons1.is_some() || self.singletons2.is_some()
            || self.r1_scan_reverse || self.r2_scan_reverse
            || self.pe_ambiguous_smallest
            || self.pe_overlap != PeOverlapMode::Auto)
        {
            bail!("--singletons*, --r1-scan-reverse, --r2-scan-reverse, --pe-overlap and --pe-ambiguous-smallest apply to paired-end mode only");
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
            // 0 disables the rule (see the field docs): the scan then owns
            // every leading run, which is the only way to avoid trimming a
            // leading G that is real biology.
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

        // --force: refuse to clobber any existing output. A dry run writes
        // nothing, so it never has to refuse (CLI-SPEC §1). The first existing
        // path is named — `if let` rather than a `for` because the body always
        // diverges (clippy::never_loop, deny-by-default).
        if !self.force && !self.dry_run {
            if let Some(p) = self
                .out_paths()
                .iter()
                .find(|p| p.as_path().exists())
            {
                bail!(
                    "{} already exists; pass --force (-f) to overwrite it",
                    p.display()
                );
            }
        }

        // Build adapter source. `explicit` records *which branch* produced the
        // set: only the two the user named (preset, FASTA) are explicit, and
        // `--pe-overlap auto` resolves against that flag alone.
        let adapter_source = if self.no_adapter {
            AdapterSource::None
        } else if let Some(preset) = self.adapter_preset {
            let adapters: Vec<Adapter> = match preset {
                AdapterPreset::Illumina => adapter::preset_illumina(),
                AdapterPreset::Nextera  => adapter::preset_nextera(),
                AdapterPreset::SmallRna => adapter::preset_small_rna(),
            };
            AdapterSource::Fixed {
                set: Arc::new(AdapterSet::new(
                    adapters,
                    self.adapter_error_rate,
                    self.min_adapter_overlap,
                )),
                explicit: true,
            }
        } else if let Some(path) = self.adapters {
            let adapters = adapter::load_fasta(&path)?;
            AdapterSource::Fixed {
                set: Arc::new(AdapterSet::new(
                    adapters,
                    self.adapter_error_rate,
                    self.min_adapter_overlap,
                )),
                explicit: true,
            }
        } else {
            // Auto-detect. The terminal k-mer is a fingerprint against the
            // built-in catalogue; the catalogue supplies the full sequence
            // that actually gets trimmed. The k-mer is never used as an
            // adapter itself: it is 12 bases, which both over-trims (a 12-mer
            // within 2 mismatches hits ~1.6% of clean reads by chance) and
            // under-trims (it is the END of the adapter fragment present, so
            // trimming from its position leaves the earlier bases behind).
            const PEEK: usize = 10_000;
            let mut seqs = peek_sequences(&in1, PEEK)?;
            if let Some(p2) = self.in2.as_ref() {
                // Both mates: TruSeq R1 and R2 diverge after a 13-base shared
                // prefix, so a fingerprint from R1's divergent tail would
                // leave R2 untrimmed.
                seqs.extend(peek_sequences(p2, PEEK)?);
            }
            let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
            let det = adapter::identify_from_reads(&refs);

            if det.adapters.is_empty() {
                // D-6 (docs/DECISIONS.md): the safe action is to trim nothing,
                // which is what we do — loudly, with the reason in the report.
                // --strict restores the refusal for anyone who prefers it.
                if self.strict {
                    bail!(
                        "--strict: adapter auto-detection found nothing in {} sampled reads.\n\
                         No sequence was guessed, because trimming on a guess is \
                         worse than not trimming. Pass one of:\n  \
                         --adapter-preset <illumina|nextera|small-rna>\n  \
                         --adapters <FASTA>   (sequence IDs become adapter names)\n  \
                         --no-adapter         (disable adapter trimming)",
                        det.sampled
                    );
                }
                eprintln!(
                    "fastprep: assumption — adapter auto-detection found nothing in {} sampled \
                     reads, so NO adapter trimming will be applied (the safe action; pass \
                     --adapter-preset/--adapters/--no-adapter to be explicit, or --strict to \
                     make this fatal)",
                    det.sampled
                );
                AdapterSource::None
            } else {
                // A run commits every read to this one decision, so say what it
                // was and on what evidence. Informational, so -q gates it.
                if !self.quiet {
                    eprintln!(
                        "fastprep: adapter auto-detection ({} reads sampled)",
                        det.sampled
                    );
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
                }

                AdapterSource::Fixed {
                    set: Arc::new(AdapterSet::new(
                        det.adapters,
                        self.adapter_error_rate,
                        self.min_adapter_overlap,
                    )),
                    // Detection is not configuration: `auto` stays off.
                    explicit: false,
                }
            }
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
            scan_interior: self.interior || !self.no_interior,
            // default ON: the anchoring cut is the measured winner, and the
            // positive flag survives so bench invocations written before 0.4.0
            // keep parsing. The negation wins if both are somehow given.
            boundary_anchored: match (self.boundary_anchored,
                                      self.no_boundary_anchored) {
                (_, true) => false,
                (Some(v), false) => v,
                (None, false) => true,
            },
            overtrim_margin: self.overtrim_margin,
            direction: Direction::Forward,
            mean_q: if self.no_mean_q { None } else { self.mean_q },
            trim_5: self.trim_5,
            trim_3: self.trim_3,
            adapter: match &adapter_source {
                AdapterSource::Fixed { set, .. } => Some(Arc::clone(set)),
                _ => None,
            },
            adapter_g_run_retry: self.adapter_g_run_retry,
            adapter_scored_retry: self.adapter_scored_retry,
            short_tail: self.adapter_short_tail.map(|min| {
                let max_len = min.max(self.min_adapter_overlap.saturating_sub(1));
                // The dataset-level corroboration signal (UT-3). It needs a
                // sample of the reads, so it is computed here rather than
                // inside the per-read path.
                let enrichment = match &adapter_source {
                    AdapterSource::Fixed { set, .. } => {
                        match peek_for_enrichment(&in1, self.in2.as_deref(), SAMPLE) {
                            Ok(sample) => crate::adapter::prefix_enrichment(
                                &sample.iter().map(|s| s.as_slice()).collect::<Vec<_>>(),
                                set,
                                min,
                                max_len,
                                self.adapter_min_enrichment as f64,
                            ),
                            Err(_) => Default::default(),
                        }
                    }
                    _ => Default::default(),
                };
                crate::adapter::ShortTail {
                    min_len: min,
                    // Never reach the ordinary floor: above it the ordinary
                    // matcher already had its chance.
                    max_len,
                    q_drop: self.adapter_q_drop,
                    n_before: 4,
                    n_ref: 10,
                    enrichment,
                    min_factor: self.adapter_min_enrichment as f64,
                    longest_wins: self.short_tail_longest,
                }
            }),
            // Set per pair in `pair_chunks` when `--pe-overlap` fires; `None`
            // here keeps the single-end path and every fallback byte-identical.
            pe_adapter: None,
            extra_committed: 0,
            original_len: None,
            gate_lut: [0.0; 256],
            // The configured offset; `run_single`/`run_paired` re-set it from
            // auto-detection before the first chunk and rebuild the LUT.
            phred_offset,
        };
        base.prepare();
        // Deduplication. Resolved here so the run refuses before allocating
        // rather than after filling the table.
        let dedup = if !self.dedup || self.no_dedup {
            None
        } else {
            let read_len = peek_sequences(&in1, 1)?
                .first()
                .map_or(150, |s| s.len());
            let est = estimate_records(&in1, read_len)?;
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

            let anchor = crate::dedup::AnchorSpec::new(&base, read_len, base.adapter.is_some());
            if !self.quiet {
                if anchor.degenerate() {
                    eprintln!(
                        "fastprep: deduplication on, {}-bit fingerprints, table pre-sized for \
                         {} reads ({} buckets, {:.0} MB); both read ends are trimmed, so the \
                         fingerprint uses an interior window ({} of {} bases)",
                        bits,
                        capacity,
                        proj.buckets,
                        crate::memory::mib(proj.bytes),
                        anchor.len,
                        read_len,
                    );
                } else {
                    eprintln!(
                        "fastprep: deduplication on, {}-bit fingerprints, table pre-sized for \
                         {} reads ({} buckets, {:.0} MB); fingerprint window {} bases at the \
                         {} end",
                        bits,
                        capacity,
                        proj.buckets,
                        crate::memory::mib(proj.bytes),
                        anchor.len,
                        anchor.end.label(),
                    );
                }
            }

            Some(DedupConfig { bits, capacity, read_len, ceiling: cap_bytes })
        };

        let mode = if paired {
            let mut params_r1 = base.clone();
            params_r1.direction = if self.r1_scan_reverse { Direction::Reverse } else { Direction::Forward };
            let mut params_r2 = base;
            params_r2.direction = if self.r2_scan_reverse { Direction::Reverse } else { Direction::Forward };
            Mode::Paired {
                in1: in1.clone(),
                in2: self.in2.unwrap(),
                out1: out1.clone(),
                out2: self.out2.unwrap(),
                singletons1: self.singletons1,
                singletons2: self.singletons2,
                params_r1,
                params_r2,
            }
        } else {
            let mut params = base;
            params.direction = if self.scan_reverse { Direction::Reverse } else { Direction::Forward };
            Mode::Single { in1, out1, params }
        };

        // Resolve the PE-overlap tri-state against the adapter source, in the
        // same place the source is built and through the one function that
        // knows the rule: the guard needs a set to check the boundary against
        // AND the user must have named it, so `auto` runs only for an explicit
        // source. `always` keeps its explicit opt-in, and the detector still
        // refuses a boundary with no set (`detect_pe_overlap_trim(.., None)` is
        // `None`).
        let pe_overlap = if self.pe_overlap.runs_with(&adapter_source) {
            self.pe_overlap
        } else {
            PeOverlapMode::Never
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
            // An explicit --phred-offset is a decision, and auto-detection must
            // not silently overrule it (`conflicts_with` only fires for an
            // explicit --auto-phred, since --auto-phred ships on by default).
            phred_auto: self.auto_phred && !self.no_auto_phred && self.phred_offset.is_none(),
            metadata: self.metadata,
            json: self.json,
            dry_run: self.dry_run,
            quiet: self.quiet,
            allow_total_drop: self.allow_total_drop,
            adapter_source,
            pe_overlap,
            adapter_evidence: self.adapter_evidence,
            pe_ambiguous_smallest: self.pe_ambiguous_smallest,
            dedup,
            adapter_error_rate: self.adapter_error_rate,
            min_adapter_overlap_se: self.min_adapter_overlap,
            thread_pool: None,
        };

        Ok((mode, runtime))
    }
}

/// `--overtrim-margin` accepts 0, 1 and 2 only.
///
/// A hard cap rather than a range: the extension fires on reads whose end cut
/// already overshot a run, so every base of margin spends biology, and the
/// useful range was measured at 0-2. Anything larger is a different feature
/// (a run-base overtrim policy), not a bigger knob on this one.
fn parse_overtrim_margin(s: &str) -> Result<usize, String> {
    let n: usize = s
        .parse()
        .map_err(|_| format!("--overtrim-margin: {s:?} is not a non-negative integer"))?;
    if n > 2 {
        return Err(format!("--overtrim-margin must be 0, 1 or 2 (got {n})"));
    }
    Ok(n)
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
    use crate::algorithm::MinLen;
    use super::*;

    fn base_cli(extra: &[&str]) -> Cli {
        let mut argv = vec!["fastprep", "-1", "in.fq", "-o", "out.fq"];
        argv.extend_from_slice(extra);
        Cli::parse_from(argv)
    }

    #[test]
    fn auto_params_refuses_without_a_model() {
        // A valid fastcheck QC report: the refusal is about the missing
        // MODEL, not the input, so a good report must still refuse.
        let dir = std::env::temp_dir().join("fastprep_autoparams_test");
        std::fs::create_dir_all(&dir).unwrap();
        let qc = dir.join("qc.json");
        std::fs::write(&qc, r#"{"tool":"fastcheck","metrics":{"mean_len":150,"gc_fraction":0.5}}"#)
            .unwrap();
        let cli = base_cli(&["--auto-params", qc.to_str().unwrap(), "--print-args"]);
        let err = cli.validate_auto_params().unwrap_err().to_string();
        assert!(err.contains("no model is trained yet"),
                "the refusal must name the missing model, got: {err}");
        assert!(err.contains("docs/AUTO-PARAMS.md"),
                "the refusal must point at the doc, got: {err}");
    }

    #[test]
    fn auto_params_refuses_a_bad_qc_path() {
        let cli = base_cli(&["--auto-params", "/nonexistent/qc.json", "--print-args"]);
        assert!(cli.validate_auto_params().is_err(),
                "a missing QC report must refuse, not silently recommend");
    }

    #[test]
    fn auto_params_refuses_non_fastcheck_json() {
        let dir = std::env::temp_dir().join("fastprep_autoparams_test");
        std::fs::create_dir_all(&dir).unwrap();
        let qc = dir.join("notqc.json");
        std::fs::write(&qc, r#"{"hello":"world"}"#).unwrap();
        let cli = base_cli(&["--auto-params", qc.to_str().unwrap(), "--print-args"]);
        let err = cli.validate_auto_params().unwrap_err().to_string();
        assert!(err.contains("metrics"),
                "must name the missing metrics block, got: {err}");
    }

    #[test]
    fn auto_params_requires_print_args() {
        let cli = base_cli(&["--auto-params", "qc.json"]);
        assert!(cli.validate_auto_params().is_err(),
                "--auto-params without --print-args must refuse");
    }

    #[test]
    fn print_args_requires_auto_params_at_clap_level() {
        // clap's `requires = "auto_params"` rejects --print-args alone before
        // validate_auto_params ever runs. parse_from would exit the process,
        // so assert through try_parse_from, which returns the error instead.
        let r = Cli::try_parse_from(["fastprep", "-1", "in.fq", "-o", "out.fq",
                                     "--print-args"]);
        assert!(r.is_err(), "clap must reject --print-args without --auto-params");
    }

    #[test]
    fn no_auto_params_is_a_noop() {
        // The refusal path must never fire on a normal run.
        let cli = base_cli(&[]);
        assert!(cli.validate_auto_params().is_ok());
    }

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

    #[test]
    fn shipped_defaults_are_the_swept_winner() {
        // The 90-dataset campaign found one configuration best everywhere; it
        // is what a bare invocation must do. Lock it so a future edit cannot
        // silently drift the default away from the validated configuration.
        use clap::Parser;
        let cli = Cli::parse_from(["fastprep", "-1", "in.fq", "-o", "out.fq"]);
        assert_eq!(cli.adapter_error_rate, 0.2, "--adapter-error-rate");
        assert_eq!(cli.min_adapter_overlap, 14, "--min-adapter-overlap");
        assert_eq!(cli.min_run, 12, "--min-run");
        // 4, not 0: the 7:1 cost model of D-16 prices an artefact left in the
        // output at seven lost bases, and at that anchor the rule pays (its
        // break-even is 2.6). The 1:1 total-error reading prefers 0 and is
        // reachable with an explicit flag.
        assert_eq!(cli.head_dropout, Some(4), "--head-dropout");
        assert_eq!(cli.mean_q, Some(25.0), "--mean-q");
        assert!(matches!(cli.min_len, MinLen::Frac(f) if (f - 0.9).abs() < 1e-9), "--min-len");
        assert!(cli.dedup, "--dedup ships on");
        assert!(cli.auto_phred, "--auto-phred ships on");
        assert!(cli.no_interior, "the interior scan ships OFF");
    }

    #[test]
    fn negation_flags_restore_the_old_behaviour() {
        use clap::Parser;
        let cli = Cli::parse_from(["fastprep", "-1", "in.fq", "-o", "out.fq",
                                   "--no-mean-q", "--no-dedup", "--no-auto-phred", "--interior"]);
        // The negations are applied when the params are built (finalize), so
        // assert the parsed flags that carry them.
        assert!(cli.no_mean_q);
        assert!(cli.no_dedup);
        assert!(cli.no_auto_phred);
        assert!(cli.interior);
    }

    #[test]
    fn pe_overlap_defaults_to_auto() {
        // The default must be the tri-state's `auto`, not a bare bool: it is
        // what resolves against the adapter source in finalize().
        use clap::Parser;
        let cli = Cli::parse_from(["fastprep", "-1", "in.fq", "-2", "in2.fq",
                                   "-o", "out.fq", "-p", "out2.fq"]);
        assert_eq!(cli.pe_overlap, PeOverlapMode::Auto, "--pe-overlap defaults to auto");
    }

    #[test]
    fn auto_tracks_the_adapter_source() {
        // The resolution the pipeline relies on: auto + explicit ⇒ runs;
        // auto + auto-detect or --no-adapter ⇒ off. Robustness note: an adapter
        // set EXISTING is not the same thing as one being CONFIGURED. That is
        // the whole point of `explicit`.
        let set = || Arc::new(AdapterSet::new(adapter::preset_illumina(), 0.2, 14));
        let explicit = AdapterSource::Fixed { set: set(), explicit: true };
        let detected = AdapterSource::Fixed { set: set(), explicit: false };
        assert!(PeOverlapMode::Auto.runs_with(&explicit));
        assert!(
            !PeOverlapMode::Auto.runs_with(&detected),
            "a DETECTED set is not an explicit one: auto must stay off"
        );
        assert!(!PeOverlapMode::Auto.runs_with(&AdapterSource::None));
        assert!(!PeOverlapMode::Never.runs_with(&explicit));
        assert!(PeOverlapMode::Always.runs_with(&detected));
        assert!(PeOverlapMode::Always.runs_with(&AdapterSource::None));
    }

    /// The resolution `finalize` performs, end to end: with no adapter flags the
    /// source comes from auto-detection, which must leave `auto` resolved to
    /// `Never`; `--adapter-preset` is explicit and keeps `Auto`; `always` keeps
    /// today's reach regardless of the source.
    #[test]
    fn finalize_resolves_auto_against_explicit_adapters_only() {
        use clap::Parser;
        let dir = std::env::temp_dir().join(format!("fastprep-pe-auto-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let inp = dir.join("in.fq");
        let inp2 = dir.join("in2.fq");
        let out = dir.join("out.fq");
        let out2 = dir.join("out2.fq");
        // Nextera readthrough in every read, so catalogue detection fires and
        // resolves to a Fixed set that no flag asked for.
        let mut body = String::new();
        for i in 0..200 {
            body.push_str(&format!(
                "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGTCTGTCTCTTATACACATCT\n+\n{}\n",
                "I".repeat(51)
            ));
        }
        std::fs::write(&inp, &body).unwrap();
        std::fs::write(&inp2, &body).unwrap();

        let argv = |extra: &[&str]| -> Vec<String> {
            let mut v: Vec<String> = vec![
                "fastprep".into(),
                "-1".into(),
                inp.to_string_lossy().into_owned(),
                "-2".into(),
                inp2.to_string_lossy().into_owned(),
                "-o".into(),
                out.to_string_lossy().into_owned(),
                "-p".into(),
                out2.to_string_lossy().into_owned(),
                "--no-dedup".into(),
                "--no-auto-phred".into(),
                "-q".into(),
            ];
            v.extend(extra.iter().map(|s| s.to_string()));
            v
        };

        // 1. No adapter flags → detected → `auto` is off.
        let (_, rt) = Cli::parse_from(argv(&[])).finalize().unwrap();
        assert!(
            matches!(rt.adapter_source, AdapterSource::Fixed { explicit: false, .. }),
            "auto-detection must be recorded as NOT explicit"
        );
        assert_eq!(
            rt.pe_overlap,
            PeOverlapMode::Never,
            "auto must not run on a detected adapter set"
        );

        // 2. --adapter-preset is explicit → `auto` stays on.
        let (_, rt) = Cli::parse_from(argv(&["--adapter-preset", "illumina"]))
            .finalize()
            .unwrap();
        assert!(matches!(rt.adapter_source, AdapterSource::Fixed { explicit: true, .. }));
        assert_eq!(rt.pe_overlap, PeOverlapMode::Auto);

        // 3. `always` keeps today's reach, detected set or not.
        let (_, rt) = Cli::parse_from(argv(&["--pe-overlap", "always"]))
            .finalize()
            .unwrap();
        assert_eq!(rt.pe_overlap, PeOverlapMode::Always);

        std::fs::remove_file(&inp).ok();
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn the_three_new_trimming_flags_default_off() {
        use clap::Parser;
        let cli = Cli::parse_from(["fastprep", "-1", "in.fq", "-o", "out.fq"]);
        assert_eq!(cli.overtrim_margin, 0, "--overtrim-margin defaults to 0 (off)");
        assert!(!cli.short_tail_longest, "--short-tail-longest defaults off");
        assert!(
            !cli.pe_ambiguous_smallest,
            "--pe-ambiguous-smallest defaults off"
        );
    }

    #[test]
    fn boundary_anchored_defaults_on_and_has_a_negation() {
        use clap::Parser;
        // 0.4.0 flipped the default: the anchoring cut is the measured winner
        // (+0.0121/read at 7:1, plan §7), and the pre-flip cut stays reachable.
        // finalize() reads the input for adapter auto-detection, so the test
        // needs a real (tiny) FASTQ: 94 bp with a qualifying 3' G-run.
        let dir = std::env::temp_dir().join(format!("fp-ba-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let inp = dir.join("in.fq");
        let out = dir.join("out.fq");
        std::fs::write(
            &inp,
            format!("@r1\n{}\n+\n{}\n", "ACGT".repeat(20) + &"G".repeat(14), "I".repeat(94)),
        )
        .unwrap();

        let anchored = |extra: &[&str]| -> bool {
            let mut v: Vec<String> = vec![
                "fastprep".into(),
                "-1".into(), inp.to_string_lossy().into_owned(),
                "-o".into(), out.to_string_lossy().into_owned(),
                "--no-auto-phred".into(), "--no-dedup".into(), "-q".into(),
            ];
            v.extend(extra.iter().map(|s| s.to_string()));
            match Cli::parse_from(v).finalize().unwrap().0 {
                Mode::Single { params, .. } => params.boundary_anchored,
                _ => panic!("expected single-end"),
            }
        };

        assert!(anchored(&[]), "anchoring is ON by default since 0.4.0");
        assert!(
            anchored(&["--boundary-anchored"]),
            "the pre-0.4.0 spelling still parses and is now a no-op"
        );
        assert!(
            !anchored(&["--no-boundary-anchored"]),
            "--no-boundary-anchored restores the mercy-walk overshoot"
        );
        assert!(
            !anchored(&["--boundary-anchored=false"]),
            "--boundary-anchored=false is the second escape"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn overtrim_margin_accepts_zero_to_two_and_refuses_more() {
        use clap::Parser;
        for n in 0..=2usize {
            let cli = Cli::parse_from([
                "fastprep",
                "-1",
                "in.fq",
                "-o",
                "out.fq",
                "--overtrim-margin",
                &n.to_string(),
            ]);
            assert_eq!(cli.overtrim_margin, n);
        }
        for bad in ["3", "9", "-1", "x"] {
            assert!(
                Cli::try_parse_from([
                    "fastprep",
                    "-1",
                    "in.fq",
                    "-o",
                    "out.fq",
                    "--overtrim-margin",
                    bad,
                ])
                .is_err(),
                "--overtrim-margin {bad} must be refused"
            );
        }
    }

    #[test]
    fn short_tail_longest_requires_the_short_tail_path() {
        use clap::Parser;
        assert!(
            Cli::try_parse_from([
                "fastprep",
                "-1",
                "in.fq",
                "-o",
                "out.fq",
                "--short-tail-longest",
            ])
            .is_err(),
            "--short-tail-longest means nothing without --adapter-short-tail"
        );
        let cli = Cli::parse_from([
            "fastprep",
            "-1",
            "in.fq",
            "-o",
            "out.fq",
            "--adapter-short-tail",
            "6",
            "--short-tail-longest",
        ]);
        assert!(cli.short_tail_longest);
    }

    /// A flag that parses but is never wired into `Params`/`ShortTail`/
    /// `Runtime` is a silent no-op: this is the end-to-end check that each of
    /// the three new flags actually reaches the code that acts on it.
    #[test]
    fn finalize_wires_the_new_trimming_flags() {
        use clap::Parser;
        let dir = std::env::temp_dir().join(format!("fastprep-new-flags-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let inp = dir.join("in.fq");
        let inp2 = dir.join("in2.fq");
        let out = dir.join("out.fq");
        let out2 = dir.join("out2.fq");
        let body = "@r0\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        std::fs::write(&inp, body).unwrap();
        std::fs::write(&inp2, body).unwrap();

        let (mode, rt) = Cli::parse_from([
            "fastprep",
            "-1",
            inp.to_str().unwrap(),
            "-2",
            inp2.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-p",
            out2.to_str().unwrap(),
            "--no-dedup",
            "--no-auto-phred",
            "-q",
            "--adapter-preset",
            "illumina",
            "--overtrim-margin",
            "2",
            "--adapter-short-tail",
            "6",
            "--short-tail-longest",
            "--pe-ambiguous-smallest",
        ])
        .finalize()
        .unwrap();

        let params = match mode {
            Mode::Paired { params_r1, .. } => params_r1,
            _ => panic!("paired-end mode expected"),
        };
        assert_eq!(params.overtrim_margin, 2, "--overtrim-margin must reach Params");
        let st = params
            .short_tail
            .as_ref()
            .expect("--adapter-short-tail must build a ShortTail");
        assert!(st.longest_wins, "--short-tail-longest must reach the ShortTail");
        assert!(
            rt.pe_ambiguous_smallest,
            "--pe-ambiguous-smallest must reach the Runtime"
        );

        std::fs::remove_dir_all(&dir).ok();
    }
}

