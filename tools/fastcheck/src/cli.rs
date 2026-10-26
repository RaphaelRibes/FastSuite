use clap::Parser;
use std::path::PathBuf;

use anyhow::Context;

use crate::overrep::OverrepMode;

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
    /// Optional processed (e.g. trimmed) files to compare against the raw ones:
    /// the report then shows a before/after side-by-side. One file for
    /// single-end, two for paired-end.
    #[arg(long, value_name = "PATH")]
    pub processed_1: Option<PathBuf>,
    #[arg(long, value_name = "PATH")]
    pub processed_2: Option<PathBuf>,


    /// Write the machine-readable report (fastsuite.report.v1) to PATH.
    #[arg(long, value_name = "PATH")]
    pub json: Option<PathBuf>,
    /// Write a fastsim profile directly, fitted from this run's measurements:
    /// the file `fastsim --profile FILE` consumes, with no separate fit step.
    #[arg(long, value_name = "PATH")]
    pub profile: Option<PathBuf>,

    /// Write a FastQC-style HTML report (exact counts, inline SVG, no JS).
    #[arg(long, value_name = "PATH")]
    pub html: Option<PathBuf>,

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

    // NOTE: this value is bounded by MAX_READ_LEN, checked by
    // `checked_max_read_len` before anything is allocated. `--insert-hist` sizes two
    // bin vectors from it, which is how an absurd value used to reach the allocator.

    /// Refuse to start if the projected fingerprint table exceeds this many megabytes.
    /// Defaults to the enforced memory limit of the process (cgroup v2/v1, SLURM,
    /// RLIMIT_AS, then MemAvailable — the minimum of them, per tools/fastprep's model).
    #[arg(long, value_name = "MB")]
    pub max_mem: Option<u64>,

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

    /// How over-represented sequences are identified (PLAN-QC §5.9). The one honest
    /// trade-off in fastcheck: exact identity needs either the sequences themselves or a
    /// second read of the file.
    ///
    ///   none        (default, zero cost) exact duplication levels, no sequence names
    ///   fingerprint +unique x 8 B/0.8 RAM; exact multiplicity, names_available: false
    ///   second-pass one extra read of the file to emit the text of sequences above the
    ///               threshold; memory stays at the fingerprint table; passes: 2
    ///   table       2-bit packed sequences (ceil(len/4)+8 B each); refuses above --max-mem;
    ///               unique_sequences_retained: true
    ///
    /// FastQC instead analyses the first 100 000 sequences and hopes. None of these modes
    /// samples: a sample is never presented as the whole.
    #[arg(long, value_name = "none|fingerprint|second-pass|table", default_value = "none")]
    pub overrep_mode: String,

    /// Minimum multiplicity for a sequence to be reported as over-represented
    /// (default: max(10, 0.1% of reads) — FastQC's 0.1% convention, but computed from the
    /// exact read count rather than from a 100 k sample).
    #[arg(long, value_name = "N")]
    pub overrep_min_count: Option<u32>,

    /// Count, per position, the reads whose adapter match STARTS there, for every adapter
    /// in the catalogue (PLAN-QC module 11). Off by default: it costs one catalogue scan
    /// per read. The matcher semantics are fastprep's, exactly.
    #[arg(long)]
    pub adapter_content: bool,

    /// Adapter catalogue for --adapter-content: illumina (TruSeq R1/R2 + Nextera),
    /// nextera, or a FASTA path.
    #[arg(long, value_name = "illumina|nextera|PATH", default_value = "illumina")]
    pub adapter_preset: String,

    /// Maximum mismatch rate for --adapter-content matching [0, 1] (fastprep's default).
    #[arg(long, default_value_t = 0.2, value_name = "FLOAT")]
    pub adapter_error_rate: f64,

    /// Minimum adapter overlap for --adapter-content matching (fastprep's default).
    #[arg(long, default_value_t = 14, value_name = "N")]
    pub min_adapter_overlap: usize,

    /// Per-tile quality (PLAN-QC module 7 / FastQC's per-tile modules): quality by
    /// tile x position x quality bin, parsed from 7-field Illumina headers
    /// (`@instrument:run:flowcell:lane:tile:x:y`). OFF by default: it costs ~5 MB per 100
    /// tiles at R=150, and only raw Illumina runs carry the coordinates. Reads whose
    /// header has no tile field are counted and reported, never silently dropped.
    #[arg(long)]
    pub tiles: bool,

    /// Run only these modules (comma-separated): bases, quality, length, gc, heatmap,
    /// dup, polyx, kmer, adapter, tile. The parser then materialises only what the
    /// enabled modules need, so `--only dup` touches neither the quality bytes nor the
    /// heatmap grid, and skips the per-read loops entirely. Absent = every module on,
    /// which is byte for byte the pre--only behaviour; `adapter` and `tile` stay opt-in
    /// behind --adapter-content / --tiles either way.
    #[arg(long, value_name = "MODULE[,MODULE...]")]
    pub only: Option<String>,

    /// Disable the per-base quality accumulators. The quality line is still read for
    /// framing, but no quality byte is decoded: every quality-derived metric is absent
    /// from the report and the Phred offset becomes a stated default, not a measurement.
    #[arg(long)]
    pub no_quality: bool,

    /// Never build the composition grid (the 256x256 heatmap).
    #[arg(long)]
    pub no_heatmap: bool,

    /// Suppress the end-of-run summary line.
    #[arg(short = 'q', long)]
    pub quiet: bool,
}

/// Longest read `--max-read-len` accepts.
///
/// The flag sizes the per-position planes, which cost `128 x 8 B` = 1 KiB per position
/// per file, so 1 Mi positions already project 1 GiB of accumulator state; on top of
/// that `--insert-hist` allocates two bin vectors of `2 * max_read_len + 2` u64. The
/// value is therefore bounded once, here, before any of that is requested: a typo such
/// as 1099511627776 used to be taken at face value and died inside the insert
/// allocator ("memory allocation of 17592186044424 bytes failed") after the whole scan.
pub const MAX_READ_LEN: usize = 1 << 20;

impl Cli {
    /// `--max-read-len`, bounded to what the per-position planes can carry.
    pub fn checked_max_read_len(&self) -> anyhow::Result<usize> {
        if self.max_read_len > MAX_READ_LEN {
            anyhow::bail!(
                "--max-read-len {} is above the supported maximum {MAX_READ_LEN}: the \
                 per-position planes cost 1 KiB per position per file and --insert-hist \
                 allocates two vectors of 2 x max-read-len + 2 bins, so a larger value is a \
                 multi-terabyte request, not a longer read length",
                self.max_read_len
            );
        }
        Ok(self.max_read_len)
    }

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
            return self.threads;
        }
        // `FASTSUITE_THREADS` supplies the default when `-j` is not given — the suite's
        // convention (fastprep, fastmaps, fastsim all read it): an explicit `-j` wins
        // outright, and unset / unparseable / 0 all mean "not set".
        std::env::var("FASTSUITE_THREADS")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            })
    }

    pub fn overrep_mode(&self) -> anyhow::Result<OverrepMode> {
        OverrepMode::parse(&self.overrep_mode)
    }

    /// The adapter set for --adapter-content, with fastprep's presets.
    pub fn adapter_set(&self) -> anyhow::Result<crate::adapter::AdapterSet> {
        if !(0.0..=1.0).contains(&self.adapter_error_rate) {
            anyhow::bail!("--adapter-error-rate must be in [0, 1]");
        }
        if self.min_adapter_overlap == 0 {
            anyhow::bail!("--min-adapter-overlap must be >= 1");
        }
        let adapters = match self.adapter_preset.as_str() {
            "illumina" => crate::adapter::preset_illumina(),
            "nextera" => crate::adapter::preset_nextera(),
            path => {
                let text = std::fs::read_to_string(path)
                    .with_context(|| format!("reading adapter FASTA {}", path))?;
                parse_fasta_adapters(&text)?
            }
        };
        Ok(crate::adapter::AdapterSet::new(
            adapters,
            self.adapter_error_rate,
            self.min_adapter_overlap,
        ))
    }

    /// True when the run was asked to measure a subset of the modules — i.e. when a
    /// gating flag was given. Absent everywhere, the ungated path (and the unchanged
    /// report) is taken; the report only carries a `modules` array when this is true.
    pub fn gating_requested(&self) -> bool {
        self.only.is_some() || self.no_quality || self.no_heatmap
    }

    /// Which modules this run measures.
    ///
    /// `--only` is an allow-list; `--no-quality` and `--no-heatmap` are deny-list entries
    /// that also apply on top of it. An unknown name, an empty selection, or a selection
    /// whose every entry is disabled is refused: running no modules at all and writing an
    /// empty report is the one outcome this must never produce silently.
    pub fn modules(&self) -> anyhow::Result<Modules> {
        let mut m = match &self.only {
            None => {
                // The tile module rides on its own opt-in flag, like adapter content.
                let mut m = Modules::defaults(self.adapter_content);
                m.tile = self.tiles;
                m
            }
            Some(list) => {
                let names: Vec<&str> =
                    list.split(',').map(str::trim).filter(|s| !s.is_empty()).collect();
                if names.is_empty() {
                    anyhow::bail!(
                        "--only: no module named (known: {})",
                        MODULE_NAMES.join(", ")
                    );
                }
                let mut m = Modules::none();
                for name in names {
                    match name {
                        "bases" => m.bases = true,
                        "quality" => m.quality = true,
                        "length" => m.length = true,
                        "gc" => m.gc = true,
                        "heatmap" => m.heatmap = true,
                        "dup" => m.dup = true,
                        "polyx" => m.polyx = true,
                        "kmer" => m.kmer = true,
                        // `adapter` is the --adapter-content module: naming it here does not
                        // switch that module on, it only says it may run.
                        "adapter" => m.adapter = self.adapter_content,
                        // `tile` is the --tiles module, opt-in the same way.
                        "tile" => m.tile = self.tiles,
                        other => anyhow::bail!(
                            "--only: unknown module {other:?} (known: {})",
                            MODULE_NAMES.join(", ")
                        ),
                    }
                }
                m
            }
        };
        if self.no_quality {
            m.quality = false;
        }
        if self.no_heatmap {
            m.heatmap = false;
        }
        if self.adapter_content && self.only.is_some() && !m.adapter {
            anyhow::bail!(
                "--adapter-content is not selected by --only: add `adapter` to the list, \
                 or drop --adapter-content"
            );
        }
        if self.tiles && self.only.is_some() && !m.tile {
            anyhow::bail!(
                "--tiles is not selected by --only: add `tile` to the list, or drop --tiles"
            );
        }
        if !m.any() {
            anyhow::bail!(
                "--only: nothing left to run — the selection was empty, or every module in it \
                 was disabled by --no-quality/--no-heatmap, or `adapter`/`tile` was named \
                 without --adapter-content/--tiles (known: {})",
                MODULE_NAMES.join(", ")
            );
        }
        Ok(m)
    }
}

/// Every module name, in the order they are reported. An unknown name in `--only` is
/// refused against this list rather than ignored.
pub const MODULE_NAMES: &[&str] = &[
    "bases", "quality", "length", "gc", "heatmap", "dup", "polyx", "kmer", "adapter", "tile",
];

/// Which modules a run measures.
///
/// `bases` owns the per-position composition, `quality` the per-position quality and the
/// byte census, `length` the read-length histogram, `gc` the ACGT/GC histograms, `dup`
/// the exact fingerprint table, `polyx` the homopolymer runs, `kmer` the 12-mer counter,
/// `heatmap` the composition grid, `adapter` the per-position adapter content, `tile` the
/// per-tile quality. Modules left off are not measured *and* not reported: an absent field
/// never reads as a zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Modules {
    pub bases: bool,
    pub quality: bool,
    pub length: bool,
    pub gc: bool,
    pub heatmap: bool,
    pub dup: bool,
    pub polyx: bool,
    pub kmer: bool,
    pub adapter: bool,
    pub tile: bool,
}

impl Modules {
    /// No module at all — the starting point of a `--only` selection.
    pub fn none() -> Self {
        Self {
            bases: false,
            quality: false,
            length: false,
            gc: false,
            heatmap: false,
            dup: false,
            polyx: false,
            kmer: false,
            adapter: false,
            tile: false,
        }
    }

    /// Today's default: every module except the opt-in tile and adapter content, which
    /// stay behind their own flags (`modules()` adds them back when asked for).
    pub fn defaults(adapter_content: bool) -> Self {
        Self {
            bases: true,
            quality: true,
            length: true,
            gc: true,
            heatmap: true,
            dup: true,
            polyx: true,
            kmer: true,
            adapter: adapter_content,
            tile: false,
        }
    }

    /// At least one module runs. A run that measures nothing is refused up front.
    pub fn any(&self) -> bool {
        !self.enabled_names().is_empty()
    }

    /// The enabled names, in `MODULE_NAMES` order.
    pub fn enabled_names(&self) -> Vec<&'static str> {
        MODULE_NAMES
            .iter()
            .copied()
            .filter(|n| match *n {
                "bases" => self.bases,
                "quality" => self.quality,
                "length" => self.length,
                "gc" => self.gc,
                "heatmap" => self.heatmap,
                "dup" => self.dup,
                "polyx" => self.polyx,
                "kmer" => self.kmer,
                "adapter" => self.adapter,
                "tile" => self.tile,
                _ => false,
            })
            .collect()
    }
}

/// Parse an adapter FASTA: `>name` then the sequence on the following lines.
fn parse_fasta_adapters(text: &str) -> anyhow::Result<Vec<crate::adapter::Adapter>> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.push(crate::adapter::Adapter { name: n, seq: std::mem::take(&mut seq).into_bytes() });
            }
            name = Some(rest.split_whitespace().next().unwrap_or("adapter").to_string());
        } else if !line.is_empty() {
            seq.push_str(&line.to_ascii_uppercase());
        }
    }
    if let Some(n) = name {
        out.push(crate::adapter::Adapter { name: n, seq: seq.into_bytes() });
    }
    if out.is_empty() {
        anyhow::bail!("adapter FASTA contained no records");
    }
    Ok(out)
}
