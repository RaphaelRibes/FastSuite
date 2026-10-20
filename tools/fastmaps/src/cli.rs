use std::fmt;
use std::path::PathBuf;

use clap::Parser;

/// What to record in the per-read index stream that `heatmap.data` is built
/// from (FM-1).
///
/// The stream is the only structural cost in the tool: 12 B per read in RAM,
/// 8 B per read of temp disk, and a full render pass — all of it to serve one
/// output file. The counting grid, the PNG, the corner metrics and the report
/// come from `state.counts`, which needs no per-read record, so every mode
/// below is exact for them.
///
/// * `Exact`  — every read, as before. The default.
/// * `Capped(n)` — at most `n` reads per `(pixel, file)`: `heatmap.data` becomes
///   a labelled sample (pixel *membership* is sampled; pixel *counts* are not).
/// * `Off` — no per-read index buffer and no index spill, so there is nothing
///   for the render pass to read: `heatmap.data` is written as empty rows, so
///   the file and its offsets stay valid. The preview stream has its own path
///   and is unaffected by this mode: `heatmap.preview` is still written.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum DataMode {
    #[default]
    Exact,
    Capped(u8),
    Off,
}

impl DataMode {
    pub fn is_exact(self) -> bool {
        matches!(self, DataMode::Exact)
    }

    /// Reads stored per `(pixel, file)`; `None` in exact mode (unbounded).
    pub fn cap(self) -> Option<u8> {
        match self {
            DataMode::Exact => None,
            DataMode::Capped(n) => Some(n),
            DataMode::Off => Some(0),
        }
    }

    /// Value for the `data_reads_per_pixel` JSON/report field in a non-exact
    /// mode (0 in exact mode, where the field is omitted).
    pub fn reads_per_pixel(self) -> u8 {
        self.cap().unwrap_or(0)
    }
}

impl fmt::Display for DataMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DataMode::Exact => f.write_str("exact"),
            DataMode::Capped(n) => write!(f, "capped:{n}"),
            DataMode::Off => f.write_str("off"),
        }
    }
}

/// Parse `--data-mode`: `exact`, `off`, or `capped:N` (`capped=N` is accepted
/// too, as is a space-free `capped` with the count attached by the caller).
pub fn parse_data_mode(s: &str) -> Result<DataMode, String> {
    let t = s.trim();
    if t.eq_ignore_ascii_case("exact") {
        return Ok(DataMode::Exact);
    }
    if t.eq_ignore_ascii_case("off") {
        return Ok(DataMode::Off);
    }
    let lower = t.to_ascii_lowercase();
    if let Some(rest) = lower.strip_prefix("capped") {
        let rest = rest.trim_start_matches([':', '=']);
        let n: u16 = rest
            .parse()
            .map_err(|_| format!("--data-mode: {s:?} is not exact|off|capped:N"))?;
        if n == 0 {
            return Err("--data-mode capped:0 is meaningless; use --data-mode off".into());
        }
        if n > u8::MAX as u16 {
            return Err(format!("--data-mode capped max is {}", u8::MAX));
        }
        return Ok(DataMode::Capped(n as u8));
    }
    Err(format!(
        "--data-mode: {s:?} is not exact|off|capped:N (e.g. capped:4)"
    ))
}

/// 2-D heatmap of FASTQ reads by base composition.
///
/// Each read is encoded as (x = T+C count, y = T+G count) and counted
/// into a grid. Outputs heatmap.json / .data / .preview (for the marimo
/// notebook), heatmap.png, and heatmap_3d.html.
#[derive(Parser, Debug)]
#[command(name = "fastmaps", version)]
pub struct Args {
    /// Input FASTQ files (plain or .gz).
    #[arg(required = true, value_name = "FILE")]
    pub files: Vec<String>,

    /// Extract reads whose (x,y) falls inside a bounding box.
    /// Repeatable: each box produces its own output file.
    ///
    /// Format: x1,y1:x2,y2
    /// Example: --corners 0,140:5,150
    /// Output:  <output-dir>/corners_0_140_5_150.fq
    #[arg(long, value_name = "x1,y1:x2,y2", action = clap::ArgAction::Append)]
    pub corners: Vec<String>,

    /// Number of preview reads stored per pixel in heatmap.preview.
    /// Set to 0 to store no reads: heatmap.preview is still written, as empty
    /// rows, so its row offsets stay valid.
    #[arg(long, default_value = "5", value_name = "N")]
    pub preview_reads: usize,

    /// Approximate total RAM budget for the worker index buffer before it is
    /// spilled to the spiller thread (megabytes). The budget is divided across
    /// workers, so an input file smaller than a worker's share is held whole
    /// without spilling.
    ///
    /// Preview reads are streamed to disk on their own threshold too, but that
    /// threshold is derived from this one (bounded by 256 KiB) instead of being
    /// a fixed 4 MiB per worker — before that, previews imposed a memory floor
    /// `--buffer-mb 1` could not move.
    ///
    /// Measured: on a real 5M-read library ~141 MB peak RSS at 512 MB against
    /// ~36 MB at 8 MB, and on 2M real reads 48.6 MB against 22.7 MB, at
    /// unchanged wall time and byte-identical outputs — the default is 8 MB for
    /// that reason. Raising it buys nothing unless the input is large enough for
    /// a worker's share to bind.
    #[arg(long, default_value = "8", value_name = "MB")]
    pub buffer_mb: usize,

    /// Per-read index stream mode for heatmap.data: `exact` (default) stores
    /// every read, `capped:N` stores at most N reads per (pixel, file) and
    /// `off` stores none (no index buffer, no index spill, nothing to render).
    /// `heatmap.preview` is written in every mode.
    ///
    /// Only heatmap.data is affected: heatmap.json carries
    /// `"data_reads_per_pixel"` and `"data_complete": false` when the mode is
    /// not exact, and the counts, PNG, corner metrics and report stay exact.
    #[arg(
        long,
        default_value = "exact",
        value_name = "exact|capped:N|off",
        value_parser = parse_data_mode
    )]
    pub data_mode: DataMode,

    /// Write all output files into this directory (created if absent).
    #[arg(short = 'o', long, default_value = ".", value_name = "PATH")]
    pub output_dir: PathBuf,

    /// Directory for temporary per-row spill files.
    /// Defaults to ./tmp_<pid>. Cleaned up automatically on success;
    /// if the job is killed, delete the directory manually before re-running.
    #[arg(long, value_name = "PATH")]
    pub temp_dir: Option<PathBuf>,

    /// Also export a vector SVG heatmap alongside the PNG and HTML outputs.
    #[arg(long)]
    pub svg: bool,

    /// Glob pattern(s) marking a file as R1 (case-insensitive). Repeatable.
    /// Supports `*` (any sequence) and `?` (one character) wildcards.
    /// If any --r1-pattern is given, ALL defaults are overridden — pass each
    /// pattern you want explicitly.
    /// Default: *R1*
    #[arg(long, value_name = "GLOB", action = clap::ArgAction::Append)]
    pub r1_pattern: Vec<String>,

    /// Glob pattern(s) marking a file as R2 (case-insensitive). Repeatable.
    /// Same override semantics as --r1-pattern.
    /// Default: *R2*
    #[arg(long, value_name = "GLOB", action = clap::ArgAction::Append)]
    pub r2_pattern: Vec<String>,

    /// Write the machine-readable FastSuite report (fastsuite.report.v1) to
    /// PATH, alongside the grid JSON. The grid JSON keeps its default
    /// location `<output-dir>/heatmap.json`; this flag only redirects it
    /// and additionally writes the v1 wrapper. `-` is not accepted:
    /// fastmaps is file-only (D-4).
    #[arg(long, value_name = "PATH")]
    pub json: Option<std::path::PathBuf>,

    /// Suppress the end-of-run summary lines on stdout.
    /// Progress and warnings still go to stderr.
    #[arg(long, short = 'q')]
    pub quiet: bool,

    /// Refuse to overwrite an output directory that already contains
    /// fastmaps outputs; without this flag the run exits 1 naming the file.
    #[arg(long, short = 'f')]
    pub force: bool,

    /// Number of worker threads to use. Defaults to the number of available
    /// CPU cores (bounded by the file count), or $FASTSUITE_THREADS when set.
    /// Pass 0 for auto.
    #[arg(short = 'j', long, default_value_t = 0)]
    pub threads: usize,
}

impl Args {
    /// Effective worker count, mirroring fastsim's `FASTSUITE_THREADS`
    /// contract (CLI-SPEC §1, gap #6) and fastdedup's `resolve_thread_count`.
    ///
    /// `-j N` (N > 0) wins outright; otherwise `$FASTSUITE_THREADS` is read
    /// when set to a parseable non-zero value; otherwise all available cores.
    /// `max` (the input file count) always clamps the result — there is no
    /// point in workers with no files, exactly as before this env default.
    pub fn resolved_threads(&self, max: usize) -> usize {
        let base = if self.threads > 0 {
            self.threads
        } else {
            std::env::var("FASTSUITE_THREADS")
                .ok()
                .and_then(|v| v.trim().parse::<usize>().ok())
                .filter(|&n| n > 0)
                .unwrap_or_else(|| {
                    std::thread::available_parallelism()
                        .map(|n| n.get())
                        .unwrap_or(4)
                })
        };
        base.min(max.max(1))
    }
}

/// Parse "x1,y1:x2,y2" into normalised corners `((min_x,min_y),(max_x,max_y))`.
pub fn parse_corners(s: &str) -> Option<((u32, u32), (u32, u32))> {
    let (a, b) = s.split_once(':')?;
    let (x1, y1) = a.split_once(',')?;
    let (x2, y2) = b.split_once(',')?;
    let x1: u32 = x1.trim().parse().ok()?;
    let y1: u32 = y1.trim().parse().ok()?;
    let x2: u32 = x2.trim().parse().ok()?;
    let y2: u32 = y2.trim().parse().ok()?;
    Some(((x1.min(x2), y1.min(y2)), (x1.max(x2), y1.max(y2))))
}

/// Build the output filename for a corners box.
/// e.g. ((0,140),(5,150)) → "corners_0_140_5_150.fq"
pub fn corners_filename(((x1, y1), (x2, y2)): ((u32, u32), (u32, u32))) -> String {
    format!("corners_{x1}_{y1}_{x2}_{y2}.fq")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_corners_basic() {
        assert_eq!(parse_corners("0,2:4,6"), Some(((0, 2), (4, 6))));
    }

    #[test]
    fn test_parse_corners_normalizes() {
        assert_eq!(parse_corners("4,6:0,2"), Some(((0, 2), (4, 6))));
    }

    #[test]
    fn test_parse_corners_invalid() {
        assert!(parse_corners("not_valid").is_none());
        assert!(parse_corners("0,2").is_none());
    }

    #[test]
    fn test_corners_filename() {
        assert_eq!(corners_filename(((0, 140), (5, 150))), "corners_0_140_5_150.fq");
    }

    #[test]
    fn test_json_arg_default_none() {
        let a = Args::parse_from(["fastmaps", "a.fq"]);
        assert_eq!(a.json, None);
    }

    #[test]
    fn test_json_arg_parsed() {
        let a = Args::parse_from(["fastmaps", "a.fq", "--json", "/tmp/r.json"]);
        assert_eq!(a.json.as_deref(), Some("/tmp/r.json".as_ref()));
    }

    #[test]
    fn test_force_default_false() {
        let a = Args::parse_from(["fastmaps", "a.fq"]);
        assert!(!a.force);
    }

    #[test]
    fn test_force_short_and_long() {
        assert!(Args::parse_from(["fastmaps", "a.fq", "--force"]).force);
        assert!(Args::parse_from(["fastmaps", "a.fq", "-f"]).force);
    }

    #[test]
    fn test_quiet_default_false() {
        assert!(!Args::parse_from(["fastmaps", "a.fq"]).quiet);
    }

    #[test]
    fn test_quiet_short_and_long() {
        assert!(Args::parse_from(["fastmaps", "a.fq", "-q"]).quiet);
        assert!(Args::parse_from(["fastmaps", "a.fq", "--quiet"]).quiet);
    }

    #[test]
    fn test_resolved_threads_prefers_flag() {
        let a = Args::parse_from(["fastmaps", "a.fq", "b.fq", "-j", "3"]);
        assert_eq!(a.resolved_threads(4), 3);
    }

    #[test]
    fn test_resolved_threads_flag_bounded_by_files() {
        let a = Args::parse_from(["fastmaps", "a.fq", "b.fq", "-j", "8"]);
        assert_eq!(a.resolved_threads(2), 2); // caller clamps by file count
    }

    #[test]
    fn test_resolved_threads_zero_means_auto() {
        let a = Args::parse_from(["fastmaps", "a.fq", "-j", "0"]);
        // No env in tests: falls back to available cores, clamped by file count.
        let n = a.resolved_threads(1);
        assert!(n >= 1 && n <= 16, "got {n}");
    }

    // ── FM-1: --data-mode ───────────────────────────────────────────────────

    #[test]
    fn test_data_mode_default_is_exact() {
        let a = Args::parse_from(["fastmaps", "a.fq"]);
        assert_eq!(a.data_mode, DataMode::Exact);
        assert!(a.data_mode.is_exact());
        assert_eq!(a.data_mode.cap(), None);
    }

    #[test]
    fn test_data_mode_parses_all_forms() {
        let p = |s: &str| parse_data_mode(s).unwrap();
        assert_eq!(p("exact"), DataMode::Exact);
        assert_eq!(p("EXACT"), DataMode::Exact);
        assert_eq!(p("off"), DataMode::Off);
        assert_eq!(p("capped:1"), DataMode::Capped(1));
        assert_eq!(p("capped:255"), DataMode::Capped(255));
        assert_eq!(p("capped=7"), DataMode::Capped(7));
        assert_eq!(p("capped:4").to_string(), "capped:4");
        assert_eq!(p("off").to_string(), "off");
        assert_eq!(DataMode::Off.reads_per_pixel(), 0);
        assert_eq!(DataMode::Capped(6).reads_per_pixel(), 6);

        let a = Args::parse_from(["fastmaps", "a.fq", "--data-mode", "off"]);
        assert_eq!(a.data_mode, DataMode::Off);
        let a = Args::parse_from(["fastmaps", "a.fq", "--data-mode", "capped:3"]);
        assert_eq!(a.data_mode, DataMode::Capped(3));
    }

    #[test]
    fn test_data_mode_rejects_garbage() {
        assert!(parse_data_mode("capped:0").is_err(), "cap 0 is `off`, not 0 reads");
        assert!(parse_data_mode("capped:256").is_err(), "u8 cap");
        assert!(parse_data_mode("capped:").is_err());
        assert!(parse_data_mode("mostly").is_err());
        assert!(parse_data_mode("").is_err());
        assert!(Args::try_parse_from(["fastmaps", "a.fq", "--data-mode", "nope"]).is_err());
    }
}
