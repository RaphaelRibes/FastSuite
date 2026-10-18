use std::path::PathBuf;

use clap::Parser;

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
    /// Set to 0 to disable preview generation entirely.
    #[arg(long, default_value = "5", value_name = "N")]
    pub preview_reads: usize,

    /// Approximate total RAM budget for the worker index buffer before it is
    /// spilled to the spiller thread (megabytes). Preview reads are streamed
    /// to disk independently and do not count against it.
    /// Measured on a real 5M-read library: ~141 MB peak RSS at the 512 MB
    /// default, ~36 MB at 8 MB, wall time unchanged. The budget is divided
    /// across workers, so an input file smaller than a worker's share is held
    /// whole without spilling.
    #[arg(long, default_value = "512", value_name = "MB")]
    pub buffer_mb: usize,

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
}
