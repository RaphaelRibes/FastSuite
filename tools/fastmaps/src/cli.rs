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

    /// Approximate total RAM budget for worker index/preview buffers
    /// before they are spilled to the spiller thread (megabytes).
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
}
