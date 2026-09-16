//! `fasthelper concatenate --positions`: many per-position TSVs → one matrix.
//!
//! fastcheck writes one `positions.tsv` per file: a row per position along the
//! read (base counts, depth, quality percentiles). Comparing samples by eye or
//! by plot needs those transposed into positions × samples, which is this
//! module. Two shapes are offered because there are two real uses:
//!
//!   * wide (default): `pos, sampleA, sampleB, …` for one metric — the shape
//!     you plot as N curves;
//!   * long (`--all-metrics`): `sample, pos, metric, value` — tidy, for
//!     computing anything in pandas/R.
//!
//! Positions are unioned across samples: a 150 bp file and a 151 bp file are
//! both representable, and a missing cell is written as empty rather than as a
//! zero that would fake a measurement.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::report::sample_from_path;

/// Columns fastcheck writes, in order. Read from the file header at runtime
/// (not hardcoded as truth): a future fastcheck may add columns, and an
/// unrelated TSV must be rejected rather than misread.
const REQUIRED_COLUMN: &str = "pos";

pub struct SamplePositions {
    pub sample: String,
    pub source: PathBuf,
    pub metrics: Vec<String>,
    /// position (1-based, as written by fastcheck) -> metric values, aligned
    /// with `metrics`
    pub rows: BTreeMap<u32, Vec<f64>>,
}

fn looks_like_positions_header(line: &str) -> bool {
    let cols: Vec<&str> = line.split('\t').collect();
    cols.first() == Some(&REQUIRED_COLUMN) && cols.len() >= 2
}

/// Read one positions TSV. `Ok(None)` means "not a positions table" — the
/// caller skips it and says so, instead of inventing numbers from it.
pub fn read_positions(path: &Path, sample_hint: Option<&str>) -> Result<Option<SamplePositions>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = text.lines();
    let header = match lines.next() {
        Some(h) => h,
        None => return Ok(None),
    };
    if !looks_like_positions_header(header) {
        return Ok(None);
    }
    let metrics: Vec<String> = header
        .split('\t')
        .skip(1)
        .map(|c| c.trim().to_string())
        .collect();

    let sample = sample_hint
        .map(str::to_string)
        .unwrap_or_else(|| sample_from_path(path));

    let mut rows = BTreeMap::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != metrics.len() + 1 {
            bail!(
                "{}: row {:?} has {} columns, expected {} (header: {})",
                path.display(),
                line.chars().take(40).collect::<String>(),
                cells.len(),
                metrics.len() + 1,
                header
            );
        }
        let pos: u32 = cells[0]
            .parse()
            .with_context(|| format!("{}: position {:?} is not a number", path.display(), cells[0]))?;
        let values: Vec<f64> = cells[1..]
            .iter()
            .map(|c| c.parse::<f64>().unwrap_or(f64::NAN))
            .collect();
        rows.insert(pos, values);
    }
    Ok(Some(SamplePositions {
        sample,
        source: path.to_path_buf(),
        metrics,
        rows,
    }))
}

fn is_positions_tsv(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("tsv") || e.eq_ignore_ascii_case("txt"))
}

fn walk(dir: &Path, out: &mut Vec<PathBuf>, depth: usize) {
    if depth > 6 {
        return;
    }
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in rd.flatten() {
        let p = entry.path();
        let name = entry.file_name().to_string_lossy().to_string();
        if p.is_dir() {
            if matches!(
                name.as_str(),
                "target" | ".git" | ".pixi" | "node_modules" | "tmp"
            ) {
                continue;
            }
            walk(&p, out, depth + 1);
        } else if is_positions_tsv(&p) {
            // Never read our own output back as input: its header (`pos, A, B…`)
            // is a valid positions table, so it would silently become another
            // sample and duplicate every column.
            if name == "positions_matrix.tsv" {
                continue;
            }
            out.push(p);
        }
    }
}

fn fmt(v: f64) -> String {
    if v.is_nan() {
        String::new()
    } else if v.fract() == 0.0 {
        format!("{v:.0}")
    } else {
        format!("{v:.4}")
    }
}

/// Concatenate. `metric` selects the wide-matrix column; `all_metrics` selects
/// the tidy shape instead.
pub fn concatenate_positions(
    inputs: &[PathBuf],
    outdir: &Path,
    metric: &str,
    all_metrics: bool,
    to_stdout: bool,
) -> Result<()> {
    let mut files = Vec::new();
    for input in inputs {
        if input.is_dir() {
            walk(input, &mut files, 0);
        } else if input.is_file() {
            files.push(input.clone());
        } else {
            bail!("input not found: {}", input.display());
        }
    }
    files.sort();
    files.dedup();
    if files.is_empty() {
        bail!(
            "no .tsv files found in {}",
            inputs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    let mut samples = Vec::new();
    let mut skipped = 0usize;
    for f in &files {
        match read_positions(f, None)? {
            Some(s) => samples.push(s),
            None => skipped += 1,
        }
    }
    if samples.is_empty() {
        bail!(
            "none of the {} .tsv file(s) is a fastcheck positions table \
             (expected a header starting with '{}')",
            files.len(),
            REQUIRED_COLUMN
        );
    }
    samples.sort_by(|a, b| a.sample.cmp(&b.sample));

    let mut out = String::new();
    if all_metrics {
        out.push_str("sample\tpos\tmetric\tvalue\n");
        for s in &samples {
            for (pos, values) in &s.rows {
                for (name, v) in s.metrics.iter().zip(values.iter()) {
                    out.push_str(&format!("{}\t{}\t{}\t{}\n", s.sample, pos, name, fmt(*v)));
                }
            }
        }
    } else {
        let missing: Vec<&str> = samples
            .iter()
            .filter(|s| !s.metrics.iter().any(|m| m == metric))
            .map(|s| s.sample.as_str())
            .collect();
        if missing.len() == samples.len() {
            bail!(
                "no positions table has a '{}' column. Available here: {}",
                metric,
                samples
                    .first()
                    .map(|s| s.metrics.join(", "))
                    .unwrap_or_default()
            );
        }
        out.push_str("pos");
        for s in &samples {
            out.push('\t');
            out.push_str(&s.sample);
        }
        out.push('\n');

        let mut positions: Vec<u32> = Vec::new();
        for s in &samples {
            positions.extend(s.rows.keys().copied());
        }
        positions.sort_unstable();
        positions.dedup();

        for pos in positions {
            out.push_str(&pos.to_string());
            for s in &samples {
                out.push('\t');
                if let (Some(idx), Some(values)) = (
                    s.metrics.iter().position(|m| m == metric),
                    s.rows.get(&pos),
                ) {
                    out.push_str(&fmt(values[idx]));
                }
            }
            out.push('\n');
        }
    }

    let positions_count = samples
        .iter()
        .map(|s| s.rows.len())
        .max()
        .unwrap_or(0);
    if to_stdout {
        print!("{out}");
    } else {
        std::fs::create_dir_all(outdir)
            .with_context(|| format!("creating {}", outdir.display()))?;
        let path = outdir.join("positions_matrix.tsv");
        std::fs::write(&path, &out).with_context(|| format!("writing {}", path.display()))?;
        println!(
            "fasthelper: concatenated {} sample(s), up to {} position(s) -> {}",
            samples.len(),
            positions_count,
            path.display()
        );
    }
    if skipped > 0 {
        eprintln!(
            "fasthelper: {} file(s) skipped (not positions tables)",
            skipped
        );
    }
    if !all_metrics {
        eprintln!("fasthelper: metric column '{metric}' (use --metric <name> or --all-metrics)");
    }
    Ok(())
}
