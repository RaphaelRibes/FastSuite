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
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::report::sample_from_path;

/// Columns fastcheck writes, in order. Read from the file header at runtime
/// (not hardcoded as truth): a future fastcheck may add columns, and an
/// unrelated TSV must be rejected rather than misread.
const REQUIRED_COLUMN: &str = "pos";

pub struct SamplePositions {
    pub sample: String,
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

/// The tidy shape: `sample, pos, metric, value`, one row per cell, in input
/// order — no transpose, so it streams straight out.
fn write_tidy(w: &mut dyn std::io::Write, samples: &[SamplePositions]) -> Result<()> {
    w.write_all(b"sample\tpos\tmetric\tvalue\n")?;
    for s in samples {
        for (pos, values) in &s.rows {
            for (name, v) in s.metrics.iter().zip(values.iter()) {
                writeln!(w, "{}\t{}\t{}\t{}", s.sample, pos, name, fmt(*v))?;
            }
        }
    }
    Ok(())
}

/// Write a streamed table to `path`, through an 8 KB buffer.
fn write_file<F>(path: &Path, body: F) -> Result<()>
where
    F: FnOnce(&mut dyn std::io::Write) -> Result<()>,
{
    let file = std::fs::File::create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    let mut w = std::io::BufWriter::new(file);
    body(&mut w)?;
    w.flush()
        .with_context(|| format!("writing {}", path.display()))?;
    Ok(())
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

    let positions_count = samples
        .iter()
        .map(|s| s.rows.len())
        .max()
        .unwrap_or(0);
    if all_metrics {
        // Tidy shape: input order already, so it streams — the table used to be
        // one String (138 MB at 2 k samples for an 88 MB file).
        if to_stdout {
            let stdout = std::io::stdout();
            let mut w = std::io::BufWriter::new(stdout.lock());
            write_tidy(&mut w, &samples)?;
            w.flush()?;
        } else {
            std::fs::create_dir_all(outdir)
                .with_context(|| format!("creating {}", outdir.display()))?;
            let path = outdir.join("positions_matrix.tsv");
            write_file(&path, |w| write_tidy(w, &samples))?;
            println!(
                "fasthelper: concatenated {} sample(s), up to {} position(s) -> {}",
                samples.len(),
                positions_count,
                path.display()
            );
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
        // A per-sample column index, resolved once: the wide table looks the
        // metric up per cell otherwise.
        let columns: Vec<Option<usize>> = samples
            .iter()
            .map(|s| s.metrics.iter().position(|m| m == metric))
            .collect();
        let mut positions: Vec<u32> = Vec::new();
        for s in &samples {
            positions.extend(s.rows.keys().copied());
        }
        positions.sort_unstable();
        positions.dedup();

        let write_wide = |w: &mut dyn std::io::Write| -> Result<()> {
            write!(w, "pos")?;
            for s in &samples {
                write!(w, "\t{}", s.sample)?;
            }
            w.write_all(b"\n")?;
            for pos in &positions {
                write!(w, "{pos}")?;
                for (s, idx) in samples.iter().zip(&columns) {
                    w.write_all(b"\t")?;
                    if let (Some(idx), Some(values)) = (idx, s.rows.get(pos)) {
                        write!(w, "{}", fmt(values[*idx]))?;
                    }
                }
                w.write_all(b"\n")?;
            }
            Ok(())
        };
        if to_stdout {
            let stdout = std::io::stdout();
            let mut w = std::io::BufWriter::new(stdout.lock());
            write_wide(&mut w)?;
            w.flush()?;
        } else {
            std::fs::create_dir_all(outdir)
                .with_context(|| format!("creating {}", outdir.display()))?;
            let path = outdir.join("positions_matrix.tsv");
            write_file(&path, write_wide)?;
            println!(
                "fasthelper: concatenated {} sample(s), up to {} position(s) -> {}",
                samples.len(),
                positions_count,
                path.display()
            );
        }
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
