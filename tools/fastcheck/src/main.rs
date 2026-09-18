//! fastcheck — exact, one-pass FASTQ QC.
//!
//! One pass over the input, dense integer counters, no temp files. BGZF input
//! is inflated block-parallel; multiple files are scanned concurrently; a
//! single file is scanned by one thread whose BGZF pool uses all the threads.
//! Budgets (`--reads`, `--bytes`) stop the source and the report declares the
//! run partial with the exact reason.

mod accum;
mod cli;
mod heatmap;
mod html;
mod polyx;
mod profile;
mod repeated;
mod rich;
mod insert;
mod reader;
mod report;

use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use clap::Parser;

use accum::{offset_verdict, Accum};
use cli::Cli;
use reader::Limits;
use report::FileReport;

fn main() -> Result<()> {
    let cli = Cli::parse();
    let requested_offset = cli.requested_offset()?;
    let threads = cli.thread_count();
    let limits = Limits {
        max_reads: cli.reads,
        max_bytes: cli.bytes,
    };

    if cli.profile.is_some() && !(1..=2).contains(&cli.files.len()) {
        anyhow::bail!(
            "--profile needs one input file (single-end) or two (paired-end), got {}",
            cli.files.len()
        );
    }

    let started = Instant::now();
    let files: Vec<PathBuf> = cli.files.clone();
    let max_pos = cli.max_read_len;
    // The per-read profile scan (adapter prefix + fragment GC) is only paid for
    // when a profile was asked for.
    let track_profile = cli.profile.is_some();

    // A single file already has all the parallelism it can use inside the
    // BGZF inflate pool; several files are scanned concurrently.
    let reports: Vec<FileReport> = if files.len() == 1 || threads <= 1 {
        let mut out = Vec::with_capacity(files.len());
        for path in &files {
            out.push(scan_one(path, threads.max(1), limits, max_pos, track_profile)?);
        }
        out
    } else {
        let workers = threads.min(files.len());
        let per_worker_threads = (threads / workers).max(1);
        let mut handles = Vec::with_capacity(workers);
        for w in 0..workers {
            let chunk: Vec<PathBuf> =
                files.iter().skip(w).step_by(workers).cloned().collect();
            handles.push(std::thread::spawn(move || -> Result<Vec<FileReport>> {
                let mut out = Vec::with_capacity(chunk.len());
                for path in &chunk {
                    out.push(scan_one(path, per_worker_threads, limits, max_pos, track_profile)?);
                }
                Ok(out)
            }));
        }
        let mut all = Vec::with_capacity(files.len());
        for h in handles {
            let mut part = h.join().map_err(|_| anyhow::anyhow!("a worker thread panicked"))??;
            all.append(&mut part);
        }
        // Restore input order, so the report is stable across -j.
        all.sort_by_key(|r| files.iter().position(|f| f == &r.path).unwrap_or(usize::MAX));
        all
    };

    let mut total = Accum::new(cli.max_read_len);
    let mut stopped: Option<&'static str> = None;
    for r in &reports {
        total.merge(&r.acc);
        stopped = stopped.or(r.scan.stopped);
    }

    // Memory projection, printed before the numbers: what this run's analysis
    // state costs is a function of read length and file count, not of size.
    let projection = total.resident_bytes() * reports.len().max(1);
    eprintln!(
        "fastcheck: accumulator state {:.1} MB for {} file(s) at --max-read-len {}",
        projection as f64 / (1024.0 * 1024.0),
        reports.len(),
        cli.max_read_len
    );

    // Insert size is a property of a PAIR, and the main path scans files
    // independently, so this is a separate opt-in lockstep pass over a sample.
    let insert = if cli.insert_hist {
        if files.len() != 2 {
            anyhow::bail!(
                "--insert-hist needs exactly two input files (got {}): paired-end overlap is \
                 what makes the insert measurable",
                files.len()
            );
        }
        let h = insert::scan_pair(
            &files[0],
            &files[1],
            threads.max(1),
            cli.insert_pairs,
            cli.max_read_len * 2,
            cli.adapter_seq.iter().map(|s| s.as_bytes().to_ascii_uppercase()).collect(),
        )?;
        Some(h)
    } else {
        None
    };

    // Optional processed files: a second scan whose results are stored as the
    // "after" side of a before/after comparison.
    let mut processed_reports: Vec<FileReport> = Vec::new();
    let mut processed_total = Accum::new(cli.max_read_len);
    if let Some(p1) = &cli.processed_1 {
        let mut pfiles = vec![p1.clone()];
        if let Some(p2) = &cli.processed_2 {
            pfiles.push(p2.clone());
        }
        for path in &pfiles {
            processed_reports.push(scan_one(path, threads.max(1), limits, max_pos, track_profile)?);
        }
        for r in &processed_reports {
            processed_total.merge(&r.acc);
        }
    }

    let verdict = offset_verdict(&total.qual_byte, requested_offset);
    let elapsed = started.elapsed().as_secs_f64();

    if let Some(p) = &cli.json {
        report::write_report(
            p,
            &reports,
            &total,
            &verdict,
            elapsed,
            stopped,
            cli.positions_tsv.as_deref(),
            insert.as_ref().map(|h| (h, total.reads / 2)),
            if processed_reports.is_empty() { None } else { Some(&processed_total) },
        )?;
    }
    if let Some(p) = &cli.html {
        html::write_html(p, &reports, &total, verdict.offset, elapsed,
                         insert.as_ref().map(|h| (h, cli.insert_pairs)),
                         if processed_reports.is_empty() { None }
                         else { Some((&processed_reports, &processed_total)) })?;
    }
    if let Some(p) = &cli.positions_tsv {
        report::write_positions_tsv(p, &total, verdict.offset)?;
    }
    if let Some(p) = &cli.profile {
        let name = files
            .first()
            .and_then(|f| f.file_stem())
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "input".to_string());
        let accs: Vec<&Accum> = reports.iter().map(|r| &r.acc).collect();
        profile::write_profile(p, &name, &accs, verdict.offset, insert.as_ref())?;
        eprintln!("fastcheck: wrote fastsim profile to {}", p.display());
    }

    if !cli.quiet {
        let gc = 100.0 * report_gc(&total);
        let mean_q = report_mean_q(&total, verdict.offset);
        let line = format!(
            "fastcheck: {} file(s), {} reads, {} bases, mean length {:.1}, GC {:.2}%, \
             mean read Q {:.2} (Phred+{}{}), {:.2}s{}",
            reports.len(),
            total.reads,
            total.bases,
            if total.reads == 0 {
                0.0
            } else {
                total.bases as f64 / total.reads as f64
            },
            gc,
            mean_q,
            verdict.offset,
            if verdict.certain { "" } else { ", uncertain" },
            elapsed,
            stopped.map(|s| format!(" [{s}]")).unwrap_or_default()
        );
        // With `--json -` stdout carries the report, so the summary goes to
        // stderr: a consumer parsing the output must not have to skip a line.
        if cli.json.as_deref() == Some(std::path::Path::new("-")) {
            eprintln!("{line}");
        } else {
            println!("{line}");
        }
    }
    Ok(())
}

fn scan_one(
    path: &PathBuf,
    gz_threads: usize,
    limits: Limits,
    max_pos: usize,
    track_profile: bool,
) -> Result<FileReport> {
    let mut acc = Accum::new(max_pos);
    acc.track_profile = track_profile;
    let scan = reader::scan_file(path, gz_threads, limits, &mut acc)?;
    Ok(FileReport {
        path: path.clone(),
        acc,
        scan,
    })
}

// The two helpers below exist only for the summary line; the report module
// owns the same computations for the JSON, and they are tested there.
fn report_gc(a: &Accum) -> f64 {
    let acgt: u64 = (0..a.acgt_hist.len()).map(|i| i as u64 * a.acgt_hist[i]).sum();
    let gc: u64 = (0..a.gc_hist.len()).map(|i| i as u64 * a.gc_hist[i]).sum();
    if acgt == 0 {
        0.0
    } else {
        gc as f64 / acgt as f64
    }
}

fn report_mean_q(a: &Accum, offset: u8) -> f64 {
    if a.reads == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for (b, &c) in a.meanq_raw_hist.iter().enumerate() {
        sum += (b as i64 - offset as i64).max(0) as f64 * c as f64;
    }
    sum / a.reads as f64
}
