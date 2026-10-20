//! fastcheck — exact, one-pass FASTQ QC.
//!
//! One pass over the input, dense integer counters, no temp files. BGZF input
//! is inflated block-parallel; multiple files are scanned concurrently; a
//! single file is scanned by one thread whose BGZF pool uses all the threads.
//! Budgets (`--reads`, `--bytes`) stop the source and the report declares the
//! run partial with the exact reason.

mod accum;
mod cli;
mod gated;
mod heatmap;
mod html;
mod polyx;
mod profile;
mod repeated;
mod rich;
mod insert;
mod memory;
mod fingerprint;
mod overrep;
mod adapter;
mod adapter_pos;
mod tile;
mod reader;
mod report;

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;

use accum::{offset_verdict, Accum};
use cli::{Cli, Modules};
use gated::Work;
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

    // Module gating is resolved and validated before a single byte is read: an unknown
    // `--only` name, or an output whose inputs the selected modules will not
    // materialise, must not cost a scan. `None` means no gating flag was given, and the
    // whole run takes the ungated path — today's behaviour, byte for byte.
    let modules = cli.modules()?;
    let gate = if cli.gating_requested() { Some(modules) } else { None };
    let overrep = cli.overrep_mode()?;
    // The per-tile module is opt-in behind --tiles; `modules()` already folded the flag
    // (and any --only selection) into `modules.tile`, so this is the one place that decides
    // whether the tile accumulator exists at all.
    let track_tile = modules.tile;
    if let Some(m) = &gate {
        check_outputs(&cli, m, overrep)?;
    }

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
    // With the dup module off, nothing is fingerprinted: the table is neither
    // projected nor admission-checked, and the memory the check guards is simply not
    // allocated.
    let track_dup = gate.map_or(true, |m| m.dup);

    // The duplication fingerprint set is the only state whose size is a function of the
    // dataset, so it is projected and admission-checked before the scan starts. The
    // estimate is deliberately conservative on compressed input (a 4x ratio under-counts,
    // so the table is under-sized and grows rather than over-allocating): getting it wrong
    // costs a rehash, not correctness. Copied from fastprep's estimate_records (D-1).
    let fp_projected = if track_dup {
        let est_records: usize = cli
            .files
            .iter()
            .map(|p| {
                // `-` is stdin, whose size is unknowable before reading; the estimate then
                // floors at the clamp minimum and the table simply grows as reads arrive.
                let bytes = if p.as_os_str() == "-" {
                    0
                } else {
                    std::fs::metadata(p)
                        .with_context(|| format!("sizing {}", p.display()))?
                        .len()
                };
                // Read length comes from the profile fit's own length histogram later; 150 is
                // the short-read default and only scales the estimate.
                let per_record = (150u64 * 2 + 50).max(1);
                let compressed = p
                    .extension()
                    .and_then(|e| e.to_str())
                    .is_some_and(|e| matches!(e, "gz" | "zst" | "zstd" | "bz2" | "xz"));
                let raw = if compressed { bytes.saturating_mul(4) } else { bytes };
                Ok::<usize, anyhow::Error>((raw / per_record) as usize)
            })
            .sum::<anyhow::Result<usize>>()?;
        let fp_capacity = est_records.clamp(1, memory::MAX_PREALLOC_ENTRIES);
        // The table is `FxHashMap<u64, u32>`: 16 B for the (fingerprint, count) record plus
        // hashbrown's 1 control byte per bucket. The projection is the size the table
        // reaches when it holds `fp_capacity` distinct reads; it is not preallocated to it.
        let entry_bytes = std::mem::size_of::<(u64, u32)>();
        let fp_proj = memory::project(fp_capacity, entry_bytes);
        eprintln!(
            "fastcheck: duplication fingerprint table projected at {:.1} MB if ~{} distinct \
             reads are seen ({} buckets x {} B entry + 1 B control = {} B/bucket); the table \
             grows into that, it is not preallocated",
            memory::mib(fp_proj.bytes),
            fp_capacity,
            fp_proj.buckets,
            entry_bytes,
            entry_bytes + 1,
        );
        Some((fp_capacity, fp_proj))
    } else {
        eprintln!(
            "fastcheck: the dup module is off: no fingerprint table is projected, built or \
             admission-checked, and no read is fingerprinted"
        );
        None
    };

    // An explicit --max-mem overrides discovery, but only DOWNWARDS: it cannot authorise
    // exceeding a cgroup limit. (fastprep's rule, D-1 copy.)
    let discovered = memory::discover();
    let limit_bytes = cli
        .max_mem
        .map(|mb| mb.saturating_mul(1024 * 1024))
        .into_iter()
        .chain(discovered.as_ref().map(|l| l.bytes))
        .min();
    let limit = limit_bytes.map(|bytes| memory::Limit {
        bytes,
        source: discovered
            .as_ref()
            .filter(|l| l.bytes == bytes)
            .map(|l| l.source.clone())
            .unwrap_or(memory::LimitSource::RlimitAs), // an explicit --max-mem is a soft cap
    });

    // --overrep-mode table keeps every unique sequence in RAM (2-bit packed), so it is
    // admission-checked with the same rule as the fingerprint table it builds on. It is
    // checked FIRST: refusing on the mode the user asked for, with the cheaper modes
    // named, is more useful than a refusal about the table they get either way.
    let adapter_set = cli.adapter_set()?;
    let table_proj = if overrep == overrep::OverrepMode::Table {
        let fp_capacity = fp_projected.map_or(0, |(c, _)| c);
        let proj = overrep::table_projection(fp_capacity, 150);
        eprintln!(
            "fastcheck: over-representation table projected at {:.1} MB for ~{} unique sequences",
            memory::mib(proj),
            fp_capacity,
        );
        match memory::admit(proj, limit.as_ref()) {
            memory::Admission::Refuse { projected, limit, source } => anyhow::bail!(
                "--overrep-mode table needs ~{:.1} MB of packed unique sequences\n\
                 but this process is limited to {:.1} MB ({})\n\
\n\
                 --overrep-mode fingerprint   exact multiplicities at 8 B per unique read\n\
                 --overrep-mode second-pass   same memory as fingerprint, plus one extra read\n\
                 --overrep-mode none         duplication levels only, zero cost",
                memory::mib(projected),
                memory::mib(limit),
                source,
            ),
            memory::Admission::ProceedWithWarning { headroom_bytes } => eprintln!(
                "fastcheck: warning: unique-sequence table projected at {:.1} MB, leaving \
                 {:.1} MB headroom",
                memory::mib(proj),
                memory::mib(headroom_bytes),
            ),
            memory::Admission::Proceed => {}
        }
        Some(proj)
    } else {
        None
    };

    // The fingerprint table is only admission-checked when it is going to exist.
    if let Some((_, fp_proj)) = &fp_projected {
        match memory::admit(fp_proj.bytes, limit.as_ref()) {
            memory::Admission::Refuse { projected, limit, source } => {
                // In table mode the sequence table is the larger allocation, so name it.
                if let Some(tp) = table_proj {
                    anyhow::bail!(
                        "the duplication fingerprint table needs ~{:.1} MB ({} buckets x 17 B: \
                         (u64, u32) entry + control byte)\n\
                         and --overrep-mode table needs ~{:.1} MB of packed unique sequences on top\n\
                         but this process is limited to {:.1} MB ({})\n\
\n\
                         --max-mem <MB>                    raise the ceiling (it cannot exceed the cgroup limit)\n\
                         --overrep-mode fingerprint        exact multiplicities at 8 B per unique read\n\
                         --overrep-mode none               duplication levels only, zero cost",
                        memory::mib(projected),
                        fp_proj.buckets,
                        memory::mib(tp),
                        memory::mib(limit),
                        source,
                    )
                }
                anyhow::bail!(
                    "the duplication fingerprint table needs ~{:.1} MB ({} buckets x 17 B: \
                     (u64, u32) entry + control byte)\n\
                     but this process is limited to {:.1} MB ({})\n\
\n\
                     --max-mem <MB>           raise the ceiling (it cannot exceed the cgroup limit)\n\
                     --overrep-mode none      run without the duplication table entirely\n\
                     reread the file with a smaller --reads budget to bound the estimate",
                    memory::mib(projected),
                    fp_proj.buckets,
                    memory::mib(limit),
                    source,
                )
            }
            memory::Admission::ProceedWithWarning { headroom_bytes } => eprintln!(
                "fastcheck: warning: fingerprint table projected at {:.1} MB, leaving only {:.1} MB headroom",
                memory::mib(fp_proj.bytes),
                memory::mib(headroom_bytes),
            ),
            memory::Admission::Proceed => {}
        }
    }

    // A single file already has all the parallelism it can use inside the
    // BGZF inflate pool; several files are scanned concurrently.
    let reports: Vec<FileReport> = if files.len() == 1 || threads <= 1 {
        let mut out = Vec::with_capacity(files.len());
        for path in &files {
            out.push(scan_one(path, threads.max(1), limits, max_pos, track_profile, overrep, &adapter_set, cli.adapter_content, track_tile, gate)?);
        }
        out
    } else {
        let workers = threads.min(files.len());
        let per_worker_threads = (threads / workers).max(1);
        let mut handles = Vec::with_capacity(workers);
        for w in 0..workers {
            let chunk: Vec<PathBuf> =
                files.iter().skip(w).step_by(workers).cloned().collect();
            let adapter_set = adapter_set.clone();
            handles.push(std::thread::spawn(move || -> Result<Vec<FileReport>> {
                let mut out = Vec::with_capacity(chunk.len());
                for path in &chunk {
                    out.push(scan_one(path, per_worker_threads, limits, max_pos, track_profile, overrep, &adapter_set, cli.adapter_content, track_tile, gate)?);
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
    let mut work = Work::default();
    for r in &reports {
        total.merge(&r.acc);
        stopped = stopped.or(r.scan.stopped);
        work.add(&r.scan.work);
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

    // Only when --tiles: the per-tile module's memory, printed so a run's cost is known
    // rather than discovered. Off, there is no accumulator and no line.
    if let Some(t) = total.tile.as_ref() {
        eprintln!(
            "fastcheck: per-tile quality for {} tile(s), {:.1} MB",
            t.tiles.len(),
            t.resident_bytes() as f64 / (1024.0 * 1024.0),
        );
    }

    // Only when a gating flag was given: what ran, and the work the run actually
    // performed. A module that is off contributes nothing to these counters because
    // its loop never ran — that is the claim, and this is the measurement of it.
    if let Some(m) = &gate {
        eprintln!("fastcheck: modules: {}", m.enabled_names().join(", "));
        eprintln!("fastcheck: measurement work: {}", work.summary());
    }

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
            processed_reports.push(scan_one(path, threads.max(1), limits, max_pos, track_profile, overrep, &adapter_set, cli.adapter_content, track_tile, gate)?);
        }
        for r in &processed_reports {
            processed_total.merge(&r.acc);
        }
    }

    let verdict = offset_verdict(&total.qual_byte, requested_offset);
    let elapsed = started.elapsed().as_secs_f64();

    let min_count = cli
        .overrep_min_count
        .unwrap_or_else(|| ((total.reads as f64 * 0.001) as u32).max(10));
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
            overrep,
            min_count,
            threads,
            total.adapter_pos.as_ref().map(|ap| {
                (ap, cli.adapter_error_rate, cli.min_adapter_overlap)
            }),
            gate.as_ref(),
            total.tile.as_ref(),
        )?;
    }
    if let Some(p) = &cli.html {
        html::write_html(p, &reports, &total, verdict.offset, elapsed,
                         insert.as_ref(), stopped,
                         if processed_reports.is_empty() { None }
                         else { Some((&processed_reports, &processed_total)) })?;
    }
    if let Some(p) = &cli.positions_tsv {
        report::write_positions_tsv(
            p,
            &total,
            verdict.offset,
            total.adapter_pos.as_ref(),
            gate.as_ref(),
        )?;
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
        // The summary line carries only what was measured: with the quality module off
        // there is no mean read Q to print, and with gc off no GC. Building it from
        // parts keeps the ungated line byte for byte what it was.
        let mut parts = vec![
            format!("{} file(s)", reports.len()),
            format!("{} reads", total.reads),
            format!("{} bases", total.bases),
            format!(
                "mean length {:.1}",
                if total.reads == 0 {
                    0.0
                } else {
                    total.bases as f64 / total.reads as f64
                }
            ),
        ];
        if gate.map_or(true, |m| m.gc) {
            parts.push(format!("GC {:.2}%", 100.0 * report_gc(&total)));
        }
        if gate.map_or(true, |m| m.quality) {
            parts.push(format!(
                "mean read Q {:.2} (Phred+{}{})",
                report_mean_q(&total, verdict.offset),
                verdict.offset,
                if verdict.certain { "" } else { ", uncertain" },
            ));
        }
        let line = format!(
            "fastcheck: {}, {:.2}s{}",
            parts.join(", "),
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
        // Citation request (MIT stays permissive; this is the community-norm ask).
        eprintln!(
            "If you use fastcheck in published research, please cite: \
             Ribes R. FastSuite — fastcheck. https://github.com/RaphaelRibes/FastSuite"
        );
    }
    Ok(())
}

fn scan_one(
    path: &PathBuf,
    gz_threads: usize,
    limits: Limits,
    max_pos: usize,
    track_profile: bool,
    overrep: overrep::OverrepMode,
    adapter_set: &crate::adapter::AdapterSet,
    adapter_content: bool,
    tiles: bool,
    gate: Option<Modules>,
) -> Result<FileReport> {
    let mut acc = Accum::new(max_pos);
    acc.track_profile = track_profile;
    acc.track_overrep_table = overrep == overrep::OverrepMode::Table;
    acc.adapter_set = adapter_set.clone();
    if adapter_content {
        acc.adapter_pos = Some(crate::adapter_pos::AdapterPos::new(adapter_set));
    }
    // Only when --tiles: the per-tile accumulator exists iff the module runs, so "off"
    // costs nothing and leaves the report unchanged.
    if tiles {
        acc.tile = Some(crate::tile::TileQual::default());
    }
    let scan = reader::scan_file_gated(path, gz_threads, limits, &mut acc, gate)?;
    Ok(FileReport {
        path: path.clone(),
        acc,
        scan,
    })
}

/// Refuse an output whose inputs the selected modules will not materialise.
///
/// Gating means "not measured", and an artifact rendered from an accumulator that was
/// never filled would read as a measurement of zero. The JSON report can say so field
/// by field (an absent key is not a zero) and the per-position TSV travels with its
/// columns, but a rendered HTML page and a fitted profile cannot: those are refused
/// while they are still cheap to refuse, naming the module that would have to run.
/// (`--insert-hist` needs no module — it is a separate lockstep pass over the pair of
/// files, and it builds no accumulator.)
fn check_outputs(cli: &Cli, m: &Modules, overrep: overrep::OverrepMode) -> Result<()> {
    let need = |what: &str, needed: &[(&str, bool)]| -> Result<()> {
        let missing: Vec<&str> =
            needed.iter().filter(|(_, on)| !*on).map(|(name, _)| *name).collect();
        if missing.is_empty() {
            return Ok(());
        }
        anyhow::bail!(
            "{} needs the {} module{}: --only/--no-quality/--no-heatmap disabled {} \
             (add {} to --only, or drop the gating flag)",
            what,
            missing.join(", "),
            if missing.len() > 1 { "s" } else { "" },
            if missing.len() > 1 { "them" } else { "it" },
            if missing.len() > 1 { "them" } else { "it" },
        )
    };
    if cli.html.is_some() {
        need(
            "--html",
            &[
                ("bases", m.bases),
                ("quality", m.quality),
                ("length", m.length),
                ("gc", m.gc),
                ("heatmap", m.heatmap),
                ("dup", m.dup),
                ("polyx", m.polyx),
                ("kmer", m.kmer),
            ],
        )?;
    }
    if cli.profile.is_some() {
        need(
            "--profile",
            &[
                ("bases", m.bases),
                ("quality", m.quality),
                ("length", m.length),
                ("polyx", m.polyx),
                ("dup", m.dup),
            ],
        )?;
    }
    if cli.positions_tsv.is_some() && !(m.bases || m.quality) {
        anyhow::bail!(
            "--positions-tsv needs one of the bases/quality modules: --only/--no-quality \
             selected neither, so no per-position row was materialised"
        );
    }
    if overrep != overrep::OverrepMode::None {
        need("--overrep-mode", &[("dup", m.dup)])?;
    }
    Ok(())
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
