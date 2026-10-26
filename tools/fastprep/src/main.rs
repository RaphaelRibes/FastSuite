//! CLI entry point.

use anyhow::{Context, Result};
use std::time::Instant;

use fastprep::cli::{Cli, Mode};
use fastprep::pipeline::{run_paired, run_single};
use fastprep::stats::{report, report_runtime, pct, to_report_paired_v1, to_report_v1, ReportContext};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// The exact argv, space-joined: enough to reproduce the run without
/// guessing which flags produced the numbers (REPORT-SCHEMA.md rule 5).
fn command_line() -> String {
    std::env::args().collect::<Vec<_>>().join(" ")
}

/// Write the fastsuite.report.v1 wrapper. `--dry-run` writes nothing
/// (CLI-SPEC §1), which is what the `.filter` at both call sites preserves.
fn write_json(
    path: &std::path::Path,
    doc: &serde_json::Value,
) -> anyhow::Result<()> {
    std::fs::write(path, serde_json::to_string_pretty(doc)? + "\n")
        .with_context(|| format!("writing JSON report to {}", path.display()))
}

/// A run that keeps nothing is a configuration error, not a result.
///
/// Every cause of a total drop (a Phred offset read the wrong way, a floor
/// above the read length, a mean-Q gate on a degraded library) leaves an empty
/// output file behind, and an empty file plus exit 0 is indistinguishable from
/// a successful run to every downstream step. Fail loudly instead; `--allow-total-drop`
/// is the explicit opt-out for a caller that really wants the empty file.
///
/// `reads_written` is every record the run actually wrote, not just the pairs:
/// a paired run whose R2 all fails the quality gate still writes R1 singletons
/// when `--singletons1` is given, and calling that an empty output was wrong.
/// Create the parent directory of an output path.
///
/// The suite's other tools disagree about this (fastmaps creates it, the rest
/// failed with ENOENT), and the documented chain assumes `-o out/x.fq.gz` works
/// on a fresh tree. Creating it is the friendlier half of the choice; a dry run
/// touches nothing.
fn ensure_parent(path: &std::path::Path) -> anyhow::Result<()> {
    use anyhow::Context;
    let Some(parent) = path.parent() else {
        return Ok(());
    };
    if parent.as_os_str().is_empty() || parent.exists() {
        return Ok(());
    }
    std::fs::create_dir_all(parent)
        .with_context(|| format!("creating output directory {}", parent.display()))
}

fn refuse_total_drop(
    reads_in: u64,
    reads_written: u64,
    allow: bool,
    mode: &str,
    outputs: &str,
) -> Result<()> {
    if allow || reads_in == 0 || reads_written > 0 {
        return Ok(());
    }
    anyhow::bail!("{}", fastprep::pipeline::total_drop_message(reads_in, mode, outputs))
}

/// The `dedup` block of the report: which anchor window the fingerprint used.
///
/// Recorded because the anchor is a *choice*: when both ends are trimmed the
/// window moves inside the read, and a consumer comparing duplicate rates
/// across runs needs to know which anchor produced them.
fn dedup_report(rt: &fastprep::cli::Runtime, params: &fastprep::algorithm::Params) -> Option<serde_json::Value> {
    let cfg = rt.dedup.as_ref()?;
    let spec = fastprep::dedup::AnchorSpec::new(params, cfg.read_len, params.adapter.is_some());
    Some(serde_json::json!({
        "anchor_end": spec.end.label(),
        // The length of the window that was ACTUALLY fingerprinted. When no
        // window is usable `slice()` falls back to the whole read, and the
        // report used to keep printing the nominal (unused) width.
        "anchor_len": if spec.is_usable() { spec.len } else { cfg.read_len },
        "usable": spec.is_usable(),
        "whole_read": !spec.is_usable(),
        "degenerate": spec.degenerate(),
    }))
}

/// The `--min-len` floor in bases, for the read length the run actually saw.
///
/// The floor and the trim budget are one number: everything the floor does not
/// reserve is trimmable, and a read needing more trim than that is dropped
/// rather than trimmed. Printing it in bases is the difference between
/// "0.67L" and a user discovering after the run that a 30 bp artefact on a
/// 150 bp read cost every read in the file.
fn length_floor_line(params: &fastprep::algorithm::Params, s: &fastprep::stats::Stats) -> String {
    if s.reads_in == 0 {
        return String::new();
    }
    let mean_len = (s.bases_in / s.reads_in) as usize;
    let floor = params.min_len.resolve(mean_len);
    format!(
        "    length floor:            {:>12}  ({} bp of a {} bp mean read; {} bp of trim allowed)",
        floor,
        floor,
        mean_len,
        mean_len.saturating_sub(floor),
    )
}

fn main() -> Result<()> {
    // Keep clap's matches alongside the parsed struct: `--auto-params` must
    // know which knobs the user actually typed (ValueSource::CommandLine) so
    // the recommendation cannot overwrite an explicit choice.
    let matches = <Cli as clap::CommandFactory>::command().get_matches();
    let mut cli = <Cli as clap::FromArgMatches>::from_arg_matches(&matches)?;

    // --auto-params is resolved BEFORE finalize(): it either prints the argv
    // and exits (a query) or patches the knobs the run will use.
    if cli.auto_params.is_some() {
        cli.validate_auto_params()?;
        let (model_path, model, rec) = cli.recommend()?;
        if cli.print_args {
            cli.print_recommendation(&model_path, &model, &rec, true)?;
            return Ok(());
        }
        let explicit = fastprep::cli::Explicit::from_matches(&matches);
        cli.print_recommendation(&model_path, &model, &rec, false)?;
        cli.apply_recommendation(&rec, &explicit);
    }

    let (mode, rt) = cli.finalize()?;

    // Configure rayon's global thread pool.
    if rt.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(rt.threads)
            .build_global()
            .ok(); // ignore error if already initialised
    }

    // -q/--quiet gates the progress and end-of-run summary stream only.
    // Warnings and errors are never gated (CLI-SPEC §1 rule 5), so the
    // quiet flag is threaded only into the reporting below.
    let quiet = rt.quiet;

    let start = Instant::now();

    match mode {
        Mode::Single { in1, out1, params } => {
            if !rt.dry_run {
                ensure_parent(&out1)?;
                if let Some(j) = rt.json.as_deref() {
                    ensure_parent(j)?;
                }
                if let Some(m) = rt.metadata.as_deref() {
                    ensure_parent(m)?;
                }
            }
            let s = run_single(&in1, &out1, &params, &rt)?;
            let elapsed = start.elapsed();
            if !quiet {
                eprintln!("fastprep: single-end complete");
                report("R1", &s);
                let floor_line = length_floor_line(&params, &s);
                if !floor_line.is_empty() {
                    eprintln!("{floor_line}");
                }
                report_runtime(elapsed, s.reads_in);
            }
            if let Some(path) = rt.json.as_ref() {
                let ctx = ReportContext {
                    mode: "single",
                    dry_run: rt.dry_run,
                    elapsed,
                    inputs: serde_json::json!({ "r1": in1 }),
                    outputs: serde_json::json!({ "r1": out1 }),
                    dedup: dedup_report(&rt, &params),
                    threads: rt.threads,
                };
                let doc = to_report_v1(env!("CARGO_PKG_VERSION"), &command_line(), &ctx, &s);
                write_json(path, &doc)?;
            }
            refuse_total_drop(
                s.reads_in,
                s.reads_kept,
                rt.allow_total_drop,
                "single-end",
                "the output file",
            )?;
        }
        Mode::Paired {
            in1,
            in2,
            out1,
            out2,
            singletons1,
            singletons2,
            params_r1,
            params_r2,
        } => {
            if !rt.dry_run {
                ensure_parent(&out1)?;
                ensure_parent(&out2)?;
                for p in [singletons1.as_deref(), singletons2.as_deref()].into_iter().flatten() {
                    ensure_parent(p)?;
                }
                if let Some(j) = rt.json.as_deref() {
                    ensure_parent(j)?;
                }
                if let Some(m) = rt.metadata.as_deref() {
                    ensure_parent(m)?;
                }
            }
            let s = run_paired(
                &in1,
                &in2,
                &out1,
                &out2,
                singletons1.as_deref(),
                singletons2.as_deref(),
                &params_r1,
                &params_r2,
                &rt,
            )?;
            let elapsed = start.elapsed();
            if !quiet {
                eprintln!("fastprep: paired-end complete");
                report("R1", &s.r1);
                report("R2", &s.r2);
                let floor_line = length_floor_line(&params_r1, &s.r1);
                if !floor_line.is_empty() {
                    eprintln!("{floor_line}");
                }
                eprintln!("  [pairs]");
                eprintln!(
                    "    both kept:              {:>12}  ({})",
                    s.pairs_out,
                    pct(s.pairs_out, s.r1.reads_in)
                );
                eprintln!(
                    "    both dropped:           {:>12}  ({})",
                    s.pairs_dropped,
                    pct(s.pairs_dropped, s.r1.reads_in)
                );
                eprintln!("    singletons written (R1):{:>12}", s.singletons_r1);
                eprintln!("    singletons written (R2):{:>12}", s.singletons_r2);
                if s.discarded_singletons_r1 > 0 {
                    eprintln!(
                        "    discarded R1 singletons:{:>12}  (pass --singletons1 to save)",
                        s.discarded_singletons_r1
                    );
                }
                if s.discarded_singletons_r2 > 0 {
                    eprintln!(
                        "    discarded R2 singletons:{:>12}  (pass --singletons2 to save)",
                        s.discarded_singletons_r2
                    );
                }
                if s.bases_corrected > 0 {
                    eprintln!(
                        "    bases corrected:        {:>12}  (--pe-correct)",
                        s.bases_corrected
                    );
                }
                report_runtime(elapsed, s.r1.reads_in + s.r2.reads_in);
            }
            if let Some(path) = rt.json.as_ref() {
                let ctx = ReportContext {
                    mode: "paired",
                    dry_run: rt.dry_run,
                    elapsed,
                    inputs: serde_json::json!({ "r1": in1, "r2": in2 }),
                    outputs: serde_json::json!({
                        "r1": out1,
                        "r2": out2,
                        "singletons1": singletons1,
                        "singletons2": singletons2,
                    }),
                    dedup: dedup_report(&rt, &params_r1),
                    threads: rt.threads,
                };
                let doc = to_report_paired_v1(env!("CARGO_PKG_VERSION"), &command_line(), &ctx, &s);
                write_json(path, &doc)?;
            }
            // Records actually written: pairs plus any salvaged singletons.
            // Counting only the pairs called a run that wrote 100 singletons
            // an empty output and exited 1.
            refuse_total_drop(
                s.r1.reads_in,
                s.pairs_out + s.singletons_r1 + s.singletons_r2,
                rt.allow_total_drop,
                "paired-end",
                "the R1/R2 outputs or either singleton file",
            )?;
        }
    }

    Ok(())
}
