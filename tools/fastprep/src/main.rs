//! CLI entry point.

use anyhow::{Context, Result};
use clap::Parser;
use std::time::Instant;

use fastprep::cli::{Cli, Mode};
use fastprep::pipeline::{run_paired, run_single};
use fastprep::stats::{report, report_runtime, pct};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> Result<()> {
    let cli = Cli::parse();
    let (mode, rt) = cli.finalize()?;

    // Configure rayon's global thread pool.
    if rt.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(rt.threads)
            .build_global()
            .ok(); // ignore error if already initialised
    }

    let start = Instant::now();

    match mode {
        Mode::Single { in1, out1, params } => {
            let s = run_single(&in1, &out1, &params, &rt)?;
            let elapsed = start.elapsed();
            eprintln!("fastprep: single-end complete");
            report("R1", &s);
            report_runtime(elapsed, s.reads_in);
            if let Some(path) = rt.json.as_ref().filter(|_| !rt.dry_run) {
                std::fs::write(path, serde_json::to_string_pretty(&s)? + "\n")
                    .with_context(|| format!("writing JSON stats to {}", path.display()))?;
            }
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
            eprintln!("fastprep: paired-end complete");
            report("R1", &s.r1);
            report("R2", &s.r2);
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
            report_runtime(elapsed, s.r1.reads_in + s.r2.reads_in);
            if let Some(path) = rt.json.as_ref().filter(|_| !rt.dry_run) {
                std::fs::write(path, serde_json::to_string_pretty(&s)? + "\n")
                    .with_context(|| format!("writing JSON stats to {}", path.display()))?;
            }
        }
    }

    Ok(())
}
