mod cli;
mod hasher;
mod parallel;
mod processor;
mod stream;
mod utils;
mod writer;
use crate::hasher::HashType;
use crate::parallel::{execute_deduplication_parallel, execute_paired_deduplication_parallel};
use crate::processor::{
    DedupConfig, execute_deduplication, execute_deduplication_with_abundance,
    execute_paired_deduplication, execute_paired_deduplication_with_abundance,
};
use crate::utils::{
    birthday_problem_square_approximation, count_sequences, estimate_sequence_capacity,
    get_hash_method,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use cli::Cli;
use std::path::Path;
use std::time::Instant;

/// Dispatches single-end deduplication based on the selected hash type.
fn dispatch(
    input_path: &str,
    output_path: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Single-End Mode: 64-bit Hash");
            }
            execute_deduplication::<u64>(input_path, output_path, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Single-End Mode: 128-bit Hash");
            }
            execute_deduplication::<u128>(input_path, output_path, cfg)
        }
    }
}

/// Dispatches single-end abundance-tracking based on the selected hash type.
fn dispatch_abundance(
    input_path: &str,
    output_path: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Single-End Mode (Abundance): 64-bit Hash");
            }
            execute_deduplication_with_abundance::<u64>(input_path, output_path, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Single-End Mode (Abundance): 128-bit Hash");
            }
            execute_deduplication_with_abundance::<u128>(input_path, output_path, cfg)
        }
    }
}

/// Dispatches paired-end deduplication based on the selected hash type.
fn dispatch_paired(
    input_r1: &str,
    input_r2: &str,
    output_r1: &str,
    output_r2: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Paired-End Mode: 64-bit Combined Hash");
            }
            execute_paired_deduplication::<u64>(input_r1, input_r2, output_r1, output_r2, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Paired-End Mode: 128-bit Combined Hash");
            }
            execute_paired_deduplication::<u128>(input_r1, input_r2, output_r1, output_r2, cfg)
        }
    }
}

/// Dispatches paired-end abundance-tracking based on the selected hash type.
fn dispatch_paired_abundance(
    input_r1: &str,
    input_r2: &str,
    output_r1: &str,
    output_r2: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Paired-End Mode (Abundance): 64-bit Combined Hash");
            }
            execute_paired_deduplication_with_abundance::<u64>(
                input_r1, input_r2, output_r1, output_r2, cfg,
            )
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Paired-End Mode (Abundance): 128-bit Combined Hash");
            }
            execute_paired_deduplication_with_abundance::<u128>(
                input_r1, input_r2, output_r1, output_r2, cfg,
            )
        }
    }
}

/// Dispatches single-end parallel deduplication based on the selected hash type.
fn dispatch_threaded(
    input_path: &str,
    output_path: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Single-End Mode (Multi-threaded): 64-bit Hash");
            }
            execute_deduplication_parallel::<u64>(input_path, output_path, cfg)
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Single-End Mode (Multi-threaded): 128-bit Hash");
            }
            execute_deduplication_parallel::<u128>(input_path, output_path, cfg)
        }
    }
}

/// Dispatches paired-end parallel deduplication based on the selected hash type.
fn dispatch_paired_threaded(
    input_r1: &str,
    input_r2: &str,
    output_r1: &str,
    output_r2: &str,
    hash_type: &HashType,
    cfg: &DedupConfig,
) -> Result<(usize, usize)> {
    match hash_type {
        HashType::XXH3_64 => {
            if cfg.verbose {
                println!("Paired-End Mode (Multi-threaded): 64-bit Combined Hash");
            }
            execute_paired_deduplication_parallel::<u64>(
                input_r1, input_r2, output_r1, output_r2, cfg,
            )
        }
        HashType::XXH3_128 => {
            if cfg.verbose {
                println!("Paired-End Mode (Multi-threaded): 128-bit Combined Hash");
            }
            execute_paired_deduplication_parallel::<u128>(
                input_r1, input_r2, output_r1, output_r2, cfg,
            )
        }
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // --- INPUT VALIDATION ---

    if !Path::new(&args.input).exists() {
        bail!("Input file not found: {}", args.input);
    }

    if let Some(ref r2) = args.input_r2 {
        if !Path::new(r2).exists() {
            bail!("R2 input file not found: {}", r2);
        }
    }

    // Resolve output paths: required unless --dry-run
    let output = match (&args.output, args.dry_run) {
        (Some(path), _) => path.clone(),
        (None, true) => String::new(), // dummy, never written in dry-run
        (None, false) => bail!(
            "An output file is required. Use --output (-o) to specify one, or --dry-run to only compute stats."
        ),
    };

    if args.dry_run && args.abundance {
        bail!("--abundance is not compatible with --dry-run");
    }

    // In dry-run mode, output_r2 is optional even when input_r2 is provided
    if args.input_r2.is_some() && !args.dry_run && args.output_r2.is_none() {
        bail!(
            "The --output-r2 (-p) argument is required when --input-r2 (-2) is used (unless --dry-run is set)."
        );
    }

    if args.hash.is_some() && args.threshold != 0.001 {
        eprintln!(
            "Warning: --hash specifies a hash size, so the automatic selection threshold ({}) is ignored.",
            args.threshold
        );
    }

    if args.verbose {
        println!("Primary input file (R1): {}", args.input);
        if args.dry_run && args.output.is_none() {
            println!("Primary output file (R1): (none, dry-run)");
        } else {
            println!("Primary output file (R1): {output}");
        }
        println!("Expected read length: {} bp", args.read_length);
    }

    let cap_input = estimate_sequence_capacity(&args.input, args.read_length)
        .context("Failed to estimate input file capacity")?;

    let cap_output = if args.force || (args.dry_run && args.output.is_none()) {
        0
    } else {
        estimate_sequence_capacity(&output, args.read_length).unwrap_or(0)
    };

    // If --pre-count is enabled, replace the file-size estimate with an exact count
    let pre_count = if args.pre_count {
        if args.verbose {
            println!("Pre-counting sequences in: {}", args.input);
        }
        let exact = count_sequences(&args.input)?;
        if args.verbose {
            println!("  Found {} sequences", exact);
        }
        // For PE, R2 must match R1 so we just count the primary input
        exact
    } else {
        0
    };
    let cap_input = if args.pre_count { pre_count } else { cap_input };

    let total_capacity = cap_input + cap_output;

    let selected_hash_type = match args.hash {
        Some(mode) => HashType::from(mode),
        None => get_hash_method(total_capacity, args.threshold),
    };

    if args.verbose {
        println!(
            "Total estimated hash table capacity: {} fragments",
            total_capacity
        );
    }

    let cfg = DedupConfig {
        force: args.force,
        verbose: args.verbose,
        dry_run: args.dry_run,
        compression_level: args.compression,
        estimated_capacity: total_capacity,
        read_length: args.read_length,
        abundance_sep: args.abundance_sep.clone(),
        threads: args.threads,
        temp_dir: args.temp_dir.clone(),
    };

    if args.abundance && args.threads != 1 {
        eprintln!(
            "Warning: --abundance does not support multi-threading. Running single-threaded."
        );
    }

    let start = Instant::now();

    let (processed, duplicates) = if let Some(input_r2) = &args.input_r2 {
        let output_r2 = args.output_r2.as_deref().unwrap_or("");

        if args.verbose {
            println!("Secondary input file (R2): {input_r2}");
            if output_r2.is_empty() {
                println!("Secondary output file (R2): (none, dry-run)");
            } else {
                println!("Secondary output file (R2): {output_r2}");
            }
            println!("--- Starting Paired-End Processing ---");
        }

        if args.abundance {
            dispatch_paired_abundance(
                &args.input,
                input_r2,
                &output,
                output_r2,
                &selected_hash_type,
                &cfg,
            )?
        } else if cfg.threads == 0 || cfg.threads > 1 {
            dispatch_paired_threaded(
                &args.input,
                input_r2,
                &output,
                output_r2,
                &selected_hash_type,
                &cfg,
            )?
        } else {
            dispatch_paired(
                &args.input,
                input_r2,
                &output,
                output_r2,
                &selected_hash_type,
                &cfg,
            )?
        }
    } else {
        if args.verbose {
            println!("--- Starting Single-End Processing ---");
        }

        if args.abundance {
            dispatch_abundance(&args.input, &output, &selected_hash_type, &cfg)?
        } else if cfg.threads == 0 || cfg.threads > 1 {
            dispatch_threaded(&args.input, &output, &selected_hash_type, &cfg)?
        } else {
            dispatch(&args.input, &output, &selected_hash_type, &cfg)?
        }
    };

    if args.verbose {
        println!(
            "Processed fragments: {}\nDuplicates removed: {:.2}%",
            processed,
            if processed > 0 {
                duplicates as f64 / processed as f64 * 100.0
            } else {
                0.0
            }
        );
        println!(
            "Estimated collisions: {:.2e}%",
            birthday_problem_square_approximation(processed, &selected_hash_type) * 100.0
        );
        println!("Total execution time: {:.2?}", start.elapsed());
    }

    // --- FastSuite common interface: summary line + machine-readable report ---
    let elapsed = start.elapsed();
    let dup_pct = if processed > 0 {
        duplicates as f64 / processed as f64 * 100.0
    } else {
        0.0
    };
    let collisions_pct = birthday_problem_square_approximation(processed, &selected_hash_type) * 100.0;
    let reads_out = processed.saturating_sub(duplicates);

    if !args.quiet {
        println!(
            "fastdedup: {} reads in, {} out, {} duplicates ({:.2}%), {} threads, {:.2?}",
            processed,
            reads_out,
            duplicates,
            dup_pct,
            args.threads,
            elapsed
        );
        // Citation request (MIT stays permissive; this is the community-norm ask).
        eprintln!(
            "If you use fastdedup in published research, please cite: \
             Ribes R., Mandier C., Baniel A. FastDedup: a fast and memory-efficient tool for \
             read deduplication. bioRxiv 2026, doi:10.64898/2026.04.29.721745"
        );
    }

    if let Some(path) = &args.json {
        let report = serde_json::json!({
            "schema": "fastsuite.report.v1",
            "tool": "fastdedup",
            "tool_version": env!("CARGO_PKG_VERSION"),
            "command": std::env::args().collect::<Vec<_>>().join(" "),
            "mode": if args.input_r2.is_some() { "paired" } else { "single" },
            "dry_run": args.dry_run,
            "inputs": { "r1": args.input, "r2": args.input_r2 },
            "outputs": {
                "r1": if args.dry_run { serde_json::Value::Null } else { serde_json::Value::String(output.clone()) },
                "r2": args.output_r2
            },
            "elapsed_s": elapsed.as_secs_f64(),
            "metrics": {
                "reads_in": processed,
                "reads_out": reads_out,
                "reads_duplicate": duplicates,
                "duplicate_rate": dup_pct / 100.0,
                "hash_bits": selected_hash_type.to_num(),
                "hash_capacity_estimate": total_capacity,
                "collision_estimate": collisions_pct / 100.0,
                "threads": args.threads,
                "gz_level": args.compression,
                "abundance": args.abundance
            }
        });
        std::fs::write(path, serde_json::to_string_pretty(&report)? + "\n")
            .with_context(|| format!("writing JSON report to {path}"))?;
    }

    Ok(())
}
