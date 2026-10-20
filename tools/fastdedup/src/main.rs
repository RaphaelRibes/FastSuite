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
    get_hash_method, path_identity,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use cli::Cli;
use crate::cli::{OutCompression, OutFormat};
use crate::writer::{OutputFormat, StreamCompression};
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

/// Refuses two CLI roles that resolve to the same file.
///
/// `-` is not a path, so a stream side is always accepted here (the stream
/// conflicts are refused separately). Everything else is resolved through the
/// filesystem with [`path_identity`] before the comparison: comparing the
/// strings let `-o ./c.fq -p c.fq` through.
fn refuse_same_file(a: &str, a_role: &str, b: &str, b_role: &str, why: &str) -> Result<()> {
    if crate::stream::is_stream(a) || crate::stream::is_stream(b) {
        return Ok(());
    }
    if path_identity(a) == path_identity(b) {
        bail!("{a_role} ({a}) and {b_role} ({b}) name the same file: {why}");
    }
    Ok(())
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // --- INPUT VALIDATION ---

    let streaming_input = crate::stream::is_stream(&args.input);
    if !streaming_input && !Path::new(&args.input).exists() {
        bail!("Input file not found: {}", args.input);
    }

    if let Some(ref r2) = args.input_r2 {
        if !crate::stream::is_stream(r2) && !Path::new(r2).exists() {
            bail!("R2 input file not found: {r2}");
        }
    }

    // --- D-4: what a single stream cannot carry ---
    //
    // These MUST be evaluated before the "an output file is required" bail further
    // down: the refusals below are about what a stream can carry, and they were
    // previously shadowed by it, so `-1 - -2 -` reported a missing output instead
    // of the real problem (two files through one stdin). Ordered refusals beat a
    // generic one.
    let out1_stream = args.output.as_deref().is_some_and(crate::stream::is_stream);
    let out2_stream = args.output_r2.as_deref().is_some_and(crate::stream::is_stream);

    if streaming_input && args.input_r2.as_deref().is_some_and(crate::stream::is_stream) {
        bail!(
            "Both --input-r1 (-1) and --input-r2 (-2) cannot be stdin: one stream cannot carry two files. Pipe one mate and name the other."
        );
    }
    if streaming_input && args.pre_count {
        bail!(
            "--pre-count needs a second pass over the input, impossible on a stream. Drop --pre-count or pass a file path."
        );
    }
    if streaming_input && args.abundance && args.temp_dir.is_none() {
        bail!(
            "--abundance needs a second pass over the input. With stdin, pass --temp-dir so unique reads are spilled to disk and the stream is read once."
        );
    }
    if out1_stream && out2_stream {
        bail!(
            "--out1 (-o) and --out2 (-p) cannot both be stdout: R1 and R2 would interleave into one corrupted stream. Send one mate to a file."
        );
    }
    // Two independent writers on one path interleave R1 and R2 into a file that
    // is neither valid FASTQ nor FASTA, and exit 0. Refuse before any write.
    // The comparison resolves both paths: `-o ./coll2.fq -p coll2.fq` is one
    // file in two spellings and used to slip past a string comparison.
    if let (Some(out1), Some(out2)) = (args.output.as_deref(), args.output_r2.as_deref()) {
        refuse_same_file(
            out1,
            "--out1 (-o)",
            out2,
            "--out2 (-p)",
            "two writers would interleave R1 and R2 into one file that is neither \
             valid FASTQ nor FASTA. Give each mate its own path.",
        )?;
    }
    // The report is written last, so `--json j.fq -o j.fq` used to replace a
    // finished FASTQ with JSON and exit 0.
    if let Some(json) = args.json.as_deref() {
        for (role, out) in [
            ("--out1 (-o)", args.output.as_deref()),
            ("--out2 (-p)", args.output_r2.as_deref()),
        ] {
            if let Some(out) = out {
                refuse_same_file(
                    json,
                    "--json",
                    out,
                    role,
                    "the report would overwrite the sequences (or the reverse). \
                     Give the report its own path.",
                )?;
            }
        }
    }
    // Reading and writing one file cannot produce correct output: `--force`
    // truncates it before the first read (`-1 f.fq -o f.fq --force` left a
    // 0-byte input), and a resume appends to the file it is reading.
    // `--dry-run` writes nothing, so an output naming an input is harmless.
    if !args.dry_run {
        for (in_role, in_path) in [
            ("--input (-1)", Some(args.input.as_str())),
            ("--input-r2 (-2)", args.input_r2.as_deref()),
        ] {
            let Some(in_path) = in_path else { continue };
            for (out_role, out) in [
                ("--out1 (-o)", args.output.as_deref()),
                ("--out2 (-p)", args.output_r2.as_deref()),
            ] {
                if let Some(out) = out {
                    refuse_same_file(
                        in_path,
                        in_role,
                        out,
                        out_role,
                        "reading and writing one file truncates or appends to the \
                         input. Give the output its own path.",
                    )?;
                }
            }
        }
    }
    if out1_stream && args.json.as_deref() == Some("-") {
        bail!(
            "--out1 (-) and --json (-) cannot both write to stdout. Send the report to a file, or drop one of them."
        );
    }

    // Resolve output paths: required unless --dry-run. Deliberately AFTER the
    // stream refusals above, so a stream-specific error is never reported as a
    // missing output.
    let output = match (&args.output, args.dry_run) {
        (Some(path), _) => path.clone(),
        (None, true) => String::new(), // dummy, never written in dry-run
        (None, false) => bail!(
            "An output file is required. Use --output (-o) to specify one, or --dry-run to only compute stats."
        ),
    };

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
        println!(
            "Primary input file (R1): {}",
            if streaming_input { "(stdin)" } else { &args.input }
        );
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

    let resolved_threads = crate::parallel::resolve_thread_count(args.threads);

    // `--out-compression` is the stream counterpart of the `.gz` extension: a
    // named path used to ignore it and write plain records in silence.
    if args.out_compression == OutCompression::Gz {
        let named_output = [args.output.as_deref(), args.output_r2.as_deref()]
            .into_iter()
            .flatten()
            .any(|path| !crate::stream::is_stream(path));
        if named_output {
            eprintln!(
                "Warning: --out-compression only applies to stream output (-o -/-p -); \
                 a named path infers compression from its extension (.gz). Writing plain records."
            );
        }
    }

    let cfg = DedupConfig {
        force: args.force,
        verbose: args.verbose,
        dry_run: args.dry_run,
        compression_level: args.compression,
        estimated_capacity: total_capacity,
        read_length: args.read_length,
        abundance_sep: args.abundance_sep.clone(),
        // `-j 0` means "auto": resolve it once here so the engine, the summary
        // and the report all announce the pool that actually runs (they used to
        // report the request, so `-j 0` printed "0 threads").
        threads: resolved_threads,
        temp_dir: args.temp_dir.clone(),
        out_format: match args.out_format {
            OutFormat::Fastq => OutputFormat::Fastq,
            OutFormat::Fasta => OutputFormat::Fasta,
        },
        out_compression: match args.out_compression {
            OutCompression::None => StreamCompression::None,
            OutCompression::Gz => StreamCompression::Gz,
        },
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
        // D-4 / CLI-SPEC §0: with `-o -` or `--json -` stdout carries data,
        // so the summary goes to stderr and never pollutes the pipe.
        let summary = format!(
            "fastdedup: {} reads in, {} out, {} duplicates ({:.2}%), {} threads, {:.2?}",
            processed, reads_out, duplicates, dup_pct, resolved_threads, elapsed
        );
        if out1_stream || args.json.as_deref() == Some("-") {
            eprintln!("{summary}");
        } else {
            println!("{summary}");
        }
        // Citation request (MIT stays permissive; this is the community-norm ask).
        eprintln!(
            "If you use fastdedup in published research, please cite: \
             Ribes R., Mandier C., Baniel A. FastDedup: a fast and memory-efficient tool for \
             read deduplication. bioRxiv 2026, doi:10.64898/2026.04.29.721745"
        );
    }

    if let Some(path) = &args.json {
        // Units. `processed`/`reads_out`/`duplicates` above count records: one
        // read in single-end mode, but one PAIR in paired mode. The report's
        // read metrics are reads, so paired mode doubles them — 4 pairs is
        // 8 reads, not 4 — and also publishes the pair counts it used to
        // report under a read name. Keys are added, never renamed: fasthelper
        // aggregates `reads_in`/`reads_out`/`reads_duplicate`.
        let paired = args.input_r2.is_some();
        let pairs = paired.then_some((processed, reads_out));
        let (report_reads_in, report_reads_out, report_reads_duplicate) = if paired {
            (processed * 2, reads_out * 2, duplicates * 2)
        } else {
            (processed, reads_out, duplicates)
        };

        let mut metrics = serde_json::json!({
            "reads_in": report_reads_in,
            "reads_out": report_reads_out,
            "reads_duplicate": report_reads_duplicate,
            "duplicate_rate": dup_pct / 100.0,
            "hash_bits": selected_hash_type.to_num(),
            "hash_capacity_estimate": total_capacity,
            "collision_estimate": collisions_pct / 100.0,
            "threads": resolved_threads,
            "gz_level": args.compression,
            "abundance": args.abundance
        });
        if let Some((pairs_in, pairs_out)) = pairs {
            metrics["pairs_in"] = pairs_in.into();
            metrics["pairs_out"] = pairs_out.into();
        }

        let report = serde_json::json!({
            "schema": "fastsuite.report.v1",
            "tool": "fastdedup",
            "tool_version": env!("CARGO_PKG_VERSION"),
            "command": std::env::args().collect::<Vec<_>>().join(" "),
            "mode": if paired { "paired" } else { "single" },
            "dry_run": args.dry_run,
            "inputs": { "r1": args.input, "r2": args.input_r2 },
            "outputs": {
                // --dry-run writes neither mate: report both as absent.
                "r1": if args.dry_run { serde_json::Value::Null } else { serde_json::Value::String(output.clone()) },
                "r2": if args.dry_run { serde_json::Value::Null } else { serde_json::json!(args.output_r2) }
            },
            "elapsed_s": elapsed.as_secs_f64(),
            "metrics": metrics
        });
        let body = serde_json::to_string_pretty(&report)? + "\n";
        if crate::stream::is_stream(path) {
            // D-4: `-` means the stream, for the report too.
            use std::io::Write;
            let mut out = std::io::stdout();
            out.write_all(body.as_bytes())
                .and_then(|_| out.flush())
                .context("writing the JSON report to stdout")?;
        } else {
            std::fs::write(path, body)
                .with_context(|| format!("writing JSON report to {path}"))?;
        }
    }

    Ok(())
}
