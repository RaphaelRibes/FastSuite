mod cli;
mod encode;
mod file_kind;
mod format;
mod output;
mod spiller;
mod varint;
mod worker;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

use clap::Parser;

use cli::{corners_filename, parse_corners, Args};
use file_kind::{classify_all, default_r1_patterns, default_r2_patterns};
use format::write_format;
use output::{write_html, write_png, write_svg};
use spiller::spiller_thread;
use worker::{flush_local, CornerExtractor, SpillerMsg, WorkerCfg, WorkerLocal};

fn main() {
    let args = Args::parse();

    // ── Validate --preview-reads ──────────────────────────────────────────
    if args.preview_reads > u8::MAX as usize {
        eprintln!("Error: --preview-reads max is {}", u8::MAX);
        std::process::exit(1);
    }

    // ── Output directory ──────────────────────────────────────────────────
    let output_dir = args.output_dir;
    std::fs::create_dir_all(&output_dir).unwrap_or_else(|e| {
        eprintln!("Cannot create output dir {}: {e}", output_dir.display());
        std::process::exit(1);
    });

    // ── Temp directory ────────────────────────────────────────────────────
    let temp_dir: PathBuf = args.temp_dir.unwrap_or_else(|| {
        PathBuf::from(format!("./tmp_{}", std::process::id()))
    });
    std::fs::create_dir_all(&temp_dir).unwrap_or_else(|e| {
        eprintln!("Cannot create temp dir {}: {e}", temp_dir.display());
        std::process::exit(1);
    });
    let temp_dir_cleanup = temp_dir.clone();

    // ── Input files ───────────────────────────────────────────────────────
    let files = args.files;
    if files.len() > u16::MAX as usize {
        eprintln!("Error: too many input files (max {})", u16::MAX);
        std::process::exit(1);
    }
    let total_files = files.len();

    // ── R1 / R2 classification ────────────────────────────────────────────
    // If user supplied patterns, they fully override the defaults.
    // `customised` controls the warning behaviour when nothing matches.
    let customised = !args.r1_pattern.is_empty() || !args.r2_pattern.is_empty();
    let r1_patterns = if args.r1_pattern.is_empty() {
        default_r1_patterns()
    } else {
        args.r1_pattern
    };
    let r2_patterns = if args.r2_pattern.is_empty() {
        default_r2_patterns()
    } else {
        args.r2_pattern
    };
    let file_kinds = Arc::new(classify_all(
        &files,
        &r1_patterns,
        &r2_patterns,
        customised,
    ));

    // ── Corner extractors ─────────────────────────────────────────────────
    // Each --corners spec produces a dedicated output file.
    let corner_extractors: Vec<CornerExtractor> = args
        .corners
        .iter()
        .map(|spec| {
            let bounds = parse_corners(spec).unwrap_or_else(|| {
                eprintln!(
                    "Error: invalid --corners value {:?}, expected x1,y1:x2,y2",
                    spec
                );
                std::process::exit(1);
            });
            let fname = corners_filename(bounds);
            let path = output_dir.join(&fname);
            let f = File::create(&path).unwrap_or_else(|e| {
                eprintln!("Cannot create {}: {e}", path.display());
                std::process::exit(1);
            });
            let ((x1, y1), (x2, y2)) = bounds;
            eprintln!(
                "Corner extraction x=[{x1},{x2}] y=[{y1},{y2}] → {fname}"
            );
            CornerExtractor {
                x1, y1, x2, y2,
                writer: Arc::new(Mutex::new(BufWriter::new(f))),
            }
        })
        .collect();

    // ── Worker configuration ──────────────────────────────────────────────
    let num_workers = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(total_files);

    let threshold_bytes = (args.buffer_mb * 1024 * 1024) / num_workers.max(1);

    let cfg = WorkerCfg {
        preview_n: args.preview_reads,
        threshold_bytes,
    };

    // ── Work queue ────────────────────────────────────────────────────────
    let work_queue: Vec<(u16, String)> = files
        .iter()
        .enumerate()
        .map(|(i, p)| (i as u16, p.clone()))
        .collect();
    let queue = Arc::new(Mutex::new(work_queue));
    let done  = Arc::new(AtomicUsize::new(0));

    // ── Spiller thread ────────────────────────────────────────────────────
    let (spill_tx, spill_rx) = mpsc::channel::<SpillerMsg>();
    let spiller_handle = thread::spawn(move || spiller_thread(spill_rx, temp_dir));

    // ── Worker threads ────────────────────────────────────────────────────
    let mut handles = Vec::with_capacity(num_workers);

    for _ in 0..num_workers {
        let queue       = Arc::clone(&queue);
        let done        = Arc::clone(&done);
        let tx          = spill_tx.clone();
        let corners     = corner_extractors.clone();
        let kinds       = Arc::clone(&file_kinds);
        let cfg         = cfg;

        handles.push(thread::spawn(move || {
            let mut state = WorkerLocal::new();

            loop {
                let item = queue.lock().unwrap().pop();
                let Some((file_id, path)) = item else { break };

                if let Err(e) = worker::process_fastq(
                    &path, file_id, &mut state, &tx, &cfg, &corners, &kinds,
                ) {
                    eprintln!("Error reading {path}: {e}");
                    std::process::exit(1);
                }

                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                eprintln!("[{n}/{total_files}] {path}");
            }

            flush_local(&mut state, &tx);
            tx.send(SpillerMsg::Done {
                counts:          state.counts,
                counts_r1:       state.counts_r1,
                counts_r2:       state.counts_r2,
                per_file_counts: state.per_file_counts,
                max_seq_len:     state.max_seq_len,
                n_reads:         state.total_reads,
            })
            .expect("spiller died");
        }));
    }

    // Drop the main sender so the channel closes once all workers finish.
    drop(spill_tx);

    for h in handles {
        h.join().unwrap();
    }

    let state = spiller_handle.join().unwrap();

    if state.total_reads == 0 {
        eprintln!("No reads found.");
        std::process::exit(1);
    }

    println!(
        "{} unique coordinates, {} total reads",
        state.counts.len(),
        state.total_reads
    );

    // Build u32-keyed map for PNG/HTML before write_format consumes state.
    let counts_u32: HashMap<(u32, u32), (u64, f64)> = state
        .counts
        .iter()
        .map(|(&(x, y), &v)| ((x as u32, y as u32), v))
        .collect();

    if let Err(e) = write_format(state, &files, &file_kinds, args.preview_reads, &output_dir) {
        eprintln!("Error writing heatmap files: {e}");
        std::process::exit(1);
    }

    write_png(&counts_u32, &output_dir);
    write_html(&counts_u32, &output_dir);

    if args.svg {
        write_svg(&counts_u32, &output_dir);
    }

    // Report corner extraction counts.
    for ce in &corner_extractors {
        let n: u64 = counts_u32
            .iter()
            .filter(|(&(x, y), _)| x >= ce.x1 && x <= ce.x2 && y >= ce.y1 && y <= ce.y2)
            .map(|(_, (c, _))| c)
            .sum();
        let fname = corners_filename(((ce.x1, ce.y1), (ce.x2, ce.y2)));
        println!("Saved {fname}  ({n} reads extracted)");
    }

    // Best-effort temp dir cleanup (only removes when empty, which it is on success).
    let _ = std::fs::remove_dir(&temp_dir_cleanup);
}
