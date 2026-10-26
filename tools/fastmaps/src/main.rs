mod cli;
mod encode;
mod file_kind;
mod format;
mod output;
mod report;
mod spiller;
mod varint;
mod worker;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::time::Instant;

use clap::Parser;

use cli::{corners_filename, parse_corners, Args};
use file_kind::{classify_all, default_r1_patterns, default_r2_patterns};
use format::{write_format, FormatOpts};
use output::{write_html, write_png, write_svg};
use spiller::spiller_thread;
use worker::{
    channel_capacity, flush_local, preview_flush_bytes, CornerExtractor, PixelAcc, SpillerMsg,
    WorkerCfg, WorkerLocal,
};

// ── Path identity ────────────────────────────────────────────────────────────
//
// Two questions in this file are really "do these paths name the same file?":
// the `--json` parent test (is the report inside `-o`?) and the output-target
// collision test (would one output of this run overwrite another?). Both used
// to compare `Path` values directly, which is component-wise and therefore
// sees `out` and `./out` as different — and a collision as no collision.

/// Lexically normalise `p` against the current directory.
///
/// Resolves `.` and `..` components without touching the filesystem:
/// `fs::canonicalize` needs the path to exist (every output target is created
/// by this run) and resolving symlinks is not the question being asked.
fn normalise_path(p: &Path) -> PathBuf {
    let abs = if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(p)
    };
    let mut out = PathBuf::new();
    for c in abs.components() {
        match c {
            Component::CurDir => {}
            Component::ParentDir => {
                // `..` above the root stays at the root.
                if !out.pop() {
                    out.push(c.as_os_str());
                }
            }
            other => out.push(other.as_os_str()),
        }
    }
    out
}

/// First pair of output targets that resolve to the same file, if any.
///
/// `targets` is `(knob, path)` so the refusal can name both sides: a report
/// written over `heatmap.data` is a collision between `--json` and
/// `heatmap.data`, not a mystery about one path.
fn first_collision(targets: &[(String, PathBuf)]) -> Option<(String, String, PathBuf)> {
    let mut seen: Vec<(String, PathBuf)> = Vec::with_capacity(targets.len());
    for (label, path) in targets {
        let norm = normalise_path(path);
        if let Some((other, _)) = seen.iter().find(|(_, n)| *n == norm) {
            return Some((other.clone(), label.clone(), norm));
        }
        seen.push((label.clone(), norm));
    }
    None
}

// ── Per-run temp directory ───────────────────────────────────────────────────

/// A unique directory under `--temp-dir` holding this run's spill rows.
///
/// Rows used to be written straight into `--temp-dir` and only removed by the
/// render pass, so a SIGKILLed run left a full index stream behind and the
/// next run with the same directory rendered *those* rows as its own: 650
/// reads produced 1940 membership entries and 193 count-zero pixels at exit 0
/// (M2). A per-run subdirectory makes leftovers unreachable instead of
/// merely unlikely, and `Drop` removes it on success and on error.
struct TempRunDir {
    root: PathBuf,
    dir: PathBuf,
}

/// Best-effort removal of a run's temp directory, and of `--temp-dir` itself
/// when this run left it empty.
///
/// Separate from `TempRunDir` because it has to be cloned into the worker and
/// spiller threads: those exit through `process::exit`, which runs no
/// destructors, so the cleanup they need must be callable explicitly.
#[derive(Clone)]
struct TempCleanup {
    root: PathBuf,
    dir: PathBuf,
}

impl TempCleanup {
    fn run(&self) {
        let _ = std::fs::remove_dir_all(&self.dir);
        // Only when empty: a `--temp-dir` the user already owned stays theirs.
        let _ = std::fs::remove_dir(&self.root);
    }
}

impl TempRunDir {
    fn new(root: PathBuf) -> std::io::Result<Self> {
        std::fs::create_dir_all(&root)?;
        let pid = std::process::id();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        for attempt in 0..64u32 {
            let dir = root.join(format!("run-{pid}-{nanos}-{attempt}"));
            match std::fs::create_dir(&dir) {
                Ok(()) => return Ok(Self { root, dir }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
        Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("cannot create a unique run directory under {}", root.display()),
        ))
    }

    fn path(&self) -> &Path {
        &self.dir
    }

    fn cleanup(&self) -> TempCleanup {
        TempCleanup {
            root: self.root.clone(),
            dir: self.dir.clone(),
        }
    }
}

impl Drop for TempRunDir {
    fn drop(&mut self) {
        self.cleanup().run();
    }
}

fn main() {
    let started = Instant::now();
    let args = Args::parse();
    // Worker count is read up front: `args` is partially moved below (pattern
    // vectors are consumed by the classifier), so the resolution must happen
    // before those moves.
    let num_workers = args.resolved_threads(args.files.len());

    // ── Validate --preview-reads ──────────────────────────────────────────
    if args.preview_reads > u8::MAX as usize {
        eprintln!("Error: --preview-reads max is {}", u8::MAX);
        std::process::exit(1);
    }

    // ── Validate --buffer-mb ──────────────────────────────────────────────
    // `--buffer-mb 2^44` is 2^64 bytes and wrapped to a zero budget, i.e. a run
    // that "succeeds" with a threshold nothing can stay under (M10).
    let Some(budget_bytes) = args.buffer_mb.checked_mul(1024 * 1024) else {
        eprintln!(
            "Error: --buffer-mb {} is too large: the budget is multiplied by 1 MiB and must \
             fit in a usize (max {} MB)",
            args.buffer_mb,
            usize::MAX / (1024 * 1024)
        );
        std::process::exit(1);
    };

    // ── Validate --json parent directory ───────────────────────────────────
    if let Some(p) = &args.json {
        if p.as_os_str() == "-" {
            eprintln!("Error: fastmaps is file-only (D-4); --json needs a PATH, not '-'");
            std::process::exit(1);
        }
        let parent = p.parent().filter(|d| !d.as_os_str().is_empty());
        if let Some(dir) = parent {
            // A report inside the output directory is fine even when that
            // directory does not exist yet: this run creates it, then the
            // report's parent (review 2026-10-21, item 5 — `-o fresh --json
            // fresh/r.json` used to fail here, before create_dir_all ran).
            //
            // Compared lexically: `-o outB --json ./outB/q/r.json` and
            // `-o ./outC --json outC/q/r.json` are inside `-o` too, but a bare
            // component-wise `starts_with` counts the `.` and refused them
            // (M5).
            if !dir.exists() && !normalise_path(dir).starts_with(normalise_path(&args.output_dir))
            {
                eprintln!(
                    "Error: --json directory {} does not exist",
                    dir.display()
                );
                std::process::exit(1);
            }
        }
    }

    // ── Output paths ──────────────────────────────────────────────────────
    let output_dir = args.output_dir;
    // The grid JSON always lands in the output dir under its historical name:
    // `--json` names the *report*, not the grid (README, CLI-SPEC §1), so -o
    // never loses heatmap.json (review 2026-10-21, item 1).
    let grid_path = output_dir.join(format::GRID_JSON_NAME);
    // `--json PATH` is the fastsuite.report.v1 destination — nothing else.
    let report_path = args.json.clone();

    // Each --corners spec produces a dedicated output file. Parsed up here,
    // before the clobber guard, because that file is one of the outputs the
    // guard must see before anything is created.
    let corner_bounds: Vec<((u32, u32), (u32, u32))> = args
        .corners
        .iter()
        .map(|spec| {
            parse_corners(spec).unwrap_or_else(|| {
                eprintln!(
                    "Error: invalid --corners value {:?}, expected x1,y1:x2,y2",
                    spec
                );
                std::process::exit(1);
            })
        })
        .collect();

    // ── Output targets ────────────────────────────────────────────────────
    // Every file this run would write, built from the same inputs the writers
    // use — the grid JSON, the data/preview streams, PNG and HTML always; the
    // SVG on --svg; one file per --corners box; and the `--json` report
    // wherever it points, inside the output dir or not (review 2026-10-21,
    // item 6: a report outside -o used to be invisible here).
    let mut targets: Vec<(String, PathBuf)> = vec![
        (format::GRID_JSON_NAME.to_string(), grid_path.clone()),
        ("heatmap.data".into(), output_dir.join("heatmap.data")),
        ("heatmap.preview".into(), output_dir.join("heatmap.preview")),
        ("heatmap.png".into(), output_dir.join("heatmap.png")),
        ("heatmap_3d.html".into(), output_dir.join("heatmap_3d.html")),
    ];
    if args.svg {
        targets.push(("heatmap.svg".into(), output_dir.join("heatmap.svg")));
    }
    targets.extend(
        corner_bounds
            .iter()
            .map(|&bounds| {
                let name = corners_filename(bounds);
                (name.clone(), output_dir.join(name))
            }),
    );
    if let Some(p) = &report_path {
        targets.push(("--json".into(), p.clone()));
    }

    // ── Collision guard ───────────────────────────────────────────────────
    // Two of this run's own outputs may not name the same file: `--json
    // oc/heatmap.data` used to exit 0 having replaced the index stream with
    // the report (M1), and two `--corners` boxes that normalise to one file
    // silently dropped the first extraction. Checked before the first
    // `File::create`, and never bypassed by --force: --force is about files
    // that were already there, not about one output of this run eating
    // another.
    if let Some((a, b, path)) = first_collision(&targets) {
        eprintln!(
            "Error: output collision: {} and {} resolve to the same file ({}). Every output of \
             a run needs its own file; change the --json path or the -o directory.",
            a,
            b,
            path.display()
        );
        std::process::exit(1);
    }

    // ── Clobber guard ─────────────────────────────────────────────────────
    // Checked before the first `File::create`, so a refused run leaves the
    // directory untouched.
    if !args.force {
        if let Some((_, p)) = targets.iter().find(|(_, p)| p.exists()) {
            eprintln!(
                "Error: {} already exists. Pass --force to overwrite.",
                p.display()
            );
            std::process::exit(1);
        }
    }
    std::fs::create_dir_all(&output_dir).unwrap_or_else(|e| {
        eprintln!("Cannot create output dir {}: {e}", output_dir.display());
        std::process::exit(1);
    });
    // A report parent inside the output dir may be a nested directory this run
    // creates (`-o out --json out/qc/r.json`).
    if let Some(dir) = report_path
        .as_deref()
        .and_then(|p| p.parent())
        .filter(|d| !d.as_os_str().is_empty() && !d.exists())
    {
        std::fs::create_dir_all(dir).unwrap_or_else(|e| {
            eprintln!("Cannot create --json directory {}: {e}", dir.display());
            std::process::exit(1);
        });
    }

    // ── Temp directory ────────────────────────────────────────────────────
    let temp_root: PathBuf = args
        .temp_dir
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("./tmp_{}", std::process::id())));
    let temp_run = TempRunDir::new(temp_root.clone()).unwrap_or_else(|e| {
        eprintln!("Cannot create temp dir under {}: {e}", temp_root.display());
        std::process::exit(1);
    });
    // Everything below writes into the per-run subdirectory, never the root,
    // so a `--temp-dir` reused after a hard kill cannot feed this run the dead
    // run's rows (M2). `temp_run`'s Drop removes it on success and on error.
    let temp_dir: PathBuf = temp_run.path().to_path_buf();
    let temp_cleanup = temp_run.cleanup();

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
    // One output file per box (the boxes were parsed before the clobber
    // guard); created only now that the guard has passed.
    let corner_extractors: Vec<CornerExtractor> = corner_bounds
        .iter()
        .map(|&bounds| {
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
    // Resolved in main()'s prologue (see `num_workers` above): -j/--threads
    // overrides everything; 0 (or absent) reads $FASTSUITE_THREADS, else all
    // cores. Either way the pool never exceeds the file count (no point in
    // workers with no files).

    let threshold_bytes = budget_bytes / num_workers.max(1);

    let cfg = WorkerCfg {
        preview_n: args.preview_reads,
        threshold_bytes,
        preview_flush_bytes: preview_flush_bytes(threshold_bytes),
        data_mode: args.data_mode,
    };

    // FM-1: a non-exact data mode must be visible in the run summary, not only
    // in the JSON — the user has to know heatmap.data is a sample even when
    // they never open it. Printed unconditionally (like the citation note):
    // `--quiet` silences the progress lines, not a caveat about the outputs.
    if !args.data_mode.is_exact() {
        eprintln!(
            "fastmaps: --data-mode {} → heatmap.data holds {} reads per (pixel, file) \
             (data_reads_per_pixel={}, data_complete=false); counts, PNG and corner \
             metrics stay exact",
            args.data_mode,
            args.data_mode.reads_per_pixel(),
            args.data_mode.reads_per_pixel(),
        );
    }

    // ── Work queue ────────────────────────────────────────────────────────
    let work_queue: Vec<(u16, String)> = files
        .iter()
        .enumerate()
        .map(|(i, p)| (i as u16, p.clone()))
        .collect();
    let queue = Arc::new(Mutex::new(work_queue));
    let done  = Arc::new(AtomicUsize::new(0));

    // ── Spiller thread ────────────────────────────────────────────────────
    // Bounded: one message costs at most the index threshold plus the preview
    // flush threshold, and the queue may hold one budget's worth of them. An
    // unbounded channel let the queue grow with however far the spiller fell
    // behind, so peak RSS scaled with the *file count* (~190 MB for 64 files
    // against the documented 22 MB for 2) whatever --buffer-mb said (M9).
    let per_message = threshold_bytes.saturating_add(cfg.preview_flush_bytes);
    let cap = channel_capacity(budget_bytes, per_message);
    let (spill_tx, spill_rx) = mpsc::sync_channel::<SpillerMsg>(cap);
    let spiller_cleanup = temp_cleanup.clone();
    let spiller_handle = thread::spawn(move || match spiller_thread(spill_rx, temp_dir) {
        Ok(state) => state,
        // The spiller is the channel's only receiver: returning an error would
        // leave every worker blocked on a full queue (or panicking on a dead
        // sender), so a spill failure stops the run here — reported, not a
        // panic, and with the per-run temp directory removed (M8/M2).
        Err(e) => {
            eprintln!("Error writing spill files: {e}");
            spiller_cleanup.run();
            std::process::exit(1);
        }
    });

    // ── Worker threads ────────────────────────────────────────────────────
    let mut handles = Vec::with_capacity(num_workers);

    for _ in 0..num_workers {
        let queue       = Arc::clone(&queue);
        let done        = Arc::clone(&done);
        let tx          = spill_tx.clone();
        let corners     = corner_extractors.clone();
        let kinds       = Arc::clone(&file_kinds);
        let cleanup     = temp_cleanup.clone();
        let cfg         = cfg;

        handles.push(thread::spawn(move || {
            let mut state = WorkerLocal::new();

            loop {
                let item = queue.lock().unwrap().pop();
                let Some((file_id, path)) = item else { break };

                if let Err(e) = worker::process_fastq(
                    &path, file_id, &mut state, &tx, &cfg, &corners, &kinds,
                ) {
                    // `process::exit` skips destructors, so the per-run temp
                    // directory is removed here rather than left for the next
                    // run to find (M2).
                    eprintln!("Error reading {path}: {e}");
                    cleanup.run();
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
                no_call:         state.no_call,
                no_call_r1:      state.no_call_r1,
                no_call_r2:      state.no_call_r2,
                poly_tail:       state.poly_tail,
                min_seq_len:     state.min_seq_len,
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
        drop(temp_run);
        std::process::exit(1);
    }

    if !args.quiet {
        println!(
            "{} unique coordinates, {} total reads",
            state.counts.len(),
            state.total_reads
        );
    }
    // Citation request (MIT stays permissive; this is the community-norm ask).
    eprintln!(
        "If you use fastmaps in published research, please cite: \
         Ribes R. FastSuite — fastmaps. https://github.com/RaphaelRibes/FastSuite"
    );

    // Build u32-keyed map for PNG/HTML before write_format consumes state.
    let counts_u32: HashMap<(u32, u32), PixelAcc> = state
        .counts
        .iter()
        .map(|(&(x, y), &v)| ((x as u32, y as u32), v))
        .collect();

    let fmt_opts = FormatOpts {
        files: &files,
        file_kinds: &file_kinds,
        preview_n: args.preview_reads,
        data_mode: args.data_mode,
        output_dir: &output_dir,
        report_path: report_path.as_deref(),
        elapsed_s: started.elapsed().as_secs_f64(),
        threads: num_workers,
        quiet: args.quiet,
    };

    if let Err(e) = write_format(state, &fmt_opts) {
        eprintln!("Error writing heatmap files: {e}");
        drop(temp_run);
        std::process::exit(1);
    }

    if let Err(e) = write_png(&counts_u32, &output_dir, args.quiet) {
        eprintln!(
            "Error writing {}: {e}",
            output_dir.join("heatmap.png").display()
        );
        drop(temp_run);
        std::process::exit(1);
    }
    if let Err(e) = write_html(&counts_u32, &output_dir, args.quiet) {
        eprintln!(
            "Error writing {}: {e}",
            output_dir.join("heatmap_3d.html").display()
        );
        drop(temp_run);
        std::process::exit(1);
    }

    if args.svg {
        if let Err(e) = write_svg(&counts_u32, &output_dir, args.quiet) {
            eprintln!(
                "Error writing {}: {e}",
                output_dir.join("heatmap.svg").display()
            );
            drop(temp_run);
            std::process::exit(1);
        }
    }

    // Report corner extraction counts.
    for ce in &corner_extractors {
        let n: u64 = counts_u32
            .iter()
            .filter(|(&(x, y), _)| x >= ce.x1 && x <= ce.x2 && y >= ce.y1 && y <= ce.y2)
            .map(|(_, acc)| acc.count)
            .sum();
        let fname = corners_filename(((ce.x1, ce.y1), (ce.x2, ce.y2)));
        if !args.quiet {
            println!("Saved {fname}  ({n} reads extracted)");
        }
    }

    // `temp_run` drops here: the per-run subdirectory is removed with it (and
    // `--temp-dir` itself only when this run left it empty).
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalise_folds_dot_and_dotdot() {
        let cwd = std::env::current_dir().unwrap();
        assert_eq!(normalise_path(Path::new("./a/b")), cwd.join("a/b"));
        assert_eq!(normalise_path(Path::new("a/./b")), cwd.join("a/b"));
        assert_eq!(normalise_path(Path::new("a/x/../b")), cwd.join("a/b"));
        // Absolute paths never depend on the cwd.
        assert_eq!(normalise_path(Path::new("/x/./y")), PathBuf::from("/x/y"));
        assert_eq!(normalise_path(Path::new("/a/b/../c")), PathBuf::from("/a/c"));
    }

    #[test]
    fn collision_sees_through_dot_prefixes() {
        let t = |pairs: &[(&str, &str)]| -> Vec<(String, PathBuf)> {
            pairs
                .iter()
                .map(|(l, p)| (l.to_string(), PathBuf::from(p)))
                .collect()
        };
        let hit = first_collision(&t(&[
            ("heatmap.data", "oc/heatmap.data"),
            ("--json", "./oc/heatmap.data"),
        ]))
        .expect("the same file spelled two ways is a collision");
        assert_eq!(hit.0, "heatmap.data");
        assert_eq!(hit.1, "--json");

        assert!(
            first_collision(&t(&[
                ("heatmap.json", "oc/heatmap.json"),
                ("--json", "other/r.json"),
            ]))
            .is_none(),
            "distinct files are not a collision"
        );
    }
}
