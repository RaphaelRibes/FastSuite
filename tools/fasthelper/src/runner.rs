//! Runs the FastSuite chain from one TOML config and leaves a manifest behind.
//!
//! Order and data flow:
//!
//!   inputs (or fastsim output)
//!     ├─ fastdedup  --json <work>/fastdedup.json
//!     ├─ fastprep   --json <work>/fastprep.json
//!     └─ fastmaps   -o    <work>/maps/            (heatmap.json)
//!   then fasthelper's own reporter over <work>.
//!
//! Each step writes its stdout/stderr to <outdir>/logs/<sample>.<tool>.log and
//! its own JSON report into the run directory, so a run is fully auditable.

use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use crate::adapter_guard::{self, Policy, Verdict};
use crate::cli::{DoctorArgs, RunArgs};
use crate::config::Config;
use crate::report;

/// Steps that trim adapter off the reads. The residual-adapter guard runs only
/// around these, so a chain with no trimming step pays nothing.
const TRIMMING_STEPS: &[&str] = &["fastprep"];

fn is_trimming(tool: &str) -> bool {
    TRIMMING_STEPS.contains(&tool)
}

/// Every step name the runner can execute. `check` validates configs against
/// this same list, so a config `run` accepts is never rejected by `check`.
pub const KNOWN_STEPS: &[&str] = &["fastdedup", "fastprep", "fastmaps", "fastcheck", "fastsim"];

/// Tools whose CLI accepts `--force` (fastcheck has none), so `run --force` can
/// be forwarded to them instead of only deleting their outputs.
const FORCE_STEPS: &[&str] = &["fastdedup", "fastprep", "fastmaps", "fastsim"];

/// The one predicate `run` and `check` share for a step name.
pub fn is_known_step(step: &str) -> bool {
    KNOWN_STEPS.contains(&step)
}

pub struct ToolSet {
    pub check: Option<PathBuf>,
    pub dedup: Option<PathBuf>,
    pub prep: Option<PathBuf>,
    pub maps: Option<PathBuf>,
    pub sim: Option<PathBuf>,
}

impl ToolSet {
    pub fn get(&self, tool: &str) -> Option<&PathBuf> {
        match tool {
            "fastcheck" => self.check.as_ref(),
            "fastdedup" => self.dedup.as_ref(),
            "fastprep" => self.prep.as_ref(),
            "fastmaps" => self.maps.as_ref(),
            "fastsim" => self.sim.as_ref(),
            _ => None,
        }
    }
}

/// Find a tool binary: explicit tools_dir (release then debug), then $PATH.
pub fn resolve(tool: &str, tools_dir: Option<&Path>) -> Option<PathBuf> {
    if let Some(dir) = tools_dir {
        for profile in ["release", "debug"] {
            let p = dir.join(tool).join("target").join(profile).join(tool);
            if p.is_file() {
                return Some(p);
            }
            let p2 = dir.join(tool).join("target").join(profile).join(format!("{tool}.exe"));
            if p2.is_file() {
                return Some(p2);
            }
        }
        // tools_dir may itself be a bin directory
        let p = dir.join(tool);
        if p.is_file() {
            return Some(p);
        }
    }
    if let Ok(path) = std::env::var("PATH") {
        for entry in path.split(':') {
            if entry.is_empty() {
                continue;
            }
            let p = Path::new(entry).join(tool);
            if p.is_file() {
                return Some(p);
            }
        }
    }
    None
}

pub fn toolset(tools_dir: Option<PathBuf>) -> ToolSet {
    let d = tools_dir.as_deref();
    ToolSet {
        check: resolve("fastcheck", d),
        dedup: resolve("fastdedup", d),
        prep: resolve("fastprep", d),
        maps: resolve("fastmaps", d),
        sim: resolve("fastsim", d),
    }
}

struct StepOutcome {
    argv: Vec<String>,
    exit: i32,
    elapsed: f64,
    log: PathBuf,
}

/// What one sample contributed to the run: its manifest entries (each tagged
/// with its step index for deterministic merging), the JSON reports it wrote,
/// plus any failures, plus the names of any samples whose trimming step tripped
/// the residual-adapter guard.
struct SampleOutcome {
    entries: Vec<(usize, serde_json::Value)>,
    /// (sample it belongs to, report path). The sample name travels with the
    /// path because the path cannot be trusted to carry it: every tool's report
    /// file is named after the tool, not the sample.
    reports: Vec<(String, PathBuf)>,
    failed: Vec<String>,
    adapter_trips: Vec<String>,
}

/// Run one sample through the configured steps, strictly in order.
///
/// Returns the manifest entries (tagged with the step index in `cfg.run.steps`)
/// and the failure labels (`sample` or `sample/tool`). Safe to call from
/// several threads: it only touches paths under `<outdir>/<sample.name>` and
/// its own log files.
fn run_sample(
    sample: &crate::config::Sample,
    cfg: &Config,
    ts: &ToolSet,
    outdir: &Path,
    logs: &Path,
    config_dir: &Path,
    dry_run: bool,
    force: bool,
    threads: usize,
    guard: &Policy,
    resume: &BTreeMap<String, ResumeEntry>,
    versions: &BTreeMap<String, Option<String>>,
    from: Option<&str>,
) -> Result<SampleOutcome> {
    let mut entries: Vec<(usize, serde_json::Value)> = Vec::new();
    let mut reports: Vec<(String, PathBuf)> = Vec::new();
    let mut failed: Vec<String> = Vec::new();
    let mut adapter_trips: Vec<String> = Vec::new();

    let work = outdir.join(&sample.name);
    std::fs::create_dir_all(&work)?;
    let in1 = resolve_rel(&sample.r1, config_dir);
    let in2 = sample.r2.as_ref().map(|p| resolve_rel(p, config_dir));
    if !dry_run && !in1.is_file() {
        eprintln!(
            "fasthelper: skipping {}: input {} not found",
            sample.name,
            in1.display()
        );
        failed.push(sample.name.clone());
        return Ok(SampleOutcome {
            entries,
            reports,
            failed,
            adapter_trips,
        });
    }

    println!(
        "fasthelper: sample {} ({} steps)",
        sample.name,
        cfg.run.steps.len()
    );
    let mut cur1 = in1.clone();
    let mut cur2 = in2.clone();
    let mut to_delete: Vec<PathBuf> = Vec::new();
    let mut from_reached = false;

    for (step_idx, tool) in cfg.run.steps.iter().enumerate() {
        let step = cfg.step(tool);
        if !step.enabled {
            continue;
        }
        let Some(bin) = ts.get(tool).cloned() else {
            // A step whose binary does not resolve must not vanish silently:
            // the run used to continue on the previous step's input and exit 0
            // while `check` rejected the very same config.
            eprintln!(
                "fasthelper: {tool} not found — step cannot run (build tools/{tool} or pass --tools-dir)"
            );
            failed.push(format!("{}/{}", sample.name, tool));
            break;
        };
        // Where a fastsim step generates. A step-private prefix, never the
        // [simulate] control's own <outdir>/<simulate.name>/sim.
        let step_prefix = work.join("fastsim").join("sim");
        let (out1, out2, json_path, maps_dir) = match tool.as_str() {
            "fastdedup" => (
                work.join(format!("{}_dedup_R1.fastq.gz", sample.name)),
                cur2.as_ref()
                    .map(|_| work.join(format!("{}_dedup_R2.fastq.gz", sample.name))),
                work.join("fastdedup.json"),
                work.join("maps"),
            ),
            "fastprep" => (
                work.join(format!("{}_clean_R1.fastq.gz", sample.name)),
                cur2.as_ref()
                    .map(|_| work.join(format!("{}_clean_R2.fastq.gz", sample.name))),
                work.join("fastprep.json"),
                work.join("maps"),
            ),
            "fastmaps" => (
                PathBuf::new(),
                None,
                work.join("maps").join("heatmap.json"),
                work.join("maps"),
            ),
            "fastcheck" => (
                PathBuf::new(), // no FASTQ output
                None,
                work.join("fastcheck.json"),
                work.join("maps"), // unused; the aggregator reads the JSON
            ),
            // A generator step writes a fresh pair of FASTQs, so they are this
            // step's outputs (and the next step's inputs). They live in a
            // step-private directory: the [simulate] control already owns
            // <outdir>/<simulate.name>/sim, and a step told to generate into
            // that same prefix tripped fastsim's clobber guard on every run.
            "fastsim" => (
                work.join("fastsim").join(sim_file("R1")),
                Some(work.join("fastsim").join(sim_file("R2"))),
                work.join("fastsim").join("fastsim.json"),
                work.join("fastsim"), // unused
            ),
            // Unreachable: run() validates the step list against KNOWN_STEPS
            // before any sample starts, so this only keeps the match total.
            _ => continue,
        };
        // fastsim does not create its own output directory (it only writes the
        // three artefacts), so the runner must — but never on a dry run.
        if tool == "fastsim" && !dry_run {
            if let Some(dir) = out1.parent() {
                std::fs::create_dir_all(dir)?;
            }
        }
        // Step-private prefix for a fastsim step: the [simulate] control owns
        // <outdir>/<name>/sim, and a step told to write there tripped the
        // clobber guard on every run.
        let sim_prefix = work.join("fastsim").join("sim");
        let argv = argv_for(
            tool,
            &cur1,
            cur2.as_deref(),
            &out1,
            out2.as_deref(),
            &json_path,
            threads,
            cfg.run.gz_level,
            &step.args,
            &maps_dir,
            if tool == "fastsim" {
                Some(&cfg.simulate)
            } else {
                None
            },
            &sim_prefix,
            config_dir,
            force,
        );
        let log = logs.join(format!("{}.{}.log", sample.name, tool));
        // Identity of the inputs this step will consume, and the files it is
        // expected to produce. Both are recorded so a later run can tell
        // whether a resume would reuse a stale output.
        let cur_sig = input_signatures(&cur1, cur2.as_deref());
        let expect_outputs: Vec<PathBuf> = match tool.as_str() {
            "fastdedup" | "fastprep" => std::iter::once(out1.clone())
                .chain(out2.clone())
                .collect(),
            // A generator step's FASTQs ARE its outputs: they are the next
            // step's input, and a resume that skipped the step while they were
            // gone (intermediate cleanup removes them) left the next tool with
            // a missing file. The GT is its third artefact.
            "fastsim" => {
                let mut v = vec![out1.clone(), json_path.clone()];
                if let Some(o2) = &out2 {
                    v.push(o2.clone());
                }
                v.push(PathBuf::from(format!("{}.gt.tsv", sim_prefix.display())));
                v
            }
            // fastmaps/fastcheck write only their JSON report.
            _ => vec![json_path.clone()],
        };
        // Re-running into the same outdir: every tool refuses to clobber an
        // output it is about to rewrite (fastprep and fastmaps guard their own
        // --json too, fastcheck has no --force at all) while cleanup keeps none
        // of them. Clearing the step's own stale artifacts first is the one
        // mechanism that works for every tool, chain position and --from.
        let mut stale: Vec<PathBuf> = expect_outputs.clone();
        if !stale.contains(&json_path) {
            stale.push(json_path.clone());
        }
        if tool == "fastmaps" {
            // The maps dir is this step's own output directory, so everything in
            // it is stale by construction. An explicit filename list is what
            // broke the second run of the documented chain: fastmaps grew a
            // clobber guard for heatmap.preview (and --svg, --corners, the v1
            // wrapper), the list did not, so the guard found a file the cleanup
            // had left behind and refused with exit 1.
            if let Ok(dir) = std::fs::read_dir(maps_dir) {
                for entry in dir.flatten() {
                    let p = entry.path();
                    if p.is_file() {
                        stale.push(p);
                    }
                }
            }
        }
        // Resume: skip this step only when a previous manifest recorded a
        // SUCCESSFUL run (exit 0) of the exact same (argv, binary, version) on
        // the SAME inputs, and the outputs it produced are still byte-for-byte
        // the ones it recorded. --from <step> forces a rerun from that step
        // onward, and --force reruns everything, so once either is in play
        // nothing is resumed.
        let forced = from.map(|f| f == tool).unwrap_or(false);
        if forced {
            from_reached = true;
        }
        let key = format!("{}/{}", sample.name, tool);
        // Resolved once per run (see resolve_versions): a version string cannot
        // change inside one run, and it used to fork once per sample per step.
        let version = versions.get(tool.as_str()).cloned().flatten();
        let prev_out_sig = output_signatures(&expect_outputs);
        if !force && !forced && !from_reached {
            if let Some(prev) = resume.get(&key) {
                let inputs_match = prev.inputs == cur_sig;
                let outputs_match = prev.outputs == prev_out_sig;
                let outputs_present = expect_outputs.iter().all(|p| p.is_file());
                if prev.exit == 0
                    && prev.argv == argv
                    && prev.binary == s(&bin)
                    && prev.version == version
                    && inputs_match
                    && outputs_match
                    && outputs_present
                {
                    println!("fasthelper: resumed {key}");
                    entries.push((
                        step_idx,
                        serde_json::json!({
                            "sample": sample.name,
                            "tool": tool,
                            "binary": prev.binary.clone(),
                            "argv": prev.argv.clone(),
                            "exit": prev.exit,
                            "elapsed_s": 0.0,
                            "log": s(&log),
                            "resumed": true,
                            "tool_version": prev.version.clone(),
                            "inputs": inputs_json(&prev.inputs),
                            "outputs": inputs_json(&prev.outputs),
                        }),
                    ));
                    reports.push((sample.name.clone(), json_path.clone()));
                    // The chain continues from this step's output: the file
                    // the resumed step produced (same argv, same inputs) is
                    // still there — we only get here when it is.
                    if !out1.as_os_str().is_empty() {
                        if !cfg.run.keep_intermediate {
                            if cur1 != in1 {
                                to_delete.push(cur1.clone());
                            }
                            if let Some(c2) = &cur2 {
                                if Some(c2) != in2.as_ref() {
                                    to_delete.push(c2.clone());
                                }
                            }
                        }
                        cur1 = out1;
                        cur2 = out2;
                    }
                    continue;
                }
            }
        }
        if !dry_run {
            for p in &stale {
                if p.is_file() {
                    let _ = std::fs::remove_file(p);
                }
            }
        }
        let outcome = run_step(tool, &bin, &argv, &log, dry_run, None)?;
        println!(
            "  {tool:<10} exit {:<3} {:>7.2}s  log {}",
            outcome.exit,
            outcome.elapsed,
            outcome.log.display()
        );
        // A step that exits 0 without writing the report it was told to write
        // (--json) is a failure: otherwise the run reports success with a
        // header-only table and no measurement behind it.
        if !dry_run && outcome.exit == 0 {
            if let Some(missing) = expect_outputs.iter().find(|p| !p.is_file()) {
                eprintln!(
                    "fasthelper: {} for {} exited 0 but did not write {} — treating the step as failed",
                    tool,
                    sample.name,
                    missing.display()
                );
                failed.push(format!("{}/{}", sample.name, tool));
                break;
            }
        }
        // What this step actually left on disk, recorded so the next run can
        // tell a completed step from a partial output of a failed one.
        let out_sig = output_signatures(&expect_outputs);
        entries.push((
            step_idx,
            serde_json::json!({
                "sample": sample.name,
                "tool": tool,
                "binary": s(&bin),
                "argv": outcome.argv,
                "exit": outcome.exit,
                "elapsed_s": outcome.elapsed,
                "log": s(&outcome.log),
                "tool_version": version.clone(),
                "inputs": inputs_json(&cur_sig),
                "outputs": inputs_json(&out_sig),
            }),
        ));
        reports.push((sample.name.clone(), json_path.clone()));
        if outcome.exit != 0 {
            eprintln!(
                "fasthelper: {} failed for {} (exit {}) — see {}",
                tool,
                sample.name,
                outcome.exit,
                outcome.log.display()
            );
            failed.push(format!("{}/{}", sample.name, tool));
            break;
        }
        // Post-trim residual-adapter guard: measure adapter on the reads entering
        // this trimming step (cur1/cur2) and on the reads it wrote (out1/out2),
        // and warn when it survived (see adapter_guard.rs). Only a trimming step
        // pays for it, so a chain without one is unaffected.
        if guard.enabled && is_trimming(tool) && !dry_run {
            if let Some(check) = ts.check.as_ref() {
                let inputs: Vec<PathBuf> =
                    std::iter::once(cur1.clone()).chain(cur2.clone()).collect();
                let outputs: Vec<PathBuf> =
                    std::iter::once(out1.clone()).chain(out2.clone()).collect();
                let g = adapter_guard::guard_step(
                    check,
                    &inputs,
                    &outputs,
                    logs,
                    &sample.name,
                    guard,
                    &step.args,
                    1,
                )?;
                if let Some(entry) = entries.last_mut() {
                    entry.1["adapter_guard"] = g.to_json();
                }
                match &g.verdict {
                    Verdict::Warn(m) | Verdict::Fail(m) => {
                        eprintln!("fasthelper: adapter guard: sample {}: {m}", sample.name);
                        // Recorded on the manifest either way; only a Fail (under
                        // --fail-on-adapter-residual) makes the run exit 2.
                        adapter_trips.push(sample.name.clone());
                    }
                    Verdict::Silent(_) => {}
                }
            }
        }
        // Only a step that actually wrote FASTQ becomes the next step's input.
        // fastcheck/fastsim/fastmaps write reports, and taking their empty
        // output path forward blanked the chain's input (and dropped the FASTQ
        // the next step still needed).
        if !out1.as_os_str().is_empty() {
            if !cfg.run.keep_intermediate {
                if cur1 != in1 {
                    to_delete.push(cur1.clone());
                }
                if let Some(c2) = &cur2 {
                    if Some(c2) != in2.as_ref() {
                        to_delete.push(c2.clone());
                    }
                }
            }
            cur1 = out1;
            cur2 = out2;
        }
    }
    for p in to_delete {
        if p.is_file() {
            let _ = std::fs::remove_file(&p);
        }
    }
    Ok(SampleOutcome {
        entries,
        reports,
        failed,
        adapter_trips,
    })
}

/// A fastsim step's output file name under its step-private prefix.
/// `-j 0` means "all cores" for the tools that document it, but fastprep
/// requires an explicit count (`--threads must be >= 1`). Resolve the sentinel
/// for the tools that refuse it instead of rewriting it for everyone.
fn threads_for(tool: &str, requested: usize) -> usize {
    if requested == 0 && tool == "fastprep" {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        requested
    }
}

fn sim_file(mate: &str) -> String {
    format!("sim_{mate}.fastq.gz")
}

fn run_step(
    tool: &str,
    bin: &Path,
    argv: &[String],
    log_path: &Path,
    dry_run: bool,
    cwd: Option<&Path>,
) -> Result<StepOutcome> {
    if dry_run {
        println!("[dry-run] {} {}", bin.display(), argv.join(" "));
        return Ok(StepOutcome {
            argv: argv.to_vec(),
            exit: 0,
            elapsed: 0.0,
            log: log_path.to_path_buf(),
        });
    }
    if let Some(dir) = log_path.parent() {
        std::fs::create_dir_all(dir).ok();
    }
    let out = std::fs::File::create(log_path)
        .with_context(|| format!("creating log {}", log_path.display()))?;
    let err = out.try_clone()?;

    let mut cmd = Command::new(bin);
    cmd.args(argv).stdout(out).stderr(err);
    if let Some(dir) = cwd {
        cmd.current_dir(dir);
    }
    let start = Instant::now();
    let status = cmd
        .status()
        .with_context(|| format!("spawning {} ({})", tool, bin.display()))?;
    let elapsed = start.elapsed().as_secs_f64();
    Ok(StepOutcome {
        argv: argv.to_vec(),
        exit: status.code().unwrap_or(-1),
        elapsed,
        log: log_path.to_path_buf(),
    })
}

fn s(p: &Path) -> String {
    p.display().to_string()
}

/// Run `<bin> --version` and return its trimmed stdout (None if it fails).
fn tool_version(bin: &Path) -> Option<String> {
    std::process::Command::new(bin)
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
}

/// `--version` for every tool the run will use, resolved once.
///
/// A version string cannot change inside one run, and the value is needed both
/// by the resume comparison and by every manifest entry: resolving it at the
/// call site forked once per sample per step (30 000 forks ≈ 25 s at 10 k
/// samples, paid on `--dry-run` too).
fn resolve_versions(cfg: &Config, ts: &ToolSet) -> BTreeMap<String, Option<String>> {
    let mut out: BTreeMap<String, Option<String>> = BTreeMap::new();
    for tool in &cfg.run.steps {
        if !cfg.step(tool).enabled || out.contains_key(tool.as_str()) {
            continue;
        }
        let version = ts.get(tool).map(|bin| tool_version(bin)).unwrap_or(None);
        out.insert(tool.clone(), version);
    }
    out
}

/// What a previous run recorded for one `sample/tool` step: the argv it ran,
/// the binary path, that binary's `--version` AT THE TIME, the exit code it
/// returned, and the identity (size + mtime) of each input it consumed and each
/// output it left behind. A step is only resumed when all of these still match
/// today — in particular the previous attempt must have SUCCEEDED: a tool that
/// wrote a partial output and exited 1 used to be resumed as if it had worked.
#[derive(Debug, Clone)]
struct ResumeEntry {
    argv: Vec<String>,
    binary: String,
    version: Option<String>,
    /// The exit code this step returned on the run that recorded it.
    exit: i32,
    /// (path, size, mtime_secs, mtime_nanos) for r1, then r2 when paired.
    inputs: Vec<(String, u64, i64, u32)>,
    /// The same signature for every file the step was expected to produce.
    outputs: Vec<(String, u64, i64, u32)>,
}

/// (path, size, mtime_secs, mtime_nanos) for each of `paths`. Cheap, and an
/// overwrite in place changes size or mtime, so a stale resume is not possible.
/// A path that does not exist contributes nothing, so its absence shows up as a
/// shorter list that can never equal a recorded one.
fn sig_of<'a>(paths: impl Iterator<Item = &'a Path>) -> Vec<(String, u64, i64, u32)> {
    paths
        .filter_map(|p| {
            let md = std::fs::metadata(p).ok()?;
            let mtime = md.modified().ok()?;
            let d = mtime.duration_since(std::time::UNIX_EPOCH).ok()?;
            Some((s(p), md.len(), d.as_secs() as i64, d.subsec_nanos()))
        })
        .collect()
}

/// Identity of one step's inputs, in order.
fn input_signatures(r1: &Path, r2: Option<&Path>) -> Vec<(String, u64, i64, u32)> {
    sig_of(std::iter::once(r1).chain(r2))
}

/// Identity of the files one step is expected to produce.
fn output_signatures(paths: &[PathBuf]) -> Vec<(String, u64, i64, u32)> {
    sig_of(paths.iter().map(|p| p.as_path()))
}

/// Serialise a signature list into a manifest `inputs`/`outputs` array.
fn inputs_json(sig: &[(String, u64, i64, u32)]) -> Vec<serde_json::Value> {
    sig.iter()
        .map(|(p, size, secs, nanos)| {
            serde_json::json!({
                "path": p.clone(),
                "size": *size,
                "mtime_secs": *secs,
                "mtime_nanos": *nanos,
            })
        })
        .collect()
}

/// Read a signature list back from a manifest `inputs`/`outputs` array.
fn sig_from_json(v: &serde_json::Value) -> Vec<(String, u64, i64, u32)> {
    v.as_array()
        .map(|arr| {
            arr.iter()
                .filter_map(|e| {
                    Some((
                        e["path"].as_str()?.to_string(),
                        e["size"].as_u64()?,
                        e["mtime_secs"].as_i64()?,
                        e["mtime_nanos"].as_u64()? as u32,
                    ))
                })
                .collect()
        })
        .unwrap_or_default()
}

/// Read `<outdir>/manifest.json` and build the `sample/tool -> ResumeEntry`
/// map. A missing or unreadable manifest yields an empty map (fresh run).
fn load_resume(outdir: &Path) -> BTreeMap<String, ResumeEntry> {
    let mut map = BTreeMap::new();
    let path = outdir.join("manifest.json");
    let Ok(text) = std::fs::read_to_string(&path) else {
        return map;
    };
    let Ok(doc) = serde_json::from_str::<serde_json::Value>(&text) else {
        return map;
    };
    if let Some(cmds) = doc["commands"].as_array() {
        for c in cmds {
            let sample = c["sample"].as_str().unwrap_or("");
            let tool = c["tool"].as_str().unwrap_or("");
            if sample.is_empty() || tool.is_empty() {
                continue;
            }
            let argv: Vec<String> = serde_json::from_value(c["argv"].clone()).unwrap_or_default();
            let binary = c["binary"].as_str().unwrap_or("").to_string();
            // The version is the one recorded AT RUN TIME, never re-resolved
            // here: re-resolving live made the later comparison tautological,
            // so a binary that changed under the same path was still resumed.
            let version = c["tool_version"].as_str().map(str::to_string);
            // Absent (an older manifest) reads as a failure: the safe
            // direction is to rerun, never to adopt.
            let exit = c["exit"].as_i64().unwrap_or(-1) as i32;
            let inputs = sig_from_json(&c["inputs"]);
            let outputs = sig_from_json(&c["outputs"]);
            map.insert(
                format!("{sample}/{tool}"),
                ResumeEntry {
                    argv,
                    binary,
                    version,
                    exit,
                    inputs,
                    outputs,
                },
            );
        }
    }
    map
}

/// Build the argv for one tool. Every tool receives the shared FastSuite flags
/// (`-j/--threads`, `--gz-level`) plus the per-step extras from the config.
fn argv_for(
    tool: &str,
    in1: &Path,
    in2: Option<&Path>,
    out1: &Path,
    out2: Option<&Path>,
    json: &Path,
    threads: usize,
    gz: u8,
    extra: &[String],
    maps_dir: &Path,
    simulate: Option<&crate::config::Simulate>,
    sim_prefix: &Path,
    config_dir: &Path,
    force: bool,
) -> Vec<String> {
    let mut a: Vec<String> = Vec::new();
    match tool {
        "fastdedup" => {
            a.push("-1".into());
            a.push(s(in1));
            if let Some(r2) = in2 {
                a.push("-2".into());
                a.push(s(r2));
            }
            a.push("-o".into());
            a.push(s(out1));
            if let Some(o2) = out2 {
                a.push("-p".into());
                a.push(s(o2));
            }
            a.push("--json".into());
            a.push(s(json));
            a.push("-j".into());
            a.push(threads_for(tool, threads).to_string());
            a.push("--gz-level".into());
            a.push(gz.to_string());
        }
        "fastprep" => {
            a.push("-1".into());
            a.push(s(in1));
            if let Some(r2) = in2 {
                a.push("-2".into());
                a.push(s(r2));
            }
            a.push("-o".into());
            a.push(s(out1));
            if let Some(o2) = out2 {
                a.push("-p".into());
                a.push(s(o2));
            }
            a.push("--json".into());
            a.push(s(json));
            a.push("-j".into());
            a.push(threads_for(tool, threads).to_string());
            a.push("--gz-level".into());
            a.push(gz.to_string());
        }
        "fastmaps" => {
            a.push("-o".into());
            a.push(s(maps_dir));
            a.push("-j".into());
            a.push(threads_for(tool, threads).to_string());
            a.push(s(in1));
            if let Some(r2) = in2 {
                a.push(s(r2));
            }
        }
        "fastcheck" => {
            // fastcheck takes positional FASTQ files, never -1/-2/-o; it
            // writes only reports. --json is mandatory so the run directory
            // always carries a fastsuite.report.v1 for the aggregator.
            a.push(s(in1));
            if let Some(r2) = in2 {
                a.push(s(r2));
            }
            a.push("--json".into());
            a.push(s(json));
            a.push("-j".into());
            a.push(threads_for(tool, threads).to_string());
        }
        "fastsim" => {
            // simulate-only: the generated FASTQs are inputs, not outputs.
            if let Some(sim) = simulate {
                a.push("--pairs".into());
                a.push(sim.pairs.to_string());
                a.push("--out".into());
                a.push(s(&sim_prefix));
                a.push("--json".into());
                a.push(s(json));
                if let Some(prof) = &sim.profile {
                    a.push("--profile".into());
                    a.push(s(&resolve_rel(prof, config_dir)));
                }
                if let Some(seed) = sim.seed {
                    a.push("--seed".into());
                    a.push(seed.to_string());
                }
                a.push("-j".into());
                a.push(threads_for(tool, threads).to_string());
                a.extend(sim.args.iter().cloned());
            }
        }
        _ => {}
    }
    // `run --force` reaches the tools that can take it: they then overwrite
    // instead of refusing, which is the documented way to redo a run.
    if force && FORCE_STEPS.contains(&tool) {
        a.push("--force".into());
    }
    a.extend(extra.iter().cloned());
    a
}

pub fn run(args: &RunArgs) -> Result<()> {
    let mut cfg = Config::load(&args.config)?;
    if let Some(o) = &args.outdir {
        cfg.run.outdir = o.clone();
    }
    if let Some(t) = args.threads {
        cfg.run.threads = t;
    }
    if args.keep_intermediate {
        cfg.run.keep_intermediate = true;
    }
    let force = args.force;
    // Minute zero: a step name that `check` would reject must fail here too,
    // before a sample directory exists or a tool is forked. A typo used to be
    // skipped with a warning (and, before that, fail only after the earlier
    // steps had already run and written reports).
    for step in &cfg.run.steps {
        if !is_known_step(step) {
            bail!(
                "unknown step '{step}' — supported: {}",
                KNOWN_STEPS.join(", ")
            );
        }
    }
    // --from names a step of THIS run. A typo used to fall through as "no
    // --from", so the run resumed everything and exited 0 while the user
    // believed a step had been re-executed.
    if let Some(from) = &args.from {
        if !cfg.run.steps.iter().any(|s| s == from) {
            bail!(
                "--from '{from}' is not a step of this run — steps: {}",
                cfg.run.steps.join(", ")
            );
        }
    }
    let tools_dir = cfg.tools_dir(args.tools_dir.clone());
    let ts = toolset(tools_dir.clone());
    // Only the tools this run actually needs. Demanding the three FASTQ tools
    // unconditionally made a QC-only chain (`steps = ["fastcheck"]`) fail with
    // "no FastSuite binaries found" while the one binary it needed was there;
    // and a missing one is named before any work starts.
    for step in &cfg.run.steps {
        if cfg.step(step).enabled && ts.get(step).is_none() {
            bail!(
                "step '{step}' is enabled but its binary was not found (looked in {} and $PATH) — \
                 build it with `cargo build --release` in tools/{step}, or pass --tools-dir",
                tools_dir
                    .as_ref()
                    .map(|d| d.display().to_string())
                    .unwrap_or_else(|| "$PATH".into())
            );
        }
    }

    let versions = resolve_versions(&cfg, &ts);

    let guard = Policy {
        enabled: !args.no_adapter_guard,
        threshold: args.adapter_residual_threshold,
        min_reduction: args.adapter_min_reduction,
        fail: args.fail_on_adapter_residual,
    };
    // The guard costs two fastcheck passes per trimming step; say so once when
    // it is enabled but its measurement tool is missing, rather than silently
    // doing nothing.
    let has_enabled_trimming = cfg
        .run
        .steps
        .iter()
        .any(|s| is_trimming(s) && cfg.step(s).enabled);
    if guard.enabled && has_enabled_trimming && ts.check.is_none() && !args.dry_run {
        eprintln!(
            "fasthelper: adapter guard: fastcheck not found — residual-adapter check skipped \
             (build tools/fastcheck or pass --tools-dir)"
        );
    }

    let outdir = cfg.run.outdir.clone();
    let logs = outdir.join("logs");
    std::fs::create_dir_all(&logs)?;
    let config_dir = args
        .config
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    let mut samples = cfg.samples.clone();
    // sim_prefix is only meaningful when fastsim is a step; the argv builder
    // needs it in scope either way. fastsim does not create its own output
    // directory, so when it runs as a step the runner must make it first.
    let sim_prefix = outdir.join(&cfg.simulate.name).join("sim");
    if cfg.run.steps.iter().any(|s| s == "fastsim") {
        if let Some(dir) = sim_prefix.parent() {
            std::fs::create_dir_all(dir)?;
        }
    }

    // Optional synthetic control: generate a known-truth sample and push it
    // through the same chain.
    let mut sim_report: Option<(String, PathBuf)> = None;
    if cfg.simulate.enabled {
        let bin = ts
            .sim
            .clone()
            .context("simulate.enabled is set but no fastsim binary was found")?;
        let simdir = if outdir.is_absolute() {
            outdir.join(&cfg.simulate.name)
        } else {
            std::env::current_dir()
                .context("resolving the working directory")?
                .join(&outdir)
                .join(&cfg.simulate.name)
        };
        std::fs::create_dir_all(&simdir)?;
        let prefix = simdir.join("sim");
        let sim_json = simdir.join("fastsim.json");
        // fastsim refuses to overwrite its own output (sim.gt.tsv and the two
        // FASTQ) and nothing else ever regenerates them, so the second run of
        // the documented workflow died here with "fastsim failed" and exit 1.
        // Reuse what a previous run left; --force regenerates it.
        let sim_outputs: Vec<PathBuf> = vec![
            PathBuf::from(format!("{}.gt.tsv", s(&prefix))),
            PathBuf::from(format!("{}_R1.fastq.gz", s(&prefix))),
            PathBuf::from(format!("{}_R2.fastq.gz", s(&prefix))),
        ];
        let reuse = !force && !args.dry_run && sim_outputs.iter().all(|p| p.is_file());
        if reuse {
            println!(
                "fasthelper: fastsim: reusing {} (pass --force to regenerate)",
                sim_outputs[0].display()
            );
        } else {
            let mut argv: Vec<String> = vec![
                "--pairs".into(),
                cfg.simulate.pairs.to_string(),
                "-o".into(),
                s(&prefix),
                "--json".into(),
                s(&sim_json),
            ];
            if let Some(prof) = &cfg.simulate.profile {
                argv.push("--profile".into());
                argv.push(s(&resolve_rel(prof, &config_dir)));
            }
            if let Some(seed) = cfg.simulate.seed {
                argv.push("--seed".into());
                argv.push(seed.to_string());
            }
            argv.extend(cfg.simulate.args.iter().cloned());
            if force {
                argv.push("--force".into());
            }
            let log = logs.join("fastsim.log");
            let outcome = run_step("fastsim", &bin, &argv, &log, args.dry_run, None)?;
            println!(
                "fasthelper: fastsim {} pairs -> {} (exit {}, {:.2}s)",
                cfg.simulate.pairs,
                prefix.display(),
                outcome.exit,
                outcome.elapsed
            );
            if outcome.exit != 0 {
                bail!("fastsim failed (see {})", log.display());
            }
        }
        sim_report = Some((cfg.simulate.name.clone(), sim_json));
        samples.push(crate::config::Sample {
            name: cfg.simulate.name.clone(),
            r1: PathBuf::from(format!("{}_R1.fastq.gz", s(&prefix))),
            r2: Some(PathBuf::from(format!("{}_R2.fastq.gz", s(&prefix)))),
        });
    }

    // Resume state from the previous manifest (empty on a fresh run). Shared
    // read-only across the worker threads.
    let resume = load_resume(&outdir);

    // Samples run through a shared work queue: min(P, samples) threads each
    // pull the next sample until the queue is empty. Per-sample step order
    // stays strict; only samples are concurrent. The merged manifest is
    // sorted by (sample index, step index) so it is byte-identical to the
    // sequential run's.
    let per_sample_threads = if cfg.run.threads == 0 {
        // 0 is the "all cores" sentinel every tool understands; dividing it
        // would silently turn it into 1.
        0
    } else {
        (cfg.run.threads / args.jobs.max(1)).max(1)
    };
    let queue = std::sync::Mutex::new(std::collections::VecDeque::from(
        samples.iter().enumerate().collect::<Vec<_>>(),
    ));
    let entries: std::sync::Mutex<Vec<(usize, usize, serde_json::Value)>> =
        std::sync::Mutex::new(Vec::new());
    let failures: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());
    let guard_trips: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());
    let report_paths: std::sync::Mutex<Vec<(String, PathBuf)>> = std::sync::Mutex::new(Vec::new());

    std::thread::scope(|scope| {
        for _ in 0..args.jobs.max(1) {
            scope.spawn(|| loop {
                let job = { queue.lock().unwrap().pop_front() };
                let Some((sample_idx, sample)) = job else { break };
                match run_sample(
                    sample,
                    &cfg,
                    &ts,
                    &outdir,
                    &logs,
                    &config_dir,
                    args.dry_run,
                    force,
                    per_sample_threads,
                    &guard,
                    &resume,
                    &versions,
                    args.from.as_deref(),
                ) {
                    Ok(mut o) => {
                        entries
                            .lock()
                            .unwrap()
                            .extend(o.entries.into_iter().map(|(step_idx, m)| (sample_idx, step_idx, m)));
                        report_paths.lock().unwrap().append(&mut o.reports);
                        failures.lock().unwrap().append(&mut o.failed);
                        guard_trips.lock().unwrap().append(&mut o.adapter_trips);
                    }
                    Err(err) => {
                        eprintln!("fasthelper: {} failed: {err}", sample.name);
                        failures.lock().unwrap().push(sample.name.clone());
                    }
                }
            });
        }
    });

    let mut tagged = entries.into_inner().unwrap();
    tagged.sort_by_key(|(sample_idx, step_idx, _)| (*sample_idx, *step_idx));
    let manifest: Vec<serde_json::Value> =
        tagged.into_iter().map(|(_, _, m)| m).collect();
    let mut failed = failures.into_inner().unwrap();
    failed.sort();
    failed.dedup();
    let mut trips = guard_trips.into_inner().unwrap();
    trips.sort();
    trips.dedup();

    // A dry run must not touch the outdir's provenance: it used to overwrite a
    // previous real manifest with entries claiming exit 0 and elapsed 0.0, so
    // `--dry-run` after a paper run destroyed the record of what actually ran.
    if args.dry_run {
        println!("fasthelper: dry-run — nothing executed");
        return Ok(());
    }

    // Manifest + aggregate report over everything THIS run produced.
    let manifest_path = outdir.join("manifest.json");
    // Merge with what is already on disk instead of overwriting it: an entry
    // this run did not produce (a step it never reached, or a rerun that failed
    // before getting there) still describes artifacts sitting in this outdir,
    // and destroying it lost the previous run's provenance. Entries this run
    // produced replace their previous counterpart by (sample, tool); the rest
    // are carried over in their original order. `failed` stays this run's.
    let mut commands = manifest;
    if let Ok(text) = std::fs::read_to_string(&manifest_path) {
        if let Ok(prev) = serde_json::from_str::<serde_json::Value>(&text) {
            let key = |m: &serde_json::Value| {
                (
                    m["sample"].as_str().unwrap_or_default().to_string(),
                    m["tool"].as_str().unwrap_or_default().to_string(),
                )
            };
            let fresh: std::collections::BTreeSet<(String, String)> =
                commands.iter().map(key).collect();
            if let Some(old) = prev["commands"].as_array() {
                commands.extend(old.iter().filter(|m| !fresh.contains(&key(m))).cloned());
            }
        }
    }
    std::fs::write(
        &manifest_path,
        serde_json::to_string_pretty(&serde_json::json!({
            "schema": "fastsuite.manifest.v1",
            "run": cfg.run.name,
            "generated_by": format!("fasthelper {}", env!("CARGO_PKG_VERSION")),
            "threads": cfg.run.threads,
            "gz_level": cfg.run.gz_level,
            "steps": cfg.run.steps,
            "tools": {
                "fastcheck": ts.check.as_ref().map(|p| s(p)),
                "fastdedup": ts.dedup.as_ref().map(|p| s(p)),
                "fastprep": ts.prep.as_ref().map(|p| s(p)),
                "fastmaps": ts.maps.as_ref().map(|p| s(p)),
                "fastsim": ts.sim.as_ref().map(|p| s(p)),
            },
            // The residual-adapter guard: its policy, and which samples tripped it.
            "adapter_guard": {
                "enabled": guard.enabled,
                "threshold": guard.threshold,
                "min_reduction": guard.min_reduction,
                "fail_on_trip": guard.fail,
                "tripped_samples": trips.clone(),
            },
            "commands": commands,
            "failed": failed,
        }))? + "\n",
    )
    .with_context(|| format!("writing {}", manifest_path.display()))?;

    // Aggregate exactly the reports this run produced (the JSON each step wrote
    // plus the simulated control's own report). Walking the outdir instead — as
    // this did — ingested whatever an earlier, differently configured run had
    // left behind, so a table could carry a measurement of a step that did not
    // run.
    let mut report_files = report_paths.into_inner().unwrap();
    // When fastsim runs as a step it already reports that simulation; adding
    // the control's report too would publish the same sample twice.
    if !cfg.run.steps.iter().any(|s| s == "fastsim") {
        report_files.extend(sim_report);
    }
    // Sort/dedupe on the path: the sample label must survive, and a step can
    // legitimately be reached twice (a resumed step records its report again).
    report_files.sort_by(|a, b| a.1.cmp(&b.1));
    report_files.dedup_by(|a, b| a.1 == b.1);
    let reports = report::collect_named(&report_files)?;
    let formats: Vec<String> = vec!["json".into(), "tsv".into(), "md".into(), "html".into()];
    let title = format!("FastSuite report — {}", cfg.run.name);
    let written = report::write_reports(&title, &reports, &outdir, &formats, None)?;
    println!(
        "fasthelper: {} report(s) aggregated -> {}",
        reports.len(),
        written
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
    if !failed.is_empty() {
        bail!("{} step(s) failed: {}", failed.len(), failed.join(", "));
    }
    // The residual-adapter guard asked to be fatal (--fail-on-adapter-residual):
    // the report and manifest are written, then exit 2 — the --fail-on convention
    // ("ran fine, but a threshold tripped", CLI-SPEC §5).
    if guard.fail && !trips.is_empty() {
        for t in &trips {
            eprintln!("fasthelper: adapter guard: {t} tripped --fail-on-adapter-residual");
        }
        eprintln!(
            "fasthelper: {} sample(s) still carry adapter after trimming — report written, exiting 2",
            trips.len()
        );
        std::process::exit(2);
    }
    // Citation request (MIT stays permissive; this is the community-norm ask).
    eprintln!(
        "If you use FastSuite in published research, please cite: \
         Ribes R. FastSuite. https://github.com/RaphaelRibes/FastSuite \
         (fastprep and fastdedup: Mandrillus Project; bioRxiv 2026, doi:10.64898/2026.04.29.721745)"
    );
    Ok(())
}

fn resolve_rel(p: &Path, base: &Path) -> PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        base.join(p)
    }
}

pub fn doctor(args: &DoctorArgs) -> Result<()> {
    let tools_dir = args.tools_dir.clone().or_else(|| {
        std::env::var("FASTSUITE_TOOLS_DIR")
            .ok()
            .filter(|s| !s.is_empty())
            .map(PathBuf::from)
    });
    let ts = toolset(tools_dir.clone());
    println!(
        "fasthelper doctor — tools_dir: {}\n",
        tools_dir
            .map(|d| d.display().to_string())
            .unwrap_or_else(|| "(none: $PATH only)".into())
    );
    let mut table: BTreeMap<&str, (String, String)> = BTreeMap::new();
    for tool in [
        "fastcheck",
        "fastdedup",
        "fastprep",
        "fastmaps",
        "fastsim",
        "fasthelper",
    ] {
        let bin = ts.get(tool).cloned().or_else(|| resolve(tool, None)).or_else(|| {
            // fasthelper is the binary running right now.
            if tool == "fasthelper" {
                std::env::current_exe().ok()
            } else {
                None
            }
        });
        match bin {
            Some(p) => {
                let ver = Command::new(&p)
                    .arg("--version")
                    .output()
                    .ok()
                    .filter(|o| o.status.success())
                    .and_then(|o| String::from_utf8(o.stdout).ok())
                    .map(|s| s.trim().to_string())
                    .unwrap_or_else(|| "no --version".into());
                table.insert(tool, (p.display().to_string(), ver));
            }
            None => {
                table.insert(tool, ("NOT FOUND".into(), "-".into()));
            }
        }
    }
    println!("{:<10} {:<60} {}", "tool", "path", "version");
    for (tool, (path, ver)) in &table {
        println!("{tool:<10} {path:<60} {ver}");
    }
    let missing: Vec<&str> = table
        .iter()
        .filter(|(t, (p, _))| **t != "fasthelper" && p == "NOT FOUND")
        .map(|(t, _)| *t)
        .collect();
    if !missing.is_empty() {
        println!(
            "\nmissing: {}\nbuild with: for t in fastcheck fastdedup fastprep fastmaps fastsim; do (cd tools/$t && cargo build --release); done",
            missing.join(", ")
        );
    }
    Ok(())
}
