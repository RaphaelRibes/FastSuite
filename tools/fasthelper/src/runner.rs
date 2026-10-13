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

use crate::cli::{DoctorArgs, RunArgs};
use crate::config::Config;
use crate::report;

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
/// with its step index for deterministic merging) plus any failures.
struct SampleOutcome {
    entries: Vec<(usize, serde_json::Value)>,
    failed: Vec<String>,
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
    threads: usize,
    resume: &BTreeMap<String, ResumeEntry>,
    from: Option<&str>,
) -> Result<SampleOutcome> {
    let mut entries: Vec<(usize, serde_json::Value)> = Vec::new();
    let mut failed: Vec<String> = Vec::new();

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
        return Ok(SampleOutcome { entries, failed });
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
            eprintln!("fasthelper: {tool} not found, step skipped");
            continue;
        };
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
            "fastsim" => (
                PathBuf::new(), // simulate-only: the generated FASTQs are inputs
                None,
                work.join("fastsim.json"),
                work.join("maps"), // unused
            ),
            other => {
                eprintln!("fasthelper: unsupported step '{other}' — skipped");
                continue;
            }
        };
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
            &outdir.join(&cfg.simulate.name).join("sim"),
            config_dir,
        );
        let log = logs.join(format!("{}.{}.log", sample.name, tool));
        // Resume: skip this step when a previous manifest recorded the exact
        // same (argv, binary, version). --from <step> forces a rerun from
        // that step onward, so once the named step is reached everything is
        // forced regardless of what the manifest says.
        let forced = from.map(|f| f == tool).unwrap_or(false);
        if forced {
            from_reached = true;
        }
        let key = format!("{}/{}", sample.name, tool);
        if !forced && !from_reached {
            if let Some((prev_argv, prev_bin, prev_ver)) = resume.get(&key) {
                if prev_argv == &argv
                    && prev_bin == &s(&bin)
                    && *prev_ver == tool_version(&bin)
                {
                    println!("fasthelper: resumed {key}");
                    entries.push((
                        step_idx,
                        serde_json::json!({
                            "sample": sample.name,
                            "tool": tool,
                            "binary": prev_bin,
                            "argv": prev_argv,
                            "exit": 0,
                            "elapsed_s": 0.0,
                            "log": s(&log),
                            "resumed": true,
                        }),
                    ));
                    // The chain continues from this step's output: the file
                    // the resumed step produced (same argv) is still there.
                    if tool != "fastmaps" {
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
        let outcome = run_step(tool, &bin, &argv, &log, dry_run, None)?;
        println!(
            "  {tool:<10} exit {:<3} {:>7.2}s  log {}",
            outcome.exit,
            outcome.elapsed,
            outcome.log.display()
        );
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
            }),
        ));
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
        if tool != "fastmaps" {
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
    Ok(SampleOutcome { entries, failed })
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

/// What a previous run recorded for one `sample/tool` step: the argv it ran,
/// the binary path, and that binary's `--version` at the time. A step is only
/// resumed when all three match today.
type ResumeEntry = (Vec<String>, String, Option<String>);

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
            // Version is re-resolved live from the recorded binary; None when
            // that binary no longer answers --version (then never resume).
            let ver = if binary.is_empty() {
                None
            } else {
                tool_version(Path::new(&binary))
            };
            map.insert(format!("{sample}/{tool}"), (argv, binary, ver));
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
            a.push(threads.to_string());
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
            a.push(threads.to_string());
            a.push("--gz-level".into());
            a.push(gz.to_string());
        }
        "fastmaps" => {
            a.push("-o".into());
            a.push(s(maps_dir));
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
                a.extend(sim.args.iter().cloned());
            }
        }
        _ => {}
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
    let tools_dir = cfg.tools_dir(args.tools_dir.clone());
    let ts = toolset(tools_dir.clone());
    if ts.dedup.is_none() && ts.prep.is_none() && ts.maps.is_none() {
        bail!(
            "no FastSuite binaries found (looked in {} and $PATH) — build them with \
             `cargo build --release` in each tools/<name> dir, or pass --tools-dir",
            tools_dir
                .map(|d| d.display().to_string())
                .unwrap_or_else(|| "$PATH".into())
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
        let mut argv: Vec<String> = vec![
            "--pairs".into(),
            cfg.simulate.pairs.to_string(),
            "-o".into(),
            s(&prefix),
            "--json".into(),
            s(&simdir.join("fastsim.json")),
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
    let per_sample_threads = (cfg.run.threads / args.jobs.max(1)).max(1);
    let queue = std::sync::Mutex::new(std::collections::VecDeque::from(
        samples.iter().enumerate().collect::<Vec<_>>(),
    ));
    let entries: std::sync::Mutex<Vec<(usize, usize, serde_json::Value)>> =
        std::sync::Mutex::new(Vec::new());
    let failures: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());

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
                    per_sample_threads,
                    &resume,
                    args.from.as_deref(),
                ) {
                    Ok(mut o) => {
                        entries
                            .lock()
                            .unwrap()
                            .extend(o.entries.into_iter().map(|(step_idx, m)| (sample_idx, step_idx, m)));
                        failures.lock().unwrap().append(&mut o.failed);
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

    // Manifest + aggregate report over everything the run produced.
    let manifest_path = outdir.join("manifest.json");
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
            "commands": manifest,
            "failed": failed,
        }))? + "\n",
    )
    .with_context(|| format!("writing {}", manifest_path.display()))?;

    if args.dry_run {
        println!("fasthelper: dry-run — nothing executed");
        return Ok(());
    }

    let reports = report::collect(&[outdir.clone()], None)?;
    let formats: Vec<String> = vec!["json".into(), "tsv".into(), "md".into(), "html".into()];
    let title = format!("FastSuite report — {}", cfg.run.name);
    let written = report::write_reports(&title, &reports, &outdir, &formats)?;
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
