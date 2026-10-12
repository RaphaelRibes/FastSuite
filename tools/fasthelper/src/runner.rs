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

    let mut manifest = Vec::new();
    let mut failed = Vec::new();

    for sample in &samples {
        let work = outdir.join(&sample.name);
        std::fs::create_dir_all(&work)?;
        let in1 = resolve_rel(&sample.r1, &config_dir);
        let in2 = sample.r2.as_ref().map(|p| resolve_rel(p, &config_dir));
        if !args.dry_run && !in1.is_file() {
            eprintln!("fasthelper: skipping {}: input {} not found", sample.name, in1.display());
            failed.push(sample.name.clone());
            continue;
        }

        println!("fasthelper: sample {} ({} steps)", sample.name, cfg.run.steps.len());
        let mut cur1 = in1.clone();
        let mut cur2 = in2.clone();
        let mut to_delete: Vec<PathBuf> = Vec::new();

        for tool in &cfg.run.steps {
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
                cfg.run.threads,
                cfg.run.gz_level,
                &step.args,
                &maps_dir,
            );
            let log = logs.join(format!("{}.{}.log", sample.name, tool));
            let outcome = run_step(tool, &bin, &argv, &log, args.dry_run, None)?;
            println!(
                "  {tool:<10} exit {:<3} {:>7.2}s  log {}",
                outcome.exit,
                outcome.elapsed,
                outcome.log.display()
            );
            manifest.push(serde_json::json!({
                "sample": sample.name,
                "tool": tool,
                "binary": s(&bin),
                "argv": outcome.argv,
                "exit": outcome.exit,
                "elapsed_s": outcome.elapsed,
                "log": s(&outcome.log),
            }));
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
    }

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
