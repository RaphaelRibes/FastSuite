//! `fasthelper check -c config.toml` — validate a run config before it starts.
//!
//! The point is to fail at minute zero, not at hour three: every [[samples]]
//! input exists, every requested binary resolves, every step name is known,
//! and every flag passed under `[steps.<tool>] args` is accepted by that
//! tool's `--help` text.

use anyhow::{bail, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::cli::CheckArgs;
use crate::config::Config;
use crate::runner;

const KNOWN_STEPS: &[&str] = &["fastdedup", "fastprep", "fastmaps", "fastcheck"];

/// Run `<bin> --help` and return the combined output text.
fn help_text(bin: &Path) -> Option<String> {
    let out = Command::new(bin).arg("--help").output().ok()?;
    let mut s = String::from_utf8_lossy(&out.stdout).to_string();
    s.push_str(&String::from_utf8_lossy(&out.stderr));
    Some(s)
}

/// Is `flag` (a token starting with `--`) accepted by the tool's help text?
fn flag_accepted(help: &str, flag: &str) -> bool {
    // --insert-hist and --insert-hist=... both count; --foo bar never does.
    help.contains(flag)
}

pub fn cmd_check(args: &CheckArgs) -> Result<()> {
    let cfg = Config::load(&args.config)?;
    let config_dir = args
        .config
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));
    let tools_dir = cfg.tools_dir(args.tools_dir.clone());
    let ts = runner::toolset(tools_dir);

    let mut errors: Vec<String> = Vec::new();

    // 1. Step names.
    for step in &cfg.run.steps {
        if !KNOWN_STEPS.contains(&step.as_str()) {
            errors.push(format!(
                "unknown step '{step}' — supported: {}",
                KNOWN_STEPS.join(", ")
            ));
        }
        if let Some(extra) = cfg.step(step).args.iter().find(|a| a.starts_with("--")) {
            if cfg.step(step).enabled {
                if let Some(bin) = ts.get(step) {
                    if let Some(h) = help_text(bin) {
                        if !flag_accepted(&h, extra) {
                            errors.push(format!(
                                "{step} does not accept '{extra}' (see {})",
                                bin.display()
                            ));
                        }
                    }
                }
            }
        }
    }

    // 2. Sample inputs.
    for s in &cfg.samples {
        for f in [Some(&s.r1), s.r2.as_ref()].into_iter().flatten() {
            let p = if f.is_absolute() {
                f.clone()
            } else {
                config_dir.join(f)
            };
            if !p.is_file() {
                errors.push(format!("sample {}: input {} not found", s.name, p.display()));
            }
        }
    }

    // 3. Tools.
    for step in &cfg.run.steps {
        if KNOWN_STEPS.contains(&step.as_str()) && cfg.step(step).enabled && ts.get(step).is_none()
        {
            errors.push(format!("tool {step} not found in the tools dir or $PATH"));
        }
    }

    if errors.is_empty() {
        println!(
            "fasthelper check: OK — {} sample(s), {} step(s), tools resolved",
            cfg.samples.len(),
            cfg.run.steps.len()
        );
        return Ok(());
    }
    for e in &errors {
        eprintln!("fasthelper check: {e}");
    }
    let n = errors.len();
    bail!("{n} problem(s) — fix them before running")
}
