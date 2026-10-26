//! `fasthelper install --prefix DIR` — one bin/ with every suite binary.
//!
//! Legacy binary names (rpolyx, fdedup, novaseq-sim) are installed as
//! symlinks alongside the canonical ones, per docs/DECISIONS.md D-8: cargo and
//! bioconda never install a binary name, fasthelper does.

use anyhow::{bail, Context, Result};
use std::os::unix::fs::symlink;
use std::path::{Path, PathBuf};

use crate::cli::InstallArgs;
use crate::runner::{self, ToolSet};

/// canonical name -> legacy symlink name (D-8)
const LEGACY: &[(&str, &str)] = &[
    ("fastprep", "rpolyx"),
    ("fastdedup", "fdedup"),
    ("fastsim", "novaseq-sim"),
];

fn link(dst: &Path, target: &Path) -> Result<()> {
    if let Ok(md) = std::fs::symlink_metadata(dst) {
        if md.file_type().is_symlink() {
            std::fs::remove_file(dst).ok();
        }
    }
    symlink(target, dst)
        .with_context(|| format!("linking {} -> {}", dst.display(), target.display()))
}

pub fn cmd_install(args: &InstallArgs) -> Result<()> {
    if cfg!(not(unix)) {
        bail!("install is only implemented for unix platforms");
    }
    let ts: ToolSet = runner::toolset(args.tools_dir.clone());
    let mut targets: Vec<(String, PathBuf)> = Vec::new();
    // The suite a prefix is supposed to hold. `fasthelper` is this binary, so it
    // always resolves; the other five may not, and leaving them out silently is
    // what made an install of 1 of 6 look like success.
    let mut missing: Vec<&str> = Vec::new();
    for tool in [
        "fastcheck",
        "fastdedup",
        "fastprep",
        "fastsim",
        "fastmaps",
        "fasthelper",
    ] {
        let bin = if tool == "fasthelper" {
            std::env::current_exe().ok()
        } else {
            ts.get(tool).cloned()
        };
        match bin {
            Some(b) => targets.push((tool.to_string(), b)),
            None => missing.push(tool),
        }
    }
    if targets.is_empty() {
        bail!("no FastSuite binaries found — build them first or pass --tools-dir");
    }
    let bin = args.prefix.join("bin");
    std::fs::create_dir_all(&bin).with_context(|| format!("creating {}", bin.display()))?;
    let mut made = Vec::new();
    for (name, target) in &targets {
        link(&bin.join(name), target)?;
        made.push(format!("{name} -> {}", target.display()));
    }
    if !args.no_legacy {
        for (canonical, legacy) in LEGACY {
            if let Some((_, target)) = targets.iter().find(|(n, _)| n == canonical) {
                link(&bin.join(legacy), target)?;
                made.push(format!("{legacy} -> {}", target.display()));
            }
        }
    }
    println!("fasthelper install: {} symlink(s) in {}", made.len(), bin.display());
    for m in &made {
        println!("  {m}");
    }
    if !missing.is_empty() {
        bail!(
            "{} tool(s) not found, so they were not linked: {} — build them (cargo build \
             --release in tools/<tool>) or pass --tools-dir",
            missing.len(),
            missing.join(", ")
        );
    }
    Ok(())
}
