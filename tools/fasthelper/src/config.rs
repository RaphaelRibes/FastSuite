//! Run-config model: `fastsuite.toml`.
//!
//! One file describes a whole preprocessing run: which samples, which steps,
//! and the shared runtime knobs every FastSuite tool understands
//! (`--threads`, `--gz-level`). Per-step escape hatches live under
//! `[steps.<tool>] args = [...]` so a tool-specific flag never blocks the
//! common path.

use anyhow::{Context, Result};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Deserialize, Default)]
pub struct Config {
    #[serde(default)]
    pub run: Run,
    #[serde(default)]
    pub samples: Vec<Sample>,
    #[serde(default)]
    pub steps: BTreeMap<String, Step>,
    #[serde(default)]
    pub simulate: Simulate,
}

#[derive(Debug, Deserialize)]
pub struct Run {
    /// Run label, used in the report title and the manifest.
    #[serde(default = "default_name")]
    pub name: String,
    /// Where every output of the run lands.
    #[serde(default = "default_outdir")]
    pub outdir: PathBuf,
    /// Threads handed to each tool via its own --threads/-j flag.
    #[serde(default = "default_threads")]
    pub threads: usize,
    /// gzip level handed to each tool via --gz-level.
    #[serde(default = "default_gz")]
    pub gz_level: u8,
    /// Steps to execute, in order. Supported: fastdedup, fastprep, fastmaps, fastcheck.
    #[serde(default = "default_steps")]
    pub steps: Vec<String>,
    /// Directory holding the tool checkouts (…/fastdedup, …/fastprep, …).
    #[serde(default)]
    pub tools_dir: Option<PathBuf>,
    /// Keep intermediate FASTQ instead of deleting it once consumed.
    #[serde(default)]
    pub keep_intermediate: bool,
}

fn default_name() -> String {
    "fastsuite".to_string()
}
fn default_outdir() -> PathBuf {
    PathBuf::from("fastsuite_run")
}
fn default_threads() -> usize {
    4
}
fn default_gz() -> u8 {
    6
}
fn default_steps() -> Vec<String> {
    vec![
        "fastdedup".to_string(),
        "fastprep".to_string(),
        "fastmaps".to_string(),
    ]
}

impl Default for Run {
    fn default() -> Self {
        Self {
            name: default_name(),
            outdir: default_outdir(),
            threads: default_threads(),
            gz_level: default_gz(),
            steps: default_steps(),
            tools_dir: None,
            keep_intermediate: false,
        }
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct Sample {
    pub name: String,
    pub r1: PathBuf,
    #[serde(default)]
    pub r2: Option<PathBuf>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct Step {
    #[serde(default = "yes")]
    pub enabled: bool,
    /// Extra arguments appended verbatim to the tool's command line.
    #[serde(default)]
    pub args: Vec<String>,
}

fn yes() -> bool {
    true
}

impl Default for Step {
    fn default() -> Self {
        Self {
            enabled: true,
            args: Vec::new(),
        }
    }
}

/// Optional synthetic-control step: generate reads with the simulator and push
/// them through the same chain, so every run carries a known-truth sample.
#[derive(Debug, Deserialize)]
pub struct Simulate {
    #[serde(default)]
    pub enabled: bool,
    #[serde(default = "default_sim_name")]
    pub name: String,
    #[serde(default = "default_pairs")]
    pub pairs: usize,
    #[serde(default)]
    pub profile: Option<PathBuf>,
    #[serde(default)]
    pub seed: Option<u64>,
    #[serde(default)]
    pub args: Vec<String>,
}

fn default_sim_name() -> String {
    "sim".to_string()
}
fn default_pairs() -> usize {
    20_000
}

impl Default for Simulate {
    fn default() -> Self {
        Self {
            enabled: false,
            name: default_sim_name(),
            pairs: default_pairs(),
            profile: None,
            seed: None,
            args: Vec::new(),
        }
    }
}

impl Config {
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading config {}", path.display()))?;
        let cfg: Config =
            toml::from_str(&text).with_context(|| format!("parsing TOML {}", path.display()))?;
        if cfg.samples.is_empty() && !cfg.simulate.enabled {
            anyhow::bail!(
                "config {} has no [[samples]] and simulate.enabled is false — nothing to do",
                path.display()
            );
        }
        // The work dir and the resume key derive from the sample name alone
        // (runner.rs), so two entries with one name would have two workers
        // create, unlink and write the same paths. Refuse it at minute zero.
        // The simulated control is appended to this list under simulate.name,
        // so its name has to be in the same set: a [[samples]] entry called
        // sim_E put two workers on one directory and the control's report
        // overwrote the sample's.
        let mut seen = std::collections::BTreeSet::new();
        for s in &cfg.samples {
            if !seen.insert(s.name.as_str()) {
                anyhow::bail!(
                    "config {}: duplicate sample name '{}' — every [[samples]] name must be \
                     unique, it names the sample's output directory",
                    path.display(),
                    s.name
                );
            }
        }
        if cfg.simulate.enabled && !seen.insert(cfg.simulate.name.as_str()) {
            anyhow::bail!(
                "config {}: duplicate sample name '{}' — [simulate] name collides with a \
                 [[samples]] entry, and it names the simulated control's output directory",
                path.display(),
                cfg.simulate.name
            );
        }
        Ok(cfg)
    }

    pub fn step(&self, tool: &str) -> Step {
        self.steps.get(tool).cloned().unwrap_or_default()
    }

    /// Resolve the tools directory: explicit flag/config, then $FASTSUITE_TOOLS_DIR,
    /// then the parent of this binary's workspace checkout.
    pub fn tools_dir(&self, cli_override: Option<PathBuf>) -> Option<PathBuf> {
        if let Some(d) = cli_override {
            return Some(d);
        }
        if let Some(d) = &self.run.tools_dir {
            return Some(d.clone());
        }
        if let Ok(d) = std::env::var("FASTSUITE_TOOLS_DIR") {
            if !d.is_empty() {
                return Some(PathBuf::from(d));
            }
        }
        // tools/fasthelper/target/release/fasthelper -> tools/
        let exe = std::env::current_exe().ok()?;
        let tools = exe
            .parent()? // release
            .parent()? // target
            .parent()? // fasthelper
            .to_path_buf();
        if tools.is_dir() {
            Some(tools)
        } else {
            None
        }
    }
}
