//! `fasthelper check -c config.toml` — validate a run config before it starts.
//!
//! The point is to fail at minute zero, not at hour three: every [[samples]]
//! input exists, every requested binary resolves, every step name is known,
//! and every flag passed under `[steps.<tool>] args` is accepted by that
//! tool's `--help` text.

use anyhow::{bail, Result};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::cli::CheckArgs;
use crate::config::Config;
use crate::runner;

/// Run `<bin> --help` and return the combined output text.
fn help_text(bin: &Path) -> Option<String> {
    let out = Command::new(bin).arg("--help").output().ok()?;
    let mut s = String::from_utf8_lossy(&out.stdout).to_string();
    s.push_str(&String::from_utf8_lossy(&out.stderr));
    Some(s)
}

/// Every long-flag token the help text advertises, as WHOLE tokens.
///
/// The old test was `help.contains(flag)`, so `--min-le` was accepted because
/// `--min-len` appears in the help. Usage lines wrap a flag in punctuation and
/// may attach a value (`--flag`, `--flag,`, `[--flag`, `--flag=VALUE`), so each
/// word is cut at `=` and trimmed of the punctuation around it.
fn help_flags(help: &str) -> BTreeSet<String> {
    help.split_whitespace()
        .map(|w| w.split('=').next().unwrap_or(w))
        .map(|w| w.trim_matches(|c: char| !(c.is_ascii_alphanumeric() || c == '-' || c == '_')))
        .filter(|w| w.starts_with("--") && w.len() > 2)
        .map(str::to_string)
        .collect()
}

/// The flag a config argument names: `--min-len=30` is `--min-len`.
fn flag_token(arg: &str) -> &str {
    arg.split('=').next().unwrap_or(arg)
}

/// A config value the TOML types accept but the tools do not.
///
/// `check` reports these next to everything else it finds, and `run` bails on
/// them at minute zero: `check -c` used to bless `gz_level = 42` (every tool
/// requires 1..=9), an absurd thread count, and a `steps = []` config that
/// produced an empty report and exited 0.
pub fn config_problems(cfg: &Config) -> Vec<String> {
    /// 0 is the documented "all cores" sentinel; an explicit count is a worker
    /// count. The ceiling is a typo gate, not a policy on oversubscription.
    const MAX_THREADS: usize = 4096;
    let mut out = Vec::new();
    if cfg.run.steps.is_empty() {
        out.push("run.steps is empty — there is nothing to run (list at least one step)".to_string());
    }
    if !(1..=9).contains(&cfg.run.gz_level) {
        out.push(format!(
            "run.gz_level is {} — every tool requires --gz-level in [1, 9]",
            cfg.run.gz_level
        ));
    }
    if cfg.run.threads > MAX_THREADS {
        out.push(format!(
            "run.threads is {} — 0 means all cores; an explicit count must be in [1, {MAX_THREADS}]",
            cfg.run.threads
        ));
    }
    out
}

pub fn cmd_check(args: &CheckArgs) -> Result<()> {
    if let Some(path) = args.report.as_ref() {
        return check_report(path);
    }
    let Some(config) = args.config.as_ref() else {
        // clap's `required_unless_present` catches this; the guard is for a
        // future signature change.
        bail!("pass -c <config.toml> or --report <report.json>");
    };
    let cfg = Config::load(config)?;
    let config_dir = config
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));
    let tools_dir = cfg.tools_dir(args.tools_dir.clone());
    let ts = runner::toolset(tools_dir);

    let mut errors: Vec<String> = config_problems(&cfg);

    // 1. Step names, then every flag the step passes.
    for step in &cfg.run.steps {
        if !runner::is_known_step(step) {
            errors.push(format!(
                "unknown step '{step}' — supported: {}",
                runner::KNOWN_STEPS.join(", ")
            ));
        }
    }
    for step in &cfg.run.steps {
        let s = cfg.step(step);
        if !s.enabled {
            continue;
        }
        let Some(bin) = ts.get(step) else { continue };
        let Some(h) = help_text(bin) else { continue };
        let accepted = help_flags(&h);
        // EVERY token, not only the first: `.find()` stopped at the first
        // `--flag`, so a typo in the second one was never looked at.
        for arg in s.args.iter().filter(|a| a.starts_with("--")) {
            let token = flag_token(arg);
            if !accepted.contains(token) {
                errors.push(format!(
                    "{step} does not accept '{token}' (see {})",
                    bin.display()
                ));
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
        if runner::is_known_step(step) && cfg.step(step).enabled && ts.get(step).is_none()
        {
            errors.push(format!("tool {step} not found in the tools dir or $PATH"));
        }
    }

    // 4. The simulated control. `run` bails at minute zero on both of these, so
    // `check` returning OK for them was a config the validator blessed and the
    // runner refused.
    if cfg.simulate.enabled {
        if ts.get("fastsim").is_none() {
            errors.push(
                "simulate.enabled is set but the fastsim binary was not found in the tools dir \
                 or $PATH"
                    .into(),
            );
        }
        if let Some(prof) = &cfg.simulate.profile {
            let p = if prof.is_absolute() {
                prof.clone()
            } else {
                config_dir.join(prof)
            };
            if !p.is_file() {
                errors.push(format!("[simulate] profile {} not found", p.display()));
            }
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

/// `fasthelper check --report <json>` — validate a tool's `fastsuite.report.v1`.
///
/// The suite has no shared crate (D-1): the report schema is the only contract
/// between the tools, and until this existed nothing could REJECT a report. The
/// rules live in [`validate_report`], which the runner calls too.
pub fn check_report(path: &Path) -> Result<()> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("reading {}: {e}", path.display()))?;
    let doc: serde_json::Value = serde_json::from_str(&text)
        .map_err(|e| anyhow::anyhow!("parsing {}: {e}", path.display()))?;

    // fasthelper's own `report.json` is a different document: a collection of
    // per-tool reports, stamped `fastsuite.aggregate.v1`. Checking it against
    // the v1 report rules rejected it on four keys it never had.
    if doc.get("schema").and_then(|v| v.as_str())
        == Some(crate::report::AGGREGATE_SCHEMA_VERSION)
    {
        return check_aggregate(path, &doc);
    }

    let verdict = validate_report(&doc);
    for w in &verdict.warnings {
        eprintln!("fasthelper check: warning: {w}");
    }
    if verdict.errors.is_empty() {
        println!(
            "fasthelper check: OK — {} report from {}, schema fastsuite.report.v1, {} metric(s), {} warning(s)",
            path.display(),
            doc.get("tool").and_then(|v| v.as_str()).unwrap_or("?"),
            verdict.metrics,
            verdict.warnings.len()
        );
        return Ok(());
    }
    for e in &verdict.errors {
        eprintln!("fasthelper check: {e}");
    }
    let n = verdict.errors.len();
    bail!("{n} problem(s) in the report")
}

/// The verdict of the `fastsuite.report.v1` rule set.
pub struct ReportCheck {
    /// Reasons the report cannot be trusted: the caller must refuse it.
    pub errors: Vec<String>,
    /// Contract advice (a missing recommended key); never fatal.
    pub warnings: Vec<String>,
    /// Metric count in `metrics` (0 when the object is missing).
    pub metrics: usize,
}

/// The ONE implementation of the `fastsuite.report.v1` rules.
///
/// `check --report` and the runner both call it, so a report the validator
/// rejects can never be ingested by the production path:
///
/// 1. the `schema` tag is `fastsuite.report.v1` and the identifying keys are
///    present (`tool`, `tool_version`, `command`);
/// 2. every `metrics` value is a number (rule 1: a flat numeric namespace);
/// 3. the recommended shared keys are reported — a WARNING, not an error,
///    because the schema says none are mandatory;
/// 4. the record identity. `reads_in == reads_out + reads_dropped` when those
///    three operands are present: a report whose `reads_out` exceeds its
///    `reads_in` is not reconcilable. A run that adopted a previous output
///    states it as `reads_preloaded`, and its identity is then
///    `reads_in + reads_preloaded == reads_out + reads_duplicate`;
/// 5. the base identity `bases_in == bases_out + bases_trimmed + bases_dropped`
///    once `bases_out` and `bases_dropped` are present, a missing operand
///    reading as 0 (REPORT-SCHEMA).
pub fn validate_report(doc: &serde_json::Value) -> ReportCheck {
    let mut errors: Vec<String> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();

    match doc.get("schema").and_then(|v| v.as_str()) {
        Some("fastsuite.report.v1") => {}
        Some(other) => errors.push(format!("schema is {other:?}, expected \"fastsuite.report.v1\"")),
        None => errors.push("no `schema` key".into()),
    }
    for key in ["tool", "tool_version", "command"] {
        if doc.get(key).and_then(|v| v.as_str()).is_none() {
            errors.push(format!("missing or non-string `{key}`"));
        }
    }
    let metrics = doc.get("metrics").and_then(|v| v.as_object());
    let empty = serde_json::Map::new();
    let metrics = match metrics {
        Some(m) => m,
        None => {
            errors.push("missing `metrics` object".into());
            &empty
        }
    };

    // Rule 1: a flat numeric namespace. A nested object here is how a tool
    // quietly breaks every aggregator downstream.
    for (k, v) in metrics {
        if !v.is_number() {
            errors.push(format!("metric `{k}` is not a number ({v})"));
        }
    }

    const RECOMMENDED: [&str; 8] = [
        "reads_in",
        "reads_out",
        "reads_dropped",
        "bases_in",
        "bases_out",
        "bases_trimmed",
        "bases_dropped",
        "reads_adapter_trimmed",
    ];
    for key in RECOMMENDED {
        if !metrics.contains_key(key) {
            warnings.push(format!("no `{key}` (recommended shared key, not mandatory)"));
        }
    }

    let num = |k: &str| metrics.get(k).and_then(|v| v.as_f64());
    // The record identity. A resumed run publishes `reads_preloaded` — the
    // reads it adopted instead of reprocessing — and its budget is the
    // three-term one; the two-term rule would reject a legitimate resume
    // (`reads_out` alone can exceed `reads_in`).
    if num("reads_preloaded").is_some() {
        if let (Some(rin), Some(pre), Some(rout), Some(rdup)) = (
            num("reads_in"),
            num("reads_preloaded"),
            num("reads_out"),
            num("reads_duplicate"),
        ) {
            if rin + pre != rout + rdup {
                errors.push(format!(
                    "record identity broken: reads_in {rin} + reads_preloaded {pre} != \
                     reads_out {rout} + reads_duplicate {rdup} (sum {})",
                    rout + rdup
                ));
            }
        }
    } else if let (Some(rin), Some(rout), Some(rdrop)) =
        (num("reads_in"), num("reads_out"), num("reads_dropped"))
    {
        if rin != rout + rdrop {
            errors.push(format!(
                "record identity broken: reads_in {rin} != reads_out {rout} + \
                 reads_dropped {rdrop} (sum {})",
                rout + rdrop
            ));
        }
    }

    // The base identity, whenever the report states what left and what was
    // dropped: `bases_trimmed` (and `bases_in`) may be absent and then read as
    // 0, which is what REPORT-SCHEMA documents. Requiring all four operands
    // skipped the check for every report that named only three.
    if let (Some(bout), Some(bdrop)) = (num("bases_out"), num("bases_dropped")) {
        let bin = num("bases_in").unwrap_or(0.0);
        let btrim = num("bases_trimmed").unwrap_or(0.0);
        if bin != bout + btrim + bdrop {
            errors.push(format!(
                "base identity broken: bases_in {bin} != bases_out {bout} + bases_trimmed \
                 {btrim} + bases_dropped {bdrop} (sum {})",
                bout + btrim + bdrop
            ));
        }
    }

    ReportCheck {
        errors,
        warnings,
        metrics: metrics.len(),
    }
}

/// `fasthelper check --report` for the aggregate document
/// (`fastsuite.aggregate.v1`): `generator`, `generator_version`, `title`, an
/// optional numeric `threads`, and a `reports` array in which every row names
/// its sample, its tool and a numeric metric map. That is the shape
/// `concatenate`/`run` write, and the shape a consumer can actually read.
fn check_aggregate(path: &Path, doc: &serde_json::Value) -> Result<()> {
    let mut errors: Vec<String> = Vec::new();

    for key in ["generator", "generator_version", "title"] {
        if doc.get(key).and_then(|v| v.as_str()).is_none() {
            errors.push(format!("missing or non-string `{key}`"));
        }
    }
    if doc.get("threads").is_some_and(|t| !t.is_number()) {
        errors.push("`threads` is not a number".into());
    }

    let rows = doc.get("reports").and_then(|v| v.as_array());
    let rows = match rows {
        Some(r) => r,
        None => {
            errors.push("missing `reports` array".into());
            for e in &errors {
                eprintln!("fasthelper check: {e}");
            }
            bail!("{} problem(s) in the aggregate", errors.len());
        }
    };

    for (i, row) in rows.iter().enumerate() {
        for key in ["sample", "tool"] {
            if row.get(key).and_then(|v| v.as_str()).is_none() {
                errors.push(format!("reports[{i}]: missing or non-string `{key}`"));
            }
        }
        match row.get("metrics").and_then(|m| m.as_object()) {
            Some(m) => {
                for (k, v) in m {
                    if !v.is_number() {
                        errors.push(format!("reports[{i}]: metric `{k}` is not a number ({v})"));
                    }
                }
            }
            None => errors.push(format!("reports[{i}]: missing `metrics` object")),
        }
    }

    if errors.is_empty() {
        println!(
            "fasthelper check: OK — {} aggregate from {}, schema {}, {} report(s)",
            path.display(),
            doc.get("generator").and_then(|v| v.as_str()).unwrap_or("?"),
            crate::report::AGGREGATE_SCHEMA_VERSION,
            rows.len()
        );
        return Ok(());
    }
    for e in &errors {
        eprintln!("fasthelper check: {e}");
    }
    let n = errors.len();
    bail!("{n} problem(s) in the aggregate")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str, body: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!(
            "fasthelper-check-{}-{}-{name}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    const GOOD: &str = r#"{
      "schema": "fastsuite.report.v1",
      "tool": "fastprep",
      "tool_version": "0.4.1",
      "command": "fastprep -1 a.fq -o b.fq",
      "metrics": {
        "reads_in": 100, "reads_out": 90, "reads_dropped": 10,
        "reads_adapter_trimmed": 0,
        "bases_in": 15000, "bases_out": 13230,
        "bases_trimmed": 1200, "bases_dropped": 570
      }
    }"#;

    #[test]
    fn a_reconciling_report_passes() {
        let p = tmp("good.json", GOOD);
        assert!(check_report(&p).is_ok());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn a_wrong_schema_tag_is_refused() {
        let p = tmp("schema.json", &GOOD.replace("fastsuite.report.v1", "fastsuite.report.v0"));
        assert!(check_report(&p).is_err(), "a v0 tag must be refused");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn a_broken_base_identity_is_refused() {
        // The pre-fix fastprep report: trims on dropped reads were not counted,
        // so bases_in did not close.
        let p = tmp("identity.json", &GOOD.replace("\"bases_trimmed\": 1200", "\"bases_trimmed\": 0"));
        let err = check_report(&p).unwrap_err().to_string();
        assert!(err.contains("problem"), "{err}");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn a_report_where_reads_out_exceeds_reads_in_is_refused() {
        // The reproduced defect: reads_in 10, reads_out 1e6 passed as OK while
        // REPORT-SCHEMA states reads_in == reads_out + reads_dropped. Nothing
        // checked the record identity at all, so a report could invent reads.
        let p = tmp(
            "reads_identity.json",
            &GOOD.replace(
                "\"reads_out\": 90, \"reads_dropped\": 10",
                "\"reads_out\": 1000000, \"reads_dropped\": 0",
            ),
        );
        let err = check_report(&p).unwrap_err().to_string();
        assert!(err.contains("problem"), "{err}");
        std::fs::remove_file(p).ok();

        // The rule itself names the broken identity, not just "a problem".
        let doc: serde_json::Value = serde_json::from_str(
            &GOOD.replace(
                "\"reads_out\": 90, \"reads_dropped\": 10",
                "\"reads_out\": 1000000, \"reads_dropped\": 0",
            ),
        )
        .unwrap();
        let v = validate_report(&doc);
        assert!(
            v.errors.iter().any(|e| e.contains("record identity broken")),
            "the record identity must be the reason: {:?}",
            v.errors
        );
    }

    #[test]
    fn a_resumed_run_reconciles_through_reads_preloaded() {
        // fastdedup publishes `reads_preloaded` when it adopts an existing
        // output: reads_out alone can legitimately exceed reads_in, so the
        // two-term rule must not be applied to that shape.
        let resumed = r#"{
          "schema": "fastsuite.report.v1",
          "tool": "fastdedup",
          "tool_version": "1.3.0",
          "command": "fastdedup -1 a.fq -o b.fq",
          "metrics": {
            "reads_in": 4, "reads_preloaded": 19999,
            "reads_out": 20000, "reads_duplicate": 3
          }
        }"#;
        let doc: serde_json::Value = serde_json::from_str(resumed).unwrap();
        let v = validate_report(&doc);
        assert!(v.errors.is_empty(), "{:?}", v.errors);

        // A resumed report whose three terms do not close is still refused.
        let broken: serde_json::Value =
            serde_json::from_str(&resumed.replace("\"reads_duplicate\": 3", "\"reads_duplicate\": 4"))
                .unwrap();
        let v = validate_report(&broken);
        assert!(
            v.errors.iter().any(|e| e.contains("reads_preloaded")),
            "the three-term identity must be reported: {:?}",
            v.errors
        );
    }

    #[test]
    fn the_base_identity_is_checked_without_bases_trimmed() {
        // REPORT-SCHEMA: the identity applies when bases_out and bases_dropped
        // are present, a missing operand reading as 0. Requiring all four meant
        // a report that simply omitted bases_trimmed skipped the check.
        let p = tmp("no_trim.json", &GOOD.replace("\"bases_trimmed\": 1200, ", ""));
        let err = check_report(&p).unwrap_err().to_string();
        assert!(err.contains("problem"), "{err}");
        std::fs::remove_file(p).ok();

        // The same shape that DOES close (bases_in == bases_out + bases_dropped)
        // must pass: the check is the identity, not the presence of a key.
        let closing = GOOD
            .replace("\"bases_trimmed\": 1200, ", "")
            .replace("\"bases_in\": 15000", "\"bases_in\": 13800");
        let p = tmp("no_trim_ok.json", &closing);
        assert!(check_report(&p).is_ok(), "a closing three-term identity must pass");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn a_nested_metric_is_refused() {
        let p = tmp(
            "nested.json",
            &GOOD.replace("\"bases_in\": 15000", "\"bases_in\": {\"r1\": 15000}"),
        );
        assert!(check_report(&p).is_err(), "a nested metric must be refused");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn a_broken_base_identity_written_as_floats_is_refused_too() {
        // H4: every operand was read with as_u64, so a report that wrote its
        // bases as JSON floats skipped the identity check and passed — while
        // the integer spelling of the same report was refused.
        let broken = GOOD
            .replace("\"bases_in\": 15000", "\"bases_in\": 15000.0")
            .replace("\"bases_out\": 13230", "\"bases_out\": 1.0")
            .replace("\"bases_trimmed\": 1200", "\"bases_trimmed\": 2.0")
            .replace("\"bases_dropped\": 570", "\"bases_dropped\": 3.0");
        let p = tmp("float_identity.json", &broken);
        let err = check_report(&p).unwrap_err().to_string();
        assert!(err.contains("problem"), "{err}");
        std::fs::remove_file(p).ok();

        // The float spelling that DOES close must still pass.
        let good = GOOD
            .replace("\"bases_in\": 15000", "\"bases_in\": 15000.0")
            .replace("\"bases_out\": 13230", "\"bases_out\": 13230.0")
            .replace("\"bases_trimmed\": 1200", "\"bases_trimmed\": 1200.0")
            .replace("\"bases_dropped\": 570", "\"bases_dropped\": 570.0");
        let p = tmp("float_identity_ok.json", &good);
        assert!(check_report(&p).is_ok(), "a closing float identity must pass");
        std::fs::remove_file(p).ok();
    }
}
