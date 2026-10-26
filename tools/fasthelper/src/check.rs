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

    let mut errors: Vec<String> = Vec::new();

    // 1. Step names.
    for step in &cfg.run.steps {
        if !runner::is_known_step(step) {
            errors.push(format!(
                "unknown step '{step}' — supported: {}",
                runner::KNOWN_STEPS.join(", ")
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
/// rules checked here are the ones the schema states:
///
/// 1. the `schema` tag is `fastsuite.report.v1` and the identifying keys are
///    present (`tool`, `tool_version`, `command`, `metrics`);
/// 2. every `metrics` value is a number (rule 1: flat numeric namespace);
/// 3. the recommended shared keys are reported — a WARNING, not an error,
///    because the schema says none are mandatory;
/// 4. when `bases_out` and `bases_dropped` are both present, the base identity
///    `bases_in == bases_out + bases_trimmed + bases_dropped` must hold. That
///    identity is what makes a run reconcilable, and no tool could be checked
///    against it before.
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
    let Some(metrics) = doc.get("metrics").and_then(|v| v.as_object()) else {
        errors.push("missing `metrics` object".into());
        for e in &errors {
            eprintln!("fasthelper check: {e}");
        }
        bail!("{} problem(s) in the report", errors.len());
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
    if let (Some(bin), Some(bout), Some(btrim), Some(bdrop)) = (
        num("bases_in"),
        num("bases_out"),
        num("bases_trimmed"),
        num("bases_dropped"),
    ) {
        if bin != bout + btrim + bdrop {
            errors.push(format!(
                "base identity broken: bases_in {bin} != bases_out {bout} + bases_trimmed \
                 {btrim} + bases_dropped {bdrop} (sum {})",
                bout + btrim + bdrop
            ));
        }
    }

    for w in &warnings {
        eprintln!("fasthelper check: warning: {w}");
    }
    if errors.is_empty() {
        println!(
            "fasthelper check: OK — {} report from {}, schema fastsuite.report.v1, {} metric(s), {} warning(s)",
            path.display(),
            doc.get("tool").and_then(|v| v.as_str()).unwrap_or("?"),
            metrics.len(),
            warnings.len()
        );
        return Ok(());
    }
    for e in &errors {
        eprintln!("fasthelper check: {e}");
    }
    let n = errors.len();
    bail!("{n} problem(s) in the report")
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
