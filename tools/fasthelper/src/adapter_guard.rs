//! Post-trimming residual-adapter guard.
//!
//! fastcheck already measures, on any FASTQ it is given, how many reads still
//! carry adapter (`--adapter-content`, PLAN-QC module 11). Nothing measured that
//! on the *output* of a trimming run — exactly the failure mode a trimmer can
//! hide: under-trimmed reads carry adapter into every downstream
//! assembler/mapper, and the run still exits 0.
//!
//! This module runs that measurement twice around a trimming step — once on the
//! reads entering it, once on the reads leaving it — and warns when adapter
//! survives, naming the likely cause so the user can act.
//!
//! Two deliberate choices:
//!
//! * It measures with fastcheck's OWN permissive defaults (illumina catalogue,
//!   `--adapter-error-rate 0.2`, `--min-adapter-overlap 14`), NOT the trimmer's
//!   settings. A trimmer configured too strictly (error rate 0, overlap floor
//!   30) is precisely the case this must catch; re-using its settings would make
//!   the guard blind in the same way the trimmer is. The one exception is a
//!   custom `--adapters FILE` on the trimming step, which is passed through so
//!   the catalogue matches — an adapter missing from the catalogue would be
//!   invisible to both.
//!
//! * It is silent when the input carries no detectable adapter at all, so a run
//!   on adapter-free data (amplicons, simulations) never nags. And it is only
//!   called for a trimming step, so a run with none pays nothing (see runner.rs).

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Residual fraction of trimmed reads that trips the guard.
///
/// 5% is conservative in the sense that matters here: a working trim leaves
/// reads with essentially no adapter (measured fractions are 0.0–0.5%), so a
/// real failure sits well above this and a passing run sits well below it.
pub const DEFAULT_THRESHOLD: f64 = 0.05;

/// Required relative fall from the raw fraction: the trimmed fraction must be at
/// most `(1 - min_reduction) * raw`. 0.5 = the trim must at least halve it.
pub const DEFAULT_MIN_REDUCTION: f64 = 0.5;

/// One adapter measurement, from a fastcheck report's `adapter_content` block.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Measurement {
    /// Total reads measured (`metrics.reads`).
    pub reads: u64,
    /// Reads carrying at least one adapter match: the largest per-adapter count,
    /// so a read matching two adapters is not counted twice.
    pub reads_with_adapter: u64,
    /// Total match observations across the catalogue (a read may match twice).
    pub total_matches: u64,
}

impl Measurement {
    /// Fraction of reads carrying adapter, in [0, 1]. Zero reads => 0.
    pub fn fraction(&self) -> f64 {
        if self.reads == 0 {
            0.0
        } else {
            self.reads_with_adapter as f64 / self.reads as f64
        }
    }

    /// Parse a `fastsuite.report.v1` fastcheck document.
    ///
    /// `None` means "not measured" — the `adapter_content` block is absent (the
    /// module did not run) or the read count is zero. That is deliberately
    /// different from a measured zero (D-6): a caller must not read "absent" as
    /// "clean".
    pub fn from_report(doc: &serde_json::Value) -> Option<Measurement> {
        let blk = doc.get("adapter_content")?.as_object()?;
        let reads = doc.get("metrics")?.get("reads")?.as_u64()?;
        if reads == 0 {
            return None;
        }
        let adapters = blk.get("adapters")?.as_array()?;
        let mut reads_with_adapter = 0u64;
        for a in adapters {
            let n = a
                .get("reads_with_match")
                .and_then(|v| v.as_u64())
                .unwrap_or(0);
            reads_with_adapter = reads_with_adapter.max(n);
        }
        let total_matches = blk
            .get("total_matches")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        Some(Measurement {
            reads,
            reads_with_adapter,
            total_matches,
        })
    }
}

/// The guard's knobs, built from the `run` CLI flags.
#[derive(Debug, Clone)]
pub struct Policy {
    pub enabled: bool,
    /// Residual fraction above which the guard trips.
    pub threshold: f64,
    /// Relative fall required from the raw fraction (see [`DEFAULT_MIN_REDUCTION`]).
    pub min_reduction: f64,
    /// Turn a trip into a hard failure (exit 2, the `--fail-on` convention)
    /// instead of a warning.
    pub fail: bool,
}

impl Default for Policy {
    fn default() -> Self {
        Policy {
            enabled: true,
            threshold: DEFAULT_THRESHOLD,
            min_reduction: DEFAULT_MIN_REDUCTION,
            fail: false,
        }
    }
}

/// What the guard has to say about one trimming step.
#[derive(Debug, Clone, PartialEq)]
pub enum Verdict {
    /// Nothing to report (no adapter in the input, or the trim cleared it).
    Silent(&'static str),
    /// Adapter survived; report it and carry on.
    Warn(String),
    /// Adapter survived and the run asked for a hard failure.
    Fail(String),
}

impl Verdict {
    /// True for both `Warn` and `Fail`.
    pub fn is_trip(&self) -> bool {
        matches!(self, Verdict::Warn(_) | Verdict::Fail(_))
    }

    /// The actionable message, when there is one.
    pub fn message(&self) -> Option<&str> {
        match self {
            Verdict::Warn(m) | Verdict::Fail(m) => Some(m.as_str()),
            Verdict::Silent(_) => None,
        }
    }

    /// Why the guard stayed silent (for the manifest).
    pub fn silent_reason(&self) -> Option<&'static str> {
        match self {
            Verdict::Silent(r) => Some(r),
            _ => None,
        }
    }
}

/// The likely reason adapter survived, inferred from how the trimming step was
/// configured. The three candidates the task names, in the order they are most
/// often the cause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cause {
    /// `fastprep --no-adapter`: trimming was switched off.
    NoAdapterSet,
    /// No adapter set was given and auto-detection found none.
    AutoDetect,
    /// A set was configured but the match failed — error rate too strict or the
    /// overlap floor too long.
    StrictMatching,
}

impl Cause {
    /// Infer the cause from the trimming step's extra args (see `[[steps.fastprep]]`).
    pub fn from_step_args(args: &[String]) -> Cause {
        let has = |flag: &str| {
            args.iter()
                .any(|a| a == flag || a.starts_with(&format!("{flag}=")))
        };
        if has("--no-adapter") {
            Cause::NoAdapterSet
        } else if has("--adapter-preset") || has("--adapters") {
            Cause::StrictMatching
        } else {
            Cause::AutoDetect
        }
    }

    /// One-line headline naming the cause.
    pub fn headline(self) -> &'static str {
        match self {
            Cause::NoAdapterSet => {
                "adapter trimming is disabled for the step (fastprep --no-adapter)"
            }
            Cause::AutoDetect => {
                "the step had no adapter set configured and auto-detection found none in the sampled reads"
            }
            Cause::StrictMatching => {
                "an adapter set is configured but the match failed: the error rate may be too strict or the minimum overlap too long"
            }
        }
    }

    /// The FASTA catalogue to measure with: the step's own custom set when it has
    /// one, otherwise fastcheck's default illumina preset.
    pub fn catalogue_arg(self, args: &[String]) -> Option<String> {
        if self != Cause::StrictMatching {
            return None;
        }
        // `--adapters FILE` or `--adapters=FILE`.
        let mut it = args.iter();
        while let Some(a) = it.next() {
            if a == "--adapters" {
                return it.next().cloned();
            }
            if let Some(rest) = a.strip_prefix("--adapters=") {
                return Some(rest.to_string());
            }
        }
        None
    }
}

/// The remedies every trip message carries, so a user can act without reading
/// the docs.
pub const REMEDIES: &str = "try: fastprep --adapter-preset illumina|nextera (or --adapters FILE) \
if no adapter set was configured; --adapter-error-rate 0.3 (default 0.2) if matching is too \
strict; --min-adapter-overlap 8 (default 14) if the overlap floor is too long";

/// Compare the adapter measurement on the reads entering and leaving a trimming
/// step and decide. Pure: no I/O, so the whole decision is unit-testable from
/// inline reports.
pub fn evaluate(
    raw: &Measurement,
    trimmed: &Measurement,
    policy: &Policy,
    step_args: &[String],
) -> Verdict {
    // Silent when the input has no detectable adapter at all.
    if raw.reads_with_adapter == 0 {
        return Verdict::Silent("the trimming step's input carries no detectable adapter");
    }
    let raw_frac = raw.fraction();
    let trim_frac = trimmed.fraction();

    let still_present = trim_frac > policy.threshold;
    let not_reduced = trim_frac > (1.0 - policy.min_reduction) * raw_frac;
    if !(still_present || not_reduced) {
        return Verdict::Silent("adapter content fell below the residual threshold");
    }

    let cause = Cause::from_step_args(step_args);
    let reason = if still_present {
        format!(
            "residual adapter above the {:.1}% threshold",
            policy.threshold * 100.0
        )
    } else {
        format!(
            "adapter content did not fall from the raw {:.1}% (it must at least halve)",
            raw_frac * 100.0
        )
    };
    let msg = format!(
        "adapter survived the trim: {:.1}% of trimmed reads still carry adapter \
(raw {:.1}%, {} match(es) left) — {reason}\n  \
likely cause: {}\n  {REMEDIES}",
        trim_frac * 100.0,
        raw_frac * 100.0,
        trimmed.total_matches,
        cause.headline(),
    );
    if policy.fail {
        Verdict::Fail(msg)
    } else {
        Verdict::Warn(msg)
    }
}

/// Run `fastcheck <files...> --adapter-content --json <json>` and parse the
/// measurement out of the report. Returns `Ok(None)` when the tool fails or the
/// report carries no measurement (never a fabricated zero).
pub fn measure(
    check_bin: &Path,
    files: &[PathBuf],
    json: &Path,
    log: &Path,
    catalogue: Option<&str>,
    threads: usize,
) -> Result<Option<Measurement>> {
    if let Some(dir) = log.parent() {
        std::fs::create_dir_all(dir).ok();
    }
    let out = std::fs::File::create(log)
        .with_context(|| format!("creating adapter-guard log {}", log.display()))?;
    let err = out.try_clone()?;

    let mut cmd = Command::new(check_bin);
    cmd.args(files);
    cmd.arg("--adapter-content");
    if let Some(cat) = catalogue {
        cmd.arg("--adapter-preset").arg(cat);
    }
    cmd.arg("--json").arg(json);
    cmd.arg("-j").arg(threads.max(1).to_string());
    cmd.arg("-q");
    cmd.stdout(out).stderr(err);

    let status = cmd
        .status()
        .with_context(|| format!("spawning {} for the adapter guard", check_bin.display()))?;
    if !status.success() {
        return Ok(None);
    }
    let text = match std::fs::read_to_string(json) {
        Ok(t) => t,
        Err(_) => return Ok(None),
    };
    let doc: serde_json::Value = match serde_json::from_str(&text) {
        Ok(v) => v,
        Err(_) => return Ok(None),
    };
    Ok(Measurement::from_report(&doc))
}

/// The guard's result for one trimming step, for the message and the manifest.
#[derive(Debug, Clone)]
pub struct GuardOutcome {
    pub verdict: Verdict,
    pub raw: Option<Measurement>,
    pub trimmed: Option<Measurement>,
}

impl GuardOutcome {
    /// The manifest view of this result.
    pub fn to_json(&self) -> serde_json::Value {
        let f = |m: &Option<Measurement>| match m {
            Some(m) => serde_json::json!({
                "reads": m.reads,
                "reads_with_adapter": m.reads_with_adapter,
                "fraction": m.fraction(),
            }),
            None => serde_json::Value::Null,
        };
        serde_json::json!({
            "raw": f(&self.raw),
            "trimmed": f(&self.trimmed),
            "tripped": self.verdict.is_trip(),
            "silent_reason": self.verdict.silent_reason(),
        })
    }
}

/// Guard one trimming step: measure the reads entering it and the reads leaving
/// it and decide. `tag` names the scratch report files (the sample name).
///
/// The scratch reports are written under `log_dir` (the run's `logs/`): that
/// directory is skipped by the report aggregator, so the guard's extra fastcheck
/// reports never turn into bogus samples in the run report.
///
/// Cheap fast-path: the trimmed output is only measured when the input actually
/// carried adapter, so the common adapter-free case costs one fastcheck pass.
pub fn guard_step(
    check_bin: &Path,
    inputs: &[PathBuf],
    outputs: &[PathBuf],
    log_dir: &Path,
    tag: &str,
    policy: &Policy,
    step_args: &[String],
    threads: usize,
) -> Result<GuardOutcome> {
    let cause = Cause::from_step_args(step_args);
    let catalogue = cause.catalogue_arg(step_args);
    let raw_json = log_dir.join(format!("{tag}.adapter-guard.raw.json"));
    let trim_json = log_dir.join(format!("{tag}.adapter-guard.trim.json"));
    let raw_log = log_dir.join(format!("{tag}.adapter-guard.raw.log"));
    let trim_log = log_dir.join(format!("{tag}.adapter-guard.trim.log"));

    let raw = measure(
        check_bin,
        inputs,
        &raw_json,
        &raw_log,
        catalogue.as_deref(),
        threads,
    )?;
    let Some(raw) = raw else {
        return Ok(GuardOutcome {
            verdict: Verdict::Silent("the trimming step's input could not be measured"),
            raw: None,
            trimmed: None,
        });
    };
    // No adapter entering the step: nothing to guard, no second pass.
    if raw.reads_with_adapter == 0 {
        return Ok(GuardOutcome {
            verdict: Verdict::Silent("the trimming step's input carries no detectable adapter"),
            raw: Some(raw),
            trimmed: None,
        });
    }

    let trimmed = measure(
        check_bin,
        outputs,
        &trim_json,
        &trim_log,
        catalogue.as_deref(),
        threads,
    )?;
    let Some(trimmed) = trimmed else {
        return Ok(GuardOutcome {
            verdict: Verdict::Silent("the trimmed output could not be measured"),
            raw: Some(raw),
            trimmed: None,
        });
    };

    let verdict = evaluate(&raw, &trimmed, policy, step_args);
    Ok(GuardOutcome {
        verdict,
        raw: Some(raw),
        trimmed: Some(trimmed),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report(reads: u64, with_match: u64, matches: u64) -> serde_json::Value {
        serde_json::json!({
            "schema": "fastsuite.report.v1",
            "tool": "fastcheck",
            "metrics": { "reads": reads },
            "adapter_content": {
                "error_rate": 0.2,
                "min_overlap": 14,
                "total_matches": matches,
                "adapters": [ { "name": "Nextera", "reads_with_match": with_match, "per_position": {} } ],
            }
        })
    }

    #[test]
    fn measurement_reads_the_block_and_the_read_count() {
        let m = Measurement::from_report(&report(200, 40, 41)).unwrap();
        assert_eq!(m.reads, 200);
        assert_eq!(m.reads_with_adapter, 40);
        assert_eq!(m.total_matches, 41);
        assert!((m.fraction() - 0.2).abs() < 1e-12);
    }

    #[test]
    fn absent_block_is_not_measured_and_not_a_zero() {
        let doc = serde_json::json!({ "metrics": { "reads": 100 } });
        assert!(Measurement::from_report(&doc).is_none(), "absent block = not measured");
        // A block that ran but measured zero reads is also "not measured".
        let empty = serde_json::json!({ "metrics": { "reads": 0 }, "adapter_content": { "adapters": [] } });
        assert!(Measurement::from_report(&empty).is_none());
    }

    #[test]
    fn residual_adapter_warns_and_names_the_causes() {
        let raw = Measurement::from_report(&report(100, 50, 50)).unwrap();
        let trimmed = Measurement::from_report(&report(100, 47, 47)).unwrap();
        let v = evaluate(&raw, &trimmed, &Policy::default(), &[]);
        assert!(v.is_trip(), "adapter survived: {v:?}");
        let msg = v.message().unwrap();
        assert!(msg.contains("survived the trim"), "{msg}");
        // The actionable remedies the task asks for, by name.
        assert!(msg.contains("--adapter-error-rate"), "{msg}");
        assert!(msg.contains("--min-adapter-overlap"), "{msg}");
        assert!(msg.contains("--adapter-preset"), "{msg}");
        assert!(msg.contains("likely cause"), "{msg}");
    }

    #[test]
    fn clean_output_is_silent() {
        let raw = Measurement::from_report(&report(100, 50, 50)).unwrap();
        let trimmed = Measurement::from_report(&report(100, 0, 0)).unwrap();
        let v = evaluate(&raw, &trimmed, &Policy::default(), &[]);
        assert_eq!(v, Verdict::Silent("adapter content fell below the residual threshold"));
        assert!(v.message().is_none());
    }

    #[test]
    fn input_without_adapter_is_silent_even_if_output_would_trip() {
        // The guard must never nag a run whose data had no adapter to begin with.
        let raw = Measurement::from_report(&report(100, 0, 0)).unwrap();
        let trimmed = Measurement::from_report(&report(100, 90, 90)).unwrap();
        let v = evaluate(&raw, &trimmed, &Policy::default(), &[]);
        assert!(!v.is_trip(), "{v:?}");
        assert!(v.silent_reason().unwrap().contains("no detectable adapter"));
    }

    #[test]
    fn fail_policy_turns_a_trip_into_fail() {
        let raw = Measurement::from_report(&report(100, 50, 50)).unwrap();
        let trimmed = Measurement::from_report(&report(100, 40, 40)).unwrap();
        let policy = Policy { fail: true, ..Policy::default() };
        assert!(matches!(evaluate(&raw, &trimmed, &policy, &[]), Verdict::Fail(_)));
    }

    #[test]
    fn tiny_residual_below_threshold_with_a_big_fall_is_silent() {
        // 40% -> 2%: below the 5% floor and a large fall. Silent.
        let raw = Measurement::from_report(&report(1000, 400, 400)).unwrap();
        let trimmed = Measurement::from_report(&report(1000, 20, 20)).unwrap();
        let v = evaluate(&raw, &trimmed, &Policy::default(), &[]);
        assert!(!v.is_trip(), "{v:?}");
    }

    #[test]
    fn a_fall_that_is_not_material_trips_even_below_the_threshold() {
        // 20% -> 18%: under the 5% floor? no, but the point of this case is the
        // relative rule; use 4% -> 3.6% so the absolute floor would NOT trip.
        let raw = Measurement::from_report(&report(1000, 40, 40)).unwrap();
        let trimmed = Measurement::from_report(&report(1000, 36, 36)).unwrap();
        let v = evaluate(&raw, &trimmed, &Policy::default(), &[]);
        assert!(v.is_trip(), "a bare 10% reduction is not material: {v:?}");
    }

    #[test]
    fn threshold_is_configurable() {
        let raw = Measurement::from_report(&report(100, 10, 10)).unwrap();
        let trimmed = Measurement::from_report(&report(100, 10, 10)).unwrap();
        // A threshold above any fraction and no required fall: nothing trips.
        let lenient = Policy { threshold: 1.0, min_reduction: 0.0, ..Policy::default() };
        assert!(!evaluate(&raw, &trimmed, &lenient, &[]).is_trip(), "1.0/0.0 disables the trip");
        let strict = Policy { threshold: 0.05, min_reduction: 0.5, ..Policy::default() };
        assert!(evaluate(&raw, &trimmed, &strict, &[]).is_trip());
    }

    #[test]
    fn cause_is_inferred_from_the_step_args() {
        let no_adapter = vec!["--no-adapter".to_string()];
        assert_eq!(Cause::from_step_args(&no_adapter), Cause::NoAdapterSet);
        let explicit = vec!["--adapter-preset".into(), "nextera".into()];
        assert_eq!(Cause::from_step_args(&explicit), Cause::StrictMatching);
        let fasta = vec!["--adapters".into(), "/tmp/a.fa".into()];
        assert_eq!(Cause::from_step_args(&fasta), Cause::StrictMatching);
        assert_eq!(
            Cause::from_step_args(&fasta).catalogue_arg(&fasta).as_deref(),
            Some("/tmp/a.fa")
        );
        assert_eq!(Cause::from_step_args(&[]), Cause::AutoDetect);
        assert_eq!(Cause::from_step_args(&[]).catalogue_arg(&[]), None);
    }

    #[test]
    fn warns_even_when_only_position_counts_would_show_it() {
        // Multiple adapters: the fraction is the max per-adapter count, never the
        // sum, so a read matching two adapters is not double-counted.
        let doc = serde_json::json!({
            "metrics": { "reads": 100 },
            "adapter_content": {
                "total_matches": 60,
                "adapters": [
                    { "name": "TruSeq/R1", "reads_with_match": 30, "per_position": {} },
                    { "name": "Nextera", "reads_with_match": 25, "per_position": {} },
                ]
            }
        });
        let m = Measurement::from_report(&doc).unwrap();
        assert_eq!(m.reads_with_adapter, 30, "max, not sum");
        assert!((m.fraction() - 0.30).abs() < 1e-12);
    }
}
