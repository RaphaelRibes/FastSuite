//! `--auto-params`: a decision list over a fastcheck QC report.
//!
//! The recommendation is a small JSON file of ordered rules, each an
//! `(field, op, value)` conjunction producing a patch of trimming knobs. The
//! FIRST rule whose conditions all hold wins; if none does, the model's
//! `default` patch applies, and a knob that no patch sets keeps its shipped
//! default. Walking it is a handful of float comparisons — no inference
//! runtime, no new dependency (PLAN-ML-STEPS F1).
//!
//! Why a decision list and not a tree ensemble: the recommendation has to be
//! inspectable and auditable. `--print-args` prints which rule fired and what
//! it set, so a user can read the reasoning and disagree with it. A model whose
//! recommendation cannot be explained is not usable on a trimming tool whose
//! failure mode is silently wrong reads.
//!
//! The features are the fields `fastcheck` measures on the INPUT before
//! trimming (`QC_FEATURES` in `bench/scripts/labels.py`): the predictor sees
//! exactly what a user could see by running `fastcheck` first. Nothing measured
//! on the trimmed output is allowed in, or the model would be predicting from
//! the answer.

use anyhow::{bail, Context, Result};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use crate::algorithm::MinLen;

/// Tag of the model format. A model with any other tag is refused rather than
/// half-understood: the rules of v1 and a future v2 do not share semantics.
pub const SCHEMA: &str = "fastsuite.auto-params.v1";

/// Default file name, searched next to the binary's `models/` directory and in
/// the working directory, and overridable with `--model` / `$FASTPREP_MODEL`.
pub const MODEL_FILE: &str = "auto-params.v1.json";

/// Environment override for the model path.
pub const MODEL_ENV: &str = "FASTPREP_MODEL";

/// The model baked into the binary: what a bare `--auto-params` uses when no
/// file is found. Embedding it (rather than only shipping the file) is what
/// makes the flag work out of the box — a `cargo install` copies one binary,
/// not a tidy tree — while the same bytes stay on disk
/// (`models/auto-params.v1.json`) so the rules can be read and audited without
/// a build. `the_embedded_model_is_the_file_on_disk` keeps the two identical.
pub const EMBEDDED_MODEL: &str = include_str!("../models/auto-params.v1.json");

/// Display label for the embedded model, used where a path is printed.
pub const EMBEDDED_LABEL: &str = "<built-in>";

/// The QC fields a rule may test. These are fastcheck's `metrics` keys on the
/// real input — the same list `bench/scripts/labels.py` builds the training
/// matrix from. A rule naming anything else is refused at load time: a
/// misspelled field would otherwise be a rule that silently never fires.
pub const KNOWN_FEATURES: [&str; 12] = [
    "read_len",
    "gc_fraction",
    "mean_read_mean_q",
    "reads_q20_frac",
    "reads_q30_frac",
    "polyx5_frac_gt_1",
    "polyx3_frac_gt_5",
    "head_g_frac",
    "head_dropout",
    "duplication_level",
    "kmer_top_frac",
    "bases_non_acgt_frac",
];

/// The knobs a patch may set: every axis the Phase D grid swept that is also a
/// fastprep flag. `r2_reverse` is deliberately absent — the sweep reversed R2
/// in the BENCH pipeline to simulate the other strand convention; it is not a
/// fastprep option, so a model that "recommended" it would be recommending a
/// flag that does not exist. `head_dropout` is here but was never swept (D-16);
/// a model may set it only if the trainer measured it.
pub const KNOBS: [&str; 9] = [
    "min_run",
    "mercy",
    "phred_sat",
    "gate_shape",
    "min_len",
    "mean_q",
    "adapter_error_rate",
    "min_adapter_overlap",
    "interior",
];

// ---------------------------------------------------------------------------
// Comparison operators.
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
}

impl Op {
    /// Whether `lhs <op> rhs` holds. `eq`/`ne` compare exactly, which is what a
    /// trained threshold on an integer-valued field (a read length, a count)
    /// means; on a float field the trainer emits inequalities.
    pub fn holds(self, lhs: f64, rhs: f64) -> bool {
        match self {
            Op::Lt => lhs < rhs,
            Op::Le => lhs <= rhs,
            Op::Gt => lhs > rhs,
            Op::Ge => lhs >= rhs,
            Op::Eq => lhs == rhs,
            Op::Ne => lhs != rhs,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Op::Lt => "<",
            Op::Le => "<=",
            Op::Gt => ">",
            Op::Ge => ">=",
            Op::Eq => "==",
            Op::Ne => "!=",
        }
    }
}

impl FromStr for Op {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "<" | "lt" => Ok(Op::Lt),
            "<=" | "le" => Ok(Op::Le),
            ">" | "gt" => Ok(Op::Gt),
            ">=" | "ge" => Ok(Op::Ge),
            "==" | "=" | "eq" => Ok(Op::Eq),
            "!=" | "ne" => Ok(Op::Ne),
            other => Err(format!(
                "unknown operator {other:?} (expected one of < <= > >= == !=)"
            )),
        }
    }
}

impl<'de> Deserialize<'de> for Op {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        Op::from_str(&s).map_err(serde::de::Error::custom)
    }
}

// ---------------------------------------------------------------------------
// Model.
// ---------------------------------------------------------------------------

/// One `(field, op, value)` test.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Condition {
    pub field: String,
    pub op: Op,
    pub value: f64,
}

/// One rule: a conjunction of conditions and the patch it applies.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    /// Empty `when` means "always", which is only legitimate for the last rule;
    /// the trainer emits no such rule (its fallback is `default`), and a model
    /// that does gets a warning in `--print-args` rather than an error.
    #[serde(default)]
    pub when: Vec<Condition>,
    #[serde(default)]
    pub params: ParamsPatch,
    /// How many training datasets the rule fired on. Reported, never used for
    /// the decision: a rule with coverage 1 still fires at predict time.
    #[serde(default)]
    pub coverage: Option<u32>,
    /// Mean regret of the rule's config on those datasets (0.0 = the winner
    /// every time). Reported by `--print-args` so a thin rule is visible.
    #[serde(default)]
    pub regret: Option<f64>,
    #[serde(default)]
    pub note: Option<String>,
}

impl Rule {
    /// Whether every condition holds. A field the profile does not have makes
    /// the condition FALSE — absence is not agreement. That is what keeps a QC
    /// report missing `polyx3_frac_gt_5` from selecting a rule trained on it.
    pub fn holds(&self, profile: &Profile) -> bool {
        self.when.iter().all(|c| match profile.get(&c.field) {
            Some(v) => c.op.holds(v, c.value),
            None => false,
        })
    }

    /// `mean_read_mean_q < 28 && polyx3_frac_gt_5 >= 0.05`, for `--print-args`.
    pub fn describe(&self) -> String {
        if self.when.is_empty() {
            return "always".to_string();
        }
        self.when
            .iter()
            .map(|c| format!("{} {} {}", c.field, c.op.as_str(), fmt_num(c.value)))
            .collect::<Vec<_>>()
            .join(" && ")
    }
}

/// A patch of trimming knobs. Every field is optional: absent means "keep what
/// the CLI (or the shipped default) says", which is what makes a partial rule
/// readable.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParamsPatch {
    #[serde(default)]
    pub min_run: Option<usize>,
    #[serde(default)]
    pub mercy: Option<f64>,
    #[serde(default)]
    pub phred_sat: Option<u8>,
    #[serde(default)]
    pub gate_shape: Option<f64>,
    #[serde(default)]
    pub min_len: Option<MinLenSpec>,
    /// `null` means `--no-mean-q`; absent means "leave the gate alone".
    #[serde(default, deserialize_with = "double_option")]
    pub mean_q: Option<Option<f64>>,
    #[serde(default)]
    pub adapter_error_rate: Option<f64>,
    #[serde(default)]
    pub min_adapter_overlap: Option<usize>,
    #[serde(default)]
    pub head_dropout: Option<usize>,
    #[serde(default)]
    pub interior: Option<bool>,
}

impl ParamsPatch {
    /// Whether the patch sets nothing (a rule that recommends the defaults).
    pub fn is_empty(&self) -> bool {
        self.min_run.is_none()
            && self.mercy.is_none()
            && self.phred_sat.is_none()
            && self.gate_shape.is_none()
            && self.min_len.is_none()
            && self.mean_q.is_none()
            && self.adapter_error_rate.is_none()
            && self.min_adapter_overlap.is_none()
            && self.head_dropout.is_none()
            && self.interior.is_none()
    }

    /// The patch as a space-joined argv fragment, in a stable order.
    pub fn to_argv(&self) -> Vec<String> {
        let mut out = Vec::new();
        if let Some(v) = self.min_run {
            out.push(format!("--min-run {v}"));
        }
        if let Some(v) = self.mercy {
            out.push(format!("--mercy {}", fmt_num(v)));
        }
        if let Some(v) = self.phred_sat {
            out.push(format!("--phred-sat {v}"));
        }
        if let Some(v) = self.gate_shape {
            out.push(format!("--gate-shape {}", fmt_num(v)));
        }
        if let Some(v) = &self.min_len {
            out.push(format!("--min-len {}", v.to_flag()));
        }
        match self.mean_q {
            Some(Some(q)) => out.push(format!("--mean-q {}", fmt_num(q))),
            Some(None) => out.push("--no-mean-q".to_string()),
            None => {}
        }
        if let Some(v) = self.adapter_error_rate {
            out.push(format!("--adapter-error-rate {}", fmt_num(v)));
        }
        if let Some(v) = self.min_adapter_overlap {
            out.push(format!("--min-adapter-overlap {v}"));
        }
        if let Some(v) = self.head_dropout {
            out.push(format!("--head-dropout {v}"));
        }
        if let Some(true) = self.interior {
            out.push("--interior".to_string());
        }
        out
    }

    /// Resolve the patch's `min_len` spec, if any.
    pub fn min_len_value(&self) -> Result<Option<MinLen>> {
        self.min_len.as_ref().map(|s| s.resolve()).transpose()
    }

    /// Reject values that the CLI itself would refuse.
    pub fn validate(&self, where_: &str) -> Result<()> {
        if let Some(v) = self.min_run {
            if v == 0 {
                bail!("{where_}: min_run must be > 0");
            }
        }
        if let Some(v) = self.mercy {
            if !(v.is_finite() && v > 0.0) {
                bail!("{where_}: mercy must be finite and > 0, got {v}");
            }
        }
        if let Some(v) = self.gate_shape {
            if !v.is_finite() {
                bail!("{where_}: gate_shape must be finite, got {v}");
            }
        }
        if let Some(v) = self.adapter_error_rate {
            if !(0.0..1.0).contains(&v) {
                bail!("{where_}: adapter_error_rate must be in [0, 1), got {v}");
            }
        }
        if let Some(v) = &self.min_len {
            v.resolve().with_context(|| format!("{where_}: bad min_len"))?;
        }
        Ok(())
    }
}

/// `min_len` in a model file: `"0.9L"` (a fraction of the insert) or `100` (an
/// absolute floor). Both forms are what `--min-len` accepts.
#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
pub enum MinLenSpec {
    Fraction(String),
    Absolute(usize),
}

impl MinLenSpec {
    pub fn resolve(&self) -> Result<MinLen> {
        match self {
            MinLenSpec::Fraction(s) => {
                MinLen::from_str(s).map_err(|e| anyhow::anyhow!("{e}"))
            }
            MinLenSpec::Absolute(n) => {
                if *n == 0 {
                    bail!("min_len must be > 0");
                }
                Ok(MinLen::Abs(*n))
            }
        }
    }

    /// Canonical flag spelling: `0.9L` or `100`.
    pub fn to_flag(&self) -> String {
        match self {
            MinLenSpec::Fraction(s) => s.clone(),
            MinLenSpec::Absolute(n) => n.to_string(),
        }
    }
}

/// The fallback when no rule matches: a patch, and an optional note.
///
/// A wrapper (rather than a bare patch) so the file reads the same way as a
/// rule: `{"params": {...}}`. A model with no `default` block at all is legal
/// and means "the shipped defaults".
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DefaultRule {
    #[serde(default)]
    pub params: ParamsPatch,
    #[serde(default)]
    pub note: Option<String>,
}

/// A trained model file.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub schema: String,
    /// Free-form provenance, printed verbatim by `--print-args`. The trainer
    /// writes the corpus, the dataset count, the domains and the date here.
    #[serde(default)]
    pub trained_on: Option<serde_json::Value>,
    /// The feature contract, in the trainer's column order. Informational, but
    /// checked non-empty: a model file with no feature list is a model file
    /// whose walk was never meant to be reproduced.
    pub features: Vec<String>,
    #[serde(default)]
    pub rules: Vec<Rule>,
    #[serde(default)]
    pub default: Option<DefaultRule>,
    /// The measured acceptance block (held-out regret, top-5 agreement, the
    /// fastp comparison). Printed by `--print-args`, because a recommendation
    /// without its measured accuracy is a guess with a schema.
    #[serde(default)]
    pub acceptance: Option<serde_json::Value>,
}

impl Model {
    pub fn load(path: &Path) -> Result<Model> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading model {}", path.display()))?;
        let model: Model = serde_json::from_str(&text)
            .with_context(|| format!("parsing model {}", path.display()))?;
        model
            .validate()
            .with_context(|| format!("model {} is not usable", path.display()))?;
        Ok(model)
    }

    /// Structural checks that catch a malformed or stale model file.
    pub fn validate(&self) -> Result<()> {
        if self.schema != SCHEMA {
            bail!("schema is {:?}, expected {:?}", self.schema, SCHEMA);
        }
        if self.features.is_empty() {
            bail!("`features` is empty: the file does not say what it was trained on");
        }
        for (i, r) in self.rules.iter().enumerate() {
            let where_ = format!("rule {i}");
            for c in &r.when {
                if !KNOWN_FEATURES.contains(&c.field.as_str()) {
                    bail!(
                        "{where_}: unknown feature {:?} — the walk could never test it \
                         (known: {})",
                        c.field,
                        KNOWN_FEATURES.join(", ")
                    );
                }
                if !c.value.is_finite() {
                    bail!("{where_}: threshold for {:?} is not finite", c.field);
                }
            }
            if r.when.is_empty() && i + 1 != self.rules.len() {
                bail!(
                    "{where_}: an unconditional rule before the last one makes every \
                     later rule dead — move it to the end or give it a condition"
                );
            }
            r.params.validate(&where_)?;
        }
        if let Some(d) = &self.default {
            d.params.validate("default")?;
        }
        Ok(())
    }

    /// The provenance line(s) for `--print-args`.
    pub fn provenance(&self) -> String {
        match &self.trained_on {
            Some(v) => serde_json::to_string(v).unwrap_or_else(|_| v.to_string()),
            None => "(no provenance block in the model file)".to_string(),
        }
    }

    /// Walk the decision list. First matching rule wins; otherwise `default`.
    pub fn walk(&self, profile: &Profile) -> Recommendation {
        for (i, rule) in self.rules.iter().enumerate() {
            if rule.holds(profile) {
                return Recommendation {
                    patch: rule.params.clone(),
                    fired: Some(i),
                    rule: Some(RuleInfo {
                        describe: rule.describe(),
                        coverage: rule.coverage,
                        regret: rule.regret,
                        note: rule.note.clone(),
                    }),
                    default_fired: false,
                };
            }
        }
        Recommendation {
            patch: self
                .default
                .as_ref()
                .map(|d| d.params.clone())
                .unwrap_or_default(),
            fired: None,
            rule: None,
            default_fired: true,
        }
    }
}

/// Where a recommendation came from, for the printed explanation.
#[derive(Debug, Clone)]
pub struct RuleInfo {
    pub describe: String,
    pub coverage: Option<u32>,
    pub regret: Option<f64>,
    pub note: Option<String>,
}

/// The result of a walk.
#[derive(Debug, Clone)]
pub struct Recommendation {
    pub patch: ParamsPatch,
    pub fired: Option<usize>,
    pub rule: Option<RuleInfo>,
    /// True when no rule matched and the model's `default` patch was used.
    pub default_fired: bool,
}

impl Recommendation {
    /// One line saying which rule fired and what it set.
    pub fn explain(&self, n_rules: usize) -> String {
        match (&self.rule, self.fired) {
            (Some(info), Some(i)) => {
                let mut s = format!("rule {}/{}: {}", i + 1, n_rules, info.describe);
                if let Some(c) = info.coverage {
                    s.push_str(&format!(" [n={c}"));
                    if let Some(r) = info.regret {
                        s.push_str(&format!(", regret {r:.4}"));
                    }
                    s.push(']');
                }
                if let Some(n) = &info.note {
                    s.push_str(&format!(" — {n}"));
                }
                s
            }
            _ => format!("no rule matched {} rule(s): using the model's default", n_rules),
        }
    }
}

// ---------------------------------------------------------------------------
// Profile: the QC report's numbers.
// ---------------------------------------------------------------------------

/// The QC numbers a rule can test, keyed by feature name.
#[derive(Debug, Clone, Default)]
pub struct Profile {
    fields: BTreeMap<String, f64>,
    /// Keys in the report that are NOT known features: reported so a QC report
    /// that names a field this build does not know is visible rather than
    /// silently ignored.
    pub unknown: Vec<String>,
}

impl Profile {
    pub fn get(&self, field: &str) -> Option<f64> {
        self.fields.get(field).copied()
    }

    pub fn len(&self) -> usize {
        self.fields.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    pub fn set(&mut self, field: &str, value: f64) {
        self.fields.insert(field.to_string(), value);
    }

    /// Build a profile from a fastcheck report.
    ///
    /// Accepts both shapes in use: the `fastsuite.report.v1` wrapper with a
    /// `metrics` block, and a bare flat object of numbers. Non-numeric values
    /// inside `metrics` are refused (they would be silently dropped otherwise,
    /// and a dropped field turns a rule off without saying so).
    pub fn from_qc_json(text: &str) -> Result<Profile> {
        let doc: serde_json::Value = serde_json::from_str(text).context("parsing the QC JSON")?;
        let obj = doc
            .get("metrics")
            .and_then(|m| m.as_object())
            .or_else(|| doc.as_object())
            .context("the QC JSON has no `metrics` object and is not a flat object")?;
        let mut p = Profile::default();
        for (k, v) in obj {
            if let Some(n) = v.as_f64() {
                p.fields.insert(k.clone(), n);
            } else if v.is_string() || v.is_object() || v.is_array() || v.is_boolean() {
                // `tool`, `schema`, nested blocks: not features, not errors.
                if !v.is_object() && !v.is_array() {
                    p.unknown.push(k.clone());
                }
            } else if v.is_null() {
                // A null metric is a field the tool could not measure. Leaving
                // it out is the documented "absence is not agreement" path.
                p.unknown.push(k.clone());
            }
        }
        if p.fields.is_empty() {
            bail!("the QC JSON contains no numeric metrics");
        }
        // Two features the training matrix builds from OTHER report keys, and
        // the alias must live here or the model's `read_len` rules could never
        // fire on a real report — the walk would parse fine, match nothing, and
        // fall through to the default without saying so (train/serve skew, the
        // one failure mode that looks like success):
        //   read_len            <- mean_len            (labels.py:85)
        //   bases_non_acgt_frac <- bases_non_acgt / bases  (labels.py:162)
        if !p.fields.contains_key("read_len") {
            if let Some(v) = p.fields.get("mean_len").copied() {
                p.fields.insert("read_len".to_string(), v);
            }
        }
        if !p.fields.contains_key("bases_non_acgt_frac") {
            if let (Some(n), Some(b)) = (
                p.fields.get("bases_non_acgt").copied(),
                p.fields.get("bases").copied(),
            ) {
                if b > 0.0 {
                    p.fields.insert("bases_non_acgt_frac".to_string(), n / b);
                }
            }
        }
        p.fields.retain(|k, _| KNOWN_FEATURES.contains(&k.as_str()));
        Ok(p)
    }

    pub fn from_file(path: &Path) -> Result<Profile> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading QC report {}", path.display()))?;
        Profile::from_qc_json(&text)
            .with_context(|| format!("QC report {}", path.display()))
    }
}

// ---------------------------------------------------------------------------
// Model discovery.
// ---------------------------------------------------------------------------

/// The paths searched for a model file, in order. `--model` is checked by the
/// caller (it is explicit); this is the implicit search that a bare
/// `--auto-params` uses.
pub fn search_paths() -> Vec<PathBuf> {
    let mut out = Vec::new();
    if let Ok(p) = std::env::var(MODEL_ENV) {
        if !p.is_empty() {
            out.push(PathBuf::from(p));
        }
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            out.push(dir.join("models").join(MODEL_FILE));
            out.push(dir.join(MODEL_FILE));
        }
    }
    out.push(PathBuf::from(MODEL_FILE));
    out
}

/// The embedded model, parsed and validated.
pub fn embedded_model() -> Result<Model> {
    let model: Model = serde_json::from_str(EMBEDDED_MODEL)
        .context("the built-in model does not parse — the binary is corrupt")?;
    model
        .validate()
        .context("the built-in model is not usable")?;
    Ok(model)
}

/// First readable model in [`search_paths`], else the embedded one.
///
/// A file on disk always wins: a user who drops their own model next to the
/// binary has said which model they mean.
pub fn find_model() -> Result<(PathBuf, Model)> {
    for p in search_paths() {
        if p.is_file() {
            return Ok((p.clone(), Model::load(&p)?));
        }
    }
    Ok((PathBuf::from(EMBEDDED_LABEL), embedded_model()?))
}

/// `4.0` -> `4`, `0.05` -> `0.05`. A flag value the user would have typed.
pub fn fmt_num(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

fn double_option<'de, D, T>(de: D) -> Result<Option<Option<T>>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    Option::<T>::deserialize(de).map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model_json() -> &'static str {
        r#"{
          "schema": "fastsuite.auto-params.v1",
          "trained_on": {"corpus": "unit-test", "datasets": 3},
          "features": ["mean_read_mean_q", "polyx3_frac_gt_5", "read_len"],
          "rules": [
            {"when": [{"field": "mean_read_mean_q", "op": "<", "value": 28.0},
                      {"field": "polyx3_frac_gt_5", "op": ">=", "value": 0.05}],
             "params": {"mercy": 5.0, "min_run": 16, "min_len": "0.85L"},
             "coverage": 12, "regret": 0.0031, "note": "degraded, tailed"},
            {"when": [{"field": "read_len", "op": ">=", "value": 140.0}],
             "params": {"min_len": 100},
             "coverage": 40},
            {"when": [{"field": "mean_read_mean_q", "op": ">", "value": 34.0}],
             "params": {"mean_q": null, "interior": true}}
          ],
          "default": {"params": {"mercy": 4.0}}
        }"#
    }

    fn model() -> Model {
        let m: Model = serde_json::from_str(model_json()).expect("fixture parses");
        m.validate().expect("fixture validates");
        m
    }

    #[test]
    fn ops_read_both_spellings_and_apply() {
        for (s, op) in [
            ("<", Op::Lt),
            ("lt", Op::Lt),
            ("<=", Op::Le),
            (">=", Op::Ge),
            ("==", Op::Eq),
            ("!=", Op::Ne),
            (">", Op::Gt),
        ] {
            assert_eq!(Op::from_str(s).unwrap(), op);
        }
        assert!(Op::Lt.holds(1.0, 2.0));
        assert!(!Op::Lt.holds(2.0, 2.0));
        assert!(Op::Ge.holds(2.0, 2.0));
        assert!(Op::Ne.holds(1.0, 2.0));
        assert!(Op::from_str("~").is_err());
    }

    #[test]
    fn the_first_matching_rule_wins() {
        let m = model();
        let mut p = Profile::default();
        p.set("mean_read_mean_q", 26.0);
        p.set("polyx3_frac_gt_5", 0.2);
        p.set("read_len", 150.0);
        let r = m.walk(&p);
        assert_eq!(r.fired, Some(0), "rule 0 matches and must win over rule 1");
        assert_eq!(r.patch.mercy, Some(5.0));
        assert_eq!(r.patch.min_run, Some(16));
        assert_eq!(r.patch.min_len_value().unwrap(), Some(MinLen::Frac(0.85)));
    }

    #[test]
    fn a_later_rule_applies_when_the_earlier_one_does_not() {
        let m = model();
        let mut p = Profile::default();
        p.set("mean_read_mean_q", 36.0);
        p.set("polyx3_frac_gt_5", 0.0);
        p.set("read_len", 150.0);
        let r = m.walk(&p);
        // Rule 1 (read_len >= 140) matches before rule 2 (mean Q > 34).
        assert_eq!(r.fired, Some(1));
        assert_eq!(r.patch.min_len_value().unwrap(), Some(MinLen::Abs(100)));
    }

    #[test]
    fn a_missing_field_is_not_agreement() {
        let m = model();
        let mut p = Profile::default();
        // No mean_read_mean_q at all: rule 0 and rule 2 cannot fire.
        p.set("polyx3_frac_gt_5", 0.9);
        let r = m.walk(&p);
        assert_eq!(r.fired, None, "a rule testing a field the profile lacks must not fire");
        assert!(r.default_fired);
        assert_eq!(r.patch.mercy, Some(4.0), "the model default applies");
    }

    #[test]
    fn mean_q_null_means_no_mean_q_and_absent_means_leave_it() {
        let m = model();
        let mut p = Profile::default();
        p.set("mean_read_mean_q", 40.0);
        let r = m.walk(&p);
        assert_eq!(r.fired, Some(2));
        assert_eq!(r.patch.mean_q, Some(None), "null in the model = --no-mean-q");
        assert!(r.patch.to_argv().contains(&"--no-mean-q".to_string()));
        let second = &m.rules[1].params;
        assert_eq!(second.mean_q, None, "absent = leave the gate alone");
        assert!(!second.to_argv().iter().any(|a| a.starts_with("--mean-q")));
    }

    #[test]
    fn argv_is_stable_and_round_trips_through_the_cli_parser() {
        let m = model();
        let mut p = Profile::default();
        p.set("mean_read_mean_q", 26.0);
        p.set("polyx3_frac_gt_5", 0.2);
        let argv = m.walk(&p).patch.to_argv();
        assert_eq!(
            argv,
            vec!["--min-run 16", "--mercy 5", "--min-len 0.85L"],
            "order and formatting are part of the contract"
        );
        // Every fragment must be accepted by the real CLI: split on spaces and
        // parse with clap, which is what a user copy-pasting the output does.
        let mut args: Vec<String> = vec!["fastprep".into()];
        for frag in &argv {
            args.extend(frag.split(' ').map(str::to_string));
        }
        args.extend(["-1", "in.fq", "-o", "out.fq"].map(str::to_string));
        let cli = <crate::cli::Cli as clap::Parser>::try_parse_from(&args)
            .expect("the printed argv must parse");
        assert_eq!(cli.min_run, 16);
        assert_eq!(cli.mercy, 5.0);
        assert_eq!(cli.min_len, MinLen::Frac(0.85));
    }

    #[test]
    fn an_unknown_feature_is_refused_at_load_time() {
        let bad = model_json().replace("\"read_len\"", "\"read_lenght\"");
        let m: Model = serde_json::from_str(&bad).unwrap();
        let err = m.validate().unwrap_err().to_string();
        assert!(err.contains("unknown feature"), "{err}");
        assert!(err.contains("read_lenght"), "{err}");
    }

    #[test]
    fn a_wrong_schema_tag_is_refused() {
        let bad = model_json().replace("fastsuite.auto-params.v1", "fastsuite.auto-params.v0");
        let m: Model = serde_json::from_str(&bad).unwrap();
        assert!(m.validate().unwrap_err().to_string().contains("schema"));
    }

    #[test]
    fn an_unconditional_rule_before_the_last_is_refused() {
        let bad = model_json().replace(
            r#"{"when": [{"field": "read_len", "op": ">=", "value": 140.0}],"#,
            r#"{"when": [],"#,
        );
        let m: Model = serde_json::from_str(&bad).unwrap();
        let err = m.validate().unwrap_err().to_string();
        assert!(err.contains("dead"), "{err}");
    }

    #[test]
    fn out_of_range_knobs_are_refused_like_the_cli_would() {
        let bad = model_json().replace("\"min_run\": 16", "\"min_run\": 0");
        let m: Model = serde_json::from_str(&bad).unwrap();
        assert!(m.validate().unwrap_err().to_string().contains("min_run"));
        let bad = model_json().replace("\"mercy\": 5.0", "\"mercy\": -1.0");
        let m: Model = serde_json::from_str(&bad).unwrap();
        assert!(m.validate().unwrap_err().to_string().contains("mercy"));
    }

    #[test]
    fn a_qc_report_in_either_shape_yields_the_same_profile() {
        let wrapped = r#"{"schema":"fastsuite.report.v1","tool":"fastcheck",
            "metrics":{"read_len":150,"mean_read_mean_q":31.5,"polyx3_frac_gt_5":0.02}}"#;
        let flat = r#"{"read_len":150,"mean_read_mean_q":31.5,"polyx3_frac_gt_5":0.02}"#;
        let a = Profile::from_qc_json(wrapped).unwrap();
        let b = Profile::from_qc_json(flat).unwrap();
        assert_eq!(a.len(), 3);
        assert_eq!(a.get("read_len"), b.get("read_len"));
        assert_eq!(a.get("read_len"), Some(150.0));
    }

    #[test]
    fn a_null_metric_is_absent_not_zero() {
        let text = r#"{"metrics":{"read_len":150,"head_g_frac":null}}"#;
        let p = Profile::from_qc_json(text).unwrap();
        assert_eq!(p.get("head_g_frac"), None, "a null must not become 0.0");
        assert!(p.unknown.contains(&"head_g_frac".to_string()));
    }

    #[test]
    fn fastcheck_aliases_are_derived_like_the_training_matrix() {
        // The real report has mean_len/bases/bases_non_acgt, not the training
        // column names. Without the derivation, `read_len` rules never fire.
        let text = r#"{"schema":"fastsuite.report.v1","tool":"fastcheck","metrics":{
            "mean_len":101,"bases":10100,"bases_non_acgt":101,
            "mean_read_mean_q":31.2,"polyx3_frac_gt_5":0.02,"polyx5_frac_gt_1":0.1}}"#;
        let p = Profile::from_qc_json(text).unwrap();
        assert_eq!(p.get("read_len"), Some(101.0), "read_len must alias mean_len");
        let f = p.get("bases_non_acgt_frac").unwrap();
        assert!((f - 0.01).abs() < 1e-9, "derived fraction was {f}");
        assert_eq!(p.get("mean_len"), None, "the raw alias is not a feature");
    }

    #[test]
    fn a_report_with_no_numbers_is_refused() {
        let err = Profile::from_qc_json(r#"{"metrics":{"tool":"fastcheck"}}"#).unwrap_err();
        assert!(err.to_string().contains("no numeric metrics"), "{err}");
    }

    #[test]
    fn the_embedded_model_is_the_file_on_disk() {
        // The binary embeds models/auto-params.v1.json at compile time; the
        // committed file must be the same bytes, or `--print-args` would show
        // rules a reader cannot find in the repo.
        let on_disk = include_str!("../models/auto-params.v1.json");
        assert_eq!(on_disk, EMBEDDED_MODEL);
        let m = embedded_model().expect("the shipped model must parse and validate");
        assert!(!m.rules.is_empty(), "the shipped model must have rules");
        assert_eq!(m.schema, SCHEMA);
        // Every rule must be satisfiable and must not be shadowed by a LATER
        // rule (first match wins): walking its own conditions lands on this
        // rule, or on an earlier one.
        for (i, rule) in m.rules.iter().enumerate() {
            let mut p = Profile::default();
            for c in &rule.when {
                let v = match c.op {
                    Op::Le | Op::Lt | Op::Eq => c.value,
                    Op::Gt | Op::Ge | Op::Ne => c.value + 1.0,
                };
                p.set(&c.field, v);
            }
            let rec = m.walk(&p);
            let fired = rec.fired.expect("a rule's own conditions must fire something");
            assert!(
                fired <= i,
                "rule {i} is shadowed by a later rule {fired}: the list is not in order"
            );
        }
    }

    #[test]
    fn model_discovery_falls_back_to_the_embedded_model() {
        // No model file in this test's CWD: the embedded model must be used,
        // never an error and never an empty recommendation.
        std::env::remove_var(MODEL_ENV);
        let (p, m) = find_model().expect("there is always a model: the embedded one");
        assert!(
            p.is_file() || p == PathBuf::from(EMBEDDED_LABEL),
            "the source must be a real file or the built-in label, got {}",
            p.display()
        );
        assert!(!m.rules.is_empty());
    }
}
