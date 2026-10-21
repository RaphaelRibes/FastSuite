//! Pipeline statistics.
//!
//! A single [`Stats`] accumulates per read on the ordered (main) thread as
//! results arrive from the compute thread — there is no per-chunk
//! accumulation and no merge step. At program end, rendered to stderr.

use std::collections::HashMap;
use std::time::Duration;

use crate::algorithm::{DropReason, TrimResult};

/// The schema string every FastSuite tool writes with `--json`.
///
/// Copied from the canonical v1 emitter in `tools/fastdedup/src/main.rs`
/// (the `serde_json::json!({...})` block); `tools/fastsim/src/main.rs`
/// repeats the same shape with `outputs.gt` (DECISIONS D-1: copy the
/// reference implementation and name the source).
pub const REPORT_SCHEMA: &str = "fastsuite.report.v1";

// ---------------------------------------------------------------------------
// Per-adapter stats (counts by adapter name).
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct AdapterStats {
    /// Trimmed-read count per adapter index, indexed exactly like the
    /// configured set (`crate::adapter::PE_OVERLAP_ADAPTER` for the PE-overlap
    /// boundary). A flat array, not a hash of names: this is written once per
    /// trimmed read on the ordered thread, and hashing a `String` there was
    /// pure cost — the name is only needed when the report is printed.
    pub counts: [u64; 256],
    /// Names by adapter index, filled ONCE per run by [`AdapterStats::set_names`].
    /// Index `255` is `"pe-overlap"`.
    pub names: Vec<String>,
    pub reads_trimmed: u64,
    pub bases_trimmed_5: u64,
    pub bases_trimmed_3: u64,
}

impl Default for AdapterStats {
    fn default() -> Self {
        Self {
            counts: [0; 256],
            names: Vec::new(),
            reads_trimmed: 0,
            bases_trimmed_5: 0,
            bases_trimmed_3: 0,
        }
    }
}

impl AdapterStats {
    /// Install the adapter names for the run. Called once from the pipeline,
    /// with the configured set's names in set order and the reserved
    /// `"pe-overlap"` slot filled in.
    pub fn set_names(&mut self, mut names: Vec<String>) {
        let reserved = crate::adapter::PE_OVERLAP_ADAPTER as usize;
        if names.len() <= reserved {
            names.resize_with(reserved + 1, String::new);
        }
        names[reserved] = "pe-overlap".to_string();
        self.names = names;
    }

    /// Count one trimmed read against adapter `idx`.
    pub fn record(&mut self, idx: u8, t5: usize, t3: usize) {
        self.counts[idx as usize] += 1;
        self.reads_trimmed += 1;
        self.bases_trimmed_5 += t5 as u64;
        self.bases_trimmed_3 += t3 as u64;
    }

    /// Name for an adapter index; empty when the index was never named (only
    /// possible in tests that record without calling `set_names`).
    pub fn name(&self, idx: usize) -> &str {
        self.names.get(idx).map(String::as_str).unwrap_or("")
    }

    /// The legacy `by_name` map, built at report time: index → count for every
    /// adapter that actually trimmed a read.
    pub fn by_name(&self) -> HashMap<String, u64> {
        let mut m = HashMap::new();
        for (idx, &n) in self.counts.iter().enumerate() {
            if n > 0 {
                m.insert(self.name(idx).to_string(), n);
            }
        }
        m
    }
}

/// Serialised exactly as the old `{by_name, reads_trimmed, bases_trimmed_5,
/// bases_trimmed_3}` object, so the `--json` adapter block does not move. Keys
/// come out sorted because `serde_json::to_value` builds a `BTreeMap` by
/// default, which is also what the old `HashMap` field produced.
impl serde::Serialize for AdapterStats {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut st = s.serialize_struct("AdapterStats", 4)?;
        st.serialize_field("by_name", &self.by_name())?;
        st.serialize_field("reads_trimmed", &self.reads_trimmed)?;
        st.serialize_field("bases_trimmed_5", &self.bases_trimmed_5)?;
        st.serialize_field("bases_trimmed_3", &self.bases_trimmed_3)?;
        st.end()
    }
}

// ---------------------------------------------------------------------------
// Main Stats struct.
// ---------------------------------------------------------------------------

/// One assumption the run made about its input, in the shape D-6 promises:
/// nothing is ever guessed silently, and every assumption appears in the JSON
/// and in the summary line.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Assumption {
    pub what: &'static str,
    pub value: serde_json::Value,
    /// The evidence the assumption rests on.
    pub basis: String,
    /// What goes wrong if the assumption is false.
    pub consequence: String,
}

#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct Stats {
    pub reads_in: u64,
    pub reads_kept: u64,
    pub reads_dropped: u64,
    pub bases_in: u64,
    pub bases_trimmed_5: u64,
    pub bases_trimmed_3: u64,
    pub bases_interior_mass: u64,
    // Per-cause drop counters (sum equals reads_dropped).
    pub reads_dropped_poly5: u64,
    pub reads_dropped_poly3: u64,
    pub reads_dropped_poly_int: u64,
    pub reads_dropped_hard: u64,
    pub reads_dropped_mean_q: u64,
    pub reads_dropped_adapter: u64,
    pub reads_dropped_duplicate: u64,
    /// Post-trim length below `--min-len` (distinct from the trim budget: this
    /// one fires on reads that were simply too short to begin with).
    pub reads_dropped_length: u64,
    /// Phred offset actually used to decode this mate's quality bytes.
    pub phred_offset: u8,
    /// Assumptions the run made. Reported at the top level of the JSON and in
    /// the summary, so a +64 file read as +33 cannot be silent again.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub assumptions: Vec<Assumption>,
    /// Adapter-specific stats; separate to isolate HashMap allocation.
    pub adapter: AdapterStats,
}

impl Stats {
    pub fn record(&mut self, original_len: usize, res: &TrimResult, adapter: Option<u8>) {
        self.reads_in += 1;
        self.bases_in += original_len as u64;
        match res {
            TrimResult::Drop { reason, t5_adapter, t3_adapter, .. } => {
                // Match ergonomics: res is &TrimResult, so fields are &T.
                let (t5a, t3a) = (*t5_adapter, *t3_adapter);
                self.reads_dropped += 1;
                match reason {
                    DropReason::Poly5    => self.reads_dropped_poly5 += 1,
                    DropReason::Poly3    => self.reads_dropped_poly3 += 1,
                    DropReason::PolyInt  => self.reads_dropped_poly_int += 1,
                    DropReason::Hard     => self.reads_dropped_hard += 1,
                    DropReason::MeanQ    => self.reads_dropped_mean_q += 1,
                    DropReason::Adapter  => self.reads_dropped_adapter += 1,
                    DropReason::Duplicate => self.reads_dropped_duplicate += 1,
                    DropReason::Length   => self.reads_dropped_length += 1,
                }
                if (t5a + t3a) > 0 {
                    if let Some(idx) = adapter {
                        self.adapter.record(idx, t5a, t3a);
                    }
                }
            }
            TrimResult::Emit { t5_adapter, t3_adapter, t5_poly, t3_poly, runs } => {
                let (t5a, t3a, t5p, t3p) = (*t5_adapter, *t3_adapter, *t5_poly, *t3_poly);
                self.reads_kept += 1;
                self.bases_trimmed_5 += (t5a + t5p) as u64;
                self.bases_trimmed_3 += (t3a + t3p) as u64;
                let h_int: usize = runs.iter().map(|(_, l)| l).sum();
                self.bases_interior_mass += h_int as u64;
                if (t5a + t3a) > 0 {
                    if let Some(idx) = adapter {
                        self.adapter.record(idx, t5a, t3a);
                    }
                }
            }
        }
    }

    pub fn bases_trimmed_total(&self) -> u64 {
        self.bases_trimmed_5 + self.bases_trimmed_3
    }

    /// The flat `metrics` block for one mate (single-end, or one side of a
    /// paired run). Nothing structured goes in here: the per-adapter counts
    /// stay under the `legacy` block, per REPORT-SCHEMA.md rule 1.
    fn metrics_flat(&self) -> serde_json::Value {
        serde_json::json!({
            "reads_in":        self.reads_in,
            "reads_out":       self.reads_kept,
            "reads_dropped":   self.reads_dropped,
            "bases_in":        self.bases_in,
            "bases_trimmed":   self.bases_trimmed_total(),
            "bases_trimmed_5": self.bases_trimmed_5,
            "bases_trimmed_3": self.bases_trimmed_3,
            "reads_dropped_poly5":  self.reads_dropped_poly5,
            "reads_dropped_poly3":  self.reads_dropped_poly3,
            "reads_dropped_poly_int": self.reads_dropped_poly_int,
            "reads_dropped_hard":   self.reads_dropped_hard,
            "reads_dropped_mean_q": self.reads_dropped_mean_q,
            "reads_dropped_adapter": self.reads_dropped_adapter,
            "reads_dropped_duplicate": self.reads_dropped_duplicate,
            "reads_dropped_length": self.reads_dropped_length,
            "reads_adapter_trimmed": self.adapter.reads_trimmed,
            "phred_offset": self.phred_offset,
        })
    }
}

/// Every piece of run context the wrapper records but the counters do not.
pub struct ReportContext<'a> {
    /// `single` or `paired`.
    pub mode: &'a str,
    pub dry_run: bool,
    pub elapsed: Duration,
    pub inputs: serde_json::Value,
    pub outputs: serde_json::Value,
    /// `dedup` block: which anchor window the fingerprint used, and whether it
    /// had to degenerate to an interior window. `None` when dedup is off.
    pub dedup: Option<serde_json::Value>,
}

/// The `fastsuite.report.v1` wrapper for a single-end run (PLAN.md §1.1).
///
/// Shape copied from the canonical emitter in `tools/fastdedup/src/main.rs`
/// (`serde_json::json!({...})` at the `--json` write site); `tools/fastsim/
/// src/main.rs` uses the same pattern. `metrics` is a flat numeric namespace
/// (REPORT-SCHEMA.md rule 1) and the raw `Stats` object is preserved under
/// `legacy` for one release, so a consumer written against the old shape
/// keeps working.
pub fn to_report_v1(
    tool_version: &str,
    command: &str,
    ctx: &ReportContext,
    s: &Stats,
) -> serde_json::Value {
    let mut doc = serde_json::json!({
        "schema":       REPORT_SCHEMA,
        "tool":         "fastprep",
        "tool_version": tool_version,
        "command":      command,
        "mode":         ctx.mode,
        "dry_run":      ctx.dry_run,
        "inputs":       ctx.inputs,
        "outputs":      ctx.outputs,
        "elapsed_s":    ctx.elapsed.as_secs_f64(),
        "metrics":      s.metrics_flat(),
        // Every assumption the run made, in D-6's shape.
        "assumptions":  &s.assumptions,
        // Raw stats object, exactly as fastprep used to write it.
        "legacy":       serde_json::to_value(s).unwrap_or(serde_json::Value::Null),
    });
    if let Some(d) = &ctx.dedup {
        doc["dedup"] = d.clone();
    }
    doc
}

/// The `fastsuite.report.v1` wrapper for a paired-end run.
///
/// `metrics` is the SUM over both mates for every read-level counter: a
/// read budget taken from one mate alone silently halves the run
/// (the same trap the legacy fasthelper adapter documents). Where a name
/// would be ambiguous across mates, both are kept, `r2_`-prefixed
/// (REPORT-SCHEMA.md rule 2).
pub fn to_report_paired_v1(
    tool_version: &str,
    command: &str,
    ctx: &ReportContext,
    p: &crate::pipeline::PairedStats,
) -> serde_json::Value {
    let mut metrics = serde_json::Map::new();
    let (a, b) = (p.r1.metrics_flat(), p.r2.metrics_flat());
    let (a, b) = (a.as_object().unwrap(), b.as_object().unwrap());
    // Every counter is a u64: sum as u64 so the JSON stays integral
    // (`8`, never `8.0`), matching what the single-end path emits.
    // `phred_offset` is an identity, not a count: summing two mates would
    // report 66 for a +33 run, so it is inserted once below.
    for key in a.keys() {
        if key == "phred_offset" {
            continue;
        }
        let sum = a[key].as_u64().unwrap_or(0) + b[key].as_u64().unwrap_or(0);
        metrics.insert(key.clone(), serde_json::json!(sum));
    }
    metrics.insert("phred_offset".into(), serde_json::json!(p.r1.phred_offset));
    // Per-mate breakdown for the keys that are genuinely mate-specific:
    // r1_* and r2_*, so neither mate is lost when the sum is what a
    // consumer wants.
    for (key, val) in a {
        metrics.insert(format!("r1_{key}"), val.clone());
    }
    for (key, val) in b {
        metrics.insert(format!("r2_{key}"), val.clone());
    }
    // Pair-level counters, from the PairedStats header.
    metrics.insert("pairs_out".into(), serde_json::json!(p.pairs_out));
    metrics.insert("pairs_dropped".into(), serde_json::json!(p.pairs_dropped));
    metrics.insert("singletons_r1".into(), serde_json::json!(p.singletons_r1));
    metrics.insert("singletons_r2".into(), serde_json::json!(p.singletons_r2));

    let mut doc = serde_json::json!({
        "schema":       REPORT_SCHEMA,
        "tool":         "fastprep",
        "tool_version": tool_version,
        "command":      command,
        "mode":         ctx.mode,
        "dry_run":      ctx.dry_run,
        "inputs":       ctx.inputs,
        "outputs":      ctx.outputs,
        "elapsed_s":    ctx.elapsed.as_secs_f64(),
        "metrics":      serde_json::Value::Object(metrics),
        // The offset census is run once, from R1, and applies to both mates.
        "assumptions":  &p.r1.assumptions,
        // Raw PairedStats object, exactly as fastprep used to write it.
        "legacy":       serde_json::to_value(p).unwrap_or(serde_json::Value::Null),
    });
    if let Some(d) = &ctx.dedup {
        doc["dedup"] = d.clone();
    }
    doc
}

// ---------------------------------------------------------------------------
// Reporting.
// ---------------------------------------------------------------------------

pub fn report(label: &str, s: &Stats) {
    eprintln!("  [{}]", label);
    eprintln!("    reads in:                {:>12}", s.reads_in);
    eprintln!(
        "    kept:                    {:>12}  ({})",
        s.reads_kept, pct(s.reads_kept, s.reads_in)
    );
    eprintln!(
        "    dropped:                 {:>12}  ({})",
        s.reads_dropped, pct(s.reads_dropped, s.reads_in)
    );
    let causes: [(&str, u64); 8] = [
        ("polyX 5′",       s.reads_dropped_poly5),
        ("polyX 3′",       s.reads_dropped_poly3),
        ("polyX interior", s.reads_dropped_poly_int),
        ("hard trim",      s.reads_dropped_hard),
        ("mean quality",   s.reads_dropped_mean_q),
        ("adapter",        s.reads_dropped_adapter),
        ("duplicate",      s.reads_dropped_duplicate),
        ("below --min-len", s.reads_dropped_length),
    ];
    let last_nonzero = causes.iter().rposition(|(_, n)| *n > 0);
    for (i, (lbl, count)) in causes.iter().enumerate() {
        if *count == 0 { continue; }
        let prefix = if Some(i) == last_nonzero { "└─" } else { "├─" };
        eprintln!("      {} {:<18} {:>12}  ({})", prefix, lbl, count, pct(*count, s.reads_in));
    }

    if s.adapter.reads_trimmed > 0 {
        eprintln!("    adapter trimmed:         {:>12}  ({})",
            s.adapter.reads_trimmed, pct(s.adapter.reads_trimmed, s.reads_in));
        eprintln!("      bases trimmed (5′):    {:>12}", s.adapter.bases_trimmed_5);
        eprintln!("      bases trimmed (3′):    {:>12}", s.adapter.bases_trimmed_3);
        // Per-adapter breakdown.
        let by_name_map = s.adapter.by_name();
        let mut by_name: Vec<(&String, &u64)> = by_name_map.iter().collect();
        by_name.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
        let last = by_name.len().saturating_sub(1);
        for (i, (name, count)) in by_name.iter().enumerate() {
            let prefix = if i == last { "└─" } else { "├─" };
            eprintln!("      {} {:<18} {:>12}  ({})", prefix, name, count, pct(**count, s.reads_in));
        }
    }

    // D-6: nothing is guessed silently. The offset assumption is the one that
    // turns a +64 file into a 31-point quality error, so it is stated here as
    // well as in the JSON.
    for a in &s.assumptions {
        eprintln!(
            "    assumption — {} = {}: {}; {}",
            a.what, a.value, a.basis, a.consequence
        );
    }

    eprintln!("    bases in:                {:>12}", s.bases_in);
    eprintln!("    bases trimmed (5′):      {:>12}", s.bases_trimmed_5);
    eprintln!("    bases trimmed (3′):      {:>12}", s.bases_trimmed_3);
    eprintln!("    interior poly mass (budget contrib):{:>12}", s.bases_interior_mass);
    let tot = s.bases_trimmed_total();
    eprintln!("    bases trimmed (total):   {:>12}  ({})", tot, pct(tot, s.bases_in));
}

pub fn report_runtime(elapsed: Duration, total_reads: u64) {
    let secs = elapsed.as_secs_f64();
    eprintln!(
        "  elapsed: {:.2}s  ({} reads/s)",
        secs,
        if secs > 0.0 { format!("{:.0}", total_reads as f64 / secs) } else { "—".to_string() }
    );
    cite();
}

/// The one-line citation request, printed at the end of every run (not under -q).
/// fastprep belongs to the Mandrillus Project (DECISIONS D-9).
pub fn cite() {
    eprintln!(
        "If you use fastprep in published research, please cite: \
         Ribes R., Mandier C., Baniel A. FastSuite. https://github.com/RaphaelRibes/FastSuite \
         (fastprep and fastdedup: Mandrillus Project; bioRxiv 2026, doi:10.64898/2026.04.29.721745)"
    );
}

pub fn pct(num: u64, den: u64) -> String {
    if den == 0 { "n/a".into() } else { format!("{:.2}%", 100.0 * num as f64 / den as f64) }
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::{DropReason, TrimResult};

    #[test]
    fn record_emit_and_drop() {
        let mut s = Stats::default();
        s.record(150, &TrimResult::Emit { t5_adapter: 0, t3_adapter: 0, t5_poly: 10, t3_poly: 2, runs: vec![] }, None);
        s.record(150, &TrimResult::Drop { reason: DropReason::Poly5, t5_adapter: 0, t3_adapter: 0, t5_poly: 0, t3_poly: 0, runs: vec![] }, None);
        s.record(100, &TrimResult::Emit { t5_adapter: 0, t3_adapter: 0, t5_poly: 0, t3_poly: 0, runs: vec![(10, 5)] }, None);
        assert_eq!(s.reads_in, 3);
        assert_eq!(s.reads_kept, 2);
        assert_eq!(s.reads_dropped, 1);
        assert_eq!(s.reads_dropped_poly5, 1);
        assert_eq!(s.bases_in, 400);
        assert_eq!(s.bases_trimmed_5, 10);
        assert_eq!(s.bases_trimmed_3, 2);
        assert_eq!(s.bases_interior_mass, 5);
        assert_eq!(s.bases_trimmed_total(), 12);
    }

    #[test]
    fn record_adapter_trim_updates_adapter_stats() {
        let mut s = Stats::default();
        s.adapter.set_names(vec!["TruSeq/R1".to_string()]);
        s.record(
            150,
            &TrimResult::Emit { t5_adapter: 5, t3_adapter: 10, t5_poly: 0, t3_poly: 0, runs: vec![] },
            Some(0),
        );
        assert_eq!(s.adapter.reads_trimmed, 1);
        assert_eq!(s.adapter.bases_trimmed_5, 5);
        assert_eq!(s.adapter.bases_trimmed_3, 10);
        assert_eq!(*s.adapter.by_name().get("TruSeq/R1").unwrap(), 1);
    }

    /// Counts are by index and names are materialised only at report time —
    /// and the serialised `--json` block keeps the shape and the values it
    /// always had, including the reserved `pe-overlap` slot.
    #[test]
    fn adapter_names_are_materialised_at_report_time() {
        let mut s = Stats::default();
        s.adapter.set_names(vec!["TruSeq/R1".to_string(), "Nextera".to_string()]);
        let emit = |t5: usize, t3: usize| TrimResult::Emit {
            t5_adapter: t5,
            t3_adapter: t3,
            t5_poly: 0,
            t3_poly: 0,
            runs: vec![],
        };
        s.record(150, &emit(5, 10), Some(0));
        s.record(150, &emit(2, 0), Some(1));
        // A boundary the mate pair pinned, and an untrimmed read.
        s.record(150, &emit(0, 7), Some(crate::adapter::PE_OVERLAP_ADAPTER));
        s.record(150, &emit(0, 0), Some(0));

        assert_eq!(s.adapter.reads_trimmed, 3);
        assert_eq!(s.adapter.bases_trimmed_5, 7);
        assert_eq!(s.adapter.bases_trimmed_3, 17);
        assert_eq!(s.adapter.counts[0], 1);
        assert_eq!(s.adapter.counts[1], 1);
        assert_eq!(s.adapter.counts[crate::adapter::PE_OVERLAP_ADAPTER as usize], 1);

        let v = serde_json::to_value(&s.adapter).unwrap();
        assert_eq!(v["by_name"]["TruSeq/R1"], 1);
        assert_eq!(v["by_name"]["Nextera"], 1);
        assert_eq!(v["by_name"]["pe-overlap"], 1);
        assert_eq!(v["reads_trimmed"], 3);
        assert_eq!(v["bases_trimmed_5"], 7);
        assert_eq!(v["bases_trimmed_3"], 17);
    }
}
