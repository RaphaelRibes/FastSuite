//! Pipeline statistics.
//!
//! Accumulated per-chunk in the parallel path (no shared mutex), then
//! merged once per chunk on the main thread. At program end, rendered to
//! stderr.

use std::collections::HashMap;
use std::time::Duration;

use crate::algorithm::{DropReason, TrimResult};

// ---------------------------------------------------------------------------
// Per-adapter stats (counts by adapter name).
// ---------------------------------------------------------------------------

#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct AdapterStats {
    /// reads → adapter name → count of reads where that adapter was trimmed.
    pub by_name: HashMap<String, u64>,
    pub reads_trimmed: u64,
    pub bases_trimmed_5: u64,
    pub bases_trimmed_3: u64,
}

impl AdapterStats {
    pub fn record(&mut self, name: &str, t5: usize, t3: usize) {
        *self.by_name.entry(name.to_string()).or_insert(0) += 1;
        self.reads_trimmed += 1;
        self.bases_trimmed_5 += t5 as u64;
        self.bases_trimmed_3 += t3 as u64;
    }

    pub fn merge(&mut self, other: &AdapterStats) {
        for (name, count) in &other.by_name {
            *self.by_name.entry(name.clone()).or_insert(0) += count;
        }
        self.reads_trimmed += other.reads_trimmed;
        self.bases_trimmed_5 += other.bases_trimmed_5;
        self.bases_trimmed_3 += other.bases_trimmed_3;
    }
}

// ---------------------------------------------------------------------------
// Main Stats struct.
// ---------------------------------------------------------------------------

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
    /// Adapter-specific stats; separate to isolate HashMap allocation.
    pub adapter: AdapterStats,
}

impl Stats {
    pub fn record(&mut self, original_len: usize, res: &TrimResult, adapter_name: Option<&str>) {
        self.reads_in += 1;
        self.bases_in += original_len as u64;
        match res {
            TrimResult::Drop { reason, t5_adapter, t3_adapter, t5_poly, t3_poly, .. } => {
                // Match ergonomics: res is &TrimResult, so fields are &T.
                let (t5a, t3a) = (*t5_adapter, *t3_adapter);
                let _ = (*t5_poly, *t3_poly);
                self.reads_dropped += 1;
                match reason {
                    DropReason::Poly5    => self.reads_dropped_poly5 += 1,
                    DropReason::Poly3    => self.reads_dropped_poly3 += 1,
                    DropReason::PolyInt  => self.reads_dropped_poly_int += 1,
                    DropReason::Hard     => self.reads_dropped_hard += 1,
                    DropReason::MeanQ    => self.reads_dropped_mean_q += 1,
                    DropReason::Adapter  => self.reads_dropped_adapter += 1,
                    DropReason::Duplicate => self.reads_dropped_duplicate += 1,
                }
                if (t5a + t3a) > 0 {
                    if let Some(name) = adapter_name {
                        self.adapter.record(name, t5a, t3a);
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
                    if let Some(name) = adapter_name {
                        self.adapter.record(name, t5a, t3a);
                    }
                }
            }
        }
    }

    pub fn merge(&mut self, other: &Stats) {
        self.reads_in             += other.reads_in;
        self.reads_kept           += other.reads_kept;
        self.reads_dropped        += other.reads_dropped;
        self.reads_dropped_poly5  += other.reads_dropped_poly5;
        self.reads_dropped_poly3  += other.reads_dropped_poly3;
        self.reads_dropped_poly_int += other.reads_dropped_poly_int;
        self.reads_dropped_hard   += other.reads_dropped_hard;
        self.reads_dropped_mean_q += other.reads_dropped_mean_q;
        self.reads_dropped_adapter += other.reads_dropped_adapter;
        self.reads_dropped_duplicate += other.reads_dropped_duplicate;
        self.bases_in             += other.bases_in;
        self.bases_trimmed_5      += other.bases_trimmed_5;
        self.bases_trimmed_3      += other.bases_trimmed_3;
        self.bases_interior_mass += other.bases_interior_mass;
        self.adapter.merge(&other.adapter);
    }

    pub fn bases_trimmed_total(&self) -> u64 {
        self.bases_trimmed_5 + self.bases_trimmed_3
    }
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
    let causes: [(&str, u64); 7] = [
        ("polyX 5′",       s.reads_dropped_poly5),
        ("polyX 3′",       s.reads_dropped_poly3),
        ("polyX interior", s.reads_dropped_poly_int),
        ("hard trim",      s.reads_dropped_hard),
        ("mean quality",   s.reads_dropped_mean_q),
        ("adapter",        s.reads_dropped_adapter),
        ("duplicate",      s.reads_dropped_duplicate),
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
        let mut by_name: Vec<(&String, &u64)> = s.adapter.by_name.iter().collect();
        by_name.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
        let last = by_name.len().saturating_sub(1);
        for (i, (name, count)) in by_name.iter().enumerate() {
            let prefix = if i == last { "└─" } else { "├─" };
            eprintln!("      {} {:<18} {:>12}  ({})", prefix, name, count, pct(**count, s.reads_in));
        }
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
        s.record(
            150,
            &TrimResult::Emit { t5_adapter: 5, t3_adapter: 10, t5_poly: 0, t3_poly: 0, runs: vec![] },
            Some("TruSeq/R1"),
        );
        assert_eq!(s.adapter.reads_trimmed, 1);
        assert_eq!(s.adapter.bases_trimmed_5, 5);
        assert_eq!(s.adapter.bases_trimmed_3, 10);
        assert_eq!(*s.adapter.by_name.get("TruSeq/R1").unwrap(), 1);
    }

    #[test]
    fn merge_adds_fields() {
        let mut a = Stats::default();
        a.record(100, &TrimResult::Emit { t5_adapter: 0, t3_adapter: 0, t5_poly: 5, t3_poly: 5, runs: vec![] }, None);
        let mut b = Stats::default();
        b.record(100, &TrimResult::Drop { reason: DropReason::Poly3, t5_adapter: 0, t3_adapter: 0, t5_poly: 0, t3_poly: 0, runs: vec![] }, None);
        a.merge(&b);
        assert_eq!(a.reads_in, 2);
        assert_eq!(a.reads_kept, 1);
        assert_eq!(a.reads_dropped, 1);
        assert_eq!(a.reads_dropped_poly3, 1);
    }

    #[test]
    fn merge_adapter_stats_accumulates_by_name() {
        let mut a = Stats::default();
        a.record(100, &TrimResult::Emit { t5_adapter: 0, t3_adapter: 5, t5_poly: 0, t3_poly: 0, runs: vec![] }, Some("Nextera"));
        let mut b = Stats::default();
        b.record(100, &TrimResult::Emit { t5_adapter: 0, t3_adapter: 8, t5_poly: 0, t3_poly: 0, runs: vec![] }, Some("Nextera"));
        a.merge(&b);
        assert_eq!(*a.adapter.by_name.get("Nextera").unwrap(), 2);
        assert_eq!(a.adapter.bases_trimmed_3, 13);
    }
}
