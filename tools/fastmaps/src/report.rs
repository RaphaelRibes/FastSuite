//! The `fastsuite.report.v1` wrapper and the corner metrics it carries.

use std::collections::BTreeMap;

/// Fraction of reads in each heatmap corner, per mate window.
///
/// Reproduces the convention fasthelper historically derived from
/// `heatmap.json` (tools/fasthelper/src/report.rs): a (T+C, T+G) grid of
/// side `size`, corner window `c = max(size/50, 2)`:
///   low_x/low_y  → polyA corner, high_x/high_y → polyT,
///   low_x/high_y → polyG,      high_x/low_y  → polyC,
/// and off-diagonal mass where `|x - y| > size/10`.
///
/// `grids` is `[("all", &counts, no_call), ("r1", &counts_r1, no_call_r1),
/// ("r2", &counts_r2, no_call_r2)]` with flat row-major index `y*size + x`.
/// Keys use the "" / "r1_" / "r2_" prefixes, exactly as fasthelper's adapter
/// did, so consumers see no difference in metric names.
///
/// `no_call` is the number of reads in that window with zero ACGT bases. Such
/// a read encodes to `(0,0)` — the polyA corner — but it has no measurable
/// composition (an all-N read is not polyA-rich), so it is excluded from every
/// fraction's denominator and reported as `<pre>no_call_reads` instead. The
/// grid itself is untouched: the pixel keeps the read.
pub fn corner_metrics(size: u32, grids: &[(&str, &[u64], u64)]) -> BTreeMap<String, f64> {
    let c = corner_window(size) as usize;
    let diag_lim = size / 10;
    let mut out = BTreeMap::new();
    for &(label, counts, no_call) in grids {
        let (mut poly_a, mut poly_t, mut poly_g, mut poly_c, mut off_diag) =
            (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut total = 0.0f64;
        for y in 0..size as usize {
            for x in 0..size as usize {
                let n = counts[y * size as usize + x] as f64;
                if n == 0.0 {
                    continue;
                }
                total += n;
                // High corner bound: `size - 1 - c`, saturating because a grid
                // no wider than the window (`size <= 1 + c`, i.e. reads of a
                // couple of bases) would otherwise wrap in release and panic in
                // debug — a 1 bp all-T read reported polyA instead of polyT
                // (review 2026-10-21, item 2).
                let high = (size as usize).saturating_sub(1 + c);
                let (low_x, high_x) = (x <= c, x >= high);
                let (low_y, high_y) = (y <= c, y >= high);
                if low_x && low_y { poly_a += n; }
                if high_x && high_y { poly_t += n; }
                if low_x && high_y { poly_g += n; }
                if high_x && low_y { poly_c += n; }
                if (x as i64 - y as i64).abs() as u32 > diag_lim { off_diag += n; }
            }
        }
        if total == 0.0 {
            continue;
        }
        // No-call reads share the (0,0) pixel with real polyA reads: take them
        // out of that corner and out of the denominator, so an all-N read can
        // neither raise polyA_frac nor dilute the other corners. `total` keeps
        // every read, so total_reads stays the window's true read count.
        let denom = total - no_call as f64;
        let poly_a = (poly_a - no_call as f64).max(0.0);
        let p = |n: f64| if denom > 0.0 { n / denom } else { 0.0 };
        let pre = if label.is_empty() { String::new() } else { format!("{label}_") };
        out.insert(format!("{pre}polyA_frac"), p(poly_a));
        out.insert(format!("{pre}polyT_frac"), p(poly_t));
        out.insert(format!("{pre}polyG_frac"), p(poly_g));
        out.insert(format!("{pre}polyC_frac"), p(poly_c));
        out.insert(format!("{pre}off_diagonal_frac"), p(off_diag));
        out.insert(format!("{pre}no_call_reads"), no_call as f64);
        if label.is_empty() {
            out.insert("total_reads".into(), total);
        }
    }
    out
}

/// Side of the corner window on a grid of side `size`: `max(size/50, 2)`.
///
/// Published as `corner_window` in both outputs, and derived here only, because
/// it is the number that decides whether a read can reach a corner at all: the
/// grid is sized for the *longest* read, so a read shorter than `size - 1 - c`
/// never gets there.
pub fn corner_window(size: u32) -> u32 {
    (size / 50).max(2)
}

// ── 3' polyX tails ───────────────────────────────────────────────────────────
//
// The corner fractions describe *whole-read* composition, so a read with a
// random 131 bp body and a 20 bp polyG tail encodes to roughly (65, 85) on a
// 2048-wide plane — mid-plane, no corner, `polyG_frac` 0. The tail is real
// contamination and the corners cannot see it (M3). The detector is the same
// corner test applied to the last `POLY_TAIL_WINDOW` bases only: an all-G
// window is the polyG corner of a local plane.

/// Length of the 3' window the terminal-homopolymer detector inspects.
pub const POLY_TAIL_WINDOW: usize = 20;

/// Shortest 3' window the detector will judge. A read shorter than this has
/// too few bases for "the read ends in a run of one base" to mean anything.
pub const POLY_TAIL_MIN_WINDOW: usize = 8;

/// Which base (index into `[A, C, G, T]`) the read ends in a run of, if any.
///
/// `None` for a read that does not end in a homopolymer, for a read shorter
/// than `POLY_TAIL_MIN_WINDOW`, and for a run of non-ACGT bytes (an N-tail has
/// no measurable composition either).
pub fn tail_base(seq: &[u8]) -> Option<usize> {
    let n = seq.len().min(POLY_TAIL_WINDOW);
    if n < POLY_TAIL_MIN_WINDOW {
        return None;
    }
    let win = &seq[seq.len() - n..];
    let base = match win[0] {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return None,
    };
    if win.iter().all(|&b| b == win[0]) {
        Some(base)
    } else {
        None
    }
}

/// Notes for 3'-end-concentrated libraries the corner fractions cannot see.
///
/// `tail` counts, for each of `[A, C, G, T]`, the reads that end in a run of at
/// least `POLY_TAIL_WINDOW` copies of that base. A note is emitted when such
/// reads exist and the matching whole-library corner fraction is 0 — the exact
/// blind spot M3 reports: `r2_polyG_frac 0.0` on 200 reads that all end in
/// 20 G. The corner definition is not changed; the run says what it cannot see.
pub fn poly_tail_notes(
    tail: &[u64; 4],
    metrics: &BTreeMap<String, f64>,
    total_reads: u64,
) -> Vec<String> {
    let mut notes = Vec::new();
    for (base, name) in [(0usize, "A"), (1, "C"), (2, "G"), (3, "T")] {
        let n = tail[base];
        if n == 0 {
            continue;
        }
        let key = format!("poly{name}_frac");
        let frac = metrics.get(&key).copied().unwrap_or(0.0);
        if frac > 0.0 {
            continue;
        }
        let total = total_reads.max(1);
        notes.push(format!(
            "3' poly{name} tail: {n} of {total} reads ({pct:.1}%) end in a run of at least \
             {POLY_TAIL_WINDOW} identical {name} bases, but {key} is 0 — the corner fractions \
             measure whole-read composition, so a 3' tail on a random body lands mid-plane and \
             never reaches a corner. The tail is contamination the corners cannot see: trim it \
             upstream (fastprep) or inspect 3' ends directly.",
            pct = 100.0 * n as f64 / total as f64,
        ));
    }
    notes
}

/// Write the `fastsuite.report.v1` wrapper (docs/REPORT-SCHEMA.md) to `path`.
///
/// This *is* the `--json PATH` file; the grid JSON keeps its own shape and its
/// own home, `<output-dir>/heatmap.json`. `grid_name` is that file's name, so
/// the `grid` block names a file the consumer can actually find (it used to
/// hard-code `heatmap.json` while `--json` was redirecting the grid).
///
/// `notes` carries what the metrics cannot: a run-level caveat (a 3' polyX
/// tail the corner fractions are blind to, M3) is a string, and `metrics` must
/// stay a flat numeric namespace. Always emitted, empty or not, so a consumer
/// can tell "no notes" from "notes not supported".
pub fn write_report_v1(
    path: &std::path::Path,
    command: &str,
    elapsed_s: f64,
    files: &[String],
    metrics: &BTreeMap<String, f64>,
    notes: &[String],
    grid_name: &str,
) -> std::io::Result<()> {
    use crate::encode::json_escape;

    let metrics_block = metrics_json(metrics);
    let files_json = files
        .iter()
        .map(|f| json_escape(f))
        .collect::<Vec<_>>()
        .join(", ");
    let notes_json = notes
        .iter()
        .map(|n| json_escape(n))
        .collect::<Vec<_>>()
        .join(", ");

    let doc = format!(
        concat!(
            "{{\n",
            "  \"schema\": \"fastsuite.report.v1\",\n",
            "  \"tool\": \"fastmaps\",\n",
            "  \"tool_version\": \"{}\",\n",
            "  \"command\": {},\n",
            "  \"mode\": \"composition\",\n",
            "  \"inputs\": {{\"files\": [{files}]}},\n",
            "  \"elapsed_s\": {},\n",
            "  \"metrics\": {{\n{metrics}\n  }},\n",
            "  \"notes\": [{notes}],\n",
            "  \"grid\": {{\"format\": {}, \"see\": \"docs/PLAN-QC.md\"}}\n",
            "}}\n"
        ),
        env!("CARGO_PKG_VERSION"),
        json_escape(command),
        fmt_f64(elapsed_s),
        grid = json_escape(grid_name),
        files = files_json,
        metrics = metrics_block,
        notes = notes_json,
    );

    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    use std::io::Write;
    f.write_all(doc.as_bytes())?;
    f.flush()
}

/// Render a metric map as the flat numeric JSON object both shapes carry.
///
/// One definition for both files. `fastsuite.report.v1` puts it under
/// `metrics`, and `<output-dir>/heatmap.json` carries the same object, because
/// the two described the same run and did not agree: the grid had no
/// `no_call_reads`, so a consumer that recomputed the corner fractions from
/// `counts` — fasthelper's documented `concatenate` adapter does that — read
/// `polyA_frac 1.000` for a run of all-N reads where the report said 0. A
/// no-call read has no measurable composition but still lands on the (0,0)
/// pixel, so only this number can take it back out of the polyA corner.
pub fn metrics_json(metrics: &BTreeMap<String, f64>) -> String {
    metrics
        .iter()
        .map(|(k, v)| format!("    {}: {}", crate::encode::json_escape(k), fmt_f64(*v)))
        .collect::<Vec<_>>()
        .join(",\n")
}

/// Compact float rendering: integers print without a `.0`, everything else
/// with 6 significant decimals. Keeps the hand-written JSON small.
fn fmt_f64(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        let s = format!("{x:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// size=11, c=2 fixture: 4 reads of length 10.
    /// R1 (file A): TTTTTTTTTT (x=10,y=10) and GGGGGGGGGG (x=0,y=10);
    /// R2 (file B): GGGGGGGGGG (x=0,y=10) and AAAAAAAAAA (x=0,y=0).
    /// Matches the verified fastmaps numbers in the convergence plan:
    /// all: polyA .25, polyT .25, polyG .5, polyC 0, off_diag .5, total 4;
    /// r1: polyT .5, polyG .5; r2: polyA .5, polyG .5.
    fn grids() -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        let size = 11;
        let mut all = vec![0u64; size * size];
        let mut r1 = vec![0u64; size * size];
        let mut r2 = vec![0u64; size * size];
        all[10 * size + 10] += 1; // T at (10,10), from file A's R1
        all[10 * size] += 2;      // G at (0,10), once per file
        r1[10 * size + 10] += 1;
        r1[10 * size] += 1;
        r2[10 * size] += 1;
        all[0] += 1;              // A at (0,0), from file B's R2
        r2[0] += 1;
        (all, r1, r2)
    }

    #[test]
    fn corner_metrics_matches_fasthelper_fixture() {
        let (all, r1, r2) = grids();
        let m = corner_metrics(11, &[("", &all, 0), ("r1", &r1, 0), ("r2", &r2, 0)]);
        let near = |k: &str, want: f64| {
            assert!((m[k] - want).abs() < 1e-9, "{k}: got {}, want {want}", m[k]);
        };
        near("polyA_frac", 0.25);
        near("polyT_frac", 0.25);
        near("polyG_frac", 0.5);
        near("polyC_frac", 0.0);
        near("off_diagonal_frac", 0.5);
        near("total_reads", 4.0);
        near("r1_polyT_frac", 0.5);
        near("r1_polyG_frac", 0.5);
        near("r1_polyA_frac", 0.0);
        near("r2_polyA_frac", 0.5);
        near("r2_polyT_frac", 0.0);
        near("r2_polyG_frac", 0.5);
    }

    #[test]
    fn corner_metrics_empty_grid_is_omitted() {
        // The two empty windows contribute no metrics at all; the populated
        // one still yields the full set. A banded grid keeps all corners and
        // the off-diagonal region empty — the opposite of the plan's original
        // "uniform grid" premise, which is corner-enriched by construction
        // (one pixel per corner vs. nine per interior pixel).
        let size = 11;
        let banded = || -> Vec<u64> {
            let mut v = vec![0u64; size * size];
            for y in 0..size {
                for x in 0..size {
                    if (x as i64 - y as i64).abs() <= 1 {
                        v[y * size + x] = 1;
                    }
                }
            }
            v
        };
        let all = vec![0u64; size * size];
        let r1 = vec![0u64; size * size];
        let r2 = banded();
        let m = corner_metrics(size as u32, &[("", &all, 0), ("r1", &r1, 0), ("r2", &r2, 0)]);
        assert!(!m.contains_key("polyT_frac"), "empty window must be skipped");
        assert!(!m.contains_key("r1_polyT_frac"));
        assert!(!m.contains_key("total_reads"), "total_reads lives in the all window");
        assert!(m.contains_key("r2_polyT_frac"));
        // The two high-y/low-x corners (polyG, polyC) and the off-diagonal
        // band are empty for a banded grid. The polyA corner holds 7 banded
        // pixels (a 3x3 window with the two |x-y|=2 pixels absent) out of 31
        // total; polyT holds 7 likewise.
        assert!((m["r2_polyG_frac"] - 0.0).abs() < 1e-9);
        assert!((m["r2_polyC_frac"] - 0.0).abs() < 1e-9);
        assert!((m["r2_off_diagonal_frac"] - 0.0).abs() < 1e-9);
        assert!((m["r2_polyT_frac"] - 7.0 / 31.0).abs() < 1e-9);
        assert!((m["r2_polyA_frac"] - 7.0 / 31.0).abs() < 1e-9);
    }

    #[test]
    fn report_v1_parses_with_expected_values() {
        let dir = std::env::temp_dir().join(format!("maps_v1_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("r.json");
        let mut metrics = std::collections::BTreeMap::new();
        metrics.insert("total_reads".into(), 4.0);
        metrics.insert("polyT_frac".into(), 0.25);
        metrics.insert("r2_polyG_frac".into(), 0.5);
        super::write_report_v1(
            &path,
            "fastmaps a_R1.fq.gz a_R2.fastq.gz -j 2",
            0.42,
            &["a_R1.fq.gz".into(), "a_R2.fastq.gz".into()],
            &metrics,
            &[],
            "heatmap.json",
        )
        .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();

        // Actually parse the emitted document: a stray brace, a trailing comma
        // or a duplicate key would all make this fail, which the old
        // substring checks could not catch.
        let v: serde_json::Value = serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("report JSON must parse: {e}\n---\n{text}"));

        assert_eq!(v["schema"], "fastsuite.report.v1");
        assert_eq!(v["tool"], "fastmaps");
        assert_eq!(v["tool_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(v["mode"], "composition");
        assert_eq!(v["command"], "fastmaps a_R1.fq.gz a_R2.fastq.gz -j 2");
        assert!((v["elapsed_s"].as_f64().unwrap() - 0.42).abs() < 1e-12);
        assert_eq!(v["inputs"]["files"], serde_json::json!(["a_R1.fq.gz", "a_R2.fastq.gz"]));
        assert_eq!(v["metrics"]["total_reads"].as_f64().unwrap(), 4.0);
        assert_eq!(v["metrics"]["polyT_frac"].as_f64().unwrap(), 0.25);
        assert_eq!(v["metrics"]["r2_polyG_frac"].as_f64().unwrap(), 0.5);
        assert_eq!(v["grid"]["format"], "heatmap.json");

        // Still-meaningful substring contract: fmt_f64 renders integral
        // floats compactly, so total_reads is `4`, never `4.0`.
        assert!(text.contains("\"total_reads\": 4"), "{text}");
        assert!(!text.contains("\"total_reads\": 4.0"), "{text}");

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The `grid` block names the grid file the writer was handed, not a
    /// hard-coded `heatmap.json`.
    #[test]
    fn report_names_the_grid_file_it_is_given() {
        let dir = std::env::temp_dir().join(format!("maps_grid_name_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("r.json");
        let mut metrics = std::collections::BTreeMap::new();
        metrics.insert("total_reads".into(), 1.0);
        super::write_report_v1(&path, "fastmaps a.fq", 0.0, &["a.fq".into()], &metrics, &[], "corner.json")
            .unwrap();
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(v["grid"]["format"], "corner.json", "{v}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A note is a string, so it goes in `notes` — `metrics` stays numeric.
    #[test]
    fn notes_are_emitted_as_a_string_array() {
        let dir = std::env::temp_dir().join(format!("maps_notes_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("r.json");
        let mut metrics = std::collections::BTreeMap::new();
        metrics.insert("total_reads".into(), 1.0);
        let notes = vec!["3' polyG tail: 1 of 1 reads ends in G".to_string()];
        super::write_report_v1(&path, "fastmaps a.fq", 0.0, &["a.fq".into()], &metrics, &notes, "heatmap.json")
            .unwrap();
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(v["notes"][0], notes[0].as_str(), "{v}");
        assert!(
            v["metrics"]["total_reads"].is_number(),
            "metrics must stay a flat numeric namespace: {v}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ── M3: the 3' tail detector ────────────────────────────────────────────

    #[test]
    fn tail_base_finds_a_terminal_homopolymer() {
        let mut seq = b"ACGTACGTACGTACGTACGT".to_vec(); // random-ish body
        seq.extend(std::iter::repeat_n(b'G', POLY_TAIL_WINDOW));
        assert_eq!(super::tail_base(&seq), Some(2), "G");

        let mut seq = b"ACGTTGCAACGTTGCAACGT".to_vec();
        seq.extend(std::iter::repeat_n(b'T', POLY_TAIL_WINDOW));
        assert_eq!(super::tail_base(&seq), Some(3), "T");
    }

    #[test]
    fn tail_base_ignores_reads_without_a_terminal_run() {
        // Random 3' end: no run.
        assert_eq!(super::tail_base(b"ACGTACGTACGTACGTACGTACGTACGTACGT"), None);
        // A run that does not reach the 3' end.
        assert_eq!(super::tail_base(b"GGGGGGGGGGGGGGGGGGGGACGT"), None);
        // Too short to judge.
        assert_eq!(super::tail_base(b"GGGGGGG"), None);
        // A run of N: no measurable composition.
        assert_eq!(super::tail_base(b"ACGTACGTACGTACGTACGTNNNNNNNNNNNNNNNNNNNN"), None);
        // A read no longer than the minimum window, entirely one base, is
        // judged: the whole read is the 3' end.
        assert_eq!(super::tail_base(b"GGGGGGGG"), Some(2));
    }

    #[test]
    fn tail_notes_fire_only_when_the_corner_is_blind() {
        let tail = [0u64, 0, 200, 0];
        let mut metrics = std::collections::BTreeMap::new();
        metrics.insert("polyG_frac".into(), 0.0);
        let notes = super::poly_tail_notes(&tail, &metrics, 200);
        assert_eq!(notes.len(), 1, "{notes:?}");
        assert!(notes[0].contains("polyG") && notes[0].contains("200") && notes[0].contains("3'"));

        // The corner already sees them: nothing to add.
        metrics.insert("polyG_frac".into(), 0.5);
        assert!(super::poly_tail_notes(&tail, &metrics, 200).is_empty());

        // No tails at all.
        assert!(super::poly_tail_notes(&[0; 4], &metrics, 200).is_empty());
    }
}
