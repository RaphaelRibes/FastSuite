//! Report model, adapters and writers.
//!
//! The suite writes one JSON per tool per sample. fasthelper normalises all of
//! them into a single metric table so a run can be summarised in one file.
//!
//! Two shapes are understood, deliberately:
//!   * `fastsuite.report.v1` — the canonical wrapper (schema/tool/metrics)
//!   * legacy raw JSONs — fastmaps' heatmap.json grid and pre-wrapper
//!     fastprep stats objects, so FastHelper works with binaries (and
//!     report files) that predate the wrapper.

use anyhow::{Context, Result};
use serde::de::{Deserialize, Deserializer, IgnoredAny, MapAccess, SeqAccess, Visitor};
use serde::Serialize;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::ReportArgs;

pub const SCHEMA_VERSION: &str = "fastsuite.report.v1";

pub const SCHEMA_DOC: &str = r#"fastsuite.report.v1 — per-tool report contract
=======================================================================
Every FastSuite tool writes exactly one JSON per run with --json <PATH>:

{
  "schema":       "fastsuite.report.v1",
  "tool":         "fastdedup",              // binary name
  "tool_version": "1.3.0",                  // CARGO_PKG_VERSION
  "command":      "fastdedup -1 a.fq ...",  // exact argv, space-joined
  "mode":         "paired",                 // optional, tool-defined
  "dry_run":      false,
  "inputs":       {"r1": "...", "r2": "..."},
  "outputs":      {"r1": "...", "r2": "..."},
  "elapsed_s":    12.5,                     // wall clock, seconds
  "metrics": {                              // flat numeric namespace
    "reads_in": 1000000,
    "reads_out": 998000,
    "reads_dropped": 2000
  }
}

Rules
-----
1. `metrics` values are numbers only (booleans as 0/1). Nested structures
   belong in tool-specific blocks outside `metrics`.
2. Metric names are snake_case and tool-prefixed only when ambiguous
   (`reads_in`, not `fastdedup_reads_in`).
3. Recommended shared keys, honoured by every tool that can produce them:
   reads_in, reads_out, reads_dropped, bases_in, bases_out, bases_trimmed,
   elapsed_s is top-level, none of the metric keys are mandatory.
4. Adding a metric is a MINOR change; renaming/removing one is a MAJOR change
   (bump the schema string, keep the old key one release as an alias).
5. FastHelper also accepts legacy raw JSON (pre-wrapper fastprep stats
   objects, fastmaps heatmap.json) and derives the keys above; every new
   tool emits v1.

Producers today
---------------
  fastdedup   --json  -> v1 (native)
  fastprep    --json  -> v1 (native)
  fastmaps    --json  -> heatmap.json in --output-dir (v1 wrapper planned)
  fastsim     -       -> no JSON yet (planned: --json with truth summary)
"#;

/// One tool's report for one sample, normalised.
#[derive(Debug, Clone)]
pub struct ToolReport {
    pub sample: String,
    pub tool: String,
    pub version: Option<String>,
    pub source: PathBuf,
    pub metrics: BTreeMap<String, f64>,
    pub notes: Vec<String>,
    /// True when this came from a `fastsuite.report.v1` wrapper. A v1 report
    /// wins over a legacy/raw file describing the same (sample, tool): the
    /// wrapper is the tool's own statement, the raw file a derivation of it.
    pub from_v1: bool,
}

fn num(v: &serde_json::Value) -> Option<f64> {
    match v {
        serde_json::Value::Number(n) => n.as_f64(),
        serde_json::Value::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
        _ => None,
    }
}

fn get_num(o: &serde_json::Value, key: &str) -> f64 {
    o.get(key).and_then(num).unwrap_or(0.0)
}

/// Sample name from a path: file stem, unless that stem is a generic report
/// name — then the parent directory carries the sample identity.
///
/// A `fastsuite.report.v1` wrapper written *beside* its grid
/// (`<stem>.report.json`, what `fastmaps --json <dir>/heatmap.json` leaves
/// behind) inherits the grid's stem, so `heatmap.report` is tested as
/// `heatmap`: generic, and the directory supplies the name. The strip is
/// gated on the remainder being generic, so a library genuinely called
/// `align.report` keeps its own name.
pub fn sample_from_path(path: &Path) -> String {
    const GENERIC: &[&str] = &[
        "report", "reports", "stats", "json", "out", "output", "heatmap", "fastprep",
        "fastdedup", "fastmaps", "fastsim", "fasthelper", "dedup", "prep", "maps", "result",
        "results", "summary", "manifest", "positions", "position", "pos", "matrix",
    ];
    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
        .to_string();
    let stem = match stem.strip_suffix(".report") {
        Some(base) if GENERIC.contains(&base.to_ascii_lowercase().as_str()) => base.to_string(),
        _ => stem,
    };
    if GENERIC.contains(&stem.to_ascii_lowercase().as_str()) {
        // Walk up until a directory that is not itself a generic report name:
        // runs/demo/<sample>/maps/heatmap.json must yield "<sample>", not "maps".
        let mut cur = path.parent();
        while let Some(dir) = cur {
            match dir.file_name().and_then(|s| s.to_str()) {
                Some(name)
                    if !name.is_empty()
                        && name != "."
                        && !GENERIC.contains(&name.to_ascii_lowercase().as_str()) =>
                {
                    return name.to_string();
                }
                Some(_) => cur = dir.parent(),
                None => break,
            }
        }
    }
    stem
}

// ------------------------------------------------------ fastmaps grid ----
//
// `heatmap.json` is a `size × size` raster of (T+C, T+G) counts plus metadata:
// 140 KB at 151 bp against 3.9 KB for an ordinary report, one per sample. It
// is read through a typed visitor rather than a `Value` tree, because the old
// reader's cost was one enum construction per cell (22 801 per file). The five
// corner sums and the total are accumulated in the same x-major order with the
// same f64 operations as before, so every fraction stays bit-identical.

/// Accumulator for one grid window: the five corner sums and the total.
#[derive(Debug, Default, Clone, Copy)]
struct WindowSums {
    /// polyA, polyT, polyG, polyC, off-diagonal
    corner: [f64; 5],
    total: f64,
    /// cells seen, including any past `size²`
    cells: usize,
    /// cells a complete window has
    needed: usize,
}

impl WindowSums {
    fn new(size: i64) -> Self {
        Self {
            needed: (size * size) as usize,
            ..Self::default()
        }
    }

    /// Fold one cell in, at linear index `i` (x-major: x = i % size,
    /// y = i / size). The `n == 0.0` shortcut and the order of the five tests
    /// are what keep the sums bit-identical to the previous reader.
    fn push(&mut self, size: i64, c: i64, i: usize, n: f64) {
        self.cells += 1;
        if size <= 0 || i >= self.needed || n == 0.0 {
            return;
        }
        let (x, y) = ((i as i64) % size, (i as i64) / size);
        self.total += n;
        let low_x = x <= c;
        let high_x = x >= size - 1 - c;
        let low_y = y <= c;
        let high_y = y >= size - 1 - c;
        if low_x && low_y {
            self.corner[0] += n;
        }
        if high_x && high_y {
            self.corner[1] += n;
        }
        if low_x && high_y {
            self.corner[2] += n;
        }
        if high_x && low_y {
            self.corner[3] += n;
        }
        if (x - y).abs() > size / 10 {
            self.corner[4] += n;
        }
    }

    /// A window shorter than the grid is not a measurement of that grid: the
    /// previous reader skipped it entirely, and so does this one.
    fn complete(&self) -> bool {
        self.cells >= self.needed
    }
}

/// One window: a bare array of counts (`counts`), or an object carrying one
/// under a nested `counts` key. `bare_array` is false for `counts_r1`/
/// `counts_r2`, which fastmaps writes as bare arrays — the previous reader
/// looked for a nested key there and so never used them; that is preserved.
struct WindowVisitor {
    size: i64,
    c: i64,
    bare_array: bool,
}

impl<'de> Visitor<'de> for WindowVisitor {
    type Value = Option<WindowSums>;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a grid of counts (an array, or an object with \"counts\")")
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        if !self.bare_array {
            while seq.next_element::<IgnoredAny>()?.is_some() {}
            return Ok(None);
        }
        let mut sums = WindowSums::new(self.size);
        let mut i = 0usize;
        while let Some(n) = seq.next_element::<f64>()? {
            sums.push(self.size, self.c, i, n);
            i += 1;
        }
        Ok(Some(sums))
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut out = None;
        while let Some(key) = map.next_key::<Cow<'_, str>>()? {
            if key == "counts" {
                out = map.next_value_seed(WindowSeed {
                    size: self.size,
                    c: self.c,
                    bare_array: true,
                })?;
            } else {
                map.next_value::<IgnoredAny>()?;
            }
        }
        Ok(out)
    }

    fn visit_unit<E>(self) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_none<E>(self) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_bool<E>(self, _: bool) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_i64<E>(self, _: i64) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_u64<E>(self, _: u64) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_f64<E>(self, _: f64) -> Result<Self::Value, E> {
        Ok(None)
    }
    fn visit_str<E>(self, _: &str) -> Result<Self::Value, E> {
        Ok(None)
    }
}

struct WindowSeed {
    size: i64,
    c: i64,
    bare_array: bool,
}

impl<'de> serde::de::DeserializeSeed<'de> for WindowSeed {
    type Value = Option<WindowSums>;

    fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_any(WindowVisitor {
            size: self.size,
            c: self.c,
            bare_array: self.bare_array,
        })
    }
}

/// The keys of `heatmap.json` fasthelper uses. Everything else (`files`,
/// `per_file_counts`, `data_row_offsets`, …) is skipped without being built.
struct GridDoc {
    size: Option<i64>,
    total_reads: Option<f64>,
    version: Option<f64>,
    /// windows in metric order: all, r1, r2
    windows: [Option<WindowSums>; 3],
    /// `counts` numbers buffered because `size` had not been seen yet — JSON
    /// object order is not guaranteed, so a grid may list them in any order.
    deferred: [Option<Vec<f64>>; 3],
    has_counts: bool,
}

impl Default for GridDoc {
    fn default() -> Self {
        Self {
            size: None,
            total_reads: None,
            version: None,
            windows: [None, None, None],
            deferred: [None, None, None],
            has_counts: false,
        }
    }
}

fn read_window<'de, A: MapAccess<'de>>(
    map: &mut A,
    doc: &mut GridDoc,
    idx: usize,
    top_level: bool,
) -> Result<(), A::Error> {
    match doc.size.filter(|s| *s > 0) {
        Some(size) => {
            let c = (size / 50).max(2); // corner window: 2% of the read length
            doc.windows[idx] = map.next_value_seed(WindowSeed {
                size,
                c,
                bare_array: top_level,
            })?;
        }
        None if top_level => {
            doc.deferred[idx] = Some(map.next_value::<Vec<f64>>()?);
        }
        None => {
            map.next_value::<IgnoredAny>()?;
        }
    }
    Ok(())
}

impl<'de> Deserialize<'de> for GridDoc {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        d.deserialize_map(GridDocVisitor)
    }
}

struct GridDocVisitor;

impl<'de> Visitor<'de> for GridDocVisitor {
    type Value = GridDoc;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a fastmaps heatmap object")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<GridDoc, A::Error> {
        let mut doc = GridDoc::default();
        while let Some(key) = map.next_key::<Cow<'_, str>>()? {
            match key.as_ref() {
                "size" => doc.size = Some(map.next_value::<f64>()? as i64),
                "total_reads" => doc.total_reads = Some(map.next_value::<f64>()?),
                "version" => doc.version = Some(map.next_value::<f64>()?),
                "counts" => {
                    doc.has_counts = true;
                    read_window(&mut map, &mut doc, 0, true)?;
                }
                "counts_r1" => read_window(&mut map, &mut doc, 1, false)?,
                "counts_r2" => read_window(&mut map, &mut doc, 2, false)?,
                _ => {
                    map.next_value::<IgnoredAny>()?;
                }
            }
        }
        Ok(doc)
    }
}

impl GridDoc {
    /// `counts` and `size` both present: this is a grid, not some other JSON
    /// object that happens to carry a `counts` key.
    fn is_grid(&self) -> bool {
        self.has_counts && self.size.is_some()
    }

    fn into_report(mut self, sample: String, path: &Path) -> ToolReport {
        self.fold_deferred();
        let size = self.size.unwrap_or(0);
        let mut metrics = BTreeMap::new();
        metrics.insert("total_reads".into(), self.total_reads.unwrap_or(0.0));
        metrics.insert("grid_size".into(), size as f64);
        for (idx, label) in ["all", "r1", "r2"].into_iter().enumerate() {
            let Some(w) = self.windows[idx] else { continue };
            if !w.complete() || w.total == 0.0 {
                continue;
            }
            let p = |n: f64| n / w.total;
            let pre = if label == "all" {
                String::new()
            } else {
                format!("{label}_")
            };
            metrics.insert(format!("{pre}polyA_frac"), p(w.corner[0]));
            metrics.insert(format!("{pre}polyT_frac"), p(w.corner[1]));
            metrics.insert(format!("{pre}polyG_frac"), p(w.corner[2]));
            metrics.insert(format!("{pre}polyC_frac"), p(w.corner[3]));
            metrics.insert(format!("{pre}off_diagonal_frac"), p(w.corner[4]));
            if label == "all" {
                metrics.insert("total_reads".into(), w.total);
            }
        }
        ToolReport {
            sample,
            tool: "fastmaps".into(),
            version: self.version.map(|n| format!("grid v{n}")),
            source: path.to_path_buf(),
            metrics,
            notes: Vec::new(),
            from_v1: false,
        }
    }

    /// Fold a `counts` array that arrived before `size` into its window.
    fn fold_deferred(&mut self) {
        let Some(size) = self.size.filter(|s| *s > 0) else {
            return;
        };
        let c = (size / 50).max(2);
        let deferred = std::mem::replace(&mut self.deferred, [None, None, None]);
        for (idx, cells) in deferred.into_iter().enumerate() {
            if let Some(cells) = cells {
                let mut sums = WindowSums::new(size);
                for (i, n) in cells.into_iter().enumerate() {
                    sums.push(size, c, i, n);
                }
                self.windows[idx] = Some(sums);
            }
        }
    }
}

/// Normalise one JSON file. Returns Ok(None) when the file is not a report we
/// recognise (caller decides whether that is a warning).
pub fn adapt_file(path: &Path, sample_hint: Option<&str>) -> Result<Option<ToolReport>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;
    let sample = sample_hint
        .map(str::to_string)
        .unwrap_or_else(|| sample_from_path(path));
    // A fastmaps grid is 36× the bytes of an ordinary report and every sample
    // gets one, so it is read through a typed visitor instead of a `Value`
    // tree (one enum construction per cell used to be the cost). The `Value`
    // branch further down still handles any grid this reader cannot claim.
    if text.contains("\"counts\"") {
        if let Ok(doc) = serde_json::from_str::<GridDoc>(&text) {
            if doc.is_grid() {
                return Ok(Some(doc.into_report(sample, path)));
            }
        }
    }
    let v: serde_json::Value = match serde_json::from_str(&text) {
        Ok(v) => v,
        Err(e) => return Ok(Some(ToolReport {
            sample,
            tool: "unparsable".into(),
            version: None,
            source: path.to_path_buf(),
            metrics: BTreeMap::new(),
            notes: vec![format!("not valid JSON: {e}")],
            from_v1: false,
        })),
    };
    if !v.is_object() {
        return Ok(None);
    }
    // fasthelper's own aggregate (report.json) carries the same schema string
    // as a per-tool report but is a collection, not a measurement. Skip it by
    // its generator/reports shape, so re-concatenating a run directory never
    // turns fasthelper's report.json into a phantom sample.
    if v.get("generator").and_then(|g| g.as_str()) == Some("fasthelper")
        || v.get("reports").is_some_and(|r| r.is_array())
    {
        return Ok(None);
    }

    // --- canonical v1 wrapper ---
    if v.get("schema").and_then(|s| s.as_str()) == Some(SCHEMA_VERSION) {
        let tool = v
            .get("tool")
            .and_then(|t| t.as_str())
            .unwrap_or("unknown")
            .to_string();
        let mut metrics = BTreeMap::new();
        if let Some(m) = v.get("metrics").and_then(|m| m.as_object()) {
            for (k, val) in m {
                if let Some(f) = num(val) {
                    metrics.insert(k.clone(), f);
                }
            }
        }
        // Derived convenience keys, so every sample has a comparable pair.
        let rin = metrics.get("reads_in").copied();
        let rout = metrics.get("reads_out").copied();
        if let (Some(rin), Some(rout)) = (rin, rout) {
            if rin > 0.0 {
                metrics.insert("read_retention".into(), rout / rin);
            }
        }
        if let Some(e) = v.get("elapsed_s").and_then(num) {
            metrics.insert("elapsed_s".into(), e);
        }
        return Ok(Some(ToolReport {
            sample,
            tool,
            version: v
                .get("tool_version")
                .and_then(|s| s.as_str())
                .map(str::to_string),
            source: path.to_path_buf(),
            metrics,
            notes: Vec::new(),
            from_v1: true,
        }));
    }

    // --- fastmaps heatmap.json (grid of (T+C, T+G) counts) ---
    if v.get("counts").is_some() && v.get("size").is_some() {
        let size = get_num(&v, "size") as i64;
        let windows: [(&str, &serde_json::Value); 3] = [
            ("all", &v),
            ("r1", v.get("counts_r1").unwrap_or(&serde_json::Value::Null)),
            ("r2", v.get("counts_r2").unwrap_or(&serde_json::Value::Null)),
        ];
        let mut metrics = BTreeMap::new();
        metrics.insert("total_reads".into(), get_num(&v, "total_reads"));
        metrics.insert("grid_size".into(), size as f64);
        let c = (size / 50).max(2); // corner window: 2% of the read length
        for (label, obj) in windows {
            let counts = match obj.get("counts") {
                Some(serde_json::Value::Array(a)) => a,
                _ => continue,
            };
            if counts.len() < (size * size) as usize {
                continue;
            }
            let (mut poly_a, mut poly_t, mut poly_g, mut poly_c, mut off_diag) =
                (0.0f64, 0.0, 0.0, 0.0, 0.0);
            let mut total = 0.0f64;
            for y in 0..size {
                for x in 0..size {
                    let n = counts[(y * size + x) as usize].as_f64().unwrap_or(0.0);
                    if n == 0.0 {
                        continue;
                    }
                    total += n;
                    let low_x = x <= c;
                    let high_x = x >= size - 1 - c;
                    let low_y = y <= c;
                    let high_y = y >= size - 1 - c;
                    if low_x && low_y {
                        poly_a += n;
                    }
                    if high_x && high_y {
                        poly_t += n;
                    }
                    if low_x && high_y {
                        poly_g += n;
                    }
                    if high_x && low_y {
                        poly_c += n;
                    }
                    if (x - y).abs() > size / 10 {
                        off_diag += n;
                    }
                }
            }
            if total == 0.0 {
                continue;
            }
            let p = |n: f64| n / total;
            let pre = if label == "all" { String::new() } else { format!("{label}_") };
            metrics.insert(format!("{pre}polyA_frac"), p(poly_a));
            metrics.insert(format!("{pre}polyT_frac"), p(poly_t));
            metrics.insert(format!("{pre}polyG_frac"), p(poly_g));
            metrics.insert(format!("{pre}polyC_frac"), p(poly_c));
            metrics.insert(format!("{pre}off_diagonal_frac"), p(off_diag));
            if label == "all" {
                metrics.insert("total_reads".into(), total);
            }
        }
        return Ok(Some(ToolReport {
            sample,
            tool: "fastmaps".into(),
            version: v.get("version").and_then(num).map(|n| format!("grid v{n}")),
            source: path.to_path_buf(),
            metrics,
            notes: Vec::new(),
            from_v1: false,
        }));
    }

    // --- fastprep legacy stats: {r1: {...}, r2: {...}, pairs_out, ...} ---
    // Only reached by report files written before fastprep emitted the v1
    // wrapper (PLAN.md 1.1). A v1 report is claimed by the branch above, so
    // this note never lands on a current fastprep --json file.
    if let Some(r1) = v.get("r1").filter(|x| x.is_object()) {
        let r2 = v.get("r2").filter(|x| x.is_object());
        let mut metrics = BTreeMap::new();
        let mut reads_in = 0.0;
        let mut reads_kept = 0.0;
        let mut reads_dropped = 0.0;
        let mut bases_in = 0.0;
        let mut bases_trimmed = 0.0;
        let mut adapter_trimmed = 0.0;
        let mut poly_dropped = 0.0;
        let mut hard = 0.0;
        let mut mean_q = 0.0;
        let mut adapter_dropped = 0.0;
        let mut duplicate_dropped = 0.0;
        for mate in [Some(r1), r2].into_iter().flatten() {
            reads_in += get_num(mate, "reads_in");
            reads_kept += get_num(mate, "reads_kept");
            reads_dropped += get_num(mate, "reads_dropped");
            bases_in += get_num(mate, "bases_in");
            bases_trimmed += get_num(mate, "bases_trimmed_5") + get_num(mate, "bases_trimmed_3");
            if let Some(a) = mate.get("adapter") {
                adapter_trimmed += get_num(a, "reads_trimmed");
            }
            for cause in [
                "reads_dropped_poly5",
                "reads_dropped_poly3",
                "reads_dropped_poly_int",
            ] {
                poly_dropped += get_num(mate, cause);
            }
            hard += get_num(mate, "reads_dropped_hard");
            mean_q += get_num(mate, "reads_dropped_mean_q");
            adapter_dropped += get_num(mate, "reads_dropped_adapter");
            duplicate_dropped += get_num(mate, "reads_dropped_duplicate");
        }
        // Per-cause drop counts are summed across mates: this is a read-level
        // budget, and reporting only the last mate silently halves it.
        metrics.insert("reads_dropped_hard".into(), hard);
        metrics.insert("reads_dropped_mean_q".into(), mean_q);
        metrics.insert("reads_dropped_adapter".into(), adapter_dropped);
        metrics.insert("reads_dropped_duplicate".into(), duplicate_dropped);
        metrics.insert("reads_in".into(), reads_in);
        metrics.insert("reads_out".into(), reads_kept);
        metrics.insert("reads_dropped".into(), reads_dropped);
        metrics.insert("bases_in".into(), bases_in);
        metrics.insert("bases_trimmed".into(), bases_trimmed);
        metrics.insert("reads_adapter_trimmed".into(), adapter_trimmed);
        metrics.insert("reads_dropped_polyx".into(), poly_dropped);
        if reads_in > 0.0 {
            metrics.insert("read_retention".into(), reads_kept / reads_in);
            metrics.insert("bases_trimmed_frac".into(), bases_trimmed / bases_in.max(1.0));
        }
        if v.get("pairs_out").is_some() {
            metrics.insert("pairs_out".into(), get_num(&v, "pairs_out"));
            metrics.insert("pairs_dropped".into(), get_num(&v, "pairs_dropped"));
            metrics.insert("singletons_r1".into(), get_num(&v, "singletons_r1"));
            metrics.insert("singletons_r2".into(), get_num(&v, "singletons_r2"));
        }
        return Ok(Some(ToolReport {
            sample,
            tool: "fastprep".into(),
            version: None,
            source: path.to_path_buf(),
            metrics,
            notes: vec!["legacy raw stats (no schema wrapper)".into()],
            from_v1: false,
        }));
    }

    // --- anything else: take top-level numbers if it smells like a report ---
    let obj = v.as_object().unwrap();
    let numeric: BTreeMap<String, f64> = obj
        .iter()
        .filter_map(|(k, val)| num(val).map(|n| (k.clone(), n)))
        .collect();
    if numeric.is_empty() {
        return Ok(None);
    }
    Ok(Some(ToolReport {
        sample,
        tool: obj
            .get("tool")
            .and_then(|t| t.as_str())
            .unwrap_or("unknown")
            .to_string(),
        version: obj
            .get("tool_version")
            .and_then(|t| t.as_str())
            .map(str::to_string),
        source: path.to_path_buf(),
        metrics: numeric,
        notes: vec!["unrecognised schema: top-level numbers only".into()],
        from_v1: false,
    }))
}

fn walk_json(dir: &Path, out: &mut Vec<PathBuf>, depth: usize) {
    if depth > 6 {
        return;
    }
    let Ok(rd) = std::fs::read_dir(dir) else { return };
    for entry in rd.flatten() {
        let p = entry.path();
        let name = entry.file_name().to_string_lossy().to_string();
        if p.is_dir() {
            if matches!(
                name.as_str(),
                "target" | ".git" | ".pixi" | "node_modules" | "logs" | "tmp"
            ) {
                continue;
            }
            walk_json(&p, out, depth + 1);
        } else if name.ends_with(".json")
            && !matches!(name.as_str(), "pixi.lock" | "Cargo.lock" | "manifest.json")
            && !p.to_string_lossy().contains(".rustc_info")
        {
            out.push(p);
        }
    }
}

/// Collect and normalise every report reachable from the given inputs.
pub fn collect(inputs: &[PathBuf], sample: Option<&str>) -> Result<Vec<ToolReport>> {
    let mut files = Vec::new();
    for input in inputs {
        if input.is_dir() {
            walk_json(input, &mut files, 0);
        } else if input.is_file() {
            files.push(input.clone());
        } else {
            anyhow::bail!("input not found: {}", input.display());
        }
    }
    files.sort();
    files.dedup();

    let mut reports = Vec::new();
    for f in files.iter() {
        let hint = if files.len() == 1 { sample } else { None };
        match adapt_file(f, hint) {
            Ok(Some(r)) => reports.push(r),
            Ok(None) => {}
            Err(e) => eprintln!("fasthelper: skipping {}: {e}", f.display()),
        }
    }
    // One measurement, one row: a `fastsuite.report.v1` wrapper and the raw
    // file it was derived from (fastmaps' grid beside its `<stem>.report.json`)
    // describe the same (sample, tool). Keep the v1 report — it is the tool's
    // own statement — and drop the derivation, so nothing double-counts.
    let v1: BTreeSet<(String, String)> = reports
        .iter()
        .filter(|r| r.from_v1)
        .map(|r| (r.sample.clone(), r.tool.clone()))
        .collect();
    if !v1.is_empty() {
        reports.retain(|r| r.from_v1 || !v1.contains(&(r.sample.clone(), r.tool.clone())));
    }
    Ok(reports)
}

/// Wide metric matrix: one row per (sample, tool), union of metric columns.
///
/// Borrowed, not cloned: at 10 k samples the owned copy cost ~61 MB of row
/// maps that existed only to be printed once.
pub struct Matrix<'a> {
    pub columns: Vec<&'a str>,
    pub rows: Vec<(&'a str, &'a str, &'a BTreeMap<String, f64>)>,
}

pub fn build_matrix(reports: &[ToolReport]) -> Matrix<'_> {
    let mut columns: BTreeSet<&str> = BTreeSet::new();
    let mut rows: Vec<(&str, &str, &BTreeMap<String, f64>)> = Vec::with_capacity(reports.len());
    for r in reports {
        columns.extend(r.metrics.keys().map(|k| k.as_str()));
        rows.push((r.sample.as_str(), r.tool.as_str(), &r.metrics));
    }
    Matrix {
        columns: columns.into_iter().collect(),
        rows,
    }
}

fn fmt_val(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else if v.abs() >= 1000.0 {
        format!("{v:.0}")
    } else if v.abs() >= 1.0 {
        format!("{v:.3}")
    } else {
        format!("{v:.5}")
    }
}

/// Write the wide matrix as TSV. Streaming: the 10 k table is 2 MB of cells
/// looked up in a borrowed row, so nothing but the 8 KB buffer is held.
pub fn write_matrix_tsv<W: Write>(w: &mut W, m: &Matrix<'_>) -> std::io::Result<()> {
    w.write_all(b"sample\ttool")?;
    for c in &m.columns {
        w.write_all(b"\t")?;
        w.write_all(c.as_bytes())?;
    }
    w.write_all(b"\n")?;
    for (sample, tool, metrics) in &m.rows {
        w.write_all(sample.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(tool.as_bytes())?;
        for c in &m.columns {
            w.write_all(b"\t")?;
            if let Some(v) = metrics.get(*c) {
                write!(w, "{}", fmt_val(*v))?;
            }
        }
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// In-memory convenience wrapper over [`write_matrix_tsv`] (tests, compare).
pub fn matrix_tsv(m: &Matrix<'_>) -> String {
    let mut buf = Vec::new();
    write_matrix_tsv(&mut buf, m).expect("writing a TSV into a Vec cannot fail");
    String::from_utf8(buf).expect("fmt_val and the header are ASCII")
}

pub fn matrix_md(reports: &[ToolReport]) -> String {
    let by_tool = group_by_tool(reports);
    let mut s = String::new();
    for (tool, reports) in &by_tool {
        let mut columns: BTreeSet<&str> = BTreeSet::new();
        for r in reports {
            columns.extend(r.metrics.keys().map(|k| k.as_str()));
        }
        s.push_str(&format!("## {tool}\n\n| metric |"));
        for r in reports {
            s.push_str(&format!(" {} |", r.sample));
        }
        s.push_str("\n|---|");
        for _ in reports {
            s.push_str("---|");
        }
        s.push('\n');
        for c in &columns {
            s.push_str(&format!("| `{c}` |"));
            for r in reports {
                s.push_str(&format!(
                    " {} |",
                    r.metrics.get(*c).map(|v| fmt_val(*v)).unwrap_or_default()
                ));
            }
            s.push('\n');
        }
        s.push('\n');
    }
    s
}

fn group_by_tool(reports: &[ToolReport]) -> BTreeMap<String, Vec<&ToolReport>> {
    let mut m: BTreeMap<String, Vec<&ToolReport>> = BTreeMap::new();
    for r in reports {
        m.entry(r.tool.clone()).or_default().push(r);
    }
    m
}

// ---------------------------------------------------------------- HTML ----

fn bar(value: f64, max: f64, colour: &str) -> String {
    let w = if max > 0.0 {
        ((value / max) * 180.0).clamp(0.0, 180.0)
    } else {
        0.0
    };
    format!(
        "<svg width='200' height='14' style='vertical-align:middle'>\
         <rect x='0' y='2' width='{w:.1}' height='10' fill='{colour}' rx='2'/></svg>"
    )
}

/// Self-contained HTML report: no scripts, no external assets.
pub fn render_html(title: &str, reports: &[ToolReport]) -> String {
    let by_tool = group_by_tool(reports);
    let mut out = String::new();
    out.push_str(&format!(
        "<!doctype html><html><head><meta charset='utf-8'><title>{title}</title>\
<style>\
body{{font-family:system-ui,sans-serif;margin:2rem;color:#1c1e21;background:#fff}}\
h1{{font-size:1.4rem;margin:0 0 .2rem}}h2{{font-size:1.05rem;margin:1.6rem 0 .4rem;border-bottom:1px solid #ddd;padding-bottom:.2rem}}\
table{{border-collapse:collapse;font-size:.82rem;margin:.4rem 0}}\
th,td{{border:1px solid #e3e3e3;padding:.25rem .5rem;text-align:right}}\
th:first-child,td:first-child{{text-align:left}}\
th{{background:#f6f7f9;font-weight:600}}\
tr:nth-child(even) td{{background:#fbfbfc}}\
.meta{{color:#616770;font-size:.85rem;margin-bottom:1rem}}\
.note{{color:#8a6d3b;font-size:.8rem}}\
code{{background:#f2f3f5;padding:0 .25rem;border-radius:3px}}\
</style></head><body>"
    ));
    out.push_str(&format!("<h1>{title}</h1>"));
    out.push_str(&format!(
        "<div class='meta'>samples: {} &middot; reports: {} &middot; tools: {}</div>",
        reports
            .iter()
            .map(|r| r.sample.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        reports.len(),
        by_tool.keys().cloned().collect::<Vec<_>>().join(", ")
    ));

    for (tool, reps) in &by_tool {
        let versions: BTreeSet<String> = reps
            .iter()
            .filter_map(|r| r.version.clone())
            .collect();
        out.push_str(&format!(
            "<h2>{tool}{}</h2>",
            if versions.is_empty() {
                String::new()
            } else {
                format!(" <span class='meta'>{}</span>", versions.into_iter().collect::<Vec<_>>().join(", "))
            }
        ));
        let mut columns: BTreeSet<&str> = BTreeSet::new();
        for r in reps {
            columns.extend(r.metrics.keys().map(|k| k.as_str()));
        }
        out.push_str("<table><tr><th>metric</th>");
        for r in reps {
            out.push_str(&format!("<th>{}</th>", r.sample));
        }
        out.push_str("</tr>");
        for c in &columns {
            let numeric: Vec<Option<f64>> = reps.iter().map(|r| r.metrics.get(*c).copied()).collect();
            let max = numeric.iter().flatten().cloned().fold(0.0f64, f64::max);
            let is_rate = c.ends_with("_frac") || c.ends_with("rate") || c.ends_with("retention");
            out.push_str(&format!("<tr><td><code>{c}</code></td>"));
            for v in &numeric {
                match v {
                    Some(v) => {
                        if is_rate && max > 0.0 {
                            out.push_str(&format!(
                                "<td>{} {}</td>",
                                fmt_val(*v),
                                bar(*v, max, "#4c78a8")
                            ));
                        } else {
                            out.push_str(&format!("<td>{}</td>", fmt_val(*v)));
                        }
                    }
                    None => out.push_str("<td></td>"),
                }
            }
            out.push_str("</tr>");
        }
        out.push_str("</table>");
        for r in reps {
            for n in &r.notes {
                out.push_str(&format!(
                    "<div class='note'>{} — {}: {}</div>",
                    r.sample, r.tool, n
                ));
            }
        }
    }
    out.push_str("</body></html>");
    out
}

/// One report as it appears inside the aggregate document. Field order is
/// declaration order, and it is alphabetical — the aggregate is a `BTreeMap`
/// document today and its bytes are part of the contract.
#[derive(Serialize)]
struct JsonRow<'a> {
    metrics: &'a BTreeMap<String, f64>,
    notes: &'a [String],
    sample: &'a str,
    source: Cow<'a, str>,
    tool: &'a str,
    tool_version: Option<&'a str>,
}

#[derive(Serialize)]
struct JsonDoc<'a> {
    generator: &'a str,
    generator_version: &'a str,
    reports: Vec<JsonRow<'a>>,
    schema: &'a str,
    title: &'a str,
}

/// Serialise the aggregate straight into `w`, borrowing the reports.
///
/// The document used to be a `serde_json::Value` tree plus one pretty `String`
/// before the first byte reached the file: at 10 k samples (20 k reports) that
/// peaked at 273 MB of RSS for a 23 MB file, ~20 k objects each with ~30 key
/// `String`s. Only the 20 k-row view below is built now, and the buffer belongs
/// to the caller's `BufWriter`.
pub fn write_json<W: Write>(w: &mut W, title: &str, reports: &[ToolReport]) -> Result<()> {
    let rows: Vec<JsonRow<'_>> = reports
        .iter()
        .map(|r| JsonRow {
            metrics: &r.metrics,
            notes: &r.notes,
            sample: &r.sample,
            source: r.source.to_string_lossy(),
            tool: &r.tool,
            tool_version: r.version.as_deref(),
        })
        .collect();
    let doc = JsonDoc {
        generator: "fasthelper",
        generator_version: env!("CARGO_PKG_VERSION"),
        reports: rows,
        schema: SCHEMA_VERSION,
        title,
    };
    serde_json::to_writer_pretty(&mut *w, &doc).context("serialising report.json")?;
    w.write_all(b"\n")?;
    Ok(())
}

/// Write the aggregate report in the requested formats. Returns written paths.
///
/// Every format is streamed into its own `BufWriter`; `matrix` is the wide
/// table built once by the caller (so `--stdout` does not build it twice).
pub fn write_reports(
    title: &str,
    reports: &[ToolReport],
    outdir: &Path,
    formats: &[String],
    matrix: Option<&Matrix<'_>>,
) -> Result<Vec<PathBuf>> {
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("creating {}", outdir.display()))?;
    let mut written = Vec::new();
    let mut owned_matrix: Option<Matrix<'_>> = None;
    for f in formats {
        let name = match f.trim() {
            "json" => "report.json",
            "tsv" => "report.tsv",
            "md" => "report.md",
            "html" => "report.html",
            other => anyhow::bail!("unknown format '{other}' (json,tsv,md,html)"),
        };
        let path = outdir.join(name);
        let file = std::fs::File::create(&path)
            .with_context(|| format!("creating {}", path.display()))?;
        let mut w = BufWriter::new(file);
        match f.trim() {
            "json" => write_json(&mut w, title, reports)?,
            "tsv" => {
                if matrix.is_none() && owned_matrix.is_none() {
                    owned_matrix = Some(build_matrix(reports));
                }
                write_matrix_tsv(&mut w, matrix.or(owned_matrix.as_ref()).expect("built above"))?;
            }
            "md" => w.write_all(matrix_md(reports).as_bytes())?,
            _ => w.write_all(render_html(title, reports).as_bytes())?,
        }
        w.flush()
            .with_context(|| format!("writing {}", path.display()))?;
        written.push(path);
    }
    Ok(written)
}

/// One `metric OP value` threshold, e.g. `polyG_frac>0.05`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Op {
    Gt,
    Lt,
    Ge,
    Le,
    Eq,
    Ne,
}

impl Op {
    fn as_str(&self) -> &'static str {
        match self {
            Op::Gt => ">",
            Op::Lt => "<",
            Op::Ge => ">=",
            Op::Le => "<=",
            Op::Eq => "==",
            Op::Ne => "!=",
        }
    }
}

/// One `metric OP value` threshold, e.g. `polyG_frac>0.05`.
#[derive(Debug, Clone)]
pub struct Threshold {
    pub metric: String,
    pub op: Op,
    pub value: f64,
}

impl Threshold {
    /// Parse `metric>0.05` / `reads<=100` / ... Ops: > < >= <= == !=.
    pub fn parse(expr: &str) -> Result<Self> {
        // Two-char operators first, so `>=` is never read as `>`.
        const OPS: &[(&str, Op)] = &[
            (">=", Op::Ge),
            ("<=", Op::Le),
            ("==", Op::Eq),
            ("!=", Op::Ne),
            (">", Op::Gt),
            ("<", Op::Lt),
        ];
        for (op, kind) in OPS {
            if let Some(i) = expr.find(op) {
                let metric = expr[..i].trim().to_string();
                let value: f64 = expr[i + op.len()..]
                    .trim()
                    .parse()
                    .with_context(|| format!("parsing threshold '{expr}': value is not a number"))?;
                if metric.is_empty() {
                    anyhow::bail!("threshold '{expr}' has no metric name");
                }
                return Ok(Threshold {
                    metric,
                    op: kind.clone(),
                    value,
                });
            }
        }
        anyhow::bail!("threshold '{expr}' is not METRIC OP VALUE (ops: > < >= <= == !=)");
    }

    pub fn tripped(&self, actual: f64) -> bool {
        match self.op {
            Op::Gt => actual > self.value,
            Op::Lt => actual < self.value,
            Op::Ge => actual >= self.value,
            Op::Le => actual <= self.value,
            Op::Eq => (actual - self.value).abs() < 1e-12,
            Op::Ne => (actual - self.value).abs() >= 1e-12,
        }
    }
}

/// Evaluate thresholds against EVERY report, not a flattened map.
///
/// Each tripped metric is returned labelled with the sample it came from, so a
/// violation in one sample of a multi-sample run is never masked by another
/// sample's value ("last wins" hid exactly that). Returns the tripped
/// descriptions, e.g. `a: duplicate_rate=0.5 > 0.01`.
pub fn evaluate_thresholds(reports: &[ToolReport], thresholds: &[Threshold]) -> Vec<String> {
    let mut out = Vec::new();
    for r in reports {
        for t in thresholds {
            if let Some(actual) = r.metrics.get(&t.metric) {
                if t.tripped(*actual) {
                    out.push(format!(
                        "{}: {}={} {} {}",
                        r.sample,
                        t.metric,
                        actual,
                        t.op.as_str(),
                        t.value
                    ));
                }
            }
        }
    }
    out
}

/// True when `path` is the aggregate `run` writes: `generator` names
/// fasthelper and the payload is a `reports` array of per-tool entries.
///
/// `collect` only recognises per-tool reports, so feeding it an aggregate used
/// to produce "no recognisable reports found" — a message that denies the file
/// is a report at all, when it is one, of a different kind.
fn is_aggregate(path: &Path) -> bool {
    let Ok(text) = std::fs::read_to_string(path) else {
        return false;
    };
    let Ok(v) = serde_json::from_str::<serde_json::Value>(&text) else {
        return false;
    };
    v["generator"] == "fasthelper" && v["reports"].is_array()
}

pub fn cmd_concatenate(a: &ReportArgs) -> Result<()> {
    let reports = collect(&a.inputs, a.sample.as_deref())?;
    if reports.is_empty() {
        let aggregates: Vec<String> = a
            .inputs
            .iter()
            .filter(|p| p.is_file() && is_aggregate(p))
            .map(|p| p.display().to_string())
            .collect();
        if !aggregates.is_empty() {
            anyhow::bail!(
                "{} is a fasthelper aggregate (\"generator\": \"fasthelper\", a \"reports\" array), \
                 not a per-tool report: concatenate merges per-tool reports such as \
                 <sample>/fastprep.json. To combine aggregates, pass the report paths each one \
                 lists in its manifest (or run the tools into a single --outdir and let `run` \
                 aggregate them).",
                aggregates.join(", ")
            );
        }
        anyhow::bail!(
            "no recognisable reports found in {}",
            a.inputs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    let formats: Vec<String> = a.format.split(',').map(|s| s.trim().to_string()).collect();
    let title = "FastSuite report";
    // Built once: the TSV writer and --stdout are the same table, and the
    // matrix borrows the reports instead of cloning them.
    let matrix = if a.stdout || formats.iter().any(|f| f == "tsv") {
        Some(build_matrix(&reports))
    } else {
        None
    };
    let written = write_reports(title, &reports, &a.outdir, &formats, matrix.as_ref())?;
    if !a.fail_on.is_empty() {
        let thresholds: Vec<Threshold> = a
            .fail_on
            .iter()
            .map(|s| Threshold::parse(s))
            .collect::<Result<Vec<_>>>()?;
        let tripped = evaluate_thresholds(&reports, &thresholds);
        if !tripped.is_empty() {
            for t in &tripped {
                eprintln!("fasthelper: threshold tripped: {t}");
            }
            eprintln!(
                "fasthelper: {} threshold(s) tripped — report written, exiting 2",
                tripped.len()
            );
            std::process::exit(2);
        }
    }
    if a.stdout {
        let stdout = std::io::stdout();
        let mut w = BufWriter::new(stdout.lock());
        write_matrix_tsv(&mut w, matrix.as_ref().expect("built when --stdout"))
            .context("printing report.tsv")?;
        w.flush()?;
    }
    println!(
        "fasthelper: concatenated {} report(s) from {} sample(s) -> {}",
        reports.len(),
        reports
            .iter()
            .map(|r| r.sample.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        written
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
    Ok(())
}
