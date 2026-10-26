//! Regression tests for the report adapters — the part of fasthelper that
//! silently wrong-numbers a report if it breaks.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use fasthelper::report::{
    adapt_file, build_matrix, collect, matrix_tsv, sample_from_path, write_json, ToolReport,
};

fn tmp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fasthelper_test_{}_{}", std::process::id(), name));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn sample_name_walks_up_past_generic_dirs() {
    assert_eq!(sample_from_path(Path::new("R1.fastq")), "R1");
    // <run>/E94/maps/heatmap.json -> E94, not "maps" and not "heatmap"
    assert_eq!(sample_from_path(Path::new("runs/demo/E94/maps/heatmap.json")), "E94");
    assert_eq!(sample_from_path(Path::new("runs/demo/E94/fastprep.json")), "E94");
    // A non-generic stem wins over the parent directory.
    assert_eq!(sample_from_path(Path::new("runs/demo/stats.json")), "demo");
    // Every tool's own report stem is generic: a QC-only chain writes
    // <sample>/fastcheck.json, so the directory carries the identity.
    assert_eq!(sample_from_path(Path::new("runs/demo/E94/fastcheck.json")), "E94");
}

/// Task 2.21 defect 1: the `.report` wrapper written beside a grid must resolve
/// to the SAMPLE, not to a phantom one.
///
/// `heatmap.report.json` has the stem `heatmap.report`. `heatmap` is generic, so
/// the wrapper has to be stripped to `heatmap` and then fall through to the
/// parent directory. Before the fix the stem `heatmap.report` matched nothing,
/// so three samples produced four rows and the summary read
/// "4 report(s) from 3 sample(s)".
#[test]
fn the_report_wrapper_beside_a_grid_resolves_to_its_sample() {
    // The exact shape the defect was reported against.
    assert_eq!(
        sample_from_path(Path::new("runs/demo/E94/heatmap.report.json")),
        "E94"
    );
    // A bare wrapper with no sample directory above it cannot invent one: the
    // suffix is still stripped, so it lands on the generic "heatmap" rather than
    // the literal stem "heatmap.report". That is the degenerate case (a report
    // at the filesystem root with nothing to walk up to); the interesting one is
    // the line above, where the directory supplies the identity.
    assert_eq!(sample_from_path(Path::new("heatmap.report.json")), "heatmap");
    // The strip is GATED on the remainder being generic, so a library genuinely
    // named `align.report` keeps its own name instead of losing the suffix.
    assert_eq!(
        sample_from_path(Path::new("runs/demo/align.report.json")),
        "align.report"
    );
    // And the non-wrapper spelling still works, which is the case that was
    // already covered: both spellings must agree on the sample.
    assert_eq!(
        sample_from_path(Path::new("runs/demo/E94/heatmap.report.json")),
        sample_from_path(Path::new("runs/demo/E94/heatmap.json"))
    );
}

#[test]
fn html_and_md_escape_the_sample_name() {
    // A sample directory may legally be named `a<b&c|d`; unescaped it closed the
    // HTML tag and broke the table. TSV stays raw — it is not markup.
    let reports = vec![ToolReport {
        sample: "a<b&c|d".into(),
        tool: "fastcheck".into(),
        version: Some("0.1.0".into()),
        source: PathBuf::from("/tmp/a<b&c|d/fastcheck.json"),
        metrics: BTreeMap::from([("reads".to_string(), 2.0)]),
        notes: vec![],
        from_v1: true,
    }];

    let html = fasthelper::report::render_html("T", &reports);
    assert!(
        html.contains("<th>a&lt;b&amp;c|d</th>"),
        "the sample column must be escaped: {html}"
    );
    assert!(!html.contains("<b&c|d"), "raw markup leaked into the HTML: {html}");

    let md = fasthelper::report::matrix_md(&reports);
    assert!(
        md.contains(" a&lt;b&amp;c|d |"),
        "the markdown header must be escaped: {md}"
    );

    let mut tsv = Vec::new();
    fasthelper::report::write_matrix_tsv(&mut tsv, &fasthelper::report::build_matrix(&reports))
        .unwrap();
    assert!(
        String::from_utf8(tsv).unwrap().contains("a<b&c|d\t"),
        "TSV is not markup and must stay raw"
    );
}

#[test]
fn heatmap_corner_fractions_match_ground_truth() {
    // Grid 151x151 (read length 150), index = y*size + x with x = T+C, y = T+G.
    // 3 polyT reads -> (150,150); 2 polyG reads -> (0,150); 1 read at (75,75).
    let size = 151usize;
    let mut counts = vec![0u64; size * size];
    counts[150 * size + 150] = 3;
    counts[150 * size + 0] = 2;
    counts[75 * size + 75] = 1;
    let doc = serde_json::json!({
        "version": 2,
        "size": size,
        "total_reads": 6,
        "counts": counts,
    });
    let path = tmp("heatmap").join("heatmap.json");
    std::fs::write(&path, serde_json::to_string(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(r.tool, "fastmaps");
    let m = &r.metrics;
    assert!((m["polyT_frac"] - 3.0 / 6.0).abs() < 1e-9, "{m:?}");
    assert!((m["polyG_frac"] - 2.0 / 6.0).abs() < 1e-9, "{m:?}");
    assert!((m["off_diagonal_frac"] - 2.0 / 6.0).abs() < 1e-9, "{m:?}");
    assert_eq!(m["polyA_frac"], 0.0);
}

#[test]
fn legacy_fastprep_stats_are_summed_across_mates() {
    let doc = serde_json::json!({
        "r1": {"reads_in": 100, "reads_kept": 90, "reads_dropped": 10, "bases_in": 15000,
               "bases_trimmed_5": 100, "bases_trimmed_3": 50,
               "reads_dropped_adapter": 4, "adapter": {"reads_trimmed": 7}},
        "r2": {"reads_in": 100, "reads_kept": 80, "reads_dropped": 20, "bases_in": 15000,
               "bases_trimmed_5": 0, "bases_trimmed_3": 50,
               "reads_dropped_adapter": 6, "adapter": {"reads_trimmed": 3}}
    });
    let path = tmp("prep").join("fastprep.json");
    std::fs::write(&path, serde_json::to_string(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, Some("S1")).unwrap().unwrap();
    let m = &r.metrics;
    assert_eq!(r.sample, "S1");
    assert_eq!(m["reads_in"], 200.0);
    assert_eq!(m["reads_out"], 170.0);
    // summed across mates, not last-mate-wins
    assert_eq!(m["reads_dropped_adapter"], 10.0);
    assert_eq!(m["reads_adapter_trimmed"], 10.0);
    assert_eq!(m["bases_trimmed"], 200.0);
    assert!((m["read_retention"] - 0.85).abs() < 1e-9);
}

#[test]
fn fastprep_v1_report_goes_through_the_v1_path_without_the_legacy_note() {
    // The shape fastprep now writes (PLAN.md 1.1): the wrapper, with the raw
    // stats object parked under `legacy`. The legacy adapter must NOT claim
    // it — that note belongs to files written before the wrapper existed.
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastprep",
        "tool_version": "0.1.0",
        "mode": "paired",
        "dry_run": false,
        "command": "fastprep -1 r1.fq -2 r2.fq -o o1.fq -p o2.fq -j 1",
        "inputs": {"r1": "r1.fq", "r2": "r2.fq"},
        "outputs": {"r1": "o1.fq", "r2": "o2.fq"},
        "elapsed_s": 1.0,
        "metrics": {
            "reads_in": 200, "reads_out": 170, "reads_dropped": 30,
            "r1_reads_in": 100, "r2_reads_in": 100,
            "pairs_out": 85
        },
        "legacy": {
            "r1": {"reads_in": 100, "reads_kept": 90, "reads_dropped": 10},
            "r2": {"reads_in": 100, "reads_kept": 80, "reads_dropped": 20},
            "pairs_out": 85
        }
    });
    let path = tmp("prep_v1").join("fastprep.json");
    std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(r.tool, "fastprep");
    assert_eq!(r.version.as_deref(), Some("0.1.0"));
    // The whole point of the wrapper: no legacy note, metrics straight from
    // the report.
    assert!(
        !r.notes.iter().any(|n| n.contains("legacy raw stats")),
        "fastprep's v1 report must not be noted as legacy: {:?}",
        r.notes
    );
    assert_eq!(r.metrics["reads_in"], 200.0);
    assert_eq!(r.metrics["reads_out"], 170.0);
    assert_eq!(r.metrics["r2_reads_in"], 100.0);
    assert_eq!(r.metrics["pairs_out"], 85.0);
    assert!((r.metrics["read_retention"] - 0.85).abs() < 1e-9);
}

#[test]
fn legacy_fastprep_file_still_gets_the_legacy_note() {
    // Files written before the wrapper existed keep the note: the branch is
    // still reachable, and still warns that the numbers were derived.
    let doc = serde_json::json!({
        "r1": {"reads_in": 100, "reads_kept": 90, "reads_dropped": 10},
        "r2": {"reads_in": 100, "reads_kept": 80, "reads_dropped": 20}
    });
    let path = tmp("prep_legacy").join("fastprep.json");
    std::fs::write(&path, serde_json::to_string(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(r.tool, "fastprep");
    assert!(
        r.notes.iter().any(|n| n.contains("legacy raw stats")),
        "the pre-wrapper file must still be noted: {:?}",
        r.notes
    );
    assert_eq!(r.metrics["reads_in"], 200.0);
}

#[test]
fn v1_wrapper_is_read_without_the_legacy_path() {
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastdedup",
        "tool_version": "1.3.0",
        "elapsed_s": 1.5,
        "metrics": {"reads_in": 200, "reads_out": 190, "duplicate_rate": 0.05}
    });
    let path = tmp("v1").join("fastdedup.json");
    std::fs::write(&path, serde_json::to_string(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(r.tool, "fastdedup");
    assert_eq!(r.version.as_deref(), Some("1.3.0"));
    assert!(r.notes.is_empty(), "v1 must not go through the legacy adapter");
    assert!((r.metrics["read_retention"] - 0.95).abs() < 1e-9);
    assert_eq!(r.metrics["elapsed_s"], 1.5);
}

#[test]
fn matrix_has_one_row_per_sample_tool_pair_and_union_columns() {
    let mk = |sample: &str, tool: &str, k: &str, v: f64| ToolReport {
        sample: sample.into(),
        tool: tool.into(),
        version: None,
        source: PathBuf::from("x.json"),
        metrics: BTreeMap::from([(k.to_string(), v)]),
        notes: vec![],
        from_v1: true,
    };
    let reports = vec![
        mk("A", "fastdedup", "reads_in", 10.0),
        mk("A", "fastprep", "reads_out", 9.0),
        mk("B", "fastdedup", "reads_in", 20.0),
    ];
    let tsv = matrix_tsv(&build_matrix(&reports));
    let lines: Vec<_> = tsv.lines().collect();
    assert_eq!(lines.len(), 4);
    assert_eq!(lines[0], "sample\ttool\treads_in\treads_out");
    assert!(tsv.contains("A\tfastprep\t\t9"));
}

// ---------------------------------------------------------------- FH-1 ----

/// A grid exactly as `fastmaps -o <dir> …` writes it (counts + metadata).
fn write_grid(path: &Path, size: u64) {
    let mut counts = vec![0u64; (size * size) as usize];
    counts[0] = 4;
    counts[(size * size - 1) as usize] = 2;
    let doc = serde_json::json!({
        "version": 2,
        "size": size,
        "files": ["r1.fastq.gz", "r2.fastq.gz"],
        "file_kinds": {"r1": 1, "r2": 1},
        "total_reads": 6,
        "per_file_counts": [4, 2],
        "preview_reads_per_pixel": 1,
        "counts": counts,
        "counts_r1": counts,
        "counts_r2": counts,
        "data_row_offsets": [0],
        "preview_row_offsets": [0],
    });
    std::fs::write(path, serde_json::to_string(&doc).unwrap()).unwrap();
}

/// The v1 wrapper `fastmaps --json <dir>/heatmap.json` leaves beside the grid
/// (`format.rs`: `<stem>.report.json`).
fn write_maps_wrapper(path: &Path) {
    let doc = serde_json::json!({
        "schema": "fastsuite.report.v1",
        "tool": "fastmaps",
        "tool_version": "0.1.0",
        "mode": "composition",
        "elapsed_s": 0.06,
        "inputs": {"files": ["r1.fastq.gz", "r2.fastq.gz"]},
        "metrics": {
            "no_call_reads": 0,
            "off_diagonal_frac": 0.111275,
            "polyA_frac": 0.001175,
            "polyC_frac": 0.0,
            "polyG_frac": 0.0,
            "polyT_frac": 0.0,
            "total_reads": 6
        }
    });
    std::fs::write(path, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
}

#[test]
fn fastmaps_wrapper_beside_its_grid_is_not_a_phantom_sample() {
    // FH-1 regression: `fastmaps --json <sample>/maps/heatmap.json` writes the
    // grid AND `<sample>/maps/heatmap.report.json`. Only an exact stem was
    // tested against GENERIC, so `heatmap.report` became a sample name: both
    // libraries collapsed into one phantom sample and each appeared twice.
    let root = tmp("fh1_phantom");
    let _ = std::fs::remove_dir_all(&root);
    for sample in ["alpha", "beta"] {
        let maps = root.join(sample).join("maps");
        std::fs::create_dir_all(&maps).unwrap();
        write_grid(&maps.join("heatmap.json"), 4);
        write_maps_wrapper(&maps.join("heatmap.report.json"));
    }

    // The name resolution itself: the `.report` suffix is stripped before the
    // GENERIC test, so the directory carries the identity — but only when the
    // stripped stem is generic (a library named `align.report` keeps its name).
    assert_eq!(
        sample_from_path(Path::new("run/alpha/maps/heatmap.report.json")),
        "alpha"
    );
    assert_eq!(
        sample_from_path(Path::new("libs/align.report.json")),
        "align.report"
    );

    let reports = collect(&[root.clone()], None).unwrap();
    let rows: Vec<(String, String)> = reports
        .iter()
        .map(|r| (r.sample.clone(), r.tool.clone()))
        .collect();
    assert_eq!(
        rows,
        vec![
            ("alpha".to_string(), "fastmaps".to_string()),
            ("beta".to_string(), "fastmaps".to_string())
        ],
        "one row per (sample, tool) — the grid must not survive beside its wrapper"
    );
    // The survivor is the v1 wrapper: it carries the tool version and the
    // metrics only the wrapper knows (`grid_size` is the grid adapter's).
    for r in &reports {
        assert!(r.from_v1, "the v1 report must win the dedup: {r:?}");
        assert_eq!(r.version.as_deref(), Some("0.1.0"));
        assert_eq!(r.metrics["total_reads"], 6.0);
        assert!(r.metrics.contains_key("no_call_reads"), "{:?}", r.metrics);
        assert!(!r.metrics.contains_key("grid_size"), "{:?}", r.metrics);
    }
}

// ---------------------------------------------------------------- FH-5 ----

/// The pre-streaming reader, verbatim, as an oracle: x-major accumulation,
/// `n == 0.0` shortcut, pitch as f64.
///
/// One change (2026-10-26): the row now also states the corner window the
/// fractions were built with (`max(size/50, 2)`), because copying the grid's
/// declared value published a window the numbers had not used. The fractions
/// themselves are still the old reader's, bit for bit.
fn reference_grid_metrics(
    size: i64,
    counts: &[u64],
    header_total: f64,
) -> BTreeMap<String, f64> {
    let mut metrics = BTreeMap::new();
    metrics.insert("total_reads".to_string(), header_total);
    metrics.insert("grid_size".to_string(), size as f64);
    let c = (size / 50).max(2);
    if size > 0 {
        metrics.insert("corner_window".to_string(), c as f64);
    }
    if counts.len() < (size * size) as usize {
        return metrics;
    }
    let (mut poly_a, mut poly_t, mut poly_g, mut poly_c, mut off_diag) = (0.0f64, 0.0, 0.0, 0.0, 0.0);
    let mut total = 0.0f64;
    for y in 0..size {
        for x in 0..size {
            let n = counts[(y * size + x) as usize] as f64;
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
        return metrics;
    }
    let p = |n: f64| n / total;
    metrics.insert("polyA_frac".to_string(), p(poly_a));
    metrics.insert("polyT_frac".to_string(), p(poly_t));
    metrics.insert("polyG_frac".to_string(), p(poly_g));
    metrics.insert("polyC_frac".to_string(), p(poly_c));
    metrics.insert("off_diagonal_frac".to_string(), p(off_diag));
    metrics.insert("total_reads".to_string(), total);
    metrics
}

#[test]
fn streamed_grid_matches_the_value_based_reference_bit_for_bit() {
    // Every cell of a raster is visited in the same order and the same f64
    // operations are applied, so the fractions must be equal — not "close".
    let size = 151i64;
    let mut counts = vec![0u64; (size * size) as usize];
    let mut seed = 12345u64;
    for (i, c) in counts.iter_mut().enumerate() {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let r = (seed >> 33) as u64 % 97;
        // Sparse, like a real heatmap: most cells are empty.
        *c = if r < 3 { r * 7 + (i % 5) as u64 } else { 0 };
    }
    counts[0] = 11;
    counts[(size * size - 1) as usize] = 13;
    let header_total: u64 = counts.iter().sum();

    let doc = serde_json::json!({
        "version": 2,
        "size": size,
        "total_reads": header_total,
        "counts": counts,
    });
    let path = tmp("fh5_grid").join("heatmap.json");
    std::fs::write(&path, serde_json::to_string(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(r.metrics, reference_grid_metrics(size, &counts, header_total as f64));
    assert_eq!(r.version.as_deref(), Some("grid v2"));
}

#[test]
fn a_grid_publishes_the_corner_window_it_actually_used() {
    // H8: the grid's declared `corner_window` was copied into the metrics while
    // the fractions were computed with `size/50`, so a row named a window it
    // had not used. A cell at distance 3 from the corner is inside a window of
    // 4 and outside a window of 2, so the two readings are distinguishable.
    let size = 9i64;
    let mut counts = vec![0u64; (size * size) as usize];
    counts[(8 * size + 8) as usize] = 5; // the far corner
    counts[(5 * size + 5) as usize] = 7; // inside a window of 4, outside of 2
    let doc = serde_json::json!({
        "version": 2,
        "size": size,
        "total_reads": 12,
        "corner_window": 4, // declared by the producer, NOT what is used
        "counts": counts,
    });
    let path = tmp("fh8_window").join("heatmap.json");
    std::fs::write(&path, serde_json::to_string(&doc).unwrap()).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(
        r.metrics["corner_window"], 2.0,
        "the published window must be the one the fractions were built with: {:?}",
        r.metrics
    );
    let mut expected = reference_grid_metrics(size, &counts, 12.0);
    expected.insert("corner_window".to_string(), 2.0);
    assert_eq!(r.metrics, expected, "fractions and the published window must agree");
}

#[test]
fn streamed_grid_keeps_short_windows_and_key_order_independent() {
    // `counts_r1`/`counts_r2` are bare arrays in fastmaps' output, which the
    // reader has always ignored (it looks for a nested "counts" key); a window
    // shorter than size*size contributes nothing at all.
    let size = 5i64;
    let full: Vec<u64> = (0..(size * size)).map(|i| (i % 3) as u64).collect();
    let short: Vec<u64> = full[..3].to_vec();
    // `counts` BEFORE `size`: the reader must not depend on field order.
    let json = format!(
        "{{\"counts\": {}, \"counts_r1\": {}, \"total_reads\": 9, \"size\": {size}, \"version\": 3}}",
        serde_json::to_string(&full).unwrap(),
        serde_json::to_string(&short).unwrap(),
    );
    let path = tmp("fh5_order").join("heatmap.json");
    std::fs::write(&path, json).unwrap();

    let r = adapt_file(&path, None).unwrap().unwrap();
    assert_eq!(r.metrics, reference_grid_metrics(size, &full, 9.0));
    assert!(
        !r.metrics.contains_key("r1_polyA_frac"),
        "a bare counts_r1 array is not a window: {:?}",
        r.metrics
    );
}

// ---------------------------------------------------------------- FH-3 ----

#[test]
fn streamed_report_json_is_byte_identical_to_the_value_tree() {
    // The streaming writer must reproduce the `serde_json::Value` +
    // `to_string_pretty` document exactly, including the alphabetical key
    // order that a `json!` BTreeMap gives for free.
    let reports = vec![
        ToolReport {
            sample: "A".into(),
            tool: "fastcheck".into(),
            version: Some("0.1.0".into()),
            source: PathBuf::from("/tmp/a fastcheck.json"),
            metrics: BTreeMap::from([
                ("reads_in".to_string(), 10.0),
                ("gc_fraction".to_string(), 0.4125),
                ("huge".to_string(), 1.0e300),
                ("negzero".to_string(), -0.0),
                ("integral".to_string(), 40000.0),
                ("nan".to_string(), f64::NAN),
            ]),
            notes: vec!["a note".into()],
            from_v1: true,
        },
        ToolReport {
            sample: "B".into(),
            tool: "fastprep".into(),
            version: None,
            source: PathBuf::from("b.json"),
            metrics: BTreeMap::new(),
            notes: vec![],
            from_v1: false,
        },
    ];
    let expected = {
        let samples: Vec<serde_json::Value> = reports
            .iter()
            .map(|r| {
                serde_json::json!({
                    "sample": r.sample,
                    "tool": r.tool,
                    "tool_version": r.version,
                    "source": r.source.to_string_lossy(),
                    "metrics": r.metrics,
                    "notes": r.notes,
                })
            })
            .collect();
        // The aggregate is a COLLECTION of reports, not a measurement: it is
        // stamped with its own schema tag, because `fastsuite.report.v1`
        // promises `tool` + `metrics` that an aggregate does not have (and
        // `check --report` used to reject fasthelper's own output over it).
        serde_json::to_string_pretty(&serde_json::json!({
            "schema": fasthelper::report::AGGREGATE_SCHEMA_VERSION,
            "generator": "fasthelper",
            "generator_version": env!("CARGO_PKG_VERSION"),
            "title": "T",
            "reports": samples,
        }))
        .unwrap() + "\n"
    };
    let mut out = Vec::new();
    write_json(&mut out, "T", &reports).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

