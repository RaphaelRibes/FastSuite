//! Regression tests for the report adapters — the part of fasthelper that
//! silently wrong-numbers a report if it breaks.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use fasthelper::report::{adapt_file, build_matrix, matrix_tsv, sample_from_path, ToolReport};

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
