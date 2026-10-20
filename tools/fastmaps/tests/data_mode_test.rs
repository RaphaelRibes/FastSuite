//! FM-1 — `--data-mode exact|capped:N|off`.
//!
//! The per-read index stream exists for exactly one output: `heatmap.data`.
//! Making it optional must not move a single number anywhere else, and a
//! non-exact mode must say so in both the grid JSON and the v1 wrapper.
//!
//! Fixture: 14 reads of length 20, so pixel membership is deep enough for a
//! cap of 2 to bite (4 reads land on (0,0), 4 on (20,20), 4 on (0,20),
//! 2 on (20,0); `size` = 21).

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

const SIZE: usize = 21;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn fixture(dir: &Path) -> String {
    let fq = dir.join("mix_R1.fastq");
    let mut body = String::new();
    for (seq, n) in [("A", 4), ("T", 4), ("G", 4), ("C", 2)] {
        for i in 0..n {
            let s = seq.repeat(20);
            body.push_str(&format!("@p{seq}{i}\n{s}\n+\n{}\n", "I".repeat(s.len())));
        }
    }
    std::fs::write(&fq, body).unwrap();
    fq.to_string_lossy().to_string()
}

struct Run {
    out: PathBuf,
    report: Option<PathBuf>,
    stdout: String,
    stderr: String,
}

/// One binary run. `extra` is appended verbatim; `--json` is added by the
/// caller when it wants the v1 wrapper.
fn run(tag: &str, fq: &str, extra: &[&str]) -> Run {
    let dir = std::env::temp_dir().join(format!("maps_dm_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("out");
    let mut c = Command::new(bin());
    c.arg("-o").arg(&out).arg("--force").arg(fq).args(extra);
    let o = c.output().unwrap();
    assert!(
        o.status.success(),
        "--data-mode run {tag} failed\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    // `--json PATH` is the v1 report itself (the grid stays in -o).
    let report = extra
        .windows(2)
        .find(|w| w[0] == "--json")
        .map(|w| PathBuf::from(w[1]));
    Run {
        out,
        report,
        stdout: String::from_utf8_lossy(&o.stdout).to_string(),
        stderr: String::from_utf8_lossy(&o.stderr).to_string(),
    }
}

fn json(r: &Run) -> Value {
    let text = std::fs::read_to_string(r.out.join("heatmap.json")).unwrap();
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("heatmap.json must parse: {e}\n{text}"))
}

fn data(r: &Run) -> Vec<u8> {
    std::fs::read(r.out.join("heatmap.data")).unwrap()
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_dm_fx_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Everything except the per-read stream must be identical between modes:
/// the grid is built from `state.counts`, which the index stream never feeds.
#[test]
fn off_and_capped_leave_counts_and_png_exact() {
    let fx = tmpdir("counts");
    let fq = fixture(&fx);

    let exact = run("exact", &fq, &["--data-mode", "exact"]);
    let capped = run("capped", &fq, &["--data-mode", "capped:2"]);
    let off = run("off", &fq, &["--data-mode", "off"]);

    let je = json(&exact);
    for (name, r) in [("capped:2", &capped), ("off", &off)] {
        let j = json(r);
        for key in ["size", "total_reads", "counts", "counts_r1", "counts_r2", "per_file_counts"] {
            assert_eq!(j[key], je[key], "{name}: {key} must stay exact");
        }
        assert_eq!(j["preview_reads_per_pixel"], je["preview_reads_per_pixel"]);
        assert_eq!(
            std::fs::read(r.out.join("heatmap.png")).unwrap(),
            std::fs::read(exact.out.join("heatmap.png")).unwrap(),
            "{name}: PNG must stay byte-identical"
        );
        assert_eq!(
            std::fs::read(r.out.join("heatmap.preview")).unwrap(),
            std::fs::read(exact.out.join("heatmap.preview")).unwrap(),
            "{name}: previews must stay byte-identical"
        );
    }

    // Sanity: the fixture really does put 4 reads on one pixel.
    assert_eq!(je["counts"][0], 4);
}

/// `off` deletes the stream without breaking the file's row structure: every
/// row's encoding is the empty pixel list, so `data_row_offsets` still points
/// at a valid offset in a `size`-byte file.
#[test]
fn off_writes_only_empty_rows() {
    let fx = tmpdir("offrows");
    let fq = fixture(&fx);
    let r = run("offrows", &fq, &["--data-mode", "off"]);
    let j = json(&r);

    let d = data(&r);
    assert_eq!(d.len(), SIZE, "off: one empty varint(0) row per grid row");
    assert!(d.iter().all(|&b| b == 0), "off: every row must be empty");

    let offs: Vec<u64> = j["data_row_offsets"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap())
        .collect();
    assert_eq!(offs.len(), SIZE + 1);
    assert_eq!(offs, (0..=SIZE as u64).collect::<Vec<_>>());
}

/// `capped:N` keeps at most N reads per (pixel, file) and shrinks the stream.
#[test]
fn capped_limits_reads_per_pixel_and_labels_json() {
    let fx = tmpdir("caprows");
    let fq = fixture(&fx);

    let exact = run("capex", &fq, &["--data-mode", "exact"]);
    let capped = run("capcap", &fq, &["--data-mode", "capped:2"]);

    let je = json(&exact);
    assert!(
        je.get("data_complete").is_none(),
        "exact mode must not carry a data-completeness label: {je}"
    );
    assert!(
        je.get("data_reads_per_pixel").is_none(),
        "exact mode must not carry a cap label: {je}"
    );

    let jc = json(&capped);
    assert_eq!(jc["data_complete"], Value::Bool(false));
    assert_eq!(jc["data_reads_per_pixel"], 2);

    let de = data(&exact);
    let dc = data(&capped);
    assert!(
        dc.len() < de.len(),
        "capped:2 must shrink heatmap.data ({} vs {})",
        dc.len(),
        de.len()
    );
    assert_ne!(dc, de, "capped:2 must not reproduce the exact stream");

    // 14 reads, 4 pixels, cap 2 → 4 pixels × 2 reads kept = 8 of 14.
    // Decode `n_stored` per pixel and check the cap held everywhere.
    let mut pos = 0usize;
    let mut total = 0usize;
    let n_rows = SIZE;
    for _ in 0..n_rows {
        let n_pixels = varint(&dc, &mut pos);
        for _ in 0..n_pixels {
            let _x = varint(&dc, &mut pos);
            let n_groups = varint(&dc, &mut pos);
            // groups: varint(file_id) varint(n_reads) varint*n_reads
            let mut stored = 0u64;
            for _ in 0..n_groups {
                let _fid = varint(&dc, &mut pos);
                let n_reads = varint(&dc, &mut pos);
                assert!(n_reads <= 2, "capped:2 stored {n_reads} reads in one group");
                for _ in 0..n_reads {
                    let _ = varint(&dc, &mut pos);
                }
                stored += n_reads;
            }
            assert!(stored <= 2, "capped:2 stored {stored} reads at one pixel");
            total += stored as usize;
        }
    }
    assert_eq!(pos, dc.len(), "capped:2 stream must decode exactly");
    assert_eq!(total, 8, "4 pixels × cap 2 = 8 reads kept");
}

fn varint(buf: &[u8], pos: &mut usize) -> u64 {
    let mut shift = 0u32;
    let mut out = 0u64;
    loop {
        let b = buf[*pos];
        *pos += 1;
        out |= ((b & 0x7f) as u64) << shift;
        if b & 0x80 == 0 {
            return out;
        }
        shift += 7;
    }
}

/// The v1 wrapper carries the same two flags as metrics, so a downstream
/// consumer cannot mistake an approximation for the exact stream.
#[test]
fn non_exact_modes_label_the_v1_wrapper() {
    let fx = tmpdir("wrapper");
    let fq = fixture(&fx);

    for (tag, mode, cap) in [("w_off", "off", 0.0), ("w_cap", "capped:2", 2.0)] {
        let dir = std::env::temp_dir().join(format!("maps_dm_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let jsonp = dir.join("grid.json");
        let r = run(
            tag,
            &fq,
            &[
                "--data-mode",
                mode,
                "--json",
                jsonp.to_str().unwrap(),
            ],
        );
        let rp = r.report.expect("report path recorded");
        let text = std::fs::read_to_string(&rp).unwrap();
        let v: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(v["schema"], "fastsuite.report.v1");
        assert_eq!(v["metrics"]["data_complete"].as_f64().unwrap(), 0.0);
        assert_eq!(v["metrics"]["data_reads_per_pixel"].as_f64().unwrap(), cap);
    }

    // Exact mode: no metrics, no keys — the wrapper shape is untouched.
    let dir = std::env::temp_dir().join(format!("maps_dm_w_ex_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let jsonp = dir.join("grid.json");
    let r = run(
        "w_exact",
        &fq,
        &["--data-mode", "exact", "--json", jsonp.to_str().unwrap()],
    );
    let text = std::fs::read_to_string(r.report.unwrap()).unwrap();
    let v: Value = serde_json::from_str(&text).unwrap();
    assert!(v["metrics"].get("data_complete").is_none(), "{text}");
    assert!(v["metrics"].get("data_reads_per_pixel").is_none(), "{text}");
}

/// The run summary has to say which stream the user got.
#[test]
fn run_summary_names_a_non_exact_mode() {
    let fx = tmpdir("summary");
    let fq = fixture(&fx);

    let off = run("sum_off", &fq, &["--data-mode", "off"]);
    assert!(
        off.stderr.contains("data-mode off") || off.stderr.contains("data mode = off"),
        "off: run summary must name the mode; stderr:\n{}",
        off.stderr
    );
    assert!(
        off.stderr.contains("data_complete=false"),
        "off: run summary must flag the approximation; stderr:\n{}",
        off.stderr
    );

    let cap = run("sum_cap", &fq, &["--data-mode", "capped:3"]);
    assert!(
        cap.stderr.contains("capped:3"),
        "capped: run summary must name the cap; stderr:\n{}",
        cap.stderr
    );

    let ex = run("sum_ex", &fq, &[]);
    assert!(
        !ex.stderr.contains("data_complete")
            && !ex.stdout.contains("data_complete")
            && !ex.stdout.contains("data-mode"),
        "exact (default): the summary must stay silent about data mode"
    );
}
