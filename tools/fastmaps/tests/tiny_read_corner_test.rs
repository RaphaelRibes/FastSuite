//! Review 2026-10-21, item 2 — unsigned underflow in the corner pass.
//!
//! `corner_metrics` computed the high corner as `x >= size as usize - 1 - c`.
//! On a grid no wider than the corner window (`size <= 1 + c`, i.e. reads of a
//! couple of bases) that subtraction wraps in release and panics in debug: a
//! 1 bp all-T read reported `polyA_frac 1, polyT_frac 0` at exit 0 — the exact
//! opposite of the truth — instead of placing the read in the polyT corner.

use std::path::PathBuf;
use std::process::Command;

use serde_json::Value;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

#[test]
fn one_base_all_t_read_lands_in_the_polyt_corner() {
    let dir: PathBuf = std::env::temp_dir().join(format!("maps_tiny_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fq = dir.join("tiny_R1.fastq");
    std::fs::write(&fq, "@t\nT\n+\nI\n").unwrap();
    let out = dir.join("out");
    let report = dir.join("r.json");

    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--json")
        .arg(&report)
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        o.status.success(),
        "a 1 bp read must not panic or fail (status {:?})\nstderr:\n{stderr}",
        o.status.code()
    );
    assert!(
        !stderr.contains("panicked"),
        "the corner pass must not panic on a 1 bp read; stderr:\n{stderr}"
    );

    let text = std::fs::read_to_string(&report).unwrap();
    let v: Value = serde_json::from_str(&text).unwrap();
    let m = &v["metrics"];
    assert!(
        (m["polyT_frac"].as_f64().unwrap() - 1.0).abs() < 1e-9,
        "an all-T read is polyT, got polyT_frac = {} ({m})",
        m["polyT_frac"]
    );
    assert_eq!(m["total_reads"].as_f64().unwrap(), 1.0);
    // size = 2, corner window = 2: every pixel of a plane narrower than the
    // window is inside every corner, so polyA counts the same read too. That
    // is the published window (`max(size/50, 2)`), unchanged by this fix.
    assert!((m["polyA_frac"].as_f64().unwrap() - 1.0).abs() < 1e-9, "{m}");

    let _ = std::fs::remove_dir_all(&dir);
}
