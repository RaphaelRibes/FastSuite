//! A read with zero ACGT bases is unmeasurable, not polyA-rich.
//!
//! Non-ACGT bytes are skipped from the (T+C, T+G) coordinates, so an all-N
//! read lands on the (0,0) pixel exactly like an all-A read. fastsim emits
//! N-runs by design, so before this the suite reported polyA contamination on
//! its own simulated data. The read keeps its pixel; the corner fractions must
//! leave it out and report it as `no_call_reads`.

use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

#[test]
fn all_n_read_is_not_poly_a() {
    let dir = std::env::temp_dir().join(format!("maps_nocall_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // 150 bp reads: one all-A (real polyA), one mid-plane (ACGT repeat, not in
    // any corner), one all-N (no call). size = 151 gives a 4-pixel corner
    // window, so (0,0) is polyA and the mid-plane read is neither corner nor
    // off-diagonal.
    let mid = "ACGT".repeat(37) + "AC";
    let seqs = ["A".repeat(150), mid, "N".repeat(150)];
    let fq = dir.join("mix_R1.fastq");
    let mut body = String::new();
    for (i, s) in seqs.iter().enumerate() {
        body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(s.len())));
    }
    std::fs::write(&fq, body).unwrap();

    let out = dir.join("out");
    let json = dir.join("r.json");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--force")
        .arg("--json")
        .arg(&json)
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let text = std::fs::read_to_string(&json).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    let m = &v["metrics"];

    // One polyA read out of two measurable reads (the all-A and the mid-plane
    // one); the all-N read must not be the third.
    assert!(
        (m["polyA_frac"].as_f64().unwrap() - 0.5).abs() < 1e-9,
        "polyA_frac must be 1/2, got {}",
        m["polyA_frac"]
    );
    assert_eq!(m["no_call_reads"].as_f64().unwrap(), 1.0);
    assert_eq!(m["r1_no_call_reads"].as_f64().unwrap(), 1.0);
    // total_reads still counts every read in the window.
    assert_eq!(m["total_reads"].as_f64().unwrap(), 3.0);

    let _ = std::fs::remove_dir_all(&dir);
}
