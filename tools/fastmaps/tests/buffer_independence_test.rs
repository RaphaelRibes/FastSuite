//! FM-2 / FM-3 / FM-9 — the spill path must be invisible in the outputs.
//!
//! `--buffer-mb` is a memory knob, not an output knob. This test drives the
//! same input through three budgets that take three different code paths:
//! `--buffer-mb 1` (one worker, 1 MiB share, so the index buffer spills several
//! times and rows get reopened by the spiller), `--buffer-mb 2` (two spills)
//! and `--buffer-mb 64` (nothing ever spills).
//!
//! All four gated outputs must hash the same. It also exercises the preview
//! flush while spilling (previews are flushed on their own, much smaller,
//! threshold), so a budget-derived preview threshold that split previews
//! differently would show up here.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

/// Deterministic 150 bp reads with a spread of compositions (so the grid has
/// many pixels and the preview buffer actually fills) and varied qualities.
fn fixture(dir: &Path, n_reads: usize) -> String {
    let fq = dir.join("buf_R1.fastq");
    let mut s = String::with_capacity(n_reads * 320);
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for i in 0..n_reads {
        let mut seq = String::with_capacity(150);
        for _ in 0..150 {
            seq.push(match next() % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            });
        }
        let mut qual = String::with_capacity(150);
        for _ in 0..150 {
            // '!'..'J' → Q0..Q41, so mean quality varies per read.
            qual.push((33 + (next() % 42) as u8) as char);
        }
        s.push_str(&format!("@b{i}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(&fq, s).unwrap();
    fq.to_string_lossy().to_string()
}

fn run(tag: &str, fq: &str, extra: &[&str]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("maps_buf_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let o = Command::new(bin())
        .arg("-o")
        .arg(&dir)
        .arg("--force")
        .arg(fq)
        .args(extra)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "run {tag} failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    dir
}

#[test]
fn output_is_identical_across_buffer_budgets() {
    let fx = std::env::temp_dir().join(format!("maps_buf_fx_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&fx);
    std::fs::create_dir_all(&fx).unwrap();
    // 70 000 reads × 8 B spilled = 560 KB per record stream; with one worker
    // the 1 MiB share holds 65 536 of them, so this spills at least once.
    let fq = fixture(&fx, 70_000);

    let a = run("mb1", &fq, &["--buffer-mb", "1", "-j", "1"]);
    let b = run("mb2", &fq, &["--buffer-mb", "2", "-j", "1"]);
    let c = run("mb64", &fq, &["--buffer-mb", "64", "-j", "1"]);

    for f in ["heatmap.json", "heatmap.png", "heatmap.data", "heatmap.preview"] {
        let ra = std::fs::read(a.join(f)).unwrap();
        for (tag, dir) in [("--buffer-mb 2", &b), ("--buffer-mb 64", &c)] {
            let rb = std::fs::read(dir.join(f)).unwrap();
            assert!(
                ra == rb,
                "{f} differs between --buffer-mb 1 and {tag} ({} vs {} bytes)",
                ra.len(),
                rb.len()
            );
        }
    }

    // Sanity: the grid is not degenerate, so the comparison is meaningful.
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(a.join("heatmap.json")).unwrap()).unwrap();
    assert_eq!(j["total_reads"], 70_000);
    let occupied = j["counts"].as_array().unwrap().iter().filter(|c| c.as_u64().unwrap() > 0).count();
    assert!(occupied > 500, "fixture should populate many pixels, got {occupied}");

    let _ = std::fs::remove_dir_all(&fx);
    let _ = std::fs::remove_dir_all(&a);
    let _ = std::fs::remove_dir_all(&b);
    let _ = std::fs::remove_dir_all(&c);
}
