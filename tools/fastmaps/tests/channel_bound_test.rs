//! M9 — the worker→spiller channel is bounded by the memory budget.
//!
//! The channel was unbounded: a spiller that fell behind let every worker keep
//! sending, so peak RSS grew with the *file count* — 64 files × 100 k reads
//! peaked at ~190 MB against the documented ~22 MB for 2 files, and
//! `--buffer-mb` could not move it. The capacity is now derived from the buffer
//! budget (`channel_capacity`), which makes the queue's size a function of the
//! budget instead of the number of inputs.
//!
//! Byte-identity is the other half of the contract: bounding the channel must
//! not reorder or drop anything. This drives 64 files through a tiny budget at
//! 8 workers (so the channel is the binding constraint) and compares every
//! output against the same run at `-j 1`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_chan_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// `files` x `reads` deterministic 150 bp reads with varied composition and
/// quality, so previews are kept at many pixels and the buffers really fill.
fn library(dir: &Path, files: usize, reads: usize) -> Vec<PathBuf> {
    let mut out = Vec::with_capacity(files);
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    for f in 0..files {
        let mut body = String::new();
        for i in 0..reads {
            let mut seq = String::with_capacity(150);
            for _ in 0..150 {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                seq.push(match state % 4 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                });
            }
            let qual = "I".repeat(150);
            body.push_str(&format!("@f{f}r{i}\n{seq}\n+\n{qual}\n"));
        }
        let p = dir.join(format!("ch{f}_R1.fastq"));
        std::fs::write(&p, body).unwrap();
        out.push(p);
    }
    out
}

fn run(tag: &str, files: &[PathBuf], extra: &[&str]) -> PathBuf {
    let out = std::env::temp_dir().join(format!("maps_chan_out_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--buffer-mb")
        .arg("1")
        .args(extra)
        .args(files)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "run {tag} failed\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    out
}

#[test]
fn many_files_through_a_tiny_budget_are_byte_identical() {
    let dir = tmpdir("fx");
    let files = library(&dir, 64, 400);

    let wide = run("j8", &files, &["-j", "8"]);
    let narrow = run("j1", &files, &["-j", "1"]);

    for f in [
        "heatmap.data",
        "heatmap.preview",
        "heatmap.json",
        "heatmap.png",
        "heatmap_3d.html",
    ] {
        let a = std::fs::read(wide.join(f)).unwrap();
        let b = std::fs::read(narrow.join(f)).unwrap();
        assert!(
            a == b,
            "{f} differs between -j 8 ({} B) and -j 1 ({} B) with a bounded channel",
            a.len(),
            b.len()
        );
    }

    // Sanity: every read was counted, so the comparison is meaningful.
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(wide.join("heatmap.json")).unwrap()).unwrap();
    assert_eq!(j["total_reads"], 64 * 400);

    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_dir_all(&wide);
    let _ = std::fs::remove_dir_all(&narrow);
}
