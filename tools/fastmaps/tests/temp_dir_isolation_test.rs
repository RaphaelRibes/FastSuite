//! M2 — a stale `--temp-dir` must never be read as this run's index stream.
//!
//! Spill rows live at `<temp-dir>/heatmap_row_NNNNN.tmp` and are only deleted
//! by the render pass, so a SIGKILLed run leaves a full set behind. Reusing the
//! same `--temp-dir` then fed the dead run's records into this run's render
//! pass: 650 reads produced 1940 membership entries and 193 count-zero pixels,
//! exit 0. The fix is a unique per-run subdirectory: leftovers can never be
//! opened, whatever state the reuse left them in.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_tmp_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// 650 reads of 60 bp, spread over a handful of pixels.
fn fixture(dir: &Path) -> PathBuf {
    let fq = dir.join("stale_R1.fastq");
    let mut body = String::new();
    for i in 0..650usize {
        let seq = match i % 4 {
            0 => "A".repeat(60),
            1 => "T".repeat(60),
            2 => "G".repeat(60),
            _ => "C".repeat(60),
        };
        body.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

fn run(out: &Path, fq: &Path, temp: Option<&Path>) -> std::process::Output {
    let mut c = Command::new(bin());
    c.arg("-o").arg(out);
    if let Some(t) = temp {
        c.arg("--temp-dir").arg(t);
    }
    c.arg(fq).output().unwrap()
}

fn entries(dir: &Path) -> Vec<String> {
    match std::fs::read_dir(dir) {
        Ok(rd) => rd
            .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
            .collect(),
        Err(_) => Vec::new(),
    }
}

/// A run that finds a dead run's row files in `--temp-dir` must produce exactly
/// the outputs of a clean run, and must leave the directory clean.
#[test]
fn stale_row_files_are_never_consumed() {
    let dir = tmpdir("stale");
    let fq = fixture(&dir);

    // A clean run in its own output dir, for byte comparison.
    let clean_out = dir.join("clean");
    let o = run(&clean_out, &fq, None);
    assert!(
        o.status.success(),
        "clean run failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    // The reuse: `--temp-dir` already holds what a SIGKILLed run leaves.
    let temp = dir.join("temp");
    std::fs::create_dir_all(&temp).unwrap();
    for y in [0u16, 1, 2, 3] {
        // 10 well-formed 8-byte records for a row nothing will write.
        std::fs::write(
            temp.join(format!("heatmap_row_{y:05}.tmp")),
            vec![0x41u8; 10 * 8],
        )
        .unwrap();
    }
    std::fs::write(temp.join("preview_row_00000.tmp"), b"stale-preview").unwrap();

    let reused_out = dir.join("reused");
    let o = run(&reused_out, &fq, Some(&temp));
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        o.status.success(),
        "a reused --temp-dir must not fail the run\nstderr:\n{stderr}"
    );

    for f in ["heatmap.json", "heatmap.data", "heatmap.preview"] {
        let clean = std::fs::read(clean_out.join(f)).unwrap();
        let reused = std::fs::read(reused_out.join(f)).unwrap();
        assert_eq!(
            clean, reused,
            "{f} differs when --temp-dir holds a dead run's row files: the stale index stream was rendered as this run's"
        );
    }

    // The stale files themselves are untouched (they are not ours), but the
    // per-run subdirectory is gone: no row of ours is left behind.
    let mut left = entries(&temp);
    left.retain(|n| n.starts_with("run-"));
    assert!(
        left.is_empty(),
        "the per-run temp subdirectory must be removed on success, found {left:?}"
    );
}

/// The cleanup also has to happen on an error path: a refused run must not
/// leave a row directory behind for the next one to trip over.
#[test]
fn a_failed_run_removes_its_per_run_directory() {
    let dir = tmpdir("failed");
    let fq = dir.join("bad_R1.fastq");
    // One good record, then a qual/seq mismatch: the worker aborts mid-file.
    std::fs::write(
        &fq,
        format!(
            "@ok\n{}\n+\n{}\n@bad\n{}\n+\n{}\n",
            "A".repeat(60),
            "I".repeat(60),
            "A".repeat(60),
            "I".repeat(10)
        ),
    )
    .unwrap();
    let temp = dir.join("temp");
    let o = run(&dir.join("out"), &fq, Some(&temp));
    assert!(
        !o.status.success(),
        "the mismatch must be refused: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let mut left = entries(&temp);
    left.retain(|n| n.starts_with("run-"));
    assert!(
        left.is_empty(),
        "a failed run must remove its per-run temp subdirectory, found {left:?}"
    );
    for e in entries(&temp) {
        assert!(
            !e.starts_with("heatmap_row_") && !e.starts_with("preview_row_"),
            "a failed run left a spill row at the top of --temp-dir: {e}"
        );
    }
}
