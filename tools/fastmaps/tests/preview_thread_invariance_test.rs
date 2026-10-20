//! Review 2026-10-21, item 7 — `heatmap.preview` must not depend on `-j`.
//!
//! `encode_preview_row` documents (and implements) "the first `preview_n` reads
//! per (pixel, kind) by `(file_id, read_num)`". But the *candidate* set is
//! truncated earlier, in the worker, by a quota keyed `((x, y), kind)` — one
//! quota **per worker**. With `-j 1` the single worker sees the files in queue
//! order and keeps whichever file it met first, so the encode-time rule only
//! ever sees that file's reads; with `-j 3` every file's prefix survives and
//! the documented first-N rule wins. Same input, different preview bytes.
//!
//! Fixture: three R1 files, ten 60 bp all-A reads each — every read on pixel
//! (0,0), well past the default `preview_n = 5`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_pv_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn files(dir: &Path) -> Vec<PathBuf> {
    let s = "A".repeat(60);
    (0..3)
        .map(|f| {
            let p = dir.join(format!("f{f}_R1.fastq"));
            let mut body = String::new();
            for i in 0..10 {
                body.push_str(&format!("@f{f}r{i}\n{s}\n+\n{}\n", "I".repeat(60)));
            }
            std::fs::write(&p, body).unwrap();
            p
        })
        .collect()
}

fn run(dir: &Path, tag: &str, threads: &str, inputs: &[PathBuf]) -> PathBuf {
    let out = dir.join(tag);
    let o = Command::new(bin())
        .arg("-q")
        .arg("-o")
        .arg(&out)
        .arg("-j")
        .arg(threads)
        .args(inputs)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "-j {threads} run failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    out
}

#[test]
fn preview_is_identical_at_one_and_three_workers() {
    let dir = tmpdir("inv");
    let inputs = files(&dir);

    let one = run(&dir, "j1", "1", &inputs);
    let three = run(&dir, "j3", "3", &inputs);

    for name in ["heatmap.preview", "heatmap.json", "heatmap.data"] {
        let a = std::fs::read(one.join(name)).unwrap();
        let b = std::fs::read(three.join(name)).unwrap();
        assert_eq!(
            a, b,
            "{name} differs between -j 1 ({} B) and -j 3 ({} B): the preview quota is per worker, so the documented first-N-by-(file_id, read_num) rule never sees the same candidates",
            a.len(),
            b.len()
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}
