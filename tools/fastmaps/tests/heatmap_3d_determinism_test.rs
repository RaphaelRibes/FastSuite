//! FM-12 — `heatmap_3d.html` must be reproducible.
//!
//! Two independent sources of non-determinism fed this file:
//!
//!   1. the per-pixel mean quality was an `f64` partial per worker, merged in
//!      *thread-completion* order (`Done` messages race), and
//!   2. the series was emitted in `HashMap` iteration order, which is seeded
//!      per process.
//!
//! So the file differed between two runs of the identical command. It is now
//! built from integer `(sum_q, sum_len)` accumulators and sorted by `(x, y)`.
//! Every other output was already deterministic (integer counts + a sorted
//! encoder); this test pins the last one.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn write_fq(path: &Path, seed: u64, n: usize) {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut s = String::with_capacity(n * 330);
    for i in 0..n {
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
            qual.push((33 + (next() % 42) as u8) as char);
        }
        s.push_str(&format!("@d{i}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(path, s).unwrap();
}

fn run(files: &[String], j: &str, tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("maps_3d_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let mut c = Command::new(bin());
    c.arg("-o").arg(&dir).arg("--force").args(files).args(["-j", j]);
    let o = c.output().unwrap();
    assert!(
        o.status.success(),
        "-j {j} run failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    dir
}

#[test]
fn heatmap_3d_html_is_stable_across_runs_and_threads() {
    let fx = std::env::temp_dir().join(format!("maps_3d_fx_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&fx);
    std::fs::create_dir_all(&fx).unwrap();

    // One mate pair: two files, so -j 2 gives two workers with four files worth
    // of nothing — the file→worker partition is forced (no race over which
    // worker owns which file), which is what makes the *other* outputs
    // comparable run to run. The html comparison still spans the thread split.
    let mut files = Vec::new();
    for (i, seed) in [11u64, 23].iter().enumerate() {
        let p = fx.join(format!("s{i}_R{}.fastq", 1 + i));
        write_fq(&p, *seed, 6_000);
        files.push(p.to_string_lossy().to_string());
    }

    // Two runs per -j: the first pair isolates the per-process `HashMap` seed,
    // the rest isolate thread scheduling.
    let mut dirs = Vec::new();
    for j in ["1", "2"] {
        for rep in ["a", "b"] {
            dirs.push((j, run(&files, j, &format!("j{j}{rep}"))));
        }
    }

    // `heatmap_3d.html` must be identical across *every* run, including across
    // -j: its series is sorted and its quality values are exact integers now.
    let first = std::fs::read(dirs[0].1.join("heatmap_3d.html")).unwrap();
    assert!(!first.is_empty());
    for (j, d) in dirs.iter().skip(1) {
        let other = std::fs::read(d.join("heatmap_3d.html")).unwrap();
        assert!(
            first == other,
            "heatmap_3d.html differs between -j {} run a and -j {j} ({} vs {} bytes)",
            dirs[0].0,
            first.len(),
            other.len()
        );
    }

    // The other outputs must be identical between two runs *at the same -j*
    // (`heatmap.preview` is not -j invariant: the per-(pixel, kind) quota is
    // per worker by design, and that predates this change).
    for pair in dirs.chunks(2) {
        let ((j, a), (_, b)) = (&pair[0], &pair[1]);
        for f in ["heatmap.json", "heatmap.data", "heatmap.preview", "heatmap.png"] {
            let ra = std::fs::read(a.join(f)).unwrap();
            let rb = std::fs::read(b.join(f)).unwrap();
            assert!(ra == rb, "{f} differs between two -j {j} runs");
        }
    }

    let _ = std::fs::remove_dir_all(&fx);
    for (_, d) in &dirs {
        let _ = std::fs::remove_dir_all(d);
    }
}
