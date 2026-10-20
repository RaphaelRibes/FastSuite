//! P2 — the `--force` guard must cover every file the run writes.
//!
//! The guard used to list `heatmap.json`, `heatmap.png`, `heatmap.data` and
//! `heatmap_3d.html` only, so an existing `heatmap.svg`, `heatmap.preview`, a
//! `--corners` extraction or a `--json` target inside the output dir was
//! overwritten without `--force` — the corners file before the guard could even
//! fire. Each output below is pre-created with sentinel bytes and the same
//! command line is re-run without `--force`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_guard_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fixture(dir: &Path) -> PathBuf {
    let fq = dir.join("reads_R1.fastq");
    let s = "A".repeat(60);
    std::fs::write(&fq, format!("@r1\n{s}\n+\n{}\n", "I".repeat(60))).unwrap();
    fq
}

fn run(out: &Path, fq: &Path, extra: &[String], force: bool) -> std::process::Output {
    let mut c = Command::new(bin());
    c.arg("-o").arg(out);
    if force {
        c.arg("--force");
    }
    c.arg(fq).args(extra);
    c.output().unwrap()
}

/// Pre-create one output file, re-run the same command line without `--force`,
/// and require a non-zero exit that names the file, with the sentinel intact.
fn assert_guarded(tag: &str, victim: &str, extra: &dyn Fn(&Path) -> Vec<String>) {
    let dir = tmpdir(tag);
    let fq = fixture(&dir);
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();

    let sentinel: &[u8] = b"sentinel: must survive a refused run\n";
    let path = out.join(victim);
    std::fs::write(&path, sentinel).unwrap();

    let o = run(&out, &fq, &extra(&out), false);
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        !o.status.success(),
        "{tag}: run must refuse without --force\nstdout: {}\nstderr: {stderr}",
        String::from_utf8_lossy(&o.stdout)
    );
    assert!(
        stderr.contains(victim),
        "{tag}: the refusal must name {victim}; stderr:\n{stderr}"
    );
    assert_eq!(
        std::fs::read(&path).unwrap(),
        sentinel,
        "{tag}: {victim} was overwritten despite the refusal"
    );

    // And the guard must not over-block: --force overwrites it.
    let o = run(&out, &fq, &extra(&out), true);
    assert!(
        o.status.success(),
        "{tag}: --force must still overwrite; stderr:\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert_ne!(
        std::fs::read(&path).unwrap(),
        sentinel,
        "{tag}: --force must actually rewrite {victim}"
    );
}

#[test]
fn guard_covers_svg_preview_and_corners() {
    assert_guarded("svg", "heatmap.svg", &|_| vec!["--svg".into()]);
    assert_guarded("preview", "heatmap.preview", &|_| Vec::new());
    assert_guarded("corners", "corners_0_0_5_5.fq", &|_| {
        vec!["--corners".into(), "0,0:5,5".into()]
    });
}

#[test]
fn guard_covers_the_json_target_and_its_wrapper() {
    let json_args = |out: &Path| {
        vec![
            "--json".into(),
            out.join("grid.json").to_string_lossy().into_owned(),
        ]
    };
    assert_guarded("json", "grid.json", &json_args);
    assert_guarded("wrapper", "grid.report.json", &json_args);
}
