//! M1 — two outputs of one run must never resolve to the same file.
//! M5 — `--json` inside a fresh `-o` must be recognised across `./` spellings.
//!
//! `--json` used to be checked only for "does this file already exist?", so
//! `-o oc --json oc/heatmap.data` ran happily and replaced the 8-byte-per-read
//! index stream with 1 KB of report JSON — exit 0, silently. The check needs to
//! be structural: no two targets of one run may name the same file, whatever
//! the pre-existing-file state is and whether or not `--force` was passed.
//!
//! M5 is the same comparison in the other direction: `Path::starts_with` is
//! component-wise and counts `CurDir`, so `-o outB --json ./outB/q/r.json` and
//! `-o ./outC --json outC/q/r.json` were both refused as "directory does not
//! exist" even though the report lands inside the output dir this run creates.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_tgt_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fixture(dir: &Path) -> PathBuf {
    let fq = dir.join("a_R1.fastq");
    let s = "A".repeat(60);
    let mut body = String::new();
    for i in 0..3 {
        body.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

fn stderr_of(o: &std::process::Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

/// Every fixed output name, asked for as the report target.
#[test]
fn json_target_colliding_with_a_fixed_output_is_refused() {
    for name in [
        "heatmap.data",
        "heatmap.json",
        "heatmap.preview",
        "heatmap.png",
        "heatmap_3d.html",
    ] {
        let dir = tmpdir(&format!("fixed_{}", name.replace('.', "_")));
        let fq = fixture(&dir);
        let out = dir.join("oc");
        let target = out.join(name);

        let o = Command::new(bin())
            .arg("-o")
            .arg(&out)
            .arg("--json")
            .arg(&target)
            .arg(&fq)
            .output()
            .unwrap();
        let stderr = stderr_of(&o);

        assert!(
            !o.status.success(),
            "--json {name} must be refused, not run\nstderr:\n{stderr}"
        );
        assert_ne!(
            o.status.code(),
            Some(101),
            "--json {name}: must be a refusal, not a panic\nstderr:\n{stderr}"
        );
        assert!(
            stderr.contains("collision") && stderr.contains(name),
            "--json {name}: the refusal must name the collision\nstderr:\n{stderr}"
        );
        // Refused before the first create: nothing was written at all.
        for f in [
            "heatmap.data",
            "heatmap.json",
            "heatmap.preview",
            "heatmap.png",
            "heatmap_3d.html",
        ] {
            assert!(
                !out.join(f).exists(),
                "--json {name}: a refused run must write nothing, found {f}"
            );
        }
    }
}

/// The SVG is only an output on `--svg`, so it is only a collision then.
#[test]
fn json_target_colliding_with_the_svg_is_refused() {
    let dir = tmpdir("svg");
    let fq = fixture(&dir);
    let out = dir.join("oc");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--svg")
        .arg("--json")
        .arg(out.join("heatmap.svg"))
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = stderr_of(&o);
    assert!(!o.status.success(), "stderr:\n{stderr}");
    assert!(
        stderr.contains("collision") && stderr.contains("heatmap.svg"),
        "stderr:\n{stderr}"
    );
    assert!(!out.join("heatmap.svg").exists());
}

/// Two `--corners` specs that normalise to the same box write the same file:
/// the second extraction would silently replace the first.
#[test]
fn two_corners_specs_naming_one_file_are_refused() {
    let dir = tmpdir("corners");
    let fq = fixture(&dir);
    let out = dir.join("oc");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--corners")
        .arg("0,0:5,5")
        .arg("--corners")
        .arg("5,5:0,0")
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = stderr_of(&o);
    assert!(!o.status.success(), "stderr:\n{stderr}");
    assert!(
        stderr.contains("collision") && stderr.contains("corners_0_0_5_5.fq"),
        "stderr:\n{stderr}"
    );
    assert!(!out.join("corners_0_0_5_5.fq").exists());
}

/// `--force` means "overwrite what is already there", not "let one output of
/// this run eat another".
#[test]
fn force_does_not_bypass_the_collision_check() {
    let dir = tmpdir("force");
    let fq = fixture(&dir);
    let out = dir.join("oc");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--force")
        .arg("--json")
        .arg(out.join("heatmap.data"))
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = stderr_of(&o);
    assert!(!o.status.success(), "stderr:\n{stderr}");
    assert!(stderr.contains("collision"), "stderr:\n{stderr}");
    assert!(
        !out.join("heatmap.data").exists(),
        "heatmap.data must not exist after a refused run"
    );
}

/// Distinct targets are still fine, including a report inside `-o`.
#[test]
fn distinct_targets_still_run() {
    let dir = tmpdir("distinct");
    let fq = fixture(&dir);
    let out = dir.join("oc");
    let report = out.join("qc").join("r.json");
    let o = Command::new(bin())
        .arg("-o")
        .arg(&out)
        .arg("--json")
        .arg(&report)
        .arg(&fq)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr:\n{}",
        stderr_of(&o)
    );
    assert!(out.join("heatmap.data").is_file());
    assert!(report.is_file());
}

/// M5: a `./` on either side must not hide that the report is inside `-o`.
#[test]
fn json_inside_a_fresh_output_dir_accepts_dot_prefixes() {
    // (tag, -o spelling, --json spelling relative to the run directory)
    for (tag, out_arg, json_arg) in [
        ("json_dot", "outB", "./outB/q/r.json"),
        ("out_dot", "./outC", "outC/q/r.json"),
        ("both_dot", "./outD", "./outD/q/r.json"),
        ("dotdot", "outE", "./outE/q/../q/r.json"),
    ] {
        let dir = tmpdir(tag);
        let fq = fixture(&dir);
        let o = Command::new(bin())
            .current_dir(&dir)
            .arg("-o")
            .arg(out_arg)
            .arg("--json")
            .arg(json_arg)
            .arg(&fq)
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "{tag}: -o {out_arg} --json {json_arg} must be accepted\nstderr:\n{}",
            stderr_of(&o)
        );
        let report = dir.join(json_arg);
        assert!(
            report.is_file(),
            "{tag}: the report must be written at {json_arg}"
        );
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&report).unwrap()).unwrap();
        assert_eq!(v["schema"], "fastsuite.report.v1");
        assert!(dir.join(out_arg.trim_start_matches("./")).join("heatmap.json").is_file());
    }
}

/// Normalising the comparison must not turn the safety check into a free pass:
/// a report parent outside `-o` that does not exist is still refused.
#[test]
fn a_missing_report_parent_outside_the_output_dir_is_still_refused() {
    let dir = tmpdir("outside_missing");
    let fq = fixture(&dir);
    let o = Command::new(bin())
        .current_dir(&dir)
        .arg("-o")
        .arg("outF")
        .arg("--json")
        .arg("../missing/r.json")
        .arg(&fq)
        .output()
        .unwrap();
    let stderr = stderr_of(&o);
    assert!(
        !o.status.success(),
        "--json outside -o with a missing parent must stay refused\nstderr:\n{stderr}"
    );
    assert!(
        stderr.contains("does not exist"),
        "the refusal must say the directory is missing; stderr:\n{stderr}"
    );
}
