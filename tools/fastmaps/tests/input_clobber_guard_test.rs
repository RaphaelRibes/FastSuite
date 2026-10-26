//! An output target must never be one of the inputs.
//!
//! `fastmaps -o maps b.fq --json b.fq --force` exited 0 with `b.fq` replaced by
//! the 753-byte `fastsuite.report.v1` document. The guards compared this run's
//! outputs with each other (`first_collision`) and with "does this file already
//! exist", but never with `args.files` — and `--force` answered the last
//! question for it. Nothing re-reads an input after the report is written, so
//! the run looks successful while the FASTQ is gone.
//!
//! Every spelling of the same file has to be caught before the first
//! `File::create`: a plain path, a `./` spelling, a hardlink (another name for
//! the same inode), a fixed output name inside `-o`, and `-o` itself.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_in_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A three-read FASTQ whose bytes are recognisable in an assertion message.
fn fixture(dir: &Path) -> PathBuf {
    let fq = dir.join("b_R1.fastq");
    let mut body = String::new();
    for (h, s) in [
        ("r0", "A".repeat(60)),
        ("r1", "ACGT".repeat(15)),
        ("r2", "N".repeat(60)),
    ] {
        body.push_str(&format!("@{h}\n{s}\n+\n{}\n", "I".repeat(s.len())));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

/// Re-run the reported command line and require a refusal that leaves the input
/// byte for byte intact and writes nothing.
fn assert_refused(tag: &str, dir: &Path, args: &[&str], input: &Path, before: &[u8]) {
    let o = Command::new(bin())
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        !o.status.success(),
        "{tag}: the run must be refused\nargs: {args:?}\nstdout: {}\nstderr: {stderr}",
        String::from_utf8_lossy(&o.stdout)
    );
    assert_ne!(
        o.status.code(),
        Some(101),
        "{tag}: must be a refusal, not a panic\nstderr: {stderr}"
    );
    assert!(
        stderr.contains("input"),
        "{tag}: the refusal must say an output names an input; stderr: {stderr}"
    );
    assert_eq!(
        std::fs::read(input).unwrap(),
        before,
        "{tag}: {} was modified by a refused run",
        input.display()
    );
    assert!(
        !dir.join("maps").exists(),
        "{tag}: a refused run must not create its output directory"
    );
}

/// The reported command line: `--json` and the input are the same plain path.
#[test]
fn json_target_naming_an_input_is_refused() {
    let dir = tmpdir("plain");
    let fq = fixture(&dir);
    let before = std::fs::read(&fq).unwrap();

    assert_refused(
        "plain",
        &dir,
        &[
            "-o", "maps", "--json", "b_R1.fastq", "--force", "b_R1.fastq",
        ],
        &fq,
        &before,
    );
}

/// `./b_R1.fastq` names the same file as `b_R1.fastq`, and so does a `..`
/// detour through a directory that exists.
#[test]
fn dot_spelled_json_target_naming_an_input_is_refused() {
    let dir = tmpdir("dot");
    let fq = fixture(&dir);
    let before = std::fs::read(&fq).unwrap();
    // The `..` spelling needs its intermediate directory to exist, or the run is
    // refused earlier for a missing report parent instead.
    std::fs::create_dir(dir.join("x")).unwrap();

    for target in ["./b_R1.fastq", "./x/../b_R1.fastq"] {
        assert_refused(
            "dot",
            &dir,
            &["-o", "maps", "--json", target, "--force", "b_R1.fastq"],
            &fq,
            &before,
        );
    }
}

/// A hardlink is the same file under another name: rewriting one destroys the
/// data behind the other, so the path comparison alone is not enough.
#[test]
fn json_target_that_is_a_hardlink_to_an_input_is_refused() {
    let dir = tmpdir("hardlink");
    let fq = fixture(&dir);
    let before = std::fs::read(&fq).unwrap();
    let link = dir.join("link.fastq");
    std::fs::hard_link(&fq, &link).unwrap();

    let o = Command::new(bin())
        .current_dir(&dir)
        .args(["-o", "maps", "--json", "link.fastq", "--force", "b_R1.fastq"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        !o.status.success(),
        "a hardlink to an input must be refused\nstderr: {stderr}"
    );
    assert!(
        stderr.contains("input"),
        "the refusal must say an output names an input; stderr: {stderr}"
    );
    assert_eq!(
        std::fs::read(&fq).unwrap(),
        before,
        "the input was modified through its hardlink"
    );
    assert_eq!(std::fs::read(&link).unwrap(), before, "the hardlink was modified");
}

/// A symlink to an input is the same file too (`canonicalize`, since the names
/// differ).
#[cfg(unix)]
#[test]
fn symlinked_json_target_naming_an_input_is_refused() {
    let dir = tmpdir("symlink");
    let fq = fixture(&dir);
    let before = std::fs::read(&fq).unwrap();
    let link = dir.join("alias.fastq");
    std::os::unix::fs::symlink("b_R1.fastq", &link).unwrap();

    let o = Command::new(bin())
        .current_dir(&dir)
        .args(["-o", "maps", "--json", "alias.fastq", "--force", "b_R1.fastq"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(!o.status.success(), "stderr: {stderr}");
    assert!(
        stderr.contains("input"),
        "the refusal must say an output names an input; stderr: {stderr}"
    );
    assert_eq!(std::fs::read(&fq).unwrap(), before, "the input was modified");
}

/// The fixed outputs are outputs too: an input that happens to sit at
/// `<output-dir>/heatmap.json` would be replaced by the grid.
#[test]
fn input_at_a_fixed_output_name_is_refused() {
    let dir = tmpdir("fixed");
    fixture(&dir);
    let out = dir.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let victim = out.join("heatmap.json");
    std::fs::write(&victim, b"@r0\nAAAA\n+\nIIII\n").unwrap();
    let before = std::fs::read(&victim).unwrap();

    let o = Command::new(bin())
        .current_dir(&dir)
        .args(["-o", "out", "--force", "out/heatmap.json"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        !o.status.success(),
        "the grid must not overwrite an input\nstderr: {stderr}"
    );
    assert!(
        stderr.contains("heatmap.json") && stderr.contains("input"),
        "the refusal must name the output and the input; stderr: {stderr}"
    );
    assert_eq!(
        std::fs::read(&victim).unwrap(),
        before,
        "heatmap.json was overwritten by the grid"
    );
}

/// `-o` may not be an input either: `create_dir_all` on an existing file used to
/// be the only thing standing between the run and the input's own name.
#[test]
fn output_dir_that_is_an_input_is_refused() {
    let dir = tmpdir("outdir");
    let fq = fixture(&dir);
    let before = std::fs::read(&fq).unwrap();

    let o = Command::new(bin())
        .current_dir(&dir)
        .args(["-o", "b_R1.fastq", "--force", "b_R1.fastq"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(!o.status.success(), "stderr: {stderr}");
    assert!(
        stderr.contains("input"),
        "the refusal must say the output directory is an input; stderr: {stderr}"
    );
    assert_eq!(std::fs::read(&fq).unwrap(), before, "-o touched the input");
}

/// The guard must not over-block: a report beside the input still runs.
#[test]
fn a_report_next_to_the_input_still_runs() {
    let dir = tmpdir("ok");
    let fq = fixture(&dir);

    let o = Command::new(bin())
        .current_dir(&dir)
        .args(["-o", "maps", "--json", "r.json", "--force", "b_R1.fastq"])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "a distinct report target must still run\nstderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(dir.join("r.json").is_file());
    assert!(dir.join("maps/heatmap.json").is_file());
    assert!(fq.is_file(), "the input must still be there");
}
