//! A corner extraction must be a whole number of FASTQ records, or not exist.
//!
//! `corners_*.fq` used to be created at its final name before the first read
//! was parsed, and the workers wrote straight into it through a `BufWriter`. A
//! run that failed half way — a truncated FASTQ record in a later file, a spill
//! error — exited 1 with the buffer's tail lost, leaving a file that ends in the
//! middle of a record: 510 matching reads produced 505 records plus a dangling
//! `+`, and nothing on stderr mentioned the file. A reader cannot tell that
//! apart from a short extraction.
//!
//! The extraction is now staged under a temporary name and renamed into place
//! once every read has been seen, so the final name only ever exists complete.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastmaps")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("maps_corner_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// 510 reads of 60 bp, all `A`: every one lands on (0, 60), inside the box, and
/// 510 records are more than one `BufWriter` flush (so a lost tail is visible).
fn good(dir: &Path) -> PathBuf {
    let fq = dir.join("good_R1.fastq");
    let s = "A".repeat(60);
    let mut body = String::new();
    for i in 0..510 {
        body.push_str(&format!("@g{i}\n{s}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&fq, body).unwrap();
    fq
}

/// A record that starts and stops after the `+` line: the mid-run error.
fn truncated(dir: &Path) -> PathBuf {
    let fq = dir.join("bad_R1.fastq");
    std::fs::write(&fq, format!("@t0\n{}\n+\n", "A".repeat(60))).unwrap();
    fq
}

fn run(dir: &Path, out: &Path, files: &[&Path]) -> std::process::Output {
    let mut c = Command::new(bin());
    c.current_dir(dir)
        .arg("-o")
        .arg(out)
        .args(["-j", "1", "--force", "--temp-dir"])
        .arg(dir.join("tmp"))
        .args(["--corners", "0,0:5,65"]);
    for f in files {
        c.arg(f);
    }
    c.output().unwrap()
}

fn names(dir: &Path) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    v.sort();
    v
}

/// The reported case: a good file whose extraction exceeds one buffer flush,
/// then a file that fails mid-run. No `corners_*.fq` may survive — not a partial
/// one, and not an empty one either — and neither may the staging file.
#[test]
fn mid_run_failure_leaves_no_corner_file() {
    let dir = tmpdir("fail");
    let good = good(&dir);
    let bad = truncated(&dir);
    let out = dir.join("out");

    // The work queue pops from the end, so the truncated file is read last.
    let o = run(&dir, &out, &[&bad, &good]);
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(!o.status.success(), "the truncated input must fail the run");
    assert!(
        stderr.contains("bad_R1.fastq") && stderr.contains("truncated"),
        "the error must name the failing input; stderr: {stderr}"
    );

    let got = names(&out);
    assert!(
        got.iter().all(|n| !n.starts_with("corners_")),
        "a failed run must not leave a corner file behind: {got:?}"
    );
    assert!(
        got.iter().all(|n| !n.contains(".part")),
        "a failed run must not leave a staged corner file behind: {got:?}"
    );
    assert!(
        got.is_empty(),
        "a failed run writes no output at all, found {got:?}"
    );
}

/// The success path still produces the extraction, in full: 510 records of four
/// lines each, ending on a newline.
#[test]
fn successful_run_writes_only_whole_records() {
    let dir = tmpdir("ok");
    let good = good(&dir);
    let out = dir.join("out");

    let o = run(&dir, &out, &[&good]);
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let path = out.join("corners_0_0_5_65.fq");
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.ends_with('\n'), "the last record must be terminated");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 510 * 4, "510 reads, four lines each");
    assert_eq!(lines.len() % 4, 0, "records must not be cut");
    for rec in lines.chunks(4) {
        assert!(rec[0].starts_with('@'), "header: {:?}", rec[0]);
        assert_eq!(rec[1].len(), 60, "sequence: {:?}", rec[1]);
        assert_eq!(rec[2], "+", "separator: {:?}", rec[2]);
        assert_eq!(rec[3].len(), 60, "quality: {:?}", rec[3]);
    }
    assert!(
        names(&out).iter().all(|n| !n.contains(".part")),
        "the staging name must be renamed away: {:?}",
        names(&out)
    );
}
