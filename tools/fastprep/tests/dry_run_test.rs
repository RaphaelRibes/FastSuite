//! `--dry-run` (CLI-SPEC §1): read the input, compute every counter, write
//! nothing — no output FASTQ, no metadata, no --json file. The summary still
//! goes to stderr, and the counters must match a real run exactly.
//!
//! fastprep had no stats-only mode (CLI-SPEC §4, gap 2), so there was no way
//! to check what a trimming configuration would do to a sample before
//! spending the I/O on it.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_dryrun_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// One clean read, one read with a long 5′ polyA tail that the default
/// polyX scan drops, one that survives with a 3′ polyT trim. Enough to make
/// every counter non-trivial.
fn write_fastq(dir: &Path) -> PathBuf {
    let p = dir.join("in.fastq");
    // Quality strings are derived from the sequence length, so a fixture
    // edit can never produce a seq/qual mismatch.
    let rec = |id: &str, seq: &str| {
        format!("@{id}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
    };
    let body = [
        rec("clean", "ACGTACGTACGTACGTACGTACGTACGTACGT"),
        rec("polyA_head", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTACGTACGT"),
        rec(
            "polyT_tail",
            "ACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
    ]
    .join("");
    std::fs::write(&p, body).unwrap();
    p
}

fn run(args: &[&str]) -> (String, i32) {
    let out = Command::new(fastprep())
        .args(args)
        .output()
        .expect("spawn fastprep");
    (
        String::from_utf8_lossy(&out.stderr).to_string(),
        out.status.code().unwrap_or(-1),
    )
}

/// Pull one counter out of the stderr summary table.
fn counter(stderr: &str, label: &str) -> u64 {
    stderr
        .lines()
        .find(|l| l.contains(label))
        .and_then(|l| {
            l.split_whitespace()
                .find_map(|w| w.replace(',', "").parse::<u64>().ok())
        })
        .unwrap_or_else(|| panic!("counter {label:?} not found in:\n{stderr}"))
}

#[test]
fn dry_run_writes_nothing_but_reports_the_same_counters() {
    let d = tmpdir("single");
    let fq = write_fastq(&d);

    // Baseline: a real run, so the dry run has something to be equal to.
    let real_out = d.join("real.fastq");
    let base = [
        "-1",
        fq.to_str().unwrap(),
        "--min-len",
        "0",
        "--no-adapter",
    ];
    let (real_err, code) = {
        let mut a = base.to_vec();
        a.extend(["-o", real_out.to_str().unwrap()]);
        run(&a)
    };
    assert_eq!(code, 0, "the real run must succeed: {real_err}");

    // The dry run points at paths that do not exist yet.
    let dry_out = d.join("dry.fastq");
    let dry_meta = d.join("dry.pxmt");
    let dry_json = d.join("dry.json");
    let (dry_err, code) = {
        let mut a = base.to_vec();
        a.extend([
            "--dry-run",
            "-o",
            dry_out.to_str().unwrap(),
            "--metadata",
            dry_meta.to_str().unwrap(),
            "--json",
            dry_json.to_str().unwrap(),
        ]);
        run(&a)
    };
    assert_eq!(code, 0, "--dry-run must succeed: {dry_err}");

    // Nothing was written. Not the FASTQ, not the metadata, not the report.
    for p in [&dry_out, &dry_meta, &dry_json] {
        assert!(!p.exists(), "--dry-run must not write {}", p.display());
    }

    // The summary is still printed, and every counter matches the real run.
    assert!(
        dry_err.contains("reads in"),
        "--dry-run must still report the summary: {dry_err}"
    );
    for label in ["reads in", "kept:", "dropped:", "bases in"] {
        assert_eq!(
            counter(&dry_err, label),
            counter(&real_err, label),
            "--dry-run's {label} must equal the real run's"
        );
    }
    assert!(
        dry_err.contains("elapsed"),
        "--dry-run must still report the runtime: {dry_err}"
    );

    // A dry run over a path that already exists must not trip --force: it
    // writes nothing, so it can clobber nothing.
    std::fs::write(&dry_out, "sentinel\n").unwrap();
    let (err, code) = {
        let mut a = base.to_vec();
        a.extend(["--dry-run", "-o", dry_out.to_str().unwrap()]);
        run(&a)
    };
    assert_eq!(code, 0, "--dry-run must not need --force: {err}");
    assert_eq!(
        std::fs::read_to_string(&dry_out).unwrap(),
        "sentinel\n",
        "--dry-run must leave an existing file untouched"
    );

    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn dry_run_paired_computes_every_pair_counter_without_writing() {
    let d = tmpdir("paired");
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    std::fs::write(&r1, "@p1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
                        @p2\nACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n").unwrap();
    std::fs::write(&r2, "@p1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
                        @p2\nACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIII\n").unwrap();

    let o1 = d.join("o1.fastq");
    let o2 = d.join("o2.fastq");
    let s1 = d.join("s1.fastq");
    let (err, code) = run(&[
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-o",
        o1.to_str().unwrap(),
        "-p",
        o2.to_str().unwrap(),
        "--singletons1",
        s1.to_str().unwrap(),
        "--dry-run",
        "--min-len",
        "0",
        "--no-adapter",
    ]);
    assert_eq!(code, 0, "paired --dry-run must succeed: {err}");

    for p in [&o1, &o2, &s1] {
        assert!(!p.exists(), "--dry-run must not write {}", p.display());
    }

    // Both mates' reports plus the pairs block are computed.
    assert!(err.contains("[R1]"), "R1 report missing:\n{err}");
    assert!(err.contains("[R2]"), "R2 report missing:\n{err}");
    for label in ["both kept:", "both dropped:", "singletons written (R1):"] {
        assert!(
            err.contains(label),
            "paired counter {label:?} missing from the dry run:\n{err}"
        );
    }
    assert_eq!(counter(&err, "both kept:"), 2, "both pairs must be kept");
    assert_eq!(counter(&err, "reads in"), 2, "R1 reads in");

    std::fs::remove_dir_all(&d).unwrap();
}
