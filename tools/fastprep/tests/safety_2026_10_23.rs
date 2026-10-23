//! Data-loss regressions found in the 2026-10-23 review. Each fails on the
//! pre-fix 0.4.1 binary:
//!
//! 1. `--adapter-scored-retry` set `t5` to the window offset as well as
//!    `t3 = l - p`, so `t5 + t3 == l` and every read the retry rescued was
//!    emitted at length 0, exit 0;
//! 2. `--json PATH` naming `-o PATH` overwrote the FASTQ with the report and
//!    exited 0; `-o x -p x` wrote two streams into one `.part` and failed at
//!    the second rename, after the first commit;
//! 3. `.bz2` / `.xz` outputs were written as PLAIN text under a compressed
//!    name, exit 0;
//! 4. a write error at the final flush was discarded, so a truncated file was
//!    renamed into place and the run reported success;
//! 5. `--pe-correct-min-overlap` above `l1 + l2` underflowed and the run spun
//!    forever on the first pair.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_safety_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_fastq(path: &Path, n: usize, len: usize, adapter: Option<&str>) {
    let mut s = String::new();
    for i in 0..n {
        let insert: String = (0..len)
            .map(|j| b"ACGT"[(i * 31 + j * 7 + i * j * 13) % 4] as char)
            .collect();
        let seq = match adapter {
            Some(a) => format!("{insert}{a}"),
            None => insert,
        };
        let q = "F".repeat(seq.len());
        s.push_str(&format!("@r{i}\n{seq}\n+\n{q}\n"));
    }
    std::fs::write(path, s).unwrap();
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

fn seq_lens(path: &Path) -> Vec<usize> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let lines: Vec<&str> = text.lines().collect();
    lines
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 4 == 1)
        .map(|(_, l)| l.len())
        .collect()
}

/// The retry must trim only the 3' end. Before the fix every rescued read was
/// written at length 0 (the whole read "trimmed"), with exit 0.
#[test]
fn scored_retry_keeps_the_insert() {
    let d = tmpdir("retry");
    let tru = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    let nrich: String = tru
        .chars()
        .enumerate()
        .map(|(i, c)| if i % 4 == 3 && i < 32 { 'N' } else { c })
        .collect();
    let read = format!("{}{}", "TTTTCCCCAAAAGGGG".repeat(4), nrich);
    let mut s = String::new();
    for i in 0..20 {
        s.push_str(&format!("@s{i}\n{read}\n+\n{}\n", "F".repeat(read.len())));
    }
    let inp = d.join("in.fq");
    std::fs::write(&inp, s).unwrap();
    let outp = d.join("out.fq");
    let json = d.join("r.json");

    let (err, rc) = run(&[
        "-1", inp.to_str().unwrap(),
        "-o", outp.to_str().unwrap(),
        "--adapter-preset", "illumina",
        "--adapter-scored-retry", "2",
        "--no-dedup", "--no-mean-q",
        "--json", json.to_str().unwrap(),
    ]);
    assert_eq!(rc, 0, "run failed: {err}");
    let lens = seq_lens(&outp);
    assert_eq!(lens.len(), 20, "every read must be written");
    let ins = 64;
    assert!(
        lens.iter().all(|&l| l == ins),
        "the retry must leave the insert at {ins} bp, got {lens:?}"
    );
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(doc["metrics"]["reads_adapter_trimmed"], 20);
}

/// A path appearing twice in the argv is refused before anything is written.
#[test]
fn colliding_output_paths_are_refused() {
    let d = tmpdir("collide");
    let inp = d.join("in.fq");
    write_fastq(&inp, 5, 60, None);
    let x = d.join("x.fq");
    let xs = x.to_str().unwrap();

    for (tag, args) in [
        ("json==out1", vec!["-1", inp.to_str().unwrap(), "-o", xs, "--json", xs, "--force"]),
        ("metadata==out1", vec!["-1", inp.to_str().unwrap(), "-o", xs, "--metadata", xs, "--force"]),
        ("out2==out1", vec![
            "-1", inp.to_str().unwrap(), "-2", inp.to_str().unwrap(),
            "-o", xs, "-p", xs, "--force",
        ]),
        ("out1==in1", vec!["-1", inp.to_str().unwrap(), "-o", inp.to_str().unwrap(), "--force"]),
    ] {
        let (err, rc) = run(&args);
        assert_ne!(rc, 0, "{tag}: a colliding path must be refused, got rc=0 ({err})");
        assert!(
            err.contains("name the same file"),
            "{tag}: the refusal must name the collision, got: {err}"
        );
        assert!(!x.exists(), "{tag}: nothing may be created on a refusal");
    }
}

/// Output compression that is not implemented must fail loudly instead of
/// writing plain text under a compressed extension.
#[test]
fn unsupported_output_compression_is_refused() {
    let d = tmpdir("compress");
    let inp = d.join("in.fq");
    write_fastq(&inp, 5, 60, None);
    for ext in ["bz2", "xz"] {
        let outp = d.join(format!("out.fq.{ext}"));
        let (err, rc) = run(&[
            "-1", inp.to_str().unwrap(), "-o", outp.to_str().unwrap(), "--no-adapter",
        ]);
        assert_ne!(rc, 0, ".{ext} output must be refused, not written plain ({err})");
        assert!(
            err.contains("not supported"),
            ".{ext}: the refusal must say so, got: {err}"
        );
        assert!(!outp.exists(), ".{ext}: no file may be left behind");
    }
}

/// A short write (RLIMIT_FSIZE) must surface as a failure and leave no output.
#[cfg(unix)]
#[test]
fn a_write_error_is_not_a_success() {
    let d = tmpdir("writer");
    let inp = d.join("in.fq");
    write_fastq(&inp, 400, 60, None);
    let outp = d.join("out.fq");
    let cmd = format!(
        "ulimit -f 1; trap '' XFSZ; exec '{}' -1 '{}' -o '{}' --no-adapter --no-dedup --force -q",
        fastprep(),
        inp.display(),
        outp.display()
    );
    let out = Command::new("bash").arg("-c").arg(&cmd).output().unwrap();
    let rc = out.status.code().unwrap_or(-1);
    assert_ne!(rc, 0, "a truncated write must not exit 0");
    assert!(
        !outp.exists(),
        "the .part must not be renamed over the destination on a failed flush"
    );
}

/// `--pe-correct-min-overlap` above `l1 + l2` used to wrap and spin forever.
#[test]
fn pe_correct_min_overlap_above_the_read_sum_terminates() {
    let d = tmpdir("pecorrect");
    let r1 = d.join("r1.fq");
    let r2 = d.join("r2.fq");
    write_fastq(&r1, 4, 150, None);
    write_fastq(&r2, 4, 150, None);
    let json = d.join("r.json");
    let (err, rc) = run(&[
        "-1", r1.to_str().unwrap(),
        "-2", r2.to_str().unwrap(),
        "-o", d.join("o1.fq").to_str().unwrap(),
        "-p", d.join("o2.fq").to_str().unwrap(),
        "--pe-correct", "2.0",
        "--pe-correct-min-overlap", "400",
        "--no-adapter", "--no-dedup", "--no-mean-q",
        "--json", json.to_str().unwrap(),
    ]);
    assert_eq!(rc, 0, "a too-large overlap floor must abstain, not hang: {err}");
    let n1 = seq_lens(&d.join("o1.fq")).len();
    let n2 = seq_lens(&d.join("o2.fq")).len();
    assert_eq!(n1, n2, "mates must stay in lockstep");
    assert!(n1 > 0, "the reads must be written");
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(
        doc["metrics"]["bases_corrected"], 0,
        "a floor above l1 + l2 must abstain, not correct"
    );
}
