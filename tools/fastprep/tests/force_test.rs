//! `--force` (CLI-SPEC §1): without it fastprep must refuse to clobber an
//! existing `-o`/`-p`, name the file, and suggest the flag. With it, the run
//! overwrites as before. fastprep used to overwrite silently (CLI-SPEC §4,
//! gap 2).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_force_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_fastq(dir: &Path) -> PathBuf {
    let p = dir.join("in.fastq");
    std::fs::write(
        &p,
        "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n\
         @r2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
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

#[test]
fn refuses_to_clobber_existing_out1_without_force() {
    let d = tmpdir("single");
    let fq = write_fastq(&d);
    let out = d.join("out.fastq");
    let args = [
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-len",
        "0",
        "--no-adapter",
    ];

    let (err, code) = run(&args);
    assert_eq!(code, 0, "first run must succeed: {err}");
    assert!(out.exists(), "first run must write the output");

    // A sentinel in place of the real output: the refusal must leave it alone.
    std::fs::write(&out, "sentinel\n").unwrap();

    let (err, code) = run(&args);
    assert_ne!(code, 0, "clobbered an existing -o without --force");
    assert!(
        err.contains("--force"),
        "the refusal must suggest --force: {err}"
    );
    assert!(
        err.contains("out.fastq"),
        "the refusal must name the file: {err}"
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "sentinel\n",
        "the refusal must not touch the existing output"
    );

    // --force overwrites; -f behaves identically.
    for flag in ["--force", "-f"] {
        let (err, code) = {
            let mut a = args.to_vec();
            a.push(flag);
            run(&a)
        };
        assert_eq!(code, 0, "{flag} must allow the overwrite: {err}");
        assert_ne!(
            std::fs::read_to_string(&out).unwrap(),
            "sentinel\n",
            "{flag} must actually overwrite"
        );
    }

    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn refuses_to_clobber_existing_out2_without_force() {
    let d = tmpdir("paired");
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    let body = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    std::fs::write(&r1, body).unwrap();
    std::fs::write(&r2, body).unwrap();

    let o1 = d.join("o1.fastq");
    let o2 = d.join("o2.fastq");
    let args = [
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-o",
        o1.to_str().unwrap(),
        "-p",
        o2.to_str().unwrap(),
        "--min-len",
        "0",
        "--no-adapter",
    ];

    let (err, code) = run(&args);
    assert_eq!(code, 0, "first paired run must succeed: {err}");

    // Only -p pre-exists: the refusal must name it, not -o.
    std::fs::remove_file(&o1).unwrap();
    let (err, code) = run(&args);
    assert_ne!(code, 0, "clobbered an existing -p without --force");
    assert!(err.contains("--force"), "must suggest --force: {err}");
    assert!(err.contains("o2.fastq"), "must name -p: {err}");

    let (err, code) = {
        let mut a = args.to_vec();
        a.push("--force");
        run(&a)
    };
    assert_eq!(code, 0, "--force must allow the paired overwrite: {err}");
    assert!(o1.exists() && o2.exists(), "--force must write both mates");

    std::fs::remove_dir_all(&d).unwrap();
}
