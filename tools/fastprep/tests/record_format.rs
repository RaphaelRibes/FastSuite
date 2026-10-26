//! Regression: the destination extension must agree with the record format.
//!
//! fastprep preserves the input record format and infers only the compression
//! from the output extension. At 0.5.0 a FASTA input written to `out.fq`
//! produced a FASTA file with a FASTQ name (exit 0); the next tool in the
//! chain then refused it with a quality-length error that pointed at the wrong
//! step. CLI-SPEC §4 says named paths use the extension.

use std::io::Write;
use std::process::Command;

fn write(path: &std::path::Path, s: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(s.as_bytes()).unwrap();
}

fn dir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep-fmt-{tag}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fastprep"))
        .args(args)
        .output()
        .expect("fastprep ran")
}

#[test]
fn fasta_input_refuses_a_fastq_named_output() {
    let d = dir("fasta-to-fq");
    let inp = d.join("in.fa");
    let out = d.join("out.fq");
    write(&inp, "@r0\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n");
    // now a real FASTA input
    write(&inp, ">r0\nACGTACGTACGTACGTACGT\n");

    let o = run(&[
        "-1",
        inp.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "-q",
    ]);
    assert!(
        !o.status.success(),
        "a FASTA input named .fq must be refused, not silently written"
    );
    let msg = String::from_utf8_lossy(&o.stderr);
    assert!(
        msg.contains("FASTA") && msg.contains("FASTQ"),
        "the refusal must name both formats, got: {msg}"
    );
    assert!(!out.exists(), "nothing may be committed on the refusal");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn fasta_input_with_a_fasta_name_is_kept_as_fasta() {
    let d = dir("fasta-to-fa");
    let inp = d.join("in.fa");
    let out = d.join("out.fa");
    write(&inp, ">r0\nACGTACGTACGTACGTACGTACGTACGTACGT\n");

    let o = run(&[
        "-1",
        inp.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "-q",
    ]);
    assert!(
        o.status.success(),
        "FASTA in / .fa out must work: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.starts_with('>'), "output must still be FASTA: {text}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn fastq_input_refuses_a_fasta_named_output() {
    let d = dir("fq-to-fa");
    let inp = d.join("in.fq");
    let out = d.join("out.fa");
    write(&inp, "@r0\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n");

    let o = run(&[
        "-1",
        inp.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "-q",
    ]);
    assert!(
        !o.status.success(),
        "a FASTQ input named .fa must be refused (the file would not be FASTA)"
    );
    assert!(!out.exists(), "nothing may be committed on the refusal");
    let _ = std::fs::remove_dir_all(&d);
}
