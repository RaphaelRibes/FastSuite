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

/// The singleton streams are outputs too: a FASTA input written into a `.fq`
/// singleton path is the same footgun `check_output_format` closed for `-o`/`-p`
/// (a FASTA file with a FASTQ name, exit 0), and it skipped the check.
///
/// The mate that survives is the one the other's `--min-len` floor rejects, so
/// each direction writes a singleton: R1 in the first run, R2 in the second.
#[test]
fn singletons_refuse_a_format_contradicting_their_extension() {
    let d = dir("singleton-fmt");
    let long = d.join("long.fa");
    let short = d.join("short.fa");
    write(&long, &format!(">long\n{}\n", "GTAC".repeat(38)));
    write(&short, &format!(">short\n{}\n", "GTAC".repeat(15)));

    // (in1, in2, singleton flag, singleton path, tag)
    let cases = [
        (&long, &short, "--singletons1", d.join("s1.fq"), "s1"),
        (&short, &long, "--singletons2", d.join("s2.fq"), "s2"),
    ];
    for (in1, in2, flag, singleton, tag) in cases {
        let out1 = d.join(format!("o1-{tag}.fa"));
        let out2 = d.join(format!("o2-{tag}.fa"));
        let o = run(&[
            "-1",
            in1.to_str().unwrap(),
            "-2",
            in2.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            "-p",
            out2.to_str().unwrap(),
            flag,
            singleton.to_str().unwrap(),
            "--no-adapter",
            "--no-mean-q",
            "--min-len",
            "100",
            "-q",
        ]);
        assert!(
            !o.status.success(),
            "{flag}: a FASTA input named .fq must be refused, not silently written"
        );
        let msg = String::from_utf8_lossy(&o.stderr);
        assert!(
            msg.contains("FASTA") && msg.contains("FASTQ"),
            "{flag}: the refusal must name both formats, got: {msg}"
        );
        assert!(
            msg.contains(singleton.file_name().unwrap().to_str().unwrap()),
            "{flag}: the refusal must name the singleton path, got: {msg}"
        );
        assert!(
            !singleton.exists(),
            "{flag}: nothing may be committed on the refusal"
        );
        assert!(
            !out1.exists() && !out2.exists(),
            "{flag}: nor the main outputs"
        );
    }
    let _ = std::fs::remove_dir_all(&d);
}
