//! Regression: the total-drop refusal must fire BEFORE the output is
//! committed.
//!
//! At 0.5.0 the refusal ran after the atomic rename, so a run that dropped
//! every read still replaced the destination with a 0-byte file and then said
//! "nothing was written to the output file" — with `--force` that destroyed a
//! good previous output (reproduced: `keep.fq` 10 B -> 0 B, exit 1).

use std::io::Write;
use std::process::Command;

fn tmpdir() -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep-totaldrop-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write(path: &std::path::Path, s: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(s.as_bytes()).unwrap();
}

const GOOD: &str = "@keep\nACGT\n+\nIIII\n";

#[test]
fn a_total_drop_does_not_destroy_an_existing_output() {
    let d = tmpdir();
    let input = d.join("drop.fq");
    let dest = d.join("keep.fq");
    // One read so poor that --mean-q 60 drops it.
    write(&input, "@drop\nACGTACGTACGT\n+\n!!!!!!!!!!!!\n");
    write(&dest, GOOD);

    let out = Command::new(env!("CARGO_BIN_EXE_fastprep"))
        .args([
            "-1",
            input.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
            "--mean-q",
            "60",
            "--no-adapter",
            "--force",
            "-q",
        ])
        .output()
        .expect("fastprep ran");

    assert!(
        !out.status.success(),
        "a total drop must exit non-zero (stderr: {})",
        String::from_utf8_lossy(&out.stderr)
    );
    let after = std::fs::read_to_string(&dest).unwrap_or_default();
    assert_eq!(
        after, GOOD,
        "the destination must be left untouched when the run refuses"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn allow_total_drop_still_writes_the_empty_output() {
    let d = tmpdir();
    let input = d.join("drop.fq");
    let dest = d.join("out.fq");
    write(&input, "@drop\nACGTACGTACGT\n+\n!!!!!!!!!!!!\n");

    let out = Command::new(env!("CARGO_BIN_EXE_fastprep"))
        .args([
            "-1",
            input.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
            "--mean-q",
            "60",
            "--no-adapter",
            "--allow-total-drop",
            "-q",
        ])
        .output()
        .expect("fastprep ran");

    assert!(
        out.status.success(),
        "--allow-total-drop must still succeed (stderr: {})",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(std::fs::read_to_string(&dest).unwrap_or_default(), "");
    let _ = std::fs::remove_dir_all(&d);
}
