//! Regression tests for the `fastprep-meta` PXMT viewer binary.
//!
//! Three defects are pinned here:
//!
//! 1. The viewer rejected the very files `fastprep --metadata` writes: the
//!    writer stamps `metadata::VERSION` (3) into the header, while the reader
//!    hard-coded `version[0] != 2`.  Since both constants date from the
//!    initial commit the viewer has never read a single file.
//! 2. `fastprep-meta --help` panicked in a debug build because the clap
//!    conflict lists named arguments by their kebab-case spellings
//!    (`"export-fasta"`, `"out-folder"`) instead of the Rust field names
//!    (`export_fasta`, `out_folder`).  clap's debug assertions abort on the
//!    unknown id, and in a release build the conflict is silently inert so an
//!    invalid flag combination is accepted.
//!
//! Both are exercised through the built binary, which is the only surface the
//! bug lives on (`src/bin/pxmt.rs` is not part of the library).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn fastprep_meta() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep-meta")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_meta_viewer_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(bin: &str, dir: &Path, args: &[&str]) -> Output {
    Command::new(bin)
        .current_dir(dir)
        .args(args)
        .output()
        .expect("spawn binary")
}

/// Write a FASTQ large enough that its metadata exceeds one BGZF member.
///
/// BGZF blocks are capped at 64 KiB, so a file under that fits in a SINGLE
/// gzip member — and `GzDecoder` reads one member fine. A regression test for
/// the multi-member read must therefore exceed 64 KiB of metadata, or it
/// passes against the broken decoder too. At ~40 bytes per record that is a
/// few thousand reads.
fn write_bulk_fastq(dir: &Path, n: usize) -> PathBuf {
    let p = dir.join("bulk.fq");
    let mut body = String::with_capacity(n * 80);
    for i in 0..n {
        // A 5' polyA run on every read, so every read is *touched* and gets a
        // metadata record (the writer skips untouched reads by design).
        let seq = format!("{}CGTACGTACGTACGTACGTACGT", "A".repeat(30));
        body.push_str(&format!(
            "@bulk{i}\n{seq}\n+\n{}\n",
            "I".repeat(seq.len()),
        ));
    }
    std::fs::write(&p, body).unwrap();
    p
}

/// Write a small FASTQ with a read that is guaranteed to be *touched*, so the
/// metadata file carries at least one record (the writer skips untouched
/// reads by design).
fn write_fastq(dir: &Path) -> PathBuf {
    let p = dir.join("in.fq");
    // r1: a long 5' polyA run that polyX trimming must cut → t5_poly > 0.
    // r2: a low-quality read, so the file has a second record type too.
    let r1_seq = format!("{}CGTACGTACGTACGTACGTACGT", "A".repeat(30));
    let r2_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string();
    let body = format!(
        "@r1\n{r1_seq}\n+\n{}\n@r2\n{r2_seq}\n+\n{}\n",
        "I".repeat(r1_seq.len()),
        "!".repeat(r2_seq.len()),
    );
    std::fs::write(&p, body).unwrap();
    p
}

/// The full `fastprep --metadata` → `fastprep-meta` round trip.
#[test]
fn metadata_written_by_fastprep_is_readable_by_the_viewer() {
    let d = tmpdir("roundtrip");
    let fq = write_fastq(&d);
    let out = d.join("out.fq");
    let meta = d.join("m.pxmt");

    let w = run(
        fastprep(),
        &d,
        &[
            "-1",
            fq.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--metadata",
            meta.to_str().unwrap(),
            "--min-len",
            "0",
            "--no-adapter",
            "-q",
        ],
    );
    assert!(
        w.status.success(),
        "fastprep --metadata must succeed, stderr: {}",
        String::from_utf8_lossy(&w.stderr)
    );

    // The header must carry exactly what the writer claims it writes.
    let bytes = std::fs::read(&meta).expect("metadata file written");
    assert_eq!(&bytes[..4], fastprep::metadata::MAGIC, "PXMT magic");
    assert_eq!(
        bytes[4],
        fastprep::metadata::VERSION,
        "header version must be fastprep's own metadata::VERSION"
    );
    assert!(
        bytes.len() > 5,
        "the file must contain at least one record"
    );

    // The viewer must read that file back: exit 0 and print records.
    let v = run(fastprep_meta(), &d, &[meta.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&v.stdout);
    let stderr = String::from_utf8_lossy(&v.stderr);
    assert!(
        v.status.success(),
        "the viewer must accept fastprep's own metadata file; exit {:?}, stderr: {stderr}",
        v.status.code()
    );
    assert!(
        !stderr.contains("unsupported PXMT version"),
        "the viewer must not reject the written version; stderr: {stderr}"
    );
    assert!(
        stdout.contains("r1"),
        "the viewer must print the recorded read IDs; stdout: {stdout}"
    );
    assert!(
        !stdout.contains("(no matching records)"),
        "records exist in the file, so the viewer must print them; stdout: {stdout}"
    );

    let _ = std::fs::remove_dir_all(&d);
}

/// A debug build must not abort on its own clap definition.
/// A `.pxmt.gz` must be GENUINELY gzip, and the viewer must read all of it.
///
/// Two defects met here. First, `open_output` inferred the codec from the
/// ".part" temp path, whose extension is not ".gz", so a file named
/// `.pxmt.gz` was written UNCOMPRESSED — the name lied. Second, once it did
/// compress, the output is BGZF (a chain of ~64 KiB gzip members) and the
/// viewer used `GzDecoder`, which stops after the FIRST member and reports
/// "failed to fill whole buffer" past 64 KiB. Each bug hid the other.
#[test]
fn a_gzipped_metadata_file_is_really_gzip_and_reads_back_whole() {
    let d = tmpdir("gz");
    // 6000 reads -> ~240 KiB of metadata -> several 64 KiB BGZF members.
    let fq = write_bulk_fastq(&d, 6000);
    let out = d.join("out.fq");
    let meta = d.join("m.pxmt.gz");

    let w = run(
        fastprep(),
        &d,
        &[
            "-1",
            fq.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--metadata",
            meta.to_str().unwrap(),
            "--min-len",
            "0",
            "--no-adapter",
            "-q",
        ],
    );
    assert!(
        w.status.success(),
        "fastprep --metadata must succeed, stderr: {}",
        String::from_utf8_lossy(&w.stderr)
    );

    // 1. The name must not lie: a .gz file is a gzip stream.
    let bytes = std::fs::read(&meta).expect("metadata file written");
    assert_eq!(
        &bytes[..2],
        &[0x1f, 0x8b],
        "a file named .pxmt.gz must actually be gzip, not raw PXMT"
    );

    // 2. The viewer must read the whole stream, not just the first member.
    //
    // The oracle is the RECORD COUNT, not the file size: BGZF members are
    // bounded by 64 KiB of UNCOMPRESSED data, and this fixture's records are
    // so repetitive that 6000 of them compress to ~14 KiB. A count of "6.0K"
    // proves every member was read; stopping after the first member would
    // report roughly 1.9K. That also keeps the test non-vacuous — an earlier
    // version used a fixture under 64 KiB and passed against the broken
    // decoder, which is exactly the trap this guard exists to catch.
    let v = run(fastprep_meta(), &d, &[meta.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&v.stdout);
    let stderr = String::from_utf8_lossy(&v.stderr);
    assert!(
        v.status.success(),
        "the viewer must read a BGZF metadata file; exit {:?}, stderr: {stderr}",
        v.status.code()
    );
    assert!(
        !stderr.contains("failed to fill whole buffer"),
        "GzDecoder stops after the first BGZF member; the viewer needs \
         MultiGzDecoder. stderr: {stderr}"
    );
    assert!(
        stdout.contains("bulk"),
        "the viewer must print the recorded read IDs; stdout: {stdout}"
    );
    // The LAST record id is the unambiguous oracle: it lives in the final
    // BGZF member, so a decoder that stops after the first member cannot
    // print it, however the viewer chooses to format its output.
    assert!(
        stdout.contains("bulk5999"),
        "the viewer must reach the LAST record, which sits in the final BGZF \
         member; its absence means the read stopped at the first member"
    );
}

#[test]
fn help_does_not_panic() {
    let d = tmpdir("help");
    let h = run(fastprep_meta(), &d, &["--help"]);
    let stdout = String::from_utf8_lossy(&h.stdout);
    let stderr = String::from_utf8_lossy(&h.stderr);
    assert!(
        h.status.success(),
        "--help must exit 0, got {:?}; stderr: {stderr}",
        h.status.code()
    );
    assert!(
        !stderr.contains("panicked"),
        "--help must not panic on clap's debug assertions; stderr: {stderr}"
    );
    assert!(
        stdout.contains("--view"),
        "--help must describe the real flags; stdout: {stdout}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// `--view` conflicts with the output modes; the conflict list must name real
/// fields, so clap rejects the combination instead of ignoring it.
#[test]
fn conflicting_output_modes_are_rejected() {
    let d = tmpdir("conflict");
    let meta = d.join("m.pxmt");
    std::fs::write(&meta, b"PXMT\x03").unwrap();
    let path = meta.to_str().unwrap();

    let cases: [&[&str]; 3] = [
        &["--view", "--summary"],
        &["--view", "--json"],
        &["--export-fasta", "x.fa", "--out-folder", "outdir"],
    ];
    for extra in cases {
        let mut args = vec![path];
        args.extend_from_slice(extra);

        let o = run(fastprep_meta(), &d, &args);
        assert_eq!(
            o.status.code(),
            Some(2),
            "`{}` must be a clap usage error (exit 2), got {:?}; stderr: {}",
            extra.join(" "),
            o.status.code(),
            String::from_utf8_lossy(&o.stderr)
        );
        let stderr = String::from_utf8_lossy(&o.stderr);
        assert!(
            stderr.contains("cannot be used with"),
            "the conflict must be reported, got: {stderr}"
        );
    }

    let _ = std::fs::remove_dir_all(&d);
}
