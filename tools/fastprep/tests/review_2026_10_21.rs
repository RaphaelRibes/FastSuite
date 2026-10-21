//! Tests pinning the fastprep defects reproduced in the 2026-10-21 review.
//!
//! Every fixture goes through the SHIPPED binary, not `Params::default()`:
//! each defect is a difference between the documented defaults and what the
//! CLI actually resolves.
//!
//! * adapter seed prefilter / mismatch tolerance (items 1 and 2)
//! * dedup anchor degenerating to the full read (item 3)
//! * silent Phred+64 read as +33 (item 4)
//! * `--min-len` overridden by the adapter budget (item 5)
//! * metadata viewer not knowing drop cause 8 = length (item 6)
//! * total-drop refusal false-firing on singleton salvage (item 7)

use std::path::{Path, PathBuf};
use std::process::Command;

const TRUSEQ: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

fn fastprep() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep")
}

fn fastprep_meta() -> &'static str {
    env!("CARGO_BIN_EXE_fastprep-meta")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "fastprep_review_20261021_{}_{}",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
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

fn records(path: &Path) -> Vec<(String, String)> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut lines = text.lines();
    while let Some(h) = lines.next() {
        let seq = lines.next().unwrap_or_default().to_string();
        if h.starts_with('@') {
            lines.next();
            lines.next();
        }
        out.push((h.to_string(), seq));
    }
    out
}

fn write_fastq(path: &Path, recs: &[(String, Vec<u8>, u8)]) {
    let mut body = String::new();
    for (name, seq, q) in recs {
        body.push_str(&format!(
            "@{name}\n{}\n+\n{}\n",
            String::from_utf8_lossy(seq),
            (*q as char).to_string().repeat(seq.len())
        ));
    }
    std::fs::write(path, body).unwrap();
}

/// A homopolymer-free, adapter-hostile body: an alternating AC repeat has no
/// run the polyX scan can see, and every 14-base window of the TruSeq adapter
/// contains at least five G/T bases, so an AC-only read cannot match it within
/// the documented tolerance by chance.
fn clean_body(len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
        .collect()
}

/// `prefix` with `positions` flipped to a base that differs from the adapter.
fn mutated(prefix_len: usize, positions: &[usize]) -> Vec<u8> {
    let mut v = TRUSEQ[..prefix_len].to_vec();
    for &p in positions {
        v[p] = if v[p] == b'C' { b'A' } else { b'C' };
    }
    v
}

fn adapter_fixture(dir: &Path, prefix_len: usize, positions: &[usize]) -> PathBuf {
    let mut seq = clean_body(100);
    seq.extend_from_slice(&mutated(prefix_len, positions));
    let fq = dir.join("insert.fastq");
    write_fastq(&fq, &[("r1".into(), seq, b'I')]);
    fq
}

fn trim_len(dir: &Path, fq: &Path, extra: &[&str]) -> (usize, String, i32) {
    let out = dir.join("out.fastq");
    let mut args = vec![
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--adapter-preset",
        "illumina",
        "--no-dedup",
        "-q",
    ];
    args.extend_from_slice(extra);
    let (err, code) = run(&args);
    let kept = records(&out);
    (kept.first().map(|(_, s)| s.len()).unwrap_or(0), err, code)
}

// ---------------------------------------------------------------------------
// Items 1 and 2 — adapter seed prefilter and mismatch tolerance.
// ---------------------------------------------------------------------------

/// Four mismatches in a 21-base TruSeq prefix at the shipped 20% rate
/// (floor(21 * 0.2) = 4) is a real match: trim 21, emit 100.
#[test]
fn a_21bp_prefix_with_4_mismatches_is_trimmed_at_the_documented_rate() {
    let d = tmpdir("adapter_4mm");
    let fq = adapter_fixture(&d, 21, &[1, 5, 9, 13]);
    let (len, err, code) = trim_len(&d, &fq, &[]);
    assert_eq!(code, 0, "run must succeed: {err}");
    assert_eq!(len, 100, "the adapter-trimmed insert must be emitted: {err}");
    std::fs::remove_dir_all(&d).ok();
}

/// Five mismatches in 21 bases is 23.8% — above the shipped 20% rate, so it is
/// not a match and must be left alone. The pre-fix `ceil` tolerance admitted 5
/// and silently trimmed 21 bases of real insert.
#[test]
fn a_21bp_prefix_with_5_mismatches_is_not_a_match_at_the_shipped_rate() {
    let d = tmpdir("adapter_5mm_reject");
    let fq = adapter_fixture(&d, 21, &[1, 5, 9, 13, 17]);
    let (len, err, code) = trim_len(&d, &fq, &[]);
    assert_eq!(code, 0, "run must succeed: {err}");
    assert_eq!(
        len, 121,
        "5/21 = 23.8% exceeds --adapter-error-rate 0.2; the read must not be trimmed: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}

/// The same five-mismatch read IS a match once the rate admits it. At 0.26 the
/// seed bound cannot cover a 21-base overlap, so the direct scan has to find
/// it; with the pre-fix non-monotone bound the position was declared
/// seed-covered and never examined.
#[test]
fn a_21bp_prefix_with_5_mismatches_is_trimmed_when_the_rate_admits_it() {
    let d = tmpdir("adapter_5mm_accept");
    let fq = adapter_fixture(&d, 21, &[1, 5, 9, 13, 17]);
    let (len, err, code) = trim_len(&d, &fq, &["--adapter-error-rate", "0.26"]);
    assert_eq!(code, 0, "run must succeed: {err}");
    assert_eq!(
        len, 100,
        "floor(21 * 0.26) = 5, so this is a real match and must be trimmed: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}

/// The `ceil` tolerance rounded 3.2 up to 4 mismatches in 16 bases — 25% at a
/// documented 20% rate. That read is not a match under `floor(16 * 0.2) = 3`
/// and must be left alone.
#[test]
fn a_16bp_prefix_with_4_mismatches_is_not_a_match_at_the_shipped_rate() {
    let d = tmpdir("adapter_k16");
    let fq = adapter_fixture(&d, 16, &[1, 5, 9, 13]);
    let (len, err, code) = trim_len(&d, &fq, &[]);
    assert_eq!(code, 0, "run must succeed: {err}");
    assert_eq!(
        len, 116,
        "4/16 = 25% exceeds --adapter-error-rate 0.2; the read must not be trimmed: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}

/// The rounded-up tolerance also fired at the shipped 14-base overlap:
/// `ceil(14 * 0.2) = 3` admits 21.4%. Three mismatches in 14 bases is not a
/// match.
#[test]
fn a_14bp_prefix_with_3_mismatches_is_not_a_match_at_the_shipped_rate() {
    let d = tmpdir("adapter_k14");
    let fq = adapter_fixture(&d, 14, &[1, 5, 9]);
    let (len, err, code) = trim_len(&d, &fq, &[]);
    assert_eq!(code, 0, "run must succeed: {err}");
    assert_eq!(
        len, 114,
        "3/14 = 21.4% exceeds --adapter-error-rate 0.2; the read must not be trimmed: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}

// ---------------------------------------------------------------------------
// Item 3 — dedup anchor when both ends are unsafe.
// ---------------------------------------------------------------------------

/// Three amplification copies of one 120-base fragment, with polyG tails of 8,
/// 12 and 14 bases and an adapter set configured: both ends are untrusted, so
/// the fingerprint must come from an interior window, not the whole read.
#[test]
fn divergent_polyg_copies_dedup_with_an_interior_anchor() {
    let d = tmpdir("dedup_interior");
    let insert = clean_body(120);
    let mut recs = Vec::new();
    for (i, tail) in [8usize, 12, 14].iter().enumerate() {
        let mut seq = insert.clone();
        seq.extend(std::iter::repeat_n(b'G', *tail));
        recs.push((format!("copy{i}"), seq, b'I'));
    }
    let fq = d.join("copies.fastq");
    write_fastq(&fq, &recs);
    let out = d.join("out.fastq");
    let json = d.join("r.json");

    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--adapter-preset",
        "illumina",
        "--json",
        json.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");
    let kept = records(&out);
    assert_eq!(
        kept.len(),
        1,
        "three copies of one fragment must dedup to one, kept {}: {err}",
        kept.len()
    );
    assert!(
        err.contains("interior"),
        "the summary line must record the degenerate anchor: {err}"
    );
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(
        doc["dedup"]["degenerate"], true,
        "the report must record the degeneracy: {doc}"
    );
    assert_eq!(doc["dedup"]["anchor_end"], "interior");
    std::fs::remove_dir_all(&d).ok();
}

/// The control: a tail-free fragment still dedups, so the fix cannot be
/// "stop deduplicating".
#[test]
fn tail_free_copies_still_dedup() {
    let d = tmpdir("dedup_clean");
    let seq = clean_body(120);
    let recs: Vec<(String, Vec<u8>, u8)> = (0..3)
        .map(|i| (format!("copy{i}"), seq.clone(), b'I'))
        .collect();
    let fq = d.join("copies.fastq");
    write_fastq(&fq, &recs);
    let out = d.join("out.fastq");
    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--adapter-preset",
        "illumina",
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");
    assert_eq!(records(&out).len(), 1, "clean copies must dedup: {err}");
    std::fs::remove_dir_all(&d).ok();
}

// ---------------------------------------------------------------------------
// Item 4 — the assumed Phred offset must be in the report and the summary.
// ---------------------------------------------------------------------------

/// Legacy Phred+64 low-quality data: byte 'B' is Q2 under +64 and Q33 under
/// +33. At the shipped defaults the census cannot exclude +64, so +33 is
/// assumed — and that assumption must be stated, not silent.
#[test]
fn an_assumed_phred_offset_is_stated_in_the_json_and_the_summary() {
    let d = tmpdir("phred_assumed");
    let recs: Vec<(String, Vec<u8>, u8)> = (0..30)
        .map(|i| (format!("r{i}"), clean_body(100), b'B'))
        .collect();
    let fq = d.join("in.fastq");
    write_fastq(&fq, &recs);
    let out = d.join("out.fastq");
    let json = d.join("r.json");

    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "--no-dedup",
        "--json",
        json.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "a +33 reading keeps these reads: {err}");
    assert_eq!(records(&out).len(), 30, "reads must be kept: {err}");
    assert!(
        err.contains("33") && err.contains("assum"),
        "the summary must state the assumed offset: {err}"
    );
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(
        doc["metrics"]["phred_offset"], 33,
        "the report must carry the offset it used: {doc}"
    );
    let a = &doc["assumptions"][0];
    assert_eq!(a["what"], "phred_offset");
    assert_eq!(a["value"], 33);
    assert!(
        a["basis"].as_str().is_some_and(|s| s.contains("census")),
        "the basis must say what the assumption rests on: {a}"
    );
    assert!(
        a["consequence"].as_str().is_some_and(|s| s.contains("31")),
        "the consequence must quantify the error: {a}"
    );
    std::fs::remove_dir_all(&d).ok();
}

/// The same data read as +64 is Q2 everywhere and the mean-quality gate drops
/// it — which is the fix the assumption's consequence points at.
#[test]
fn an_explicit_phred64_offset_drops_the_legacy_data() {
    let d = tmpdir("phred64");
    let recs: Vec<(String, Vec<u8>, u8)> = (0..30)
        .map(|i| (format!("r{i}"), clean_body(100), b'B'))
        .collect();
    let fq = d.join("in.fastq");
    write_fastq(&fq, &recs);
    let out = d.join("out.fastq");
    let json = d.join("r.json");

    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "--phred-offset",
        "64",
        "--allow-total-drop",
        "--json",
        json.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "--allow-total-drop keeps the exit code at 0: {err}");
    assert_eq!(records(&out).len(), 0, "Q2 reads fail --mean-q 25: {err}");
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(doc["metrics"]["phred_offset"], 64);
    assert_eq!(doc["assumptions"][0]["basis"], "explicit --phred-offset 64");
    std::fs::remove_dir_all(&d).ok();
}

// ---------------------------------------------------------------------------
// Item 5 — the min-len floor must not be overridden by the adapter budget.
// ---------------------------------------------------------------------------

/// A 150-base read carrying a 120-base insert plus 30 bases of TruSeq. The
/// adapter trim leaves 120 bases — the whole insert — so the read must be
/// emitted as 120 bp at the shipped `--min-len 0.9L`. The floor is a floor on
/// the insert, and the adapter budget must not veto it.
#[test]
fn a_30bp_adapter_readthrough_does_not_veto_the_min_len_floor() {
    let d = tmpdir("minlen_adapter");
    let mut seq = clean_body(120);
    seq.extend_from_slice(&TRUSEQ[..30]);
    assert_eq!(seq.len(), 150);
    let fq = d.join("in.fastq");
    write_fastq(&fq, &[("r1".into(), seq, b'I')]);
    let out = d.join("out.fastq");

    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--adapter-preset",
        "illumina",
        "--no-dedup",
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");
    let kept = records(&out);
    assert_eq!(
        kept.len(),
        1,
        "the read must be emitted, not dropped as adapter overtrim: {err}"
    );
    assert_eq!(
        kept[0].1.len(),
        120,
        "the adapter-trimmed insert must be emitted whole: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}

// ---------------------------------------------------------------------------
// Item 6 — the metadata viewer must know drop cause 8.
// ---------------------------------------------------------------------------

/// `--min-len 100` drops 40-base reads for length; the viewer must label that
/// cause, filter on it, and count it, instead of reporting "Dropped N" with
/// every cause row at zero.
#[test]
fn the_metadata_viewer_knows_the_length_drop_cause() {
    let d = tmpdir("meta_length");
    let recs: Vec<(String, Vec<u8>, u8)> = (0..40)
        .map(|i| (format!("r{i}"), clean_body(40), b'I'))
        .collect();
    let fq = d.join("in.fastq");
    write_fastq(&fq, &recs);
    let out = d.join("out.fastq");
    let meta = d.join("m.pxmt");

    let (err, code) = run(&[
        "-1",
        fq.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-adapter",
        "--no-dedup",
        "--min-len",
        "100",
        "--metadata",
        meta.to_str().unwrap(),
        "--allow-total-drop",
        "-q",
    ]);
    assert_eq!(code, 0, "run must succeed: {err}");

    let s = Command::new(fastprep_meta())
        .args([meta.to_str().unwrap(), "--status", "length", "--summary"])
        .output()
        .expect("spawn fastprep-meta");
    // The human summary goes to stderr; --json is the stdout stream.
    let summary = String::from_utf8_lossy(&s.stderr).to_string();
    assert!(
        s.status.success(),
        "--status length must be accepted, got: {summary}"
    );
    assert!(
        summary.contains("below --min-len"),
        "the summary must carry the cause row: {summary}"
    );
    assert!(
        !summary.contains("unknown"),
        "no record may be labelled unknown: {summary}"
    );
    // The cause row must account for all 40 drops, not zero.
    let line = summary
        .lines()
        .find(|l| l.contains("below --min-len"))
        .unwrap();
    assert!(
        line.contains("40"),
        "the length row must count the drops: {summary}"
    );

    let j = Command::new(fastprep_meta())
        .args([meta.to_str().unwrap(), "--status", "length", "--json"])
        .output()
        .expect("spawn fastprep-meta");
    let jd: serde_json::Value =
        serde_json::from_slice(&j.stdout).expect("viewer --json must parse");
    let recs = jd["records"].as_array().expect("records array");
    assert_eq!(recs.len(), 40, "--status length must select every dropped read");
    assert_eq!(recs[0]["status"], "length");
    std::fs::remove_dir_all(&d).ok();
}

// ---------------------------------------------------------------------------
// Item 7 — a total drop is not total when singletons were written.
// ---------------------------------------------------------------------------

/// 100 pairs where R2 is all-Q2: with `--singletons1` the run writes 100 R1
/// singletons, so it did not produce an empty output and must exit 0.
#[test]
fn singleton_salvage_is_not_a_total_drop() {
    let d = tmpdir("singletons");
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for i in 0..100usize {
        // Distinct bodies so dedup cannot collapse them.
        let mut a = clean_body(100);
        a[0] = if i % 2 == 0 { b'A' } else { b'C' };
        let mut b = clean_body(100);
        b[1] = if (i / 2) % 2 == 0 { b'A' } else { b'C' };
        r1.push((format!("p{i}"), a, b'I'));
        r2.push((format!("p{i}"), b, b'#'));
    }
    let f1 = d.join("r1.fastq");
    let f2 = d.join("r2.fastq");
    write_fastq(&f1, &r1);
    write_fastq(&f2, &r2);
    let o1 = d.join("o1.fastq");
    let o2 = d.join("o2.fastq");
    let sg1 = d.join("sg1.fastq");

    let (err, code) = run(&[
        "-1",
        f1.to_str().unwrap(),
        "-2",
        f2.to_str().unwrap(),
        "-o",
        o1.to_str().unwrap(),
        "-p",
        o2.to_str().unwrap(),
        "--singletons1",
        sg1.to_str().unwrap(),
        "--no-dedup",
        "-q",
    ]);
    assert_eq!(
        code, 0,
        "singletons were written, so this is not a total drop: {err}"
    );
    assert_eq!(
        records(&sg1).len(),
        100,
        "the 100 R1 singletons must be written: {err}"
    );
    std::fs::remove_dir_all(&d).ok();
}
