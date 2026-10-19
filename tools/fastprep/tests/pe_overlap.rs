//! `--pe-overlap` (paper §`sec:adp-pe`, opt-in): paired-end overlap detection.
//!
//! When the insert is shorter than the read, R1's 3′ suffix and the
//! reverse-complement of R2's prefix describe the same fragment, so the mate
//! pair pins where the adapter readthrough begins on each mate. The detector
//! must be *conservative*: it returns nothing when the mates cannot overlap,
//! when the alignment is ambiguous, when no anchor is found, or when only a
//! chance anchor hit exists. The chance-match case copies the discipline of
//! `tools/fastcheck/src/insert.rs:651`
//! (`a_chance_anchor_hit_does_not_produce_a_measurement`).
//!
//! **Default off, byte-identical when off**: every detection step is gated
//! behind `rt.pe_overlap`, and a pair whose detection returns nothing must run
//! exactly as the flag being absent. The last test of this file is that
//! property, asserted on the full pipeline output (FASTQ and PXMT).
//!
//! Test layout: the first block is the detector's unit contract (pure
//! functions, no I/O); the second is the wiring contract end-to-end.

use fastprep::pe_overlap::{detect_pe_overlap, revcomp};
use fastprep::pipeline::{run_paired, run_single};

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

// ---------------------------------------------------------------------------
// Fixtures and helpers.
// ---------------------------------------------------------------------------

/// A deterministic ACGT generator: the same sequence for a given seed on every
/// machine, so a test that finds a "chance" alignment stays reproducible.
fn dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed | 1;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        _ => b'N',
    }
}

/// A fragment of `insert` biology bases followed by `a` adapter bases, the
/// standard Illumina layout: R1 reads `insert` 5′→3′, R2 reads the same
/// fragment from the opposite strand, so R2 is the reverse complement of the
/// fragment read 5′→3′.
///
/// Both mates are the same length and both may carry readthrough: with
/// `insert < read_len` the sequencer runs past the fragment into the adapter
/// ligated at the far end, on BOTH mates.
fn pair(read_len: usize, insert: usize, seed: u64) -> (Vec<u8>, Vec<u8>) {
    let frag = dna(insert, seed);
    // R1 reads the fragment 5'->3'; R2 reads the same fragment from the
    // opposite strand, so R2 = revcomp(fragment) + its own adapter. Both mates
    // then agree on their first `insert` bases once R2 is put in the fragment's
    // own frame (the paper's R1[i] == complement(R2[I-1-i])).
    let mut r1 = frag.clone();
    r1.extend_from_slice(&dna(read_len - insert, seed ^ 0x9E37_79B9));
    let mut r2: Vec<u8> = frag.iter().rev().map(|&b| comp(b)).collect();
    r2.extend_from_slice(&dna(read_len - insert, seed ^ 0x85EB_CA6B));
    (r1, r2)
}

// ---------------------------------------------------------------------------
// Detector unit contract.
// ---------------------------------------------------------------------------

#[test]
fn detects_the_readthrough_length_on_both_mates() {
    // insert = 60, read_len = 100 ⇒ readthrough = 40 on each mate, and the
    // reverse-complement prefix of R2 equals R1's first 60 bases exactly.
    let (r1, r2) = pair(100, 60, 5);
    let got = detect_pe_overlap(&r1, &r2, 12).expect("a 60-base insert on 100-base reads must overlap");
    assert_eq!(got, (1, 40), "R1 carries 40 bases of adapter readthrough");
    let got = detect_pe_overlap(&r2, &r1, 12).expect("symmetric in the mate order");
    assert_eq!(got, (1, 40), "R2 carries the same 40 bases of readthrough");
}

#[test]
fn mates_too_long_to_overlap_return_none() {
    // insert = 100 on 100-base reads: the mates touch but do not overlap, so
    // there is no shared suffix/prefix and no boundary to pin.
    let (r1, r2) = pair(100, 100, 7);
    assert_eq!(detect_pe_overlap(&r1, &r2, 12), None);
    // And a genuinely longer insert: unrelated random tails, nothing shared.
    let a = dna(60, 11);
    let b = dna(60, 22);
    assert_eq!(detect_pe_overlap(&a, &b, 12), None);
}

#[test]
fn ambiguous_alignment_returns_none() {
    // Two admissible alignments at different offsets: the pair is a repeated
    // fragment (a 3-base unit), so the anchor matches at several offsets and
    // each offset verifies over the full shared window. The detector must
    // refuse rather than pick one.
    let unit = b"ACG";
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for i in 0..25 {
        r1.push(unit[i % 3]);
        r2.push(unit[i % 3]);
    }
    assert_eq!(detect_pe_overlap(&r1, &r2, 12), None,
        "a repeated fragment admits several alignments: refuse, do not guess");
}

#[test]
fn no_anchor_returns_none() {
    // The mates share a long region, but neither of the 12-base anchors at the
    // two ends of the shared region matches the other read — a chance anchor
    // hit is not required to be absent for this to happen, so make it certain:
    // R1's 3′ 12-mer occurs nowhere in R2's frame.
    let (mut r1, mut r2) = pair(80, 40, 13);
    // Overwrite the anchors with a fixed, non-matching sequence.
    let n = r1.len();
    r1[n - 12..].copy_from_slice(b"ACGTACGTACGT");
    r2[..12].copy_from_slice(b"TTTTTTTTTTTT");
    assert_eq!(detect_pe_overlap(&r1, &r2, 12), None,
        "no anchor match anywhere in the other read: nothing to verify");
}

#[test]
fn a_chance_anchor_hit_is_never_accepted() {
    // The discipline of fastcheck's `a_chance_anchor_hit_does_not_produce_a_measurement`:
    // the anchor matches a window of R1 by construction, but the rest of the
    // alignment is not the complement of R1 — the full verification must reject
    // it, so the pair returns nothing rather than a fabricated boundary.
    let s1 = dna(60, 17);
    let p = 40;
    // The anchor (revcomp(R2)'s first biology bases) matches s1 at p, and the
    // remainder of R2's frame is built from s1's own bases — so no alignment
    // between R1's suffix and R2's reverse complement can verify.
    let mut frame = s1[p..p + 12].to_vec();
    for &b in &s1[p + 12..] {
        frame.push(match b {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            b'T' => b'G',
            _ => b'N',
        });
    }
    frame.resize(60, b'N');
    let r2 = revcomp(&frame);
    assert_eq!(detect_pe_overlap(&s1, &r2, 12), None,
        "a chance anchor hit must not produce a measurement");
}

#[test]
fn revcomp_is_an_involution() {
    let s = dna(37, 29);
    let once = revcomp(&s);
    assert_eq!(revcomp(&once), s);
}

#[test]
fn an_exact_boundary_is_not_readthrough() {
    // insert == read_len: the mates share the whole read, so the boundary is
    // at 0 and there is no adapter to trim.
    let (r1, r2) = pair(100, 100, 31);
    let shared = r1[..100].iter().zip(r2[..100].iter().rev()).all(|(&a, &b)| {
        comp(a) == b
    });
    assert!(shared, "the fixture must actually share the whole read");
    assert_eq!(detect_pe_overlap(&r1, &r2, 12), None,
        "a boundary of 0 is no boundary at all");
}

#[test]
fn a_single_error_inside_the_overlap_still_detects() {
    // One substitution error in the shared region: the pair is still the same
    // fragment, so the boundary must be found (the paper's ε > 0 exists for
    // this), and it must land on the same readthrough length.
    let (mut r1, r2) = pair(100, 60, 33);
    let i = 30;
    r1[i] = if r1[i] == b'A' { b'C' } else { b'A' };
    assert_eq!(detect_pe_overlap(&r1, &r2, 12), Some((1, 40)));
}

// ---------------------------------------------------------------------------
// Pipeline plumbing (shared with the end-to-end tests below).
// ---------------------------------------------------------------------------

fn tmp_path(name: &str) -> PathBuf {
    let id = {
        static C: AtomicU64 = AtomicU64::new(0);
        C.fetch_add(1, Ordering::Relaxed)
    };
    let mut p = std::env::temp_dir();
    p.push(format!(
        "pe-overlap-it-{}-{}-{}",
        std::process::id(),
        id,
        name
    ));
    p
}

fn write_plain(path: &Path, body: &str) {
    let mut f = File::create(path).unwrap();
    f.write_all(body.as_bytes()).unwrap();
}

fn read_plain(path: &Path) -> String {
    let mut s = String::new();
    File::open(path).unwrap().read_to_string(&mut s).unwrap();
    s
}

/// One pair as FASTQ text, with quality `I` (Phred 40) throughout.
fn pair_fastq(p1: &[u8], p2: &[u8], name: &str) -> (String, String) {
    let rec = |id: &str, s: &[u8]| {
        format!(
            "@{id}\n{}\n+\n{}\n",
            String::from_utf8_lossy(s),
            "I".repeat(s.len())
        )
    };
    (rec(name, p1), rec(name, p2))
}

fn runtime(pe_overlap: bool) -> fastprep::cli::Runtime {
    fastprep::cli::Runtime {
        threads: 1,
        chunk: 8,
        compression: fastprep::io::CompressionLevels::default(),
        phred_offset: 33,
        phred_auto: false,
        metadata: None,
        json: None,
        dry_run: false,
        quiet: false,
        adapter_source: fastprep::cli::AdapterSource::None,
        adapter_error_rate: 0.2,
        min_adapter_overlap_se: 14,
        pe_overlap,
        dedup: None,
        thread_pool: None,
    }
}

fn params() -> fastprep::algorithm::Params {
    // The shipped winner's own config (CLI-SPEC): adapter matching OFF so the
    // only possible adapter trim is the PE one, mean-q OFF so the trim cannot
    // be masked by a drop, and Δ_L = L (min-len 0) so nothing is dropped for
    // length. What is left is the boundary the detector must find.
    fastprep::algorithm::Params {
        min_len: fastprep::algorithm::MinLen::Abs(0),
        mean_q: None,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// THE BYTE-IDENTICAL TEST — the critical no-regression requirement.
// ---------------------------------------------------------------------------

/// With `--pe-overlap` off, and on any pair for which detection returns
/// nothing, the pipeline output must be byte-identical to the flag-off run.
///
/// This is not "close": the FASTQ and the PXMT records are compared byte for
/// byte, and the set of pairs here is deliberately the set the mode CANNOT
/// fire on (mates too long to overlap, plus unrelated random mates).
#[test]
fn pairs_without_overlap_are_byte_identical_with_and_without_the_flag() {
    // Long-insert pairs (insert == read_len ⇒ no overlap at any I < L) and
    // unrelated mates (nothing shared). Neither can produce a boundary.
    let mut r1_body = String::new();
    let mut r2_body = String::new();
    for i in 0..8 {
        let (a, b) = pair(100, 100, 1000 + i as u64);
        let (qa, qb) = pair_fastq(&a, &b, &format!("long{i}"));
        r1_body.push_str(&qa);
        r2_body.push_str(&qb);
    }
    for i in 0..8 {
        let a = dna(60, 2000 + i as u64);
        let b = dna(60, 3000 + i as u64);
        let (qa, qb) = pair_fastq(&a, &b, &format!("unrelated{i}"));
        r1_body.push_str(&qa);
        r2_body.push_str(&qb);
    }

    let in1 = tmp_path("in_r1.fq");
    let in2 = tmp_path("in_r2.fq");
    write_plain(&in1, &r1_body);
    write_plain(&in2, &r2_body);

    let p = params();
    let mut runs: Vec<(PathBuf, PathBuf, PathBuf)> = Vec::new();
    for (tag, on) in [("off", false), ("on", true)] {
        let out1 = tmp_path(&format!("{tag}_1.fq"));
        let out2 = tmp_path(&format!("{tag}_2.fq"));
        let meta = tmp_path(&format!("{tag}.pxmt"));
        let mut rt = runtime(on);
        rt.metadata = Some(meta.clone());
        let st = run_paired(
            &in1, &in2, &out1, &out2, None, None, &p, &p, &rt,
        )
        .expect("paired run must succeed");
        assert_eq!(st.pairs_out, 16, "{tag}: every pair is kept (Δ_L = L, mean-q off)");
        runs.push((out1, out2, meta));
    }

    let (o1_off, o2_off, m_off) = &runs[0];
    let (o1_on, o2_on, m_on) = &runs[1];
    assert_eq!(
        read_plain(o1_off),
        read_plain(o1_on),
        "R1 output must be byte-identical when no pair can overlap"
    );
    assert_eq!(
        read_plain(o2_off),
        read_plain(o2_on),
        "R2 output must be byte-identical when no pair can overlap"
    );
    assert_eq!(
        std::fs::read(m_off).unwrap(),
        std::fs::read(m_on).unwrap(),
        "PXMT records must be byte-identical when no pair can overlap"
    );

    for p in [&in1, &in2, o1_off, o2_off, o1_on, o2_on, m_off, m_on] {
        std::fs::remove_file(p).ok();
    }
}

/// The same property on single-end input: the flag has no paired path there,
/// so single-end output must be unchanged as well.
#[test]
fn single_end_output_is_untouched_by_the_flag() {
    let a = dna(60, 41);
    let (qa, _) = pair_fastq(&a, &a, "s1");
    let in1 = tmp_path("se.fq");
    write_plain(&in1, &qa);

    let p = params();
    let mut outs = Vec::new();
    for (tag, on) in [("off", false), ("on", true)] {
        let out = tmp_path(&format!("se_{tag}.fq"));
        run_single(&in1, &out, &p, &runtime(on)).expect("single-end run must succeed");
        outs.push((out, tag));
    }
    assert_eq!(
        read_plain(&outs[0].0),
        read_plain(&outs[1].0),
        "single-end output must be identical: the flag is a paired-end mode"
    );

    for (p, _) in outs {
        std::fs::remove_file(p).ok();
    }
    std::fs::remove_file(&in1).ok();
}

// ---------------------------------------------------------------------------
// THE DETECTOR-REACHES-THE-TRIM TEST.
// ---------------------------------------------------------------------------

/// With the flag on, a pair whose insert is shorter than the read must be
/// trimmed at the boundary the overlap reveals — the readthrough length on
/// both mates — while the same pair with the flag off is left alone (no
/// adapter set is configured, so there is nothing else to trim).
#[test]
fn a_pair_with_overlap_trims_differently_with_the_flag_on() {
    let (r1, r2) = pair(100, 60, 5);
    let (qa, qb) = pair_fastq(&r1, &r2, "short_insert");
    let in1 = tmp_path("ov_r1.fq");
    let in2 = tmp_path("ov_r2.fq");
    write_plain(&in1, &qa);
    write_plain(&in2, &qb);

    let p = params();
    let mut runs: Vec<(String, String)> = Vec::new();
    for (tag, on) in [("off", false), ("on", true)] {
        let out1 = tmp_path(&format!("ov_{tag}_1.fq"));
        let out2 = tmp_path(&format!("ov_{tag}_2.fq"));
        let st = run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &runtime(on))
            .expect("paired run must succeed");
        assert_eq!(st.pairs_out, 1, "{tag}: the pair survives (Δ_L = L)");
        runs.push((read_plain(&out1), read_plain(&out2)));
    }

    let (off1, off2) = &runs[0];
    let (on1, on2) = &runs[1];

    // Off: nothing trimmed (60 insert + 40 adapter, all of it written).
    assert_eq!(off1, &qa, "without the flag the read is untouched");
    assert_eq!(off2, &qb, "without the flag the mate is untouched");

    // On: exactly the readthrough is removed, on both mates.
    let kept1 = on1.lines().nth(1).unwrap();
    let kept2 = on2.lines().nth(1).unwrap();
    assert_eq!(kept1.len(), 60, "R1 must be cut to the 60-base insert, got {}", kept1.len());
    assert_eq!(kept2.len(), 60, "R2 must be cut to the 60-base insert, got {}", kept2.len());
    assert_ne!(on1, off1, "the flag must change R1");
    assert_ne!(on2, off2, "the flag must change R2");
    // And the kept bases are the fragment, not a truncated tail of it.
    assert_eq!(kept1, String::from_utf8_lossy(&r1[..60]));
    assert_eq!(kept2, String::from_utf8_lossy(&r2[..60]));

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
    for tag in ["off", "on"] {
        for m in ["1", "2"] {
            std::fs::remove_file(tmp_path(&format!("ov_{tag}_{m}.fq"))).ok();
        }
    }
}
