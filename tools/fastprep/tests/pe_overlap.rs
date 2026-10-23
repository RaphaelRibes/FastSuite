//! `--pe-overlap <never|auto|always>` (paper §`sec:adp-pe`, default `auto`):
//! paired-end overlap detection.
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
//! # The tri-state, and the guard it exists to serve
//!
//! The mode runs on `Auto` (the default): only when the adapter source is
//! EXPLICIT (`--adapter-preset` / `--adapters`). With auto-detection or
//! `--no-adapter` the mode stays off, because the safety guard needs a set to
//! check the boundary against. `Always` is the explicit opt-in to run
//! regardless, but the guard is structural: **no set ⇒ no boundary**, so even
//! `always` with no set refuses to trim biology (see
//! `no_adapter_set_means_no_boundary`).
//!
//! **Byte-identical when off**: every detection step is gated behind
//! `rt.pe_overlap`, and a pair whose detection returns nothing must run
//! exactly as the flag being absent. The byte-identical test below is that
//! property, asserted on the full pipeline output (FASTQ and PXMT).
//!
//! Test layout: the first block is the detector's unit contract (pure
//! functions, no I/O); the second is the wiring contract end-to-end; the
//! third is the mode's semantics (the tri-state) as the CLI resolves it.

use fastprep::cli::{AdapterSource, PeOverlapMode, Runtime};
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

fn runtime(pe_overlap: PeOverlapMode) -> Runtime {
    Runtime {
        threads: 1,
        chunk: 8,
        compression: fastprep::io::CompressionLevels::default(),
        phred_offset: 33,
        phred_auto: false,
        metadata: None,
        json: None,
        dry_run: false,
        quiet: false,
        allow_total_drop: false,
        adapter_source: AdapterSource::None,
        adapter_error_rate: 0.2,
        min_adapter_overlap_se: 14,
        pe_overlap,
        adapter_evidence: fastprep::adapter::AdapterEvidence::Floor,
        pe_ambiguous_smallest: false,
        pe_correct: None,
        pe_correct_min_overlap: 30,
        dedup: None,
        thread_pool: None,
    }
}

/// An explicit adapter set that shares no 6-mer with the fixtures: the ordinary
/// sequence matcher finds nothing, so it isolates the PE decision.
fn inert_params() -> fastprep::algorithm::Params {
    let set = fastprep::adapter::AdapterSet::new(
        vec![fastprep::adapter::Adapter {
            name: "inert".into(),
            seq: b"AGTCAGTCAGTCAGTCAGTCAGTCAGTC".to_vec(),
        }],
        0.2,
        14,
    );
    fastprep::algorithm::Params {
        adapter: Some(std::sync::Arc::new(set)),
        ..params()
    }
}

/// The shipped adapter set, as `finalize()` builds it for `--adapter-preset`.
/// Both mates in a real TruSeq library carry TruSeq readthrough past the
/// insert, so this is the set the PE boundary is checked against.
fn explicit_params() -> fastprep::algorithm::Params {
    let set = fastprep::adapter::AdapterSet::new(
        fastprep::adapter::preset_illumina(),
        0.2,
        14,
    );
    fastprep::algorithm::Params {
        adapter: Some(std::sync::Arc::new(set)),
        ..params()
    }
}

/// A pair whose insert is shorter than the read, with a REAL adapter readthrough
/// past the fragment on BOTH mates. The overlap pins the boundary at `insert`,
/// and the sequence past it is adapter — the mode's real target.
///
/// A real TruSeq readthrough is longer than the 33-base adapter: the read runs
/// through the adapter into the downstream primer, so the tail is the adapter
/// followed by primer sequence, not a repeat of the adapter (which would
/// self-match the matcher). That is real adapter sequence by construction,
/// accepted by `--adapters FASTA`, and what `readthrough_is_adapter(r1, …)`
/// must recognise.
///
/// `both_mates` picks the geometry: `true` gives each mate its own adapter
/// (TruSeq R1 on R1, TruSeq R2 on R2), `false` gives R1 the adapter and R2 only
/// biology — the case that proves the guard checks the MATE, not a shared
/// boundary.
fn adapter_pair(read_len: usize, insert: usize, seed: u64, both_mates: bool) -> (Vec<u8>, Vec<u8>) {
    let frag = dna(insert, seed);
    let a1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let a2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
    // Downstream primer (non-adapter, non-repeating) to lengthen the tail.
    let p1 = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";
    let p2 = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";
    let n = read_len - insert;
    let tail = |a: &[u8], p: &[u8]| -> Vec<u8> {
        let mut v = a.to_vec();
        let mut i = 0;
        while v.len() < n {
            v.push(p[i % p.len()]);
            i += 1;
        }
        v.truncate(n);
        v
    };
    let mut r1 = frag.clone();
    r1.extend_from_slice(&tail(a1, p1));
    let mut r2: Vec<u8> = frag.iter().rev().map(|&b| comp(b)).collect();
    if both_mates {
        r2.extend_from_slice(&tail(a2, p2));
    } else {
        r2.extend_from_slice(&dna(n, seed ^ 0x85EB_CA6B));
    }
    (r1, r2)
}

/// Write a FASTQ body from `(seq, name)` pairs.
fn fastq_body(pairs: &[(&[u8], &str)]) -> String {
    let mut body = String::new();
    for (s, name) in pairs {
        body.push_str(&pair_fastq(s, s, name).0);
    }
    body
}

/// Run one paired-end pass and return the R1 output text.
fn run_r1(in1: &Path, in2: &Path, p: &fastprep::algorithm::Params, rt: &Runtime) -> String {
    let out1 = tmp_path("run_1.fq");
    let out2 = tmp_path("run_2.fq");
    run_paired(in1, in2, &out1, &out2, None, None, p, p, rt).expect("paired run must succeed");
    let body = read_plain(&out1);
    std::fs::remove_file(&out1).ok();
    std::fs::remove_file(&out2).ok();
    body
}

/// The shipped winner's own config (CLI-SPEC): adapter matching OFF so the
/// only possible adapter trim is the PE one, mean-q OFF so the trim cannot
/// be masked by a drop, and Δ_L = L (min-len 0) so nothing is dropped for
/// length. What is left is the boundary the detector must find.
fn params() -> fastprep::algorithm::Params {
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
    for (tag, on) in [("off", PeOverlapMode::Never), ("on", PeOverlapMode::Always)] {
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
    for (tag, on) in [("off", PeOverlapMode::Never), ("on", PeOverlapMode::Always)] {
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

/// With the mode on and a real adapter set, a pair whose insert is shorter than
/// the read is trimmed at the boundary the overlap reveals. R1 carries adapter
/// (so the ordinary matcher trims it too); R2 carries only biology, so the ONLY
/// thing that can trim R2 is the PE boundary — which is what makes this a test
/// of the mode rather than of the matcher.
#[test]
fn a_pair_with_overlap_trims_differently_with_the_flag_on() {
    let (r1, r2) = adapter_pair(150, 100, 5, false);
    let (qa, qb) = pair_fastq(&r1, &r2, "short_insert");
    let in1 = tmp_path("ov_r1.fq");
    let in2 = tmp_path("ov_r2.fq");
    write_plain(&in1, &qa);
    write_plain(&in2, &qb);

    let p = explicit_params();
    let mut runs: Vec<(String, String)> = Vec::new();
    for (tag, on) in [("off", PeOverlapMode::Never), ("on", PeOverlapMode::Auto)] {
        let out1 = tmp_path(&format!("ov_{tag}_1.fq"));
        let out2 = tmp_path(&format!("ov_{tag}_2.fq"));
        let st = run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &runtime(on))
            .expect("paired run must succeed");
        assert_eq!(st.pairs_out, 1, "{tag}: the pair survives (Δ_L = L)");
        runs.push((read_plain(&out1), read_plain(&out2)));
    }

    let (off1, off2) = &runs[0];
    let (on1, on2) = &runs[1];

    // Off: R1 is trimmed by the ordinary matcher (adapter present); R2 is
    // biology, so nothing trims it.
    assert_eq!(off1.lines().nth(1).unwrap().len(), 100, "off: R1 trimmed by the matcher");
    assert_eq!(off2.lines().nth(1).unwrap().len(), 150, "off: R2 (biology) is untouched");

    // On: the PE boundary is applied on BOTH mates. R2 can only be trimmed by
    // the mode — its tail is not adapter — so 100 proves the boundary fired.
    let kept1 = on1.lines().nth(1).unwrap();
    let kept2 = on2.lines().nth(1).unwrap();
    assert_eq!(kept1.len(), 100, "R1 must be cut to the 100-base insert, got {}", kept1.len());
    assert_eq!(kept2.len(), 100, "R2 must be cut to the 100-base insert, got {}", kept2.len());
    assert_ne!(on2, off2, "the mode must change R2, whose tail is not adapter");
    assert_eq!(kept1, String::from_utf8_lossy(&r1[..100]));
    assert_eq!(kept2, String::from_utf8_lossy(&r2[..100]));

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
    for tag in ["off", "on"] {
        for m in ["1", "2"] {
            std::fs::remove_file(tmp_path(&format!("ov_{tag}_{m}.fq"))).ok();
        }
    }
}

// ---------------------------------------------------------------------------
// THE TRI-STATE, AND THE GUARD IT EXISTS TO SERVE.
// ---------------------------------------------------------------------------

/// The CLI-level resolution of `Auto`: the mode runs only when the adapter
/// source is EXPLICIT. This is what `finalize()` computes, asserted through the
/// public `Cli` API (not a private field), so it cannot be edited away by a
/// refactor of the flag alone.
#[test]
fn auto_resolves_on_only_for_an_explicit_adapter_source() {
    use clap::Parser;
    let dir = std::env::temp_dir().join("fastprep_pe_tristate");
    std::fs::create_dir_all(&dir).unwrap();
    let adapters = dir.join("adapters.fa");
    write_plain(&adapters, ">tru\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n");
    // Auto-detection peeks at these; they must exist and hold no catalogue
    // fingerprint, so detection finds nothing.
    let clean1 = dir.join("clean1.fq");
    let clean2 = dir.join("clean2.fq");
    let body = fastq_body(&[(&dna(60, 3), "n1")]);
    write_plain(&clean1, &body);
    write_plain(&clean2, &body);

    // Explicit preset ⇒ Auto is ON (the guard has a set to check against).
    let cli = fastprep::cli::Cli::parse_from([
        "fastprep", "-1", clean1.to_str().unwrap(), "-2", clean2.to_str().unwrap(),
        "-o", "o1.fq", "-p", "o2.fq", "--adapter-preset", "illumina",
    ]);
    assert_eq!(cli.pe_overlap, PeOverlapMode::Auto, "the default is auto");
    let (_, rt) = cli.finalize().expect("an explicit adapter source finalizes");
    assert_eq!(
        rt.pe_overlap,
        PeOverlapMode::Auto,
        "an explicit source keeps auto, and the guard then has a set"
    );
    assert!(rt.pe_overlap.runs_with(&rt.adapter_source), "auto + explicit preset must run the mode");

    // Explicit FASTA ⇒ Auto is ON.
    let cli = fastprep::cli::Cli::parse_from([
        "fastprep", "-1", clean1.to_str().unwrap(), "-2", clean2.to_str().unwrap(),
        "-o", "o1.fq", "-p", "o2.fq", "--adapters", adapters.to_str().unwrap(),
    ]);
    let (_, rt) = cli.finalize().expect("an explicit FASTA finalizes");
    assert!(rt.pe_overlap.runs_with(&rt.adapter_source), "auto + explicit FASTA must run the mode");

    // Auto-detection ⇒ the mode stays OFF (the reproduction's exact config).
    let cli = fastprep::cli::Cli::parse_from([
        "fastprep", "-1", clean1.to_str().unwrap(), "-2", clean2.to_str().unwrap(),
        "-o", "o1.fq", "-p", "o2.fq", "--force",
    ]);
    let (_, rt) = cli.finalize().expect("auto-detection finalizes");
    assert_eq!(rt.pe_overlap, PeOverlapMode::Never, "auto + auto-detect resolves to never");

    // --no-adapter ⇒ the mode stays OFF.
    let cli = fastprep::cli::Cli::parse_from([
        "fastprep", "-1", clean1.to_str().unwrap(), "-2", clean2.to_str().unwrap(),
        "-o", "o1.fq", "-p", "o2.fq", "--no-adapter",
    ]);
    let (_, rt) = cli.finalize().expect("--no-adapter finalizes");
    assert_eq!(rt.pe_overlap, PeOverlapMode::Never, "auto + --no-adapter resolves to never");

    std::fs::remove_file(&adapters).ok();
    std::fs::remove_file(&clean1).ok();
    std::fs::remove_file(&clean2).ok();
}

/// The three spellings parse, and `never`/`always` are respected regardless of
/// the adapter source.
#[test]
fn the_tri_state_parses() {
    use clap::Parser;
    let parse = |mode: &str| {
        let cli = fastprep::cli::Cli::parse_from([
            "fastprep", "-1", "in.fq", "-2", "in2.fq", "-o", "o1.fq", "-p", "o2.fq",
            "--pe-overlap", mode, "--adapter-preset", "illumina",
        ]);
        cli.pe_overlap
    };
    assert_eq!(parse("never"), PeOverlapMode::Never);
    assert_eq!(parse("auto"), PeOverlapMode::Auto);
    assert_eq!(parse("always"), PeOverlapMode::Always);
    assert!(!PeOverlapMode::Never.runs(), "never never runs");
    assert!(PeOverlapMode::Always.runs(), "always runs");
    // The bare spelling (no value) is the pre-tri-state alias for `always`, so
    // existing invocations that asked for the mode keep working.
    let bare = fastprep::cli::Cli::parse_from([
        "fastprep", "-1", "in.fq", "-2", "in2.fq", "-o", "o1.fq", "-p", "o2.fq",
        "--pe-overlap",
    ]);
    assert_eq!(bare.pe_overlap, PeOverlapMode::Always, "bare --pe-overlap == always");
    // A bad mode is a clap error, not a silent fallback.
    let bad = fastprep::cli::Cli::try_parse_from([
        "fastprep", "-1", "in.fq", "-2", "in2.fq", "-o", "o1.fq", "-p", "o2.fq",
        "--pe-overlap", "sometimes",
    ]);
    assert!(bad.is_err(), "--pe-overlap must reject an unknown mode");
}

/// The detector's guard, asserted directly: **no set ⇒ no boundary**. This is
/// the structural invariant that makes `--pe-overlap=always` safe to expose:
/// even the explicit opt-in cannot apply a boundary it has no set to check.
#[test]
fn no_adapter_set_means_no_boundary() {
    // A pair whose mates genuinely overlap (insert 100 of 150): detection
    // succeeds and pins the boundary.
    let (r1, r2) = pair(150, 100, 5);
    assert_eq!(detect_pe_overlap(&r1, &r2, 12), Some((1, 50)));
    // Without a set the boundary is refused — the guard cannot be bypassed by
    // the flag, because the detector has nothing to check the tail against.
    assert_eq!(
        fastprep::pe_overlap::detect_pe_overlap_trim(&r1, &r2, 12, None),
        None,
        "no adapter set must never yield a PE boundary"
    );
}

/// The reproduction, end-to-end: auto-detection finds nothing (a random tail
/// holds no catalogue fingerprint), so the mode is OFF at its default and the
/// read is left exactly as it came in.
#[test]
fn auto_detect_with_no_adapter_leaves_the_read_untouched_at_the_default() {
    // Random tails: no TruSeq/Nextera/SmallRNA fingerprint, so auto-detection
    // finds nothing and `Params.adapter` is None.
    let (r1, r2) = pair(150, 100, 5);
    let in1 = tmp_path("repro_r1.fq");
    let in2 = tmp_path("repro_r2.fq");
    write_plain(&in1, &fastq_body(&[(&r1, "s1")]));
    write_plain(&in2, &fastq_body(&[(&r2, "s1")]));

    let p = params();
    let body = run_r1(&in1, &in2, &p, &runtime(PeOverlapMode::Auto));
    let kept = body.lines().nth(1).unwrap();
    assert_eq!(kept.len(), 150, "the default must not trim a random tail, got {}", kept.len());
    assert_eq!(kept, String::from_utf8_lossy(&r1), "the read must be untouched");

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
}

/// The same data, explicit `--adapter-preset`, mode at its default: the guard
/// refuses the boundary because the tail past it is biology, not adapter.
#[test]
fn explicit_preset_with_a_random_tail_is_not_trimmed() {
    let (r1, r2) = pair(150, 100, 5);
    let in1 = tmp_path("preset_r1.fq");
    let in2 = tmp_path("preset_r2.fq");
    write_plain(&in1, &fastq_body(&[(&r1, "s1")]));
    write_plain(&in2, &fastq_body(&[(&r2, "s1")]));

    let p = inert_params();
    let body = run_r1(&in1, &in2, &p, &runtime(PeOverlapMode::Auto));
    let kept = body.lines().nth(1).unwrap();
    assert_eq!(kept.len(), 150, "a boundary whose tail is not adapter must not be applied, got {}", kept.len());
    assert_eq!(kept, String::from_utf8_lossy(&r1));

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
}

/// Run the built `fastprep` binary and return its R1 output text.
///
/// The binary is what a user runs, and it is the only place the whole chain
/// (clap → finalize → Runtime → pipeline → trim) is exercised together, so the
/// end-to-end cases below use it rather than the library entry points.
fn run_cli(dir: &Path, in1: &Path, in2: &Path, tag: &str, extra: &[&str]) -> String {
    let out1 = tmp_path(&format!("cli_{tag}_1.fq"));
    let out2 = tmp_path(&format!("cli_{tag}_2.fq"));
    let mut args: Vec<&str> = vec![
        "-1", in1.to_str().unwrap(), "-2", in2.to_str().unwrap(),
        "-o", out1.to_str().unwrap(), "-p", out2.to_str().unwrap(),
        "--min-len", "0", "--force",
    ];
    args.extend_from_slice(extra);
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fastprep"))
        .args(&args)
        .current_dir(dir)
        .status()
        .expect("spawn fastprep");
    assert!(status.success(), "fastprep must exit 0 for {tag}");
    let body = read_plain(&out1);
    std::fs::remove_file(&out1).ok();
    std::fs::remove_file(&out2).ok();
    body
}

/// The reproduction, through the CLI binary: auto-detection finds nothing (a
/// random tail holds no catalogue fingerprint), so the mode is OFF at its
/// default and the read is left exactly as it came in.
#[test]
fn cli_auto_detect_with_no_adapter_leaves_the_read_untouched() {
    let (r1, r2) = pair(150, 100, 5);
    let in1 = tmp_path("cli_repro_r1.fq");
    let in2 = tmp_path("cli_repro_r2.fq");
    write_plain(&in1, &fastq_body(&[(&r1, "s1")]));
    write_plain(&in2, &fastq_body(&[(&r2, "s1")]));

    // No --pe-overlap at all: the default is auto, which is off here.
    let body = run_cli(std::path::Path::new("."), &in1, &in2, "repro", &[]);
    let kept = body.lines().nth(1).unwrap();
    assert_eq!(kept.len(), 150, "the default must not trim a random tail, got {}", kept.len());

    // Explicitly `--pe-overlap auto`: same result.
    let body = run_cli(std::path::Path::new("."), &in1, &in2, "auto", &["--pe-overlap", "auto"]);
    assert_eq!(body.lines().nth(1).unwrap().len(), 150, "auto + auto-detect must not trim");

    // `--pe-overlap always`: the explicit opt-in, still refused (no set).
    let body = run_cli(std::path::Path::new("."), &in1, &in2, "always", &["--pe-overlap", "always"]);
    assert_eq!(body.lines().nth(1).unwrap().len(), 150, "always + no set must still refuse");

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
}

/// Explicit `--adapter-preset` on the same adapter-free data: the guard refuses
/// the boundary because the tail past it is biology, not adapter.
#[test]
fn cli_explicit_preset_with_a_random_tail_is_not_trimmed() {
    let (r1, r2) = pair(150, 100, 5);
    let in1 = tmp_path("cli_preset_r1.fq");
    let in2 = tmp_path("cli_preset_r2.fq");
    write_plain(&in1, &fastq_body(&[(&r1, "s1")]));
    write_plain(&in2, &fastq_body(&[(&r2, "s1")]));

    let body = run_cli(std::path::Path::new("."), &in1, &in2, "preset", &["--adapter-preset", "illumina"]);
    let kept = body.lines().nth(1).unwrap();
    assert_eq!(kept.len(), 150, "a boundary whose tail is not adapter must not be applied, got {}", kept.len());

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
}

/// The mode's real target, through the CLI: a real TruSeq readthrough past a
/// real overlap, with the explicit preset. Both mates are cut to the 100-base
/// insert — and with the mode off the same pair is untouched by the PE path.
#[test]
fn cli_explicit_preset_with_real_adapter_readthrough_trims_at_the_boundary() {
    let (r1, r2) = adapter_pair(150, 100, 5, true);
    let in1 = tmp_path("cli_adp_r1.fq");
    let in2 = tmp_path("cli_adp_r2.fq");
    write_plain(&in1, &fastq_body(&[(&r1, "s1")]));
    write_plain(&in2, &fastq_body(&[(&r2, "s1")]));

    // auto (the default) with an explicit preset: the mode runs and fires.
    let on = run_cli(std::path::Path::new("."), &in1, &in2, "adp_auto", &["--adapter-preset", "illumina"]);
    let kept1 = on.lines().nth(1).unwrap();
    assert_eq!(kept1.len(), 100, "R1 must be cut to the 100-base insert, got {}", kept1.len());
    assert_eq!(kept1, String::from_utf8_lossy(&r1[..100]));

    // never: the mode is off, so the ordinary matcher alone decides.
    let off = run_cli(std::path::Path::new("."), &in1, &in2, "adp_never",
        &["--adapter-preset", "illumina", "--pe-overlap", "never"]);
    assert_eq!(off.lines().nth(1).unwrap().len(), 100,
        "the matcher alone trims the readthrough to the same boundary");

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
}

/// `--pe-overlap=always` with auto-detection and adapter-free data: the mode
/// is deliberately opted into, but the guard still refuses — **no set, no
/// boundary**. Detection may run and report; it may not trim biology. This is
/// the conservative reading of `always`, pinned here so a future change cannot
/// silently turn it into an opt-in to trims the default forbids.
#[test]
fn always_with_no_adapter_set_still_refuses_the_boundary() {
    let (r1, r2) = pair(150, 100, 5);
    let in1 = tmp_path("always_r1.fq");
    let in2 = tmp_path("always_r2.fq");
    write_plain(&in1, &fastq_body(&[(&r1, "s1")]));
    write_plain(&in2, &fastq_body(&[(&r2, "s1")]));

    let p = params();
    let body = run_r1(&in1, &in2, &p, &runtime(PeOverlapMode::Always));
    let kept = body.lines().nth(1).unwrap();
    assert_eq!(kept.len(), 150, "always + no set must still refuse the boundary, got {}", kept.len());
    assert_eq!(kept, String::from_utf8_lossy(&r1));

    for p in [&in1, &in2] {
        std::fs::remove_file(p).ok();
    }
}


// ---------------------------------------------------------------------------
// UT-1 end to end: the evidence policy through the whole pipeline.
//
// The unit tests cover the guard in isolation. This one drives the real
// `run_paired` over a corpus whose readthrough is BELOW the overlap floor,
// and asserts the two properties that matter:
//
//   1. with `Boundary`, the pinned boundary is applied (the readthrough is
//      removed);
//   2. with `Floor` (the default), the output is BYTE-IDENTICAL to a run
//      with the policy unset -- today's behaviour is untouched.
//
// The corpus: 100-base mates with a 90-base insert, so each mate carries 10
// bases of readthrough -- below the default floor of 14. The readthrough is
// the real TruSeq/R1 prefix on R1 and the real TruSeq/R2 prefix on R2, so
// the ordinary matcher could trim it if the floor allowed.
// ---------------------------------------------------------------------------

const TRUSEQ_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const TRUSEQ_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

/// A mate pair with `readthrough` bases of adapter past the insert, the
/// adapter being the real TruSeq prefix for that mate.
fn pair_with_readthrough(
    read_len: usize,
    readthrough: usize,
    seed: u64,
) -> (Vec<u8>, Vec<u8>) {
    let insert = read_len - readthrough;
    let frag = dna(insert, seed);
    let mut r1 = frag.clone();
    r1.extend_from_slice(&TRUSEQ_R1[..readthrough]);
    let mut r2: Vec<u8> = frag.iter().rev().map(|&b| comp(b)).collect();
    r2.extend_from_slice(&TRUSEQ_R2[..readthrough]);
    (r1, r2)
}

/// The adapter set the two new tests run with: the two real TruSeq
/// prefixes, and the default floor of 14.
///
/// NOTE: the pipeline reads the set from `Params::adapter` (what the CLI
/// fills from the adapter source), not from `Runtime::adapter_source` — the
/// detector is handed `params_r1.adapter`.
fn params_with_adapters() -> fastprep::algorithm::Params {
    fastprep::algorithm::Params {
        min_len: fastprep::algorithm::MinLen::Abs(0),
        mean_q: None,
        adapter: Some(std::sync::Arc::new(fastprep::adapter::AdapterSet::new(
            vec![
                fastprep::adapter::Adapter {
                    name: "TruSeq/R1".into(),
                    seq: TRUSEQ_R1.to_vec(),
                },
                fastprep::adapter::Adapter {
                    name: "TruSeq/R2".into(),
                    seq: TRUSEQ_R2.to_vec(),
                },
            ],
            0.2,
            14,
        ))),
        ..Default::default()
    }
}

#[test]
fn the_boundary_policy_trims_a_readthrough_below_the_floor() {
    let mut r1_body = String::new();
    let mut r2_body = String::new();
    for i in 0..4u64 {
        let (a, b) = pair_with_readthrough(100, 10, 5000 + i);
        let (qa, qb) = pair_fastq(&a, &b, &format!("short{i}"));
        r1_body.push_str(&qa);
        r2_body.push_str(&qb);
    }

    let in1 = tmp_path("ut1_r1.fq");
    let in2 = tmp_path("ut1_r2.fq");
    write_plain(&in1, &r1_body);
    write_plain(&in2, &r2_body);

    let p = params_with_adapters();
    let mut outs: Vec<(PathBuf, PathBuf)> = Vec::new();
    for (tag, ev) in [
        ("floor", fastprep::adapter::AdapterEvidence::Floor),
        ("boundary", fastprep::adapter::AdapterEvidence::Boundary),
    ] {
        let out1 = tmp_path(&format!("ut1_{tag}_1.fq"));
        let out2 = tmp_path(&format!("ut1_{tag}_2.fq"));
        let mut rt = runtime(PeOverlapMode::Always);
        rt.adapter_evidence = ev;
        run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &rt)
            .expect("paired run must succeed");
        outs.push((out1, out2));
    }

    let floor_len = read_plain(&outs[0].0).lines().count();
    let boundary_len = read_plain(&outs[1].0).lines().count();
    assert_eq!(
        floor_len, boundary_len,
        "both policies keep every record; only the TRIM differs"
    );

    // The decisive difference: with `Boundary` the 10-base readthrough is
    // removed, so the emitted R1 records are 10 bases shorter.
    let floor_seq_len = read_plain(&outs[0].0)
        .lines()
        .nth(1)
        .map(|s| s.len())
        .unwrap_or(0);
    let boundary_seq_len = read_plain(&outs[1].0)
        .lines()
        .nth(1)
        .map(|s| s.len())
        .unwrap_or(0);
    assert_eq!(
        floor_seq_len - boundary_seq_len,
        10,
        "the boundary policy removes the 10-base readthrough the floor refused"
    );

    for p in [
        &in1, &in2, &outs[0].0, &outs[0].1, &outs[1].0, &outs[1].1,
    ] {
        std::fs::remove_file(p).ok();
    }
}

#[test]
fn the_floor_policy_is_byte_identical_to_no_policy_at_all() {
    // The default must not move a single byte.
    let mut r1_body = String::new();
    let mut r2_body = String::new();
    for i in 0..4u64 {
        let (a, b) = pair_with_readthrough(100, 10, 6000 + i);
        let (qa, qb) = pair_fastq(&a, &b, &format!("def{i}"));
        r1_body.push_str(&qa);
        r2_body.push_str(&qb);
    }

    let in1 = tmp_path("ut1d_r1.fq");
    let in2 = tmp_path("ut1d_r2.fq");
    write_plain(&in1, &r1_body);
    write_plain(&in2, &r2_body);

    let p = params_with_adapters();
    let mut outs: Vec<(PathBuf, PathBuf)> = Vec::new();
    for tag in ["default", "explicit_floor"] {
        let out1 = tmp_path(&format!("ut1d_{tag}_1.fq"));
        let out2 = tmp_path(&format!("ut1d_{tag}_2.fq"));
        let mut rt = runtime(PeOverlapMode::Always);
        if tag == "explicit_floor" {
            rt.adapter_evidence = fastprep::adapter::AdapterEvidence::Floor;
        }
        run_paired(&in1, &in2, &out1, &out2, None, None, &p, &p, &rt)
            .expect("paired run must succeed");
        outs.push((out1, out2));
    }

    assert_eq!(
        read_plain(&outs[0].0),
        read_plain(&outs[1].0),
        "the default policy and an explicit `floor` must produce identical R1"
    );
    assert_eq!(
        read_plain(&outs[0].1),
        read_plain(&outs[1].1),
        "the default policy and an explicit `floor` must produce identical R2"
    );

    for p in [
        &in1, &in2, &outs[0].0, &outs[0].1, &outs[1].0, &outs[1].1,
    ] {
        std::fs::remove_file(p).ok();
    }
}

// ---------------------------------------------------------------------------
// --pe-correct: posterior consensus correction over the doubly-retained mate
// overlap (src/correct.rs). At overlap disagreements the two-hypothesis
// posterior decides: when |log10 odds| >= gate the LOSING mate's base and Q
// are rewritten (the winner is untouched), the written Q becoming the
// posterior Q capped at 40. The pipeline contract pinned here: flag-off is
// byte-identical, a decisive disagreement flips exactly one base, a base one
// mate trimmed away can never reach (nor corrupt) its mate, and an insert the
// unique verifier refuses leaves everything exactly as sequenced.
// ---------------------------------------------------------------------------

/// FASTQ text with explicit per-base qualities: the rule acts on quality, so
/// the fixtures need to control it (`pair_fastq` writes Q40 throughout).
fn pair_fastq_q(p1: &[u8], q1: &[u8], p2: &[u8], q2: &[u8], name: &str) -> (String, String) {
    let rec = |id: &str, s: &[u8], q: &[u8]| {
        assert_eq!(s.len(), q.len(), "seq/qual length mismatch in fixture");
        format!(
            "@{id}\n{}\n+\n{}\n",
            String::from_utf8_lossy(s),
            String::from_utf8_lossy(q)
        )
    };
    (rec(name, p1, q1), rec(name, p2, q2))
}

/// One FASTQ record body back to its `(seq, qual)` halves.
fn seq_qual(body: &str) -> (Vec<u8>, Vec<u8>) {
    let lines: Vec<&str> = body.lines().collect();
    assert_eq!(lines.len(), 4, "one FASTQ record expected, got: {body:?}");
    (lines[1].as_bytes().to_vec(), lines[3].as_bytes().to_vec())
}

#[test]
fn pe_correct_rewrites_the_low_q_error_and_rescores_it() {
    // insert 100 on 150-base mates: the shared fragment window is [0, 100).
    // At fragment position 50 R1 reports the WRONG base at Q2 ('#') while R2
    // reports the correct one at Q40 ('I') -- the posterior is decisive for
    // R2, so R1's base is rewritten and its Q rescored to the capped posterior
    // Q ('I' = Phred 40). Everything else must be untouched. The posterior
    // acts on QUALITY, not truth: the wrong call must be the low-Q one.
    let (mut r1, r2) = pair(150, 100, 77);
    let fpos = 50usize;
    let j2 = 100 - 1 - fpos; // 49: R2's read-orientation partner of fpos
    let truth = r1[fpos];
    r1[fpos] = if truth == b'A' { b'C' } else { b'A' };
    assert_eq!(comp(r2[j2]), truth, "fixture: R2 must carry the true base");
    let mut q1 = vec![b'I'; 150];
    q1[fpos] = b'#'; // Phred 2: confidently wrong
    let q2 = vec![b'I'; 150]; // R2 correct at Phred 40

    let in1 = tmp_path("pc_r1.fq");
    let in2 = tmp_path("pc_r2.fq");
    let (a, b) = pair_fastq_q(&r1, &q1, &r2, &q2, "pc");
    write_plain(&in1, &a);
    write_plain(&in2, &b);
    let p = params();

    // Flag OFF: byte-identical to the input (nothing trims on this fixture).
    let off1 = tmp_path("pc_off_1.fq");
    let off2 = tmp_path("pc_off_2.fq");
    run_paired(
        &in1, &in2, &off1, &off2, None, None, &p, &p,
        &runtime(PeOverlapMode::Never),
    )
    .expect("paired run must succeed");
    assert_eq!(read_plain(&off1), a, "flag off: R1 must be byte-identical");
    assert_eq!(read_plain(&off2), b, "flag off: R2 must be byte-identical");

    // Flag ON: exactly one base + one quality byte change, on R1 only.
    let mut rt = runtime(PeOverlapMode::Never);
    rt.pe_correct = Some(2.0);
    rt.pe_correct_min_overlap = 30;
    let on1 = tmp_path("pc_on_1.fq");
    let on2 = tmp_path("pc_on_2.fq");
    let st = run_paired(&in1, &in2, &on1, &on2, None, None, &p, &p, &rt)
        .expect("paired run must succeed");
    assert_eq!(st.pairs_out, 1, "the pair survives");
    assert_eq!(st.bases_corrected, 1, "exactly one base was rewritten");

    let (o1s, o1q) = seq_qual(&read_plain(&on1));
    let (o2s, o2q) = seq_qual(&read_plain(&on2));
    assert_eq!(o1s[fpos], truth, "the low-Q error must become the true base");
    assert_eq!(
        o1q[fpos], b'I',
        "the rewritten Q is the capped posterior Q (Phred 40 -> 'I')"
    );
    let mut r1_exp = r1.clone();
    r1_exp[fpos] = truth;
    let mut q1_exp = q1.clone();
    q1_exp[fpos] = b'I';
    assert_eq!(o1s, r1_exp, "no other R1 base may move");
    assert_eq!(o1q, q1_exp, "no other R1 quality may move");
    assert_eq!(o2s, r2, "the winning mate's bases are untouched");
    assert_eq!(o2q, q2, "the winning mate's qualities are untouched");

    for f in [&in1, &in2, &off1, &off2, &on1, &on2] {
        std::fs::remove_file(f).ok();
    }
}

#[test]
fn a_base_one_mate_trims_away_must_not_rewrite_its_mate() {
    // Both mates hard-trim their 3' 60 bases (trim_3 = 60), so with insert 100
    // the doubly-retained overlap is fragment positions 10..90. Position 95 is
    // SHARED but only R2 keeps its partner (R2 index 4): R1's base there is
    // trimmed-away junk. Give that junk the WRONG call and the HIGHER quality
    // -- unmasked it would win the posterior and corrupt R2's kept base. The
    // masking contract: output is byte-identical to the flag-off run.
    let (mut r1, r2) = pair(150, 100, 91);
    let fpos = 95usize;
    let j2 = 100 - 1 - fpos; // 4
    let truth = r1[fpos];
    assert_eq!(
        comp(r2[j2]),
        truth,
        "fixture: R2's kept base must be the true call"
    );
    r1[fpos] = if truth == b'A' { b'C' } else { b'A' }; // wrong, and trimmed away on R1
    let q1 = vec![b'I'; 150]; // Q40 junk: would win the posterior if it counted
    let mut q2 = vec![b'I'; 150];
    q2[j2] = b'+'; // '+' = Phred 10: R2's correct call is the weaker voice

    let p = fastprep::algorithm::Params {
        trim_3: Some(60),
        ..params()
    };
    let in1 = tmp_path("pc2_r1.fq");
    let in2 = tmp_path("pc2_r2.fq");
    let (a, b) = pair_fastq_q(&r1, &q1, &r2, &q2, "pc2");
    write_plain(&in1, &a);
    write_plain(&in2, &b);

    let mut rt = runtime(PeOverlapMode::Never);
    rt.pe_correct = Some(2.0);
    let on1 = tmp_path("pc2_on_1.fq");
    let on2 = tmp_path("pc2_on_2.fq");
    let st = run_paired(&in1, &in2, &on1, &on2, None, None, &p, &p, &rt)
        .expect("paired run must succeed");
    assert_eq!(st.bases_corrected, 0, "nothing may be rewritten here");

    let off1 = tmp_path("pc2_off_1.fq");
    let off2 = tmp_path("pc2_off_2.fq");
    run_paired(
        &in1, &in2, &off1, &off2, None, None, &p, &p,
        &runtime(PeOverlapMode::Never),
    )
    .expect("paired run must succeed");

    // Self-check the fixture geometry before the real assertions.
    assert_eq!(
        read_plain(&off1).lines().nth(1).unwrap().len(),
        90,
        "fixture: trim_3 = 60 must leave 90 bases"
    );
    assert_eq!(
        read_plain(&on1),
        read_plain(&off1),
        "R1: nothing to correct inside its kept window"
    );
    assert_eq!(
        read_plain(&on2),
        read_plain(&off2),
        "R2: the trimmed mate's base must never reach it"
    );
    let (o2s, o2q) = seq_qual(&read_plain(&on2));
    assert_eq!(o2s[j2], comp(truth), "R2's kept base stays the true call");
    assert_eq!(o2q[j2], b'+', "and keeps its honest quality");

    for f in [&in1, &in2, &on1, &on2, &off1, &off2] {
        std::fs::remove_file(f).ok();
    }
}

#[test]
fn an_uninferrable_insert_is_left_untouched() {
    // A period-3 fragment fits MANY inserts (every candidate whose shift is a
    // multiple of 3 verifies on this geometry), so the unique verifier refuses
    // rather than guessing. Even a decisive disagreement must then stay
    // exactly as sequenced -- censoring beats fabrication.
    let unit = b"ACG";
    let r1_clean: Vec<u8> = (0..60).map(|i| unit[i % 3]).collect();
    let r2: Vec<u8> = r1_clean.iter().rev().map(|&b| comp(b)).collect();
    let mut r1 = r1_clean;
    let truth = r1[30];
    assert_eq!(truth, b'A', "fixture: position 30 is unit[0]");
    r1[30] = b'T'; // decisive: wrong call, and no new homopolymer run
    let mut q1 = vec![b'I'; 60];
    q1[30] = b'#'; // Phred 2: would lose the posterior outright
    let q2 = vec![b'I'; 60];

    let in1 = tmp_path("pc3_r1.fq");
    let in2 = tmp_path("pc3_r2.fq");
    let (a, b) = pair_fastq_q(&r1, &q1, &r2, &q2, "pc3");
    write_plain(&in1, &a);
    write_plain(&in2, &b);

    let mut rt = runtime(PeOverlapMode::Never);
    rt.pe_correct = Some(2.0);
    rt.pe_correct_min_overlap = 30;
    let on1 = tmp_path("pc3_on_1.fq");
    let on2 = tmp_path("pc3_on_2.fq");
    let st = run_paired(&in1, &in2, &on1, &on2, None, None, &params(), &params(), &rt)
        .expect("paired run must succeed");
    assert_eq!(
        st.bases_corrected, 0,
        "an ambiguous insert must never be corrected"
    );
    assert_eq!(read_plain(&on1), a, "R1 stays exactly as sequenced");
    assert_eq!(read_plain(&on2), b, "R2 stays exactly as sequenced");

    for f in [&in1, &in2, &on1, &on2] {
        std::fs::remove_file(f).ok();
    }
}
