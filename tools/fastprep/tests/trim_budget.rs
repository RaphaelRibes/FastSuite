//! Regression: a qualifying 3′ polyX run must be trimmed even when the
//! bounded reverse scan runs out of budget on the flank instead of on the run.
//!
//! Reproduced on the CLI at the shipped defaults (`--min-len 0.67L`,
//! `--min-run 12`, `--mercy 1`, `--head-dropout 4`): a 151 bp read whose 3′
//! end carries a 48-base polyA artefact was emitted **untrimmed** (151 bp,
//! `bases_trimmed_3 = 0`) while a 47-base artefact was trimmed — non-monotone
//! in artefact length. Cause: when the read's first base is a 1-base G-run,
//! `head_dropout` spends one base of the trim budget, so the reverse scan's
//! cap is exactly the run length; the walk then stops on `StopMode::Budget`
//! after the run ended, and `trim()` treated every Budget stop as "no run".

use fastprep::algorithm::*;

/// The shipped CLI defaults, copied from `fastprep --help`
/// (`[default: 3]` for `--mercy`, `[default: 0.67L]` for `--min-len`, ...).
/// NOT `Params::default()`, which is the library's own 0.6L / mercy 3 /
/// min_run 6 / no head dropout / not boundary-anchored. When a shipped default
/// changes, change it here too — this fixture exists to test the CLI path.
fn shipped() -> Params {
    let mut p = Params::default();
    p.min_len = MinLen::Frac(0.67);
    p.mercy = 3.0;
    p.phred_sat = 20;
    p.gate_shape = 1.0;
    p.match_leniency = 0.1;
    p.mismatch_strength = 0.5;
    p.min_run = 12;
    p.head_dropout = Some(4);
    p.scan_interior = false;
    p.boundary_anchored = true;
    p.prepare();
    p
}

/// 151 bp read, leading `G` (one-base dropout: head_dropout trims it and
/// spends one base of the budget), `run` trailing `A`s, Q40 everywhere.
/// The base immediately before the tail is a `C`, so the artefact is exactly
/// `run` bases long (a body ending in `A` would silently extend it).
fn read_with_tail(run: usize) -> (Vec<u8>, Vec<u8>) {
    let l = 151usize;
    let body = b"GTAC";
    let mut seq: Vec<u8> = Vec::with_capacity(l);
    seq.push(b'G'); // the dropout base, alone
    let mut i = 1usize;
    while seq.len() < l - run {
        seq.push(body[i % 4]);
        i += 1;
    }
    *seq.last_mut().expect("non-empty") = b'C';
    seq.extend(std::iter::repeat(b'A').take(run));
    assert_eq!(seq.len(), l);
    (seq, vec![40u8; l])
}

#[test]
fn a_run_at_the_budget_boundary_is_trimmed_not_kept() {
    let p = shipped();
    for run in [20usize, 47, 48, 49] {
        let (seq, q) = read_with_tail(run);
        let res = polyx_trim(&seq, &q, &p);
        match res {
            TrimResult::Emit { t3_poly, .. } => assert_eq!(
                t3_poly, run,
                "run {run}: expected the whole artefact trimmed, got {t3_poly} \
                 (0 means the artefact was emitted intact)"
            ),
            TrimResult::Drop { reason, .. } => assert_eq!(
                reason,
                DropReason::Poly3,
                "run {run}: only a budget-exhausting run may drop"
            ),
        }
    }
}

#[test]
fn a_tail_below_min_run_is_still_kept() {
    let p = shipped();
    let (seq, q) = read_with_tail(10);
    match polyx_trim(&seq, &q, &p) {
        TrimResult::Emit { t3_poly, .. } => assert_eq!(
            t3_poly, 0,
            "a 10-base run is below min_run 12 and must not be trimmed"
        ),
        TrimResult::Drop { reason, .. } => panic!("10-base tail dropped as {reason:?}"),
    }
}

/// 151 bp read: a `five`-base polyA 5′ run, a GTAC body, a `three`-base polyA
/// 3′ tail. The body's first base is C, so the 5′ run is exactly `five`; its
/// last base is A, so the 3′ artefact is `three + 1` bases long — that is what
/// makes the reviewed 48-base control trim 49.
fn read_151(five: usize, three: usize) -> (Vec<u8>, Vec<u8>) {
    let l = 151usize;
    let body_len = l - five - three;
    let pat = b"GTAC";
    let mut seq: Vec<u8> = Vec::with_capacity(l);
    seq.resize(five, b'A');
    for i in 0..body_len {
        seq.push(pat[i % 4]);
    }
    seq[five] = b'C';
    seq[five + body_len - 1] = b'A';
    seq.resize(l, b'A');
    assert_eq!(seq.len(), l);
    (seq, vec![40u8; l])
}

/// `Some((emitted, t5, t3))` for an emitted read, `None` when it dropped.
fn outcome(p: &Params, five: usize, three: usize) -> Option<(usize, usize, usize)> {
    let (seq, q) = read_151(five, three);
    match polyx_trim(&seq, &q, p) {
        TrimResult::Emit {
            t5_poly, t3_poly, ..
        } => Some((seq.len() - t5_poly - t3_poly, t5_poly, t3_poly)),
        TrimResult::Drop { .. } => None,
    }
}

/// The four rows reproduced on the CLI at the shipped defaults (151 bp,
/// Δ_L = 50), which the library must reproduce because `shipped()` is exactly
/// that configuration:
///
/// | read            | emitted | t5 | t3   |
/// |-----------------|---------|----|------|
/// | 5′40 only       | 111     | 40 | 0    |
/// | 3′60 only       | drop    |    |      |
/// | 5′40 + 3′60     | 111     | 40 | **0**|  ← the defect: 60 bp kept whole
/// | 3′48 (control)  | 102     | 0  | 49   |
///
/// The 5′ cut spends 40 of the 50-base budget, leaving the reverse walk 10
/// steps for a 60-base tail. It stops on `StopMode::Budget` after 10 matches,
/// and the old qualifying test (`matches + 1 >= min_run`) is unsatisfiable
/// below `min_run - 1` reach, so the real run resolved to "no run".
#[test]
fn a_three_prime_artefact_is_not_kept_when_the_five_prime_one_spent_the_budget() {
    let p = shipped();
    assert_eq!(outcome(&p, 40, 0), Some((111, 40, 0)), "5′40 only");
    assert_eq!(outcome(&p, 0, 60), None, "3′60 only must drop");
    assert_eq!(
        outcome(&p, 40, 60),
        None,
        "5′40 + 3′60: the 60-base tail must not be emitted whole"
    );
    assert_eq!(outcome(&p, 0, 48), Some((102, 0, 49)), "3′48 control");

    // The chosen outcome for the two-artefact row is the single-sided 3′60
    // row's fate: a real tail the remaining budget cannot cover drops the read
    // as Poly3, never emits it with the tail intact.
    let (seq, q) = read_151(40, 60);
    match polyx_trim(&seq, &q, &p) {
        TrimResult::Drop { reason, .. } => assert_eq!(reason, DropReason::Poly3),
        TrimResult::Emit { t3_poly, .. } => {
            panic!("emitted with t3 = {t3_poly}; the tail was kept whole")
        }
    }
}

/// The reverse pass owns the 5′ end under `--r1-scan-reverse`. Same defect,
/// other physical end: the streaming pass spends 41 of the 50-base budget on
/// the 3′ run, leaving the reverse pass 9 steps for a 60-base 5′ run.
#[test]
fn the_reverse_pass_cuts_a_five_prime_run_the_same_way() {
    let mut p = shipped();
    p.direction = Direction::Reverse;
    let (seq, q) = read_151(60, 40);
    match polyx_trim(&seq, &q, &p) {
        // Poly5 is the reverse pass reporting a run it could not cover; before
        // the fix this emitted the 151 bp read with t5 = 0.
        TrimResult::Drop { reason, .. } => assert_eq!(reason, DropReason::Poly5),
        TrimResult::Emit { t5_poly, .. } => {
            panic!("5′ tail emitted with t5 = {t5_poly}")
        }
    }
}

#[test]
fn a_tail_longer_than_the_budget_still_drops() {
    let p = shipped();
    let (seq, q) = read_with_tail(60);
    match polyx_trim(&seq, &q, &p) {
        TrimResult::Drop { reason, .. } => assert_eq!(reason, DropReason::Poly3),
        TrimResult::Emit { t3_poly, .. } => {
            panic!("a 60-base tail on a 151 bp read must drop, trimmed {t3_poly}")
        }
    }
}
