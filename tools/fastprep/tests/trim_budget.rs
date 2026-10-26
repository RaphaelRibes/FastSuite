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
