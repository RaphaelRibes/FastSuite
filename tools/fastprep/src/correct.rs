//! Posterior PE-overlap consensus correction.
//!
//! When the two mates of a pair overlap, they read the *same* fragment from
//! opposite strands, so at every shared fragment position they carry two
//! independent observations of one true base. A disagreement between them is
//! almost always one mate being wrong. This module picks the better-supported
//! base and rewrites **only the losing mate**, keeping the paired-end record
//! structure intact (unlike a merge).
//!
//! The rule is the two-hypothesis posterior from `bench/scripts/
//! consensus_headroom.py`, which was measured against `fastp --correction`
//! (its rigid "one mate ≥Q30 and the other ≤Q14" rule) on a truth-exact
//! corpus and dominated it on both axes at every gate — and, decisively, did
//! **not** degrade faster than the rigid rule when the quality strings lie
//! (the `--err-scale 2.0` stress). See `docs/plans/` for the decision record.
//!
//! # Model
//!
//! At a shared position the mates report bases `b1` (R1, already in fragment
//! orientation) and `b2 = complement(R2[j])` (R2 flipped into fragment
//! orientation). Let `p_i = 10^(-q_i/10)` be mate *i*'s per-base error
//! probability. The posterior odds that `b1` is the true base rather than `b2`
//! is
//!
//! ```text
//!   odds = (1-p1)·p2 / ((1-p2)·p1)        (the uniform wrong-base /3 factors cancel)
//!   log10(odds) = log10((1-p1)·p2) - log10((1-p2)·p1)
//! ```
//!
//! We act only when `|log10 odds|` reaches `gate`, rewriting the losing mate's
//! base and setting its quality to the posterior quality of the consensus,
//! `min(40, floor(10·log10(1 + 10^|log10 odds|)))` — an evidence-based rescore,
//! never inflation.
//!
//! # Guarantees
//!
//! * **Ties do nothing.** `log10 odds == 0` (equal quality evidence on both
//!   mates) never fires, so a symmetric disagreement is never resolved by fiat.
//! * **At most one base changes per position** — only the loser — so the
//!   corruption ledger stays auditable.
//! * **`N` and other non-ACGT bases are skipped**; there is no evidence to
//!   arbitrate on.
//! * This is a pure function over the two mates. It never touches trim
//!   coordinates; the pipeline calls it *after* trimming and adapter/overlap
//!   decisions, so "same trimming decisions, fewer wrong bases" holds by
//!   construction.

/// The consensus quality is capped here. A capped score keeps the rescore
/// honest — the posterior can be extremely confident on a decisive two-mate
/// disagreement, but a single capped column is the convention downstream
/// tools expect and it never overstates the evidence.
pub const MAX_Q: u8 = 40;

/// Complement a single base, uppercase in → uppercase out. `N` (and any
/// non-ACGT byte) maps to `N`, which can never satisfy a comparison against a
/// real base.
#[inline]
pub fn comp_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        _ => b'N',
    }
}

#[inline]
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't')
}

/// Case-insensitive base equality, so a soft-masked (lowercase) call of the
/// same base agrees with its uppercase mate. `N` never equals a real base.
#[inline]
fn base_eq(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// Infer the fragment insert `I` from a mate pair, for consensus correction.
///
/// For a candidate `I` the shared fragment window is `[max(0, I - l2), min(I,
/// l1))`; fragment position `fpos` there pairs R1 index `fpos` with R2 index
/// `I - 1 - fpos`. The true `I` makes the two mates agree over that window up
/// to sequencing errors (a low mismatch count); a wrong `I` pairs unrelated
/// bases (a high mismatch count). Returns the insert only when **exactly one**
/// candidate clears the mismatch budget — a repeated fragment that verifies at
/// several inserts is refused rather than guessed (the suite's censoring
/// discipline). `None` means "do not correct", which leaves the trim path
/// byte-for-byte untouched.
///
/// This handles both regimes the correction must cover: `I <= L` (full
/// overlap plus adapter readthrough) and `L < I < 2L` (partial overlap, each
/// mate reading past the fragment midpoint). The existing
/// [`crate::pe_overlap`] `detect_pe_overlap` finds only the readthrough case
/// (`I <= min(l1, l2)`), so it cannot on its own reach the partial-overlap
/// error mass.
pub fn infer_insert(r1: &[u8], r2: &[u8], min_ovl: usize) -> Option<usize> {
    let (l1, l2) = (r1.len(), r2.len());
    if l1 == 0 || l2 == 0 || min_ovl == 0 {
        return None;
    }
    // `min_ovl` is user-controlled (`--pe-correct-min-overlap`) and the guard
    // above only rejects zero: when it exceeds `l1 + l2` this subtraction
    // wrapped to ~1.8e19 in release (the loop then spun forever on the first
    // pair) and panicked under `overflow-checks` in debug.
    let max_i = (l1 + l2).checked_sub(min_ovl)?;
    let mut found: Option<usize> = None;
    for i in min_ovl..=max_i {
        let lo = i.saturating_sub(l2);
        let hi = std::cmp::min(i, l1);
        if hi <= lo || hi - lo < min_ovl {
            continue;
        }
        let shared = hi - lo;
        // A wrong insert pairs unrelated bases (≈ 3/4 mismatch); the true one
        // mismatches only where a sequencing error split the mates (≈ 1/70).
        // shared/6 sits in the wide gap between them.
        let tol = std::cmp::max(2, shared / 6);
        let mut mm = 0usize;
        let mut compared = 0usize;
        let mut over = false;
        for fpos in lo..hi {
            let j2 = i - 1 - fpos;
            if j2 >= l2 {
                over = true;
                break;
            }
            let b1 = r1[fpos];
            let r2b = r2[j2];
            if !is_acgt(b1) || !is_acgt(r2b) {
                continue; // ambiguous base: no evidence either way
            }
            compared += 1;
            if !base_eq(b1, comp_base(r2b)) {
                mm += 1;
                if mm > tol {
                    over = true;
                    break;
                }
            }
        }
        // A candidate needs at least `min_ovl` *informative* (both-ACGT)
        // positions. Without this an all-`N` window trivially "verifies" with
        // zero mismatches and would collide with the true insert, which the
        // ambiguity refusal would then throw away — turning a confident call
        // into `None`.
        if over || compared < min_ovl || mm > tol {
            continue;
        }
        if found.is_some() {
            return None; // ambiguous: more than one insert verifies
        }
        found = Some(i);
    }
    found
}

/// Posterior-gated consensus correction over the mate overlap.
///
/// `s1`/`q1` are R1's bases and raw quality **bytes**; `s2`/`q2` are R2's, in
/// R2's own read orientation exactly as FASTQ stores it. `phred_offset` is the
/// encoding offset of the quality bytes (33 for standard Phred+33, 0 when the
/// caller passes already-decoded Phred scores). The loser's quality byte is
/// rewritten as `qn + phred_offset`.
///
/// `insert` is the fragment length `I`. The shared fragment positions are
/// `fpos ∈ [max(0, I - l2), min(I, l1))`; `fpos` maps to R1 index `fpos` and R2
/// index `I - 1 - fpos`. `min_ovl` is the minimum shared length below which the
/// pair is left untouched (same filter the validated census applied).
///
/// `gate` is the minimum `|log10 odds|` required to act (`>= 0`). A tie — both
/// mates carrying exactly equal evidence — always has `log10 odds == 0` and is
/// skipped regardless of the gate, so a symmetric disagreement is never
/// resolved by fiat; the caller may pass `0` to resolve every non-tie
/// disagreement.
///
/// Returns the number of bases rewritten (≤ one per shared position).
#[allow(clippy::too_many_arguments)]
pub fn correct_overlap_posterior(
    s1: &mut [u8],
    q1: &mut [u8],
    s2: &mut [u8],
    q2: &mut [u8],
    phred_offset: u8,
    insert: usize,
    gate: f64,
    min_ovl: usize,
) -> usize {
    let (l1, l2) = (s1.len(), s2.len());
    debug_assert_eq!(l1, q1.len(), "R1 seq/qual length mismatch");
    debug_assert_eq!(l2, q2.len(), "R2 seq/qual length mismatch");
    debug_assert!(gate >= 0.0, "gate must be >= 0");

    // Shared fragment window: both mates have a base here.
    let lo = insert.saturating_sub(l2);
    let hi = std::cmp::min(insert, l1);
    if hi <= lo || hi - lo < min_ovl {
        return 0;
    }

    let mut edits = 0usize;
    for fpos in lo..hi {
        let j1 = fpos;
        let j2 = insert - 1 - fpos;
        if j2 >= l2 {
            continue; // belt and suspenders; `lo` already keeps j2 in range
        }
        let b1 = s1[j1];
        let r2 = s2[j2];
        // No evidence to arbitrate on: an ambiguous base is never rewritten.
        if !is_acgt(b1) || !is_acgt(r2) {
            continue;
        }
        let b2 = comp_base(r2); // mate 2 into fragment orientation
        if base_eq(b1, b2) {
            continue; // agreement is left untouched
        }

        // Q0 means "no confidence", not "certain to be right": at Q0 the
        // decoded probability is 1.0, which makes num or den zero and the
        // position is skipped — the most decisive disagreement (a Q0 loser
        // against a high-Q winner) was the one never corrected. Clamp the
        // decoded Q at 1 so the posterior stays positive.
        let dq1 = q1[j1].saturating_sub(phred_offset).max(1);
        let dq2 = q2[j2].saturating_sub(phred_offset).max(1);
        let p1 = 10f64.powf(-(dq1 as f64) / 10.0);
        let p2 = 10f64.powf(-(dq2 as f64) / 10.0);
        let num = (1.0 - p1) * p2;
        let den = (1.0 - p2) * p1;
        if num <= 0.0 || den <= 0.0 {
            continue; // degenerate posterior; no basis to act
        }
        let lodds = num.log10() - den.log10();
        let a = lodds.abs();
        // A tie (equal evidence) is a no-op, and anything under the gate is a
        // no-op. Both are deliberate: the gate does all the safety work.
        if lodds == 0.0 || a < gate {
            continue;
        }

        let qn = std::cmp::min(
            MAX_Q,
            (10.0f64 * (1.0 + 10f64.powf(a)).log10()).floor() as u8,
        );
        if lodds > 0.0 {
            // Mate 1 wins; mate 2 loses. Rewrite mate 2 in *its* read
            // orientation, i.e. the complement of the winning fragment base.
            s2[j2] = comp_base(b1);
            q2[j2] = qn.saturating_add(phred_offset);
        } else {
            // Mate 2 wins; mate 1 loses. R1 is already fragment orientation.
            s1[j1] = b2;
            q1[j1] = qn.saturating_add(phred_offset);
        }
        edits += 1;
    }
    edits
}

/// Mask every position outside `keep` to `N` (base) and `0` (quality), in
/// place. A masked base is non-ACGT, so [`correct_overlap_posterior`] skips the
/// position outright — a trimmed-away (adapter / poly-X / low-Q) base can
/// therefore never supply evidence to rewrite its mate's *retained* base, and
/// is never itself rewritten.
fn mask_outside(seq: &mut [u8], qual: &mut [u8], keep: &std::ops::Range<usize>) {
    for i in 0..seq.len() {
        if !keep.contains(&i) {
            seq[i] = b'N';
            qual[i] = 0;
        }
    }
}

/// Consensus-correct the mate overlap over each mate's *retained* window.
///
/// `keep_i` is the half-open range of read positions mate `i` actually emits
/// (its `t5..len-t3` window). Positions outside `keep_i` are masked to `N`
/// first (see [`mask_outside`]), so the correction acts only where *both* mates
/// keep a base — exactly the doubly-retained overlap — and trimmed-away junk
/// can never corrupt a retained base. Coordinates are full-read (index 0 is the
/// read's 5′ end), so `insert` maps fragment position `fpos` to `s1[fpos]` and
/// `s2[insert-1-fpos]` exactly regardless of how much was trimmed; the caller
/// slices `keep_i` out afterwards.
///
/// Returns the number of retained bases rewritten (≤ one per shared position).
#[allow(clippy::too_many_arguments)]
pub fn correct_overlap_retained(
    s1: &mut [u8],
    q1: &mut [u8],
    s2: &mut [u8],
    q2: &mut [u8],
    keep1: std::ops::Range<usize>,
    keep2: std::ops::Range<usize>,
    phred_offset: u8,
    insert: usize,
    gate: f64,
    min_ovl: usize,
) -> usize {
    // The floor belongs on the DOUBLY-RETAINED evidence, not on the geometric
    // overlap. With a 32-base geometric window but one retained base on each
    // side, `min_ovl = 30` still rewrote that base's mate — a single retained
    // base cannot outvote anything, and the caller passed 30 because it wanted
    // 30 bases of evidence.
    let lo = insert.saturating_sub(s2.len());
    let hi = std::cmp::min(insert, s1.len());
    let retained = (lo..hi)
        .filter(|&fpos| {
            let j1 = fpos;
            let j2 = insert - 1 - fpos;
            j2 < s2.len() && keep1.contains(&j1) && keep2.contains(&j2)
        })
        .count();
    if retained < min_ovl {
        return 0;
    }
    mask_outside(s1, q1, &keep1);
    mask_outside(s2, q2, &keep2);
    correct_overlap_posterior(s1, q1, s2, q2, phred_offset, insert, gate, min_ovl)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a full-overlap pair: insert = read length, R1 = fragment,
    // R2 = reverse-complement(fragment). Then force a single disagreement.
    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| comp_base(b)).collect()
    }

    fn pair_from(fragment: &[u8]) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, usize) {
        let s1 = fragment.to_vec();
        let s2 = rc(fragment);
        let n = fragment.len();
        let q1 = vec![30u8; n];
        let q2 = vec![30u8; n];
        (s1, q1, s2, q2, n)
    }

    #[test]
    fn complement_mapping() {
        assert_eq!(comp_base(b'A'), b'T');
        assert_eq!(comp_base(b'T'), b'A');
        assert_eq!(comp_base(b'G'), b'C');
        assert_eq!(comp_base(b'C'), b'G');
        assert_eq!(comp_base(b'N'), b'N');
        assert_eq!(comp_base(b'a'), b'T'); // case-insensitive, uppercase out
    }

    #[test]
    fn agreement_is_left_untouched() {
        let frag = b"ACGTACGTACGTACGTACGTACGTACGTAC";
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(frag);
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 0.5, 10);
        assert_eq!(e, 0, "an agreeing pair must never be edited");
        assert_eq!(s1, frag.to_vec());
    }

    #[test]
    fn tie_does_nothing_even_at_zero_gate() {
        // Force a disagreement at a shared position; give both mates equal
        // quality so the posterior is exactly 1 (log-odds 0). Must not fire.
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A'; // R1 will report A here
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        // Corrupt R2's corresponding base so the mates disagree, keeping equal Q.
        // R2 index for fpos=10 is ins-1-10 = 18.
        s2[19] = comp_base(b'T'); // R2 now reads fragment base T at fpos 10
        q1[10] = 30;
        q2[19] = 30;
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 0.0, 1);
        assert_eq!(e, 0, "a tie must never fire, even at gate 0");
    }

    #[test]
    fn higher_quality_mate_writes_lower_quality_mate() {
        // R1 says A at Q40, R2 (fragment orientation) says T at Q2. Mate 1 wins;
        // mate 2 is rewritten to A, with a capped posterior Q.
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        // Make the mates disagree at fpos=10 (R2 index 18).
        s2[19] = comp_base(b'T'); // fragment T at fpos 10
        q1[10] = 40;
        q2[19] = 2;
        let before_s1 = s1.clone();
        let before_q1 = q1.clone();
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 0.5, 1);
        assert_eq!(e, 1);
        // Winner untouched.
        assert_eq!(s1, before_s1);
        assert_eq!(q1, before_q1);
        // Loser rewritten to the winning base (complemented back to R2 frame)
        // with a capped, evidence-based quality.
        assert_eq!(
            s2[19],
            comp_base(b'A'),
            "loser base must become the winner's"
        );
        assert!(
            q2[19] > 30 && q2[19] <= MAX_Q,
            "rescored Q out of range: {}",
            q2[19]
        );
    }

    #[test]
    fn q_out_is_capped() {
        // Two extremely confident but disagreeing mates: posterior Q must cap.
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        s2[19] = comp_base(b'T');
        q1[10] = 93;
        q2[19] = 1;
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 0.5, 1);
        assert_eq!(e, 1);
        let loser_q = q2[19];
        assert!(
            loser_q <= MAX_Q,
            "Q must be capped at {}: got {}",
            MAX_Q,
            loser_q
        );
    }

    #[test]
    fn n_is_skipped() {
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'N';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        s2[19] = comp_base(b'T'); // disagree with N
        q1[10] = 40;
        q2[19] = 2;
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 0.5, 1);
        assert_eq!(e, 0, "an N on either mate must be skipped");
    }

    #[test]
    fn below_gate_is_a_no_op() {
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        s2[19] = comp_base(b'T');
        q1[10] = 25;
        q2[19] = 24; // near-equal evidence → tiny |log-odds| → below any real gate
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 2.0, 1);
        assert_eq!(e, 0, "a weak disagreement below the gate must not fire");
    }

    #[test]
    fn min_overlap_filters_short_inserts() {
        // Insert 20 with min_ovl 30: the pair is left alone regardless.
        let frag = b"ACGTACGTACGTACGTACGT";
        let (mut s1, mut q1, mut s2, mut q2, _) = pair_from(frag);
        // Seed a strong disagreement at fpos=5 (R2 index = 20-1-5 = 14).
        s2[14] = comp_base(b'T');
        q1[5] = 40;
        q2[14] = 2;
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, 20, 0.5, 30);
        assert_eq!(e, 0, "short overlaps must be left untouched");
    }

    #[test]
    fn gate_symmetry_swapping_mates_is_mirrored() {
        // Mate 1 strong-A vs mate 2 weak-T, then the mirror image: the winner is
        // always the strong mate, and the edited base always lands on the weak one.
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A';

        // Case 1: mate1 strong.
        let (mut s1a, mut q1a, mut s2a, mut q2a, ins) = pair_from(&frag);
        s2a[19] = comp_base(b'T');
        q1a[10] = 40;
        q2a[19] = 2;
        correct_overlap_posterior(&mut s1a, &mut q1a, &mut s2a, &mut q2a, 0, ins, 0.5, 1);
        assert_eq!(
            s2a[19],
            comp_base(b'A'),
            "mate1 strong ⇒ mate2 rewritten to A"
        );

        // Case 2: mate2 strong (mirrored). Mate 2's fragment base is T (via
        // comp of its read base), mate 1 says A at low Q; mate2 wins, mate1
        // rewritten to T.
        let (mut s1b, mut q1b, mut s2b, mut q2b, ins) = pair_from(&frag);
        s2b[19] = comp_base(b'T'); // fragment T at fpos 10
        q1b[10] = 2;
        q2b[19] = 40;
        correct_overlap_posterior(&mut s1b, &mut q1b, &mut s2b, &mut q2b, 0, ins, 0.5, 1);
        assert_eq!(s1b[10], b'T', "mate2 strong ⇒ mate1 rewritten to T");
    }

    #[test]
    fn at_most_one_base_changes_per_position() {
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        s2[19] = comp_base(b'T');
        q1[10] = 40;
        q2[19] = 2;
        let edits = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 0, ins, 0.5, 1);
        // Exactly one base differs from the input across both mates.
        let mut changed = 0;
        if s1[10] != frag[10] {
            changed += 1;
        }
        if s2[19] != comp_base(b'T') {
            changed += 1;
        }
        assert_eq!(edits, 1);
        assert_eq!(
            changed, 1,
            "exactly one base must change per fired position"
        );
    }

    #[test]
    fn phred_offset_is_honoured() {
        // Same scenario but quality stored as Phred+33 bytes (offset 33).
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'A';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        s2[19] = comp_base(b'T');
        q1[10] = 40 + 33; // Q40 as a Phred+33 byte
        q2[19] = 2 + 33; // Q2 as a Phred+33 byte
        let e = correct_overlap_posterior(&mut s1, &mut q1, &mut s2, &mut q2, 33, ins, 0.5, 1);
        assert_eq!(e, 1);
        assert_eq!(s2[19], comp_base(b'A'));
        // The rewritten quality byte must be a valid Phred+33 encoding (≥ offset).
        assert!(q2[19] >= 33 && q2[19] <= MAX_Q + 33, "q byte {}", q2[19]);
    }

    // --- trimmed-base influence guard ---

    #[test]
    fn trimmed_away_base_never_rewrites_retained_mate() {
        // Mate 1's high-Q base at a shared fragment position is trimmed away
        // (outside keep1); mate 2's retained base there is correct. The
        // discarded mate-1 base must NOT be used to rewrite mate 2's retained
        // base — masking it to N makes `is_acgt` skip the position outright.
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'T'; // true fragment base at fpos 10 is T
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        // R2 index for fpos=10 is ins-1-10 = 19. pair_from already leaves mate 2
        // reading fragment base T there (correct). Give it low Q but keep it.
        q2[19] = 2;
        // Mate 1 reports A (wrong) at high Q — unmasked, the posterior would
        // pick it and rewrite mate 2's correct base (the corruption vector).
        s1[10] = b'A';
        q1[10] = 40;
        // Mate 1's position 10 is trimmed away; mate 2 is fully retained.
        let edits = correct_overlap_retained(
            &mut s1,
            &mut q1,
            &mut s2,
            &mut q2,
            0..10,
            0..30,
            0,
            ins,
            0.5,
            10,
        );
        // Mate 2's retained base must still read fragment T (unchanged).
        assert_eq!(
            s2[19],
            comp_base(b'T'),
            "trimmed-away mate-1 base rewrote retained mate-2 base"
        );
        // Mate 1's discarded base is masked to N so it can never act as evidence.
        assert_eq!(
            s1[10], b'N',
            "trimmed-away base must be masked, not trusted"
        );
        // No edit fired: the only disagreement involved a trimmed-away base.
        assert_eq!(edits, 0, "a trimmed-away base must not trigger any rewrite");
    }

    #[test]
    fn retained_mate_still_corrects_with_masking() {
        // Control for the influence guard: when BOTH mates retain a disagreeing
        // position, correction still fires exactly as before.
        let mut frag = *b"ACGTACGTACGTACGTACGTACGTACGTAC";
        frag[10] = b'T';
        let (mut s1, mut q1, mut s2, mut q2, ins) = pair_from(&frag);
        s1[10] = b'A'; // mate 1 wrong, LOW Q — must lose the posterior
        q1[10] = 2;
        q2[19] = 40; // mate 2 correct (fragment T), high Q — wins
        let edits = correct_overlap_retained(
            &mut s1,
            &mut q1,
            &mut s2,
            &mut q2,
            0..30,
            0..30,
            0,
            ins,
            0.5,
            10,
        );
        assert_eq!(edits, 1, "a doubly-retained disagreement must be corrected");
        // Mate 2 won (its fragment base T is correct) → mate 1 is rewritten to T.
        assert_eq!(
            s1[10], b'T',
            "retained wrong base must be corrected to the winner"
        );
    }

    // --- insert inference ---

    fn rand_frag(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    }

    /// Build a mate pair from a fragment and an insert size: R1 reads `F[0..]`,
    /// R2 reads the opposite strand so `R2[j] = complement(F[insert-1-j])`.
    /// Any read length beyond the insert is left as `N` (adapter / unsequenced).
    fn make_pair(frag: &[u8], insert: usize, l1: usize, l2: usize) -> (Vec<u8>, Vec<u8>) {
        let mut r1 = vec![b'N'; l1];
        let mut r2 = vec![b'N'; l2];
        for i in 0..l1.min(insert) {
            r1[i] = frag[i];
        }
        for j in 0..l2.min(insert) {
            r2[j] = comp_base(frag[insert - 1 - j]);
        }
        (r1, r2)
    }

    #[test]
    fn infer_insert_full_overlap() {
        let frag = rand_frag(60, 123);
        let (r1, r2) = make_pair(&frag, 60, 60, 60);
        assert_eq!(infer_insert(&r1, &r2, 30), Some(60));
    }

    #[test]
    fn infer_insert_partial_overlap() {
        // L < I < 2L: each mate reads past the fragment midpoint.
        let frag = rand_frag(90, 456);
        let (r1, r2) = make_pair(&frag, 90, 60, 60);
        assert_eq!(infer_insert(&r1, &r2, 30), Some(90));
    }

    #[test]
    fn infer_insert_readthrough_overlap() {
        // I < L: mates overlap fully and carry N padding past the insert.
        let frag = rand_frag(40, 789);
        let (r1, r2) = make_pair(&frag, 40, 60, 60);
        assert_eq!(infer_insert(&r1, &r2, 30), Some(40));
    }

    #[test]
    fn infer_insert_unrelated_is_none() {
        let r1 = rand_frag(60, 1);
        let r2 = rand_frag(60, 2);
        assert_eq!(infer_insert(&r1, &r2, 30), None);
    }

    #[test]
    fn infer_insert_repeated_fragment_is_refused() {
        // Periodic fragment: several inserts verify → refuse rather than guess.
        let frag: Vec<u8> = (0..60).map(|i| b"ACGT"[i % 4]).collect();
        let (r1, r2) = make_pair(&frag, 60, 60, 60);
        assert_eq!(infer_insert(&r1, &r2, 30), None);
    }
}
