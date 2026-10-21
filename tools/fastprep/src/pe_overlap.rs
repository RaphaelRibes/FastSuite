//! Paired-end overlap detection (paper §`sec:adp-pe`).
//!
//! When the insert is shorter than the read, R1 reads the fragment 5′→3′ while
//! R2 reads the same fragment from the opposite strand. Writing `I` for the
//! insert size and `L1`/`L2` for the read lengths:
//!
//! ```text
//!   R1 = F[0..I]        + adapter_top[0..L1-I]
//!   R2 = revcomp(F[0..I]) + adapter_bottom[0..L2-I]
//! ```
//!
//! so the two mates agree on their first `I` bases once R2 is put in the
//! fragment's own frame:
//!
//! ```text
//!   R1[i] == complement(R2[I-1-i])      for i < I
//! ```
//!
//! That overlap pins the insert size `I`, and therefore the boundary at which
//! each mate leaves the insert and enters its adapter: the boundary falls at
//! position `I` on **both** mates (each carries its own adapter ligated at the
//! far end, read 3′ on both). The mate pair thus determines the readthrough
//! length that single-read sequence matching can only guess at, and it needs
//! no adapter sequence at all.
//!
//! # Why this is anchored, not the paper's scan
//!
//! The paper's `@alg:adp-pe` scans insert sizes `I` from `min(L1, L2)` down to
//! `o_PE`, tests `R1[i] == complement(R2[I-1-i])` for `i < I`, and accepts the
//! first `I` whose mismatch *fraction* is within `ε`. Taken literally on real
//! data that fires constantly and wrongly: at `I = min(L1, L2)` the comparison
//! is between unrelated sequence, and at `ε = 0.1` a 150-base window admits 15
//! mismatches, so a chance alignment passes at essentially every offset. The
//! paper itself flags the cost of that scan and notes that a probe gate would
//! be "not sound".
//!
//! This module therefore follows the suite's own discipline —
//! `tools/fastcheck/src/insert.rs` (`ANCHOR = 12`, full-window verification,
//! censoring instead of guessing) — with the 12-base anchor the task names as
//! the natural minimum. The contract:
//!
//! * a candidate insert size is admitted only when the 12-base anchor at the
//!   3′ end of its shared window matches (`R1[I-1-j] == complement(R2[j])`
//!   for `j < 12`) — a chance hit in 250 bases is ~1e-5;
//! * the whole shared window is then verified base for base, up to
//!   [`MAX_OVERLAP_MISMATCHES`];
//! * if more than one insert size verifies, the alignment is **ambiguous** and
//!   the detector returns nothing rather than choosing (a repeated fragment is
//!   exactly this case);
//! * if the mates cannot overlap, no anchor matches, or the only matches are
//!   chance hits, it returns nothing. A chance match is never accepted (the
//!   discipline of fastcheck's
//!   `a_chance_anchor_hit_does_not_produce_a_measurement`).
//!
//! `None` is the caller's signal to fall back to the ordinary per-mate
//! behaviour — which is exactly what happens when `--pe-overlap` is off.

/// Bases of the mate used to locate the overlap. Twelve keeps a spurious hit
/// in a 250-base read at roughly 1e-5 while staying short enough that real
/// overlaps are still found. Same value as
/// `tools/fastcheck/src/insert.rs::ANCHOR`.
pub const ANCHOR: usize = 12;

/// Mismatches tolerated inside the anchor. One base of slack lets a sequencing
/// error (or the artefact boundary) inside the 12 anchor bases survive without
/// dropping a pair the full verification would have accepted.
pub const MAX_ANCHOR_MISMATCHES: usize = 1;

/// Mismatches tolerated across the whole verified window. Not zero, and
/// deliberately so: a real pair disagrees at the base adjacent to an artefact
/// by construction (each mate's boundary base is called independently), and
/// sequencing errors add a few more. A chance alignment differs at a quarter of
/// positions, so 3 rejects it over any window of 12 bases or more while
/// absorbing the real disagreements. Same value as
/// `tools/fastcheck/src/insert.rs::MAX_OVERLAP_MISMATCHES`.
pub const MAX_OVERLAP_MISMATCHES: usize = 3;

#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        _ => b'N',
    }
}

/// Case-insensitive base equality. Bases are compared as bytes; the comparison
/// is what makes a lowercase or upper-case call of the same base agree, while
/// `N` (which is what `comp` returns for anything outside ACGT) can never equal
/// a real base.
#[inline]
fn base_eq(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// Reverse complement, uppercase output. `N` (and any non-ACGT byte) maps to
/// `N`, which can never satisfy a comparison against a real base, so an
/// ambiguous base can only cost a pair a detection, never create a false one.
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

/// True when two equal-length slices agree within [`MAX_ANCHOR_MISMATCHES`]
/// mismatches, case-insensitively.
#[inline]
fn anchor_matches(a: &[u8], b: &[u8]) -> bool {
    debug_assert_eq!(a.len(), b.len());
    let mut mm = 0usize;
    for i in 0..a.len() {
        if !base_eq(a[i], b[i]) {
            mm += 1;
            if mm > MAX_ANCHOR_MISMATCHES {
                return false;
            }
        }
    }
    true
}

/// Mismatches over the candidate shared window.
///
/// The shared region is R1's PREFIX of length `o` (`r1[0..o]`) against the
/// SUFFIX of `revcomp(R2)` (`rc2[L2-o..L2]`): R1 reads the fragment 5′→3′
/// while `revcomp(R2)` is the same fragment read from the other strand, so the
/// two agree on their first `o` bases. This is the paper's identity
/// `R1[i] == complement(R2[I-1-i])` for `i < I`.
#[inline]
fn window_mismatches(r1: &[u8], rc2: &[u8], o: usize) -> usize {
    let l2 = rc2.len();
    let mut mm = 0usize;
    for j in 0..o {
        if !base_eq(r1[j], rc2[l2 - o + j]) {
            mm += 1;
            if mm > MAX_OVERLAP_MISMATCHES {
                return mm; // saturated: caller only tests `<=`
            }
        }
    }
    mm
}

/// Detect the paired-end overlap and the readthrough it implies.
///
/// `r1`, `r2` are the two mates as read (R2 in its own 5′→3′ orientation, as
/// FASTQ stores it). `anchor` is the minimum overlap the caller is willing to
/// act on; it is clamped up to [`ANCHOR`], the length actually used to locate
/// the alignment. Returns `Some((mate, readthrough))` where `mate` is 1 or 2
/// (the mate whose boundary was pinned) and `readthrough` is the number of
/// adapter bases that mate carries past the insert, `L - I`.
///
/// Returns `None` when the mates cannot overlap, when no anchor matches, when
/// the only matches are chance hits, or when more than one insert size
/// verifies (ambiguous: refuse rather than guess).
pub fn detect_pe_overlap(
    r1: &[u8],
    r2: &[u8],
    anchor: usize,
) -> Option<(usize /* mate */, usize /* readthrough */)> {
    if anchor < 2 {
        return None;
    }
    let anchor = anchor.max(ANCHOR);
    let l1 = r1.len();
    let l2 = r2.len();
    if l1 < anchor || l2 < anchor {
        return None;
    }
    // The insert cannot exceed the shorter mate: at I = min(l1, l2) the mates
    // touch without sharing anything, so there is no boundary to pin.
    let max_insert = l1.min(l2);
    if max_insert < anchor {
        return None;
    }

    let rc2 = revcomp(r2);
    let l2 = rc2.len();
    // Candidate overlap `o`: R1's prefix of length `o` against the SUFFIX of
    // `revcomp(R2)`, which is the 5' end of R2 read backwards. That is the
    // paper's identity `R1[i] == complement(R2[I-1-i])` for `i < I`: R1 reads
    // the fragment 5'->3' while R2 is the same fragment on the other strand.
    //
    // Scan in DECREASING order (largest insert first, smallest trim first) —
    // the paper's conservative scan. Collect every admissible `o` and refuse
    // unless exactly one exists: the task's ambiguity rule, and the case a
    // repeated fragment produces. A single admissible overlap is the paper's
    // "first admissible" too, since the scan runs largest-first and the window
    // verification is monotone in `o` for a genuine fragment.
    //
    // The anchor is checked first (cheap): the 3' end of the shared window is
    // R1[o-anchor..o] against `rc2[L2-anchor..L2]`. A single 12-base match is
    // ~1e-5, so this is the locality device that keeps the scan from verifying
    // every offset. A chance hit additionally has to agree base for base over
    // the window it implies, which is ~1e-5 * (3/4)^o.
    let mut admissible: Vec<usize> = Vec::new();
    for o in (anchor..=max_insert).rev() {
        if !anchor_matches(&r1[o - anchor..o], &rc2[l2 - anchor..]) {
            continue;
        }
        if window_mismatches(r1, &rc2, o) > MAX_OVERLAP_MISMATCHES {
            continue;
        }
        admissible.push(o);
    }

    // Ambiguity: more than one insert size verifies. Refuse rather than choose.
    admissible.sort_unstable();
    admissible.dedup();
    let o = match admissible.as_slice() {
        [single] => *single,
        _ => return None,
    };
    // A boundary of 0 is no boundary at all: the read is entirely insert.
    if o >= l1 {
        return None;
    }
    Some((1, l1 - o))
}

/// The per-mate adapter trim the detected overlap implies, as
/// `((t5, t3), (t5, t3))` for R1 and R2.
///
/// The boundary falls at the insert length `I` on both mates, and each mate's
/// adapter is read 3′, so mate *i* receives `(0, l_i - I)`. Equal-length mates
/// give the same number twice — the published case, and the only one the
/// original single-value return could express. Unequal mates do not: the
/// single value (derived from R1) over-trimmed R2 by `l1 - l2` bases of real
/// insert, and if that pushed R2 past its budget the mate was dropped as
/// adapter instead of being emitted.
///
/// `adapter_set` is the adapter set the run actually has, or `None` when
/// adapter trimming is disabled — **and `None` refuses the boundary outright**.
/// The safety guard is structural, not a flag default: an overlap proves where
/// the fragment ends, never that anything was ligated past it, so with nothing
/// to check the tail against the mode must not apply a boundary. The check runs
/// on R1, the mate whose boundary the detector pinned.
///
/// Returns `None` when nothing is detected, or when the boundary cannot be
/// checked (no set) or does not hold adapter behind it.
pub fn detect_pe_overlap_trim(
    r1: &[u8],
    r2: &[u8],
    anchor: usize,
    adapter_set: Option<&crate::adapter::AdapterSet>,
) -> Option<((usize, usize), (usize, usize))> {
    detect_pe_overlap_trim_with(
        r1,
        r2,
        anchor,
        adapter_set,
        crate::adapter::AdapterEvidence::Floor,
    )
}

/// [`detect_pe_overlap_trim`] with an explicit evidence policy.
///
/// The policy decides what the guard accepts once the boundary is pinned:
/// [`crate::adapter::AdapterEvidence::Floor`] is today's behaviour,
/// [`crate::adapter::AdapterEvidence::Boundary`] also accepts a tail shorter
/// than the overlap floor, because the mate geometry already fixed the
/// boundary (UT-1).
pub fn detect_pe_overlap_trim_with(
    r1: &[u8],
    r2: &[u8],
    anchor: usize,
    adapter_set: Option<&crate::adapter::AdapterSet>,
    evidence: crate::adapter::AdapterEvidence,
) -> Option<((usize, usize), (usize, usize))> {
    let (_, t3_1) = detect_pe_overlap(r1, r2, anchor)?;
    let set = adapter_set?;
    // R1's readthrough is `l1 - I`, so the insert is `l1 - t3_1`.
    let insert = r1.len().checked_sub(t3_1)?;
    let t3_2 = r2.len().saturating_sub(insert);
    if !readthrough_is_adapter(r1, t3_1, set, evidence) {
        // The mates overlap, but the sequence past the boundary is not
        // adapter: there is nothing to trim. Falling through to the
        // ordinary path (which will also find no adapter) is the safe
        // action, and it is what keeps the mode inert on adapter-free
        // libraries -- the Cell A guard.
        return None;
    }
    Some(((0, t3_1), (0, t3_2)))
}

/// Is the `readthrough` bases at the 3′ end of `seq` an adapter that *starts at*
/// the overlap boundary?
///
/// The overlap gives the boundary; this asks whether the adapter begins exactly
/// there. That distinction is the whole guard: a tail that merely *contains* an
/// adapter fragment further in is ordinary sequence as far as this mate is
/// concerned, and applying the PE boundary to it trims biology. (A real tail can
/// legitimately contain a partial adapter internally -- a 38-base random tail hits
/// a 14-base adapter prefix by chance -- which is why "contains" is not enough.)
///
/// Requires the adapter's own start to coincide with the boundary: the matched
/// position within the tail must be 0.
///
/// # The evidence policy (`evidence`)
///
/// [`AdapterEvidence::Floor`] is today's behaviour: the general matcher,
/// with its `min_overlap_se` floor. A tail shorter than the floor is
/// provably unmatchable there.
///
/// [`AdapterEvidence::Boundary`] lifts the length floor, because the caller
/// has already fixed the boundary from the mate geometry — the tail only
/// has to BE the adapter at that boundary. Without this the detector pins a
/// boundary exactly and then throws it away (UT-1).
///
/// The two rules that protect biology are unchanged in both policies: the
/// adapter must start AT the boundary, and the tail must actually be the
/// adapter prefix.
pub fn readthrough_is_adapter(
    seq: &[u8],
    readthrough: usize,
    set: &crate::adapter::AdapterSet,
    evidence: crate::adapter::AdapterEvidence,
) -> bool {
    if readthrough == 0 || readthrough > seq.len() {
        return false;
    }
    let tail = &seq[seq.len() - readthrough..];
    // PRIMARY (boundary evidence): a boundary-anchored prefix check. The
    // caller already pinned this boundary from the mate geometry, so the
    // length floor is the wrong guard here.
    if evidence == crate::adapter::AdapterEvidence::Boundary
        && crate::adapter::any_prefix_matches(tail, set)
    {
        return true;
    }
    // FALLBACK: the general matcher, unchanged. It also covers tails longer
    // than any adapter, where the prefix check cannot apply.
    match crate::adapter::match_adapters(tail, set) {
        Some(m) => m.t3 == readthrough,
        None => false,
    }
}

// ---------------------------------------------------------------------------
// Tests. The integration contract (including the byte-identical no-regression
// property) lives in tests/pe_overlap.rs.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{Adapter, AdapterEvidence, AdapterSet};

    fn dna(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        (0..len)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((x >> 33) & 3) as usize]
            })
            .collect()
    }

    /// A pair whose insert is `insert` bases on `read_len`-base reads, both
    /// mates carrying adapter readthrough past the fragment.
    fn pair(read_len: usize, insert: usize, seed: u64) -> (Vec<u8>, Vec<u8>) {
        let frag = dna(insert, seed);
        // R1 reads the fragment 5'->3'; R2 reads the same fragment from the
        // opposite strand, so R2 = revcomp(fragment) + its own adapter. Both
        // mates then agree on their first `insert` bases once R2 is put in the
        // fragment's own frame (the paper's R1[i] == complement(R2[I-1-i])).
        let mut r1 = frag.clone();
        r1.extend_from_slice(&dna(read_len - insert, seed ^ 0x9E37_79B9));
        let mut r2: Vec<u8> = frag.iter().rev().map(|&b| comp(b)).collect();
        r2.extend_from_slice(&dna(read_len - insert, seed ^ 0x85EB_CA6B));
        (r1, r2)
    }

    #[test]
    fn recovers_the_readthrough_on_both_mates() {
        let (r1, r2) = pair(100, 60, 5);
        assert_eq!(detect_pe_overlap(&r1, &r2, ANCHOR), Some((1, 40)));
        assert_eq!(detect_pe_overlap(&r2, &r1, ANCHOR), Some((1, 40)));
        // The boundary is applied only with a set to check it against; with
        // none, the structural guard refuses (no set ⇒ no boundary).
        let set = AdapterSet::new(
            vec![Adapter {
                name: "test".into(),
                seq: b"AGATCGGAAGAGC".to_vec(),
            }],
            0.1,
            4,
        );
        assert_eq!(
            detect_pe_overlap_trim(&r1, &r2, ANCHOR, Some(&set)),
            None,
            "this pair's tail is not that adapter: refuse"
        );
    }

    #[test]
    fn unequal_mates_get_their_own_readthrough() {
        // Insert 60, R1 150, R2 100: R1 carries 90 adapter bases past the
        // insert, R2 carries 40. The boundary is the same fragment position on
        // both mates, but the tail length is per-mate. Handing R1's 90 to R2
        // over-trimmed 50 bases of R2's biology (and dropped the mate outright
        // when that exceeded its budget).
        let frag = dna(60, 17);
        let adapter = dna(120, 99);
        let mut r1 = frag.clone();
        r1.extend_from_slice(&adapter[..90]);
        let mut r2: Vec<u8> = frag.iter().rev().map(|&b| comp(b)).collect();
        r2.extend_from_slice(&adapter[..40]);
        assert_eq!((r1.len(), r2.len()), (150, 100));

        let set = AdapterSet::new(
            vec![Adapter {
                name: "long".into(),
                seq: adapter.clone(),
            }],
            0.2,
            4,
        );
        assert_eq!(
            detect_pe_overlap_trim(&r1, &r2, ANCHOR, Some(&set)),
            Some(((0, 90), (0, 40))),
            "each mate's tail is l_i - I: 150-60 and 100-60"
        );
    }

    #[test]
    fn a_boundary_without_adapter_behind_it_is_not_applied() {
        // The mates overlap (insert 60 of 100), but the sequence past the
        // boundary is biology, not adapter. With an adapter set configured the
        // boundary must NOT be applied -- this is the Cell A guard, and it is
        // what keeps the mode from trimming clean reads.
        let (r1, r2) = pair(100, 60, 5);
        let set = AdapterSet::new(
            vec![Adapter {
                name: "test".into(),
                seq: b"AGATCGGAAGAGC".to_vec(),
            }],
            0.1,
            4,
        );
        assert!(
            detect_pe_overlap_trim(&r1, &r2, ANCHOR, Some(&set)).is_none(),
            "an overlap with no adapter behind it must not produce a trim"
        );
        // And with no adapter set at all the boundary is refused too: the
        // guard is structural, not a flag default, so the mode can never
        // apply a boundary it has no set to check.
        assert_eq!(
            detect_pe_overlap_trim(&r1, &r2, ANCHOR, None),
            None,
            "no adapter set must never yield a PE boundary"
        );
    }

    #[test]
    fn too_long_to_overlap_is_none() {
        let (r1, r2) = pair(100, 100, 7);
        assert_eq!(detect_pe_overlap(&r1, &r2, ANCHOR), None);
        assert_eq!(detect_pe_overlap(&dna(60, 11), &dna(60, 22), ANCHOR), None);
    }

    #[test]
    fn a_boundary_of_zero_is_not_a_boundary() {
        // insert == read_len: the mates share the whole read.
        let (r1, r2) = pair(100, 100, 31);
        assert_eq!(detect_pe_overlap(&r1, &r2, ANCHOR), None);
    }

    /// A repeated fragment admits several insert sizes at once: the pair is one
    /// fragment, but the geometry cannot say which of the candidate overlaps is
    /// the true boundary. Refuse rather than guess.
    #[test]
    fn ambiguity_is_refused() {
        // "ACGT" repeated is periodic with period 4, so every insert size that
        // is a multiple of 4 away from the true one verifies as well.
        let unit = b"ACGT";
        let frag: Vec<u8> = (0..24).map(|i| unit[i % 4]).collect();
        let mut r1 = frag.clone();
        r1.extend_from_slice(b"TTTTTTTTTTTT");
        let mut r2: Vec<u8> = frag.iter().rev().map(|&b| comp(b)).collect();
        r2.extend_from_slice(b"GGGGGGGGGGGG");
        assert_eq!(detect_pe_overlap(&r1, &r2, ANCHOR), None,
            "a repeated fragment admits several overlaps: refuse");
    }

    #[test]
    fn a_chance_anchor_hit_is_refused() {
        // The discipline of fastcheck's
        // `a_chance_anchor_hit_does_not_produce_a_measurement`: the anchor
        // matches a window of R1 by construction, but the rest of the
        // alignment is not the complement of R1, so the full verification
        // rejects it.
        let s1 = dna(60, 17);
        let p = 40;
        let mut frame = s1[p..p + ANCHOR].to_vec();
        for &b in &s1[p + ANCHOR..] {
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
        assert_eq!(detect_pe_overlap(&s1, &r2, ANCHOR), None);
    }

    #[test]
    fn one_error_inside_the_overlap_still_detects() {
        let (mut r1, r2) = pair(100, 60, 33);
        let i = 30;
        r1[i] = if r1[i] == b'A' { b'C' } else { b'A' };
        assert_eq!(detect_pe_overlap(&r1, &r2, ANCHOR), Some((1, 40)));
    }

    #[test]
    fn revcomp_is_an_involution() {
        let s = dna(37, 29);
        assert_eq!(revcomp(&revcomp(&s)), s);
    }

    #[test]
    fn n_in_either_mate_never_fabricates_a_hit() {
        let n = vec![b'N'; 60];
        let r = dna(60, 61);
        assert_eq!(detect_pe_overlap(&n, &r, ANCHOR), None);
        assert_eq!(detect_pe_overlap(&r, &n, ANCHOR), None);
    }

    /// A corrupted anchor region is refused: the anchor is part of the shared
    /// window, so breaking it breaks the candidate.
    #[test]
    fn a_corrupted_anchor_is_refused() {
        let (mut r1, mut r2) = pair(100, 60, 13);
        // The anchor pairs R1[48..60] with R2[0..12]; break both.
        r1[48..60].copy_from_slice(b"ACGTACGTACGT");
        r2[..12].copy_from_slice(b"TTTTTTTTTTTT");
        assert_eq!(detect_pe_overlap(&r1, &r2, ANCHOR), None);
    }

    // ------------------------------------------------------------------
    // UT-1: the pinned boundary must survive its own guard.
    //
    // `readthrough_is_adapter` re-checks the tail with the general matcher,
    // which applies `min_overlap_se` to the tail it is given. A 10-base
    // tail is therefore PROVABLY unmatchable there, and the boundary the
    // detector pinned exactly is thrown away.
    //
    // The fix is a boundary-anchored prefix check that does not consult the
    // floor, used only where the mate geometry already fixed the boundary.
    // ------------------------------------------------------------------

    const TRUSEQ: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

    fn illumina_set(min_overlap: usize) -> AdapterSet {
        AdapterSet::new(
            vec![Adapter {
                name: "TruSeq/R1".into(),
                seq: TRUSEQ.to_vec(),
            }],
            0.2,
            min_overlap,
        )
    }

    /// A 10-base adapter readthrough: below the default floor of 14.
    fn short_readthrough_seq() -> Vec<u8> {
        let mut seq = dna(40, 4242);
        seq.extend_from_slice(&TRUSEQ[..10]);
        seq
    }

    #[test]
    fn a_pinned_boundary_survives_a_tail_below_the_overlap_floor() {
        // The floor the guard must NOT re-impose on a pinned boundary.
        let set = illumina_set(14);
        let seq = short_readthrough_seq();
        assert!(
            readthrough_is_adapter(&seq, 10, &set, AdapterEvidence::Boundary),
            "a boundary the detector pinned must survive the tail check"
        );
    }

    #[test]
    fn the_floor_policy_still_refuses_a_short_tail() {
        // The byte-identical property: today's behaviour is untouched.
        let set = illumina_set(14);
        let seq = short_readthrough_seq();
        assert!(
            !readthrough_is_adapter(&seq, 10, &set, AdapterEvidence::Floor),
            "without independent evidence the floor still refuses a short tail"
        );
    }

    #[test]
    fn a_tail_that_is_not_the_adapter_prefix_is_refused() {
        // Boundary evidence does not make a random tail acceptable: the tail
        // must still BE adapter at the boundary.
        let set = illumina_set(14);
        let mut seq = dna(40, 4242);
        seq.extend_from_slice(b"TTTTTTTTTT");
        assert!(
            !readthrough_is_adapter(&seq, 10, &set, AdapterEvidence::Boundary),
            "a tail that is not the adapter prefix must be refused"
        );
    }

    #[test]
    fn the_floor_policy_still_accepts_a_long_tail() {
        // The long-tail path is unchanged in both policies.
        let set = illumina_set(14);
        let mut seq = dna(40, 4242);
        seq.extend_from_slice(&TRUSEQ[..20]);
        for ev in [AdapterEvidence::Floor, AdapterEvidence::Boundary] {
            assert!(
                readthrough_is_adapter(&seq, 20, &set, ev),
                "a 20-base adapter tail passes the floor in both policies"
            );
        }
    }
}

