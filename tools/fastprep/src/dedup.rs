//! Deduplication: the first stage of the pipeline, ahead of every other
//! operation on the read. See the specification, §6.1.
//!
//! Each read (or pair) is reduced to a fingerprint over an *anchor window*
//! chosen to lie outside the region the configured pipeline will trim. That
//! matters on two-colour chemistry: polyG tails come from per-cluster signal
//! dropout, which is stochastic, so two amplification copies of one fragment
//! can read `INSERT+GGGGGGGG` and `INSERT+GGGGGGGGGGGG`. Fingerprinting the
//! full read before trimming would keep both.
//!
//! The fingerprint set is the only unbounded state in the pipeline, so it is
//! admission-checked against [`crate::memory`] before allocation and
//! *saturates* rather than being OOM-killed if the projection was wrong.

use rustc_hash::FxHashSet;
use std::hash::Hash;
use xxhash_rust::xxh3::{xxh3_128, xxh3_64};

use crate::algorithm::Params;
use crate::memory;

// ---------------------------------------------------------------------------
// Fingerprints.
// ---------------------------------------------------------------------------

/// A fingerprint width. Both widths use xxh3.
///
/// `hash_pair` is order-dependent by construction: R1 and R2 carry distinct
/// geometric meaning, so swapping them must not collide.
pub trait Fingerprint: Hash + Eq + Copy + Send + Sync {
    const BYTES: usize;
    fn of(seq: &[u8]) -> Self;
    fn of_pair(r1: &[u8], r2: &[u8]) -> Self;
}

impl Fingerprint for u64 {
    const BYTES: usize = 8;

    #[inline(always)]
    fn of(seq: &[u8]) -> Self {
        xxh3_64(seq)
    }

    /// A plain XOR would vanish identically whenever `r1 == r2`, mapping every
    /// pair with identical mates onto zero. Rotating one operand by half the
    /// width removes that: rotation is bijective, so it preserves entropy,
    /// and `h ^ rot(h, 32) == 0` only for hashes with two identical halves.
    #[inline(always)]
    fn of_pair(r1: &[u8], r2: &[u8]) -> Self {
        xxh3_64(r1) ^ xxh3_64(r2).rotate_left(32)
    }
}

impl Fingerprint for u128 {
    const BYTES: usize = 16;

    #[inline(always)]
    fn of(seq: &[u8]) -> Self {
        xxh3_128(seq)
    }

    #[inline(always)]
    fn of_pair(r1: &[u8], r2: &[u8]) -> Self {
        xxh3_128(r1) ^ xxh3_128(r2).rotate_left(64)
    }
}

/// Probability that *any* collision exists across `x` fingerprints of `bits`
/// width, by the birthday square approximation `x^2 / 2^(bits+1)`.
///
/// This bounds the existence of one collision over the whole run, not a
/// per-read false-duplicate rate: a collision removes exactly one read.
pub fn collision_probability(x: u64, bits: u32) -> f64 {
    let x = x as f64;
    x * x / 2f64.powi(bits as i32 + 1)
}

/// Narrowest width whose collision probability stays under `threshold`.
pub fn select_width(estimated_reads: u64, threshold: f64) -> u32 {
    if collision_probability(estimated_reads, 64) < threshold { 64 } else { 128 }
}

// ---------------------------------------------------------------------------
// Anchor window.
// ---------------------------------------------------------------------------

/// Fingerprinting fewer than this many bases is not worth doing; the caller
/// falls back to the whole available read.
const MIN_ANCHOR: usize = 20;

/// Which end of the read the fingerprint is taken from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorEnd {
    /// 5′ prefix. Correct when the 5′ end is not scanned.
    Five,
    /// 3′ suffix. Used when the 5′ end is scanned instead.
    Three,
    /// Both ends are targeted, so no artefact-free window exists.
    None,
}

/// The anchor window, resolved **once** for the run.
///
/// The width must not depend on the individual read's length. Deriving it
/// per read (`k = min(L_min, l - ...)` with a fractional `L_min`) defeats the
/// purpose: two copies of one fragment that differ in length get different
/// widths, hence different windows, hence different fingerprints, which is
/// exactly the divergence the anchor exists to absorb. Illumina output is
/// fixed-length so the two agree there, but already-trimmed input, ONT and
/// FASTA are not, and the property should not be conditional on the platform.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnchorSpec {
    pub end: AnchorEnd,
    /// Fixed window width in bases.
    pub len: usize,
    t5: usize,
    t3: usize,
}

impl AnchorSpec {
    /// `reference_len` is the expected read length (the `--read-length` hint,
    /// or the length of the first record). It fixes the window width for the
    /// whole run.
    pub fn new(params: &Params, reference_len: usize, adapters_enabled: bool) -> Self {
        let end = anchor_end(params, adapters_enabled);
        let t5 = params.trim_5.unwrap_or(0).min(reference_len);
        let t3 = params.trim_3.unwrap_or(0).min(reference_len - t5);
        let usable = reference_len - t5 - t3;
        let len = params.min_len.resolve(reference_len).min(usable);
        Self { end, len, t5, t3 }
    }

    /// Whether a usable artefact-free window exists at all.
    pub fn is_usable(&self) -> bool {
        self.end != AnchorEnd::None && self.len >= MIN_ANCHOR
    }

    /// The slice to fingerprint.
    ///
    /// Falls back to the whole read when no window is usable, or when this
    /// read is too short to carry the fixed window. Reads shorter than the
    /// reference therefore lose the invariance; that is unavoidable and
    /// affects only reads already below `L_min`.
    #[inline]
    pub fn slice<'a>(&self, seq: &'a [u8]) -> &'a [u8] {
        if !self.is_usable() {
            return seq;
        }
        let l = seq.len();
        if l < self.t5 + self.t3 + self.len {
            return seq;
        }
        match self.end {
            AnchorEnd::Five => &seq[self.t5..self.t5 + self.len],
            AnchorEnd::Three => &seq[l - self.t3 - self.len..l - self.t3],
            AnchorEnd::None => seq,
        }
    }
}

/// Picks the anchor end: the end the configured pipeline will *not* touch.
///
/// Adapter readthrough is a 3′ event, so an enabled adapter set makes the 3′
/// end unsafe just as `scan_3` does. Hard trims are deterministic and do not
/// cause copies to diverge, so they only shift the window's origin.
pub fn anchor_end(params: &Params, adapters_enabled: bool) -> AnchorEnd {
    let five_unsafe = params.scan_5;
    let three_unsafe = params.scan_3 || adapters_enabled;
    match (five_unsafe, three_unsafe) {
        (false, _) => AnchorEnd::Five,
        (true, false) => AnchorEnd::Three,
        (true, true) => AnchorEnd::None,
    }
}

// ---------------------------------------------------------------------------
// Verifier.
// ---------------------------------------------------------------------------

/// Outcome of offering a fingerprint to the set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Seen {
    /// Not seen before, and recorded. Emit the read.
    New,
    /// Already present. Drop the read as a duplicate.
    Duplicate,
    /// Not seen before, but the set is saturated so it was *not* recorded.
    /// Emit the read; later copies of it will not be detected.
    NewUnrecorded,
}

/// How many insertions between memory polls. Polling reads a file, so it is
/// far too expensive per read, and far too coarse per chunk.
const POLL_INTERVAL: u64 = 1 << 16;

/// The fingerprint set, with the saturation behaviour of §6.1.2.
///
/// Once resident usage crosses `ceiling`, the set stops accepting insertions
/// but keeps answering queries. This is one-sided: the saturated set is a
/// subset of the unsaturated one at every point, so a saturated run can fail
/// to remove duplicates but can never remove a read that the unsaturated
/// procedure would have kept.
pub struct Verifier<T: Fingerprint> {
    seen: FxHashSet<T>,
    ceiling: Option<u64>,
    since_poll: u64,
    saturated_at: Option<u64>,
    offered: u64,
}

impl<T: Fingerprint> Verifier<T> {
    /// `ceiling` is the resident-usage bound above which the set saturates;
    /// `None` disables the runtime check.
    pub fn with_capacity(capacity: usize, ceiling: Option<u64>) -> Self {
        Self {
            seen: FxHashSet::with_capacity_and_hasher(capacity, Default::default()),
            ceiling,
            since_poll: 0,
            saturated_at: None,
            offered: 0,
        }
    }

    /// Offers a fingerprint. `&mut self`, so this runs on the ordered thread;
    /// computing the fingerprint itself is a pure function and belongs in the
    /// parallel stage.
    pub fn offer(&mut self, fp: T) -> Seen {
        self.offered += 1;
        if self.seen.contains(&fp) {
            return Seen::Duplicate;
        }
        if self.saturated() {
            return Seen::NewUnrecorded;
        }
        self.since_poll += 1;
        if self.since_poll >= POLL_INTERVAL {
            self.since_poll = 0;
            if self.over_ceiling() {
                self.saturated_at = Some(self.offered);
                return Seen::NewUnrecorded;
            }
        }
        self.seen.insert(fp);
        Seen::New
    }

    fn over_ceiling(&self) -> bool {
        match (self.ceiling, memory::current_usage()) {
            (Some(c), Some(used)) => used >= c,
            _ => false,
        }
    }

    #[inline]
    pub fn saturated(&self) -> bool {
        self.saturated_at.is_some()
    }

    /// Read index at which saturation began, for the statistics block.
    pub fn saturated_at(&self) -> Option<u64> {
        self.saturated_at
    }

    pub fn len(&self) -> usize {
        self.seen.len()
    }
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::{MatchRule, MinLen};

    fn params(min_len: MinLen, scan_5: bool, scan_3: bool) -> Params {
        let mut p = Params { min_len, scan_5, scan_3, ..Default::default() };
        p.prepare();
        let _ = MatchRule::Homopolymer;
        p
    }

    // --- fingerprints ---

    #[test]
    fn fingerprint_is_deterministic() {
        assert_eq!(u64::of(b"ACGTACGT"), u64::of(b"ACGTACGT"));
        assert_eq!(u128::of(b"ACGTACGT"), u128::of(b"ACGTACGT"));
    }

    #[test]
    fn pair_fingerprint_is_order_dependent() {
        let (a, b) = (b"ACGTACGT".as_slice(), b"TGCATGCA".as_slice());
        assert_ne!(u64::of_pair(a, b), u64::of_pair(b, a));
        assert_ne!(u128::of_pair(a, b), u128::of_pair(b, a));
    }

    /// The failure a plain XOR would produce: identical mates mapping to zero.
    #[test]
    fn identical_mates_do_not_cancel() {
        let s = b"ACGTACGTACGTACGT".as_slice();
        assert_ne!(u64::of_pair(s, s), 0);
        assert_ne!(u128::of_pair(s, s), 0);
    }

    #[test]
    fn width_selection_follows_the_birthday_bound() {
        assert_eq!(select_width(1_000_000, 0.001), 64);
        assert_eq!(select_width(2_000_000_000, 0.001), 128);
    }

    #[test]
    fn collision_probability_at_the_documented_point() {
        // ~0.6e9 sequences at 64 bits is where p reaches ~1%.
        let p = collision_probability(600_000_000, 64);
        assert!((0.009..0.011).contains(&p), "p = {p}");
    }

    // --- anchor window ---

    #[test]
    fn anchors_on_the_end_that_is_not_scanned() {
        assert_eq!(anchor_end(&params(MinLen::Frac(0.6), false, true), true), AnchorEnd::Five);
        assert_eq!(anchor_end(&params(MinLen::Frac(0.6), true, false), false), AnchorEnd::Three);
        assert_eq!(anchor_end(&params(MinLen::Frac(0.6), true, true), false), AnchorEnd::None);
        // Adapters make the 3' end unsafe even with scan_3 off.
        assert_eq!(anchor_end(&params(MinLen::Frac(0.6), true, false), true), AnchorEnd::None);
    }

    #[test]
    fn window_is_lmin_wide_from_the_chosen_end() {
        let p = params(MinLen::Frac(0.6), false, true);
        let five = AnchorSpec::new(&p, 150, true);
        assert_eq!((five.end, five.len), (AnchorEnd::Five, 90));
        assert_eq!(five.slice(&vec![b'A'; 150]).len(), 90);

        let q = params(MinLen::Frac(0.6), true, false);
        let three = AnchorSpec::new(&q, 150, false);
        assert_eq!((three.end, three.len), (AnchorEnd::Three, 90));
    }

    #[test]
    fn window_starts_after_a_hard_trim() {
        let mut p = params(MinLen::Abs(50), false, true);
        p.trim_5 = Some(10);
        p.trim_3 = Some(5);
        let spec = AnchorSpec::new(&p, 150, true);
        assert_eq!(spec.len, 50);
        let seq: Vec<u8> = (0..150u8).collect();
        assert_eq!(spec.slice(&seq), &seq[10..60]);
    }

    #[test]
    fn short_reads_fall_back_to_the_whole_read() {
        let p = params(MinLen::Abs(10), false, true);
        let spec = AnchorSpec::new(&p, 30, true);
        assert!(!spec.is_usable(), "10-base window is below MIN_ANCHOR");
        assert_eq!(spec.slice(b"ACGTACGTAC"), b"ACGTACGTAC");
    }

    /// The property the anchor exists for. This failed against a per-read
    /// window width, which is why `AnchorSpec` fixes the width for the run.
    #[test]
    fn anchor_is_invariant_to_divergent_polyg_tails() {
        let p = params(MinLen::Frac(0.6), false, true);
        let spec = AnchorSpec::new(&p, 150, true);

        let insert: Vec<u8> = b"ACGTTGCA".iter().cycle().take(120).copied().collect();
        let mut short_tail = insert.clone();
        short_tail.extend(std::iter::repeat_n(b'G', 8));
        let mut long_tail = insert.clone();
        long_tail.extend(std::iter::repeat_n(b'G', 30));

        assert_ne!(short_tail.len(), long_tail.len(), "copies differ in length");
        assert_eq!(
            u64::of(spec.slice(&short_tail)),
            u64::of(spec.slice(&long_tail)),
            "copies with divergent polyG tails must fingerprint identically"
        );
        // And the full-read fingerprint, which is what a naive implementation
        // would use, does not have the property.
        assert_ne!(u64::of(&short_tail), u64::of(&long_tail));
    }

    #[test]
    fn anchor_still_separates_genuinely_different_fragments() {
        let p = params(MinLen::Frac(0.6), false, true);
        let spec = AnchorSpec::new(&p, 150, true);
        let a: Vec<u8> = b"ACGTTGCA".iter().cycle().take(150).copied().collect();
        let mut b = a.clone();
        // Pick a base that differs from whatever is there. Hardcoding the
        // replacement made this a no-op: position 3 of "ACGTTGCA" is already
        // a T, so the test compared a read with itself and the property went
        // untested until the assertion below was added.
        b[3] = if a[3] == b'A' { b'C' } else { b'A' };
        assert_ne!(a, b, "the mutation must actually change the read");
        assert_ne!(u64::of(spec.slice(&a)), u64::of(spec.slice(&b)));
    }

    // --- verifier ---

    #[test]
    fn detects_duplicates() {
        let mut v = Verifier::<u64>::with_capacity(16, None);
        assert_eq!(v.offer(u64::of(b"ACGT")), Seen::New);
        assert_eq!(v.offer(u64::of(b"ACGT")), Seen::Duplicate);
        assert_eq!(v.offer(u64::of(b"TGCA")), Seen::New);
        assert_eq!(v.len(), 2);
    }

    #[test]
    fn saturated_verifier_still_detects_known_duplicates() {
        let mut v = Verifier::<u64>::with_capacity(16, None);
        assert_eq!(v.offer(u64::of(b"ACGT")), Seen::New);
        v.saturated_at = Some(1);
        // Known fingerprint: still caught.
        assert_eq!(v.offer(u64::of(b"ACGT")), Seen::Duplicate);
        // New fingerprint: passed through, and not recorded.
        assert_eq!(v.offer(u64::of(b"TGCA")), Seen::NewUnrecorded);
        assert_eq!(v.offer(u64::of(b"TGCA")), Seen::NewUnrecorded);
        assert_eq!(v.len(), 1);
    }

    /// The one-sidedness of §6.1.2, tested directly rather than assumed: the
    /// saturated run's survivors are a superset of the unsaturated run's.
    #[test]
    fn saturation_is_one_sided() {
        let reads: Vec<Vec<u8>> = (0..200u32)
            .map(|i| format!("READ{}", i % 70).into_bytes())
            .collect();

        let mut full = Verifier::<u64>::with_capacity(256, None);
        let unsaturated: Vec<bool> =
            reads.iter().map(|r| full.offer(u64::of(r)) != Seen::Duplicate).collect();

        for freeze_at in [1usize, 5, 37, 199] {
            let mut v = Verifier::<u64>::with_capacity(256, None);
            let saturated: Vec<bool> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    if i == freeze_at {
                        v.saturated_at = Some(i as u64);
                    }
                    v.offer(u64::of(r)) != Seen::Duplicate
                })
                .collect();

            for (i, (&sat, &unsat)) in saturated.iter().zip(&unsaturated).enumerate() {
                assert!(
                    !unsat || sat,
                    "freeze_at={freeze_at}: read {i} kept unsaturated but dropped saturated"
                );
            }
        }
    }
}
