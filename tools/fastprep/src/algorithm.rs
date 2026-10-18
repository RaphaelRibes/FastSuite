//! Core polyX trimming algorithm.
//!
//! Implements the specification faithfully: asymmetric quality-weighted
//! per-step score, separated end-scan and interior-scan, fused forward
//! pass + budget-capped reverse pass + reconciliation.
//!
//! This module is I/O-free: `polyx_trim` takes bases and qualities as
//! slices and returns a decision. The caller wraps it in needletail +
//! flate2 and a rayon pool.

use std::str::FromStr;
use std::sync::Arc;

use crate::adapter::AdapterSet;

// ---------------------------------------------------------------------------
// MinLen — either an absolute base count or a fraction of read length L.
// ---------------------------------------------------------------------------

/// Minimum post-trim length, either absolute or as a fraction of `L`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MinLen {
    Abs(usize),
    Frac(f64),
}

impl MinLen {
    /// Resolve against a concrete read length `l`.
    pub fn resolve(self, l: usize) -> usize {
        match self {
            MinLen::Abs(n) => n.min(l),
            MinLen::Frac(f) => ((f * l as f64).floor() as usize).min(l),
        }
    }
}

impl FromStr for MinLen {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(num) = s.strip_suffix('L').or_else(|| s.strip_suffix('l')) {
            let f: f64 = num
                .parse()
                .map_err(|_| format!("cannot parse `{}` as fraction", num))?;
            if !(0.0 < f && f <= 1.0) {
                return Err(format!("fraction must be in (0, 1], got {}", f));
            }
            Ok(MinLen::Frac(f))
        } else {
            let n: usize = s
                .parse()
                .map_err(|_| format!("cannot parse `{}` as integer", s))?;
            Ok(MinLen::Abs(n))
        }
    }
}

// ---------------------------------------------------------------------------
// MatchRule — homopolymer or anchored to a target set.
// ---------------------------------------------------------------------------

/// What counts as a match at a given step.
///
/// `Homopolymer`: match iff `b[k] == b[k+1]` (the generic polyX case).
/// `Anchored(mask)`: match iff `b[k+1]` is in the target set (per spec
/// remark §1.4; note this does not check `b[k]` at `k = 0`).
#[derive(Debug, Clone)]
pub enum MatchRule {
    Homopolymer,
    Anchored(Box<[bool; 256]>),
}

impl MatchRule {
    /// Construct an anchored rule from a set of ASCII letters.
    /// Accepts both upper and lowercase for each letter.
    pub fn anchored(letters: &[u8]) -> Self {
        let mut mask = Box::new([false; 256]);
        for &b in letters {
            mask[b.to_ascii_uppercase() as usize] = true;
            mask[b.to_ascii_lowercase() as usize] = true;
        }
        MatchRule::Anchored(mask)
    }
}

// ---------------------------------------------------------------------------
// Direction and stop modes.
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Reverse,
}

/// Why a read was dropped. Stored in [`TrimResult::Drop`] so callers can
/// attribute drops to specific causes without a separate channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    /// Hard-trim (--trim-5 / --trim-3) alone exhausted the budget or the read.
    Hard,
    /// 5′-end polyX scan exceeded the trim budget.
    Poly5,
    /// 3′-end polyX scan exceeded the trim budget.
    Poly3,
    /// Interior polyX accumulation pushed the total over the trim budget.
    PolyInt,
    /// Whole-read mean quality fell below --mean-q.
    MeanQ,
    /// Adapter trim alone exhausted the min-len budget.
    Adapter,
    /// Removed by deduplication, before any trimming ran. Carries no trim
    /// counts: the read was never scanned.
    Duplicate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopMode {
    Mercy,
    Boundary,
    Budget,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanResult {
    pub k_star: usize,
    pub mode: StopMode,
    /// Number of matching steps in the scan (excludes mismatches).
    /// Used to threshold against `m_h`: a run qualifies iff
    /// `matches + 1 ≥ m_h` (the +1 counts the initial anchor base).
    pub matches: usize,
    /// Boundary-anchored run extent in bases: `last_match_step + 2`, i.e.
    /// the number of bases from `start` up to and including the far base of
    /// the last matching step. This is the homopolymer run's true span — it
    /// keeps error bases inside the run but excludes the flanking overshoot
    /// the mercy walk absorbs. 0 when the scan saw no matching step.
    pub run_span: usize,
}

// ---------------------------------------------------------------------------
// Params — all knobs from the spec.
// ---------------------------------------------------------------------------

/// All parameters configuring the pipeline for one read.
///
/// Defaults match table (2) of the specification. The `direction` field
/// is applied uniformly to all scans: `Forward` corresponds to 5′→3′
/// primary end-scan; `Reverse` flips the 5′ and 3′ roles (useful for
/// Nextera R2 polyG).
#[derive(Debug, Clone)]
pub struct Params {
    pub min_len: MinLen,
    pub mercy: f64,             // Δ_m
    pub phred_sat: u8,          // t
    pub gate_shape: f64,        // x
    pub match_leniency: f64,    // μ
    pub mismatch_strength: f64, // ν
    pub min_run: usize,         // m_h
    /// NovaSeq two-colour first-cycle G-dropout: a short (1–4 bp) leading
    /// G-run at the 5′ end with *normal* quality, so the quality-gated
    /// scan never sees it (it looks for the quality collapse of a long
    /// run, and `min_run` rejects short ones). `Some(n)` trims exactly
    /// the leading G-run when its length is in `[1, n]` — a hard trim,
    /// no mercy walk — before the scan runs. `n < min_run` keeps the
    /// scan in charge of longer runs. `None` disables the rule.
    pub head_dropout: Option<usize>,
    pub rule: MatchRule,
    pub scan_5: bool,
    pub scan_3: bool,
    pub scan_interior: bool,
    /// Boundary-anchored cut (new algorithm, opt-in). When `true`, a polyX
    /// scan that terminates by mercy trims only the homopolymer *run extent*
    /// (start ..= last matching base) instead of `k_star`, the full mercy
    /// walk. The mercy walk absorbs a few flanking non-run bases past the run
    /// before the gap mass trips Δ_m — on a clean run followed by ordinary
    /// biology that overshoots the boundary by ~2 bp (measured). Cutting the
    /// run extent places the cut exactly at the run/flanking boundary, keeping
    /// error bases *inside* the run (the last match sits past them) while
    /// dropping the flanking overshoot. Default `false` preserves the shipped
    /// whole-run/mercy cut.
    pub boundary_anchored: bool,
    pub direction: Direction,
    pub mean_q: Option<f64>,
    pub trim_5: Option<usize>,
    pub trim_3: Option<usize>,

    /// Adapter configuration. `None` disables adapter trimming.
    /// Wrapped in `Arc` so per-mate `Params` clones share the allocation.
    pub adapter: Option<Arc<AdapterSet>>,

    /// Trim budget already consumed by a prior step (e.g. adapter trim).
    /// Subtracted from the computed Δ_L inside [`polyx_trim`] so the
    /// combined trim stays within the min-len contract.
    pub extra_committed: usize,

    /// If `Some(l)`, use `l` instead of the actual slice length when
    /// computing Δ_L. Set by [`trim`] before calling [`polyx_trim`] on
    /// an adapter-sliced sub-read so that min-len is enforced relative
    /// to the original read length.
    pub original_len: Option<usize>,

    /// Precomputed quality-gate values: `gate_lut[q] = g(q)`.
    /// Built by [`Params::prepare`] (called automatically by
    /// `Default::default()` and `cli::Cli::finalize()`). On the hot
    /// path, `per_step_gap_and_match` reads from this table instead
    /// of calling `powf`, eliminating the dominant micro-cost.
    #[doc(hidden)]
    pub gate_lut: [f64; 256],
}

impl Default for Params {
    fn default() -> Self {
        let mut p = Self {
            min_len: MinLen::Frac(0.6),
            mercy: 3.0,
            phred_sat: 20,
            gate_shape: 1.0,
            match_leniency: 0.1,
            mismatch_strength: 0.5,
            min_run: 6,
            head_dropout: None,
            rule: MatchRule::Homopolymer,
            scan_5: true,
            scan_3: true,
            scan_interior: true,
            boundary_anchored: false,
            direction: Direction::Forward,
            mean_q: None,
            trim_5: None,
            trim_3: None,
            adapter: None,
            extra_committed: 0,
            original_len: None,
            gate_lut: [0.0; 256],
        };
        p.prepare();
        p
    }
}

impl Params {
    /// Trim budget Δ_L = L_eff − L_min − extra_committed.
    /// `L_eff` is `original_len` when set (adapter-trim sub-read path),
    /// otherwise `l` (the actual slice length).
    pub fn trim_budget(&self, l: usize) -> usize {
        let l_eff = self.original_len.unwrap_or(l);
        l_eff
            .saturating_sub(self.min_len.resolve(l_eff))
            .saturating_sub(self.extra_committed)
    }

    /// Populate `gate_lut` from the current `phred_sat` and `gate_shape`.
    /// Cheap (~256 `powf` calls, ~10 µs) but should not be done per
    /// read.
    ///
    /// `Default::default()` and `cli::Cli::finalize()` call this for
    /// you. Tests that mutate `phred_sat` or `gate_shape` after
    /// construction must call `prepare()` again to keep the LUT in
    /// sync; tests that mutate other fields don't need to.
    pub fn prepare(&mut self) {
        for q in 0u8..=255 {
            self.gate_lut[q as usize] = quality_gate(q, self.phred_sat, self.gate_shape);
        }
    }
}

// ---------------------------------------------------------------------------
// Scan primitive (§1.7 / @alg:scan).
// ---------------------------------------------------------------------------

/// Quality gate g(q) from spec §1.3.
#[inline]
fn quality_gate(q: u8, t: u8, x: f64) -> f64 {
    if q == 0 {
        0.0
    } else if q >= t {
        1.0
    } else {
        (q as f64 / t as f64).powf(x)
    }
}

/// Per-step gap contribution + match flag, for a single pair.
/// Per-step gap contribution + match flag for a single pair.
///
/// Hot-path function: called per (b1, b2, q1, q2) quadruple inside the
/// scan loop. Cost was once dominated by two `powf` calls in
/// `quality_gate`; we now pull `g(q)` from a 256-entry lookup table
/// on `Params`, populated once via `Params::prepare()`. Per-call cost:
/// two LUT loads + one fadd + a branch.
#[inline]
fn per_step_gap_and_match(
    b1: u8,
    b2: u8,
    q1: u8,
    q2: u8,
    params: &Params,
    rule: &MatchRule,
) -> (f64, bool) {
    let matches = match rule {
        MatchRule::Homopolymer => b1 == b2,
        MatchRule::Anchored(mask) => mask[b2 as usize],
    };
    // Quality bytes are raw Phred values: pipeline::ScanInput pre-decodes
    // the FASTQ ASCII offset on ingestion. Tests construct raw values
    // directly. So no offset subtraction here.
    let w = (params.gate_lut[q1 as usize] + params.gate_lut[q2 as usize]) * 0.5;
    let gap = if matches {
        // σ = (1 − μ) + μW  ⇒  gap = μ(1 − W)
        params.match_leniency * (1.0 - w)
    } else {
        // σ = −νW  ⇒  gap = 1 + νW
        1.0 + params.mismatch_strength * w
    };
    (gap, matches)
}

/// Direction-aware index translation.
#[inline]
fn idx(i: usize, n: usize, dir: Direction) -> usize {
    match dir {
        Direction::Forward => i,
        Direction::Reverse => n - 1 - i,
    }
}

/// Walk from offset `start` toward `upper`, accumulating gap mass.
///
/// `upper` is an inclusive upper bound on `start + K` (mirrors spec's
/// boundary clause `s + K = u`). `budget`, if present, caps K at that
/// value and returns `StopMode::Budget` — used by the 5′ short-circuit.
pub fn scan(
    bases: &[u8],
    quals: &[u8],
    start: usize,
    upper: usize,
    direction: Direction,
    params: &Params,
    budget: Option<usize>,
) -> ScanResult {
    debug_assert_eq!(bases.len(), quals.len());
    debug_assert!(upper < bases.len());
    debug_assert!(start <= upper);

    let n = bases.len();
    let mut gamma = 0.0f64;
    let mut k = 0usize;
    let mut matches = 0usize;
    let mut run_span = 0usize;

    loop {
        // Boundary check: next pair would need index start + k + 1 ≤ upper.
        if start + k + 1 > upper {
            return ScanResult {
                k_star: k,
                mode: StopMode::Boundary,
                matches,
                run_span,
            };
        }
        // Budget check (optional).
        if let Some(b) = budget {
            if k >= b {
                return ScanResult {
                    k_star: k,
                    mode: StopMode::Budget,
                    matches,
                    run_span,
                };
            }
        }
        let i1 = idx(start + k, n, direction);
        let i2 = idx(start + k + 1, n, direction);
        let (gap, is_match) = per_step_gap_and_match(
            bases[i1],
            bases[i2],
            quals[i1],
            quals[i2],
            params,
            &params.rule,
        );
        gamma += gap;
        if is_match {
            matches += 1;
            // The matching pair (start+k, start+k+1) confirms the run reaches
            // offset k+1; the run extent from `start` is therefore k+2 bases.
            // Only tracked when the boundary-anchored cut will consume it, so
            // the default whole-run path does no extra work.
            if params.boundary_anchored {
                run_span = k + 2;
            }
        }
        k += 1;
        if gamma > params.mercy {
            return ScanResult {
                k_star: k,
                mode: StopMode::Mercy,
                matches,
                run_span,
            };
        }
    }
}

// ---------------------------------------------------------------------------
// Main pipeline (@alg:main).
// ---------------------------------------------------------------------------

/// Decision for one read.
///
/// Both variants carry separate adapter and polyX trim amounts so the
/// pipeline can attribute bases precisely. `Emit.runs` lists surviving
/// interior-run positions as `(start, len)` pairs in absolute 0-based
/// coordinates of the original read (after adapter trim is added back).
/// `h_int` is derivable as `runs.iter().map(|(_, l)| l).sum()`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrimResult {
    Drop {
        reason: DropReason,
        t5_adapter: usize,
        t3_adapter: usize,
        t5_poly: usize,
        t3_poly: usize,
        runs: Vec<(usize, usize)>,
    },
    Emit {
        t5_adapter: usize,
        t3_adapter: usize,
        t5_poly: usize,
        t3_poly: usize,
        runs: Vec<(usize, usize)>,
    },
}

impl TrimResult {
    /// Total 5′ trim (adapter + polyX).
    pub fn t5(&self) -> usize {
        match self {
            TrimResult::Drop {
                t5_adapter,
                t5_poly,
                ..
            } => t5_adapter + t5_poly,
            TrimResult::Emit {
                t5_adapter,
                t5_poly,
                ..
            } => t5_adapter + t5_poly,
        }
    }
    /// Total 3′ trim (adapter + polyX).
    pub fn t3(&self) -> usize {
        match self {
            TrimResult::Drop {
                t3_adapter,
                t3_poly,
                ..
            } => t3_adapter + t3_poly,
            TrimResult::Emit {
                t3_adapter,
                t3_poly,
                ..
            } => t3_adapter + t3_poly,
        }
    }
    /// Total bases removed on a kept read (0 on a drop).
    pub fn total_trim(&self) -> usize {
        match self {
            TrimResult::Drop { .. } => 0,
            TrimResult::Emit {
                t5_adapter,
                t3_adapter,
                t5_poly,
                t3_poly,
                runs,
            } => {
                t5_adapter
                    + t3_adapter
                    + t5_poly
                    + t3_poly
                    + runs.iter().map(|(_, l)| l).sum::<usize>()
            }
        }
    }
}

/// Resolve trim length from a scan result, applying the m_h threshold
/// and the boundary +1 correction.
///
/// Must not be called with `StopMode::Budget` — callers short-circuit
/// that case to `TrimResult::Drop`.
///
/// The homopolymer-length threshold is applied to the match count, not
/// to `K*`: a scan can have a large `K*` when mercy fires after only a
/// few matches (absorbed by high-quality mismatches at the scan's end),
/// and that should not count as a qualifying run. The actual run length
/// is `matches + 1` (the anchor base plus one matching neighbour per
/// matching step).
#[inline]
fn resolve_trim(res: ScanResult, min_run: usize, boundary_anchored: bool) -> usize {
    debug_assert!(res.mode != StopMode::Budget);
    if res.matches + 1 < min_run {
        0
    } else if boundary_anchored {
        // Cut the homopolymer run extent (start ..= last matching base),
        // dropping the flanking overshoot the mercy walk absorbs.
        res.run_span
    } else if res.mode == StopMode::Mercy {
        res.k_star
    } else {
        res.k_star + 1 // boundary
    }
}

/// Output of the streaming walk: state needed by the rest of
/// `polyx_trim`. Bundled so we can return it from a const-generic
/// helper without polluting the call site with tuples.
struct StreamOut {
    quality_sum: u64,
    first_trim: usize,
    h_int: usize,
    runs: Vec<(usize, usize)>,
    /// True if the streaming pass exhausted the trim budget on the
    /// first sub-scan (5′ short-circuit equivalent of `StopMode::Budget`).
    dropped_by_budget: bool,
}

/// Inner streaming walk over `[0, l)`. Monomorphized on direction so
/// that `i_phys` and `prev_phys` collapse to constant-time index
/// expressions and the per-iteration direction check disappears.
///
/// `FORWARD = true` means the streaming walk runs in the natural
/// 5′→3′ direction; `FORWARD = false` walks the read in reverse and
/// is used when the user asked for `Direction::Reverse` (typical
/// for R2 polyG short-circuiting).
#[inline]
fn streaming_walk<const FORWARD: bool>(
    bases: &[u8],
    quals: &[u8],
    params: &Params,
    win_start: usize,
    win_end: usize,
    rem_budget: usize,
) -> StreamOut {
    let l = bases.len();

    let mut quality_sum: u64 = 0;
    let mut gamma = 0.0f64;
    let mut k_curr: usize = 0;
    let mut m_curr: usize = 0;
    let mut run_span_curr: usize = 0;
    let mut sub_start_logical: usize = win_start;
    let mut first_subscan = true;

    let mut first_trim: usize = 0;
    let mut h_int: usize = 0;
    let mut runs: Vec<(usize, usize)> = Vec::new();

    #[inline]
    fn resolve_subscan(
        k_star: usize,
        m_star: usize,
        mode: StopMode,
        min_run: usize,
        run_span: usize,
        boundary_anchored: bool,
    ) -> usize {
        debug_assert!(mode != StopMode::Budget);
        if m_star + 1 < min_run {
            0
        } else if boundary_anchored {
            run_span
        } else if mode == StopMode::Mercy {
            k_star
        } else {
            k_star + 1
        }
    }

    // Compile-time-resolved index translation. Because FORWARD is a
    // const generic parameter, the optimizer monomorphizes this
    // helper into either an identity or `n - 1 - i`, eliminating the
    // per-step branch we previously paid 10% of total runtime for.
    #[inline(always)]
    fn idx_const<const FWD: bool>(i: usize, n: usize) -> usize {
        if FWD {
            i
        } else {
            n - 1 - i
        }
    }

    for i_logical in 0..l {
        let i_phys = idx_const::<FORWARD>(i_logical, l);

        // Always sum quality (mean-q is over the whole read).
        // Safety: i_phys < l by construction (idx_const returns either
        // i_logical < l or l - 1 - i_logical, both < l). The compiler
        // should hoist the bounds check; if not, we can switch to
        // get_unchecked once the safety argument is mechanically clear.
        quality_sum += quals[i_phys] as u64;

        // Polyx window check: pair (b[i-1], b[i]) requires
        // i_logical in [win_start + 1, win_end - 1].
        if i_logical < win_start + 1 || i_logical >= win_end {
            continue;
        }

        let prev_phys = idx_const::<FORWARD>(i_logical - 1, l);
        let (gap, is_match) = per_step_gap_and_match(
            bases[prev_phys],
            bases[i_phys],
            quals[prev_phys],
            quals[i_phys],
            params,
            &params.rule,
        );

        gamma += gap;
        if is_match {
            m_curr += 1;
            // Run reaches the pair's far base `i_logical`; extent from the
            // sub-scan anchor is (i_logical - sub_start_logical + 1) bases.
            // Only tracked for the boundary-anchored cut.
            if params.boundary_anchored {
                run_span_curr = i_logical - sub_start_logical + 1;
            }
        }
        k_curr += 1;

        // Streaming 5′ budget short-circuit.
        if first_subscan && rem_budget > 0 && k_curr >= rem_budget {
            return StreamOut {
                quality_sum, // partial; caller still drops, doesn't use
                first_trim: 0,
                h_int: 0,
                runs: Vec::new(),
                dropped_by_budget: true,
            };
        }

        if gamma > params.mercy {
            let sub_trim = resolve_subscan(
                k_curr,
                m_curr,
                StopMode::Mercy,
                params.min_run,
                run_span_curr,
                params.boundary_anchored,
            );
            if first_subscan {
                first_trim = sub_trim;
                first_subscan = false;
            } else if sub_trim > 0 && params.scan_interior {
                let sub_pos = sub_start_logical - win_start;
                h_int += sub_trim;
                runs.push((sub_pos, sub_trim));
                if first_trim + h_int > rem_budget {
                    return StreamOut {
                        quality_sum,
                        first_trim,
                        h_int,
                        runs,
                        dropped_by_budget: true,
                    };
                }
            }
            gamma = 0.0;
            k_curr = 0;
            m_curr = 0;
            run_span_curr = 0;
            sub_start_logical = i_logical + 1;
        }
    }

    // End of walk: handle in-flight sub-scan. Boundary mode for the
    // 5′-end scan resolves into first_trim. Interior boundary-mode
    // sub-scans are not recorded (spec §1.9).
    if k_curr > 0 && first_subscan {
        first_trim = resolve_subscan(
            k_curr,
            m_curr,
            StopMode::Boundary,
            params.min_run,
            run_span_curr,
            params.boundary_anchored,
        );
    }

    StreamOut {
        quality_sum,
        first_trim,
        h_int,
        runs,
        dropped_by_budget: false,
    }
}

/// Two-tier cascade, tier 1 (paper §11, Prop "rejection test").
///
/// Quality-blind per-end rejection: the end-scan of a window visits its
/// first `w = m_h + floor(Δ_m) - 1` positions first. If those contain
/// fewer than `m_h - 1` adjacent matches (the generic-case anchor count,
/// paper Def "anchor count"), the scan terminates by mercy within the
/// window (Prop "prefilter soundness") with `m* ≤ c_X < m_h - 1`, so it
/// resolves to `T = 0` (Prop "rejection test"). O(1): only the end's
/// terminal window is examined, not the whole read.
///
/// Returns `true` when the end-scan starting at `lo` (walking toward
/// `hi`) can be skipped with `T = 0`. One-sided: may pass windows the
/// Γ-scan later rejects, never rejects one it would trim.
#[inline]
fn prefilter_end_rejects(bases: &[u8], lo: usize, hi: usize, w: usize, min_run: usize) -> bool {
    if w < 2 || hi.saturating_sub(lo) < w + 1 {
        return false; // degenerate window or too short to matter
    }
    let need = min_run - 1;
    let mut matches = 0usize;
    for i in lo..lo + w {
        if bases[i] == bases[i + 1] {
            matches += 1;
        }
    }
    matches < need
}

/// Two-tier cascade, tier 1, interior gate (paper §11 Rem "one-sided
/// error"). Interior sub-scans start wherever a mercy exit lands, so any
/// window of `w` adjacent pairs with `≥ m_h - 1` matches could host one.
/// Rejecting the interior therefore requires the whole-read check: every
/// window of `w` pairs (spanning `w + 1` positions) in `[lo, hi)` has
/// fewer than `m_h - 1` matches. O(L) sliding pass, byte compares only.
#[inline]
fn prefilter_interior_rejects(
    bases: &[u8],
    lo: usize,
    hi: usize,
    w: usize,
    min_run: usize,
) -> bool {
    if w < 2 || hi.saturating_sub(lo) < w + 1 {
        return false;
    }
    let need = min_run - 1;
    let mut matches = 0usize;
    for i in lo..lo + w {
        if bases[i] == bases[i + 1] {
            matches += 1;
        }
    }
    if matches >= need {
        return false;
    }
    for s in lo + 1..=hi - (w + 1) {
        // Slide one position: drop pair (s-1, s), add pair (s+w-1, s+w).
        if bases[s - 1] == bases[s] {
            matches -= 1;
        }
        if bases[s + w - 1] == bases[s + w] {
            matches += 1;
        }
        if matches >= need {
            return false;
        }
    }
    true
}

/// Test-only kill switch for the two-tier cascade's tier 1. The conformance
/// test runs the full walk (prefilter disabled) and the cascade side by side
/// over randomized reads and asserts identical TrimResults — the paper's
/// `cascade(R,Q) = PolyXTrim(R,Q)` property.
#[cfg(test)]
static FORCE_NO_PREFILTER: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

/// Whether tier 1 may skip the walk. Production is always enabled (inlined
/// to a constant `true`); the test build can disable it via
/// [`FORCE_NO_PREFILTER`] for the conformance comparison.
#[inline(always)]
fn prefilter_enabled() -> bool {
    #[cfg(test)]
    {
        use std::sync::atomic::Ordering;
        !FORCE_NO_PREFILTER.load(Ordering::Relaxed)
    }
    #[cfg(not(test))]
    {
        true
    }
}

/// Run the full pipeline on a single read.
///
/// `bases` and `quals` must have equal length. `t5_poly` and `t3_poly` in the
/// returned [`TrimResult::Emit`] are the polyX 5′ and 3′ trim counts,
/// including any hard-trim contributions (adapter trims are always 0 here).
///
/// # Streaming structure
///
/// The forward direction is a single walk over `[0, L)` that does two
/// things at once: it sums the full-read quality (for the mean-q
/// filter) and feeds pairs into the polyx scan state machine, but only
/// for indices inside the polyx window `[t_front, l - t_tail]`. This
/// fuses what used to be two distinct passes (mean-q pre-pass + scan)
/// into a single byte walk. The walk is monomorphized on direction
/// (Forward / Reverse) via [`streaming_walk`], which eliminates a
/// per-iteration branch that profiling showed at ~10% of runtime.
///
/// The reverse pass runs unchanged after, and is bounded by the trim
/// budget so it walks at most `Δ_L'` bases.
///
/// `Direction::Reverse` flips which physical end the streaming pass
/// owns: in `Forward` mode the streaming pass produces `T_5` (and any
/// interior runs); in `Reverse` mode it produces `T_3` (with the
/// streaming walk happening in the reversed frame).
pub fn polyx_trim(bases: &[u8], quals: &[u8], params: &Params) -> TrimResult {
    polyx_trim_inner(
        bases,
        quals,
        params,
        params.extra_committed,
        params.original_len,
    )
}

/// Inner implementation. `extra_committed` and `original_len` are passed explicitly
/// so `trim()` can invoke it without cloning `Params` (avoiding an `Arc` bump on
/// every adapter-trimmed read — ~15 M atomic ops on a typical Illumina lane).
fn polyx_trim_inner(
    bases: &[u8],
    quals: &[u8],
    params: &Params,
    extra_committed: usize,
    original_len: Option<usize>,
) -> TrimResult {
    assert_eq!(bases.len(), quals.len());
    let l = bases.len();

    // Preserve existing behavior for very short inputs.
    if l < 2 {
        return TrimResult::Emit {
            t5_adapter: 0,
            t3_adapter: 0,
            t5_poly: 0,
            t3_poly: 0,
            runs: vec![],
        };
    }

    let delta_l = {
        let l_eff = original_len.unwrap_or(l);
        l_eff
            .saturating_sub(params.min_len.resolve(l_eff))
            .saturating_sub(extra_committed)
    };

    // Hard trims, capped to the read length.
    let t_front = params.trim_5.unwrap_or(0).min(l);
    let t_tail = params.trim_3.unwrap_or(0).min(l - t_front);

    // NovaSeq first-cycle 5′ G-dropout: a short leading G-run with normal
    // quality. The quality-gated scan can't see it (it looks for the
    // quality collapse of a long run, and min_run rejects short ones), so
    // it is trimmed here as a hard trim — exactly the run, no mercy walk.
    // Only applied when the 5′ scan would otherwise be active and the run
    // is below the scan's threshold (params.head_dropout < min_run is
    // enforced by the CLI).
    let head_trim = if params.scan_5 {
        match params.head_dropout {
            Some(n) => {
                let run = bases[..l].iter().take_while(|&&b| b == b'G').count();
                if (1..=n).contains(&run) {
                    run
                } else {
                    0
                }
            }
            None => 0,
        }
    } else {
        0
    };
    let t_front = t_front + head_trim;
    let sub_l = l - t_front - t_tail;

    // Hard trim alone consumes the read, or leaves nothing scannable.
    if sub_l < 2 {
        return TrimResult::Drop {
            reason: DropReason::Hard,
            t5_adapter: 0,
            t3_adapter: 0,
            t5_poly: t_front,
            t3_poly: t_tail,
            runs: vec![],
        };
    }
    if t_front + t_tail > delta_l {
        return TrimResult::Drop {
            reason: DropReason::Hard,
            t5_adapter: 0,
            t3_adapter: 0,
            t5_poly: t_front,
            t3_poly: t_tail,
            runs: vec![],
        };
    }
    let rem_budget = delta_l - (t_front + t_tail);

    // ----- Two-tier cascade, tier 1 (paper §11, Prop "rejection test") -----
    //
    // Quality-blind rejection. The 5′ end-scan can be skipped when its
    // terminal window (the first w positions it visits) has < m_h - 1
    // matches — O(1). The interior needs the whole-read sliding check
    // (interior sub-scans start at mercy exits anywhere). When both the
    // 5′ end and the interior reject, the streaming walk resolves to
    // T_5 = H_int = 0 and is skipped entirely. The 3′ end-scan (reverse
    // pass) is gated separately on its own terminal window.
    //
    // Guard: only skip when rem_budget > w. A rejected scan terminates by
    // mercy within w steps; with a budget cap ≤ w the walk could instead
    // fire the Budget short-circuit and DROP the read, so the skip must
    // not change that outcome.
    let w = params
        .min_run
        .saturating_add(params.mercy.floor() as usize)
        .saturating_sub(1);
    let scan_range_lo = t_front;
    let scan_range_hi = l - t_tail;
    let five_end_rejects = params.scan_5
        && prefilter_enabled()
        && prefilter_end_rejects(bases, scan_range_lo, scan_range_hi, w, params.min_run);
    let interior_rejects = !params.scan_interior
        || prefilter_enabled()
            && prefilter_interior_rejects(bases, scan_range_lo, scan_range_hi, w, params.min_run);
    let three_end_rejects = params.scan_3
        && prefilter_enabled()
        && prefilter_end_rejects(
            bases,
            scan_range_hi.saturating_sub(w + 1),
            scan_range_hi,
            w,
            params.min_run,
        );

    // Direction::Forward  — streaming pass = physical 5′,  reverse pass = physical 3′
    // Direction::Reverse  — streaming pass = physical 3′,  reverse pass = physical 5′
    let dir = params.direction;
    let (stream_enabled, reverse_enabled) = match dir {
        Direction::Forward => (params.scan_5, params.scan_3),
        Direction::Reverse => (params.scan_3, params.scan_5),
    };

    // Which end each pass actually scans, resolved against the physical
    // frame: the streaming pass's terminal window is the 5′ end in
    // Forward mode and the 3′ end in Reverse mode; the reverse pass is
    // the other end.
    let (stream_rejects, reverse_rejects) = match dir {
        Direction::Forward => (five_end_rejects, three_end_rejects),
        Direction::Reverse => (three_end_rejects, five_end_rejects),
    };

    // The tier-1 pass is pure overhead when nothing will be scanned (e.g.
    // the bench `floor` stage: --no-5 --no-3 --no-interior). Only compute
    // it when at least one scan is enabled. Verified by profile: with all
    // scans off, `composition_prefilter_rejects` was 38% of floor runtime.
    let prefiltered =
        stream_enabled && w >= 2 && rem_budget > w && stream_rejects && interior_rejects;

    // ----- Streaming forward pass -----
    //
    // The streaming walk is dispatched to one of two specialized
    // functions, monomorphized on the `FORWARD` const parameter. This
    // eliminates the per-iteration direction check that perf annotation
    // showed at >10% of total runtime: in `Forward` mode `i_phys = i`
    // and the bounds checks on `bases[i]`/`quals[i]` simplify to
    // single-comparison patterns the optimizer can hoist; in `Reverse`
    // mode `i_phys = n - 1 - i` is constant-folded the same way.
    let stream_out = if prefiltered {
        // quality_sum still needed for the mean-q filter when configured;
        // the bench never sets --mean-q, so this is usually a plain 0.
        StreamOut {
            quality_sum: if params.mean_q.is_some() {
                quals.iter().map(|&q| q as u64).sum()
            } else {
                0
            },
            first_trim: 0,
            h_int: 0,
            runs: Vec::new(),
            dropped_by_budget: false,
        }
    } else if !stream_enabled {
        // Quality-only walk (no scan). Still cheap; one path with no
        // direction branch.
        StreamOut {
            quality_sum: quals.iter().map(|&q| q as u64).sum(),
            first_trim: 0,
            h_int: 0,
            runs: Vec::new(),
            dropped_by_budget: false,
        }
    } else if matches!(dir, Direction::Forward) {
        streaming_walk::<true>(bases, quals, params, t_front, l - t_tail, rem_budget)
    } else {
        streaming_walk::<false>(bases, quals, params, t_front, l - t_tail, rem_budget)
    };

    if stream_out.dropped_by_budget {
        let reason = match dir {
            Direction::Forward => DropReason::Poly5,
            Direction::Reverse => DropReason::Poly3,
        };
        return TrimResult::Drop {
            reason,
            t5_adapter: 0,
            t3_adapter: 0,
            t5_poly: t_front,
            t3_poly: t_tail,
            runs: vec![],
        };
    }
    let first_trim = stream_out.first_trim;
    let h_int = stream_out.h_int;
    let mut runs = stream_out.runs;
    let quality_sum = stream_out.quality_sum;

    // ----- Mean-quality filter -----
    //
    // Computed over the whole read, after the streaming pass. Doing it
    // here instead of up front means we may have done some "wasted"
    // scan work on a soon-to-be-dropped read, but mean-q failure rates
    // are typically 1-5% on Illumina data so this is comfortably a win.
    if let Some(threshold) = params.mean_q {
        let mean = quality_sum as f64 / l as f64;
        if mean < threshold {
            return TrimResult::Drop {
                reason: DropReason::MeanQ,
                t5_adapter: 0,
                t3_adapter: 0,
                t5_poly: 0,
                t3_poly: 0,
                runs: vec![],
            };
        }
    }

    // ----- Reverse pass (unchanged single bounded call) -----
    //
    // Operates on the same sub-slice as the streaming forward pass but
    // in the opposite frame. The budget cap ensures it walks at most
    // (rem_budget - first_trim) bases.
    let sub_bases = &bases[t_front..l - t_tail];
    let sub_quals = &quals[t_front..l - t_tail];
    let upper = sub_l - 1;

    let reverse_budget = rem_budget.saturating_sub(first_trim);
    let second_trim = if reverse_enabled && reverse_rejects && reverse_budget > w {
        // Tier-1 rejection covers the reverse end-scan: its terminal window
        // has < m_h - 1 matches, so the pass would resolve to 0. Guard: the
        // pass's own budget cap is (rem_budget - first_trim), and a cap
        // ≤ w could fire the Budget short-circuit and DROP the read (the
        // rejection test only guarantees mercy fires within w steps when
        // the cap does not preempt it).
        0
    } else if reverse_enabled {
        let second_dir = match dir {
            Direction::Forward => Direction::Reverse,
            Direction::Reverse => Direction::Forward,
        };
        let budget = rem_budget.saturating_sub(first_trim);
        let budget_opt = if budget > 0 { Some(budget) } else { None };
        let res = scan(
            sub_bases, sub_quals, 0, upper, second_dir, params, budget_opt,
        );
        if res.mode == StopMode::Budget {
            let reason = match dir {
                Direction::Forward => DropReason::Poly3,
                Direction::Reverse => DropReason::Poly5,
            };
            let (scan_t5_poly, scan_t3_poly) = match dir {
                Direction::Forward => (first_trim, 0),
                Direction::Reverse => (0, first_trim),
            };
            return TrimResult::Drop {
                reason,
                t5_adapter: 0,
                t3_adapter: 0,
                t5_poly: t_front + scan_t5_poly,
                t3_poly: t_tail + scan_t3_poly,
                runs: vec![],
            };
        }
        resolve_trim(res, params.min_run, params.boundary_anchored)
    } else {
        0
    };

    // ----- Reconciliation -----
    //
    // Remove interior runs whose recorded position overlaps the 3'
    // (post-hard-trim) trim zone. Positions are in the sub-slice frame
    // (relative to win_start); the cut is at sub_l - second_trim.
    let mut h_int_adj = h_int;
    if second_trim > 0 {
        let cut = sub_l.saturating_sub(second_trim);
        let mut drop_count = 0usize;
        for &(start, len) in runs.iter().rev() {
            if start + len > cut {
                h_int_adj = h_int_adj.saturating_sub(len);
                drop_count += 1;
            } else {
                break;
            }
        }
        runs.truncate(runs.len() - drop_count);
    }

    if first_trim + second_trim + h_int_adj > rem_budget {
        let (scan_t5_poly, scan_t3_poly) = match dir {
            Direction::Forward => (first_trim, second_trim),
            Direction::Reverse => (second_trim, first_trim),
        };
        let hard = t_front + t_tail;
        let max_poly = scan_t5_poly.max(scan_t3_poly).max(h_int_adj);
        let reason = if max_poly >= hard {
            if scan_t5_poly >= scan_t3_poly && scan_t5_poly >= h_int_adj {
                DropReason::Poly5
            } else if scan_t3_poly >= h_int_adj {
                DropReason::Poly3
            } else {
                DropReason::PolyInt
            }
        } else {
            DropReason::Hard
        };
        return TrimResult::Drop {
            reason,
            t5_adapter: 0,
            t3_adapter: 0,
            t5_poly: t_front + scan_t5_poly,
            t3_poly: t_tail + scan_t3_poly,
            runs,
        };
    }

    let (scan_t5_poly, scan_t3_poly) = match dir {
        Direction::Forward => (first_trim, second_trim),
        Direction::Reverse => (second_trim, first_trim),
    };

    let abs_runs: Vec<(usize, usize)> = runs
        .into_iter()
        .map(|(sub_pos, len)| match dir {
            Direction::Forward => (t_front + sub_pos, len),
            Direction::Reverse => (l - (t_front + sub_pos) - len, len),
        })
        .collect();

    TrimResult::Emit {
        t5_adapter: 0,
        t3_adapter: 0,
        t5_poly: t_front + scan_t5_poly,
        t3_poly: t_tail + scan_t3_poly,
        runs: abs_runs,
    }
}

// ---------------------------------------------------------------------------
// Public trim() — adapter + polyX in one pass.
// ---------------------------------------------------------------------------

/// Output of [`trim`]: the trim decision plus the name of the matched
/// adapter (if any), used by the pipeline to update per-adapter stats.
#[derive(Debug)]
pub struct TrimOutput {
    pub result: TrimResult,
    /// Name of the adapter that was trimmed, if any.
    pub adapter_name: Option<String>,
}

/// Full single-read processing: adapter detection then polyX trimming.
///
/// Steps (single byte-stream pass design):
/// 1. Adapter detection by sequence matching from `params.adapter` — the only
///    mode. The PE-overlap override path and the `AdapterHint` dispatcher that
///    carried its result were removed with the feature (adapter.rs).
/// 2. Budget check: if adapter trim alone exceeds Δ_L, return
///    `Drop { reason: Adapter }`.
/// 3. polyX scan on the adapter-trimmed sub-slice with the remaining
///    budget enforced via `polyx_trim_inner` (no `Params` clone — the
///    `extra_committed` and `original_len` overrides are passed explicitly).
/// 4. Reconstruct absolute coordinates (add `t5_adapter` back to all
///    interior-run positions).
pub fn trim(bases: &[u8], quals: &[u8], params: &Params) -> TrimOutput {
    let l = bases.len();

    // --- Step 1: adapter detection (sequence matching) ---
    let (t5_adapter, t3_adapter, adapter_name) = match params.adapter.as_deref() {
        None => (0, 0, None),
        Some(set) => match crate::adapter::match_adapters(bases, set) {
            None => (0, 0, None),
            Some(m) => (m.t5, m.t3, Some(m.name)),
        },
    };

    // --- Step 2: adapter-only budget check ---
    let committed = t5_adapter + t3_adapter;
    if committed > 0 {
        let delta_l = params.trim_budget(l); // uses original l
        if committed > delta_l {
            return TrimOutput {
                result: TrimResult::Drop {
                    reason: DropReason::Adapter,
                    t5_adapter,
                    t3_adapter,
                    t5_poly: 0,
                    t3_poly: 0,
                    runs: vec![],
                },
                adapter_name,
            };
        }
    }

    // --- Step 3: polyX scan on adapter-trimmed sub-slice ---
    // Guard against pathological adapter trim consuming the whole read.
    let safe_t5 = t5_adapter.min(l);
    let safe_t3 = t3_adapter.min(l.saturating_sub(safe_t5));
    let sub_bases = &bases[safe_t5..l - safe_t3];
    let sub_quals = &quals[safe_t5..l - safe_t3];

    // Call the inner implementation directly, passing `committed` and `Some(l)` so
    // the budget is enforced relative to the original read length — no Params clone.
    let poly_result = polyx_trim_inner(sub_bases, sub_quals, params, committed, Some(l));

    // --- Step 4: reconstruct absolute coordinates ---
    let result = match poly_result {
        TrimResult::Drop {
            reason,
            t5_poly,
            t3_poly,
            runs,
            ..
        } => TrimResult::Drop {
            reason,
            t5_adapter,
            t3_adapter,
            t5_poly,
            t3_poly,
            runs: runs
                .into_iter()
                .map(|(pos, len)| (pos + safe_t5, len))
                .collect(),
        },
        TrimResult::Emit {
            t5_poly,
            t3_poly,
            runs,
            ..
        } => {
            // Interior runs were recorded relative to the polyX sub-slice;
            // add t5_adapter to bring them into original-read coordinates.
            let abs_runs = runs
                .into_iter()
                .map(|(pos, len)| (pos + safe_t5, len))
                .collect();
            TrimResult::Emit {
                t5_adapter,
                t3_adapter,
                t5_poly,
                t3_poly,
                runs: abs_runs,
            }
        }
    };

    TrimOutput {
        result,
        adapter_name,
    }
}
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a read with constant Phred quality `q` over bases `b`.
    fn read(b: &[u8], q: u8) -> (Vec<u8>, Vec<u8>) {
        (b.to_vec(), vec![q; b.len()])
    }

    // ------- MinLen parsing -----------------------------------------------

    #[test]
    fn minlen_parse_int() {
        assert_eq!(MinLen::from_str("100"), Ok(MinLen::Abs(100)));
        assert_eq!(MinLen::from_str("0"), Ok(MinLen::Abs(0)));
    }

    #[test]
    fn minlen_parse_frac() {
        assert_eq!(MinLen::from_str("0.65L"), Ok(MinLen::Frac(0.65)));
        assert_eq!(MinLen::from_str("1.0L"), Ok(MinLen::Frac(1.0)));
        assert_eq!(MinLen::from_str("0.5l"), Ok(MinLen::Frac(0.5)));
    }

    #[test]
    fn minlen_parse_rejects_bad() {
        assert!(MinLen::from_str("-5").is_err());
        assert!(MinLen::from_str("2.0L").is_err());
        assert!(MinLen::from_str("0L").is_err());
        assert!(MinLen::from_str("abc").is_err());
        assert!(MinLen::from_str("0.5X").is_err()); // wrong suffix
    }

    #[test]
    fn minlen_resolve() {
        assert_eq!(MinLen::Abs(100).resolve(150), 100);
        assert_eq!(MinLen::Abs(200).resolve(150), 150); // clamped
        assert_eq!(MinLen::Frac(0.6).resolve(150), 90);
        assert_eq!(MinLen::Frac(0.65).resolve(150), 97); // floor(97.5)
        assert_eq!(MinLen::Frac(1.0).resolve(150), 150);
    }

    // ------- Quality gate -------------------------------------------------

    #[test]
    fn quality_gate_boundaries() {
        assert_eq!(quality_gate(0, 20, 1.0), 0.0);
        assert_eq!(quality_gate(20, 20, 1.0), 1.0);
        assert_eq!(quality_gate(40, 20, 1.0), 1.0); // saturated
        assert!((quality_gate(10, 20, 1.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn quality_gate_shape() {
        // x = 2 concentrates weight near t
        let g1 = quality_gate(10, 20, 2.0);
        assert!((g1 - 0.25).abs() < 1e-9);
        // x = 0.5 admits sub-threshold
        let g2 = quality_gate(10, 20, 0.5);
        assert!((g2 - (0.5f64).sqrt()).abs() < 1e-9);
    }

    #[test]
    fn gate_lut_matches_direct_quality_gate() {
        // Regression test for the LUT optimisation: after prepare(),
        // every entry must equal a direct quality_gate call. If this
        // ever drifts, hot-path scan results would silently change.
        let params = Params::default();
        for q in 0u8..=255 {
            let direct = quality_gate(q, params.phred_sat, params.gate_shape);
            assert_eq!(
                params.gate_lut[q as usize], direct,
                "LUT mismatch at q={}",
                q
            );
        }
    }

    #[test]
    fn gate_lut_updates_after_prepare_with_changed_shape() {
        // Mutating phred_sat/gate_shape requires prepare() to refresh.
        let mut params = Params::default();
        params.phred_sat = 30;
        params.gate_shape = 2.0;
        params.prepare();
        for q in 0u8..=255 {
            assert_eq!(
                params.gate_lut[q as usize],
                quality_gate(q, 30, 2.0),
                "LUT not updated after prepare() at q={}",
                q
            );
        }
    }

    // ------- Scan primitive -----------------------------------------------

    #[test]
    fn scan_pure_polya_boundary() {
        // All-A, high quality: matches all the way, mercy never fires.
        let (b, q) = read(b"AAAAAAAAAA", 40);
        let params = Params::default();
        let res = scan(&b, &q, 0, 9, Direction::Forward, &params, None);
        assert_eq!(res.mode, StopMode::Boundary);
        assert_eq!(res.k_star, 9);
    }

    #[test]
    fn scan_random_mercy() {
        // Highly variable sequence, high quality: mercy fires fast.
        let (b, q) = read(b"ACGTACGTAC", 40);
        let params = Params::default();
        let res = scan(&b, &q, 0, 9, Direction::Forward, &params, None);
        assert_eq!(res.mode, StopMode::Mercy);
        // Δ_m = 3, ν = 0.5, W = 1 ⇒ each mismatch adds 1.5.
        // Γ crosses 3 after step 3 (Γ = 4.5). K* = 3.
        assert_eq!(res.k_star, 3);
    }

    #[test]
    fn scan_polya_then_random() {
        // Seven A's then a mismatch to C, rest random.
        let (b, q) = read(b"AAAAAAACGTACGTAC", 40);
        let params = Params::default();
        let res = scan(&b, &q, 0, b.len() - 1, Direction::Forward, &params, None);
        assert_eq!(res.mode, StopMode::Mercy);
        // Steps 0..5 are matches (A==A, gap ≈ 0 each).
        // Step 6 mismatch A→C adds 1.5.
        // Step 7 mismatch C→G adds 1.5. Γ = 3.0, not > 3.
        // Step 8 mismatch G→T adds 1.5. Γ = 4.5 > 3 ⇒ mercy at K* = 9.
        assert_eq!(res.k_star, 9);
    }

    #[test]
    fn scan_budget_fires() {
        // All-A of length 20, budget 5 ⇒ Budget mode after 5 matches.
        let (b, q) = read(b"AAAAAAAAAAAAAAAAAAAA", 40);
        let params = Params::default();
        let res = scan(&b, &q, 0, b.len() - 1, Direction::Forward, &params, Some(5));
        assert_eq!(res.mode, StopMode::Budget);
        assert_eq!(res.k_star, 5);
    }

    #[test]
    fn scan_reverse_direction() {
        // polyG at the 3′ end only.
        let (b, q) = read(b"ACGTACGTGGGGGGGG", 40);
        let params = Params::default();
        let res = scan(&b, &q, 0, b.len() - 1, Direction::Reverse, &params, None);
        // Reverse frame sees 8 G's at the start, then random.
        assert_eq!(res.mode, StopMode::Mercy);
        // Similar math: 7 matches (gap 0), then mismatches accumulate.
        assert!(res.k_star >= 7);
    }

    #[test]
    fn scan_low_quality_mismatch_weak() {
        // Low-quality mismatch should not fire mercy as fast as hi-Q.
        let (b, q) = read(b"ACGTACGTACGT", 2); // q=2 well below t=20
        let params = Params::default();
        let res = scan(&b, &q, 0, b.len() - 1, Direction::Forward, &params, None);
        // W ≈ 0.1, mismatch gap ≈ 1.05 each. Still fires mercy, but later-ish;
        // main point: reaches mercy within reasonable K.
        assert_eq!(res.mode, StopMode::Mercy);
    }

    // ------- polyx_trim: full pipeline ------------------------------------

    #[test]
    fn trim_clean_read_unchanged() {
        // Random sequence, no homopolymer at either end.
        let (b, q) = read(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", 40);
        let params = Params::default();
        let res = polyx_trim(&b, &q, &params);
        match res {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 0);
                assert_eq!(t3_poly, 0);
                assert!(runs.is_empty(), "expected no interior runs");
            }
            _ => panic!("should not drop clean read, got {:?}", res),
        }
    }

    #[test]
    fn trim_5prime_polyg_trimmed() {
        // 10 G's then random tail, length 52.
        let mut b = vec![b'G'; 10];
        b.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 40 more
        let q = vec![40u8; b.len()];
        let params = Params::default();
        let res = polyx_trim(&b, &q, &params);
        match res {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                // 10 matches, then mismatches fire mercy. T_5 ≈ 10 + a few.
                assert!(t5_poly >= 10, "t5_poly = {}", t5_poly);
                assert!(t5_poly < 15);
                assert_eq!(t3_poly, 0);
                assert!(runs.is_empty());
            }
            _ => panic!("unexpected drop"),
        }
    }

    #[test]
    fn trim_all_polya_dropped() {
        // All-A of length 150, Δ_L ≈ 60: budget fires on 5′ end-scan.
        let (b, q) = read(&vec![b'A'; 150], 40);
        let params = Params::default();
        assert!(matches!(
            polyx_trim(&b, &q, &params),
            TrimResult::Drop { .. }
        ));
    }

    #[test]
    fn trim_3prime_polyg() {
        // Random prefix then polyG tail.
        let mut b = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(); // 40
        b.extend(std::iter::repeat(b'G').take(12));
        let q = vec![40u8; b.len()];
        let params = Params::default();
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 0);
                assert!(t3_poly >= 12, "t3_poly = {}", t3_poly);
                assert!(t3_poly < 17);
                assert!(runs.is_empty());
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn trim_interior_run_accumulates() {
        // Random + interior polyA(12) + random.
        let mut b = b"ACGTACGTACGTACGT".to_vec(); // 16
        b.extend(std::iter::repeat(b'A').take(12)); // interior run
        b.extend_from_slice(b"CGTACGTACGTACGT"); // 15 more — ends with T so no 3' run
        let q = vec![40u8; b.len()];
        let params = Params::default();
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 0);
                assert_eq!(t3_poly, 0);
                let h_int: usize = runs.iter().map(|(_, l)| l).sum();
                assert!(h_int >= 11, "h_int = {}", h_int);
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn trim_reconciliation_preserves_non_overlapping_run() {
        // An interior run that sits entirely *before* the 3' trim zone
        // must survive reconciliation. This is the load-bearing invariant:
        // a false positive (removing a non-overlapping run) would inflate
        // the drop rate.
        //
        // Layout: 16 clean | 10 A interior | 23 clean middle | 10 G tail.
        let mut b = vec![];
        b.extend_from_slice(b"ACGTACGTACGTACGT"); // 16
        b.extend(std::iter::repeat(b'A').take(10));
        b.extend_from_slice(b"CGTACGTACGTACGTACGTACGT"); // 23
        b.extend(std::iter::repeat(b'G').take(10));
        let q = vec![40u8; b.len()];
        let params = Params::default();
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 0);
                assert!(t3_poly >= 10, "t3_poly = {}, len = {}", t3_poly, b.len());
                // Interior A-run is well before the 3' zone → must be kept.
                let h_int: usize = runs.iter().map(|(_, l)| l).sum();
                assert!(h_int >= 9, "h_int = {} (A-run should survive)", h_int);
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn trim_reconciliation_removes_overlapping_run() {
        // An interior run whose last position lies in [L - T_3, L - 1]
        // must be removed from H_int.
        //
        // Forcing the scenario requires both:
        //   (a) the forward interior scan mercy-terminates on the run
        //       (needs a mismatch *after* the run so it doesn't hit EOR),
        //   (b) the reverse scan's T_3 reaches back through that mismatch.
        //
        // We use a larger min_len fraction so Δ_L is small but T_3 can
        // still land in the right place, and a small enough m_h so short
        // runs qualify.
        //
        // Layout: 30 clean | 8 G | 1 T | clean tail of length 1.
        //   L = 40. With Δ_L = 20 (min_len = 0.5) this is plenty.
        //   Forward 5': exits in ~4 mismatch steps, T_5 = 0.
        //   Forward interior: mercy-terminates on the G-run at pos ~30,
        //     records (30, ~10). Advances past mercy-end.
        //   Reverse 3': starts at the final clean base, absorbs the T
        //     mismatch, walks back through 8 G's, keeps going into clean
        //     region, eventually mercy. T_3 straddles or covers the G-run.
        let mut b = vec![];
        b.extend_from_slice(b"CGTACGTACGTACGTACGTACGTACGTACG"); // 30 clean
        b.extend(std::iter::repeat(b'G').take(8));
        b.push(b'T');
        b.push(b'A'); // final clean base, total L = 40
        let q = vec![40u8; b.len()];
        let params = Params {
            min_len: MinLen::Frac(0.5),
            min_run: 4, // allow 8-G run to qualify with some margin
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 0);
                // T_3 is large enough that the G-run's last position
                // (≥ 38) falls into [L - T_3, L - 1] = [40 - T_3, 39].
                // If T_3 ≥ 2 (trivially satisfied) reconciliation triggers.
                assert!(t3_poly >= 2, "t3_poly = {}", t3_poly);
                // Either: reconciliation removed the run (H_int == 0), or
                // the forward scan never recorded it (also H_int == 0).
                // Both satisfy the no-double-counting invariant.
                let h_int: usize = runs.iter().map(|(_, l)| l).sum();
                assert!(
                    t3_poly + h_int <= 11,
                    "no double counting: t3_poly={}, h_int={}",
                    t3_poly,
                    h_int
                );
            }
            TrimResult::Drop { .. } => {
                // Dropping is also acceptable here — what matters is that
                // double-counting doesn't inflate the drop rate when the
                // physical trim is within budget.
            }
        }
    }

    #[test]
    fn trim_anchored_polyg_only() {
        // Anchored to G: a polyA prefix should NOT be trimmed.
        let mut b = vec![b'A'; 12];
        b.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 40
        let q = vec![40u8; b.len()];
        let params = Params {
            rule: MatchRule::anchored(b"G"),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                // Anchored-to-G mode ignores A runs.
                assert_eq!(t5_poly, 0, "anchored-G should not trim polyA");
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn trim_anchored_multi_target() {
        // Anchored to {A, G}: a polyA prefix SHOULD be trimmed.
        let mut b = vec![b'A'; 12];
        b.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"); // 40 C's
        let q = vec![40u8; b.len()];
        let params = Params {
            rule: MatchRule::anchored(b"AG"),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                assert!(t5_poly >= 12, "t5_poly = {}", t5_poly);
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn head_dropout_trims_short_leading_g_run() {
        // NovaSeq first-cycle dropout: 1-base G head, then clean template.
        // The quality-gated scan (min_run 6) leaves it; --head-dropout 5
        // must trim exactly the G run and nothing more.
        let mut b = b"GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let q = vec![40u8; b.len()];
        let params = Params {
            min_run: 6,
            head_dropout: Some(5),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, runs, .. } => {
                assert_eq!(t5_poly, 1, "exactly the 1-base G head");
                assert!(runs.is_empty(), "no interior runs: {:?}", runs);
            }
            r => panic!("unexpected {:?}", r),
        }

        // A 3-base G head is trimmed exactly; longer runs (>= min_run)
        // are left to the scan, and a 6-base run must NOT be touched by
        // the dropout rule.
        let mut b = b"GGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let q = vec![40u8; b.len()];
        let params = Params {
            min_run: 6,
            head_dropout: Some(5),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                assert_eq!(t5_poly, 3, "exactly the 3-base G head");
            }
            r => panic!("unexpected {:?}", r),
        }

        let mut b = b"GGGGGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let q = vec![40u8; b.len()];
        let params = Params {
            min_run: 6,
            head_dropout: Some(5),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                // 6-base run >= min_run: dropout rule must not fire, the
                // scan owns it (t5_poly >= 6).
                assert!(t5_poly >= 6, "scan owns long run, t5_poly = {}", t5_poly);
            }
            r => panic!("unexpected {:?}", r),
        }

        // Without the flag, a 1-base G head is left untouched.
        let mut b = b"GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let q = vec![40u8; b.len()];
        let params = Params {
            min_run: 6,
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                assert_eq!(t5_poly, 0, "no dropout flag: head left alone");
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn boundary_anchored_places_cut_at_run_boundary() {
        // A clean homopolymer run followed by non-run biology: the mercy walk
        // absorbs a couple of flanking bases past the run before Δ_m trips, so
        // the default cut overshoots by ~2. `boundary_anchored` must cut
        // exactly the run extent.
        for (t, min_run) in [(12usize, 12usize), (16, 12), (20, 12)] {
            let l = 150usize;
            let mut seq: Vec<u8> = Vec::new();
            while seq.len() < l - t {
                seq.extend_from_slice(b"ACGT");
            }
            seq.truncate(l - t);
            for _ in 0..t {
                seq.push(b'G');
            }
            let q = vec![40u8; l];
            // default (whole-run/mercy): overshoots the clean run by 2.
            let p = Params {
                min_run,
                scan_interior: false,
                ..Default::default()
            };
            let over = match polyx_trim(&seq, &q, &p) {
                TrimResult::Emit { t3_poly, .. } => t3_poly,
                r => panic!("unexpected {:?}", r),
            };
            assert!(over >= t + 2, "default should overshoot: {} vs {}", over, t);
            // boundary-anchored: exactly the run, no flanking overshoot.
            let p = Params {
                min_run,
                scan_interior: false,
                boundary_anchored: true,
                ..Default::default()
            };
            let cut = match polyx_trim(&seq, &q, &p) {
                TrimResult::Emit { t3_poly, .. } => t3_poly,
                r => panic!("unexpected {:?}", r),
            };
            assert_eq!(cut, t, "boundary-anchored must cut exactly the run");
        }
    }

    #[test]
    fn boundary_anchored_keeps_internal_error_bases() {
        // A run with a substitution error inside it is still one artefact run;
        // the boundary-anchored cut reaches past the internal mismatch (the
        // last match sits beyond it) and removes the whole run.
        let l = 150usize;
        let t = 14usize;
        let mut seq: Vec<u8> = Vec::new();
        while seq.len() < l - t {
            seq.extend_from_slice(b"ACGT");
        }
        seq.truncate(l - t);
        for i in 0..t {
            seq.push(if i == 5 { b'A' } else { b'G' });
        }
        let q = vec![40u8; l];
        let p = Params {
            min_run: 12,
            scan_interior: false,
            boundary_anchored: true,
            ..Default::default()
        };
        match polyx_trim(&seq, &q, &p) {
            TrimResult::Emit { t3_poly, .. } => {
                assert_eq!(t3_poly, t, "must include the internal error base");
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn boundary_anchored_respects_min_run() {
        // A run shorter than min_run is not polyX: no trim either way.
        let l = 150usize;
        let t = 8usize;
        let mut seq: Vec<u8> = Vec::new();
        while seq.len() < l - t {
            seq.extend_from_slice(b"ACGT");
        }
        seq.truncate(l - t);
        for _ in 0..t {
            seq.push(b'G');
        }
        let q = vec![40u8; l];
        let p = Params {
            min_run: 12,
            scan_interior: false,
            boundary_anchored: true,
            ..Default::default()
        };
        match polyx_trim(&seq, &q, &p) {
            TrimResult::Emit { t3_poly, .. } => assert_eq!(t3_poly, 0),
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn trim_scan_disable_toggles() {
        // Disabling 5' scan should leave a 5' polyA run untouched.
        let mut b = vec![b'A'; 15];
        b.extend_from_slice(b"CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA");
        let q = vec![40u8; b.len()];
        let params = Params {
            scan_5: false,
            scan_interior: false, // avoid the head being picked up as interior
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, runs, .. } => {
                assert_eq!(t5_poly, 0);
                assert!(runs.is_empty(), "expected no interior runs, got {:?}", runs);
            }
            r => panic!("unexpected {:?}", r),
        }
    }

    #[test]
    fn trim_very_short_read() {
        // Reads shorter than 2 bases should pass through unchanged.
        let params = Params::default();
        assert!(matches!(
            polyx_trim(b"A", &[40], &params),
            TrimResult::Emit {
                t5_poly: 0,
                t3_poly: 0,
                ..
            }
        ));
        assert!(matches!(
            polyx_trim(b"", &[], &params),
            TrimResult::Emit {
                t5_poly: 0,
                t3_poly: 0,
                ..
            }
        ));
    }

    #[test]
    fn trim_min_run_threshold_rejects_short() {
        // A 4-base run is below default m_h = 6, should not trigger T_5.
        let mut b = vec![b'A'; 4];
        b.extend_from_slice(b"CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA");
        let q = vec![40u8; b.len()];
        let params = Params::default();
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                assert_eq!(t5_poly, 0, "4-base run should be below m_h = 6");
            }
            r => panic!("unexpected {:?}", r),
        }
        // Fractional min_len means Δ_L scales with L.
        let params = Params {
            min_len: MinLen::Frac(0.5),
            ..Default::default()
        };
        // L = 40, Δ_L = 20: a 40-A read exceeds budget ⇒ drop.
        let (b, q) = read(&vec![b'A'; 40], 40);
        assert!(matches!(
            polyx_trim(&b, &q, &params),
            TrimResult::Drop { .. }
        ));
        // L = 96, Δ_L = 48: a 12-A prefix then clean tail ⇒ emit.
        let mut b2 = vec![b'A'; 12];
        b2.extend_from_slice(
            b"CGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCG",
        );
        let q2 = vec![40u8; b2.len()];
        match polyx_trim(&b2, &q2, &params) {
            TrimResult::Emit { t5_poly, .. } => assert!((12..20).contains(&t5_poly)),
            r => panic!("should emit, got {:?}", r),
        }
    }

    #[test]
    fn trim_reverse_direction_swaps_priority() {
        // Construct a read with a huge 3′ polyG and a small 5′ polyA.
        // In Forward mode the 5′ scan runs first on a short run (no trim);
        // the 3′ polyG then exceeds budget ⇒ drop.
        // In Reverse mode the first scan runs from the 3′ end and budget
        // short-circuits on the long polyG immediately ⇒ drop (same
        // outcome), but via a different path.
        let mut b = b"AAACGTACGTACGTACGTACGTACGTAC".to_vec(); // 28
        b.extend(std::iter::repeat(b'G').take(100));
        let q = vec![40u8; b.len()];

        // With default params (Δ_L = 0.4L = 51), 100 G's exceeds budget.
        let res_fwd = polyx_trim(&b, &q, &Params::default());
        let res_rev = polyx_trim(
            &b,
            &q,
            &Params {
                direction: Direction::Reverse,
                ..Default::default()
            },
        );
        assert!(matches!(res_fwd, TrimResult::Drop { .. }));
        assert!(matches!(res_rev, TrimResult::Drop { .. }));

        // Smaller 3′ run that fits: both modes should emit with the same
        // physical t_3 value (trim is end-local, direction only changes
        // order of evaluation).
        let mut b = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(); // 52
        b.extend(std::iter::repeat(b'G').take(10));
        let q = vec![40u8; b.len()];
        let fwd = polyx_trim(&b, &q, &Params::default());
        let rev = polyx_trim(
            &b,
            &q,
            &Params {
                direction: Direction::Reverse,
                ..Default::default()
            },
        );
        match (fwd, rev) {
            (
                TrimResult::Emit {
                    t5_poly: f5,
                    t3_poly: f3,
                    ..
                },
                TrimResult::Emit {
                    t5_poly: r5,
                    t3_poly: r3,
                    ..
                },
            ) => {
                assert_eq!(f5, r5, "t5_poly should be identical across directions");
                assert_eq!(f3, r3, "t3_poly should be identical across directions");
                assert!(f3 >= 10);
            }
            r => panic!("both should emit: {:?}", r),
        }
    }

    #[test]
    fn trim_mean_quality_filter() {
        let b = b"ACGTACGTAC"; // Length 10

        // `polyx_trim` receives raw Phred values (pipeline decodes the
        // ASCII+33 offset from FASTQ input before calling). So these
        // tests use raw Phred directly.
        let q_20 = vec![20u8; b.len()];
        let q_40 = vec![40u8; b.len()];

        let mut params = Params::default();
        params.mean_q = Some(30.0);

        // Threshold 30, read is Phred 20 → drop.
        assert!(matches!(
            polyx_trim(b, &q_20, &params),
            TrimResult::Drop {
                reason: DropReason::MeanQ,
                ..
            }
        ));

        // Threshold 30, read is Phred 40 → keep (clean random-looking seq).
        assert!(matches!(
            polyx_trim(b, &q_40, &params),
            TrimResult::Emit {
                t5_poly: 0,
                t3_poly: 0,
                ..
            }
        ));
    }

    #[test]
    fn trim_mean_q_uses_full_read_with_hard_trims() {
        // Streaming-specific regression test: mean-q must be computed
        // over ALL L bases, not just the polyx window. We construct a
        // read where the polyx-window quality is high but the
        // hard-trimmed flanks have low quality, dragging the whole-read
        // mean below threshold.
        //
        // Length 20 read: 5 low-Q (Phred 5) bases | 10 high-Q (Phred 40)
        // bases | 5 low-Q (Phred 5) bases. Hard-trim 5 from each end.
        //
        // - Window-only mean would be 40 → passes threshold 30.
        // - Whole-read mean is (5*5 + 40*10 + 5*5) / 20 = 22.5 → fails 30.
        //
        // We use min_len = 10 (absolute) so Δ_L = 10 and the 5+5 hard
        // trims fit exactly (rem_budget = 0). The hard-trim guards
        // therefore pass and the streaming walk runs, which is what
        // exercises the mean-q computation.
        let b = b"ACGTACGTACGTACGTACGT"; // L = 20
        let mut q = vec![5u8; b.len()];
        for i in 5..15 {
            q[i] = 40;
        }
        let mut params = Params::default();
        params.min_len = MinLen::Abs(10);
        params.mean_q = Some(30.0);
        params.trim_5 = Some(5);
        params.trim_3 = Some(5);

        // Whole-read mean = 22.5 < 30 → drop with DropMeanQ.
        // The streaming walk visits all 20 quality bytes (including the
        // low-Q flanks outside the polyx window) before checking mean-q.
        assert!(matches!(
            polyx_trim(b, &q, &params),
            TrimResult::Drop {
                reason: DropReason::MeanQ,
                ..
            }
        ));

        // With threshold 20: 22.5 ≥ 20 → keep. Polyx window content is
        // ACGTACGTAC (no qualifying run), so the only trim is the hard
        // cuts.
        params.mean_q = Some(20.0);
        match polyx_trim(b, &q, &params) {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 5, "5' hard trim only");
                assert_eq!(t3_poly, 5, "3' hard trim only");
                assert!(runs.is_empty());
            }
            r => panic!("expected Emit, got {:?}", r),
        }
    }

    #[test]
    fn trim_hard_trim_5_and_3() {
        // Hard trim N bases from either end, no dynamic trimming needed.
        let (b, q) = read(b"ACGTACGTACGTACGTACGTACGT", 40);
        let params = Params {
            trim_5: Some(3),
            trim_3: Some(2),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit {
                t5_poly,
                t3_poly,
                runs,
                ..
            } => {
                assert_eq!(t5_poly, 3);
                assert_eq!(t3_poly, 2);
                assert!(runs.is_empty());
            }
            r => panic!("expected Emit, got {:?}", r),
        }
    }

    #[test]
    fn trim_hard_trim_plus_scan() {
        // Hard trim consumes 5' adapter-like region; 5' scan then trims
        // any leftover polyX at the new 5' end.
        let mut b = vec![b'N'; 5]; // fake adapter
        b.extend(std::iter::repeat(b'A').take(10));
        b.extend_from_slice(b"CGTACGTACGTACGTACGTACGTACGTACGTACGTA");
        let q = vec![40u8; b.len()];
        let params = Params {
            trim_5: Some(5),
            ..Default::default()
        };
        match polyx_trim(&b, &q, &params) {
            TrimResult::Emit { t5_poly, .. } => {
                // Hard trim 5 + scan-trimmed polyA ≥ 10.
                assert!(t5_poly >= 15, "t5_poly = {}", t5_poly);
            }
            r => panic!("expected Emit, got {:?}", r),
        }
    }

    #[test]
    fn trim_hard_trim_consumes_whole_read_drops() {
        // trim_5 + trim_3 = L ⇒ sub_l = 0 < 2 ⇒ drop.
        let (b, q) = read(b"ACGTACGTACGT", 40);
        let params = Params {
            trim_5: Some(6),
            trim_3: Some(6),
            ..Default::default()
        };
        assert!(matches!(
            polyx_trim(&b, &q, &params),
            TrimResult::Drop {
                reason: DropReason::Hard,
                ..
            }
        ));
    }

    #[test]
    fn trim_hard_trim_exceeds_budget_drops() {
        // With min_len = 0.6L, L=10 ⇒ Δ_L = 4. trim_5=3 + trim_3=3 = 6 > 4.
        let (b, q) = read(b"ACGTACGTAC", 40);
        let params = Params {
            trim_5: Some(3),
            trim_3: Some(3),
            ..Default::default()
        };
        assert!(matches!(
            polyx_trim(&b, &q, &params),
            TrimResult::Drop {
                reason: DropReason::Hard,
                ..
            }
        ));
    }

    /// Conformance property from the paper (§11): `cascade(R,Q) =
    /// PolyXTrim(R,Q)`. With tier 1 enabled the walk is skipped for reads
    /// whose every w-window has < m_h - 1 matches; the result must be
    /// byte-identical to running the full Γ walk on the same read.
    ///
    /// Randomized over reads that exercise the boundary: mostly clean
    /// random sequence (the >99% case the prefilter exists for), plus
    /// reads with homopolymer runs of length m_h - 2 .. m_h + 4 at either
    /// end and in the interior, plus hard-trim and head-dropout corners.
    #[test]
    fn prefilter_conformance_cascade_equals_full_walk() {
        // Tiny deterministic xorshift64* PRNG — no external dep needed.
        struct Xs(u64);
        impl Xs {
            fn next(&mut self) -> u64 {
                let mut x = self.0;
                x ^= x >> 12;
                x ^= x << 25;
                x ^= x >> 27;
                self.0 = x;
                x.wrapping_mul(0x2545F4914F6CDD1D)
            }
            fn below(&mut self, n: usize) -> usize {
                (self.next() % n as u64) as usize
            }
            fn chance(&mut self, pct: u64) -> bool {
                self.next() % 100 < pct
            }
        }
        let mut rng = Xs(0x72706F6C7978); // "rpolyx"
        let bases = [b'A', b'C', b'G', b'T'];

        let mut param_sets: Vec<Params> = vec![Params::default()];
        // Parameter corners from the QC sweep + defaults.
        for &min_run in &[4usize, 6, 9] {
            for &mercy in &[1.5f64, 3.0, 3.5] {
                for &scan_interior in &[true, false] {
                    param_sets.push(Params {
                        min_run,
                        mercy,
                        scan_interior,
                        ..Default::default()
                    });
                }
            }
        }
        // Interior + mean-q corner: prefilter must still honor mean-q drops.
        param_sets.push(Params {
            scan_interior: true,
            mean_q: Some(25.0),
            ..Default::default()
        });
        // Head-dropout corner.
        param_sets.push(Params {
            head_dropout: Some(5),
            min_run: 6,
            ..Default::default()
        });

        let mut n_checked = 0usize;
        let mut n_rejected = 0usize;
        for params in &param_sets {
            for _ in 0..800 {
                let l = 2 + rng.below(178);
                let mut b: Vec<u8> = (0..l).map(|_| bases[rng.below(4)]).collect();
                // Sprinkle homopolymers of varied length at random positions.
                if rng.chance(35) {
                    let run_len = 2 + rng.below(10);
                    let start = rng.below(l.saturating_sub(run_len) + 1);
                    let base = bases[rng.below(4)];
                    for i in start..(start + run_len).min(l) {
                        b[i] = base;
                    }
                }
                let q: Vec<u8> = (0..l).map(|_| (2 + rng.below(39)) as u8).collect();

                // Count how often tier 1 would actually skip (engage check).
                let w = params
                    .min_run
                    .saturating_add(params.mercy.floor() as usize)
                    .saturating_sub(1);
                if w >= 2
                    && l >= w + 1
                    && (prefilter_end_rejects(&b, 0, l, w, params.min_run)
                        || prefilter_end_rejects(&b, l - w - 1, l, w, params.min_run))
                {
                    n_rejected += 1;
                }

                // Reference: full walk, prefilter disabled.
                FORCE_NO_PREFILTER.store(true, std::sync::atomic::Ordering::Relaxed);
                let full = polyx_trim(&b, &q, params);
                // Cascade: prefilter enabled (production path).
                FORCE_NO_PREFILTER.store(false, std::sync::atomic::Ordering::Relaxed);
                let casc = polyx_trim(&b, &q, params);

                match (&full, &casc) {
                    (
                        TrimResult::Drop {
                            reason: r1,
                            t5_poly: a,
                            t3_poly: c,
                            ..
                        },
                        TrimResult::Drop {
                            reason: r2,
                            t5_poly: d,
                            t3_poly: e,
                            ..
                        },
                    ) => {
                        assert_eq!(r1, r2, "drop reason mismatch");
                        assert_eq!(a, d, "t5_poly mismatch on drop");
                        assert_eq!(c, e, "t3_poly mismatch on drop");
                    }
                    (
                        TrimResult::Emit {
                            t5_poly: a,
                            t3_poly: c,
                            runs: ra,
                            ..
                        },
                        TrimResult::Emit {
                            t5_poly: d,
                            t3_poly: e,
                            runs: rb,
                            ..
                        },
                    ) => {
                        assert_eq!(a, d, "t5_poly mismatch");
                        assert_eq!(c, e, "t3_poly mismatch");
                        assert_eq!(ra, rb, "interior runs mismatch");
                    }
                    (f, c) => {
                        panic!("result kind mismatch: full={f:?} cascade={c:?}");
                    }
                }
                n_checked += 1;
            }
        }

        FORCE_NO_PREFILTER.store(false, std::sync::atomic::Ordering::Relaxed);
        assert!(n_checked > 0);
        // The prefilter must actually engage on this corpus — otherwise the
        // conformance test above would be comparing two identical walks.
        assert!(
            n_rejected > 0,
            "tier-1 prefilter never engaged on the randomized corpus"
        );
    }
}
