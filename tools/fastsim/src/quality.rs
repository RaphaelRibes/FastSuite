//! Quality model: per-position Phred curve, per-read baseline variance,
//! quality collapse inside G-runs, and base substitution errors.
//!
//! For each mate's read:
//! 1. Sample a per-read baseline `q_base` from `read_meanq[mate]` (the
//!    measured per-read mean-Q histogram). Offset δ = q_base − mean(pos_q).
//! 2. Per-position curve: q_i = pos_q[mate][i] + δ, clamped to [0, 41].
//! 3. Bases inside a G-run (head [0,h) or tail [L−t, L)) lose
//!    `tail_q_collapse[mate]` Phred (two-colour chemistry signature).
//! 4. Adapter bases keep the position curve (no collapse — adapter isn't a
//!    G-run).
//! 5. Substitution error per base with p = 10^(−q/10) (uniform over the
//!    other three bases). The emitted base reflects the error, so the FASTQ
//!    is self-consistent: a Q30 base has a 1/1000 chance of being wrong.
//! 6. Ambiguous-base (N) overlay, profile-driven: the instrument's no-call /
//!    padding behaviour, placed inside the biology window only (see
//!    `inject_ambiguous`). Absent from a profile -> nothing is drawn.
//!
//! Full-length (entirely artefact) reads use a flat `full_length_q`, no
//! position curve.

use rand::Rng;

use crate::artefact::Read;
use crate::profile::{Mate, Profile};
use crate::randutil::{sample_cdf, sample_gauss};

const MAX_Q: u8 = 41;
const MAX_Q_F: f64 = MAX_Q as f64;

/// Expected Phred lost to the hard cap at 41 when emitting q + z (z ~ N(0,sd)).
/// E[(q + z - 41)+] = sd·(φ(a) − a·Φ(−a)), a = (41−q)/sd.
fn cap_loss(q: f64, sd: f64) -> f64 {
    if q >= MAX_Q_F {
        return q - MAX_Q_F + sd * 0.3989; // always capped; residual ~ sd·φ(0)
    }
    let a = (MAX_Q_F - q) / sd;
    let phi = (-0.5 * a * a).exp() / 2.5066;
    let phi_neg = 0.5 * (1.0 - erf(a / std::f64::consts::SQRT_2)); // Φ(−a) = 1 − Φ(a)
    sd * (phi - a * phi_neg)
}

/// Standard normal CDF via the error function.
fn erf(x: f64) -> f64 {
    // Abramowitz-Stegun 7.1.26 (max error 1.5e-7) — plenty for cap loss.
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
        + 0.254829592)
        * t;
    sign * (1.0 - poly * (-x * x).exp())
}

/// One simulated mate: emitted (possibly substituted) sequence and Phreds.
#[derive(Debug, Clone)]
pub struct MateRead {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>, // Phred scores (0..=41), NOT ASCII
}

/// Assign quality + errors to one read.
pub fn qualify<R: Rng + ?Sized>(p: &Profile, mate: Mate, read: &Read, rng: &mut R) -> MateRead {
    let l = read.seq.len();
    let mut qual = vec![0u8; l];

    if read.full {
        // Entirely artefact: flat quality, errors at that quality.
        let q = p.full_length_q[Profile::mate_idx(mate)];
        qual.iter_mut().for_each(|v| *v = q);
        let seq = apply_errors(&read.seq, &qual, rng);
        return MateRead { seq, qual };
    }

    // 1. Per-read baseline.
    let q_base = sample_cdf(&p.read_meanq[Profile::mate_idx(mate)], rng);
    // Center on the mean of the SAME distribution q_base is drawn from:
    // E[offset] == 0 by construction, so per-position means reproduce the
    // measurement exactly ((1-p_i)*pos_q[i] + p_i*tail_q_mean == measured[i]
    // by the deconvolution in build_profile.py). Do NOT center on
    // mean(pos_q_uncond): the profiler floors per-read means into bins, so
    // mean(read_meanq) < mean(pos_q_uncond) by ~0.5 Phred (R1: 37.36 vs
    // 37.85) — that would bias every read low by the same amount.
    let offset = q_base as f64 - p.read_meanq_mean(mate);

    // 2-4. Position curve + run collapse.
    // NOTE: measured data: R1 heads are 1bp and normal (38.3), R2 heads are
    // long real polyG runs that collapse (34.5); 3' tails collapse on both
    // mates (R1 26.4, R2 37.2). pos_q is the deconvolved OUT-OF-RUN position
    // curve. Head/tail bases use the flat measured run means, EXCEPT the
    // FIRST head base (position h of the run), which is normal quality (R1's
    // 1bp heads prove it: 38.3 == outside). Per-base Gaussian noise is added
    // so the final integer rounding does not bias the mean.
    let tail_q = p.tail_q_mean(mate);
    let head_q = p.head_q_mean(mate);
    let pos_q = p.pos_q(mate);
    let h = read.h as usize;
    let tail_start = l - read.t as usize;
    for i in 0..l {
        let q = if i < h {
            // first head base is normal; deeper ones collapse
            if i == 0 {
                pos_q[i] + offset
            } else {
                head_q + offset
            }
        } else if i >= tail_start {
            tail_q + offset
        } else {
            pos_q[i] + offset
        };
        // Per-base Gaussian noise + cap compensation. Without noise, rounding
        // q biases the mean (38.689 -> 39 -> +0.31). With zero-mean noise,
        // E[round(q+z)] == q, BUT the hard cap at 41 truncates the noise for
        // high-quality bases (E[min(41, q+z)] < q). Compensate: emit q + the
        // expected cap loss, so E[round(clamp(q+loss+z, 41))] == q.
        let sd = p.qual_noise;
        let mut v = q;
        if sd > 0.0 {
            let z = sample_gauss(0.0, sd, -20.0, 20.0, rng);
            v += cap_loss(q, sd) + z;
        }
        qual[i] = v.round().clamp(0.0, MAX_Q_F) as u8;
    }

    // 5. Substitutions.
    let mut seq = apply_errors(&read.seq, &qual, rng);
    // 6. Ambiguous bases: the instrument's no-call / padding overlay.
    inject_ambiguous(p, read, &mut seq, rng);
    MateRead { seq, qual }
}

/// Overlay the profile's measured ambiguous-base (N) pattern on an emitted read.
///
/// The run is placed inside the BIOLOGY window `[h, h+bio_len)` only, never over
/// an injected 5′/3′ G-run. The two are different artefacts: the G-runs are the
/// two-colour dropout the trim ground truth's `[h, min(I, L−t))` contract is
/// built on, while the N overlay is the instrument failing to make a call. An N
/// inside the window therefore cannot change `h`, `t`, or the trim window, and
/// the GT records the injected count in its own columns (see truth.rs) rather
/// than letting the overlay silently rewrite the artefact truth.
///
/// The placement reproduces the measured shapes: a run touching either edge of
/// the window gives the `NNNN...` / `GNNNN...` padding, and a run strictly
/// inside gives the sparse no-calls. When the read has no biology window (a
/// full-length artefact) or the profile carries no model, nothing is drawn, so
/// the RNG stream of every profile without the field is unchanged.
///
/// Returns the number of N bases injected.
fn inject_ambiguous<R: Rng + ?Sized>(p: &Profile, read: &Read, seq: &mut [u8], rng: &mut R) -> u32 {
    let Some(n) = &p.n_run else {
        return 0;
    };
    if read.full || read.bio_len == 0 || n.read_rate <= 0.0 {
        return 0;
    }
    if rng.gen::<f64>() >= n.read_rate {
        return 0;
    }
    let bio_len = read.bio_len as usize;
    let len = (sample_cdf(&n.len, rng) as usize).clamp(1, bio_len);
    let lo = read.h as usize; // first biology base
    let hi = lo + bio_len; // one past the last biology base
    let u: f64 = rng.gen();
    let start = if u < n.leading {
        lo
    } else if u < n.leading + n.trailing {
        hi - len
    } else if bio_len > len + 1 {
        lo + 1 + rng.gen_range(0..(bio_len - len - 1))
    } else {
        lo
    };
    for b in &mut seq[start..start + len] {
        *b = b'N';
    }
    len as u32
}

/// Apply per-base substitution errors at probability 10^(−q/10).
/// N bases are left untouched (no base to substitute to — keep N).
fn apply_errors<R: Rng + ?Sized>(seq: &[u8], qual: &[u8], rng: &mut R) -> Vec<u8> {
    let mut out = seq.to_vec();
    for i in 0..seq.len() {
        let q = qual[i];
        if q == 0 {
            // Q0: error probability 1 — always substitute.
            if seq[i] != b'N' {
                out[i] = substitute(seq[i], rng);
            }
            continue;
        }
        let p = 10f64.powf(-(q as f64) / 10.0);
        let u: f64 = rng.gen();
        if u < p && seq[i] != b'N' {
            out[i] = substitute(seq[i], rng);
        }
    }
    out
}

/// Substitute to a random one of the three other bases.
fn substitute<R: Rng + ?Sized>(base: u8, rng: &mut R) -> u8 {
    let options = match base {
        b'A' => *b"CGT",
        b'C' => *b"AGT",
        b'G' => *b"ACT",
        b'T' => *b"ACG",
        other => return other, // N etc.
    };
    options[rng.gen_range(0..3)]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::artefact::{build_read, Artefact};
    use crate::profile::RawProfile;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(1234)
    }

    /// Profile: flat pos_q = 30 everywhere, read_meanq all at 30,
    /// tail-run mean 20 (R1) / 30 (R2), adapter none.
    fn q_profile(tail_q_r1: f64, tail_q_r2: f64) -> Profile {
        let q150 = vec!["30"; 150].join(",");
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "qual_noise": 0.0,
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 10.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "pos_q_uncond": {{"R1": [{q}], "R2": [{q}]}},
            "run_q_mean": {{
                "R1": {{"head": 30.0, "tail": {tq1}}},
                "R2": {{"head": 30.0, "tail": {tq2}}}
            }},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            tq1 = tail_q_r1,
            tq2 = tail_q_r2,
            q = q150
        );
        Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap()
    }

    /// The q_profile shape plus an N model, so the overlay can be exercised.
    fn n_profile(read_rate: f64, run_len: u32, leading: f64, trailing: f64) -> Profile {
        let q150 = vec!["30"; 150].join(",");
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "qual_noise": 0.0,
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}},
            "n_run": {{"read_rate": {rr}, "len": {{"{rl}": 100}}, "leading": {ld}, "trailing": {tr}}}
        }}"#,
            rr = read_rate,
            rl = run_len,
            ld = leading,
            tr = trailing,
            q = q150
        );
        Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap()
    }

    #[test]
    fn ambiguous_overlay_is_absent_without_the_profile_field() {
        // q_profile carries no n_run: not a single N may be emitted, which is
        // also why the RNG stream of every pre-existing profile is unchanged.
        let p = q_profile(10.0, 0.0);
        assert!(p.n_run.is_none());
        let mut rng = rng();
        let r = read_with_tail(0);
        for _ in 0..200 {
            let mr = qualify(&p, Mate::R1, &r, &mut rng);
            assert!(!mr.seq.contains(&b'N'), "no N without the model");
        }
    }

    #[test]
    fn ambiguous_overlay_rate_and_placement_follow_the_profile() {
        // every read carries a 5 bp run at the 3' edge of the biology window
        let p = n_profile(1.0, 5, 0.0, 1.0);
        let mut rng = rng();
        let r = read_with_tail(0); // h=0, t=0 -> window is the whole read
        for _ in 0..100 {
            let mr = qualify(&p, Mate::R1, &r, &mut rng);
            assert_eq!(mr.seq.iter().filter(|&&b| b == b'N').count(), 5);
            assert!(
                mr.seq[145..150].iter().all(|&b| b == b'N'),
                "a trailing run touches the 3' edge"
            );
        }
        // a rate of 0.5 fires on about half the reads
        let p = n_profile(0.5, 5, 0.0, 1.0);
        let mut with_n = 0;
        for _ in 0..2000 {
            let mr = qualify(&p, Mate::R1, &r, &mut rng);
            if mr.seq.contains(&b'N') {
                with_n += 1;
            }
        }
        let frac = with_n as f64 / 2000.0;
        assert!((frac - 0.5).abs() < 0.05, "read rate {frac}");
    }

    #[test]
    fn ambiguous_overlay_never_touches_the_artefact_runs() {
        // h=5 head, t=10 tail: a leading run must sit between them, so the
        // measured h/t (and the trim ground truth built on them) is unchanged
        let p = n_profile(1.0, 5, 1.0, 0.0);
        let mut rng = rng();
        let r = Read {
            seq: vec![b'A'; 150],
            h: 5,
            t: 10,
            bio_len: 135,
            adapter_len: 0,
            full: false,
        };
        let mr = qualify(&p, Mate::R1, &r, &mut rng);
        assert!(mr.seq[0..5].iter().all(|&b| b != b'N'), "head run untouched");
        assert!(
            mr.seq[140..150].iter().all(|&b| b != b'N'),
            "tail run untouched"
        );
        assert!(
            mr.seq[5..10].iter().all(|&b| b == b'N'),
            "the leading run starts at the first biology base"
        );
    }

    fn read_with_tail(t: u32) -> Read {
        Read {
            seq: vec![b'A'; 150],
            h: 0,
            t,
            bio_len: 150 - t,
            adapter_len: 0,
            full: false,
        }
    }

    #[test]
    fn tail_bases_use_flat_run_mean() {
        let p = q_profile(20.0, 30.0);
        let mut rng = rng();
        // tail of 10 Gs; R1 tail-run mean 20 -> tail bases Q20, others Q30
        let r = read_with_tail(10);
        let mr = qualify(&p, Mate::R1, &r, &mut rng);
        assert!(mr.qual[0..140].iter().all(|&q| q == 30));
        assert!(mr.qual[140..150].iter().all(|&q| q == 20));
        // R2 tail-run mean 30 -> all 30
        let mr2 = qualify(&p, Mate::R2, &r, &mut rng);
        assert!(mr2.qual.iter().all(|&q| q == 30));
    }

    #[test]
    fn baseline_offset_centers_on_unconditional_mean() {
        // Deconvolved out-of-run curve (pos_q) = 35 everywhere; 50% of reads
        // carry an 80-base tail (T = 15) covering positions 70..149, so the
        // measured UNCONDITIONAL mean is 35 for i<70 and 25 for i>=70
        // (mean 29.67). mean(read_meanq) == 29.5 ≈ mean(pos_q_uncond).
        // New centering (E[offset]=0): overall sim mean ~29.5, pos0 ~35,
        // pos100 ~25. Old centering on mean(pos_q)=35 would give overall
        // ~24.2 (pos0 ~29.5, pos100 ~19.5) — a ~5.3 Phred bias.
        let q35 = vec!["35"; 150].join(",");
        let q_uncond = {
            let mut v = vec!["35"; 70];
            v.extend(vec!["25"; 80]);
            v.join(",")
        };
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "qual_noise": 0.0,
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 50, "80": 50}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 20.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q35}], "R2": [{q35}]}},
            "pos_q_uncond": {{"R1": [{q_uncond}], "R2": [{q35}]}},
            "run_q_mean": {{
                "R1": {{"head": 30.0, "tail": 15.0}},
                "R2": {{"head": 30.0, "tail": 30.0}}
            }},
            "read_meanq": {{"R1": {{"35": 50, "24": 50}}, "R2": {{"35": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            q35 = q35,
            q_uncond = q_uncond
        );
        let p = Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap();
        let mut rng = rng();
        // 50% no-tail reads + 50% 80-base-tail reads.
        let mut total = 0f64;
        let mut pos0_total = 0f64;
        let mut pos100_total = 0f64;
        let n = 2000;
        for i in 0..n {
            let t = if i % 2 == 0 { 0 } else { 80 };
            let r = read_with_tail(t);
            let mr = qualify(&p, Mate::R1, &r, &mut rng);
            total += mr.qual.iter().map(|&q| q as f64).sum::<f64>();
            pos0_total += mr.qual[0] as f64;
            pos100_total += mr.qual[100] as f64;
        }
        let mean = total / (n as f64 * 150.0);
        // Correct centering (on read_meanq mean 29.5, what q_base is drawn
        // from, E[offset]=0): overall ~30.2, pos0 ~35.5, pos100 ~25.5.
        // Buggy centering (on pos_q_uncond mean 29.67): overall ~29.2,
        // pos0 ~34.5, pos100 ~24.5 — a ~1 Phred bias on every position.
        // Expected values include integer-Phred rounding.
        assert!(
            (mean - 30.2).abs() < 0.5,
            "overall sim mean {mean:.2} must be ~30.2 (read_meanq centering), not ~29.2 (pos_q_uncond centering)"
        );
        assert!(
            (pos0_total / n as f64 - 35.5).abs() < 0.5,
            "5' out-of-run bases must be ~35.5, got {:.2}",
            pos0_total / n as f64
        );
        assert!(
            (pos100_total / n as f64 - 25.5).abs() < 0.5,
            "3' tail-covered bases must be ~25.5 unconditional, got {:.2}",
            pos100_total / n as f64
        );
    }

    #[test]
    fn head_bases_use_flat_run_mean() {
        // Regression: heads were once excluded from collapse entirely, which
        // broke R2's 5' end (real R2 head runs collapse to 34.5). Head bases
        // must use the measured head-run mean, not the out-of-run curve.
        let mut rng = rng();
        let r = Read {
            seq: vec![b'G'; 150],
            h: 5,
            t: 0,
            bio_len: 145,
            adapter_len: 0,
            full: false,
        };
        let q30 = vec!["30"; 150].join(",");
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "qual_noise": 0.0,
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 10.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q30}], "R2": [{q30}]}},
            "pos_q_uncond": {{"R1": [{q30}], "R2": [{q30}]}},
            "run_q_mean": {{
                "R1": {{"head": 15.0, "tail": 20.0}},
                "R2": {{"head": 30.0, "tail": 30.0}}
            }},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            q30 = q30
        );
        let p = Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap();
        let mr = qualify(&p, Mate::R1, &r, &mut rng);
        // FIRST head base (position 0) is normal (pos_q 30); deeper head
        // bases (1..4) collapse to head_q_mean 15.
        assert_eq!(mr.qual[0], 30, "first head base is normal quality");
        assert!(
            mr.qual[1..5].iter().all(|&q| q == 15),
            "deep head bases use head_q_mean 15"
        );
        assert!(
            mr.qual[5..150].iter().all(|&q| q == 30),
            "body uses pos_q 30"
        );
    }

    #[test]
    fn substitution_rate_matches_q30() {
        // apply_errors directly on 1e6 bases at Q30 -> ~0.1% substitutions
        let seq = vec![b'A'; 1_000_000];
        let qual = vec![30u8; 1_000_000];
        let mut rng = rng();
        let out = apply_errors(&seq, &qual, &mut rng);
        let n: u64 = out.len() as u64;
        let mismatches = out.iter().zip(seq.iter()).filter(|(a, b)| a != b).count() as u64;
        let rate = mismatches as f64 / n as f64;
        assert!((rate - 0.001).abs() < 0.0005, "substitution rate {rate}");
        assert_eq!(out.len(), qual.len());
    }

    #[test]
    fn determinism_same_seed() {
        let p = q_profile(10.0, 0.0);
        let r = read_with_tail(10);
        let a = qualify(&p, Mate::R1, &r, &mut ChaCha8Rng::seed_from_u64(5));
        let b = qualify(&p, Mate::R1, &r, &mut ChaCha8Rng::seed_from_u64(5));
        assert_eq!(a.seq, b.seq);
        assert_eq!(a.qual, b.qual);
    }

    #[test]
    fn per_read_variance_bimodal() {
        // read_meanq bimodal: half the reads at Q20, half at Q40.
        let q150 = vec!["30"; 150].join(",");
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "qual_noise": 0.0,
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"20": 50, "40": 50}}, "R2": {{"20": 50, "40": 50}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            q = q150
        );
        let p = Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap();
        let mut rng = rng();
        let r = read_with_tail(0);
        let mut means = Vec::new();
        for _ in 0..2000 {
            let mr = qualify(&p, Mate::R1, &r, &mut rng);
            let m: f64 = mr.qual.iter().map(|&q| q as f64).sum::<f64>() / mr.qual.len() as f64;
            means.push(m);
        }
        // two clusters around 20 and 40
        let low = means.iter().filter(|&&m| (m - 20.0).abs() < 3.0).count();
        let high = means.iter().filter(|&&m| (m - 40.0).abs() < 3.0).count();
        assert!(low > 800, "low cluster {low}");
        assert!(high > 800, "high cluster {high}");
    }

    #[test]
    fn full_length_read_uses_flat_quality() {
        let p = q_profile(10.0, 0.0);
        let mut rng = rng();
        let r = Read {
            seq: vec![b'G'; 150],
            h: 0,
            t: 150,
            bio_len: 0,
            adapter_len: 0,
            full: true,
        };
        let mr = qualify(&p, Mate::R1, &r, &mut rng);
        assert!(mr.qual.iter().all(|&q| q == 20));
    }

    #[test]
    fn quality_round_trip_with_build_read() {
        // integration: build a real read, qualify it, check lengths + bounds
        let p = q_profile(5.0, 2.0);
        let mut rng = rng();
        let insert = vec![b'A'; 300];
        let a = Artefact {
            h1: 0,
            t1: 5,
            h2: 0,
            t2: 3,
            a1: 0,
            a2: 0,
            full: false,
        };
        let r1 = build_read(&p, &insert, 300, Mate::R1, &a, &mut rng);
        let r2 = build_read(&p, &insert, 300, Mate::R2, &a, &mut rng);
        let m1 = qualify(&p, Mate::R1, &r1, &mut rng);
        let m2 = qualify(&p, Mate::R2, &r2, &mut rng);
        assert_eq!(m1.seq.len(), 150);
        assert_eq!(m1.qual.len(), 150);
        assert_eq!(m2.seq.len(), 150);
        assert!(m1.qual.iter().all(|&q| q <= 41));
        assert!(m2.qual.iter().all(|&q| q <= 41));
        assert!(m1
            .seq
            .iter()
            .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')));
    }
}
