//! PolyG artefact injection: heads, tails, the full-length spike, and
//! adapter readthrough.
//!
//! Everything is driven by the profile's CDFs — the `geometry` string is
//! provenance metadata, NOT branching logic. A profile with mass in R2's
//! head CDF (the old paper's 5′-R2 regime) produces 5′ R2 artefacts; the
//! mandril profile (mass in both tail CDFs) produces 3′ tails on both
//! mates. Same code path, no hardcoded geometry.
//!
//! Read layout (5′ → 3′): `[head Gs] [biology] [adapter?] [tail Gs]`.
//! The adapter appears only when the insert is shorter than the read
//! window (readthrough). Adapter bases are NOT biology — the GT window
//! `[h, min(I, L−t))` already excludes them.

use rand::Rng;

use crate::biology::revcomp;
use crate::profile::{Mate, Profile};
use crate::randutil::{sample_cdf, sample_joint};

/// Per-pair artefact decision (one per mate: head/tail G-run lengths).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Artefact {
    pub h1: u32,
    pub t1: u32,
    pub h2: u32,
    pub t2: u32,
    /// adapter readthrough length per mate (0 = no adapter on that read).
    /// Filled in by the caller after build_read; the GT writer reports it so
    /// "where is the adapter" is explicit, not derived.
    pub a1: u32,
    pub a2: u32,
    /// Entirely-artefact pair (full-length spike or biology window < 1 bp):
    /// the correct answer for a trimmer is to DROP the pair.
    pub full: bool,
}

/// A constructed read: sequence + the boundaries needed by the quality
/// model and the ground-truth writer.
#[derive(Debug, Clone)]
pub struct Read {
    pub seq: Vec<u8>,
    /// 5′ G-run length
    pub h: u32,
    /// 3′ G-run length
    pub t: u32,
    /// Number of biology bases in the read (from the insert)
    pub bio_len: u32,
    /// Number of adapter bases (readthrough), 0 if insert >= window
    pub adapter_len: u32,
    /// Entirely artefact -> drop
    pub full: bool,
}

impl Read {
    /// Biology window in read coordinates: [h, h + bio_len).
    pub fn biology_window(&self) -> (u32, u32) {
        (self.h, self.h + self.bio_len)
    }
}

/// Sample one pair's artefact decision from the profile.
pub fn decide<R: Rng + ?Sized>(p: &Profile, rng: &mut R) -> Artefact {
    let l = p.read_len;

    // NOTE: full-length (entirely-polyG) reads are NOT sampled here. They
    // appear in the profile's head AND tail histograms (an all-G read has a
    // 150bp head and a 150bp tail), so the marginal CDFs already contain
    // them, and the near-full check below (h+t >= L) turns them into the
    // all-G class. A separate pair-level spike would DOUBLE them and, worse,
    // fabricate full-length R1 reads that the measured data does not have
    // (mandril: full_length_rate R1 = 0.00000, R2 = 0.00217 — R2 only).

    // 1. Correlated tails from the joint distribution, then marginal heads.
    let (t1, t2) = sample_joint(&p.joint_tail, rng);
    let h1 = sample_cdf(p.head_cdf(Mate::R1), rng);
    let h2 = sample_cdf(p.head_cdf(Mate::R2), rng);

    // 2. Near-full: biology window collapsed to < 1 bp -> drop class.
    let full = h1 + t1 >= l || h2 + t2 >= l;
    Artefact {
            h1,
            t1,
            h2,
            t2,
            a1: 0,
            a2: 0,
            full,
        }
    }

/// Build one mate's read from the insert.
///
/// `insert` is the top strand of the fragment; R1 takes its 5′ end, R2 the
/// reverse complement of its 3′ end. The read window `L − h − t` is filled
/// with biology first, then adapter (if the insert is shorter), then the
/// tail Gs; head Gs go at the 5′ end.
///
/// Boundary rule (important for distribution fidelity): the profile's tail
/// histogram is the TOTAL measured 3′ G-run (artefact + any biological Gs
/// at the boundary). To reproduce it exactly, the base immediately before
/// the injected tail (and after the head) must be non-G; if biology ends/
/// starts with G, that single base is replaced with a random non-G. Without
/// this, simulated runs would be `t + boundary_Gs`, over-convolving the
/// distribution (measured 51% vs 35% in the 20M-pair closed loop).
pub fn build_read<R: Rng + ?Sized>(
    p: &Profile,
    insert: &[u8],
    insert_len: u32,
    mate: Mate,
    a: &Artefact,
    rng: &mut R,
) -> Read {
    let l = p.read_len as usize;
    let (h, t) = match mate {
        Mate::R1 => (a.h1, a.t1),
        Mate::R2 => (a.h2, a.t2),
    };

    // Degenerate: biology window < 1 bp -> the read is all artefact.
    // NOTE: full-length is PER MATE (measured R1 rate 0.00000, R2 0.00217):
    // a pair where only R2 is full must not turn R1 into an all-G read.
    if (h + t) >= p.read_len {
        return Read {
            seq: vec![b'G'; l],
            h,
            t,
            bio_len: 0,
            adapter_len: 0,
            full: true,
        };
    }

    let window = p.read_len - h - t; // bases available for biology+adapter
    let bio = window.min(insert_len);
    let adapter_len = window - bio;

    // Biology slice of the insert, oriented for this mate.
    let mut seq: Vec<u8> = Vec::with_capacity(l);
    seq.extend(std::iter::repeat_n(b'G', h as usize));
    match mate {
        Mate::R1 => seq.extend_from_slice(&insert[0..bio as usize]),
        Mate::R2 => {
            let start = insert_len - bio;
            let far_end = &insert[start as usize..insert_len as usize];
            seq.extend_from_slice(&revcomp(far_end));
        }
    }
    // Adapter readthrough: fill the rest of the window with adapter bases
    // (cycled if the adapter is shorter than the gap).
    if adapter_len > 0 && p.adapter.mode == crate::profile::AdapterMode::Readthrough {
        let ad = &p.adapter.seq;
        for i in 0..adapter_len as usize {
            seq.push(ad[i % ad.len()]);
        }
    } else if adapter_len > 0 {
        // mode None: pad with Ns so the read still has the right length
        seq.extend(std::iter::repeat_n(b'N', adapter_len as usize));
    }
    seq.extend(std::iter::repeat_n(b'G', t as usize));

    // Boundary rule: make the measured G-runs exactly h and t.
    // - 5′: base at index h (first biology/adapter base) must be non-G.
    // - 3′: base at index L−t−1 (last base before the tail) must be non-G.
    // This includes t=0 / h=0: the measured run must be exactly the sampled
    // length, so a read with t=0 must not end in G either.
    // Replace with a random non-G when biology happens to start/end in G.
    let non_g = *b"ACT";
    // 5': base at index h (first base after the head; index 0 when h=0)
    // must be non-G, so a read sampled with h=0 never starts with a
    // biological G (otherwise the measured head rate inflates by ~GC/2).
    if (h as usize) < l {
        let i = h as usize;
        if seq[i] == b'G' {
            seq[i] = non_g[rng.gen_range(0..3)];
        }
    }
    if t < p.read_len {
        let i = l - t as usize - 1;
        if seq[i] == b'G' {
            seq[i] = non_g[rng.gen_range(0..3)];
        }
    }

    debug_assert_eq!(seq.len(), l, "read length mismatch");
    Read {
        seq,
        h,
        t,
        bio_len: bio,
        adapter_len,
        full: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profile::{AdapterMode, RawProfile, JOINT_CAP};
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(99)
    }

    fn q150() -> Vec<u8> {
        vec![38; 150]
    }

    /// Profile with tails only from joint (0:0.72, 1:0.24, 2:0.04 diagonal),
    /// heads all 0, no full-length spike, adapter off.
    fn tail_profile() -> Profile {
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 60, "1": 30, "2": 10}}, "R2": {{"0": 60, "1": 30, "2": 10}}}},
            "joint_tail": {{"0,0": 36, "1,1": 12, "2,2": 2}},
            "tail_q_collapse": {{"R1": 11.7, "R2": 1.7}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            q = q150()
                .iter()
                .map(|b| b.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap()
    }

    fn fixed_joint_profile(pairs: &[(u32, u32, u64)]) -> Profile {
        let mut joint = std::collections::BTreeMap::new();
        for (a, b, c) in pairs {
            joint.insert(format!("{a},{b}"), *c);
        }
        let q = q150();
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {j},
            "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            j = serde_json::to_string(&joint).unwrap(),
            q = q
                .iter()
                .map(|b| b.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap()
    }

    #[test]
    fn joint_tail_5_5_produces_exact_tails() {
        // all mass at (5,5): every pair gets t1=t2=5
        let p = fixed_joint_profile(&[(5, 5, 100)]);
        let mut rng = rng();
        for _ in 0..200 {
            let a = decide(&p, &mut rng);
            assert_eq!((a.t1, a.t2), (5, 5), "tail pair must be (5,5)");
            assert!(!a.full);
        }
        // build a read and check the tail is exactly 5 Gs preceded by non-G
        let insert = vec![b'A'; 300];
        let a = Artefact {
            h1: 0,
            t1: 5,
            h2: 0,
            t2: 5,
            a1: 0,
            a2: 0,
            full: false,
        };
        let r1 = build_read(&p, &insert, 300, Mate::R1, &a, &mut rng);
        assert_eq!(r1.seq.len(), 150);
        assert!(r1.seq[144] != b'G', "base before tail must be non-G");
        assert!(r1.seq[145..150].iter().all(|&b| b == b'G'));
        assert_eq!(r1.bio_len, 145);
        assert_eq!(r1.adapter_len, 0);
    }

    #[test]
    fn full_length_comes_from_near_full_path_not_spike() {
        // Regression: decide() used to sample a pair-level full-length spike
        // that set BOTH mates all-G. Real data has full-length only on R2
        // (mandril R1 rate 0.00000), so the spike fabricated R1 all-G reads.
        // Full-length must emerge from the near-full check (h+t >= L).
        // Profile with heads up to 150 on R2 only.
        let q150 = vec!["30"; 150].join(",");
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "tail3p_both",
            "qual_noise": 0.0,
            "head": {{
                "R1": {{"0": 90, "1": 10}},
                "R2": {{"0": 90, "1": 5, "150": 5}}
            }},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "pos_q_uncond": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            q = q150
        );
        let p = Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap();
        let mut rng = rng();
        let insert = vec![b'A'; 300];
        let mut r1_full = 0;
        let mut r2_full = 0;
        for _ in 0..2000 {
            let a = decide(&p, &mut rng);
            let r1 = build_read(&p, &insert, 300, Mate::R1, &a, &mut rng);
            let r2 = build_read(&p, &insert, 300, Mate::R2, &a, &mut rng);
            if r1.full {
                r1_full += 1;
            }
            if r2.full {
                r2_full += 1;
            }
        }
        // R1 head max is 1 -> never full; R2 has 5% heads of length 150.
        assert_eq!(r1_full, 0, "R1 must never be full-length (no long heads)");
        assert!(
            r2_full > 50 && r2_full < 250,
            "R2 full-length rate ~5% (got {r2_full}/2000)"
        );
    }

    #[test]
    fn zero_head_read_never_starts_with_g() {
        // Regression: the 5' boundary rule used to skip h=0 reads, so a read
        // sampled with h=0 could start with a biological G, inflating the
        // measured head rate by ~GC/2. The base at index 0 must be non-G.
        let p = fixed_joint_profile(&[(0, 0, 100)]);
        let mut rng = rng();
        let insert: Vec<u8> = vec![b'G'; 300]; // biology starts with G
        let a = decide(&p, &mut rng);
        assert_eq!((a.h1, a.t1), (0, 0));
        let r1 = build_read(&p, &insert, 300, Mate::R1, &a, &mut rng);
        assert!(r1.seq[0] != b'G', "h=0 read must not start with G");
    }

    #[test]
    fn no_artefact_read_equals_insert_slice() {
        // head=0, tail=0 -> R1 must equal insert[0..150]
        let p = fixed_joint_profile(&[(0, 0, 100)]);
        let mut rng = rng();
        let insert: Vec<u8> = (0..300).map(|i| b"ACGT"[i % 4]).collect();
        let a = decide(&p, &mut rng);
        assert_eq!((a.t1, a.t2), (0, 0));
        let r1 = build_read(&p, &insert, 300, Mate::R1, &a, &mut rng);
        assert_eq!(r1.seq, insert[0..150]);
        assert_eq!(r1.bio_len, 150);
        // R2 = revcomp(insert[150..300])
        let r2 = build_read(&p, &insert, 300, Mate::R2, &a, &mut rng);
        assert_eq!(r2.seq, revcomp(&insert[150..300]));
    }

    #[test]
    fn distributional_marginals_match() {
        let p = tail_profile();
        let mut rng = rng();
        let n = 200_000;
        let mut t1_count = [0u64; 3];
        let mut t2_count = [0u64; 3];
        for _ in 0..n {
            let a = decide(&p, &mut rng);
            t1_count[a.t1 as usize] += 1;
            t2_count[a.t2 as usize] += 1;
        }
        // joint diagonal P(0)=0.72, P(1)=0.24, P(2)=0.04
        for (i, pv) in [0.72, 0.24, 0.04].iter().enumerate() {
            let f1 = t1_count[i] as f64 / n as f64;
            assert!((f1 - pv).abs() < 0.005, "R1 bin {i}: {f1}");
            let f2 = t2_count[i] as f64 / n as f64;
            assert!((f2 - pv).abs() < 0.005, "R2 bin {i}: {f2}");
        }
    }

    #[test]
    fn near_full_marks_drop() {
        // h+t >= L -> full
        let p = tail_profile();
        let a = Artefact {
            h1: 100,
            t1: 60,
            h2: 0,
            t2: 0,
            a1: 0,
            a2: 0,
            full: false,
        };
        let mut rng = rng();
        let insert = vec![b'A'; 300];
        let r1 = build_read(&p, &insert, 300, Mate::R1, &a, &mut rng);
        assert!(r1.full);
        assert_eq!(r1.bio_len, 0);
        assert_eq!(r1.seq.len(), 150);
    }

    #[test]
    fn adapter_readthrough_fills_short_insert() {
        // insert 120 < read 150, no head/tail -> positions 120..150 = adapter
        let p = fixed_joint_profile(&[(0, 0, 100)]);
        // adapter mode on: make a profile with readthrough
        let mut p = p;
        p.adapter.mode = AdapterMode::Readthrough;
        p.adapter.seq = b"CTGTCTCTTATACACATCT".to_vec();
        let mut rng = rng();
        let insert: Vec<u8> = (0..120).map(|i| b"ACGT"[i % 4]).collect();
        let a = decide(&p, &mut rng);
        let r1 = build_read(&p, &insert, 120, Mate::R1, &a, &mut rng);
        assert_eq!(r1.seq.len(), 150);
        assert_eq!(r1.bio_len, 120);
        assert_eq!(r1.adapter_len, 30);
        let ad = &p.adapter.seq;
        for (i, &b) in r1.seq[120..150].iter().enumerate() {
            assert_eq!(b, ad[i % ad.len()], "adapter base at 120+{i}");
        }
    }

    #[test]
    fn head5p_r2_geometry_from_profile() {
        // old paper regime: mass in R2's head CDF only -> 5′ R2 artefacts
        let q = q150();
        let json = format!(
            r#"{{
            "read_len": 150, "gc": 0.5, "geometry": "head5p_r2",
            "head": {{"R1": {{"0": 100}}, "R2": {{"0": 60, "20": 40}}}},
            "tail": {{"R1": {{"0": 100}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}}
        }}"#,
            q = q
                .iter()
                .map(|b| b.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        let p = Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap();
        let mut rng = rng();
        let insert = vec![b'A'; 300];
        let mut n20 = 0u64;
        let n = 10_000;
        for _ in 0..n {
            let a = decide(&p, &mut rng);
            assert_eq!(a.t2, 0, "no tail in the 5′-R2 regime");
            if a.h2 == 20 {
                n20 += 1;
            }
            let r2 = build_read(&p, &insert, 300, Mate::R2, &a, &mut rng);
            // 20-G head at the 5′ end of R2
            if a.h2 == 20 {
                assert!(r2.seq[0..20].iter().all(|&b| b == b'G'));
                assert!(r2.seq[20] != b'G');
            }
        }
        let frac = n20 as f64 / n as f64;
        assert!((frac - 0.4).abs() < 0.02, "h2=20 rate {frac}");
    }

    #[test]
    fn joint_cap_bins_collapse() {
        // tails >= 64 collapse into the cap bin
        let p = fixed_joint_profile(&[(JOINT_CAP, JOINT_CAP, 100)]);
        let mut rng = rng();
        for _ in 0..50 {
            let a = decide(&p, &mut rng);
            assert_eq!(a.t1, JOINT_CAP);
            assert_eq!(a.t2, JOINT_CAP);
        }
    }

    #[test]
    fn boundary_rule_keeps_measured_run_exact() {
        // regression: biology ending in G must not extend the injected tail
        let p = fixed_joint_profile(&[(0, 0, 100)]);
        let mut rng = rng();
        // insert ends in a long G-run (biology boundary is G)
        let insert: Vec<u8> = [vec![b'A'; 140], vec![b'G'; 60]].concat();
        for _ in 0..500 {
            let a = decide(&p, &mut rng);
            assert_eq!((a.t1, a.t2), (0, 0));
            let r1 = build_read(&p, &insert, 200, Mate::R1, &a, &mut rng);
            // no tail injected -> measured run must be 0 (boundary replaced)
            let run = r1.seq.iter().rev().take_while(|&&b| b == b'G').count();
            assert_eq!(run, 0, "read must not end in G when t=0");
        }
        // with an injected tail t=5, measured run is exactly 5
        let a = Artefact {
            h1: 0,
            t1: 5,
            h2: 0,
            t2: 5,
            a1: 0,
            a2: 0,
            full: false,
        };
        let r1 = build_read(&p, &insert, 200, Mate::R1, &a, &mut rng);
        let run = r1.seq.iter().rev().take_while(|&&b| b == b'G').count();
        assert_eq!(run, 5, "measured run must equal injected tail exactly");
    }
}
