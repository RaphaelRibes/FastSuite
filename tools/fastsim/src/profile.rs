//! Typed profile: the JSON file that drives every simulation parameter.
//!
//! The profile is built from real measurements (see `scripts/build_profile.py`):
//! - `head` / `tail`: per-mate G-run length histograms `{"len": count}`
//! - `joint_tail`: R1×R2 tail-length counts `{"r1,r2": count}`, lengths capped
//! - `pos_q`: per-position mean Phred curves (R1/R2)
//! - `read_meanq`: per-read mean-Q histograms (for per-read variance)
//! - `tail_q_collapse`: Phred lost inside G-runs (two-colour chemistry)
//! - `full_length_*`: the rare entirely-artefact read class
//! - `adapter`: sequence + readthrough mode (e.g. Nextera)
//! - `geometry`: where the artefact sits (`tail3p_both`, `head5p_r2`, ...)
//!
//! Histograms are stored as counts (exactly like the profiler emits) and
//! converted to cumulative CDFs once at load time for O(log n) sampling.

use anyhow::{bail, Context, Result};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::Path;

pub const JOINT_CAP: u32 = 64; // tail lengths >= 64 collapse into the cap bin

fn default_qual_noise() -> f64 {
    0.5
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mate {
    R1,
    R2,
}

impl Mate {
    pub fn both() -> [Mate; 2] {
        [Mate::R1, Mate::R2]
    }
    pub fn as_str(&self) -> &'static str {
        match self {
            Mate::R1 => "R1",
            Mate::R2 => "R2",
        }
    }
}

/// A length histogram `{len: count}`, converted to a cumulative CDF
/// `[(len, cumulative_fraction)]` (monotonic, ends at 1.0).
#[derive(Debug, Clone)]
pub struct Cdf {
    pub bins: Vec<(u32, f64)>,
    /// Total count the CDF was built from (for diagnostics).
    pub total: u64,
}

#[derive(Debug, Clone)]
pub struct Profile {
    pub read_len: u32,
    /// Insert-length histogram `{len: count}`, measured by `fastcheck
    /// --insert-hist`. When present it replaces the mean/sd model, because real
    /// insert distributions are not Gaussian: a tagmentation library is broad
    /// and heavy-tailed, and the spread is what sets the readthrough rate.
    ///
    /// Lengths above the measurement's resolution go in ONE saturated bin (the
    /// profile builder uses 2*read_len+1). That is not a modelling choice: for a
    /// read shorter than the insert, every such length produces the same reads —
    /// no adapter, no overlap — so a bin is exactly as faithful as a fitted tail
    /// and does not invent a shape the data cannot show.
    pub insert: Option<Cdf>,
    pub gc: f64,
    pub geometry: String,
    /// per-base Gaussian noise sd added before integer rounding (Phred).
    /// Without it, rounding the per-position mean biases emitted bases.
    pub qual_noise: f64,
    /// head CDF per mate: 5′ G-run length (0 = read does not start with G)
    pub head: [Cdf; 2],
    /// tail CDF per mate: 3′ G-run length (0 = read does not end with G)
    pub tail: [Cdf; 2],
    /// joint R1×R2 tail CDF: cumulative over 65×65 grid, cap JOINT_CAP
    pub joint_tail: JointCdf,
    /// Phred lost for bases inside a G-run, per mate (derived, kept for compat)
    pub tail_q_collapse: [f64; 2],
    /// per-position mean Phred OUTSIDE G-runs, per mate (deconvolved).
    /// Kept FRACTIONAL: rounding the curve early biases every emitted base
    /// (e.g. 38.689 -> 39 -> +0.31 on average); only the final per-base
    /// quality is rounded, after the per-read offset.
    pub pos_q: [Vec<f64>; 2],
    /// per-position mean Phred unconditional (measured), per mate
    pub pos_q_uncond: [Vec<f64>; 2],
    /// flat mean Phred INSIDE 3' tail runs, per mate (measured run_q_mean.tail)
    pub tail_q_mean: [f64; 2],
    /// flat mean Phred INSIDE 5' head runs EXCLUDING the first base, per mate
    /// (deep head mean, `deep_head_q` in the profile; falls back to
    /// run_q_mean.head). The first head base is normal quality (a 1bp head is
    /// a biological G); only deeper bases collapse. R1: heads 1-2bp, Q ~38.3
    /// (no collapse); R2: long real polyG runs, deep Q ~32.
    pub head_q_mean: [f64; 2],
    /// per-read mean-Q histogram per mate (for baseline variance)
    pub read_meanq: [Cdf; 2],
    pub full_length_rate: f64,
    pub full_length_q: [u8; 2],
    pub adapter: Adapter,
}

#[derive(Debug, Clone)]
pub struct Adapter {
    pub seq: Vec<u8>,
    pub mode: AdapterMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterMode {
    None,
    Readthrough,
}

/// 65×65 cumulative distribution (R1 length × R2 length, capped).
#[derive(Debug, Clone)]
pub struct JointCdf {
    /// cumulative over the flat grid; index = r1 * 65 + r2
    pub cum: Vec<f64>,
}

impl JointCdf {
    pub const N: usize = (JOINT_CAP as usize) + 1;

    pub fn get(&self, r1: u32, r2: u32) -> f64 {
        self.cum[r1 as usize * Self::N + r2 as usize]
    }
}

impl Cdf {
    /// Build a CDF from a `{len: count}` histogram. Missing bins are 0.
    pub fn from_map(m: &BTreeMap<String, u64>) -> Result<Self> {
        let total: u64 = m.values().sum();
        if total == 0 {
            bail!("histogram has zero total count");
        }
        // Expand to a dense vec over [0, max_len].
        let max_len = m
            .keys()
            .map(|k| {
                k.parse::<u32>()
                    .with_context(|| format!("bad histogram key {k:?}"))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .max()
            .unwrap_or(0);
        let mut dense = vec![0u64; max_len as usize + 1];
        for (k, v) in m {
            let i = k
                .parse::<usize>()
                .with_context(|| format!("bad histogram key {k:?}"))?;
            dense[i] += v;
        }
        Ok(Self::from_dense(&dense))
    }

    pub fn from_dense(dense: &[u64]) -> Self {
        let total: u64 = dense.iter().sum();
        let mut cum = 0f64;
        let mut bins = Vec::with_capacity(dense.len());
        for (len, &c) in dense.iter().enumerate() {
            cum += c as f64 / total as f64;
            bins.push((len as u32, cum));
        }
        // Ensure last bin == 1.0 (float safety).
        if let Some(last) = bins.last_mut() {
            last.1 = 1.0;
        }
        Cdf { bins, total }
    }

    pub fn from_vec(v: &[u64]) -> Self {
        Self::from_dense(v)
    }
}

impl Profile {
    /// Insert length for the next pair. The profile's measured distribution wins
    /// when it has one; the CLI's Gaussian is the fallback for profiles fitted
    /// before the insert module existed.
    pub fn sample_insert<R: rand::Rng + ?Sized>(
        &self,
        mean: f64,
        sd: f64,
        floor: u32,
        rng: &mut R,
    ) -> u32 {
        match &self.insert {
            Some(cdf) => crate::randutil::sample_cdf(cdf, rng),
            None => crate::biology::sample_insert_len(mean, sd, floor, rng),
        }
    }

    /// Load and validate a profile JSON.
    pub fn from_json(path: &Path) -> Result<Self> {
        let raw: RawProfile = serde_json::from_slice(
            &std::fs::read(path).with_context(|| format!("read {}", path.display()))?,
        )
        .with_context(|| format!("parse {}", path.display()))?;
        Self::from_raw(raw)
    }

    pub fn from_raw(raw: RawProfile) -> Result<Self> {
        let head = [Cdf::from_map(&raw.head.R1)?, Cdf::from_map(&raw.head.R2)?];
        let tail = [Cdf::from_map(&raw.tail.R1)?, Cdf::from_map(&raw.tail.R2)?];

        // Joint: parse "r1,r2" -> count, cap at JOINT_CAP, build flat cumulative.
        let mut grid = vec![0f64; JointCdf::N * JointCdf::N];
        let mut jtotal = 0u64;
        for (k, v) in &raw.joint_tail {
            let (a, b) = k
                .split_once(',')
                .with_context(|| format!("bad joint key {k:?}"))?;
            let r1 = a.parse::<u32>()?.min(JOINT_CAP);
            let r2 = b.parse::<u32>()?.min(JOINT_CAP);
            grid[r1 as usize * JointCdf::N + r2 as usize] += *v as f64;
            jtotal += v;
        }
        if jtotal == 0 {
            bail!("joint_tail has zero total count");
        }
        let mut cum = Vec::with_capacity(grid.len());
        let mut acc = 0f64;
        for &c in &grid {
            acc += c / jtotal as f64;
            cum.push(acc);
        }
        let n = cum.len();
        cum[n - 1] = 1.0;

        let read_len = raw.read_len as usize;
        if raw.pos_q.R1.len() != read_len || raw.pos_q.R2.len() != read_len {
            bail!(
                "pos_q length {} must equal read_len {}",
                raw.pos_q.R1.len(),
                read_len
            );
        }
        // Measured per-position means are fractional; keep them fractional —
        // rounding the curve here would bias every emitted base.
        let pos_q = [raw.pos_q.R1.clone(), raw.pos_q.R2.clone()];

        // Unconditional (measured) per-position curve. Profiles built before
        // the deconvolution (build_profile.py) lack pos_q_uncond: fall back
        // to pos_q itself (then out-of-run == unconditional, correct for
        // reads without tails).
        let pos_q_uncond = match raw.pos_q_uncond {
            Some(u) => [u.R1.clone(), u.R2.clone()],
            None => pos_q.clone(),
        };

        // Flat in-tail mean quality (measured run_q_mean.tail). Fall back to
        // deriving it from the collapse: in_tail = mean(pos_q) - collapse.
        let mean_of = |v: &[f64]| -> f64 { v.iter().sum::<f64>() / v.len().max(1) as f64 };
        let (tail_q_mean, head_q_mean) = match &raw.run_q_mean {
            Some(rm) => (
                [
                    rm.R1
                        .tail
                        .unwrap_or_else(|| mean_of(&pos_q[0]) - raw.tail_q_collapse.R1),
                    rm.R2
                        .tail
                        .unwrap_or_else(|| mean_of(&pos_q[1]) - raw.tail_q_collapse.R2),
                ],
                [
                    raw.deep_head_q
                        .as_ref()
                        .map(|d| d.R1)
                        .or(rm.R1.head)
                        .unwrap_or_else(|| mean_of(&pos_q[0])),
                    raw.deep_head_q
                        .as_ref()
                        .map(|d| d.R2)
                        .or(rm.R2.head)
                        .unwrap_or_else(|| mean_of(&pos_q[1])),
                ],
            ),
            None => (
                [
                    mean_of(&pos_q[0]) - raw.tail_q_collapse.R1,
                    mean_of(&pos_q[1]) - raw.tail_q_collapse.R2,
                ],
                [mean_of(&pos_q[0]), mean_of(&pos_q[1])],
            ),
        };

        let mode = match raw.adapter.mode.as_deref() {
            None | Some("readthrough") => AdapterMode::Readthrough,
            Some("none") => AdapterMode::None,
            Some(other) => bail!("unknown adapter mode {other:?}"),
        };

        let insert = match &raw.insert {
            Some(m) => Some(Cdf::from_map(m)?),
            None => None,
        };

        Ok(Profile {
            read_len: raw.read_len,
            gc: raw.gc,
            geometry: raw.geometry,
            qual_noise: raw.qual_noise,
            insert,
            head,
            tail,
            joint_tail: JointCdf { cum },
            tail_q_collapse: [raw.tail_q_collapse.R1, raw.tail_q_collapse.R2],
            pos_q,
            pos_q_uncond,
            tail_q_mean,
            head_q_mean,
            read_meanq: [
                Cdf::from_map(&raw.read_meanq.R1)?,
                Cdf::from_map(&raw.read_meanq.R2)?,
            ],
            full_length_rate: raw.full_length_rate,
            full_length_q: [
                raw.full_length_q.R1.round().clamp(0.0, 41.0) as u8,
                raw.full_length_q.R2.round().clamp(0.0, 41.0) as u8,
            ],
            adapter: Adapter {
                seq: raw.adapter.seq.into_bytes(),
                mode,
            },
        })
    }

    pub fn mate_idx(mate: Mate) -> usize {
        match mate {
            Mate::R1 => 0,
            Mate::R2 => 1,
        }
    }

    pub fn pos_q(&self, mate: Mate) -> &[f64] {
        &self.pos_q[Self::mate_idx(mate)]
    }

    pub fn pos_q_uncond(&self, mate: Mate) -> &[f64] {
        &self.pos_q_uncond[Self::mate_idx(mate)]
    }

    pub fn tail_cdf(&self, mate: Mate) -> &Cdf {
        &self.tail[Self::mate_idx(mate)]
    }

    pub fn head_cdf(&self, mate: Mate) -> &Cdf {
        &self.head[Self::mate_idx(mate)]
    }

    pub fn tail_collapse(&self, mate: Mate) -> f64 {
        self.tail_q_collapse[Self::mate_idx(mate)]
    }

    pub fn tail_q_mean(&self, mate: Mate) -> f64 {
        self.tail_q_mean[Self::mate_idx(mate)]
    }

    pub fn head_q_mean(&self, mate: Mate) -> f64 {
        self.head_q_mean[Self::mate_idx(mate)]
    }

    /// Mean of the per-read mean-Q distribution (the CDF the baseline is
    /// sampled from). Centering the per-read offset on THIS mean gives
    /// E[offset] == 0 by construction. NOTE: this is NOT the same as
    /// mean(pos_q_uncond): the profiler floors per-read means into histogram
    /// bins (e.g. R1: 37.36 vs 37.85 — a 0.49 floor bias), so centering on
    /// pos_q_uncond would bias every read low.
    pub fn read_meanq_mean(&self, mate: Mate) -> f64 {
        let cdf = &self.read_meanq[Self::mate_idx(mate)];
        let mut prev = 0.0;
        let mut mean = 0.0;
        for (len, cum) in &cdf.bins {
            mean += *len as f64 * (cum - prev);
            prev = *cum;
        }
        mean
    }
}

// ---------------------------------------------------------------------------
// Raw JSON schema (serde). Histograms use string keys because JSON.
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
pub struct RawProfile {
    pub read_len: u32,
    pub gc: f64,
    pub geometry: String,
    #[serde(default = "default_qual_noise")]
    pub qual_noise: f64,
    pub head: RawMateHists,
    pub tail: RawMateHists,
    pub joint_tail: BTreeMap<String, u64>,
    pub tail_q_collapse: RawMateFloats,
    pub pos_q: RawMatePosQ,
    /// measured unconditional per-position curve (deconvolution input);
    /// absent in profiles built before the deconvolution
    #[serde(default)]
    pub pos_q_uncond: Option<RawMatePosQ>,
    /// measured mean Phred inside head/tail runs; absent in old profiles
    #[serde(default)]
    pub run_q_mean: Option<RawRunQMean>,
    /// deep head mean Phred (head runs minus their first base); absent in
    /// profiles built before the first-base-aware head model
    #[serde(default)]
    pub deep_head_q: Option<RawMateFloats>,
    pub read_meanq: RawMateHists,
    #[serde(default)]
    pub full_length_rate: f64,
    #[serde(default)]
    pub full_length_q: RawMateU8s,
    #[serde(default)]
    pub adapter: RawAdapter,
    /// measured insert-length histogram; absent in profiles fitted before the
    /// insert module existed, and then `--insert-mean/--insert-sd` are used
    #[serde(default)]
    pub insert: Option<BTreeMap<String, u64>>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)] // field names mirror the JSON keys (R1/R2)
pub struct RawRunQMean {
    pub R1: RawRunQMeanMate,
    pub R2: RawRunQMeanMate,
}

#[derive(Debug, Deserialize)]
pub struct RawRunQMeanMate {
    pub head: Option<f64>,
    pub tail: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)] // field names mirror the JSON keys (R1/R2)
pub struct RawMateHists {
    pub R1: BTreeMap<String, u64>,
    pub R2: BTreeMap<String, u64>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
pub struct RawMateFloats {
    pub R1: f64,
    pub R2: f64,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
pub struct RawMatePosQ {
    pub R1: Vec<f64>,
    pub R2: Vec<f64>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
pub struct RawMateU8s {
    #[serde(default)]
    pub R1: f64,
    #[serde(default)]
    pub R2: f64,
}

#[derive(Debug, Deserialize)]
pub struct RawAdapter {
    #[serde(default = "default_adapter_seq")]
    pub seq: String,
    #[serde(default)]
    pub mode: Option<String>,
}

fn default_adapter_seq() -> String {
    "CTGTCTCTTATACACATCT".to_string() // Nextera
}

impl Default for RawMateU8s {
    fn default() -> Self {
        RawMateU8s { R1: 20.0, R2: 20.0 }
    }
}

impl Default for RawAdapter {
    fn default() -> Self {
        RawAdapter {
            seq: default_adapter_seq(),
            mode: Some("readthrough".into()),
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_profile_json() -> String {
        // pos_q arrays: 150 elements of 38 (Phred)
        let q150 = vec!["38"; 150].join(",");
        format!(
            r#"{{
            "read_len": 150, "gc": 0.49, "geometry": "tail3p_both",
            "head": {{"R1": {{"0": 100, "1": 100}}, "R2": {{"0": 100, "1": 100}}}},
            "tail": {{"R1": {{"0": 60, "1": 30, "2": 10}}, "R2": {{"0": 60, "1": 30, "2": 10}}}},
            "joint_tail": {{"0,0": 36, "1,1": 12, "2,2": 2}},
            "tail_q_collapse": {{"R1": 11.7, "R2": 1.7}},
            "pos_q": {{"R1": [{q150}], "R2": [{q150}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}}
        }}"#
        )
    }

    fn test_profile() -> Profile {
        let raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        Profile::from_raw(raw).unwrap()
    }

    #[test]
    fn a_profile_insert_distribution_replaces_the_gaussian() {
        use rand::SeedableRng as _;
        // all mass at 40: every draw must be 40, whatever the CLI says
        let mut raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        raw.insert = Some(BTreeMap::from([("40".into(), 100u64)]));
        let p = Profile::from_raw(raw).unwrap();
        // a real RNG: StepRng's counter lands on the CDF's empty low bins at
        // r = 0, which is a property of the mock, not of the sampler
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        for _ in 0..200 {
            assert_eq!(p.sample_insert(300.0, 50.0, 150, &mut rng), 40);
        }

        // and with no distribution it falls back to the Gaussian around the mean
        let raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        let p = Profile::from_raw(raw).unwrap();
        assert!(p.insert.is_none());
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let l = p.sample_insert(300.0, 50.0, 150, &mut rng);
        assert!((150..=450).contains(&l), "gaussian draw {l} out of range");
    }

    #[test]
    fn cdf_from_map_total_and_ends_at_one() {
        let c = Cdf::from_map(&BTreeMap::from([
            ("0".into(), 60u64),
            ("1".into(), 30),
            ("2".into(), 10),
        ]))
        .unwrap();
        assert_eq!(c.total, 100);
        assert_eq!(c.bins.len(), 3);
        assert!((c.bins.last().unwrap().1 - 1.0).abs() < 1e-12);
        // cumulative at len 1 = 0.9
        assert!((c.bins[1].1 - 0.9).abs() < 1e-12);
    }

    #[test]
    fn cdf_zero_total_is_error() {
        let err = Cdf::from_map(&BTreeMap::new()).unwrap_err();
        assert!(err.to_string().contains("zero total"));
    }

    #[test]
    fn profile_parses_all_fields() {
        let p = test_profile();
        assert_eq!(p.read_len, 150);
        assert_eq!(p.geometry, "tail3p_both");
        assert_eq!(p.tail_cdf(Mate::R1).total, 100);
        assert_eq!(p.pos_q(Mate::R1).len(), 150);
        assert_eq!(p.tail_collapse(Mate::R1), 11.7);
        assert_eq!(p.adapter.seq, b"CTGTCTCTTATACACATCT");
        assert_eq!(p.adapter.mode, AdapterMode::Readthrough);
    }

    #[test]
    fn joint_cdf_marginals_match() {
        let p = test_profile();
        // P(R1 tail == 0) = (36 + 0) / 50 = 0.72
        // joint table: {0,0:36, 1,1:12, 2,2:2} -> total 50
        assert!((p.joint_tail.get(0, 0) - 36.0 / 50.0).abs() < 1e-12);
        // cumulative after {0,0} is 0.72; the flat order is row-major (r1 outer)
        let cum_00 = p.joint_tail.cum[0]; // index r1=0,r2=0
        assert!((cum_00 - 0.72).abs() < 1e-12);
        // after {1,1} (flat index N + 1) cumulative = 0.72+0.24 = 0.96
        let cum_11 = p.joint_tail.cum[JointCdf::N + 1];
        assert!((cum_11 - 0.96).abs() < 1e-12);
    }

    #[test]
    fn missing_profile_file_is_error() {
        let err = Profile::from_json(Path::new("/nonexistent/profile.json")).unwrap_err();
        assert!(err.to_string().contains("/nonexistent"));
    }

    #[test]
    fn bad_pos_q_length_is_error() {
        let mut raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        raw.pos_q.R1 = vec![38.0; 10];
        let err = Profile::from_raw(raw).unwrap_err();
        assert!(err.to_string().contains("pos_q"));
    }

    #[test]
    fn fractional_pos_q_accepted_and_kept_fractional() {
        // regression: build_profile.py emits fractional means (e.g. 38.784);
        // the curve must stay fractional so the final per-base rounding
        // (after the per-read offset) does not bias the emitted qualities.
        let mut raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        let mut q1 = vec![38.0; 150];
        q1[0] = 38.6;
        q1[1] = 30.4;
        q1[2] = 41.2;
        q1[3] = 55.0;
        q1[4] = -3.0;
        raw.pos_q.R1 = q1;
        raw.pos_q.R2 = raw.pos_q.R1.clone();
        raw.full_length_q.R1 = 26.431;
        raw.full_length_q.R2 = 36.93;
        let p = Profile::from_raw(raw).unwrap();
        assert_eq!(p.pos_q(Mate::R1)[0], 38.6); // kept fractional
        assert_eq!(p.pos_q(Mate::R1)[1], 30.4);
        // full_length_q is a single scalar and IS rounded at load
        assert_eq!(p.full_length_q[0], 26); // 26.431 rounds to 26
    }
}
