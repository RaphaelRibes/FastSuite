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
    /// Measured duplication level from the fitter's exact fingerprint counter
    /// (`duplication_level`). When present it is the duplicate-class rate; the
    /// CLI's --dup-rate is only the fallback for profiles without it.
    pub dup_rate: Option<f64>,
    pub geometry: String,
    /// per-position [a,c,g,t] composition at the fragment 5′ end (from
    /// fastcheck `base_position`). When non-empty, the fragment's leading bases
    /// are sampled from this instead of uniform GC, reproducing *located* base
    /// biases (a 5′ dropout that starts mid-read, an enriched position) that a
    /// G-run-length head model cannot express.
    pub head_bias: Vec<[f32; 4]>,
    /// per-position [a,c,g,t] composition at the fragment 3′ end, in the
    /// fragment's own 5′→3′ frame (R2's measured head, reverse-complemented).
    pub tail_bias: Vec<[f32; 4]>,
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
    /// flat mean Phred INSIDE 3' tail runs, per mate (measured run_q_mean.tail)
    pub tail_q_mean: [f64; 2],
    /// E[q_base] per mate: the mean of `read_meanq`, precomputed at load.
    ///
    /// `qualify` needs it once per read to centre the per-read offset, and it
    /// is a pure function of the profile — recomputing the CDF sum per read
    /// was ~40 bins of work for a constant (FS-4). The value is byte-identical
    /// to the per-read recomputation: same loop, same order, stored once.
    pub read_meanq_mean: [f64; 2],
    /// flat mean Phred INSIDE 5' head runs EXCLUDING the first base, per mate
    /// (deep head mean, `deep_head_q` in the profile; falls back to
    /// run_q_mean.head). The first head base is normal quality (a 1bp head is
    /// a biological G); only deeper bases collapse. R1: heads 1-2bp, Q ~38.3
    /// (no collapse); R2: long real polyG runs, deep Q ~32.
    pub head_q_mean: [f64; 2],
    /// per-read mean-Q histogram per mate (for baseline variance)
    pub read_meanq: [Cdf; 2],
    pub full_length_q: [u8; 2],
    pub adapter: Adapter,
    /// measured ambiguous-base (N) emission; absent -> emit none (see NModel)
    pub n_run: Option<NModel>,
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

/// Ambiguous-base (N) model, measured from the library's own output.
///
/// Real Illumina output carries `N` wherever the instrument could not make a
/// base call. Two regimes coexist and their mix is strongly library-dependent:
/// isolated no-calls sprinkled through the read, and — in a library with a
/// failed-cluster / padding share — a long poly-N run at a read end
/// (`NNNN...`, `GNNNN...`) that is often the single most-duplicated sequence in
/// the file. `random_dna_bias` drew only ACGT, so the simulator emitted none of
/// them and the padding-driven duplicates had no counterpart in the simulation.
///
/// * `read_rate`   — measured share of reads carrying at least one `N`.
/// * `len`         — measured run-length distribution (a run's length).
/// * `leading`     — share of runs touching the read's first base.
/// * `trailing`    — share of runs touching the read's last base.
///                   The remainder sit interior.
///
/// When the field is absent (a profile fitted before the model existed, or the
/// hand-authored test profiles) the simulator draws NOTHING extra, so the RNG
/// stream of every existing profile is unchanged.
#[derive(Debug, Clone)]
pub struct NModel {
    pub read_rate: f64,
    pub len: Cdf,
    pub leading: f64,
    pub trailing: f64,
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
}


/// P(read position p is covered by an injected 5' G-run), for p in 0..read_len.
///
/// `m` is the measured histogram of run LENGTHS {len: count}. Position p is
/// inside the run exactly when the run length h satisfies h > p, so
/// P(p covered) = (sum of counts for h > p) / (total counts).
fn head_cover(m: &BTreeMap<String, u64>, read_len: usize) -> Vec<f32> {
    let total: u64 = m.values().sum();
    if total == 0 {
        return vec![0.0; read_len];
    }
    let maxlen = m.keys().filter_map(|k| k.parse::<usize>().ok()).max().unwrap_or(0);
    let mut counts = vec![0u64; maxlen + 2];
    for (k, v) in m.iter() {
        if let Ok(len) = k.parse::<usize>() {
            counts[len] += *v;
        }
    }
    // greater[p] = number of reads with h > p
    let mut greater = vec![0u64; maxlen + 2];
    let mut acc = 0u64;
    for p in (0..=maxlen + 1).rev() {
        if p <= maxlen {
            greater[p] = acc;
        }
        if p <= maxlen {
            acc += counts[p];
        }
    }
    (0..read_len)
        .map(|p| if p <= maxlen { greater[p] as f32 / total as f32 } else { 0.0 })
        .collect()
}

/// P(read position p is covered by an injected 3' G-run), for p in 0..read_len.
///
/// A run of length t fills the LAST t positions of the read (read coordinates
/// L-t..L-1), so it covers read position p exactly when t >= L - p.
fn tail_cover(m: &BTreeMap<String, u64>, read_len: usize) -> Vec<f32> {
    let total: u64 = m.values().sum();
    if total == 0 {
        return vec![0.0; read_len];
    }
    let maxlen = m.keys().filter_map(|k| k.parse::<usize>().ok()).max().unwrap_or(0);
    let mut counts = vec![0u64; maxlen + 1];
    for (k, v) in m.iter() {
        if let Ok(len) = k.parse::<usize>() {
            counts[len] += *v;
        }
    }
    // ge[k] = number of reads with t >= k
    let mut ge = vec![0u64; maxlen + 2];
    let mut acc = 0u64;
    for k in (0..=maxlen).rev() {
        ge[k] = acc + counts[k];
        acc += counts[k];
    }
    (0..read_len)
        .map(|p| {
            let need = read_len - p; // >= 1 for every valid read position
            if need <= maxlen { ge[need] as f32 / total as f32 } else { 0.0 }
        })
        .collect()
}

/// Remove one (plane × artefact) pair's contribution from a single measured
/// `[a,c,g,t]` row already normalised to sum to 1.
///
/// `cover` is P(this position carries the injected G). A covered position is
/// forced G, so `g_measured = cover * 1 + (1 - cover) * g_bio`; invert for the
/// biology G fraction and renormalise the remaining three channels over the
/// non-G mass. A no-op when the pair does not reach this position.
fn deconvolve_row(a: &mut [f32; 4], cover: f32) {
    if cover <= 1e-6 || cover >= 1.0 - 1e-6 {
        return;
    }
    let g_bio = ((a[2] - cover) / (1.0 - cover)).clamp(0.0, 1.0);
    let rest_total: f32 = a[0] + a[1] + a[3];
    let rest_bio = (1.0 - g_bio).max(0.0);
    if rest_total > 1e-6 && rest_bio > 1e-6 {
        let scale = rest_bio / rest_total;
        a[0] *= scale;
        a[1] *= scale;
        a[3] *= scale;
    }
    a[2] = g_bio;
    let t: f32 = a.iter().sum();
    if t > 1e-6 {
        for x in a.iter_mut() {
            *x /= t;
        }
    }
}

/// Build a fragment-frame per-position composition from one measured `base`
/// plane, deconvolving BOTH injected runs that contaminate it.
///
/// The plane is the MEASURED read-position composition, so every position an
/// injected head or tail run reaches already contains that artefact. The
/// consumer draws biology from this plane and injects the runs separately, so
/// each pair must be divided out or the artefact is counted twice
/// (docs/GFRAC-MECHANISM.md §3).
fn biased_plane(plane: &[Vec<f32>], head_cov: &[f32], tail_cov: &[f32]) -> Vec<[f32; 4]> {
    plane
        .iter()
        .enumerate()
        .map(|(p, v)| {
            let mut a = [0.0f32; 4];
            for (i, x) in v.iter().take(4).enumerate() {
                if x.is_finite() && *x >= 0.0 {
                    a[i] = *x;
                }
            }
            let s: f32 = a.iter().sum();
            if s > 0.0 {
                for x in a.iter_mut() {
                    *x /= s;
                }
            } else {
                a = [0.25; 4];
            }
            deconvolve_row(&mut a, head_cov.get(p).copied().unwrap_or(0.0));
            deconvolve_row(&mut a, tail_cov.get(p).copied().unwrap_or(0.0));
            a
        })
        .collect()
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
        // A zero-length read is not a library: it emitted empty sequence and
        // quality lines, which is not valid FASTQ, and exited 0.
        if raw.read_len == 0 {
            bail!("profile read_len must be >= 1");
        }
        // SE profiles omit R2: mirror R1 so the parser accepts them (in --single
        // mode the simulator only ever reads R1).
        let r2_head = if raw.head.R2.is_empty() { &raw.head.R1 } else { &raw.head.R2 };
        let head = [Cdf::from_map(&raw.head.R1)?, Cdf::from_map(r2_head)?];
        let r2_tail = if raw.tail.R2.is_empty() { &raw.tail.R1 } else { &raw.tail.R2 };
        let tail = [Cdf::from_map(&raw.tail.R1)?, Cdf::from_map(r2_tail)?];

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
        let r2_posq = if raw.pos_q.R2.is_empty() { raw.pos_q.R1.clone() } else { raw.pos_q.R2.clone() };
        if raw.pos_q.R1.len() != read_len || r2_posq.len() != read_len {
            bail!(
                "pos_q length {} must equal read_len {}",
                raw.pos_q.R1.len(),
                read_len
            );
        }
        // Measured per-position means are fractional; keep them fractional —
        // rounding the curve here would bias every emitted base.
        let pos_q = [raw.pos_q.R1.clone(), r2_posq];

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
                        .as_ref()
                        .and_then(|m| m.tail)
                        .unwrap_or_else(|| mean_of(&pos_q[1]) - raw.tail_q_collapse.R2.unwrap_or(raw.tail_q_collapse.R1)),
                ],
                [
                    raw.deep_head_q
                        .as_ref()
                        .map(|d| d.R1)
                        .or(rm.R1.head)
                        .unwrap_or_else(|| mean_of(&pos_q[0])),
                    raw.deep_head_q
                        .as_ref()
                        .and_then(|d| d.R2)
                        .or_else(|| rm.R2.as_ref().and_then(|m| m.head))
                        .unwrap_or_else(|| mean_of(&pos_q[1])),
                ],
            ),
            None => (
                [
                    mean_of(&pos_q[0]) - raw.tail_q_collapse.R1,
                    mean_of(&pos_q[1]) - raw.tail_q_collapse.R2.unwrap_or(raw.tail_q_collapse.R1),
                ],
                [mean_of(&pos_q[0]), mean_of(&pos_q[1])],
            ),
        };

        let mut mode = match raw.adapter.mode.as_deref() {
            None | Some("readthrough") => AdapterMode::Readthrough,
            Some("none") => AdapterMode::None,
            Some(other) => bail!("unknown adapter mode {other:?}"),
        };
        // An empty adapter sequence cannot be read through: cycling a zero-length
        // string would panic on the modulo. The CLI path already forces `None`
        // for `--adapter-seq ""`; the profile path must do the same.
        let adapter_seq = raw.adapter.seq.into_bytes();
        if adapter_seq.is_empty() {
            mode = AdapterMode::None;
        }
        // The profile's GC is the authority for random DNA (the CLI default is
        // only the fallback), so an out-of-range value must not run clean.
        if !(0.0..=1.0).contains(&raw.gc) {
            bail!("profile gc {} must be in [0,1]", raw.gc);
        }

        let insert = match &raw.insert {
            Some(m) => Some(Cdf::from_map(m)?),
            None => None,
        };
        // N model: an all-zero rate (or an empty histogram) means "emit none",
        // NOT an error — a library with no ambiguous bases is a real library
        // (measured: human_gut/DRR425489 duplicates carry no N at all), and a
        // profile without the field must keep the simulator's stream unchanged.
        let n_run = match &raw.n_run {
            Some(n) if n.read_rate > 0.0 && !n.len.is_empty() => {
                let lead = n.leading.clamp(0.0, 1.0);
                let trail = n.trailing.clamp(0.0, 1.0 - lead);
                Some(NModel {
                    read_rate: n.read_rate.clamp(0.0, 1.0),
                    len: Cdf::from_map(&n.len)?,
                    leading: lead,
                    trailing: trail,
                })
            }
            _ => None,
        };
        // R1 shows the fragment's 5′ end, so its head composition anchors the fragment head;
        // R2 shows the fragment's 3′ end reverse-complemented, so its head composition
        // becomes the fragment tail after A<->T, C<->G.
        //
        // `base.R1[pos]` and `base.R2[pos]` are the MEASURED read-position
        // compositions. Each is a MIXTURE: reads carrying an injected run
        // contribute a forced G at every position the run covers, and reads with
        // no run there contribute biology. The simulator draws biology from this
        // plane and separately injects the runs, so feeding it the mixed value
        // double-counts the artefact. All FOUR (plane × artefact) pairs must be
        // deconvolved out before the plane is used as a biology draw:
        //
        //   base.R1 × R1 head run   (h > p)      -- handled historically
        //   base.R1 × R1 tail run   (t >= L - p) -- was NOT handled
        //   base.R2 × R2 head run   (h > p)      -- was NOT handled
        //   base.R2 × R2 tail run   (t >= L - p) -- was NOT handled
        //
        // Measured (docs/GFRAC-MECHANISM.md §3): with only the first pair
        // handled, DRR083326's R1 terminal G was 0.5270 against a real 0.4245,
        // and the 0.5285 model is exactly the double count. Deconvolving the
        // other three takes it to 0.4140.
        //
        // The head run's mass is P(h > p) -- a run of length h fills positions
        // 0..h-1, including position 0 for an h = 1 run (the boundary rule
        // rewrites index h, the FIRST BIOLOGY position, not inside the run).
        let r1_head = head_cover(&raw.head.R1, read_len);
        let r1_tail = tail_cover(&raw.tail.R1, read_len);
        let r2_head_hist = if raw.head.R2.is_empty() { &raw.head.R1 } else { &raw.head.R2 };
        let r2_tail_hist = if raw.tail.R2.is_empty() { &raw.tail.R1 } else { &raw.tail.R2 };
        let r2_head = head_cover(r2_head_hist, read_len);
        let r2_tail = tail_cover(r2_tail_hist, read_len);

        let head_bias: Vec<[f32; 4]> = match &raw.base {
            Some(b) => match b.get("R1") {
                Some(plane) => biased_plane(plane, &r1_head, &r1_tail),
                None => Vec::new(),
            },
            None => Vec::new(),
        };
        let tail_bias: Vec<[f32; 4]> = match &raw.base {
            Some(b) => match b.get("R2") {
                Some(plane) => biased_plane(plane, &r2_head, &r2_tail)
                    .into_iter()
                    // revcomp: A<->T (0<->3), C<->G (1<->2)
                    .map(|a| [a[3], a[2], a[1], a[0]])
                    .collect(),
                None => Vec::new(),
            },
            None => Vec::new(),
        };

        let read_meanq = [
            Cdf::from_map(&raw.read_meanq.R1)?,
            Cdf::from_map(if raw.read_meanq.R2.is_empty() { &raw.read_meanq.R1 } else { &raw.read_meanq.R2 })?,
        ];
        let read_meanq_mean = [cdf_mean(&read_meanq[0]), cdf_mean(&read_meanq[1])];

        Ok(Profile {
            read_len: raw.read_len,
            gc: raw.gc,
            dup_rate: raw.dup_rate,
            geometry: raw.geometry,
            head_bias,
            tail_bias,
            qual_noise: raw.qual_noise,
            insert,
            head,
            tail,
            joint_tail: JointCdf { cum },
            tail_q_collapse: [raw.tail_q_collapse.R1,
                              raw.tail_q_collapse.R2.unwrap_or(raw.tail_q_collapse.R1)],
            pos_q,
            tail_q_mean,
            head_q_mean,
            read_meanq,
            read_meanq_mean,
            full_length_q: [
                raw.full_length_q.R1.round().clamp(0.0, 41.0) as u8,
                if raw.full_length_q.R2 > 0.0 { raw.full_length_q.R2 } else { raw.full_length_q.R1 }
                    .round().clamp(0.0, 41.0) as u8,
            ],
            adapter: Adapter {
                seq: adapter_seq,
                mode,
            },
            n_run,
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
    ///
    /// Precomputed at load (see the field): this used to walk the CDF for
    /// every read.
    pub fn read_meanq_mean(&self, mate: Mate) -> f64 {
        self.read_meanq_mean[Self::mate_idx(mate)]
    }

    /// Does the profile carry a measured per-position base composition?
    ///
    /// When it does, that plane governs the emitted bases position by position;
    /// the single `gc` fraction only reaches an interior the plane does not
    /// cover, so for a read whose two plane windows meet (`insert <=
    /// 2 x read_len`) `--gc` cannot move a single base. `main` warns and reports
    /// the profile's own gc instead of a value the run never used.
    pub fn has_base_plane(&self) -> bool {
        !self.head_bias.is_empty() || !self.tail_bias.is_empty()
    }
}

/// E[Q] of a CDF stored as `(value, cumulative fraction)` bins.
fn cdf_mean(cdf: &Cdf) -> f64 {
    let mut prev = 0.0;
    let mut mean = 0.0;
    for (len, cum) in &cdf.bins {
        mean += *len as f64 * (cum - prev);
        prev = *cum;
    }
    mean
}

// ---------------------------------------------------------------------------
// Raw JSON schema (serde). Histograms use string keys because JSON.
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
pub struct RawProfile {
    pub read_len: u32,
    pub gc: f64,
    /// measured duplication level; absent in profiles fitted before the field
    #[serde(default)]
    pub dup_rate: Option<f64>,
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
    /// per-position [a,c,g,t] fractions per mate (fastcheck `base_position`)
    #[serde(default)]
    pub base: Option<std::collections::HashMap<String, Vec<Vec<f32>>>>,
    /// measured ambiguous-base (N) emission; absent in profiles fitted before
    /// the N model existed, and then the simulator emits no N at all
    #[serde(default)]
    pub n_run: Option<RawNRun>,
}

/// Raw JSON shape of the N model (`len` is a `{length: count}` histogram like
/// every other profile histogram).
#[derive(Debug, Deserialize)]
pub struct RawNRun {
    pub read_rate: f64,
    pub len: BTreeMap<String, u64>,
    #[serde(default)]
    pub leading: f64,
    #[serde(default)]
    pub trailing: f64,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)] // field names mirror the JSON keys (R1/R2)
pub struct RawRunQMean {
    pub R1: RawRunQMeanMate,
    #[serde(default)]
    pub R2: Option<RawRunQMeanMate>,
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
    /// absent for a single-end profile; the parser mirrors R1
    #[serde(default)]
    pub R2: BTreeMap<String, u64>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
pub struct RawMateFloats {
    pub R1: f64,
    #[serde(default)]
    pub R2: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
pub struct RawMatePosQ {
    pub R1: Vec<f64>,
    /// absent for a single-end profile; the parser mirrors R1
    #[serde(default)]
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

    // ---- consumer-side plane deconvolution (docs/GFRAC-MECHANISM.md §3) ----

    #[test]
    fn run_coverage_masses_match_the_geometry() {
        // {"0":60,"1":30,"2":10}: P(h>p) and P(t >= L-p) for L = 5.
        let m = BTreeMap::from([("0".into(), 60u64), ("1".into(), 30), ("2".into(), 10)]);
        let h = head_cover(&m, 5);
        assert!((h[0] - 0.4).abs() < 1e-6, "P(h>0)=0.4, got {}", h[0]);
        assert!((h[1] - 0.1).abs() < 1e-6, "P(h>1)=0.1, got {}", h[1]);
        assert_eq!(h[2], 0.0);
        let t = tail_cover(&m, 5);
        // p=4 -> need 1 -> P(t>=1)=0.4; p=3 -> need 2 -> 0.1; p=2 -> 0
        assert!((t[4] - 0.4).abs() < 1e-6, "P(t>=1)=0.4, got {}", t[4]);
        assert!((t[3] - 0.1).abs() < 1e-6, "P(t>=2)=0.1, got {}", t[3]);
        assert_eq!(t[2], 0.0);
        // an empty histogram covers nothing
        assert!(head_cover(&BTreeMap::new(), 3).iter().all(|&x| x == 0.0));
        assert!(tail_cover(&BTreeMap::new(), 3).iter().all(|&x| x == 0.0));
    }

    #[test]
    fn deconvolve_row_inverts_the_mixture() {
        // measured = cover*1 + (1-cover)*bio  ->  deconvolve must recover bio
        for &(cover, bio) in &[(0.2f32, 0.1f32), (0.4, 0.3), (0.25, 0.6)] {
            let measured = cover + (1.0 - cover) * bio;
            let mut a = [(1.0 - measured) / 3.0, (1.0 - measured) / 3.0, measured, (1.0 - measured) / 3.0];
            deconvolve_row(&mut a, cover);
            assert!((a[2] - bio).abs() < 1e-5, "cover {cover}: got {}, want {bio}", a[2]);
            let s: f32 = a.iter().sum();
            assert!((s - 1.0).abs() < 1e-5, "row must stay normalised, sum {s}");
        }
        // a pair that does not reach the position is a no-op
        let mut a = [0.2, 0.2, 0.4, 0.2];
        deconvolve_row(&mut a, 0.0);
        assert_eq!(a, [0.2, 0.2, 0.4, 0.2]);
    }

    #[test]
    fn profile_deconvolves_both_planes_for_both_injected_runs() {
        // read_len 4; R1 head {"0":50,"2":50} and tail {"0":70,"2":30};
        // R2 has neither a head nor a tail run.
        // The measured plane's terminal R1 base is G = 0.5. At p = 3 the head
        // run does not reach (P(h>3)=0) and the tail run covers with
        // P(t>=1) = 0.3, so the biology value must be (0.5-0.3)/0.7 = 0.2857.
        // Before this fix that terminal base was used as 0.5 and the injected
        // tail was counted on top of it -- the DRR083326 double count.
        let q = vec!["30"; 4].join(",");
        let json = format!(
            r#"{{
            "read_len": 4, "gc": 0.5, "geometry": "tail3p_both",
            "head": {{"R1": {{"0": 50, "2": 50}}, "R2": {{"0": 100}}}},
            "tail": {{"R1": {{"0": 70, "2": 30}}, "R2": {{"0": 100}}}},
            "joint_tail": {{"0,0": 100}},
            "tail_q_collapse": {{"R1": 0.0, "R2": 0.0}},
            "pos_q": {{"R1": [{q}], "R2": [{q}]}},
            "read_meanq": {{"R1": {{"30": 100}}, "R2": {{"30": 100}}}},
            "full_length_rate": 0.0,
            "full_length_q": {{"R1": 20, "R2": 20}},
            "adapter": {{"seq": "CTGTCTCTTATACACATCT", "mode": "none"}},
            "base": {{
                "R1": [[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],
                       [0.25,0.25,0.25,0.25],[0.1667,0.1667,0.5,0.1666]],
                "R2": [[0.1,0.1,0.7,0.1],[0.25,0.25,0.25,0.25],
                       [0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]]
            }}
        }}"#,
            q = q
        );
        let p = Profile::from_raw(serde_json::from_str::<RawProfile>(&json).unwrap()).unwrap();
        assert_eq!(p.head_bias.len(), 4);
        assert_eq!(p.tail_bias.len(), 4);
        // R1 terminal: tail pair removed, not just the head pair
        assert!((p.head_bias[3][2] - 0.2857).abs() < 1e-3,
                "R1 terminal biology G should be 0.2857, got {}", p.head_bias[3][2]);
        // R1 position 0: head covers with P(h>0) = 0.5, tail does not reach
        assert!((p.head_bias[0][2] - 0.0).abs() < 1e-6,
                "R1 pos0 is fully covered by the head run, got {}", p.head_bias[0][2]);
        // R2's plane is revcomp'd into the fragment tail: measured G 0.7 at R2
        // position 0 becomes C (index 1) at tail_bias[0]; no run reaches R2.
        assert!((p.tail_bias[0][1] - 0.7).abs() < 1e-6,
                "revcomp must map G->C, got {:?}", p.tail_bias[0]);
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
    fn n_model_parses_and_absent_or_zero_is_none() {
        // absent -> None (the simulator must then draw nothing)
        let p = test_profile();
        assert!(p.n_run.is_none());

        // present -> the measured rate and run-length CDF are carried
        let mut raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        raw.n_run = Some(RawNRun {
            read_rate: 0.024,
            len: BTreeMap::from([("1".into(), 90u64), ("150".into(), 10)]),
            leading: 0.1,
            trailing: 0.6,
        });
        let p = Profile::from_raw(raw).unwrap();
        let n = p.n_run.as_ref().expect("model present");
        assert!((n.read_rate - 0.024).abs() < 1e-12);
        assert_eq!(n.len.total, 100);
        assert!((n.leading - 0.1).abs() < 1e-12);
        assert!((n.trailing - 0.6).abs() < 1e-12);

        // a zero rate is "no N", not an error: some libraries genuinely have none
        let mut raw: RawProfile = serde_json::from_str(&sample_profile_json()).unwrap();
        raw.n_run = Some(RawNRun {
            read_rate: 0.0,
            len: BTreeMap::from([("1".into(), 100u64)]),
            leading: 0.0,
            trailing: 0.0,
        });
        assert!(Profile::from_raw(raw).unwrap().n_run.is_none());
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
