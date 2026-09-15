//! Random sampling utilities: CDF sampling, joint (R1×R2) sampling, gaussian.
//!
//! All sampling uses a `rand::Rng`; determinism comes from the caller passing
//! a seeded `ChaCha8Rng` (see `cli.rs` / main).

use rand::Rng;

use crate::profile::{Cdf, JointCdf};

/// Sample a length from a cumulative distribution.
///
/// `cdf.bins` is monotonic in cumulative fraction; binary search for the
/// first bin whose cumulative >= u ∈ [0, 1).
#[inline]
pub fn sample_cdf<R: Rng + ?Sized>(cdf: &Cdf, rng: &mut R) -> u32 {
    let u: f64 = rng.gen();
    // binary search: first bin with cum >= u
    match cdf
        .bins
        .binary_search_by(|(_, c)| c.partial_cmp(&u).unwrap())
    {
        Ok(i) => cdf.bins[i].0,
        Err(i) => cdf.bins[i.min(cdf.bins.len() - 1)].0,
    }
}

/// Sample an (R1, R2) tail-length pair from the joint cumulative grid.
#[inline]
pub fn sample_joint<R: Rng + ?Sized>(joint: &JointCdf, rng: &mut R) -> (u32, u32) {
    let u: f64 = rng.gen();
    let idx = joint
        .cum
        .binary_search_by(|c| c.partial_cmp(&u).unwrap())
        .unwrap_or_else(|i| i.min(joint.cum.len() - 1));
    let r1 = (idx / JointCdf::N) as u32;
    let r2 = (idx % JointCdf::N) as u32;
    (r1, r2)
}

/// Sample a gaussian N(mean, sd), clamped to [min, max].
pub fn sample_gauss<R: Rng + ?Sized>(mean: f64, sd: f64, min: f64, max: f64, rng: &mut R) -> f64 {
    use rand_distr::{Distribution, Normal};
    // normal() panics on invalid params; sd > 0 checked by caller config
    let n = Normal::new(mean, sd).expect("valid gaussian params");
    n.sample(rng).clamp(min, max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profile::{Cdf, JOINT_CAP};
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(42)
    }

    fn cdf_from(probs: &[(u32, f64)]) -> Cdf {
        // probs: (len, probability) — build dense counts
        let max = probs.iter().map(|(l, _)| *l).max().unwrap() as usize;
        let mut dense = vec![0u64; max + 1];
        for (l, p) in probs {
            dense[*l as usize] = (*p * 1_000_000.0) as u64;
        }
        Cdf::from_dense(&dense)
    }

    #[test]
    fn sample_cdf_matches_probabilities() {
        // P(0)=0.5, P(1)=0.25, P(2)=0.25
        let c = cdf_from(&[(0, 0.5), (1, 0.25), (2, 0.25)]);
        let mut rng = rng();
        let n = 1_000_000;
        let mut counts = [0u64; 3];
        for _ in 0..n {
            counts[sample_cdf(&c, &mut rng) as usize] += 1;
        }
        for (i, p) in [0.5, 0.25, 0.25].iter().enumerate() {
            let frac = counts[i] as f64 / n as f64;
            assert!((frac - p).abs() < 0.005, "bin {i}: {frac} vs {p}");
        }
    }

    #[test]
    fn sample_cdf_deterministic() {
        let c = cdf_from(&[(0, 0.5), (1, 0.25), (2, 0.25)]);
        let mut a = rng();
        let mut b = ChaCha8Rng::seed_from_u64(42);
        for _ in 0..100 {
            assert_eq!(sample_cdf(&c, &mut a), sample_cdf(&c, &mut b));
        }
    }

    #[test]
    fn sample_joint_marginals_match() {
        // joint: {0,0}=36, {1,1}=12, {2,2}=2  (diagonal only, total 50)
        let mut grid = vec![0f64; JointCdf::N * JointCdf::N];
        let put = |grid: &mut Vec<f64>, r1: u32, r2: u32, v: f64| {
            grid[r1 as usize * JointCdf::N + r2 as usize] = v;
        };
        put(&mut grid, 0, 0, 36.0);
        put(&mut grid, 1, 1, 12.0);
        put(&mut grid, 2, 2, 2.0);
        let mut cum = Vec::with_capacity(grid.len());
        let mut acc = 0.0;
        for &c in &grid {
            acc += c / 50.0;
            cum.push(acc);
        }
        let joint = JointCdf { cum };

        let mut rng = rng();
        let n = 500_000;
        let mut marginal1 = [0u64; 3];
        let mut marginal2 = [0u64; 3];
        for _ in 0..n {
            let (a, b) = sample_joint(&joint, &mut rng);
            assert_eq!(a, b, "diagonal-only joint must return equal pair");
            marginal1[a as usize] += 1;
            marginal2[b as usize] += 1;
        }
        // P(0)=0.72, P(1)=0.24, P(2)=0.04
        for (i, p) in [0.72, 0.24, 0.04].iter().enumerate() {
            let f = marginal1[i] as f64 / n as f64;
            assert!((f - p).abs() < 0.005, "R1 bin {i}: {f} vs {p}");
            let f = marginal2[i] as f64 / n as f64;
            assert!((f - p).abs() < 0.005, "R2 bin {i}: {f} vs {p}");
        }
    }

    #[test]
    fn sample_joint_caps_at_64() {
        // all mass at (100, 100) -> capped to (64, 64)
        let mut grid = vec![0f64; JointCdf::N * JointCdf::N];
        grid[100.min(JOINT_CAP) as usize * JointCdf::N + 100.min(JOINT_CAP) as usize] = 1.0;
        let mut cum = Vec::with_capacity(grid.len());
        let mut acc = 0.0;
        for &c in &grid {
            acc += c;
            cum.push(acc);
        }
        let joint = JointCdf { cum };
        let mut rng = rng();
        for _ in 0..100 {
            let (a, b) = sample_joint(&joint, &mut rng);
            assert_eq!((a, b), (JOINT_CAP, JOINT_CAP));
        }
    }
}
