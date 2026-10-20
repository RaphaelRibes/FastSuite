//! Instrument counters (FS-11): the primitives that make "same stream" and
//! "same amount of work" checkable instead of asserted.
//!
//! [`TrackedRng`] is an `Rng` wrapper that counts the words drawn. It is a pure
//! pass-through — every method forwards to the inner generator — so the byte
//! stream is identical with the counter on or off. Counting is disabled by
//! default (`--stats` turns it on) and the increment sits behind a branch that
//! never fires, which is what keeps the default path free.
//!
//! The count is the number of *underlying words* (`next_u64`/`next_u32`/
//! `fill_bytes` calls), not of `gen()` calls: that is the quantity a
//! determinism argument actually needs (FS-2's per-pair streams, FS-5's
//! geometric gaps), and it is cheap to maintain exactly.
//!
//! The other counters (substitutions, mean Q, gzip bytes) live where their
//! state does — `QCache` for the quality model, the writer's files for the
//! byte count — so a run's counters cannot leak into another run's.

use rand::RngCore;

/// An `Rng` that forwards every draw and counts the underlying words.
///
/// `enabled == false` is the default: the counter is then a no-op behind an
/// always-false branch, and the draws are exactly the inner generator's.
pub struct TrackedRng<R> {
    inner: R,
    draws: u64,
    enabled: bool,
}

impl<R: RngCore> TrackedRng<R> {
    pub fn new(inner: R, enabled: bool) -> Self {
        TrackedRng {
            inner,
            draws: 0,
            enabled,
        }
    }

    /// Underlying words drawn so far (0 while the counter is disabled).
    pub fn draws(&self) -> u64 {
        self.draws
    }

    /// The wrapped generator.
    pub fn inner(&self) -> &R {
        &self.inner
    }

    #[inline]
    fn tick(&mut self) {
        if self.enabled {
            self.draws += 1;
        }
    }
}

impl<R: RngCore> RngCore for TrackedRng<R> {
    #[inline]
    fn next_u32(&mut self) -> u32 {
        let v = self.inner.next_u32();
        self.tick();
        v
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let v = self.inner.next_u64();
        self.tick();
        v
    }

    #[inline]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.inner.fill_bytes(dest);
        self.tick();
    }

    #[inline]
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        let r = self.inner.try_fill_bytes(dest);
        self.tick();
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    /// The counter must not consume or reorder a single draw: same seed, same
    /// values, counted or not.
    #[test]
    fn counting_does_not_change_the_stream() {
        let mut plain = ChaCha8Rng::seed_from_u64(7);
        let mut off = TrackedRng::new(ChaCha8Rng::seed_from_u64(7), false);
        let mut on = TrackedRng::new(ChaCha8Rng::seed_from_u64(7), true);

        let mut a = Vec::new();
        let mut b = Vec::new();
        let mut c = Vec::new();
        for _ in 0..1000 {
            a.push(plain.gen::<f64>());
            b.push(off.gen::<f64>());
            c.push(on.gen::<f64>());
        }
        assert_eq!(a, b, "the disabled wrapper must be a pass-through");
        assert_eq!(a, c, "the enabled wrapper must be a pass-through");
        assert_eq!(off.draws(), 0, "counting is off by default");
        assert!(on.draws() >= 1000, "one word per f64 draw: {}", on.draws());
    }
}
