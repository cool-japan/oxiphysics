// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Deterministic random-number toolkit of the self-contained MD drivers.
//!
//! The stochastic parts of [`MdSimulationPlain`](super::MdSimulationPlain) --
//! the canonical-sampling velocity-rescaling thermostat and the Monte-Carlo
//! barostat -- draw uniform, Gaussian and chi-squared deviates from a
//! [`SimRng`].  The generator is seeded explicitly, so a simulation is
//! reproducible bit for bit and independent of the version of any external
//! random-number crate.

use std::f64::consts::PI;

/// Maximum number of rejection rounds of the Marsaglia-Tsang gamma sampler.
///
/// The acceptance probability of one round is above 95 % for every shape
/// `>= 1`, so exhausting the budget has a probability below `10^-1000`; the
/// bound only makes the loop provably finite.
const MAX_GAMMA_ROUNDS: usize = 10_000;

/// SplitMix64 pseudo-random generator with Gaussian and gamma samplers.
///
/// SplitMix64 (Steele, Lea & Flood 2014) passes the BigCrush battery, has a
/// 64-bit state, and every seed (including zero) yields a full-period
/// sequence.
#[derive(Debug, Clone)]
pub struct SimRng {
    state: u64,
}

impl SimRng {
    /// Create a generator from `seed`.
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Restart the generator from `seed`.
    pub fn reseed(&mut self, seed: u64) {
        self.state = seed;
    }

    /// Next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform deviate in `[0, 1)` built from the upper 53 bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform deviate in `(0, 1]`; never zero, so its logarithm is finite.
    pub fn uniform_open(&mut self) -> f64 {
        1.0 - self.uniform()
    }

    /// Standard normal deviate `N(0, 1)` (Box-Muller transform, cosine branch).
    pub fn normal(&mut self) -> f64 {
        let u1 = self.uniform_open();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }

    /// Draw from `Gamma(shape, 1)` (mean and variance both equal to `shape`).
    ///
    /// Uses the Marsaglia-Tsang squeeze/rejection method (ACM TOMS 26, 363
    /// (2000)) for `shape >= 1` and the boost
    /// `Gamma(a) = Gamma(a + 1) U^(1/a)` for `0 < shape < 1`.  A non-positive
    /// or non-finite `shape` yields `0.0`.
    pub fn gamma(&mut self, shape: f64) -> f64 {
        if !shape.is_finite() || shape <= 0.0 {
            return 0.0;
        }
        if shape < 1.0 {
            let boosted = self.gamma(shape + 1.0);
            let u = self.uniform_open();
            return boosted * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        for _ in 0..MAX_GAMMA_ROUNDS {
            let x = self.normal();
            let t = 1.0 + c * x;
            if t <= 0.0 {
                continue;
            }
            let v = t * t * t;
            let u = self.uniform_open();
            if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
                return d * v;
            }
        }
        shape
    }

    /// Draw a chi-squared deviate with `dof` (real, `>= 0`) degrees of freedom,
    /// `chi2(dof) = 2 Gamma(dof / 2)`.  `dof <= 0` gives `0.0`.
    pub fn chi_squared(&mut self, dof: f64) -> f64 {
        2.0 * self.gamma(0.5 * dof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sample mean and (population) variance of `n` draws of `draw`.
    fn moments(n: usize, mut draw: impl FnMut() -> f64) -> (f64, f64) {
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for _ in 0..n {
            let x = draw();
            sum += x;
            sum_sq += x * x;
        }
        let mean = sum / n as f64;
        (mean, sum_sq / n as f64 - mean * mean)
    }

    #[test]
    fn test_same_seed_gives_identical_stream() {
        let mut a = SimRng::new(42);
        let mut b = SimRng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        let mut c = SimRng::new(43);
        let mut d = SimRng::new(42);
        let differs = (0..16).any(|_| c.next_u64() != d.next_u64());
        assert!(differs, "different seeds must give different streams");
    }

    #[test]
    fn test_reseed_restarts_the_sequence() {
        let mut rng = SimRng::new(7);
        let first: Vec<u64> = (0..5).map(|_| rng.next_u64()).collect();
        rng.reseed(7);
        let again: Vec<u64> = (0..5).map(|_| rng.next_u64()).collect();
        assert_eq!(first, again);
    }

    #[test]
    fn test_uniform_ranges_and_moments() {
        let mut rng = SimRng::new(1);
        for _ in 0..10_000 {
            let u = rng.uniform();
            assert!((0.0..1.0).contains(&u), "uniform out of range: {u}");
            let v = rng.uniform_open();
            assert!(v > 0.0 && v <= 1.0, "uniform_open out of range: {v}");
        }
        let (mean, var) = moments(200_000, || rng.uniform());
        assert!((mean - 0.5).abs() < 0.005, "mean = {mean}");
        assert!((var - 1.0 / 12.0).abs() < 0.002, "variance = {var}");
    }

    #[test]
    fn test_normal_moments() {
        let mut rng = SimRng::new(2);
        let (mean, var) = moments(400_000, || rng.normal());
        assert!(mean.abs() < 0.01, "mean = {mean}");
        assert!((var - 1.0).abs() < 0.02, "variance = {var}");
    }

    #[test]
    fn test_gamma_moments_for_several_shapes() {
        // Gamma(a, 1) has mean a and variance a, including shapes below one.
        let mut rng = SimRng::new(3);
        for shape in [0.3, 1.0, 2.5, 20.0, 150.0] {
            let (mean, var) = moments(200_000, || rng.gamma(shape));
            assert!(
                (mean - shape).abs() < 0.02 * shape.max(1.0),
                "shape {shape}: mean = {mean}"
            );
            assert!(
                (var - shape).abs() < 0.04 * shape.max(1.0),
                "shape {shape}: variance = {var}"
            );
        }
    }

    #[test]
    fn test_chi_squared_moments_and_degenerate_dof() {
        let mut rng = SimRng::new(4);
        // chi2(k) has mean k and variance 2k.
        let (mean, var) = moments(200_000, || rng.chi_squared(9.0));
        assert!((mean - 9.0).abs() < 0.1, "mean = {mean}");
        assert!((var - 18.0).abs() < 0.6, "variance = {var}");
        assert_eq!(rng.chi_squared(0.0), 0.0);
        assert_eq!(rng.chi_squared(-3.0), 0.0);
        assert_eq!(rng.gamma(f64::NAN), 0.0);
    }
}
