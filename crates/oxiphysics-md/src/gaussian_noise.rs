// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Gaussian random variates for the Langevin-type thermostats of the
//! path-integral and polarizable-model modules.
//!
//! Langevin and Ornstein-Uhlenbeck (OU) integrators need *standard normal*
//! noise.  A uniform draw on `[-1, 1)` has variance 1/3 and therefore makes a
//! thermostat equilibrate at one third of its target temperature.  The PILE-L
//! ring-polymer thermostat and the extended-Lagrangian Drude Langevin kick draw
//! their noise through [`standard_normal`] instead; the unit tests of the Drude
//! thermostat in `polarizable_md` build three-component noise vectors with the
//! test-only helper `standard_normal_vec3`.

use rand::{Rng, RngExt};
use std::f64::consts::TAU;

/// Draw one standard normal variate `N(0, 1)` with the Box-Muller transform.
///
/// `u1` is drawn from `(0, 1]` so that `ln(u1)` is always finite, and `u2`
/// from `[0, 1)`.  The sine partner of the Box-Muller pair is discarded so
/// that the generator keeps no hidden state and a seeded `rng` fully
/// determines the stream.
pub(crate) fn standard_normal<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    let u1: f64 = 1.0 - rng.random::<f64>();
    let u2: f64 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

/// Draw three independent standard normal variates (one per Cartesian axis).
///
/// Used by the unit tests that exercise the thermostats with caller-supplied
/// noise vectors.
#[cfg(test)]
pub(crate) fn standard_normal_vec3<R: Rng + ?Sized>(rng: &mut R) -> [f64; 3] {
    [
        standard_normal(rng),
        standard_normal(rng),
        standard_normal(rng),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    /// Sample mean, variance and kurtosis of `n` draws.
    fn moments(n: usize, seed: u64) -> (f64, f64, f64) {
        let mut rng = StdRng::seed_from_u64(seed);
        let samples: Vec<f64> = (0..n).map(|_| standard_normal(&mut rng)).collect();
        let nf = n as f64;
        let mean = samples.iter().sum::<f64>() / nf;
        let var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / nf;
        let m4 = samples.iter().map(|x| (x - mean).powi(4)).sum::<f64>() / nf;
        (mean, var, m4 / (var * var))
    }

    #[test]
    fn standard_normal_has_unit_variance_and_gaussian_kurtosis() {
        let (mean, var, kurtosis) = moments(400_000, 7);
        assert!(mean.abs() < 0.01, "mean = {mean}");
        // A uniform draw on [-1, 1) would give variance 1/3 and kurtosis 1.8.
        assert!((var - 1.0).abs() < 0.01, "variance = {var}");
        assert!((kurtosis - 3.0).abs() < 0.1, "kurtosis = {kurtosis}");
    }

    #[test]
    fn standard_normal_is_reproducible_for_a_fixed_seed() {
        let mut a = StdRng::seed_from_u64(99);
        let mut b = StdRng::seed_from_u64(99);
        for _ in 0..16 {
            assert_eq!(standard_normal(&mut a), standard_normal(&mut b));
        }
    }

    #[test]
    fn standard_normal_is_always_finite() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..100_000 {
            assert!(standard_normal(&mut rng).is_finite());
        }
    }

    #[test]
    fn standard_normal_vec3_components_are_independent_draws() {
        let mut rng = StdRng::seed_from_u64(5);
        let v = standard_normal_vec3(&mut rng);
        assert!(v[0] != v[1] && v[1] != v[2]);
        let n = 200_000;
        let mut cross = 0.0;
        for _ in 0..n {
            let w = standard_normal_vec3(&mut rng);
            cross += w[0] * w[1];
        }
        assert!((cross / n as f64).abs() < 0.01, "cross-correlation");
    }
}
