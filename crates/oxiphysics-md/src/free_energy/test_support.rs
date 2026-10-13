// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared helpers for the statistical unit tests of the free-energy estimators.

use rand::{Rng, RngExt};
use std::f64::consts::TAU;

/// Draw one standard normal variate with the Box-Muller transform.
///
/// `u1` is drawn from `(0, 1]` so that `ln(u1)` is always finite.
pub(crate) fn standard_normal<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    let u1: f64 = 1.0 - rng.random::<f64>();
    let u2: f64 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

/// Sample mean and (unbiased) sample standard deviation of `values`.
pub(crate) fn mean_and_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, var.sqrt())
}
