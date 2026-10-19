// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared kernels of the Bussi-Donadio-Parrinello stochastic velocity
//! rescaling thermostats.
//!
//! The canonical-sampling velocity-rescaling (CSVR) thermostat evolves the
//! total kinetic energy `K` of `N_f` degrees of freedom with the stochastic
//! differential equation
//!
//! ```text
//! dK = (K̄ − K)·dt/τ + 2·sqrt(K·K̄/N_f)·dW/√τ ,        K̄ = N_f·k_B·T / 2 ,
//! ```
//!
//! whose stationary distribution is the canonical gamma distribution
//! `P(K) ∝ K^(N_f/2 − 1)·exp(−K / k_B T)` with `⟨K⟩ = K̄` and
//! `Var(K) = 2·K̄² / N_f`.  Integrating it exactly over a step `Δt` gives
//! (Bussi, Donadio & Parrinello, J. Chem. Phys. 126, 014101 (2007), eq. A7)
//!
//! ```text
//! K' = c·K + (1 − c)·K̄/N_f·(R₁² + S) + 2·R₁·sqrt(c·(1 − c)·K·K̄ / N_f) ,
//! ```
//!
//! where `c = exp(−Δt/τ)`, `R₁ ~ N(0, 1)` and `S ~ χ²(N_f − 1)` is
//! independent of `R₁`.  Every velocity is then multiplied by
//! `α = sqrt(K'/K)`.
//!
//! [`CsvrThermostat`](super::CsvrThermostat) and
//! [`SvrThermostat`](super::SvrThermostat) both call [`sample_alpha2`], so the
//! two thermostats cannot drift apart.  The `χ²` deviate is sampled exactly
//! (as `2·Gamma((N_f − 1)/2)` with the Marsaglia-Tsang method), for any
//! number of degrees of freedom.

use std::f64::consts::PI;

/// Maximum number of rejection rounds of the Marsaglia-Tsang sampler.
///
/// The acceptance probability is above 95 % for every shape `≥ 1`, so the
/// chance of exhausting this budget is below `10^-1000`; the bound only
/// exists so that the loop is provably finite.
const MAX_GAMMA_ROUNDS: usize = 10_000;

/// Advance the 64-bit linear congruential generator and return a uniform
/// deviate in `[0, 1)` built from its upper 53 bits.
pub(super) fn lcg_uniform(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

/// Uniform deviate in `(0, 1]`; never zero, so `ln` of it is always finite.
fn lcg_uniform_open(state: &mut u64) -> f64 {
    1.0 - lcg_uniform(state)
}

/// Standard normal deviate `N(0, 1)` (Box-Muller transform, cosine branch).
pub(super) fn lcg_normal(state: &mut u64) -> f64 {
    let u1 = lcg_uniform_open(state);
    let u2 = lcg_uniform(state);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// Draw from `Gamma(shape, 1)` (mean and variance both equal to `shape`).
///
/// Uses the Marsaglia-Tsang squeeze/rejection method (ACM TOMS 26, 363
/// (2000)) for `shape ≥ 1` and the boost `Gamma(a) = Gamma(a + 1)·U^(1/a)`
/// for `shape < 1`.  A non-positive or non-finite `shape` yields `0.0`.
pub(super) fn sample_gamma(state: &mut u64, shape: f64) -> f64 {
    if !shape.is_finite() || shape <= 0.0 {
        return 0.0;
    }
    if shape < 1.0 {
        let boosted = sample_gamma(state, shape + 1.0);
        let u = lcg_uniform_open(state);
        return boosted * u.powf(1.0 / shape);
    }
    let d = shape - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    for _ in 0..MAX_GAMMA_ROUNDS {
        let x = lcg_normal(state);
        let t = 1.0 + c * x;
        if t <= 0.0 {
            continue;
        }
        let v = t * t * t;
        let u = lcg_uniform_open(state);
        if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
            return d * v;
        }
    }
    shape
}

/// Draw a chi-squared deviate with `dof` (real, `≥ 0`) degrees of freedom,
/// `χ²(dof) = 2·Gamma(dof / 2)`.  `dof ≤ 0` gives `0.0`.
pub(super) fn sample_chi2(state: &mut u64, dof: f64) -> f64 {
    2.0 * sample_gamma(state, 0.5 * dof)
}

/// Per-step decay factor `c = exp(−dt/τ)` of the kinetic-energy relaxation.
///
/// No elapsed time (`dt ≤ 0` or NaN) means no coupling (`c = 1`); a
/// non-positive coupling time means instantaneous coupling, i.e. a complete
/// resampling of the kinetic energy (`c = 0`).
pub(super) fn bussi_decay(dt: f64, tau: f64) -> f64 {
    if dt.is_nan() || dt <= 0.0 {
        return 1.0;
    }
    if tau.is_nan() || tau <= 0.0 {
        return 0.0;
    }
    (-dt / tau).exp()
}

/// Squared velocity-scaling factor `α² = K'/K` of one exact Bussi step.
///
/// * `k`     – current kinetic energy `K > 0`.
/// * `k_bar` – target kinetic energy `K̄ = N_f·k_B·T / 2 ≥ 0`.
/// * `nf`    – number of degrees of freedom `N_f > 0`.
/// * `c`     – decay factor `exp(−Δt/τ)` in `[0, 1]`.
/// * `r1`    – a standard normal deviate `R₁`.
/// * `chi2`  – a `χ²(N_f − 1)` deviate.
///
/// The update is evaluated in the sum-of-squares form
///
/// ```text
/// α² = (√c + R₁·b)² + b²·χ² ,        b² = (1 − c)·K̄ / (K·N_f) ,
/// ```
///
/// which is algebraically identical to eq. A7 of Bussi *et al.* and is
/// non-negative by construction, so no clamping is required.
pub(super) fn bussi_alpha2(k: f64, k_bar: f64, nf: f64, c: f64, r1: f64, chi2: f64) -> f64 {
    let b2 = (1.0 - c) * k_bar / (k * nf);
    let lead = c.sqrt() + r1 * b2.sqrt();
    lead * lead + b2 * chi2
}

/// Sample `α²` for one thermostat step, drawing `R₁` and `χ²(N_f − 1)` from
/// the LCG `state`.  See [`bussi_alpha2`] for the argument meaning.
pub(super) fn sample_alpha2(state: &mut u64, k: f64, k_bar: f64, nf: f64, c: f64) -> f64 {
    let r1 = lcg_normal(state);
    let chi2 = sample_chi2(state, nf - 1.0);
    bussi_alpha2(k, k_bar, nf, c, r1, chi2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomSet;
    use crate::thermostat::{CsvrThermostat, SvrThermostat, Thermostat};
    use oxiphysics_core::Vec3;

    /// Sample mean and population variance (two-pass).
    fn mean_and_variance(samples: &[f64]) -> (f64, f64) {
        let n = samples.len() as f64;
        let mean = samples.iter().sum::<f64>() / n;
        let var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        (mean, var)
    }

    /// `n` unit-mass atoms whose instantaneous temperature is `temp` (`k_B = 1`).
    fn ideal_gas(n: usize, temp: f64) -> AtomSet {
        let v = temp.sqrt();
        let mut atoms = AtomSet::new();
        for i in 0..n {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            atoms.add_atom(
                Vec3::zeros(),
                Vec3::new(sign * v, sign * v, sign * v),
                1.0,
                0.0,
                0,
            );
        }
        atoms
    }

    /// Run `thermostat` on a force-free gas and return the mean and variance of
    /// the total kinetic energy over `samples` steps after `burn_in` steps.
    fn kinetic_energy_statistics(
        thermostat: &mut dyn Thermostat,
        atoms: &mut AtomSet,
        target_temp: f64,
        dt: f64,
        burn_in: usize,
        samples: usize,
    ) -> (f64, f64) {
        for _ in 0..burn_in {
            thermostat.apply(atoms, target_temp, dt, 1.0);
        }
        let mut series = Vec::with_capacity(samples);
        for _ in 0..samples {
            thermostat.apply(atoms, target_temp, dt, 1.0);
            series.push(atoms.kinetic_energy());
        }
        mean_and_variance(&series)
    }

    #[test]
    fn lcg_uniform_stays_in_unit_interval_with_half_mean() {
        let mut state = 17_u64;
        let draws: Vec<f64> = (0..200_000).map(|_| lcg_uniform(&mut state)).collect();
        assert!(draws.iter().all(|u| (0.0..1.0).contains(u)));
        let (mean, var) = mean_and_variance(&draws);
        assert!((mean - 0.5).abs() < 5e-3, "mean = {mean}");
        assert!((var - 1.0 / 12.0).abs() < 2e-3, "variance = {var}");
    }

    #[test]
    fn lcg_normal_has_unit_variance_and_gaussian_kurtosis() {
        let mut state = 99_u64;
        let draws: Vec<f64> = (0..400_000).map(|_| lcg_normal(&mut state)).collect();
        let (mean, var) = mean_and_variance(&draws);
        let m4 = draws.iter().map(|x| (x - mean).powi(4)).sum::<f64>() / draws.len() as f64;
        assert!(mean.abs() < 1e-2, "mean = {mean}");
        assert!((var - 1.0).abs() < 1.5e-2, "variance = {var}");
        assert!(
            (m4 / (var * var) - 3.0).abs() < 0.1,
            "kurtosis = {}",
            m4 / (var * var)
        );
        assert!(draws.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn chi2_moments_are_exact_for_every_degree_of_freedom() {
        // Mean = dof and variance = 2 dof, including the shape < 1 branch
        // (dof = 0.5, 1) and a large dof where Wilson-Hilferty would be used.
        for (seed, dof) in [
            (1_u64, 0.5_f64),
            (2, 1.0),
            (3, 2.0),
            (4, 5.0),
            (5, 29.0),
            (6, 299.0),
            (7, 2999.0),
        ] {
            let mut state = seed;
            let draws: Vec<f64> = (0..200_000).map(|_| sample_chi2(&mut state, dof)).collect();
            let (mean, var) = mean_and_variance(&draws);
            assert!(draws.iter().all(|x| x.is_finite() && *x >= 0.0));
            assert!((mean / dof - 1.0).abs() < 0.02, "chi2({dof}) mean = {mean}");
            assert!(
                (var / (2.0 * dof) - 1.0).abs() < 0.05,
                "chi2({dof}) variance = {var}, expected {}",
                2.0 * dof
            );
        }
    }

    #[test]
    fn chi2_with_no_degrees_of_freedom_is_zero() {
        let mut state = 5_u64;
        assert_eq!(sample_chi2(&mut state, 0.0), 0.0);
        assert_eq!(sample_chi2(&mut state, -3.0), 0.0);
        assert_eq!(sample_chi2(&mut state, f64::NAN), 0.0);
        assert_eq!(sample_gamma(&mut state, f64::INFINITY), 0.0);
    }

    #[test]
    fn decay_factor_limits() {
        assert_eq!(bussi_decay(0.0, 1.0), 1.0);
        assert_eq!(bussi_decay(-1.0, 1.0), 1.0);
        assert_eq!(bussi_decay(f64::NAN, 1.0), 1.0);
        assert_eq!(bussi_decay(0.5, 0.0), 0.0);
        assert_eq!(bussi_decay(0.5, -2.0), 0.0);
        assert_eq!(bussi_decay(0.5, f64::INFINITY), 1.0);
        assert!((bussi_decay(0.5, 2.0) - (-0.25_f64).exp()).abs() < 1e-15);
    }

    /// The sum-of-squares kernel must reproduce eq. A7 of Bussi et al. verbatim.
    #[test]
    fn alpha2_matches_equation_a7() {
        let cases: [(f64, f64, f64, f64, f64, f64); 4] = [
            (3.0, 5.0, 12.0, 0.8, 0.7, 9.0),
            (10.0, 2.0, 30.0, 0.25, -1.3, 31.0),
            (0.4, 6.0, 3.0, 0.95, 2.1, 1.7),
            (7.5, 7.5, 300.0, 0.0, 0.2, 280.0),
        ];
        for (k, k_bar, nf, c, r1, chi2) in cases {
            // Eq. A7 as printed: K' = K + (1-c)(K̄ (R1² + S)/N_f - K) + 2 R1 sqrt(K K̄/N_f (1-c) c).
            let k_new = k
                + (1.0 - c) * (k_bar * (chi2 + r1 * r1) / nf - k)
                + 2.0 * r1 * (k * k_bar / nf * (1.0 - c) * c).sqrt();
            let alpha2 = bussi_alpha2(k, k_bar, nf, c, r1, chi2);
            assert!(
                (alpha2 * k - k_new).abs() < 1e-12 * k_new.abs().max(1.0),
                "K' mismatch: kernel {} vs A7 {k_new}",
                alpha2 * k
            );
        }
    }

    #[test]
    fn alpha2_is_nonnegative_and_vanishes_at_the_boundary_draw() {
        let (k, k_bar, nf, c) = (4.0_f64, 9.0_f64, 20.0_f64, 0.64_f64);
        // R1 chosen so the squared bracket cancels, chi2 = 0: K' = 0 exactly.
        let r1 = -(c.sqrt()) / ((1.0 - c) * k_bar / (k * nf)).sqrt();
        assert!(bussi_alpha2(k, k_bar, nf, c, r1, 0.0).abs() < 1e-12);
        let mut state = 4242_u64;
        for _ in 0..10_000 {
            assert!(sample_alpha2(&mut state, k, k_bar, nf, c) >= 0.0);
        }
    }

    #[test]
    fn alpha2_without_coupling_is_exactly_one() {
        let mut state = 11_u64;
        for _ in 0..100 {
            assert_eq!(sample_alpha2(&mut state, 3.0, 8.0, 30.0, 1.0), 1.0);
        }
    }

    /// One step from a fixed `K` has the analytic mean and variance of eq. A7:
    /// `E[K'] = cK + (1-c)K̄` and
    /// `Var[K'] = 2(1-c)²K̄²/N_f + 4c(1-c)K K̄/N_f`.
    #[test]
    fn one_step_moments_match_the_bussi_formulas() {
        let (k, k_bar, nf, c) = (25.0_f64, 15.0_f64, 30.0_f64, 0.8_f64);
        let mut state = 20_240_607_u64;
        let k_new: Vec<f64> = (0..400_000)
            .map(|_| k * sample_alpha2(&mut state, k, k_bar, nf, c))
            .collect();
        let (mean, var) = mean_and_variance(&k_new);
        let mean_exact = c * k + (1.0 - c) * k_bar;
        let var_exact =
            2.0 * (1.0 - c).powi(2) * k_bar * k_bar / nf + 4.0 * c * (1.0 - c) * k * k_bar / nf;
        assert!(
            (mean / mean_exact - 1.0).abs() < 2e-3,
            "E[K'] = {mean}, expected {mean_exact}"
        );
        assert!(
            (var / var_exact - 1.0).abs() < 0.03,
            "Var[K'] = {var}, expected {var_exact}"
        );
    }

    /// With `c = 0` the kinetic energy is resampled from the canonical
    /// distribution regardless of its old value: mean `K̄`, variance `2K̄²/N_f`.
    #[test]
    fn complete_resampling_is_canonical() {
        let (k_bar, nf) = (45.0_f64, 30.0_f64);
        for (seed, k_old) in [(31_u64, 1.0_f64), (32, 500.0)] {
            let mut state = seed;
            let k_new: Vec<f64> = (0..300_000)
                .map(|_| k_old * sample_alpha2(&mut state, k_old, k_bar, nf, 0.0))
                .collect();
            let (mean, var) = mean_and_variance(&k_new);
            assert!((mean / k_bar - 1.0).abs() < 4e-3, "mean = {mean}");
            assert!(
                (var / (2.0 * k_bar * k_bar / nf) - 1.0).abs() < 0.03,
                "variance = {var}"
            );
        }
    }

    #[test]
    fn csvr_reproduces_canonical_kinetic_energy_statistics() {
        // N_f = 300, dt/tau = 0.1: ~2500 independent samples per run, so the
        // standard errors are ~0.2 % on <K> and ~3 % on Var(K).
        let (n_atoms, target, dt, tau) = (100, 300.0, 0.01, 0.1);
        let nf = (3 * n_atoms) as f64;
        let k_bar = 0.5 * nf * target;
        for seed in [3_u64, 11, 42, 2024] {
            let mut atoms = ideal_gas(n_atoms, 600.0);
            let mut thermostat = CsvrThermostat::with_seed(tau, seed);
            let (mean, var) =
                kinetic_energy_statistics(&mut thermostat, &mut atoms, target, dt, 5_000, 45_000);
            assert!(
                (mean / k_bar - 1.0).abs() < 0.02,
                "seed {seed}: <T> = {} K, target {target} K",
                2.0 * mean / nf
            );
            let var_exact = 2.0 * k_bar * k_bar / nf;
            assert!(
                (var / var_exact - 1.0).abs() < 0.15,
                "seed {seed}: Var(K) = {var}, canonical value {var_exact} (ratio {})",
                var / var_exact
            );
        }
    }

    #[test]
    fn csvr_small_system_is_gamma_distributed() {
        // One atom: N_f = 3, K ~ Gamma(3/2, k_B T) => <K> = 1.5, Var(K) = 1.5.
        let target = 1.0;
        let mut atoms = ideal_gas(1, 1.0);
        let mut thermostat = CsvrThermostat::with_seed(0.02, 77);
        let (mean, var) =
            kinetic_energy_statistics(&mut thermostat, &mut atoms, target, 0.01, 2_000, 200_000);
        assert!((mean / 1.5 - 1.0).abs() < 0.03, "<K> = {mean}");
        assert!((var / 1.5 - 1.0).abs() < 0.10, "Var(K) = {var}");
    }

    #[test]
    fn svr_reproduces_canonical_kinetic_energy_statistics() {
        let (n_atoms, target, dt, tau) = (100, 300.0, 0.01, 0.1);
        let nf = 3 * n_atoms;
        let k_bar = 0.5 * nf as f64 * target;
        for seed in [5_u64, 9, 1234, 987_654] {
            let mut atoms = ideal_gas(n_atoms, 150.0);
            let mut thermostat = SvrThermostat::with_seed(tau, nf, seed);
            let (mean, var) =
                kinetic_energy_statistics(&mut thermostat, &mut atoms, target, dt, 5_000, 45_000);
            assert!(
                (mean / k_bar - 1.0).abs() < 0.02,
                "seed {seed}: <T> = {} K, target {target} K (the old SVR equilibrated near T/2)",
                2.0 * mean / nf as f64
            );
            let var_exact = 2.0 * k_bar * k_bar / nf as f64;
            assert!(
                (var / var_exact - 1.0).abs() < 0.15,
                "seed {seed}: Var(K) = {var}, canonical value {var_exact} (ratio {})",
                var / var_exact
            );
        }
    }

    #[test]
    fn svr_rescale_factor_one_step_moments() {
        // c = exp(-dt/tau) = 0.8 with dt = 1.
        let tau = -1.0 / 0.8_f64.ln();
        let mut svr = SvrThermostat::with_seed(tau, 30, 8_675_309);
        let (k, k_bar, c, nf) = (25.0_f64, 15.0_f64, 0.8_f64, 30.0_f64);
        let k_new: Vec<f64> = (0..400_000)
            .map(|_| {
                let s = svr.rescale_factor(k, k_bar, 1.0);
                s * s * k
            })
            .collect();
        let (mean, var) = mean_and_variance(&k_new);
        let mean_exact = c * k + (1.0 - c) * k_bar;
        let var_exact =
            2.0 * (1.0 - c).powi(2) * k_bar * k_bar / nf + 4.0 * c * (1.0 - c) * k * k_bar / nf;
        assert!((mean / mean_exact - 1.0).abs() < 2e-3, "E[K'] = {mean}");
        assert!((var / var_exact - 1.0).abs() < 0.03, "Var[K'] = {var}");
    }

    #[test]
    fn svr_rescale_factor_guards_degenerate_input() {
        let mut svr = SvrThermostat::with_seed(0.1, 30, 1);
        assert_eq!(svr.rescale_factor(0.0, 10.0, 0.01), 1.0);
        assert_eq!(svr.rescale_factor(-1.0, 10.0, 0.01), 1.0);
        assert_eq!(svr.rescale_factor(f64::NAN, 10.0, 0.01), 1.0);
        assert_eq!(svr.rescale_factor(10.0, -1.0, 0.01), 1.0);
        assert_eq!(svr.rescale_factor(10.0, f64::NAN, 0.01), 1.0);
        let mut empty = SvrThermostat::with_seed(0.1, 0, 1);
        assert_eq!(empty.rescale_factor(10.0, 10.0, 0.01), 1.0);
        // No elapsed time: no coupling, velocities untouched.
        assert_eq!(svr.rescale_factor(10.0, 20.0, 0.0), 1.0);
    }

    #[test]
    fn svr_cools_towards_zero_temperature_without_noise() {
        // Target K̄ = 0: the update collapses to pure friction, α² = c.
        let mut svr = SvrThermostat::with_seed(2.0, 30, 3);
        let dt = 0.5;
        let expected = (-dt / 2.0_f64).exp().sqrt();
        let s = svr.rescale_factor(12.0, 0.0, dt);
        assert!((s - expected).abs() < 1e-14, "s = {s}, expected {expected}");
    }

    #[test]
    fn csvr_apply_edge_cases_leave_the_system_alone() {
        let mut thermostat = CsvrThermostat::with_seed(0.1, 5);
        let mut empty = AtomSet::new();
        thermostat.apply(&mut empty, 300.0, 0.01, 1.0);
        assert!(empty.is_empty());

        let mut at_rest = ideal_gas(4, 0.0);
        thermostat.apply(&mut at_rest, 300.0, 0.01, 1.0);
        assert!(at_rest.velocities.iter().all(|v| v.norm() == 0.0));

        let mut atoms = ideal_gas(4, 100.0);
        let before = atoms.velocities.clone();
        thermostat.apply(&mut atoms, -5.0, 0.01, 1.0);
        assert_eq!(atoms.velocities, before, "negative target must be ignored");
        thermostat.apply(&mut atoms, 300.0, 0.0, 1.0);
        assert_eq!(atoms.velocities, before, "dt = 0 means no coupling");
    }
}
