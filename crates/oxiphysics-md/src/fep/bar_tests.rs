// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the BAR estimator [`Bar`].
//!
//! The statistical tests draw Gaussian energy differences from a seeded
//! generator.  For forward work `w = β ΔU ~ N(μ, s²)` sampled in state A, the
//! Crooks relation forces the reverse samples drawn in state B to follow
//! `N(μ − s², s²)`, and the exact free energy is `β ΔF = μ − s²/2`.

use super::*;
use crate::gaussian_noise::standard_normal;
use rand::SeedableRng;
use rand::rngs::StdRng;

/// BAR estimator in reduced units (`β = 1`): samples and `ΔF` are in kT.
fn reduced_bar(samples_a: Vec<f64>, samples_b: Vec<f64>) -> Bar {
    Bar {
        beta: 1.0,
        samples_a,
        samples_b,
    }
}

/// `n` draws from `N(mean, sigma²)`.
fn gaussian_samples(rng: &mut StdRng, n: usize, mean: f64, sigma: f64) -> Vec<f64> {
    (0..n)
        .map(|_| mean + sigma * standard_normal(rng))
        .collect()
}

/// Gaussian forward/reverse work samples with exact `β ΔF = mu − sigma²/2`.
fn gaussian_pair(
    rng: &mut StdRng,
    n_a: usize,
    n_b: usize,
    mu: f64,
    sigma: f64,
) -> (Vec<f64>, Vec<f64>) {
    (
        gaussian_samples(rng, n_a, mu, sigma),
        gaussian_samples(rng, n_b, mu - sigma * sigma, sigma),
    )
}

#[test]
fn constant_energy_shift_is_returned_for_equal_and_unequal_sample_sizes() {
    // Every sample of both states has ΔU = 2 kT: the only consistent free energy
    // is ΔF = 2 kT.  The old iteration (flipped Fermi arguments, no +C) moved by
    // -ΔU per step and never converged; with n_A != n_B a missing ln(n_A/n_B)
    // shows up as a bias of exactly that size.
    for (n_a, n_b) in [(50, 50), (30, 10), (10, 30), (1, 1), (1, 400), (400, 1)] {
        let mut bar = Bar::new(300.0);
        let du = 2.0 / bar.beta;
        for _ in 0..n_a {
            bar.add_sample_a(du);
        }
        for _ in 0..n_b {
            bar.add_sample_b(du);
        }
        let df = bar.delta_f_estimate();
        assert!(
            (df - du).abs() < 1e-9 * du,
            "n_A = {n_a}, n_B = {n_b}: ΔF = {df} J, expected {du} J"
        );
    }
}

#[test]
fn gaussian_work_with_unequal_sample_sizes_recovers_the_analytic_free_energy() {
    // σ = 1 kT, μ = 6 kT  ⇒  β ΔF = 5.5.  n_A = 2 n_B makes ln(n_A/n_B) = 0.69 kT
    // matter: dropping it biases the result by 12 %.
    let (mu, sigma) = (6.0, 1.0);
    let expected = mu - 0.5 * sigma * sigma;
    for seed in [3239_u64, 17] {
        let mut rng = StdRng::seed_from_u64(seed);
        let (w_a, w_b) = gaussian_pair(&mut rng, 4000, 2000, mu, sigma);
        let df = reduced_bar(w_a, w_b).delta_f_estimate();
        assert!(
            (df - expected).abs() < 0.02 * expected,
            "seed {seed}: β ΔF = {df}, analytic {expected}"
        );
    }
}

#[test]
fn gaussian_work_with_poor_overlap_still_converges_to_the_analytic_value() {
    // s = 3 dissipates s²/2 = 4.5 kT: the forward and reverse distributions only
    // overlap in their tails, where a fixed-point iteration contracts very slowly.
    let (mu, sigma) = (20.0, 3.0);
    let expected = mu - 0.5 * sigma * sigma;
    let mut rng = StdRng::seed_from_u64(2024);
    let (w_a, w_b) = gaussian_pair(&mut rng, 5000, 5000, mu, sigma);
    let df = reduced_bar(w_a, w_b).delta_f_estimate();
    assert!(
        (df - expected).abs() < 0.25,
        "β ΔF = {df}, analytic {expected}"
    );
}

#[test]
fn bar_beats_the_one_sided_zwanzig_estimators_on_the_same_data() {
    // With 4.5 kT of dissipation the forward Zwanzig estimate is dominated by a few
    // rare samples and is off by about half a kT on average, while BAR, which uses
    // both ends, is off by less than a tenth of a kT.  Averaged over 20 data sets.
    let (mu, sigma) = (20.0, 3.0);
    let expected = mu - 0.5 * sigma * sigma;
    let (mut bar_error, mut forward_error) = (0.0, 0.0);
    let n_sets = 20_u64;
    for seed in 0..n_sets {
        let mut rng = StdRng::seed_from_u64(7700 + seed);
        let (w_a, w_b) = gaussian_pair(&mut rng, 1000, 1000, mu, sigma);
        let forward_only = reduced_bar(w_a.clone(), Vec::new()).delta_f_estimate();
        forward_error += (forward_only - expected).abs();
        bar_error += (reduced_bar(w_a, w_b).delta_f_estimate() - expected).abs();
    }
    let (bar_error, forward_error) = (bar_error / n_sets as f64, forward_error / n_sets as f64);
    assert!(
        bar_error < 0.15 && forward_error > 2.5 * bar_error,
        "mean |error|: BAR {bar_error}, forward Zwanzig {forward_error}"
    );
}

#[test]
fn symmetric_samples_give_zero_free_energy() {
    let mut bar = Bar::new(300.0);
    let kt = 1.0 / bar.beta;
    for du in [1.0, -1.0, 2.0, -2.0] {
        bar.add_sample_a(du * kt);
        bar.add_sample_b(-du * kt);
    }
    let df = bar.delta_f_estimate();
    assert!(df.abs() < 1e-12 * kt, "ΔF = {df} J");
}

#[test]
fn joule_sized_energy_differences_do_not_overflow() {
    // ΔU = ±1 J, ±2 J gives β ΔU ~ 1e20 at 300 K, which overflows any direct
    // evaluation of exp(β ΔU); the log-space solver returns the symmetric answer.
    let mut bar = Bar::new(300.0);
    for du in [1.0_f64, -1.0, 2.0, -2.0] {
        bar.add_sample_a(du);
        bar.add_sample_b(-du);
    }
    let df = bar.delta_f_estimate();
    assert!(df.is_finite() && df.abs() < 1e-6, "ΔF = {df} J");
}

#[test]
fn free_energy_follows_a_common_shift_of_all_samples() {
    // ΔU → ΔU + d on both sides moves ΔF by exactly d (U_B → U_B + d).
    let mut rng = StdRng::seed_from_u64(5);
    let (w_a, w_b) = gaussian_pair(&mut rng, 300, 500, 4.0, 1.5);
    let base = reduced_bar(w_a.clone(), w_b.clone()).delta_f_estimate();
    for d in [-7.25_f64, 0.5, 31.0] {
        let shifted = reduced_bar(
            w_a.iter().map(|w| w + d).collect(),
            w_b.iter().map(|w| w + d).collect(),
        )
        .delta_f_estimate();
        assert!(
            (shifted - (base + d)).abs() < 1e-9 * (1.0 + d.abs()),
            "d = {d}: {shifted} vs {}",
            base + d
        );
    }
}

#[test]
fn exchanging_the_two_states_negates_the_free_energy() {
    // Relabelling A <-> B flips ΔU → −ΔU and swaps the sample sets.
    let mut rng = StdRng::seed_from_u64(11);
    let (w_a, w_b) = gaussian_pair(&mut rng, 700, 300, 3.0, 1.2);
    let forward = reduced_bar(w_a.clone(), w_b.clone()).delta_f_estimate();
    let swapped = reduced_bar(
        w_b.iter().map(|w| -w).collect(),
        w_a.iter().map(|w| -w).collect(),
    )
    .delta_f_estimate();
    assert!(
        (forward + swapped).abs() < 1e-9 * (1.0 + forward.abs()),
        "{forward} and {swapped} are not opposite"
    );
}

#[test]
fn sample_order_does_not_matter() {
    let mut rng = StdRng::seed_from_u64(23);
    let (w_a, w_b) = gaussian_pair(&mut rng, 200, 200, 2.0, 1.0);
    let forward = reduced_bar(w_a.clone(), w_b.clone()).delta_f_estimate();
    let reversed = reduced_bar(
        w_a.into_iter().rev().collect(),
        w_b.into_iter().rev().collect(),
    )
    .delta_f_estimate();
    assert!(
        (forward - reversed).abs() < 1e-11,
        "{forward} vs {reversed}"
    );
}

#[test]
fn one_sided_data_reduces_to_the_zwanzig_estimators() {
    let mut bar = Bar::new(300.0);
    let mut zwanzig = ExponentialAveraging::new(300.0);
    let kt = 1.0 / bar.beta;
    for i in 0..200 {
        let du = (0.3 + 0.01 * f64::from(i)) * kt;
        bar.add_sample_a(du);
        zwanzig.add_sample(du);
    }
    let forward = bar.delta_f_estimate();
    assert!(
        (forward - zwanzig.delta_f()).abs() < 1e-12 * kt,
        "{forward} vs {}",
        zwanzig.delta_f()
    );

    // State B only: ΔF = kT ln ⟨exp(β ΔU)⟩_B.
    let mut reverse_only = Bar::new(300.0);
    let w_b = [-1.5_f64, -0.5, -2.5, -1.0];
    for w in w_b {
        reverse_only.add_sample_b(w * kt);
    }
    let expected = kt * (w_b.iter().map(|w| w.exp()).sum::<f64>() / 4.0).ln();
    assert!((reverse_only.delta_f_estimate() - expected).abs() < 1e-12 * kt);
}

#[test]
fn one_sided_estimates_survive_energies_that_overflow_a_direct_exponential() {
    // exp(−β ΔU) underflows to zero for β ΔU = 2e3, which made the old fallback
    // return +∞; the log-sum-exp form returns the finite, exact answer.
    let mut bar = Bar::new(300.0);
    let kt = 1.0 / bar.beta;
    bar.add_sample_a(2.0e3 * kt);
    bar.add_sample_a(2.0e3 * kt);
    let df = bar.delta_f_estimate();
    assert!((df - 2.0e3 * kt).abs() < 1e-9 * kt * 2.0e3, "ΔF = {df}");
}

#[test]
fn empty_and_invalid_input() {
    assert_eq!(Bar::new(300.0).delta_f_estimate(), 0.0);

    let mut nan_bar = Bar::new(300.0);
    nan_bar.add_sample_a(1.0e-21);
    nan_bar.add_sample_a(f64::NAN);
    nan_bar.add_sample_b(2.0e-21);
    assert!(nan_bar.delta_f_estimate().is_nan());

    // One sample per side at ΔU = 0: no free-energy difference.
    let mut single = Bar::new(300.0);
    single.add_sample_a(0.0);
    single.add_sample_b(0.0);
    assert!(single.delta_f_estimate().abs() < 1e-30);
}

#[test]
fn bar_agrees_with_two_state_mbar() {
    // For two states MBAR reduces to BAR: both estimators must return the same
    // free-energy difference from the same samples, to solver precision.
    let kappa = [1.0_f64, 4.0];
    let mut rng = StdRng::seed_from_u64(1976);
    let mut mbar = Mbar::new(2, 300.0);
    let (mut w_a, mut w_b) = (Vec::new(), Vec::new());
    for (state, n) in [(0_usize, 60_usize), (1, 40)] {
        for _ in 0..n {
            let x = standard_normal(&mut rng) / kappa[state].sqrt();
            let u = [0.5 * kappa[0] * x * x, 0.5 * kappa[1] * x * x];
            mbar.add_sample(state, &u);
            let w = u[1] - u[0];
            if state == 0 {
                w_a.push(w);
            } else {
                w_b.push(w);
            }
        }
    }
    let f = mbar.estimate_free_energies(20_000);
    let bar = reduced_bar(w_a, w_b).delta_f_estimate();
    assert!(
        (f[1] - f[0] - bar).abs() < 1e-7,
        "MBAR {} vs BAR {bar}",
        f[1] - f[0]
    );
}
