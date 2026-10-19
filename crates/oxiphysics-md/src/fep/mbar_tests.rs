// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`Mbar`] and [`mbar_jackknife_uncertainty`].
//!
//! The reference system is a pair of one-dimensional harmonic states with
//! reduced potentials `u_k(x) = κ_k x² / 2`.  Their partition functions are
//! `Z_k = sqrt(2π / κ_k)`, so the exact reduced free-energy difference is
//! `f_1 − f_0 = ½ ln(κ_1 / κ_0)`, and a configuration drawn from state `k` is
//! `x ~ N(0, 1/κ_k)`.

use super::*;
use crate::gaussian_noise::standard_normal;
use rand::SeedableRng;
use rand::rngs::StdRng;

const KAPPA: [f64; 2] = [1.0, 4.0];

/// Exact `f_1 − f_0` of the two harmonic states.
fn exact_free_energy_difference() -> f64 {
    0.5 * (KAPPA[1] / KAPPA[0]).ln()
}

/// MBAR data set with `counts[k]` configurations drawn from harmonic state `k`.
fn harmonic_dataset(rng: &mut StdRng, counts: [usize; 2]) -> Mbar {
    let mut mbar = Mbar::new(2, 300.0);
    for (state, &count) in counts.iter().enumerate() {
        for _ in 0..count {
            let x = standard_normal(rng) / KAPPA[state].sqrt();
            let u = [0.5 * KAPPA[0] * x * x, 0.5 * KAPPA[1] * x * x];
            mbar.add_sample(state, &u);
        }
    }
    mbar
}

/// A copy of `mbar` built through the public API without sample `leave_out`.
fn without_sample(mbar: &Mbar, leave_out: usize) -> Mbar {
    let mut reduced = Mbar::new(mbar.n_states, 300.0);
    for sample in (0..mbar.n_samples()).filter(|&s| s != leave_out) {
        let column: Vec<f64> = mbar.u_kn.iter().map(|row| row[sample]).collect();
        reduced.add_sample(mbar.sampled_from[sample], &column);
    }
    reduced
}

/// Jackknife standard errors by brute force: fully converged cold-start solves
/// of every leave-one-out data set, then the textbook formula.
fn brute_force_jackknife(mbar: &Mbar) -> Vec<f64> {
    let n = mbar.n_samples();
    let replicates: Vec<Vec<f64>> = (0..n)
        .map(|s| without_sample(mbar, s).estimate_free_energies(5000))
        .collect();
    (0..mbar.n_states)
        .map(|k| {
            let mean = replicates.iter().map(|r| r[k]).sum::<f64>() / n as f64;
            let sum_sq: f64 = replicates.iter().map(|r| (r[k] - mean).powi(2)).sum();
            ((n as f64 - 1.0) / n as f64 * sum_sq).sqrt()
        })
        .collect()
}

#[test]
fn add_sample_records_the_state_of_every_sample() {
    let mut mbar = Mbar::new(3, 300.0);
    mbar.add_sample(2, &[0.1, 0.2, 0.3]);
    mbar.add_sample(0, &[0.4, 0.5, 0.6]);
    mbar.add_sample(2, &[0.7, 0.8, 0.9]);
    assert_eq!(mbar.sampled_from, vec![2, 0, 2]);
    assert_eq!(mbar.n_k, vec![1, 0, 2]);
    assert_eq!(mbar.n_samples(), 3);
    assert_eq!(mbar.u_kn[1], vec![0.2, 0.5, 0.8]);
    assert!(mbar.is_consistent());

    // Direct edits that break the bookkeeping are detected.
    let mut broken = Mbar::new(3, 300.0);
    broken.add_sample(2, &[0.1, 0.2, 0.3]);
    broken.sampled_from.clear();
    assert!(!broken.is_consistent());
    let mut wrong_count = Mbar::new(3, 300.0);
    wrong_count.add_sample(1, &[0.1, 0.2, 0.3]);
    wrong_count.n_k[1] = 5;
    assert!(!wrong_count.is_consistent());
    let mut bad_state = Mbar::new(3, 300.0);
    bad_state.add_sample(1, &[0.1, 0.2, 0.3]);
    bad_state.sampled_from[0] = 9;
    assert!(!bad_state.is_consistent());
}

#[test]
fn two_harmonic_states_give_the_exact_free_energy_difference() {
    let mut rng = StdRng::seed_from_u64(42);
    let mbar = harmonic_dataset(&mut rng, [600, 600]);
    let f = mbar.estimate_free_energies(10_000);
    assert_eq!(f[0], 0.0);
    assert!(
        (f[1] - exact_free_energy_difference()).abs() < 0.1,
        "f_1 = {}, exact {}",
        f[1],
        exact_free_energy_difference()
    );
}

#[test]
fn warm_start_from_the_solution_is_a_fixed_point() {
    let mut rng = StdRng::seed_from_u64(9);
    let mbar = harmonic_dataset(&mut rng, [80, 50]);
    let converged = mbar.estimate_free_energies(10_000);
    // One sweep from the converged answer leaves it in place ...
    let again = mbar.estimate_free_energies_from(&converged, 1);
    for (a, b) in converged.iter().zip(&again) {
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }
    // ... while one sweep from zero is still far away from it.
    let cold = mbar.estimate_free_energies(1);
    assert!(
        (cold[1] - converged[1]).abs() > 1e-3,
        "{cold:?} vs {converged:?}"
    );
    // A wrong start still converges to the same solution.
    let from_far = mbar.estimate_free_energies_from(&[3.0, -8.0], 10_000);
    assert!((from_far[1] - converged[1]).abs() < 1e-8);
    // Missing start entries count as zero.
    let from_short = mbar.estimate_free_energies_from(&[], 10_000);
    assert!((from_short[1] - converged[1]).abs() < 1e-8);
}

#[test]
fn states_without_samples_do_not_disturb_the_solution() {
    // State 2 has no samples (N_2 = 0): its free energy follows by reweighting and
    // the other two must not change.
    let mut rng = StdRng::seed_from_u64(31);
    let two = harmonic_dataset(&mut rng, [70, 70]);
    let mut three = Mbar::new(3, 300.0);
    for sample in 0..two.n_samples() {
        let u0 = two.u_kn[0][sample];
        let u1 = two.u_kn[1][sample];
        // State 2 is harmonic with kappa = 2: u_2 = x^2 = 2 u_0 (kappa_0 = 1).
        three.add_sample(two.sampled_from[sample], &[u0, u1, 2.0 * u0]);
    }
    let f2 = two.estimate_free_energies(10_000);
    let f3 = three.estimate_free_energies(10_000);
    assert!((f3[1] - f2[1]).abs() < 1e-8, "{} vs {}", f3[1], f2[1]);
    assert!(f3[2].is_finite());
}

#[test]
fn nan_reduced_potentials_are_reported_not_hidden() {
    let mut mbar = Mbar::new(2, 300.0);
    mbar.add_sample(0, &[0.1, 0.4]);
    mbar.add_sample(1, &[f64::NAN, 0.2]);
    let f = mbar.estimate_free_energies(100);
    assert!(f.iter().all(|v| v.is_nan()), "{f:?}");
}

#[test]
fn jackknife_equals_the_brute_force_leave_one_out_definition() {
    // Unequal state populations, so that returning a sample to the wrong state
    // (or forgetting to decrement N_k) changes the answer.
    let mut rng = StdRng::seed_from_u64(3242);
    let mbar = harmonic_dataset(&mut rng, [14, 9]);
    let expected = brute_force_jackknife(&mbar);
    let sigma = mbar_jackknife_uncertainty(&mbar, 5000);
    assert_eq!(sigma.len(), 2);
    assert_eq!(sigma[0], 0.0, "state 0 is the reference");
    assert!(expected[1] > 0.05, "degenerate reference {expected:?}");
    assert!(
        (sigma[1] - expected[1]).abs() < 1e-6,
        "jackknife {} vs brute force {}",
        sigma[1],
        expected[1]
    );
}

#[test]
fn jackknife_depends_on_the_data_not_just_on_their_number() {
    // The same number of samples from the same states but different energies must
    // give a different uncertainty (the old code returned a function of the first
    // sample only).
    let mut rng = StdRng::seed_from_u64(1);
    let first = harmonic_dataset(&mut rng, [40, 40]);
    let second = harmonic_dataset(&mut rng, [40, 40]);
    let sigma_first = mbar_jackknife_uncertainty(&first, 5000)[1];
    let sigma_second = mbar_jackknife_uncertainty(&second, 5000)[1];
    assert!(sigma_first > 0.0 && sigma_second > 0.0);
    assert!(
        (sigma_first - sigma_second).abs() > 1e-4,
        "{sigma_first} vs {sigma_second}"
    );
}

#[test]
fn jackknife_error_shrinks_like_one_over_root_n() {
    let mut rng = StdRng::seed_from_u64(2025);
    let small = harmonic_dataset(&mut rng, [100, 100]);
    let large = harmonic_dataset(&mut rng, [400, 400]);
    let sigma_small = mbar_jackknife_uncertainty(&small, 2000)[1];
    let sigma_large = mbar_jackknife_uncertainty(&large, 2000)[1];
    let ratio = sigma_small / sigma_large;
    assert!(
        ratio > 1.5 && ratio < 2.7,
        "σ(100) = {sigma_small}, σ(400) = {sigma_large}, ratio {ratio} (expected ~2)"
    );
}

#[test]
fn jackknife_error_matches_the_spread_over_independent_data_sets() {
    // 100 independent data sets of the two harmonic states: the standard deviation of
    // the MBAR estimates is the true statistical error.  The jackknife error of a
    // single data set (averaged over a few of them to tame its own noise) must
    // reproduce it to within 30 %.
    let n_replicates = 100;
    let n_jackknife = 8;
    for seed in [3243_u64, 8] {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut estimates = Vec::new();
        let mut jackknife = Vec::new();
        for replicate in 0..n_replicates {
            let mbar = harmonic_dataset(&mut rng, [100, 100]);
            estimates.push(mbar.estimate_free_energies(2000)[1]);
            if replicate < n_jackknife {
                jackknife.push(mbar_jackknife_uncertainty(&mbar, 2000)[1]);
            }
        }
        let n = n_replicates as f64;
        let mean = estimates.iter().sum::<f64>() / n;
        let spread = (estimates.iter().map(|e| (e - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        let jackknife_mean = jackknife.iter().sum::<f64>() / jackknife.len() as f64;
        assert!(
            (jackknife_mean / spread - 1.0).abs() < 0.30,
            "seed {seed}: jackknife σ = {jackknife_mean}, spread of {n_replicates} data sets = {spread}"
        );
        // The estimates scatter around the exact value (3 standard errors of the mean).
        assert!(
            (mean - exact_free_energy_difference()).abs() < 3.0 * spread / n.sqrt() + 0.01,
            "seed {seed}: mean estimate {mean}, exact {}",
            exact_free_energy_difference()
        );
    }
}

#[test]
fn identical_states_have_zero_uncertainty() {
    let mut mbar = Mbar::new(3, 300.0);
    for i in 0..12 {
        let u = 0.1 * f64::from(i);
        mbar.add_sample(i as usize % 3, &[u, u, u]);
    }
    let sigma = mbar_jackknife_uncertainty(&mbar, 500);
    assert_eq!(sigma.len(), 3);
    for (k, s) in sigma.iter().enumerate() {
        assert!(s.abs() < 1e-9, "σ[{k}] = {s}");
    }
}

#[test]
fn jackknife_is_undefined_without_enough_or_consistent_data() {
    // No state at all.
    assert!(mbar_jackknife_uncertainty(&Mbar::new(0, 300.0), 10).is_empty());
    // No sample, one sample: no leave-one-out replicates exist.
    let empty = Mbar::new(2, 300.0);
    assert!(
        mbar_jackknife_uncertainty(&empty, 10)
            .iter()
            .all(|s| s.is_nan())
    );
    let mut single = Mbar::new(2, 300.0);
    single.add_sample(0, &[0.0, 1.0]);
    assert!(
        mbar_jackknife_uncertainty(&single, 10)
            .iter()
            .all(|s| s.is_nan())
    );
    // Without a single self-consistency sweep there is no spread to measure.
    let mut rng = StdRng::seed_from_u64(4);
    let enough = harmonic_dataset(&mut rng, [10, 10]);
    assert!(
        mbar_jackknife_uncertainty(&enough, 0)
            .iter()
            .all(|s| s.is_nan())
    );
    // Bookkeeping that disagrees with the data must not index out of range.
    let mut mbar = harmonic_dataset(&mut rng, [10, 10]);
    mbar.sampled_from.pop();
    assert!(
        mbar_jackknife_uncertainty(&mbar, 10)
            .iter()
            .all(|s| s.is_nan())
    );
    let mut out_of_range = harmonic_dataset(&mut rng, [10, 10]);
    out_of_range.sampled_from[3] = 7;
    assert!(
        mbar_jackknife_uncertainty(&out_of_range, 10)
            .iter()
            .all(|s| s.is_nan())
    );
}

#[test]
fn samples_from_a_single_state_still_give_finite_errors() {
    // All configurations come from state 0; the other state is reweighted.
    let mut rng = StdRng::seed_from_u64(12);
    let mbar = harmonic_dataset(&mut rng, [60, 0]);
    let sigma = mbar_jackknife_uncertainty(&mbar, 5000);
    assert_eq!(sigma[0], 0.0);
    assert!(sigma[1].is_finite() && sigma[1] > 0.0, "{sigma:?}");
}

// ─────────────────────────────────────────────────────────────────────────────
// Analytic (Shirts–Chodera) uncertainty
// ─────────────────────────────────────────────────────────────────────────────

/// MBAR data set of one-dimensional harmonic states with force constants
/// `kappas`; `counts[k]` configurations are drawn from state `k`.
fn harmonic_dataset_k(rng: &mut StdRng, kappas: &[f64], counts: &[usize]) -> Mbar {
    let mut mbar = Mbar::new(kappas.len(), 300.0);
    for (state, &count) in counts.iter().enumerate() {
        for _ in 0..count {
            let x = standard_normal(rng) / kappas[state].sqrt();
            let u: Vec<f64> = kappas.iter().map(|kappa| 0.5 * kappa * x * x).collect();
            mbar.add_sample(state, &u);
        }
    }
    mbar
}

/// Sample standard deviation of `values`.
fn sample_std(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
}

#[test]
fn analytic_uncertainty_matches_the_spread_over_independent_data_sets() {
    // The standard deviation of the MBAR estimates over 100 independent data sets
    // is the true statistical error; the asymptotic formula, averaged over the
    // data sets to tame its own noise, must reproduce it to within 30 %.
    let n_replicates = 100;
    for seed in [3244_u64, 11] {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut estimates = Vec::new();
        let mut analytic = Vec::new();
        for _ in 0..n_replicates {
            let mbar = harmonic_dataset(&mut rng, [100, 100]);
            estimates.push(mbar.estimate_free_energies(2000)[1]);
            let sigma = mbar_analytic_uncertainty(&mbar, 2000).expect("well-formed data set");
            assert_eq!(sigma.len(), 2);
            assert_eq!(sigma[0], 0.0, "state 0 is the reference");
            analytic.push(sigma[1]);
        }
        let spread = sample_std(&estimates);
        let analytic_mean = analytic.iter().sum::<f64>() / n_replicates as f64;
        assert!(
            (analytic_mean / spread - 1.0).abs() < 0.30,
            "seed {seed}: analytic σ = {analytic_mean}, spread of {n_replicates} data sets = {spread}"
        );
    }
}

#[test]
fn analytic_uncertainty_of_three_states_matches_the_spread() {
    // Three harmonic states (κ = 1, 2, 4) with 50 configurations each: the K × K
    // pseudo-inverse of the covariance formula has a genuine two-dimensional
    // range here, not just the single difference of the two-state case.
    let kappas = [1.0_f64, 2.0, 4.0];
    let n_replicates = 100;
    let mut rng = StdRng::seed_from_u64(3246);
    let mut estimates = vec![Vec::new(); 3];
    let mut analytic = vec![Vec::new(); 3];
    for _ in 0..n_replicates {
        let mbar = harmonic_dataset_k(&mut rng, &kappas, &[50, 50, 50]);
        let f = mbar.estimate_free_energies(2000);
        let sigma = mbar_analytic_uncertainty(&mbar, 2000).expect("well-formed data set");
        assert_eq!(sigma.len(), 3);
        for (k, (estimate, error)) in estimates.iter_mut().zip(analytic.iter_mut()).enumerate() {
            estimate.push(f[k]);
            error.push(sigma[k]);
        }
    }
    assert!(
        analytic[0].iter().all(|&s| s == 0.0),
        "state 0 is the reference"
    );
    for k in 1..3 {
        let spread = sample_std(&estimates[k]);
        let analytic_mean = analytic[k].iter().sum::<f64>() / n_replicates as f64;
        assert!(
            (analytic_mean / spread - 1.0).abs() < 0.30,
            "state {k}: analytic σ = {analytic_mean}, spread of {n_replicates} data sets = {spread}"
        );
        // The estimates scatter around the exact value ½ ln(κ_k / κ_0).
        let mean = estimates[k].iter().sum::<f64>() / n_replicates as f64;
        let exact = 0.5 * (kappas[k] / kappas[0]).ln();
        assert!(
            (mean - exact).abs() < 3.0 * spread / (n_replicates as f64).sqrt() + 0.01,
            "state {k}: mean estimate {mean}, exact {exact}"
        );
    }
}

#[test]
fn analytic_and_jackknife_errors_agree_on_one_data_set() {
    let mut rng = StdRng::seed_from_u64(3245);
    let mbar = harmonic_dataset(&mut rng, [300, 300]);
    let analytic = mbar_analytic_uncertainty(&mbar, 5000).expect("well-formed data set");
    let jackknife = mbar_jackknife_uncertainty(&mbar, 5000);
    assert_eq!(analytic[0], 0.0);
    assert!(analytic[1] > 0.0 && jackknife[1] > 0.0);
    // Two independent routes to the same asymptotic error: the covariance formula
    // and the resampled leave-one-out spread agree to about 1 % at this size.
    assert!(
        (analytic[1] / jackknife[1] - 1.0).abs() < 0.05,
        "analytic σ = {}, jackknife σ = {}",
        analytic[1],
        jackknife[1]
    );
}

#[test]
fn analytic_uncertainty_shrinks_like_one_over_root_n() {
    // The error of a single data set fluctuates by ~8 %, so compare the averages
    // over 20 data sets of each size: quadrupling the samples halves the error.
    let mut rng = StdRng::seed_from_u64(2026);
    let mut mean_sigma = |counts: [usize; 2]| {
        let total: f64 = (0..20)
            .map(|_| {
                let mbar = harmonic_dataset(&mut rng, counts);
                mbar_analytic_uncertainty(&mbar, 2000).expect("well-formed data set")[1]
            })
            .sum();
        total / 20.0
    };
    let sigma_small = mean_sigma([100, 100]);
    let sigma_large = mean_sigma([400, 400]);
    let ratio = sigma_small / sigma_large;
    assert!(
        ratio > 1.8 && ratio < 2.2,
        "σ(100) = {sigma_small}, σ(400) = {sigma_large}, ratio {ratio} (expected ~2)"
    );
}

#[test]
fn identical_states_have_zero_analytic_uncertainty() {
    let mut mbar = Mbar::new(3, 300.0);
    for i in 0..12 {
        let u = 0.1 * f64::from(i);
        mbar.add_sample(i as usize % 3, &[u, u, u]);
    }
    let sigma = mbar_analytic_uncertainty(&mbar, 500).expect("well-formed data set");
    assert_eq!(sigma.len(), 3);
    for (k, s) in sigma.iter().enumerate() {
        assert!(s.abs() < 1e-6, "σ[{k}] = {s}");
    }
}

#[test]
fn analytic_uncertainty_reports_undefined_input_as_errors() {
    // No state at all: nothing to report.
    let no_states = mbar_analytic_uncertainty(&Mbar::new(0, 300.0), 10).expect("no states");
    assert!(no_states.is_empty());
    // States without a sample, and a request for no solver sweep at all.
    assert!(mbar_analytic_uncertainty(&Mbar::new(2, 300.0), 10).is_err());
    let mut rng = StdRng::seed_from_u64(4);
    let enough = harmonic_dataset(&mut rng, [10, 10]);
    assert!(mbar_analytic_uncertainty(&enough, 0).is_err());
    // Bookkeeping that disagrees with the data.
    let mut broken = harmonic_dataset(&mut rng, [10, 10]);
    broken.sampled_from.pop();
    assert!(mbar_analytic_uncertainty(&broken, 10).is_err());
    // A NaN reduced potential leaves the free energies undefined.
    let mut nan = Mbar::new(2, 300.0);
    nan.add_sample(0, &[0.1, 0.4]);
    nan.add_sample(1, &[f64::NAN, 0.2]);
    assert!(mbar_analytic_uncertainty(&nan, 100).is_err());
    // A single state is its own reference.
    let mut single = Mbar::new(1, 300.0);
    single.add_sample(0, &[0.3]);
    single.add_sample(0, &[0.9]);
    let sigma = mbar_analytic_uncertainty(&single, 10).expect("single state");
    assert_eq!(sigma, vec![0.0]);
}
