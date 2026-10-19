// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Bennett acceptance ratio (BAR) and multistate BAR
//! (MBAR) free-energy estimators.
//!
//! * `bar_iteration` evaluated the Fermi weights at the shift `ΔA + M` but
//!   returned `ΔA + k_B T ln(⟨f⟩_B/⟨f⟩_A)` without the shift, so for unequal sample
//!   sizes the fixed point was biased by `−k_B T ln(n_A/n_B)`.  `bar_uncertainty`
//!   used a variance formula without the `⟨f⟩²` normalisation.
//! * `mbar_free_energies` revisited the pooled samples `0..n_j` for every state
//!   `j` (never using the later ones), put `f_j + ln N_j` instead of `1` in the
//!   numerator and recomputed the denominators per state.
//!
//! All random samples come from a seeded generator, so the statistical checks are
//! deterministic.

use super::*;
use rand::SeedableRng;
use rand::rngs::StdRng;

/// Temperature of the tests (K).
const TEMP: f64 = 300.0;

/// Standard normal variate (Box–Muller, cosine branch).
fn standard_normal(rng: &mut StdRng) -> f64 {
    let u1: f64 = 1.0 - rng.random::<f64>();
    let u2: f64 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// `n` samples of `N(mean, std²)`.
fn gaussian(rng: &mut StdRng, mean: f64, std: f64, n: usize) -> Vec<f64> {
    (0..n).map(|_| mean + std * standard_normal(rng)).collect()
}

/// Forward and reverse energy differences with Gaussian statistics that obey
/// Crooks' theorem: `ΔU_AB ~ N(μ, σ²)` and `ΔU_BA ~ N(−μ + βσ², σ²)`, for which
/// `ΔA = μ − βσ²/2`.
fn crooks_gaussian_samples(
    rng: &mut StdRng,
    mu: f64,
    sigma: f64,
    n_a: usize,
    n_b: usize,
) -> (Vec<f64>, Vec<f64>) {
    let beta = beta_mol_kjmol(TEMP);
    let du_ab = gaussian(rng, mu, sigma, n_a);
    let du_ba = gaussian(rng, -mu + beta * sigma * sigma, sigma, n_b);
    (du_ab, du_ba)
}

/// Pooled samples of one-dimensional harmonic states `U_k(x) = ½ κ_k x²`
/// (kJ mol⁻¹ nm⁻²) drawn exactly from their Boltzmann distributions.
struct HarmonicEnsemble {
    /// `energies[k][n]`: energy of the pooled sample `n` in state `k` (kJ mol⁻¹).
    energies: Vec<Vec<f64>>,
    /// Number of samples drawn from each state.
    n_k: Vec<usize>,
    /// Pooled coordinates, the samples of state 0 first.
    samples: Vec<f64>,
}

impl HarmonicEnsemble {
    fn draw(rng: &mut StdRng, force_constants: &[f64], counts: &[usize]) -> Self {
        let beta = beta_mol_kjmol(TEMP);
        let mut samples = Vec::new();
        for (&kappa, &count) in force_constants.iter().zip(counts.iter()) {
            samples.extend(gaussian(rng, 0.0, (1.0 / (beta * kappa)).sqrt(), count));
        }
        let energies = force_constants
            .iter()
            .map(|&kappa| samples.iter().map(|&x| 0.5 * kappa * x * x).collect())
            .collect();
        Self {
            energies,
            n_k: counts.to_vec(),
            samples,
        }
    }

    /// Exact reduced free energy differences `f_k − f_0 = ½ ln(κ_k / κ_0)`.
    fn exact_reduced_free_energies(force_constants: &[f64]) -> Vec<f64> {
        let kappa_0 = force_constants[0];
        force_constants
            .iter()
            .map(|&kappa| 0.5 * (kappa / kappa_0).ln())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Scalar helpers
// ---------------------------------------------------------------------------

#[test]
fn log_fermi_matches_the_direct_evaluation_and_saturates() {
    for &x in &[-30.0, -2.5, -0.1, 0.0, 0.1, 2.5, 30.0] {
        let direct = fermi(x).ln();
        assert!(
            (log_fermi(x) - direct).abs() < 1.0e-12,
            "x = {x}: {} vs {direct}",
            log_fermi(x)
        );
    }
    assert_eq!(log_fermi(f64::INFINITY), f64::NEG_INFINITY);
    assert_eq!(log_fermi(f64::NEG_INFINITY), 0.0);
    // The naive 1/(1+exp(x)) is exactly zero here.
    assert!((log_fermi(1000.0) + 1000.0).abs() < 1.0e-12);
    assert!(log_fermi(f64::NAN).is_nan());
}

#[test]
fn log_sum_exp_handles_degenerate_inputs() {
    assert_eq!(log_sum_exp(&[]), f64::NEG_INFINITY);
    assert_eq!(
        log_sum_exp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]),
        f64::NEG_INFINITY
    );
    assert!((log_sum_exp(&[f64::NEG_INFINITY, 2.0]) - 2.0).abs() < 1.0e-15);
    assert_eq!(log_sum_exp(&[1.0, f64::INFINITY]), f64::INFINITY);
    assert!(log_sum_exp(&[1.0, f64::NAN]).is_nan());
    assert!(log_sum_exp(&[f64::NAN]).is_nan());
    // ln(e + e²) = 2 + ln(1 + e⁻¹), also far outside the range of exp.
    let expected = 2.0 + (1.0 + (-1.0_f64).exp()).ln();
    assert!((log_sum_exp(&[1.0, 2.0]) - expected).abs() < 1.0e-14);
    assert!((log_sum_exp(&[1000.0, 1001.0]) - (1001.0 + (-1.0_f64).exp().ln_1p())).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// BAR
// ---------------------------------------------------------------------------

#[test]
fn bar_recovers_the_crooks_free_energy_for_unequal_sample_counts() {
    // nA = 40000 and nB = 5000: the former implementation was biased by
    // k_B T ln(nA/nB) = 5.2 kJ/mol here.
    let (mu, sigma) = (6.0, 1.5);
    let mut rng = StdRng::seed_from_u64(2024);
    let (du_ab, du_ba) = crooks_gaussian_samples(&mut rng, mu, sigma, 40_000, 5_000);
    let expected = mu - 0.5 * beta_mol_kjmol(TEMP) * sigma * sigma;
    let (da, iterations) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-10, 1000);
    assert!(iterations < 1000, "BAR did not converge");
    assert!(
        (da - expected).abs() < 0.1,
        "BAR gave {da}, the exact value is {expected}"
    );
    let sigma_da = bar_uncertainty(&du_ab, &du_ba, da, TEMP);
    assert!(
        (da - expected).abs() < 5.0 * sigma_da,
        "error {} exceeds five times the predicted sigma {sigma_da}",
        (da - expected).abs()
    );
}

#[test]
fn bar_recovers_the_crooks_free_energy_when_the_reverse_set_is_larger() {
    let (mu, sigma) = (-3.0, 2.0);
    let mut rng = StdRng::seed_from_u64(77);
    let (du_ab, du_ba) = crooks_gaussian_samples(&mut rng, mu, sigma, 3_000, 24_000);
    let expected = mu - 0.5 * beta_mol_kjmol(TEMP) * sigma * sigma;
    let (da, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-10, 1000);
    assert!(
        (da - expected).abs() < 0.1,
        "BAR gave {da}, the exact value is {expected}"
    );
}

#[test]
fn bar_recovers_a_constant_shift_with_unbalanced_sample_counts() {
    // U_B = U_A + 1 kJ/mol everywhere: ΔU_AB = +1 and ΔU_BA = −1 whatever the
    // configuration, so ΔA = 1 exactly, with any sample sizes.
    for &(n_a, n_b) in &[(100_usize, 50_usize), (50, 100), (10, 1000), (700, 700)] {
        let du_ab = vec![1.0; n_a];
        let du_ba = vec![-1.0; n_b];
        let (da, iterations) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-12, 10_000);
        assert!(iterations < 10_000, "no convergence for ({n_a}, {n_b})");
        assert!(
            (da - 1.0).abs() < 1.0e-9,
            "n_A = {n_a}, n_B = {n_b}: ΔA = {da}, expected 1"
        );
    }
}

#[test]
fn bar_iteration_leaves_the_exact_solution_unchanged() {
    let (n_a, n_b) = (400_usize, 100_usize);
    let du_ab = vec![2.5; n_a];
    let du_ba = vec![-2.5; n_b];
    let m = (n_a as f64 / n_b as f64).ln() / beta_mol_kjmol(TEMP);
    let da = bar_iteration(&du_ab, &du_ba, 2.5, TEMP, m);
    assert!((da - 2.5).abs() < 1.0e-12, "fixed point moved to {da}");
    // Away from the solution the iteration moves towards it.
    let below = bar_iteration(&du_ab, &du_ba, 0.0, TEMP, m);
    assert!(below > 0.0 && below <= 2.5 + 1.0e-12, "step to {below}");
}

#[test]
fn bar_is_antisymmetric_under_exchange_of_the_states() {
    let mut rng = StdRng::seed_from_u64(5);
    let (du_ab, du_ba) = crooks_gaussian_samples(&mut rng, 4.0, 1.5, 3_000, 500);
    let (forward, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-13, 2000);
    let (backward, _) = bar_converge(&du_ba, &du_ab, TEMP, 1.0e-13, 2000);
    assert!(
        (forward + backward).abs() < 1.0e-8,
        "ΔA(A→B) = {forward}, ΔA(B→A) = {backward}"
    );
}

#[test]
fn bar_handles_energy_differences_of_thousands_of_kt() {
    // βΔU = 400: the Fermi weight f(400) = 4e-175 and f(1e5 · β) underflow a plain
    // evaluation, which then stalled at the starting estimate.
    let du_ab = vec![1000.0; 20];
    let du_ba = vec![-1000.0; 20];
    let (da, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-9, 1000);
    assert!((da - 1000.0).abs() < 1.0e-6, "ΔA = {da}");
    let (huge, _) = bar_converge(&[5.0e6; 8], &[-5.0e6; 8], TEMP, 1.0e-6, 1000);
    assert!((huge - 5.0e6).abs() < 1.0e-3, "ΔA = {huge}");
}

#[test]
fn bar_iteration_keeps_the_estimate_when_a_side_has_no_weight() {
    let du_ab = [f64::INFINITY, f64::INFINITY];
    let du_ba = [0.0, 0.0];
    assert_eq!(bar_iteration(&du_ab, &du_ba, 1.25, TEMP, 0.0), 1.25);
    // Empty sample sets also leave the estimate alone.
    assert_eq!(bar_iteration(&[], &du_ba, 1.25, TEMP, 0.0), 1.25);
    assert_eq!(bar_iteration(&du_ab, &[], 1.25, TEMP, 0.0), 1.25);
    assert_eq!(bar_converge(&[], &du_ba, TEMP, 1.0e-9, 10), (0.0, 0));
}

#[test]
fn bar_uncertainty_follows_the_normalised_variance_formula() {
    let mut rng = StdRng::seed_from_u64(31);
    let (du_ab, du_ba) = crooks_gaussian_samples(&mut rng, 3.0, 2.0, 600, 150);
    let (da, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-12, 2000);

    // Direct evaluation of σ² = (1/β²)[(⟨f²⟩/⟨f⟩² − 1)/n_A + (⟨f²⟩/⟨f⟩² − 1)/n_B].
    let beta = beta_mol_kjmol(TEMP);
    let m = (du_ab.len() as f64 / du_ba.len() as f64).ln() / beta;
    let weights_a: Vec<f64> = du_ab
        .iter()
        .map(|&du| fermi(beta * (du - da + m)))
        .collect();
    let weights_b: Vec<f64> = du_ba
        .iter()
        .map(|&du| fermi(beta * (du + da - m)))
        .collect();
    let mean = |values: &[f64]| values.iter().sum::<f64>() / values.len() as f64;
    let relative_variance = |values: &[f64]| {
        let squares: Vec<f64> = values.iter().map(|v| v * v).collect();
        mean(&squares) / mean(values).powi(2) - 1.0
    };
    let expected = ((relative_variance(&weights_a) / du_ab.len() as f64
        + relative_variance(&weights_b) / du_ba.len() as f64)
        .sqrt())
        / beta;
    let sigma = bar_uncertainty(&du_ab, &du_ba, da, TEMP);
    assert!(sigma > 0.0 && sigma.is_finite());
    assert!(
        (sigma - expected).abs() < 1.0e-9 * expected,
        "σ = {sigma}, direct evaluation {expected}"
    );
}

#[test]
fn bar_uncertainty_vanishes_for_constant_weights() {
    // A pure shift of the potential: every Fermi weight of a side is identical.
    let du_ab = vec![1.0; 200];
    let du_ba = vec![-1.0; 50];
    let (da, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-12, 10_000);
    let sigma = bar_uncertainty(&du_ab, &du_ba, da, TEMP);
    assert!(sigma < 1.0e-6, "σ = {sigma}");
}

#[test]
fn bar_uncertainty_degenerate_inputs() {
    assert_eq!(bar_uncertainty(&[], &[1.0], 0.0, TEMP), 0.0);
    assert_eq!(bar_uncertainty(&[1.0], &[], 0.0, TEMP), 0.0);
    // Zero weight on one side: the uncertainty is undefined.
    assert!(bar_uncertainty(&[f64::INFINITY; 3], &[0.0; 3], 0.0, TEMP).is_nan());
    // A single sample per side has no spread.
    assert!(bar_uncertainty(&[1.0], &[-1.0], 1.0, TEMP) < 1.0e-6);
}

#[test]
fn bar_uncertainty_predicts_the_scatter_between_replicates() {
    let (mu, sigma) = (3.0, 2.0);
    let (n_a, n_b) = (400_usize, 100_usize);
    let mut rng = StdRng::seed_from_u64(99);
    let mut estimates = Vec::new();
    let mut predicted = Vec::new();
    for _ in 0..200 {
        let (du_ab, du_ba) = crooks_gaussian_samples(&mut rng, mu, sigma, n_a, n_b);
        let (da, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-9, 500);
        estimates.push(da);
        predicted.push(bar_uncertainty(&du_ab, &du_ba, da, TEMP));
    }
    let n = estimates.len() as f64;
    let mean = estimates.iter().sum::<f64>() / n;
    let scatter = (estimates.iter().map(|e| (e - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
    let mean_predicted = predicted.iter().sum::<f64>() / n;
    let exact = mu - 0.5 * beta_mol_kjmol(TEMP) * sigma * sigma;
    assert!(
        (mean - exact).abs() < 3.0 * scatter / n.sqrt() + 0.02,
        "replicate mean {mean}, exact {exact}"
    );
    let ratio = mean_predicted / scatter;
    assert!(
        (0.8..1.25).contains(&ratio),
        "predicted σ = {mean_predicted}, observed scatter {scatter} (ratio {ratio})"
    );
}

// ---------------------------------------------------------------------------
// MBAR
// ---------------------------------------------------------------------------

#[test]
fn mbar_recovers_the_free_energy_of_two_harmonic_states() {
    // κ₁ = 4 κ₀: f₁ − f₀ = ½ ln 4.
    let force_constants = [2.0, 8.0];
    let mut rng = StdRng::seed_from_u64(11);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &force_constants, &[20_000, 5_000]);
    let f = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 20_000, 1.0e-11);
    let exact = HarmonicEnsemble::exact_reduced_free_energies(&force_constants);
    assert_eq!(f.len(), 2);
    assert_eq!(f[0], 0.0);
    assert!(
        (f[1] - exact[1]).abs() < 0.05,
        "f₁ − f₀ = {}, exact {}",
        f[1],
        exact[1]
    );
}

#[test]
fn mbar_recovers_the_free_energies_of_three_states_with_unequal_sample_counts() {
    let force_constants = [2.0, 4.0, 8.0];
    let mut rng = StdRng::seed_from_u64(12);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &force_constants, &[12_000, 6_000, 3_000]);
    let f = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 20_000, 1.0e-11);
    let exact = HarmonicEnsemble::exact_reduced_free_energies(&force_constants);
    assert_eq!(f.len(), 3);
    for k in 0..3 {
        assert!(
            (f[k] - exact[k]).abs() < 0.05,
            "state {k}: f = {}, exact {}",
            f[k],
            exact[k]
        );
    }
}

#[test]
fn mbar_estimates_a_state_without_samples_by_reweighting() {
    // Only the outer states are sampled; the middle one is reached through the
    // pooled configurations of its neighbours.
    let force_constants = [2.0, 4.0, 8.0];
    let mut rng = StdRng::seed_from_u64(13);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &force_constants, &[10_000, 0, 10_000]);
    assert_eq!(ensemble.n_k, vec![10_000, 0, 10_000]);
    let f = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 20_000, 1.0e-11);
    let exact = HarmonicEnsemble::exact_reduced_free_energies(&force_constants);
    for k in 0..3 {
        assert!(
            (f[k] - exact[k]).abs() < 0.05,
            "state {k}: f = {}, exact {}",
            f[k],
            exact[k]
        );
    }
}

#[test]
fn two_state_mbar_equals_bar_on_the_same_samples() {
    // For two states MBAR is BAR (Shirts & Chodera 2008), so the two estimators must
    // agree to the convergence tolerance with unequal sample sizes.
    let force_constants = [2.0, 8.0];
    let (n_a, n_b) = (3_000_usize, 800_usize);
    let mut rng = StdRng::seed_from_u64(21);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &force_constants, &[n_a, n_b]);
    let f = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 50_000, 1.0e-13);
    let mbar_da = mbar_delta_a(&f, TEMP)[1];

    // ΔU_AB = U₁ − U₀ on the samples from state 0, ΔU_BA = U₀ − U₁ on those from state 1.
    let du: Vec<f64> = ensemble.energies[1]
        .iter()
        .zip(ensemble.energies[0].iter())
        .map(|(u1, u0)| u1 - u0)
        .collect();
    let du_ab: Vec<f64> = du[..n_a].to_vec();
    let du_ba: Vec<f64> = du[n_a..].iter().map(|d| -d).collect();
    let (bar_da, _) = bar_converge(&du_ab, &du_ba, TEMP, 1.0e-13, 50_000);

    assert!(
        (mbar_da - bar_da).abs() < 1.0e-6,
        "MBAR ΔA = {mbar_da}, BAR ΔA = {bar_da}"
    );
    // And both are close to the exact value ½ k_B T ln 4.
    let exact =
        HarmonicEnsemble::exact_reduced_free_energies(&force_constants)[1] / beta_mol_kjmol(TEMP);
    assert!(
        (mbar_da - exact).abs() < 0.15,
        "ΔA = {mbar_da}, exact {exact}"
    );
}

#[test]
fn mbar_does_not_depend_on_the_order_of_the_pooled_samples() {
    let force_constants = [2.0, 8.0];
    let mut rng = StdRng::seed_from_u64(31);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &force_constants, &[1_500, 500]);
    let reference = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 50_000, 1.0e-13);

    // Reverse the pooled sample order in every row: the estimator only sees the
    // multiset of configurations.
    let reversed: Vec<Vec<f64>> = ensemble
        .energies
        .iter()
        .map(|row| row.iter().rev().copied().collect())
        .collect();
    let permuted = mbar_free_energies(&reversed, &ensemble.n_k, TEMP, 50_000, 1.0e-13);
    assert!(
        (reference[1] - permuted[1]).abs() < 1.0e-9,
        "{} vs {}",
        reference[1],
        permuted[1]
    );
    assert!(ensemble.samples.len() == 2_000);
}

#[test]
fn mbar_shifts_a_state_free_energy_with_its_potential() {
    // Adding c to U₁ everywhere multiplies Z₁ by exp(−βc): f₁ grows by βc exactly.
    let force_constants = [2.0, 8.0];
    let mut rng = StdRng::seed_from_u64(41);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &force_constants, &[1_200, 800]);
    let reference = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 50_000, 1.0e-13);
    let shift_kj_mol = 5.0;
    let mut shifted = ensemble.energies.clone();
    for energy in shifted[1].iter_mut() {
        *energy += shift_kj_mol;
    }
    let result = mbar_free_energies(&shifted, &ensemble.n_k, TEMP, 50_000, 1.0e-13);
    let expected_shift = beta_mol_kjmol(TEMP) * shift_kj_mol;
    assert!(
        ((result[1] - reference[1]) - expected_shift).abs() < 1.0e-6,
        "f₁ moved by {}, expected {expected_shift}",
        result[1] - reference[1]
    );
}

#[test]
fn mbar_of_identical_states_is_zero() {
    let mut rng = StdRng::seed_from_u64(51);
    let ensemble = HarmonicEnsemble::draw(&mut rng, &[3.0, 3.0, 3.0], &[300, 200, 100]);
    let f = mbar_free_energies(&ensemble.energies, &ensemble.n_k, TEMP, 1_000, 1.0e-12);
    assert_eq!(f.len(), 3);
    for (k, value) in f.iter().enumerate() {
        assert!(value.abs() < 1.0e-9, "state {k}: f = {value}");
    }
}

#[test]
fn mbar_single_state_and_empty_pool() {
    let f = mbar_free_energies(&[vec![0.5, 1.5, 2.5]], &[3], TEMP, 100, 1.0e-12);
    assert_eq!(f, vec![0.0]);
    // No samples at all: nothing to estimate, all zeros.
    let empty = mbar_free_energies(&[vec![], vec![]], &[0, 0], TEMP, 100, 1.0e-12);
    assert_eq!(empty, vec![0.0, 0.0]);
}

#[test]
fn mbar_rejects_malformed_input() {
    let rows = vec![vec![0.0, 1.0, 2.0], vec![0.5, 1.5, 2.5]];
    // No states.
    assert!(mbar_free_energies(&[], &[], TEMP, 10, 1.0e-9).is_empty());
    // `n_k` does not have one entry per state.
    assert!(mbar_free_energies(&rows, &[3], TEMP, 10, 1.0e-9).is_empty());
    assert!(mbar_free_energies(&rows, &[1, 1, 1], TEMP, 10, 1.0e-9).is_empty());
    // A row is shorter or longer than the pooled sample count.
    assert!(mbar_free_energies(&rows, &[2, 2], TEMP, 10, 1.0e-9).is_empty());
    assert!(mbar_free_energies(&rows, &[1, 1], TEMP, 10, 1.0e-9).is_empty());
    let ragged = vec![vec![0.0, 1.0, 2.0], vec![0.5, 1.5]];
    assert!(mbar_free_energies(&ragged, &[2, 1], TEMP, 10, 1.0e-9).is_empty());
    // Invalid temperatures.
    assert!(mbar_free_energies(&rows, &[2, 1], 0.0, 10, 1.0e-9).is_empty());
    assert!(mbar_free_energies(&rows, &[2, 1], -300.0, 10, 1.0e-9).is_empty());
    assert!(mbar_free_energies(&rows, &[2, 1], f64::NAN, 10, 1.0e-9).is_empty());
    // A NaN energy cannot be interpreted.
    let nan_rows = vec![vec![0.0, f64::NAN, 2.0], vec![0.5, 1.5, 2.5]];
    assert!(mbar_free_energies(&nan_rows, &[2, 1], TEMP, 10, 1.0e-9).is_empty());
    // The well-formed counterpart is accepted.
    assert_eq!(
        mbar_free_energies(&rows, &[2, 1], TEMP, 10, 1.0e-9).len(),
        2
    );
}

#[test]
fn mbar_rejects_a_configuration_that_no_state_can_generate() {
    let rows = vec![vec![0.0, f64::INFINITY, 2.0], vec![0.5, f64::INFINITY, 2.5]];
    assert!(mbar_free_energies(&rows, &[2, 1], TEMP, 10, 1.0e-9).is_empty());
}

#[test]
fn mbar_delta_a_converts_reduced_free_energies_to_kj_per_mol() {
    let reduced = [0.0, 0.5, -1.25];
    let delta_a = mbar_delta_a(&reduced, TEMP);
    let kt = R_GAS * 1.0e-3 * TEMP;
    assert_eq!(delta_a.len(), 3);
    assert_eq!(delta_a[0], 0.0);
    assert!((delta_a[1] - 0.5 * kt).abs() < 1.0e-12);
    assert!((delta_a[2] + 1.25 * kt).abs() < 1.0e-12);
}
