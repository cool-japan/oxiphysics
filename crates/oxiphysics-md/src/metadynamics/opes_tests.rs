// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the reweighted OPES estimator ([`OpesState`]).
//!
//! The samples of an OPES run come from the *biased* ensemble, so each kernel
//! must carry the weight `exp(V(s_k) / k_B T)` of the bias that produced it.
//! Kernels of unit height estimate the biased distribution, which flattens the
//! bias instead of letting it reach the well-tempered target.

use super::functions::KB;
use super::types::{GaussianHill, OpesState};
use crate::error::MdError;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

const TEMPERATURE: f64 = 300.0;

fn kt() -> f64 {
    KB * TEMPERATURE
}

fn assert_close(actual: f64, expected: f64, rel: f64, what: &str) {
    let scale = expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        (actual - expected).abs() <= rel * scale,
        "{what}: {actual} vs {expected}"
    );
}

/// `Z` recomputed as the plain average of the density estimate over the centres.
fn brute_force_normalisation(opes: &OpesState) -> f64 {
    let sum: f64 = opes
        .kernels
        .iter()
        .map(|k| opes.kernel_density(&k.center))
        .sum();
    sum / opes.kernels.len() as f64
}

#[test]
fn first_kernel_has_unit_weight_and_unit_normalisation() {
    let mut opes = OpesState::new(10.0, TEMPERATURE, 1, vec![0.5]);
    assert_eq!(opes.reweighting_factor(&[0.3]), 1.0);
    opes.update(&[0.3]);
    assert_eq!(opes.kernels.len(), 1);
    assert_eq!(opes.kernels[0].height, 1.0);
    assert_close(opes.z, 1.0, 1e-14, "Z after one kernel");
}

#[test]
fn second_kernel_at_the_same_cv_carries_the_reweighting_factor() {
    let gamma = 10.0;
    let mut opes = OpesState::new(gamma, TEMPERATURE, 1, vec![0.5]);
    opes.update(&[0.3]);
    let bias_before = opes.bias_potential(&[0.3]);
    assert!(bias_before > 0.0, "bias at the only kernel: {bias_before}");
    opes.update(&[0.3]);
    let weight = opes.kernels[1].height;
    // w = exp(V / kT) = (1 + epsilon)^(1 - 1/gamma) since P_hat / Z = 1 there.
    let expected = (1.0 + opes.epsilon()).powf(1.0 - 1.0 / gamma);
    assert_close(weight, expected, 1e-12, "closed-form weight");
    assert_close(weight, (bias_before / kt()).exp(), 1e-12, "exp(V/kT)");
    assert!(weight > 1.0, "weight {weight}");
    // Both centres coincide, so Z is the density 1 + w.
    assert_close(
        opes.z,
        1.0 + weight,
        1e-12,
        "Z after two coincident kernels",
    );
}

#[test]
fn kernel_far_from_the_sampled_region_gets_the_barrier_weight() {
    let gamma = 8.0;
    let mut opes = OpesState::new(gamma, TEMPERATURE, 1, vec![0.5]);
    opes.update(&[0.0]);
    // 100 bandwidths away the Gaussian tail underflows, so P_hat / Z = 0 and the
    // ratio is epsilon: V = -barrier and w = exp(-barrier / kT).
    let far = [50.0];
    let barrier = opes.barrier;
    assert_close(barrier, gamma * kt(), 1e-14, "default barrier");
    assert_close(opes.bias_potential(&far), -barrier, 1e-9, "bias floor");
    opes.update(&far);
    assert_close(
        opes.kernels[1].height,
        (-barrier / kt()).exp(),
        1e-9,
        "weight at the floor",
    );
}

#[test]
fn explicit_barrier_sets_epsilon_and_the_bias_floor() {
    let gamma = 6.0;
    let barrier = 12.0;
    let mut opes = OpesState::new(gamma, TEMPERATURE, 1, vec![0.4])
        .with_barrier(barrier)
        .expect("a positive barrier is valid");
    assert_eq!(opes.barrier, barrier);
    let expected_epsilon = (-barrier / ((1.0 - 1.0 / gamma) * kt())).exp();
    assert_close(opes.epsilon(), expected_epsilon, 1e-14, "epsilon");
    opes.update(&[1.0]);
    assert_close(
        opes.bias_potential(&[100.0]),
        -barrier,
        1e-9,
        "floor = -barrier",
    );
    // A deeper barrier lowers the floor of the bias.
    let deeper = OpesState::new(gamma, TEMPERATURE, 1, vec![0.4])
        .with_barrier(2.0 * barrier)
        .expect("a positive barrier is valid");
    assert!(deeper.epsilon() < opes.epsilon());
}

#[test]
fn invalid_barrier_is_rejected_and_infinite_barrier_removes_epsilon() {
    let state = || OpesState::new(10.0, TEMPERATURE, 1, vec![0.5]);
    assert!(matches!(
        state().with_barrier(-1.0),
        Err(MdError::General(_))
    ));
    assert!(matches!(
        state().with_barrier(f64::NAN),
        Err(MdError::General(_))
    ));
    let mut unregularised = state()
        .with_barrier(f64::INFINITY)
        .expect("an infinite barrier is allowed");
    assert_eq!(unregularised.epsilon(), 0.0);
    unregularised.update(&[0.0]);
    let far = unregularised.bias_potential(&[1.0e3]);
    assert!(far.is_finite(), "bias far from all kernels: {far}");
    assert!(far < 0.0);
}

#[test]
fn normalisation_is_the_mean_density_over_the_kernel_centres() {
    let mut opes = OpesState::new(7.0, TEMPERATURE, 2, vec![0.3, 0.6]);
    let mut rng = StdRng::seed_from_u64(11);
    for _ in 0..80 {
        let cv = [rng.random::<f64>() * 2.0 - 1.0, rng.random::<f64>() * 3.0];
        opes.update(&cv);
    }
    assert_eq!(opes.kernels.len(), 40, "one kernel every second step");
    assert_close(
        opes.z,
        brute_force_normalisation(&opes),
        1e-12,
        "incremental Z",
    );
    let incremental = opes.z;
    opes.recompute_normalisation();
    assert_close(opes.z, incremental, 1e-12, "recomputed Z");
    // The incremental update also holds when the kernels have different widths.
    let mut mixed = OpesState::new(7.0, TEMPERATURE, 1, vec![0.2]);
    mixed.update(&[0.0]);
    mixed.bandwidth = vec![0.9];
    mixed.update(&[0.4]);
    mixed.bandwidth = vec![0.05];
    mixed.update(&[0.1]);
    assert_close(
        mixed.z,
        brute_force_normalisation(&mixed),
        1e-12,
        "mixed widths",
    );
}

#[test]
fn no_bias_and_unit_weights_for_gamma_at_most_one() {
    let mut opes = OpesState::new(1.0, TEMPERATURE, 1, vec![0.5]);
    for s in [0.0, 0.2, 0.2, 0.4] {
        opes.update(&[s]);
    }
    assert_eq!(opes.bias_potential(&[0.2]), 0.0);
    assert!(opes.kernels.iter().all(|k| k.height == 1.0));
    assert_close(
        opes.z,
        brute_force_normalisation(&opes),
        1e-12,
        "Z without bias",
    );
}

#[test]
fn compress_merges_kernels_and_refreshes_the_normalisation() {
    let mut opes = OpesState::new(10.0, TEMPERATURE, 1, vec![0.5]);
    for s in [0.0, 0.01, 0.02, 3.0, 3.01, 6.0] {
        opes.update(&[s]);
    }
    let total_weight: f64 = opes.kernels.iter().map(|k| k.height).sum();
    opes.compress(0.5);
    assert_eq!(opes.kernels.len(), 3, "three clusters remain");
    let merged_weight: f64 = opes.kernels.iter().map(|k| k.height).sum();
    assert_close(merged_weight, total_weight, 1e-12, "weights are conserved");
    assert_close(
        opes.z,
        brute_force_normalisation(&opes),
        1e-12,
        "Z after merge",
    );
}

#[test]
fn free_energy_estimate_handles_empty_and_uncovered_points() {
    let mut opes = OpesState::new(10.0, TEMPERATURE, 1, vec![0.5]);
    let grid: Vec<Vec<f64>> = vec![vec![0.0], vec![1.0], vec![1.0e3]];
    assert_eq!(opes.free_energy_estimate(&grid), vec![0.0; 3]);
    opes.update(&[0.0]);
    let free_energy = opes.free_energy_estimate(&grid);
    assert_eq!(free_energy[0], 0.0, "highest density has F = 0");
    // A unit Gaussian of width 0.5 has density exp(-2) one unit away from its centre.
    assert_close(free_energy[1], 2.0 * kt(), 1e-12, "Gaussian tail");
    assert!(free_energy[2].is_infinite() && free_energy[2] > 0.0);
}

/// Double well `U(s) = h ((s/L)^2 - 1)^2` with minima at `±L`.
fn double_well(s: f64, barrier: f64, half_width: f64) -> f64 {
    let x = s / half_width;
    barrier * (x * x - 1.0).powi(2)
}

/// Metropolis Monte Carlo in `U(s) + V(s)` with OPES depositing every `stride`
/// moves.  Returns the final state.
fn run_opes_double_well(gamma: f64, n_steps: usize, stride: usize, seed: u64) -> OpesState {
    let (barrier, half_width, bandwidth, max_move) = (5.0 * kt(), 3.0, 0.15, 0.5);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut opes = OpesState::new(gamma, TEMPERATURE, stride, vec![bandwidth]);
    let mut s = -half_width;
    for _ in 0..n_steps {
        opes.update(&[s]);
        let proposal = s + max_move * (2.0 * rng.random::<f64>() - 1.0);
        let delta = (double_well(proposal, barrier, half_width) + opes.bias_potential(&[proposal]))
            - (double_well(s, barrier, half_width) + opes.bias_potential(&[s]));
        if delta <= 0.0 || rng.random::<f64>() < (-delta / kt()).exp() {
            s = proposal;
        }
    }
    opes
}

/// `(max, rms)` deviation (units of kT) between `F_est(s) - offset` and `U(s)`
/// over the region `U - U_min < window`, the offset being the mean difference.
fn free_energy_error(estimate: &[f64], grid: &[f64], window: f64) -> (f64, f64) {
    let (barrier, half_width) = (5.0 * kt(), 3.0);
    let reference: Vec<f64> = grid
        .iter()
        .map(|&s| double_well(s, barrier, half_width))
        .collect();
    let inside: Vec<usize> = (0..grid.len())
        .filter(|&i| reference[i] < window * kt() && estimate[i].is_finite())
        .collect();
    let offset = inside
        .iter()
        .map(|&i| estimate[i] - reference[i])
        .sum::<f64>()
        / inside.len() as f64;
    let deviations: Vec<f64> = inside
        .iter()
        .map(|&i| (estimate[i] - reference[i] - offset).abs() / kt())
        .collect();
    let max = deviations.iter().copied().fold(0.0_f64, f64::max);
    let rms = (deviations.iter().map(|d| d * d).sum::<f64>() / deviations.len() as f64).sqrt();
    (max, rms)
}

fn evaluation_grid() -> (Vec<f64>, Vec<Vec<f64>>) {
    let points: Vec<f64> = (0..=1200)
        .map(|i| -6.0 + 12.0 * i as f64 / 1200.0)
        .collect();
    let as_cv = points.iter().map(|&s| vec![s]).collect();
    (points, as_cv)
}

#[test]
fn reweighted_kernels_recover_the_unbiased_double_well_free_energy() {
    let (points, cv_grid) = evaluation_grid();
    for (gamma, seed) in [(5.0, 1), (10.0, 2)] {
        let opes = run_opes_double_well(gamma, 40_000, 50, seed);
        assert_eq!(opes.kernels.len(), 800);
        let (max, rms) = free_energy_error(&opes.free_energy_estimate(&cv_grid), &points, 4.0);
        // 800 kernels estimate the density only to a statistical accuracy that
        // is worst at the 4 kT edge of the window: over independent seeds the
        // maximum error is 0.3-0.65 kT and the rms error 0.1-0.26 kT, so the
        // bounds below leave a margin without admitting a biased estimate (the
        // unit-weight estimate of the test further down is off by 2.5 kT).
        assert!(
            max < 1.0 && rms < 0.35,
            "gamma {gamma}: max error {max} kT, rms {rms} kT"
        );
    }
}

#[test]
fn run_visits_both_wells_and_flattens_the_biased_distribution() {
    let opes = run_opes_double_well(10.0, 40_000, 50, 3);
    let left = opes.kernels.iter().filter(|k| k.center[0] < 0.0).count();
    let right = opes.kernels.len() - left;
    assert!(
        left > 100 && right > 100,
        "the bias must let the walker cross the barrier: {left} / {right}"
    );
}

#[test]
fn unit_weights_would_estimate_the_biased_distribution_not_the_target() {
    let (points, cv_grid) = evaluation_grid();
    let opes = run_opes_double_well(10.0, 40_000, 50, 2);
    let weighted = opes.free_energy_estimate(&cv_grid);
    let (max_weighted, _) = free_energy_error(&weighted, &points, 4.0);
    // The same sample positions with unit kernel heights.
    let mut unweighted = opes.clone();
    for kernel in &mut unweighted.kernels {
        *kernel = GaussianHill::new(kernel.center.clone(), 1.0, kernel.widths.clone());
    }
    let (max_unweighted, rms_unweighted) =
        free_energy_error(&unweighted.free_energy_estimate(&cv_grid), &points, 4.0);
    assert!(
        max_unweighted > 2.0 * max_weighted && rms_unweighted > 0.8,
        "unit weights: max {max_unweighted} kT, rms {rms_unweighted} kT; weighted max {max_weighted} kT"
    );
}
