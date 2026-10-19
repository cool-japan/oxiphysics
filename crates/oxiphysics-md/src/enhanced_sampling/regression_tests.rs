// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the well-tempered metadynamics hill rule, its free
//! energy estimate, and the WHAM-based umbrella-sampling PMF.
//!
//! The metadynamics tests pin the Barducci–Bussi–Parrinello rule
//! `w = w₀ exp(−V / (k_B ΔT))` against the independent
//! [`WellTemperedMetadynamics`] engine and against a converged run on a model
//! free-energy surface.  The umbrella tests sample the biased distributions of
//! known free-energy profiles exactly (they are Gaussian for linear and
//! quadratic profiles) and require WHAM to stitch the windows together.

use super::*;
use crate::enhanced_sampling_md::WellTemperedMetadynamics;
use crate::gaussian_noise::standard_normal;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// Simulation temperature (K) of the metadynamics tests.
const TEMPERATURE: f64 = 300.0;

// ─────────────────────────────────────────────────────────────────────────────
// Well-tempered metadynamics
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn second_well_tempered_hill_is_scaled_by_exp_minus_h_over_kb_delta_t() {
    // After one hill of height h at s the bias at s equals h, so the next hill
    // deposited there has height h exp(−h / (kB ΔT)).  The old rule used
    // kB (T + ΔT) and deposited a visibly taller hill.
    let (height, width, delta_t) = (1.5, 0.3, 450.0);
    let mut meta = Metadynamics::new(height, width, 1, TEMPERATURE, delta_t);
    meta.maybe_deposit(0.7);
    meta.maybe_deposit(0.7);
    assert_eq!(
        meta.hills[0].height, height,
        "empty bias deposits the full w0"
    );
    let expected = height * (-height / (KB * delta_t)).exp();
    assert!(
        (meta.hills[1].height - expected).abs() < 1e-14 * height,
        "second hill {} vs h exp(-h/(kB dT)) = {expected}",
        meta.hills[1].height
    );
    let old_rule = height * (-height / (KB * (TEMPERATURE + delta_t))).exp();
    assert!(
        (meta.hills[1].height - old_rule).abs() > 0.05 * height,
        "the kB (T + dT) rule would give {old_rule}"
    );
}

#[test]
fn plain_metadynamics_keeps_a_constant_hill_height() {
    let mut meta = Metadynamics::new(1.5, 0.3, 1, TEMPERATURE, 0.0);
    for _ in 0..5 {
        meta.maybe_deposit(0.7);
    }
    assert!(meta.hills.iter().all(|hill| hill.height == 1.5));
    assert_eq!(meta.n_hills(), 5);
}

#[test]
fn metadynamics_matches_the_independent_well_tempered_engine() {
    // Same inputs, two implementations: the hill heights, the total bias and the
    // free-energy surface (scale γ/(γ−1) = (T+ΔT)/ΔT) must agree.
    let (height, width, delta_t) = (1.2, 0.25, 600.0);
    let bias_factor = (TEMPERATURE + delta_t) / TEMPERATURE;
    let mut meta = Metadynamics::new(height, width, 1, TEMPERATURE, delta_t);
    let mut reference = WellTemperedMetadynamics::new(height, width, bias_factor, TEMPERATURE, 1);
    for i in 0..40 {
        let s = 0.8 * (f64::from(i) * 0.37).sin();
        meta.maybe_deposit(s);
        let reference_height = reference.deposit_hill(s);
        assert!(
            (meta.hills[i as usize].height - reference_height).abs() < 1e-10 * height,
            "hill {i}: {} vs {reference_height}",
            meta.hills[i as usize].height
        );
    }
    assert_eq!(meta.n_hills(), reference.n_hills());
    // Later hills are strictly smaller than the first: the bias has built up.
    assert!(meta.hills[39].height < 0.5 * height);

    let (grid, reference_fes) = reference.free_energy_surface(-1.0, 1.0, 21);
    let bias_scale = 1.0 / (1.0 - 1.0 / bias_factor);
    assert!((bias_scale - (TEMPERATURE + delta_t) / delta_t).abs() < 1e-12);
    let fes = meta.free_energy_estimate(&grid);
    for (i, (&s, (&ours, &theirs))) in grid.iter().zip(fes.iter().zip(&reference_fes)).enumerate() {
        assert!(
            (meta.compute_bias(s) - reference.bias_at(s)).abs() < 1e-9,
            "bias at grid point {i}"
        );
        assert!(
            (ours - theirs).abs() < 1e-9 * (1.0 + theirs.abs()),
            "FES at grid point {i}: {ours} vs {theirs}"
        );
    }
}

#[test]
fn well_tempered_free_energy_estimate_undoes_the_bias_scaling() {
    let delta_t = 900.0;
    let mut meta = Metadynamics::new(1.0, 0.4, 1, TEMPERATURE, delta_t);
    for i in 0..30 {
        meta.maybe_deposit(0.5 * (f64::from(i) * 0.9).sin());
    }
    let grid: Vec<f64> = (0..=10).map(|i| -0.6 + 0.12 * f64::from(i)).collect();
    let scale = (TEMPERATURE + delta_t) / delta_t;
    let fes = meta.free_energy_estimate(&grid);
    for (&s, &f) in grid.iter().zip(&fes) {
        let expected = -scale * meta.compute_bias(s);
        assert!(
            (f - expected).abs() < 1e-12 * (1.0 + expected.abs()),
            "s = {s}"
        );
        // The bias alone (scale 1) would under-estimate the depth by ΔT/(T+ΔT).
        assert!(f.abs() > 1.2 * meta.compute_bias(s) - 1e-12);
    }
}

#[test]
fn adaptive_metadynamics_uses_the_same_well_tempered_rule() {
    let (height, delta_t) = (2.0, 300.0);
    let mut meta = AdaptiveBiasMetadynamics::new(height, 0.1, 0.4, 0.5, 1, TEMPERATURE, delta_t);
    // First hill: empty bias, no velocity estimate yet, so the widest Gaussian.
    meta.maybe_deposit(0.0);
    assert_eq!(meta.hills[0].height, height);
    assert_eq!(meta.hills[0].width, 0.4);
    // Second hill 0.3 away: the bias of the first one at that point is
    // h exp(−0.3² / (2 · 0.4²)).
    meta.maybe_deposit(0.3);
    let bias = height * (-0.3_f64 * 0.3 / (2.0 * 0.4 * 0.4)).exp();
    let expected = height * (-bias / (KB * delta_t)).exp();
    assert!(
        (meta.hills[1].height - expected).abs() < 1e-14 * height,
        "second hill {} vs {expected}",
        meta.hills[1].height
    );
    // Adaptive width: sigma_min + (sigma_max - sigma_min) exp(-v/vel_scale), v = 0.3.
    let width = 0.1 + 0.3 * (-0.3_f64 / 0.5).exp();
    assert!((meta.hills[1].width - width).abs() < 1e-14);
    let old_rule = height * (-bias / (KB * (TEMPERATURE + delta_t))).exp();
    assert!((meta.hills[1].height - old_rule).abs() > 0.02 * height);
}

/// Asymmetric double well (kJ/mol): minima near −1 and +1, barrier ~ 10 kJ/mol.
fn double_well(s: f64) -> f64 {
    10.0 * (s * s - 1.0).powi(2) + 2.0 * s
}

/// Draw a grid index with probability proportional to `exp(log_weights[i])`.
fn draw_index(rng: &mut StdRng, log_weights: &[f64]) -> usize {
    let max = log_weights
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let weights: Vec<f64> = log_weights.iter().map(|&lw| (lw - max).exp()).collect();
    let total: f64 = weights.iter().sum();
    let target = rng.random::<f64>() * total;
    let mut cumulative = 0.0;
    for (index, weight) in weights.iter().enumerate() {
        cumulative += weight;
        if cumulative >= target {
            return index;
        }
    }
    weights.len() - 1
}

#[test]
fn well_tempered_metadynamics_converges_to_the_free_energy_surface() {
    // Quasi-static metadynamics on a known free-energy profile: at every
    // deposition the collective variable is drawn from the exact biased
    // distribution p(s) ∝ exp(−(F(s) + V(s)) / kT).  In the long-time limit the
    // bias becomes V = −ΔT/(T+ΔT) F + c, so −(T+ΔT)/ΔT V recovers F up to a
    // constant.  With ΔT = T the bias only compensates half of the surface: the
    // unscaled estimate −V would have the barrier at half of its true height.
    let k_t = KB * TEMPERATURE;
    let delta_t = TEMPERATURE;
    let mut meta = Metadynamics::new(1.0, 0.2, 1, TEMPERATURE, delta_t);
    let grid: Vec<f64> = (0..=220).map(|i| -2.2 + 0.02 * f64::from(i)).collect();
    let mut bias = vec![0.0_f64; grid.len()];
    let mut rng = StdRng::seed_from_u64(3215);
    for _ in 0..2500 {
        let log_weights: Vec<f64> = grid
            .iter()
            .zip(&bias)
            .map(|(&s, &v)| -(double_well(s) + v) / k_t)
            .collect();
        let s = grid[draw_index(&mut rng, &log_weights)];
        meta.maybe_deposit(s);
        let hill = &meta.hills[meta.hills.len() - 1];
        for (&g, v) in grid.iter().zip(bias.iter_mut()) {
            *v += hill.evaluate(g);
        }
    }

    let estimate = meta.free_energy_estimate(&grid);
    let (minimum, _) = grid
        .iter()
        .enumerate()
        .filter(|&(_, &s)| s.abs() < 1.5)
        .map(|(i, &s)| (i, double_well(s)))
        .fold((0, f64::INFINITY), |best, item| {
            if item.1 < best.1 { item } else { best }
        });
    let mut worst = 0.0_f64;
    // Compare where the surface is actually sampled: beyond |s| = 1.2 it climbs
    // steeply (about 20 kJ/mol above the minimum at |s| = 1.5) and is hardly visited.
    for (i, &s) in grid.iter().enumerate().filter(|&(_, &s)| s.abs() <= 1.2) {
        let truth = double_well(s) - double_well(grid[minimum]);
        let found = estimate[i] - estimate[minimum];
        worst = worst.max((found - truth).abs());
        assert!(
            (found - truth).abs() < 1.5,
            "s = {s}: estimated F = {found}, exact {truth}"
        );
    }
    assert!(worst < 1.5, "largest deviation {worst} kJ/mol");

    // The barrier at s ~ 0 is recovered to 15 %; −V alone gives only about half of it.
    let barrier_index = grid
        .iter()
        .position(|&s| s.abs() < 1e-9)
        .expect("grid contains s = 0");
    let true_barrier = double_well(0.0) - double_well(grid[minimum]);
    let estimated_barrier = estimate[barrier_index] - estimate[minimum];
    assert!(
        (estimated_barrier - true_barrier).abs() < 0.15 * true_barrier,
        "barrier {estimated_barrier} vs exact {true_barrier}"
    );
    let unscaled_barrier = -(meta.compute_bias(0.0) - meta.compute_bias(grid[minimum]));
    assert!(
        (unscaled_barrier - true_barrier).abs() > 0.3 * true_barrier,
        "unscaled barrier {unscaled_barrier} would be accepted"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Umbrella sampling / WHAM
// ─────────────────────────────────────────────────────────────────────────────

/// Thermal energy (kJ/mol) of the umbrella tests (T ≈ 298 K).
const K_T: f64 = 2.479;

/// Umbrella windows on the profile `U(ξ) = ½ stiffness ξ²  +  slope ξ`.
///
/// The biased density `exp(−(U + V_i)/kT)` with the window potential
/// `V_i = k/2 (ξ − c_i)²` is exactly Gaussian with precision
/// `(k + stiffness)/kT` and mean `(k c_i − slope)/(k + stiffness)`.
fn sample_windows(
    rng: &mut StdRng,
    centres: &[f64],
    force_constant: f64,
    stiffness: f64,
    slope: f64,
    samples_per_window: usize,
) -> UmbrellaSampling {
    let mut umbrella = UmbrellaSampling::from_centers(centres, force_constant);
    let precision = (force_constant + stiffness) / K_T;
    for window in &mut umbrella.windows {
        let mean = (force_constant * window.center - slope) / (force_constant + stiffness);
        for _ in 0..samples_per_window {
            window.add_sample(mean + standard_normal(rng) / precision.sqrt());
        }
    }
    umbrella
}

/// Least-squares slope and intercept of `y(x)` over the points with
/// `x ∈ [lo, hi]` and finite `y`.
fn fit_line(x: &[f64], y: &[f64], lo: f64, hi: f64) -> (f64, f64) {
    let points: Vec<(f64, f64)> = x
        .iter()
        .zip(y)
        .filter(|&(&xi, &yi)| (lo..=hi).contains(&xi) && yi.is_finite())
        .map(|(&xi, &yi)| (xi, yi))
        .collect();
    let n = points.len() as f64;
    let mean_x = points.iter().map(|p| p.0).sum::<f64>() / n;
    let mean_y = points.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = points.iter().map(|p| (p.0 - mean_x).powi(2)).sum();
    let sxy: f64 = points.iter().map(|p| (p.0 - mean_x) * (p.1 - mean_y)).sum();
    let slope = sxy / sxx;
    (slope, mean_y - slope * mean_x)
}

#[test]
fn wham_recovers_a_linear_free_energy_slope_across_ten_windows() {
    // U(ξ) = 2 kT ξ sampled with ten harmonic windows spaced 0.3 apart.  Shifting
    // every window's PMF by its own minimum (the old code) discards the relative
    // window offsets and flattens the slope towards zero.
    let slope = 2.0 * K_T;
    let centres: Vec<f64> = (0..10).map(|i| 0.3 * f64::from(i)).collect();
    for seed in [3216_u64, 99] {
        let mut rng = StdRng::seed_from_u64(seed);
        let umbrella = sample_windows(&mut rng, &centres, 40.0 * K_T, 0.0, slope, 16_000);
        let (bins, pmf) = umbrella.compute_pmf(60, K_T);
        assert_eq!(bins.len(), 60);
        // Interior bins only: the outermost ones are sampled by a single window's tail.
        let (fitted, _) = fit_line(&bins, &pmf, 0.1, 2.5);
        assert!(
            (fitted / slope - 1.0).abs() < 0.05,
            "seed {seed}: fitted slope {fitted} kJ/mol per unit, exact {slope}"
        );
        // The profile is straight: no interior bin leaves the fitted line by more than 0.3 kT.
        let (fitted, intercept) = fit_line(&bins, &pmf, 0.1, 2.5);
        for (&x, &f) in bins
            .iter()
            .zip(&pmf)
            .filter(|&(&x, _)| (0.1..=2.5).contains(&x))
        {
            assert!(
                (f - (intercept + fitted * x)).abs() < 0.3 * K_T,
                "seed {seed}: PMF({x}) = {f} deviates from the line"
            );
        }
        let finite_minimum = pmf
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::INFINITY, f64::min);
        assert_eq!(finite_minimum, 0.0, "the PMF is shifted to a zero minimum");
    }
}

#[test]
fn wham_recovers_the_curvature_of_a_harmonic_profile() {
    // U(ξ) = ½ a ξ² with a = 30 kT per unit² (the windows' own force constant is 40 kT):
    // −kT ln P(ξ) must reproduce a quadratic with that curvature.
    let stiffness = 30.0 * K_T;
    let centres: Vec<f64> = (-6..=6).map(|i| 0.1 * f64::from(i)).collect();
    let mut rng = StdRng::seed_from_u64(17);
    let umbrella = sample_windows(&mut rng, &centres, 40.0 * K_T, stiffness, 0.0, 3000);
    let (bins, pmf) = umbrella.compute_pmf(50, K_T);
    // PMF = c + (a/2) ξ² is linear in ξ²: regress over the well-sampled interior.
    let squares: Vec<f64> = bins.iter().map(|x| x * x).collect();
    let n_interior = squares
        .iter()
        .zip(&pmf)
        .filter(|&(&s, &f)| s <= 0.16 && f.is_finite())
        .count();
    assert!(n_interior > 15, "only {n_interior} interior bins");
    let (half_curvature, _) = fit_line(&squares, &pmf, 0.0, 0.16);
    let curvature = 2.0 * half_curvature;
    assert!(
        (curvature / stiffness - 1.0).abs() < 0.06,
        "curvature {curvature} vs exact {stiffness}"
    );
}

#[test]
fn wham_is_independent_of_window_order_and_handles_unequal_sample_counts() {
    let slope = 2.0 * K_T;
    let centres: Vec<f64> = (0..8).map(|i| 0.3 * f64::from(i)).collect();
    let mut rng = StdRng::seed_from_u64(5);
    let mut umbrella = sample_windows(&mut rng, &centres, 40.0 * K_T, 0.0, slope, 12_000);
    // Thin out some windows: WHAM weights windows by their sample number.
    umbrella.windows[1].samples.truncate(2000);
    umbrella.windows[4].samples.truncate(3600);
    umbrella.windows[6].samples.truncate(1200);
    let (bins, pmf) = umbrella.compute_pmf(40, K_T);
    let (fitted, _) = fit_line(&bins, &pmf, 0.1, 1.9);
    assert!(
        (fitted / slope - 1.0).abs() < 0.08,
        "fitted slope {fitted}, exact {slope}"
    );

    let mut reversed = umbrella.clone();
    reversed.windows.reverse();
    let (bins_rev, pmf_rev) = reversed.compute_pmf(40, K_T);
    assert_eq!(bins, bins_rev);
    for (a, b) in pmf.iter().zip(&pmf_rev) {
        assert!(
            (a.is_infinite() && b.is_infinite()) || (a - b).abs() < 1e-6,
            "{a} vs {b}"
        );
    }
}

#[test]
fn flat_histogram_under_a_stiff_window_unbiases_to_minus_the_window_potential() {
    // One window (k = 200) and one sample in each of ten bins.  The unbiased density
    // is n(b) exp(+V(b)/kT), so PMF(b) = −V(b) + const: a histogram that is flat
    // under a stiff bias belongs to a profile that FALLS away from the window
    // centre, whose bins therefore hold the PMF maximum.
    let mut umbrella = UmbrellaSampling::new();
    umbrella.add_window(1.0, 200.0);
    for j in 0..10 {
        umbrella.windows[0].add_sample(0.55 + 0.1 * f64::from(j));
    }
    let (bins, pmf) = umbrella.compute_pmf(10, K_T);
    assert_eq!(bins.len(), 10);
    let bias: Vec<f64> = bins
        .iter()
        .map(|&b| umbrella.windows[0].bias_potential(b))
        .collect();
    let bias_max = bias.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    for (b, (&p, &v)) in pmf.iter().zip(&bias).enumerate() {
        assert!(
            (p - (bias_max - v)).abs() < 1e-9,
            "bin {b}: PMF {p}, expected {}",
            bias_max - v
        );
    }
    // The outermost bins (largest bias, equal by symmetry) have the lowest PMF, the
    // bins next to the window centre the highest.
    assert!(pmf[0].abs() < 1e-9 && pmf[9].abs() < 1e-9, "{pmf:?}");
    assert!(pmf[4] > pmf[2] && pmf[2] > pmf[0]);
    assert!(pmf[5] > pmf[7] && pmf[7] > pmf[9]);
}

#[test]
fn wham_pmf_is_flat_for_a_flat_profile() {
    let centres: Vec<f64> = (0..8).map(|i| 0.3 * f64::from(i)).collect();
    let mut rng = StdRng::seed_from_u64(31);
    let umbrella = sample_windows(&mut rng, &centres, 40.0 * K_T, 0.0, 0.0, 4000);
    let (bins, pmf) = umbrella.compute_pmf(30, K_T);
    // Over the interior the PMF of U = 0 stays within 0.3 kT of its mean.
    let interior: Vec<f64> = bins
        .iter()
        .zip(&pmf)
        .filter(|&(&x, &f)| (0.1..=1.9).contains(&x) && f.is_finite())
        .map(|(_, &f)| f)
        .collect();
    let mean = interior.iter().sum::<f64>() / interior.len() as f64;
    for f in interior {
        assert!(
            (f - mean).abs() < 0.3 * K_T,
            "flat PMF value {f}, mean {mean}"
        );
    }
}

#[test]
fn compute_pmf_handles_degenerate_input() {
    let mut umbrella = UmbrellaSampling::from_centers(&[0.0, 1.0], 100.0);
    // No samples at all.
    assert_eq!(umbrella.compute_pmf(10, K_T), (vec![], vec![]));
    umbrella.windows[0].add_sample(0.2);
    umbrella.windows[0].add_sample(f64::NAN);
    umbrella.windows[0].add_sample(f64::INFINITY);
    // One finite sample: a single degenerate bin at that value.
    assert_eq!(umbrella.compute_pmf(10, K_T), (vec![0.2], vec![0.0]));
    umbrella.windows[1].add_sample(0.9);
    // Invalid binning and temperature.
    assert_eq!(umbrella.compute_pmf(0, K_T), (vec![], vec![]));
    assert_eq!(umbrella.compute_pmf(10, 0.0), (vec![], vec![]));
    assert_eq!(umbrella.compute_pmf(10, -1.0), (vec![], vec![]));
    assert_eq!(umbrella.compute_pmf(10, f64::NAN), (vec![], vec![]));
    // A valid call ignores the non-finite samples and a sample-free window.
    let mut sparse = UmbrellaSampling::from_centers(&[0.0, 0.5, 5.0], 100.0);
    sparse.windows[0].add_sample(0.0);
    sparse.windows[0].add_sample(0.1);
    sparse.windows[1].add_sample(0.4);
    sparse.windows[1].add_sample(f64::NAN);
    let (bins, pmf) = sparse.compute_pmf(4, K_T);
    assert_eq!(bins.len(), 4);
    assert!(pmf.contains(&0.0), "{pmf:?}");
    assert!(pmf.iter().all(|f| *f >= 0.0));
}
