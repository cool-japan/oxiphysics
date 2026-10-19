// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Monte Carlo routines of the surface-science models.
//!
//! `SurfaceAdsorption::metropolis_step` and `SurfaceDiffusion::monte_carlo_step`
//! drew their "uniform" numbers from `(state >> 33) / u32::MAX`, a 31-bit value
//! divided by 2³² − 1, so every number lay in `[0, 0.5)`:
//!
//! * the Metropolis test `u < exp(−ΔE/kT)` accepted uphill moves with probability
//!   `min(1, 2 exp(−ΔE/kT))`,
//! * the direction test `u < 0.5` was always true, so every hop went to the right
//!   and the surface diffusion was a pure drift.
//!
//! The routines also reseeded their generator on every call, so a loop passing a
//! constant seed drew the same random numbers on every call: with the same hop
//! decisions in every step the diffusing adatoms of a flat surface moved
//! ballistically (mean square displacement `n²` instead of `n`).
//!
//! The generators are seeded, so the statistical checks are deterministic; their
//! tolerances are several standard errors wide.

use super::*;

/// k_B T at 300 K (J).
const KBT: f64 = 4.141_947e-21;

/// Binding energy of the test sites (J, negative = attractive).
const SITE_ENERGY: f64 = -1.0e-20;

/// `n` empty or occupied sites of energy [`SITE_ENERGY`].
fn sites(n: usize, occupied: bool) -> Vec<AdsorptionSite> {
    (0..n)
        .map(|i| {
            let mut site = AdsorptionSite::new([i as f64, 0.0, 0.0], SITE_ENERGY);
            site.occupied = occupied;
            site
        })
        .collect()
}

/// Gas chemical potential for which `(E_site − μ) / k_B T = x`.
fn chemical_potential_for(x: f64) -> f64 {
    SITE_ENERGY - x * KBT
}

/// Equilibrium occupation of an independent site, `1 / (1 + exp(x))`.
fn fermi_dirac_occupation(x: f64) -> f64 {
    1.0 / (1.0 + x.exp())
}

/// A flat-barrier diffusion model whose adatoms hop on every call: the barrier is
/// zero, so `k dt = ν₀ dt = 10 > 1` is clamped to a hop probability of one.
fn always_hopping(n_adatoms: usize) -> SurfaceDiffusion {
    SurfaceDiffusion::new(vec![500.0; n_adatoms], vec![0.0; 1000], 300.0, KBT, 1.0e13)
}

// ---------------------------------------------------------------------------
// SurfaceAdsorption
// ---------------------------------------------------------------------------

#[test]
fn uphill_adsorption_is_accepted_with_probability_exp_minus_delta_e_over_kt() {
    // E_site − μ = +kT: adsorbing costs one kT, so a sweep over empty sites fills
    // the fraction exp(−1) = 0.368.  A generator on [0, 0.5) gave min(1, 2/e) = 0.736.
    for (x, seed) in [(1.0, 3_u64), (2.0, 4)] {
        let mut model = SurfaceAdsorption::new(sites(40_000, false), 300.0, KBT);
        model.metropolis_step(1.0e5, chemical_potential_for(x), seed);
        let expected = (-x).exp();
        let coverage = model.coverage();
        let sigma = (expected * (1.0 - expected) / 40_000.0).sqrt();
        assert!(
            (coverage - expected).abs() < 5.0 * sigma,
            "x = {x}: coverage {coverage}, expected {expected} (σ = {sigma})"
        );
    }
}

#[test]
fn uphill_desorption_is_accepted_with_probability_exp_minus_delta_e_over_kt() {
    // E_site − μ = −kT: leaving an occupied site costs one kT.
    let mut model = SurfaceAdsorption::new(sites(40_000, true), 300.0, KBT);
    model.metropolis_step(1.0e5, chemical_potential_for(-1.0), 11);
    let expected_remaining = 1.0 - (-1.0_f64).exp();
    let coverage = model.coverage();
    let sigma = (expected_remaining * (1.0 - expected_remaining) / 40_000.0).sqrt();
    assert!(
        (coverage - expected_remaining).abs() < 5.0 * sigma,
        "coverage {coverage}, expected {expected_remaining}"
    );
}

#[test]
fn downhill_moves_are_always_accepted() {
    // Adsorbing at E_site − μ = −5 kT lowers the energy: all empty sites fill ...
    let mut filling = SurfaceAdsorption::new(sites(1_000, false), 300.0, KBT);
    filling.metropolis_step(1.0e5, chemical_potential_for(-5.0), 1);
    assert_eq!(filling.coverage(), 1.0);
    // ... and desorbing at E_site − μ = +5 kT empties all occupied sites.
    let mut emptying = SurfaceAdsorption::new(sites(1_000, true), 300.0, KBT);
    emptying.metropolis_step(1.0e5, chemical_potential_for(5.0), 1);
    assert_eq!(emptying.coverage(), 0.0);
}

#[test]
fn metropolis_coverage_relaxes_to_the_fermi_dirac_occupation() {
    // Independent sites: the Metropolis chain samples 1 / (1 + exp((E − μ)/kT)).  A
    // constant seed on every call must not replay the sweep: the sweeps are averaged
    // after a burn-in, with the same seed each time.
    for x in [0.5, -0.5, -3.0] {
        let mut model = SurfaceAdsorption::new(sites(2_000, false), 300.0, KBT);
        let mu = chemical_potential_for(x);
        for _ in 0..50 {
            model.metropolis_step(1.0e5, mu, 5);
        }
        let mut accumulated = 0.0;
        let sweeps = 100;
        for _ in 0..sweeps {
            model.metropolis_step(1.0e5, mu, 5);
            accumulated += model.coverage();
        }
        let mean = accumulated / sweeps as f64;
        let expected = fermi_dirac_occupation(x);
        assert!(
            (mean - expected).abs() < 0.01,
            "(E − μ)/kT = {x}: mean coverage {mean}, Fermi-Dirac {expected}"
        );
    }
}

#[test]
fn a_repeated_seed_continues_the_stream_and_a_new_seed_restarts_it() {
    let mut slot = None;
    let first: f64 = stream_for_seed(&mut slot, 5).random();
    let second: f64 = stream_for_seed(&mut slot, 5).random();
    assert_ne!(first, second, "the same seed must continue the stream");

    // A different seed starts that seed's stream from its beginning ...
    let other: f64 = stream_for_seed(&mut slot, 6).random();
    let mut fresh = None;
    assert_eq!(other, stream_for_seed(&mut fresh, 6).random::<f64>());
    // ... and so does returning to an earlier seed.
    let back: f64 = stream_for_seed(&mut slot, 5).random();
    assert_eq!(back, first);

    // The draws are uniform on [0, 1) with a mean of one half.
    let mut stream = None;
    let draws: Vec<f64> = (0..100_000)
        .map(|_| stream_for_seed(&mut stream, 1).random::<f64>())
        .collect();
    let mean = draws.iter().sum::<f64>() / draws.len() as f64;
    let maximum = draws.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let minimum = draws.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!((mean - 0.5).abs() < 0.01, "mean = {mean}");
    assert!(maximum > 0.99 && minimum < 0.01 && draws.iter().all(|&u| (0.0..1.0).contains(&u)));
}

#[test]
fn metropolis_runs_are_reproducible_for_a_given_sequence_of_seeds() {
    let run = |seeds: &[u64]| -> Vec<bool> {
        let mut model = SurfaceAdsorption::new(sites(500, false), 300.0, KBT);
        for &seed in seeds {
            model.metropolis_step(1.0e5, chemical_potential_for(0.2), seed);
        }
        model.sites.iter().map(|s| s.occupied).collect()
    };
    assert_eq!(run(&[1, 1, 2, 3]), run(&[1, 1, 2, 3]));
    assert_ne!(run(&[1, 1, 2, 3]), run(&[7, 7, 8, 9]));
}

// ---------------------------------------------------------------------------
// SurfaceDiffusion
// ---------------------------------------------------------------------------

#[test]
fn hops_are_equally_likely_in_both_directions() {
    // Every adatom hops on every call: a symmetric random walk with unit steps.  A
    // generator on [0, 0.5) sent all of them to the right (mean displacement +n).
    let n_adatoms = 5_000;
    let steps = 200;
    let mut model = always_hopping(n_adatoms);
    let initial = model.positions.clone();
    for _ in 0..steps {
        let hops = model.monte_carlo_step(1.0e-12, 42);
        assert_eq!(hops, n_adatoms);
    }
    let mean_displacement = model
        .positions
        .iter()
        .zip(initial.iter())
        .map(|(x, x0)| x - x0)
        .sum::<f64>()
        / n_adatoms as f64;
    let msd = model.mean_square_displacement(&initial);
    // Standard error of the mean displacement of a symmetric walk: sqrt(MSD / N).
    let standard_error = (msd / n_adatoms as f64).sqrt();
    assert!(
        mean_displacement.abs() < 4.0 * standard_error,
        "drift {mean_displacement} exceeds four standard errors ({standard_error})"
    );
    // Diffusive growth: MSD = n steps for unit hops (a replayed step would give n²).
    let expected = steps as f64;
    assert!(
        (msd - expected).abs() < 0.1 * expected,
        "MSD = {msd}, expected {expected}"
    );
}

#[test]
fn hop_probability_equals_rate_times_time_step() {
    // Zero barrier: k = ν₀ = 1e13 / s and dt = 0.3e-13 s give a hop probability of
    // 0.3.  A generator on [0, 0.5) hopped with probability 0.6.
    let n_adatoms = 20_000;
    let mut model = always_hopping(n_adatoms);
    let hops = model.monte_carlo_step(0.3e-13, 7);
    let fraction = hops as f64 / n_adatoms as f64;
    let sigma = (0.3_f64 * 0.7 / n_adatoms as f64).sqrt();
    assert!(
        (fraction - 0.3).abs() < 5.0 * sigma,
        "hopping fraction {fraction}, expected 0.3 (σ = {sigma})"
    );
}

#[test]
fn barrier_suppresses_hops_with_the_arrhenius_factor() {
    // One kT of barrier: k dt = exp(−1) with ν₀ dt = 1.
    let n_adatoms = 20_000;
    let mut model = SurfaceDiffusion::new(vec![0.0; n_adatoms], vec![KBT], 300.0, KBT, 1.0e13);
    let hops = model.monte_carlo_step(1.0e-13, 21);
    let expected = (-1.0_f64).exp();
    let fraction = hops as f64 / n_adatoms as f64;
    let sigma = (expected * (1.0 - expected) / n_adatoms as f64).sqrt();
    assert!(
        (fraction - expected).abs() < 5.0 * sigma,
        "hopping fraction {fraction}, expected {expected}"
    );
}

#[test]
fn repeated_calls_with_the_same_seed_do_not_replay_the_step() {
    let n_adatoms = 64;
    let mut model = always_hopping(n_adatoms);
    let start = model.positions.clone();
    model.monte_carlo_step(1.0e-12, 5);
    let after_first = model.positions.clone();
    model.monte_carlo_step(1.0e-12, 5);
    let first_step: Vec<f64> = after_first.iter().zip(&start).map(|(a, b)| a - b).collect();
    let second_step: Vec<f64> = model
        .positions
        .iter()
        .zip(&after_first)
        .map(|(a, b)| a - b)
        .collect();
    assert_ne!(
        first_step, second_step,
        "the second call replayed the hop directions of the first"
    );
    assert!(first_step.iter().all(|d| d.abs() == 1.0));
    assert!(second_step.iter().all(|d| d.abs() == 1.0));
}

#[test]
fn diffusion_runs_are_reproducible_for_a_given_sequence_of_seeds() {
    let run = |seeds: &[u64]| -> Vec<f64> {
        let mut model = always_hopping(100);
        for &seed in seeds {
            model.monte_carlo_step(1.0e-12, seed);
        }
        model.positions
    };
    assert_eq!(run(&[1, 1, 2, 3]), run(&[1, 1, 2, 3]));
    assert_ne!(run(&[1, 1, 2, 3]), run(&[7, 7, 8, 9]));
}
