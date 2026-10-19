// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the entropic root search and for the viscosity
//! dependence of the entropic, KBC and EOR collisions.

use super::*;
use std::f64::consts::PI;

/// Density and velocity of a population vector.
fn moments(f: &[f64; 9]) -> (f64, [f64; 2]) {
    let rho: f64 = f.iter().sum();
    let mut j = [0.0; 2];
    for i in 0..9 {
        j[0] += f[i] * C[i].0 as f64;
        j[1] += f[i] * C[i].1 as f64;
    }
    (rho, [j[0] / rho, j[1] / rho])
}

/// Entropic equilibrium with the density and momentum of `f`.
fn local_equilibrium(f: &[f64; 9]) -> [f64; 9] {
    let (rho, u) = moments(f);
    EntropicD2Q9::equilibrium(rho, u)
}

/// Equilibrium at rho = 1, u = (0.05, 0.02) plus a mass- and momentum-free
/// stress perturbation of relative size `scale` (`scale = 1` is a 90 % shear
/// and 60 % normal-stress perturbation, i.e. far from equilibrium).
fn perturbed_state(scale: f64) -> [f64; 9] {
    let mut f = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        f[i] += W[i] * scale * (0.9 * cx * cy + 0.6 * (cx * cx - cy * cy));
    }
    f
}

/// Entropy change H(f_post) - H(f_pre).
fn entropy_gain(f_pre: &[f64; 9], f_post: &[f64; 9]) -> f64 {
    h_function(f_post) - h_function(f_pre)
}

// ---------------------------------------------------------------------------
// Equilibrium
// ---------------------------------------------------------------------------

#[test]
fn entropic_equilibrium_conserves_mass_and_momentum() {
    for &(rho, u) in &[
        (1.0, [0.0, 0.0]),
        (1.2, [0.05, -0.03]),
        (0.8, [0.1, 0.06]),
        (1.0, [-0.2, 0.15]),
    ] {
        let f = d2q9_equilibrium(rho, u);
        let (rho_f, u_f) = moments(&f);
        assert!((rho_f - rho).abs() < 1e-14, "rho {rho_f} != {rho}");
        assert!((u_f[0] - u[0]).abs() < 1e-14, "ux {} != {}", u_f[0], u[0]);
        assert!((u_f[1] - u[1]).abs() < 1e-14, "uy {} != {}", u_f[1], u[1]);
        assert!(f.iter().all(|&v| v > 0.0), "equilibrium must be positive");
    }
}

#[test]
fn entropic_equilibrium_matches_polynomial_to_third_order() {
    let rest = d2q9_equilibrium(1.0, [0.0, 0.0]);
    for i in 0..9 {
        assert!((rest[i] - W[i]).abs() < 1e-15, "rest state must be w_i");
    }
    let mut previous = f64::INFINITY;
    for &mach in &[0.1, 0.05, 0.025] {
        let u = [mach, 0.6 * mach];
        let exact = d2q9_equilibrium(1.0, u);
        let poly = d2q9_polynomial_equilibrium(1.0, u);
        let diff = exact
            .iter()
            .zip(poly.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            diff > 0.0 && diff < 1e-3 * (mach / 0.1).powi(3) + 1e-12,
            "diff {diff}"
        );
        assert!(
            diff < previous,
            "difference must shrink with the Mach number"
        );
        previous = diff;
    }
}

#[test]
fn entropic_equilibrium_minimises_the_h_function() {
    let (rho, u) = (1.0, [0.1, 0.06]);
    let exact = d2q9_equilibrium(rho, u);
    let poly = d2q9_polynomial_equilibrium(rho, u);
    assert!(
        h_function(&exact) < h_function(&poly),
        "the entropic equilibrium must have lower H than the polynomial one"
    );
    // Any mass- and momentum-preserving perturbation raises H.
    for scale in [1e-3, -1e-3, 1e-2, -1e-2] {
        let mut f = exact;
        for i in 0..9 {
            let cx = C[i].0 as f64;
            let cy = C[i].1 as f64;
            f[i] += W[i] * scale * (0.9 * cx * cy + 0.6 * (cx * cx - cy * cy));
        }
        assert!(h_function(&f) > h_function(&exact), "scale {scale}");
    }
}

// ---------------------------------------------------------------------------
// Entropic root search
// ---------------------------------------------------------------------------

#[test]
fn find_entropic_alpha_is_not_identically_two() {
    // 45 % shear / 30 % normal-stress perturbation of the distribution.
    let f = perturbed_state(0.5);
    let f_eq = local_equilibrium(&f);
    let d = mirror_direction(&f, &f_eq);
    let alpha = find_entropic_alpha(&f, &f_eq, 50);
    assert!(
        (alpha - 2.0).abs() > 1e-3,
        "strongly perturbed state must not get alpha = 2, got {alpha}"
    );
    assert!(alpha > 1.0, "root lies above 1, got {alpha}");
    // alpha is the iso-entropy point: H(f + alpha d) == H(f) and the H-theorem holds.
    let mut f_mirror = f;
    for i in 0..9 {
        f_mirror[i] += alpha * d[i];
    }
    let gain = entropy_gain(&f, &f_mirror);
    assert!(gain <= 1e-12, "H must not increase, gain {gain}");
    assert!(
        gain.abs() < 1e-9,
        "alpha must be the iso-entropy root, gain {gain}"
    );
}

#[test]
fn find_entropic_alpha_stays_inside_the_positive_orthant() {
    // At this perturbation f_7 would hit zero (alpha = 1.8967) before the
    // entropy returns to H(f): there is no iso-entropy root and alpha is the
    // positivity limit, backed off so that no population vanishes.
    let f = perturbed_state(1.0);
    let f_eq = local_equilibrium(&f);
    let d = mirror_direction(&f, &f_eq);
    let positivity_limit = (0..9)
        .filter(|&i| d[i] < 0.0)
        .map(|i| f[i] / -d[i])
        .fold(f64::INFINITY, f64::min);
    for alpha in [
        find_entropic_alpha(&f, &f_eq, 60),
        find_entropic_alpha_newton(&f, &f_eq, 40),
    ] {
        assert!(
            alpha < positivity_limit,
            "alpha {alpha} >= {positivity_limit}"
        );
        assert!(
            alpha > positivity_limit * (1.0 - 1e-5),
            "alpha {alpha} should sit just below the positivity limit {positivity_limit}"
        );
        let mut f_mirror = f;
        for i in 0..9 {
            f_mirror[i] += alpha * d[i];
        }
        assert!(
            f_mirror.iter().all(|&v| v > 0.0),
            "mirror state not positive"
        );
        assert!(entropy_gain(&f, &f_mirror) < 0.0, "H must decrease");
    }
}

#[test]
fn find_entropic_alpha_newton_agrees_with_bisection() {
    for &scale in &[0.1, 0.3, 0.6, 0.9, 1.0] {
        let f = perturbed_state(scale);
        let f_eq = local_equilibrium(&f);
        let bisect = find_entropic_alpha(&f, &f_eq, 60);
        let newton = find_entropic_alpha_newton(&f, &f_eq, 40);
        assert!(
            (bisect - newton).abs() < 1e-9,
            "scale {scale}: bisection {bisect} vs Newton {newton}"
        );
        assert!((bisect - 2.0).abs() > 1e-5, "scale {scale}: alpha {bisect}");
    }
}

#[test]
fn solve_alpha_newton_raphson_finds_the_nontrivial_root() {
    let f = perturbed_state(0.5);
    let f_eq = local_equilibrium(&f);
    let mirror = mirror_state(&f, &f_eq);
    let alpha = solve_alpha_newton_raphson(&f, &f_eq, &mirror, 1e-14, 50);
    let reference = find_entropic_alpha(&f, &f_eq, 60);
    assert!(
        (alpha - 2.0).abs() > 1e-3,
        "mirror-start Newton must not stay at 2, got {alpha}"
    );
    assert!(
        (alpha - reference).abs() < 1e-9,
        "Newton {alpha} vs bisection {reference}"
    );
}

#[test]
fn entropic_alpha_is_two_close_to_equilibrium() {
    let f = perturbed_state(1e-7);
    let f_eq = local_equilibrium(&f);
    assert_eq!(find_entropic_alpha(&f, &f_eq, 50), 2.0);
    assert_eq!(find_entropic_alpha_newton(&f, &f_eq, 20), 2.0);
    // Exactly at equilibrium as well.
    let f_eq = d2q9_equilibrium(1.0, [0.03, -0.02]);
    assert_eq!(find_entropic_alpha(&f_eq, &f_eq, 50), 2.0);
}

#[test]
fn entropic_alpha_falls_back_to_two_not_one_without_a_root() {
    // f is the entropy minimiser, so the (polynomial) "equilibrium" handed in
    // has higher entropy: phi(1) > 0 and there is no entropic correction.
    // The fall-back must be BGK over-relaxation (alpha = 2), not alpha = 1,
    // which would double the viscosity.
    let u = [0.1, 0.06];
    let f = d2q9_equilibrium(1.0, u);
    let f_eq = d2q9_polynomial_equilibrium(1.0, u);
    assert!(h_function(&f_eq) > h_function(&f), "test setup");
    assert!(!is_near_equilibrium(&f, &f_eq), "test setup");
    assert_eq!(find_entropic_alpha(&f, &f_eq, 50), 2.0);
    assert_eq!(find_entropic_alpha_newton(&f, &f_eq, 20), 2.0);
}

// ---------------------------------------------------------------------------
// Viscosity enters through beta
// ---------------------------------------------------------------------------

#[test]
fn relaxation_beta_follows_the_viscosity() {
    // tau = 3 nu + 1/2, beta = 1/(2 tau)
    assert!((relaxation_beta(1.0 / 6.0) - 0.5).abs() < 1e-15);
    assert!((relaxation_beta(0.05) - 1.0 / 1.3).abs() < 1e-15);
    assert!((relaxation_beta(1e-9) - 1.0).abs() < 1e-8);
    for &invalid in &[0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert_eq!(relaxation_beta(invalid), 1.0, "nu = {invalid}");
    }
    assert!((EntropicD2Q9::new(2, 2, 0.05).beta() - 1.0 / 1.3).abs() < 1e-15);
    assert!((KbcD2Q9::new(2, 2, 0.05).beta() - 1.0 / 1.3).abs() < 1e-15);
    assert!((EntropyOverRelaxation::new(0.05).beta() - 1.0 / 1.3).abs() < 1e-15);
}

#[test]
fn try_new_validates_the_viscosity_and_grid() {
    assert!(EntropicD2Q9::try_new(4, 4, 0.1).is_ok());
    assert!(KbcD2Q9::try_new(4, 4, 0.1).is_ok());
    assert!(EntropyOverRelaxation::try_new(0.1).is_ok());

    // nu = 0 means tau = 1/2, which Error::check_tau rejects.
    assert!(matches!(
        EntropicD2Q9::try_new(4, 4, 0.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    assert!(matches!(
        KbcD2Q9::try_new(4, 4, 0.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    assert!(matches!(
        EntropyOverRelaxation::try_new(0.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    for &invalid in &[-0.05, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            matches!(
                EntropicD2Q9::try_new(4, 4, invalid),
                Err(Error::ParameterOutOfBounds { .. })
            ),
            "nu = {invalid}"
        );
        assert!(KbcD2Q9::try_new(4, 4, invalid).is_err(), "nu = {invalid}");
        assert!(
            EntropyOverRelaxation::try_new(invalid).is_err(),
            "nu = {invalid}"
        );
    }
    assert!(matches!(
        EntropicD2Q9::try_new(0, 4, 0.1),
        Err(Error::InvalidGridDimensions { .. })
    ));
    assert!(matches!(
        KbcD2Q9::try_new(4, 0, 0.1),
        Err(Error::InvalidGridDimensions { .. })
    ));
}

/// Single-cell solver holding `f`, collided once (no streaming).
fn collide_single_cell(f: &[f64; 9], nu: f64, newton: bool) -> [f64; 9] {
    let mut solver = EntropicD2Q9::new(1, 1, nu);
    solver.f[0] = *f;
    if newton {
        solver.collide_newton();
    } else {
        solver.collide();
    }
    solver.f[0]
}

fn collide_single_cell_kbc(f: &[f64; 9], nu: f64) -> [f64; 9] {
    let mut solver = KbcD2Q9::new(1, 1, nu);
    solver.f[0] = *f;
    solver.collide_kbc();
    solver.f[0]
}

/// Close to equilibrium every variant is exactly BGK with tau = 3 nu + 1/2:
/// the non-equilibrium part is multiplied by (1 - 1/tau).
#[test]
fn collisions_reduce_to_bgk_with_the_requested_viscosity_near_equilibrium() {
    let f = perturbed_state(1e-7);
    let f_eq = local_equilibrium(&f);
    for &nu in &[0.02, 0.05, 1.0 / 6.0, 0.3] {
        let tau = 3.0 * nu + 0.5;
        let eor = EntropyOverRelaxation::new(nu);
        let (rho, u) = moments(&f);
        let results = [
            ("collide", collide_single_cell(&f, nu, false)),
            ("collide_newton", collide_single_cell(&f, nu, true)),
            ("kbc", collide_single_cell_kbc(&f, nu)),
            ("eor", eor.collide(&f, rho, u)),
        ];
        for (name, post) in &results {
            for i in 0..9 {
                let expected = f_eq[i] + (1.0 - 1.0 / tau) * (f[i] - f_eq[i]);
                assert!(
                    (post[i] - expected).abs() < 1e-14,
                    "{name}, nu = {nu}, i = {i}: {} vs BGK {expected}",
                    post[i]
                );
            }
        }
    }
}

/// Far from equilibrium the collisions still conserve mass and momentum and
/// obey the discrete H-theorem, for every viscosity including the inviscid
/// limit beta -> 1.
#[test]
fn collisions_conserve_moments_and_satisfy_the_h_theorem_far_from_equilibrium() {
    let f = perturbed_state(1.0);
    let (rho0, u0) = moments(&f);
    for &nu in &[1e-6, 0.01, 0.05, 1.0 / 6.0, 0.4] {
        let eor = EntropyOverRelaxation::new(nu);
        let results = [
            ("collide", collide_single_cell(&f, nu, false)),
            ("collide_newton", collide_single_cell(&f, nu, true)),
            ("kbc", collide_single_cell_kbc(&f, nu)),
            ("eor", eor.collide(&f, rho0, u0)),
        ];
        for (name, post) in &results {
            let (rho, u) = moments(post);
            assert!((rho - rho0).abs() < 1e-13, "{name}, nu = {nu}: rho {rho}");
            assert!(
                (u[0] - u0[0]).abs() < 1e-13,
                "{name}, nu = {nu}: ux {}",
                u[0]
            );
            assert!(
                (u[1] - u0[1]).abs() < 1e-13,
                "{name}, nu = {nu}: uy {}",
                u[1]
            );
            // The entropic and EOR collisions stay on the admissible segment
            // between f and the (positive) over-relaxed state, so they obey the
            // H-theorem and keep every population positive; the KBC
            // stabiliser only minimises a quadratic entropy estimate and does
            // not promise that for such extreme states.
            if *name != "kbc" {
                assert!(
                    post.iter().all(|&v| v > 0.0),
                    "{name}, nu = {nu}: negative f"
                );
                let gain = entropy_gain(&f, post);
                assert!(gain <= 1e-12, "{name}, nu = {nu}: H grew by {gain}");
            }
        }
    }
}

/// The entropic collision is not simply BGK far from equilibrium: alpha
/// deviates from 2, so the result differs from f + 2 beta (f_eq - f).
#[test]
fn entropic_collision_differs_from_bgk_far_from_equilibrium() {
    let f = perturbed_state(1.0);
    let f_eq = local_equilibrium(&f);
    let nu = 0.05;
    let beta = relaxation_beta(nu);
    let post = collide_single_cell(&f, nu, false);
    let max_dev = (0..9)
        .map(|i| (post[i] - (f[i] + 2.0 * beta * (f_eq[i] - f[i]))).abs())
        .fold(0.0_f64, f64::max);
    assert!(max_dev > 1e-3, "entropic result equals BGK (dev {max_dev})");
}

// ---------------------------------------------------------------------------
// Shear-wave decay: the macroscopic viscosity is the requested one
// ---------------------------------------------------------------------------

const SHEAR_NX: usize = 2;
const SHEAR_NY: usize = 64;

fn initial_shear_wave(amplitude: f64) -> Vec<[f64; 9]> {
    let mut cells = Vec::with_capacity(SHEAR_NX * SHEAR_NY);
    for y in 0..SHEAR_NY {
        let ux = amplitude * (2.0 * PI * y as f64 / SHEAR_NY as f64).sin();
        let f = EntropicD2Q9::equilibrium(1.0, [ux, 0.0]);
        for _ in 0..SHEAR_NX {
            cells.push(f);
        }
    }
    cells
}

/// First Fourier amplitude of u_x(y) = A sin(2 pi y / ny).
fn shear_wave_amplitude(cells: &[[f64; 9]]) -> f64 {
    let mut acc = 0.0;
    for y in 0..SHEAR_NY {
        let (_, u) = moments(&cells[y * SHEAR_NX]);
        acc += u[0] * (2.0 * PI * y as f64 / SHEAR_NY as f64).sin();
    }
    2.0 * acc / SHEAR_NY as f64
}

/// Kinematic viscosity inferred from the decay A(T) = A(0) exp(-nu k^2 T).
fn viscosity_from_decay(a_start: f64, a_end: f64, steps: usize) -> f64 {
    let k = 2.0 * PI / SHEAR_NY as f64;
    -(a_end / a_start).ln() / (k * k * steps as f64)
}

const SHEAR_STEPS: usize = 200;
const SHEAR_AMPLITUDE: f64 = 0.01;

/// Pure BGK reference (tau = 3 nu + 1/2) on the same lattice and initial
/// condition: measures the lattice's own decay rate, including its O(k^2)
/// discretisation error.
fn bgk_shear_viscosity(nu: f64) -> f64 {
    let omega = 1.0 / (3.0 * nu + 0.5);
    let mut solver = EntropicD2Q9::new(SHEAR_NX, SHEAR_NY, nu);
    solver.f = initial_shear_wave(SHEAR_AMPLITUDE);
    let a0 = shear_wave_amplitude(&solver.f);
    for _ in 0..SHEAR_STEPS {
        for cell in solver.f.iter_mut() {
            let (rho, u) = moments(cell);
            let f_eq = EntropicD2Q9::equilibrium(rho, u);
            for i in 0..9 {
                cell[i] += omega * (f_eq[i] - cell[i]);
            }
        }
        solver.stream();
    }
    viscosity_from_decay(a0, shear_wave_amplitude(&solver.f), SHEAR_STEPS)
}

fn entropic_shear_viscosity(nu: f64, newton: bool) -> f64 {
    let mut solver = EntropicD2Q9::new(SHEAR_NX, SHEAR_NY, nu);
    solver.f = initial_shear_wave(SHEAR_AMPLITUDE);
    let a0 = shear_wave_amplitude(&solver.f);
    for _ in 0..SHEAR_STEPS {
        if newton {
            solver.collide_newton();
        } else {
            solver.collide();
        }
        solver.stream();
    }
    viscosity_from_decay(a0, shear_wave_amplitude(&solver.f), SHEAR_STEPS)
}

/// The measured viscosity must match the requested one up to the lattice's
/// own discretisation error (visible in the BGK reference), and the entropic
/// variant must reproduce the BGK reference almost exactly.
fn assert_matches_requested_viscosity(name: &str, nu: f64, measured: f64) {
    let reference = bgk_shear_viscosity(nu);
    assert!(
        (reference - nu).abs() < 0.03 * nu,
        "{name}: BGK reference nu = {reference} vs requested {nu}"
    );
    assert!(
        (measured - nu).abs() < 0.02 * nu + (reference - nu).abs(),
        "{name}: measured nu = {measured}, requested {nu} (BGK lattice value {reference})"
    );
    assert!(
        (measured - reference).abs() < 0.005 * nu,
        "{name}: measured nu = {measured} deviates from BGK {reference}"
    );
}

#[test]
fn entropic_shear_wave_decays_with_the_requested_viscosity() {
    let nu = 0.05;
    assert_matches_requested_viscosity("bisection", nu, entropic_shear_viscosity(nu, false));
}

#[test]
fn entropic_newton_shear_wave_decays_with_the_requested_viscosity() {
    for &nu in &[0.02, 0.05, 0.08] {
        assert_matches_requested_viscosity("newton", nu, entropic_shear_viscosity(nu, true));
    }
}

#[test]
fn shear_wave_decay_rate_scales_with_the_viscosity() {
    let low = entropic_shear_viscosity(0.02, true);
    let high = entropic_shear_viscosity(0.08, true);
    assert!(
        (high / low - 4.0).abs() < 0.2,
        "viscosity ratio {} should be about 4",
        high / low
    );
}

#[test]
fn kbc_shear_wave_decays_with_the_requested_viscosity() {
    let nu = 0.05;
    let mut solver = KbcD2Q9::new(SHEAR_NX, SHEAR_NY, nu);
    solver.f = initial_shear_wave(SHEAR_AMPLITUDE);
    let a0 = shear_wave_amplitude(&solver.f);
    for _ in 0..SHEAR_STEPS {
        solver.step();
    }
    let measured = viscosity_from_decay(a0, shear_wave_amplitude(&solver.f), SHEAR_STEPS);
    assert_matches_requested_viscosity("kbc", nu, measured);
}

#[test]
fn eor_shear_wave_decays_with_the_requested_viscosity() {
    let nu = 0.05;
    let eor = EntropyOverRelaxation::new(nu);
    let mut solver = EntropicD2Q9::new(SHEAR_NX, SHEAR_NY, nu);
    solver.f = initial_shear_wave(SHEAR_AMPLITUDE);
    let a0 = shear_wave_amplitude(&solver.f);
    // EOR acts cell-wise; the solver provides the pull streaming.
    for _ in 0..SHEAR_STEPS {
        for cell in solver.f.iter_mut() {
            let (rho, u) = moments(cell);
            *cell = eor.collide(cell, rho, u);
        }
        solver.stream();
    }
    let measured = viscosity_from_decay(a0, shear_wave_amplitude(&solver.f), SHEAR_STEPS);
    assert_matches_requested_viscosity("eor", nu, measured);
}
