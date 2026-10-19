//! Regression tests for the entropic (ELBM) mirror-state search and collision.
//!
//! `entropic_alpha_d3q19` bisected the wrong way: it moved the bracket away from the root and
//! returned the lower end (`1`) or the upper end (`2`), and `collide_entropic_stabilized` applied
//! `alpha (feq - f)` without `beta`, so the viscosity could not be set. These tests pin the
//! mirror-state condition `H(f + alpha (feq - f)) = H(f)`, the BGK limit near equilibrium, the
//! positivity cap and the viscosity-controlled collision.

use super::test_support::{
    assert_conserved, double_shear_layer_survival, perturbed_populations,
    perturbed_rest_populations,
};
use super::*;
use crate::lattice::{D3Q19_VELOCITIES, D3Q19_WEIGHTS};

fn mirror(f: &[f64; 19], feq: &[f64; 19], alpha: f64) -> [f64; 19] {
    std::array::from_fn(|i| f[i] + alpha * (feq[i] - f[i]))
}

/// `(f, feq)` of a deterministic perturbed state with the equilibrium of its own `(rho, u)`.
fn state(seed: u64, amplitude: f64) -> ([f64; 19], [f64; 19]) {
    let f = perturbed_populations(seed, amplitude);
    let (rho, u) = compute_macroscopic(&f);
    (f, compute_equilibrium_d3q19(rho, u))
}

#[test]
fn mirror_state_has_the_entropy_of_the_input() {
    let mut interior = 0;
    let mut beyond_two = 0;
    for seed in 1..=60 {
        let (f, feq) = state(seed, 0.3);
        let alpha = entropic_alpha_d3q19(&f, &feq);
        assert!(alpha > 1.0 && alpha <= 2.0, "seed {seed}: alpha = {alpha}");
        let h_gap = discrete_h_function(&mirror(&f, &feq, alpha)) - discrete_h_function(&f);
        if alpha < 2.0 {
            interior += 1;
            // The root is not one of the bracket ends the inverted bisection used to return.
            assert!(alpha > 1.9, "seed {seed}: alpha = {alpha}");
            assert!(
                h_gap.abs() < 1e-10,
                "seed {seed}: |H(mirror) - H(f)| = {}",
                h_gap.abs()
            );
        } else {
            beyond_two += 1;
            // The root lies beyond the BGK limit: the entropy has not returned to H(f) yet.
            assert!(
                h_gap <= 1e-12,
                "seed {seed}: H(f + 2 (feq - f)) - H(f) = {h_gap}"
            );
        }
    }
    assert!(
        interior >= 10,
        "only {interior} states with an interior root"
    );
    assert!(beyond_two >= 10, "only {beyond_two} states clamped to 2");
}

#[test]
fn root_is_found_to_round_off() {
    for seed in 1..=60 {
        let (f, feq) = state(seed, 0.3);
        let alpha = entropic_alpha_d3q19(&f, &feq);
        if alpha >= 2.0 {
            continue;
        }
        let h = |a: f64| discrete_h_function(&mirror(&f, &feq, a)) - discrete_h_function(&f);
        assert!(h(alpha).abs() < 1e-13, "seed {seed}: gap {}", h(alpha));
        // H is negative strictly before the root and positive strictly after it.
        assert!(
            h(alpha - 1e-3) < 0.0 && h(alpha + 1e-3) > 0.0,
            "seed {seed}"
        );
    }
}

#[test]
fn near_equilibrium_populations_give_alpha_close_to_two() {
    // For the entropy minimiser alpha* = 2 - O(perturbation), so the deviation shrinks with the
    // non-equilibrium. Around the rest state the polynomial equilibrium is the minimiser up to
    // O(u^3) ...
    for seed in 1..=20 {
        for amplitude in [1e-4, 1e-3] {
            let f = perturbed_rest_populations(seed, amplitude);
            let (rho, u) = compute_macroscopic(&f);
            let alpha = entropic_alpha_d3q19(&f, &compute_equilibrium_d3q19(rho, u));
            assert!(
                (alpha - 2.0).abs() < 0.5 * amplitude,
                "rest state, seed {seed}, perturbation {amplitude}: alpha = {alpha}"
            );
        }
    }
    // ... and for a moving planar flow the entropic (product-form) equilibrium is the exact
    // minimiser. The low-Mach polynomial equilibrium of a moving fluid is not: it leaves alpha*
    // at 2 + O(u^3 / perturbation), which exits (1, 2] and is clamped to the BGK value.
    let base = entropic_equilibrium_d3q19(1.0, [0.05, -0.03, 0.0]);
    for seed in 1..=20 {
        for amplitude in [1e-4, 1e-3] {
            let noise = perturbed_rest_populations(seed, 1.0);
            let f: [f64; 19] = std::array::from_fn(|i| {
                base[i] * (1.0 + amplitude * (noise[i] / D3Q19_WEIGHTS[i] - 1.0))
            });
            let (rho, u) = compute_macroscopic(&f);
            let alpha = entropic_alpha_d3q19(&f, &entropic_equilibrium_d3q19(rho, u));
            assert!(
                (alpha - 2.0).abs() < 0.5 * amplitude,
                "moving state, seed {seed}, perturbation {amplitude}: alpha = {alpha}"
            );
        }
    }
}

#[test]
fn equilibrium_gives_the_bgk_value() {
    let feq = compute_equilibrium_d3q19(1.3, [0.04, -0.02, 0.0]);
    assert_eq!(entropic_alpha_d3q19(&feq, &feq), 2.0);
}

#[test]
fn populations_that_are_not_positive_fall_back_to_the_bgk_value() {
    let feq = compute_equilibrium_d3q19(1.0, [0.02, 0.0, 0.0]);
    let (f, _) = state(4, 0.2);
    let mut negative_f = f;
    negative_f[5] = -0.01;
    assert_eq!(entropic_alpha_d3q19(&negative_f, &feq), 2.0);
    let mut zero_f = f;
    zero_f[2] = 0.0;
    assert_eq!(entropic_alpha_d3q19(&zero_f, &feq), 2.0);
    let mut negative_feq = feq;
    negative_feq[7] = -0.1;
    assert_eq!(entropic_alpha_d3q19(&f, &negative_feq), 2.0);
    let mut nan_f = f;
    nan_f[0] = f64::NAN;
    assert_eq!(entropic_alpha_d3q19(&nan_f, &feq), 2.0);
}

#[test]
fn positivity_limit_caps_alpha() {
    // Rest-state weights, except that population 5 sits far above its equilibrium value, so the
    // line f + alpha (feq - f) reaches zero in that population at alpha = 0.1 / 0.07 < 2.
    let mut f = D3Q19_WEIGHTS;
    let mut feq = D3Q19_WEIGHTS;
    f[5] = 0.1;
    feq[5] = 0.03;
    let limit = 0.1 / 0.07;
    let alpha = entropic_alpha_d3q19(&f, &feq);
    assert!(
        alpha > 1.0 && alpha <= limit + 1e-12,
        "alpha = {alpha}, limit = {limit}"
    );
    assert!(mirror(&f, &feq, alpha).iter().all(|&x| x >= -1e-15));
    // The entropy is still below H(f) at the limit, so the limit itself is returned.
    assert!((alpha - limit).abs() < 1e-12);
}

/// Index of the velocity `(cx, cy, -cz)` for every D3Q19 direction.
fn z_mirror() -> [usize; 19] {
    std::array::from_fn(|i| {
        let c = D3Q19_VELOCITIES[i];
        (0..19)
            .find(|&j| {
                let d = D3Q19_VELOCITIES[j];
                d[0] == c[0] && d[1] == c[1] && d[2] == -c[2]
            })
            .unwrap_or(i)
    })
}

/// Non-equilibrium state of a flow in the xy plane: noise symmetric under `z -> -z`, so that
/// `u_z = 0` and the product-form equilibrium is exactly conservative.
fn planar_state(seed: u64, amplitude: f64) -> [f64; 19] {
    let noisy = perturbed_populations(seed, amplitude);
    let reference = compute_equilibrium_d3q19(1.15, [0.03, -0.02, 0.025]);
    let base = compute_equilibrium_d3q19(1.1, [0.03, -0.02, 0.0]);
    let mirrored = z_mirror();
    std::array::from_fn(|i| {
        let j = mirrored[i];
        base[i] * 0.5 * (noisy[i] / reference[i] + noisy[j] / reference[j])
    })
}

#[test]
fn collision_matches_the_explicit_elbm_update() {
    for seed in 1..=10 {
        let f = planar_state(seed, 0.2);
        let (rho, u) = compute_macroscopic(&f);
        for omega in [0.6, 1.2, 1.9] {
            let feq = entropic_equilibrium_d3q19(rho, u);
            let alpha = entropic_alpha_d3q19(&f, &feq);
            let out = collide_entropic_stabilized(&f, rho, u, omega);
            for i in 0..19 {
                let expected = f[i] + alpha * 0.5 * omega * (feq[i] - f[i]);
                assert!(
                    (out[i] - expected).abs() < 1e-15,
                    "seed {seed}, omega {omega}"
                );
            }
        }
    }
}

#[test]
fn collision_conserves_mass_and_momentum_for_planar_flows() {
    for seed in 1..=10 {
        let f = planar_state(seed, 0.2);
        let (rho, u) = compute_macroscopic(&f);
        assert!(u[2].abs() < 1e-14, "planar state has u_z = {}", u[2]);
        for omega in [0.6, 1.2, 1.9] {
            let out = collide_entropic_stabilized(&f, rho, u, omega);
            assert_conserved(&f, &out, 1e-13, &format!("seed {seed}, omega {omega}"));
        }
    }
}

#[test]
fn collision_never_increases_the_entropy_for_beta_up_to_one() {
    for seed in 1..=30 {
        let f = planar_state(seed, 0.3);
        let (rho, u) = compute_macroscopic(&f);
        for omega in [0.5, 1.0, 1.7, 2.0] {
            let out = collide_entropic_stabilized(&f, rho, u, omega);
            assert!(
                discrete_h_function(&out) <= discrete_h_function(&f) + 1e-12,
                "seed {seed}, omega {omega}"
            );
            assert!(
                out.iter().all(|&x| x > 0.0),
                "seed {seed}, omega {omega}: negative population"
            );
        }
    }
}

#[test]
fn elbm_survives_an_under_resolved_shear_layer_where_bgk_breaks_down() {
    // The entropic collision adapts alpha (and with it the dissipation) locally: on a double
    // shear layer far thinner than the lattice spacing at nu = 1e-4 BGK produces negative
    // populations after ~230 steps whereas ELBM stays positive.
    let (n, nu, u0, steps) = (12, 1e-4, 0.1, 320);
    let omega = 1.0 / (3.0 * nu + 0.5);
    let bgk = BgkCollision::from_viscosity(nu);
    let bgk_run = double_shear_layer_survival(|f, rho, u| bgk.collide(f, rho, u), n, u0, steps);
    assert!(
        bgk_run.completed_steps < steps,
        "BGK completed {} steps without a negative population",
        bgk_run.completed_steps
    );
    let elbm_run = double_shear_layer_survival(
        |f, rho, u| collide_entropic_stabilized(f, rho, u, omega),
        n,
        u0,
        steps,
    );
    assert_eq!(elbm_run.completed_steps, steps);
    assert!(elbm_run.mass_drift.abs() < 1e-8, "{}", elbm_run.mass_drift);
    for (k, drift) in elbm_run.momentum_drift.iter().enumerate() {
        assert!(drift.abs() < 1e-8, "momentum component {k}: {drift}");
    }
}

#[test]
fn collision_is_bgk_with_the_requested_rate_close_to_equilibrium() {
    // alpha -> 2 as f -> feq, so beta = omega / 2 turns the update into BGK with rate omega
    // towards the entropic equilibrium and the viscosity nu = cs^2 (1/omega - 1/2) is controlled
    // by omega.
    let (rho, u) = (1.0, [0.05, -0.03, 0.0]);
    for omega in [0.8, 1.5, 1.9] {
        let feq = entropic_equilibrium_d3q19(rho, u);
        let f: [f64; 19] = std::array::from_fn(|i| feq[i] * (1.0 + 1e-6 * (i as f64 - 9.0)));
        let (rho_f, u_f) = compute_macroscopic(&f);
        let feq_f = entropic_equilibrium_d3q19(rho_f, u_f);
        let elbm = collide_entropic_stabilized(&f, rho_f, u_f, omega);
        let alpha = entropic_alpha_d3q19(&f, &feq_f);
        assert!((alpha - 2.0).abs() < 1e-3, "omega {omega}: alpha = {alpha}");
        for i in 0..19 {
            let bgk = f[i] + omega * (feq_f[i] - f[i]);
            assert!(
                (elbm[i] - bgk).abs() < 1e-9,
                "omega {omega}, population {i}: {} vs {bgk}",
                elbm[i]
            );
        }
    }
}
