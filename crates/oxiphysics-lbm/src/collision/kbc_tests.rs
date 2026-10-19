//! Regression tests for the KBC collision operator.
//!
//! The previous operator applied `f + gamma (feq - f)` with `gamma` either `omega_s` or `2`:
//! there was no shear / higher-order split (the mirror state was computed and never used), and
//! `gamma = 2` relaxed the shear stress with `omega = 2` (zero viscosity). These tests pin the
//! decomposition `f - feq = ds + dh`, the closed-form stabiliser, the conservation laws and the
//! viscosity of the shear relaxation.

use super::test_support::{
    assert_conserved, double_shear_layer_survival, mass_and_momentum, perturbed_populations,
    shear_wave_amplitude_ratio, shear_wave_decay,
};
use super::*;

/// Second moment `sum_i x_i c_i c_i` as `[xx, yy, zz, xy, xz, yz]`.
fn second_moment(x: &[f64; 19]) -> [f64; 6] {
    compute_pi1_tensor(x, &[0.0; 19])
}

/// Quadratic entropy estimate `<x|x>` of the deviation from equilibrium.
fn entropy_estimate(f: &[f64; 19], feq: &[f64; 19]) -> f64 {
    let dev: [f64; 19] = std::array::from_fn(|i| f[i] - feq[i]);
    kbc_entropic_product(&dev, &dev, feq)
}

#[test]
fn decomposition_splits_the_non_equilibrium_part_into_shear_and_higher_order() {
    for seed in 1..=10 {
        let f = perturbed_populations(seed, 0.1);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let (delta_s, delta_h) = kbc_decompose_d3q19(&f, &feq);
        for i in 0..19 {
            assert!(
                (delta_s[i] + delta_h[i] - (f[i] - feq[i])).abs() < 1e-15,
                "seed {seed}: ds + dh != f - feq at {i}"
            );
        }
        // Neither part carries mass or momentum.
        for (name, part) in [("ds", &delta_s), ("dh", &delta_h)] {
            for (k, value) in mass_and_momentum(part).iter().enumerate() {
                assert!(
                    value.abs() < 1e-14,
                    "seed {seed}: {name} moment {k}: {value}"
                );
            }
        }
        // The second moment of ds is the traceless stress; the trace stays in dh.
        let pi = compute_pi1_tensor(&f, &feq);
        let mean_normal = (pi[0] + pi[1] + pi[2]) / 3.0;
        let expected_s = [
            pi[0] - mean_normal,
            pi[1] - mean_normal,
            pi[2] - mean_normal,
            pi[3],
            pi[4],
            pi[5],
        ];
        let expected_h = [mean_normal, mean_normal, mean_normal, 0.0, 0.0, 0.0];
        for (k, (&s, &h)) in second_moment(&delta_s)
            .iter()
            .zip(second_moment(&delta_h).iter())
            .enumerate()
        {
            assert!(
                (s - expected_s[k]).abs() < 1e-14,
                "seed {seed}: ds stress {k}"
            );
            assert!(
                (h - expected_h[k]).abs() < 1e-14,
                "seed {seed}: dh stress {k}"
            );
        }
    }
}

#[test]
fn collision_conserves_mass_and_momentum() {
    for nu in [0.3, 0.1, 0.02, 0.001] {
        let kbc = KbcCollision::new(nu);
        for seed in 1..=8 {
            let f = perturbed_populations(seed, 0.1);
            let (rho, u) = compute_macroscopic(&f);
            let out = kbc.collide(&f, rho, u);
            assert_conserved(&f, &out, 1e-14, &format!("nu = {nu}, seed {seed}"));
            let max_change = f
                .iter()
                .zip(out.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            assert!(max_change > 1e-4, "the collision must change the state");
        }
    }
}

#[test]
fn equilibrium_is_a_fixed_point() {
    let (rho, u) = (0.9, [0.04, -0.01, 0.02]);
    let feq = compute_equilibrium_d3q19(rho, u);
    for nu in [0.3, 0.02, 0.0005] {
        let out = KbcCollision::new(nu).collide(&feq, rho, u);
        for i in 0..19 {
            assert!((out[i] - feq[i]).abs() < 1e-15, "nu = {nu}, population {i}");
        }
    }
}

#[test]
fn vanishing_higher_order_part_gives_plain_bgk() {
    let (rho, u) = (1.1, [0.03, -0.02, 0.01]);
    let feq = compute_equilibrium_d3q19(rho, u);
    // A traceless stress and nothing else: dh = 0.
    let delta_s = kbc_shear_part_d3q19(&[2e-3, -1e-3, -1e-3, 5e-4, -3e-4, 2e-4]);
    let f: [f64; 19] = std::array::from_fn(|i| feq[i] + delta_s[i]);
    let (_, delta_h) = kbc_decompose_d3q19(&f, &feq);
    assert!(delta_h.iter().all(|v| v.abs() < 1e-14), "{delta_h:?}");
    let (rho_f, u_f) = compute_macroscopic(&f);
    assert!((rho_f - rho).abs() < 1e-14);
    for nu in [0.3, 0.05, 0.002] {
        let kbc = KbcCollision::new(nu);
        let bgk = collide_bgk(&f, rho_f, u_f, kbc.omega_s);
        let out = kbc.collide(&f, rho_f, u_f);
        for i in 0..19 {
            assert!(
                (out[i] - bgk[i]).abs() < 1e-14,
                "nu = {nu}, population {i}: {} vs {}",
                out[i],
                bgk[i]
            );
        }
    }
}

#[test]
fn round_off_level_higher_order_part_is_ignored() {
    // The ratio <ds|dh> / <dh|dh> does not vanish with dh, so a dh made of round-off noise must
    // not be allowed to steer gamma: the state is treated as having no higher-order part.
    let (rho, u) = (1.1, [0.03, -0.02, 0.01]);
    let feq = compute_equilibrium_d3q19(rho, u);
    let delta_s = kbc_shear_part_d3q19(&[2e-3, -1e-3, -1e-3, 5e-4, -3e-4, 2e-4]);
    let f: [f64; 19] =
        std::array::from_fn(|i| feq[i] + delta_s[i] + 1e-13 * feq[i] * (i as f64 * 1.7).sin());
    let (rho_f, u_f) = compute_macroscopic(&f);
    let kbc = KbcCollision::new(0.002);
    assert_eq!(kbc.stabilizer(&f, rho_f, u_f), 2.0);
    let bgk = collide_bgk(&f, rho_f, u_f, kbc.omega_s);
    let out = kbc.collide(&f, rho_f, u_f);
    for (a, b) in out.iter().zip(bgk.iter()) {
        assert!((a - b).abs() < 1e-14, "{a} vs {b}");
    }
}

#[test]
fn unit_shear_rate_is_plain_bgk_for_every_state() {
    // omega_s = 1 means beta = 1/2 and gamma = 2 whatever dh is: f* = f - (f - feq) = feq.
    let kbc = KbcCollision::new(1.0 / 6.0);
    assert!((kbc.omega_s - 1.0).abs() < 1e-15);
    for seed in 1..=10 {
        let f = perturbed_populations(seed, 0.2);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        assert!((kbc.stabilizer(&f, rho, u) - 2.0).abs() < 1e-12);
        let out = kbc.collide(&f, rho, u);
        for i in 0..19 {
            assert!(
                (out[i] - feq[i]).abs() < 1e-13,
                "seed {seed}, population {i}"
            );
        }
    }
}

#[test]
fn stabilizer_follows_the_closed_form_for_synthetic_parts() {
    let feq = compute_equilibrium_d3q19(1.0, [0.0; 3]);
    let beta = 0.8;
    let unit = |k: usize| -> [f64; 19] { std::array::from_fn(|i| if i == k { 0.01 } else { 0.0 }) };
    // Orthogonal parts: gamma = 1 / beta (the higher-order part is reset to equilibrium).
    let gamma = kbc_stabilizer(beta, &unit(3), &unit(9), &feq);
    assert!((gamma - 1.0 / beta).abs() < 1e-14);
    // Parallel parts: gamma = 2 / beta - 2.
    let gamma = kbc_stabilizer(beta, &unit(9), &unit(9), &feq);
    assert!((gamma - (2.0 / beta - 2.0)).abs() < 1e-14);
    // Anti-parallel parts: gamma = 2, whatever beta is.
    let negative: [f64; 19] = std::array::from_fn(|i| -unit(9)[i]);
    let gamma = kbc_stabilizer(beta, &unit(9), &negative, &feq);
    assert!((gamma - 2.0).abs() < 1e-14);
}

#[test]
fn stabilizer_falls_back_to_bgk_for_degenerate_input() {
    let feq = compute_equilibrium_d3q19(1.0, [0.01, 0.0, 0.0]);
    let zero = [0.0_f64; 19];
    let part: [f64; 19] = std::array::from_fn(|i| 1e-3 * (i as f64 - 9.0));
    // dh = 0: the value of gamma is immaterial and defined as the BGK value.
    assert_eq!(kbc_stabilizer(0.4, &part, &zero, &feq), 2.0);
    // Non-positive or NaN beta.
    assert_eq!(kbc_stabilizer(0.0, &part, &part, &feq), 2.0);
    assert_eq!(kbc_stabilizer(-0.2, &part, &part, &feq), 2.0);
    assert_eq!(kbc_stabilizer(f64::NAN, &part, &part, &feq), 2.0);
    // The entropic product is undefined without a positive equilibrium.
    let mut bad = feq;
    bad[4] = -0.01;
    assert_eq!(kbc_stabilizer(0.4, &part, &part, &bad), 2.0);
}

#[test]
fn stabilizer_minimises_the_quadratic_entropy_estimate() {
    for nu in [0.05, 0.005] {
        let kbc = KbcCollision::new(nu);
        let beta = 0.5 * kbc.omega_s;
        for seed in 1..=10 {
            let f = perturbed_populations(seed, 0.15);
            let (rho, u) = compute_macroscopic(&f);
            let feq = compute_equilibrium_d3q19(rho, u);
            let (delta_s, delta_h) = kbc_decompose_d3q19(&f, &feq);
            let post = |gamma: f64| -> [f64; 19] {
                std::array::from_fn(|i| f[i] - beta * (2.0 * delta_s[i] + gamma * delta_h[i]))
            };
            let gamma = kbc.stabilizer(&f, rho, u);
            // The estimate is a parabola in gamma and gamma is its vertex.
            let (below, above) = (
                entropy_estimate(&post(gamma - 0.5), &feq),
                entropy_estimate(&post(gamma + 0.5), &feq),
            );
            assert!(
                (above - below).abs() < 1e-9 * (above + below),
                "nu = {nu}, seed {seed}: not the vertex ({below} vs {above})"
            );
            let at_gamma = entropy_estimate(&post(gamma), &feq);
            assert!(at_gamma < below && at_gamma < above);
            // In particular it beats BGK (gamma = 2) and the KBC output is what `collide` returns.
            assert!(at_gamma <= entropy_estimate(&post(2.0), &feq) * (1.0 + 1e-12));
            let out = kbc.collide(&f, rho, u);
            let expected = post(gamma);
            for (a, b) in out.iter().zip(expected.iter()) {
                assert!((a - b).abs() < 1e-15);
            }
        }
    }
}

#[test]
fn stabilizer_reduces_the_entropy_estimate_below_bgk_for_strong_non_equilibrium() {
    let nu = 0.0005;
    let kbc = KbcCollision::new(nu);
    let mut strictly_better = 0;
    for seed in 1..=20 {
        let f = perturbed_populations(seed, 0.4);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let bgk = collide_bgk(&f, rho, u, kbc.omega_s);
        let out = kbc.collide(&f, rho, u);
        let (q_kbc, q_bgk) = (entropy_estimate(&out, &feq), entropy_estimate(&bgk, &feq));
        assert!(
            q_kbc <= q_bgk * (1.0 + 1e-12),
            "seed {seed}: {q_kbc} > {q_bgk}"
        );
        if q_kbc < 0.999 * q_bgk {
            strictly_better += 1;
        }
    }
    assert!(
        strictly_better >= 10,
        "only {strictly_better} of 20 improved"
    );
}

#[test]
fn kbc_survives_an_under_resolved_shear_layer_where_bgk_breaks_down() {
    // Double shear layer far thinner than the lattice spacing at nu = 1e-4: BGK produces negative
    // populations after ~230 steps, the entropic stabiliser keeps every population positive and
    // conserves mass and momentum to round-off.
    let (n, nu, u0, steps) = (12, 1e-4, 0.1, 320);
    let bgk = BgkCollision::from_viscosity(nu);
    let kbc = KbcCollision::new(nu);
    let bgk_run = double_shear_layer_survival(|f, rho, u| bgk.collide(f, rho, u), n, u0, steps);
    assert!(
        bgk_run.completed_steps < steps,
        "BGK completed {} steps without a negative population",
        bgk_run.completed_steps
    );
    let kbc_run = double_shear_layer_survival(|f, rho, u| kbc.collide(f, rho, u), n, u0, steps);
    assert_eq!(kbc_run.completed_steps, steps);
    assert!(kbc_run.mass_drift.abs() < 1e-10, "{}", kbc_run.mass_drift);
    for (k, drift) in kbc_run.momentum_drift.iter().enumerate() {
        assert!(drift.abs() < 1e-10, "momentum component {k}: {drift}");
    }
}

#[test]
fn shear_wave_decays_with_the_requested_viscosity() {
    let (nu, ny, steps, amplitude) = (0.02, 32, 400, 0.02);
    let kbc = KbcCollision::new(nu);
    assert!((kbc.effective_viscosity() - nu).abs() < 1e-14);
    let ratio =
        shear_wave_amplitude_ratio(|f, rho, u| kbc.collide(f, rho, u), ny, amplitude, steps);
    let expected = shear_wave_decay(nu, ny, steps);
    assert!(
        (ratio - expected).abs() < 0.01,
        "decay {ratio} vs exp(-nu k^2 t) = {expected}"
    );
    // The previous operator relaxed the shear stress at omega = 2 (nu = 0) and barely damped it.
    assert!(ratio < 0.8, "the wave must decay, ratio {ratio}");
    // BGK with the same viscosity decays the same way.
    let bgk = BgkCollision::from_viscosity(nu);
    let bgk_ratio =
        shear_wave_amplitude_ratio(|f, rho, u| bgk.collide(f, rho, u), ny, amplitude, steps);
    assert!(
        (ratio - bgk_ratio).abs() < 2e-3,
        "KBC {ratio} vs BGK {bgk_ratio}"
    );
}
