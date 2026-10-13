//! Regression tests for the central-moment and cumulant collision operators.
//!
//! Both operators used to scatter moment `k` onto population `k` (central moments) or to
//! reconstruct from an incomplete set of moments (cumulants), so neither conserved mass or
//! momentum. They are now built on the exact inverse central-moment transform; these tests pin
//! the conservation laws, the relaxation of the stress moments and the defining properties of
//! each operator.

use super::*;
use crate::lattice::D3Q19_VELOCITIES;

/// Deterministic randomly perturbed populations around a moving equilibrium.
fn perturbed_populations(seed: u64, amplitude: f64) -> [f64; 19] {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let feq = compute_equilibrium_d3q19(1.15, [0.03, -0.02, 0.025]);
    std::array::from_fn(|i| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
        feq[i] * (1.0 + amplitude * (2.0 * unit - 1.0))
    })
}

/// Mass and momentum `[rho, jx, jy, jz]`.
fn mass_and_momentum(f: &[f64; 19]) -> [f64; 4] {
    let mut m = [0.0_f64; 4];
    for (&fi, c) in f.iter().zip(D3Q19_VELOCITIES.iter()) {
        m[0] += fi;
        m[1] += fi * f64::from(c[0]);
        m[2] += fi * f64::from(c[1]);
        m[3] += fi * f64::from(c[2]);
    }
    m
}

fn assert_conserved(before: &[f64; 19], after: &[f64; 19], context: &str) {
    let a = mass_and_momentum(before);
    let b = mass_and_momentum(after);
    for (c, name) in ["rho", "jx", "jy", "jz"].iter().enumerate() {
        assert!(
            (a[c] - b[c]).abs() < 1e-12,
            "{context}: {name} changed from {} to {}",
            a[c],
            b[c]
        );
    }
}

/// Non-equilibrium stress `Pi_ab = sum_i (f_i - feq_i) c_a c_b` of a population set.
fn neq_stress(f: &[f64; 19], rho: f64, u: [f64; 3]) -> [[f64; 3]; 3] {
    non_equilibrium_stress(f, &compute_equilibrium_d3q19(rho, u))
}

fn trace(pi: &[[f64; 3]; 3]) -> f64 {
    pi[0][0] + pi[1][1] + pi[2][2]
}

#[test]
fn central_moment_collision_conserves_mass_and_momentum() {
    for seed in 1..=8 {
        let f = perturbed_populations(seed, 0.2);
        let (rho, u) = compute_macroscopic(&f);
        for nu in [0.02, 0.1, 1.0 / 6.0, 0.4] {
            let out = CentralMomentCollision::from_viscosity(nu).collide(&f, rho, u);
            assert_conserved(&f, &out, &format!("seed {seed}, nu {nu}"));
        }
    }
}

#[test]
fn central_moment_collision_with_uniform_rate_equals_bgk() {
    for seed in 1..=6 {
        let f = perturbed_populations(seed, 0.2);
        let (rho, u) = compute_macroscopic(&f);
        for omega in [0.6, 1.0, 1.45, 1.9] {
            let operator = CentralMomentCollision { s: [omega; 19] };
            let out = operator.collide(&f, rho, u);
            let bgk = collide_bgk(&f, rho, u, omega);
            for i in 0..19 {
                assert!(
                    (out[i] - bgk[i]).abs() < 1e-12,
                    "seed {seed}, omega {omega}, population {i}: {} vs BGK {}",
                    out[i],
                    bgk[i]
                );
            }
        }
    }
}

#[test]
fn central_moment_collision_keeps_the_equilibrium_fixed() {
    let rho = 1.3;
    let u = [0.04, -0.03, 0.02];
    let feq = compute_equilibrium_d3q19(rho, u);
    let out = CentralMomentCollision::from_viscosity(0.07).collide(&feq, rho, u);
    for i in 0..19 {
        assert!((out[i] - feq[i]).abs() < 1e-14, "population {i}");
    }
}

#[test]
fn central_moment_collision_relaxes_each_moment_with_its_own_rate() {
    let f = perturbed_populations(5, 0.25);
    let (rho, u) = compute_macroscopic(&f);
    let nu = 0.12;
    let omega = 1.0 / (3.0 * nu + 0.5);
    let out = CentralMomentCollision::from_viscosity(nu).collide(&f, rho, u);

    // Second-order non-equilibrium stress relaxes by (1 - omega), trace included.
    let before = neq_stress(&f, rho, u);
    let after = neq_stress(&out, rho, u);
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (after[a][b] - (1.0 - omega) * before[a][b]).abs() < 1e-13,
                "Pi[{a}][{b}]: {} vs {}",
                after[a][b],
                (1.0 - omega) * before[a][b]
            );
        }
    }

    // Third- and fourth-order central moments are reset to their equilibrium values.
    let kappa = compute_central_moments_cumulant(&out, u);
    let kappa_eq = compute_central_moments_cumulant(&compute_equilibrium_d3q19(rho, u), u);
    for k in 10..19 {
        assert!(
            (kappa[k] - kappa_eq[k]).abs() < 1e-13,
            "central moment {k}: {} vs equilibrium {}",
            kappa[k],
            kappa_eq[k]
        );
    }
}

#[test]
fn cumulant_collision_conserves_mass_and_momentum() {
    for seed in 1..=8 {
        let f = perturbed_populations(seed, 0.2);
        let (rho, u) = compute_macroscopic(&f);
        for (nu, zeta) in [(0.02, 0.01), (0.1, 0.0), (1.0 / 6.0, 0.3), (0.4, 0.05)] {
            let operator = CumulantCollision::from_viscosities(nu, zeta);
            let out = operator.collide(&f, rho, u);
            assert_conserved(&f, &out, &format!("seed {seed}, nu {nu}, zeta {zeta}"));
        }
        // Unequal shear and bulk rates must not leak mass through the trace mode.
        let operator = CumulantCollision {
            omega: 1.4,
            omega_b: 0.6,
        };
        assert_conserved(
            &f,
            &operator.collide(&f, rho, u),
            &format!("seed {seed}, omega_b"),
        );
    }
}

#[test]
fn cumulant_collision_keeps_the_equilibrium_fixed() {
    for u in [[0.0; 3], [0.04, -0.03, 0.02], [0.1, 0.0, 0.0]] {
        let rho = 1.2;
        let feq = compute_equilibrium_d3q19(rho, u);
        let operator = CumulantCollision::from_viscosities(0.05, 0.02);
        let out = operator.collide(&feq, rho, u);
        for i in 0..19 {
            assert!((out[i] - feq[i]).abs() < 1e-14, "u = {u:?}, population {i}");
        }
    }
}

#[test]
fn cumulant_collision_with_unit_rates_returns_the_equilibrium() {
    let f = perturbed_populations(9, 0.3);
    let (rho, u) = compute_macroscopic(&f);
    let operator = CumulantCollision {
        omega: 1.0,
        omega_b: 1.0,
    };
    let out = operator.collide(&f, rho, u);
    let feq = compute_equilibrium_d3q19(rho, u);
    for i in 0..19 {
        assert!(
            (out[i] - feq[i]).abs() < 1e-13,
            "population {i}: {} vs {}",
            out[i],
            feq[i]
        );
    }
}

#[test]
fn cumulant_collision_separates_shear_and_bulk_relaxation() {
    let f = perturbed_populations(4, 0.25);
    let (rho, u) = compute_macroscopic(&f);
    let (omega, omega_b) = (1.3, 0.7);
    let out = CumulantCollision { omega, omega_b }.collide(&f, rho, u);

    let before = neq_stress(&f, rho, u);
    let after = neq_stress(&out, rho, u);
    // Trace (bulk) relaxes with omega_b.
    assert!(
        (trace(&after) - (1.0 - omega_b) * trace(&before)).abs() < 1e-13,
        "trace: {} vs {}",
        trace(&after),
        (1.0 - omega_b) * trace(&before)
    );
    // Shear and the trace-free normal stresses relax with omega.
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        assert!(
            (after[a][b] - (1.0 - omega) * before[a][b]).abs() < 1e-13,
            "Pi[{a}][{b}]: {} vs {}",
            after[a][b],
            (1.0 - omega) * before[a][b]
        );
    }
    let deviator = |pi: &[[f64; 3]; 3]| (pi[0][0] - pi[1][1], pi[0][0] - pi[2][2]);
    let (d_after, d_before) = (deviator(&after), deviator(&before));
    assert!((d_after.0 - (1.0 - omega) * d_before.0).abs() < 1e-13);
    assert!((d_after.1 - (1.0 - omega) * d_before.1).abs() < 1e-13);
}

#[test]
fn cumulant_collision_resets_higher_cumulants_to_equilibrium() {
    let f = perturbed_populations(2, 0.3);
    let (rho, u) = compute_macroscopic(&f);
    let out = CumulantCollision::from_viscosity(0.1).collide(&f, rho, u);
    let cumulants_out =
        central_moments_to_cumulants_d3q19(&compute_central_moments_cumulant(&out, u));
    let cumulants_eq = central_moments_to_cumulants_d3q19(&compute_central_moments_cumulant(
        &compute_equilibrium_d3q19(rho, u),
        u,
    ));
    for k in 10..19 {
        assert!(
            (cumulants_out[k] - cumulants_eq[k]).abs() < 1e-12,
            "cumulant {k}: {} vs equilibrium {}",
            cumulants_out[k],
            cumulants_eq[k]
        );
    }
}

#[test]
fn cumulant_collision_differs_from_independent_central_moment_relaxation() {
    // With the same stress relaxation the two operators agree to second order but differ in the
    // fourth-order moments: the cumulant reset follows the relaxed second-order moments.
    let f = perturbed_populations(7, 0.3);
    let (rho, u) = compute_macroscopic(&f);
    let nu = 0.05;
    let omega = 1.0 / (3.0 * nu + 0.5);
    let cumulant = CumulantCollision {
        omega,
        omega_b: omega,
    }
    .collide(&f, rho, u);
    let central = CentralMomentCollision::from_viscosity(nu).collide(&f, rho, u);
    let diff = cumulant
        .iter()
        .zip(central.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        diff > 1e-8,
        "operators should not coincide, max difference {diff}"
    );
    assert_conserved(&f, &cumulant, "cumulant");
    assert_conserved(&f, &central, "central moment");
}

#[test]
fn cumulant_viscosities_round_trip() {
    for (nu, zeta) in [(0.1, 0.05), (1.0 / 6.0, 0.2), (0.3, 0.0)] {
        let operator = CumulantCollision::from_viscosities(nu, zeta);
        assert!((operator.shear_viscosity() - nu).abs() < 1e-14, "nu = {nu}");
        assert!(
            (operator.bulk_viscosity() - zeta).abs() < 1e-14,
            "zeta = {zeta}"
        );
    }
    // The BGK relation: tau-relaxed trace gives zeta = 2 nu / 3.
    let nu = 0.12;
    let bgk = CumulantCollision::from_viscosities(nu, 2.0 * nu / 3.0);
    assert!((bgk.omega - bgk.omega_b).abs() < 1e-14);
}

#[test]
fn central_moments_to_f_cumulant_is_the_exact_inverse() {
    let f = perturbed_populations(3, 0.2);
    let (rho, u) = compute_macroscopic(&f);
    let kappa = compute_central_moments_cumulant(&f, u);
    let back = central_moments_to_f_cumulant(&kappa, rho, u);
    for i in 0..19 {
        assert!((back[i] - f[i]).abs() < 1e-13, "population {i}");
    }
}
