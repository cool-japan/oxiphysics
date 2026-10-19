//! Regression tests for the regularized collision operators.
//!
//! The reconstruction of `f^(1)` from the non-equilibrium stress used the Hermite-free
//! `c_a c_b` polynomial, so `sum_i f^(1)_i = 3/2 tr(Pi)` instead of zero and every regularized
//! operator created or destroyed mass whenever the relaxation rate differed from one. The
//! former mass test ran with `nu = 1/6` (`omega = 1`), which multiplies `f^(1)` by zero and hid
//! the defect. These tests use `omega != 1` and check mass and momentum to round-off for every
//! operator that reconstructs a non-equilibrium part.

use super::test_support::{
    assert_conserved, face_scaled_equilibrium, mass_and_momentum, perturbed_populations,
    shear_wave_amplitude_ratio, shear_wave_decay,
};
use super::*;

/// Second moment `sum_i x_i c_i c_i` as `[xx, yy, zz, xy, xz, yz]`.
fn second_moment(x: &[f64; 19]) -> [f64; 6] {
    compute_pi1_tensor(x, &[0.0; 19])
}

const STRESSES: [[f64; 6]; 4] = [
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [0.3, -0.1, 0.2, 0.05, -0.07, 0.11],
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [-0.4, 0.9, 0.25, -0.3, 0.2, 0.6],
];

const CORRECTION: [f64; 6] = [0.2, 0.1, -0.3, 0.02, 0.04, -0.06];

#[test]
fn reconstructed_f1_carries_no_mass_and_no_momentum() {
    for pi in &STRESSES {
        for (name, f1) in [
            ("pi only", regularized_f1_d3q19(pi, 1.0)),
            ("pi + q", regularized_f1_with_q_d3q19(pi, &CORRECTION)),
        ] {
            for (k, value) in mass_and_momentum(&f1).iter().enumerate() {
                assert!(
                    value.abs() < 1e-14,
                    "{name}, stress {pi:?}: moment {k} = {value}"
                );
            }
        }
    }
}

#[test]
fn pure_trace_stress_is_not_counted_as_mass() {
    // The old reconstruction produced sum f1 = 3/2 tr(Pi) = 4.5 for this tensor.
    let f1 = regularized_f1_d3q19(&[1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 1.0);
    assert!(f1.iter().sum::<f64>().abs() < 1e-14);
    // The trace is nevertheless represented in the second moment.
    let pi_back = second_moment(&f1);
    for (k, &value) in pi_back.iter().enumerate() {
        let expected = if k < 3 { 1.0 } else { 0.0 };
        assert!((value - expected).abs() < 1e-14, "component {k}: {value}");
    }
}

#[test]
fn second_moment_of_f1_is_the_stress_it_was_built_from() {
    for pi in &STRESSES {
        let back = second_moment(&regularized_f1_d3q19(pi, 1.0));
        for (k, (&b, &p)) in back.iter().zip(pi.iter()).enumerate() {
            assert!(
                (b - p).abs() < 1e-14,
                "stress {pi:?}, component {k}: {b} vs {p}"
            );
        }
    }
}

#[test]
fn f1_with_q_is_the_reconstruction_of_the_summed_tensor() {
    for pi in &STRESSES {
        let summed: [f64; 6] = std::array::from_fn(|k| pi[k] + CORRECTION[k]);
        let combined = regularized_f1_with_q_d3q19(pi, &CORRECTION);
        let direct = regularized_f1_d3q19(&summed, 1.0);
        for (a, b) in combined.iter().zip(direct.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
        let back = second_moment(&combined);
        for (k, (&b, &s)) in back.iter().zip(summed.iter()).enumerate() {
            assert!((b - s).abs() < 1e-14, "component {k}");
        }
        // Without a correction both reconstructions coincide.
        let plain = regularized_f1_with_q_d3q19(pi, &[0.0; 6]);
        let base = regularized_f1_d3q19(pi, 1.0);
        assert!(
            plain
                .iter()
                .zip(base.iter())
                .all(|(a, b)| (a - b).abs() < 1e-15)
        );
    }
}

#[test]
fn operators_conserve_mass_and_momentum_with_non_unit_relaxation() {
    let f = face_scaled_equilibrium();
    let (rho, u) = compute_macroscopic(&f);
    let nu = 0.1;
    // nu = 0.1 gives omega = 1.25, so 1 - omega = -0.25 weights the reconstructed part.
    let outputs = [
        (
            "regularized full",
            RegularizedCollisionFull::from_viscosity(nu).collide(&f, rho, u),
        ),
        (
            "recursive regularized",
            RecursiveRegularized::from_viscosity(nu).collide(&f, rho, u),
        ),
        (
            "hybrid recursive regularized (recursive branch)",
            HybridRecursiveRegularized::new(nu, 1.0).collide(&f, rho, u),
        ),
        (
            "hybrid recursive regularized (BGK branch)",
            HybridRecursiveRegularized::new(nu, 1e-9).collide(&f, rho, u),
        ),
        (
            "cumulant",
            CumulantCollision::from_viscosity(nu).collide(&f, rho, u),
        ),
        (
            "raw moment",
            RawMomentCollision::from_viscosity(nu, 0.05).collide(&f, rho, u),
        ),
    ];
    for (name, out) in &outputs {
        assert_conserved(&f, out, 1e-14, name);
        let max_change = f
            .iter()
            .zip(out.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_change > 1e-5,
            "{name}: the collision must change the state"
        );
    }
}

#[test]
fn regularized_operators_conserve_for_random_non_equilibrium_states() {
    for nu in [0.4, 0.1, 0.01] {
        let regularized = RegularizedCollisionFull::from_viscosity(nu);
        let recursive = RecursiveRegularized::from_viscosity(nu);
        for seed in 1..=10 {
            let f = perturbed_populations(seed, 0.2);
            let (rho, u) = compute_macroscopic(&f);
            assert_conserved(
                &f,
                &regularized.collide(&f, rho, u),
                1e-14,
                &format!("regularized, nu = {nu}, seed {seed}"),
            );
            assert_conserved(
                &f,
                &recursive.collide(&f, rho, u),
                1e-14,
                &format!("recursive, nu = {nu}, seed {seed}"),
            );
        }
    }
}

#[test]
fn regularized_collision_relaxes_the_stress_with_the_viscosity() {
    let nu = 0.1;
    let omega = 1.0 / (3.0 * nu + 0.5);
    let coll = RegularizedCollisionFull::from_viscosity(nu);
    for seed in 1..=6 {
        let f = perturbed_populations(seed, 0.05);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let before = compute_pi1_tensor(&f, &feq);
        let after = compute_pi1_tensor(&coll.collide(&f, rho, u), &feq);
        for (k, (&a, &b)) in after.iter().zip(before.iter()).enumerate() {
            assert!(
                (a - (1.0 - omega) * b).abs() < 1e-14,
                "seed {seed}, stress component {k}: {a} vs {}",
                (1.0 - omega) * b
            );
        }
    }
}

#[test]
fn regularized_collision_keeps_the_equilibrium_fixed() {
    let (rho, u) = (1.2, [0.03, -0.04, 0.02]);
    let feq = compute_equilibrium_d3q19(rho, u);
    for out in [
        RegularizedCollisionFull::from_viscosity(0.1).collide(&feq, rho, u),
        RecursiveRegularized::from_viscosity(0.1).collide(&feq, rho, u),
        HybridRecursiveRegularized::new(0.1, 0.05).collide(&feq, rho, u),
    ] {
        for (a, b) in out.iter().zip(feq.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }
}

#[test]
fn shear_wave_decays_with_the_requested_viscosity() {
    let (nu, ny, steps, amplitude) = (0.02, 32, 400, 0.02);
    let coll = RegularizedCollisionFull::from_viscosity(nu);
    let ratio =
        shear_wave_amplitude_ratio(|f, rho, u| coll.collide(f, rho, u), ny, amplitude, steps);
    let expected = shear_wave_decay(nu, ny, steps);
    assert!(
        (ratio - expected).abs() < 0.01,
        "decay {ratio} vs exp(-nu k^2 t) = {expected}"
    );
}
