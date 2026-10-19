// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the private stability analysis and validation helpers of
//! [`crate::dynamics`].

use super::*;

/// Helper: build a 1-DOF system with K = k, M = m.
fn single_dof(k: f64, m: f64) -> DynamicsSystem {
    DynamicsSystem::new(1, vec![m], vec![vec![k]])
}
// ── Stability-analysis helpers ─────────────────────────────────────────

/// Companion matrix of the monic cubic `λ³ + c2·λ² + c1·λ + c0`.
fn companion(c2: f64, c1: f64, c0: f64) -> [[f64; 3]; 3] {
    [[0.0, 0.0, -c0], [1.0, 0.0, -c1], [0.0, 1.0, -c2]]
}

/// Coefficients `(c2, c1, c0)` of `(λ − r0)(λ² − 2·m·cos θ·λ + m²)` for a real root `r0` and
/// a complex pair of modulus `m` and argument `θ`.
fn cubic_with_complex_pair(r0: f64, m: f64, theta: f64) -> (f64, f64, f64) {
    let s = 2.0 * m * theta.cos();
    let q = m * m;
    (-(r0 + s), q + r0 * s, -r0 * q)
}

#[test]
fn test_schur_cohn_test_matches_the_roots_of_a_cubic() {
    // Real roots: stable iff all |r| ≤ 1 (roots on the circle count as stable).
    let real_roots: [([f64; 3], bool); 6] = [
        ([0.5, -0.3, 0.9], true),
        ([0.0, 0.0, 0.0], true),
        ([1.0, 0.5, -0.5], true),
        ([-1.0, 0.2, 0.3], true),
        ([1.1, 0.2, 0.3], false),
        ([-1.05, 0.0, 0.9], false),
    ];
    for (r, stable) in real_roots {
        let c2 = -(r[0] + r[1] + r[2]);
        let c1 = r[0] * r[1] + r[0] * r[2] + r[1] * r[2];
        let c0 = -(r[0] * r[1] * r[2]);
        assert_eq!(
            schur_cohn_cubic_stable(&companion(c2, c1, c0)),
            stable,
            "roots {r:?}"
        );
    }
    // Complex pair with a real root.
    for (r0, m, theta, stable) in [
        (0.5, 0.9, 0.7, true),
        (-0.4, 0.99, 2.5, true),
        (0.5, 1.0, 1.3, true),
        (0.5, 1.1, 0.7, false),
        (0.95, 1.02, 3.0, false),
        (-1.2, 0.5, 1.0, false),
    ] {
        let (c2, c1, c0) = cubic_with_complex_pair(r0, m, theta);
        assert_eq!(
            schur_cohn_cubic_stable(&companion(c2, c1, c0)),
            stable,
            "real root {r0}, pair modulus {m}, argument {theta}"
        );
    }
    let nan = [[f64::NAN; 3]; 3];
    assert!(!schur_cohn_cubic_stable(&nan));
}

#[test]
fn test_scanned_critical_step_matches_the_closed_form_of_classical_newmark() {
    for (beta, gamma) in [
        (0.0, 0.5),
        (1.0 / 6.0, 0.5),
        (0.2, 0.5),
        (0.1, 0.7),
        (0.0, 0.6),
        (0.05, 0.6),
        (0.2, 0.6),
    ] {
        let integrator = NewmarkIntegrator {
            beta,
            gamma,
            hht_alpha: 0.0,
        };
        let closed = integrator.critical_omega_dt();
        let scanned = integrator.scan_critical_omega_dt();
        assert!(closed.is_finite());
        assert!(
            (closed - scanned).abs() < 1e-6 * closed,
            "β = {beta}, γ = {gamma}: closed form {closed}, scan {scanned}"
        );
    }
    // Unconditionally stable and unstable classical parameter sets.
    let stable = NewmarkIntegrator::new_average_acceleration();
    assert_eq!(stable.scan_critical_omega_dt(), f64::INFINITY);
    let unstable = NewmarkIntegrator {
        beta: 0.25,
        gamma: 0.4,
        hht_alpha: 0.0,
    };
    assert_eq!(unstable.scan_critical_omega_dt(), 0.0);
    assert_eq!(unstable.critical_omega_dt(), 0.0);
}

#[test]
fn test_hht_amplification_matches_the_asymptotic_spectral_radius() {
    // For the optimal HHT-α parameters the spectral radius at ω·Δt → ∞ is (1+α)/(1−α):
    // the amplification matrix has two eigenvalues of that modulus.  With one real
    // eigenvalue the determinant identity |det| = ρ²·|λ3| gives a quick check through the
    // Schur–Cohn test: scaling the matrix by 1/ρ∞ keeps it stable, by 1.01/ρ∞ does not.
    let alpha = -0.2;
    let hht = NewmarkIntegrator::new_hht_alpha(alpha).expect("alpha in range");
    let rho_inf = (1.0 + alpha) / (1.0 - alpha);
    let a = hht.modal_amplification(1.0e6);
    let scaled = |factor: f64| {
        let mut m = a;
        for row in m.iter_mut() {
            for x in row.iter_mut() {
                *x *= factor;
            }
        }
        m
    };
    assert!(schur_cohn_cubic_stable(&scaled(0.999 / rho_inf)));
    assert!(!schur_cohn_cubic_stable(&scaled(1.01 / rho_inf)));
}

// ── System validation and helpers ──────────────────────────────────────

#[test]
fn test_mass_stiffness_combination() {
    let sys = DynamicsSystem::new(2, vec![2.0, 3.0], vec![vec![4.0, -1.0], vec![-1.0, 5.0]]);
    let c = sys.mass_stiffness_combination(0.5, 2.0);
    assert!((c[0][0] - (1.0 + 8.0)).abs() < 1e-14);
    assert!((c[0][1] + 2.0).abs() < 1e-14);
    assert!((c[1][0] + 2.0).abs() < 1e-14);
    assert!((c[1][1] - (1.5 + 10.0)).abs() < 1e-14);
}

#[test]
fn test_system_validation_rejects_inconsistent_data() {
    let mut sys = single_dof(1.0, 1.0);
    assert!(sys.validate().is_ok());
    sys.mass_diag = vec![1.0, 2.0];
    assert!(matches!(sys.validate(), Err(Error::InvalidInput(_))));
    let mut sys = single_dof(1.0, 1.0);
    sys.stiffness = vec![vec![f64::NAN]];
    assert!(matches!(sys.validate(), Err(Error::InvalidInput(_))));
    let mut sys = single_dof(1.0, 1.0);
    sys.damping_beta = f64::INFINITY;
    assert!(matches!(sys.validate(), Err(Error::InvalidInput(_))));
    let mut sys = single_dof(1.0, 1.0);
    sys.stiffness = vec![vec![1.0, 2.0]];
    assert!(matches!(sys.validate(), Err(Error::InvalidInput(_))));
}

#[test]
fn test_integrator_parameter_validation() {
    let sys = single_dof(1.0, 1.0);
    for bad in [
        NewmarkIntegrator {
            beta: -0.1,
            gamma: 0.5,
            hht_alpha: 0.0,
        },
        NewmarkIntegrator {
            beta: 0.25,
            gamma: f64::NAN,
            hht_alpha: 0.0,
        },
        NewmarkIntegrator {
            beta: 0.25,
            gamma: 0.5,
            hht_alpha: -1.0,
        },
    ] {
        assert!(matches!(
            bad.stepper(&sys, 0.1),
            Err(Error::InvalidInput(_))
        ));
    }
}
