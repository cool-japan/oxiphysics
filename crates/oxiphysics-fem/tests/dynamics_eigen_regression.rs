// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the structural-dynamics integrators in
//! [`oxiphysics_fem::dynamics`] and the dense subspace eigen-solvers in
//! [`oxiphysics_fem::eigen`].
//!
//! What used to be wrong (and is pinned down here):
//!
//! * `NewmarkIntegrator::step` and `BatheIntegrator::step` were *explicit* predictor
//!   corrector updates dressed up as implicit schemes: they diverged for `ω·Δt` beyond
//!   about 2 although `critical_dt` advertised unconditional stability.
//! * `NewmarkIntegrator::new_hht_alpha` used `γ = 1/2 + α < 1/2` (negative numerical damping)
//!   and had no `(1 + α)` force weighting.
//! * The Bathe step used an invented `0.9/0.1` blend and overwrote the velocity with a
//!   first-order difference.
//! * `lump_mass_hrz` divided the total mass by `n_dim` although the matrix already carries the
//!   mass once per direction, and `estimate_omega_max` used `K_ii/M_ii`, which underestimates
//!   `ω_max` of every coupled system, so `safe_explicit_dt` was not safe.
//! * `subspace_iteration` iterated with `A` (largest eigenpairs) without ever
//!   re-orthonormalising, and `shift_invert_power` returned `num_modes` copies of the single
//!   eigenpair nearest to the shift.

use oxiphysics_fem::Error;
use oxiphysics_fem::dynamic::element_mass_matrix;
use oxiphysics_fem::dynamics::{
    BatheIntegrator, CentralDifference, DynamicsSystem, NewmarkIntegrator, estimate_omega_max,
    lump_mass_hrz, safe_explicit_dt,
};
use oxiphysics_fem::eigen::{
    eigenpair_residuals, jacobi_eigen_dense, shift_invert_power, subspace_iteration,
};
use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Undamped single-DOF oscillator with stiffness `k` and mass `m`.
fn sdof(k: f64, m: f64) -> DynamicsSystem {
    DynamicsSystem::new(1, vec![m], vec![vec![k]])
}

/// Three coupled, damped DOFs with unequal masses.
fn coupled_damped_system() -> DynamicsSystem {
    let mut sys = DynamicsSystem::new(
        3,
        vec![2.0, 1.0, 0.5],
        vec![
            vec![300.0, -150.0, 0.0],
            vec![-150.0, 250.0, -100.0],
            vec![0.0, -100.0, 100.0],
        ],
    );
    sys.damping_alpha = 0.7;
    sys.damping_beta = 0.002;
    sys
}

/// Time-varying load of [`coupled_damped_system`].
fn coupled_load(t: f64) -> Vec<f64> {
    vec![10.0 * (3.0 * t).sin(), 0.0, 5.0 * (7.0 * t).cos()]
}

/// Largest absolute component of `M·a + C·v + K·u − f`.
fn motion_residual(sys: &DynamicsSystem, u: &[f64], v: &[f64], a: &[f64], f: &[f64]) -> f64 {
    let k_u = sys.internal_force(u);
    let c_v = sys.rayleigh_damping_force(v);
    (0..sys.n_dof)
        .map(|i| (sys.mass_diag[i] * a[i] + c_v[i] + k_u[i] - f[i]).abs())
        .fold(0.0, f64::max)
}

/// `½(u² + v²)` of a unit-frequency, unit-mass oscillator.
fn sdof_energy(u: f64, v: f64) -> f64 {
    0.5 * (u * u + v * v)
}

/// Euclidean inner product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Dense symmetric matrix `Q·diag(eigenvalues)·Qᵀ` with `Q` a product of Householder
/// reflectors; returns the matrix and the columns of `Q` (the eigenvectors).
fn matrix_with_spectrum(eigenvalues: &[f64]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = eigenvalues.len();
    // Orthogonal Q = H1·H2·H3 with H = I − 2·w·wᵀ/(wᵀw).
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for seed in 1..=3_usize {
        let w: Vec<f64> = (0..n)
            .map(|i| ((i + 1) as f64 * seed as f64 * 0.731).sin() + 0.3 * seed as f64)
            .collect();
        let ww = dot(&w, &w);
        // q ← H·q, one column at a time: column_j ← column_j − 2·w·(wᵀ·column_j)/(wᵀw).
        let projections: Vec<f64> = (0..n)
            .map(|col| {
                let column: Vec<f64> = q.iter().map(|row| row[col]).collect();
                2.0 * dot(&w, &column) / ww
            })
            .collect();
        for (row, w_row) in q.iter_mut().zip(w.iter()) {
            for (entry, proj) in row.iter_mut().zip(projections.iter()) {
                *entry -= proj * w_row;
            }
        }
    }
    let columns: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| q[i][j]).collect()).collect();
    let a: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    columns
                        .iter()
                        .zip(eigenvalues.iter())
                        .map(|(c, lambda)| lambda * c[i] * c[j])
                        .sum()
                })
                .collect()
        })
        .collect();
    (a, columns)
}

/// Assert that the vectors are orthonormal.
fn assert_orthonormal(vectors: &[Vec<f64>], tol: f64) {
    for (i, vi) in vectors.iter().enumerate() {
        for (j, vj) in vectors.iter().enumerate() {
            let expected = if i == j { 1.0 } else { 0.0 };
            let d = dot(vi, vj);
            assert!(
                (d - expected).abs() < tol,
                "<v{i}, v{j}> = {d}, expected {expected}"
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Newmark-β: the step is genuinely implicit (F1141, F1142)
// ─────────────────────────────────────────────────────────────────────────────

/// Average acceleration is unconditionally stable and conserves the quadratic energy: even at
/// `ω·Δt = 4π` (two periods per step), 10 and 100 the oscillator neither grows nor loses
/// energy.  The old "predictor-corrector" step reached 3·10¹¹⁸ for `ω·Δt ≳ 2.5`.
#[test]
fn average_acceleration_is_unconditionally_stable_and_conserves_energy() {
    let sys = sdof(1.0, 1.0);
    let integrator = NewmarkIntegrator::new_average_acceleration();
    for (omega_dt, steps) in [(4.0 * PI, 100_usize), (10.0, 1000), (100.0, 1000)] {
        let (mut u, mut v, mut a) = (vec![1.0], vec![0.0], vec![-1.0]);
        for step in 1..=steps {
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &[0.0], omega_dt)
                .expect("implicit step");
            assert!(
                u[0].abs() <= 1.0 + 1e-9,
                "ω·Δt = {omega_dt}, step {step}: |u| = {}",
                u[0].abs()
            );
            let drift = (sdof_energy(u[0], v[0]) - 0.5).abs();
            assert!(
                drift < 1e-9,
                "ω·Δt = {omega_dt}, step {step}: energy drift {drift}"
            );
        }
    }
}

/// The exact discrete modal recurrence of the trapezoidal rule, for one mode `(ω, c, ċ)`.
fn trapezoid_modal_step(omega: f64, dt: f64, c: f64, dc: f64) -> (f64, f64) {
    let q = 0.25 * (omega * dt).powi(2);
    let d = 1.0 + q;
    (
        ((1.0 - q) * c + dt * dc) / d,
        ((1.0 - q) * dc - omega * omega * dt * c) / d,
    )
}

/// A coupled two-DOF system must follow the trapezoidal rule applied to its *modes*.  The old
/// decoupled predictor treated each DOF independently and therefore missed the stiffness
/// coupling in the effective matrix.
#[test]
fn average_acceleration_reproduces_the_modal_recurrence_of_a_coupled_system() {
    // K = [[2, -1], [-1, 1]], M = I: ω² = (3 ∓ √5)/2 with mode shapes [1, 2 − ω²].
    let sys = DynamicsSystem::new(2, vec![1.0, 1.0], vec![vec![2.0, -1.0], vec![-1.0, 1.0]]);
    let lambdas = [(3.0 - 5.0_f64.sqrt()) / 2.0, (3.0 + 5.0_f64.sqrt()) / 2.0];
    let shapes: Vec<[f64; 2]> = lambdas.iter().map(|l| [1.0, 2.0 - l]).collect();
    assert!(
        dot(&shapes[0], &shapes[1]).abs() < 1e-14,
        "modes must be orthogonal"
    );

    let integrator = NewmarkIntegrator::new_average_acceleration();
    for dt in [0.2, 1.0, 2.5] {
        let u0 = vec![1.0, 0.3];
        let v0 = vec![-0.2, 0.5];
        let mut modal: Vec<(f64, f64)> = shapes
            .iter()
            .map(|s| {
                let norm = dot(s, s);
                (dot(s, &u0) / norm, dot(s, &v0) / norm)
            })
            .collect();
        let (mut u, mut v) = (u0.clone(), v0.clone());
        let mut a = sys
            .initial_acceleration(&u, &v, &[0.0, 0.0])
            .expect("consistent initial acceleration");
        for step in 1..=40 {
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &[0.0, 0.0], dt)
                .expect("implicit step");
            for (state, lambda) in modal.iter_mut().zip(lambdas.iter()) {
                *state = trapezoid_modal_step(lambda.sqrt(), dt, state.0, state.1);
            }
            for i in 0..2 {
                let u_ref: f64 = modal
                    .iter()
                    .zip(shapes.iter())
                    .map(|(m, s)| m.0 * s[i])
                    .sum();
                let v_ref: f64 = modal
                    .iter()
                    .zip(shapes.iter())
                    .map(|(m, s)| m.1 * s[i])
                    .sum();
                assert!(
                    (u[i] - u_ref).abs() < 1e-11 && (v[i] - v_ref).abs() < 1e-11,
                    "dt = {dt}, step {step}, dof {i}: u = {}, expected {u_ref}; v = {}, expected {v_ref}",
                    u[i],
                    v[i]
                );
            }
        }
    }
}

/// After every step the discrete equation `M·a + C·v + K·u = f` must hold exactly at the new
/// time level for a coupled, damped, loaded system and any Newmark parameters.
#[test]
fn newmark_satisfies_the_equation_of_motion_at_the_end_of_each_step() {
    let sys = coupled_damped_system();
    let dt = 0.02;
    for integrator in [
        NewmarkIntegrator::new_average_acceleration(),
        NewmarkIntegrator::new_linear_acceleration(),
        NewmarkIntegrator {
            beta: 0.3,
            gamma: 0.6,
            hht_alpha: 0.0,
        },
    ] {
        let mut u = vec![0.01, -0.02, 0.03];
        let mut v = vec![0.1, 0.0, -0.1];
        let mut a = sys
            .initial_acceleration(&u, &v, &coupled_load(0.0))
            .expect("consistent initial acceleration");
        for step in 1..=200 {
            let f = coupled_load(step as f64 * dt);
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &f, dt)
                .expect("implicit step");
            let residual = motion_residual(&sys, &u, &v, &a, &f);
            assert!(
                residual < 1e-9,
                "{integrator:?}, step {step}: equilibrium residual {residual}"
            );
        }
    }
}

/// Damped free vibration follows the analytic response with second-order accuracy.
#[test]
fn newmark_with_rayleigh_damping_converges_to_the_analytic_decay() {
    let mut sys = sdof(1.0, 1.0);
    sys.damping_alpha = 0.2;
    sys.damping_beta = 0.01;
    let zeta = sys.damping_alpha / 2.0 + sys.damping_beta / 2.0; // ω = 1
    let omega_d = (1.0 - zeta * zeta).sqrt();
    let exact =
        |t: f64| (-zeta * t).exp() * ((omega_d * t).cos() + zeta / omega_d * (omega_d * t).sin());
    let t_end = 4.0 * PI;
    let integrator = NewmarkIntegrator::new_average_acceleration();
    let errors: Vec<f64> = [100_usize, 200, 400]
        .iter()
        .map(|&steps| {
            let dt = t_end / steps as f64;
            let (mut u, mut v) = (vec![1.0], vec![0.0]);
            let mut a = sys.initial_acceleration(&u, &v, &[0.0]).expect("a0");
            for _ in 0..steps {
                integrator
                    .step(&sys, &mut u, &mut v, &mut a, &[0.0], dt)
                    .expect("implicit step");
            }
            (u[0] - exact(t_end)).abs()
        })
        .collect();
    assert!(errors[2] < 1e-4, "errors {errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.6..4.4).contains(&ratio),
            "halving dt should cut the error by about 4, got {ratio} ({errors:?})"
        );
    }
}

/// Reusing a prepared stepper gives the same result as the one-shot step.
#[test]
fn stepper_matches_the_one_shot_step() {
    let sys = coupled_damped_system();
    let integrator = NewmarkIntegrator::new_average_acceleration();
    let dt = 0.01;
    let stepper = integrator.stepper(&sys, dt).expect("stepper");
    let (mut u1, mut v1) = (vec![0.01, 0.0, 0.0], vec![0.0; 3]);
    let mut a1 = sys
        .initial_acceleration(&u1, &v1, &coupled_load(0.0))
        .expect("a0");
    let (mut u2, mut v2, mut a2) = (u1.clone(), v1.clone(), a1.clone());
    for step in 1..=50 {
        let f = coupled_load(step as f64 * dt);
        stepper
            .step(&mut u1, &mut v1, &mut a1, &f)
            .expect("stepper step");
        integrator
            .step(&sys, &mut u2, &mut v2, &mut a2, &f, dt)
            .expect("one-shot step");
    }
    for i in 0..3 {
        assert!((u1[i] - u2[i]).abs() < 1e-14 && (v1[i] - v2[i]).abs() < 1e-14);
    }
}

/// Invalid input and singular systems are reported, and a failed step leaves the state alone.
#[test]
fn newmark_reports_invalid_input_and_singular_systems() {
    let integrator = NewmarkIntegrator::new_average_acceleration();
    let sys = sdof(1.0, 1.0);
    let (mut u, mut v, mut a) = (vec![1.0], vec![0.5], vec![-1.0]);
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        let result = integrator.step(&sys, &mut u, &mut v, &mut a, &[0.0], dt);
        assert!(
            matches!(result, Err(Error::InvalidInput(_))),
            "dt = {dt}: {result:?}"
        );
    }
    let wrong = integrator.step(&sys, &mut u, &mut v, &mut a, &[0.0, 0.0], 0.1);
    assert!(matches!(wrong, Err(Error::InvalidInput(_))), "{wrong:?}");
    let mut short = vec![1.0, 2.0];
    let wrong = integrator.step(&sys, &mut short, &mut v, &mut a, &[0.0], 0.1);
    assert!(matches!(wrong, Err(Error::InvalidInput(_))), "{wrong:?}");

    let negative_mass = DynamicsSystem::new(1, vec![-1.0], vec![vec![1.0]]);
    let result = integrator.step(&negative_mass, &mut u, &mut v, &mut a, &[0.0], 0.1);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");

    // A massless DOF without stiffness leaves the effective matrix singular.
    let singular = DynamicsSystem::new(2, vec![0.0, 1.0], vec![vec![0.0, 0.0], vec![0.0, 1.0]]);
    let (mut u2, mut v2, mut a2) = (vec![0.0; 2], vec![0.0; 2], vec![0.0; 2]);
    let result = integrator.step(&singular, &mut u2, &mut v2, &mut a2, &[0.0; 2], 0.1);
    assert!(matches!(result, Err(Error::SingularSystem)), "{result:?}");

    assert_eq!(
        (u[0], v[0], a[0]),
        (1.0, 0.5, -1.0),
        "failed steps must not touch the state"
    );

    // A non-finite load is caught instead of silently poisoning the state.
    let result = integrator.step(&sys, &mut u, &mut v, &mut a, &[f64::NAN], 0.1);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");
    assert_eq!((u[0], v[0], a[0]), (1.0, 0.5, -1.0));
}

/// `initial_acceleration` returns `M⁻¹(f − C·v − K·u)` and refuses massless DOFs.
#[test]
fn initial_acceleration_balances_the_equation_of_motion() {
    let sys = coupled_damped_system();
    let (u, v, f) = (
        vec![0.01, -0.02, 0.03],
        vec![0.1, 0.0, -0.1],
        coupled_load(0.3),
    );
    let a = sys.initial_acceleration(&u, &v, &f).expect("a0");
    assert!(motion_residual(&sys, &u, &v, &a, &f) < 1e-12);
    let massless = DynamicsSystem::new(1, vec![0.0], vec![vec![1.0]]);
    assert!(matches!(
        massless.initial_acceleration(&[0.0], &[0.0], &[0.0]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        sys.initial_acceleration(&u[..2], &v, &f),
        Err(Error::InvalidInput(_))
    ));
}

// ─────────────────────────────────────────────────────────────────────────────
// critical_dt (F1141)
// ─────────────────────────────────────────────────────────────────────────────

/// Largest `|u|` of a unit oscillator integrated for `steps` steps of `dt`.
fn peak_amplitude(integrator: &NewmarkIntegrator, dt: f64, steps: usize) -> f64 {
    let sys = sdof(1.0, 1.0);
    let (mut u, mut v, mut a) = (vec![1.0], vec![0.0], vec![-1.0]);
    let mut peak = 1.0_f64;
    for _ in 0..steps {
        integrator
            .step(&sys, &mut u, &mut v, &mut a, &[0.0], dt)
            .expect("implicit step");
        peak = peak.max(u[0].abs());
        if !peak.is_finite() || peak > 1.0e100 {
            break;
        }
    }
    peak
}

/// `critical_dt` is infinite only for `2β ≥ γ ≥ 1/2`; otherwise it equals
/// `1/(ω·√(γ/2 − β))`, and a time step just inside/outside it is stable/unstable.
#[test]
fn critical_dt_follows_newmark_stability_theory() {
    let omega = 2.5;
    assert_eq!(
        NewmarkIntegrator::new_average_acceleration().critical_dt(omega),
        f64::INFINITY
    );
    for alpha in [-1.0 / 3.0, -0.2, -0.05, 0.0] {
        let hht = NewmarkIntegrator::new_hht_alpha(alpha).expect("alpha in range");
        assert_eq!(hht.critical_dt(omega), f64::INFINITY, "alpha = {alpha}");
    }
    let central = NewmarkIntegrator {
        beta: 0.0,
        gamma: 0.5,
        hht_alpha: 0.0,
    };
    assert!((central.critical_dt(omega) - 2.0 / omega).abs() < 1e-12);
    let linear = NewmarkIntegrator::new_linear_acceleration();
    assert!((linear.critical_dt(omega) - 2.0 * 3.0_f64.sqrt() / omega).abs() < 1e-12);
    let mixed = NewmarkIntegrator {
        beta: 0.1,
        gamma: 0.7,
        hht_alpha: 0.0,
    };
    assert!((mixed.critical_dt(1.0) - 2.0).abs() < 1e-12);
    // γ < 1/2 is unstable for every step size, however 2β ≥ γ looks.
    let anti_damped = NewmarkIntegrator {
        beta: 0.25,
        gamma: 0.4,
        hht_alpha: 0.0,
    };
    assert_eq!(anti_damped.critical_dt(1.0), 0.0);
    assert_eq!(linear.critical_dt(0.0), f64::INFINITY);
    assert!(linear.critical_dt(-1.0).is_nan());
    assert!(linear.critical_dt(f64::NAN).is_nan());

    for (integrator, label) in [(central, "central"), (linear, "linear"), (mixed, "mixed")] {
        let dt_crit = integrator.critical_dt(1.0);
        let inside = peak_amplitude(&integrator, 0.95 * dt_crit, 3000);
        let outside = peak_amplitude(&integrator, 1.05 * dt_crit, 3000);
        assert!(
            inside <= 1.0 + 1e-9,
            "{label}: inside the limit |u| = {inside}"
        );
        assert!(
            outside > 1.0e6,
            "{label}: outside the limit |u| = {outside}"
        );
    }
}

/// For HHT-α with parameters other than the optimal ones the plain Newmark criterion does not
/// apply (e.g. `γ = 1/2` is unstable once α < 0); the limit follows from the amplification
/// matrix, and unconditionally stable combinations still report no limit.
#[test]
fn critical_dt_of_non_canonical_hht_parameters_uses_the_amplification_matrix() {
    // 2β ≥ γ ≥ 1/2 holds, but with α = −0.1 the scheme grows for every step.
    let unstable = NewmarkIntegrator {
        beta: 0.25,
        gamma: 0.5,
        hht_alpha: -0.1,
    };
    assert!(
        unstable.critical_dt(1.0) < 0.05,
        "{}",
        unstable.critical_dt(1.0)
    );
    assert!(peak_amplitude(&unstable, 0.5, 200) > 1.0);
    // γ = 1/2 − α and β = γ/2 is still unconditionally stable for HHT.
    let stable = NewmarkIntegrator {
        beta: 0.3,
        gamma: 0.6,
        hht_alpha: -0.1,
    };
    assert_eq!(stable.critical_dt(1.0), f64::INFINITY);
    assert!(peak_amplitude(&stable, 50.0, 500) <= 1.0 + 1e-9);
}

// ─────────────────────────────────────────────────────────────────────────────
// HHT-α (F1140)
// ─────────────────────────────────────────────────────────────────────────────

/// Hilber–Hughes–Taylor parameters: `γ = 1/2 − α ≥ 1/2` and `β = (1 − α)²/4`.  The old code
/// used `γ = 1/2 + α < 1/2`, which amplifies energy.
#[test]
fn hht_alpha_uses_the_hilber_hughes_taylor_parameters() {
    let hht = NewmarkIntegrator::new_hht_alpha(-0.1).expect("alpha in range");
    assert!((hht.gamma - 0.6).abs() < 1e-14, "gamma = {}", hht.gamma);
    assert!((hht.beta - 0.3025).abs() < 1e-14, "beta = {}", hht.beta);
    assert!(hht.gamma >= 0.5 && 2.0 * hht.beta >= hht.gamma);
    for alpha in [-1.0 / 3.0, -0.25, -0.1] {
        let hht = NewmarkIntegrator::new_hht_alpha(alpha).expect("alpha in range");
        assert!((hht.gamma - (0.5 - alpha)).abs() < 1e-14);
        assert!((hht.beta - (1.0 - alpha).powi(2) / 4.0).abs() < 1e-14);
    }
    assert!(matches!(
        NewmarkIntegrator::new_hht_alpha(0.1),
        Err(Error::InvalidInput(_))
    ));
}

/// HHT-α dissipates when the step is far too large for the mode (`ω·Δt = π`, α = −0.3) and
/// never grows.  The old sign-flipped parameters amplified the response.
#[test]
fn hht_alpha_dissipates_poorly_resolved_modes_without_growth() {
    let sys = sdof(1.0, 1.0);
    let hht = NewmarkIntegrator::new_hht_alpha(-0.3).expect("alpha in range");
    let (mut u, mut v) = (vec![1.0], vec![0.0]);
    let mut a = sys.initial_acceleration(&u, &v, &[0.0]).expect("a0");
    let mut peak = 1.0_f64;
    for _ in 0..100 {
        hht.step(&sys, &mut u, &mut v, &mut a, &[0.0], PI)
            .expect("implicit step");
        peak = peak.max(u[0].abs());
    }
    assert!(peak <= 1.0 + 1e-9, "peak |u| = {peak}");
    assert!(
        u[0].abs() < 1e-5 && v[0].abs() < 1e-5,
        "the oscillation should have decayed, got u = {}, v = {}",
        u[0],
        v[0]
    );
}

/// The `(1 + α)` weighting is what keeps HHT accurate at low frequencies.  Plain Newmark with
/// the very same `β`, `γ` has the same asymptotic dissipation (its spectral radius at large
/// `ω·Δt` is also `(1+α)/(1−α)`) but loses a quarter of the energy per period at 40 steps per
/// period, whereas HHT loses less than 0.1 %.
#[test]
fn hht_alpha_keeps_low_frequency_energy_unlike_plain_newmark_with_the_same_beta_gamma() {
    let sys = sdof(1.0, 1.0);
    let hht = NewmarkIntegrator::new_hht_alpha(-0.3).expect("alpha in range");
    let plain = NewmarkIntegrator {
        hht_alpha: 0.0,
        ..hht
    };
    let dt = 2.0 * PI / 40.0;
    let energy_after_one_period = |integrator: &NewmarkIntegrator| {
        let (mut u, mut v) = (vec![1.0], vec![0.0]);
        let mut a = sys.initial_acceleration(&u, &v, &[0.0]).expect("a0");
        for _ in 0..40 {
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &[0.0], dt)
                .expect("implicit step");
        }
        sdof_energy(u[0], v[0]) / 0.5
    };
    let hht_ratio = energy_after_one_period(&hht);
    let plain_ratio = energy_after_one_period(&plain);
    assert!(hht_ratio > 0.998, "HHT energy ratio {hht_ratio}");
    assert!(
        plain_ratio < 0.8,
        "plain Newmark energy ratio {plain_ratio}"
    );
}

/// HHT with a time-varying load needs the load at both ends of the step; with it the forced
/// response of a SDOF oscillator converges with second order, for α = 0 and α < 0.
#[test]
fn hht_alpha_forced_response_is_second_order_with_both_loads() {
    let (omega_f, amplitude) = (0.7, 1.0);
    let load = |t: f64| amplitude * (omega_f * t).sin();
    let exact =
        |t: f64| amplitude / (1.0 - omega_f * omega_f) * ((omega_f * t).sin() - omega_f * t.sin());
    let sys = sdof(1.0, 1.0);
    let t_end = 20.0;
    for alpha in [0.0, -0.2] {
        let integrator = NewmarkIntegrator::new_hht_alpha(alpha).expect("alpha in range");
        let errors: Vec<f64> = [200_usize, 400, 800]
            .iter()
            .map(|&steps| {
                let dt = t_end / steps as f64;
                let stepper = integrator.stepper(&sys, dt).expect("stepper");
                let (mut u, mut v) = (vec![0.0], vec![0.0]);
                let mut a = sys.initial_acceleration(&u, &v, &[load(0.0)]).expect("a0");
                for n in 0..steps {
                    let (t0, t1) = (n as f64 * dt, (n + 1) as f64 * dt);
                    stepper
                        .step_with_loads(&mut u, &mut v, &mut a, &[load(t0)], &[load(t1)])
                        .expect("stepper step");
                }
                (u[0] - exact(t_end)).abs()
            })
            .collect();
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.8..4.4).contains(&ratio),
                "alpha = {alpha}: error ratio {ratio} ({errors:?})"
            );
        }
    }
}

/// With α ≠ 0 the end-of-step residual is the α-weighted equation
/// `M·a₊ + (1+α)(C·v₊ + K·u₊) − α(C·v + K·u) = (1+α)·f₊ − α·f`.
#[test]
fn hht_alpha_satisfies_the_alpha_weighted_equation_of_motion() {
    let sys = coupled_damped_system();
    let alpha = -0.2;
    let integrator = NewmarkIntegrator::new_hht_alpha(alpha).expect("alpha in range");
    let dt = 0.02;
    let stepper = integrator.stepper(&sys, dt).expect("stepper");
    let mut u = vec![0.01, -0.02, 0.03];
    let mut v = vec![0.1, 0.0, -0.1];
    let mut a = sys
        .initial_acceleration(&u, &v, &coupled_load(0.0))
        .expect("a0");
    for step in 0..100 {
        let f_prev = coupled_load(step as f64 * dt);
        let f_next = coupled_load((step + 1) as f64 * dt);
        let (u_old, v_old) = (u.clone(), v.clone());
        stepper
            .step_with_loads(&mut u, &mut v, &mut a, &f_prev, &f_next)
            .expect("stepper step");
        let (ku_new, cv_new) = (sys.internal_force(&u), sys.rayleigh_damping_force(&v));
        let (ku_old, cv_old) = (
            sys.internal_force(&u_old),
            sys.rayleigh_damping_force(&v_old),
        );
        for i in 0..3 {
            let lhs = sys.mass_diag[i] * a[i] + (1.0 + alpha) * (cv_new[i] + ku_new[i])
                - alpha * (cv_old[i] + ku_old[i]);
            let rhs = (1.0 + alpha) * f_next[i] - alpha * f_prev[i];
            assert!(
                (lhs - rhs).abs() < 1e-9,
                "step {step}, dof {i}: {lhs} vs {rhs}"
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bathe composite method (F1143)
// ─────────────────────────────────────────────────────────────────────────────

/// Mean period of a unit oscillator from the downward zero crossings of `u`.
fn measured_period(trace: &[f64], dt: f64) -> f64 {
    let crossings: Vec<f64> = trace
        .windows(2)
        .enumerate()
        .filter(|(_, w)| w[0] > 0.0 && w[1] <= 0.0)
        .map(|(i, w)| (i as f64 + w[0] / (w[0] - w[1])) * dt)
        .collect();
    assert!(
        crossings.len() >= 10,
        "too few zero crossings: {}",
        crossings.len()
    );
    (crossings[crossings.len() - 1] - crossings[0]) / (crossings.len() - 1) as f64
}

/// At 20 steps per period the Bathe–Noh scheme (p = 1/2) elongates the period by about 0.41 %
/// (below 0.5 %).  The old invented 0.9/0.1 blend with a first-order velocity overwrite is
/// both less accurate and not the Bathe method.
#[test]
fn bathe_period_error_is_below_half_a_percent_at_twenty_steps_per_period() {
    let sys = sdof(1.0, 1.0);
    let bathe = BatheIntegrator::new();
    let dt = 2.0 * PI / 20.0;
    let (mut u, mut v) = (vec![1.0], vec![0.0]);
    let mut a = sys.initial_acceleration(&u, &v, &[0.0]).expect("a0");
    let mut trace = vec![1.0];
    for _ in 0..(20 * 20) {
        bathe
            .step(&sys, &mut u, &mut v, &mut a, &[0.0], dt)
            .expect("Bathe step");
        trace.push(u[0]);
    }
    let period_error = measured_period(&trace, dt) / (2.0 * PI) - 1.0;
    assert!(
        period_error > 0.0 && period_error < 0.005,
        "relative period error {period_error}"
    );
}

/// The composite scheme is unconditionally stable and annihilates unresolved modes: with
/// `ω·Δt = 100` the oscillation never exceeds its initial amplitude and is gone in a few
/// steps.  The old explicit sub-steps diverge there.
#[test]
fn bathe_is_stable_and_dissipative_for_huge_steps() {
    let sys = sdof(1.0, 1.0);
    for split_ratio in [0.2, 0.5, 0.8] {
        let bathe = BatheIntegrator::with_split_ratio(split_ratio);
        let (mut u, mut v, mut a) = (vec![1.0], vec![0.0], vec![-1.0]);
        let mut peak = 1.0_f64;
        for _ in 0..200 {
            bathe
                .step(&sys, &mut u, &mut v, &mut a, &[0.0], 100.0)
                .expect("Bathe step");
            peak = peak.max(u[0].abs());
        }
        assert!(peak <= 1.0 + 1e-9, "p = {split_ratio}: peak |u| = {peak}");
        assert!(
            u[0].abs() < 1e-6 && v[0].abs() < 1e-6,
            "p = {split_ratio}: u = {}, v = {}",
            u[0],
            v[0]
        );
    }
}

/// Second-order accuracy for free vibration, for several split ratios.
#[test]
fn bathe_free_vibration_converges_with_second_order() {
    let sys = sdof(1.0, 1.0);
    for split_ratio in [0.3, 0.5, 0.7] {
        let bathe = BatheIntegrator::with_split_ratio(split_ratio);
        let errors: Vec<f64> = [20_usize, 40, 80, 160]
            .iter()
            .map(|&per_period| {
                let dt = 2.0 * PI / per_period as f64;
                let (mut u, mut v, mut a) = (vec![1.0], vec![0.0], vec![-1.0]);
                for _ in 0..(2 * per_period) {
                    bathe
                        .step(&sys, &mut u, &mut v, &mut a, &[0.0], dt)
                        .expect("Bathe step");
                }
                (u[0] - 1.0).abs() + v[0].abs() // exact after two periods: u = 1, v = 0
            })
            .collect();
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.8..4.3).contains(&ratio),
                "p = {split_ratio}: error ratio {ratio} ({errors:?})"
            );
        }
    }
}

/// The new state satisfies `M·a + C·v + K·u = f(t_{n+1})` for a coupled, damped, loaded system,
/// with the intermediate load taken at `t_n + p·Δt`.
#[test]
fn bathe_satisfies_the_equation_of_motion_at_the_end_of_each_step() {
    let sys = coupled_damped_system();
    let dt = 0.02;
    for split_ratio in [0.3, 0.5, 0.8] {
        let bathe = BatheIntegrator::with_split_ratio(split_ratio);
        let stepper = bathe.stepper(&sys, dt).expect("stepper");
        let mut u = vec![0.01, -0.02, 0.03];
        let mut v = vec![0.1, 0.0, -0.1];
        let mut a = sys
            .initial_acceleration(&u, &v, &coupled_load(0.0))
            .expect("a0");
        for step in 0..200 {
            let t0 = step as f64 * dt;
            let f_mid = coupled_load(t0 + split_ratio * dt);
            let f_end = coupled_load(t0 + dt);
            stepper
                .step_with_loads(&mut u, &mut v, &mut a, &f_mid, &f_end)
                .expect("Bathe step");
            let residual = motion_residual(&sys, &u, &v, &a, &f_end);
            assert!(
                residual < 1e-9,
                "p = {split_ratio}, step {step}: equilibrium residual {residual}"
            );
        }
    }
}

/// Forced response with the two loads of the composite step converges with second order.
#[test]
fn bathe_forced_response_is_second_order() {
    let omega_f = 0.7;
    let load = |t: f64| (omega_f * t).sin();
    let exact = |t: f64| ((omega_f * t).sin() - omega_f * t.sin()) / (1.0 - omega_f * omega_f);
    let sys = sdof(1.0, 1.0);
    let t_end = 20.0;
    for split_ratio in [0.3, 0.5] {
        let bathe = BatheIntegrator::with_split_ratio(split_ratio);
        let errors: Vec<f64> = [200_usize, 400, 800]
            .iter()
            .map(|&steps| {
                let dt = t_end / steps as f64;
                let stepper = bathe.stepper(&sys, dt).expect("stepper");
                let (mut u, mut v) = (vec![0.0], vec![0.0]);
                let mut a = sys.initial_acceleration(&u, &v, &[load(0.0)]).expect("a0");
                for n in 0..steps {
                    let t0 = n as f64 * dt;
                    stepper
                        .step_with_loads(
                            &mut u,
                            &mut v,
                            &mut a,
                            &[load(t0 + split_ratio * dt)],
                            &[load(t0 + dt)],
                        )
                        .expect("Bathe step");
                }
                (u[0] - exact(t_end)).abs()
            })
            .collect();
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.8..4.3).contains(&ratio),
                "p = {split_ratio}: error ratio {ratio} ({errors:?})"
            );
        }
    }
}

/// Invalid split ratios and time steps are rejected without panicking.
#[test]
fn bathe_reports_invalid_input() {
    let sys = sdof(1.0, 1.0);
    let (mut u, mut v, mut a) = (vec![1.0], vec![0.0], vec![-1.0]);
    for split_ratio in [0.0, 1.0, -0.5, 1.5, f64::NAN] {
        let bathe = BatheIntegrator { split_ratio };
        let result = bathe.step(&sys, &mut u, &mut v, &mut a, &[0.0], 0.1);
        assert!(
            matches!(result, Err(Error::InvalidInput(_))),
            "p = {split_ratio}: {result:?}"
        );
    }
    let result = BatheIntegrator::new().step(&sys, &mut u, &mut v, &mut a, &[0.0], 0.0);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");
    assert_eq!((u[0], v[0], a[0]), (1.0, 0.0, -1.0));
}

// ─────────────────────────────────────────────────────────────────────────────
// HRZ mass lumping (F1144)
// ─────────────────────────────────────────────────────────────────────────────

/// Consistent mass matrix of a linear tetrahedron with volume 6.
fn tetrahedron_mass(density: f64) -> (Vec<Vec<f64>>, f64) {
    let nodes = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 6.0],
    ];
    let matrix = element_mass_matrix(&nodes, density);
    (matrix.iter().map(|row| row.to_vec()).collect(), 6.0)
}

/// Lumping a tetrahedral consistent mass gives every DOF `ρV/4`; the lumped masses sum to
/// `3ρV` (the element mass once per direction).  The old code divided the total by `n_dim` and
/// returned `ρV/12`.
#[test]
fn hrz_lumping_of_a_tetrahedron_gives_a_quarter_of_the_element_mass_per_dof() {
    let density = 2.5;
    let (matrix, volume) = tetrahedron_mass(density);
    let lumped = lump_mass_hrz(&matrix, 3);
    assert_eq!(lumped.len(), 12);
    for (i, value) in lumped.iter().enumerate() {
        assert!(
            (value - density * volume / 4.0).abs() < 1e-12,
            "lumped[{i}] = {value}, expected {}",
            density * volume / 4.0
        );
    }
    let total: f64 = lumped.iter().sum();
    assert!(
        (total - 3.0 * density * volume).abs() < 1e-12,
        "total {total}"
    );
    for direction in 0..3 {
        let carried: f64 = lumped.iter().skip(direction).step_by(3).sum();
        assert!((carried - density * volume).abs() < 1e-12);
    }
}

/// HRZ scaling is done per direction, so directions with different mass keep their own total.
#[test]
fn hrz_lumping_preserves_the_mass_of_every_direction() {
    // Two nodes with (x, y) DOFs ordered [n0x, n0y, n1x, n1y]; x carries 6a, y carries 6b.
    let (a, b) = (1.2, 0.4);
    let matrix = vec![
        vec![2.0 * a, 0.0, a, 0.0],
        vec![0.0, 2.0 * b, 0.0, b],
        vec![a, 0.0, 2.0 * a, 0.0],
        vec![0.0, b, 0.0, 2.0 * b],
    ];
    let lumped = lump_mass_hrz(&matrix, 2);
    let expected = [3.0 * a, 3.0 * b, 3.0 * a, 3.0 * b];
    for (value, want) in lumped.iter().zip(expected.iter()) {
        assert!((value - want).abs() < 1e-12, "{lumped:?}");
    }
}

/// Scalar problems, empty input, an all-zero diagonal and `n_dim = 0` are handled gracefully.
#[test]
fn hrz_lumping_edge_cases() {
    let matrix = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
    assert_eq!(lump_mass_hrz(&matrix, 1), vec![3.0, 3.0]);
    assert_eq!(lump_mass_hrz(&matrix, 0), vec![3.0, 3.0]);
    assert!(lump_mass_hrz(&[], 3).is_empty());
    let zero = vec![vec![0.0; 3]; 3];
    assert_eq!(lump_mass_hrz(&zero, 3), vec![0.0; 3]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Explicit time step estimate (F1146)
// ─────────────────────────────────────────────────────────────────────────────

/// Stiffness of a fixed-free chain of unit springs: `tridiag(−1, 2, −1)` with a free end.
fn fixed_free_chain(n: usize) -> Vec<Vec<f64>> {
    let mut k = vec![vec![0.0; n]; n];
    for i in 0..n {
        k[i][i] = 2.0;
        if i > 0 {
            k[i][i - 1] = -1.0;
        }
        if i + 1 < n {
            k[i][i + 1] = -1.0;
        }
    }
    k[n - 1][n - 1] = 1.0;
    k
}

/// Exact `ω_max` of `K·φ = ω²·M·φ` for a diagonal `M`.
fn exact_omega_max(mass: &[f64], stiffness: &[Vec<f64>]) -> f64 {
    let n = mass.len();
    let scaled: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| stiffness[i][j] / (mass[i] * mass[j]).sqrt())
                .collect()
        })
        .collect();
    let (values, _) = jacobi_eigen_dense(&scaled, n);
    values[n - 1].sqrt()
}

/// The estimate is an upper bound of the true `ω_max` for coupled systems (the old `K_ii/M_ii`
/// ratio was `√2` too small for a uniform chain) and stays reasonably tight.
#[test]
fn estimate_omega_max_is_an_upper_bound_of_the_true_frequency() {
    for n in [3_usize, 8, 20] {
        let stiffness = fixed_free_chain(n);
        let masses: Vec<Vec<f64>> = vec![
            vec![1.0; n],
            (0..n).map(|i| 1.0 + 0.5 * (i as f64).sin().abs()).collect(),
        ];
        for mass in masses {
            let estimate = estimate_omega_max(&mass, &stiffness);
            let exact = exact_omega_max(&mass, &stiffness);
            assert!(
                estimate >= exact - 1e-12,
                "n = {n}: estimate {estimate} below the exact {exact}"
            );
            assert!(
                estimate < 1.5 * exact,
                "n = {n}: estimate {estimate} vs exact {exact}"
            );
        }
    }
    // A decoupled system is exact.
    let estimate = estimate_omega_max(&[1.0, 2.0], &[vec![100.0, 0.0], vec![0.0, 200.0]]);
    assert!((estimate - 10.0).abs() < 1e-12);
}

/// Central difference at `safe_explicit_dt(.., 0.95)` of a 20-DOF chain stays bounded for 10⁴
/// steps for every excitation.  The old estimate gave `Δt ≈ 1.34` against a critical 1.003 and
/// blew up within about 200 steps.
#[test]
fn safe_explicit_dt_keeps_central_difference_bounded_on_a_spring_chain() {
    let n = 20;
    let stiffness = fixed_free_chain(n);
    let mass = vec![1.0; n];
    let dt = safe_explicit_dt(&mass, &stiffness, 0.95);
    let critical = 2.0 / exact_omega_max(&mass, &stiffness);
    assert!(
        dt < critical,
        "dt {dt} must be below the critical step {critical}"
    );
    assert!(
        (dt - 0.95).abs() < 1e-12,
        "Gershgorin gives ω_max ≤ 2, so dt = 0.95, got {dt}"
    );

    let sys = DynamicsSystem::new(n, mass, stiffness);
    let excitations: Vec<Vec<f64>> = vec![
        (0..n).map(|i| if i == n - 1 { 1.0 } else { 0.0 }).collect(),
        (0..n)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect(),
    ];
    for u0 in excitations {
        let a0 = sys
            .initial_acceleration(&u0, &vec![0.0; n], &vec![0.0; n])
            .expect("a0");
        let mut u = u0.clone();
        let mut u_prev: Vec<f64> = (0..n).map(|i| u0[i] + 0.5 * dt * dt * a0[i]).collect();
        let mut v = vec![0.0; n];
        let f = vec![0.0; n];
        let mut peak = 0.0_f64;
        for _ in 0..10_000 {
            CentralDifference::step(&sys, &mut u, &mut u_prev, &mut v, &f, dt);
            peak = u.iter().fold(peak, |m, x| m.max(x.abs()));
            assert!(
                peak.is_finite() && peak < 10.0,
                "central difference diverged: peak {peak}"
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Dense subspace eigen-solvers (F1148, F1149)
// ─────────────────────────────────────────────────────────────────────────────

/// `tridiag(−1, 2, −1)` of order `n`; eigenvalues `4·sin²(kπ / (2(n+1)))`, `k = 1..n`.
fn dirichlet_laplacian(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        2.0
                    } else if i.abs_diff(j) == 1 {
                        -1.0
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect()
}

fn dirichlet_eigenvalue(n: usize, k: usize) -> f64 {
    4.0 * (k as f64 * PI / (2.0 * (n as f64 + 1.0))).sin().powi(2)
}

fn diagonal_matrix(d: &[f64]) -> Vec<Vec<f64>> {
    let n = d.len();
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { d[i] } else { 0.0 }).collect())
        .collect()
}

/// `diag(1..10)` with `p = 3` must give `{1, 2, 3}`, not the largest `{8, 9, 10}` the old
/// forward iteration converged to.
#[test]
fn subspace_iteration_returns_the_smallest_eigenpairs() {
    let a = diagonal_matrix(&(1..=10).map(f64::from).collect::<Vec<_>>());
    let pairs = subspace_iteration(&a, 10, 3, 300, 1e-10);
    assert_eq!(pairs.len(), 3);
    for (k, (lambda, vector)) in pairs.iter().enumerate() {
        let expected = (k + 1) as f64;
        assert!((lambda - expected).abs() < 1e-9, "lambda[{k}] = {lambda}");
        assert!(vector[k].abs() > 1.0 - 1e-9, "eigenvector {k} = {vector:?}");
    }
    let residuals = eigenpair_residuals(&a, 10, &pairs).residuals;
    assert!(residuals.iter().all(|r| *r < 1e-8), "{residuals:?}");
    let vectors: Vec<Vec<f64>> = pairs.iter().map(|(_, v)| v.clone()).collect();
    assert_orthonormal(&vectors, 1e-10);
}

/// A non-diagonal matrix whose block (p + guard vectors) is smaller than the space: the
/// 30×30 second-difference matrix with known eigenvalues.
#[test]
fn subspace_iteration_matches_the_analytic_spectrum_of_a_tridiagonal_matrix() {
    let n = 30;
    let a = dirichlet_laplacian(n);
    let pairs = subspace_iteration(&a, n, 4, 500, 1e-10);
    assert_eq!(pairs.len(), 4);
    for (k, (lambda, _)) in pairs.iter().enumerate() {
        let expected = dirichlet_eigenvalue(n, k + 1);
        assert!(
            (lambda - expected).abs() < 1e-9,
            "lambda[{k}] = {lambda}, expected {expected}"
        );
    }
    for window in pairs.windows(2) {
        assert!(window[0].0 <= window[1].0, "eigenvalues must be ascending");
    }
    let residuals = eigenpair_residuals(&a, n, &pairs).residuals;
    assert!(residuals.iter().all(|r| *r < 1e-7), "{residuals:?}");
    let vectors: Vec<Vec<f64>> = pairs.iter().map(|(_, v)| v.clone()).collect();
    assert_orthonormal(&vectors, 1e-9);
    // The k-th Dirichlet mode is sin(k·π·j/(n+1)).
    for (k, vector) in vectors.iter().enumerate() {
        let mode: Vec<f64> = (1..=n)
            .map(|j| ((k + 1) as f64 * PI * j as f64 / (n as f64 + 1.0)).sin())
            .collect();
        let norm = dot(&mode, &mode).sqrt();
        let overlap = dot(vector, &mode).abs() / norm;
        assert!(overlap > 1.0 - 1e-9, "mode {k}: overlap {overlap}");
    }
}

/// A dense matrix with a prescribed spectrum and rotated eigenvectors.
#[test]
fn subspace_iteration_recovers_rotated_eigenvectors() {
    let spectrum: Vec<f64> = (1..=12).map(|i| 0.5 * (i * i) as f64).collect();
    let (a, q) = matrix_with_spectrum(&spectrum);
    let pairs = subspace_iteration(&a, 12, 3, 500, 1e-10);
    assert_eq!(pairs.len(), 3);
    for (k, (lambda, vector)) in pairs.iter().enumerate() {
        assert!(
            (lambda - spectrum[k]).abs() < 1e-8,
            "lambda[{k}] = {lambda}"
        );
        assert!(dot(vector, &q[k]).abs() > 1.0 - 1e-8, "mode {k}");
    }
}

/// Inverse iteration finds the eigenvalues smallest in *magnitude*, which for an indefinite
/// matrix are not the algebraically smallest ones.
#[test]
fn subspace_iteration_finds_the_eigenvalues_closest_to_zero_of_an_indefinite_matrix() {
    let spectrum = [-3.0, -0.5, 0.25, 2.0, 4.0, 7.0, 9.0, 11.0];
    let (a, q) = matrix_with_spectrum(&spectrum);
    let pairs = subspace_iteration(&a, 8, 2, 500, 1e-10);
    assert_eq!(pairs.len(), 2);
    assert!(
        (pairs[0].0 + 0.5).abs() < 1e-8 && (pairs[1].0 - 0.25).abs() < 1e-8,
        "{pairs:?}"
    );
    assert!(dot(&pairs[0].1, &q[1]).abs() > 1.0 - 1e-8);
    assert!(dot(&pairs[1].1, &q[2]).abs() > 1.0 - 1e-8);
}

/// A singular matrix (a free-free chain has a rigid-body mode) is not rejected: the zero
/// eigenpair is returned together with the first flexible mode.
#[test]
fn subspace_iteration_handles_a_singular_matrix_with_a_rigid_body_mode() {
    let n = 12;
    let mut a = dirichlet_laplacian(n);
    a[0][0] = 1.0;
    a[n - 1][n - 1] = 1.0;
    let pairs = subspace_iteration(&a, n, 2, 500, 1e-10);
    assert_eq!(pairs.len(), 2);
    assert!(
        pairs[0].0.abs() < 1e-9,
        "rigid-body eigenvalue {}",
        pairs[0].0
    );
    let flexible = 2.0 - 2.0 * (PI / n as f64).cos();
    assert!(
        (pairs[1].0 - flexible).abs() < 1e-8,
        "{} vs {flexible}",
        pairs[1].0
    );
    let expected = 1.0 / (n as f64).sqrt();
    assert!(
        pairs[0].1.iter().all(|x| (x.abs() - expected).abs() < 1e-6),
        "{:?}",
        pairs[0].1
    );
    let residuals = eigenpair_residuals(&a, n, &pairs).residuals;
    assert!(residuals.iter().all(|r| *r < 1e-7), "{residuals:?}");
}

/// `diag(1, 2, 3, 4)` with `σ = 0` and three modes must return `{1, 2, 3}` with mutually
/// orthogonal vectors; the old implementation returned the eigenpair `1` three times.
#[test]
fn shift_invert_returns_distinct_orthogonal_eigenpairs() {
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);
    let pairs = shift_invert_power(&a, 4, 0.0, 3, 300, 1e-10);
    assert_eq!(pairs.len(), 3);
    for (k, (lambda, vector)) in pairs.iter().enumerate() {
        assert!(
            (lambda - (k + 1) as f64).abs() < 1e-9,
            "lambda[{k}] = {lambda}"
        );
        assert!(vector[k].abs() > 1.0 - 1e-9, "eigenvector {k} = {vector:?}");
    }
    let vectors: Vec<Vec<f64>> = pairs.iter().map(|(_, v)| v.clone()).collect();
    assert_orthonormal(&vectors, 1e-10);
}

/// An interior shift returns the eigenpairs nearest to it, nearest first.
#[test]
fn shift_invert_returns_the_eigenpairs_nearest_to_an_interior_shift() {
    let n = 30;
    let a = dirichlet_laplacian(n);
    let sigma = 1.0;
    let mut all: Vec<f64> = (1..=n).map(|k| dirichlet_eigenvalue(n, k)).collect();
    all.sort_by(|x, y| (x - sigma).abs().total_cmp(&(y - sigma).abs()));
    let pairs = shift_invert_power(&a, n, sigma, 3, 500, 1e-10);
    assert_eq!(pairs.len(), 3);
    for (k, (lambda, _)) in pairs.iter().enumerate() {
        assert!(
            (lambda - all[k]).abs() < 1e-8,
            "lambda[{k}] = {lambda}, expected {}",
            all[k]
        );
    }
    let residuals = eigenpair_residuals(&a, n, &pairs).residuals;
    assert!(residuals.iter().all(|r| *r < 1e-7), "{residuals:?}");
    let vectors: Vec<Vec<f64>> = pairs.iter().map(|(_, v)| v.clone()).collect();
    assert_orthonormal(&vectors, 1e-9);
}

/// A shift that sits exactly on an eigenvalue makes `A − σI` singular; that eigenpair is
/// returned first and the iteration does not fail.
#[test]
fn shift_invert_with_the_shift_on_an_eigenvalue() {
    let a = diagonal_matrix(&[1.0, 2.0, 4.0, 8.0]);
    let pairs = shift_invert_power(&a, 4, 2.0, 2, 300, 1e-10);
    assert_eq!(pairs.len(), 2);
    assert!(
        (pairs[0].0 - 2.0).abs() < 1e-9 && pairs[0].1[1].abs() > 1.0 - 1e-9,
        "{pairs:?}"
    );
    assert!(
        (pairs[1].0 - 1.0).abs() < 1e-9 && pairs[1].1[0].abs() > 1.0 - 1e-9,
        "{pairs:?}"
    );
}

/// A repeated eigenvalue yields an orthonormal basis of its eigenspace, not one vector twice.
#[test]
fn shift_invert_resolves_a_repeated_eigenvalue() {
    let spectrum = [1.0, 1.0, 2.0, 3.0, 5.0, 8.0];
    let (a, q) = matrix_with_spectrum(&spectrum);
    let pairs = shift_invert_power(&a, 6, 0.0, 2, 500, 1e-10);
    assert_eq!(pairs.len(), 2);
    for (lambda, vector) in &pairs {
        assert!((lambda - 1.0).abs() < 1e-8, "lambda = {lambda}");
        let in_eigenspace = dot(vector, &q[0]).powi(2) + dot(vector, &q[1]).powi(2);
        assert!(
            in_eigenspace > 1.0 - 1e-8,
            "vector leaves the eigenspace: {in_eigenspace}"
        );
    }
    assert!(
        dot(&pairs[0].1, &pairs[1].1).abs() < 1e-8,
        "the two vectors must be orthogonal"
    );
}

/// Degenerate requests return nothing, and an exhausted iteration budget still returns the
/// current Ritz pairs.
#[test]
fn subspace_solvers_handle_degenerate_requests() {
    let a = dirichlet_laplacian(6);
    assert!(subspace_iteration(&a, 6, 0, 10, 1e-8).is_empty());
    assert!(subspace_iteration(&a, 0, 3, 10, 1e-8).is_empty());
    assert!(
        subspace_iteration(&a, 9, 3, 10, 1e-8).is_empty(),
        "matrix smaller than n"
    );
    assert!(shift_invert_power(&a, 6, 0.0, 0, 10, 1e-8).is_empty());
    assert!(shift_invert_power(&a, 6, f64::NAN, 2, 10, 1e-8).is_empty());
    let mut bad = a.clone();
    bad[2][3] = f64::NAN;
    assert!(subspace_iteration(&bad, 6, 2, 10, 1e-8).is_empty());

    let rough = subspace_iteration(&a, 6, 3, 1, 0.0);
    assert_eq!(rough.len(), 3);
    assert!(
        rough
            .iter()
            .all(|(lambda, v)| lambda.is_finite() && v.iter().all(|x| x.is_finite()))
    );
    // Asking for more modes than the dimension returns all of them.
    let all = subspace_iteration(&a, 6, 10, 100, 1e-10);
    assert_eq!(all.len(), 6);
    for (k, (lambda, _)) in all.iter().enumerate() {
        assert!((lambda - dirichlet_eigenvalue(6, k + 1)).abs() < 1e-9);
    }
}
