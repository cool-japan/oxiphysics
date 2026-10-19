// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Newmark-β and Wilson-θ steppers of
//! [`DynamicAnalysis`].
//!
//! Both steppers used to solve the *displacement-form* effective system but
//! treated the solution as an acceleration, so the returned acceleration was
//! wrong by a factor `1/(β·dt²)` (Newmark) resp. the Wilson right-hand side
//! lacked the `2·C·v`, `τ/2·C·a` and `f_prev` terms.  The checks below compare
//! against closed-form responses and against the discrete equation of motion.

use oxiphysics_fem::analysis::DynamicAnalysis;

/// Single-DOF oscillator `m·ü + c·u̇ + k·u = f`.
fn sdof(m: f64, c: f64, k: f64) -> DynamicAnalysis {
    DynamicAnalysis::new(vec![m], vec![c], vec![k], 1)
}

/// Dense symmetric-matrix · vector product used to evaluate residuals.
fn matvec(a: &[f64], x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
        .collect()
}

/// Undamped free vibration `m = 1`, `k = 100` (ω = 10 rad/s), `u₀ = 1`,
/// `a₀ = −100`.  The average-acceleration scheme must track `cos(10 t)`
/// (phase error ≈ 8·10⁻⁵ rad after one second) and conserve the energy
/// `½v² + 50u² = 50` to round-off.
#[test]
fn test_newmark_free_vibration_matches_analytic_solution() {
    let da = sdof(1.0, 0.0, 100.0);
    let dt = 1.0e-3;
    let mut u = vec![1.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![-100.0_f64];
    for step in 1..=1000 {
        let (un, vn, an) = da
            .newmark_beta_step(&u, &v, &a, &[0.0], dt, 0.25, 0.5)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        let t = step as f64 * dt;
        assert!(
            (u[0] - (10.0 * t).cos()).abs() < 1e-3,
            "step {step}: u = {}, exact {}",
            u[0],
            (10.0 * t).cos()
        );
        assert!(
            (v[0] + 10.0 * (10.0 * t).sin()).abs() < 1e-2,
            "step {step}: v = {}, exact {}",
            v[0],
            -10.0 * (10.0 * t).sin()
        );
        let energy = 0.5 * v[0] * v[0] + 50.0 * u[0] * u[0];
        assert!(
            (energy - 50.0).abs() < 1e-8,
            "step {step}: energy {energy} drifted from 50"
        );
        assert!(
            (a[0] + 100.0 * u[0]).abs() < 1e-8,
            "step {step}: a = {} violates a = -k u / m for u = {}",
            a[0],
            u[0]
        );
    }
}

/// The new-time-level state must satisfy `M·a + C·v + K·u = f` for a coupled,
/// Rayleigh-damped, forced 3-DOF system and for several (β, γ) — this pins
/// the scaling of the effective system (the pre-fix scheme returned
/// `β·dt²·a`, i.e. an acceleration that was `1/(β·dt²)` too small).
#[test]
fn test_newmark_step_satisfies_the_equation_of_motion() {
    let m = vec![2.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.0];
    let k = vec![
        600.0, -200.0, 0.0, //
        -200.0, 400.0, -200.0, //
        0.0, -200.0, 200.0,
    ];
    let c: Vec<f64> = (0..9).map(|i| 0.02 * m[i] + 0.003 * k[i]).collect();
    let da = DynamicAnalysis::new(m.clone(), c.clone(), k.clone(), 3);
    let u = [0.01, -0.02, 0.03];
    let v = [0.1, 0.0, -0.1];
    let a = [1.0, -2.0, 0.5];
    let f_ext = [5.0, -3.0, 2.0];
    let dt = 0.01;
    for &(beta, gamma) in &[(0.25, 0.5), (1.0 / 6.0, 0.5), (0.3025, 0.6)] {
        let (un, vn, an) = da
            .newmark_beta_step(&u, &v, &a, &f_ext, dt, beta, gamma)
            .expect("integration step");
        let ma = matvec(&m, &an);
        let cv = matvec(&c, &vn);
        let ku = matvec(&k, &un);
        for i in 0..3 {
            let residual = ma[i] + cv[i] + ku[i] - f_ext[i];
            assert!(
                residual.abs() < 1e-9,
                "β = {beta}, γ = {gamma}, DOF {i}: residual {residual}"
            );
            // Newmark kinematic relations.
            let u_expected = u[i] + dt * v[i] + dt * dt * ((0.5 - beta) * a[i] + beta * an[i]);
            let v_expected = v[i] + dt * ((1.0 - gamma) * a[i] + gamma * an[i]);
            assert!((un[i] - u_expected).abs() < 1e-12);
            assert!((vn[i] - v_expected).abs() < 1e-12);
        }
    }
}

/// Halving `dt` must cut the end-of-run error by about four (second order).
#[test]
fn test_newmark_is_second_order_accurate() {
    let error_at_one_second = |dt: f64| {
        let da = sdof(1.0, 0.0, 100.0);
        let steps = (1.0 / dt).round() as usize;
        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-100.0_f64];
        for _ in 0..steps {
            let (un, vn, an) = da
                .newmark_beta_step(&u, &v, &a, &[0.0], dt, 0.25, 0.5)
                .expect("integration step");
            u = un;
            v = vn;
            a = an;
        }
        (u[0] - 10.0_f64.cos()).abs()
    };
    let coarse = error_at_one_second(2.0e-3);
    let fine = error_at_one_second(1.0e-3);
    let ratio = coarse / fine;
    assert!(
        (3.6..4.4).contains(&ratio),
        "error ratio {ratio} (coarse {coarse:e}, fine {fine:e}) is not ≈ 4"
    );
}

/// A step load on a damped oscillator must settle on the static solution
/// `u = f / k = 0.5` (m = 1, c = 2, k = 10 ⇒ decay rate ζω = 1/s).
#[test]
fn test_newmark_static_equilibrium() {
    let da = sdof(1.0, 2.0, 10.0);
    let mut u = vec![0.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![5.0_f64]; // consistent initial acceleration: m·a₀ = f
    let f = vec![5.0_f64];
    for _ in 0..5000 {
        let (un, vn, an) = da
            .newmark_beta_step(&u, &v, &a, &f, 0.01, 0.25, 0.5)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
    }
    assert!((u[0] - 0.5).abs() < 1e-3, "u = {}", u[0]);
    assert!(v[0].abs() < 1e-3, "v = {}", v[0]);
    assert!(a[0].abs() < 1e-2, "a = {}", a[0]);
}

/// Average-acceleration Newmark is unconditionally stable and conserves the
/// energy even when ω·dt = 10 ≫ 1.
#[test]
fn test_newmark_average_acceleration_is_unconditionally_stable() {
    let da = sdof(1.0, 0.0, 1.0e4); // ω = 100 rad/s
    let dt = 0.1;
    let mut u = vec![1.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![-1.0e4_f64];
    for step in 0..200 {
        let (un, vn, an) = da
            .newmark_beta_step(&u, &v, &a, &[0.0], dt, 0.25, 0.5)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        let energy = 0.5 * v[0] * v[0] + 0.5 * 1.0e4 * u[0] * u[0];
        assert!(
            (energy - 5.0e3).abs() < 1e-6 * 5.0e3,
            "step {step}: energy {energy}"
        );
    }
}

/// `m = k = 1`, constant unit load, rest initial state with `a₀ = 1`:
/// the exact response is `u = 1 − cos t`, `v = sin t`.
#[test]
fn test_wilson_theta_forced_oscillator_matches_analytic_solution() {
    let da = sdof(1.0, 0.0, 1.0);
    let dt = 0.01;
    let mut u = vec![0.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![1.0_f64];
    for step in 1..=100 {
        let (un, vn, an) = da
            .wilson_theta_step(&u, &v, &a, &[1.0], &[1.0], dt, 1.4)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        let t = step as f64 * dt;
        assert!(
            (u[0] - (1.0 - t.cos())).abs() < 1e-3,
            "step {step}: u = {}, exact {}",
            u[0],
            1.0 - t.cos()
        );
        assert!(
            (v[0] - t.sin()).abs() < 1e-3,
            "step {step}: v = {}, exact {}",
            v[0],
            t.sin()
        );
        assert!(
            (a[0] - t.cos()).abs() < 1e-2,
            "step {step}: a = {}, exact {}",
            a[0],
            t.cos()
        );
    }
}

/// With `K = C = 0` and constant acceleration the linear-acceleration
/// assumption is exact, so Wilson-θ must reproduce `u = t²/2` to round-off.
/// The pre-fix scheme produced `a ≈ 2·10⁻⁵` instead of `1`.
#[test]
fn test_wilson_theta_integrates_constant_acceleration_exactly() {
    let da = sdof(1.0, 0.0, 0.0);
    let dt = 0.01;
    let mut u = vec![0.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![1.0_f64];
    for step in 1..=10 {
        let (un, vn, an) = da
            .wilson_theta_step(&u, &v, &a, &[1.0], &[1.0], dt, 1.4)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        let t = step as f64 * dt;
        assert!((a[0] - 1.0).abs() < 1e-10, "a = {}", a[0]);
        assert!((v[0] - t).abs() < 1e-12, "v = {}", v[0]);
        assert!((u[0] - 0.5 * t * t).abs() < 1e-12, "u = {}", u[0]);
    }
}

/// A linearly increasing load `f(t) = t` on a free mass gives `a = t`,
/// `v = t²/2`, `u = t³/6`, which the linear-acceleration assumption captures
/// exactly *provided the load increment `f_ext − f_prev` is used*.
#[test]
fn test_wilson_theta_uses_the_load_increment() {
    let da = sdof(1.0, 0.0, 0.0);
    let dt = 0.02;
    let mut u = vec![0.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![0.0_f64];
    for step in 0..50 {
        let t0 = step as f64 * dt;
        let (un, vn, an) = da
            .wilson_theta_step(&u, &v, &a, &[t0], &[t0 + dt], dt, 1.4)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        let t = t0 + dt;
        assert!((a[0] - t).abs() < 1e-10, "a = {}, exact {t}", a[0]);
        assert!((v[0] - t * t / 2.0).abs() < 1e-10, "v = {}", v[0]);
        assert!((u[0] - t * t * t / 6.0).abs() < 1e-10, "u = {}", u[0]);
    }
}

/// 2-DOF chain `K = [[2, −1], [−1, 2]]`, `M = I` started in its first mode
/// `(1, 1)` (ω₁ = 1) must keep both DOFs equal and follow `cos t`.
#[test]
fn test_wilson_theta_two_dof_mode_shape() {
    let da = DynamicAnalysis::new(
        vec![1.0, 0.0, 0.0, 1.0],
        vec![0.0; 4],
        vec![2.0, -1.0, -1.0, 2.0],
        2,
    );
    let dt = 0.01;
    let mut u = vec![1.0_f64, 1.0];
    let mut v = vec![0.0_f64, 0.0];
    let mut a = vec![-1.0_f64, -1.0];
    let zero = [0.0, 0.0];
    for step in 1..=200 {
        let (un, vn, an) = da
            .wilson_theta_step(&u, &v, &a, &zero, &zero, dt, 1.4)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        let t = step as f64 * dt;
        for (dof, &disp) in u.iter().enumerate() {
            assert!(
                (disp - t.cos()).abs() < 5e-3,
                "step {step}, DOF {dof}: u = {disp}, exact {}",
                t.cos()
            );
        }
        assert!((u[0] - u[1]).abs() < 1e-9, "mode shape distorted");
    }
}

/// Step load on a damped oscillator: settles on `u = f / k = 0.5`.
#[test]
fn test_wilson_theta_static_equilibrium() {
    let da = sdof(1.0, 2.0, 10.0);
    let dt = 0.01;
    let mut u = vec![0.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![5.0_f64];
    for _ in 0..5000 {
        let (un, vn, an) = da
            .wilson_theta_step(&u, &v, &a, &[5.0], &[5.0], dt, 1.4)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
    }
    assert!((u[0] - 0.5).abs() < 1e-3, "u = {}", u[0]);
}

/// Free vibration of an undamped oscillator (`ω = 100 rad/s`, `u₀ = 1`) with the
/// large step `ω·dt = 10`; returns `(max |u|, final u)`.
fn large_step_free_vibration(theta: f64, steps: usize) -> (f64, f64) {
    let da = sdof(1.0, 0.0, 1.0e4);
    let dt = 0.1;
    let mut u = vec![1.0_f64];
    let mut v = vec![0.0_f64];
    let mut a = vec![-1.0e4_f64];
    let mut max_abs_u = 0.0_f64;
    for _ in 0..steps {
        let (un, vn, an) = da
            .wilson_theta_step(&u, &v, &a, &[0.0], &[0.0], dt, theta)
            .expect("integration step");
        u = un;
        v = vn;
        a = an;
        max_abs_u = max_abs_u.max(u[0].abs());
    }
    (max_abs_u, u[0])
}

/// θ = 1.4 is unconditionally stable: even with `ω·dt = 10` the response stays
/// bounded (a finite start-up overshoot, ≈ 14 here) and is numerically damped
/// out, whereas θ = 1 (linear acceleration) is only conditionally stable
/// (`ω·dt < 3.46`) and blows up at the same step.
#[test]
fn test_wilson_theta_is_unconditionally_stable() {
    let (max_abs_u, final_u) = large_step_free_vibration(1.4, 400);
    assert!(
        max_abs_u.is_finite() && max_abs_u < 20.0,
        "overshoot {max_abs_u}"
    );
    assert!(final_u.abs() < 1e-6, "response did not decay: {final_u}");

    let (unstable_max, _) = large_step_free_vibration(1.0, 40);
    assert!(
        unstable_max > 1.0e10,
        "θ = 1 should be unstable at ω·dt = 10, got max |u| = {unstable_max}"
    );
}

// ---------------------------------------------------------------------------
// Error reporting and factorisation reuse
// ---------------------------------------------------------------------------

/// The step matrices used to be solved with a Gaussian elimination that
/// `continue`d past a vanishing pivot and then divided by it, returning NaN /
/// garbage; the factorisation was also rebuilt from scratch on every step.
mod step_matrix_handling {
    use super::*;
    use oxiphysics_fem::analysis::FemError;
    use std::f64::consts::PI;

    /// `m = 1`, `k = (2π)²`: a unit-period oscillator started at `u₀ = 1`.
    /// After one second the average-acceleration scheme (dt = 1 ms) is back at
    /// `u = 1`; its phase error is only ≈ 2·10⁻⁵.
    #[test]
    fn unit_period_oscillator_returns_to_its_start_after_one_second() {
        let stiffness = (2.0 * PI) * (2.0 * PI);
        let da = sdof(1.0, 0.0, stiffness);
        let dt = 1.0e-3;
        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-stiffness];
        for _ in 0..1000 {
            let (un, vn, an) = da
                .newmark_beta_step(&u, &v, &a, &[0.0], dt, 0.25, 0.5)
                .expect("regular step matrix");
            u = un;
            v = vn;
            a = an;
        }
        assert!((u[0] - 1.0).abs() < 1.0e-3, "u(1 s) = {}", u[0]);
        assert!(v[0].abs() < 2.0e-2, "v(1 s) = {}", v[0]);
        // One factorisation served all 1000 steps.
        assert_eq!(da.factorization_count(), 1);
    }

    /// Central difference (β = 0) is an explicit scheme: its step matrix is
    /// just `M + γ·dt·C`, which the solver must accept.
    #[test]
    fn explicit_central_difference_step_is_supported() {
        let da = sdof(1.0, 0.0, 100.0);
        let dt = 1.0e-3;
        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-100.0_f64];
        for step in 1..=1000 {
            let (un, vn, an) = da
                .newmark_beta_step(&u, &v, &a, &[0.0], dt, 0.0, 0.5)
                .expect("explicit step");
            u = un;
            v = vn;
            a = an;
            let t = step as f64 * dt;
            assert!(
                (u[0] - (10.0 * t).cos()).abs() < 5.0e-3,
                "step {step}: u = {}",
                u[0]
            );
        }
    }

    /// `M = C = 0` and β = 0 give the zero step matrix.
    #[test]
    fn singular_newmark_matrix_is_an_error_not_nan() {
        let da = sdof(0.0, 0.0, 1.0);
        let result = da.newmark_beta_step(&[0.0], &[0.0], &[0.0], &[1.0], 0.01, 0.0, 0.5);
        assert!(
            matches!(result, Err(FemError::SingularSystem)),
            "{result:?}"
        );
        // A mass cancelled exactly by the stiffness term: with dt = 1/2 and
        // β = 1/4 (exact in binary) M + β·dt²·K = 1 + (1/16)·(−16) = 0.
        let dt = 0.5_f64;
        let beta = 0.25_f64;
        let cancelling = sdof(1.0, 0.0, -1.0 / (beta * dt * dt));
        let result = cancelling.newmark_beta_step(&[0.0], &[0.0], &[0.0], &[1.0], dt, beta, 0.5);
        assert!(
            matches!(result, Err(FemError::SingularSystem)),
            "{result:?}"
        );
        // Nothing was factorised for a singular matrix.
        assert_eq!(da.factorization_count(), 0);
    }

    /// `K̂ = K + 6/τ²·M` vanishes for `K = −6/τ²` (τ = θ·dt): with dt = 1/2 and
    /// θ = 2, τ = 1 and `K = −6` cancel exactly.
    #[test]
    fn singular_wilson_matrix_is_an_error_not_nan() {
        let dt = 0.5_f64;
        let theta = 2.0_f64;
        let da = sdof(1.0, 0.0, -6.0);
        let result = da.wilson_theta_step(&[0.0], &[0.0], &[0.0], &[0.0], &[1.0], dt, theta);
        assert!(
            matches!(result, Err(FemError::SingularSystem)),
            "{result:?}"
        );
        assert_eq!(da.factorization_count(), 0);
    }

    #[test]
    fn invalid_arguments_are_rejected() {
        let da = sdof(1.0, 0.0, 10.0);
        let zero = [0.0_f64];
        let bad_dt = |dt: f64| da.newmark_beta_step(&zero, &zero, &zero, &zero, dt, 0.25, 0.5);
        for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(bad_dt(dt), Err(FemError::InvalidInput(_))),
                "dt = {dt}"
            );
        }
        assert!(matches!(
            da.newmark_beta_step(&zero, &zero, &zero, &zero, 0.01, -0.1, 0.5),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            da.newmark_beta_step(&zero, &zero, &zero, &zero, 0.01, 0.25, f64::NAN),
            Err(FemError::InvalidInput(_))
        ));
        // Wrong vector length and a non-finite load.
        assert!(matches!(
            da.newmark_beta_step(&[0.0, 0.0], &zero, &zero, &zero, 0.01, 0.25, 0.5),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            da.newmark_beta_step(&zero, &zero, &zero, &[f64::NAN], 0.01, 0.25, 0.5),
            Err(FemError::InvalidInput(_))
        ));
        // θ below 1 is not a Wilson scheme.
        assert!(matches!(
            da.wilson_theta_step(&zero, &zero, &zero, &zero, &zero, 0.01, 0.9),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            da.wilson_theta_step(&zero, &zero, &zero, &zero, &zero, 0.01, f64::NAN),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            da.wilson_theta_step(&zero, &zero, &zero, &[], &zero, 0.01, 1.4),
            Err(FemError::InvalidInput(_))
        ));
        // A system whose matrices do not match `n`.
        let broken = DynamicAnalysis::new(vec![1.0; 3], vec![0.0; 4], vec![1.0; 4], 2);
        assert!(matches!(
            broken.newmark_beta_step(&[0.0; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2], 0.01, 0.25, 0.5),
            Err(FemError::InvalidInput(_))
        ));
    }

    /// The step matrix is factorised once per (scheme, parameters) and reused.
    #[test]
    fn factorisation_is_reused_while_the_step_matrix_is_constant() {
        let da = sdof(1.0, 0.1, 50.0);
        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-50.0_f64];
        assert_eq!(da.factorization_count(), 0);
        for _ in 0..40 {
            let (un, vn, an) = da
                .newmark_beta_step(&u, &v, &a, &[0.0], 0.01, 0.25, 0.5)
                .expect("step");
            u = un;
            v = vn;
            a = an;
        }
        assert_eq!(da.factorization_count(), 1);
        // A different time step changes the matrix ...
        let (un, vn, an) = da
            .newmark_beta_step(&u, &v, &a, &[0.0], 0.02, 0.25, 0.5)
            .expect("step");
        assert_eq!(da.factorization_count(), 2);
        // ... and so do β and γ.
        da.newmark_beta_step(&un, &vn, &an, &[0.0], 0.02, 1.0 / 6.0, 0.5)
            .expect("step");
        assert_eq!(da.factorization_count(), 3);
        // Wilson-θ has its own slot, so interleaving the schemes never refactors.
        let zero = [0.0_f64];
        for _ in 0..5 {
            da.wilson_theta_step(&un, &vn, &an, &zero, &zero, 0.02, 1.4)
                .expect("step");
            da.newmark_beta_step(&un, &vn, &an, &zero, 0.02, 1.0 / 6.0, 0.5)
                .expect("step");
        }
        assert_eq!(da.factorization_count(), 4);
    }

    /// `m`, `c`, `k` are public: mutating them between steps must be honoured,
    /// not answered from a stale factorisation of the old matrices.
    #[test]
    fn mutating_the_system_matrices_invalidates_the_cached_factors() {
        let mut da = sdof(1.0, 0.0, 100.0);
        let (u, v, a) = ([1.0_f64], [0.0_f64], [-100.0_f64]);
        let first = da
            .newmark_beta_step(&u, &v, &a, &[0.0], 0.01, 0.25, 0.5)
            .expect("step");
        assert_eq!(da.factorization_count(), 1);

        da.k[0] = 400.0;
        let after_change = da
            .newmark_beta_step(&u, &v, &a, &[0.0], 0.01, 0.25, 0.5)
            .expect("step");
        assert_eq!(da.factorization_count(), 2);

        let fresh = sdof(1.0, 0.0, 400.0)
            .newmark_beta_step(&u, &v, &a, &[0.0], 0.01, 0.25, 0.5)
            .expect("step");
        assert_eq!(after_change, fresh);
        assert_ne!(after_change, first);

        // Changing it back reproduces the very first answer.
        da.k[0] = 100.0;
        let restored = da
            .newmark_beta_step(&u, &v, &a, &[0.0], 0.01, 0.25, 0.5)
            .expect("step");
        assert_eq!(restored, first);
    }

    /// Reusing the factors must not change the numbers.
    #[test]
    fn cached_steps_equal_steps_of_a_fresh_instance() {
        let m = vec![2.0, 0.0, 0.0, 1.0];
        let c = vec![0.3, -0.1, -0.1, 0.2];
        let k = vec![60.0, -20.0, -20.0, 30.0];
        let reused = DynamicAnalysis::new(m.clone(), c.clone(), k.clone(), 2);
        let mut state = (vec![0.01, -0.02], vec![0.1, 0.0], vec![1.0, -1.0]);
        for step in 0..25 {
            let load = [(0.3 * step as f64).sin(), 0.5];
            let from_reused = reused
                .newmark_beta_step(&state.0, &state.1, &state.2, &load, 0.01, 0.25, 0.5)
                .expect("step");
            let fresh = DynamicAnalysis::new(m.clone(), c.clone(), k.clone(), 2);
            let from_fresh = fresh
                .newmark_beta_step(&state.0, &state.1, &state.2, &load, 0.01, 0.25, 0.5)
                .expect("step");
            assert_eq!(from_reused, from_fresh, "step {step}");
            state = from_reused;
        }
        assert_eq!(reused.factorization_count(), 1);
    }

    #[test]
    fn dynamic_analysis_can_be_shared_between_threads() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<DynamicAnalysis>();

        let da = sdof(1.0, 0.05, 80.0);
        let reference = {
            let mut state = (vec![1.0_f64], vec![0.0_f64], vec![-80.0_f64]);
            for _ in 0..50 {
                state = da
                    .newmark_beta_step(&state.0, &state.1, &state.2, &[0.0], 0.005, 0.25, 0.5)
                    .expect("step");
            }
            state
        };
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..4)
                .map(|_| {
                    scope.spawn(|| {
                        let mut state = (vec![1.0_f64], vec![0.0_f64], vec![-80.0_f64]);
                        for _ in 0..50 {
                            state = da
                                .newmark_beta_step(
                                    &state.0,
                                    &state.1,
                                    &state.2,
                                    &[0.0],
                                    0.005,
                                    0.25,
                                    0.5,
                                )
                                .expect("step");
                        }
                        state
                    })
                })
                .collect();
            for handle in handles {
                let state = handle.join().expect("worker thread");
                assert_eq!(state, reference);
            }
        });
        assert_eq!(da.factorization_count(), 1);
    }
}
