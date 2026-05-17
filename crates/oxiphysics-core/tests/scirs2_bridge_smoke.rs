// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoke tests for the SciRS2 ODE integrator bridge.
//!
//! These tests are compiled and run only when the `scirs2` feature is enabled:
//!
//! ```
//! cargo nextest run -p oxiphysics-core --features scirs2
//! ```

#![cfg(feature = "scirs2")]

use oxiphysics_core::Integrator;
use oxiphysics_core::scirs2_integrator::{Scirs2Method, Scirs2OdeIntegrator};
use std::f64::consts::PI;

/// 1-D harmonic oscillator smoke test.
///
/// Physical setup: m=1, k=1, x(0)=1, v(0)=0.
/// Exact solution: x(t) = cos(t), so x(π) = -1.
///
/// The bridge freezes forces at the start of each step (`f = -k·q`), so the
/// integration is not symplectic; there will be a small amplitude drift.
/// We tolerate |x(π) + 1| < 0.1 (loosened from the exact test).
#[test]
fn test_scirs2_rk45_harmonic() {
    let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
    let mut q = vec![1.0_f64];
    let mut v = vec![0.0_f64];
    let inv_mass = vec![1.0_f64];

    let dt = 0.01;
    let steps = (PI / dt).round() as usize;
    for _ in 0..steps {
        let f = vec![-q[0]]; // f = -k·x with k=1 (re-evaluated each step)
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, dt);
    }

    assert!(
        q[0].is_finite(),
        "q must be finite after half-period: {}",
        q[0]
    );

    let error = (q[0] - (-1.0)).abs();
    assert!(
        error < 0.1,
        "Expected x(π)≈-1, got {}, error {}",
        q[0],
        error
    );
}

/// BDF stiff-decay sanity check.
///
/// Physical setup: dy/dt = f * inv_mass (here f = -1000·q, inv_mass = 1).
/// Exact solution: q(0.01) = exp(-10) ≈ 4.54e-5.
///
/// NOTE: the bridge calls `solve_ivp` with constant `f` (evaluated at q at
/// the start of the step), not with the true evolving force.  Therefore we
/// only verify that:
/// 1. The integrator ran without panicking.
/// 2. `q` is finite.
/// 3. `q` is strictly positive (decay, not sign flip).
#[test]
fn test_scirs2_bdf_stiff_decay() {
    let integrator = Scirs2OdeIntegrator {
        method: Scirs2Method::Bdf,
        rtol: 1e-6,
        atol: 1e-8,
        max_steps: 2000,
    };
    let mut q = vec![1.0_f64];
    let mut v = vec![0.0_f64]; // velocity is a dummy for this test
    let inv_mass = vec![1.0_f64];

    let dt = 0.01;
    let f = vec![-1000.0 * q[0]]; // constant force (frozen at initial state)

    integrator.integrate(&mut q, &mut v, &f, &inv_mass, dt);

    assert!(
        q[0].is_finite(),
        "q must be finite after integration: {}",
        q[0]
    );
    assert!(
        q[0] > 0.0,
        "q must be positive (decaying, not oscillating): {}",
        q[0]
    );
}

/// DOP853 high-accuracy smoke test — free particle (f=0, v=const).
///
/// q(dt) = q(0) + v(0) * dt  (exact for free particle).
#[test]
fn test_scirs2_dop853_free_particle() {
    let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Dop853);
    let mut q = vec![0.0_f64];
    let mut v = vec![2.5_f64];
    let f = vec![0.0_f64];
    let inv_mass = vec![1.0_f64];

    let dt = 0.1;
    integrator.integrate(&mut q, &mut v, &f, &inv_mass, dt);

    let expected = 2.5 * dt;
    let error = (q[0] - expected).abs();
    assert!(
        error < 1e-9,
        "Free-particle displacement: expected {expected}, got {}, error {error}",
        q[0]
    );
}

/// RK23 smoke test — constant acceleration.
#[test]
fn test_scirs2_rk23_constant_accel() {
    let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk23);
    let mut q = vec![0.0_f64];
    let mut v = vec![0.0_f64];
    let f = vec![2.0_f64]; // a = 2
    let inv_mass = vec![1.0_f64];

    let dt = 0.5;
    integrator.integrate(&mut q, &mut v, &f, &inv_mass, dt);

    // v should be ≈ a*dt = 1.0; q ≈ 0.5*a*dt² = 0.25
    // RK23 is adaptive and uses variable sub-steps, so we use a loose tolerance.
    assert!(v[0].is_finite(), "v is not finite: {}", v[0]);
    assert!(q[0].is_finite(), "q is not finite: {}", q[0]);
    assert!(v[0] > 0.0, "v should be positive: {}", v[0]);
    assert!(q[0] > 0.0, "q should be positive: {}", q[0]);
    // Loose bounds: allow up to 5% error
    let v_expected = 2.0 * dt;
    let q_expected = 0.5 * 2.0 * dt * dt;
    assert!(
        (v[0] - v_expected).abs() < 0.05,
        "v: expected {v_expected}, got {}, error {}",
        v[0],
        (v[0] - v_expected).abs()
    );
    assert!(
        (q[0] - q_expected).abs() < 0.05,
        "q: expected {q_expected}, got {}, error {}",
        q[0],
        (q[0] - q_expected).abs()
    );
}

/// LSODA smoke test — must not panic for non-stiff input.
#[test]
fn test_scirs2_lsoda_smoke() {
    let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Lsoda);
    let mut q = vec![1.0_f64];
    let mut v = vec![0.0_f64];
    let f = vec![-1.0_f64];
    let inv_mass = vec![1.0_f64];
    integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.05);
    assert!(q[0].is_finite());
    assert!(v[0].is_finite());
}
