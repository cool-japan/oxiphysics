// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

#![cfg(feature = "scirs2")]
//! SciRS2 ODE integrator bridge for the [`Integrator`] trait.
//!
//! This module adapts `scirs2_integrate::ode::solve_ivp` so that any of
//! SciRS2's ODE methods (RK45, RK23, BDF, Radau, LSODA, DOP853) can be used
//! as a drop-in replacement for the built-in integrators.
//!
//! # Feature gate
//!
//! The entire module is compiled only when the crate feature `scirs2` is
//! enabled:
//!
//! ```toml
//! [dependencies.oxiphysics-core]
//! features = ["scirs2"]
//! ```
//!
//! # Design
//!
//! `integrate(q, v, f, inv_mass, dt)` is the interface the engine calls every
//! time step.  We adapt it to `solve_ivp` as follows:
//!
//! 1. Pack the combined state `y = [q₀, …, qₙ₋₁, v₀, …, vₙ₋₁]` into an
//!    `Array1<f64>` of length `2n`.
//! 2. Build a right-hand-side closure: `ẏ = [vᵢ, aᵢ]` where `aᵢ = f[i] *
//!    inv_mass[i]`.  Forces are pre-computed by the caller and held constant
//!    over the step.
//! 3. Call `solve_ivp(rhs, [0.0, dt], y0, options)`.
//! 4. Extract the terminal state from `result.y.last()` and unpack back into
//!    `q` and `v`.
//! 5. On integration failure, log a warning and fall back to a single
//!    explicit-Euler step so that the simulation keeps running.

use ndarray::{Array1, ArrayView1};
use scirs2_integrate::ode::{ODEMethod, ODEOptions, solve_ivp};

use crate::math::Real;
use crate::traits::Integrator;

// ---------------------------------------------------------------------------
// Public method enum
// ---------------------------------------------------------------------------

/// ODE method variants available through the SciRS2 bridge.
///
/// These map 1-to-1 to [`scirs2_integrate::ode::ODEMethod`] variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Scirs2Method {
    /// Dormand-Prince RK4(5) — adaptive, non-stiff (default).
    #[default]
    Rk45,
    /// Bogacki-Shampine RK2(3) — adaptive, non-stiff, low accuracy.
    Rk23,
    /// Backward differentiation formula — implicit, stiff.
    Bdf,
    /// Radau IIA — implicit, L-stable, stiff.
    Radau,
    /// LSODA — auto-switching Adams / BDF.
    Lsoda,
    /// Dormand-Prince DOP853 — 8th-order explicit, high accuracy.
    Dop853,
    /// Enhanced BDF — improved Jacobian handling and adaptive order.
    EnhancedBdf,
    /// Enhanced LSODA — improved stiffness detection and method switching.
    EnhancedLsoda,
}

impl Scirs2Method {
    /// Convert to the scirs2-integrate [`ODEMethod`].
    fn to_ode_method(self) -> ODEMethod {
        match self {
            Scirs2Method::Rk45 => ODEMethod::RK45,
            Scirs2Method::Rk23 => ODEMethod::RK23,
            Scirs2Method::Bdf => ODEMethod::Bdf,
            Scirs2Method::Radau => ODEMethod::Radau,
            Scirs2Method::Lsoda => ODEMethod::LSODA,
            Scirs2Method::Dop853 => ODEMethod::DOP853,
            Scirs2Method::EnhancedBdf => ODEMethod::EnhancedBDF,
            Scirs2Method::EnhancedLsoda => ODEMethod::EnhancedLSODA,
        }
    }
}

// ---------------------------------------------------------------------------
// Integrator struct
// ---------------------------------------------------------------------------

/// An [`Integrator`] that delegates to SciRS2's `solve_ivp`.
///
/// # Example
///
/// ```rust,ignore
/// use oxiphysics_core::scirs2_integrator::{Scirs2Method, Scirs2OdeIntegrator};
/// use oxiphysics_core::traits::Integrator;
///
/// let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
/// let mut q = vec![1.0_f64];
/// let mut v = vec![0.0_f64];
/// let f = vec![-1.0_f64];       // restoring force
/// let inv_mass = vec![1.0_f64];
/// integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.01);
/// ```
#[derive(Debug, Clone)]
pub struct Scirs2OdeIntegrator {
    /// ODE method to use.
    pub method: Scirs2Method,
    /// Relative tolerance for adaptive methods.
    pub rtol: f64,
    /// Absolute tolerance for adaptive methods.
    pub atol: f64,
    /// Maximum number of solver steps per `integrate` call.
    pub max_steps: usize,
}

impl Scirs2OdeIntegrator {
    /// Create a new integrator with default tolerances (`rtol=1e-6`, `atol=1e-9`).
    pub fn new(method: Scirs2Method) -> Self {
        Self {
            method,
            rtol: 1e-6,
            atol: 1e-9,
            max_steps: 1000,
        }
    }

    /// Set relative tolerance (builder-style).
    pub fn with_rtol(mut self, rtol: f64) -> Self {
        self.rtol = rtol;
        self
    }

    /// Set absolute tolerance (builder-style).
    pub fn with_atol(mut self, atol: f64) -> Self {
        self.atol = atol;
        self
    }

    /// Set maximum steps per call (builder-style).
    pub fn with_max_steps(mut self, max_steps: usize) -> Self {
        self.max_steps = max_steps;
        self
    }
}

impl Default for Scirs2OdeIntegrator {
    fn default() -> Self {
        Self::new(Scirs2Method::Rk45)
    }
}

// ---------------------------------------------------------------------------
// Integrator trait implementation
// ---------------------------------------------------------------------------

impl Integrator for Scirs2OdeIntegrator {
    /// Advance positions `q` and velocities `v` by `dt` using SciRS2's
    /// `solve_ivp`.
    ///
    /// # State encoding
    ///
    /// The combined state passed to the ODE solver is
    /// `y = [q[0], …, q[n-1], v[0], …, v[n-1]]` (length `2n`).
    ///
    /// The right-hand side is the first-order system:
    /// ```text
    /// dq/dt = v
    /// dv/dt = f * inv_mass    (forces are pre-computed and constant over dt)
    /// ```
    ///
    /// # Fallback
    ///
    /// If `solve_ivp` returns an error, a warning is emitted to `stderr` and
    /// a single explicit-Euler step is performed instead.
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        if n == 0 || dt <= 0.0 {
            return;
        }

        // Pre-compute accelerations a[i] = f[i] * inv_mass[i] — constant over step.
        let accel: Vec<f64> = (0..n).map(|i| f[i] * inv_mass[i]).collect();

        // Build y0 = [q; v] (length 2n).
        let mut y0_vec = Vec::with_capacity(2 * n);
        y0_vec.extend_from_slice(&q[..n]);
        y0_vec.extend_from_slice(&v[..n]);
        let y0 = Array1::from_vec(y0_vec);

        // RHS closure: reads current v from solver state, applies constant accel.
        // `accel` is moved in; closure must be Clone for solve_ivp.
        let accel_clone = accel.clone();
        let rhs = move |_t: f64, y: ArrayView1<f64>| -> Array1<f64> {
            let mut dy = Array1::zeros(2 * n);
            for i in 0..n {
                // dq/dt = v  (read from solver state, not from outer slice)
                dy[i] = y[n + i];
                // dv/dt = accel (constant over this step)
                dy[n + i] = accel_clone[i];
            }
            dy
        };

        let options = ODEOptions {
            method: self.method.to_ode_method(),
            rtol: self.rtol,
            atol: self.atol,
            max_steps: self.max_steps,
            max_step: Some(dt),
            ..Default::default()
        };

        match solve_ivp(rhs, [0.0_f64, dt], y0, Some(options)) {
            Ok(result) => {
                // Extract terminal state (last entry in result.y).
                if let Some(y_final) = result.y.last() {
                    // Unpack q and v from y_final.
                    for i in 0..n {
                        q[i] = y_final[i];
                        v[i] = y_final[n + i];
                    }
                } else {
                    // solve_ivp returned no steps — fall back to explicit Euler.
                    eprintln!(
                        "[oxiphysics-core/scirs2_integrator] solve_ivp returned empty trajectory; \
                         falling back to explicit Euler."
                    );
                    explicit_euler_step(q, v, &accel, n, dt);
                }
            }
            Err(e) => {
                eprintln!(
                    "[oxiphysics-core/scirs2_integrator] solve_ivp failed ({e}); \
                     falling back to explicit Euler."
                );
                explicit_euler_step(q, v, &accel, n, dt);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Fallback helper
// ---------------------------------------------------------------------------

/// Single explicit-Euler step used as fallback on `solve_ivp` failure.
#[inline]
fn explicit_euler_step(q: &mut [Real], v: &mut [Real], accel: &[Real], n: usize, dt: Real) {
    for i in 0..n {
        q[i] += v[i] * dt;
        v[i] += accel[i] * dt;
    }
}

// ---------------------------------------------------------------------------
// Unit tests (compiled only with the feature)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Sanity check: integrating zero force leaves state unchanged.
    #[test]
    fn test_zero_force_no_change() {
        let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
        let mut q = vec![5.0_f64];
        let mut v = vec![0.0_f64];
        let f = vec![0.0_f64];
        let inv_mass = vec![1.0_f64];
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.1);
        assert!(
            (q[0] - 5.0).abs() < 1e-10,
            "q changed unexpectedly: {}",
            q[0]
        );
        assert!((v[0]).abs() < 1e-10, "v changed unexpectedly: {}", v[0]);
    }

    /// Constant force: verify position grows roughly as q ≈ 0.5*a*dt²
    /// (bridge is not symplectic, but should give bounded results).
    #[test]
    fn test_constant_force_finite() {
        let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
        let mut q = vec![0.0_f64];
        let mut v = vec![0.0_f64];
        let f = vec![1.0_f64];
        let inv_mass = vec![1.0_f64];
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.1);
        assert!(q[0].is_finite(), "q must be finite: {}", q[0]);
        assert!(v[0].is_finite(), "v must be finite: {}", v[0]);
        assert!(
            q[0] > 0.0,
            "q must be positive under positive force: {}",
            q[0]
        );
    }

    /// Harmonic oscillator: x(0)=1, v(0)=0 under f=-x (k=1, m=1).
    /// After half-period (t=π), x ≈ -1.  The bridge uses constant-force
    /// approximation per step so a small amplitude drift is expected; we
    /// tolerate error < 0.1.
    #[test]
    fn test_harmonic_oscillator_rk45() {
        let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
        let mut q = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let inv_mass = vec![1.0_f64];

        let dt = 0.01;
        let steps = (PI / dt).round() as usize;
        for _ in 0..steps {
            let fx = vec![-q[0]]; // restoring force (re-evaluated each step)
            integrator.integrate(&mut q, &mut v, &fx, &inv_mass, dt);
        }

        // After one half-period the oscillator should be near -1.
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

    /// Empty slice: must not panic.
    #[test]
    fn test_empty_state() {
        let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
        let mut q: Vec<f64> = vec![];
        let mut v: Vec<f64> = vec![];
        integrator.integrate(&mut q, &mut v, &[], &[], 0.1);
    }

    /// Zero dt: must not panic and must not change state.
    #[test]
    fn test_zero_dt() {
        let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
        let mut q = vec![3.0_f64];
        let mut v = vec![2.0_f64];
        let f = vec![1.0_f64];
        let inv_mass = vec![1.0_f64];
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.0);
        assert!((q[0] - 3.0).abs() < 1e-15);
        assert!((v[0] - 2.0).abs() < 1e-15);
    }

    /// Multi-DOF: 3-component state integrates without panic.
    #[test]
    fn test_multidof() {
        let integrator = Scirs2OdeIntegrator::new(Scirs2Method::Rk45);
        let mut q = vec![1.0_f64, 2.0, 3.0];
        let mut v = vec![0.1_f64, 0.2, 0.3];
        let f = vec![0.0_f64, -1.0, 0.5];
        let inv_mass = vec![1.0_f64, 0.5, 2.0];
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.05);
        for (i, qi) in q.iter().enumerate() {
            assert!(qi.is_finite(), "q[{i}] is not finite: {qi}");
        }
        for (i, vi) in v.iter().enumerate() {
            assert!(vi.is_finite(), "v[{i}] is not finite: {vi}");
        }
    }
}
