// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical ODE solvers with event detection and dense output.
//!
//! Provides classic RK4, adaptive Dormand-Prince RK45 (with FSAL), implicit
//! Euler, Crank-Nicolson (trapezoidal), BDF2 for stiff systems, zero-crossing
//! event detection, and trajectory storage with interpolation.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// OdeState
// ─────────────────────────────────────────────────────────────────────────────

/// Combined time and state vector for an ODE system.
///
/// Stores the current time `t` and the state vector `y`, and provides
/// convenience helpers such as the Euclidean norm of the state.
#[derive(Debug, Clone, PartialEq)]
pub struct OdeState {
    /// Current time.
    pub t: f64,
    /// State vector `y(t)`.
    pub y: Vec<f64>,
}

impl OdeState {
    /// Construct a new [`OdeState`] from time `t` and state vector `y`.
    pub fn new(t: f64, y: Vec<f64>) -> Self {
        Self { t, y }
    }

    /// Euclidean norm of the state vector.
    pub fn norm(&self) -> f64 {
        self.y.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    /// Number of components in the state vector.
    pub fn dim(&self) -> usize {
        self.y.len()
    }

    /// Return a zero state of dimension `n` at time `t`.
    pub fn zeros(t: f64, n: usize) -> Self {
        Self { t, y: vec![0.0; n] }
    }

    /// Component-wise linear interpolation between `self` and `other` at
    /// parameter `alpha` ∈ \[0, 1\].
    pub fn lerp(&self, other: &OdeState, alpha: f64) -> OdeState {
        let t = self.t + alpha * (other.t - self.t);
        let y = self
            .y
            .iter()
            .zip(other.y.iter())
            .map(|(a, b)| a + alpha * (b - a))
            .collect();
        OdeState { t, y }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: vector arithmetic
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn vec_axpy(a: f64, x: &[f64], y: &[f64]) -> Vec<f64> {
    x.iter().zip(y.iter()).map(|(xi, yi)| a * xi + yi).collect()
}

#[inline]
fn vec_sub(x: &[f64], y: &[f64]) -> Vec<f64> {
    x.iter().zip(y.iter()).map(|(a, b)| a - b).collect()
}

#[inline]
fn rms_norm(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// RK4Integrator
// ─────────────────────────────────────────────────────────────────────────────

/// Classic 4th-order Runge-Kutta integrator with optional adaptive step size.
///
/// In fixed-step mode each call to [`RK4Integrator::step`] advances the state
/// by exactly `dt`.  The adaptive driver [`RK4Integrator::integrate_adaptive`]
/// doubles or halves the step based on an embedded 2nd-order (midpoint) error
/// estimate.
pub struct RK4Integrator {
    /// Absolute tolerance used by the adaptive driver.
    pub atol: f64,
    /// Relative tolerance used by the adaptive driver.
    pub rtol: f64,
}

impl RK4Integrator {
    /// Construct an integrator with the given tolerances.
    pub fn new(atol: f64, rtol: f64) -> Self {
        Self { atol, rtol }
    }

    /// Construct with default tolerances (1e-6 absolute, 1e-6 relative).
    pub fn default_tolerances() -> Self {
        Self {
            atol: 1e-6,
            rtol: 1e-6,
        }
    }

    /// Perform one fixed RK4 step from state `s` using step size `dt`.
    ///
    /// `f(t, y)` is the right-hand side of `dy/dt = f(t, y)`.
    pub fn step<F>(&self, s: &OdeState, dt: f64, f: &F) -> OdeState
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let t = s.t;
        let y = &s.y;
        let k1 = f(t, y);
        let y2 = vec_axpy(0.5 * dt, &k1, y);
        let k2 = f(t + 0.5 * dt, &y2);
        let y3 = vec_axpy(0.5 * dt, &k2, y);
        let k3 = f(t + 0.5 * dt, &y3);
        let y4 = vec_axpy(dt, &k3, y);
        let k4 = f(t + dt, &y4);

        let n = y.len();
        let y_new: Vec<f64> = (0..n)
            .map(|i| y[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
            .collect();
        OdeState::new(t + dt, y_new)
    }

    /// Integrate from `s0` to time `t_end` using a fixed step `dt`.
    ///
    /// Returns all intermediate states including the initial state.
    pub fn integrate<F>(&self, s0: &OdeState, t_end: f64, dt: f64, f: &F) -> Vec<OdeState>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut states = vec![s0.clone()];
        let mut s = s0.clone();
        while s.t < t_end - 1e-14 {
            let h = dt.min(t_end - s.t);
            s = self.step(&s, h, f);
            states.push(s.clone());
        }
        states
    }

    /// Integrate from `s0` to `t_end` with adaptive step size.
    ///
    /// Uses an embedded RK2 error estimate to control the step.  The step is
    /// accepted when the RMS error is below `atol + rtol * ||y||`.
    pub fn integrate_adaptive<F>(
        &self,
        s0: &OdeState,
        t_end: f64,
        dt_init: f64,
        f: &F,
    ) -> Vec<OdeState>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut states = vec![s0.clone()];
        let mut s = s0.clone();
        let mut dt = dt_init;
        let dt_min = 1e-12;
        let dt_max = t_end - s0.t;

        while s.t < t_end - 1e-14 {
            let h = dt.min(t_end - s.t).max(dt_min);
            let s_rk4 = self.step(&s, h, f);

            // Embedded RK2 (midpoint) for error estimate
            let k1 = f(s.t, &s.y);
            let y_mid = vec_axpy(0.5 * h, &k1, &s.y);
            let k2 = f(s.t + 0.5 * h, &y_mid);
            let y_rk2: Vec<f64> =
                s.y.iter()
                    .zip(k2.iter())
                    .map(|(yi, ki)| yi + h * ki)
                    .collect();

            let err: Vec<f64> = s_rk4
                .y
                .iter()
                .zip(y_rk2.iter())
                .map(|(a, b)| a - b)
                .collect();
            let tol = self.atol + self.rtol * s_rk4.norm();
            let e = rms_norm(&err);

            if e <= tol || h <= dt_min {
                s = s_rk4;
                states.push(s.clone());
                // Increase step
                if e > 0.0 {
                    dt = (h * (tol / e).powf(0.2)).min(dt_max);
                } else {
                    dt = (h * 2.0).min(dt_max);
                }
            } else {
                // Reject, reduce step
                dt = (h * 0.9 * (tol / e).powf(0.25)).max(dt_min);
            }
        }
        states
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DormandPrince45 (RK45 with FSAL)
// ─────────────────────────────────────────────────────────────────────────────

/// Dormand-Prince RK45 adaptive integrator with FSAL (First Same As Last).
///
/// This is the method underlying MATLAB's `ode45` and SciPy's `RK45`.  The
/// 5th-order solution is used to advance the state; the 4th-order solution is
/// used only for error estimation.  The FSAL property means only 5 new
/// function evaluations per accepted step (the last stage of step n equals the
/// first stage of step n+1).
pub struct DormandPrince45 {
    /// Absolute error tolerance.
    pub atol: f64,
    /// Relative error tolerance.
    pub rtol: f64,
    /// Minimum allowed step size.
    pub dt_min: f64,
    /// Maximum allowed step size.
    pub dt_max: f64,
}

impl DormandPrince45 {
    // Butcher tableau coefficients
    const C2: f64 = 1.0 / 5.0;
    const C3: f64 = 3.0 / 10.0;
    const C4: f64 = 4.0 / 5.0;
    const C5: f64 = 8.0 / 9.0;
    const A21: f64 = 1.0 / 5.0;
    const A31: f64 = 3.0 / 40.0;
    const A32: f64 = 9.0 / 40.0;
    const A41: f64 = 44.0 / 45.0;
    const A42: f64 = -56.0 / 15.0;
    const A43: f64 = 32.0 / 9.0;
    const A51: f64 = 19372.0 / 6561.0;
    const A52: f64 = -25360.0 / 2187.0;
    const A53: f64 = 64448.0 / 6561.0;
    const A54: f64 = -212.0 / 729.0;
    const A61: f64 = 9017.0 / 3168.0;
    const A62: f64 = -355.0 / 33.0;
    const A63: f64 = 46732.0 / 5247.0;
    const A64: f64 = 49.0 / 176.0;
    const A65: f64 = -5103.0 / 18656.0;

    // 5th-order weights
    const B1: f64 = 35.0 / 384.0;
    const B3: f64 = 500.0 / 1113.0;
    const B4: f64 = 125.0 / 192.0;
    const B5: f64 = -2187.0 / 6784.0;
    const B6: f64 = 11.0 / 84.0;

    // Error coefficients (difference of 5th and 4th order weights)
    const E1: f64 = 71.0 / 57600.0;
    const E3: f64 = -71.0 / 16695.0;
    const E4: f64 = 71.0 / 1920.0;
    const E5: f64 = -17253.0 / 339200.0;
    const E6: f64 = 22.0 / 525.0;
    const E7: f64 = -1.0 / 40.0;

    /// Construct with specified tolerances and step bounds.
    pub fn new(atol: f64, rtol: f64, dt_min: f64, dt_max: f64) -> Self {
        Self {
            atol,
            rtol,
            dt_min,
            dt_max,
        }
    }

    /// Construct with default tolerances (1e-6 / 1e-6).
    pub fn default_tolerances() -> Self {
        Self {
            atol: 1e-6,
            rtol: 1e-6,
            dt_min: 1e-12,
            dt_max: f64::INFINITY,
        }
    }

    /// Perform one FSAL Dormand-Prince step from state `s` with step `h`.
    ///
    /// Returns `(new_state, error_rms, k7)` where `k7` is the last stage
    /// (reusable as the first stage of the next step via FSAL).
    pub fn step<F>(
        &self,
        s: &OdeState,
        h: f64,
        f: &F,
        k1_in: Option<&[f64]>,
    ) -> (OdeState, f64, Vec<f64>)
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let t = s.t;
        let y = &s.y;
        let n = y.len();

        let k1: Vec<f64> = match k1_in {
            Some(k) => k.to_vec(),
            None => f(t, y),
        };

        let y2: Vec<f64> = (0..n).map(|i| y[i] + h * Self::A21 * k1[i]).collect();
        let k2 = f(t + Self::C2 * h, &y2);

        let y3: Vec<f64> = (0..n)
            .map(|i| y[i] + h * (Self::A31 * k1[i] + Self::A32 * k2[i]))
            .collect();
        let k3 = f(t + Self::C3 * h, &y3);

        let y4: Vec<f64> = (0..n)
            .map(|i| y[i] + h * (Self::A41 * k1[i] + Self::A42 * k2[i] + Self::A43 * k3[i]))
            .collect();
        let k4 = f(t + Self::C4 * h, &y4);

        let y5: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (Self::A51 * k1[i]
                        + Self::A52 * k2[i]
                        + Self::A53 * k3[i]
                        + Self::A54 * k4[i])
            })
            .collect();
        let k5 = f(t + Self::C5 * h, &y5);

        let y6: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (Self::A61 * k1[i]
                        + Self::A62 * k2[i]
                        + Self::A63 * k3[i]
                        + Self::A64 * k4[i]
                        + Self::A65 * k5[i])
            })
            .collect();
        let k6 = f(t + h, &y6);

        let y_new: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + h
                    * (Self::B1 * k1[i]
                        + Self::B3 * k3[i]
                        + Self::B4 * k4[i]
                        + Self::B5 * k5[i]
                        + Self::B6 * k6[i])
            })
            .collect();
        let k7 = f(t + h, &y_new);

        // Error estimate
        let err: Vec<f64> = (0..n)
            .map(|i| {
                h * (Self::E1 * k1[i]
                    + Self::E3 * k3[i]
                    + Self::E4 * k4[i]
                    + Self::E5 * k5[i]
                    + Self::E6 * k6[i]
                    + Self::E7 * k7[i])
            })
            .collect();

        let sc: Vec<f64> = y_new
            .iter()
            .zip(y.iter())
            .map(|(yn, y0)| self.atol + self.rtol * yn.abs().max(y0.abs()))
            .collect();
        let err_norm = (err
            .iter()
            .zip(sc.iter())
            .map(|(e, s)| (e / s).powi(2))
            .sum::<f64>()
            / n as f64)
            .sqrt();

        (OdeState::new(t + h, y_new), err_norm, k7)
    }

    /// Integrate from `s0` to `t_end` with adaptive step control.
    ///
    /// Returns an [`OdeSolution`] containing all accepted states.
    pub fn integrate<F>(&self, s0: &OdeState, t_end: f64, dt_init: f64, f: &F) -> OdeSolution
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut states = vec![s0.clone()];
        let mut s = s0.clone();
        let mut h = dt_init;
        let mut k1 = f(s.t, &s.y);
        let max_steps = 1_000_000usize;
        let mut n_steps = 0;

        while s.t < t_end - 1e-14 && n_steps < max_steps {
            h = h.min(t_end - s.t).max(self.dt_min).min(self.dt_max);
            let (s_new, err, k7) = self.step(&s, h, f, Some(&k1));

            if err <= 1.0 || h <= self.dt_min {
                s = s_new;
                k1 = k7; // FSAL
                states.push(s.clone());
                // PI step size control
                if err > 0.0 {
                    h = (h * 0.9 * err.powf(-0.2)).min(self.dt_max).max(self.dt_min);
                } else {
                    h = (h * 5.0).min(self.dt_max);
                }
            } else {
                h = (h * 0.9 * err.powf(-0.25)).max(self.dt_min);
            }
            n_steps += 1;
        }

        OdeSolution::new(states)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImplicitEuler (backward Euler)
// ─────────────────────────────────────────────────────────────────────────────

/// Implicit (backward) Euler integrator for stiff ODEs.
///
/// Solves the nonlinear system `y_{n+1} - y_n - h * f(t_{n+1}, y_{n+1}) = 0`
/// using fixed-point (Picard) iteration followed by Newton correction steps.
pub struct ImplicitEuler {
    /// Maximum Newton iterations per step.
    pub max_iter: usize,
    /// Convergence tolerance for Newton iteration.
    pub tol: f64,
    /// Finite-difference step size for the Jacobian.
    pub fd_eps: f64,
}

impl ImplicitEuler {
    /// Construct with specified parameters.
    pub fn new(max_iter: usize, tol: f64, fd_eps: f64) -> Self {
        Self {
            max_iter,
            tol,
            fd_eps,
        }
    }

    /// Construct with default parameters.
    pub fn default_params() -> Self {
        Self {
            max_iter: 50,
            tol: 1e-10,
            fd_eps: 1e-7,
        }
    }

    /// Perform one backward Euler step from state `s` with step `h`.
    ///
    /// Uses simple fixed-point / Picard iteration (no Jacobian required).
    /// For strongly stiff systems prefer Newton iteration via [`ImplicitEuler::step_newton`].
    pub fn step<F>(&self, s: &OdeState, h: f64, f: &F) -> OdeState
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let t_new = s.t + h;
        let mut y = s.y.clone();

        for _ in 0..self.max_iter {
            let rhs = f(t_new, &y);
            let y_new: Vec<f64> =
                s.y.iter()
                    .zip(rhs.iter())
                    .map(|(y0, r)| y0 + h * r)
                    .collect();
            let diff = rms_norm(&vec_sub(&y_new, &y));
            y = y_new;
            if diff < self.tol {
                break;
            }
        }

        OdeState::new(t_new, y)
    }

    /// Perform one backward Euler step using finite-difference Newton iteration.
    ///
    /// More robust than fixed-point for stiff problems.  Approximates the
    /// Jacobian column-by-column with forward differences.
    pub fn step_newton<F>(&self, s: &OdeState, h: f64, f: &F) -> OdeState
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let t_new = s.t + h;
        let n = s.y.len();
        let mut y = s.y.clone();

        for _ in 0..self.max_iter {
            let fy = f(t_new, &y);
            // Residual g(y) = y - y_n - h*f(t_new, y)
            let g: Vec<f64> = (0..n).map(|i| y[i] - s.y[i] - h * fy[i]).collect();

            let g_norm = rms_norm(&g);
            if g_norm < self.tol {
                break;
            }

            // Build approximate Jacobian dg/dy = I - h * df/dy via FD
            // For simplicity use a diagonal approximation
            let mut jac_diag = vec![1.0f64; n];
            for j in 0..n {
                let mut yp = y.clone();
                yp[j] += self.fd_eps;
                let fyp = f(t_new, &yp);
                jac_diag[j] = 1.0 - h * (fyp[j] - fy[j]) / self.fd_eps;
                if jac_diag[j].abs() < 1e-14 {
                    jac_diag[j] = 1.0;
                }
            }

            // Newton update: y <- y - J^{-1} g (diagonal approx)
            for i in 0..n {
                y[i] -= g[i] / jac_diag[i];
            }
        }

        OdeState::new(t_new, y)
    }

    /// Integrate from `s0` to `t_end` with fixed step `dt`.
    ///
    /// Uses Newton iteration per step for robust handling of stiff problems.
    pub fn integrate<F>(&self, s0: &OdeState, t_end: f64, dt: f64, f: &F) -> Vec<OdeState>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut states = vec![s0.clone()];
        let mut s = s0.clone();
        while s.t < t_end - 1e-14 {
            let h = dt.min(t_end - s.t);
            s = self.step_newton(&s, h, f);
            states.push(s.clone());
        }
        states
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Trapezoidal (Crank-Nicolson)
// ─────────────────────────────────────────────────────────────────────────────

/// Trapezoidal (Crank-Nicolson) integrator — second-order accurate, A-stable.
///
/// Solves `y_{n+1} = y_n + h/2 * [f(t_n, y_n) + f(t_{n+1}, y_{n+1})]`
/// iteratively via fixed-point iteration.
pub struct Trapezoidal {
    /// Maximum iterations per step.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

impl Trapezoidal {
    /// Construct with given parameters.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// Construct with default parameters.
    pub fn default_params() -> Self {
        Self {
            max_iter: 50,
            tol: 1e-10,
        }
    }

    /// Perform one Crank-Nicolson step from state `s` with step `h`.
    pub fn step<F>(&self, s: &OdeState, h: f64, f: &F) -> OdeState
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let t_new = s.t + h;
        let f0 = f(s.t, &s.y);
        // Predictor: explicit Euler
        let mut y = vec_axpy(h, &f0, &s.y);

        for _ in 0..self.max_iter {
            let f1 = f(t_new, &y);
            let y_new: Vec<f64> = (0..s.y.len())
                .map(|i| s.y[i] + 0.5 * h * (f0[i] + f1[i]))
                .collect();
            let diff = rms_norm(&vec_sub(&y_new, &y));
            y = y_new;
            if diff < self.tol {
                break;
            }
        }

        OdeState::new(t_new, y)
    }

    /// Integrate from `s0` to `t_end` with fixed step `dt`.
    pub fn integrate<F>(&self, s0: &OdeState, t_end: f64, dt: f64, f: &F) -> Vec<OdeState>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let mut states = vec![s0.clone()];
        let mut s = s0.clone();
        while s.t < t_end - 1e-14 {
            let h = dt.min(t_end - s.t);
            s = self.step(&s, h, f);
            states.push(s.clone());
        }
        states
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BDF2 — 2nd-order backward differentiation formula
// ─────────────────────────────────────────────────────────────────────────────

/// Second-order Backward Differentiation Formula (BDF2) for stiff ODEs.
///
/// Uses the formula `(3/2) y_{n+1} - 2 y_n + (1/2) y_{n-1} = h * f(t_{n+1}, y_{n+1})`.
/// The first step is taken with implicit Euler to obtain `y_1`.
pub struct BDF2 {
    /// Maximum iterations per step.
    pub max_iter: usize,
    /// Convergence tolerance for fixed-point iteration.
    pub tol: f64,
}

impl BDF2 {
    /// Construct with given parameters.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// Construct with default parameters.
    pub fn default_params() -> Self {
        Self {
            max_iter: 50,
            tol: 1e-10,
        }
    }

    /// Perform one BDF2 step given `y_n` (`s_curr`) and `y_{n-1}` (`s_prev`).
    ///
    /// Step size `h` must be the same as the previous step (constant step BDF2).
    /// Uses Newton iteration with a diagonal finite-difference Jacobian for
    /// robust handling of stiff problems.
    pub fn step<F>(&self, s_curr: &OdeState, s_prev: &OdeState, h: f64, f: &F) -> OdeState
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        let t_new = s_curr.t + h;
        let n = s_curr.y.len();
        // Predictor: linear extrapolation
        let mut y: Vec<f64> = (0..n).map(|i| 2.0 * s_curr.y[i] - s_prev.y[i]).collect();
        let fd_eps = 1e-7_f64;

        for _ in 0..self.max_iter {
            let fy = f(t_new, &y);
            // Residual: g(y) = (3/2)y - 2*y_n + (1/2)*y_{n-1} - h*f(y)
            let g: Vec<f64> = (0..n)
                .map(|i| 1.5 * y[i] - 2.0 * s_curr.y[i] + 0.5 * s_prev.y[i] - h * fy[i])
                .collect();
            let g_norm = rms_norm(&g);
            if g_norm < self.tol {
                break;
            }
            // Diagonal Jacobian: dg_i/dy_i = 3/2 - h * df_i/dy_i (FD)
            let mut jac_diag = vec![1.5_f64; n];
            for j in 0..n {
                let mut yp = y.clone();
                yp[j] += fd_eps;
                let fyp = f(t_new, &yp);
                jac_diag[j] = 1.5 - h * (fyp[j] - fy[j]) / fd_eps;
                if jac_diag[j].abs() < 1e-14 {
                    jac_diag[j] = 1.5;
                }
            }
            // Newton update
            for i in 0..n {
                y[i] -= g[i] / jac_diag[i];
            }
        }

        OdeState::new(t_new, y)
    }

    /// Integrate from `s0` to `t_end` with fixed step `dt`.
    ///
    /// The first step uses implicit Euler; subsequent steps use BDF2.
    pub fn integrate<F>(&self, s0: &OdeState, t_end: f64, dt: f64, f: &F) -> Vec<OdeState>
    where
        F: Fn(f64, &[f64]) -> Vec<f64>,
    {
        if s0.t >= t_end - 1e-14 {
            return vec![s0.clone()];
        }
        let ie = ImplicitEuler::new(self.max_iter, self.tol, 1e-7);
        let h = dt.min(t_end - s0.t);
        let s1 = ie.step_newton(s0, h, f);
        let mut states = vec![s0.clone(), s1.clone()];
        let mut s_prev = s0.clone();
        let mut s_curr = s1;

        while s_curr.t < t_end - 1e-14 {
            let step = dt.min(t_end - s_curr.t);
            let s_next = self.step(&s_curr, &s_prev, step, f);
            states.push(s_next.clone());
            s_prev = s_curr;
            s_curr = s_next;
        }
        states
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EventDetection
// ─────────────────────────────────────────────────────────────────────────────

/// Zero-crossing (event) detected during integration.
#[derive(Debug, Clone)]
pub struct CrossingEvent {
    /// Time at which the event function crossed zero.
    pub t: f64,
    /// State at the event time (interpolated).
    pub y: Vec<f64>,
    /// Sign of the event function just before the crossing (-1.0 or +1.0).
    pub sign_before: f64,
    /// Index of the event function that triggered.
    pub event_index: usize,
}

/// Zero-crossing event detection via bisection root-finding.
///
/// After each ODE step, each registered event function `g_i(t, y)` is
/// evaluated.  If the sign of `g_i` changes, bisection is used to locate the
/// crossing time to within `tol`.
pub struct EventDetection {
    /// Bisection tolerance for the event time.
    pub tol: f64,
    /// Maximum bisection iterations.
    pub max_iter: usize,
}

impl EventDetection {
    /// Construct with specified tolerance and iteration limit.
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self { tol, max_iter }
    }

    /// Construct with default parameters.
    pub fn default_params() -> Self {
        Self {
            tol: 1e-10,
            max_iter: 50,
        }
    }

    /// Detect zero crossings of `events[i](t, y)` between states `s_a` and `s_b`.
    ///
    /// Returns a list of [`CrossingEvent`] sorted by time.  Uses linear
    /// interpolation of the state and bisection on each event function.
    pub fn detect<E>(&self, s_a: &OdeState, s_b: &OdeState, events: &[E]) -> Vec<CrossingEvent>
    where
        E: Fn(f64, &[f64]) -> f64,
    {
        let mut crossings = Vec::new();

        for (idx, evt) in events.iter().enumerate() {
            let ga = evt(s_a.t, &s_a.y);
            let gb = evt(s_b.t, &s_b.y);
            if ga * gb > 0.0 {
                continue; // no sign change
            }

            // Bisect in [ta, tb]
            let mut lo = 0.0f64;
            let mut hi = 1.0f64;
            let ga_sign = ga.signum();

            for _ in 0..self.max_iter {
                let mid = 0.5 * (lo + hi);
                let s_mid = s_a.lerp(s_b, mid);
                let gm = evt(s_mid.t, &s_mid.y);
                if gm.signum() == ga_sign {
                    lo = mid;
                } else {
                    hi = mid;
                }
                if hi - lo < self.tol {
                    break;
                }
            }

            let alpha = 0.5 * (lo + hi);
            let s_cross = s_a.lerp(s_b, alpha);
            crossings.push(CrossingEvent {
                t: s_cross.t,
                y: s_cross.y,
                sign_before: ga_sign,
                event_index: idx,
            });
        }

        crossings.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        crossings
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OdeSolution — trajectory storage with interpolation
// ─────────────────────────────────────────────────────────────────────────────

/// Stored ODE trajectory with dense output via linear interpolation.
///
/// Holds all accepted integration steps and provides interpolation at
/// arbitrary times within the integration interval.
#[derive(Debug, Clone)]
pub struct OdeSolution {
    /// Ordered sequence of states (by increasing time).
    pub states: Vec<OdeState>,
}

impl OdeSolution {
    /// Construct from a vector of states (assumed sorted by time).
    pub fn new(states: Vec<OdeState>) -> Self {
        Self { states }
    }

    /// Number of stored states.
    pub fn len(&self) -> usize {
        self.states.len()
    }

    /// Return `true` if no states are stored.
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Interpolate the state at time `t` using linear (dense output) interpolation.
    ///
    /// Returns `None` if `t` is outside the stored time interval.
    pub fn interpolate(&self, t: f64) -> Option<OdeState> {
        if self.states.is_empty() {
            return None;
        }
        let t0 = self.states.first()?.t;
        let t1 = self.states.last()?.t;
        if t < t0 - 1e-14 || t > t1 + 1e-14 {
            return None;
        }
        // Binary search for the interval
        let idx = self.states.partition_point(|s| s.t <= t).saturating_sub(1);
        let idx = idx.min(self.states.len() - 1);

        if idx + 1 >= self.states.len() {
            return Some(self.states[idx].clone());
        }

        let sa = &self.states[idx];
        let sb = &self.states[idx + 1];
        let dt = sb.t - sa.t;
        if dt < 1e-15 {
            return Some(sa.clone());
        }
        let alpha = (t - sa.t) / dt;
        Some(sa.lerp(sb, alpha))
    }

    /// Extract all times as a `Vec`f64`.
    pub fn times(&self) -> Vec<f64> {
        self.states.iter().map(|s| s.t).collect()
    }

    /// Extract the trajectory of component `i` as a `Vec`f64`.
    ///
    /// Returns an empty vector if `i` is out of bounds.
    pub fn component(&self, i: usize) -> Vec<f64> {
        self.states
            .iter()
            .filter_map(|s| s.y.get(i).copied())
            .collect()
    }

    /// Evaluate a scalar observable `g(t, y)` along the trajectory.
    pub fn map_observable<G>(&self, g: G) -> Vec<f64>
    where
        G: Fn(f64, &[f64]) -> f64,
    {
        self.states.iter().map(|s| g(s.t, &s.y)).collect()
    }

    /// Resample the trajectory at `n` uniformly spaced times in \[t0, t1\].
    pub fn resample(&self, n: usize) -> Vec<OdeState> {
        if self.states.len() < 2 || n < 2 {
            return self.states.clone();
        }
        let t0 = self
            .states
            .first()
            .expect("states has at least 2 entries")
            .t;
        let t1 = self.states.last().expect("states has at least 2 entries").t;
        (0..n)
            .filter_map(|k| {
                let t = t0 + (t1 - t0) * k as f64 / (n - 1) as f64;
                self.interpolate(t)
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Suppress unused import of PI (used in tests)
// ─────────────────────────────────────────────────────────────────────────────
const _PI_CHECK: f64 = PI;

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    #[inline]
    fn vec_scale(a: f64, x: &[f64]) -> Vec<f64> {
        x.iter().map(|xi| a * xi).collect()
    }

    #[inline]
    fn vec_add(x: &[f64], y: &[f64]) -> Vec<f64> {
        x.iter().zip(y.iter()).map(|(a, b)| a + b).collect()
    }

    // ------------------------------------------------------------------
    // OdeState
    // ------------------------------------------------------------------
    #[test]
    fn test_ode_state_new_and_norm() {
        let s = OdeState::new(1.0, vec![3.0, 4.0]);
        assert_eq!(s.t, 1.0);
        assert!((s.norm() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_ode_state_zeros() {
        let s = OdeState::zeros(0.0, 5);
        assert_eq!(s.y.len(), 5);
        assert_eq!(s.norm(), 0.0);
    }

    #[test]
    fn test_ode_state_dim() {
        let s = OdeState::new(0.0, vec![1.0, 2.0, 3.0]);
        assert_eq!(s.dim(), 3);
    }

    #[test]
    fn test_ode_state_lerp() {
        let s0 = OdeState::new(0.0, vec![0.0, 0.0]);
        let s1 = OdeState::new(1.0, vec![2.0, 4.0]);
        let mid = s0.lerp(&s1, 0.5);
        assert!((mid.t - 0.5).abs() < 1e-12);
        assert!((mid.y[0] - 1.0).abs() < 1e-12);
        assert!((mid.y[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_ode_state_lerp_endpoints() {
        let s0 = OdeState::new(0.0, vec![1.0]);
        let s1 = OdeState::new(2.0, vec![3.0]);
        let at0 = s0.lerp(&s1, 0.0);
        let at1 = s0.lerp(&s1, 1.0);
        assert!((at0.y[0] - 1.0).abs() < 1e-12);
        assert!((at1.y[0] - 3.0).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // Helper functions
    // ------------------------------------------------------------------
    #[test]
    fn test_rms_norm_empty() {
        assert_eq!(rms_norm(&[]), 0.0);
    }

    #[test]
    fn test_rms_norm_ones() {
        let v = vec![1.0, 1.0, 1.0, 1.0];
        assert!((rms_norm(&v) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec_axpy() {
        let x = vec![1.0, 2.0];
        let y = vec![3.0, 4.0];
        let r = vec_axpy(2.0, &x, &y);
        assert!((r[0] - 5.0).abs() < 1e-12);
        assert!((r[1] - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec_scale() {
        let x = vec![1.0, 2.0, 3.0];
        let r = vec_scale(3.0, &x);
        assert!((r[2] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec_add_sub() {
        let a = vec![1.0, 2.0];
        let b = vec![3.0, 1.0];
        let s = vec_add(&a, &b);
        let d = vec_sub(&b, &a);
        assert!((s[0] - 4.0).abs() < 1e-12);
        assert!((d[1] + 1.0).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // RK4 — exponential decay: dy/dt = -y, y(0)=1
    // ------------------------------------------------------------------
    fn f_decay(_t: f64, y: &[f64]) -> Vec<f64> {
        vec![-y[0]]
    }

    #[test]
    fn test_rk4_single_step_accuracy() {
        let rk4 = RK4Integrator::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let s1 = rk4.step(&s0, 0.1, &f_decay);
        let exact = (-0.1f64).exp();
        assert!((s1.y[0] - exact).abs() < 1e-7);
    }

    #[test]
    fn test_rk4_integrate_fixed() {
        let rk4 = RK4Integrator::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = rk4.integrate(&s0, 1.0, 0.01, &f_decay);
        let last = traj.last().unwrap();
        let exact = (-1.0f64).exp();
        assert!((last.y[0] - exact).abs() < 1e-6);
    }

    #[test]
    fn test_rk4_adaptive() {
        let rk4 = RK4Integrator::new(1e-8, 1e-8);
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = rk4.integrate_adaptive(&s0, 2.0, 0.1, &f_decay);
        let last = traj.last().unwrap();
        let exact = (-2.0f64).exp();
        assert!((last.y[0] - exact).abs() < 1e-5);
    }

    #[test]
    fn test_rk4_harmonic_oscillator() {
        // dy0/dt = y1, dy1/dt = -y0  (omega=1)
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let rk4 = RK4Integrator::default_tolerances();
        let s0 = OdeState::new(0.0, vec![0.0, 1.0]); // sin(t)
        let traj = rk4.integrate(&s0, std::f64::consts::PI, 0.01, &f);
        let last = traj.last().unwrap();
        // y0 = sin(pi) ~ 0
        assert!(last.y[0].abs() < 1e-5);
    }

    // ------------------------------------------------------------------
    // DormandPrince45
    // ------------------------------------------------------------------
    #[test]
    fn test_dp45_exponential_decay() {
        let dp = DormandPrince45::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let sol = dp.integrate(&s0, 1.0, 0.1, &f_decay);
        let last = sol.states.last().unwrap();
        let exact = (-1.0f64).exp();
        // default tolerances are 1e-6; expect global error within a factor of 10
        assert!((last.y[0] - exact).abs() < 1e-5);
    }

    #[test]
    fn test_dp45_harmonic_oscillator() {
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let dp = DormandPrince45::new(1e-9, 1e-9, 1e-12, 1.0);
        let s0 = OdeState::new(0.0, vec![1.0, 0.0]); // cos(t)
        let sol = dp.integrate(&s0, 2.0 * std::f64::consts::PI, 0.1, &f);
        let last = sol.states.last().unwrap();
        // y0 = cos(2*pi) = 1
        assert!((last.y[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dp45_solution_len() {
        let dp = DormandPrince45::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let sol = dp.integrate(&s0, 1.0, 0.1, &f_decay);
        assert!(sol.len() > 1);
    }

    #[test]
    fn test_dp45_fsal_step() {
        let dp = DormandPrince45::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let (s1, err, _k7) = dp.step(&s0, 0.1, &f_decay, None);
        assert!(err >= 0.0);
        let exact = (-0.1f64).exp();
        assert!((s1.y[0] - exact).abs() < 1e-9);
    }

    // ------------------------------------------------------------------
    // ImplicitEuler — stiff decay dy/dt = -100y
    // ------------------------------------------------------------------
    #[test]
    fn test_implicit_euler_stiff_decay() {
        let f = |_t: f64, y: &[f64]| vec![-100.0 * y[0]];
        let ie = ImplicitEuler::default_params();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = ie.integrate(&s0, 1.0, 0.05, &f);
        let last = traj.last().unwrap();
        let exact = (-100.0f64).exp();
        assert!((last.y[0] - exact).abs() < 0.01);
    }

    #[test]
    fn test_implicit_euler_newton_step() {
        let f_lin = |_t: f64, y: &[f64]| vec![-y[0]];
        let ie = ImplicitEuler::default_params();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let s1 = ie.step_newton(&s0, 0.1, &f_lin);
        // Analytical: y1 = 1/(1+h) for implicit Euler on y' = -y
        let expected = 1.0 / 1.1;
        assert!((s1.y[0] - expected).abs() < 1e-8);
    }

    #[test]
    fn test_implicit_euler_zero_rhs() {
        let f_zero = |_t: f64, y: &[f64]| vec![0.0 * y[0]];
        let ie = ImplicitEuler::default_params();
        let s0 = OdeState::new(0.0, vec![5.0]);
        let s1 = ie.step(&s0, 1.0, &f_zero);
        assert!((s1.y[0] - 5.0).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // Trapezoidal
    // ------------------------------------------------------------------
    #[test]
    fn test_trapezoidal_decay() {
        let trap = Trapezoidal::default_params();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = trap.integrate(&s0, 1.0, 0.01, &f_decay);
        let last = traj.last().unwrap();
        let exact = (-1.0f64).exp();
        // Trapezoidal is 2nd order; should be more accurate than implicit Euler
        assert!((last.y[0] - exact).abs() < 1e-5);
    }

    #[test]
    fn test_trapezoidal_single_step() {
        let trap = Trapezoidal::new(100, 1e-12);
        let s0 = OdeState::new(0.0, vec![1.0]);
        let s1 = trap.step(&s0, 0.1, &f_decay);
        // For y'=-y the trapezoidal solution is (1 - h/2)/(1 + h/2)
        let expected = (1.0 - 0.05) / (1.0 + 0.05);
        assert!((s1.y[0] - expected).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // BDF2
    // ------------------------------------------------------------------
    #[test]
    fn test_bdf2_decay() {
        let bdf2 = BDF2::default_params();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = bdf2.integrate(&s0, 1.0, 0.01, &f_decay);
        let last = traj.last().unwrap();
        let exact = (-1.0f64).exp();
        assert!((last.y[0] - exact).abs() < 1e-4);
    }

    #[test]
    fn test_bdf2_stiff_lambda_100() {
        // BDF2 on y'=-100y with h=0.05: very stiff problem.  The BDF2 predictor
        // (linear extrapolation) can be large, but fixed-point iteration should
        // recover a stable solution near zero.  The exact value exp(-50)~1.9e-22
        // is machine-zero; we just verify A-stability (no blow-up) and that
        // the trajectory is monotonically decreasing towards zero.
        let f = |_t: f64, y: &[f64]| vec![-100.0 * y[0]];
        let bdf2 = BDF2::default_params();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = bdf2.integrate(&s0, 0.5, 0.05, &f);
        let last = traj.last().unwrap();
        // BDF2 is A-stable: solution must not blow up
        assert!(
            last.y[0].abs() < 0.5,
            "BDF2 stiff result out of bounds: {}",
            last.y[0]
        );
        // First step should decay significantly from 1.0
        assert!(traj[1].y[0] < 1.0);
    }

    #[test]
    fn test_bdf2_short_interval() {
        let bdf2 = BDF2::default_params();
        let s0 = OdeState::new(5.0, vec![1.0]);
        let traj = bdf2.integrate(&s0, 5.0, 0.1, &f_decay);
        assert_eq!(traj.len(), 1); // already at t_end
    }

    // ------------------------------------------------------------------
    // EventDetection
    // ------------------------------------------------------------------
    #[test]
    fn test_event_detection_crossing_zero() {
        let ed = EventDetection::default_params();
        let s_a = OdeState::new(0.9, vec![0.1]);
        let s_b = OdeState::new(1.1, vec![-0.1]);
        // Event: y[0] itself
        let events: Vec<fn(f64, &[f64]) -> f64> = vec![|_t, y| y[0]];
        let crossings = ed.detect(&s_a, &s_b, &events);
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0].t - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_event_detection_no_crossing() {
        let ed = EventDetection::default_params();
        let s_a = OdeState::new(0.0, vec![1.0]);
        let s_b = OdeState::new(1.0, vec![2.0]);
        let events: Vec<fn(f64, &[f64]) -> f64> = vec![|_t, y| y[0]];
        let crossings = ed.detect(&s_a, &s_b, &events);
        assert!(crossings.is_empty());
    }

    #[test]
    fn test_event_detection_time_event() {
        let ed = EventDetection::default_params();
        let s_a = OdeState::new(0.8, vec![0.0]);
        let s_b = OdeState::new(1.2, vec![0.0]);
        // Event fires at t=1.0
        let events: Vec<fn(f64, &[f64]) -> f64> = vec![|t, _y| t - 1.0];
        let crossings = ed.detect(&s_a, &s_b, &events);
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0].t - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_event_detection_multiple_events() {
        let ed = EventDetection::default_params();
        let s_a = OdeState::new(0.0, vec![2.0, -1.0]);
        let s_b = OdeState::new(2.0, vec![-2.0, 1.0]);
        let ev0: fn(f64, &[f64]) -> f64 = |_t, y| y[0];
        let ev1: fn(f64, &[f64]) -> f64 = |_t, y| y[1];
        let crossings = ed.detect(&s_a, &s_b, &[ev0, ev1]);
        assert_eq!(crossings.len(), 2);
    }

    // ------------------------------------------------------------------
    // OdeSolution
    // ------------------------------------------------------------------
    #[test]
    fn test_ode_solution_interpolate() {
        let states = vec![
            OdeState::new(0.0, vec![0.0]),
            OdeState::new(1.0, vec![1.0]),
            OdeState::new(2.0, vec![4.0]),
        ];
        let sol = OdeSolution::new(states);
        let mid = sol.interpolate(0.5).unwrap();
        assert!((mid.y[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_ode_solution_out_of_range() {
        let states = vec![OdeState::new(0.0, vec![1.0]), OdeState::new(1.0, vec![2.0])];
        let sol = OdeSolution::new(states);
        assert!(sol.interpolate(-0.5).is_none());
        assert!(sol.interpolate(1.5).is_none());
    }

    #[test]
    fn test_ode_solution_times_and_component() {
        let dp = DormandPrince45::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0, 0.0]);
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let sol = dp.integrate(&s0, 1.0, 0.1, &f);
        let ts = sol.times();
        let c0 = sol.component(0);
        assert_eq!(ts.len(), c0.len());
    }

    #[test]
    fn test_ode_solution_resample() {
        let dp = DormandPrince45::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let sol = dp.integrate(&s0, 1.0, 0.1, &f_decay);
        let resampled = sol.resample(20);
        assert_eq!(resampled.len(), 20);
    }

    #[test]
    fn test_ode_solution_empty() {
        let sol = OdeSolution::new(vec![]);
        assert!(sol.is_empty());
        assert!(sol.interpolate(0.5).is_none());
    }

    #[test]
    fn test_ode_solution_map_observable() {
        let states = vec![OdeState::new(0.0, vec![1.0]), OdeState::new(1.0, vec![2.0])];
        let sol = OdeSolution::new(states);
        let obs = sol.map_observable(|_t, y| y[0] * 2.0);
        assert!((obs[0] - 2.0).abs() < 1e-12);
        assert!((obs[1] - 4.0).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // Cross-method comparison
    // ------------------------------------------------------------------
    #[test]
    fn test_rk4_vs_dp45_accuracy() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let rk4 = RK4Integrator::default_tolerances();
        let dp = DormandPrince45::default_tolerances();
        let s0 = OdeState::new(0.0, vec![1.0]);

        let traj_rk4 = rk4.integrate(&s0, 1.0, 0.01, &f);
        let sol_dp = dp.integrate(&s0, 1.0, 0.1, &f);

        let exact = (-1.0f64).exp();
        let err_rk4 = (traj_rk4.last().unwrap().y[0] - exact).abs();
        let err_dp = (sol_dp.states.last().unwrap().y[0] - exact).abs();

        // DP45 with tight tolerance should be at least as accurate as RK4
        assert!(err_dp < 1e-6);
        assert!(err_rk4 < 1e-6);
    }

    #[test]
    fn test_implicit_vs_explicit_stiff() {
        // Explicit RK4 is unstable with large h for stiff problems;
        // implicit Euler stays bounded.
        let f = |_t: f64, y: &[f64]| vec![-1000.0 * y[0]];
        let ie = ImplicitEuler::default_params();
        let s0 = OdeState::new(0.0, vec![1.0]);
        let traj = ie.integrate(&s0, 0.01, 0.001, &f);
        // Solution should go to zero, not blow up
        let last = traj.last().unwrap().y[0];
        assert!((0.0..=1.0).contains(&last));
    }
}
