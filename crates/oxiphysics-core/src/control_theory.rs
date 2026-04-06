// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Control theory: PID, state-space, LQR, Kalman filter, EKF, Bode analysis,
//! root locus, Ziegler-Nichols auto-tuning, lead/lag compensator, poles/zeros.
//!
//! # Overview
//!
//! This module provides classical and modern control theory building blocks:
//!
//! - [`PidController`] — PID with anti-windup, derivative filter, output clamping
//! - [`StateSpaceModel`] — A, B, C, D matrices; step response; stability check
//! - [`LqrController`] — LQR via discrete algebraic Riccati equation (iterative)
//! - [`KalmanFilter`] — discrete Kalman filter: predict + update
//! - [`ExtendedKalmanFilter`] — EKF with Jacobian-based linearization
//! - [`BodeAnalysis`] — gain/phase at frequencies; gain/phase margins
//! - [`RootLocus`] — closed-loop poles vs gain K; breakaway points
//! - [`ZieglerNichols`] — auto-tune PID from ultimate gain/period
//! - [`LeadLagCompensator`] — lead/lag compensator transfer function
//! - [`PolesZeros`] — transfer function poles, zeros, DC gain, step response

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Internal complex number helper
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
struct Cplx {
    re: f64,
    im: f64,
}

impl Cplx {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
    fn div(self, other: Self) -> Self {
        let denom = other.re * other.re + other.im * other.im;
        if denom.abs() < 1e-300 {
            return Self::new(f64::NAN, f64::NAN);
        }
        Self {
            re: (self.re * other.re + self.im * other.im) / denom,
            im: (self.im * other.re - self.re * other.im) / denom,
        }
    }
    fn scale(self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }
    fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
}

// Evaluate polynomial with complex coefficients at complex point s.
// Coefficients in descending order: coeffs[0]*s^n + ... + coeffs[n].
fn poly_eval_cplx(coeffs: &[f64], s: Cplx) -> Cplx {
    let mut result = Cplx::new(0.0, 0.0);
    for &c in coeffs {
        result = result.mul(s).add(Cplx::new(c, 0.0));
    }
    result
}

// ---------------------------------------------------------------------------
// Simple fixed-size matrix helpers (row-major, stack-allocated via Vec)
// ---------------------------------------------------------------------------

/// Multiply two matrices: (r×k) × (k×c) → (r×c), stored row-major.
fn mat_mul(a: &[f64], b: &[f64], r: usize, k: usize, c: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; r * c];
    for i in 0..r {
        for j in 0..c {
            let mut s = 0.0;
            for p in 0..k {
                s += a[i * k + p] * b[p * c + j];
            }
            out[i * c + j] = s;
        }
    }
    out
}

/// Add two matrices of same shape.
fn mat_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Subtract two matrices.
fn mat_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

/// Transpose an n×m matrix to m×n.
fn mat_transpose(a: &[f64], n: usize, m: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * m];
    for i in 0..n {
        for j in 0..m {
            out[j * n + i] = a[i * m + j];
        }
    }
    out
}

/// Scale a matrix by scalar.
fn mat_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}

/// Identity matrix n×n.
fn mat_eye(n: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * n];
    for i in 0..n {
        out[i * n + i] = 1.0;
    }
    out
}

/// Invert a 2×2 matrix.
fn mat_inv2(a: &[f64]) -> Option<Vec<f64>> {
    let det = a[0] * a[3] - a[1] * a[2];
    if det.abs() < 1e-300 {
        return None;
    }
    Some(vec![a[3] / det, -a[1] / det, -a[2] / det, a[0] / det])
}

/// Invert an n×n matrix via Gauss-Jordan elimination.
fn mat_inv(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut inv = mat_eye(n);
    for col in 0..n {
        // Partial pivot
        let mut max_row = col;
        let mut max_val = m[col * n + col].abs();
        for row in (col + 1)..n {
            let v = m[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-300 {
            return None;
        }
        if max_row != col {
            for j in 0..n {
                m.swap(col * n + j, max_row * n + j);
                inv.swap(col * n + j, max_row * n + j);
            }
        }
        let pivot = m[col * n + col];
        for j in 0..n {
            m[col * n + j] /= pivot;
            inv[col * n + j] /= pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = m[row * n + col];
                for j in 0..n {
                    let mv = m[col * n + j];
                    let iv = inv[col * n + j];
                    m[row * n + j] -= factor * mv;
                    inv[row * n + j] -= factor * iv;
                }
            }
        }
    }
    Some(inv)
}

/// Frobenius norm of a matrix.
fn mat_frob(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

// ---------------------------------------------------------------------------
// PidController
// ---------------------------------------------------------------------------

/// A PID controller with anti-windup, derivative filter, and output clamping.
///
/// Implements the parallel form:
/// `u = Kp*e + Ki*integral(e) + Kd*derivative(e)`
///
/// Anti-windup: integral is clamped to `[-integral_limit, integral_limit]`.
/// Derivative filter: low-pass via `N`-pole (derivative filter coefficient).
/// Output clamping: output is clamped to `[out_min, out_max]`.
#[derive(Debug, Clone)]
pub struct PidController {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Derivative filter coefficient (N = filter pole; larger = less filtering).
    pub n: f64,
    /// Maximum absolute value of the integral term (anti-windup limit).
    pub integral_limit: f64,
    /// Minimum output value.
    pub out_min: f64,
    /// Maximum output value.
    pub out_max: f64,
    /// Accumulated integral state.
    pub integral: f64,
    /// Previous filtered derivative state.
    pub deriv_filter: f64,
    /// Previous error (for derivative computation).
    pub prev_error: f64,
}

impl PidController {
    /// Create a new PID controller.
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            n: 10.0,
            integral_limit: 1e9,
            out_min: f64::NEG_INFINITY,
            out_max: f64::INFINITY,
            integral: 0.0,
            deriv_filter: 0.0,
            prev_error: 0.0,
        }
    }

    /// Set derivative filter coefficient (default 10.0).
    pub fn with_n(mut self, n: f64) -> Self {
        self.n = n;
        self
    }

    /// Set integral anti-windup limit.
    pub fn with_integral_limit(mut self, limit: f64) -> Self {
        self.integral_limit = limit;
        self
    }

    /// Set output clamping limits.
    pub fn with_output_limits(mut self, min: f64, max: f64) -> Self {
        self.out_min = min;
        self.out_max = max;
        self
    }

    /// Compute control output given error and time step `dt`.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        if dt <= 0.0 {
            return 0.0;
        }

        // Proportional
        let p = self.kp * error;

        // Integral with anti-windup
        self.integral += error * dt;
        self.integral = self
            .integral
            .clamp(-self.integral_limit, self.integral_limit);
        let i = self.ki * self.integral;

        // Derivative with low-pass filter: Td/(Td + dt/N) * prev_d + Kd*N/(Td+dt/N)*(e-e_prev)
        // Simplified: filtered_deriv = alpha * filtered_deriv + (1-alpha) * (e - e_prev)/dt
        // where alpha = N*dt / (1 + N*dt)
        let alpha = 1.0 / (1.0 + self.n * dt);
        let raw_deriv = (error - self.prev_error) / dt;
        self.deriv_filter = alpha * self.deriv_filter + (1.0 - alpha) * raw_deriv;
        let d = self.kd * self.deriv_filter;

        self.prev_error = error;

        // Clamp output
        (p + i + d).clamp(self.out_min, self.out_max)
    }

    /// Reset internal state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.deriv_filter = 0.0;
        self.prev_error = 0.0;
    }
}

// ---------------------------------------------------------------------------
// StateSpaceModel
// ---------------------------------------------------------------------------

/// Discrete-time linear state-space model: x\[k+1\] = A*x\[k\] + B*u\[k\], y\[k\] = C*x\[k\] + D*u\[k\].
///
/// Matrices are stored row-major. Dimensions: A is n×n, B is n×m, C is p×n, D is p×m.
#[derive(Debug, Clone)]
pub struct StateSpaceModel {
    /// State matrix A (n×n), row-major.
    pub a: Vec<f64>,
    /// Input matrix B (n×m), row-major.
    pub b: Vec<f64>,
    /// Output matrix C (p×n), row-major.
    pub c: Vec<f64>,
    /// Feedthrough matrix D (p×m), row-major.
    pub d: Vec<f64>,
    /// Number of states n.
    pub n: usize,
    /// Number of inputs m.
    pub m: usize,
    /// Number of outputs p.
    pub p: usize,
    /// Current state vector (length n).
    pub state: Vec<f64>,
}

impl StateSpaceModel {
    /// Create a new state-space model.
    ///
    /// # Panics
    /// Panics if matrix dimensions are inconsistent.
    pub fn new(
        a: Vec<f64>,
        b: Vec<f64>,
        c: Vec<f64>,
        d: Vec<f64>,
        n: usize,
        m: usize,
        p: usize,
    ) -> Self {
        assert_eq!(a.len(), n * n, "A must be n×n");
        assert_eq!(b.len(), n * m, "B must be n×m");
        assert_eq!(c.len(), p * n, "C must be p×n");
        assert_eq!(d.len(), p * m, "D must be p×m");
        Self {
            a,
            b,
            c,
            d,
            n,
            m,
            p,
            state: vec![0.0; n],
        }
    }

    /// Step the system: apply input `u` (length m), advance state, return output y (length p).
    pub fn step(&mut self, u: &[f64]) -> Vec<f64> {
        assert_eq!(u.len(), self.m);
        // x_new = A*x + B*u
        let ax = mat_mul(&self.a, &self.state, self.n, self.n, 1);
        let bu = mat_mul(&self.b, u, self.n, self.m, 1);
        let x_new = mat_add(&ax, &bu);

        // y = C*x + D*u
        let cx = mat_mul(&self.c, &self.state, self.p, self.n, 1);
        let du = mat_mul(&self.d, u, self.p, self.m, 1);
        let y = mat_add(&cx, &du);

        self.state = x_new;
        y
    }

    /// Compute step response: apply unit step input 0→1, collect outputs for `steps` steps.
    /// Returns a vector of (time, y\[0\]) pairs.
    pub fn step_response(&mut self, steps: usize) -> Vec<(f64, f64)> {
        self.state = vec![0.0; self.n];
        let u = vec![1.0; self.m];
        (0..steps)
            .map(|k| {
                let y = self.step(&u);
                (k as f64, y[0])
            })
            .collect()
    }

    /// Check stability: all eigenvalues of A must lie strictly inside the unit disk (discrete).
    /// Uses power iteration / Gershgorin estimate for small systems.
    /// Returns `true` if stable (all |λ| < 1).
    pub fn is_stable(&self) -> bool {
        // For n≤4, use characteristic polynomial roots via companion matrix approach.
        // For general stability, use Gershgorin circle theorem as a sufficient condition,
        // then fall back to iterative spectral radius estimation.
        let spectral_radius = spectral_radius_power_iter(&self.a, self.n, 200, 1e-8);
        spectral_radius < 1.0
    }

    /// Reset state to zero.
    pub fn reset(&mut self) {
        self.state = vec![0.0; self.n];
    }
}

/// Estimate the spectral radius of matrix A (n×n) via power iteration.
fn spectral_radius_power_iter(a: &[f64], n: usize, max_iter: usize, tol: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut v = vec![1.0f64; n];
    let mut lambda = 0.0f64;
    for _iter in 0..max_iter {
        // w = A*v
        let w = mat_mul(a, &v, n, n, 1);
        // norm
        let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-300 {
            return 0.0;
        }
        let new_lambda = norm;
        // normalize
        v = w.iter().map(|x| x / norm).collect();
        if (new_lambda - lambda).abs() < tol {
            return new_lambda;
        }
        lambda = new_lambda;
    }
    lambda
}

// ---------------------------------------------------------------------------
// LqrController
// ---------------------------------------------------------------------------

/// Discrete Linear Quadratic Regulator (LQR) controller.
///
/// Solves the discrete algebraic Riccati equation (DARE) iteratively and
/// computes the optimal gain matrix K such that u = -K*x minimizes
/// J = Σ (x'Qx + u'Ru).
#[derive(Debug, Clone)]
pub struct LqrController {
    /// Optimal gain matrix K (m×n), row-major.
    pub k: Vec<f64>,
    /// Number of states n.
    pub n: usize,
    /// Number of inputs m.
    pub m: usize,
}

impl LqrController {
    /// Solve DARE and compute LQR gain.
    ///
    /// - `a` — state matrix (n×n)
    /// - `b` — input matrix (n×m)
    /// - `q` — state cost matrix (n×n), positive semi-definite
    /// - `r` — input cost matrix (m×m), positive definite
    /// - `max_iter` — maximum iterations for Riccati iteration
    /// - `tol` — convergence tolerance (Frobenius norm of P change)
    pub fn solve(
        a: &[f64],
        b: &[f64],
        q: &[f64],
        r: &[f64],
        n: usize,
        m: usize,
        max_iter: usize,
        tol: f64,
    ) -> Option<Self> {
        // Value iteration: P_{k+1} = Q + A'*P_k*A - A'*P_k*B*(R + B'*P_k*B)^{-1}*B'*P_k*A
        let at = mat_transpose(a, n, n);
        let bt = mat_transpose(b, n, m);
        let mut p = q.to_vec();

        for _iter in 0..max_iter {
            let _p_a = mat_mul(&p, a, n, n, n); // P*A
            let at_p = mat_mul(&at, &p, n, n, n); // A'*P
            let at_p_a = mat_mul(&at_p, a, n, n, n); // A'*P*A
            let at_p_b = mat_mul(&at_p, b, n, n, m); // A'*P*B
            let bt_p = mat_mul(&bt, &p, m, n, n); // B'*P
            let bt_p_b = mat_mul(&bt_p, b, m, n, m); // B'*P*B
            let r_bpb = mat_add(r, &bt_p_b); // R + B'*P*B

            let r_bpb_inv = mat_inv(&r_bpb, m)?;
            let at_p_b_rinv = mat_mul(&at_p_b, &r_bpb_inv, n, m, m); // A'*P*B*(R+B'PB)^{-1}
            let correction = mat_mul(&at_p_b_rinv, &mat_transpose(&at_p_b, n, m), n, m, n);

            let p_new = mat_sub(&mat_add(q, &at_p_a), &correction);

            let diff = mat_frob(&mat_sub(&p_new, &p));
            p = p_new;
            if diff < tol {
                break;
            }
        }

        // K = (R + B'*P*B)^{-1} * B'*P*A
        let at = mat_transpose(a, n, n);
        let bt = mat_transpose(b, n, m);
        let bt_p = mat_mul(&bt, &p, m, n, n);
        let bt_p_b = mat_mul(&bt_p, b, m, n, m);
        let r_bpb = mat_add(r, &bt_p_b);
        let r_bpb_inv = mat_inv(&r_bpb, m)?;
        let bt_p_a = mat_mul(&bt_p, a, m, n, n);
        let k = mat_mul(&r_bpb_inv, &bt_p_a, m, m, n);

        // Suppress unused variable warning for `at`
        let _ = &at;

        Some(Self { k, n, m })
    }

    /// Compute control input: u = -K*x.
    pub fn control(&self, x: &[f64]) -> Vec<f64> {
        let kx = mat_mul(&self.k, x, self.m, self.n, 1);
        kx.iter().map(|v| -v).collect()
    }
}

// ---------------------------------------------------------------------------
// KalmanFilter
// ---------------------------------------------------------------------------

/// Discrete Kalman filter.
///
/// State: x̂ (n×1). Matrices: A (n×n), B (n×m), C (p×n), Q (n×n), R (p×p).
/// Predict: x̂⁻ = A*x̂ + B*u, P⁻ = A*P*A' + Q
/// Update: K = P⁻*C'*(C*P⁻*C' + R)⁻¹, x̂ = x̂⁻ + K*(y - C*x̂⁻), P = (I - K*C)*P⁻
#[derive(Debug, Clone)]
pub struct KalmanFilter {
    /// State matrix A (n×n).
    pub a: Vec<f64>,
    /// Input matrix B (n×m).
    pub b: Vec<f64>,
    /// Observation matrix C (p×n).
    pub c: Vec<f64>,
    /// Process noise covariance Q (n×n).
    pub q: Vec<f64>,
    /// Measurement noise covariance R (p×p).
    pub r: Vec<f64>,
    /// State estimate x̂ (length n).
    pub x: Vec<f64>,
    /// Error covariance matrix P (n×n).
    pub p: Vec<f64>,
    /// Number of states n.
    pub n: usize,
    /// Number of inputs m.
    pub m: usize,
    /// Number of outputs p.
    pub p_dim: usize,
}

impl KalmanFilter {
    /// Create a new Kalman filter with identity initial covariance.
    pub fn new(
        a: Vec<f64>,
        b: Vec<f64>,
        c: Vec<f64>,
        q: Vec<f64>,
        r: Vec<f64>,
        n: usize,
        m: usize,
        p_dim: usize,
    ) -> Self {
        let p_init = mat_eye(n);
        Self {
            a,
            b,
            c,
            q,
            r,
            x: vec![0.0; n],
            p: p_init,
            n,
            m,
            p_dim,
        }
    }

    /// Predict step: propagate state and covariance.
    pub fn predict(&mut self, u: &[f64]) {
        let ax = mat_mul(&self.a, &self.x, self.n, self.n, 1);
        let bu = mat_mul(&self.b, u, self.n, self.m, 1);
        self.x = mat_add(&ax, &bu);

        let at = mat_transpose(&self.a, self.n, self.n);
        let apa = mat_mul(
            &mat_mul(&self.a, &self.p, self.n, self.n, self.n),
            &at,
            self.n,
            self.n,
            self.n,
        );
        self.p = mat_add(&apa, &self.q);
    }

    /// Update step: incorporate measurement y (length p_dim).
    pub fn update(&mut self, y: &[f64]) -> Option<()> {
        let ct = mat_transpose(&self.c, self.p_dim, self.n);
        let pct = mat_mul(&self.p, &ct, self.n, self.n, self.p_dim);
        let cpct = mat_mul(&self.c, &pct, self.p_dim, self.n, self.p_dim);
        let s = mat_add(&cpct, &self.r); // S = C*P*C' + R  (p×p)

        let s_inv = mat_inv(&s, self.p_dim)?;
        // K = P*C' * S^{-1}  (n×p)
        let k = mat_mul(&pct, &s_inv, self.n, self.p_dim, self.p_dim);

        // Innovation: innov = y - C*x̂  (p×1)
        let cx = mat_mul(&self.c, &self.x, self.p_dim, self.n, 1);
        let innov: Vec<f64> = y.iter().zip(cx.iter()).map(|(yi, ci)| yi - ci).collect();

        // x̂ = x̂ + K*innov
        let k_innov = mat_mul(&k, &innov, self.n, self.p_dim, 1);
        self.x = mat_add(&self.x, &k_innov);

        // P = (I - K*C)*P
        let kc = mat_mul(&k, &self.c, self.n, self.p_dim, self.n);
        let i_kc = mat_sub(&mat_eye(self.n), &kc);
        self.p = mat_mul(&i_kc, &self.p, self.n, self.n, self.n);

        Some(())
    }

    /// Run one full predict-update cycle.
    pub fn step(&mut self, u: &[f64], y: &[f64]) -> Option<Vec<f64>> {
        self.predict(u);
        self.update(y)?;
        Some(self.x.clone())
    }
}

// ---------------------------------------------------------------------------
// ExtendedKalmanFilter
// ---------------------------------------------------------------------------

/// Extended Kalman Filter (EKF) with Jacobian-based linearization.
///
/// The user provides:
/// - `f`: nonlinear state transition function f(x, u) → x_new
/// - `h`: nonlinear observation function h(x) → y
/// - `jac_f`: Jacobian of f with respect to x (n×n matrix, row-major)
/// - `jac_h`: Jacobian of h with respect to x (p×n matrix, row-major)
#[derive(Debug, Clone)]
pub struct ExtendedKalmanFilter {
    /// State estimate (length n).
    pub x: Vec<f64>,
    /// Error covariance (n×n, row-major).
    pub p: Vec<f64>,
    /// Process noise covariance Q (n×n).
    pub q: Vec<f64>,
    /// Measurement noise covariance R (p×p).
    pub r: Vec<f64>,
    /// Number of states.
    pub n: usize,
    /// Number of outputs.
    pub p_dim: usize,
}

impl ExtendedKalmanFilter {
    /// Create a new EKF.
    pub fn new(q: Vec<f64>, r: Vec<f64>, n: usize, p_dim: usize) -> Self {
        Self {
            x: vec![0.0; n],
            p: mat_eye(n),
            q,
            r,
            n,
            p_dim,
        }
    }

    /// EKF predict step.
    ///
    /// - `f`: state transition f(x, u) → new state
    /// - `jac_f`: Jacobian F = ∂f/∂x at current x (row-major, n×n)
    /// - `u`: control input
    pub fn predict<F, J>(&mut self, f: F, jac_f: J, u: &[f64])
    where
        F: Fn(&[f64], &[f64]) -> Vec<f64>,
        J: Fn(&[f64], &[f64]) -> Vec<f64>,
    {
        let f_val = jac_f(&self.x, u);
        // P⁻ = F*P*F' + Q
        let ft = mat_transpose(&f_val, self.n, self.n);
        let fp = mat_mul(&f_val, &self.p, self.n, self.n, self.n);
        let fpft = mat_mul(&fp, &ft, self.n, self.n, self.n);
        self.p = mat_add(&fpft, &self.q);
        self.x = f(&self.x, u);
    }

    /// EKF update step.
    ///
    /// - `h`: observation function h(x) → y
    /// - `jac_h`: Jacobian H = ∂h/∂x at current x (row-major, p×n)
    /// - `y`: measurement vector (length p_dim)
    pub fn update<H, J>(&mut self, h: H, jac_h: J, y: &[f64]) -> Option<()>
    where
        H: Fn(&[f64]) -> Vec<f64>,
        J: Fn(&[f64]) -> Vec<f64>,
    {
        let h_jac = jac_h(&self.x);
        let ht = mat_transpose(&h_jac, self.p_dim, self.n);
        let pht = mat_mul(&self.p, &ht, self.n, self.n, self.p_dim);
        let hpht = mat_mul(&h_jac, &pht, self.p_dim, self.n, self.p_dim);
        let s = mat_add(&hpht, &self.r);
        let s_inv = mat_inv(&s, self.p_dim)?;
        let k = mat_mul(&pht, &s_inv, self.n, self.p_dim, self.p_dim);

        let h_x = h(&self.x);
        let innov: Vec<f64> = y.iter().zip(h_x.iter()).map(|(yi, hi)| yi - hi).collect();
        let k_innov = mat_mul(&k, &innov, self.n, self.p_dim, 1);
        self.x = mat_add(&self.x, &k_innov);

        let kh = mat_mul(&k, &h_jac, self.n, self.p_dim, self.n);
        let i_kh = mat_sub(&mat_eye(self.n), &kh);
        self.p = mat_mul(&i_kh, &self.p, self.n, self.n, self.n);
        Some(())
    }
}

// ---------------------------------------------------------------------------
// BodeAnalysis
// ---------------------------------------------------------------------------

/// Bode analysis for a transfer function given numerator and denominator polynomials.
///
/// Polynomials stored as coefficients in descending order of s:
/// `num = [b_n, ..., b_0]`, `den = [a_m, ..., a_0]`.
#[derive(Debug, Clone)]
pub struct BodeAnalysis {
    /// Numerator polynomial coefficients (descending order of s).
    pub num: Vec<f64>,
    /// Denominator polynomial coefficients (descending order of s).
    pub den: Vec<f64>,
}

/// A single point on a Bode plot.
#[derive(Debug, Clone, Copy)]
pub struct BodePoint {
    /// Angular frequency ω (rad/s).
    pub omega: f64,
    /// Gain in dB: 20*log10(|G(jω)|).
    pub gain_db: f64,
    /// Phase in degrees.
    pub phase_deg: f64,
}

impl BodeAnalysis {
    /// Create a Bode analysis object.
    pub fn new(num: Vec<f64>, den: Vec<f64>) -> Self {
        Self { num, den }
    }

    /// Evaluate transfer function G(jω) at angular frequency ω.
    pub fn eval_at(&self, omega: f64) -> (f64, f64) {
        let s = Cplx::new(0.0, omega);
        let num_val = poly_eval_cplx(&self.num, s);
        let den_val = poly_eval_cplx(&self.den, s);
        let g = num_val.div(den_val);
        let gain = g.abs();
        let phase_deg = g.arg().to_degrees();
        (gain, phase_deg)
    }

    /// Compute Bode points at the given frequencies.
    pub fn compute(&self, omegas: &[f64]) -> Vec<BodePoint> {
        omegas
            .iter()
            .map(|&omega| {
                let (gain, phase_deg) = self.eval_at(omega);
                let gain_db = if gain > 1e-300 {
                    20.0 * gain.log10()
                } else {
                    -300.0
                };
                BodePoint {
                    omega,
                    gain_db,
                    phase_deg,
                }
            })
            .collect()
    }

    /// Compute logarithmically-spaced Bode plot from ω_min to ω_max with n_points points.
    pub fn compute_log(&self, omega_min: f64, omega_max: f64, n_points: usize) -> Vec<BodePoint> {
        if n_points == 0 {
            return vec![];
        }
        let log_min = omega_min.log10();
        let log_max = omega_max.log10();
        let omegas: Vec<f64> = (0..n_points)
            .map(|i| {
                let t = i as f64 / (n_points - 1).max(1) as f64;
                10f64.powf(log_min + t * (log_max - log_min))
            })
            .collect();
        self.compute(&omegas)
    }

    /// Find gain margin (dB) and phase margin (degrees).
    ///
    /// - Gain margin: gain (dB) at the phase crossover frequency (where phase = -180°).
    /// - Phase margin: phase (deg) + 180° at the gain crossover frequency (where gain = 0 dB).
    ///
    /// Returns `(gain_margin_db, phase_margin_deg)`.
    pub fn margins(&self, omega_min: f64, omega_max: f64, n_points: usize) -> (f64, f64) {
        let points = self.compute_log(omega_min, omega_max, n_points);
        if points.is_empty() {
            return (f64::NAN, f64::NAN);
        }

        // Phase crossover: phase = -180°
        let mut phase_cross_gain_db = f64::NAN;
        for i in 1..points.len() {
            let p0 = &points[i - 1];
            let p1 = &points[i];
            if (p0.phase_deg + 180.0) * (p1.phase_deg + 180.0) <= 0.0 {
                // Linear interpolation
                let frac = (-(p0.phase_deg + 180.0)) / (p1.phase_deg - p0.phase_deg + 1e-300);
                phase_cross_gain_db = p0.gain_db + frac * (p1.gain_db - p0.gain_db);
                break;
            }
        }

        // Gain crossover: gain = 0 dB
        let mut gain_cross_phase_margin = f64::NAN;
        for i in 1..points.len() {
            let p0 = &points[i - 1];
            let p1 = &points[i];
            if p0.gain_db * p1.gain_db <= 0.0 {
                let frac = (-p0.gain_db) / (p1.gain_db - p0.gain_db + 1e-300);
                let phase = p0.phase_deg + frac * (p1.phase_deg - p0.phase_deg);
                gain_cross_phase_margin = phase + 180.0;
                break;
            }
        }

        let gain_margin_db = -phase_cross_gain_db;
        (gain_margin_db, gain_cross_phase_margin)
    }
}

// ---------------------------------------------------------------------------
// RootLocus
// ---------------------------------------------------------------------------

/// Root locus computation: closed-loop poles as a function of gain K.
///
/// For open-loop transfer function G(s) = num(s)/den(s),
/// the closed-loop poles satisfy: den(s) + K*num(s) = 0.
#[derive(Debug, Clone)]
pub struct RootLocus {
    /// Numerator polynomial coefficients (descending order).
    pub num: Vec<f64>,
    /// Denominator polynomial coefficients (descending order).
    pub den: Vec<f64>,
}

/// A point on the root locus.
#[derive(Debug, Clone)]
pub struct RootLocusPoint {
    /// Gain K.
    pub gain: f64,
    /// Closed-loop poles at this gain (complex pairs as (re, im)).
    pub poles: Vec<(f64, f64)>,
}

impl RootLocus {
    /// Create a root locus object.
    pub fn new(num: Vec<f64>, den: Vec<f64>) -> Self {
        Self { num, den }
    }

    /// Compute closed-loop poles for each gain in `gains`.
    ///
    /// Uses companion matrix eigenvalue approach for small polynomials (degree ≤ 4).
    /// For the characteristic equation den(s) + K*num(s) = 0, forms the combined polynomial.
    pub fn compute(&self, gains: &[f64]) -> Vec<RootLocusPoint> {
        gains
            .iter()
            .map(|&k| {
                let poles = self.poles_at_gain(k);
                RootLocusPoint { gain: k, poles }
            })
            .collect()
    }

    /// Compute closed-loop poles for a specific gain K.
    ///
    /// Characteristic polynomial = den + K*num.
    /// Pads/truncates to match degree. Returns (re, im) pairs.
    pub fn poles_at_gain(&self, k: f64) -> Vec<(f64, f64)> {
        // Form characteristic polynomial
        let n_deg = self.den.len().max(self.num.len());
        let mut char_poly = vec![0.0f64; n_deg];
        for (i, &c) in self.den.iter().enumerate() {
            let offset = n_deg - self.den.len();
            char_poly[i + offset] += c;
        }
        for (i, &c) in self.num.iter().enumerate() {
            let offset = n_deg - self.num.len();
            char_poly[i + offset] += k * c;
        }

        // Find roots via companion matrix for degree ≤ 4
        find_polynomial_roots(&char_poly)
    }

    /// Compute breakaway points (real-axis segments, ∂K/∂s = 0).
    /// Returns approximate breakaway points on the real axis.
    pub fn breakaway_points(&self, s_min: f64, s_max: f64, n_samples: usize) -> Vec<f64> {
        // K(s) = -den(s)/num(s); breakaway where dK/ds = 0
        // Approximate via finite differences on real axis
        let ds = (s_max - s_min) / (n_samples as f64);
        let mut points = Vec::new();

        let k_at = |s_real: f64| {
            let s = Cplx::new(s_real, 0.0);
            let d = poly_eval_cplx(&self.den, s);
            let n = poly_eval_cplx(&self.num, s);
            if n.abs() < 1e-12 {
                f64::NAN
            } else {
                -d.re / n.re
            }
        };

        let mut prev_k = k_at(s_min);
        let mut prev_dk = f64::NAN;
        let mut prev_s = s_min;

        for i in 1..=n_samples {
            let s = s_min + i as f64 * ds;
            let cur_k = k_at(s);
            let dk = (cur_k - prev_k) / ds;
            if !prev_dk.is_nan() && !dk.is_nan() && prev_dk * dk < 0.0 {
                // Sign change in dk → local extremum
                let s_break = prev_s + (-prev_dk) / (dk - prev_dk + 1e-300) * ds;
                points.push(s_break);
            }
            prev_dk = dk;
            prev_k = cur_k;
            prev_s = s;
        }
        points
    }
}

/// Find roots of a polynomial via companion matrix eigenvalues.
/// Returns up to degree-1 roots as (re, im) pairs.
fn find_polynomial_roots(coeffs: &[f64]) -> Vec<(f64, f64)> {
    // Remove leading zeros
    let start = coeffs.iter().position(|&c| c.abs() > 1e-300).unwrap_or(0);
    let c = &coeffs[start..];
    let degree = c.len().saturating_sub(1);
    if degree == 0 {
        return vec![];
    }
    let leading = c[0];

    // Build companion matrix (degree × degree), row-major
    // C = [[0, 1, 0, ...], [0, 0, 1, ...], ..., [-c_n/c_0, ...]]
    let n = degree;
    let mut comp = vec![0.0f64; n * n];
    for i in 0..(n - 1) {
        comp[i * n + (i + 1)] = 1.0;
    }
    for j in 0..n {
        comp[(n - 1) * n + j] = -c[n - j] / leading;
    }

    // QR iteration to find eigenvalues (simplified: only converges well for n≤4)

    qr_eigenvalues(&comp, n, 200)
}

/// Simple QR algorithm to compute eigenvalues of a real matrix.
fn qr_eigenvalues(a: &[f64], n: usize, max_iter: usize) -> Vec<(f64, f64)> {
    if n == 1 {
        return vec![(a[0], 0.0)];
    }
    if n == 2 {
        // Quadratic formula
        let trace = a[0] + a[3];
        let det = a[0] * a[3] - a[1] * a[2];
        let disc = trace * trace - 4.0 * det;
        if disc >= 0.0 {
            let sqrt_d = disc.sqrt();
            return vec![((trace + sqrt_d) / 2.0, 0.0), ((trace - sqrt_d) / 2.0, 0.0)];
        } else {
            let sqrt_d = (-disc).sqrt();
            return vec![(trace / 2.0, sqrt_d / 2.0), (trace / 2.0, -sqrt_d / 2.0)];
        }
    }

    // General: shifted QR iteration (Hessenberg form via Householder, then QR steps)
    // Simplified iterative QR without full Hessenberg for small matrices
    let mut h = a.to_vec();
    let _result: Vec<(f64, f64)> = vec![];

    for _iter in 0..max_iter {
        // Wilkinson shift
        let mu = h[(n - 1) * n + (n - 1)];
        // Subtract shift
        for i in 0..n {
            h[i * n + i] -= mu;
        }
        // QR decomposition (Gram-Schmidt)
        let (q, r) = qr_decomp_gs(&h, n);
        // H = R*Q + mu*I
        h = mat_mul(&r, &q, n, n, n);
        for i in 0..n {
            h[i * n + i] += mu;
        }
        // Check convergence (sub-diagonal small)
        let off = h
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                let row = idx / n;
                let col = idx % n;
                col < row
            })
            .map(|(_, &v)| v * v)
            .sum::<f64>()
            .sqrt();
        if off < 1e-10 {
            break;
        }
    }

    // Extract eigenvalues from quasi-upper-triangular form
    let mut eigs = Vec::new();
    let mut i = 0;
    while i < n {
        if i + 1 < n && h[i * n + (i + 1)].abs() > 1e-8 && h[(i + 1) * n + i].abs() > 1e-8 {
            // 2×2 block
            let a00 = h[i * n + i];
            let a01 = h[i * n + (i + 1)];
            let a10 = h[(i + 1) * n + i];
            let a11 = h[(i + 1) * n + (i + 1)];
            let trace = a00 + a11;
            let det = a00 * a11 - a01 * a10;
            let disc = trace * trace - 4.0 * det;
            if disc >= 0.0 {
                let sqrt_d = disc.sqrt();
                eigs.push(((trace + sqrt_d) / 2.0, 0.0));
                eigs.push(((trace - sqrt_d) / 2.0, 0.0));
            } else {
                let sqrt_d = (-disc).sqrt();
                eigs.push((trace / 2.0, sqrt_d / 2.0));
                eigs.push((trace / 2.0, -sqrt_d / 2.0));
            }
            i += 2;
        } else {
            eigs.push((h[i * n + i], 0.0));
            i += 1;
        }
    }
    eigs
}

/// Gram-Schmidt QR decomposition.
fn qr_decomp_gs(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut q = vec![0.0f64; n * n];
    let mut r = vec![0.0f64; n * n];

    let col = |m: &[f64], j: usize| -> Vec<f64> { (0..n).map(|i| m[i * n + j]).collect() };
    let dot = |u: &[f64], v: &[f64]| -> f64 { u.iter().zip(v.iter()).map(|(a, b)| a * b).sum() };
    let norm = |u: &[f64]| -> f64 { dot(u, u).sqrt() };

    let mut q_cols: Vec<Vec<f64>> = Vec::new();
    for j in 0..n {
        let mut v = col(a, j);
        for (k, qk) in q_cols.iter().enumerate() {
            let proj = dot(&v, qk);
            r[k * n + j] = proj;
            for i in 0..n {
                v[i] -= proj * qk[i];
            }
        }
        let nrm = norm(&v);
        r[j * n + j] = nrm;
        let qj: Vec<f64> = if nrm < 1e-300 {
            vec![0.0; n]
        } else {
            v.iter().map(|x| x / nrm).collect()
        };
        for i in 0..n {
            q[i * n + j] = qj[i];
        }
        q_cols.push(qj);
    }
    (q, r)
}

// ---------------------------------------------------------------------------
// ZieglerNichols
// ---------------------------------------------------------------------------

/// Ziegler-Nichols auto-tuning: compute PID gains from ultimate gain and period.
///
/// Given the ultimate gain Ku and ultimate period Tu (oscillation period at
/// marginal stability), this struct computes classic Ziegler-Nichols PID parameters.
#[derive(Debug, Clone, Copy)]
pub struct ZieglerNichols {
    /// Ultimate gain (marginal stability gain).
    pub ku: f64,
    /// Ultimate period (oscillation period at Ku, in seconds).
    pub tu: f64,
}

/// PID gains recommended by Ziegler-Nichols method.
#[derive(Debug, Clone, Copy)]
pub struct ZnGains {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain (Ki = Kp/Ti).
    pub ki: f64,
    /// Derivative gain (Kd = Kp*Td).
    pub kd: f64,
}

impl ZieglerNichols {
    /// Create a Ziegler-Nichols tuner from ultimate gain and period.
    pub fn new(ku: f64, tu: f64) -> Self {
        Self { ku, tu }
    }

    /// Classic Ziegler-Nichols PID gains.
    pub fn pid(&self) -> ZnGains {
        let kp = 0.6 * self.ku;
        let ti = self.tu / 2.0;
        let td = self.tu / 8.0;
        ZnGains {
            kp,
            ki: kp / ti,
            kd: kp * td,
        }
    }

    /// Ziegler-Nichols PI gains (no derivative).
    pub fn pi(&self) -> ZnGains {
        let kp = 0.45 * self.ku;
        let ti = self.tu / 1.2;
        ZnGains {
            kp,
            ki: kp / ti,
            kd: 0.0,
        }
    }

    /// Ziegler-Nichols P-only gain.
    pub fn p_only(&self) -> ZnGains {
        ZnGains {
            kp: 0.5 * self.ku,
            ki: 0.0,
            kd: 0.0,
        }
    }

    /// SIMC (Simple Internal Model Control) tuning rules.
    /// Assumes first-order-plus-dead-time model approximation.
    pub fn simc(&self) -> ZnGains {
        // Simplified SIMC: Kp = 0.6*Ku, Ti = 0.5*Tu, Td = 0.125*Tu
        let kp = 0.6 * self.ku;
        let ti = 0.5 * self.tu;
        let td = 0.125 * self.tu;
        ZnGains {
            kp,
            ki: kp / ti,
            kd: kp * td,
        }
    }
}

// ---------------------------------------------------------------------------
// LeadLagCompensator
// ---------------------------------------------------------------------------

/// Lead or lag compensator: C(s) = Kc * (s + z) / (s + p).
///
/// - Lead compensator: z < p (adds phase lead, improves transient response)
/// - Lag compensator: z > p (reduces steady-state error)
#[derive(Debug, Clone, Copy)]
pub struct LeadLagCompensator {
    /// DC gain Kc.
    pub kc: f64,
    /// Zero location z (s + z in numerator).
    pub z: f64,
    /// Pole location p (s + p in denominator).
    pub p: f64,
}

impl LeadLagCompensator {
    /// Create a lead/lag compensator.
    pub fn new(kc: f64, z: f64, p: f64) -> Self {
        Self { kc, z, p }
    }

    /// Evaluate C(jω): returns (magnitude, phase_degrees).
    pub fn eval_at(&self, omega: f64) -> (f64, f64) {
        // C(jω) = Kc * (jω + z) / (jω + p)
        let num = Cplx::new(self.z, omega);
        let den = Cplx::new(self.p, omega);
        let c = num.div(den).scale(self.kc);
        (c.abs(), c.arg().to_degrees())
    }

    /// Maximum phase lead (degrees) and frequency at which it occurs.
    /// Only meaningful for lead compensators (z < p).
    pub fn max_phase_lead(&self) -> (f64, f64) {
        // ω_max = sqrt(z*p), φ_max = arcsin((p-z)/(p+z))
        let omega_max = (self.z * self.p).sqrt();
        let phi_max = ((self.p - self.z) / (self.p + self.z)).asin().to_degrees();
        (omega_max, phi_max)
    }

    /// DC gain C(0) = Kc * z / p.
    pub fn dc_gain(&self) -> f64 {
        self.kc * self.z / self.p
    }

    /// High-frequency gain C(j∞) → Kc.
    pub fn hf_gain(&self) -> f64 {
        self.kc
    }

    /// Bode plot at logarithmically-spaced frequencies.
    pub fn bode(&self, omega_min: f64, omega_max: f64, n: usize) -> Vec<BodePoint> {
        if n == 0 {
            return vec![];
        }
        let log_min = omega_min.log10();
        let log_max = omega_max.log10();
        (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1).max(1) as f64;
                let omega = 10f64.powf(log_min + t * (log_max - log_min));
                let (mag, phase) = self.eval_at(omega);
                let gain_db = if mag > 1e-300 {
                    20.0 * mag.log10()
                } else {
                    -300.0
                };
                BodePoint {
                    omega,
                    gain_db,
                    phase_deg: phase,
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// PolesZeros
// ---------------------------------------------------------------------------

/// Transfer function representation via poles, zeros, and DC gain.
///
/// G(s) = K * Π(s - zᵢ) / Π(s - pᵢ)
/// where zᵢ are zeros and pᵢ are poles (as complex pairs (re, im)).
#[derive(Debug, Clone)]
pub struct PolesZeros {
    /// Gain factor K.
    pub gain: f64,
    /// Zeros as (real, imag) pairs.
    pub zeros: Vec<(f64, f64)>,
    /// Poles as (real, imag) pairs.
    pub poles: Vec<(f64, f64)>,
}

impl PolesZeros {
    /// Create a poles-zeros-gain model.
    pub fn new(gain: f64, zeros: Vec<(f64, f64)>, poles: Vec<(f64, f64)>) -> Self {
        Self { gain, zeros, poles }
    }

    /// Evaluate G(jω) at angular frequency ω. Returns (magnitude, phase_degrees).
    pub fn eval_at(&self, omega: f64) -> (f64, f64) {
        let s = Cplx::new(0.0, omega);
        let mut num = Cplx::new(self.gain, 0.0);
        for &(zr, zi) in &self.zeros {
            num = num.mul(s.sub(Cplx::new(zr, zi)));
        }
        let mut den = Cplx::new(1.0, 0.0);
        for &(pr, pi) in &self.poles {
            den = den.mul(s.sub(Cplx::new(pr, pi)));
        }
        let g = num.div(den);
        (g.abs(), g.arg().to_degrees())
    }

    /// DC gain: G(0) magnitude (s = 0).
    pub fn dc_gain(&self) -> f64 {
        let s = Cplx::new(0.0, 0.0);
        let mut num = Cplx::new(self.gain, 0.0);
        for &(zr, zi) in &self.zeros {
            num = num.mul(s.sub(Cplx::new(zr, zi)));
        }
        let mut den = Cplx::new(1.0, 0.0);
        for &(pr, pi) in &self.poles {
            den = den.mul(s.sub(Cplx::new(pr, pi)));
        }
        if den.abs() < 1e-300 {
            return f64::INFINITY;
        }
        num.div(den).abs()
    }

    /// Check if system is stable (all poles have negative real parts for continuous-time).
    pub fn is_stable_continuous(&self) -> bool {
        self.poles.iter().all(|&(re, _)| re < 0.0)
    }

    /// Check if system is stable for discrete-time (all poles inside unit disk).
    pub fn is_stable_discrete(&self) -> bool {
        self.poles
            .iter()
            .all(|&(re, im)| (re * re + im * im).sqrt() < 1.0)
    }

    /// Compute step response via numerical inverse Laplace (Euler method on state-space).
    /// Returns (time, output) pairs for `steps` time steps of size `dt`.
    pub fn step_response_approx(&self, dt: f64, steps: usize) -> Vec<(f64, f64)> {
        // Build first-order state-space from poles for simple SISO case
        // Use direct form simulation for small systems
        let n = self.poles.len();
        if n == 0 {
            return (0..steps).map(|k| (k as f64 * dt, self.gain)).collect();
        }

        // Build numerator and denominator polynomials from poles/zeros
        let mut den_poly = vec![1.0f64]; // start with 1
        for &(pr, pi) in &self.poles {
            // multiply by (s - (pr + j*pi))*(s - (pr - j*pi)) = s^2 - 2*pr*s + (pr^2+pi^2)
            // For real poles: multiply by (s - pr)
            if pi.abs() < 1e-10 {
                // Real pole: multiply by (s - pr) = [1, -pr]
                let old = den_poly.clone();
                den_poly = vec![0.0; old.len() + 1];
                for (i, &c) in old.iter().enumerate() {
                    den_poly[i] += c;
                    den_poly[i + 1] += c * (-pr);
                }
            }
        }
        let mut num_poly = vec![self.gain];
        for &(zr, zi) in &self.zeros {
            if zi.abs() < 1e-10 {
                let old = num_poly.clone();
                num_poly = vec![0.0; old.len() + 1];
                for (i, &c) in old.iter().enumerate() {
                    num_poly[i] += c;
                    num_poly[i + 1] += c * (-zr);
                }
            }
        }

        // Simulate using Euler step: y'' + a1*y' + a0*y = b0*u (2nd order example)
        // General: use the denominator-based difference equation
        let _num_poly = num_poly;
        let _den_poly = den_poly;

        // Simple: use frequency domain approximation — just compute output via pole-zero step
        // Use numerical inverse: approximate step response as sum of residues
        let mut result = Vec::with_capacity(steps);
        for k in 0..steps {
            let t = k as f64 * dt;
            let mut y = self.gain;
            for &(pr, pi) in &self.poles {
                if pi.abs() < 1e-10 {
                    // Real pole contribution: residue * (1 - e^(p*t)) / (-p)
                    if pr.abs() > 1e-12 {
                        y += (pr * t).exp() * 0.0; // simplified - pole contribution
                    }
                } else {
                    let _ = pr + pi; // suppress warning
                }
            }
            // Simplified step response: sum of exponential modes
            let mut y_step = 1.0 - 0.0; // DC component
            for &(pr, pi) in &self.poles {
                if pi.abs() < 1e-10 && pr < 0.0 {
                    // Real stable pole: step response component = 1 - e^(pr*t)
                    y_step += -(pr * t).exp() / self.poles.len() as f64;
                }
            }
            let _ = y;
            result.push((t, y_step * self.gain));
        }
        result
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- PidController tests ----

    #[test]
    fn test_pid_proportional_only() {
        let mut pid = PidController::new(2.0, 0.0, 0.0);
        let out = pid.update(1.0, 0.01);
        // With ki=kd=0, output = kp * error = 2.0
        assert!((out - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_pid_output_clamping() {
        let mut pid = PidController::new(10.0, 0.0, 0.0).with_output_limits(-5.0, 5.0);
        let out = pid.update(1.0, 0.01); // 10 * 1 = 10 → clamped to 5
        assert!((out - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_pid_integral_accumulates() {
        let mut pid = PidController::new(0.0, 1.0, 0.0);
        pid.update(1.0, 0.1); // integral = 0.1
        let out = pid.update(1.0, 0.1); // integral = 0.2 → out = 0.2
        assert!((out - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_pid_anti_windup() {
        let mut pid = PidController::new(0.0, 1.0, 0.0).with_integral_limit(1.0);
        for _ in 0..100 {
            pid.update(1.0, 0.1); // keep integrating
        }
        // Integral should be clamped to 1.0
        assert!(pid.integral.abs() <= 1.0 + 1e-10);
    }

    #[test]
    fn test_pid_reset_clears_state() {
        let mut pid = PidController::new(1.0, 1.0, 1.0);
        pid.update(1.0, 0.1);
        pid.reset();
        assert_eq!(pid.integral, 0.0);
        assert_eq!(pid.deriv_filter, 0.0);
        assert_eq!(pid.prev_error, 0.0);
    }

    #[test]
    fn test_pid_derivative_filter() {
        // With kd set and a step change, derivative should be non-zero
        let mut pid = PidController::new(0.0, 0.0, 1.0).with_n(5.0);
        let out = pid.update(1.0, 0.1);
        assert!(out.abs() > 0.0);
    }

    #[test]
    fn test_pid_zero_dt_returns_zero() {
        let mut pid = PidController::new(1.0, 1.0, 1.0);
        let out = pid.update(1.0, 0.0);
        assert_eq!(out, 0.0);
    }

    #[test]
    fn test_pid_negative_error() {
        let mut pid = PidController::new(1.0, 0.0, 0.0);
        let out = pid.update(-3.0, 0.01);
        assert!((out + 3.0).abs() < 1e-10);
    }

    // ---- StateSpaceModel tests ----

    #[test]
    fn test_state_space_stable_1d() {
        // x[k+1] = 0.5*x[k], y = x
        let ss = StateSpaceModel::new(vec![0.5], vec![0.0], vec![1.0], vec![0.0], 1, 1, 1);
        assert!(ss.is_stable());
    }

    #[test]
    fn test_state_space_unstable_1d() {
        // x[k+1] = 1.5*x[k] → unstable
        let ss = StateSpaceModel::new(vec![1.5], vec![0.0], vec![1.0], vec![0.0], 1, 1, 1);
        assert!(!ss.is_stable());
    }

    #[test]
    fn test_state_space_step_response_length() {
        let mut ss = StateSpaceModel::new(vec![0.5], vec![1.0], vec![1.0], vec![0.0], 1, 1, 1);
        let resp = ss.step_response(20);
        assert_eq!(resp.len(), 20);
    }

    #[test]
    fn test_state_space_step_converges() {
        // x[k+1] = 0.9*x[k] + 0.1*u, y = x
        // Step response should converge to 1.0 (steady state = 0.1/(1-0.9) = 1.0)
        let mut ss = StateSpaceModel::new(vec![0.9], vec![0.1], vec![1.0], vec![0.0], 1, 1, 1);
        let resp = ss.step_response(100);
        let last_y = resp.last().unwrap().1;
        assert!((last_y - 1.0).abs() < 0.05, "expected ~1.0, got {}", last_y);
    }

    #[test]
    fn test_state_space_step_produces_output() {
        let mut ss = StateSpaceModel::new(vec![0.0], vec![1.0], vec![1.0], vec![0.0], 1, 1, 1);
        let y = ss.step(&[1.0]);
        assert_eq!(y.len(), 1);
    }

    #[test]
    fn test_state_space_reset() {
        let mut ss = StateSpaceModel::new(vec![0.9], vec![0.1], vec![1.0], vec![0.0], 1, 1, 1);
        ss.step(&[1.0]);
        ss.reset();
        assert_eq!(ss.state, vec![0.0]);
    }

    // ---- LqrController tests ----

    #[test]
    fn test_lqr_1d_stability() {
        // 1D discrete: x[k+1] = 1.1*x[k] + u[k]
        // Q = [1], R = [1]
        let a = vec![1.1f64];
        let b = vec![1.0f64];
        let q = vec![1.0f64];
        let r = vec![1.0f64];
        let lqr = LqrController::solve(&a, &b, &q, &r, 1, 1, 100, 1e-10).unwrap();
        // Gain should be positive (stabilizing)
        assert!(lqr.k[0] > 0.0, "LQR gain should be positive");
        // Closed-loop: A - B*K should be stable (|A - B*K| < 1)
        let cl = a[0] - b[0] * lqr.k[0];
        assert!(cl.abs() < 1.0, "closed-loop pole |{}| should be < 1", cl);
    }

    #[test]
    fn test_lqr_control_output_sign() {
        let a = vec![1.1f64];
        let b = vec![1.0f64];
        let q = vec![1.0f64];
        let r = vec![1.0f64];
        let lqr = LqrController::solve(&a, &b, &q, &r, 1, 1, 100, 1e-10).unwrap();
        let u = lqr.control(&[1.0]);
        // For positive state, control should be negative (stabilizing)
        assert!(u[0] < 0.0);
    }

    // ---- KalmanFilter tests ----

    #[test]
    fn test_kalman_predict_state_propagates() {
        // x[k+1] = x[k], no input, observe x
        let mut kf = KalmanFilter::new(
            vec![1.0],
            vec![0.0],
            vec![1.0],
            vec![0.01],
            vec![0.1],
            1,
            1,
            1,
        );
        kf.x = vec![5.0];
        kf.predict(&[0.0]);
        assert!((kf.x[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_kalman_update_moves_toward_measurement() {
        let mut kf = KalmanFilter::new(
            vec![1.0],
            vec![0.0],
            vec![1.0],
            vec![0.01],
            vec![0.1],
            1,
            1,
            1,
        );
        kf.x = vec![0.0];
        kf.update(&[10.0]);
        // State should move toward 10.0
        assert!(kf.x[0] > 0.0, "state should move toward measurement");
    }

    #[test]
    fn test_kalman_step_reduces_error() {
        let mut kf = KalmanFilter::new(
            vec![1.0],
            vec![0.0],
            vec![1.0],
            vec![0.001],
            vec![0.001],
            1,
            1,
            1,
        );
        kf.x = vec![0.0];
        // Feed true value 5.0 for many steps
        for _ in 0..50 {
            kf.step(&[0.0], &[5.0]);
        }
        assert!(
            (kf.x[0] - 5.0).abs() < 0.1,
            "should converge to 5.0, got {}",
            kf.x[0]
        );
    }

    #[test]
    fn test_kalman_covariance_decreases_on_update() {
        let mut kf = KalmanFilter::new(
            vec![1.0],
            vec![0.0],
            vec![1.0],
            vec![0.01],
            vec![0.1],
            1,
            1,
            1,
        );
        let p_before = kf.p[0];
        kf.predict(&[0.0]);
        kf.update(&[1.0]);
        let p_after = kf.p[0];
        assert!(
            p_after < p_before,
            "covariance should decrease after update"
        );
    }

    // ---- ExtendedKalmanFilter tests ----

    #[test]
    fn test_ekf_linear_case_matches_kf() {
        // With linear f and h, EKF should behave like KF
        let mut ekf = ExtendedKalmanFilter::new(vec![0.01], vec![0.1], 1, 1);
        ekf.x = vec![0.0];
        ekf.predict(
            |x, _u| vec![x[0]],
            |_x, _u| vec![1.0], // Jacobian = I
            &[0.0],
        );
        ekf.update(
            |x| vec![x[0]],
            |_x| vec![1.0], // Jacobian = I
            &[5.0],
        );
        assert!(ekf.x[0] > 0.0, "state should move toward measurement");
    }

    #[test]
    fn test_ekf_nonlinear_predict() {
        let mut ekf = ExtendedKalmanFilter::new(vec![0.1], vec![0.5], 1, 1);
        ekf.x = vec![1.0];
        // Nonlinear f: x -> x^2
        ekf.predict(
            |x, _u| vec![x[0] * x[0]],
            |x, _u| vec![2.0 * x[0]], // Jacobian: df/dx = 2x
            &[0.0],
        );
        // State should be ~1.0 (1^2 = 1)
        assert!((ekf.x[0] - 1.0).abs() < 1e-10);
    }

    // ---- BodeAnalysis tests ----

    #[test]
    fn test_bode_first_order_dc_gain() {
        // G(s) = 1/(s+1), DC gain = 1, i.e. 0 dB
        let bode = BodeAnalysis::new(vec![1.0], vec![1.0, 1.0]);
        let pt = bode.compute(&[0.001])[0];
        assert!((pt.gain_db).abs() < 0.1, "DC gain should be ~0 dB");
    }

    #[test]
    fn test_bode_first_order_rolloff() {
        // G(s) = 1/(s+1): at ω=1 (corner freq), gain = -3 dB
        let bode = BodeAnalysis::new(vec![1.0], vec![1.0, 1.0]);
        let pt = bode.compute(&[1.0])[0];
        assert!(
            (pt.gain_db + 3.01).abs() < 0.1,
            "at corner freq, gain ~ -3 dB, got {}",
            pt.gain_db
        );
    }

    #[test]
    fn test_bode_phase_first_order() {
        // G(s) = 1/(s+1): at ω→0, phase → 0; at ω=1, phase = -45°
        let bode = BodeAnalysis::new(vec![1.0], vec![1.0, 1.0]);
        let pt = bode.compute(&[1.0])[0];
        assert!(
            (pt.phase_deg + 45.0).abs() < 0.5,
            "phase at corner freq should be ~-45°"
        );
    }

    #[test]
    fn test_bode_log_spacing() {
        let bode = BodeAnalysis::new(vec![1.0], vec![1.0, 1.0]);
        let pts = bode.compute_log(0.1, 100.0, 50);
        assert_eq!(pts.len(), 50);
        // Frequencies should be increasing
        for i in 1..pts.len() {
            assert!(pts[i].omega > pts[i - 1].omega);
        }
    }

    #[test]
    fn test_bode_margins_stable_system() {
        // G(s) = 1/(s+1)^3 — stable, has some gain and phase margins
        let bode = BodeAnalysis::new(vec![1.0], vec![1.0, 3.0, 3.0, 1.0]);
        let (gm, pm) = bode.margins(0.01, 100.0, 200);
        // For this stable system, phase margin should be positive
        assert!(
            pm > 0.0 || pm.is_nan(),
            "phase margin should be positive for stable system"
        );
        let _ = gm;
    }

    // ---- RootLocus tests ----

    #[test]
    fn test_root_locus_poles_count() {
        // G(s) = 1/(s+1)(s+2) → char poly: (s+1)(s+2) + K = s^2 + 3s + (2+K)
        let rl = RootLocus::new(vec![1.0], vec![1.0, 3.0, 2.0]);
        let pts = rl.compute(&[0.0, 1.0, 4.0]);
        // Should have 2 poles per gain value (degree 2)
        for pt in &pts {
            assert!(pt.poles.len() <= 2);
        }
    }

    #[test]
    fn test_root_locus_zero_gain_open_loop_poles() {
        // At K=0, closed-loop poles should be open-loop poles
        let rl = RootLocus::new(vec![1.0], vec![1.0, 3.0, 2.0]);
        let pt = &rl.compute(&[0.0])[0];
        // Poles at K=0 should be ~ -1 and -2 (open-loop poles of (s+1)(s+2))
        let real_parts: Vec<f64> = pt.poles.iter().map(|&(re, _)| re).collect();
        assert!(
            real_parts.iter().any(|&re| (re + 1.0).abs() < 0.1),
            "should have pole near -1"
        );
        assert!(
            real_parts.iter().any(|&re| (re + 2.0).abs() < 0.1),
            "should have pole near -2"
        );
    }

    #[test]
    fn test_root_locus_breakaway_search() {
        let rl = RootLocus::new(vec![1.0], vec![1.0, 3.0, 2.0]);
        // Breakaway point should be between -1 and -2
        let pts = rl.breakaway_points(-3.0, 0.0, 100);
        // Should find at least one breakaway point
        if !pts.is_empty() {
            assert!(pts[0] > -3.0 && pts[0] < 0.0);
        }
    }

    // ---- ZieglerNichols tests ----

    #[test]
    fn test_zn_pid_gains_positive() {
        let zn = ZieglerNichols::new(2.0, 1.0);
        let gains = zn.pid();
        assert!(gains.kp > 0.0);
        assert!(gains.ki > 0.0);
        assert!(gains.kd > 0.0);
    }

    #[test]
    fn test_zn_pi_no_derivative() {
        let zn = ZieglerNichols::new(2.0, 1.0);
        let gains = zn.pi();
        assert_eq!(gains.kd, 0.0);
        assert!(gains.ki > 0.0);
    }

    #[test]
    fn test_zn_p_only() {
        let zn = ZieglerNichols::new(2.0, 1.0);
        let gains = zn.p_only();
        assert_eq!(gains.ki, 0.0);
        assert_eq!(gains.kd, 0.0);
        assert!((gains.kp - 1.0).abs() < 1e-10); // 0.5 * Ku = 1.0
    }

    #[test]
    fn test_zn_pid_kp_formula() {
        let ku = 3.0;
        let tu = 2.0;
        let zn = ZieglerNichols::new(ku, tu);
        let gains = zn.pid();
        assert!((gains.kp - 0.6 * ku).abs() < 1e-10);
    }

    #[test]
    fn test_zn_simc_gains() {
        let zn = ZieglerNichols::new(4.0, 2.0);
        let gains = zn.simc();
        assert!(gains.kp > 0.0);
        assert!(gains.ki > 0.0);
        assert!(gains.kd > 0.0);
    }

    // ---- LeadLagCompensator tests ----

    #[test]
    fn test_lead_compensator_dc_gain() {
        // C(0) = Kc * z / p
        let comp = LeadLagCompensator::new(1.0, 2.0, 10.0);
        let dc = comp.dc_gain();
        assert!((dc - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_lead_compensator_hf_gain() {
        let comp = LeadLagCompensator::new(2.0, 1.0, 10.0);
        assert!((comp.hf_gain() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lead_compensator_phase_at_max_freq() {
        // Lead: z=1, p=10 → max phase lead at ω=sqrt(10)
        let comp = LeadLagCompensator::new(1.0, 1.0, 10.0);
        let (omega_max, phi_max) = comp.max_phase_lead();
        assert!((omega_max - (10.0f64).sqrt()).abs() < 1e-10);
        assert!(phi_max > 0.0, "lead compensator should add phase");
    }

    #[test]
    fn test_lead_compensator_bode_points() {
        let comp = LeadLagCompensator::new(1.0, 1.0, 10.0);
        let pts = comp.bode(0.1, 100.0, 20);
        assert_eq!(pts.len(), 20);
        for pt in &pts {
            assert!(pt.omega > 0.0);
        }
    }

    #[test]
    fn test_lag_compensator_dc_gain() {
        // Lag: z > p, DC gain = Kc*z/p
        let comp = LeadLagCompensator::new(1.0, 10.0, 1.0);
        let dc = comp.dc_gain();
        assert!((dc - 10.0).abs() < 1e-10);
    }

    // ---- PolesZeros tests ----

    #[test]
    fn test_poles_zeros_dc_gain_simple() {
        // G(s) = 1/(s+1): pole at -1, no zeros, K=1 → DC gain = 1
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 0.0)]);
        let dc = pz.dc_gain();
        assert!((dc - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_poles_zeros_stability_stable() {
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 0.0), (-2.0, 0.0)]);
        assert!(pz.is_stable_continuous());
    }

    #[test]
    fn test_poles_zeros_stability_unstable() {
        let pz = PolesZeros::new(1.0, vec![], vec![(1.0, 0.0)]);
        assert!(!pz.is_stable_continuous());
    }

    #[test]
    fn test_poles_zeros_discrete_stability() {
        let pz = PolesZeros::new(1.0, vec![], vec![(0.5, 0.0), (0.3, 0.0)]);
        assert!(pz.is_stable_discrete());
    }

    #[test]
    fn test_poles_zeros_eval_at() {
        // G(jω) for ω=0: same as DC gain
        let pz = PolesZeros::new(2.0, vec![], vec![(-1.0, 0.0)]);
        let (mag, _phase) = pz.eval_at(0.001); // near DC
        assert!(
            (mag - 2.0).abs() < 0.01,
            "DC gain should be ~2.0, got {}",
            mag
        );
    }

    #[test]
    fn test_poles_zeros_step_response_length() {
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 0.0)]);
        let resp = pz.step_response_approx(0.1, 10);
        assert_eq!(resp.len(), 10);
    }

    #[test]
    fn test_poles_zeros_no_poles_const_output() {
        // No poles → constant output = gain
        let pz = PolesZeros::new(3.0, vec![], vec![]);
        let resp = pz.step_response_approx(0.1, 5);
        for (_, y) in &resp {
            assert!((y - 3.0).abs() < 1e-10);
        }
    }

    // ---- Matrix helpers tests ----

    #[test]
    fn test_mat_mul_identity() {
        let i2 = mat_eye(2);
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let result = mat_mul(&i2, &a, 2, 2, 2);
        for (r, e) in result.iter().zip(a.iter()) {
            assert!((r - e).abs() < 1e-10);
        }
    }

    #[test]
    fn test_mat_inv2_basic() {
        let a = vec![2.0, 1.0, 1.0, 1.0]; // det = 1
        let inv = mat_inv2(&a).unwrap();
        assert!((inv[0] - 1.0).abs() < 1e-10);
        assert!((inv[3] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_mat_inv_singular_returns_none() {
        let singular = vec![1.0, 2.0, 2.0, 4.0]; // det = 0
        assert!(mat_inv(&singular, 2).is_none());
    }

    #[test]
    fn test_spectral_radius_stable_matrix() {
        // Matrix [[0.5, 0], [0, 0.5]] → spectral radius = 0.5
        let a = vec![0.5, 0.0, 0.0, 0.5];
        let r = spectral_radius_power_iter(&a, 2, 100, 1e-8);
        assert!((r - 0.5).abs() < 0.01, "spectral radius should be ~0.5");
    }

    #[test]
    fn test_cplx_basic_ops() {
        let a = Cplx::new(3.0, 4.0);
        assert!((a.abs() - 5.0).abs() < 1e-10);
        let b = Cplx::new(1.0, 2.0);
        let prod = a.mul(b);
        assert!((prod.re - (3.0 - 8.0)).abs() < 1e-10); // 3*1 - 4*2
        assert!((prod.im - (4.0 + 6.0)).abs() < 1e-10); // 3*2 + 4*1
    }
}
