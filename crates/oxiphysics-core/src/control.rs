#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Control theory and systems: PID, state-space, transfer functions, Bode analysis,
//! LQR, Kalman filter, MPC, adaptive control, fuzzy control, and neural controller.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Helper math utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Matrix-vector multiply: C = A * b, where A is n×n and b is length n.
fn mat_vec_mul(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mut out = vec![0.0; n];
    for i in 0..n {
        for j in 0..b.len() {
            out[i] += a[i][j] * b[j];
        }
    }
    out
}

/// Matrix-matrix multiply: C = A * B.
fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let rows = a.len();
    let cols = b[0].len();
    let inner = b.len();
    let mut c = vec![vec![0.0; cols]; rows];
    for i in 0..rows {
        for j in 0..cols {
            for k in 0..inner {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Transpose of a matrix.
fn mat_transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return vec![];
    }
    let rows = a.len();
    let cols = a[0].len();
    let mut t = vec![vec![0.0; rows]; cols];
    for i in 0..rows {
        for j in 0..cols {
            t[j][i] = a[i][j];
        }
    }
    t
}

/// Add two vectors element-wise.
fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Scale a vector.
fn vec_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}

/// Dot product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Identity matrix n×n.
fn eye(n: usize) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        m[i][i] = 1.0;
    }
    m
}

/// Add two matrices element-wise.
fn mat_add(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let rows = a.len();
    let cols = a[0].len();
    let mut c = vec![vec![0.0; cols]; rows];
    for i in 0..rows {
        for j in 0..cols {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Scale a matrix.
fn mat_scale(a: &[Vec<f64>], s: f64) -> Vec<Vec<f64>> {
    a.iter()
        .map(|row| row.iter().map(|x| x * s).collect())
        .collect()
}

/// Simple matrix inverse using Gauss-Jordan elimination (n×n).
fn mat_inv(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut aug: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            for j in 0..n {
                r.push(if i == j { 1.0 } else { 0.0 });
            }
            r
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| {
            aug[i][col]
                .abs()
                .partial_cmp(&aug[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        aug.swap(col, pivot);
        let diag = aug[col][col];
        if diag.abs() < 1e-14 {
            return None;
        }
        for j in 0..(2 * n) {
            aug[col][j] /= diag;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row][col];
                for j in 0..(2 * n) {
                    let val = aug[col][j] * factor;
                    aug[row][j] -= val;
                }
            }
        }
    }
    Some(aug.iter().map(|row| row[n..].to_vec()).collect())
}

// ─────────────────────────────────────────────────────────────────────────────
// PID Controller
// ─────────────────────────────────────────────────────────────────────────────

/// Proportional-integral-derivative (PID) controller with anti-windup and
/// derivative filter.
#[derive(Debug, Clone)]
pub struct Pid {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Integral accumulator.
    pub integral: f64,
    /// Previous error for derivative computation.
    pub prev_error: f64,
    /// Filtered derivative state.
    pub derivative_filter_state: f64,
    /// Derivative filter coefficient N (cutoff factor).
    pub filter_coeff: f64,
    /// Anti-windup: maximum integral magnitude.
    pub integral_limit: f64,
    /// Output saturation maximum.
    pub output_max: f64,
    /// Output saturation minimum.
    pub output_min: f64,
}

impl Pid {
    /// Create a new PID controller with given gains.
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            derivative_filter_state: 0.0,
            filter_coeff: 10.0,
            integral_limit: f64::MAX,
            output_max: f64::MAX,
            output_min: f64::MIN,
        }
    }

    /// Set gains.
    pub fn set_gains(&mut self, kp: f64, ki: f64, kd: f64) {
        self.kp = kp;
        self.ki = ki;
        self.kd = kd;
    }

    /// Set anti-windup integral limit.
    pub fn set_integral_limit(&mut self, limit: f64) {
        self.integral_limit = limit.abs();
    }

    /// Set output saturation limits.
    pub fn set_output_limits(&mut self, min: f64, max: f64) {
        self.output_min = min;
        self.output_max = max;
    }

    /// Compute PID output given current `error` and time step `dt`.
    pub fn compute(&mut self, error: f64, dt: f64) -> f64 {
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

        // Derivative with low-pass filter (Tustin bilinear)
        // Filtered derivative: d_f = N*(e - x_f), x_f += dt * d_f
        let n = self.filter_coeff;
        let d_f = n * (error - self.derivative_filter_state);
        self.derivative_filter_state += dt * d_f;
        let d = self.kd * d_f;

        self.prev_error = error;

        let output = p + i + d;
        output.clamp(self.output_min, self.output_max)
    }

    /// Reset internal state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
        self.derivative_filter_state = 0.0;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// State-Space Model
// ─────────────────────────────────────────────────────────────────────────────

/// Continuous or discrete linear state-space model: ẋ = Ax + Bu, y = Cx + Du.
#[derive(Debug, Clone)]
pub struct StateSpaceModel {
    /// State matrix (n×n).
    pub a: Vec<Vec<f64>>,
    /// Input matrix (n×m).
    pub b: Vec<Vec<f64>>,
    /// Output matrix (p×n).
    pub c: Vec<Vec<f64>>,
    /// Feedthrough matrix (p×m).
    pub d: Vec<Vec<f64>>,
    /// Current state vector (length n).
    pub state: Vec<f64>,
    /// Whether this model is discrete-time.
    pub is_discrete: bool,
}

impl StateSpaceModel {
    /// Create a new state-space model.
    pub fn new(a: Vec<Vec<f64>>, b: Vec<Vec<f64>>, c: Vec<Vec<f64>>, d: Vec<Vec<f64>>) -> Self {
        let n = a.len();
        Self {
            a,
            b,
            c,
            d,
            state: vec![0.0; n],
            is_discrete: false,
        }
    }

    /// Set initial state.
    pub fn set_state(&mut self, x0: Vec<f64>) {
        self.state = x0;
    }

    /// Step the continuous system using Euler integration; returns output y.
    pub fn step(&mut self, u: &[f64], dt: f64) -> Vec<f64> {
        if self.is_discrete {
            // x[k+1] = A*x[k] + B*u[k]
            let ax = mat_vec_mul(&self.a, &self.state);
            let bu = mat_vec_mul(&self.b, u);
            self.state = vec_add(&ax, &bu);
        } else {
            // Euler: x += dt*(A*x + B*u)
            let ax = mat_vec_mul(&self.a, &self.state);
            let bu = mat_vec_mul(&self.b, u);
            let dx = vec_add(&ax, &bu);
            let dx_dt = vec_scale(&dx, dt);
            self.state = vec_add(&self.state, &dx_dt);
        }
        // y = C*x + D*u
        let cx = mat_vec_mul(&self.c, &self.state);
        let du = mat_vec_mul(&self.d, u);
        vec_add(&cx, &du)
    }

    /// Discretize using zero-order hold (ZOH) approximation.
    pub fn discretize_zoh(&self, dt: f64) -> Self {
        let ad = zoh_discretize(&self.a, dt);
        // B_d ≈ dt * B (first-order ZOH approximation)
        let bd = mat_scale(&self.b, dt);
        let mut disc = Self::new(ad, bd, self.c.clone(), self.d.clone());
        disc.is_discrete = true;
        disc
    }

    /// Number of states.
    pub fn n_states(&self) -> usize {
        self.a.len()
    }

    /// Number of inputs.
    pub fn n_inputs(&self) -> usize {
        if self.b.is_empty() {
            0
        } else {
            self.b[0].len()
        }
    }

    /// Number of outputs.
    pub fn n_outputs(&self) -> usize {
        self.c.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Transfer Function
// ─────────────────────────────────────────────────────────────────────────────

/// SISO transfer function defined by numerator and denominator polynomials
/// (descending powers).
#[derive(Debug, Clone)]
pub struct TransferFunction {
    /// Numerator polynomial coefficients (highest power first).
    pub num: Vec<f64>,
    /// Denominator polynomial coefficients (highest power first).
    pub den: Vec<f64>,
}

impl TransferFunction {
    /// Create a new transfer function.
    pub fn new(num: Vec<f64>, den: Vec<f64>) -> Self {
        Self { num, den }
    }

    /// Create from PID parameters: H(s) = kd*s^2 + kp*s + ki / s.
    pub fn from_pid(kp: f64, ki: f64, kd: f64) -> Self {
        // C(s) = (kd*s^2 + kp*s + ki) / s
        Self::new(vec![kd, kp, ki], vec![1.0, 0.0])
    }

    /// Evaluate the polynomial at complex point s = (re, im).
    fn eval_poly_complex(coeffs: &[f64], re: f64, im: f64) -> (f64, f64) {
        let mut re_acc = 0.0_f64;
        let mut im_acc = 0.0_f64;
        for &c in coeffs {
            // multiply by s = re + j*im and add c
            let new_re = re_acc * re - im_acc * im + c;
            let new_im = re_acc * im + im_acc * re;
            re_acc = new_re;
            im_acc = new_im;
        }
        (re_acc, im_acc)
    }

    /// Evaluate H(jω) and return (magnitude, phase_deg).
    pub fn frequency_response(&self, omega: f64) -> (f64, f64) {
        let (nr, ni) = Self::eval_poly_complex(&self.num, 0.0, omega);
        let (dr, di) = Self::eval_poly_complex(&self.den, 0.0, omega);
        let num_mag = (nr * nr + ni * ni).sqrt();
        let den_mag = (dr * dr + di * di).sqrt();
        if den_mag < 1e-30 {
            return (f64::INFINITY, 0.0);
        }
        let mag = num_mag / den_mag;
        let phase = ni.atan2(nr) - di.atan2(dr);
        (mag, phase.to_degrees())
    }

    /// DC gain (ω → 0), i.e., H(0).
    pub fn dc_gain(&self) -> f64 {
        let n_val = self.num.last().copied().unwrap_or(0.0);
        let d_val = self.den.last().copied().unwrap_or(1.0);
        if d_val.abs() < 1e-30 {
            f64::INFINITY
        } else {
            n_val / d_val
        }
    }

    /// Compute poles (roots of denominator) using companion matrix eigenvalues
    /// (only real poles for now via Durand-Kerner).
    pub fn poles(&self) -> Vec<f64> {
        poly_roots_real(&self.den)
    }

    /// Compute zeros (roots of numerator).
    pub fn zeros(&self) -> Vec<f64> {
        poly_roots_real(&self.num)
    }
}

/// Find real roots of a polynomial with descending-power coefficients using
/// numerical deflation (simple Newton-Raphson approximation).
fn poly_roots_real(coeffs: &[f64]) -> Vec<f64> {
    let n = coeffs.len();
    if n <= 1 {
        return vec![];
    }
    // Normalize
    let leading = coeffs[0];
    if leading.abs() < 1e-30 {
        return poly_roots_real(&coeffs[1..]);
    }
    let p: Vec<f64> = coeffs.iter().map(|&c| c / leading).collect();
    let degree = p.len() - 1;
    if degree == 0 {
        return vec![];
    }
    if degree == 1 {
        return vec![-p[1]];
    }
    if degree == 2 {
        let disc = p[1] * p[1] - 4.0 * p[2];
        if disc >= 0.0 {
            return vec![(-p[1] + disc.sqrt()) / 2.0, (-p[1] - disc.sqrt()) / 2.0];
        }
        return vec![];
    }
    // For higher degree, use Newton on the polynomial
    let mut roots = Vec::new();
    let mut remaining = p.clone();
    while remaining.len() > 2 {
        let mut x = 0.1_f64;
        for _ in 0..200 {
            let (fx, dfx) = eval_poly_and_deriv(&remaining, x);
            if dfx.abs() < 1e-30 {
                x += 0.01;
                continue;
            }
            let dx = fx / dfx;
            x -= dx;
            if dx.abs() < 1e-12 {
                break;
            }
        }
        roots.push(x);
        // deflate
        remaining = poly_deflate(&remaining, x);
    }
    if remaining.len() == 2 {
        roots.push(-remaining[1] / remaining[0]);
    }
    roots
}

fn eval_poly_and_deriv(p: &[f64], x: f64) -> (f64, f64) {
    let mut f = p[0];
    let mut df = 0.0_f64;
    for &c in &p[1..] {
        df = df * x + f;
        f = f * x + c;
    }
    (f, df)
}

fn poly_deflate(p: &[f64], root: f64) -> Vec<f64> {
    let n = p.len();
    let mut q = vec![0.0; n - 1];
    q[0] = p[0];
    for i in 1..n - 1 {
        q[i] = p[i] + q[i - 1] * root;
    }
    q
}

// ─────────────────────────────────────────────────────────────────────────────
// Bode Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Bode frequency response analysis for a SISO transfer function.
#[derive(Debug, Clone)]
pub struct BodeAnalysis {
    /// The transfer function being analyzed.
    pub tf: TransferFunction,
}

impl BodeAnalysis {
    /// Create from a transfer function.
    pub fn new(tf: TransferFunction) -> Self {
        Self { tf }
    }

    /// Compute magnitude in dB at angular frequency ω.
    pub fn magnitude_db(&self, omega: f64) -> f64 {
        let (mag, _) = self.tf.frequency_response(omega);
        20.0 * mag.log10()
    }

    /// Compute phase in degrees at angular frequency ω.
    pub fn phase_deg(&self, omega: f64) -> f64 {
        let (_, phase) = self.tf.frequency_response(omega);
        phase
    }

    /// Sweep frequencies \[ω_min, ω_max\] with `n` points (log-spaced).
    /// Returns `(omegas, magnitudes_db, phases_deg)`.
    pub fn sweep(
        &self,
        omega_min: f64,
        omega_max: f64,
        n: usize,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let log_min = omega_min.ln();
        let log_max = omega_max.ln();
        let omegas: Vec<f64> = (0..n)
            .map(|i| (log_min + (log_max - log_min) * i as f64 / (n - 1) as f64).exp())
            .collect();
        let mags: Vec<f64> = omegas.iter().map(|&w| self.magnitude_db(w)).collect();
        let phases: Vec<f64> = omegas.iter().map(|&w| self.phase_deg(w)).collect();
        (omegas, mags, phases)
    }

    /// Estimate gain margin and phase margin.
    /// Returns `(gain_margin_db, phase_crossover_freq, phase_margin_deg, gain_crossover_freq)`.
    pub fn stability_margins(
        &self,
        omega_min: f64,
        omega_max: f64,
        n: usize,
    ) -> (f64, f64, f64, f64) {
        let (omegas, mags, phases) = self.sweep(omega_min, omega_max, n);

        // Phase crossover: where phase ≈ -180°
        let mut phase_crossover_freq = 0.0;
        let mut gain_at_phase_cross = 0.0;
        for i in 1..n {
            if phases[i - 1] > -180.0 && phases[i] <= -180.0 {
                let t = (-180.0 - phases[i - 1]) / (phases[i] - phases[i - 1]);
                phase_crossover_freq = omegas[i - 1] + t * (omegas[i] - omegas[i - 1]);
                gain_at_phase_cross = mags[i - 1] + t * (mags[i] - mags[i - 1]);
                break;
            }
        }
        let gain_margin = -gain_at_phase_cross;

        // Gain crossover: where magnitude ≈ 0 dB
        let mut gain_crossover_freq = 0.0;
        let mut phase_at_gain_cross = 0.0;
        for i in 1..n {
            if mags[i - 1] > 0.0 && mags[i] <= 0.0 {
                let t = (0.0 - mags[i - 1]) / (mags[i] - mags[i - 1]);
                gain_crossover_freq = omegas[i - 1] + t * (omegas[i] - omegas[i - 1]);
                phase_at_gain_cross = phases[i - 1] + t * (phases[i] - phases[i - 1]);
                break;
            }
        }
        let phase_margin = 180.0 + phase_at_gain_cross;

        (
            gain_margin,
            phase_crossover_freq,
            phase_margin,
            gain_crossover_freq,
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LQR Controller
// ─────────────────────────────────────────────────────────────────────────────

/// Linear Quadratic Regulator (LQR): computes optimal gain matrix K via
/// Discrete Algebraic Riccati Equation (DARE).
#[derive(Debug, Clone)]
pub struct LqrController {
    /// Gain matrix K (m×n).
    pub k: Vec<Vec<f64>>,
    /// State weighting matrix Q (n×n).
    pub q: Vec<Vec<f64>>,
    /// Input weighting matrix R (m×m).
    pub r: Vec<Vec<f64>>,
}

impl LqrController {
    /// Solve DARE and compute LQR gain K for discrete system (A, B).
    /// Returns `None` if the Riccati equation does not converge.
    pub fn solve(
        a: &[Vec<f64>],
        b: &[Vec<f64>],
        q: Vec<Vec<f64>>,
        r: Vec<Vec<f64>>,
    ) -> Option<Self> {
        let k = riccati_solve_dare(a, b, &q, &r, 1000, 1e-10)?;
        Some(Self { k, q, r })
    }

    /// Compute optimal control: u = -K * x.
    pub fn compute(&self, x: &[f64]) -> Vec<f64> {
        let kx = mat_vec_mul(&self.k, x);
        vec_scale(&kx, -1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Kalman Filter
// ─────────────────────────────────────────────────────────────────────────────

/// Full linear Kalman filter with predict and update steps.
#[derive(Debug, Clone)]
pub struct KalmanFilter {
    /// State estimate (length n).
    pub x_hat: Vec<f64>,
    /// Error covariance (n×n).
    pub p: Vec<Vec<f64>>,
}

impl KalmanFilter {
    /// Create with initial state estimate and covariance.
    pub fn new(x0: Vec<f64>, p0: Vec<Vec<f64>>) -> Self {
        Self { x_hat: x0, p: p0 }
    }

    /// Prediction step: x̂⁻ = A*x̂ + B*u, P⁻ = A*P*Aᵀ + Q.
    pub fn predict(&mut self, a: &[Vec<f64>], b: &[Vec<f64>], u: &[f64], q: &[Vec<f64>]) {
        let ax = mat_vec_mul(a, &self.x_hat);
        let bu = mat_vec_mul(b, u);
        self.x_hat = vec_add(&ax, &bu);

        let at = mat_transpose(a);
        let ap = mat_mul(a, &self.p);
        let apat = mat_mul(&ap, &at);
        self.p = mat_add(&apat, q);
    }

    /// Update step: K = P⁻Hᵀ(H*P⁻*Hᵀ+R)⁻¹, x̂ = x̂⁻ + K*(z - H*x̂⁻), P = (I-KH)*P⁻.
    pub fn update(&mut self, h: &[Vec<f64>], r: &[Vec<f64>], z: &[f64]) -> Option<()> {
        let n = self.x_hat.len();
        let ht = mat_transpose(h);
        let ph_t = mat_mul(&self.p, &ht);
        let hph_t = mat_mul(h, &ph_t);
        let s = mat_add(&hph_t, r);
        let s_inv = mat_inv(&s)?;
        // Kalman gain K = P*Hᵀ * S⁻¹
        let k = mat_mul(&ph_t, &s_inv);
        // Innovation: z - H*x̂
        let hx = mat_vec_mul(h, &self.x_hat);
        let innov: Vec<f64> = z.iter().zip(hx.iter()).map(|(zi, hxi)| zi - hxi).collect();
        // State update
        let k_innov = mat_vec_mul(&k, &innov);
        self.x_hat = vec_add(&self.x_hat, &k_innov);
        // Covariance update: P = (I - K*H) * P
        let kh = mat_mul(&k, h);
        let i = eye(n);
        let i_minus_kh = mat_add(&i, &mat_scale(&kh, -1.0));
        self.p = mat_mul(&i_minus_kh, &self.p);
        Some(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Model Predictive Control
// ─────────────────────────────────────────────────────────────────────────────

/// Simple linear MPC with prediction horizon N, tracking and input penalties.
#[derive(Debug, Clone)]
pub struct ModelPredictiveControl {
    /// System A matrix (n×n).
    pub a: Vec<Vec<f64>>,
    /// System B matrix (n×m).
    pub b: Vec<Vec<f64>>,
    /// Prediction horizon.
    pub horizon: usize,
    /// State tracking weight (diagonal Q).
    pub q_weight: f64,
    /// Input weight (diagonal R).
    pub r_weight: f64,
    /// Input constraint: max magnitude.
    pub u_max: f64,
}

impl ModelPredictiveControl {
    /// Create a new MPC controller.
    pub fn new(
        a: Vec<Vec<f64>>,
        b: Vec<Vec<f64>>,
        horizon: usize,
        q_weight: f64,
        r_weight: f64,
    ) -> Self {
        Self {
            a,
            b,
            horizon,
            q_weight,
            r_weight,
            u_max: f64::MAX,
        }
    }

    /// Set input constraint.
    pub fn set_input_constraint(&mut self, u_max: f64) {
        self.u_max = u_max;
    }

    /// Compute first optimal control action via unconstrained MPC (gradient step).
    /// `x` is current state, `x_ref` is reference trajectory (constant horizon×n).
    pub fn compute(&self, x: &[f64], x_ref: &[f64]) -> Vec<f64> {
        let n = x.len();
        let m = if self.b.is_empty() {
            0
        } else {
            self.b[0].len()
        };
        if m == 0 {
            return vec![];
        }

        // Receding horizon: predict states and compute gradient numerically
        // Simple steepest-descent for unconstrained linear MPC
        let mut u_seq: Vec<Vec<f64>> = vec![vec![0.0; m]; self.horizon];
        let alpha = 0.01; // step size

        for _iter in 0..50 {
            // Evaluate cost gradient
            let mut x_pred = x.to_vec();
            let mut grads: Vec<Vec<f64>> = vec![vec![0.0; m]; self.horizon];

            for k in 0..self.horizon {
                let ax = mat_vec_mul(&self.a, &x_pred);
                let bu = mat_vec_mul(&self.b, &u_seq[k]);
                x_pred = vec_add(&ax, &bu);
                // State error
                let err: Vec<f64> = x_pred
                    .iter()
                    .zip(x_ref.iter().take(n))
                    .map(|(xi, ri)| xi - ri)
                    .collect();
                let cost_grad_x = vec_scale(&err, 2.0 * self.q_weight);
                // Propagate gradient back to u: dJ/du_k ≈ Bᵀ * dJ/dx_{k+1}
                let bt = mat_transpose(&self.b);
                let gu = mat_vec_mul(&bt, &cost_grad_x);
                // Add input penalty gradient
                let gu_r: Vec<f64> = u_seq[k]
                    .iter()
                    .zip(gu.iter())
                    .map(|(ui, gi)| gi + 2.0 * self.r_weight * ui)
                    .collect();
                grads[k] = gu_r;
            }
            // Update u_seq
            for k in 0..self.horizon {
                for j in 0..m {
                    u_seq[k][j] -= alpha * grads[k][j];
                    u_seq[k][j] = u_seq[k][j].clamp(-self.u_max, self.u_max);
                }
            }
        }

        u_seq[0].clone()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Adaptive Controller (MRAS / MIT rule)
// ─────────────────────────────────────────────────────────────────────────────

/// Model Reference Adaptive System (MRAS) controller using the MIT rule for
/// parameter adaptation.
#[derive(Debug, Clone)]
pub struct AdaptiveController {
    /// Reference model numerator (first order: \[b_m\]).
    pub bm: f64,
    /// Reference model denominator pole (am).
    pub am: f64,
    /// Adaptive parameter θ (scalar for SISO first order).
    pub theta: f64,
    /// Adaptation gain γ.
    pub gamma: f64,
    /// Reference model state.
    pub ym: f64,
    /// Plant output.
    pub y: f64,
    /// Integral of sensitivity signal.
    pub sens_int: f64,
}

impl AdaptiveController {
    /// Create with reference model parameters and adaptation gain.
    pub fn new(bm: f64, am: f64, gamma: f64) -> Self {
        Self {
            bm,
            am,
            theta: 1.0,
            gamma,
            ym: 0.0,
            y: 0.0,
            sens_int: 0.0,
        }
    }

    /// Update the adaptive controller.
    /// `r` is reference input, `y_plant` is measured plant output, `dt` is time step.
    /// Returns the control signal u.
    pub fn update(&mut self, r: f64, y_plant: f64, dt: f64) -> f64 {
        // Reference model: ẏ_m = -a_m * y_m + b_m * r
        let dym = -self.am * self.ym + self.bm * r;
        self.ym += dym * dt;

        // Error
        let e = y_plant - self.ym;

        // MIT rule: dθ/dt = -γ * e * ∂e/∂θ ≈ -γ * e * y_m
        let dtheta = -self.gamma * e * self.ym;
        self.theta += dtheta * dt;

        self.y = y_plant;
        // Control signal
        self.theta * r
    }

    /// Reset state.
    pub fn reset(&mut self) {
        self.ym = 0.0;
        self.y = 0.0;
        self.theta = 1.0;
        self.sens_int = 0.0;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fuzzy Controller (Mamdani)
// ─────────────────────────────────────────────────────────────────────────────

/// Triangular membership function: zero outside \[a, c\], peak at b.
#[derive(Debug, Clone)]
pub struct TriangleMf {
    /// Left foot.
    pub a: f64,
    /// Peak.
    pub b: f64,
    /// Right foot.
    pub c: f64,
}

impl TriangleMf {
    /// Create a triangular membership function.
    pub fn new(a: f64, b: f64, c: f64) -> Self {
        Self { a, b, c }
    }

    /// Evaluate μ(x).
    pub fn eval(&self, x: f64) -> f64 {
        if x <= self.a || x >= self.c {
            0.0
        } else if x <= self.b {
            (x - self.a) / (self.b - self.a)
        } else {
            (self.c - x) / (self.c - self.b)
        }
    }
}

/// Fuzzy rule: (input_mf_index) → (output_mf_index).
#[derive(Debug, Clone)]
pub struct FuzzyRule {
    /// Index into input membership functions.
    pub input_idx: usize,
    /// Index into output membership functions.
    pub output_idx: usize,
}

/// Mamdani fuzzy controller with triangular MFs and centroid defuzzification.
#[derive(Debug, Clone)]
pub struct FuzzyController {
    /// Input membership functions.
    pub input_mfs: Vec<TriangleMf>,
    /// Output membership functions.
    pub output_mfs: Vec<TriangleMf>,
    /// Fuzzy rules.
    pub rules: Vec<FuzzyRule>,
    /// Output universe lower bound.
    pub out_min: f64,
    /// Output universe upper bound.
    pub out_max: f64,
    /// Number of discretization points for defuzzification.
    pub n_points: usize,
}

impl FuzzyController {
    /// Create a new fuzzy controller.
    pub fn new(
        input_mfs: Vec<TriangleMf>,
        output_mfs: Vec<TriangleMf>,
        rules: Vec<FuzzyRule>,
        out_min: f64,
        out_max: f64,
    ) -> Self {
        Self {
            input_mfs,
            output_mfs,
            rules,
            out_min,
            out_max,
            n_points: 100,
        }
    }

    /// Compute output via Mamdani inference and centroid defuzzification.
    pub fn compute(&self, input: f64) -> f64 {
        // Fuzzify input
        let memberships: Vec<f64> = self.input_mfs.iter().map(|mf| mf.eval(input)).collect();

        // Activate output MFs (min operation for each rule)
        let mut activated = vec![0.0_f64; self.output_mfs.len()];
        for rule in &self.rules {
            if rule.input_idx < memberships.len() && rule.output_idx < activated.len() {
                activated[rule.output_idx] =
                    activated[rule.output_idx].max(memberships[rule.input_idx]);
            }
        }

        // Centroid defuzzification
        let step = (self.out_max - self.out_min) / self.n_points as f64;
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for i in 0..self.n_points {
            let y = self.out_min + (i as f64 + 0.5) * step;
            // Aggregate: max of clipped output MFs
            let mu = self
                .output_mfs
                .iter()
                .enumerate()
                .fold(0.0_f64, |acc, (j, mf)| {
                    let clipped = mf.eval(y).min(activated[j]);
                    acc.max(clipped)
                });
            num += y * mu;
            den += mu;
        }
        if den < 1e-12 {
            (self.out_min + self.out_max) / 2.0
        } else {
            num / den
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Neural Controller (MLP)
// ─────────────────────────────────────────────────────────────────────────────

/// Simple experience replay buffer for neural controller training.
#[derive(Debug, Clone)]
pub struct ReplayBuffer {
    /// Stored (state, action, reward) tuples.
    pub buffer: Vec<(Vec<f64>, Vec<f64>, f64)>,
    /// Maximum capacity.
    pub capacity: usize,
    /// Write index.
    pub idx: usize,
}

impl ReplayBuffer {
    /// Create a new replay buffer with given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            capacity,
            idx: 0,
        }
    }

    /// Push a new experience.
    pub fn push(&mut self, state: Vec<f64>, action: Vec<f64>, reward: f64) {
        let entry = (state, action, reward);
        if self.buffer.len() < self.capacity {
            self.buffer.push(entry);
        } else {
            self.buffer[self.idx % self.capacity] = entry;
        }
        self.idx += 1;
    }

    /// Number of stored experiences.
    pub fn len(&self) -> usize {
        self.buffer.len()
    }

    /// Check if buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }
}

/// Simple multi-layer perceptron (MLP) neural network controller.
#[derive(Debug, Clone)]
pub struct NeuralController {
    /// Layer weights: each element is a weight matrix (out×in).
    pub weights: Vec<Vec<Vec<f64>>>,
    /// Layer biases: each element is a bias vector.
    pub biases: Vec<Vec<f64>>,
    /// Learning rate.
    pub lr: f64,
    /// Replay buffer.
    pub replay: ReplayBuffer,
}

impl NeuralController {
    /// Create an MLP controller with given layer sizes.
    /// `layer_sizes` includes input and output dimensions.
    pub fn new(layer_sizes: &[usize], lr: f64) -> Self {
        let n_layers = layer_sizes.len() - 1;
        let mut weights = Vec::with_capacity(n_layers);
        let mut biases = Vec::with_capacity(n_layers);
        for i in 0..n_layers {
            let rows = layer_sizes[i + 1];
            let cols = layer_sizes[i];
            // Xavier initialization
            let scale = (2.0 / (rows + cols) as f64).sqrt();
            let w: Vec<Vec<f64>> = (0..rows)
                .map(|r| {
                    (0..cols)
                        .map(|c| ((r * cols + c) as f64 * 0.1).sin() * scale)
                        .collect()
                })
                .collect();
            weights.push(w);
            biases.push(vec![0.0; rows]);
        }
        Self {
            weights,
            biases,
            lr,
            replay: ReplayBuffer::new(10000),
        }
    }

    /// ReLU activation.
    fn relu(x: f64) -> f64 {
        x.max(0.0)
    }

    /// Tanh activation.
    fn tanh(x: f64) -> f64 {
        x.tanh()
    }

    /// Forward pass: returns output and all layer activations (for backprop).
    pub fn forward(&self, input: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>) {
        let mut activations = vec![input.to_vec()];
        let n_layers = self.weights.len();
        for layer in 0..n_layers {
            let w = &self.weights[layer];
            let b = &self.biases[layer];
            let rows = w.len();
            let z: Vec<f64> = (0..rows)
                .map(|i| {
                    let s: f64 = w[i]
                        .iter()
                        .zip(activations.last().expect("activations is non-empty").iter())
                        .map(|(wi, ai)| wi * ai)
                        .sum();
                    s + b[i]
                })
                .collect();
            // Use tanh on hidden, linear on output
            let a = if layer < n_layers - 1 {
                z.iter().map(|&zi| Self::tanh(zi)).collect()
            } else {
                z
            };
            activations.push(a);
        }
        let out = activations
            .last()
            .expect("activations is non-empty")
            .clone();
        (out, activations)
    }

    /// Gradient update step (simple policy gradient / supervised).
    /// `target` is the desired output, `input` is the current state.
    pub fn gradient_update(&mut self, input: &[f64], target: &[f64]) {
        let n_layers = self.weights.len();
        let (_, activations) = self.forward(input);

        // MSE loss gradient at output layer
        let out = &activations[n_layers];
        let mut delta: Vec<f64> = out.iter().zip(target.iter()).map(|(o, t)| o - t).collect();

        for layer in (0..n_layers).rev() {
            let a_prev = &activations[layer];
            let rows = self.weights[layer].len();
            // Weight gradient
            let new_delta: Vec<f64> = if layer > 0 {
                let cols = self.weights[layer][0].len();
                let mut nd = vec![0.0; cols];
                for i in 0..rows {
                    for j in 0..cols {
                        nd[j] += self.weights[layer][i][j] * delta[i];
                    }
                }
                // Backprop through tanh: (1 - a²)
                nd.iter()
                    .zip(activations[layer].iter())
                    .map(|(d, a)| d * (1.0 - a * a))
                    .collect()
            } else {
                vec![0.0; a_prev.len()]
            };
            // Update weights and biases
            for i in 0..rows {
                for j in 0..a_prev.len() {
                    self.weights[layer][i][j] -= self.lr * delta[i] * a_prev[j];
                }
                self.biases[layer][i] -= self.lr * delta[i];
            }
            delta = new_delta;
        }
    }

    /// Store experience in replay buffer.
    pub fn remember(&mut self, state: Vec<f64>, action: Vec<f64>, reward: f64) {
        self.replay.push(state, action, reward);
    }

    /// Compute control output for given state.
    pub fn compute(&self, state: &[f64]) -> Vec<f64> {
        let (out, _) = self.forward(state);
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Zero-order hold (ZOH) discretization of continuous-time A matrix.
/// Uses Padé approximation: e^(A*dt) ≈ I + A*dt + (A*dt)²/2 + ...
pub fn zoh_discretize(a: &[Vec<f64>], dt: f64) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut result = eye(n);
    let mut term = eye(n);
    let a_dt = mat_scale(a, dt);
    for k in 1..=10 {
        term = mat_mul(&term, &a_dt);
        let scaled = mat_scale(&term, 1.0 / factorial(k) as f64);
        result = mat_add(&result, &scaled);
    }
    result
}

fn factorial(n: u64) -> u64 {
    (1..=n).product()
}

/// Solve the Discrete Algebraic Riccati Equation (DARE) iteratively.
/// Returns gain matrix K = (R + BᵀPB)⁻¹ * BᵀPA.
pub fn riccati_solve_dare(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    q: &[Vec<f64>],
    r: &[Vec<f64>],
    max_iter: usize,
    tol: f64,
) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let m = if b.is_empty() {
        return None;
    } else {
        b[0].len()
    };
    let _ = m;

    let at = mat_transpose(a);
    let bt = mat_transpose(b);

    let mut p = q.to_vec();
    for _ in 0..max_iter {
        // K = (R + Bᵀ P B)⁻¹ Bᵀ P A
        let pb = mat_mul(&p, b);
        let btpb = mat_mul(&bt, &pb);
        let r_btpb = mat_add(r, &btpb);
        let r_btpb_inv = mat_inv(&r_btpb)?;
        let btp = mat_mul(&bt, &p);
        let btpa = mat_mul(&btp, a);
        let k = mat_mul(&r_btpb_inv, &btpa);

        // P_new = Q + Aᵀ P A - Aᵀ P B K
        let atpa = mat_mul(&at, &mat_mul(&p, a));
        let atpbk = mat_mul(&at, &mat_mul(&pb, &k));
        let p_new = mat_add(q, &mat_add(&atpa, &mat_scale(&atpbk, -1.0)));

        // Check convergence
        let mut max_diff = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                max_diff = max_diff.max((p_new[i][j] - p[i][j]).abs());
            }
        }
        p = p_new;
        if max_diff < tol {
            // Compute final K
            let pb2 = mat_mul(&p, b);
            let btpb2 = mat_mul(&bt, &pb2);
            let r_btpb2 = mat_add(r, &btpb2);
            let r_btpb2_inv = mat_inv(&r_btpb2)?;
            let btp2 = mat_mul(&bt, &p);
            let btpa2 = mat_mul(&btp2, a);
            let k_final = mat_mul(&r_btpb2_inv, &btpa2);
            return Some(k_final);
        }
    }
    None
}

/// Compute controllability matrix: \[B, AB, A²B, ..., A^(n-1)B\].
pub fn controllability_matrix(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = if b.is_empty() {
        return vec![];
    } else {
        b[0].len()
    };
    let total_cols = n * m;
    let mut c_mat = vec![vec![0.0; total_cols]; n];

    let mut ak_b = b.to_vec();
    for k in 0..n {
        for i in 0..n {
            for j in 0..m {
                c_mat[i][k * m + j] = ak_b[i][j];
            }
        }
        if k < n - 1 {
            ak_b = mat_mul(a, &ak_b);
        }
    }
    c_mat
}

/// Pole placement: compute state feedback gain K such that eigenvalues of
/// (A - B*K) match desired poles (Ackermann's formula, SISO only).
pub fn pole_placement(a: &[Vec<f64>], b: &[Vec<f64>], poles: &[f64]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if b.is_empty() || b[0].len() != 1 || poles.len() != n {
        return None; // Only SISO supported
    }

    // Desired characteristic polynomial
    let mut p_coeffs = vec![1.0_f64];
    for &pole in poles {
        // Multiply by (s - pole): shift and subtract
        let mut new_coeffs = vec![0.0; p_coeffs.len() + 1];
        for (i, &c) in p_coeffs.iter().enumerate() {
            new_coeffs[i] += c;
            new_coeffs[i + 1] -= c * pole;
        }
        p_coeffs = new_coeffs;
    }

    // Controllability matrix
    let c_mat = controllability_matrix(a, b);

    // Characteristic polynomial of A (Cayley-Hamilton)
    let a_poly = char_poly(a);

    // Ackermann's formula: K = e_n * C^-1 * phi(A)
    // phi(A) = A^n + alpha_1 * A^(n-1) + ... + alpha_n * I
    let c_inv = mat_inv(&c_mat)?;

    // Evaluate phi(A)
    let mut phi_a = mat_scale(&eye(n), p_coeffs[n]);
    let mut ak = eye(n);
    for k in (0..n).rev() {
        ak = mat_mul(&ak, a);
        phi_a = mat_add(&phi_a, &mat_scale(&ak, p_coeffs[k]));
    }
    let _ = a_poly;

    // e_n = last row of identity
    let en: Vec<Vec<f64>> = vec![(0..n).map(|j| if j == n - 1 { 1.0 } else { 0.0 }).collect()];
    let en_cinv = mat_mul(&en, &c_inv);
    let k = mat_mul(&en_cinv, &phi_a);
    Some(k)
}

/// Compute characteristic polynomial coefficients of matrix A (length n+1, leading 1).
fn char_poly(a: &[Vec<f64>]) -> Vec<f64> {
    let n = a.len();
    // Use Faddeev-LeVerrier algorithm
    let mut m = eye(n);
    let mut coeffs = vec![1.0_f64];
    for k in 1..=n {
        let am = mat_mul(a, &m);
        let tr: f64 = (0..n).map(|i| am[i][i]).sum();
        let c = -tr / k as f64;
        coeffs.push(c);
        m = mat_add(&am, &mat_scale(&eye(n), c));
    }
    coeffs
}

/// Constant π exported for convenience.
pub const TWO_PI: f64 = 2.0 * PI;

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── PID tests ────────────────────────────────────────────────────────────

    #[test]
    fn pid_proportional_only() {
        let mut pid = Pid::new(2.0, 0.0, 0.0);
        let out = pid.compute(1.0, 0.01);
        assert!((out - 2.0).abs() < 1e-10, "P-only: out={out}");
    }

    #[test]
    fn pid_integral_accumulates() {
        let mut pid = Pid::new(0.0, 1.0, 0.0);
        pid.compute(1.0, 0.1);
        pid.compute(1.0, 0.1);
        let out = pid.compute(1.0, 0.1);
        // integral ≈ 0.3, output = ki * integral = 0.3
        assert!((out - 0.3).abs() < 1e-6, "I-only: out={out}");
    }

    #[test]
    fn pid_anti_windup_clamps_integral() {
        let mut pid = Pid::new(0.0, 1.0, 0.0);
        pid.set_integral_limit(1.0);
        for _ in 0..100 {
            pid.compute(1.0, 0.1);
        }
        assert!(pid.integral <= 1.0 + 1e-10);
    }

    #[test]
    fn pid_output_saturation() {
        let mut pid = Pid::new(100.0, 0.0, 0.0);
        pid.set_output_limits(-1.0, 1.0);
        let out = pid.compute(10.0, 0.01);
        assert!((out - 1.0).abs() < 1e-10, "saturated output: {out}");
    }

    #[test]
    fn pid_reset_clears_state() {
        let mut pid = Pid::new(1.0, 1.0, 1.0);
        pid.compute(5.0, 0.1);
        pid.reset();
        assert!(pid.integral.abs() < 1e-12);
        assert!(pid.prev_error.abs() < 1e-12);
    }

    #[test]
    fn pid_zero_dt_returns_zero() {
        let mut pid = Pid::new(1.0, 1.0, 1.0);
        let out = pid.compute(1.0, 0.0);
        assert!((out).abs() < 1e-12);
    }

    #[test]
    fn pid_set_gains() {
        let mut pid = Pid::new(1.0, 0.0, 0.0);
        pid.set_gains(2.0, 3.0, 4.0);
        assert!((pid.kp - 2.0).abs() < 1e-12);
        assert!((pid.ki - 3.0).abs() < 1e-12);
        assert!((pid.kd - 4.0).abs() < 1e-12);
    }

    // ── StateSpaceModel tests ────────────────────────────────────────────────

    #[test]
    fn state_space_step_first_order() {
        // ẋ = -x + u, y = x (integrator approx)
        let a = vec![vec![-1.0]];
        let b = vec![vec![1.0]];
        let c = vec![vec![1.0]];
        let d = vec![vec![0.0]];
        let mut ss = StateSpaceModel::new(a, b, c, d);
        let y = ss.step(&[1.0], 0.01);
        // x += dt*(-x+u) = 0 + 0.01*(0+1) = 0.01
        assert!((y[0] - 0.01).abs() < 1e-10, "y={}", y[0]);
    }

    #[test]
    fn state_space_n_states() {
        let a = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        let b = vec![vec![0.0], vec![1.0]];
        let c = vec![vec![1.0, 0.0]];
        let d = vec![vec![0.0]];
        let ss = StateSpaceModel::new(a, b, c, d);
        assert_eq!(ss.n_states(), 2);
        assert_eq!(ss.n_inputs(), 1);
        assert_eq!(ss.n_outputs(), 1);
    }

    #[test]
    fn state_space_discrete_step() {
        let a = vec![vec![0.5]];
        let b = vec![vec![1.0]];
        let c = vec![vec![1.0]];
        let d = vec![vec![0.0]];
        let mut ss = StateSpaceModel::new(a, b, c, d);
        ss.is_discrete = true;
        ss.set_state(vec![1.0]);
        let y = ss.step(&[0.0], 0.01);
        // x_new = 0.5 * 1 + 1.0 * 0 = 0.5
        assert!((y[0] - 0.5).abs() < 1e-10, "y={}", y[0]);
    }

    // ── TransferFunction tests ───────────────────────────────────────────────

    #[test]
    fn tf_dc_gain_unity() {
        let tf = TransferFunction::new(vec![1.0], vec![1.0]);
        assert!((tf.dc_gain() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn tf_dc_gain_ratio() {
        let tf = TransferFunction::new(vec![2.0], vec![4.0]);
        assert!((tf.dc_gain() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn tf_from_pid_structure() {
        let tf = TransferFunction::from_pid(1.0, 2.0, 0.5);
        // num = [0.5, 1, 2], den = [1, 0]
        assert_eq!(tf.num.len(), 3);
        assert_eq!(tf.den.len(), 2);
    }

    #[test]
    fn tf_frequency_response_at_zero() {
        let tf = TransferFunction::new(vec![3.0], vec![1.0]);
        let (mag, phase) = tf.frequency_response(0.001);
        assert!((mag - 3.0).abs() < 0.01, "mag={mag}");
        assert!(phase.abs() < 1.0, "phase={phase}");
    }

    // ── BodeAnalysis tests ───────────────────────────────────────────────────

    #[test]
    fn bode_magnitude_db_at_low_freq() {
        let tf = TransferFunction::new(vec![10.0], vec![1.0]);
        let bode = BodeAnalysis::new(tf);
        let db = bode.magnitude_db(0.001);
        // Should be ~20 dB
        assert!((db - 20.0).abs() < 1.0, "db={db}");
    }

    #[test]
    fn bode_sweep_returns_correct_length() {
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]);
        let bode = BodeAnalysis::new(tf);
        let (omegas, mags, phases) = bode.sweep(0.1, 100.0, 50);
        assert_eq!(omegas.len(), 50);
        assert_eq!(mags.len(), 50);
        assert_eq!(phases.len(), 50);
    }

    // ── KalmanFilter tests ───────────────────────────────────────────────────

    #[test]
    fn kalman_predict_zero_input() {
        let mut kf = KalmanFilter::new(vec![1.0], vec![vec![1.0]]);
        let a = vec![vec![1.0]];
        let b = vec![vec![0.0]];
        let q = vec![vec![0.01]];
        kf.predict(&a, &b, &[0.0], &q);
        // x_hat should still be 1.0
        assert!((kf.x_hat[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn kalman_update_corrects_state() {
        let mut kf = KalmanFilter::new(vec![0.0], vec![vec![1.0]]);
        let h = vec![vec![1.0]];
        let r = vec![vec![1.0]];
        kf.update(&h, &r, &[10.0]).unwrap();
        // Should move toward measurement
        assert!(
            kf.x_hat[0] > 0.0,
            "state should increase toward measurement"
        );
    }

    #[test]
    fn kalman_update_reduces_covariance() {
        let mut kf = KalmanFilter::new(vec![0.0], vec![vec![10.0]]);
        let h = vec![vec![1.0]];
        let r = vec![vec![1.0]];
        let p_before = kf.p[0][0];
        kf.update(&h, &r, &[0.0]).unwrap();
        assert!(
            kf.p[0][0] < p_before,
            "covariance should decrease after update"
        );
    }

    // ── LQR tests ────────────────────────────────────────────────────────────

    #[test]
    fn lqr_scalar_system() {
        // Scalar: a=1.1, b=1, q=1, r=0.1 → should converge
        let a = vec![vec![0.9]];
        let b = vec![vec![1.0]];
        let q = vec![vec![1.0]];
        let r = vec![vec![0.1]];
        let lqr = LqrController::solve(&a, &b, q, r);
        assert!(
            lqr.is_some(),
            "LQR should converge for stable scalar system"
        );
    }

    #[test]
    fn lqr_control_output_sign() {
        let a = vec![vec![0.9]];
        let b = vec![vec![1.0]];
        let q = vec![vec![1.0]];
        let r = vec![vec![0.1]];
        if let Some(lqr) = LqrController::solve(&a, &b, q, r) {
            let u = lqr.compute(&[1.0]);
            // Stabilizing control should push state negative
            assert!(u[0] < 0.0, "LQR should apply negative feedback: u={}", u[0]);
        }
    }

    // ── AdaptiveController tests ─────────────────────────────────────────────

    #[test]
    fn adaptive_controller_produces_output() {
        let mut ctrl = AdaptiveController::new(1.0, 1.0, 0.1);
        let u = ctrl.update(1.0, 0.5, 0.01);
        assert!(u.is_finite(), "u={u}");
    }

    #[test]
    fn adaptive_controller_reset() {
        let mut ctrl = AdaptiveController::new(1.0, 1.0, 0.1);
        ctrl.update(5.0, 3.0, 0.1);
        ctrl.reset();
        assert!(ctrl.ym.abs() < 1e-12);
        assert!((ctrl.theta - 1.0).abs() < 1e-12);
    }

    // ── FuzzyController tests ────────────────────────────────────────────────

    #[test]
    fn triangle_mf_peak_value() {
        let mf = TriangleMf::new(0.0, 1.0, 2.0);
        assert!((mf.eval(1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn triangle_mf_zero_outside() {
        let mf = TriangleMf::new(0.0, 1.0, 2.0);
        assert!(mf.eval(-0.1).abs() < 1e-10);
        assert!(mf.eval(2.1).abs() < 1e-10);
    }

    #[test]
    fn triangle_mf_half_left() {
        let mf = TriangleMf::new(0.0, 1.0, 2.0);
        assert!((mf.eval(0.5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn fuzzy_controller_center_input() {
        let imfs = vec![
            TriangleMf::new(-2.0, -1.0, 0.0),
            TriangleMf::new(-1.0, 0.0, 1.0),
            TriangleMf::new(0.0, 1.0, 2.0),
        ];
        let omfs = vec![
            TriangleMf::new(-2.0, -1.0, 0.0),
            TriangleMf::new(-1.0, 0.0, 1.0),
            TriangleMf::new(0.0, 1.0, 2.0),
        ];
        let rules = vec![
            FuzzyRule {
                input_idx: 0,
                output_idx: 0,
            },
            FuzzyRule {
                input_idx: 1,
                output_idx: 1,
            },
            FuzzyRule {
                input_idx: 2,
                output_idx: 2,
            },
        ];
        let ctrl = FuzzyController::new(imfs, omfs, rules, -2.0, 2.0);
        let out = ctrl.compute(0.0);
        // With input = 0, rule 2 fires with membership 1, output should be ~0
        assert!(out.abs() < 0.5, "out={out}");
    }

    // ── NeuralController tests ───────────────────────────────────────────────

    #[test]
    fn neural_controller_forward_pass() {
        let nn = NeuralController::new(&[2, 4, 1], 0.01);
        let (out, _acts) = nn.forward(&[1.0, 0.5]);
        assert_eq!(out.len(), 1);
        assert!(out[0].is_finite());
    }

    #[test]
    fn neural_controller_gradient_update() {
        let mut nn = NeuralController::new(&[2, 4, 1], 0.01);
        let (out_before, _) = nn.forward(&[1.0, 0.5]);
        nn.gradient_update(&[1.0, 0.5], &[1.0]);
        let (out_after, _) = nn.forward(&[1.0, 0.5]);
        // Loss should decrease after one step toward target=1.0
        let loss_before = (out_before[0] - 1.0).powi(2);
        let loss_after = (out_after[0] - 1.0).powi(2);
        assert!(loss_after <= loss_before + 1e-6, "loss should not increase");
    }

    #[test]
    fn replay_buffer_capacity() {
        let mut buf = ReplayBuffer::new(3);
        for i in 0..5 {
            buf.push(vec![i as f64], vec![0.0], 0.0);
        }
        assert_eq!(buf.len(), 3, "buffer should not exceed capacity");
    }

    #[test]
    fn replay_buffer_is_empty() {
        let buf = ReplayBuffer::new(10);
        assert!(buf.is_empty());
    }

    // ── Helper function tests ────────────────────────────────────────────────

    #[test]
    fn zoh_discretize_identity_zero_dt() {
        let a = vec![vec![1.0, 0.0], vec![0.0, -1.0]];
        let ad = zoh_discretize(&a, 0.0);
        // e^0 ≈ I
        assert!((ad[0][0] - 1.0).abs() < 1e-6);
        assert!((ad[1][1] - 1.0).abs() < 1e-6);
        assert!(ad[0][1].abs() < 1e-6);
    }

    #[test]
    fn controllability_matrix_size() {
        let a = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        let b = vec![vec![0.0], vec![1.0]];
        let cm = controllability_matrix(&a, &b);
        assert_eq!(cm.len(), 2);
        assert_eq!(cm[0].len(), 2); // 2×(n*m)
    }

    #[test]
    fn mat_inv_identity() {
        let i = eye(3);
        let inv = mat_inv(&i).unwrap();
        for r in 0..3 {
            for c in 0..3 {
                let expected = if r == c { 1.0 } else { 0.0 };
                assert!((inv[r][c] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn mat_mul_identity() {
        let a = vec![vec![2.0, 3.0], vec![1.0, 4.0]];
        let i = eye(2);
        let ai = mat_mul(&a, &i);
        assert!((ai[0][0] - 2.0).abs() < 1e-12);
        assert!((ai[1][1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn poly_roots_linear() {
        // x + 2 = 0 → root = -2
        let roots = poly_roots_real(&[1.0, 2.0]);
        assert_eq!(roots.len(), 1);
        assert!((roots[0] + 2.0).abs() < 1e-10);
    }

    #[test]
    fn poly_roots_quadratic() {
        // x^2 - 5x + 6 = (x-2)(x-3) → roots 2, 3
        let roots = poly_roots_real(&[1.0, -5.0, 6.0]);
        assert_eq!(roots.len(), 2);
        let mut sorted = roots.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 2.0).abs() < 1e-6, "root0={}", sorted[0]);
        assert!((sorted[1] - 3.0).abs() < 1e-6, "root1={}", sorted[1]);
    }

    #[test]
    fn mpc_compute_output_length() {
        let a = vec![vec![1.0, 0.1], vec![0.0, 1.0]];
        let b = vec![vec![0.0], vec![0.1]];
        let mpc = ModelPredictiveControl::new(a, b, 5, 1.0, 0.1);
        let u = mpc.compute(&[1.0, 0.0], &[0.0, 0.0]);
        assert_eq!(u.len(), 1);
    }

    #[test]
    fn tf_poles_first_order() {
        // H(s) = 1 / (s + 3) → pole at -3
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 3.0]);
        let poles = tf.poles();
        assert_eq!(poles.len(), 1);
        assert!((poles[0] + 3.0).abs() < 1e-6, "pole={}", poles[0]);
    }

    #[test]
    fn bode_phase_first_order_system() {
        // First-order: H(s) = 1/(s+1), at ω=1: phase = -45°
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]);
        let bode = BodeAnalysis::new(tf);
        let phase = bode.phase_deg(1.0);
        assert!((phase + 45.0).abs() < 1.0, "phase={phase}");
    }

    #[test]
    fn nn_compute_shape() {
        let nn = NeuralController::new(&[3, 8, 2], 0.001);
        let out = nn.compute(&[1.0, 2.0, 3.0]);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn pi_only_settles() {
        let mut pid = Pid::new(1.0, 0.5, 0.0);
        // Simulate simple integrating plant: y += u * dt
        let mut y = 0.0_f64;
        let setpoint = 1.0;
        let dt = 0.01;
        for _ in 0..2000 {
            let err = setpoint - y;
            let output = pid.compute(err, dt);
            y += output * dt;
        }
        assert!((y - setpoint).abs() < 0.05, "PI settled y={y}");
    }

    #[test]
    fn state_space_zoh_discrete() {
        let a = vec![vec![-1.0]];
        let b = vec![vec![1.0]];
        let c = vec![vec![1.0]];
        let d = vec![vec![0.0]];
        let ss = StateSpaceModel::new(a, b, c, d);
        let disc = ss.discretize_zoh(0.01);
        assert!(disc.is_discrete);
        assert!((disc.a[0][0] - 1.0).abs() < 0.1); // e^(-0.01) ≈ 0.99
    }
}
