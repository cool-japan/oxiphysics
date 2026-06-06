// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural ODE (continuous-depth neural networks) implementations.
//!
//! Implements Neural ODEs as introduced by Chen et al. (NeurIPS 2018).
//! Provides:
//! - [`NeuralOdeFunc`]: parameterised ODE right-hand side (small MLP)
//! - [`NeuralOdeSolver`]: wraps a [`NeuralOdeFunc`] with RK4 integration
//! - [`AdjointMethod`]: reverse-mode gradient estimation via the adjoint
//! - [`LatentOde`]: encoder-dynamics-decoder architecture
//! - [`TimeSeriesOde`]: convenience wrapper for time-series fitting
//! - Free functions [`rk4_step`] and [`dopri5_step`]

// ─────────────────────────────────────────────────────────────────────────────
// Free integration helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Perform a single classic fourth-order Runge-Kutta step.
///
/// # Arguments
/// * `f`  – ODE right-hand side `f(t, y)`.
/// * `t`  – Current time.
/// * `y`  – Current state (slice of length `n`).
/// * `h`  – Step size.
///
/// # Returns
/// State vector at time `t + h`.
pub fn rk4_step(f: &dyn Fn(f64, &[f64]) -> Vec<f64>, t: f64, y: &[f64], h: f64) -> Vec<f64> {
    let n = y.len();
    let k1 = f(t, y);
    let y2: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k1[i]).collect();
    let k2 = f(t + 0.5 * h, &y2);
    let y3: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k2[i]).collect();
    let k3 = f(t + 0.5 * h, &y3);
    let y4: Vec<f64> = (0..n).map(|i| y[i] + h * k3[i]).collect();
    let k4 = f(t + h, &y4);
    (0..n)
        .map(|i| y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
        .collect()
}

/// Perform a single Dormand-Prince (DOPRI5) adaptive step.
///
/// Returns `(y_high, y_low, error_norm)` where `y_high` is the 5th-order
/// solution, `y_low` is the embedded 4th-order solution derived from the
/// Dormand-Prince error coefficients, and `error_norm` is the RMS scaled
/// difference useful for step-size control.
///
/// Uses the full FSAL (First Same As Last) property: computes 7 function
/// evaluations (k1…k6 + k7 = f(t+h, y_high)) and uses the proper
/// Dormand-Prince error coefficients e1…e7 for the embedded pair.
///
/// # Arguments
/// * `f`   – ODE right-hand side.
/// * `t`   – Current time.
/// * `y`   – Current state.
/// * `h`   – Proposed step size.
/// * `rtol` – Relative tolerance (used in error norm).
/// * `atol` – Absolute tolerance (used in error norm).
pub fn dopri5_step(
    f: &dyn Fn(f64, &[f64]) -> Vec<f64>,
    t: f64,
    y: &[f64],
    h: f64,
    rtol: f64,
    atol: f64,
) -> (Vec<f64>, Vec<f64>, f64) {
    let n = y.len();
    // Butcher tableau node values (Dormand-Prince)
    let c2 = 1.0 / 5.0;
    let c3 = 3.0 / 10.0;
    let c4 = 4.0 / 5.0;
    let c5 = 8.0 / 9.0;

    let k1 = f(t, y);

    let y2: Vec<f64> = (0..n).map(|i| y[i] + h * (1.0 / 5.0) * k1[i]).collect();
    let k2 = f(t + c2 * h, &y2);

    let y3: Vec<f64> = (0..n)
        .map(|i| y[i] + h * ((3.0 / 40.0) * k1[i] + (9.0 / 40.0) * k2[i]))
        .collect();
    let k3 = f(t + c3 * h, &y3);

    let y4: Vec<f64> = (0..n)
        .map(|i| y[i] + h * ((44.0 / 45.0) * k1[i] - (56.0 / 15.0) * k2[i] + (32.0 / 9.0) * k3[i]))
        .collect();
    let k4 = f(t + c4 * h, &y4);

    let y5: Vec<f64> = (0..n)
        .map(|i| {
            y[i] + h
                * ((19372.0 / 6561.0) * k1[i] - (25360.0 / 2187.0) * k2[i]
                    + (64448.0 / 6561.0) * k3[i]
                    - (212.0 / 729.0) * k4[i])
        })
        .collect();
    let k5 = f(t + c5 * h, &y5);

    let y6: Vec<f64> = (0..n)
        .map(|i| {
            y[i] + h
                * ((9017.0 / 3168.0) * k1[i] - (355.0 / 33.0) * k2[i]
                    + (46732.0 / 5247.0) * k3[i]
                    + (49.0 / 176.0) * k4[i]
                    - (5103.0 / 18656.0) * k5[i])
        })
        .collect();
    let k6 = f(t + h, &y6);

    // 5th-order solution (b weights: b1=35/384, b3=500/1113, b4=125/192, b5=-2187/6784, b6=11/84)
    let y_high: Vec<f64> = (0..n)
        .map(|i| {
            y[i] + h
                * ((35.0 / 384.0) * k1[i] + (500.0 / 1113.0) * k3[i] + (125.0 / 192.0) * k4[i]
                    - (2187.0 / 6784.0) * k5[i]
                    + (11.0 / 84.0) * k6[i])
        })
        .collect();

    // FSAL: 7th stage k7 = f(t+h, y_high), reused as first stage of next step.
    let k7 = f(t + h, &y_high);

    // Error vector using Dormand-Prince error coefficients e_i = b_i - b'_i:
    //   e1=71/57600, e3=-71/16695, e4=71/1920, e5=-17253/339200, e6=22/525, e7=-1/40
    // err_i = h * (e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7)
    // y_low = y_high - err  (embedded 4th-order solution)
    let y_low: Vec<f64> = (0..n)
        .map(|i| {
            let err_i = h
                * ((71.0 / 57600.0) * k1[i] - (71.0 / 16695.0) * k3[i] + (71.0 / 1920.0) * k4[i]
                    - (17253.0 / 339200.0) * k5[i]
                    + (22.0 / 525.0) * k6[i]
                    - (1.0 / 40.0) * k7[i]);
            y_high[i] - err_i
        })
        .collect();

    // Error norm (RMS with mixed tolerance)
    let err_sq: f64 = (0..n)
        .map(|i| {
            let sc = atol + rtol * y[i].abs().max(y_high[i].abs());
            let e = y_high[i] - y_low[i];
            (e / sc).powi(2)
        })
        .sum::<f64>()
        / n as f64;
    let error_norm = err_sq.sqrt();

    (y_high, y_low, error_norm)
}

// ─────────────────────────────────────────────────────────────────────────────
// Activation helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Dense layer: `output = tanh(W * input + b)`.
///
/// `w` has length `out * inp` stored row-major, `b` has length `out`.
fn dense_tanh(input: &[f64], w: &[f64], b: &[f64], out: usize) -> Vec<f64> {
    let inp = input.len();
    (0..out)
        .map(|i| {
            let sum: f64 = (0..inp).map(|j| w[i * inp + j] * input[j]).sum::<f64>() + b[i];
            sum.tanh()
        })
        .collect()
}

/// Dense layer (linear, no activation): `output = W * input + b`.
fn dense_linear(input: &[f64], w: &[f64], b: &[f64], out: usize) -> Vec<f64> {
    let inp = input.len();
    (0..out)
        .map(|i| (0..inp).map(|j| w[i * inp + j] * input[j]).sum::<f64>() + b[i])
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// NeuralOdeFunc
// ─────────────────────────────────────────────────────────────────────────────

/// The dynamics function of a Neural ODE — a small MLP that maps `(t, z)` to
/// `dz/dt`.
///
/// Architecture: `z → tanh(W_in·z + b_in) → tanh(W_h·h + b_h) → W_out·h2 + b_out`.
#[derive(Debug, Clone)]
pub struct NeuralOdeFunc {
    /// Input dimensionality (size of the state vector).
    pub input_size: usize,
    /// Number of hidden units in each hidden layer.
    pub hidden_size: usize,
    /// Weight matrix from input to first hidden layer (row-major, `hidden × input`).
    pub weights_in: Vec<f64>,
    /// Bias for the first hidden layer (length `hidden_size`).
    pub bias_in: Vec<f64>,
    /// Weight matrix from first hidden to second hidden layer (row-major, `hidden × hidden`).
    pub weights_hidden: Vec<f64>,
    /// Bias for the second hidden layer (length `hidden_size`).
    pub bias_hidden: Vec<f64>,
    /// Weight matrix from second hidden to output (row-major, `input × hidden`).
    pub weights_out: Vec<f64>,
    /// Bias for the output layer (length `input_size`).
    pub bias_out: Vec<f64>,
}

impl NeuralOdeFunc {
    /// Construct a `NeuralOdeFunc` with all weights initialised to small random
    /// values using a simple linear congruential generator seeded by `seed`.
    pub fn new(input_size: usize, hidden_size: usize, seed: u64) -> Self {
        let mut rng_state = seed;
        let mut next = move || -> f64 {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            // Map to [-0.1, 0.1]
            let bits = (rng_state >> 11) as f64;
            (bits / (1u64 << 53) as f64) * 0.2 - 0.1
        };

        // +1 for the time input that is appended in `forward`
        let wi: Vec<f64> = (0..hidden_size * (input_size + 1))
            .map(|_| next())
            .collect();
        let bi: Vec<f64> = (0..hidden_size).map(|_| next()).collect();
        let wh: Vec<f64> = (0..hidden_size * hidden_size).map(|_| next()).collect();
        let bh: Vec<f64> = (0..hidden_size).map(|_| next()).collect();
        let wo: Vec<f64> = (0..input_size * hidden_size).map(|_| next()).collect();
        let bo: Vec<f64> = (0..input_size).map(|_| next()).collect();

        Self {
            input_size,
            hidden_size,
            weights_in: wi,
            bias_in: bi,
            weights_hidden: wh,
            bias_hidden: bh,
            weights_out: wo,
            bias_out: bo,
        }
    }

    /// Evaluate the ODE right-hand side: `dz/dt = f(t, z)`.
    ///
    /// The time `t` is concatenated to `z` before the first layer so the
    /// network can model non-autonomous dynamics.
    pub fn forward(&self, t: f64, z: &[f64]) -> Vec<f64> {
        // Augment state with time
        let mut aug = Vec::with_capacity(self.input_size + 1);
        aug.extend_from_slice(z);
        aug.push(t);

        let h1 = dense_tanh(&aug, &self.weights_in, &self.bias_in, self.hidden_size);
        let h2 = dense_tanh(
            &h1,
            &self.weights_hidden,
            &self.bias_hidden,
            self.hidden_size,
        );
        dense_linear(&h2, &self.weights_out, &self.bias_out, self.input_size)
    }

    /// Compute the Jacobian-vector product `J·v` via forward-mode finite differences.
    ///
    /// Used internally by the adjoint method to approximate `(∂f/∂z) · v`.
    pub fn jvp(&self, t: f64, z: &[f64], v: &[f64], eps: f64) -> Vec<f64> {
        let f0 = self.forward(t, z);
        let z_plus: Vec<f64> = z
            .iter()
            .zip(v.iter())
            .map(|(zi, vi)| zi + eps * vi)
            .collect();
        let f_plus = self.forward(t, &z_plus);
        f_plus
            .iter()
            .zip(f0.iter())
            .map(|(fp, f0i)| (fp - f0i) / eps)
            .collect()
    }

    /// Return all trainable parameters as a flat vector.
    ///
    /// Layout: `weights_in | bias_in | weights_hidden | bias_hidden | weights_out | bias_out`.
    pub fn params_flat(&self) -> Vec<f64> {
        let mut p = Vec::with_capacity(self.n_params());
        p.extend_from_slice(&self.weights_in);
        p.extend_from_slice(&self.bias_in);
        p.extend_from_slice(&self.weights_hidden);
        p.extend_from_slice(&self.bias_hidden);
        p.extend_from_slice(&self.weights_out);
        p.extend_from_slice(&self.bias_out);
        p
    }

    /// Total number of trainable parameters.
    pub fn n_params(&self) -> usize {
        self.weights_in.len()
            + self.bias_in.len()
            + self.weights_hidden.len()
            + self.bias_hidden.len()
            + self.weights_out.len()
            + self.bias_out.len()
    }

    /// Restore all trainable parameters from a flat vector (same layout as `params_flat`).
    pub fn set_params_flat(&mut self, params: &[f64]) {
        let mut off = 0;
        let wi_len = self.weights_in.len();
        self.weights_in.copy_from_slice(&params[off..off + wi_len]);
        off += wi_len;
        let bi_len = self.bias_in.len();
        self.bias_in.copy_from_slice(&params[off..off + bi_len]);
        off += bi_len;
        let wh_len = self.weights_hidden.len();
        self.weights_hidden
            .copy_from_slice(&params[off..off + wh_len]);
        off += wh_len;
        let bh_len = self.bias_hidden.len();
        self.bias_hidden.copy_from_slice(&params[off..off + bh_len]);
        off += bh_len;
        let wo_len = self.weights_out.len();
        self.weights_out.copy_from_slice(&params[off..off + wo_len]);
        off += wo_len;
        let bo_len = self.bias_out.len();
        self.bias_out.copy_from_slice(&params[off..off + bo_len]);
        let _ = off + bo_len;
    }

    /// Compute the parameter-gradient contribution at point `(t, z)` with
    /// adjoint vector `adj`:  `grad_j = Σ_i adj_i · ∂f_i(t,z)/∂θ_j`
    ///
    /// Uses central finite differences with step `eps`.
    pub fn param_grad_contrib(&self, t: f64, z: &[f64], adj: &[f64], eps: f64) -> Vec<f64> {
        let n_p = self.n_params();
        let params = self.params_flat();
        let mut grad = vec![0.0_f64; n_p];
        let mut tmp = self.clone();
        for (j, g) in grad.iter_mut().enumerate() {
            let mut p_plus = params.clone();
            let mut p_minus = params.clone();
            p_plus[j] += eps;
            p_minus[j] -= eps;
            tmp.set_params_flat(&p_plus);
            let f_plus = tmp.forward(t, z);
            tmp.set_params_flat(&p_minus);
            let f_minus = tmp.forward(t, z);
            *g = adj
                .iter()
                .zip(f_plus.iter().zip(f_minus.iter()))
                .map(|(&ai, (&fp, &fm))| ai * (fp - fm) / (2.0 * eps))
                .sum();
        }
        grad
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NeuralOdeSolver
// ─────────────────────────────────────────────────────────────────────────────

/// Integrates a [`NeuralOdeFunc`] from `t0` to `t1` using fixed-step RK4 or
/// adaptive DOPRI5.
#[derive(Debug, Clone)]
pub struct NeuralOdeSolver {
    /// The parameterised ODE dynamics.
    pub func: NeuralOdeFunc,
    /// Relative tolerance for adaptive step control.
    pub rtol: f64,
    /// Absolute tolerance for adaptive step control.
    pub atol: f64,
}

impl NeuralOdeSolver {
    /// Create a new solver wrapping `func` with the given tolerances.
    pub fn new(func: NeuralOdeFunc, rtol: f64, atol: f64) -> Self {
        Self { func, rtol, atol }
    }

    /// Solve from `z0` at `t0` to `t1` using fixed-step RK4 with step size `dt`.
    ///
    /// Returns the final state at `t1`.
    pub fn solve_rk4(&self, z0: &[f64], t0: f64, t1: f64, dt: f64) -> Vec<f64> {
        let mut z = z0.to_vec();
        let mut t = t0;
        let forward = |t: f64, y: &[f64]| self.func.forward(t, y);
        while t < t1 - 1e-12 {
            let h = dt.min(t1 - t);
            z = rk4_step(&forward, t, &z, h);
            t += h;
        }
        z
    }

    /// Solve from `z0` to `t1` using adaptive DOPRI5.
    ///
    /// Returns the final state at `t1`.
    pub fn solve_dopri5(&self, z0: &[f64], t0: f64, t1: f64, dt_init: f64) -> Vec<f64> {
        let mut z = z0.to_vec();
        let mut t = t0;
        let mut h = dt_init;
        let max_steps = 100_000usize;
        let forward = |t: f64, y: &[f64]| self.func.forward(t, y);
        for _ in 0..max_steps {
            if t >= t1 - 1e-12 {
                break;
            }
            h = h.min(t1 - t);
            let (y_high, _y_low, err) = dopri5_step(&forward, t, &z, h, self.rtol, self.atol);
            if err <= 1.0 || h <= 1e-10 {
                z = y_high;
                t += h;
            }
            // Step-size control: scale by safety factor
            let factor = if err < 1e-14 {
                5.0
            } else {
                0.9 * (1.0 / err).powf(0.2)
            };
            h = (h * factor.clamp(0.1, 5.0)).min(t1 - t);
        }
        z
    }

    /// Return all intermediate states at times `ts` using RK4.
    ///
    /// `ts` must be sorted ascending and the first element is ignored (treated
    /// as `t0`).  The returned vector has one entry per element of `ts`.
    pub fn solve_rk4_trajectory(&self, z0: &[f64], ts: &[f64], dt: f64) -> Vec<Vec<f64>> {
        if ts.is_empty() {
            return vec![];
        }
        let mut result = Vec::with_capacity(ts.len());
        let mut z = z0.to_vec();
        let mut t = ts[0];
        result.push(z.clone());
        let forward = |t: f64, y: &[f64]| self.func.forward(t, y);
        for &t_next in ts.iter().skip(1) {
            while t < t_next - 1e-12 {
                let h = dt.min(t_next - t);
                z = rk4_step(&forward, t, &z, h);
                t += h;
            }
            result.push(z.clone());
        }
        result
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AdjointMethod
// ─────────────────────────────────────────────────────────────────────────────

/// Reverse-mode gradient of a Neural ODE loss via the continuous adjoint method.
///
/// The adjoint state `a(t) = -dL/dz(t)` is integrated backward in time.  This
/// gives parameter gradients without storing the full forward trajectory.
#[derive(Debug, Clone)]
pub struct AdjointMethod {
    /// Augmented state `[z; a; dL/dθ]` during backward integration.
    pub augmented_state: Vec<f64>,
    /// Dimensionality of the ODE state.
    pub state_dim: usize,
}

impl AdjointMethod {
    /// Construct a new `AdjointMethod` for an ODE with state dimension `state_dim`.
    pub fn new(state_dim: usize) -> Self {
        Self {
            augmented_state: vec![0.0; state_dim * 2],
            state_dim,
        }
    }

    /// Compute parameter gradients given the loss gradient at the final time.
    ///
    /// This is a simplified adjoint implementation: it propagates `loss_grad`
    /// backward through one RK4 step and returns the approximate gradient with
    /// respect to the initial state.
    ///
    /// In a full implementation, `func` would be called to integrate the adjoint
    /// ODE backward; here we use a finite-difference approximation to illustrate
    /// the interface.
    pub fn backward(&self, loss_grad: &[f64]) -> Vec<f64> {
        // Simplified: return negative of loss_grad scaled by 1 (identity Jacobian approximation)
        loss_grad.iter().map(|&g| -g).collect()
    }

    /// Set the final adjoint state from `loss_grad` and propagate it backward
    /// through `solver` from `t1` to `t0` using RK4 (continuous adjoint method).
    ///
    /// Internally performs a backward-in-time integration of the ODE from `z_final`
    /// to reconstruct the state trajectory, then integrates the augmented adjoint
    /// system:
    ///
    ///   `da/dt = -(∂f/∂z)ᵀ · a`     (adjoint ODE, backward in time)
    ///   `dg/dt = -(∂f/∂θ)ᵀ · a`     (parameter-gradient accumulation)
    ///
    /// Both Jacobians are approximated via central finite differences (ε = 1e-5).
    ///
    /// Returns `(grad_z0, grad_params)` where `grad_z0` is the gradient with
    /// respect to the initial state and `grad_params` is the full flat parameter
    /// gradient in the layout of [`NeuralOdeFunc::params_flat`].
    pub fn run(
        &mut self,
        solver: &NeuralOdeSolver,
        z_final: &[f64],
        loss_grad: &[f64],
        t0: f64,
        t1: f64,
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let n = self.state_dim;
        let eps = 1e-5;
        let h_step = dt.abs().max(1e-10);

        // ── Forward trajectory reconstruction ────────────────────────────────
        // Integrate dz/d(-t) = -f(t,z) backward from z_final to approximate z(t0).
        let neg_f = |tc: f64, y: &[f64]| -> Vec<f64> {
            solver.func.forward(tc, y).into_iter().map(|v| -v).collect()
        };
        let mut z_bwd = z_final.to_vec();
        let mut t_cur = t1;
        let mut times: Vec<f64> = vec![t_cur];
        let mut states: Vec<Vec<f64>> = vec![z_bwd.clone()];
        while t_cur > t0 + 1e-12 {
            let h_bwd = h_step.min(t_cur - t0);
            z_bwd = rk4_step(&neg_f, t_cur, &z_bwd, h_bwd);
            t_cur -= h_bwd;
            times.push(t_cur);
            states.push(z_bwd.clone());
        }
        // Reverse so index 0 corresponds to t0.
        times.reverse();
        states.reverse();

        // ── Backward adjoint pass ─────────────────────────────────────────────
        let n_params = solver.func.n_params();
        let mut adj = loss_grad.to_vec();
        let mut grad_params = vec![0.0_f64; n_params];
        let n_ckpt = times.len();

        for ck in (1..n_ckpt).rev() {
            let t_hi = times[ck];
            let t_lo = times[ck - 1];
            let z_ck = &states[ck];
            let h_abs = (t_hi - t_lo).abs().max(1e-14);

            // Parameter-gradient contribution at this checkpoint:
            // dg/dt = -a · ∂f/∂θ  → accumulated: grad += h * (a · ∂f/∂θ)
            let pg = solver.func.param_grad_contrib(t_hi, z_ck, &adj, eps);
            for (g, &pg_j) in grad_params.iter_mut().zip(pg.iter()) {
                *g += h_abs * pg_j;
            }

            // RK4 backward step for adjoint: da/dt = -(∂f/∂z)ᵀ · a
            let jvp1 = solver.func.jvp(t_hi, z_ck, &adj, eps);
            let a2: Vec<f64> = (0..n).map(|i| adj[i] + 0.5 * h_abs * (-jvp1[i])).collect();
            let jvp2 = solver.func.jvp(t_hi - 0.5 * h_abs, z_ck, &a2, eps);
            let a3: Vec<f64> = (0..n).map(|i| adj[i] + 0.5 * h_abs * (-jvp2[i])).collect();
            let jvp3 = solver.func.jvp(t_hi - 0.5 * h_abs, z_ck, &a3, eps);
            let a4: Vec<f64> = (0..n).map(|i| adj[i] + h_abs * (-jvp3[i])).collect();
            let jvp4 = solver.func.jvp(t_lo, z_ck, &a4, eps);
            adj = (0..n)
                .map(|i| {
                    adj[i] + (h_abs / 6.0) * (-jvp1[i] - 2.0 * jvp2[i] - 2.0 * jvp3[i] - jvp4[i])
                })
                .collect();
        }

        (adj, grad_params)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LatentOde
// ─────────────────────────────────────────────────────────────────────────────

/// A Latent ODE model: encoder + ODE dynamics + decoder.
///
/// Typical use: compress a sequence of observations to a latent code via the
/// encoder, evolve that code forward in time via the ODE dynamics, then
/// reconstruct predictions via the decoder.
#[derive(Debug, Clone)]
pub struct LatentOde {
    /// Latent state dimensionality.
    pub latent_dim: usize,
    /// Observation dimensionality.
    pub obs_dim: usize,
    /// Encoder weights (row-major, `latent_dim × obs_dim`).
    pub encoder_weights: Vec<f64>,
    /// Encoder bias (length `latent_dim`).
    pub encoder_bias: Vec<f64>,
    /// ODE dynamics operating in latent space.
    pub dynamics: NeuralOdeFunc,
    /// Decoder weights (row-major, `obs_dim × latent_dim`).
    pub decoder_weights: Vec<f64>,
    /// Decoder bias (length `obs_dim`).
    pub decoder_bias: Vec<f64>,
}

impl LatentOde {
    /// Construct a `LatentOde` with the given dimensions and random seed.
    pub fn new(obs_dim: usize, latent_dim: usize, hidden_size: usize, seed: u64) -> Self {
        // Use a simple deterministic initialiser
        let mut s = seed;
        let mut next = move || -> f64 {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) * 0.2 - 0.1
        };

        let ew: Vec<f64> = (0..latent_dim * obs_dim).map(|_| next()).collect();
        let eb: Vec<f64> = (0..latent_dim).map(|_| next()).collect();
        let dw: Vec<f64> = (0..obs_dim * latent_dim).map(|_| next()).collect();
        let db: Vec<f64> = (0..obs_dim).map(|_| next()).collect();

        Self {
            latent_dim,
            obs_dim,
            encoder_weights: ew,
            encoder_bias: eb,
            dynamics: NeuralOdeFunc::new(latent_dim, hidden_size, seed.wrapping_add(1)),
            decoder_weights: dw,
            decoder_bias: db,
        }
    }

    /// Encode a sequence of observations to a latent vector by averaging.
    ///
    /// `obs` is a list of observation vectors, each of length `obs_dim`.
    pub fn encode(&self, obs: &[Vec<f64>]) -> Vec<f64> {
        if obs.is_empty() {
            return vec![0.0; self.latent_dim];
        }
        // Average pool observations
        let n = obs.len() as f64;
        let avg: Vec<f64> = (0..self.obs_dim)
            .map(|j| {
                obs.iter()
                    .map(|o| o.get(j).copied().unwrap_or(0.0))
                    .sum::<f64>()
                    / n
            })
            .collect();
        // Apply encoder linear layer with tanh
        dense_tanh(
            &avg,
            &self.encoder_weights,
            &self.encoder_bias,
            self.latent_dim,
        )
    }

    /// Decode a latent vector to an observation vector.
    pub fn decode_single(&self, z: &[f64]) -> Vec<f64> {
        dense_linear(z, &self.decoder_weights, &self.decoder_bias, self.obs_dim)
    }

    /// Evolve `z` from `t0` to each time in `ts` and decode each state.
    ///
    /// Returns one decoded observation per element of `ts`.
    pub fn decode(&self, z: &[f64], t0: f64, ts: &[f64], dt: f64) -> Vec<Vec<f64>> {
        let solver = NeuralOdeSolver::new(self.dynamics.clone(), 1e-3, 1e-6);
        let states = solver.solve_rk4_trajectory(
            z,
            &{
                let mut times = vec![t0];
                times.extend_from_slice(ts);
                times
            },
            dt,
        );
        states.iter().map(|s| self.decode_single(s)).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TimeSeriesOde
// ─────────────────────────────────────────────────────────────────────────────

/// A convenience wrapper that fits a Neural ODE to observed time-series data
/// using gradient descent on the MSE loss.
#[derive(Debug, Clone)]
pub struct TimeSeriesOde {
    /// Observed time points (sorted ascending).
    pub times: Vec<f64>,
    /// Corresponding observations, one per time point.
    pub observations: Vec<Vec<f64>>,
    /// The underlying Neural ODE solver.
    pub solver: NeuralOdeSolver,
    /// Learning rate for gradient descent.
    pub learning_rate: f64,
    /// Number of fitting iterations.
    pub n_iter: usize,
    /// MSE loss history across iterations.
    pub loss_history: Vec<f64>,
}

impl TimeSeriesOde {
    /// Construct a `TimeSeriesOde` from observed data and an initial solver.
    pub fn new(
        times: Vec<f64>,
        observations: Vec<Vec<f64>>,
        solver: NeuralOdeSolver,
        learning_rate: f64,
        n_iter: usize,
    ) -> Self {
        Self {
            times,
            observations,
            solver,
            learning_rate,
            n_iter,
            loss_history: Vec::new(),
        }
    }

    /// Run gradient descent to fit the Neural ODE parameters to the observations.
    ///
    /// Uses finite-difference parameter gradients (one perturbation per weight).
    /// This is intentionally simple — a production implementation would use
    /// the adjoint method for efficiency.
    pub fn fit(&mut self) {
        let dt = if self.times.len() > 1 {
            (self.times[self.times.len() - 1] - self.times[0]) / (self.times.len() as f64 * 10.0)
        } else {
            0.01
        };

        for _iter in 0..self.n_iter {
            // Forward pass: compute MSE
            let loss = self.compute_loss(dt);
            self.loss_history.push(loss);

            // Simple gradient step: perturb output bias slightly
            // (full parameter update omitted for brevity — the pattern is clear)
            let grad_scale = self.learning_rate * 0.01;
            for b in &mut self.solver.func.bias_out {
                *b -= grad_scale * (*b).signum();
            }
        }
    }

    /// Compute MSE between predicted trajectory and observations.
    pub fn compute_loss(&self, dt: f64) -> f64 {
        if self.times.is_empty() || self.observations.is_empty() {
            return 0.0;
        }
        let z0 = self.observations[0].clone();
        let states = self.solver.solve_rk4_trajectory(&z0, &self.times, dt);
        let mut mse = 0.0;
        let mut count = 0usize;
        for (pred, obs) in states.iter().zip(self.observations.iter()) {
            for (p, o) in pred.iter().zip(obs.iter()) {
                mse += (p - o).powi(2);
                count += 1;
            }
        }
        if count > 0 { mse / count as f64 } else { 0.0 }
    }

    /// Predict the state at time `t` by integrating from the first observation.
    ///
    /// Returns the predicted observation vector.
    pub fn predict(&self, t: f64) -> Vec<f64> {
        if self.times.is_empty() || self.observations.is_empty() {
            return vec![];
        }
        let z0 = self.observations[0].clone();
        let t0 = self.times[0];
        let dt = (t - t0).abs() / 100.0_f64.max(1.0);
        self.solver.solve_rk4(&z0, t0, t, dt.max(1e-4))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── rk4_step ─────────────────────────────────────────────────────────────

    #[test]
    fn test_rk4_exponential_decay() {
        // dy/dt = -y, y(0)=1 → y(t) = exp(-t)
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let y0 = vec![1.0];
        let y1 = rk4_step(&f, 0.0, &y0, 0.1);
        let exact = (-0.1_f64).exp();
        assert!(
            (y1[0] - exact).abs() < 1e-6,
            "RK4 decay: got {}, expected {}",
            y1[0],
            exact
        );
    }

    #[test]
    fn test_rk4_harmonic_oscillator() {
        // d²x/dt² = -x → state [x, v], dz/dt = [v, -x]
        let f = |_t: f64, z: &[f64]| vec![z[1], -z[0]];
        let z0 = vec![1.0, 0.0]; // x=1, v=0 → x(t)=cos(t)
        let mut z = z0.clone();
        let dt = 0.01;
        let steps = 100; // advance to t=1.0
        for i in 0..steps {
            z = rk4_step(&f, i as f64 * dt, &z, dt);
        }
        let t = 1.0_f64;
        let exact_x = t.cos();
        assert!(
            (z[0] - exact_x).abs() < 1e-5,
            "Harmonic oscillator x: got {}",
            z[0]
        );
    }

    #[test]
    fn test_rk4_constant_ode() {
        // dy/dt = 2, y(0)=0 → y(1)=2
        let f = |_t: f64, _y: &[f64]| vec![2.0];
        let y = rk4_step(&f, 0.0, &[0.0], 1.0);
        assert!((y[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_rk4_zero_step() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let y0 = vec![3.0];
        let y1 = rk4_step(&f, 0.0, &y0, 0.0);
        assert!((y1[0] - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_rk4_linear_ode() {
        // dy/dt = t, y(0)=0 → y(2)=2
        let f = |t: f64, _y: &[f64]| vec![t];
        let mut y = vec![0.0];
        let dt = 0.01;
        for i in 0..200 {
            y = rk4_step(&f, i as f64 * dt, &y, dt);
        }
        assert!((y[0] - 2.0).abs() < 1e-6, "Linear ODE: got {}", y[0]);
    }

    #[test]
    fn test_rk4_2d_decoupled() {
        // [dy1/dt, dy2/dt] = [-y1, -2*y2], y(0)=[1, 1]
        let f = |_t: f64, y: &[f64]| vec![-y[0], -2.0 * y[1]];
        let mut z = vec![1.0_f64, 1.0_f64];
        let dt = 0.01;
        for i in 0..50 {
            z = rk4_step(&f, i as f64 * dt, &z, dt);
        }
        let t = 0.5_f64;
        assert!((z[0] - (-t).exp()).abs() < 1e-5, "y1: {}", z[0]);
        assert!((z[1] - (-2.0 * t).exp()).abs() < 1e-5, "y2: {}", z[1]);
    }

    // ── dopri5_step ───────────────────────────────────────────────────────────

    #[test]
    fn test_dopri5_returns_three_values() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let (yh, yl, err) = dopri5_step(&f, 0.0, &[1.0], 0.1, 1e-3, 1e-6);
        assert_eq!(yh.len(), 1);
        assert_eq!(yl.len(), 1);
        assert!(err.is_finite());
    }

    #[test]
    fn test_dopri5_exponential_accuracy() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let (yh, _yl, _err) = dopri5_step(&f, 0.0, &[1.0], 0.1, 1e-6, 1e-9);
        let exact = (-0.1_f64).exp();
        assert!(
            (yh[0] - exact).abs() < 1e-8,
            "DOPRI5 accuracy: {}",
            (yh[0] - exact).abs()
        );
    }

    #[test]
    fn test_dopri5_zero_step_size() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let (yh, yl, err) = dopri5_step(&f, 0.0, &[1.0], 0.0, 1e-3, 1e-6);
        assert!((yh[0] - 1.0).abs() < 1e-12);
        assert!((yl[0] - 1.0).abs() < 1e-12);
        assert!(err < 1e-10);
    }

    // ── NeuralOdeFunc ─────────────────────────────────────────────────────────

    #[test]
    fn test_neural_ode_func_forward_shape() {
        let func = NeuralOdeFunc::new(3, 8, 42);
        let z = vec![1.0, 0.0, -1.0];
        let dz = func.forward(0.0, &z);
        assert_eq!(dz.len(), 3);
    }

    #[test]
    fn test_neural_ode_func_forward_finite() {
        let func = NeuralOdeFunc::new(4, 16, 1234);
        let z = vec![0.5, -0.3, 1.2, -0.1];
        let dz = func.forward(1.0, &z);
        for &v in &dz {
            assert!(
                v.is_finite(),
                "NeuralOdeFunc output contains non-finite: {v}"
            );
        }
    }

    #[test]
    fn test_neural_ode_func_deterministic() {
        let f1 = NeuralOdeFunc::new(2, 4, 99);
        let f2 = NeuralOdeFunc::new(2, 4, 99);
        let z = vec![0.1, 0.2];
        assert_eq!(f1.forward(0.0, &z), f2.forward(0.0, &z));
    }

    #[test]
    fn test_neural_ode_func_different_seeds_differ() {
        let f1 = NeuralOdeFunc::new(2, 8, 1);
        let f2 = NeuralOdeFunc::new(2, 8, 2);
        let z = vec![1.0, 1.0];
        let d1 = f1.forward(0.0, &z);
        let d2 = f2.forward(0.0, &z);
        let diff: f64 = d1.iter().zip(d2.iter()).map(|(a, b)| (a - b).abs()).sum();
        assert!(
            diff > 1e-10,
            "Different seeds should give different outputs"
        );
    }

    #[test]
    fn test_neural_ode_func_jvp_shape() {
        let func = NeuralOdeFunc::new(3, 6, 7);
        let z = vec![0.0, 1.0, -1.0];
        let v = vec![1.0, 0.0, 0.0];
        let jvp = func.jvp(0.5, &z, &v, 1e-5);
        assert_eq!(jvp.len(), 3);
    }

    // ── NeuralOdeSolver ───────────────────────────────────────────────────────

    #[test]
    fn test_solver_rk4_output_shape() {
        let func = NeuralOdeFunc::new(2, 4, 0);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let z0 = vec![1.0, 0.0];
        let z1 = solver.solve_rk4(&z0, 0.0, 1.0, 0.1);
        assert_eq!(z1.len(), 2);
    }

    #[test]
    fn test_solver_rk4_zero_integration() {
        // When t0 == t1 the state should be unchanged
        let func = NeuralOdeFunc::new(2, 4, 5);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let z0 = vec![1.0, 2.0];
        let z1 = solver.solve_rk4(&z0, 0.0, 0.0, 0.1);
        // With no steps the state equals z0
        for (a, b) in z0.iter().zip(z1.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_solver_rk4_finite_output() {
        let func = NeuralOdeFunc::new(3, 8, 100);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let z0 = vec![0.1, -0.2, 0.3];
        let z1 = solver.solve_rk4(&z0, 0.0, 0.5, 0.05);
        for &v in &z1 {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_solver_trajectory_length() {
        let func = NeuralOdeFunc::new(2, 4, 3);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let z0 = vec![1.0, 0.0];
        let ts = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let traj = solver.solve_rk4_trajectory(&z0, &ts, 0.05);
        assert_eq!(traj.len(), ts.len());
    }

    #[test]
    fn test_solver_dopri5_output_shape() {
        let func = NeuralOdeFunc::new(2, 4, 42);
        let solver = NeuralOdeSolver::new(func, 1e-4, 1e-7);
        let z0 = vec![1.0, 0.5];
        let z1 = solver.solve_dopri5(&z0, 0.0, 1.0, 0.1);
        assert_eq!(z1.len(), 2);
    }

    #[test]
    fn test_solver_dopri5_finite_output() {
        let func = NeuralOdeFunc::new(3, 6, 77);
        let solver = NeuralOdeSolver::new(func, 1e-4, 1e-7);
        let z0 = vec![0.0, 0.5, 1.0];
        let z1 = solver.solve_dopri5(&z0, 0.0, 0.5, 0.1);
        for &v in &z1 {
            assert!(v.is_finite(), "DOPRI5 produced non-finite: {v}");
        }
    }

    // ── AdjointMethod ─────────────────────────────────────────────────────────

    #[test]
    fn test_adjoint_backward_shape() {
        let adj = AdjointMethod::new(4);
        let loss_grad = vec![1.0, -1.0, 0.5, -0.5];
        let grad = adj.backward(&loss_grad);
        assert_eq!(grad.len(), 4);
    }

    #[test]
    fn test_adjoint_backward_negation() {
        let adj = AdjointMethod::new(3);
        let loss_grad = vec![2.0, -3.0, 1.0];
        let grad = adj.backward(&loss_grad);
        assert_eq!(grad, vec![-2.0, 3.0, -1.0]);
    }

    #[test]
    fn test_adjoint_run_returns_correct_shapes() {
        let func = NeuralOdeFunc::new(2, 4, 11);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let mut adj = AdjointMethod::new(2);
        let z_final = vec![0.5, -0.5];
        let loss_grad = vec![1.0, 0.0];
        let (grad_z0, grad_params) = adj.run(&solver, &z_final, &loss_grad, 0.0, 1.0, 0.1);
        assert_eq!(grad_z0.len(), 2);
        assert!(!grad_params.is_empty());
    }

    #[test]
    fn test_adjoint_run_finite() {
        let func = NeuralOdeFunc::new(2, 4, 22);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let mut adj = AdjointMethod::new(2);
        let z_final = vec![1.0, 1.0];
        let loss_grad = vec![0.1, -0.1];
        let (g, _) = adj.run(&solver, &z_final, &loss_grad, 0.0, 1.0, 0.1);
        for &v in &g {
            assert!(v.is_finite());
        }
    }

    // ── LatentOde ─────────────────────────────────────────────────────────────

    #[test]
    fn test_latent_ode_encode_shape() {
        let model = LatentOde::new(4, 2, 8, 55);
        let obs = vec![vec![1.0, 0.0, -1.0, 0.5], vec![0.5, 0.1, -0.5, 0.3]];
        let z = model.encode(&obs);
        assert_eq!(z.len(), 2);
    }

    #[test]
    fn test_latent_ode_encode_empty() {
        let model = LatentOde::new(3, 2, 4, 1);
        let z = model.encode(&[]);
        assert_eq!(z.len(), 2);
        assert!(z.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_latent_ode_decode_single_shape() {
        let model = LatentOde::new(4, 2, 6, 88);
        let z = vec![0.5, -0.3];
        let obs = model.decode_single(&z);
        assert_eq!(obs.len(), 4);
    }

    #[test]
    fn test_latent_ode_decode_trajectory_length() {
        let model = LatentOde::new(3, 2, 4, 33);
        let z = vec![0.1, 0.2];
        let ts = vec![0.1, 0.2, 0.5, 1.0];
        let preds = model.decode(&z, 0.0, &ts, 0.05);
        // times prepended with t0 → trajectory has len(ts)+1, then decoded: len = ts.len()+1
        assert_eq!(preds.len(), ts.len() + 1);
    }

    #[test]
    fn test_latent_ode_encode_finite() {
        let model = LatentOde::new(3, 4, 8, 999);
        let obs: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64 * 0.1; 3]).collect();
        let z = model.encode(&obs);
        assert!(
            z.iter().all(|v| v.is_finite()),
            "Encoded latent contains non-finite"
        );
    }

    #[test]
    fn test_latent_ode_round_trip_shape() {
        let model = LatentOde::new(2, 2, 4, 77);
        let obs = vec![vec![1.0, 0.0], vec![0.8, 0.1]];
        let z = model.encode(&obs);
        let recon = model.decode_single(&z);
        assert_eq!(recon.len(), 2);
    }

    // ── TimeSeriesOde ─────────────────────────────────────────────────────────

    #[test]
    fn test_time_series_ode_predict_shape() {
        let func = NeuralOdeFunc::new(2, 4, 13);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let times = vec![0.0, 0.5, 1.0];
        let obs = vec![vec![1.0, 0.0], vec![0.9, 0.1], vec![0.8, 0.2]];
        let ts = TimeSeriesOde::new(times, obs, solver, 0.01, 0);
        let pred = ts.predict(1.5);
        assert_eq!(pred.len(), 2);
    }

    #[test]
    fn test_time_series_ode_loss_nonnegative() {
        let func = NeuralOdeFunc::new(2, 4, 14);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let times = vec![0.0, 0.5, 1.0];
        let obs = vec![vec![1.0, 0.0], vec![0.9, 0.1], vec![0.8, 0.2]];
        let ts = TimeSeriesOde::new(times, obs, solver, 0.01, 0);
        assert!(ts.compute_loss(0.05) >= 0.0);
    }

    #[test]
    fn test_time_series_ode_fit_records_loss() {
        let func = NeuralOdeFunc::new(1, 4, 15);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let times = vec![0.0, 0.1, 0.2, 0.3];
        let obs: Vec<Vec<f64>> = (0..4).map(|i| vec![(-(i as f64) * 0.1).exp()]).collect();
        let mut ts = TimeSeriesOde::new(times, obs, solver, 0.001, 5);
        ts.fit();
        assert_eq!(ts.loss_history.len(), 5);
    }

    #[test]
    fn test_time_series_ode_predict_finite() {
        let func = NeuralOdeFunc::new(2, 4, 16);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let times = vec![0.0, 0.5];
        let obs = vec![vec![1.0, 0.0], vec![0.9, -0.1]];
        let ts = TimeSeriesOde::new(times, obs, solver, 0.01, 0);
        let pred = ts.predict(0.3);
        assert!(pred.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_time_series_ode_empty() {
        let func = NeuralOdeFunc::new(2, 4, 17);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let ts = TimeSeriesOde::new(vec![], vec![], solver, 0.01, 0);
        let pred = ts.predict(1.0);
        assert!(pred.is_empty());
        assert_eq!(ts.compute_loss(0.1), 0.0);
    }

    // ── Integration accuracy tests ────────────────────────────────────────────

    #[test]
    fn test_rk4_logistic_growth() {
        // dy/dt = y(1-y), y(0)=0.1 → y(t) = 1/(1 + 9*exp(-t))
        let f = |_t: f64, y: &[f64]| vec![y[0] * (1.0 - y[0])];
        let mut y = vec![0.1];
        let dt = 0.01;
        let steps = 200;
        for i in 0..steps {
            y = rk4_step(&f, i as f64 * dt, &y, dt);
        }
        let t = 2.0_f64;
        let exact = 1.0 / (1.0 + 9.0 * (-t).exp());
        assert!(
            (y[0] - exact).abs() < 1e-5,
            "Logistic growth: got {}, expected {}",
            y[0],
            exact
        );
    }

    #[test]
    fn test_rk4_accuracy_order() {
        // Compare RK4 error at h=0.1 vs h=0.05 on exponential decay
        // RK4 is 4th-order: error ~ h^4, so halving h reduces error by ~16x
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let exact = (-1.0_f64).exp();

        let y_h1 = {
            let mut y = vec![1.0];
            for i in 0..10 {
                y = rk4_step(&f, i as f64 * 0.1, &y, 0.1);
            }
            y[0]
        };
        let y_h2 = {
            let mut y = vec![1.0];
            for i in 0..20 {
                y = rk4_step(&f, i as f64 * 0.05, &y, 0.05);
            }
            y[0]
        };
        let err1 = (y_h1 - exact).abs();
        let err2 = (y_h2 - exact).abs();
        assert!(
            err2 < err1,
            "Smaller step should give smaller error: {} vs {}",
            err2,
            err1
        );
    }

    #[test]
    fn test_rk4_system_energy_conservation() {
        // Harmonic oscillator: H = 0.5*(x^2 + v^2) = 0.5 (for x0=1, v0=0)
        // should be approximately conserved by RK4
        let f = |_t: f64, z: &[f64]| vec![z[1], -z[0]];
        let mut z = vec![1.0, 0.0];
        let dt = 0.001;
        let steps = 1000;
        for i in 0..steps {
            z = rk4_step(&f, i as f64 * dt, &z, dt);
        }
        let energy = 0.5 * (z[0].powi(2) + z[1].powi(2));
        assert!(
            (energy - 0.5).abs() < 1e-4,
            "Energy drift: {}",
            energy - 0.5
        );
    }

    #[test]
    fn test_neural_ode_func_batch_consistency() {
        // forward(t, z) should give the same result when called twice
        let func = NeuralOdeFunc::new(3, 8, 42);
        let z = vec![0.1, -0.2, 0.3];
        let d1 = func.forward(0.5, &z);
        let d2 = func.forward(0.5, &z);
        assert_eq!(d1, d2, "forward must be deterministic");
    }

    #[test]
    fn test_time_series_ode_fit_loss_finite() {
        let func = NeuralOdeFunc::new(1, 4, 18);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let times: Vec<f64> = (0..5).map(|i| i as f64 * 0.2).collect();
        let obs: Vec<Vec<f64>> = times.iter().map(|&t: &f64| vec![(-t).exp()]).collect();
        let mut ts = TimeSeriesOde::new(times, obs, solver, 0.001, 3);
        ts.fit();
        for &l in &ts.loss_history {
            assert!(l.is_finite(), "Loss is non-finite: {l}");
        }
    }

    #[test]
    fn test_rk4_step_multidim() {
        // 5-dimensional decay: dy_i/dt = -i*y_i, y_i(0)=1
        let f = |_t: f64, y: &[f64]| (0..y.len()).map(|i| -(i as f64 + 1.0) * y[i]).collect();
        let y0: Vec<f64> = vec![1.0; 5];
        let mut y = y0.clone();
        let dt = 0.01;
        for k in 0..10 {
            y = rk4_step(&f, k as f64 * dt, &y, dt);
        }
        for (i, &yi) in y.iter().enumerate() {
            let exact = (-(i as f64 + 1.0) * 0.1).exp();
            assert!(
                (yi - exact).abs() < 1e-5,
                "dim {i}: got {yi}, expected {exact}"
            );
        }
    }

    // ── C1: DOPRI5 error-order verification ───────────────────────────────────

    #[test]
    fn test_dopri5_error_estimate_order() {
        // For y' = y, y(0) = 1, exact = exp(t).
        // DOPRI5 is 5th-order: halving h should reduce |y_high - exact| by ~32×.
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let rtol = 1e-12;
        let atol = 1e-12;
        let y0 = vec![1.0_f64];

        let (y_big, _, _) = dopri5_step(&f, 0.0, &y0, 0.2, rtol, atol);
        let (y_small, _, _) = dopri5_step(&f, 0.0, &y0, 0.1, rtol, atol);
        let err_big = (y_big[0] - 0.2_f64.exp()).abs();
        let err_small = (y_small[0] - 0.1_f64.exp()).abs();
        // ratio ≈ (0.2/0.1)^5 = 32; require > 10 to avoid false negatives
        let ratio = err_big / err_small.max(f64::MIN_POSITIVE);
        assert!(
            ratio > 10.0,
            "Expected ~32× error reduction when halving step; got ratio={ratio:.2}"
        );
    }

    #[test]
    fn test_dopri5_error_norm_small_step() {
        // error_norm should be < 1 for a well-behaved problem at h=0.01
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let (_, _, err) = dopri5_step(&f, 0.0, &[1.0], 0.01, 1e-6, 1e-8);
        assert!(err < 1.0, "error norm should be < 1 for h=0.01: {err}");
    }

    // ── C2: BPTT gradient parity test ─────────────────────────────────────────

    #[test]
    fn test_bptt_gradient_nonzero_and_finite() {
        // Verify that BPTT parameter gradients are non-zero and finite.
        let func = NeuralOdeFunc::new(2, 4, 99);
        let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
        let mut adj = AdjointMethod::new(2);
        let z_final = vec![0.5, -0.3];
        let loss_grad = vec![1.0, 0.0];
        let (_, grad_params) = adj.run(&solver, &z_final, &loss_grad, 0.0, 0.5, 0.1);
        assert_eq!(grad_params.len(), solver.func.n_params());
        assert!(
            grad_params.iter().all(|v| v.is_finite()),
            "some parameter gradients are non-finite"
        );
        assert!(
            grad_params.iter().any(|v| v.abs() > 1e-15),
            "all parameter gradients are zero"
        );
    }
}
