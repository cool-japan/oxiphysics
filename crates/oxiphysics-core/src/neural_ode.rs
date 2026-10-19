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
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Returns `base + scale · direction` as a new vector.
fn add_scaled(base: &[f64], scale: f64, direction: &[f64]) -> Vec<f64> {
    base.iter()
        .zip(direction.iter())
        .map(|(b, d)| b + scale * d)
        .collect()
}

/// Performs `target += scale · direction` in place.
fn add_scaled_in_place(target: &mut [f64], scale: f64, direction: &[f64]) {
    for (t, d) in target.iter_mut().zip(direction.iter()) {
        *t += scale * d;
    }
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
    /// Approximates the directional derivative `(∂f/∂z) · v`.  Note that the
    /// continuous adjoint needs the *transposed* product `(∂f/∂z)ᵀ · a`; see
    /// [`NeuralOdeFunc::vjp`] and [`NeuralOdeFunc::backprop`].
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

    /// Compute the vector-Jacobian product `(∂f/∂z)ᵀ · adj` via central finite
    /// differences with step `eps`.
    ///
    /// Component `j` of the result is `Σ_i adj_i · (f_i(z + ε e_j) − f_i(z − ε e_j)) / 2ε`.
    /// This is the derivative the adjoint equation `da/dt = −(∂f/∂z)ᵀ a` needs
    /// (the Jacobian is transposed, unlike [`NeuralOdeFunc::jvp`]).  For the
    /// exact value use [`NeuralOdeFunc::backprop`].
    pub fn vjp(&self, t: f64, z: &[f64], adj: &[f64], eps: f64) -> Vec<f64> {
        (0..z.len())
            .map(|j| {
                let mut z_plus = z.to_vec();
                let mut z_minus = z.to_vec();
                z_plus[j] += eps;
                z_minus[j] -= eps;
                let f_plus = self.forward(t, &z_plus);
                let f_minus = self.forward(t, &z_minus);
                adj.iter()
                    .zip(f_plus.iter().zip(f_minus.iter()))
                    .map(|(&a, (&fp, &fm))| a * (fp - fm) / (2.0 * eps))
                    .sum()
            })
            .collect()
    }

    /// Exact reverse-mode derivatives of the network output at `(t, z)` for the
    /// output cotangent `adj`.
    ///
    /// Back-propagates `adj` through
    /// `f = W_out·tanh(W_h·tanh(W_in·[z; t] + b_in) + b_h) + b_out` and returns
    /// `(vjp_z, grad_theta)`:
    /// * `vjp_z[j] = Σ_i adj_i · ∂f_i/∂z_j`, the product `(∂f/∂z)ᵀ · adj`,
    /// * `grad_theta[p] = Σ_i adj_i · ∂f_i/∂θ_p` in the layout of
    ///   [`NeuralOdeFunc::params_flat`].
    ///
    /// The result is exact up to rounding error and costs one forward and one
    /// backward pass, whereas [`NeuralOdeFunc::vjp`] and
    /// [`NeuralOdeFunc::param_grad_contrib`] need a forward pass per perturbed
    /// coordinate.
    pub fn backprop(&self, t: f64, z: &[f64], adj: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let out = self.input_size;
        let hidden = self.hidden_size;
        let mut aug = Vec::with_capacity(z.len() + 1);
        aug.extend_from_slice(z);
        aug.push(t);
        let inp = aug.len();
        let h1 = dense_tanh(&aug, &self.weights_in, &self.bias_in, hidden);
        let h2 = dense_tanh(&h1, &self.weights_hidden, &self.bias_hidden, hidden);

        // Output layer: f = W_out h2 + b_out.
        let cotangent: Vec<f64> = (0..out)
            .map(|i| adj.get(i).copied().unwrap_or(0.0))
            .collect();
        let grad_weights_out: Vec<f64> = (0..out * hidden)
            .map(|idx| cotangent[idx / hidden] * h2[idx % hidden])
            .collect();
        let grad_bias_out = cotangent.clone();

        // Second hidden layer: h2 = tanh(W_h h1 + b_h).
        let d_pre2: Vec<f64> = (0..hidden)
            .map(|k| {
                let d_h2: f64 = (0..out)
                    .map(|i| self.weights_out[i * hidden + k] * cotangent[i])
                    .sum();
                d_h2 * (1.0 - h2[k] * h2[k])
            })
            .collect();
        let grad_weights_hidden: Vec<f64> = (0..hidden * hidden)
            .map(|idx| d_pre2[idx / hidden] * h1[idx % hidden])
            .collect();
        let grad_bias_hidden = d_pre2.clone();

        // First hidden layer: h1 = tanh(W_in [z; t] + b_in).
        let d_pre1: Vec<f64> = (0..hidden)
            .map(|l| {
                let d_h1: f64 = (0..hidden)
                    .map(|k| self.weights_hidden[k * hidden + l] * d_pre2[k])
                    .sum();
                d_h1 * (1.0 - h1[l] * h1[l])
            })
            .collect();
        let grad_weights_in: Vec<f64> = (0..hidden * inp)
            .map(|idx| d_pre1[idx / inp] * aug[idx % inp])
            .collect();
        let grad_bias_in = d_pre1.clone();
        let vjp_z: Vec<f64> = (0..z.len())
            .map(|j| {
                (0..hidden)
                    .map(|l| self.weights_in[l * inp + j] * d_pre1[l])
                    .sum()
            })
            .collect();

        let mut grad_theta = Vec::with_capacity(self.n_params());
        grad_theta.extend_from_slice(&grad_weights_in);
        grad_theta.extend_from_slice(&grad_bias_in);
        grad_theta.extend_from_slice(&grad_weights_hidden);
        grad_theta.extend_from_slice(&grad_bias_hidden);
        grad_theta.extend_from_slice(&grad_weights_out);
        grad_theta.extend_from_slice(&grad_bias_out);
        (vjp_z, grad_theta)
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

    /// Exact reverse-mode derivative of one fixed RK4 step.
    ///
    /// Let `y₁ = Φ_θ(y₀)` be the classical RK4 step of size `h` from `(t, y)`.
    /// Given the cotangent `y1_bar = ∂L/∂y₁`, this back-propagates through the
    /// four stages (using [`NeuralOdeFunc::backprop`]) and returns
    /// `(∂L/∂y₀, ∂L/∂θ)`, where `∂L/∂θ` follows the layout of
    /// [`NeuralOdeFunc::params_flat`].  Chaining it over the steps of a solve
    /// gives the exact gradient of the *discretised* solution map.
    pub fn rk4_step_vjp(&self, t: f64, y: &[f64], h: f64, y1_bar: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let func = &self.func;
        // Forward stages (needed again for the backward sweep).
        let k1 = func.forward(t, y);
        let y2 = add_scaled(y, 0.5 * h, &k1);
        let k2 = func.forward(t + 0.5 * h, &y2);
        let y3 = add_scaled(y, 0.5 * h, &k2);
        let k3 = func.forward(t + 0.5 * h, &y3);
        let y4 = add_scaled(y, h, &k3);

        // y₁ = y₀ + h/6 (k₁ + 2k₂ + 2k₃ + k₄)
        let mut y_bar = y1_bar.to_vec();
        let mut theta_bar = vec![0.0_f64; func.n_params()];
        let mut k1_bar: Vec<f64> = y1_bar.iter().map(|v| v * h / 6.0).collect();
        let mut k2_bar: Vec<f64> = y1_bar.iter().map(|v| v * h / 3.0).collect();
        let mut k3_bar: Vec<f64> = y1_bar.iter().map(|v| v * h / 3.0).collect();
        let k4_bar: Vec<f64> = y1_bar.iter().map(|v| v * h / 6.0).collect();

        // k₄ = f(t + h, y₀ + h k₃)
        let (j4, g4) = func.backprop(t + h, &y4, &k4_bar);
        add_scaled_in_place(&mut theta_bar, 1.0, &g4);
        add_scaled_in_place(&mut y_bar, 1.0, &j4);
        add_scaled_in_place(&mut k3_bar, h, &j4);
        // k₃ = f(t + h/2, y₀ + h/2 k₂)
        let (j3, g3) = func.backprop(t + 0.5 * h, &y3, &k3_bar);
        add_scaled_in_place(&mut theta_bar, 1.0, &g3);
        add_scaled_in_place(&mut y_bar, 1.0, &j3);
        add_scaled_in_place(&mut k2_bar, 0.5 * h, &j3);
        // k₂ = f(t + h/2, y₀ + h/2 k₁)
        let (j2, g2) = func.backprop(t + 0.5 * h, &y2, &k2_bar);
        add_scaled_in_place(&mut theta_bar, 1.0, &g2);
        add_scaled_in_place(&mut y_bar, 1.0, &j2);
        add_scaled_in_place(&mut k1_bar, 0.5 * h, &j2);
        // k₁ = f(t, y₀)
        let (j1, g1) = func.backprop(t, y, &k1_bar);
        add_scaled_in_place(&mut theta_bar, 1.0, &g1);
        add_scaled_in_place(&mut y_bar, 1.0, &j1);
        (y_bar, theta_bar)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AdjointMethod
// ─────────────────────────────────────────────────────────────────────────────

/// Reverse-mode gradient of a Neural ODE loss via the continuous adjoint method.
///
/// The adjoint state `a(t) = ∂L/∂z(t)` obeys `da/dt = −(∂f/∂z)ᵀ a` and is
/// integrated backward in time together with the state itself.  This gives
/// parameter gradients without storing the forward trajectory.
#[derive(Debug, Clone)]
pub struct AdjointMethod {
    /// Augmented state `[z; a; dL/dθ]` of the backward integration.  After
    /// [`AdjointMethod::run`] it holds the state at `t0`.
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

    /// One-step adjoint VJP: propagate `loss_grad` backward through a single RK4
    /// step of `solver`'s dynamics.
    ///
    /// Given the loss gradient `g = ∂L/∂z₁` at the *output* state `z₁ = Φ(z₀)`
    /// (where `Φ` is one fixed-step RK4 integration of `dz/dt = f(t,z)` over
    /// `[t0, t0 + dt]`), this returns the gradient with respect to the *input*
    /// state:
    ///
    ///   `∂L/∂z₀ = (∂Φ/∂z₀)ᵀ · g`
    ///
    /// The transpose-Jacobian–vector product `(∂Φ/∂z₀)ᵀ · g` is formed column by
    /// column via central finite differences of the actual RK4 step map: for each
    /// input coordinate `j`, the `j`-th component of the result is
    /// `Σ_i g_i · [Φ(z₀ + ε e_j) − Φ(z₀ − ε e_j)]_i / (2ε)`.
    ///
    /// This is the honest reverse-mode gradient of one integration step — it
    /// invokes the real dynamics, so its value depends on `z0` and `f`, not just
    /// on the sign of `loss_grad`.
    pub fn backward(
        &self,
        solver: &NeuralOdeSolver,
        z0: &[f64],
        loss_grad: &[f64],
        t0: f64,
        dt: f64,
    ) -> Vec<f64> {
        let n = self.state_dim.min(z0.len()).min(loss_grad.len());
        if n == 0 {
            return vec![0.0; self.state_dim];
        }
        // Relative finite-difference step, robust to the scale of z0.
        let eps_base = 1e-6;
        let forward = |t: f64, y: &[f64]| solver.func.forward(t, y);

        let mut grad = vec![0.0_f64; z0.len()];
        for j in 0..n {
            let eps = eps_base * (1.0 + z0[j].abs());
            let mut z_plus = z0.to_vec();
            let mut z_minus = z0.to_vec();
            z_plus[j] += eps;
            z_minus[j] -= eps;
            // Φ(z₀ ± ε e_j): a single RK4 step of the real dynamics.
            let phi_plus = rk4_step(&forward, t0, &z_plus, dt);
            let phi_minus = rk4_step(&forward, t0, &z_minus, dt);
            // (column j of ∂Φ/∂z₀) ⋅ loss_grad  →  component j of (∂Φ/∂z₀)ᵀ g.
            let mut acc = 0.0;
            for i in 0..n {
                acc += loss_grad[i] * (phi_plus[i] - phi_minus[i]) / (2.0 * eps);
            }
            grad[j] = acc;
        }
        grad
    }

    /// Propagate the final adjoint state `loss_grad = ∂L/∂z(t1)` backward through
    /// `solver` from `t1` to `t0` (continuous adjoint method) and return
    /// `(grad_z0, grad_params)`.
    ///
    /// `z_final = z(t1)` is the end state of the forward solve.  In the reversed
    /// time `s = t1 − t` the augmented state `[z, a, g]` obeys
    ///
    /// ```text
    /// dz/ds = −f(t1 − s, z)             (the trajectory is re-integrated backward)
    /// da/ds = +(∂f/∂z)ᵀ · a             (adjoint, a(0) = ∂L/∂z(t1))
    /// dg/ds = +aᵀ · ∂f/∂θ               (parameter gradient, g(0) = 0)
    /// ```
    ///
    /// and is advanced by classical RK4 with all four stages evaluated at the
    /// stage states and stage times of the *augmented* system.  The step grid is
    /// the one [`NeuralOdeSolver::solve_rk4`] uses for the same `(t0, t1, dt)`,
    /// traversed in reverse.  The transposed Jacobian and the parameter
    /// derivatives are exact ([`NeuralOdeFunc::backprop`]).
    ///
    /// `grad_z0` is `∂L/∂z(t0)`; `grad_params` is `∂L/∂θ` in the layout of
    /// [`NeuralOdeFunc::params_flat`].  The result agrees with the gradient of the
    /// discretised forward solve up to the RK4 truncation error, which vanishes
    /// as `dt → 0`.  Afterwards [`AdjointMethod::augmented_state`] holds
    /// `[z(t0); grad_z0; grad_params]`, with `z(t0)` the backward reconstruction
    /// of the initial state.
    pub fn run(
        &mut self,
        solver: &NeuralOdeSolver,
        z_final: &[f64],
        loss_grad: &[f64],
        t0: f64,
        t1: f64,
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let func = &solver.func;
        let n_params = func.n_params();
        let n = self.state_dim.min(z_final.len()).min(loss_grad.len());
        if n == 0 {
            return (vec![0.0; self.state_dim], vec![0.0; n_params]);
        }
        let h_step = dt.abs().max(1e-10);

        // Step grid of the forward solve (same loop as `NeuralOdeSolver::solve_rk4`).
        let mut times = vec![t0];
        let mut t_cur = t0;
        while t_cur < t1 - 1e-12 {
            let next = t_cur + h_step.min(t1 - t_cur);
            if next <= t_cur {
                break;
            }
            t_cur = next;
            times.push(t_cur);
        }

        // Right-hand side of the augmented system in reversed time.
        let rhs = |t: f64, z: &[f64], a: &[f64]| -> (Vec<f64>, Vec<f64>, Vec<f64>) {
            let neg_f: Vec<f64> = func.forward(t, z).into_iter().map(|v| -v).collect();
            let (vjp_z, grad_theta) = func.backprop(t, z, a);
            (neg_f, vjp_z, grad_theta)
        };

        let mut z = z_final[..n].to_vec();
        let mut adj = loss_grad[..n].to_vec();
        let mut grad_params = vec![0.0_f64; n_params];
        for window in times.windows(2).rev() {
            let (t_lo, t_hi) = (window[0], window[1]);
            let h = t_hi - t_lo;
            let t_mid = 0.5 * (t_lo + t_hi);
            let (kz1, ka1, kg1) = rhs(t_hi, &z, &adj);
            let (kz2, ka2, kg2) = rhs(
                t_mid,
                &add_scaled(&z, 0.5 * h, &kz1),
                &add_scaled(&adj, 0.5 * h, &ka1),
            );
            let (kz3, ka3, kg3) = rhs(
                t_mid,
                &add_scaled(&z, 0.5 * h, &kz2),
                &add_scaled(&adj, 0.5 * h, &ka2),
            );
            let (kz4, ka4, kg4) = rhs(t_lo, &add_scaled(&z, h, &kz3), &add_scaled(&adj, h, &ka3));
            for i in 0..n {
                z[i] += h / 6.0 * (kz1[i] + 2.0 * kz2[i] + 2.0 * kz3[i] + kz4[i]);
                adj[i] += h / 6.0 * (ka1[i] + 2.0 * ka2[i] + 2.0 * ka3[i] + ka4[i]);
            }
            for (j, g) in grad_params.iter_mut().enumerate() {
                *g += h / 6.0 * (kg1[j] + 2.0 * kg2[j] + 2.0 * kg3[j] + kg4[j]);
            }
        }

        let mut augmented = z;
        augmented.extend_from_slice(&adj);
        augmented.extend_from_slice(&grad_params);
        self.augmented_state = augmented;
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

/// Start of one RK4 step recorded by [`TimeSeriesOde::loss_gradient`].
struct Rk4StepRecord {
    /// Time at the start of the step.
    t: f64,
    /// Step size.
    h: f64,
    /// State at the start of the step.
    state: Vec<f64>,
}

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
    /// Every iteration records the current MSE in [`TimeSeriesOde::loss_history`]
    /// and then moves all parameters along the negative loss gradient,
    /// `θ ← θ − learning_rate · ∇L`, where `∇L` comes from
    /// [`TimeSeriesOde::loss_gradient`].  The integration step is a tenth of the
    /// mean spacing of the observation times.  If the gradient becomes
    /// non-finite (the learning rate is too large and the model diverged) the
    /// run stops early and the parameters keep their last finite values, so
    /// `loss_history` may then hold fewer than `n_iter` entries.
    pub fn fit(&mut self) {
        let span_dt = if self.times.len() > 1 {
            (self.times[self.times.len() - 1] - self.times[0]) / (self.times.len() as f64 * 10.0)
        } else {
            0.01
        };
        let dt = if span_dt.is_finite() && span_dt > 0.0 {
            span_dt
        } else {
            0.01
        };

        for _iter in 0..self.n_iter {
            let loss = self.compute_loss(dt);
            self.loss_history.push(loss);

            let gradient = self.loss_gradient(dt);
            if gradient.iter().any(|g| !g.is_finite()) {
                break;
            }
            let mut params = self.solver.func.params_flat();
            for (p, g) in params.iter_mut().zip(gradient.iter()) {
                *p -= self.learning_rate * g;
            }
            self.solver.func.set_params_flat(&params);
        }
    }

    /// Gradient of [`TimeSeriesOde::compute_loss`] with respect to all trainable
    /// parameters, in the layout of [`NeuralOdeFunc::params_flat`].
    ///
    /// This is the exact derivative of the discretised loss: the RK4 solve of
    /// `compute_loss(dt)` is replayed once forward (remembering the start state
    /// of every step) and then differentiated in reverse through every step with
    /// [`NeuralOdeSolver::rk4_step_vjp`], injecting the residual cotangent
    /// `2 (z_k − o_k) / count` at each observation time.  The initial state is
    /// the first observation, which does not depend on the parameters.  Returns
    /// zeros when there is nothing to fit.
    pub fn loss_gradient(&self, dt: f64) -> Vec<f64> {
        let mut grad = vec![0.0_f64; self.solver.func.n_params()];
        if self.times.is_empty() || self.observations.is_empty() {
            return grad;
        }
        let forward = |t: f64, y: &[f64]| self.solver.func.forward(t, y);

        // Forward pass identical to `solve_rk4_trajectory`, with a record per step.
        let mut state = self.observations[0].clone();
        let mut t = self.times[0];
        let mut steps: Vec<Rk4StepRecord> = Vec::new();
        // (number of steps taken so far, state) for every observation time.
        let mut snapshots: Vec<(usize, Vec<f64>)> = vec![(0, state.clone())];
        for &t_next in self.times.iter().skip(1) {
            while t < t_next - 1e-12 {
                let h = dt.min(t_next - t);
                steps.push(Rk4StepRecord {
                    t,
                    h,
                    state: state.clone(),
                });
                state = rk4_step(&forward, t, &state, h);
                t += h;
            }
            snapshots.push((steps.len(), state.clone()));
        }

        // Number of squared residuals in the MSE (mirrors `compute_loss`).
        let paired = snapshots.len().min(self.observations.len());
        let count: usize = snapshots
            .iter()
            .zip(self.observations.iter())
            .map(|((_, pred), obs)| pred.len().min(obs.len()))
            .sum();
        if count == 0 {
            return grad;
        }
        let scale = 2.0 / count as f64;

        // Reverse sweep over the observation intervals.
        let mut cotangent = vec![0.0_f64; snapshots[0].1.len()];
        for k in (0..paired).rev() {
            let (step_end, pred) = &snapshots[k];
            for ((c, p), o) in cotangent
                .iter_mut()
                .zip(pred.iter())
                .zip(self.observations[k].iter())
            {
                *c += scale * (p - o);
            }
            if k == 0 {
                break;
            }
            let step_start = snapshots[k - 1].0;
            for record in steps[step_start..*step_end].iter().rev() {
                let (state_bar, theta_bar) =
                    self.solver
                        .rk4_step_vjp(record.t, &record.state, record.h, &cotangent);
                cotangent = state_bar;
                add_scaled_in_place(&mut grad, 1.0, &theta_bar);
            }
        }
        grad
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
mod tests;
