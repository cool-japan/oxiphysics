#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural network-accelerated Lattice Boltzmann methods.
//!
//! This module provides deep-learning augmented LBM components:
//! - Fully-connected neural networks for collision closures
//! - Physics-informed LBM with neural force corrections
//! - DeepONet-based solution operators
//! - Data-driven reduced-order models via DMD

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Activation functions
// ─────────────────────────────────────────────────────────────────────────────

/// Activation function used in neural network layers.
#[derive(Debug, Clone, PartialEq)]
pub enum ActivationFn {
    /// Rectified linear unit: max(0, x).
    Relu,
    /// Hyperbolic tangent.
    Tanh,
    /// Logistic sigmoid: 1 / (1 + e^{-x}).
    Sigmoid,
    /// Identity (no nonlinearity).
    Linear,
}

impl ActivationFn {
    /// Apply this activation function to a single value.
    pub fn apply(&self, x: f64) -> f64 {
        match self {
            ActivationFn::Relu => x.max(0.0),
            ActivationFn::Tanh => x.tanh(),
            ActivationFn::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            ActivationFn::Linear => x,
        }
    }

    /// Apply the element-wise derivative of this activation function.
    pub fn derivative(&self, x: f64) -> f64 {
        match self {
            ActivationFn::Relu => {
                if x > 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
            ActivationFn::Tanh => {
                let t = x.tanh();
                1.0 - t * t
            }
            ActivationFn::Sigmoid => {
                let s = 1.0 / (1.0 + (-x).exp());
                s * (1.0 - s)
            }
            ActivationFn::Linear => 1.0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Single dense layer
// ─────────────────────────────────────────────────────────────────────────────

/// A single fully-connected (dense) layer in a neural network.
///
/// `weights[i][j]` is the weight from input neuron `j` to output neuron `i`.
#[derive(Debug, Clone)]
pub struct NeuralLbmLayer {
    /// Weight matrix stored as `[output_neurons][input_neurons]`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector with one entry per output neuron.
    pub biases: Vec<f64>,
    /// Activation function applied after the linear transform.
    pub activation: ActivationFn,
}

impl NeuralLbmLayer {
    /// Create a new layer with given weight matrix, biases, and activation.
    ///
    /// # Panics
    /// Panics if `weights` and `biases` have different lengths.
    pub fn new(weights: Vec<Vec<f64>>, biases: Vec<f64>, activation: ActivationFn) -> Self {
        assert_eq!(
            weights.len(),
            biases.len(),
            "weights rows must match biases length"
        );
        Self {
            weights,
            biases,
            activation,
        }
    }

    /// Number of input features expected by this layer.
    pub fn input_size(&self) -> usize {
        self.weights.first().map_or(0, |row| row.len())
    }

    /// Number of output neurons in this layer.
    pub fn output_size(&self) -> usize {
        self.weights.len()
    }

    /// Forward pass: compute `activation(W * x + b)`.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        self.weights
            .iter()
            .zip(self.biases.iter())
            .map(|(row, &b)| {
                let pre = row
                    .iter()
                    .zip(input.iter())
                    .map(|(&w, &x)| w * x)
                    .sum::<f64>()
                    + b;
                self.activation.apply(pre)
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-layer perceptron
// ─────────────────────────────────────────────────────────────────────────────

/// A multi-layer perceptron (MLP) used to accelerate LBM computations.
///
/// The network is a stack of [`NeuralLbmLayer`] instances evaluated
/// sequentially during the forward pass.
#[derive(Debug, Clone)]
pub struct NeuralLbm {
    /// Ordered list of dense layers.
    pub layers: Vec<NeuralLbmLayer>,
    /// Number of inputs to the first layer.
    pub input_size: usize,
    /// Number of outputs from the final layer.
    pub output_size: usize,
}

impl NeuralLbm {
    /// Build a `NeuralLbm` from pre-constructed layers.
    ///
    /// `input_size` and `output_size` are inferred from the first and last
    /// layer if `layers` is non-empty; otherwise the supplied values are used.
    pub fn new(layers: Vec<NeuralLbmLayer>, input_size: usize, output_size: usize) -> Self {
        Self {
            layers,
            input_size,
            output_size,
        }
    }

    /// Create a small default network with random-like (deterministic) weights.
    ///
    /// Architecture: `input_size` → 16 (Tanh) → 16 (Tanh) → `output_size` (Linear).
    pub fn default_network(input_size: usize, output_size: usize) -> Self {
        let hidden = 16usize;
        // Use a simple deterministic pseudo-random initialisation so we stay
        // dependency-free (no external rand crate needed here).
        let layer1_w: Vec<Vec<f64>> = (0..hidden)
            .map(|i| {
                (0..input_size)
                    .map(|j| simple_weight(i, j, input_size) * 0.5)
                    .collect()
            })
            .collect();
        let layer1_b = vec![0.0f64; hidden];
        let layer2_w: Vec<Vec<f64>> = (0..hidden)
            .map(|i| {
                (0..hidden)
                    .map(|j| simple_weight(i, j, hidden) * 0.3)
                    .collect()
            })
            .collect();
        let layer2_b = vec![0.0f64; hidden];
        let layer3_w: Vec<Vec<f64>> = (0..output_size)
            .map(|i| {
                (0..hidden)
                    .map(|j| simple_weight(i, j, hidden) * 0.2)
                    .collect()
            })
            .collect();
        let layer3_b = vec![0.0f64; output_size];

        let layers = vec![
            NeuralLbmLayer::new(layer1_w, layer1_b, ActivationFn::Tanh),
            NeuralLbmLayer::new(layer2_w, layer2_b, ActivationFn::Tanh),
            NeuralLbmLayer::new(layer3_w, layer3_b, ActivationFn::Linear),
        ];
        Self::new(layers, input_size, output_size)
    }

    /// Run the forward pass through all layers.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let mut current = input.to_vec();
        for layer in &self.layers {
            current = layer.forward(&current);
        }
        current
    }

    /// Predict the equilibrium distribution `f_eq` from density `rho` and
    /// velocity components `u` using the neural network.
    ///
    /// The feature vector is `[rho, u[0\], u[1], ...]`.
    pub fn predict_equilibrium(&self, rho: f64, u: &[f64]) -> Vec<f64> {
        let mut features = vec![rho];
        features.extend_from_slice(u);
        // Pad or truncate to `input_size`.
        features.resize(self.input_size, 0.0);
        self.forward(&features)
    }

    /// Evaluate a turbulence closure model returning the eddy viscosity
    /// correction for a given `strain_rate` magnitude.
    pub fn closure_model(&self, strain_rate: f64) -> f64 {
        let features: Vec<f64> = std::iter::once(strain_rate)
            .chain(std::iter::repeat_n(0.0, self.input_size.saturating_sub(1)))
            .take(self.input_size)
            .collect();
        let out = self.forward(&features);
        *out.first().unwrap_or(&0.0)
    }

    /// Predict the turbulent (sub-grid) viscosity from a feature vector.
    ///
    /// `features` should contain flow-field descriptors such as local velocity
    /// gradients, vorticity, etc.
    pub fn turbulent_viscosity(&self, features: &[f64]) -> f64 {
        let mut padded = features.to_vec();
        padded.resize(self.input_size, 0.0);
        let out = self.forward(&padded);
        out.first().copied().unwrap_or(0.0).abs()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Physics-informed LBM
// ─────────────────────────────────────────────────────────────────────────────

/// Physics-informed LBM that couples a neural network correction term with a
/// standard LBM grid.
///
/// The neural network learns the residual between the coarse-grained LBM
/// solution and the reference (DNS / experimental) solution.
#[derive(Debug, Clone)]
pub struct PhysicsInformedLbm {
    /// Neural network that provides force corrections.
    pub neural_correction: NeuralLbm,
    /// Flattened distribution-function grid (`Q * Nx * Ny`).
    pub lbm_grid: Vec<f64>,
    /// Relaxation time τ.
    pub tau: f64,
    /// Number of velocity directions `Q`.
    pub q: usize,
    /// Number of grid nodes.
    pub n_nodes: usize,
}

impl PhysicsInformedLbm {
    /// Create a new `PhysicsInformedLbm` with a given grid size and
    /// default neural network.
    pub fn new(n_nodes: usize, q: usize, tau: f64) -> Self {
        let lbm_grid = vec![1.0 / q as f64; n_nodes * q];
        let neural_correction = NeuralLbm::default_network(q + 1, q);
        Self {
            neural_correction,
            lbm_grid,
            tau,
            q,
            n_nodes,
        }
    }

    /// Advance the simulation by one LBM time step with neural correction.
    pub fn step(&mut self) {
        let q = self.q;
        let tau = self.tau;
        for node in 0..self.n_nodes {
            let slice = &self.lbm_grid[node * q..(node + 1) * q];
            // Compute macroscopic density.
            let rho: f64 = slice.iter().sum();
            let correction = self.neural_correction.neural_force_correction_inner(slice);
            // BGK update with neural correction.
            let fi_new: Vec<f64> = slice
                .iter()
                .zip(correction.iter())
                .map(|(&fi, &corr)| {
                    let fi_eq = rho / q as f64;
                    fi - (fi - fi_eq) / tau + corr
                })
                .collect();
            self.lbm_grid[node * q..(node + 1) * q].copy_from_slice(&fi_new);
        }
    }

    /// Compute the L2 residual of the BGK equation across the grid.
    pub fn compute_residual(&self) -> f64 {
        let q = self.q;
        let tau = self.tau;
        let mut residual = 0.0f64;
        for node in 0..self.n_nodes {
            let slice = &self.lbm_grid[node * q..(node + 1) * q];
            let rho: f64 = slice.iter().sum();
            for &fi in slice {
                let fi_eq = rho / q as f64;
                let r = fi - (fi - fi_eq) / tau;
                residual += r * r;
            }
        }
        (residual / (self.n_nodes * q) as f64).sqrt()
    }

    /// Compute the neural force correction for a distribution slice.
    pub fn neural_force_correction(&self, fi: &[f64]) -> Vec<f64> {
        self.neural_correction.neural_force_correction_inner(fi)
    }
}

// Internal helper on NeuralLbm.
impl NeuralLbm {
    fn neural_force_correction_inner(&self, fi: &[f64]) -> Vec<f64> {
        let rho: f64 = fi.iter().sum();
        let mut features = vec![rho];
        features.extend_from_slice(fi);
        features.resize(self.input_size, 0.0);
        let raw = self.forward(&features);
        raw.into_iter().take(fi.len()).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DeepONet LBM
// ─────────────────────────────────────────────────────────────────────────────

/// Deep operator network (DeepONet) for LBM solution operators.
///
/// A DeepONet consists of a **branch net** encoding the input function (e.g.
/// initial/boundary conditions) and a **trunk net** encoding the query
/// location (space-time coordinates).  The output is their inner product.
#[derive(Debug, Clone)]
pub struct DeepONetLbm {
    /// Branch network — encodes the input function.
    pub branch_net: NeuralLbm,
    /// Trunk network — encodes the query coordinates.
    pub trunk_net: NeuralLbm,
}

impl DeepONetLbm {
    /// Create a `DeepONetLbm` with given branch and trunk networks.
    pub fn new(branch_net: NeuralLbm, trunk_net: NeuralLbm) -> Self {
        Self {
            branch_net,
            trunk_net,
        }
    }

    /// Create a default `DeepONetLbm` with small networks.
    pub fn default_network(branch_input: usize, trunk_input: usize, basis_size: usize) -> Self {
        let branch_net = NeuralLbm::default_network(branch_input, basis_size);
        let trunk_net = NeuralLbm::default_network(trunk_input, basis_size);
        Self::new(branch_net, trunk_net)
    }

    /// Evaluate the DeepONet: compute `branch(u) · trunk(y)`.
    ///
    /// - `branch_input`: the discrete input function values.
    /// - `trunk_input`: the query coordinates (e.g. `[x, y, t]`).
    pub fn evaluate(&self, branch_input: &[f64], trunk_input: &[f64]) -> f64 {
        let b = self.branch_net.forward(branch_input);
        let t = self.trunk_net.forward(trunk_input);
        b.iter().zip(t.iter()).map(|(&bi, &ti)| bi * ti).sum()
    }

    /// Predict the LBM solution at multiple query points given the input
    /// function encoded in `branch_input`.
    ///
    /// Returns one scalar per query point.
    pub fn predict_solution(&self, branch_input: &[f64], query_points: &[Vec<f64>]) -> Vec<f64> {
        let b = self.branch_net.forward(branch_input);
        query_points
            .iter()
            .map(|qp| {
                let t = self.trunk_net.forward(qp);
                b.iter().zip(t.iter()).map(|(&bi, &ti)| bi * ti).sum()
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Neural BGK collision
// ─────────────────────────────────────────────────────────────────────────────

/// Neural-network augmented BGK collision operator.
///
/// Adds a learned correction term `correction[i]` to the standard BGK update:
///
/// ```text
/// f_i* = f_i - (f_i - f_i^eq) / τ + correction_i
/// ```
pub fn neural_bgk_collision(fi: &[f64], fi_eq: &[f64], tau: f64, correction: &[f64]) -> Vec<f64> {
    assert_eq!(
        fi.len(),
        fi_eq.len(),
        "fi and fi_eq must have the same length"
    );
    let corr_len = correction.len();
    fi.iter()
        .zip(fi_eq.iter())
        .enumerate()
        .map(|(i, (&f, &feq))| {
            let corr = if i < corr_len { correction[i] } else { 0.0 };
            f - (f - feq) / tau + corr
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Data-driven LBM (DMD / ROM)
// ─────────────────────────────────────────────────────────────────────────────

/// Data-driven LBM using Dynamic Mode Decomposition (DMD) for reduced-order
/// modelling.
///
/// Snapshots of the distribution-function field are stored and used to extract
/// coherent spatial modes via DMD.  The resulting ROM can predict future states
/// without running the full LBM.
#[derive(Debug, Clone)]
pub struct DataDrivenLbm {
    /// Collection of solution snapshots; each snapshot is a flattened
    /// distribution-function field.
    pub snapshots: Vec<Vec<f64>>,
    /// DMD / POD spatial modes extracted from `snapshots`.
    pub modes: Vec<Vec<f64>>,
    /// DMD eigenvalues (real part only for simplicity).
    pub eigenvalues: Vec<f64>,
}

impl DataDrivenLbm {
    /// Create an empty `DataDrivenLbm`.
    pub fn new() -> Self {
        Self {
            snapshots: Vec::new(),
            modes: Vec::new(),
            eigenvalues: Vec::new(),
        }
    }

    /// Append a new snapshot to the data set.
    pub fn add_snapshot(&mut self, snapshot: Vec<f64>) {
        self.snapshots.push(snapshot);
    }

    /// Compute DMD modes from the stored snapshots using a simplified
    /// rank-1 approximation.
    ///
    /// After calling this method, `self.modes` and `self.eigenvalues` are
    /// populated.
    pub fn dmd_modes(&mut self) {
        if self.snapshots.len() < 2 {
            return;
        }
        let n = self.snapshots.len();
        let dim = self.snapshots[0].len();

        // Build X (all but last) and Y (all but first) matrices.
        // DMD approximates Y ≈ A X, so A ≈ Y * pinv(X).
        // We use the power iteration / rank-1 approx for simplicity.
        let mut modes = Vec::new();
        let mut eigenvalues = Vec::new();

        for k in 0..(n - 1) {
            let x = &self.snapshots[k];
            let y = &self.snapshots[k + 1];
            // Rayleigh quotient as scalar eigenvalue approximation.
            let xnorm2: f64 = x.iter().map(|v| v * v).sum();
            if xnorm2 < 1e-30 {
                eigenvalues.push(0.0);
                modes.push(vec![0.0; dim]);
                continue;
            }
            let lambda: f64 = x
                .iter()
                .zip(y.iter())
                .map(|(&xi, &yi)| xi * yi)
                .sum::<f64>()
                / xnorm2;
            eigenvalues.push(lambda);
            // Mode is the normalised snapshot.
            let xnorm = xnorm2.sqrt();
            let mode: Vec<f64> = x.iter().map(|v| v / xnorm).collect();
            modes.push(mode);
        }

        self.modes = modes;
        self.eigenvalues = eigenvalues;
    }

    /// Predict the state at future time step `steps_ahead` using the ROM.
    ///
    /// The prediction is a linear combination of modes evolved by their
    /// respective eigenvalues.
    pub fn rom_predict(&self, steps_ahead: usize) -> Vec<f64> {
        if self.modes.is_empty() {
            return Vec::new();
        }
        let dim = self.modes[0].len();
        let mut prediction = vec![0.0f64; dim];

        for (mode, &lambda) in self.modes.iter().zip(self.eigenvalues.iter()) {
            let amp = lambda.powi(steps_ahead as i32);
            for (p, &m) in prediction.iter_mut().zip(mode.iter()) {
                *p += amp * m;
            }
        }
        prediction
    }

    /// Return the number of snapshots currently stored.
    pub fn snapshot_count(&self) -> usize {
        self.snapshots.len()
    }

    /// Compute the reconstruction error for snapshot `k` using the stored modes.
    pub fn reconstruction_error(&self, k: usize) -> f64 {
        if k >= self.snapshots.len() || self.modes.is_empty() {
            return f64::NAN;
        }
        let snap = &self.snapshots[k];
        let mode = &self.modes[k.min(self.modes.len() - 1)];
        let dot: f64 = snap.iter().zip(mode.iter()).map(|(&s, &m)| s * m).sum();
        let reconstructed: Vec<f64> = mode.iter().map(|&m| dot * m).collect();
        snap.iter()
            .zip(reconstructed.iter())
            .map(|(&s, &r)| (s - r).powi(2))
            .sum::<f64>()
            .sqrt()
    }
}

impl Default for DataDrivenLbm {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Simple deterministic weight initialiser (no external dependencies).
///
/// Produces values in `[-1, 1]` based on row/column indices.
fn simple_weight(row: usize, col: usize, _n_cols: usize) -> f64 {
    ((row * 7 + col * 13) % 17) as f64 / 8.5 - 1.0
}

/// Compute the macroscopic density from a distribution function slice.
///
/// `rho = sum_i f_i`
pub fn compute_density(fi: &[f64]) -> f64 {
    fi.iter().sum()
}

/// Compute the macroscopic velocity component from distributions.
///
/// Assumes D2Q9 conventions with `e[i]` being the x-component of the
/// `i`-th velocity direction.
pub fn compute_velocity(fi: &[f64], e: &[f64]) -> f64 {
    let rho = compute_density(fi);
    if rho.abs() < 1e-15 {
        return 0.0;
    }
    fi.iter().zip(e.iter()).map(|(&f, &ei)| f * ei).sum::<f64>() / rho
}

/// Compute the Maxwell-Boltzmann equilibrium distribution for a single
/// velocity direction with weight `w`, lattice velocity `e·u`, and speed of
/// sound `cs`.
pub fn equilibrium_fi(rho: f64, eu: f64, u2: f64, w: f64, cs: f64) -> f64 {
    let cs2 = cs * cs;
    rho * w * (1.0 + eu / cs2 + eu * eu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- ActivationFn tests ---

    #[test]
    fn test_relu_positive() {
        assert_eq!(ActivationFn::Relu.apply(3.0), 3.0);
    }

    #[test]
    fn test_relu_negative() {
        assert_eq!(ActivationFn::Relu.apply(-2.0), 0.0);
    }

    #[test]
    fn test_tanh_zero() {
        assert!((ActivationFn::Tanh.apply(0.0)).abs() < 1e-15);
    }

    #[test]
    fn test_sigmoid_zero() {
        let s = ActivationFn::Sigmoid.apply(0.0);
        assert!((s - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_linear_passthrough() {
        assert_eq!(ActivationFn::Linear.apply(42.0), 42.0);
    }

    #[test]
    fn test_relu_derivative() {
        assert_eq!(ActivationFn::Relu.derivative(1.0), 1.0);
        assert_eq!(ActivationFn::Relu.derivative(-1.0), 0.0);
    }

    #[test]
    fn test_sigmoid_derivative_zero() {
        let d = ActivationFn::Sigmoid.derivative(0.0);
        assert!((d - 0.25).abs() < 1e-10);
    }

    // --- NeuralLbmLayer tests ---

    #[test]
    fn test_layer_forward_identity() {
        // 2×2 identity weight, zero bias, Linear activation.
        let layer = NeuralLbmLayer::new(
            vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            vec![0.0, 0.0],
            ActivationFn::Linear,
        );
        let out = layer.forward(&[3.0, 5.0]);
        assert!((out[0] - 3.0).abs() < 1e-12);
        assert!((out[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_layer_sizes() {
        let layer = NeuralLbmLayer::new(
            vec![vec![1.0, 2.0, 3.0]; 4],
            vec![0.0; 4],
            ActivationFn::Relu,
        );
        assert_eq!(layer.input_size(), 3);
        assert_eq!(layer.output_size(), 4);
    }

    #[test]
    fn test_layer_relu_clamps_negative() {
        let layer = NeuralLbmLayer::new(vec![vec![-1.0]], vec![0.0], ActivationFn::Relu);
        let out = layer.forward(&[1.0]);
        assert_eq!(out[0], 0.0);
    }

    // --- NeuralLbm tests ---

    #[test]
    fn test_neural_lbm_forward_output_size() {
        let net = NeuralLbm::default_network(4, 9);
        let out = net.forward(&[0.1, 0.2, 0.05, 0.0]);
        assert_eq!(out.len(), 9);
    }

    #[test]
    fn test_predict_equilibrium_output_size() {
        let net = NeuralLbm::default_network(3, 9);
        let out = net.predict_equilibrium(1.0, &[0.1, 0.05]);
        assert_eq!(out.len(), 9);
    }

    #[test]
    fn test_closure_model_returns_scalar() {
        let net = NeuralLbm::default_network(1, 1);
        let nu_t = net.closure_model(0.5);
        assert!(nu_t.is_finite());
    }

    #[test]
    fn test_turbulent_viscosity_non_negative() {
        let net = NeuralLbm::default_network(4, 1);
        let nu_t = net.turbulent_viscosity(&[0.1, 0.2, 0.3, 0.4]);
        assert!(nu_t >= 0.0);
    }

    #[test]
    fn test_neural_lbm_empty_layers() {
        let net = NeuralLbm::new(vec![], 3, 3);
        let out = net.forward(&[1.0, 2.0, 3.0]);
        assert_eq!(out, vec![1.0, 2.0, 3.0]);
    }

    // --- PhysicsInformedLbm tests ---

    #[test]
    fn test_pilbm_residual_initial() {
        let pil = PhysicsInformedLbm::new(10, 9, 1.0);
        let res = pil.compute_residual();
        assert!(res.is_finite());
    }

    #[test]
    fn test_pilbm_step_preserves_size() {
        let mut pil = PhysicsInformedLbm::new(4, 9, 0.6);
        let before = pil.lbm_grid.len();
        pil.step();
        assert_eq!(pil.lbm_grid.len(), before);
    }

    #[test]
    fn test_pilbm_neural_correction_size() {
        let pil = PhysicsInformedLbm::new(2, 9, 0.6);
        let fi = vec![1.0 / 9.0; 9];
        let corr = pil.neural_force_correction(&fi);
        assert_eq!(corr.len(), 9);
    }

    #[test]
    fn test_pilbm_multiple_steps() {
        let mut pil = PhysicsInformedLbm::new(4, 9, 1.0);
        for _ in 0..5 {
            pil.step();
        }
        assert!(pil.compute_residual().is_finite());
    }

    // --- DeepONetLbm tests ---

    #[test]
    fn test_deeponet_evaluate_finite() {
        let net = DeepONetLbm::default_network(9, 3, 8);
        let branch_in = vec![1.0 / 9.0; 9];
        let trunk_in = vec![0.5, 0.5, 0.0];
        let val = net.evaluate(&branch_in, &trunk_in);
        assert!(val.is_finite());
    }

    #[test]
    fn test_deeponet_predict_solution_count() {
        let net = DeepONetLbm::default_network(9, 2, 8);
        let branch_in = vec![0.1; 9];
        let queries = vec![vec![0.0, 0.0], vec![0.5, 0.5], vec![1.0, 1.0]];
        let preds = net.predict_solution(&branch_in, &queries);
        assert_eq!(preds.len(), 3);
    }

    // --- neural_bgk_collision tests ---

    #[test]
    fn test_bgk_no_correction_equilibrium() {
        let feq: Vec<f64> = vec![1.0 / 9.0; 9];
        let fi = feq.clone();
        let out = neural_bgk_collision(&fi, &feq, 1.0, &[]);
        for (a, b) in out.iter().zip(fi.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_bgk_correction_applied() {
        let fi = vec![0.2, 0.1];
        let feq = vec![0.15, 0.15];
        let corr = vec![0.01, -0.01];
        let out = neural_bgk_collision(&fi, &feq, 1.0, &corr);
        // out[0] = 0.2 - (0.2-0.15)/1.0 + 0.01 = 0.2 - 0.05 + 0.01 = 0.16
        assert!((out[0] - 0.16).abs() < 1e-12);
    }

    #[test]
    fn test_bgk_output_length() {
        let fi = vec![0.1; 5];
        let feq = vec![0.1; 5];
        let out = neural_bgk_collision(&fi, &feq, 1.5, &[0.0; 5]);
        assert_eq!(out.len(), 5);
    }

    // --- DataDrivenLbm tests ---

    #[test]
    fn test_dmd_modes_empty() {
        let mut ddl = DataDrivenLbm::new();
        ddl.dmd_modes(); // should not panic
        assert!(ddl.modes.is_empty());
    }

    #[test]
    fn test_dmd_modes_two_snapshots() {
        let mut ddl = DataDrivenLbm::new();
        ddl.add_snapshot(vec![1.0, 0.0, 0.0]);
        ddl.add_snapshot(vec![0.9, 0.1, 0.0]);
        ddl.dmd_modes();
        assert_eq!(ddl.modes.len(), 1);
        assert_eq!(ddl.eigenvalues.len(), 1);
    }

    #[test]
    fn test_rom_predict_size() {
        let mut ddl = DataDrivenLbm::new();
        for k in 0..5 {
            let snap: Vec<f64> = (0..9).map(|i| (i + k) as f64 * 0.1).collect();
            ddl.add_snapshot(snap);
        }
        ddl.dmd_modes();
        let pred = ddl.rom_predict(2);
        assert_eq!(pred.len(), 9);
    }

    #[test]
    fn test_snapshot_count() {
        let mut ddl = DataDrivenLbm::new();
        assert_eq!(ddl.snapshot_count(), 0);
        ddl.add_snapshot(vec![1.0, 2.0]);
        assert_eq!(ddl.snapshot_count(), 1);
    }

    #[test]
    fn test_reconstruction_error_out_of_bounds() {
        let ddl = DataDrivenLbm::new();
        assert!(ddl.reconstruction_error(0).is_nan());
    }

    // --- Helper utilities tests ---

    #[test]
    fn test_compute_density() {
        let fi = vec![0.1, 0.2, 0.3, 0.4];
        assert!((compute_density(&fi) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_velocity_zero_rho() {
        let fi = vec![0.0; 4];
        let e = vec![1.0, 0.0, -1.0, 0.0];
        assert_eq!(compute_velocity(&fi, &e), 0.0);
    }

    #[test]
    fn test_equilibrium_fi_uniform() {
        // At rest (u=0), equilibrium should be rho * w.
        let val = equilibrium_fi(1.0, 0.0, 0.0, 4.0 / 9.0, (1.0f64 / 3.0).sqrt());
        assert!((val - 4.0 / 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_simple_weight_range() {
        for r in 0..8 {
            for c in 0..8 {
                let w = simple_weight(r, c, 8);
                assert!(w >= -1.0 && w <= 1.0, "weight {w} out of range");
            }
        }
    }
}
