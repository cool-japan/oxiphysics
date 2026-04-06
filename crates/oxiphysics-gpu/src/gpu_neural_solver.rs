// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU neural network solver for physics (CPU mock backend).
//!
//! Provides a multi-layer perceptron (MLP) framework and physics-informed
//! neural network (PINN) utilities:
//! - [`NeuralLayer`]: single dense layer with several activation modes.
//! - [`GpuNeuralSolver`]: stacked MLP forward pass.
//! - [`PhysicsNeuralNet`]: PINN wrapper with PDE + boundary-condition loss.
//! - Activation functions: [`ns_relu`], [`ns_sigmoid`], [`ns_softmax`].
//! - Loss functions: [`ns_mse_loss`], [`ns_mae_loss`].
//! - PINN residuals: [`pinn_residual`], [`pinn_boundary_loss`].

#![allow(dead_code)]

// ── Activation functions ──────────────────────────────────────────────────────

/// Rectified linear unit activation: `max(0, x)`.
pub fn ns_relu(x: f64) -> f64 {
    x.max(0.0)
}

/// Logistic sigmoid activation: `1 / (1 + e^{-x})`.
pub fn ns_sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Softmax of a slice: normalised exponentials `exp(x_i) / Σ exp(x_j)`.
///
/// Uses the max-subtraction trick for numerical stability.
/// Returns an empty `Vec` when `x` is empty.
pub fn ns_softmax(x: &[f64]) -> Vec<f64> {
    if x.is_empty() {
        return Vec::new();
    }
    let max_val = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = x.iter().map(|&v| (v - max_val).exp()).collect();
    let sum: f64 = exps.iter().sum();
    exps.iter().map(|&e| e / sum).collect()
}

// ── Loss functions ────────────────────────────────────────────────────────────

/// Mean-squared error loss: `mean((predicted_i − target_i)²)`.
///
/// Returns `0.0` when `predicted` is empty.
pub fn ns_mse_loss(predicted: &[f64], target: &[f64]) -> f64 {
    if predicted.is_empty() {
        return 0.0;
    }
    let n = predicted.len().min(target.len());
    predicted[..n]
        .iter()
        .zip(target[..n].iter())
        .map(|(p, t)| (p - t).powi(2))
        .sum::<f64>()
        / n as f64
}

/// Mean absolute error loss: `mean(|predicted_i − target_i|)`.
///
/// Returns `0.0` when `predicted` is empty.
pub fn ns_mae_loss(predicted: &[f64], target: &[f64]) -> f64 {
    if predicted.is_empty() {
        return 0.0;
    }
    let n = predicted.len().min(target.len());
    predicted[..n]
        .iter()
        .zip(target[..n].iter())
        .map(|(p, t)| (p - t).abs())
        .sum::<f64>()
        / n as f64
}

// ── PINN residuals ────────────────────────────────────────────────────────────

/// Physics-informed residual for the 1-D Poisson equation `−u_xx = source`.
///
/// Returns `−u_xx − source`; this should be driven to zero during training.
pub fn pinn_residual(u: f64, u_xx: f64, source: f64) -> f64 {
    let _ = u; // u itself is not needed for the Poisson residual
    -u_xx - source
}

/// Boundary-condition loss: MSE between the predicted boundary values and the
/// prescribed Dirichlet data.
pub fn pinn_boundary_loss(u_boundary: &[f64], u_target: &[f64]) -> f64 {
    ns_mse_loss(u_boundary, u_target)
}

// ── NeuralLayer ───────────────────────────────────────────────────────────────

/// A single fully-connected (dense) neural network layer.
///
/// Weights are stored in row-major order: `weights[i * n_in + j]` is the
/// weight from input `j` to output neuron `i`.
#[derive(Debug, Clone)]
pub struct NeuralLayer {
    /// Flattened weight matrix of shape `[n_out × n_in]`.
    pub weights: Vec<f64>,
    /// Bias vector of length `n_out`.
    pub biases: Vec<f64>,
    /// Number of input features.
    pub n_in: usize,
    /// Number of output neurons.
    pub n_out: usize,
}

impl NeuralLayer {
    /// Create a new layer with all weights and biases initialised to `0.1`.
    pub fn new(n_in: usize, n_out: usize) -> Self {
        Self {
            weights: vec![0.1; n_out * n_in],
            biases: vec![0.0; n_out],
            n_in,
            n_out,
        }
    }

    /// Linear forward pass (no activation): `W·x + b`.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let n = self.n_in.min(input.len());
        (0..self.n_out)
            .map(|i| {
                let base = i * self.n_in;
                let dot: f64 = (0..n).map(|j| self.weights[base + j] * input[j]).sum();
                dot + self.biases[i]
            })
            .collect()
    }

    /// Forward pass with ReLU activation applied element-wise.
    pub fn relu_forward(&self, input: &[f64]) -> Vec<f64> {
        self.forward(input).into_iter().map(ns_relu).collect()
    }

    /// Forward pass with tanh activation applied element-wise.
    pub fn tanh_forward(&self, input: &[f64]) -> Vec<f64> {
        self.forward(input).into_iter().map(|v| v.tanh()).collect()
    }

    /// Number of output neurons.
    pub fn output_size(&self) -> usize {
        self.n_out
    }

    /// Number of input features.
    pub fn input_size(&self) -> usize {
        self.n_in
    }
}

// ── GpuNeuralSolver ───────────────────────────────────────────────────────────

/// Multi-layer perceptron running on the CPU mock GPU backend.
///
/// Layers are stacked so that the output of layer `i` feeds into layer `i+1`.
/// All hidden activations use ReLU; the final layer is linear.
#[derive(Debug, Clone)]
pub struct GpuNeuralSolver {
    /// Ordered list of dense layers.
    pub layers: Vec<NeuralLayer>,
    /// Learning rate (stored for future back-prop use).
    pub learning_rate: f64,
}

impl GpuNeuralSolver {
    /// Build a network from an ordered list of layer widths.
    ///
    /// `layer_sizes` must contain at least two elements: the first is the
    /// input dimension and the last is the output dimension.
    pub fn new(layer_sizes: &[usize], lr: f64) -> Self {
        assert!(
            layer_sizes.len() >= 2,
            "Need at least input and output sizes"
        );
        let layers = layer_sizes
            .windows(2)
            .map(|w| NeuralLayer::new(w[0], w[1]))
            .collect();
        Self {
            layers,
            learning_rate: lr,
        }
    }

    /// Run a full forward pass through all layers (ReLU hidden, linear output).
    pub fn forward_pass(&self, input: &[f64]) -> Vec<f64> {
        let mut x: Vec<f64> = input.to_vec();
        let last = self.layers.len().saturating_sub(1);
        for (i, layer) in self.layers.iter().enumerate() {
            x = if i < last {
                layer.relu_forward(&x)
            } else {
                layer.forward(&x)
            };
        }
        x
    }

    /// Number of layers (= number of weight matrices).
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }

    /// Alias for [`forward_pass`](GpuNeuralSolver::forward_pass).
    pub fn predict(&self, input: &[f64]) -> Vec<f64> {
        self.forward_pass(input)
    }
}

// ── PhysicsNeuralNet ──────────────────────────────────────────────────────────

/// Physics-informed neural network (PINN).
///
/// Wraps a [`GpuNeuralSolver`] and provides a composite loss that combines
/// a PDE residual term with a boundary-condition MSE term.
#[derive(Debug, Clone)]
pub struct PhysicsNeuralNet {
    /// Underlying neural solver.
    pub solver: GpuNeuralSolver,
    /// Weight applied to the PDE residual in the total loss.
    pub pde_weight: f64,
    /// Weight applied to the boundary-condition loss in the total loss.
    pub bc_weight: f64,
}

impl PhysicsNeuralNet {
    /// Create a new PINN from layer sizes, PDE weight, and BC weight.
    pub fn new(layer_sizes: &[usize], pde_weight: f64, bc_weight: f64) -> Self {
        Self {
            solver: GpuNeuralSolver::new(layer_sizes, 1e-3),
            pde_weight,
            bc_weight,
        }
    }

    /// Compute the weighted total loss:
    /// `pde_weight * |pde_residual| + bc_weight * bc_loss`.
    pub fn total_loss(&self, pde_residual: f64, bc_loss: f64) -> f64 {
        self.pde_weight * pde_residual.abs() + self.bc_weight * bc_loss
    }

    /// Predict the solution at a 1-D coordinate `x`.
    ///
    /// Wraps the scalar input into a slice and extracts the first output.
    pub fn predict(&self, x: f64) -> f64 {
        let out = self.solver.predict(&[x]);
        out.first().copied().unwrap_or(0.0)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ns_relu ──────────────────────────────────────────────────────────

    #[test]
    fn relu_negative_is_zero() {
        assert!((ns_relu(-1.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn relu_positive_identity() {
        assert!((ns_relu(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn relu_zero_boundary() {
        assert!((ns_relu(0.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn relu_large_positive() {
        assert!((ns_relu(1000.0) - 1000.0).abs() < 1e-8);
    }

    // ── ns_sigmoid ───────────────────────────────────────────────────────

    #[test]
    fn sigmoid_at_zero_is_half() {
        assert!((ns_sigmoid(0.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn sigmoid_large_positive_near_one() {
        assert!(ns_sigmoid(100.0) > 0.999);
    }

    #[test]
    fn sigmoid_large_negative_near_zero() {
        assert!(ns_sigmoid(-100.0) < 0.001);
    }

    #[test]
    fn sigmoid_symmetry() {
        let s = ns_sigmoid(2.0);
        assert!((ns_sigmoid(-2.0) - (1.0 - s)).abs() < 1e-12);
    }

    // ── ns_softmax ────────────────────────────────────────────────────────

    #[test]
    fn softmax_sums_to_one() {
        let x = [1.0, 2.0, 3.0];
        let s = ns_softmax(&x);
        let total: f64 = s.iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn softmax_empty_input() {
        let s = ns_softmax(&[]);
        assert!(s.is_empty());
    }

    #[test]
    fn softmax_single_element() {
        let s = ns_softmax(&[42.0]);
        assert!((s[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn softmax_uniform_input() {
        let x = [1.0f64; 4];
        let s = ns_softmax(&x);
        for &v in &s {
            assert!((v - 0.25).abs() < 1e-12);
        }
    }

    #[test]
    fn softmax_all_non_negative() {
        let x = [-3.0, 0.0, 1.0, 5.0];
        let s = ns_softmax(&x);
        for &v in &s {
            assert!(v >= 0.0);
        }
    }

    // ── ns_mse_loss ──────────────────────────────────────────────────────

    #[test]
    fn mse_zero_for_identical() {
        let v = [1.0, 2.0, 3.0];
        assert!((ns_mse_loss(&v, &v) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn mse_known_value() {
        let pred = [3.0];
        let target = [1.0];
        // (3-1)^2 / 1 = 4
        assert!((ns_mse_loss(&pred, &target) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn mse_empty_returns_zero() {
        assert!((ns_mse_loss(&[], &[]) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn mse_positive_values() {
        let pred = [1.0, 2.0];
        let target = [0.0, 0.0];
        let loss = ns_mse_loss(&pred, &target);
        assert!(loss > 0.0);
    }

    // ── ns_mae_loss ──────────────────────────────────────────────────────

    #[test]
    fn mae_zero_for_identical() {
        let v = [1.0, 2.0, 3.0];
        assert!((ns_mae_loss(&v, &v) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn mae_known_value() {
        let pred = [3.0, 1.0];
        let target = [1.0, 1.0];
        // |3-1| + |1-1| = 2; / 2 = 1
        assert!((ns_mae_loss(&pred, &target) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mae_empty_returns_zero() {
        assert!((ns_mae_loss(&[], &[]) - 0.0).abs() < 1e-12);
    }

    // ── NeuralLayer ───────────────────────────────────────────────────────

    #[test]
    fn neural_layer_output_size() {
        let layer = NeuralLayer::new(4, 3);
        assert_eq!(layer.output_size(), 3);
    }

    #[test]
    fn neural_layer_input_size() {
        let layer = NeuralLayer::new(4, 3);
        assert_eq!(layer.input_size(), 4);
    }

    #[test]
    fn neural_layer_forward_output_length() {
        let layer = NeuralLayer::new(4, 3);
        let out = layer.forward(&[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn neural_layer_relu_forward_non_negative() {
        let layer = NeuralLayer::new(2, 4);
        let out = layer.relu_forward(&[-10.0, -10.0]);
        for &v in &out {
            assert!(v >= 0.0);
        }
    }

    #[test]
    fn neural_layer_tanh_bounded() {
        let layer = NeuralLayer::new(3, 3);
        let out = layer.tanh_forward(&[1.0, 2.0, 3.0]);
        for &v in &out {
            assert!(v > -1.0 && v < 1.0);
        }
    }

    #[test]
    fn neural_layer_zero_input() {
        // With all weights=0.1 and biases=0, output should be 0
        let mut layer = NeuralLayer::new(3, 2);
        layer.weights = vec![0.0; 6];
        let out = layer.forward(&[0.0, 0.0, 0.0]);
        for &v in &out {
            assert!(v.abs() < 1e-12);
        }
    }

    // ── GpuNeuralSolver ───────────────────────────────────────────────────

    #[test]
    fn solver_layer_count() {
        let s = GpuNeuralSolver::new(&[4, 8, 8, 2], 1e-3);
        assert_eq!(s.layer_count(), 3);
    }

    #[test]
    fn solver_forward_output_shape() {
        let s = GpuNeuralSolver::new(&[3, 5, 2], 1e-3);
        let out = s.forward_pass(&[1.0, 0.0, -1.0]);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn solver_predict_same_as_forward() {
        let s = GpuNeuralSolver::new(&[2, 4, 1], 1e-3);
        let input = [0.5, -0.5];
        let a = s.forward_pass(&input);
        let b = s.predict(&input);
        assert_eq!(a, b);
    }

    #[test]
    fn solver_single_layer() {
        let s = GpuNeuralSolver::new(&[2, 1], 1e-3);
        let out = s.forward_pass(&[1.0, 1.0]);
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn solver_deep_network_no_panic() {
        let s = GpuNeuralSolver::new(&[10, 20, 20, 20, 5], 1e-4);
        let input = vec![0.1; 10];
        let out = s.forward_pass(&input);
        assert_eq!(out.len(), 5);
    }

    // ── pinn_residual ────────────────────────────────────────────────────

    #[test]
    fn pinn_residual_formula() {
        // residual = -u_xx - source
        let r = pinn_residual(0.0, 2.0, 1.0);
        assert!((r - (-3.0)).abs() < 1e-12);
    }

    #[test]
    fn pinn_residual_zero_when_satisfied() {
        // If -u_xx == source then residual == 0
        let u_xx = -1.0;
        let source = 1.0;
        let r = pinn_residual(0.0, u_xx, source);
        assert!(r.abs() < 1e-12);
    }

    #[test]
    fn pinn_boundary_loss_zero_for_equal() {
        let v = [1.0, 0.0, -1.0];
        assert!((pinn_boundary_loss(&v, &v) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn pinn_boundary_loss_positive_for_different() {
        let u_boundary = [1.0, 2.0];
        let u_target = [0.0, 0.0];
        assert!(pinn_boundary_loss(&u_boundary, &u_target) > 0.0);
    }

    // ── PhysicsNeuralNet ─────────────────────────────────────────────────

    #[test]
    fn pinn_total_loss_formula() {
        let net = PhysicsNeuralNet::new(&[1, 4, 1], 2.0, 3.0);
        // pde_residual = 1, bc_loss = 1 => 2*1 + 3*1 = 5
        let loss = net.total_loss(1.0, 1.0);
        assert!((loss - 5.0).abs() < 1e-12);
    }

    #[test]
    fn pinn_total_loss_zero_when_both_zero() {
        let net = PhysicsNeuralNet::new(&[1, 4, 1], 1.0, 1.0);
        assert!((net.total_loss(0.0, 0.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn pinn_predict_returns_scalar() {
        let net = PhysicsNeuralNet::new(&[1, 8, 1], 1.0, 1.0);
        let _v = net.predict(0.5); // should not panic
    }

    #[test]
    fn pinn_total_loss_pde_only() {
        let net = PhysicsNeuralNet::new(&[1, 4, 1], 5.0, 0.0);
        assert!((net.total_loss(2.0, 100.0) - 10.0).abs() < 1e-12);
    }

    #[test]
    fn pinn_total_loss_bc_only() {
        let net = PhysicsNeuralNet::new(&[1, 4, 1], 0.0, 4.0);
        assert!((net.total_loss(100.0, 3.0) - 12.0).abs() < 1e-12);
    }
}
