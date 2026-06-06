// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural network-based physics acceleration (CPU mock).
//!
//! Provides a simple feed-forward neural network with multiple activation
//! functions, forward-pass inference, MSE loss, ML force potentials, and
//! collision probability prediction — all on CPU as a GPU mock backend.

// ── Activation type ──────────────────────────────────────────────────────────

/// Activation function used in a neural layer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ActivationType {
    /// Rectified linear unit: `max(0, x)`.
    Relu,
    /// Hyperbolic tangent.
    Tanh,
    /// Logistic sigmoid: `1 / (1 + exp(-x))`.
    Sigmoid,
    /// Identity / no activation.
    Linear,
}

// ── Layer and network structs ─────────────────────────────────────────────────

/// A single fully-connected neural network layer.
#[derive(Debug, Clone)]
pub struct NeuralLayer {
    /// Weight matrix: `weights[out][in]`.
    pub weights: Vec<Vec<f32>>,
    /// Bias vector, one entry per output neuron.
    pub biases: Vec<f32>,
    /// Activation function applied after the linear transform.
    pub activation: ActivationType,
}

/// A feed-forward neural network composed of stacked [`NeuralLayer`]s.
#[derive(Debug, Clone)]
pub struct NeuralNet {
    /// The ordered list of layers.
    pub layers: Vec<NeuralLayer>,
    /// Expected size of the input vector.
    pub input_size: usize,
    /// Size of the network's output.
    pub output_size: usize,
}

// ── Activation functions ──────────────────────────────────────────────────────

/// Evaluate activation function `act` at scalar `x`.
pub fn activate(x: f32, act: &ActivationType) -> f32 {
    match act {
        ActivationType::Relu => x.max(0.0),
        ActivationType::Tanh => x.tanh(),
        ActivationType::Sigmoid => 1.0 / (1.0 + (-x).exp()),
        ActivationType::Linear => x,
    }
}

/// Evaluate the derivative of activation function `act` at *pre-activation* `x`.
pub fn activate_derivative(x: f32, act: &ActivationType) -> f32 {
    match act {
        ActivationType::Relu => {
            if x > 0.0 {
                1.0
            } else {
                0.0
            }
        }
        ActivationType::Tanh => {
            let t = x.tanh();
            1.0 - t * t
        }
        ActivationType::Sigmoid => {
            let s = 1.0 / (1.0 + (-x).exp());
            s * (1.0 - s)
        }
        ActivationType::Linear => 1.0,
    }
}

// ── Inference ─────────────────────────────────────────────────────────────────

/// Run a forward pass through `net`, returning the output vector.
///
/// # Panics
/// Panics in debug mode if `input.len() != net.input_size`.
pub fn forward_pass(net: &NeuralNet, input: &[f32]) -> Vec<f32> {
    debug_assert_eq!(input.len(), net.input_size);
    let mut current: Vec<f32> = input.to_vec();
    for layer in &net.layers {
        let n_out = layer.biases.len();
        let mut next = Vec::with_capacity(n_out);
        for o in 0..n_out {
            let mut sum = layer.biases[o];
            for (i, &inp) in current.iter().enumerate() {
                if i < layer.weights[o].len() {
                    sum += layer.weights[o][i] * inp;
                }
            }
            next.push(activate(sum, &layer.activation));
        }
        current = next;
    }
    current
}

// ── Loss ──────────────────────────────────────────────────────────────────────

/// Mean squared error between `predicted` and `target` vectors.
///
/// Returns 0 if either slice is empty or lengths differ.
pub fn mse_loss(predicted: &[f32], target: &[f32]) -> f32 {
    if predicted.is_empty() || predicted.len() != target.len() {
        return 0.0;
    }
    let n = predicted.len() as f32;
    predicted
        .iter()
        .zip(target.iter())
        .map(|(p, t)| (p - t) * (p - t))
        .sum::<f32>()
        / n
}

// ── Physics applications ──────────────────────────────────────────────────────

/// Predict interatomic forces using a neural network potential.
///
/// For each atom, concatenates its position `[x, y, z]` with its type index,
/// runs a forward pass, and interprets the first three output components as the
/// predicted force `[fx, fy, fz]`.
pub fn neural_force_prediction(
    net: &NeuralNet,
    positions: &[[f32; 3]],
    types: &[u32],
) -> Vec<[f32; 3]> {
    positions
        .iter()
        .zip(types.iter())
        .map(|(pos, &atom_type)| {
            let mut inp = Vec::with_capacity(net.input_size);
            inp.push(pos[0]);
            inp.push(pos[1]);
            inp.push(pos[2]);
            inp.push(atom_type as f32);
            // Pad or truncate to net.input_size
            inp.resize(net.input_size, 0.0);
            let out = forward_pass(net, &inp);
            let fx = out.first().copied().unwrap_or(0.0);
            let fy = out.get(1).copied().unwrap_or(0.0);
            let fz = out.get(2).copied().unwrap_or(0.0);
            [fx, fy, fz]
        })
        .collect()
}

/// Predict collision probability between two spheres using a neural network.
///
/// Input features: relative displacement `[dx, dy, dz]`, radii `[ra, rb]`.
/// Returns a scalar in `[0, 1]`.
pub fn neural_collision_check(
    net: &NeuralNet,
    pos_a: [f32; 3],
    pos_b: [f32; 3],
    radii: [f32; 2],
) -> f32 {
    let dx = pos_b[0] - pos_a[0];
    let dy = pos_b[1] - pos_a[1];
    let dz = pos_b[2] - pos_a[2];
    let mut inp = vec![dx, dy, dz, radii[0], radii[1]];
    inp.resize(net.input_size, 0.0);
    let out = forward_pass(net, &inp);
    // Clamp output to [0, 1]
    out.first().copied().unwrap_or(0.0).clamp(0.0, 1.0)
}

/// Run a batched GPU-style forward pass for multiple input vectors.
pub fn gpu_neural_batch_forward(net: &NeuralNet, batch: &[Vec<f32>]) -> Vec<Vec<f32>> {
    batch.iter().map(|inp| forward_pass(net, inp)).collect()
}

// ── Network construction ──────────────────────────────────────────────────────

/// Create a fully-connected network with the given layer sizes and random weights.
///
/// `layer_sizes` must contain at least 2 entries (input + output).
/// All hidden layers use `activation`; the output layer uses `Linear`.
pub fn create_network(layer_sizes: &[usize], activation: ActivationType) -> NeuralNet {
    use rand::RngExt;
    assert!(
        layer_sizes.len() >= 2,
        "Need at least input and output sizes"
    );

    let mut rng = rand::rng();
    let mut layers = Vec::new();

    for i in 0..layer_sizes.len() - 1 {
        let n_in = layer_sizes[i];
        let n_out = layer_sizes[i + 1];
        let is_last = i == layer_sizes.len() - 2;
        let act = if is_last {
            ActivationType::Linear
        } else {
            activation
        };

        let scale = (2.0_f32 / n_in as f32).sqrt();
        let weights: Vec<Vec<f32>> = (0..n_out)
            .map(|_| (0..n_in).map(|_| rng.random_range(-scale..scale)).collect())
            .collect();
        let biases: Vec<f32> = (0..n_out).map(|_| 0.0_f32).collect();
        layers.push(NeuralLayer {
            weights,
            biases,
            activation: act,
        });
    }

    NeuralNet {
        input_size: layer_sizes[0],
        output_size: *layer_sizes.last().expect("collection should not be empty"),
        layers,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_net() -> NeuralNet {
        // 2 → 3 → 1
        create_network(&[2, 3, 1], ActivationType::Relu)
    }

    #[test]
    fn test_activate_relu_positive() {
        assert!((activate(2.0, &ActivationType::Relu) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_activate_relu_negative() {
        assert!((activate(-1.0, &ActivationType::Relu)).abs() < 1e-6);
    }

    #[test]
    fn test_activate_relu_zero() {
        assert!((activate(0.0, &ActivationType::Relu)).abs() < 1e-6);
    }

    #[test]
    fn test_activate_tanh_zero() {
        assert!((activate(0.0, &ActivationType::Tanh)).abs() < 1e-6);
    }

    #[test]
    fn test_activate_sigmoid_zero() {
        assert!((activate(0.0, &ActivationType::Sigmoid) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_activate_linear() {
        assert!((activate(3.125, &ActivationType::Linear) - 3.125).abs() < 1e-6);
    }

    #[test]
    fn test_activate_derivative_relu_positive() {
        assert!((activate_derivative(1.0, &ActivationType::Relu) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_activate_derivative_relu_negative() {
        assert!((activate_derivative(-1.0, &ActivationType::Relu)).abs() < 1e-6);
    }

    #[test]
    fn test_activate_derivative_tanh_zero() {
        assert!((activate_derivative(0.0, &ActivationType::Tanh) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_activate_derivative_sigmoid_zero() {
        assert!((activate_derivative(0.0, &ActivationType::Sigmoid) - 0.25).abs() < 1e-5);
    }

    #[test]
    fn test_activate_derivative_linear() {
        assert!((activate_derivative(99.0, &ActivationType::Linear) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_mse_loss_zero() {
        let a = vec![1.0, 2.0, 3.0];
        assert!((mse_loss(&a, &a)).abs() < 1e-6);
    }

    #[test]
    fn test_mse_loss_known() {
        let p = vec![0.0, 0.0];
        let t = vec![1.0, 1.0];
        assert!((mse_loss(&p, &t) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_mse_loss_empty() {
        assert!((mse_loss(&[], &[])).abs() < 1e-6);
    }

    #[test]
    fn test_mse_loss_length_mismatch() {
        assert!((mse_loss(&[1.0], &[1.0, 2.0])).abs() < 1e-6);
    }

    #[test]
    fn test_create_network_sizes() {
        let net = create_network(&[4, 8, 8, 3], ActivationType::Relu);
        assert_eq!(net.input_size, 4);
        assert_eq!(net.output_size, 3);
        assert_eq!(net.layers.len(), 3);
    }

    #[test]
    fn test_create_network_layer_dims() {
        let net = create_network(&[3, 5, 2], ActivationType::Tanh);
        assert_eq!(net.layers[0].weights.len(), 5);
        assert_eq!(net.layers[0].weights[0].len(), 3);
        assert_eq!(net.layers[1].weights.len(), 2);
        assert_eq!(net.layers[1].weights[0].len(), 5);
    }

    #[test]
    fn test_create_network_output_activation_linear() {
        let net = create_network(&[2, 4, 1], ActivationType::Relu);
        assert_eq!(
            net.layers.last().unwrap().activation,
            ActivationType::Linear
        );
    }

    #[test]
    fn test_forward_pass_output_size() {
        let net = simple_net();
        let out = forward_pass(&net, &[0.5, -0.3]);
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn test_forward_pass_deterministic() {
        let net = simple_net();
        let a = forward_pass(&net, &[1.0, 0.0]);
        let b = forward_pass(&net, &[1.0, 0.0]);
        assert_eq!(a, b);
    }

    #[test]
    fn test_forward_pass_zero_input() {
        let net = simple_net();
        let out = forward_pass(&net, &[0.0, 0.0]);
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn test_forward_pass_sigmoid_net() {
        let net = create_network(&[2, 2, 1], ActivationType::Sigmoid);
        let out = forward_pass(&net, &[0.0, 0.0]);
        // Output of sigmoid net on zero input should be in range [0,1] roughly
        assert!(out[0].is_finite());
    }

    #[test]
    fn test_neural_force_prediction_shape() {
        let net = create_network(&[4, 8, 3], ActivationType::Relu);
        let positions = vec![[1.0_f32, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let types = vec![0u32, 1];
        let forces = neural_force_prediction(&net, &positions, &types);
        assert_eq!(forces.len(), 2);
    }

    #[test]
    fn test_neural_force_prediction_finite() {
        let net = create_network(&[4, 6, 3], ActivationType::Tanh);
        let positions = vec![[0.0_f32; 3]];
        let types = vec![0u32];
        let forces = neural_force_prediction(&net, &positions, &types);
        assert!(forces[0][0].is_finite());
        assert!(forces[0][1].is_finite());
        assert!(forces[0][2].is_finite());
    }

    #[test]
    fn test_neural_collision_check_range() {
        let net = create_network(&[5, 4, 1], ActivationType::Sigmoid);
        let prob = neural_collision_check(&net, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5]);
        assert!((0.0..=1.0).contains(&prob));
    }

    #[test]
    fn test_neural_collision_check_zero_sep() {
        let net = create_network(&[5, 4, 1], ActivationType::Sigmoid);
        let prob = neural_collision_check(&net, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0]);
        assert!((0.0..=1.0).contains(&prob));
    }

    #[test]
    fn test_gpu_neural_batch_forward_shape() {
        let net = create_network(&[3, 4, 2], ActivationType::Relu);
        let batch: Vec<Vec<f32>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![0.0, 0.0, 0.0],
            vec![-1.0, 0.5, 0.1],
        ];
        let results = gpu_neural_batch_forward(&net, &batch);
        assert_eq!(results.len(), 3);
        for r in &results {
            assert_eq!(r.len(), 2);
        }
    }

    #[test]
    fn test_gpu_neural_batch_forward_empty() {
        let net = create_network(&[2, 2, 1], ActivationType::Linear);
        let results = gpu_neural_batch_forward(&net, &[]);
        assert!(results.is_empty());
    }

    #[test]
    fn test_create_network_two_layers() {
        let net = create_network(&[1, 1], ActivationType::Linear);
        assert_eq!(net.layers.len(), 1);
        assert_eq!(net.input_size, 1);
        assert_eq!(net.output_size, 1);
    }

    #[test]
    fn test_network_weights_finite() {
        let net = create_network(&[5, 10, 3], ActivationType::Relu);
        for layer in &net.layers {
            for row in &layer.weights {
                for &w in row {
                    assert!(w.is_finite());
                }
            }
        }
    }

    #[test]
    fn test_forward_pass_tanh_bounded() {
        let net = create_network(&[2, 4, 1], ActivationType::Tanh);
        let out = forward_pass(&net, &[100.0, -100.0]);
        // tanh saturates; linear output should still be finite
        assert!(out[0].is_finite());
    }

    #[test]
    fn test_mse_loss_asymmetric() {
        let p = vec![2.0_f32, 0.0];
        let t = vec![0.0_f32, 2.0];
        // (4 + 4) / 2 = 4
        assert!((mse_loss(&p, &t) - 4.0).abs() < 1e-5);
    }

    #[test]
    fn test_neural_force_empty_input() {
        let net = create_network(&[4, 4, 3], ActivationType::Linear);
        let forces = neural_force_prediction(&net, &[], &[]);
        assert!(forces.is_empty());
    }

    #[test]
    fn test_batch_forward_single_item() {
        let net = create_network(&[2, 3, 1], ActivationType::Relu);
        let batch = vec![vec![0.5_f32, -0.5]];
        let out = gpu_neural_batch_forward(&net, &batch);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].len(), 1);
    }
}
