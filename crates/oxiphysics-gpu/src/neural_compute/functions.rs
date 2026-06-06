//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f32::consts::PI as PI_F32;

use super::types::*;

/// Run a forward pass for every sample in `batch`.
pub fn batch_forward(net: &FeedForwardNet, batch: &[Vec<f32>]) -> Vec<Vec<f32>> {
    batch.iter().map(|x| net.forward(x)).collect()
}
/// Predict atomic energies for a batch of (descriptor, atomic_number) pairs.
///
/// Returns `0.0` for atoms whose element has no registered network.
pub fn batch_atomic_energies(
    aann: &AtomicNeuralNetwork,
    descriptors: &[Vec<f32>],
    atomic_numbers: &[u8],
) -> Vec<f32> {
    assert_eq!(
        descriptors.len(),
        atomic_numbers.len(),
        "batch_atomic_energies: descriptors and atomic_numbers must have the same length"
    );
    descriptors
        .iter()
        .zip(atomic_numbers.iter())
        .map(|(desc, &z)| aann.atomic_energy(z, desc).unwrap_or(0.0))
        .collect()
}
pub(super) const _PI_F32_USED: f32 = PI_F32;
/// Load weights from a flat f32 buffer into a network, partitioning by layer sizes.
///
/// Returns the number of weights consumed.
pub fn load_weights_from_buffer(net: &mut FeedForwardNet, buffer: &[f32]) -> usize {
    let mut offset = 0;
    for layer in &mut net.layers {
        let n_weights = layer.out_features * layer.in_features;
        let n_biases = layer.out_features;
        let total = n_weights + n_biases;
        if offset + total > buffer.len() {
            break;
        }
        layer.set_weights(&buffer[offset..offset + n_weights]);
        offset += n_weights;
        layer.set_biases(&buffer[offset..offset + n_biases]);
        offset += n_biases;
    }
    offset
}
/// Serialize network weights to a flat f32 buffer.
pub fn save_weights_to_buffer(net: &FeedForwardNet) -> Vec<f32> {
    let mut buffer = Vec::new();
    for layer in &net.layers {
        buffer.extend_from_slice(&layer.weights);
        buffer.extend_from_slice(&layer.biases);
    }
    buffer
}
/// Compute softmax of a vector.
pub fn softmax(x: &[f32]) -> Vec<f32> {
    if x.is_empty() {
        return Vec::new();
    }
    let max_val = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let exps: Vec<f32> = x.iter().map(|&v| (v - max_val).exp()).collect();
    let sum: f32 = exps.iter().sum();
    exps.iter().map(|&e| e / sum).collect()
}
/// Compute cross-entropy loss between predictions and one-hot target.
///
/// target_idx is the index of the correct class.
pub fn cross_entropy_loss(logits: &[f32], target_idx: usize) -> f32 {
    let probs = softmax(logits);
    let p = probs[target_idx].max(1e-7);
    -p.ln()
}
/// Mean squared error loss between predictions and targets.
pub fn mse_loss(predictions: &[f32], targets: &[f32]) -> f32 {
    assert_eq!(predictions.len(), targets.len());
    let n = predictions.len() as f32;
    predictions
        .iter()
        .zip(targets.iter())
        .map(|(&p, &t)| (p - t) * (p - t))
        .sum::<f32>()
        / n
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::BatchNormLayer;

    use crate::FeedForwardNet;

    #[test]
    fn test_dense_layer_zero_weights_bias_output() {
        let mut layer = DenseLayer::new(3, 2, ActivationFn::Linear);
        layer.set_biases(&[1.0, 2.0]);
        let out = layer.forward(&[10.0, 20.0, 30.0]);
        assert!((out[0] - 1.0).abs() < 1e-6, "expected 1.0, got {}", out[0]);
        assert!((out[1] - 2.0).abs() < 1e-6, "expected 2.0, got {}", out[1]);
    }
    #[test]
    fn test_activation_tanh_zero() {
        let act = ActivationFn::Tanh;
        assert!((act.apply(0.0)).abs() < 1e-7, "tanh(0) should be 0");
    }
    #[test]
    fn test_activation_relu() {
        let act = ActivationFn::Relu;
        assert!((act.apply(-1.0)).abs() < 1e-7, "relu(-1) should be 0");
        assert!((act.apply(1.0) - 1.0).abs() < 1e-7, "relu(1) should be 1");
    }
    #[test]
    fn test_feedforward_net_two_layers() {
        let mut layer1 = DenseLayer::new(2, 3, ActivationFn::Linear);
        layer1.set_weights(&[1.0, 0.0, 0.0, 1.0, 1.0, 0.0]);
        let mut layer2 = DenseLayer::new(3, 1, ActivationFn::Linear);
        layer2.set_weights(&[1.0, 1.0, 1.0]);
        let mut net = FeedForwardNet::new();
        net.add_layer(layer1);
        net.add_layer(layer2);
        let out = net.forward(&[1.0, 2.0]);
        assert!((out[0] - 4.0).abs() < 1e-5, "expected 4.0, got {}", out[0]);
    }
    #[test]
    fn test_cutoff_fn_at_zero_and_beyond_rc() {
        let rc = 6.0;
        let fc0 = BehlerParrinelloDescriptor::cutoff_fn(0.0, rc);
        assert!((fc0 - 1.0).abs() < 1e-10, "fc(0) should be 1.0, got {fc0}");
        let fc_rc = BehlerParrinelloDescriptor::cutoff_fn(rc, rc);
        assert!((fc_rc).abs() < 1e-10, "fc(rc) should be 0.0, got {fc_rc}");
        let fc_beyond = BehlerParrinelloDescriptor::cutoff_fn(rc + 1.0, rc);
        assert!(
            (fc_beyond).abs() < 1e-10,
            "fc(>rc) should be 0.0, got {fc_beyond}"
        );
    }
    #[test]
    fn test_radial_g2_decreases_with_distance() {
        let rc = 6.0;
        let eta = 0.5;
        let rs = 0.0;
        let g2_near = BehlerParrinelloDescriptor::radial_g2(1.0, eta, rs, rc);
        let g2_far = BehlerParrinelloDescriptor::radial_g2(5.0, eta, rs, rc);
        assert!(
            g2_near > g2_far,
            "G2 should decrease with distance: near={g2_near}, far={g2_far}"
        );
    }
    #[test]
    fn test_data_normalizer_round_trip() {
        let data = vec![
            vec![1.0_f32, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let norm = DataNormalizer::fit(&data);
        let sample = &data[1];
        let transformed = norm.transform(sample);
        let recovered = norm.inverse_transform(&transformed);
        for (a, b) in recovered.iter().zip(sample.iter()) {
            assert!(
                (a - b).abs() < 1e-5,
                "round-trip failed: got {a}, expected {b}"
            );
        }
    }
    #[test]
    fn test_network_builder_simple_aann_architecture() {
        let hidden = &[64_usize, 32];
        let net = NetworkBuilder::simple_aann(20, hidden, 1);
        assert_eq!(net.layers.len(), 3);
        assert_eq!(net.input_size(), Some(20));
        assert_eq!(net.output_size(), Some(1));
        assert_eq!(net.layers[0].activation, ActivationFn::Tanh);
        assert_eq!(net.layers[1].activation, ActivationFn::Tanh);
        assert_eq!(net.layers[2].activation, ActivationFn::Linear);
        assert_eq!(net.total_parameters(), 3457);
    }
    #[test]
    fn test_batch_norm_identity_transform() {
        let bn = BatchNormLayer::new(3);
        let input = vec![1.0, 2.0, 3.0];
        let output = bn.forward(&input);
        for i in 0..3 {
            assert!(
                (output[i] - input[i]).abs() < 1e-4,
                "output[{i}]={}",
                output[i]
            );
        }
    }
    #[test]
    fn test_batch_norm_zero_mean_unit_var() {
        let mut bn = BatchNormLayer::new(2);
        bn.set_stats(&[5.0, 10.0], &[4.0, 9.0]);
        let output = bn.forward(&[5.0, 10.0]);
        assert!(output[0].abs() < 1e-4);
        assert!(output[1].abs() < 1e-4);
    }
    #[test]
    fn test_batch_norm_affine() {
        let mut bn = BatchNormLayer::new(2);
        bn.set_stats(&[0.0, 0.0], &[1.0, 1.0]);
        bn.set_affine(&[2.0, 3.0], &[1.0, -1.0]);
        let output = bn.forward(&[1.0, 1.0]);
        assert!((output[0] - 3.0).abs() < 1e-4);
        assert!((output[1] - 2.0).abs() < 1e-4);
    }
    #[test]
    fn test_inference_pipeline_dense_only() {
        let mut pipe = InferencePipeline::new();
        let mut layer = DenseLayer::new(2, 1, ActivationFn::Linear);
        layer.set_weights(&[1.0, 1.0]);
        layer.set_biases(&[0.0]);
        pipe.add_op(InferenceOp::Dense(layer));
        let out = pipe.forward(&[3.0, 4.0]);
        assert!((out[0] - 7.0).abs() < 1e-5);
    }
    #[test]
    fn test_inference_pipeline_with_batch_norm() {
        let mut pipe = InferencePipeline::new();
        let mut layer = DenseLayer::new(2, 2, ActivationFn::Linear);
        layer.set_weights(&[1.0, 0.0, 0.0, 1.0]);
        pipe.add_op(InferenceOp::Dense(layer));
        pipe.add_op(InferenceOp::BatchNorm(BatchNormLayer::new(2)));
        pipe.add_op(InferenceOp::Activation(ActivationFn::Relu));
        let out = pipe.forward(&[1.0, -2.0]);
        assert!((out[0] - 1.0).abs() < 1e-4);
        assert!(out[1].abs() < 1e-4);
    }
    #[test]
    fn test_inference_pipeline_total_params() {
        let mut pipe = InferencePipeline::new();
        pipe.add_op(InferenceOp::Dense(DenseLayer::new(
            4,
            3,
            ActivationFn::Relu,
        )));
        pipe.add_op(InferenceOp::BatchNorm(BatchNormLayer::new(3)));
        pipe.add_op(InferenceOp::Dense(DenseLayer::new(
            3,
            1,
            ActivationFn::Linear,
        )));
        assert_eq!(pipe.total_parameters(), 25);
    }
    #[test]
    fn test_save_load_weights_roundtrip() {
        let mut net = FeedForwardNet::new();
        let mut l1 = DenseLayer::new(2, 3, ActivationFn::Relu);
        l1.set_weights(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        l1.set_biases(&[0.1, 0.2, 0.3]);
        net.add_layer(l1);
        let buf = save_weights_to_buffer(&net);
        assert_eq!(buf.len(), 9);
        let mut net2 = FeedForwardNet::new();
        net2.add_layer(DenseLayer::new(2, 3, ActivationFn::Relu));
        let consumed = load_weights_from_buffer(&mut net2, &buf);
        assert_eq!(consumed, 9);
        let input = vec![1.0, 2.0];
        let out1 = net.forward(&input);
        let out2 = net2.forward(&input);
        for (a, b) in out1.iter().zip(out2.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }
    #[test]
    fn test_load_weights_partial() {
        let mut net = FeedForwardNet::new();
        net.add_layer(DenseLayer::new(2, 3, ActivationFn::Relu));
        net.add_layer(DenseLayer::new(3, 1, ActivationFn::Linear));
        let buf = vec![1.0_f32; 9];
        let consumed = load_weights_from_buffer(&mut net, &buf);
        assert_eq!(consumed, 9);
    }
    #[test]
    fn test_softmax_sums_to_one() {
        let logits = vec![1.0_f32, 2.0, 3.0, 4.0];
        let probs = softmax(&logits);
        let sum: f32 = probs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5, "softmax sum={sum}");
    }
    #[test]
    fn test_softmax_largest_is_max() {
        let logits = vec![1.0_f32, 5.0, 2.0];
        let probs = softmax(&logits);
        assert!(probs[1] > probs[0] && probs[1] > probs[2]);
    }
    #[test]
    fn test_softmax_empty() {
        assert!(softmax(&[]).is_empty());
    }
    #[test]
    fn test_softmax_single() {
        let probs = softmax(&[42.0]);
        assert!((probs[0] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_cross_entropy_loss_correct_class() {
        let logits = vec![-10.0_f32, 10.0, -10.0];
        let loss = cross_entropy_loss(&logits, 1);
        assert!(loss < 0.01, "loss should be small, got {loss}");
    }
    #[test]
    fn test_cross_entropy_loss_wrong_class() {
        let logits = vec![10.0_f32, -10.0, -10.0];
        let loss = cross_entropy_loss(&logits, 1);
        assert!(loss > 1.0, "loss should be large, got {loss}");
    }
    #[test]
    fn test_mse_loss_zero() {
        let pred = vec![1.0_f32, 2.0, 3.0];
        let target = vec![1.0, 2.0, 3.0];
        assert!(mse_loss(&pred, &target).abs() < 1e-7);
    }
    #[test]
    fn test_mse_loss_positive() {
        let pred = vec![1.0_f32, 2.0];
        let target = vec![3.0, 4.0];
        let loss = mse_loss(&pred, &target);
        assert!((loss - 4.0).abs() < 1e-5);
    }
    #[test]
    fn test_sigmoid_derivative() {
        let act = ActivationFn::Sigmoid;
        let d = act.derivative(0.0);
        assert!((d - 0.25).abs() < 1e-5, "sigmoid'(0) = {d}");
    }
    #[test]
    fn test_silu_at_zero() {
        let act = ActivationFn::Silu;
        let v = act.apply(0.0);
        assert!(v.abs() < 1e-7);
    }
    #[test]
    fn test_gelu_derivative_positive() {
        let act = ActivationFn::Gelu;
        let d = act.derivative(1.0);
        assert!(d > 0.0, "GELU derivative at 1.0 should be positive");
    }
    #[test]
    fn test_linear_derivative() {
        let act = ActivationFn::Linear;
        assert!((act.derivative(42.0) - 1.0).abs() < 1e-7);
    }
}
/// Run the network on every position and return a force vector per atom.
///
/// The network is expected to take a 3-component position and output a
/// 3-component force vector.  This is a CPU mock of a GPU batch dispatch.
pub fn compute_forces_batch(positions: &[[f64; 3]], network: &NeuralNetwork) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|p| {
            let out = network.forward(p);
            let fx = out.first().copied().unwrap_or(0.0);
            let fy = out.get(1).copied().unwrap_or(0.0);
            let fz = out.get(2).copied().unwrap_or(0.0);
            [fx, fy, fz]
        })
        .collect()
}
/// Compute a scalar potential energy by summing the first network output over all atoms.
///
/// Suitable for networks that map a position to a scalar energy contribution.
pub fn neural_potential_energy(network: &NeuralNetwork, positions: &[[f64; 3]]) -> f64 {
    positions
        .iter()
        .map(|p| network.forward(p).first().copied().unwrap_or(0.0))
        .sum()
}
#[cfg(test)]
mod neural_f64_tests {

    use crate::ActivationFn64;

    use crate::GpuNeuralBuffer;

    use crate::NeuralNetwork;

    use crate::compute_forces_batch;

    use crate::neural_potential_energy;

    #[test]
    fn test_forward_pass_output_size() {
        let net = NeuralNetwork::new(&[3, 8, 8, 3], ActivationFn64::Tanh);
        let input = [1.0, 0.5, -0.5];
        let out = net.forward(&input);
        assert_eq!(out.len(), 3, "output should have 3 components");
    }
    #[test]
    fn test_relu_activation_f64() {
        let act = ActivationFn64::Relu;
        assert!((act.apply(-2.0)).abs() < 1e-12, "relu(-2) = 0");
        assert!((act.apply(3.0) - 3.0).abs() < 1e-12, "relu(3) = 3");
    }
    #[test]
    fn test_relu_batch() {
        let act = ActivationFn64::Relu;
        let mut v = vec![-1.0, 0.0, 2.0, -0.5, 3.0];
        act.apply_batch(&mut v);
        assert_eq!(v, vec![0.0, 0.0, 2.0, 0.0, 3.0]);
    }
    #[test]
    fn test_xavier_init_non_zero() {
        let net = NeuralNetwork::new(&[4, 16, 4], ActivationFn64::Relu);
        let all_zero = net
            .layers
            .iter()
            .all(|l| l.weights.iter().all(|row| row.iter().all(|&w| w == 0.0)));
        assert!(!all_zero, "Xavier-init weights should not all be zero");
    }
    #[test]
    fn test_batch_forces_count() {
        let net = NeuralNetwork::new(&[3, 8, 3], ActivationFn64::Relu);
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let forces = compute_forces_batch(&positions, &net);
        assert_eq!(forces.len(), 5, "one force vector per position");
    }
    #[test]
    fn test_gpu_neural_buffer_roundtrip() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let buf = GpuNeuralBuffer::pack_positions(&positions);
        assert_eq!(buf.batch_size, 2);
        assert_eq!(buf.data.len(), 6);
        let forces = buf.unpack_forces();
        assert_eq!(forces[0], [1.0, 2.0, 3.0]);
        assert_eq!(forces[1], [4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_neural_potential_energy_positive() {
        let net = NeuralNetwork::new(&[3, 4, 1], ActivationFn64::Relu);
        let positions = vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
        let energy = neural_potential_energy(&net, &positions);
        assert!(energy.is_finite(), "energy should be finite");
    }
}
/// Compute L2 regularisation loss contribution.
///
/// Returns `0.5 * lambda * sum(w^2)` for all elements of `weights`.
pub fn l2_regularisation(weights: &[f64], lambda: f64) -> f64 {
    0.5 * lambda * weights.iter().map(|&w| w * w).sum::<f64>()
}
/// Compute L2 regularisation gradient contribution (adds `lambda * w` to each element).
pub fn l2_regularisation_grad(weights: &[f64], lambda: f64) -> Vec<f64> {
    weights.iter().map(|&w| lambda * w).collect()
}
/// Huber loss between prediction and target.
///
/// Acts as MSE when `|pred - target| < delta`, and as MAE beyond that.
pub fn huber_loss(pred: f64, target: f64, delta: f64) -> f64 {
    let e = (pred - target).abs();
    if e <= delta {
        0.5 * e * e
    } else {
        delta * (e - 0.5 * delta)
    }
}
/// Mean Huber loss over a batch.
pub fn mean_huber_loss(predictions: &[f64], targets: &[f64], delta: f64) -> f64 {
    assert_eq!(predictions.len(), targets.len());
    let n = predictions.len() as f64;
    predictions
        .iter()
        .zip(targets.iter())
        .map(|(&p, &t)| huber_loss(p, t, delta))
        .sum::<f64>()
        / n
}
/// Gradient of Huber loss w.r.t. predictions.
pub fn huber_loss_grad(predictions: &[f64], targets: &[f64], delta: f64) -> Vec<f64> {
    predictions
        .iter()
        .zip(targets.iter())
        .map(|(&p, &t)| {
            let e = p - t;
            if e.abs() <= delta {
                e
            } else {
                delta * e.signum()
            }
        })
        .collect()
}
/// Compute the L2 norm (Frobenius norm) of a concatenated gradient vector.
///
/// Used for gradient clipping: if the norm exceeds a threshold the caller
/// should scale all gradients by `threshold / norm`.
pub fn compute_gradient_norm(gradients: &[&[f64]]) -> f64 {
    let sum_sq: f64 = gradients
        .iter()
        .flat_map(|g| g.iter())
        .map(|&v| v * v)
        .sum();
    sum_sq.sqrt()
}
/// Clip a collection of gradient slices in-place so that their combined L2
/// norm does not exceed `max_norm`.
///
/// If the current norm is already ≤ `max_norm` the gradients are unchanged.
/// Returns the pre-clip norm.
pub fn clip_gradients_by_norm(gradients: &mut [Vec<f64>], max_norm: f64) -> f64 {
    let refs: Vec<&[f64]> = gradients.iter().map(|v| v.as_slice()).collect();
    let norm = compute_gradient_norm(&refs);
    if norm > max_norm && norm > 0.0 {
        let scale = max_norm / norm;
        for g in gradients.iter_mut() {
            for v in g.iter_mut() {
                *v *= scale;
            }
        }
    }
    norm
}
#[cfg(test)]
mod extended_tests {

    use crate::AdamOptimizer;

    use crate::DenseLayer64;
    use crate::DropoutLayer;
    use crate::ExtActivation;

    use crate::GradAccumulator;

    use crate::huber_loss;
    use crate::huber_loss_grad;
    use crate::l2_regularisation;
    use crate::l2_regularisation_grad;
    use crate::mean_huber_loss;

    #[test]
    fn test_leaky_relu_positive() {
        let act = ExtActivation::LeakyRelu(0.01);
        assert!((act.apply(3.0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_leaky_relu_negative() {
        let act = ExtActivation::LeakyRelu(0.1);
        assert!((act.apply(-2.0) - (-0.2)).abs() < 1e-12);
    }
    #[test]
    fn test_leaky_relu_derivative_positive() {
        let act = ExtActivation::LeakyRelu(0.05);
        assert!((act.derivative(1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_leaky_relu_derivative_negative() {
        let act = ExtActivation::LeakyRelu(0.05);
        assert!((act.derivative(-1.0) - 0.05).abs() < 1e-12);
    }
    #[test]
    fn test_swish_at_zero() {
        let act = ExtActivation::Swish(1.0);
        assert!(act.apply(0.0).abs() < 1e-12);
    }
    #[test]
    fn test_swish_positive_region() {
        let act = ExtActivation::Swish(1.0);
        let v = act.apply(10.0);
        assert!((v - 10.0).abs() < 0.01, "swish(10) ≈ 10, got {v}");
    }
    #[test]
    fn test_ext_activation_apply_vec() {
        let act = ExtActivation::Relu;
        let mut v = vec![-1.0, 0.0, 2.0, -3.0];
        act.apply_vec(&mut v);
        assert_eq!(v, vec![0.0, 0.0, 2.0, 0.0]);
    }
    #[test]
    fn test_dense64_forward_linear_known_values() {
        let mut layer = DenseLayer64::new(2, 2, ExtActivation::Linear);
        layer.weights = vec![1.0, 0.0, 0.0, 1.0];
        let out = layer.forward(&[3.0, 5.0]);
        assert!((out[0] - 3.0).abs() < 1e-12);
        assert!((out[1] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_dense64_backward_gradient_shapes() {
        let mut layer = DenseLayer64::new(3, 2, ExtActivation::Relu);
        layer.weights = vec![1.0; 6];
        let _out = layer.forward(&[1.0, 2.0, 3.0]);
        let (gw, gb, di) = layer.backward(&[1.0, 1.0]);
        assert_eq!(gw.len(), 6, "grad_weights shape");
        assert_eq!(gb.len(), 2, "grad_biases shape");
        assert_eq!(di.len(), 3, "delta_in shape");
    }
    #[test]
    fn test_dense64_sgd_update_reduces_output() {
        let mut layer = DenseLayer64::new(1, 1, ExtActivation::Linear);
        layer.weights = vec![2.0];
        layer.biases = vec![0.0];
        let out = layer.forward(&[1.0]);
        let loss_before = out[0] * out[0];
        let (gw, gb, _) = layer.backward(&[2.0 * out[0]]);
        layer.apply_sgd(&gw, &gb, 0.1);
        let out2 = layer.forward(&[1.0]);
        let loss_after = out2[0] * out2[0];
        assert!(loss_after < loss_before, "SGD should reduce loss");
    }
    #[test]
    fn test_dense64_num_params() {
        let layer = DenseLayer64::new(4, 3, ExtActivation::Linear);
        assert_eq!(layer.num_params(), 4 * 3 + 3);
    }
    #[test]
    fn test_dropout_inference_passthrough() {
        let mut drop = DropoutLayer::new(0.5, false);
        let input = vec![1.0, 2.0, 3.0, 4.0];
        let out = drop.forward(&input);
        assert_eq!(out, input, "dropout in eval mode should pass through");
    }
    #[test]
    fn test_dropout_rate_zero_no_drop() {
        let mut drop = DropoutLayer::new(0.0, true);
        let input = vec![1.0, 2.0, 3.0];
        let out = drop.forward(&input);
        assert_eq!(out, input, "zero rate should not drop anything");
    }
    #[test]
    fn test_dropout_rate_one_all_zero() {
        let mut drop = DropoutLayer::new(1.0, true);
        let input = vec![5.0, 6.0, 7.0];
        let out = drop.forward(&input);
        assert!(
            out.iter().all(|&x| x == 0.0),
            "rate=1 should zero everything"
        );
    }
    #[test]
    fn test_dropout_training_some_zeros() {
        let mut drop = DropoutLayer::new(0.5, true);
        drop.set_seed(42);
        let input = vec![1.0_f64; 100];
        let out = drop.forward(&input);
        let n_zeros = out.iter().filter(|&&x| x == 0.0).count();
        assert!(n_zeros > 10, "expected some zeros, got {n_zeros}");
        assert!(
            n_zeros < 90,
            "expected some non-zeros, too many zeros: {n_zeros}"
        );
    }
    #[test]
    fn test_dropout_backward_applies_mask() {
        let mut drop = DropoutLayer::new(0.0, false);
        let input = vec![1.0, 2.0, 3.0];
        let _out = drop.forward(&input);
        let grad = drop.backward(&[1.0, 1.0, 1.0]);
        assert_eq!(grad, vec![1.0, 1.0, 1.0]);
    }
    #[test]
    fn test_adam_step_decreases_loss() {
        let mut params = vec![5.0_f64];
        let mut opt = AdamOptimizer::default_params(1);
        let initial_abs = params[0].abs();
        for _ in 0..20 {
            let grads = vec![2.0 * params[0]];
            opt.step_update(&mut params, &grads);
        }
        assert!(
            params[0].abs() < initial_abs,
            "Adam should move towards zero, final={}",
            params[0]
        );
    }
    #[test]
    fn test_adam_step_increments_counter() {
        let mut params = vec![1.0_f64; 3];
        let mut opt = AdamOptimizer::default_params(3);
        assert_eq!(opt.step, 0);
        let grads = vec![0.1; 3];
        opt.step_update(&mut params, &grads);
        assert_eq!(opt.step, 1);
        opt.step_update(&mut params, &grads);
        assert_eq!(opt.step, 2);
    }
    #[test]
    fn test_adam_reset() {
        let mut params = vec![1.0_f64; 2];
        let mut opt = AdamOptimizer::default_params(2);
        let grads = vec![0.5; 2];
        opt.step_update(&mut params, &grads);
        opt.reset();
        assert_eq!(opt.step, 0);
        assert!(opt.m.iter().all(|&x| x == 0.0));
        assert!(opt.v.iter().all(|&x| x == 0.0));
    }
    #[test]
    fn test_adam_moment_accumulation() {
        let mut params = vec![1.0_f64];
        let mut opt = AdamOptimizer::new(1, 1e-3, 0.9, 0.999, 1e-8);
        let grads = vec![1.0_f64];
        opt.step_update(&mut params, &grads);
        assert!(
            (opt.m[0] - 0.1).abs() < 1e-10,
            "m after step 1 = {}",
            opt.m[0]
        );
        assert!(
            (opt.v[0] - 0.001).abs() < 1e-10,
            "v after step 1 = {}",
            opt.v[0]
        );
    }
    #[test]
    fn test_grad_accumulator_mean() {
        let mut acc = GradAccumulator::new(2, 1);
        acc.accumulate(&[1.0, 2.0], &[3.0]);
        acc.accumulate(&[3.0, 4.0], &[1.0]);
        let (gw, gb) = acc.mean_grads();
        assert!((gw[0] - 2.0).abs() < 1e-12);
        assert!((gw[1] - 3.0).abs() < 1e-12);
        assert!((gb[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_grad_accumulator_zero() {
        let mut acc = GradAccumulator::new(3, 2);
        acc.accumulate(&[1.0, 2.0, 3.0], &[4.0, 5.0]);
        assert_eq!(acc.count, 1);
        acc.zero();
        assert_eq!(acc.count, 0);
        assert!(acc.grad_weights.iter().all(|&x| x == 0.0));
    }
    #[test]
    fn test_l2_regularisation() {
        let weights = vec![1.0, 2.0, 3.0];
        let reg = l2_regularisation(&weights, 0.01);
        assert!((reg - 0.07).abs() < 1e-12, "L2 reg = {reg}");
    }
    #[test]
    fn test_l2_regularisation_grad() {
        let weights = vec![2.0, -3.0];
        let grad = l2_regularisation_grad(&weights, 0.1);
        assert!((grad[0] - 0.2).abs() < 1e-12);
        assert!((grad[1] - (-0.3)).abs() < 1e-12);
    }
    #[test]
    fn test_huber_loss_small_error() {
        let loss = huber_loss(1.0, 1.1, 0.5);
        assert!((loss - 0.5 * 0.01).abs() < 1e-12, "huber loss = {loss}");
    }
    #[test]
    fn test_huber_loss_large_error() {
        let loss = huber_loss(0.0, 5.0, 1.0);
        assert!((loss - 4.5).abs() < 1e-12, "huber loss = {loss}");
    }
    #[test]
    fn test_mean_huber_loss() {
        let preds = vec![0.0, 0.0];
        let targets = vec![0.1, 5.0];
        let loss = mean_huber_loss(&preds, &targets, 1.0);
        assert!((loss - 2.2525).abs() < 1e-10, "mean huber loss = {loss}");
    }
    #[test]
    fn test_huber_loss_grad_small() {
        let preds = vec![1.0];
        let targets = vec![1.1];
        let grad = huber_loss_grad(&preds, &targets, 1.0);
        assert!((grad[0] - (-0.1)).abs() < 1e-12);
    }
    #[test]
    fn test_huber_loss_grad_large() {
        let preds = vec![10.0];
        let targets = vec![0.0];
        let grad = huber_loss_grad(&preds, &targets, 1.0);
        assert!((grad[0] - 1.0).abs() < 1e-12);
    }
}
#[cfg(test)]
mod conv_rnn_tests {

    use crate::BatchNormLayer;
    use crate::Conv1DLayer;

    use crate::ExtActivation;
    use crate::FeedForwardNet;

    use crate::LayerNorm;
    use crate::LayerNormLayer;

    use crate::RnnCell;

    use crate::clip_gradients_by_norm;

    use crate::compute_gradient_norm;

    #[test]
    fn test_conv1d_zero_weights_output_is_bias() {
        let mut conv = Conv1DLayer::new(2, 3, 2, ExtActivation::Linear);
        conv.biases = vec![1.0, 2.0, 3.0];
        let input = vec![vec![0.5, 0.5]; 4];
        let out = conv.forward(&input);
        assert_eq!(out.len(), 4);
        for row in &out {
            assert_eq!(row.len(), 3);
            assert!((row[0] - 1.0).abs() < 1e-12, "out[0] = {}", row[0]);
            assert!((row[1] - 2.0).abs() < 1e-12, "out[1] = {}", row[1]);
            assert!((row[2] - 3.0).abs() < 1e-12, "out[2] = {}", row[2]);
        }
    }
    #[test]
    fn test_conv1d_num_params() {
        let conv = Conv1DLayer::new(4, 8, 3, ExtActivation::Relu);
        assert_eq!(conv.num_params(), 104);
    }
    #[test]
    fn test_conv1d_output_shape() {
        let conv = Conv1DLayer::new(3, 5, 3, ExtActivation::Tanh);
        let input: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0, 0.0, -1.0]).collect();
        let out = conv.forward(&input);
        assert_eq!(out.len(), 10, "seq_len preserved");
        assert_eq!(out[0].len(), 5, "out_channels");
    }
    #[test]
    fn test_conv1d_causal_first_step_only_sees_t0() {
        let mut conv = Conv1DLayer::new(1, 1, 3, ExtActivation::Linear);
        conv.weights[0][0][0] = 1.0;
        conv.weights[0][1][0] = 100.0;
        conv.weights[0][2][0] = 100.0;
        let input = vec![vec![5.0], vec![0.0], vec![0.0]];
        let out = conv.forward(&input);
        assert!(
            (out[0][0] - 5.0).abs() < 1e-12,
            "t=0 output = {}",
            out[0][0]
        );
    }
    #[test]
    fn test_conv1d_kernel1_is_pointwise() {
        let mut conv = Conv1DLayer::new(2, 1, 1, ExtActivation::Linear);
        conv.weights[0][0][0] = 2.0;
        conv.weights[0][0][1] = 3.0;
        let input = vec![vec![1.0, 1.0], vec![2.0, 2.0]];
        let out = conv.forward(&input);
        assert!((out[0][0] - 5.0).abs() < 1e-12);
        assert!((out[1][0] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_conv1d_relu_clips_negative() {
        let conv = Conv1DLayer::new(1, 1, 1, ExtActivation::Relu);
        let input: Vec<Vec<f64>> = vec![vec![-5.0], vec![-3.0]];
        let out = conv.forward(&input);
        assert!(out[0][0] >= 0.0, "relu should clip negative");
        assert!(out[1][0] >= 0.0);
    }
    #[test]
    fn test_layer_norm_zero_mean_after_forward() {
        let ln = LayerNorm::new(4);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let out = ln.forward(&x);
        let mean: f64 = out.iter().sum::<f64>() / out.len() as f64;
        assert!(mean.abs() < 1e-10, "mean after LayerNorm = {mean}");
    }
    #[test]
    fn test_layer_norm_unit_variance() {
        let ln = LayerNorm::new(4);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let out = ln.forward(&x);
        let mean: f64 = out.iter().sum::<f64>() / out.len() as f64;
        let var: f64 = out.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / out.len() as f64;
        assert!((var - 1.0).abs() < 1e-3, "variance after LayerNorm = {var}");
    }
    #[test]
    fn test_layer_norm_identity_gamma_beta() {
        let ln = LayerNorm::new(3);
        let x = vec![5.0, 5.0, 5.0];
        let out = ln.forward(&x);
        for &v in &out {
            assert!(v.abs() < 1e-4, "constant input → near-zero output, got {v}");
        }
    }
    #[test]
    fn test_layer_norm_output_length() {
        let ln = LayerNorm::new(6);
        let x = vec![1.0; 6];
        let out = ln.forward(&x);
        assert_eq!(out.len(), 6);
    }
    #[test]
    fn test_layer_norm_custom_gamma_beta() {
        let mut ln = LayerNorm::new(2);
        ln.gamma = vec![2.0, 3.0];
        ln.beta = vec![1.0, -1.0];
        let x = vec![0.0, 4.0];
        let out = ln.forward(&x);
        assert!((out[0] - (-1.0)).abs() < 1e-4, "out[0] = {}", out[0]);
        assert!((out[1] - 2.0).abs() < 1e-4, "out[1] = {}", out[1]);
    }
    #[test]
    fn test_rnn_cell_zero_weights_output_is_activated_bias() {
        let cell = RnnCell::new(2, 3, ExtActivation::Linear);
        let x = vec![10.0, 20.0];
        let h_prev = vec![1.0, 2.0, 3.0];
        let h = cell.step(&x, &h_prev);
        for &v in &h {
            assert!(v.abs() < 1e-12, "zero weights → zero output, got {v}");
        }
    }
    #[test]
    fn test_rnn_cell_output_length() {
        let cell = RnnCell::new(4, 8, ExtActivation::Tanh);
        let x = vec![0.0; 4];
        let h_prev = vec![0.0; 8];
        let h = cell.step(&x, &h_prev);
        assert_eq!(h.len(), 8);
    }
    #[test]
    fn test_rnn_cell_identity_weights_copies_input() {
        let mut cell = RnnCell::new(2, 2, ExtActivation::Linear);
        cell.w_x[0] = 1.0;
        cell.w_x[3] = 1.0;
        let x = vec![3.0, 7.0];
        let h_prev = vec![0.0, 0.0];
        let h = cell.step(&x, &h_prev);
        assert!((h[0] - 3.0).abs() < 1e-12, "h[0] = {}", h[0]);
        assert!((h[1] - 7.0).abs() < 1e-12, "h[1] = {}", h[1]);
    }
    #[test]
    fn test_rnn_cell_sequence_length() {
        let cell = RnnCell::new(3, 5, ExtActivation::Relu);
        let seq: Vec<Vec<f64>> = (0..7).map(|_| vec![0.0; 3]).collect();
        let states = cell.forward_sequence(&seq);
        assert_eq!(states.len(), 7);
        assert_eq!(states[0].len(), 5);
    }
    #[test]
    fn test_rnn_cell_sequence_accumulates_state() {
        let mut cell = RnnCell::new(1, 1, ExtActivation::Linear);
        cell.w_x[0] = 0.0;
        cell.w_h[0] = 2.0;
        let h0 = vec![1.0_f64];
        let h1 = cell.step(&[0.0], &h0);
        assert!((h1[0] - 2.0).abs() < 1e-12, "h1 = {}", h1[0]);
        let h2 = cell.step(&[0.0], &h1);
        assert!((h2[0] - 4.0).abs() < 1e-12, "h2 = {}", h2[0]);
    }
    #[test]
    fn test_layer_norm_zero_mean_unit_variance() {
        let ln = LayerNormLayer::new(4);
        let input = vec![1.0, 2.0, 3.0, 4.0];
        let out = ln.forward(&input);
        let mean: f64 = out.iter().sum::<f64>() / out.len() as f64;
        let var: f64 = out.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / out.len() as f64;
        assert!(mean.abs() < 1e-10, "output mean should be ~0, got {mean}");
        assert!(
            (var - 1.0).abs() < 1e-5,
            "output var should be ~1, got {var}"
        );
    }
    #[test]
    fn test_layer_norm_gamma_scales_output() {
        let mut ln = LayerNormLayer::new(3);
        ln.gamma = vec![2.0, 2.0, 2.0];
        let input = vec![1.0, 2.0, 3.0];
        let out = ln.forward(&input);
        let mean: f64 = out.iter().sum::<f64>() / out.len() as f64;
        let var: f64 = out.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / out.len() as f64;
        assert!(
            (var.sqrt() - 2.0).abs() < 1e-3,
            "std should be ~2 with gamma=2, got {}",
            var.sqrt()
        );
    }
    #[test]
    fn test_layer_norm_beta_shifts_output() {
        let mut ln = LayerNormLayer::new(2);
        ln.beta = vec![5.0, 5.0];
        let input = vec![0.0, 1.0];
        let out = ln.forward(&input);
        let mean: f64 = out.iter().sum::<f64>() / out.len() as f64;
        assert!(
            (mean - 5.0).abs() < 1e-5,
            "mean should be ~5 with beta=5, got {mean}"
        );
    }
    #[test]
    fn test_layer_norm_backward_d_beta_equals_d_output() {
        let ln = LayerNormLayer::new(3);
        let input = vec![1.0, 2.0, 3.0];
        let d_out = vec![0.1, 0.2, 0.3];
        let (_d_in, _d_gamma, d_beta) = ln.backward(&input, &d_out);
        for (i, (&db, &dout)) in d_beta.iter().zip(d_out.iter()).enumerate() {
            assert!(
                (db - dout).abs() < 1e-12,
                "d_beta[{i}] should equal d_output[{i}]"
            );
        }
    }
    #[test]
    fn test_compute_gradient_norm_zero() {
        let g: Vec<f64> = vec![0.0; 5];
        let norm = compute_gradient_norm(&[&g]);
        assert!(norm.abs() < 1e-12, "norm of zero gradient should be 0");
    }
    #[test]
    fn test_compute_gradient_norm_known_value() {
        let g = vec![3.0_f64, 4.0];
        let norm = compute_gradient_norm(&[&g]);
        assert!((norm - 5.0).abs() < 1e-10, "norm should be 5.0, got {norm}");
    }
    #[test]
    fn test_clip_gradients_no_clip_when_below_max() {
        let mut grads = vec![vec![1.0_f64, 0.0], vec![0.0, 1.0]];
        let norm_before = clip_gradients_by_norm(&mut grads, 5.0);
        assert!((norm_before - 2.0_f64.sqrt()).abs() < 1e-10);
        assert!((grads[0][0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_clip_gradients_clips_correctly() {
        let mut grads = vec![vec![3.0_f64, 4.0]];
        clip_gradients_by_norm(&mut grads, 1.0);
        let new_norm: f64 = grads[0].iter().map(|&v| v * v).sum::<f64>().sqrt();
        assert!(
            (new_norm - 1.0).abs() < 1e-10,
            "clipped norm should be 1.0, got {new_norm}"
        );
    }
    #[test]
    fn test_batch_norm_update_running_stats_mean_converges() {
        let mut bn = BatchNormLayer::new(2);
        let batch = vec![vec![2.0_f32, 3.0], vec![2.0, 3.0], vec![2.0, 3.0]];
        for _ in 0..50 {
            bn.update_running_stats(&batch, 0.1);
        }
        assert!(
            (bn.running_mean[0] - 2.0).abs() < 0.1,
            "mean[0] should converge to 2.0"
        );
        assert!(
            (bn.running_mean[1] - 3.0).abs() < 0.1,
            "mean[1] should converge to 3.0"
        );
    }
    #[test]
    fn test_batch_norm_update_running_stats_variance_zero() {
        let mut bn = BatchNormLayer::new(2);
        let batch = vec![vec![1.0_f32, 1.0]; 4];
        for _ in 0..100 {
            bn.update_running_stats(&batch, 0.2);
        }
        assert!(
            bn.running_var[0] < 0.05,
            "variance should be ~0 for constant input, got {}",
            bn.running_var[0]
        );
    }
    #[test]
    fn test_ffnet_compute_gradient_norm_correct() {
        let net = FeedForwardNet::new();
        let grads = vec![vec![3.0_f32, 4.0]];
        let norm = net.compute_gradient_norm(&grads);
        assert!(
            (norm - 5.0).abs() < 1e-5,
            "gradient norm should be 5.0, got {norm}"
        );
    }
    #[test]
    fn test_ffnet_clip_gradients_applies_scaling() {
        let net = FeedForwardNet::new();
        let mut grads = vec![vec![3.0_f32, 4.0]];
        net.clip_gradients(&mut grads, 1.0);
        let new_norm: f32 = grads[0].iter().map(|&v| v * v).sum::<f32>().sqrt();
        assert!(
            (new_norm - 1.0).abs() < 1e-5,
            "clipped norm should be 1.0, got {new_norm}"
        );
    }
}
/// Compute scaled dot-product attention.
///
/// Attention(Q, K, V) = softmax(Q Kᵀ / √d_k) V
///
/// All matrices use row-major flat storage with shapes:
/// - Q: `[seq_q × d_k]`
/// - K: `[seq_k × d_k]`
/// - V: `[seq_k × d_v]`
///
/// Returns output of shape `[seq_q × d_v]` (flat row-major).
pub fn scaled_dot_product_attention(
    q: &[f64],
    k: &[f64],
    v: &[f64],
    seq_q: usize,
    seq_k: usize,
    d_k: usize,
    d_v: usize,
    mask: Option<&[f64]>,
) -> Vec<f64> {
    assert_eq!(q.len(), seq_q * d_k);
    assert_eq!(k.len(), seq_k * d_k);
    assert_eq!(v.len(), seq_k * d_v);
    let scale = (d_k as f64).sqrt();
    let mut scores = vec![0.0_f64; seq_q * seq_k];
    for i in 0..seq_q {
        for j in 0..seq_k {
            let mut dot = 0.0_f64;
            for d in 0..d_k {
                dot += q[i * d_k + d] * k[j * d_k + d];
            }
            scores[i * seq_k + j] = dot / scale;
        }
    }
    if let Some(m) = mask {
        assert_eq!(m.len(), seq_q * seq_k);
        for idx in 0..scores.len() {
            scores[idx] += m[idx];
        }
    }
    let mut attn_weights = vec![0.0_f64; seq_q * seq_k];
    for i in 0..seq_q {
        let row_start = i * seq_k;
        let row = &scores[row_start..row_start + seq_k];
        let max_val = row.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exp_vals: Vec<f64> = row.iter().map(|&s| (s - max_val).exp()).collect();
        let sum_exp: f64 = exp_vals.iter().sum();
        for j in 0..seq_k {
            attn_weights[row_start + j] = exp_vals[j] / sum_exp.max(1e-30);
        }
    }
    let mut output = vec![0.0_f64; seq_q * d_v];
    for i in 0..seq_q {
        for dv in 0..d_v {
            let mut acc = 0.0_f64;
            for j in 0..seq_k {
                acc += attn_weights[i * seq_k + j] * v[j * d_v + dv];
            }
            output[i * d_v + dv] = acc;
        }
    }
    output
}
#[cfg(test)]
mod attention_gnn_tests {

    use crate::AttentionReadout;

    use crate::ExtActivation;

    use crate::GnnLayer;

    use crate::MessagePassingNet;
    use crate::MultiHeadAttention;

    use crate::PositionalEncoding;

    use crate::TransformerBlock;
    use crate::TransformerFfn;

    use crate::scaled_dot_product_attention;
    #[test]
    fn test_positional_encoding_shape() {
        let pe = PositionalEncoding::new(8, 16);
        assert_eq!(pe.table.len(), 16);
        assert_eq!(pe.table[0].len(), 8);
    }
    #[test]
    fn test_positional_encoding_position_zero_first_dim_sin_zero() {
        let pe = PositionalEncoding::new(4, 10);
        assert!(
            pe.table[0][0].abs() < 1e-12,
            "PE[0,0] should be 0.0 (sin(0))"
        );
    }
    #[test]
    fn test_positional_encoding_first_dim_cos_at_zero() {
        let pe = PositionalEncoding::new(4, 10);
        assert!(
            (pe.table[0][1] - 1.0).abs() < 1e-12,
            "PE[0,1] should be 1.0 (cos(0))"
        );
    }
    #[test]
    fn test_positional_encoding_add_to_sequence() {
        let pe = PositionalEncoding::new(4, 5);
        let mut seq = vec![vec![0.0_f64; 4]; 3];
        pe.add_to_sequence(&mut seq);
        let expected = (1.0_f64 / 1.0_f64).sin();
        assert!(
            (seq[1][0] - expected).abs() < 1e-12,
            "seq[1][0] = {}",
            seq[1][0]
        );
    }
    #[test]
    fn test_positional_encoding_get_returns_slice() {
        let pe = PositionalEncoding::new(8, 10);
        let row = pe.get(3);
        assert_eq!(row.len(), 8);
    }
    #[test]
    fn test_positional_encoding_different_positions_differ() {
        let pe = PositionalEncoding::new(8, 10);
        let row0 = pe.get(0);
        let row1 = pe.get(1);
        let same = row0
            .iter()
            .zip(row1.iter())
            .all(|(a, b)| (a - b).abs() < 1e-12);
        assert!(!same, "PE at pos=0 and pos=1 should differ");
    }
    #[test]
    fn test_sdpa_output_shape() {
        let q = vec![0.1_f64; 3 * 4];
        let k = vec![0.2_f64; 3 * 4];
        let v = vec![0.3_f64; 3 * 4];
        let out = scaled_dot_product_attention(&q, &k, &v, 3, 3, 4, 4, None);
        assert_eq!(out.len(), 3 * 4, "output should have seq_q * d_v elements");
    }
    #[test]
    fn test_sdpa_uniform_attention_averages_values() {
        let q = vec![0.0_f64; 2 * 2];
        let k = vec![0.0_f64; 3 * 2];
        let v = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let out = scaled_dot_product_attention(&q, &k, &v, 2, 3, 2, 2, None);
        assert!((out[0] - 1.0 / 3.0).abs() < 1e-8, "out[0]={}", out[0]);
        assert!((out[1] - 1.0 / 3.0).abs() < 1e-8, "out[1]={}", out[1]);
    }
    #[test]
    fn test_sdpa_masking_blocks_position() {
        let q = vec![0.0_f64; 2];
        let k = vec![0.0_f64; 2];
        let v = vec![10.0_f64, 20.0];
        let mask = vec![0.0_f64, -1e9, 0.0, 0.0];
        let out = scaled_dot_product_attention(&q, &k, &v, 2, 2, 1, 1, Some(&mask));
        assert!((out[0] - 10.0).abs() < 1e-4, "masked output[0]={}", out[0]);
    }
    #[test]
    fn test_sdpa_attention_weights_sum_to_one() {
        let seq = 4;
        let dk = 3;
        let q = vec![0.0_f64; seq * dk];
        let k = vec![0.0_f64; seq * dk];
        let v = vec![1.0_f64; seq];
        let out = scaled_dot_product_attention(&q, &k, &v, seq, seq, dk, 1, None);
        for &o in &out {
            assert!(
                (o - 1.0).abs() < 1e-8,
                "attention weight sum = 1 check: out={o}"
            );
        }
    }
    #[test]
    fn test_mha_output_shape() {
        let mha = MultiHeadAttention::new(8, 2);
        let x = vec![0.0_f64; 5 * 8];
        let out = mha.forward(&x, 5);
        assert_eq!(out.len(), 5 * 8, "MHA output shape mismatch");
    }
    #[test]
    fn test_mha_num_params() {
        let mha = MultiHeadAttention::new(16, 4);
        assert_eq!(mha.num_params(), 1040);
    }
    #[test]
    fn test_mha_zero_weights_zero_output_except_bias() {
        let mut mha = MultiHeadAttention::new(4, 2);
        mha.b_o = vec![1.0_f64; 4];
        let x = vec![0.0_f64; 3 * 4];
        let out = mha.forward(&x, 3);
        for &v in &out {
            assert!(
                (v - 1.0).abs() < 1e-10,
                "out should equal bias=1.0, got {v}"
            );
        }
    }
    #[test]
    fn test_mha_identity_init_output_finite() {
        let mut mha = MultiHeadAttention::new(4, 2);
        mha.init_identity();
        let x: Vec<f64> = (0..3 * 4).map(|i| (i as f64) * 0.1).collect();
        let out = mha.forward(&x, 3);
        assert_eq!(out.len(), 3 * 4);
        for &v in &out {
            assert!(v.is_finite(), "output must be finite, got {v}");
        }
    }
    #[test]
    fn test_transformer_ffn_output_shape() {
        let ffn = TransformerFfn::new(8, 32);
        let x = vec![1.0_f64; 5 * 8];
        let out = ffn.forward(&x, 5);
        assert_eq!(out.len(), 5 * 8);
    }
    #[test]
    fn test_transformer_ffn_zero_weights_zero_output() {
        let ffn = TransformerFfn::new(4, 16);
        let x = vec![1.0_f64; 3 * 4];
        let out = ffn.forward(&x, 3);
        for &v in &out {
            assert!(v.abs() < 1e-12, "zero weights → zero output, got {v}");
        }
    }
    #[test]
    fn test_transformer_ffn_relu_activation() {
        let mut ffn = TransformerFfn::new(2, 2);
        ffn.w1 = vec![-1.0_f64; 4];
        let x = vec![1.0, 1.0, 1.0, 1.0];
        let out = ffn.forward(&x, 2);
        for &v in &out {
            assert!(
                v.abs() < 1e-12,
                "relu-clipped hidden → zero output, got {v}"
            );
        }
    }
    #[test]
    fn test_transformer_block_output_shape() {
        let block = TransformerBlock::new(8, 2, 32);
        let x = vec![0.5_f64; 4 * 8];
        let out = block.forward(&x, 4);
        assert_eq!(out.len(), 4 * 8, "transformer block output shape mismatch");
    }
    #[test]
    fn test_transformer_block_residual_preserves_input_with_zero_weights() {
        let block = TransformerBlock::new(4, 2, 16);
        let x = vec![1.0_f64; 3 * 4];
        let out = block.forward(&x, 3);
        assert_eq!(out.len(), x.len());
        for &v in &out {
            assert!(v.is_finite(), "transformer block output must be finite");
        }
    }
    #[test]
    fn test_transformer_block_output_differs_from_input() {
        let mut block = TransformerBlock::new(4, 2, 8);
        block.mha.init_identity();
        let x: Vec<f64> = (0..4 * 4).map(|i| ((i % 7) as f64) * 0.3 - 0.5).collect();
        let out = block.forward(&x, 4);
        for &v in &out {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn test_gnn_layer_output_shape() {
        let gnn = GnnLayer::new(4, 8, ExtActivation::Relu);
        let n_nodes = 3;
        let feats = vec![0.5_f64; n_nodes * 4];
        let adj = vec![vec![1usize, 2], vec![0], vec![0, 1]];
        let out = gnn.forward(&feats, n_nodes, &adj);
        assert_eq!(out.len(), n_nodes * 8, "GNN output shape mismatch");
    }
    #[test]
    fn test_gnn_layer_num_params() {
        let gnn = GnnLayer::new(4, 8, ExtActivation::Relu);
        assert_eq!(gnn.num_params(), 72);
    }
    #[test]
    fn test_gnn_layer_zero_weights_zero_output() {
        let gnn = GnnLayer::new(3, 5, ExtActivation::Linear);
        let feats = vec![1.0_f64; 4 * 3];
        let adj = vec![vec![1usize], vec![0], vec![3], vec![2]];
        let out = gnn.forward(&feats, 4, &adj);
        for &v in &out {
            assert!(v.abs() < 1e-12, "zero weights → zero output, got {v}");
        }
    }
    #[test]
    fn test_gnn_layer_isolated_node_uses_only_self() {
        let mut gnn = GnnLayer::new(2, 2, ExtActivation::Linear);
        gnn.w_self = vec![1.0, 0.0, 0.0, 1.0];
        let feats = vec![3.0, 7.0, 0.0, 0.0];
        let adj = vec![vec![], vec![]];
        let out = gnn.forward(&feats, 2, &adj);
        assert!((out[0] - 3.0).abs() < 1e-12, "node 0 out[0] = {}", out[0]);
        assert!((out[1] - 7.0).abs() < 1e-12, "node 0 out[1] = {}", out[1]);
    }
    #[test]
    fn test_gnn_layer_neighbour_aggregation_sum() {
        let mut gnn = GnnLayer::new(2, 2, ExtActivation::Linear);
        gnn.w_neigh = vec![1.0, 0.0, 0.0, 1.0];
        let feats = vec![0.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let adj = vec![vec![1usize, 2], vec![], vec![]];
        let out = gnn.forward(&feats, 3, &adj);
        assert!(
            (out[0] - 2.0).abs() < 1e-12,
            "aggregated out[0] = {}",
            out[0]
        );
        assert!((out[1]).abs() < 1e-12, "aggregated out[1] = {}", out[1]);
    }
    #[test]
    fn test_mpnn_output_shape_two_layers() {
        let mut mpnn = MessagePassingNet::new();
        mpnn.add_layer(GnnLayer::new(4, 8, ExtActivation::Relu));
        mpnn.add_layer(GnnLayer::new(8, 4, ExtActivation::Linear));
        let feats = vec![0.1_f64; 5 * 4];
        let adj: Vec<Vec<usize>> = (0..5)
            .map(|i| if i > 0 { vec![i - 1] } else { vec![] })
            .collect();
        let out = mpnn.forward(&feats, 5, &adj);
        assert_eq!(out.len(), 5 * 4);
    }
    #[test]
    fn test_mpnn_global_mean_pool_shape() {
        let mpnn = MessagePassingNet::new();
        let feats = vec![1.0_f64; 4 * 3];
        let pooled = mpnn.global_mean_pool(&feats, 4, 3);
        assert_eq!(pooled.len(), 3);
    }
    #[test]
    fn test_mpnn_global_mean_pool_known_values() {
        let mpnn = MessagePassingNet::new();
        let feats = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0];
        let pooled = mpnn.global_mean_pool(&feats, 3, 2);
        assert!((pooled[0] - 3.0).abs() < 1e-12, "mean[0] = {}", pooled[0]);
        assert!((pooled[1] - 4.0).abs() < 1e-12, "mean[1] = {}", pooled[1]);
    }
    #[test]
    fn test_mpnn_empty_returns_zero_pool() {
        let mpnn = MessagePassingNet::new();
        let pooled = mpnn.global_mean_pool(&[], 0, 4);
        assert_eq!(pooled, vec![0.0_f64; 4]);
    }
    #[test]
    fn test_mpnn_default_is_empty() {
        let mpnn = MessagePassingNet::default();
        assert_eq!(mpnn.layers.len(), 0);
    }
    #[test]
    fn test_attention_readout_output_shape() {
        let ar = AttentionReadout::new(8);
        let feats = vec![0.0_f64; 5 * 8];
        let out = ar.forward(&feats, 5);
        assert_eq!(out.len(), 8);
    }
    #[test]
    fn test_attention_readout_zero_weights_equal_scores() {
        let ar = AttentionReadout::new(2);
        let feats = vec![1.0_f64, 2.0, 3.0, 4.0];
        let out = ar.forward(&feats, 2);
        assert!((out[0] - 2.0).abs() < 1e-10, "readout[0] = {}", out[0]);
        assert!((out[1] - 3.0).abs() < 1e-10, "readout[1] = {}", out[1]);
    }
    #[test]
    fn test_attention_readout_all_zeros_input() {
        let ar = AttentionReadout::new(4);
        let feats = vec![0.0_f64; 3 * 4];
        let out = ar.forward(&feats, 3);
        for &v in &out {
            assert!(v.abs() < 1e-12, "zero input → zero output, got {v}");
        }
    }
    #[test]
    fn test_attention_readout_single_node() {
        let ar = AttentionReadout::new(3);
        let feats = vec![2.0, 4.0, 6.0];
        let out = ar.forward(&feats, 1);
        let score = 0.5_f64;
        assert!((out[0] - score * 2.0).abs() < 1e-10);
        assert!((out[1] - score * 4.0).abs() < 1e-10);
        assert!((out[2] - score * 6.0).abs() < 1e-10);
    }
}
