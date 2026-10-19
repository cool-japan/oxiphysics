// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the neural-network compute module: activations, forward
//! passes, exact backpropagation (finite-difference gradient checks for every
//! layer type), Conv1D geometry and SGD / Adam training.

use super::*;

// ── Activation function tests ────────────────────────────────────────

#[test]
fn test_relu_positive() {
    assert!((relu(3.0) - 3.0).abs() < 1e-12);
}

#[test]
fn test_relu_negative() {
    assert!((relu(-5.0)).abs() < 1e-12);
}

#[test]
fn test_relu_zero() {
    assert!((relu(0.0)).abs() < 1e-12);
}

#[test]
fn test_sigmoid_zero() {
    assert!((sigmoid(0.0) - 0.5).abs() < 1e-12);
}

#[test]
fn test_sigmoid_large_positive() {
    assert!((sigmoid(100.0) - 1.0).abs() < 1e-6);
}

#[test]
fn test_sigmoid_large_negative() {
    assert!(sigmoid(-100.0) < 1e-6);
}

#[test]
fn test_sigmoid_symmetry() {
    let x = 2.5;
    assert!((sigmoid(x) + sigmoid(-x) - 1.0).abs() < 1e-12);
}

#[test]
fn test_softmax_sums_to_one() {
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let s = softmax(&x);
    let sum: f64 = s.iter().sum();
    assert!((sum - 1.0).abs() < 1e-12);
}

#[test]
fn test_softmax_monotone() {
    let x = vec![1.0, 2.0, 3.0];
    let s = softmax(&x);
    assert!(s[0] < s[1] && s[1] < s[2]);
}

#[test]
fn test_softmax_uniform() {
    let x = vec![0.0, 0.0, 0.0];
    let s = softmax(&x);
    for &v in &s {
        assert!((v - 1.0 / 3.0).abs() < 1e-12);
    }
}

#[test]
fn test_softmax_empty() {
    let s = softmax(&[]);
    assert!(s.is_empty());
}

#[test]
fn test_softmax_single() {
    let s = softmax(&[42.0]);
    assert!((s[0] - 1.0).abs() < 1e-12);
}

#[test]
fn test_softmax_numerical_stability() {
    // Large values shouldn't overflow
    let x = vec![1000.0, 1001.0, 1002.0];
    let s = softmax(&x);
    let sum: f64 = s.iter().sum();
    assert!((sum - 1.0).abs() < 1e-10);
}

// ── MSE loss tests ───────────────────────────────────────────────────

#[test]
fn test_mse_loss_perfect() {
    let pred = vec![1.0, 2.0, 3.0];
    assert!((mse_loss(&pred, &pred)).abs() < 1e-12);
}

#[test]
fn test_mse_loss_known() {
    // mse([0, 0], [1, 1]) = 1.0
    let pred = vec![0.0, 0.0];
    let target = vec![1.0, 1.0];
    assert!((mse_loss(&pred, &target) - 1.0).abs() < 1e-12);
}

#[test]
fn test_mse_loss_empty() {
    assert!((mse_loss(&[], &[])).abs() < 1e-12);
}

#[test]
fn test_mse_loss_positive() {
    let pred = vec![1.0, 2.0, 3.0];
    let target = vec![0.0, 0.0, 0.0];
    assert!(mse_loss(&pred, &target) > 0.0);
}

// ── LayerType / NeuralLayer forward pass tests ───────────────────────

#[test]
fn test_relu_layer_forward() {
    let layer = NeuralLayer::new(3, 3, LayerType::ReLU);
    let out = layer.forward(&[-1.0, 0.0, 2.0]);
    assert_eq!(out, vec![0.0, 0.0, 2.0]);
}

#[test]
fn test_sigmoid_layer_forward() {
    let layer = NeuralLayer::new(1, 1, LayerType::Sigmoid);
    let out = layer.forward(&[0.0]);
    assert!((out[0] - 0.5).abs() < 1e-12);
}

#[test]
fn test_tanh_layer_forward() {
    let layer = NeuralLayer::new(1, 1, LayerType::Tanh);
    let out = layer.forward(&[0.0]);
    assert!((out[0]).abs() < 1e-12);
}

#[test]
fn test_softmax_layer_forward() {
    let layer = NeuralLayer::new(3, 3, LayerType::Softmax);
    let out = layer.forward(&[1.0, 2.0, 3.0]);
    let sum: f64 = out.iter().sum();
    assert!((sum - 1.0).abs() < 1e-12);
}

#[test]
fn test_dropout_layer_passthrough() {
    let layer = NeuralLayer::new(4, 4, LayerType::Dropout);
    let input = vec![1.0, 2.0, 3.0, 4.0];
    let out = layer.forward(&input);
    assert_eq!(out, input);
}

#[test]
fn test_dense_layer_identity() {
    // Single neuron with weight=1, bias=0 → identity
    let mut layer = NeuralLayer::new(1, 1, LayerType::Dense);
    layer.weights[0] = 1.0;
    let out = layer.forward(&[5.0]);
    assert!((out[0] - 5.0).abs() < 1e-12);
}

#[test]
fn test_dense_layer_known_output() {
    // 2-input → 1-output: w=[1,2], b=0.5
    let mut layer = NeuralLayer::new(2, 1, LayerType::Dense);
    layer.weights = vec![1.0, 2.0];
    layer.biases = vec![0.5];
    // out = 1*3 + 2*4 + 0.5 = 11.5
    let out = layer.forward(&[3.0, 4.0]);
    assert!((out[0] - 11.5).abs() < 1e-12);
}

#[test]
fn test_dense_layer_multi_out() {
    let mut layer = NeuralLayer::new(2, 2, LayerType::Dense);
    // Row 0: w=[1,0], b=0 → out0 = x0
    // Row 1: w=[0,1], b=0 → out1 = x1
    layer.weights = vec![1.0, 0.0, 0.0, 1.0];
    layer.biases = vec![0.0, 0.0];
    let out = layer.forward(&[7.0, 3.0]);
    assert!((out[0] - 7.0).abs() < 1e-12);
    assert!((out[1] - 3.0).abs() < 1e-12);
}

#[test]
fn test_batchnorm_layer_zero_mean() {
    let mut layer = NeuralLayer::new(4, 4, LayerType::BatchNorm);
    layer.weights = vec![1.0; 4]; // gamma = 1
    layer.biases = vec![0.0; 4]; // beta = 0
    let input = vec![1.0, 2.0, 3.0, 4.0];
    let out = layer.forward(&input);
    let mean_out: f64 = out.iter().sum::<f64>() / out.len() as f64;
    assert!(mean_out.abs() < 1e-10);
}

// ── GpuNeuralNet tests ───────────────────────────────────────────────

#[test]
fn test_empty_net_passthrough() {
    let net = GpuNeuralNet::new();
    let input = vec![1.0, 2.0, 3.0];
    let out = net.forward_pass(&input);
    assert_eq!(out, input);
}

#[test]
fn test_single_relu_net() {
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::new(3, 3, LayerType::ReLU));
    let out = net.forward_pass(&[-1.0, 0.0, 2.0]);
    assert_eq!(out, vec![0.0, 0.0, 2.0]);
}

#[test]
fn test_net_dense_then_relu() {
    let mut net = GpuNeuralNet::new();
    let mut dense = NeuralLayer::new(2, 2, LayerType::Dense);
    dense.weights = vec![1.0, 0.0, 0.0, -1.0];
    dense.biases = vec![0.0, 0.0];
    net.add_layer(dense);
    net.add_layer(NeuralLayer::new(2, 2, LayerType::ReLU));
    let out = net.forward_pass(&[3.0, 4.0]);
    // dense: [3, -4], relu: [3, 0]
    assert!((out[0] - 3.0).abs() < 1e-12);
    assert!((out[1]).abs() < 1e-12);
}

#[test]
fn test_batch_forward() {
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::new(2, 2, LayerType::ReLU));
    let inputs = vec![vec![-1.0, 2.0], vec![3.0, -4.0]];
    let outs = net.batch_forward(&inputs);
    assert_eq!(outs.len(), 2);
    assert_eq!(outs[0], vec![0.0, 2.0]);
    assert_eq!(outs[1], vec![3.0, 0.0]);
}

#[test]
fn test_net_default() {
    let net = GpuNeuralNet::default();
    assert!(net.layers.is_empty());
}

// ── BackpropGpu tests ────────────────────────────────────────────────

#[test]
fn test_backprop_gradient_shape() {
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::new(3, 2, LayerType::Dense));
    let bp = BackpropGpu::new(&net);
    assert_eq!(bp.gradients.len(), 1);
    assert_eq!(bp.gradients[0].len(), 6); // 3*2
    assert_eq!(bp.bias_gradients.len(), 1);
    assert_eq!(bp.bias_gradients[0].len(), 2);
}

/// The gradients of a dense layer are the outer product `δ xᵀ` (weights),
/// `δ` (biases) and `Wᵀ δ` (input) - not a constant broadcast of
/// `loss_grad[0]`.
#[test]
fn test_backprop_dense_gradients_are_exact() {
    let mut net = GpuNeuralNet::new();
    let mut dense = NeuralLayer::new(2, 2, LayerType::Dense);
    dense.weights = vec![1.0, 2.0, 3.0, 4.0];
    dense.biases = vec![0.5, -0.5];
    net.add_layer(dense);
    let cache = net.forward_with_cache(&[1.0, -1.0]);
    let mut bp = BackpropGpu::new(&net);
    let grad_input = bp
        .backward_pass(&net, &cache, &[0.1, 0.2])
        .expect("matching shapes");
    let expected_w = [0.1, -0.1, 0.2, -0.2];
    for (g, e) in bp.gradients[0].iter().zip(expected_w) {
        assert!((g - e).abs() < 1e-12, "dW {g} vs {e}");
    }
    assert!((bp.bias_gradients[0][0] - 0.1).abs() < 1e-12);
    assert!((bp.bias_gradients[0][1] - 0.2).abs() < 1e-12);
    // Wᵀ δ = [1·0.1 + 3·0.2, 2·0.1 + 4·0.2]
    assert!((grad_input[0] - 0.7).abs() < 1e-12);
    assert!((grad_input[1] - 1.0).abs() < 1e-12);
}

#[test]
fn test_backprop_zero_loss_grad_gives_zero_gradients() {
    let mut net = GpuNeuralNet::new();
    let mut dense = NeuralLayer::new(2, 2, LayerType::Dense);
    dense.weights = vec![0.3, -0.7, 1.1, 0.4];
    net.add_layer(dense);
    net.add_layer(NeuralLayer::new(2, 2, LayerType::Tanh));
    let cache = net.forward_with_cache(&[0.8, -1.3]);
    let mut bp = BackpropGpu::new(&net);
    let grad_input = bp
        .backward_pass(&net, &cache, &[0.0, 0.0])
        .expect("matching shapes");
    assert!(bp.gradients.iter().flatten().all(|&g| g == 0.0));
    assert!(bp.bias_gradients.iter().flatten().all(|&g| g == 0.0));
    assert!(grad_input.iter().all(|&g| g == 0.0));
}

/// Distinct output gradients must produce distinct weight gradients; the
/// old implementation filled every entry with `loss_grad[0]`.
#[test]
fn test_backprop_gradients_depend_on_every_output_gradient_entry() {
    let mut net = GpuNeuralNet::new();
    let mut dense = NeuralLayer::new(3, 2, LayerType::Dense);
    dense.weights = vec![0.2, -0.4, 0.6, 0.1, 0.3, -0.5];
    net.add_layer(dense);
    let cache = net.forward_with_cache(&[1.0, 2.0, 3.0]);
    let mut bp = BackpropGpu::new(&net);
    bp.backward_pass(&net, &cache, &[1.0, 0.0])
        .expect("matching shapes");
    // Row 1 receives no gradient, row 0 gets the input vector.
    assert_eq!(bp.gradients[0], vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);
    assert_eq!(bp.bias_gradients[0], vec![1.0, 0.0]);
}

#[test]
fn test_backprop_accumulates_over_a_minibatch_and_zero_grad_resets() {
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::new(2, 1, LayerType::Dense));
    let samples = [[1.0, 2.0], [3.0, -1.0]];
    let mut single = Vec::new();
    for s in &samples {
        let cache = net.forward_with_cache(s);
        let mut bp = BackpropGpu::new(&net);
        bp.backward_pass(&net, &cache, &[0.5]).expect("shapes");
        single.push(bp.gradients[0].clone());
    }
    let mut batch = BackpropGpu::new(&net);
    for s in &samples {
        let cache = net.forward_with_cache(s);
        batch
            .accumulate_backward(&net, &cache, &[0.5])
            .expect("shapes");
    }
    let summed: Vec<f64> = single[0]
        .iter()
        .zip(single[1].iter())
        .map(|(a, b)| a + b)
        .collect();
    assert_eq!(batch.gradients[0].len(), summed.len());
    for (accumulated, expected) in batch.gradients[0].iter().zip(summed.iter()) {
        assert!((accumulated - expected).abs() < 1e-12);
    }
    batch.zero_grad();
    assert!(batch.gradients[0].iter().all(|&g| g == 0.0));
    assert!(batch.bias_gradients[0].iter().all(|&g| g == 0.0));
}

#[test]
fn test_backprop_rejects_mismatched_shapes() {
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::new(2, 2, LayerType::Dense));
    let cache = net.forward_with_cache(&[1.0, 2.0]);
    let mut bp = BackpropGpu::new(&net);
    // Wrong loss-gradient length.
    let err = bp.backward_pass(&net, &cache, &[1.0]).expect_err("length");
    assert!(
        matches!(err, NeuralNetError::ShapeMismatch { .. }),
        "{err:?}"
    );
    // Backprop state built for a different network.
    let mut other = GpuNeuralNet::new();
    other.add_layer(NeuralLayer::new(2, 2, LayerType::Dense));
    other.add_layer(NeuralLayer::new(2, 2, LayerType::ReLU));
    let mut other_bp = BackpropGpu::new(&other);
    assert!(other_bp.backward_pass(&net, &cache, &[1.0, 1.0]).is_err());
    // Cache recorded for a different network.
    let other_cache = other.forward_with_cache(&[1.0, 2.0]);
    assert!(bp.backward_pass(&net, &other_cache, &[1.0, 1.0]).is_err());
    // A layer-level mismatch (gradient buffer of the wrong size).
    let layer = &net.layers[0];
    let mut short = vec![0.0; 3];
    let mut bias = vec![0.0; 2];
    assert!(
        layer
            .backward(
                &[1.0, 2.0],
                cache.output(),
                &[1.0, 1.0],
                &mut short,
                &mut bias
            )
            .is_err()
    );
}

#[test]
fn test_forward_cache_matches_forward_pass() {
    let (net, input, _) = every_layer_net(11);
    let cache = net.forward_with_cache(&input);
    assert_eq!(cache.output(), net.forward_pass(&input).as_slice());
    assert_eq!(cache.input(), input.as_slice());
    assert_eq!(cache.layer_outputs().len(), net.layers.len());
    assert_eq!(cache.layer_input(0), Some(input.as_slice()));
    for l in 1..net.layers.len() {
        assert_eq!(
            cache.layer_input(l),
            Some(cache.layer_outputs()[l - 1].as_slice())
        );
    }
    assert!(cache.layer_input(net.layers.len() + 1).is_none());
    // A network without layers passes its input through.
    let empty = GpuNeuralNet::new().forward_with_cache(&[1.0, 2.0]);
    assert_eq!(empty.output(), &[1.0, 2.0]);
}

// ── Finite-difference gradient checks ────────────────────────────────

/// Deterministic generator of values in `[-1, 1)`; independent of any
/// external RNG so the tests never change with a dependency upgrade.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64) * 2.0 - 1.0
    }
}

/// Fill a layer's weights and biases with values in `[-scale, scale)`.
fn randomize(layer: &mut NeuralLayer, rng: &mut Lcg, scale: f64) {
    for w in layer.weights.iter_mut() {
        *w = scale * rng.next_unit();
    }
    for b in layer.biases.iter_mut() {
        *b = 0.5 * scale * rng.next_unit();
    }
}

/// A 9-layer network containing every layer type:
/// `Conv1D(8 -> 2x3) ReLU BatchNorm Tanh Dense(6 -> 4) Sigmoid Dropout
/// Dense(4 -> 3) Softmax`, plus an input and a target.
fn every_layer_net(seed: u64) -> (GpuNeuralNet, Vec<f64>, Vec<f64>) {
    let mut rng = Lcg(seed);
    let mut net = GpuNeuralNet::new();
    let mut conv = NeuralLayer::conv1d(1, 2, 8, 3, 2).expect("valid geometry");
    randomize(&mut conv, &mut rng, 0.9);
    net.add_layer(conv);
    net.add_layer(NeuralLayer::new(6, 6, LayerType::ReLU));
    let mut norm = NeuralLayer::new(6, 6, LayerType::BatchNorm);
    for i in 0..6 {
        norm.weights[i] = 1.0 + 0.3 * rng.next_unit();
        norm.biases[i] = 0.2 * rng.next_unit();
    }
    net.add_layer(norm);
    net.add_layer(NeuralLayer::new(6, 6, LayerType::Tanh));
    let mut dense1 = NeuralLayer::new(6, 4, LayerType::Dense);
    randomize(&mut dense1, &mut rng, 0.8);
    net.add_layer(dense1);
    net.add_layer(NeuralLayer::new(4, 4, LayerType::Sigmoid));
    net.add_layer(NeuralLayer::new(4, 4, LayerType::Dropout));
    let mut dense2 = NeuralLayer::new(4, 3, LayerType::Dense);
    randomize(&mut dense2, &mut rng, 0.8);
    net.add_layer(dense2);
    net.add_layer(NeuralLayer::new(3, 3, LayerType::Softmax));
    let input: Vec<f64> = (0..8).map(|_| rng.next_unit()).collect();
    let target = vec![0.2, 0.5, 0.3];
    (net, input, target)
}

/// MSE loss of one sample.
fn sample_loss(net: &GpuNeuralNet, input: &[f64], target: &[f64]) -> f64 {
    mse_loss(&net.forward_pass(input), target)
}

/// Analytic gradients (weights, biases, input) of the MSE loss.
fn analytic_gradients(
    net: &GpuNeuralNet,
    input: &[f64],
    target: &[f64],
) -> (BackpropGpu, Vec<f64>) {
    let cache = net.forward_with_cache(input);
    let n = cache.output().len().min(target.len()) as f64;
    let loss_grad: Vec<f64> = cache
        .output()
        .iter()
        .zip(target.iter())
        .map(|(p, t)| 2.0 * (p - t) / n)
        .collect();
    let mut bp = BackpropGpu::new(net);
    let grad_input = bp
        .backward_pass(net, &cache, &loss_grad)
        .expect("backward pass");
    (bp, grad_input)
}

/// Largest deviation between the analytic gradient and a central finite
/// difference over every weight, bias and input entry.
fn max_gradient_error(net: &GpuNeuralNet, input: &[f64], target: &[f64]) -> f64 {
    let h = 1e-6;
    let (bp, grad_input) = analytic_gradients(net, input, target);
    let mut worst = 0.0_f64;
    let mut probe = net.clone();
    for l in 0..net.layers.len() {
        for i in 0..net.layers[l].weights.len() {
            let original = probe.layers[l].weights[i];
            probe.layers[l].weights[i] = original + h;
            let plus = sample_loss(&probe, input, target);
            probe.layers[l].weights[i] = original - h;
            let minus = sample_loss(&probe, input, target);
            probe.layers[l].weights[i] = original;
            let numeric = (plus - minus) / (2.0 * h);
            worst = worst.max((numeric - bp.gradients[l][i]).abs());
        }
        for i in 0..net.layers[l].biases.len() {
            let original = probe.layers[l].biases[i];
            probe.layers[l].biases[i] = original + h;
            let plus = sample_loss(&probe, input, target);
            probe.layers[l].biases[i] = original - h;
            let minus = sample_loss(&probe, input, target);
            probe.layers[l].biases[i] = original;
            let numeric = (plus - minus) / (2.0 * h);
            worst = worst.max((numeric - bp.bias_gradients[l][i]).abs());
        }
    }
    let mut shifted = input.to_vec();
    for i in 0..input.len() {
        shifted[i] = input[i] + h;
        let plus = sample_loss(net, &shifted, target);
        shifted[i] = input[i] - h;
        let minus = sample_loss(net, &shifted, target);
        shifted[i] = input[i];
        let numeric = (plus - minus) / (2.0 * h);
        worst = worst.max((numeric - grad_input[i]).abs());
    }
    worst
}

/// Every weight, bias and input gradient of a network containing all layer
/// types matches a central finite difference to 1e-6.
#[test]
fn test_gradient_check_network_with_every_layer_type() {
    for seed in [11_u64, 29, 1234] {
        let (net, input, target) = every_layer_net(seed);
        // Keep the ReLU away from its kink and BatchNorm away from a
        // vanishing variance, where finite differences lose accuracy.
        let pre_relu = net.layers[0].forward(&input);
        assert!(
            pre_relu.iter().all(|z| z.abs() > 1e-3),
            "seed {seed}: ReLU input too close to zero: {pre_relu:?}"
        );
        let post_relu = net.layers[1].forward(&pre_relu);
        let mean = post_relu.iter().sum::<f64>() / 6.0;
        let var = post_relu.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / 6.0;
        assert!(
            var > 0.01,
            "seed {seed}: BatchNorm variance {var} too small"
        );

        let error = max_gradient_error(&net, &input, &target);
        assert!(error < 1e-6, "seed {seed}: gradient error {error:e}");
    }
}

/// The same check layer type by layer, so a failure names its culprit.
#[test]
fn test_gradient_check_each_layer_type_alone() {
    let mut rng = Lcg(77);
    let input5: Vec<f64> = vec![0.7, -0.4, 1.3, -1.1, 0.25];
    let target5: Vec<f64> = vec![0.1, -0.2, 0.3, 0.0, 0.5];
    let target3: Vec<f64> = vec![0.2, 0.3, 0.5];
    let layers: Vec<(&str, NeuralLayer, Vec<f64>)> = {
        let mut dense = NeuralLayer::new(5, 3, LayerType::Dense);
        randomize(&mut dense, &mut rng, 0.9);
        let mut conv = NeuralLayer::conv1d(1, 2, 5, 2, 1).expect("valid");
        randomize(&mut conv, &mut rng, 0.9);
        let mut norm = NeuralLayer::new(5, 5, LayerType::BatchNorm);
        for i in 0..5 {
            norm.weights[i] = 1.0 + 0.4 * rng.next_unit();
            norm.biases[i] = 0.3 * rng.next_unit();
        }
        vec![
            ("dense", dense, target3.clone()),
            (
                "conv1d",
                conv,
                vec![0.1, 0.2, 0.3, 0.4, -0.1, -0.2, -0.3, 0.0],
            ),
            (
                "relu",
                NeuralLayer::new(5, 5, LayerType::ReLU),
                target5.clone(),
            ),
            (
                "sigmoid",
                NeuralLayer::new(5, 5, LayerType::Sigmoid),
                target5.clone(),
            ),
            (
                "tanh",
                NeuralLayer::new(5, 5, LayerType::Tanh),
                target5.clone(),
            ),
            (
                "softmax",
                NeuralLayer::new(5, 5, LayerType::Softmax),
                target5.clone(),
            ),
            ("batchnorm", norm, target5.clone()),
            (
                "dropout",
                NeuralLayer::new(5, 5, LayerType::Dropout),
                target5.clone(),
            ),
        ]
    };
    for (name, layer, target) in layers {
        let mut net = GpuNeuralNet::new();
        net.add_layer(layer);
        let error = max_gradient_error(&net, &input5, &target);
        assert!(error < 1e-6, "{name}: gradient error {error:e}");
    }
}

/// Multi-channel / strided / longer convolutions (the geometry that makes
/// Conv1D differ from Dense) are differentiated correctly as well.
#[test]
fn test_gradient_check_conv1d_geometries() {
    let mut rng = Lcg(5);
    for &(in_ch, out_ch, len, kernel, stride) in &[
        (1_usize, 1_usize, 9_usize, 3_usize, 1_usize),
        (2, 3, 10, 4, 2),
        (3, 2, 7, 2, 3),
        (1, 2, 6, 6, 1),
    ] {
        let mut conv = NeuralLayer::conv1d(in_ch, out_ch, len, kernel, stride).expect("valid");
        randomize(&mut conv, &mut rng, 0.8);
        let out_size = conv.output_size;
        let mut net = GpuNeuralNet::new();
        net.add_layer(conv);
        let input: Vec<f64> = (0..in_ch * len).map(|_| rng.next_unit()).collect();
        let target: Vec<f64> = (0..out_size).map(|_| rng.next_unit()).collect();
        let error = max_gradient_error(&net, &input, &target);
        assert!(
            error < 1e-6,
            "conv1d({in_ch},{out_ch},{len},{kernel},{stride}): gradient error {error:e}"
        );
    }
}

// ── Layer-level backward values ──────────────────────────────────────

#[test]
fn test_activation_backward_known_values() {
    // ReLU: gradient passes only where the input is positive.
    let layer = NeuralLayer::new(3, 3, LayerType::ReLU);
    let x = [-1.0, 0.0, 2.0];
    let y = layer.forward(&x);
    let mut gw = vec![0.0; layer.weights.len()];
    let mut gb = vec![0.0; layer.biases.len()];
    let gx = layer
        .backward(&x, &y, &[1.0, 1.0, 1.0], &mut gw, &mut gb)
        .expect("shapes");
    assert_eq!(gx, vec![0.0, 0.0, 1.0]);

    // Sigmoid at 0: s = 1/2, ds = 1/4.
    let layer = NeuralLayer::new(1, 1, LayerType::Sigmoid);
    let y = layer.forward(&[0.0]);
    let mut gw = vec![0.0; 1];
    let mut gb = vec![0.0; 1];
    let gx = layer
        .backward(&[0.0], &y, &[2.0], &mut gw, &mut gb)
        .expect("shapes");
    assert!((gx[0] - 0.5).abs() < 1e-12);

    // Tanh at 0: derivative 1.
    let layer = NeuralLayer::new(1, 1, LayerType::Tanh);
    let y = layer.forward(&[0.0]);
    let gx = layer
        .backward(&[0.0], &y, &[3.0], &mut gw, &mut gb)
        .expect("shapes");
    assert!((gx[0] - 3.0).abs() < 1e-12);

    // Softmax of equal logits: dx_i = s_i (g_i - sum_j g_j s_j).
    let layer = NeuralLayer::new(3, 3, LayerType::Softmax);
    let x = [0.0, 0.0, 0.0];
    let y = layer.forward(&x);
    let mut gw = vec![0.0; 9];
    let mut gb = vec![0.0; 3];
    let gx = layer
        .backward(&x, &y, &[1.0, 0.0, 0.0], &mut gw, &mut gb)
        .expect("shapes");
    let expected = [2.0 / 9.0, -1.0 / 9.0, -1.0 / 9.0];
    for (g, e) in gx.iter().zip(expected) {
        assert!((g - e).abs() < 1e-12, "{g} vs {e}");
    }
    // The softmax Jacobian-vector product of any g sums to zero.
    assert!(gx.iter().sum::<f64>().abs() < 1e-12);

    // Dropout is the identity in both directions.
    let layer = NeuralLayer::new(2, 2, LayerType::Dropout);
    let mut gw = vec![0.0; 4];
    let mut gb = vec![0.0; 2];
    let gx = layer
        .backward(&[1.0, 2.0], &[1.0, 2.0], &[0.3, -0.4], &mut gw, &mut gb)
        .expect("shapes");
    assert_eq!(gx, vec![0.3, -0.4]);
}

/// BatchNorm: `dγ = δ x̂`, `dβ = δ`, and the input gradient of a
/// normalisation sums to zero (shifting every input leaves `x̂` fixed).
#[test]
fn test_batchnorm_backward_has_trainable_gamma_and_beta() {
    let mut layer = NeuralLayer::new(4, 4, LayerType::BatchNorm);
    layer.weights = vec![1.5, 0.5, 2.0, 1.0];
    layer.biases = vec![0.1, 0.2, 0.3, 0.4];
    let x = [1.0, 2.0, 4.0, 7.0];
    let y = layer.forward(&x);
    let delta = [0.5, -1.0, 0.25, 2.0];
    let mut gw = vec![0.0; 4];
    let mut gb = vec![0.0; 4];
    let gx = layer
        .backward(&x, &y, &delta, &mut gw, &mut gb)
        .expect("shapes");

    let mean = x.iter().sum::<f64>() / 4.0;
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / 4.0;
    let std = (var + 1e-5).sqrt();
    for i in 0..4 {
        let x_hat = (x[i] - mean) / std;
        assert!((gw[i] - delta[i] * x_hat).abs() < 1e-12, "dgamma[{i}]");
        assert!((gb[i] - delta[i]).abs() < 1e-12, "dbeta[{i}]");
    }
    assert!(gx.iter().sum::<f64>().abs() < 1e-12, "dx must sum to zero");
    assert!(gx.iter().any(|g| g.abs() > 1e-3));
}

// ── Conv1D ───────────────────────────────────────────────────────────

#[test]
fn test_conv1d_spec_validation() {
    assert!(Conv1dSpec::new(1, 1, 8, 3, 1).is_ok());
    for bad in [
        Conv1dSpec::new(0, 1, 8, 3, 1),
        Conv1dSpec::new(1, 0, 8, 3, 1),
        Conv1dSpec::new(1, 1, 8, 0, 1),
        Conv1dSpec::new(1, 1, 8, 3, 0),
        Conv1dSpec::new(1, 1, 2, 3, 1),
    ] {
        assert!(
            matches!(bad, Err(NeuralNetError::InvalidLayer(_))),
            "{bad:?}"
        );
    }
    assert!(NeuralLayer::conv1d(1, 1, 2, 3, 1).is_err());
}

#[test]
fn test_conv1d_geometry_follows_kernel_and_stride() {
    let layer = NeuralLayer::conv1d(2, 3, 11, 4, 3).expect("valid");
    let spec = layer.conv_spec().expect("conv layer");
    // positions = (11 - 4) / 3 + 1 = 3
    assert_eq!(spec.output_len(), 3);
    assert_eq!(layer.input_size, 22);
    assert_eq!(layer.output_size, 9);
    assert_eq!(layer.weights.len(), 3 * 2 * 4);
    assert_eq!(layer.biases.len(), 3);
    assert_eq!(layer.forward(&[0.0; 22]).len(), 9);
    assert!(
        NeuralLayer::new(3, 3, LayerType::ReLU)
            .conv_spec()
            .is_none()
    );
}

/// A sliding-window correlation, not a dense matrix product: a width-3
/// difference kernel over five samples yields three outputs.
#[test]
fn test_conv1d_layer_forward_slides_the_kernel() {
    let mut layer = NeuralLayer::conv1d(1, 1, 5, 3, 1).expect("valid");
    layer.weights = vec![1.0, 0.0, -1.0];
    layer.biases = vec![0.5];
    let out = layer.forward(&[1.0, 2.0, 4.0, 8.0, 16.0]);
    // y[p] = x[p] - x[p+2] + 0.5
    assert_eq!(out, vec![-2.5, -5.5, -11.5]);
}

#[test]
fn test_conv1d_stride_skips_positions() {
    let mut layer = NeuralLayer::conv1d(1, 1, 6, 2, 2).expect("valid");
    layer.weights = vec![1.0, 1.0];
    let out = layer.forward(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(out, vec![3.0, 7.0, 11.0]);
}

/// Translation equivariance - the defining property a dense layer lacks:
/// shifting the input by one sample shifts the output by one position.
#[test]
fn test_conv1d_is_translation_equivariant() {
    let mut layer = NeuralLayer::conv1d(1, 1, 7, 3, 1).expect("valid");
    layer.weights = vec![0.5, -1.0, 2.0];
    let a = layer.forward(&[0.0, 1.0, 3.0, 2.0, 0.0, 0.0, 0.0]);
    let b = layer.forward(&[0.0, 0.0, 1.0, 3.0, 2.0, 0.0, 0.0]);
    assert_eq!(&a[..4], &b[1..5]);
}

/// Multi-channel forward pass against an independent reference.
#[test]
fn test_conv1d_multichannel_matches_reference() {
    let (in_ch, out_ch, len, kernel, stride) = (2_usize, 2_usize, 6_usize, 3_usize, 1_usize);
    let mut rng = Lcg(3);
    let mut layer = NeuralLayer::conv1d(in_ch, out_ch, len, kernel, stride).expect("valid");
    randomize(&mut layer, &mut rng, 1.0);
    let input: Vec<f64> = (0..in_ch * len).map(|_| rng.next_unit()).collect();
    let out = layer.forward(&input);
    let positions = (len - kernel) / stride + 1;
    assert_eq!(out.len(), out_ch * positions);
    for oc in 0..out_ch {
        for p in 0..positions {
            let mut expected = layer.biases[oc];
            for ic in 0..in_ch {
                for k in 0..kernel {
                    expected += layer.weights[oc * in_ch * kernel + ic * kernel + k]
                        * input[ic * len + p + k];
                }
            }
            assert!(
                (out[oc * positions + p] - expected).abs() < 1e-12,
                "oc {oc} p {p}"
            );
        }
    }
}

/// `NeuralLayer::new(.., Conv1D)` is the full-width convolution: one dot
/// product per filter (the shape the legacy constructor always had).
#[test]
fn test_conv1d_layer_forward() {
    let mut layer = NeuralLayer::new(3, 1, LayerType::Conv1D);
    layer.weights = vec![1.0, 1.0, 1.0];
    layer.biases = vec![0.0];
    let out = layer.forward(&[1.0, 2.0, 3.0]);
    assert!((out[0] - 6.0).abs() < 1e-12);
    let spec = layer.conv_spec().expect("conv layer");
    assert_eq!(
        (spec.kernel_size, spec.stride, spec.output_len()),
        (3, 1, 1)
    );
}

/// A hand-constructed `Conv1D` layer without a geometry falls back to the
/// full-width reading instead of panicking.
#[test]
fn test_conv1d_without_spec_reads_as_full_width() {
    let layer = NeuralLayer {
        weights: vec![1.0, 2.0, 3.0, 4.0],
        biases: vec![0.0, 1.0],
        layer_type: LayerType::Conv1D,
        input_size: 2,
        output_size: 2,
        conv: None,
    };
    assert_eq!(layer.forward(&[1.0, 1.0]), vec![3.0, 8.0]);
}

// ── AdamOptimizer tests ──────────────────────────────────────────────

#[test]
fn test_adam_decreases_loss() {
    let mut params = vec![1.0, -1.0, 2.0];
    let mut adam = AdamOptimizer::new(3, 0.1, 0.9, 0.999, 1e-8);
    // Target: params = 0, grad = 2*params
    for _ in 0..500 {
        let grads: Vec<f64> = params.iter().map(|&p| 2.0 * p).collect();
        adam.update(&mut params, &grads);
    }
    for &p in &params {
        assert!(p.abs() < 0.1, "param={p}");
    }
}

#[test]
fn test_adam_timestep_increments() {
    let mut adam = AdamOptimizer::new(2, 0.01, 0.9, 0.999, 1e-8);
    let mut params = vec![1.0, 1.0];
    let grads = vec![0.1, 0.1];
    adam.update(&mut params, &grads);
    assert_eq!(adam.t, 1);
    adam.update(&mut params, &grads);
    assert_eq!(adam.t, 2);
}

#[test]
fn test_adam_moment_buffers_update() {
    let mut adam = AdamOptimizer::new(1, 0.01, 0.9, 0.999, 1e-8);
    let mut params = vec![1.0];
    adam.update(&mut params, &[0.5]);
    assert!((adam.m[0] - 0.1 * 0.5).abs() < 1e-12); // (1-0.9)*0.5
    assert!(adam.v[0] > 0.0);
}

// ── GpuTrainer tests ─────────────────────────────────────────────────

/// A `1 -> 1` dense network with the given weight and bias.
fn scalar_net(weight: f64, bias: f64) -> GpuNeuralNet {
    let mut net = GpuNeuralNet::new();
    let mut layer = NeuralLayer::new(1, 1, LayerType::Dense);
    layer.weights = vec![weight];
    layer.biases = vec![bias];
    net.add_layer(layer);
    net
}

/// One SGD step: `w -= lr dL/dw` and `b -= lr dL/db`.  For `w = b = 0`,
/// `x = 1`, `t = 1` the loss gradient is `-2`, so both move by `+2 lr`.
#[test]
fn test_trainer_sgd_step_updates_weights_and_biases() {
    let mut trainer = GpuTrainer::new(scalar_net(0.0, 0.0), 0.1, OptimizerType::Sgd);
    let loss = trainer.train_step(&[1.0], &[1.0]).expect("shapes agree");
    assert!((loss - 1.0).abs() < 1e-12);
    assert!((trainer.net.layers[0].weights[0] - 0.2).abs() < 1e-12);
    assert!((trainer.net.layers[0].biases[0] - 0.2).abs() < 1e-12);
}

/// The first Adam step moves every parameter by about `lr * sign(grad)`;
/// the bias must move too.
#[test]
fn test_trainer_adam_step_updates_weights_and_biases() {
    let mut trainer = GpuTrainer::new(scalar_net(0.0, 0.0), 0.01, OptimizerType::Adam);
    let adam_len = trainer.adam.as_ref().map(|a| a.m.len());
    assert_eq!(adam_len, Some(2), "Adam must track the weight and the bias");
    trainer.train_step(&[1.0], &[1.0]).expect("shapes agree");
    assert!((trainer.net.layers[0].weights[0] - 0.01).abs() < 1e-6);
    assert!((trainer.net.layers[0].biases[0] - 0.01).abs() < 1e-6);
}

/// Learn `y = 2x + 1`: the slope needs the weight gradient, the offset
/// needs the bias gradient.
#[test]
fn test_trainer_sgd_fits_a_line() {
    let mut trainer = GpuTrainer::new(scalar_net(0.0, 0.0), 0.05, OptimizerType::Sgd);
    let xs = [-1.0, 0.0, 1.0, 2.0];
    let mut previous = f64::INFINITY;
    for epoch in 0..400 {
        let mut epoch_loss = 0.0;
        for &x in &xs {
            epoch_loss += trainer
                .train_step(&[x], &[2.0 * x + 1.0])
                .expect("shapes agree");
        }
        if epoch < 20 {
            assert!(
                epoch_loss < previous,
                "loss must fall early on (epoch {epoch})"
            );
        }
        previous = epoch_loss;
    }
    let layer = &trainer.net.layers[0];
    assert!(
        (layer.weights[0] - 2.0).abs() < 1e-3,
        "w = {}",
        layer.weights[0]
    );
    assert!(
        (layer.biases[0] - 1.0).abs() < 1e-3,
        "b = {}",
        layer.biases[0]
    );
    assert!(previous < 1e-6);
}

#[test]
fn test_trainer_adam_fits_a_line() {
    let mut trainer = GpuTrainer::new(scalar_net(0.0, 0.0), 0.05, OptimizerType::Adam);
    let xs = [-1.0, 0.0, 1.0, 2.0];
    for _ in 0..1500 {
        for &x in &xs {
            trainer
                .train_step(&[x], &[2.0 * x + 1.0])
                .expect("shapes agree");
        }
    }
    let layer = &trainer.net.layers[0];
    assert!(
        (layer.weights[0] - 2.0).abs() < 1e-2,
        "w = {}",
        layer.weights[0]
    );
    assert!(
        (layer.biases[0] - 1.0).abs() < 1e-2,
        "b = {}",
        layer.biases[0]
    );
}

/// Hidden-layer XOR network `2 -> 4 -> 1` (tanh hidden, sigmoid output) with
/// a deterministic non-symmetric initialisation.
fn xor_net(seed: u64) -> GpuNeuralNet {
    let mut rng = Lcg(seed);
    let mut net = GpuNeuralNet::new();
    let mut hidden = NeuralLayer::new(2, 4, LayerType::Dense);
    randomize(&mut hidden, &mut rng, 1.0);
    net.add_layer(hidden);
    net.add_layer(NeuralLayer::new(4, 4, LayerType::Tanh));
    let mut output = NeuralLayer::new(4, 1, LayerType::Dense);
    randomize(&mut output, &mut rng, 1.0);
    net.add_layer(output);
    net.add_layer(NeuralLayer::new(1, 1, LayerType::Sigmoid));
    net
}

/// XOR is not linearly separable: it only converges when the gradients
/// flow correctly through the hidden layer and its nonlinearity.
#[test]
fn test_trainer_learns_xor() {
    let samples: [([f64; 2], f64); 4] = [
        ([0.0, 0.0], 0.0),
        ([0.0, 1.0], 1.0),
        ([1.0, 0.0], 1.0),
        ([1.0, 1.0], 0.0),
    ];
    let mut trainer = GpuTrainer::new(xor_net(2024), 0.05, OptimizerType::Adam);
    let mut last_epoch_loss = f64::INFINITY;
    for _ in 0..1500 {
        last_epoch_loss = 0.0;
        for (x, t) in &samples {
            last_epoch_loss += trainer.train_step(x, &[*t]).expect("shapes agree");
        }
    }
    assert!(last_epoch_loss < 0.01, "XOR loss {last_epoch_loss}");
    for (x, t) in &samples {
        let y = trainer.net.forward_pass(x)[0];
        assert!((y - t).abs() < 0.2, "XOR({x:?}) = {y}, expected about {t}");
    }
}

/// Plain SGD solves XOR as well (weights and biases both trained).
#[test]
fn test_trainer_learns_xor_with_sgd() {
    let samples: [([f64; 2], f64); 4] = [
        ([0.0, 0.0], 0.0),
        ([0.0, 1.0], 1.0),
        ([1.0, 0.0], 1.0),
        ([1.0, 1.0], 0.0),
    ];
    let mut trainer = GpuTrainer::new(xor_net(2024), 0.5, OptimizerType::Sgd);
    let mut last_epoch_loss = f64::INFINITY;
    for _ in 0..6000 {
        last_epoch_loss = 0.0;
        for (x, t) in &samples {
            last_epoch_loss += trainer.train_step(x, &[*t]).expect("shapes agree");
        }
    }
    assert!(last_epoch_loss < 0.02, "XOR loss {last_epoch_loss}");
    for (x, t) in &samples {
        let y = trainer.net.forward_pass(x)[0];
        assert!((y - t).abs() < 0.3, "XOR({x:?}) = {y}, expected about {t}");
    }
}

/// A convolution trained with SGD recovers a known filter and offset.
#[test]
fn test_trainer_learns_a_conv1d_filter() {
    let true_kernel = [1.0, -1.0, 0.5];
    let true_bias = 0.3;
    let mut rng = Lcg(99);
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::conv1d(1, 1, 8, 3, 1).expect("valid"));
    let mut trainer = GpuTrainer::new(net, 0.05, OptimizerType::Sgd);
    let inputs: Vec<Vec<f64>> = (0..6)
        .map(|_| (0..8).map(|_| rng.next_unit()).collect())
        .collect();
    let reference = |x: &[f64]| -> Vec<f64> {
        (0..6)
            .map(|p| {
                true_bias
                    + true_kernel
                        .iter()
                        .enumerate()
                        .map(|(k, w)| w * x[p + k])
                        .sum::<f64>()
            })
            .collect()
    };
    for _ in 0..800 {
        for x in &inputs {
            trainer.train_step(x, &reference(x)).expect("shapes agree");
        }
    }
    let layer = &trainer.net.layers[0];
    for (w, t) in layer.weights.iter().zip(true_kernel) {
        assert!((w - t).abs() < 1e-3, "kernel {w} vs {t}");
    }
    assert!((layer.biases[0] - true_bias).abs() < 1e-3);
}

#[test]
fn test_trainer_adam_is_sized_for_weights_and_biases() {
    let (net, _, _) = every_layer_net(11);
    let params: usize = net
        .layers
        .iter()
        .map(|l| l.weights.len() + l.biases.len())
        .sum();
    let trainer = GpuTrainer::new(net, 0.01, OptimizerType::Adam);
    let adam = trainer.adam.as_ref().expect("Adam optimiser");
    assert_eq!(adam.m.len(), params);
    assert_eq!(adam.v.len(), params);
}

/// Changing the network behind the trainer's back is reported, not
/// silently half-applied.
#[test]
fn test_trainer_reports_inconsistent_state() {
    let mut trainer = GpuTrainer::new(scalar_net(0.0, 0.0), 0.01, OptimizerType::Adam);
    trainer
        .net
        .add_layer(NeuralLayer::new(1, 1, LayerType::Dense));
    let result = trainer.train_step(&[1.0], &[1.0]);
    assert!(
        matches!(result, Err(NeuralNetError::ShapeMismatch { .. })),
        "{result:?}"
    );
}

#[test]
fn test_mse_loss_without_overlap_is_zero() {
    assert_eq!(mse_loss(&[1.0, 2.0], &[]), 0.0);
    assert_eq!(mse_loss(&[], &[1.0]), 0.0);
}

#[test]
fn test_softmax_net_output_probabilities() {
    let mut net = GpuNeuralNet::new();
    net.add_layer(NeuralLayer::new(3, 3, LayerType::Softmax));
    let out = net.forward_pass(&[0.0, 1.0, 2.0]);
    let sum: f64 = out.iter().sum();
    assert!((sum - 1.0).abs() < 1e-12);
    for &p in &out {
        assert!((0.0..=1.0).contains(&p));
    }
}

#[test]
fn test_mse_symmetric() {
    let a = vec![1.0, 2.0];
    let b = vec![3.0, 4.0];
    assert!((mse_loss(&a, &b) - mse_loss(&b, &a)).abs() < 1e-12);
}

#[test]
fn test_layer_type_debug() {
    let lt = LayerType::Dense;
    let s = format!("{lt:?}");
    assert!(s.contains("Dense"));
}

#[test]
fn test_optimizer_type_eq() {
    assert_eq!(OptimizerType::Sgd, OptimizerType::Sgd);
    assert_ne!(OptimizerType::Sgd, OptimizerType::Adam);
}

#[test]
fn test_sigmoid_vs_exp() {
    // Verify sigmoid matches the direct formula
    let x = 1.0_f64;
    assert!((sigmoid(x) - 1.0 / (1.0 + (-x).exp())).abs() < 1e-12);
}
