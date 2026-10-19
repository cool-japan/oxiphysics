// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated neural network compute (CPU mock backend).
//!
//! Provides layer-wise forward passes, exact backpropagation (weight, bias and
//! input gradients for every [`LayerType`]), SGD and Adam training — all
//! running on the CPU as a mock GPU backend.
//!
//! # Backpropagation
//!
//! [`GpuNeuralNet::forward_with_cache`] records every layer's input and output
//! in a [`ForwardCache`]; [`BackpropGpu::backward_pass`] then applies the chain
//! rule layer by layer (Dense `dW = δ xᵀ`, `δ_prev = Wᵀ δ`; ReLU `δ·1[x>0]`;
//! Sigmoid `δ·s(1−s)`; Tanh `δ·(1−t²)`; Softmax via its Jacobian-vector
//! product; BatchNorm with trainable γ/β; Conv1D with kernel size and stride).

use thiserror::Error;

// ---------------------------------------------------------------------------
// Activation helpers (public, free functions)
// ---------------------------------------------------------------------------

/// Rectified linear unit: `max(0, x)`.
pub fn relu(x: f64) -> f64 {
    x.max(0.0)
}

/// Logistic sigmoid: `1 / (1 + e^{-x})`.
pub fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Softmax of a slice: `exp(x_i) / sum(exp(x_j))`.
///
/// Numerically stable implementation via max-subtraction.
pub fn softmax(x: &[f64]) -> Vec<f64> {
    if x.is_empty() {
        return Vec::new();
    }
    let max_val = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = x.iter().map(|&v| (v - max_val).exp()).collect();
    let sum: f64 = exps.iter().sum();
    exps.iter().map(|&e| e / sum).collect()
}

/// Mean-squared error: `mean((pred_i - target_i)^2)`.
///
/// Compares the first `min(pred.len(), target.len())` entries; returns `0.0`
/// when there is nothing to compare.
pub fn mse_loss(pred: &[f64], target: &[f64]) -> f64 {
    let n = pred.len().min(target.len());
    if n == 0 {
        return 0.0;
    }
    let sum: f64 = pred[..n]
        .iter()
        .zip(target[..n].iter())
        .map(|(p, t)| (p - t).powi(2))
        .sum();
    sum / n as f64
}

// ---------------------------------------------------------------------------
// NeuralNetError
// ---------------------------------------------------------------------------

/// Errors reported by the network construction and backpropagation routines.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum NeuralNetError {
    /// A buffer does not have the length the operation requires.
    #[error("shape mismatch in {context}: expected {expected}, got {got}")]
    ShapeMismatch {
        /// What was being checked (e.g. `"layer 2 weight gradient"`).
        context: String,
        /// Required length.
        expected: usize,
        /// Actual length.
        got: usize,
    },
    /// A layer was configured with impossible parameters.
    #[error("invalid layer configuration: {0}")]
    InvalidLayer(String),
}

/// Check that `got == expected`, naming `context` in the error otherwise.
fn expect_len(context: &str, expected: usize, got: usize) -> Result<(), NeuralNetError> {
    if expected == got {
        Ok(())
    } else {
        Err(NeuralNetError::ShapeMismatch {
            context: context.to_string(),
            expected,
            got,
        })
    }
}

// ---------------------------------------------------------------------------
// LayerType
// ---------------------------------------------------------------------------

/// The computational type of a single neural network layer.
#[derive(Debug, Clone, PartialEq)]
pub enum LayerType {
    /// Fully-connected (dense) layer.
    Dense,
    /// 1-D convolution layer (see [`Conv1dSpec`] for the geometry).
    Conv1D,
    /// Rectified linear unit activation.
    ReLU,
    /// Sigmoid activation.
    Sigmoid,
    /// Hyperbolic tangent activation.
    Tanh,
    /// Softmax activation.
    Softmax,
    /// Batch normalisation layer.
    BatchNorm,
    /// Dropout regularisation layer.
    Dropout,
}

// ---------------------------------------------------------------------------
// Conv1dSpec
// ---------------------------------------------------------------------------

/// Geometry of a [`LayerType::Conv1D`] layer.
///
/// The layer cross-correlates `in_channels` input channels of `input_len`
/// samples (channel-major: `x[ic * input_len + t]`) with `out_channels`
/// filters of `kernel_size` taps, advancing `stride` samples between output
/// positions and applying no padding:
///
/// `y[oc * positions + p] = b[oc] + Σ_ic Σ_k w[(oc * in_channels + ic) * kernel_size + k] · x[ic * input_len + p * stride + k]`
///
/// with `positions = (input_len − kernel_size) / stride + 1`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Conv1dSpec {
    /// Number of input channels.
    pub in_channels: usize,
    /// Number of output channels (filters).
    pub out_channels: usize,
    /// Samples per input channel.
    pub input_len: usize,
    /// Taps per filter and channel.
    pub kernel_size: usize,
    /// Step between successive kernel positions.
    pub stride: usize,
}

impl Conv1dSpec {
    /// Validate and build a convolution geometry.
    ///
    /// # Errors
    ///
    /// Returns [`NeuralNetError::InvalidLayer`] when a channel count, the
    /// kernel size or the stride is zero, or when the kernel is longer than
    /// the input.
    pub fn new(
        in_channels: usize,
        out_channels: usize,
        input_len: usize,
        kernel_size: usize,
        stride: usize,
    ) -> Result<Self, NeuralNetError> {
        if in_channels == 0 || out_channels == 0 {
            return Err(NeuralNetError::InvalidLayer(format!(
                "Conv1D needs positive channel counts, got {in_channels} -> {out_channels}"
            )));
        }
        if kernel_size == 0 || stride == 0 {
            return Err(NeuralNetError::InvalidLayer(format!(
                "Conv1D needs a positive kernel size and stride, got kernel {kernel_size}, \
                 stride {stride}"
            )));
        }
        if kernel_size > input_len {
            return Err(NeuralNetError::InvalidLayer(format!(
                "Conv1D kernel size {kernel_size} exceeds the input length {input_len}"
            )));
        }
        Ok(Self {
            in_channels,
            out_channels,
            input_len,
            kernel_size,
            stride,
        })
    }

    /// Single-channel geometry whose kernel spans the whole input: one output
    /// position per output channel (the geometry of [`NeuralLayer::new`]).
    fn full_width(input_len: usize, out_channels: usize) -> Self {
        Self {
            in_channels: 1,
            out_channels,
            input_len,
            kernel_size: input_len,
            stride: 1,
        }
    }

    /// Output positions per channel; `0` for a degenerate geometry (zero
    /// kernel or stride, or a kernel longer than the input).
    pub fn output_len(&self) -> usize {
        if self.stride == 0 || self.kernel_size == 0 || self.kernel_size > self.input_len {
            0
        } else {
            (self.input_len - self.kernel_size) / self.stride + 1
        }
    }

    /// Length of the flattened input, `in_channels * input_len`.
    pub fn input_size(&self) -> usize {
        self.in_channels * self.input_len
    }

    /// Length of the flattened output, `out_channels * output_len()`.
    pub fn output_size(&self) -> usize {
        self.out_channels * self.output_len()
    }

    /// Number of filter weights, `out_channels * in_channels * kernel_size`.
    pub fn weight_len(&self) -> usize {
        self.out_channels * self.in_channels * self.kernel_size
    }
}

// ---------------------------------------------------------------------------
// NeuralLayer
// ---------------------------------------------------------------------------

/// Variance regulariser of the normalisation layer.
const BATCH_NORM_EPSILON: f64 = 1e-5;

/// A single layer in a neural network, carrying weights, biases and a type.
#[derive(Debug, Clone)]
pub struct NeuralLayer {
    /// Flattened weight matrix (row-major: `[out, in]`).  For a
    /// [`LayerType::Conv1D`] layer: the filter bank `[out_ch, in_ch, kernel]`.
    /// For [`LayerType::BatchNorm`]: the scale γ in the first `input_size`
    /// entries.
    pub weights: Vec<f64>,
    /// Bias vector (length = number of output neurons; one per output channel
    /// for `Conv1D`; the shift β for `BatchNorm`).
    pub biases: Vec<f64>,
    /// Computational type of this layer.
    pub layer_type: LayerType,
    /// Number of input neurons / features.
    pub input_size: usize,
    /// Number of output neurons.
    pub output_size: usize,
    /// Convolution geometry of a [`LayerType::Conv1D`] layer, `None` for every
    /// other layer type (see [`NeuralLayer::conv_spec`]).
    pub conv: Option<Conv1dSpec>,
}

impl NeuralLayer {
    /// Create a new layer with given dimensions and type.
    ///
    /// Weights and biases are zero-initialised; assign the public `weights`
    /// and `biases` fields to set custom values.
    ///
    /// A [`LayerType::Conv1D`] layer created here has `output_size` filters
    /// whose single-channel kernel spans the whole input (one output position
    /// per filter); use [`NeuralLayer::conv1d`] for a sliding kernel with an
    /// explicit kernel size, stride and channel counts.
    pub fn new(input_size: usize, output_size: usize, layer_type: LayerType) -> Self {
        let conv = (layer_type == LayerType::Conv1D)
            .then(|| Conv1dSpec::full_width(input_size, output_size));
        Self {
            weights: vec![0.0; input_size * output_size],
            biases: vec![0.0; output_size],
            layer_type,
            input_size,
            output_size,
            conv,
        }
    }

    /// Create a zero-initialised 1-D convolution layer.
    ///
    /// The input is `in_channels` channels of `input_len` samples
    /// (channel-major); the output is `out_channels` channels of
    /// `(input_len − kernel_size) / stride + 1` positions.
    ///
    /// # Errors
    ///
    /// Returns [`NeuralNetError::InvalidLayer`] for a zero channel count,
    /// kernel size or stride, or a kernel longer than the input.
    pub fn conv1d(
        in_channels: usize,
        out_channels: usize,
        input_len: usize,
        kernel_size: usize,
        stride: usize,
    ) -> Result<Self, NeuralNetError> {
        let spec = Conv1dSpec::new(in_channels, out_channels, input_len, kernel_size, stride)?;
        Ok(Self {
            weights: vec![0.0; spec.weight_len()],
            biases: vec![0.0; out_channels],
            layer_type: LayerType::Conv1D,
            input_size: spec.input_size(),
            output_size: spec.output_size(),
            conv: Some(spec),
        })
    }

    /// Convolution geometry of a `Conv1D` layer; `None` for other layer types.
    ///
    /// A `Conv1D` layer whose [`conv`](Self::conv) field is `None` is read as
    /// the full-width single-channel geometry of [`NeuralLayer::new`].
    pub fn conv_spec(&self) -> Option<Conv1dSpec> {
        if self.layer_type != LayerType::Conv1D {
            return None;
        }
        Some(
            self.conv
                .unwrap_or_else(|| Conv1dSpec::full_width(self.input_size, self.output_size)),
        )
    }

    /// Execute the forward pass of this layer on `input`.
    ///
    /// For activation layers (`ReLU`, `Sigmoid`, `Tanh`, `Softmax`) the
    /// weights/biases are ignored and the input is transformed element-wise.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        match self.layer_type {
            LayerType::ReLU => input.iter().map(|&x| relu(x)).collect(),
            LayerType::Sigmoid => input.iter().map(|&x| sigmoid(x)).collect(),
            LayerType::Tanh => input.iter().map(|&x| x.tanh()).collect(),
            LayerType::Softmax => softmax(input),
            LayerType::BatchNorm => {
                // Normalise the feature vector to zero-mean / unit-variance and
                // scale/shift it with the stored (gamma, beta) pairs.
                let n = input.len();
                if n == 0 {
                    return Vec::new();
                }
                let mean = input.iter().sum::<f64>() / n as f64;
                let var = input.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
                let std = (var + BATCH_NORM_EPSILON).sqrt();
                input
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| {
                        let gamma = self.weights.get(i).copied().unwrap_or(1.0);
                        let beta = self.biases.get(i).copied().unwrap_or(0.0);
                        gamma * (x - mean) / std + beta
                    })
                    .collect()
            }
            LayerType::Dropout => {
                // At inference time Dropout is a no-op.
                input.to_vec()
            }
            LayerType::Dense => {
                // Matrix–vector multiply: out[j] = sum_i w[j*in + i]*input[i] + bias[j]
                let in_sz = input.len();
                let out_sz = self.output_size;
                let mut out = vec![0.0; out_sz];
                for (j, out_val) in out.iter_mut().enumerate() {
                    let mut acc = self.biases.get(j).copied().unwrap_or(0.0);
                    for (i, &inp) in input.iter().enumerate() {
                        let w = self.weights.get(j * in_sz + i).copied().unwrap_or(0.0);
                        acc += w * inp;
                    }
                    *out_val = acc;
                }
                out
            }
            LayerType::Conv1D => self.conv1d_forward(input),
        }
    }

    /// Sliding-window cross-correlation of a `Conv1D` layer (see
    /// [`Conv1dSpec`]).  Missing input samples or weights read as zero.
    fn conv1d_forward(&self, input: &[f64]) -> Vec<f64> {
        let Some(spec) = self.conv_spec() else {
            return Vec::new();
        };
        let positions = spec.output_len();
        let mut out = vec![0.0; spec.out_channels * positions];
        for oc in 0..spec.out_channels {
            let bias = self.biases.get(oc).copied().unwrap_or(0.0);
            for p in 0..positions {
                let mut acc = bias;
                for ic in 0..spec.in_channels {
                    for k in 0..spec.kernel_size {
                        let w_idx = (oc * spec.in_channels + ic) * spec.kernel_size + k;
                        let x_idx = ic * spec.input_len + p * spec.stride + k;
                        let w = self.weights.get(w_idx).copied().unwrap_or(0.0);
                        let x = input.get(x_idx).copied().unwrap_or(0.0);
                        acc += w * x;
                    }
                }
                if let Some(slot) = out.get_mut(oc * positions + p) {
                    *slot = acc;
                }
            }
        }
        out
    }

    /// Backward pass of this layer.
    ///
    /// `input` and `output` are the values the layer received and produced in
    /// the forward pass, `grad_output` is `∂L/∂output`.  The weight and bias
    /// gradients are *added* to `grad_weights` / `grad_biases` (so several
    /// samples can be accumulated); the returned vector is `∂L/∂input`.
    ///
    /// * Dense: `dW[j·in + i] += δ_j x_i`, `db_j += δ_j`, `δ_prev = Wᵀ δ`.
    /// * Conv1D: the correlation of `δ` with the input per filter tap, with
    ///   `δ_prev` scattered back through the same windows.
    /// * ReLU: `δ · 1[x > 0]`; Sigmoid: `δ · s (1 − s)`; Tanh: `δ · (1 − t²)`.
    /// * Softmax: `s ⊙ (δ − ⟨δ, s⟩)` (the Jacobian-vector product).
    /// * BatchNorm: `dγ_i += δ_i x̂_i`, `dβ_i += δ_i` and the full
    ///   normalisation Jacobian for the input gradient.
    /// * Dropout: the identity (it is a no-op in [`NeuralLayer::forward`]).
    ///
    /// # Errors
    ///
    /// Returns [`NeuralNetError::ShapeMismatch`] when `output`, `grad_output`,
    /// `grad_weights` or `grad_biases` do not have the lengths this layer
    /// produces.
    pub fn backward(
        &self,
        input: &[f64],
        output: &[f64],
        grad_output: &[f64],
        grad_weights: &mut [f64],
        grad_biases: &mut [f64],
    ) -> Result<Vec<f64>, NeuralNetError> {
        expect_len("output gradient", output.len(), grad_output.len())?;
        expect_len("weight gradient", self.weights.len(), grad_weights.len())?;
        expect_len("bias gradient", self.biases.len(), grad_biases.len())?;
        match self.layer_type {
            LayerType::Dense => {
                expect_len("dense output", self.output_size, output.len())?;
                Ok(self.dense_backward(input, grad_output, grad_weights, grad_biases))
            }
            LayerType::Conv1D => {
                let expected = self.conv_spec().map_or(0, |spec| spec.output_size());
                expect_len("conv1d output", expected, output.len())?;
                Ok(self.conv1d_backward(input, grad_output, grad_weights, grad_biases))
            }
            LayerType::ReLU => {
                expect_len("activation input", output.len(), input.len())?;
                Ok(input
                    .iter()
                    .zip(grad_output)
                    .map(|(&x, &g)| if x > 0.0 { g } else { 0.0 })
                    .collect())
            }
            LayerType::Sigmoid => {
                expect_len("activation input", output.len(), input.len())?;
                Ok(output
                    .iter()
                    .zip(grad_output)
                    .map(|(&s, &g)| g * s * (1.0 - s))
                    .collect())
            }
            LayerType::Tanh => {
                expect_len("activation input", output.len(), input.len())?;
                Ok(output
                    .iter()
                    .zip(grad_output)
                    .map(|(&t, &g)| g * (1.0 - t * t))
                    .collect())
            }
            LayerType::Softmax => {
                expect_len("activation input", output.len(), input.len())?;
                let dot: f64 = output.iter().zip(grad_output).map(|(&s, &g)| s * g).sum();
                Ok(output
                    .iter()
                    .zip(grad_output)
                    .map(|(&s, &g)| s * (g - dot))
                    .collect())
            }
            LayerType::Dropout => {
                expect_len("dropout input", output.len(), input.len())?;
                Ok(grad_output.to_vec())
            }
            LayerType::BatchNorm => {
                expect_len("batch-norm input", output.len(), input.len())?;
                Ok(self.batch_norm_backward(input, grad_output, grad_weights, grad_biases))
            }
        }
    }

    /// Dense-layer backward pass (strides follow [`NeuralLayer::forward`]).
    fn dense_backward(
        &self,
        input: &[f64],
        grad_output: &[f64],
        grad_weights: &mut [f64],
        grad_biases: &mut [f64],
    ) -> Vec<f64> {
        let in_sz = input.len();
        let mut grad_input = vec![0.0; in_sz];
        for (j, &g) in grad_output.iter().enumerate() {
            if let Some(gb) = grad_biases.get_mut(j) {
                *gb += g;
            }
            for (i, (&x, gx)) in input.iter().zip(grad_input.iter_mut()).enumerate() {
                let idx = j * in_sz + i;
                if let Some(gw) = grad_weights.get_mut(idx) {
                    *gw += g * x;
                }
                if let Some(&w) = self.weights.get(idx) {
                    *gx += w * g;
                }
            }
        }
        grad_input
    }

    /// Convolution backward pass: the exact adjoint of
    /// [`NeuralLayer::conv1d_forward`].
    fn conv1d_backward(
        &self,
        input: &[f64],
        grad_output: &[f64],
        grad_weights: &mut [f64],
        grad_biases: &mut [f64],
    ) -> Vec<f64> {
        let mut grad_input = vec![0.0; input.len()];
        let Some(spec) = self.conv_spec() else {
            return grad_input;
        };
        let positions = spec.output_len();
        for oc in 0..spec.out_channels {
            for p in 0..positions {
                let g = grad_output.get(oc * positions + p).copied().unwrap_or(0.0);
                if let Some(gb) = grad_biases.get_mut(oc) {
                    *gb += g;
                }
                for ic in 0..spec.in_channels {
                    for k in 0..spec.kernel_size {
                        let w_idx = (oc * spec.in_channels + ic) * spec.kernel_size + k;
                        let x_idx = ic * spec.input_len + p * spec.stride + k;
                        let Some(&x) = input.get(x_idx) else {
                            continue;
                        };
                        if let Some(gw) = grad_weights.get_mut(w_idx) {
                            *gw += g * x;
                        }
                        if let (Some(&w), Some(gx)) =
                            (self.weights.get(w_idx), grad_input.get_mut(x_idx))
                        {
                            *gx += w * g;
                        }
                    }
                }
            }
        }
        grad_input
    }

    /// Backward pass of the feature-wise normalisation layer.
    ///
    /// With `x̂ = (x − μ)/σ`, `y = γ x̂ + β` and `σ² = var + ε`:
    /// `dγ = δ x̂`, `dβ = δ` and
    /// `dx = (1/σ) (δγ − mean(δγ) − x̂ · mean(δγ ⊙ x̂))`.
    fn batch_norm_backward(
        &self,
        input: &[f64],
        grad_output: &[f64],
        grad_weights: &mut [f64],
        grad_biases: &mut [f64],
    ) -> Vec<f64> {
        let n = input.len();
        if n == 0 {
            return Vec::new();
        }
        let count = n as f64;
        let mean = input.iter().sum::<f64>() / count;
        let var = input.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / count;
        let std = (var + BATCH_NORM_EPSILON).sqrt();
        let normalised: Vec<f64> = input.iter().map(|&x| (x - mean) / std).collect();

        let mut grad_normalised = Vec::with_capacity(n);
        for (i, (&g, &x_hat)) in grad_output.iter().zip(normalised.iter()).enumerate() {
            let gamma = self.weights.get(i).copied().unwrap_or(1.0);
            if let Some(gw) = grad_weights.get_mut(i) {
                *gw += g * x_hat;
            }
            if let Some(gb) = grad_biases.get_mut(i) {
                *gb += g;
            }
            grad_normalised.push(g * gamma);
        }
        let mean_grad = grad_normalised.iter().sum::<f64>() / count;
        let mean_grad_x = grad_normalised
            .iter()
            .zip(normalised.iter())
            .map(|(&g, &x_hat)| g * x_hat)
            .sum::<f64>()
            / count;
        grad_normalised
            .iter()
            .zip(normalised.iter())
            .map(|(&g, &x_hat)| (g - mean_grad - x_hat * mean_grad_x) / std)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// GpuNeuralNet
// ---------------------------------------------------------------------------

/// A sequential neural network backed by a CPU mock GPU context.
#[derive(Debug, Clone)]
pub struct GpuNeuralNet {
    /// Ordered list of layers in the network.
    pub layers: Vec<NeuralLayer>,
}

impl GpuNeuralNet {
    /// Create an empty network with no layers.
    pub fn new() -> Self {
        Self { layers: Vec::new() }
    }

    /// Append a layer to the network.
    pub fn add_layer(&mut self, layer: NeuralLayer) {
        self.layers.push(layer);
    }

    /// Run a single forward pass through all layers.
    pub fn forward_pass(&self, input: &[f64]) -> Vec<f64> {
        let mut current = input.to_vec();
        for layer in &self.layers {
            current = layer.forward(&current);
        }
        current
    }

    /// Run a forward pass that records every layer's input and output, as
    /// needed by [`BackpropGpu::backward_pass`].
    ///
    /// [`ForwardCache::output`] equals [`forward_pass`](Self::forward_pass).
    pub fn forward_with_cache(&self, input: &[f64]) -> ForwardCache {
        let mut layer_outputs: Vec<Vec<f64>> = Vec::with_capacity(self.layers.len());
        for layer in &self.layers {
            let next = layer.forward(layer_outputs.last().map_or(input, Vec::as_slice));
            layer_outputs.push(next);
        }
        ForwardCache {
            input: input.to_vec(),
            layer_outputs,
        }
    }

    /// Run forward passes for a batch of inputs.
    pub fn batch_forward(&self, inputs: &[Vec<f64>]) -> Vec<Vec<f64>> {
        inputs.iter().map(|inp| self.forward_pass(inp)).collect()
    }
}

impl Default for GpuNeuralNet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ForwardCache
// ---------------------------------------------------------------------------

/// Activations recorded by [`GpuNeuralNet::forward_with_cache`].
///
/// Layer `l` received [`layer_input(l)`](Self::layer_input) and produced
/// `layer_outputs()[l]`; the network output is the last layer's output (the
/// input itself for a network without layers).
#[derive(Debug, Clone)]
pub struct ForwardCache {
    /// The network input.
    input: Vec<f64>,
    /// `layer_outputs[l]` is the output of layer `l`.
    layer_outputs: Vec<Vec<f64>>,
}

impl ForwardCache {
    /// The network input the cache was recorded for.
    pub fn input(&self) -> &[f64] {
        &self.input
    }

    /// The output of every layer, in network order.
    pub fn layer_outputs(&self) -> &[Vec<f64>] {
        &self.layer_outputs
    }

    /// The input that layer `layer` received, or `None` past the last layer.
    pub fn layer_input(&self, layer: usize) -> Option<&[f64]> {
        if layer == 0 {
            Some(&self.input)
        } else {
            self.layer_outputs.get(layer - 1).map(Vec::as_slice)
        }
    }

    /// The network output (the input itself for a network without layers).
    pub fn output(&self) -> &[f64] {
        self.layer_outputs.last().map_or(&self.input, Vec::as_slice)
    }
}

// ---------------------------------------------------------------------------
// BackpropGpu
// ---------------------------------------------------------------------------

/// Backpropagation state: stores per-layer weight and bias gradients.
#[derive(Debug, Clone)]
pub struct BackpropGpu {
    /// Per-layer weight gradients (same shape as the layer's weight vector).
    pub gradients: Vec<Vec<f64>>,
    /// Per-layer bias gradients (same shape as the layer's bias vector).
    pub bias_gradients: Vec<Vec<f64>>,
}

impl BackpropGpu {
    /// Initialise gradient buffers matching the network's weight and bias
    /// shapes.
    pub fn new(net: &GpuNeuralNet) -> Self {
        let gradients = net
            .layers
            .iter()
            .map(|l| vec![0.0; l.weights.len()])
            .collect();
        let bias_gradients = net
            .layers
            .iter()
            .map(|l| vec![0.0; l.biases.len()])
            .collect();
        Self {
            gradients,
            bias_gradients,
        }
    }

    /// Reset every gradient to zero.
    pub fn zero_grad(&mut self) {
        for buffer in self
            .gradients
            .iter_mut()
            .chain(self.bias_gradients.iter_mut())
        {
            buffer.fill(0.0);
        }
    }

    /// Backpropagate the loss gradient `loss_grad` (`∂L/∂output`) through the
    /// network, replacing the stored gradients with those of this sample.
    ///
    /// `cache` must come from [`GpuNeuralNet::forward_with_cache`] on the same
    /// `net`.  Returns `∂L/∂input`.
    ///
    /// # Errors
    ///
    /// Returns [`NeuralNetError::ShapeMismatch`] when the cache, the gradient
    /// buffers or `loss_grad` do not match the network.
    pub fn backward_pass(
        &mut self,
        net: &GpuNeuralNet,
        cache: &ForwardCache,
        loss_grad: &[f64],
    ) -> Result<Vec<f64>, NeuralNetError> {
        self.zero_grad();
        self.accumulate_backward(net, cache, loss_grad)
    }

    /// Like [`backward_pass`](Self::backward_pass) but *adds* to the stored
    /// gradients, which accumulates a mini-batch (call
    /// [`zero_grad`](Self::zero_grad) before the first sample).
    ///
    /// # Errors
    ///
    /// Returns [`NeuralNetError::ShapeMismatch`] when the cache, the gradient
    /// buffers or `loss_grad` do not match the network.
    pub fn accumulate_backward(
        &mut self,
        net: &GpuNeuralNet,
        cache: &ForwardCache,
        loss_grad: &[f64],
    ) -> Result<Vec<f64>, NeuralNetError> {
        let layer_count = net.layers.len();
        expect_len("gradient buffers", layer_count, self.gradients.len())?;
        expect_len(
            "bias gradient buffers",
            layer_count,
            self.bias_gradients.len(),
        )?;
        expect_len(
            "cached layer outputs",
            layer_count,
            cache.layer_outputs.len(),
        )?;
        expect_len("loss gradient", cache.output().len(), loss_grad.len())?;

        let mut grad = loss_grad.to_vec();
        for index in (0..layer_count).rev() {
            let (Some(layer), Some(input), Some(output), Some(grad_w), Some(grad_b)) = (
                net.layers.get(index),
                cache.layer_input(index),
                cache.layer_outputs.get(index),
                self.gradients.get_mut(index),
                self.bias_gradients.get_mut(index),
            ) else {
                return Err(NeuralNetError::ShapeMismatch {
                    context: format!("layer {index} backward buffers"),
                    expected: layer_count,
                    got: index,
                });
            };
            grad = layer.backward(input, output, &grad, grad_w, grad_b)?;
        }
        Ok(grad)
    }
}

// ---------------------------------------------------------------------------
// Optimizer type
// ---------------------------------------------------------------------------

/// Optimiser selection for the GPU trainer.
#[derive(Debug, Clone, PartialEq)]
pub enum OptimizerType {
    /// Stochastic gradient descent.
    Sgd,
    /// Adaptive moment estimation (Adam).
    Adam,
}

// ---------------------------------------------------------------------------
// AdamOptimizer
// ---------------------------------------------------------------------------

/// Adam adaptive moment estimator.
///
/// Reference: Kingma & Ba (2015) — "Adam: A Method for Stochastic Optimization".
#[derive(Debug, Clone)]
pub struct AdamOptimizer {
    /// First moment decay rate (typically 0.9).
    pub beta1: f64,
    /// Second moment decay rate (typically 0.999).
    pub beta2: f64,
    /// Numerical stability constant (typically 1e-8).
    pub eps: f64,
    /// Learning rate.
    pub lr: f64,
    /// First moment (mean) buffer.
    pub m: Vec<f64>,
    /// Second moment (variance) buffer.
    pub v: Vec<f64>,
    /// Current time-step (number of update calls so far).
    pub t: u64,
}

impl AdamOptimizer {
    /// Create a new Adam optimiser for `n` parameters.
    pub fn new(n: usize, lr: f64, beta1: f64, beta2: f64, eps: f64) -> Self {
        Self {
            beta1,
            beta2,
            eps,
            lr,
            m: vec![0.0; n],
            v: vec![0.0; n],
            t: 0,
        }
    }

    /// Apply one Adam update step.
    ///
    /// * `params` — mutable slice of parameter values.
    /// * `grads`  — gradient slice of the same length.
    pub fn update(&mut self, params: &mut [f64], grads: &[f64]) {
        self.t += 1;
        let t = self.t as f64;
        let lr_t = self.lr * (1.0 - self.beta2.powf(t)).sqrt() / (1.0 - self.beta1.powf(t));
        let n = params
            .len()
            .min(grads.len())
            .min(self.m.len())
            .min(self.v.len());
        for i in 0..n {
            self.m[i] = self.beta1 * self.m[i] + (1.0 - self.beta1) * grads[i];
            self.v[i] = self.beta2 * self.v[i] + (1.0 - self.beta2) * grads[i].powi(2);
            params[i] -= lr_t * self.m[i] / (self.v[i].sqrt() + self.eps);
        }
    }
}

// ---------------------------------------------------------------------------
// GpuTrainer
// ---------------------------------------------------------------------------

/// Combines a network, a backprop context, and an optimiser for training.
#[derive(Debug)]
pub struct GpuTrainer {
    /// Neural network being trained.
    pub net: GpuNeuralNet,
    /// Backpropagation state.
    pub backprop: BackpropGpu,
    /// Learning rate.
    pub learning_rate: f64,
    /// Which optimiser to use.
    pub optimizer: OptimizerType,
    /// Adam optimiser instance (only active when `optimizer == Adam`); sized
    /// for every weight and bias of the network, layer by layer
    /// (`w₀, b₀, w₁, b₁, …`).
    pub adam: Option<AdamOptimizer>,
}

impl GpuTrainer {
    /// Create a new trainer wrapping `net` with the given optimiser.
    pub fn new(net: GpuNeuralNet, learning_rate: f64, optimizer: OptimizerType) -> Self {
        let backprop = BackpropGpu::new(&net);
        let total_params: usize = net
            .layers
            .iter()
            .map(|l| l.weights.len() + l.biases.len())
            .sum();
        let adam = if optimizer == OptimizerType::Adam {
            Some(AdamOptimizer::new(
                total_params,
                learning_rate,
                0.9,
                0.999,
                1e-8,
            ))
        } else {
            None
        };
        Self {
            net,
            backprop,
            learning_rate,
            optimizer,
            adam,
        }
    }

    /// Execute one training step: forward pass → loss → backward → update.
    ///
    /// * `input`  — network input.
    /// * `target` — ground-truth target.
    ///
    /// The mean-squared-error gradient `2 (pred − target) / n` is
    /// backpropagated through the whole network and every weight *and bias*
    /// is updated with SGD or Adam.  Returns the MSE loss before the update.
    ///
    /// # Errors
    ///
    /// Returns [`NeuralNetError::ShapeMismatch`] when the network, the
    /// gradient buffers or the optimiser were changed after the trainer was
    /// created so that their shapes no longer agree.
    pub fn train_step(&mut self, input: &[f64], target: &[f64]) -> Result<f64, NeuralNetError> {
        // Forward
        let cache = self.net.forward_with_cache(input);
        let pred = cache.output();
        let loss = mse_loss(pred, target);

        // Output gradient of the MSE: 2 (pred - target) / n
        let n = pred.len().min(target.len());
        let mut loss_grad = vec![0.0; pred.len()];
        for ((g, &p), &t) in loss_grad.iter_mut().zip(pred.iter()).zip(target.iter()) {
            *g = 2.0 * (p - t) / n as f64;
        }

        // Backward
        self.backprop.backward_pass(&self.net, &cache, &loss_grad)?;

        // Update weights and biases with SGD or Adam
        match self.optimizer {
            OptimizerType::Sgd => {
                let lr = self.learning_rate;
                for ((layer, grad_w), grad_b) in self
                    .net
                    .layers
                    .iter_mut()
                    .zip(self.backprop.gradients.iter())
                    .zip(self.backprop.bias_gradients.iter())
                {
                    for (w, &g) in layer.weights.iter_mut().zip(grad_w.iter()) {
                        *w -= lr * g;
                    }
                    for (b, &g) in layer.biases.iter_mut().zip(grad_b.iter()) {
                        *b -= lr * g;
                    }
                }
            }
            OptimizerType::Adam => {
                if let Some(adam) = &mut self.adam {
                    // Flatten all weights and biases into a single buffer,
                    // update, scatter back.
                    let mut params: Vec<f64> = Vec::with_capacity(adam.m.len());
                    let mut grads: Vec<f64> = Vec::with_capacity(adam.m.len());
                    for ((layer, grad_w), grad_b) in self
                        .net
                        .layers
                        .iter()
                        .zip(self.backprop.gradients.iter())
                        .zip(self.backprop.bias_gradients.iter())
                    {
                        params.extend_from_slice(&layer.weights);
                        params.extend_from_slice(&layer.biases);
                        grads.extend_from_slice(grad_w);
                        grads.extend_from_slice(grad_b);
                    }
                    expect_len("Adam parameter buffer", adam.m.len(), params.len())?;
                    adam.update(&mut params, &grads);
                    // Scatter back in the same weights-then-biases order.
                    let mut updated = params.into_iter();
                    for layer in &mut self.net.layers {
                        for value in layer.weights.iter_mut().chain(layer.biases.iter_mut()) {
                            if let Some(next) = updated.next() {
                                *value = next;
                            }
                        }
                    }
                }
            }
        }

        Ok(loss)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests;
