//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use std::collections::HashMap;
use std::f64::consts::PI as PI_F64;

use super::functions::scaled_dot_product_attention;
#[allow(unused_imports)]
use super::functions::*;

/// Batch normalization layer (inference mode).
///
/// Normalizes input features using stored running mean and variance,
/// then applies learned scale (gamma) and shift (beta).
#[derive(Debug, Clone)]
pub struct BatchNormLayer {
    /// Running mean for each feature.
    pub running_mean: Vec<f32>,
    /// Running variance for each feature.
    pub running_var: Vec<f32>,
    /// Learned scale parameter (gamma).
    pub gamma: Vec<f32>,
    /// Learned shift parameter (beta).
    pub beta: Vec<f32>,
    /// Small constant for numerical stability.
    pub epsilon: f32,
    /// Number of features.
    pub n_features: usize,
}
impl BatchNormLayer {
    /// Create a new batch norm layer with identity transform (gamma=1, beta=0).
    pub fn new(n_features: usize) -> Self {
        Self {
            running_mean: vec![0.0; n_features],
            running_var: vec![1.0; n_features],
            gamma: vec![1.0; n_features],
            beta: vec![0.0; n_features],
            epsilon: 1e-5,
            n_features,
        }
    }
    /// Apply batch normalization in inference mode.
    ///
    /// output\[i\] = gamma\[i\] * (input\[i\] - mean\[i\]) / sqrt(var\[i\] + eps) + beta\[i\]
    pub fn forward(&self, input: &[f32]) -> Vec<f32> {
        assert_eq!(input.len(), self.n_features);
        let mut output = Vec::with_capacity(self.n_features);
        for i in 0..self.n_features {
            let normalized =
                (input[i] - self.running_mean[i]) / (self.running_var[i] + self.epsilon).sqrt();
            output.push(self.gamma[i] * normalized + self.beta[i]);
        }
        output
    }
    /// Set the running statistics.
    pub fn set_stats(&mut self, mean: &[f32], var: &[f32]) {
        assert_eq!(mean.len(), self.n_features);
        assert_eq!(var.len(), self.n_features);
        self.running_mean.copy_from_slice(mean);
        self.running_var.copy_from_slice(var);
    }
    /// Set the affine parameters.
    pub fn set_affine(&mut self, gamma: &[f32], beta: &[f32]) {
        assert_eq!(gamma.len(), self.n_features);
        assert_eq!(beta.len(), self.n_features);
        self.gamma.copy_from_slice(gamma);
        self.beta.copy_from_slice(beta);
    }
}
impl BatchNormLayer {
    /// Update running statistics from a mini-batch (training mode).
    ///
    /// Uses exponential moving average:
    /// `running_mean = (1-momentum) * running_mean + momentum * batch_mean`
    ///
    /// # Panics
    /// Panics if `batch` is empty or if any sample has the wrong feature count.
    pub fn update_running_stats(&mut self, batch: &[Vec<f32>], momentum: f32) {
        assert!(
            !batch.is_empty(),
            "update_running_stats: batch must not be empty"
        );
        let n = batch.len() as f32;
        let mut batch_mean = vec![0.0_f32; self.n_features];
        for sample in batch {
            assert_eq!(sample.len(), self.n_features, "sample length mismatch");
            for (k, &v) in sample.iter().enumerate() {
                batch_mean[k] += v;
            }
        }
        for m in &mut batch_mean {
            *m /= n;
        }
        let mut batch_var = vec![0.0_f32; self.n_features];
        for sample in batch {
            for (k, &v) in sample.iter().enumerate() {
                let d = v - batch_mean[k];
                batch_var[k] += d * d;
            }
        }
        for v in &mut batch_var {
            *v /= n;
        }
        for k in 0..self.n_features {
            self.running_mean[k] =
                (1.0 - momentum) * self.running_mean[k] + momentum * batch_mean[k];
            self.running_var[k] = (1.0 - momentum) * self.running_var[k] + momentum * batch_var[k];
        }
    }
}
/// A single-step Elman RNN cell:
/// `h_t = activation(W_x * x_t + W_h * h_{t-1} + b)`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RnnCell {
    /// Input-to-hidden weight matrix `[hidden_size × input_size]`.
    pub w_x: Vec<f64>,
    /// Hidden-to-hidden weight matrix `[hidden_size × hidden_size]`.
    pub w_h: Vec<f64>,
    /// Bias vector `[hidden_size]`.
    pub b: Vec<f64>,
    /// Input dimensionality.
    pub input_size: usize,
    /// Hidden state dimensionality.
    pub hidden_size: usize,
    /// Activation function.
    pub activation: ExtActivation,
}
impl RnnCell {
    /// Create a new RNN cell with zero weights.
    pub fn new(input_size: usize, hidden_size: usize, activation: ExtActivation) -> Self {
        Self {
            w_x: vec![0.0_f64; hidden_size * input_size],
            w_h: vec![0.0_f64; hidden_size * hidden_size],
            b: vec![0.0_f64; hidden_size],
            input_size,
            hidden_size,
            activation,
        }
    }
    /// One forward step.
    ///
    /// Returns the new hidden state `h_t` of length `hidden_size`.
    pub fn step(&self, x: &[f64], h_prev: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.input_size);
        assert_eq!(h_prev.len(), self.hidden_size);
        let mut h = Vec::with_capacity(self.hidden_size);
        for o in 0..self.hidden_size {
            let mut acc = self.b[o];
            for i in 0..self.input_size {
                acc += self.w_x[o * self.input_size + i] * x[i];
            }
            for i in 0..self.hidden_size {
                acc += self.w_h[o * self.hidden_size + i] * h_prev[i];
            }
            h.push(self.activation.apply(acc));
        }
        h
    }
    /// Run the RNN over a full sequence `[seq_len][input_size]`.
    ///
    /// Returns all hidden states `[seq_len][hidden_size]`.
    pub fn forward_sequence(&self, sequence: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let mut h = vec![0.0_f64; self.hidden_size];
        let mut hidden_states = Vec::with_capacity(sequence.len());
        for x in sequence {
            h = self.step(x, &h);
            hidden_states.push(h.clone());
        }
        hidden_states
    }
}
/// An inference pipeline that chains DenseLayer and BatchNormLayer operations.
#[derive(Debug, Clone)]
pub struct InferencePipeline {
    /// Ordered list of operations.
    pub ops: Vec<InferenceOp>,
}
impl InferencePipeline {
    /// Create an empty pipeline.
    pub fn new() -> Self {
        Self { ops: Vec::new() }
    }
    /// Add an operation to the pipeline.
    pub fn add_op(&mut self, op: InferenceOp) {
        self.ops.push(op);
    }
    /// Run forward pass through all operations.
    pub fn forward(&self, input: &[f32]) -> Vec<f32> {
        let mut current = input.to_vec();
        for op in &self.ops {
            current = match op {
                InferenceOp::Dense(layer) => layer.forward(&current),
                InferenceOp::BatchNorm(bn) => bn.forward(&current),
                InferenceOp::Activation(act) => current.iter().map(|&x| act.apply(x)).collect(),
            };
        }
        current
    }
    /// Total number of trainable parameters.
    pub fn total_parameters(&self) -> usize {
        self.ops
            .iter()
            .map(|op| match op {
                InferenceOp::Dense(layer) => layer.parameter_count(),
                InferenceOp::BatchNorm(bn) => 2 * bn.n_features,
                InferenceOp::Activation(_) => 0,
            })
            .sum()
    }
}
/// Convenience builder for standard AANN architectures.
pub struct NetworkBuilder;
impl NetworkBuilder {
    /// Build a simple element-specific network:
    ///
    /// `n_descriptors → hidden[0] (Tanh) → hidden[1] (Tanh) → … → 1 (Linear)`
    ///
    /// At least one hidden size must be provided.
    pub fn simple_aann(
        n_descriptors: usize,
        hidden_sizes: &[usize],
        _element: u8,
    ) -> FeedForwardNet {
        let mut net = FeedForwardNet::new();
        let mut prev = n_descriptors;
        for &h in hidden_sizes {
            net.add_layer(DenseLayer::new(prev, h, ActivationFn::Tanh));
            prev = h;
        }
        net.add_layer(DenseLayer::new(prev, 1, ActivationFn::Linear));
        net
    }
}
/// Position-wise feed-forward network used inside a transformer block.
///
/// FFN(x) = max(0, x W1 + b1) W2 + b2
///
/// Applied identically to each position in the sequence.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TransformerFfn {
    /// Input / output dimensionality.
    pub d_model: usize,
    /// Inner (hidden) dimensionality.
    pub d_ff: usize,
    /// W1: \[d_ff × d_model\]
    pub w1: Vec<f64>,
    /// b1: \[d_ff\]
    pub b1: Vec<f64>,
    /// W2: \[d_model × d_ff\]
    pub w2: Vec<f64>,
    /// b2: \[d_model\]
    pub b2: Vec<f64>,
}
impl TransformerFfn {
    /// Create with zero weights.
    pub fn new(d_model: usize, d_ff: usize) -> Self {
        Self {
            d_model,
            d_ff,
            w1: vec![0.0_f64; d_ff * d_model],
            b1: vec![0.0_f64; d_ff],
            w2: vec![0.0_f64; d_model * d_ff],
            b2: vec![0.0_f64; d_model],
        }
    }
    /// Forward pass over a sequence `[seq_len × d_model]` (flat row-major).
    pub fn forward(&self, x: &[f64], seq_len: usize) -> Vec<f64> {
        let dm = self.d_model;
        let df = self.d_ff;
        let mut out = vec![0.0_f64; seq_len * dm];
        for t in 0..seq_len {
            let mut hidden = vec![0.0_f64; df];
            for j in 0..df {
                let mut acc = self.b1[j];
                for i in 0..dm {
                    acc += x[t * dm + i] * self.w1[j * dm + i];
                }
                hidden[j] = acc.max(0.0);
            }
            for j in 0..dm {
                let mut acc = self.b2[j];
                for i in 0..df {
                    acc += hidden[i] * self.w2[j * df + i];
                }
                out[t * dm + j] = acc;
            }
        }
        out
    }
}
/// A 1-D convolutional layer operating on a sequence of feature vectors.
///
/// Applies a set of `out_channels` filters, each of length `kernel_size`
/// spanning `in_channels` input channels, using causal (left) padding so the
/// output length equals the input length.
///
/// Layout:
/// - `weights[o][k][c]` = weight for output channel `o`, kernel position `k`,
///   input channel `c`.
/// - `biases[o]` = bias for output channel `o`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Conv1DLayer {
    /// Number of input channels per time step.
    pub in_channels: usize,
    /// Number of output channels per time step.
    pub out_channels: usize,
    /// Kernel (filter) length along the time axis.
    pub kernel_size: usize,
    /// Filter weights: `weights[out_ch][kernel_pos][in_ch]`.
    pub weights: Vec<Vec<Vec<f64>>>,
    /// Bias per output channel.
    pub biases: Vec<f64>,
    /// Activation function applied after convolution.
    pub activation: ExtActivation,
}
impl Conv1DLayer {
    /// Create a new Conv1D layer with zero-initialised weights.
    pub fn new(
        in_channels: usize,
        out_channels: usize,
        kernel_size: usize,
        activation: ExtActivation,
    ) -> Self {
        let weights = vec![vec![vec![0.0_f64; in_channels]; kernel_size]; out_channels];
        let biases = vec![0.0_f64; out_channels];
        Self {
            in_channels,
            out_channels,
            kernel_size,
            weights,
            biases,
            activation,
        }
    }
    /// Forward pass.
    ///
    /// `input` has shape `[seq_len][in_channels]`.  Returns a tensor of shape
    /// `[seq_len][out_channels]` using causal (left-zero) padding.
    pub fn forward(&self, input: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let seq_len = input.len();
        let mut output = vec![vec![0.0_f64; self.out_channels]; seq_len];
        for t in 0..seq_len {
            for o in 0..self.out_channels {
                let mut acc = self.biases[o];
                for k in 0..self.kernel_size {
                    let src_t = t as isize - k as isize;
                    if src_t < 0 {
                        continue;
                    }
                    let src_t = src_t as usize;
                    for c in 0..self.in_channels {
                        acc += self.weights[o][k][c] * input[src_t][c];
                    }
                }
                output[t][o] = self.activation.apply(acc);
            }
        }
        output
    }
    /// Total number of trainable parameters.
    pub fn num_params(&self) -> usize {
        self.out_channels * self.kernel_size * self.in_channels + self.out_channels
    }
}
/// A single operation in the inference pipeline.
#[derive(Debug, Clone)]
pub enum InferenceOp {
    /// Dense (fully-connected) layer.
    Dense(DenseLayer),
    /// Batch normalization layer.
    BatchNorm(BatchNormLayer),
    /// Activation function (standalone).
    Activation(ActivationFn),
}
/// Adam optimizer for a flat parameter vector.
///
/// Reference: Kingma & Ba (2015) "Adam: A Method for Stochastic Optimization".
#[derive(Debug, Clone)]
pub struct AdamOptimizer {
    /// Learning rate α.
    pub lr: f64,
    /// Exponential decay rate for first moment estimates.
    pub beta1: f64,
    /// Exponential decay rate for second moment estimates.
    pub beta2: f64,
    /// Small constant for numerical stability.
    pub epsilon: f64,
    /// First moment vector (m).
    pub m: Vec<f64>,
    /// Second moment vector (v).
    pub v: Vec<f64>,
    /// Current step count (t).
    pub step: u64,
}
impl AdamOptimizer {
    /// Create a new Adam optimizer for a parameter vector of length `n_params`.
    pub fn new(n_params: usize, lr: f64, beta1: f64, beta2: f64, epsilon: f64) -> Self {
        Self {
            lr,
            beta1,
            beta2,
            epsilon,
            m: vec![0.0; n_params],
            v: vec![0.0; n_params],
            step: 0,
        }
    }
    /// Create an Adam optimizer with default hyperparameters (lr=1e-3, β1=0.9, β2=0.999, ε=1e-8).
    pub fn default_params(n_params: usize) -> Self {
        Self::new(n_params, 1e-3, 0.9, 0.999, 1e-8)
    }
    /// Apply one Adam update step to `params` using `grads`.
    ///
    /// Updates `params` in-place and increments the step counter.
    pub fn step_update(&mut self, params: &mut [f64], grads: &[f64]) {
        assert_eq!(
            params.len(),
            self.m.len(),
            "AdamOptimizer::step_update: params/m length mismatch"
        );
        assert_eq!(
            grads.len(),
            self.m.len(),
            "AdamOptimizer::step_update: grads/m length mismatch"
        );
        self.step += 1;
        let t = self.step as f64;
        let bias_corr1 = 1.0 - self.beta1.powf(t);
        let bias_corr2 = 1.0 - self.beta2.powf(t);
        for i in 0..params.len() {
            self.m[i] = self.beta1 * self.m[i] + (1.0 - self.beta1) * grads[i];
            self.v[i] = self.beta2 * self.v[i] + (1.0 - self.beta2) * grads[i] * grads[i];
            let m_hat = self.m[i] / bias_corr1;
            let v_hat = self.v[i] / bias_corr2;
            params[i] -= self.lr * m_hat / (v_hat.sqrt() + self.epsilon);
        }
    }
    /// Reset state (moments and step counter) to zero.
    pub fn reset(&mut self) {
        self.m.iter_mut().for_each(|x| *x = 0.0);
        self.v.iter_mut().for_each(|x| *x = 0.0);
        self.step = 0;
    }
}
/// A single graph neural network layer implementing the sum-aggregation
/// message passing update:
///
/// h_i^(l+1) = σ(W_self * h_i^(l) + W_neigh * Σ_{j ∈ N(i)} h_j^(l) + b)
///
/// All nodes share the same weight matrices.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GnnLayer {
    /// Input feature dimension.
    pub in_dim: usize,
    /// Output feature dimension.
    pub out_dim: usize,
    /// Self-loop weight matrix W_self \[out_dim × in_dim\].
    pub w_self: Vec<f64>,
    /// Neighbour aggregation weight matrix W_neigh \[out_dim × in_dim\].
    pub w_neigh: Vec<f64>,
    /// Bias vector \[out_dim\].
    pub bias: Vec<f64>,
    /// Activation function.
    pub activation: ExtActivation,
}
impl GnnLayer {
    /// Create a new GNN layer with zero-initialised weights.
    pub fn new(in_dim: usize, out_dim: usize, activation: ExtActivation) -> Self {
        Self {
            in_dim,
            out_dim,
            w_self: vec![0.0_f64; out_dim * in_dim],
            w_neigh: vec![0.0_f64; out_dim * in_dim],
            bias: vec![0.0_f64; out_dim],
            activation,
        }
    }
    /// Forward pass.
    ///
    /// * `node_feats` – `[n_nodes × in_dim]` flat row-major node feature matrix.
    /// * `adj`        – adjacency list: `adj[i]` contains the neighbour indices of node `i`.
    ///
    /// Returns `[n_nodes × out_dim]` flat row-major.
    pub fn forward(&self, node_feats: &[f64], n_nodes: usize, adj: &[Vec<usize>]) -> Vec<f64> {
        assert_eq!(node_feats.len(), n_nodes * self.in_dim);
        assert_eq!(adj.len(), n_nodes);
        let in_d = self.in_dim;
        let out_d = self.out_dim;
        let mut out = vec![0.0_f64; n_nodes * out_d];
        for i in 0..n_nodes {
            let h_self = &node_feats[i * in_d..(i + 1) * in_d];
            let mut agg = vec![0.0_f64; in_d];
            for &j in &adj[i] {
                let h_j = &node_feats[j * in_d..(j + 1) * in_d];
                for d in 0..in_d {
                    agg[d] += h_j[d];
                }
            }
            for o in 0..out_d {
                let mut acc = self.bias[o];
                for d in 0..in_d {
                    acc += self.w_self[o * in_d + d] * h_self[d];
                    acc += self.w_neigh[o * in_d + d] * agg[d];
                }
                out[i * out_d + o] = self.activation.apply(acc);
            }
        }
        out
    }
    /// Total number of trainable parameters.
    pub fn num_params(&self) -> usize {
        2 * self.out_dim * self.in_dim + self.out_dim
    }
}
/// A multi-layer message passing neural network stacking `GnnLayer`s.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MessagePassingNet {
    /// Ordered list of GNN layers.
    pub layers: Vec<GnnLayer>,
}
impl MessagePassingNet {
    /// Create an empty MPNN.
    pub fn new() -> Self {
        Self { layers: Vec::new() }
    }
    /// Add a GNN layer to the stack.
    pub fn add_layer(&mut self, layer: GnnLayer) {
        self.layers.push(layer);
    }
    /// Run all layers in sequence over a fixed graph.
    ///
    /// Returns the final node feature matrix `[n_nodes × last_out_dim]`.
    pub fn forward(&self, node_feats: &[f64], n_nodes: usize, adj: &[Vec<usize>]) -> Vec<f64> {
        let mut h = node_feats.to_vec();
        for layer in &self.layers {
            h = layer.forward(&h, n_nodes, adj);
        }
        h
    }
    /// Aggregate node features to a single graph-level representation (mean pooling).
    pub fn global_mean_pool(&self, node_feats: &[f64], n_nodes: usize, out_dim: usize) -> Vec<f64> {
        if n_nodes == 0 {
            return vec![0.0_f64; out_dim];
        }
        let mut pooled = vec![0.0_f64; out_dim];
        for i in 0..n_nodes {
            for d in 0..out_dim {
                pooled[d] += node_feats[i * out_dim + d];
            }
        }
        let inv_n = 1.0 / n_nodes as f64;
        for v in &mut pooled {
            *v *= inv_n;
        }
        pooled
    }
}
/// Accumulates gradients from multiple backward passes for mini-batch training.
#[derive(Debug, Clone)]
pub struct GradAccumulator {
    /// Accumulated weight gradients.
    pub grad_weights: Vec<f64>,
    /// Accumulated bias gradients.
    pub grad_biases: Vec<f64>,
    /// Number of samples accumulated.
    pub count: usize,
}
impl GradAccumulator {
    /// Create a new accumulator sized for `n_weights` weights and `n_biases` biases.
    pub fn new(n_weights: usize, n_biases: usize) -> Self {
        Self {
            grad_weights: vec![0.0; n_weights],
            grad_biases: vec![0.0; n_biases],
            count: 0,
        }
    }
    /// Add a set of gradients (accumulate without dividing).
    pub fn accumulate(&mut self, gw: &[f64], gb: &[f64]) {
        assert_eq!(gw.len(), self.grad_weights.len());
        assert_eq!(gb.len(), self.grad_biases.len());
        for (acc, &g) in self.grad_weights.iter_mut().zip(gw.iter()) {
            *acc += g;
        }
        for (acc, &g) in self.grad_biases.iter_mut().zip(gb.iter()) {
            *acc += g;
        }
        self.count += 1;
    }
    /// Compute mean gradients (divide by count) and return them.
    pub fn mean_grads(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.count.max(1) as f64;
        let gw: Vec<f64> = self.grad_weights.iter().map(|&g| g / n).collect();
        let gb: Vec<f64> = self.grad_biases.iter().map(|&g| g / n).collect();
        (gw, gb)
    }
    /// Zero all accumulated gradients and reset count.
    pub fn zero(&mut self) {
        self.grad_weights.iter_mut().for_each(|x| *x = 0.0);
        self.grad_biases.iter_mut().for_each(|x| *x = 0.0);
        self.count = 0;
    }
}
/// Multi-head attention module.
///
/// Projects Q, K, V with learned linear projections, runs `n_heads`
/// parallel attention heads, then concatenates and projects the output.
///
/// All weight matrices are stored flat row-major.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MultiHeadAttention {
    /// Model dimensionality.
    pub d_model: usize,
    /// Number of attention heads.
    pub n_heads: usize,
    /// Dimensionality per head: `d_model / n_heads`.
    pub d_head: usize,
    /// W_Q projection \[d_model × d_model\].
    pub w_q: Vec<f64>,
    /// W_K projection \[d_model × d_model\].
    pub w_k: Vec<f64>,
    /// W_V projection \[d_model × d_model\].
    pub w_v: Vec<f64>,
    /// W_O output projection \[d_model × d_model\].
    pub w_o: Vec<f64>,
    /// Output bias \[d_model\].
    pub b_o: Vec<f64>,
}
impl MultiHeadAttention {
    /// Create a new MHA module with zero-initialised projections.
    pub fn new(d_model: usize, n_heads: usize) -> Self {
        assert_eq!(d_model % n_heads, 0, "d_model must be divisible by n_heads");
        let d_head = d_model / n_heads;
        let dm2 = d_model * d_model;
        Self {
            d_model,
            n_heads,
            d_head,
            w_q: vec![0.0_f64; dm2],
            w_k: vec![0.0_f64; dm2],
            w_v: vec![0.0_f64; dm2],
            w_o: vec![0.0_f64; dm2],
            b_o: vec![0.0_f64; d_model],
        }
    }
    /// Initialise W_Q, W_K, W_V, W_O with identity-like weights for testing.
    pub fn init_identity(&mut self) {
        let dm = self.d_model;
        for row in 0..dm {
            self.w_q[row * dm + row] = 1.0;
            self.w_k[row * dm + row] = 1.0;
            self.w_v[row * dm + row] = 1.0;
            self.w_o[row * dm + row] = 1.0;
        }
    }
    /// Linear projection: `output = input @ W^T`  where W is `[out × in]`.
    fn project(
        input: &[f64],
        w: &[f64],
        seq_len: usize,
        in_dim: usize,
        out_dim: usize,
    ) -> Vec<f64> {
        let mut out = vec![0.0_f64; seq_len * out_dim];
        for t in 0..seq_len {
            for o in 0..out_dim {
                let mut acc = 0.0_f64;
                for i in 0..in_dim {
                    acc += input[t * in_dim + i] * w[o * in_dim + i];
                }
                out[t * out_dim + o] = acc;
            }
        }
        out
    }
    /// Forward pass.
    ///
    /// `x` has shape `[seq_len × d_model]` (flat row-major).
    /// Returns output of shape `[seq_len × d_model]`.
    pub fn forward(&self, x: &[f64], seq_len: usize) -> Vec<f64> {
        let dm = self.d_model;
        let dh = self.d_head;
        let nh = self.n_heads;
        let q_full = Self::project(x, &self.w_q, seq_len, dm, dm);
        let k_full = Self::project(x, &self.w_k, seq_len, dm, dm);
        let v_full = Self::project(x, &self.w_v, seq_len, dm, dm);
        let mut concat = vec![0.0_f64; seq_len * dm];
        for h in 0..nh {
            let mut q_h = vec![0.0_f64; seq_len * dh];
            let mut k_h = vec![0.0_f64; seq_len * dh];
            let mut v_h = vec![0.0_f64; seq_len * dh];
            for t in 0..seq_len {
                for d in 0..dh {
                    q_h[t * dh + d] = q_full[t * dm + h * dh + d];
                    k_h[t * dh + d] = k_full[t * dm + h * dh + d];
                    v_h[t * dh + d] = v_full[t * dm + h * dh + d];
                }
            }
            let head_out =
                scaled_dot_product_attention(&q_h, &k_h, &v_h, seq_len, seq_len, dh, dh, None);
            for t in 0..seq_len {
                for d in 0..dh {
                    concat[t * dm + h * dh + d] = head_out[t * dh + d];
                }
            }
        }
        let projected = Self::project(&concat, &self.w_o, seq_len, dm, dm);
        let mut output = projected;
        for t in 0..seq_len {
            for d in 0..dm {
                output[t * dm + d] += self.b_o[d];
            }
        }
        output
    }
    /// Total number of trainable parameters.
    pub fn num_params(&self) -> usize {
        4 * self.d_model * self.d_model + self.d_model
    }
}
/// A single fully-connected layer with f64 weights.
#[derive(Debug, Clone)]
pub struct NeuralLayer {
    /// Weight matrix: `weights[out][in]`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector of length `out_features`.
    pub biases: Vec<f64>,
    /// Activation function applied after the affine transform.
    pub activation: ActivationFn64,
}
impl NeuralLayer {
    /// Create a new layer with Xavier-uniform initialised weights.
    ///
    /// Xavier uniform: U(-limit, limit) where limit = sqrt(6 / (fan_in + fan_out)).
    pub fn new_xavier(in_features: usize, out_features: usize, activation: ActivationFn64) -> Self {
        let limit = (6.0_f64 / (in_features + out_features) as f64).sqrt();
        let mut state: u64 = 0x123456789abcdef0;
        let lcg_next = |s: &mut u64| -> f64 {
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let bits = (*s >> 33) as f64;
            bits / (u64::MAX as f64) * 2.0 * limit - limit
        };
        let weights: Vec<Vec<f64>> = (0..out_features)
            .map(|_| (0..in_features).map(|_| lcg_next(&mut state)).collect())
            .collect();
        let biases = vec![0.0_f64; out_features];
        Self {
            weights,
            biases,
            activation,
        }
    }
    /// Forward pass: activation(W * input + b).
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let out_features = self.weights.len();
        let mut output = Vec::with_capacity(out_features);
        for o in 0..out_features {
            let mut acc = self.biases[o];
            for (i, &x) in input.iter().enumerate() {
                acc += self.weights[o][i] * x;
            }
            output.push(self.activation.apply(acc));
        }
        output
    }
}
/// Attention-based graph readout that computes a weighted sum of node features.
///
/// For each node i, computes a scalar attention score a_i = sigmoid(w · h_i + b),
/// then returns Σ_i a_i * h_i.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AttentionReadout {
    /// Feature dimensionality.
    pub d_feat: usize,
    /// Attention weight vector \[d_feat\].
    pub w_attn: Vec<f64>,
    /// Attention bias (scalar).
    pub b_attn: f64,
}
impl AttentionReadout {
    /// Create with zero attention weights.
    pub fn new(d_feat: usize) -> Self {
        Self {
            d_feat,
            w_attn: vec![0.0_f64; d_feat],
            b_attn: 0.0,
        }
    }
    /// Compute attention-weighted sum over nodes.
    ///
    /// `node_feats`: `[n_nodes × d_feat]` flat row-major.
    pub fn forward(&self, node_feats: &[f64], n_nodes: usize) -> Vec<f64> {
        let df = self.d_feat;
        let mut out = vec![0.0_f64; df];
        let mut attn_scores = Vec::with_capacity(n_nodes);
        for i in 0..n_nodes {
            let h = &node_feats[i * df..(i + 1) * df];
            let raw: f64 = h
                .iter()
                .zip(self.w_attn.iter())
                .map(|(&x, &w)| x * w)
                .sum::<f64>()
                + self.b_attn;
            let score = 1.0 / (1.0 + (-raw).exp());
            attn_scores.push(score);
        }
        for i in 0..n_nodes {
            let h = &node_feats[i * df..(i + 1) * df];
            for d in 0..df {
                out[d] += attn_scores[i] * h[d];
            }
        }
        out
    }
}
/// A single transformer encoder block:
/// x → MHA(LayerNorm(x)) + x → FFN(LayerNorm(·)) + ·
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TransformerBlock {
    /// Multi-head self-attention module.
    pub mha: MultiHeadAttention,
    /// Feed-forward network.
    pub ffn: TransformerFfn,
    /// Layer norm before MHA.
    pub ln1: LayerNorm,
    /// Layer norm before FFN.
    pub ln2: LayerNorm,
    /// Model dimensionality.
    pub d_model: usize,
}
impl TransformerBlock {
    /// Create a new transformer block with zero weights.
    pub fn new(d_model: usize, n_heads: usize, d_ff: usize) -> Self {
        Self {
            mha: MultiHeadAttention::new(d_model, n_heads),
            ffn: TransformerFfn::new(d_model, d_ff),
            ln1: LayerNorm::new(d_model),
            ln2: LayerNorm::new(d_model),
            d_model,
        }
    }
    /// Forward pass with pre-norm style residual connections.
    ///
    /// `x` is flat row-major `[seq_len × d_model]`.
    pub fn forward(&self, x: &[f64], seq_len: usize) -> Vec<f64> {
        let dm = self.d_model;
        let mut normed1 = vec![0.0_f64; seq_len * dm];
        for t in 0..seq_len {
            let row = &x[t * dm..(t + 1) * dm];
            let n = self.ln1.forward(row);
            normed1[t * dm..(t + 1) * dm].copy_from_slice(&n);
        }
        let attn_out = self.mha.forward(&normed1, seq_len);
        let mut x1 = vec![0.0_f64; seq_len * dm];
        for i in 0..x1.len() {
            x1[i] = x[i] + attn_out[i];
        }
        let mut normed2 = vec![0.0_f64; seq_len * dm];
        for t in 0..seq_len {
            let row = &x1[t * dm..(t + 1) * dm];
            let n = self.ln2.forward(row);
            normed2[t * dm..(t + 1) * dm].copy_from_slice(&n);
        }
        let ffn_out = self.ffn.forward(&normed2, seq_len);
        let mut x2 = vec![0.0_f64; seq_len * dm];
        for i in 0..x2.len() {
            x2[i] = x1[i] + ffn_out[i];
        }
        x2
    }
}
/// Activation functions for neural network layers.
#[derive(Debug, Clone, PartialEq)]
pub enum ActivationFn {
    /// Hyperbolic tangent.
    Tanh,
    /// Rectified linear unit.
    Relu,
    /// Logistic sigmoid.
    Sigmoid,
    /// Sigmoid-weighted linear unit.
    Silu,
    /// Gaussian error linear unit (approximation).
    Gelu,
    /// Identity / no activation.
    Linear,
}
impl ActivationFn {
    /// Evaluate the activation function at `x`.
    pub fn apply(&self, x: f32) -> f32 {
        match self {
            ActivationFn::Tanh => x.tanh(),
            ActivationFn::Relu => x.max(0.0),
            ActivationFn::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            ActivationFn::Silu => x / (1.0 + (-x).exp()),
            ActivationFn::Gelu => {
                let cdf = 0.5
                    * (1.0
                        + (std::f32::consts::FRAC_2_SQRT_PI.sqrt() * (x + 0.044715 * x * x * x))
                            .tanh());
                x * cdf
            }
            ActivationFn::Linear => x,
        }
    }
    /// Evaluate the derivative of the activation function at `x`.
    pub fn derivative(&self, x: f32) -> f32 {
        match self {
            ActivationFn::Tanh => {
                let t = x.tanh();
                1.0 - t * t
            }
            ActivationFn::Relu => {
                if x > 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
            ActivationFn::Sigmoid => {
                let s = 1.0 / (1.0 + (-x).exp());
                s * (1.0 - s)
            }
            ActivationFn::Silu => {
                let s = 1.0 / (1.0 + (-x).exp());
                s + x * s * (1.0 - s)
            }
            ActivationFn::Gelu => {
                let eps = 1e-5_f32;
                (self.apply(x + eps) - self.apply(x - eps)) / (2.0 * eps)
            }
            ActivationFn::Linear => 1.0,
        }
    }
}
/// Layer normalisation (Ba et al., 2016).
///
/// Normalises the *entire* feature vector of a single sample to zero mean and
/// unit variance, then applies learned scale (gamma) and shift (beta).
#[derive(Debug, Clone)]
pub struct LayerNormLayer {
    /// Number of features (last dimension size).
    pub n_features: usize,
    /// Learned scale parameter (gamma), initialised to 1.
    pub gamma: Vec<f64>,
    /// Learned shift parameter (beta), initialised to 0.
    pub beta: Vec<f64>,
    /// Numerical stability constant.
    pub epsilon: f64,
}
impl LayerNormLayer {
    /// Create a new LayerNorm with identity transform (gamma=1, beta=0).
    pub fn new(n_features: usize) -> Self {
        Self {
            n_features,
            gamma: vec![1.0; n_features],
            beta: vec![0.0; n_features],
            epsilon: 1e-5,
        }
    }
    /// Apply layer normalisation to one sample vector.
    ///
    /// output\[i\] = gamma\[i\] * (input\[i\] - mean) / sqrt(var + eps) + beta\[i\]
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        assert_eq!(
            input.len(),
            self.n_features,
            "LayerNorm: input size mismatch"
        );
        let n = self.n_features as f64;
        let mean: f64 = input.iter().sum::<f64>() / n;
        let var: f64 = input.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / n;
        let std_inv = 1.0 / (var + self.epsilon).sqrt();
        (0..self.n_features)
            .map(|i| self.gamma[i] * (input[i] - mean) * std_inv + self.beta[i])
            .collect()
    }
    /// Compute gradient of the layer norm output with respect to the input.
    ///
    /// Returns `(d_input, d_gamma, d_beta)` given upstream gradient `d_output`.
    #[allow(non_snake_case)]
    pub fn backward(&self, input: &[f64], d_output: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        assert_eq!(input.len(), self.n_features);
        assert_eq!(d_output.len(), self.n_features);
        let n = self.n_features as f64;
        let mean: f64 = input.iter().sum::<f64>() / n;
        let var: f64 = input.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / n;
        let std_inv = 1.0 / (var + self.epsilon).sqrt();
        let x_hat: Vec<f64> = input.iter().map(|&x| (x - mean) * std_inv).collect();
        let d_gamma: Vec<f64> = (0..self.n_features)
            .map(|i| d_output[i] * x_hat[i])
            .collect();
        let d_beta: Vec<f64> = d_output.to_vec();
        let d_x_hat: Vec<f64> = (0..self.n_features)
            .map(|i| d_output[i] * self.gamma[i])
            .collect();
        let sum_d_x_hat: f64 = d_x_hat.iter().sum();
        let sum_d_x_hat_xhat: f64 = d_x_hat.iter().zip(x_hat.iter()).map(|(&a, &b)| a * b).sum();
        let d_input: Vec<f64> = (0..self.n_features)
            .map(|i| std_inv * (d_x_hat[i] - (sum_d_x_hat + x_hat[i] * sum_d_x_hat_xhat) / n))
            .collect();
        (d_input, d_gamma, d_beta)
    }
}
/// Feature-wise Z-score normalizer.
///
/// Stores per-feature mean and standard deviation fitted on a training corpus.
#[derive(Debug, Clone)]
pub struct DataNormalizer {
    /// Per-feature mean.
    pub mean: Vec<f32>,
    /// Per-feature standard deviation.
    pub std_dev: Vec<f32>,
}
impl DataNormalizer {
    /// Fit normalizer statistics from a collection of sample vectors.
    ///
    /// # Panics
    /// Panics if `data` is empty or if sample vectors have inconsistent lengths.
    pub fn fit(data: &[Vec<f32>]) -> Self {
        assert!(
            !data.is_empty(),
            "DataNormalizer::fit: data must be non-empty"
        );
        let n_features = data[0].len();
        let n = data.len() as f32;
        let mut mean = vec![0.0_f32; n_features];
        for sample in data {
            assert_eq!(
                sample.len(),
                n_features,
                "DataNormalizer::fit: inconsistent sample length"
            );
            for (k, &v) in sample.iter().enumerate() {
                mean[k] += v;
            }
        }
        for m in &mut mean {
            *m /= n;
        }
        let mut variance = vec![0.0_f32; n_features];
        for sample in data {
            for (k, &v) in sample.iter().enumerate() {
                let diff = v - mean[k];
                variance[k] += diff * diff;
            }
        }
        let std_dev: Vec<f32> = variance
            .iter()
            .map(|&v| {
                let s = (v / n).sqrt();
                if s < 1e-8 { 1.0 } else { s }
            })
            .collect();
        DataNormalizer { mean, std_dev }
    }
    /// Standardise a single sample: `(x - mean) / std`.
    pub fn transform(&self, x: &[f32]) -> Vec<f32> {
        x.iter()
            .zip(self.mean.iter())
            .zip(self.std_dev.iter())
            .map(|((&xi, &m), &s)| (xi - m) / s)
            .collect()
    }
    /// Invert standardisation: `x * std + mean`.
    pub fn inverse_transform(&self, x: &[f32]) -> Vec<f32> {
        x.iter()
            .zip(self.mean.iter())
            .zip(self.std_dev.iter())
            .map(|((&xi, &m), &s)| xi * s + m)
            .collect()
    }
}
/// Activation function for f64-precision neural network layers.
#[derive(Debug, Clone, PartialEq)]
pub enum ActivationFn64 {
    /// Rectified linear unit: max(0, x).
    Relu,
    /// Logistic sigmoid: 1 / (1 + exp(-x)).
    Sigmoid,
    /// Hyperbolic tangent.
    Tanh,
    /// Identity (no activation).
    Linear,
}
impl ActivationFn64 {
    /// Evaluate the activation at a single value.
    pub fn apply(&self, x: f64) -> f64 {
        match self {
            ActivationFn64::Relu => x.max(0.0),
            ActivationFn64::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            ActivationFn64::Tanh => x.tanh(),
            ActivationFn64::Linear => x,
        }
    }
    /// Apply the activation in-place to every element of a vector.
    pub fn apply_batch(&self, v: &mut Vec<f64>) {
        for x in v.iter_mut() {
            *x = self.apply(*x);
        }
    }
}
/// A sequential feed-forward neural network using f64 precision.
#[derive(Debug, Clone)]
pub struct NeuralNetwork {
    /// Ordered list of layers.
    pub layers: Vec<NeuralLayer>,
}
impl NeuralNetwork {
    /// Build a network with Xavier-initialised weights from a list of layer sizes.
    ///
    /// All hidden layers use `activation`; the final layer uses `ActivationFn64::Linear`.
    pub fn new(layer_sizes: &[usize], activation: ActivationFn64) -> Self {
        assert!(
            layer_sizes.len() >= 2,
            "need at least input and output size"
        );
        let mut layers = Vec::new();
        for i in 0..layer_sizes.len() - 1 {
            let act = if i == layer_sizes.len() - 2 {
                ActivationFn64::Linear
            } else {
                activation.clone()
            };
            layers.push(NeuralLayer::new_xavier(
                layer_sizes[i],
                layer_sizes[i + 1],
                act,
            ));
        }
        Self { layers }
    }
    /// Run a forward pass through all layers.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let mut current: Vec<f64> = input.to_vec();
        for layer in &self.layers {
            current = layer.forward(&current);
        }
        current
    }
    /// Expected input dimension (size of the first layer's input).
    pub fn input_dim(&self) -> usize {
        self.layers
            .first()
            .map_or(0, |l| l.weights.first().map_or(0, |r| r.len()))
    }
    /// Expected output dimension (size of the last layer's output).
    pub fn output_dim(&self) -> usize {
        self.layers.last().map_or(0, |l| l.biases.len())
    }
}
/// A sequential feed-forward neural network.
#[derive(Debug, Clone)]
pub struct FeedForwardNet {
    /// Ordered list of dense layers.
    pub layers: Vec<DenseLayer>,
}
impl FeedForwardNet {
    /// Create an empty network.
    pub fn new() -> Self {
        FeedForwardNet { layers: Vec::new() }
    }
    /// Append a layer to the end of the network.
    pub fn add_layer(&mut self, layer: DenseLayer) {
        self.layers.push(layer);
    }
    /// Run a forward pass through all layers.
    pub fn forward(&self, input: &[f32]) -> Vec<f32> {
        let mut current: Vec<f32> = input.to_vec();
        for layer in &self.layers {
            current = layer.forward(&current);
        }
        current
    }
    /// Returns the expected input width (`None` if the network has no layers).
    pub fn input_size(&self) -> Option<usize> {
        self.layers.first().map(|l| l.in_features)
    }
    /// Returns the output width of the last layer (`None` if empty).
    pub fn output_size(&self) -> Option<usize> {
        self.layers.last().map(|l| l.out_features)
    }
    /// Sum of parameters across all layers.
    pub fn total_parameters(&self) -> usize {
        self.layers.iter().map(|l| l.parameter_count()).sum()
    }
}
impl FeedForwardNet {
    /// Compute the total gradient norm across all layers.
    ///
    /// `layer_grads[i]` is the concatenated `[grad_weights, grad_biases]` for
    /// layer `i`.  Returns the L2 norm of all gradients combined.
    pub fn compute_gradient_norm(&self, layer_grads: &[Vec<f32>]) -> f32 {
        let sum_sq: f32 = layer_grads
            .iter()
            .flat_map(|g| g.iter())
            .map(|&v| v * v)
            .sum();
        sum_sq.sqrt()
    }
    /// Clip per-layer gradient vectors in-place so their combined norm ≤ `max_norm`.
    /// Returns the pre-clip norm.
    pub fn clip_gradients(&self, layer_grads: &mut Vec<Vec<f32>>, max_norm: f32) -> f32 {
        let norm = self.compute_gradient_norm(layer_grads);
        if norm > max_norm && norm > 0.0 {
            let scale = max_norm / norm;
            for g in layer_grads.iter_mut() {
                for v in g.iter_mut() {
                    *v *= scale;
                }
            }
        }
        norm
    }
}
/// Layer normalisation applied to each time step independently.
///
/// Normalises a feature vector of length `n_features` to zero mean and unit
/// variance, then applies learnable scale (gamma) and bias (beta).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LayerNorm {
    /// Number of features.
    pub n_features: usize,
    /// Learnable scale parameter.
    pub gamma: Vec<f64>,
    /// Learnable shift parameter.
    pub beta: Vec<f64>,
    /// Numerical stability constant.
    pub epsilon: f64,
}
impl LayerNorm {
    /// Create a new layer norm with identity initialisation (gamma=1, beta=0).
    pub fn new(n_features: usize) -> Self {
        Self {
            n_features,
            gamma: vec![1.0_f64; n_features],
            beta: vec![0.0_f64; n_features],
            epsilon: 1e-5,
        }
    }
    /// Normalise a single feature vector.
    pub fn forward(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.n_features);
        let n = self.n_features as f64;
        let mean = x.iter().sum::<f64>() / n;
        let var = x.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n;
        let std = (var + self.epsilon).sqrt();
        x.iter()
            .enumerate()
            .map(|(i, &v)| self.gamma[i] * (v - mean) / std + self.beta[i])
            .collect()
    }
}
/// A mock GPU buffer for batched neural network inference.
///
/// In a real GPU backend this would be a device-side buffer; here we
/// store a flat f64 array in host memory.
#[derive(Debug, Clone)]
pub struct GpuNeuralBuffer {
    /// Batch size (number of samples).
    pub batch_size: usize,
    /// Dimensionality of each input sample.
    pub input_dim: usize,
    /// Dimensionality of each output sample.
    pub output_dim: usize,
    /// Flat data storage: `batch_size * max(input_dim, output_dim)` elements.
    pub data: Vec<f64>,
}
impl GpuNeuralBuffer {
    /// Pack a slice of 3-D positions into a buffer suitable for network input.
    ///
    /// Each position `[x, y, z]` becomes three consecutive f64 values.
    pub fn pack_positions(positions: &[[f64; 3]]) -> Self {
        let batch_size = positions.len();
        let input_dim = 3;
        let output_dim = 3;
        let mut data = Vec::with_capacity(batch_size * input_dim);
        for p in positions {
            data.push(p[0]);
            data.push(p[1]);
            data.push(p[2]);
        }
        Self {
            batch_size,
            input_dim,
            output_dim,
            data,
        }
    }
    /// Unpack the buffer contents as a list of 3-D force vectors.
    ///
    /// Assumes `data` has been filled with `batch_size * 3` values by
    /// a prior inference step.
    pub fn unpack_forces(&self) -> Vec<[f64; 3]> {
        self.data.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
    }
}
/// Sinusoidal positional encoding (Vaswani et al., 2017).
///
/// For each position `pos` and dimension `i` in a `d_model`-dimensional
/// embedding space:
///   PE\[pos, 2i\]   = sin(pos / 10000^(2i/d_model))
///   PE\[pos, 2i+1\] = cos(pos / 10000^(2i/d_model))
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PositionalEncoding {
    /// Embedding dimensionality.
    pub d_model: usize,
    /// Maximum sequence length supported.
    pub max_len: usize,
    /// Pre-computed encoding table: `table[pos][dim]`.
    pub table: Vec<Vec<f64>>,
}
impl PositionalEncoding {
    /// Build the positional encoding table up to `max_len` positions.
    pub fn new(d_model: usize, max_len: usize) -> Self {
        let mut table = vec![vec![0.0_f64; d_model]; max_len];
        for pos in 0..max_len {
            for i in 0..(d_model / 2) {
                let angle = (pos as f64) / (10000.0_f64).powf(2.0 * i as f64 / d_model as f64);
                table[pos][2 * i] = angle.sin();
                if 2 * i + 1 < d_model {
                    table[pos][2 * i + 1] = angle.cos();
                }
            }
        }
        Self {
            d_model,
            max_len,
            table,
        }
    }
    /// Add positional encoding to a sequence of embeddings in-place.
    ///
    /// `embeddings[t]` is a feature vector of length `d_model`.
    pub fn add_to_sequence(&self, embeddings: &mut Vec<Vec<f64>>) {
        for (t, emb) in embeddings.iter_mut().enumerate() {
            if t >= self.max_len {
                break;
            }
            for d in 0..emb.len().min(self.d_model) {
                emb[d] += self.table[t][d];
            }
        }
    }
    /// Return the positional encoding vector for position `pos`.
    pub fn get(&self, pos: usize) -> &[f64] {
        &self.table[pos.min(self.max_len - 1)]
    }
}
/// A fully-connected layer with f64 weights supporting forward pass and
/// gradient computation for backpropagation.
#[derive(Debug, Clone)]
pub struct DenseLayer64 {
    /// Weight matrix in row-major layout: `weights[out * in_features + in]`.
    pub weights: Vec<f64>,
    /// Bias vector of length `out_features`.
    pub biases: Vec<f64>,
    /// Number of input features.
    pub in_features: usize,
    /// Number of output features.
    pub out_features: usize,
    /// Activation function.
    pub activation: ExtActivation,
    /// Pre-activation outputs from the last forward pass (z = W*x + b).
    pub last_pre_act: Vec<f64>,
    /// Post-activation outputs from the last forward pass.
    pub last_output: Vec<f64>,
    /// Last input fed to this layer.
    pub last_input: Vec<f64>,
}
impl DenseLayer64 {
    /// Create a new layer with zero-initialised weights and biases.
    pub fn new(in_features: usize, out_features: usize, activation: ExtActivation) -> Self {
        Self {
            weights: vec![0.0_f64; out_features * in_features],
            biases: vec![0.0_f64; out_features],
            in_features,
            out_features,
            activation,
            last_pre_act: Vec::new(),
            last_output: Vec::new(),
            last_input: Vec::new(),
        }
    }
    /// Forward pass: computes `activation(W * input + b)`.
    /// Caches `pre_act`, `output`, and `input` for backprop.
    pub fn forward(&mut self, input: &[f64]) -> Vec<f64> {
        assert_eq!(
            input.len(),
            self.in_features,
            "DenseLayer64::forward: input size mismatch"
        );
        self.last_input = input.to_vec();
        let mut pre_act = Vec::with_capacity(self.out_features);
        for o in 0..self.out_features {
            let row = o * self.in_features;
            let mut acc = self.biases[o];
            for i in 0..self.in_features {
                acc += self.weights[row + i] * input[i];
            }
            pre_act.push(acc);
        }
        let output: Vec<f64> = pre_act.iter().map(|&z| self.activation.apply(z)).collect();
        self.last_pre_act = pre_act;
        self.last_output = output.clone();
        output
    }
    /// Backward pass: computes gradients w.r.t. weights, biases, and input.
    ///
    /// `delta_out` is the gradient of the loss w.r.t. this layer's output
    /// (same shape as `last_output`).
    ///
    /// Returns `(grad_weights, grad_biases, delta_in)` where `delta_in` is the
    /// gradient passed to the previous layer.
    #[allow(clippy::too_many_arguments)]
    pub fn backward(&self, delta_out: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        assert_eq!(
            delta_out.len(),
            self.out_features,
            "DenseLayer64::backward: delta_out size mismatch"
        );
        let delta_pre: Vec<f64> = delta_out
            .iter()
            .zip(self.last_pre_act.iter())
            .map(|(&d, &z)| d * self.activation.derivative(z))
            .collect();
        let mut grad_weights = vec![0.0_f64; self.out_features * self.in_features];
        for o in 0..self.out_features {
            let row = o * self.in_features;
            for i in 0..self.in_features {
                grad_weights[row + i] = delta_pre[o] * self.last_input[i];
            }
        }
        let grad_biases = delta_pre.clone();
        let mut delta_in = vec![0.0_f64; self.in_features];
        for o in 0..self.out_features {
            let row = o * self.in_features;
            for i in 0..self.in_features {
                delta_in[i] += self.weights[row + i] * delta_pre[o];
            }
        }
        (grad_weights, grad_biases, delta_in)
    }
    /// Apply gradient updates using a simple SGD step.
    pub fn apply_sgd(&mut self, grad_weights: &[f64], grad_biases: &[f64], lr: f64) {
        for (w, &gw) in self.weights.iter_mut().zip(grad_weights.iter()) {
            *w -= lr * gw;
        }
        for (b, &gb) in self.biases.iter_mut().zip(grad_biases.iter()) {
            *b -= lr * gb;
        }
    }
    /// Total number of parameters.
    pub fn num_params(&self) -> usize {
        self.out_features * self.in_features + self.out_features
    }
}
/// Atomic neural network potential (NNP) with one sub-network per element.
///
/// Architecture follows Behler (2011): each atom contributes an atomic energy
/// predicted by an element-specific feed-forward network whose input is the
/// Behler-Parrinello descriptor vector.
#[derive(Debug)]
pub struct AtomicNeuralNetwork {
    /// Element-specific networks keyed by atomic number.
    pub networks: HashMap<u8, FeedForwardNet>,
    /// Symmetry-function descriptor shared by all elements.
    pub descriptor: BehlerParrinelloDescriptor,
}
impl AtomicNeuralNetwork {
    /// Create a new AANN with the given descriptor.
    pub fn new(descriptor: BehlerParrinelloDescriptor) -> Self {
        AtomicNeuralNetwork {
            networks: HashMap::new(),
            descriptor,
        }
    }
    /// Register a sub-network for the given atomic number.
    pub fn add_element_network(&mut self, atomic_number: u8, net: FeedForwardNet) {
        self.networks.insert(atomic_number, net);
    }
    /// Predict the atomic energy for one atom given its descriptor.
    ///
    /// Returns `None` if no network is registered for `atomic_number`.
    pub fn atomic_energy(&self, atomic_number: u8, descriptor: &[f32]) -> Option<f32> {
        self.networks
            .get(&atomic_number)
            .map(|net| net.forward(descriptor)[0])
    }
    /// Sum of atomic energies over all atoms.
    ///
    /// Atoms whose element has no registered network contribute 0.
    pub fn total_energy(&self, positions: &[[f64; 3]], atomic_numbers: &[u8]) -> f64 {
        assert_eq!(
            positions.len(),
            atomic_numbers.len(),
            "total_energy: positions and atomic_numbers must have the same length"
        );
        let mut e_total = 0.0_f64;
        for (i, &z) in atomic_numbers.iter().enumerate() {
            let desc_f64 = self.descriptor.descriptor_vector(positions, i);
            let desc_f32: Vec<f32> = desc_f64.iter().map(|&v| v as f32).collect();
            if let Some(e) = self.atomic_energy(z, &desc_f32) {
                e_total += e as f64;
            }
        }
        e_total
    }
}
/// A single fully-connected (dense) layer with an activation function.
///
/// Weights are stored in row-major order: `weights[out * in_features + in]`.
#[derive(Debug, Clone)]
pub struct DenseLayer {
    /// Weight matrix in row-major layout `[out_features × in_features]`.
    pub weights: Vec<f32>,
    /// Bias vector of length `out_features`.
    pub biases: Vec<f32>,
    /// Number of input features.
    pub in_features: usize,
    /// Number of output features.
    pub out_features: usize,
    /// Activation function applied after the affine transform.
    pub activation: ActivationFn,
}
impl DenseLayer {
    /// Create a new layer with zero-initialised weights and biases.
    pub fn new(in_features: usize, out_features: usize, activation: ActivationFn) -> Self {
        DenseLayer {
            weights: vec![0.0_f32; out_features * in_features],
            biases: vec![0.0_f32; out_features],
            in_features,
            out_features,
            activation,
        }
    }
    /// Compute `activation(W * input + b)`.
    ///
    /// # Panics
    /// Panics if `input.len() != self.in_features`.
    pub fn forward(&self, input: &[f32]) -> Vec<f32> {
        assert_eq!(
            input.len(),
            self.in_features,
            "DenseLayer::forward: input length {} != in_features {}",
            input.len(),
            self.in_features
        );
        let mut output = Vec::with_capacity(self.out_features);
        for o in 0..self.out_features {
            let row_offset = o * self.in_features;
            let mut acc = self.biases[o];
            for i in 0..self.in_features {
                acc += self.weights[row_offset + i] * input[i];
            }
            output.push(self.activation.apply(acc));
        }
        output
    }
    /// Replace the weight matrix (must have length `out_features * in_features`).
    ///
    /// # Panics
    /// Panics if `w.len()` does not match.
    pub fn set_weights(&mut self, w: &[f32]) {
        assert_eq!(
            w.len(),
            self.out_features * self.in_features,
            "set_weights: expected {} elements, got {}",
            self.out_features * self.in_features,
            w.len()
        );
        self.weights.copy_from_slice(w);
    }
    /// Replace the bias vector (must have length `out_features`).
    ///
    /// # Panics
    /// Panics if `b.len()` does not match.
    pub fn set_biases(&mut self, b: &[f32]) {
        assert_eq!(
            b.len(),
            self.out_features,
            "set_biases: expected {} elements, got {}",
            self.out_features,
            b.len()
        );
        self.biases.copy_from_slice(b);
    }
    /// Total number of trainable parameters (weights + biases).
    pub fn parameter_count(&self) -> usize {
        self.out_features * self.in_features + self.out_features
    }
}
/// Dropout regularisation layer.
///
/// During training, each neuron is set to zero with probability `rate`.
/// During inference (training=false) the layer passes inputs through unchanged
/// but scales by `1 - rate` to maintain expected magnitude.
#[derive(Debug, Clone)]
pub struct DropoutLayer {
    /// Probability of dropping a unit (0.0 = no dropout, 1.0 = all dropped).
    pub rate: f64,
    /// Whether the layer is in training mode.
    pub training: bool,
    /// Mask used in the last forward pass (1.0 = kept, 0.0 = dropped).
    pub last_mask: Vec<f64>,
    /// Seed for the deterministic LCG used to generate the mask.
    pub(super) seed: u64,
}
impl DropoutLayer {
    /// Create a new dropout layer.
    pub fn new(rate: f64, training: bool) -> Self {
        assert!(
            (0.0..=1.0).contains(&rate),
            "dropout rate must be in [0, 1]"
        );
        Self {
            rate,
            training,
            last_mask: Vec::new(),
            seed: 0xdeadbeefcafe1234,
        }
    }
    /// Set the seed for reproducible mask generation.
    pub fn set_seed(&mut self, seed: u64) {
        self.seed = seed;
    }
    /// Apply dropout to `input`.
    ///
    /// In training mode, randomly zeros elements with probability `rate`
    /// and scales the rest by `1 / (1 - rate)` (inverted dropout).
    /// In eval mode, passes through unchanged.
    pub fn forward(&mut self, input: &[f64]) -> Vec<f64> {
        if !self.training || self.rate == 0.0 {
            self.last_mask = vec![1.0; input.len()];
            return input.to_vec();
        }
        if self.rate == 1.0 {
            self.last_mask = vec![0.0; input.len()];
            return vec![0.0; input.len()];
        }
        let scale = 1.0 / (1.0 - self.rate);
        let mut mask = Vec::with_capacity(input.len());
        let mut output = Vec::with_capacity(input.len());
        for &x in input {
            self.seed = self
                .seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u = (self.seed >> 11) as f64 / (1u64 << 53) as f64;
            let m = if u >= self.rate { scale } else { 0.0 };
            mask.push(m);
            output.push(x * m);
        }
        self.last_mask = mask;
        output
    }
    /// Backward pass: applies the stored mask to the upstream gradient.
    pub fn backward(&self, delta_out: &[f64]) -> Vec<f64> {
        delta_out
            .iter()
            .zip(self.last_mask.iter())
            .map(|(&d, &m)| d * m)
            .collect()
    }
}
/// Behler-Parrinello symmetry functions for constructing atomic descriptors.
///
/// Reference: J. Behler and M. Parrinello, PRL 98, 146401 (2007).
#[derive(Debug, Clone)]
pub struct BehlerParrinelloDescriptor {
    /// Radial decay parameters η for G2 functions.
    pub eta: Vec<f64>,
    /// Shift parameters R_s for G2 functions.
    pub rs: Vec<f64>,
    /// Cutoff radius R_c in Ångström (or same units as positions).
    pub cutoff: f64,
}
impl BehlerParrinelloDescriptor {
    /// Smooth cutoff function.
    ///
    /// f_c(r) = 0.5 * (cos(π r / R_c) + 1)  for r < R_c, else 0.
    pub fn cutoff_fn(r: f64, rc: f64) -> f64 {
        if r < rc {
            0.5 * ((PI_F64 * r / rc).cos() + 1.0)
        } else {
            0.0
        }
    }
    /// G1 radial symmetry function: G1 = f_c(r).
    pub fn radial_g1(r: f64, rc: f64) -> f64 {
        Self::cutoff_fn(r, rc)
    }
    /// G2 radial symmetry function: G2 = exp(-η (r - R_s)²) * f_c(r).
    pub fn radial_g2(r: f64, eta: f64, rs: f64, rc: f64) -> f64 {
        let diff = r - rs;
        (-eta * diff * diff).exp() * Self::cutoff_fn(r, rc)
    }
    /// G4 angular symmetry function (two-body factor for a triplet i-j-k).
    ///
    /// G4 = 2^(1-ζ) * (1 + λ cos θ)^ζ * exp(-η (r_ij² + r_ik² + r_jk²)) * f_c(r_ij) f_c(r_ik) f_c(r_jk)
    #[allow(clippy::too_many_arguments)]
    pub fn angular_g4(
        r_ij: f64,
        r_ik: f64,
        r_jk: f64,
        cos_theta: f64,
        eta: f64,
        zeta: f64,
        lambda: f64,
        rc: f64,
    ) -> f64 {
        let angular = (1.0 + lambda * cos_theta).powf(zeta);
        let radial = (-eta * (r_ij * r_ij + r_ik * r_ik + r_jk * r_jk)).exp();
        let fc = Self::cutoff_fn(r_ij, rc) * Self::cutoff_fn(r_ik, rc) * Self::cutoff_fn(r_jk, rc);
        2.0_f64.powf(1.0 - zeta) * angular * radial * fc
    }
    /// Compute a single G2 descriptor value (convenience wrapper).
    pub fn compute(r_ij: f64, eta: f64, rs: f64, cutoff: f64) -> f64 {
        Self::radial_g2(r_ij, eta, rs, cutoff)
    }
    /// Build a full G2 descriptor vector for atom `center_idx`.
    ///
    /// For every (η, R_s) pair the function sums G2(r_ij, η, R_s, R_c) over all
    /// neighbours j ≠ center_idx that lie within the cutoff radius.
    pub fn descriptor_vector(&self, positions: &[[f64; 3]], center_idx: usize) -> Vec<f64> {
        let n_descriptors = self.eta.len();
        let mut desc = vec![0.0_f64; n_descriptors];
        let ci = positions[center_idx];
        for (j, pos_j) in positions.iter().enumerate() {
            if j == center_idx {
                continue;
            }
            let dx = pos_j[0] - ci[0];
            let dy = pos_j[1] - ci[1];
            let dz = pos_j[2] - ci[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r >= self.cutoff {
                continue;
            }
            for k in 0..n_descriptors {
                desc[k] += Self::radial_g2(r, self.eta[k], self.rs[k], self.cutoff);
            }
        }
        desc
    }
}
/// Extended activation functions with additional variants for f64 paths.
#[derive(Debug, Clone, PartialEq)]
pub enum ExtActivation {
    /// Leaky ReLU: max(alpha * x, x) where alpha is the negative slope.
    LeakyRelu(f64),
    /// Swish: x * sigmoid(beta * x).  beta=1 recovers SiLU.
    Swish(f64),
    /// Standard ReLU.
    Relu,
    /// Logistic sigmoid.
    Sigmoid,
    /// Hyperbolic tangent.
    Tanh,
    /// Identity.
    Linear,
}
impl ExtActivation {
    /// Evaluate the activation function at `x`.
    pub fn apply(&self, x: f64) -> f64 {
        match self {
            ExtActivation::LeakyRelu(alpha) => {
                if x >= 0.0 {
                    x
                } else {
                    alpha * x
                }
            }
            ExtActivation::Swish(beta) => x / (1.0 + (-beta * x).exp()),
            ExtActivation::Relu => x.max(0.0),
            ExtActivation::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            ExtActivation::Tanh => x.tanh(),
            ExtActivation::Linear => x,
        }
    }
    /// Evaluate the derivative of the activation function at `x`.
    pub fn derivative(&self, x: f64) -> f64 {
        match self {
            ExtActivation::LeakyRelu(alpha) => {
                if x >= 0.0 {
                    1.0
                } else {
                    *alpha
                }
            }
            ExtActivation::Swish(beta) => {
                let sig = 1.0 / (1.0 + (-beta * x).exp());
                sig + beta * x * sig * (1.0 - sig)
            }
            ExtActivation::Relu => {
                if x > 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
            ExtActivation::Sigmoid => {
                let s = 1.0 / (1.0 + (-x).exp());
                s * (1.0 - s)
            }
            ExtActivation::Tanh => {
                let t = x.tanh();
                1.0 - t * t
            }
            ExtActivation::Linear => 1.0,
        }
    }
    /// Apply elementwise to a vector in-place.
    pub fn apply_vec(&self, v: &mut [f64]) {
        for x in v.iter_mut() {
            *x = self.apply(*x);
        }
    }
}
