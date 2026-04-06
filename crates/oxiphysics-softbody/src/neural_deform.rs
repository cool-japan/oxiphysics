// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Learned / neural deformation models for soft bodies.
//!
//! Provides:
//!
//! - A simple multi-layer perceptron ([`NeuralDeformNet`]) that maps rest
//!   positions to displacement vectors.
//! - Positional (Fourier feature) encoding ([`PositionalEncoding`]).
//! - PCA-based shape deformation basis ([`ShapeDeformationBasis`]).
//! - Latent-space deformer that combines the basis with an MLP decoder
//!   ([`LatentSpaceDeformer`]).
//! - Blend-shapes ([`blend_shapes`]), delta-mush Laplacian smoothing
//!   ([`delta_mush_smooth`]), and linear blend skinning ([`linear_blend_skinning`]).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Small 3-D helpers (plain arrays, no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-D vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors: a − b.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by a scalar.
#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// 3×3 matrix–vector product.
#[inline]
fn mat3_mul_vec3(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ---------------------------------------------------------------------------
// ActivationFn
// ---------------------------------------------------------------------------

/// Activation function for neural-network layers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ActivationFn {
    /// Rectified linear unit: max(0, x).
    Relu,
    /// Hyperbolic tangent.
    Tanh,
    /// Logistic sigmoid: 1 / (1 + e^(−x)).
    Sigmoid,
    /// Identity (no activation).
    Linear,
}

impl ActivationFn {
    /// Apply the activation function to a scalar value.
    pub fn apply(&self, x: f64) -> f64 {
        match self {
            ActivationFn::Relu => x.max(0.0),
            ActivationFn::Tanh => x.tanh(),
            ActivationFn::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            ActivationFn::Linear => x,
        }
    }

    /// Derivative of the activation function with respect to x.
    ///
    /// For `Relu` the subgradient at 0 is taken as 0.
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

// ---------------------------------------------------------------------------
// NeuralLayer
// ---------------------------------------------------------------------------

/// A single fully-connected (dense) layer of a neural network.
///
/// Stores a weight matrix `weights[out][in]` and a bias vector `bias[out]`.
#[derive(Debug, Clone)]
pub struct NeuralLayer {
    /// Weight matrix, indexed as `weights[output_neuron][input_neuron]`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector, one entry per output neuron.
    pub bias: Vec<f64>,
    /// Activation function applied element-wise to the pre-activation output.
    pub activation: ActivationFn,
}

impl NeuralLayer {
    /// Create a new layer with zero-initialised weights and biases.
    ///
    /// # Arguments
    /// * `in_size`  – number of input features.
    /// * `out_size` – number of output neurons.
    pub fn new(in_size: usize, out_size: usize) -> Self {
        Self {
            weights: vec![vec![0.0; in_size]; out_size],
            bias: vec![0.0; out_size],
            activation: ActivationFn::Relu,
        }
    }

    /// Create a layer with a specified activation function.
    pub fn with_activation(in_size: usize, out_size: usize, activation: ActivationFn) -> Self {
        Self {
            weights: vec![vec![0.0; in_size]; out_size],
            bias: vec![0.0; out_size],
            activation,
        }
    }

    /// Forward pass: computes output = activation(W·x + b).
    ///
    /// Panics in debug mode if `x.len()` does not match the layer's input size.
    pub fn forward(&self, x: &[f64]) -> Vec<f64> {
        let out_size = self.weights.len();
        let mut out = Vec::with_capacity(out_size);
        for i in 0..out_size {
            let mut sum = self.bias[i];
            let row = &self.weights[i];
            for (j, &xj) in x.iter().enumerate() {
                if j < row.len() {
                    sum += row[j] * xj;
                }
            }
            out.push(self.activation.apply(sum));
        }
        out
    }

    /// Number of output neurons.
    pub fn out_size(&self) -> usize {
        self.weights.len()
    }

    /// Number of input features.
    pub fn in_size(&self) -> usize {
        self.weights.first().map(|r| r.len()).unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// NeuralDeformNet
// ---------------------------------------------------------------------------

/// Multi-layer perceptron (MLP) that predicts per-vertex displacements from
/// rest positions.
///
/// Given rest-pose vertex positions (each a 3-D point), the network predicts
/// a displacement Δp for each vertex.  The deformed position is then
/// rest + Δp.
#[derive(Debug, Clone)]
pub struct NeuralDeformNet {
    /// Sequence of fully-connected layers.
    pub layers: Vec<NeuralLayer>,
}

impl NeuralDeformNet {
    /// Construct an MLP with the given layer sizes.
    ///
    /// `layer_sizes` must have at least two entries (input and output).
    /// Hidden layers use `Relu`; the final layer uses `Linear`.
    pub fn new(layer_sizes: &[usize]) -> Self {
        assert!(layer_sizes.len() >= 2, "need at least input + output layer");
        let n = layer_sizes.len();
        let layers = (0..n - 1)
            .map(|i| {
                let activation = if i + 2 == n {
                    ActivationFn::Linear
                } else {
                    ActivationFn::Relu
                };
                NeuralLayer::with_activation(layer_sizes[i], layer_sizes[i + 1], activation)
            })
            .collect();
        Self { layers }
    }

    /// Run the network on a single flat input vector and return the flat output.
    pub fn forward_flat(&self, x: &[f64]) -> Vec<f64> {
        let mut current = x.to_vec();
        for layer in &self.layers {
            current = layer.forward(&current);
        }
        current
    }

    /// Predict deformed positions for all vertices.
    ///
    /// Each vertex's 3-D rest position is flattened and fed through the
    /// network.  The output is interpreted as a displacement vector which is
    /// added to the rest position.
    ///
    /// Assumes the network input size is 3 and output size is 3.
    pub fn forward(&self, rest_pos: &[[f64; 3]]) -> Vec<[f64; 3]> {
        rest_pos
            .iter()
            .map(|p| {
                let disp = self.forward_flat(&[p[0], p[1], p[2]]);
                let dx = disp.first().copied().unwrap_or(0.0);
                let dy = disp.get(1).copied().unwrap_or(0.0);
                let dz = disp.get(2).copied().unwrap_or(0.0);
                [p[0] + dx, p[1] + dy, p[2] + dz]
            })
            .collect()
    }

    /// Overwrite the weights and bias of a specific layer.
    ///
    /// `layer` is a 0-based index.  Panics if the index is out of range.
    pub fn set_weights_layer(&mut self, layer: usize, weights: Vec<Vec<f64>>, bias: Vec<f64>) {
        self.layers[layer].weights = weights;
        self.layers[layer].bias = bias;
    }
}

// ---------------------------------------------------------------------------
// PositionalEncoding
// ---------------------------------------------------------------------------

/// Fourier feature positional encoding for NeRF-style networks.
///
/// Encodes a 3-D point as a vector of sinusoidal features:
/// \[sin(2^0·π·x), cos(2^0·π·x), …, sin(2^(k-1)·π·x), cos(2^(k-1)·π·x)\]
/// for each coordinate.
#[derive(Debug, Clone)]
pub struct PositionalEncoding {
    /// Number of frequency bands k.
    pub n_freqs: usize,
}

impl PositionalEncoding {
    /// Create a `PositionalEncoding` with `n_freqs` frequency bands.
    pub fn new(n_freqs: usize) -> Self {
        Self { n_freqs }
    }

    /// Encode a 3-D point into a sinusoidal feature vector.
    ///
    /// Output length = `output_size()`.
    pub fn encode(&self, x: &[f64; 3]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.output_size());
        for coord in x.iter() {
            for k in 0..self.n_freqs {
                let freq = (2_u64.pow(k as u32) as f64) * PI;
                out.push((freq * coord).sin());
                out.push((freq * coord).cos());
            }
        }
        out
    }

    /// Size of the encoded output vector = 3 × 2 × n_freqs.
    pub fn output_size(&self) -> usize {
        3 * 2 * self.n_freqs
    }
}

// ---------------------------------------------------------------------------
// ShapeDeformationBasis
// ---------------------------------------------------------------------------

/// PCA / learned shape deformation basis.
///
/// Represents a set of shape modes (principal components) and a mean shape.
/// New shapes can be reconstructed as: `shape = mean + Σ coeffs[i] * modes[i]`.
#[derive(Debug, Clone)]
pub struct ShapeDeformationBasis {
    /// Shape basis modes (principal components), each with the same number of
    /// vertices as `mean_shape`.
    pub modes: Vec<Vec<[f64; 3]>>,
    /// Mean (rest) shape.
    pub mean_shape: Vec<[f64; 3]>,
}

impl ShapeDeformationBasis {
    /// Create a `ShapeDeformationBasis` from a mean shape and a set of modes.
    pub fn new(mean: Vec<[f64; 3]>, modes: Vec<Vec<[f64; 3]>>) -> Self {
        Self {
            modes,
            mean_shape: mean,
        }
    }

    /// Reconstruct a shape from latent coefficients.
    ///
    /// `coeffs[i]` scales `modes[i]`.  Extra coefficients are ignored;
    /// missing coefficients default to 0.
    pub fn reconstruct(&self, coeffs: &[f64]) -> Vec<[f64; 3]> {
        let n = self.mean_shape.len();
        let mut out = self.mean_shape.clone();
        for (i, coeff) in coeffs.iter().enumerate() {
            if i >= self.modes.len() {
                break;
            }
            let mode = &self.modes[i];
            for (j, p) in out.iter_mut().enumerate() {
                if j < mode.len() {
                    p[0] += coeff * mode[j][0];
                    p[1] += coeff * mode[j][1];
                    p[2] += coeff * mode[j][2];
                }
            }
        }
        let _ = n;
        out
    }

    /// Project a shape onto the basis to obtain latent coefficients.
    ///
    /// Computes `coeffs[i] = <shape − mean, modes[i]>` (assuming orthonormal
    /// modes, which is the case for PCA).
    pub fn project(&self, shape: &[[f64; 3]]) -> Vec<f64> {
        self.modes
            .iter()
            .map(|mode| {
                let mut dot = 0.0;
                for (j, p) in shape.iter().enumerate() {
                    if j < self.mean_shape.len() && j < mode.len() {
                        let d = sub3(*p, self.mean_shape[j]);
                        dot += d[0] * mode[j][0] + d[1] * mode[j][1] + d[2] * mode[j][2];
                    }
                }
                dot
            })
            .collect()
    }

    /// Compute a shape basis from sample shapes using PCA (via covariance).
    ///
    /// # Arguments
    /// * `shapes`  – slice of example shapes, each a `Vec<[f64;3]>`.
    /// * `n_modes` – number of principal modes to retain.
    ///
    /// Uses a simple power-iteration approach for the covariance eigenvectors
    /// (suitable for small vertex counts; replace with a proper SVD for
    /// large meshes).
    pub fn from_sample_shapes(shapes: &[Vec<[f64; 3]>], n_modes: usize) -> Self {
        assert!(!shapes.is_empty(), "need at least one sample shape");
        let n_verts = shapes[0].len();
        let n_shapes = shapes.len();

        // Compute mean shape.
        let mut mean = vec![[0.0_f64; 3]; n_verts];
        for shape in shapes {
            for (j, p) in shape.iter().enumerate() {
                mean[j][0] += p[0];
                mean[j][1] += p[1];
                mean[j][2] += p[2];
            }
        }
        let inv_n = 1.0 / n_shapes as f64;
        for p in mean.iter_mut() {
            p[0] *= inv_n;
            p[1] *= inv_n;
            p[2] *= inv_n;
        }

        // Build centred data matrix X: shape (n_shapes, 3*n_verts).
        let dim = 3 * n_verts;
        let mut x_data: Vec<Vec<f64>> = Vec::with_capacity(n_shapes);
        for shape in shapes {
            let mut row = Vec::with_capacity(dim);
            for (j, p) in shape.iter().enumerate() {
                row.push(p[0] - mean[j][0]);
                row.push(p[1] - mean[j][1]);
                row.push(p[2] - mean[j][2]);
            }
            x_data.push(row);
        }

        // Compute n_modes principal modes via deflation / power iteration on X^T X.
        let actual_modes = n_modes.min(n_shapes).min(dim);
        let mut modes: Vec<Vec<[f64; 3]>> = Vec::with_capacity(actual_modes);
        let mut residuals = x_data.clone();

        for _ in 0..actual_modes {
            // Initialise vector randomly (deterministic: use first residual).
            let mut v: Vec<f64> = residuals[0].clone();
            // Normalise.
            let len: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if len < 1e-30 {
                break;
            }
            for vi in v.iter_mut() {
                *vi /= len;
            }

            // Power iteration.
            for _ in 0..30 {
                // u = X v  (shape n_shapes)
                let u: Vec<f64> = residuals
                    .iter()
                    .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
                    .collect();
                // v_new = X^T u / ||X^T u||
                let mut v_new = vec![0.0_f64; dim];
                for (row, &ui) in residuals.iter().zip(u.iter()) {
                    for (k, &rk) in row.iter().enumerate() {
                        v_new[k] += ui * rk;
                    }
                }
                let nrm: f64 = v_new.iter().map(|x| x * x).sum::<f64>().sqrt();
                if nrm < 1e-30 {
                    break;
                }
                for vi in v_new.iter_mut() {
                    *vi /= nrm;
                }
                v = v_new;
            }

            // Convert flat vector v back to per-vertex format.
            let mode_verts: Vec<[f64; 3]> = (0..n_verts)
                .map(|j| [v[3 * j], v[3 * j + 1], v[3 * j + 2]])
                .collect();
            modes.push(mode_verts);

            // Deflate residuals: subtract projection onto v.
            for row in residuals.iter_mut() {
                let proj: f64 = row.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
                for (k, rk) in row.iter_mut().enumerate() {
                    *rk -= proj * v[k];
                }
            }
        }

        Self {
            modes,
            mean_shape: mean,
        }
    }
}

// ---------------------------------------------------------------------------
// LatentSpaceDeformer
// ---------------------------------------------------------------------------

/// Combines a `ShapeDeformationBasis` with a neural MLP decoder.
///
/// A latent vector is first decoded by the MLP into shape coefficients which
/// are then used to reconstruct a mesh via the basis.
#[derive(Debug, Clone)]
pub struct LatentSpaceDeformer {
    /// Shape basis used for mesh reconstruction.
    pub basis: ShapeDeformationBasis,
    /// MLP decoder mapping latent codes to shape-basis coefficients.
    pub decoder: NeuralDeformNet,
}

impl LatentSpaceDeformer {
    /// Create a `LatentSpaceDeformer`.
    pub fn new(basis: ShapeDeformationBasis, decoder: NeuralDeformNet) -> Self {
        Self { basis, decoder }
    }

    /// Decode a latent code into a set of deformed vertex positions.
    ///
    /// The latent code is passed through the MLP to produce shape-basis
    /// coefficients, which are then used to reconstruct the mesh.
    pub fn decode(&self, latent: &[f64]) -> Vec<[f64; 3]> {
        let coeffs = self.decoder.forward_flat(latent);
        self.basis.reconstruct(&coeffs)
    }
}

// ---------------------------------------------------------------------------
// Blend shapes
// ---------------------------------------------------------------------------

/// Blend-shape deformation: weighted sum of target shapes.
///
/// Computes `Σ weights[i] * shapes[i]` element-wise.  All shapes must have the
/// same vertex count.  Weights need not sum to 1.
pub fn blend_shapes(shapes: &[Vec<[f64; 3]>], weights: &[f64]) -> Vec<[f64; 3]> {
    assert!(!shapes.is_empty(), "need at least one shape");
    let n = shapes[0].len();
    let mut out = vec![[0.0_f64; 3]; n];
    for (shape, &w) in shapes.iter().zip(weights.iter()) {
        for (j, p) in shape.iter().enumerate() {
            if j < n {
                out[j] = add3(out[j], scale3(*p, w));
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Delta-mush smoothing
// ---------------------------------------------------------------------------

/// Delta-mush Laplacian smoothing.
///
/// Applies `iterations` rounds of umbrella-operator smoothing with weight
/// `lambda` ∈ (0, 1] to the vertex positions, using the supplied adjacency
/// list.  The delta-mush correction (adding back the original delta) is *not*
/// included here; this function implements the core Laplacian pass used inside
/// the algorithm.
pub fn delta_mush_smooth(
    positions: &[[f64; 3]],
    adjacency: &[Vec<usize>],
    lambda: f64,
    iterations: usize,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut current = positions.to_vec();
    let mut next = current.clone();

    for _ in 0..iterations {
        for i in 0..n {
            let neighbours = &adjacency[i];
            if neighbours.is_empty() {
                next[i] = current[i];
                continue;
            }
            let inv_k = 1.0 / neighbours.len() as f64;
            let mut avg = [0.0_f64; 3];
            for &j in neighbours {
                avg = add3(avg, current[j]);
            }
            avg = scale3(avg, inv_k);
            // Blend: p_new = (1 − λ) * p + λ * avg
            let blended = add3(scale3(current[i], 1.0 - lambda), scale3(avg, lambda));
            next[i] = blended;
        }
        current.clone_from(&next);
    }
    current
}

// ---------------------------------------------------------------------------
// Linear blend skinning
// ---------------------------------------------------------------------------

/// Linear Blend Skinning (LBS) deformation.
///
/// Deforms `rest` vertex positions using a set of bone transforms.
///
/// # Arguments
/// * `rest`    – rest-pose vertex positions.
/// * `bones`   – slice of `(translation, rotation_matrix)` bone transforms.
///   The rotation matrix is row-major `[[f64;3\];3]`.
/// * `weights` – per-vertex skin weights as `(bone_index, weight)` pairs.
///
/// The deformed position of vertex *i* is:
/// Σ_b  w_ib · (R_b · p_i + t_b)
pub fn linear_blend_skinning(
    rest: &[[f64; 3]],
    bones: &[([f64; 3], [[f64; 3]; 3])],
    weights: &[Vec<(usize, f64)>],
) -> Vec<[f64; 3]> {
    rest.iter()
        .zip(weights.iter())
        .map(|(p, vw)| {
            let mut deformed = [0.0_f64; 3];
            for &(bone_idx, w) in vw {
                if bone_idx >= bones.len() {
                    continue;
                }
                let (translation, rotation) = &bones[bone_idx];
                let rotated = mat3_mul_vec3(rotation, *p);
                let transformed = add3(rotated, *translation);
                deformed = add3(deformed, scale3(transformed, w));
            }
            deformed
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── ActivationFn ─────────────────────────────────────────────────────────

    #[test]
    fn relu_positive_input() {
        assert!((ActivationFn::Relu.apply(2.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn relu_negative_input_is_zero() {
        assert_eq!(ActivationFn::Relu.apply(-1.0), 0.0);
    }

    #[test]
    fn relu_derivative_positive() {
        assert!((ActivationFn::Relu.derivative(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn relu_derivative_negative_is_zero() {
        assert_eq!(ActivationFn::Relu.derivative(-1.0), 0.0);
    }

    #[test]
    fn tanh_at_zero_is_zero() {
        assert!(ActivationFn::Tanh.apply(0.0).abs() < 1e-12);
    }

    #[test]
    fn tanh_derivative_at_zero_is_one() {
        assert!((ActivationFn::Tanh.derivative(0.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sigmoid_at_zero_is_half() {
        assert!((ActivationFn::Sigmoid.apply(0.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn sigmoid_derivative_at_zero() {
        // σ'(0) = 0.25
        assert!((ActivationFn::Sigmoid.derivative(0.0) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn linear_apply_identity() {
        assert!((ActivationFn::Linear.apply(3.7) - 3.7).abs() < 1e-12);
    }

    #[test]
    fn linear_derivative_is_one() {
        assert!((ActivationFn::Linear.derivative(99.0) - 1.0).abs() < 1e-12);
    }

    // ── NeuralLayer ───────────────────────────────────────────────────────────

    #[test]
    fn neural_layer_zero_init_forward() {
        let layer = NeuralLayer::new(3, 2);
        // Zero weights + zero bias → ReLU(0) = 0
        let out = layer.forward(&[1.0, 2.0, 3.0]);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], 0.0);
        assert_eq!(out[1], 0.0);
    }

    #[test]
    fn neural_layer_size_accessors() {
        let layer = NeuralLayer::new(4, 6);
        assert_eq!(layer.in_size(), 4);
        assert_eq!(layer.out_size(), 6);
    }

    #[test]
    fn neural_layer_custom_weights() {
        let mut layer = NeuralLayer::with_activation(2, 1, ActivationFn::Linear);
        layer.weights[0] = vec![1.0, 2.0];
        layer.bias[0] = 0.5;
        let out = layer.forward(&[1.0, 1.0]);
        // 1*1 + 2*1 + 0.5 = 3.5
        assert!((out[0] - 3.5).abs() < 1e-12);
    }

    #[test]
    fn neural_layer_relu_clamps_negative() {
        let mut layer = NeuralLayer::new(1, 1);
        layer.bias[0] = -5.0; // output pre-activation = -5
        let out = layer.forward(&[0.0]);
        assert_eq!(out[0], 0.0, "ReLU should clamp to 0");
    }

    // ── NeuralDeformNet ───────────────────────────────────────────────────────

    #[test]
    fn neural_deform_net_zero_displacement() {
        // All-zero weights → zero displacement → deformed = rest.
        let net = NeuralDeformNet::new(&[3, 8, 3]);
        let rest = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let deformed = net.forward(&rest);
        for (d, r) in deformed.iter().zip(rest.iter()) {
            for k in 0..3 {
                assert!(
                    (d[k] - r[k]).abs() < 1e-12,
                    "zero net should not displace vertices"
                );
            }
        }
    }

    #[test]
    fn neural_deform_net_set_weights() {
        let mut net = NeuralDeformNet::new(&[3, 3]);
        net.set_weights_layer(0, vec![vec![1.0, 0.0, 0.0]; 3], vec![0.0; 3]);
        let out = net.forward_flat(&[5.0, 0.0, 0.0]);
        // Linear layer: W*x = [5, 5, 5] (weight row [1,0,0] * [5,0,0])
        assert!((out[0] - 5.0).abs() < 1e-12, "out[0]={}", out[0]);
    }

    #[test]
    fn neural_deform_net_forward_output_count() {
        let net = NeuralDeformNet::new(&[3, 16, 16, 3]);
        let rest: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let deformed = net.forward(&rest);
        assert_eq!(deformed.len(), rest.len());
    }

    // ── PositionalEncoding ────────────────────────────────────────────────────

    #[test]
    fn positional_encoding_output_size() {
        let pe = PositionalEncoding::new(4);
        assert_eq!(pe.output_size(), 3 * 2 * 4);
        let encoded = pe.encode(&[0.1, 0.2, 0.3]);
        assert_eq!(encoded.len(), pe.output_size());
    }

    #[test]
    fn positional_encoding_zero_point() {
        // sin(k·π·0) = 0, cos(k·π·0) = 1 for all k
        let pe = PositionalEncoding::new(3);
        let enc = pe.encode(&[0.0, 0.0, 0.0]);
        for i in 0..3 {
            for k in 0..3 {
                let sin_val = enc[i * 6 + 2 * k];
                let cos_val = enc[i * 6 + 2 * k + 1];
                assert!(sin_val.abs() < 1e-12, "sin at 0 should be 0, got {sin_val}");
                assert!(
                    (cos_val - 1.0).abs() < 1e-12,
                    "cos at 0 should be 1, got {cos_val}"
                );
            }
        }
    }

    #[test]
    fn positional_encoding_different_freqs() {
        let pe1 = PositionalEncoding::new(1);
        let pe2 = PositionalEncoding::new(4);
        let pt = [0.5, 0.5, 0.5];
        assert_ne!(pe1.encode(&pt).len(), pe2.encode(&pt).len());
    }

    // ── ShapeDeformationBasis ─────────────────────────────────────────────────

    #[test]
    fn basis_reconstruct_zero_coeffs_gives_mean() {
        let mean = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mode0 = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let basis = ShapeDeformationBasis::new(mean.clone(), vec![mode0]);
        let rec = basis.reconstruct(&[0.0]);
        for (r, m) in rec.iter().zip(mean.iter()) {
            for k in 0..3 {
                assert!((r[k] - m[k]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn basis_reconstruct_unit_coeff_adds_mode() {
        let mean = vec![[0.0, 0.0, 0.0]];
        let mode0 = vec![[1.0, 2.0, 3.0]];
        let basis = ShapeDeformationBasis::new(mean, vec![mode0]);
        let rec = basis.reconstruct(&[2.0]);
        assert!((rec[0][0] - 2.0).abs() < 1e-12);
        assert!((rec[0][1] - 4.0).abs() < 1e-12);
        assert!((rec[0][2] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn basis_project_recovers_coeff() {
        let mean = vec![[0.0, 0.0, 0.0]];
        // Orthonormal mode: unit vector in x
        let mode0 = vec![[1.0, 0.0, 0.0]];
        let basis = ShapeDeformationBasis::new(mean.clone(), vec![mode0]);
        // shape displaced by 3 in x
        let shape = vec![[3.0, 0.0, 0.0]];
        let coeffs = basis.project(&shape);
        assert!((coeffs[0] - 3.0).abs() < 1e-12, "coeff={}", coeffs[0]);
    }

    #[test]
    fn basis_from_samples_returns_correct_n_modes() {
        let shapes: Vec<Vec<[f64; 3]>> = (0..5)
            .map(|i| vec![[i as f64, 0.0, 0.0], [0.0, i as f64, 0.0]])
            .collect();
        let basis = ShapeDeformationBasis::from_sample_shapes(&shapes, 2);
        assert!(basis.modes.len() <= 2);
        assert_eq!(basis.mean_shape.len(), 2);
    }

    #[test]
    fn basis_from_samples_mean_is_average() {
        // Three identical shapes: mean must equal that shape.
        let shape = vec![[2.0_f64, 0.0, 0.0]];
        let shapes = vec![shape.clone(), shape.clone(), shape.clone()];
        let basis = ShapeDeformationBasis::from_sample_shapes(&shapes, 1);
        assert!((basis.mean_shape[0][0] - 2.0).abs() < 1e-10);
    }

    // ── LatentSpaceDeformer ───────────────────────────────────────────────────

    #[test]
    fn latent_space_deformer_decode_returns_correct_count() {
        let mean = vec![[0.0_f64; 3]; 4];
        let modes = vec![vec![[1.0_f64, 0.0, 0.0]; 4]];
        let basis = ShapeDeformationBasis::new(mean, modes);
        // MLP: 2 inputs → 1 output (number of modes)
        let decoder = NeuralDeformNet::new(&[2, 1]);
        let ld = LatentSpaceDeformer::new(basis, decoder);
        let decoded = ld.decode(&[0.5, 0.5]);
        assert_eq!(decoded.len(), 4);
    }

    // ── blend_shapes ─────────────────────────────────────────────────────────

    #[test]
    fn blend_shapes_equal_weights_is_average() {
        let s0 = vec![[0.0_f64, 0.0, 0.0]];
        let s1 = vec![[2.0_f64, 0.0, 0.0]];
        let blended = blend_shapes(&[s0, s1], &[0.5, 0.5]);
        assert!(
            (blended[0][0] - 1.0).abs() < 1e-12,
            "blended={}",
            blended[0][0]
        );
    }

    #[test]
    fn blend_shapes_single_weight_one() {
        let s0 = vec![[3.0_f64, 4.0, 5.0]];
        let blended = blend_shapes(std::slice::from_ref(&s0), &[1.0]);
        assert_eq!(blended[0], s0[0]);
    }

    #[test]
    fn blend_shapes_zero_weights_gives_zero() {
        let s0 = vec![[1.0_f64, 2.0, 3.0]];
        let blended = blend_shapes(&[s0], &[0.0]);
        assert_eq!(blended[0], [0.0, 0.0, 0.0]);
    }

    // ── delta_mush_smooth ─────────────────────────────────────────────────────

    #[test]
    fn delta_mush_smooth_isolated_vertex_unchanged() {
        // Vertex 0 has no neighbours → position unchanged.
        let positions = vec![[1.0, 2.0, 3.0], [5.0, 0.0, 0.0]];
        let adjacency = vec![vec![], vec![0_usize]];
        let smoothed = delta_mush_smooth(&positions, &adjacency, 0.5, 1);
        assert_eq!(smoothed[0], positions[0]);
    }

    #[test]
    fn delta_mush_smooth_chain_converges() {
        // Three vertices in a line: 0–1–2.
        // After many iterations the interior vertex (1) should approach the
        // average of its neighbours.
        let positions = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let adjacency = vec![vec![1usize], vec![0usize, 2usize], vec![1usize]];
        let smoothed = delta_mush_smooth(&positions, &adjacency, 0.5, 20);
        // Vertex 1 should move toward average of (0, 2) = (0, 0, 0)
        assert!(
            smoothed[1][0] < 9.0,
            "vertex 1 should have moved, x={}",
            smoothed[1][0]
        );
    }

    #[test]
    fn delta_mush_smooth_zero_lambda_no_change() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let adjacency = vec![vec![1_usize], vec![0_usize]];
        let smoothed = delta_mush_smooth(&positions, &adjacency, 0.0, 5);
        for (s, p) in smoothed.iter().zip(positions.iter()) {
            for k in 0..3 {
                assert!((s[k] - p[k]).abs() < 1e-12);
            }
        }
    }

    // ── linear_blend_skinning ─────────────────────────────────────────────────

    #[test]
    fn lbs_identity_transform_no_change() {
        let rest = vec![[1.0, 2.0, 3.0]];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let bones = vec![([0.0_f64, 0.0, 0.0], identity)];
        let weights = vec![vec![(0usize, 1.0)]];
        let deformed = linear_blend_skinning(&rest, &bones, &weights);
        for k in 0..3 {
            assert!((deformed[0][k] - rest[0][k]).abs() < 1e-12);
        }
    }

    #[test]
    fn lbs_translation_only() {
        let rest = vec![[0.0, 0.0, 0.0]];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let bones = vec![([5.0_f64, 3.0, 1.0], identity)];
        let weights = vec![vec![(0usize, 1.0)]];
        let deformed = linear_blend_skinning(&rest, &bones, &weights);
        assert!((deformed[0][0] - 5.0).abs() < 1e-12);
        assert!((deformed[0][1] - 3.0).abs() < 1e-12);
        assert!((deformed[0][2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn lbs_two_bones_weighted() {
        let rest = vec![[1.0, 0.0, 0.0]];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Bone 0: translate +2 in x; Bone 1: translate -2 in x.
        let bones = vec![
            ([2.0_f64, 0.0, 0.0], identity),
            ([-2.0_f64, 0.0, 0.0], identity),
        ];
        // Equal weights → net translation = 0.
        let weights = vec![vec![(0usize, 0.5), (1usize, 0.5)]];
        let deformed = linear_blend_skinning(&rest, &bones, &weights);
        // Deformed: 0.5*(R*p + [2,0,0]) + 0.5*(R*p + [-2,0,0]) = R*p = [1,0,0]
        assert!((deformed[0][0] - 1.0).abs() < 1e-12, "x={}", deformed[0][0]);
    }

    #[test]
    fn lbs_rotation_90_deg_around_z() {
        // Rotation matrix R_z(90°): x→y, y→-x
        let rot90z = [[0.0_f64, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let rest = vec![[1.0, 0.0, 0.0]];
        let bones = vec![([0.0_f64, 0.0, 0.0], rot90z)];
        let weights = vec![vec![(0usize, 1.0)]];
        let deformed = linear_blend_skinning(&rest, &bones, &weights);
        // Expected: [0, 1, 0]
        assert!((deformed[0][0] - 0.0).abs() < 1e-12, "x={}", deformed[0][0]);
        assert!((deformed[0][1] - 1.0).abs() < 1e-12, "y={}", deformed[0][1]);
    }

    #[test]
    fn lbs_out_of_range_bone_index_ignored() {
        let rest = vec![[1.0, 1.0, 1.0]];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let bones = vec![([0.0_f64, 0.0, 0.0], identity)];
        // Weight pointing to bone index 99 (out of range) → should be ignored.
        let weights = vec![vec![(99usize, 1.0)]];
        let deformed = linear_blend_skinning(&rest, &bones, &weights);
        // No valid bone → deformed should be zero.
        assert_eq!(deformed[0], [0.0, 0.0, 0.0]);
    }

    // ── edge cases ───────────────────────────────────────────────────────────

    #[test]
    fn blend_shapes_multiple_verts() {
        let s0 = vec![[0.0_f64, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let s1 = vec![[2.0_f64, 2.0, 2.0], [3.0, 3.0, 3.0]];
        let blended = blend_shapes(&[s0, s1], &[0.5, 0.5]);
        assert!((blended[0][0] - 1.0).abs() < 1e-12);
        assert!((blended[1][0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn positional_encoding_known_value() {
        let pe = PositionalEncoding::new(1);
        // For x = 0.5: sin(π * 0.5) = sin(π/2) = 1, cos = 0
        let enc = pe.encode(&[0.5, 0.0, 0.0]);
        assert!((enc[0] - 1.0).abs() < 1e-10, "sin(π*0.5)≈1, got {}", enc[0]);
        assert!(enc[1].abs() < 1e-10, "cos(π*0.5)≈0, got {}", enc[1]);
    }

    #[test]
    fn neural_layer_tanh_output_bounded() {
        let mut layer = NeuralLayer::with_activation(1, 3, ActivationFn::Tanh);
        layer.weights[0][0] = 100.0;
        layer.weights[1][0] = -100.0;
        layer.weights[2][0] = 0.0;
        let out = layer.forward(&[1.0]);
        assert!((out[0] - 1.0).abs() < 1e-6, "tanh(100)≈1");
        assert!((out[1] + 1.0).abs() < 1e-6, "tanh(-100)≈-1");
    }

    #[test]
    fn basis_from_samples_reconstruct_close_to_sample() {
        // With 2 identical samples the mean == sample and 0 modes (no variance).
        let shapes = vec![vec![[1.0_f64, 0.0, 0.0]], vec![[1.0_f64, 0.0, 0.0]]];
        let basis = ShapeDeformationBasis::from_sample_shapes(&shapes, 2);
        // Reconstruct with zero coefficients → should give mean = [1,0,0].
        let rec = basis.reconstruct(&[]);
        assert!((rec[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn delta_mush_iterations_zero_is_identity() {
        let positions = vec![[3.0, 1.0, 2.0], [0.0, 0.0, 0.0]];
        let adjacency = vec![vec![1_usize], vec![0_usize]];
        let smoothed = delta_mush_smooth(&positions, &adjacency, 0.8, 0);
        assert_eq!(smoothed, positions);
    }
}
