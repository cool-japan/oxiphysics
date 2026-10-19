// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural network-based softbody simulation.
//!
//! This module implements a physics-informed neural network (PINN) approach to
//! softbody deformation, including:
//!
//! - **Learned constitutive model**: neural network predicts stress from strain
//! - **PINN deformation**: physics-informed loss penalises equilibrium residual
//! - **Learned collision response**: data-driven impulse prediction
//! - **Material parameter identification**: gradient descent on material params
//! - **Neural elasticity graph**: graph-neural-network-style node/edge updates
//! - **Training data generation**: forward physics simulation to build datasets
//! - **Inference-time correction**: residual correction at simulation time
//! - **Neural pose-driven deformation**: blend-shape + neural residual
//! - **Latent space deformation**: encoder/decoder over vertex positions
//! - **Real-time neural softbody**: low-latency inference loop
//!
//! Strains and stresses are Voigt 6-vectors `[xx, yy, zz, xy, xz, yz]` with
//! tensor (not engineering) shear components.

use crate::{Error, Result};

mod pde_loss;
mod rest_shape;
#[cfg(test)]
mod tests;

pub use rest_shape::TetRestShape;

// ---------------------------------------------------------------------------
// Internal math helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// In-place vector addition: `a += b`.
#[inline]
fn vec_add_inplace(a: &mut [f64], b: &[f64]) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x += y;
    }
}

/// ReLU activation.
#[inline]
fn relu(x: f64) -> f64 {
    x.max(0.0)
}

/// Derivative of ReLU.
#[inline]
fn relu_grad(x: f64) -> f64 {
    if x > 0.0 { 1.0 } else { 0.0 }
}

/// Tanh activation.
#[inline]
fn tanh_act(x: f64) -> f64 {
    x.tanh()
}

/// Dense forward pass: `y = W x + b` (row-major weights, size `n_out × n_in`).
fn dense_forward(w: &[f64], b: &[f64], x: &[f64], n_in: usize, n_out: usize) -> Vec<f64> {
    assert_eq!(w.len(), n_out * n_in);
    assert_eq!(b.len(), n_out);
    assert_eq!(x.len(), n_in);
    (0..n_out)
        .map(|i| b[i] + (0..n_in).map(|j| w[i * n_in + j] * x[j]).sum::<f64>())
        .collect()
}

/// Dense backward pass (gradient w.r.t. input).
/// Returns `grad_x` of shape `[n_in]`.
fn dense_backward_input(w: &[f64], grad_y: &[f64], n_in: usize, n_out: usize) -> Vec<f64> {
    let mut grad_x = vec![0.0f64; n_in];
    for i in 0..n_out {
        for j in 0..n_in {
            grad_x[j] += w[i * n_in + j] * grad_y[i];
        }
    }
    grad_x
}

/// Compute gradient of weights from `grad_y` and `x`.
fn dense_weight_grad(x: &[f64], grad_y: &[f64], n_in: usize, n_out: usize) -> Vec<f64> {
    let mut gw = vec![0.0f64; n_out * n_in];
    for i in 0..n_out {
        for j in 0..n_in {
            gw[i * n_in + j] = grad_y[i] * x[j];
        }
    }
    gw
}

// ---------------------------------------------------------------------------
// NeuralConstitutiveModel
// ---------------------------------------------------------------------------

/// Seed of the deterministic He initialisation used by
/// [`NeuralConstitutiveModel::new`] ("OXPHYSIC" in ASCII).
pub const DEFAULT_INIT_SEED: u64 = 0x4F58_5048_5953_4943;

/// SplitMix64 pseudo-random generator (Steele, Lea & Flood, 2014).
///
/// Only used for the weight initialisation, so that a model built from a given
/// seed is bit-identical on every platform and independent of the `rand` version.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform sample in `[-limit, limit)` with 53 random bits.
    fn next_symmetric(&mut self, limit: f64) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        limit * (2.0 * unit - 1.0)
    }
}

/// Half-width `sqrt(6 / fan_in)` of the He-uniform weight distribution.
fn he_uniform_limit(fan_in: usize) -> f64 {
    if fan_in == 0 {
        0.0
    } else {
        (6.0 / fan_in as f64).sqrt()
    }
}

/// A 2-layer ReLU network that maps strain (Voigt 6-vector) to stress (6-vector).
///
/// Architecture: `6 → n_hidden (ReLU) → 6`.
#[derive(Debug, Clone, PartialEq)]
pub struct NeuralConstitutiveModel {
    /// Hidden layer size.
    pub n_hidden: usize,
    /// Layer 1 weights (row-major `n_hidden × 6`).
    pub w1: Vec<f64>,
    /// Layer 1 biases (`n_hidden`).
    pub b1: Vec<f64>,
    /// Layer 2 weights (row-major `6 × n_hidden`).
    pub w2: Vec<f64>,
    /// Layer 2 biases (`6`).
    pub b2: Vec<f64>,
}

impl NeuralConstitutiveModel {
    /// Construct a model with all weights and biases zero.
    ///
    /// Such a model predicts zero stress for every strain, and gradient descent
    /// cannot leave that state through the hidden layer: with `w1 = w2 = 0`
    /// every weight gradient except the one of the output bias vanishes.  Use
    /// [`new`](Self::new) for a trainable network.
    pub fn zeros(n_hidden: usize) -> Self {
        Self {
            n_hidden,
            w1: vec![0.0; n_hidden * 6],
            b1: vec![0.0; n_hidden],
            w2: vec![0.0; 6 * n_hidden],
            b2: vec![0.0; 6],
        }
    }

    /// Construct a trainable model with deterministic He-uniform weights.
    ///
    /// Equivalent to [`with_seed`](Self::with_seed) with [`DEFAULT_INIT_SEED`].
    pub fn new(n_hidden: usize) -> Self {
        Self::with_seed(n_hidden, DEFAULT_INIT_SEED)
    }

    /// Construct a model with He-uniform weights drawn from a SplitMix64
    /// generator seeded with `seed`: every weight is uniform in
    /// `±sqrt(6 / fan_in)` (fan-in 6 for the first layer, `n_hidden` for the
    /// second) and all biases are zero.
    ///
    /// The random weights break the symmetry that keeps a zero-initialised
    /// network from learning, and the zero biases make the stress vanish at
    /// zero strain, so the undeformed body starts force free.
    pub fn with_seed(n_hidden: usize, seed: u64) -> Self {
        let mut rng = SplitMix64::new(seed);
        let limit1 = he_uniform_limit(6);
        let limit2 = he_uniform_limit(n_hidden);
        let w1 = (0..n_hidden * 6)
            .map(|_| rng.next_symmetric(limit1))
            .collect();
        let w2 = (0..6 * n_hidden)
            .map(|_| rng.next_symmetric(limit2))
            .collect();
        Self {
            n_hidden,
            w1,
            b1: vec![0.0; n_hidden],
            w2,
            b2: vec![0.0; 6],
        }
    }

    /// Construct a "linear-elastic" model with `stress = scale * strain`.
    ///
    /// A ReLU unit cannot pass a negative input, so every Voigt component uses
    /// a pair of hidden units, `relu(x) - relu(-x) = x`: the network is then
    /// exactly linear for tension and compression alike.  Six components need
    /// `n_hidden >= 12`; with fewer hidden units only the first `n_hidden / 2`
    /// components get a pair and the remaining stress components are zero.
    /// Hidden units beyond twelve have zero weights.
    pub fn linear_elastic(n_hidden: usize, scale: f64) -> Self {
        let mut model = Self::zeros(n_hidden);
        for component in 0..(n_hidden / 2).min(6) {
            let positive = 2 * component;
            let negative = positive + 1;
            model.w1[positive * 6 + component] = 1.0;
            model.w1[negative * 6 + component] = -1.0;
            model.w2[component * n_hidden + positive] = scale;
            model.w2[component * n_hidden + negative] = -scale;
        }
        model
    }

    /// Forward inference: given a Voigt strain vector (length 6), return stress (length 6).
    pub fn forward(&self, strain: &[f64]) -> Vec<f64> {
        assert_eq!(strain.len(), 6);
        let h = dense_forward(&self.w1, &self.b1, strain, 6, self.n_hidden);
        let h_act: Vec<f64> = h.iter().map(|&x| relu(x)).collect();
        dense_forward(&self.w2, &self.b2, &h_act, self.n_hidden, 6)
    }

    /// Forward inference on a fixed-size Voigt strain; returns the Voigt stress.
    pub fn forward_voigt(&self, strain: &[f64; 6]) -> [f64; 6] {
        let mut stress = [0.0f64; 6];
        for (out, value) in stress.iter_mut().zip(self.forward(strain)) {
            *out = value;
        }
        stress
    }

    /// Backward pass: returns `(grad_w1, grad_b1, grad_w2, grad_b2)` for the
    /// squared-error loss `||stress_pred - stress_target||^2` of one sample.
    pub fn backward(
        &self,
        strain: &[f64],
        stress_target: &[f64],
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let stress_pred = self.forward(strain);
        // dL/d_stress_pred = 2*(pred - target)
        let grad_out: Vec<f64> = stress_pred
            .iter()
            .zip(stress_target.iter())
            .map(|(p, t)| 2.0 * (p - t))
            .collect();
        self.backward_from_output_grad(strain, &grad_out)
    }

    /// Backward pass for an arbitrary upstream gradient `grad_out = dL/d(stress)`
    /// (length 6) of the stress predicted for `strain`; returns
    /// `(grad_w1, grad_b1, grad_w2, grad_b2)`.
    ///
    /// This is how losses that are not a plain per-sample error, such as the
    /// equilibrium residual of a stress field, reach the network weights.
    pub fn backward_from_output_grad(
        &self,
        strain: &[f64],
        grad_out: &[f64],
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let h_pre = dense_forward(&self.w1, &self.b1, strain, 6, self.n_hidden);
        let h_act: Vec<f64> = h_pre.iter().map(|&x| relu(x)).collect();

        let gw2 = dense_weight_grad(&h_act, grad_out, self.n_hidden, 6);
        let gb2 = grad_out.to_vec();

        // Backprop through ReLU
        let grad_h_pre_act = dense_backward_input(&self.w2, grad_out, self.n_hidden, 6);
        let grad_h_pre: Vec<f64> = grad_h_pre_act
            .iter()
            .zip(h_pre.iter())
            .map(|(g, &pre)| g * relu_grad(pre))
            .collect();

        let gw1 = dense_weight_grad(strain, &grad_h_pre, 6, self.n_hidden);
        let gb1 = grad_h_pre;

        (gw1, gb1, gw2, gb2)
    }

    /// SGD update step with learning rate `lr`.
    pub fn sgd_step(&mut self, gw1: &[f64], gb1: &[f64], gw2: &[f64], gb2: &[f64], lr: f64) {
        for (w, g) in self.w1.iter_mut().zip(gw1) {
            *w -= lr * g;
        }
        for (b, g) in self.b1.iter_mut().zip(gb1) {
            *b -= lr * g;
        }
        for (w, g) in self.w2.iter_mut().zip(gw2) {
            *w -= lr * g;
        }
        for (b, g) in self.b2.iter_mut().zip(gb2) {
            *b -= lr * g;
        }
    }
}

// ---------------------------------------------------------------------------
// PhysicsInformedLoss
// ---------------------------------------------------------------------------

/// Physics-informed loss combining data MSE with an equilibrium residual penalty.
///
/// `L = L_data + lambda_pde * L_pde`
///
/// `L_data` is the mean squared error between predicted and target stress
/// ([`data_loss`](Self::data_loss)).  `L_pde` is the mean squared residual
/// `||div(sigma) + f||^2` of the static equilibrium equation of the predicted
/// stress field, evaluated with central differences on a regular grid of
/// spacing [`fd_step`](Self::fd_step) (see
/// [`pde_loss`](Self::pde_loss) and [`compute_field`](Self::compute_field)).
/// The divergence needs spatial neighbours, so the PDE term only exists for a
/// stress field sampled on a grid.
#[derive(Debug, Clone)]
pub struct PhysicsInformedLoss {
    /// Weight on the PDE residual term.
    pub lambda_pde: f64,
    /// Grid spacing `h` of the stress field, the finite-difference step of
    /// the stress divergence.
    pub fd_step: f64,
    /// Regular-grid layout `[nx, ny, nz]` (flat index `i + nx * (j + ny * k)`)
    /// of the batches given to [`NeuralPinnSolver::train_step`].  With
    /// `Some(dims)` the batch is the stress field of the grid nodes and the PDE
    /// term takes part in training; with `None` the samples are independent and
    /// only the data term is trained.
    pub grid_dims: Option<[usize; 3]>,
}

impl PhysicsInformedLoss {
    /// Create a PINN loss with given `lambda_pde` and grid spacing `fd_step`,
    /// without a grid layout (data term only until [`with_grid`](Self::with_grid)).
    pub fn new(lambda_pde: f64, fd_step: f64) -> Self {
        Self {
            lambda_pde,
            fd_step,
            grid_dims: None,
        }
    }

    /// Declare that training batches are stress fields on a regular grid of
    /// `dims = [nx, ny, nz]` nodes, which activates the PDE term.
    pub fn with_grid(mut self, dims: [usize; 3]) -> Self {
        self.grid_dims = Some(dims);
        self
    }

    /// Mean squared error between predicted and target stress components.
    ///
    /// Returns 0 for empty slices (never NaN).  The comparison covers the
    /// common length of the two slices.
    pub fn data_loss(&self, stress_pred: &[f64], stress_target: &[f64]) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for (p, t) in stress_pred.iter().zip(stress_target.iter()) {
            sum += (p - t) * (p - t);
            count += 1;
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
}

// ---------------------------------------------------------------------------
// LearnedCollisionResponse
// ---------------------------------------------------------------------------

/// Learned collision response: given contact normal, relative velocity, and
/// penetration depth, predicts the corrective impulse via a small network.
#[derive(Debug, Clone)]
pub struct LearnedCollisionResponse {
    /// Weights for the single hidden layer (row-major `n_hidden × 5`).
    pub w1: Vec<f64>,
    /// Biases for hidden layer.
    pub b1: Vec<f64>,
    /// Output weights (row-major `3 × n_hidden`).
    pub w2: Vec<f64>,
    /// Output biases.
    pub b2: Vec<f64>,
    /// Hidden layer size.
    pub n_hidden: usize,
}

impl LearnedCollisionResponse {
    /// Construct with given hidden size; weights default to zero.
    pub fn new(n_hidden: usize) -> Self {
        Self {
            w1: vec![0.0; n_hidden * 5],
            b1: vec![0.0; n_hidden],
            w2: vec![0.0; 3 * n_hidden],
            b2: vec![0.0; 3],
            n_hidden,
        }
    }

    /// Predict the impulse vector `[ix, iy, iz]` from:
    ///
    /// - `contact_normal`: `[nx, ny, nz]`
    /// - `rel_vel`: relative normal velocity (scalar)
    /// - `penetration`: penetration depth (scalar, positive means overlap)
    pub fn predict(&self, contact_normal: &[f64; 3], rel_vel: f64, penetration: f64) -> [f64; 3] {
        let input = [
            contact_normal[0],
            contact_normal[1],
            contact_normal[2],
            rel_vel,
            penetration,
        ];
        let h = dense_forward(&self.w1, &self.b1, &input, 5, self.n_hidden);
        let h_act: Vec<f64> = h.iter().map(|&x| relu(x)).collect();
        let out = dense_forward(&self.w2, &self.b2, &h_act, self.n_hidden, 3);
        [out[0], out[1], out[2]]
    }
}

// ---------------------------------------------------------------------------
// MaterialParameterIdentifier
// ---------------------------------------------------------------------------

/// Data-driven material parameter identification via gradient descent.
///
/// Optimises Young's modulus `E` and Poisson ratio `nu` to minimise the
/// mean-squared error between simulated and observed vertex displacements.
#[derive(Debug, Clone)]
pub struct MaterialParameterIdentifier {
    /// Current estimate of Young's modulus (Pa).
    pub young_modulus: f64,
    /// Current estimate of Poisson ratio.
    pub poisson_ratio: f64,
    /// Learning rate for gradient descent.
    pub lr: f64,
    /// Regularisation weight (L2 on parameters).
    pub reg_weight: f64,
}

impl MaterialParameterIdentifier {
    /// Create a new identifier with initial guesses and learning rate.
    pub fn new(young_modulus: f64, poisson_ratio: f64, lr: f64, reg_weight: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            lr,
            reg_weight,
        }
    }

    /// Compute Lamé first parameter `lambda` from current E and nu.
    pub fn lame_lambda(&self) -> f64 {
        self.young_modulus * self.poisson_ratio
            / ((1.0 + self.poisson_ratio) * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Compute shear modulus `mu` from current E and nu.
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Perform one gradient-descent step given simulated and observed displacements.
    ///
    /// Uses finite differences to estimate the gradient of the MSE loss with
    /// respect to `E` and `nu`.
    pub fn step(
        &mut self,
        sim_disp: &[[f64; 3]],
        obs_disp: &[[f64; 3]],
        sim_fn: &dyn Fn(f64, f64) -> Vec<[f64; 3]>,
    ) {
        let loss = |e: f64, nu: f64| -> f64 {
            let disp = sim_fn(e, nu);
            disp.iter()
                .zip(obs_disp.iter())
                .map(|(s, o)| (s[0] - o[0]).powi(2) + (s[1] - o[1]).powi(2) + (s[2] - o[2]).powi(2))
                .sum::<f64>()
                / disp.len() as f64
        };

        let _ = sim_disp; // used via sim_fn
        let eps = 1e-4;
        let l0 = loss(self.young_modulus, self.poisson_ratio);
        let de = (loss(self.young_modulus + eps, self.poisson_ratio) - l0) / eps;
        let dnu = (loss(self.young_modulus, self.poisson_ratio + eps) - l0) / eps;

        self.young_modulus -= self.lr * (de + self.reg_weight * self.young_modulus);
        self.poisson_ratio -= self.lr * (dnu + self.reg_weight * self.poisson_ratio);

        // Clamp to physically valid range
        self.young_modulus = self.young_modulus.max(1.0);
        self.poisson_ratio = self.poisson_ratio.clamp(0.0, 0.499);
    }
}

// ---------------------------------------------------------------------------
// NeuralElasticityGraph
// ---------------------------------------------------------------------------

/// A node in the neural elasticity graph.
#[derive(Debug, Clone)]
pub struct NeuralElasticNode {
    /// 3-D position of the node.
    pub position: [f64; 3],
    /// Node feature vector (latent state).
    pub features: Vec<f64>,
    /// Mass of this node.
    pub mass: f64,
}

impl NeuralElasticNode {
    /// Create a new node at `position` with given feature dimensionality and mass.
    pub fn new(position: [f64; 3], feature_dim: usize, mass: f64) -> Self {
        Self {
            position,
            features: vec![0.0; feature_dim],
            mass,
        }
    }
}

/// An edge in the neural elasticity graph.
#[derive(Debug, Clone)]
pub struct NeuralElasticEdge {
    /// Index of the source node.
    pub src: usize,
    /// Index of the destination node.
    pub dst: usize,
    /// Rest length of this edge.
    pub rest_length: f64,
}

impl NeuralElasticEdge {
    /// Create a new edge between `src` and `dst` with given `rest_length`.
    pub fn new(src: usize, dst: usize, rest_length: f64) -> Self {
        Self {
            src,
            dst,
            rest_length,
        }
    }
}

/// Graph neural network-style elasticity model operating on a particle mesh.
///
/// Each message-passing step aggregates edge features and updates node states,
/// then uses the updated states to predict nodal forces.
#[derive(Debug, Clone)]
pub struct NeuralElasticityGraph {
    /// Nodes of the graph.
    pub nodes: Vec<NeuralElasticNode>,
    /// Edges of the graph.
    pub edges: Vec<NeuralElasticEdge>,
    /// Edge message MLP weights (row-major `n_msg × (2*n_feat+1)`).
    pub edge_w: Vec<f64>,
    /// Edge message biases.
    pub edge_b: Vec<f64>,
    /// Node update MLP weights (row-major `n_feat × (n_feat + n_msg)`).
    pub node_w: Vec<f64>,
    /// Node update biases.
    pub node_b: Vec<f64>,
    /// Force prediction weights (row-major `3 × n_feat`).
    pub force_w: Vec<f64>,
    /// Force prediction biases.
    pub force_b: Vec<f64>,
    /// Feature dimensionality.
    pub n_feat: usize,
    /// Message dimensionality.
    pub n_msg: usize,
}

impl NeuralElasticityGraph {
    /// Construct a graph with given nodes, edges, feature dim, and message dim.
    ///
    /// All network weights are initialised to zero.
    pub fn new(
        nodes: Vec<NeuralElasticNode>,
        edges: Vec<NeuralElasticEdge>,
        n_feat: usize,
        n_msg: usize,
    ) -> Self {
        let edge_in = 2 * n_feat + 1; // features of both endpoints + edge_length
        Self {
            nodes,
            edges,
            edge_w: vec![0.0; n_msg * edge_in],
            edge_b: vec![0.0; n_msg],
            node_w: vec![0.0; n_feat * (n_feat + n_msg)],
            node_b: vec![0.0; n_feat],
            force_w: vec![0.0; 3 * n_feat],
            force_b: vec![0.0; 3],
            n_feat,
            n_msg,
        }
    }

    /// Compute current edge length between two nodes.
    fn edge_length(&self, e: &NeuralElasticEdge) -> f64 {
        let a = &self.nodes[e.src].position;
        let b = &self.nodes[e.dst].position;
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Run one message-passing step and return predicted nodal forces.
    pub fn predict_forces(&self) -> Vec<[f64; 3]> {
        let n = self.nodes.len();
        let edge_in = 2 * self.n_feat + 1;

        // Aggregate messages for each node
        let mut agg = vec![vec![0.0f64; self.n_msg]; n];
        for e in &self.edges {
            let fi = &self.nodes[e.src].features;
            let fj = &self.nodes[e.dst].features;
            let len = self.edge_length(e);
            let mut inp = Vec::with_capacity(edge_in);
            inp.extend_from_slice(fi);
            inp.extend_from_slice(fj);
            inp.push(len - e.rest_length);

            let msg = dense_forward(&self.edge_w, &self.edge_b, &inp, edge_in, self.n_msg);
            let msg_act: Vec<f64> = msg.iter().map(|&x| relu(x)).collect();
            for k in 0..self.n_msg {
                agg[e.src][k] += msg_act[k];
                agg[e.dst][k] += msg_act[k];
            }
        }

        // Node update and force prediction
        let node_in = self.n_feat + self.n_msg;
        let mut forces = vec![[0.0f64; 3]; n];
        for (i, node) in self.nodes.iter().enumerate() {
            let mut inp = Vec::with_capacity(node_in);
            inp.extend_from_slice(&node.features);
            inp.extend_from_slice(&agg[i]);
            let updated = dense_forward(&self.node_w, &self.node_b, &inp, node_in, self.n_feat);
            let updated_act: Vec<f64> = updated.iter().map(|&x| relu(x)).collect();
            let fout = dense_forward(&self.force_w, &self.force_b, &updated_act, self.n_feat, 3);
            forces[i] = [fout[0], fout[1], fout[2]];
        }
        forces
    }
}

// ---------------------------------------------------------------------------
// TrainingDataset
// ---------------------------------------------------------------------------

/// A single training sample: strain-stress pair from physics simulation.
#[derive(Debug, Clone)]
pub struct StrainStressSample {
    /// Voigt strain vector (length 6).
    pub strain: [f64; 6],
    /// Corresponding stress vector (length 6).
    pub stress: [f64; 6],
}

/// Dataset of strain-stress samples generated from forward physics simulation.
#[derive(Debug, Clone)]
pub struct TrainingDataset {
    /// Collected samples.
    pub samples: Vec<StrainStressSample>,
    /// Maximum dataset size before oldest samples are evicted.
    pub max_size: usize,
}

impl TrainingDataset {
    /// Create an empty dataset with given maximum capacity.
    pub fn new(max_size: usize) -> Self {
        Self {
            samples: Vec::new(),
            max_size,
        }
    }

    /// Add a sample. If the dataset is full, evict the oldest sample.
    ///
    /// A dataset with `max_size == 0` stores nothing.
    pub fn push(&mut self, sample: StrainStressSample) {
        if self.max_size == 0 {
            return;
        }
        if self.samples.len() >= self.max_size {
            self.samples.remove(0);
        }
        self.samples.push(sample);
    }

    /// Generate synthetic training data from a linear-elastic constitutive law.
    ///
    /// Creates `n_samples` strain-stress pairs using the isotropic stiffness
    /// matrix `C` (represented by `E` and `nu`).
    pub fn generate_linear_elastic(
        &mut self,
        young: f64,
        nu: f64,
        n_samples: usize,
        rng: &mut impl rand::Rng,
    ) {
        let mu = young / (2.0 * (1.0 + nu));
        let lam = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        use rand::RngExt as _;
        for _ in 0..n_samples {
            let mut strain = [0.0f64; 6];
            for s in strain.iter_mut() {
                *s = rng.random_range(-0.01f64..0.01f64);
            }
            // Isotropic C * strain in Voigt notation
            let (e11, e22, e33, e12, e13, e23) = (
                strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
            );
            let tr = e11 + e22 + e33;
            let stress = [
                lam * tr + 2.0 * mu * e11,
                lam * tr + 2.0 * mu * e22,
                lam * tr + 2.0 * mu * e33,
                2.0 * mu * e12,
                2.0 * mu * e13,
                2.0 * mu * e23,
            ];
            self.push(StrainStressSample { strain, stress });
        }
    }

    /// Return the number of stored samples.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Return `true` if there are no samples.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
}

// ---------------------------------------------------------------------------
// InferenceTimeCorrector
// ---------------------------------------------------------------------------

/// Inference-time residual corrector: adjusts network-predicted displacements
/// using a fast linear correction step.
///
/// Given a predicted displacement `u_pred` and the actual equilibrium residual
/// `r = K * u_pred - f`, this module applies a gradient step to reduce `||r||`.
#[derive(Debug, Clone)]
pub struct InferenceTimeCorrector {
    /// Step size for the correction.
    pub alpha: f64,
    /// Maximum number of correction iterations.
    pub max_iters: usize,
    /// Convergence tolerance on `||r||`.
    pub tol: f64,
}

impl InferenceTimeCorrector {
    /// Create a new corrector with given step size, max iterations, and tolerance.
    pub fn new(alpha: f64, max_iters: usize, tol: f64) -> Self {
        Self {
            alpha,
            max_iters,
            tol,
        }
    }

    /// Apply iterative correction to `u` given stiffness-matrix-vector product
    /// function `kv(u) = K * u` and external force `f`.
    ///
    /// Returns the number of iterations performed.
    pub fn correct(&self, u: &mut [f64], f: &[f64], kv: &dyn Fn(&[f64]) -> Vec<f64>) -> usize {
        let n = u.len();
        for iter in 0..self.max_iters {
            let ku = kv(u);
            let mut r = vec![0.0f64; n];
            let mut r_norm = 0.0f64;
            for i in 0..n {
                r[i] = ku[i] - f[i];
                r_norm += r[i] * r[i];
            }
            r_norm = r_norm.sqrt();
            if r_norm < self.tol {
                return iter;
            }
            for i in 0..n {
                u[i] -= self.alpha * r[i];
            }
        }
        self.max_iters
    }
}

// ---------------------------------------------------------------------------
// NeuralPoseDeformer
// ---------------------------------------------------------------------------

/// Neural pose-driven deformation: blend-shape skinning plus a learned
/// neural residual correction.
///
/// Given `n_poses` blend-shape vectors of dimension `n_verts * 3` and a
/// weight vector, produces the deformed vertex positions.
#[derive(Debug, Clone)]
pub struct NeuralPoseDeformer {
    /// Number of vertices.
    pub n_verts: usize,
    /// Number of blend shapes.
    pub n_poses: usize,
    /// Blend-shape matrix (row-major `n_verts*3 × n_poses`).
    pub blend_shapes: Vec<f64>,
    /// Neural residual MLP weights (`n_verts*3 × n_poses` → `n_hidden` → `n_verts*3`).
    pub res_w1: Vec<f64>,
    /// Residual MLP hidden biases.
    pub res_b1: Vec<f64>,
    /// Residual MLP output weights.
    pub res_w2: Vec<f64>,
    /// Residual MLP output biases.
    pub res_b2: Vec<f64>,
    /// Hidden layer size for residual MLP.
    pub n_hidden: usize,
}

impl NeuralPoseDeformer {
    /// Create a new deformer for `n_verts` vertices with `n_poses` blend shapes.
    pub fn new(n_verts: usize, n_poses: usize, n_hidden: usize) -> Self {
        let dim = n_verts * 3;
        Self {
            n_verts,
            n_poses,
            blend_shapes: vec![0.0; dim * n_poses],
            res_w1: vec![0.0; n_hidden * n_poses],
            res_b1: vec![0.0; n_hidden],
            res_w2: vec![0.0; dim * n_hidden],
            res_b2: vec![0.0; dim],
            n_hidden,
        }
    }

    /// Compute deformed positions for a given pose-weight vector.
    ///
    /// Result is a flat vector of length `n_verts * 3`.
    pub fn deform(&self, pose_weights: &[f64]) -> Vec<f64> {
        assert_eq!(pose_weights.len(), self.n_poses);
        let dim = self.n_verts * 3;

        // Linear blend-shape contribution
        let blend = dense_forward(
            &self.blend_shapes,
            &vec![0.0; dim],
            pose_weights,
            self.n_poses,
            dim,
        );

        // Neural residual
        let h = dense_forward(
            &self.res_w1,
            &self.res_b1,
            pose_weights,
            self.n_poses,
            self.n_hidden,
        );
        let h_act: Vec<f64> = h.iter().map(|&x| tanh_act(x)).collect();
        let residual = dense_forward(&self.res_w2, &self.res_b2, &h_act, self.n_hidden, dim);

        blend
            .iter()
            .zip(residual.iter())
            .map(|(b, r)| b + r)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// LatentSpaceDeformer
// ---------------------------------------------------------------------------

/// Latent-space deformation: encode vertex positions to a low-dimensional
/// latent code and decode back, enabling compact neural deformation.
#[derive(Debug, Clone)]
pub struct LatentSpaceDeformer {
    /// Number of vertices.
    pub n_verts: usize,
    /// Latent dimensionality.
    pub latent_dim: usize,
    /// Encoder weights (row-major `latent_dim × n_verts*3`).
    pub enc_w: Vec<f64>,
    /// Encoder biases.
    pub enc_b: Vec<f64>,
    /// Decoder weights (row-major `n_verts*3 × latent_dim`).
    pub dec_w: Vec<f64>,
    /// Decoder biases.
    pub dec_b: Vec<f64>,
}

impl LatentSpaceDeformer {
    /// Create a new latent-space deformer with zero-initialised weights.
    pub fn new(n_verts: usize, latent_dim: usize) -> Self {
        let dim = n_verts * 3;
        Self {
            n_verts,
            latent_dim,
            enc_w: vec![0.0; latent_dim * dim],
            enc_b: vec![0.0; latent_dim],
            dec_w: vec![0.0; dim * latent_dim],
            dec_b: vec![0.0; dim],
        }
    }

    /// Encode a flat vertex position vector to latent code.
    pub fn encode(&self, positions: &[f64]) -> Vec<f64> {
        let dim = self.n_verts * 3;
        assert_eq!(positions.len(), dim);
        let pre = dense_forward(&self.enc_w, &self.enc_b, positions, dim, self.latent_dim);
        pre.iter().map(|&x| tanh_act(x)).collect()
    }

    /// Decode a latent code back to flat vertex positions.
    pub fn decode(&self, latent: &[f64]) -> Vec<f64> {
        assert_eq!(latent.len(), self.latent_dim);
        let dim = self.n_verts * 3;
        dense_forward(&self.dec_w, &self.dec_b, latent, self.latent_dim, dim)
    }

    /// Round-trip encode then decode.
    pub fn round_trip(&self, positions: &[f64]) -> Vec<f64> {
        let z = self.encode(positions);
        self.decode(&z)
    }

    /// Reconstruction loss: MSE between `positions` and `round_trip(positions)`.
    pub fn reconstruction_loss(&self, positions: &[f64]) -> f64 {
        let recon = self.round_trip(positions);
        let n = positions.len() as f64;
        positions
            .iter()
            .zip(recon.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            / n
    }
}

// ---------------------------------------------------------------------------
// RealTimeNeuralSoftBody
// ---------------------------------------------------------------------------

/// A vertex in a real-time neural softbody mesh.
#[derive(Debug, Clone)]
pub struct NeuralVertex {
    /// Current position.
    pub position: [f64; 3],
    /// Current velocity.
    pub velocity: [f64; 3],
    /// Mass.
    pub mass: f64,
    /// `true` if the vertex is kinematically fixed.
    pub is_static: bool,
}

impl NeuralVertex {
    /// Create a new dynamic vertex at `position` with given `mass`.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            is_static: false,
        }
    }

    /// Create a static (pinned) vertex.
    pub fn new_static(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass: f64::INFINITY,
            is_static: true,
        }
    }
}

/// Nodal displacements `x - X` of the four corners of `elem`, or `None` if a
/// corner index is out of range.
fn element_displacements(
    positions: &[[f64; 3]],
    rest: &[[f64; 3]],
    elem: &[usize; 4],
) -> Option<[[f64; 3]; 4]> {
    let mut u = [[0.0f64; 3]; 4];
    for (slot, &corner) in u.iter_mut().zip(elem.iter()) {
        let x = positions.get(corner)?;
        let x_rest = rest.get(corner)?;
        *slot = [x[0] - x_rest[0], x[1] - x_rest[1], x[2] - x_rest[2]];
    }
    Some(u)
}

/// Real-time neural softbody: combines a learned constitutive model with
/// a simple explicit Euler integrator for high-speed soft simulation.
///
/// The mesh consists of linear tetrahedra.  The strain of an element is the
/// small strain `eps = sym(grad u)` of its displacement field relative to the
/// rest configuration recorded when the body is built (see [`TetRestShape`]);
/// the learned model maps it to a Voigt stress, and the stress becomes the
/// nodal forces `f_a = -V0 * sigma * grad N_a`.  The forces of an element sum
/// to zero and vanish for a rigid translation.  Being linearised, the strain
/// is not invariant under large rigid rotations.
#[derive(Debug, Clone)]
pub struct RealTimeNeuralSoftBody {
    /// Vertices of the mesh.
    pub vertices: Vec<NeuralVertex>,
    /// Tetrahedral element connectivity (4 vertex indices per element).
    ///
    /// Call [`refresh_rest_shapes`](Self::refresh_rest_shapes) after changing it.
    pub elements: Vec<[usize; 4]>,
    /// Neural constitutive model for stress prediction.
    pub constitutive: NeuralConstitutiveModel,
    /// Mass density.
    pub density: f64,
    /// Rayleigh damping coefficient.
    pub damping: f64,
    /// Rest positions of the vertices.
    rest_positions: Vec<[f64; 3]>,
    /// Per-element rest-state data (`None` for invalid or degenerate elements).
    rest_shapes: Vec<Option<TetRestShape>>,
}

impl RealTimeNeuralSoftBody {
    /// Construct a new real-time neural softbody.
    ///
    /// The current vertex positions are the rest configuration.  The inverse
    /// rest edge matrix and the rest volume of every element are precomputed
    /// here, once.  Elements with an out-of-range vertex index or a degenerate
    /// rest shape exert no force; they are listed by
    /// [`invalid_elements`](Self::invalid_elements).
    pub fn new(
        vertices: Vec<NeuralVertex>,
        elements: Vec<[usize; 4]>,
        constitutive: NeuralConstitutiveModel,
        density: f64,
        damping: f64,
    ) -> Self {
        let rest_positions: Vec<[f64; 3]> = vertices.iter().map(|v| v.position).collect();
        let rest_shapes = Self::build_rest_shapes(&rest_positions, &elements);
        Self {
            vertices,
            elements,
            constitutive,
            density,
            damping,
            rest_positions,
            rest_shapes,
        }
    }

    fn build_rest_shapes(
        rest_positions: &[[f64; 3]],
        elements: &[[usize; 4]],
    ) -> Vec<Option<TetRestShape>> {
        elements
            .iter()
            .map(|elem| TetRestShape::from_rest(rest_positions, elem))
            .collect()
    }

    /// Rest positions of the vertices.
    pub fn rest_positions(&self) -> &[[f64; 3]] {
        &self.rest_positions
    }

    /// Replace the rest configuration and recompute the element rest data.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] unless there is exactly one rest
    /// position per vertex.
    pub fn set_rest_positions(&mut self, rest_positions: Vec<[f64; 3]>) -> Result<()> {
        if rest_positions.len() != self.vertices.len() {
            return Err(Error::InvalidParameter(format!(
                "expected {} rest positions (one per vertex), got {}",
                self.vertices.len(),
                rest_positions.len()
            )));
        }
        self.rest_positions = rest_positions;
        self.refresh_rest_shapes();
        Ok(())
    }

    /// Recompute the per-element rest data from the rest positions, e.g.
    /// after modifying [`elements`](Self::elements).
    pub fn refresh_rest_shapes(&mut self) {
        self.rest_shapes = Self::build_rest_shapes(&self.rest_positions, &self.elements);
    }

    /// Indices of the elements that exert no force because a vertex index is
    /// out of range or the rest tetrahedron is degenerate.
    pub fn invalid_elements(&self) -> Vec<usize> {
        (0..self.elements.len())
            .filter(|&e| self.rest_shapes.get(e).copied().flatten().is_none())
            .collect()
    }

    /// Compute the Voigt strain for a tetrahedral element given its current
    /// and rest positions.
    ///
    /// The strain is the linearised strain tensor `eps = 0.5*(grad u + grad u^T)`
    /// with the displacement gradient `grad u = Du * Dm^-1`, where `Du` holds
    /// the displacement differences `u_i - u_0` and `Dm` the rest edge vectors
    /// `X_i - X_0` of the element.  Returns zero strain for a degenerate rest
    /// tetrahedron or an out-of-range vertex index.
    pub fn element_strain(
        positions: &[[f64; 3]],
        rest: &[[f64; 3]],
        elem: &[usize; 4],
    ) -> [f64; 6] {
        match (
            TetRestShape::from_rest(rest, elem),
            element_displacements(positions, rest, elem),
        ) {
            (Some(shape), Some(u)) => shape.voigt_strain(&u),
            _ => [0.0; 6],
        }
    }

    /// Elastic force on every vertex for the current positions.
    ///
    /// The learned stress of each element is converted to nodal forces
    /// `f_a = -V0 * sigma * grad N_a` (columns of `-V0 sigma Dm^-T`, and
    /// `f_0 = -(f_1 + f_2 + f_3)`) and accumulated per vertex.  Reaction forces
    /// on static vertices are included.
    pub fn elastic_forces(&self) -> Vec<[f64; 3]> {
        let positions: Vec<[f64; 3]> = self.vertices.iter().map(|v| v.position).collect();
        let mut forces = vec![[0.0f64; 3]; self.vertices.len()];
        for (elem, shape) in self.elements.iter().zip(self.rest_shapes.iter()) {
            let Some(shape) = shape else {
                continue;
            };
            let Some(u) = element_displacements(&positions, &self.rest_positions, elem) else {
                continue;
            };
            let strain = shape.voigt_strain(&u);
            let stress = self.constitutive.forward_voigt(&strain);
            let nodal = shape.nodal_forces(&stress);
            for (&vertex, force) in elem.iter().zip(nodal.iter()) {
                if let Some(total) = forces.get_mut(vertex) {
                    for (t, f) in total.iter_mut().zip(force.iter()) {
                        *t += f;
                    }
                }
            }
        }
        forces
    }

    /// Advance the simulation by one time step `dt` given external gravity vector.
    pub fn step(&mut self, dt: f64, gravity: &[f64; 3]) {
        let mut forces = self.elastic_forces();

        // Apply gravity
        for (force, v) in forces.iter_mut().zip(self.vertices.iter()) {
            if !v.is_static {
                for (f, g) in force.iter_mut().zip(gravity.iter()) {
                    *f += v.mass * g;
                }
            }
        }

        // Explicit Euler integration
        for (v, force) in self.vertices.iter_mut().zip(forces.iter()) {
            if v.is_static {
                continue;
            }
            let inv_m = 1.0 / v.mass;
            for (k, &fk) in force.iter().enumerate() {
                let acc = fk * inv_m - self.damping * v.velocity[k];
                v.velocity[k] += dt * acc;
                v.position[k] += dt * v.velocity[k];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// NeuralPinnSolver
// ---------------------------------------------------------------------------

/// Physics-informed neural network solver for the constitutive model.
///
/// Minimises `L = L_data + lambda_pde * L_pde` over the network weights by
/// gradient descent.  `L_data` is the mean squared error between predicted
/// and target stress.  `L_pde` is the equilibrium residual
/// `||div(sigma) + f||^2` of the *predicted* stress field on a regular grid;
/// it takes part in training when the solver's loss has a grid layout
/// ([`PhysicsInformedLoss::grid_dims`]) and the batch holds the samples of
/// the grid nodes in x-fastest order.  Without a grid layout the samples are
/// independent, carry no spatial neighbours, and only the data term is trained.
#[derive(Debug, Clone)]
pub struct NeuralPinnSolver {
    /// The neural constitutive model used to predict stress from strain.
    pub model: NeuralConstitutiveModel,
    /// PINN loss calculator.
    pub pinn_loss: PhysicsInformedLoss,
    /// Learning rate.
    pub lr: f64,
    /// Training iteration count.
    pub iteration: usize,
}

impl NeuralPinnSolver {
    /// Create a new PINN solver with a He-initialised model.
    ///
    /// The loss starts without a grid layout (data term only) and with unit
    /// grid spacing; set [`pinn_loss`](Self::pinn_loss) to
    /// `PhysicsInformedLoss::new(lambda_pde, spacing).with_grid(dims)` to train
    /// on a stress field.
    pub fn new(n_hidden: usize, lambda_pde: f64, lr: f64) -> Self {
        Self {
            model: NeuralConstitutiveModel::new(n_hidden),
            pinn_loss: PhysicsInformedLoss::new(lambda_pde, 1.0),
            lr,
            iteration: 0,
        }
    }

    /// Run one gradient-descent step on a batch of strain-stress samples.
    ///
    /// Returns the total loss `L_data + lambda_pde * L_pde` of the batch
    /// evaluated at the weights *before* the update, with
    /// `L_data = 1/(6 N) sum_k ||sigma_k - target_k||^2`.  The PDE term uses
    /// `body_force` and is only present for a gridded batch (see
    /// [`PhysicsInformedLoss::grid_dims`]); its gradient
    /// `lambda_pde * d L_pde / d sigma` is added to the output gradient of the
    /// data term before backpropagation through the network.
    ///
    /// An empty batch changes nothing and returns 0.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] (leaving the model untouched) if a
    /// grid layout is set but the batch size does not match it, or the grid
    /// is invalid (see [`PhysicsInformedLoss::equilibrium_residuals`]).
    pub fn train_step(
        &mut self,
        batch: &[StrainStressSample],
        body_force: [f64; 3],
    ) -> Result<f64> {
        if batch.is_empty() {
            return Ok(0.0);
        }
        let predictions: Vec<[f64; 6]> = batch
            .iter()
            .map(|s| self.model.forward_voigt(&s.strain))
            .collect();

        // Data term L_data = 1/(6 N) sum ||pred - target||^2 and its gradient.
        let inv_count = 1.0 / (6.0 * batch.len() as f64);
        let mut loss = 0.0;
        let mut grad_out: Vec<[f64; 6]> = Vec::with_capacity(batch.len());
        for (pred, sample) in predictions.iter().zip(batch.iter()) {
            let mut grad = [0.0f64; 6];
            for ((g, p), t) in grad.iter_mut().zip(pred.iter()).zip(sample.stress.iter()) {
                let diff = p - t;
                loss += diff * diff * inv_count;
                *g = 2.0 * diff * inv_count;
            }
            grad_out.push(grad);
        }

        // Equilibrium term of the predicted stress field.
        if let Some(dims) = self.pinn_loss.grid_dims {
            let lambda = self.pinn_loss.lambda_pde;
            let pde_loss = self.pinn_loss.pde_loss(&predictions, dims, body_force)?;
            let pde_grad = self
                .pinn_loss
                .pde_gradient(&predictions, dims, body_force)?;
            loss += lambda * pde_loss;
            for (grad, pde) in grad_out.iter_mut().zip(pde_grad.iter()) {
                for (g, d) in grad.iter_mut().zip(pde.iter()) {
                    *g += lambda * d;
                }
            }
        }

        let mut acc_gw1 = vec![0.0f64; self.model.w1.len()];
        let mut acc_gb1 = vec![0.0f64; self.model.b1.len()];
        let mut acc_gw2 = vec![0.0f64; self.model.w2.len()];
        let mut acc_gb2 = vec![0.0f64; self.model.b2.len()];
        for (sample, grad) in batch.iter().zip(grad_out.iter()) {
            let (gw1, gb1, gw2, gb2) = self.model.backward_from_output_grad(&sample.strain, grad);
            vec_add_inplace(&mut acc_gw1, &gw1);
            vec_add_inplace(&mut acc_gb1, &gb1);
            vec_add_inplace(&mut acc_gw2, &gw2);
            vec_add_inplace(&mut acc_gb2, &gb2);
        }

        self.model
            .sgd_step(&acc_gw1, &acc_gb1, &acc_gw2, &acc_gb2, self.lr);
        self.iteration += 1;
        Ok(loss)
    }
}

// ---------------------------------------------------------------------------
// NeuralSoftBodySimulator  (high-level orchestrator)
// ---------------------------------------------------------------------------

/// High-level simulator that orchestrates training data generation, online
/// training of the neural constitutive model, and real-time stepping.
///
/// The solver's model is the single source of truth for the learned
/// constitutive law: after every training step it is copied into the soft
/// body, so training changes the simulated material.
#[derive(Debug, Clone)]
pub struct NeuralSoftBodySimulator {
    /// The real-time soft body.
    pub body: RealTimeNeuralSoftBody,
    /// The PINN solver used for offline/online training.
    pub solver: NeuralPinnSolver,
    /// Training dataset.
    pub dataset: TrainingDataset,
    /// Inference-time corrector.
    pub corrector: InferenceTimeCorrector,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

impl NeuralSoftBodySimulator {
    /// Construct a new simulator from a vertex/element mesh.
    ///
    /// The initial vertex positions are the rest configuration; the soft body
    /// starts with a copy of the solver's (He-initialised) model.
    pub fn new(
        vertices: Vec<NeuralVertex>,
        elements: Vec<[usize; 4]>,
        n_hidden: usize,
        density: f64,
        damping: f64,
        gravity: [f64; 3],
    ) -> Self {
        let solver = NeuralPinnSolver::new(n_hidden, 0.01, 1e-4);
        let body =
            RealTimeNeuralSoftBody::new(vertices, elements, solver.model.clone(), density, damping);
        let dataset = TrainingDataset::new(1024);
        let corrector = InferenceTimeCorrector::new(0.01, 10, 1e-6);
        Self {
            body,
            solver,
            dataset,
            corrector,
            gravity,
        }
    }

    /// Advance the simulation by `dt`, optionally performing one online training step.
    ///
    /// With `do_train` and a non-empty dataset, the whole dataset is one
    /// training batch (no body force); afterwards the trained model replaces
    /// the soft body's constitutive model and the batch loss is returned.
    /// Otherwise returns `None`.
    ///
    /// # Errors
    ///
    /// Propagates the error of [`NeuralPinnSolver::train_step`] (for example a
    /// grid layout of the solver that does not match the dataset size).
    pub fn step(&mut self, dt: f64, do_train: bool) -> Result<Option<f64>> {
        self.body.step(dt, &self.gravity);
        if !do_train || self.dataset.is_empty() {
            return Ok(None);
        }
        let loss = self.solver.train_step(&self.dataset.samples, [0.0; 3])?;
        self.body.constitutive = self.solver.model.clone();
        Ok(Some(loss))
    }

    /// Rest positions of the vertices.
    pub fn rest_positions(&self) -> &[[f64; 3]] {
        self.body.rest_positions()
    }

    /// Return the current positions of all vertices.
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.body.vertices.iter().map(|v| v.position).collect()
    }
}

// ---------------------------------------------------------------------------
// Utility: simple 3-vector operations
// ---------------------------------------------------------------------------

/// Compute the norm (magnitude) of a 3-vector.
#[inline]
pub fn vec3_norm(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalise a 3-vector. Returns `[0,0,0]` if the norm is less than `1e-12`.
#[inline]
pub fn vec3_normalise(v: &[f64; 3]) -> [f64; 3] {
    let n = vec3_norm(v);
    if n < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Dot product of two 3-vectors.
#[inline]
pub fn vec3_dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn vec3_cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
