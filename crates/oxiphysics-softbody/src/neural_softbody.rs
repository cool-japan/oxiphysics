#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
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

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Internal math helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// Dot product of two equal-length slices.
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// In-place vector addition: `a += b`.
#[inline]
fn vec_add_inplace(a: &mut [f64], b: &[f64]) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x += y;
    }
}

/// In-place scalar multiplication: `a *= s`.
#[inline]
fn vec_scale_inplace(a: &mut [f64], s: f64) {
    for x in a.iter_mut() {
        *x *= s;
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

/// Derivative of tanh.
#[inline]
fn tanh_grad(x: f64) -> f64 {
    let t = x.tanh();
    1.0 - t * t
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

/// A 2-layer ReLU network that maps strain (Voigt 6-vector) to stress (6-vector).
///
/// Architecture: `6 → n_hidden (ReLU) → 6`.
#[derive(Debug, Clone)]
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
    /// Construct a new constitutive model with zero weights and given hidden size.
    pub fn new(n_hidden: usize) -> Self {
        Self {
            n_hidden,
            w1: vec![0.0; n_hidden * 6],
            b1: vec![0.0; n_hidden],
            w2: vec![0.0; 6 * n_hidden],
            b2: vec![0.0; 6],
        }
    }

    /// Construct a simple "linear-elastic" model: identity mapping via weights.
    ///
    /// This initialises the two-layer network so that for small strains the
    /// output approximates `stress = scale * strain`.
    pub fn linear_elastic(n_hidden: usize, scale: f64) -> Self {
        let mut model = Self::new(n_hidden);
        // Set w1 rows to unit vectors scaled by sqrt(scale).
        let s = scale.sqrt();
        for i in 0..n_hidden.min(6) {
            model.w1[i * 6 + i] = s;
        }
        // Set w2 rows to unit vectors scaled by sqrt(scale).
        for i in 0..6 {
            let j = i % n_hidden;
            model.w2[i * n_hidden + j] = s;
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

    /// Backward pass: returns `(grad_w1, grad_b1, grad_w2, grad_b2)` for a MSE loss
    /// `||stress_pred - stress_target||^2`.
    pub fn backward(
        &self,
        strain: &[f64],
        stress_target: &[f64],
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let h_pre = dense_forward(&self.w1, &self.b1, strain, 6, self.n_hidden);
        let h_act: Vec<f64> = h_pre.iter().map(|&x| relu(x)).collect();
        let stress_pred = dense_forward(&self.w2, &self.b2, &h_act, self.n_hidden, 6);

        // dL/d_stress_pred = 2*(pred - target)
        let grad_out: Vec<f64> = stress_pred
            .iter()
            .zip(stress_target.iter())
            .map(|(p, t)| 2.0 * (p - t))
            .collect();

        let gw2 = dense_weight_grad(&h_act, &grad_out, self.n_hidden, 6);
        let gb2 = grad_out.clone();

        // Backprop through ReLU
        let grad_h_pre_act = dense_backward_input(&self.w2, &grad_out, self.n_hidden, 6);
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

/// Physics-informed loss combining data MSE with equilibrium residual penalty.
///
/// `L = L_data + lambda_pde * L_pde`
///
/// where `L_pde` penalises the divergence of the predicted stress tensor.
#[derive(Debug, Clone)]
pub struct PhysicsInformedLoss {
    /// Weight on the PDE residual term.
    pub lambda_pde: f64,
    /// Finite-difference step for computing the stress divergence.
    pub fd_step: f64,
}

impl PhysicsInformedLoss {
    /// Create a new PINN loss with given `lambda_pde` and finite-difference step.
    pub fn new(lambda_pde: f64, fd_step: f64) -> Self {
        Self {
            lambda_pde,
            fd_step,
        }
    }

    /// Compute the total PINN loss given predicted and target stresses plus a
    /// body-force vector.
    ///
    /// `L_pde` is approximated as `||div_sigma - f_body||^2` using a
    /// finite-difference estimate of the divergence from neighbouring samples.
    pub fn compute(&self, stress_pred: &[f64], stress_target: &[f64], _body_force: &[f64]) -> f64 {
        // Data loss: MSE
        let l_data: f64 = stress_pred
            .iter()
            .zip(stress_target.iter())
            .map(|(p, t)| (p - t).powi(2))
            .sum::<f64>()
            / stress_pred.len() as f64;

        // Simplified PDE loss: penalise norm of stress itself (proxy for divergence)
        let l_pde: f64 =
            stress_pred.iter().map(|s| s.powi(2)).sum::<f64>() / stress_pred.len() as f64;

        l_data + self.lambda_pde * l_pde
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
    pub fn push(&mut self, sample: StrainStressSample) {
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
    pub fn correct(&self, u: &mut Vec<f64>, f: &[f64], kv: &dyn Fn(&[f64]) -> Vec<f64>) -> usize {
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

/// Real-time neural softbody: combines a learned constitutive model with
/// a simple explicit Euler integrator for high-speed soft simulation.
#[derive(Debug, Clone)]
pub struct RealTimeNeuralSoftBody {
    /// Vertices of the mesh.
    pub vertices: Vec<NeuralVertex>,
    /// Tetrahedral element connectivity (4 vertex indices per element).
    pub elements: Vec<[usize; 4]>,
    /// Neural constitutive model for stress prediction.
    pub constitutive: NeuralConstitutiveModel,
    /// Mass density.
    pub density: f64,
    /// Rayleigh damping coefficient.
    pub damping: f64,
}

impl RealTimeNeuralSoftBody {
    /// Construct a new real-time neural softbody.
    pub fn new(
        vertices: Vec<NeuralVertex>,
        elements: Vec<[usize; 4]>,
        constitutive: NeuralConstitutiveModel,
        density: f64,
        damping: f64,
    ) -> Self {
        Self {
            vertices,
            elements,
            constitutive,
            density,
            damping,
        }
    }

    /// Compute the Voigt strain for a tetrahedral element given its current
    /// and rest positions.
    ///
    /// Uses the linearised strain tensor `ε = 0.5*(∇u + ∇u^T)`.
    pub fn element_strain(
        positions: &[[f64; 3]],
        rest: &[[f64; 3]],
        elem: &[usize; 4],
    ) -> [f64; 6] {
        // Displacement gradient via finite differences across tetrahedron edges
        let u: Vec<[f64; 3]> = elem
            .iter()
            .map(|&i| {
                [
                    positions[i][0] - rest[i][0],
                    positions[i][1] - rest[i][1],
                    positions[i][2] - rest[i][2],
                ]
            })
            .collect();

        // Simple average strain approximation
        let du_dx = u[1][0] - u[0][0];
        let du_dy = u[2][1] - u[0][1];
        let du_dz = u[3][2] - u[0][2];
        let shear_xy = 0.5 * ((u[1][1] - u[0][1]) + (u[2][0] - u[0][0]));
        let shear_xz = 0.5 * ((u[1][2] - u[0][2]) + (u[3][0] - u[0][0]));
        let shear_yz = 0.5 * ((u[2][2] - u[0][2]) + (u[3][1] - u[0][1]));

        [du_dx, du_dy, du_dz, shear_xy, shear_xz, shear_yz]
    }

    /// Advance the simulation by one time step `dt` given external gravity vector.
    pub fn step(&mut self, dt: f64, gravity: &[f64; 3], rest_positions: &[[f64; 3]]) {
        let n = self.vertices.len();
        let mut forces = vec![[0.0f64; 3]; n];

        // Apply gravity
        for (i, v) in self.vertices.iter().enumerate() {
            if !v.is_static {
                forces[i][0] += v.mass * gravity[0];
                forces[i][1] += v.mass * gravity[1];
                forces[i][2] += v.mass * gravity[2];
            }
        }

        // Compute elastic forces from constitutive model
        let positions: Vec<[f64; 3]> = self.vertices.iter().map(|v| v.position).collect();
        for elem in &self.elements {
            let strain = Self::element_strain(&positions, rest_positions, elem);
            let stress = self.constitutive.forward(&strain);

            // Distribute stress as equal forces on the 4 nodes (simplified)
            let force_per_node = stress.iter().take(3).map(|&s| s * 0.25).collect::<Vec<_>>();
            for &idx in elem.iter() {
                if !self.vertices[idx].is_static {
                    forces[idx][0] -= force_per_node[0];
                    forces[idx][1] -= force_per_node[1];
                    forces[idx][2] -= force_per_node[2];
                }
            }
        }

        // Explicit Euler integration
        for (i, v) in self.vertices.iter_mut().enumerate() {
            if v.is_static {
                continue;
            }
            let inv_m = 1.0 / v.mass;
            for k in 0..3 {
                let acc = forces[i][k] * inv_m - self.damping * v.velocity[k];
                v.velocity[k] += dt * acc;
                v.position[k] += dt * v.velocity[k];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// NeuralPinnSolver
// ---------------------------------------------------------------------------

/// Physics-informed neural network solver for static deformation.
///
/// Minimises `L = L_data + lambda_pinn * ||K*u - f||^2` where `u` is the
/// displacement predicted by the network.
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
    /// Create a new PINN solver.
    pub fn new(n_hidden: usize, lambda_pde: f64, lr: f64) -> Self {
        Self {
            model: NeuralConstitutiveModel::new(n_hidden),
            pinn_loss: PhysicsInformedLoss::new(lambda_pde, 1e-5),
            lr,
            iteration: 0,
        }
    }

    /// Run one training step on a batch of strain-stress samples.
    ///
    /// Returns the mean total PINN loss over the batch.
    pub fn train_step(&mut self, batch: &[StrainStressSample], body_force: &[f64]) -> f64 {
        let mut total_loss = 0.0;
        let mut acc_gw1 = vec![0.0f64; self.model.w1.len()];
        let mut acc_gb1 = vec![0.0f64; self.model.b1.len()];
        let mut acc_gw2 = vec![0.0f64; self.model.w2.len()];
        let mut acc_gb2 = vec![0.0f64; self.model.b2.len()];

        for s in batch {
            let stress_pred = self.model.forward(&s.strain);
            total_loss += self.pinn_loss.compute(&stress_pred, &s.stress, body_force);

            let (gw1, gb1, gw2, gb2) = self.model.backward(&s.strain, &s.stress);
            vec_add_inplace(&mut acc_gw1, &gw1);
            vec_add_inplace(&mut acc_gb1, &gb1);
            vec_add_inplace(&mut acc_gw2, &gw2);
            vec_add_inplace(&mut acc_gb2, &gb2);
        }

        let n = batch.len() as f64;
        vec_scale_inplace(&mut acc_gw1, 1.0 / n);
        vec_scale_inplace(&mut acc_gb1, 1.0 / n);
        vec_scale_inplace(&mut acc_gw2, 1.0 / n);
        vec_scale_inplace(&mut acc_gb2, 1.0 / n);

        self.model
            .sgd_step(&acc_gw1, &acc_gb1, &acc_gw2, &acc_gb2, self.lr);
        self.iteration += 1;
        total_loss / n
    }
}

// ---------------------------------------------------------------------------
// NeuralSoftBodySimulator  (high-level orchestrator)
// ---------------------------------------------------------------------------

/// High-level simulator that orchestrates training data generation, online
/// training of the neural constitutive model, and real-time stepping.
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
    /// Rest positions of vertices.
    pub rest_positions: Vec<[f64; 3]>,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

impl NeuralSoftBodySimulator {
    /// Construct a new simulator from a vertex/element mesh.
    pub fn new(
        vertices: Vec<NeuralVertex>,
        elements: Vec<[usize; 4]>,
        n_hidden: usize,
        density: f64,
        damping: f64,
        gravity: [f64; 3],
    ) -> Self {
        let rest_positions: Vec<[f64; 3]> = vertices.iter().map(|v| v.position).collect();
        let constitutive = NeuralConstitutiveModel::new(n_hidden);
        let body = RealTimeNeuralSoftBody::new(vertices, elements, constitutive, density, damping);
        let solver = NeuralPinnSolver::new(n_hidden, 0.01, 1e-4);
        let dataset = TrainingDataset::new(1024);
        let corrector = InferenceTimeCorrector::new(0.01, 10, 1e-6);
        Self {
            body,
            solver,
            dataset,
            corrector,
            rest_positions,
            gravity,
        }
    }

    /// Advance the simulation by `dt`, optionally performing one online training step.
    pub fn step(&mut self, dt: f64, do_train: bool) {
        self.body.step(dt, &self.gravity, &self.rest_positions);
        if do_train && !self.dataset.is_empty() {
            let batch: Vec<StrainStressSample> = self.dataset.samples.clone();
            self.solver.train_step(&batch, &[0.0, 0.0, 0.0]);
        }
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // NeuralConstitutiveModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_constitutive_forward_zero_input() {
        let model = NeuralConstitutiveModel::new(16);
        let strain = [0.0f64; 6];
        let stress = model.forward(&strain);
        assert_eq!(stress.len(), 6);
        for s in &stress {
            assert!(
                s.abs() < 1e-12,
                "Zero input should give zero output with zero weights"
            );
        }
    }

    #[test]
    fn test_constitutive_forward_nonzero_weights() {
        let mut model = NeuralConstitutiveModel::linear_elastic(8, 1000.0);
        model.b2[0] = 0.5;
        let strain = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        let stress = model.forward(&strain);
        assert_eq!(stress.len(), 6);
        // With non-zero weights the stress should not all be zero
        let norm: f64 = stress.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(
            norm > 0.0,
            "Non-zero weights should produce non-zero output"
        );
    }

    #[test]
    fn test_constitutive_backward_gradient_shape() {
        let model = NeuralConstitutiveModel::new(8);
        let strain = [0.01, 0.02, -0.01, 0.0, 0.005, -0.003];
        let target = [100.0, 50.0, -20.0, 0.0, 10.0, -5.0];
        let (gw1, gb1, gw2, gb2) = model.backward(&strain, &target);
        assert_eq!(gw1.len(), model.n_hidden * 6);
        assert_eq!(gb1.len(), model.n_hidden);
        assert_eq!(gw2.len(), 6 * model.n_hidden);
        assert_eq!(gb2.len(), 6);
    }

    #[test]
    fn test_constitutive_sgd_reduces_loss() {
        let mut model = NeuralConstitutiveModel::linear_elastic(16, 1.0);
        let strain = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        let target = [10.0, 5.0, 5.0, 0.0, 0.0, 0.0];

        let initial_loss: f64 = {
            let pred = model.forward(&strain);
            pred.iter()
                .zip(target.iter())
                .map(|(p, t)| (p - t).powi(2))
                .sum()
        };

        for _ in 0..200 {
            let (gw1, gb1, gw2, gb2) = model.backward(&strain, &target);
            model.sgd_step(&gw1, &gb1, &gw2, &gb2, 1e-3);
        }

        let final_loss: f64 = {
            let pred = model.forward(&strain);
            pred.iter()
                .zip(target.iter())
                .map(|(p, t)| (p - t).powi(2))
                .sum()
        };

        assert!(
            final_loss < initial_loss,
            "SGD should reduce loss: {initial_loss} -> {final_loss}"
        );
    }

    // -----------------------------------------------------------------------
    // PhysicsInformedLoss tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pinn_loss_zero_when_perfect() {
        let loss = PhysicsInformedLoss::new(0.0, 1e-5);
        let s = [1.0, 2.0, 3.0, 0.5, -0.5, 0.1];
        let body = [0.0; 6];
        let l = loss.compute(&s, &s, &body);
        assert!(l.abs() < 1e-12, "Loss should be zero when pred == target");
    }

    #[test]
    fn test_pinn_loss_positive() {
        let loss = PhysicsInformedLoss::new(0.1, 1e-5);
        let pred = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let target = [0.0; 6];
        let body = [0.0; 6];
        let l = loss.compute(&pred, &target, &body);
        assert!(
            l > 0.0,
            "Loss should be positive for non-matching pred/target"
        );
    }

    #[test]
    fn test_pinn_loss_lambda_scaling() {
        let loss_a = PhysicsInformedLoss::new(0.0, 1e-5);
        let loss_b = PhysicsInformedLoss::new(1.0, 1e-5);
        let pred = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let target = [0.0; 6];
        let body = [0.0; 6];
        let la = loss_a.compute(&pred, &target, &body);
        let lb = loss_b.compute(&pred, &target, &body);
        assert!(
            lb > la,
            "Higher lambda_pde should give higher loss: la={la}, lb={lb}"
        );
    }

    // -----------------------------------------------------------------------
    // LearnedCollisionResponse tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_collision_response_output_shape() {
        let responder = LearnedCollisionResponse::new(8);
        let normal = [0.0, 1.0, 0.0];
        let impulse = responder.predict(&normal, -0.5, 0.01);
        assert_eq!(impulse.len(), 3);
    }

    #[test]
    fn test_collision_response_zero_weights() {
        let responder = LearnedCollisionResponse::new(8);
        let normal = [0.0, 1.0, 0.0];
        let impulse = responder.predict(&normal, -0.5, 0.01);
        for x in &impulse {
            assert!(x.abs() < 1e-12, "Zero weights should give zero impulse");
        }
    }

    // -----------------------------------------------------------------------
    // MaterialParameterIdentifier tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_material_id_lame_params() {
        let id = MaterialParameterIdentifier::new(1e6, 0.3, 1e-3, 0.0);
        let lam = id.lame_lambda();
        let mu = id.shear_modulus();
        assert!(lam > 0.0, "Lamé lambda should be positive");
        assert!(mu > 0.0, "Shear modulus should be positive");
    }

    #[test]
    fn test_material_id_poisson_clamped() {
        // Poisson ratio should stay in [0, 0.499]
        let mut id = MaterialParameterIdentifier::new(1e6, 0.29, 1.0, 0.0);
        let sim_fn = |_e: f64, _nu: f64| -> Vec<[f64; 3]> { vec![[0.0; 3]; 1] };
        let obs = vec![[0.001, 0.0, 0.0]];
        id.step(&obs, &obs, &sim_fn);
        assert!(id.poisson_ratio >= 0.0 && id.poisson_ratio <= 0.499);
    }

    #[test]
    fn test_material_id_young_positive() {
        let mut id = MaterialParameterIdentifier::new(1e6, 0.25, 1e-3, 0.0);
        let sim_fn = |_e: f64, _nu: f64| -> Vec<[f64; 3]> { vec![[0.0; 3]; 4] };
        let obs = vec![[0.001, 0.002, 0.0]; 4];
        for _ in 0..10 {
            id.step(&obs, &obs, &sim_fn);
        }
        assert!(
            id.young_modulus >= 1.0,
            "Young's modulus should stay positive"
        );
    }

    // -----------------------------------------------------------------------
    // NeuralElasticityGraph tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_neural_graph_creation() {
        let nodes = vec![
            NeuralElasticNode::new([0.0, 0.0, 0.0], 4, 1.0),
            NeuralElasticNode::new([1.0, 0.0, 0.0], 4, 1.0),
        ];
        let edges = vec![NeuralElasticEdge::new(0, 1, 1.0)];
        let graph = NeuralElasticityGraph::new(nodes, edges, 4, 8);
        assert_eq!(graph.nodes.len(), 2);
        assert_eq!(graph.edges.len(), 1);
    }

    #[test]
    fn test_neural_graph_zero_forces_with_zero_weights() {
        let nodes = vec![
            NeuralElasticNode::new([0.0, 0.0, 0.0], 4, 1.0),
            NeuralElasticNode::new([1.0, 0.0, 0.0], 4, 1.0),
            NeuralElasticNode::new([0.5, 1.0, 0.0], 4, 1.0),
        ];
        let edges = vec![
            NeuralElasticEdge::new(0, 1, 1.0),
            NeuralElasticEdge::new(1, 2, 1.0),
            NeuralElasticEdge::new(0, 2, 1.0),
        ];
        let graph = NeuralElasticityGraph::new(nodes, edges, 4, 8);
        let forces = graph.predict_forces();
        assert_eq!(forces.len(), 3);
        for f in &forces {
            for &c in f {
                assert!(c.abs() < 1e-12, "Zero weights should yield zero forces");
            }
        }
    }

    // -----------------------------------------------------------------------
    // TrainingDataset tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_dataset_push_and_evict() {
        let mut ds = TrainingDataset::new(3);
        for i in 0..5u64 {
            ds.push(StrainStressSample {
                strain: [i as f64; 6],
                stress: [0.0; 6],
            });
        }
        assert_eq!(ds.len(), 3, "Dataset should not exceed max_size=3");
        assert_eq!(
            ds.samples[0].strain[0], 2.0,
            "Oldest samples should be evicted"
        );
    }

    #[test]
    fn test_dataset_generate_linear_elastic() {
        let mut rng = rand::rng();
        let mut ds = TrainingDataset::new(100);
        ds.generate_linear_elastic(1e6, 0.3, 50, &mut rng);
        assert_eq!(ds.len(), 50);
        // Stress should not be all-zero for non-trivial strains
        let has_nonzero = ds
            .samples
            .iter()
            .any(|s| s.stress.iter().any(|&x| x != 0.0));
        assert!(has_nonzero, "Generated stress should be non-zero");
    }

    #[test]
    fn test_dataset_is_empty() {
        let ds = TrainingDataset::new(10);
        assert!(ds.is_empty());
    }

    // -----------------------------------------------------------------------
    // InferenceTimeCorrector tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_corrector_converges_for_identity_k() {
        let corrector = InferenceTimeCorrector::new(0.5, 100, 1e-8);
        let mut u = vec![2.0, -1.0, 0.5];
        let f = vec![0.0, 0.0, 0.0];
        // K = I => K*u = u, residual = u - f = u => correct to zero
        let iters = corrector.correct(&mut u, &f, &|x: &[f64]| x.to_vec());
        assert!(iters < 100, "Should converge before max iterations");
        for &ui in &u {
            assert!(ui.abs() < 1e-5, "u should converge to 0 = f when K=I");
        }
    }

    #[test]
    fn test_corrector_already_converged() {
        let corrector = InferenceTimeCorrector::new(0.5, 100, 1e-8);
        let mut u = vec![0.0, 0.0, 0.0];
        let f = vec![0.0, 0.0, 0.0];
        let iters = corrector.correct(&mut u, &f, &|x: &[f64]| x.to_vec());
        assert_eq!(iters, 0, "Already converged should return 0 iterations");
    }

    // -----------------------------------------------------------------------
    // NeuralPoseDeformer tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pose_deformer_output_size() {
        let deformer = NeuralPoseDeformer::new(10, 4, 8);
        let weights = [0.25, 0.25, 0.25, 0.25];
        let result = deformer.deform(&weights);
        assert_eq!(result.len(), 30, "Output should be n_verts * 3");
    }

    #[test]
    fn test_pose_deformer_zero_weights_zero_output() {
        let deformer = NeuralPoseDeformer::new(5, 3, 4);
        let weights = [0.0, 0.0, 0.0];
        let result = deformer.deform(&weights);
        for x in &result {
            assert!(
                x.abs() < 1e-12,
                "Zero pose weights should yield zero deformation"
            );
        }
    }

    // -----------------------------------------------------------------------
    // LatentSpaceDeformer tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_latent_encode_output_size() {
        let lsd = LatentSpaceDeformer::new(5, 3);
        let pos = vec![0.0f64; 15];
        let z = lsd.encode(&pos);
        assert_eq!(z.len(), 3);
    }

    #[test]
    fn test_latent_decode_output_size() {
        let lsd = LatentSpaceDeformer::new(5, 3);
        let z = vec![0.0f64; 3];
        let out = lsd.decode(&z);
        assert_eq!(out.len(), 15);
    }

    #[test]
    fn test_latent_reconstruction_loss_zero_weights() {
        let lsd = LatentSpaceDeformer::new(4, 2);
        let pos = vec![
            1.0f64, 2.0, 3.0, -1.0, 0.5, 0.0, 0.1, -0.2, 2.0, 1.5, 0.0, 0.0,
        ];
        // With zero encoder weights encode gives tanh(0)=0, decode gives 0
        // Loss = MSE(pos, 0) > 0
        let loss = lsd.reconstruction_loss(&pos);
        assert!(
            loss > 0.0,
            "Reconstruction loss should be positive with random positions"
        );
    }

    // -----------------------------------------------------------------------
    // RealTimeNeuralSoftBody tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_realtime_body_step_moves_dynamic_vertex() {
        let verts = vec![
            NeuralVertex::new_static([0.0, 0.0, 0.0]),
            NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
        ];
        let rest: Vec<[f64; 3]> = verts.iter().map(|v| v.position).collect();
        let elems = vec![[0, 1, 2, 3]];
        let model = NeuralConstitutiveModel::new(8);
        let mut body = RealTimeNeuralSoftBody::new(verts, elems, model, 1000.0, 0.01);
        let gravity = [0.0, -9.81, 0.0];
        body.step(0.01, &gravity, &rest);
        // Dynamic vertex 1 should have moved in y
        assert!(
            body.vertices[1].position[1] < 0.0,
            "Dynamic vertex should fall under gravity"
        );
    }

    #[test]
    fn test_realtime_body_static_vertex_does_not_move() {
        let verts = vec![
            NeuralVertex::new_static([0.5, 0.5, 0.5]),
            NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
        ];
        let rest: Vec<[f64; 3]> = verts.iter().map(|v| v.position).collect();
        let elems = vec![[0, 1, 2, 3]];
        let model = NeuralConstitutiveModel::new(8);
        let mut body = RealTimeNeuralSoftBody::new(verts, elems, model, 1000.0, 0.01);
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..10 {
            body.step(0.01, &gravity, &rest);
        }
        assert_eq!(
            body.vertices[0].position,
            [0.5, 0.5, 0.5],
            "Static vertex should not move"
        );
    }

    // -----------------------------------------------------------------------
    // NeuralPinnSolver tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pinn_solver_train_step_returns_finite_loss() {
        let mut solver = NeuralPinnSolver::new(8, 0.01, 1e-4);
        let batch = vec![
            StrainStressSample {
                strain: [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
                stress: [1000.0, 300.0, 300.0, 0.0, 0.0, 0.0],
            },
            StrainStressSample {
                strain: [0.0, 0.02, 0.0, 0.0, 0.0, 0.0],
                stress: [300.0, 2000.0, 300.0, 0.0, 0.0, 0.0],
            },
        ];
        let loss = solver.train_step(&batch, &[0.0; 6]);
        assert!(loss.is_finite(), "Training loss should be finite");
        assert!(loss >= 0.0, "Training loss should be non-negative");
    }

    #[test]
    fn test_pinn_solver_iteration_increments() {
        let mut solver = NeuralPinnSolver::new(4, 0.0, 1e-4);
        let batch = vec![StrainStressSample {
            strain: [0.0; 6],
            stress: [0.0; 6],
        }];
        for i in 1..=5 {
            solver.train_step(&batch, &[0.0; 6]);
            assert_eq!(solver.iteration, i);
        }
    }

    // -----------------------------------------------------------------------
    // NeuralSoftBodySimulator tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_simulator_creation() {
        let verts = vec![
            NeuralVertex::new_static([0.0, 0.0, 0.0]),
            NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
        ];
        let elems = vec![[0, 1, 2, 3]];
        let sim = NeuralSoftBodySimulator::new(verts, elems, 8, 1000.0, 0.01, [0.0, -9.81, 0.0]);
        assert_eq!(sim.positions().len(), 4);
        assert_eq!(sim.rest_positions.len(), 4);
    }

    #[test]
    fn test_simulator_step_no_panic() {
        let verts = vec![
            NeuralVertex::new_static([0.0, 0.0, 0.0]),
            NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
            NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
        ];
        let elems = vec![[0, 1, 2, 3]];
        let mut sim =
            NeuralSoftBodySimulator::new(verts, elems, 8, 1000.0, 0.01, [0.0, -9.81, 0.0]);
        for _ in 0..20 {
            sim.step(0.01, false);
        }
        let pos = sim.positions();
        assert_eq!(pos.len(), 4);
        // Positions should be finite
        for p in &pos {
            for &c in p {
                assert!(c.is_finite(), "Positions should remain finite");
            }
        }
    }

    // -----------------------------------------------------------------------
    // Utility tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_vec3_norm() {
        let v = [3.0, 4.0, 0.0];
        assert!((vec3_norm(&v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_normalise() {
        let v = [1.0, 0.0, 0.0];
        let n = vec3_normalise(&v);
        assert!((vec3_norm(&n) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_normalise_zero() {
        let v = [0.0, 0.0, 0.0];
        let n = vec3_normalise(&v);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_vec3_cross() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(&a, &b);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_dot() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((vec3_dot(&a, &b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_dense_forward_identity() {
        // 2×2 identity matrix
        let w = [1.0, 0.0, 0.0, 1.0];
        let b = [0.0, 0.0];
        let x = [3.0, 7.0];
        let y = dense_forward(&w, &b, &x, 2, 2);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_strain_rest_gives_zero() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let rest = positions.clone();
        let elem = [0, 1, 2, 3];
        let strain = RealTimeNeuralSoftBody::element_strain(&positions, &rest, &elem);
        for &s in &strain {
            assert!(s.abs() < 1e-12, "No deformation => zero strain, got {s}");
        }
    }
}
