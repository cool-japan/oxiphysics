//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{angular_fourier_basis, triplet_cos_angle};
use crate::ml_potential::{MlPotential, SymmetryFunction, SymmetryFunctionSet};
use oxiphysics_core::math::Vec3;

/// Ensemble of feedforward neural network potentials.
///
/// Uses multiple independently-trained models to estimate mean energy and
/// epistemic uncertainty (standard deviation).
#[derive(Debug, Clone)]
pub struct EnsembleNNP {
    /// The individual NNP models.
    pub models: Vec<FeedForwardPotential>,
}
impl EnsembleNNP {
    /// Create an ensemble from a list of models.
    pub fn new(models: Vec<FeedForwardPotential>) -> Self {
        assert!(!models.is_empty(), "ensemble must have at least one model");
        Self { models }
    }
    /// Compute the mean energy prediction across the ensemble.
    pub fn mean_energy(&self, positions: &[Vec3], species: &[u32]) -> f64 {
        let n = self.models.len() as f64;
        let sum: f64 = self
            .models
            .iter()
            .map(|m| {
                let (e, _) = m.energy_and_forces(positions, species);
                e
            })
            .sum();
        sum / n
    }
    /// Compute the uncertainty (standard deviation) of energy predictions.
    pub fn uncertainty(&self, positions: &[Vec3], species: &[u32]) -> f64 {
        let energies: Vec<f64> = self
            .models
            .iter()
            .map(|m| {
                let (e, _) = m.energy_and_forces(positions, species);
                e
            })
            .collect();
        let n = energies.len() as f64;
        let mean = energies.iter().sum::<f64>() / n;
        let var = energies
            .iter()
            .map(|&e| (e - mean) * (e - mean))
            .sum::<f64>()
            / n;
        var.sqrt()
    }
    /// Number of models in the ensemble.
    pub fn n_models(&self) -> usize {
        self.models.len()
    }
}
/// A simple multi-layer perceptron interatomic potential.
///
/// Descriptor → hidden layers → single scalar (atomic energy).
/// Total energy = sum of atomic energies.
/// Forces are computed via central finite differences (proof of concept).
#[derive(Debug, Clone)]
pub struct FeedForwardPotential {
    /// Network layers.
    pub layers: Vec<DenseLayer>,
    /// Symmetry function set used to build input descriptors.
    pub symmetry_functions: SymmetryFunctionSet,
    /// Cutoff radius for the descriptor computation.
    pub cutoff: f64,
    /// Finite-difference step size for force computation.
    pub fd_step: f64,
}
impl FeedForwardPotential {
    /// Create a new feedforward potential with random (deterministic) weights.
    ///
    /// `layer_sizes` describes the full architecture including input and
    /// output dimensions, e.g. `[8, 16, 16, 1]` for 8 descriptor inputs,
    /// two hidden layers of 16, and one energy output.
    ///
    /// The symmetry function set and cutoff must be provided separately.
    pub fn new(
        layer_sizes: &[usize],
        activation: Activation,
        symmetry_functions: SymmetryFunctionSet,
        cutoff: f64,
    ) -> Self {
        assert!(
            layer_sizes.len() >= 2,
            "Need at least input and output sizes"
        );
        let mut rng = Lcg::new(42);
        let mut layers = Vec::with_capacity(layer_sizes.len() - 1);
        for l in 0..layer_sizes.len() - 1 {
            let n_in = layer_sizes[l];
            let n_out = layer_sizes[l + 1];
            let scale = (2.0 / (n_in + n_out) as f64).sqrt();
            let mut weights = vec![vec![0.0; n_in]; n_out];
            for row in &mut weights {
                for v in row.iter_mut() {
                    *v = rng.next_f64(scale);
                }
            }
            let biases = vec![0.0; n_out];
            let act = if l == layer_sizes.len() - 2 {
                Activation::Linear
            } else {
                activation
            };
            layers.push(DenseLayer::new(weights, biases, act));
        }
        Self {
            layers,
            symmetry_functions,
            cutoff,
            fd_step: 1e-5,
        }
    }
    /// Evaluate the neural network on a descriptor vector, returning the
    /// scalar atomic energy.
    fn evaluate_network(&self, descriptor: &[f64]) -> f64 {
        let mut x = descriptor.to_vec();
        for layer in &self.layers {
            x = layer.forward(&x);
        }
        x.iter().sum()
    }
    /// Compute total energy for the given positions and species.
    pub(super) fn total_energy(&self, positions: &[Vec3], species: &[u32]) -> f64 {
        let mut energy = 0.0;
        for i in 0..positions.len() {
            let desc =
                self.symmetry_functions
                    .compute_descriptor(positions, i, species, self.cutoff);
            energy += self.evaluate_network(&desc);
        }
        energy
    }
    /// Compute the gradient of the descriptor with respect to atom positions.
    ///
    /// Uses central finite differences to estimate ∂G_k / ∂r_{i,α} for
    /// each descriptor component `k`, atom index `i`, and Cartesian component
    /// `α ∈ {0, 1, 2}`.
    ///
    /// The result is a flattened matrix of shape `[N_desc × N_atoms × 3]`
    /// stored as `grad[k * N_atoms * 3 + i * 3 + alpha]`.
    ///
    /// # Arguments
    /// * `positions` – Cartesian positions (Vec3) of all atoms.
    /// * `center`    – index of the atom whose descriptor is differentiated.
    /// * `species`   – species labels for all atoms.
    /// * `h`         – finite-difference step size (default 1e-5).
    ///
    /// # Returns
    /// Gradient tensor flattened to `Vec`f64` of length
    /// `n_descriptors × N_atoms × 3`.
    pub fn compute_descriptor_gradient(
        &self,
        positions: &[Vec3],
        center: usize,
        species: &[u32],
        h: f64,
    ) -> Vec<f64> {
        let n = positions.len();
        let desc0 =
            self.symmetry_functions
                .compute_descriptor(positions, center, species, self.cutoff);
        let n_desc = desc0.len();
        let mut grad = vec![0.0_f64; n_desc * n * 3];
        let mut pos_pert = positions.to_vec();
        for i in 0..n {
            for alpha in 0..3_usize {
                pos_pert[i][alpha] += h;
                let d_plus = self.symmetry_functions.compute_descriptor(
                    &pos_pert,
                    center,
                    species,
                    self.cutoff,
                );
                pos_pert[i][alpha] -= 2.0 * h;
                let d_minus = self.symmetry_functions.compute_descriptor(
                    &pos_pert,
                    center,
                    species,
                    self.cutoff,
                );
                pos_pert[i][alpha] += h;
                for k in 0..n_desc {
                    grad[k * n * 3 + i * 3 + alpha] = (d_plus[k] - d_minus[k]) / (2.0 * h);
                }
            }
        }
        grad
    }
    /// Build a default symmetry function set with a few G2 functions.
    ///
    /// This is a convenience helper for testing and prototyping.
    pub fn default_symmetry_set(cutoff: f64) -> SymmetryFunctionSet {
        let mut set = SymmetryFunctionSet::new();
        for &eta in &[0.1, 0.5, 1.0, 2.0] {
            set.push(SymmetryFunction::G2 {
                eta,
                rs: 0.0,
                rc: cutoff,
            });
        }
        set
    }
}
/// A single message passing step: aggregates neighbour embeddings.
///
/// Each node `i` collects messages from all neighbours `j` within cutoff and
/// updates its hidden state:
/// ```text
/// h_i^(t+1) = activation( W_self * h_i^t + sum_{j ∈ N(i)} W_msg * h_j^t + b )
/// ```
#[derive(Debug, Clone)]
pub struct MpnnLayer {
    /// Weight matrix for self (n_hidden × n_hidden).
    pub w_self: Vec<Vec<f64>>,
    /// Weight matrix for messages from neighbours (n_hidden × n_hidden).
    pub w_msg: Vec<Vec<f64>>,
    /// Bias vector (n_hidden).
    pub bias: Vec<f64>,
    /// Activation applied after aggregation.
    pub activation: Activation,
}
impl MpnnLayer {
    /// Create a new MPNN layer with the given weight matrices.
    pub fn new(
        w_self: Vec<Vec<f64>>,
        w_msg: Vec<Vec<f64>>,
        bias: Vec<f64>,
        activation: Activation,
    ) -> Self {
        Self {
            w_self,
            w_msg,
            bias,
            activation,
        }
    }
    /// Create a randomly-initialised MPNN layer (deterministic seed).
    pub fn random(n_hidden: usize, activation: Activation, seed: u64) -> Self {
        let mut lcg = Lcg::new(seed);
        let scale = (2.0 / (2 * n_hidden) as f64).sqrt();
        let make_mat = |lcg: &mut Lcg| -> Vec<Vec<f64>> {
            (0..n_hidden)
                .map(|_| (0..n_hidden).map(|_| lcg.next_f64(scale)).collect())
                .collect()
        };
        let w_self = make_mat(&mut lcg);
        let w_msg = make_mat(&mut lcg);
        let bias = vec![0.0; n_hidden];
        Self {
            w_self,
            w_msg,
            bias,
            activation,
        }
    }
    /// Forward pass for a single node.
    ///
    /// `h_self` – hidden state of node i.
    /// `neighbour_sum` – pre-aggregated sum of neighbour hidden states.
    pub fn forward_node(&self, h_self: &[f64], neighbour_sum: &[f64]) -> Vec<f64> {
        let n = self.bias.len();
        let mut out = self.bias.clone();
        for (i, out_i) in out.iter_mut().enumerate().take(n) {
            for j in 0..h_self.len() {
                *out_i += self.w_self[i][j] * h_self[j];
                *out_i += self.w_msg[i][j] * neighbour_sum[j];
            }
        }
        self.activation.apply(&mut out);
        out
    }
}
/// Batch evaluator for Behler-Parrinello symmetry functions.
///
/// This struct wraps a [`SymmetryFunctionSet`] and a cutoff radius and
/// provides a `compute_symmetry_functions_batch` method that evaluates
/// descriptors for every atom in a structure in a single call, returning
/// all atom-centered descriptors as a 2-D matrix (one row per atom).
///
/// # Usage
/// ```text
/// let bp = BehlerParrinello::new(sym_set, 6.0);
/// let batch = bp.compute_symmetry_functions_batch(&positions, &species);
/// // batch[i] is the descriptor vector for atom i
/// ```
#[derive(Debug, Clone)]
pub struct BehlerParrinello {
    /// Symmetry function set defining the descriptor features.
    pub symmetry_functions: SymmetryFunctionSet,
    /// Cutoff radius (Å).
    pub cutoff: f64,
}
impl BehlerParrinello {
    /// Create a new batch evaluator.
    ///
    /// # Arguments
    /// * `symmetry_functions` – the set of symmetry functions to evaluate.
    /// * `cutoff`             – cutoff radius for all symmetry functions.
    pub fn new(symmetry_functions: SymmetryFunctionSet, cutoff: f64) -> Self {
        Self {
            symmetry_functions,
            cutoff,
        }
    }
    /// Create a BehlerParrinello evaluator with a standard set of G2 radial
    /// functions (convenience constructor for prototyping).
    ///
    /// The set contains G2 functions with `eta` values `\[0.1, 0.5, 1.0, 2.0\]`
    /// and `rs = 0`.
    pub fn default_radial(cutoff: f64) -> Self {
        let mut set = SymmetryFunctionSet::new();
        for &eta in &[0.1_f64, 0.5, 1.0, 2.0] {
            set.push(SymmetryFunction::G2 {
                eta,
                rs: 0.0,
                rc: cutoff,
            });
        }
        Self::new(set, cutoff)
    }
    /// Evaluate symmetry functions for all atoms in a structure.
    ///
    /// Computes the atom-centered descriptor vector for every atom
    /// independently and returns the result as a `Vec<Vec`f64`>` where
    /// `result\[i\]` is the descriptor for atom `i`.
    ///
    /// # Arguments
    /// * `positions` – Cartesian positions `\[f64; 3\]` for each atom.
    /// * `species`   – integer species labels (same length as `positions`).
    ///
    /// # Returns
    /// Matrix of descriptors: `result\[i\]\[k\]` is the k-th descriptor
    /// value for atom `i`.  Shape: `N_atoms × N_desc`.
    pub fn compute_symmetry_functions_batch(
        &self,
        positions: &[[f64; 3]],
        species: &[u32],
    ) -> Vec<Vec<f64>> {
        let n = positions.len();
        let pos_vec3: Vec<Vec3> = positions
            .iter()
            .map(|p| Vec3::new(p[0], p[1], p[2]))
            .collect();
        (0..n)
            .map(|i| {
                self.symmetry_functions
                    .compute_descriptor(&pos_vec3, i, species, self.cutoff)
            })
            .collect()
    }
    /// Compute the descriptor for a single atom.
    ///
    /// Convenience wrapper around [`SymmetryFunctionSet::compute_descriptor`].
    ///
    /// # Arguments
    /// * `positions` – all atom positions as `\[f64; 3\]`.
    /// * `center`    – index of the atom whose descriptor is computed.
    /// * `species`   – species labels.
    pub fn compute_single(
        &self,
        positions: &[[f64; 3]],
        center: usize,
        species: &[u32],
    ) -> Vec<f64> {
        let pos_vec3: Vec<Vec3> = positions
            .iter()
            .map(|p| Vec3::new(p[0], p[1], p[2]))
            .collect();
        self.symmetry_functions
            .compute_descriptor(&pos_vec3, center, species, self.cutoff)
    }
    /// Descriptor dimensionality (number of symmetry functions).
    pub fn n_descriptors(&self) -> usize {
        self.symmetry_functions.functions.len()
    }
}
/// Angular triplet (i, j, k) for angular message passing.
#[derive(Debug, Clone, Copy)]
pub struct AngularTriplet {
    /// Central atom index.
    pub i: usize,
    /// First neighbour.
    pub j: usize,
    /// Second neighbour.
    pub k: usize,
}
/// Minimal LCG for deterministic weight initialisation.
pub(super) struct Lcg {
    pub(super) state: u64,
}
impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    /// Return a pseudo-random f64 in `\[-scale, scale\]`.
    fn next_f64(&mut self, scale: f64) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let u = ((self.state >> 33) as f64) / (2.0_f64.powi(31));
        (u - 0.5) * 2.0 * scale
    }
}
/// Ensemble of neural network potentials exposing a `predict_energy` /
/// `predict_forces` API that returns both the mean prediction and an
/// uncertainty estimate (standard deviation across models).
#[derive(Debug, Clone)]
pub struct EnsembleNnp {
    /// The individual NNP models in the ensemble.
    pub models: Vec<FeedForwardPotential>,
}
impl EnsembleNnp {
    /// Create an ensemble from a list of models.
    ///
    /// At least one model must be provided.
    pub fn new(models: Vec<FeedForwardPotential>) -> Self {
        assert!(
            !models.is_empty(),
            "EnsembleNnp must have at least one model"
        );
        Self { models }
    }
    /// Predict energy for a given descriptor.
    ///
    /// Returns `(mean, std)` across all models in the ensemble.
    /// Each model evaluates the descriptor through its network layers.
    pub fn predict_energy(&self, descriptor: &[f64]) -> (f64, f64) {
        let energies: Vec<f64> = self
            .models
            .iter()
            .map(|m| m.evaluate_network(descriptor))
            .collect();
        let n = energies.len() as f64;
        let mean = energies.iter().sum::<f64>() / n;
        let var = energies
            .iter()
            .map(|&e| (e - mean) * (e - mean))
            .sum::<f64>()
            / n;
        (mean, var.sqrt())
    }
    /// Predict forces for a set of atom positions.
    ///
    /// Returns the mean force vector per atom across all ensemble members.
    /// Forces are obtained by finite differences on the mean energy.
    pub fn predict_forces(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let species: Vec<u32> = vec![1; n];
        let pos_vec3: Vec<oxiphysics_core::math::Vec3> = positions
            .iter()
            .map(|p| oxiphysics_core::math::Vec3::new(p[0], p[1], p[2]))
            .collect();
        let n_models = self.models.len() as f64;
        let mut force_sum = vec![[0.0_f64; 3]; n];
        for model in &self.models {
            let (_, forces) = model.energy_and_forces(&pos_vec3, &species);
            for (i, f) in forces.iter().enumerate() {
                force_sum[i][0] += f[0];
                force_sum[i][1] += f[1];
                force_sum[i][2] += f[2];
            }
        }
        force_sum.iter_mut().for_each(|f| {
            f[0] /= n_models;
            f[1] /= n_models;
            f[2] /= n_models;
        });
        force_sum
    }
    /// Number of models in the ensemble.
    pub fn n_models(&self) -> usize {
        self.models.len()
    }
    /// Compute the per-descriptor disagreement (standard deviation) across
    /// ensemble members for a single atom descriptor vector.
    ///
    /// For each feature index `k` of the descriptor, the disagreement is
    /// defined as the standard deviation of the descriptor values predicted
    /// by evaluating the *network output sensitivity* (first-layer response)
    /// across models.  As a practical approximation this computes the
    /// standard deviation of the network outputs for a set of perturbed
    /// descriptors `G + δ_k · e_k`:
    ///
    /// ```text
    /// disagreement_k = std_dev_m { NN_m(G + δ_k * e_k) }
    /// ```
    ///
    /// This gives a per-descriptor uncertainty measure that reflects how
    /// differently the ensemble members respond to small changes in each
    /// feature.
    ///
    /// # Arguments
    /// * `descriptor` – atom-centered descriptor vector G.
    /// * `delta`      – finite-difference step δ_k (typically 1e-3).
    ///
    /// # Returns
    /// A vector of length `descriptor.len()` with the ensemble disagreement
    /// for each descriptor component.
    pub fn compute_disagreement(&self, descriptor: &[f64], delta: f64) -> Vec<f64> {
        let dim = descriptor.len();
        let mut result = vec![0.0_f64; dim];
        let n_models = self.models.len() as f64;
        for k in 0..dim {
            let mut desc_pert = descriptor.to_vec();
            desc_pert[k] += delta;
            let energies: Vec<f64> = self
                .models
                .iter()
                .map(|m| m.evaluate_network(&desc_pert))
                .collect();
            let mean = energies.iter().sum::<f64>() / n_models;
            let var = energies
                .iter()
                .map(|&e| (e - mean) * (e - mean))
                .sum::<f64>()
                / n_models;
            result[k] = var.sqrt();
        }
        result
    }
}
/// A simple message-passing neural network for molecular property prediction.
///
/// Architecture: initial embedding → T message-passing layers → readout sum.
#[derive(Debug, Clone)]
pub struct Mpnn {
    /// Dimensionality of hidden node embeddings.
    pub n_hidden: usize,
    /// Message-passing layers.
    pub mp_layers: Vec<MpnnLayer>,
    /// Readout layer: maps final hidden state to scalar energy per atom.
    pub readout: DenseLayer,
    /// Cutoff distance for neighbour graph construction.
    pub cutoff: f64,
}
impl Mpnn {
    /// Build a random MPNN with `n_layers` message-passing steps.
    pub fn new(n_hidden: usize, n_layers: usize, cutoff: f64) -> Self {
        let mut mp_layers = Vec::with_capacity(n_layers);
        for l in 0..n_layers {
            mp_layers.push(MpnnLayer::random(
                n_hidden,
                Activation::Tanh,
                100 + l as u64,
            ));
        }
        let scale = (2.0 / (n_hidden + 1) as f64).sqrt();
        let mut lcg = Lcg::new(999);
        let w_read = vec![(0..n_hidden).map(|_| lcg.next_f64(scale)).collect()];
        let b_read = vec![0.0];
        let readout = DenseLayer::new(w_read, b_read, Activation::Linear);
        Self {
            n_hidden,
            mp_layers,
            readout,
            cutoff,
        }
    }
    /// Compute per-atom energies using message passing.
    ///
    /// Initial node embeddings are set to a constant vector of ones (proof-of-concept).
    pub fn atom_energies(&self, positions: &[[f64; 3]]) -> Vec<f64> {
        let n = positions.len();
        let mut h: Vec<Vec<f64>> = (0..n).map(|_| vec![0.1; self.n_hidden]).collect();
        for layer in &self.mp_layers {
            let mut h_new = h.clone();
            for i in 0..n {
                let mut neigh_sum = vec![0.0; self.n_hidden];
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let dx = positions[j][0] - positions[i][0];
                    let dy = positions[j][1] - positions[i][1];
                    let dz = positions[j][2] - positions[i][2];
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    if r < self.cutoff {
                        for k in 0..self.n_hidden {
                            neigh_sum[k] += h[j][k];
                        }
                    }
                }
                h_new[i] = layer.forward_node(&h[i], &neigh_sum);
            }
            h = h_new;
        }
        h.iter().map(|hi| self.readout.forward(hi)[0]).collect()
    }
    /// Total energy (sum of atomic contributions).
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        self.atom_energies(positions).iter().sum()
    }
}
/// SchNet interaction block: a continuous-filter convolution layer.
///
/// For each atom i, collects contributions from all neighbours j:
///   h_i^(l+1) = h_i^l + sum_j  h_j^l ⊙ W(r_ij)
/// where W(r) = DenseNetwork(basis(r)).
#[derive(Clone)]
pub struct SchNetInteraction {
    /// Filter network: maps basis vector → feature vector.
    pub filter_network: Vec<DenseLayer>,
    /// Number of features.
    pub n_features: usize,
    /// Gaussian basis for distance expansion.
    pub basis: GaussianBasis,
}
impl SchNetInteraction {
    /// Build a SchNet interaction layer.
    pub fn new(n_features: usize, basis: GaussianBasis) -> Self {
        let n_basis = basis.n_basis();
        let mut lcg = Lcg::new(555);
        let scale = (2.0 / (n_basis + n_features) as f64).sqrt();
        let w1: Vec<Vec<f64>> = (0..n_features)
            .map(|_| (0..n_basis).map(|_| lcg.next_f64(scale)).collect())
            .collect();
        let layer1 = DenseLayer::new(w1, vec![0.0; n_features], Activation::Tanh);
        let scale2 = (2.0 / (2 * n_features) as f64).sqrt();
        let w2: Vec<Vec<f64>> = (0..n_features)
            .map(|_| (0..n_features).map(|_| lcg.next_f64(scale2)).collect())
            .collect();
        let layer2 = DenseLayer::new(w2, vec![0.0; n_features], Activation::Linear);
        Self {
            filter_network: vec![layer1, layer2],
            n_features,
            basis,
        }
    }
    fn run_filter(&self, r: f64) -> Vec<f64> {
        let mut x = self.basis.evaluate(r);
        for layer in &self.filter_network {
            x = layer.forward(&x);
        }
        x
    }
    /// Apply one interaction step to all atom features.
    pub fn forward(&self, features: &[Vec<f64>], positions: &[[f64; 3]]) -> Vec<Vec<f64>> {
        let n = features.len();
        let cutoff = self.basis.cutoff;
        let mut new_features = features.to_vec();
        for i in 0..n {
            let mut delta = vec![0.0_f64; self.n_features];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r >= cutoff {
                    continue;
                }
                let w = self.run_filter(r);
                for k in 0..self.n_features {
                    delta[k] += features[j][k] * w[k];
                }
            }
            for k in 0..self.n_features {
                new_features[i][k] += delta[k];
            }
        }
        new_features
    }
}
/// Activation function for a dense layer.
#[derive(Debug, Clone, Copy)]
pub enum Activation {
    /// Hyperbolic tangent.
    Tanh,
    /// Logistic sigmoid 1/(1+exp(-x)).
    Sigmoid,
    /// Rectified linear unit max(0, x).
    ReLU,
    /// Identity (no activation).
    Linear,
}
impl Activation {
    /// Apply the activation element-wise to a mutable slice.
    fn apply(&self, v: &mut [f64]) {
        for x in v.iter_mut() {
            *x = match self {
                Activation::Tanh => x.tanh(),
                Activation::Sigmoid => 1.0 / (1.0 + (-*x).exp()),
                Activation::ReLU => x.max(0.0),
                Activation::Linear => *x,
            };
        }
    }
}
/// Graph Attention Network Potential.
#[derive(Clone)]
pub struct GatPotential {
    /// GAT interaction layers.
    pub gat_layers: Vec<GraphAttentionLayer>,
    /// Readout dense layer.
    pub readout: DenseLayer,
    /// Cutoff (Å).
    pub cutoff: f64,
    /// Number of node features.
    pub n_features: usize,
}
impl GatPotential {
    /// Create a GAT potential.
    pub fn new(n_features: usize, n_layers: usize, cutoff: f64) -> Self {
        let gat_layers = (0..n_layers)
            .map(|_| GraphAttentionLayer::new(n_features, n_features, cutoff))
            .collect();
        let mut lcg = Lcg::new(1234);
        let scale = (2.0 / (n_features + 1) as f64).sqrt();
        let w = vec![
            (0..n_features)
                .map(|_| lcg.next_f64(scale))
                .collect::<Vec<f64>>(),
        ];
        let readout = DenseLayer::new(w, vec![0.0], Activation::Linear);
        Self {
            gat_layers,
            readout,
            cutoff,
            n_features,
        }
    }
    /// Total energy.
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let n = positions.len();
        let mut features: Vec<Vec<f64>> = (0..n).map(|_| vec![0.1; self.n_features]).collect();
        for layer in &self.gat_layers {
            features = layer.forward(&features, positions);
        }
        let mut energy = 0.0;
        for f in &features {
            energy += self.readout.forward(f).iter().sum::<f64>();
        }
        energy
    }
}
/// Validator that checks energy conservation along a trajectory.
pub struct EnergyConservationValidator {
    /// Tolerance for relative energy drift.
    pub tolerance: f64,
    /// Energy values recorded at each step.
    pub energy_history: Vec<f64>,
}
impl EnergyConservationValidator {
    /// Create a new validator with a given relative tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            energy_history: Vec::new(),
        }
    }
    /// Record the energy at the current step.
    pub fn record(&mut self, energy: f64) {
        self.energy_history.push(energy);
    }
    /// Check if energy drift is within tolerance.
    ///
    /// Returns `Ok(max_drift)` if conserved, `Err(max_drift)` otherwise.
    pub fn validate(&self) -> Result<f64, f64> {
        if self.energy_history.len() < 2 {
            return Ok(0.0);
        }
        let e0 = self.energy_history[0];
        let max_drift = self
            .energy_history
            .iter()
            .map(|&e| {
                if e0.abs() > 1e-30 {
                    (e - e0).abs() / e0.abs()
                } else {
                    (e - e0).abs()
                }
            })
            .fold(0.0_f64, f64::max);
        if max_drift <= self.tolerance {
            Ok(max_drift)
        } else {
            Err(max_drift)
        }
    }
    /// Maximum observed drift.
    pub fn max_drift(&self) -> f64 {
        match self.validate() {
            Ok(d) | Err(d) => d,
        }
    }
    /// Clear recorded history.
    pub fn reset(&mut self) {
        self.energy_history.clear();
    }
    /// Number of steps recorded.
    pub fn n_steps(&self) -> usize {
        self.energy_history.len()
    }
}
/// Descriptor type for the DeePMD-style environment matrix.
#[derive(Debug, Clone)]
pub struct DeepPotDescriptor {
    /// Maximum number of neighbours per atom.
    pub max_neighbours: usize,
    /// Cutoff radius (Å).
    pub cutoff: f64,
    /// Inner cutoff for smooth transition (Å).
    pub inner_cutoff: f64,
}
impl DeepPotDescriptor {
    /// Create a new DeePMD descriptor builder.
    pub fn new(max_neighbours: usize, cutoff: f64, inner_cutoff: f64) -> Self {
        assert!(inner_cutoff < cutoff, "inner_cutoff must be < cutoff");
        Self {
            max_neighbours,
            cutoff,
            inner_cutoff,
        }
    }
    /// Smooth switching function s(r) used in DeePMD:
    ///   s(r) = 1/r  if r < r_inner
    ///   s(r) = (1/r) * smooth(r, r_inner, r_cut)  otherwise
    pub fn smooth_fn(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        if r <= 0.0 {
            return 0.0;
        }
        let inv_r = 1.0 / r;
        if r <= self.inner_cutoff {
            return inv_r;
        }
        let u = (r - self.inner_cutoff) / (self.cutoff - self.inner_cutoff);
        let envelope = 1.0 - 3.0 * u * u + 2.0 * u * u * u;
        inv_r * envelope
    }
    /// Compute the DeePMD environment matrix for atom `center`.
    ///
    /// Returns a flattened row-major matrix of shape [M × 4]:
    ///   column 0: s(r_ij)
    ///   columns 1-3: s(r_ij) * r_hat_ij[0..2]
    pub fn environment_matrix(&self, positions: &[[f64; 3]], center: usize) -> Vec<f64> {
        let _n = positions.len();
        let mut rows: Vec<[f64; 4]> = Vec::new();
        let rc = positions[center];
        for (j, pos_j) in positions.iter().enumerate() {
            if j == center {
                continue;
            }
            let dx = pos_j[0] - rc[0];
            let dy = pos_j[1] - rc[1];
            let dz = pos_j[2] - rc[2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            if r >= self.cutoff {
                continue;
            }
            let s = self.smooth_fn(r);
            let inv_r = if r > 1e-15 { 1.0 / r } else { 0.0 };
            rows.push([s, s * dx * inv_r, s * dy * inv_r, s * dz * inv_r]);
        }
        rows.truncate(self.max_neighbours);
        while rows.len() < self.max_neighbours {
            rows.push([0.0; 4]);
        }
        rows.iter().flat_map(|row| row.iter().copied()).collect()
    }
    /// Descriptor dimension = max_neighbours × 4.
    pub fn descriptor_dim(&self) -> usize {
        self.max_neighbours * 4
    }
}
/// Neural network pair potential: E_pair(r) = NN(r).
///
/// A simple 1D neural network that maps an interatomic distance to
/// a pair energy contribution.
#[derive(Clone)]
pub struct NnPairPotential {
    /// The underlying dense layers.
    pub layers: Vec<DenseLayer>,
    /// Cutoff radius for the pair interaction.
    pub cutoff: f64,
}
impl NnPairPotential {
    /// Create a pair potential with a single hidden layer.
    ///
    /// Architecture: 1 → n_hidden → 1.
    pub fn new(n_hidden: usize, cutoff: f64) -> Self {
        let mut rng_state: u64 = 12345;
        let mut next_w = |scale: f64| -> f64 {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u = ((rng_state >> 33) as f64) / (2.0_f64.powi(31));
            (u - 0.5) * 2.0 * scale
        };
        let scale1 = (2.0 / (1 + n_hidden) as f64).sqrt();
        let w1: Vec<Vec<f64>> = (0..n_hidden).map(|_| vec![next_w(scale1)]).collect();
        let b1 = vec![0.0; n_hidden];
        let scale2 = (2.0 / (n_hidden + 1) as f64).sqrt();
        let w2 = vec![(0..n_hidden).map(|_| next_w(scale2)).collect::<Vec<f64>>()];
        let b2 = vec![0.0];
        let layers = vec![
            DenseLayer::new(w1, b1, Activation::Tanh),
            DenseLayer::new(w2, b2, Activation::Linear),
        ];
        Self { layers, cutoff }
    }
    /// Evaluate the NN for a single distance r.
    pub fn evaluate(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        let mut x = vec![r];
        for layer in &self.layers {
            x = layer.forward(&x);
        }
        x[0]
    }
    /// Compute pair energy and force for atoms i and j at distance r.
    ///
    /// Returns (energy, force_magnitude) where force is the scalar dE/dr.
    pub fn energy_and_force(&self, r: f64) -> (f64, f64) {
        let e = self.evaluate(r);
        let h = 1e-5;
        let de_dr = if r >= self.cutoff {
            0.0
        } else {
            (self.evaluate(r + h) - self.evaluate(r - h)) / (2.0 * h)
        };
        (e, -de_dr)
    }
    /// Total energy for a set of atom positions (sum over all pairs within cutoff).
    pub fn total_energy(&self, positions: &[Vec3]) -> f64 {
        let n = positions.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = (positions[i] - positions[j]).norm();
                energy += self.evaluate(r);
            }
        }
        energy
    }
}
/// Statistics of a descriptor over a training set.
#[derive(Debug, Clone)]
pub struct DescriptorStats {
    /// Per-feature minimum.
    pub min: Vec<f64>,
    /// Per-feature maximum.
    pub max: Vec<f64>,
    /// Per-feature mean.
    pub mean: Vec<f64>,
    /// Per-feature standard deviation.
    pub std: Vec<f64>,
}
impl DescriptorStats {
    /// Compute statistics from a collection of descriptor vectors.
    pub fn compute(descriptors: &[Vec<f64>]) -> Self {
        assert!(!descriptors.is_empty());
        let dim = descriptors[0].len();
        let n = descriptors.len() as f64;
        let mut min = vec![f64::INFINITY; dim];
        let mut max = vec![f64::NEG_INFINITY; dim];
        let mut mean = vec![0.0f64; dim];
        for d in descriptors {
            for (i, &v) in d.iter().enumerate() {
                if v < min[i] {
                    min[i] = v;
                }
                if v > max[i] {
                    max[i] = v;
                }
                mean[i] += v;
            }
        }
        for m in mean.iter_mut() {
            *m /= n;
        }
        let mut var = vec![0.0f64; dim];
        for d in descriptors {
            for (i, &v) in d.iter().enumerate() {
                var[i] += (v - mean[i]) * (v - mean[i]);
            }
        }
        let std: Vec<f64> = var.iter().map(|&v| (v / n).sqrt()).collect();
        Self {
            min,
            max,
            mean,
            std,
        }
    }
    /// Range (max - min) per feature.
    pub fn range(&self) -> Vec<f64> {
        self.min
            .iter()
            .zip(self.max.iter())
            .map(|(lo, hi)| hi - lo)
            .collect()
    }
}
/// A single dense (fully connected) layer: y = activation(W * x + b).
#[derive(Debug, Clone)]
pub struct DenseLayer {
    /// Weight matrix stored row-major: `weights\[i\]\[j\]` is the weight from
    /// input `j` to output `i`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector (length = number of outputs).
    pub biases: Vec<f64>,
    /// Activation function applied after the affine transform.
    pub activation: Activation,
}
impl DenseLayer {
    /// Create a new dense layer with the given weights, biases, and
    /// activation.
    pub fn new(weights: Vec<Vec<f64>>, biases: Vec<f64>, activation: Activation) -> Self {
        Self {
            weights,
            biases,
            activation,
        }
    }
    /// Forward pass: compute `activation(W * input + b)`.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let n_out = self.biases.len();
        let mut output = self.biases.clone();
        for (i, out_i) in output.iter_mut().enumerate().take(n_out) {
            for (j, &inp) in input.iter().enumerate() {
                *out_i += self.weights[i][j] * inp;
            }
        }
        self.activation.apply(&mut output);
        output
    }
}
/// Gaussian basis function expansion for interatomic distances.
///
/// Maps r to a vector of Gaussian basis values:
///   phi_k(r) = exp(-0.5 * ((r - mu_k) / sigma)^2)
#[derive(Debug, Clone)]
pub struct GaussianBasis {
    /// Centers of the Gaussians.
    pub centers: Vec<f64>,
    /// Width (sigma) of each Gaussian.
    pub sigma: f64,
    /// Cutoff radius.
    pub cutoff: f64,
}
impl GaussianBasis {
    /// Create an evenly-spaced Gaussian basis from `r_min` to `r_max`.
    pub fn new(r_min: f64, r_max: f64, n_basis: usize, sigma: f64, cutoff: f64) -> Self {
        let centers = if n_basis <= 1 {
            vec![r_min]
        } else {
            let step = (r_max - r_min) / (n_basis - 1) as f64;
            (0..n_basis).map(|k| r_min + k as f64 * step).collect()
        };
        Self {
            centers,
            sigma,
            cutoff,
        }
    }
    /// Evaluate all Gaussian basis values at distance `r`.
    pub fn evaluate(&self, r: f64) -> Vec<f64> {
        if r >= self.cutoff {
            return vec![0.0; self.centers.len()];
        }
        self.centers
            .iter()
            .map(|&mu| {
                let x = (r - mu) / self.sigma;
                (-0.5 * x * x).exp()
            })
            .collect()
    }
    /// Number of basis functions.
    pub fn n_basis(&self) -> usize {
        self.centers.len()
    }
}
/// SchNet-style neural network potential.
#[derive(Clone)]
pub struct SchNetPotential {
    /// Number of atom features.
    pub n_features: usize,
    /// Interaction layers.
    pub interactions: Vec<SchNetInteraction>,
    /// Readout network: features → atomic energy.
    pub readout: Vec<DenseLayer>,
    /// Cutoff (Å).
    pub cutoff: f64,
}
impl SchNetPotential {
    /// Create a SchNet potential.
    pub fn new(n_features: usize, n_interactions: usize, n_basis: usize, cutoff: f64) -> Self {
        let basis = GaussianBasis::new(0.5, cutoff, n_basis, 0.5, cutoff);
        let interactions = (0..n_interactions)
            .map(|_| SchNetInteraction::new(n_features, basis.clone()))
            .collect();
        let mut lcg = Lcg::new(888);
        let scale = (2.0 / (n_features + 1) as f64).sqrt();
        let w: Vec<Vec<f64>> = vec![(0..n_features).map(|_| lcg.next_f64(scale)).collect()];
        let readout = vec![DenseLayer::new(w, vec![0.0], Activation::Linear)];
        Self {
            n_features,
            interactions,
            readout,
            cutoff,
        }
    }
    /// Total energy for a set of positions.
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let n = positions.len();
        let mut features: Vec<Vec<f64>> = (0..n).map(|_| vec![0.1; self.n_features]).collect();
        for interaction in &self.interactions {
            features = interaction.forward(&features, positions);
        }
        let mut energy = 0.0;
        for f in &features {
            let mut x = f.clone();
            for layer in &self.readout {
                x = layer.forward(&x);
            }
            energy += x.iter().sum::<f64>();
        }
        energy
    }
}
/// A single snapshot of an atomistic system using a neural-network potential.
#[derive(Debug, Clone)]
pub struct NnAtomisticSystem {
    /// Atom positions (Cartesian, Å).
    pub positions: Vec<oxiphysics_core::math::Vec3>,
    /// Atom velocities (Å/ps).
    pub velocities: Vec<oxiphysics_core::math::Vec3>,
    /// Atom masses (amu).
    pub masses: Vec<f64>,
    /// Species labels.
    pub species: Vec<u32>,
}
impl NnAtomisticSystem {
    /// Create a new system at rest.
    pub fn new(
        positions: Vec<oxiphysics_core::math::Vec3>,
        masses: Vec<f64>,
        species: Vec<u32>,
    ) -> Self {
        let n = positions.len();
        Self {
            positions,
            velocities: vec![oxiphysics_core::math::Vec3::zeros(); n],
            masses,
            species,
        }
    }
    /// Perform one velocity-Verlet step.
    ///
    /// `potential` provides forces via `energy_and_forces`.
    /// `dt` is the time step in ps.
    pub fn step(&mut self, potential: &FeedForwardPotential, dt: f64) {
        let (_, forces) = potential.energy_and_forces(&self.positions, &self.species);
        let n = self.positions.len();
        let mut half_vel = self.velocities.clone();
        for i in 0..n {
            let a = forces[i] * (1.0 / self.masses[i]);
            half_vel[i] = self.velocities[i] + a * (0.5 * dt);
            self.positions[i] += half_vel[i] * dt;
        }
        let (_, forces_new) = potential.energy_and_forces(&self.positions, &self.species);
        for i in 0..n {
            let a_new = forces_new[i] * (1.0 / self.masses[i]);
            self.velocities[i] = half_vel[i] + a_new * (0.5 * dt);
        }
    }
    /// Compute instantaneous kinetic energy (amu·Å²/ps²).
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * v.norm_squared())
            .sum()
    }
    /// Compute total energy = kinetic + potential.
    pub fn total_energy(&self, potential: &FeedForwardPotential) -> f64 {
        let (e_pot, _) = potential.energy_and_forces(&self.positions, &self.species);
        e_pot + self.kinetic_energy()
    }
}
/// Normalises symmetry-function descriptors to zero mean and unit variance.
///
/// Equivalent to [`DescriptorNormalizer`] but exposes a `fit` / `transform`
/// naming convention that matches common ML-library conventions.
#[derive(Debug, Clone)]
pub struct NnpNormalizer {
    /// Per-feature mean (length = descriptor dimension).
    pub mean: Vec<f64>,
    /// Per-feature standard deviation (clamped ≥ 1e-12).
    pub std: Vec<f64>,
}
impl NnpNormalizer {
    /// Fit the normalizer on a collection of descriptor vectors.
    ///
    /// Computes per-feature mean and standard deviation; the floor ensures
    /// no division-by-zero when a feature is constant.
    pub fn fit(data: &[Vec<f64>]) -> Self {
        assert!(!data.is_empty(), "NnpNormalizer::fit called on empty data");
        let dim = data[0].len();
        let n = data.len() as f64;
        let mut mean = vec![0.0_f64; dim];
        for row in data {
            assert_eq!(
                row.len(),
                dim,
                "descriptor dimension mismatch in NnpNormalizer::fit"
            );
            for (j, &v) in row.iter().enumerate() {
                mean[j] += v;
            }
        }
        for m in mean.iter_mut() {
            *m /= n;
        }
        let mut var = vec![0.0_f64; dim];
        for row in data {
            for (j, &v) in row.iter().enumerate() {
                let d = v - mean[j];
                var[j] += d * d;
            }
        }
        let std: Vec<f64> = var.iter().map(|&v| (v / n).sqrt().max(1e-12)).collect();
        Self { mean, std }
    }
    /// Transform a descriptor vector to zero mean / unit variance.
    ///
    /// Returns a new vector; the input is not modified.
    pub fn transform(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(
            x.len(),
            self.mean.len(),
            "descriptor dimension mismatch in NnpNormalizer::transform"
        );
        x.iter()
            .enumerate()
            .map(|(i, &v)| (v - self.mean[i]) / self.std[i])
            .collect()
    }
    /// Inverse transform (denormalize) back to original scale.
    pub fn inverse_transform(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(
            x.len(),
            self.mean.len(),
            "dimension mismatch in NnpNormalizer::inverse_transform"
        );
        x.iter()
            .enumerate()
            .map(|(i, &v)| v * self.std[i] + self.mean[i])
            .collect()
    }
}
/// DeePMD-style neural network potential.
///
/// Uses the environment matrix as descriptor input to a feedforward network.
#[derive(Clone)]
pub struct DeepPotPotential {
    /// Descriptor builder.
    pub descriptor: DeepPotDescriptor,
    /// Embedding network layers.
    pub embedding: Vec<DenseLayer>,
    /// Fitting network layers.
    pub fitting: Vec<DenseLayer>,
}
impl DeepPotPotential {
    /// Create a DeePMD potential with random weights.
    ///
    /// # Arguments
    /// * `max_neigh`    – Max neighbours in environment matrix.
    /// * `embed_sizes`  – Hidden sizes for embedding network (e.g., [8, 16]).
    /// * `fit_sizes`    – Hidden sizes for fitting network (e.g., [64, 64, 1]).
    /// * `cutoff`       – Outer cutoff (Å).
    /// * `inner_cutoff` – Inner cutoff (Å).
    pub fn new(
        max_neigh: usize,
        embed_sizes: &[usize],
        fit_sizes: &[usize],
        cutoff: f64,
        inner_cutoff: f64,
    ) -> Self {
        let descriptor = DeepPotDescriptor::new(max_neigh, cutoff, inner_cutoff);
        let input_dim = descriptor.descriptor_dim();
        let mut lcg = Lcg::new(777);
        let make_layers = |dims: &[usize], lcg: &mut Lcg| -> Vec<DenseLayer> {
            let mut layers = Vec::new();
            for w in dims.windows(2) {
                let n_in = w[0];
                let n_out = w[1];
                let scale = (2.0 / (n_in + n_out) as f64).sqrt();
                let weights: Vec<Vec<f64>> = (0..n_out)
                    .map(|_| (0..n_in).map(|_| lcg.next_f64(scale)).collect())
                    .collect();
                let biases = vec![0.0; n_out];
                layers.push(DenseLayer::new(weights, biases, Activation::Tanh));
            }
            layers
        };
        let mut embed_full = vec![input_dim];
        embed_full.extend_from_slice(embed_sizes);
        let embedding = make_layers(&embed_full, &mut lcg);
        let embed_out = *embed_sizes.last().unwrap_or(&input_dim);
        let mut fit_full = vec![embed_out];
        fit_full.extend_from_slice(fit_sizes);
        let fitting = make_layers(&fit_full, &mut lcg);
        Self {
            descriptor,
            embedding,
            fitting,
        }
    }
    fn run_layers(layers: &[DenseLayer], input: &[f64]) -> Vec<f64> {
        let mut x = input.to_vec();
        for layer in layers {
            x = layer.forward(&x);
        }
        x
    }
    /// Compute atomic energy for a single atom.
    pub fn atomic_energy(&self, positions: &[[f64; 3]], center: usize) -> f64 {
        let desc = self.descriptor.environment_matrix(positions, center);
        let embedded = Self::run_layers(&self.embedding, &desc);
        let out = Self::run_layers(&self.fitting, &embedded);
        out.iter().sum()
    }
    /// Total energy = sum of atomic energies.
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        (0..positions.len())
            .map(|i| self.atomic_energy(positions, i))
            .sum()
    }
    /// Forces via finite differences.
    pub fn forces(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let h = 1e-5;
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut pos = positions.to_vec();
        for i in 0..n {
            for a in 0..3 {
                pos[i][a] += h;
                let ep = self.total_energy(&pos);
                pos[i][a] -= 2.0 * h;
                let em = self.total_energy(&pos);
                pos[i][a] += h;
                forces[i][a] = -(ep - em) / (2.0 * h);
            }
        }
        forces
    }
}
/// Normalises descriptor vectors to zero mean and unit variance.
///
/// Given a training set of descriptor vectors, computes per-feature mean and
/// standard deviation, then normalises new descriptors via:
///
/// ```text
/// x_norm = (x - mean) / std
/// ```
#[derive(Debug, Clone)]
pub struct DescriptorNormalizer {
    /// Per-feature mean.
    pub mean: Vec<f64>,
    /// Per-feature standard deviation (clamped to a floor to avoid division by zero).
    pub std: Vec<f64>,
}
impl DescriptorNormalizer {
    /// Create a normalizer with explicit mean and std.
    pub fn new(mean: Vec<f64>, std: Vec<f64>) -> Self {
        assert_eq!(mean.len(), std.len(), "mean and std length mismatch");
        Self { mean, std }
    }
    /// Fit the normalizer to a batch of descriptor vectors.
    ///
    /// `data` is a slice of descriptor vectors, each of length `dim`.
    pub fn fit(data: &[Vec<f64>]) -> Self {
        assert!(!data.is_empty(), "cannot fit on empty data");
        let dim = data[0].len();
        let n = data.len() as f64;
        let mut mean = vec![0.0; dim];
        for row in data {
            for (j, &val) in row.iter().enumerate() {
                mean[j] += val;
            }
        }
        for v in mean.iter_mut() {
            *v /= n;
        }
        let mut variance = vec![0.0; dim];
        for row in data {
            for (j, &val) in row.iter().enumerate() {
                let diff = val - mean[j];
                variance[j] += diff * diff;
            }
        }
        let std: Vec<f64> = variance
            .iter()
            .map(|&v| (v / n).sqrt().max(1e-12))
            .collect();
        Self { mean, std }
    }
    /// Normalize a single descriptor vector in-place.
    pub fn normalize(&self, descriptor: &mut [f64]) {
        assert_eq!(descriptor.len(), self.mean.len(), "dimension mismatch");
        for (i, x) in descriptor.iter_mut().enumerate() {
            *x = (*x - self.mean[i]) / self.std[i];
        }
    }
    /// Normalize a single descriptor, returning a new vector.
    pub fn normalize_vec(&self, descriptor: &[f64]) -> Vec<f64> {
        let mut result = descriptor.to_vec();
        self.normalize(&mut result);
        result
    }
    /// Inverse-transform a normalized descriptor back to original scale.
    pub fn denormalize(&self, normalized: &[f64]) -> Vec<f64> {
        normalized
            .iter()
            .enumerate()
            .map(|(i, &x)| x * self.std[i] + self.mean[i])
            .collect()
    }
}
/// Transfer learning wrapper for a feedforward neural network potential.
///
/// Freezes base layers and allows fine-tuning only the final layers.
#[derive(Debug, Clone)]
pub struct TransferLearning {
    /// The base model (pre-trained).
    pub base: FeedForwardPotential,
    /// Indices of layers that are trainable (not frozen).
    pub fine_tune_layers: Vec<usize>,
    /// Whether the base layers are frozen.
    pub frozen: bool,
}
impl TransferLearning {
    /// Create a transfer learning wrapper.
    ///
    /// `fine_tune_layers` specifies which layer indices remain trainable.
    pub fn new(base: FeedForwardPotential, fine_tune_layers: Vec<usize>) -> Self {
        Self {
            base,
            fine_tune_layers,
            frozen: false,
        }
    }
    /// Freeze all base layers (mark as non-trainable).
    pub fn freeze_base(&mut self) {
        self.frozen = true;
    }
    /// Unfreeze all layers.
    pub fn unfreeze(&mut self) {
        self.frozen = false;
    }
    /// Check if a given layer index is trainable.
    pub fn is_trainable(&self, layer_idx: usize) -> bool {
        if !self.frozen {
            return true;
        }
        self.fine_tune_layers.contains(&layer_idx)
    }
    /// Forward pass through the full network (base + fine-tune layers).
    pub fn forward(&self, descriptor: &[f64]) -> f64 {
        self.base.evaluate_network(descriptor)
    }
    /// Number of total layers.
    pub fn n_layers(&self) -> usize {
        self.base.layers.len()
    }
    /// Number of trainable layers (when frozen).
    pub fn n_trainable(&self) -> usize {
        if !self.frozen {
            return self.base.layers.len();
        }
        self.fine_tune_layers.len()
    }
}
impl TransferLearning {
    /// Fine-tune the trainable layers using stochastic gradient descent
    /// (SGD) on a small labelled dataset.
    ///
    /// `data` is a slice of `(descriptor, target_energy)` pairs.
    /// Only layers listed in `fine_tune_layers` are updated when the model
    /// is frozen; all layers are updated when unfrozen.
    ///
    /// Uses a simple finite-difference gradient estimate with learning rate
    /// `lr` for `epochs` passes over the dataset (proof-of-concept SGD).
    pub fn fine_tune(&mut self, data: &[(Vec<f64>, f64)], lr: f64, epochs: usize) {
        let h = 1e-5_f64;
        for _epoch in 0..epochs {
            for (desc, target) in data {
                let target = *target;
                let e_pred = self.forward(desc);
                let loss = (e_pred - target) * (e_pred - target);
                let _ = loss;
                let n_layers = self.base.layers.len();
                for layer_idx in 0..n_layers {
                    if self.frozen && !self.fine_tune_layers.contains(&layer_idx) {
                        continue;
                    }
                    let n_out = self.base.layers[layer_idx].biases.len();
                    let n_in = self.base.layers[layer_idx].weights[0].len();
                    for i in 0..n_out {
                        for j in 0..n_in {
                            let w0 = self.base.layers[layer_idx].weights[i][j];
                            self.base.layers[layer_idx].weights[i][j] = w0 + h;
                            let ep = self.forward(desc);
                            self.base.layers[layer_idx].weights[i][j] = w0 - h;
                            let em = self.forward(desc);
                            self.base.layers[layer_idx].weights[i][j] = w0;
                            let grad_e = (ep - em) / (2.0 * h);
                            let grad_loss = 2.0 * (e_pred - target) * grad_e;
                            self.base.layers[layer_idx].weights[i][j] -= lr * grad_loss;
                        }
                    }
                    for i in 0..n_out {
                        let b0 = self.base.layers[layer_idx].biases[i];
                        self.base.layers[layer_idx].biases[i] = b0 + h;
                        let ep = self.forward(desc);
                        self.base.layers[layer_idx].biases[i] = b0 - h;
                        let em = self.forward(desc);
                        self.base.layers[layer_idx].biases[i] = b0;
                        let grad_e = (ep - em) / (2.0 * h);
                        let grad_loss = 2.0 * (e_pred - target) * grad_e;
                        self.base.layers[layer_idx].biases[i] -= lr * grad_loss;
                    }
                }
            }
        }
    }
}
/// DimeNet-style directional message passing: computes interaction
/// features for edge (i→j) based on angular information from all
/// other edges (k→j) impinging on atom j.
#[derive(Clone)]
pub struct DimeNetMessageLayer {
    /// Envelope function cutoff.
    pub cutoff: f64,
    /// Dense layer for combining angular features.
    pub dense: DenseLayer,
    /// Number of angular basis functions.
    pub n_angular: usize,
    /// Number of output features.
    pub n_out: usize,
}
impl DimeNetMessageLayer {
    /// Create a DimeNet message layer.
    pub fn new(n_angular: usize, n_out: usize, cutoff: f64) -> Self {
        let mut lcg = Lcg::new(321);
        let scale = (2.0 / (n_angular + n_out) as f64).sqrt();
        let weights: Vec<Vec<f64>> = (0..n_out)
            .map(|_| (0..n_angular).map(|_| lcg.next_f64(scale)).collect())
            .collect();
        let dense = DenseLayer::new(weights, vec![0.0; n_out], Activation::Tanh);
        Self {
            cutoff,
            dense,
            n_angular,
            n_out,
        }
    }
    /// Compute outgoing message from atom i to atom j.
    ///
    /// Aggregates angular features from all atoms k≠i that are neighbours of j.
    pub fn message(&self, positions: &[[f64; 3]], i: usize, j: usize) -> Vec<f64> {
        let n = positions.len();
        let mut agg = vec![0.0_f64; self.n_angular];
        for k in 0..n {
            if k == i || k == j {
                continue;
            }
            let dk = {
                let dx = positions[k][0] - positions[j][0];
                let dy = positions[k][1] - positions[j][1];
                let dz = positions[k][2] - positions[j][2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            if dk >= self.cutoff {
                continue;
            }
            let triplet = AngularTriplet { i: j, j: k, k: i };
            let cos_theta = triplet_cos_angle(positions, triplet);
            let basis = angular_fourier_basis(cos_theta, self.n_angular);
            for (m, &v) in basis.iter().enumerate() {
                agg[m] += v;
            }
        }
        self.dense.forward(&agg)
    }
}
/// A graph attention layer for molecular potentials.
///
/// For each atom i, computes attention weights over its neighbours j:
///   alpha_ij = softmax( LeakyReLU(a^T [Wh_i || Wh_j]) )
///   h_i^new  = sigma( sum_j alpha_ij * W * h_j )
#[derive(Clone)]
pub struct GraphAttentionLayer {
    /// Weight matrix W (n_out × n_in).
    pub w: Vec<Vec<f64>>,
    /// Attention vector a (1 × 2*n_out), stored as flat Vec of length 2*n_out.
    pub a: Vec<f64>,
    /// Number of input features.
    pub n_in: usize,
    /// Number of output features.
    pub n_out: usize,
    /// Cutoff distance.
    pub cutoff: f64,
    /// LeakyReLU slope for negative inputs.
    pub leaky_slope: f64,
}
impl GraphAttentionLayer {
    /// Create a randomly-initialised GAT layer.
    pub fn new(n_in: usize, n_out: usize, cutoff: f64) -> Self {
        let mut lcg = Lcg::new(9001);
        let scale = (2.0 / (n_in + n_out) as f64).sqrt();
        let w: Vec<Vec<f64>> = (0..n_out)
            .map(|_| (0..n_in).map(|_| lcg.next_f64(scale)).collect())
            .collect();
        let a: Vec<f64> = (0..2 * n_out).map(|_| lcg.next_f64(0.1)).collect();
        Self {
            w,
            a,
            n_in,
            n_out,
            cutoff,
            leaky_slope: 0.2,
        }
    }
    fn transform(&self, h: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0_f64; self.n_out];
        for (i, out_i) in out.iter_mut().enumerate().take(self.n_out) {
            let n_w = h.len().min(self.n_in);
            for (&hj, &wij) in h.iter().zip(self.w[i].iter()).take(n_w) {
                *out_i += wij * hj;
            }
        }
        out
    }
    fn leaky_relu(&self, x: f64) -> f64 {
        if x >= 0.0 { x } else { self.leaky_slope * x }
    }
    fn attention_score(&self, hi: &[f64], hj: &[f64]) -> f64 {
        let concat: Vec<f64> = hi.iter().chain(hj.iter()).copied().collect();
        let dot: f64 = concat.iter().zip(self.a.iter()).map(|(x, a)| x * a).sum();
        self.leaky_relu(dot)
    }
    /// Forward pass: update all node features using graph attention.
    pub fn forward(&self, features: &[Vec<f64>], positions: &[[f64; 3]]) -> Vec<Vec<f64>> {
        let n = features.len();
        let transformed: Vec<Vec<f64>> = features.iter().map(|h| self.transform(h)).collect();
        let mut new_features = transformed.clone();
        for i in 0..n {
            let mut neigh_scores: Vec<(usize, f64)> = Vec::new();
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r >= self.cutoff {
                    continue;
                }
                let score = self.attention_score(&transformed[i], &transformed[j]);
                neigh_scores.push((j, score));
            }
            if neigh_scores.is_empty() {
                continue;
            }
            let max_s = neigh_scores
                .iter()
                .map(|&(_, s)| s)
                .fold(f64::NEG_INFINITY, f64::max);
            let exp_scores: Vec<f64> = neigh_scores
                .iter()
                .map(|&(_, s)| (s - max_s).exp())
                .collect();
            let sum_exp: f64 = exp_scores.iter().sum();
            let alphas: Vec<f64> = exp_scores.iter().map(|e| e / sum_exp.max(1e-30)).collect();
            let mut agg = vec![0.0_f64; self.n_out];
            for (idx, &(j, _)) in neigh_scores.iter().enumerate() {
                for k in 0..self.n_out {
                    agg[k] += alphas[idx] * transformed[j][k];
                }
            }
            new_features[i] = agg
                .iter()
                .map(|&x| if x >= 0.0 { x } else { x.exp() - 1.0 })
                .collect();
        }
        new_features
    }
}
/// A stitched PES combining multiple sub-potentials via smooth switching.
///
/// E_total = sum_k w_k(x) * E_k(x)
/// where w_k are Gaussian weights in some reaction coordinate.
pub struct PesStiching {
    /// Energy functions (closures stored as `Box<dyn Fn>`).
    /// For simplicity we store them as trait objects.
    pub n_regions: usize,
    /// Centers of Gaussian switching weights in the reaction coordinate.
    pub centers: Vec<f64>,
    /// Widths of Gaussian switching weights.
    pub widths: Vec<f64>,
    /// Fixed energy offsets for each region (proxy for sub-potential energies).
    pub energy_offsets: Vec<f64>,
}
impl PesStiching {
    /// Create a PES stitching with `n_regions` sub-potentials.
    pub fn new(centers: Vec<f64>, widths: Vec<f64>, energy_offsets: Vec<f64>) -> Self {
        let n = centers.len();
        assert_eq!(widths.len(), n, "widths length mismatch");
        assert_eq!(energy_offsets.len(), n, "energy_offsets length mismatch");
        Self {
            n_regions: n,
            centers,
            widths,
            energy_offsets,
        }
    }
    /// Compute Gaussian weight for region k at reaction coordinate s.
    pub fn weight(&self, k: usize, s: f64) -> f64 {
        let d = (s - self.centers[k]) / self.widths[k];
        (-0.5 * d * d).exp()
    }
    /// Normalised weight (partition of unity).
    pub fn normalised_weight(&self, k: usize, s: f64) -> f64 {
        let w_k = self.weight(k, s);
        let total: f64 = (0..self.n_regions).map(|i| self.weight(i, s)).sum();
        if total < 1e-30 {
            return 1.0 / self.n_regions as f64;
        }
        w_k / total
    }
    /// Blended energy at reaction coordinate s (using energy offsets as proxy).
    pub fn blended_energy(&self, s: f64) -> f64 {
        (0..self.n_regions)
            .map(|k| self.normalised_weight(k, s) * self.energy_offsets[k])
            .sum()
    }
    /// Gradient of blended energy with respect to s.
    pub fn blended_energy_gradient(&self, s: f64) -> f64 {
        let h = 1e-6;
        (self.blended_energy(s + h) - self.blended_energy(s - h)) / (2.0 * h)
    }
}
