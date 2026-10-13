// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Machine-learning potential interface and Behler-Parrinello symmetry
//! function descriptors.
//!
//! This module defines the [`MlPotential`] trait that any ML-based
//! interatomic potential must implement, and provides
//! [`SymmetryFunctionSet`] for computing atom-centered descriptors in
//! the style of Behler and Parrinello (2007).
//!
//! Additionally, this module provides:
//! - [`compute_g2`]: standalone Gaussian radial symmetry function.
//! - [`compute_g4`]: standalone angular symmetry function.
//! - [`BehlerParrinelloDescriptor`]: BP descriptor builder.
//! - [`NeuralNetworkLayer`]: a single dense layer with function-pointer activation.
//! - [`BpNeuralNetworkPotential`]: a complete BP neural network potential.
//! - Activation functions: [`tanh_activation`], [`relu_activation`], [`identity_activation`].

use crate::error::MdError;
use oxiphysics_core::math::Vec3;
use oxiphysics_core::numerical_methods::cholesky;
use std::f64::consts::PI;
use std::sync::OnceLock;

// ---------------------------------------------------------------------------
// ML potential trait
// ---------------------------------------------------------------------------

/// Trait for machine-learning interatomic potentials.
///
/// Implementations receive the full set of atomic positions and species
/// and return the total potential energy together with the per-atom
/// force vectors.
pub trait MlPotential: Send + Sync {
    /// Compute the total energy and per-atom forces.
    ///
    /// # Arguments
    /// * `positions` - Cartesian coordinates of every atom.
    /// * `species`   - Integer species label for every atom (same length as
    ///   `positions`).
    ///
    /// # Returns
    /// A tuple `(energy, forces)` where `forces[i]` is the force on atom `i`.
    fn energy_and_forces(&self, positions: &[Vec3], species: &[u32]) -> (f64, Vec<Vec3>);
}

// ---------------------------------------------------------------------------
// Cutoff function
// ---------------------------------------------------------------------------

/// Smooth cutoff function (cosine taper).
///
/// Returns `0.5 * (1 + cos(pi * r / rc))` when `r < rc`, and `0`
/// otherwise.
#[inline]
pub fn cutoff_function(r: f64, rc: f64) -> f64 {
    if r < rc {
        0.5 * (1.0 + (PI * r / rc).cos())
    } else {
        0.0
    }
}

/// Derivative of [`cutoff_function`] with respect to `r`.
///
/// Returns `-0.5 * (pi / rc) * sin(pi * r / rc)` when `r < rc`, and `0`
/// otherwise.  The taper has zero slope at `rc`, so the derivative is
/// continuous there.
#[inline]
pub fn cutoff_function_derivative(r: f64, rc: f64) -> f64 {
    if r < rc {
        -0.5 * (PI / rc) * (PI * r / rc).sin()
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Symmetry functions (enum)
// ---------------------------------------------------------------------------

/// A single Behler-Parrinello symmetry function.
#[derive(Debug, Clone)]
pub enum SymmetryFunction {
    /// Radial (two-body) symmetry function G2.
    ///
    /// G2 = sum_j exp(-eta * (r_ij - rs)^2) * fc(r_ij)
    G2 {
        /// Width parameter.
        eta: f64,
        /// Shift parameter.
        rs: f64,
        /// Cutoff radius.
        rc: f64,
    },
    /// Angular (three-body) symmetry function G4.
    ///
    /// G4 = 2^(1-zeta) * sum_{j,k!=j} (1 + lambda*cos(theta))^zeta
    ///        * exp(-eta*(r_ij^2 + r_ik^2 + r_jk^2)) * fc(r_ij)*fc(r_ik)*fc(r_jk)
    G4 {
        /// Width parameter.
        eta: f64,
        /// Angular exponent.
        zeta: f64,
        /// +1 or -1 selecting cos / -cos.
        lambda: f64,
        /// Cutoff radius.
        rc: f64,
    },
}

/// A collection of symmetry functions used to build an atom-centered
/// descriptor vector.
#[derive(Debug, Clone)]
pub struct SymmetryFunctionSet {
    /// The symmetry functions in this set.
    pub functions: Vec<SymmetryFunction>,
}

impl SymmetryFunctionSet {
    /// Create a new, empty symmetry function set.
    pub fn new() -> Self {
        Self {
            functions: Vec::new(),
        }
    }

    /// Add a symmetry function to the set.
    pub fn push(&mut self, sf: SymmetryFunction) {
        self.functions.push(sf);
    }

    /// Compute the descriptor vector for atom `center`.
    ///
    /// The returned vector has the same length as `self.functions`.
    pub fn compute_descriptor(
        &self,
        positions: &[Vec3],
        center: usize,
        _species: &[u32],
        cutoff: f64,
    ) -> Vec<f64> {
        let n = positions.len();
        let ri = positions[center];
        let mut descriptor = vec![0.0; self.functions.len()];

        for (idx, sf) in self.functions.iter().enumerate() {
            match sf {
                SymmetryFunction::G2 { eta, rs, rc } => {
                    let rc_use = rc.min(cutoff);
                    for (j, &pos_j) in positions.iter().enumerate() {
                        if j == center {
                            continue;
                        }
                        let rij = (pos_j - ri).norm();
                        if rij < rc_use {
                            let fc = cutoff_function(rij, rc_use);
                            descriptor[idx] += (-eta * (rij - rs).powi(2)).exp() * fc;
                        }
                    }
                }
                SymmetryFunction::G4 {
                    eta,
                    zeta,
                    lambda,
                    rc,
                } => {
                    let rc_use = rc.min(cutoff);
                    let prefactor = 2.0_f64.powf(1.0 - zeta);
                    for j in 0..n {
                        if j == center {
                            continue;
                        }
                        let rij_vec = positions[j] - ri;
                        let rij = rij_vec.norm();
                        if rij >= rc_use {
                            continue;
                        }
                        for k in (j + 1)..n {
                            if k == center {
                                continue;
                            }
                            let rik_vec = positions[k] - ri;
                            let rik = rik_vec.norm();
                            if rik >= rc_use {
                                continue;
                            }
                            let rjk = (positions[k] - positions[j]).norm();
                            if rjk >= rc_use {
                                continue;
                            }
                            // Rounding can push |cos θ| a hair above 1 for (anti)parallel
                            // bonds; clamp so `1 + λ cos θ` never becomes slightly negative
                            // (a NaN under a fractional exponent).
                            let cos_theta = (rij_vec.dot(&rik_vec) / (rij * rik)).clamp(-1.0, 1.0);
                            let angular = (1.0 + lambda * cos_theta).powf(*zeta);
                            let radial = (-eta * (rij.powi(2) + rik.powi(2) + rjk.powi(2))).exp();
                            let fc = cutoff_function(rij, rc_use)
                                * cutoff_function(rik, rc_use)
                                * cutoff_function(rjk, rc_use);
                            descriptor[idx] += prefactor * angular * radial * fc;
                        }
                    }
                }
            }
        }
        descriptor
    }
}

impl Default for SymmetryFunctionSet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Standalone G2 and G4 functions
// ---------------------------------------------------------------------------

/// Compute the G2 radial symmetry function value for one neighbor.
///
/// ```text
/// G2 = exp(-eta * (r_ij - rs)^2) * fc(r_ij, rc)
/// ```
///
/// Returns 0.0 when `r_ij >= rc`.
pub fn compute_g2(r_ij: f64, eta: f64, rs: f64, rc: f64) -> f64 {
    if r_ij >= rc {
        return 0.0;
    }
    (-eta * (r_ij - rs).powi(2)).exp() * cutoff_function(r_ij, rc)
}

/// Compute the G4 angular symmetry function contribution for one triplet.
///
/// ```text
/// G4 = 2^(1-zeta) * (1 + lambda*cos_theta)^zeta
///       * exp(-eta*(r_ij^2 + r_ik^2 + r_jk^2))
///       * fc(r_ij) * fc(r_ik) * fc(r_jk)
/// ```
///
/// Returns 0.0 if any distance is >= rc.
pub fn compute_g4(
    r_ij: f64,
    r_ik: f64,
    r_jk: f64,
    cos_theta: f64,
    eta: f64,
    zeta: f64,
    lambda: f64,
    rc: f64,
) -> f64 {
    if r_ij >= rc || r_ik >= rc || r_jk >= rc {
        return 0.0;
    }
    let prefactor = 2.0_f64.powf(1.0 - zeta);
    let angular = (1.0 + lambda * cos_theta).powf(zeta);
    let radial = (-eta * (r_ij * r_ij + r_ik * r_ik + r_jk * r_jk)).exp();
    let fc = cutoff_function(r_ij, rc) * cutoff_function(r_ik, rc) * cutoff_function(r_jk, rc);
    prefactor * angular * radial * fc
}

// ---------------------------------------------------------------------------
// BehlerParrinelloDescriptor
// ---------------------------------------------------------------------------

/// Parameters for a G2 symmetry function.
#[derive(Debug, Clone)]
pub struct G2Params {
    /// Width parameter η.
    pub eta: f64,
    /// Shift parameter r_s.
    pub rs: f64,
}

/// Parameters for a G4 symmetry function.
#[derive(Debug, Clone)]
pub struct G4Params {
    /// Width parameter η.
    pub eta: f64,
    /// Angular exponent ζ.
    pub zeta: f64,
    /// ±1 angular selector λ.
    pub lambda: f64,
}

/// Builder for Behler-Parrinello atom-centered descriptors.
#[derive(Debug, Clone)]
pub struct BehlerParrinelloDescriptor {
    /// G2 parameter sets.
    pub g2_params: Vec<G2Params>,
    /// G4 parameter sets.
    pub g4_params: Vec<G4Params>,
    /// Cutoff radius (same units as positions).
    pub r_cut: f64,
}

impl BehlerParrinelloDescriptor {
    /// Create a new descriptor builder.
    pub fn new(r_cut: f64) -> Self {
        Self {
            g2_params: Vec::new(),
            g4_params: Vec::new(),
            r_cut,
        }
    }

    /// Add a G2 parameter set.
    pub fn add_g2(&mut self, eta: f64, rs: f64) {
        self.g2_params.push(G2Params { eta, rs });
    }

    /// Add a G4 parameter set.
    pub fn add_g4(&mut self, eta: f64, zeta: f64, lambda: f64) {
        self.g4_params.push(G4Params { eta, zeta, lambda });
    }

    /// Expected length of the descriptor vector.
    pub fn descriptor_len(&self) -> usize {
        self.g2_params.len() + self.g4_params.len()
    }

    /// Compute the descriptor vector for atom `center_idx`.
    ///
    /// `positions` is a slice of `[f64;3]` arrays; `neighbor_idxs` lists
    /// the neighbors to consider (typically all atoms except center).
    pub fn compute(
        &self,
        positions: &[[f64; 3]],
        center_idx: usize,
        neighbor_idxs: &[usize],
    ) -> Vec<f64> {
        let rc = self.r_cut;
        let ri = positions[center_idx];
        let mut desc = vec![0.0f64; self.descriptor_len()];
        let n_g2 = self.g2_params.len();

        // G2 contributions
        for &j in neighbor_idxs {
            if j == center_idx {
                continue;
            }
            let rj = positions[j];
            let r_ij = dist3(ri, rj);
            for (p_idx, p) in self.g2_params.iter().enumerate() {
                desc[p_idx] += compute_g2(r_ij, p.eta, p.rs, rc);
            }
        }

        // G4 contributions (triplets)
        let neighbors: Vec<usize> = neighbor_idxs
            .iter()
            .copied()
            .filter(|&j| j != center_idx)
            .collect();
        for (ni, &j) in neighbors.iter().enumerate() {
            let rj = positions[j];
            let r_ij = dist3(ri, rj);
            if r_ij >= rc {
                continue;
            }
            for &k in &neighbors[ni + 1..] {
                let rk = positions[k];
                let r_ik = dist3(ri, rk);
                if r_ik >= rc {
                    continue;
                }
                let r_jk = dist3(rj, rk);
                if r_jk >= rc {
                    continue;
                }
                // cos theta_ijk
                let cos_theta = {
                    let v_ij = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
                    let v_ik = [rk[0] - ri[0], rk[1] - ri[1], rk[2] - ri[2]];
                    let dot = v_ij[0] * v_ik[0] + v_ij[1] * v_ik[1] + v_ij[2] * v_ik[2];
                    dot / (r_ij * r_ik)
                };
                for (p_idx, p) in self.g4_params.iter().enumerate() {
                    desc[n_g2 + p_idx] +=
                        compute_g4(r_ij, r_ik, r_jk, cos_theta, p.eta, p.zeta, p.lambda, rc);
                }
            }
        }
        desc
    }
}

// ---------------------------------------------------------------------------
// Activation functions and their derivatives
// ---------------------------------------------------------------------------

/// Hyperbolic tangent activation.
pub fn tanh_activation(x: f64) -> f64 {
    x.tanh()
}

/// Derivative of the hyperbolic tangent: 1 - tanh(x)².
pub fn tanh_derivative(x: f64) -> f64 {
    let t = x.tanh();
    1.0 - t * t
}

/// Rectified linear unit activation.
pub fn relu_activation(x: f64) -> f64 {
    x.max(0.0)
}

/// Derivative of the ReLU: 1 if x > 0, else 0.
pub fn relu_derivative(x: f64) -> f64 {
    if x > 0.0 { 1.0 } else { 0.0 }
}

/// Identity (linear) activation.
pub fn identity_activation(x: f64) -> f64 {
    x
}

/// Derivative of the identity activation: always 1.
pub fn identity_derivative(_x: f64) -> f64 {
    1.0
}

/// Automatically derive the analytical derivative for a known activation
/// function pointer. Falls back to `identity_derivative` for unknown functions.
fn derive_activation_for(f: fn(f64) -> f64) -> fn(f64) -> f64 {
    // Cast through raw pointer to avoid the `function_casts_as_integer` lint.
    let f_addr = f as *const () as usize;
    if f_addr == tanh_activation as *const () as usize {
        tanh_derivative
    } else if f_addr == relu_activation as *const () as usize {
        relu_derivative
    } else {
        // Handles identity_activation and any unknown function.
        identity_derivative
    }
}

// ---------------------------------------------------------------------------
// NeuralNetworkLayer
// ---------------------------------------------------------------------------

/// A single fully-connected layer: `y = activation(W * x + b)`.
///
/// Uses a function pointer for the activation, enabling zero-overhead dispatch.
#[derive(Clone)]
pub struct NeuralNetworkLayer {
    /// Weight matrix, row-major: `weights[i][j]` = weight from input `j` to output `i`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector.
    pub biases: Vec<f64>,
    /// Activation function.
    pub activation: fn(f64) -> f64,
    /// Derivative of the activation function.
    pub dactivation: fn(f64) -> f64,
}

impl NeuralNetworkLayer {
    /// Create a new layer. The derivative is auto-detected from the activation
    /// function pointer for the three built-in activations (tanh, relu, identity).
    /// Unknown activation functions fall back to `identity_derivative`.
    pub fn new(weights: Vec<Vec<f64>>, biases: Vec<f64>, activation: fn(f64) -> f64) -> Self {
        let dactivation = derive_activation_for(activation);
        Self {
            weights,
            biases,
            activation,
            dactivation,
        }
    }

    /// Create a new layer with an explicitly provided derivative function.
    /// Use this for custom activation functions not in the built-in set.
    pub fn with_explicit_derivative(
        weights: Vec<Vec<f64>>,
        biases: Vec<f64>,
        activation: fn(f64) -> f64,
        dactivation: fn(f64) -> f64,
    ) -> Self {
        Self {
            weights,
            biases,
            activation,
            dactivation,
        }
    }

    /// Forward pass: `output[i] = activation(sum_j W[i][j] * input[j] + b[i])`.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let n_out = self.biases.len();
        let mut output = self.biases.clone();
        for (i, out_i) in output.iter_mut().enumerate().take(n_out) {
            let n_w = input.len().min(self.weights[i].len());
            for (&inp_j, &w_ij) in input.iter().zip(self.weights[i].iter()).take(n_w) {
                *out_i += w_ij * inp_j;
            }
            *out_i = (self.activation)(*out_i);
        }
        output
    }
}

// ---------------------------------------------------------------------------
// BpNeuralNetworkPotential
// ---------------------------------------------------------------------------

/// A Behler-Parrinello neural network potential.
///
/// Maps BP descriptors through a stack of [`NeuralNetworkLayer`]s to a
/// scalar atomic energy.  Total energy = sum of atomic energies.
pub struct BpNeuralNetworkPotential {
    /// Network layers.
    pub layers: Vec<NeuralNetworkLayer>,
    /// Descriptor builder.
    pub descriptor: BehlerParrinelloDescriptor,
}

impl BpNeuralNetworkPotential {
    /// Create a new potential with the given layers and descriptor.
    pub fn new(layers: Vec<NeuralNetworkLayer>, descriptor: BehlerParrinelloDescriptor) -> Self {
        Self { layers, descriptor }
    }

    /// Evaluate the network on a descriptor vector → scalar energy.
    pub fn energy(&self, descriptor_vec: &[f64]) -> f64 {
        let mut x = descriptor_vec.to_vec();
        for layer in &self.layers {
            x = layer.forward(&x);
        }
        x.iter().sum()
    }

    /// Compute forces on atom `center` by central finite differences.
    ///
    /// Returns a vector of force contributions (length = `neighbors.len() + 1`
    /// for center and each neighbor; here we return forces on the center atom
    /// for each spatial dimension as a `Vec<[f64;3]>` with one entry).
    pub fn forces_by_finite_diff(
        &self,
        positions: &[[f64; 3]],
        center: usize,
        neighbors: &[usize],
        dx: f64,
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0f64; 3]; n];
        let mut pos = positions.to_vec();

        for a in 0..3 {
            pos[center][a] += dx;
            let desc_plus = self.descriptor.compute(&pos, center, neighbors);
            let e_plus = self.energy(&desc_plus);

            pos[center][a] -= 2.0 * dx;
            let desc_minus = self.descriptor.compute(&pos, center, neighbors);
            let e_minus = self.energy(&desc_minus);

            pos[center][a] += dx; // restore
            forces[center][a] = -(e_plus - e_minus) / (2.0 * dx);
        }
        forces
    }
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ---------------------------------------------------------------------------
// Gaussian Approximation Potential (GAP) kernel
// ---------------------------------------------------------------------------

/// Squared-exponential (RBF) kernel for Gaussian Approximation Potentials.
///
/// k(x, x') = σ_f² · exp(-‖x - x'‖² / (2 l²))
pub struct GapKernel {
    /// Signal variance σ_f².
    pub signal_variance: f64,
    /// Length scale l.
    pub length_scale: f64,
}

impl GapKernel {
    /// Create a new GAP kernel.
    pub fn new(signal_variance: f64, length_scale: f64) -> Self {
        Self {
            signal_variance,
            length_scale,
        }
    }

    /// Evaluate the kernel k(x, x').
    pub fn evaluate(&self, x: &[f64], x_prime: &[f64]) -> f64 {
        debug_assert_eq!(x.len(), x_prime.len());
        let sq_dist: f64 = x
            .iter()
            .zip(x_prime.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum();
        self.signal_variance * (-sq_dist / (2.0 * self.length_scale * self.length_scale)).exp()
    }

    /// Build the full kernel matrix K\[i\]\[j\] = k(X\[i\], X\[j\]).
    pub fn kernel_matrix(&self, training_data: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = training_data.len();
        let mut k = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                k[i][j] = self.evaluate(&training_data[i], &training_data[j]);
            }
        }
        k
    }

    /// Compute the kernel vector k*(x*) = \[k(x*, x_1), ..., k(x*, x_n)\].
    pub fn kernel_vector(&self, x_star: &[f64], training_data: &[Vec<f64>]) -> Vec<f64> {
        training_data
            .iter()
            .map(|xi| self.evaluate(x_star, xi))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Kernel Ridge Regression (KRR) ML potential
// ---------------------------------------------------------------------------

/// Kernel Ridge Regression potential (the posterior mean of a Gaussian process).
///
/// Fits: E(x) = Σ_i α_i k(x, x_i)
/// where α are the regression coefficients solved via
/// (K + λI) α = y.
///
/// `K + λI` is symmetric positive definite for a positive-definite kernel and
/// `λ > 0`; it is factorised once as `L Lᵀ` (Cholesky), which gives both the
/// coefficients α and the posterior variance of
/// [`KrrPotential::predict_uncertainty`].
///
/// The Cholesky factor is computed from the hyper-parameters in effect when it
/// is first needed ([`KrrPotential::fit`] computes it eagerly,
/// [`KrrPotential::new`] lazily on the first uncertainty query).  Build a new
/// potential instead of mutating `kernel`, `lambda` or `training_descriptors`
/// of an existing one.
pub struct KrrPotential {
    /// Training descriptors.
    pub training_descriptors: Vec<Vec<f64>>,
    /// Regression coefficients α.
    pub alpha: Vec<f64>,
    /// Kernel function.
    pub kernel: GapKernel,
    /// Regularisation parameter λ.
    pub lambda: f64,
    /// Lower-triangular Cholesky factor `L` of `K + λI`, `None` when the matrix
    /// is not (numerically) symmetric positive definite.  Filled on first use.
    chol: OnceLock<Option<Vec<Vec<f64>>>>,
}

impl KrrPotential {
    /// Create a new KRR potential with pre-computed α coefficients.
    ///
    /// The Cholesky factor of `K + λI` needed for
    /// [`KrrPotential::predict_uncertainty`] is computed lazily on the first
    /// uncertainty query.
    pub fn new(
        training_descriptors: Vec<Vec<f64>>,
        alpha: Vec<f64>,
        kernel: GapKernel,
        lambda: f64,
    ) -> Self {
        Self {
            training_descriptors,
            alpha,
            kernel,
            lambda,
            chol: OnceLock::new(),
        }
    }

    /// Fit a KRR potential to descriptors and energies.
    ///
    /// Solves (K + λI) α = y with a Cholesky factorisation
    /// (O(n³), intended for small training sets) and keeps the factor for
    /// uncertainty predictions.
    ///
    /// # Errors
    /// Returns an error if the number of descriptors and energies differ, the
    /// descriptors have different lengths, the energies, `lambda` or the
    /// kernel hyper-parameters are not finite (or `lambda < 0`, or the kernel
    /// parameters are not positive), or if `K + λI` is not numerically
    /// symmetric positive definite (for example duplicate descriptors with
    /// `lambda = 0`; increase `lambda`).
    pub fn fit(
        descriptors: Vec<Vec<f64>>,
        energies: &[f64],
        kernel: GapKernel,
        lambda: f64,
    ) -> Result<Self, MdError> {
        let n = descriptors.len();
        if n != energies.len() {
            return Err(MdError::General(format!(
                "KRR fit: {n} descriptors but {} energies",
                energies.len()
            )));
        }
        if let Some(first) = descriptors.first()
            && let Some(bad) = descriptors.iter().position(|d| d.len() != first.len())
        {
            return Err(MdError::General(format!(
                "KRR fit: descriptor {bad} has length {} but descriptor 0 has length {}",
                descriptors[bad].len(),
                first.len()
            )));
        }
        if !lambda.is_finite() || lambda < 0.0 {
            return Err(MdError::General(format!(
                "KRR fit: lambda must be finite and non-negative, got {lambda}"
            )));
        }
        let kernel_ok = kernel.signal_variance.is_finite()
            && kernel.signal_variance > 0.0
            && kernel.length_scale.is_finite()
            && kernel.length_scale > 0.0;
        if !kernel_ok {
            return Err(MdError::General(format!(
                "KRR fit: kernel signal variance and length scale must be finite and positive, \
                 got {} and {}",
                kernel.signal_variance, kernel.length_scale
            )));
        }
        if energies.iter().any(|e| !e.is_finite()) {
            return Err(MdError::General(
                "KRR fit: training energies must be finite".to_string(),
            ));
        }

        // Build (K + λI)
        let mut km = kernel.kernel_matrix(&descriptors);
        for (i, row) in km.iter_mut().enumerate() {
            row[i] += lambda;
        }

        let factor = factorise_spd(&km).ok_or_else(|| {
            MdError::General(
                "KRR fit: K + lambda*I is not symmetric positive definite \
                 (duplicate descriptors or lambda too small); increase lambda"
                    .to_string(),
            )
        })?;

        // (K + λI) α = y  <=>  L z = y, Lᵀ α = z
        let z = forward_substitution(&factor, energies);
        let alpha = backward_substitution_transposed(&factor, &z);

        Ok(Self {
            training_descriptors: descriptors,
            alpha,
            kernel,
            lambda,
            chol: OnceLock::from(Some(factor)),
        })
    }

    /// Predict energy for a new descriptor.
    pub fn predict_energy(&self, descriptor: &[f64]) -> f64 {
        let kv = self
            .kernel
            .kernel_vector(descriptor, &self.training_descriptors);
        kv.iter().zip(self.alpha.iter()).map(|(k, a)| k * a).sum()
    }

    /// Cholesky factor of `K + λI`, factorised on first use.
    fn cholesky_factor(&self) -> Option<&Vec<Vec<f64>>> {
        self.chol
            .get_or_init(|| {
                let mut km = self.kernel.kernel_matrix(&self.training_descriptors);
                for (i, row) in km.iter_mut().enumerate() {
                    row[i] += self.lambda;
                }
                factorise_spd(&km)
            })
            .as_ref()
    }

    /// Predict uncertainty (posterior variance of the latent energy).
    ///
    /// σ²(x*) = k(x*,x*) - k*(K + λI)^{-1} k*^T
    ///
    /// evaluated through the Cholesky factor `K + λI = L Lᵀ` as
    /// `k(x*,x*) - ‖v‖²` with `L v = k*`.  The variance does not depend on the
    /// training energies: it equals the signal variance far from the data and
    /// shrinks towards `σ_f² λ / (σ_f² + λ)` at a training point.
    ///
    /// If `K + λI` cannot be factorised (only possible for a potential built
    /// with [`KrrPotential::new`] from inconsistent parameters) nothing can be
    /// learned from the data and the prior variance `k(x*,x*)` is returned.
    pub fn predict_uncertainty(&self, descriptor: &[f64]) -> f64 {
        let k_star_star = self.kernel.evaluate(descriptor, descriptor);
        let kv = self
            .kernel
            .kernel_vector(descriptor, &self.training_descriptors);
        match self.cholesky_factor() {
            Some(factor) if factor.len() == kv.len() => {
                let v = forward_substitution(factor, &kv);
                let explained: f64 = v.iter().map(|x| x * x).sum();
                (k_star_star - explained).max(0.0)
            }
            _ => k_star_star.max(0.0),
        }
    }
}

/// Cholesky factor `L` (`A = L Lᵀ`) of a symmetric positive definite matrix.
///
/// Returns `None` if `a` is not numerically positive definite: the core
/// factorisation fails, or a pivot `L_ii²` is not larger than
/// `n · ε · max_i A_ii` (a pivot at rounding level means the matrix is
/// singular to working precision, and dividing by it would amplify noise),
/// or the factor contains non-finite entries.
fn factorise_spd(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let factor = cholesky(a)?;
    let largest_diagonal = a
        .iter()
        .enumerate()
        .map(|(i, row)| row[i])
        .fold(0.0_f64, f64::max);
    let pivot_floor = n as f64 * f64::EPSILON * largest_diagonal;
    let sound = factor
        .iter()
        .enumerate()
        .all(|(i, row)| row[..=i].iter().all(|x| x.is_finite()) && row[i] * row[i] > pivot_floor);
    sound.then_some(factor)
}

/// Solve `L x = b` for lower-triangular `L` by forward substitution.
fn forward_substitution(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut sum = b[i];
        for j in 0..i {
            sum -= l[i][j] * x[j];
        }
        x[i] = sum / l[i][i];
    }
    x
}

/// Solve `Lᵀ x = b` for lower-triangular `L` by backward substitution.
fn backward_substitution_transposed(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= l[j][i] * x[j];
        }
        x[i] = sum / l[i][i];
    }
    x
}

// ---------------------------------------------------------------------------
// Active learning / uncertainty-based sampling
// ---------------------------------------------------------------------------

/// Active learning criterion: select the sample with the highest uncertainty.
pub struct ActiveLearningSampler {
    /// Minimum uncertainty threshold to trigger new calculation.
    pub uncertainty_threshold: f64,
}

impl ActiveLearningSampler {
    /// Create a new active learning sampler.
    pub fn new(threshold: f64) -> Self {
        Self {
            uncertainty_threshold: threshold,
        }
    }

    /// Select the index with the highest uncertainty from a list of
    /// (descriptor, uncertainty) pairs.
    ///
    /// Returns `None` if all uncertainties are below the threshold.
    pub fn select_most_uncertain(&self, uncertainties: &[f64]) -> Option<usize> {
        let (idx, &max_unc) = uncertainties
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        if max_unc >= self.uncertainty_threshold {
            Some(idx)
        } else {
            None
        }
    }

    /// Query-by-committee uncertainty: standard deviation across committee predictions.
    pub fn committee_uncertainty(predictions: &[f64]) -> f64 {
        if predictions.is_empty() {
            return 0.0;
        }
        let n = predictions.len() as f64;
        let mean = predictions.iter().sum::<f64>() / n;
        let var = predictions
            .iter()
            .map(|&p| (p - mean) * (p - mean))
            .sum::<f64>()
            / n;
        var.sqrt()
    }
}

// ---------------------------------------------------------------------------
// SOAP-like descriptor (simplified Smooth Overlap of Atomic Positions)
// ---------------------------------------------------------------------------

/// Simplified SOAP power spectrum descriptor for a single atom.
///
/// This implements a reduced version: the radial basis is a set of Gaussian
/// functions centred at distances r_n, and the angular part uses spherical
/// harmonics of degree l=0 (just the monopole), giving a 1D radial profile.
pub struct SoapDescriptor {
    /// Radial basis centres r_n.
    pub r_centres: Vec<f64>,
    /// Gaussian width σ.
    pub sigma: f64,
    /// Cutoff radius.
    pub r_cut: f64,
}

impl SoapDescriptor {
    /// Create a new SOAP descriptor.
    pub fn new(r_centres: Vec<f64>, sigma: f64, r_cut: f64) -> Self {
        Self {
            r_centres,
            sigma,
            r_cut,
        }
    }

    /// Build a default SOAP descriptor with evenly spaced radial centres.
    pub fn default_descriptor(r_cut: f64, n_radial: usize) -> Self {
        let dr = r_cut / n_radial as f64;
        let centres: Vec<f64> = (1..=n_radial).map(|i| i as f64 * dr).collect();
        Self::new(centres, 0.5 * dr, r_cut)
    }

    /// Compute the SOAP power spectrum for atom `center_idx`.
    ///
    /// Returns a vector of length `n_radial` (one element per radial basis).
    pub fn compute(&self, positions: &[[f64; 3]], center_idx: usize) -> Vec<f64> {
        let rc = self.r_cut;
        let ri = positions[center_idx];
        let n_r = self.r_centres.len();
        let mut spectrum = vec![0.0f64; n_r];

        for (j, &rj) in positions.iter().enumerate() {
            if j == center_idx {
                continue;
            }
            let r_ij = dist3(ri, rj);
            if r_ij >= rc {
                continue;
            }
            let fc = cutoff_function(r_ij, rc);
            for (n, &r_n) in self.r_centres.iter().enumerate() {
                let exponent = -(r_ij - r_n) * (r_ij - r_n) / (2.0 * self.sigma * self.sigma);
                spectrum[n] += exponent.exp() * fc;
            }
        }

        // Power spectrum: square element-wise
        for s in spectrum.iter_mut() {
            *s = *s * *s;
        }
        spectrum
    }

    /// Number of features in the descriptor.
    pub fn n_features(&self) -> usize {
        self.r_centres.len()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Cutoff tests ───────────────────────────────────────────────────────

    #[test]
    fn test_cutoff_at_zero() {
        let fc = cutoff_function(0.0, 5.0);
        assert!((fc - 1.0).abs() < 1e-12, "fc(0) should be 1, got {fc}");
    }

    #[test]
    fn test_cutoff_at_rc() {
        let fc = cutoff_function(5.0, 5.0);
        assert!(fc.abs() < 1e-12, "fc(rc) should be 0, got {fc}");
    }

    #[test]
    fn test_cutoff_beyond_rc() {
        let fc = cutoff_function(6.0, 5.0);
        assert_eq!(fc, 0.0);
    }

    // ── SymmetryFunctionSet tests ──────────────────────────────────────────

    #[test]
    fn test_g2_two_atoms() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let species = vec![1, 1];

        let mut set = SymmetryFunctionSet::new();
        set.push(SymmetryFunction::G2 {
            eta: 0.5,
            rs: 0.0,
            rc: 5.0,
        });

        let desc = set.compute_descriptor(&positions, 0, &species, 5.0);
        let expected = (-0.5_f64).exp() * cutoff_function(1.0, 5.0);
        assert!(
            (desc[0] - expected).abs() < 1e-12,
            "G2 mismatch: {} vs {}",
            desc[0],
            expected
        );
    }

    #[test]
    fn test_g4_three_atoms() {
        let s = 3.0_f64.sqrt() / 2.0;
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, s, 0.0),
        ];
        let species = vec![1, 1, 1];

        let mut set = SymmetryFunctionSet::new();
        set.push(SymmetryFunction::G4 {
            eta: 0.1,
            zeta: 1.0,
            lambda: 1.0,
            rc: 5.0,
        });

        let desc = set.compute_descriptor(&positions, 0, &species, 5.0);
        assert!(desc[0].is_finite(), "G4 should be finite");
        assert!(desc[0] > 0.0, "G4 should be positive for this geometry");
    }

    #[test]
    fn test_descriptor_dimension() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
        let species = vec![1, 2];

        let mut set = SymmetryFunctionSet::new();
        set.push(SymmetryFunction::G2 {
            eta: 0.5,
            rs: 0.0,
            rc: 6.0,
        });
        set.push(SymmetryFunction::G2 {
            eta: 1.0,
            rs: 1.0,
            rc: 6.0,
        });
        set.push(SymmetryFunction::G4 {
            eta: 0.1,
            zeta: 2.0,
            lambda: 1.0,
            rc: 6.0,
        });

        let desc = set.compute_descriptor(&positions, 0, &species, 6.0);
        assert_eq!(desc.len(), 3);
    }

    #[test]
    fn test_g4_descriptor_is_finite_for_collinear_atoms_with_fractional_zeta() {
        // Atoms on a skew line: about a quarter of the computed cos θ values exceed 1 by
        // rounding, which made `1 - cos θ` slightly negative and `powf(1.5)` NaN.
        let direction = Vec3::new(1.0, 2.0, 3.0).normalize();
        let positions: Vec<Vec3> = (0..12)
            .map(|k| direction * (0.5 + 0.37 * k as f64))
            .collect();
        let species = vec![1; positions.len()];
        let mut set = SymmetryFunctionSet::new();
        set.push(SymmetryFunction::G4 {
            eta: 0.1,
            zeta: 1.5,
            lambda: -1.0,
            rc: 6.0,
        });
        set.push(SymmetryFunction::G4 {
            eta: 0.1,
            zeta: 1.5,
            lambda: 1.0,
            rc: 6.0,
        });
        let desc = set.compute_descriptor(&positions, 0, &species, 6.0);
        assert!(desc.iter().all(|x| x.is_finite()), "descriptor = {desc:?}");
        // All neighbors lie on one side: cos θ = 1, so the λ = -1 term vanishes and the
        // λ = +1 term is strictly positive.
        assert!(desc[0].abs() < 1e-12, "lambda = -1 term = {}", desc[0]);
        assert!(desc[1] > 0.0, "lambda = +1 term = {}", desc[1]);
    }

    // ── compute_g2 / compute_g4 standalone tests ───────────────────────────

    #[test]
    fn test_standalone_g2_zero_at_cutoff() {
        let v = compute_g2(5.0, 0.5, 0.0, 5.0);
        assert!(v.abs() < 1e-10, "G2 should be ~0 at rc, got {v}");
    }

    #[test]
    fn test_standalone_g2_positive_inside() {
        let v = compute_g2(1.0, 0.5, 0.0, 5.0);
        assert!(v > 0.0, "G2 should be positive inside cutoff, got {v}");
    }

    #[test]
    fn test_standalone_g2_zero_beyond_cutoff() {
        let v = compute_g2(6.0, 0.5, 0.0, 5.0);
        assert_eq!(v, 0.0, "G2 should be exactly 0 beyond cutoff");
    }

    #[test]
    fn test_standalone_g4_angular_sensitivity() {
        // Same distances, different cos_theta → different G4 values with lambda=1
        let v_aligned = compute_g4(1.0, 1.0, 1.0, 1.0, 0.1, 1.0, 1.0, 5.0);
        let v_opposing = compute_g4(1.0, 1.0, 1.0, -1.0, 0.1, 1.0, 1.0, 5.0);
        assert!(
            v_aligned > v_opposing,
            "G4 should be larger for aligned (cos=1) than opposing (cos=-1): {v_aligned} vs {v_opposing}"
        );
    }

    #[test]
    fn test_standalone_g4_zero_beyond_cutoff() {
        let v = compute_g4(6.0, 1.0, 1.0, 0.5, 0.1, 1.0, 1.0, 5.0);
        assert_eq!(v, 0.0, "G4 should be 0 when r_ij >= rc");
    }

    // ── BehlerParrinelloDescriptor tests ──────────────────────────────────

    #[test]
    fn test_bp_descriptor_length() {
        let mut desc = BehlerParrinelloDescriptor::new(5.0);
        desc.add_g2(0.5, 0.0);
        desc.add_g2(1.0, 1.0);
        desc.add_g4(0.1, 1.0, 1.0);
        assert_eq!(desc.descriptor_len(), 3);
    }

    #[test]
    fn test_bp_descriptor_positive_g2() {
        let mut d = BehlerParrinelloDescriptor::new(5.0);
        d.add_g2(0.5, 0.0);
        let positions = [[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let neighbors = vec![0, 1];
        let desc = d.compute(&positions, 0, &neighbors);
        assert_eq!(desc.len(), 1);
        assert!(
            desc[0] > 0.0,
            "G2 descriptor should be positive inside cutoff"
        );
    }

    // ── NeuralNetworkLayer tests ───────────────────────────────────────────

    #[test]
    fn test_nn_layer_identity_activation() {
        let layer = NeuralNetworkLayer::new(
            vec![vec![1.0, 2.0], vec![3.0, 4.0]],
            vec![0.5, -0.5],
            identity_activation,
        );
        let out = layer.forward(&[1.0, 1.0]);
        // [1+2+0.5, 3+4-0.5] = [3.5, 6.5]
        assert!((out[0] - 3.5).abs() < 1e-12, "Expected 3.5, got {}", out[0]);
        assert!((out[1] - 6.5).abs() < 1e-12, "Expected 6.5, got {}", out[1]);
    }

    #[test]
    fn test_nn_layer_relu_activation() {
        let layer =
            NeuralNetworkLayer::new(vec![vec![1.0], vec![-1.0]], vec![0.0, 0.0], relu_activation);
        let out = layer.forward(&[2.0]);
        assert!((out[0] - 2.0).abs() < 1e-12, "ReLU(2) should be 2");
        assert!(out[1].abs() < 1e-12, "ReLU(-2) should be 0");
    }

    #[test]
    fn test_nn_layer_tanh_activation() {
        let layer = NeuralNetworkLayer::new(vec![vec![1.0]], vec![0.0], tanh_activation);
        let out = layer.forward(&[0.0]);
        assert!(out[0].abs() < 1e-12, "tanh(0) should be 0, got {}", out[0]);
    }

    // ── Activation function tests ──────────────────────────────────────────

    #[test]
    fn test_activation_functions() {
        assert!((tanh_activation(0.0)).abs() < 1e-12);
        assert!((tanh_activation(1.0) - 1.0_f64.tanh()).abs() < 1e-12);
        assert!((relu_activation(-3.0)).abs() < 1e-12);
        assert!((relu_activation(3.0) - 3.0).abs() < 1e-12);
        assert!((identity_activation(42.0) - 42.0).abs() < 1e-12);
    }

    // ── GAP kernel tests ───────────────────────────────────────────────────

    #[test]
    fn test_gap_kernel_self_similarity() {
        let kernel = GapKernel::new(1.0, 1.0);
        let x = vec![0.5, 1.0, 0.3];
        let k = kernel.evaluate(&x, &x);
        assert!(
            (k - 1.0).abs() < 1e-10,
            "k(x,x) should equal sigma_f^2, got {k}"
        );
    }

    #[test]
    fn test_gap_kernel_decreases_with_distance() {
        let kernel = GapKernel::new(1.0, 1.0);
        let x = vec![0.0, 0.0];
        let y_near = vec![0.1, 0.0];
        let y_far = vec![2.0, 0.0];
        let k_near = kernel.evaluate(&x, &y_near);
        let k_far = kernel.evaluate(&x, &y_far);
        assert!(
            k_near > k_far,
            "kernel should decrease with distance: {k_near} > {k_far}"
        );
    }

    #[test]
    fn test_gap_kernel_matrix_symmetric() {
        let kernel = GapKernel::new(1.0, 0.5);
        let data = vec![vec![0.0f64], vec![1.0], vec![2.0]];
        let km = kernel.kernel_matrix(&data);
        for (i, row) in km.iter().enumerate() {
            for (j, &kmij) in row.iter().enumerate() {
                assert!(
                    (kmij - km[j][i]).abs() < 1e-12,
                    "kernel matrix not symmetric"
                );
            }
        }
    }

    #[test]
    fn test_gap_kernel_vector_length() {
        let kernel = GapKernel::new(1.0, 1.0);
        let training = vec![vec![0.0f64], vec![1.0], vec![2.0]];
        let kv = kernel.kernel_vector(&[0.5], &training);
        assert_eq!(kv.len(), 3);
    }

    // ── KRR potential tests ───────────────────────────────────────────────

    #[test]
    fn test_krr_fit_trivial() {
        // Fit to a constant: all energies equal to 1.0
        let descriptors = vec![vec![0.0f64], vec![1.0], vec![2.0]];
        let energies = vec![1.0, 1.0, 1.0];
        let kernel = GapKernel::new(1.0, 1.0);
        let krr = KrrPotential::fit(descriptors, &energies, kernel, 1e-6).unwrap();
        let pred = krr.predict_energy(&[1.5]);
        assert!((pred - 1.0).abs() < 0.5, "constant fit prediction = {pred}");
    }

    #[test]
    fn test_krr_alpha_length() {
        let descriptors = vec![vec![0.0f64], vec![1.0]];
        let energies = vec![-1.0, -2.0];
        let krr =
            KrrPotential::fit(descriptors, &energies, GapKernel::new(1.0, 1.0), 1e-3).unwrap();
        assert_eq!(krr.alpha.len(), 2);
    }

    #[test]
    fn test_krr_uncertainty_non_negative() {
        let descriptors = vec![vec![0.0f64], vec![1.0], vec![2.0]];
        let energies = vec![-1.0, -2.0, -1.5];
        let krr =
            KrrPotential::fit(descriptors, &energies, GapKernel::new(1.0, 1.0), 1e-3).unwrap();
        let unc = krr.predict_uncertainty(&[1.5]);
        assert!(unc >= 0.0, "uncertainty should be non-negative, got {unc}");
    }

    /// Training set on a 1-D grid with smooth, deterministic energies.
    fn grid_training_set(n: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
        let descriptors: Vec<Vec<f64>> = (0..n).map(|i| vec![0.6 * i as f64]).collect();
        let energies: Vec<f64> = descriptors
            .iter()
            .map(|d| d[0].sin() - 0.3 * d[0])
            .collect();
        (descriptors, energies)
    }

    #[test]
    fn test_krr_fit_solves_the_regularised_normal_equations() {
        let (descriptors, energies) = grid_training_set(8);
        let kernel = GapKernel::new(1.7, 0.9);
        let lambda = 1e-3;
        let km = kernel.kernel_matrix(&descriptors);
        let krr = KrrPotential::fit(descriptors, &energies, kernel, lambda).unwrap();
        for (i, row) in km.iter().enumerate() {
            let lhs: f64 =
                row.iter().zip(&krr.alpha).map(|(k, a)| k * a).sum::<f64>() + lambda * krr.alpha[i];
            assert!(
                (lhs - energies[i]).abs() < 1e-9,
                "row {i}: (K + lambda I) alpha = {lhs}, y = {}",
                energies[i]
            );
        }
    }

    #[test]
    fn test_krr_fit_interpolates_training_energies_for_tiny_lambda() {
        let (descriptors, energies) = grid_training_set(6);
        let krr = KrrPotential::fit(
            descriptors.clone(),
            &energies,
            GapKernel::new(1.0, 1.0),
            1e-12,
        )
        .unwrap();
        for (d, e) in descriptors.iter().zip(&energies) {
            let p = krr.predict_energy(d);
            assert!(
                (p - e).abs() < 1e-5,
                "prediction {p} vs training energy {e}"
            );
        }
    }

    #[test]
    fn test_krr_uncertainty_matches_closed_form_for_one_training_point() {
        // One training point x0 = 0: σ²(x*) = k** - k(x*,x0)² / (k00 + λ).
        let kernel = GapKernel::new(2.0, 1.0);
        let lambda = 0.5;
        let krr = KrrPotential::fit(vec![vec![0.0]], &[3.0], kernel, lambda).unwrap();
        let reference = GapKernel::new(2.0, 1.0);
        let k_star = reference.evaluate(&[1.0], &[0.0]);
        let expected = 2.0 - k_star * k_star / (2.0 + lambda);
        let got = krr.predict_uncertainty(&[1.0]);
        assert!(
            (got - expected).abs() < 1e-12,
            "got {got}, expected {expected}"
        );
        // At the training point itself: σ_f² λ / (σ_f² + λ) = 0.4.
        let at_data = krr.predict_uncertainty(&[0.0]);
        assert!((at_data - 0.4).abs() < 1e-12, "at the data: {at_data}");
    }

    #[test]
    fn test_krr_uncertainty_is_independent_of_the_energy_scale() {
        // Regression: the "variance" was k** - |k*·alpha| = k** - |predicted energy|, which
        // scales with the training targets.  The GP posterior variance does not depend on y.
        let (descriptors, energies) = grid_training_set(7);
        let scaled: Vec<f64> = energies.iter().map(|e| 10.0 * e).collect();
        let kernel = || GapKernel::new(1.3, 0.8);
        let a = KrrPotential::fit(descriptors.clone(), &energies, kernel(), 1e-3).unwrap();
        let b = KrrPotential::fit(descriptors, &scaled, kernel(), 1e-3).unwrap();
        for &x in &[-0.7, 0.3, 1.1, 2.9, 4.4] {
            let ua = a.predict_uncertainty(&[x]);
            let ub = b.predict_uncertainty(&[x]);
            assert!((ua - ub).abs() < 1e-12, "x = {x}: {ua} vs {ub}");
            let ea = a.predict_energy(&[x]);
            let eb = b.predict_energy(&[x]);
            assert!(
                (eb - 10.0 * ea).abs() < 1e-8,
                "predicted energies must scale with the targets"
            );
        }
    }

    #[test]
    fn test_krr_uncertainty_far_from_data_is_the_signal_variance() {
        let (descriptors, energies) = grid_training_set(5);
        let krr =
            KrrPotential::fit(descriptors, &energies, GapKernel::new(2.5, 0.7), 1e-4).unwrap();
        let far = krr.predict_uncertainty(&[1.0e3]);
        assert!((far - 2.5).abs() < 1e-12, "far-field variance = {far}");
    }

    #[test]
    fn test_krr_uncertainty_is_small_near_data_and_grows_away_from_it() {
        let (descriptors, energies) = grid_training_set(6);
        let krr =
            KrrPotential::fit(descriptors, &energies, GapKernel::new(1.0, 0.5), 1e-6).unwrap();
        let on_data = krr.predict_uncertainty(&[0.6]);
        let between = krr.predict_uncertainty(&[0.9]);
        let outside = krr.predict_uncertainty(&[6.5]);
        assert!(on_data < 1e-5, "variance at a training point = {on_data}");
        assert!(between > on_data, "{between} should exceed {on_data}");
        assert!(outside > between, "{outside} should exceed {between}");
        assert!(outside <= 1.0 + 1e-12);
    }

    #[test]
    fn test_krr_new_computes_the_same_uncertainty_lazily() {
        let (descriptors, energies) = grid_training_set(6);
        let fitted = KrrPotential::fit(
            descriptors.clone(),
            &energies,
            GapKernel::new(1.1, 0.9),
            1e-3,
        )
        .unwrap();
        let rebuilt = KrrPotential::new(
            descriptors,
            fitted.alpha.clone(),
            GapKernel::new(1.1, 0.9),
            1e-3,
        );
        for &x in &[0.2, 1.7, 3.3, 9.0] {
            let a = fitted.predict_uncertainty(&[x]);
            let b = rebuilt.predict_uncertainty(&[x]);
            assert!((a - b).abs() < 1e-14, "x = {x}: {a} vs {b}");
            assert!((fitted.predict_energy(&[x]) - rebuilt.predict_energy(&[x])).abs() < 1e-14);
        }
    }

    #[test]
    fn test_krr_unfactorisable_model_reports_the_prior_variance() {
        // Duplicate descriptors with lambda = 0: K is singular, so nothing can be conditioned.
        let krr = KrrPotential::new(
            vec![vec![0.0], vec![0.0]],
            vec![0.0, 0.0],
            GapKernel::new(1.5, 1.0),
            0.0,
        );
        assert!((krr.predict_uncertainty(&[0.3]) - 1.5).abs() < 1e-14);
    }

    #[test]
    fn test_krr_fit_rejects_inconsistent_input() {
        let kernel = || GapKernel::new(1.0, 1.0);
        let d = || vec![vec![0.0], vec![1.0]];
        // length mismatch
        assert!(KrrPotential::fit(d(), &[1.0], kernel(), 1e-3).is_err());
        // ragged descriptors
        assert!(
            KrrPotential::fit(vec![vec![0.0], vec![1.0, 2.0]], &[1.0, 2.0], kernel(), 1e-3)
                .is_err()
        );
        // non-finite energy / lambda
        assert!(KrrPotential::fit(d(), &[1.0, f64::NAN], kernel(), 1e-3).is_err());
        assert!(KrrPotential::fit(d(), &[1.0, 2.0], kernel(), f64::NAN).is_err());
        assert!(KrrPotential::fit(d(), &[1.0, 2.0], kernel(), -1.0).is_err());
        // degenerate kernel hyper-parameters
        assert!(KrrPotential::fit(d(), &[1.0, 2.0], GapKernel::new(1.0, 0.0), 1e-3).is_err());
        assert!(KrrPotential::fit(d(), &[1.0, 2.0], GapKernel::new(-1.0, 1.0), 1e-3).is_err());
    }

    #[test]
    fn test_krr_fit_rejects_singular_kernel_matrix() {
        // Two identical descriptors and no regularisation: K + 0*I is singular.
        let result = KrrPotential::fit(
            vec![vec![0.5], vec![0.5], vec![1.5]],
            &[1.0, 1.0, 2.0],
            GapKernel::new(1.0, 1.0),
            0.0,
        );
        assert!(result.is_err());
        // The same data with a ridge term is fine.
        let result = KrrPotential::fit(
            vec![vec![0.5], vec![0.5], vec![1.5]],
            &[1.0, 1.0, 2.0],
            GapKernel::new(1.0, 1.0),
            1e-3,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_krr_fit_empty_training_set_is_the_prior() {
        let krr = KrrPotential::fit(Vec::new(), &[], GapKernel::new(1.2, 1.0), 1e-3).unwrap();
        assert!(krr.alpha.is_empty());
        assert_eq!(krr.predict_energy(&[0.4]), 0.0);
        assert!((krr.predict_uncertainty(&[0.4]) - 1.2).abs() < 1e-14);
    }

    #[test]
    fn test_cholesky_triangular_solves_recover_the_solution() {
        let a = vec![
            vec![4.0, 1.0, 0.5],
            vec![1.0, 3.0, 0.2],
            vec![0.5, 0.2, 2.0],
        ];
        let factor = factorise_spd(&a).unwrap();
        let b = [1.0, -2.0, 0.5];
        let z = forward_substitution(&factor, &b);
        let x = backward_substitution_transposed(&factor, &z);
        for (i, row) in a.iter().enumerate() {
            let lhs: f64 = row.iter().zip(&x).map(|(aij, xj)| aij * xj).sum();
            assert!((lhs - b[i]).abs() < 1e-12, "row {i}: {lhs} vs {}", b[i]);
        }
        assert!(factorise_spd(&[vec![1.0, 2.0], vec![2.0, 1.0]]).is_none());
        assert!(factorise_spd(&[vec![1.0, 1.0], vec![1.0, 1.0]]).is_none());
    }

    // ── ActiveLearningSampler tests ───────────────────────────────────────

    #[test]
    fn test_active_learning_select_most_uncertain() {
        let sampler = ActiveLearningSampler::new(0.1);
        let uncertainties = vec![0.05, 0.3, 0.15];
        let idx = sampler.select_most_uncertain(&uncertainties);
        assert_eq!(idx, Some(1), "should select index 1 (highest uncertainty)");
    }

    #[test]
    fn test_active_learning_all_below_threshold() {
        let sampler = ActiveLearningSampler::new(0.5);
        let uncertainties = vec![0.1, 0.2, 0.3];
        let idx = sampler.select_most_uncertain(&uncertainties);
        assert_eq!(idx, None, "all below threshold, should return None");
    }

    #[test]
    fn test_active_learning_committee_uncertainty_zero() {
        let preds = vec![1.0, 1.0, 1.0];
        let unc = ActiveLearningSampler::committee_uncertainty(&preds);
        assert!(
            unc.abs() < 1e-10,
            "identical predictions should have 0 uncertainty"
        );
    }

    #[test]
    fn test_active_learning_committee_uncertainty_nonzero() {
        let preds = vec![1.0, 3.0];
        let unc = ActiveLearningSampler::committee_uncertainty(&preds);
        assert!(
            unc > 0.0,
            "committee uncertainty should be positive for different predictions"
        );
    }

    // ── SoapDescriptor tests ─────────────────────────────────────────────

    #[test]
    fn test_soap_descriptor_length() {
        let soap = SoapDescriptor::default_descriptor(5.0, 8);
        let positions = vec![[0.0f64; 3], [1.0, 0.0, 0.0]];
        let desc = soap.compute(&positions, 0);
        assert_eq!(desc.len(), 8, "SOAP descriptor should have 8 features");
    }

    #[test]
    fn test_soap_descriptor_non_negative() {
        let soap = SoapDescriptor::default_descriptor(5.0, 4);
        let positions = vec![[0.0f64; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let desc = soap.compute(&positions, 0);
        for &v in &desc {
            assert!(v >= 0.0, "SOAP power spectrum should be non-negative");
        }
    }

    #[test]
    fn test_soap_descriptor_zero_outside_cutoff() {
        let soap = SoapDescriptor::default_descriptor(1.0, 4);
        // Neighbor at r = 5.0, far beyond cutoff = 1.0
        let positions = vec![[0.0f64; 3], [5.0, 0.0, 0.0]];
        let desc = soap.compute(&positions, 0);
        let total: f64 = desc.iter().sum();
        assert!(
            total < 1e-10,
            "SOAP should be zero when no neighbors in cutoff"
        );
    }

    #[test]
    fn test_soap_n_features() {
        let soap = SoapDescriptor::default_descriptor(3.0, 5);
        assert_eq!(soap.n_features(), 5);
    }

    // ── BpNeuralNetworkPotential tests ────────────────────────────────────

    #[test]
    fn test_bp_nn_energy_single_atom() {
        let mut d = BehlerParrinelloDescriptor::new(5.0);
        d.add_g2(0.5, 0.0);
        d.add_g2(1.0, 1.0);
        // 2-input → 4 → 1 network
        let l1 = NeuralNetworkLayer::new(
            vec![
                vec![0.1, -0.2],
                vec![0.3, 0.1],
                vec![-0.1, 0.4],
                vec![0.2, 0.2],
            ],
            vec![0.0; 4],
            tanh_activation,
        );
        let l2 = NeuralNetworkLayer::new(
            vec![vec![0.1, 0.2, -0.1, 0.3]],
            vec![0.0],
            identity_activation,
        );
        let pot = BpNeuralNetworkPotential::new(vec![l1, l2], d);
        let positions = [[0.0f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let neighbors = vec![0, 1];
        let desc = pot.descriptor.compute(&positions, 0, &neighbors);
        let e = pot.energy(&desc);
        assert!(e.is_finite(), "Energy should be finite, got {e}");
    }

    #[test]
    fn test_bp_nn_forces_finite() {
        let mut d = BehlerParrinelloDescriptor::new(5.0);
        d.add_g2(0.5, 0.0);
        let l1 = NeuralNetworkLayer::new(vec![vec![0.5]], vec![0.0], tanh_activation);
        let l2 = NeuralNetworkLayer::new(vec![vec![1.0]], vec![0.0], identity_activation);
        let pot = BpNeuralNetworkPotential::new(vec![l1, l2], d);
        let positions = [[0.0f64, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let neighbors = vec![0, 1];
        let forces = pot.forces_by_finite_diff(&positions, 0, &neighbors, 1e-5);
        for &f in forces[0].iter() {
            assert!(f.is_finite(), "Force component should be finite, got {f}");
        }
    }
}
