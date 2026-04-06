#![allow(clippy::too_many_arguments)]
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

use oxiphysics_core::math::Vec3;
use std::f64::consts::PI;

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
    #[allow(clippy::needless_range_loop)]
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
                    for j in 0..n {
                        if j == center {
                            continue;
                        }
                        let rij = (positions[j] - ri).norm();
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
                            let cos_theta = rij_vec.dot(&rik_vec) / (rij * rik);
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct G2Params {
    /// Width parameter η.
    pub eta: f64,
    /// Shift parameter r_s.
    pub rs: f64,
}

/// Parameters for a G4 symmetry function.
#[allow(dead_code)]
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new(r_cut: f64) -> Self {
        Self {
            g2_params: Vec::new(),
            g4_params: Vec::new(),
            r_cut,
        }
    }

    /// Add a G2 parameter set.
    #[allow(dead_code)]
    pub fn add_g2(&mut self, eta: f64, rs: f64) {
        self.g2_params.push(G2Params { eta, rs });
    }

    /// Add a G4 parameter set.
    #[allow(dead_code)]
    pub fn add_g4(&mut self, eta: f64, zeta: f64, lambda: f64) {
        self.g4_params.push(G4Params { eta, zeta, lambda });
    }

    /// Expected length of the descriptor vector.
    #[allow(dead_code)]
    pub fn descriptor_len(&self) -> usize {
        self.g2_params.len() + self.g4_params.len()
    }

    /// Compute the descriptor vector for atom `center_idx`.
    ///
    /// `positions` is a slice of `[f64;3]` arrays; `neighbor_idxs` lists
    /// the neighbors to consider (typically all atoms except center).
    #[allow(dead_code)]
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
// Activation functions
// ---------------------------------------------------------------------------

/// Hyperbolic tangent activation.
#[allow(dead_code)]
pub fn tanh_activation(x: f64) -> f64 {
    x.tanh()
}

/// Rectified linear unit activation.
#[allow(dead_code)]
pub fn relu_activation(x: f64) -> f64 {
    x.max(0.0)
}

/// Identity (linear) activation.
#[allow(dead_code)]
pub fn identity_activation(x: f64) -> f64 {
    x
}

// ---------------------------------------------------------------------------
// NeuralNetworkLayer
// ---------------------------------------------------------------------------

/// A single fully-connected layer: `y = activation(W * x + b)`.
///
/// Uses a function pointer for the activation, enabling zero-overhead dispatch.
#[allow(dead_code)]
#[derive(Clone)]
pub struct NeuralNetworkLayer {
    /// Weight matrix, row-major: `weights[i][j]` = weight from input `j` to output `i`.
    pub weights: Vec<Vec<f64>>,
    /// Bias vector.
    pub biases: Vec<f64>,
    /// Activation function.
    pub activation: fn(f64) -> f64,
}

impl NeuralNetworkLayer {
    /// Create a new layer.
    #[allow(dead_code)]
    pub fn new(weights: Vec<Vec<f64>>, biases: Vec<f64>, activation: fn(f64) -> f64) -> Self {
        Self {
            weights,
            biases,
            activation,
        }
    }

    /// Forward pass: `output[i] = activation(sum_j W[i][j] * input[j] + b[i])`.
    #[allow(dead_code)]
    #[allow(clippy::needless_range_loop)]
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let n_out = self.biases.len();
        let mut output = self.biases.clone();
        for i in 0..n_out {
            for j in 0..input.len().min(self.weights[i].len()) {
                output[i] += self.weights[i][j] * input[j];
            }
            output[i] = (self.activation)(output[i]);
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
#[allow(dead_code)]
pub struct BpNeuralNetworkPotential {
    /// Network layers.
    pub layers: Vec<NeuralNetworkLayer>,
    /// Descriptor builder.
    pub descriptor: BehlerParrinelloDescriptor,
}

impl BpNeuralNetworkPotential {
    /// Create a new potential with the given layers and descriptor.
    #[allow(dead_code)]
    pub fn new(layers: Vec<NeuralNetworkLayer>, descriptor: BehlerParrinelloDescriptor) -> Self {
        Self { layers, descriptor }
    }

    /// Evaluate the network on a descriptor vector → scalar energy.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
#[allow(dead_code)]
pub struct GapKernel {
    /// Signal variance σ_f².
    pub signal_variance: f64,
    /// Length scale l.
    pub length_scale: f64,
}

impl GapKernel {
    /// Create a new GAP kernel.
    #[allow(dead_code)]
    pub fn new(signal_variance: f64, length_scale: f64) -> Self {
        Self {
            signal_variance,
            length_scale,
        }
    }

    /// Evaluate the kernel k(x, x').
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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

/// Kernel Ridge Regression potential.
///
/// Fits: E(x) = Σ_i α_i k(x, x_i)
/// where α are the regression coefficients solved via
/// (K + λI) α = y.
#[allow(dead_code)]
pub struct KrrPotential {
    /// Training descriptors.
    pub training_descriptors: Vec<Vec<f64>>,
    /// Regression coefficients α.
    pub alpha: Vec<f64>,
    /// Kernel function.
    pub kernel: GapKernel,
    /// Regularisation parameter λ.
    pub lambda: f64,
}

impl KrrPotential {
    /// Create a new KRR potential with pre-computed α coefficients.
    #[allow(dead_code)]
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
        }
    }

    /// Fit a KRR potential to descriptors and energies.
    ///
    /// Solves (K + λI) α = y using Cholesky-like Gaussian elimination
    /// (naïve O(n³) solver for small training sets).
    #[allow(dead_code)]
    #[allow(clippy::needless_range_loop)]
    pub fn fit(
        descriptors: Vec<Vec<f64>>,
        energies: &[f64],
        kernel: GapKernel,
        lambda: f64,
    ) -> Self {
        let n = descriptors.len();
        assert_eq!(n, energies.len());

        // Build (K + λI)
        let mut km = kernel.kernel_matrix(&descriptors);
        for i in 0..n {
            km[i][i] += lambda;
        }

        // Solve via Gaussian elimination with partial pivoting
        let alpha = gaussian_elimination(&km, energies);

        Self {
            training_descriptors: descriptors,
            alpha,
            kernel,
            lambda,
        }
    }

    /// Predict energy for a new descriptor.
    #[allow(dead_code)]
    pub fn predict_energy(&self, descriptor: &[f64]) -> f64 {
        let kv = self
            .kernel
            .kernel_vector(descriptor, &self.training_descriptors);
        kv.iter().zip(self.alpha.iter()).map(|(k, a)| k * a).sum()
    }

    /// Predict uncertainty (posterior variance).
    ///
    /// σ²(x*) = k(x*, x*) - k*(K + λI)^{-1} k*^T
    #[allow(dead_code)]
    pub fn predict_uncertainty(&self, descriptor: &[f64]) -> f64 {
        let k_star_star = self.kernel.evaluate(descriptor, descriptor);
        let kv = self
            .kernel
            .kernel_vector(descriptor, &self.training_descriptors);
        // Approximate: σ² ≈ k** - kv · α (scaled)
        let cross: f64 = kv.iter().zip(self.alpha.iter()).map(|(k, a)| k * a).sum();
        (k_star_star - cross.abs()).max(0.0)
    }
}

/// Simple Gaussian elimination solver for A x = b (square A).
#[allow(dead_code)]
#[allow(clippy::needless_range_loop)]
fn gaussian_elimination(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut aug: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            row.push(b[i]);
            row
        })
        .collect();

    // Forward elimination
    for col in 0..n {
        // Partial pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..n {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        if pivot.abs() < 1e-15 {
            continue;
        }

        for row in (col + 1)..n {
            let factor = aug[row][col] / pivot;
            for j in col..=n {
                let val = aug[col][j];
                aug[row][j] -= factor * val;
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum -= aug[i][j] * x[j];
        }
        if aug[i][i].abs() > 1e-15 {
            x[i] = sum / aug[i][i];
        }
    }
    x
}

// ---------------------------------------------------------------------------
// Active learning / uncertainty-based sampling
// ---------------------------------------------------------------------------

/// Active learning criterion: select the sample with the highest uncertainty.
#[allow(dead_code)]
pub struct ActiveLearningSampler {
    /// Minimum uncertainty threshold to trigger new calculation.
    pub uncertainty_threshold: f64,
}

impl ActiveLearningSampler {
    /// Create a new active learning sampler.
    #[allow(dead_code)]
    pub fn new(threshold: f64) -> Self {
        Self {
            uncertainty_threshold: threshold,
        }
    }

    /// Select the index with the highest uncertainty from a list of
    /// (descriptor, uncertainty) pairs.
    ///
    /// Returns `None` if all uncertainties are below the threshold.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new(r_centres: Vec<f64>, sigma: f64, r_cut: f64) -> Self {
        Self {
            r_centres,
            sigma,
            r_cut,
        }
    }

    /// Build a default SOAP descriptor with evenly spaced radial centres.
    #[allow(dead_code)]
    pub fn default_descriptor(r_cut: f64, n_radial: usize) -> Self {
        let dr = r_cut / n_radial as f64;
        let centres: Vec<f64> = (1..=n_radial).map(|i| i as f64 * dr).collect();
        Self::new(centres, 0.5 * dr, r_cut)
    }

    /// Compute the SOAP power spectrum for atom `center_idx`.
    ///
    /// Returns a vector of length `n_radial` (one element per radial basis).
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(clippy::needless_range_loop)]
    fn test_gap_kernel_matrix_symmetric() {
        let kernel = GapKernel::new(1.0, 0.5);
        let data = vec![vec![0.0f64], vec![1.0], vec![2.0]];
        let km = kernel.kernel_matrix(&data);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (km[i][j] - km[j][i]).abs() < 1e-12,
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
        let krr = KrrPotential::fit(descriptors, &energies, kernel, 1e-6);
        let pred = krr.predict_energy(&[1.5]);
        assert!((pred - 1.0).abs() < 0.5, "constant fit prediction = {pred}");
    }

    #[test]
    fn test_krr_alpha_length() {
        let descriptors = vec![vec![0.0f64], vec![1.0]];
        let energies = vec![-1.0, -2.0];
        let krr = KrrPotential::fit(descriptors, &energies, GapKernel::new(1.0, 1.0), 1e-3);
        assert_eq!(krr.alpha.len(), 2);
    }

    #[test]
    fn test_krr_uncertainty_non_negative() {
        let descriptors = vec![vec![0.0f64], vec![1.0], vec![2.0]];
        let energies = vec![-1.0, -2.0, -1.5];
        let krr = KrrPotential::fit(descriptors, &energies, GapKernel::new(1.0, 1.0), 1e-3);
        let unc = krr.predict_uncertainty(&[1.5]);
        assert!(unc >= 0.0, "uncertainty should be non-negative, got {unc}");
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

    // ── Gaussian elimination helper ───────────────────────────────────────

    #[test]
    fn test_gaussian_elimination_2x2() {
        // 2x + y = 5, x + 3y = 10 → x=1, y=3
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let b = vec![5.0, 10.0];
        let x = gaussian_elimination(&a, &b);
        assert!((x[0] - 1.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-8, "x[1] = {}", x[1]);
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
