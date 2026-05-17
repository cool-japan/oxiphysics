// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Analytical gradient computation (backpropagation) for the
//! Behler-Parrinello neural network potential.
//!
//! This module provides:
//! - [`DifferentiableForceField`]: trait for potentials that expose analytical
//!   parameter gradients alongside energy and forces.
//! - [`NnTape`] / [`nn_forward_with_tape`] / [`NnGrads`] / [`nn_backward`]:
//!   forward-with-tape and backward pass for the dense NN.
//! - [`bp_descriptor_forces`]: chain rule through the BP G2/G4 descriptors to
//!   produce per-atom force contributions.
//! - `impl DifferentiableForceField for BpNeuralNetworkPotential`.

use std::f64::consts::PI;

use oxiphysics_core::math::Vec3;

use crate::ml_potential::{
    BehlerParrinelloDescriptor, BpNeuralNetworkPotential, MlPotential, NeuralNetworkLayer,
};

// ---------------------------------------------------------------------------
// DifferentiableForceField trait
// ---------------------------------------------------------------------------

/// Extension of [`MlPotential`] that exposes analytical gradients.
///
/// Implementors must return the total energy, per-atom analytical forces, and
/// a flattened vector of gradients w.r.t. every trainable parameter (weights
/// first, then biases, layer by layer).
pub trait DifferentiableForceField: MlPotential {
    /// Compute energy, per-atom forces (analytical), and flattened parameter gradients.
    ///
    /// The parameter-gradient vector is ordered as: for each layer in order,
    /// all weights row-major, followed by all biases.
    fn energy_forces_param_grads(
        &self,
        positions: &[Vec3],
        species: &[u32],
    ) -> (f64, Vec<Vec3>, Vec<f64>);

    /// Total number of trainable parameters.
    fn n_params(&self) -> usize;

    /// Apply an in-place parameter update: `param[i] += step[i]`.
    fn apply_param_step(&mut self, step: &[f64]);
}

// ---------------------------------------------------------------------------
// NN forward-with-tape
// ---------------------------------------------------------------------------

/// Recorded intermediate values from a NN forward pass, needed for backprop.
pub struct NnTape {
    /// Pre-activation values `z[ℓ]` at layer ℓ (0-indexed).
    pub z: Vec<Vec<f64>>,
    /// Post-activation values `a[ℓ]`; `a[0]` is the input descriptor,
    /// `a[l+1]` is the output after layer `l`.
    pub a: Vec<Vec<f64>>,
}

/// Run the neural network forward pass, recording `z` and `a` at every layer.
///
/// Returns `(energy, tape)` where `energy = Σ a_last`.
#[allow(clippy::needless_range_loop)]
pub fn nn_forward_with_tape(layers: &[NeuralNetworkLayer], input: &[f64]) -> (f64, NnTape) {
    let mut tape = NnTape {
        z: Vec::with_capacity(layers.len()),
        a: Vec::with_capacity(layers.len() + 1),
    };
    tape.a.push(input.to_vec());
    let mut x = input.to_vec();

    for layer in layers {
        let n_out = layer.biases.len();
        let n_in = x.len();
        let mut z_l = layer.biases.clone();
        for i in 0..n_out {
            let n_weights = layer.weights[i].len().min(n_in);
            for j in 0..n_weights {
                z_l[i] += layer.weights[i][j] * x[j];
            }
        }
        let a_l: Vec<f64> = z_l.iter().map(|&zj| (layer.activation)(zj)).collect();
        tape.z.push(z_l);
        tape.a.push(a_l.clone());
        x = a_l;
    }

    let energy: f64 = x.iter().sum();
    (energy, tape)
}

// ---------------------------------------------------------------------------
// NN backward pass
// ---------------------------------------------------------------------------

/// Gradients produced by the NN backward pass.
pub struct NnGrads {
    /// `dw[ℓ][i][j]` = ∂E/∂W\[ℓ\]\[i\]\[j\].
    pub dw: Vec<Vec<Vec<f64>>>,
    /// `db[ℓ][i]` = ∂E/∂b\[ℓ\]\[i\].
    pub db: Vec<Vec<f64>>,
    /// Gradient w.r.t. the input descriptor: ∂E/∂desc\[k\].
    pub ddesc: Vec<f64>,
}

/// Backpropagate through the neural network.
///
/// `tape` must have been produced by [`nn_forward_with_tape`] on the same
/// `layers`.  Returns weight/bias gradients and the descriptor gradient.
#[allow(clippy::needless_range_loop)]
pub fn nn_backward(layers: &[NeuralNetworkLayer], tape: &NnTape) -> NnGrads {
    let n_layers = layers.len();
    let mut dw: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_layers);
    let mut db: Vec<Vec<f64>> = Vec::with_capacity(n_layers);

    // Output gradient: dE/da[L] = [1, 1, ..., 1]  (energy = sum of outputs)
    let last_output_size = tape.a[n_layers].len();
    let mut delta: Vec<f64> = vec![1.0; last_output_size];

    for l in (0..n_layers).rev() {
        let layer = &layers[l];
        let z_l = &tape.z[l];
        let a_prev = &tape.a[l]; // activation output of prev layer (= input to layer l)

        // Element-wise: delta_local[i] = delta[i] * f'(z_l[i])
        let delta_local: Vec<f64> = delta
            .iter()
            .zip(z_l.iter())
            .map(|(&d, &z)| d * (layer.dactivation)(z))
            .collect();

        let n_out = delta_local.len();
        let n_in = a_prev.len();

        // Weight gradients: dW[i][j] = delta_local[i] * a_prev[j]
        let dw_l: Vec<Vec<f64>> = (0..n_out)
            .map(|i| (0..n_in).map(|j| delta_local[i] * a_prev[j]).collect())
            .collect();

        // Bias gradients: db[i] = delta_local[i]
        let db_l = delta_local.clone();

        // Propagate delta to previous layer: delta_prev[j] = Σ_i W[i][j] * delta_local[i]
        let mut delta_prev = vec![0.0f64; n_in];
        for i in 0..n_out {
            let n_w = layer.weights[i].len().min(n_in);
            for j in 0..n_w {
                delta_prev[j] += layer.weights[i][j] * delta_local[i];
            }
        }

        dw.push(dw_l);
        db.push(db_l);
        delta = delta_prev;
    }

    // Built in reverse; flip back to forward order.
    dw.reverse();
    db.reverse();

    NnGrads {
        dw,
        db,
        ddesc: delta,
    }
}

// ---------------------------------------------------------------------------
// Cutoff derivative helper
// ---------------------------------------------------------------------------

/// Derivative of the cosine cutoff function: `dfc/dr`.
///
/// Returns `-(π / (2·rc)) · sin(π·r/rc)` for `r < rc`, else `0.0`.
#[inline]
fn dcutoff_dr(r: f64, rc: f64) -> f64 {
    if r < rc {
        -(PI / (2.0 * rc)) * (PI * r / rc).sin()
    } else {
        0.0
    }
}

/// Cosine cutoff: `0.5 * (1 + cos(π·r/rc))` for `r < rc`, else `0.0`.
#[inline]
fn cutoff(r: f64, rc: f64) -> f64 {
    if r < rc {
        0.5 * (1.0 + (PI * r / rc).cos())
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// BP descriptor backward pass
// ---------------------------------------------------------------------------

/// Compute per-atom force contributions from the descriptor gradient of atom
/// `center`.
///
/// `ddesc[k]` = ∂E/∂G_k (from the NN backward pass).
/// Returns forces `f_a` for every atom `a` in `positions` such that
/// `f_a[α] += -∂E_center/∂r_a[α]`.
pub fn bp_descriptor_forces(
    descriptor: &BehlerParrinelloDescriptor,
    ddesc: &[f64],
    positions: &[[f64; 3]],
    center: usize,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let rc = descriptor.r_cut;
    let ri = positions[center];
    let n_g2 = descriptor.g2_params.len();

    // --- G2 contributions ---
    for (p_idx, p) in descriptor.g2_params.iter().enumerate() {
        let de_dg = ddesc[p_idx];
        if de_dg == 0.0 {
            continue;
        }
        for j in 0..n {
            if j == center {
                continue;
            }
            let rj = positions[j];
            let dr = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
            let r_ij = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
            if r_ij >= rc || r_ij < 1e-15 {
                continue;
            }
            let r_hat = [dr[0] / r_ij, dr[1] / r_ij, dr[2] / r_ij];

            // dG2/dr_ij = exp(-η(r-rs)²) * [dfc/dr + (-2η(r-rs)) * fc]
            let eta = p.eta;
            let rs = p.rs;
            let diff = r_ij - rs;
            let gauss = (-eta * diff * diff).exp();
            let fc_val = cutoff(r_ij, rc);
            let dfc = dcutoff_dr(r_ij, rc);
            let dg2_drij = gauss * (dfc + (-2.0 * eta * diff) * fc_val);

            // Force from chain rule: f_i_α += -dE/dG2 * dG2/dr_ij * (-r_hat)
            // f_j_α += -dE/dG2 * dG2/dr_ij * (+r_hat)
            // (dr_ij/d(r_i) = -r_hat, dr_ij/d(r_j) = +r_hat)
            let coeff = de_dg * dg2_drij;
            for alpha in 0..3 {
                forces[center][alpha] += coeff * r_hat[alpha]; // -(-r_hat) = +r_hat
                forces[j][alpha] -= coeff * r_hat[alpha]; // -(+r_hat)
            }
        }
    }

    // --- G4 contributions ---
    // Collect neighbors once
    let neighbors: Vec<usize> = (0..n).filter(|&j| j != center).collect();

    for (g4_idx, p) in descriptor.g4_params.iter().enumerate() {
        let desc_idx = n_g2 + g4_idx;
        let de_dg = ddesc[desc_idx];
        if de_dg == 0.0 {
            continue;
        }
        let eta = p.eta;
        let zeta = p.zeta;
        let lambda = p.lambda;
        let prefactor = 2.0_f64.powf(1.0 - zeta);

        for (ni, &j) in neighbors.iter().enumerate() {
            let rj = positions[j];
            let dr_ij = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
            let r_ij = (dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2]).sqrt();
            if r_ij >= rc || r_ij < 1e-15 {
                continue;
            }

            for &k in &neighbors[ni + 1..] {
                let rk = positions[k];
                let dr_ik = [rk[0] - ri[0], rk[1] - ri[1], rk[2] - ri[2]];
                let r_ik = (dr_ik[0] * dr_ik[0] + dr_ik[1] * dr_ik[1] + dr_ik[2] * dr_ik[2]).sqrt();
                if r_ik >= rc || r_ik < 1e-15 {
                    continue;
                }

                let dr_jk = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
                let r_jk = (dr_jk[0] * dr_jk[0] + dr_jk[1] * dr_jk[1] + dr_jk[2] * dr_jk[2]).sqrt();
                if r_jk >= rc || r_jk < 1e-15 {
                    continue;
                }

                // Dot product r_ij · r_ik  (for cos_theta)
                let dot_ij_ik = dr_ij[0] * dr_ik[0] + dr_ij[1] * dr_ik[1] + dr_ij[2] * dr_ik[2];
                let cos_theta = dot_ij_ik / (r_ij * r_ik);
                let one_plus_lcos = 1.0 + lambda * cos_theta;

                // Guard against zero / negative base for non-integer zeta
                if one_plus_lcos <= 0.0 {
                    continue;
                }

                let angular = one_plus_lcos.powf(zeta);
                let exp_val = (-eta * (r_ij * r_ij + r_ik * r_ik + r_jk * r_jk)).exp();
                let fc_ij = cutoff(r_ij, rc);
                let fc_ik = cutoff(r_ik, rc);
                let fc_jk = cutoff(r_jk, rc);
                let fc_prod = fc_ij * fc_ik * fc_jk;

                let dfc_ij = dcutoff_dr(r_ij, rc);
                let dfc_ik = dcutoff_dr(r_ik, rc);
                let dfc_jk = dcutoff_dr(r_jk, rc);

                // ∂T/∂cosθ = prefactor * zeta * (1+λcosθ)^(ζ-1) * λ * exp_val * fc_prod
                let dt_dcostheta =
                    prefactor * zeta * one_plus_lcos.powf(zeta - 1.0) * lambda * exp_val * fc_prod;

                // ∂T/∂r_ij  (holding r_ik, r_jk, cosθ constant — radial + cutoff only)
                // = prefactor * angular * exp_val * (-2η r_ij) * fc_prod
                // + prefactor * angular * exp_val * dfc_ij/dr * fc_ik * fc_jk
                let dt_drij = prefactor
                    * angular
                    * (exp_val * (-2.0 * eta * r_ij) * fc_prod + exp_val * dfc_ij * fc_ik * fc_jk);

                // ∂T/∂r_ik  (holding r_ij, r_jk, cosθ constant)
                let dt_drik = prefactor
                    * angular
                    * (exp_val * (-2.0 * eta * r_ik) * fc_prod + exp_val * fc_ij * dfc_ik * fc_jk);

                // ∂T/∂r_jk  (holding r_ij, r_ik, cosθ constant)
                let dt_drjk = prefactor
                    * angular
                    * (exp_val * (-2.0 * eta * r_jk) * fc_prod + exp_val * fc_ij * fc_ik * dfc_jk);

                // Scale by dE/dG4
                let coeff_rij = de_dg * dt_drij;
                let coeff_rik = de_dg * dt_drik;
                let coeff_rjk = de_dg * dt_drjk;
                let coeff_cos = de_dg * dt_dcostheta;

                // r_hat vectors
                let r_hat_ij = [dr_ij[0] / r_ij, dr_ij[1] / r_ij, dr_ij[2] / r_ij];
                let r_hat_ik = [dr_ik[0] / r_ik, dr_ik[1] / r_ik, dr_ik[2] / r_ik];
                let r_hat_jk = [dr_jk[0] / r_jk, dr_jk[1] / r_jk, dr_jk[2] / r_jk];

                // Cartesian force via chain rule: f_a = -∂E/∂x_a
                //
                // x_a contributes through:
                //   r_ij (if a=i or a=j), r_ik (if a=i or a=k), r_jk (if a=j or a=k), cosθ
                //
                // ∂r_ij/∂x_i = -r_hat_ij, ∂r_ij/∂x_j = +r_hat_ij
                // ∂r_ik/∂x_i = -r_hat_ik, ∂r_ik/∂x_k = +r_hat_ik
                // ∂r_jk/∂x_j = -r_hat_jk, ∂r_jk/∂x_k = +r_hat_jk
                //
                // ∂cosθ/∂x_i = -(r_ij_vec/(r_ij r_ik) + r_ik_vec/(r_ij r_ik))    (vector)
                //             = -(dr_ij/(r_ij r_ik) + dr_ik/(r_ij r_ik))
                // ∂cosθ/∂x_j = +r_ik_vec/(r_ij r_ik) - cosθ * r_ij_vec/r_ij²
                //             = dr_ik/(r_ij r_ik) - cosθ * dr_ij/r_ij²
                //   Wait: more precisely:
                //   cosθ = (r_j - r_i)·(r_k - r_i) / (r_ij * r_ik)
                //   ∂cosθ/∂r_j_α = (r_k - r_i)_α / (r_ij * r_ik)
                //                 - [(r_j - r_i)·(r_k - r_i) / (r_ij³ r_ik)] * (r_j - r_i)_α
                //                = dr_ik[α] / (r_ij * r_ik) - cosθ * dr_ij[α] / r_ij²
                //   ∂cosθ/∂r_k_α = dr_ij[α] / (r_ij * r_ik) - cosθ * dr_ik[α] / r_ik²
                //   ∂cosθ/∂r_i_α = -(∂cosθ/∂r_j_α + ∂cosθ/∂r_k_α)

                for alpha in 0..3 {
                    // --- Radial contributions ---
                    // r_ij: ∂/∂x_i = -r_hat_ij, ∂/∂x_j = +r_hat_ij
                    forces[center][alpha] += coeff_rij * r_hat_ij[alpha]; // -(-r_hat)
                    forces[j][alpha] -= coeff_rij * r_hat_ij[alpha];

                    // r_ik: ∂/∂x_i = -r_hat_ik, ∂/∂x_k = +r_hat_ik
                    forces[center][alpha] += coeff_rik * r_hat_ik[alpha];
                    forces[k][alpha] -= coeff_rik * r_hat_ik[alpha];

                    // r_jk: ∂/∂x_j = -r_hat_jk, ∂/∂x_k = +r_hat_jk
                    forces[j][alpha] += coeff_rjk * r_hat_jk[alpha];
                    forces[k][alpha] -= coeff_rjk * r_hat_jk[alpha];

                    // --- Angular (cos_theta) contributions ---
                    // ∂cosθ/∂x_j_α = dr_ik[α]/(r_ij*r_ik) - cosθ*dr_ij[α]/r_ij²
                    let dcostheta_xj =
                        dr_ik[alpha] / (r_ij * r_ik) - cos_theta * dr_ij[alpha] / (r_ij * r_ij);
                    // ∂cosθ/∂x_k_α = dr_ij[α]/(r_ij*r_ik) - cosθ*dr_ik[α]/r_ik²
                    let dcostheta_xk =
                        dr_ij[alpha] / (r_ij * r_ik) - cos_theta * dr_ik[alpha] / (r_ik * r_ik);
                    // ∂cosθ/∂x_i_α = -(∂cosθ/∂x_j_α + ∂cosθ/∂x_k_α)
                    let dcostheta_xi = -(dcostheta_xj + dcostheta_xk);

                    // Force = -∂E/∂x = -coeff_cos * ∂cosθ/∂x
                    forces[center][alpha] -= coeff_cos * dcostheta_xi;
                    forces[j][alpha] -= coeff_cos * dcostheta_xj;
                    forces[k][alpha] -= coeff_cos * dcostheta_xk;
                }
            }
        }
    }

    forces
}

// ---------------------------------------------------------------------------
// impl DifferentiableForceField for BpNeuralNetworkPotential
// ---------------------------------------------------------------------------

impl MlPotential for BpNeuralNetworkPotential {
    fn energy_and_forces(&self, positions: &[Vec3], species: &[u32]) -> (f64, Vec<Vec3>) {
        let (e, forces, _) = self.energy_forces_param_grads(positions, species);
        (e, forces)
    }
}

impl DifferentiableForceField for BpNeuralNetworkPotential {
    fn energy_forces_param_grads(
        &self,
        positions: &[Vec3],
        species: &[u32],
    ) -> (f64, Vec<Vec3>, Vec<f64>) {
        let n = positions.len();
        let _ = species; // species not used by this descriptor

        // Convert Vec3 positions to [f64;3] arrays
        let pos_arr: Vec<[f64; 3]> = positions.iter().map(|p| [p[0], p[1], p[2]]).collect();
        let neighbors: Vec<usize> = (0..n).collect();

        let mut total_energy = 0.0_f64;
        let mut forces = vec![Vec3::zeros(); n];

        // We accumulate parameter gradients.  They are summed across all
        // center atoms (one NN evaluation per center).
        let n_params = self.n_params();
        let mut all_param_grads: Vec<f64> = vec![0.0; n_params];

        for center in 0..n {
            // Forward: descriptor → NN with tape
            let desc = self.descriptor.compute(&pos_arr, center, &neighbors);
            let (e_atom, tape) = nn_forward_with_tape(&self.layers, &desc);
            total_energy += e_atom;

            // NN backward: get dE/d(descriptor) and dE/d(parameters)
            let nn_grads = nn_backward(&self.layers, &tape);

            // Descriptor backward: per-atom force contributions
            let atom_forces =
                bp_descriptor_forces(&self.descriptor, &nn_grads.ddesc, &pos_arr, center);
            for i in 0..n {
                // forces already carry the sign: f = -∂E/∂x
                forces[i] += Vec3::new(atom_forces[i][0], atom_forces[i][1], atom_forces[i][2]);
            }

            // Accumulate (sum) parameter gradients across atoms
            let mut idx = 0;
            for (dw_l, db_l) in nn_grads.dw.iter().zip(nn_grads.db.iter()) {
                for row in dw_l {
                    for &dw in row {
                        all_param_grads[idx] += dw;
                        idx += 1;
                    }
                }
                for &db in db_l {
                    all_param_grads[idx] += db;
                    idx += 1;
                }
            }
        }

        (total_energy, forces, all_param_grads)
    }

    fn n_params(&self) -> usize {
        self.layers
            .iter()
            .map(|l| l.weights.iter().map(|r| r.len()).sum::<usize>() + l.biases.len())
            .sum()
    }

    fn apply_param_step(&mut self, step: &[f64]) {
        let mut idx = 0;
        for layer in &mut self.layers {
            for row in &mut layer.weights {
                for w in row.iter_mut() {
                    *w += step[idx];
                    idx += 1;
                }
            }
            for b in &mut layer.biases {
                *b += step[idx];
                idx += 1;
            }
        }
    }
}
