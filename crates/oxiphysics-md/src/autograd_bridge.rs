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
//!   accumulate per-atom force contributions from a neighbour list.
//! - `impl DifferentiableForceField for BpNeuralNetworkPotential`.

use std::collections::HashMap;
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
        for (i, z_val) in z_l.iter_mut().enumerate().take(n_out) {
            let n_weights = layer.weights[i].len().min(n_in);
            for (j, x_val) in x.iter().enumerate().take(n_weights) {
                *z_val += layer.weights[i][j] * x_val;
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
        for (i, &dl) in delta_local.iter().enumerate().take(n_out) {
            let n_w = layer.weights[i].len().min(n_in);
            for (j, dp) in delta_prev.iter_mut().enumerate().take(n_w) {
                *dp += layer.weights[i][j] * dl;
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
// Neighbour lists
// ---------------------------------------------------------------------------

/// Full neighbour lists of an open-boundary (non-periodic) configuration.
///
/// `lists[i]` holds, in ascending order, every atom `j != i` with
/// `|r_ij| < r_cut`; the relation is symmetric.  The atoms are binned into
/// cubic cells of edge `r_cut`, so only the 27 surrounding cells of an atom are
/// searched and the cost is linear in the number of atoms at bounded density
/// (testing all pairs would be quadratic).
///
/// A non-positive or NaN `r_cut` gives empty lists, an infinite `r_cut` makes
/// every atom a neighbour of every other one, and an atom with a non-finite
/// coordinate has no neighbours.
fn open_boundary_neighbor_lists(positions: &[[f64; 3]], r_cut: f64) -> Vec<Vec<usize>> {
    let n = positions.len();
    if r_cut.is_nan() || r_cut <= 0.0 {
        return vec![Vec::new(); n];
    }
    if !r_cut.is_finite() {
        return (0..n)
            .map(|i| (0..n).filter(|&j| j != i).collect())
            .collect();
    }
    let cell_of = |p: &[f64; 3]| -> Option<[i64; 3]> {
        let coordinate = |x: f64| {
            let cell = (x / r_cut).floor();
            // `as` saturates for huge values; NaN and infinities are rejected.
            cell.is_finite().then_some(cell as i64)
        };
        Some([coordinate(p[0])?, coordinate(p[1])?, coordinate(p[2])?])
    };
    let cells: Vec<Option<[i64; 3]>> = positions.iter().map(cell_of).collect();
    let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (atom, cell) in cells.iter().enumerate() {
        if let Some(cell) = cell {
            grid.entry(*cell).or_default().push(atom);
        }
    }
    let cutoff_sq = r_cut * r_cut;
    let mut lists: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (atom, (cell, list)) in cells.iter().zip(lists.iter_mut()).enumerate() {
        let Some(cell) = cell else {
            continue;
        };
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    let key = [
                        cell[0].saturating_add(dx),
                        cell[1].saturating_add(dy),
                        cell[2].saturating_add(dz),
                    ];
                    let Some(members) = grid.get(&key) else {
                        continue;
                    };
                    for &other in members {
                        if other != atom {
                            let d = sub3(positions[other], positions[atom]);
                            if dot3(d, d) < cutoff_sq {
                                list.push(other);
                            }
                        }
                    }
                }
            }
        }
        // Cells beyond the `i64` range saturate and may be visited twice.
        list.sort_unstable();
        list.dedup();
    }
    lists
}

// ---------------------------------------------------------------------------
// Small 3-vector helpers
// ---------------------------------------------------------------------------

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// `acc += s * v`.
#[inline]
fn add_scaled3(acc: &mut [f64; 3], v: [f64; 3], s: f64) {
    acc[0] += s * v[0];
    acc[1] += s * v[1];
    acc[2] += s * v[2];
}

// ---------------------------------------------------------------------------
// BP descriptor backward pass
// ---------------------------------------------------------------------------

/// One G4 parameter set with a non-zero `∂E/∂G4` and its prefactor `2^(1-ζ)`.
struct ActiveG4 {
    de_dg: f64,
    eta: f64,
    zeta: f64,
    lambda: f64,
    prefactor: f64,
}

/// Accumulate the per-atom force contributions from the descriptor gradient of
/// atom `center` into `forces`.
///
/// `ddesc[k]` = ∂E/∂G_k (from the NN backward pass), ordered G2 sets first and
/// G4 sets after them.  For every atom `a` the force
/// `-∂E_center/∂r_a[α]` is **added** to `forces[a][α]`, so the contributions of
/// all centres can be summed into one caller-owned buffer.
///
/// `neighbors` lists, each atom at most once, the atoms that can lie within
/// the descriptor cutoff of `center`, typically the neighbour list of `center`;
/// atoms at or beyond the cutoff, `center` itself and out-of-range indices
/// contribute nothing, so a complete neighbour list gives the same result as
/// passing every atom.  The
/// cost per call is `O(m·P₂ + m²·P₄)` for `m` neighbours with `P₂`/`P₄` G2/G4
/// parameter sets.  The geometry of each triplet (distances, cutoff factors and
/// their derivatives, `cos θ` and its gradient) is evaluated once and shared by
/// all G4 sets.
///
/// Nothing is done if `center` is out of range or `forces` is shorter than
/// `positions`.
pub fn bp_descriptor_forces(
    descriptor: &BehlerParrinelloDescriptor,
    ddesc: &[f64],
    positions: &[[f64; 3]],
    center: usize,
    neighbors: &[usize],
    forces: &mut [[f64; 3]],
) {
    let Some(&ri) = positions.get(center) else {
        return;
    };
    if forces.len() < positions.len() {
        return;
    }
    let rc = descriptor.r_cut;
    let n_g2 = descriptor.g2_params.len();

    // Neighbours inside the cutoff: (index, displacement r_j - r_i, distance).
    let near: Vec<(usize, [f64; 3], f64)> = neighbors
        .iter()
        .filter(|&&j| j != center)
        .filter_map(|&j| {
            let rj = *positions.get(j)?;
            let dr = sub3(rj, ri);
            let r = dot3(dr, dr).sqrt();
            (r < rc && r >= 1e-15).then_some((j, dr, r))
        })
        .collect();

    // --- G2 contributions ---
    // dG2/dr_ij = exp(-η(r-rs)²) · [dfc/dr + (-2η(r-rs)) · fc]; with
    // dr_ij/dr_i = -r_hat and dr_ij/dr_j = +r_hat the force is
    // f_i += (∂E/∂G2 · dG2/dr_ij) r_hat, f_j -= (∂E/∂G2 · dG2/dr_ij) r_hat.
    for &(j, dr, r_ij) in &near {
        let r_hat = scale3(dr, 1.0 / r_ij);
        let fc_val = cutoff(r_ij, rc);
        let dfc = dcutoff_dr(r_ij, rc);
        let mut coeff = 0.0_f64;
        for (p_idx, p) in descriptor.g2_params.iter().enumerate() {
            let de_dg = ddesc.get(p_idx).copied().unwrap_or(0.0);
            if de_dg == 0.0 {
                continue;
            }
            let diff = r_ij - p.rs;
            let gauss = (-p.eta * diff * diff).exp();
            coeff += de_dg * gauss * (dfc + (-2.0 * p.eta * diff) * fc_val);
        }
        add_scaled3(&mut forces[center], r_hat, coeff);
        add_scaled3(&mut forces[j], r_hat, -coeff);
    }

    // --- G4 contributions ---
    let active: Vec<ActiveG4> = descriptor
        .g4_params
        .iter()
        .enumerate()
        .filter_map(|(g4_idx, p)| {
            let de_dg = ddesc.get(n_g2 + g4_idx).copied().unwrap_or(0.0);
            (de_dg != 0.0).then(|| ActiveG4 {
                de_dg,
                eta: p.eta,
                zeta: p.zeta,
                lambda: p.lambda,
                prefactor: 2.0_f64.powf(1.0 - p.zeta),
            })
        })
        .collect();
    if active.is_empty() {
        return;
    }

    for (ni, &(j, dr_ij, r_ij)) in near.iter().enumerate() {
        for &(k, dr_ik, r_ik) in &near[ni + 1..] {
            let dr_jk = sub3(positions[k], positions[j]);
            let r_jk = dot3(dr_jk, dr_jk).sqrt();
            if r_jk >= rc || r_jk < 1e-15 {
                continue;
            }

            // Geometry shared by all G4 parameter sets of this triplet.
            let cos_theta = dot3(dr_ij, dr_ik) / (r_ij * r_ik);
            let fc_ij = cutoff(r_ij, rc);
            let fc_ik = cutoff(r_ik, rc);
            let fc_jk = cutoff(r_jk, rc);
            let fc_prod = fc_ij * fc_ik * fc_jk;
            let dfc_ij_term = dcutoff_dr(r_ij, rc) * fc_ik * fc_jk;
            let dfc_ik_term = fc_ij * dcutoff_dr(r_ik, rc) * fc_jk;
            let dfc_jk_term = fc_ij * fc_ik * dcutoff_dr(r_jk, rc);
            let sum_r2 = r_ij * r_ij + r_ik * r_ik + r_jk * r_jk;

            // Unit vectors of the three pair distances.
            let r_hat_ij = scale3(dr_ij, 1.0 / r_ij);
            let r_hat_ik = scale3(dr_ik, 1.0 / r_ik);
            let r_hat_jk = scale3(dr_jk, 1.0 / r_jk);

            // cosθ = (r_j - r_i)·(r_k - r_i) / (r_ij r_ik), hence
            //   ∂cosθ/∂r_j = dr_ik/(r_ij r_ik) - cosθ dr_ij/r_ij²,
            //   ∂cosθ/∂r_k = dr_ij/(r_ij r_ik) - cosθ dr_ik/r_ik²,
            //   ∂cosθ/∂r_i = -(∂cosθ/∂r_j + ∂cosθ/∂r_k).
            let inv_rr = 1.0 / (r_ij * r_ik);
            let dcos_j: [f64; 3] =
                std::array::from_fn(|a| dr_ik[a] * inv_rr - cos_theta * dr_ij[a] / (r_ij * r_ij));
            let dcos_k: [f64; 3] =
                std::array::from_fn(|a| dr_ij[a] * inv_rr - cos_theta * dr_ik[a] / (r_ik * r_ik));
            let dcos_i: [f64; 3] = std::array::from_fn(|a| -(dcos_j[a] + dcos_k[a]));

            // Sum ∂E/∂G4 · ∂G4/∂(r_ij, r_ik, r_jk, cosθ) over the parameter
            // sets; the forces are linear in these four coefficients.
            let mut coeff_rij = 0.0_f64;
            let mut coeff_rik = 0.0_f64;
            let mut coeff_rjk = 0.0_f64;
            let mut coeff_cos = 0.0_f64;
            for p in &active {
                let one_plus_lcos = 1.0 + p.lambda * cos_theta;
                // Guard against a zero / negative base for non-integer zeta.
                if one_plus_lcos <= 0.0 {
                    continue;
                }
                let angular = one_plus_lcos.powf(p.zeta);
                let exp_val = (-p.eta * sum_r2).exp();
                let radial = p.prefactor * angular * exp_val;
                // ∂T/∂cosθ = prefactor·ζ·(1+λcosθ)^(ζ-1)·λ·exp_val·fc_prod
                coeff_cos += p.de_dg
                    * p.prefactor
                    * p.zeta
                    * one_plus_lcos.powf(p.zeta - 1.0)
                    * p.lambda
                    * exp_val
                    * fc_prod;
                // ∂T/∂r_xy holding the other distances and cosθ fixed: the
                // Gaussian factor and the cutoff of that distance.
                coeff_rij += p.de_dg * radial * (-2.0 * p.eta * r_ij * fc_prod + dfc_ij_term);
                coeff_rik += p.de_dg * radial * (-2.0 * p.eta * r_ik * fc_prod + dfc_ik_term);
                coeff_rjk += p.de_dg * radial * (-2.0 * p.eta * r_jk * fc_prod + dfc_jk_term);
            }

            // Cartesian forces f_a = -∂E/∂x_a with
            //   ∂r_ij/∂x_i = -r_hat_ij, ∂r_ij/∂x_j = +r_hat_ij,
            //   ∂r_ik/∂x_i = -r_hat_ik, ∂r_ik/∂x_k = +r_hat_ik,
            //   ∂r_jk/∂x_j = -r_hat_jk, ∂r_jk/∂x_k = +r_hat_jk.
            let f_center = &mut forces[center];
            add_scaled3(f_center, r_hat_ij, coeff_rij);
            add_scaled3(f_center, r_hat_ik, coeff_rik);
            add_scaled3(f_center, dcos_i, -coeff_cos);
            let f_j = &mut forces[j];
            add_scaled3(f_j, r_hat_ij, -coeff_rij);
            add_scaled3(f_j, r_hat_jk, coeff_rjk);
            add_scaled3(f_j, dcos_j, -coeff_cos);
            let f_k = &mut forces[k];
            add_scaled3(f_k, r_hat_ik, -coeff_rik);
            add_scaled3(f_k, r_hat_jk, -coeff_rjk);
            add_scaled3(f_k, dcos_k, -coeff_cos);
        }
    }
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
        // Only atoms within the cutoff of a centre enter its descriptor, so one
        // cell-list pass replaces the all-atoms neighbour list (and with it the
        // cubic scaling of the triplet loops).
        let neighbor_lists = open_boundary_neighbor_lists(&pos_arr, self.descriptor.r_cut);

        let mut total_energy = 0.0_f64;
        let mut atom_forces = vec![[0.0_f64; 3]; n];

        // We accumulate parameter gradients.  They are summed across all
        // center atoms (one NN evaluation per center).
        let n_params = self.n_params();
        let mut all_param_grads: Vec<f64> = vec![0.0; n_params];

        for (center, neighbors) in neighbor_lists.iter().enumerate() {
            // Forward: descriptor → NN with tape
            let desc = self.descriptor.compute(&pos_arr, center, neighbors);
            let (e_atom, tape) = nn_forward_with_tape(&self.layers, &desc);
            total_energy += e_atom;

            // NN backward: get dE/d(descriptor) and dE/d(parameters)
            let nn_grads = nn_backward(&self.layers, &tape);

            // Descriptor backward: per-atom force contributions, accumulated
            // in place (they already carry the sign f = -∂E/∂x).
            bp_descriptor_forces(
                &self.descriptor,
                &nn_grads.ddesc,
                &pos_arr,
                center,
                neighbors,
                &mut atom_forces,
            );

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

        let forces = atom_forces
            .iter()
            .map(|f| Vec3::new(f[0], f[1], f[2]))
            .collect();
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic generator of uniform numbers in `[0, 1)`.
    struct Lcg(u64);

    impl Lcg {
        fn next_in(&mut self, lo: f64, hi: f64) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let unit = (self.0 >> 11) as f64 / (1u64 << 53) as f64;
            lo + (hi - lo) * unit
        }
    }

    /// All-pairs reference: ascending indices of the atoms closer than `r_cut`.
    fn brute_force_lists(positions: &[[f64; 3]], r_cut: f64) -> Vec<Vec<usize>> {
        (0..positions.len())
            .map(|i| {
                (0..positions.len())
                    .filter(|&j| {
                        let d = sub3(positions[j], positions[i]);
                        j != i && dot3(d, d) < r_cut * r_cut
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn cell_lists_match_all_pairs_for_random_clouds() {
        // Negative coordinates, dense and sparse clouds, cutoffs larger and
        // smaller than the cloud.
        for (seed, n, half_edge, r_cut) in [
            (1_u64, 150_usize, 5.0_f64, 2.5_f64),
            (2, 200, 3.0, 1.0),
            (3, 120, 15.0, 4.0),
            (4, 80, 1.5, 5.0),
        ] {
            let mut rng = Lcg(seed);
            let positions: Vec<[f64; 3]> = (0..n)
                .map(|_| {
                    [
                        rng.next_in(-half_edge, half_edge),
                        rng.next_in(-half_edge, half_edge),
                        rng.next_in(-half_edge, half_edge),
                    ]
                })
                .collect();
            let lists = open_boundary_neighbor_lists(&positions, r_cut);
            assert_eq!(lists, brute_force_lists(&positions, r_cut), "seed {seed}");
            let pairs: usize = lists.iter().map(Vec::len).sum();
            assert!(pairs > 0, "seed {seed}: the cloud must have neighbours");
        }
    }

    #[test]
    fn cell_lists_are_symmetric_sorted_and_exclude_the_atom_itself() {
        let mut rng = Lcg(9);
        let positions: Vec<[f64; 3]> = (0..100)
            .map(|_| {
                [
                    rng.next_in(0.0, 8.0),
                    rng.next_in(0.0, 8.0),
                    rng.next_in(0.0, 8.0),
                ]
            })
            .collect();
        let lists = open_boundary_neighbor_lists(&positions, 2.0);
        for (i, list) in lists.iter().enumerate() {
            assert!(list.windows(2).all(|w| w[0] < w[1]), "atom {i}: {list:?}");
            assert!(!list.contains(&i));
            for &j in list {
                assert!(lists[j].contains(&i), "{j} must list {i} as well");
            }
        }
    }

    #[test]
    fn cell_lists_use_a_strict_cutoff_and_do_not_wrap_around() {
        // A pair at exactly r_cut is outside (the descriptor cutoff is
        // open), a pair just inside is a neighbour, and atoms on opposite
        // faces of the bounding box are not periodic images of each other.
        let positions = [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 2.999_999, 0.0],
            [-3.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
        ];
        let lists = open_boundary_neighbor_lists(&positions, 3.0);
        assert_eq!(lists[0], vec![2]);
        assert_eq!(lists[1], Vec::<usize>::new());
        assert_eq!(lists[2], vec![0]);
        assert_eq!(lists[3], Vec::<usize>::new());
        assert_eq!(lists[4], Vec::<usize>::new());
    }

    #[test]
    fn cell_lists_handle_degenerate_input() {
        assert!(open_boundary_neighbor_lists(&[], 3.0).is_empty());
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.7, 0.0]];
        for r_cut in [0.0, -1.0, f64::NAN] {
            let lists = open_boundary_neighbor_lists(&positions, r_cut);
            assert_eq!(lists, vec![Vec::<usize>::new(); 3], "r_cut {r_cut}");
        }
        let everyone = open_boundary_neighbor_lists(&positions, f64::INFINITY);
        assert_eq!(everyone, vec![vec![1, 2], vec![0, 2], vec![0, 1]]);

        // Non-finite coordinates never pair up; the others are unaffected.
        let with_bad = [
            [0.0, 0.0, 0.0],
            [f64::NAN, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [f64::INFINITY, 0.0, 0.0],
        ];
        let lists = open_boundary_neighbor_lists(&with_bad, 1.0);
        assert_eq!(lists, vec![vec![2], vec![], vec![0], vec![]]);

        // Huge coordinates saturate the cell index instead of overflowing.
        let huge = [[1e300, 0.0, 0.0], [1e300, 0.5, 0.0], [-1e300, 0.0, 0.0]];
        let lists = open_boundary_neighbor_lists(&huge, 1.0);
        assert_eq!(lists, vec![vec![1], vec![0], vec![]]);
    }

    #[test]
    fn descriptor_forces_ignore_atoms_beyond_the_cutoff() {
        // Two atoms inside the cutoff of the centre and one outside: listing
        // the outside atom changes nothing.
        let mut descriptor = BehlerParrinelloDescriptor::new(2.0);
        descriptor.add_g2(0.4, 0.0);
        descriptor.add_g4(0.1, 2.0, 1.0);
        let positions = [
            [0.0, 0.0, 0.0],
            [1.1, 0.2, 0.0],
            [0.3, 1.2, 0.4],
            [5.0, 0.0, 0.0],
        ];
        let ddesc = [0.7, -0.4];
        let mut inside_only = vec![[0.0_f64; 3]; 4];
        bp_descriptor_forces(
            &descriptor,
            &ddesc,
            &positions,
            0,
            &[1, 2],
            &mut inside_only,
        );
        let mut with_far = vec![[0.0_f64; 3]; 4];
        bp_descriptor_forces(
            &descriptor,
            &ddesc,
            &positions,
            0,
            &[1, 2, 3],
            &mut with_far,
        );
        assert_eq!(inside_only, with_far);
        assert_eq!(inside_only[3], [0.0; 3]);
        // Momentum conservation of the descriptor forces.
        for a in 0..3 {
            let total: f64 = inside_only.iter().map(|f| f[a]).sum();
            assert!(total.abs() < 1e-14, "net force component {a}: {total}");
        }
        assert!(inside_only[1].iter().any(|c| c.abs() > 1e-6));
    }
}
