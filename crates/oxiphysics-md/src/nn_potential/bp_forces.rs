// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Analytic energies and forces of the Behler-Parrinello feed-forward potential.
//!
//! The total energy is `E = Σ_i NN(G_i)` where `G_i` is the symmetry-function
//! descriptor of atom `i`.  The force on atom `a` is
//!
//! ```text
//! F_a = − Σ_i Σ_k (∂NN/∂G_{i,k}) (∂G_{i,k}/∂r_a)
//! ```
//!
//! with `∂NN/∂G` from back-propagation through the dense layers
//! ([`NetworkTrace`]) and `∂G/∂r` from the closed-form derivatives of the G2
//! and G4 symmetry functions below.  Only atoms within the cutoff of `i`
//! contribute to `G_i`, so the work is organised around neighbor lists
//! ([`neighbor_lists`]) instead of full-system energy re-evaluations.

use super::backprop::NetworkTrace;
use super::locality::{Neighbor, neighbor_lists};
use super::types::{DenseLayer, FeedForwardPotential};
use crate::ml_potential::{
    SymmetryFunction, SymmetryFunctionSet, cutoff_function, cutoff_function_derivative,
};
use oxiphysics_core::math::Vec3;

/// Pairs closer than this have no defined direction; their terms are skipped.
const MIN_DISTANCE: f64 = 1e-12;

/// Geometry of the triplet `(i; j, k)` entering a G4 symmetry function.
struct Triplet {
    /// `r_j - r_i`.
    a: Vec3,
    /// `r_k - r_i`.
    b: Vec3,
    /// `r_k - r_j`.
    c: Vec3,
    rij: f64,
    rik: f64,
    rjk: f64,
    /// `cos θ_jik`, clamped to `[-1, 1]` so rounding cannot push `1 + λ cos θ` below zero.
    cos_theta: f64,
}

impl Triplet {
    /// The triplet of centre `i` with neighbors `nj` and `nk`, or `None` if
    /// any of the three distances reaches `rc` or the centre coincides with a
    /// neighbor (the angle is undefined).
    fn new(nj: &Neighbor, nk: &Neighbor, positions: &[Vec3], rc: f64) -> Option<Self> {
        if nj.r >= rc || nk.r >= rc || nj.r < MIN_DISTANCE || nk.r < MIN_DISTANCE {
            return None;
        }
        let a = Vec3::from(nj.dr);
        let b = Vec3::from(nk.dr);
        let c = positions[nk.j] - positions[nj.j];
        let rjk = c.norm();
        if rjk >= rc {
            return None;
        }
        let cos_theta = (a.dot(&b) / (nj.r * nk.r)).clamp(-1.0, 1.0);
        Some(Self {
            a,
            b,
            c,
            rij: nj.r,
            rik: nk.r,
            rjk,
            cos_theta,
        })
    }

    /// G4 term `2^(1-ζ) (1 + λ cos θ)^ζ exp(-η (r_ij² + r_ik² + r_jk²)) f_c(r_ij) f_c(r_ik) f_c(r_jk)`.
    fn value(&self, eta: f64, zeta: f64, lambda: f64, rc: f64) -> f64 {
        let prefactor = 2.0_f64.powf(1.0 - zeta);
        let angular = (1.0 + lambda * self.cos_theta).powf(zeta);
        let radial = (-eta * (self.rij.powi(2) + self.rik.powi(2) + self.rjk.powi(2))).exp();
        let fc = cutoff_function(self.rij, rc)
            * cutoff_function(self.rik, rc)
            * cutoff_function(self.rjk, rc);
        prefactor * angular * radial * fc
    }

    /// Gradient of the G4 term with respect to the positions of `(i, j, k)`.
    fn gradient(&self, eta: f64, zeta: f64, lambda: f64, rc: f64) -> [Vec3; 3] {
        let prefactor = 2.0_f64.powf(1.0 - zeta);
        let base = 1.0 + lambda * self.cos_theta;
        let angular = base.powf(zeta);
        // d/d(cos θ) of (1 + λ cos θ)^ζ; the term vanishes (for ζ > 1) or is
        // not differentiable (ζ < 1) where the base reaches zero.
        let d_angular = if base > 0.0 {
            zeta * lambda * base.powf(zeta - 1.0)
        } else {
            0.0
        };
        let radial = (-eta * (self.rij.powi(2) + self.rik.powi(2) + self.rjk.powi(2))).exp();
        let (f_ij, f_ik, f_jk) = (
            cutoff_function(self.rij, rc),
            cutoff_function(self.rik, rc),
            cutoff_function(self.rjk, rc),
        );
        let (d_ij, d_ik, d_jk) = (
            cutoff_function_derivative(self.rij, rc),
            cutoff_function_derivative(self.rik, rc),
            cutoff_function_derivative(self.rjk, rc),
        );
        let common = prefactor * radial;
        let t_rij = common * angular * (-2.0 * eta * self.rij * f_ij + d_ij) * f_ik * f_jk;
        let t_rik = common * angular * f_ij * (-2.0 * eta * self.rik * f_ik + d_ik) * f_jk;
        let t_rjk = common * angular * f_ij * f_ik * (-2.0 * eta * self.rjk * f_jk + d_jk);
        let t_cos = common * d_angular * f_ij * f_ik * f_jk;

        let a_hat = self.a / self.rij;
        let b_hat = self.b / self.rik;
        let c_hat = if self.rjk > MIN_DISTANCE {
            self.c / self.rjk
        } else {
            // dT/dr_jk vanishes with r_jk (f_c'(0) = 0), so the direction is irrelevant.
            Vec3::zeros()
        };
        let dcos_da = (b_hat - a_hat * self.cos_theta) / self.rij;
        let dcos_db = (a_hat - b_hat * self.cos_theta) / self.rik;
        let g_j = a_hat * t_rij + dcos_da * t_cos - c_hat * t_rjk;
        let g_k = b_hat * t_rik + dcos_db * t_cos + c_hat * t_rjk;
        // The term depends only on position differences.
        let g_i = -(g_j + g_k);
        [g_i, g_j, g_k]
    }
}

/// Derivative with respect to `r` of the G2 term `exp(-η (r - r_s)²) f_c(r)`.
fn g2_radial_derivative(r: f64, eta: f64, rs: f64, rc: f64) -> f64 {
    let gauss = (-eta * (r - rs).powi(2)).exp();
    gauss * (-2.0 * eta * (r - rs) * cutoff_function(r, rc) + cutoff_function_derivative(r, rc))
}

/// Descriptor `G_i` of one atom from its neighbor list.
///
/// Accumulates exactly the terms of
/// [`SymmetryFunctionSet::compute_descriptor`], in the same order, but only
/// over the atoms within the cutoff.
fn descriptor_from_neighbors(
    set: &SymmetryFunctionSet,
    positions: &[Vec3],
    neighbors: &[Neighbor],
    cutoff: f64,
) -> Vec<f64> {
    let mut descriptor = vec![0.0; set.functions.len()];
    for (idx, sf) in set.functions.iter().enumerate() {
        match sf {
            SymmetryFunction::G2 { eta, rs, rc } => {
                let rc_use = rc.min(cutoff);
                for nb in neighbors {
                    if nb.r < rc_use {
                        descriptor[idx] +=
                            (-eta * (nb.r - rs).powi(2)).exp() * cutoff_function(nb.r, rc_use);
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
                for (a, nj) in neighbors.iter().enumerate() {
                    for nk in &neighbors[a + 1..] {
                        if let Some(triplet) = Triplet::new(nj, nk, positions, rc_use) {
                            descriptor[idx] += triplet.value(*eta, *zeta, *lambda, rc_use);
                        }
                    }
                }
            }
        }
    }
    descriptor
}

/// Add `-Σ_k weights[k] ∂G_k/∂r_a` to `forces[a]` for the descriptor of atom `i`.
///
/// `weights[k]` is `∂E_i/∂G_{i,k}` from back-propagation.
fn accumulate_descriptor_forces(
    set: &SymmetryFunctionSet,
    positions: &[Vec3],
    i: usize,
    neighbors: &[Neighbor],
    cutoff: f64,
    weights: &[f64],
    forces: &mut [Vec3],
) {
    for (idx, sf) in set.functions.iter().enumerate() {
        let w = weights.get(idx).copied().unwrap_or(0.0);
        if w == 0.0 {
            continue;
        }
        match sf {
            SymmetryFunction::G2 { eta, rs, rc } => {
                let rc_use = rc.min(cutoff);
                for nb in neighbors {
                    if nb.r < rc_use && nb.r > MIN_DISTANCE {
                        let slope = g2_radial_derivative(nb.r, *eta, *rs, rc_use);
                        // ∂G/∂r_j = slope r̂_ij and ∂G/∂r_i = -slope r̂_ij.
                        let pull = Vec3::from(nb.dr) * (w * slope / nb.r);
                        forces[nb.j] -= pull;
                        forces[i] += pull;
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
                for (a, nj) in neighbors.iter().enumerate() {
                    for nk in &neighbors[a + 1..] {
                        if let Some(triplet) = Triplet::new(nj, nk, positions, rc_use) {
                            let [g_i, g_j, g_k] = triplet.gradient(*eta, *zeta, *lambda, rc_use);
                            forces[i] -= g_i * w;
                            forces[nj.j] -= g_j * w;
                            forces[nk.j] -= g_k * w;
                        }
                    }
                }
            }
        }
    }
}

/// Atomic energy `Σ NN(G)` and its gradient `∂E/∂G` with respect to the descriptor.
fn atomic_energy_and_descriptor_gradient(
    layers: &[DenseLayer],
    descriptor: &[f64],
) -> (f64, Vec<f64>) {
    let trace = NetworkTrace::forward(layers, descriptor);
    let output = trace.output();
    let energy = output.iter().sum();
    let grad = trace.backward(layers, &vec![1.0; output.len()]);
    (energy, grad)
}

impl FeedForwardPotential {
    /// Radius beyond which no symmetry function of this potential has support.
    fn interaction_radius(&self) -> f64 {
        self.symmetry_functions
            .functions
            .iter()
            .map(|sf| match sf {
                SymmetryFunction::G2 { rc, .. } | SymmetryFunction::G4 { rc, .. } => {
                    rc.min(self.cutoff)
                }
            })
            .fold(0.0_f64, f64::max)
    }

    /// Total energy and analytic forces in one pass over cutoff neighbor lists.
    ///
    /// The cost is `O(N · k²)` for `k` neighbors per atom (descriptors,
    /// back-propagation and force assembly), instead of the `O(N³)` of `6N`
    /// full-system energy evaluations.  The `species` labels are ignored, as in
    /// [`SymmetryFunctionSet::compute_descriptor`].
    pub(super) fn energy_and_forces_analytic(
        &self,
        positions: &[Vec3],
        _species: &[u32],
    ) -> (f64, Vec<Vec3>) {
        let n = positions.len();
        let raw: Vec<[f64; 3]> = positions.iter().map(|p| [p[0], p[1], p[2]]).collect();
        let lists = neighbor_lists(&raw, self.interaction_radius());
        let mut energy = 0.0;
        let mut forces = vec![Vec3::zeros(); n];
        for (i, neighbors) in lists.iter().enumerate() {
            let descriptor = descriptor_from_neighbors(
                &self.symmetry_functions,
                positions,
                neighbors,
                self.cutoff,
            );
            let (atomic_energy, weights) =
                atomic_energy_and_descriptor_gradient(&self.layers, &descriptor);
            energy += atomic_energy;
            accumulate_descriptor_forces(
                &self.symmetry_functions,
                positions,
                i,
                neighbors,
                self.cutoff,
                &weights,
                &mut forces,
            );
        }
        (energy, forces)
    }

    /// Forces from central finite differences of the total energy (step `fd_step`).
    ///
    /// This displaces every coordinate by `±fd_step` and re-evaluates the
    /// whole-system energy (`6N` evaluations, `O(N³)` overall).  It is
    /// intended as an independent reference for validating the analytic forces
    /// returned by [`MlPotential::energy_and_forces`](crate::ml_potential::MlPotential::energy_and_forces),
    /// not for production dynamics.
    pub fn finite_difference_forces(&self, positions: &[Vec3], species: &[u32]) -> Vec<Vec3> {
        let n = positions.len();
        let h = self.fd_step;
        let mut forces = vec![Vec3::zeros(); n];
        let mut pos_plus = positions.to_vec();
        let mut pos_minus = positions.to_vec();
        for i in 0..n {
            for a in 0..3_usize {
                pos_plus[i][a] += h;
                pos_minus[i][a] -= h;
                let e_plus = self.total_energy(&pos_plus, species);
                let e_minus = self.total_energy(&pos_minus, species);
                forces[i][a] = -(e_plus - e_minus) / (2.0 * h);
                pos_plus[i] = positions[i];
                pos_minus[i] = positions[i];
            }
        }
        forces
    }
}

#[cfg(test)]
mod tests {
    use super::super::types::Activation;
    use super::*;
    use crate::ml_potential::{MlPotential, SymmetryFunction};
    use rand::{RngExt, SeedableRng};

    /// G2 + G4 descriptor set with 3 radial and 3 angular functions.
    fn mixed_symmetry_set(cutoff: f64) -> SymmetryFunctionSet {
        let mut set = SymmetryFunctionSet::new();
        for &(eta, rs) in &[(0.3, 0.0), (0.8, 1.0), (1.5, 2.0)] {
            set.push(SymmetryFunction::G2 {
                eta,
                rs,
                rc: cutoff,
            });
        }
        for &(eta, zeta, lambda) in &[(0.2, 1.0, 1.0), (0.1, 2.0, -1.0), (0.05, 4.0, 1.0)] {
            set.push(SymmetryFunction::G4 {
                eta,
                zeta,
                lambda,
                rc: cutoff,
            });
        }
        set
    }

    fn mixed_potential(activation: Activation) -> FeedForwardPotential {
        let cutoff = 4.5;
        FeedForwardPotential::new(
            &[6, 10, 8, 1],
            activation,
            mixed_symmetry_set(cutoff),
            cutoff,
        )
    }

    /// Random cluster with a minimum pair distance (no overlapping atoms).
    fn random_cluster(n: usize, extent: f64, min_dist: f64, seed: u64) -> Vec<Vec3> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut atoms: Vec<Vec3> = Vec::new();
        while atoms.len() < n {
            let p = Vec3::new(
                rng.random_range(-extent..extent),
                rng.random_range(-extent..extent),
                rng.random_range(-extent..extent),
            );
            if atoms.iter().all(|q| (q - p).norm() > min_dist) {
                atoms.push(p);
            }
        }
        atoms
    }

    fn max_force_error(analytic: &[Vec3], reference: &[Vec3]) -> (f64, f64) {
        let mut max_diff = 0.0_f64;
        let mut max_force = 0.0_f64;
        for (fa, fr) in analytic.iter().zip(reference.iter()) {
            for a in 0..3 {
                max_diff = max_diff.max((fa[a] - fr[a]).abs());
                max_force = max_force.max(fr[a].abs());
            }
        }
        (max_diff, max_force)
    }

    #[test]
    fn analytic_forces_match_finite_differences_for_a_ten_atom_cluster() {
        let pot = mixed_potential(Activation::Tanh);
        let positions = random_cluster(10, 2.2, 0.9, 42);
        let species = vec![1; 10];
        let (energy, forces) = pot.energy_and_forces(&positions, &species);
        let reference = pot.finite_difference_forces(&positions, &species);
        let (max_diff, max_force) = max_force_error(&forces, &reference);
        assert!(max_force > 1e-3, "test cluster must feel real forces");
        assert!(
            max_diff / max_force < 1e-5,
            "relative force error {} (max |F| = {max_force})",
            max_diff / max_force
        );
        let e_ref = pot.total_energy(&positions, &species);
        assert!((energy - e_ref).abs() < 1e-12 * e_ref.abs().max(1.0));
    }

    #[test]
    fn analytic_forces_match_finite_differences_for_sigmoid_and_relu_networks() {
        for activation in [Activation::Sigmoid, Activation::ReLU, Activation::Linear] {
            let pot = mixed_potential(activation);
            let positions = random_cluster(8, 2.0, 0.9, 7);
            let species = vec![1; 8];
            let (_, forces) = pot.energy_and_forces(&positions, &species);
            let reference = pot.finite_difference_forces(&positions, &species);
            let (max_diff, max_force) = max_force_error(&forces, &reference);
            assert!(
                max_diff / max_force.max(1e-3) < 1e-4,
                "{activation:?}: error {max_diff} vs max |F| {max_force}"
            );
        }
    }

    #[test]
    fn forces_obey_newton_and_are_torque_free() {
        let pot = mixed_potential(Activation::Tanh);
        let positions = random_cluster(12, 2.5, 0.9, 5);
        let (_, forces) = pot.energy_and_forces(&positions, &[1; 12]);
        let total: Vec3 = forces.iter().sum();
        assert!(total.norm() < 1e-10, "net force = {total:?}");
        let torque: Vec3 = positions
            .iter()
            .zip(forces.iter())
            .map(|(r, f)| r.cross(f))
            .sum();
        assert!(torque.norm() < 1e-8, "net torque = {torque:?}");
    }

    #[test]
    fn cutoff_crossing_is_smooth_in_energy_and_force() {
        // Two atoms straddling the cutoff: the cosine taper makes both the energy and the
        // force vanish continuously there.
        let cutoff = 3.0;
        let pot = FeedForwardPotential::new(
            &[4, 6, 1],
            Activation::Tanh,
            FeedForwardPotential::default_symmetry_set(cutoff),
            cutoff,
        );
        let at = |d: f64| vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(d, 0.0, 0.0)];
        let (e_in, f_in) = pot.energy_and_forces(&at(cutoff - 1e-7), &[1, 1]);
        let (e_out, f_out) = pot.energy_and_forces(&at(cutoff + 1e-7), &[1, 1]);
        assert!((e_in - e_out).abs() < 1e-9, "energy jump {}", e_in - e_out);
        assert!(f_in[0].norm() < 1e-6 && f_out[0].norm() < 1e-12);
    }

    #[test]
    fn collinear_atoms_with_fractional_zeta_stay_finite() {
        // λ = -1 on a straight line puts 1 + λ cos θ at 0; rounding must not make it negative.
        let cutoff = 5.0;
        let mut set = SymmetryFunctionSet::new();
        set.push(SymmetryFunction::G4 {
            eta: 0.1,
            zeta: 1.5,
            lambda: -1.0,
            rc: cutoff,
        });
        set.push(SymmetryFunction::G4 {
            eta: 0.1,
            zeta: 1.5,
            lambda: 1.0,
            rc: cutoff,
        });
        let pot = FeedForwardPotential::new(&[2, 4, 1], Activation::Tanh, set, cutoff);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.9, 0.0, 0.0),
            Vec3::new(1.7, 0.0, 0.0),
            Vec3::new(2.9, 0.0, 0.0),
        ];
        let (energy, forces) = pot.energy_and_forces(&positions, &[1; 4]);
        assert!(energy.is_finite());
        assert!(forces.iter().all(|f| f.iter().all(|x| x.is_finite())));
    }

    #[test]
    fn forces_are_strictly_local() {
        // Moving a distant, non-interacting group of atoms leaves the forces elsewhere unchanged.
        let pot = mixed_potential(Activation::Tanh);
        let near = random_cluster(6, 1.5, 0.9, 3);
        let mut positions = near.clone();
        positions.extend(near.iter().map(|p| p + Vec3::new(50.0, 0.0, 0.0)));
        let species = vec![1; 12];
        let (_, forces_a) = pot.energy_and_forces(&positions, &species);
        for p in positions.iter_mut().skip(6) {
            *p += Vec3::new(3.0, 7.0, -11.0);
        }
        let (_, forces_b) = pot.energy_and_forces(&positions, &species);
        for i in 0..12 {
            assert!((forces_a[i] - forces_b[i]).norm() < 1e-12, "atom {i}");
        }
        for i in 0..6 {
            assert!(
                (forces_a[i] - forces_a[i + 6]).norm() < 1e-10,
                "copies agree"
            );
        }
    }

    #[test]
    fn a_thousand_atom_cluster_is_cheap_and_consistent() {
        // The previous 6N full-energy finite differences would need ~6000 O(N²) energy sweeps.
        let cutoff = 4.0;
        let pot = FeedForwardPotential::new(
            &[4, 8, 1],
            Activation::Tanh,
            FeedForwardPotential::default_symmetry_set(cutoff),
            cutoff,
        );
        let mut positions = Vec::new();
        for ix in 0..10 {
            for iy in 0..10 {
                for iz in 0..10 {
                    positions.push(Vec3::new(ix as f64, iy as f64, iz as f64) * 1.5);
                }
            }
        }
        let (energy, forces) = pot.energy_and_forces(&positions, &vec![1; positions.len()]);
        assert!(energy.is_finite());
        let total: Vec3 = forces.iter().sum();
        assert!(total.norm() < 1e-8, "net force = {total:?}");
        // The central atom of a symmetric lattice shell feels no net force.
        let centre = 5 * 100 + 5 * 10 + 5;
        assert!(
            forces[centre].norm() < 1e-9,
            "centre force {:?}",
            forces[centre]
        );
    }

    #[test]
    fn nve_dynamics_conserves_energy_with_analytic_forces() {
        use super::super::types::NnAtomisticSystem;
        let cutoff = 5.0;
        let pot = FeedForwardPotential::new(
            &[4, 6, 1],
            Activation::Tanh,
            FeedForwardPotential::default_symmetry_set(cutoff),
            cutoff,
        );
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.8, 0.1, 0.0),
            Vec3::new(0.7, 1.6, 0.2),
        ];
        let mut system = NnAtomisticSystem::new(positions, vec![12.0; 3], vec![1; 3]);
        system.velocities[0] = Vec3::new(0.02, 0.0, 0.01);
        system.velocities[2] = Vec3::new(-0.01, 0.015, 0.0);
        let e0 = system.total_energy(&pot);
        let mut max_drift = 0.0_f64;
        for _ in 0..400 {
            system.step(&pot, 0.02);
            max_drift = max_drift.max((system.total_energy(&pot) - e0).abs());
        }
        assert!(
            max_drift < 1e-5 * e0.abs().max(1.0),
            "energy drift {max_drift} (E0 = {e0})"
        );
    }

    #[test]
    fn g2_derivative_matches_numerical_derivative() {
        let (eta, rs, rc) = (0.7, 0.8, 4.0);
        for &r in &[0.5, 1.3, 2.9, 3.9] {
            let h = 1e-6;
            let value = |x: f64| (-eta * (x - rs).powi(2)).exp() * cutoff_function(x, rc);
            let numeric = (value(r + h) - value(r - h)) / (2.0 * h);
            let analytic = g2_radial_derivative(r, eta, rs, rc);
            assert!((numeric - analytic).abs() < 1e-8, "r = {r}");
        }
    }

    #[test]
    fn g4_gradient_matches_numerical_gradient_and_is_translation_invariant() {
        let (eta, zeta, lambda, rc) = (0.15, 2.0, -1.0, 4.0);
        let base = [
            Vec3::new(0.1, -0.2, 0.0),
            Vec3::new(1.4, 0.3, -0.5),
            Vec3::new(-0.6, 1.2, 0.8),
        ];
        let term = |p: &[Vec3; 3]| -> f64 {
            let nj = Neighbor {
                j: 1,
                dr: (p[1] - p[0]).into(),
                r: (p[1] - p[0]).norm(),
            };
            let nk = Neighbor {
                j: 2,
                dr: (p[2] - p[0]).into(),
                r: (p[2] - p[0]).norm(),
            };
            match Triplet::new(&nj, &nk, p, rc) {
                Some(t) => t.value(eta, zeta, lambda, rc),
                None => 0.0,
            }
        };
        let nj = Neighbor {
            j: 1,
            dr: (base[1] - base[0]).into(),
            r: (base[1] - base[0]).norm(),
        };
        let nk = Neighbor {
            j: 2,
            dr: (base[2] - base[0]).into(),
            r: (base[2] - base[0]).norm(),
        };
        let triplet = Triplet::new(&nj, &nk, &base, rc).expect("triplet inside the cutoff");
        let grads = triplet.gradient(eta, zeta, lambda, rc);
        let total = grads[0] + grads[1] + grads[2];
        assert!(total.norm() < 1e-14, "translation invariance: {total:?}");
        let h = 1e-6;
        for (atom, grad) in grads.iter().enumerate() {
            for axis in 0..3 {
                let mut plus = base;
                let mut minus = base;
                plus[atom][axis] += h;
                minus[atom][axis] -= h;
                let numeric = (term(&plus) - term(&minus)) / (2.0 * h);
                assert!(
                    (numeric - grad[axis]).abs() < 1e-8,
                    "atom {atom} axis {axis}: {numeric} vs {}",
                    grad[axis]
                );
            }
        }
    }
}
