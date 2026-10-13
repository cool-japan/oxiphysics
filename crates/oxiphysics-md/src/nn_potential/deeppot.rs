// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! DeePMD-SE (`se_e2_a`) descriptor, energy and analytic forces.
//!
//! For atom `i` with neighbors `j` (sorted by distance, nearest `N_c` kept):
//!
//! ```text
//! R̃_j = [s, s x̂, s ŷ, s ẑ]_ij                environment-matrix row
//! G_j  = embedding(s(r_ij))      ∈ R^{M1}
//! A    = (1/N_c) Σ_j G_j ⊗ R̃_j   ∈ R^{M1×4}
//! D    = A · A[:M2]ᵀ             ∈ R^{M1×M2}   (rotation/permutation invariant)
//! e_i  = fitting(D)
//! ```
//!
//! The same quantities are differentiated analytically for the forces:
//! `∂e/∂D` by back-propagation through the fitting network, then through
//! `D = A A[:M2]ᵀ`, `A`, the embedding network and the geometry of
//! `s(r) R̂`.

use super::backprop::NetworkTrace;
use super::locality::{Neighbor, neighbor_lists, neighbors_of_atom};
use super::types::{DeepPotDescriptor, DeepPotPotential};

/// Pairs closer than this have no defined direction and are ignored.
const MIN_DISTANCE: f64 = 1e-12;

impl DeepPotDescriptor {
    /// Keep the neighbors that enter the descriptor: inside the cutoff, nearest
    /// first (ties broken by atom index), at most `max_neighbours`.
    fn select_nearest(&self, mut neighbors: Vec<Neighbor>) -> Vec<Neighbor> {
        neighbors.retain(|nb| nb.r < self.cutoff && nb.r > MIN_DISTANCE);
        neighbors.sort_by(|a, b| a.r.total_cmp(&b.r).then(a.j.cmp(&b.j)));
        neighbors.truncate(self.max_neighbours);
        neighbors
    }

    /// Compute the DeePMD environment matrix for atom `center`.
    ///
    /// Returns a flattened row-major matrix of shape `[M × 4]` with
    /// `M = max_neighbours`:
    ///   column 0: s(r_ij)
    ///   columns 1-3: s(r_ij) * r_hat_ij[0..2]
    ///
    /// The rows are the neighbors inside the cutoff ordered by distance
    /// (nearest first), truncated to the nearest `max_neighbours` and padded
    /// with zero rows, so the matrix does not depend on the atom ordering.
    pub fn environment_matrix(&self, positions: &[[f64; 3]], center: usize) -> Vec<f64> {
        let neighbors = self.select_nearest(neighbors_of_atom(positions, center, self.cutoff));
        let mut matrix = Vec::with_capacity(4 * self.max_neighbours);
        for nb in &neighbors {
            let s = self.smooth_fn(nb.r);
            matrix.extend_from_slice(&[
                s,
                s * nb.dr[0] / nb.r,
                s * nb.dr[1] / nb.r,
                s * nb.dr[2] / nb.r,
            ]);
        }
        matrix.resize(4 * self.max_neighbours, 0.0);
        matrix
    }
}

/// Forward-pass quantities of one atom, kept for back-propagation.
struct AtomForward {
    neighbors: Vec<Neighbor>,
    /// `s(r_j)` per neighbor.
    switch: Vec<f64>,
    /// Unit vectors `r̂_ij` per neighbor.
    unit: Vec<[f64; 3]>,
    /// Embedding-network trace per neighbor (input `s(r_j)`, output `G_j`).
    embedding: Vec<NetworkTrace>,
    /// `A`, row-major `M1 × 4`.
    a: Vec<f64>,
    /// Fitting-network trace (input `D`, output the atomic energy).
    fitting: NetworkTrace,
}

impl AtomForward {
    fn energy(&self) -> f64 {
        self.fitting.output().iter().sum()
    }
}

impl DeepPotPotential {
    /// Embedding feature dimension `M1`.
    fn feature_dim(&self) -> usize {
        self.embedding.last().map_or(1, |layer| layer.biases.len())
    }

    /// Number of axis neurons actually used, `M2 ≤ M1`.
    fn axis_dim(&self) -> usize {
        self.axis_neurons.min(self.feature_dim())
    }

    /// Normalisation `N_c` of the neighbor sum (the fixed, padded neighbor count).
    fn norm_count(&self) -> f64 {
        self.descriptor.max_neighbours.max(1) as f64
    }

    /// Forward pass for one atom from its raw neighbor list.
    fn forward(&self, raw_neighbors: Vec<Neighbor>) -> AtomForward {
        let m1 = self.feature_dim();
        let m2 = self.axis_dim();
        let inv_nc = 1.0 / self.norm_count();
        let neighbors = self.descriptor.select_nearest(raw_neighbors);

        let mut a = vec![0.0; m1 * 4];
        let mut switch = Vec::with_capacity(neighbors.len());
        let mut unit = Vec::with_capacity(neighbors.len());
        let mut embedding = Vec::with_capacity(neighbors.len());
        for nb in &neighbors {
            let s = self.descriptor.smooth_fn(nb.r);
            let u = [nb.dr[0] / nb.r, nb.dr[1] / nb.r, nb.dr[2] / nb.r];
            let row = [s, s * u[0], s * u[1], s * u[2]];
            let trace = NetworkTrace::forward(&self.embedding, &[s]);
            for (m, &g) in trace.output().iter().enumerate().take(m1) {
                for (c, &r) in row.iter().enumerate() {
                    a[m * 4 + c] += g * r * inv_nc;
                }
            }
            switch.push(s);
            unit.push(u);
            embedding.push(trace);
        }

        // D = A · A[:M2]ᵀ, flattened row-major (M1 × M2).
        let mut d = vec![0.0; m1 * m2];
        for m in 0..m1 {
            for q in 0..m2 {
                d[m * m2 + q] = (0..4).map(|c| a[m * 4 + c] * a[q * 4 + c]).sum();
            }
        }
        let fitting = NetworkTrace::forward(&self.fitting, &d);
        AtomForward {
            neighbors,
            switch,
            unit,
            embedding,
            a,
            fitting,
        }
    }

    /// `∂e/∂d_j` (with `d_j = r_j - r_i`) for every selected neighbor.
    fn backward(&self, fwd: &AtomForward) -> Vec<[f64; 3]> {
        let m1 = self.feature_dim();
        let m2 = self.axis_dim();
        let inv_nc = 1.0 / self.norm_count();
        let a = &fwd.a;

        // Γ = ∂e/∂D from the fitting network.
        let seed = vec![1.0; fwd.fitting.output().len()];
        let gamma = fwd.fitting.backward(&self.fitting, &seed);

        // Ā[p][c] = ∂e/∂A[p][c]: D[m][q] = Σ_c A[m][c] A[q][c] enters A through
        // the first factor (row p = m) and, for p < M2, the second (q = p).
        let mut a_bar = vec![0.0; m1 * 4];
        for p in 0..m1 {
            for c in 0..4 {
                let mut v: f64 = (0..m2).map(|q| gamma[p * m2 + q] * a[q * 4 + c]).sum();
                if p < m2 {
                    v += (0..m1)
                        .map(|m| gamma[m * m2 + p] * a[m * 4 + c])
                        .sum::<f64>();
                }
                a_bar[p * 4 + c] = v;
            }
        }

        let mut grads = Vec::with_capacity(fwd.neighbors.len());
        for (n, nb) in fwd.neighbors.iter().enumerate() {
            let s = fwd.switch[n];
            let u = fwd.unit[n];
            let row = [s, s * u[0], s * u[1], s * u[2]];
            let features = fwd.embedding[n].output();

            // ∂e/∂G_j[m] and ∂e/∂R̃_j[c] from A = (1/N_c) Σ_j G_j ⊗ R̃_j.
            let d_features: Vec<f64> = (0..m1)
                .map(|m| (0..4).map(|c| a_bar[m * 4 + c] * row[c]).sum::<f64>() * inv_nc)
                .collect();
            let mut d_row = [0.0; 4];
            for (c, d_c) in d_row.iter_mut().enumerate() {
                *d_c = (0..m1)
                    .map(|m| a_bar[m * 4 + c] * features.get(m).copied().unwrap_or(0.0))
                    .sum::<f64>()
                    * inv_nc;
            }

            // ∂e/∂s: through the embedding input, and explicitly through R̃.
            let through_embedding = fwd.embedding[n]
                .backward(&self.embedding, &d_features)
                .first()
                .copied()
                .unwrap_or(0.0);
            let d_switch =
                through_embedding + d_row[0] + d_row[1] * u[0] + d_row[2] * u[1] + d_row[3] * u[2];
            // ∂e/∂r̂ (explicit): R̃[1+α] = s r̂_α.
            let d_unit = [d_row[1] * s, d_row[2] * s, d_row[3] * s];
            let along = d_unit[0] * u[0] + d_unit[1] * u[1] + d_unit[2] * u[2];
            let slope = self.descriptor.smooth_fn_derivative(nb.r);

            // ∂r/∂d = r̂ and ∂r̂/∂d = (I - r̂ r̂ᵀ)/r.
            let mut grad = [0.0; 3];
            for (alpha, g) in grad.iter_mut().enumerate() {
                *g = d_switch * slope * u[alpha] + (d_unit[alpha] - along * u[alpha]) / nb.r;
            }
            grads.push(grad);
        }
        grads
    }

    /// Atomic energy of atom `center` (neighbors inside the cutoff, nearest
    /// `max_neighbours`, enter through the invariant DeePMD-SE descriptor).
    pub fn atomic_energy(&self, positions: &[[f64; 3]], center: usize) -> f64 {
        let raw = neighbors_of_atom(positions, center, self.descriptor.cutoff);
        self.forward(raw).energy()
    }

    /// Total energy = sum of atomic energies.
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let lists = neighbor_lists(positions, self.descriptor.cutoff);
        lists
            .into_iter()
            .map(|list| self.forward(list).energy())
            .sum()
    }

    /// Total energy and analytic forces `F = -∂E/∂r`.
    ///
    /// One pass over cutoff neighbor lists, `O(N · N_c · M1)` for the
    /// descriptors and back-propagation; the previous `6N` full-system energy
    /// evaluations are only available as
    /// [`DeepPotPotential::forces_finite_difference`].
    pub fn energy_and_forces(&self, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let lists = neighbor_lists(positions, self.descriptor.cutoff);
        let mut energy = 0.0;
        let mut forces = vec![[0.0_f64; 3]; positions.len()];
        for (i, list) in lists.into_iter().enumerate() {
            let fwd = self.forward(list);
            energy += fwd.energy();
            for (nb, grad) in fwd.neighbors.iter().zip(self.backward(&fwd)) {
                for (alpha, g) in grad.iter().enumerate() {
                    forces[nb.j][alpha] -= g;
                    forces[i][alpha] += g;
                }
            }
        }
        (energy, forces)
    }

    /// Analytic forces on every atom (see [`DeepPotPotential::energy_and_forces`]).
    pub fn forces(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.energy_and_forces(positions).1
    }

    /// Forces from central finite differences of the total energy with step `h`.
    ///
    /// `6N` full-system energy evaluations: an independent reference for
    /// validating [`DeepPotPotential::forces`], not for production dynamics.
    pub fn forces_finite_difference(&self, positions: &[[f64; 3]], h: f64) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut pos = positions.to_vec();
        for i in 0..n {
            for a in 0..3 {
                pos[i][a] += h;
                let e_plus = self.total_energy(&pos);
                pos[i][a] -= 2.0 * h;
                let e_minus = self.total_energy(&pos);
                pos[i][a] += h;
                forces[i][a] = -(e_plus - e_minus) / (2.0 * h);
            }
        }
        forces
    }
}

#[cfg(test)]
mod tests {
    use super::super::types::Activation;
    use super::*;
    use rand::{RngExt, SeedableRng};

    fn potential() -> DeepPotPotential {
        DeepPotPotential::new(6, &[8, 12], &[16, 16, 1], 5.0, 2.0)
    }

    /// Random cluster with a minimum pair distance.
    fn cluster(n: usize, extent: f64, min_dist: f64, seed: u64) -> Vec<[f64; 3]> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut atoms: Vec<[f64; 3]> = Vec::new();
        while atoms.len() < n {
            let p = [
                rng.random_range(-extent..extent),
                rng.random_range(-extent..extent),
                rng.random_range(-extent..extent),
            ];
            let far_enough = atoms.iter().all(|q| {
                let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() > min_dist
            });
            if far_enough {
                atoms.push(p);
            }
        }
        atoms
    }

    /// Random proper rotation from a unit quaternion.
    fn random_rotation(seed: u64) -> [[f64; 3]; 3] {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut q = [0.0_f64; 4];
        for x in q.iter_mut() {
            *x = crate::gaussian_noise::standard_normal(&mut rng);
        }
        let norm = q.iter().map(|x| x * x).sum::<f64>().sqrt();
        let (w, x, y, z) = (q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm);
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - z * w),
                2.0 * (x * z + y * w),
            ],
            [
                2.0 * (x * y + z * w),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - x * w),
            ],
            [
                2.0 * (x * z - y * w),
                2.0 * (y * z + x * w),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    }

    fn rotate(rot: &[[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
        [
            rot[0][0] * p[0] + rot[0][1] * p[1] + rot[0][2] * p[2],
            rot[1][0] * p[0] + rot[1][1] * p[1] + rot[1][2] * p[2],
            rot[2][0] * p[0] + rot[2][1] * p[1] + rot[2][2] * p[2],
        ]
    }

    fn rel_close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
    }

    #[test]
    fn energy_is_invariant_under_rotation() {
        // Regression: the raw environment matrix was fed into an MLP, so rotating the
        // structure changed the energy.
        let pot = potential();
        let positions = cluster(7, 1.8, 0.8, 11);
        let e0 = pot.total_energy(&positions);
        for seed in 0..5 {
            let rot = random_rotation(100 + seed);
            let rotated: Vec<[f64; 3]> = positions.iter().map(|&p| rotate(&rot, p)).collect();
            let e1 = pot.total_energy(&rotated);
            assert!(rel_close(e0, e1, 1e-10), "rotation {seed}: {e0} vs {e1}");
        }
        // Improper rotation (inversion) as well.
        let inverted: Vec<[f64; 3]> = positions.iter().map(|p| [-p[0], -p[1], -p[2]]).collect();
        assert!(rel_close(e0, pot.total_energy(&inverted), 1e-10));
    }

    #[test]
    fn energy_is_invariant_under_translation() {
        let pot = potential();
        let positions = cluster(7, 1.8, 0.8, 12);
        let e0 = pot.total_energy(&positions);
        let shifted: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| [p[0] + 13.7, p[1] - 4.2, p[2] + 0.31])
            .collect();
        assert!(rel_close(e0, pot.total_energy(&shifted), 1e-10));
    }

    #[test]
    fn energy_is_invariant_under_atom_permutation() {
        // Regression: neighbors were truncated and ordered by atom index.
        let pot = potential();
        let positions = cluster(7, 1.8, 0.8, 13);
        let e0 = pot.total_energy(&positions);
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        for _ in 0..5 {
            let mut order: Vec<usize> = (0..positions.len()).collect();
            for k in (1..order.len()).rev() {
                let j = rng.random_range(0..=k);
                order.swap(k, j);
            }
            let permuted: Vec<[f64; 3]> = order.iter().map(|&k| positions[k]).collect();
            assert!(rel_close(e0, pot.total_energy(&permuted), 1e-10));
            // Atomic energies follow the relabelling.
            for (new_idx, &old_idx) in order.iter().enumerate() {
                let a = pot.atomic_energy(&positions, old_idx);
                let b = pot.atomic_energy(&permuted, new_idx);
                assert!(rel_close(a, b, 1e-10), "atom {old_idx}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn truncation_to_nearest_neighbors_keeps_the_invariances() {
        // 3 neighbors kept out of 8 atoms within the cutoff.
        let pot = DeepPotPotential::new(3, &[6, 10], &[12, 1], 6.0, 2.5);
        let positions = cluster(9, 1.6, 0.7, 21);
        let e0 = pot.total_energy(&positions);
        let rot = random_rotation(9);
        let rotated: Vec<[f64; 3]> = positions.iter().rev().map(|&p| rotate(&rot, p)).collect();
        assert!(rel_close(e0, pot.total_energy(&rotated), 1e-10));
    }

    #[test]
    fn environment_matrix_is_sorted_by_distance_and_padded() {
        let pot = DeepPotPotential::new(5, &[4], &[4, 1], 6.0, 2.0);
        let d = &pot.descriptor;
        // Atoms deliberately listed far-to-near.
        let positions = [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 1.0],
            [9.0, 0.0, 0.0], // outside the cutoff
        ];
        let m = d.environment_matrix(&positions, 0);
        assert_eq!(m.len(), 20);
        // rows: r = 1, 2, 3 then two zero rows
        assert!((m[0] - 1.0).abs() < 1e-12 && (m[3] - 1.0).abs() < 1e-12); // along +z at r = 1
        assert!((m[4] - d.smooth_fn(2.0)).abs() < 1e-12 && m[5].abs() < 1e-12);
        assert!((m[8] - d.smooth_fn(3.0)).abs() < 1e-12);
        assert!(m[12..].iter().all(|&x| x == 0.0), "padding rows are zero");
        // Reordering the atoms does not change the matrix.
        let reordered = [
            positions[0],
            positions[4],
            positions[3],
            positions[1],
            positions[2],
        ];
        assert_eq!(d.environment_matrix(&reordered, 0), m);
    }

    #[test]
    fn environment_matrix_keeps_only_the_nearest_max_neighbours() {
        let pot = DeepPotPotential::new(2, &[4], &[4, 1], 6.0, 2.0);
        let positions = [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ];
        let m = pot.descriptor.environment_matrix(&positions, 0);
        assert_eq!(m.len(), 8);
        assert!((m[0] - 1.0).abs() < 1e-12);
        assert!((m[4] - pot.descriptor.smooth_fn(2.0)).abs() < 1e-12);
    }

    #[test]
    fn smooth_fn_derivative_matches_numerical_derivative() {
        let d = DeepPotDescriptor::new(4, 5.0, 2.0);
        for &r in &[0.7, 1.5, 2.0001, 2.8, 3.9, 4.9] {
            let h = 1e-6;
            let numeric = (d.smooth_fn(r + h) - d.smooth_fn(r - h)) / (2.0 * h);
            assert!(
                (numeric - d.smooth_fn_derivative(r)).abs() < 1e-7,
                "r = {r}: {numeric} vs {}",
                d.smooth_fn_derivative(r)
            );
        }
        assert_eq!(d.smooth_fn_derivative(5.0), 0.0);
        assert_eq!(d.smooth_fn_derivative(7.0), 0.0);
        assert_eq!(d.smooth_fn_derivative(0.0), 0.0);
    }

    #[test]
    fn analytic_forces_match_finite_differences() {
        let pot = potential();
        let positions = cluster(7, 1.8, 0.8, 31);
        let (energy, forces) = pot.energy_and_forces(&positions);
        assert!((energy - pot.total_energy(&positions)).abs() < 1e-12);
        let reference = pot.forces_finite_difference(&positions, 1e-5);
        let max_force = reference
            .iter()
            .flat_map(|f| f.iter().map(|x| x.abs()))
            .fold(0.0_f64, f64::max);
        assert!(max_force > 1e-4, "need non-trivial forces, got {max_force}");
        for (fa, fr) in forces.iter().zip(reference.iter()) {
            for a in 0..3 {
                assert!(
                    (fa[a] - fr[a]).abs() < 1e-6 * max_force,
                    "analytic {} vs numeric {}",
                    fa[a],
                    fr[a]
                );
            }
        }
    }

    #[test]
    fn analytic_forces_match_finite_differences_with_truncated_neighbors() {
        let pot = DeepPotPotential::new(3, &[6, 10], &[12, 12, 1], 6.0, 2.5);
        let positions = cluster(8, 1.7, 0.8, 41);
        let forces = pot.forces(&positions);
        let reference = pot.forces_finite_difference(&positions, 1e-5);
        let max_force = reference
            .iter()
            .flat_map(|f| f.iter().map(|x| x.abs()))
            .fold(0.0_f64, f64::max);
        for (fa, fr) in forces.iter().zip(reference.iter()) {
            for a in 0..3 {
                assert!((fa[a] - fr[a]).abs() < 1e-6 * max_force.max(1e-3));
            }
        }
    }

    #[test]
    fn forces_obey_newton_and_are_torque_free() {
        let pot = potential();
        let positions = cluster(8, 1.8, 0.8, 51);
        let forces = pot.forces(&positions);
        for axis in 0..3 {
            let total: f64 = forces.iter().map(|f| f[axis]).sum();
            assert!(total.abs() < 1e-10, "net force along {axis} = {total}");
        }
        let mut torque = [0.0_f64; 3];
        for (r, f) in positions.iter().zip(forces.iter()) {
            torque[0] += r[1] * f[2] - r[2] * f[1];
            torque[1] += r[2] * f[0] - r[0] * f[2];
            torque[2] += r[0] * f[1] - r[1] * f[0];
        }
        for (axis, t) in torque.iter().enumerate() {
            assert!(t.abs() < 1e-8, "net torque along {axis} = {t}");
        }
    }

    #[test]
    fn energy_and_force_are_continuous_at_the_cutoff() {
        let pot = DeepPotPotential::new(4, &[6, 8], &[10, 1], 4.0, 2.0);
        let at = |d: f64| vec![[0.0, 0.0, 0.0], [d, 0.0, 0.0]];
        let inside = at(4.0 - 1e-7);
        let outside = at(4.0 + 1e-7);
        assert!((pot.total_energy(&inside) - pot.total_energy(&outside)).abs() < 1e-8);
        let f_in = pot.forces(&inside);
        assert!(
            f_in[0][0].abs() < 1e-5,
            "force at the cutoff = {}",
            f_in[0][0]
        );
    }

    #[test]
    fn isolated_atom_has_a_constant_finite_energy() {
        let mut pot = potential();
        // Non-zero biases make the isolated-atom energy a non-trivial constant.
        pot.fitting[0].biases[0] = 0.3;
        if let Some(last) = pot.fitting.last_mut() {
            last.biases[0] = 0.7;
        }
        let e1 = pot.total_energy(&[[0.0; 3]]);
        let e2 = pot.total_energy(&[[0.0; 3], [50.0, 0.0, 0.0]]);
        assert!(e1.is_finite() && e1.abs() > 1e-3, "e1 = {e1}");
        assert!(
            (e2 - 2.0 * e1).abs() < 1e-12,
            "two isolated atoms = 2 x single"
        );
        assert_eq!(pot.forces(&[[0.0; 3]])[0], [0.0; 3]);
    }

    #[test]
    fn zero_neighbor_capacity_is_well_defined() {
        let pot = DeepPotPotential::new(0, &[4], &[4, 1], 5.0, 2.0);
        let positions = [[0.0; 3], [1.0, 0.0, 0.0]];
        let e = pot.total_energy(&positions);
        assert!(e.is_finite());
        assert!(pot.forces(&positions).iter().all(|f| *f == [0.0; 3]));
    }

    #[test]
    fn constructor_builds_the_se_architecture() {
        let pot = DeepPotPotential::new(10, &[8, 16], &[32, 32, 1], 6.0, 3.0);
        // The embedding network maps the scalar s(r) to M1 = 16 features.
        assert_eq!(pot.embedding[0].weights[0].len(), 1);
        assert_eq!(pot.embedding.last().map(|l| l.biases.len()), Some(16));
        // M2 = min(M1, 4) axis neurons -> the fitting input is M1 * M2.
        assert_eq!(pot.axis_neurons, 4);
        assert_eq!(pot.fitting[0].weights[0].len(), 16 * 4);
        // The energy output is linear, hidden layers are tanh.
        assert!(matches!(
            pot.fitting.last().map(|l| l.activation),
            Some(Activation::Linear)
        ));
        assert!(matches!(pot.fitting[0].activation, Activation::Tanh));
        // Explicit axis neurons are clamped to 1..=M1.
        let wide = DeepPotPotential::new_with_axis_neurons(6, &[8], 99, &[8, 1], 5.0, 2.0);
        assert_eq!(wide.axis_neurons, 8);
        let narrow = DeepPotPotential::new_with_axis_neurons(6, &[8], 0, &[8, 1], 5.0, 2.0);
        assert_eq!(narrow.axis_neurons, 1);
        assert_eq!(narrow.fitting[0].weights[0].len(), 8);
    }

    #[test]
    fn atomic_energy_is_not_confined_to_the_tanh_range() {
        // With a linear output layer the energy is unbounded; the previous all-tanh network
        // clamped every atomic energy into (-1, 1).
        let mut pot = DeepPotPotential::new(4, &[4], &[4, 1], 5.0, 2.0);
        if let Some(last) = pot.fitting.last_mut() {
            last.biases[0] = 25.0;
        }
        let e = pot.atomic_energy(&[[0.0; 3], [1.0, 0.0, 0.0]], 0);
        assert!(e > 20.0, "e = {e}");
    }
}
