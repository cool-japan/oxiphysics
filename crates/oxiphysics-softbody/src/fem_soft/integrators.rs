// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM integrators, boundary conditions, and assembled soft body types.
//!
//! Contains LumpedMassMatrix, RayleighDamping, VelocityVerletIntegrator,
//! BoundaryConditions, CorotFemElement4, FemSoftBodyVerlet, and
//! assemble_internal_forces.

use super::corot_sim::{CorotFemNode, CorotFemTet};
use super::math_helpers::{det3x3, edge_matrix_raw, inv3x3, mul3x3, transpose3x3};

// ---------------------------------------------------------------------------
// Lumped mass matrix
// ---------------------------------------------------------------------------

/// Lumped (diagonal) mass matrix for a FEM mesh of tetrahedral elements.
///
/// Each element distributes its mass equally to its four nodes (1/4 each).
/// This avoids inverting a consistent mass matrix and is widely used in
/// explicit FEM.
#[derive(Debug, Clone)]
pub struct LumpedMassMatrix {
    /// Per-node lumped mass (kg).  Index matches the node array.
    pub masses: Vec<f64>,
}

impl LumpedMassMatrix {
    /// Build a lumped mass matrix from element rest volumes and material density.
    ///
    /// `n_nodes`  -- total number of nodes in the mesh.
    /// `elements` -- slice of `(node_indices, rest_volume)` pairs.
    /// `density`  -- mass density (kg/m^3).
    pub fn from_elements(n_nodes: usize, elements: &[([usize; 4], f64)], density: f64) -> Self {
        let mut masses = vec![0.0_f64; n_nodes];
        for (indices, vol) in elements {
            let m_elem = density * vol;
            let m_per_node = m_elem * 0.25;
            for &idx in indices {
                masses[idx] += m_per_node;
            }
        }
        Self { masses }
    }

    /// Return the inverse mass for node `i` (0 if mass is zero or tiny).
    #[inline]
    pub fn inv_mass(&self, i: usize) -> f64 {
        let m = self.masses[i];
        if m > 1e-30 { 1.0 / m } else { 0.0 }
    }

    /// Apply the inverse mass matrix: return `M^{-1} * f` (element-wise).
    pub fn apply_inv(&self, forces: &[[f64; 3]]) -> Vec<[f64; 3]> {
        forces
            .iter()
            .enumerate()
            .map(|(i, &f)| {
                let inv_m = self.inv_mass(i);
                [f[0] * inv_m, f[1] * inv_m, f[2] * inv_m]
            })
            .collect()
    }

    /// Total mass of the mesh.
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// Rayleigh damping (standalone, f64-array)
// ---------------------------------------------------------------------------

/// Rayleigh damping: `C = alpha * M + beta * K`.
///
/// The damping force on node `i` is `f_d,i = -alpha * m_i * v_i - beta_K * (K v)_i`.
///
/// The stiffness-proportional term needs the product of the tangent stiffness
/// `K = -d f_el / d x` with the nodal velocities.  It is a function of the
/// **elastic** forces only; it must never be derived from an array that already
/// holds gravity or other external loads (a velocity-independent force is not
/// damping).  [`stiffness_velocity_product`](Self::stiffness_velocity_product)
/// evaluates `K v` matrix-free from an elastic force callback, and
/// [`apply`](Self::apply) adds both damping contributions to the force array.
#[derive(Debug, Clone, Copy)]
pub struct RayleighDamping {
    /// Mass-proportional coefficient alpha (1/s).
    pub alpha: f64,
    /// Stiffness-proportional coefficient beta (s).
    pub beta_k: f64,
}

impl RayleighDamping {
    /// Create new Rayleigh damping coefficients.
    pub fn new(alpha: f64, beta_k: f64) -> Self {
        Self { alpha, beta_k }
    }

    /// Compute the mass-proportional damping force for one node.
    ///
    /// `f_mass = -alpha * mass * velocity`
    #[inline]
    pub fn mass_damping_force(&self, mass: f64, velocity: [f64; 3]) -> [f64; 3] {
        let c = -self.alpha * mass;
        [c * velocity[0], c * velocity[1], c * velocity[2]]
    }

    /// Compute the stiffness-proportional damping force for one node.
    ///
    /// `kv` is the node's entry of the stiffness-velocity product `K v` (see
    /// [`stiffness_velocity_product`](Self::stiffness_velocity_product)); the force
    /// is `f_stiff = -beta_K * (K v)`.
    #[inline]
    pub fn stiffness_damping_force(&self, kv: [f64; 3]) -> [f64; 3] {
        let c = -self.beta_k;
        [c * kv[0], c * kv[1], c * kv[2]]
    }

    /// `true` when the stiffness-proportional term is active (`beta_K != 0`), i.e.
    /// when callers have to supply a real `K v` to [`apply`](Self::apply).
    #[inline]
    pub fn uses_stiffness_term(&self) -> bool {
        self.beta_k != 0.0
    }

    /// Matrix-free stiffness-velocity product `K v` with `K = -d f_el / d x`.
    ///
    /// Evaluates the directional derivative of the elastic force field along the
    /// nodal velocities by a central difference,
    /// `K v ~ -(f_el(x + eps v) - f_el(x - eps v)) / (2 eps)`, which costs two
    /// evaluations of `elastic_forces`.
    ///
    /// * `positions` -- current nodal positions `x`.
    /// * `velocities` -- nodal velocities `v` (pass zero for fixed nodes).
    /// * `elastic_forces` -- returns the **restoring elastic** nodal forces of all
    ///   nodes for the given positions: no gravity, no external loads, no damping.
    ///
    /// `eps` is chosen so that the largest probe displacement is `1e-6` times the
    /// extent of the node cloud.  The result is all zeros when every velocity is
    /// zero.
    pub fn stiffness_velocity_product<F>(
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        mut elastic_forces: F,
    ) -> Vec<[f64; 3]>
    where
        F: FnMut(&[[f64; 3]]) -> Vec<[f64; 3]>,
    {
        let mut kv = vec![[0.0_f64; 3]; positions.len()];
        let v_max = velocities
            .iter()
            .flatten()
            .fold(0.0_f64, |acc, c| acc.max(c.abs()));
        if positions.is_empty() || !v_max.is_finite() || v_max <= 0.0 {
            return kv;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in positions {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        let extent = (0..3).map(|d| hi[d] - lo[d]).fold(0.0_f64, f64::max);
        let length_scale = if extent > 0.0 { extent } else { 1.0 };
        let eps = 1.0e-6 * length_scale / v_max;

        let probe = |sign: f64| -> Vec<[f64; 3]> {
            let mut shifted = positions.to_vec();
            for (x, v) in shifted.iter_mut().zip(velocities.iter()) {
                for d in 0..3 {
                    x[d] += sign * eps * v[d];
                }
            }
            shifted
        };
        let f_plus = elastic_forces(&probe(1.0));
        let f_minus = elastic_forces(&probe(-1.0));
        let scale = -1.0 / (2.0 * eps);
        for ((out, fp), fm) in kv.iter_mut().zip(f_plus.iter()).zip(f_minus.iter()) {
            for d in 0..3 {
                out[d] = scale * (fp[d] - fm[d]);
            }
        }
        kv
    }

    /// Apply Rayleigh damping to the full force array.
    ///
    /// Adds `-alpha * m_i * v_i - beta_K * kv_i` to `forces[i]`.  `kv` is the
    /// stiffness-velocity product `K v` computed from the elastic forces only (see
    /// [`stiffness_velocity_product`](Self::stiffness_velocity_product)); it may be
    /// all zeros when [`uses_stiffness_term`](Self::uses_stiffness_term) is `false`.
    pub fn apply(
        &self,
        forces: &mut [[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        kv: &[[f64; 3]],
    ) {
        debug_assert_eq!(forces.len(), velocities.len());
        debug_assert_eq!(forces.len(), masses.len());
        debug_assert_eq!(forces.len(), kv.len());
        for (((force, velocity), mass), kv_i) in forces
            .iter_mut()
            .zip(velocities.iter())
            .zip(masses.iter())
            .zip(kv.iter())
        {
            let f_stiff = self.stiffness_damping_force(*kv_i);
            let f_mass = self.mass_damping_force(*mass, *velocity);
            for d in 0..3 {
                force[d] += f_stiff[d] + f_mass[d];
            }
        }
    }

    /// Critical damping ratio for a given natural frequency omega (rad/s).
    ///
    /// xi = (alpha/(2*omega)) + (beta_K * omega / 2)
    pub fn damping_ratio(&self, omega: f64) -> f64 {
        if omega.abs() < 1e-30 {
            return 0.0;
        }
        self.alpha / (2.0 * omega) + self.beta_k * omega / 2.0
    }

    /// Compute alpha and beta that achieve critical damping ratios xi1, xi2 at two
    /// reference frequencies omega1, omega2 (rad/s).
    ///
    /// Solves: \[ 1/(2*omega1)  omega1/2 \] \[alpha\]   \[xi1\]
    ///          \[ 1/(2*omega2)  omega2/2 \] \[beta \] = \[xi2\]
    pub fn from_two_modes(xi1: f64, omega1: f64, xi2: f64, omega2: f64) -> Self {
        // det = (omega2 - omega1) / (2 * omega1 * omega2)
        let denom = omega2 * omega2 - omega1 * omega1;
        if denom.abs() < 1e-30 {
            return Self::new(0.0, 0.0);
        }
        let alpha = 2.0 * omega1 * omega2 * (xi1 * omega2 - xi2 * omega1) / denom;
        let beta_k = 2.0 * (xi2 * omega2 - xi1 * omega1) / denom;
        Self { alpha, beta_k }
    }
}

// ---------------------------------------------------------------------------
// Velocity Verlet integrator (standalone)
// ---------------------------------------------------------------------------

/// Velocity Verlet time integrator for FEM nodes.
///
/// The classic "leap-frog" Velocity Verlet scheme:
/// ```text
/// x_{n+1} = x_n + dt * v_n + 0.5 * dt^2 * a_n
/// a_{n+1} = M^{-1} * f(x_{n+1})
/// v_{n+1} = v_n + 0.5 * dt * (a_n + a_{n+1})
/// ```
/// This is a second-order, time-reversible, symplectic integrator.
#[derive(Debug, Clone)]
pub struct VelocityVerletIntegrator {
    /// Acceleration from the previous step.
    pub prev_accel: Vec<[f64; 3]>,
    /// Time step (s).
    pub dt: f64,
}

impl VelocityVerletIntegrator {
    /// Create a new integrator for `n_nodes` nodes with time step `dt`.
    pub fn new(n_nodes: usize, dt: f64) -> Self {
        Self {
            prev_accel: vec![[0.0; 3]; n_nodes],
            dt,
        }
    }

    /// Predictor step: update positions using previous acceleration.
    ///
    /// `x_{n+1} = x_n + dt * v + 0.5 * dt^2 * a_prev`
    pub fn predict_positions(
        &self,
        positions: &mut [[f64; 3]],
        velocities: &[[f64; 3]],
        pinned: &[bool],
    ) {
        let dt = self.dt;
        for i in 0..positions.len() {
            if pinned[i] {
                continue;
            }
            let a = self.prev_accel[i];
            let v = velocities[i];
            for d in 0..3 {
                positions[i][d] += dt * v[d] + 0.5 * dt * dt * a[d];
            }
        }
    }

    /// Corrector step: update velocities using average of old and new acceleration.
    ///
    /// `v_{n+1} = v_n + 0.5 * dt * (a_prev + a_new)`
    pub fn correct_velocities(
        &mut self,
        velocities: &mut [[f64; 3]],
        new_accel: &[[f64; 3]],
        pinned: &[bool],
    ) {
        let dt = self.dt;
        for i in 0..velocities.len() {
            if pinned[i] {
                continue;
            }
            let a_old = self.prev_accel[i];
            let a_new = new_accel[i];
            for d in 0..3 {
                velocities[i][d] += 0.5 * dt * (a_old[d] + a_new[d]);
            }
        }
        self.prev_accel.copy_from_slice(new_accel);
    }

    /// Full Velocity Verlet step for a `CorotFemSim`-like simulation.
    ///
    /// Applies elastic forces computed from the tet elements and optional
    /// Rayleigh damping.
    pub fn step_sim(
        &mut self,
        positions: &mut [[f64; 3]],
        velocities: &mut [[f64; 3]],
        masses: &[f64],
        tets: &[CorotFemTet],
        mu: f64,
        lambda: f64,
        gravity: [f64; 3],
        rayleigh: Option<&RayleighDamping>,
        pinned: &[bool],
    ) {
        let n = positions.len();

        // 1. Predict positions
        for i in 0..n {
            if pinned[i] {
                continue;
            }
            let a = self.prev_accel[i];
            let v = velocities[i];
            let dt = self.dt;
            for d in 0..3 {
                positions[i][d] += dt * v[d] + 0.5 * dt * dt * a[d];
            }
        }

        // 2. Compute forces at new positions
        let mut forces = vec![[0.0_f64; 3]; n];

        // Gravity
        for i in 0..n {
            if pinned[i] {
                continue;
            }
            for d in 0..3 {
                forces[i][d] += masses[i] * gravity[d];
            }
        }

        // Elastic (restoring) forces from the tets.
        let elastic = tet_elastic_forces(positions, masses, tets, mu, lambda);
        for (force, f_el) in forces.iter_mut().zip(elastic.iter()) {
            for d in 0..3 {
                force[d] += f_el[d];
            }
        }

        // Rayleigh damping.  The stiffness term uses K v of the elastic forces
        // alone, never the gravity-inclusive force array.
        if let Some(rd) = rayleigh {
            let kv = if rd.uses_stiffness_term() {
                let moving: Vec<[f64; 3]> = velocities
                    .iter()
                    .zip(pinned.iter())
                    .map(|(v, &fixed)| if fixed { [0.0; 3] } else { *v })
                    .collect();
                RayleighDamping::stiffness_velocity_product(positions, &moving, |x| {
                    tet_elastic_forces(x, masses, tets, mu, lambda)
                })
            } else {
                vec![[0.0_f64; 3]; n]
            };
            rd.apply(&mut forces, velocities, masses, &kv);
        }

        // 3. Compute new accelerations
        let new_accel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                if pinned[i] || masses[i] < 1e-30 {
                    return [0.0; 3];
                }
                let inv_m = 1.0 / masses[i];
                [
                    forces[i][0] * inv_m,
                    forces[i][1] * inv_m,
                    forces[i][2] * inv_m,
                ]
            })
            .collect();

        // 4. Correct velocities
        let dt = self.dt;
        for i in 0..n {
            if pinned[i] {
                continue;
            }
            let a_old = self.prev_accel[i];
            let a_new = new_accel[i];
            for d in 0..3 {
                velocities[i][d] += 0.5 * dt * (a_old[d] + a_new[d]);
            }
        }

        self.prev_accel = new_accel;
    }
}

/// Restoring elastic forces of all `tets` at the given nodal `positions`.
fn tet_elastic_forces(
    positions: &[[f64; 3]],
    masses: &[f64],
    tets: &[CorotFemTet],
    mu: f64,
    lambda: f64,
) -> Vec<[f64; 3]> {
    let nodes: Vec<CorotFemNode> = positions
        .iter()
        .zip(masses.iter())
        .map(|(position, mass)| CorotFemNode {
            position: *position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: *mass,
        })
        .collect();
    let mut forces = vec![[0.0_f64; 3]; positions.len()];
    for tet in tets {
        let fs = tet.compute_elastic_forces(&nodes, mu, lambda);
        for (k, fs_k) in fs.iter().enumerate() {
            let idx = tet.node_indices[k];
            for (force_d, fs_d) in forces[idx].iter_mut().zip(fs_k.iter()) {
                *force_d += fs_d;
            }
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// Fixed boundary condition manager
// ---------------------------------------------------------------------------

/// Manager for Dirichlet (fixed / prescribed) boundary conditions.
///
/// Stores a list of node indices that are held at fixed positions.
/// During each time step, these nodes are reset to their prescribed positions
/// and their velocities are zeroed.
#[derive(Debug, Clone)]
pub struct BoundaryConditions {
    /// Fixed nodes: (node_index, fixed_position).
    fixed: Vec<(usize, [f64; 3])>,
}

impl BoundaryConditions {
    /// Create a new (empty) boundary condition set.
    pub fn new() -> Self {
        Self { fixed: Vec::new() }
    }

    /// Pin node `idx` at its current position.
    pub fn pin_at(&mut self, idx: usize, position: [f64; 3]) {
        // Replace if already pinned
        for entry in &mut self.fixed {
            if entry.0 == idx {
                entry.1 = position;
                return;
            }
        }
        self.fixed.push((idx, position));
    }

    /// Unpin a node, removing its constraint.
    pub fn unpin(&mut self, idx: usize) {
        self.fixed.retain(|(i, _)| *i != idx);
    }

    /// Return true if node `idx` is pinned.
    pub fn is_pinned(&self, idx: usize) -> bool {
        self.fixed.iter().any(|(i, _)| *i == idx)
    }

    /// Number of pinned nodes.
    pub fn len(&self) -> usize {
        self.fixed.len()
    }

    /// True if no nodes are pinned.
    pub fn is_empty(&self) -> bool {
        self.fixed.is_empty()
    }

    /// Enforce constraints: reset pinned nodes to their prescribed positions
    /// and zero their velocities.
    pub fn enforce(&self, positions: &mut [[f64; 3]], velocities: &mut [[f64; 3]]) {
        for &(idx, pos) in &self.fixed {
            positions[idx] = pos;
            velocities[idx] = [0.0; 3];
        }
    }

    /// Build a `pinned` boolean mask of length `n_nodes`.
    pub fn pinned_mask(&self, n_nodes: usize) -> Vec<bool> {
        let mut mask = vec![false; n_nodes];
        for &(idx, _) in &self.fixed {
            if idx < n_nodes {
                mask[idx] = true;
            }
        }
        mask
    }
}

impl Default for BoundaryConditions {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CorotFemElement4: 4-node tetrahedral element (explicit struct, new name)
// ---------------------------------------------------------------------------

/// A 4-node tetrahedral corotational FEM element storing Young's modulus and
/// Poisson's ratio directly (so callers don't need to pre-compute Lame params).
#[derive(Debug, Clone)]
pub struct CorotFemElement4 {
    /// Vertex indices (into a flat node array).
    pub indices: [usize; 4],
    /// Rest volume (m^3).
    pub rest_volume: f64,
    /// Inverse of the rest-shape edge matrix Dm.
    dm_inv: [[f64; 3]; 3],
    /// Young's modulus (Pa).
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Derived shear modulus mu.
    mu: f64,
    /// Derived first Lame parameter lambda.
    lambda: f64,
}

impl CorotFemElement4 {
    /// Create a new element from node positions, Young's modulus, and Poisson's ratio.
    pub fn new(indices: [usize; 4], positions: &[[f64; 3]], young: f64, poisson: f64) -> Self {
        let p0 = positions[indices[0]];
        let p1 = positions[indices[1]];
        let p2 = positions[indices[2]];
        let p3 = positions[indices[3]];
        let dm = edge_matrix_raw(p0, p1, p2, p3);
        let rest_volume = det3x3(dm).abs() / 6.0;
        let dm_inv = inv3x3(dm);
        let mu = young / (2.0 * (1.0 + poisson));
        let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        Self {
            indices,
            rest_volume,
            dm_inv,
            young,
            poisson,
            mu,
            lambda,
        }
    }

    /// Compute the deformation gradient F = Ds * Dm^{-1}.
    pub fn deformation_gradient(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let p0 = positions[self.indices[0]];
        let p1 = positions[self.indices[1]];
        let p2 = positions[self.indices[2]];
        let p3 = positions[self.indices[3]];
        let ds = edge_matrix_raw(p0, p1, p2, p3);
        mul3x3(ds, self.dm_inv)
    }

    /// Compute corotational elastic forces on the four nodes.
    pub fn elastic_forces(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 4] {
        let f = self.deformation_gradient(positions);
        let r = CorotFemTet::polar_decompose_rotation(f);
        let rt = transpose3x3(r);
        let rtf = mul3x3(rt, f);
        let mut strain = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                strain[i][j] = rtf[i][j] - delta;
            }
        }
        let trace = strain[0][0] + strain[1][1] + strain[2][2];
        let mut stress = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                stress[i][j] = 2.0 * self.mu * strain[i][j] + self.lambda * trace * delta;
            }
        }
        let piola = mul3x3(r, stress);
        let dm_inv_t = transpose3x3(self.dm_inv);
        let h = mul3x3(piola, dm_inv_t);
        let mut forces = [[0.0f64; 3]; 4];
        for i in 0..3 {
            forces[1][i] = -self.rest_volume * h[i][0];
            forces[2][i] = -self.rest_volume * h[i][1];
            forces[3][i] = -self.rest_volume * h[i][2];
            forces[0][i] = -(forces[1][i] + forces[2][i] + forces[3][i]);
        }
        forces
    }

    /// Compute corotational stiffness matrix using finite differences (12x12).
    pub fn stiffness_matrix(&self, positions: &[[f64; 3]]) -> Vec<f64> {
        let h = 1e-6;
        let dof = 12;
        let mut k = vec![0.0_f64; dof * dof];
        let f0 = self.elastic_forces(positions);
        for col_node in 0..4 {
            for col_dim in 0..3 {
                let col = col_node * 3 + col_dim;
                let mut pos_p = positions.to_vec();
                pos_p[self.indices[col_node]][col_dim] += h;
                let f1 = self.elastic_forces(&pos_p);
                for row_node in 0..4 {
                    for row_dim in 0..3 {
                        let row = row_node * 3 + row_dim;
                        k[row * dof + col] = (f1[row_node][row_dim] - f0[row_node][row_dim]) / h;
                    }
                }
            }
        }
        k
    }

    /// Elastic strain energy of this element.
    pub fn strain_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let f = self.deformation_gradient(positions);
        let r = CorotFemTet::polar_decompose_rotation(f);
        let rt = transpose3x3(r);
        let rtf = mul3x3(rt, f);
        let mut strain = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                strain[i][j] = rtf[i][j] - delta;
            }
        }
        let trace = strain[0][0] + strain[1][1] + strain[2][2];
        let mut frob_sq = 0.0;
        for row in strain.iter() {
            for &x in row.iter() {
                frob_sq += x * x;
            }
        }
        self.rest_volume * (self.mu * frob_sq + 0.5 * self.lambda * trace * trace)
    }
}

// ---------------------------------------------------------------------------
// FemSoftBodyVerlet: complete FEM soft body with Velocity Verlet + Rayleigh
// ---------------------------------------------------------------------------

/// A complete FEM soft body simulation using:
/// - `CorotFemElement4` for corotational tetrahedral elements
/// - Velocity Verlet time integration
/// - Rayleigh damping
/// - Fixed boundary conditions
#[derive(Debug, Clone)]
pub struct FemSoftBodyVerlet {
    /// Node positions.
    pub positions: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Node masses (lumped).
    pub masses: Vec<f64>,
    /// Tetrahedral elements.
    pub elements: Vec<CorotFemElement4>,
    /// Velocity Verlet integrator state.
    pub integrator: VelocityVerletIntegrator,
    /// Rayleigh damping parameters.
    pub rayleigh: RayleighDamping,
    /// Boundary conditions.
    pub bc: BoundaryConditions,
    /// Gravitational acceleration.
    pub gravity: [f64; 3],
}

impl FemSoftBodyVerlet {
    /// Create a new FEM soft body.
    pub fn new(
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
        elements: Vec<CorotFemElement4>,
        gravity: [f64; 3],
        dt: f64,
        rayleigh_alpha: f64,
        rayleigh_beta_k: f64,
    ) -> Self {
        let n = positions.len();
        Self {
            velocities: vec![[0.0; 3]; n],
            integrator: VelocityVerletIntegrator::new(n, dt),
            positions,
            masses,
            elements,
            rayleigh: RayleighDamping::new(rayleigh_alpha, rayleigh_beta_k),
            bc: BoundaryConditions::new(),
            gravity,
        }
    }

    /// Advance the simulation by one time step.
    pub fn step(&mut self) {
        let n = self.positions.len();
        let pinned = self.bc.pinned_mask(n);
        let dt = self.integrator.dt;

        // 1. Predict positions: x_{n+1} = x_n + dt*v + 0.5*dt^2*a_prev
        for (i, (pos, (vel, (pinned_i, a)))) in self
            .positions
            .iter_mut()
            .zip(
                self.velocities
                    .iter()
                    .zip(pinned.iter().zip(self.integrator.prev_accel.iter())),
            )
            .enumerate()
        {
            let _ = i;
            if *pinned_i {
                continue;
            }
            for (p_d, (v_d, a_d)) in pos.iter_mut().zip(vel.iter().zip(a.iter())) {
                *p_d += dt * v_d + 0.5 * dt * dt * a_d;
            }
        }

        // Enforce boundary conditions after position prediction
        self.bc.enforce(&mut self.positions, &mut self.velocities);

        // 2. Compute forces at new positions
        let mut forces = vec![[0.0_f64; 3]; n];

        // Gravity
        let gravity = self.gravity;
        for (i, (force, (pinned_i, mass))) in forces
            .iter_mut()
            .zip(pinned.iter().zip(self.masses.iter()))
            .enumerate()
        {
            let _ = i;
            if *pinned_i {
                continue;
            }
            for (f_d, g_d) in force.iter_mut().zip(gravity.iter()) {
                *f_d += mass * g_d;
            }
        }

        // Elastic (restoring) forces
        let elastic = assemble_internal_forces(&self.positions, &self.elements);
        for (force, f_el) in forces.iter_mut().zip(elastic.iter()) {
            for d in 0..3 {
                force[d] += f_el[d];
            }
        }

        // Rayleigh damping.  The stiffness term uses K v of the elastic forces
        // alone, never the gravity-inclusive force array.
        let kv = if self.rayleigh.uses_stiffness_term() {
            let moving: Vec<[f64; 3]> = self
                .velocities
                .iter()
                .zip(pinned.iter())
                .map(|(v, &fixed)| if fixed { [0.0; 3] } else { *v })
                .collect();
            let elements = &self.elements;
            RayleighDamping::stiffness_velocity_product(&self.positions, &moving, |x| {
                assemble_internal_forces(x, elements)
            })
        } else {
            vec![[0.0_f64; 3]; n]
        };
        self.rayleigh
            .apply(&mut forces, &self.velocities, &self.masses, &kv);

        // Zero forces on pinned nodes
        for i in 0..n {
            if pinned[i] {
                forces[i] = [0.0; 3];
            }
        }

        // 3. New accelerations
        let new_accel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                if pinned[i] || self.masses[i] < 1e-30 {
                    return [0.0; 3];
                }
                let inv_m = 1.0 / self.masses[i];
                [
                    forces[i][0] * inv_m,
                    forces[i][1] * inv_m,
                    forces[i][2] * inv_m,
                ]
            })
            .collect();

        // 4. Correct velocities: v_{n+1} = v_n + 0.5*dt*(a_prev + a_new)
        for i in 0..n {
            if pinned[i] {
                continue;
            }
            let a_old = self.integrator.prev_accel[i];
            let a_new = new_accel[i];
            for d in 0..3 {
                self.velocities[i][d] += 0.5 * dt * (a_old[d] + a_new[d]);
            }
        }

        self.integrator.prev_accel = new_accel;

        // Enforce BCs again after velocity update
        self.bc.enforce(&mut self.positions, &mut self.velocities);
    }

    /// Pin a node at its current position.
    pub fn pin_node(&mut self, idx: usize) {
        let pos = self.positions[idx];
        self.bc.pin_at(idx, pos);
        self.velocities[idx] = [0.0; 3];
    }

    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.positions.len() {
            let v = self.velocities[i];
            ke += 0.5 * self.masses[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        ke
    }

    /// Total elastic strain energy.
    pub fn strain_energy(&self) -> f64 {
        self.elements
            .iter()
            .map(|e| e.strain_energy(&self.positions))
            .sum()
    }

    /// Total mechanical energy (kinetic + elastic).
    pub fn mechanical_energy(&self) -> f64 {
        self.kinetic_energy() + self.strain_energy()
    }
}

// ---------------------------------------------------------------------------
// Internal force vector (assembled from all elements)
// ---------------------------------------------------------------------------

/// Assemble the global internal force vector from all tetrahedral elements.
///
/// Returns a `Vec<[f64; 3]>` of length `n_nodes`.
pub fn assemble_internal_forces(
    positions: &[[f64; 3]],
    elements: &[CorotFemElement4],
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut forces = vec![[0.0_f64; 3]; n];
    for elem in elements {
        let fs = elem.elastic_forces(positions);
        for (k, fs_k) in fs.iter().enumerate() {
            let idx = elem.indices[k];
            for (force_d, fs_d) in forces[idx].iter_mut().zip(fs_k.iter()) {
                *force_d += fs_d;
            }
        }
    }
    forces
}
