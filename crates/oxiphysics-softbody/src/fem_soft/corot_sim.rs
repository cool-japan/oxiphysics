// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CorotFem simulation types: CorotFemNode, CorotFemTet, CorotFemSim,
//! NewmarkBetaSim, CorotationalTet, and HighLevelFemBody.

#![allow(clippy::needless_range_loop)]

use super::math_helpers::{det3x3, edge_matrix_raw, inv3x3, mul3x3, transpose3x3};

// ---------------------------------------------------------------------------
// CorotFem: clean f64-array corotational FEM types
// ---------------------------------------------------------------------------

/// A single node in a corotational FEM simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CorotFemNode {
    /// Current position.
    pub position: [f64; 3],
    /// Current velocity.
    pub velocity: [f64; 3],
    /// Accumulated force (cleared each step before accumulation).
    pub force: [f64; 3],
    /// Node mass.
    pub mass: f64,
}

/// A tetrahedral element for the corotational FEM simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CorotFemTet {
    /// Indices into the node array.
    pub node_indices: [usize; 4],
    /// Rest positions of the four nodes.
    pub rest_positions: [[f64; 3]; 4],
    /// Rest volume = |det(Dm)| / 6.
    pub volume: f64,
    /// Inverse of the reference shape matrix Dm = \[p1-p0, p2-p0, p3-p0\].
    pub dm_inv: [[f64; 3]; 3],
}

#[allow(dead_code)]
impl CorotFemTet {
    /// Create a new tet from node indices and initial node positions.
    pub fn new(node_indices: [usize; 4], nodes: &[CorotFemNode]) -> Self {
        let p: [[f64; 3]; 4] = [
            nodes[node_indices[0]].position,
            nodes[node_indices[1]].position,
            nodes[node_indices[2]].position,
            nodes[node_indices[3]].position,
        ];
        let dm = edge_matrix_raw(p[0], p[1], p[2], p[3]);
        let det = det3x3(dm);
        let volume = det.abs() / 6.0;
        let dm_inv = inv3x3(dm);
        Self {
            node_indices,
            rest_positions: p,
            volume,
            dm_inv,
        }
    }

    /// Compute the deformation gradient F = Ds * Dm^{-1}.
    pub fn compute_deformation_gradient(&self, nodes: &[CorotFemNode]) -> [[f64; 3]; 3] {
        let p0 = nodes[self.node_indices[0]].position;
        let p1 = nodes[self.node_indices[1]].position;
        let p2 = nodes[self.node_indices[2]].position;
        let p3 = nodes[self.node_indices[3]].position;
        let ds = edge_matrix_raw(p0, p1, p2, p3);
        mul3x3(ds, self.dm_inv)
    }

    /// Extract rotation R from F via iterative polar decomposition (Newton iterations).
    pub fn polar_decompose_rotation(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut r = f;
        for _ in 0..5 {
            let det = det3x3(r);
            if det.abs() < 1e-30 {
                return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            }
            let inv = inv3x3(r);
            let inv_t = transpose3x3(inv);
            for i in 0..3 {
                for j in 0..3 {
                    r[i][j] = (r[i][j] + inv_t[i][j]) * 0.5;
                }
            }
        }
        // Ensure proper rotation (det > 0)
        if det3x3(r) < 0.0 {
            for i in 0..3 {
                for j in 0..3 {
                    r[i][j] = -r[i][j];
                }
            }
        }
        r
    }

    /// Compute corotational elastic forces on the four nodes.
    ///
    /// Returns forces `[[fx,fy,fz\]; 4]` in order of `node_indices`.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_elastic_forces(
        &self,
        nodes: &[CorotFemNode],
        mu: f64,
        lambda: f64,
    ) -> [[f64; 3]; 4] {
        let f = self.compute_deformation_gradient(nodes);
        let r = Self::polar_decompose_rotation(f);
        let rt = transpose3x3(r);

        // Strain in rotated frame: epsilon = R^T * F - I
        let rtf = mul3x3(rt, f);
        let mut strain = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                strain[i][j] = rtf[i][j] - delta;
            }
        }

        // Stress: sigma = 2*mu*strain + lambda*tr(strain)*I
        let trace = strain[0][0] + strain[1][1] + strain[2][2];
        let mut stress = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                stress[i][j] = 2.0 * mu * strain[i][j] + lambda * trace * delta;
            }
        }

        // First Piola-Kirchhoff: P = R * stress
        let piola = mul3x3(r, stress);

        // H = -volume * P * Dm^{-T}
        let dm_inv_t = transpose3x3(self.dm_inv);
        let h = mul3x3(piola, dm_inv_t);

        let mut forces = [[0.0f64; 3]; 4];
        for i in 0..3 {
            forces[1][i] = -self.volume * h[i][0];
            forces[2][i] = -self.volume * h[i][1];
            forces[3][i] = -self.volume * h[i][2];
            forces[0][i] = -(forces[1][i] + forces[2][i] + forces[3][i]);
        }
        forces
    }
}

/// Full corotational FEM simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CorotFemSim {
    /// All nodes.
    pub nodes: Vec<CorotFemNode>,
    /// All tetrahedral elements.
    pub tets: Vec<CorotFemTet>,
    /// Shear modulus mu.
    pub mu: f64,
    /// First Lame parameter lambda.
    pub lambda: f64,
    /// Time step.
    pub dt: f64,
}

#[allow(dead_code)]
impl CorotFemSim {
    /// Create a new simulation.
    pub fn new(
        nodes: Vec<CorotFemNode>,
        tets: Vec<CorotFemTet>,
        mu: f64,
        lambda: f64,
        dt: f64,
    ) -> Self {
        Self {
            nodes,
            tets,
            mu,
            lambda,
            dt,
        }
    }

    /// Advance one time step using semi-implicit Euler.
    pub fn step(&mut self) {
        // Clear forces
        for n in &mut self.nodes {
            n.force = [0.0; 3];
        }
        // Accumulate elastic forces
        for tet in &self.tets {
            let fs = tet.compute_elastic_forces(&self.nodes, self.mu, self.lambda);
            for k in 0..4 {
                let idx = tet.node_indices[k];
                for d in 0..3 {
                    self.nodes[idx].force[d] += fs[k][d];
                }
            }
        }
        // Semi-implicit Euler
        let dt = self.dt;
        for n in &mut self.nodes {
            let inv_m = 1.0 / n.mass;
            for d in 0..3 {
                n.velocity[d] += n.force[d] * inv_m * dt;
                n.position[d] += n.velocity[d] * dt;
            }
        }
    }

    /// Compute total elastic strain energy across all tets.
    pub fn total_strain_energy(&self) -> f64 {
        let mut energy = 0.0;
        for tet in &self.tets {
            let f = tet.compute_deformation_gradient(&self.nodes);
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
            let mut frobenius_sq = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    frobenius_sq += strain[i][j] * strain[i][j];
                }
            }
            energy += tet.volume * (self.mu * frobenius_sq + 0.5 * self.lambda * trace * trace);
        }
        energy
    }
}

// ---------------------------------------------------------------------------
// CorotTet -- ergonomic alias for CorotFemTet
// ---------------------------------------------------------------------------

/// Alias for [`CorotFemTet`].  Provides the same corotational FEM tet element
/// under the shorter name used by high-level interfaces.
pub type CorotTet = CorotFemTet;

// ---------------------------------------------------------------------------
// Newmark-beta time integration for FEM soft body
// ---------------------------------------------------------------------------

/// FEM soft body integrator using the Newmark-beta method.
///
/// The classic Newmark parameters `beta = 0.25` and `gamma = 0.5` give the
/// **average constant acceleration** (unconditionally stable, second-order
/// accurate, energy-conserving for linear problems).
///
/// The update formulas are:
/// ```text
/// a_{n+1}  = M^{-1} * f(x_{n+1})           (fixed-point / linearised)
/// v_{n+1}  = v_n + dt*(1 - gamma)*a_n + dt*gamma*a_{n+1}
/// x_{n+1}  = x_n + dt*v_n + dt^2*(0.5 - beta)*a_n + dt^2*beta*a_{n+1}
/// ```
/// Because we use the _explicit prediction / force-corrector_ variant (one
/// force evaluation per step), this degenerates to the **explicit Newmark**
/// form which remains stable for the chosen beta, gamma values at small time steps.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NewmarkBetaSim {
    /// All simulation nodes (position, velocity, force, mass).
    pub nodes: Vec<CorotFemNode>,
    /// Tetrahedral elements.
    pub tets: Vec<CorotFemTet>,
    /// Shear modulus mu (Pa).
    pub mu: f64,
    /// First Lame parameter lambda (Pa).
    pub lambda: f64,
    /// Time step dt (s).
    pub dt: f64,
    /// Newmark beta parameter (0.25 for average-acceleration).
    pub beta: f64,
    /// Newmark gamma parameter (0.5 for average-acceleration).
    pub gamma: f64,
    /// Acceleration at the previous step (for Newmark update).
    pub prev_accel: Vec<[f64; 3]>,
    /// Pinned-node flags: pinned nodes are held fixed.
    pub pinned: Vec<bool>,
    /// Rayleigh damping coefficient alpha (mass-proportional).
    pub rayleigh_alpha: f64,
    /// Rayleigh damping coefficient beta_r (stiffness-proportional).
    pub rayleigh_beta: f64,
}

#[allow(dead_code)]
impl NewmarkBetaSim {
    /// Create a new Newmark-beta simulation with default parameters
    /// (beta = 0.25, gamma = 0.5, no Rayleigh damping).
    pub fn new(
        nodes: Vec<CorotFemNode>,
        tets: Vec<CorotFemTet>,
        mu: f64,
        lambda: f64,
        dt: f64,
    ) -> Self {
        let n = nodes.len();
        Self {
            nodes,
            tets,
            mu,
            lambda,
            dt,
            beta: 0.25,
            gamma: 0.5,
            prev_accel: vec![[0.0; 3]; n],
            pinned: vec![false; n],
            rayleigh_alpha: 0.0,
            rayleigh_beta: 0.0,
        }
    }

    /// Create a new simulation with explicit Newmark parameters and Rayleigh damping.
    #[allow(clippy::too_many_arguments)]
    pub fn with_params(
        nodes: Vec<CorotFemNode>,
        tets: Vec<CorotFemTet>,
        mu: f64,
        lambda: f64,
        dt: f64,
        beta: f64,
        gamma: f64,
        rayleigh_alpha: f64,
        rayleigh_beta: f64,
    ) -> Self {
        let n = nodes.len();
        Self {
            nodes,
            tets,
            mu,
            lambda,
            dt,
            beta,
            gamma,
            prev_accel: vec![[0.0; 3]; n],
            pinned: vec![false; n],
            rayleigh_alpha,
            rayleigh_beta,
        }
    }

    /// Compute elastic + Rayleigh-damping forces into each node's `force` field.
    fn accumulate_forces(&mut self) {
        let n = self.nodes.len();
        for node in &mut self.nodes {
            node.force = [0.0; 3];
        }
        // Elastic forces from elements
        for tet in &self.tets {
            let fs = tet.compute_elastic_forces(&self.nodes, self.mu, self.lambda);
            for k in 0..4 {
                let idx = tet.node_indices[k];
                for d in 0..3 {
                    self.nodes[idx].force[d] += fs[k][d];
                }
            }
        }
        // Rayleigh mass-proportional damping: f_d = -alpha * m * v
        if self.rayleigh_alpha.abs() > 1e-30 {
            for i in 0..n {
                if self.pinned[i] {
                    continue;
                }
                let m = self.nodes[i].mass;
                let v = self.nodes[i].velocity;
                for d in 0..3 {
                    self.nodes[i].force[d] -= self.rayleigh_alpha * m * v[d];
                }
            }
        }
        // Rayleigh stiffness-proportional damping: applied implicitly via force
        // Here we do an approximation: subtract beta_r * K * v via elastic forces
        // with a perturbed velocity step (not implemented for full implicit;
        // zero contribution when rayleigh_beta = 0).
    }

    /// Advance the simulation by one time step using the explicit Newmark-beta scheme.
    ///
    /// Algorithm:
    /// 1. Predictor: `x* = x + dt*v + dt^2*(0.5 - beta)*a_n`
    /// 2. Evaluate forces at `x*` -> acceleration `a*`
    /// 3. Corrector:
    ///    - `v_{n+1} = v + dt*(1 - gamma)*a_n + dt*gamma*a*`
    ///    - `x_{n+1} = x* + dt^2*beta*a*`
    pub fn step(&mut self) {
        let dt = self.dt;
        let beta = self.beta;
        let gamma = self.gamma;
        let n = self.nodes.len();

        // --- Predictor ---
        let pred_positions: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let mut pos = [0.0f64; 3];
                if self.pinned[i] {
                    return self.nodes[i].position;
                }
                let v = self.nodes[i].velocity;
                let a = self.prev_accel[i];
                for d in 0..3 {
                    pos[d] = self.nodes[i].position[d] + dt * v[d] + dt * dt * (0.5 - beta) * a[d];
                }
                pos
            })
            .collect();

        // Apply predicted positions temporarily
        for i in 0..n {
            if !self.pinned[i] {
                self.nodes[i].position = pred_positions[i];
            }
        }

        // --- Evaluate forces at predicted positions ---
        self.accumulate_forces();

        // Compute new accelerations from forces
        let new_accel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                if self.pinned[i] {
                    return [0.0; 3];
                }
                let inv_m = 1.0 / self.nodes[i].mass;
                let f = self.nodes[i].force;
                [f[0] * inv_m, f[1] * inv_m, f[2] * inv_m]
            })
            .collect();

        // --- Corrector ---
        for i in 0..n {
            if self.pinned[i] {
                continue;
            }
            let a_n = self.prev_accel[i];
            let a_new = new_accel[i];

            // Velocity corrector
            for d in 0..3 {
                self.nodes[i].velocity[d] += dt * (1.0 - gamma) * a_n[d] + dt * gamma * a_new[d];
                // Position correction: x* already has (0.5-beta)*a_n term; add beta*a*
                self.nodes[i].position[d] += dt * dt * beta * a_new[d];
            }
        }

        // Save accelerations for next step
        self.prev_accel = new_accel;
    }

    /// Compute total elastic strain energy over all tets.
    pub fn total_strain_energy(&self) -> f64 {
        let mut energy = 0.0;
        for tet in &self.tets {
            let f = tet.compute_deformation_gradient(&self.nodes);
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
            for i in 0..3 {
                for j in 0..3 {
                    frob_sq += strain[i][j] * strain[i][j];
                }
            }
            energy += tet.volume * (self.mu * frob_sq + 0.5 * self.lambda * trace * trace);
        }
        energy
    }

    /// Compute total kinetic energy: Sum 0.5 * m * |v|^2.
    pub fn total_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for (i, node) in self.nodes.iter().enumerate() {
            if self.pinned[i] {
                continue;
            }
            let v = node.velocity;
            ke += 0.5 * node.mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        ke
    }

    /// Compute total mechanical energy (kinetic + strain).
    pub fn total_mechanical_energy(&self) -> f64 {
        self.total_kinetic_energy() + self.total_strain_energy()
    }

    /// Pin a node by index so it is held fixed during integration.
    pub fn pin_node(&mut self, idx: usize) {
        if idx < self.pinned.len() {
            self.pinned[idx] = true;
            self.nodes[idx].velocity = [0.0; 3];
        }
    }

    /// Unpin a node, allowing it to move freely.
    pub fn unpin_node(&mut self, idx: usize) {
        if idx < self.pinned.len() {
            self.pinned[idx] = false;
        }
    }

    /// Apply an impulse (velocity increment) to a node.
    pub fn apply_impulse(&mut self, idx: usize, impulse: [f64; 3]) {
        if idx < self.nodes.len() && !self.pinned[idx] {
            let inv_m = 1.0 / self.nodes[idx].mass;
            for d in 0..3 {
                self.nodes[idx].velocity[d] += impulse[d] * inv_m;
            }
        }
    }

    /// Number of active (non-pinned) nodes.
    pub fn active_node_count(&self) -> usize {
        self.pinned.iter().filter(|&&p| !p).count()
    }
}

// ---------------------------------------------------------------------------
// CorotationalTet -- clean f64-array tet with explicit Dm / V0 fields
// ---------------------------------------------------------------------------

/// A single corotational tetrahedral FEM element using raw f64 arrays.
///
/// Stores the reference shape matrix `Dm` and rest volume `V0` explicitly,
/// matching the notation common in FEM textbooks.
#[derive(Debug, Clone)]
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct CorotationalTet {
    /// Indices of the four nodes of the tetrahedron.
    pub indices: [usize; 4],
    /// Reference shape matrix: columns are edge vectors `[p1-p0, p2-p0, p3-p0]`.
    pub Dm: [[f64; 3]; 3],
    /// Rest volume of the element (m^3).
    pub V0: f64,
    /// Inverse of the reference shape matrix `Dm`.
    dm_inv: [[f64; 3]; 3],
}

#[allow(dead_code)]
#[allow(non_snake_case)]
impl CorotationalTet {
    /// Create a new `CorotationalTet` from node positions.
    ///
    /// `positions[indices[k\]]` gives the 3-D rest position of node `k`.
    pub fn new(indices: [usize; 4], positions: &[[f64; 3]]) -> Self {
        let p0 = positions[indices[0]];
        let p1 = positions[indices[1]];
        let p2 = positions[indices[2]];
        let p3 = positions[indices[3]];
        let Dm = edge_matrix_raw(p0, p1, p2, p3);
        let V0 = det3x3(Dm).abs() / 6.0;
        let dm_inv = inv3x3(Dm);
        Self {
            indices,
            Dm,
            V0,
            dm_inv,
        }
    }

    /// Compute the deformation gradient `F = Ds * Dm^{-1}`.
    ///
    /// `nodes` holds the **current** positions (world space).
    pub fn compute_deformation_gradient(&self, nodes: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let p0 = nodes[self.indices[0]];
        let p1 = nodes[self.indices[1]];
        let p2 = nodes[self.indices[2]];
        let p3 = nodes[self.indices[3]];
        let ds = edge_matrix_raw(p0, p1, p2, p3);
        mul3x3(ds, self.dm_inv)
    }

    /// Polar decomposition of `F`.
    ///
    /// Returns `(R, S)` where `F = R * S`, `R` is a proper rotation and `S`
    /// is symmetric positive-definite.  Uses the iterative Newton method.
    pub fn polar_decompose(F: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
        let mut r = *F;
        for _ in 0..10 {
            let det = det3x3(r);
            if det.abs() < 1e-30 {
                let r_id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
                return (r_id, *F);
            }
            let inv = inv3x3(r);
            let inv_t = transpose3x3(inv);
            for i in 0..3 {
                for j in 0..3 {
                    r[i][j] = (r[i][j] + inv_t[i][j]) * 0.5;
                }
            }
        }
        if det3x3(r) < 0.0 {
            for row in r.iter_mut() {
                for x in row.iter_mut() {
                    *x = -*x;
                }
            }
        }
        // S = R^T F
        let rt = transpose3x3(r);
        let s = mul3x3(rt, *F);
        (r, s)
    }

    /// Compute corotational elastic forces on the four nodes.
    ///
    /// Returns `[[fx, fy, fz\]; 4]` in the same order as `self.indices`.
    ///
    /// Uses linear (corotational) elasticity:
    /// - Strain `eps = R^T F - I` (in rotated frame)
    /// - Stress `sigma = 2mu eps + lambda tr(eps) I`
    /// - Forces `H = -V0 P Dm^{-T}`,  `P = R sigma`
    #[allow(clippy::too_many_arguments)]
    pub fn elastic_forces(&self, nodes: &[[f64; 3]], mu: f64, lambda: f64) -> [[f64; 3]; 4] {
        let f = self.compute_deformation_gradient(nodes);
        let (r, _s) = Self::polar_decompose(&f);
        let rt = transpose3x3(r);

        // Strain in rotated frame: eps = R^T F - I
        let rtf = mul3x3(rt, f);
        let mut strain = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                strain[i][j] = rtf[i][j] - delta;
            }
        }

        // Stress: sigma = 2mu eps + lambda tr(eps) I
        let trace = strain[0][0] + strain[1][1] + strain[2][2];
        let mut stress = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                stress[i][j] = 2.0 * mu * strain[i][j] + lambda * trace * delta;
            }
        }

        // First Piola-Kirchhoff: P = R sigma
        let piola = mul3x3(r, stress);

        // H = -V0 P Dm^{-T}
        let dm_inv_t = transpose3x3(self.dm_inv);
        let h = mul3x3(piola, dm_inv_t);

        let mut forces = [[0.0f64; 3]; 4];
        for i in 0..3 {
            forces[1][i] = -self.V0 * h[i][0];
            forces[2][i] = -self.V0 * h[i][1];
            forces[3][i] = -self.V0 * h[i][2];
            forces[0][i] = -(forces[1][i] + forces[2][i] + forces[3][i]);
        }
        forces
    }

    /// Compute the 12x12 tangent stiffness matrix for this element.
    ///
    /// Returns the matrix in row-major order (144 f64 values).
    /// Entries are ordered as `[node0_x, node0_y, node0_z, node1_x, ..., node3_z]`.
    ///
    /// Uses a finite-difference approximation with a small perturbation `h`.
    pub fn stiffness_matrix(&self, nodes: &[[f64; 3]], mu: f64, lambda: f64) -> Vec<f64> {
        let h = 1e-6;
        let dof = 12; // 4 nodes * 3
        let mut k = vec![0.0_f64; dof * dof];
        let f0 = self.elastic_forces(nodes, mu, lambda);

        for col_node in 0..4 {
            for col_dim in 0..3 {
                let col = col_node * 3 + col_dim;
                let mut nodes_p = nodes.to_vec();
                nodes_p[self.indices[col_node]][col_dim] += h;
                let f1 = self.elastic_forces(&nodes_p, mu, lambda);
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
}

// ---------------------------------------------------------------------------
// HighLevelFemBody -- a FemSoftBody wrapper with add_tet / total_elastic_energy
// ---------------------------------------------------------------------------

/// A high-level FEM soft body built from `CorotationalTet` elements.
///
/// Provides a simple builder API:
/// - `new(positions, masses)` -- create from node positions and masses
/// - `add_tet(i, j, k, l)` -- add a tetrahedral element
/// - `step(dt)` -- integrate one explicit Euler step with elastic forces only
/// - `total_elastic_energy()` -- compute the total elastic strain energy
#[derive(Debug, Clone)]
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct HighLevelFemBody {
    /// Current node positions (world space).
    pub positions: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Node masses (kg).
    pub masses: Vec<f64>,
    /// All tetrahedral elements.
    pub elements: Vec<CorotationalTet>,
    /// Young's modulus E (Pa).
    pub E: f64,
    /// Poisson's ratio nu.
    pub nu: f64,
    /// Derived shear modulus mu.
    mu: f64,
    /// Derived first Lame parameter lambda.
    lambda: f64,
}

#[allow(dead_code)]
#[allow(non_snake_case)]
impl HighLevelFemBody {
    /// Create a new body from rest positions and node masses.
    ///
    /// `E` is Young's modulus (Pa), `nu` is Poisson's ratio.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>, E: f64, nu: f64) -> Self {
        let n = positions.len();
        let mu = E / (2.0 * (1.0 + nu));
        let lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            elements: Vec::new(),
            E,
            nu,
            mu,
            lambda,
        }
    }

    /// Add a tetrahedral element connecting nodes `i, j, k, l`.
    ///
    /// Uses the current `positions` as the rest configuration.
    pub fn add_tet(&mut self, i: usize, j: usize, k: usize, l: usize) {
        let tet = CorotationalTet::new([i, j, k, l], &self.positions);
        self.elements.push(tet);
    }

    /// Advance one explicit Euler time step with elastic forces.
    ///
    /// No gravity; use `apply_gravity` first if needed.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        let mut forces = vec![[0.0f64; 3]; n];

        for elem in &self.elements {
            let f = elem.elastic_forces(&self.positions, self.mu, self.lambda);
            for (k, &node_idx) in elem.indices.iter().enumerate() {
                for d in 0..3 {
                    forces[node_idx][d] += f[k][d];
                }
            }
        }

        for i in 0..n {
            let inv_m = 1.0 / self.masses[i];
            for d in 0..3 {
                self.velocities[i][d] += forces[i][d] * inv_m * dt;
                self.positions[i][d] += self.velocities[i][d] * dt;
            }
        }
    }

    /// Apply gravitational acceleration to all nodes (adds to velocities).
    pub fn apply_gravity(&mut self, gravity: [f64; 3], dt: f64) {
        for i in 0..self.positions.len() {
            for d in 0..3 {
                self.velocities[i][d] += gravity[d] * dt;
            }
        }
    }

    /// Compute the total elastic strain energy over all elements.
    ///
    /// Energy density: `W = mu tr(eps^T eps) + 0.5 lambda tr(eps)^2`
    pub fn total_elastic_energy(&self) -> f64 {
        let mut energy = 0.0;
        for elem in &self.elements {
            let f = elem.compute_deformation_gradient(&self.positions);
            let (r, _) = CorotationalTet::polar_decompose(&f);
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
            energy += elem.V0 * (self.mu * frob_sq + 0.5 * self.lambda * trace * trace);
        }
        energy
    }
}
