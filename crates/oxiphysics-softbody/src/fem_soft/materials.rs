// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Material models for FEM soft bodies: SoftMaterial, NeoHookean, HyperelasticBody.
//!
//! Uses raw f64 arrays (no nalgebra dependency).

#![allow(clippy::needless_range_loop)]

use super::math_helpers::{det3x3, edge_matrix_raw, inv3x3, mul3x3, transpose3x3};

// ---------------------------------------------------------------------------
// Standalone FEM types (f64 array-based, no nalgebra dependency)
// ---------------------------------------------------------------------------

/// Material properties for a soft body element.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct SoftMaterial {
    /// Young's modulus (stiffness).
    pub youngs_modulus: f64,
    /// Poisson's ratio (0..0.5).
    pub poisson_ratio: f64,
    /// Damping coefficient.
    pub damping: f64,
}

/// A corotational tetrahedral element using raw f64 arrays.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CorotationalElementRaw {
    /// Indices of the four nodes of the tetrahedron.
    pub node_indices: [usize; 4],
    /// Rest volume of the element.
    pub rest_volume: f64,
    /// Rest-shape edge matrix columns: `[e1, e2, e3]` where `ei = pi - p0`.
    pub rest_shape: [[f64; 3]; 3],
    /// Material properties.
    pub material: SoftMaterial,
}

#[allow(dead_code)]
impl CorotationalElementRaw {
    /// Create a new element from node positions.
    pub fn new(node_indices: [usize; 4], positions: &[[f64; 3]], material: SoftMaterial) -> Self {
        let p0 = positions[node_indices[0]];
        let p1 = positions[node_indices[1]];
        let p2 = positions[node_indices[2]];
        let p3 = positions[node_indices[3]];
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let e3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
        let det = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
            + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
        let rest_volume = det.abs() / 6.0;
        Self {
            node_indices,
            rest_volume,
            rest_shape: [e1, e2, e3],
            material,
        }
    }

    /// Compute a 3x3 rotation matrix from the deformation gradient via
    /// iterative polar decomposition.
    pub fn compute_rotation(current_positions: &[[f64; 3]], indices: &[usize; 4]) -> [[f64; 3]; 3] {
        let p0 = current_positions[indices[0]];
        let p1 = current_positions[indices[1]];
        let p2 = current_positions[indices[2]];
        let p3 = current_positions[indices[3]];

        let ds = [
            [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]],
            [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]],
            [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]],
        ];

        // Start with ds as initial guess for R.
        let mut r = ds;
        for _ in 0..10 {
            let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
                - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
                + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
            if det.abs() < 1e-30 {
                return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            }
            let inv_det = 1.0 / det;
            let inv = [
                [
                    (r[1][1] * r[2][2] - r[1][2] * r[2][1]) * inv_det,
                    (r[0][2] * r[2][1] - r[0][1] * r[2][2]) * inv_det,
                    (r[0][1] * r[1][2] - r[0][2] * r[1][1]) * inv_det,
                ],
                [
                    (r[1][2] * r[2][0] - r[1][0] * r[2][2]) * inv_det,
                    (r[0][0] * r[2][2] - r[0][2] * r[2][0]) * inv_det,
                    (r[0][2] * r[1][0] - r[0][0] * r[1][2]) * inv_det,
                ],
                [
                    (r[1][0] * r[2][1] - r[1][1] * r[2][0]) * inv_det,
                    (r[0][1] * r[2][0] - r[0][0] * r[2][1]) * inv_det,
                    (r[0][0] * r[1][1] - r[0][1] * r[1][0]) * inv_det,
                ],
            ];
            // R = (R + inv^T) / 2
            for i in 0..3 {
                for j in 0..3 {
                    r[i][j] = (r[i][j] + inv[j][i]) * 0.5;
                }
            }
        }
        r
    }

    /// Compute elastic forces on the four nodes.
    pub fn compute_forces(&self, current_positions: &[[f64; 3]]) -> [[f64; 3]; 4] {
        let r = Self::compute_rotation(current_positions, &self.node_indices);
        // Simplified: force proportional to displacement from rest shape
        let p0 = current_positions[self.node_indices[0]];
        let mut forces = [[0.0; 3]; 4];
        let stiffness = self.material.youngs_modulus * self.rest_volume;

        for k in 1..4 {
            let pk = current_positions[self.node_indices[k]];
            let deformed = [pk[0] - p0[0], pk[1] - p0[1], pk[2] - p0[2]];
            // Rotated rest edge
            let rest_e = self.rest_shape[k - 1];
            let rotated = [
                r[0][0] * rest_e[0] + r[0][1] * rest_e[1] + r[0][2] * rest_e[2],
                r[1][0] * rest_e[0] + r[1][1] * rest_e[1] + r[1][2] * rest_e[2],
                r[2][0] * rest_e[0] + r[2][1] * rest_e[1] + r[2][2] * rest_e[2],
            ];
            let diff = [
                deformed[0] - rotated[0],
                deformed[1] - rotated[1],
                deformed[2] - rotated[2],
            ];
            let f = [
                -stiffness * diff[0],
                -stiffness * diff[1],
                -stiffness * diff[2],
            ];
            forces[k] = f;
            forces[0][0] -= f[0];
            forces[0][1] -= f[1];
            forces[0][2] -= f[2];
        }
        forces
    }
}

/// A FEM soft body using raw f64 arrays.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FemSoftBodyRaw {
    /// Node positions.
    pub nodes: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Tetrahedral elements.
    pub elements: Vec<CorotationalElementRaw>,
    /// Mass per node.
    pub masses: Vec<f64>,
}

#[allow(dead_code)]
impl FemSoftBodyRaw {
    /// Create a new FEM soft body.
    pub fn new(nodes: Vec<[f64; 3]>, elements: Vec<CorotationalElementRaw>, mass: f64) -> Self {
        let n = nodes.len();
        Self {
            velocities: vec![[0.0; 3]; n],
            nodes,
            elements,
            masses: vec![mass; n],
        }
    }

    /// Perform one explicit Euler time step.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.nodes.len();
        let mut forces = vec![[0.0f64; 3]; n];

        // Gravity
        for i in 0..n {
            forces[i][0] += self.masses[i] * gravity[0];
            forces[i][1] += self.masses[i] * gravity[1];
            forces[i][2] += self.masses[i] * gravity[2];
        }

        // Element forces
        for elem in &self.elements {
            let f = elem.compute_forces(&self.nodes);
            for k in 0..4 {
                let idx = elem.node_indices[k];
                forces[idx][0] += f[k][0];
                forces[idx][1] += f[k][1];
                forces[idx][2] += f[k][2];
            }
        }

        // Integrate
        for i in 0..n {
            let inv_m = 1.0 / self.masses[i];
            for d in 0..3 {
                self.velocities[i][d] += forces[i][d] * inv_m * dt;
                self.nodes[i][d] += self.velocities[i][d] * dt;
            }
        }
    }

    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.nodes.len() {
            let v = self.velocities[i];
            ke += 0.5 * self.masses[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        ke
    }

    /// Total gravitational potential energy.
    pub fn potential_energy(&self, gravity: [f64; 3]) -> f64 {
        let g_mag =
            (gravity[0] * gravity[0] + gravity[1] * gravity[1] + gravity[2] * gravity[2]).sqrt();
        if g_mag < 1e-30 {
            return 0.0;
        }
        let mut pe = 0.0;
        for i in 0..self.nodes.len() {
            let dot = self.nodes[i][0] * gravity[0]
                + self.nodes[i][1] * gravity[1]
                + self.nodes[i][2] * gravity[2];
            pe += -self.masses[i] * dot;
        }
        pe
    }
}

// ---------------------------------------------------------------------------
// Hyperelastic soft body (Neo-Hookean, raw f64)
// ---------------------------------------------------------------------------

/// Neo-Hookean hyperelastic material model.
///
/// Strain energy density: W = (mu/2)(I1 - 3) - mu ln(J) + (lambda/2)(ln J)^2
/// where I1 = tr(F^T F), J = det(F).
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct NeoHookeanMaterial {
    /// First Lame parameter mu (shear modulus).
    pub mu: f64,
    /// Second Lame parameter lambda.
    pub lambda: f64,
}

#[allow(dead_code)]
impl NeoHookeanMaterial {
    /// Create from Young's modulus and Poisson's ratio.
    pub fn from_young_poisson(young: f64, poisson: f64) -> Self {
        let mu = young / (2.0 * (1.0 + poisson));
        let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        Self { mu, lambda }
    }

    /// Compute the first Piola-Kirchhoff stress tensor P = dW/dF.
    ///
    /// For Neo-Hookean: P = mu(F - F^{-T}) + lambda ln(J) F^{-T}
    ///
    /// Returns P as `[[P00,P01,P02\],[P10,P11,P12],[P20,P21,P22]]`.
    pub fn piola_kirchhoff(&self, deformation_gradient: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let f = deformation_gradient;
        let det = det3x3(f);
        if det.abs() < 1e-30 {
            return [[0.0; 3]; 3];
        }
        let f_inv_t = super::math_helpers::inv3x3_transpose(f);
        let ln_j = det.abs().ln();

        let mut p = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                p[i][j] = self.mu * (f[i][j] - f_inv_t[i][j]) + self.lambda * ln_j * f_inv_t[i][j];
            }
        }
        p
    }

    /// Compute strain energy density W for a given deformation gradient.
    pub fn strain_energy(&self, f: [[f64; 3]; 3]) -> f64 {
        let det = det3x3(f);
        if det.abs() < 1e-30 {
            return 0.0;
        }
        // I1 = tr(F^T F)
        let mut i1 = 0.0;
        for i in 0..3 {
            for k in 0..3 {
                i1 += f[k][i] * f[k][i];
            }
        }
        let ln_j = det.abs().ln();
        0.5 * self.mu * (i1 - 3.0) - self.mu * ln_j + 0.5 * self.lambda * ln_j * ln_j
    }
}

/// A Neo-Hookean tetrahedral element using raw f64 arrays.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NeoHookeanElement {
    /// Indices of the four nodes.
    pub node_indices: [usize; 4],
    /// Rest-shape inverse of the edge matrix (Dm^{-1}).
    pub rest_inv: [[f64; 3]; 3],
    /// Rest volume.
    pub rest_volume: f64,
    /// Material.
    pub material: NeoHookeanMaterial,
}

#[allow(dead_code)]
impl NeoHookeanElement {
    /// Create a new element from node positions.
    pub fn new(
        node_indices: [usize; 4],
        positions: &[[f64; 3]],
        material: NeoHookeanMaterial,
    ) -> Self {
        let p0 = positions[node_indices[0]];
        let p1 = positions[node_indices[1]];
        let p2 = positions[node_indices[2]];
        let p3 = positions[node_indices[3]];
        let dm = edge_matrix_raw(p0, p1, p2, p3);
        let det = det3x3(dm);
        let rest_volume = det.abs() / 6.0;
        let rest_inv = inv3x3(dm);
        Self {
            node_indices,
            rest_inv,
            rest_volume,
            material,
        }
    }

    /// Compute the deformation gradient F = Ds * Dm^{-1}.
    pub fn deformation_gradient(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let p0 = positions[self.node_indices[0]];
        let p1 = positions[self.node_indices[1]];
        let p2 = positions[self.node_indices[2]];
        let p3 = positions[self.node_indices[3]];
        let ds = edge_matrix_raw(p0, p1, p2, p3);
        mul3x3(ds, self.rest_inv)
    }

    /// Compute elastic forces on the four nodes.
    pub fn compute_forces(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 4] {
        let f = self.deformation_gradient(positions);
        let p = self.material.piola_kirchhoff(f);
        // H = -V * P * Dm^{-T}
        let dm_inv_t = transpose3x3(self.rest_inv);
        let h = mul3x3(p, dm_inv_t);

        let mut forces = [[0.0; 3]; 4];
        for i in 0..3 {
            forces[1][i] = -self.rest_volume * h[i][0];
            forces[2][i] = -self.rest_volume * h[i][1];
            forces[3][i] = -self.rest_volume * h[i][2];
            forces[0][i] = -(forces[1][i] + forces[2][i] + forces[3][i]);
        }
        forces
    }

    /// Compute the strain energy of this element.
    pub fn strain_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let f = self.deformation_gradient(positions);
        self.rest_volume * self.material.strain_energy(f)
    }
}

/// A hyperelastic FEM body using Neo-Hookean elements.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HyperelasticBody {
    /// Node positions.
    pub nodes: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Elements.
    pub elements: Vec<NeoHookeanElement>,
    /// Mass per node.
    pub masses: Vec<f64>,
    /// Damping coefficient.
    pub damping: f64,
}

#[allow(dead_code)]
impl HyperelasticBody {
    /// Create a new hyperelastic body.
    pub fn new(
        nodes: Vec<[f64; 3]>,
        elements: Vec<NeoHookeanElement>,
        mass: f64,
        damping: f64,
    ) -> Self {
        let n = nodes.len();
        Self {
            velocities: vec![[0.0; 3]; n],
            nodes,
            elements,
            masses: vec![mass; n],
            damping,
        }
    }

    /// Perform one explicit Euler time step.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.nodes.len();
        let mut forces = vec![[0.0f64; 3]; n];

        for i in 0..n {
            for d in 0..3 {
                forces[i][d] += self.masses[i] * gravity[d];
            }
        }

        for elem in &self.elements {
            let f = elem.compute_forces(&self.nodes);
            for k in 0..4 {
                let idx = elem.node_indices[k];
                for d in 0..3 {
                    forces[idx][d] += f[k][d];
                }
            }
        }

        for i in 0..n {
            let inv_m = 1.0 / self.masses[i];
            for d in 0..3 {
                self.velocities[i][d] += forces[i][d] * inv_m * dt;
                self.velocities[i][d] *= 1.0 - self.damping;
                self.nodes[i][d] += self.velocities[i][d] * dt;
            }
        }
    }

    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.nodes.len() {
            let v = self.velocities[i];
            ke += 0.5 * self.masses[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        ke
    }

    /// Total strain energy.
    pub fn strain_energy(&self) -> f64 {
        self.elements
            .iter()
            .map(|e| e.strain_energy(&self.nodes))
            .sum()
    }
}
