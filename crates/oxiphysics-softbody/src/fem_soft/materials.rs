// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Material models for FEM soft bodies: SoftMaterial, NeoHookean, HyperelasticBody.
//!
//! Uses raw f64 arrays (no nalgebra dependency).

use super::corot_functions::{corot_cauchy_stress, corotational_strain, polar_decompose_r};
use super::math_helpers::{det3x3, edge_matrix_raw, inv3x3, mul3x3, transpose3x3};

// ---------------------------------------------------------------------------
// Standalone FEM types (f64 array-based, no nalgebra dependency)
// ---------------------------------------------------------------------------

/// Material properties for a soft body element.
#[derive(Debug, Clone, Copy)]
pub struct SoftMaterial {
    /// Young's modulus (stiffness).
    pub youngs_modulus: f64,
    /// Poisson's ratio (0..0.5).
    pub poisson_ratio: f64,
    /// Kelvin-Voigt damping time \[s\].
    ///
    /// The viscous stress is `damping` times the elastic stress of the strain
    /// *rate* (stiffness-proportional damping), so rigid motions are not damped.
    pub damping: f64,
}

impl SoftMaterial {
    /// Lame parameters `(lambda, mu)` of the linear-elastic law of this material:
    /// `lambda = E nu / ((1 + nu)(1 - 2 nu))`, `mu = E / (2 (1 + nu))`.
    ///
    /// The Poisson ratio is clamped to `[-0.99, 0.4999]` so that `lambda` stays
    /// finite as the material approaches the incompressible limit `nu = 1/2`.
    pub fn lame_parameters(&self) -> (f64, f64) {
        let nu = self.poisson_ratio.clamp(-0.99, 0.4999);
        let e = self.youngs_modulus;
        (
            e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)),
            e / (2.0 * (1.0 + nu)),
        )
    }
}

/// A corotational tetrahedral element using raw f64 arrays.
///
/// With the deformation gradient `F = Ds Dm^-1` (current and rest edge matrices
/// have the edge vectors as columns) and its polar decomposition `F = R S`, the
/// element stores the corotational linear-elastic energy of the strain
/// `eps = R^T F - I`.  It is force free in its rest shape and under every rigid
/// motion, also for skewed rest geometries.
#[derive(Debug, Clone)]
pub struct CorotationalElementRaw {
    /// Indices of the four nodes of the tetrahedron.
    pub node_indices: [usize; 4],
    /// Rest volume of the element.
    pub rest_volume: f64,
    /// Rest-shape edge vectors `[e1, e2, e3]` with `ei = pi - p0`, one edge per
    /// **row** (the rest edge matrix `Dm` holds them as columns).
    pub rest_shape: [[f64; 3]; 3],
    /// Inverse `Dm^-1` of the rest edge matrix, stored at construction.
    pub rest_inv: [[f64; 3]; 3],
    /// Material properties.
    pub material: SoftMaterial,
}

impl CorotationalElementRaw {
    /// Create a new element from node positions.
    pub fn new(node_indices: [usize; 4], positions: &[[f64; 3]], material: SoftMaterial) -> Self {
        let p0 = positions[node_indices[0]];
        let p1 = positions[node_indices[1]];
        let p2 = positions[node_indices[2]];
        let p3 = positions[node_indices[3]];
        let dm = edge_matrix_raw(p0, p1, p2, p3);
        Self {
            node_indices,
            rest_volume: det3x3(dm).abs() / 6.0,
            rest_shape: transpose3x3(dm),
            rest_inv: inv3x3(dm),
            material,
        }
    }

    /// Deformation gradient `F = Ds Dm^-1` for the given node positions.
    pub fn deformation_gradient(&self, current_positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let [i0, i1, i2, i3] = self.node_indices;
        let ds = edge_matrix_raw(
            current_positions[i0],
            current_positions[i1],
            current_positions[i2],
            current_positions[i3],
        );
        mul3x3(ds, self.rest_inv)
    }

    /// Rotation `R` of the polar decomposition `F = R S` of the deformation
    /// gradient `F = Ds Dm^-1`.
    ///
    /// The rotation is a property of `F`, not of the raw current edges: in a
    /// skewed rest shape the raw edge matrix already differs from a rotation at
    /// rest.
    pub fn compute_rotation(&self, current_positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
        polar_decompose_r(self.deformation_gradient(current_positions))
    }

    /// Nodal forces `f_k` of the first Piola-Kirchhoff stress `P = R sigma`,
    /// `sigma = lambda tr(eps) I + 2 mu eps`, `eps = R^T F - I`.
    ///
    /// `H = -V0 P Dm^-T` holds the forces of nodes 1..3 as its columns and the
    /// force of node 0 balances them: `f_0 = -(f_1 + f_2 + f_3)`.
    pub fn compute_forces(&self, current_positions: &[[f64; 3]]) -> [[f64; 3]; 4] {
        let f = self.deformation_gradient(current_positions);
        let r = polar_decompose_r(f);
        let strain = corotational_strain(r, f);
        let (lambda, mu) = self.material.lame_parameters();
        let sigma = corot_cauchy_stress(strain, lambda, mu);
        self.nodal_forces_from_stress(mul3x3(r, sigma), -self.rest_volume)
    }

    /// Viscous (Kelvin-Voigt) forces `-V0 R (damping C : eps_dot) Dm^-T`, where
    /// `C` is the elastic law and `eps_dot` the strain rate of the corotated
    /// element for the nodal `velocities`.
    ///
    /// The strain rate is `S_dot` of the polar decomposition `F = R S`,
    /// `S_dot = R^T F_dot - [w]x S`, with the rotation rate `w` solving
    /// `(tr(S) I - S) w = 2 axial(skew(R^T F_dot))`, so a rigid motion carries
    /// no strain rate at any deformation.  The forces are zero for
    /// `material.damping == 0`.
    pub fn damping_forces(
        &self,
        current_positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> [[f64; 3]; 4] {
        let eta = self.material.damping;
        if eta == 0.0 || self.rest_volume == 0.0 {
            return [[0.0; 3]; 4];
        }
        let f = self.deformation_gradient(current_positions);
        let r = polar_decompose_r(f);
        let rt = transpose3x3(r);
        let s = mul3x3(rt, f);
        let [i0, i1, i2, i3] = self.node_indices;
        let dv = edge_matrix_raw(
            velocities[i0],
            velocities[i1],
            velocities[i2],
            velocities[i3],
        );
        let a = mul3x3(rt, mul3x3(dv, self.rest_inv));
        let axial = [
            0.5 * (a[2][1] - a[1][2]),
            0.5 * (a[0][2] - a[2][0]),
            0.5 * (a[1][0] - a[0][1]),
        ];
        let trace_s = s[0][0] + s[1][1] + s[2][2];
        let mut m = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] = -s[i][j];
            }
            m[i][i] += trace_s;
        }
        let w = if det3x3(m).abs() > 1e-12 {
            let m_inv = inv3x3(m);
            let mut w = [0.0_f64; 3];
            for i in 0..3 {
                w[i] = 2.0
                    * (m_inv[i][0] * axial[0] + m_inv[i][1] * axial[1] + m_inv[i][2] * axial[2]);
            }
            w
        } else {
            [0.0; 3]
        };
        let w_cross = [[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]];
        let ws = mul3x3(w_cross, s);
        let mut strain_rate = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                strain_rate[i][j] = 0.5 * ((a[i][j] - ws[i][j]) + (a[j][i] - ws[j][i]));
            }
        }
        let (lambda, mu) = self.material.lame_parameters();
        let sigma_rate = corot_cauchy_stress(strain_rate, lambda, mu);
        self.nodal_forces_from_stress(mul3x3(r, sigma_rate), -self.rest_volume * eta)
    }

    /// Nodal forces `scale * P Dm^-T` (columns for nodes 1..3, balancing node 0).
    fn nodal_forces_from_stress(&self, piola: [[f64; 3]; 3], scale: f64) -> [[f64; 3]; 4] {
        let h = mul3x3(piola, transpose3x3(self.rest_inv));
        let mut forces = [[0.0_f64; 3]; 4];
        for i in 0..3 {
            forces[1][i] = scale * h[i][0];
            forces[2][i] = scale * h[i][1];
            forces[3][i] = scale * h[i][2];
            forces[0][i] = -(forces[1][i] + forces[2][i] + forces[3][i]);
        }
        forces
    }
}

/// A FEM soft body using raw f64 arrays.
#[derive(Debug, Clone)]
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

    /// Perform one explicit (semi-implicit Euler) time step.
    ///
    /// Both the elastic and the Kelvin-Voigt force are evaluated at the start of
    /// the step, so the step is only conditionally stable: for a mode of angular
    /// frequency `omega` of the elements, with `eta` the
    /// [`SoftMaterial::damping`] time, `omega^2 (dt^2 + 2 eta dt) < 4` has to
    /// hold. Without damping this is the usual `dt < 2 / omega`; a large damping
    /// time lowers the admissible step to about `dt < 2 / (eta omega^2)`.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.nodes.len();
        let mut forces = vec![[0.0f64; 3]; n];

        // Gravity
        for (force, mass) in forces.iter_mut().zip(self.masses.iter()) {
            for (f_d, g_d) in force.iter_mut().zip(gravity.iter()) {
                *f_d += mass * g_d;
            }
        }

        // Element forces: elastic restoring forces plus Kelvin-Voigt damping.
        for elem in &self.elements {
            let f = elem.compute_forces(&self.nodes);
            let f_damp = elem.damping_forces(&self.nodes, &self.velocities);
            for (k, (f_k, fd_k)) in f.iter().zip(f_damp.iter()).enumerate() {
                let idx = elem.node_indices[k];
                for (fd, (fkd, fdd)) in forces[idx].iter_mut().zip(f_k.iter().zip(fd_k.iter())) {
                    *fd += fkd + fdd;
                }
            }
        }

        // Integrate
        for (i, (node, (vel, (force, mass)))) in self
            .nodes
            .iter_mut()
            .zip(
                self.velocities
                    .iter_mut()
                    .zip(forces.iter().zip(self.masses.iter())),
            )
            .enumerate()
        {
            let _ = i;
            let inv_m = 1.0 / mass;
            for (n_d, (v_d, f_d)) in node.iter_mut().zip(vel.iter_mut().zip(force.iter())) {
                *v_d += f_d * inv_m * dt;
                *n_d += *v_d * dt;
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
pub struct NeoHookeanMaterial {
    /// First Lame parameter mu (shear modulus).
    pub mu: f64,
    /// Second Lame parameter lambda.
    pub lambda: f64,
}

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
        let i1: f64 = f.iter().flat_map(|row| row.iter()).map(|x| x * x).sum();
        let ln_j = det.abs().ln();
        0.5 * self.mu * (i1 - 3.0) - self.mu * ln_j + 0.5 * self.lambda * ln_j * ln_j
    }
}

/// A Neo-Hookean tetrahedral element using raw f64 arrays.
#[derive(Debug, Clone)]
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

        for (force, mass) in forces.iter_mut().zip(self.masses.iter()) {
            for (f_d, g_d) in force.iter_mut().zip(gravity.iter()) {
                *f_d += mass * g_d;
            }
        }

        for elem in &self.elements {
            let f = elem.compute_forces(&self.nodes);
            for (k, f_k) in f.iter().enumerate() {
                let idx = elem.node_indices[k];
                for (fd, fkd) in forces[idx].iter_mut().zip(f_k.iter()) {
                    *fd += fkd;
                }
            }
        }

        let damping = self.damping;
        for (node, (vel, (force, mass))) in self.nodes.iter_mut().zip(
            self.velocities
                .iter_mut()
                .zip(forces.iter().zip(self.masses.iter())),
        ) {
            let inv_m = 1.0 / mass;
            for (n_d, (v_d, f_d)) in node.iter_mut().zip(vel.iter_mut().zip(force.iter())) {
                *v_d += f_d * inv_m * dt;
                *v_d *= 1.0 - damping;
                *n_d += *v_d * dt;
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

// ---------------------------------------------------------------------------
// Stable Neo-Hookean (Smith et al. 2018)
// ---------------------------------------------------------------------------

/// Cofactor matrix of a 3x3 matrix (row-major).
///
/// `cofactor[i][j] = (-1)^(i+j) * minor_ij`. This equals `dJ/dF` for
/// `J = det(F)` and is computed from 2x2 minors so it stays finite as
/// `J -> 0` (unlike `J * F^{-T}`).
fn cofactor3(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            f[1][1] * f[2][2] - f[1][2] * f[2][1],
            -(f[1][0] * f[2][2] - f[1][2] * f[2][0]),
            f[1][0] * f[2][1] - f[1][1] * f[2][0],
        ],
        [
            -(f[0][1] * f[2][2] - f[0][2] * f[2][1]),
            f[0][0] * f[2][2] - f[0][2] * f[2][0],
            -(f[0][0] * f[2][1] - f[0][1] * f[2][0]),
        ],
        [
            f[0][1] * f[1][2] - f[0][2] * f[1][1],
            -(f[0][0] * f[1][2] - f[0][2] * f[1][0]),
            f[0][0] * f[1][1] - f[0][1] * f[1][0],
        ],
    ]
}

/// Stable Neo-Hookean material (Smith et al. 2018, "Stable Neo-Hookean Flesh
/// Simulation"). Free of volumetric locking at high Poisson ratios and stable
/// through inversion.
#[derive(Debug, Clone, Copy)]
pub struct StableNeoHookeanMaterial {
    /// Effective shear constant used internally. NOTE: `from_young_poisson`
    /// stores the Smith-2018 remapped value `(4/3)*mu`, not the textbook Lame
    /// `mu`, so the stable model matches classic linear elasticity at small
    /// strain. Construct via `from_young_poisson` unless you have already
    /// applied the remap yourself.
    pub mu: f64,
    /// Effective volumetric constant used internally. NOTE: `from_young_poisson`
    /// stores the Smith-2018 remapped value `lambda + (5/6)*mu`, not the
    /// textbook Lame `lambda`.
    pub lambda: f64,
}

impl StableNeoHookeanMaterial {
    /// Create from Young's modulus and Poisson's ratio.
    ///
    /// The Smith et al. (2018) "Stable Neo-Hookean Flesh Simulation" energy does
    /// not linearize to the standard Lame constants when fed the classic
    /// `mu`/`lambda`. The paper derives a remapping so the stable model
    /// reproduces the same small-strain (linear-elastic) limit as the classic
    /// Neo-Hookean: `mu_hat = (4/3) mu` and `lambda_hat = lambda + (5/6) mu`.
    /// These hat values are stored so that around `F = I` the stable stress
    /// matches classic linear elasticity.
    pub fn from_young_poisson(young: f64, poisson: f64) -> Self {
        let mu = young / (2.0 * (1.0 + poisson));
        let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        let mu_hat = 4.0 / 3.0 * mu;
        let lambda_hat = lambda + 5.0 / 6.0 * mu;
        Self {
            mu: mu_hat,
            lambda: lambda_hat,
        }
    }

    /// Rest-stability constant `alpha = 1 + mu/lambda - mu/(4 lambda)`.
    ///
    /// Chosen so the undeformed state `F = I` is a stress-free energy minimum.
    pub fn alpha(&self) -> f64 {
        1.0 + self.mu / self.lambda - self.mu / (4.0 * self.lambda)
    }

    /// Strain energy density `Psi(F)`.
    pub fn strain_energy(&self, f: [[f64; 3]; 3]) -> f64 {
        let i_c: f64 = f.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
        let j = det3x3(f);
        let a = self.alpha();
        0.5 * self.mu * (i_c - 3.0) + 0.5 * self.lambda * (j - a) * (j - a)
            - 0.5 * self.mu * (i_c + 1.0).ln()
    }

    /// First Piola-Kirchhoff stress `P = dPsi/dF`.
    pub fn piola_kirchhoff(&self, f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let i_c: f64 = f.iter().flat_map(|r| r.iter()).map(|x| x * x).sum();
        let j = det3x3(f);
        let a = self.alpha();
        let cof = cofactor3(f);
        let mu_term = self.mu * (1.0 - 1.0 / (i_c + 1.0));
        let vol_term = self.lambda * (j - a);
        let mut p = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                p[i][k] = mu_term * f[i][k] + vol_term * cof[i][k];
            }
        }
        p
    }
}

/// A stable Neo-Hookean tetrahedral element (raw f64 arrays).
#[derive(Debug, Clone)]
pub struct StableNeoHookeanElement {
    /// Node indices.
    pub node_indices: [usize; 4],
    /// Rest-shape inverse of the edge matrix (Dm^{-1}), row-major.
    pub rest_inv: [[f64; 3]; 3],
    /// Rest volume.
    pub rest_volume: f64,
    /// Material.
    pub material: StableNeoHookeanMaterial,
}

impl StableNeoHookeanElement {
    /// Create a new element from node positions.
    pub fn new(
        node_indices: [usize; 4],
        positions: &[[f64; 3]],
        material: StableNeoHookeanMaterial,
    ) -> Self {
        let p0 = positions[node_indices[0]];
        let p1 = positions[node_indices[1]];
        let p2 = positions[node_indices[2]];
        let p3 = positions[node_indices[3]];
        let dm = edge_matrix_raw(p0, p1, p2, p3);
        let rest_volume = det3x3(dm).abs() / 6.0;
        let rest_inv = inv3x3(dm);
        Self {
            node_indices,
            rest_inv,
            rest_volume,
            material,
        }
    }

    /// Deformation gradient `F = Ds * Dm^{-1}`.
    pub fn deformation_gradient(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let p0 = positions[self.node_indices[0]];
        let p1 = positions[self.node_indices[1]];
        let p2 = positions[self.node_indices[2]];
        let p3 = positions[self.node_indices[3]];
        let ds = edge_matrix_raw(p0, p1, p2, p3);
        mul3x3(ds, self.rest_inv)
    }

    /// Elastic forces on the four nodes.
    pub fn compute_forces(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 4] {
        let f = self.deformation_gradient(positions);
        let p = self.material.piola_kirchhoff(f);
        let dm_inv_t = transpose3x3(self.rest_inv);
        let h = mul3x3(p, dm_inv_t);
        let mut forces = [[0.0f64; 3]; 4];
        for i in 0..3 {
            forces[1][i] = -self.rest_volume * h[i][0];
            forces[2][i] = -self.rest_volume * h[i][1];
            forces[3][i] = -self.rest_volume * h[i][2];
            forces[0][i] = -(forces[1][i] + forces[2][i] + forces[3][i]);
        }
        forces
    }

    /// Strain energy of this element.
    pub fn strain_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let f = self.deformation_gradient(positions);
        self.rest_volume * self.material.strain_energy(f)
    }
}

#[cfg(test)]
mod tests_stable_neohookean {
    use super::*;

    fn identity() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    #[test]
    fn test_piola_finite_difference() {
        let mat = StableNeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        // A generic non-symmetric deformation gradient.
        let f = [[1.1, 0.05, -0.02], [0.03, 0.95, 0.04], [-0.01, 0.02, 1.05]];
        let analytic = mat.piola_kirchhoff(f);
        let eps = 1e-6;
        let mut max_err = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                let mut fp = f;
                let mut fm = f;
                fp[i][j] += eps;
                fm[i][j] -= eps;
                let num = (mat.strain_energy(fp) - mat.strain_energy(fm)) / (2.0 * eps);
                let err = (num - analytic[i][j]).abs();
                max_err = max_err.max(err);
                assert!(
                    err < 1e-4,
                    "P[{i}][{j}] fd={num} analytic={} err={err}",
                    analytic[i][j]
                );
            }
        }
        assert!(max_err < 1e-4, "max P fd error {max_err}");
    }

    #[test]
    fn test_rest_state_stress_free() {
        // At F = I the stable model is (near) stress-free by construction of alpha.
        let mat = StableNeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        let p = mat.piola_kirchhoff(identity());
        for (i, row) in p.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1e-9, "rest stress P[{i}][{j}]={val}");
            }
        }
    }

    #[test]
    fn test_stable_vs_classic_small_strain() {
        // Small stretch: both Neo-Hookean variants converge to identical linear
        // elasticity, so their stresses must agree to ~1%.
        let young = 1000.0;
        let poisson = 0.3;
        let stable = StableNeoHookeanMaterial::from_young_poisson(young, poisson);
        let classic = NeoHookeanMaterial::from_young_poisson(young, poisson);
        let eps = 1e-4;
        let f = [
            [1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps],
        ];
        let ps = stable.piola_kirchhoff(f);
        let pc = classic.piola_kirchhoff(f);
        // Compare Frobenius norms of the stress (both ~ O(eps * young)).
        let mut ns = 0.0;
        let mut nc = 0.0;
        let mut nd = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                ns += ps[i][j] * ps[i][j];
                nc += pc[i][j] * pc[i][j];
                let d = ps[i][j] - pc[i][j];
                nd += d * d;
            }
        }
        let rel = nd.sqrt() / nc.sqrt().max(1e-30);
        assert!(
            rel < 0.01,
            "small-strain stress disagreement {rel} (stable|{}| classic|{}|)",
            ns.sqrt(),
            nc.sqrt()
        );
    }

    #[test]
    fn test_stable_neohookean_element_forces() {
        // Forces must equal the negative gradient of total strain energy.
        let rest = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mat = StableNeoHookeanMaterial::from_young_poisson(500.0, 0.3);
        let elem = StableNeoHookeanElement::new([0, 1, 2, 3], &rest, mat);
        // Apply a stretch.
        let mut pos = rest.clone();
        pos[1][0] = 1.3;
        pos[2][1] = 1.15;
        pos[3][2] = 0.9;

        let forces = elem.compute_forces(&pos);
        let eps = 1e-6;
        let mut max_err = 0.0f64;
        for node in 0..4 {
            for d in 0..3 {
                let mut pp = pos.clone();
                let mut pm = pos.clone();
                pp[node][d] += eps;
                pm[node][d] -= eps;
                // Force = -dEnergy/dx.
                let num = -(elem.strain_energy(&pp) - elem.strain_energy(&pm)) / (2.0 * eps);
                let err = (num - forces[node][d]).abs();
                max_err = max_err.max(err);
                assert!(
                    err < 1e-2,
                    "force node {node} dim {d}: fd={num} analytic={} err={err}",
                    forces[node][d]
                );
            }
        }
        assert!(max_err < 1e-2);
    }

    #[test]
    fn test_stable_no_locking() {
        // Near-incompressible (nu=0.499) should deflect comparably to nu=0.45
        // under the same shear/stretch (no volumetric locking), unlike the
        // classic compressible Neo-Hookean which stiffens (locks).
        let young = 1000.0;
        let stable_high = StableNeoHookeanMaterial::from_young_poisson(young, 0.499);
        let stable_low = StableNeoHookeanMaterial::from_young_poisson(young, 0.45);
        // An isochoric-ish shear deformation (det ~ 1).
        let f = [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Deviatoric (shear) stress magnitude should be governed by mu, which is
        // nearly equal for both nu, so the stable model's shear response barely
        // changes with nu (no locking).
        let ph = stable_high.piola_kirchhoff(f);
        let pl = stable_low.piola_kirchhoff(f);
        let shear_h = ph[0][1].abs();
        let shear_l = pl[0][1].abs();
        let rel = (shear_h - shear_l).abs() / shear_l.max(1e-30);
        assert!(
            rel < 0.2,
            "stable shear response changed too much with nu (locking): {rel}"
        );
    }
}
