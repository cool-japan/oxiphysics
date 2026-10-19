// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Force model and explicit integrator of the electroactive tetrahedral FEM body.
//!
//! Each [`ElectroacitiveTetraElement`] carries two kinds of nodal forces:
//!
//! - **Elastic forces** of a compressible neo-Hookean solid,
//!   `Ψ = μ/2 (tr FᵀF − 3) − μ ln J + λ/2 (ln J)²` with the Lamé constant
//!   `λ = κ − 2μ/3` (so that the small-strain bulk modulus is the stored `κ` and
//!   the shear modulus the stored `μ`).  The deformation gradient of the linear
//!   tetrahedron is `F = Σ_k x_k ⊗ ∇N_k`, the first Piola–Kirchhoff stress
//!   `P = μ (F − F⁻ᵀ) + λ ln J F⁻ᵀ` and the nodal forces `f_k = −V₀ P ∇N_k`.
//! - **Maxwell forces** of the electric field, `f_k = −v T ∇N_k` with the Maxwell
//!   stress `T = ε (E ⊗ E − |E|²/2 I)`.  The *deformed* variant evaluates the field
//!   and the gradients on the current geometry (a thinner film has a stronger
//!   field, the origin of the pull-in instability of dielectric elastomers); at
//!   fixed nodal potentials it is the exact derivative of the stored electric
//!   energy `½ ε v |E|²` with respect to the node positions.
//!
//! [`ElectroacitiveFemBody::step`] accumulates both force sets and integrates the
//! equations of motion with semi-implicit (symplectic) Euler.

use super::{
    EPSILON_0, ElectroacitiveFemBody, ElectroacitiveFemNode, ElectroacitiveTetraElement,
    invert_3x3, v3_dot, v3_sub,
};

/// Row-major 3×3 matrix.
type Mat3 = [[f64; 3]; 3];

/// Smallest volume ratio `J = det F` at which the neo-Hookean law is evaluated.
///
/// Below it (a nearly flat or inverted element) `ln J` and `F⁻ᵀ = cof(F) / J` use
/// this value instead, which keeps the forces finite and still pushes the element
/// back towards positive volume (the cofactor is the gradient of `J`).
const MIN_JACOBIAN: f64 = 1.0e-3;

/// Cofactor matrix `cof(M) = det(M) M⁻ᵀ`.
fn cofactor(m: &Mat3) -> Mat3 {
    [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
        ],
        [
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
        ],
        [
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ]
}

/// Determinant of a 3×3 matrix.
fn determinant(m: &Mat3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// `M v` for a 3×3 matrix and a 3-vector.
fn mat_vec(m: &Mat3, v: [f64; 3]) -> [f64; 3] {
    [v3_dot(m[0], v), v3_dot(m[1], v), v3_dot(m[2], v)]
}

/// Volume and shape-function gradients `∇N₀ … ∇N₃` of a tetrahedron with the
/// given corner positions, or `None` for a (nearly) degenerate tetrahedron whose
/// volume is below `min_volume`.
fn tet_geometry(corners: [[f64; 3]; 4], min_volume: f64) -> Option<(f64, [[f64; 3]; 4])> {
    let edges = [
        v3_sub(corners[1], corners[0]),
        v3_sub(corners[2], corners[0]),
        v3_sub(corners[3], corners[0]),
    ];
    let volume = determinant(&edges).abs() / 6.0;
    if volume.is_nan() || volume <= min_volume {
        return None;
    }
    // Rows of `edges` are the edge vectors, so column k - 1 of its inverse is ∇N_k.
    let inverse = invert_3x3(edges);
    let mut grads = [[0.0_f64; 3]; 4];
    for k in 1..4 {
        for d in 0..3 {
            grads[k][d] = inverse[d][k - 1];
            grads[0][d] -= inverse[d][k - 1];
        }
    }
    Some((volume, grads))
}

impl ElectroacitiveTetraElement {
    /// First Lamé constant `λ = κ − 2μ/3` of the neo-Hookean law.
    pub fn lame_lambda(&self) -> f64 {
        self.bulk_modulus - 2.0 * self.shear_modulus / 3.0
    }

    /// Deformation gradient `F = Σ_k x_k ⊗ ∇N_k` of the element for the current
    /// node positions (`F = I` in the rest configuration).
    pub fn deformation_gradient(&self, nodes: &[ElectroacitiveFemNode]) -> [[f64; 3]; 3] {
        let grads = self.shape_gradients();
        let mut f = [[0.0_f64; 3]; 3];
        for (k, grad) in grads.iter().enumerate() {
            let x = nodes[self.node_indices[k]].position;
            for i in 0..3 {
                for j in 0..3 {
                    f[i][j] += x[i] * grad[j];
                }
            }
        }
        f
    }

    /// First Piola–Kirchhoff stress `P = μ F + (λ ln J − μ) cof(F) / J` of the
    /// compressible neo-Hookean law (`= μ (F − F⁻ᵀ) + λ ln J F⁻ᵀ`) for the
    /// deformation gradient `f`; `J` is clamped from below by [`MIN_JACOBIAN`].
    pub fn first_piola_stress(&self, f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mu = self.shear_modulus;
        let j = determinant(&f).max(MIN_JACOBIAN);
        let cof = cofactor(&f);
        let coefficient = (self.lame_lambda() * j.ln() - mu) / j;
        let mut p = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                p[i][k] = mu * f[i][k] + coefficient * cof[i][k];
            }
        }
        p
    }

    /// Strain energy `V₀ (μ/2 (tr FᵀF − 3) − μ ln J + λ/2 (ln J)²)` \[J\] of the
    /// element for the current node positions.
    pub fn elastic_energy(&self, nodes: &[ElectroacitiveFemNode]) -> f64 {
        let f = self.deformation_gradient(nodes);
        let j = determinant(&f).max(MIN_JACOBIAN);
        let ln_j = j.ln();
        let trace_ftf: f64 = f.iter().flat_map(|row| row.iter()).map(|x| x * x).sum();
        let mu = self.shear_modulus;
        self.volume_0
            * (0.5 * mu * (trace_ftf - 3.0) - mu * ln_j + 0.5 * self.lame_lambda() * ln_j * ln_j)
    }

    /// Elastic nodal forces `f_k = −V₀ P ∇N_k` of the neo-Hookean law.  They are the
    /// exact negative gradient of [`elastic_energy`](Self::elastic_energy), sum to
    /// zero and vanish for rigid motions and in the rest configuration.
    pub fn elastic_nodal_forces(&self, nodes: &[ElectroacitiveFemNode]) -> [[f64; 3]; 4] {
        let p = self.first_piola_stress(self.deformation_gradient(nodes));
        let grads = self.shape_gradients();
        let mut forces = [[0.0_f64; 3]; 4];
        for (force, grad) in forces.iter_mut().zip(grads.iter()) {
            let p_grad = mat_vec(&p, *grad);
            for d in 0..3 {
                force[d] = -self.volume_0 * p_grad[d];
            }
        }
        forces
    }

    /// Electric field `E = −Σ_k φ_k ∇N_k` and current volume of the element and the
    /// current shape-function gradients, or `None` if the element is degenerate.
    fn current_field(
        &self,
        nodes: &[ElectroacitiveFemNode],
    ) -> Option<([f64; 3], f64, [[f64; 3]; 4])> {
        let corners = [
            nodes[self.node_indices[0]].position,
            nodes[self.node_indices[1]].position,
            nodes[self.node_indices[2]].position,
            nodes[self.node_indices[3]].position,
        ];
        let (volume, grads) = tet_geometry(corners, 1.0e-9 * self.volume_0)?;
        let mut field = [0.0_f64; 3];
        for (k, grad) in grads.iter().enumerate() {
            let potential = nodes[self.node_indices[k]].potential;
            for d in 0..3 {
                field[d] -= potential * grad[d];
            }
        }
        Some((field, volume, grads))
    }

    /// Maxwell nodal forces `f_k = −v T ∇N_k` evaluated on the **current**
    /// geometry: field, volume and gradients follow the deformed element.  At the
    /// rest configuration this is [`maxwell_nodal_forces`](Self::maxwell_nodal_forces);
    /// under deformation it captures that a film that thins at fixed potentials sees
    /// a stronger field.  At fixed nodal potentials the forces are the exact
    /// derivative `+∂/∂x_k` of [`electric_energy`](Self::electric_energy).  A
    /// degenerate (flat) element carries no force.
    pub fn maxwell_nodal_forces_deformed(&self, nodes: &[ElectroacitiveFemNode]) -> [[f64; 3]; 4] {
        let Some((field, volume, grads)) = self.current_field(nodes) else {
            return [[0.0; 3]; 4];
        };
        let eps = EPSILON_0 * self.epsilon_r;
        let e2 = v3_dot(field, field);
        let mut stress = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                stress[i][j] = eps * (field[i] * field[j] - 0.5 * e2 * delta);
            }
        }
        let mut forces = [[0.0_f64; 3]; 4];
        for (force, grad) in forces.iter_mut().zip(grads.iter()) {
            let t_grad = mat_vec(&stress, *grad);
            for d in 0..3 {
                force[d] = -volume * t_grad[d];
            }
        }
        forces
    }

    /// Stored electric energy `½ ε v |E|²` \[J\] of the element at the current
    /// geometry and nodal potentials (zero for a degenerate element).
    pub fn electric_energy(&self, nodes: &[ElectroacitiveFemNode]) -> f64 {
        match self.current_field(nodes) {
            Some((field, volume, _)) => {
                0.5 * EPSILON_0 * self.epsilon_r * volume * v3_dot(field, field)
            }
            None => 0.0,
        }
    }

    /// Total nodal forces of the element: neo-Hookean elastic forces plus the
    /// Maxwell forces of the deformed geometry.
    pub fn nodal_forces(&self, nodes: &[ElectroacitiveFemNode]) -> [[f64; 3]; 4] {
        let elastic = self.elastic_nodal_forces(nodes);
        let maxwell = self.maxwell_nodal_forces_deformed(nodes);
        let mut total = [[0.0_f64; 3]; 4];
        for k in 0..4 {
            for d in 0..3 {
                total[k][d] = elastic[k][d] + maxwell[k][d];
            }
        }
        total
    }
}

impl ElectroacitiveFemBody {
    /// Add the per-element nodal forces given by `element_forces` to `forces`
    /// (one entry per node; entries beyond `forces.len()` are skipped).
    fn accumulate_forces(
        &self,
        forces: &mut [[f64; 3]],
        element_forces: impl Fn(&ElectroacitiveTetraElement) -> [[f64; 3]; 4],
    ) {
        for element in &self.elements {
            let nodal = element_forces(element);
            for (k, &node) in element.node_indices.iter().enumerate() {
                if let Some(slot) = forces.get_mut(node) {
                    for d in 0..3 {
                        slot[d] += nodal[k][d];
                    }
                }
            }
        }
    }

    /// Add the elastic (neo-Hookean) nodal forces \[N\] of all elements to `forces`
    /// (one entry per node).  Pure accumulation: nodes are neither accelerated nor
    /// moved, [`step`](Self::step) integrates the accumulated forces.
    pub fn apply_elastic_forces(&self, forces: &mut [[f64; 3]]) {
        self.accumulate_forces(forces, |element| element.elastic_nodal_forces(&self.nodes));
    }

    /// Add the Maxwell nodal forces \[N\] of all elements to `forces` (one entry per
    /// node), evaluated on the current geometry (see
    /// [`ElectroacitiveTetraElement::maxwell_nodal_forces_deformed`]).  Pure
    /// accumulation: nodes are neither accelerated nor moved, [`step`](Self::step)
    /// integrates the accumulated forces with the time step.
    pub fn apply_maxwell_forces(&self, forces: &mut [[f64; 3]]) {
        self.accumulate_forces(forces, |element| {
            element.maxwell_nodal_forces_deformed(&self.nodes)
        });
    }

    /// Elastic (neo-Hookean) force \[N\] on every node.
    pub fn elastic_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0_f64; 3]; self.nodes.len()];
        self.apply_elastic_forces(&mut forces);
        forces
    }

    /// Maxwell force \[N\] of the electric field on every node.  A constant field
    /// exerts no net force: the forces of the body sum to zero.
    pub fn maxwell_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0_f64; 3]; self.nodes.len()];
        self.apply_maxwell_forces(&mut forces);
        forces
    }

    /// Elastic plus Maxwell force \[N\] on every node (gravity and damping are
    /// added by [`step`](Self::step)).
    pub fn net_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0_f64; 3]; self.nodes.len()];
        self.apply_elastic_forces(&mut forces);
        self.apply_maxwell_forces(&mut forces);
        forces
    }

    /// Advance the body by `dt` with semi-implicit Euler.
    ///
    /// The elastic and Maxwell forces are accumulated, the velocity of every free
    /// node changes by `(F / m + g) dt` and is damped by the mass-proportional
    /// Rayleigh coefficient, `v ← v (1 − α dt)`, and the positions advance with the
    /// new velocity.  Pinned nodes do not move.  The explicit scheme needs a time
    /// step well below the transit time of the elastic wave through the smallest
    /// element.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let forces = self.net_forces();
        let damping = 1.0 - self.damping * dt;
        for (node, force) in self.nodes.iter_mut().zip(forces.iter()) {
            if node.pinned {
                continue;
            }
            let inv_mass = 1.0 / node.mass.max(1.0e-30);
            for d in 0..3 {
                let acceleration = force[d] * inv_mass + gravity[d];
                node.velocity[d] = (node.velocity[d] + acceleration * dt) * damping;
            }
        }
        for node in &mut self.nodes {
            if node.pinned {
                continue;
            }
            for d in 0..3 {
                node.position[d] += node.velocity[d] * dt;
            }
        }
    }

    /// Total strain energy of the neo-Hookean elements \[J\].
    pub fn elastic_energy(&self) -> f64 {
        self.elements
            .iter()
            .map(|element| element.elastic_energy(&self.nodes))
            .sum()
    }

    /// Total stored electric energy `Σ ½ ε v |E|²` \[J\] at the current geometry.
    pub fn electric_energy(&self) -> f64 {
        self.elements
            .iter()
            .map(|element| element.electric_energy(&self.nodes))
            .sum()
    }

    /// Kinetic energy of the free nodes \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        self.nodes
            .iter()
            .filter(|node| !node.pinned)
            .map(|node| 0.5 * node.mass * v3_dot(node.velocity, node.velocity))
            .sum()
    }
}
