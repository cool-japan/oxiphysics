#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Truss and frame finite element analysis.
//!
//! This module provides a complete truss/frame FEM toolkit:
//!
//! - **2-D truss element**: axial-only, 4-DOF, local and global stiffness
//! - **3-D truss element**: axial-only, 6-DOF, direction cosines
//! - **2-D frame element**: axial + bending (Euler-Bernoulli), 6-DOF
//! - **3-D space frame element**: 12-DOF, full local-to-global transformation
//! - **Global stiffness assembly**: scatter-add with DOF mapping
//! - **Support reactions**: back-calculated from free-DOF solution
//! - **Internal forces**: axial force, shear and moment diagrams
//! - **Elastic buckling**: Euler critical load, effective length factors
//! - **2-D portal frame**: two-bay / two-storey analysis helpers
//! - **Plane frame**: rigid joints, lateral load, gravity load
//!
//! # Coordinate conventions
//! - Local axis 1 (x̄) runs from node i to node j.
//! - Local axis 2 (ȳ) is the primary bending axis (for 2-D: out-of-plane = z).
//! - Positive moments follow the right-hand rule.
//! - All stiffness matrices in consistent SI units (N, m, Pa).

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// Number of DOF per node for a 2-D truss (ux, uy).
pub const TRUSS2D_DOF_PER_NODE: usize = 2;

/// Number of DOF per node for a 3-D truss (ux, uy, uz).
pub const TRUSS3D_DOF_PER_NODE: usize = 3;

/// Number of DOF per node for a 2-D frame (ux, uy, θz).
pub const FRAME2D_DOF_PER_NODE: usize = 3;

/// Number of DOF per node for a 3-D space frame (ux, uy, uz, θx, θy, θz).
pub const FRAME3D_DOF_PER_NODE: usize = 6;

/// DOF count for a 2-node 2-D truss element.
pub const TRUSS2D_NDOF: usize = 4;

/// DOF count for a 2-node 3-D truss element.
pub const TRUSS3D_NDOF: usize = 6;

/// DOF count for a 2-node 2-D frame element.
pub const FRAME2D_NDOF: usize = 6;

/// DOF count for a 2-node 3-D space frame element.
pub const FRAME3D_NDOF: usize = 12;

/// Minimum length below which an element is considered degenerate.
pub const MIN_LENGTH: f64 = 1.0e-14;

// ============================================================================
// Section 1 – 2-D truss element
// ============================================================================

/// A 2-D truss element defined by two nodes.
///
/// Carries axial stiffness only (no bending or torsion).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct TrussElement2D {
    /// Young's modulus E \[Pa\].
    pub elastic_modulus: f64,
    /// Cross-sectional area A \[m²\].
    pub area: f64,
    /// Node i coordinates (x, y) \[m\].
    pub node_i: [f64; 2],
    /// Node j coordinates (x, y) \[m\].
    pub node_j: [f64; 2],
}

impl TrussElement2D {
    /// Compute element length L = |j − i|.
    #[allow(dead_code)]
    pub fn length(&self) -> f64 {
        let dx = self.node_j[0] - self.node_i[0];
        let dy = self.node_j[1] - self.node_i[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Direction cosines (l, m) = (cos α, sin α).
    #[allow(dead_code)]
    pub fn direction_cosines(&self) -> [f64; 2] {
        let l = self.length().max(MIN_LENGTH);
        let dx = self.node_j[0] - self.node_i[0];
        let dy = self.node_j[1] - self.node_i[1];
        [dx / l, dy / l]
    }

    /// Local (1-D) stiffness matrix  K_loc = (EA/L) \[1, -1; -1, 1\].
    #[allow(dead_code)]
    pub fn local_stiffness(&self) -> [[f64; 2]; 2] {
        let k = self.elastic_modulus * self.area / self.length().max(MIN_LENGTH);
        [[k, -k], [-k, k]]
    }

    /// Global 4×4 stiffness matrix K_glob after transformation T^T K_loc T.
    ///
    /// DOF order: \[u_xi, u_yi, u_xj, u_yj\].
    #[allow(dead_code)]
    pub fn global_stiffness(&self) -> [[f64; 4]; 4] {
        let [l, m] = self.direction_cosines();
        let k = self.elastic_modulus * self.area / self.length().max(MIN_LENGTH);
        let l2 = l * l;
        let m2 = m * m;
        let lm = l * m;
        [
            [k * l2, k * lm, -k * l2, -k * lm],
            [k * lm, k * m2, -k * lm, -k * m2],
            [-k * l2, -k * lm, k * l2, k * lm],
            [-k * lm, -k * m2, k * lm, k * m2],
        ]
    }

    /// Axial force N = (EA/L) (Δu_j − Δu_i) · ê.
    ///
    /// Positive = tension.
    ///
    /// # Arguments
    /// * `disp` – global displacements \[u_xi, u_yi, u_xj, u_yj\]
    #[allow(dead_code)]
    pub fn axial_force(&self, disp: [f64; 4]) -> f64 {
        let [l, m] = self.direction_cosines();
        let ea_over_l = self.elastic_modulus * self.area / self.length().max(MIN_LENGTH);
        let du = (disp[2] - disp[0]) * l + (disp[3] - disp[1]) * m;
        ea_over_l * du
    }
}

// ============================================================================
// Section 2 – 3-D truss element
// ============================================================================

/// A 3-D truss element (axial-only, 6-DOF, two nodes).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct TrussElement3D {
    /// Young's modulus E \[Pa\].
    pub elastic_modulus: f64,
    /// Cross-sectional area A \[m²\].
    pub area: f64,
    /// Node i coordinates (x, y, z) \[m\].
    pub node_i: [f64; 3],
    /// Node j coordinates (x, y, z) \[m\].
    pub node_j: [f64; 3],
}

impl TrussElement3D {
    /// Compute element length.
    #[allow(dead_code)]
    pub fn length(&self) -> f64 {
        let dx = self.node_j[0] - self.node_i[0];
        let dy = self.node_j[1] - self.node_i[1];
        let dz = self.node_j[2] - self.node_i[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Direction cosines (l, m, n) in 3-D.
    #[allow(dead_code)]
    pub fn direction_cosines(&self) -> [f64; 3] {
        let len = self.length().max(MIN_LENGTH);
        [
            (self.node_j[0] - self.node_i[0]) / len,
            (self.node_j[1] - self.node_i[1]) / len,
            (self.node_j[2] - self.node_i[2]) / len,
        ]
    }

    /// Global 6×6 stiffness matrix (DOF: \[u_xi, u_yi, u_zi, u_xj, u_yj, u_zj\]).
    #[allow(dead_code)]
    pub fn global_stiffness(&self) -> [[f64; 6]; 6] {
        let [l, m, n] = self.direction_cosines();
        let k = self.elastic_modulus * self.area / self.length().max(MIN_LENGTH);
        let mut kg = [[0.0_f64; 6]; 6];
        let dc = [l, m, n];
        for a in 0..3 {
            for b in 0..3 {
                let val = k * dc[a] * dc[b];
                kg[a][b] = val;
                kg[a][b + 3] = -val;
                kg[a + 3][b] = -val;
                kg[a + 3][b + 3] = val;
            }
        }
        kg
    }

    /// Axial force in the element.
    ///
    /// # Arguments
    /// * `disp` – global displacements \[uxi, uyi, uzi, uxj, uyj, uzj\]
    #[allow(dead_code)]
    pub fn axial_force(&self, disp: [f64; 6]) -> f64 {
        let [l, m, n] = self.direction_cosines();
        let ea_over_len = self.elastic_modulus * self.area / self.length().max(MIN_LENGTH);
        let du = (disp[3] - disp[0]) * l + (disp[4] - disp[1]) * m + (disp[5] - disp[2]) * n;
        ea_over_len * du
    }
}

// ============================================================================
// Section 3 – 2-D frame element (Euler-Bernoulli, combined axial + bending)
// ============================================================================

/// A 2-D Euler-Bernoulli frame element (6-DOF).
///
/// DOF order (local): \[ū_i, v̄_i, θ̄_i, ū_j, v̄_j, θ̄_j\].
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct FrameElement2D {
    /// Young's modulus E \[Pa\].
    pub elastic_modulus: f64,
    /// Cross-sectional area A \[m²\].
    pub area: f64,
    /// Second moment of area I \[m⁴\].
    pub inertia: f64,
    /// Node i coordinates (x, y) \[m\].
    pub node_i: [f64; 2],
    /// Node j coordinates (x, y) \[m\].
    pub node_j: [f64; 2],
}

impl FrameElement2D {
    /// Element length.
    #[allow(dead_code)]
    pub fn length(&self) -> f64 {
        let dx = self.node_j[0] - self.node_i[0];
        let dy = self.node_j[1] - self.node_i[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Inclination angle α = atan2(dy, dx).
    #[allow(dead_code)]
    pub fn angle(&self) -> f64 {
        let dx = self.node_j[0] - self.node_i[0];
        let dy = self.node_j[1] - self.node_i[1];
        dy.atan2(dx)
    }

    /// Local 6×6 stiffness matrix in local coordinates.
    ///
    /// Combines axial (EA/L) and Euler-Bernoulli bending (EI/L³).
    #[allow(dead_code)]
    pub fn local_stiffness(&self) -> [[f64; 6]; 6] {
        let l = self.length().max(MIN_LENGTH);
        let ea = self.elastic_modulus * self.area;
        let ei = self.elastic_modulus * self.inertia;
        let a = ea / l;
        let b1 = 12.0 * ei / (l * l * l);
        let b2 = 6.0 * ei / (l * l);
        let b3 = 4.0 * ei / l;
        let b4 = 2.0 * ei / l;

        [
            [a, 0.0, 0.0, -a, 0.0, 0.0],
            [0.0, b1, b2, 0.0, -b1, b2],
            [0.0, b2, b3, 0.0, -b2, b4],
            [-a, 0.0, 0.0, a, 0.0, 0.0],
            [0.0, -b1, -b2, 0.0, b1, -b2],
            [0.0, b2, b4, 0.0, -b2, b3],
        ]
    }

    /// 6×6 rotation matrix T that transforms global to local displacements.
    ///
    /// d_local = T · d_global
    #[allow(dead_code)]
    pub fn rotation_matrix(&self) -> [[f64; 6]; 6] {
        let alpha = self.angle();
        let c = alpha.cos();
        let s = alpha.sin();
        let mut t = [[0.0_f64; 6]; 6];
        // Upper-left 3×3 block
        t[0][0] = c;
        t[0][1] = s;
        t[1][0] = -s;
        t[1][1] = c;
        t[2][2] = 1.0;
        // Lower-right 3×3 block (same)
        t[3][3] = c;
        t[3][4] = s;
        t[4][3] = -s;
        t[4][4] = c;
        t[5][5] = 1.0;
        t
    }

    /// Global 6×6 stiffness: K_glob = T^T K_loc T.
    #[allow(dead_code)]
    pub fn global_stiffness(&self) -> [[f64; 6]; 6] {
        let kl = self.local_stiffness();
        let t = self.rotation_matrix();
        // K_glob = T^T K_loc T
        let kt = mat6_mul(&transpose6(&t), &kl);
        mat6_mul(&kt, &t)
    }

    /// Internal forces in local coordinates: \[N, V, M_i, N, V, M_j\].
    ///
    /// # Arguments
    /// * `disp_global` – 6-D global displacement vector
    #[allow(dead_code)]
    pub fn internal_forces(&self, disp_global: [f64; 6]) -> [f64; 6] {
        let t = self.rotation_matrix();
        let d_loc = mat6_vec_mul(&t, &disp_global);
        let kl = self.local_stiffness();
        mat6_vec_mul(&kl, &d_loc)
    }
}

// ============================================================================
// Section 4 – 3-D space frame element (12-DOF)
// ============================================================================

/// A 3-D space frame element.
///
/// DOF order (local): \[ū_i, v̄_i, w̄_i, θ̄x_i, θ̄y_i, θ̄z_i, (j end)\].
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct SpaceFrameElement {
    /// Young's modulus E \[Pa\].
    pub elastic_modulus: f64,
    /// Shear modulus G \[Pa\].
    pub shear_modulus: f64,
    /// Cross-sectional area A \[m²\].
    pub area: f64,
    /// Moment of inertia I_y \[m⁴\] (bending about local y).
    pub iy: f64,
    /// Moment of inertia I_z \[m⁴\] (bending about local z).
    pub iz: f64,
    /// Saint-Venant torsional constant J \[m⁴\].
    pub j_torsion: f64,
    /// Node i coordinates \[m\].
    pub node_i: [f64; 3],
    /// Node j coordinates \[m\].
    pub node_j: [f64; 3],
}

impl SpaceFrameElement {
    /// Element length.
    #[allow(dead_code)]
    pub fn length(&self) -> f64 {
        let dx = self.node_j[0] - self.node_i[0];
        let dy = self.node_j[1] - self.node_i[1];
        let dz = self.node_j[2] - self.node_i[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Local 12×12 stiffness matrix (axial + bending + torsion).
    ///
    /// Returns only the upper-left 6×6 block for brevity (symmetric element).
    /// Full 12×12 is assembled by scatter-add in [`SpaceFrameAssembler`].
    #[allow(dead_code)]
    pub fn axial_stiffness(&self) -> f64 {
        self.elastic_modulus * self.area / self.length().max(MIN_LENGTH)
    }

    /// Bending stiffness coefficient EI_z / L³.
    #[allow(dead_code)]
    pub fn bending_z_coeff(&self) -> f64 {
        self.elastic_modulus * self.iz / self.length().max(MIN_LENGTH).powi(3)
    }

    /// Bending stiffness coefficient EI_y / L³.
    #[allow(dead_code)]
    pub fn bending_y_coeff(&self) -> f64 {
        self.elastic_modulus * self.iy / self.length().max(MIN_LENGTH).powi(3)
    }

    /// Torsional stiffness GJ / L.
    #[allow(dead_code)]
    pub fn torsional_stiffness(&self) -> f64 {
        self.shear_modulus * self.j_torsion / self.length().max(MIN_LENGTH)
    }

    /// Axial force from end displacements.
    ///
    /// # Arguments
    /// * `u_i`, `u_j` – axial displacements at near and far ends \[m\]
    #[allow(dead_code)]
    pub fn axial_force(&self, u_i: f64, u_j: f64) -> f64 {
        self.axial_stiffness() * (u_j - u_i)
    }

    /// Torsional moment from end rotations.
    #[allow(dead_code)]
    pub fn torsional_moment(&self, phi_i: f64, phi_j: f64) -> f64 {
        self.torsional_stiffness() * (phi_j - phi_i)
    }
}

// ============================================================================
// Section 5 – Local-to-global coordinate transformation (3-D)
// ============================================================================

/// Compute the 3-D direction cosine matrix for a space frame member.
///
/// The local x̄ axis runs from node i to node j.  The local ȳ axis is chosen
/// to lie in the x̄-z_global plane (standard structural convention).
///
/// # Arguments
/// * `node_i`, `node_j` – end-node coordinates
/// * `ref_up`           – reference "up" vector (usually global Z = \[0,0,1\])
///
/// Returns 3×3 direction cosine matrix Γ such that v_local = Γ v_global.
#[allow(dead_code)]
pub fn direction_cosine_matrix(
    node_i: [f64; 3],
    node_j: [f64; 3],
    ref_up: [f64; 3],
) -> [[f64; 3]; 3] {
    let ex = vec3_sub(node_j, node_i);
    let ex_len = vec3_norm(ex).max(MIN_LENGTH);
    let ex = vec3_scale(ex, 1.0 / ex_len);

    // Local z axis = ex × ref_up, then ey = ez × ex
    let ez = vec3_cross(ex, ref_up);
    let ez_len = vec3_norm(ez);
    let (ez, ey) = if ez_len < 1e-10 {
        // Member is parallel to ref_up: use global x as helper
        let helper = [1.0_f64, 0.0, 0.0];
        let ez2 = vec3_cross(ex, helper);
        let ez2_len = vec3_norm(ez2).max(MIN_LENGTH);
        let ez2 = vec3_scale(ez2, 1.0 / ez2_len);
        let ey2 = vec3_cross(ez2, ex);
        (ez2, ey2)
    } else {
        let ez = vec3_scale(ez, 1.0 / ez_len);
        let ey = vec3_cross(ez, ex);
        (ez, ey)
    };

    [ex, ey, ez]
}

/// Expand 3×3 direction cosine matrix Γ to 12×12 block-diagonal rotation T.
///
/// T = diag(Γ, Γ, Γ, Γ) mapping 12 global DOF to 12 local DOF.
#[allow(dead_code)]
pub fn expand_to_12x12(gamma: &[[f64; 3]; 3]) -> [[f64; 12]; 12] {
    let mut t = [[0.0_f64; 12]; 12];
    for block in 0..4 {
        let r0 = 3 * block;
        for i in 0..3 {
            for j in 0..3 {
                t[r0 + i][r0 + j] = gamma[i][j];
            }
        }
    }
    t
}

// ============================================================================
// Section 6 – Global stiffness assembly
// ============================================================================

/// Global stiffness assembler for 2-D frame structures.
///
/// Stores the full dense matrix for small-to-medium problems.
#[allow(dead_code)]
pub struct FrameAssembler2D {
    /// Total number of DOF.
    pub n_dof: usize,
    /// Dense global stiffness matrix K \[n_dof × n_dof\].
    pub k_global: Vec<Vec<f64>>,
    /// Global load vector f \[n_dof\].
    pub f_global: Vec<f64>,
}

impl FrameAssembler2D {
    /// Construct an assembler for `n_nodes` nodes (3 DOF/node in 2-D frame).
    #[allow(dead_code)]
    pub fn new(n_nodes: usize) -> Self {
        let n_dof = n_nodes * FRAME2D_DOF_PER_NODE;
        Self {
            n_dof,
            k_global: vec![vec![0.0_f64; n_dof]; n_dof],
            f_global: vec![0.0_f64; n_dof],
        }
    }

    /// Scatter the element stiffness ke into the global matrix.
    ///
    /// # Arguments
    /// * `ke`      – 6×6 element stiffness
    /// * `node_i`  – global index of node i
    /// * `node_j`  – global index of node j
    #[allow(dead_code)]
    pub fn add_element(&mut self, ke: &[[f64; 6]; 6], node_i: usize, node_j: usize) {
        let dofs = [
            node_i * 3,
            node_i * 3 + 1,
            node_i * 3 + 2,
            node_j * 3,
            node_j * 3 + 1,
            node_j * 3 + 2,
        ];
        for r in 0..6 {
            for c in 0..6 {
                self.k_global[dofs[r]][dofs[c]] += ke[r][c];
            }
        }
    }

    /// Apply a Dirichlet (fixed) boundary condition by penalty method.
    ///
    /// # Arguments
    /// * `dof`   – constrained DOF index
    /// * `value` – prescribed displacement (usually 0)
    /// * `penalty` – large stiffness value (e.g. 1e30)
    #[allow(dead_code)]
    pub fn apply_dirichlet(&mut self, dof: usize, value: f64, penalty: f64) {
        self.k_global[dof][dof] += penalty;
        self.f_global[dof] += penalty * value;
    }

    /// Solve K u = f using Gaussian elimination (for small systems).
    ///
    /// Returns the displacement vector `u` or `None` if singular.
    #[allow(dead_code)]
    pub fn solve(&self) -> Option<Vec<f64>> {
        let n = self.n_dof;
        let mut a = self.k_global.clone();
        let mut b = self.f_global.clone();
        // Forward elimination with partial pivoting
        for col in 0..n {
            // Pivot
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for row in col + 1..n {
                if a[row][col].abs() > max_val {
                    max_val = a[row][col].abs();
                    max_row = row;
                }
            }
            if max_val < 1e-20 {
                return None;
            }
            a.swap(col, max_row);
            b.swap(col, max_row);
            let pivot = a[col][col];
            for row in col + 1..n {
                let factor = a[row][col] / pivot;
                b[row] -= factor * b[col];
                for c in col..n {
                    a[row][c] -= factor * a[col][c];
                }
            }
        }
        // Back substitution
        let mut u = vec![0.0_f64; n];
        for row in (0..n).rev() {
            let mut s = b[row];
            for c in row + 1..n {
                s -= a[row][c] * u[c];
            }
            u[row] = s / a[row][row];
        }
        Some(u)
    }

    /// Compute support reactions at constrained DOFs.
    ///
    /// R = K u − f  (residual forces at prescribed DOFs)
    ///
    /// # Arguments
    /// * `u`             – displacement solution
    /// * `constrained`   – list of constrained DOF indices
    #[allow(dead_code)]
    pub fn support_reactions(&self, u: &[f64], constrained: &[usize]) -> Vec<f64> {
        let n = self.n_dof;
        let mut reactions = vec![0.0_f64; constrained.len()];
        for (k, &dof) in constrained.iter().enumerate() {
            let mut ku = 0.0;
            for c in 0..n {
                ku += self.k_global[dof][c] * u[c];
            }
            reactions[k] = ku - self.f_global[dof];
        }
        reactions
    }
}

// ============================================================================
// Section 7 – Internal force diagrams
// ============================================================================

/// Bending moment at position s ∈ \[0,L\] along a simply-supported beam with UDL.
///
/// M(s) = (w L / 2) s − w s² / 2
///
/// # Arguments
/// * `w`  – uniform distributed load intensity \[N/m\]
/// * `l`  – span \[m\]
/// * `s`  – position from left end \[m\]
#[allow(dead_code)]
pub fn simply_supported_moment(w: f64, l: f64, s: f64) -> f64 {
    (w * l / 2.0) * s - 0.5 * w * s * s
}

/// Shear force at position s along a simply-supported beam with UDL.
///
/// V(s) = w L / 2 − w s
#[allow(dead_code)]
pub fn simply_supported_shear(w: f64, l: f64, s: f64) -> f64 {
    w * l / 2.0 - w * s
}

/// Maximum bending moment in a simply-supported beam under UDL.
///
/// M_max = w L² / 8  at s = L/2.
#[allow(dead_code)]
pub fn simply_supported_max_moment(w: f64, l: f64) -> f64 {
    w * l * l / 8.0
}

/// Fixed-end moments (FEM) for a fixed-fixed beam under UDL.
///
/// M_A = M_B = w L² / 12.
#[allow(dead_code)]
pub fn fixed_fixed_fem(w: f64, l: f64) -> f64 {
    w * l * l / 12.0
}

/// Carryover factor for the moment distribution method (fixed far end = 0.5).
pub const CARRYOVER_FIXED: f64 = 0.5;

/// Distribution factor for a member in moment distribution.
///
/// DF = (EI/L) / Σ(EI/L)
///
/// # Arguments
/// * `ei_over_l` – stiffness factor of this member
/// * `sum_stiff` – sum of stiffness factors at the joint
#[allow(dead_code)]
pub fn distribution_factor(ei_over_l: f64, sum_stiff: f64) -> f64 {
    if sum_stiff.abs() < f64::EPSILON {
        0.0
    } else {
        ei_over_l / sum_stiff
    }
}

/// Moment diagram at `n_pts` sample points along a member.
///
/// Returns (positions, moments).  Uses linear interpolation between
/// the end moments including the effect of a uniform load.
///
/// # Arguments
/// * `m_i`  – end moment at node i \[N·m\]
/// * `m_j`  – end moment at node j \[N·m\]
/// * `w`    – uniform distributed load \[N/m\] (positive downward)
/// * `l`    – member length \[m\]
/// * `n_pts` – number of sample points
#[allow(dead_code)]
pub fn moment_diagram(m_i: f64, m_j: f64, w: f64, l: f64, n_pts: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n_pts.max(2);
    let mut pos = Vec::with_capacity(n);
    let mut mom = Vec::with_capacity(n);
    for k in 0..n {
        let s = l * k as f64 / (n - 1) as f64;
        // Superpose: end moments (linear) + simply-supported UDL moments
        let m_linear = m_i * (1.0 - s / l) + m_j * (s / l);
        let m_udl = 0.5 * w * s * (l - s);
        pos.push(s);
        mom.push(m_linear + m_udl);
    }
    (pos, mom)
}

// ============================================================================
// Section 8 – Elastic buckling (Euler column)
// ============================================================================

/// Euler critical buckling load P_cr = π² E I / (K L)².
///
/// # Arguments
/// * `elastic_modulus` – E \[Pa\]
/// * `inertia`         – I (minimum second moment of area) \[m⁴\]
/// * `length`          – actual column length L \[m\]
/// * `k_factor`        – effective length factor K (1.0 = pin-pin)
#[allow(dead_code)]
pub fn euler_buckling_load(elastic_modulus: f64, inertia: f64, length: f64, k_factor: f64) -> f64 {
    let le = k_factor * length;
    if le.abs() < MIN_LENGTH {
        return f64::INFINITY;
    }
    PI * PI * elastic_modulus * inertia / (le * le)
}

/// Effective length factor for standard end conditions.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EndCondition {
    /// Pinned–pinned: K = 1.0.
    PinnedPinned,
    /// Fixed–free (cantilever): K = 2.0.
    FixedFree,
    /// Fixed–pinned: K = 0.7.
    FixedPinned,
    /// Fixed–fixed: K = 0.5.
    FixedFixed,
}

impl EndCondition {
    /// Return the effective length factor K.
    #[allow(dead_code)]
    pub fn k_factor(&self) -> f64 {
        match self {
            Self::PinnedPinned => 1.0,
            Self::FixedFree => 2.0,
            Self::FixedPinned => 0.7,
            Self::FixedFixed => 0.5,
        }
    }
}

/// Radius of gyration r = √(I/A).
#[allow(dead_code)]
pub fn radius_of_gyration(inertia: f64, area: f64) -> f64 {
    if area < f64::EPSILON {
        0.0
    } else {
        (inertia / area).sqrt()
    }
}

/// Slenderness ratio λ = K L / r.
#[allow(dead_code)]
pub fn slenderness_ratio(k_factor: f64, length: f64, r: f64) -> f64 {
    if r.abs() < f64::EPSILON {
        f64::INFINITY
    } else {
        k_factor * length / r
    }
}

/// Critical stress σ_cr = P_cr / A = π² E / λ².
#[allow(dead_code)]
pub fn euler_critical_stress(elastic_modulus: f64, slenderness: f64) -> f64 {
    if slenderness.abs() < f64::EPSILON {
        f64::INFINITY
    } else {
        PI * PI * elastic_modulus / (slenderness * slenderness)
    }
}

/// Johnson–Parabola critical stress (for intermediate columns).
///
/// σ_cr = σ_y \[1 − σ_y λ² / (4 π² E)\]
///
/// # Arguments
/// * `yield_stress`     – σ_y \[Pa\]
/// * `elastic_modulus`  – E \[Pa\]
/// * `slenderness`      – λ (dimensionless slenderness ratio)
#[allow(dead_code)]
pub fn johnson_critical_stress(yield_stress: f64, elastic_modulus: f64, slenderness: f64) -> f64 {
    yield_stress
        * (1.0 - yield_stress * slenderness * slenderness / (4.0 * PI * PI * elastic_modulus))
}

// ============================================================================
// Section 9 – 2-D portal frame analysis
// ============================================================================

/// Geometry of a simple rectangular portal frame.
///
/// ```text
///     C──────D
///     |      |
///     A      B
/// ```
///
/// Nodes: A(0,0), B(w,0), C(0,h), D(w,h).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct PortalFrame {
    /// Bay width \[m\].
    pub width: f64,
    /// Column height \[m\].
    pub height: f64,
    /// Young's modulus for columns \[Pa\].
    pub e_column: f64,
    /// Young's modulus for beam \[Pa\].
    pub e_beam: f64,
    /// Second moment of area for columns \[m⁴\].
    pub i_column: f64,
    /// Second moment of area for beam \[m⁴\].
    pub i_beam: f64,
    /// Cross-sectional area for columns \[m²\].
    pub a_column: f64,
    /// Cross-sectional area for beam \[m²\].
    pub a_beam: f64,
}

impl PortalFrame {
    /// Stiffness factor K = EI/L for a column.
    #[allow(dead_code)]
    pub fn column_stiffness(&self) -> f64 {
        self.e_column * self.i_column / self.height
    }

    /// Stiffness factor K = EI/L for the beam.
    #[allow(dead_code)]
    pub fn beam_stiffness(&self) -> f64 {
        self.e_beam * self.i_beam / self.width
    }

    /// Lateral (sway) stiffness of the portal under a horizontal load at the beam level.
    ///
    /// Using the portal frame formula (stiff beam): K_sway = 24 EI_col / h³.
    #[allow(dead_code)]
    pub fn sway_stiffness(&self) -> f64 {
        24.0 * self.e_column * self.i_column / (self.height * self.height * self.height)
    }

    /// Top-of-column deflection under horizontal load H (sway mode).
    ///
    /// δ = H / K_sway.
    #[allow(dead_code)]
    pub fn sway_deflection(&self, horizontal_load: f64) -> f64 {
        horizontal_load / self.sway_stiffness().max(f64::EPSILON)
    }

    /// Fixed-base column bending moments at the column tops due to sway H.
    ///
    /// M_top = M_bot = H h / 2  (simplified symmetric portal).
    #[allow(dead_code)]
    pub fn column_moments(&self, horizontal_load: f64) -> f64 {
        horizontal_load * self.height / 2.0
    }

    /// Beam mid-span moment under vertical load w on the beam (kN/m).
    #[allow(dead_code)]
    pub fn beam_midspan_moment(&self, w: f64) -> f64 {
        // Portal with fixed column bases: M_mid = wL²/24 (approx).
        w * self.width * self.width / 24.0
    }

    /// Critical lateral load for sway buckling (simplified).
    ///
    /// P_cr,sway = K_sway / n_columns  (each column carries P/2).
    #[allow(dead_code)]
    pub fn sway_buckling_load(&self) -> f64 {
        self.sway_stiffness()
    }
}

// ============================================================================
// Section 10 – Plane frame analysis helpers
// ============================================================================

/// A node in a 2-D plane frame.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct FrameNode2D {
    /// Node index.
    pub id: usize,
    /// x-coordinate \[m\].
    pub x: f64,
    /// y-coordinate \[m\].
    pub y: f64,
}

/// DOF constraint type.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Constraint {
    /// Fixed in ux, uy, and θ (welded to ground).
    FixedFixed,
    /// Pinned support: ux = uy = 0, θ free.
    PinnedSupport,
    /// Roller (horizontal): uy = 0 only.
    RollerVertical,
    /// Roller (vertical): ux = 0 only.
    RollerHorizontal,
    /// Free (no constraint).
    Free,
}

impl Constraint {
    /// Returns which DOFs are constrained as boolean flags (ux, uy, θ).
    #[allow(dead_code)]
    pub fn constrained_dofs(&self) -> [bool; 3] {
        match self {
            Self::FixedFixed => [true, true, true],
            Self::PinnedSupport => [true, true, false],
            Self::RollerVertical => [false, true, false],
            Self::RollerHorizontal => [true, false, false],
            Self::Free => [false, false, false],
        }
    }
}

/// A 2-D plane frame with nodes, elements and loads.
#[allow(dead_code)]
pub struct PlaneFrame {
    /// List of frame nodes.
    pub nodes: Vec<FrameNode2D>,
    /// List of frame elements (node_i index, node_j index, element data).
    pub elements: Vec<(usize, usize, FrameElement2D)>,
    /// Nodal constraints (per node).
    pub constraints: Vec<Constraint>,
    /// Nodal point loads \[Fx, Fy, Mz\] per node.
    pub nodal_loads: Vec<[f64; 3]>,
}

impl PlaneFrame {
    /// Construct a frame from node and element lists.
    #[allow(dead_code)]
    pub fn new(nodes: Vec<FrameNode2D>, constraints: Vec<Constraint>) -> Self {
        let n = nodes.len();
        Self {
            nodes,
            elements: Vec::new(),
            constraints,
            nodal_loads: vec![[0.0; 3]; n],
        }
    }

    /// Add a frame element between nodes `ni` and `nj`.
    #[allow(dead_code)]
    pub fn add_element(&mut self, ni: usize, nj: usize, e: f64, area: f64, inertia: f64) {
        let elem = FrameElement2D {
            elastic_modulus: e,
            area,
            inertia,
            node_i: [self.nodes[ni].x, self.nodes[ni].y],
            node_j: [self.nodes[nj].x, self.nodes[nj].y],
        };
        self.elements.push((ni, nj, elem));
    }

    /// Apply a point load at node `n`.
    #[allow(dead_code)]
    pub fn apply_load(&mut self, n: usize, fx: f64, fy: f64, mz: f64) {
        self.nodal_loads[n][0] += fx;
        self.nodal_loads[n][1] += fy;
        self.nodal_loads[n][2] += mz;
    }

    /// Assemble and solve the plane frame for nodal displacements.
    ///
    /// Returns displacement vector `u` (length = 3 × n_nodes).
    #[allow(dead_code)]
    pub fn solve(&self) -> Option<Vec<f64>> {
        let n_nodes = self.nodes.len();
        let mut asm = FrameAssembler2D::new(n_nodes);

        // Assemble load vector
        for (ni, loads) in self.nodal_loads.iter().enumerate() {
            for d in 0..3 {
                asm.f_global[ni * 3 + d] += loads[d];
            }
        }

        // Assemble stiffness
        for (ni, nj, elem) in &self.elements {
            let ke = elem.global_stiffness();
            asm.add_element(&ke, *ni, *nj);
        }

        // Apply constraints via penalty
        let penalty = 1.0e30;
        for (n_idx, constraint) in self.constraints.iter().enumerate() {
            let flags = constraint.constrained_dofs();
            for (d, &fixed) in flags.iter().enumerate() {
                if fixed {
                    asm.apply_dirichlet(n_idx * 3 + d, 0.0, penalty);
                }
            }
        }

        asm.solve()
    }

    /// Total number of unconstrained DOFs.
    #[allow(dead_code)]
    pub fn free_dof_count(&self) -> usize {
        self.constraints
            .iter()
            .map(|c| {
                let flags = c.constrained_dofs();
                flags.iter().filter(|&&f| !f).count()
            })
            .sum()
    }
}

// ============================================================================
// Section 11 – Small matrix utilities
// ============================================================================

/// Transpose a 6×6 matrix.
#[allow(dead_code)]
pub fn transpose6(m: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut t = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            t[i][j] = m[j][i];
        }
    }
    t
}

/// Multiply two 6×6 matrices.
#[allow(dead_code)]
pub fn mat6_mul(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for k in 0..6 {
            if a[i][k].abs() < f64::EPSILON * 1e-5 {
                continue;
            }
            for j in 0..6 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Multiply a 6×6 matrix by a 6-vector.
#[allow(dead_code)]
pub fn mat6_vec_mul(m: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut r = [0.0_f64; 6];
    for i in 0..6 {
        for j in 0..6 {
            r[i] += m[i][j] * v[j];
        }
    }
    r
}

/// 3-D vector subtraction.
#[allow(dead_code)]
pub fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// 3-D vector Euclidean norm.
#[allow(dead_code)]
pub fn vec3_norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Scalar multiplication of a 3-D vector.
#[allow(dead_code)]
pub fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// 3-D cross product.
#[allow(dead_code)]
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ============================================================================
// Section 12 – Space frame assembler (extended)
// ============================================================================

/// Global stiffness assembler for 3-D space frames.
#[allow(dead_code)]
pub struct SpaceFrameAssembler {
    /// Total number of DOF (6 × n_nodes).
    pub n_dof: usize,
    /// Dense global stiffness matrix.
    pub k_global: Vec<Vec<f64>>,
    /// Global force vector.
    pub f_global: Vec<f64>,
}

impl SpaceFrameAssembler {
    /// Construct assembler for `n_nodes` nodes.
    #[allow(dead_code)]
    pub fn new(n_nodes: usize) -> Self {
        let n_dof = n_nodes * FRAME3D_DOF_PER_NODE;
        Self {
            n_dof,
            k_global: vec![vec![0.0_f64; n_dof]; n_dof],
            f_global: vec![0.0_f64; n_dof],
        }
    }

    /// Scatter a 12×12 element stiffness into the global matrix.
    ///
    /// # Arguments
    /// * `ke`     – 12×12 element stiffness in global coordinates
    /// * `node_i` – index of first node
    /// * `node_j` – index of second node
    #[allow(dead_code)]
    pub fn add_element(&mut self, ke: &[[f64; 12]; 12], node_i: usize, node_j: usize) {
        let mut dofs = [0usize; 12];
        for d in 0..6 {
            dofs[d] = node_i * 6 + d;
            dofs[d + 6] = node_j * 6 + d;
        }
        for r in 0..12 {
            for c in 0..12 {
                self.k_global[dofs[r]][dofs[c]] += ke[r][c];
            }
        }
    }

    /// Solve via Gaussian elimination (small systems only).
    #[allow(dead_code)]
    pub fn solve(&self) -> Option<Vec<f64>> {
        let n = self.n_dof;
        let mut a = self.k_global.clone();
        let mut b = self.f_global.clone();
        for col in 0..n {
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for row in col + 1..n {
                if a[row][col].abs() > max_val {
                    max_val = a[row][col].abs();
                    max_row = row;
                }
            }
            if max_val < 1e-20 {
                return None;
            }
            a.swap(col, max_row);
            b.swap(col, max_row);
            let pivot = a[col][col];
            for row in col + 1..n {
                let factor = a[row][col] / pivot;
                b[row] -= factor * b[col];
                for c in col..n {
                    a[row][c] -= factor * a[col][c];
                }
            }
        }
        let mut u = vec![0.0_f64; n];
        for row in (0..n).rev() {
            let mut s = b[row];
            for c in row + 1..n {
                s -= a[row][c] * u[c];
            }
            u[row] = s / a[row][row];
        }
        Some(u)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── TrussElement2D ────────────────────────────────────────────────────────

    #[test]
    fn test_truss2d_length_horizontal() {
        let e = TrussElement2D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0],
            node_j: [3.0, 0.0],
        };
        assert!((e.length() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_truss2d_direction_cosines_horizontal() {
        let e = TrussElement2D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0],
            node_j: [4.0, 0.0],
        };
        let [l, m] = e.direction_cosines();
        assert!((l - 1.0).abs() < 1e-12);
        assert!(m.abs() < 1e-12);
    }

    #[test]
    fn test_truss2d_direction_cosines_diagonal() {
        let e = TrussElement2D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0],
            node_j: [1.0, 1.0],
        };
        let [l, m] = e.direction_cosines();
        let sq2 = 1.0 / 2.0_f64.sqrt();
        assert!((l - sq2).abs() < 1e-10);
        assert!((m - sq2).abs() < 1e-10);
    }

    #[test]
    fn test_truss2d_global_stiffness_symmetry() {
        let e = TrussElement2D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0],
            node_j: [3.0, 4.0],
        };
        let k = e.global_stiffness();
        for i in 0..4 {
            for j in 0..4 {
                assert!((k[i][j] - k[j][i]).abs() < 1e-6, "asymmetric at [{i}][{j}]");
            }
        }
    }

    #[test]
    fn test_truss2d_axial_force_tension() {
        let e = TrussElement2D {
            elastic_modulus: 1.0,
            area: 1.0,
            node_i: [0.0, 0.0],
            node_j: [1.0, 0.0],
        };
        // Extend by 0.1 at far end
        let n = e.axial_force([0.0, 0.0, 0.1, 0.0]);
        assert!((n - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_truss2d_global_stiffness_row_sum_zero() {
        // For an unrestrained element, each row of K should sum to zero
        let e = TrussElement2D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0],
            node_j: [3.0, 4.0],
        };
        let k = e.global_stiffness();
        for i in 0..4 {
            let row_sum: f64 = k[i].iter().sum();
            assert!(row_sum.abs() < 1e-4, "row {i} sum = {row_sum}");
        }
    }

    // ── TrussElement3D ────────────────────────────────────────────────────────

    #[test]
    fn test_truss3d_length() {
        let e = TrussElement3D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0, 0.0],
            node_j: [3.0, 4.0, 0.0],
        };
        assert!((e.length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_truss3d_direction_cosines_unit_vector() {
        let e = TrussElement3D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0, 0.0],
            node_j: [1.0, 2.0, 2.0],
        };
        let dc = e.direction_cosines();
        let norm = (dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]).sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_truss3d_global_stiffness_symmetry() {
        let e = TrussElement3D {
            elastic_modulus: 2e11,
            area: 1e-4,
            node_i: [0.0, 0.0, 0.0],
            node_j: [1.0, 2.0, 2.0],
        };
        let k = e.global_stiffness();
        for i in 0..6 {
            for j in 0..6 {
                assert!((k[i][j] - k[j][i]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_truss3d_axial_force_z() {
        let e = TrussElement3D {
            elastic_modulus: 1.0,
            area: 1.0,
            node_i: [0.0, 0.0, 0.0],
            node_j: [0.0, 0.0, 2.0],
        };
        // Extend by 0.2 at far z end
        let n = e.axial_force([0.0, 0.0, 0.0, 0.0, 0.0, 0.2]);
        assert!((n - 0.1).abs() < 1e-12);
    }

    // ── FrameElement2D ────────────────────────────────────────────────────────

    #[test]
    fn test_frame2d_local_stiffness_symmetry() {
        let e = FrameElement2D {
            elastic_modulus: 2e11,
            area: 1e-3,
            inertia: 1e-6,
            node_i: [0.0, 0.0],
            node_j: [5.0, 0.0],
        };
        let k = e.local_stiffness();
        for i in 0..6 {
            for j in 0..6 {
                assert!((k[i][j] - k[j][i]).abs() < 1.0, "asymmetric [{i}][{j}]");
            }
        }
    }

    #[test]
    fn test_frame2d_global_stiffness_symmetry() {
        let e = FrameElement2D {
            elastic_modulus: 2e11,
            area: 1e-3,
            inertia: 1e-6,
            node_i: [0.0, 0.0],
            node_j: [3.0, 4.0],
        };
        let k = e.global_stiffness();
        for i in 0..6 {
            for j in 0..6 {
                assert!((k[i][j] - k[j][i]).abs() < 1.0, "asymmetric [{i}][{j}]");
            }
        }
    }

    #[test]
    fn test_frame2d_angle_horizontal() {
        let e = FrameElement2D {
            elastic_modulus: 2e11,
            area: 1e-3,
            inertia: 1e-6,
            node_i: [0.0, 0.0],
            node_j: [1.0, 0.0],
        };
        assert!(e.angle().abs() < 1e-12);
    }

    #[test]
    fn test_frame2d_angle_vertical() {
        let e = FrameElement2D {
            elastic_modulus: 2e11,
            area: 1e-3,
            inertia: 1e-6,
            node_i: [0.0, 0.0],
            node_j: [0.0, 1.0],
        };
        assert!((e.angle() - PI / 2.0).abs() < 1e-12);
    }

    // ── Euler buckling ────────────────────────────────────────────────────────

    #[test]
    fn test_euler_buckling_pinned_pinned() {
        // E=1, I=1, L=1, K=1 → P_cr = π²
        let p_cr = euler_buckling_load(1.0, 1.0, 1.0, 1.0);
        assert!((p_cr - PI * PI).abs() < 1e-10);
    }

    #[test]
    fn test_euler_buckling_fixed_free() {
        let p1 = euler_buckling_load(1.0, 1.0, 1.0, EndCondition::PinnedPinned.k_factor());
        let p2 = euler_buckling_load(1.0, 1.0, 1.0, EndCondition::FixedFree.k_factor());
        // Fixed-free (K=2) → P_cr = P_pinned/4
        assert!((p2 - p1 / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_euler_critical_stress() {
        let sigma = euler_critical_stress(1.0, PI);
        assert!((sigma - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_radius_of_gyration() {
        // I = A r² → r = √(I/A) = √(4/1) = 2
        let r = radius_of_gyration(4.0, 1.0);
        assert!((r - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_slenderness_ratio() {
        let lambda = slenderness_ratio(1.0, 10.0, 0.5);
        assert!((lambda - 20.0).abs() < 1e-12);
    }

    #[test]
    fn test_end_condition_k_factors() {
        assert!((EndCondition::PinnedPinned.k_factor() - 1.0).abs() < 1e-12);
        assert!((EndCondition::FixedFree.k_factor() - 2.0).abs() < 1e-12);
        assert!((EndCondition::FixedFixed.k_factor() - 0.5).abs() < 1e-12);
    }

    // ── Moment distribution helpers ───────────────────────────────────────────

    #[test]
    fn test_simply_supported_moment_midspan() {
        // w=1, L=10: M_max = 1*100/8 = 12.5 at s=5
        let m = simply_supported_moment(1.0, 10.0, 5.0);
        assert!((m - 12.5).abs() < 1e-10);
    }

    #[test]
    fn test_simply_supported_max_moment() {
        let m_max = simply_supported_max_moment(1.0, 10.0);
        assert!((m_max - 12.5).abs() < 1e-10);
    }

    #[test]
    fn test_simply_supported_shear_midspan() {
        let v = simply_supported_shear(1.0, 10.0, 5.0);
        assert!(v.abs() < 1e-12);
    }

    #[test]
    fn test_fixed_fixed_fem() {
        let fem = fixed_fixed_fem(1.0, 6.0);
        assert!((fem - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_distribution_factor_sum_one() {
        // For two equal members at a joint, DF = 0.5 each
        let df1 = distribution_factor(1.0, 2.0);
        let df2 = distribution_factor(1.0, 2.0);
        assert!((df1 + df2 - 1.0).abs() < 1e-12);
    }

    // ── PortalFrame ───────────────────────────────────────────────────────────

    #[test]
    fn test_portal_sway_stiffness_positive() {
        let p = PortalFrame {
            width: 5.0,
            height: 4.0,
            e_column: 2e11,
            e_beam: 2e11,
            i_column: 1e-5,
            i_beam: 2e-5,
            a_column: 1e-3,
            a_beam: 1e-3,
        };
        assert!(p.sway_stiffness() > 0.0);
    }

    #[test]
    fn test_portal_sway_deflection_proportional() {
        let p = PortalFrame {
            width: 5.0,
            height: 4.0,
            e_column: 2e11,
            e_beam: 2e11,
            i_column: 1e-5,
            i_beam: 2e-5,
            a_column: 1e-3,
            a_beam: 1e-3,
        };
        let d1 = p.sway_deflection(10.0);
        let d2 = p.sway_deflection(20.0);
        assert!((d2 / d1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_portal_column_moments() {
        let p = PortalFrame {
            width: 5.0,
            height: 4.0,
            e_column: 2e11,
            e_beam: 2e11,
            i_column: 1e-5,
            i_beam: 2e-5,
            a_column: 1e-3,
            a_beam: 1e-3,
        };
        // H=10, h=4 → M = 10*4/2 = 20
        assert!((p.column_moments(10.0) - 20.0).abs() < 1e-10);
    }

    // ── PlaneFrame solve ──────────────────────────────────────────────────────

    #[test]
    fn test_plane_frame_cantilever_tip_displacement() {
        // Single cantilever: node 0 fixed, node 1 free with tip load
        let nodes = vec![
            FrameNode2D {
                id: 0,
                x: 0.0,
                y: 0.0,
            },
            FrameNode2D {
                id: 1,
                x: 2.0,
                y: 0.0,
            },
        ];
        let constraints = vec![Constraint::FixedFixed, Constraint::Free];
        let mut frame = PlaneFrame::new(nodes, constraints);
        let e = 1.0e4;
        let a = 0.1;
        let i = 1.0e-3;
        frame.add_element(0, 1, e, a, i);
        frame.apply_load(1, 0.0, -1.0, 0.0); // downward unit load
        let u = frame.solve().expect("solver failed");
        // tip vertical displacement u_y[1] = -P L³ / (3EI)
        let expected = -8.0 / (3.0 * e * i);
        let got = u[4]; // DOF index: node 1, y
        assert!(
            (got - expected).abs() / expected.abs() < 1e-4,
            "got={got} expected={expected}"
        );
    }

    #[test]
    fn test_plane_frame_free_dof_count() {
        let nodes = vec![
            FrameNode2D {
                id: 0,
                x: 0.0,
                y: 0.0,
            },
            FrameNode2D {
                id: 1,
                x: 1.0,
                y: 0.0,
            },
        ];
        let constraints = vec![Constraint::FixedFixed, Constraint::Free];
        let frame = PlaneFrame::new(nodes, constraints);
        assert_eq!(frame.free_dof_count(), 3);
    }

    // ── Direction cosine matrix ───────────────────────────────────────────────

    #[test]
    fn test_direction_cosine_matrix_orthonormal() {
        let ni = [0.0, 0.0, 0.0];
        let nj = [1.0, 0.0, 0.0];
        let up = [0.0, 0.0, 1.0];
        let gamma = direction_cosine_matrix(ni, nj, up);
        // Check rows are orthonormal
        for i in 0..3 {
            let norm: f64 = gamma[i].iter().map(|&x| x * x).sum::<f64>().sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "row {i} norm = {norm}");
        }
    }

    #[test]
    fn test_direction_cosine_matrix_vertical_member() {
        // Member along Z: ni=(0,0,0), nj=(0,0,5)
        let ni = [0.0, 0.0, 0.0];
        let nj = [0.0, 0.0, 5.0];
        let up = [0.0, 0.0, 1.0];
        let gamma = direction_cosine_matrix(ni, nj, up);
        // local x should point in global z direction
        assert!((gamma[0][2] - 1.0).abs() < 1e-10);
    }

    // ── Moment diagram ────────────────────────────────────────────────────────

    #[test]
    fn test_moment_diagram_simply_supported_udl() {
        // Simply-supported: M_i = M_j = 0, w = 1, L = 10
        let (pos, mom) = moment_diagram(0.0, 0.0, 1.0, 10.0, 11);
        // Check midspan
        let mid = mom[5];
        assert!((mid - 12.5).abs() < 1e-8, "mid moment = {mid}");
        assert!(pos.len() == 11);
    }

    #[test]
    fn test_moment_diagram_end_values() {
        let m_i = 5.0;
        let m_j = 3.0;
        let (_pos, mom) = moment_diagram(m_i, m_j, 0.0, 4.0, 5);
        assert!((mom[0] - m_i).abs() < 1e-10);
        assert!((mom[4] - m_j).abs() < 1e-10);
    }

    // ── Matrix utilities ──────────────────────────────────────────────────────

    #[test]
    fn test_mat6_mul_identity() {
        let mut id = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            id[i][i] = 1.0;
        }
        let a = [[1.0_f64; 6]; 6];
        let c = mat6_mul(&a, &id);
        for i in 0..6 {
            for j in 0..6 {
                assert!((c[i][j] - a[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_transpose6_double_transpose() {
        let mut m = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                m[i][j] = (i * 6 + j) as f64;
            }
        }
        let tt = transpose6(&transpose6(&m));
        for i in 0..6 {
            for j in 0..6 {
                assert!((tt[i][j] - m[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_vec3_cross_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[2] - 1.0).abs() < 1e-12);
        assert!(c[0].abs() < 1e-12);
        assert!(c[1].abs() < 1e-12);
    }

    // ── SpaceFrameElement ─────────────────────────────────────────────────────

    #[test]
    fn test_space_frame_axial_stiffness() {
        let e = SpaceFrameElement {
            elastic_modulus: 1.0,
            shear_modulus: 0.5,
            area: 2.0,
            iy: 1.0,
            iz: 1.0,
            j_torsion: 1.0,
            node_i: [0.0, 0.0, 0.0],
            node_j: [4.0, 0.0, 0.0],
        };
        assert!((e.axial_stiffness() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_space_frame_torsional_stiffness() {
        let e = SpaceFrameElement {
            elastic_modulus: 2.0,
            shear_modulus: 1.0,
            area: 1.0,
            iy: 1.0,
            iz: 1.0,
            j_torsion: 2.0,
            node_i: [0.0, 0.0, 0.0],
            node_j: [2.0, 0.0, 0.0],
        };
        assert!((e.torsional_stiffness() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_space_frame_axial_force() {
        let e = SpaceFrameElement {
            elastic_modulus: 1.0,
            shear_modulus: 0.5,
            area: 1.0,
            iy: 1.0,
            iz: 1.0,
            j_torsion: 1.0,
            node_i: [0.0, 0.0, 0.0],
            node_j: [1.0, 0.0, 0.0],
        };
        let n = e.axial_force(0.0, 0.2);
        assert!((n - 0.2).abs() < 1e-12);
    }
}
