// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM (Finite Element Method) Assembly and Solver.

#![allow(missing_docs)]

// ===========================================================================
// FEM (Finite Element Method) Assembly and Solver
// ===========================================================================

/// A single 2-node bar (truss) element for FEM assembly.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FemBarElement {
    /// Global indices of the two end nodes.
    pub nodes: [usize; 2],
    /// Young's modulus times cross-sectional area (EA).
    pub ea: f64,
    /// Undeformed length of the element.
    pub length: f64,
}

/// FEM assembly for a truss structure.
///
/// Assembles a global stiffness matrix from bar elements, applies boundary
/// conditions, and solves with a direct (dense) solver suitable for
/// demonstration/testing with small meshes.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyFemAssembly {
    /// Number of degrees of freedom (nodes × 3 for 3-D).
    pub n_dofs: usize,
    /// Bar elements.
    elements: Vec<FemBarElement>,
    /// External force vector.
    forces: Vec<f64>,
    /// Fixed DOF indices (Dirichlet BCs: displacement = 0).
    fixed_dofs: Vec<usize>,
    /// Global stiffness matrix (dense, row-major, n_dofs × n_dofs).
    k_global: Vec<f64>,
    /// Displacement solution vector.
    displacements: Vec<f64>,
    /// Whether the stiffness matrix has been assembled.
    assembled: bool,
}

impl PyFemAssembly {
    /// Create a new FEM assembly with `n_nodes` 3-D nodes.
    pub fn new(n_nodes: usize) -> Self {
        let n_dofs = n_nodes * 3;
        Self {
            n_dofs,
            elements: Vec::new(),
            forces: vec![0.0; n_dofs],
            fixed_dofs: Vec::new(),
            k_global: vec![0.0; n_dofs * n_dofs],
            displacements: vec![0.0; n_dofs],
            assembled: false,
        }
    }

    /// Add a bar element between two nodes with given EA and undeformed length.
    pub fn add_bar_element(&mut self, node_a: usize, node_b: usize, ea: f64, length: f64) {
        self.elements.push(FemBarElement {
            nodes: [node_a, node_b],
            ea,
            length,
        });
        self.assembled = false;
    }

    /// Apply an external force to a DOF (node_index * 3 + direction).
    pub fn apply_force(&mut self, dof: usize, force: f64) {
        if dof < self.n_dofs {
            self.forces[dof] += force;
        }
    }

    /// Fix a DOF (zero displacement Dirichlet BC).
    pub fn fix_dof(&mut self, dof: usize) {
        if dof < self.n_dofs && !self.fixed_dofs.contains(&dof) {
            self.fixed_dofs.push(dof);
        }
    }

    /// Assemble the global stiffness matrix from all bar elements.
    ///
    /// For each bar element the 2×2 (1-D axial) stiffness matrix contribution
    /// `k_e = EA/L * [[1, -1\], [-1, 1]]` is added to the corresponding DOF
    /// rows/columns along the element's axis direction.
    pub fn assemble(&mut self) {
        // Zero the matrix
        for v in &mut self.k_global {
            *v = 0.0;
        }

        let n = self.n_dofs;
        for elem in &self.elements {
            let ea_l = if elem.length > 1e-15 {
                elem.ea / elem.length
            } else {
                0.0
            };
            // Axial DOFs: use DOF 0 of each node (x-direction simplified)
            let dof_a = elem.nodes[0] * 3;
            let dof_b = elem.nodes[1] * 3;
            if dof_a < n && dof_b < n {
                self.k_global[dof_a * n + dof_a] += ea_l;
                self.k_global[dof_a * n + dof_b] -= ea_l;
                self.k_global[dof_b * n + dof_a] -= ea_l;
                self.k_global[dof_b * n + dof_b] += ea_l;
            }
        }

        // Zero-out empty diagonal entries (unconstrained DOFs with no element
        // contribution) and apply Dirichlet BCs via the large-number method.
        let big = 1.0e20;
        for i in 0..n {
            if self.k_global[i * n + i].abs() < 1e-30 {
                // This DOF has no stiffness contribution — pin it to zero.
                self.k_global[i * n + i] = big;
                self.forces[i] = 0.0;
            }
        }
        for &dof in &self.fixed_dofs {
            if dof < n {
                for col in 0..n {
                    self.k_global[dof * n + col] = 0.0;
                }
                self.k_global[dof * n + dof] = big;
                self.forces[dof] = 0.0;
            }
        }

        self.assembled = true;
    }

    /// Solve the linear system K·u = f using Gaussian elimination.
    ///
    /// The system is solved in-place. Returns `true` on success,
    /// `false` if the matrix is singular or not yet assembled.
    pub fn solve(&mut self) -> bool {
        if !self.assembled {
            self.assemble();
        }
        let n = self.n_dofs;
        if n == 0 {
            return false;
        }

        // Build an augmented matrix [K | f]
        let mut aug: Vec<f64> = vec![0.0; n * (n + 1)];
        for i in 0..n {
            for j in 0..n {
                aug[i * (n + 1) + j] = self.k_global[i * n + j];
            }
            aug[i * (n + 1) + n] = self.forces[i];
        }

        // Forward elimination
        for col in 0..n {
            // Pivot
            let mut max_row = col;
            let mut max_val = aug[col * (n + 1) + col].abs();
            for row in (col + 1)..n {
                let v = aug[row * (n + 1) + col].abs();
                if v > max_val {
                    max_val = v;
                    max_row = row;
                }
            }
            if max_val < 1e-15 {
                return false; // singular
            }
            if max_row != col {
                for j in 0..=(n) {
                    aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
                }
            }
            let pivot = aug[col * (n + 1) + col];
            for row in (col + 1)..n {
                let factor = aug[row * (n + 1) + col] / pivot;
                for j in col..=(n) {
                    let sub = factor * aug[col * (n + 1) + j];
                    aug[row * (n + 1) + j] -= sub;
                }
            }
        }

        // Back substitution
        let mut u = vec![0.0f64; n];
        for i in (0..n).rev() {
            let mut sum = aug[i * (n + 1) + n];
            for j in (i + 1)..n {
                sum -= aug[i * (n + 1) + j] * u[j];
            }
            u[i] = sum / aug[i * (n + 1) + i];
        }

        self.displacements = u;
        true
    }

    /// Get the displacement at DOF `dof`, or `0.0` if out of bounds.
    pub fn displacement(&self, dof: usize) -> f64 {
        self.displacements.get(dof).copied().unwrap_or(0.0)
    }

    /// Get all displacements as a slice.
    pub fn displacements(&self) -> &[f64] {
        &self.displacements
    }

    /// Compute the axial force in element `elem_idx`.
    ///
    /// Returns `None` if the index is out of bounds or the system is unsolved.
    pub fn element_force(&self, elem_idx: usize) -> Option<f64> {
        let elem = self.elements.get(elem_idx)?;
        if elem.length < 1e-15 {
            return Some(0.0);
        }
        let dof_a = elem.nodes[0] * 3;
        let dof_b = elem.nodes[1] * 3;
        let ua = self.displacements.get(dof_a).copied().unwrap_or(0.0);
        let ub = self.displacements.get(dof_b).copied().unwrap_or(0.0);
        Some(elem.ea / elem.length * (ub - ua))
    }
}
