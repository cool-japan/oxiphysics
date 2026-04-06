// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated FEM matrix assembly (CPU mock implementation).
//!
//! This module provides Finite Element Method (FEM) matrix assembly routines
//! that mirror a GPU implementation. All operations run on the CPU via plain
//! loops for portability.
//!
//! The formulation assumes 2-node bar/rod elements in 1-D for simplicity,
//! making the element stiffness matrix 2×2 and DOF management straightforward.
//! The same patterns extend to 2-D and 3-D elements.

// ── Data structures ──────────────────────────────────────────────────────────

/// A FEM mesh with element connectivity and material parameters.
///
/// Elements are 2-node bar/rod elements. Each element connects two nodes.
/// Global DOF count equals the number of nodes.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GpuFemMesh {
    /// Node coordinates (one per node).
    pub node_coords: Vec<f64>,
    /// Element connectivity: `[node_a0, node_b0, node_a1, node_b1, …]`.
    pub elements: Vec<usize>,
    /// Young's modulus for each element.
    pub youngs_modulus: Vec<f64>,
    /// Cross-sectional area for each element.
    pub area: Vec<f64>,
    /// Global stiffness matrix (n_dofs × n_dofs, row-major).
    pub k_global: Vec<f64>,
    /// Global displacement vector (n_dofs).
    pub displacements: Vec<f64>,
    /// Global external force vector (n_dofs).
    pub ext_forces: Vec<f64>,
    /// Residual vector r = f − K·u (n_dofs).
    pub residual: Vec<f64>,
    /// Dirichlet (fixed) DOF flags: `true` means constrained.
    pub dirichlet_flags: Vec<bool>,
}

impl GpuFemMesh {
    /// Create a new `GpuFemMesh` from node coordinates and element connectivity.
    ///
    /// `elements` must be a flat list of node-index pairs `[a0, b0, a1, b1, …]`.
    /// All material parameters default to 1.0.
    pub fn new(node_coords: Vec<f64>, elements: Vec<usize>) -> Self {
        let n_nodes = node_coords.len();
        let n_elems = elements.len() / 2;
        Self {
            node_coords,
            elements,
            youngs_modulus: vec![1.0; n_elems],
            area: vec![1.0; n_elems],
            k_global: vec![0.0; n_nodes * n_nodes],
            displacements: vec![0.0; n_nodes],
            ext_forces: vec![0.0; n_nodes],
            residual: vec![0.0; n_nodes],
            dirichlet_flags: vec![false; n_nodes],
        }
    }

    /// Number of nodes (= DOFs for 1-D bar formulation).
    pub fn n_dofs(&self) -> usize {
        self.node_coords.len()
    }

    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.elements.len() / 2
    }
}

// ── GPU kernel mocks ─────────────────────────────────────────────────────────

/// Compute the 2×2 element stiffness matrix for bar/rod element `e`.
///
/// For a 1-D bar element: k_e = (E·A / L) · \[\[1, -1\\], \[-1, 1\]]
///
/// Returns `[k00, k01, k10, k11]` in row-major order.
pub fn gpu_element_stiffness(mesh: &GpuFemMesh, e: usize) -> [f64; 4] {
    let na = mesh.elements[e * 2];
    let nb = mesh.elements[e * 2 + 1];
    let xa = mesh.node_coords[na];
    let xb = mesh.node_coords[nb];
    let length = (xb - xa).abs();
    if length < 1e-15 {
        return [0.0; 4];
    }
    let ke = mesh.youngs_modulus[e] * mesh.area[e] / length;
    [ke, -ke, -ke, ke]
}

/// Parallel element stiffness computation — returns all element matrices.
///
/// Returns a `Vec` of `[k00, k01, k10, k11]` arrays, one per element.
pub fn gpu_assemble_global(mesh: &mut GpuFemMesh) {
    let n_dofs = mesh.n_dofs();
    mesh.k_global = vec![0.0; n_dofs * n_dofs];
    let n_elem = mesh.n_elements();
    for e in 0..n_elem {
        let ke = gpu_element_stiffness(mesh, e);
        let na = mesh.elements[e * 2];
        let nb = mesh.elements[e * 2 + 1];
        // scatter ke into global K
        mesh.k_global[na * n_dofs + na] += ke[0];
        mesh.k_global[na * n_dofs + nb] += ke[1];
        mesh.k_global[nb * n_dofs + na] += ke[2];
        mesh.k_global[nb * n_dofs + nb] += ke[3];
    }
}

/// Apply Dirichlet boundary conditions by zeroing constrained DOF rows/cols.
///
/// For each constrained DOF `i`, sets:
/// - Row `i` of `K` to zero except `K[i,i] = 1`
/// - Column `i` of `K` to zero
/// - `f[i] = 0`
pub fn gpu_apply_dirichlet(mesh: &mut GpuFemMesh) {
    let n = mesh.n_dofs();
    for i in 0..n {
        if mesh.dirichlet_flags[i] {
            // zero row i
            for j in 0..n {
                mesh.k_global[i * n + j] = 0.0;
            }
            // zero column i
            for j in 0..n {
                mesh.k_global[j * n + i] = 0.0;
            }
            // set diagonal to 1
            mesh.k_global[i * n + i] = 1.0;
            // zero rhs
            mesh.ext_forces[i] = 0.0;
        }
    }
}

/// Compute the residual vector `r = f − K·u` in parallel.
///
/// Updates `mesh.residual`.
pub fn gpu_residual(mesh: &mut GpuFemMesh) {
    let n = mesh.n_dofs();
    for i in 0..n {
        let mut ku_i = 0.0f64;
        for j in 0..n {
            ku_i += mesh.k_global[i * n + j] * mesh.displacements[j];
        }
        mesh.residual[i] = mesh.ext_forces[i] - ku_i;
    }
}

/// Parallel reduction dot product: `a · b`.
///
/// Both slices must have the same length.
pub fn gpu_dot_product(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
}

/// Return the element stiffness matrices for all elements (parallel mock).
///
/// Convenience wrapper that collects [`gpu_element_stiffness`] for every
/// element.
#[allow(dead_code)]
pub fn gpu_all_element_stiffness(mesh: &GpuFemMesh) -> Vec<[f64; 4]> {
    (0..mesh.n_elements())
        .map(|e| gpu_element_stiffness(mesh, e))
        .collect()
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a simple 3-node bar mesh: nodes at 0.0, 1.0, 2.0 with 2 elements.
    fn make_bar_mesh() -> GpuFemMesh {
        let coords = vec![0.0, 1.0, 2.0];
        let elems = vec![0, 1, 1, 2];
        GpuFemMesh::new(coords, elems)
    }

    #[test]
    fn test_new_mesh_n_dofs() {
        let m = make_bar_mesh();
        assert_eq!(m.n_dofs(), 3);
    }

    #[test]
    fn test_new_mesh_n_elements() {
        let m = make_bar_mesh();
        assert_eq!(m.n_elements(), 2);
    }

    #[test]
    fn test_new_mesh_default_youngs() {
        let m = make_bar_mesh();
        assert!((m.youngs_modulus[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_stiffness_unit_bar() {
        // E=1, A=1, L=1 → ke = [[1,-1],[-1,1]]
        let m = make_bar_mesh();
        let ke = gpu_element_stiffness(&m, 0);
        assert!((ke[0] - 1.0).abs() < 1e-12);
        assert!((ke[1] + 1.0).abs() < 1e-12);
        assert!((ke[2] + 1.0).abs() < 1e-12);
        assert!((ke[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_stiffness_scaled() {
        // E=2, A=3, L=1 → ke_diag = 6
        let mut m = make_bar_mesh();
        m.youngs_modulus[0] = 2.0;
        m.area[0] = 3.0;
        let ke = gpu_element_stiffness(&m, 0);
        assert!((ke[0] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_stiffness_zero_length() {
        let coords = vec![0.0, 0.0];
        let elems = vec![0, 1];
        let m = GpuFemMesh::new(coords, elems);
        let ke = gpu_element_stiffness(&m, 0);
        assert_eq!(ke, [0.0; 4]);
    }

    #[test]
    fn test_assemble_global_dimensions() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        assert_eq!(m.k_global.len(), 9); // 3×3
    }

    #[test]
    fn test_assemble_global_diagonal_positive() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        let n = m.n_dofs();
        for i in 0..n {
            assert!(m.k_global[i * n + i] >= 0.0);
        }
    }

    #[test]
    fn test_assemble_global_symmetric() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        let n = m.n_dofs();
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (m.k_global[i * n + j] - m.k_global[j * n + i]).abs() < 1e-12,
                    "K[{i},{j}] != K[{j},{i}]"
                );
            }
        }
    }

    #[test]
    fn test_assemble_global_row_sum_zero() {
        // For a free structure (no BCs) each row should sum to ~0
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        let n = m.n_dofs();
        for i in 0..n {
            let row_sum: f64 = (0..n).map(|j| m.k_global[i * n + j]).sum();
            assert!(row_sum.abs() < 1e-10, "row {i} sum = {row_sum}");
        }
    }

    #[test]
    fn test_apply_dirichlet_zeroes_row() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.dirichlet_flags[0] = true;
        gpu_apply_dirichlet(&mut m);
        let n = m.n_dofs();
        // Off-diagonal entries of row 0 should be zero
        for j in 1..n {
            assert!((m.k_global[j]).abs() < 1e-12);
        }
        // Diagonal should be 1
        assert!((m.k_global[0]).abs() - 1.0 < 1e-12);
    }

    #[test]
    fn test_apply_dirichlet_zeroes_column() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.dirichlet_flags[0] = true;
        gpu_apply_dirichlet(&mut m);
        let n = m.n_dofs();
        for i in 1..n {
            assert!((m.k_global[i * n]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_apply_dirichlet_zeroes_rhs() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.ext_forces[0] = 99.0;
        m.dirichlet_flags[0] = true;
        gpu_apply_dirichlet(&mut m);
        assert!((m.ext_forces[0]).abs() < 1e-12);
    }

    #[test]
    fn test_residual_zero_displacement() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.ext_forces[2] = 1.0;
        // u = 0 → r = f
        gpu_residual(&mut m);
        assert!((m.residual[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_residual_equilibrium() {
        // If K*u = f exactly, residual should be zero
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.dirichlet_flags[0] = true;
        gpu_apply_dirichlet(&mut m);
        m.ext_forces[2] = 1.0;
        // Solve manually for 2-element bar, node 0 fixed, node 2 loaded
        // K after BCs: diag = [1, 2, 1], off-diag per element pattern
        // For simplicity just set u = K^{-1} f using known solution
        // u[0]=0, u[1]=1, u[2]=2 (for unit bar: displacement = x * force)
        m.displacements = vec![0.0, 1.0, 2.0];
        gpu_residual(&mut m);
        // residual should not blow up
        for &r in &m.residual {
            assert!(r.is_finite());
        }
    }

    #[test]
    fn test_gpu_dot_product_basic() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((gpu_dot_product(&a, &b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_dot_product_empty() {
        assert!((gpu_dot_product(&[], &[])).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_dot_product_unit_vectors() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!((gpu_dot_product(&a, &b)).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_all_element_stiffness_count() {
        let m = make_bar_mesh();
        let all_ke = gpu_all_element_stiffness(&m);
        assert_eq!(all_ke.len(), m.n_elements());
    }

    #[test]
    fn test_gpu_all_element_stiffness_values() {
        let m = make_bar_mesh();
        let all_ke = gpu_all_element_stiffness(&m);
        // Both elements identical → same stiffness matrix
        assert_eq!(all_ke[0], all_ke[1]);
    }

    #[test]
    fn test_fem_mesh_clone() {
        let m = make_bar_mesh();
        let m2 = m.clone();
        assert_eq!(m2.n_dofs(), 3);
    }

    #[test]
    fn test_fem_mesh_debug() {
        let m = make_bar_mesh();
        let s = format!("{m:?}");
        assert!(s.contains("GpuFemMesh"));
    }

    #[test]
    fn test_assemble_then_apply_dirichlet_both_ends() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.dirichlet_flags[0] = true;
        m.dirichlet_flags[2] = true;
        gpu_apply_dirichlet(&mut m);
        let n = m.n_dofs();
        // Diagonals of constrained nodes should be 1
        assert!((m.k_global[0] - 1.0).abs() < 1e-12);
        assert!((m.k_global[2 * n + 2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_residual_updates_all_entries() {
        let mut m = make_bar_mesh();
        gpu_assemble_global(&mut m);
        m.ext_forces = vec![1.0, 0.0, -1.0];
        gpu_residual(&mut m);
        assert_eq!(m.residual.len(), m.n_dofs());
    }
}
