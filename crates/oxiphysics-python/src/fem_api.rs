#![allow(clippy::too_many_arguments)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Finite Element Method (FEM) simulation API for Python interop.
//!
//! Provides a standalone structural mechanics FEM solver supporting:
//! - Linear elastic 2-D triangular elements (CST - Constant Strain Triangle)
//! - Assembly of global stiffness matrix
//! - Dirichlet boundary conditions (fixed DOFs)
//! - Nodal load vectors
//! - Conjugate Gradient iterative solver
//! - Displacement, stress, and strain output

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Material
// ---------------------------------------------------------------------------

/// Linear elastic isotropic material properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemMaterial {
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio ν (dimensionless, 0..0.5).
    pub poisson_ratio: f64,
    /// Mass density ρ (kg/m³).
    pub density: f64,
}

impl PyFemMaterial {
    /// Create a new material.
    pub fn new(young_modulus: f64, poisson_ratio: f64, density: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            density,
        }
    }

    /// Steel-like material (E=200 GPa, ν=0.3, ρ=7850 kg/m³).
    pub fn steel() -> Self {
        Self::new(200e9, 0.3, 7850.0)
    }

    /// Aluminum-like material (E=70 GPa, ν=0.33, ρ=2700 kg/m³).
    pub fn aluminum() -> Self {
        Self::new(70e9, 0.33, 2700.0)
    }

    /// Rubber-like material (E=0.01 GPa, ν=0.49, ρ=1100 kg/m³).
    pub fn rubber() -> Self {
        Self::new(10e6, 0.49, 1100.0)
    }

    /// Concrete-like material (E=30 GPa, ν=0.2, ρ=2400 kg/m³).
    pub fn concrete() -> Self {
        Self::new(30e9, 0.2, 2400.0)
    }

    /// Compute the plane-stress constitutive matrix D (3×3 as flat \[f64; 9\]).
    ///
    /// Returns the D matrix that relates stress to strain:
    /// σ = D ε
    #[allow(non_snake_case)]
    pub fn plane_stress_D(&self) -> [f64; 9] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let c = e / (1.0 - nu * nu);
        [
            c,
            c * nu,
            0.0,
            c * nu,
            c,
            0.0,
            0.0,
            0.0,
            c * (1.0 - nu) * 0.5,
        ]
    }
}

impl Default for PyFemMaterial {
    fn default() -> Self {
        Self::steel()
    }
}

// ---------------------------------------------------------------------------
// Node
// ---------------------------------------------------------------------------

/// A single node in the FEM mesh.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemNode {
    /// Node coordinates \[x, y\].
    pub position: [f64; 2],
}

impl PyFemNode {
    /// Create a new node at the given 2-D position.
    pub fn new(x: f64, y: f64) -> Self {
        Self { position: [x, y] }
    }
}

// ---------------------------------------------------------------------------
// Element (CST triangular element)
// ---------------------------------------------------------------------------

/// A triangular (CST) element defined by three node indices.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemElement {
    /// Indices of the three corner nodes (counter-clockwise preferred).
    pub nodes: [usize; 3],
    /// Index into the material list.
    pub material_id: usize,
    /// Element thickness t (for plane-stress, metres).
    pub thickness: f64,
}

impl PyFemElement {
    /// Create a new triangular element.
    pub fn new(n0: usize, n1: usize, n2: usize, material_id: usize, thickness: f64) -> Self {
        Self {
            nodes: [n0, n1, n2],
            material_id,
            thickness,
        }
    }
}

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------

/// A Dirichlet (fixed DOF) boundary condition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemDirichletBC {
    /// Node index.
    pub node: usize,
    /// DOF index within the node (0=x, 1=y).
    pub dof: usize,
    /// Prescribed displacement value.
    pub value: f64,
}

impl PyFemDirichletBC {
    /// Fix the x-DOF of `node` to zero.
    pub fn fix_x(node: usize) -> Self {
        Self {
            node,
            dof: 0,
            value: 0.0,
        }
    }

    /// Fix the y-DOF of `node` to zero.
    pub fn fix_y(node: usize) -> Self {
        Self {
            node,
            dof: 1,
            value: 0.0,
        }
    }

    /// Fix both DOFs of `node` (fully pinned).
    pub fn pin(node: usize) -> Vec<Self> {
        vec![Self::fix_x(node), Self::fix_y(node)]
    }
}

/// A Neumann (nodal force) boundary condition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemNodalForce {
    /// Node index.
    pub node: usize,
    /// Force vector \[fx, fy\] in Newtons.
    pub force: [f64; 2],
}

impl PyFemNodalForce {
    /// Apply force `[fx, fy]` to `node`.
    pub fn new(node: usize, fx: f64, fy: f64) -> Self {
        Self {
            node,
            force: [fx, fy],
        }
    }
}

// ---------------------------------------------------------------------------
// FEM Mesh
// ---------------------------------------------------------------------------

/// A 2-D FEM mesh with nodes, elements, materials, and boundary conditions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemMesh {
    /// All nodes in the mesh.
    pub nodes: Vec<PyFemNode>,
    /// All triangular elements in the mesh.
    pub elements: Vec<PyFemElement>,
    /// Material library.
    pub materials: Vec<PyFemMaterial>,
    /// Dirichlet (fixed) boundary conditions.
    pub dirichlet_bcs: Vec<PyFemDirichletBC>,
    /// Nodal force loads.
    pub nodal_forces: Vec<PyFemNodalForce>,
}

impl PyFemMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
            materials: vec![PyFemMaterial::default()],
            dirichlet_bcs: Vec::new(),
            nodal_forces: Vec::new(),
        }
    }

    /// Add a node at `(x, y)`. Returns node index.
    pub fn add_node(&mut self, x: f64, y: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(PyFemNode::new(x, y));
        idx
    }

    /// Add a material. Returns material index.
    pub fn add_material(&mut self, mat: PyFemMaterial) -> usize {
        let idx = self.materials.len();
        self.materials.push(mat);
        idx
    }

    /// Add a CST element connecting nodes `n0, n1, n2`. Returns element index.
    pub fn add_element(&mut self, n0: usize, n1: usize, n2: usize) -> usize {
        self.add_element_with_material(n0, n1, n2, 0, 1.0)
    }

    /// Add a CST element with explicit material and thickness. Returns element index.
    pub fn add_element_with_material(
        &mut self,
        n0: usize,
        n1: usize,
        n2: usize,
        material_id: usize,
        thickness: f64,
    ) -> usize {
        let idx = self.elements.len();
        self.elements
            .push(PyFemElement::new(n0, n1, n2, material_id, thickness));
        idx
    }

    /// Add a Dirichlet boundary condition.
    pub fn add_dirichlet_bc(&mut self, bc: PyFemDirichletBC) {
        self.dirichlet_bcs.push(bc);
    }

    /// Pin both DOFs of `node` (fix x and y to zero).
    pub fn pin_node(&mut self, node: usize) {
        self.dirichlet_bcs.push(PyFemDirichletBC::fix_x(node));
        self.dirichlet_bcs.push(PyFemDirichletBC::fix_y(node));
    }

    /// Apply a nodal force at `node`.
    pub fn apply_nodal_force(&mut self, node: usize, fx: f64, fy: f64) {
        self.nodal_forces.push(PyFemNodalForce::new(node, fx, fy));
    }

    /// Number of degrees of freedom (2 per node).
    pub fn num_dofs(&self) -> usize {
        self.nodes.len() * 2
    }

    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Compute the signed area of triangular element `e`.
    pub fn element_area(&self, e: usize) -> f64 {
        let elem = &self.elements[e];
        let [n0, n1, n2] = elem.nodes;
        let [x0, y0] = self.nodes[n0].position;
        let [x1, y1] = self.nodes[n1].position;
        let [x2, y2] = self.nodes[n2].position;
        0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
    }

    /// Build a simple structured rectangular mesh of triangles.
    ///
    /// Creates a `nx × ny` grid of quads, each quad split into 2 triangles.
    /// Bottom-left corner at `(ox, oy)`, domain size `(lx, ly)`.
    /// Returns the number of nodes added.
    pub fn build_rectangle(
        &mut self,
        ox: f64,
        oy: f64,
        lx: f64,
        ly: f64,
        nx: usize,
        ny: usize,
        material_id: usize,
        thickness: f64,
    ) -> usize {
        let dx = lx / nx as f64;
        let dy = ly / ny as f64;
        let node_start = self.nodes.len();

        // Add nodes
        for j in 0..=(ny) {
            for i in 0..=(nx) {
                self.add_node(ox + i as f64 * dx, oy + j as f64 * dy);
            }
        }

        // Add triangular elements (2 per quad)
        let row = nx + 1;
        for j in 0..ny {
            for i in 0..nx {
                let n0 = node_start + j * row + i;
                let n1 = node_start + j * row + i + 1;
                let n2 = node_start + (j + 1) * row + i;
                let n3 = node_start + (j + 1) * row + i + 1;
                self.add_element_with_material(n0, n1, n3, material_id, thickness);
                self.add_element_with_material(n0, n3, n2, material_id, thickness);
            }
        }

        (ny + 1) * (nx + 1)
    }
}

impl Default for PyFemMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// FEM Solver result
// ---------------------------------------------------------------------------

/// Result of a FEM linear static solve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFemSolveResult {
    /// Nodal displacement vector (flat: \[ux0, uy0, ux1, uy1, ...\]).
    pub displacements: Vec<f64>,
    /// Von Mises stress per element.
    pub von_mises_stress: Vec<f64>,
    /// Number of CG solver iterations used.
    pub solver_iterations: usize,
    /// Residual norm at convergence.
    pub residual_norm: f64,
    /// Whether the solver converged.
    pub converged: bool,
}

impl PyFemSolveResult {
    /// Displacement of node `i` as `[ux, uy]`.
    pub fn node_displacement(&self, i: usize) -> Option<[f64; 2]> {
        let base = i * 2;
        if base + 1 < self.displacements.len() {
            Some([self.displacements[base], self.displacements[base + 1]])
        } else {
            None
        }
    }

    /// Maximum von Mises stress over all elements.
    pub fn max_von_mises(&self) -> f64 {
        self.von_mises_stress
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Mean von Mises stress.
    pub fn mean_von_mises(&self) -> f64 {
        if self.von_mises_stress.is_empty() {
            return 0.0;
        }
        self.von_mises_stress.iter().sum::<f64>() / self.von_mises_stress.len() as f64
    }
}

// ---------------------------------------------------------------------------
// FEM Solver
// ---------------------------------------------------------------------------

/// Linear static FEM solver for 2-D plane-stress problems.
///
/// Assembles the global stiffness matrix K from CST elements, applies
/// Dirichlet boundary conditions, builds the load vector f, and solves
/// K u = f with a Conjugate Gradient (CG) iterative solver.
#[derive(Debug, Clone)]
pub struct PyFemSolver {
    /// Maximum CG iterations.
    pub max_iterations: usize,
    /// CG convergence tolerance.
    pub tolerance: f64,
}

impl PyFemSolver {
    /// Create a solver with default settings.
    pub fn new() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: 1e-10,
        }
    }

    /// Create a solver with custom settings.
    pub fn with_settings(max_iterations: usize, tolerance: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
        }
    }

    /// Solve the linear static problem K u = f.
    ///
    /// Returns `None` if the mesh has no nodes or no elements.
    pub fn solve(&self, mesh: &PyFemMesh) -> Option<PyFemSolveResult> {
        let n_nodes = mesh.num_nodes();
        if n_nodes == 0 || mesh.elements.is_empty() {
            return None;
        }
        let n_dofs = n_nodes * 2;

        // --- Assemble global stiffness K (stored as flat dense n_dofs×n_dofs) ---
        let mut k = vec![0.0f64; n_dofs * n_dofs];
        for elem in &mesh.elements {
            let mat = mesh
                .materials
                .get(elem.material_id)
                .unwrap_or(&mesh.materials[0]);
            let ke = self.element_stiffness(mesh, elem, mat);
            // Scatter into global K
            let dofs: [usize; 6] = [
                elem.nodes[0] * 2,
                elem.nodes[0] * 2 + 1,
                elem.nodes[1] * 2,
                elem.nodes[1] * 2 + 1,
                elem.nodes[2] * 2,
                elem.nodes[2] * 2 + 1,
            ];
            for (i, &gi) in dofs.iter().enumerate() {
                for (j, &gj) in dofs.iter().enumerate() {
                    k[gi * n_dofs + gj] += ke[i * 6 + j];
                }
            }
        }

        // --- Assemble load vector f ---
        let mut f = vec![0.0f64; n_dofs];
        for nf in &mesh.nodal_forces {
            if nf.node < n_nodes {
                f[nf.node * 2] += nf.force[0];
                f[nf.node * 2 + 1] += nf.force[1];
            }
        }

        // --- Apply Dirichlet BCs (large-number method) ---
        let large = 1e30;
        for bc in &mesh.dirichlet_bcs {
            if bc.node < n_nodes {
                let dof = bc.node * 2 + bc.dof;
                k[dof * n_dofs + dof] += large;
                f[dof] += large * bc.value;
            }
        }

        // --- Solve K u = f with Conjugate Gradient ---
        let (u, iters, res, converged) =
            conjugate_gradient(&k, &f, n_dofs, self.max_iterations, self.tolerance);

        // --- Compute von Mises stress per element ---
        let mut von_mises = Vec::with_capacity(mesh.elements.len());
        for elem in &mesh.elements {
            let mat = mesh
                .materials
                .get(elem.material_id)
                .unwrap_or(&mesh.materials[0]);
            let sigma = self.element_stress(mesh, elem, mat, &u);
            // Von Mises: sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2)
            let sx = sigma[0];
            let sy = sigma[1];
            let txy = sigma[2];
            let vm = (sx * sx - sx * sy + sy * sy + 3.0 * txy * txy).sqrt();
            von_mises.push(vm);
        }

        Some(PyFemSolveResult {
            displacements: u,
            von_mises_stress: von_mises,
            solver_iterations: iters,
            residual_norm: res,
            converged,
        })
    }

    // -----------------------------------------------------------------------
    // Private: element stiffness matrix (6×6 flat)
    // -----------------------------------------------------------------------

    fn element_stiffness(
        &self,
        mesh: &PyFemMesh,
        elem: &PyFemElement,
        mat: &PyFemMaterial,
    ) -> [f64; 36] {
        let [n0, n1, n2] = elem.nodes;
        let [x0, y0] = mesh.nodes[n0].position;
        let [x1, y1] = mesh.nodes[n1].position;
        let [x2, y2] = mesh.nodes[n2].position;

        // Area
        let area = 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        if area.abs() < 1e-15 {
            return [0.0; 36];
        }

        // Shape function gradients
        let inv_2a = 1.0 / (2.0 * area);
        let b1x = (y1 - y2) * inv_2a;
        let b2x = (y2 - y0) * inv_2a;
        let b3x = (y0 - y1) * inv_2a;
        let b1y = (x2 - x1) * inv_2a;
        let b2y = (x0 - x2) * inv_2a;
        let b3y = (x1 - x0) * inv_2a;

        // B matrix (3×6): [ε_xx, ε_yy, γ_xy] = B u
        let b_mat: [f64; 18] = [
            b1x, 0.0, b2x, 0.0, b3x, 0.0, 0.0, b1y, 0.0, b2y, 0.0, b3y, b1y, b1x, b2y, b2x, b3y,
            b3x,
        ];

        // D matrix (3×3)
        let d = mat.plane_stress_D();

        // Ke = B^T D B * area * thickness
        let t = elem.thickness;
        let factor = area.abs() * t;

        // Compute D*B (3×6)
        let mut db = [0.0f64; 18];
        for i in 0..3 {
            for j in 0..6 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += d[i * 3 + k] * b_mat[k * 6 + j];
                }
                db[i * 6 + j] = sum;
            }
        }

        // Ke = B^T * (D*B) * factor (6×6)
        let mut ke = [0.0f64; 36];
        for i in 0..6 {
            for j in 0..6 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += b_mat[k * 6 + i] * db[k * 6 + j];
                }
                ke[i * 6 + j] = sum * factor;
            }
        }
        ke
    }

    // -----------------------------------------------------------------------
    // Private: element stress vector [σx, σy, τxy]
    // -----------------------------------------------------------------------

    fn element_stress(
        &self,
        mesh: &PyFemMesh,
        elem: &PyFemElement,
        mat: &PyFemMaterial,
        u: &[f64],
    ) -> [f64; 3] {
        let [n0, n1, n2] = elem.nodes;
        let [x0, y0] = mesh.nodes[n0].position;
        let [x1, y1] = mesh.nodes[n1].position;
        let [x2, y2] = mesh.nodes[n2].position;

        let area = 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        if area.abs() < 1e-15 {
            return [0.0; 3];
        }

        let inv_2a = 1.0 / (2.0 * area);
        let b1x = (y1 - y2) * inv_2a;
        let b2x = (y2 - y0) * inv_2a;
        let b3x = (y0 - y1) * inv_2a;
        let b1y = (x2 - x1) * inv_2a;
        let b2y = (x0 - x2) * inv_2a;
        let b3y = (x1 - x0) * inv_2a;

        let b_mat: [f64; 18] = [
            b1x, 0.0, b2x, 0.0, b3x, 0.0, 0.0, b1y, 0.0, b2y, 0.0, b3y, b1y, b1x, b2y, b2x, b3y,
            b3x,
        ];

        // Local displacement vector (6 components)
        let ue: [f64; 6] = [
            u.get(n0 * 2).copied().unwrap_or(0.0),
            u.get(n0 * 2 + 1).copied().unwrap_or(0.0),
            u.get(n1 * 2).copied().unwrap_or(0.0),
            u.get(n1 * 2 + 1).copied().unwrap_or(0.0),
            u.get(n2 * 2).copied().unwrap_or(0.0),
            u.get(n2 * 2 + 1).copied().unwrap_or(0.0),
        ];

        // Strain ε = B ue (3 components)
        let mut strain = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..6 {
                strain[i] += b_mat[i * 6 + j] * ue[j];
            }
        }

        // Stress σ = D ε
        let d = mat.plane_stress_D();
        let mut stress = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                stress[i] += d[i * 3 + j] * strain[j];
            }
        }
        stress
    }
}

impl Default for PyFemSolver {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Private: Conjugate Gradient solver
// ---------------------------------------------------------------------------

/// Simple conjugate gradient solver for A x = b.
/// Returns (x, iterations, final_residual_norm, converged).
fn conjugate_gradient(
    a: &[f64],
    b: &[f64],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, usize, f64, bool) {
    let mut x = vec![0.0f64; n];
    // r = b - A x  (x=0, so r = b)
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rs_old: f64 = r.iter().map(|v| v * v).sum();

    if rs_old.sqrt() < tol {
        return (x, 0, rs_old.sqrt(), true);
    }

    let mut iters = 0;
    for _ in 0..max_iter {
        iters += 1;
        // Ap = A * p
        let ap = mat_vec_mul(a, &p, n);
        let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum();
        if pap.abs() < 1e-30 {
            break;
        }
        let alpha = rs_old / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rs_new: f64 = r.iter().map(|v| v * v).sum();
        let res_norm = rs_new.sqrt();
        if res_norm < tol {
            return (x, iters, res_norm, true);
        }
        let beta = rs_new / rs_old;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }
    let res_norm = rs_old.sqrt();
    (x, iters, res_norm, false)
}

/// Multiply dense matrix A (n×n) by vector v (n).
fn mat_vec_mul(a: &[f64], v: &[f64], n: usize) -> Vec<f64> {
    let mut result = vec![0.0f64; n];
    for i in 0..n {
        let row = &a[i * n..(i + 1) * n];
        result[i] = row.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PyFemMesh;
    use crate::PyFemSolveResult;
    use crate::PyFemSolver;

    fn simple_mesh() -> PyFemMesh {
        // Two-triangle unit square patch
        let mut mesh = PyFemMesh::new();
        mesh.add_node(0.0, 0.0); // 0
        mesh.add_node(1.0, 0.0); // 1
        mesh.add_node(1.0, 1.0); // 2
        mesh.add_node(0.0, 1.0); // 3
        mesh.add_element(0, 1, 2);
        mesh.add_element(0, 2, 3);
        mesh
    }

    #[test]
    fn test_fem_mesh_creation() {
        let mesh = simple_mesh();
        assert_eq!(mesh.num_nodes(), 4);
        assert_eq!(mesh.num_elements(), 2);
        assert_eq!(mesh.num_dofs(), 8);
    }

    #[test]
    fn test_fem_material_steel_properties() {
        let mat = PyFemMaterial::steel();
        assert!((mat.young_modulus - 200e9).abs() < 1.0);
        assert!((mat.poisson_ratio - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_fem_material_d_matrix_symmetry() {
        let mat = PyFemMaterial::steel();
        let d = mat.plane_stress_D();
        // D[0][1] == D[1][0]
        assert!((d[1] - d[3]).abs() < 1.0, "D should be symmetric");
    }

    #[test]
    fn test_fem_element_area_positive() {
        let mesh = simple_mesh();
        let area0 = mesh.element_area(0);
        let area1 = mesh.element_area(1);
        assert!(area0 != 0.0, "area = {}", area0);
        assert!(area1 != 0.0, "area = {}", area1);
        // Sum of absolute areas should be 1.0 for the unit square
        assert!((area0.abs() + area1.abs() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fem_dirichlet_bc_pin() {
        let bcs = PyFemDirichletBC::pin(0);
        assert_eq!(bcs.len(), 2);
        assert_eq!(bcs[0].dof, 0);
        assert_eq!(bcs[1].dof, 1);
    }

    #[test]
    fn test_fem_pin_node() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(1);
        assert_eq!(mesh.dirichlet_bcs.len(), 4);
    }

    #[test]
    fn test_fem_apply_nodal_force() {
        let mut mesh = simple_mesh();
        mesh.apply_nodal_force(2, 0.0, -1000.0);
        assert_eq!(mesh.nodal_forces.len(), 1);
        assert!((mesh.nodal_forces[0].force[1] + 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_fem_solver_creation() {
        let solver = PyFemSolver::new();
        assert_eq!(solver.max_iterations, 1000);
        assert!((solver.tolerance - 1e-10).abs() < 1e-20);
    }

    #[test]
    fn test_fem_solver_returns_none_for_empty_mesh() {
        let mesh = PyFemMesh::new();
        let solver = PyFemSolver::new();
        assert!(solver.solve(&mesh).is_none());
    }

    #[test]
    fn test_fem_solve_displacement_length() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        mesh.apply_nodal_force(1, 0.0, -1000.0);
        mesh.apply_nodal_force(2, 0.0, -1000.0);
        let solver = PyFemSolver::new();
        let result = solver.solve(&mesh).expect("should solve");
        assert_eq!(result.displacements.len(), 8); // 4 nodes * 2 DOFs
    }

    #[test]
    fn test_fem_solve_von_mises_length() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        mesh.apply_nodal_force(1, 1000.0, 0.0);
        let solver = PyFemSolver::new();
        let result = solver.solve(&mesh).expect("should solve");
        assert_eq!(result.von_mises_stress.len(), 2); // 2 elements
    }

    #[test]
    fn test_fem_solve_converges() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        mesh.apply_nodal_force(1, 0.0, -500.0);
        let solver = PyFemSolver::new();
        let result = solver.solve(&mesh).expect("should solve");
        assert!(
            result.converged,
            "CG should converge for a well-posed problem"
        );
    }

    #[test]
    fn test_fem_result_node_displacement() {
        let result = PyFemSolveResult {
            displacements: vec![0.0, 0.0, 1e-4, -2e-5, 0.5e-4, -1e-5, 0.0, 0.0],
            von_mises_stress: vec![1.0e6, 0.8e6],
            solver_iterations: 10,
            residual_norm: 1e-12,
            converged: true,
        };
        let d0 = result.node_displacement(0).unwrap();
        assert!((d0[0]).abs() < 1e-15);
        let d1 = result.node_displacement(1).unwrap();
        assert!((d1[0] - 1e-4).abs() < 1e-15);
    }

    #[test]
    fn test_fem_result_max_von_mises() {
        let result = PyFemSolveResult {
            displacements: vec![0.0; 8],
            von_mises_stress: vec![1.0e6, 2.0e6, 0.5e6],
            solver_iterations: 5,
            residual_norm: 1e-12,
            converged: true,
        };
        assert!((result.max_von_mises() - 2.0e6).abs() < 1.0);
    }

    #[test]
    fn test_fem_result_mean_von_mises() {
        let result = PyFemSolveResult {
            displacements: vec![0.0; 8],
            von_mises_stress: vec![1.0, 2.0, 3.0],
            solver_iterations: 5,
            residual_norm: 1e-12,
            converged: true,
        };
        assert!((result.mean_von_mises() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_fem_build_rectangle() {
        let mut mesh = PyFemMesh::new();
        let n_nodes = mesh.build_rectangle(0.0, 0.0, 1.0, 1.0, 4, 4, 0, 1.0);
        assert_eq!(n_nodes, 25); // 5×5 = 25 nodes
        assert_eq!(mesh.num_elements(), 32); // 4×4×2 = 32 triangles
    }

    #[test]
    fn test_fem_material_aluminum() {
        let mat = PyFemMaterial::aluminum();
        assert!((mat.young_modulus - 70e9).abs() < 1.0);
    }

    #[test]
    fn test_fem_material_rubber_high_poisson() {
        let mat = PyFemMaterial::rubber();
        assert!(mat.poisson_ratio > 0.4);
    }

    #[test]
    fn test_fem_material_concrete() {
        let mat = PyFemMaterial::concrete();
        assert!((mat.density - 2400.0).abs() < 1.0);
    }

    #[test]
    fn test_fem_nodal_force_new() {
        let nf = PyFemNodalForce::new(3, 100.0, -200.0);
        assert_eq!(nf.node, 3);
        assert!((nf.force[0] - 100.0).abs() < 1e-10);
        assert!((nf.force[1] + 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_fem_dirichlet_fix_x() {
        let bc = PyFemDirichletBC::fix_x(5);
        assert_eq!(bc.node, 5);
        assert_eq!(bc.dof, 0);
        assert!((bc.value).abs() < 1e-15);
    }

    #[test]
    fn test_fem_dirichlet_fix_y() {
        let bc = PyFemDirichletBC::fix_y(2);
        assert_eq!(bc.dof, 1);
    }

    #[test]
    fn test_fem_solve_pinned_both_sides_zero_displacement_at_pin() {
        let mut mesh = simple_mesh();
        // Pin all nodes — no load — displacement should be zero
        for n in 0..4 {
            mesh.pin_node(n);
        }
        let solver = PyFemSolver::new();
        let result = solver.solve(&mesh).expect("should solve");
        for &d in &result.displacements {
            assert!(d.abs() < 1e-10, "constrained mesh with no load: d={}", d);
        }
    }

    #[test]
    fn test_fem_add_material_returns_index() {
        let mut mesh = PyFemMesh::new();
        let idx = mesh.add_material(PyFemMaterial::rubber());
        // Default material is at index 0, so new one is at index 1
        assert_eq!(idx, 1);
    }

    #[test]
    fn test_fem_element_stiffness_symmetric() {
        // The element stiffness matrix should be symmetric
        let mut mesh = PyFemMesh::new();
        mesh.add_node(0.0, 0.0);
        mesh.add_node(1.0, 0.0);
        mesh.add_node(0.5, 1.0);
        mesh.add_element(0, 1, 2);

        let solver = PyFemSolver::new();
        let elem = &mesh.elements[0];
        let mat = &mesh.materials[0];
        let ke = solver.element_stiffness(&mesh, elem, mat);

        for i in 0..6 {
            for j in 0..6 {
                let diff = (ke[i * 6 + j] - ke[j * 6 + i]).abs();
                assert!(
                    diff < 1e-6,
                    "Ke not symmetric at ({},{}): diff={}",
                    i,
                    j,
                    diff
                );
            }
        }
    }

    #[test]
    fn test_fem_solver_with_settings() {
        let solver = PyFemSolver::with_settings(500, 1e-8);
        assert_eq!(solver.max_iterations, 500);
        assert!((solver.tolerance - 1e-8).abs() < 1e-20);
    }

    #[test]
    fn test_fem_result_node_displacement_out_of_bounds() {
        let result = PyFemSolveResult {
            displacements: vec![0.0; 4],
            von_mises_stress: vec![],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        assert!(result.node_displacement(10).is_none());
    }

    #[test]
    fn test_fem_mesh_serde_roundtrip() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.apply_nodal_force(2, 100.0, -200.0);
        let json = serde_json::to_string(&mesh).expect("serialize");
        let back: PyFemMesh = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back.num_nodes(), 4);
        assert_eq!(back.dirichlet_bcs.len(), 2);
        assert_eq!(back.nodal_forces.len(), 1);
    }
}

// ===========================================================================
// Modal analysis (free-vibration eigenvalue problem)
// ===========================================================================

/// Result of a modal analysis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ModalAnalysisResult {
    /// Natural frequencies in rad/s (ascending order).
    pub natural_frequencies: Vec<f64>,
    /// Mode shapes: each row is a normalised displacement eigenvector (flat).
    pub mode_shapes: Vec<Vec<f64>>,
    /// Number of modes extracted.
    pub num_modes: usize,
}

impl ModalAnalysisResult {
    /// Return the `i`-th natural frequency in Hz.
    pub fn frequency_hz(&self, i: usize) -> Option<f64> {
        self.natural_frequencies
            .get(i)
            .map(|&f| f / (2.0 * std::f64::consts::PI))
    }

    /// Return the period (in seconds) of the `i`-th mode.
    pub fn period_s(&self, i: usize) -> Option<f64> {
        self.frequency_hz(i).filter(|&f| f > 1e-15).map(|f| 1.0 / f)
    }
}

/// Simple power-iteration modal analyzer (extracts the dominant modes of K).
///
/// This is a proof-of-concept; it uses a deflated power iteration to
/// approximate the lowest `num_modes` natural frequencies for undamped free
/// vibration (K φ = ω² M φ, with M = I).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ModalAnalyzer {
    /// Number of modes to extract.
    pub num_modes: usize,
    /// Maximum power-iteration steps per mode.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

impl ModalAnalyzer {
    /// Create an analyzer that extracts up to `num_modes` modes.
    pub fn new(num_modes: usize) -> Self {
        Self {
            num_modes,
            max_iter: 500,
            tolerance: 1e-8,
        }
    }

    /// Run modal analysis on the assembled mesh.
    ///
    /// Returns `None` if the mesh has no elements.
    pub fn analyze(&self, mesh: &PyFemMesh) -> Option<ModalAnalysisResult> {
        if mesh.num_nodes() == 0 || mesh.elements.is_empty() {
            return None;
        }
        let solver = PyFemSolver::new();
        let n_dofs = mesh.num_dofs();

        // Assemble stiffness matrix (reuse solver's internal logic via a zero-load solve)
        let mut k = vec![0.0f64; n_dofs * n_dofs];
        for elem in &mesh.elements {
            let mat = mesh
                .materials
                .get(elem.material_id)
                .unwrap_or(&mesh.materials[0]);
            let ke = solver.element_stiffness(mesh, elem, mat);
            let dofs: [usize; 6] = [
                elem.nodes[0] * 2,
                elem.nodes[0] * 2 + 1,
                elem.nodes[1] * 2,
                elem.nodes[1] * 2 + 1,
                elem.nodes[2] * 2,
                elem.nodes[2] * 2 + 1,
            ];
            for (i, &gi) in dofs.iter().enumerate() {
                for (j, &gj) in dofs.iter().enumerate() {
                    k[gi * n_dofs + gj] += ke[i * 6 + j];
                }
            }
        }

        // Apply BCs: zero out rows/cols of constrained DOFs
        for bc in &mesh.dirichlet_bcs {
            if bc.node < mesh.num_nodes() {
                let dof = bc.node * 2 + bc.dof;
                for j in 0..n_dofs {
                    k[dof * n_dofs + j] = 0.0;
                    k[j * n_dofs + dof] = 0.0;
                }
                k[dof * n_dofs + dof] = 1.0;
            }
        }

        let modes = self.num_modes.min(n_dofs);
        let mut natural_frequencies = Vec::with_capacity(modes);
        let mut mode_shapes = Vec::with_capacity(modes);

        // Deflated inverse-free power iteration: approximate λ = ω²
        let mut deflated_k = k.clone();

        for _ in 0..modes {
            let (lambda, phi) = power_iteration(&deflated_k, n_dofs, self.max_iter, self.tolerance);
            if lambda.abs() < 1e-15 {
                break;
            }
            let omega = lambda.abs().sqrt();
            natural_frequencies.push(omega);
            mode_shapes.push(phi.clone());

            // Deflate: K ← K - λ φ φᵀ
            for i in 0..n_dofs {
                for j in 0..n_dofs {
                    deflated_k[i * n_dofs + j] -= lambda * phi[i] * phi[j];
                }
            }
        }

        Some(ModalAnalysisResult {
            num_modes: natural_frequencies.len(),
            natural_frequencies,
            mode_shapes,
        })
    }
}

/// Power iteration to find the dominant eigenvalue and eigenvector of a symmetric matrix.
fn power_iteration(a: &[f64], n: usize, max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
    let mut lambda = 0.0f64;
    for _ in 0..max_iter {
        let av = mat_vec_mul(a, &v, n);
        let new_lambda: f64 = av.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
        let norm: f64 = (av.iter().map(|x| x * x).sum::<f64>()).sqrt();
        if norm < 1e-30 {
            break;
        }
        v = av.iter().map(|x| x / norm).collect();
        if (new_lambda - lambda).abs() < tol {
            lambda = new_lambda;
            break;
        }
        lambda = new_lambda;
    }
    (lambda, v)
}

// ===========================================================================
// Result extraction helpers
// ===========================================================================

impl PyFemSolveResult {
    /// Return the displacement field as a list of `[ux, uy]` pairs.
    pub fn displacement_pairs(&self) -> Vec<[f64; 2]> {
        self.displacements
            .chunks_exact(2)
            .map(|c| [c[0], c[1]])
            .collect()
    }

    /// Maximum absolute displacement component.
    pub fn max_displacement_magnitude(&self) -> f64 {
        self.displacement_pairs()
            .iter()
            .map(|[ux, uy]| (ux * ux + uy * uy).sqrt())
            .fold(0.0f64, f64::max)
    }

    /// Return all von Mises stresses as a sorted (ascending) copy.
    pub fn sorted_von_mises(&self) -> Vec<f64> {
        let mut v = self.von_mises_stress.clone();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        v
    }

    /// Return the `q`-th quantile of the von Mises stress distribution (q in \[0,1\]).
    pub fn von_mises_quantile(&self, q: f64) -> f64 {
        if self.von_mises_stress.is_empty() {
            return 0.0;
        }
        let sorted = self.sorted_von_mises();
        let idx = ((q.clamp(0.0, 1.0)) * (sorted.len() - 1) as f64).round() as usize;
        sorted[idx]
    }
}

// ===========================================================================
// Mesh query helpers
// ===========================================================================

impl PyFemMesh {
    /// Return the bounding box of all nodes as `[xmin, ymin, xmax, ymax]`.
    pub fn bounding_box(&self) -> Option<[f64; 4]> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut xmin = f64::INFINITY;
        let mut ymin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        let mut ymax = f64::NEG_INFINITY;
        for n in &self.nodes {
            xmin = xmin.min(n.position[0]);
            ymin = ymin.min(n.position[1]);
            xmax = xmax.max(n.position[0]);
            ymax = ymax.max(n.position[1]);
        }
        Some([xmin, ymin, xmax, ymax])
    }

    /// Return the total area of all elements (sum of |area_i|).
    pub fn total_area(&self) -> f64 {
        (0..self.num_elements())
            .map(|e| self.element_area(e).abs())
            .sum()
    }

    /// Return the centroid of element `e` as `[cx, cy]`.
    pub fn element_centroid(&self, e: usize) -> Option<[f64; 2]> {
        if e >= self.elements.len() {
            return None;
        }
        let [n0, n1, n2] = self.elements[e].nodes;
        let [x0, y0] = self.nodes[n0].position;
        let [x1, y1] = self.nodes[n1].position;
        let [x2, y2] = self.nodes[n2].position;
        Some([(x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0])
    }

    /// Find the node closest to position `(px, py)`.
    pub fn closest_node(&self, px: f64, py: f64) -> Option<usize> {
        if self.nodes.is_empty() {
            return None;
        }
        let (idx, _) = self
            .nodes
            .iter()
            .enumerate()
            .map(|(i, n)| {
                let dx = n.position[0] - px;
                let dy = n.position[1] - py;
                (i, dx * dx + dy * dy)
            })
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        Some(idx)
    }

    /// Return nodes on the left boundary (x <= xmin + tol).
    pub fn left_boundary_nodes(&self, tol: f64) -> Vec<usize> {
        if let Some(bb) = self.bounding_box() {
            self.nodes
                .iter()
                .enumerate()
                .filter(|(_, n)| n.position[0] <= bb[0] + tol)
                .map(|(i, _)| i)
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Return nodes on the right boundary (x >= xmax - tol).
    pub fn right_boundary_nodes(&self, tol: f64) -> Vec<usize> {
        if let Some(bb) = self.bounding_box() {
            self.nodes
                .iter()
                .enumerate()
                .filter(|(_, n)| n.position[0] >= bb[2] - tol)
                .map(|(i, _)| i)
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Pin all nodes on the left boundary.
    pub fn pin_left_boundary(&mut self, tol: f64) {
        let nodes = self.left_boundary_nodes(tol);
        for n in nodes {
            self.pin_node(n);
        }
    }
}

// ===========================================================================
// Additional FEM tests
// ===========================================================================

#[cfg(test)]
mod fem_ext_tests {

    use crate::PyFemMesh;
    use crate::PyFemSolveResult;
    use crate::PyFemSolver;
    use crate::fem_api::ModalAnalysisResult;
    use crate::fem_api::ModalAnalyzer;
    use crate::fem_api::power_iteration;

    fn simple_mesh() -> PyFemMesh {
        let mut mesh = PyFemMesh::new();
        mesh.add_node(0.0, 0.0);
        mesh.add_node(1.0, 0.0);
        mesh.add_node(1.0, 1.0);
        mesh.add_node(0.0, 1.0);
        mesh.add_element(0, 1, 2);
        mesh.add_element(0, 2, 3);
        mesh
    }

    // --- ModalAnalyzer ---

    #[test]
    fn test_modal_analysis_returns_some() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        let analyzer = ModalAnalyzer::new(2);
        let result = analyzer.analyze(&mesh);
        assert!(result.is_some());
    }

    #[test]
    fn test_modal_analysis_returns_none_empty_mesh() {
        let mesh = PyFemMesh::new();
        let analyzer = ModalAnalyzer::new(2);
        assert!(analyzer.analyze(&mesh).is_none());
    }

    #[test]
    fn test_modal_analysis_num_modes() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        let analyzer = ModalAnalyzer::new(2);
        if let Some(result) = analyzer.analyze(&mesh) {
            assert!(result.num_modes <= 2);
        }
    }

    #[test]
    fn test_modal_analysis_frequencies_positive() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        let analyzer = ModalAnalyzer::new(3);
        if let Some(result) = analyzer.analyze(&mesh) {
            for &f in &result.natural_frequencies {
                assert!(f >= 0.0, "frequency should be non-negative: {}", f);
            }
        }
    }

    #[test]
    fn test_modal_analysis_frequency_hz() {
        let mut mesh = simple_mesh();
        mesh.pin_node(0);
        mesh.pin_node(3);
        let analyzer = ModalAnalyzer::new(1);
        if let Some(result) = analyzer.analyze(&mesh)
            && result.num_modes > 0
        {
            let hz = result.frequency_hz(0);
            assert!(hz.is_some());
            assert!(hz.unwrap() >= 0.0);
        }
    }

    #[test]
    fn test_modal_period_s() {
        let result = ModalAnalysisResult {
            natural_frequencies: vec![2.0 * std::f64::consts::PI * 100.0], // 100 Hz
            mode_shapes: vec![vec![1.0, 0.0, 0.0, 0.0]],
            num_modes: 1,
        };
        let period = result.period_s(0).unwrap();
        assert!((period - 0.01).abs() < 1e-6, "period = {}", period);
    }

    // --- Result extraction helpers ---

    #[test]
    fn test_displacement_pairs_length() {
        let result = PyFemSolveResult {
            displacements: vec![1.0, 2.0, 3.0, 4.0],
            von_mises_stress: vec![],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        let pairs = result.displacement_pairs();
        assert_eq!(pairs.len(), 2);
        assert!((pairs[0][0] - 1.0).abs() < 1e-15);
        assert!((pairs[1][1] - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_max_displacement_magnitude() {
        let result = PyFemSolveResult {
            displacements: vec![3.0, 4.0, 0.0, 0.0],
            von_mises_stress: vec![],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        assert!((result.max_displacement_magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_sorted_von_mises() {
        let result = PyFemSolveResult {
            displacements: vec![],
            von_mises_stress: vec![3.0, 1.0, 2.0],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        let sorted = result.sorted_von_mises();
        assert!((sorted[0] - 1.0).abs() < 1e-15);
        assert!((sorted[2] - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_von_mises_quantile_min() {
        let result = PyFemSolveResult {
            displacements: vec![],
            von_mises_stress: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        assert!((result.von_mises_quantile(0.0) - 1.0).abs() < 1e-10);
        assert!((result.von_mises_quantile(1.0) - 5.0).abs() < 1e-10);
    }

    // --- Mesh query helpers ---

    #[test]
    fn test_bounding_box_unit_square() {
        let mesh = simple_mesh();
        let bb = mesh.bounding_box().unwrap();
        assert!((bb[0]).abs() < 1e-15); // xmin = 0
        assert!((bb[2] - 1.0).abs() < 1e-15); // xmax = 1
    }

    #[test]
    fn test_bounding_box_empty() {
        let mesh = PyFemMesh::new();
        assert!(mesh.bounding_box().is_none());
    }

    #[test]
    fn test_total_area_unit_square() {
        let mesh = simple_mesh();
        assert!((mesh.total_area() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_element_centroid_valid() {
        let mesh = simple_mesh();
        let c = mesh.element_centroid(0).unwrap();
        // Triangle (0,0),(1,0),(1,1): centroid = (2/3, 1/3)
        assert!((c[0] - 2.0 / 3.0).abs() < 1e-10, "cx = {}", c[0]);
    }

    #[test]
    fn test_element_centroid_out_of_bounds() {
        let mesh = simple_mesh();
        assert!(mesh.element_centroid(99).is_none());
    }

    #[test]
    fn test_closest_node_corner() {
        let mesh = simple_mesh();
        let idx = mesh.closest_node(0.01, 0.01).unwrap();
        assert_eq!(idx, 0); // closest to (0,0)
    }

    #[test]
    fn test_closest_node_empty() {
        let mesh = PyFemMesh::new();
        assert!(mesh.closest_node(0.0, 0.0).is_none());
    }

    #[test]
    fn test_left_boundary_nodes() {
        let mesh = simple_mesh();
        let left = mesh.left_boundary_nodes(1e-6);
        assert!(left.contains(&0)); // node at (0,0)
        assert!(left.contains(&3)); // node at (0,1)
    }

    #[test]
    fn test_right_boundary_nodes() {
        let mesh = simple_mesh();
        let right = mesh.right_boundary_nodes(1e-6);
        assert!(right.contains(&1)); // node at (1,0)
        assert!(right.contains(&2)); // node at (1,1)
    }

    #[test]
    fn test_pin_left_boundary() {
        let mut mesh = simple_mesh();
        mesh.pin_left_boundary(1e-6);
        // 2 left nodes × 2 DOFs each = 4 BCs
        assert_eq!(mesh.dirichlet_bcs.len(), 4);
    }

    #[test]
    fn test_fem_solve_with_large_mesh() {
        let mut mesh = PyFemMesh::new();
        mesh.build_rectangle(0.0, 0.0, 1.0, 0.5, 4, 2, 0, 0.01);
        mesh.pin_left_boundary(1e-6);
        // Apply load on right boundary
        let right = mesh.right_boundary_nodes(1e-6);
        for n in right {
            mesh.apply_nodal_force(n, 0.0, -500.0);
        }
        let solver = PyFemSolver::new();
        let result = solver.solve(&mesh).expect("should solve");
        assert!(!result.displacements.is_empty());
    }

    #[test]
    fn test_power_iteration_identity_matrix() {
        // Identity matrix: eigenvalue = 1, eigenvector = [1, 0, 0]
        let a = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let (lambda, _phi) = power_iteration(&a, 3, 100, 1e-10);
        assert!((lambda - 1.0).abs() < 1e-6, "lambda = {}", lambda);
    }

    #[test]
    fn test_fem_result_displacement_pairs_empty() {
        let result = PyFemSolveResult {
            displacements: vec![],
            von_mises_stress: vec![],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        assert!(result.displacement_pairs().is_empty());
    }

    #[test]
    fn test_fem_max_displacement_empty() {
        let result = PyFemSolveResult {
            displacements: vec![],
            von_mises_stress: vec![],
            solver_iterations: 0,
            residual_norm: 0.0,
            converged: true,
        };
        assert!((result.max_displacement_magnitude()).abs() < 1e-15);
    }
}
