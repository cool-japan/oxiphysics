// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Static, modal, and dynamic FEM analysis drivers.
//!
//! Provides:
//! - [`LinearStaticAnalysis`] — classic K·u = f solve
//! - [`StaticAnalysis`] — standalone dense-matrix static analysis
//! - [`ModalAnalysis`] — natural frequencies and mode shapes of `K·φ = ω²·M·φ`
//!   (Cholesky reduction + Jacobi, shifted inverse iteration for single modes)
//! - [`DynamicAnalysis`] — Newmark-β and Wilson-θ time integration

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use oxiphysics_core::math::Vec3;

use crate::assembly::{assemble_load_vector, assemble_stiffness};
use crate::boundary::{DirichletBc, NeumannBc, apply_dirichlet, apply_neumann};
use crate::constitutive::LinearElasticMaterial;
use crate::element::LinearTetrahedron;
use crate::error::Error;
use crate::mesh::TetrahedralMesh;
use crate::solvers::{jacobi_preconditioner, pcg_solve};

pub(crate) mod dense;

// ---------------------------------------------------------------------------
// FEM error type
// ---------------------------------------------------------------------------

/// Errors that can arise during FEM analysis.
#[derive(Debug, Clone, PartialEq)]
pub enum FemError {
    /// The linear system is singular or ill-conditioned.
    SingularSystem,
    /// An input argument has an invalid size or value.
    InvalidInput(String),
    /// The iterative solver did not converge within the allowed iterations.
    SolverNotConverged,
}

impl std::fmt::Display for FemError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FemError::SingularSystem => write!(f, "singular or ill-conditioned system"),
            FemError::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            FemError::SolverNotConverged => write!(f, "iterative solver did not converge"),
        }
    }
}

impl std::error::Error for FemError {}

impl From<Error> for FemError {
    /// Map the crate-level [`Error`] onto the analysis error type so that the
    /// fallible helpers (boundary-condition application, ...) can be used with `?`
    /// inside functions that return [`FemError`].
    fn from(err: Error) -> Self {
        match err {
            Error::SingularSystem => FemError::SingularSystem,
            Error::NotConverged { .. } => FemError::SolverNotConverged,
            Error::InvalidInput(message) | Error::General(message) => {
                FemError::InvalidInput(message)
            }
            other => FemError::InvalidInput(other.to_string()),
        }
    }
}

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

/// Result of a finite element analysis.
#[derive(Debug, Clone)]
pub struct FemResult {
    /// Nodal displacement vectors.
    pub displacements: Vec<Vec3>,
    /// Stress tensor per element (Voigt: \[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz\]).
    pub stress: Vec<[f64; 6]>,
    /// Strain tensor per element (Voigt: \[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz\]).
    pub strain: Vec<[f64; 6]>,
}

// ---------------------------------------------------------------------------
// LinearStaticAnalysis (existing, preserved)
// ---------------------------------------------------------------------------

/// Linear static analysis driver.
///
/// Assembles the global stiffness matrix and load vector, applies boundary
/// conditions, solves the linear system, and post-processes stresses and strains.
/// Use [`try_solve`](Self::try_solve) to run it: it validates the problem and
/// reports a solve that did not converge instead of returning garbage.
pub struct LinearStaticAnalysis {
    /// Maximum iterations for the iterative solver.
    pub max_iter: usize,
    /// Convergence tolerance for the iterative solver.
    pub tolerance: f64,
}

impl Default for LinearStaticAnalysis {
    fn default() -> Self {
        Self {
            max_iter: 10000,
            tolerance: 1e-10,
        }
    }
}

impl LinearStaticAnalysis {
    /// Create a new linear static analysis with default solver parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Run the analysis, validating the input and the numerical result.
    ///
    /// Assembles the global stiffness matrix and load vector, applies the
    /// boundary conditions, solves `K·u = f` with Jacobi-preconditioned conjugate
    /// gradients and recovers the element strains and stresses.
    ///
    /// # Errors
    ///
    /// - [`FemError::InvalidInput`] when the problem is ill-posed before any
    ///   solve is attempted: an empty mesh, an element that references a missing
    ///   node, a degenerate element (`|det J| < 1e-12·h³`, `h` being its longest
    ///   edge), a Young's modulus that is not positive and finite, a Poisson ratio
    ///   outside `(-1, ½)`, non-finite coordinates / loads / prescribed values /
    ///   solver tolerance, a Dirichlet BC whose node or local DOF (0, 1, 2) does
    ///   not exist, conflicting prescribed values for one DOF, or a Neumann BC on
    ///   a missing node.
    /// - [`FemError::SolverNotConverged`] when the conjugate-gradient iteration
    ///   does not reach the relative residual `tolerance` within `max_iter`
    ///   iterations — typically because the structure is not constrained against
    ///   rigid-body motion, so that `K` is singular.
    /// - [`FemError::SingularSystem`] when the displacement, strain or stress
    ///   field contains a non-finite number.
    pub fn try_solve(
        &self,
        mesh: &TetrahedralMesh,
        material: &LinearElasticMaterial,
        dirichlet_bcs: &[DirichletBc],
        neumann_bcs: &[NeumannBc],
        body_force: &Vec3,
    ) -> Result<FemResult, FemError> {
        self.validate_input(mesh, material, dirichlet_bcs, neumann_bcs, body_force)?;

        let mut k = assemble_stiffness(mesh, material);
        let mut f = assemble_load_vector(mesh, body_force);

        apply_neumann(&mut f, neumann_bcs);
        apply_dirichlet(&mut k, &mut f, dirichlet_bcs)?;

        let precond = jacobi_preconditioner(&k);
        let (u, stats) = pcg_solve(&k, &f, &precond, self.max_iter, self.tolerance);
        if !stats.converged {
            return Err(FemError::SolverNotConverged);
        }
        if u.iter().any(|value| !value.is_finite()) {
            return Err(FemError::SingularSystem);
        }

        let num_nodes = mesh.num_nodes();
        let mut displacements = Vec::with_capacity(num_nodes);
        for i in 0..num_nodes {
            displacements.push(Vec3::new(u[i * 3], u[i * 3 + 1], u[i * 3 + 2]));
        }

        let d_matrix = material.constitutive_matrix();
        let num_elements = mesh.num_elements();
        let mut stress = Vec::with_capacity(num_elements);
        let mut strain = Vec::with_capacity(num_elements);

        for elem in &mesh.elements {
            let nodes = [
                mesh.nodes[elem[0]],
                mesh.nodes[elem[1]],
                mesh.nodes[elem[2]],
                mesh.nodes[elem[3]],
            ];

            let (b_mat, _vol) = LinearTetrahedron::b_matrix(&nodes);

            let mut u_e = [0.0f64; 12];
            for (local, &global) in elem.iter().enumerate() {
                u_e[local * 3] = u[global * 3];
                u_e[local * 3 + 1] = u[global * 3 + 1];
                u_e[local * 3 + 2] = u[global * 3 + 2];
            }

            let mut eps = [0.0f64; 6];
            for i in 0..6 {
                for j in 0..12 {
                    eps[i] += b_mat[i * 12 + j] * u_e[j];
                }
            }

            let mut sig = [0.0f64; 6];
            for i in 0..6 {
                for j in 0..6 {
                    sig[i] += d_matrix[i][j] * eps[j];
                }
            }

            strain.push(eps);
            stress.push(sig);
        }

        let all_finite = strain
            .iter()
            .chain(stress.iter())
            .all(|voigt| voigt.iter().all(|component| component.is_finite()));
        if !all_finite {
            return Err(FemError::SingularSystem);
        }

        Ok(FemResult {
            displacements,
            stress,
            strain,
        })
    }

    /// Run the analysis and return the results, without reporting failures.
    ///
    /// This is the infallible convenience form of [`try_solve`](Self::try_solve)
    /// and delegates to it.  When `try_solve` fails (invalid input, a CG solve
    /// that did not converge, a non-finite result, ...) every displacement,
    /// strain and stress component of the returned [`FemResult`] is `NaN`, so a
    /// failed analysis can never be mistaken for a valid one.  Prefer
    /// [`try_solve`](Self::try_solve): it tells you *why* the analysis failed.
    pub fn solve(
        &self,
        mesh: &TetrahedralMesh,
        material: &LinearElasticMaterial,
        dirichlet_bcs: &[DirichletBc],
        neumann_bcs: &[NeumannBc],
        body_force: &Vec3,
    ) -> FemResult {
        self.try_solve(mesh, material, dirichlet_bcs, neumann_bcs, body_force)
            .unwrap_or_else(|_| FemResult::not_a_number(mesh.num_nodes(), mesh.num_elements()))
    }

    /// Check everything that would otherwise panic or silently corrupt the
    /// system before it is assembled.
    fn validate_input(
        &self,
        mesh: &TetrahedralMesh,
        material: &LinearElasticMaterial,
        dirichlet_bcs: &[DirichletBc],
        neumann_bcs: &[NeumannBc],
        body_force: &Vec3,
    ) -> Result<(), FemError> {
        let invalid = |message: String| Err(FemError::InvalidInput(message));
        let num_nodes = mesh.num_nodes();
        if num_nodes == 0 || mesh.num_elements() == 0 {
            return invalid(format!(
                "mesh must contain nodes and elements (nodes: {num_nodes}, elements: {})",
                mesh.num_elements()
            ));
        }
        if !(self.tolerance.is_finite() && self.tolerance > 0.0) {
            return invalid(format!(
                "solver tolerance must be finite and positive, got {}",
                self.tolerance
            ));
        }
        let young = material.young_modulus;
        if !(young.is_finite() && young > 0.0) {
            return invalid(format!(
                "Young's modulus must be finite and positive, got {young}"
            ));
        }
        let poisson = material.poisson_ratio;
        if !(poisson.is_finite() && poisson > -1.0 && poisson < 0.5) {
            return invalid(format!(
                "Poisson's ratio must lie in (-1, 0.5), got {poisson}"
            ));
        }
        if !body_force.iter().all(|c| c.is_finite()) {
            return invalid("body force must be finite".to_string());
        }
        for (index, node) in mesh.nodes.iter().enumerate() {
            if !node.iter().all(|c| c.is_finite()) {
                return invalid(format!("node {index} has non-finite coordinates"));
            }
        }
        for (index, elem) in mesh.elements.iter().enumerate() {
            if let Some(&missing) = elem.iter().find(|&&node| node >= num_nodes) {
                return invalid(format!(
                    "element {index} references node {missing} but the mesh has {num_nodes} nodes"
                ));
            }
            let corner = [
                mesh.nodes[elem[0]],
                mesh.nodes[elem[1]],
                mesh.nodes[elem[2]],
                mesh.nodes[elem[3]],
            ];
            if Self::is_degenerate_tetrahedron(&corner) {
                return invalid(format!(
                    "element {index} is degenerate (zero or numerically zero volume)"
                ));
            }
        }
        for (index, bc) in dirichlet_bcs.iter().enumerate() {
            if bc.node >= num_nodes || bc.dof >= 3 {
                return invalid(format!(
                    "Dirichlet BC #{index} addresses node {} DOF {}, but the mesh has {num_nodes} \
                     nodes with DOFs 0..3",
                    bc.node, bc.dof
                ));
            }
            if !bc.value.is_finite() {
                return invalid(format!("Dirichlet BC #{index} has a non-finite value"));
            }
        }
        for (index, bc) in neumann_bcs.iter().enumerate() {
            if bc.node >= num_nodes {
                return invalid(format!(
                    "Neumann BC #{index} addresses node {} but the mesh has {num_nodes} nodes",
                    bc.node
                ));
            }
            if !bc.force.iter().all(|c| c.is_finite()) {
                return invalid(format!("Neumann BC #{index} has a non-finite force"));
            }
        }
        Ok(())
    }

    /// `true` when `|det J| < 1e-12·h³` (`h` = longest edge) or the geometry is
    /// not finite, i.e. the B-matrix would divide by a (numerically) zero volume.
    fn is_degenerate_tetrahedron(nodes: &[Vec3; 4]) -> bool {
        let x10 = nodes[1] - nodes[0];
        let x20 = nodes[2] - nodes[0];
        let x30 = nodes[3] - nodes[0];
        let det_j = x10.dot(&x20.cross(&x30));
        let longest_edge = [
            x10,
            x20,
            x30,
            nodes[2] - nodes[1],
            nodes[3] - nodes[1],
            nodes[3] - nodes[2],
        ]
        .iter()
        .map(|edge| edge.norm())
        .fold(0.0_f64, f64::max);
        // A NaN determinant or edge length makes the comparison false, so a
        // non-finite geometry is reported as degenerate as well.
        let well_shaped = det_j.abs() > 1e-12 * longest_edge.powi(3);
        !well_shaped
    }
}

impl FemResult {
    /// A result whose every component is `NaN`, used to flag a failed analysis.
    fn not_a_number(num_nodes: usize, num_elements: usize) -> Self {
        Self {
            displacements: vec![Vec3::new(f64::NAN, f64::NAN, f64::NAN); num_nodes],
            stress: vec![[f64::NAN; 6]; num_elements],
            strain: vec![[f64::NAN; 6]; num_elements],
        }
    }
}

// ---------------------------------------------------------------------------
// StaticAnalysis — standalone dense-matrix static solver
// ---------------------------------------------------------------------------

/// Standalone static analysis using a dense global stiffness matrix.
///
/// Suitable for small problems where CSR overhead is not needed.
pub struct StaticAnalysis {
    /// Number of degrees of freedom.
    pub n_dofs: usize,
    /// Dense global stiffness matrix stored row-major (n_dofs × n_dofs).
    pub stiffness_matrix: Vec<f64>,
    /// Global force vector (length n_dofs).
    pub force_vector: Vec<f64>,
}

impl StaticAnalysis {
    /// Create a new `StaticAnalysis` with zeros.
    pub fn new(n_dofs: usize) -> Self {
        Self {
            n_dofs,
            stiffness_matrix: vec![0.0; n_dofs * n_dofs],
            force_vector: vec![0.0; n_dofs],
        }
    }

    /// Assemble the global stiffness matrix from element stiffness matrices.
    ///
    /// `elements` is a slice of `(dof_indices, ke)` pairs where `ke` is stored
    /// row-major for the element's local DOFs.
    pub fn assemble_global_stiffness(&mut self, elements: &[(Vec<usize>, Vec<f64>)]) {
        let n = self.n_dofs;
        for (dofs, ke) in elements {
            let ndof_e = dofs.len();
            for i in 0..ndof_e {
                for j in 0..ndof_e {
                    let gi = dofs[i];
                    let gj = dofs[j];
                    self.stiffness_matrix[gi * n + gj] += ke[i * ndof_e + j];
                }
            }
        }
    }

    /// Solve the linear system K·u = f using Gaussian elimination with partial pivoting.
    ///
    /// Returns the displacement vector or a [`FemError`] if the system is singular.
    pub fn solve(&self) -> Result<Vec<f64>, FemError> {
        let n = self.n_dofs;
        if n == 0 {
            return Ok(Vec::new());
        }
        // Augmented matrix [K | f]
        let mut aug: Vec<f64> = Vec::with_capacity(n * (n + 1));
        for i in 0..n {
            for j in 0..n {
                aug.push(self.stiffness_matrix[i * n + j]);
            }
            aug.push(self.force_vector[i]);
        }

        let cols = n + 1;

        for col in 0..n {
            // Find pivot row
            let mut max_val = aug[col * cols + col].abs();
            let mut max_row = col;
            for row in (col + 1)..n {
                let v = aug[row * cols + col].abs();
                if v > max_val {
                    max_val = v;
                    max_row = row;
                }
            }
            if max_val < 1e-14 {
                return Err(FemError::SingularSystem);
            }
            // Swap rows
            if max_row != col {
                for k in 0..cols {
                    aug.swap(col * cols + k, max_row * cols + k);
                }
            }
            let pivot = aug[col * cols + col];
            // Eliminate below
            for row in (col + 1)..n {
                let factor = aug[row * cols + col] / pivot;
                for k in col..cols {
                    let sub = factor * aug[col * cols + k];
                    aug[row * cols + k] -= sub;
                }
            }
        }

        // Back substitution
        let mut u = vec![0.0; n];
        for i in (0..n).rev() {
            let mut sum = aug[i * cols + n];
            for j in (i + 1)..n {
                sum -= aug[i * cols + j] * u[j];
            }
            u[i] = sum / aug[i * cols + i];
        }
        Ok(u)
    }

    /// Compute engineering strains ε = B·u for each element.
    ///
    /// `elements` is a slice of `(dof_indices, B_matrix)` where `B` is 6×n_dof_e
    /// stored row-major. Returns one `[f64;6]` per element.
    pub fn compute_strains(
        &self,
        displacements: &[f64],
        elements: &[(Vec<usize>, Vec<f64>)],
    ) -> Vec<[f64; 6]> {
        let mut strains = Vec::with_capacity(elements.len());
        for (dofs, b_mat) in elements {
            let ndof_e = dofs.len();
            let mut u_e = vec![0.0; ndof_e];
            for (k, &g) in dofs.iter().enumerate() {
                u_e[k] = displacements[g];
            }
            let mut eps = [0.0f64; 6];
            for i in 0..6 {
                for j in 0..ndof_e {
                    eps[i] += b_mat[i * ndof_e + j] * u_e[j];
                }
            }
            strains.push(eps);
        }
        strains
    }

    /// Compute stresses σ = D·ε from the constitutive matrix and strains.
    ///
    /// `d_matrix` is the 6×6 constitutive matrix (row-major flat).
    pub fn compute_stresses(strains: &[[f64; 6]], d_matrix: &[f64; 36]) -> Vec<[f64; 6]> {
        strains
            .iter()
            .map(|eps| {
                let mut sig = [0.0f64; 6];
                for i in 0..6 {
                    for j in 0..6 {
                        sig[i] += d_matrix[i * 6 + j] * eps[j];
                    }
                }
                sig
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// ModalAnalysis — generalized symmetric eigenproblem K·φ = ω²·M·φ
// ---------------------------------------------------------------------------

/// Convergence tolerance (change of the normalised iterate) of the shifted
/// inverse iteration behind [`ModalAnalysis::compute_mode_shapes`].
const MODE_SHAPE_TOLERANCE: f64 = 1.0e-10;

/// Relative size of the first shift perturbation applied by the mode-shape
/// iteration when `K − ω²·M` is singular (the shift sits exactly on an
/// eigenvalue); it is grown a hundredfold per retry.
const SHIFT_PERTURBATION: f64 = 1.0e-10;

/// Number of perturbed retries before `K − ω²·M` is declared singular.
const SHIFT_RETRIES: usize = 3;

/// A negative eigenvalue larger in magnitude than this fraction of the
/// spectral radius is a genuine instability, not round-off of a rigid-body mode.
const NEGATIVE_EIGENVALUE_TOLERANCE: f64 = 1.0e-10;

/// Modal analysis of the dense undamped system `K·φ = ω²·M·φ`.
///
/// Natural frequencies and mode shapes come from the exact generalized
/// symmetric-definite eigen-decomposition: `M = L·Lᵀ` (Cholesky), the reduced
/// matrix `C = L⁻¹·K·L⁻ᵀ` is diagonalised by cyclic Jacobi sweeps and the mode
/// shapes follow from `φ = L⁻ᵀ·y`.  The `n_modes` lowest modes are therefore
/// returned **in ascending order of frequency** and none can be skipped —
/// something a Rayleigh-quotient iteration started from an arbitrary vector
/// cannot guarantee.  [`compute_mode_shapes`](Self::compute_mode_shapes)
/// recovers the shape of the mode nearest a prescribed frequency by shifted
/// inverse iteration.
pub struct ModalAnalysis {
    /// Iteration cap: Jacobi sweeps for the frequencies / modes and
    /// inverse-iteration steps for a single mode shape.
    pub max_iter: usize,
    /// Convergence tolerance on the change of the normalised mode shape in the
    /// shifted inverse iteration (used by [`mode_shape`](Self::mode_shape)).
    pub tol: f64,
}

impl Default for ModalAnalysis {
    fn default() -> Self {
        Self {
            max_iter: 500,
            tol: MODE_SHAPE_TOLERANCE,
        }
    }
}

impl ModalAnalysis {
    /// Create with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Dot product helper.
    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    /// Dense matrix-vector product (n×n matrix stored row-major).
    fn matvec_dense(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
        let mut y = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                y[i] += a[i * n + j] * x[j];
            }
        }
        y
    }

    /// The `n_modes` lowest natural frequencies (rad/s) of `K·φ = ω²·M·φ` in
    /// ascending order, using this instance's [`max_iter`](Self::max_iter).
    ///
    /// See [`compute_natural_frequencies`](Self::compute_natural_frequencies).
    pub fn frequencies(
        &self,
        k: &[f64],
        m: &[f64],
        n: usize,
        n_modes: usize,
    ) -> Result<Vec<f64>, FemError> {
        Self::compute_natural_frequencies(k, m, n, n_modes, self.max_iter)
    }

    /// Mode shape nearest `omega`, using this instance's
    /// [`max_iter`](Self::max_iter) and [`tol`](Self::tol).
    ///
    /// See [`compute_mode_shapes`](Self::compute_mode_shapes).
    pub fn mode_shape(
        &self,
        k: &[f64],
        m: &[f64],
        n: usize,
        omega: f64,
    ) -> Result<Vec<f64>, FemError> {
        Self::shifted_inverse_iteration(k, m, n, omega, self.max_iter, self.tol)
    }

    /// The `n_modes` lowest natural frequencies (rad/s) of the dense
    /// generalized eigenproblem `K·φ = ω²·M·φ`, in **ascending** order.
    ///
    /// `k` and `m` are row-major `n × n` arrays; `max_iter` caps the Jacobi
    /// sweeps (about ten are needed in practice).  At most `n` frequencies are
    /// returned.  Eigenvalues that are negative only by round-off (rigid-body
    /// modes of a free structure) are reported as `ω = 0`.
    ///
    /// # Errors
    ///
    /// - [`FemError::InvalidInput`] — wrongly sized or non-finite `k` / `m`, a
    ///   non-symmetric `k` or `m`, a mass matrix that is not positive definite
    ///   (e.g. a massless DOF), or a stiffness matrix with a genuinely negative
    ///   eigenvalue (an unstable structure has no real frequency);
    /// - [`FemError::SolverNotConverged`] — the Jacobi iteration needs more than
    ///   `max_iter` sweeps.
    pub fn compute_natural_frequencies(
        k: &[f64],
        m: &[f64],
        n: usize,
        n_modes: usize,
        max_iter: usize,
    ) -> Result<Vec<f64>, FemError> {
        Ok(Self::compute_modes(k, m, n, n_modes, max_iter)?
            .into_iter()
            .map(|(omega, _)| omega)
            .collect())
    }

    /// The `n_modes` lowest modes `(ω, φ)` of `K·φ = ω²·M·φ` in ascending order
    /// of the angular frequency `ω` (rad/s).
    ///
    /// The shapes are mass-orthonormal (`φᵢᵀ·M·φⱼ = δᵢⱼ`) and carry a
    /// deterministic sign: the first component whose magnitude is (within
    /// `1e-9` relative) the largest one is positive.
    /// Arguments and errors are those of
    /// [`compute_natural_frequencies`](Self::compute_natural_frequencies).
    pub fn compute_modes(
        k: &[f64],
        m: &[f64],
        n: usize,
        n_modes: usize,
        max_iter: usize,
    ) -> Result<Vec<(f64, Vec<f64>)>, FemError> {
        let eigen = dense::generalized_symmetric_eigen(
            ("stiffness matrix K", k),
            ("mass matrix M", m),
            n,
            max_iter,
        )?;
        let count = n_modes.min(n);
        let spectral_radius = eigen.values.iter().fold(0.0_f64, |r, v| r.max(v.abs()));
        let mut modes = Vec::with_capacity(count);
        for (index, &lambda) in eigen.values.iter().enumerate().take(count) {
            if lambda < -NEGATIVE_EIGENVALUE_TOLERANCE * spectral_radius {
                return Err(FemError::InvalidInput(format!(
                    "stiffness matrix K has a negative eigenvalue {lambda:e} (relative to the \
                     spectral radius {spectral_radius:e}): the structure is unstable and has no \
                     real natural frequency"
                )));
            }
            let mut shape = eigen.eigenvector(index);
            dense::normalise_sign(&mut shape);
            modes.push((lambda.max(0.0).sqrt(), shape));
        }
        Ok(modes)
    }

    /// Mode shape of the mode whose eigenvalue is nearest `omega²`, by shifted
    /// inverse iteration: repeatedly solve `(K − ω²·M)·y = M·φ`.
    ///
    /// `k` and `m` are row-major `n × n` arrays and `max_iter` caps the
    /// iteration steps.  The returned shape is mass-normalised
    /// (`φᵀ·M·φ = 1`) with a deterministic sign.
    ///
    /// When `omega` is itself a natural frequency, `K − ω²·M` is singular; the
    /// shift is then perturbed by `1e-10·max(ω², 1)` (grown a hundredfold per
    /// retry, three retries) so that the iteration converges to that mode in a
    /// step or two instead of failing.
    ///
    /// # Errors
    ///
    /// - [`FemError::InvalidInput`] — `max_iter == 0`, a non-finite `omega`,
    ///   wrongly sized or non-finite `k` / `m`, or a `M` that is not positive
    ///   along the computed shape;
    /// - [`FemError::SingularSystem`] — `K − ω²·M` stays singular under every
    ///   perturbation (a singular matrix pencil);
    /// - [`FemError::SolverNotConverged`] — the iteration did not settle within
    ///   `max_iter` steps (for instance `ω²` midway between two eigenvalues).
    pub fn compute_mode_shapes(
        k: &[f64],
        m: &[f64],
        n: usize,
        omega: f64,
        max_iter: usize,
    ) -> Result<Vec<f64>, FemError> {
        Self::shifted_inverse_iteration(k, m, n, omega, max_iter, MODE_SHAPE_TOLERANCE)
    }

    /// Shared implementation of [`compute_mode_shapes`](Self::compute_mode_shapes)
    /// and [`mode_shape`](Self::mode_shape).
    fn shifted_inverse_iteration(
        k: &[f64],
        m: &[f64],
        n: usize,
        omega: f64,
        max_iter: usize,
        tol: f64,
    ) -> Result<Vec<f64>, FemError> {
        if max_iter == 0 {
            return Err(FemError::InvalidInput(
                "max_iter must be at least 1".to_string(),
            ));
        }
        if !omega.is_finite() {
            return Err(FemError::InvalidInput(format!(
                "angular frequency must be finite, got {omega}"
            )));
        }
        if !(tol.is_finite() && tol > 0.0) {
            return Err(FemError::InvalidInput(format!(
                "convergence tolerance must be finite and positive, got {tol}"
            )));
        }
        dense::validate_square("stiffness matrix K", k, n)?;
        dense::validate_square("mass matrix M", m, n)?;
        if n == 0 {
            return Ok(Vec::new());
        }

        let base_shift = omega * omega;
        let mut factors = None;
        for attempt in 0..=SHIFT_RETRIES {
            let shift = if attempt == 0 {
                base_shift
            } else {
                let growth = 100.0_f64.powi(attempt as i32 - 1);
                base_shift + SHIFT_PERTURBATION * base_shift.max(1.0) * growth
            };
            let shifted: Vec<f64> = k
                .iter()
                .zip(m.iter())
                .map(|(k_ij, m_ij)| k_ij - shift * m_ij)
                .collect();
            match dense::LuFactors::factor(&shifted, n) {
                Ok(lu) => {
                    factors = Some(lu);
                    break;
                }
                Err(Error::SingularSystem) => continue,
                Err(other) => return Err(other.into()),
            }
        }
        let factors = factors.ok_or(FemError::SingularSystem)?;

        // Deterministic start vector with no special symmetry (additive
        // recurrence of the golden ratio), so that it is not accidentally
        // M-orthogonal to a mode of a symmetric structure.
        let mut q: Vec<f64> = (0..n)
            .map(|i| ((i + 1) as f64 * 0.618_033_988_749_894_9).fract() - 0.5)
            .collect();
        Self::scale_to_unit_norm(&mut q)?;

        let mut converged = false;
        for _ in 0..max_iter {
            let rhs = Self::matvec_dense(m, &q, n);
            let mut y = factors.solve(&rhs)?;
            Self::scale_to_unit_norm(&mut y)?;
            // The direction matters, not the sign of the iterate.
            let same: f64 = q.iter().zip(&y).map(|(a, b)| (a - b) * (a - b)).sum();
            let flipped: f64 = q.iter().zip(&y).map(|(a, b)| (a + b) * (a + b)).sum();
            q = y;
            if same.min(flipped).sqrt() < tol {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(FemError::SolverNotConverged);
        }

        let mass_q = Self::matvec_dense(m, &q, n);
        let modal_mass = Self::dot(&q, &mass_q);
        if !(modal_mass.is_finite() && modal_mass > 0.0) {
            return Err(FemError::InvalidInput(format!(
                "mass matrix M is not positive along the computed mode shape (φᵀMφ = {modal_mass:e})"
            )));
        }
        let scale = modal_mass.sqrt();
        for value in q.iter_mut() {
            *value /= scale;
        }
        dense::normalise_sign(&mut q);
        Ok(q)
    }

    /// Scale `v` to unit Euclidean norm; a zero or non-finite vector means the
    /// shifted system overflowed or collapsed.
    fn scale_to_unit_norm(v: &mut [f64]) -> Result<(), FemError> {
        let norm = Self::dot(v, v).sqrt();
        if !(norm.is_finite() && norm > 0.0) {
            return Err(FemError::SingularSystem);
        }
        for value in v.iter_mut() {
            *value /= norm;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// DynamicAnalysis — Newmark-β and Wilson-θ
// ---------------------------------------------------------------------------

/// State `(displacement, velocity, acceleration)` produced by one step of
/// [`DynamicAnalysis::newmark_beta_step`] or [`DynamicAnalysis::wilson_theta_step`].
pub type StepState = (Vec<f64>, Vec<f64>, Vec<f64>);

/// The two direct time-integration schemes; each has its own constant step
/// matrix and therefore its own slot in the factorisation cache.
#[derive(Debug, Clone, Copy)]
enum StepScheme {
    /// Newmark-β: `M + γ·dt·C + β·dt²·K`.
    Newmark,
    /// Wilson-θ: `K + 6/τ²·M + 3/τ·C`.
    Wilson,
}

impl StepScheme {
    /// Slot of this scheme in the factorisation cache.
    fn slot(self) -> usize {
        match self {
            StepScheme::Newmark => 0,
            StepScheme::Wilson => 1,
        }
    }
}

/// An LU factorisation together with the exact matrix it was computed from.
struct CachedFactorization {
    /// The factorised matrix.  A later step reuses the factors only when it
    /// assembles exactly this matrix, so changing `dt`, `β`, `γ`, `θ` — or
    /// mutating the public `m`, `c`, `k` — can never return a stale solution.
    matrix: Vec<f64>,
    /// Shared factors, so the substitution can run outside the cache lock.
    factors: Arc<dense::LuFactors>,
}

/// Direct time-integration for structural dynamics: M·ü + C·u̇ + K·u = f(t).
///
/// Provides Newmark-β and Wilson-θ stepping methods.  The step matrix of each
/// scheme is constant while the time step and the scheme parameters stay
/// unchanged, so it is LU-factorised once and the factors are reused by every
/// following step (`O(n²)` per step instead of `O(n³)`); see
/// [`factorization_count`](Self::factorization_count).
pub struct DynamicAnalysis {
    /// Mass matrix (row-major, n×n).
    pub m: Vec<f64>,
    /// Damping matrix (row-major, n×n).
    pub c: Vec<f64>,
    /// Stiffness matrix (row-major, n×n).
    pub k: Vec<f64>,
    /// Number of DOFs.
    pub n: usize,
    /// Latest factorisation per scheme, see [`CachedFactorization`].
    step_factors: Mutex<[Option<CachedFactorization>; 2]>,
    /// Number of LU factorisations performed so far.
    factorizations: AtomicUsize,
}

impl DynamicAnalysis {
    /// Construct a new `DynamicAnalysis`.
    pub fn new(m: Vec<f64>, c: Vec<f64>, k: Vec<f64>, n: usize) -> Self {
        Self {
            m,
            c,
            k,
            n,
            step_factors: Mutex::new([None, None]),
            factorizations: AtomicUsize::new(0),
        }
    }

    /// Number of LU factorisations of a step matrix performed so far.
    ///
    /// Stepping repeatedly with the same `dt` and scheme parameters leaves this
    /// count unchanged after the first step; it grows only when a step matrix
    /// actually changes.
    pub fn factorization_count(&self) -> usize {
        self.factorizations.load(Ordering::Relaxed)
    }

    /// Dense matrix-vector product.
    fn mv(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
        let mut y = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                y[i] += a[i * n + j] * x[j];
            }
        }
        y
    }

    /// Check the system matrices against `n`.
    fn validate_system(&self) -> Result<(), FemError> {
        dense::validate_square("mass matrix M", &self.m, self.n)?;
        dense::validate_square("damping matrix C", &self.c, self.n)?;
        dense::validate_square("stiffness matrix K", &self.k, self.n)?;
        Ok(())
    }

    /// Check that a state / load vector has `n` finite entries.
    fn validate_vector(&self, name: &str, values: &[f64]) -> Result<(), FemError> {
        if values.len() != self.n {
            return Err(FemError::InvalidInput(format!(
                "{name} has {} entries but the system has {} DOFs",
                values.len(),
                self.n
            )));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(FemError::InvalidInput(format!(
                "{name} contains a non-finite entry"
            )));
        }
        Ok(())
    }

    /// Factors of `matrix` for `scheme`, reusing the cached factorisation when
    /// `matrix` is exactly the one it was computed from.
    fn step_factors_for(
        &self,
        scheme: StepScheme,
        matrix: Vec<f64>,
    ) -> Result<Arc<dense::LuFactors>, FemError> {
        // The cache holds whole `Option` values that are replaced atomically, so
        // the data is consistent even if another thread panicked while holding
        // the lock; a poisoned lock is therefore safe to recover.
        let mut cache = match self.step_factors.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        };
        let slot = scheme.slot();
        if let Some(entry) = &cache[slot]
            && entry.matrix == matrix
        {
            return Ok(Arc::clone(&entry.factors));
        }
        let factors = Arc::new(dense::LuFactors::factor(&matrix, self.n)?);
        self.factorizations.fetch_add(1, Ordering::Relaxed);
        cache[slot] = Some(CachedFactorization {
            matrix,
            factors: Arc::clone(&factors),
        });
        Ok(factors)
    }

    /// Reject a step result that contains a non-finite number.
    fn finite_state(state: StepState) -> Result<StepState, FemError> {
        let all_finite = state
            .0
            .iter()
            .chain(state.1.iter())
            .chain(state.2.iter())
            .all(|value| value.is_finite());
        if all_finite {
            Ok(state)
        } else {
            Err(FemError::SingularSystem)
        }
    }

    /// One step of the Newmark-β method.
    ///
    /// Updates displacement `u`, velocity `v`, and acceleration `a` from time `t`
    /// to time `t + dt`.
    ///
    /// The Newmark updates are `u₊ = ũ + β·dt²·a₊` and `v₊ = ṽ + γ·dt·a₊` with
    /// the predictors `ũ = u + dt·v + dt²·(½ − β)·a` and `ṽ = v + dt·(1 − γ)·a`.
    /// Substituting them into the equation of motion at `t + dt`,
    /// `M·a₊ + C·v₊ + K·u₊ = f₊`, gives a system for the *acceleration*
    ///
    /// `(M + γ·dt·C + β·dt²·K) · a₊ = f₊ − K·ũ − C·ṽ`.
    ///
    /// The matrix on the left is constant while `dt`, `β`, `γ` and the system
    /// matrices are unchanged; its LU factors are cached and reused.
    ///
    /// Typical parameter choices:
    /// - Constant average acceleration (unconditionally stable): β = 0.25, γ = 0.5
    /// - Linear acceleration (conditionally stable): β = 1/6, γ = 0.5
    /// - Central difference (explicit, conditionally stable): β = 0, γ = 0.5
    ///
    /// # Returns
    /// `(u_{n+1}, v_{n+1}, a_{n+1})`
    ///
    /// # Errors
    ///
    /// - [`FemError::InvalidInput`] — a system matrix or a state / load vector
    ///   with the wrong size or a non-finite entry, `dt` not finite and
    ///   positive, `β` negative or non-finite, or a non-finite `γ`;
    /// - [`FemError::SingularSystem`] — the step matrix
    ///   `M + γ·dt·C + β·dt²·K` is singular to working precision, or the step
    ///   overflowed.
    pub fn newmark_beta_step(
        &self,
        u: &[f64],
        v: &[f64],
        a: &[f64],
        f_ext: &[f64],
        dt: f64,
        beta: f64,
        gamma: f64,
    ) -> Result<StepState, FemError> {
        let n = self.n;
        self.validate_system()?;
        self.validate_vector("displacement u", u)?;
        self.validate_vector("velocity v", v)?;
        self.validate_vector("acceleration a", a)?;
        self.validate_vector("external load f_ext", f_ext)?;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(FemError::InvalidInput(format!(
                "time step must be finite and positive, got {dt}"
            )));
        }
        if !(beta.is_finite() && beta >= 0.0 && gamma.is_finite()) {
            return Err(FemError::InvalidInput(format!(
                "Newmark parameters must be finite with β >= 0, got β = {beta}, γ = {gamma}"
            )));
        }

        // Predictors
        let u_pred: Vec<f64> = (0..n)
            .map(|i| u[i] + dt * v[i] + dt * dt * (0.5 - beta) * a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n).map(|i| v[i] + dt * (1.0 - gamma) * a[i]).collect();

        // System matrix for the acceleration: A = M + γ·dt·C + β·dt²·K
        let c_damping = gamma * dt;
        let c_stiffness = beta * dt * dt;
        let a_eff: Vec<f64> = (0..n * n)
            .map(|idx| self.m[idx] + c_damping * self.c[idx] + c_stiffness * self.k[idx])
            .collect();
        let factors = self.step_factors_for(StepScheme::Newmark, a_eff)?;

        // Right-hand side: f₊ − K·ũ − C·ṽ
        let ku_pred = Self::mv(&self.k, &u_pred, n);
        let cv_pred = Self::mv(&self.c, &v_pred, n);
        let rhs: Vec<f64> = (0..n).map(|i| f_ext[i] - ku_pred[i] - cv_pred[i]).collect();

        let a_new = factors.solve(&rhs)?;
        let u_new: Vec<f64> = (0..n)
            .map(|i| u_pred[i] + beta * dt * dt * a_new[i])
            .collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + gamma * dt * a_new[i]).collect();

        Self::finite_state((u_new, v_new, a_new))
    }

    /// One step of the Wilson-θ method (Bathe, *Finite Element Procedures*, ch. 9).
    ///
    /// The acceleration is assumed to vary linearly over the extended interval
    /// `[t, t + τ]`, `τ = θ·dt`.  The displacement at `t + τ` follows from
    ///
    /// `(K + 6/τ²·M + 3/τ·C) · u_τ = f̂`
    ///
    /// with `f̂ = f_prev + θ·(f_ext − f_prev) + M·(6/τ²·u + 6/τ·v + 2a)
    /// + C·(3/τ·u + 2v + τ/2·a)`, and the state at `t + dt` is recovered from
    /// `a₊ = 6/(θ·τ²)·(u_τ − u) − 6/(θ·τ)·v + (1 − 3/θ)·a`,
    /// `v₊ = v + dt/2·(a + a₊)`, `u₊ = u + dt·v + dt²/6·(2a + a₊)`.
    ///
    /// The matrix on the left is constant while `dt`, `θ` and the system
    /// matrices are unchanged; its LU factors are cached and reused.
    ///
    /// `theta` must be at least 1.366 for unconditional stability (1.37–1.42 is
    /// typical); `theta = 1` reduces to the linear-acceleration method.
    ///
    /// # Arguments
    /// * `u`, `v`, `a` – displacement, velocity and acceleration at time `t`
    /// * `f_prev`      – external load at time `t`
    /// * `f_ext`       – external load at time `t + dt`
    ///
    /// # Returns
    /// `(u_{n+1}, v_{n+1}, a_{n+1})`
    ///
    /// # Errors
    ///
    /// - [`FemError::InvalidInput`] — a system matrix or a state / load vector
    ///   with the wrong size or a non-finite entry, `dt` not finite and
    ///   positive, or `theta` non-finite or below 1;
    /// - [`FemError::SingularSystem`] — the effective stiffness
    ///   `K + 6/τ²·M + 3/τ·C` is singular to working precision, or the step
    ///   overflowed.
    pub fn wilson_theta_step(
        &self,
        u: &[f64],
        v: &[f64],
        a: &[f64],
        f_prev: &[f64],
        f_ext: &[f64],
        dt: f64,
        theta: f64,
    ) -> Result<StepState, FemError> {
        let n = self.n;
        self.validate_system()?;
        self.validate_vector("displacement u", u)?;
        self.validate_vector("velocity v", v)?;
        self.validate_vector("acceleration a", a)?;
        self.validate_vector("previous external load f_prev", f_prev)?;
        self.validate_vector("external load f_ext", f_ext)?;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(FemError::InvalidInput(format!(
                "time step must be finite and positive, got {dt}"
            )));
        }
        if !(theta.is_finite() && theta >= 1.0) {
            return Err(FemError::InvalidInput(format!(
                "Wilson parameter θ must be finite and at least 1, got {theta}"
            )));
        }
        let tau = theta * dt;

        // Effective stiffness K̂ = K + 6/τ²·M + 3/τ·C
        let c_mass = 6.0 / (tau * tau);
        let c_damping = 3.0 / tau;
        let k_hat: Vec<f64> = (0..n * n)
            .map(|idx| self.k[idx] + c_mass * self.m[idx] + c_damping * self.c[idx])
            .collect();
        let factors = self.step_factors_for(StepScheme::Wilson, k_hat)?;

        // Effective load at t + τ:
        //   f̂ = f_prev + θ(f_ext − f_prev) + M(6/τ²·u + 6/τ·v + 2a) + C(3/τ·u + 2v + τ/2·a)
        let mass_arg: Vec<f64> = (0..n)
            .map(|i| c_mass * u[i] + 6.0 / tau * v[i] + 2.0 * a[i])
            .collect();
        let damping_arg: Vec<f64> = (0..n)
            .map(|i| c_damping * u[i] + 2.0 * v[i] + 0.5 * tau * a[i])
            .collect();
        let mass_term = Self::mv(&self.m, &mass_arg, n);
        let damping_term = Self::mv(&self.c, &damping_arg, n);
        let rhs: Vec<f64> = (0..n)
            .map(|i| f_prev[i] + theta * (f_ext[i] - f_prev[i]) + mass_term[i] + damping_term[i])
            .collect();

        // Displacement at t + τ
        let u_tau = factors.solve(&rhs)?;

        // Acceleration, velocity and displacement at t + dt
        let a_new: Vec<f64> = (0..n)
            .map(|i| {
                6.0 / (theta * tau * tau) * (u_tau[i] - u[i]) - 6.0 / (theta * tau) * v[i]
                    + (1.0 - 3.0 / theta) * a[i]
            })
            .collect();
        let v_new: Vec<f64> = (0..n)
            .map(|i| v[i] + dt / 2.0 * (a[i] + a_new[i]))
            .collect();
        let u_new: Vec<f64> = (0..n)
            .map(|i| u[i] + dt * v[i] + dt * dt / 6.0 * (2.0 * a[i] + a_new[i]))
            .collect();

        Self::finite_state((u_new, v_new, a_new))
    }
}

// ---------------------------------------------------------------------------
// Stress state and invariants
// ---------------------------------------------------------------------------

/// Full 3-D stress state in Voigt notation: \[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz\].
#[derive(Debug, Clone, Copy)]
pub struct StressState {
    /// Voigt stress components.
    pub sigma: [f64; 6],
}

impl StressState {
    /// Construct from a Voigt array.
    pub fn from_voigt(sigma: [f64; 6]) -> Self {
        Self { sigma }
    }

    /// Construct from individual components.
    pub fn new(sxx: f64, syy: f64, szz: f64, txy: f64, tyz: f64, txz: f64) -> Self {
        Self {
            sigma: [sxx, syy, szz, txy, tyz, txz],
        }
    }

    /// Hydrostatic (mean) stress: p = (σ_xx + σ_yy + σ_zz) / 3.
    pub fn hydrostatic(&self) -> f64 {
        (self.sigma[0] + self.sigma[1] + self.sigma[2]) / 3.0
    }

    /// First stress invariant I₁ = σ_xx + σ_yy + σ_zz.
    pub fn i1(&self) -> f64 {
        self.sigma[0] + self.sigma[1] + self.sigma[2]
    }

    /// Second deviatoric stress invariant J₂.
    ///
    /// J₂ = ½ s_ij s_ij  where  s_ij = σ_ij − p δ_ij.
    pub fn j2(&self) -> f64 {
        let p = self.hydrostatic();
        let sx = self.sigma[0] - p;
        let sy = self.sigma[1] - p;
        let sz = self.sigma[2] - p;
        let txy = self.sigma[3];
        let tyz = self.sigma[4];
        let txz = self.sigma[5];
        0.5 * (sx * sx + sy * sy + sz * sz + 2.0 * (txy * txy + tyz * tyz + txz * txz))
    }

    /// Third deviatoric stress invariant J₃.
    pub fn j3(&self) -> f64 {
        let p = self.hydrostatic();
        let sx = self.sigma[0] - p;
        let sy = self.sigma[1] - p;
        let sz = self.sigma[2] - p;
        let txy = self.sigma[3];
        let tyz = self.sigma[4];
        let txz = self.sigma[5];
        sx * (sy * sz - tyz * tyz) - txy * (txy * sz - tyz * txz) + txz * (txy * tyz - sy * txz)
    }

    /// Von Mises equivalent stress: σ_vm = √(3 J₂).
    pub fn von_mises(&self) -> f64 {
        (3.0 * self.j2()).max(0.0).sqrt()
    }

    /// Tresca equivalent stress: max principal − min principal.
    pub fn tresca_stress(&self) -> f64 {
        let [s1, _s2, s3] = self.principal_stresses();
        s1 - s3
    }

    /// Principal stresses sorted in descending order \[σ₁ ≥ σ₂ ≥ σ₃\].
    ///
    /// Uses analytical solution for the characteristic equation of the
    /// symmetric 3×3 stress tensor.
    pub fn principal_stresses(&self) -> [f64; 3] {
        let p = self.hydrostatic();
        let j2 = self.j2();
        let j3 = self.j3();
        let r = (j2 / 3.0).max(0.0).sqrt();
        if r < 1e-30 {
            return [p, p, p];
        }
        let cos_arg = (j3 / (2.0 * r * r * r)).clamp(-1.0, 1.0);
        let theta = cos_arg.acos() / 3.0;
        use std::f64::consts::PI;
        let s1 = p + 2.0 * r * theta.cos();
        let s2 = p + 2.0 * r * (theta - 2.0 * PI / 3.0).cos();
        let s3 = p + 2.0 * r * (theta - 4.0 * PI / 3.0).cos();
        let mut arr = [s1, s2, s3];
        arr.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        arr
    }

    /// Lode angle in degrees θ ∈ \[0°, 60°\].
    ///
    /// θ = 0°: axisymmetric tension; θ = 60°: axisymmetric compression.
    pub fn lode_angle_deg(&self) -> f64 {
        let j2 = self.j2();
        let j3 = self.j3();
        if j2 < 1e-30 {
            return 0.0;
        }
        let cos_arg = (3.0_f64.sqrt() * 3.0 * j3 / (2.0 * j2.powf(1.5))).clamp(-1.0, 1.0);
        let theta_rad = cos_arg.acos() / 3.0;
        theta_rad.to_degrees()
    }

    /// Safety factor by the von Mises criterion.
    ///
    /// `sf = yield_stress / sigma_vm`.  Returns `f64::INFINITY` if von Mises stress is zero.
    pub fn safety_factor_von_mises(&self, yield_stress: f64) -> f64 {
        let vm = self.von_mises();
        if vm < 1e-30 {
            return f64::INFINITY;
        }
        yield_stress / vm
    }

    /// Safety factor by the Tresca criterion.
    pub fn safety_factor_tresca(&self, yield_stress: f64) -> f64 {
        let tresca = self.tresca_stress();
        if tresca < 1e-30 {
            return f64::INFINITY;
        }
        yield_stress / tresca
    }

    /// Drucker-Prager equivalent stress: σ_dp = α·I₁ + √J₂.
    ///
    /// Useful for granular materials where hydrostatic pressure affects yielding.
    pub fn drucker_prager(&self, alpha: f64) -> f64 {
        alpha * self.i1() + self.j2().max(0.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Fatigue life estimator (Basquin's law)
// ---------------------------------------------------------------------------

/// Estimates high-cycle fatigue life using Basquin's power law:
///
/// N = N_ref * (σ_ref / σ_a)^m
///
/// where σ_a is stress amplitude, σ_ref is reference stress at N_ref cycles,
/// and m = log(N_ref) / log(σ_uts / σ_ref) (slope exponent).
#[derive(Debug, Clone)]
pub struct FatigueLifeEstimator {
    /// Ultimate tensile strength (UTS) \[Pa\].
    pub sigma_uts: f64,
    /// Reference number of cycles (e.g., 10⁶).
    pub n_ref: f64,
    /// Endurance limit \[Pa\]. Stress amplitudes below this never fail.
    pub endurance_limit: f64,
}

impl FatigueLifeEstimator {
    /// Construct a new Basquin fatigue model.
    pub fn new(sigma_uts: f64, n_ref: f64, endurance_limit: f64) -> Self {
        assert!(sigma_uts > 0.0);
        assert!(n_ref > 0.0);
        assert!(endurance_limit >= 0.0 && endurance_limit < sigma_uts);
        Self {
            sigma_uts,
            n_ref,
            endurance_limit,
        }
    }

    /// Compute the number of cycles to failure for a given stress amplitude.
    ///
    /// Returns `f64::INFINITY` if the stress amplitude is at or below the endurance limit.
    pub fn cycles_to_failure(&self, sigma_a: f64) -> f64 {
        if sigma_a <= self.endurance_limit {
            return f64::INFINITY;
        }
        if sigma_a >= self.sigma_uts {
            return 1.0;
        }
        // Basquin exponent: assume reference stress = σ_uts / 2 at N_ref
        // N = N_ref * (sigma_ref / sigma_a)^m, where sigma_ref = sigma_uts at N=1
        let m = self.basquin_exponent();
        self.n_ref * (self.sigma_uts / sigma_a).powf(m)
    }

    /// Basquin exponent m from N_ref and σ_uts.
    ///
    /// Derived from the two-point definition: (σ_uts, 1) and (σ_ref, N_ref).
    pub fn basquin_exponent(&self) -> f64 {
        // N_ref = (σ_uts / σ_ref)^m
        // Assuming σ_ref = 0.5 * σ_uts:
        self.n_ref.log10() / 2.0_f64.log10()
    }

    /// Goodman correction for mean stress σ_mean.
    ///
    /// Reduces the allowable stress amplitude: σ_a_eff = σ_a / (1 - σ_mean/σ_uts).
    pub fn goodman_corrected_amplitude(&self, sigma_a: f64, sigma_mean: f64) -> f64 {
        let denom = 1.0 - sigma_mean / self.sigma_uts;
        if denom <= 0.0 {
            return f64::INFINITY;
        }
        sigma_a / denom
    }
}

/// Compute Palmgren-Miner cumulative damage from a slice of cycle ratios.
///
/// D = Σ (n_i / N_i)  where each element is the ratio n_i/N_i.
/// Failure is predicted when D ≥ 1.
pub fn miner_damage(cycle_ratios: &[f64]) -> f64 {
    cycle_ratios.iter().sum()
}

// ---------------------------------------------------------------------------
// Post-processing utilities
// ---------------------------------------------------------------------------

/// Compute element-averaged strains from a displacement field and B-matrices.
///
/// For each element, ε = B · u_local.  Stores result as `[f64; 6]` in Voigt order.
pub fn compute_element_strains(
    displacements: &[f64],
    element_dofs: &[Vec<usize>],
    b_matrices: &[Vec<f64>],
) -> Vec<[f64; 6]> {
    element_dofs
        .iter()
        .zip(b_matrices.iter())
        .map(|(dofs, b)| {
            let mut eps = [0.0f64; 6];
            let n_dof = dofs.len();
            for row in 0..6 {
                let mut sum = 0.0;
                for col in 0..n_dof {
                    if dofs[col] < displacements.len() {
                        sum += b[row * n_dof + col] * displacements[dofs[col]];
                    }
                }
                eps[row] = sum;
            }
            eps
        })
        .collect()
}

/// Compute element von Mises stresses from a stress field.
pub fn von_mises_field(stresses: &[[f64; 6]]) -> Vec<f64> {
    stresses
        .iter()
        .map(|&s| StressState::from_voigt(s).von_mises())
        .collect()
}

/// Compute element safety factors by von Mises criterion.
pub fn safety_factor_field(stresses: &[[f64; 6]], yield_stress: f64) -> Vec<f64> {
    stresses
        .iter()
        .map(|&s| StressState::from_voigt(s).safety_factor_von_mises(yield_stress))
        .collect()
}

/// Smooth a field over neighboring elements by averaging with weights.
///
/// Useful for stress smoothing after solving (superconvergent patch recovery
/// approximation via simple neighbor average).
pub fn smooth_field(values: &[f64], neighbors: &[Vec<usize>]) -> Vec<f64> {
    values
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            let nbrs = &neighbors[i];
            if nbrs.is_empty() {
                return v;
            }
            let sum: f64 = nbrs
                .iter()
                .map(|&j| if j < values.len() { values[j] } else { 0.0 })
                .sum();
            (v + sum) / (1 + nbrs.len()) as f64
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Rayleigh damping
// ---------------------------------------------------------------------------

/// Rayleigh damping: C = α·M + β·K.
///
/// α and β are derived from two modal damping ratios at given frequencies.
#[derive(Debug, Clone, Copy)]
pub struct RayleighDamping {
    /// Mass-proportional coefficient.
    pub alpha: f64,
    /// Stiffness-proportional coefficient.
    pub beta: f64,
}

impl RayleighDamping {
    /// Compute Rayleigh coefficients from two target modes.
    ///
    /// At ω_i and ω_j with damping ratios ζ_i and ζ_j:
    /// \[α; β\] = 2·\[ω_i, ω_j; 1/ω_i, 1/ω_j\]^{-1} · \[ζ_i; ζ_j\]
    pub fn from_two_modes(omega_i: f64, omega_j: f64, zeta_i: f64, zeta_j: f64) -> Self {
        // Solve 2×2 system:
        // ζ_i = α/(2ω_i) + β·ω_i/2
        // ζ_j = α/(2ω_j) + β·ω_j/2
        let det = (omega_j * omega_j - omega_i * omega_i) / (2.0 * omega_i * omega_j);
        if det.abs() < 1e-30 {
            return Self {
                alpha: 0.0,
                beta: 0.0,
            };
        }
        // Rearranged from 2x2 linear system
        let alpha = 2.0 * omega_i * omega_j * (zeta_i * omega_j - zeta_j * omega_i)
            / (omega_j * omega_j - omega_i * omega_i);
        let beta =
            2.0 * (zeta_j * omega_j - zeta_i * omega_i) / (omega_j * omega_j - omega_i * omega_i);
        Self { alpha, beta }
    }

    /// Damping ratio at a given angular frequency ω.
    pub fn damping_ratio(&self, omega: f64) -> f64 {
        if omega < 1e-30 {
            return 0.0;
        }
        self.alpha / (2.0 * omega) + self.beta * omega / 2.0
    }

    /// Assemble dense damping matrix C = α·M + β·K.
    pub fn assemble(&self, mass: &[f64], stiffness: &[f64], n: usize) -> Vec<f64> {
        assert_eq!(mass.len(), n * n);
        assert_eq!(stiffness.len(), n * n);
        (0..n * n)
            .map(|k| self.alpha * mass[k] + self.beta * stiffness[k])
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- StaticAnalysis tests ----

    #[test]
    fn test_static_solve_2dof() {
        // K = [[4, -1],[-1, 3]], f = [1, 2]
        // Solution: 4u0 - u1 = 1, -u0 + 3u1 = 2 → u0 = 5/11, u1 = 9/11
        let mut sa = StaticAnalysis::new(2);
        sa.stiffness_matrix = vec![4.0, -1.0, -1.0, 3.0];
        sa.force_vector = vec![1.0, 2.0];
        let u = sa.solve().unwrap();
        let expected0 = 5.0 / 11.0;
        let expected1 = 9.0 / 11.0;
        assert!(
            (u[0] - expected0).abs() < 1e-12,
            "u[0] = {}, expected {expected0}",
            u[0]
        );
        assert!(
            (u[1] - expected1).abs() < 1e-12,
            "u[1] = {}, expected {expected1}",
            u[1]
        );
    }

    #[test]
    fn test_static_solve_singular() {
        let mut sa = StaticAnalysis::new(2);
        sa.stiffness_matrix = vec![1.0, 1.0, 1.0, 1.0];
        sa.force_vector = vec![1.0, 1.0];
        let result = sa.solve();
        assert!(result.is_err(), "expected SingularSystem error");
    }

    #[test]
    fn test_static_assemble_and_solve() {
        // Single 2-DOF element: ke = [[2,-2],[-2,2]], dofs=[0,1]
        // Fix DOF 0: set K[0,0]=1, K[0,1]=0, f[0]=0
        let mut sa = StaticAnalysis::new(2);
        let elements = vec![(vec![0usize, 1], vec![2.0, -2.0, -2.0, 2.0])];
        sa.assemble_global_stiffness(&elements);
        // Apply Dirichlet BC at DOF 0
        sa.stiffness_matrix[0] = 1.0;
        sa.stiffness_matrix[1] = 0.0;
        sa.stiffness_matrix[2] = 0.0;
        sa.force_vector[0] = 0.0;
        sa.force_vector[1] = 1.0; // unit force on DOF 1
        let u = sa.solve().unwrap();
        assert!((u[0]).abs() < 1e-12, "u[0] should be 0, got {}", u[0]);
        assert!(u[1] > 0.0, "u[1] should be positive, got {}", u[1]);
    }

    #[test]
    fn test_compute_strains() {
        // 1-DOF trivial: B = [1.0, 0, 0, 0, 0, 0], u=[0.5]
        let sa = StaticAnalysis::new(1);
        let b_mat: Vec<f64> = {
            let mut v = vec![0.0; 6];
            v[0] = 1.0;
            v
        };
        let elements = vec![(vec![0usize], b_mat)];
        let u = vec![0.5_f64];
        let strains = sa.compute_strains(&u, &elements);
        assert_eq!(strains.len(), 1);
        assert!((strains[0][0] - 0.5).abs() < 1e-15);
        for (k, &val) in strains[0].iter().enumerate().skip(1) {
            assert_eq!(val, 0.0, "strains[0][{k}] should be 0");
        }
    }

    #[test]
    fn test_compute_stresses_identity_d() {
        // D = I_6 → σ = ε
        let mut d = [0.0f64; 36];
        for i in 0..6 {
            d[i * 6 + i] = 1.0;
        }
        let eps = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let stresses = StaticAnalysis::compute_stresses(&[eps], &d);
        assert_eq!(stresses[0], eps);
    }

    // ---- ModalAnalysis tests ----

    #[test]
    fn test_modal_1dof() {
        // K = [k], M = [m] → ω = sqrt(k/m)
        let k_val = 400.0f64;
        let m_val = 1.0f64;
        let k = vec![k_val];
        let m = vec![m_val];
        let _ma = ModalAnalysis::default();
        let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 1, 1, 200)
            .expect("1-DOF modal analysis");
        let expected = k_val.sqrt();
        assert!(
            (freqs[0] - expected).abs() / expected < 1e-5,
            "ω = {}, expected {expected}",
            freqs[0]
        );
    }

    #[test]
    fn test_modal_2dof_diagonal() {
        // K = diag(1, 9), M = I → ω₁=1, ω₂=3
        let k = vec![1.0, 0.0, 0.0, 9.0];
        let m = vec![1.0, 0.0, 0.0, 1.0];
        let _ma = ModalAnalysis::default();
        // No sorting here: the frequencies must already come out ascending.
        let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 2, 2, 500)
            .expect("2-DOF modal analysis");
        assert!((freqs[0] - 1.0).abs() < 1e-12, "ω₁ = {}", freqs[0]);
        assert!((freqs[1] - 3.0).abs() < 1e-12, "ω₂ = {}", freqs[1]);
    }

    #[test]
    fn test_mode_shape_1dof() {
        let k = vec![25.0_f64];
        let m = vec![1.0_f64];
        // ω = 5 is exactly the natural frequency, so K − ω²·M is singular.
        let phi =
            ModalAnalysis::compute_mode_shapes(&k, &m, 1, 5.0, 100).expect("1-DOF mode shape");
        // Mass-normalised: φᵀ·M·φ = 1 → |φ| = 1
        let norm: f64 = phi.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-6, "mass norm = {norm}");
    }

    // ---- LinearStaticAnalysis tests ----

    #[test]
    fn test_cantilever_deflection() {
        let length: f64 = 0.5;
        let width: f64 = 0.1;
        let height: f64 = 0.1;
        let e: f64 = 200.0e9;
        let nu: f64 = 0.3;
        let p: f64 = -500.0;

        let mesh = TetrahedralMesh::generate_beam(length, width, height, 4, 1, 1);
        let material = LinearElasticMaterial::new(e, nu);

        let mut dirichlet_bcs = Vec::new();
        for (i, node) in mesh.nodes.iter().enumerate() {
            if node.x.abs() < 1e-10 {
                dirichlet_bcs.push(DirichletBc::new(i, 0, 0.0));
                dirichlet_bcs.push(DirichletBc::new(i, 1, 0.0));
                dirichlet_bcs.push(DirichletBc::new(i, 2, 0.0));
            }
        }

        let free_end_nodes: Vec<usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| (n.x - length).abs() < 1e-10)
            .map(|(i, _)| i)
            .collect();
        let num_free = free_end_nodes.len() as f64;
        let neumann_bcs: Vec<NeumannBc> = free_end_nodes
            .iter()
            .map(|&i| NeumannBc::new(i, [0.0, p / num_free, 0.0]))
            .collect();

        let analysis = LinearStaticAnalysis {
            max_iter: 20000,
            tolerance: 1e-10,
        };
        let result = analysis
            .try_solve(
                &mesh,
                &material,
                &dirichlet_bcs,
                &neumann_bcs,
                &Vec3::new(0.0, 0.0, 0.0),
            )
            .expect("cantilever analysis must converge");

        let max_y_disp = free_end_nodes
            .iter()
            .map(|&i| result.displacements[i].y)
            .fold(f64::INFINITY, f64::min);

        assert!(
            max_y_disp < 0.0,
            "Tip deflection should be downward (negative y), got {max_y_disp:.6e}"
        );
    }

    // ---- StressInvariants ----

    #[test]
    fn test_hydrostatic_pure_hydrostatic() {
        let sigma = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        assert!(
            (s.hydrostatic() - 1.0).abs() < 1e-12,
            "hydrostatic={}",
            s.hydrostatic()
        );
    }

    #[test]
    fn test_von_mises_zero_for_hydrostatic() {
        let sigma = [5.0, 5.0, 5.0, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        assert!(
            s.von_mises() < 1e-10,
            "von Mises should be 0 for hydrostatic: {}",
            s.von_mises()
        );
    }

    #[test]
    fn test_von_mises_uniaxial() {
        let s0 = 300e6;
        let sigma = [s0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        assert!(
            (s.von_mises() - s0).abs() < 1.0,
            "von Mises uniaxial: {}",
            s.von_mises()
        );
    }

    #[test]
    fn test_principal_stresses_sorted() {
        let sigma = [200e6, 100e6, 50e6, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        let [s1, s2, s3] = s.principal_stresses();
        assert!(
            s1 >= s2 && s2 >= s3,
            "principal stresses not sorted: {s1} {s2} {s3}"
        );
    }

    #[test]
    fn test_safety_factor_von_mises() {
        let sigma = [200e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        let sf = s.safety_factor_von_mises(400e6);
        assert!((sf - 2.0).abs() < 1e-6, "safety factor={sf}");
    }

    #[test]
    fn test_lode_angle_uniaxial_tension() {
        let sigma = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        // Lode angle for uniaxial tension = 30° (π/6)
        let theta = s.lode_angle_deg();
        assert!(theta.is_finite(), "Lode angle should be finite: {theta}");
    }

    #[test]
    fn test_tresca_criterion_uniaxial() {
        let sigma = [200e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s = StressState::from_voigt(sigma);
        // Tresca = max principal - min principal
        let tresca = s.tresca_stress();
        assert!((tresca - 200e6).abs() < 1.0, "Tresca={tresca}");
    }

    // ---- FatigueLifeEstimator ----

    #[test]
    fn test_fatigue_life_below_endurance() {
        let est = FatigueLifeEstimator::new(500e6, 1e6, 200e6);
        let life = est.cycles_to_failure(100e6);
        assert!(
            life == f64::INFINITY,
            "below endurance should return infinity: {life}"
        );
    }

    #[test]
    fn test_fatigue_life_at_ultimate() {
        let est = FatigueLifeEstimator::new(500e6, 1e6, 200e6);
        let life = est.cycles_to_failure(500e6);
        assert!(life <= 1e6 + 1.0, "at UTS life should be <= N_ref: {life}");
    }

    #[test]
    fn test_damage_accumulation_miner() {
        let d = miner_damage(&[0.3, 0.5, 0.2]);
        assert!((d - 1.0).abs() < 1e-12, "Miner sum = {d}");
    }

    // ---- CantileverBeam ----

    #[test]
    fn test_cantilever_beam() {
        let length: f64 = 1.0;
        let width: f64 = 0.1;
        let height: f64 = 0.1;
        let e: f64 = 200.0e9;
        let nu: f64 = 0.3;
        let p: f64 = -1000.0;

        let inertia = width * height.powi(3) / 12.0;
        let analytical_deflection = p * length.powi(3) / (3.0 * e * inertia);

        let nx = 10;
        let ny = 2;
        let nz = 2;
        let mesh = TetrahedralMesh::generate_beam(length, width, height, nx, ny, nz);
        let material = LinearElasticMaterial::new(e, nu);

        let mut dirichlet_bcs = Vec::new();
        for (i, node) in mesh.nodes.iter().enumerate() {
            if node.x.abs() < 1e-10 {
                dirichlet_bcs.push(DirichletBc::new(i, 0, 0.0));
                dirichlet_bcs.push(DirichletBc::new(i, 1, 0.0));
                dirichlet_bcs.push(DirichletBc::new(i, 2, 0.0));
            }
        }

        let free_end_nodes: Vec<usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| (n.x - length).abs() < 1e-10)
            .map(|(i, _)| i)
            .collect();
        let num_free = free_end_nodes.len() as f64;
        let neumann_bcs: Vec<NeumannBc> = free_end_nodes
            .iter()
            .map(|&i| NeumannBc::new(i, [0.0, p / num_free, 0.0]))
            .collect();

        let analysis = LinearStaticAnalysis {
            max_iter: 50000,
            tolerance: 1e-12,
        };

        let result = analysis
            .try_solve(
                &mesh,
                &material,
                &dirichlet_bcs,
                &neumann_bcs,
                &Vec3::new(0.0, 0.0, 0.0),
            )
            .expect("cantilever analysis must converge");

        let max_y_disp = free_end_nodes
            .iter()
            .map(|&i| result.displacements[i].y)
            .fold(f64::INFINITY, f64::min);

        let error = ((max_y_disp - analytical_deflection) / analytical_deflection).abs();

        assert!(
            max_y_disp < 0.0,
            "Deflection should be negative (downward), got {max_y_disp:.6e}",
        );
        assert!(
            error < 1.0,
            "Cantilever deflection error too large: {:.1}%. \
             FEM = {max_y_disp:.6e}, analytical = {analytical_deflection:.6e}",
            error * 100.0,
        );
    }
}
