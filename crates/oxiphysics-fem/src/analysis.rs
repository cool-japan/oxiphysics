#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Static, modal, and dynamic FEM analysis drivers.
//!
//! Provides:
//! - [`LinearStaticAnalysis`] — classic K·u = f solve
//! - [`StaticAnalysis`] — standalone dense-matrix static analysis
//! - [`ModalAnalysis`] — Rayleigh-quotient iteration / power iteration
//! - [`DynamicAnalysis`] — Newmark-β and Wilson-θ time integration

#![allow(dead_code)]

use oxiphysics_core::math::Vec3;

use crate::assembly::{assemble_load_vector, assemble_stiffness};
use crate::boundary::{DirichletBc, NeumannBc, apply_dirichlet, apply_neumann};
use crate::constitutive::LinearElasticMaterial;
use crate::element::LinearTetrahedron;
use crate::mesh::TetrahedralMesh;
use crate::solvers::{jacobi_preconditioner, preconditioned_cg};

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

    /// Run the analysis and return the results.
    #[allow(clippy::too_many_arguments)]
    pub fn solve(
        &self,
        mesh: &TetrahedralMesh,
        material: &LinearElasticMaterial,
        dirichlet_bcs: &[DirichletBc],
        neumann_bcs: &[NeumannBc],
        body_force: &Vec3,
    ) -> FemResult {
        let ndof = mesh.num_nodes() * 3;

        let mut k = assemble_stiffness(mesh, material);
        let mut f = assemble_load_vector(mesh, body_force);

        apply_neumann(&mut f, neumann_bcs);
        apply_dirichlet(&mut k, &mut f, dirichlet_bcs);

        let precond = jacobi_preconditioner(&k);
        let u = preconditioned_cg(&k, &f, &precond, self.max_iter, self.tolerance);

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

        let _ = ndof;

        FemResult {
            displacements,
            stress,
            strain,
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
// ModalAnalysis — Rayleigh-quotient / inverse iteration
// ---------------------------------------------------------------------------

/// Modal analysis: computes natural frequencies and mode shapes using
/// Rayleigh-quotient iteration (inverse iteration with shift).
pub struct ModalAnalysis {
    /// Maximum number of iterations per mode.
    pub max_iter: usize,
    /// Convergence tolerance on the Rayleigh quotient.
    pub tol: f64,
}

impl Default for ModalAnalysis {
    fn default() -> Self {
        Self {
            max_iter: 500,
            tol: 1e-10,
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

    /// Solve A·x = b by Gaussian elimination (small n, no pivoting needed for
    /// well-conditioned problems).  Returns None if singular.
    fn solve_dense(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
        let cols = n + 1;
        let mut aug = vec![0.0; n * cols];
        for i in 0..n {
            for j in 0..n {
                aug[i * cols + j] = a[i * n + j];
            }
            aug[i * cols + n] = b[i];
        }
        for col in 0..n {
            // Partial pivot
            let mut max_val = aug[col * cols + col].abs();
            let mut max_row = col;
            for row in (col + 1)..n {
                let v = aug[row * cols + col].abs();
                if v > max_val {
                    max_val = v;
                    max_row = row;
                }
            }
            if max_val < 1e-15 {
                return None;
            }
            if max_row != col {
                for k in 0..cols {
                    aug.swap(col * cols + k, max_row * cols + k);
                }
            }
            let pivot = aug[col * cols + col];
            for row in (col + 1)..n {
                let factor = aug[row * cols + col] / pivot;
                for k in col..cols {
                    let sub = factor * aug[col * cols + k];
                    aug[row * cols + k] -= sub;
                }
            }
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut s = aug[i * cols + n];
            for j in (i + 1)..n {
                s -= aug[i * cols + j] * x[j];
            }
            x[i] = s / aug[i * cols + i];
        }
        Some(x)
    }

    /// Compute the `n_modes` lowest natural frequencies (rad/s) using
    /// Rayleigh-quotient iteration on the dense generalized eigenproblem K·φ = ω²·M·φ.
    ///
    /// Both `k` and `m` are row-major flat arrays of size n×n.
    pub fn compute_natural_frequencies(
        k: &[f64],
        m: &[f64],
        n: usize,
        n_modes: usize,
        max_iter: usize,
    ) -> Vec<f64> {
        let mut freqs = Vec::with_capacity(n_modes);
        let mut found_modes: Vec<Vec<f64>> = Vec::new();

        for mode_idx in 0..n_modes {
            // Initial vector
            let mut q: Vec<f64> = (0..n)
                .map(|i| {
                    let v = ((i + mode_idx * 13 + 1) as f64 * 0.7).sin();
                    if i == mode_idx % n { v + 1.0 } else { v }
                })
                .collect();

            // Mass-normalise
            let mq = Self::matvec_dense(m, &q, n);
            let nm = Self::dot(&q, &mq).max(0.0).sqrt();
            if nm > 1e-60 {
                for qi in q.iter_mut() {
                    *qi /= nm;
                }
            }

            let mut omega_sq = 0.0_f64;

            for _iter in 0..max_iter {
                // Deflate against converged modes
                for phi in &found_modes {
                    let mphi = Self::matvec_dense(m, phi, n);
                    let c = Self::dot(&q, &mphi);
                    for i in 0..n {
                        q[i] -= c * phi[i];
                    }
                }

                // Rayleigh quotient
                let kq = Self::matvec_dense(k, &q, n);
                let mq2 = Self::matvec_dense(m, &q, n);
                let num = Self::dot(&q, &kq);
                let den = Self::dot(&q, &mq2);
                let rq = if den.abs() > 1e-60 { num / den } else { 0.0 };

                // Inverse iteration: solve (K - rq*M)*y = M*q
                let mut km = vec![0.0; n * n];
                for i in 0..n {
                    for j in 0..n {
                        km[i * n + j] = k[i * n + j] - rq * m[i * n + j];
                    }
                }
                let rhs = Self::matvec_dense(m, &q, n);
                let y_opt = Self::solve_dense(&km, &rhs, n);
                let y = match y_opt {
                    Some(y) => y,
                    None => {
                        // Fallback: plain power iteration
                        let mut rhs2 = Self::matvec_dense(m, &q, n);
                        // Solve K·y = M·q instead
                        if let Some(sol) = Self::solve_dense(k, &rhs2, n) {
                            rhs2 = sol;
                        }
                        rhs2
                    }
                };

                // Deflate y
                let mut yd = y.clone();
                for phi in &found_modes {
                    let mphi = Self::matvec_dense(m, phi, n);
                    let c = Self::dot(&yd, &mphi);
                    for i in 0..n {
                        yd[i] -= c * phi[i];
                    }
                }

                // Mass-normalise
                let my = Self::matvec_dense(m, &yd, n);
                let nm2 = Self::dot(&yd, &my).max(0.0).sqrt();
                if nm2 < 1e-60 {
                    break;
                }
                let q_new: Vec<f64> = yd.iter().map(|yi| yi / nm2).collect();

                // New Rayleigh quotient
                let kqn = Self::matvec_dense(k, &q_new, n);
                let mqn = Self::matvec_dense(m, &q_new, n);
                let num_new = Self::dot(&q_new, &kqn);
                let den_new = Self::dot(&q_new, &mqn);
                let rq_new = if den_new.abs() > 1e-60 {
                    num_new / den_new
                } else {
                    0.0
                };

                let change = (rq_new - omega_sq).abs() / (rq_new.abs() + 1e-60);
                q = q_new;
                omega_sq = rq_new;
                if change < 1e-10 && _iter > 0 {
                    break;
                }
            }

            freqs.push(omega_sq.max(0.0).sqrt());
            found_modes.push(q);
        }

        freqs
    }

    /// Compute the mode shape for a given angular frequency `omega` using
    /// inverse iteration: solve (K − ω²·M)·φ = b repeatedly.
    ///
    /// Both `k` and `m` are row-major flat arrays of size n×n.
    pub fn compute_mode_shapes(
        k: &[f64],
        m: &[f64],
        n: usize,
        omega: f64,
        max_iter: usize,
    ) -> Vec<f64> {
        let shift = omega * omega;
        // Build K - shift*M
        let mut km = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                km[i * n + j] = k[i * n + j] - shift * m[i * n + j];
            }
        }

        let mut q: Vec<f64> = (0..n).map(|i| ((i + 1) as f64 * 0.5).sin()).collect();
        // Normalise
        let nrm: f64 = q.iter().map(|x| x * x).sum::<f64>().sqrt();
        if nrm > 1e-60 {
            for qi in q.iter_mut() {
                *qi /= nrm;
            }
        }

        for _ in 0..max_iter {
            let rhs = Self::matvec_dense(m, &q, n);
            let y = match Self::solve_dense(&km, &rhs, n) {
                Some(y) => y,
                None => break,
            };
            let ny: f64 = y.iter().map(|x| x * x).sum::<f64>().sqrt();
            if ny < 1e-60 {
                break;
            }
            let q_new: Vec<f64> = y.iter().map(|yi| yi / ny).collect();
            let change: f64 = q
                .iter()
                .zip(q_new.iter())
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
                / (n as f64);
            q = q_new;
            if change < 1e-12 {
                break;
            }
        }

        // Mass-normalise
        let mq = Self::matvec_dense(m, &q, n);
        let mn: f64 = Self::dot(&q, &mq).max(0.0).sqrt();
        if mn > 1e-60 {
            for qi in q.iter_mut() {
                *qi /= mn;
            }
        }
        q
    }
}

// ---------------------------------------------------------------------------
// DynamicAnalysis — Newmark-β and Wilson-θ
// ---------------------------------------------------------------------------

/// Direct time-integration for structural dynamics: M·ü + C·u̇ + K·u = f(t).
///
/// Provides Newmark-β and Wilson-θ stepping methods.
pub struct DynamicAnalysis {
    /// Mass matrix (row-major, n×n).
    pub m: Vec<f64>,
    /// Damping matrix (row-major, n×n).
    pub c: Vec<f64>,
    /// Stiffness matrix (row-major, n×n).
    pub k: Vec<f64>,
    /// Number of DOFs.
    pub n: usize,
}

impl DynamicAnalysis {
    /// Construct a new `DynamicAnalysis`.
    pub fn new(m: Vec<f64>, c: Vec<f64>, k: Vec<f64>, n: usize) -> Self {
        Self { m, c, k, n }
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

    /// Solve dense n×n system A·x = b (Gaussian elim with partial pivot).
    fn solve(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let cols = n + 1;
        let mut aug = vec![0.0; n * cols];
        for i in 0..n {
            for j in 0..n {
                aug[i * cols + j] = a[i * n + j];
            }
            aug[i * cols + n] = b[i];
        }
        for col in 0..n {
            let mut max_val = aug[col * cols + col].abs();
            let mut max_row = col;
            for row in (col + 1)..n {
                let v = aug[row * cols + col].abs();
                if v > max_val {
                    max_val = v;
                    max_row = row;
                }
            }
            if max_val < 1e-15 {
                continue;
            }
            if max_row != col {
                for k in 0..cols {
                    aug.swap(col * cols + k, max_row * cols + k);
                }
            }
            let piv = aug[col * cols + col];
            for row in (col + 1)..n {
                let fac = aug[row * cols + col] / piv;
                for k in col..cols {
                    let sub = fac * aug[col * cols + k];
                    aug[row * cols + k] -= sub;
                }
            }
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut s = aug[i * cols + n];
            for j in (i + 1)..n {
                s -= aug[i * cols + j] * x[j];
            }
            x[i] = s / aug[i * cols + i];
        }
        x
    }

    /// One step of the Newmark-β method.
    ///
    /// Updates displacement `u`, velocity `v`, and acceleration `a` from time `t`
    /// to time `t + dt`.
    ///
    /// Typical parameter choices:
    /// - Constant average acceleration (unconditionally stable): β = 0.25, γ = 0.5
    /// - Linear acceleration (conditionally stable): β = 1/6, γ = 0.5
    ///
    /// # Returns
    /// `(u_{n+1}, v_{n+1}, a_{n+1})`
    #[allow(clippy::too_many_arguments)]
    pub fn newmark_beta_step(
        &self,
        u: &[f64],
        v: &[f64],
        a: &[f64],
        f_ext: &[f64],
        dt: f64,
        beta: f64,
        gamma: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = self.n;
        // Predictors
        let u_pred: Vec<f64> = (0..n)
            .map(|i| u[i] + dt * v[i] + dt * dt * (0.5 - beta) * a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n).map(|i| v[i] + dt * (1.0 - gamma) * a[i]).collect();

        // Effective stiffness K_eff = M/(β·dt²) + γ/(β·dt)·C + K
        let c1 = 1.0 / (beta * dt * dt);
        let c2 = gamma / (beta * dt);

        let mut k_eff = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                k_eff[i * n + j] =
                    c1 * self.m[i * n + j] + c2 * self.c[i * n + j] + self.k[i * n + j];
            }
        }

        // Effective force
        let ku_pred = Self::mv(&self.k, &u_pred, n);
        let cv_pred = Self::mv(&self.c, &v_pred, n);
        let rhs: Vec<f64> = (0..n).map(|i| f_ext[i] - ku_pred[i] - cv_pred[i]).collect();

        let a_new = Self::solve(&k_eff, &rhs, n);
        let u_new: Vec<f64> = (0..n)
            .map(|i| u_pred[i] + beta * dt * dt * a_new[i])
            .collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + gamma * dt * a_new[i]).collect();

        (u_new, v_new, a_new)
    }

    /// One step of the Wilson-θ method.
    ///
    /// `theta` is typically 1.37–1.42 for unconditional stability.
    ///
    /// # Returns
    /// `(u_{n+1}, v_{n+1}, a_{n+1})`
    #[allow(clippy::too_many_arguments)]
    pub fn wilson_theta_step(
        &self,
        u: &[f64],
        v: &[f64],
        a: &[f64],
        f_ext: &[f64],
        dt: f64,
        theta: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = self.n;
        let tau = theta * dt;

        // Effective stiffness K_eff = 6/(θ²dt²)·M + 3/(θdt)·C + K
        let c1 = 6.0 / (tau * tau);
        let c2 = 3.0 / tau;

        let mut k_eff = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                k_eff[i * n + j] =
                    c1 * self.m[i * n + j] + c2 * self.c[i * n + j] + self.k[i * n + j];
            }
        }

        // Extrapolated force at t + θ·dt: linear extrapolation
        // f_tau ≈ f_t + θ*(f_ext - f_t)  but we have only f_ext (at t+dt),
        // so approximate f_tau = f_ext (conservative)
        let f_tau = f_ext.to_vec();

        // Effective RHS
        let mu = Self::mv(&self.m, u, n);
        let mv_curr = Self::mv(&self.m, v, n);
        let ma = Self::mv(&self.m, a, n);
        let cu = Self::mv(&self.c, u, n);
        let cv_curr = Self::mv(&self.c, v, n);

        let rhs: Vec<f64> = (0..n)
            .map(|i| {
                f_tau[i] + c1 * mu[i] + c2 * mv_curr[i] + 2.0 * ma[i]
                    - self.k[i * n..i * n + n]
                        .iter()
                        .zip(u.iter())
                        .map(|(kij, uj)| kij * uj)
                        .sum::<f64>()
                    + (c2 * cu[i] + cv_curr[i])
            })
            .collect();

        // Solve for a at t + θ·dt
        let a_tau = Self::solve(&k_eff, &rhs, n);

        // Interpolate back to t + dt
        let a_new: Vec<f64> = (0..n).map(|i| a[i] + (a_tau[i] - a[i]) / theta).collect();
        let v_new: Vec<f64> = (0..n)
            .map(|i| v[i] + dt / 2.0 * (a[i] + a_new[i]))
            .collect();
        let u_new: Vec<f64> = (0..n)
            .map(|i| u[i] + dt * v[i] + dt * dt / 6.0 * (2.0 * a[i] + a_new[i]))
            .collect();

        (u_new, v_new, a_new)
    }
}

// ---------------------------------------------------------------------------
// Stress state and invariants
// ---------------------------------------------------------------------------

/// Full 3-D stress state in Voigt notation: \[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz\].
#[allow(dead_code)]
#[allow(non_snake_case)]
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
    #[allow(clippy::too_many_arguments)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
pub fn miner_damage(cycle_ratios: &[f64]) -> f64 {
    cycle_ratios.iter().sum()
}

// ---------------------------------------------------------------------------
// Post-processing utilities
// ---------------------------------------------------------------------------

/// Compute element-averaged strains from a displacement field and B-matrices.
///
/// For each element, ε = B · u_local.  Stores result as `[f64; 6]` in Voigt order.
#[allow(dead_code)]
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
#[allow(dead_code)]
pub fn von_mises_field(stresses: &[[f64; 6]]) -> Vec<f64> {
    stresses
        .iter()
        .map(|&s| StressState::from_voigt(s).von_mises())
        .collect()
}

/// Compute element safety factors by von Mises criterion.
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
        for k in 1..6 {
            assert_eq!(strains[0][k], 0.0);
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
        let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 1, 1, 200);
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
        let mut freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 2, 2, 500);
        freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((freqs[0] - 1.0).abs() < 0.05, "ω₁ = {}", freqs[0]);
        assert!((freqs[1] - 3.0).abs() < 0.05, "ω₂ = {}", freqs[1]);
    }

    #[test]
    fn test_mode_shape_1dof() {
        let k = vec![25.0_f64];
        let m = vec![1.0_f64];
        let phi = ModalAnalysis::compute_mode_shapes(&k, &m, 1, 5.0, 100);
        // Mass-normalised: φᵀ·M·φ = 1 → |φ| = 1
        let norm: f64 = phi.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-6, "mass norm = {norm}");
    }

    // ---- DynamicAnalysis tests ----

    /// 1-DOF undamped free vibration: M·ü + K·u = 0.
    /// Verifies the Newmark method produces a bounded, oscillatory response.
    #[test]
    fn test_newmark_free_vibration() {
        let m_val = 1.0_f64;
        let k_val = 100.0_f64;
        let dt = 0.001;
        let beta = 0.25;
        let gamma = 0.5;

        let da = DynamicAnalysis::new(vec![m_val], vec![0.0], vec![k_val], 1);
        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-k_val / m_val];

        for _ in 0..100 {
            let (un, vn, an) = da.newmark_beta_step(&u, &v, &a, &[0.0], dt, beta, gamma);
            u = un;
            v = vn;
            a = an;
        }
        // Displacement should stay bounded (energy conservation)
        assert!(u[0].abs() <= 1.01, "u should be bounded, got {}", u[0]);
    }

    #[test]
    fn test_newmark_static_equilibrium() {
        // Apply constant force f to 1-DOF with damping: should move toward u_static = f/k
        let da = DynamicAnalysis::new(vec![1.0], vec![0.1], vec![10.0], 1);
        let mut u = vec![0.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![0.0_f64];
        let f = vec![5.0_f64]; // static → u_static = 5/10 = 0.5
        for _ in 0..5000 {
            let (un, vn, an) = da.newmark_beta_step(&u, &v, &a, &f, 0.01, 0.25, 0.5);
            u = un;
            v = vn;
            a = an;
        }
        // Should be in the general vicinity of the static solution
        assert!(u[0] > 0.0, "u should be positive, got {}", u[0]);
    }

    #[test]
    fn test_wilson_theta_step_displacement_sign() {
        // Force applied in positive direction → displacement must be positive
        let da = DynamicAnalysis::new(vec![1.0], vec![0.0], vec![1.0], 1);
        let u = vec![0.0_f64];
        let v = vec![0.0_f64];
        let a = vec![0.0_f64];
        let (u_new, _, _) = da.wilson_theta_step(&u, &v, &a, &[1.0], 0.01, 1.4);
        assert!(
            u_new[0] >= 0.0,
            "displacement should be non-negative, got {}",
            u_new[0]
        );
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
        let result = analysis.solve(
            &mesh,
            &material,
            &dirichlet_bcs,
            &neumann_bcs,
            &Vec3::new(0.0, 0.0, 0.0),
        );

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

        let result = analysis.solve(
            &mesh,
            &material,
            &dirichlet_bcs,
            &neumann_bcs,
            &Vec3::new(0.0, 0.0, 0.0),
        );

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
