// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Poromechanics and Biot's consolidation theory.
//!
//! Implements coupled solid-fluid mechanics for saturated porous media:
//!
//! - [`BiotCoeff`]: Biot effective stress coefficient and modulus
//! - [`PorousElement`]: Finite element for porous media with coupling matrix
//! - [`ConsolidationSolver`]: Sparse backward-Euler u-p solver for Biot consolidation
//! - [`TerzaghiConsolidation`]: 1-D Terzaghi consolidation (analytical)
//! - [`DarcyFlow`]: Darcy seepage velocity computation
//! - [`terzaghi_u`]: Settlement degree of consolidation series expansion
//! - [`consolidation_time`]: Real time from dimensionless time factor
//! - [`effective_stress`]: Terzaghi effective stress principle

use std::f64::consts::PI;

use crate::error::{Error, Result};
use crate::sparse::CsrMatrix;

mod spd_solve;
use spd_solve::{SpdPreconditioner, csr_triplets_checked, mask_matrix, preconditioned_cg};

// ============================================================================
// BiotCoeff
// ============================================================================

/// Biot effective-stress coefficient and Biot modulus for a saturated porous medium.
///
/// The Biot coefficient α relates total stress changes to pore-pressure changes:
/// σ′ = σ - α p
///
/// The Biot modulus M controls the compressibility of the fluid-solid mixture.
pub struct BiotCoeff {
    /// Biot effective-stress coefficient α (dimensionless, 0 < α ≤ 1).
    pub alpha: f64,
    /// Biot modulus M (Pa) – inverse of the storage coefficient.
    pub modulus: f64,
}

impl BiotCoeff {
    /// Create a `BiotCoeff` with the given α and M values.
    pub fn new(alpha: f64, modulus: f64) -> Self {
        Self { alpha, modulus }
    }

    /// Compute Biot parameters from drained bulk modulus, fluid bulk modulus and porosity.
    ///
    /// Uses the Biot-Gassmann relations:
    /// α = 1 - K_dry / K_s
    /// 1/M = φ/K_f + (α - φ)/K_s
    ///
    /// # Arguments
    /// * `k_dry` – drained (skeleton) bulk modulus \[Pa\]
    /// * `k_s`   – solid grain bulk modulus \[Pa\]
    /// * `k_f`   – pore fluid bulk modulus \[Pa\]
    /// * `phi`   – porosity (0 < φ < 1)
    pub fn compute_from_bulk(k_dry: f64, k_s: f64, k_f: f64, phi: f64) -> Self {
        let alpha = 1.0 - k_dry / k_s;
        let inv_m = phi / k_f + (alpha - phi) / k_s;
        let modulus = if inv_m.abs() < 1e-300 {
            f64::INFINITY
        } else {
            1.0 / inv_m
        };
        Self { alpha, modulus }
    }

    /// Storage coefficient S = α² / M + φ / K_f (Pa⁻¹).
    pub fn storage_coefficient(&self, phi: f64, k_f: f64) -> f64 {
        self.alpha * self.alpha / self.modulus + phi / k_f
    }

    /// Undrained bulk modulus K_u = K_dry + α² M.
    pub fn undrained_bulk(k_dry: f64, alpha: f64, modulus: f64) -> f64 {
        k_dry + alpha * alpha * modulus
    }

    /// Skempton's B coefficient B = α M / K_u.
    pub fn skempton_b(k_dry: f64, alpha: f64, modulus: f64) -> f64 {
        let k_u = Self::undrained_bulk(k_dry, alpha, modulus);
        if k_u.abs() < 1e-300 {
            0.0
        } else {
            alpha * modulus / k_u
        }
    }
}

// ============================================================================
// PorousElement
// ============================================================================

/// Finite element representing a porous solid-fluid element.
///
/// Stores connectivity, permeability tensor, porosity, and solid stiffness
/// for assembly of the coupled u-p system.
pub struct PorousElement {
    /// Global node IDs for this element.
    pub node_ids: Vec<usize>,
    /// Permeability tensor k \[m²\] stored in row-major 3×3 order.
    pub permeability: [f64; 9],
    /// Porosity φ (dimensionless).
    pub porosity: f64,
    /// Drained Young's modulus E \[Pa\].
    pub solid_stiffness: f64,
}

impl PorousElement {
    /// Create a new `PorousElement`.
    ///
    /// # Arguments
    /// * `node_ids`       – global DOF node indices
    /// * `permeability`   – 3×3 permeability tensor in row-major order \[m²\]
    /// * `porosity`       – porosity φ
    /// * `solid_stiffness`– drained Young's modulus E \[Pa\]
    pub fn new(
        node_ids: Vec<usize>,
        permeability: [f64; 9],
        porosity: f64,
        solid_stiffness: f64,
    ) -> Self {
        Self {
            node_ids,
            permeability,
            porosity,
            solid_stiffness,
        }
    }

    /// Compute a simplified coupling matrix L = α B^T m N_p dV.
    ///
    /// Returns a flat vector of length `3 * n_nodes` representing coupling
    /// contributions for a single Gauss point with unit volume and Biot α = 1.
    pub fn compute_coupling_matrix(&self) -> Vec<f64> {
        let n = self.node_ids.len();
        // Each displacement node has 3 DOFs; each pressure node has 1 DOF.
        // Simplified: return B^T m shape for unit Jacobian, α=1.
        let mut coupling = vec![0.0_f64; 3 * n];
        let n_f64 = n as f64;
        for i in 0..n {
            // Uniform shape function gradient contribution (1/n per direction)
            coupling[3 * i] = 1.0 / n_f64;
            coupling[3 * i + 1] = 1.0 / n_f64;
            coupling[3 * i + 2] = 1.0 / n_f64;
        }
        coupling
    }

    /// Permeability in the x-direction k_xx \[m²\].
    pub fn k_xx(&self) -> f64 {
        self.permeability[0]
    }

    /// Isotropic permeability scalar (mean of diagonal entries).
    pub fn isotropic_permeability(&self) -> f64 {
        (self.permeability[0] + self.permeability[4] + self.permeability[8]) / 3.0
    }
}

// ============================================================================
// ConsolidationSolver
// ============================================================================

/// Block of the discrete Biot system addressed by [`ConsolidationSolver`]
/// assembly calls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BiotBlock {
    /// Drained skeleton stiffness `K_uu = ∫ Bᵀ D B dV` (`n_u × n_u`).
    ///
    /// Must be symmetric and positive definite on the free displacement DOFs.
    Stiffness,
    /// Coupling matrix `Q = ∫ α Bᵀ m N_p dV` (`n_u × n_p`).
    Coupling,
    /// Storage matrix `S = ∫ N_pᵀ (1/M) N_p dV` (`n_p × n_p`, symmetric positive
    /// semi-definite; zero for incompressible constituents).
    Storage,
    /// Flow matrix `H = ∫ ∇N_pᵀ (k/μ) ∇N_p dV` (`n_p × n_p`, symmetric positive
    /// semi-definite).
    Permeability,
}

/// Iteration statistics of one [`ConsolidationSolver::step`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConsolidationStepReport {
    /// Iterations of the outer preconditioned CG on the pressure Schur complement.
    pub outer_iterations: usize,
    /// Total iterations of the inner preconditioned CG solves with `K_uu`.
    pub inner_iterations: usize,
    /// Final relative residual of the pressure Schur-complement system.
    pub relative_residual: f64,
}

/// Time-domain solver for Biot's quasi-static consolidation equations.
///
/// The semi-discrete u-p system is
///
/// ```text
/// K_uu u        − Q p        = f_u      (equilibrium)
/// Qᵀ u̇ + S ṗ + H p           = f_p      (fluid mass balance)
/// ```
///
/// with skeleton stiffness `K_uu`, coupling `Q = ∫ α Bᵀ m N_p dV`, storage `S`
/// and flow (permeability) matrix `H`.  Backward-Euler integration gives, per
/// step of size `Δt`,
///
/// ```text
/// [ K_uu   −Q        ] { u_{n+1} }   { f_u                           }
/// [ Qᵀ     S + Δt·H  ] { p_{n+1} } = { Qᵀ u_n + S p_n + Δt·f_p        }
/// ```
///
/// The matrices are stored sparse ([`CsrMatrix`]) and assembled through
/// [`ConsolidationSolver::add_triplets`], [`ConsolidationSolver::add_element_block`]
/// or [`ConsolidationSolver::set_block`].  Dirichlet conditions are imposed
/// exactly by elimination ([`ConsolidationSolver::fix_displacement`],
/// [`ConsolidationSolver::fix_pressure`]).
///
/// Each step eliminates the displacements and solves the symmetric positive
/// definite pressure Schur complement
/// `Σ = S + Δt·H + Qᵀ K_uu⁻¹ Q` with preconditioned conjugate gradients
/// (diagonal approximation of `Σ` as preconditioner); every application of
/// `K_uu⁻¹` is an inner preconditioned CG with an incomplete-Cholesky (IC(0),
/// diagonally shifted on breakdown) preconditioner.  `Δt = 0` is admissible and
/// yields the instantaneous undrained response.
#[derive(Debug, Clone)]
pub struct ConsolidationSolver {
    /// Displacement DOF vector `u_n`.
    pub u_dofs: Vec<f64>,
    /// Pore-pressure DOF vector `p_n`.
    pub p_dofs: Vec<f64>,
    /// Current simulation time \[s\].
    pub time: f64,
    /// Relative residual tolerance of the outer Schur-complement CG.
    pub tolerance: f64,
    /// Maximum iterations of each (outer and inner) CG solve.
    pub max_iterations: usize,
    stiffness: CsrMatrix,
    coupling: CsrMatrix,
    storage: CsrMatrix,
    permeability: CsrMatrix,
    u_constraints: Vec<Option<f64>>,
    p_constraints: Vec<Option<f64>>,
    /// Coordinate entries added since the last flush, per [`BiotBlock`].
    pending: [Vec<(usize, usize, f64)>; 4],
    prepared: Option<PreparedBiot>,
}

impl ConsolidationSolver {
    /// Create a new `ConsolidationSolver` with empty (all-zero) sparse matrices.
    ///
    /// # Arguments
    /// * `n_u`  – number of displacement DOFs
    /// * `n_p`  – number of pressure DOFs
    pub fn new(n_u: usize, n_p: usize) -> Self {
        Self {
            u_dofs: vec![0.0; n_u],
            p_dofs: vec![0.0; n_p],
            time: 0.0,
            tolerance: 1e-10,
            max_iterations: 10_000,
            stiffness: CsrMatrix::new(n_u, n_u),
            coupling: CsrMatrix::new(n_u, n_p),
            storage: CsrMatrix::new(n_p, n_p),
            permeability: CsrMatrix::new(n_p, n_p),
            u_constraints: vec![None; n_u],
            p_constraints: vec![None; n_p],
            pending: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            prepared: None,
        }
    }

    /// Number of displacement DOFs.
    pub fn n_u(&self) -> usize {
        self.u_dofs.len()
    }

    /// Number of pressure DOFs.
    pub fn n_p(&self) -> usize {
        self.p_dofs.len()
    }

    /// Required `(rows, cols)` of a block.
    pub fn block_shape(&self, block: BiotBlock) -> (usize, usize) {
        let n_u = self.n_u();
        let n_p = self.n_p();
        match block {
            BiotBlock::Stiffness => (n_u, n_u),
            BiotBlock::Coupling => (n_u, n_p),
            BiotBlock::Storage | BiotBlock::Permeability => (n_p, n_p),
        }
    }

    /// The assembled sparse matrix of a block (pending entries added through
    /// [`ConsolidationSolver::add_triplets`] are merged first).
    pub fn block(&mut self, block: BiotBlock) -> &CsrMatrix {
        self.flush_block(block);
        match block {
            BiotBlock::Stiffness => &self.stiffness,
            BiotBlock::Coupling => &self.coupling,
            BiotBlock::Storage => &self.storage,
            BiotBlock::Permeability => &self.permeability,
        }
    }

    fn block_index(block: BiotBlock) -> usize {
        match block {
            BiotBlock::Stiffness => 0,
            BiotBlock::Coupling => 1,
            BiotBlock::Storage => 2,
            BiotBlock::Permeability => 3,
        }
    }

    /// Merge the pending coordinate entries of `block` into its CSR storage
    /// (one `O(nnz log nnz)` rebuild, duplicates summed).
    fn flush_block(&mut self, block: BiotBlock) {
        let pending = std::mem::take(&mut self.pending[Self::block_index(block)]);
        if pending.is_empty() {
            return;
        }
        let (rows, cols) = self.block_shape(block);
        let current = self.block_mut(block);
        let mut merged = Vec::with_capacity(current.nnz() + pending.len());
        for row in 0..current.nrows {
            for idx in current.row_ptr[row]..current.row_ptr[row + 1] {
                merged.push((row, current.col_indices[idx], current.values[idx]));
            }
        }
        merged.extend(pending);
        *current = CsrMatrix::from_triplets(rows, cols, &merged);
    }

    fn block_mut(&mut self, block: BiotBlock) -> &mut CsrMatrix {
        match block {
            BiotBlock::Stiffness => &mut self.stiffness,
            BiotBlock::Coupling => &mut self.coupling,
            BiotBlock::Storage => &mut self.storage,
            BiotBlock::Permeability => &mut self.permeability,
        }
    }

    /// Replace a block by `matrix`.
    ///
    /// The matrix is validated (shape, CSR structure, finite values) and
    /// normalised to sorted, duplicate-free rows.  Entries of this block still
    /// pending from earlier [`ConsolidationSolver::add_triplets`] calls are
    /// discarded.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] on a shape mismatch, malformed CSR arrays or
    /// non-finite values.
    pub fn set_block(&mut self, block: BiotBlock, matrix: CsrMatrix) -> Result<()> {
        let (rows, cols) = self.block_shape(block);
        if matrix.nrows != rows || matrix.ncols != cols {
            return Err(Error::InvalidInput(format!(
                "{block:?} block must be {rows}×{cols}, got {}×{}",
                matrix.nrows, matrix.ncols
            )));
        }
        let triplets = csr_triplets_checked(&matrix)?;
        self.pending[Self::block_index(block)].clear();
        *self.block_mut(block) = CsrMatrix::from_triplets(rows, cols, &triplets);
        self.prepared = None;
        Ok(())
    }

    /// Add coordinate entries `(row, col, value)` to a block (duplicates are
    /// summed).  Entries are validated and buffered; they are merged into the
    /// CSR storage once, by the next [`ConsolidationSolver::step`] or
    /// [`ConsolidationSolver::block`] call, so element-by-element assembly costs
    /// `O(nnz log nnz)` in total.
    ///
    /// # Errors
    /// [`Error::DofOutOfRange`] for an index outside the block,
    /// [`Error::InvalidInput`] for a non-finite value.
    pub fn add_triplets(
        &mut self,
        block: BiotBlock,
        triplets: &[(usize, usize, f64)],
    ) -> Result<()> {
        if triplets.is_empty() {
            return Ok(());
        }
        let (rows, cols) = self.block_shape(block);
        for &(r, c, v) in triplets {
            if r >= rows {
                return Err(Error::DofOutOfRange { dof: r, n: rows });
            }
            if c >= cols {
                return Err(Error::DofOutOfRange { dof: c, n: cols });
            }
            if !v.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "non-finite value {v} at ({r}, {c}) of the {block:?} block"
                )));
            }
        }
        self.pending[Self::block_index(block)].extend_from_slice(triplets);
        self.prepared = None;
        Ok(())
    }

    /// Scatter-add a dense element matrix into a block.
    ///
    /// `local` is row-major with `row_dofs.len()` rows and `col_dofs.len()`
    /// columns; `row_dofs` / `col_dofs` are the global indices of its rows and
    /// columns in the block.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `local` has the wrong length, plus the errors
    /// of [`ConsolidationSolver::add_triplets`].
    pub fn add_element_block(
        &mut self,
        block: BiotBlock,
        row_dofs: &[usize],
        col_dofs: &[usize],
        local: &[f64],
    ) -> Result<()> {
        let n_cols = col_dofs.len();
        if local.len() != row_dofs.len() * n_cols {
            return Err(Error::InvalidInput(format!(
                "element matrix has {} entries, expected {}×{}",
                local.len(),
                row_dofs.len(),
                n_cols
            )));
        }
        let mut triplets = Vec::with_capacity(local.len());
        for (a, &r) in row_dofs.iter().enumerate() {
            for (b, &c) in col_dofs.iter().enumerate() {
                let v = local[a * n_cols + b];
                if v != 0.0 {
                    triplets.push((r, c, v));
                }
            }
        }
        self.add_triplets(block, &triplets)
    }

    /// Prescribe displacement DOF `dof` to `value` (Dirichlet condition).
    ///
    /// # Errors
    /// [`Error::DofOutOfRange`] if `dof >= n_u`, [`Error::InvalidInput`] if
    /// `value` is not finite.
    pub fn fix_displacement(&mut self, dof: usize, value: f64) -> Result<()> {
        let n = self.n_u();
        Self::fix_dof(&mut self.u_constraints, &mut self.prepared, dof, value, n)
    }

    /// Prescribe pressure DOF `dof` to `value` (drained boundary for `0.0`).
    ///
    /// # Errors
    /// [`Error::DofOutOfRange`] if `dof >= n_p`, [`Error::InvalidInput`] if
    /// `value` is not finite.
    pub fn fix_pressure(&mut self, dof: usize, value: f64) -> Result<()> {
        let n = self.n_p();
        Self::fix_dof(&mut self.p_constraints, &mut self.prepared, dof, value, n)
    }

    fn fix_dof(
        constraints: &mut [Option<f64>],
        prepared: &mut Option<PreparedBiot>,
        dof: usize,
        value: f64,
        n: usize,
    ) -> Result<()> {
        if !value.is_finite() {
            return Err(Error::InvalidInput(format!(
                "prescribed value {value} for DOF {dof} is not finite"
            )));
        }
        let slot = constraints
            .get_mut(dof)
            .ok_or(Error::DofOutOfRange { dof, n })?;
        if slot.is_none() {
            // The set of free DOFs changed: the eliminated operators are stale.
            *prepared = None;
        }
        *slot = Some(value);
        Ok(())
    }

    /// Remove every displacement and pressure Dirichlet condition.
    pub fn clear_constraints(&mut self) {
        self.u_constraints.iter_mut().for_each(|c| *c = None);
        self.p_constraints.iter_mut().for_each(|c| *c = None);
        self.prepared = None;
    }

    /// Perform one backward-Euler consolidation step of size `dt` \[s\].
    ///
    /// Solves the coupled system documented on [`ConsolidationSolver`] for
    /// `{u, p}_{n+1}` with external loads `f_u` (length `n_u`) and fluid
    /// sources `f_p` (length `n_p`), updates [`Self::u_dofs`],
    /// [`Self::p_dofs`] and advances [`Self::time`].  `dt = 0` gives the
    /// instantaneous (undrained) response to a load change.
    ///
    /// # Errors
    /// * [`Error::InvalidInput`] for a negative / non-finite `dt`, wrong load
    ///   vector lengths, an invalid tolerance or iteration limit, or a stiffness
    ///   block that is not positive on a free displacement DOF.
    /// * [`Error::SingularSystem`] when a CG solve detects a matrix that is not
    ///   positive definite on the free DOFs (e.g. unconstrained rigid-body
    ///   modes).
    /// * [`Error::NotConverged`] when a CG solve exhausts `max_iterations`.
    pub fn step(&mut self, dt: f64, f_u: &[f64], f_p: &[f64]) -> Result<ConsolidationStepReport> {
        if !dt.is_finite() || dt < 0.0 {
            return Err(Error::InvalidInput(format!(
                "time step must be finite and non-negative, got {dt}"
            )));
        }
        let n_u = self.n_u();
        let n_p = self.n_p();
        if f_u.len() != n_u {
            return Err(Error::InvalidInput(format!(
                "f_u has length {}, expected {n_u}",
                f_u.len()
            )));
        }
        if f_p.len() != n_p {
            return Err(Error::InvalidInput(format!(
                "f_p has length {}, expected {n_p}",
                f_p.len()
            )));
        }
        if !(self.tolerance.is_finite() && self.tolerance > 0.0) {
            return Err(Error::InvalidInput(format!(
                "tolerance must be finite and positive, got {}",
                self.tolerance
            )));
        }
        if self.max_iterations == 0 {
            return Err(Error::InvalidInput(
                "max_iterations must be at least 1".to_string(),
            ));
        }
        for block in [
            BiotBlock::Stiffness,
            BiotBlock::Coupling,
            BiotBlock::Storage,
            BiotBlock::Permeability,
        ] {
            self.flush_block(block);
        }
        if self.prepared.is_none() {
            self.prepared = Some(PreparedBiot::new(
                &self.stiffness,
                &self.coupling,
                &self.storage,
                &self.permeability,
                &self.u_constraints,
                &self.p_constraints,
            )?);
        }
        let prepared = self
            .prepared
            .as_ref()
            .ok_or_else(|| Error::General("Biot operators were not prepared".to_string()))?;

        // Prescribed parts u_c, p_c (zero on free DOFs).
        let u_c: Vec<f64> = self
            .u_constraints
            .iter()
            .map(|c| c.unwrap_or(0.0))
            .collect();
        let p_c: Vec<f64> = self
            .p_constraints
            .iter()
            .map(|c| c.unwrap_or(0.0))
            .collect();

        // Mass-balance history term r_p = Qᵀ u_n + S p_n + dt f_p.
        let qt_un = prepared.coupling_t.mul_vec(&self.u_dofs);
        let s_pn = self.storage.mul_vec(&self.p_dofs);
        // Lifting of the prescribed values.
        let k_uc = self.stiffness.mul_vec(&u_c);
        let q_pc = self.coupling.mul_vec(&p_c);
        let qt_uc = prepared.coupling_t.mul_vec(&u_c);
        let s_pc = self.storage.mul_vec(&p_c);
        let h_pc = self.permeability.mul_vec(&p_c);

        let b_u: Vec<f64> = (0..n_u)
            .map(|i| {
                if prepared.u_free[i] {
                    f_u[i] - k_uc[i] + q_pc[i]
                } else {
                    0.0
                }
            })
            .collect();
        let b_p: Vec<f64> = (0..n_p)
            .map(|i| {
                if prepared.p_free[i] {
                    qt_un[i] + s_pn[i] + dt * f_p[i] - qt_uc[i] - s_pc[i] - dt * h_pc[i]
                } else {
                    0.0
                }
            })
            .collect();

        let inner_tolerance = (self.tolerance * 1e-2).max(1e-14);
        let max_iterations = self.max_iterations;
        let mut inner_iterations = 0_usize;

        // w = K̃⁻¹ b_u ;  g = b_p − Q̃ᵀ w.
        let mut w = vec![0.0; n_u];
        inner_iterations +=
            prepared.solve_stiffness(&b_u, &mut w, inner_tolerance, max_iterations)?;
        let qt_w = prepared.coupling_masked_t.mul_vec(&w);
        let g: Vec<f64> = b_p.iter().zip(&qt_w).map(|(b, q)| b - q).collect();

        // Diagonal preconditioner of Σ̃ = S̃ + dt H̃ + Q̃ᵀ K̃⁻¹ Q̃ (+ I on fixed rows).
        let schur_diag_inv: Vec<f64> = (0..n_p)
            .map(|i| {
                let d = if prepared.p_free[i] {
                    prepared.storage_diag[i]
                        + dt * prepared.permeability_diag[i]
                        + prepared.coupling_schur_diag[i]
                } else {
                    1.0
                };
                if d.is_finite() && d > 0.0 {
                    1.0 / d
                } else {
                    1.0
                }
            })
            .collect();

        let mut p_free_part: Vec<f64> = (0..n_p)
            .map(|i| {
                if prepared.p_free[i] {
                    self.p_dofs[i]
                } else {
                    0.0
                }
            })
            .collect();
        let mut schur_work = vec![0.0; n_u];
        let outer = preconditioned_cg(
            |x, y| {
                let qx = prepared.coupling_masked.mul_vec(x);
                schur_work.iter_mut().for_each(|v| *v = 0.0);
                inner_iterations += prepared.solve_stiffness(
                    &qx,
                    &mut schur_work,
                    inner_tolerance,
                    max_iterations,
                )?;
                let qt_kq = prepared.coupling_masked_t.mul_vec(&schur_work);
                let sx = prepared.storage_masked.mul_vec(x);
                let hx = prepared.permeability_masked.mul_vec(x);
                for i in 0..y.len() {
                    y[i] = if prepared.p_free[i] {
                        qt_kq[i] + sx[i] + dt * hx[i]
                    } else {
                        x[i]
                    };
                }
                Ok(())
            },
            |r, z| {
                for ((zi, ri), di) in z.iter_mut().zip(r).zip(&schur_diag_inv) {
                    *zi = ri * di;
                }
            },
            &g,
            &mut p_free_part,
            self.tolerance,
            max_iterations,
            "Biot pressure Schur-complement CG",
        )?;

        // u_f = K̃⁻¹ (b_u + Q̃ p_f).
        let q_pf = prepared.coupling_masked.mul_vec(&p_free_part);
        let rhs_u: Vec<f64> = b_u.iter().zip(&q_pf).map(|(b, q)| b + q).collect();
        let mut u_free_part: Vec<f64> = (0..n_u)
            .map(|i| {
                if prepared.u_free[i] {
                    self.u_dofs[i]
                } else {
                    0.0
                }
            })
            .collect();
        inner_iterations +=
            prepared.solve_stiffness(&rhs_u, &mut u_free_part, inner_tolerance, max_iterations)?;

        for i in 0..n_u {
            self.u_dofs[i] = if prepared.u_free[i] {
                u_free_part[i]
            } else {
                u_c[i]
            };
        }
        for i in 0..n_p {
            self.p_dofs[i] = if prepared.p_free[i] {
                p_free_part[i]
            } else {
                p_c[i]
            };
        }
        self.time += dt;
        Ok(ConsolidationStepReport {
            outer_iterations: outer.iterations,
            inner_iterations,
            relative_residual: outer.relative_residual,
        })
    }
}

/// Biot operators with the Dirichlet DOFs eliminated, cached between steps.
#[derive(Debug, Clone)]
struct PreparedBiot {
    u_free: Vec<bool>,
    p_free: Vec<bool>,
    /// `M_u K M_u + (I − M_u)`.
    stiffness_masked: CsrMatrix,
    stiffness_preconditioner: SpdPreconditioner,
    /// Full `Qᵀ` (history and lifting terms).
    coupling_t: CsrMatrix,
    /// `M_u Q M_p`.
    coupling_masked: CsrMatrix,
    /// `(M_u Q M_p)ᵀ`.
    coupling_masked_t: CsrMatrix,
    /// `M_p S M_p`.
    storage_masked: CsrMatrix,
    /// `M_p H M_p`.
    permeability_masked: CsrMatrix,
    storage_diag: Vec<f64>,
    permeability_diag: Vec<f64>,
    /// `diag(Q̃ᵀ diag(K̃)⁻¹ Q̃)`.
    coupling_schur_diag: Vec<f64>,
}

impl PreparedBiot {
    fn new(
        stiffness: &CsrMatrix,
        coupling: &CsrMatrix,
        storage: &CsrMatrix,
        permeability: &CsrMatrix,
        u_constraints: &[Option<f64>],
        p_constraints: &[Option<f64>],
    ) -> Result<Self> {
        let u_free: Vec<bool> = u_constraints.iter().map(Option::is_none).collect();
        let p_free: Vec<bool> = p_constraints.iter().map(Option::is_none).collect();
        let stiffness_masked = mask_matrix(stiffness, &u_free, &u_free, true);
        let stiffness_diag = stiffness_masked.diagonal_vec();
        if let Some((dof, value)) = stiffness_diag
            .iter()
            .enumerate()
            .find(|&(_, &d)| !(d.is_finite() && d > 0.0))
        {
            return Err(Error::InvalidInput(format!(
                "stiffness block has non-positive diagonal {value} at free displacement DOF {dof}; \
                 assemble K_uu or constrain the DOF"
            )));
        }
        let stiffness_preconditioner = SpdPreconditioner::new(&stiffness_masked, &stiffness_diag);
        let coupling_masked = mask_matrix(coupling, &u_free, &p_free, false);
        let coupling_masked_t = coupling_masked.transpose();
        let mut coupling_schur_diag = vec![0.0; coupling.ncols];
        for (row, &k_rr) in stiffness_diag.iter().enumerate() {
            for idx in coupling_masked.row_ptr[row]..coupling_masked.row_ptr[row + 1] {
                let v = coupling_masked.values[idx];
                coupling_schur_diag[coupling_masked.col_indices[idx]] += v * v / k_rr;
            }
        }
        let storage_masked = mask_matrix(storage, &p_free, &p_free, false);
        let permeability_masked = mask_matrix(permeability, &p_free, &p_free, false);
        Ok(Self {
            storage_diag: storage_masked.diagonal_vec(),
            permeability_diag: permeability_masked.diagonal_vec(),
            u_free,
            p_free,
            stiffness_masked,
            stiffness_preconditioner,
            coupling_t: coupling.transpose(),
            coupling_masked,
            coupling_masked_t,
            storage_masked,
            permeability_masked,
            coupling_schur_diag,
        })
    }

    /// Solve `K̃ x = b` by preconditioned CG starting from the given `x`;
    /// returns the iteration count.
    fn solve_stiffness(
        &self,
        b: &[f64],
        x: &mut [f64],
        tolerance: f64,
        max_iterations: usize,
    ) -> Result<usize> {
        let outcome = preconditioned_cg(
            |v, y| {
                let kv = self.stiffness_masked.mul_vec(v);
                y.copy_from_slice(&kv);
                Ok(())
            },
            |r, z| self.stiffness_preconditioner.apply(r, z),
            b,
            x,
            tolerance,
            max_iterations,
            "Biot stiffness CG",
        )?;
        Ok(outcome.iterations)
    }
}

// ============================================================================
// TerzaghiConsolidation
// ============================================================================

/// Analytical 1-D Terzaghi consolidation model.
///
/// Describes the time-dependent settlement of a saturated clay layer of
/// thickness H drained at the top (single drainage) or both surfaces (double
/// drainage) under a suddenly applied load.
pub struct TerzaghiConsolidation {
    /// Drainage height H \[m\] (half-thickness for double drainage).
    pub drainage_height: f64,
    /// Coefficient of consolidation c_v \[m²/s\].
    pub cv: f64,
}

impl TerzaghiConsolidation {
    /// Create a new `TerzaghiConsolidation`.
    ///
    /// # Arguments
    /// * `drainage_height` – drainage path length H \[m\]
    /// * `cv`              – coefficient of consolidation c_v \[m²/s\]
    pub fn new(drainage_height: f64, cv: f64) -> Self {
        Self {
            drainage_height,
            cv,
        }
    }

    /// Dimensionless time factor T_v = c_v t / H².
    pub fn time_factor(&self, time: f64) -> f64 {
        self.cv * time / (self.drainage_height * self.drainage_height)
    }

    /// Average degree of consolidation U (0–1) using Terzaghi series.
    ///
    /// Uses `terzaghi_u` with T_v = c_v t / H².
    pub fn degree_of_consolidation(&self, time: f64) -> f64 {
        let tv = self.time_factor(time);
        terzaghi_u(tv)
    }

    /// Settlement at time `time` given final settlement `s_final` \[m\].
    pub fn settlement(&self, time: f64, s_final: f64) -> f64 {
        self.degree_of_consolidation(time) * s_final
    }

    /// Time \[s\] required to reach degree of consolidation `u_target` (0–1).
    ///
    /// Uses the classical Terzaghi/Casagrande closed-form fits:
    ///
    /// * `T_v = (π/4)·U²` for `U ≤ 0.6`,
    /// * `T_v = 1.781 − 0.933·log10(100 − U%) = −0.9332·log10(1 − U) − 0.0851`
    ///   for `U > 0.6`.
    ///
    /// Both branches agree to within 0.004 at `U = 0.6`.  A non-positive (or
    /// NaN) target returns `0.0`; a target `≥ 1` is never reached in finite
    /// time and returns `f64::INFINITY`.
    pub fn time_to_consolidation(&self, u_target: f64) -> f64 {
        if u_target.is_nan() || u_target <= 0.0 {
            return 0.0;
        }
        if u_target >= 1.0 {
            return f64::INFINITY;
        }
        let tv = if u_target <= 0.6 {
            PI / 4.0 * u_target * u_target
        } else {
            -0.9332 * (1.0 - u_target).log10() - 0.0851
        };
        tv * self.drainage_height * self.drainage_height / self.cv
    }
}

// ============================================================================
// DarcyFlow
// ============================================================================

/// Darcy seepage flow model.
///
/// Computes the specific discharge (Darcy velocity) q = -k/μ * ∇h
/// for laminar flow through a porous medium.
pub struct DarcyFlow {
    /// Intrinsic permeability k \[m²\].
    pub permeability: f64,
    /// Dynamic viscosity μ \[Pa·s\].
    pub viscosity: f64,
    /// Hydraulic gradient ∇h \[m/m\] (dimensionless).
    pub hydraulic_gradient: [f64; 3],
}

impl DarcyFlow {
    /// Create a new `DarcyFlow`.
    pub fn new(permeability: f64, viscosity: f64, hydraulic_gradient: [f64; 3]) -> Self {
        Self {
            permeability,
            viscosity,
            hydraulic_gradient,
        }
    }

    /// Compute the Darcy velocity q = -k/μ * ∇h \[m/s\].
    pub fn darcy_velocity(&self) -> [f64; 3] {
        let factor = -self.permeability / self.viscosity;
        [
            factor * self.hydraulic_gradient[0],
            factor * self.hydraulic_gradient[1],
            factor * self.hydraulic_gradient[2],
        ]
    }

    /// Magnitude of the Darcy velocity |q| \[m/s\].
    pub fn darcy_speed(&self) -> f64 {
        let q = self.darcy_velocity();
        (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt()
    }

    /// Hydraulic conductivity K = k ρ g / μ \[m/s\].
    ///
    /// # Arguments
    /// * `density` – fluid density ρ \[kg/m³\]
    /// * `g`       – gravitational acceleration \[m/s²\]
    pub fn hydraulic_conductivity(&self, density: f64, g: f64) -> f64 {
        self.permeability * density * g / self.viscosity
    }

    /// Reynolds number for Darcy flow Re = ρ |q| d / μ.
    ///
    /// # Arguments
    /// * `density`          – fluid density ρ \[kg/m³\]
    /// * `particle_diameter`– mean grain diameter d \[m\]
    pub fn reynolds_number(&self, density: f64, particle_diameter: f64) -> f64 {
        density * self.darcy_speed() * particle_diameter / self.viscosity
    }
}

// ============================================================================
// Free functions
// ============================================================================

/// Average degree of consolidation U using Terzaghi Fourier series.
///
/// U(T_v) = 1 - Σ_{m=0}^{∞} (8 / (2m+1)² π²) exp(-(2m+1)² π² T_v / 4)
///
/// The series is truncated after 50 terms for practical accuracy.
pub fn terzaghi_u(tv: f64) -> f64 {
    if tv <= 0.0 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for m in 0_u32..50 {
        let two_m1 = (2 * m + 1) as f64;
        let coeff = 8.0 / (two_m1 * two_m1 * PI * PI);
        let exponent = -two_m1 * two_m1 * PI * PI * tv / 4.0;
        sum += coeff * exponent.exp();
    }
    (1.0 - sum).clamp(0.0, 1.0)
}

/// Real time t \[s\] from dimensionless time factor T_v, c_v \[m²/s\] and H \[m\].
///
/// t = T_v * H² / c_v
pub fn consolidation_time(tv: f64, cv: f64, h: f64) -> f64 {
    tv * h * h / cv
}

/// Effective stress σ′ = σ - p (Terzaghi effective stress principle) \[Pa\].
///
/// # Arguments
/// * `total`         – total stress σ \[Pa\]
/// * `pore_pressure` – pore fluid pressure p \[Pa\]
pub fn effective_stress(total: f64, pore_pressure: f64) -> f64 {
    total - pore_pressure
}

/// Consolidation coefficient c_v from permeability, compressibility and fluid properties.
///
/// c_v = k / (μ m_v)  where m_v is the coefficient of volume compressibility.
///
/// # Arguments
/// * `k`   – hydraulic conductivity k \[m/s\] (K = k / (ρ g) style input in m/s)
/// * `mu`  – dynamic viscosity \[Pa·s\]
/// * `m_v` – coefficient of volume compressibility \[1/Pa\]
pub fn consolidation_coefficient(k: f64, mu: f64, m_v: f64) -> f64 {
    if mu.abs() < 1e-300 || m_v.abs() < 1e-300 {
        return 0.0;
    }
    k / (mu * m_v)
}

/// Pore pressure dissipation ratio at depth z in a layer of height H at time factor T_v.
///
/// u(z, T_v) / u_0 = Σ_{m=0}^{∞} (4 / (2m+1) π) sin((2m+1) π z / (2H)) exp(...)
pub fn pore_pressure_ratio(z: f64, h: f64, tv: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for m in 0_u32..50 {
        let two_m1 = (2 * m + 1) as f64;
        let coeff = 4.0 / (two_m1 * PI);
        let spatial = (two_m1 * PI * z / (2.0 * h)).sin();
        let temporal = (-two_m1 * two_m1 * PI * PI * tv / 4.0).exp();
        sum += coeff * spatial * temporal;
    }
    sum.clamp(0.0, 1.0)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- BiotCoeff ---

    #[test]
    fn test_biot_new() {
        let b = BiotCoeff::new(0.8, 1e9);
        assert!((b.alpha - 0.8).abs() < 1e-12);
        assert!((b.modulus - 1e9).abs() < 1.0);
    }

    #[test]
    fn test_biot_compute_from_bulk_incompressible_solid() {
        // If K_dry ≈ 0 (very soft skeleton) and K_s >> K_f, α → 1
        let b = BiotCoeff::compute_from_bulk(0.0, 1e12, 2.2e9, 0.4);
        assert!((b.alpha - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_biot_compute_from_bulk_typical_sandstone() {
        // Typical sandstone: K_dry=10 GPa, K_s=36 GPa, K_f=2.2 GPa, phi=0.2
        let b = BiotCoeff::compute_from_bulk(10e9, 36e9, 2.2e9, 0.2);
        assert!(b.alpha > 0.0 && b.alpha < 1.0);
        assert!(b.modulus > 0.0);
    }

    #[test]
    fn test_biot_alpha_range() {
        let b = BiotCoeff::compute_from_bulk(5e9, 36e9, 2.2e9, 0.3);
        assert!((0.0..=1.0).contains(&b.alpha));
    }

    #[test]
    fn test_biot_storage_coefficient() {
        let b = BiotCoeff::new(1.0, 1e9);
        let s = b.storage_coefficient(0.4, 2.2e9);
        assert!(s > 0.0);
    }

    #[test]
    fn test_biot_undrained_bulk() {
        let k_u = BiotCoeff::undrained_bulk(10e9, 0.8, 5e9);
        // K_u > K_dry
        assert!(k_u > 10e9);
    }

    #[test]
    fn test_biot_skempton_b() {
        let b_val = BiotCoeff::skempton_b(10e9, 0.8, 5e9);
        assert!(b_val > 0.0 && b_val <= 1.0);
    }

    // --- PorousElement ---

    #[test]
    fn test_porous_element_new() {
        let k = [1e-12, 0.0, 0.0, 0.0, 1e-12, 0.0, 0.0, 0.0, 1e-12];
        let elem = PorousElement::new(vec![0, 1, 2, 3], k, 0.3, 1e7);
        assert_eq!(elem.node_ids.len(), 4);
        assert!((elem.porosity - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_porous_element_coupling_matrix_length() {
        let k = [1e-12; 9];
        let elem = PorousElement::new(vec![0, 1, 2, 3], k, 0.3, 1e7);
        let coupling = elem.compute_coupling_matrix();
        assert_eq!(coupling.len(), 12); // 3 * 4
    }

    #[test]
    fn test_porous_element_coupling_matrix_values() {
        let k = [1e-12; 9];
        let elem = PorousElement::new(vec![0, 1, 2, 3], k, 0.3, 1e7);
        let coupling = elem.compute_coupling_matrix();
        // Each entry should be 1/4 = 0.25
        for &v in &coupling {
            assert!((v - 0.25).abs() < 1e-12);
        }
    }

    #[test]
    fn test_porous_element_k_xx() {
        let mut k = [0.0_f64; 9];
        k[0] = 2e-12;
        let elem = PorousElement::new(vec![0], k, 0.2, 1e8);
        assert!((elem.k_xx() - 2e-12).abs() < 1e-300);
    }

    #[test]
    fn test_porous_element_isotropic_permeability() {
        let k = [1e-12, 0.0, 0.0, 0.0, 2e-12, 0.0, 0.0, 0.0, 3e-12];
        let elem = PorousElement::new(vec![0], k, 0.2, 1e8);
        assert!((elem.isotropic_permeability() - 2e-12).abs() < 1e-24);
    }

    // --- ConsolidationSolver ---

    #[test]
    fn test_consolidation_solver_new() {
        let solver = ConsolidationSolver::new(6, 2);
        assert_eq!(solver.n_u(), 6);
        assert_eq!(solver.n_p(), 2);
        assert!((solver.time - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_consolidation_solver_is_sparse_for_large_systems() {
        // The former dense layout allocated n_u² + n_u·n_p doubles (≈ 240 GB here).
        let mut solver = ConsolidationSolver::new(150_000, 50_000);
        for block in [
            BiotBlock::Stiffness,
            BiotBlock::Coupling,
            BiotBlock::Storage,
            BiotBlock::Permeability,
        ] {
            assert_eq!(solver.block(block).nnz(), 0);
        }
        assert_eq!(solver.block_shape(BiotBlock::Coupling), (150_000, 50_000));
    }

    #[test]
    fn test_consolidation_solver_element_assembly_scales() {
        // 20 000 element-by-element insertions must stay cheap and sparse.
        let n_el = 20_000;
        let mut solver = biot_column(n_el, 1.0, 1.0, 1.0, 1.0);
        let n = n_el + 1;
        assert_eq!(solver.block(BiotBlock::Stiffness).nnz(), 3 * n - 2);
        assert_eq!(solver.block(BiotBlock::Permeability).nnz(), 3 * n - 2);
        // Q has the tridiagonal pattern too; its interior diagonal cancels to 0.
        assert_eq!(solver.block(BiotBlock::Coupling).nnz(), 3 * n - 2);
        assert_eq!(solver.block(BiotBlock::Coupling).get(7, 7), 0.0);
        assert_eq!(solver.block(BiotBlock::Storage).nnz(), 0);
        // Interior stiffness row: (−M/h, 2M/h, −M/h) with h = 1/n_el.
        let k = solver.block(BiotBlock::Stiffness);
        let scale = n_el as f64;
        assert!((k.get(7, 7) - 2.0 * scale).abs() < 1e-9 * scale);
        assert!((k.get(7, 6) + scale).abs() < 1e-9 * scale);
    }

    /// Equal-order linear 1-D Biot column on `[0, height]` with `n_el` elements:
    /// base (node 0) fixed, top (node `n_el`) drained, incompressible
    /// constituents (`S = 0`).
    fn biot_column(
        n_el: usize,
        height: f64,
        modulus: f64,
        mobility: f64,
        alpha: f64,
    ) -> ConsolidationSolver {
        let n = n_el + 1;
        let h = height / n_el as f64;
        let mut solver = ConsolidationSolver::new(n, n);
        for e in 0..n_el {
            let dofs = [e, e + 1];
            let ke = [1.0, -1.0, -1.0, 1.0].map(|v| v * modulus / h);
            let qe = [-1.0, -1.0, 1.0, 1.0].map(|v| v * alpha / 2.0);
            let he = [1.0, -1.0, -1.0, 1.0].map(|v| v * mobility / h);
            solver
                .add_element_block(BiotBlock::Stiffness, &dofs, &dofs, &ke)
                .expect("stiffness");
            solver
                .add_element_block(BiotBlock::Coupling, &dofs, &dofs, &qe)
                .expect("coupling");
            solver
                .add_element_block(BiotBlock::Permeability, &dofs, &dofs, &he)
                .expect("permeability");
        }
        solver.fix_displacement(0, 0.0).expect("fix base");
        solver.fix_pressure(n_el, 0.0).expect("drain top");
        solver
    }

    #[test]
    fn test_consolidation_solver_terzaghi_column_matches_analytical() {
        let (height, modulus, mobility, load) = (2.0, 5.0e6, 2.0e-9, 1.0e4);
        let n_el = 20;
        let mut solver = biot_column(n_el, height, modulus, mobility, 1.0);
        let mut f_u = vec![0.0; n_el + 1];
        f_u[n_el] = -load;
        let f_p = vec![0.0; n_el + 1];
        // α = 1, S = 0 → c_v = (k/μ)·M; single drainage → H_dr = height.
        let cv = mobility * modulus;
        let tv_end = 0.2;
        let t_end = tv_end * height * height / cv;
        let n_steps = 400;
        let dt = t_end / n_steps as f64;
        let mut outer_total = 0;
        for _ in 0..n_steps {
            let report = solver.step(dt, &f_u, &f_p).expect("Biot step");
            assert!(report.relative_residual <= solver.tolerance);
            outer_total += report.outer_iterations;
        }
        assert!(outer_total > 0);
        assert!((solver.time - t_end).abs() < 1e-9 * t_end);

        let s_final = load * height / modulus;
        let degree = -solver.u_dofs[n_el] / s_final;
        let exact = terzaghi_u(tv_end);
        assert!(
            ((degree - exact) / exact).abs() < 0.03,
            "U_num = {degree}, U_exact = {exact}"
        );
        // Excess pore pressure at the impermeable base (distance H from the drain).
        let p_base = solver.p_dofs[0] / load;
        let p_exact = pore_pressure_ratio(height, height, tv_end);
        assert!(
            (p_base - p_exact).abs() < 0.03,
            "p_base/q = {p_base}, exact = {p_exact}"
        );
        assert_eq!(solver.u_dofs[0], 0.0);
        assert_eq!(solver.p_dofs[n_el], 0.0);
    }

    #[test]
    fn test_consolidation_solver_column_reaches_drained_state() {
        let (height, modulus, mobility, load) = (1.0, 1.0, 1.0, 2.0);
        let n_el = 10;
        let mut solver = biot_column(n_el, height, modulus, mobility, 1.0);
        let mut f_u = vec![0.0; n_el + 1];
        f_u[n_el] = -load;
        let f_p = vec![0.0; n_el + 1];
        for _ in 0..60 {
            solver.step(0.25, &f_u, &f_p).expect("Biot step");
        }
        // Fully drained: p → 0, linear settlement profile u(z) = −q z / M.
        for (i, (&u, &p)) in solver.u_dofs.iter().zip(&solver.p_dofs).enumerate() {
            let z = height * i as f64 / n_el as f64;
            assert!((u + load * z / modulus).abs() < 1e-6, "u[{i}] = {u}");
            assert!(p.abs() < 1e-6, "p[{i}] = {p}");
        }
    }

    #[test]
    fn test_consolidation_solver_undrained_then_drained_single_dof() {
        // K = 4, Q = 2, S = 0.5, H = 1.
        let mut solver = ConsolidationSolver::new(1, 1);
        solver
            .add_triplets(BiotBlock::Stiffness, &[(0, 0, 4.0)])
            .expect("K");
        solver
            .add_triplets(BiotBlock::Coupling, &[(0, 0, 2.0)])
            .expect("Q");
        solver
            .add_triplets(BiotBlock::Storage, &[(0, 0, 0.5)])
            .expect("S");
        solver
            .add_triplets(BiotBlock::Permeability, &[(0, 0, 1.0)])
            .expect("H");
        // Undrained (dt = 0): 4u − 2p = 12, 2u + 0.5p = 0 → u = 1, p = −4.
        solver.step(0.0, &[12.0], &[0.0]).expect("undrained");
        assert!(
            (solver.u_dofs[0] - 1.0).abs() < 1e-9,
            "u = {}",
            solver.u_dofs[0]
        );
        assert!(
            (solver.p_dofs[0] + 4.0).abs() < 1e-9,
            "p = {}",
            solver.p_dofs[0]
        );
        assert_eq!(solver.time, 0.0);
        // One backward-Euler step, dt = 1, solved by hand:
        // 4u − 2p = 12 ; 2u + 1.5p = 2·1 + 0.5·(−4) = 0 → u = 1.8, p = −2.4.
        solver.step(1.0, &[12.0], &[0.0]).expect("step");
        assert!(
            (solver.u_dofs[0] - 1.8).abs() < 1e-9,
            "u = {}",
            solver.u_dofs[0]
        );
        assert!(
            (solver.p_dofs[0] + 2.4).abs() < 1e-9,
            "p = {}",
            solver.p_dofs[0]
        );
        // Long-term drained limit: p → 0, u → f/K = 3.
        for _ in 0..200 {
            solver.step(1.0, &[12.0], &[0.0]).expect("step");
        }
        assert!((solver.u_dofs[0] - 3.0).abs() < 1e-8);
        assert!(solver.p_dofs[0].abs() < 1e-8);
        assert!((solver.time - 201.0).abs() < 1e-9);
    }

    #[test]
    fn test_consolidation_solver_pure_pressure_diffusion() {
        // No displacement DOFs: (S + dt H) p_{n+1} = S p_n.
        let mut solver = ConsolidationSolver::new(0, 1);
        solver.p_dofs[0] = 100.0;
        solver
            .add_triplets(BiotBlock::Storage, &[(0, 0, 1.0)])
            .expect("S");
        solver
            .add_triplets(BiotBlock::Permeability, &[(0, 0, 1.0)])
            .expect("H");
        solver.step(1.0, &[], &[0.0]).expect("step");
        assert!((solver.p_dofs[0] - 50.0).abs() < 1e-8);
        assert!((solver.time - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_consolidation_solver_prescribed_pressures_give_linear_profile() {
        // Steady seepage between p = 0 and p = 1 through a uniform 1-D chain.
        let n = 5;
        let mut solver = ConsolidationSolver::new(0, n);
        let mut h = Vec::new();
        for e in 0..n - 1 {
            h.extend_from_slice(&[
                (e, e, 1.0),
                (e, e + 1, -1.0),
                (e + 1, e, -1.0),
                (e + 1, e + 1, 1.0),
            ]);
        }
        solver.add_triplets(BiotBlock::Permeability, &h).expect("H");
        solver.fix_pressure(0, 0.0).expect("p0");
        solver.fix_pressure(n - 1, 1.0).expect("p1");
        solver.step(1.0, &[], &vec![0.0; n]).expect("step");
        for (i, &p) in solver.p_dofs.iter().enumerate() {
            let expected = i as f64 / (n - 1) as f64;
            assert!((p - expected).abs() < 1e-8, "p[{i}] = {p}");
        }
    }

    #[test]
    fn test_consolidation_solver_elastic_only_with_prescribed_displacement() {
        // Spring chain (k = 2) without pressure DOFs.
        let mut solver = ConsolidationSolver::new(3, 0);
        let k = [
            (0, 0, 2.0),
            (0, 1, -2.0),
            (1, 0, -2.0),
            (1, 1, 4.0),
            (1, 2, -2.0),
            (2, 1, -2.0),
            (2, 2, 2.0),
        ];
        solver.add_triplets(BiotBlock::Stiffness, &k).expect("K");
        solver.fix_displacement(0, 0.0).expect("fix");
        solver.step(0.0, &[0.0, 0.0, 1.0], &[]).expect("load");
        assert!((solver.u_dofs[1] - 0.5).abs() < 1e-9);
        assert!((solver.u_dofs[2] - 1.0).abs() < 1e-9);
        // Prescribe the free end instead of loading it.
        solver.fix_displacement(2, 0.3).expect("fix end");
        solver.step(0.0, &[0.0; 3], &[]).expect("prescribed");
        assert!((solver.u_dofs[1] - 0.15).abs() < 1e-9);
        assert_eq!(solver.u_dofs[2], 0.3);
    }

    #[test]
    fn test_consolidation_solver_set_block_roundtrip_and_shape_check() {
        let mut solver = ConsolidationSolver::new(2, 1);
        let q = CsrMatrix::from_triplets(2, 1, &[(0, 0, 0.5), (1, 0, -0.5)]);
        solver.set_block(BiotBlock::Coupling, q).expect("Q");
        assert_eq!(solver.block(BiotBlock::Coupling).get(1, 0), -0.5);
        let wrong = CsrMatrix::from_triplets(1, 2, &[(0, 1, 1.0)]);
        assert!(matches!(
            solver.set_block(BiotBlock::Coupling, wrong),
            Err(Error::InvalidInput(_))
        ));
        solver
            .add_triplets(BiotBlock::Coupling, &[(0, 0, 0.25)])
            .expect("add");
        assert!((solver.block(BiotBlock::Coupling).get(0, 0) - 0.75).abs() < 1e-15);
    }

    #[test]
    fn test_consolidation_solver_input_errors() {
        let mut solver = ConsolidationSolver::new(2, 2);
        assert!(matches!(
            solver.add_triplets(BiotBlock::Coupling, &[(2, 0, 1.0)]),
            Err(Error::DofOutOfRange { dof: 2, n: 2 })
        ));
        assert!(matches!(
            solver.add_triplets(BiotBlock::Storage, &[(0, 0, f64::NAN)]),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            solver.add_element_block(BiotBlock::Stiffness, &[0, 1], &[0, 1], &[1.0; 3]),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            solver.fix_pressure(5, 0.0),
            Err(Error::DofOutOfRange { dof: 5, n: 2 })
        ));
        assert!(matches!(
            solver.step(-1.0, &[0.0; 2], &[0.0; 2]),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            solver.step(1.0, &[0.0; 3], &[0.0; 2]),
            Err(Error::InvalidInput(_))
        ));
        // Free displacement DOFs without stiffness are rejected, state untouched.
        assert!(matches!(
            solver.step(1.0, &[1.0, 0.0], &[0.0; 2]),
            Err(Error::InvalidInput(_))
        ));
        assert_eq!(solver.time, 0.0);
    }

    #[test]
    fn test_not_converged_error_maps_to_fem_error() {
        let err = Error::NotConverged {
            solver: "Biot stiffness CG".to_string(),
            iterations: 7,
            residual: 0.5,
        };
        assert!(err.to_string().contains("7 iterations"));
        let fem_err: crate::analysis::FemError = err.into();
        assert!(matches!(
            fem_err,
            crate::analysis::FemError::SolverNotConverged
        ));
    }

    #[test]
    fn test_consolidation_solver_floating_body_is_reported() {
        // Unconstrained spring: K is singular and the load is not balanced.
        let mut solver = ConsolidationSolver::new(2, 0);
        solver
            .add_triplets(
                BiotBlock::Stiffness,
                &[(0, 0, 1.0), (0, 1, -1.0), (1, 0, -1.0), (1, 1, 1.0)],
            )
            .expect("K");
        solver.max_iterations = 50;
        let result = solver.step(0.0, &[1.0, 0.0], &[]);
        assert!(
            matches!(
                result,
                Err(Error::SingularSystem | Error::NotConverged { .. })
            ),
            "{result:?}"
        );
        // A failed step leaves the state untouched.
        assert_eq!(solver.time, 0.0);
        assert_eq!(solver.u_dofs, vec![0.0, 0.0]);
        assert!(solver.p_dofs.is_empty());
    }

    // --- TerzaghiConsolidation ---

    #[test]
    fn test_terzaghi_time_factor_zero() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        assert!((tc.time_factor(0.0)).abs() < 1e-12);
    }

    #[test]
    fn test_terzaghi_time_factor() {
        let cv = 1e-6;
        let h = 2.0;
        let tc = TerzaghiConsolidation::new(h, cv);
        let t = 1e5;
        let tv = cv * t / (h * h);
        assert!((tc.time_factor(t) - tv).abs() < 1e-20);
    }

    #[test]
    fn test_terzaghi_degree_zero_time() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        assert!((tc.degree_of_consolidation(0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_terzaghi_degree_large_time() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        // t → ∞ should give U → 1
        let u = tc.degree_of_consolidation(1e20);
        assert!((u - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_terzaghi_degree_monotonic() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-7);
        let u1 = tc.degree_of_consolidation(1e4);
        let u2 = tc.degree_of_consolidation(2e4);
        assert!(u2 >= u1);
    }

    #[test]
    fn test_terzaghi_settlement() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        let s = tc.settlement(1e20, 0.1);
        assert!((s - 0.1).abs() < 1e-5);
    }

    #[test]
    fn test_terzaghi_time_to_consolidation_50pct() {
        let cv = 1e-6;
        let h = 1.0;
        let tc = TerzaghiConsolidation::new(h, cv);
        let t50 = tc.time_to_consolidation(0.5);
        // Verify degree at computed time is near 0.5
        let u_check = tc.degree_of_consolidation(t50);
        assert!((u_check - 0.5).abs() < 0.05);
    }

    #[test]
    fn test_terzaghi_time_to_consolidation_90pct_uses_log10() {
        // Casagrande: T_90 = −0.9332·log10(0.1) − 0.0851 = 0.8481.
        let tc = TerzaghiConsolidation::new(1.0, 1.0);
        let t90 = tc.time_to_consolidation(0.9);
        assert!((t90 - 0.848).abs() < 0.002, "T_90 = {t90}");
    }

    #[test]
    fn test_terzaghi_time_to_consolidation_continuous_at_60pct() {
        let tc = TerzaghiConsolidation::new(1.0, 1.0);
        let below = tc.time_to_consolidation(0.6);
        let above = tc.time_to_consolidation(0.6 + 1e-9);
        assert!((above - below).abs() < 0.01, "jump {below} -> {above}");
    }

    #[test]
    fn test_terzaghi_time_to_consolidation_roundtrip() {
        let tc = TerzaghiConsolidation::new(2.5, 3e-7);
        for &u in &[0.3, 0.7, 0.95] {
            let t = tc.time_to_consolidation(u);
            let u_back = terzaghi_u(tc.time_factor(t));
            assert!(
                (u_back - u).abs() < 0.01,
                "U = {u}: round trip gave {u_back}"
            );
        }
    }

    #[test]
    fn test_terzaghi_time_to_consolidation_guards() {
        let tc = TerzaghiConsolidation::new(1.0, 1.0);
        assert_eq!(tc.time_to_consolidation(0.0), 0.0);
        assert_eq!(tc.time_to_consolidation(-0.5), 0.0);
        assert_eq!(tc.time_to_consolidation(f64::NAN), 0.0);
        assert!(tc.time_to_consolidation(1.0).is_infinite());
        assert!(tc.time_to_consolidation(1.5).is_infinite());
    }

    // --- DarcyFlow ---

    #[test]
    fn test_darcy_velocity_1d() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [1.0, 0.0, 0.0]);
        let v = flow.darcy_velocity();
        // q = -k/mu * grad = -1e-12 / 1e-3 = -1e-9 m/s
        assert!((v[0] - (-1e-9)).abs() < 1e-21);
        assert!((v[1]).abs() < 1e-30);
        assert!((v[2]).abs() < 1e-30);
    }

    #[test]
    fn test_darcy_velocity_zero_gradient() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [0.0, 0.0, 0.0]);
        let v = flow.darcy_velocity();
        assert!(v[0].abs() < 1e-30 && v[1].abs() < 1e-30 && v[2].abs() < 1e-30);
    }

    #[test]
    fn test_darcy_speed() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [1.0, 0.0, 0.0]);
        let speed = flow.darcy_speed();
        assert!((speed - 1e-9).abs() < 1e-21);
    }

    #[test]
    fn test_darcy_hydraulic_conductivity() {
        // k=1e-12 m², rho=1000, g=9.81, mu=1e-3
        let flow = DarcyFlow::new(1e-12, 1e-3, [0.0; 3]);
        let k_cond = flow.hydraulic_conductivity(1000.0, 9.81);
        let expected = 1e-12 * 1000.0 * 9.81 / 1e-3;
        assert!((k_cond - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_darcy_reynolds_number() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [1.0, 0.0, 0.0]);
        let re = flow.reynolds_number(1000.0, 1e-3);
        assert!(re >= 0.0);
    }

    // --- Free functions ---

    #[test]
    fn test_terzaghi_u_zero() {
        assert!((terzaghi_u(0.0)).abs() < 1e-12);
    }

    #[test]
    fn test_terzaghi_u_large() {
        assert!((terzaghi_u(100.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_terzaghi_u_half() {
        // At T_v ≈ 0.197, U ≈ 0.5
        let u = terzaghi_u(0.197);
        assert!((u - 0.5).abs() < 0.02);
    }

    #[test]
    fn test_terzaghi_u_monotonic() {
        let u1 = terzaghi_u(0.1);
        let u2 = terzaghi_u(0.3);
        let u3 = terzaghi_u(1.0);
        assert!(u1 <= u2 && u2 <= u3);
    }

    #[test]
    fn test_consolidation_time_roundtrip() {
        let cv = 2e-8;
        let h = 3.0;
        let tv = 0.5;
        let t = consolidation_time(tv, cv, h);
        let tv_back = cv * t / (h * h);
        assert!((tv_back - tv).abs() < 1e-12);
    }

    #[test]
    fn test_effective_stress_positive() {
        let sigma_eff = effective_stress(100e3, 30e3);
        assert!((sigma_eff - 70e3).abs() < 1e-8);
    }

    #[test]
    fn test_effective_stress_zero_pore() {
        assert!((effective_stress(50e3, 0.0) - 50e3).abs() < 1e-10);
    }

    #[test]
    fn test_effective_stress_negative() {
        // Tensile effective stress
        let s = effective_stress(10e3, 50e3);
        assert!((s - (-40e3)).abs() < 1e-8);
    }

    #[test]
    fn test_consolidation_coefficient_basic() {
        let cc = consolidation_coefficient(1e-9, 1e-3, 1e-7);
        assert!(cc > 0.0);
    }

    #[test]
    fn test_consolidation_coefficient_zero_viscosity() {
        let cc = consolidation_coefficient(1e-9, 0.0, 1e-7);
        assert!((cc).abs() < 1e-12);
    }

    #[test]
    fn test_pore_pressure_ratio_surface() {
        // At z=0 (drainage surface) ratio → 0
        let r = pore_pressure_ratio(0.0, 1.0, 1.0);
        assert!(r.abs() < 0.01);
    }

    #[test]
    fn test_pore_pressure_ratio_zero_time() {
        // At T_v=0, ratio should be ~1 at mid-layer
        let r = pore_pressure_ratio(1.0, 1.0, 0.0);
        assert!(r > 0.9);
    }

    #[test]
    fn test_pore_pressure_ratio_range() {
        for i in 0..10 {
            let tv = i as f64 * 0.1;
            let r = pore_pressure_ratio(0.5, 1.0, tv);
            assert!((0.0..=1.0).contains(&r));
        }
    }
}
