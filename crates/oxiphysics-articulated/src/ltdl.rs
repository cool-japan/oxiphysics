// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sparse `LᵀDL` factorization of the joint-space mass matrix `H`.
//!
//! The Composite Rigid Body Algorithm (CRBA) produces the symmetric
//! positive-definite joint-space inertia matrix `H(q)`. For a kinematic *tree*
//! this matrix is **branch-induced sparse**: entry `H[i][j]` is structurally
//! nonzero only when DOF `i` and DOF `j` lie on a common root-to-leaf path,
//! i.e. one is an ancestor of the other in the assembly tree. Featherstone's
//! `LᵀDL` factorization exploits exactly this structure: it touches only the
//! ancestor entries of each column, so the unit-lower-triangular factor `L`
//! inherits the sparsity pattern of `H` with **zero fill-in**.
//!
//! Concretely we factor `H = LᵀDL`, where `L` is unit lower-triangular (its
//! implicit diagonal is 1) and `D` is diagonal. A solve `H·x = b` then costs
//! `O(N·d)` operations (`N` = number of DOF, `d` = tree depth) instead of the
//! `O(N³)` of a dense Cholesky, which is a large win for branched mechanisms
//! such as humanoids and quadrupeds where `d ≪ N`.
//!
//! The sparsity bookkeeping is driven by the *expanded parent array* `λ`
//! (lambda): `λ[d]` is the DOF index of the parent of DOF `d` in the tree
//! (1-based), or `0` if `d` is a root DOF. Walking the `λ` chain from any DOF
//! enumerates its ancestors cheaply.
//!
//! Reference: Roy Featherstone, *Rigid Body Dynamics Algorithms*, Springer
//! 2008, chapters 6 and 9 (in particular Algorithm 9.5, the `LᵀDL`
//! factorization of a branched joint-space inertia matrix).

use crate::model::ArticulatedModel;

/// Build the expanded parent array `λ` (1-based) for an articulated model.
///
/// `λ[d]` is the DOF index of the parent of DOF `d` in the assembly tree
/// (1-based), or `0` if `d` is a root DOF. The valid range of `d` is
/// `1..=n` with `n = model.total_dof()`; index `0` is unused padding so that
/// the 1-based indexing used by the factorization lines up directly.
///
/// Within a multi-DOF joint the DOFs form a chain, so every DOF after the
/// first points at its immediate predecessor. The first DOF of a body points
/// at the *last* DOF of its nearest ancestor that actually has DOFs: fixed
/// (0-DOF) joints are transparent and are skipped by walking `parent_id`
/// upward until a body with at least one DOF is found (or the root is reached,
/// in which case `λ[d] = 0`).
///
/// The returned array always satisfies `λ[d] < d` for every non-root DOF `d`,
/// which is the invariant the `LᵀDL` recursion relies on.
pub fn build_lambda(model: &ArticulatedModel) -> Vec<usize> {
    let n = model.total_dof();
    let mut lambda = vec![0usize; n + 1];
    for b in 0..model.num_bodies() {
        let ds = model.dof_start(b);
        let dc = model.dof_count(b);
        for local_d in 0..dc {
            let d = ds + local_d + 1; // 1-based DOF index
            if local_d > 0 {
                // Chain within a multi-DOF joint: point at the predecessor.
                lambda[d] = d - 1;
            } else {
                // First DOF of this body: find the nearest ancestor that has
                // at least one DOF, skipping fixed (0-DOF) bodies.
                let mut ancestor = model.bodies[b].parent_id;
                let mut resolved = 0usize;
                while let Some(a) = ancestor {
                    if model.dof_count(a) > 0 {
                        // 1-based index of the last DOF of ancestor `a`.
                        resolved = model.dof_start(a) + model.dof_count(a);
                        break;
                    }
                    ancestor = model.bodies[a].parent_id;
                }
                lambda[d] = resolved;
            }
        }
    }
    lambda
}

/// Factor a full symmetric SPD matrix `H` in place as `H = LᵀDL`.
///
/// On entry `h` is the full symmetric SPD `n×n` CRBA matrix (both triangles
/// filled). On exit the storage is overwritten so that `h[k-1][k-1] = D[k]`
/// (the diagonal factor) and `h[k-1][i-1] = L[k][i]` for each ancestor `i` of
/// DOF `k`. `L` is unit lower-triangular with an *implicit* unit diagonal, and
/// only ancestor entries are ever touched, so the factorization produces zero
/// fill-in beyond the original branch-induced sparsity pattern.
///
/// `lambda` must be the expanded parent array from [`build_lambda`] and `n`
/// must equal its logical DOF count (`lambda.len() == n + 1`).
///
/// Reference: Featherstone 2008, Algorithm 9.5.
pub fn ltdl_factor_inplace(h: &mut [Vec<f64>], lambda: &[usize], n: usize) {
    for k in (1..=n).rev() {
        let mut i = lambda[k];
        while i != 0 {
            let a = h[k - 1][i - 1] / h[k - 1][k - 1];
            let mut j = i;
            while j != 0 {
                h[i - 1][j - 1] -= a * h[k - 1][j - 1];
                j = lambda[j];
            }
            h[k - 1][i - 1] = a;
            i = lambda[i];
        }
    }
}

/// Solve `(LᵀDL)·x = b` using the in-place `LᵀDL` factor produced by
/// [`ltdl_factor_inplace`].
///
/// The solve runs in three sparse passes over the `λ` chains:
///
/// 1. forward elimination against `Lᵀ` (ascending `k`),
/// 2. diagonal scaling by `D⁻¹`,
/// 3. back substitution against `L` (descending `k`).
///
/// Each pass only visits ancestor entries, so the whole solve is `O(N·d)`.
/// `h_factored` is the factored storage, `b` the right-hand side (length `n`),
/// and `lambda` the expanded parent array. The returned vector has length `n`.
pub fn ltdl_solve(h_factored: &[Vec<f64>], b: &[f64], lambda: &[usize], n: usize) -> Vec<f64> {
    let mut x = b.to_vec();
    // Convention note: the factor satisfies H = LᵀDL (verified by exact
    // reconstruction), so solving H·x = b means applying (Lᵀ)⁻¹, then D⁻¹,
    // then L⁻¹ in that order.
    //
    // The k-iteration *direction* of each elimination pass is dictated by which
    // triangle is being inverted, not by the H = LᵀDL vs H = LDLᵀ choice:
    //   * Lᵀ is upper-triangular ⇒ solving Lᵀ z = b is a *back*-substitution,
    //     so the dependent entry z[k] must be finalized before it updates its
    //     ancestors ⇒ k must run DESCENDING.
    //   * L is unit lower-triangular ⇒ solving L x = y is a *forward*-
    //     substitution ⇒ k must run ASCENDING.
    // (A short 2-level branch happens to tolerate either direction, which is why
    // a shallow test can mask the wrong order; a deep, multiply-branched tree
    // does not — hence the directions below.)

    // Phase 1: Lᵀ z = b  (back-substitution ⇒ k DESCENDING)
    for k in (1..=n).rev() {
        let mut i = lambda[k];
        while i != 0 {
            x[i - 1] -= h_factored[k - 1][i - 1] * x[k - 1];
            i = lambda[i];
        }
    }
    // Phase 2: D y = z
    for k in 1..=n {
        x[k - 1] /= h_factored[k - 1][k - 1];
    }
    // Phase 3: L x = y  (forward-substitution ⇒ k ASCENDING)
    for k in 1..=n {
        let mut i = lambda[k];
        while i != 0 {
            x[k - 1] -= h_factored[k - 1][i - 1] * x[i - 1];
            i = lambda[i];
        }
    }
    x
}

/// Error returned by [`SparseMassFactorization`] operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LtdlError {
    /// A diagonal pivot was ≤ 0 ⇒ matrix not positive-definite.
    NotPositiveDefinite {
        /// Zero-based pivot index where the non-positive pivot was found.
        index: usize,
    },
    /// A matrix/vector dimension did not match the expected DOF count.
    DimensionMismatch {
        /// Expected length / dimension.
        expected: usize,
        /// Actual length / dimension supplied.
        got: usize,
    },
}

impl std::fmt::Display for LtdlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotPositiveDefinite { index } => {
                write!(f, "matrix is not positive-definite at pivot index {index}")
            }
            Self::DimensionMismatch { expected, got } => {
                write!(f, "dimension mismatch: expected {expected}, got {got}")
            }
        }
    }
}

impl std::error::Error for LtdlError {}

/// Sparse LTDL factorization of the joint-space mass matrix H, exploiting the
/// branch-induced sparsity of a kinematic tree for O(Nd) solves (N=#DOF, d=tree depth).
#[derive(Clone, Debug)]
pub struct SparseMassFactorization {
    h_factored: Vec<Vec<f64>>,
    lambda: Vec<usize>,
    n: usize,
}

impl SparseMassFactorization {
    /// Factor an owned full symmetric SPD matrix `h` as `H = LᵀDL`.
    ///
    /// `h` must be `n×n` (checked) and `lambda` the expanded parent array from
    /// [`build_lambda`]. The diagonal is validated to be strictly positive both
    /// before and after the in-place factorization so that the result is a
    /// usable factor of an SPD matrix.
    ///
    /// # Errors
    /// Returns [`LtdlError::DimensionMismatch`] if `h` is not `n×n`, or
    /// [`LtdlError::NotPositiveDefinite`] if any diagonal pivot is non-positive
    /// before or after factorization.
    pub fn factor(h: Vec<Vec<f64>>, lambda: Vec<usize>, n: usize) -> Result<Self, LtdlError> {
        if h.len() != n {
            return Err(LtdlError::DimensionMismatch {
                expected: n,
                got: h.len(),
            });
        }
        for row in &h {
            if row.len() != n {
                return Err(LtdlError::DimensionMismatch {
                    expected: n,
                    got: row.len(),
                });
            }
        }
        // Pre-factor positive-definiteness check on the diagonal.
        for (k, row) in h.iter().enumerate() {
            if row[k] <= 0.0 {
                return Err(LtdlError::NotPositiveDefinite { index: k });
            }
        }
        let mut h_factored = h;
        ltdl_factor_inplace(&mut h_factored, &lambda, n);
        // Post-factor check: the D entries (diagonal) must stay positive.
        for (k, row) in h_factored.iter().enumerate() {
            if row[k] <= 0.0 {
                return Err(LtdlError::NotPositiveDefinite { index: k });
            }
        }
        Ok(Self {
            h_factored,
            lambda,
            n,
        })
    }

    /// Solve `H·x = b` for a single right-hand side using the stored factor.
    ///
    /// # Errors
    /// Returns [`LtdlError::DimensionMismatch`] if `b.len()` differs from the
    /// factored DOF count.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, LtdlError> {
        if b.len() != self.n {
            return Err(LtdlError::DimensionMismatch {
                expected: self.n,
                got: b.len(),
            });
        }
        Ok(ltdl_solve(&self.h_factored, b, &self.lambda, self.n))
    }

    /// Solve `H·X = B` for several right-hand-side columns at once.
    ///
    /// Each element of `b_cols` is one right-hand side; the result holds the
    /// corresponding solution columns in the same order.
    ///
    /// # Errors
    /// Returns [`LtdlError::DimensionMismatch`] if any column length differs
    /// from the factored DOF count.
    pub fn solve_columns(&self, b_cols: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, LtdlError> {
        b_cols.iter().map(|c| self.solve(c)).collect()
    }

    /// Return the expanded parent array `λ` used by this factorization.
    pub fn lambda(&self) -> &[usize] {
        &self.lambda
    }

    /// Return the number of DOF (the dimension `n` of the factored matrix).
    pub fn n_dof(&self) -> usize {
        self.n
    }

    /// Read a single entry of the in-place `L`/`D` storage (0-based indices).
    ///
    /// This exposes the raw factored storage: the diagonal holds `D[k]` and the
    /// strictly-lower entries hold the unit-lower factor `L`. Out-of-range
    /// indices return `0.0` rather than panicking.
    pub fn factored_entry(&self, row: usize, col: usize) -> f64 {
        if row >= self.n || col >= self.n {
            return 0.0;
        }
        self.h_factored[row][col]
    }
}

/// Compute the CRBA mass matrix at configuration `q` and immediately LTDL-factor it.
///
/// # Errors
/// Returns [`LtdlError`] if the resulting matrix is not SPD.
pub fn factor_crba(
    model: &mut ArticulatedModel,
    q: &[f64],
) -> Result<SparseMassFactorization, LtdlError> {
    use crate::crba::compute_mass_matrix_crba;
    let h = compute_mass_matrix_crba(model, q);
    let n = h.len();
    let lambda = build_lambda(model);
    SparseMassFactorization::factor(h, lambda, n)
}
