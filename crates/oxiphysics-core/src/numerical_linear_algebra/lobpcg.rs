//! LOBPCG: Locally Optimal Block Preconditioned Conjugate Gradient (Knyazev, 2001)
//! for the smallest eigenpairs of a symmetric sparse matrix.
//!
//! Each iteration builds the search space `S = [X, W, P]` from
//!
//! * `X` -- the current Ritz vectors (a block slightly larger than the number
//!   of wanted pairs, which speeds up convergence of the last wanted pair),
//! * `W` -- the Jacobi-preconditioned residuals `D^-1 (A X - X Lambda)` of the
//!   pairs that have not converged yet,
//! * `P` -- the previous search directions,
//!
//! orthonormalises it in that order (dropping numerically dependent `W` and `P`
//! columns, never `X` columns), solves the Rayleigh-Ritz problem `S^T A S` with a
//! dense symmetric eigen-solver and keeps the smallest Ritz pairs.  The new `P`
//! is the `W`/`P` part of `S Y`.  Convergence is always confirmed on residuals
//! recomputed from scratch, so drift in the recurrences cannot cause a false
//! `converged = true`.

use super::functions::{axpy, dot, norm2};
use super::krylov_eigen::{VectorGenerator, dot_blocked};
use super::types::{EigenResult, SparseCSR};
use crate::numerical_methods::symmetric_eigen_n;

/// Seed of the generator that produces the initial block.
const INITIAL_BLOCK_SEED: u64 = 0x1F2E_3D4C_5B6A_7988;
/// A direction whose norm falls below this fraction of its original norm during
/// orthogonalisation is numerically dependent and is dropped.
const DEPENDENCE_TOLERANCE: f64 = 1.0e-10;
/// Upper limit on the number of extra block vectors beyond the wanted ones.
const MAX_EXTRA_BLOCK_VECTORS: usize = 5;

/// Residual data of a block of Ritz pairs.
struct Residuals {
    /// Residual vectors `A x_i - lambda_i x_i`.
    vectors: Vec<Vec<f64>>,
    /// Euclidean norms of the residual vectors.
    norms: Vec<f64>,
    /// Acceptance thresholds `tol * max(1, |lambda_i|)`.
    thresholds: Vec<f64>,
}

impl Residuals {
    fn new(x: &[Vec<f64>], ax: &[Vec<f64>], lambda: &[f64], tol: f64) -> Self {
        let vectors: Vec<Vec<f64>> = x
            .iter()
            .zip(ax.iter())
            .zip(lambda.iter())
            .map(|((xi, axi), &li)| {
                let mut r = axi.clone();
                axpy(-li, xi, &mut r);
                r
            })
            .collect();
        let norms = vectors.iter().map(|r| norm2(r)).collect();
        let thresholds = lambda.iter().map(|l| tol * l.abs().max(1.0)).collect();
        Self {
            vectors,
            norms,
            thresholds,
        }
    }

    /// Whether the first `wanted` pairs satisfy their thresholds (a NaN residual
    /// never does).
    fn wanted_converged(&self, wanted: usize) -> bool {
        self.norms
            .iter()
            .zip(self.thresholds.iter())
            .take(wanted)
            .all(|(norm, threshold)| norm <= threshold)
            && self.norms.len() >= wanted
    }

    /// Indices of the pairs that still need work.
    fn active(&self) -> Vec<usize> {
        (0..self.norms.len())
            .filter(|&i| {
                let converged = self.norms[i] <= self.thresholds[i];
                !converged
            })
            .collect()
    }
}

/// `1 / |a_ii|` for every row (1 where the diagonal entry vanishes): the Jacobi
/// preconditioner, kept positive so that it stays a valid preconditioner for
/// indefinite matrices.
fn jacobi_inverse_diagonal(a: &SparseCSR) -> Vec<f64> {
    a.diagonal()
        .iter()
        .map(|&d| {
            if d.abs() > 1.0e-300 {
                1.0 / d.abs()
            } else {
                1.0
            }
        })
        .collect()
}

/// Applies `a` to every vector of `block`.
fn matvec_block(a: &SparseCSR, block: &[Vec<f64>]) -> Vec<Vec<f64>> {
    block.iter().map(|v| a.matvec(v)).collect()
}

/// Orthonormalises `candidates` against the orthonormal vectors in `against`
/// and against each other (two Gram-Schmidt passes), dropping every candidate
/// whose norm collapses to [`DEPENDENCE_TOLERANCE`] of its original value.
fn orthonormalise_block(candidates: Vec<Vec<f64>>, against: &[&[Vec<f64>]]) -> Vec<Vec<f64>> {
    let mut kept: Vec<Vec<f64>> = Vec::with_capacity(candidates.len());
    for mut v in candidates {
        let original = norm2(&v);
        if original == 0.0 || !original.is_finite() {
            continue;
        }
        for _pass in 0..2 {
            for block in against {
                for q in block.iter() {
                    let c = dot_blocked(q, &v);
                    axpy(-c, q, &mut v);
                }
            }
            for q in &kept {
                let c = dot_blocked(q, &v);
                axpy(-c, q, &mut v);
            }
        }
        let remaining = norm2(&v);
        if remaining <= DEPENDENCE_TOLERANCE * original {
            continue;
        }
        v.iter_mut().for_each(|e| *e /= remaining);
        kept.push(v);
    }
    kept
}

/// Rayleigh-Ritz on the orthonormal columns `s` with images `a_s = A s`: returns
/// the ascending Ritz values and the eigenvectors `y[j][i]` of `S^T A S`.
fn rayleigh_ritz(s: &[Vec<f64>], a_s: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let size = s.len();
    let mut gram = vec![vec![0.0_f64; size]; size];
    for i in 0..size {
        for j in i..size {
            // The matrix is symmetric, so only the upper triangle is formed and
            // mirrored; this keeps the Gram matrix exactly symmetric.
            let value = dot_blocked(&s[i], &a_s[j]);
            gram[i][j] = value;
            gram[j][i] = value;
        }
    }
    symmetric_eigen_n(&gram)
}

/// `count` linear combinations of `columns[range]` with the coefficients
/// `y[j][i]` (column `i` of `S`, combination `j`).
fn combine(
    columns: &[Vec<f64>],
    y: &[Vec<f64>],
    first_column: usize,
    count: usize,
) -> Vec<Vec<f64>> {
    let n = columns.first().map_or(0, Vec::len);
    y.iter()
        .take(count)
        .map(|coefficients| {
            let mut combined = vec![0.0_f64; n];
            for (column, &c) in columns.iter().zip(coefficients.iter()).skip(first_column) {
                axpy(c, column, &mut combined);
            }
            combined
        })
        .collect()
}

/// Smallest `wanted` eigenpairs by a dense symmetric eigen-decomposition.
fn dense_smallest(a: &SparseCSR, wanted: usize, tol: f64) -> EigenResult {
    let n = a.nrows;
    let mut dense = vec![vec![0.0_f64; n]; n];
    for (i, bounds) in a.row_ptr.windows(2).enumerate() {
        for idx in bounds[0]..bounds[1] {
            let j = a.col_idx[idx];
            let half = 0.5 * a.values[idx];
            dense[i][j] += half;
            dense[j][i] += half;
        }
    }
    let (values, vectors) = symmetric_eigen_n(&dense);
    let eigenvalues: Vec<f64> = values.into_iter().take(wanted).collect();
    let eigenvectors: Vec<Vec<f64>> = vectors.into_iter().take(wanted).collect();
    let ax = matvec_block(a, &eigenvectors);
    let residuals = Residuals::new(&eigenvectors, &ax, &eigenvalues, tol);
    EigenResult {
        converged: residuals.wanted_converged(wanted),
        eigenvalues,
        eigenvectors,
        iterations: 1,
    }
}

/// LOBPCG for the `k` smallest eigenpairs; see [`super::EigenSolver::lobpcg`].
pub(super) fn lobpcg(a: &SparseCSR, k: usize, tol: f64, max_iter: usize) -> EigenResult {
    let n = a.nrows;
    if a.ncols != n {
        // Eigenpairs are only defined for square matrices.
        return EigenResult {
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
            iterations: 0,
            converged: false,
        };
    }
    let wanted = k.min(n);
    if wanted == 0 {
        return EigenResult {
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
            iterations: 0,
            converged: true,
        };
    }
    let block = (wanted + wanted.clamp(1, MAX_EXTRA_BLOCK_VECTORS)).min(n);
    if 3 * block >= n {
        return dense_smallest(a, wanted, tol);
    }

    let inv_diag = jacobi_inverse_diagonal(a);
    let mut generator = VectorGenerator::new(INITIAL_BLOCK_SEED);
    let mut x: Vec<Vec<f64>> = Vec::with_capacity(block);
    while x.len() < block {
        let missing = block - x.len();
        let candidates: Vec<Vec<f64>> = (0..missing).map(|_| generator.unit_vector(n)).collect();
        let fresh = orthonormalise_block(candidates, &[&x]);
        x.extend(fresh);
    }
    let mut ax = matvec_block(a, &x);
    let (initial_values, initial_y) = rayleigh_ritz(&x, &ax);
    let mut lambda: Vec<f64> = initial_values.into_iter().take(block).collect();
    x = combine(&x, &initial_y, 0, block);
    ax = combine(&ax, &initial_y, 0, block);
    let mut directions: Vec<Vec<f64>> = Vec::new();

    for iteration in 0..max_iter {
        let mut residuals = Residuals::new(&x, &ax, &lambda, tol);
        if residuals.wanted_converged(wanted) {
            // Confirm on residuals computed from scratch before declaring victory.
            ax = matvec_block(a, &x);
            lambda = x
                .iter()
                .zip(ax.iter())
                .map(|(xi, axi)| dot(xi, axi))
                .collect();
            residuals = Residuals::new(&x, &ax, &lambda, tol);
            if residuals.wanted_converged(wanted) {
                return finish(x, lambda, wanted, iteration, true);
            }
        }
        let active = residuals.active();

        let preconditioned: Vec<Vec<f64>> = active
            .iter()
            .map(|&i| {
                residuals.vectors[i]
                    .iter()
                    .zip(inv_diag.iter())
                    .map(|(r, d)| r * d)
                    .collect()
            })
            .collect();
        let w_block = orthonormalise_block(preconditioned, &[&x]);
        let previous: Vec<Vec<f64>> = if directions.is_empty() {
            Vec::new()
        } else {
            active.iter().map(|&i| directions[i].clone()).collect()
        };
        let p_block = orthonormalise_block(previous, &[&x, &w_block]);

        let aw = matvec_block(a, &w_block);
        let ap = matvec_block(a, &p_block);
        let mut s = x.clone();
        s.extend(w_block);
        s.extend(p_block);
        let mut a_s = ax.clone();
        a_s.extend(aw);
        a_s.extend(ap);

        let (theta, y) = rayleigh_ritz(&s, &a_s);
        x = combine(&s, &y, 0, block);
        ax = combine(&a_s, &y, 0, block);
        directions = combine(&s, &y, block, block);
        lambda = theta.into_iter().take(block).collect();
    }

    // Out of iterations: report the honest state.
    ax = matvec_block(a, &x);
    lambda = x
        .iter()
        .zip(ax.iter())
        .map(|(xi, axi)| dot(xi, axi))
        .collect();
    let residuals = Residuals::new(&x, &ax, &lambda, tol);
    let converged = residuals.wanted_converged(wanted);
    finish(x, lambda, wanted, max_iter, converged)
}

/// Packages the first `wanted` Ritz pairs (the block is sorted by eigenvalue).
fn finish(
    x: Vec<Vec<f64>>,
    lambda: Vec<f64>,
    wanted: usize,
    iterations: usize,
    converged: bool,
) -> EigenResult {
    EigenResult {
        eigenvalues: lambda.into_iter().take(wanted).collect(),
        eigenvectors: x.into_iter().take(wanted).collect(),
        iterations,
        converged,
    }
}
