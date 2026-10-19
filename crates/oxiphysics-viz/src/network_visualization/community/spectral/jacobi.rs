// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eigen-decomposition of small dense symmetric matrices by cyclic Jacobi
//! rotations.
//!
//! The Lanczos iteration reduces the Laplacian to a tridiagonal matrix of at
//! most a few dozen rows; Jacobi's method is simple, unconditionally stable and
//! accurate to working precision for matrices of that size.

/// Bound on the sweeps over all off-diagonal entries.  Jacobi converges
/// quadratically, so a handful of sweeps suffice for any realistic input.
const MAX_SWEEPS: usize = 64;

/// The off-diagonal part counts as zero once its squared Frobenius norm is
/// below this fraction of the squared norm of the whole matrix, i.e. once the
/// off-diagonal entries are at rounding level (`1e-15` relative).
const NEGLIGIBLE_OFF_DIAGONAL: f64 = 1e-30;

/// Eigenvalues (ascending) and unit eigenvectors of the symmetric matrix `a`.
///
/// `vectors[i]` is the eigenvector that belongs to `values[i]`.  Only the
/// symmetric part of `a` is meaningful; `a` must be square.
pub(super) fn symmetric_eigen(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    // Accumulated rotations; column `i` converges to the eigenvector of `a[i][i]`.
    let mut rotations: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _ in 0..MAX_SWEEPS {
        let mut off = 0.0_f64;
        let mut diagonal = 0.0_f64;
        for (i, row) in a.iter().enumerate() {
            for (j, &entry) in row.iter().enumerate() {
                if i == j {
                    diagonal += entry * entry;
                } else {
                    off += entry * entry;
                }
            }
        }
        let total = diagonal + off;
        if off <= NEGLIGIBLE_OFF_DIAGONAL * total {
            break;
        }
        for p in 0..n.saturating_sub(1) {
            for q in p + 1..n {
                // Entries already at rounding level are not worth a rotation.
                if a[p][q] * a[p][q] > NEGLIGIBLE_OFF_DIAGONAL * total * 1e-6 {
                    rotate(&mut a, &mut rotations, p, q);
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.iter().map(|&i| a[i][i]).collect();
    let vectors = order
        .iter()
        .map(|&i| rotations.iter().map(|row| row[i]).collect())
        .collect();
    (values, vectors)
}

/// Apply the Jacobi rotation in the `(p, q)` plane that annihilates `a[p][q]`
/// to `a` (from both sides) and accumulate it in `rotations`.
fn rotate(a: &mut [Vec<f64>], rotations: &mut [Vec<f64>], p: usize, q: usize) {
    let apq = a[p][q];
    if apq == 0.0 {
        return;
    }
    // tan(theta) is the smaller root of t^2 + 2 t * cot(2 theta) - 1 = 0.
    let cot_2theta = (a[q][q] - a[p][p]) / (2.0 * apq);
    let t = cot_2theta.signum() / (cot_2theta.abs() + cot_2theta.hypot(1.0));
    let c = 1.0 / t.hypot(1.0);
    let s = t * c;
    // Columns p and q of `a`, then rows p and q (A <- J^T A J), then the
    // accumulated eigenvector matrix (V <- V J).
    for row in a.iter_mut() {
        let (x, y) = (row[p], row[q]);
        row[p] = c * x - s * y;
        row[q] = s * x + c * y;
    }
    let row_p = a[p].clone();
    let row_q = a[q].clone();
    for (k, (&x, &y)) in row_p.iter().zip(&row_q).enumerate() {
        a[p][k] = c * x - s * y;
        a[q][k] = s * x + c * y;
    }
    for row in rotations.iter_mut() {
        let (x, y) = (row[p], row[q]);
        row[p] = c * x - s * y;
        row[q] = s * x + c * y;
    }
}
