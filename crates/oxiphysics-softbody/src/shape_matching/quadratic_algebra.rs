// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dense 9×9 algebra for quadratic shape matching.
//!
//! Quadratic shape matching (Müller et al., 2005, §5.2) needs the inverse of
//! the moment matrix `Σ mᵢ q̃ᵢ q̃ᵢᵀ` of the nine-component feature vectors
//! `q̃ = [qx, qy, qz, qx², qy², qz², qx·qy, qy·qz, qz·qx]`. That matrix is only
//! positive *semi*-definite whenever the rest points do not excite all nine
//! quadratic modes: a body of fewer than nine particles, planar bodies, and
//! bodies whose particles sit on the corners of a box (`qx² = qy² = qz²`) all
//! give a singular matrix. The Moore–Penrose pseudo-inverse used here is the
//! exact inverse on the span of the feature vectors and zero elsewhere, which
//! is precisely what the least-squares fit of the goal positions requires.

/// Nine-component feature vector of a rest position.
pub(super) type Feature = [f64; 9];

/// Dense 9×9 matrix, row-major.
pub(super) type Mat9 = [[f64; 9]; 9];

/// Dense 3×9 matrix, row-major: maps a [`Feature`] to a position.
pub(super) type Mat3x9 = [[f64; 9]; 3];

/// Eigenvalues below this fraction of the largest one count as zero.
const RANK_TOLERANCE: f64 = 1.0e-12;

/// Maximum number of cyclic Jacobi sweeps (convergence takes 6–10 for 9×9).
const MAX_SWEEPS: usize = 64;

/// Rotate columns `p` and `q` of `m` by the plane rotation `(c, s)`:
/// `m[:, p] ← c m[:, p] − s m[:, q]`, `m[:, q] ← s m[:, p] + c m[:, q]`.
fn rotate_columns(m: &mut Mat9, p: usize, q: usize, c: f64, s: f64) {
    for row in m.iter_mut() {
        let (x, y) = (row[p], row[q]);
        row[p] = c * x - s * y;
        row[q] = s * x + c * y;
    }
}

/// Rotate rows `p < q` of `m` by the plane rotation `(c, s)`.
fn rotate_rows(m: &mut Mat9, p: usize, q: usize, c: f64, s: f64) {
    let (head, tail) = m.split_at_mut(q);
    let (row_p, row_q) = (&mut head[p], &mut tail[0]);
    for (x_ref, y_ref) in row_p.iter_mut().zip(row_q.iter_mut()) {
        let (x, y) = (*x_ref, *y_ref);
        *x_ref = c * x - s * y;
        *y_ref = s * x + c * y;
    }
}

/// Eigen-decomposition `m = V diag(λ) Vᵀ` of a symmetric 9×9 matrix by the
/// cyclic Jacobi method. Returns the eigenvalues and the matrix `V` whose
/// *columns* are the corresponding orthonormal eigenvectors.
fn symmetric_eigen(m: &Mat9) -> ([f64; 9], Mat9) {
    // Enforce exact symmetry so that rounding in the input cannot make the
    // rotations inconsistent.
    let mut a = [[0.0_f64; 9]; 9];
    for (i, row) in a.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = 0.5 * (m[i][j] + m[j][i]);
        }
    }
    let mut v = [[0.0_f64; 9]; 9];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let frobenius_sq: f64 = a.iter().flatten().map(|x| x * x).sum();
    if frobenius_sq == 0.0 {
        return ([0.0; 9], v);
    }
    for _ in 0..MAX_SWEEPS {
        let off_sq: f64 = 2.0
            * a.iter()
                .enumerate()
                .flat_map(|(p, row)| row.iter().skip(p + 1).map(|x| x * x))
                .sum::<f64>();
        if off_sq <= 1.0e-30 * frobenius_sq {
            break;
        }
        for p in 0..8 {
            for q in (p + 1)..9 {
                let a_pq = a[p][q];
                if a_pq == 0.0 {
                    continue;
                }
                // Rotation that annihilates a[p][q]: tan θ is the smaller root
                // of t² + 2 t ϑ − 1 = 0 with ϑ = (a_qq − a_pp) / (2 a_pq).
                let theta = (a[q][q] - a[p][p]) / (2.0 * a_pq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                rotate_columns(&mut a, p, q, c, s);
                rotate_rows(&mut a, p, q, c, s);
                rotate_columns(&mut v, p, q, c, s);
            }
        }
    }
    let mut eigenvalues = [0.0_f64; 9];
    for (i, value) in eigenvalues.iter_mut().enumerate() {
        *value = a[i][i];
    }
    (eigenvalues, v)
}

/// Moore–Penrose pseudo-inverse of a symmetric positive-semidefinite 9×9
/// matrix.
///
/// Eigenvalues not larger than `1e-12` times the largest one are treated as
/// zero. The result `M⁺` satisfies `M M⁺ M = M` and `M⁺ M M⁺ = M⁺`; `M M⁺` is
/// the orthogonal projector onto the range of `M`, so `M M⁺ x = x` for every
/// `x` in the span of the data the moment matrix was assembled from.
pub(super) fn symmetric_pseudo_inverse(m: &Mat9) -> Mat9 {
    let (eigenvalues, v) = symmetric_eigen(m);
    let largest = eigenvalues.iter().copied().fold(0.0_f64, f64::max);
    let cutoff = RANK_TOLERANCE * largest;
    let mut inverse = [[0.0_f64; 9]; 9];
    for (k, &lambda) in eigenvalues.iter().enumerate() {
        if lambda <= cutoff || lambda <= 0.0 {
            continue;
        }
        let weight = 1.0 / lambda;
        for (i, row) in inverse.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry += weight * v[i][k] * v[j][k];
            }
        }
    }
    inverse
}

/// Mass-weighted moment matrix `Σ mᵢ q̃ᵢ q̃ᵢᵀ` of feature vectors.
pub(super) fn feature_moment_matrix(features: &[Feature], masses: &[f64]) -> Mat9 {
    let mut moment = [[0.0_f64; 9]; 9];
    for (feature, &mass) in features.iter().zip(masses.iter()) {
        for (row, &f_row) in moment.iter_mut().zip(feature.iter()) {
            for (entry, &f_col) in row.iter_mut().zip(feature.iter()) {
                *entry += mass * f_row * f_col;
            }
        }
    }
    moment
}

/// Product `a · b` of a 3×9 and a 9×9 matrix.
pub(super) fn mul_3x9_9x9(a: &Mat3x9, b: &Mat9) -> Mat3x9 {
    let mut product = [[0.0_f64; 9]; 3];
    for (out_row, a_row) in product.iter_mut().zip(a.iter()) {
        for (k, &a_ik) in a_row.iter().enumerate() {
            for (out, &b_kj) in out_row.iter_mut().zip(b[k].iter()) {
                *out += a_ik * b_kj;
            }
        }
    }
    product
}

/// Product `a · x` of a 3×9 matrix and a feature vector.
pub(super) fn mul_3x9_vec(a: &Mat3x9, x: &Feature) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for (out, row) in result.iter_mut().zip(a.iter()) {
        *out = row.iter().zip(x.iter()).map(|(a_ij, x_j)| a_ij * x_j).sum();
    }
    result
}
