//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::math::{Mat3, Real, Vec3};

pub(super) const EPS: Real = 1e-14;
/// Determinant of a 3x3 matrix.
pub fn det3(m: Mat3) -> Real {
    m.m11 * (m.m22 * m.m33 - m.m23 * m.m32) - m.m12 * (m.m21 * m.m33 - m.m23 * m.m31)
        + m.m13 * (m.m21 * m.m32 - m.m22 * m.m31)
}
/// Trace of a 3x3 matrix.
pub fn trace3(m: Mat3) -> Real {
    m.m11 + m.m22 + m.m33
}
/// Frobenius norm of a 3x3 matrix.
pub fn frobenius_norm3(m: Mat3) -> Real {
    (m.m11 * m.m11
        + m.m12 * m.m12
        + m.m13 * m.m13
        + m.m21 * m.m21
        + m.m22 * m.m22
        + m.m23 * m.m23
        + m.m31 * m.m31
        + m.m32 * m.m32
        + m.m33 * m.m33)
        .sqrt()
}
/// 1-norm of a 3x3 matrix (maximum column sum of absolute values).
pub fn norm1_3(m: Mat3) -> Real {
    let c0 = m.m11.abs() + m.m21.abs() + m.m31.abs();
    let c1 = m.m12.abs() + m.m22.abs() + m.m32.abs();
    let c2 = m.m13.abs() + m.m23.abs() + m.m33.abs();
    c0.max(c1).max(c2)
}
/// Infinity-norm of a 3x3 matrix (maximum row sum of absolute values).
pub fn norm_inf3(m: Mat3) -> Real {
    let r0 = m.m11.abs() + m.m12.abs() + m.m13.abs();
    let r1 = m.m21.abs() + m.m22.abs() + m.m23.abs();
    let r2 = m.m31.abs() + m.m32.abs() + m.m33.abs();
    r0.max(r1).max(r2)
}
/// Check if a matrix is approximately symmetric (max |M - M^T| < tol).
pub fn is_symmetric3(m: Mat3, tol: Real) -> bool {
    (m.m12 - m.m21).abs() < tol && (m.m13 - m.m31).abs() < tol && (m.m23 - m.m32).abs() < tol
}
/// Check if a matrix is approximately the identity.
pub fn is_identity3(m: Mat3, tol: Real) -> bool {
    (m - Mat3::identity()).abs().max() < tol
}
/// Check if a matrix is approximately orthogonal: M^T * M ~ I.
pub fn is_orthogonal3(m: Mat3, tol: Real) -> bool {
    is_identity3(m.transpose() * m, tol)
}
/// Inverse of a 3x3 matrix. Returns `None` if the matrix is singular.
pub fn inv3(m: Mat3) -> Option<Mat3> {
    let d = det3(m);
    if d.abs() < EPS {
        return None;
    }
    let inv_d = 1.0 / d;
    Some(Mat3::new(
        (m.m22 * m.m33 - m.m23 * m.m32) * inv_d,
        (m.m13 * m.m32 - m.m12 * m.m33) * inv_d,
        (m.m12 * m.m23 - m.m13 * m.m22) * inv_d,
        (m.m23 * m.m31 - m.m21 * m.m33) * inv_d,
        (m.m11 * m.m33 - m.m13 * m.m31) * inv_d,
        (m.m13 * m.m21 - m.m11 * m.m23) * inv_d,
        (m.m21 * m.m32 - m.m22 * m.m31) * inv_d,
        (m.m12 * m.m31 - m.m11 * m.m32) * inv_d,
        (m.m11 * m.m22 - m.m12 * m.m21) * inv_d,
    ))
}
/// Solve a 3x3 linear system `M * x = b` via Cramer's rule.
/// Returns `None` if the matrix is singular.
pub fn solve3(m: Mat3, b: Vec3) -> Option<Vec3> {
    let d = det3(m);
    if d.abs() < EPS {
        return None;
    }
    let inv_d = 1.0 / d;
    let d0 = b.x * (m.m22 * m.m33 - m.m23 * m.m32) - m.m12 * (b.y * m.m33 - m.m23 * b.z)
        + m.m13 * (b.y * m.m32 - m.m22 * b.z);
    let d1 = m.m11 * (b.y * m.m33 - m.m23 * b.z) - b.x * (m.m21 * m.m33 - m.m23 * m.m31)
        + m.m13 * (m.m21 * b.z - b.y * m.m31);
    let d2 = m.m11 * (m.m22 * b.z - b.y * m.m32) - m.m12 * (m.m21 * b.z - b.y * m.m31)
        + b.x * (m.m21 * m.m32 - m.m22 * m.m31);
    Some(Vec3::new(d0 * inv_d, d1 * inv_d, d2 * inv_d))
}
/// Characteristic polynomial coefficients of a 3x3 matrix.
///
/// det(lambda*I - M) = lambda^3 - c2*lambda^2 + c1*lambda - c0
/// Returns `(c0, c1, c2)`.
pub fn characteristic_poly3(m: Mat3) -> (Real, Real, Real) {
    let c2 = trace3(m);
    let c1 = (m.m11 * m.m22 - m.m12 * m.m21)
        + (m.m11 * m.m33 - m.m13 * m.m31)
        + (m.m22 * m.m33 - m.m23 * m.m32);
    let c0 = det3(m);
    (c0, c1, c2)
}
/// Analytical eigenvalues of a symmetric 3x3 matrix using Cardano's formula.
/// Returns eigenvalues in ascending order.
pub fn symmetric_eigenvalues3(m: Mat3) -> Vec3 {
    let (c0, c1, c2) = characteristic_poly3(m);
    let p_coeff = c1 - c2 * c2 / 3.0;
    let q_coeff = 2.0 * c2 * c2 * c2 / 27.0 - c1 * c2 / 3.0 + c0;
    let discriminant = -(4.0 * p_coeff * p_coeff * p_coeff + 27.0 * q_coeff * q_coeff);
    let shift = c2 / 3.0;
    let mut eigs = if discriminant >= 0.0 && p_coeff < 0.0 {
        let m_val = 2.0 * (-p_coeff / 3.0).sqrt();
        let arg = (3.0 * q_coeff / (p_coeff * m_val)).clamp(-1.0, 1.0);
        let theta = arg.acos() / 3.0;
        use std::f64::consts::PI;
        let l0 = m_val * theta.cos() + shift;
        let l1 = m_val * (theta - 2.0 * PI / 3.0).cos() + shift;
        let l2 = m_val * (theta - 4.0 * PI / 3.0).cos() + shift;
        [l0, l1, l2]
    } else {
        [shift, shift, shift]
    };
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Vec3::new(eigs[0], eigs[1], eigs[2])
}
/// Symmetric 3x3 matrix eigendecomposition via Jacobi iterations.
///
/// Returns `(eigenvalues, eigenvectors)` where columns of `eigenvectors` are
/// the corresponding unit eigenvectors. Input `m` should be symmetric.
/// Eigenvalues are returned in ascending order.
///
/// Algorithm: cyclic Jacobi -- for each off-diagonal element (p,q), compute a
/// Givens rotation J that zeroes it, then update A <- J^T A J and accumulate
/// V <- V J.  After convergence the diagonal of A holds the eigenvalues and
/// the columns of V are the eigenvectors.
pub fn symmetric_eigen3(m: Mat3) -> (Vec3, Mat3) {
    pub(super) const MAX_SWEEPS: usize = 50;
    let mut a = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = m[(i, j)];
        }
    }
    let mut v = [[0.0_f64; 3]; 3];
    v[0][0] = 1.0;
    v[1][1] = 1.0;
    v[2][2] = 1.0;
    let pairs = [(0, 1), (0, 2), (1, 2)];
    'outer: for _ in 0..MAX_SWEEPS {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off < EPS {
            break 'outer;
        }
        for &(p, q) in &pairs {
            let apq = a[p][q];
            if apq.abs() < EPS {
                continue;
            }
            let app = a[p][p];
            let aqq = a[q][q];
            let tau = (aqq - app) / (2.0 * apq);
            let t = if tau >= 0.0 {
                1.0 / (tau + (1.0 + tau * tau).sqrt())
            } else {
                1.0 / (tau - (1.0 + tau * tau).sqrt())
            };
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = t * c;
            a[p][p] = app - t * apq;
            a[q][q] = aqq + t * apq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            let r = 3 - p - q;
            let arp = a[r][p];
            let arq = a[r][q];
            let new_arp = c * arp - s * arq;
            let new_arq = s * arp + c * arq;
            a[r][p] = new_arp;
            a[p][r] = new_arp;
            a[r][q] = new_arq;
            a[q][r] = new_arq;
            for k in 0..3 {
                let vkp = v[k][p];
                let vkq = v[k][q];
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    let diag = [a[0][0], a[1][1], a[2][2]];
    let mut pairs_sorted: [(Real, usize); 3] = [(diag[0], 0), (diag[1], 1), (diag[2], 2)];
    pairs_sorted.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let eigenvalues = Vec3::new(pairs_sorted[0].0, pairs_sorted[1].0, pairs_sorted[2].0);
    let col = |i: usize| -> Vec3 { Vec3::new(v[0][i], v[1][i], v[2][i]) };
    let eigenvectors = Mat3::from_columns(&[
        col(pairs_sorted[0].1),
        col(pairs_sorted[1].1),
        col(pairs_sorted[2].1),
    ]);
    (eigenvalues, eigenvectors)
}
/// LU decomposition of a 3x3 matrix with partial pivoting.
///
/// Returns `(L, U, perm)` where `P * M = L * U`.
/// `perm` is the row permutation \[p0, p1, p2\] such that row `perm[i]` of `M`
/// becomes row `i` of `P * M`.
///
/// Returns `None` if the matrix is singular.
pub fn lu_decomp3(m: Mat3) -> Option<(Mat3, Mat3, [usize; 3])> {
    let mut a = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = m[(i, j)];
        }
    }
    let mut perm = [0usize, 1, 2];
    for k in 0..3 {
        let mut max_val = a[k][k].abs();
        let mut max_row = k;
        for i in (k + 1)..3 {
            if a[i][k].abs() > max_val {
                max_val = a[i][k].abs();
                max_row = i;
            }
        }
        if max_val < EPS {
            return None;
        }
        if max_row != k {
            a.swap(k, max_row);
            perm.swap(k, max_row);
        }
        for i in (k + 1)..3 {
            let factor = a[i][k] / a[k][k];
            a[i][k] = factor;
            for j in (k + 1)..3 {
                a[i][j] -= factor * a[k][j];
            }
        }
    }
    let l = Mat3::new(1.0, 0.0, 0.0, a[1][0], 1.0, 0.0, a[2][0], a[2][1], 1.0);
    let u = Mat3::new(
        a[0][0], a[0][1], a[0][2], 0.0, a[1][1], a[1][2], 0.0, 0.0, a[2][2],
    );
    Some((l, u, perm))
}
/// Solve `M * x = b` using LU decomposition.
///
/// Returns `None` if the matrix is singular.
pub fn lu_solve3(m: Mat3, b: Vec3) -> Option<Vec3> {
    let (l, u, perm) = lu_decomp3(m)?;
    let pb = Vec3::new(b[perm[0]], b[perm[1]], b[perm[2]]);
    let y = forward_substitute3(l, pb)?;
    backward_substitute3(u, y)
}
/// QR decomposition of a 3x3 matrix via modified Gram-Schmidt.
///
/// Returns `(Q, R)` where `Q` is orthogonal and `R` is upper triangular.
/// If a column becomes near-zero, a fallback direction is used.
pub fn qr_decomp3(m: Mat3) -> (Mat3, Mat3) {
    let cols: [Vec3; 3] = [
        Vec3::new(m[(0, 0)], m[(1, 0)], m[(2, 0)]),
        Vec3::new(m[(0, 1)], m[(1, 1)], m[(2, 1)]),
        Vec3::new(m[(0, 2)], m[(1, 2)], m[(2, 2)]),
    ];
    let mut q_cols = [Vec3::zeros(); 3];
    let mut r = Mat3::zeros();
    for j in 0..3 {
        let mut u = cols[j];
        for i in 0..j {
            let proj = q_cols[i].dot(&cols[j]);
            r[(i, j)] = proj;
            u -= q_cols[i] * proj;
        }
        let norm = u.norm();
        r[(j, j)] = norm;
        if norm > EPS {
            q_cols[j] = u / norm;
        } else {
            q_cols[j] = fallback_orthogonal(j, &q_cols);
        }
    }
    let q = Mat3::from_columns(&q_cols);
    (q, r)
}
/// QR decomposition via Householder reflections.
///
/// Returns `(Q, R)` where `Q` is orthogonal and `R` is upper triangular.
/// Householder is more numerically stable than Gram-Schmidt for ill-conditioned
/// matrices.
pub fn qr_householder3(m: Mat3) -> (Mat3, Mat3) {
    let mut a = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = m[(i, j)];
        }
    }
    let mut q = [[0.0_f64; 3]; 3];
    q[0][0] = 1.0;
    q[1][1] = 1.0;
    q[2][2] = 1.0;
    for k in 0..2 {
        let mut x = [0.0_f64; 3];
        let mut norm_x = 0.0_f64;
        for i in k..3 {
            x[i] = a[i][k];
            norm_x += x[i] * x[i];
        }
        norm_x = norm_x.sqrt();
        if norm_x < EPS {
            continue;
        }
        let sign = if a[k][k] >= 0.0 { 1.0 } else { -1.0 };
        let mut v = [0.0_f64; 3];
        v[k..3].copy_from_slice(&x[k..3]);
        v[k] += sign * norm_x;
        let mut v_norm = 0.0_f64;
        for i in k..3 {
            v_norm += v[i] * v[i];
        }
        if v_norm < EPS {
            continue;
        }
        for j in k..3 {
            let mut dot_val = 0.0_f64;
            for i in k..3 {
                dot_val += v[i] * a[i][j];
            }
            let factor = 2.0 * dot_val / v_norm;
            for i in k..3 {
                a[i][j] -= factor * v[i];
            }
        }
        for j in 0..3 {
            let mut dot_val = 0.0_f64;
            for i in k..3 {
                dot_val += q[j][i] * v[i];
            }
            let factor = 2.0 * dot_val / v_norm;
            for i in k..3 {
                q[j][i] -= factor * v[i];
            }
        }
    }
    let q_mat = Mat3::new(
        q[0][0], q[0][1], q[0][2], q[1][0], q[1][1], q[1][2], q[2][0], q[2][1], q[2][2],
    );
    let r_mat = Mat3::new(
        a[0][0], a[0][1], a[0][2], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2],
    );
    (q_mat, r_mat)
}
/// Return a unit vector orthogonal to all *filled* columns of `cols`.
pub(super) fn fallback_orthogonal_filled(
    skip_j: usize,
    cols: &[Vec3; 3],
    filled: &[bool; 3],
) -> Vec3 {
    let candidates = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    for c in &candidates {
        let mut u = *c;
        for k in 0..3 {
            if k == skip_j || !filled[k] {
                continue;
            }
            u -= cols[k] * cols[k].dot(c);
        }
        let n = u.norm();
        if n > EPS {
            return u / n;
        }
    }
    Vec3::new(0.0, 0.0, 1.0)
}
/// Return a unit vector orthogonal to the first `j` columns of `q_cols`.
pub(super) fn fallback_orthogonal(j: usize, q_cols: &[Vec3; 3]) -> Vec3 {
    let candidates = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    for c in &candidates {
        let mut u = *c;
        for k in 0..j {
            u -= q_cols[k] * q_cols[k].dot(c);
        }
        let n = u.norm();
        if n > EPS {
            return u / n;
        }
    }
    Vec3::new(0.0, 0.0, 1.0)
}
/// Forward substitution: solve L * x = b where L is lower triangular 3x3.
///
/// Returns `None` if a diagonal element is zero.
pub fn forward_substitute3(l: Mat3, b: Vec3) -> Option<Vec3> {
    let mut x = [0.0_f64; 3];
    for i in 0..3 {
        let mut sum = b[i];
        for j in 0..i {
            sum -= l[(i, j)] * x[j];
        }
        if l[(i, i)].abs() < EPS {
            return None;
        }
        x[i] = sum / l[(i, i)];
    }
    Some(Vec3::new(x[0], x[1], x[2]))
}
/// Backward substitution: solve U * x = b where U is upper triangular 3x3.
///
/// Returns `None` if a diagonal element is zero.
pub fn backward_substitute3(u: Mat3, b: Vec3) -> Option<Vec3> {
    let mut x = [0.0_f64; 3];
    for i in (0..3).rev() {
        let mut sum = b[i];
        for j in (i + 1)..3 {
            sum -= u[(i, j)] * x[j];
        }
        if u[(i, i)].abs() < EPS {
            return None;
        }
        x[i] = sum / u[(i, i)];
    }
    Some(Vec3::new(x[0], x[1], x[2]))
}
/// 3x3 SVD: M = U Sigma V^T.
///
/// Returns `(U, singular_values, V)` where `singular_values` are in
/// **descending** order. Both U and V are orthogonal matrices.
///
/// Algorithm: form B = M^T M (symmetric PSD), eigen-decompose B to get V and
/// sigma_i = sqrt(lambda_i), then compute U via u_i = M v_i / sigma_i.
pub fn svd3(m: Mat3) -> (Mat3, Vec3, Mat3) {
    let b = m.transpose() * m;
    let (eigenvalues, v) = symmetric_eigen3(b);
    let sv_asc = [
        eigenvalues.x.max(0.0).sqrt(),
        eigenvalues.y.max(0.0).sqrt(),
        eigenvalues.z.max(0.0).sqrt(),
    ];
    let mut idx = [0usize, 1, 2];
    idx.sort_by(|&a, &b| {
        sv_asc[b]
            .partial_cmp(&sv_asc[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let sigma = Vec3::new(sv_asc[idx[0]], sv_asc[idx[1]], sv_asc[idx[2]]);
    let v_col = |i: usize| Vec3::new(v[(0, i)], v[(1, i)], v[(2, i)]);
    let v_sorted = Mat3::from_columns(&[v_col(idx[0]), v_col(idx[1]), v_col(idx[2])]);
    let mut u_cols = [Vec3::zeros(); 3];
    let mut filled = [false; 3];
    for i in 0..3 {
        let vi = Vec3::new(v_sorted[(0, i)], v_sorted[(1, i)], v_sorted[(2, i)]);
        let si = sigma[i];
        if si > EPS {
            u_cols[i] = (m * vi) / si;
            filled[i] = true;
        }
    }
    let has_zero = filled.iter().any(|&f| !f);
    let u = if has_zero {
        let mut u_orth = u_cols;
        for j in 0..3 {
            if !filled[j] {
                u_orth[j] = fallback_orthogonal_filled(j, &u_orth, &filled);
            }
        }
        Mat3::from_columns(&u_orth)
    } else {
        Mat3::from_columns(&u_cols)
    };
    (u, sigma, v_sorted)
}
/// Estimate the condition number of a 3x3 matrix using the 1-norm.
///
/// kappa(M) = ||M||_1 * ||M^{-1}||_1
/// Returns `f64::INFINITY` if the matrix is singular.
pub fn condition_number_1(m: Mat3) -> Real {
    let n = norm1_3(m);
    match inv3(m) {
        Some(inv) => n * norm1_3(inv),
        None => f64::INFINITY,
    }
}
/// Estimate the condition number using the Frobenius norm.
///
/// kappa_F(M) = ||M||_F * ||M^{-1}||_F
/// Returns `f64::INFINITY` if the matrix is singular.
pub fn condition_number_frobenius(m: Mat3) -> Real {
    let n = frobenius_norm3(m);
    match inv3(m) {
        Some(inv) => n * frobenius_norm3(inv),
        None => f64::INFINITY,
    }
}
/// Estimate the 2-norm condition number using SVD.
///
/// kappa_2(M) = sigma_max / sigma_min
/// Returns `f64::INFINITY` if the matrix is singular (sigma_min ~ 0).
pub fn condition_number_2(m: Mat3) -> Real {
    let (_, sigma, _) = svd3(m);
    if sigma.z < EPS {
        f64::INFINITY
    } else {
        sigma.x / sigma.z
    }
}
/// Polar decomposition: M = R * S.
///
/// R is orthogonal (a rotation/reflection) and S is symmetric positive
/// semi-definite (stretch).  Uses SVD: M = U Sigma V^T -> R = U V^T, S = V Sigma V^T.
///
/// Used in corotational FEM and shape matching.
pub fn polar_decomp3(m: Mat3) -> (Mat3, Mat3) {
    let (u, sigma, v) = svd3(m);
    let r = u * v.transpose();
    let sigma_mat = Mat3::from_diagonal(&sigma);
    let s = v * sigma_mat * v.transpose();
    (r, s)
}
/// Matrix exponential of a symmetric 3x3 matrix via eigendecomposition.
///
/// exp(M) = V * diag(exp(lambda_i)) * V^T
#[allow(dead_code)]
pub fn symmetric_matrix_exp3(m: Mat3) -> Mat3 {
    let (evals, evecs) = symmetric_eigen3(m);
    let exp_evals = Vec3::new(evals.x.exp(), evals.y.exp(), evals.z.exp());
    let d = Mat3::from_diagonal(&exp_evals);
    evecs * d * evecs.transpose()
}
/// Matrix logarithm of a symmetric positive-definite 3x3 matrix.
///
/// log(M) = V * diag(ln(lambda_i)) * V^T
/// Returns `None` if any eigenvalue is non-positive.
#[allow(dead_code)]
pub fn symmetric_matrix_log3(m: Mat3) -> Option<Mat3> {
    let (evals, evecs) = symmetric_eigen3(m);
    if evals.x <= 0.0 || evals.y <= 0.0 || evals.z <= 0.0 {
        return None;
    }
    let log_evals = Vec3::new(evals.x.ln(), evals.y.ln(), evals.z.ln());
    let d = Mat3::from_diagonal(&log_evals);
    Some(evecs * d * evecs.transpose())
}
/// Matrix square root of a symmetric positive-semidefinite 3x3 matrix.
///
/// sqrt(M) = V * diag(sqrt(lambda_i)) * V^T
#[allow(dead_code)]
pub fn symmetric_matrix_sqrt3(m: Mat3) -> Mat3 {
    let (evals, evecs) = symmetric_eigen3(m);
    let sqrt_evals = Vec3::new(
        evals.x.max(0.0).sqrt(),
        evals.y.max(0.0).sqrt(),
        evals.z.max(0.0).sqrt(),
    );
    let d = Mat3::from_diagonal(&sqrt_evals);
    evecs * d * evecs.transpose()
}
/// LU factorisation with partial pivoting for an n×n matrix (stored row-major).
///
/// Returns `(lu, perm)` where `lu` is the in-place LU factorisation (L below
/// diagonal, U on and above) and `perm[i]` is the row permuted to position `i`.
/// Returns `None` if the matrix is (near-)singular.
pub fn lu_factor_n(a: &[Vec<f64>]) -> Option<(Vec<Vec<f64>>, Vec<usize>)> {
    let n = a.len();
    let mut lu: Vec<Vec<f64>> = a.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut max_val = lu[k][k].abs();
        let mut max_row = k;
        for i in (k + 1)..n {
            if lu[i][k].abs() > max_val {
                max_val = lu[i][k].abs();
                max_row = i;
            }
        }
        if max_val < EPS {
            return None;
        }
        if max_row != k {
            lu.swap(k, max_row);
            perm.swap(k, max_row);
        }
        for i in (k + 1)..n {
            let factor = lu[i][k] / lu[k][k];
            lu[i][k] = factor;
            for j in (k + 1)..n {
                let lkj = lu[k][j];
                lu[i][j] -= factor * lkj;
            }
        }
    }
    Some((lu, perm))
}
/// Forward substitution for n×n lower-triangular system L * x = b.
///
/// `lu` is the compact LU array (diagonal of L is implicitly 1).
pub fn forward_sub_n(lu: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = lu.len();
    let mut x = vec![0.0f64; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= lu[i][j] * x[j];
        }
        x[i] = s;
    }
    x
}
/// Backward substitution for n×n upper-triangular system U * x = b.
///
/// `lu` is the compact LU array; the diagonal and upper part hold U.
/// Returns `None` if a diagonal element is zero.
pub fn backward_sub_n(lu: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = lu.len();
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        if lu[i][i].abs() < EPS {
            return None;
        }
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= lu[i][j] * x[j];
        }
        x[i] = s / lu[i][i];
    }
    Some(x)
}
/// Solve A * x = b for general n×n system using LU with partial pivoting.
///
/// Returns `None` if the matrix is singular.
pub fn lu_solve_n(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    assert_eq!(b.len(), n, "lu_solve_n: dimension mismatch");
    let (lu, perm) = lu_factor_n(a)?;
    let pb: Vec<f64> = (0..n).map(|i| b[perm[i]]).collect();
    let y = forward_sub_n(&lu, &pb);
    backward_sub_n(&lu, &y)
}
/// QR factorisation of an m×n matrix (m ≥ n) via Householder reflections.
///
/// Returns `(q, r)` where `q` is m×m orthogonal and `r` is m×n upper triangular,
/// satisfying A = Q * R.
/// Stored as `Vec<Vec`f64`>` (row-major).
pub fn qr_factor_n(a: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let m = a.len();
    let n = if m > 0 { a[0].len() } else { 0 };
    assert!(m >= n, "qr_factor_n: requires m >= n");
    let mut r: Vec<Vec<f64>> = a.to_vec();
    let mut q: Vec<Vec<f64>> = (0..m)
        .map(|i| {
            let mut row = vec![0.0f64; m];
            row[i] = 1.0;
            row
        })
        .collect();
    for k in 0..n {
        let norm_x: f64 = (k..m).map(|i| r[i][k] * r[i][k]).sum::<f64>().sqrt();
        if norm_x < EPS {
            continue;
        }
        let sign = if r[k][k] >= 0.0 { 1.0 } else { -1.0 };
        let mut v = vec![0.0f64; m];
        for i in k..m {
            v[i] = r[i][k];
        }
        v[k] += sign * norm_x;
        let v_norm_sq: f64 = (k..m).map(|i| v[i] * v[i]).sum();
        if v_norm_sq < EPS {
            continue;
        }
        for j in k..n {
            let dot: f64 = (k..m).map(|i| v[i] * r[i][j]).sum();
            let factor = 2.0 * dot / v_norm_sq;
            for i in k..m {
                r[i][j] -= factor * v[i];
            }
        }
        for i in 0..m {
            let dot: f64 = (k..m).map(|l| q[i][l] * v[l]).sum();
            let factor = 2.0 * dot / v_norm_sq;
            for l in k..m {
                q[i][l] -= factor * v[l];
            }
        }
    }
    (q, r)
}
/// Solve the least-squares problem min ||A x - b||_2 via QR factorisation.
///
/// `a` is m×n (m ≥ n), `b` has length m.
/// Returns the length-n solution vector.
pub fn qr_solve_n(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let m = a.len();
    let n = if m > 0 { a[0].len() } else { 0 };
    assert_eq!(b.len(), m, "qr_solve_n: b length must equal number of rows");
    let (q, r) = qr_factor_n(a);
    let qtb: Vec<f64> = (0..m)
        .map(|i| (0..m).map(|j| q[j][i] * b[j]).sum())
        .collect();
    let r_sq: Vec<Vec<f64>> = (0..n).map(|i| r[i][..n].to_vec()).collect();
    let rhs: Vec<f64> = qtb[..n].to_vec();
    backward_sub_n(&r_sq, &rhs)
}
/// Frobenius norm of an n×n (or m×n) matrix stored as `Vec<Vec`f64`>`.
pub fn frobenius_norm_n(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .flat_map(|row| row.iter())
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
}
/// Matrix-vector multiply for a dense m×n matrix: y = A * x.
pub fn matvec_n(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(a, b)| a * b).sum())
        .collect()
}
/// Compute A^T * A for a dense m×n matrix.
pub fn ata_n(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let n = if m > 0 { a[0].len() } else { 0 };
    let mut ata = vec![vec![0.0f64; n]; n];
    for row in a {
        for i in 0..n {
            for j in 0..n {
                ata[i][j] += row[i] * row[j];
            }
        }
    }
    ata
}
/// Power iteration to find the dominant eigenvalue/eigenvector of an n×n matrix.
///
/// Returns `(eigenvalue, eigenvector)`. Converges to the largest-magnitude eigenvalue.
/// Returns `(0.0, zero_vec)` if the starting vector is zero or n=0.
#[allow(dead_code)]
pub fn power_iteration_n(a: &[Vec<f64>], max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let n = a.len();
    if n == 0 {
        return (0.0, Vec::new());
    }
    let mut v = vec![1.0f64; n];
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    for vi in v.iter_mut() {
        *vi /= norm;
    }
    let mut lambda = 0.0f64;
    for _ in 0..max_iter {
        let av = matvec_n(a, &v);
        let new_lambda: f64 = av.iter().zip(v.iter()).map(|(ai, vi)| ai * vi).sum();
        let norm: f64 = av.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < EPS {
            break;
        }
        let new_v: Vec<f64> = av.iter().map(|x| x / norm).collect();
        if (new_lambda - lambda).abs() < tol {
            lambda = new_lambda;
            v = new_v;
            break;
        }
        lambda = new_lambda;
        v = new_v;
    }
    (lambda, v)
}
/// Inverse iteration (shift-invert) to find eigenvalue closest to `shift`.
///
/// Solves `(A - shift*I) * v = v_old` at each step.
/// Returns `(eigenvalue, eigenvector)`.
#[allow(dead_code)]
pub fn inverse_iteration_n(
    a: &[Vec<f64>],
    shift: f64,
    max_iter: usize,
    tol: f64,
) -> (f64, Vec<f64>) {
    let n = a.len();
    if n == 0 {
        return (0.0, Vec::new());
    }
    let mut as_ = a.to_vec();
    for i in 0..n {
        as_[i][i] -= shift;
    }
    let mut v = vec![1.0f64; n];
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    for vi in v.iter_mut() {
        *vi /= norm;
    }
    let mut lambda = shift;
    for _ in 0..max_iter {
        let new_v_opt = lu_solve_n(&as_, &v);
        let new_v_raw = match new_v_opt {
            Some(x) => x,
            None => break,
        };
        let norm: f64 = new_v_raw.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < EPS {
            break;
        }
        let new_v: Vec<f64> = new_v_raw.iter().map(|x| x / norm).collect();
        let av = matvec_n(a, &new_v);
        let new_lambda: f64 = av.iter().zip(new_v.iter()).map(|(ai, vi)| ai * vi).sum();
        if (new_lambda - lambda).abs() < tol {
            lambda = new_lambda;
            v = new_v;
            break;
        }
        lambda = new_lambda;
        v = new_v;
    }
    (lambda, v)
}
/// Arnoldi iteration: builds an orthonormal Krylov basis Q and upper Hessenberg H.
///
/// Given matrix `a` (n×n) and starting vector `b`, runs `k` steps.
/// Returns `(Q, H)` where Q is n×k orthonormal and H is (k+1)×k upper Hessenberg.
#[allow(dead_code)]
pub fn arnoldi_n(a: &[Vec<f64>], b: &[f64], k: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    let k = k.min(n);
    let mut q: Vec<Vec<f64>> = vec![vec![0.0f64; k]; n];
    let mut h: Vec<Vec<f64>> = vec![vec![0.0f64; k]; k + 1];
    let norm: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm < EPS {
        return (q, h);
    }
    for i in 0..n {
        q[i][0] = b[i] / norm;
    }
    for j in 0..k {
        let qj: Vec<f64> = (0..n).map(|i| q[i][j]).collect();
        let mut z = matvec_n(a, &qj);
        for i in 0..=j {
            let qi: Vec<f64> = (0..n).map(|row| q[row][i]).collect();
            let dot: f64 = z.iter().zip(qi.iter()).map(|(zi, qi)| zi * qi).sum();
            h[i][j] = dot;
            for row in 0..n {
                z[row] -= dot * qi[row];
            }
        }
        let norm: f64 = z.iter().map(|x| x * x).sum::<f64>().sqrt();
        h[j + 1][j] = norm;
        if j + 1 < k && norm > EPS {
            for row in 0..n {
                q[row][j + 1] = z[row] / norm;
            }
        }
    }
    (q, h)
}
/// GMRES (Generalised Minimal Residual) solver for A x = b.
///
/// Restarts after `max_iter` iterations if not converged.
/// Returns `None` if the algorithm stagnates or fails.
#[allow(dead_code)]
pub fn gmres_n(a: &[Vec<f64>], b: &[f64], max_iter: usize, tol: f64) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut x = vec![0.0f64; n];
    let k = max_iter.min(n);
    let ax = matvec_n(a, &x);
    let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, axi)| bi - axi).collect();
    let beta: f64 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
    if beta < tol {
        return Some(x);
    }
    let mut q_cols: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
    q_cols.push(r.iter().map(|v| v / beta).collect());
    let mut h: Vec<Vec<f64>> = vec![vec![0.0f64; k]; k + 1];
    let mut cs = vec![0.0f64; k];
    let mut sn = vec![0.0f64; k];
    let mut e1 = vec![0.0f64; k + 1];
    e1[0] = beta;
    let mut j_final = k;
    for j in 0..k {
        let qj = q_cols[j].clone();
        let mut w = matvec_n(a, &qj);
        for i in 0..=j {
            let dot: f64 = w.iter().zip(q_cols[i].iter()).map(|(wi, qi)| wi * qi).sum();
            h[i][j] = dot;
            for row in 0..n {
                w[row] -= dot * q_cols[i][row];
            }
        }
        let norm: f64 = w.iter().map(|v| v * v).sum::<f64>().sqrt();
        h[j + 1][j] = norm;
        if norm > EPS {
            q_cols.push(w.iter().map(|v| v / norm).collect());
        } else {
            q_cols.push(vec![0.0f64; n]);
        }
        for i in 0..j {
            let t = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
            h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
            h[i][j] = t;
        }
        let r_val = (h[j][j] * h[j][j] + h[j + 1][j] * h[j + 1][j]).sqrt();
        if r_val < EPS {
            j_final = j;
            break;
        }
        cs[j] = h[j][j] / r_val;
        sn[j] = h[j + 1][j] / r_val;
        h[j][j] = r_val;
        h[j + 1][j] = 0.0;
        e1[j + 1] = -sn[j] * e1[j];
        e1[j] *= cs[j];
        if e1[j + 1].abs() < tol {
            j_final = j + 1;
            break;
        }
    }
    let m = j_final.min(k);
    if m == 0 {
        return Some(x);
    }
    let h_sq: Vec<Vec<f64>> = (0..m).map(|i| h[i][..m].to_vec()).collect();
    let rhs: Vec<f64> = e1[..m].to_vec();
    let y = backward_sub_n(&h_sq, &rhs)?;
    for i in 0..m {
        for row in 0..n {
            x[row] += y[i] * q_cols[i][row];
        }
    }
    Some(x)
}
/// Diagonal (Jacobi) preconditioner: returns the diagonal of `a`.
#[allow(dead_code)]
pub fn diagonal_preconditioner_n(a: &[Vec<f64>]) -> Vec<f64> {
    (0..a.len()).map(|i| a[i][i]).collect()
}
/// Apply diagonal preconditioner: z = D^{-1} r.
#[allow(dead_code)]
pub fn apply_diagonal_preconditioner(diag: &[f64], r: &[f64]) -> Vec<f64> {
    r.iter()
        .zip(diag.iter())
        .map(|(ri, di)| if di.abs() > EPS { ri / di } else { *ri })
        .collect()
}
/// ILU(0) incomplete LU factorisation preserving the sparsity pattern.
///
/// Returns `(L, U)` where L is strictly lower triangular with unit diagonal
/// and U is upper triangular, both stored as dense `Vec<Vec`f64`>`.
#[allow(dead_code)]
pub fn ilu0_n(a: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut l: Vec<Vec<f64>> = vec![vec![0.0f64; n]; n];
    let mut u: Vec<Vec<f64>> = a.to_vec();
    for i in 0..n {
        l[i][i] = 1.0;
    }
    for k in 0..n {
        for i in (k + 1)..n {
            if a[i][k].abs() < EPS {
                continue;
            }
            let factor = u[i][k] / u[k][k].max(EPS);
            l[i][k] = factor;
            for j in k..n {
                if a[i][j].abs() < EPS && j != k {
                    continue;
                }
                u[i][j] -= factor * u[k][j];
            }
            u[i][k] = 0.0;
        }
    }
    (l, u)
}
/// Real Schur decomposition of an n×n matrix via the QR algorithm.
///
/// Returns `(Q, T)` where Q is orthogonal and T is quasi-upper triangular
/// (upper triangular for matrices with only real eigenvalues),
/// satisfying `A = Q * T * Q^T`.
///
/// Uses the unshifted QR iteration: A_0 = A, Q_0 = I; A_{k+1} = R_k Q_k.
#[allow(dead_code)]
pub fn schur_n(a: &[Vec<f64>], max_iter: usize, tol: f64) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut t = a.to_vec();
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0f64; n];
            row[i] = 1.0;
            row
        })
        .collect();
    for _ in 0..max_iter {
        let off: f64 = (1..n)
            .flat_map(|i| (0..i).map(move |j| (i, j)))
            .map(|(i, j)| t[i][j].abs())
            .sum();
        if off < tol {
            break;
        }
        let (qk, rk) = qr_factor_n(&t);
        t = matmul_n_sq(&rk, &qk);
        q = matmul_n_sq(&q, &qk);
    }
    (q, t)
}
/// Dense n×n matrix multiplication (square).
pub(super) fn matmul_n_sq(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut c = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for k in 0..n {
            for j in 0..n {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
#[cfg(test)]
mod tests {
    use super::*;
    pub(super) const TOL: Real = 1e-9;
    fn mat_approx_eq(a: Mat3, b: Mat3, tol: Real) -> bool {
        (a - b).abs().max() < tol
    }
    fn vec_approx_eq(a: Vec3, b: Vec3, tol: Real) -> bool {
        (a - b).abs().max() < tol
    }
    #[test]
    fn test_det3_identity() {
        let d = det3(Mat3::identity());
        assert!((d - 1.0).abs() < TOL, "det(I) = {}", d);
    }
    #[test]
    fn test_det3_singular() {
        let m = Mat3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        let d = det3(m);
        assert!(d.abs() < 1e-10, "det of singular matrix = {}", d);
    }
    #[test]
    fn test_inv3_identity() {
        let inv = inv3(Mat3::identity()).expect("should invert I");
        assert!(mat_approx_eq(inv, Mat3::identity(), TOL));
    }
    #[test]
    fn test_inv3_singular() {
        let m = Mat3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        assert!(inv3(m).is_none());
    }
    #[test]
    fn test_inv3_roundtrip() {
        let m = Mat3::new(2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0);
        let inv = inv3(m).expect("invertible");
        let product = m * inv;
        assert!(
            mat_approx_eq(product, Mat3::identity(), 1e-10),
            "M * M^-1 not identity"
        );
    }
    #[test]
    fn test_solve3_basic() {
        let x = solve3(Mat3::identity(), Vec3::new(1.0, 2.0, 3.0)).expect("should solve");
        assert!(vec_approx_eq(x, Vec3::new(1.0, 2.0, 3.0), TOL));
    }
    #[test]
    fn test_solve3_singular() {
        let m = Mat3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        assert!(solve3(m, Vec3::new(1.0, 0.0, 0.0)).is_none());
    }
    #[test]
    fn test_eigen_diagonal() {
        let m = Mat3::from_diagonal(&Vec3::new(3.0, 1.0, 2.0));
        let (evals, _evecs) = symmetric_eigen3(m);
        assert!(
            (evals.x - 1.0).abs() < 1e-8
                && (evals.y - 2.0).abs() < 1e-8
                && (evals.z - 3.0).abs() < 1e-8,
            "eigenvalues: {:?}",
            evals
        );
    }
    #[test]
    fn test_eigen_symmetric() {
        let m = Mat3::new(4.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 2.0);
        let (evals, evecs) = symmetric_eigen3(m);
        for i in 0..3 {
            let vi = Vec3::new(evecs[(0, i)], evecs[(1, i)], evecs[(2, i)]);
            let mv = m * vi;
            let lv = vi * evals[i];
            let residual = (mv - lv).norm();
            assert!(
                residual < 1e-8,
                "eigenpair {}: M*v - lambda*v residual = {}",
                i,
                residual
            );
        }
    }
    #[test]
    fn test_eigen_symmetric_reconstruction() {
        let m = Mat3::new(4.0, 1.0, 0.5, 1.0, 3.0, 0.25, 0.5, 0.25, 2.0);
        let (evals, evecs) = symmetric_eigen3(m);
        let sigma = Mat3::from_diagonal(&evals);
        let reconstructed = evecs * sigma * evecs.transpose();
        assert!(
            mat_approx_eq(reconstructed, m, 1e-8),
            "reconstruction error: max = {}",
            (reconstructed - m).abs().max()
        );
    }
    #[test]
    fn test_polar_decomp_rotation_only() {
        use std::f64::consts::FRAC_PI_2;
        let (c, s) = (FRAC_PI_2.cos(), FRAC_PI_2.sin());
        let rot = Mat3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
        let (r, stretch) = polar_decomp3(rot);
        let r_diff = (r - rot).abs().max();
        assert!(r_diff < 1e-8, "R != M for pure rotation; diff = {}", r_diff);
        let s_diff = (stretch - Mat3::identity()).abs().max();
        assert!(s_diff < 1e-8, "S != I for pure rotation; diff = {}", s_diff);
    }
    #[test]
    fn test_polar_decomp_r_orthogonal() {
        let m = Mat3::new(2.0, 1.0, 0.0, 0.5, 3.0, 0.5, 0.0, 0.25, 1.5);
        let (r, _s) = polar_decomp3(m);
        let rt_r = r.transpose() * r;
        let diff = (rt_r - Mat3::identity()).abs().max();
        assert!(diff < 1e-8, "R^T * R != I; diff = {}", diff);
    }
    #[test]
    fn test_polar_decomp_reconstruction() {
        let m = Mat3::new(2.0, 1.0, 0.0, 0.5, 3.0, 0.5, 0.0, 0.25, 1.5);
        let (r, s) = polar_decomp3(m);
        let reconstructed = r * s;
        let diff = (reconstructed - m).abs().max();
        assert!(diff < 1e-8, "R * S != M; diff = {}", diff);
    }
    #[test]
    fn test_svd3_reconstruction() {
        let m = Mat3::new(1.0, 2.0, 0.0, 0.0, 3.0, 1.0, 1.0, 0.0, 4.0);
        let (u, sigma, v) = svd3(m);
        let sigma_mat = Mat3::from_diagonal(&sigma);
        let reconstructed = u * sigma_mat * v.transpose();
        let diff = (reconstructed - m).abs().max();
        assert!(diff < 1e-8, "U * Sigma * V^T != M; diff = {}", diff);
    }
    #[test]
    fn test_svd3_orthogonal() {
        let m = Mat3::new(1.0, 2.0, 0.0, 0.0, 3.0, 1.0, 1.0, 0.0, 4.0);
        let (u, _sigma, v) = svd3(m);
        let ut_u = u.transpose() * u;
        let diff_u = (ut_u - Mat3::identity()).abs().max();
        assert!(diff_u < 1e-8, "U not orthogonal; diff = {}", diff_u);
        let vt_v = v.transpose() * v;
        let diff_v = (vt_v - Mat3::identity()).abs().max();
        assert!(diff_v < 1e-8, "V not orthogonal; diff = {}", diff_v);
    }
    #[test]
    fn test_frobenius_identity() {
        let f = frobenius_norm3(Mat3::identity());
        let expected = (3.0_f64).sqrt();
        assert!(
            (f - expected).abs() < TOL,
            "||I||_F = {}, expected {}",
            f,
            expected
        );
    }
    #[test]
    fn test_norm1_identity() {
        let n = norm1_3(Mat3::identity());
        assert!((n - 1.0).abs() < TOL, "||I||_1 = {}, expected 1", n);
    }
    #[test]
    fn test_norm_inf_identity() {
        let n = norm_inf3(Mat3::identity());
        assert!((n - 1.0).abs() < TOL, "||I||_inf = {}, expected 1", n);
    }
    #[test]
    fn test_norm1_known() {
        let m = Mat3::new(1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0);
        let n = norm1_3(m);
        assert!((n - 18.0).abs() < TOL, "||M||_1 = {}", n);
    }
    #[test]
    fn test_norm_inf_known() {
        let m = Mat3::new(1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0);
        let n = norm_inf3(m);
        assert!((n - 24.0).abs() < TOL, "||M||_inf = {}", n);
    }
    #[test]
    fn test_trace3() {
        let m = Mat3::new(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0);
        assert!((trace3(m) - 6.0).abs() < TOL);
    }
    #[test]
    fn test_is_symmetric3() {
        let sym = Mat3::new(1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0);
        assert!(is_symmetric3(sym, 1e-10));
        let asym = Mat3::new(1.0, 2.0, 3.0, 0.0, 4.0, 5.0, 3.0, 5.0, 6.0);
        assert!(!is_symmetric3(asym, 1e-10));
    }
    #[test]
    fn test_is_orthogonal3() {
        assert!(is_orthogonal3(Mat3::identity(), 1e-10));
        let m = Mat3::new(2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        assert!(!is_orthogonal3(m, 1e-10));
    }
    #[test]
    fn test_characteristic_poly3_identity() {
        let (c0, c1, c2) = characteristic_poly3(Mat3::identity());
        assert!((c0 - 1.0).abs() < TOL, "c0 = {}", c0);
        assert!((c1 - 3.0).abs() < TOL, "c1 = {}", c1);
        assert!((c2 - 3.0).abs() < TOL, "c2 = {}", c2);
    }
    #[test]
    fn test_symmetric_eigenvalues3_diagonal() {
        let m = Mat3::from_diagonal(&Vec3::new(5.0, 2.0, 8.0));
        let evals = symmetric_eigenvalues3(m);
        assert!(
            (evals.x - 2.0).abs() < 1e-8
                && (evals.y - 5.0).abs() < 1e-8
                && (evals.z - 8.0).abs() < 1e-8,
            "evals = {:?}",
            evals
        );
    }
    #[test]
    fn test_lu_decomp3_identity() {
        let (l, u, perm) = lu_decomp3(Mat3::identity()).expect("should decompose identity");
        assert!(mat_approx_eq(l, Mat3::identity(), TOL));
        assert!(mat_approx_eq(u, Mat3::identity(), TOL));
        assert_eq!(perm, [0, 1, 2]);
    }
    #[test]
    fn test_lu_decomp3_reconstruction() {
        let m = Mat3::new(2.0, 1.0, 0.0, 4.0, 3.0, 1.0, 0.0, 1.0, 4.0);
        let (l, u, perm) = lu_decomp3(m).expect("should decompose");
        let lu = l * u;
        let pm = Mat3::new(
            m[(perm[0], 0)],
            m[(perm[0], 1)],
            m[(perm[0], 2)],
            m[(perm[1], 0)],
            m[(perm[1], 1)],
            m[(perm[1], 2)],
            m[(perm[2], 0)],
            m[(perm[2], 1)],
            m[(perm[2], 2)],
        );
        assert!(mat_approx_eq(lu, pm, 1e-10), "L*U != P*M");
    }
    #[test]
    fn test_lu_decomp3_singular() {
        let m = Mat3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        assert!(lu_decomp3(m).is_none());
    }
    #[test]
    fn test_lu_solve3() {
        let m = Mat3::new(2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0);
        let b = Vec3::new(1.0, 2.0, 3.0);
        let x = lu_solve3(m, b).expect("should solve");
        let residual = (m * x - b).norm();
        assert!(residual < 1e-10, "LU solve residual = {}", residual);
    }
    #[test]
    fn test_qr_householder_reconstruction() {
        let m = Mat3::new(2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0);
        let (q, r) = qr_householder3(m);
        let reconstructed = q * r;
        assert!(
            mat_approx_eq(reconstructed, m, 1e-8),
            "Q*R != M; max diff = {}",
            (reconstructed - m).abs().max()
        );
    }
    #[test]
    fn test_qr_householder_orthogonal() {
        let m = Mat3::new(2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0);
        let (q, _r) = qr_householder3(m);
        let qt_q = q.transpose() * q;
        assert!(
            mat_approx_eq(qt_q, Mat3::identity(), 1e-8),
            "Q not orthogonal; max diff = {}",
            (qt_q - Mat3::identity()).abs().max()
        );
    }
    #[test]
    fn test_qr_householder_r_upper_triangular() {
        let m = Mat3::new(2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 4.0);
        let (_q, r) = qr_householder3(m);
        assert!(r[(1, 0)].abs() < 1e-10, "R[1,0] = {}", r[(1, 0)]);
        assert!(r[(2, 0)].abs() < 1e-10, "R[2,0] = {}", r[(2, 0)]);
        assert!(r[(2, 1)].abs() < 1e-10, "R[2,1] = {}", r[(2, 1)]);
    }
    #[test]
    fn test_forward_substitute() {
        let l = Mat3::new(1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 3.0, 4.0, 1.0);
        let b = Vec3::new(1.0, 4.0, 15.0);
        let x = forward_substitute3(l, b).expect("should solve");
        let residual = (l * x - b).norm();
        assert!(
            residual < 1e-10,
            "forward substitution residual = {}",
            residual
        );
    }
    #[test]
    fn test_backward_substitute() {
        let u = Mat3::new(2.0, 1.0, 3.0, 0.0, 4.0, 5.0, 0.0, 0.0, 6.0);
        let b = Vec3::new(14.0, 29.0, 18.0);
        let x = backward_substitute3(u, b).expect("should solve");
        let residual = (u * x - b).norm();
        assert!(
            residual < 1e-10,
            "backward substitution residual = {}",
            residual
        );
    }
    #[test]
    fn test_condition_number_identity() {
        let k1 = condition_number_1(Mat3::identity());
        assert!((k1 - 1.0).abs() < 1e-8, "kappa_1(I) = {}", k1);
    }
    #[test]
    fn test_condition_number_singular() {
        let m = Mat3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        assert!(condition_number_1(m).is_infinite());
        assert!(condition_number_frobenius(m).is_infinite());
        let k2 = condition_number_2(m);
        assert!(
            k2 > 1e6 || k2.is_infinite(),
            "kappa_2 should be huge or inf, got {}",
            k2
        );
    }
    #[test]
    fn test_condition_number_2_well_conditioned() {
        let m = Mat3::from_diagonal(&Vec3::new(1.0, 2.0, 3.0));
        let k2 = condition_number_2(m);
        assert!((k2 - 3.0).abs() < 1e-6, "kappa_2 = {}", k2);
    }
    #[test]
    fn test_condition_number_frobenius() {
        let k = condition_number_frobenius(Mat3::identity());
        assert!((k - 3.0).abs() < 1e-8, "kappa_F(I) = {}", k);
    }
    #[test]
    fn test_symmetric_matrix_exp_identity() {
        let result = symmetric_matrix_exp3(Mat3::zeros());
        assert!(mat_approx_eq(result, Mat3::identity(), 1e-8), "exp(0) != I");
    }
    #[test]
    fn test_symmetric_matrix_log_identity() {
        let result = symmetric_matrix_log3(Mat3::identity()).expect("should compute log(I)");
        assert!(mat_approx_eq(result, Mat3::zeros(), 1e-8), "log(I) != 0");
    }
    #[test]
    fn test_symmetric_matrix_sqrt() {
        let m = Mat3::from_diagonal(&Vec3::new(4.0, 9.0, 16.0));
        let sq = symmetric_matrix_sqrt3(m);
        let expected = Mat3::from_diagonal(&Vec3::new(2.0, 3.0, 4.0));
        assert!(
            mat_approx_eq(sq, expected, 1e-8),
            "sqrt(diag(4,9,16)) != diag(2,3,4)"
        );
    }
    #[test]
    fn test_symmetric_matrix_exp_log_roundtrip() {
        let m = Mat3::from_diagonal(&Vec3::new(1.0, 2.0, 3.0));
        let exp_m = symmetric_matrix_exp3(m);
        let log_exp_m = symmetric_matrix_log3(exp_m).expect("should compute log(exp(M))");
        assert!(mat_approx_eq(log_exp_m, m, 1e-6), "log(exp(M)) != M");
    }
    #[test]
    fn test_qr_gram_schmidt_reconstruction() {
        let m = Mat3::new(1.0, 2.0, 0.0, 0.0, 3.0, 1.0, 1.0, 0.0, 4.0);
        let (q, r) = qr_decomp3(m);
        let reconstructed = q * r;
        assert!(
            mat_approx_eq(reconstructed, m, 1e-8),
            "QR reconstruction failed"
        );
    }
    #[test]
    fn test_qr_gram_schmidt_orthogonal() {
        let m = Mat3::new(1.0, 2.0, 0.0, 0.0, 3.0, 1.0, 1.0, 0.0, 4.0);
        let (q, _r) = qr_decomp3(m);
        let qt_q = q.transpose() * q;
        assert!(
            mat_approx_eq(qt_q, Mat3::identity(), 1e-8),
            "Q not orthogonal"
        );
    }
    fn make_4x4() -> Vec<Vec<f64>> {
        vec![
            vec![4.0, 3.0, 2.0, 1.0],
            vec![3.0, 4.0, 3.0, 2.0],
            vec![2.0, 3.0, 4.0, 3.0],
            vec![1.0, 2.0, 3.0, 4.0],
        ]
    }
    fn matmul_n(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = a.len();
        let p = b[0].len();
        let mut c = vec![vec![0.0f64; p]; n];
        for i in 0..n {
            for k in 0..a[i].len() {
                for j in 0..p {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    }
    #[test]
    fn test_lu_factor_n_identity() {
        let id: Vec<Vec<f64>> = (0..4)
            .map(|i| {
                let mut row = vec![0.0f64; 4];
                row[i] = 1.0;
                row
            })
            .collect();
        let (lu, perm) = lu_factor_n(&id).expect("should factor identity");
        assert_eq!(
            perm,
            vec![0, 1, 2, 3],
            "identity permutation should be trivial"
        );
        for i in 0..4 {
            assert!((lu[i][i] - 1.0).abs() < 1e-10, "LU diagonal[{i}] != 1");
        }
    }
    #[test]
    fn test_lu_factor_n_singular() {
        let singular: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
        ];
        assert!(
            lu_factor_n(&singular).is_none(),
            "singular matrix should return None"
        );
    }
    #[test]
    fn test_lu_solve_n_basic() {
        let a = make_4x4();
        let b = vec![10.0, 12.0, 12.0, 10.0];
        let x = lu_solve_n(&a, &b).expect("should solve");
        let ax = matvec_n(&a, &x);
        for i in 0..4 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-8,
                "LU solve 4×4: residual at [{i}]: {} vs {}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_lu_solve_n_identity() {
        let id: Vec<Vec<f64>> = (0..5)
            .map(|i| {
                let mut row = vec![0.0f64; 5];
                row[i] = 1.0;
                row
            })
            .collect();
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let x = lu_solve_n(&id, &b).expect("should solve I*x=b");
        for i in 0..5 {
            assert!(
                (x[i] - b[i]).abs() < 1e-10,
                "I*x=b should give x=b at [{i}]"
            );
        }
    }
    #[test]
    fn test_lu_solve_n_singular() {
        let a: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let b = vec![1.0, 2.0, 3.0];
        assert!(
            lu_solve_n(&a, &b).is_none(),
            "singular system should return None"
        );
    }
    #[test]
    fn test_forward_sub_n_identity_l() {
        let l: Vec<Vec<f64>> = (0..4)
            .map(|i| {
                let mut row = vec![0.0f64; 4];
                row[i] = 1.0;
                row
            })
            .collect();
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x = forward_sub_n(&l, &b);
        for i in 0..4 {
            assert!(
                (x[i] - b[i]).abs() < 1e-12,
                "forward sub with L=I: x[{i}]={} expected {}",
                x[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_backward_sub_n_basic() {
        let u: Vec<Vec<f64>> = vec![vec![2.0, 1.0], vec![0.0, 3.0]];
        let b = vec![5.0, 6.0];
        let x = backward_sub_n(&u, &b).expect("should solve");
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1] = {} expected 2", x[1]);
        assert!((x[0] - 1.5).abs() < 1e-10, "x[0] = {} expected 1.5", x[0]);
    }
    #[test]
    fn test_qr_factor_n_reconstruction_4x4() {
        let a = make_4x4();
        let (q, r) = qr_factor_n(&a);
        let qr = matmul_n(&q, &r);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (qr[i][j] - a[i][j]).abs() < 1e-8,
                    "QR reconstruction failed at [{i}][{j}]: {} vs {}",
                    qr[i][j],
                    a[i][j]
                );
            }
        }
    }
    #[test]
    fn test_qr_factor_n_q_orthogonal() {
        let a = make_4x4();
        let (q, _r) = qr_factor_n(&a);
        let n = q.len();
        for i in 0..n {
            for j in 0..n {
                let dot: f64 = (0..n).map(|k| q[k][i] * q[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-8,
                    "Q not orthogonal at ({i},{j}): {} vs {}",
                    dot,
                    expected
                );
            }
        }
    }
    #[test]
    fn test_qr_factor_n_r_upper_triangular() {
        let a = make_4x4();
        let (_q, r) = qr_factor_n(&a);
        for i in 0..4 {
            for j in 0..i {
                assert!(
                    r[i][j].abs() < 1e-8,
                    "R should be upper triangular: R[{i}][{j}] = {}",
                    r[i][j]
                );
            }
        }
    }
    #[test]
    fn test_qr_solve_n_square() {
        let a = make_4x4();
        let b = vec![10.0, 12.0, 12.0, 10.0];
        let x = qr_solve_n(&a, &b).expect("should solve");
        let ax = matvec_n(&a, &x);
        for i in 0..4 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-6,
                "QR solve 4×4: residual at [{i}]: {} vs {}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_qr_solve_n_overdetermined() {
        let a: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
            vec![0.0, 1.0, 1.0],
        ];
        let b = vec![1.0, 2.0, 3.0, 3.5, 4.5];
        let x = qr_solve_n(&a, &b).expect("should solve least squares");
        assert_eq!(x.len(), 3, "QR solve should return 3-vector");
        let ax = matvec_n(&a, &x);
        let res: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(ai, bi)| (ai - bi).powi(2))
            .sum();
        assert!(res < 1.0, "QR LS residual should be small: {res}");
    }
    #[test]
    fn test_frobenius_norm_n_identity() {
        let id: Vec<Vec<f64>> = (0..3)
            .map(|i| {
                let mut row = vec![0.0f64; 3];
                row[i] = 1.0;
                row
            })
            .collect();
        let f = frobenius_norm_n(&id);
        let expected = (3.0_f64).sqrt();
        assert!(
            (f - expected).abs() < 1e-10,
            "||I_3||_F = {} expected sqrt(3)",
            f
        );
    }
    #[test]
    fn test_matvec_n_identity() {
        let id: Vec<Vec<f64>> = (0..4)
            .map(|i| {
                let mut row = vec![0.0f64; 4];
                row[i] = 1.0;
                row
            })
            .collect();
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = matvec_n(&id, &x);
        for i in 0..4 {
            assert!(
                (y[i] - x[i]).abs() < 1e-12,
                "matvec_n with I: y[{i}] != x[{i}]"
            );
        }
    }
    #[test]
    fn test_ata_n_identity() {
        let id: Vec<Vec<f64>> = (0..3)
            .map(|i| {
                let mut row = vec![0.0f64; 3];
                row[i] = 1.0;
                row
            })
            .collect();
        let ata = ata_n(&id);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (ata[i][j] - expected).abs() < 1e-10,
                    "I^T * I should be I at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    fn test_ata_n_column_vector() {
        let a: Vec<Vec<f64>> = vec![vec![1.0], vec![2.0], vec![3.0]];
        let ata = ata_n(&a);
        assert!(
            (ata[0][0] - 14.0).abs() < 1e-10,
            "A^T A for column vector = {}",
            ata[0][0]
        );
    }
    #[test]
    fn test_polar_decomp3_stretch_positive_semidefinite() {
        let m = Mat3::new(3.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 2.0);
        let (_r, s) = polar_decomp3(m);
        let diff = (s - s.transpose()).abs().max();
        assert!(diff < 1e-8, "S should be symmetric; max diff = {diff}");
        let evals = symmetric_eigenvalues3(s);
        assert!(
            evals.x >= -1e-8,
            "S eigenvalue x should be >= 0: {}",
            evals.x
        );
        assert!(
            evals.y >= -1e-8,
            "S eigenvalue y should be >= 0: {}",
            evals.y
        );
        assert!(
            evals.z >= -1e-8,
            "S eigenvalue z should be >= 0: {}",
            evals.z
        );
    }
    #[test]
    fn test_polar_decomp3_det_r_is_pm1() {
        let m = Mat3::new(2.0, 0.5, 0.0, 0.5, 2.0, 0.0, 0.0, 0.0, 1.5);
        let (r, _s) = polar_decomp3(m);
        let d = det3(r);
        assert!((d.abs() - 1.0).abs() < 1e-8, "det(R) should be ±1, got {d}");
    }
    #[test]
    fn test_symmetric_eigenvalues3_stress_hydrostatic() {
        let p = 5.0_f64;
        let sigma = Mat3::from_diagonal(&Vec3::new(p, p, p));
        let evals = symmetric_eigenvalues3(sigma);
        assert!(
            (evals.x - p).abs() < 1e-8,
            "hydrostatic eigenvalue x = {}",
            evals.x
        );
        assert!(
            (evals.y - p).abs() < 1e-8,
            "hydrostatic eigenvalue y = {}",
            evals.y
        );
        assert!(
            (evals.z - p).abs() < 1e-8,
            "hydrostatic eigenvalue z = {}",
            evals.z
        );
    }
    #[test]
    fn test_symmetric_eigenvalues3_sum_equals_trace() {
        let m = Mat3::new(3.0, 1.0, 0.5, 1.0, 4.0, 0.2, 0.5, 0.2, 2.0);
        let evals = symmetric_eigenvalues3(m);
        let sum = evals.x + evals.y + evals.z;
        let tr = trace3(m);
        assert!(
            (sum - tr).abs() < 1e-6,
            "sum of eigenvalues should equal trace: {sum} vs {tr}"
        );
    }
    #[test]
    fn test_symmetric_eigenvalues3_product_equals_det() {
        let m = Mat3::new(3.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 2.0);
        let (evals, _) = symmetric_eigen3(m);
        let product = evals.x * evals.y * evals.z;
        let d = det3(m);
        assert!(
            (product - d).abs() < 1e-6,
            "product of eigenvalues should equal det: {product} vs {d}"
        );
    }
    #[test]
    fn test_frobenius_norm3_scaling() {
        let m = Mat3::identity();
        let scaled_m = m * 3.0;
        let f = frobenius_norm3(scaled_m);
        let expected = 3.0 * (3.0_f64).sqrt();
        assert!(
            (f - expected).abs() < 1e-8,
            "||3*I||_F = {f}, expected {expected}"
        );
    }
    #[test]
    fn test_norm1_3_column_dominant() {
        let m = Mat3::new(1.0, 0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 0.0, 2.0);
        let n = norm1_3(m);
        assert!((n - 6.0).abs() < 1e-10, "1-norm should be 6, got {n}");
    }
    #[test]
    fn test_norm_inf3_row_dominant() {
        let m = Mat3::new(1.0, 0.0, 0.0, 2.0, 2.0, 0.0, 3.0, 3.0, 3.0);
        let n = norm_inf3(m);
        assert!((n - 9.0).abs() < 1e-10, "inf-norm should be 9, got {n}");
    }
    #[test]
    fn test_power_iteration_diagonal() {
        let a = vec![
            vec![5.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let (lambda, _) = power_iteration_n(&a, 200, 1e-10);
        assert!((lambda - 5.0).abs() < 1e-6, "dominant eigenvalue={lambda}");
    }
    #[test]
    fn test_power_iteration_2x2() {
        let a = vec![vec![3.0, 1.0], vec![1.0, 3.0]];
        let (lambda, _) = power_iteration_n(&a, 200, 1e-10);
        assert!((lambda - 4.0).abs() < 1e-6, "dominant eigenvalue={lambda}");
    }
    #[test]
    fn test_inverse_iteration_smallest_eigenvalue() {
        let a = vec![
            vec![5.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let (lambda, _) = inverse_iteration_n(&a, 0.0, 200, 1e-10);
        assert!((lambda - 1.0).abs() < 1e-5, "smallest eigenvalue={lambda}");
    }
    #[test]
    fn test_inverse_iteration_shift() {
        let a = vec![
            vec![5.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let (lambda, _) = inverse_iteration_n(&a, 2.1, 200, 1e-8);
        assert!(
            (lambda - 2.0).abs() < 1e-4,
            "eigenvalue near 2, got {lambda}"
        );
    }
    #[test]
    fn test_arnoldi_orthonormality() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let b = vec![1.0, 0.0, 0.0];
        let (q, _h) = arnoldi_n(&a, &b, 3);
        let n = q.len();
        for i in 0..n {
            for j in 0..n {
                let dot: f64 = (0..a.len()).map(|k| q[k][i] * q[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-8,
                    "Arnoldi Q not orthonormal at ({i},{j}): {dot}"
                );
            }
        }
    }
    #[test]
    fn test_arnoldi_krylov_relation() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let b = vec![1.0, 0.0, 0.0];
        let (q, h) = arnoldi_n(&a, &b, 2);
        let m = a.len();
        let aq0: Vec<f64> = (0..m)
            .map(|i| a[i].iter().zip(q.iter()).map(|(aij, qj)| aij * qj[0]).sum())
            .collect();
        let rhs: Vec<f64> = (0..m)
            .map(|i| h[0][0] * q[i][0] + h[1][0] * q[i][1])
            .collect();
        for i in 0..m {
            assert!(
                (aq0[i] - rhs[i]).abs() < 1e-8,
                "Krylov relation at {i}: {} vs {}",
                aq0[i],
                rhs[i]
            );
        }
    }
    #[test]
    fn test_gmres_identity_system() {
        let a: Vec<Vec<f64>> = (0..3)
            .map(|i| {
                let mut row = vec![0.0; 3];
                row[i] = 1.0;
                row
            })
            .collect();
        let b = vec![1.0, 2.0, 3.0];
        let x = gmres_n(&a, &b, 50, 1e-10).expect("GMRES on identity should solve");
        for i in 0..3 {
            assert!(
                (x[i] - b[i]).abs() < 1e-8,
                "GMRES identity: x[{i}]={} expected {}",
                x[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_gmres_spd_system() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let b = vec![1.0, 2.0, 3.0];
        let x = gmres_n(&a, &b, 50, 1e-10).expect("GMRES on SPD");
        let ax = matvec_n(&a, &x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-6,
                "GMRES residual at {i}: {} vs {}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_diagonal_preconditioner() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let d = diagonal_preconditioner_n(&a);
        assert!((d[0] - 4.0).abs() < 1e-12, "d[0]={}", d[0]);
        assert!((d[1] - 3.0).abs() < 1e-12, "d[1]={}", d[1]);
        assert!((d[2] - 2.0).abs() < 1e-12, "d[2]={}", d[2]);
    }
    #[test]
    fn test_diagonal_preconditioning_apply() {
        let diag = vec![2.0, 4.0, 8.0];
        let r = vec![4.0, 8.0, 16.0];
        let z = apply_diagonal_preconditioner(&diag, &r);
        assert!((z[0] - 2.0).abs() < 1e-12);
        assert!((z[1] - 2.0).abs() < 1e-12);
        assert!((z[2] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_ilu0_factorization() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let (l, u) = ilu0_n(&a);
        for i in 0..3 {
            assert!((l[i][i] - 1.0).abs() < 1e-12, "L diagonal = 1");
            for j in (i + 1)..3 {
                assert!(l[i][j].abs() < 1e-12, "L upper triangular part should be 0");
            }
        }
        for i in 0..3 {
            for j in 0..i {
                assert!(u[i][j].abs() < 1e-12, "U lower triangular part should be 0");
            }
        }
    }
    #[test]
    fn test_ilu0_approximate_factorization() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let (l, u) = ilu0_n(&a);
        let lu = matmul_n(&l, &u);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (lu[i][j] - a[i][j]).abs() < 1e-8,
                    "ILU0: L*U[{i}][{j}]={} expected {}",
                    lu[i][j],
                    a[i][j]
                );
            }
        }
    }
    #[test]
    fn test_schur_decomp_diagonal_matrix() {
        let a = vec![
            vec![3.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let (q, t) = schur_n(&a, 100, 1e-10);
        for i in 1..3 {
            for j in 0..i {
                assert!(
                    t[i][j].abs() < 1e-6,
                    "Schur T[{i}][{j}]={} should be near 0",
                    t[i][j]
                );
            }
        }
        let n = q.len();
        for i in 0..n {
            for j in 0..n {
                let dot: f64 = (0..n).map(|k| q[k][i] * q[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-6,
                    "Q not orthogonal at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_schur_decomp_reconstruction() {
        let a = vec![
            vec![2.0, 1.0, 0.0],
            vec![0.0, 3.0, 1.0],
            vec![0.0, 0.0, 4.0],
        ];
        let (q, t) = schur_n(&a, 100, 1e-10);
        let qt: Vec<Vec<f64>> = (0..3).map(|i| (0..3).map(|j| q[j][i]).collect()).collect();
        let qt_a = matmul_n(&qt, &a);
        let qt_a_q = matmul_n(&qt_a, &q);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (qt_a_q[i][j] - t[i][j]).abs() < 1e-6,
                    "Q^T A Q != T at [{i}][{j}]: {} vs {}",
                    qt_a_q[i][j],
                    t[i][j]
                );
            }
        }
    }
    #[test]
    fn test_frobenius_norm_n_basic() {
        let a = vec![vec![3.0, 4.0], vec![0.0, 0.0]];
        let n = frobenius_norm_n(&a);
        assert!((n - 5.0).abs() < 1e-10, "||[3,4;0,0]||_F = 5, got {n}");
    }
    #[test]
    fn test_matvec_n_basic() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let x = vec![1.0, 1.0];
        let y = matvec_n(&a, &x);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_ata_n_basic() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 2.0]];
        let ata = ata_n(&a);
        assert!((ata[0][0] - 1.0).abs() < 1e-12);
        assert!((ata[1][1] - 4.0).abs() < 1e-12);
        assert!(ata[0][1].abs() < 1e-12);
    }
}
