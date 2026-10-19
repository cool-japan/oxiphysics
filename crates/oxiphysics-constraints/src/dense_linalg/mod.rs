//! Small dense linear-algebra kernels for the control and robotics modules.
//!
//! Matrices are row-major `f64` slices.  The kernels target the tiny systems
//! (a handful to a few hundred unknowns) met in controller design and robot
//! kinematics, where robustness and accuracy matter far more than raw speed:
//!
//! * [`cholesky_factor`] / [`cholesky_solve`] for symmetric positive definite
//!   systems (they double as a positive-definiteness test),
//! * [`lu_solve`] for general square systems (partial pivoting),
//! * [`singular_values`] by the one-sided Jacobi method, which keeps full
//!   relative accuracy for the tiny singular values of a nearly singular matrix,
//! * [`real_eigenvalues`] for the eigenvalues of a general real matrix
//!   (balancing, Householder reduction to Hessenberg form and the Francis
//!   double-shift QR iteration).

mod eigen;
mod svd;

#[cfg(test)]
mod tests;

pub(crate) use eigen::real_eigenvalues;
pub(crate) use svd::singular_values;

/// Relative pivot threshold below which a Cholesky pivot counts as not positive.
const CHOLESKY_PIVOT_TOL: f64 = 1e-13;

/// Lower-triangular Cholesky factor `L` (row-major `n x n`) with `A = L Lᵀ`.
///
/// Only the lower triangle of `a` is read.  Returns `None` when `a` is not
/// (numerically) positive definite: a pivot that is not larger than
/// `1e-13 * a[i][i]` — which includes zero, negative and non-finite pivots —
/// is rejected, so the call doubles as a positive-definiteness test.
pub(crate) fn cholesky_factor(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if !s.is_finite() || s <= CHOLESKY_PIVOT_TOL * a[i * n + i].abs() {
                    return None;
                }
                l[i * n + i] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    Some(l)
}

/// Solve `L Lᵀ X = B` for `X` given the factor from [`cholesky_factor`].
///
/// `b` is a row-major `n x nrhs` block of right-hand sides and the result has
/// the same layout.
pub(crate) fn cholesky_solve(l: &[f64], n: usize, b: &[f64], nrhs: usize) -> Vec<f64> {
    let mut x = b.to_vec();
    for c in 0..nrhs {
        for i in 0..n {
            let mut s = x[i * nrhs + c];
            for k in 0..i {
                s -= l[i * n + k] * x[k * nrhs + c];
            }
            x[i * nrhs + c] = s / l[i * n + i];
        }
        for i in (0..n).rev() {
            let mut s = x[i * nrhs + c];
            for k in (i + 1)..n {
                s -= l[k * n + i] * x[k * nrhs + c];
            }
            x[i * nrhs + c] = s / l[i * n + i];
        }
    }
    x
}

/// Solve `A X = B` for a general square `A` (row-major `n x n`) by Gaussian
/// elimination with scaled partial pivoting.
///
/// `b` is a row-major `n x nrhs` block of right-hand sides and the result has
/// the same layout.  Each row is judged against its own largest entry, so a
/// stiff matrix whose rows differ by many orders of magnitude is not mistaken
/// for a singular one.  Returns `None` when `A` is singular to working
/// precision (a scaled pivot not larger than `n · ε`), has a zero row, or
/// contains non-finite entries.
pub(crate) fn lu_solve(a: &[f64], n: usize, b: &[f64], nrhs: usize) -> Option<Vec<f64>> {
    let mut row_scale: Vec<f64> = (0..n)
        .map(|i| (0..n).fold(0.0_f64, |m, j| m.max(a[i * n + j].abs())))
        .collect();
    if row_scale.iter().any(|s| !s.is_finite() || *s == 0.0) {
        return None;
    }
    let tol = n as f64 * f64::EPSILON;
    let mut lu = a.to_vec();
    let mut x = b.to_vec();
    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_scaled = lu[col * n + col].abs() / row_scale[col];
        for row in (col + 1)..n {
            let v = lu[row * n + col].abs() / row_scale[row];
            if v > pivot_scaled {
                pivot_scaled = v;
                pivot_row = row;
            }
        }
        if pivot_scaled.is_nan() || pivot_scaled <= tol {
            return None;
        }
        if pivot_row != col {
            for k in 0..n {
                lu.swap(col * n + k, pivot_row * n + k);
            }
            for c in 0..nrhs {
                x.swap(col * nrhs + c, pivot_row * nrhs + c);
            }
            row_scale.swap(col, pivot_row);
        }
        let diag = lu[col * n + col];
        for row in (col + 1)..n {
            let factor = lu[row * n + col] / diag;
            if factor == 0.0 {
                continue;
            }
            for k in (col + 1)..n {
                lu[row * n + k] -= factor * lu[col * n + k];
            }
            for c in 0..nrhs {
                x[row * nrhs + c] -= factor * x[col * nrhs + c];
            }
        }
    }
    for c in 0..nrhs {
        for i in (0..n).rev() {
            let mut s = x[i * nrhs + c];
            for k in (i + 1)..n {
                s -= lu[i * n + k] * x[k * nrhs + c];
            }
            x[i * nrhs + c] = s / lu[i * n + i];
        }
    }
    Some(x)
}
