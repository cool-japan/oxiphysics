//! Singular values by the one-sided (Hestenes) Jacobi method.

/// Maximum number of sweeps over all column pairs.
const MAX_SWEEPS: usize = 60;

/// Columns whose normalised inner product is below this are treated as orthogonal.
const ORTHOGONALITY_TOL: f64 = 1e-14;

/// Singular values of the row-major `rows x cols` matrix `a`, in descending
/// order (`min(rows, cols)` values).
///
/// The one-sided Jacobi iteration rotates pairs of columns until all columns are
/// mutually orthogonal; the singular values are then the column norms.  Unlike
/// the eigenvalues of `AᵀA`, which square the condition number and lose every
/// singular value below `√ε · σ_max`, the method resolves tiny singular values
/// with high relative accuracy, so the product of the singular values of a
/// rank-deficient matrix vanishes to rounding level.  A wide matrix is
/// transposed first (the iteration needs at least as many rows as columns).
/// Non-finite input yields non-finite or meaningless values.
pub(crate) fn singular_values(a: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    if rows == 0 || cols == 0 {
        return Vec::new();
    }
    let (m, n, mut b) = if rows >= cols {
        (rows, cols, a.to_vec())
    } else {
        let mut t = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                t[c * rows + r] = a[r * cols + c];
            }
        }
        (cols, rows, t)
    };
    for _ in 0..MAX_SWEEPS {
        let mut rotated = false;
        for p in 0..n {
            for q in (p + 1)..n {
                let mut alpha = 0.0;
                let mut beta = 0.0;
                let mut gamma = 0.0;
                for i in 0..m {
                    let bp = b[i * n + p];
                    let bq = b[i * n + q];
                    alpha += bp * bp;
                    beta += bq * bq;
                    gamma += bp * bq;
                }
                if gamma == 0.0 || gamma.abs() <= ORTHOGONALITY_TOL * (alpha * beta).sqrt() {
                    continue;
                }
                rotated = true;
                // Rotation that makes columns p and q orthogonal:
                //   b_p' = c b_p - s b_q,   b_q' = s b_p + c b_q.
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = zeta.signum() / (zeta.abs() + zeta.hypot(1.0));
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for i in 0..m {
                    let bp = b[i * n + p];
                    let bq = b[i * n + q];
                    b[i * n + p] = c * bp - s * bq;
                    b[i * n + q] = s * bp + c * bq;
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let mut values: Vec<f64> = (0..n)
        .map(|j| {
            (0..m)
                .map(|i| b[i * n + j] * b[i * n + j])
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    values.sort_by(|x, y| y.total_cmp(x));
    values
}
