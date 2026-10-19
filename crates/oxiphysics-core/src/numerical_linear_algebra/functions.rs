//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute the dot product of two vectors.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Compute the L2 norm of a vector.
pub fn norm2(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}
/// Subtract: out = a - b
pub fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
/// Add: out = a + b
pub fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
/// Scale: out = alpha * a
pub fn vec_scale(alpha: f64, a: &[f64]) -> Vec<f64> {
    a.iter().map(|x| alpha * x).collect()
}
/// axpy: y += alpha * x
pub fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}
/// Dense matrix-vector product y = A * x, A stored row-major with shape (m, n).
pub fn matvec_dense(a: &[f64], m: usize, n: usize, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0f64; m];
    for i in 0..m {
        for j in 0..n {
            y[i] += a[i * n + j] * x[j];
        }
    }
    y
}
/// Transpose-vector product y = A^T * x.
pub fn matvec_dense_t(a: &[f64], m: usize, n: usize, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0f64; n];
    for i in 0..m {
        for j in 0..n {
            y[j] += a[i * n + j] * x[i];
        }
    }
    y
}
/// Trait for preconditioners: given a residual vector r, return M^{-1} r.
pub trait Preconditioner {
    /// Apply the preconditioner: solve M z = r and return z.
    fn apply(&self, r: &[f64]) -> Vec<f64>;
}
/// Block matrix-vector product: compute `C = A * B` where `A` is `m×n` and `B` is `n×k`.
pub fn matvec_block(a: &[f64], m: usize, n: usize, b: &[f64], k: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; m * k];
    for i in 0..m {
        for l in 0..n {
            if a[i * n + l] == 0.0 {
                continue;
            }
            for j in 0..k {
                c[i * k + j] += a[i * n + l] * b[l * k + j];
            }
        }
    }
    c
}
/// Block transposed matrix-vector product: compute `C = A^T * B` where `A` is `m×n` and `B` is `m×k`.
pub fn matvec_block_t(a: &[f64], m: usize, n: usize, b: &[f64], k: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; n * k];
    for i in 0..m {
        for l in 0..k {
            for j in 0..n {
                c[j * k + l] += a[i * n + j] * b[i * k + l];
            }
        }
    }
    c
}
/// Rectangular matrix multiply: compute `C = A * B` where `A` is `m×k` and `B` is `k×p`.
pub fn matmul_rect(a: &[f64], m: usize, k: usize, b: &[f64], _k2: usize, p: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; m * p];
    for i in 0..m {
        for l in 0..k {
            for j in 0..p {
                c[i * p + j] += a[i * k + l] * b[l * p + j];
            }
        }
    }
    c
}
/// Rectangular transposed matrix multiply: compute `C = A^T * B` where `A` is `m×k` and `B` is `m×n`.
pub fn matmul_rect_t(a: &[f64], m: usize, k: usize, b: &[f64], _m2: usize, n: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; k * n];
    for i in 0..m {
        for l in 0..k {
            for j in 0..n {
                c[l * n + j] += a[i * k + l] * b[i * n + j];
            }
        }
    }
    c
}
/// Relative remainder (after orthogonalisation) below which [`qr_thin`] treats a
/// column as linearly dependent on its predecessors.
const QR_DEPENDENCE_TOL: f64 = 64.0 * f64::EPSILON;
/// Maximum number of cyclic sweeps of the one-sided Jacobi iteration in
/// [`svd_small`].  Convergence is quadratic, so a handful of sweeps suffices.
const SVD_MAX_SWEEPS: usize = 60;
/// Unit vector orthogonal to every vector of `basis`, or `None` when `basis`
/// already spans `R^dim`.
///
/// `basis` must hold orthonormal vectors of length `dim`.  The canonical basis
/// vector with the largest component outside their span is orthogonalised twice
/// and normalised; its component is at least `sqrt((dim - basis.len()) / dim)`,
/// so the result is deterministic and well conditioned.
fn orthogonal_unit_vector(basis: &[Vec<f64>], dim: usize) -> Option<Vec<f64>> {
    if basis.len() >= dim {
        return None;
    }
    let mut best: Vec<f64> = Vec::new();
    let mut best_norm = 0.0f64;
    for candidate in 0..dim {
        let mut r = vec![0.0f64; dim];
        r[candidate] = 1.0;
        for _pass in 0..2 {
            for q in basis {
                let projection = dot(q, &r);
                axpy(-projection, q, &mut r);
            }
        }
        let nrm = norm2(&r);
        if nrm > best_norm {
            best_norm = nrm;
            best = r;
        }
    }
    if best_norm > 0.0 {
        for value in best.iter_mut() {
            *value /= best_norm;
        }
        Some(best)
    } else {
        None
    }
}
/// Compute the thin QR factorisation of an `m×k` matrix stored row-major.
///
/// Returns the orthonormal `Q` factor as a flat `m×k` row-major vector.  Columns
/// are orthogonalised against their predecessors twice (Gram-Schmidt with
/// re-orthogonalisation).  A column that is numerically dependent on the earlier
/// ones — as happens for a rank-deficient sketch — is replaced by a unit vector
/// orthogonal to all previous columns, so `Q` has orthonormal columns whenever
/// `m >= k`.
pub fn qr_thin(a: &[f64], m: usize, k: usize) -> Vec<f64> {
    let mut columns: Vec<Vec<f64>> = Vec::with_capacity(k);
    for j in 0..k {
        let mut v: Vec<f64> = (0..m).map(|i| a[i * k + j]).collect();
        let original = norm2(&v);
        for _pass in 0..2 {
            for q in &columns {
                let projection = dot(q, &v);
                axpy(-projection, q, &mut v);
            }
        }
        let remaining = norm2(&v);
        if remaining <= QR_DEPENDENCE_TOL * original {
            let replacement = orthogonal_unit_vector(&columns, m).unwrap_or_else(|| vec![0.0; m]);
            columns.push(replacement);
        } else {
            for value in v.iter_mut() {
                *value /= remaining;
            }
            columns.push(v);
        }
    }
    let mut q = vec![0.0f64; m * k];
    for (j, column) in columns.iter().enumerate() {
        for (i, &value) in column.iter().enumerate() {
            q[i * k + j] = value;
        }
    }
    q
}
/// Transpose a row-major `rows×cols` matrix into a row-major `cols×rows` matrix.
pub(super) fn transpose_dense(a: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut t = vec![0.0f64; rows * cols];
    for i in 0..rows {
        for j in 0..cols {
            t[j * rows + i] = a[i * cols + j];
        }
    }
    t
}
/// Rotate columns `p < q` of `columns` in place: `(c_p, c_q) ← (c·c_p − s·c_q, s·c_p + c·c_q)`.
fn rotate_column_pair(columns: &mut [Vec<f64>], p: usize, q: usize, c: f64, s: f64) {
    let (head, tail) = columns.split_at_mut(q);
    for (xp, xq) in head[p].iter_mut().zip(tail[0].iter_mut()) {
        let (vp, vq) = (*xp, *xq);
        *xp = c * vp - s * vq;
        *xq = s * vp + c * vq;
    }
}
/// One-sided Jacobi (Hestenes) SVD of a tall row-major `rows×cols` matrix
/// (`rows >= cols >= 1`).
///
/// Returns `(U, s, V)`: `U` is `rows×cols` with orthonormal columns, `s` holds the
/// singular values in descending order and `V` is the orthogonal `cols×cols`
/// matrix of right singular vectors (stored as `V`, not `Vᵀ`), so that
/// `A = U · diag(s) · Vᵀ`.  Pairs of columns are rotated until they are mutually
/// orthogonal; the column norms are then the singular values.  Columns belonging
/// to a zero singular value are completed to an orthonormal set.
fn jacobi_svd_tall(a: &[f64], rows: usize, cols: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // Scale by the largest entry so that squared column norms cannot overflow.
    let scale = a
        .iter()
        .take(rows * cols)
        .fold(0.0f64, |acc, v| acc.max(v.abs()));
    let divisor = if scale > 0.0 && scale.is_finite() {
        scale
    } else {
        1.0
    };
    let mut work: Vec<Vec<f64>> = (0..cols)
        .map(|j| (0..rows).map(|i| a[i * cols + j] / divisor).collect())
        .collect();
    let mut right: Vec<Vec<f64>> = (0..cols)
        .map(|j| {
            let mut unit = vec![0.0f64; cols];
            unit[j] = 1.0;
            unit
        })
        .collect();
    let tol = f64::EPSILON * (rows as f64).sqrt();
    for _sweep in 0..SVD_MAX_SWEEPS {
        let mut rotated = false;
        for p in 0..cols.saturating_sub(1) {
            for q in (p + 1)..cols {
                let alpha = dot(&work[p], &work[p]);
                let beta = dot(&work[q], &work[q]);
                let gamma = dot(&work[p], &work[q]);
                if gamma == 0.0 || gamma.abs() <= tol * alpha.sqrt() * beta.sqrt() {
                    continue;
                }
                rotated = true;
                // Rotation angle that makes the two columns orthogonal.
                let zeta = (beta - alpha) / (2.0 * gamma);
                let sign = if zeta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (zeta.abs() + zeta.hypot(1.0));
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                rotate_column_pair(&mut work, p, q, c, s);
                rotate_column_pair(&mut right, p, q, c, s);
            }
        }
        if !rotated {
            break;
        }
    }
    let norms: Vec<f64> = work.iter().map(|column| norm2(column)).collect();
    let mut order: Vec<usize> = (0..cols).collect();
    order.sort_by(|&x, &y| norms[y].total_cmp(&norms[x]));
    let mut u_columns: Vec<Vec<f64>> = Vec::with_capacity(cols);
    let mut v_columns: Vec<Vec<f64>> = Vec::with_capacity(cols);
    let mut singular_values: Vec<f64> = Vec::with_capacity(cols);
    let mut missing: Vec<usize> = Vec::new();
    for (position, &source) in order.iter().enumerate() {
        let nrm = norms[source];
        singular_values.push(nrm * divisor);
        v_columns.push(right[source].clone());
        if nrm > f64::MIN_POSITIVE {
            u_columns.push(work[source].iter().map(|x| x / nrm).collect());
        } else {
            u_columns.push(Vec::new());
            missing.push(position);
        }
    }
    for position in missing {
        let basis: Vec<Vec<f64>> = u_columns
            .iter()
            .filter(|column| !column.is_empty())
            .cloned()
            .collect();
        u_columns[position] =
            orthogonal_unit_vector(&basis, rows).unwrap_or_else(|| vec![0.0; rows]);
    }
    let mut u = vec![0.0f64; rows * cols];
    for (j, column) in u_columns.iter().enumerate() {
        for (i, &value) in column.iter().enumerate() {
            u[i * cols + j] = value;
        }
    }
    let mut v = vec![0.0f64; cols * cols];
    for (j, column) in v_columns.iter().enumerate() {
        for (i, &value) in column.iter().enumerate() {
            v[i * cols + j] = value;
        }
    }
    (u, singular_values, v)
}
/// Compute the thin SVD of a small dense `m×n` matrix.
///
/// `a` is stored row-major with `a.len() == m * n`.  With `k = min(m, n)` the
/// result `(U, s, Vt)` satisfies `A = U · diag(s) · Vt`, where
/// * `U` is `m×k` (row-major, stride `k`) with orthonormal columns,
/// * `s` holds the `k` singular values in descending order,
/// * `Vt` is `k×n` (row-major, stride `n`) with orthonormal rows.
///
/// The factorisation uses one-sided Jacobi rotations (Hestenes), which computes
/// small singular values to high relative accuracy.  A wide matrix (`m < n`) is
/// factorised through its transpose.  Columns of `U` and rows of `Vt` that belong
/// to zero singular values are completed so that both factors stay orthonormal.
/// Empty input (`m == 0` or `n == 0`) yields empty factors.
pub fn svd_small(a: &[f64], m: usize, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    if m == 0 || n == 0 {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    if m >= n {
        let (u, s, v) = jacobi_svd_tall(a, m, n);
        (u, s, transpose_dense(&v, n, n))
    } else {
        // A^T = U' S V'^T  =>  A = V' S U'^T.
        let a_t = transpose_dense(a, m, n);
        let (u_t, s, v_t) = jacobi_svd_tall(&a_t, n, m);
        (v_t, s, transpose_dense(&u_t, n, m))
    }
}
/// Compute the Gram matrix G = A^T A of shape (r x r) for A of shape (m x r).
pub fn gram(a: &[f64], m: usize, r: usize) -> Vec<f64> {
    let mut g = vec![0.0f64; r * r];
    for i in 0..m {
        for p in 0..r {
            for q in 0..r {
                g[p * r + q] += a[i * r + p] * a[i * r + q];
            }
        }
    }
    g
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::numerical_linear_algebra::EigenSolver;
    use crate::numerical_linear_algebra::ILU0Preconditioner;
    use crate::numerical_linear_algebra::ILUTPreconditioner;
    use crate::numerical_linear_algebra::IterativeSolver;
    use crate::numerical_linear_algebra::JacobiPreconditioner;
    use crate::numerical_linear_algebra::LowRankApproximation;
    use crate::numerical_linear_algebra::MatrixFunctions;
    use crate::numerical_linear_algebra::SPAIPreconditioner;
    use crate::numerical_linear_algebra::SSORPreconditioner;
    use crate::numerical_linear_algebra::SparseCSR;
    use crate::numerical_linear_algebra::TensorDecomposition;
    fn laplacian(n: usize) -> SparseCSR {
        SparseCSR::laplacian_1d(n)
    }
    fn rhs(n: usize) -> Vec<f64> {
        (0..n).map(|i| (i + 1) as f64).collect()
    }
    #[test]
    fn test_csr_matvec_identity() {
        let n = 4;
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for i in 0..n {
            rows.push(i);
            cols.push(i);
            vals.push(1.0);
        }
        let id = SparseCSR::from_triplets(n, n, &rows, &cols, &vals);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(id.matvec(&x), x);
    }
    #[test]
    fn test_csr_diagonal() {
        let lap = laplacian(5);
        let d = lap.diagonal();
        assert!(d.iter().all(|&v| (v - 2.0).abs() < 1e-12));
    }
    #[test]
    fn test_csr_laplacian_matvec() {
        let lap = laplacian(3);
        let x = vec![1.0, 1.0, 1.0];
        let y = lap.matvec(&x);
        assert!((y[0] - 1.0).abs() < 1e-12);
        assert!((y[1]).abs() < 1e-12);
        assert!((y[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_cg_small() {
        let n = 5;
        let lap = laplacian(n);
        let b = rhs(n);
        let res = IterativeSolver::conjugate_gradient(&lap, &b, None, 1e-8, 200);
        assert!(res.converged, "CG did not converge");
        let ax = lap.matvec(&res.x);
        for i in 0..n {
            assert!((ax[i] - b[i]).abs() < 1e-6, "CG residual too large at {i}");
        }
    }
    #[test]
    fn test_cg_residual_decreases() {
        let lap = laplacian(10);
        let b = rhs(10);
        let res = IterativeSolver::conjugate_gradient(&lap, &b, None, 1e-10, 100);
        assert!(res.residual < 1e-6);
    }
    #[test]
    fn test_cg_with_initial_guess() {
        let n = 6;
        let lap = laplacian(n);
        let b = rhs(n);
        let x0 = vec![0.1; n];
        let res = IterativeSolver::conjugate_gradient(&lap, &b, Some(&x0), 1e-8, 200);
        assert!(res.converged);
    }
    #[test]
    fn test_minres_runs() {
        let lap = laplacian(5);
        let b = rhs(5);
        let res = IterativeSolver::minres(&lap, &b, None, 1e-6, 100);
        assert!(res.iterations > 0);
    }
    #[test]
    fn test_gmres_small() {
        let n = 4;
        let lap = laplacian(n);
        let b = rhs(n);
        let res = IterativeSolver::gmres(&lap, &b, None, 1e-8, 50, 10);
        let ax = lap.matvec(&res.x);
        for i in 0..n {
            assert!((ax[i] - b[i]).abs() < 1e-4, "GMRES residual at {i}");
        }
    }
    #[test]
    fn test_gmres_restart() {
        let lap = laplacian(8);
        let b = rhs(8);
        let res = IterativeSolver::gmres(&lap, &b, None, 1e-6, 100, 4);
        assert!(res.iterations > 0);
    }
    #[test]
    fn test_bicgstab_converges() {
        let n = 6;
        let lap = laplacian(n);
        let b = rhs(n);
        let res = IterativeSolver::bicgstab(&lap, &b, None, 1e-8, 200);
        assert!(res.converged, "BiCGSTAB did not converge");
        let ax = lap.matvec(&res.x);
        for i in 0..n {
            assert!((ax[i] - b[i]).abs() < 1e-5, "BiCGSTAB residual at {i}");
        }
    }
    #[test]
    fn test_tfqmr_runs() {
        let lap = laplacian(5);
        let b = rhs(5);
        let res = IterativeSolver::tfqmr(&lap, &b, None, 1e-6, 100);
        assert!(res.iterations > 0);
        assert!(res.converged, "TFQMR residual {}", res.residual);
        let ax = lap.matvec(&res.x);
        for i in 0..5 {
            assert!((ax[i] - b[i]).abs() < 1e-5, "TFQMR residual at {i}");
        }
    }
    #[test]
    fn test_jacobi_preconditioner() {
        let lap = laplacian(5);
        let prec = JacobiPreconditioner::new(&lap);
        let r = vec![2.0; 5];
        let z = prec.apply(&r);
        for &v in &z {
            assert!((v - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_ssor_apply() {
        let lap = laplacian(4);
        let prec = SSORPreconditioner::new(lap, 1.0);
        let r = rhs(4);
        let z = prec.apply(&r);
        assert_eq!(z.len(), 4);
    }
    #[test]
    fn test_ilu0_apply() {
        let lap = laplacian(5);
        let prec = ILU0Preconditioner::new(&lap);
        let r = rhs(5);
        let z = prec.apply(&r);
        assert_eq!(z.len(), 5);
        assert!(z.iter().all(|v| v.is_finite()));
    }
    #[test]
    fn test_ilut_apply() {
        let lap = laplacian(5);
        let prec = ILUTPreconditioner::new(&lap, 0.01, 5);
        let r = rhs(5);
        let z = prec.apply(&r);
        assert_eq!(z.len(), 5);
    }
    #[test]
    fn test_spai_apply() {
        let lap = laplacian(5);
        let prec = SPAIPreconditioner::new(&lap);
        let r = rhs(5);
        let z = prec.apply(&r);
        assert_eq!(z.len(), 5);
    }
    #[test]
    fn test_power_iteration_dominant() {
        let lap = laplacian(5);
        let res = EigenSolver::power_iteration(&lap, 1e-8, 500);
        assert!(
            res.eigenvalues[0] > 3.0 && res.eigenvalues[0] < 4.0,
            "Power iteration eigenvalue = {}",
            res.eigenvalues[0]
        );
    }
    #[test]
    fn test_inverse_iteration_smallest() {
        let lap = laplacian(5);
        let res = EigenSolver::inverse_iteration(&lap, 0.0, 1e-6, 50);
        assert!(res.eigenvalues[0].is_finite());
    }
    #[test]
    fn test_arnoldi_returns_k_values() {
        let lap = laplacian(8);
        let res = EigenSolver::arnoldi(&lap, 3, 6, 1e-8);
        assert_eq!(res.eigenvalues.len(), 3);
        assert_eq!(res.eigenvectors.len(), 3);
        for w in res.eigenvalues.windows(2) {
            assert!(w[0] <= w[1], "eigenvalues not sorted");
        }
        // Ritz values are enclosed by the spectrum of the 8-point Laplacian.
        let upper = 2.0 - 2.0 * (8.0 * std::f64::consts::PI / 9.0).cos();
        assert!(
            res.eigenvalues
                .iter()
                .all(|&v| v > 0.0 && v < upper + 1e-12)
        );
    }
    #[test]
    fn test_hessenberg_eigenvalues_symmetric_2x2() {
        // [[2,-1],[-1,2]] has eigenvalues 1 and 3.  The old fabrication returned
        // the diagonal {2,2}.
        let h = vec![2.0, -1.0, -1.0, 2.0];
        let ev = EigenSolver::hessenberg_eigenvalues(&h, 2);
        assert_eq!(ev.len(), 2);
        assert!((ev[0] - 1.0).abs() < 1e-9, "ev0 = {}", ev[0]);
        assert!((ev[1] - 3.0).abs() < 1e-9, "ev1 = {}", ev[1]);
    }
    #[test]
    fn test_hessenberg_eigenvalues_tridiagonal_3x3() {
        // 1-D Laplacian (n=3): eigenvalues 2-√2, 2, 2+√2.  Old code: {2,2,2}.
        let h = vec![2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0];
        let ev = EigenSolver::hessenberg_eigenvalues(&h, 3);
        let s2 = std::f64::consts::SQRT_2;
        let expected = [2.0 - s2, 2.0, 2.0 + s2];
        assert_eq!(ev.len(), 3);
        for (got, want) in ev.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-8, "got {got}, want {want}");
        }
    }
    #[test]
    fn test_hessenberg_eigenvalues_nonsymmetric_companion() {
        // Upper-Hessenberg companion matrix of (x-1)(x-2)(x-4) =
        // x³ - 7x² + 14x - 8, eigenvalues {1, 2, 4}.  This is non-symmetric, so
        // the diagonal {7,0,0} bears no resemblance to the real spectrum.
        let h = vec![
            7.0, -14.0, 8.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0,
        ];
        let ev = EigenSolver::hessenberg_eigenvalues(&h, 3);
        let expected = [1.0, 2.0, 4.0];
        assert_eq!(ev.len(), 3);
        for (got, want) in ev.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-6, "got {got}, want {want}");
        }
    }
    #[test]
    fn test_lanczos_symmetric() {
        let lap = laplacian(8);
        let res = EigenSolver::lanczos(&lap, 4, 8, 1e-8);
        assert_eq!(res.eigenvalues.len(), 4);
        assert!(res.converged);
        for w in res.eigenvalues.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "eigenvalues not sorted");
        }
        // The four largest eigenvalues 2 - 2 cos(j pi / 9), j = 5..8.
        for (idx, j) in (5..=8).enumerate() {
            let want = 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / 9.0).cos();
            assert!(
                (res.eigenvalues[idx] - want).abs() < 1e-8,
                "j = {j}: {} vs {want}",
                res.eigenvalues[idx]
            );
        }
    }
    #[test]
    fn test_lobpcg_finds_eigenvalues() {
        let lap = laplacian(6);
        let res = EigenSolver::lobpcg(&lap, 2, 1e-6, 100);
        assert_eq!(res.eigenvalues.len(), 2);
        assert!(res.converged);
        // The two smallest eigenvalues 2 - 2 cos(j pi / 7), j = 1, 2.
        for j in 1..=2 {
            let want = 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / 7.0).cos();
            assert!((res.eigenvalues[j - 1] - want).abs() < 1e-6, "j = {j}");
        }
    }
    #[test]
    fn test_expm_identity() {
        let id = MatrixFunctions::eye(3);
        let exp_id = MatrixFunctions::expm(&id, 3);
        let e = std::f64::consts::E;
        for i in 0..3 {
            assert!((exp_id[i * 3 + i] - e).abs() < 1e-6, "exp(I) diagonal");
        }
    }
    #[test]
    fn test_expm_zero() {
        let z = vec![0.0f64; 9];
        let exp_z = MatrixFunctions::expm(&z, 3);
        let id = MatrixFunctions::eye(3);
        for (a, b) in exp_z.iter().zip(&id) {
            assert!((a - b).abs() < 1e-8);
        }
    }
    #[test]
    fn test_sqrtm_positive_definite() {
        let a = vec![4.0, 0.0, 0.0, 9.0];
        let sq = MatrixFunctions::sqrtm(&a, 2);
        assert!((sq[0] - 2.0).abs() < 1e-6, "sqrt(4) = {}", sq[0]);
        assert!((sq[3] - 3.0).abs() < 1e-6, "sqrt(9) = {}", sq[3]);
    }
    #[test]
    fn test_matrix_poly() {
        let id = MatrixFunctions::eye(3);
        let a: Vec<f64> = id.iter().map(|v| v * 2.0).collect();
        let res = MatrixFunctions::matrix_poly(&a, 3, &[1.0, 1.0]);
        for i in 0..3 {
            assert!((res[i * 3 + i] - 3.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_inv_dense_2x2() {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let inv = MatrixFunctions::inv_dense(&a, 2);
        assert!((inv[0] - 0.6).abs() < 1e-10);
        assert!((inv[1] + 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_randomized_svd_rank_1() {
        let m = 4;
        let n = 4;
        let a: Vec<f64> = (0..m * n)
            .map(|idx| {
                let i = idx / n;
                let j = idx % n;
                (i + 1) as f64 * (j + 1) as f64
            })
            .collect();
        let res = LowRankApproximation::randomized_svd(&a, m, n, 1, 2);
        assert_eq!(res.s.len(), 1);
        assert!(res.s[0] > 0.0);
    }
    #[test]
    fn test_randomized_svd_shape() {
        let res = LowRankApproximation::randomized_svd(&vec![1.0f64; 6 * 8], 6, 8, 3, 1);
        assert_eq!(res.k, 3);
        assert_eq!(res.m, 6);
        assert_eq!(res.n, 8);
    }
    #[test]
    fn test_nystrom_output_shape() {
        let n = 6;
        let k = 3;
        let a: Vec<f64> = (0..n * n)
            .map(|idx| {
                let i = idx / n;
                let j = idx % n;
                if i == j { (i + 1) as f64 } else { 0.0 }
            })
            .collect();
        let res = LowRankApproximation::nystrom(&a, n, k);
        assert_eq!(res.k, k);
    }
    #[test]
    fn test_cur_decomposition_shapes() {
        let m = 4;
        let n = 5;
        let k = 2;
        let a: Vec<f64> = (0..m * n).map(|i| i as f64).collect();
        let (c, u, r) = LowRankApproximation::cur_decomposition(&a, m, n, k);
        assert_eq!(c.len(), m * k);
        assert_eq!(u.len(), k * k);
        assert_eq!(r.len(), k * n);
    }
    #[test]
    fn test_cp_als_runs() {
        let shape = [3, 3, 3];
        let x: Vec<f64> = (0..27).map(|i| i as f64).collect();
        let res = TensorDecomposition::cp_als(&x, &shape, 2, 20, 1e-6);
        assert_eq!(res.factors.len(), 3);
        assert_eq!(res.core.len(), 2);
    }
    #[test]
    fn test_cp_als_factor_shapes() {
        let shape = [4, 3, 2];
        let x: Vec<f64> = (0..24).map(|i| i as f64 * 0.1).collect();
        let res = TensorDecomposition::cp_als(&x, &shape, 3, 10, 1e-4);
        assert_eq!(res.factors[0].len(), 4 * 3);
        assert_eq!(res.factors[1].len(), 3 * 3);
        assert_eq!(res.factors[2].len(), 2 * 3);
    }
    #[test]
    fn test_tucker_hosvd_runs() {
        let shape = [4, 3, 2];
        let x: Vec<f64> = (0..24).map(|i| i as f64).collect();
        let res = TensorDecomposition::tucker_hosvd(&x, &shape, &[2, 2, 2]);
        assert_eq!(res.factors.len(), 3);
        assert!(!res.core.is_empty());
    }
    #[test]
    fn test_tucker_core_shape() {
        let shape = [3, 3, 3];
        let x: Vec<f64> = (0..27).map(|i| i as f64).collect();
        let ranks = [2, 2, 2];
        let res = TensorDecomposition::tucker_hosvd(&x, &shape, &ranks);
        assert_eq!(res.core.len(), 2 * 2 * 2);
    }
    #[test]
    fn test_tensor_train_runs() {
        let shape = [4, 3, 2];
        let x: Vec<f64> = (0..24).map(|i| (i as f64).sin()).collect();
        let res = TensorDecomposition::tensor_train(&x, &shape, 4, 1e-6);
        assert_eq!(res.shape, shape);
        assert_eq!(res.factors.len(), 3);
    }
    #[test]
    fn test_tensor_train_shape_preserved() {
        let shape = [2, 3, 4];
        let x: Vec<f64> = (0..24).map(|i| i as f64).collect();
        let res = TensorDecomposition::tensor_train(&x, &shape, 2, 1e-4);
        assert_eq!(res.shape, vec![2, 3, 4]);
    }
    #[test]
    fn test_preconditioned_cg_jacobi() {
        let n = 8;
        let lap = laplacian(n);
        let prec = JacobiPreconditioner::new(&lap);
        let b = rhs(n);
        let pb: Vec<f64> = prec.apply(&b);
        let res = IterativeSolver::conjugate_gradient(&lap, &pb, None, 1e-6, 200);
        assert!(res.residual.is_finite());
    }
    #[test]
    fn test_gram_matrix_symmetry() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let g = gram(&a, 3, 2);
        assert!((g[1] - g[2]).abs() < 1e-12);
    }
    #[test]
    fn test_matvec_block_transpose() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.0, 0.0, 0.0, 1.0];
        let c = matvec_block(&a, 2, 2, &b, 2);
        for i in 0..4 {
            assert!((c[i] - a[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_svd_small_singular_values_nonneg() {
        let a = vec![3.0, 0.0, 4.0, 0.0, 2.0, 0.0];
        let (_, s, _) = svd_small(&a, 2, 3);
        for &sv in &s {
            assert!(sv >= 0.0);
        }
    }
    #[test]
    fn test_qr_thin_orthonormality() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let q = qr_thin(&a, 3, 2);
        let g = gram(&q, 3, 2);
        assert!((g[0] - 1.0).abs() < 1e-8, "g[0] = {}", g[0]);
        assert!((g[3] - 1.0).abs() < 1e-8, "g[3] = {}", g[3]);
        assert!(g[1].abs() < 1e-8, "g[1] = {}", g[1]);
    }
    #[test]
    fn test_logm_exp_inverse() {
        let a = vec![4.0, 0.0, 0.0, 9.0];
        let log_a = MatrixFunctions::logm(&a, 2);
        assert!(log_a[0].is_finite(), "log_a[0] = {}", log_a[0]);
        assert!(log_a[3].is_finite(), "log_a[3] = {}", log_a[3]);
        assert!(
            (log_a[0] - 4.0_f64.ln()).abs() < 1e-9,
            "log(4) = {}, got {}",
            4.0_f64.ln(),
            log_a[0]
        );
        assert!(
            (log_a[3] - 9.0_f64.ln()).abs() < 1e-9,
            "log(9) = {}, got {}",
            9.0_f64.ln(),
            log_a[3]
        );
    }
    #[test]
    fn test_dot_and_norm() {
        let a = vec![3.0, 4.0];
        assert!((norm2(&a) - 5.0).abs() < 1e-12);
        assert!((dot(&a, &a) - 25.0).abs() < 1e-12);
    }
    #[test]
    fn test_axpy() {
        let x = vec![1.0, 2.0, 3.0];
        let mut y = vec![10.0, 20.0, 30.0];
        axpy(2.0, &x, &mut y);
        assert_eq!(y, vec![12.0, 24.0, 36.0]);
    }
}
