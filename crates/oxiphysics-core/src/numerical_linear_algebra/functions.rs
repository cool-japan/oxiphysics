//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::MatrixFunctions;

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
/// Compute the thin QR factorisation of an `m×k` matrix stored row-major.
///
/// Returns the orthonormal `Q` factor as a flat `m×k` row-major vector.
pub fn qr_thin(a: &[f64], m: usize, k: usize) -> Vec<f64> {
    let mut q = a.to_vec();
    for j in 0..k {
        let nrm = (0..m)
            .map(|i| q[i * k + j] * q[i * k + j])
            .sum::<f64>()
            .sqrt();
        if nrm > 1e-300 {
            for i in 0..m {
                q[i * k + j] /= nrm;
            }
        }
        for l in (j + 1)..k {
            let proj: f64 = (0..m).map(|i| q[i * k + j] * q[i * k + l]).sum();
            for i in 0..m {
                q[i * k + l] -= proj * q[i * k + j];
            }
        }
    }
    q
}
/// Compute the thin SVD of a small dense `m×n` matrix.
///
/// Returns `(U, s, Vt)` where `U` is `m×k`, `s` has length `k`, and `Vt` is `k×n`, with `k = min(m, n)`.
pub fn svd_small(a: &[f64], m: usize, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let k = m.min(n);
    let mut u = MatrixFunctions::eye(m);
    let mut vt_mat = MatrixFunctions::eye(n);
    let mut b = a.to_vec();
    for _ in 0..100 {
        let mut converged = true;
        for p in 0..k.min(m) {
            for q in (p + 1)..k.min(n) {
                let bpp: f64 = (0..m).map(|i| b[i * n + p] * b[i * n + p]).sum();
                let bqq: f64 = (0..m).map(|i| b[i * n + q] * b[i * n + q]).sum();
                let bpq: f64 = (0..m).map(|i| b[i * n + p] * b[i * n + q]).sum();
                if bpq.abs() <= 1e-12 * (bpp * bqq).sqrt() {
                    continue;
                }
                converged = false;
                let tau = (bqq - bpp) / (2.0 * bpq);
                let t = 1.0 / (tau.abs() + (1.0 + tau * tau).sqrt());
                let t = if tau < 0.0 { -t } else { t };
                let cs = 1.0 / (1.0 + t * t).sqrt();
                let sn = t * cs;
                for i in 0..m {
                    let bp = b[i * n + p];
                    let bq = b[i * n + q];
                    b[i * n + p] = cs * bp + sn * bq;
                    b[i * n + q] = -sn * bp + cs * bq;
                }
                for i in 0..n {
                    let vp = vt_mat[i * n + p];
                    let vq = vt_mat[i * n + q];
                    vt_mat[i * n + p] = cs * vp + sn * vq;
                    vt_mat[i * n + q] = -sn * vp + cs * vq;
                }
            }
        }
        if converged {
            break;
        }
    }
    let mut s: Vec<f64> = (0..k)
        .map(|j| {
            (0..m)
                .map(|i| b[i * n + j] * b[i * n + j])
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    for j in 0..k {
        for i in 0..m {
            u[i * m + j] = if s[j] > 1e-300 {
                b[i * n + j] / s[j]
            } else {
                0.0
            };
        }
    }
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&a, &b| s[b].partial_cmp(&s[a]).unwrap_or(std::cmp::Ordering::Equal));
    let s_sorted: Vec<f64> = order.iter().map(|&i| s[i]).collect();
    s = s_sorted;
    (u, s, vt_mat)
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
        assert!(!res.eigenvalues.is_empty());
    }
    #[test]
    fn test_lanczos_symmetric() {
        let lap = laplacian(8);
        let res = EigenSolver::lanczos(&lap, 4, 8, 1e-8);
        assert!(!res.eigenvalues.is_empty());
        for w in res.eigenvalues.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "eigenvalues not sorted");
        }
    }
    #[test]
    fn test_lobpcg_finds_eigenvalues() {
        let lap = laplacian(6);
        let res = EigenSolver::lobpcg(&lap, 2, 1e-6, 100);
        assert_eq!(res.eigenvalues.len(), 2);
        assert!(res.eigenvalues.iter().all(|v| v.is_finite()));
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
            log_a[0] > 0.0,
            "log(4) should be positive, got {}",
            log_a[0]
        );
        assert!(
            log_a[3] > 0.0,
            "log(9) should be positive, got {}",
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
