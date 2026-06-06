//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
#[cfg(test)]
use super::types::*;
use crate::sparse::CsrMatrix;
/// GMRES solver for sparse systems (Arnoldi-based, restart = `max_iter`).
///
/// Solves `A x = b` iteratively using the Generalized Minimal Residual method.
///
/// # Arguments
/// * `a`        – sparse CSR matrix
/// * `b`        – right-hand side
/// * `max_iter` – maximum Krylov subspace size (= restart dimension)
/// * `tol`      – convergence tolerance on relative residual ‖r‖/‖b‖
pub fn fem_gmres(a: &CsrMatrix, b: &[f64], max_iter: usize, tol: f64) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    let b_norm: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    if b_norm < 1e-60 {
        return x;
    }
    let mut r: Vec<f64> = b.to_vec();
    for _restart in 0..=max_iter {
        let r_norm: f64 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if r_norm / b_norm < tol {
            break;
        }
        let m = max_iter.min(n);
        let mut v_basis: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
        let mut h = vec![vec![0.0_f64; m]; m + 1];
        let mut g = vec![0.0_f64; m + 1];
        let v0: Vec<f64> = r.iter().map(|&v| v / r_norm).collect();
        v_basis.push(v0);
        g[0] = r_norm;
        let mut cos_vals = vec![0.0_f64; m];
        let mut sin_vals = vec![0.0_f64; m];
        let mut k = 0usize;
        while k < m {
            let mut w = vec![0.0_f64; n];
            for (row, w_row) in w.iter_mut().enumerate() {
                for idx in a.row_ptr[row]..a.row_ptr[row + 1] {
                    *w_row += a.values[idx] * v_basis[k][a.col_indices[idx]];
                }
            }
            for j in 0..=k {
                h[j][k] = w
                    .iter()
                    .zip(v_basis[j].iter())
                    .map(|(wi, vj)| wi * vj)
                    .sum();
                for (i, w_i) in w.iter_mut().enumerate() {
                    *w_i -= h[j][k] * v_basis[j][i];
                }
            }
            let w_norm: f64 = w.iter().map(|v| v * v).sum::<f64>().sqrt();
            h[k + 1][k] = w_norm;
            for j in 0..k {
                let tmp = cos_vals[j] * h[j][k] + sin_vals[j] * h[j + 1][k];
                h[j + 1][k] = -sin_vals[j] * h[j][k] + cos_vals[j] * h[j + 1][k];
                h[j][k] = tmp;
            }
            let denom = (h[k][k] * h[k][k] + h[k + 1][k] * h[k + 1][k]).sqrt();
            if denom < 1e-60 {
                k += 1;
                break;
            }
            cos_vals[k] = h[k][k] / denom;
            sin_vals[k] = h[k + 1][k] / denom;
            h[k][k] = cos_vals[k] * h[k][k] + sin_vals[k] * h[k + 1][k];
            h[k + 1][k] = 0.0;
            g[k + 1] = -sin_vals[k] * g[k];
            g[k] *= cos_vals[k];
            if w_norm > 1e-60 {
                v_basis.push(w.iter().map(|&v| v / w_norm).collect());
            } else {
                v_basis.push(vec![0.0; n]);
            }
            k += 1;
            if g[k].abs() / b_norm < tol {
                break;
            }
        }
        let mut y = vec![0.0_f64; k];
        for i in (0..k).rev() {
            let mut s = g[i];
            for j in (i + 1)..k {
                s -= h[i][j] * y[j];
            }
            y[i] = if h[i][i].abs() > 1e-60 {
                s / h[i][i]
            } else {
                0.0
            };
        }
        for j in 0..k {
            for i in 0..n {
                x[i] += y[j] * v_basis[j][i];
            }
        }
        r = b.to_vec();
        for (row, r_row) in r.iter_mut().enumerate() {
            let ax_row: f64 = (a.row_ptr[row]..a.row_ptr[row + 1])
                .map(|idx| a.values[idx] * x[a.col_indices[idx]])
                .sum();
            *r_row -= ax_row;
        }
        if r.iter().map(|v| v * v).sum::<f64>().sqrt() / b_norm < tol {
            break;
        }
    }
    x
}
/// AMG smoother step: apply `n_smooth` Gauss-Seidel sweeps as a simple
/// Algebraic Multigrid smoother.
///
/// In a full AMG cycle this would alternate with coarse-grid corrections.
/// Here we expose just the smoothing step for testing and composition.
///
/// # Arguments
/// * `a`        – system matrix (CSR)
/// * `b`        – right-hand side
/// * `x_in`     – initial guess
/// * `n_smooth` – number of Gauss-Seidel sweeps
pub fn amg_smooth_step(a: &CsrMatrix, b: &[f64], x_in: &[f64], n_smooth: usize) -> Vec<f64> {
    let n = b.len();
    let mut x = x_in.to_vec();
    for _sweep in 0..n_smooth {
        for i in 0..n {
            let row_start = a.row_ptr[i];
            let row_end = a.row_ptr[i + 1];
            let mut diag = 0.0_f64;
            let mut off = 0.0_f64;
            for idx in row_start..row_end {
                let j = a.col_indices[idx];
                if j == i {
                    diag = a.values[idx];
                } else {
                    off += a.values[idx] * x[j];
                }
            }
            if diag.abs() > 1e-60 {
                x[i] = (b[i] - off) / diag;
            }
        }
    }
    x
}
/// ILU(0) preconditioned iterative solver.
///
/// Computes an ILU(0) factorisation of the dense matrix `a` (row-major, n×n),
/// then uses it as a preconditioner in a preconditioned Richardson iteration
/// until convergence or `max_iter` steps.
///
/// # Arguments
/// * `a`        – dense n×n matrix (row-major)
/// * `b`        – right-hand side
/// * `n`        – system size
/// * `max_iter` – maximum iterations
/// * `tol`      – convergence tolerance
pub fn ilu0_precond_solve(a: &[f64], b: &[f64], n: usize, max_iter: usize, tol: f64) -> Vec<f64> {
    let (l_mat, u_mat) = ilu0_dense(a, n);
    let mut x = vec![0.0_f64; n];
    for _iter in 0..max_iter {
        let mut r = b.to_vec();
        for i in 0..n {
            let ax_i: f64 = (0..n).map(|j| a[i * n + j] * x[j]).sum();
            r[i] -= ax_i;
        }
        let r_norm: f64 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if r_norm < tol {
            break;
        }
        let z = ilu0_solve(&l_mat, &u_mat, &r, n);
        for i in 0..n {
            x[i] += z[i];
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cg_small_system() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let x0 = [0.0, 0.0];
        let x = conjugate_gradient(&a, &b, &x0, 100, 1e-12);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-10, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-10, "x[1] = {}", x[1]);
    }
    #[test]
    fn test_cg_3x3() {
        let triplets = vec![
            (0, 0, 6.0),
            (0, 1, 2.0),
            (0, 2, 1.0),
            (1, 0, 2.0),
            (1, 1, 5.0),
            (1, 2, 1.0),
            (2, 0, 1.0),
            (2, 1, 1.0),
            (2, 2, 4.0),
        ];
        let a = CsrMatrix::from_triplets(3, 3, &triplets);
        let b = [9.0, 8.0, 6.0];
        let x0 = [0.0, 0.0, 0.0];
        let x = conjugate_gradient(&a, &b, &x0, 100, 1e-12);
        let ax = a.mul_vec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-8,
                "Ax[{i}] = {}, b[{i}] = {}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_preconditioned_cg() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let precond = jacobi_preconditioner(&a);
        let x = preconditioned_cg(&a, &b, &precond, 100, 1e-12);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-10, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-10, "x[1] = {}", x[1]);
    }
    #[test]
    fn test_sparse_matrix_cg_simple() {
        let triplets = vec![(0, 0, 4.0), (1, 1, 4.0)];
        let k = CsrMatrix::from_triplets(2, 2, &triplets);
        let f = [4.0, 8.0];
        let x0 = [0.0, 0.0];
        let x = conjugate_gradient(&k, &f, &x0, 100, 1e-12);
        assert!((x[0] - 1.0).abs() < 1e-10, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1] = {}", x[1]);
    }
    #[test]
    fn cg_struct_dense_2x2() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let b = vec![1.0, 2.0];
        let cg = ConjugateGradient::new(100, 1e-12);
        let x = cg.solve(&a, &b, 2);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-9, "x[0]={}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-9, "x[1]={}", x[1]);
    }
    #[test]
    fn cholesky_2x2_factorize_solve() {
        let chol = CholeskyDense::new(2);
        let a = vec![4.0, 2.0, 2.0, 3.0];
        let l = chol.factorize(&a).expect("should factorize");
        assert!((l[0] - 2.0).abs() < 1e-12, "L[0,0]={}", l[0]);
        assert!((l[2] - 1.0).abs() < 1e-12, "L[1,0]={}", l[2]);
        let b = vec![8.0, 7.0];
        let x = chol.solve(&l, &b);
        let ax0 = 4.0 * x[0] + 2.0 * x[1];
        let ax1 = 2.0 * x[0] + 3.0 * x[1];
        assert!((ax0 - 8.0).abs() < 1e-10, "A*x[0]={ax0}");
        assert!((ax1 - 7.0).abs() < 1e-10, "A*x[1]={ax1}");
    }
    #[test]
    fn cholesky_not_positive_definite_returns_error() {
        let chol = CholeskyDense::new(2);
        let a = vec![-1.0, 0.0, 0.0, 1.0];
        let result = chol.factorize(&a);
        assert_eq!(result, Err(SolverError::NotPositiveDefinite));
    }
    #[test]
    fn gauss_seidel_converges_diagonally_dominant() {
        let a = vec![10.0, 1.0, 1.0, 10.0];
        let b = vec![11.0, 11.0];
        let gs = GaussSeidel::new(1000, 1e-12, 1.0);
        let x = gs.solve(&a, &b, 2);
        assert!((x[0] - 1.0).abs() < 1e-8, "x[0]={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-8, "x[1]={}", x[1]);
    }
    #[test]
    fn gauss_seidel_stats_converged() {
        let a = vec![10.0, 1.0, 1.0, 10.0];
        let b = vec![11.0, 11.0];
        let gs = GaussSeidel::new(1000, 1e-10, 1.0);
        let (x, stats) = gs.solve_with_stats(&a, &b, 2);
        assert!(stats.converged, "should converge");
        assert!(stats.residual < 1e-8, "residual={}", stats.residual);
        assert!((x[0] - 1.0).abs() < 1e-8);
    }
    #[test]
    fn sparse_matvec_identity() {
        let values = vec![1.0, 1.0, 1.0];
        let col_idx = vec![0usize, 1, 2];
        let row_ptr = vec![0usize, 1, 2, 3];
        let x = vec![3.0, 4.0, 5.0];
        let y = sparse_matvec(&values, &col_idx, &row_ptr, &x);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 4.0).abs() < 1e-12);
        assert!((y[2] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn sparse_matvec_general() {
        let values = vec![2.0, 1.0, 1.0, 3.0];
        let col_idx = vec![0usize, 1, 0, 1];
        let row_ptr = vec![0usize, 2, 4];
        let x = vec![1.0, 2.0];
        let y = sparse_matvec(&values, &col_idx, &row_ptr, &x);
        assert!((y[0] - 4.0).abs() < 1e-12, "y[0]={}", y[0]);
        assert!((y[1] - 7.0).abs() < 1e-12, "y[1]={}", y[1]);
    }
    #[test]
    fn pcg_solve_converges() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let precond = jacobi_preconditioner(&a);
        let (x, stats) = pcg_solve(&a, &b, &precond, 100, 1e-12);
        assert!(stats.converged, "PCG did not converge");
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-9);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-9);
    }
    #[test]
    fn bicgstab_solves_2x2() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let b = vec![1.0, 2.0];
        let (x, stats) = bicgstab_dense(&a, &b, 2, 100, 1e-12);
        assert!(stats.converged, "BiCGSTAB did not converge");
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0]={}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1]={}", x[1]);
    }
    #[test]
    fn bicgstab_asymmetric_3x3() {
        let a = vec![3.0, 1.0, 0.0, -1.0, 4.0, 1.0, 0.0, -1.0, 3.0];
        let b = vec![4.0, 5.0, 3.0];
        let (x, stats) = bicgstab_dense(&a, &b, 3, 200, 1e-10);
        assert!(stats.converged, "BiCGSTAB 3x3 did not converge");
        for i in 0..3 {
            let ax_i: f64 = (0..3).map(|j| a[i * 3 + j] * x[j]).sum();
            assert!((ax_i - b[i]).abs() < 1e-7, "row {i}: Ax={ax_i}, b={}", b[i]);
        }
    }
    #[test]
    fn gmres_solves_2x2_spd() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let b = vec![1.0, 2.0];
        let (x, stats) = gmres_dense(&a, &b, 2, 10, 100, 1e-12);
        assert!(stats.converged, "GMRES did not converge");
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0]={}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1]={}", x[1]);
    }
    #[test]
    fn gmres_restart_3x3() {
        let a = vec![6.0, 2.0, 1.0, 2.0, 5.0, 1.0, 1.0, 1.0, 4.0];
        let b = vec![9.0, 8.0, 6.0];
        let (x, stats) = gmres_dense(&a, &b, 3, 3, 50, 1e-10);
        assert!(
            stats.converged || stats.residual < 1e-6,
            "GMRES restart: converged={}, res={}",
            stats.converged,
            stats.residual
        );
        for i in 0..3 {
            let ax_i: f64 = (0..3).map(|j| a[i * 3 + j] * x[j]).sum();
            assert!(
                (ax_i - b[i]).abs() < 1e-5,
                "GMRES row {i}: Ax={ax_i} b={}",
                b[i]
            );
        }
    }
    #[test]
    fn ilu0_dense_produces_l_u() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let (l, u) = ilu0_dense(&a, 2);
        assert!((l[0] - 1.0).abs() < 1e-12);
        assert!((l[1]).abs() < 1e-12);
        assert!((u[1] - 1.0).abs() < 1e-12);
        for i in 0..2 {
            for j in 0..2 {
                let lu_ij: f64 = (0..2).map(|k| l[i * 2 + k] * u[k * 2 + j]).sum();
                assert!(
                    (lu_ij - a[i * 2 + j]).abs() < 1e-10,
                    "LU[{i}][{j}] = {lu_ij}, A[{i}][{j}] = {}",
                    a[i * 2 + j]
                );
            }
        }
    }
    #[test]
    fn incomplete_cholesky_2x2() {
        let a = vec![4.0, 2.0, 2.0, 3.0];
        let l = incomplete_cholesky_dense(&a, 2);
        for i in 0..2 {
            for j in 0..2 {
                let llt_ij: f64 = (0..2).map(|k| l[i * 2 + k] * l[j * 2 + k]).sum();
                assert!(
                    (llt_ij - a[i * 2 + j]).abs() < 1e-8,
                    "LL^T[{i}][{j}]={llt_ij} vs A={}",
                    a[i * 2 + j]
                );
            }
        }
    }
    #[test]
    fn minres_solves_spd() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let b = vec![1.0, 2.0];
        let (_x, stats) = minres_dense(&a, &b, 2, 100, 1e-12);
        assert!(
            stats.residual < 1e-6 || stats.converged,
            "MINRES residual too large: res={}, converged={}",
            stats.residual,
            stats.converged
        );
    }
    #[test]
    fn schur_complement_2x2_plus_2x2() {
        let a = vec![4.0, 0.0, 0.0, 4.0];
        let b = vec![1.0, 0.0, 0.0, 1.0];
        let c = vec![1.0, 0.0, 0.0, 1.0];
        let d = vec![2.0, 0.0, 0.0, 2.0];
        let s = schur_complement_dense(&a, &b, &c, &d, 2);
        assert!((s[0] - 1.75).abs() < 1e-10, "S[0,0]={}", s[0]);
        assert!((s[3] - 1.75).abs() < 1e-10, "S[1,1]={}", s[3]);
    }
    #[test]
    fn multigrid_vcycle_reduces_residual() {
        let n = 4;
        let mut a = vec![0.0_f64; n * n];
        for i in 0..n {
            a[i * n + i] = 2.0;
            if i > 0 {
                a[i * n + (i - 1)] = -1.0;
            }
            if i < n - 1 {
                a[i * n + (i + 1)] = -1.0;
            }
        }
        let b = vec![1.0; n];
        let x0 = vec![0.0; n];
        let (x, _stats) = multigrid_vcycle_dense(&a, &b, &x0, n, 3, 3);
        let ax = dense_matvec(&a, &x, n);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let r_norm = dot_slice(&r, &r).sqrt();
        assert!(
            r_norm < 2.0,
            "multigrid residual should be reduced: {r_norm}"
        );
    }
    #[test]
    fn convergence_monitor_records_history() {
        let mut mon = ConvergenceMonitor::new(1e-10, 100);
        assert_eq!(mon.iterations(), 0);
        mon.record(1.0);
        mon.record(0.5);
        mon.record(1e-12);
        assert_eq!(mon.iterations(), 3);
        assert!(mon.converged(), "should converge after 1e-12 entry");
    }
    #[test]
    fn convergence_monitor_not_converged() {
        let mut mon = ConvergenceMonitor::new(1e-10, 100);
        mon.record(1.0);
        mon.record(0.1);
        assert!(!mon.converged(), "should not yet converge");
    }
    #[test]
    fn convergence_monitor_rate() {
        let mut mon = ConvergenceMonitor::new(1e-10, 100);
        mon.record(1.0);
        mon.record(0.5);
        let rate = mon.convergence_rate().unwrap();
        assert!((rate - 0.5).abs() < 1e-12, "rate = {rate}");
    }
    #[test]
    fn convergence_monitor_stagnation() {
        let mut mon = ConvergenceMonitor::new(1e-10, 100);
        mon.record(1.0);
        mon.record(1.0);
        assert!(mon.stagnated(), "should detect stagnation");
    }
    #[test]
    fn pcg_icc_solves_spd() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let (x, stats) = pcg_icc(&a, &b, 100, 1e-12);
        assert!(
            stats.converged || stats.residual < 1e-9,
            "PCG-ICC did not converge"
        );
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1] = {}", x[1]);
    }
    #[test]
    fn pcg_icc_tridiagonal() {
        let n = 5;
        let mut triplets = vec![];
        for i in 0..n {
            triplets.push((i, i, 4.0));
            if i > 0 {
                triplets.push((i, i - 1, -1.0));
            }
            if i < n - 1 {
                triplets.push((i, i + 1, -1.0));
            }
        }
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let b = vec![1.0; n];
        let (x, stats) = pcg_icc(&a, &b, 200, 1e-10);
        assert!(
            stats.converged || stats.residual < 1e-8,
            "PCG-ICC tridiagonal: res={}, converged={}",
            stats.residual,
            stats.converged
        );
        let ax = a.mul_vec(&x);
        for i in 0..n {
            assert!(
                (ax[i] - b[i]).abs() < 1e-7,
                "Ax[{i}] = {}, b = {}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn pcg_with_monitor_records_decreasing_residuals() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let precond = jacobi_preconditioner(&a);
        let (_x, monitor) = pcg_with_monitor(&a, &b, &precond, 50, 1e-12);
        assert!(monitor.converged(), "should converge");
        assert!(monitor.iterations() > 0);
    }
    #[test]
    fn bicgstab_sparse_solves_spd() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let (x, stats) = bicgstab_sparse(&a, &b, 100, 1e-12);
        assert!(
            stats.converged,
            "BiCGSTAB sparse did not converge: res={}",
            stats.residual
        );
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1] = {}", x[1]);
    }
    #[test]
    fn bicgstab_sparse_identity_system() {
        let triplets = vec![(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)];
        let a = CsrMatrix::from_triplets(3, 3, &triplets);
        let b = [3.0, 5.0, 7.0];
        let (x, stats) = bicgstab_sparse(&a, &b, 100, 1e-12);
        assert!(stats.converged, "BiCGSTAB identity: res={}", stats.residual);
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-10, "x[{i}] = {}", x[i]);
        }
    }
    #[test]
    fn minres_sparse_solves_spd() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let (_x, stats) = minres_sparse(&a, &b, 100, 1e-10);
        assert!(
            stats.residual < 1e-6 || stats.converged,
            "MINRES sparse: res={}, converged={}",
            stats.residual,
            stats.converged
        );
    }
    #[test]
    fn minres_sparse_diagonal_system() {
        let triplets = vec![(0, 0, 2.0), (1, 1, 3.0), (2, 2, 5.0)];
        let a = CsrMatrix::from_triplets(3, 3, &triplets);
        let b = [4.0, 9.0, 15.0];
        let (_x, stats) = minres_sparse(&a, &b, 200, 1e-8);
        assert!(
            stats.residual < 1.0 || stats.converged,
            "MINRES diagonal: res={}, converged={}",
            stats.residual,
            stats.converged
        );
    }
    #[test]
    fn gmres_sparse_solves_spd() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let (x, stats) = gmres_sparse(&a, &b, 10, 50, 1e-12);
        assert!(
            stats.converged,
            "GMRES sparse did not converge: res={}",
            stats.residual
        );
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1] = {}", x[1]);
    }
    #[test]
    fn gmres_sparse_3x3() {
        let triplets = vec![
            (0, 0, 6.0),
            (0, 1, 2.0),
            (0, 2, 1.0),
            (1, 0, 2.0),
            (1, 1, 5.0),
            (1, 2, 1.0),
            (2, 0, 1.0),
            (2, 1, 1.0),
            (2, 2, 4.0),
        ];
        let a = CsrMatrix::from_triplets(3, 3, &triplets);
        let b = [9.0, 8.0, 6.0];
        let (x, stats) = gmres_sparse(&a, &b, 3, 20, 1e-10);
        assert!(
            stats.converged || stats.residual < 1e-6,
            "GMRES sparse 3×3: converged={}, res={}",
            stats.converged,
            stats.residual
        );
        let ax = a.mul_vec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-5,
                "row {i}: ax={}, b={}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn amg_pcg_solves_1d_poisson() {
        let n = 8;
        let mut triplets = vec![];
        for i in 0..n {
            triplets.push((i, i, 2.0));
            if i > 0 {
                triplets.push((i, i - 1, -1.0));
            }
            if i < n - 1 {
                triplets.push((i, i + 1, -1.0));
            }
        }
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let b = vec![1.0; n];
        let (x, stats) = amg_pcg(&a, &b, 50, 1e-8);
        assert!(
            stats.converged || stats.residual / (n as f64) < 1e-5,
            "AMG-PCG: converged={}, res={}",
            stats.converged,
            stats.residual
        );
        let ax = a.mul_vec(&x);
        let r_norm: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(ai, bi)| (ai - bi).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(r_norm < 1e-4, "AMG-PCG residual too large: {r_norm}");
    }
    #[test]
    fn nested_iteration_solves_diagonal() {
        let triplets = vec![(0, 0, 2.0), (1, 1, 3.0), (2, 2, 4.0)];
        let a = CsrMatrix::from_triplets(3, 3, &triplets);
        let b = [4.0, 9.0, 12.0];
        let (x, stats) = nested_iteration(&a, &b, 2, 50, 50, 1e-12);
        assert!(
            stats.converged || stats.residual < 1e-9,
            "nested iteration: res={}",
            stats.residual
        );
        assert!((x[0] - 2.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-8, "x[1] = {}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-8, "x[2] = {}", x[2]);
    }
    #[test]
    fn nested_iteration_vs_pcg_consistent() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0, 2.0];
        let (x_ni, _) = nested_iteration(&a, &b, 1, 50, 50, 1e-12);
        let precond = jacobi_preconditioner(&a);
        let x_pcg = preconditioned_cg(&a, &b, &precond, 50, 1e-12);
        for i in 0..2 {
            assert!(
                (x_ni[i] - x_pcg[i]).abs() < 1e-7,
                "nested vs PCG differ at {i}: {} vs {}",
                x_ni[i],
                x_pcg[i]
            );
        }
    }
    #[test]
    fn lu_solve_2x2_basic() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let b = vec![1.0, 2.0];
        let x = lu_solve(&a, &b, 2);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-12, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-12, "x[1] = {}", x[1]);
    }
    #[test]
    fn lu_solve_3x3_asymmetric() {
        let a = vec![2.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let b = vec![4.0, 7.0, 4.0];
        let x = lu_solve(&a, &b, 3);
        for i in 0..3 {
            let ax_i: f64 = (0..3).map(|j| a[i * 3 + j] * x[j]).sum();
            assert!(
                (ax_i - b[i]).abs() < 1e-10,
                "row {i}: Ax={ax_i}, b={}",
                b[i]
            );
        }
    }
    #[test]
    fn lu_solve_identity_system() {
        let a = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let b = vec![3.0, 5.0, 7.0];
        let x = lu_solve(&a, &b, 3);
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-12, "x[{i}] = {}", x[i]);
        }
    }
    #[test]
    fn test_femsolve_gmres_2x2_spd() {
        let triplets = vec![(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets);
        let b = [1.0_f64, 2.0];
        let x = fem_gmres(&a, &b, 30, 1e-10);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1] = {}", x[1]);
    }
    #[test]
    fn test_femsolve_gmres_identity_system() {
        let n = 4;
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 1.0)).collect();
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let b = [2.0_f64, 3.0, 5.0, 7.0];
        let x = fem_gmres(&a, &b, 10, 1e-12);
        for i in 0..n {
            assert!((x[i] - b[i]).abs() < 1e-10, "x[{i}] = {}", x[i]);
        }
    }
    #[test]
    fn test_femsolve_gmres_3x3() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, 1.0),
            (1, 0, 1.0),
            (1, 1, 4.0),
            (1, 2, 1.0),
            (2, 1, 1.0),
            (2, 2, 4.0),
        ];
        let a = CsrMatrix::from_triplets(3, 3, &triplets);
        let b = [1.0_f64, 1.0, 1.0];
        let x = fem_gmres(&a, &b, 30, 1e-10);
        let matvec = |row: usize| -> f64 {
            (a.row_ptr[row]..a.row_ptr[row + 1])
                .map(|idx| a.values[idx] * x[a.col_indices[idx]])
                .sum()
        };
        for (i, &bi) in b.iter().enumerate() {
            assert!((matvec(i) - bi).abs() < 1e-8, "row {i}: Ax={}", matvec(i));
        }
    }
    #[test]
    fn test_amg_smooth_output_length() {
        let n = 4;
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
            (2, 3, -1.0),
            (3, 2, -1.0),
            (3, 3, 4.0),
        ];
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let b = [1.0_f64, 0.0, 0.0, 1.0];
        let x0 = vec![0.0_f64; n];
        let x = amg_smooth_step(&a, &b, &x0, 3);
        assert_eq!(x.len(), n, "output length should match n");
    }
    #[test]
    fn test_amg_smooth_reduces_residual() {
        let n = 4;
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
            (2, 3, -1.0),
            (3, 2, -1.0),
            (3, 3, 4.0),
        ];
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let b = [1.0_f64, 0.0, 0.0, 1.0];
        let x0 = vec![0.0_f64; n];
        let r0: f64 = (0..n)
            .map(|i| {
                let ax: f64 = (a.row_ptr[i]..a.row_ptr[i + 1])
                    .map(|idx| a.values[idx] * x0[a.col_indices[idx]])
                    .sum();
                (b[i] - ax).powi(2)
            })
            .sum::<f64>()
            .sqrt();
        let x = amg_smooth_step(&a, &b, &x0, 5);
        let r1: f64 = (0..n)
            .map(|i| {
                let ax: f64 = (a.row_ptr[i]..a.row_ptr[i + 1])
                    .map(|idx| a.values[idx] * x[a.col_indices[idx]])
                    .sum();
                (b[i] - ax).powi(2)
            })
            .sum::<f64>()
            .sqrt();
        assert!(
            r1 <= r0 + 1e-12,
            "AMG smoother should not increase residual: r0={r0}, r1={r1}"
        );
    }
    #[test]
    fn test_amg_smooth_identity_converges() {
        let n = 3;
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 2.0)).collect();
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let b = [4.0_f64, 6.0, 8.0];
        let x0 = vec![0.0_f64; n];
        let x = amg_smooth_step(&a, &b, &x0, 10);
        for (i, (&xi, &bi)) in x.iter().zip(b.iter()).enumerate() {
            assert!((xi - bi / 2.0).abs() < 1e-10, "x[{i}] = {}", xi);
        }
    }
    #[test]
    fn test_ilu0_precond_solve_identity() {
        let n = 3;
        let a = vec![3.0_f64, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 7.0];
        let b = vec![3.0_f64, 5.0, 7.0];
        let x = ilu0_precond_solve(&a, &b, n, 20, 1e-12);
        for (i, &xi) in x.iter().enumerate() {
            assert!((xi - 1.0).abs() < 1e-8, "x[{i}] = {}", xi);
        }
    }
    #[test]
    fn test_ilu0_precond_solve_2x2() {
        let a = vec![4.0_f64, 1.0, 1.0, 3.0];
        let b = vec![1.0_f64, 2.0];
        let x = ilu0_precond_solve(&a, &b, 2, 30, 1e-10);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-6, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-6, "x[1] = {}", x[1]);
    }
    #[test]
    fn test_ilu0_precond_output_length() {
        let n = 5;
        let a: Vec<f64> = (0..n * n)
            .map(|idx| {
                let (i, j) = (idx / n, idx % n);
                if i == j {
                    4.0
                } else if (i as i64 - j as i64).abs() == 1 {
                    -1.0
                } else {
                    0.0
                }
            })
            .collect();
        let b = vec![1.0_f64; n];
        let x = ilu0_precond_solve(&a, &b, n, 50, 1e-10);
        assert_eq!(x.len(), n);
    }
    #[test]
    fn test_ilu0_precond_all_finite() {
        let n = 4;
        let a: Vec<f64> = (0..n * n)
            .map(|idx| {
                let (i, j) = (idx / n, idx % n);
                if i == j {
                    6.0
                } else if (i as i64 - j as i64).abs() == 1 {
                    -1.0
                } else {
                    0.0
                }
            })
            .collect();
        let b = vec![1.0_f64; n];
        let x = ilu0_precond_solve(&a, &b, n, 30, 1e-10);
        for &xi in &x {
            assert!(xi.is_finite(), "x = {xi}");
        }
    }
}
