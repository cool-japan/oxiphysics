//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::CsrMatrix;

/// Extract the upper triangular part of a sparse matrix (including diagonal).
pub fn upper_triangular(a: &CsrMatrix) -> CsrMatrix {
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for i in 0..a.nrows {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let j = a.col_idx[k];
            if j >= i {
                rows.push(i);
                cols.push(j);
                vals.push(a.values[k]);
            }
        }
    }
    CsrMatrix::from_triplets(a.nrows, a.ncols, &rows, &cols, &vals)
}
/// Symmetrise a sparse matrix: A_sym = (A + A^T) / 2.
pub fn symmetrise_sparse(a: &CsrMatrix) -> CsrMatrix {
    let at = a.transpose();
    let sum = a.add(&at);
    sum.scale(0.5)
}
/// Compute the 1-norm (max column sum of absolute values) of a sparse matrix.
pub fn sparse_one_norm(a: &CsrMatrix) -> f64 {
    let mut col_sums = vec![0.0f64; a.ncols];
    for (&col, &val) in a.col_idx.iter().zip(a.values.iter()) {
        col_sums[col] += val.abs();
    }
    col_sums.iter().cloned().fold(0.0f64, f64::max)
}
/// Compute the infinity-norm (max row sum of absolute values) of a sparse matrix.
pub fn sparse_inf_norm(a: &CsrMatrix) -> f64 {
    (0..a.nrows)
        .map(|i| {
            (a.row_ptr[i]..a.row_ptr[i + 1])
                .map(|k| a.values[k].abs())
                .sum::<f64>()
        })
        .fold(0.0f64, f64::max)
}
/// Compute the Frobenius norm of a sparse matrix.
pub fn sparse_frobenius_norm(a: &CsrMatrix) -> f64 {
    a.values.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Set the diagonal of a CSR matrix (in-place). Only existing diagonal entries are modified.
pub fn set_diagonal(a: &mut CsrMatrix, d: &[f64]) {
    assert_eq!(d.len(), a.nrows.min(a.ncols));
    for (i, &di) in d.iter().enumerate() {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            if a.col_idx[k] == i {
                a.values[k] = di;
                break;
            }
        }
    }
}
/// Add a scalar multiple of the identity to a CSR matrix: A + alpha * I.
/// Requires the diagonal entries to already exist in the sparsity pattern.
pub fn add_identity_scaled(a: &CsrMatrix, alpha: f64) -> CsrMatrix {
    let n = a.nrows.min(a.ncols);
    let mut new_vals = a.values.clone();
    for i in 0..n {
        let row_start = a.row_ptr[i];
        let row_end = a.row_ptr[i + 1];
        for (k_off, (&col, nv)) in a.col_idx[row_start..row_end]
            .iter()
            .zip(new_vals[row_start..row_end].iter_mut())
            .enumerate()
        {
            if col == i {
                *nv += alpha;
                let _ = k_off;
                break;
            }
        }
    }
    CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values: new_vals,
    }
}
#[cfg(test)]
mod tests_extended {
    use super::super::*;
    use crate::sparse::BlockJacobi;
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    #[test]
    fn test_bicgstab_identity_system() {
        let a = identity_csr(5);
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (x, conv, _) = bicgstab_solve(&a, &b, None, 1e-10, 100);
        assert!(conv, "should converge");
        for i in 0..5 {
            assert!(approx(x[i], b[i], 1e-6));
        }
    }
    #[test]
    fn test_bicgstab_laplacian() {
        let a = laplacian_1d(5, 1.0);
        let x_true = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        let b = a.matvec(&x_true);
        let (x, conv, _) = bicgstab_solve(&a, &b, None, 1e-8, 200);
        assert!(conv, "BICGStab should converge on 1D Laplacian");
        for i in 0..5 {
            assert!(approx(x[i], x_true[i], 1e-4), "x[{i}]={}", x[i]);
        }
    }
    #[test]
    fn test_bicgstab_zero_rhs() {
        let a = laplacian_1d(4, 1.0);
        let b = vec![0.0; 4];
        let (x, conv, _) = bicgstab_solve(&a, &b, None, 1e-10, 50);
        assert!(conv);
        for &xi in x.iter() {
            assert!(approx(xi, 0.0, 1e-10));
        }
    }
    #[test]
    fn test_bicgstab_initial_guess() {
        let a = identity_csr(3);
        let b = vec![4.0, 5.0, 6.0];
        let x0 = vec![3.9, 4.9, 5.9];
        let (x, conv, iters) = bicgstab_solve(&a, &b, Some(&x0), 1e-10, 50);
        assert!(conv, "should converge from near-exact initial guess");
        for (&xi, &bi) in x.iter().zip(b.iter()) {
            assert!(approx(xi, bi, 1e-6));
        }
        assert!(iters < 10, "iters={iters}");
    }
    #[test]
    fn test_gauss_seidel_diagonal_system() {
        let a = identity_csr(4);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let (x, iters) = gauss_seidel_solve(&a, &b, 1e-10, 100);
        assert!(iters <= 3, "iters={iters}");
        for i in 0..4 {
            assert!(approx(x[i], b[i], 1e-10));
        }
    }
    #[test]
    fn test_gauss_seidel_converges_laplacian() {
        let a = laplacian_1d(5, 1.0);
        let a2 = add_identity_scaled(&a, 4.0);
        let x_true = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = a2.matvec(&x_true);
        let (x, _iters) = gauss_seidel_solve(&a2, &b, 1e-8, 500);
        for i in 0..5 {
            assert!(approx(x[i], x_true[i], 1e-4), "x[{i}]={}", x[i]);
        }
    }
    #[test]
    fn test_forward_substitution_identity() {
        let l = identity_csr(4);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x = forward_substitution(&l, &b);
        for i in 0..4 {
            assert!(approx(x[i], b[i], 1e-12));
        }
    }
    #[test]
    fn test_backward_substitution_identity() {
        let u = identity_csr(4);
        let b = vec![4.0, 3.0, 2.0, 1.0];
        let x = backward_substitution(&u, &b);
        for i in 0..4 {
            assert!(approx(x[i], b[i], 1e-12));
        }
    }
    #[test]
    fn test_forward_substitution_lower_triangular() {
        let l = CsrMatrix::from_triplets(2, 2, &[0, 1, 1], &[0, 0, 1], &[2.0, 3.0, 4.0]);
        let b = vec![4.0, 11.0];
        let x = forward_substitution(&l, &b);
        assert!(approx(x[0], 2.0, 1e-10));
        assert!(approx(x[1], 1.25, 1e-10));
    }
    #[test]
    fn test_backward_substitution_upper_triangular() {
        let u = CsrMatrix::from_triplets(2, 2, &[0, 0, 1], &[0, 1, 1], &[2.0, 3.0, 4.0]);
        let b = vec![11.0, 8.0];
        let x = backward_substitution(&u, &b);
        assert!(approx(x[1], 2.0, 1e-10));
        assert!(approx(x[0], 2.5, 1e-10));
    }
    #[test]
    fn test_banded_matrix_tridiagonal() {
        let a = banded_matrix(4, &[-1, 0, 1], &[-1.0, 2.0, -1.0]);
        assert!(approx(a.get(0, 0), 2.0, 1e-12));
        assert!(approx(a.get(0, 1), -1.0, 1e-12));
        assert!(approx(a.get(1, 0), -1.0, 1e-12));
        assert!(approx(a.get(3, 3), 2.0, 1e-12));
        assert!(approx(a.get(0, 3), 0.0, 1e-12));
    }
    #[test]
    fn test_banded_matrix_diagonal_only() {
        let a = banded_matrix(3, &[0], &[5.0]);
        for i in 0..3 {
            assert!(approx(a.get(i, i), 5.0, 1e-12));
        }
        assert!(approx(a.get(0, 1), 0.0, 1e-12));
    }
    #[test]
    fn test_row_scale_identity() {
        let a = identity_csr(3);
        let s = vec![2.0, 3.0, 5.0];
        let b = row_scale(&a, &s);
        for (i, &si) in s.iter().enumerate() {
            assert!(approx(b.get(i, i), si, 1e-12));
        }
    }
    #[test]
    fn test_col_scale_identity() {
        let a = identity_csr(3);
        let s = vec![7.0, 11.0, 13.0];
        let b = col_scale(&a, &s);
        for (i, &si) in s.iter().enumerate() {
            assert!(approx(b.get(i, i), si, 1e-12));
        }
    }
    #[test]
    fn test_equilibration_scales_identity() {
        let a = identity_csr(4);
        let (rs, cs) = equilibration_scales(&a);
        for &r in rs.iter() {
            assert!(approx(r, 1.0, 1e-10));
        }
        for &c in cs.iter() {
            assert!(approx(c, 1.0, 1e-10));
        }
    }
    #[test]
    fn test_multi_rhs_cg_identity() {
        let a = identity_csr(3);
        let b1 = vec![1.0, 2.0, 3.0];
        let b2 = vec![4.0, 5.0, 6.0];
        let xs = multi_rhs_cg(&a, &[b1.clone(), b2.clone()], 1e-10, 50);
        for (&xi, &bi) in xs[0].iter().zip(b1.iter()) {
            assert!(approx(xi, bi, 1e-6));
        }
        for (&xi, &bi) in xs[1].iter().zip(b2.iter()) {
            assert!(approx(xi, bi, 1e-6));
        }
    }
    #[test]
    fn test_lower_triangular_extraction() {
        let a = laplacian_1d(4, 1.0);
        let l = lower_triangular(&a);
        for i in 0..4 {
            for j in (i + 1)..4 {
                assert!(
                    approx(l.get(i, j), 0.0, 1e-12),
                    "L[{i},{j}]={}",
                    l.get(i, j)
                );
            }
        }
    }
    #[test]
    fn test_upper_triangular_extraction() {
        let a = laplacian_1d(4, 1.0);
        let u = upper_triangular(&a);
        for i in 1..4 {
            for j in 0..i {
                assert!(
                    approx(u.get(i, j), 0.0, 1e-12),
                    "U[{i},{j}]={}",
                    u.get(i, j)
                );
            }
        }
    }
    #[test]
    fn test_symmetrise_sparse_already_symmetric() {
        let a = laplacian_1d(4, 1.0);
        let asym = symmetrise_sparse(&a);
        for i in 0..4 {
            for j in 0..4 {
                assert!(approx(asym.get(i, j), a.get(i, j), 1e-10));
            }
        }
    }
    #[test]
    fn test_symmetrise_sparse_asymmetric_input() {
        let a = CsrMatrix::from_triplets(2, 2, &[0, 0, 1], &[0, 1, 0], &[1.0, 2.0, 0.0]);
        let asym = symmetrise_sparse(&a);
        assert!(approx(asym.get(0, 1), 1.0, 1e-10));
        assert!(approx(asym.get(1, 0), 1.0, 1e-10));
    }
    #[test]
    fn test_sparse_one_norm_identity() {
        let a = identity_csr(5);
        assert!(approx(sparse_one_norm(&a), 1.0, 1e-12));
    }
    #[test]
    fn test_sparse_inf_norm_identity() {
        let a = identity_csr(5);
        assert!(approx(sparse_inf_norm(&a), 1.0, 1e-12));
    }
    #[test]
    fn test_sparse_frobenius_norm_identity() {
        let a = identity_csr(5);
        assert!(approx(sparse_frobenius_norm(&a), (5.0f64).sqrt(), 1e-12));
    }
    #[test]
    fn test_sparse_one_norm_laplacian() {
        let a = laplacian_1d(4, 1.0);
        assert!(approx(sparse_one_norm(&a), 4.0, 1e-10));
    }
    #[test]
    fn test_add_identity_scaled_zero() {
        let a = laplacian_1d(3, 1.0);
        let b = add_identity_scaled(&a, 0.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(approx(b.get(i, j), a.get(i, j), 1e-12));
            }
        }
    }
    #[test]
    fn test_add_identity_scaled_positive() {
        let a = identity_csr(3);
        let b = add_identity_scaled(&a, 2.0);
        for i in 0..3 {
            assert!(approx(b.get(i, i), 3.0, 1e-12));
        }
    }
    #[test]
    fn test_block_jacobi_identity() {
        let a = identity_csr(4);
        let bj = BlockJacobi::new(&a, 2);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x = bj.apply(&b);
        for i in 0..4 {
            assert!(approx(x[i], b[i], 1e-10));
        }
    }
    #[test]
    fn test_block_jacobi_diagonal_matrix() {
        let a = CsrMatrix::from_triplets(4, 4, &[0, 1, 2, 3], &[0, 1, 2, 3], &[2.0, 4.0, 3.0, 6.0]);
        let bj = BlockJacobi::new(&a, 2);
        let b = vec![2.0, 4.0, 3.0, 6.0];
        let x = bj.apply(&b);
        for &xi in x.iter() {
            assert!(approx(xi, 1.0, 1e-10));
        }
    }
    #[test]
    fn test_lanczos_identity_eigenvalue() {
        let a = identity_csr(5);
        let b0 = vec![1.0, 0.0, 0.0, 0.0, 0.0];
        let (alphas, _betas, _q) = lanczos(&a, &b0, 3);
        assert!(approx(alphas[0], 1.0, 1e-10));
    }
    #[test]
    fn test_lanczos_returns_k_alphas() {
        let a = laplacian_1d(6, 1.0);
        let b0 = vec![1.0; 6];
        let (alphas, betas, q) = lanczos(&a, &b0, 4);
        assert!(alphas.len() <= 4);
        assert!(betas.len() <= 4);
        assert!(!q.is_empty());
    }
    #[test]
    fn test_lanczos_vectors_orthogonal() {
        let a = laplacian_1d(6, 1.0);
        let b0: Vec<f64> = (0..6).map(|i| (i + 1) as f64).collect();
        let (_alphas, _betas, q) = lanczos(&a, &b0, 4);
        for i in 0..q.len() {
            for j in (i + 1)..q.len() {
                let dot: f64 = q[i].iter().zip(q[j].iter()).map(|(a, b)| a * b).sum();
                assert!(dot.abs() < 1e-8, "q[{i}]·q[{j}] = {dot}");
            }
        }
    }
    #[test]
    fn test_ritz_values_1x1() {
        let alphas = vec![3.125];
        let betas: Vec<f64> = vec![];
        let rv = ritz_values(&alphas, &betas);
        assert_eq!(rv.len(), 1);
        assert!(approx(rv[0], 3.125, 1e-6));
    }
    #[test]
    fn test_ritz_values_count() {
        let alphas = vec![2.0, 5.0, 7.0];
        let betas = vec![0.5, 0.3];
        let rv = ritz_values(&alphas, &betas);
        assert_eq!(rv.len(), 3, "should return 3 Ritz values");
        for v in &rv {
            assert!(v.is_finite(), "Ritz value should be finite: {v}");
        }
    }
    #[test]
    fn test_minres_identity() {
        let a = identity_csr(4);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let (x, conv, _) = minres_solve(&a, &b, 1e-8, 100);
        assert!(
            conv || {
                let ax = a.matvec(&x);
                let res: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
                res.iter().map(|x| x * x).sum::<f64>().sqrt() < 1e-4
            }
        );
    }
    #[test]
    fn test_minres_zero_rhs() {
        let a = identity_csr(3);
        let b = vec![0.0, 0.0, 0.0];
        let (x, _conv, _) = minres_solve(&a, &b, 1e-10, 50);
        for &xi in x.iter() {
            assert!(xi.abs() < 1e-10);
        }
    }
}
