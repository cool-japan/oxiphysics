//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use crate::math::{Mat3, Real, Vec3};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::*;
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
        for (i, lu_row) in lu.iter().enumerate() {
            assert!((lu_row[i] - 1.0).abs() < 1e-10, "LU diagonal[{i}] != 1");
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
        for (i, r_row) in r.iter().enumerate() {
            for (j, cell) in r_row[..i].iter().enumerate() {
                assert!(
                    cell.abs() < 1e-8,
                    "R should be upper triangular: R[{i}][{j}] = {}",
                    cell
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
        for (i, ata_row) in ata.iter().enumerate() {
            for (j, cell) in ata_row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (cell - expected).abs() < 1e-10,
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
        for (i, l_row) in l.iter().enumerate() {
            assert!((l_row[i] - 1.0).abs() < 1e-12, "L diagonal = 1");
            for cell in l_row[(i + 1)..].iter() {
                assert!(cell.abs() < 1e-12, "L upper triangular part should be 0");
            }
        }
        for (i, u_row) in u.iter().enumerate() {
            for (j, cell) in u_row[..i].iter().enumerate() {
                assert!(
                    cell.abs() < 1e-12,
                    "U lower triangular part should be 0 at [{i}][{j}]"
                );
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
