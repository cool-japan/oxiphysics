//! Auto-generated test module (consolidated from inline `#[cfg(test)] mod` blocks)

use std::f64::consts::PI;

use super::*;

#[cfg(test)]
mod tests_2 {
    use super::*;
    #[test]
    fn test_newton_raphson_sqrt2() {
        let result = newton_raphson(|x| x * x - 2.0, |x| 2.0 * x, 1.5, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 2.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_newton_raphson_cubic() {
        let result = newton_raphson(
            |x| x * x * x - x - 2.0,
            |x| 3.0 * x * x - 1.0,
            2.0,
            1e-12,
            50,
        );
        assert!(result.converged);
        assert!(result.f_val.abs() < 1e-10);
    }
    #[test]
    fn test_secant_method() {
        let result = secant_method(|x| x * x - 9.0, 2.0, 4.0, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_bisect_simple() {
        let result = bisect(|x| x - 5.0, 0.0, 10.0, 1e-10, 100);
        assert!(result.converged);
        assert!((result.root - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_bisect_trigonometric() {
        let result = bisect(|x| x.sin(), 2.5, 3.5, 1e-10, 100);
        assert!(result.converged);
        assert!((result.root - PI).abs() < 1e-9);
    }
    #[test]
    fn test_brentq_converges() {
        let result = brentq(|x| x.exp() - 2.0, 0.0, 2.0, 1e-10, 100);
        assert!(result.converged || result.f_val.abs() < 1e-8);
    }
    #[test]
    fn test_brentq_bad_bracket() {
        let result = brentq(|x| x * x + 1.0, -1.0, 1.0, 1e-10, 50);
        assert!(!result.converged);
    }
    #[test]
    fn test_euler_exponential() {
        let result = euler_ode(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 0.001);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 0.01);
    }
    #[test]
    fn test_rk2_exponential() {
        let result = rk2(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 0.01);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 0.001);
    }
    #[test]
    fn test_rk4_exponential() {
        let result = rk4(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 0.01);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 1e-6);
    }
    #[test]
    fn test_rk4_harmonic_oscillator() {
        let result = rk4(
            |_t, y| vec![y[1], -y[0]],
            0.0,
            vec![0.0, 1.0],
            2.0 * PI,
            0.01,
        );
        let final_val = result.last().unwrap().1[0];
        assert!(final_val.abs() < 0.01);
    }
    #[test]
    fn test_rk45_step() {
        let (y_new, err) = rk45_step(&|_t, y: &[f64]| vec![-y[0]], 0.0, &[1.0], 0.1);
        assert!((y_new[0] - (-0.1_f64).exp()).abs() < 1e-5);
        assert!(err.len() == 1);
    }
    #[test]
    fn test_rk45_adaptive() {
        let result = rk45_adaptive(|_t, y| vec![-y[0]], 0.0, vec![1.0], 1.0, 1e-6, 1e-9);
        let final_val = result.last().unwrap().1[0];
        assert!((final_val - (-1.0_f64).exp()).abs() < 1e-5);
    }
    #[test]
    fn test_bvp_shoot() {
        let path = bvp_shoot(|_x, y, _dy| -PI * PI * y, 0.0, 1.0, 0.0, 0.0, PI, 100, 1e-6);
        assert!(!path.is_empty());
    }
    #[test]
    fn test_fd1d_laplacian_size() {
        let mat = fd1d_laplacian(5, 0.25);
        assert_eq!(mat.len(), 3);
        assert_eq!(mat[0].len(), 3);
    }
    #[test]
    fn test_solve_poisson_1d() {
        let n = 10;
        let dx = 1.0 / (n + 1) as f64;
        let rhs = vec![1.0; n];
        let u = solve_poisson_1d(&rhs, dx);
        let x_mid = 5.0 / (n + 1) as f64;
        let expected = x_mid * (1.0 - x_mid) / 2.0;
        assert!((u[4] - expected).abs() < 0.02);
    }
    #[test]
    fn test_fd2d_laplacian() {
        let n = 5;
        let u: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| (i * j) as f64).collect())
            .collect();
        let lap = fd2d_laplacian(&u, 0.1, 0.1);
        assert_eq!(lap.len(), n);
    }
    #[test]
    fn test_fd3d_laplacian() {
        let n = 4;
        let u = vec![vec![vec![1.0f64; n]; n]; n];
        let lap = fd3d_laplacian(&u, 0.1, 0.1, 0.1);
        assert_eq!(lap[2][2][2], 0.0);
    }
    #[test]
    fn test_deriv_central_sin() {
        let d = deriv_central(|x| x.sin(), 0.0, 1e-5);
        assert!((d - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_deriv_complex_step() {
        let d = deriv_complex_step(|x| x.powi(3), 2.0);
        assert!((d - 12.0).abs() < 1e-6);
    }
    #[test]
    fn test_jacobian_identity() {
        let jac = jacobian(|x| x.to_vec(), &[1.0, 2.0, 3.0], 1e-5);
        for (i, jac_row) in jac.iter().enumerate() {
            for (j, cell) in jac_row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((cell - expected).abs() < 1e-6);
            }
        }
    }
    #[test]
    fn test_gauss_legendre_polynomial() {
        let result = gauss_legendre_integrate(|x| x.powi(3), -1.0, 1.0, 3);
        assert!(result.abs() < 1e-12);
    }
    #[test]
    fn test_gauss_legendre_sin() {
        let result = gauss_legendre_integrate(|x| x.sin(), 0.0, PI, 5);
        assert!((result - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_gauss_hermite_gaussian() {
        let result = gauss_hermite_integrate(|_x| 1.0, 5);
        assert!((result - PI.sqrt()).abs() < 1e-6);
    }
    #[test]
    fn test_adaptive_simpson() {
        let result = adaptive_simpson(&|x: f64| x.sin(), 0.0, PI, 1e-10, 10);
        assert!((result - 2.0).abs() < 1e-8);
    }
    #[test]
    fn test_lu_solve_2x2() {
        let a = vec![vec![2.0, 1.0], vec![5.0, 7.0]];
        let b = vec![11.0, 13.0];
        let x = lu_solve(&a, &b);
        let residual: Vec<f64> = vec![
            a[0][0] * x[0] + a[0][1] * x[1] - b[0],
            a[1][0] * x[0] + a[1][1] * x[1] - b[1],
        ];
        for r in &residual {
            assert!(r.abs() < 1e-10);
        }
    }
    #[test]
    fn test_lu_solve_3x3() {
        let a = vec![
            vec![1.0, 2.0, 3.0],
            vec![0.0, 1.0, 4.0],
            vec![5.0, 6.0, 0.0],
        ];
        let b = vec![14.0, 6.0, 2.0];
        let x = lu_solve(&a, &b);
        for (i, row) in a.iter().enumerate() {
            let sum: f64 = row.iter().zip(x.iter()).map(|(a, x)| a * x).sum();
            assert!((sum - b[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_cholesky_spd() {
        let a = vec![
            vec![4.0, 2.0, 2.0],
            vec![2.0, 5.0, 2.5],
            vec![2.0, 2.5, 3.5],
        ];
        let l = cholesky(&a);
        assert!(l.is_some());
        let l = l.unwrap();
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                let entry: f64 = (0..n).map(|k| l[i][k] * l[j][k]).sum();
                assert!((entry - a[i][j]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_cholesky_not_spd() {
        let a = vec![vec![-1.0, 0.0], vec![0.0, 1.0]];
        assert!(cholesky(&a).is_none());
    }
    #[test]
    fn test_power_iteration() {
        let a = vec![vec![4.0, 1.0], vec![2.0, 3.0]];
        let (lambda, _v) = power_iteration(&a, 1000, 1e-10);
        assert!((lambda - 5.0).abs() < 0.1);
    }
    #[test]
    fn test_qr_eigenvalues_symmetric() {
        let a = vec![
            vec![4.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let eigs = qr_eigenvalues(&a, 100);
        assert_eq!(eigs.len(), 3);
        for e in &eigs {
            assert!(e.is_finite());
        }
    }
    #[test]
    fn test_svd_singular_values() {
        let a = vec![vec![3.0, 0.0], vec![0.0, 2.0]];
        let svs = svd_singular_values(&a, 100);
        assert!(!svs.is_empty());
        assert!(svs[0] > 1.0);
    }
    #[test]
    fn test_continuation() {
        let branch = pseudo_arclength_continuation(|x, lam| x * x - lam, 1.0, 1.0, 0.1, 10, 1e-8);
        assert!(branch.len() > 1);
        for (x, lam) in &branch {
            assert!((x * x - lam).abs() < 0.1);
        }
    }
    #[test]
    fn test_bifurcation_detection() {
        let branch: Vec<(f64, f64)> = (-5..=5)
            .map(|i| {
                let lam = i as f64;
                let x = if lam >= 0.0 { lam.sqrt() } else { 0.0 };
                (x, lam)
            })
            .collect();
        let bifs = detect_bifurcations(|x, lam| x * x - lam, &branch);
        let _ = bifs;
    }
    #[test]
    fn test_matrix_multiply() {
        let a = eye(3);
        let b = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        let c = mat_mul(&a, &b);
        assert_eq!(c, b);
    }
    #[test]
    fn test_frobenius_norm() {
        let a = eye(3);
        assert!((frobenius_norm(&a) - 3.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_transpose() {
        let a = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let at = transpose(&a);
        assert_eq!(at.len(), 3);
        assert_eq!(at[0].len(), 2);
        assert_eq!(at[0][0], 1.0);
        assert_eq!(at[0][1], 4.0);
    }
    #[test]
    fn test_gauss_legendre_nodes_count() {
        for n in 1..=5 {
            let (nodes, weights) = gauss_legendre_nodes(n);
            assert_eq!(nodes.len(), n);
            assert_eq!(weights.len(), n);
        }
    }
    #[test]
    fn test_infinity_norm() {
        let a = vec![vec![1.0, -2.0, 3.0], vec![0.0, 1.0, 0.0]];
        assert!((infinity_norm(&a) - 6.0).abs() < 1e-10);
    }
}
#[cfg(test)]
mod eigen_tests {
    use super::*;
    fn mat2(a: f64, b: f64, c: f64) -> Vec<Vec<f64>> {
        vec![vec![a, b], vec![b, c]]
    }
    #[test]
    fn symmetric_eigen_2x2_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let (evals, evecs) = symmetric_eigen_n(&a);
        assert!((evals[0] - 1.0).abs() < 1e-12);
        assert!((evals[1] - 1.0).abs() < 1e-12);
        let _ = evecs;
    }
    #[test]
    fn symmetric_eigen_2x2_known() {
        let a = mat2(2.0, 1.0, 2.0);
        let (evals, evecs) = symmetric_eigen_n(&a);
        assert!((evals[0] - 1.0).abs() < 1e-10, "eval0={}", evals[0]);
        assert!((evals[1] - 3.0).abs() < 1e-10, "eval1={}", evals[1]);
        let dot = evecs[0][0] * evecs[1][0] + evecs[0][1] * evecs[1][1];
        assert!(dot.abs() < 1e-10, "eigenvectors not orthogonal: dot={dot}");
    }
    #[test]
    fn symmetric_eigen_3x3_diagonal() {
        let a = vec![
            vec![3.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let (mut evals, _) = symmetric_eigen_n(&a);
        evals.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((evals[0] - 1.0).abs() < 1e-12);
        assert!((evals[1] - 2.0).abs() < 1e-12);
        assert!((evals[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn symmetric_eigen_roundtrip() {
        let a = vec![
            vec![4.0, -2.0, 1.0],
            vec![-2.0, 3.0, 0.5],
            vec![1.0, 0.5, 2.0],
        ];
        let (evals, evecs) = symmetric_eigen_n(&a);
        for i in 0..3 {
            let v = &evecs[i];
            let av: Vec<f64> = (0..3)
                .map(|r| (0..3).map(|c| a[r][c] * v[c]).sum::<f64>())
                .collect();
            for r in 0..3 {
                assert!((av[r] - evals[i] * v[r]).abs() < 1e-8, "Av ≠ λv at row {r}");
            }
        }
    }
    #[test]
    fn generalized_eigen_2x2() {
        let f = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        let s = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let (evals, _evecs) = generalized_symmetric_eigen_n(&f, &s).unwrap();
        assert!((evals[0] - 2.0).abs() < 1e-10);
        assert!((evals[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn generalized_eigen_with_overlap() {
        let f = vec![vec![1.5, 0.5], vec![0.5, 1.5]];
        let s = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        let (evals, evecs) = generalized_symmetric_eigen_n(&f, &s).unwrap();
        for i in 0..2 {
            let c = &evecs[i];
            let fc: Vec<f64> = (0..2)
                .map(|r| (0..2).map(|j| f[r][j] * c[j]).sum::<f64>())
                .collect();
            let sc: Vec<f64> = (0..2)
                .map(|r| (0..2).map(|j| s[r][j] * c[j]).sum::<f64>())
                .collect();
            for r in 0..2 {
                assert!(
                    (fc[r] - evals[i] * sc[r]).abs() < 1e-8,
                    "FC ≠ ε*SC at row {r}, evec {i}: fc={} eval*sc={}",
                    fc[r],
                    evals[i] * sc[r]
                );
            }
        }
    }
    #[test]
    fn generalized_eigen_none_on_singular_s() {
        let f = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let s = vec![vec![0.0, 0.0], vec![0.0, 0.0]];
        let result = generalized_symmetric_eigen_n(&f, &s);
        assert!(result.is_none());
    }
}
