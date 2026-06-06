//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{dot, mat_t_vec, mat_vec, norm, prox_l1, solve_spd};
use super::types::{DualDecompResult, PenaltyResult, ProxGradResult};

/// FISTA (Fast ISTA) - accelerated proximal gradient method.
///
/// Same problem as proximal_gradient_l1 but with Nesterov acceleration.
pub fn fista_l1(
    h: &[f64],
    c: &[f64],
    n: usize,
    lambda: f64,
    step_size: f64,
    max_iter: usize,
    tol: f64,
) -> ProxGradResult {
    let mut x = vec![0.0; n];
    let mut y = vec![0.0; n];
    let mut t = 1.0_f64;
    let mut converged = false;
    let mut iterations = 0;
    for iter in 0..max_iter {
        iterations = iter + 1;
        let hy = mat_vec(h, &y, n, n);
        let mut grad = vec![0.0; n];
        for i in 0..n {
            grad[i] = hy[i] + c[i];
        }
        let y_grad: Vec<f64> = (0..n).map(|i| y[i] - step_size * grad[i]).collect();
        let x_new = prox_l1(&y_grad, lambda * step_size);
        let t_new = (1.0 + (1.0 + 4.0 * t * t).sqrt()) / 2.0;
        let momentum = (t - 1.0) / t_new;
        let y_new: Vec<f64> = (0..n)
            .map(|i| x_new[i] + momentum * (x_new[i] - x[i]))
            .collect();
        let diff_norm = norm(&(0..n).map(|i| x_new[i] - x[i]).collect::<Vec<_>>());
        x = x_new;
        y = y_new;
        t = t_new;
        if diff_norm < tol {
            converged = true;
            break;
        }
    }
    let hx = mat_vec(h, &x, n, n);
    let objective = 0.5 * dot(&x, &hx) + dot(c, &x);
    ProxGradResult {
        x,
        objective,
        iterations,
        converged,
    }
}
/// Projected gradient method for box-constrained QP.
///
/// min 0.5 x^T H x + c^T x  s.t. lo <= x <= hi
pub fn projected_gradient_box(
    h: &[f64],
    c: &[f64],
    n: usize,
    lo: &[f64],
    hi: &[f64],
    step_size: f64,
    max_iter: usize,
    tol: f64,
) -> ProxGradResult {
    let mut x: Vec<f64> = (0..n).map(|i| 0.5 * (lo[i] + hi[i])).collect();
    let mut converged = false;
    let mut iterations = 0;
    for iter in 0..max_iter {
        iterations = iter + 1;
        let hx = mat_vec(h, &x, n, n);
        let mut grad = vec![0.0; n];
        for i in 0..n {
            grad[i] = hx[i] + c[i];
        }
        let x_new: Vec<f64> = (0..n)
            .map(|i| (x[i] - step_size * grad[i]).max(lo[i]).min(hi[i]))
            .collect();
        let diff_norm = norm(&(0..n).map(|i| x_new[i] - x[i]).collect::<Vec<_>>());
        x = x_new;
        if diff_norm < tol {
            converged = true;
            break;
        }
    }
    let hx = mat_vec(h, &x, n, n);
    let objective = 0.5 * dot(&x, &hx) + dot(c, &x);
    ProxGradResult {
        x,
        objective,
        iterations,
        converged,
    }
}
/// Dual decomposition for distributed optimization.
///
/// Each subproblem i solves: min 0.5 x_i^T H_i x_i + c_i^T x_i + y^T x_i
///
/// The consensus constraint is: x_i = z for all i.
///
/// `h_list` contains the Hessian for each subproblem (each n x n).
/// `c_list` contains the linear cost for each subproblem (each length n).
pub fn dual_decomposition(
    h_list: &[Vec<f64>],
    c_list: &[Vec<f64>],
    n: usize,
    step_size: f64,
    max_iter: usize,
    tol: f64,
) -> DualDecompResult {
    let num_sub = h_list.len();
    let mut y = vec![0.0; n];
    let mut x_locals: Vec<Vec<f64>> = vec![vec![0.0; n]; num_sub];
    let mut z = vec![0.0; n];
    let mut converged = false;
    let mut iterations = 0;
    for iter in 0..max_iter {
        iterations = iter + 1;
        for i in 0..num_sub {
            let mut rhs = vec![0.0; n];
            for j in 0..n {
                rhs[j] = -(c_list[i][j] + y[j]);
            }
            x_locals[i] = match solve_spd(&h_list[i], &rhs, n) {
                Some(s) => s,
                None => x_locals[i].clone(),
            };
        }
        let inv_num = 1.0 / num_sub as f64;
        for (j, z_j) in z.iter_mut().enumerate() {
            *z_j = x_locals.iter().map(|xl| xl[j]).sum::<f64>() * inv_num;
        }
        let mut residual = vec![0.0; n];
        for xl in x_locals.iter() {
            for (j, (res_j, z_j)) in residual.iter_mut().zip(z.iter()).enumerate() {
                *res_j += (xl[j] - z_j).abs();
            }
        }
        let res_norm = norm(&residual);
        for (j, y_j) in y.iter_mut().enumerate() {
            let avg_diff: f64 = x_locals.iter().map(|xl| xl[j] - z[j]).sum::<f64>();
            *y_j += step_size * avg_diff / num_sub as f64;
        }
        if res_norm < tol {
            converged = true;
            break;
        }
    }
    DualDecompResult {
        x_locals,
        z,
        dual: y,
        iterations,
        converged,
    }
}
/// Quadratic penalty method for equality constraints.
///
/// min 0.5 x^T H x + c^T x + (rho/2) ||A x - b||^2
pub fn penalty_method(
    h: &[f64],
    c: &[f64],
    a_eq: &[f64],
    b_eq: &[f64],
    n: usize,
    m: usize,
    rho_init: f64,
    rho_factor: f64,
    outer_iter: usize,
    tol: f64,
) -> PenaltyResult {
    let mut x = vec![0.0; n];
    let mut rho = rho_init;
    let mut converged = false;
    let mut violation = f64::MAX;
    for _outer in 0..outer_iter {
        let mut h_pen = h.to_vec();
        for i in 0..n {
            for j in 0..n {
                for k in 0..m {
                    h_pen[i * n + j] += rho * a_eq[k * n + i] * a_eq[k * n + j];
                }
            }
        }
        let at_b = mat_t_vec(a_eq, b_eq, m, n);
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            rhs[i] = -c[i] + rho * at_b[i];
        }
        x = match solve_spd(&h_pen, &rhs, n) {
            Some(s) => s,
            None => x,
        };
        let ax = mat_vec(a_eq, &x, m, n);
        let mut res = vec![0.0; m];
        for i in 0..m {
            res[i] = ax[i] - b_eq[i];
        }
        violation = norm(&res);
        if violation < tol {
            converged = true;
            break;
        }
        rho *= rho_factor;
    }
    let hx = mat_vec(h, &x, n, n);
    let objective = 0.5 * dot(&x, &hx) + dot(c, &x);
    PenaltyResult {
        x,
        objective,
        penalty: rho,
        violation,
        converged,
    }
}
/// Build a symmetric matrix from its upper triangle (row-major).
pub fn symmetric_from_upper(upper: &[f64], n: usize) -> Vec<f64> {
    let mut full = vec![0.0; n * n];
    let mut idx = 0;
    for i in 0..n {
        for j in i..n {
            full[i * n + j] = upper[idx];
            full[j * n + i] = upper[idx];
            idx += 1;
        }
    }
    full
}
/// Build a diagonal matrix.
pub fn diagonal_matrix(diag: &[f64]) -> Vec<f64> {
    let n = diag.len();
    let mut mat = vec![0.0; n * n];
    for i in 0..n {
        mat[i * n + i] = diag[i];
    }
    mat
}
/// Identity matrix of size n.
pub fn identity_matrix(n: usize) -> Vec<f64> {
    let mut mat = vec![0.0; n * n];
    for i in 0..n {
        mat[i * n + i] = 1.0;
    }
    mat
}
#[cfg(test)]
mod tests {
    use super::super::*;
    pub(super) const TOL: f64 = 1e-6;
    #[test]
    fn test_dot_product() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot(&a, &b) - 32.0).abs() < TOL);
    }
    #[test]
    fn test_norm() {
        let v = [3.0, 4.0];
        assert!((norm(&v) - 5.0).abs() < TOL);
    }
    #[test]
    fn test_mat_vec_identity() {
        let id = [1.0, 0.0, 0.0, 1.0];
        let x = [3.0, 7.0];
        let y = mat_vec(&id, &x, 2, 2);
        assert!((y[0] - 3.0).abs() < TOL);
        assert!((y[1] - 7.0).abs() < TOL);
    }
    #[test]
    fn test_solve_spd_2x2() {
        let h = [4.0, 2.0, 2.0, 3.0];
        let b = [1.0, 2.0];
        let x = solve_spd(&h, &b, 2).unwrap();
        let hx = mat_vec(&h, &x, 2, 2);
        assert!((hx[0] - b[0]).abs() < TOL);
        assert!((hx[1] - b[1]).abs() < TOL);
    }
    #[test]
    fn test_solve_spd_not_spd() {
        let h = [-1.0, 0.0, 0.0, -1.0];
        let b = [1.0, 1.0];
        assert!(solve_spd(&h, &b, 2).is_none());
    }
    #[test]
    fn test_qp_unconstrained() {
        let h = [1.0, 0.0, 0.0, 1.0];
        let c = [-1.0, -2.0];
        let a: Vec<f64> = vec![];
        let b: Vec<f64> = vec![];
        let result = qp_active_set(&h, &c, &a, &b, 2, 0, 100);
        assert!((result.x[0] - 1.0).abs() < TOL);
        assert!((result.x[1] - 2.0).abs() < TOL);
    }
    #[test]
    fn test_qp_single_constraint() {
        let h = [1.0, 0.0, 0.0, 1.0];
        let c = [-2.0, -2.0];
        let a = [1.0, 1.0];
        let b = [2.0];
        let result = qp_active_set(&h, &c, &a, &b, 2, 1, 100);
        assert!((result.x[0] - 1.0).abs() < 0.1, "x[0]={}", result.x[0]);
        assert!((result.x[1] - 1.0).abs() < 0.1, "x[1]={}", result.x[1]);
    }
    #[test]
    fn test_lcp_trivial() {
        let m_mat = [1.0, 0.0, 0.0, 1.0];
        let q = [1.0, 2.0];
        let result = lcp_lemke(&m_mat, &q, 2, 100);
        assert!(result.solved);
        assert!(result.z[0].abs() < TOL);
        assert!(result.z[1].abs() < TOL);
    }
    #[test]
    fn test_lcp_simple() {
        let m_mat = [1.0, 0.0, 0.0, 1.0];
        let q = [-1.0, -1.0];
        let result = lcp_lemke(&m_mat, &q, 2, 100);
        assert!(result.solved);
        assert!((result.z[0] - 1.0).abs() < TOL, "z[0]={}", result.z[0]);
        assert!((result.z[1] - 1.0).abs() < TOL, "z[1]={}", result.z[1]);
    }
    #[test]
    fn test_lcp_empty() {
        let result = lcp_lemke(&[], &[], 0, 100);
        assert!(result.solved);
        assert!(result.z.is_empty());
    }
    #[test]
    fn test_lcp_mixed() {
        let m_mat = [1.0, 0.0, 0.0, 1.0];
        let q = [-2.0, 1.0];
        let result = lcp_lemke(&m_mat, &q, 2, 100);
        assert!(result.solved);
        assert!((result.z[0] - 2.0).abs() < TOL, "z[0]={}", result.z[0]);
        assert!(result.z[1].abs() < TOL, "z[1]={}", result.z[1]);
    }
    #[test]
    fn test_barrier_simple() {
        let f = [1.0];
        let a = [-1.0];
        let b = [-1.0];
        let x0 = [2.0];
        let result = barrier_method(BarrierMethodParams {
            f: &f,
            a_ineq: &a,
            b_ineq: &b,
            n: 1,
            m: 1,
            x0: &x0,
            mu_init: 1.0,
            mu_factor: 0.5,
            outer_iter: 50,
            inner_iter: 50,
            tol: 1e-6,
        });
        assert!((result.x[0] - 1.0).abs() < 0.1, "x={}", result.x[0]);
    }
    #[test]
    fn test_aug_lag_equality() {
        let h = [1.0, 0.0, 0.0, 1.0];
        let c = [0.0, 0.0];
        let a_eq = [1.0, 1.0];
        let b_eq = [1.0];
        let result = augmented_lagrangian(AugLagParams {
            h: &h,
            c: &c,
            a_eq: &a_eq,
            b_eq: &b_eq,
            n: 2,
            m: 1,
            rho_init: 1.0,
            rho_factor: 2.0,
            outer_iter: 50,
            inner_iter: 50,
            tol: 1e-6,
        });
        assert!(result.converged, "violation={}", result.violation);
        assert!((result.x[0] - 0.5).abs() < 0.01, "x[0]={}", result.x[0]);
        assert!((result.x[1] - 0.5).abs() < 0.01, "x[1]={}", result.x[1]);
    }
    #[test]
    fn test_aug_lag_single_eq() {
        let h = [2.0];
        let c = [0.0];
        let a_eq = [1.0];
        let b_eq = [3.0];
        let result = augmented_lagrangian(AugLagParams {
            h: &h,
            c: &c,
            a_eq: &a_eq,
            b_eq: &b_eq,
            n: 1,
            m: 1,
            rho_init: 1.0,
            rho_factor: 2.0,
            outer_iter: 50,
            inner_iter: 50,
            tol: 1e-6,
        });
        assert!((result.x[0] - 3.0).abs() < 0.01, "x={}", result.x[0]);
    }
    #[test]
    fn test_admm_unconstrained() {
        let h = [1.0, 0.0, 0.0, 1.0];
        let c = [-1.0, -2.0];
        let a: Vec<f64> = vec![];
        let b: Vec<f64> = vec![];
        let result = admm_consensus(&h, &c, &a, &b, 2, 0, 1.0, 200, 1e-6);
        assert!((result.x[0] - 1.0).abs() < 0.1, "x[0]={}", result.x[0]);
        assert!((result.x[1] - 2.0).abs() < 0.1, "x[1]={}", result.x[1]);
    }
    #[test]
    fn test_prox_l1_shrink() {
        let v = [3.0, -2.0, 0.5];
        let result = prox_l1(&v, 1.0);
        assert!((result[0] - 2.0).abs() < TOL);
        assert!((result[1] + 1.0).abs() < TOL);
        assert!(result[2].abs() < TOL);
    }
    #[test]
    fn test_prox_l2_shrink() {
        let v = [3.0, 4.0];
        let result = prox_l2(&v, 2.0);
        assert!((result[0] - 1.8).abs() < TOL);
        assert!((result[1] - 2.4).abs() < TOL);
    }
    #[test]
    fn test_prox_l2_zero() {
        let v = [1.0, 0.0];
        let result = prox_l2(&v, 2.0);
        assert!(result[0].abs() < TOL);
        assert!(result[1].abs() < TOL);
    }
    #[test]
    fn test_prox_nonneg() {
        let v = [1.0, -2.0, 0.0, 3.0];
        let result = prox_nonneg(&v);
        assert!((result[0] - 1.0).abs() < TOL);
        assert!(result[1].abs() < TOL);
        assert!(result[2].abs() < TOL);
        assert!((result[3] - 3.0).abs() < TOL);
    }
    #[test]
    fn test_prox_box() {
        let v = [5.0, -3.0, 0.5];
        let result = prox_box(&v, 0.0, 1.0);
        assert!((result[0] - 1.0).abs() < TOL);
        assert!(result[1].abs() < TOL);
        assert!((result[2] - 0.5).abs() < TOL);
    }
    #[test]
    fn test_prox_soc_inside() {
        let v = [3.0, 1.0, 1.0];
        let result = prox_soc(&v);
        assert!((result[0] - 3.0).abs() < TOL);
        assert!((result[1] - 1.0).abs() < TOL);
    }
    #[test]
    fn test_prox_soc_outside() {
        let v = [0.0, 3.0, 4.0];
        let result = prox_soc(&v);
        assert!((result[0] - 2.5).abs() < TOL, "t={}", result[0]);
    }
    #[test]
    fn test_prox_elastic_net() {
        let v = [5.0, -3.0, 0.1];
        let result = prox_elastic_net(&v, 1.0, 0.5);
        assert!((result[0] - 3.0).abs() < TOL);
    }
    #[test]
    fn test_moreau_l1() {
        assert!(moreau_envelope_l1(&[0.0], 1.0).abs() < TOL);
        assert!((moreau_envelope_l1(&[2.0], 1.0) - 1.5).abs() < TOL);
    }
    #[test]
    fn test_licq_independent() {
        let grads = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(check_licq(&grads, 2, 1e-10));
    }
    #[test]
    fn test_licq_dependent() {
        let grads = vec![vec![1.0, 0.0], vec![2.0, 0.0]];
        assert!(!check_licq(&grads, 2, 1e-10));
    }
    #[test]
    fn test_licq_empty() {
        let grads: Vec<Vec<f64>> = vec![];
        assert!(check_licq(&grads, 3, 1e-10));
    }
    #[test]
    fn test_mfcq_simple() {
        let eq_grads: Vec<Vec<f64>> = vec![];
        let ineq_grads = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(check_mfcq(&eq_grads, &ineq_grads, 2, 1e-10));
    }
    #[test]
    fn test_kkt_satisfied() {
        let h = [1.0];
        let c = [-1.0];
        let a = [1.0];
        let b = [0.5];
        let x = [0.5];
        let mu = [0.5];
        let kkt = check_kkt(&h, &c, &a, &b, &x, &mu, 1, 1, 1e-6);
        assert!(
            kkt.satisfied,
            "stationarity={}, pf={}, df={}, cs={}",
            kkt.stationarity, kkt.primal_feasibility, kkt.dual_feasibility, kkt.complementarity
        );
    }
    #[test]
    fn test_kkt_not_satisfied() {
        let h = [1.0];
        let c = [-1.0];
        let a = [1.0];
        let b = [0.5];
        let x = [0.0];
        let mu = [0.0];
        let kkt = check_kkt(&h, &c, &a, &b, &x, &mu, 1, 1, 1e-6);
        assert!(!kkt.satisfied);
    }
    #[test]
    fn test_proximal_gradient_l1_sparse() {
        let h = identity_matrix(3);
        let c = [0.0, 0.0, 0.0];
        let result = proximal_gradient_l1(&h, &c, 3, 10.0, 0.5, 100, 1e-6);
        for i in 0..3 {
            assert!(result.x[i].abs() < TOL, "x[{i}]={}", result.x[i]);
        }
    }
    #[test]
    fn test_fista_convergence() {
        let h = [1.0];
        let c = [-3.0];
        let result = fista_l1(&h, &c, 1, 0.1, 0.5, 500, 1e-6);
        assert!((result.x[0] - 2.9).abs() < 0.1, "x={}", result.x[0]);
    }
    #[test]
    fn test_projected_gradient_box() {
        let h = [1.0];
        let c = [-5.0];
        let lo = [0.0];
        let hi = [3.0];
        let result = projected_gradient_box(&h, &c, 1, &lo, &hi, 0.5, 200, 1e-6);
        assert!((result.x[0] - 3.0).abs() < 0.01, "x={}", result.x[0]);
    }
    #[test]
    fn test_dual_decomposition_consensus() {
        let h1 = vec![1.0, 0.0, 0.0, 1.0];
        let h2 = vec![1.0, 0.0, 0.0, 1.0];
        let c1 = vec![-1.0, 0.0];
        let c2 = vec![0.0, -1.0];
        let result = dual_decomposition(&[h1, h2], &[c1, c2], 2, 0.1, 200, 1e-4);
        assert!((result.z[0] - 0.5).abs() < 0.2, "z[0]={}", result.z[0]);
        assert!((result.z[1] - 0.5).abs() < 0.2, "z[1]={}", result.z[1]);
    }
    #[test]
    fn test_penalty_method_equality() {
        let h = [1.0];
        let c = [0.0];
        let a_eq = [1.0];
        let b_eq = [2.0];
        let result = penalty_method(&h, &c, &a_eq, &b_eq, 1, 1, 1.0, 10.0, 20, 1e-6);
        assert!((result.x[0] - 2.0).abs() < 0.01, "x={}", result.x[0]);
    }
    #[test]
    fn test_symmetric_from_upper() {
        let upper = [1.0, 2.0, 3.0];
        let full = symmetric_from_upper(&upper, 2);
        assert!((full[0] - 1.0).abs() < TOL);
        assert!((full[1] - 2.0).abs() < TOL);
        assert!((full[2] - 2.0).abs() < TOL);
        assert!((full[3] - 3.0).abs() < TOL);
    }
    #[test]
    fn test_diagonal_matrix() {
        let d = diagonal_matrix(&[2.0, 3.0]);
        assert!((d[0] - 2.0).abs() < TOL);
        assert!(d[1].abs() < TOL);
        assert!(d[2].abs() < TOL);
        assert!((d[3] - 3.0).abs() < TOL);
    }
    #[test]
    fn test_identity_matrix() {
        let id = identity_matrix(3);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((id[i * 3 + j] - expected).abs() < TOL);
            }
        }
    }
}
