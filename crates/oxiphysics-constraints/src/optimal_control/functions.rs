//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Multiply two n×n matrices stored row-major. Returns C = A * B.
pub(super) fn mat_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            for j in 0..n {
                c[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    c
}
/// Add two n×n matrices: C = A + B.
pub(super) fn mat_add(a: &[f64], b: &[f64], _n: usize) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
/// Subtract: C = A - B.
pub(super) fn mat_sub(a: &[f64], b: &[f64], _n: usize) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
/// Transpose n×n matrix.
pub(super) fn mat_transpose(a: &[f64], n: usize) -> Vec<f64> {
    let mut t = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            t[j * n + i] = a[i * n + j];
        }
    }
    t
}
/// Identity n×n matrix.
pub(super) fn mat_eye(n: usize) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        m[i * n + i] = 1.0;
    }
    m
}
/// Matrix-vector product: y = A * x  (A is n×n, x is n).
pub(super) fn mat_vec(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            y[i] += a[i * n + j] * x[j];
        }
    }
    y
}
/// Dot product of two vectors.
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Frobenius norm of a matrix or vector.
pub(super) fn frob_norm(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Invert a general n×n matrix via Gaussian elimination with partial pivoting.
/// Returns `None` if the matrix is singular to within tolerance.
pub(super) fn mat_inv(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            aug[i * 2 * n + j] = a[i * n + j];
        }
        aug[i * 2 * n + n + i] = 1.0;
    }
    let w = 2 * n;
    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_val = aug[col * w + col].abs();
        for row in (col + 1)..n {
            let v = aug[row * w + col].abs();
            if v > pivot_val {
                pivot_val = v;
                pivot_row = row;
            }
        }
        if pivot_val < 1e-15 {
            return None;
        }
        aug.swap_with_slice_helper(col, pivot_row, w);
        let inv = 1.0 / aug[col * w + col];
        for j in 0..w {
            aug[col * w + j] *= inv;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row * w + col];
            for j in 0..w {
                let v = aug[col * w + j];
                aug[row * w + j] -= factor * v;
            }
        }
    }
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = aug[i * w + n + j];
        }
    }
    Some(inv)
}
/// Helper: swap two rows of length `w` in a flat array.
trait SwapRows {
    fn swap_with_slice_helper(&mut self, r1: usize, r2: usize, w: usize);
}
impl SwapRows for Vec<f64> {
    fn swap_with_slice_helper(&mut self, r1: usize, r2: usize, w: usize) {
        if r1 == r2 {
            return;
        }
        for j in 0..w {
            self.swap(r1 * w + j, r2 * w + j);
        }
    }
}
/// Multiply rectangular matrix A (na×ka) by B (ka×nb).
///
/// Returns C (na×nb) = A * B.
pub(super) fn mat_mul_rect(a: &[f64], na: usize, ka: usize, b: &[f64], nb: usize) -> Vec<f64> {
    let mut c = vec![0.0; na * nb];
    for i in 0..na {
        for k in 0..ka {
            let aik = a[i * ka + k];
            for j in 0..nb {
                c[i * nb + j] += aik * b[k * nb + j];
            }
        }
    }
    c
}
/// Transpose a (r×c) matrix to (c×r).
pub(super) fn mat_transpose_rect(a: &[f64], r: usize, c: usize) -> Vec<f64> {
    let mut t = vec![0.0; c * r];
    for i in 0..r {
        for j in 0..c {
            t[j * r + i] = a[i * c + j];
        }
    }
    t
}
/// Estimate spectral radius via power iteration (for stability checks).
pub(super) fn spectral_radius(a: &[f64], n: usize) -> f64 {
    let mut v = vec![1.0; n];
    let norm = frob_norm(&v);
    for x in &mut v {
        *x /= norm;
    }
    let mut lambda = 0.0;
    for _ in 0..200 {
        let av = mat_vec(a, &v, n);
        lambda = frob_norm(&av);
        if lambda < 1e-15 {
            return 0.0;
        }
        for (x, av_x) in v.iter_mut().zip(av.iter()) {
            *x = av_x / lambda;
        }
    }
    lambda
}
/// Solve the discrete algebraic Riccati equation (DARE) iteratively.
///
/// Returns (K, P) where K = (R + Bᵀ P B)⁻¹ Bᵀ P A  and
/// P satisfies P = Q + Aᵀ P A - Aᵀ P B (R + Bᵀ P B)⁻¹ Bᵀ P A.
pub fn solve_dare(
    a: &[f64],
    b: &[f64],
    q: &[f64],
    r: &[f64],
    n: usize,
    m: usize,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, Vec<f64>) {
    let at = mat_transpose(a, n);
    let bt = mat_transpose_rect(b, n, m);
    let mut p = q.to_vec();
    for _ in 0..max_iter {
        let bt_p = mat_mul_rect(&bt, m, n, &p, n);
        let bt_p_b = mat_mul_rect(&bt_p, m, n, b, m);
        let s = mat_add(r, &bt_p_b, m);
        let s_inv = mat_inv(&s, m).expect("R + BᵀPB must be invertible");
        let bt_p_a = mat_mul_rect(&bt_p, m, n, a, n);
        let k = mat_mul_rect(&s_inv, m, m, &bt_p_a, n);
        let at_p = mat_mul_rect(&at, n, n, &p, n);
        let at_p_a = mat_mul_rect(&at_p, n, n, a, n);
        let at_p_b = mat_mul_rect(&at_p, n, n, b, m);
        let at_p_b_k = mat_mul_rect(&at_p_b, n, m, &k, n);
        let p_new = mat_sub(&mat_add(q, &at_p_a, n), &at_p_b_k, n);
        let diff = frob_norm(&mat_sub(&p_new, &p, n));
        p = p_new;
        if diff < tol {
            let bt_p2 = mat_mul_rect(&bt, m, n, &p, n);
            let bt_p_b2 = mat_mul_rect(&bt_p2, m, n, b, m);
            let s2 = mat_add(r, &bt_p_b2, m);
            let s_inv2 = mat_inv(&s2, m).expect("R + BᵀPB must be invertible");
            let bt_p_a2 = mat_mul_rect(&bt_p2, m, n, a, n);
            let k_final = mat_mul_rect(&s_inv2, m, m, &bt_p_a2, n);
            return (k_final, p);
        }
    }
    let bt_p = mat_mul_rect(&bt, m, n, &p, n);
    let bt_p_b = mat_mul_rect(&bt_p, m, n, b, m);
    let s = mat_add(r, &bt_p_b, m);
    let s_inv = mat_inv(&s, m).expect("R + BᵀPB must be invertible");
    let bt_p_a = mat_mul_rect(&bt_p, m, n, a, n);
    let k = mat_mul_rect(&s_inv, m, m, &bt_p_a, n);
    (k, p)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::optimal_control::DifferentialDynamicProgramming;
    use crate::optimal_control::FiniteHorizonLqr;
    use crate::optimal_control::HamiltonJacobiBellman;
    use crate::optimal_control::IlqrSolver;
    use crate::optimal_control::KalmanFilter;
    use crate::optimal_control::LqgController;
    use crate::optimal_control::LqrController;
    use crate::optimal_control::MinimumTimeControl;
    use crate::optimal_control::ModelPredictiveControl;
    use crate::optimal_control::MpcProblem;
    use crate::optimal_control::PontryaginMinimum;
    use crate::optimal_control::RobustControl;
    use crate::optimal_control::TrajectoryOptimization;
    fn double_integrator() -> (usize, usize, Vec<f64>, Vec<f64>) {
        let n = 2;
        let m = 1;
        let dt = 0.1_f64;
        let a = vec![1.0, dt, 0.0, 1.0];
        let b = vec![dt * dt / 2.0, dt];
        (n, m, a, b)
    }
    #[test]
    fn test_dare_p_positive_definite() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let (_, p) = solve_dare(&a, &b, &q, &r, n, m, 1000, 1e-10);
        for i in 0..n {
            assert!(
                p[i * n + i] > 0.0,
                "P[{i},{i}] = {} should be positive",
                p[i * n + i]
            );
        }
    }
    #[test]
    fn test_dare_riccati_residual() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let (k, p) = solve_dare(&a, &b, &q, &r, n, m, 2000, 1e-12);
        let _bt = mat_transpose_rect(&b, n, m);
        let at = mat_transpose(&a, n);
        let at_p = mat_mul_rect(&at, n, n, &p, n);
        let at_p_a = mat_mul_rect(&at_p, n, n, &a, n);
        let at_p_b = mat_mul_rect(&at_p, n, n, &b, m);
        let at_p_b_k = mat_mul_rect(&at_p_b, n, m, &k, n);
        let rhs = mat_sub(&mat_add(&q, &at_p_a, n), &at_p_b_k, n);
        let resid = frob_norm(&mat_sub(&p, &rhs, n));
        assert!(resid < 1e-6, "DARE residual = {resid:.2e}");
    }
    #[test]
    fn test_lqr_closed_loop_stable() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let ctrl = LqrController::new(n, m, a, b, q, r);
        assert!(ctrl.is_stable(), "LQR closed-loop should be stable");
    }
    #[test]
    fn test_lqr_control_dimension() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let ctrl = LqrController::new(n, m, a, b, q, r);
        let x = vec![1.0, 0.0];
        let u = ctrl.control(&x);
        assert_eq!(u.len(), m);
    }
    #[test]
    fn test_lqr_zero_state_zero_control() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let ctrl = LqrController::new(n, m, a, b, q, r);
        let x = vec![0.0; n];
        let u = ctrl.control(&x);
        for ui in &u {
            assert!(ui.abs() < 1e-10, "u at zero state should be zero, got {ui}");
        }
    }
    #[test]
    fn test_lqr_heavier_r_smaller_control() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r1 = mat_eye(m);
        let r100 = vec![100.0];
        let x = vec![1.0, 0.0];
        let ctrl1 = LqrController::new(n, m, a.clone(), b.clone(), q.clone(), r1);
        let ctrl100 = LqrController::new(n, m, a, b, q, r100);
        let u1 = ctrl1.control(&x);
        let u100 = ctrl100.control(&x);
        assert!(
            u1[0].abs() > u100[0].abs(),
            "Heavier R should reduce control: u1={:.4}, u100={:.4}",
            u1[0].abs(),
            u100[0].abs()
        );
    }
    #[test]
    fn test_finite_lqr_gains_count() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let horizon = 10;
        let fh = FiniteHorizonLqr::new(n, m, &a, &b, &q, &r, &q, horizon);
        assert_eq!(fh.gains.len(), horizon);
    }
    #[test]
    fn test_finite_lqr_zero_state() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let fh = FiniteHorizonLqr::new(n, m, &a, &b, &q, &r, &q, 10);
        let u = fh.control(0, &vec![0.0; n]);
        for ui in &u {
            assert!(ui.abs() < 1e-10);
        }
    }
    fn make_kalman() -> KalmanFilter {
        let n = 2;
        let p = 1;
        let m = 1;
        let dt = 0.1_f64;
        let a = vec![1.0, dt, 0.0, 1.0];
        let b = vec![dt * dt / 2.0, dt];
        let c = vec![1.0, 0.0];
        let q_noise = vec![1e-4, 0.0, 0.0, 1e-4];
        let r_noise = vec![0.01];
        let x0 = vec![0.0, 0.0];
        let p0 = mat_eye(n);
        KalmanFilter::new(n, p, m, a, b, c, q_noise, r_noise, x0, p0)
    }
    #[test]
    fn test_kalman_predict_increases_covariance() {
        let mut kf = make_kalman();
        let trace_before: f64 = (0..kf.n).map(|i| kf.cov[i * kf.n + i]).sum();
        kf.predict(&[0.0]);
        let trace_after: f64 = (0..kf.n).map(|i| kf.cov[i * kf.n + i]).sum();
        assert!(
            trace_after >= trace_before - 1e-10,
            "Covariance trace should not decrease after predict: before={trace_before}, after={trace_after}"
        );
    }
    #[test]
    fn test_kalman_update_decreases_covariance() {
        let mut kf = make_kalman();
        kf.predict(&[0.0]);
        let trace_before: f64 = (0..kf.n).map(|i| kf.cov[i * kf.n + i]).sum();
        kf.update(&[0.5]);
        let trace_after: f64 = (0..kf.n).map(|i| kf.cov[i * kf.n + i]).sum();
        assert!(
            trace_after <= trace_before + 1e-10,
            "Covariance trace should decrease after update: before={trace_before}, after={trace_after}"
        );
    }
    #[test]
    fn test_kalman_tracks_constant_position() {
        let mut kf = make_kalman();
        for _ in 0..50 {
            kf.predict(&[0.0]);
            kf.update(&[1.0]);
        }
        let pos_estimate = kf.x_hat[0];
        assert!(
            (pos_estimate - 1.0).abs() < 0.1,
            "Kalman should track position 1.0, got {pos_estimate:.4}"
        );
    }
    #[test]
    fn test_kalman_steady_state_gain_dims() {
        let kf = make_kalman();
        let kg = kf.steady_state_gain(200);
        assert_eq!(kg.len(), kf.n * kf.p);
    }
    #[test]
    fn test_kalman_steady_state_gain_positive() {
        let kf = make_kalman();
        let kg = kf.steady_state_gain(200);
        assert!(kg[0] > 0.0, "Steady-state gain K[0,0]={:.6}", kg[0]);
    }
    #[test]
    fn test_lqg_step_dimension() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let lqr = LqrController::new(n, m, a.clone(), b.clone(), q, r);
        let kf = make_kalman();
        let mut lqg = LqgController::new(lqr, kf);
        let u = lqg.step(&[0.5]);
        assert_eq!(u.len(), m);
    }
    #[test]
    fn test_lqg_does_not_diverge() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let lqr = LqrController::new(n, m, a.clone(), b.clone(), q, r);
        let kf = make_kalman();
        let mut lqg = LqgController::new(lqr, kf);
        for _ in 0..100 {
            let u = lqg.step(&[0.1]);
            assert!(u[0].abs() < 1e6, "LQG control diverged: {}", u[0]);
        }
    }
    fn make_mpc_problem() -> MpcProblem {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        MpcProblem {
            horizon: 5,
            n,
            m,
            a,
            b,
            q: q.clone(),
            qf: q,
            r,
            x_lb: vec![f64::NEG_INFINITY; n],
            x_ub: vec![f64::INFINITY; n],
            u_lb: vec![-2.0],
            u_ub: vec![2.0],
        }
    }
    #[test]
    fn test_mpc_solution_dimensions() {
        let prob = make_mpc_problem();
        let h = prob.horizon;
        let m = prob.m;
        let mpc = ModelPredictiveControl::new(prob);
        let sol = mpc.solve(&[1.0, 0.0]);
        assert_eq!(sol.u_seq.len(), h * m);
    }
    #[test]
    fn test_mpc_input_constraints() {
        let prob = make_mpc_problem();
        let h = prob.horizon;
        let m = prob.m;
        let mpc = ModelPredictiveControl::new(prob);
        let sol = mpc.solve(&[1.0, 0.5]);
        for k in 0..h {
            for j in 0..m {
                let u = sol.u_seq[k * m + j];
                assert!(u >= -2.0 - 1e-6, "u[{k}][{j}]={u} < -2");
                assert!(u <= 2.0 + 1e-6, "u[{k}][{j}]={u} > 2");
            }
        }
    }
    #[test]
    fn test_mpc_reduces_state() {
        let prob = make_mpc_problem();
        let mpc = ModelPredictiveControl::new(prob);
        let sol = mpc.solve(&[1.0, 0.0]);
        assert!(
            sol.u_seq[0].abs() > 1e-6,
            "MPC should produce non-zero control"
        );
    }
    #[test]
    fn test_mpc_feasible() {
        let prob = make_mpc_problem();
        let mpc = ModelPredictiveControl::new(prob);
        let sol = mpc.solve(&[0.1, 0.0]);
        assert!(
            sol.feasible,
            "MPC should be feasible for small initial state"
        );
    }
    #[test]
    fn test_pontryagin_costate_changes() {
        let n = 2;
        let m = 1;
        let a = vec![0.0_f64, 1.0, -1.0, -0.5];
        let q = mat_eye(n);
        let x = vec![1.0, 0.0];
        let lam0 = vec![0.5, 0.5];
        let mut pmp = PontryaginMinimum::new(n, m, lam0.clone());
        pmp.costate_step(&a, &q, &x, 0.01);
        let changed = pmp
            .costate
            .iter()
            .zip(lam0.iter())
            .any(|(a, b)| (a - b).abs() > 1e-12);
        assert!(changed, "Co-state should change after backward step");
    }
    #[test]
    fn test_pontryagin_hamiltonian_finite() {
        let (n, m, a, b) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let x = vec![1.0, 0.5];
        let u = vec![-0.2];
        let lam = vec![0.3, 0.1];
        let pmp = PontryaginMinimum::new(n, m, lam);
        let h = pmp.eval_hamiltonian(&a, &b, &q, &r, &x, &u);
        assert!(h.is_finite(), "Hamiltonian should be finite, got {h}");
    }
    #[test]
    fn test_pontryagin_optimal_control_dim() {
        let (n, m, _, b) = double_integrator();
        let r = mat_eye(m);
        let lam = vec![0.3, 0.1];
        let pmp = PontryaginMinimum::new(n, m, lam);
        let u = pmp.optimal_control(&b, &r);
        assert_eq!(u.len(), m);
    }
    #[test]
    fn test_hjb_initial_value_zero() {
        let hjb = HamiltonJacobiBellman::new(10, -1.0, 1.0, 0.99);
        for v in &hjb.value {
            assert!(v.abs() < 1e-10);
        }
    }
    #[test]
    fn test_hjb_bellman_sweep_changes_values() {
        let mut hjb = HamiltonJacobiBellman::new(5, -1.0, 1.0, 0.99);
        let dynamics = |x1: f64, x2: f64, u: f64| -> (f64, f64) { (x1 + 0.1 * x2, x2 + 0.1 * u) };
        let cost = |x1: f64, x2: f64, u: f64| -> f64 { x1 * x1 + x2 * x2 + u * u };
        let u_grid: Vec<f64> = (-5..=5).map(|i| i as f64).collect();
        hjb.bellman_sweep(dynamics, cost, &u_grid);
        let any_nonzero = hjb.value.iter().any(|v| v.abs() > 1e-10);
        assert!(
            any_nonzero,
            "Value function should be updated after Bellman sweep"
        );
    }
    #[test]
    fn test_hjb_set_get_value() {
        let mut hjb = HamiltonJacobiBellman::new(10, 0.0, 1.0, 1.0);
        hjb.set_value(0.5, 0.5, 42.0);
        let v = hjb.get_value(0.5, 0.5);
        assert!((v - 42.0).abs() < 1.0, "Set/get should round-trip: got {v}");
    }
    #[test]
    fn test_ilqr_linear_system_cost() {
        let (n, m, a_mat, b_mat) = double_integrator();
        let q = mat_eye(n);
        let r = mat_eye(m);
        let dynamics = |x: &[f64], u: &[f64]| -> Vec<f64> {
            let ax = mat_vec(&a_mat, x, n);
            let bu = mat_mul_rect(&b_mat, n, m, u, 1);
            ax.iter().zip(bu.iter()).map(|(a, b)| a + b).collect()
        };
        let running_cost = |x: &[f64], u: &[f64]| -> f64 {
            dot(x, &mat_vec(&q, x, n)) + dot(u, &mat_vec(&r, u, m))
        };
        let terminal_cost = |x: &[f64]| -> f64 { dot(x, &mat_vec(&q, x, n)) };
        let solver = IlqrSolver::new(n, m, 20, 50, 1e-6, 1e-4);
        let x0 = vec![1.0, 0.0];
        let u_init = vec![0.0; 20 * m];
        let result = solver.solve(&x0, &u_init, dynamics, running_cost, terminal_cost);
        let init_cost: f64 =
            (0..20).map(|_| running_cost(&x0, &[0.0])).sum::<f64>() + terminal_cost(&x0);
        assert!(
            result.cost < init_cost,
            "iLQR should reduce cost: init={init_cost:.4}, result={:.4}",
            result.cost
        );
    }
    #[test]
    fn test_ilqr_trajectory_shape() {
        let (n, m, a_mat, b_mat) = double_integrator();
        let dynamics = |x: &[f64], u: &[f64]| -> Vec<f64> {
            let ax = mat_vec(&a_mat, x, n);
            let bu = mat_mul_rect(&b_mat, n, m, u, 1);
            ax.iter().zip(bu.iter()).map(|(a, b)| a + b).collect()
        };
        let t = 10;
        let solver = IlqrSolver::new(n, m, t, 5, 1e-4, 1e-3);
        let x0 = vec![0.5, 0.0];
        let u_init = vec![0.0; t * m];
        let result = solver.solve(
            &x0,
            &u_init,
            dynamics,
            |x, u| dot(x, x) + dot(u, u),
            |x| dot(x, x),
        );
        assert_eq!(result.x_traj.len(), (t + 1) * n);
        assert_eq!(result.u_seq.len(), t * m);
    }
    #[test]
    fn test_ddp_finite_cost() {
        let (n, m, a_mat, b_mat) = double_integrator();
        let dynamics = |x: &[f64], u: &[f64]| -> Vec<f64> {
            let ax = mat_vec(&a_mat, x, n);
            let bu = mat_mul_rect(&b_mat, n, m, u, 1);
            ax.iter().zip(bu.iter()).map(|(a, b)| a + b).collect()
        };
        let ddp = DifferentialDynamicProgramming::new(n, m, 10, 20, 1e-4);
        let result = ddp.solve(
            &[1.0, 0.0],
            &vec![0.0; 10 * m],
            dynamics,
            |x, u| dot(x, x) + dot(u, u),
            |x| dot(x, x),
        );
        assert!(
            result.cost.is_finite(),
            "DDP cost should be finite: {}",
            result.cost
        );
    }
    #[test]
    fn test_direct_collocation_defects_exact() {
        let n = 1;
        let m = 1;
        let dt = 0.1_f64;
        let traj_opt = TrajectoryOptimization::new(n, m, 3, dt);
        let dynamics = |_x: &[f64], _u: &[f64]| vec![0.0];
        let traj = vec![1.0, 1.0, 1.0];
        let u_seq = vec![0.0, 0.0];
        let defects = traj_opt.defects(&traj, &u_seq, dynamics);
        for (i, d) in defects.iter().enumerate() {
            assert!(d.abs() < 1e-10, "Defect[{i}] = {d} should be zero");
        }
    }
    #[test]
    fn test_direct_collocation_nonzero_defects() {
        let n = 1;
        let m = 1;
        let dt = 0.1_f64;
        let traj_opt = TrajectoryOptimization::new(n, m, 3, dt);
        let dynamics = |x: &[f64], _u: &[f64]| vec![1.0 * x[0]];
        let traj = vec![1.0, 5.0, 3.0];
        let u_seq = vec![0.0, 0.0];
        let defects = traj_opt.defects(&traj, &u_seq, dynamics);
        let any_nonzero = defects.iter().any(|d| d.abs() > 1e-6);
        assert!(
            any_nonzero,
            "Inconsistent trajectory should have nonzero defects"
        );
    }
    #[test]
    fn test_bangbang_magnitude() {
        let ctrl = MinimumTimeControl::new(1.0, 0.0);
        let u = ctrl.control(1.0, 0.0);
        assert!(
            (u.abs() - 1.0).abs() < 1e-10,
            "Bang-bang magnitude should equal u_max, got {u}"
        );
    }
    #[test]
    fn test_bangbang_sign() {
        let ctrl = MinimumTimeControl::new(1.0, 0.0);
        let u = ctrl.control(0.0, 0.5);
        assert!(
            u < 0.0,
            "Should decelerate when moving toward target, got {u}"
        );
    }
    #[test]
    fn test_bangbang_min_time_at_rest() {
        let ctrl = MinimumTimeControl::new(1.0, 5.0);
        let t = ctrl.min_time(5.0, 0.0);
        assert!(
            t < 1e-8,
            "Min time at target with zero vel should be ~0, got {t}"
        );
    }
    #[test]
    fn test_bangbang_min_time_positive() {
        let ctrl = MinimumTimeControl::new(1.0, 0.0);
        let t = ctrl.min_time(5.0, 0.0);
        assert!(
            t > 0.0,
            "Min time away from target should be positive, got {t}"
        );
    }
    #[test]
    fn test_robust_h_inf_finite() {
        let n = 2;
        let m = 1;
        let p = 1;
        let a = vec![-1.0, 0.0, 0.0, -2.0];
        let b = vec![1.0, 0.0];
        let c = vec![1.0, 0.0];
        let d = vec![0.0];
        let rc = RobustControl::new(n, m, p, a, b, c, d);
        let bound = rc.h_inf_upper_bound();
        assert!(bound.is_finite(), "H-inf bound should be finite: {bound}");
    }
    #[test]
    fn test_robust_small_gain_stable() {
        let n = 2;
        let m = 1;
        let p = 1;
        let a = vec![-10.0, 0.0, 0.0, -10.0];
        let b = vec![1.0, 0.0];
        let c = vec![1.0, 0.0];
        let d = vec![0.0];
        let rc = RobustControl::new(n, m, p, a, b, c, d);
        assert!(
            rc.small_gain_stable(0.01),
            "Should be stable with small gamma"
        );
    }
    #[test]
    fn test_mat_eye() {
        let i = mat_eye(3);
        for r in 0..3 {
            for c in 0..3 {
                let expected = if r == c { 1.0 } else { 0.0 };
                assert!((i[r * 3 + c] - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_mat_mul_identity() {
        let i = mat_eye(2);
        let m = vec![1.0, 2.0, 3.0, 4.0];
        let r = mat_mul(&i, &m, 2);
        for (a, b) in r.iter().zip(m.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
    #[test]
    fn test_mat_transpose_self_inverse() {
        let m = vec![1.0, 2.0, 3.0, 4.0];
        let tt = mat_transpose(&mat_transpose(&m, 2), 2);
        for (a, b) in tt.iter().zip(m.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
    #[test]
    fn test_mat_inv_2x2() {
        let m = vec![2.0, 1.0, 5.0, 3.0];
        let inv = mat_inv(&m, 2).expect("Should be invertible");
        let prod = mat_mul(&m, &inv, 2);
        let i = mat_eye(2);
        for (a, b) in prod.iter().zip(i.iter()) {
            assert!(
                (a - b).abs() < 1e-8,
                "m*inv should be identity, got {a} vs {b}"
            );
        }
    }
    #[test]
    fn test_frob_norm_identity() {
        let n = 3;
        let i = mat_eye(n);
        let norm = frob_norm(&i);
        assert!(
            (norm - (n as f64).sqrt()).abs() < 1e-10,
            "frob_norm(I_3) = {norm}"
        );
    }
    #[test]
    fn test_spectral_radius_zero() {
        let z = vec![0.0; 4];
        let r = spectral_radius(&z, 2);
        assert!(
            r < 1e-10,
            "spectral radius of zero matrix should be 0, got {r}"
        );
    }
}
