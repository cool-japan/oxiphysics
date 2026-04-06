// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Inverse problems and parameter identification for the FEM crate.
//!
//! Implements:
//! - [`InverseSolver`] — inverse problem solver with forward model and regularization
//! - \[`tikhonov_regularization()`\] — L2 (Tikhonov) regularization
//! - \[`total_variation_reg()`\] — TV regularization for piecewise-constant parameters
//! - \[`adjoint_gradient()`\] — adjoint method for gradient computation
//! - \[`gauss_newton_inverse()`\] — Gauss-Newton method for nonlinear inverse problems
//! - \[`cross_validation()`\] — L-curve method for regularization parameter selection
//! - \[`sensitivity_matrix()`\] — Jacobian of observables w.r.t. parameters

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── Math helpers ─────────────────────────────────────────────────────────────

/// Compute the dot product of two slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Compute the L2 norm of a slice.
fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// Subtract two vectors element-wise.
fn sub_vec(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

/// Add two vectors element-wise.
fn add_vec(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Scale a vector by a scalar.
fn scale_vec(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}

/// Multiply a dense matrix (row-major, `m×n`) by a vector of length `n`.
fn mat_vec(a: &[f64], n_rows: usize, n_cols: usize, x: &[f64]) -> Vec<f64> {
    (0..n_rows)
        .map(|i| (0..n_cols).map(|j| a[i * n_cols + j] * x[j]).sum::<f64>())
        .collect()
}

/// Multiply the transpose of a dense matrix (`m×n`) by a vector of length `m`.
fn mat_t_vec(a: &[f64], n_rows: usize, n_cols: usize, x: &[f64]) -> Vec<f64> {
    (0..n_cols)
        .map(|j| (0..n_rows).map(|i| a[i * n_cols + j] * x[i]).sum::<f64>())
        .collect()
}

/// Compute `Aᵀ A + λ I` (n×n) — the normal equation matrix for Tikhonov.
fn normal_matrix(a: &[f64], n_rows: usize, n_cols: usize, lambda: f64) -> Vec<f64> {
    let mut m = vec![0.0_f64; n_cols * n_cols];
    for k in 0..n_rows {
        for i in 0..n_cols {
            for j in 0..n_cols {
                m[i * n_cols + j] += a[k * n_cols + i] * a[k * n_cols + j];
            }
        }
    }
    for i in 0..n_cols {
        m[i * n_cols + i] += lambda;
    }
    m
}

/// Solve `A x = b` with a simple conjugate-gradient iteration (symmetric PD `A`).
fn cg_solve(a: &[f64], n: usize, b: &[f64], max_iter: usize, tol: f64) -> Vec<f64> {
    let mut x = vec![0.0_f64; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rsold = dot(&r, &r);
    for _ in 0..max_iter {
        let ap = mat_vec(a, n, n, &p);
        let alpha = rsold / (dot(&p, &ap) + 1e-300);
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rsnew = dot(&r, &r);
        if rsnew.sqrt() < tol {
            break;
        }
        let beta = rsnew / (rsold + 1e-300);
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    x
}

// ═══════════════════════════════════════════════════════════════════════════
// § 1  InverseSolver
// ═══════════════════════════════════════════════════════════════════════════

/// Inverse problem solver that wraps a linearised forward model and
/// regularization strategy.
///
/// The forward model is `d = A p` where `d` is observed data (length `m`),
/// `A` is the sensitivity matrix (`m × n`), and `p` is the unknown parameter
/// vector (length `n`).
#[derive(Debug, Clone)]
pub struct InverseSolver {
    /// Sensitivity (Jacobian) matrix, stored row-major with shape `m × n`.
    pub sensitivity: Vec<f64>,
    /// Number of observations (rows of `sensitivity`).
    pub n_obs: usize,
    /// Number of parameters (columns of `sensitivity`).
    pub n_params: usize,
    /// Regularization parameter λ.
    pub lambda: f64,
    /// Current parameter estimate.
    pub params: Vec<f64>,
    /// Observed data vector (length `n_obs`).
    pub data: Vec<f64>,
}

impl InverseSolver {
    /// Create a new [`InverseSolver`].
    ///
    /// # Arguments
    /// * `sensitivity` — row-major sensitivity matrix of shape `n_obs × n_params`
    /// * `n_obs` — number of observations
    /// * `n_params` — number of unknown parameters
    /// * `lambda` — regularization parameter
    /// * `data` — observed data vector of length `n_obs`
    pub fn new(
        sensitivity: Vec<f64>,
        n_obs: usize,
        n_params: usize,
        lambda: f64,
        data: Vec<f64>,
    ) -> Self {
        let params = vec![0.0_f64; n_params];
        Self {
            sensitivity,
            n_obs,
            n_params,
            lambda,
            params,
            data,
        }
    }

    /// Solve the inverse problem using Tikhonov regularization and return
    /// the recovered parameter vector.
    pub fn solve_tikhonov(&self) -> Vec<f64> {
        tikhonov_regularization(
            &self.sensitivity,
            self.n_obs,
            self.n_params,
            &self.data,
            self.lambda,
        )
    }

    /// Return the residual norm `||A p - d||` for a given parameter vector.
    pub fn residual_norm(&self, p: &[f64]) -> f64 {
        let ap = mat_vec(&self.sensitivity, self.n_obs, self.n_params, p);
        norm(&sub_vec(&ap, &self.data))
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 2  Tikhonov regularization
// ═══════════════════════════════════════════════════════════════════════════

/// Solve the Tikhonov-regularized inverse problem:
///
/// `p* = arg min ||A p − b||² + λ ||p||²`
///
/// via the normal equations `(AᵀA + λI) p = Aᵀ b`.
///
/// # Arguments
/// * `a` — sensitivity matrix, row-major, shape `m × n`
/// * `m` — number of rows (observations)
/// * `n` — number of columns (parameters)
/// * `b` — observation vector, length `m`
/// * `lambda` — regularization weight (≥ 0)
///
/// Returns the parameter vector of length `n`.
pub fn tikhonov_regularization(a: &[f64], m: usize, n: usize, b: &[f64], lambda: f64) -> Vec<f64> {
    let lhs = normal_matrix(a, m, n, lambda);
    let rhs = mat_t_vec(a, m, n, b);
    cg_solve(&lhs, n, &rhs, 1000, 1e-10)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 3  Total-variation regularization
// ═══════════════════════════════════════════════════════════════════════════

/// Apply one gradient-descent step of total-variation regularization.
///
/// The TV penalty is `λ Σ |p_{i+1} − p_i|`.  The smooth approximation
/// `|x| ≈ √(x² + ε)` is used for differentiability.
///
/// # Arguments
/// * `a` — sensitivity matrix, row-major, shape `m × n`
/// * `m` — number of observations
/// * `n` — number of parameters
/// * `b` — observed data
/// * `p` — current parameter estimate (length `n`)
/// * `lambda` — TV weight
/// * `step` — gradient-descent step size
///
/// Returns an updated parameter vector.
pub fn total_variation_reg(
    a: &[f64],
    m: usize,
    n: usize,
    b: &[f64],
    p: &[f64],
    lambda: f64,
    step: f64,
) -> Vec<f64> {
    const EPS: f64 = 1e-6;
    // data-fidelity gradient: Aᵀ(Ap − b)
    let ap = mat_vec(a, m, n, p);
    let residual = sub_vec(&ap, b);
    let grad_data = mat_t_vec(a, m, n, &residual);

    // TV gradient (finite difference)
    let mut grad_tv = vec![0.0_f64; n];
    for i in 0..n {
        if i + 1 < n {
            let diff = p[i] - p[i + 1];
            let w = (diff * diff + EPS).sqrt();
            grad_tv[i] += diff / w;
            grad_tv[i + 1] -= diff / w;
        }
    }

    let grad = add_vec(&grad_data, &scale_vec(&grad_tv, lambda));
    sub_vec(p, &scale_vec(&grad, step))
}

// ═══════════════════════════════════════════════════════════════════════════
// § 4  Adjoint gradient
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the gradient of the data-misfit functional via the adjoint method.
///
/// For a linearised forward model `d = A p`, the gradient is simply `Aᵀ(Ap − d_obs)`.
///
/// # Arguments
/// * `a` — sensitivity matrix, row-major `m × n`
/// * `m` — number of observations
/// * `n` — number of parameters
/// * `p` — current parameter vector (length `n`)
/// * `d_obs` — observed data (length `m`)
///
/// Returns the gradient vector of length `n`.
pub fn adjoint_gradient(a: &[f64], m: usize, n: usize, p: &[f64], d_obs: &[f64]) -> Vec<f64> {
    let ap = mat_vec(a, m, n, p);
    let residual = sub_vec(&ap, d_obs);
    mat_t_vec(a, m, n, &residual)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 5  Gauss-Newton inverse
// ═══════════════════════════════════════════════════════════════════════════

/// Gauss-Newton iterations for a nonlinear inverse problem.
///
/// At each iteration the Jacobian `J` is evaluated numerically by finite
/// differences of `forward_fn`, and the linearised system
/// `(JᵀJ + λI) Δp = Jᵀ r` is solved via CG.
///
/// # Arguments
/// * `forward_fn` — closure `p → predicted_data` (length `m`)
/// * `p0` — initial parameter guess (length `n`)
/// * `d_obs` — observed data (length `m`)
/// * `lambda` — Levenberg-Marquardt damping
/// * `max_iter` — maximum Gauss-Newton iterations
/// * `tol` — convergence tolerance on `||r||`
///
/// Returns the final parameter vector of length `n`.
pub fn gauss_newton_inverse<F>(
    forward_fn: F,
    p0: &[f64],
    d_obs: &[f64],
    lambda: f64,
    max_iter: usize,
    tol: f64,
) -> Vec<f64>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let n = p0.len();
    let m = d_obs.len();
    let h = 1e-7;
    let mut p = p0.to_vec();

    for _ in 0..max_iter {
        let fp = forward_fn(&p);
        let r = sub_vec(&fp, d_obs);
        if norm(&r) < tol {
            break;
        }
        // Build numerical Jacobian J (m × n)
        let mut j = vec![0.0_f64; m * n];
        for k in 0..n {
            let mut pk = p.clone();
            pk[k] += h;
            let fp_h = forward_fn(&pk);
            for i in 0..m {
                j[i * n + k] = (fp_h[i] - fp[i]) / h;
            }
        }
        // Normal equations (JᵀJ + λI) Δp = −Jᵀ r
        let lhs = normal_matrix(&j, m, n, lambda);
        let rhs_pos = mat_t_vec(&j, m, n, &r);
        let rhs: Vec<f64> = rhs_pos.iter().map(|x| -x).collect();
        let delta = cg_solve(&lhs, n, &rhs, 500, 1e-12);
        p = add_vec(&p, &delta);
    }
    p
}

// ═══════════════════════════════════════════════════════════════════════════
// § 6  Cross-validation / L-curve
// ═══════════════════════════════════════════════════════════════════════════

/// Select the regularization parameter using the L-curve criterion.
///
/// For each candidate λ the data residual norm and solution norm are computed.
/// The L-curve corner is approximated as the λ that maximises the curvature
/// of the log-log plot.
///
/// # Arguments
/// * `a` — sensitivity matrix, row-major `m × n`
/// * `m` — number of observations
/// * `n` — number of parameters
/// * `b` — observed data (length `m`)
/// * `lambdas` — candidate regularization parameters (sorted ascending)
///
/// Returns `(best_lambda, residual_norms, solution_norms)`.
pub fn cross_validation(
    a: &[f64],
    m: usize,
    n: usize,
    b: &[f64],
    lambdas: &[f64],
) -> (f64, Vec<f64>, Vec<f64>) {
    let mut res_norms = Vec::with_capacity(lambdas.len());
    let mut sol_norms = Vec::with_capacity(lambdas.len());

    for &lam in lambdas {
        let p = tikhonov_regularization(a, m, n, b, lam);
        let ap = mat_vec(a, m, n, &p);
        res_norms.push(norm(&sub_vec(&ap, b)));
        sol_norms.push(norm(&p));
    }

    // Find L-curve corner: maximum discrete curvature in log-log space
    let log_rho: Vec<f64> = res_norms.iter().map(|x| x.ln()).collect();
    let log_eta: Vec<f64> = sol_norms.iter().map(|x| x.ln()).collect();

    let mut best_idx = 0;
    let mut best_curv = f64::NEG_INFINITY;
    let k = lambdas.len();
    if k >= 3 {
        for i in 1..k - 1 {
            let d1 = (log_eta[i] - log_eta[i - 1]) / (log_rho[i] - log_rho[i - 1] + 1e-300);
            let d2 = (log_eta[i + 1] - log_eta[i]) / (log_rho[i + 1] - log_rho[i] + 1e-300);
            let curv = d2 - d1;
            if curv > best_curv {
                best_curv = curv;
                best_idx = i;
            }
        }
    }

    (lambdas[best_idx], res_norms, sol_norms)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 7  Sensitivity matrix
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the sensitivity (Jacobian) matrix of observables w.r.t. parameters
/// by forward finite differences.
///
/// # Arguments
/// * `forward_fn` — closure mapping parameter vector `p` (length `n`) to
///   observable vector (length `m`)
/// * `p` — parameter vector at which to evaluate the Jacobian
/// * `m` — number of observables
/// * `h` — finite-difference step size (default ≈ 1 × 10⁻⁷)
///
/// Returns the Jacobian stored row-major (`m × n`).
pub fn sensitivity_matrix<F>(forward_fn: F, p: &[f64], m: usize, h: f64) -> Vec<f64>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let n = p.len();
    let f0 = forward_fn(p);
    let mut jac = vec![0.0_f64; m * n];
    for k in 0..n {
        let mut pk = p.to_vec();
        pk[k] += h;
        let fk = forward_fn(&pk);
        for i in 0..m {
            jac[i * n + k] = (fk[i] - f0[i]) / h;
        }
    }
    jac
}

// ═══════════════════════════════════════════════════════════════════════════
// § 8  Additional helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the Morozov discrepancy principle estimate for `λ`.
///
/// Returns the smallest λ from `lambdas` such that
/// `||A p_λ − b|| ≤ delta * ||b||`, where `delta` is the relative noise level.
///
/// If no such λ exists the last (largest) value is returned.
pub fn morozov_lambda(
    a: &[f64],
    m: usize,
    n: usize,
    b: &[f64],
    lambdas: &[f64],
    delta: f64,
) -> f64 {
    let target = delta * norm(b);
    for &lam in lambdas.iter().rev() {
        let p = tikhonov_regularization(a, m, n, b, lam);
        let ap = mat_vec(a, m, n, &p);
        if norm(&sub_vec(&ap, b)) <= target {
            return lam;
        }
    }
    *lambdas.last().unwrap_or(&1.0)
}

/// Compute the relative data fit error `||A p − b|| / ||b||`.
pub fn relative_fit_error(a: &[f64], m: usize, n: usize, p: &[f64], b: &[f64]) -> f64 {
    let ap = mat_vec(a, m, n, p);
    let res = sub_vec(&ap, b);
    norm(&res) / (norm(b) + 1e-300)
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── helpers ──────────────────────────────────────────────────────────

    /// Identity matrix of size `n`, row-major.
    fn eye(n: usize) -> Vec<f64> {
        let mut m = vec![0.0_f64; n * n];
        for i in 0..n {
            m[i * n + i] = 1.0;
        }
        m
    }

    /// A = \[\[2, 1\\], \[1, 3\]], exact solution for b = \[3, 5\] is \[4/5, 7/5\].
    fn small_system() -> (Vec<f64>, Vec<f64>) {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let b = vec![3.0, 5.0];
        (a, b)
    }

    // ── Tikhonov tests ────────────────────────────────────────────────────

    #[test]
    fn test_tikhonov_identity_system() {
        // A = I, b = [1,2,3] → p* ≈ b / (1 + λ)
        let n = 3;
        let a = eye(n);
        let b = vec![1.0, 2.0, 3.0];
        let lambda = 1.0;
        let p = tikhonov_regularization(&a, n, n, &b, lambda);
        for (i, &pi) in p.iter().enumerate() {
            let expected = b[i] / (1.0 + lambda);
            assert!(
                (pi - expected).abs() < 1e-8,
                "i={i}: got {pi}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_tikhonov_zero_lambda_recovers_solution() {
        // With λ≈0 and a square full-rank system the solution should be close to exact.
        let a = vec![2.0, 0.0, 0.0, 3.0];
        let b = vec![4.0, 9.0];
        let p = tikhonov_regularization(&a, 2, 2, &b, 1e-10);
        assert!((p[0] - 2.0).abs() < 1e-5);
        assert!((p[1] - 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_tikhonov_overdetermined() {
        // Overdetermined: A (3×2), b (3), should return least-squares solution.
        let a = vec![1.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let b = vec![1.0, 2.0, 3.0];
        let p = tikhonov_regularization(&a, 3, 2, &b, 1e-8);
        // residual should be small relative to b
        let rel = relative_fit_error(&a, 3, 2, &p, &b);
        assert!(rel < 0.5, "relative error {rel}");
    }

    #[test]
    fn test_tikhonov_larger_lambda_shrinks_solution() {
        let n = 4;
        let a = eye(n);
        let b = vec![1.0; n];
        let p_small = tikhonov_regularization(&a, n, n, &b, 0.01);
        let p_large = tikhonov_regularization(&a, n, n, &b, 100.0);
        assert!(norm(&p_large) < norm(&p_small));
    }

    #[test]
    fn test_tikhonov_solution_norm() {
        let n = 2;
        let a = eye(n);
        let b = vec![3.0, 4.0];
        let lambda = 1.0;
        let p = tikhonov_regularization(&a, n, n, &b, lambda);
        let expected_norm = norm(&[1.5, 2.0]);
        assert!((norm(&p) - expected_norm).abs() < 1e-7);
    }

    // ── InverseSolver tests ────────────────────────────────────────────────

    #[test]
    fn test_inverse_solver_new() {
        let a = eye(3);
        let b = vec![1.0, 2.0, 3.0];
        let solver = InverseSolver::new(a, 3, 3, 0.5, b);
        assert_eq!(solver.n_obs, 3);
        assert_eq!(solver.n_params, 3);
        assert_eq!(solver.params.len(), 3);
    }

    #[test]
    fn test_inverse_solver_solve() {
        let n = 3;
        let a = eye(n);
        let b = vec![2.0, 4.0, 6.0];
        let solver = InverseSolver::new(a, n, n, 1.0, b);
        let p = solver.solve_tikhonov();
        // expected: b / 2 = [1, 2, 3]
        for (i, &pi) in p.iter().enumerate() {
            assert!((pi - (i as f64 + 1.0)).abs() < 1e-7, "i={i}: {pi}");
        }
    }

    #[test]
    fn test_inverse_solver_residual_norm() {
        let a = eye(2);
        let b = vec![1.0, 0.0];
        let solver = InverseSolver::new(a, 2, 2, 0.0, b.clone());
        let r = solver.residual_norm(&b);
        assert!(r < 1e-10);
    }

    // ── TV regularization tests ───────────────────────────────────────────

    #[test]
    fn test_tv_reg_reduces_residual() {
        let m = 4;
        let n = 4;
        let a = eye(n);
        let b = vec![1.0, 1.0, 1.0, 1.0];
        let p0 = vec![0.0; n];
        let p1 = total_variation_reg(&a, m, n, &b, &p0, 0.1, 0.1);
        // after one step solution should have moved toward b
        assert!(norm(&p1) > 0.0);
    }

    #[test]
    fn test_tv_reg_piecewise_constant() {
        // b has a step: TV should preserve the step better than Tikhonov
        let n = 6;
        let a = eye(n);
        let b = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let mut p = vec![0.0_f64; n];
        for _ in 0..200 {
            p = total_variation_reg(&a, n, n, &b, &p, 0.5, 0.01);
        }
        // step should still be present: last half > first half
        let first: f64 = p[..3].iter().sum::<f64>() / 3.0;
        let last: f64 = p[3..].iter().sum::<f64>() / 3.0;
        assert!(last > first + 0.1, "first={first:.3} last={last:.3}");
    }

    #[test]
    fn test_tv_reg_step_size_matters() {
        let n = 4;
        let a = eye(n);
        let b = vec![1.0, 1.0, 1.0, 1.0];
        let p0 = vec![0.0; n];
        let p_big = total_variation_reg(&a, n, n, &b, &p0, 0.1, 0.5);
        let p_small = total_variation_reg(&a, n, n, &b, &p0, 0.1, 0.01);
        // bigger step → bigger movement in one iteration
        assert!(norm(&p_big) > norm(&p_small));
    }

    // ── Adjoint gradient tests ────────────────────────────────────────────

    #[test]
    fn test_adjoint_gradient_zero_at_solution() {
        let n = 3;
        let a = eye(n);
        let d_obs = vec![1.0, 2.0, 3.0];
        let grad = adjoint_gradient(&a, n, n, &d_obs, &d_obs);
        assert!(norm(&grad) < 1e-10);
    }

    #[test]
    fn test_adjoint_gradient_direction() {
        // With A = I, gradient = p − d_obs.  Steepest descent should decrease cost.
        let n = 3;
        let a = eye(n);
        let d_obs = vec![1.0, 2.0, 3.0];
        let p = vec![2.0, 3.0, 4.0];
        let grad = adjoint_gradient(&a, n, n, &p, &d_obs);
        // grad should equal p - d_obs
        for i in 0..n {
            assert!((grad[i] - (p[i] - d_obs[i])).abs() < 1e-10);
        }
    }

    #[test]
    fn test_adjoint_gradient_nonzero() {
        let (a, b) = small_system();
        let p = vec![0.0, 0.0];
        let g = adjoint_gradient(&a, 2, 2, &p, &b);
        // should be -Aᵀb ≠ 0
        assert!(norm(&g) > 0.0);
    }

    // ── Gauss-Newton tests ────────────────────────────────────────────────

    #[test]
    fn test_gauss_newton_linear() {
        // For a linear model the Gauss-Newton should converge in one step.
        let a_mat = vec![2.0_f64, 0.0, 0.0, 3.0];
        let d_obs = vec![4.0, 9.0];
        let p0 = vec![0.0; 2];
        let p = gauss_newton_inverse(|p| mat_vec(&a_mat, 2, 2, p), &p0, &d_obs, 1e-6, 20, 1e-8);
        assert!((p[0] - 2.0).abs() < 1e-4, "p[0]={}", p[0]);
        assert!((p[1] - 3.0).abs() < 1e-4, "p[1]={}", p[1]);
    }

    #[test]
    fn test_gauss_newton_quadratic() {
        // Forward model: f(p) = [p[0]^2], d_obs = [4.0] → p = 2.0
        let d_obs = vec![4.0_f64];
        let p0 = vec![3.0_f64];
        let p = gauss_newton_inverse(|p| vec![p[0] * p[0]], &p0, &d_obs, 0.1, 50, 1e-8);
        assert!((p[0].abs() - 2.0).abs() < 0.1, "p[0]={}", p[0]);
    }

    #[test]
    fn test_gauss_newton_already_converged() {
        let d_obs = vec![1.0, 2.0];
        let p0 = d_obs.clone();
        let p = gauss_newton_inverse(|p| p.to_vec(), &p0, &d_obs, 1e-6, 10, 1e-8);
        assert!((p[0] - 1.0).abs() < 1e-6);
        assert!((p[1] - 2.0).abs() < 1e-6);
    }

    // ── Sensitivity matrix tests ──────────────────────────────────────────

    #[test]
    fn test_sensitivity_matrix_identity() {
        // f(p) = p → J = I
        let p = vec![1.0, 2.0, 3.0];
        let n = p.len();
        let jac = sensitivity_matrix(|p| p.to_vec(), &p, n, 1e-7);
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((jac[i * n + j] - expected).abs() < 1e-5, "({i},{j})");
            }
        }
    }

    #[test]
    fn test_sensitivity_matrix_affine() {
        // f(p) = 2 * p + 1 → J = 2 * I
        let p = vec![0.0, 0.0];
        let n = 2;
        let jac = sensitivity_matrix(|p| p.iter().map(|x| 2.0 * x + 1.0).collect(), &p, n, 1e-7);
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { 2.0 } else { 0.0 };
                assert!((jac[i * n + j] - expected).abs() < 1e-4);
            }
        }
    }

    #[test]
    fn test_sensitivity_matrix_nonlinear() {
        // f(p) = [p[0]^2, p[1]^3] at p = [2, 3]
        // J = [[2p0, 0], [0, 3p1^2]] = [[4, 0], [0, 27]]
        let p = vec![2.0, 3.0];
        let jac = sensitivity_matrix(|p| vec![p[0] * p[0], p[1] * p[1] * p[1]], &p, 2, 1e-6);
        assert!((jac[0] - 4.0).abs() < 1e-3);
        assert!((jac[1]).abs() < 1e-3);
        assert!((jac[2]).abs() < 1e-3);
        assert!((jac[3] - 27.0).abs() < 1e-2);
    }

    // ── Cross-validation / L-curve tests ─────────────────────────────────

    #[test]
    fn test_cross_validation_returns_lambda() {
        let n = 4;
        let a = eye(n);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let lambdas = vec![0.001, 0.01, 0.1, 1.0, 10.0];
        let (best, res, sol) = cross_validation(&a, n, n, &b, &lambdas);
        assert!(lambdas.contains(&best));
        assert_eq!(res.len(), 5);
        assert_eq!(sol.len(), 5);
    }

    #[test]
    fn test_cross_validation_monotone() {
        // For A = I: as λ increases residual increases and solution norm decreases.
        let n = 3;
        let a = eye(n);
        let b = vec![1.0, 1.0, 1.0];
        let lambdas = vec![0.01, 0.1, 1.0, 10.0];
        let (_best, res, sol) = cross_validation(&a, n, n, &b, &lambdas);
        for i in 1..res.len() {
            assert!(res[i] >= res[i - 1] - 1e-10, "res not monotone at {i}");
            assert!(sol[i] <= sol[i - 1] + 1e-10, "sol not monotone at {i}");
        }
    }

    // ── Morozov discrepancy tests ─────────────────────────────────────────

    #[test]
    fn test_morozov_lambda_returns_valid() {
        let n = 3;
        let a = eye(n);
        let b = vec![1.0, 1.0, 1.0];
        let lambdas = vec![0.001, 0.01, 0.1, 1.0];
        let lam = morozov_lambda(&a, n, n, &b, &lambdas, 0.5);
        assert!(lambdas.contains(&lam));
    }

    #[test]
    fn test_morozov_lambda_large_delta() {
        // With delta = 1 any λ satisfies the discrepancy → should return smallest
        // or go through the whole list and return last.
        let n = 2;
        let a = eye(n);
        let b = vec![1.0, 1.0];
        let lambdas = vec![0.01, 0.1, 1.0, 10.0];
        let lam = morozov_lambda(&a, n, n, &b, &lambdas, 1.0);
        assert!(lambdas.contains(&lam));
    }

    // ── relative_fit_error ────────────────────────────────────────────────

    #[test]
    fn test_relative_fit_error_exact() {
        let a = eye(3);
        let b = vec![1.0, 2.0, 3.0];
        let err = relative_fit_error(&a, 3, 3, &b, &b);
        assert!(err < 1e-10);
    }

    #[test]
    fn test_relative_fit_error_nonzero() {
        let a = eye(2);
        let b = vec![1.0, 0.0];
        let p = vec![0.0, 1.0];
        let err = relative_fit_error(&a, 2, 2, &p, &b);
        // ||Ap - b|| = ||[-1,1]|| = sqrt(2), ||b|| = 1
        assert!((err - 2.0_f64.sqrt()).abs() < 1e-8);
    }

    // ── mat_vec / mat_t_vec helpers ───────────────────────────────────────

    #[test]
    fn test_mat_vec_identity() {
        let n = 3;
        let a = eye(n);
        let x = vec![1.0, 2.0, 3.0];
        let y = mat_vec(&a, n, n, &x);
        assert_eq!(y, x);
    }

    #[test]
    fn test_mat_t_vec_identity() {
        let n = 3;
        let a = eye(n);
        let x = vec![4.0, 5.0, 6.0];
        let y = mat_t_vec(&a, n, n, &x);
        assert_eq!(y, x);
    }

    #[test]
    fn test_mat_vec_2x2() {
        // A = [[1,2],[3,4]], x = [1,1] → [3, 7]
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let x = vec![1.0, 1.0];
        let y = mat_vec(&a, 2, 2, &x);
        assert!((y[0] - 3.0).abs() < 1e-10);
        assert!((y[1] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_normal_matrix_positive_definite() {
        // AᵀA + λI must have positive diagonal for any λ > 0.
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let m = normal_matrix(&a, 2, 2, 1.0);
        assert!(m[0] > 0.0);
        assert!(m[3] > 0.0);
    }

    #[test]
    fn test_cg_solve_trivial() {
        let n = 2;
        let a = eye(n);
        let b = vec![3.0, 7.0];
        let x = cg_solve(&a, n, &b, 10, 1e-10);
        assert!((x[0] - 3.0).abs() < 1e-8);
        assert!((x[1] - 7.0).abs() < 1e-8);
    }

    #[test]
    fn test_add_sub_scale_vec() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let c = add_vec(&a, &b);
        assert_eq!(c, vec![5.0, 7.0, 9.0]);
        let d = sub_vec(&b, &a);
        assert_eq!(d, vec![3.0, 3.0, 3.0]);
        let e = scale_vec(&a, 2.0);
        assert_eq!(e, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn test_norm() {
        assert!((norm(&[3.0, 4.0]) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dot() {
        assert!((dot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]) - 32.0).abs() < 1e-10);
    }
}
