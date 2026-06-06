//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    AdmmResult, AugLagResult, BarrierResult, KktCheck, LcpResult, ProxGradResult, QpResult,
    SocConstraint, SocpResult,
};

/// Dot product of two slices of equal length.
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
}
/// L2 norm of a vector.
pub(super) fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}
/// Matrix-vector product: y = A * x. A is row-major with `cols` columns.
pub(super) fn mat_vec(a: &[f64], x: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut y = vec![0.0; rows];
    for i in 0..rows {
        let mut s = 0.0;
        for j in 0..cols {
            s += a[i * cols + j] * x[j];
        }
        y[i] = s;
    }
    y
}
/// Transpose matrix-vector product: y = A^T * x.
pub(super) fn mat_t_vec(a: &[f64], x: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut y = vec![0.0; cols];
    for i in 0..rows {
        for j in 0..cols {
            y[j] += a[i * cols + j] * x[i];
        }
    }
    y
}
/// Row-major index helper.
#[inline]
pub(super) fn i_rc(r: usize, c: usize, cols: usize) -> usize {
    r * cols + c
}
/// Solve a general dense linear system via Gauss elimination with partial pivoting.
pub(super) fn gauss_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col * (n + 1) + col].abs();
        for row in (col + 1)..n {
            let v = aug[row * (n + 1) + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-15 {
            return None;
        }
        if max_row != col {
            for j in 0..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }
        let pivot = aug[col * (n + 1) + col];
        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / pivot;
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            x[i] -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] /= aug[i * (n + 1) + i];
    }
    Some(x)
}
/// Solve a small dense symmetric positive definite system via Cholesky.
/// Returns None if the matrix is not SPD.
pub(super) fn solve_spd(h: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = 0.0;
            for k in 0..j {
                s += l[i * n + k] * l[j * n + k];
            }
            if i == j {
                let diag = h[i * n + i] - s;
                if diag <= 0.0 {
                    return None;
                }
                l[i * n + j] = diag.sqrt();
            } else {
                l[i * n + j] = (h[i * n + j] - s) / l[j * n + j];
            }
        }
    }
    let mut y = b.to_vec();
    for i in 0..n {
        for j in 0..i {
            y[i] -= l[i * n + j] * y[j];
        }
        y[i] /= l[i * n + i];
    }
    let mut x = y;
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            x[i] -= l[j * n + i] * x[j];
        }
        x[i] /= l[i * n + i];
    }
    Some(x)
}
/// Quadratic programming solver using the active-set method.
///
/// Solves: min  0.5 * x^T H x + c^T x
///         s.t. A_ineq * x <= b_ineq
///
/// `h` is the Hessian (n x n, row-major, symmetric positive definite).
/// `c` is the linear cost vector (length n).
/// `a_ineq` is the inequality constraint matrix (m x n, row-major).
/// `b_ineq` is the inequality right-hand side (length m).
/// `max_iter` is the maximum number of iterations.
pub fn qp_active_set(
    h: &[f64],
    c: &[f64],
    a_ineq: &[f64],
    b_ineq: &[f64],
    n: usize,
    m: usize,
    max_iter: usize,
) -> QpResult {
    let mut x = match solve_spd(h, &c.iter().map(|ci| -ci).collect::<Vec<_>>(), n) {
        Some(sol) => sol,
        None => vec![0.0; n],
    };
    let mut active: Vec<bool> = vec![false; m];
    let mut multipliers = vec![0.0; m];
    let mut converged = false;
    for iter in 0..max_iter {
        let ax = mat_vec(a_ineq, &x, m, n);
        let mut most_violated = None;
        let mut max_violation = 1e-10;
        for i in 0..m {
            if !active[i] && ax[i] > b_ineq[i] + max_violation {
                max_violation = ax[i] - b_ineq[i];
                most_violated = Some(i);
            }
        }
        if let Some(idx) = most_violated {
            active[idx] = true;
        }
        let active_indices: Vec<usize> = (0..m).filter(|&i| active[i]).collect();
        let n_active = active_indices.len();
        if n_active == 0 {
            x = match solve_spd(h, &c.iter().map(|ci| -ci).collect::<Vec<_>>(), n) {
                Some(sol) => sol,
                None => x,
            };
        } else {
            let dim = n + n_active;
            let mut kkt = vec![0.0; dim * dim];
            let mut rhs = vec![0.0; dim];
            for i in 0..n {
                for j in 0..n {
                    kkt[i * dim + j] = h[i * n + j];
                }
            }
            for (k, &ai) in active_indices.iter().enumerate() {
                for j in 0..n {
                    let val = a_ineq[ai * n + j];
                    kkt[i_rc(n + k, j, dim)] = val;
                    kkt[i_rc(j, n + k, dim)] = val;
                }
            }
            for (i, ri) in rhs[..n].iter_mut().enumerate() {
                *ri = -c[i];
            }
            for (k, &ai) in active_indices.iter().enumerate() {
                rhs[n + k] = b_ineq[ai];
            }
            if let Some(sol) = gauss_solve(&kkt, &rhs, dim) {
                x[..n].copy_from_slice(&sol[..n]);
                for (k, &ai) in active_indices.iter().enumerate() {
                    multipliers[ai] = sol[n + k];
                }
            }
        }
        let mut removed = false;
        for &ai in &active_indices {
            if multipliers[ai] < -1e-10 {
                active[ai] = false;
                multipliers[ai] = 0.0;
                removed = true;
            }
        }
        if most_violated.is_none() && !removed {
            converged = true;
            if iter > 0 {
                break;
            }
        }
    }
    QpResult {
        x,
        multipliers,
        iterations: max_iter,
        converged,
    }
}
/// Solve a linear complementarity problem using Lemke's algorithm.
///
/// Find z >= 0 such that w = M*z + q >= 0 and z^T * w = 0.
///
/// `m_mat` is the LCP matrix (n x n, row-major).
/// `q` is the constant vector (length n).
/// `max_pivots` is the maximum number of pivot operations.
pub fn lcp_lemke(m_mat: &[f64], q: &[f64], n: usize, max_pivots: usize) -> LcpResult {
    if n == 0 {
        return LcpResult {
            z: vec![],
            w: vec![],
            pivots: 0,
            solved: true,
        };
    }
    if q.iter().all(|&qi| qi >= -1e-12) {
        return LcpResult {
            z: vec![0.0; n],
            w: q.to_vec(),
            pivots: 0,
            solved: true,
        };
    }
    let total_vars = 2 * n + 1;
    let mut tableau = vec![0.0; n * (total_vars + 1)];
    let cols = total_vars + 1;
    for i in 0..n {
        tableau[i * cols + i] = 1.0;
        for j in 0..n {
            tableau[i * cols + n + j] = -m_mat[i * n + j];
        }
        tableau[i * cols + 2 * n] = -1.0;
        tableau[i * cols + total_vars] = q[i];
    }
    let mut basis: Vec<usize> = (0..n).collect();
    let mut pivot_row = 0;
    let mut min_q = q[0];
    for (i, &q_i) in q.iter().enumerate().take(n).skip(1) {
        if q_i < min_q {
            min_q = q_i;
            pivot_row = i;
        }
    }
    let mut entering = 2 * n;
    let mut pivots = 0;
    let mut solved = false;
    for _step in 0..max_pivots {
        pivots += 1;
        let pivot_val = tableau[pivot_row * cols + entering];
        if pivot_val.abs() < 1e-15 {
            break;
        }
        let inv_pivot = 1.0 / pivot_val;
        for j in 0..cols {
            tableau[pivot_row * cols + j] *= inv_pivot;
        }
        for i in 0..n {
            if i == pivot_row {
                continue;
            }
            let factor = tableau[i * cols + entering];
            if factor.abs() < 1e-15 {
                continue;
            }
            for j in 0..cols {
                tableau[i * cols + j] -= factor * tableau[pivot_row * cols + j];
            }
        }
        let leaving = basis[pivot_row];
        basis[pivot_row] = entering;
        if leaving < n {
            entering = n + leaving;
        } else if leaving < 2 * n {
            entering = leaving - n;
        } else {
            solved = true;
            break;
        }
        let mut min_ratio = f64::MAX;
        let mut next_pivot_row = None;
        for i in 0..n {
            let col_val = tableau[i * cols + entering];
            if col_val > 1e-12 {
                let ratio = tableau[i * cols + total_vars] / col_val;
                if ratio < min_ratio {
                    min_ratio = ratio;
                    next_pivot_row = Some(i);
                }
            }
        }
        match next_pivot_row {
            Some(row) => pivot_row = row,
            None => break,
        }
    }
    let mut z = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n {
        let var = basis[i];
        let val = tableau[i * cols + total_vars].max(0.0);
        if var < n {
            w[var] = val;
        } else if var < 2 * n {
            z[var - n] = val;
        }
    }
    if solved {
        let mz = mat_vec(m_mat, &z, n, n);
        for (i, wi) in w.iter_mut().enumerate() {
            *wi = mz[i] + q[i];
            if *wi < 0.0 && *wi > -1e-10 {
                *wi = 0.0;
            }
        }
    }
    LcpResult {
        z,
        w,
        pivots,
        solved,
    }
}
/// Solve a second-order cone program via interior-point method.
///
/// min  f^T x
/// s.t. ||A_i x + b_i|| <= c_i^T x + d_i  for each cone constraint
///
/// Uses a barrier approach with Newton steps.
pub fn socp_solve(
    f: &[f64],
    constraints: &[SocConstraint],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> SocpResult {
    let mut x = vec![0.0; n];
    let mut converged = false;
    let mut _best_obj = f64::MAX;
    for constraint in constraints {
        let slack = constraint
            .c_vec
            .iter()
            .zip(x.iter())
            .map(|(ci, xi)| ci * xi)
            .sum::<f64>()
            + constraint.d_val;
        if slack <= 0.0 {
            let c_norm_sq: f64 = constraint.c_vec.iter().map(|ci| ci * ci).sum();
            if c_norm_sq > 1e-15 {
                let needed = 1.0 - slack;
                for (i, xi) in x.iter_mut().enumerate() {
                    *xi += needed * constraint.c_vec[i] / c_norm_sq;
                }
            }
        }
    }
    let mut mu = 10.0;
    for iter in 0..max_iter {
        let mut grad = f.to_vec();
        let mut all_feasible = true;
        for constraint in constraints {
            let ax = mat_vec(&constraint.a_mat, &x, constraint.cone_dim, n);
            let mut residual = vec![0.0; constraint.cone_dim];
            for (k, rk) in residual.iter_mut().enumerate() {
                *rk = ax[k] + constraint.b_vec[k];
            }
            let res_norm = norm(&residual);
            let slack = dot(&constraint.c_vec, &x) + constraint.d_val;
            let gap = slack - res_norm;
            if gap <= 0.0 {
                all_feasible = false;
                for (i, gi) in grad.iter_mut().enumerate() {
                    *gi -= constraint.c_vec[i];
                }
                continue;
            }
            let barrier_denom = slack * slack - res_norm * res_norm;
            if barrier_denom.abs() < 1e-15 {
                continue;
            }
            let inv_barrier = mu / barrier_denom;
            let at_res = mat_t_vec(&constraint.a_mat, &residual, constraint.cone_dim, n);
            for (i, gi) in grad.iter_mut().enumerate() {
                *gi += inv_barrier * (2.0 * at_res[i] - 2.0 * constraint.c_vec[i] * slack);
            }
        }
        let grad_norm = norm(&grad);
        if grad_norm < tol && all_feasible {
            converged = true;
            if iter > 0 {
                break;
            }
        }
        let step_size = 0.01 / (1.0 + grad_norm);
        for (i, xi) in x.iter_mut().enumerate() {
            *xi -= step_size * grad[i];
        }
        mu *= 0.95;
    }
    let objective = dot(f, &x);
    SocpResult {
        x,
        objective,
        iterations: max_iter,
        converged,
    }
}
/// Parameters for [`barrier_method`].
#[derive(Debug, Clone)]
pub struct BarrierMethodParams<'a> {
    /// Objective cost vector `f` (length `n`).
    pub f: &'a [f64],
    /// Inequality constraint matrix `A_ineq` (row-major, `m × n`).
    pub a_ineq: &'a [f64],
    /// Inequality constraint right-hand side `b_ineq` (length `m`).
    pub b_ineq: &'a [f64],
    /// Number of primal variables.
    pub n: usize,
    /// Number of inequality constraints.
    pub m: usize,
    /// Initial primal iterate (length `n`).
    pub x0: &'a [f64],
    /// Initial barrier parameter.
    pub mu_init: f64,
    /// Barrier reduction factor per outer iteration.
    pub mu_factor: f64,
    /// Maximum outer iterations.
    pub outer_iter: usize,
    /// Maximum inner (Newton) iterations.
    pub inner_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

/// Log-barrier method for inequality-constrained optimization.
///
/// min  f^T x
/// s.t. A x <= b
///
/// Uses the log-barrier approach: min f^T x - mu * sum(log(b_i - a_i^T x)).
pub fn barrier_method(p: BarrierMethodParams<'_>) -> BarrierResult {
    let BarrierMethodParams {
        f,
        a_ineq,
        b_ineq,
        n,
        m,
        x0,
        mu_init,
        mu_factor,
        outer_iter,
        inner_iter,
        tol,
    } = p;
    let mut x = x0.to_vec();
    let mut mu = mu_init;
    let mut converged = false;
    for _outer in 0..outer_iter {
        for _inner in 0..inner_iter {
            let mut grad = f.to_vec();
            let ax = mat_vec(a_ineq, &x, m, n);
            let mut feasible = true;
            for i in 0..m {
                let slack = b_ineq[i] - ax[i];
                if slack <= 1e-15 {
                    feasible = false;
                    for j in 0..n {
                        grad[j] -= a_ineq[i * n + j] * 100.0;
                    }
                } else {
                    let inv_slack = mu / slack;
                    for j in 0..n {
                        grad[j] += a_ineq[i * n + j] * inv_slack;
                    }
                }
            }
            if !feasible {
                let grad_norm = norm(&grad);
                let step_size = 0.001 / (1.0 + grad_norm);
                for (i, xi) in x.iter_mut().enumerate() {
                    *xi -= step_size * grad[i];
                }
                continue;
            }
            let mut hess = vec![0.0; n * n];
            for i in 0..m {
                let slack = b_ineq[i] - ax[i];
                let inv_slack_sq = mu / (slack * slack);
                for j in 0..n {
                    for k in 0..n {
                        hess[j * n + k] += a_ineq[i * n + j] * a_ineq[i * n + k] * inv_slack_sq;
                    }
                }
            }
            for i in 0..n {
                hess[i * n + i] += 1e-8;
            }
            let neg_grad: Vec<f64> = grad.iter().map(|g| -g).collect();
            let step = match solve_spd(&hess, &neg_grad, n) {
                Some(s) => s,
                None => {
                    let grad_norm = norm(&grad);
                    let s = 0.01 / (1.0 + grad_norm);
                    grad.iter().map(|g| -s * g).collect()
                }
            };
            let mut alpha = 1.0;
            for _ls in 0..20 {
                let x_new: Vec<f64> = (0..n).map(|i| x[i] + alpha * step[i]).collect();
                let ax_new = mat_vec(a_ineq, &x_new, m, n);
                let all_feasible = (0..m).all(|i| ax_new[i] < b_ineq[i] - 1e-15);
                if all_feasible {
                    break;
                }
                alpha *= 0.5;
            }
            for (i, xi) in x.iter_mut().enumerate() {
                *xi += alpha * step[i];
            }
            let grad_norm = norm(&grad);
            if grad_norm < tol {
                break;
            }
        }
        if mu * (m as f64) < tol {
            converged = true;
            break;
        }
        mu *= mu_factor;
    }
    let objective = dot(f, &x);
    BarrierResult {
        x,
        objective,
        mu,
        outer_iterations: outer_iter,
        converged,
    }
}
/// Parameters for [`augmented_lagrangian`].
#[derive(Debug, Clone)]
pub struct AugLagParams<'a> {
    /// Quadratic objective matrix `H` (row-major, `n × n`).
    pub h: &'a [f64],
    /// Linear cost vector `c` (length `n`).
    pub c: &'a [f64],
    /// Equality constraint matrix `A_eq` (row-major, `m × n`).
    pub a_eq: &'a [f64],
    /// Equality constraint right-hand side `b_eq` (length `m`).
    pub b_eq: &'a [f64],
    /// Number of primal variables.
    pub n: usize,
    /// Number of equality constraints.
    pub m: usize,
    /// Initial penalty parameter.
    pub rho_init: f64,
    /// Penalty growth factor per outer iteration.
    pub rho_factor: f64,
    /// Maximum outer iterations.
    pub outer_iter: usize,
    /// Maximum inner (linear solve) iterations.
    pub inner_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

/// Augmented Lagrangian method for equality-constrained optimization.
///
/// min  0.5 * x^T H x + c^T x
/// s.t. A_eq * x = b_eq
///
/// The augmented Lagrangian is:
/// L_rho(x, lambda) = f(x) + lambda^T (Ax - b) + (rho/2) ||Ax - b||^2
pub fn augmented_lagrangian(p: AugLagParams<'_>) -> AugLagResult {
    let AugLagParams {
        h,
        c,
        a_eq,
        b_eq,
        n,
        m,
        rho_init,
        rho_factor,
        outer_iter,
        inner_iter,
        tol,
    } = p;
    let mut x = vec![0.0; n];
    let mut lambda = vec![0.0; m];
    let mut rho = rho_init;
    let mut converged = false;
    let mut violation = f64::MAX;
    for _outer in 0..outer_iter {
        for _inner in 0..inner_iter {
            let hx = mat_vec(h, &x, n, n);
            let ax = mat_vec(a_eq, &x, m, n);
            let mut residual = vec![0.0; m];
            for (i, ri) in residual.iter_mut().enumerate() {
                *ri = ax[i] - b_eq[i];
            }
            let at_lambda = mat_t_vec(a_eq, &lambda, m, n);
            let at_res = mat_t_vec(a_eq, &residual, m, n);
            let mut grad = vec![0.0; n];
            for (i, gi) in grad.iter_mut().enumerate() {
                *gi = hx[i] + c[i] + at_lambda[i] + rho * at_res[i];
            }
            let mut hess = h.to_vec();
            for i in 0..n {
                for j in 0..n {
                    for k in 0..m {
                        hess[i * n + j] += rho * a_eq[k * n + i] * a_eq[k * n + j];
                    }
                }
            }
            let neg_grad: Vec<f64> = grad.iter().map(|g| -g).collect();
            let step = match solve_spd(&hess, &neg_grad, n) {
                Some(s) => s,
                None => break,
            };
            for (i, xi) in x.iter_mut().enumerate() {
                *xi += step[i];
            }
            if norm(&grad) < tol * 0.1 {
                break;
            }
        }
        let ax = mat_vec(a_eq, &x, m, n);
        let mut residual = vec![0.0; m];
        for (i, ri) in residual.iter_mut().enumerate() {
            *ri = ax[i] - b_eq[i];
        }
        violation = norm(&residual);
        for (i, li) in lambda.iter_mut().enumerate() {
            *li += rho * residual[i];
        }
        if violation < tol {
            converged = true;
            break;
        }
        rho *= rho_factor;
    }
    AugLagResult {
        x,
        lambda,
        rho,
        outer_iterations: outer_iter,
        violation,
        converged,
    }
}
/// ADMM for consensus optimization.
///
/// Solves: min f(x) + g(z)  s.t. x = z
///
/// where f(x) = 0.5 x^T H x + c^T x and g(z) = indicator of {z : A z <= b}.
///
/// `h` is the objective Hessian (n x n).
/// `c` is the objective linear term (n).
/// `a_ineq` / `b_ineq` define inequality constraints on z (m x n / m).
/// `rho` is the penalty parameter.
pub fn admm_consensus(
    h: &[f64],
    c: &[f64],
    a_ineq: &[f64],
    b_ineq: &[f64],
    n: usize,
    m: usize,
    rho: f64,
    max_iter: usize,
    tol: f64,
) -> AdmmResult {
    let mut x = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut y = vec![0.0; n];
    let mut converged = false;
    let mut primal_res = f64::MAX;
    let mut dual_res = f64::MAX;
    let mut iterations = 0;
    let mut h_aug = h.to_vec();
    for i in 0..n {
        h_aug[i * n + i] += rho;
    }
    for iter in 0..max_iter {
        iterations = iter + 1;
        let mut rhs = vec![0.0; n];
        for (i, ri) in rhs.iter_mut().enumerate() {
            *ri = -c[i] + rho * (z[i] - y[i]);
        }
        let x_new = match solve_spd(&h_aug, &rhs, n) {
            Some(s) => s,
            None => break,
        };
        let z_old = z.clone();
        let mut z_target: Vec<f64> = (0..n).map(|i| x_new[i] + y[i]).collect();
        for _proj_iter in 0..10 {
            let az = mat_vec(a_ineq, &z_target, m, n);
            let mut violated = false;
            for i in 0..m {
                if az[i] > b_ineq[i] {
                    violated = true;
                    let row = &a_ineq[i * n..(i + 1) * n];
                    let row_norm_sq: f64 = row.iter().map(|r| r * r).sum();
                    if row_norm_sq > 1e-15 {
                        let excess = az[i] - b_ineq[i];
                        for (j, zt) in z_target.iter_mut().enumerate() {
                            *zt -= (excess / row_norm_sq) * row[j];
                        }
                    }
                }
            }
            if !violated {
                break;
            }
        }
        z = z_target;
        for (i, yi) in y.iter_mut().enumerate() {
            *yi += x_new[i] - z[i];
        }
        primal_res = norm(&(0..n).map(|i| x_new[i] - z[i]).collect::<Vec<_>>());
        dual_res = rho * norm(&(0..n).map(|i| z[i] - z_old[i]).collect::<Vec<_>>());
        x = x_new;
        if primal_res < tol && dual_res < tol {
            converged = true;
            break;
        }
    }
    AdmmResult {
        x,
        z,
        y,
        iterations,
        primal_residual: primal_res,
        dual_residual: dual_res,
        converged,
    }
}
/// Proximal operator for the L1 norm (soft thresholding).
///
/// prox_{t * ||.||_1}(v) = sign(v_i) * max(|v_i| - t, 0)
pub fn prox_l1(v: &[f64], t: f64) -> Vec<f64> {
    v.iter()
        .map(|&vi| {
            let abs_vi = vi.abs();
            if abs_vi > t {
                vi.signum() * (abs_vi - t)
            } else {
                0.0
            }
        })
        .collect()
}
/// Proximal operator for the L2 norm (block soft thresholding).
///
/// prox_{t * ||.||_2}(v) = (1 - t / ||v||)_+ * v
pub fn prox_l2(v: &[f64], t: f64) -> Vec<f64> {
    let v_norm = norm(v);
    if v_norm <= t {
        vec![0.0; v.len()]
    } else {
        let scale = 1.0 - t / v_norm;
        v.iter().map(|&vi| scale * vi).collect()
    }
}
/// Proximal operator for the indicator of the non-negative orthant.
///
/// prox_{I_{R+}}(v) = max(v, 0)
pub fn prox_nonneg(v: &[f64]) -> Vec<f64> {
    v.iter().map(|&vi| vi.max(0.0)).collect()
}
/// Proximal operator for the indicator of a box constraint \[lo, hi\].
///
/// prox_{I_{\[lo,hi\]}}(v) = clamp(v, lo, hi)
pub fn prox_box(v: &[f64], lo: f64, hi: f64) -> Vec<f64> {
    v.iter().map(|&vi| vi.max(lo).min(hi)).collect()
}
/// Proximal operator for the indicator of a second-order cone.
///
/// The SOC is {(t, x) : ||x|| <= t}.
/// Input v = (v_t, v_x...) where v_t is the first element.
pub fn prox_soc(v: &[f64]) -> Vec<f64> {
    if v.is_empty() {
        return vec![];
    }
    if v.len() == 1 {
        return vec![v[0].max(0.0)];
    }
    let vt = v[0];
    let vx = &v[1..];
    let vx_norm = norm(vx);
    if vx_norm <= vt {
        v.to_vec()
    } else if vx_norm <= -vt {
        vec![0.0; v.len()]
    } else {
        let scale = 0.5 * (1.0 + vt / vx_norm);
        let t_proj = scale * vx_norm;
        let mut result = vec![t_proj];
        for &xi in vx {
            result.push(scale * xi);
        }
        result
    }
}
/// Proximal operator for the elastic net penalty: alpha * ||x||_1 + (1-alpha)/2 * ||x||_2^2.
pub fn prox_elastic_net(v: &[f64], t: f64, alpha: f64) -> Vec<f64> {
    let l1_part = alpha * t;
    let l2_scale = 1.0 / (1.0 + (1.0 - alpha) * t);
    v.iter()
        .map(|&vi| {
            let soft = if vi.abs() > l1_part {
                vi.signum() * (vi.abs() - l1_part)
            } else {
                0.0
            };
            l2_scale * soft
        })
        .collect()
}
/// Moreau envelope evaluation: M_{t*f}(v) = min_x { f(x) + 1/(2t) ||x - v||^2 }.
/// For f = ||.||_1, the Moreau envelope at v is the Huber function.
pub fn moreau_envelope_l1(v: &[f64], t: f64) -> f64 {
    v.iter()
        .map(|&vi| {
            let abs_vi = vi.abs();
            if abs_vi <= t {
                vi * vi / (2.0 * t)
            } else {
                abs_vi - t / 2.0
            }
        })
        .sum()
}
/// Check the Linear Independence Constraint Qualification (LICQ).
///
/// LICQ holds at x if the gradients of all active inequality constraints
/// and all equality constraints are linearly independent.
///
/// `active_gradients` is a list of constraint gradient vectors (each length n).
/// Returns true if the gradients are linearly independent (via rank check).
pub fn check_licq(active_gradients: &[Vec<f64>], n: usize, tol: f64) -> bool {
    let m = active_gradients.len();
    if m == 0 {
        return true;
    }
    if m > n {
        return false;
    }
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m);
    for grad in active_gradients {
        let mut v = grad.clone();
        for b in &basis {
            let proj = dot(&v, b) / dot(b, b).max(1e-30);
            for (i, vi) in v.iter_mut().enumerate() {
                *vi -= proj * b[i];
            }
        }
        let v_norm = norm(&v);
        if v_norm > tol {
            for vi in &mut v {
                *vi /= v_norm;
            }
            basis.push(v);
        }
    }
    basis.len() == m
}
/// Check the Mangasarian-Fromovitz Constraint Qualification (MFCQ).
///
/// MFCQ holds at x if:
/// 1. The equality constraint gradients are linearly independent.
/// 2. There exists a direction d such that:
///    - For all active inequalities: grad(g_i)^T * d < 0
///    - For all equalities: grad(h_j)^T * d = 0
///
/// `eq_gradients` - equality constraint gradients (each length n).
/// `ineq_gradients` - active inequality constraint gradients (each length n).
///
/// Uses a simple feasibility check: try d in the null space of equalities
/// that has negative inner product with all active inequality gradients.
pub fn check_mfcq(
    eq_gradients: &[Vec<f64>],
    ineq_gradients: &[Vec<f64>],
    n: usize,
    tol: f64,
) -> bool {
    if !eq_gradients.is_empty() && !check_licq(eq_gradients, n, tol) {
        return false;
    }
    if ineq_gradients.is_empty() {
        return true;
    }
    let mut d = vec![0.0; n];
    for grad in ineq_gradients {
        for (i, di) in d.iter_mut().enumerate() {
            *di -= grad[i];
        }
    }
    for eq_grad in eq_gradients {
        let eq_norm_sq = dot(eq_grad, eq_grad);
        if eq_norm_sq > 1e-15 {
            let proj = dot(&d, eq_grad) / eq_norm_sq;
            for (i, di) in d.iter_mut().enumerate() {
                *di -= proj * eq_grad[i];
            }
        }
    }
    let d_norm = norm(&d);
    if d_norm < tol {
        return false;
    }
    for di in &mut d {
        *di /= d_norm;
    }
    for grad in ineq_gradients {
        if dot(&d, grad) >= -tol {
            return false;
        }
    }
    for eq_grad in eq_gradients {
        if dot(&d, eq_grad).abs() > tol {
            return false;
        }
    }
    true
}
/// Check KKT conditions for a QP with inequality constraints.
///
/// min 0.5 x^T H x + c^T x  s.t. A x <= b
///
/// KKT conditions:
/// 1. Stationarity: H*x + c + A^T * mu = 0
/// 2. Primal feasibility: A*x <= b
/// 3. Dual feasibility: mu >= 0
/// 4. Complementary slackness: mu_i * (a_i^T x - b_i) = 0
pub fn check_kkt(
    h: &[f64],
    c: &[f64],
    a_ineq: &[f64],
    b_ineq: &[f64],
    x: &[f64],
    mu: &[f64],
    n: usize,
    m: usize,
    tol: f64,
) -> KktCheck {
    let hx = mat_vec(h, x, n, n);
    let at_mu = mat_t_vec(a_ineq, mu, m, n);
    let mut stationarity_vec = vec![0.0; n];
    for (i, si) in stationarity_vec.iter_mut().enumerate() {
        *si = hx[i] + c[i] + at_mu[i];
    }
    let stationarity = norm(&stationarity_vec);
    let ax = mat_vec(a_ineq, x, m, n);
    let mut pf = 0.0_f64;
    for (&axi, &bi) in ax.iter().zip(b_ineq.iter()) {
        let violation = axi - bi;
        if violation > pf {
            pf = violation;
        }
    }
    let mut df = 0.0_f64;
    for &mui in &mu[..m] {
        if -mui > df {
            df = -mui;
        }
    }
    let mut cs = 0.0_f64;
    for ((&bi, &axi), &mui) in b_ineq.iter().zip(ax.iter()).zip(mu.iter()) {
        let slack = bi - axi;
        let violation = (mui * slack).abs();
        if violation > cs {
            cs = violation;
        }
    }
    let satisfied = stationarity < tol && pf < tol && df < tol && cs < tol;
    KktCheck {
        stationarity,
        primal_feasibility: pf,
        dual_feasibility: df,
        complementarity: cs,
        satisfied,
    }
}
/// Proximal gradient method for composite optimization.
///
/// min f(x) + g(x)
///
/// where f is smooth (quadratic: 0.5 x^T H x + c^T x) and g is the L1 norm
/// (g(x) = lambda * ||x||_1).
///
/// Uses ISTA (Iterative Shrinkage-Thresholding Algorithm).
pub fn proximal_gradient_l1(
    h: &[f64],
    c: &[f64],
    n: usize,
    lambda: f64,
    step_size: f64,
    max_iter: usize,
    tol: f64,
) -> ProxGradResult {
    let mut x = vec![0.0; n];
    let mut converged = false;
    let mut iterations = 0;
    for iter in 0..max_iter {
        iterations = iter + 1;
        let hx = mat_vec(h, &x, n, n);
        let mut grad = vec![0.0; n];
        for (i, gi) in grad.iter_mut().enumerate() {
            *gi = hx[i] + c[i];
        }
        let mut x_grad: Vec<f64> = (0..n).map(|i| x[i] - step_size * grad[i]).collect();
        let x_new = prox_l1(&x_grad, lambda * step_size);
        let diff_norm = norm(&(0..n).map(|i| x_new[i] - x[i]).collect::<Vec<_>>());
        x = x_new;
        let _unused = &mut x_grad;
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
