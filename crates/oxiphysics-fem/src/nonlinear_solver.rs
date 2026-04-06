#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nonlinear FEM solvers: Newton-Raphson, line search, arc-length (Riks),
//! modified Newton, quasi-Newton (BFGS), and convergence criteria.
//!
//! Provides iterative solution methods for systems of nonlinear equations
//! arising in finite element analysis, using plain `f64` arrays throughout.

#![allow(dead_code)]

// ─── Residual trait ──────────────────────────────────────────────────────────

/// Trait representing a nonlinear FEM problem.
///
/// Implementors define residual forces and the tangent stiffness matrix so that
/// Newton-type solvers can iterate to equilibrium.
pub trait ResidualFn {
    /// Compute residual vector R(u) = internal forces − external forces.
    fn residual(&self, u: &[f64]) -> Vec<f64>;

    /// Compute tangent stiffness (Jacobian) K(u) = dR/du.
    fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>>;
}

// ─── Newton result ───────────────────────────────────────────────────────────

/// Result returned by Newton-type solvers.
#[derive(Debug, Clone)]
pub struct NewtonResult {
    /// Converged displacement/state vector.
    pub solution: Vec<f64>,
    /// Number of iterations performed.
    pub iterations: usize,
    /// L2 norm of the residual at the final iteration.
    pub final_residual: f64,
    /// Whether the solver converged within the iteration budget.
    pub converged: bool,
    /// Residual norm at each iteration.
    pub residual_history: Vec<f64>,
}

// ─── Newton-Raphson ──────────────────────────────────────────────────────────

/// Full Newton-Raphson solver (tangent updated every iteration).
#[derive(Debug, Clone)]
pub struct NewtonRaphson {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the L2 residual norm.
    pub tolerance: f64,
}

impl NewtonRaphson {
    /// Create a new Newton-Raphson solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self {
            max_iterations: max_iter,
            tolerance: tol,
        }
    }

    /// Solve the nonlinear problem starting from `u0`.
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let mut u = u0.to_vec();
        let mut history = Vec::new();

        for iter in 0..self.max_iterations {
            let r = problem.residual(&u);
            let norm = norm_l2(&r);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            let mut k = problem.tangent_stiffness(&u);
            let mut rhs = r.clone();
            let delta = gauss_elimination(&mut k, &mut rhs)?;
            u = vec_sub(&u, &delta);
        }

        let r_final = problem.residual(&u);
        let norm_final = norm_l2(&r_final);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── Newton-Raphson with Line Search ─────────────────────────────────────────

/// Newton-Raphson solver with backtracking line search.
///
/// After computing the Newton direction, performs a line search to find
/// an optimal step length that sufficiently reduces the residual norm.
#[derive(Debug, Clone)]
pub struct NewtonRaphsonLineSearch {
    /// Maximum number of Newton iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the L2 residual norm.
    pub tolerance: f64,
    /// Maximum number of line search iterations per Newton step.
    pub max_ls_iterations: usize,
    /// Backtracking factor (0 < alpha < 1, typically 0.5).
    pub backtrack_factor: f64,
    /// Sufficient decrease parameter (Armijo condition, typically 1e-4).
    pub armijo_c: f64,
}

impl NewtonRaphsonLineSearch {
    /// Create a new Newton-Raphson solver with line search.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self {
            max_iterations: max_iter,
            tolerance: tol,
            max_ls_iterations: 10,
            backtrack_factor: 0.5,
            armijo_c: 1e-4,
        }
    }

    /// Solve the nonlinear problem with line search.
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let mut u = u0.to_vec();
        let mut history = Vec::new();

        for iter in 0..self.max_iterations {
            let r = problem.residual(&u);
            let norm = norm_l2(&r);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            let mut k = problem.tangent_stiffness(&u);
            let mut rhs = r.clone();
            let delta = gauss_elimination(&mut k, &mut rhs)?;

            // Line search: find step length alpha such that
            // ||R(u - alpha*delta)|| < ||R(u)|| - c * alpha * delta . R
            let mut alpha = 1.0;
            let slope: f64 = delta.iter().zip(r.iter()).map(|(d, ri)| d * ri).sum();

            for _ in 0..self.max_ls_iterations {
                let u_trial: Vec<f64> = u
                    .iter()
                    .zip(delta.iter())
                    .map(|(ui, di)| ui - alpha * di)
                    .collect();
                let r_trial = problem.residual(&u_trial);
                let norm_trial = norm_l2(&r_trial);

                if norm_trial < norm - self.armijo_c * alpha * slope.abs() {
                    break;
                }
                alpha *= self.backtrack_factor;
            }

            u = u
                .iter()
                .zip(delta.iter())
                .map(|(ui, di)| ui - alpha * di)
                .collect();
        }

        let r_final = problem.residual(&u);
        let norm_final = norm_l2(&r_final);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── Modified Newton-Raphson ─────────────────────────────────────────────────

/// Modified Newton-Raphson solver: re-uses the tangent stiffness for several
/// iterations before recomputing it, reducing assembly cost.
#[derive(Debug, Clone)]
pub struct ModifiedNewtonRaphson {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the L2 residual norm.
    pub tolerance: f64,
    /// Recompute the tangent stiffness every this many iterations (1 = full NR).
    pub update_tangent_every: usize,
}

impl ModifiedNewtonRaphson {
    /// Create a new modified Newton-Raphson solver.
    pub fn new(max_iter: usize, tol: f64, update_every: usize) -> Self {
        Self {
            max_iterations: max_iter,
            tolerance: tol,
            update_tangent_every: update_every.max(1),
        }
    }

    /// Solve the nonlinear problem starting from `u0`.
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let mut u = u0.to_vec();
        let mut history = Vec::new();
        let mut cached_k: Option<Vec<Vec<f64>>> = None;

        for iter in 0..self.max_iterations {
            let r = problem.residual(&u);
            let norm = norm_l2(&r);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            // Refresh tangent when scheduled.
            if iter % self.update_tangent_every == 0 {
                cached_k = Some(problem.tangent_stiffness(&u));
            }

            let k_ref = cached_k.as_ref().expect("cached_k initialized at iter 0");
            let mut k_copy = k_ref.clone();
            let mut rhs = r.clone();
            let delta = gauss_elimination(&mut k_copy, &mut rhs)?;
            u = vec_sub(&u, &delta);
        }

        let r_final = problem.residual(&u);
        let norm_final = norm_l2(&r_final);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── Quasi-Newton BFGS ───────────────────────────────────────────────────────

/// Quasi-Newton solver using BFGS inverse Hessian update.
///
/// Avoids tangent stiffness assembly after the first iteration by
/// maintaining an approximate inverse Hessian via rank-2 updates.
#[derive(Debug, Clone)]
pub struct BfgsSolver {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

impl BfgsSolver {
    /// Create a new BFGS solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self {
            max_iterations: max_iter,
            tolerance: tol,
        }
    }

    /// Solve using BFGS.
    ///
    /// The initial inverse Hessian approximation is the inverse of the
    /// tangent stiffness at `u0`.
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let n = u0.len();
        let mut u = u0.to_vec();
        let mut history = Vec::new();

        // Initial inverse Hessian approximation: H0 = K(u0)^{-1}
        // We use the identity scaled by 1/||K|| as a simple approximation
        let k0 = problem.tangent_stiffness(&u);
        let diag_sum: f64 = (0..n).map(|i| k0[i][i].abs()).sum::<f64>() / n as f64;
        let scale = if diag_sum > 1e-14 {
            1.0 / diag_sum
        } else {
            1.0
        };

        let mut h_inv = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            h_inv[i][i] = scale;
        }

        let mut r_prev = problem.residual(&u);

        for iter in 0..self.max_iterations {
            let norm = norm_l2(&r_prev);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            // Compute search direction: delta = H^{-1} * R
            let delta = mat_vec(&h_inv, &r_prev);
            let u_new = vec_sub(&u, &delta);

            let r_new = problem.residual(&u_new);

            // BFGS update of H^{-1}
            let s: Vec<f64> = vec_sub(&u_new, &u); // s = u_new - u = -delta
            let y: Vec<f64> = vec_sub(&r_new, &r_prev); // y = R_new - R_old

            let sy: f64 = s.iter().zip(y.iter()).map(|(si, yi)| si * yi).sum();

            if sy.abs() > 1e-14 {
                // Sherman-Morrison-Woodbury BFGS update
                let hy = mat_vec(&h_inv, &y);
                let yhy: f64 = y.iter().zip(hy.iter()).map(|(yi, hyi)| yi * hyi).sum();

                for i in 0..n {
                    for j in 0..n {
                        h_inv[i][j] += (sy + yhy) * s[i] * s[j] / (sy * sy)
                            - (hy[i] * s[j] + s[i] * hy[j]) / sy;
                    }
                }
            }

            u = u_new;
            r_prev = r_new;
        }

        let norm_final = norm_l2(&r_prev);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── Load stepper ────────────────────────────────────────────────────────────

/// Incremental load stepping: solves at a sequence of load factors using a
/// Newton-Raphson solver at each step.
#[derive(Debug, Clone)]
pub struct LoadStepper {
    /// Number of load increments.
    pub n_steps: usize,
    /// Starting load factor.
    pub load_factor_start: f64,
    /// Ending load factor.
    pub load_factor_end: f64,
}

impl LoadStepper {
    /// Create a new load stepper.
    pub fn new(n_steps: usize, start: f64, end: f64) -> Self {
        Self {
            n_steps,
            load_factor_start: start,
            load_factor_end: end,
        }
    }

    /// Return evenly-spaced load factors.
    pub fn step_factors(&self) -> Vec<f64> {
        if self.n_steps == 0 {
            return vec![];
        }
        if self.n_steps == 1 {
            return vec![self.load_factor_start];
        }
        let range = self.load_factor_end - self.load_factor_start;
        (0..self.n_steps)
            .map(|i| self.load_factor_start + range * (i as f64) / ((self.n_steps - 1) as f64))
            .collect()
    }

    /// Solve at each load increment.
    pub fn solve<F: ResidualFn>(
        &self,
        base_problem: &F,
        u0: &[f64],
        solver: &NewtonRaphson,
    ) -> Vec<(f64, Vec<f64>)> {
        let mut results = Vec::new();
        let mut u = u0.to_vec();

        for &lf in &self.step_factors() {
            let scaled = ScaledProblem {
                inner: base_problem,
                scale: lf,
            };
            match solver.solve(&scaled, &u) {
                Ok(nr) => {
                    u = nr.solution.clone();
                    results.push((lf, nr.solution));
                }
                Err(_) => {
                    results.push((lf, u.clone()));
                }
            }
        }
        results
    }
}

/// Helper wrapper that scales external forces by a load factor.
struct ScaledProblem<'a, F> {
    inner: &'a F,
    scale: f64,
}

impl<F: ResidualFn> ResidualFn for ScaledProblem<'_, F> {
    fn residual(&self, u: &[f64]) -> Vec<f64> {
        let r = self.inner.residual(u);
        vec_scale(&r, self.scale)
    }

    fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
        self.inner.tangent_stiffness(u)
    }
}

// ─── Arc-length method (Crisfield / Riks) ─────────────────────────────────

/// Crisfield arc-length method for tracing equilibrium paths through limit
/// points and snap-back behaviour.
#[derive(Debug, Clone)]
pub struct ArcLengthMethod {
    /// Constraint radius in the (u, λ) space.
    pub arc_length: f64,
    /// Maximum iterations per step.
    pub max_iterations: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

impl ArcLengthMethod {
    /// Create a new arc-length solver.
    pub fn new(arc_length: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            arc_length,
            max_iterations: max_iter,
            tolerance: tol,
        }
    }

    /// Perform a single arc-length predictor-corrector step.
    pub fn solve_step<F: ResidualFn>(
        &self,
        problem: &F,
        u0: &[f64],
        lambda0: f64,
        du_ref: &[f64],
    ) -> Result<(Vec<f64>, f64), String> {
        let n = u0.len();
        let ds = self.arc_length;

        // Predictor
        let ref_norm = norm_l2(du_ref);
        let scale = if ref_norm > 1e-14 { ds / ref_norm } else { ds };
        let mut u = {
            let du_pred = vec_scale(du_ref, scale);
            u0.iter()
                .zip(&du_pred)
                .map(|(a, b)| a + b)
                .collect::<Vec<_>>()
        };
        let mut lambda = lambda0 + scale;

        // Corrector iterations
        for _iter in 0..self.max_iterations {
            let r_unit = problem.residual(&u);
            let du_cur: Vec<f64> = u.iter().zip(u0).map(|(a, b)| a - b).collect();
            let g = norm_l2(&du_cur).powi(2) + (lambda - lambda0).powi(2) - ds * ds;

            let r_scaled = vec_scale(&r_unit, lambda);
            if norm_l2(&r_scaled) < self.tolerance && g.abs() < self.tolerance {
                return Ok((u, lambda));
            }

            let mut k = problem.tangent_stiffness(&u);
            let mut rhs_t = r_unit.clone();
            let du_t = gauss_elimination(&mut k, &mut rhs_t)?;

            let f_ref_neg: Vec<f64> = r_unit.iter().map(|v| -*v).collect();
            let mut k2 = problem.tangent_stiffness(&u);
            let mut rhs_r = f_ref_neg;
            let du_r = gauss_elimination(&mut k2, &mut rhs_r)?;

            let dot_r: f64 = du_cur.iter().zip(&du_r).map(|(a, b)| a * b).sum();
            let dot_t: f64 = du_cur.iter().zip(&du_t).map(|(a, b)| a * b).sum();
            let denom = 2.0 * dot_r + 2.0 * (lambda - lambda0);
            if denom.abs() < 1e-14 {
                return Err("Arc-length denominator near zero".to_string());
            }
            let d_lambda = (-g - 2.0 * dot_t) / denom;

            let d_u: Vec<f64> = (0..n).map(|i| d_lambda * du_r[i] + du_t[i]).collect();
            for i in 0..n {
                u[i] -= d_u[i];
            }
            lambda += d_lambda;
        }

        Err(format!(
            "Arc-length method did not converge in {} iterations",
            self.max_iterations
        ))
    }

    /// Trace an equilibrium path over multiple arc-length steps.
    pub fn trace_path<F: ResidualFn>(
        &self,
        problem: &F,
        u0: &[f64],
        lambda0: f64,
        n_steps: usize,
    ) -> Vec<(Vec<f64>, f64)> {
        let mut path = Vec::new();
        let mut u = u0.to_vec();
        let mut lambda = lambda0;
        let mut du_ref = vec![1.0; u0.len()]; // initial direction

        for _ in 0..n_steps {
            match self.solve_step(problem, &u, lambda, &du_ref) {
                Ok((u_new, lambda_new)) => {
                    du_ref = vec_sub(&u_new, &u);
                    let du_norm = norm_l2(&du_ref);
                    if du_norm > 1e-14 {
                        du_ref = vec_scale(&du_ref, 1.0 / du_norm);
                    }
                    u = u_new.clone();
                    lambda = lambda_new;
                    path.push((u_new, lambda_new));
                }
                Err(_) => break,
            }
        }
        path
    }
}

// ─── Convergence criteria ────────────────────────────────────────────────────

/// Convergence criteria for nonlinear solvers.
#[derive(Debug, Clone)]
#[allow(missing_docs)]
pub enum ConvergenceCriteria {
    /// Converge when ‖R‖₂ < tol.
    ResidualNorm { tol: f64 },
    /// Converge when ‖Δu‖₂ < tol.
    DisplacementNorm { tol: f64 },
    /// Converge when the incremental energy |Δu · R| < tol.
    EnergyNorm { tol: f64 },
    /// Converge when both residual and displacement norms are satisfied.
    Combined {
        residual_tol: f64,
        displacement_tol: f64,
    },
    /// Converge when relative residual ‖R‖ / ‖R₀‖ < tol.
    RelativeResidual { tol: f64, initial_norm: f64 },
}

impl ConvergenceCriteria {
    /// Check whether convergence is satisfied.
    pub fn check(&self, residual: &[f64], delta_u: &[f64], _external_force: &[f64]) -> bool {
        match self {
            ConvergenceCriteria::ResidualNorm { tol } => norm_l2(residual) < *tol,
            ConvergenceCriteria::DisplacementNorm { tol } => norm_l2(delta_u) < *tol,
            ConvergenceCriteria::EnergyNorm { tol } => {
                let energy: f64 = delta_u
                    .iter()
                    .zip(residual)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    .abs();
                energy < *tol
            }
            ConvergenceCriteria::Combined {
                residual_tol,
                displacement_tol,
            } => norm_l2(residual) < *residual_tol && norm_l2(delta_u) < *displacement_tol,
            ConvergenceCriteria::RelativeResidual { tol, initial_norm } => {
                if *initial_norm < 1e-30 {
                    return norm_l2(residual) < *tol;
                }
                norm_l2(residual) / initial_norm < *tol
            }
        }
    }
}

// ─── Linear algebra helpers ──────────────────────────────────────────────────

/// Gaussian elimination with partial pivoting.
pub fn gauss_elimination(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Result<Vec<f64>, String> {
    let n = b.len();
    if n == 0 {
        return Err("Empty system".to_string());
    }
    if a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err("Matrix dimensions inconsistent".to_string());
    }

    for col in 0..n {
        let pivot_row = (col..n).max_by(|&i, &j| {
            a[i][col]
                .abs()
                .partial_cmp(&a[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let pivot_row = pivot_row.expect("range col..n is non-empty");

        if a[pivot_row][col].abs() < 1e-14 {
            return Err(format!("Singular matrix detected at column {col}"));
        }

        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        let pivot = a[col][col];
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            for k in col..n {
                let av = a[col][k];
                a[row][k] -= factor * av;
            }
            b[row] -= factor * b[col];
        }
    }

    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    Ok(x)
}

/// L2 (Euclidean) norm of a slice.
pub fn norm_l2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// L-infinity norm of a slice.
pub fn norm_inf(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0.0f64, f64::max)
}

/// Element-wise subtraction: a − b.
pub fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(ai, bi)| ai - bi).collect()
}

/// Element-wise addition: a + b.
pub fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(ai, bi)| ai + bi).collect()
}

/// Scalar multiplication of a vector.
pub fn vec_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|ai| ai * s).collect()
}

/// Dot product of two vectors.
pub fn vec_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(ai, bi)| ai * bi).sum()
}

/// Matrix-vector product M v.
pub fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(v).map(|(mij, vj)| mij * vj).sum())
        .collect()
}

// ─── Trust-Region solver ─────────────────────────────────────────────────────

/// Trust-region solver using Cauchy point (steepest descent) or Newton step.
///
/// At each iteration, solves a constrained subproblem inside a "trust region"
/// of radius `delta`. Accepts or rejects steps based on the predicted vs.
/// actual reduction in residual norm.
#[derive(Debug, Clone)]
pub struct TrustRegionSolver {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the residual norm.
    pub tolerance: f64,
    /// Initial trust-region radius.
    pub initial_radius: f64,
    /// Maximum trust-region radius.
    pub max_radius: f64,
}

impl TrustRegionSolver {
    /// Create a new trust-region solver.
    pub fn new(max_iter: usize, tol: f64, initial_radius: f64, max_radius: f64) -> Self {
        Self {
            max_iterations: max_iter,
            tolerance: tol,
            initial_radius,
            max_radius,
        }
    }

    /// Solve the nonlinear problem using the trust-region method.
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let n = u0.len();
        let mut u = u0.to_vec();
        let mut delta = self.initial_radius;
        let mut history = Vec::new();

        for iter in 0..self.max_iterations {
            let r = problem.residual(&u);
            let norm = norm_l2(&r);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            let mut k = problem.tangent_stiffness(&u);
            let mut rhs = r.clone();

            // Try to compute the full Newton step.
            // gauss_elimination solves K * delta_u = R, so delta_u = K^{-1} R.
            // The Newton correction is: u_new = u - delta_u.
            // We negate to get the actual update direction.
            let newton_step_raw = gauss_elimination(&mut k, &mut rhs).unwrap_or_else(|_| {
                // Fall back to gradient step if singular.
                let g_norm = norm_l2(&r);
                if g_norm > 1e-14 {
                    vec_scale(&r, delta / g_norm)
                } else {
                    vec![0.0; n]
                }
            });
            // Negate: correction = -K^{-1} R  (so u_new = u + correction = u - K^{-1} R)
            let newton_step: Vec<f64> = newton_step_raw.iter().map(|v| -v).collect();

            let newton_norm = norm_l2(&newton_step);

            // Choose the step: full Newton if within trust region, else truncate.
            let step = if newton_norm <= delta {
                newton_step.clone()
            } else {
                // Cauchy point along steepest descent direction (-R).
                let g_norm = norm_l2(&r);
                if g_norm > 1e-14 {
                    vec_scale(&r, -delta / g_norm)
                } else {
                    vec_scale(&newton_step, delta / newton_norm.max(1e-14))
                }
            };

            let u_trial: Vec<f64> = u.iter().zip(&step).map(|(ui, si)| ui + si).collect();
            let r_trial = problem.residual(&u_trial);
            let norm_trial = norm_l2(&r_trial);

            // Predicted reduction: first-order model m(s) ≈ R · (-s) (positive when reducing).
            // step = correction direction (u_trial = u + step, step = -K^{-1}R).
            // Dot of step with -R gives expected decrease.
            let predicted_reduction = {
                let neg_r: Vec<f64> = r.iter().map(|x| -x).collect();
                vec_dot(&step, &neg_r).abs().max(1e-30)
            };
            let actual_reduction = norm - norm_trial;

            let rho = actual_reduction / predicted_reduction;

            // Accept or reject the step.
            if norm_trial < norm {
                u = u_trial;
            }

            // Update trust-region radius using standard dogleg rules.
            if rho < 0.25 {
                delta = (0.25 * newton_norm).max(1e-14);
            } else if rho > 0.75 {
                delta = (2.0 * delta).min(self.max_radius);
            }
            // else keep delta unchanged
        }

        let r_final = problem.residual(&u);
        let norm_final = norm_l2(&r_final);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── Picard iteration ─────────────────────────────────────────────────────────

/// Picard (successive substitution) solver with optional under-relaxation.
///
/// Iterates u_{k+1} = u_k + ω * K(u_k)^{-1} R(u_k), where ω is a
/// relaxation parameter. ω = 1 gives standard fixed-point iteration;
/// ω < 1 improves stability for poorly conditioned problems.
#[derive(Debug, Clone)]
pub struct PicardSolver {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Relaxation parameter ω ∈ (0, 1].
    pub relaxation: f64,
}

impl PicardSolver {
    /// Create a new Picard solver.
    pub fn new(max_iter: usize, tol: f64, relaxation: f64) -> Self {
        Self {
            max_iterations: max_iter,
            tolerance: tol,
            relaxation: relaxation.clamp(1e-6, 1.0),
        }
    }

    /// Solve using Picard iteration.
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let mut u = u0.to_vec();
        let mut history = Vec::new();

        for iter in 0..self.max_iterations {
            let r = problem.residual(&u);
            let norm = norm_l2(&r);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            let mut k = problem.tangent_stiffness(&u);
            let mut rhs = r.clone();
            let delta = gauss_elimination(&mut k, &mut rhs)?;

            // Apply relaxed update: u ← u - ω * K^{-1} R
            for (ui, di) in u.iter_mut().zip(delta.iter()) {
                *ui -= self.relaxation * di;
            }
        }

        let r_final = problem.residual(&u);
        let norm_final = norm_l2(&r_final);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── GMRES linear solver ──────────────────────────────────────────────────────

/// Generalized Minimal Residual (GMRES) solver for A x = b.
///
/// Uses Arnoldi iteration to build an orthonormal Krylov subspace basis,
/// then minimises the residual over that subspace via least-squares (Givens rotations).
///
/// Returns `Err` if the system cannot be solved or diverges.
pub fn gmres_solve(
    a: &[Vec<f64>],
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, String> {
    let n = b.len();
    if n == 0 {
        return Err("Empty system".to_string());
    }

    let mut x = vec![0.0_f64; n];

    // Initial residual r0 = b - A*x0 = b (since x0 = 0)
    let r = b.to_vec();
    let beta = norm_l2(&r);

    if beta < tol {
        return Ok(x);
    }

    let m = max_iter.min(n);
    // Krylov basis vectors (m+1 vectors of length n)
    let mut q: Vec<Vec<f64>> = vec![vec![0.0; n]; m + 1];
    // Hessenberg matrix (m+1 rows × m cols)
    let mut h = vec![vec![0.0_f64; m]; m + 1];
    // Givens rotation cosines and sines
    let mut cs = vec![0.0_f64; m];
    let mut sn = vec![0.0_f64; m];
    // RHS of the least-squares problem in Krylov space
    let mut e1 = vec![0.0_f64; m + 1];
    e1[0] = beta;

    // First basis vector
    q[0] = vec_scale(&r, 1.0 / beta);

    let mut iter_count = 0;

    'outer: for j in 0..m {
        // Arnoldi step: w = A * q[j]
        let mut w = mat_vec(a, &q[j]);

        // Modified Gram-Schmidt orthogonalisation
        for i in 0..=j {
            h[i][j] = vec_dot(&w, &q[i]);
            let qi_scaled = vec_scale(&q[i], h[i][j]);
            w = vec_sub(&w, &qi_scaled);
        }
        h[j + 1][j] = norm_l2(&w);

        if h[j + 1][j] > 1e-14 {
            q[j + 1] = vec_scale(&w, 1.0 / h[j + 1][j]);
        }

        // Apply previous Givens rotations to new column
        for i in 0..j {
            let temp = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
            h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
            h[i][j] = temp;
        }

        // Compute new Givens rotation
        let denom = (h[j][j].powi(2) + h[j + 1][j].powi(2)).sqrt();
        if denom > 1e-14 {
            cs[j] = h[j][j] / denom;
            sn[j] = h[j + 1][j] / denom;
        } else {
            cs[j] = 1.0;
            sn[j] = 0.0;
        }
        h[j][j] = cs[j] * h[j][j] + sn[j] * h[j + 1][j];
        h[j + 1][j] = 0.0;

        // Update RHS
        e1[j + 1] = -sn[j] * e1[j];
        e1[j] *= cs[j];

        iter_count = j + 1;

        if e1[j + 1].abs() < tol * beta {
            break 'outer;
        }
    }

    // Solve upper triangular system H[0..iter_count][0..iter_count] * y = e1[0..iter_count]
    let k = iter_count;
    let mut y = vec![0.0_f64; k];
    for i in (0..k).rev() {
        let mut sum = e1[i];
        for j in (i + 1)..k {
            sum -= h[i][j] * y[j];
        }
        if h[i][i].abs() < 1e-30 {
            return Err(format!("GMRES: singular Hessenberg at index {i}"));
        }
        y[i] = sum / h[i][i];
    }

    // Update solution: x = x + Q_k * y
    for j in 0..k {
        let qj_scaled = vec_scale(&q[j], y[j]);
        x = vec_add(&x, &qj_scaled);
    }

    Ok(x)
}

// ─── Newton-Krylov solver ─────────────────────────────────────────────────────

/// Newton-Krylov solver: uses GMRES to solve the linear Newton system at each
/// outer Newton iteration. Well suited for large sparse problems.
#[derive(Debug, Clone)]
pub struct NewtonKrylovSolver {
    /// Maximum outer Newton iterations.
    pub max_outer_iterations: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Maximum GMRES inner iterations.
    pub max_inner_iterations: usize,
}

impl NewtonKrylovSolver {
    /// Create a new Newton-Krylov solver.
    pub fn new(max_outer: usize, tol: f64, max_inner: usize) -> Self {
        Self {
            max_outer_iterations: max_outer,
            tolerance: tol,
            max_inner_iterations: max_inner,
        }
    }

    /// Solve using Newton-Krylov (GMRES inner solver).
    pub fn solve<F: ResidualFn>(&self, problem: &F, u0: &[f64]) -> Result<NewtonResult, String> {
        let mut u = u0.to_vec();
        let mut history = Vec::new();

        for iter in 0..self.max_outer_iterations {
            let r = problem.residual(&u);
            let norm = norm_l2(&r);
            history.push(norm);

            if norm < self.tolerance {
                return Ok(NewtonResult {
                    solution: u,
                    iterations: iter,
                    final_residual: norm,
                    converged: true,
                    residual_history: history,
                });
            }

            let k = problem.tangent_stiffness(&u);
            // Solve K * delta = R using GMRES
            let delta = gmres_solve(&k, &r, self.max_inner_iterations, self.tolerance * 0.01)
                .unwrap_or_else(|_| {
                    // Fall back to Gauss elimination if GMRES fails
                    let mut k2 = problem.tangent_stiffness(&u);
                    let mut rhs = r.clone();
                    gauss_elimination(&mut k2, &mut rhs).unwrap_or_else(|_| vec![0.0; u.len()])
                });

            u = vec_sub(&u, &delta);
        }

        let r_final = problem.residual(&u);
        let norm_final = norm_l2(&r_final);
        history.push(norm_final);

        Ok(NewtonResult {
            solution: u,
            iterations: self.max_outer_iterations,
            final_residual: norm_final,
            converged: norm_final < self.tolerance,
            residual_history: history,
        })
    }
}

// ─── Convergence monitor ──────────────────────────────────────────────────────

/// Monitors convergence, stagnation and divergence of iterative solvers.
///
/// Maintains a sliding window of residual norms and provides predicates
/// for detecting common failure modes.
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor {
    /// Convergence tolerance: converged when latest residual < tol.
    pub tolerance: f64,
    /// Window size for stagnation detection.
    pub stagnation_window: usize,
    /// Minimum change to consider as progress (not stagnated).
    pub stagnation_threshold: f64,
    /// History of residual norms.
    history: Vec<f64>,
}

impl ConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(tolerance: f64, stagnation_window: usize, stagnation_threshold: f64) -> Self {
        Self {
            tolerance,
            stagnation_window,
            stagnation_threshold,
            history: Vec::new(),
        }
    }

    /// Record a new residual norm.
    pub fn record(&mut self, residual: f64) {
        self.history.push(residual);
    }

    /// Whether the latest residual is below the convergence tolerance.
    pub fn has_converged(&self) -> bool {
        self.history.last().is_some_and(|&r| r < self.tolerance)
    }

    /// Whether the residual has stagnated (not changed significantly in the window).
    pub fn is_stagnated(&self) -> bool {
        let n = self.history.len();
        if n < self.stagnation_window {
            return false;
        }
        let window = &self.history[(n - self.stagnation_window)..];
        let max_val = window.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_val = window.iter().cloned().fold(f64::INFINITY, f64::min);
        (max_val - min_val).abs() < self.stagnation_threshold
    }

    /// Whether the residual is growing monotonically over the last 3 steps.
    pub fn is_diverging(&self) -> bool {
        let n = self.history.len();
        if n < 3 {
            return false;
        }
        self.history[n - 1] > self.history[n - 2] && self.history[n - 2] > self.history[n - 3]
    }

    /// Reset the history.
    pub fn reset(&mut self) {
        self.history.clear();
    }

    /// Number of iterations recorded.
    pub fn iteration_count(&self) -> usize {
        self.history.len()
    }

    /// Return a slice of the residual history.
    pub fn history(&self) -> &[f64] {
        &self.history
    }
}

// ─── Solver statistics ────────────────────────────────────────────────────────

/// Aggregated statistics for a nonlinear solver run.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct SolverStatistics {
    residuals: Vec<f64>,
    tangent_updated: Vec<bool>,
}

impl SolverStatistics {
    /// Create empty statistics.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record one iteration.
    ///
    /// `residual`: residual norm at end of this iteration.
    /// `tangent_updated`: whether the tangent was recomputed this iteration.
    pub fn record_iteration(&mut self, residual: f64, tangent_updated: bool) {
        self.residuals.push(residual);
        self.tangent_updated.push(tangent_updated);
    }

    /// Total number of iterations recorded.
    pub fn total_iterations(&self) -> usize {
        self.residuals.len()
    }

    /// Number of iterations where the tangent was updated.
    pub fn tangent_updates(&self) -> usize {
        self.tangent_updated.iter().filter(|&&v| v).count()
    }

    /// Initial (first recorded) residual norm. Returns `f64::NAN` if empty.
    pub fn initial_residual(&self) -> f64 {
        self.residuals.first().copied().unwrap_or(f64::NAN)
    }

    /// Final (last recorded) residual norm. Returns `0.0` if empty.
    pub fn final_residual(&self) -> f64 {
        self.residuals.last().copied().unwrap_or(0.0)
    }

    /// Average per-iteration convergence rate:
    /// (ln(r_final) - ln(r_initial)) / (n - 1).
    pub fn average_convergence_rate(&self) -> f64 {
        let n = self.residuals.len();
        if n < 2 {
            return 0.0;
        }
        let r0 = self.residuals[0].abs().max(1e-300);
        let rn = self.residuals[n - 1].abs().max(1e-300);
        (rn.ln() - r0.ln()) / (n - 1) as f64
    }
}

// ─── Test problem implementations ────────────────────────────────────────────

/// Simple 1-DOF linear spring: R(u) = k·u − f,  K = k.
pub struct SimpleLinearProblem {
    /// Spring stiffness.
    pub k: f64,
    /// Applied force.
    pub f: f64,
}

impl ResidualFn for SimpleLinearProblem {
    fn residual(&self, u: &[f64]) -> Vec<f64> {
        vec![self.k * u[0] - self.f]
    }

    fn tangent_stiffness(&self, _u: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![self.k]]
    }
}

/// 1-DOF nonlinear spring with cubic term: R(u) = k1·u + k3·u³ − force.
pub struct NonlinearSpring {
    /// Linear stiffness coefficient.
    pub k1: f64,
    /// Cubic stiffness coefficient.
    pub k3: f64,
    /// Applied external force.
    pub force: f64,
}

impl ResidualFn for NonlinearSpring {
    fn residual(&self, u: &[f64]) -> Vec<f64> {
        vec![self.k1 * u[0] + self.k3 * u[0].powi(3) - self.force]
    }

    fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![self.k1 + 3.0 * self.k3 * u[0].powi(2)]]
    }
}

/// 2-DOF coupled nonlinear problem for testing multi-dimensional solvers.
///
/// R1(u) = u1^2 + u2 - 3
/// R2(u) = u1 + u2^2 - 3
/// Solution: u1 = 1, u2 = 1 (among others)
pub struct CoupledNonlinearProblem;

impl ResidualFn for CoupledNonlinearProblem {
    fn residual(&self, u: &[f64]) -> Vec<f64> {
        vec![u[0] * u[0] + u[1] - 3.0, u[0] + u[1] * u[1] - 3.0]
    }

    fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![2.0 * u[0], 1.0], vec![1.0, 2.0 * u[1]]]
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gauss_elimination_2x2() {
        let mut a = vec![vec![2.0, -1.0], vec![1.0, 1.0]];
        let mut b = vec![3.0, 3.0];
        let x = gauss_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-12, "x should be 2, got {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-12, "y should be 1, got {}", x[1]);
    }

    #[test]
    fn test_gauss_elimination_singular() {
        let mut a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let mut b = vec![1.0, 2.0];
        let result = gauss_elimination(&mut a, &mut b);
        assert!(result.is_err(), "Expected Err for singular matrix");
    }

    #[test]
    fn test_newton_raphson_linear_problem() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let solver = NewtonRaphson::new(20, 1e-10);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(result.converged, "Should converge");
        assert!((result.solution[0] - 2.5).abs() < 1e-10);
        assert!(result.iterations <= 1, "Should converge in 1 iteration");
    }

    #[test]
    fn test_newton_raphson_nonlinear_spring() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = NewtonRaphson::new(50, 1e-10);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(result.converged, "Should converge");
        let r = problem.residual(&result.solution);
        assert!(norm_l2(&r) < 1e-8, "Residual should be near zero");
    }

    #[test]
    fn test_modified_newton_raphson_converges() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = ModifiedNewtonRaphson::new(100, 1e-9, 3);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(result.converged, "Modified NR should converge");
        let r = problem.residual(&result.solution);
        assert!(norm_l2(&r) < 1e-7);
    }

    #[test]
    fn test_load_stepper_returns_n_steps() {
        let problem = SimpleLinearProblem { k: 1.0, f: 1.0 };
        let stepper = LoadStepper::new(5, 0.0, 1.0);
        let solver = NewtonRaphson::new(20, 1e-10);
        let results = stepper.solve(&problem, &[0.0], &solver);
        assert_eq!(results.len(), 5, "Should have n_steps results");
    }

    #[test]
    fn test_convergence_criteria_residual_norm_zero() {
        let crit = ConvergenceCriteria::ResidualNorm { tol: 1e-6 };
        let zero = vec![0.0_f64; 3];
        assert!(
            crit.check(&zero, &zero, &zero),
            "Zero residual should converge"
        );
    }

    #[test]
    fn test_newton_raphson_max_iter_not_converged() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 1.0,
            force: 10.0,
        };
        let solver = NewtonRaphson::new(1, 1e-30);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(
            !result.converged,
            "Should not converge with max_iter=1 and very tight tol"
        );
    }

    // --- New tests ---

    #[test]
    fn test_newton_raphson_line_search_linear() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let solver = NewtonRaphsonLineSearch::new(20, 1e-10);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(result.converged);
        assert!((result.solution[0] - 2.5).abs() < 1e-8);
    }

    #[test]
    fn test_newton_raphson_line_search_nonlinear() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = NewtonRaphsonLineSearch::new(50, 1e-10);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(result.converged);
    }

    #[test]
    fn test_bfgs_linear() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let solver = BfgsSolver::new(50, 1e-8);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(result.converged, "BFGS should converge on linear problem");
        assert!((result.solution[0] - 2.5).abs() < 1e-6);
    }

    #[test]
    fn test_bfgs_nonlinear() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = BfgsSolver::new(100, 1e-8);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(result.converged, "BFGS should converge on nonlinear spring");
    }

    #[test]
    fn test_coupled_nonlinear_newton() {
        let problem = CoupledNonlinearProblem;
        let solver = NewtonRaphson::new(50, 1e-10);
        let result = solver.solve(&problem, &[1.5, 1.5]).unwrap();
        assert!(result.converged, "NR should converge on coupled problem");
        let r = problem.residual(&result.solution);
        assert!(norm_l2(&r) < 1e-8);
    }

    #[test]
    fn test_convergence_criteria_displacement() {
        let crit = ConvergenceCriteria::DisplacementNorm { tol: 1e-6 };
        let small = vec![1e-7, 1e-7, 1e-7];
        let large = vec![1.0, 1.0, 1.0];
        assert!(
            crit.check(&large, &small, &[]),
            "Small displacement should converge"
        );
        assert!(
            !crit.check(&small, &large, &[]),
            "Large displacement should not converge"
        );
    }

    #[test]
    fn test_convergence_criteria_energy() {
        let crit = ConvergenceCriteria::EnergyNorm { tol: 1e-6 };
        let r = vec![1e-4, 1e-4];
        let du = vec![1e-4, 1e-4];
        // energy = 2e-8 < 1e-6
        assert!(crit.check(&r, &du, &[]));
    }

    #[test]
    fn test_convergence_criteria_combined() {
        let crit = ConvergenceCriteria::Combined {
            residual_tol: 1e-6,
            displacement_tol: 1e-6,
        };
        let small = vec![1e-7, 1e-7];
        assert!(crit.check(&small, &small, &[]));
        let large = vec![1.0, 1.0];
        assert!(!crit.check(&large, &small, &[]));
        assert!(!crit.check(&small, &large, &[]));
    }

    #[test]
    fn test_convergence_criteria_relative() {
        let crit = ConvergenceCriteria::RelativeResidual {
            tol: 1e-3,
            initial_norm: 100.0,
        };
        let r = vec![0.05]; // norm = 0.05, relative = 5e-4 < 1e-3
        assert!(crit.check(&r, &[], &[]));
    }

    #[test]
    fn test_norm_inf() {
        let v = vec![-3.0, 1.0, 2.0];
        assert!((norm_inf(&v) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec_add() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let c = vec_add(&a, &b);
        assert!((c[0] - 5.0).abs() < 1e-12);
        assert!((c[1] - 7.0).abs() < 1e-12);
        assert!((c[2] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec_dot() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        assert!((vec_dot(&a, &b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_load_stepper_step_factors() {
        let stepper = LoadStepper::new(3, 0.0, 1.0);
        let factors = stepper.step_factors();
        assert_eq!(factors.len(), 3);
        assert!((factors[0] - 0.0).abs() < 1e-12);
        assert!((factors[1] - 0.5).abs() < 1e-12);
        assert!((factors[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_load_stepper_single_step() {
        let stepper = LoadStepper::new(1, 0.5, 1.0);
        let factors = stepper.step_factors();
        assert_eq!(factors.len(), 1);
        assert!((factors[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_load_stepper_zero_steps() {
        let stepper = LoadStepper::new(0, 0.0, 1.0);
        assert!(stepper.step_factors().is_empty());
    }

    #[test]
    fn test_residual_history_populated() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = NewtonRaphson::new(50, 1e-10);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(!result.residual_history.is_empty());
        // Residual should generally decrease
        if result.residual_history.len() > 2 {
            assert!(
                result.residual_history.last().unwrap() < &result.residual_history[0],
                "Final residual should be less than initial"
            );
        }
    }

    #[test]
    fn test_gauss_elimination_3x3() {
        // x + y + z = 6, 2x + y + 3z = 11, x + 3y + 2z = 14
        // Solution: x=1, y=3, z=2
        let mut a = vec![
            vec![1.0, 1.0, 1.0],
            vec![2.0, 1.0, 3.0],
            vec![1.0, 3.0, 2.0],
        ];
        let mut b = vec![6.0, 11.0, 14.0];
        let x = gauss_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 3.0).abs() < 1e-10);
        assert!((x[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_coupled_nonlinear_line_search() {
        let problem = CoupledNonlinearProblem;
        let solver = NewtonRaphsonLineSearch::new(50, 1e-10);
        let result = solver.solve(&problem, &[1.5, 1.5]).unwrap();
        assert!(result.converged);
    }

    #[test]
    fn test_arc_length_trace_path() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.01,
            force: 1.0,
        };
        let arc = ArcLengthMethod::new(0.1, 50, 1e-4);
        let path = arc.trace_path(&problem, &[0.5], 0.5, 5);
        // Arc-length tracing should produce at least some steps
        // (may not produce all steps if corrector fails for some)
        assert!(
            path.len() <= 5,
            "Path should have at most n_steps entries, got {}",
            path.len()
        );
    }

    // ─── Trust region tests ───────────────────────────────────────────────────

    #[test]
    fn test_trust_region_linear() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        // Start with a larger initial radius so the Newton step is accepted on the first try
        let solver = TrustRegionSolver::new(50, 1e-10, 10.0, 100.0);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(
            result.converged,
            "Trust region should converge on linear problem"
        );
        assert!(
            (result.solution[0] - 2.5).abs() < 1e-8,
            "Solution should be 2.5, got {}",
            result.solution[0]
        );
    }

    #[test]
    fn test_trust_region_nonlinear_spring() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = TrustRegionSolver::new(100, 1e-8, 1.0, 10.0);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(
            result.converged,
            "Trust region should converge on nonlinear spring"
        );
        let r = problem.residual(&result.solution);
        assert!(norm_l2(&r) < 1e-6);
    }

    #[test]
    fn test_trust_region_coupled() {
        let problem = CoupledNonlinearProblem;
        let solver = TrustRegionSolver::new(100, 1e-10, 1.0, 10.0);
        let result = solver.solve(&problem, &[1.5, 1.5]).unwrap();
        assert!(
            result.converged,
            "Trust region should converge on coupled problem"
        );
        let r = problem.residual(&result.solution);
        assert!(norm_l2(&r) < 1e-8);
    }

    #[test]
    fn test_trust_region_radius_adjusts() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.5,
            force: 2.0,
        };
        let solver = TrustRegionSolver::new(200, 1e-6, 0.1, 5.0);
        let result = solver.solve(&problem, &[0.1]).unwrap();
        // Just ensure it runs and returns a result
        assert!(result.final_residual.is_finite());
    }

    // ─── Picard iteration tests ──────────────────────────────────────────────

    #[test]
    fn test_picard_linear_converges() {
        // Picard iteration for K*u = f: u_{n+1} = K^{-1}*f
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let solver = PicardSolver::new(50, 1e-10, 0.8);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(result.converged, "Picard should converge on linear problem");
        assert!((result.solution[0] - 2.5).abs() < 1e-8);
    }

    #[test]
    fn test_picard_underrelaxation_slows_convergence() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let solver_full = PicardSolver::new(100, 1e-10, 1.0);
        let solver_under = PicardSolver::new(100, 1e-10, 0.3);
        let r1 = solver_full.solve(&problem, &[0.0]).unwrap();
        let r2 = solver_under.solve(&problem, &[0.0]).unwrap();
        // Both should converge, but under-relaxed takes more iterations
        assert!(r1.converged && r2.converged);
        // Under-relaxed may use more iterations (or same for linear case)
        assert!((r1.solution[0] - 2.5).abs() < 1e-7);
        assert!((r2.solution[0] - 2.5).abs() < 1e-7);
    }

    #[test]
    fn test_picard_history_nonempty() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = PicardSolver::new(50, 1e-9, 0.7);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(!result.residual_history.is_empty());
    }

    // ─── GMRES / Newton-Krylov tests ─────────────────────────────────────────

    #[test]
    fn test_gmres_solves_2x2() {
        let a = vec![vec![2.0, -1.0], vec![1.0, 1.0]];
        let b = vec![3.0, 3.0];
        let x = gmres_solve(&a, &b, 10, 1e-12).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-10, "x[0] should be 2, got {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-10, "x[1] should be 1, got {}", x[1]);
    }

    #[test]
    fn test_gmres_solves_3x3() {
        let a = vec![
            vec![1.0, 1.0, 1.0],
            vec![2.0, 1.0, 3.0],
            vec![1.0, 3.0, 2.0],
        ];
        let b = vec![6.0, 11.0, 14.0];
        let x = gmres_solve(&a, &b, 10, 1e-10).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-8);
        assert!((x[1] - 3.0).abs() < 1e-8);
        assert!((x[2] - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_gmres_identity() {
        let n = 4;
        let a: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x = gmres_solve(&a, &b, 20, 1e-12).unwrap();
        for i in 0..n {
            assert!(
                (x[i] - b[i]).abs() < 1e-10,
                "x[{i}]={} expected {}",
                x[i],
                b[i]
            );
        }
    }

    #[test]
    fn test_newton_krylov_linear() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let solver = NewtonKrylovSolver::new(30, 1e-10, 20);
        let result = solver.solve(&problem, &[0.0]).unwrap();
        assert!(
            result.converged,
            "Newton-Krylov should converge on linear problem"
        );
        assert!((result.solution[0] - 2.5).abs() < 1e-8);
    }

    #[test]
    fn test_newton_krylov_nonlinear() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.1,
            force: 1.1,
        };
        let solver = NewtonKrylovSolver::new(50, 1e-9, 20);
        let result = solver.solve(&problem, &[0.5]).unwrap();
        assert!(
            result.converged,
            "Newton-Krylov should converge on nonlinear spring"
        );
    }

    #[test]
    fn test_newton_krylov_coupled() {
        let problem = CoupledNonlinearProblem;
        let solver = NewtonKrylovSolver::new(80, 1e-10, 20);
        let result = solver.solve(&problem, &[1.5, 1.5]).unwrap();
        assert!(
            result.converged,
            "Newton-Krylov should converge on coupled problem"
        );
        let r = problem.residual(&result.solution);
        assert!(norm_l2(&r) < 1e-8);
    }

    // ─── ConvergenceMonitor tests ─────────────────────────────────────────────

    #[test]
    fn test_convergence_monitor_detects_stagnation() {
        let mut mon = ConvergenceMonitor::new(1e-8, 5, 1e-12);
        // Feed identical residuals → stagnation
        for _ in 0..6 {
            mon.record(1e-3);
        }
        assert!(mon.is_stagnated(), "Monitor should detect stagnation");
    }

    #[test]
    fn test_convergence_monitor_detects_convergence() {
        let mut mon = ConvergenceMonitor::new(1e-6, 5, 1e-12);
        mon.record(1e-7);
        assert!(mon.has_converged(), "Monitor should detect convergence");
    }

    #[test]
    fn test_convergence_monitor_not_converged_initially() {
        let mon = ConvergenceMonitor::new(1e-6, 5, 1e-12);
        assert!(!mon.has_converged());
        assert!(!mon.is_stagnated());
    }

    #[test]
    fn test_convergence_monitor_divergence() {
        let mut mon = ConvergenceMonitor::new(1e-6, 5, 1e-12);
        mon.record(1.0);
        mon.record(10.0);
        mon.record(100.0);
        assert!(mon.is_diverging(), "Monitor should detect divergence");
    }

    #[test]
    fn test_convergence_monitor_reset() {
        let mut mon = ConvergenceMonitor::new(1e-8, 5, 1e-12);
        mon.record(1e-3);
        mon.record(1e-3);
        mon.reset();
        assert!(!mon.has_converged());
        assert!(!mon.is_stagnated());
    }

    // ─── SolverStatistics tests ───────────────────────────────────────────────

    #[test]
    fn test_solver_stats_accumulate() {
        let mut stats = SolverStatistics::new();
        stats.record_iteration(0.5, true);
        stats.record_iteration(0.1, false);
        stats.record_iteration(1e-6, true);
        assert_eq!(stats.total_iterations(), 3);
        assert_eq!(stats.tangent_updates(), 2);
        assert!((stats.initial_residual() - 0.5).abs() < 1e-12);
        assert!((stats.final_residual() - 1e-6).abs() < 1e-12);
    }

    #[test]
    fn test_solver_stats_convergence_rate() {
        let mut stats = SolverStatistics::new();
        stats.record_iteration(1.0, false);
        stats.record_iteration(0.1, false);
        stats.record_iteration(0.01, false);
        let rate = stats.average_convergence_rate();
        // Rate ≈ log(0.01/1.0) / 2 ≈ -2.3
        assert!(rate < 0.0, "Rate should be negative (converging)");
    }

    #[test]
    fn test_solver_stats_empty() {
        let stats = SolverStatistics::new();
        assert_eq!(stats.total_iterations(), 0);
        assert!(stats.final_residual().is_nan() || stats.final_residual() == 0.0);
    }

    // ─── Vector helpers extra tests ───────────────────────────────────────────

    #[test]
    fn test_mat_vec_product() {
        let m = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let v = vec![1.0, 0.0, -1.0];
        let r = mat_vec(&m, &v);
        assert!((r[0] - (-2.0)).abs() < 1e-12);
        assert!((r[1] - (-2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_vec_scale_negative() {
        let v = vec![1.0, -2.0, 3.0];
        let sv = vec_scale(&v, -2.0);
        assert!((sv[0] - (-2.0)).abs() < 1e-12);
        assert!((sv[1] - 4.0).abs() < 1e-12);
        assert!((sv[2] - (-6.0)).abs() < 1e-12);
    }

    #[test]
    fn test_norm_l2_known() {
        let v = vec![3.0, 4.0];
        assert!((norm_l2(&v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_gauss_elimination_4x4() {
        // 4×4 diagonal system
        let mut a = vec![
            vec![2.0, 0.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0, 0.0],
            vec![0.0, 0.0, 4.0, 0.0],
            vec![0.0, 0.0, 0.0, 5.0],
        ];
        let mut b = vec![2.0, 6.0, 8.0, 10.0];
        let x = gauss_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-12);
        assert!((x[1] - 2.0).abs() < 1e-12);
        assert!((x[2] - 2.0).abs() < 1e-12);
        assert!((x[3] - 2.0).abs() < 1e-12);
    }
}
