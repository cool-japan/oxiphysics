//! Solver for the condensed linear-quadratic MPC problem
//!
//! ```text
//! min_u  J(u) = Σ_{k=0}^{N-1} (x_kᵀ Q x_k + u_kᵀ R u_k) + x_Nᵀ Qf x_N
//! s.t.   x_{k+1} = A x_k + B u_k,   x_0 given
//!        u_lb <= u_k <= u_ub,       x_lb <= x_k <= x_ub  (k = 1..N)
//! ```
//!
//! over the input sequence `u`.  The gradient of `J` is obtained exactly with
//! one backward (adjoint) pass,
//!
//! ```text
//! p_N = 2 Qf x_N,   p_k = 2 Q x_k + Aᵀ p_{k+1},   ∂J/∂u_k = 2 R u_k + Bᵀ p_{k+1},
//! ```
//!
//! The input box is handled by an accelerated projected-gradient method (FISTA
//! with gradient-scheme restart and backtracking, started from the step `1/L`
//! where `L` is a power-iteration estimate of the largest curvature) that stops
//! once the gradient-mapping norm is below a scale-invariant tolerance, see
//! [`StoppingRule`].
//!
//! State bounds are enforced by an augmented-Lagrangian outer loop (method of
//! multipliers) around that solver.  The multiplier-weighted violations enter
//! the adjoint recursion as extra state forces, so the exact solution of the
//! constrained problem is reached with a finite penalty.  The inner solves are
//! inexact, with an accuracy that tracks the violation left by the previous
//! outer iteration; the penalty is raised only while the bounds are violated
//! and the violation stalls; the loop stops when the iterate, paired with the
//! updated multipliers, is a KKT point of the original problem (stationarity
//! from the inner solve, primal violation and complementarity gap in state
//! units).  An infeasible problem ends with the least-violating input sequence
//! and `feasible = false`.

use super::types::{MpcProblem, MpcSolution};

/// Iteration limit of one box-constrained inner solve.
const MAX_INNER_ITERATIONS: usize = 5000;
/// Iteration limit of the augmented-Lagrangian outer loop.
const MAX_OUTER_ITERATIONS: usize = 80;
/// Power iterations used to estimate the largest curvature.
const POWER_ITERATIONS: usize = 100;
/// Safety factor on the power-iteration estimate (which approaches `λ_max` from below).
const LIPSCHITZ_SAFETY: f64 = 1.05;
/// Relative accuracy of the stopping rule on the projected-gradient
/// (gradient-mapping) norm, see [`StoppingRule`].
const GRADIENT_TOL: f64 = 1e-10;
/// Loosest relative accuracy of an inner solve (the first outer iteration).
const INEXACT_INNER_TOL: f64 = 1e-3;
/// State-bound violation and complementarity gap below which the outer loop
/// has converged (relative to the size of the bounds, at least `1`).
const FEASIBILITY_TOL: f64 = 1e-9;
/// Violation up to which [`MpcSolution::feasible`] still reports `true`.
const REPORT_FEASIBILITY_TOL: f64 = 1e-6;
/// Factor by which the penalty is raised when the violation stalls.
const PENALTY_GROWTH: f64 = 10.0;
/// The penalty is raised when the violation shrinks by less than this factor
/// per outer iteration.
const PENALTY_STALL_RATIO: f64 = 0.5;
/// The penalty never exceeds this multiple of its initial value.  The cost of
/// an inner solve grows like the square root of the penalty, so a larger one
/// would only pay off for problems that are infeasible anyway.
const PENALTY_MAX_FACTOR: f64 = 1e2;

/// `Aᵀ` of the row-major `rows x cols` matrix.
fn transpose(a: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut t = vec![0.0; a.len()];
    for r in 0..rows {
        for c in 0..cols {
            t[c * rows + r] = a[r * cols + c];
        }
    }
    t
}

/// `y = A x` for the row-major `rows x cols` matrix `A`, accumulated into `y`.
fn mat_vec_add(a: &[f64], rows: usize, cols: usize, x: &[f64], y: &mut [f64]) {
    for r in 0..rows {
        let row = &a[r * cols..(r + 1) * cols];
        y[r] += row.iter().zip(x.iter()).map(|(p, q)| p * q).sum::<f64>();
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm_inf(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

/// `M + Mᵀ` of a row-major `n x n` matrix (twice its symmetric part).
fn doubled_symmetric_part(a: &[f64], n: usize) -> Vec<f64> {
    let mut s = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            s[i * n + j] = a[i * n + j] + a[j * n + i];
        }
    }
    s
}

/// Bound `i` of `bounds`, or `default` when the slice is shorter or the entry is NaN.
fn bound_or(bounds: &[f64], i: usize, default: f64) -> f64 {
    match bounds.get(i) {
        Some(v) if !v.is_nan() => *v,
        _ => default,
    }
}

/// Multipliers and penalty parameter of the augmented Lagrangian for the
/// state bounds (one entry per stage `0..=N` and state; stage 0 is unused
/// because the initial state is not a decision variable).
struct Penalty {
    rho: f64,
    upper: Vec<f64>,
    lower: Vec<f64>,
}

/// The condensed MPC problem in the form the solver works on.
struct Mpc<'a> {
    problem: &'a MpcProblem,
    n: usize,
    m: usize,
    horizon: usize,
    /// `Q + Qᵀ`, `Qf + Qfᵀ`, `R + Rᵀ` (the gradient of `xᵀ M x` is `(M + Mᵀ) x`).
    q2: Vec<f64>,
    qf2: Vec<f64>,
    r2: Vec<f64>,
    a_t: Vec<f64>,
    b_t: Vec<f64>,
    u_lb: Vec<f64>,
    u_ub: Vec<f64>,
    x_lb: Vec<f64>,
    x_ub: Vec<f64>,
    /// Whether state `i` has a finite lower or upper bound.
    constrained: Vec<bool>,
}

impl<'a> Mpc<'a> {
    /// `None` when the matrix sizes do not match the declared dimensions.
    fn new(problem: &'a MpcProblem) -> Option<Self> {
        let (n, m) = (problem.n, problem.m);
        let consistent = problem.a.len() == n * n
            && problem.b.len() == n * m
            && problem.q.len() == n * n
            && problem.qf.len() == n * n
            && problem.r.len() == m * m;
        if !consistent {
            return None;
        }
        let u_lb: Vec<f64> = (0..m)
            .map(|j| bound_or(&problem.u_lb, j, f64::NEG_INFINITY))
            .collect();
        let u_ub: Vec<f64> = (0..m)
            .map(|j| bound_or(&problem.u_ub, j, f64::INFINITY))
            .collect();
        let x_lb: Vec<f64> = (0..n)
            .map(|i| bound_or(&problem.x_lb, i, f64::NEG_INFINITY))
            .collect();
        let x_ub: Vec<f64> = (0..n)
            .map(|i| bound_or(&problem.x_ub, i, f64::INFINITY))
            .collect();
        let constrained = (0..n)
            .map(|i| x_lb[i].is_finite() || x_ub[i].is_finite())
            .collect();
        Some(Self {
            problem,
            n,
            m,
            horizon: problem.horizon,
            q2: doubled_symmetric_part(&problem.q, n),
            qf2: doubled_symmetric_part(&problem.qf, n),
            r2: doubled_symmetric_part(&problem.r, m),
            a_t: transpose(&problem.a, n, n),
            b_t: transpose(&problem.b, n, m),
            u_lb,
            u_ub,
            x_lb,
            x_ub,
            constrained,
        })
    }

    /// State trajectory `(N + 1) x n` from `x0` under the inputs `u`.
    fn simulate(&self, x0: &[f64], u: &[f64]) -> Vec<f64> {
        let (n, m) = (self.n, self.m);
        let mut traj = vec![0.0; (self.horizon + 1) * n];
        traj[..n].copy_from_slice(x0);
        for k in 0..self.horizon {
            let (head, tail) = traj.split_at_mut((k + 1) * n);
            let x_k = &head[k * n..];
            let x_next = &mut tail[..n];
            mat_vec_add(&self.problem.a, n, n, x_k, x_next);
            mat_vec_add(&self.problem.b, n, m, &u[k * m..(k + 1) * m], x_next);
        }
        traj
    }

    /// `J(u)` for the trajectory `traj` generated by `u`.
    fn cost(&self, traj: &[f64], u: &[f64]) -> f64 {
        let (n, m) = (self.n, self.m);
        let quad = |matrix: &[f64], v: &[f64], dim: usize| -> f64 {
            let mut mv = vec![0.0; dim];
            mat_vec_add(matrix, dim, dim, v, &mut mv);
            dot(v, &mv)
        };
        let mut total = 0.0;
        for k in 0..self.horizon {
            total += quad(&self.problem.q, &traj[k * n..(k + 1) * n], n);
            total += quad(&self.problem.r, &u[k * m..(k + 1) * m], m);
        }
        total + quad(&self.problem.qf, &traj[self.horizon * n..], n)
    }

    /// Gradient with respect to `u` by the adjoint recursion.
    ///
    /// `cost_weight` scales the cost terms (`1.0` for the gradient of `J`,
    /// `0.0` to keep only `force`) and `force`, if given, holds an extra
    /// gradient `∂P/∂x_k` for every stage `0..=N` (length `(N + 1) n`).
    fn gradient(
        &self,
        traj: &[f64],
        u: &[f64],
        cost_weight: f64,
        force: Option<&[f64]>,
    ) -> Vec<f64> {
        let (n, m, h) = (self.n, self.m, self.horizon);
        let mut grad = vec![0.0; h * m];
        if h == 0 {
            return grad;
        }
        let mut p = vec![0.0; n];
        mat_vec_add(&self.qf2, n, n, &traj[h * n..], &mut p);
        for v in p.iter_mut() {
            *v *= cost_weight;
        }
        if let Some(f) = force {
            for (pi, fi) in p.iter_mut().zip(&f[h * n..]) {
                *pi += fi;
            }
        }
        for k in (0..h).rev() {
            let g = &mut grad[k * m..(k + 1) * m];
            mat_vec_add(&self.r2, m, m, &u[k * m..(k + 1) * m], g);
            for v in g.iter_mut() {
                *v *= cost_weight;
            }
            mat_vec_add(&self.b_t, m, n, &p, g);
            if k >= 1 {
                let mut next = vec![0.0; n];
                mat_vec_add(&self.q2, n, n, &traj[k * n..(k + 1) * n], &mut next);
                for v in next.iter_mut() {
                    *v *= cost_weight;
                }
                if let Some(f) = force {
                    for (ni, fi) in next.iter_mut().zip(&f[k * n..(k + 1) * n]) {
                        *ni += fi;
                    }
                }
                mat_vec_add(&self.a_t, n, n, &p, &mut next);
                p = next;
            }
        }
        grad
    }

    /// Clamp every input to its box.
    fn project(&self, u: &mut [f64]) {
        let m = self.m;
        for (i, v) in u.iter_mut().enumerate() {
            let j = i % m;
            *v = v.max(self.u_lb[j]).min(self.u_ub[j]);
        }
    }

    /// Largest violation of a state bound over the stages `first..=N`.
    fn state_violation(&self, traj: &[f64], first: usize) -> f64 {
        let n = self.n;
        let mut worst = 0.0_f64;
        for k in first..=self.horizon {
            for i in 0..n {
                let x = traj[k * n + i];
                worst = worst.max(self.x_lb[i] - x).max(x - self.x_ub[i]);
            }
        }
        worst
    }

    /// Penalty value `Σ max(0, s)² / (2 ρ)` and, optionally, its gradient with
    /// respect to the states, where `s = λ + ρ c` for every bound `c <= 0`.
    fn penalty(&self, pen: &Penalty, traj: &[f64], mut force: Option<&mut [f64]>) -> f64 {
        let n = self.n;
        let rho = pen.rho;
        let mut value = 0.0;
        for k in 1..=self.horizon {
            for i in 0..n {
                if !self.constrained[i] {
                    continue;
                }
                let x = traj[k * n + i];
                let idx = k * n + i;
                if self.x_ub[i].is_finite() {
                    let s = pen.upper[idx] + rho * (x - self.x_ub[i]);
                    if s > 0.0 {
                        value += s * s / (2.0 * rho);
                        if let Some(f) = force.as_deref_mut() {
                            f[idx] += s;
                        }
                    }
                }
                if self.x_lb[i].is_finite() {
                    let s = pen.lower[idx] + rho * (self.x_lb[i] - x);
                    if s > 0.0 {
                        value += s * s / (2.0 * rho);
                        if let Some(f) = force.as_deref_mut() {
                            f[idx] -= s;
                        }
                    }
                }
            }
        }
        value
    }

    /// Objective of the inner problem and its gradient (cost plus penalty).
    fn objective(
        &self,
        pen: Option<&Penalty>,
        x0: &[f64],
        u: &[f64],
        want_gradient: bool,
    ) -> (f64, Vec<f64>) {
        let traj = self.simulate(x0, u);
        let mut value = self.cost(&traj, u);
        let mut force = pen.map(|_| vec![0.0; (self.horizon + 1) * self.n]);
        if let (Some(p), Some(f)) = (pen, force.as_mut()) {
            value += self.penalty(p, &traj, Some(f.as_mut_slice()));
        }
        let grad = if want_gradient {
            self.gradient(&traj, u, 1.0, force.as_deref())
        } else {
            Vec::new()
        };
        (value, grad)
    }

    /// Dominant eigenvalue of `w H + ρ Γᵀ M Γ` by power iteration, where `H` is
    /// the Hessian of `J`, `Γ` maps the inputs to the states and `M` selects the
    /// constrained states.  `Hv` and `Γᵀ M Γ v` come from the same adjoint pass
    /// applied to the response of the zero state to the input `v`; the cost
    /// weight `w` is `1` or `0` and `rho` is `None` to leave out the state-bound
    /// term.  The estimate approaches the eigenvalue from below.
    fn largest_curvature(&self, cost_weight: f64, rho: Option<f64>) -> f64 {
        let dim = self.horizon * self.m;
        if dim == 0 {
            return 0.0;
        }
        let zero_state = vec![0.0; self.n];
        let mut v: Vec<f64> = (0..dim)
            .map(|i| 1.0 + 0.37 * ((i as f64 + 1.0) * 1.618_033_988_7).sin())
            .collect();
        let start_norm = dot(&v, &v).sqrt();
        for x in v.iter_mut() {
            *x /= start_norm;
        }
        let mut estimate = 0.0_f64;
        for _ in 0..POWER_ITERATIONS {
            let traj = self.simulate(&zero_state, &v);
            let force = rho.map(|r| {
                let mut f = vec![0.0; (self.horizon + 1) * self.n];
                for k in 1..=self.horizon {
                    for i in 0..self.n {
                        if self.constrained[i] {
                            f[k * self.n + i] = r * traj[k * self.n + i];
                        }
                    }
                }
                f
            });
            let hv = self.gradient(&traj, &v, cost_weight, force.as_deref());
            let norm = dot(&hv, &hv).sqrt();
            let rayleigh = dot(&v, &hv);
            let next_estimate = norm.max(rayleigh);
            if norm == 0.0 {
                return 0.0;
            }
            let converged = (next_estimate - estimate).abs() <= 1e-9 * next_estimate;
            estimate = next_estimate;
            for (x, h) in v.iter_mut().zip(hv.iter()) {
                *x = h / norm;
            }
            if converged {
                break;
            }
        }
        estimate
    }
}

/// Stopping rule of the inner solve.
///
/// The gradient-mapping norm at the extrapolated point `y` must not exceed
/// `relative · (gradient_scale + curvature · ‖y‖∞)`.  This is the usual
/// backward-error test of a linear system `H u = -f`: it is invariant under a
/// scaling of the cost (a problem whose weights are all multiplied by `1e-6`
/// is solved to the same relative accuracy) and it never asks for less noise
/// than the rounding error of the gradient itself, which is of the order
/// `ε (‖f‖ + ‖H‖ ‖y‖)`.
#[derive(Clone, Copy)]
struct StoppingRule {
    /// Relative accuracy.
    relative: f64,
    /// `‖∇J(u₀)‖∞` at the starting point: the size of the gradient the solve
    /// has to remove.
    gradient_scale: f64,
    /// Largest curvature of the cost `J` (without any penalty).
    curvature: f64,
}

impl StoppingRule {
    fn threshold(&self, y: &[f64]) -> f64 {
        self.relative * (self.gradient_scale + self.curvature * norm_inf(y))
    }
}

/// Result of one box-constrained inner solve.
struct InnerOutcome {
    iterations: usize,
    converged: bool,
}

/// FISTA with gradient-scheme restart and backtracking for
/// `min F(u)` over the input box, started at `u` (updated in place).
///
/// The momentum is dropped as soon as it points against the gradient mapping,
/// `(y - x⁺)ᵀ (x⁺ - x) > 0` (O'Donoghue and Candès, 2015).  Unlike a restart
/// on an increase of `F`, this test does not depend on function values, whose
/// differences sink below the rounding error long before the gradient is small
/// enough to stop: there the restart would fire at random and reduce the method
/// to plain gradient descent.
///
/// Stops when the gradient-mapping norm `L ‖y - P(y - ∇F(y)/L)‖∞` is at most
/// the threshold of the [`StoppingRule`].  `lipschitz` is the running estimate
/// of the gradient's Lipschitz constant: it only ever grows, by backtracking.
fn accelerated_projected_gradient(
    mpc: &Mpc<'_>,
    pen: Option<&Penalty>,
    x0: &[f64],
    u: &mut Vec<f64>,
    lipschitz: &mut f64,
    rule: StoppingRule,
) -> InnerOutcome {
    mpc.project(u);
    let mut x = u.clone();
    let mut y = x.clone();
    let mut t = 1.0_f64;
    for iteration in 0..MAX_INNER_ITERATIONS {
        let (f_y, g_y) = mpc.objective(pen, x0, &y, true);
        let threshold = rule.threshold(&y);
        let slack = 1e-13 * f_y.abs();
        let x_new = loop {
            let mut candidate: Vec<f64> = y
                .iter()
                .zip(g_y.iter())
                .map(|(yi, gi)| yi - gi / *lipschitz)
                .collect();
            mpc.project(&mut candidate);
            let (f_c, _) = mpc.objective(pen, x0, &candidate, false);
            let step: Vec<f64> = candidate
                .iter()
                .zip(y.iter())
                .map(|(c, yi)| c - yi)
                .collect();
            let bound = f_y + dot(&g_y, &step) + 0.5 * *lipschitz * dot(&step, &step) + slack;
            if f_c <= bound || !lipschitz.is_finite() || *lipschitz > 1e300 {
                break candidate;
            }
            *lipschitz *= 2.0;
        };
        let mapping = *lipschitz
            * x_new
                .iter()
                .zip(y.iter())
                .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
        let t_next = 0.5 * (1.0 + (1.0 + 4.0 * t * t).sqrt());
        let momentum_opposes_mapping: f64 = y
            .iter()
            .zip(x_new.iter())
            .zip(x.iter())
            .map(|((yi, xn), xo)| (yi - xn) * (xn - xo))
            .sum();
        if momentum_opposes_mapping > 0.0 {
            y = x_new.clone();
            t = 1.0;
        } else {
            let beta = (t - 1.0) / t_next;
            y = x_new
                .iter()
                .zip(x.iter())
                .map(|(xn, xo)| xn + beta * (xn - xo))
                .collect();
            t = t_next;
        }
        x = x_new;
        if mapping <= threshold {
            *u = x;
            return InnerOutcome {
                iterations: iteration + 1,
                converged: true,
            };
        }
    }
    *u = x;
    InnerOutcome {
        iterations: MAX_INNER_ITERATIONS,
        converged: false,
    }
}

/// Solve the MPC problem for the initial state `x0`.
///
/// Inconsistent matrix sizes, a wrong-sized or non-finite `x0` yield a failed
/// solution (zero inputs, NaN cost and trajectory, `feasible = false`,
/// `converged = false`) instead of a panic.
pub(super) fn solve(problem: &MpcProblem, x0: &[f64]) -> MpcSolution {
    let failure = || MpcSolution {
        u_seq: vec![0.0; problem.horizon * problem.m],
        x_traj: vec![f64::NAN; (problem.horizon + 1) * problem.n],
        cost: f64::NAN,
        feasible: false,
        iterations: 0,
        converged: false,
    };
    let Some(mpc) = Mpc::new(problem) else {
        return failure();
    };
    if x0.len() != mpc.n
        || x0.iter().any(|v| !v.is_finite())
        || [&problem.a, &problem.b, &problem.q, &problem.qf, &problem.r]
            .iter()
            .any(|matrix| matrix.iter().any(|v| !v.is_finite()))
    {
        return failure();
    }
    let dim = mpc.horizon * mpc.m;
    let mut u = vec![0.0; dim];
    mpc.project(&mut u);
    let mut iterations = 0;
    let mut converged = true;

    if dim > 0 {
        let (_, g0) = mpc.objective(None, x0, &u, true);
        let rule = StoppingRule {
            relative: GRADIENT_TOL,
            gradient_scale: norm_inf(&g0),
            curvature: mpc.largest_curvature(1.0, None),
        };
        if mpc.constrained.iter().any(|&c| c) {
            let (inner, ok) = augmented_lagrangian(&mpc, x0, &mut u, rule);
            iterations = inner;
            converged = ok;
        } else {
            let mut lipschitz = initial_lipschitz(rule.curvature);
            let outcome =
                accelerated_projected_gradient(&mpc, None, x0, &mut u, &mut lipschitz, rule);
            iterations = outcome.iterations;
            converged = outcome.converged;
        }
    }

    let traj = mpc.simulate(x0, &u);
    let cost = mpc.cost(&traj, &u);
    let violation = mpc.state_violation(&traj, 0);
    MpcSolution {
        feasible: violation <= REPORT_FEASIBILITY_TOL && traj.iter().all(|v| v.is_finite()),
        u_seq: u,
        x_traj: traj,
        cost,
        iterations,
        converged,
    }
}

/// Starting Lipschitz estimate from a power-iteration curvature estimate.
fn initial_lipschitz(curvature: f64) -> f64 {
    (LIPSCHITZ_SAFETY * curvature).max(f64::MIN_POSITIVE.sqrt())
}

/// Augmented-Lagrangian loop for the state bounds; returns the total number
/// of inner iterations and whether the bounds were met to [`FEASIBILITY_TOL`]
/// with stationary multipliers.
fn augmented_lagrangian(
    mpc: &Mpc<'_>,
    x0: &[f64],
    u: &mut Vec<f64>,
    rule: StoppingRule,
) -> (usize, bool) {
    let (n, h) = (mpc.n, mpc.horizon);
    // Initial penalty: ten times the ratio of cost to constraint curvature.
    let cost_curvature = rule.curvature;
    let constraint_curvature = mpc.largest_curvature(0.0, Some(1.0));
    let rho0 = if constraint_curvature > 0.0 && cost_curvature > 0.0 {
        10.0 * cost_curvature / constraint_curvature
    } else {
        1.0
    };
    let rho_cap = PENALTY_MAX_FACTOR * rho0;
    let mut pen = Penalty {
        rho: rho0,
        upper: vec![0.0; (h + 1) * n],
        lower: vec![0.0; (h + 1) * n],
    };
    let bound_scale = mpc
        .x_lb
        .iter()
        .chain(mpc.x_ub.iter())
        .filter(|b| b.is_finite())
        .fold(1.0_f64, |m, b| m.max(b.abs()));
    let feasibility_tol = FEASIBILITY_TOL * bound_scale;
    let mut total_iterations = 0;
    let mut previous_violation = f64::INFINITY;
    let mut penalty_just_raised = false;
    let mut converged = false;
    for _ in 0..MAX_OUTER_ITERATIONS {
        let curvature = mpc.largest_curvature(1.0, Some(pen.rho));
        let mut lipschitz = initial_lipschitz(curvature);
        // Inexact inner solves: the accuracy asked for tracks the violation
        // left by the previous outer iteration (a tenth of it, relative to the
        // size of the bounds), down to the final accuracy once the iterate is
        // feasible.  Solving a subproblem much more accurately than the
        // multipliers it will be re-solved with are known would waste the
        // iterations.
        let inner_relative = (0.1 * previous_violation / bound_scale)
            .min(INEXACT_INNER_TOL)
            .max(rule.relative);
        let inner_rule = StoppingRule {
            relative: inner_relative,
            ..rule
        };
        let outcome =
            accelerated_projected_gradient(mpc, Some(&pen), x0, u, &mut lipschitz, inner_rule);
        total_iterations += outcome.iterations;
        let traj = mpc.simulate(x0, u);
        let violation = mpc.state_violation(&traj, 1);
        // Multiplier update λ <- max(0, λ + ρ c).  The updated multiplier pairs
        // with the iterate as a KKT point of the original problem: the inner
        // gradient is its stationarity residual, `violation` its primal
        // residual, and for every bound with a positive multiplier the
        // constraint value `|c|` is the complementarity gap.  (The change of
        // the multipliers itself is `ρ |c|`, which mixes in the penalty and
        // cannot serve as a convergence measure.)
        let mut gap = 0.0_f64;
        for k in 1..=h {
            for i in 0..n {
                let idx = k * n + i;
                if mpc.x_ub[i].is_finite() {
                    let c = traj[idx] - mpc.x_ub[i];
                    let updated = (pen.upper[idx] + pen.rho * c).max(0.0);
                    if updated > 0.0 {
                        gap = gap.max(c.abs());
                    }
                    pen.upper[idx] = updated;
                }
                if mpc.x_lb[i].is_finite() {
                    let c = mpc.x_lb[i] - traj[idx];
                    let updated = (pen.lower[idx] + pen.rho * c).max(0.0);
                    if updated > 0.0 {
                        gap = gap.max(c.abs());
                    }
                    pen.lower[idx] = updated;
                }
            }
        }
        if outcome.converged
            && inner_relative <= rule.relative
            && violation <= feasibility_tol
            && gap <= feasibility_tol
        {
            converged = true;
            break;
        }
        // The penalty only has to grow while the bounds are violated and the
        // violation does not shrink fast enough.  Once the iterate is feasible
        // the multiplier updates alone converge at a fixed penalty, and a
        // larger one would just make the inner problems harder.  The first
        // iteration after a change is not judged: its multipliers were adapted
        // to the old penalty, so its contraction ratio says nothing.
        let mut raised = false;
        if violation > feasibility_tol
            && !penalty_just_raised
            && violation > PENALTY_STALL_RATIO * previous_violation
        {
            if pen.rho < rho_cap {
                pen.rho = (pen.rho * PENALTY_GROWTH).min(rho_cap);
                raised = true;
            } else if violation >= 0.99 * previous_violation {
                // The penalty is exhausted and the violation no longer shrinks:
                // the bounds cannot be met.
                break;
            }
        }
        penalty_just_raised = raised;
        previous_violation = violation;
    }
    (total_iterations, converged)
}

/// Cost and exact adjoint gradient `∇J(u)` of the condensed problem at the
/// input sequence `u` (test support).  `None` for inconsistent dimensions.
#[cfg(test)]
pub(super) fn cost_and_gradient(
    problem: &MpcProblem,
    x0: &[f64],
    u: &[f64],
) -> Option<(f64, Vec<f64>)> {
    let mpc = Mpc::new(problem)?;
    if x0.len() != mpc.n || u.len() != mpc.horizon * mpc.m {
        return None;
    }
    Some(mpc.objective(None, x0, u, true))
}

/// Power-iteration estimate of the largest eigenvalue of `w H + ρ Γᵀ M Γ`
/// (test support): `cost_weight` is `w`, `rho` is `ρ` or `None` to leave out
/// the state-bound term.  `None` for inconsistent dimensions.
#[cfg(test)]
pub(super) fn curvature_estimate(
    problem: &MpcProblem,
    cost_weight: f64,
    rho: Option<f64>,
) -> Option<f64> {
    Mpc::new(problem).map(|mpc| mpc.largest_curvature(cost_weight, rho))
}
