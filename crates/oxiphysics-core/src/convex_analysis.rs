#![allow(clippy::needless_range_loop, clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Convex analysis, proximal operators, and constrained optimisation.
//!
//! Implements convex sets (hyperplane, halfspace, ball, polytope), projection
//! operators, proximal maps, dual decomposition with subgradient updates,
//! ADMM, Frank-Wolfe / conditional gradient, and the ellipsoid method.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─── Convexity Check ─────────────────────────────────────────────────────────

/// Check whether a discretely sampled 1-D function is (approximately) convex.
///
/// Convexity is verified by checking that all finite second differences are
/// ≥ −`tol`.
///
/// Returns `true` for empty or single-element slices.
pub fn is_convex_1d(f: &[f64], tol: f64) -> bool {
    if f.len() < 3 {
        return true;
    }
    for i in 1..f.len() - 1 {
        let second_diff = f[i + 1] - 2.0 * f[i] + f[i - 1];
        if second_diff < -tol {
            return false;
        }
    }
    true
}

// ─── ConvexSet ───────────────────────────────────────────────────────────────

/// A representation of a convex set in ℝⁿ.
///
/// Supports hyperplanes, halfspaces, Euclidean balls, and polytopes in
/// H-representation (intersection of halfspaces).
#[derive(Debug, Clone)]
pub enum ConvexSet {
    /// Hyperplane `{x : aᵀx = b}`.
    Hyperplane {
        /// Normal vector `a`.
        a: Vec<f64>,
        /// Right-hand side `b`.
        b: f64,
    },
    /// Halfspace `{x : aᵀx ≤ b}`.
    Halfspace {
        /// Normal vector `a`.
        a: Vec<f64>,
        /// Right-hand side `b`.
        b: f64,
    },
    /// Euclidean ball `{x : ||x - c|| ≤ r}`.
    Ball {
        /// Centre `c`.
        center: Vec<f64>,
        /// Radius `r`.
        radius: f64,
    },
    /// Polytope in H-form: `{x : A x ≤ b}`.
    ///
    /// `a_rows[i]` is the `i`-th row of A and `b_vec[i]` is `b_i`.
    Polytope {
        /// Rows of the constraint matrix.
        a_rows: Vec<Vec<f64>>,
        /// Right-hand side vector.
        b_vec: Vec<f64>,
    },
}

impl ConvexSet {
    /// Check whether a point `x` belongs to this convex set.
    pub fn contains(&self, x: &[f64]) -> bool {
        match self {
            ConvexSet::Hyperplane { a, b } => {
                let v: f64 = dot(a, x);
                (v - b).abs() < 1e-9
            }
            ConvexSet::Halfspace { a, b } => dot(a, x) <= b + 1e-9,
            ConvexSet::Ball { center, radius } => {
                let d: f64 = x
                    .iter()
                    .zip(center.iter())
                    .map(|(xi, ci)| (xi - ci).powi(2))
                    .sum::<f64>()
                    .sqrt();
                d <= radius + 1e-9
            }
            ConvexSet::Polytope { a_rows, b_vec } => a_rows
                .iter()
                .zip(b_vec.iter())
                .all(|(row, &bi)| dot(row, x) <= bi + 1e-9),
        }
    }
}

/// Compute the inner product of two equal-length slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

// ─── ProjectionOperator ───────────────────────────────────────────────────────

/// Collection of projection operators onto standard convex sets.
///
/// Each method returns the projection of a given point onto a specified set.
pub struct ProjectionOperator;

impl ProjectionOperator {
    /// Project `x` onto the hyperplane `{z : aᵀz = b}`.
    ///
    /// `proj(x) = x - (aᵀx - b) / ||a||² · a`
    pub fn onto_hyperplane(x: &[f64], a: &[f64], b: f64) -> Vec<f64> {
        let a_norm_sq: f64 = dot(a, a);
        if a_norm_sq < 1e-14 {
            return x.to_vec();
        }
        let alpha = (dot(a, x) - b) / a_norm_sq;
        x.iter()
            .zip(a.iter())
            .map(|(xi, ai)| xi - alpha * ai)
            .collect()
    }

    /// Project `x` onto the halfspace `{z : aᵀz ≤ b}`.
    ///
    /// If `x` is already in the halfspace returns `x` unchanged.
    pub fn onto_halfspace(x: &[f64], a: &[f64], b: f64) -> Vec<f64> {
        let ax = dot(a, x);
        if ax <= b {
            return x.to_vec();
        }
        let a_norm_sq: f64 = dot(a, a);
        if a_norm_sq < 1e-14 {
            return x.to_vec();
        }
        let alpha = (ax - b) / a_norm_sq;
        x.iter()
            .zip(a.iter())
            .map(|(xi, ai)| xi - alpha * ai)
            .collect()
    }

    /// Project `x` onto the Euclidean ball `B(center, radius)`.
    pub fn onto_ball(x: &[f64], center: &[f64], radius: f64) -> Vec<f64> {
        let diff: Vec<f64> = x
            .iter()
            .zip(center.iter())
            .map(|(xi, ci)| xi - ci)
            .collect();
        let dist: f64 = diff.iter().map(|d| d * d).sum::<f64>().sqrt();
        if dist <= radius {
            return x.to_vec();
        }
        let scale = radius / dist;
        center
            .iter()
            .zip(diff.iter())
            .map(|(ci, di)| ci + scale * di)
            .collect()
    }

    /// Project `x` onto the axis-aligned box `[lo, hi]` component-wise.
    pub fn onto_box(x: &[f64], lo: &[f64], hi: &[f64]) -> Vec<f64> {
        x.iter()
            .zip(lo.iter().zip(hi.iter()))
            .map(|(&xi, (&li, &hi_i))| xi.max(li).min(hi_i))
            .collect()
    }

    /// Project `x` onto the probability simplex Δ = {z ≥ 0 : sum z = 1}.
    pub fn onto_simplex(x: &[f64]) -> Vec<f64> {
        proj_simplex(x)
    }
}

// ─── Proximal Operators ───────────────────────────────────────────────────────

/// A proximal operator parameterised by a regularisation strength `lambda`.
#[derive(Debug, Clone)]
pub struct ProximalOperator {
    /// Regularisation strength λ > 0.
    pub lambda: f64,
}

impl ProximalOperator {
    /// Create a new proximal operator with the given `lambda`.
    pub fn new(lambda: f64) -> Self {
        Self { lambda }
    }

    /// Apply the L1 (soft-thresholding) proximal operator: `prox_{λ|·|}(x)`.
    pub fn prox_l1(&self, x: f64) -> f64 {
        prox_l1(x, self.lambda)
    }

    /// Apply the L2-squared proximal operator: `prox_{λ||·||²/2}(x)`.
    pub fn prox_l2_sq(&self, x: f64) -> f64 {
        prox_l2_sq(x, self.lambda)
    }

    /// Apply the Huber proximal operator with `delta` parameter.
    pub fn prox_huber_val(&self, x: f64, delta: f64) -> f64 {
        prox_huber(x, self.lambda, delta)
    }

    /// Apply the indicator proximal operator for the box `[lo, hi]`.
    ///
    /// The proximal map of the indicator function of a convex set is just the
    /// projection onto that set.
    pub fn prox_box_indicator(&self, x: f64, lo: f64, hi: f64) -> f64 {
        x.max(lo).min(hi)
    }

    /// Apply the L1 proximal operator elementwise to a vector.
    pub fn prox_l1_vec(&self, x: &[f64]) -> Vec<f64> {
        x.iter().map(|&xi| prox_l1(xi, self.lambda)).collect()
    }

    /// Apply the L2-squared proximal operator elementwise to a vector.
    pub fn prox_l2_sq_vec(&self, x: &[f64]) -> Vec<f64> {
        x.iter().map(|&xi| prox_l2_sq(xi, self.lambda)).collect()
    }
}

/// Soft-thresholding (proximal operator of the L1 norm): `prox_{λ|·|}(x)`.
///
/// Also known as the shrinkage operator: `sign(x) * max(|x| - λ, 0)`.
pub fn prox_l1(x: f64, lambda: f64) -> f64 {
    x.signum() * (x.abs() - lambda).max(0.0)
}

/// Proximal operator of the squared L2 norm (ridge): `x / (1 + λ)`.
pub fn prox_l2_sq(x: f64, lambda: f64) -> f64 {
    x / (1.0 + lambda)
}

/// Proximal operator of the Huber loss with parameter `delta`.
///
/// The Huber loss is:
/// - `x²/2` for `|x| ≤ delta`
/// - `delta * (|x| - delta/2)` for `|x| > delta`
pub fn prox_huber(x: f64, lambda: f64, delta: f64) -> f64 {
    let threshold = delta * (1.0 + lambda);
    if x.abs() <= threshold {
        x / (1.0 + lambda)
    } else {
        x - lambda * delta * x.signum()
    }
}

// ─── Subgradient Descent ─────────────────────────────────────────────────────

/// Configuration for subgradient descent.
#[derive(Debug, Clone)]
pub struct SubgradientMethod {
    /// Step size (learning rate).
    pub step_size: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
}

impl SubgradientMethod {
    /// Create a new subgradient method configuration.
    pub fn new(step_size: f64, max_iter: usize) -> Self {
        Self {
            step_size,
            max_iter,
        }
    }
}

/// Run `iters` steps of subgradient descent on `grad_f` from `x0` with step
/// size `lr`.
///
/// Returns the final iterate.
pub fn subgradient_descent(
    grad_f: impl Fn(&[f64]) -> Vec<f64>,
    mut x0: Vec<f64>,
    lr: f64,
    iters: usize,
) -> Vec<f64> {
    for _ in 0..iters {
        let g = grad_f(&x0);
        for (xi, gi) in x0.iter_mut().zip(g.iter()) {
            *xi -= lr * gi;
        }
    }
    x0
}

// ─── DualDecomposition ───────────────────────────────────────────────────────

/// Lagrangian relaxation with subgradient updates for dual decomposition.
///
/// Minimises `f(x) + g(z)` subject to `Ax + Bz = c` by relaxing the
/// coupling constraint.  The dual variables `y` are updated using the
/// subgradient of the dual function.
pub struct DualDecomposition {
    /// Step size for dual variable updates.
    pub step_size: f64,
    /// Maximum number of dual iterations.
    pub max_iter: usize,
    /// Convergence tolerance on the primal feasibility gap ‖Ax + Bz − c‖.
    pub tol: f64,
}

impl DualDecomposition {
    /// Create a new `DualDecomposition` solver.
    pub fn new(step_size: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            step_size,
            max_iter,
            tol,
        }
    }

    /// Run dual decomposition.
    ///
    /// - `x_solve`: minimiser of the Lagrangian w.r.t. `x` given dual `y` and `z`.
    /// - `z_solve`: minimiser of the Lagrangian w.r.t. `z` given dual `y` and `x`.
    /// - `mat_a`, `mat_b`, `c`: constraint `Ax + Bz = c`.
    ///
    /// Returns `(x_best, z_best, y_final, iterations)`.
    pub fn solve(
        &self,
        x_init: Vec<f64>,
        z_init: Vec<f64>,
        y_init: Vec<f64>,
        x_solve: &dyn Fn(&[f64], &[f64]) -> Vec<f64>,
        z_solve: &dyn Fn(&[f64], &[f64]) -> Vec<f64>,
        mat_a: &[Vec<f64>],
        mat_b: &[Vec<f64>],
        c: &[f64],
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, usize) {
        let m = c.len();
        let mut x = x_init;
        let mut z = z_init;
        let mut y = y_init;
        let mut best_x = x.clone();
        let mut best_z = z.clone();
        let mut best_gap = f64::MAX;

        for iter in 0..self.max_iter {
            // x-update
            x = x_solve(&x, &y);
            // z-update
            z = z_solve(&z, &y);
            // Dual update: y += step * (Ax + Bz - c)
            let ax = matvec(mat_a, &x);
            let bz = matvec(mat_b, &z);
            let mut gap = 0.0_f64;
            for i in 0..m {
                let subgrad = ax[i] + bz[i] - c[i];
                gap += subgrad * subgrad;
                y[i] += self.step_size * subgrad;
            }
            gap = gap.sqrt();
            if gap < best_gap {
                best_gap = gap;
                best_x = x.clone();
                best_z = z.clone();
            }
            if gap < self.tol {
                return (best_x, best_z, y, iter + 1);
            }
        }
        (best_x, best_z, y, self.max_iter)
    }
}

/// Multiply matrix `a` (rows × cols) by vector `x`.
fn matvec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

// ─── ADMMSolver ──────────────────────────────────────────────────────────────

/// Result returned by the ADMM solver.
#[derive(Debug, Clone)]
pub struct AdmmResult {
    /// Primal variable `x`.
    pub x: Vec<f64>,
    /// Split variable `z`.
    pub z: Vec<f64>,
    /// Dual variable `u` (scaled form).
    pub u: Vec<f64>,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the algorithm converged within tolerance.
    pub converged: bool,
    /// Final primal residual ‖x − z‖.
    pub primal_residual: f64,
}

/// Alternating Direction Method of Multipliers (ADMM) for problems of the
/// form:
///
/// ```text
///   minimise   f(x) + g(z)
///   subject to x = z
/// ```
///
/// Uses the augmented Lagrangian with penalty parameter `rho`.
pub struct AdmmSolver {
    /// Augmented Lagrangian penalty ρ > 0.
    pub rho: f64,
    /// Maximum number of ADMM iterations.
    pub max_iter: usize,
    /// Convergence tolerance on primal and dual residuals.
    pub tol: f64,
}

impl AdmmSolver {
    /// Create a new ADMM solver.
    pub fn new(rho: f64, max_iter: usize, tol: f64) -> Self {
        Self { rho, max_iter, tol }
    }

    /// Solve the consensus ADMM problem.
    ///
    /// - `x_update`: solves `argmin_x f(x) + (ρ/2)||x − z + u||²`.
    /// - `z_update`: solves `argmin_z g(z) + (ρ/2)||x − z + u||²`.
    pub fn solve(
        &self,
        x_init: Vec<f64>,
        x_update: &dyn Fn(&[f64], &[f64], f64) -> Vec<f64>,
        z_update: &dyn Fn(&[f64], &[f64], f64) -> Vec<f64>,
    ) -> AdmmResult {
        let n = x_init.len();
        let mut x = x_init;
        let mut z = vec![0.0_f64; n];
        let mut u = vec![0.0_f64; n];

        for iter in 0..self.max_iter {
            // x-update
            let x_arg: Vec<f64> = z.iter().zip(u.iter()).map(|(zi, ui)| zi - ui).collect();
            x = x_update(&x, &x_arg, self.rho);

            // z-update
            let z_arg: Vec<f64> = x.iter().zip(u.iter()).map(|(xi, ui)| xi + ui).collect();
            let z_new = z_update(&z, &z_arg, self.rho);

            // Dual update
            let mut primal_res = 0.0_f64;
            let mut dual_res = 0.0_f64;
            for i in 0..n {
                u[i] += x[i] - z_new[i];
                primal_res += (x[i] - z_new[i]).powi(2);
                dual_res += self.rho * (z_new[i] - z[i]).powi(2);
            }
            primal_res = primal_res.sqrt();
            dual_res = dual_res.sqrt();
            z = z_new;

            if primal_res < self.tol && dual_res < self.tol {
                return AdmmResult {
                    x,
                    z,
                    u,
                    iterations: iter + 1,
                    converged: true,
                    primal_residual: primal_res,
                };
            }
        }

        let pr: f64 = x
            .iter()
            .zip(z.iter())
            .map(|(xi, zi)| (xi - zi).powi(2))
            .sum::<f64>()
            .sqrt();
        AdmmResult {
            x,
            z,
            u,
            iterations: self.max_iter,
            converged: false,
            primal_residual: pr,
        }
    }
}

// ─── FrankWolfeOptimizer ──────────────────────────────────────────────────────

/// Result from the Frank-Wolfe algorithm.
#[derive(Debug, Clone)]
pub struct FrankWolfeResult {
    /// Approximate minimiser.
    pub x: Vec<f64>,
    /// Final objective value.
    pub objective: f64,
    /// Number of iterations.
    pub iterations: usize,
    /// Whether the duality gap dropped below tolerance.
    pub converged: bool,
}

/// Frank-Wolfe (conditional gradient) method for constrained convex problems.
///
/// Minimises a smooth convex function `f` over a compact convex feasible set
/// by solving a linear minimisation oracle at each step.
pub struct FrankWolfeOptimizer {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Tolerance on the Frank-Wolfe duality gap.
    pub tol: f64,
}

impl FrankWolfeOptimizer {
    /// Create a new Frank-Wolfe optimiser.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// Run Frank-Wolfe minimisation.
    ///
    /// - `grad_f`: gradient oracle `∇f(x)`.
    /// - `lmo`: linear minimisation oracle: given `d`, return `argmin_{s ∈ C} d·s`.
    /// - `f_val`: function value oracle `f(x)`.
    /// - `x_init`: feasible starting point.
    pub fn minimize(
        &self,
        x_init: Vec<f64>,
        grad_f: &dyn Fn(&[f64]) -> Vec<f64>,
        lmo: &dyn Fn(&[f64]) -> Vec<f64>,
        f_val: &dyn Fn(&[f64]) -> f64,
    ) -> FrankWolfeResult {
        let n = x_init.len();
        let mut x = x_init;

        for iter in 0..self.max_iter {
            let g = grad_f(&x);
            let s = lmo(&g);

            // Duality gap: g · (x − s)
            let gap: f64 = g
                .iter()
                .zip(x.iter().zip(s.iter()))
                .map(|(gi, (xi, si))| gi * (xi - si))
                .sum();

            if gap < self.tol {
                return FrankWolfeResult {
                    objective: f_val(&x),
                    x,
                    iterations: iter + 1,
                    converged: true,
                };
            }

            // Step size: 2 / (iter + 2)
            let step = 2.0_f64 / (iter as f64 + 2.0);
            for i in 0..n {
                x[i] = (1.0 - step) * x[i] + step * s[i];
            }
        }

        FrankWolfeResult {
            objective: f_val(&x),
            x,
            iterations: self.max_iter,
            converged: false,
        }
    }
}

// ─── EllipsoidMethod ─────────────────────────────────────────────────────────

/// Result from the ellipsoid method.
#[derive(Debug, Clone)]
pub struct EllipsoidResult {
    /// Approximate centre of feasibility.
    pub x: Vec<f64>,
    /// Whether a feasible point was found (all constraints satisfied).
    pub feasible: bool,
    /// Number of iterations performed.
    pub iterations: usize,
}

/// Ellipsoid method for convex feasibility problems.
///
/// Given a list of constraint functions `g_i(x) ≤ 0` and a starting ellipsoid
/// `E(x0, R²I)`, iteratively cuts the ellipsoid using the most-violated
/// constraint until a feasible point is found or `max_iter` is exhausted.
pub struct EllipsoidMethod {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Initial radius of the starting ball.
    pub initial_radius: f64,
    /// Tolerance for constraint satisfaction.
    pub tol: f64,
}

impl EllipsoidMethod {
    /// Create a new ellipsoid method solver.
    pub fn new(max_iter: usize, initial_radius: f64, tol: f64) -> Self {
        Self {
            max_iter,
            initial_radius,
            tol,
        }
    }

    /// Find a feasible point satisfying all `constraints`.
    ///
    /// Each element of `constraints` is a pair `(g, subgrad)` where `g(x)`
    /// returns the constraint value and `subgrad(x)` returns a subgradient.
    pub fn find_feasible(
        &self,
        x0: Vec<f64>,
        constraints: &[(Box<dyn Fn(&[f64]) -> f64>, Box<dyn Fn(&[f64]) -> Vec<f64>>)],
    ) -> EllipsoidResult {
        let n = x0.len();
        let mut xc = x0;
        // Represent ellipsoid as E = { x : (x-xc)^T P^{-1} (x-xc) ≤ 1 }
        // Store P (n×n symmetric positive definite), initially P = R²·I.
        let r2 = self.initial_radius * self.initial_radius;
        let mut p: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let mut row = vec![0.0_f64; n];
                row[i] = r2;
                row
            })
            .collect();

        for iter in 0..self.max_iter {
            // Check feasibility
            let all_feasible = constraints.iter().all(|(g, _)| g(&xc) <= self.tol);
            if all_feasible {
                return EllipsoidResult {
                    x: xc,
                    feasible: true,
                    iterations: iter,
                };
            }

            // Find the most-violated constraint
            let (most_viol_idx, _max_viol) = constraints
                .iter()
                .enumerate()
                .map(|(i, (g, _))| (i, g(&xc)))
                .filter(|(_, v)| *v > self.tol)
                .fold((0usize, f64::NEG_INFINITY), |(bi, bv), (i, v)| {
                    if v > bv { (i, v) } else { (bi, bv) }
                });

            let sg = constraints[most_viol_idx].1(&xc);
            // Compute g_hat = P * sg / sqrt(sg^T P sg)
            let p_sg = matvec(&p, &sg);
            let sg_p_sg: f64 = dot(&sg, &p_sg);
            if sg_p_sg < 1e-20 {
                break;
            }
            let denom = sg_p_sg.sqrt();
            let g_hat: Vec<f64> = p_sg.iter().map(|x| x / denom).collect();

            let nf = n as f64;
            // Update centre: xc_new = xc - (1/(n+1)) * g_hat
            for i in 0..n {
                xc[i] -= g_hat[i] / (nf + 1.0);
            }
            // Update P: P_new = (n²/(n²-1)) * (P - (2/(n+1)) * g_hat * g_hat^T)
            let scale = nf * nf / (nf * nf - 1.0);
            let rank1_scale = 2.0 / (nf + 1.0);
            for i in 0..n {
                for j in 0..n {
                    p[i][j] = scale * (p[i][j] - rank1_scale * g_hat[i] * g_hat[j]);
                }
            }
        }

        let all_feasible = constraints.iter().all(|(g, _)| g(&xc) <= self.tol);
        EllipsoidResult {
            x: xc,
            feasible: all_feasible,
            iterations: self.max_iter,
        }
    }
}

// ─── Legendre–Fenchel Conjugate ───────────────────────────────────────────────

/// Compute the Legendre–Fenchel conjugate `f*(y) = sup_{x ∈ domain} (y·x - f(x))`
/// via a grid search over `n` equally-spaced points on `domain = (lo, hi)`.
///
/// # Panics
/// Panics if `domain.0 >= domain.1` or `n == 0`.
pub fn conjugate_function_1d(f: impl Fn(f64) -> f64, y: f64, domain: (f64, f64), n: usize) -> f64 {
    assert!(domain.0 < domain.1, "domain must be non-empty");
    assert!(n > 0, "n must be at least 1");
    let (lo, hi) = domain;
    let step = (hi - lo) / (n as f64 - 1.0).max(1.0);
    (0..n)
        .map(|i| {
            let x = lo + i as f64 * step;
            y * x - f(x)
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

// ─── Support Function ─────────────────────────────────────────────────────────

/// Support function of a convex polygon given by its vertices in ℝ².
///
/// `h_C(d) = max_{v ∈ C} d · v`.
#[derive(Debug, Clone)]
pub struct SupportFunction {
    /// Vertices of the convex set.
    pub vertices: Vec<[f64; 2]>,
}

impl SupportFunction {
    /// Construct from a list of vertices.
    pub fn new(vertices: Vec<[f64; 2]>) -> Self {
        Self { vertices }
    }

    /// Evaluate the support function in the given `direction`.
    ///
    /// Returns `f64::NEG_INFINITY` for an empty vertex set.
    pub fn evaluate(&self, direction: [f64; 2]) -> f64 {
        self.vertices
            .iter()
            .map(|v| v[0] * direction[0] + v[1] * direction[1])
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

// ─── Dykstra's Alternating Projections ───────────────────────────────────────

/// Project `x` onto the intersection of convex sets using Dykstra's algorithm.
///
/// Each element of `sets` is a projection operator `π_i: ℝⁿ → ℝⁿ`.
/// Runs `iters` passes over all sets.
pub fn dykstra_projection(
    x: &[f64],
    sets: &[&dyn Fn(&[f64]) -> Vec<f64>],
    iters: usize,
) -> Vec<f64> {
    let n = x.len();
    if sets.is_empty() {
        return x.to_vec();
    }
    let m = sets.len();
    let mut y = x.to_vec();
    let mut increments: Vec<Vec<f64>> = vec![vec![0.0; n]; m];

    for _ in 0..iters {
        for (i, proj) in sets.iter().enumerate() {
            let z: Vec<f64> = y
                .iter()
                .zip(increments[i].iter())
                .map(|(&yi, &ii)| yi + ii)
                .collect();
            let p = proj(&z);
            for k in 0..n {
                increments[i][k] = z[k] - p[k];
            }
            y = p;
        }
    }
    y
}

// ─── Box Projection ───────────────────────────────────────────────────────────

/// Project `x` onto the axis-aligned box `[lo, hi]` component-wise.
///
/// # Panics
/// Panics if `lo`, `hi`, and `x` have different lengths.
pub fn proj_box(x: &[f64], lo: &[f64], hi: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), lo.len());
    assert_eq!(x.len(), hi.len());
    x.iter()
        .zip(lo.iter().zip(hi.iter()))
        .map(|(&xi, (&li, &hi_i))| xi.max(li).min(hi_i))
        .collect()
}

// ─── Simplex Projection ───────────────────────────────────────────────────────

/// Project `x` onto the probability simplex Δ = {z ≥ 0 : sum z = 1}.
///
/// Uses the O(n log n) algorithm by Duchi et al. (2008).
pub fn proj_simplex(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let mut sorted = x.to_vec();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    let mut cumsum = 0.0_f64;
    let mut rho = 0usize;
    for (i, &si) in sorted.iter().enumerate() {
        cumsum += si;
        if si - (cumsum - 1.0) / (i as f64 + 1.0) > 0.0 {
            rho = i;
        }
    }

    let sum_rho: f64 = sorted[..=rho].iter().sum();
    let theta = (sum_rho - 1.0) / (rho as f64 + 1.0);

    x.iter().map(|&xi| (xi - theta).max(0.0)).collect()
}

// ─── Moreau Envelope ─────────────────────────────────────────────────────────

/// Moreau envelope (inf-convolution) smoothing of a convex function.
///
/// The Moreau envelope `M_μf(x) = inf_z [f(z) + ||x - z||²/(2μ)]` provides a
/// smooth approximation to `f`.
#[derive(Debug, Clone)]
pub struct MorseEnvelope {
    /// Smoothing parameter μ > 0.
    pub mu: f64,
}

impl MorseEnvelope {
    /// Create a Moreau-envelope smoother with parameter `mu`.
    pub fn new(mu: f64) -> Self {
        Self { mu }
    }

    /// Evaluate the Moreau envelope of |·| at `x`.
    pub fn eval_l1(&self, x: f64) -> f64 {
        moreau_envelope_l1(x, self.mu)
    }
}

/// Moreau envelope of the L1 norm: smooth approximation to `|x|`.
///
/// `M_μ|·|(x) = |x| - μ/2` if `|x| ≥ μ`, else `x²/(2μ)`.
pub fn moreau_envelope_l1(x: f64, mu: f64) -> f64 {
    if x.abs() >= mu {
        x.abs() - mu / 2.0
    } else {
        x * x / (2.0 * mu)
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- is_convex_1d ---

    #[test]
    fn test_convex_quadratic() {
        let f: Vec<f64> = (-10..=10).map(|i| (i as f64).powi(2)).collect();
        assert!(is_convex_1d(&f, 1e-9));
    }

    #[test]
    fn test_non_convex_neg_quadratic() {
        let f: Vec<f64> = (-10..=10).map(|i| -(i as f64).powi(2)).collect();
        assert!(!is_convex_1d(&f, 1e-9));
    }

    #[test]
    fn test_convex_empty() {
        assert!(is_convex_1d(&[], 1e-9));
    }

    #[test]
    fn test_convex_single() {
        assert!(is_convex_1d(&[3.0], 1e-9));
    }

    #[test]
    fn test_convex_linear() {
        let f: Vec<f64> = (0..20).map(|i| 3.0 * i as f64 + 1.0).collect();
        assert!(is_convex_1d(&f, 1e-9));
    }

    // --- ConvexSet::contains ---

    #[test]
    fn test_convex_set_hyperplane_contains() {
        let hs = ConvexSet::Hyperplane {
            a: vec![1.0, 0.0],
            b: 3.0,
        };
        assert!(hs.contains(&[3.0, 5.0]));
        assert!(!hs.contains(&[2.0, 5.0]));
    }

    #[test]
    fn test_convex_set_halfspace_contains() {
        let hs = ConvexSet::Halfspace {
            a: vec![1.0, 1.0],
            b: 2.0,
        };
        assert!(hs.contains(&[0.5, 0.5]));
        assert!(!hs.contains(&[2.0, 2.0]));
    }

    #[test]
    fn test_convex_set_ball_contains() {
        let ball = ConvexSet::Ball {
            center: vec![0.0, 0.0],
            radius: 1.0,
        };
        assert!(ball.contains(&[0.5, 0.5]));
        assert!(!ball.contains(&[1.0, 1.0]));
    }

    #[test]
    fn test_convex_set_polytope_contains() {
        // Unit cube [0,1]²
        let poly = ConvexSet::Polytope {
            a_rows: vec![
                vec![1.0, 0.0],
                vec![-1.0, 0.0],
                vec![0.0, 1.0],
                vec![0.0, -1.0],
            ],
            b_vec: vec![1.0, 0.0, 1.0, 0.0],
        };
        assert!(poly.contains(&[0.5, 0.5]));
        assert!(!poly.contains(&[1.5, 0.5]));
    }

    // --- ProjectionOperator ---

    #[test]
    fn test_proj_onto_hyperplane() {
        // Hyperplane: x + y = 1; project (0, 0)
        let a = vec![1.0_f64, 1.0];
        let p = ProjectionOperator::onto_hyperplane(&[0.0, 0.0], &a, 1.0);
        let val: f64 = p[0] + p[1];
        assert!(
            (val - 1.0).abs() < 1e-10,
            "projected point must satisfy a·x = b"
        );
    }

    #[test]
    fn test_proj_onto_halfspace_inside() {
        let a = vec![1.0_f64, 0.0];
        let x = vec![0.5_f64, 1.0];
        let p = ProjectionOperator::onto_halfspace(&x, &a, 1.0);
        assert!((p[0] - 0.5).abs() < 1e-12); // already inside
    }

    #[test]
    fn test_proj_onto_halfspace_outside() {
        let a = vec![1.0_f64, 0.0];
        let x = vec![2.0_f64, 1.0];
        let p = ProjectionOperator::onto_halfspace(&x, &a, 1.0);
        assert!((p[0] - 1.0).abs() < 1e-10); // projected to boundary
    }

    #[test]
    fn test_proj_onto_ball_inside() {
        let x = [0.1_f64, 0.1];
        let p = ProjectionOperator::onto_ball(&x, &[0.0, 0.0], 1.0);
        assert!((p[0] - 0.1).abs() < 1e-12);
        assert!((p[1] - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_proj_onto_ball_outside() {
        let x = [3.0_f64, 4.0];
        let p = ProjectionOperator::onto_ball(&x, &[0.0, 0.0], 1.0);
        let dist: f64 = (p[0].powi(2) + p[1].powi(2)).sqrt();
        assert!(
            (dist - 1.0).abs() < 1e-10,
            "projected point must lie on the ball boundary"
        );
    }

    #[test]
    fn test_proj_onto_box() {
        let x = vec![-1.0_f64, 0.5, 3.0];
        let lo = vec![0.0_f64, 0.0, 0.0];
        let hi = vec![1.0_f64, 1.0, 1.0];
        let p = ProjectionOperator::onto_box(&x, &lo, &hi);
        assert_eq!(p, vec![0.0, 0.5, 1.0]);
    }

    #[test]
    fn test_proj_onto_simplex_via_operator() {
        let x = vec![3.0_f64, 1.0, -1.0, 0.5];
        let p = ProjectionOperator::onto_simplex(&x);
        let sum: f64 = p.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    // --- prox_l1 ---

    #[test]
    fn test_prox_l1_positive_above_threshold() {
        assert!((prox_l1(3.0, 1.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_prox_l1_negative_above_threshold() {
        assert!((prox_l1(-3.0, 1.0) + 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_prox_l1_below_threshold() {
        assert_eq!(prox_l1(0.5, 1.0), 0.0);
    }

    #[test]
    fn test_prox_l1_at_threshold() {
        assert_eq!(prox_l1(1.0, 1.0), 0.0);
    }

    #[test]
    fn test_prox_l1_zero() {
        assert_eq!(prox_l1(0.0, 0.5), 0.0);
    }

    // --- prox_l2_sq ---

    #[test]
    fn test_prox_l2_sq_basic() {
        assert!((prox_l2_sq(3.0, 2.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_prox_l2_sq_lambda_zero() {
        assert!((prox_l2_sq(5.0, 0.0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_prox_l2_sq_negative() {
        assert!((prox_l2_sq(-4.0, 1.0) + 2.0).abs() < 1e-12);
    }

    // --- prox_huber ---

    #[test]
    fn test_prox_huber_small_x() {
        let val = prox_huber(0.5, 1.0, 1.0);
        assert!((val - 0.25).abs() < 1e-12, "val = {}", val);
    }

    #[test]
    fn test_prox_huber_large_x() {
        let val = prox_huber(5.0, 1.0, 1.0);
        assert!((val - 4.0).abs() < 1e-12, "val = {}", val);
    }

    // --- ProximalOperator struct ---

    #[test]
    fn test_proximal_operator_l1() {
        let po = ProximalOperator::new(2.0);
        assert!((po.prox_l1(5.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_proximal_operator_l2_sq() {
        let po = ProximalOperator::new(1.0);
        assert!((po.prox_l2_sq(4.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_proximal_operator_box_indicator() {
        let po = ProximalOperator::new(1.0);
        assert!((po.prox_box_indicator(3.0, 0.0, 2.0) - 2.0).abs() < 1e-12);
        assert!((po.prox_box_indicator(-1.0, 0.0, 2.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_proximal_operator_l1_vec() {
        let po = ProximalOperator::new(1.0);
        let x = vec![3.0, -3.0, 0.5, 0.0];
        let p = po.prox_l1_vec(&x);
        assert!((p[0] - 2.0).abs() < 1e-12);
        assert!((p[1] + 2.0).abs() < 1e-12);
        assert_eq!(p[2], 0.0);
        assert_eq!(p[3], 0.0);
    }

    // --- subgradient_descent ---

    #[test]
    fn test_subgradient_descent_quadratic() {
        let result = subgradient_descent(|x| vec![2.0 * (x[0] - 3.0)], vec![0.0], 0.1, 500);
        assert!((result[0] - 3.0).abs() < 0.2, "result = {}", result[0]);
    }

    #[test]
    fn test_subgradient_descent_zero_iters() {
        let result = subgradient_descent(|_| vec![1.0], vec![5.0], 0.1, 0);
        assert!((result[0] - 5.0).abs() < 1e-12);
    }

    // --- ADMM solver ---

    #[test]
    fn test_admm_lasso_converges() {
        // Minimise (1/2)||x - v||² + lambda||z||_1  s.t. x = z
        // x-update: argmin_x (1/2)||x-v||² + (rho/2)||x - (z-u)||²
        //         = (v + rho*(z-u)) / (1 + rho)
        // z-update: prox_{lambda/rho}(x + u) = soft-threshold
        // Closed-form solution: x* = z* = soft_threshold(v, lambda)
        let v = vec![3.0_f64, -2.0];
        let rho = 1.0_f64;
        let lambda = 0.5_f64;
        let v_clone1 = v.clone();
        let v_clone2 = v.clone();
        let admm = AdmmSolver::new(rho, 500, 1e-6);
        let x_upd = move |_x: &[f64], z_minus_u: &[f64], rho_val: f64| -> Vec<f64> {
            v_clone1
                .iter()
                .zip(z_minus_u.iter())
                .map(|(&vi, &zui)| (vi + rho_val * zui) / (1.0 + rho_val))
                .collect()
        };
        let z_upd = move |_z: &[f64], x_plus_u: &[f64], rho_val: f64| -> Vec<f64> {
            let _ = v_clone2;
            x_plus_u
                .iter()
                .map(|&xi| prox_l1(xi, lambda / rho_val))
                .collect()
        };
        let x_init = vec![0.0_f64, 0.0];
        let result = admm.solve(x_init, &x_upd, &z_upd);
        // Closed-form solution: soft_threshold(v, lambda)
        let expected: Vec<f64> = v.iter().map(|&vi| prox_l1(vi, lambda)).collect();
        for i in 0..2 {
            assert!(
                (result.x[i] - expected[i]).abs() < 0.1,
                "x[{}] = {}, expected {}",
                i,
                result.x[i],
                expected[i]
            );
        }
    }

    #[test]
    fn test_admm_result_fields() {
        let admm = AdmmSolver::new(1.0, 10, 1e-8);
        let x_upd = |_x: &[f64], z_u: &[f64], _rho: f64| z_u.to_vec();
        let z_upd = |_z: &[f64], x_u: &[f64], _rho: f64| x_u.to_vec();
        let res = admm.solve(vec![0.0], &x_upd, &z_upd);
        assert!(res.iterations > 0);
    }

    // --- Frank-Wolfe on simplex ---

    #[test]
    fn test_frank_wolfe_simplex_quadratic() {
        // Minimise ||x - c||² over the simplex, c = [0.3, 0.7]
        let c = [0.3_f64, 0.7];
        let grad_f = |x: &[f64]| -> Vec<f64> {
            x.iter()
                .zip(c.iter())
                .map(|(xi, ci)| 2.0 * (xi - ci))
                .collect()
        };
        let lmo = |d: &[f64]| -> Vec<f64> {
            // argmin d·s over simplex = standard basis vector for argmin d_i
            let n = d.len();
            let min_idx = d
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);
            let mut s = vec![0.0_f64; n];
            s[min_idx] = 1.0;
            s
        };
        let f_val = |x: &[f64]| -> f64 {
            x.iter()
                .zip(c.iter())
                .map(|(xi, ci)| (xi - ci).powi(2))
                .sum()
        };
        let fw = FrankWolfeOptimizer::new(500, 1e-6);
        let x0 = vec![0.5_f64, 0.5];
        let result = fw.minimize(x0, &grad_f, &lmo, &f_val);
        // Optimal point on simplex closest to c = (0.3, 0.7) is exactly c (it's already on simplex)
        let sum: f64 = result.x.iter().sum();
        assert!((sum - 1.0).abs() < 0.01, "sum = {}", sum);
        for xi in &result.x {
            assert!(*xi >= -1e-10, "simplex component negative");
        }
    }

    #[test]
    fn test_frank_wolfe_result_converged_flag() {
        let grad_f = |x: &[f64]| vec![2.0 * x[0]];
        let lmo = |_d: &[f64]| vec![0.0_f64];
        let f_val = |x: &[f64]| x[0].powi(2);
        let fw = FrankWolfeOptimizer::new(1000, 1e-6);
        let result = fw.minimize(vec![1.0], &grad_f, &lmo, &f_val);
        let _ = result.converged; // field exists and is accessible
        assert!(result.objective >= 0.0);
    }

    // --- EllipsoidMethod ---

    #[test]
    fn test_ellipsoid_ball_feasibility() {
        // Find x ∈ ℝ² with ||x - (1,1)||² ≤ 0.5 and x[0] ≥ 0
        let constraint1: Box<dyn Fn(&[f64]) -> f64> =
            Box::new(|x: &[f64]| (x[0] - 1.0).powi(2) + (x[1] - 1.0).powi(2) - 0.5);
        let sg1: Box<dyn Fn(&[f64]) -> Vec<f64>> =
            Box::new(|x: &[f64]| vec![2.0 * (x[0] - 1.0), 2.0 * (x[1] - 1.0)]);
        let constraint2: Box<dyn Fn(&[f64]) -> f64> = Box::new(|x: &[f64]| -x[0]);
        let sg2: Box<dyn Fn(&[f64]) -> Vec<f64>> = Box::new(|_x: &[f64]| vec![-1.0, 0.0]);
        let ellipsoid = EllipsoidMethod::new(200, 5.0, 1e-4);
        let result =
            ellipsoid.find_feasible(vec![0.0, 0.0], &[(constraint1, sg1), (constraint2, sg2)]);
        assert!(
            result.feasible,
            "ellipsoid method should find feasible point"
        );
    }

    #[test]
    fn test_ellipsoid_result_fields() {
        let c: Box<dyn Fn(&[f64]) -> f64> = Box::new(|x: &[f64]| x[0] - 1.0);
        let sg: Box<dyn Fn(&[f64]) -> Vec<f64>> = Box::new(|_x: &[f64]| vec![1.0]);
        let em = EllipsoidMethod::new(50, 10.0, 1e-4);
        let result = em.find_feasible(vec![0.5], &[(c, sg)]);
        assert!(result.iterations > 0 || result.feasible);
    }

    // --- conjugate_function_1d ---

    #[test]
    fn test_conjugate_x_squared_half() {
        let y = 2.0_f64;
        let conj = conjugate_function_1d(|x| x * x / 2.0, y, (-10.0, 10.0), 10_000);
        assert!((conj - y * y / 2.0).abs() < 0.01, "conj = {}", conj);
    }

    #[test]
    fn test_conjugate_abs() {
        let conj = conjugate_function_1d(|x| x.abs(), 0.5, (-5.0, 5.0), 1000);
        assert!(conj >= -0.01, "f*(0.5) should be ≥ 0, got {}", conj);
    }

    // --- SupportFunction ---

    #[test]
    fn test_support_function_unit_square() {
        let verts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let sf = SupportFunction::new(verts);
        assert!((sf.evaluate([1.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!((sf.evaluate([1.0, 1.0]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_support_function_empty() {
        let sf = SupportFunction::new(vec![]);
        assert!(sf.evaluate([1.0, 0.0]).is_infinite());
    }

    // --- proj_box ---

    #[test]
    fn test_proj_box_clamps() {
        let x = vec![-1.0, 0.5, 3.0];
        let lo = vec![0.0, 0.0, 0.0];
        let hi = vec![1.0, 1.0, 1.0];
        let p = proj_box(&x, &lo, &hi);
        assert_eq!(p, vec![0.0, 0.5, 1.0]);
    }

    #[test]
    fn test_proj_box_no_clamp() {
        let x = vec![0.3, 0.7];
        let lo = vec![0.0, 0.0];
        let hi = vec![1.0, 1.0];
        let p = proj_box(&x, &lo, &hi);
        assert!((p[0] - 0.3).abs() < 1e-12);
        assert!((p[1] - 0.7).abs() < 1e-12);
    }

    // --- proj_simplex ---

    #[test]
    fn test_proj_simplex_sums_to_one() {
        let x = vec![3.0, 1.0, -1.0, 0.5];
        let p = proj_simplex(&x);
        let sum: f64 = p.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "sum = {}", sum);
    }

    #[test]
    fn test_proj_simplex_non_negative() {
        let x = vec![-5.0, -3.0, -1.0];
        let p = proj_simplex(&x);
        for &pi in &p {
            assert!(pi >= -1e-12, "component {} is negative", pi);
        }
    }

    #[test]
    fn test_proj_simplex_already_on_simplex() {
        let x = vec![0.25, 0.25, 0.25, 0.25];
        let p = proj_simplex(&x);
        let sum: f64 = p.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
        for (&xi, &pi) in x.iter().zip(p.iter()) {
            assert!((xi - pi).abs() < 1e-10);
        }
    }

    #[test]
    fn test_proj_simplex_empty() {
        let p = proj_simplex(&[]);
        assert!(p.is_empty());
    }

    #[test]
    fn test_proj_simplex_single() {
        let p = proj_simplex(&[3.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
    }

    // --- dykstra_projection ---

    #[test]
    fn test_dykstra_no_sets() {
        let x = vec![1.0, 2.0, 3.0];
        let result = dykstra_projection(&x, &[], 10);
        assert_eq!(result, x);
    }

    #[test]
    fn test_dykstra_box_projection() {
        let x = vec![-0.5, 1.5];
        let proj_box_fn =
            |z: &[f64]| -> Vec<f64> { z.iter().map(|&xi| xi.clamp(0.0, 1.0)).collect() };
        let sets: Vec<&dyn Fn(&[f64]) -> Vec<f64>> = vec![&proj_box_fn];
        let result = dykstra_projection(&x, &sets, 5);
        assert!((result[0] - 0.0).abs() < 1e-10);
        assert!((result[1] - 1.0).abs() < 1e-10);
    }

    // --- moreau_envelope_l1 ---

    #[test]
    fn test_moreau_envelope_l1_large_x() {
        let val = moreau_envelope_l1(2.0, 0.5);
        assert!((val - (2.0 - 0.25)).abs() < 1e-12, "val = {}", val);
    }

    #[test]
    fn test_moreau_envelope_l1_small_x() {
        let val = moreau_envelope_l1(0.5, 2.0);
        assert!((val - (0.25 / 4.0)).abs() < 1e-12, "val = {}", val);
    }

    #[test]
    fn test_moreau_envelope_l1_non_negative() {
        for x in [-3.0, -1.0, 0.0, 0.5, 2.0] {
            let val = moreau_envelope_l1(x, 1.0);
            assert!(
                val >= 0.0,
                "Moreau envelope of |·| should be ≥ 0, got {}",
                val
            );
        }
    }

    #[test]
    fn test_moreau_envelope_l1_continuous_at_mu() {
        let mu = 1.5;
        let from_above = moreau_envelope_l1(mu + 1e-9, mu);
        let from_inside = moreau_envelope_l1(mu - 1e-9, mu);
        assert!(
            (from_above - from_inside).abs() < 1e-6,
            "discontinuity at mu"
        );
    }

    #[test]
    fn test_morse_envelope_struct() {
        let me = MorseEnvelope::new(1.0);
        let v = me.eval_l1(2.0);
        assert!((v - 1.5).abs() < 1e-12);
    }

    // --- SubgradientMethod struct ---

    #[test]
    fn test_subgradient_method_struct() {
        let sm = SubgradientMethod::new(0.01, 100);
        assert!((sm.step_size - 0.01).abs() < 1e-12);
        assert_eq!(sm.max_iter, 100);
    }

    // --- DualDecomposition ---

    #[test]
    fn test_dual_decomposition_basic() {
        // Simple consensus: min x² + z², s.t. x = z (A=I, B=-I, c=0)
        // Optimal: x = z = 0
        let dd = DualDecomposition::new(0.1, 200, 1e-4);
        let x_solve = |_x: &[f64], y: &[f64]| -> Vec<f64> {
            // argmin_x x² + y*x = -y/2
            vec![-y[0] / 2.0]
        };
        let z_solve = |_z: &[f64], y: &[f64]| -> Vec<f64> {
            // argmin_z z² - y*z = y/2
            vec![y[0] / 2.0]
        };
        let mat_a = vec![vec![1.0_f64]];
        let mat_b = vec![vec![-1.0_f64]];
        let c = vec![0.0_f64];
        let (x, z, _y, _iters) = dd.solve(
            vec![1.0],
            vec![1.0],
            vec![0.0],
            &x_solve,
            &z_solve,
            &mat_a,
            &mat_b,
            &c,
        );
        // Best iterate should have small constraint violation |x - z|
        assert!(
            (x[0] - z[0]).abs() < 0.5,
            "|x - z| = {}",
            (x[0] - z[0]).abs()
        );
    }
}
