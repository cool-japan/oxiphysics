// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Proportional load control: incremental load stepping and the Crisfield
//! arc-length method.
//!
//! A structure under proportional loading is described by its internal force
//! vector `f_int(u)` and a fixed reference load `f_ref`; the equilibrium states
//! at load factor `λ` solve
//!
//! ```text
//! R(u, λ) = f_int(u) − λ f_ref = 0,        K(u) = ∂f_int/∂u.
//! ```
//!
//! [`LoadedResidualFn`] exposes exactly this split, which both solvers of this
//! module need: [`LoadStepper`] solves `R(u, λ_k) = 0` for a sequence of load
//! factors with Newton-Raphson, and [`ArcLengthMethod`] follows the equilibrium
//! path `(u(s), λ(s))` through limit points (snap-through) where load control
//! fails, using the spherical constraint
//!
//! ```text
//! Δu·Δu + ψ² Δλ² = Δs².
//! ```
//!
//! Each corrector iteration factors the tangent stiffness once and solves for
//! the two displacement increments `δu_t = K⁻¹ f_ref` (load direction) and
//! `δu_r = −K⁻¹ R` (residual correction); the load increment `δλ` follows from
//! the quadratic obtained by inserting `δu = δu_r + δλ δu_t` into the
//! constraint, and the root is selected with the angle criterion of Crisfield
//! (1981).

use super::{
    CoupledNonlinearProblem, NewtonRaphson, NonlinearSpring, ResidualFn, SimpleLinearProblem,
    norm_l2, vec_dot, vec_sub,
};

/// Number of times [`ArcLengthMethod::trace_path`] halves the arc length of a
/// step that fails before it gives up (the smallest retry uses `Δs / 2¹⁰`).
const MAX_STEP_CUTS: u32 = 10;

/// `f_int − λ f_ref`, the single definition of the residual at a load factor.
fn load_residual(internal: &[f64], reference: &[f64], load_factor: f64) -> Vec<f64> {
    internal
        .iter()
        .zip(reference.iter())
        .map(|(f_int, f_ref)| f_int - load_factor * f_ref)
        .collect()
}

// ─── Loaded residual trait ───────────────────────────────────────────────────

/// A nonlinear problem under proportional loading,
/// `R(u, λ) = f_int(u) − λ f_ref`.
///
/// This is a [`ResidualFn`] whose full-load residual (`λ = 1`) is
/// `residual(u) = internal_force(u) − reference_load()` and whose
/// `tangent_stiffness(u)` is `∂ internal_force / ∂u`.  Implementors must keep
/// these three methods consistent.  A load that does not scale with the load
/// factor (dead load, prestress) belongs into the internal force as a
/// constant offset.
pub trait LoadedResidualFn: ResidualFn {
    /// Internal force vector `f_int(u)`.
    fn internal_force(&self, u: &[f64]) -> Vec<f64>;

    /// Reference load vector `f_ref` (the external load at load factor 1).
    fn reference_load(&self) -> Vec<f64>;

    /// Residual `f_int(u) − λ f_ref` at the load factor `load_factor`.
    fn residual_at(&self, u: &[f64], load_factor: f64) -> Vec<f64> {
        load_residual(&self.internal_force(u), &self.reference_load(), load_factor)
    }
}

impl LoadedResidualFn for SimpleLinearProblem {
    fn internal_force(&self, u: &[f64]) -> Vec<f64> {
        vec![self.k * u.first().copied().unwrap_or(0.0)]
    }

    fn reference_load(&self) -> Vec<f64> {
        vec![self.f]
    }
}

impl LoadedResidualFn for NonlinearSpring {
    fn internal_force(&self, u: &[f64]) -> Vec<f64> {
        let x = u.first().copied().unwrap_or(0.0);
        vec![self.k1 * x + self.k3 * x.powi(3)]
    }

    fn reference_load(&self) -> Vec<f64> {
        vec![self.force]
    }
}

impl LoadedResidualFn for CoupledNonlinearProblem {
    fn internal_force(&self, u: &[f64]) -> Vec<f64> {
        let (x, y) = (
            u.first().copied().unwrap_or(0.0),
            u.get(1).copied().unwrap_or(0.0),
        );
        vec![x * x + y, x + y * y]
    }

    fn reference_load(&self) -> Vec<f64> {
        vec![3.0, 3.0]
    }
}

// ─── Dense LU ────────────────────────────────────────────────────────────────

/// LU factorisation with partial pivoting of a dense square matrix, so that the
/// tangent stiffness is factored once per corrector iteration and reused for
/// several right-hand sides.
struct DenseLu {
    n: usize,
    /// Row-major `L` (unit lower, below the diagonal) and `U` factors.
    lu: Vec<f64>,
    /// Row exchanged with row `k` at elimination step `k`.
    pivots: Vec<usize>,
}

impl DenseLu {
    /// Factor `a`.
    ///
    /// # Errors
    /// An error message if `a` is not square, contains non-finite entries or is
    /// singular (a pivot negligible against the largest entry).
    fn factor(a: &[Vec<f64>]) -> Result<Self, String> {
        let n = a.len();
        if a.iter().any(|row| row.len() != n) {
            return Err("Matrix dimensions inconsistent".to_string());
        }
        let mut lu: Vec<f64> = a.iter().flat_map(|row| row.iter().copied()).collect();
        if lu.iter().any(|v| !v.is_finite()) {
            return Err("Tangent stiffness contains non-finite entries".to_string());
        }
        let scale = lu.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let negligible = scale * n as f64 * f64::EPSILON;
        let mut pivots = vec![0_usize; n];
        for col in 0..n {
            let mut pivot_row = col;
            let mut best = lu[col * n + col].abs();
            for row in (col + 1)..n {
                let value = lu[row * n + col].abs();
                if value > best {
                    best = value;
                    pivot_row = row;
                }
            }
            if best <= negligible {
                return Err(format!("Singular matrix detected at column {col}"));
            }
            pivots[col] = pivot_row;
            if pivot_row != col {
                for j in 0..n {
                    lu.swap(col * n + j, pivot_row * n + j);
                }
            }
            let pivot = lu[col * n + col];
            for row in (col + 1)..n {
                let factor = lu[row * n + col] / pivot;
                lu[row * n + col] = factor;
                if factor != 0.0 {
                    for j in (col + 1)..n {
                        let update = factor * lu[col * n + j];
                        lu[row * n + j] -= update;
                    }
                }
            }
        }
        Ok(Self { n, lu, pivots })
    }

    /// Solve `A x = b`.
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x = b.to_vec();
        for (col, &pivot_row) in self.pivots.iter().enumerate() {
            x.swap(col, pivot_row);
        }
        for i in 1..n {
            let row = &self.lu[i * n..i * n + i];
            let sum: f64 = row.iter().zip(x[..i].iter()).map(|(l, v)| l * v).sum();
            x[i] -= sum;
        }
        for i in (0..n).rev() {
            let row = &self.lu[i * n + i + 1..(i + 1) * n];
            let sum: f64 = row.iter().zip(x[i + 1..].iter()).map(|(u, v)| u * v).sum();
            x[i] = (x[i] - sum) / self.lu[i * n + i];
        }
        x
    }
}

// ─── Load stepper ────────────────────────────────────────────────────────────

/// Incremental load stepping: solves `f_int(u) − λ f_ref = 0` at a sequence of
/// load factors using a Newton-Raphson solver at each step.
///
/// The solution of a step is the starting point of the next one.  Load control
/// cannot pass a limit point of the equilibrium path; use [`ArcLengthMethod`]
/// for that.
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

    /// Solve the equilibrium `f_int(u) − λ f_ref = 0` at every load factor of
    /// [`Self::step_factors`], starting from `u0`, and return the pairs
    /// `(load factor, displacement)`.
    ///
    /// # Errors
    /// An error message naming the failed step if the dimensions disagree, the
    /// Newton solver reports an error (for instance a singular tangent) or a
    /// step does not converge; no partial result is returned.
    pub fn solve<F: LoadedResidualFn>(
        &self,
        problem: &F,
        u0: &[f64],
        solver: &NewtonRaphson,
    ) -> Result<Vec<(f64, Vec<f64>)>, String> {
        let reference = problem.reference_load();
        if reference.len() != u0.len() {
            return Err(format!(
                "reference load has {} entries but the start vector has {}",
                reference.len(),
                u0.len()
            ));
        }
        let factors = self.step_factors();
        let mut results = Vec::with_capacity(factors.len());
        let mut u = u0.to_vec();
        for (step, &load_factor) in factors.iter().enumerate() {
            let scaled = ScaledProblem {
                inner: problem,
                load_factor,
                reference: &reference,
            };
            let outcome = solver.solve(&scaled, &u).map_err(|error| {
                format!(
                    "load step {} of {} (load factor {load_factor}): {error}",
                    step + 1,
                    factors.len()
                )
            })?;
            if !outcome.converged {
                return Err(format!(
                    "load step {} of {} (load factor {load_factor}) did not converge: residual \
                     {:.3e} after {} Newton iterations",
                    step + 1,
                    factors.len(),
                    outcome.final_residual,
                    outcome.iterations
                ));
            }
            u = outcome.solution.clone();
            results.push((load_factor, outcome.solution));
        }
        Ok(results)
    }
}

/// The residual `f_int(u) − λ f_ref` at a fixed load factor, as a
/// [`ResidualFn`] for the Newton solvers.  Only the external load is scaled;
/// the internal force and the tangent stiffness are those of the structure.
struct ScaledProblem<'a, F> {
    inner: &'a F,
    load_factor: f64,
    reference: &'a [f64],
}

impl<F: LoadedResidualFn> ResidualFn for ScaledProblem<'_, F> {
    fn residual(&self, u: &[f64]) -> Vec<f64> {
        load_residual(
            &self.inner.internal_force(u),
            self.reference,
            self.load_factor,
        )
    }

    fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
        self.inner.tangent_stiffness(u)
    }
}

// ─── Arc-length method (Crisfield) ───────────────────────────────────────────

/// Crisfield spherical arc-length method for tracing equilibrium paths through
/// limit points (snap-through and snap-back).
///
/// A step starts from an equilibrium state `(u₀, λ₀)`, predicts along the
/// tangent `(K⁻¹ f_ref, 1)` scaled to the arc length and then iterates on the
/// constraint `Δu·Δu + ψ² Δλ² = Δs²` together with equilibrium (see the module
/// documentation).
#[derive(Debug, Clone)]
pub struct ArcLengthMethod {
    /// Constraint radius `Δs` in the `(u, ψλ)` space.
    pub arc_length: f64,
    /// Maximum corrector iterations per step.
    pub max_iterations: usize,
    /// Convergence tolerance on the Euclidean norm of the residual
    /// `f_int(u) − λ f_ref` (absolute, in force units).
    pub tolerance: f64,
    /// Weight `ψ` of the load factor in the constraint metric (default 1);
    /// `ψ = 0` gives the cylindrical (displacement-only) constraint.  Choose it
    /// so that `ψ Δλ` is of the order of the displacement increments.
    pub load_scaling: f64,
}

impl ArcLengthMethod {
    /// Create a new arc-length solver with `load_scaling = 1`.
    pub fn new(arc_length: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            arc_length,
            max_iterations: max_iter,
            tolerance: tol,
            load_scaling: 1.0,
        }
    }

    /// Return the solver with the load-factor weight `ψ` of the constraint set.
    pub fn with_load_scaling(mut self, load_scaling: f64) -> Self {
        self.load_scaling = load_scaling;
        self
    }

    /// Perform one arc-length predictor-corrector step from the equilibrium
    /// state `(u0, lambda0)` and return the new `(u, λ)`.
    ///
    /// `previous_increment` is the displacement and load increment of the
    /// preceding step; the new step proceeds in the direction that has a
    /// positive inner product (in the metric of the constraint) with it, which
    /// keeps the path from turning back at a limit point.  Without history
    /// (`None`) the load is increased.
    ///
    /// # Errors
    /// An error message for inconsistent dimensions or parameters, a singular
    /// tangent stiffness, a corrector for which the arc-length constraint has no
    /// real solution (the arc length is too large), or when `max_iterations`
    /// are exhausted.
    pub fn solve_step<F: LoadedResidualFn>(
        &self,
        problem: &F,
        u0: &[f64],
        lambda0: f64,
        previous_increment: Option<(&[f64], f64)>,
    ) -> Result<(Vec<f64>, f64), String> {
        self.advance(problem, u0, lambda0, previous_increment, self.arc_length)
    }

    /// Trace the equilibrium path over `n_steps` arc-length steps from the
    /// equilibrium state `(u0, lambda0)`, returning the converged points after
    /// each step (the start is not included).
    ///
    /// A step that does not converge is retried with the arc length halved, up
    /// to ten times, before the tracing stops; the returned path is then
    /// shorter than `n_steps`.  The next step starts again with the nominal arc
    /// length.
    pub fn trace_path<F: LoadedResidualFn>(
        &self,
        problem: &F,
        u0: &[f64],
        lambda0: f64,
        n_steps: usize,
    ) -> Vec<(Vec<f64>, f64)> {
        let mut path = Vec::with_capacity(n_steps);
        let mut u = u0.to_vec();
        let mut lambda = lambda0;
        let mut previous: Option<(Vec<f64>, f64)> = None;
        for _ in 0..n_steps {
            let mut arc_length = self.arc_length;
            let mut accepted = None;
            for _ in 0..=MAX_STEP_CUTS {
                let history = previous.as_ref().map(|(du, dl)| (du.as_slice(), *dl));
                if let Ok(next) = self.advance(problem, &u, lambda, history, arc_length) {
                    accepted = Some(next);
                    break;
                }
                arc_length *= 0.5;
            }
            let Some((u_new, lambda_new)) = accepted else {
                break;
            };
            previous = Some((vec_sub(&u_new, &u), lambda_new - lambda));
            u = u_new.clone();
            lambda = lambda_new;
            path.push((u_new, lambda_new));
        }
        path
    }

    /// One predictor-corrector step with the constraint radius `arc_length`.
    fn advance<F: LoadedResidualFn>(
        &self,
        problem: &F,
        u0: &[f64],
        lambda0: f64,
        previous_increment: Option<(&[f64], f64)>,
        arc_length: f64,
    ) -> Result<(Vec<f64>, f64), String> {
        let n = u0.len();
        if n == 0 {
            return Err("Empty system".to_string());
        }
        if !(arc_length.is_finite() && arc_length > 0.0) {
            return Err(format!(
                "arc length must be finite and positive, got {arc_length}"
            ));
        }
        if !(self.tolerance.is_finite() && self.tolerance > 0.0) {
            return Err(format!(
                "tolerance must be finite and positive, got {}",
                self.tolerance
            ));
        }
        if !(self.load_scaling.is_finite() && self.load_scaling >= 0.0) {
            return Err(format!(
                "load scaling must be finite and non-negative, got {}",
                self.load_scaling
            ));
        }
        let reference = problem.reference_load();
        if reference.len() != n {
            return Err(format!(
                "reference load has {} entries but the displacement vector has {n}",
                reference.len()
            ));
        }
        if let Some((du_previous, _)) = previous_increment
            && du_previous.len() != n
        {
            return Err(format!(
                "previous increment has {} entries but the displacement vector has {n}",
                du_previous.len()
            ));
        }
        let psi_sq = self.load_scaling * self.load_scaling;
        let step_sq = arc_length * arc_length;

        // Predictor along the tangent (K⁻¹ f_ref, 1).
        let tangent_factor = factor_tangent(problem, u0, n)?;
        let du_t = tangent_factor.solve(&reference);
        let norm = (vec_dot(&du_t, &du_t) + psi_sq).sqrt();
        if !(norm.is_finite() && norm > 0.0) {
            return Err("the reference load produces no tangent displacement".to_string());
        }
        let forward = match previous_increment {
            Some((du_previous, dlambda_previous)) => {
                vec_dot(&du_t, du_previous) + psi_sq * dlambda_previous >= 0.0
            }
            None => true,
        };
        let sign = if forward { 1.0 } else { -1.0 };
        let mut d_lambda = sign * arc_length / norm;
        let mut du: Vec<f64> = du_t.iter().map(|v| d_lambda * v).collect();
        let mut u: Vec<f64> = u0.iter().zip(du.iter()).map(|(a, b)| a + b).collect();
        let mut lambda = lambda0 + d_lambda;

        for iteration in 0..=self.max_iterations {
            let internal = problem.internal_force(&u);
            if internal.len() != n {
                return Err(format!(
                    "internal force has {} entries but the displacement vector has {n}",
                    internal.len()
                ));
            }
            let residual = load_residual(&internal, &reference, lambda);
            let residual_norm = norm_l2(&residual);
            if !residual_norm.is_finite() {
                return Err("Arc-length corrector produced a non-finite residual".to_string());
            }
            if residual_norm < self.tolerance {
                return Ok((u, lambda));
            }
            if iteration == self.max_iterations {
                break;
            }
            // One factorisation serves both right-hand sides.
            let tangent_factor = factor_tangent(problem, &u, n)?;
            let du_t = tangent_factor.solve(&reference);
            let du_r: Vec<f64> = tangent_factor.solve(&residual).iter().map(|v| -v).collect();

            // |Δu + δu_r + δλ δu_t|² + ψ² (Δλ + δλ)² = Δs²
            let shifted: Vec<f64> = du.iter().zip(du_r.iter()).map(|(a, b)| a + b).collect();
            let a1 = vec_dot(&du_t, &du_t) + psi_sq;
            let a2 = 2.0 * (vec_dot(&du_t, &shifted) + psi_sq * d_lambda);
            let a3 = vec_dot(&shifted, &shifted) + psi_sq * d_lambda * d_lambda - step_sq;
            let roots = solve_quadratic(a1, a2, a3)?;

            // Angle criterion: the root whose new increment makes the smaller
            // angle with the current one, i.e. the larger inner product.
            let alignment = |root: f64| {
                vec_dot(&shifted, &du)
                    + root * vec_dot(&du_t, &du)
                    + psi_sq * (d_lambda + root) * d_lambda
            };
            let delta_lambda = if alignment(roots[0]) >= alignment(roots[1]) {
                roots[0]
            } else {
                roots[1]
            };

            for i in 0..n {
                let delta_u = du_r[i] + delta_lambda * du_t[i];
                du[i] += delta_u;
                u[i] += delta_u;
            }
            d_lambda += delta_lambda;
            lambda += delta_lambda;
        }
        Err(format!(
            "Arc-length method did not converge in {} iterations",
            self.max_iterations
        ))
    }
}

/// Assemble the tangent stiffness at `u` and factor it.
fn factor_tangent<F: LoadedResidualFn>(
    problem: &F,
    u: &[f64],
    n: usize,
) -> Result<DenseLu, String> {
    let tangent = problem.tangent_stiffness(u);
    if tangent.len() != n {
        return Err(format!(
            "tangent stiffness has {} rows but the displacement vector has {n}",
            tangent.len()
        ));
    }
    DenseLu::factor(&tangent)
}

/// Real roots of `a x² + b x + c = 0` with `a > 0`, computed without
/// cancellation.
///
/// A discriminant that is negative only by rounding is treated as zero.
///
/// # Errors
/// An error message if `a` is not positive or the roots are complex, i.e. the
/// constraint sphere does not meet the linearised equilibrium line.
fn solve_quadratic(a: f64, b: f64, c: f64) -> Result<[f64; 2], String> {
    if !(a.is_finite() && a > 0.0 && b.is_finite() && c.is_finite()) {
        return Err("Arc-length quadratic is degenerate".to_string());
    }
    let mut discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        if discriminant >= -1.0e-10 * (b * b + (4.0 * a * c).abs()) {
            discriminant = 0.0;
        } else {
            return Err(
                "Arc-length constraint has no real solution; reduce the arc length".to_string(),
            );
        }
    }
    let q = -0.5 * (b + b.signum() * discriminant.sqrt());
    let first = q / a;
    let second = if q != 0.0 { c / q } else { first };
    Ok([first, second])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cubic snap-through spring `f_int(u) = u³ − 6u² + 9u = u (u − 3)²`.
    ///
    /// Load control from the origin reaches the limit point `(u, λ) = (1, 4)`
    /// and cannot continue; the path then returns through `(3, 0)` (stiffness
    /// `3u² − 12u + 9 < 0` between `u = 1` and `u = 3`) and rises again.
    struct CubicSnap;

    impl ResidualFn for CubicSnap {
        fn residual(&self, u: &[f64]) -> Vec<f64> {
            self.internal_force(u)
                .iter()
                .zip(self.reference_load().iter())
                .map(|(f, r)| f - r)
                .collect()
        }

        fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
            vec![vec![3.0 * u[0] * u[0] - 12.0 * u[0] + 9.0]]
        }
    }

    impl LoadedResidualFn for CubicSnap {
        fn internal_force(&self, u: &[f64]) -> Vec<f64> {
            vec![u[0] * (u[0] - 3.0) * (u[0] - 3.0)]
        }

        fn reference_load(&self) -> Vec<f64> {
            vec![1.0]
        }
    }

    /// Two-bar von Mises truss: bars from the supports `(±half_span, 0)` to the
    /// apex at height `rise`; `u = w` is the downward apex displacement and the
    /// load is a downward force at the apex.  Geometrically exact, linear elastic
    /// bars:  `f_int(w) = 2 EA (l₀ − l)(h − w) / (l₀ l)`,
    /// `l = √(L² + (h − w)²)`.
    struct VonMisesTruss {
        ea: f64,
        half_span: f64,
        rise: f64,
    }

    impl VonMisesTruss {
        fn initial_length(&self) -> f64 {
            self.half_span.hypot(self.rise)
        }

        fn force(&self, w: f64) -> f64 {
            let s = self.rise - w;
            let l = self.half_span.hypot(s);
            let l0 = self.initial_length();
            2.0 * self.ea * (l0 - l) * s / (l0 * l)
        }

        /// Stiffness `d f_int / d w = 2 EA/l₀ (1 − l₀ L² / l³)`.
        fn stiffness(&self, w: f64) -> f64 {
            let l = self.half_span.hypot(self.rise - w);
            let l0 = self.initial_length();
            2.0 * self.ea / l0 * (1.0 - l0 * self.half_span * self.half_span / l.powi(3))
        }

        /// Limit load `max_{0<w<h} f_int(w)`, from the stationary point
        /// `l³ = l₀ L²` of the force.
        fn limit_load(&self) -> f64 {
            let l0 = self.initial_length();
            let l_star = (l0 * self.half_span * self.half_span).cbrt();
            let s_star = (l_star * l_star - self.half_span * self.half_span).sqrt();
            self.force(self.rise - s_star)
        }
    }

    impl ResidualFn for VonMisesTruss {
        fn residual(&self, u: &[f64]) -> Vec<f64> {
            vec![self.force(u[0]) - 1.0]
        }

        fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
            vec![vec![self.stiffness(u[0])]]
        }
    }

    impl LoadedResidualFn for VonMisesTruss {
        fn internal_force(&self, u: &[f64]) -> Vec<f64> {
            vec![self.force(u[0])]
        }

        fn reference_load(&self) -> Vec<f64> {
            vec![1.0]
        }
    }

    /// Cubic snap-through spring coupled to a linear spring and loaded at both
    /// DOFs: `f_int = [g(u₀) + k (u₀ − u₁), k (u₁ − u₀) + c u₁]`.
    struct CoupledSnap {
        k: f64,
        c: f64,
    }

    impl ResidualFn for CoupledSnap {
        fn residual(&self, u: &[f64]) -> Vec<f64> {
            self.internal_force(u)
                .iter()
                .zip(self.reference_load().iter())
                .map(|(f, r)| f - r)
                .collect()
        }

        fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
            vec![
                vec![3.0 * u[0] * u[0] - 12.0 * u[0] + 9.0 + self.k, -self.k],
                vec![-self.k, self.k + self.c],
            ]
        }
    }

    impl LoadedResidualFn for CoupledSnap {
        fn internal_force(&self, u: &[f64]) -> Vec<f64> {
            vec![
                u[0] * (u[0] - 3.0) * (u[0] - 3.0) + self.k * (u[0] - u[1]),
                self.k * (u[1] - u[0]) + self.c * u[1],
            ]
        }

        fn reference_load(&self) -> Vec<f64> {
            vec![1.0, 0.5]
        }
    }

    /// Max-norm of `f_int(u) − λ f_ref` along a path.
    fn equilibrium_error<F: LoadedResidualFn>(problem: &F, path: &[(Vec<f64>, f64)]) -> f64 {
        path.iter()
            .map(|(u, lambda)| {
                problem
                    .residual_at(u, *lambda)
                    .iter()
                    .fold(0.0_f64, |m, r| m.max(r.abs()))
            })
            .fold(0.0, f64::max)
    }

    // ---- dense LU ------------------------------------------------------------

    #[test]
    fn lu_solves_systems_that_need_pivoting_for_several_right_hand_sides() {
        // The first pivot is zero, so a row exchange is required.
        let a = vec![
            vec![0.0, 2.0, 1.0],
            vec![1.0, 1.0, 1.0],
            vec![3.0, 0.0, 2.0],
        ];
        let lu = DenseLu::factor(&a).expect("regular matrix");
        for x_true in [[1.0, -2.0, 3.0], [0.5, 0.25, -4.0]] {
            let b: Vec<f64> = a
                .iter()
                .map(|row| row.iter().zip(x_true.iter()).map(|(m, x)| m * x).sum())
                .collect();
            let x = lu.solve(&b);
            for (xi, ti) in x.iter().zip(x_true.iter()) {
                assert!((xi - ti).abs() < 1e-13, "{xi} vs {ti}");
            }
        }
    }

    #[test]
    fn lu_rejects_singular_and_malformed_matrices() {
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(DenseLu::factor(&singular).is_err());
        assert!(DenseLu::factor(&[vec![0.0]]).is_err());
        assert!(DenseLu::factor(&[vec![1.0, 2.0]]).is_err());
        assert!(DenseLu::factor(&[vec![f64::NAN]]).is_err());
        // Scale invariance: a tiny but regular matrix is accepted.
        assert!(DenseLu::factor(&[vec![1e-30, 0.0], vec![0.0, 2e-30]]).is_ok());
    }

    #[test]
    fn quadratic_roots_are_stable_and_ordered_by_construction() {
        let [r1, r2] = solve_quadratic(1.0, -3.0, 2.0).expect("real roots");
        let (lo, hi) = (r1.min(r2), r1.max(r2));
        assert!((lo - 1.0).abs() < 1e-14 && (hi - 2.0).abs() < 1e-14);
        // Large b: no cancellation in the small root.
        let [s1, s2] = solve_quadratic(1.0, 1e8, 1.0).expect("real roots");
        let small = if s1.abs() < s2.abs() { s1 } else { s2 };
        assert!((small + 1e-8).abs() < 1e-20, "{small}");
        assert!(solve_quadratic(1.0, 0.0, 1.0).is_err());
        assert!(solve_quadratic(0.0, 1.0, 1.0).is_err());
        // Tangent circle: a discriminant of rounding size is accepted.
        let [t1, t2] = solve_quadratic(1.0, -2.0, 1.0 + 1e-17).expect("double root");
        assert!((t1 - 1.0).abs() < 1e-7 && (t2 - 1.0).abs() < 1e-7);
    }

    // ---- load stepper --------------------------------------------------------

    #[test]
    fn load_stepper_solves_the_scaled_equilibrium_not_the_full_load_problem() {
        // f_int = u + u³, f_ref = 1:  u(0.5) = 0.4239 (u + u³ = 0.5), u(1) = 0.6823.
        // The former implementation scaled the whole residual and solved the
        // full-load problem at every step.
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 1.0,
            force: 1.0,
        };
        let stepper = LoadStepper::new(3, 0.0, 1.0);
        let solver = NewtonRaphson::new(50, 1e-12);
        let results = stepper.solve(&problem, &[0.0], &solver).expect("converges");
        assert_eq!(results.len(), 3);
        assert!((results[0].0 - 0.0).abs() < 1e-15 && results[0].1[0].abs() < 1e-12);
        assert!((results[1].0 - 0.5).abs() < 1e-15);
        assert!(
            (results[1].1[0] - 0.4239).abs() < 5e-5,
            "u(0.5) = {}",
            results[1].1[0]
        );
        assert!(
            (results[2].1[0] - 0.682_327_803_8).abs() < 1e-8,
            "u(1) = {}",
            results[2].1[0]
        );
        for (load_factor, u) in &results {
            let residual = problem.residual_at(u, *load_factor);
            assert!(residual[0].abs() < 1e-10);
        }
    }

    #[test]
    fn load_stepper_follows_a_linear_problem_proportionally() {
        let problem = SimpleLinearProblem { k: 2.0, f: 5.0 };
        let stepper = LoadStepper::new(6, 0.0, 1.0);
        let solver = NewtonRaphson::new(20, 1e-12);
        let results = stepper.solve(&problem, &[0.0], &solver).expect("converges");
        assert_eq!(results.len(), 6);
        for (load_factor, u) in &results {
            assert!((u[0] - 2.5 * load_factor).abs() < 1e-12, "{load_factor}");
        }
    }

    #[test]
    fn load_stepper_handles_coupled_systems() {
        // f_int = (x² + y, x + y²), f_ref = (3, 3):  the symmetric branch has
        // x = y with x² + x = 3λ, i.e. x = (√13 - 1)/2 at λ = 1.
        let problem = CoupledNonlinearProblem;
        let stepper = LoadStepper::new(5, 0.2, 1.0);
        let solver = NewtonRaphson::new(50, 1e-12);
        let results = stepper
            .solve(&problem, &[0.1, 0.1], &solver)
            .expect("converges");
        let (load_factor, u) = results.last().expect("steps");
        assert!((load_factor - 1.0).abs() < 1e-15);
        let expected = (13.0_f64.sqrt() - 1.0) / 2.0;
        assert!((u[0] - expected).abs() < 1e-9 && (u[1] - expected).abs() < 1e-9);
        for (factor, state) in &results {
            let residual = problem.residual_at(state, *factor);
            assert!(residual.iter().all(|r| r.abs() < 1e-10));
        }
    }

    #[test]
    fn load_stepper_reports_a_failed_step_instead_of_returning_unconverged_states() {
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 1.0,
            force: 10.0,
        };
        let stepper = LoadStepper::new(4, 0.25, 1.0);
        // One Newton iteration cannot reach the tolerance.
        let starved = NewtonRaphson::new(1, 1e-14);
        let error = stepper
            .solve(&problem, &[0.0], &starved)
            .expect_err("must not converge");
        assert!(error.contains("did not converge"), "{error}");
        assert!(error.contains("load step 1 of 4"), "{error}");
        // A singular tangent is propagated with the step it occurred in.
        let singular = NonlinearSpring {
            k1: 0.0,
            k3: 1.0,
            force: 1.0,
        };
        let solver = NewtonRaphson::new(20, 1e-12);
        let error = LoadStepper::new(2, 0.5, 1.0)
            .solve(&singular, &[0.0], &solver)
            .expect_err("singular");
        assert!(error.contains("Singular"), "{error}");
        // Mismatched dimensions.
        assert!(stepper.solve(&problem, &[0.0, 0.0], &solver).is_err());
    }

    #[test]
    fn load_stepper_with_no_steps_returns_no_results() {
        let problem = SimpleLinearProblem { k: 1.0, f: 1.0 };
        let solver = NewtonRaphson::new(5, 1e-10);
        let results = LoadStepper::new(0, 0.0, 1.0)
            .solve(&problem, &[0.0], &solver)
            .expect("nothing to do");
        assert!(results.is_empty());
    }

    // ---- arc length ----------------------------------------------------------

    #[test]
    fn linear_step_lies_on_the_constraint_sphere_in_both_directions() {
        // K = 2, f_ref = 1:  u = λ/2.  Δs = 1, ψ = 1:  Δλ = Δs / √(1/4 + 1).
        let problem = SimpleLinearProblem { k: 2.0, f: 1.0 };
        let arc = ArcLengthMethod::new(1.0, 20, 1e-12);
        let expected_dl = 1.0 / (0.25_f64 + 1.0).sqrt();
        let (u, lambda) = arc
            .solve_step(&problem, &[0.0], 0.0, None)
            .expect("forward step");
        assert!((lambda - expected_dl).abs() < 1e-13, "λ = {lambda}");
        assert!((u[0] - 0.5 * expected_dl).abs() < 1e-13);
        assert!((u[0] * u[0] + lambda * lambda - 1.0).abs() < 1e-12);
        // A history pointing backwards in load reverses the step.
        let (u_back, lambda_back) = arc
            .solve_step(&problem, &[0.0], 0.0, Some((&[-0.1], -0.2)))
            .expect("backward step");
        assert!((lambda_back + expected_dl).abs() < 1e-13);
        assert!((u_back[0] + 0.5 * expected_dl).abs() < 1e-13);
    }

    #[test]
    fn corrector_uses_a_genuine_load_direction() {
        // With the former du_r = -du_t the corrector had no load direction; here the
        // converged point must be an equilibrium state of the nonlinear spring,
        // f_int(u) = lambda f_ref, at distance Δs from the start.
        let problem = NonlinearSpring {
            k1: 1.0,
            k3: 0.5,
            force: 1.0,
        };
        let arc = ArcLengthMethod::new(0.3, 30, 1e-12);
        let (u, lambda) = arc
            .solve_step(&problem, &[0.0], 0.0, None)
            .expect("step converges");
        let f_int = problem.internal_force(&u)[0];
        assert!((f_int - lambda).abs() < 1e-11, "f_int {f_int}, λ {lambda}");
        assert!(lambda > 0.0 && u[0] > 0.0);
        assert!((u[0] * u[0] + lambda * lambda - 0.09).abs() < 1e-10);
        // A hardening spring is stiffer than the tangent: u < λ.
        assert!(u[0] < lambda);
    }

    #[test]
    fn cubic_snap_through_is_traced_past_the_limit_point() {
        let problem = CubicSnap;
        let arc = ArcLengthMethod::new(0.25, 40, 1e-10);
        let path = arc.trace_path(&problem, &[0.0], 0.0, 70);
        assert_eq!(path.len(), 70, "tracing stopped early");
        assert!(equilibrium_error(&problem, &path) < 1e-9);

        // The displacement increases monotonically along the whole path ...
        for pair in path.windows(2) {
            assert!(
                pair[1].0[0] > pair[0].0[0],
                "{:?} -> {:?}",
                pair[0],
                pair[1]
            );
        }
        // ... while the load first rises to the limit load 4 at u = 1 ...
        let peak_index = path
            .windows(2)
            .position(|p| p[1].1 < p[0].1)
            .expect("the load decreases after the limit point");
        let (peak_u, peak) = (path[peak_index].0[0], path[peak_index].1);
        assert!((peak - 4.0).abs() < 0.1, "limit load {peak}");
        assert!((peak_u - 1.0).abs() < 0.3, "limit displacement {peak_u}");
        // ... falls monotonically along the unstable branch to the second
        // equilibrium (u, lambda) = (3, 0) ...
        let trough_index = (peak_index..path.len())
            .min_by(|&a, &b| path[a].1.total_cmp(&path[b].1))
            .expect("non-empty");
        assert!(
            path[peak_index..=trough_index]
                .windows(2)
                .all(|p| p[1].1 < p[0].1)
        );
        let (trough_u, trough) = (path[trough_index].0[0], path[trough_index].1);
        assert!(
            (trough_u - 3.0).abs() < 0.3,
            "second equilibrium at u = {trough_u}"
        );
        assert!((0.0..0.15).contains(&trough), "trough load {trough}");
        // ... and rises again on the stiff third branch.
        let last = path.last().expect("non-empty");
        assert!(last.0[0] > 3.5 && last.1 > 1.0, "path ended at {last:?}");
    }

    #[test]
    fn von_mises_truss_snaps_through_at_the_analytic_limit_load() {
        let truss = VonMisesTruss {
            ea: 100.0,
            half_span: 2.0,
            rise: 1.0,
        };
        let limit = truss.limit_load();
        assert!(limit > 3.0 && limit < 5.0, "limit load {limit}");
        let arc = ArcLengthMethod::new(0.1, 40, 1e-9);
        let path = arc.trace_path(&truss, &[0.0], 0.0, 170);
        assert_eq!(path.len(), 170);
        assert!(equilibrium_error(&truss, &path) < 1e-8);

        // The first local maximum of the load is the analytic limit load.
        let peak_index = path
            .windows(2)
            .position(|p| p[1].1 < p[0].1)
            .expect("load decreases after the limit point");
        let peak = path[peak_index].1;
        assert!(
            (peak - limit).abs() < 0.02 * limit,
            "peak {peak}, analytic {limit}"
        );
        assert!(
            path[peak_index].0[0] < 0.6,
            "limit state at w = {}",
            path[peak_index].0[0]
        );
        // The apex then passes the flat configuration w = h = 1 with zero load ...
        let crossing = path
            .windows(2)
            .find(|p| p[0].1 > 0.0 && p[1].1 <= 0.0)
            .expect("load changes sign");
        let (w0, w1) = (crossing[0].0[0], crossing[1].0[0]);
        assert!(
            w0 < 1.0 + 1e-9 && w1 > 1.0 - 1e-9,
            "sign change in [{w0}, {w1}]"
        );
        // ... reaches the mirrored limit load -limit (the structure is symmetric
        // about w = h) ...
        let trough = path.iter().map(|(_, l)| *l).fold(f64::MAX, f64::min);
        assert!(
            (trough + limit).abs() < 0.02 * limit,
            "trough {trough}, analytic {}",
            -limit
        );
        // ... and returns to zero load at the mirrored configuration w = 2h.
        let second_crossing = path
            .windows(2)
            .skip(peak_index + 20)
            .find(|p| p[0].1 < 0.0 && p[1].1 >= 0.0)
            .expect("load changes sign again");
        assert!(
            (second_crossing[0].0[0] - 2.0).abs() < 0.1,
            "mirrored state near w = {}",
            second_crossing[0].0[0]
        );
        assert!(path.last().expect("non-empty").0[0] > 2.0);
    }

    /// Two-DOF system whose equilibrium path is a circle in the displacement
    /// plane, `(u₀ - 1)² + u₁² = 1`, traversed with the load
    /// `λ = u₀ + u₁/2`: `u₀` rises to 2 and then **decreases** again (snap-back)
    /// while the load passes two limit points.
    struct CircleLoop;

    impl ResidualFn for CircleLoop {
        fn residual(&self, u: &[f64]) -> Vec<f64> {
            self.internal_force(u)
                .iter()
                .zip(self.reference_load().iter())
                .map(|(f, r)| f - r)
                .collect()
        }

        fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
            vec![vec![1.0, 0.5], vec![2.0 * (u[0] - 1.0), 2.0 * u[1]]]
        }
    }

    impl LoadedResidualFn for CircleLoop {
        fn internal_force(&self, u: &[f64]) -> Vec<f64> {
            vec![
                u[0] + 0.5 * u[1],
                (u[0] - 1.0) * (u[0] - 1.0) + u[1] * u[1] - 1.0,
            ]
        }

        fn reference_load(&self) -> Vec<f64> {
            vec![1.0, 0.0]
        }
    }

    #[test]
    fn snap_back_path_turns_back_in_displacement() {
        let problem = CircleLoop;
        let arc = ArcLengthMethod::new(0.1, 40, 1e-10);
        let path = arc.trace_path(&problem, &[0.0, 0.0], 0.0, 70);
        assert_eq!(path.len(), 70, "tracing stopped early");
        assert!(equilibrium_error(&problem, &path) < 1e-9);
        // Every point lies on the circle and carries the load u0 + u1 / 2.
        for (u, lambda) in &path {
            let radius = ((u[0] - 1.0).powi(2) + u[1] * u[1]).sqrt();
            assert!((radius - 1.0).abs() < 1e-9);
            assert!((lambda - (u[0] + 0.5 * u[1])).abs() < 1e-9);
        }
        // u0 grows to the rightmost point of the circle (2, 0) ...
        let (apex_index, apex) = path
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.0[0].total_cmp(&b.1.0[0]))
            .map(|(i, (u, _))| (i, u[0]))
            .expect("non-empty");
        assert!(apex > 1.98 && apex <= 2.0 + 1e-9, "largest u0 = {apex}");
        // ... and then decreases again: snap-back.
        assert!(apex_index > 5 && apex_index + 10 < path.len());
        assert!(path[apex_index + 10].0[0] < apex - 0.2);
        assert!(
            path.last().expect("non-empty").0[1] < 0.0,
            "second half of the circle"
        );
        // The load passes its first limit point on the way (a local maximum).
        assert!(path.windows(2).any(|p| p[1].1 < p[0].1));
    }

    #[test]
    fn coupled_snap_through_is_traced_with_vector_unknowns() {
        let problem = CoupledSnap { k: 2.0, c: 3.0 };
        let arc = ArcLengthMethod::new(0.2, 40, 1e-10);
        let path = arc.trace_path(&problem, &[0.0, 0.0], 0.0, 60);
        assert_eq!(path.len(), 60);
        assert!(equilibrium_error(&problem, &path) < 1e-8);
        // The load rises to a limit point, falls while the first DOF keeps growing,
        // and rises again after the second equilibrium.
        let peak_index = path
            .windows(2)
            .position(|p| p[1].1 < p[0].1)
            .expect("load decreases after the limit point");
        assert!(peak_index > 5, "limit point after {peak_index} steps");
        let trough_index = (peak_index..path.len())
            .min_by(|&a, &b| path[a].1.total_cmp(&path[b].1))
            .expect("non-empty");
        assert!(trough_index > peak_index + 3 && trough_index + 3 < path.len());
        assert!(path[trough_index].1 < path[peak_index].1 - 0.5);
        assert!(path[trough_index + 3].1 > path[trough_index].1);
        for pair in path[peak_index..=trough_index].windows(2) {
            assert!(pair[1].0[0] > pair[0].0[0], "u_0 must keep growing");
        }
        // Both unknowns follow the loading: the second DOF is driven by the first.
        assert!(path[trough_index].0[1] > path[peak_index].0[1]);
    }

    /// Hardening spring `f_int = u + u³` whose tangent is declared singular for
    /// displacements in a window, so that a large step landing in the window fails
    /// and has to be cut back.
    struct WindowedSpring {
        lower: f64,
        upper: f64,
    }

    impl ResidualFn for WindowedSpring {
        fn residual(&self, u: &[f64]) -> Vec<f64> {
            vec![u[0] + u[0].powi(3) - 1.0]
        }

        fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
            if u[0] > self.lower && u[0] < self.upper {
                vec![vec![0.0]]
            } else {
                vec![vec![1.0 + 3.0 * u[0] * u[0]]]
            }
        }
    }

    impl LoadedResidualFn for WindowedSpring {
        fn internal_force(&self, u: &[f64]) -> Vec<f64> {
            vec![u[0] + u[0].powi(3)]
        }

        fn reference_load(&self) -> Vec<f64> {
            vec![1.0]
        }
    }

    #[test]
    fn failed_steps_are_retried_with_a_smaller_arc_length() {
        let problem = WindowedSpring {
            lower: 0.4,
            upper: 0.6,
        };
        let arc = ArcLengthMethod::new(0.7, 30, 1e-12);
        // The nominal step lands at u ≈ 0.49, inside the singular window.
        let error = arc
            .solve_step(&problem, &[0.0], 0.0, None)
            .expect_err("tangent is singular in the window");
        assert!(error.contains("Singular"), "{error}");
        // The path tracer halves the arc length once (0.35 puts the corrector
        // outside the window) and then continues with the nominal length.
        let path = arc.trace_path(&problem, &[0.0], 0.0, 3);
        assert_eq!(path.len(), 3);
        assert!(equilibrium_error(&problem, &path) < 1e-10);
        let (u1, l1) = (path[0].0[0], path[0].1);
        assert!(
            ((u1 * u1 + l1 * l1).sqrt() - 0.35).abs() < 1e-9,
            "first step length {}",
            (u1 * u1 + l1 * l1).sqrt()
        );
        let (du, dl) = (path[1].0[0] - u1, path[1].1 - l1);
        assert!(((du * du + dl * dl).sqrt() - 0.7).abs() < 1e-9);
    }

    #[test]
    fn invalid_arguments_and_singular_tangents_are_reported() {
        let arc = ArcLengthMethod::new(0.1, 20, 1e-10);
        // Zero tangent stiffness at the origin of f_int = u³.
        let flat = NonlinearSpring {
            k1: 0.0,
            k3: 1.0,
            force: 1.0,
        };
        let error = arc
            .solve_step(&flat, &[0.0], 0.0, None)
            .expect_err("singular tangent");
        assert!(error.contains("Singular"), "{error}");
        let problem = CubicSnap;
        assert!(arc.solve_step(&problem, &[], 0.0, None).is_err());
        assert!(arc.solve_step(&problem, &[0.0, 0.0], 0.0, None).is_err());
        assert!(
            arc.solve_step(&problem, &[0.0], 0.0, Some((&[1.0, 2.0], 1.0)))
                .is_err()
        );
        assert!(
            ArcLengthMethod::new(0.0, 20, 1e-10)
                .solve_step(&problem, &[0.0], 0.0, None)
                .is_err()
        );
        assert!(
            ArcLengthMethod::new(0.1, 20, 0.0)
                .solve_step(&problem, &[0.0], 0.0, None)
                .is_err()
        );
        assert!(
            ArcLengthMethod::new(0.1, 20, 1e-10)
                .with_load_scaling(f64::NAN)
                .solve_step(&problem, &[0.0], 0.0, None)
                .is_err()
        );
        // No corrector iterations: the nonlinear predictor is not in equilibrium.
        let error = ArcLengthMethod::new(0.5, 0, 1e-12)
            .solve_step(&problem, &[0.0], 0.0, None)
            .expect_err("not converged");
        assert!(error.contains("did not converge"), "{error}");
    }

    #[test]
    fn load_scaling_changes_the_step_but_not_the_path() {
        // ψ = 0 is the cylindrical constraint |Δu| = Δs; the equilibrium path is the
        // same, only the metric of the step differs.  (Δs = 0.17 does not hit the
        // limit point u = 1 exactly, where the tangent is singular.)
        let problem = CubicSnap;
        let cylindrical = ArcLengthMethod::new(0.17, 40, 1e-10).with_load_scaling(0.0);
        let path = cylindrical.trace_path(&problem, &[0.0], 0.0, 30);
        assert_eq!(path.len(), 30);
        assert!(equilibrium_error(&problem, &path) < 1e-9);
        let mut previous = 0.0;
        for (u, _) in &path {
            assert!(
                (u[0] - previous - 0.17).abs() < 1e-9,
                "Δu = {}",
                u[0] - previous
            );
            previous = u[0];
        }
        // The same landscape with ψ = 1 puts the steps on the (u, λ) sphere instead.
        let spherical = ArcLengthMethod::new(0.17, 40, 1e-10);
        let path = spherical.trace_path(&problem, &[0.0], 0.0, 30);
        assert_eq!(path.len(), 30);
        let (u1, l1) = (path[0].0[0], path[0].1);
        assert!(((u1 * u1 + l1 * l1).sqrt() - 0.17).abs() < 1e-9);
    }
}
