// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! General statistical manifolds defined by a parametric log-density and a
//! quadrature rule over the sample space.
//!
//! All geometric quantities are computed from their defining expectations:
//!
//! * Fisher metric `g_ij(θ) = E_θ[∂_iℓ ∂_jℓ]`,
//! * Amari–Chentsov skewness tensor `T_ijk(θ) = E_θ[∂_iℓ ∂_jℓ ∂_kℓ]`,
//! * α-connection (first kind)
//!   `Γ^(α)_{ij,k}(θ) = E_θ[(∂_i∂_jℓ + (1−α)/2 ∂_iℓ ∂_jℓ) ∂_kℓ]`,
//!
//! where `ℓ(x; θ) = ln p(x; θ)` and `E_θ[f] ≈ Σ_n w_n p(x_n; θ) f(x_n)` uses the
//! user-supplied quadrature nodes `(x_n, w_n)`.  Parameter derivatives of `ℓ`
//! are taken by central finite differences, so only the log-density itself has
//! to be provided.  The outer-product form of the metric keeps it positive
//! semi-definite regardless of discretisation error.

use std::fmt;
use std::sync::Arc;

use super::invert_matrix;
use crate::error::{PhysicsError, PhysicsResult};

/// Parametric log-density `ℓ(x; θ)`: the first argument is the sample point
/// `x`, the second the parameter vector `θ`.
///
/// The density must be normalised with respect to the quadrature measure of the
/// owning [`StatisticalManifold`].  Return `f64::NEG_INFINITY` outside the
/// support and `f64::NAN` for parameters outside the parameter domain.
pub type LogDensityFn = Arc<dyn Fn(&[f64], &[f64]) -> f64 + Send + Sync>;

/// Relative step for first parameter derivatives (≈ ε^{1/3}).
const SCORE_STEP_SCALE: f64 = 6.0e-6;
/// Relative step for second parameter derivatives (≈ ε^{1/4}).
const HESSIAN_STEP_SCALE: f64 = 1.2e-4;
/// Maximum tolerated deviation of the total quadrature mass from one.
const MASS_TOLERANCE: f64 = 1.0e-6;
/// Number of RK4 steps used to integrate the geodesic equation.
const GEODESIC_RK4_STEPS: usize = 64;
/// Maximum number of Newton iterations of the geodesic shooting method.
const SHOOTING_MAX_ITER: usize = 60;
/// Relative tolerance of the geodesic shooting residual `‖exp_p(v) − q‖∞`.
const SHOOTING_TOL: f64 = 1.0e-7;

/// Finite-difference step for parameter component `value`.
fn fd_step(scale: f64, value: f64) -> f64 {
    scale * value.abs().max(1.0)
}

/// Reusable per-node derivative buffers of the log-density at fixed `θ`.
struct NodeJet {
    /// Score `∂_iℓ(x_n; θ)`.
    score: Vec<f64>,
    /// Row-major Hessian `∂_i∂_jℓ(x_n; θ)` (only filled when requested).
    hessian: Vec<f64>,
}

/// A statistical manifold `{p(x; θ) : θ ∈ Θ ⊂ ℝ^dim}`.
///
/// The manifold is described by its log-density `ℓ(x; θ)` and by quadrature
/// nodes `(x_n, w_n)` approximating integrals over the sample space (for a
/// discrete sample space use the support points with unit weights).  It
/// provides the Fisher information metric, the Amari–Chentsov tensor, the
/// α-connections (Levi-Civita for `α = 0`), the exponential map and geodesics.
#[derive(Clone)]
pub struct StatisticalManifold {
    dim: usize,
    log_density: LogDensityFn,
    nodes: Vec<(Vec<f64>, f64)>,
}

impl fmt::Debug for StatisticalManifold {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("StatisticalManifold")
            .field("dim", &self.dim)
            .field("log_density", &"<fn(x, theta) -> f64>")
            .field("nodes", &self.nodes.len())
            .finish()
    }
}

/// Composite trapezoid rule with `n` nodes on `[lo, hi]` (requires `n ≥ 2`).
fn trapezoid_nodes(lo: f64, hi: f64, n: usize) -> Vec<(Vec<f64>, f64)> {
    let spacing = (hi - lo) / (n - 1) as f64;
    (0..n)
        .map(|i| {
            let x = lo + spacing * i as f64;
            let w = if i == 0 || i + 1 == n {
                0.5 * spacing
            } else {
                spacing
            };
            (vec![x], w)
        })
        .collect()
}

/// Normalised log-density of `N(μ, σ²)` with `θ = (μ, σ)`.
fn gaussian_log_density(x: &[f64], theta: &[f64]) -> f64 {
    let (Some(&xv), Some(&mu), Some(&sigma)) = (x.first(), theta.first(), theta.get(1)) else {
        return f64::NAN;
    };
    if sigma.is_nan() || sigma <= 0.0 {
        return f64::NAN;
    }
    let z = (xv - mu) / sigma;
    -0.5 * z * z - sigma.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln()
}

impl StatisticalManifold {
    /// Create a statistical manifold of parameter dimension `dim`.
    ///
    /// * `log_density(x, θ)` must be the normalised log-density with respect to
    ///   the quadrature measure.
    /// * `nodes` are `(x_n, w_n)` pairs; weights must be finite and
    ///   non-negative with at least one positive weight, and every node
    ///   coordinate must be finite.
    pub fn new<F>(dim: usize, log_density: F, nodes: Vec<(Vec<f64>, f64)>) -> PhysicsResult<Self>
    where
        F: Fn(&[f64], &[f64]) -> f64 + Send + Sync + 'static,
    {
        if dim == 0 {
            return Err(invalid("dim", "parameter dimension must be positive"));
        }
        if nodes.is_empty() {
            return Err(invalid("nodes", "at least one quadrature node is required"));
        }
        let mut has_positive = false;
        for (x, w) in &nodes {
            if !w.is_finite() || *w < 0.0 {
                return Err(invalid(
                    "nodes",
                    "quadrature weights must be finite and non-negative",
                ));
            }
            if x.iter().any(|v| !v.is_finite()) {
                return Err(invalid("nodes", "node coordinates must be finite"));
            }
            has_positive |= *w > 0.0;
        }
        if !has_positive {
            return Err(invalid("nodes", "at least one weight must be positive"));
        }
        Ok(Self {
            dim,
            log_density: Arc::new(log_density),
            nodes,
        })
    }

    /// The univariate Gaussian family `N(μ, σ²)` with parameters `θ = (μ, σ)`.
    ///
    /// Integrals are evaluated with a 481-node trapezoid rule on `[−12, 12]`,
    /// which is spectrally accurate for smooth, rapidly decaying integrands.
    /// The rule is valid while the probability mass lies inside the grid and is
    /// resolved by it, roughly `|μ| + 8σ ≤ 12` and `σ ≥ 0.1`; use
    /// [`StatisticalManifold::gaussian_1d_on_grid`] otherwise.  For geodesic
    /// computations the grid must also cover the intermediate shooting
    /// iterates, which can overshoot the endpoints.
    pub fn gaussian_1d() -> Self {
        Self {
            dim: 2,
            log_density: Arc::new(gaussian_log_density),
            nodes: trapezoid_nodes(-12.0, 12.0, 481),
        }
    }

    /// The univariate Gaussian family `N(μ, σ²)`, `θ = (μ, σ)`, integrated with
    /// an `n`-node trapezoid rule on `[lo, hi]`.
    pub fn gaussian_1d_on_grid(lo: f64, hi: f64, n: usize) -> PhysicsResult<Self> {
        if !lo.is_finite() || !hi.is_finite() || lo >= hi {
            return Err(invalid("grid", "require finite lo < hi"));
        }
        if n < 3 {
            return Err(invalid("n", "the trapezoid grid needs at least 3 nodes"));
        }
        Ok(Self {
            dim: 2,
            log_density: Arc::new(gaussian_log_density),
            nodes: trapezoid_nodes(lo, hi, n),
        })
    }

    /// Dimension of the parameter space.
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Quadrature nodes `(x_n, w_n)` over the sample space.
    pub fn nodes(&self) -> &[(Vec<f64>, f64)] {
        &self.nodes
    }

    /// Evaluate the log-density `ℓ(x; θ)`.
    pub fn log_density(&self, x: &[f64], params: &[f64]) -> f64 {
        (self.log_density)(x, params)
    }

    fn check_params(&self, params: &[f64]) -> PhysicsResult<()> {
        if params.len() != self.dim {
            return Err(invalid(
                "params",
                &format!("expected {} parameters, got {}", self.dim, params.len()),
            ));
        }
        if params.iter().any(|v| !v.is_finite()) {
            return Err(invalid("params", "parameters must be finite"));
        }
        Ok(())
    }

    /// Evaluate `ℓ` at a perturbed parameter, requiring a finite value.
    fn eval_perturbed(&self, x: &[f64], theta: &[f64]) -> PhysicsResult<f64> {
        let v = (self.log_density)(x, theta);
        if v.is_finite() {
            Ok(v)
        } else {
            Err(invalid(
                "params",
                "log-density is not finite in a neighbourhood of params (parameter on the \
                 boundary of the domain or support depends on params)",
            ))
        }
    }

    /// Visit every quadrature node of positive density with its mass
    /// `w_n p(x_n; θ)` and the score (and, if `with_hessian`, the Hessian) of
    /// `ℓ` at `params`, after checking that the total mass is one.
    ///
    /// Scores use a central-difference step of `≈ ε^{1/3}` and Hessians a step
    /// of `≈ ε^{1/4}` (relative to `max(1, |θ_i|)`).
    fn visit_nodes<V>(&self, params: &[f64], with_hessian: bool, mut visit: V) -> PhysicsResult<()>
    where
        V: FnMut(f64, &NodeJet),
    {
        self.check_params(params)?;
        let n = self.dim;
        let steps: Vec<f64> = params
            .iter()
            .map(|&v| fd_step(SCORE_STEP_SCALE, v))
            .collect();
        let hess_steps: Vec<f64> = params
            .iter()
            .map(|&v| fd_step(HESSIAN_STEP_SCALE, v))
            .collect();
        let mut jet = NodeJet {
            score: vec![0.0; n],
            hessian: if with_hessian {
                vec![0.0; n * n]
            } else {
                Vec::new()
            },
        };
        let mut total_mass = 0.0f64;
        let mut work = params.to_vec();
        for (x, w) in &self.nodes {
            if *w == 0.0 {
                continue;
            }
            let l0 = (self.log_density)(x, params);
            if l0 == f64::NEG_INFINITY {
                continue;
            }
            if !l0.is_finite() {
                return Err(invalid(
                    "params",
                    "log-density is NaN or +inf at params (outside the parameter domain)",
                ));
            }
            let mass = w * l0.exp();
            total_mass += mass;
            if mass == 0.0 {
                continue;
            }
            for i in 0..n {
                let h = steps[i];
                work[i] = params[i] + h;
                let lp = self.eval_perturbed(x, &work)?;
                work[i] = params[i] - h;
                let lm = self.eval_perturbed(x, &work)?;
                work[i] = params[i];
                jet.score[i] = (lp - lm) / (2.0 * h);
            }
            if with_hessian {
                for i in 0..n {
                    let hi = hess_steps[i];
                    work[i] = params[i] + hi;
                    let lp = self.eval_perturbed(x, &work)?;
                    work[i] = params[i] - hi;
                    let lm = self.eval_perturbed(x, &work)?;
                    work[i] = params[i];
                    jet.hessian[i * n + i] = (lp - 2.0 * l0 + lm) / (hi * hi);
                    for j in (i + 1)..n {
                        let hj = hess_steps[j];
                        let mut corner = |si: f64, sj: f64| -> PhysicsResult<f64> {
                            work[i] = params[i] + si * hi;
                            work[j] = params[j] + sj * hj;
                            let v = self.eval_perturbed(x, &work);
                            work[i] = params[i];
                            work[j] = params[j];
                            v
                        };
                        let lpp = corner(1.0, 1.0)?;
                        let lpm = corner(1.0, -1.0)?;
                        let lmp = corner(-1.0, 1.0)?;
                        let lmm = corner(-1.0, -1.0)?;
                        let v = (lpp - lpm - lmp + lmm) / (4.0 * hi * hj);
                        jet.hessian[i * n + j] = v;
                        jet.hessian[j * n + i] = v;
                    }
                }
            }
            visit(mass, &jet);
        }
        if (total_mass - 1.0).abs() > MASS_TOLERANCE {
            return Err(invalid(
                "nodes",
                &format!(
                    "quadrature mass of the density is {total_mass}, not 1: the log-density is \
                     not normalised or the nodes do not cover its support"
                ),
            ));
        }
        Ok(())
    }

    /// Fisher information metric `g_ij(θ) = E_θ[∂_iℓ ∂_jℓ]` at `params`.
    ///
    /// Returns a symmetric positive semi-definite `dim × dim` matrix.
    pub fn fisher_metric(&self, params: &[f64]) -> PhysicsResult<Vec<Vec<f64>>> {
        let n = self.dim;
        let mut g = vec![0.0f64; n * n];
        self.visit_nodes(params, false, |mass, jet| {
            accumulate_metric(&mut g, mass, &jet.score);
        })?;
        Ok(unflatten_symmetric(&g, n))
    }

    /// Amari–Chentsov (skewness) tensor `T_ijk(θ) = E_θ[∂_iℓ ∂_jℓ ∂_kℓ]`,
    /// indexed `[i][j][k]` (fully symmetric).
    pub fn skewness_tensor(&self, params: &[f64]) -> PhysicsResult<Vec<Vec<Vec<f64>>>> {
        let n = self.dim;
        let mut t = vec![vec![vec![0.0f64; n]; n]; n];
        self.visit_nodes(params, false, |mass, jet| {
            let s = &jet.score;
            for (i, t_i) in t.iter_mut().enumerate() {
                for (j, t_ij) in t_i.iter_mut().enumerate() {
                    let sij = mass * s[i] * s[j];
                    for (t_ijk, &sk) in t_ij.iter_mut().zip(s) {
                        *t_ijk += sij * sk;
                    }
                }
            }
        })?;
        Ok(t)
    }

    /// α-connection coefficients of the first kind
    /// `Γ^(α)_{ij,k} = E_θ[(∂_i∂_jℓ + (1−α)/2 ∂_iℓ ∂_jℓ) ∂_kℓ]`, indexed
    /// `[i][j][k]`.
    ///
    /// `α = 0` gives the Levi-Civita connection of the Fisher metric, `α = 1`
    /// the exponential (e-) connection and `α = −1` the mixture (m-)
    /// connection.  They satisfy `Γ^(α) = Γ^(0) − (α/2) T`.
    pub fn alpha_christoffel_first_kind(
        &self,
        params: &[f64],
        alpha: f64,
    ) -> PhysicsResult<Vec<Vec<Vec<f64>>>> {
        let n = self.dim;
        let mut first = vec![0.0f64; n * n * n];
        self.visit_nodes(params, true, |mass, jet| {
            accumulate_first_kind(&mut first, mass, jet, n, alpha);
        })?;
        Ok((0..n)
            .map(|i| {
                (0..n)
                    .map(|j| (0..n).map(|k| first[(i * n + j) * n + k]).collect())
                    .collect()
            })
            .collect())
    }

    /// α-connection coefficients of the second kind `Γ^(α)k_{ij}`, indexed
    /// `[k][i][j]`, obtained by raising the last index of
    /// [`StatisticalManifold::alpha_christoffel_first_kind`] with `g⁻¹`.
    pub fn alpha_christoffel_symbols(
        &self,
        params: &[f64],
        alpha: f64,
    ) -> PhysicsResult<Vec<Vec<Vec<f64>>>> {
        let n = self.dim;
        let mut g = vec![0.0f64; n * n];
        let mut first = vec![0.0f64; n * n * n];
        self.visit_nodes(params, true, |mass, jet| {
            accumulate_metric(&mut g, mass, &jet.score);
            accumulate_first_kind(&mut first, mass, jet, n, alpha);
        })?;
        let g_inv = invert_matrix(&unflatten_symmetric(&g, n)).ok_or_else(singular_metric)?;
        let mut gamma = vec![vec![vec![0.0f64; n]; n]; n];
        for (l, gamma_l) in gamma.iter_mut().enumerate() {
            for (i, gamma_li) in gamma_l.iter_mut().enumerate() {
                for (j, gamma_lij) in gamma_li.iter_mut().enumerate() {
                    *gamma_lij = (0..n)
                        .map(|k| g_inv[l][k] * first[(i * n + j) * n + k])
                        .sum();
                }
            }
        }
        Ok(gamma)
    }

    /// Christoffel symbols `Γ^k_{ij}` of the Levi-Civita connection of the
    /// Fisher metric at `params`, indexed `[k][i][j]`.
    pub fn christoffel_symbols(&self, params: &[f64]) -> PhysicsResult<Vec<Vec<Vec<f64>>>> {
        self.alpha_christoffel_symbols(params, 0.0)
    }

    /// Geodesic acceleration `a^k = −Γ^k_{ij} u^i u^j` at `x` with velocity `u`.
    fn geodesic_acceleration(&self, x: &[f64], u: &[f64]) -> PhysicsResult<Vec<f64>> {
        let gamma = self.christoffel_symbols(x)?;
        Ok(gamma
            .iter()
            .map(|gk| {
                -gk.iter()
                    .zip(u)
                    .map(|(row, &ui)| ui * row.iter().zip(u).map(|(g, &uj)| g * uj).sum::<f64>())
                    .sum::<f64>()
            })
            .collect())
    }

    /// Exponential map of the Fisher–Rao (Levi-Civita) connection: the point
    /// `γ(t)` of the geodesic with `γ(0) = p` and `γ'(0) = v`.
    ///
    /// The geodesic equation `γ'' + Γ(γ)[γ', γ'] = 0` is integrated with a
    /// fixed-step classical Runge–Kutta scheme (64 steps over `[0, t]`).
    pub fn exponential_map(&self, p: &[f64], v: &[f64], t: f64) -> PhysicsResult<Vec<f64>> {
        self.check_params(p)?;
        if v.len() != self.dim || v.iter().any(|c| !c.is_finite()) {
            return Err(invalid(
                "v",
                "tangent vector must be finite with length dim",
            ));
        }
        if !t.is_finite() {
            return Err(invalid("t", "must be finite"));
        }
        let h = t / GEODESIC_RK4_STEPS as f64;
        let mut x = p.to_vec();
        let mut u = v.to_vec();
        let axpy = |base: &[f64], k: &[f64], s: f64| -> Vec<f64> {
            base.iter().zip(k).map(|(b, kk)| b + s * kk).collect()
        };
        for _ in 0..GEODESIC_RK4_STEPS {
            let a1 = self.geodesic_acceleration(&x, &u)?;
            let x2 = axpy(&x, &u, 0.5 * h);
            let u2 = axpy(&u, &a1, 0.5 * h);
            let a2 = self.geodesic_acceleration(&x2, &u2)?;
            let x3 = axpy(&x, &u2, 0.5 * h);
            let u3 = axpy(&u, &a2, 0.5 * h);
            let a3 = self.geodesic_acceleration(&x3, &u3)?;
            let x4 = axpy(&x, &u3, h);
            let u4 = axpy(&u, &a3, h);
            let a4 = self.geodesic_acceleration(&x4, &u4)?;
            for idx in 0..self.dim {
                x[idx] += h / 6.0 * (u[idx] + 2.0 * u2[idx] + 2.0 * u3[idx] + u4[idx]);
                u[idx] += h / 6.0 * (a1[idx] + 2.0 * a2[idx] + 2.0 * a3[idx] + a4[idx]);
            }
        }
        Ok(x)
    }

    /// Logarithmic map: the initial velocity `v` of the geodesic with
    /// `γ(0) = p` and `γ(1) = q`, found by damped quasi-Newton shooting on
    /// `exp_p(v) − q = 0` starting from `v = q − p` (finite-difference
    /// Jacobian with Broyden updates).
    ///
    /// Returns [`PhysicsError::ConvergenceFailed`] if the shooting iteration
    /// does not converge (for example when `q` lies beyond a conjugate point),
    /// and propagates [`PhysicsError::InvalidInput`] if a trial geodesic leaves
    /// the parameter domain or the region covered by the quadrature nodes.
    pub fn log_map(&self, p: &[f64], q: &[f64]) -> PhysicsResult<Vec<f64>> {
        self.check_params(p)?;
        self.check_params(q)?;
        let n = self.dim;
        let q_scale = q.iter().fold(1.0f64, |m, v| m.max(v.abs()));
        let tol = SHOOTING_TOL * q_scale;
        let residual = |v: &[f64]| -> PhysicsResult<Vec<f64>> {
            let end = self.exponential_map(p, v, 1.0)?;
            Ok(end.iter().zip(q).map(|(e, qq)| e - qq).collect())
        };
        let inf_norm = |r: &[f64]| r.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let mut v: Vec<f64> = q.iter().zip(p).map(|(qq, pp)| qq - pp).collect();
        let mut r = residual(&v)?;
        let mut r_norm = inf_norm(&r);
        // Jacobian of the shooting map: a forward-difference Jacobian that is
        // kept up to date with Broyden rank-one updates and recomputed whenever
        // the quasi-Newton step fails.
        let mut jac: Option<Vec<Vec<f64>>> = None;
        for _ in 0..SHOOTING_MAX_ITER {
            if r_norm <= tol {
                return Ok(v);
            }
            let fresh = jac.is_none();
            let current = match jac.take() {
                Some(j) => j,
                None => {
                    let v_scale = v.iter().fold(1.0f64, |m, c| m.max(c.abs()));
                    let delta = 1.0e-6 * v_scale;
                    let mut fd = vec![vec![0.0f64; n]; n];
                    for col in 0..n {
                        let mut vp = v.clone();
                        vp[col] += delta;
                        let rp = residual(&vp)?;
                        for (row, jrow) in fd.iter_mut().enumerate() {
                            jrow[col] = (rp[row] - r[row]) / delta;
                        }
                    }
                    fd
                }
            };
            let Some(jac_inv) = invert_matrix(&current) else {
                if fresh {
                    break;
                }
                continue;
            };
            let step: Vec<f64> = jac_inv
                .iter()
                .map(|row| -row.iter().zip(&r).map(|(a, b)| a * b).sum::<f64>())
                .collect();
            // Damped update: halve the step until the residual decreases.
            let mut lambda = 1.0f64;
            let mut accepted = None;
            for _ in 0..30 {
                let trial: Vec<f64> = v.iter().zip(&step).map(|(a, s)| a + lambda * s).collect();
                if let Ok(rt) = residual(&trial) {
                    let rt_norm = inf_norm(&rt);
                    if rt_norm < r_norm {
                        accepted = Some((trial, rt, rt_norm));
                        break;
                    }
                }
                lambda *= 0.5;
            }
            let Some((v_new, r_new, r_new_norm)) = accepted else {
                if fresh {
                    break;
                }
                // Stale Broyden Jacobian: retry with a finite-difference one.
                continue;
            };
            // Broyden update J += (y − J s) sᵀ / (sᵀ s).
            let s_vec: Vec<f64> = v_new.iter().zip(&v).map(|(a, b)| a - b).collect();
            let s_norm2: f64 = s_vec.iter().map(|c| c * c).sum();
            let mut updated = current;
            if s_norm2 > 0.0 {
                for (row, jrow) in updated.iter_mut().enumerate() {
                    let js: f64 = jrow.iter().zip(&s_vec).map(|(a, b)| a * b).sum();
                    let coeff = (r_new[row] - r[row] - js) / s_norm2;
                    for (cell, &sc) in jrow.iter_mut().zip(&s_vec) {
                        *cell += coeff * sc;
                    }
                }
            }
            jac = Some(updated);
            v = v_new;
            r = r_new;
            r_norm = r_new_norm;
        }
        if r_norm <= tol {
            return Ok(v);
        }
        Err(PhysicsError::ConvergenceFailed {
            iterations: SHOOTING_MAX_ITER,
            tolerance: tol,
            residual: r_norm,
        })
    }

    /// Point at time `t` of the Fisher–Rao geodesic from `p` (`t = 0`) to `q`
    /// (`t = 1`).
    ///
    /// The boundary-value problem is solved by shooting
    /// ([`StatisticalManifold::log_map`]) and the geodesic is then evaluated
    /// with [`StatisticalManifold::exponential_map`]; `t` outside `[0, 1]`
    /// extrapolates along the same geodesic.  Every call re-solves the
    /// boundary-value problem; to sample many points of one geodesic call
    /// `log_map` once and then `exponential_map` for each `t`.
    pub fn geodesic(&self, p: &[f64], q: &[f64], t: f64) -> PhysicsResult<Vec<f64>> {
        let v = self.log_map(p, q)?;
        self.exponential_map(p, &v, t)
    }

    /// Fisher–Rao geodesic distance `√(vᵀ g(p) v)` with `v = log_p(q)`.
    pub fn geodesic_distance(&self, p: &[f64], q: &[f64]) -> PhysicsResult<f64> {
        let v = self.log_map(p, q)?;
        let g = self.fisher_metric(p)?;
        let quad: f64 = g
            .iter()
            .zip(&v)
            .map(|(row, &vi)| vi * row.iter().zip(&v).map(|(a, &vj)| a * vj).sum::<f64>())
            .sum();
        Ok(quad.max(0.0).sqrt())
    }
}

/// Add `mass · s sᵀ` to the upper triangle of the row-major matrix `g`.
fn accumulate_metric(g: &mut [f64], mass: f64, score: &[f64]) {
    let n = score.len();
    for i in 0..n {
        let mi = mass * score[i];
        for j in i..n {
            g[i * n + j] += mi * score[j];
        }
    }
}

/// Expand the upper triangle of a row-major `n × n` matrix into a symmetric
/// nested matrix.
fn unflatten_symmetric(g: &[f64], n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i <= j { g[i * n + j] } else { g[j * n + i] })
                .collect()
        })
        .collect()
}

/// Add `mass · (H_ij + (1−α)/2 s_i s_j) s_k` to the flat `[i][j][k]` tensor.
fn accumulate_first_kind(first: &mut [f64], mass: f64, jet: &NodeJet, n: usize, alpha: f64) {
    let c = 0.5 * (1.0 - alpha);
    let s = &jet.score;
    for i in 0..n {
        for j in 0..n {
            let a = mass * (jet.hessian[i * n + j] + c * s[i] * s[j]);
            let base = (i * n + j) * n;
            for (slot, &sk) in first[base..base + n].iter_mut().zip(s) {
                *slot += a * sk;
            }
        }
    }
}

fn invalid(field: &str, reason: &str) -> PhysicsError {
    PhysicsError::InvalidInput {
        field: field.to_string(),
        reason: reason.to_string(),
    }
}

fn singular_metric() -> PhysicsError {
    invalid(
        "params",
        "Fisher metric is singular at params (parameterisation not identifiable)",
    )
}

#[cfg(test)]
#[path = "statistical_manifold_tests.rs"]
mod tests;
