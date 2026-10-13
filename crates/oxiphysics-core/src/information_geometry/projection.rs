// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Information (KL) projections of discrete distributions onto exponential
//! families.

use nalgebra::{DMatrix, DVector};

use crate::error::{PhysicsError, PhysicsResult};

/// Maximum number of damped Newton iterations.
const NEWTON_MAX_ITER: usize = 200;
/// Maximum number of step halvings in the backtracking line search.
const LINE_SEARCH_MAX_HALVINGS: usize = 60;
/// Armijo sufficient-decrease constant.
const ARMIJO_C: f64 = 1.0e-4;
/// Base gradient tolerance (scaled by the magnitude of the statistics).
const GRADIENT_TOL: f64 = 1.0e-12;
/// Probability mass that the forward projection may place on points outside
/// the support of `p` before `E_p[T]` is deemed to lie on the boundary of the
/// marginal polytope.
const BOUNDARY_MASS_TOL: f64 = 1.0e-10;

/// Information projection onto a discrete exponential family.
///
/// The sample space is the finite support `x_j = j` (`j = 0, …, n−1`, as the
/// one-element slice `[j as f64]`), where `n` is the length of the probability
/// vector being projected.  The target family is
///
/// `q_θ(x_j) = exp(θ · T(x_j) − A(θ))`, `A(θ) = ln Σ_j exp(θ · T(x_j))`,
///
/// with sufficient statistics `T = (T_1, …, T_k)` given by `target_family`.
/// The statistics must be affinely independent on the support (in particular
/// `k ≤ n − 1`), so that `θ` is identifiable.
pub struct InformationProjection {
    /// Sufficient statistics defining the target exponential family.
    pub target_family: Vec<fn(&[f64]) -> f64>,
}

/// Moments of `q_θ` on the support.
struct FamilyMoments {
    /// Log-partition function `A(θ)`.
    log_partition: f64,
    /// `ln q_θ(x_j)`.
    log_q: Vec<f64>,
    /// `q_θ(x_j)`.
    q: Vec<f64>,
    /// `E_θ[T]`.
    mean: DVector<f64>,
    /// `Cov_θ[T]`.
    cov: DMatrix<f64>,
}

/// Objective value with gradient, Newton Hessian and a positive-definite
/// fallback curvature matrix.
struct Objective {
    value: f64,
    grad: DVector<f64>,
    hessian: DMatrix<f64>,
    fallback: DMatrix<f64>,
}

impl InformationProjection {
    /// Create a new `InformationProjection` onto the family with the given
    /// sufficient statistics.
    pub fn new(target_family: Vec<fn(&[f64]) -> f64>) -> Self {
        Self { target_family }
    }

    /// Forward (moment-matching) projection
    /// `argmin_θ KL(p ‖ q_θ)`, characterised by `E_θ[T] = E_p[T]`.
    ///
    /// `p` is a probability vector over `x_j = j`; it is normalised by its sum.
    /// The dual problem `min_θ A(θ) − θ · E_p[T]` is strictly convex with
    /// Hessian `Cov_θ[T]` and is solved by damped Newton iteration.
    ///
    /// Returns the natural parameters `θ`.  Errors:
    /// * [`PhysicsError::InvalidInput`] for an invalid `p`, an empty family,
    ///   affinely dependent statistics, or when `E_p[T]` lies on the boundary
    ///   of the convex hull of `{T(x_j)}` (the infimum is then only approached
    ///   as `‖θ‖ → ∞`);
    /// * [`PhysicsError::ConvergenceFailed`] if the Newton iteration stalls.
    pub fn project(&self, p: &[f64]) -> PhysicsResult<Vec<f64>> {
        let p = normalise(p)?;
        let table = self.stats_table(p.len())?;
        let k = self.target_family.len();
        check_identifiable(&table, k)?;
        let eta = DVector::from_iterator(
            k,
            (0..k).map(|a| table.iter().zip(&p).map(|(t, pj)| pj * t[a]).sum::<f64>()),
        );
        let tol = GRADIENT_TOL * stats_scale(&table);
        let theta = damped_newton(
            vec![0.0; k],
            |theta| {
                let m = family_moments(&table, theta)?;
                let value = m.log_partition - dot(theta, eta.as_slice());
                Ok(Objective {
                    value,
                    grad: &m.mean - &eta,
                    hessian: m.cov.clone(),
                    fallback: m.cov,
                })
            },
            tol,
        )?;
        // A finite minimiser puts positive mass everywhere; if (numerically)
        // all mass sits on the support of p while p misses some points, E_p[T]
        // is on a face of the marginal polytope and θ has diverged.
        if p.contains(&0.0) {
            let m = family_moments(&table, &theta)?;
            let off_support: f64 =
                m.q.iter()
                    .zip(&p)
                    .filter(|(_, pj)| **pj == 0.0)
                    .map(|(qj, _)| qj)
                    .sum();
            if off_support <= BOUNDARY_MASS_TOL {
                return Err(invalid(
                    "p",
                    "E_p[T] lies on the boundary of the convex hull of the sufficient \
                     statistics; the projection is not attained at finite natural parameters",
                ));
            }
        }
        Ok(theta)
    }

    /// Reverse (information / I-) projection `argmin_θ KL(q_θ ‖ p)`.
    ///
    /// The gradient is `Cov_θ[T, ln q_θ − ln p]` and the Hessian
    /// `Cov_θ[T] + E_θ[(T − μ)(T − μ)ᵀ (ln q_θ − ln p − KL)]`; Newton steps
    /// are used where the Hessian is positive definite and natural-gradient
    /// steps (`Cov_θ[T]⁻¹ ∇`) otherwise, with Armijo backtracking.  The reverse
    /// KL is not convex in `θ`, so the result is the local minimiser reached
    /// from `init_theta`.
    ///
    /// Errors: [`PhysicsError::InvalidInput`] if `p` has a zero entry
    /// (`KL(q‖p)` is then infinite for every member of the family), is
    /// otherwise invalid, if `init_theta` has the wrong length or the
    /// statistics are affinely dependent; [`PhysicsError::ConvergenceFailed`]
    /// if the iteration stalls.
    pub fn reverse_kl_projection(&self, p: &[f64], init_theta: &[f64]) -> PhysicsResult<Vec<f64>> {
        let p = normalise(p)?;
        if p.contains(&0.0) {
            return Err(invalid(
                "p",
                "reverse KL projection requires p_j > 0 for every support point",
            ));
        }
        let k = self.target_family.len();
        if init_theta.len() != k {
            return Err(invalid(
                "init_theta",
                &format!("expected {k} natural parameters, got {}", init_theta.len()),
            ));
        }
        if init_theta.iter().any(|v| !v.is_finite()) {
            return Err(invalid("init_theta", "must be finite"));
        }
        let table = self.stats_table(p.len())?;
        check_identifiable(&table, k)?;
        let log_p: Vec<f64> = p.iter().map(|pj| pj.ln()).collect();
        let log_p_scale = log_p.iter().fold(1.0f64, |m, v| m.max(v.abs()));
        let tol = GRADIENT_TOL * stats_scale(&table) * log_p_scale;
        damped_newton(
            init_theta.to_vec(),
            |theta| {
                let m = family_moments(&table, theta)?;
                let r: Vec<f64> = m.log_q.iter().zip(&log_p).map(|(lq, lp)| lq - lp).collect();
                let kl: f64 = m.q.iter().zip(&r).map(|(qj, rj)| qj * rj).sum();
                let mut grad = DVector::zeros(k);
                let mut hessian = m.cov.clone();
                for ((t, qj), rj) in table.iter().zip(&m.q).zip(&r) {
                    let centred = DVector::from_iterator(k, (0..k).map(|a| t[a] - m.mean[a]));
                    grad += &centred * (qj * rj);
                    hessian += &centred * centred.transpose() * (qj * (rj - kl));
                }
                Ok(Objective {
                    value: kl,
                    grad,
                    hessian,
                    fallback: m.cov,
                })
            },
            tol,
        )
    }

    /// The distribution `q_θ(x_j)`, `j = 0, …, n−1`, of the family member with
    /// natural parameters `theta`.
    pub fn distribution(&self, theta: &[f64], n: usize) -> PhysicsResult<Vec<f64>> {
        if theta.len() != self.target_family.len() {
            return Err(invalid(
                "theta",
                &format!(
                    "expected {} natural parameters, got {}",
                    self.target_family.len(),
                    theta.len()
                ),
            ));
        }
        if n == 0 {
            return Err(invalid("n", "support must be non-empty"));
        }
        let table = self.stats_table(n)?;
        Ok(family_moments(&table, theta)?.q)
    }

    /// Table `T_a(x_j)` with rows indexed by support point `j`.
    fn stats_table(&self, n: usize) -> PhysicsResult<Vec<Vec<f64>>> {
        if self.target_family.is_empty() {
            return Err(invalid(
                "target_family",
                "at least one sufficient statistic is required",
            ));
        }
        let table: Vec<Vec<f64>> = (0..n)
            .map(|j| {
                let x = [j as f64];
                self.target_family.iter().map(|t| t(&x)).collect()
            })
            .collect();
        if table.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid(
                "target_family",
                "sufficient statistics must be finite on the support",
            ));
        }
        Ok(table)
    }
}

/// Validate and normalise a probability vector.
fn normalise(p: &[f64]) -> PhysicsResult<Vec<f64>> {
    if p.is_empty() {
        return Err(invalid("p", "probability vector must be non-empty"));
    }
    if p.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(invalid("p", "entries must be finite and non-negative"));
    }
    let total: f64 = p.iter().sum();
    if total <= 0.0 || !total.is_finite() {
        return Err(invalid("p", "entries must have a positive finite sum"));
    }
    Ok(p.iter().map(|v| v / total).collect())
}

/// Largest absolute value of the sufficient statistics (at least 1).
fn stats_scale(table: &[Vec<f64>]) -> f64 {
    table.iter().flatten().fold(1.0f64, |m, v| m.max(v.abs()))
}

/// Require `Cov_unif[T]` to be positive definite, i.e. the statistics are
/// affinely independent on the support.
fn check_identifiable(table: &[Vec<f64>], k: usize) -> PhysicsResult<()> {
    let m = family_moments(table, &vec![0.0; k])?;
    if m.cov.cholesky().is_none() {
        return Err(invalid(
            "target_family",
            "sufficient statistics are affinely dependent on the support (need k ≤ n − 1 \
             non-degenerate statistics); natural parameters are not identifiable",
        ));
    }
    Ok(())
}

/// Compute `ln q_θ`, `q_θ`, `E_θ[T]` and `Cov_θ[T]` with a stable
/// log-sum-exp.
fn family_moments(table: &[Vec<f64>], theta: &[f64]) -> PhysicsResult<FamilyMoments> {
    let k = theta.len();
    let s: Vec<f64> = table.iter().map(|t| dot(theta, t)).collect();
    let s_max = s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !s_max.is_finite() {
        return Err(PhysicsError::NumericalDivergence {
            message: "natural parameters produced a non-finite exponent".to_string(),
        });
    }
    let sum_exp: f64 = s.iter().map(|v| (v - s_max).exp()).sum();
    let log_partition = s_max + sum_exp.ln();
    let log_q: Vec<f64> = s.iter().map(|v| v - log_partition).collect();
    let q: Vec<f64> = log_q.iter().map(|v| v.exp()).collect();
    let mut mean = DVector::zeros(k);
    for (t, qj) in table.iter().zip(&q) {
        for a in 0..k {
            mean[a] += qj * t[a];
        }
    }
    let mut cov = DMatrix::zeros(k, k);
    for (t, qj) in table.iter().zip(&q) {
        for a in 0..k {
            let da = t[a] - mean[a];
            for b in 0..k {
                cov[(a, b)] += qj * da * (t[b] - mean[b]);
            }
        }
    }
    Ok(FamilyMoments {
        log_partition,
        log_q,
        q,
        mean,
        cov,
    })
}

/// Damped Newton minimisation with Armijo backtracking.
///
/// Uses the Newton direction when `hessian` is positive definite and the
/// `fallback` (positive-definite) metric otherwise.  Converges when
/// `‖∇f‖∞ ≤ tol`.
fn damped_newton<F>(theta0: Vec<f64>, evaluate: F, tol: f64) -> PhysicsResult<Vec<f64>>
where
    F: Fn(&[f64]) -> PhysicsResult<Objective>,
{
    let mut theta = theta0;
    let mut obj = evaluate(&theta)?;
    for iteration in 0..NEWTON_MAX_ITER {
        let grad_norm = obj.grad.amax();
        if grad_norm <= tol {
            return Ok(theta);
        }
        let direction = match obj.hessian.clone().cholesky() {
            Some(chol) => chol.solve(&(-&obj.grad)),
            None => match obj.fallback.clone().cholesky() {
                Some(chol) => chol.solve(&(-&obj.grad)),
                None => {
                    return Err(PhysicsError::ConvergenceFailed {
                        iterations: iteration,
                        tolerance: tol,
                        residual: grad_norm,
                    });
                }
            },
        };
        let slope = obj.grad.dot(&direction);
        // Rounding slack lets the iteration finish once f is flat to machine
        // precision near the optimum.
        let slack = 8.0 * f64::EPSILON * obj.value.abs().max(1.0);
        let mut lambda = 1.0f64;
        let mut next = None;
        for _ in 0..LINE_SEARCH_MAX_HALVINGS {
            let trial: Vec<f64> = theta
                .iter()
                .zip(direction.iter())
                .map(|(t, d)| t + lambda * d)
                .collect();
            if let Ok(trial_obj) = evaluate(&trial)
                && trial_obj.value <= obj.value + ARMIJO_C * lambda * slope + slack
            {
                next = Some((trial, trial_obj));
                break;
            }
            lambda *= 0.5;
        }
        match next {
            Some((t, o)) => {
                theta = t;
                obj = o;
            }
            None => {
                return Err(PhysicsError::ConvergenceFailed {
                    iterations: iteration,
                    tolerance: tol,
                    residual: grad_norm,
                });
            }
        }
    }
    if obj.grad.amax() <= tol {
        return Ok(theta);
    }
    Err(PhysicsError::ConvergenceFailed {
        iterations: NEWTON_MAX_ITER,
        tolerance: tol,
        residual: obj.grad.amax(),
    })
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn invalid(field: &str, reason: &str) -> PhysicsError {
    PhysicsError::InvalidInput {
        field: field.to_string(),
        reason: reason.to_string(),
    }
}

#[cfg(test)]
#[path = "projection_tests.rs"]
mod tests;
