// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lambert's problem: Izzo's (2015) universal-variable solver.
//!
//! Given two position vectors and a time of flight, [`lambert_solver`] finds
//! the velocities of the two-body conic that connects them.  The method
//! follows D. Izzo, "Revisiting Lambert's problem", *Celestial Mechanics and
//! Dynamical Astronomy* 121 (2015).
//!
//! # Formulation
//!
//! With the chord `c = |r2 - r1|`, the semi-perimeter `s = (r1 + r2 + c)/2`
//! and `λ² = 1 - c/s` (`λ < 0` when the transfer angle exceeds π), the
//! problem reduces to one scalar equation `T(x) = T*` for the free variable
//! `x`, where `T* = tof · √(2μ/s³)` and the transfer ellipse has semi-major
//! axis `a = s / (2 (1 - x²))`.  The branch `x < 1` is elliptic, `x = 1`
//! parabolic and `x > 1` hyperbolic, so every transfer is covered by a single
//! equation with `x ∈ (-1, ∞)`:
//!
//! ```text
//! y = √(1 - λ² (1 - x²))
//! T(x) = ( ψ / √|1 - x²| - x + λ y ) / (1 - x²)
//! ψ = atan2( (y - λx) √(1 - x²), x y + λ (1 - x²) )   (x < 1)
//! ψ = asinh( (y - λx) √(x² - 1) )                      (x > 1)
//! ```
//!
//! `T(x)` is smooth and strictly decreasing, with the closed-form
//! derivatives `T' = (3Tx - 2 + 2λ³x/y)/(1 - x²)`,
//! `T'' = (3T + 5xT' + 2(1-λ²)λ³/y³)/(1 - x²)` and
//! `T''' = (7xT'' + 8T' - 6(1-λ²)λ⁵x/y⁵)/(1 - x²)`.  All of these are
//! `0/0` at the parabolic point `x = 1`, so for `√0.6 < x < √1.4` the solver
//! switches to Battin's hypergeometric form
//! `T = (2/3) η³ ₂F₁(3, 1; 5/2; S) + 2 λ η` with `η = y - λx` and
//! `S = (1 - λ - xη)/2`, and obtains the derivatives by differentiating that
//! series with truncated Taylor arithmetic.  This keeps `T` and its
//! derivatives accurate through the parabolic transfer.
//!
//! The root is found with third-order Householder iterations started from
//! Izzo's initial guess and guarded by a bracketing interval, so the
//! iteration cannot leave the domain `x > -1` and always shrinks the
//! bracket when a Householder step is rejected.
//!
//! Only zero-revolution (direct) transfers are supported.

use super::{vec3_add, vec3_cross, vec3_norm, vec3_scale, vec3_sub};
use std::ops::{Add, Mul, Sub};

#[cfg(test)]
mod tests;

/// Lower edge `√0.6` of the interval around `x = 1` where Battin's series is used.
const BATTIN_LOWER: f64 = 0.774_596_669_241_483_4;

/// Upper edge `√1.4` of the interval around `x = 1` where Battin's series is used.
const BATTIN_UPPER: f64 = 1.183_215_956_619_923_2;

/// Iteration limit of the root finder.
const MAX_ITERATIONS: usize = 100;

/// Residual (relative to `max(1, T*)`) at which the root finder stops.
const CONVERGENCE_TOL: f64 = 1e-13;

/// Residual (relative to `max(1, T*)`) accepted when the iteration stalls.
const ACCEPTANCE_TOL: f64 = 1e-11;

/// Minimum `|r1 × r2| / (|r1||r2|)` for the transfer plane to be defined.
const COLLINEAR_TOL: f64 = 1e-9;

/// Slack on the z component of the orbit normal when deciding "prograde".
const PLANE_TOL: f64 = 1e-12;

/// Result of Lambert's problem.
#[derive(Debug, Clone)]
pub struct LambertSolution {
    /// Departure velocity vector (m/s) at `r1`.
    pub v1: [f64; 3],
    /// Arrival velocity vector (m/s) at `r2`.
    pub v2: [f64; 3],
    /// Transfer semi-major axis (m): positive for an ellipse, negative for a
    /// hyperbola and infinite for the parabolic transfer.
    pub semi_major_axis: f64,
    /// Number of solver iterations used.
    pub iterations: usize,
}

/// Solve Lambert's problem (Izzo 2015, zero revolutions).
///
/// Find the orbit connecting positions `r1` and `r2` in transfer time `tof`.
/// Elliptic, parabolic and hyperbolic transfers are all supported, and so are
/// transfer angles above 180°.
///
/// # Arguments
/// * `r1` – Departure position (m) in ECI.
/// * `r2` – Arrival position (m) in ECI.
/// * `tof` – Time of flight (s).
/// * `mu` – Gravitational parameter (m³/s²).
/// * `prograde` – True for a prograde transfer (orbital angular momentum with
///   a non-negative component along +z), false for retrograde.  The transfer
///   angle is then the angle swept in that sense, so it may exceed 180°.
///
/// Returns `None` if the problem is ill-posed (`tof` or `mu` not positive and
/// finite, a zero or non-finite position, or `r1` and `r2` collinear with the
/// centre, where the transfer plane is undefined) or if the solver fails to
/// converge within 100 iterations.
pub fn lambert_solver(
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    mu: f64,
    prograde: bool,
) -> Option<LambertSolution> {
    let finite = |v: [f64; 3]| v.iter().all(|c| c.is_finite());
    if !(tof.is_finite() && tof > 0.0 && mu.is_finite() && mu > 0.0 && finite(r1) && finite(r2)) {
        return None;
    }
    let r1_mag = vec3_norm(r1);
    let r2_mag = vec3_norm(r2);
    if !(r1_mag > 0.0 && r2_mag > 0.0 && r1_mag.is_finite() && r2_mag.is_finite()) {
        return None;
    }

    let normal = vec3_cross(r1, r2);
    let normal_mag = vec3_norm(normal);
    if normal_mag <= COLLINEAR_TOL * r1_mag * r2_mag {
        return None;
    }
    let chord = vec3_norm(vec3_sub(r2, r1));
    let semi_perimeter = 0.5 * (r1_mag + r2_mag + chord);
    let mut lambda = (1.0 - chord / semi_perimeter).max(0.0).sqrt();

    // Orbit normal in the sense of motion, tangential unit vectors from it.
    let mut h_hat = vec3_scale(normal, 1.0 / normal_mag);
    if h_hat[2] < -PLANE_TOL {
        // The prograde transfer is the long way round (angle above 180°).
        lambda = -lambda;
        h_hat = vec3_scale(h_hat, -1.0);
    }
    let r1_hat = vec3_scale(r1, 1.0 / r1_mag);
    let r2_hat = vec3_scale(r2, 1.0 / r2_mag);
    let mut t1_hat = vec3_cross(h_hat, r1_hat);
    let mut t2_hat = vec3_cross(h_hat, r2_hat);
    if !prograde {
        lambda = -lambda;
        t1_hat = vec3_scale(t1_hat, -1.0);
        t2_hat = vec3_scale(t2_hat, -1.0);
    }

    let t_target = (2.0 * mu / semi_perimeter.powi(3)).sqrt() * tof;
    let (x, iterations) = solve_for_x(lambda, t_target)?;

    // Velocity components in the radial/tangential basis (Izzo, section 3.5).
    let y = (1.0 - lambda * lambda * (1.0 - x * x)).sqrt();
    let gamma = (0.5 * mu * semi_perimeter).sqrt();
    let rho = (r1_mag - r2_mag) / chord;
    let sigma = (1.0 - rho * rho).max(0.0).sqrt();
    let v_r1 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r1_mag;
    let v_r2 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r2_mag;
    let v_t1 = gamma * sigma * (y + lambda * x) / r1_mag;
    let v_t2 = gamma * sigma * (y + lambda * x) / r2_mag;

    let v1 = vec3_add(vec3_scale(r1_hat, v_r1), vec3_scale(t1_hat, v_t1));
    let v2 = vec3_add(vec3_scale(r2_hat, v_r2), vec3_scale(t2_hat, v_t2));
    if !(finite(v1) && finite(v2)) {
        return None;
    }

    let one_minus_x2 = 1.0 - x * x;
    let semi_major_axis = if one_minus_x2 == 0.0 {
        f64::INFINITY
    } else {
        semi_perimeter / (2.0 * one_minus_x2)
    };
    Some(LambertSolution {
        v1,
        v2,
        semi_major_axis,
        iterations,
    })
}

/// Value and first three derivatives of the non-dimensional time of flight.
#[derive(Debug, Clone, Copy)]
struct TofEval {
    time: f64,
    d1: f64,
    d2: f64,
    d3: f64,
}

impl TofEval {
    const INVALID: Self = Self {
        time: f64::NAN,
        d1: f64::NAN,
        d2: f64::NAN,
        d3: f64::NAN,
    };
}

/// Non-dimensional time of flight `T(x)` and its derivatives for the
/// zero-revolution transfer with parameter `lambda`.
///
/// Not defined for `x <= -1`, where NaNs are returned.
fn tof_and_derivatives(x: f64, lambda: f64) -> TofEval {
    if !(x > -1.0 && x.is_finite()) {
        return TofEval::INVALID;
    }
    if x > BATTIN_LOWER && x < BATTIN_UPPER {
        tof_battin(x, lambda)
    } else {
        tof_lagrange(x, lambda)
    }
}

/// Lagrange's form of the time of flight with its closed-form derivatives.
///
/// Valid away from the parabolic point, where the derivative formulas are `0/0`.
fn tof_lagrange(x: f64, lambda: f64) -> TofEval {
    let lambda_sq = lambda * lambda;
    let one_minus_x2 = 1.0 - x * x;
    let y = (1.0 - lambda_sq * one_minus_x2).sqrt();
    // η = y - λx without cancellation when λ and x have the same sign
    // (y² - λ²x² = 1 - λ²).
    let eta = if x * lambda > 0.0 {
        (1.0 - lambda_sq) / (y + x * lambda)
    } else {
        y - x * lambda
    };
    let psi = if x < 1.0 {
        (eta * one_minus_x2.sqrt()).atan2(x * y + lambda * one_minus_x2)
    } else {
        (eta * (-one_minus_x2).sqrt()).asinh()
    };
    let time = (psi / one_minus_x2.abs().sqrt() - x + lambda * y) / one_minus_x2;
    let lambda3 = lambda_sq * lambda;
    let d1 = (3.0 * time * x - 2.0 + 2.0 * lambda3 * x / y) / one_minus_x2;
    let d2 =
        (3.0 * time + 5.0 * x * d1 + 2.0 * (1.0 - lambda_sq) * lambda3 / y.powi(3)) / one_minus_x2;
    let d3 = (7.0 * x * d2 + 8.0 * d1
        - 6.0 * (1.0 - lambda_sq) * lambda3 * lambda_sq * x / y.powi(5))
        / one_minus_x2;
    TofEval { time, d1, d2, d3 }
}

/// Battin's hypergeometric form of the time of flight, differentiated with
/// truncated Taylor arithmetic.  Regular at the parabolic point `x = 1`.
fn tof_battin(x: f64, lambda: f64) -> TofEval {
    let lambda_sq = lambda * lambda;
    let x_jet = Taylor3::variable(x);
    let y = (Taylor3::constant(1.0 - lambda_sq) + (x_jet * x_jet).scale(lambda_sq)).sqrt();
    let eta = y - x_jet.scale(lambda);
    let s = (Taylor3::constant(1.0 - lambda) - x_jet * eta).scale(0.5);
    let hyp = hypergeometric_3_1_five_halves(s);
    let time_jet = (eta * eta * eta * hyp).scale(2.0 / 3.0) + eta.scale(2.0 * lambda);
    let [t0, t1, t2, t3] = time_jet.0;
    TofEval {
        time: t0,
        d1: t1,
        d2: 2.0 * t2,
        d3: 6.0 * t3,
    }
}

/// Truncated Taylor series `c0 + c1 h + c2 h² + c3 h³` of a function of the
/// expansion variable `h`, closed under the operations needed above.
#[derive(Debug, Clone, Copy)]
struct Taylor3([f64; 4]);

impl Taylor3 {
    /// The constant function `value`.
    fn constant(value: f64) -> Self {
        Self([value, 0.0, 0.0, 0.0])
    }

    /// The identity function `h ↦ at + h`.
    fn variable(at: f64) -> Self {
        Self([at, 1.0, 0.0, 0.0])
    }

    /// Multiply by a scalar.
    fn scale(self, factor: f64) -> Self {
        Self(self.0.map(|c| c * factor))
    }

    /// Square root, for a positive constant term.
    fn sqrt(self) -> Self {
        let [a0, a1, a2, a3] = self.0;
        let b0 = a0.sqrt();
        let b1 = a1 / (2.0 * b0);
        let b2 = (a2 - b1 * b1) / (2.0 * b0);
        let b3 = (a3 - 2.0 * b1 * b2) / (2.0 * b0);
        Self([b0, b1, b2, b3])
    }
}

impl Add for Taylor3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let [a0, a1, a2, a3] = self.0;
        let [b0, b1, b2, b3] = rhs.0;
        Self([a0 + b0, a1 + b1, a2 + b2, a3 + b3])
    }
}

impl Sub for Taylor3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let [a0, a1, a2, a3] = self.0;
        let [b0, b1, b2, b3] = rhs.0;
        Self([a0 - b0, a1 - b1, a2 - b2, a3 - b3])
    }
}

impl Mul for Taylor3 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let [a0, a1, a2, a3] = self.0;
        let [b0, b1, b2, b3] = rhs.0;
        Self([
            a0 * b0,
            a0 * b1 + a1 * b0,
            a0 * b2 + a1 * b1 + a2 * b0,
            a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0,
        ])
    }
}

/// `C(k, m)` for `m <= 3`.
fn binomial_up_to_cubic(k: usize, m: usize) -> f64 {
    let k = k as f64;
    match m {
        0 => 1.0,
        1 => k,
        2 => k * (k - 1.0) / 2.0,
        _ => k * (k - 1.0) * (k - 2.0) / 6.0,
    }
}

/// Composition of the Gauss hypergeometric function `₂F₁(3, 1; 5/2; ·)` with
/// the Taylor series `s`, whose constant term must lie well inside `(-1, 1)`.
///
/// `₂F₁(3, 1; 5/2; z) = Σ a_k z^k` with `a_0 = 1` and
/// `a_{k+1} = a_k (3 + k)/(5/2 + k)`.  Expanding it about `z = s0` gives the
/// coefficients `F⁽ᵐ⁾(s0)/m! = Σ_{k≥m} a_k C(k, m) s0^(k-m)` that multiply the
/// powers of `s - s0`.
fn hypergeometric_3_1_five_halves(s: Taylor3) -> Taylor3 {
    const MAX_TERMS: usize = 400;
    let s0 = s.0[0];
    let mut coefficients = [0.0_f64; 4];
    for (m, coefficient) in coefficients.iter_mut().enumerate() {
        let mut a_k = 1.0_f64;
        let mut power = 1.0_f64;
        let mut sum = 0.0_f64;
        for k in 0..MAX_TERMS {
            if k >= m {
                let term = a_k * binomial_up_to_cubic(k, m) * power;
                sum += term;
                power *= s0;
                if term.abs() <= 1e-18 * sum.abs() {
                    break;
                }
            }
            a_k *= (3.0 + k as f64) / (2.5 + k as f64);
        }
        *coefficient = sum;
    }
    let delta = Taylor3([0.0, s.0[1], s.0[2], s.0[3]]);
    let delta_sq = delta * delta;
    let delta_cu = delta_sq * delta;
    Taylor3::constant(coefficients[0])
        + delta.scale(coefficients[1])
        + delta_sq.scale(coefficients[2])
        + delta_cu.scale(coefficients[3])
}

/// Izzo's initial guess for the zero-revolution root of `T(x) = t_target`.
fn initial_guess(lambda: f64, t_target: f64) -> f64 {
    let t_zero = lambda.acos() + lambda * (1.0 - lambda * lambda).sqrt();
    let t_parabolic = 2.0 * (1.0 - lambda.powi(3)) / 3.0;
    if t_target >= t_zero {
        (t_zero / t_target).powf(2.0 / 3.0) - 1.0
    } else if t_target < t_parabolic {
        2.5 * t_parabolic * (t_parabolic - t_target) / (t_target * (1.0 - lambda.powi(5))) + 1.0
    } else {
        (t_zero / t_target).powf((t_parabolic / t_zero).log2()) - 1.0
    }
}

/// Solve `T(x) = t_target` for `x ∈ (-1, ∞)`.
///
/// Returns the root and the number of iterations used, or `None` if no root
/// can be bracketed or the iteration does not converge.
///
/// `T` decreases strictly from `+∞` at `x = -1` to zero, so every iterate
/// narrows a bracket `(lo, hi)`; a Householder step that is not finite or
/// leaves the bracket is replaced by bisection (or by moving upward when only
/// the lower end is known).
fn solve_for_x(lambda: f64, t_target: f64) -> Option<(f64, usize)> {
    if !(t_target.is_finite() && t_target > 0.0 && lambda.abs() < 1.0) {
        return None;
    }
    let mut x = initial_guess(lambda, t_target);
    if !(x.is_finite() && x > -1.0) {
        return None;
    }
    let scale = t_target.max(1.0);
    let mut lo = -1.0_f64;
    let mut hi = f64::INFINITY;
    let mut used = 0;
    for iteration in 1..=MAX_ITERATIONS {
        used = iteration;
        let eval = tof_and_derivatives(x, lambda);
        let residual = eval.time - t_target;
        if !residual.is_finite() {
            return None;
        }
        if residual.abs() <= CONVERGENCE_TOL * scale {
            return Some((x, iteration));
        }
        if residual > 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let TofEval { d1, d2, d3, .. } = eval;
        let denominator = d1 * (d1 * d1 - residual * d2) + d3 * residual * residual / 6.0;
        let step = residual * (d1 * d1 - 0.5 * residual * d2) / denominator;
        let mut next = x - step;
        if !(next.is_finite() && next > lo && next < hi) {
            next = if hi.is_finite() {
                0.5 * (lo + hi)
            } else {
                lo + 1.0 + lo.abs()
            };
        }
        if next == x {
            break;
        }
        x = next;
    }
    let residual = tof_and_derivatives(x, lambda).time - t_target;
    (residual.abs() <= ACCEPTANCE_TOL * scale).then_some((x, used))
}
