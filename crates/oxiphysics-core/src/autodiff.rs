// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Forward-mode automatic differentiation using dual numbers.
//!
//! A [`Dual`] number `(v, dv)` carries a value and its derivative
//! simultaneously.  All arithmetic operations and transcendental functions
//! apply the chain rule so that the derivative tracks through expressions
//! automatically.
//!
//! # Quick start
//!
//! ```text
//! use oxiphysics_core::autodiff::{dual, grad1};
//! let f = |x: oxiphysics_core::autodiff::Dual| x * x + x;  // f(x) = x² + x
//! let df = grad1(f, 3.0);                                   // f'(3) = 2·3 + 1 = 7
//! assert!((df - 7.0).abs() < 1e-12);
//! ```

use std::ops::{Add, Div, Mul, Neg, Sub};

// ---------------------------------------------------------------------------
// Dual number
// ---------------------------------------------------------------------------

/// A dual number `(v, dv)` representing a value and its first derivative.
///
/// All arithmetic operators and elementary functions implement the standard
/// chain-rule so that derivatives propagate automatically.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual {
    /// Primal (real) value.
    pub v: f64,
    /// Derivative (dual) component.
    pub dv: f64,
}

/// Construct a [`Dual`] from a value and a derivative seed.
///
/// Set `dx = 1.0` to differentiate with respect to this variable,
/// or `dx = 0.0` to treat it as a constant.
pub fn dual(x: f64, dx: f64) -> Dual {
    Dual { v: x, dv: dx }
}

impl Dual {
    /// Create a dual number representing the constant `c` (derivative = 0).
    #[allow(dead_code)]
    pub fn constant(c: f64) -> Self {
        Dual { v: c, dv: 0.0 }
    }

    /// Create a dual number representing the variable `x` with seed `dx = 1`.
    #[allow(dead_code)]
    pub fn variable(x: f64) -> Self {
        Dual { v: x, dv: 1.0 }
    }

    /// Sine with chain rule: `sin(u)' = cos(u) * u'`.
    pub fn sin(self) -> Self {
        Dual {
            v: self.v.sin(),
            dv: self.v.cos() * self.dv,
        }
    }

    /// Cosine with chain rule: `cos(u)' = -sin(u) * u'`.
    pub fn cos(self) -> Self {
        Dual {
            v: self.v.cos(),
            dv: -self.v.sin() * self.dv,
        }
    }

    /// Natural exponential with chain rule: `exp(u)' = exp(u) * u'`.
    pub fn exp(self) -> Self {
        let ev = self.v.exp();
        Dual {
            v: ev,
            dv: ev * self.dv,
        }
    }

    /// Natural logarithm with chain rule: `ln(u)' = u'/u`.
    pub fn ln(self) -> Self {
        Dual {
            v: self.v.ln(),
            dv: self.dv / self.v,
        }
    }

    /// Square root with chain rule: `sqrt(u)' = u' / (2 * sqrt(u))`.
    pub fn sqrt(self) -> Self {
        let sv = self.v.sqrt();
        Dual {
            v: sv,
            dv: self.dv / (2.0 * sv),
        }
    }

    /// Absolute value with (sub)derivative: `|u|' = sign(u) * u'`.
    pub fn abs(self) -> Self {
        Dual {
            v: self.v.abs(),
            dv: self.v.signum() * self.dv,
        }
    }

    /// Integer power with chain rule: `u^n' = n * u^(n-1) * u'`.
    pub fn powi(self, n: i32) -> Self {
        Dual {
            v: self.v.powi(n),
            dv: (n as f64) * self.v.powi(n - 1) * self.dv,
        }
    }

    /// Floating-point power: `u^p' = p * u^(p-1) * u'`.
    pub fn powf(self, p: f64) -> Self {
        Dual {
            v: self.v.powf(p),
            dv: p * self.v.powf(p - 1.0) * self.dv,
        }
    }

    /// Tangent: `tan(u)' = u' / cos²(u)`.
    #[allow(dead_code)]
    pub fn tan(self) -> Self {
        let c = self.v.cos();
        Dual {
            v: self.v.tan(),
            dv: self.dv / (c * c),
        }
    }

    /// Hyperbolic sine: `sinh(u)' = cosh(u) * u'`.
    #[allow(dead_code)]
    pub fn sinh(self) -> Self {
        Dual {
            v: self.v.sinh(),
            dv: self.v.cosh() * self.dv,
        }
    }

    /// Hyperbolic cosine: `cosh(u)' = sinh(u) * u'`.
    #[allow(dead_code)]
    pub fn cosh(self) -> Self {
        Dual {
            v: self.v.cosh(),
            dv: self.v.sinh() * self.dv,
        }
    }
}

// ---------------------------------------------------------------------------
// Arithmetic operators
// ---------------------------------------------------------------------------

impl Add for Dual {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Dual {
            v: self.v + rhs.v,
            dv: self.dv + rhs.dv,
        }
    }
}

impl Add<f64> for Dual {
    type Output = Self;
    fn add(self, rhs: f64) -> Self {
        Dual {
            v: self.v + rhs,
            dv: self.dv,
        }
    }
}

impl Add<Dual> for f64 {
    type Output = Dual;
    fn add(self, rhs: Dual) -> Dual {
        Dual {
            v: self + rhs.v,
            dv: rhs.dv,
        }
    }
}

impl Sub for Dual {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Dual {
            v: self.v - rhs.v,
            dv: self.dv - rhs.dv,
        }
    }
}

impl Sub<f64> for Dual {
    type Output = Self;
    fn sub(self, rhs: f64) -> Self {
        Dual {
            v: self.v - rhs,
            dv: self.dv,
        }
    }
}

impl Sub<Dual> for f64 {
    type Output = Dual;
    fn sub(self, rhs: Dual) -> Dual {
        Dual {
            v: self - rhs.v,
            dv: -rhs.dv,
        }
    }
}

impl Mul for Dual {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Dual {
            v: self.v * rhs.v,
            dv: self.dv * rhs.v + self.v * rhs.dv,
        }
    }
}

impl Mul<f64> for Dual {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Dual {
            v: self.v * rhs,
            dv: self.dv * rhs,
        }
    }
}

impl Mul<Dual> for f64 {
    type Output = Dual;
    fn mul(self, rhs: Dual) -> Dual {
        Dual {
            v: self * rhs.v,
            dv: self * rhs.dv,
        }
    }
}

impl Div for Dual {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        Dual {
            v: self.v / rhs.v,
            dv: (self.dv * rhs.v - self.v * rhs.dv) / (rhs.v * rhs.v),
        }
    }
}

impl Div<f64> for Dual {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Dual {
            v: self.v / rhs,
            dv: self.dv / rhs,
        }
    }
}

impl Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        Dual {
            v: -self.v,
            dv: -self.dv,
        }
    }
}

// ---------------------------------------------------------------------------
// Convenience differentiation utilities
// ---------------------------------------------------------------------------

/// Compute the first derivative of `f` at `x` using forward-mode AD.
///
/// Sets the seed derivative to 1 and evaluates `f(dual(x, 1.0)).dv`.
pub fn grad1(f: impl Fn(Dual) -> Dual, x: f64) -> f64 {
    f(dual(x, 1.0)).dv
}

/// Compute the diagonal of the Hessian of a scalar function `f: Rⁿ → R`.
///
/// The second derivative `∂²f/∂xᵢ²` is approximated by perturbing the
/// `i`-th coordinate with a small step `h` and differencing derivatives:
///
/// ```text
/// d²f/dxᵢ² ≈ (f'(xᵢ + h) - f'(xᵢ - h)) / (2h)
/// ```
///
/// where `f'(xᵢ)` is obtained via forward-mode AD holding all other
/// coordinates constant.
pub fn hessian_diag(f: impl Fn(Dual) -> Dual, xs: &[f64]) -> Vec<f64> {
    let h = 1e-5;
    xs.iter()
        .map(|&xi| {
            let fp = f(dual(xi + h, 1.0)).dv;
            let fm = f(dual(xi - h, 1.0)).dv;
            (fp - fm) / (2.0 * h)
        })
        .collect()
}

/// Compute one row of the Jacobian: `∂f/∂xᵢ` for all `i`.
///
/// `f` maps a slice of dual numbers to a single dual output.
/// For each index `i`, we seed `xs[i]` with derivative 1 and evaluate.
pub fn jacobian_row(f: impl Fn(&[Dual]) -> Dual, xs: &[f64]) -> Vec<f64> {
    let n = xs.len();
    let mut row = Vec::with_capacity(n);
    for i in 0..n {
        let duals: Vec<Dual> = xs
            .iter()
            .enumerate()
            .map(|(j, &x)| dual(x, if j == i { 1.0 } else { 0.0 }))
            .collect();
        row.push(f(&duals).dv);
    }
    row
}

// ---------------------------------------------------------------------------
// DualVec — vector of Dual numbers
// ---------------------------------------------------------------------------

/// A vector of [`Dual`] numbers for multivariate forward-mode differentiation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DualVec {
    /// The dual components.
    pub components: Vec<Dual>,
}

impl DualVec {
    /// Construct from a slice of `(value, derivative)` pairs.
    #[allow(dead_code)]
    pub fn from_pairs(pairs: &[(f64, f64)]) -> Self {
        DualVec {
            components: pairs.iter().map(|&(v, dv)| Dual { v, dv }).collect(),
        }
    }

    /// Construct a variable vector: component `i` has seed 1, others 0.
    #[allow(dead_code)]
    pub fn variable(xs: &[f64], seed_idx: usize) -> Self {
        DualVec {
            components: xs
                .iter()
                .enumerate()
                .map(|(i, &x)| dual(x, if i == seed_idx { 1.0 } else { 0.0 }))
                .collect(),
        }
    }

    /// Length of the vector.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// Returns `true` if the vector is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }
}

/// Compute the full gradient of `f: Rⁿ → R` at `xs`.
///
/// Makes `n` forward passes, one per variable.
pub fn gradient(f: impl Fn(&[Dual]) -> Dual, xs: &[f64]) -> Vec<f64> {
    jacobian_row(f, xs)
}

// ---------------------------------------------------------------------------
// Newton-Raphson step via autodiff
// ---------------------------------------------------------------------------

/// Perform one Newton-Raphson step: `x_new = x - f(x) / f'(x)`.
///
/// Uses forward-mode autodiff to compute `f'(x)` automatically.
///
/// # Panics
///
/// Does **not** panic on zero derivative; instead returns `x` unchanged when
/// `|f'(x)| < 1e-14`.
pub fn newton_step(f: impl Fn(Dual) -> Dual, x: f64) -> f64 {
    let d = f(dual(x, 1.0));
    if d.dv.abs() < 1e-14 {
        x
    } else {
        x - d.v / d.dv
    }
}

// ---------------------------------------------------------------------------
// TaylorExpand
// ---------------------------------------------------------------------------

/// Stores the value of a function and its first three derivatives at a point,
/// enabling Taylor polynomial evaluation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TaylorExpand {
    /// The expansion center `x₀`.
    pub center: f64,
    /// `f(x₀)`.
    pub f0: f64,
    /// `f'(x₀)`.
    pub f1: f64,
    /// `f''(x₀)` (second derivative).
    pub f2: f64,
    /// `f'''(x₀)` (third derivative).
    pub f3: f64,
}

impl TaylorExpand {
    /// Build a [`TaylorExpand`] by numerically estimating derivatives via
    /// forward-mode AD and central finite differences.
    ///
    /// - `f1` — exact via autodiff.
    /// - `f2` — via central difference of `f'`.
    /// - `f3` — via central difference of `f''`.
    pub fn build<F>(f: F, x0: f64) -> Self
    where
        F: Fn(Dual) -> Dual + Copy,
    {
        let h = 1e-5;
        let f0 = f(dual(x0, 0.0)).v;
        let f1 = f(dual(x0, 1.0)).dv;
        // Second derivative: (f'(x0+h) - f'(x0-h)) / (2h)
        let f2 = (f(dual(x0 + h, 1.0)).dv - f(dual(x0 - h, 1.0)).dv) / (2.0 * h);
        // Third derivative: (f''(x0+h) - f''(x0-h)) / (2h)
        let f2p = |xv: f64| (f(dual(xv + h, 1.0)).dv - f(dual(xv - h, 1.0)).dv) / (2.0 * h);
        let f3 = (f2p(x0 + h) - f2p(x0 - h)) / (2.0 * h);
        TaylorExpand {
            center: x0,
            f0,
            f1,
            f2,
            f3,
        }
    }

    /// Evaluate the Taylor polynomial at `x`:
    /// `f0 + f1*(x-x0) + f2/2*(x-x0)² + f3/6*(x-x0)³`
    pub fn eval(&self, x: f64) -> f64 {
        let dx = x - self.center;
        self.f0 + self.f1 * dx + (self.f2 / 2.0) * dx * dx + (self.f3 / 6.0) * dx * dx * dx
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{E, PI};

    const EPS: f64 = 1e-9;
    const LOOSE: f64 = 1e-5;

    // --- Dual arithmetic ---

    #[test]
    fn test_dual_add() {
        let a = dual(3.0, 1.0);
        let b = dual(2.0, 4.0);
        let c = a + b;
        assert!((c.v - 5.0).abs() < EPS);
        assert!((c.dv - 5.0).abs() < EPS);
    }

    #[test]
    fn test_dual_sub() {
        let a = dual(5.0, 2.0);
        let b = dual(3.0, 1.0);
        let c = a - b;
        assert!((c.v - 2.0).abs() < EPS);
        assert!((c.dv - 1.0).abs() < EPS);
    }

    #[test]
    fn test_dual_mul_product_rule() {
        // d/dx [x * x] at x=3 = 2x = 6
        let x = dual(3.0, 1.0);
        let y = x * x;
        assert!((y.v - 9.0).abs() < EPS);
        assert!((y.dv - 6.0).abs() < EPS);
    }

    #[test]
    fn test_dual_div_quotient_rule() {
        // d/dx [x / (x+1)] at x=2 = 1/(x+1)^2 = 1/9
        let x = dual(2.0, 1.0);
        let y = x / (x + 1.0);
        assert!((y.v - 2.0 / 3.0).abs() < EPS);
        assert!((y.dv - 1.0 / 9.0).abs() < EPS);
    }

    #[test]
    fn test_dual_neg() {
        let x = dual(4.0, 1.0);
        let y = -x;
        assert!((y.v + 4.0).abs() < EPS);
        assert!((y.dv + 1.0).abs() < EPS);
    }

    #[test]
    fn test_dual_add_f64() {
        let x = dual(1.0, 1.0) + 5.0;
        assert!((x.v - 6.0).abs() < EPS);
        assert!((x.dv - 1.0).abs() < EPS);
    }

    #[test]
    fn test_dual_mul_f64() {
        let x = dual(3.0, 1.0) * 4.0;
        assert!((x.v - 12.0).abs() < EPS);
        assert!((x.dv - 4.0).abs() < EPS);
    }

    #[test]
    fn test_dual_sub_f64_rhs() {
        let x = 10.0_f64 - dual(3.0, 1.0);
        assert!((x.v - 7.0).abs() < EPS);
        assert!((x.dv + 1.0).abs() < EPS);
    }

    // --- Transcendental functions ---

    #[test]
    fn test_dual_sin_derivative() {
        // d/dx sin(x) at x = π/4 = cos(π/4)
        let x = dual(PI / 4.0, 1.0);
        let y = x.sin();
        assert!((y.v - (PI / 4.0).sin()).abs() < EPS);
        assert!((y.dv - (PI / 4.0).cos()).abs() < EPS);
    }

    #[test]
    fn test_dual_cos_derivative() {
        let x = dual(PI / 3.0, 1.0);
        let y = x.cos();
        assert!((y.v - (PI / 3.0).cos()).abs() < EPS);
        assert!((y.dv + (PI / 3.0).sin()).abs() < EPS);
    }

    #[test]
    fn test_dual_exp_derivative() {
        // d/dx exp(x) = exp(x)
        let x = dual(2.0, 1.0);
        let y = x.exp();
        assert!((y.v - E * E).abs() < 1e-10);
        assert!((y.dv - E * E).abs() < 1e-10);
    }

    #[test]
    fn test_dual_ln_derivative() {
        // d/dx ln(x) = 1/x
        let x = dual(3.0, 1.0);
        let y = x.ln();
        assert!((y.v - 3.0_f64.ln()).abs() < EPS);
        assert!((y.dv - 1.0 / 3.0).abs() < EPS);
    }

    #[test]
    fn test_dual_sqrt_derivative() {
        // d/dx sqrt(x) = 1/(2*sqrt(x))
        let x = dual(4.0, 1.0);
        let y = x.sqrt();
        assert!((y.v - 2.0).abs() < EPS);
        assert!((y.dv - 0.25).abs() < EPS);
    }

    #[test]
    fn test_dual_abs_positive() {
        let x = dual(3.0, 1.0);
        let y = x.abs();
        assert!((y.v - 3.0).abs() < EPS);
        assert!((y.dv - 1.0).abs() < EPS);
    }

    #[test]
    fn test_dual_abs_negative() {
        let x = dual(-3.0, 1.0);
        let y = x.abs();
        assert!((y.v - 3.0).abs() < EPS);
        assert!((y.dv + 1.0).abs() < EPS);
    }

    #[test]
    fn test_dual_powi() {
        // d/dx x^3 at x=2 = 3x^2 = 12
        let x = dual(2.0, 1.0);
        let y = x.powi(3);
        assert!((y.v - 8.0).abs() < EPS);
        assert!((y.dv - 12.0).abs() < EPS);
    }

    #[test]
    fn test_dual_powf() {
        // d/dx x^2.5 at x=4 = 2.5 * 4^1.5 = 2.5 * 8 = 20
        let x = dual(4.0, 1.0);
        let y = x.powf(2.5);
        assert!((y.v - 4.0_f64.powf(2.5)).abs() < 1e-10);
        assert!((y.dv - 2.5 * 4.0_f64.powf(1.5)).abs() < 1e-10);
    }

    // --- grad1 ---

    #[test]
    fn test_grad1_quadratic() {
        // f(x) = x² + 3x - 5 → f'(x) = 2x + 3 → f'(4) = 11
        let f = |x: Dual| x * x + dual(3.0, 0.0) * x - 5.0;
        assert!((grad1(f, 4.0) - 11.0).abs() < EPS);
    }

    #[test]
    fn test_grad1_sin() {
        // f(x) = sin(x) → f'(π/6) = cos(π/6)
        let f = |x: Dual| x.sin();
        let expected = (PI / 6.0).cos();
        assert!((grad1(f, PI / 6.0) - expected).abs() < EPS);
    }

    #[test]
    fn test_grad1_chain_rule() {
        // f(x) = exp(x^2) → f'(x) = 2x·exp(x²) at x=1 → 2e
        let f = |x: Dual| (x * x).exp();
        let expected = 2.0 * E;
        assert!((grad1(f, 1.0) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_grad1_constant_function() {
        let f = |_x: Dual| dual(42.0, 0.0);
        assert!(grad1(f, 1.0).abs() < EPS);
    }

    // --- hessian_diag ---

    #[test]
    fn test_hessian_diag_quadratic() {
        // f(x) = x² → f''(x) = 2 for all x
        let f = |x: Dual| x * x;
        let xs = [1.0, 2.0, -3.0];
        let h = hessian_diag(f, &xs);
        for hi in &h {
            assert!((hi - 2.0).abs() < LOOSE, "expected 2, got {hi}");
        }
    }

    #[test]
    fn test_hessian_diag_sin() {
        // f(x) = sin(x) → f''(x) = -sin(x)
        let f = |x: Dual| x.sin();
        let xs = [PI / 4.0];
        let h = hessian_diag(f, &xs);
        let expected = -(PI / 4.0).sin();
        assert!((h[0] - expected).abs() < LOOSE);
    }

    // --- jacobian_row / gradient ---

    #[test]
    fn test_jacobian_row_linear() {
        // f(x,y) = 2x + 3y → ∇f = [2, 3]
        let f = |xs: &[Dual]| dual(2.0, 0.0) * xs[0] + dual(3.0, 0.0) * xs[1];
        let row = jacobian_row(f, &[1.0, 1.0]);
        assert!((row[0] - 2.0).abs() < EPS);
        assert!((row[1] - 3.0).abs() < EPS);
    }

    #[test]
    fn test_gradient_quadratic_surface() {
        // f(x,y) = x² + y² → ∇f(3,4) = [6, 8]
        let f = |xs: &[Dual]| xs[0] * xs[0] + xs[1] * xs[1];
        let g = gradient(f, &[3.0, 4.0]);
        assert!((g[0] - 6.0).abs() < EPS);
        assert!((g[1] - 8.0).abs() < EPS);
    }

    #[test]
    fn test_gradient_cross_term() {
        // f(x,y,z) = x*y*z → ∂f/∂x = y*z, etc. at (2,3,4)
        let f = |xs: &[Dual]| xs[0] * xs[1] * xs[2];
        let g = gradient(f, &[2.0, 3.0, 4.0]);
        assert!((g[0] - 12.0).abs() < EPS); // y*z = 12
        assert!((g[1] - 8.0).abs() < EPS); // x*z = 8
        assert!((g[2] - 6.0).abs() < EPS); // x*y = 6
    }

    // --- newton_step ---

    #[test]
    fn test_newton_step_sqrt2() {
        // f(x) = x² - 2, root = √2 ≈ 1.4142
        let f = |x: Dual| x * x - 2.0;
        let mut x = 2.0_f64;
        for _ in 0..20 {
            x = newton_step(f, x);
        }
        assert!((x - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_newton_step_cube_root() {
        // f(x) = x³ - 8, root = 2
        let f = |x: Dual| x.powi(3) - 8.0;
        let mut x = 3.0_f64;
        for _ in 0..30 {
            x = newton_step(f, x);
        }
        assert!((x - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_newton_step_zero_derivative_safe() {
        // At the extremum of f(x) = x², f'(0) = 0; must not blow up.
        let f = |x: Dual| x * x;
        let x_new = newton_step(f, 0.0);
        assert!(x_new.is_finite());
    }

    // --- DualVec ---

    #[test]
    fn test_dualvec_len() {
        let dv = DualVec::variable(&[1.0, 2.0, 3.0], 0);
        assert_eq!(dv.len(), 3);
    }

    #[test]
    fn test_dualvec_seed() {
        let dv = DualVec::variable(&[1.0, 2.0, 3.0], 1);
        assert!((dv.components[0].dv).abs() < EPS);
        assert!((dv.components[1].dv - 1.0).abs() < EPS);
        assert!((dv.components[2].dv).abs() < EPS);
    }

    #[test]
    fn test_dualvec_from_pairs() {
        let dv = DualVec::from_pairs(&[(1.0, 0.0), (2.0, 1.0)]);
        assert_eq!(dv.len(), 2);
        assert!((dv.components[1].dv - 1.0).abs() < EPS);
    }

    // --- TaylorExpand ---

    #[test]
    fn test_taylor_expand_sin_at_zero() {
        // sin expanded at 0: f0=0, f1=1, f2=0, f3=-1
        // T3(x) = x - x³/6
        let t = TaylorExpand::build(|x| x.sin(), 0.0);
        let x = 0.3;
        let approx = t.eval(x);
        let exact = x.sin();
        assert!(
            (approx - exact).abs() < 1e-4,
            "approx={approx}, exact={exact}"
        );
    }

    #[test]
    fn test_taylor_expand_exp_at_zero() {
        // exp expanded at 0: T3(x) = 1 + x + x²/2 + x³/6
        let t = TaylorExpand::build(|x| x.exp(), 0.0);
        let x = 0.5;
        let approx = t.eval(x);
        let exact = x.exp();
        assert!(
            (approx - exact).abs() < 1e-2,
            "approx={approx}, exact={exact}"
        );
    }

    #[test]
    fn test_taylor_expand_exact_at_center() {
        let t = TaylorExpand::build(|x| x * x + x, 2.0);
        // At x = center, eval should equal f(center) = 6
        let approx = t.eval(2.0);
        assert!((approx - 6.0).abs() < 1e-9, "approx={approx}");
    }

    #[test]
    fn test_taylor_expand_f1_matches_grad1() {
        let f = |x: Dual| (x * x).sin();
        let x0 = 1.0;
        let t = TaylorExpand::build(f, x0);
        let g = grad1(f, x0);
        assert!((t.f1 - g).abs() < 1e-10);
    }

    // --- Dual::constant / variable helpers ---

    #[test]
    fn test_dual_constant_zero_derivative() {
        let c = Dual::constant(7.0);
        assert!((c.v - 7.0).abs() < EPS);
        assert!(c.dv.abs() < EPS);
    }

    #[test]
    fn test_dual_variable_unit_derivative() {
        let v = Dual::variable(5.0);
        assert!((v.v - 5.0).abs() < EPS);
        assert!((v.dv - 1.0).abs() < EPS);
    }

    // --- tan / sinh / cosh ---

    #[test]
    fn test_dual_tan_derivative() {
        // d/dx tan(x) at x = π/4 = 1/cos²(π/4) = 2
        let x = dual(PI / 4.0, 1.0);
        let y = x.tan();
        let expected_dv = 1.0 / (PI / 4.0).cos().powi(2);
        assert!((y.dv - expected_dv).abs() < 1e-10);
    }

    #[test]
    fn test_dual_sinh_derivative() {
        // d/dx sinh(x) = cosh(x)
        let x = dual(1.0, 1.0);
        let y = x.sinh();
        assert!((y.dv - 1.0_f64.cosh()).abs() < EPS);
    }

    #[test]
    fn test_dual_cosh_derivative() {
        // d/dx cosh(x) = sinh(x)
        let x = dual(1.0, 1.0);
        let y = x.cosh();
        assert!((y.dv - 1.0_f64.sinh()).abs() < EPS);
    }

    // --- operator overloads for f64 lhs ---

    #[test]
    fn test_f64_add_dual() {
        let d = 3.0_f64 + dual(2.0, 1.0);
        assert!((d.v - 5.0).abs() < EPS);
        assert!((d.dv - 1.0).abs() < EPS);
    }

    #[test]
    fn test_f64_mul_dual() {
        let d = 5.0_f64 * dual(3.0, 1.0);
        assert!((d.v - 15.0).abs() < EPS);
        assert!((d.dv - 5.0).abs() < EPS);
    }

    // --- composed functions ---

    #[test]
    fn test_composed_sin_exp() {
        // f(x) = sin(exp(x)) at x=0: f'(0) = cos(exp(0)) * exp(0) = cos(1)
        let f = |x: Dual| x.exp().sin();
        let expected = 1.0_f64.cos();
        assert!((grad1(f, 0.0) - expected).abs() < EPS);
    }

    #[test]
    fn test_composed_sqrt_ln() {
        // f(x) = sqrt(ln(x)) at x=e: f'(e) = 1/(2*sqrt(ln(e))) * 1/e = 1/(2e)
        let f = |x: Dual| x.ln().sqrt();
        let expected = 1.0 / (2.0 * E);
        assert!((grad1(f, E) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_div_f64() {
        let d = dual(6.0, 3.0) / 2.0;
        assert!((d.v - 3.0).abs() < EPS);
        assert!((d.dv - 1.5).abs() < EPS);
    }
}
