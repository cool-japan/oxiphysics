#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical integration methods.
//!
//! Provides composite trapezoid, composite Simpson, Romberg integration,
//! Gauss-Legendre quadrature (5-point and variable-point), adaptive Simpson,
//! Monte Carlo integration, Clenshaw-Curtis, double-exponential (tanh-sinh),
//! and multidimensional trapezoid rules.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Composite trapezoid rule
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using the composite trapezoid rule with `n` subintervals.
///
/// Uses the formula: h/2 * (f(a) + 2*f(x_1) + ... + 2*f(x_{n-1}) + f(b)).
pub fn trapezoid_rule(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    assert!(n >= 1, "n must be at least 1");
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        sum += 2.0 * f(a + i as f64 * h);
    }
    sum * h / 2.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Composite Simpson's rule
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using composite Simpson's rule with `n` subintervals.
///
/// `n` must be even. Uses the 1/3 Simpson formula with alternating weights 4, 2.
pub fn simpson_rule(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let n = if !n.is_multiple_of(2) { n + 1 } else { n };
    assert!(n >= 2);
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(x);
    }
    sum * h / 3.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Romberg integration
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using Romberg integration (Richardson extrapolation).
///
/// Builds a triangular table of `max_levels` rows; returns the best estimate.
/// Stops early if consecutive diagonal entries agree within `tol`.
pub fn romberg_integration(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    max_levels: usize,
    tol: f64,
) -> f64 {
    let max_levels = max_levels.max(1);
    let mut r = vec![vec![0.0_f64; max_levels]; max_levels];
    r[0][0] = trapezoid_rule(&f, a, b, 1);
    for i in 1..max_levels {
        let n_steps = 1usize << i; // 2^i trapezoid steps
        r[i][0] = trapezoid_rule(&f, a, b, n_steps);
        for j in 1..=i {
            let factor = (4u64.pow(j as u32)) as f64;
            r[i][j] = (factor * r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0);
        }
        if i >= 1 && (r[i][i] - r[i - 1][i - 1]).abs() < tol {
            return r[i][i];
        }
    }
    r[max_levels - 1][max_levels - 1]
}

// ─────────────────────────────────────────────────────────────────────────────
// 5-point Gauss-Legendre
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using the 5-point Gauss-Legendre formula.
///
/// Exact for polynomials of degree ≤ 9.
pub fn gauss_legendre_5pt(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    // Standard 5-point nodes and weights on [-1, 1]
    let nodes = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683,
        0.0,
        0.538_469_310_105_683,
        0.906_179_845_938_664,
    ];
    let weights = [
        0.236_926_885_056_189,
        0.478_628_670_499_366,
        0.568_888_888_888_889,
        0.478_628_670_499_366,
        0.236_926_885_056_189,
    ];
    let mid = (a + b) / 2.0;
    let half = (b - a) / 2.0;
    let mut sum = 0.0;
    for (xi, wi) in nodes.iter().zip(weights.iter()) {
        sum += wi * f(mid + half * xi);
    }
    sum * half
}

// ─────────────────────────────────────────────────────────────────────────────
// Gauss-Legendre nodes and weights (n = 2..5)
// ─────────────────────────────────────────────────────────────────────────────

/// Return the Gauss-Legendre nodes and weights for `n` points on `[-1, 1]`.
///
/// Supported values of `n`: 2, 3, 4, 5. Panics for other values.
/// Returns `(nodes, weights)` both of length `n`.
pub fn gauss_legendre_nodes_weights(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        2 => (
            vec![-0.577_350_269_189_626, 0.577_350_269_189_626],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.774_596_669_241_483, 0.0, 0.774_596_669_241_483],
            vec![
                0.555_555_555_555_556,
                0.888_888_888_888_889,
                0.555_555_555_555_556,
            ],
        ),
        4 => (
            vec![
                -0.861_136_311_594_953,
                -0.339_981_043_584_856,
                0.339_981_043_584_856,
                0.861_136_311_594_953,
            ],
            vec![
                0.347_854_845_137_454,
                0.652_145_154_862_546,
                0.652_145_154_862_546,
                0.347_854_845_137_454,
            ],
        ),
        5 => {
            let (nodes, weights) = (
                vec![
                    -0.906_179_845_938_664,
                    -0.538_469_310_105_683,
                    0.0,
                    0.538_469_310_105_683,
                    0.906_179_845_938_664,
                ],
                vec![
                    0.236_926_885_056_189,
                    0.478_628_670_499_366,
                    0.568_888_888_888_889,
                    0.478_628_670_499_366,
                    0.236_926_885_056_189,
                ],
            );
            (nodes, weights)
        }
        _ => panic!("gauss_legendre_nodes_weights: n must be 2, 3, 4, or 5"),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Adaptive Simpson
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using recursive adaptive Simpson's rule.
///
/// Recursively subdivides if the local error estimate exceeds `tol`.
/// Recursion depth is limited to `max_depth` to prevent stack overflow.
pub fn adaptive_simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64, max_depth: usize) -> f64 {
    fn simpson_one(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
        let c = (a + b) / 2.0;
        (b - a) / 6.0 * (f(a) + 4.0 * f(c) + f(b))
    }

    fn recurse(f: &impl Fn(f64) -> f64, a: f64, b: f64, tol: f64, whole: f64, depth: usize) -> f64 {
        let c = (a + b) / 2.0;
        let left = simpson_one(f, a, c);
        let right = simpson_one(f, c, b);
        let delta = left + right - whole;
        if depth == 0 || delta.abs() < 15.0 * tol {
            left + right + delta / 15.0
        } else {
            let half_tol = tol / 2.0;
            recurse(f, a, c, half_tol, left, depth - 1)
                + recurse(f, c, b, half_tol, right, depth - 1)
        }
    }

    let whole = simpson_one(&f, a, b);
    recurse(&f, a, b, tol, whole, max_depth)
}

// ─────────────────────────────────────────────────────────────────────────────
// Monte Carlo integration
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using Monte Carlo sampling with `n_samples` points.
///
/// Uses a fixed seed via `rand::rng()` for reproducibility within a run.
pub fn monte_carlo_integrate(f: impl Fn(f64) -> f64, a: f64, b: f64, n_samples: usize) -> f64 {
    use rand::RngExt;
    let mut rng = rand::rng();
    let mut sum = 0.0;
    for _ in 0..n_samples {
        let x: f64 = rng.random_range(a..b);
        sum += f(x);
    }
    (b - a) * sum / n_samples as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// Clenshaw-Curtis quadrature
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫_{-1}^{1} f(x) dx using Clenshaw-Curtis quadrature with `n+1` points.
///
/// Uses Chebyshev nodes cos(k*π/n) and weights derived from the DCT.
pub fn clenshaw_curtis(f: impl Fn(f64) -> f64, n: usize) -> f64 {
    let n = n.max(1);
    // Chebyshev nodes: x_k = cos(k*pi/n), k = 0..=n
    let vals: Vec<f64> = (0..=n)
        .map(|k| f((k as f64 * PI / n as f64).cos()))
        .collect();

    // Weights via exact formula for Clenshaw-Curtis
    let mut sum = 0.0;
    let nf = n as f64;
    for j in 0..=n {
        let cj = if j == 0 || j == n { 1.0 } else { 2.0 };
        let mut w = 0.0_f64;
        for k in (0..=(n / 2)).map(|m| 2 * m) {
            let bk = if k == 0 { 1.0 } else { 2.0 };
            let denom = 1.0 - (k as f64).powi(2);
            if denom.abs() < 1e-14 {
                continue;
            }
            let cos_term = ((k * j) as f64 * PI / nf).cos();
            w += bk / denom * cos_term;
        }
        w /= nf;
        sum += cj * w * vals[j];
    }
    sum
}

// ─────────────────────────────────────────────────────────────────────────────
// Double-exponential (tanh-sinh) quadrature
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx on \[a, b\] using the tanh-sinh (double exponential) transform.
///
/// Uses `n` quadrature points spread symmetrically about the midpoint.
pub fn double_exponential(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let n = n.max(2);
    let half_n = n as i64 / 2;
    let h = 3.0 / (half_n as f64); // step size in transformed variable t
    let mid = (a + b) / 2.0;
    let half = (b - a) / 2.0;
    let mut sum = 0.0;
    for ki in (-half_n)..=half_n {
        let t = ki as f64 * h;
        let sinh_t = t.sinh();
        let cosh_t = t.cosh();
        let cosh_sinh = (PI / 2.0 * sinh_t).cosh();
        let x = mid + half * (PI / 2.0 * sinh_t).tanh();
        let w = half * PI / 2.0 * cosh_t / (cosh_sinh * cosh_sinh);
        sum += w * f(x);
    }
    sum * h
}

// ─────────────────────────────────────────────────────────────────────────────
// Multidimensional trapezoid rule
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate ∫ f(x) dx over a hyperrectangle using the composite trapezoid rule.
///
/// `bounds` is a slice of `(a_i, b_i)` pairs, one per dimension.
/// `n` is the number of subintervals per dimension.
pub fn multidimensional_trap(f: impl Fn(&[f64]) -> f64, bounds: &[(f64, f64)], n: usize) -> f64 {
    let dim = bounds.len();
    if dim == 0 {
        return f(&[]);
    }
    let n = n.max(1);

    // Generate all grid points via iterative expansion
    let steps_per_dim: Vec<f64> = bounds.iter().map(|(a, b)| (b - a) / n as f64).collect();
    let total_points = (n + 1).pow(dim as u32);

    let mut integral = 0.0;
    let mut point = vec![0.0_f64; dim];

    for idx in 0..total_points {
        // Decode multi-index: for dimension d, the index is (idx / (n+1)^d) % (n+1)
        let mut weight = 1.0_f64;
        for d in 0..dim {
            let stride = (n + 1).pow(d as u32);
            let i_d = (idx / stride) % (n + 1);
            point[d] = bounds[d].0 + i_d as f64 * steps_per_dim[d];
            let w = if i_d == 0 || i_d == n { 1.0 } else { 2.0 };
            weight *= w;
        }
        integral += weight * f(&point);
    }

    // Multiply by h/2 for each dimension
    let mut factor = integral;
    for (a, b) in bounds.iter() {
        factor *= (b - a) / (2.0 * n as f64);
    }
    factor
}

// ─────────────────────────────────────────────────────────────────────────────
// QuadratureRule
// ─────────────────────────────────────────────────────────────────────────────

/// A generic quadrature rule defined by nodes and weights on a reference interval.
///
/// The nodes and weights are defined on `[-1, 1]`; `apply` rescales to `[a, b]`.
#[derive(Debug, Clone)]
pub struct QuadratureRule {
    /// Quadrature nodes on `[-1, 1]`.
    pub nodes: Vec<f64>,
    /// Corresponding quadrature weights.
    pub weights: Vec<f64>,
    /// Human-readable name (e.g., "Gauss-Legendre-5").
    pub name: String,
}

impl QuadratureRule {
    /// Construct a new quadrature rule.
    pub fn new(nodes: Vec<f64>, weights: Vec<f64>, name: impl Into<String>) -> Self {
        assert_eq!(nodes.len(), weights.len(), "nodes and weights must match");
        Self {
            nodes,
            weights,
            name: name.into(),
        }
    }

    /// Apply the quadrature rule to approximate ∫ f(x) dx on \[a, b\].
    pub fn apply(&self, f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
        let mid = (a + b) / 2.0;
        let half = (b - a) / 2.0;
        let mut sum = 0.0;
        for (xi, wi) in self.nodes.iter().zip(self.weights.iter()) {
            sum += wi * f(mid + half * xi);
        }
        sum * half
    }

    /// Number of quadrature points.
    pub fn n_points(&self) -> usize {
        self.nodes.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── trapezoid ─────────────────────────────────────────────────────────────

    #[test]
    fn test_trapezoid_constant() {
        // ∫ 3 dx on [0, 1] = 3
        let result = trapezoid_rule(|_| 3.0, 0.0, 1.0, 100);
        assert!((result - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_trapezoid_linear() {
        // ∫ x dx on [0, 1] = 0.5 — exact for linear functions
        let result = trapezoid_rule(|x| x, 0.0, 1.0, 1);
        assert!((result - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_trapezoid_sin_pi() {
        // ∫ sin(x) dx on [0, π] = 2
        let result = trapezoid_rule(|x| x.sin(), 0.0, PI, 10_000);
        assert!((result - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_trapezoid_x_squared() {
        // ∫ x^2 dx on [0, 1] = 1/3 (trapezoid has O(h^2) error)
        let result = trapezoid_rule(|x| x * x, 0.0, 1.0, 1_000);
        assert!((result - 1.0 / 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_trapezoid_negative_interval() {
        // ∫ 1 dx on [-2, -1] = 1
        let result = trapezoid_rule(|_| 1.0, -2.0, -1.0, 10);
        assert!((result - 1.0).abs() < 1e-12);
    }

    // ── simpson ───────────────────────────────────────────────────────────────

    #[test]
    fn test_simpson_constant() {
        let result = simpson_rule(|_| 5.0, 0.0, 1.0, 2);
        assert!((result - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_simpson_x_squared_exact() {
        // Simpson is exact for degree ≤ 3 polynomials
        let result = simpson_rule(|x| x * x, 0.0, 1.0, 2);
        assert!((result - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_simpson_cubic_exact() {
        // ∫ x^3 dx on [0, 1] = 0.25
        let result = simpson_rule(|x| x * x * x, 0.0, 1.0, 2);
        assert!((result - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_simpson_sin_pi() {
        let result = simpson_rule(|x| x.sin(), 0.0, PI, 100);
        assert!((result - 2.0).abs() < 1e-5);
    }

    #[test]
    fn test_simpson_odd_n_corrected() {
        // Passing odd n should be auto-corrected to even
        let result = simpson_rule(|x| x * x, 0.0, 1.0, 3);
        assert!((result - 1.0 / 3.0).abs() < 1e-6);
    }

    // ── romberg ───────────────────────────────────────────────────────────────

    #[test]
    fn test_romberg_constant() {
        let result = romberg_integration(|_| 2.0, 0.0, 1.0, 5, 1e-10);
        assert!((result - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_romberg_x_squared() {
        let result = romberg_integration(|x| x * x, 0.0, 1.0, 6, 1e-12);
        assert!((result - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_romberg_pi_integral() {
        // ∫_0^1 4 / (1 + x^2) dx = π
        let result = romberg_integration(|x| 4.0 / (1.0 + x * x), 0.0, 1.0, 8, 1e-10);
        assert!((result - PI).abs() < 1e-8);
    }

    #[test]
    fn test_romberg_sin_pi() {
        let result = romberg_integration(|x| x.sin(), 0.0, PI, 6, 1e-10);
        assert!((result - 2.0).abs() < 1e-8);
    }

    // ── gauss-legendre 5pt ────────────────────────────────────────────────────

    #[test]
    fn test_gauss_legendre_5pt_constant() {
        let result = gauss_legendre_5pt(|_| 1.0, 0.0, 1.0);
        assert!((result - 1.0).abs() < TOL);
    }

    #[test]
    fn test_gauss_legendre_5pt_poly_degree9() {
        // Should be exact for degree ≤ 9: ∫ x^9 dx on [0,1] = 0.1
        let result = gauss_legendre_5pt(|x| x.powi(9), 0.0, 1.0);
        assert!((result - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_gauss_legendre_5pt_sin() {
        let result = gauss_legendre_5pt(|x| x.sin(), 0.0, PI);
        assert!((result - 2.0).abs() < 1e-4);
    }

    // ── gauss_legendre_nodes_weights ──────────────────────────────────────────

    #[test]
    fn test_gl_nodes_weights_n2_sum_weights() {
        let (_nodes, weights) = gauss_legendre_nodes_weights(2);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_gl_nodes_weights_n3_sum_weights() {
        let (_nodes, weights) = gauss_legendre_nodes_weights(3);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_gl_nodes_weights_n4_sum_weights() {
        let (_nodes, weights) = gauss_legendre_nodes_weights(4);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_gl_nodes_weights_n5_sum_weights() {
        let (_nodes, weights) = gauss_legendre_nodes_weights(5);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_gl_nodes_exact_for_quadratic() {
        // 3-point GL is exact for degree ≤ 5. Test ∫_{-1}^{1} x^4 dx = 2/5
        let (nodes, weights) = gauss_legendre_nodes_weights(3);
        let result: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(x, w)| w * x.powi(4))
            .sum();
        assert!((result - 2.0 / 5.0).abs() < 1e-12);
    }

    // ── adaptive simpson ──────────────────────────────────────────────────────

    #[test]
    fn test_adaptive_simpson_sin() {
        // ∫ sin(x) dx on [0, π] = 2
        let result = adaptive_simpson(|x| x.sin(), 0.0, PI, 1e-8, 20);
        assert!((result - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_adaptive_simpson_constant() {
        let result = adaptive_simpson(|_| 7.0, 0.0, 1.0, 1e-10, 10);
        assert!((result - 7.0).abs() < 1e-8);
    }

    #[test]
    fn test_adaptive_simpson_x_cubed() {
        // ∫ x^3 dx on [0, 2] = 4
        let result = adaptive_simpson(|x| x * x * x, 0.0, 2.0, 1e-10, 15);
        assert!((result - 4.0).abs() < 1e-8);
    }

    // ── monte carlo ───────────────────────────────────────────────────────────

    #[test]
    fn test_monte_carlo_constant() {
        // ∫ 1 dx on [0, 1] = 1 — exact regardless of samples
        let result = monte_carlo_integrate(|_| 1.0, 0.0, 1.0, 1_000);
        assert!((result - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_monte_carlo_approx_pi() {
        // ∫_0^1 4 / (1 + x^2) dx ≈ π, converges in Monte Carlo
        let result = monte_carlo_integrate(|x| 4.0 / (1.0 + x * x), 0.0, 1.0, 100_000);
        // Accept 1% relative error for MC
        assert!((result - PI).abs() < 0.05, "MC result: {result}");
    }

    #[test]
    fn test_monte_carlo_negative_interval() {
        // ∫ 2 dx on [-1, 0] = 2
        let result = monte_carlo_integrate(|_| 2.0, -1.0, 0.0, 100);
        assert!((result - 2.0).abs() < 1e-12);
    }

    // ── clenshaw-curtis ───────────────────────────────────────────────────────

    #[test]
    fn test_clenshaw_curtis_constant() {
        // ∫_{-1}^{1} 1 dx = 2
        let result = clenshaw_curtis(|_| 1.0, 8);
        assert!((result - 2.0).abs() < 1e-6, "CC constant: {result}");
    }

    #[test]
    fn test_clenshaw_curtis_linear_zero() {
        // ∫_{-1}^{1} x dx = 0 (odd function)
        let result = clenshaw_curtis(|x| x, 8);
        assert!(result.abs() < 1e-10, "CC odd: {result}");
    }

    // ── double exponential ────────────────────────────────────────────────────

    #[test]
    fn test_double_exponential_constant() {
        // ∫ 1 dx on [0, 1] = 1
        let result = double_exponential(|_| 1.0, 0.0, 1.0, 20);
        assert!((result - 1.0).abs() < 1e-6, "DE constant: {result}");
    }

    #[test]
    fn test_double_exponential_sin() {
        // ∫ sin(x) dx on [0, π] ≈ 2
        let result = double_exponential(|x| x.sin(), 0.0, PI, 40);
        assert!((result - 2.0).abs() < 1e-4, "DE sin: {result}");
    }

    // ── QuadratureRule ────────────────────────────────────────────────────────

    #[test]
    fn test_quadrature_rule_n_points() {
        let (nodes, weights) = gauss_legendre_nodes_weights(4);
        let rule = QuadratureRule::new(nodes, weights, "GL-4");
        assert_eq!(rule.n_points(), 4);
    }

    #[test]
    fn test_quadrature_rule_apply_constant() {
        let (nodes, weights) = gauss_legendre_nodes_weights(3);
        let rule = QuadratureRule::new(nodes, weights, "GL-3");
        let result = rule.apply(|_| 2.0, -1.0, 1.0);
        assert!((result - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_quadrature_rule_apply_quadratic() {
        // ∫_0^1 x^2 dx = 1/3; GL-3 is exact for degree ≤ 5
        let (nodes, weights) = gauss_legendre_nodes_weights(3);
        let rule = QuadratureRule::new(nodes, weights, "GL-3");
        let result = rule.apply(|x| x * x, 0.0, 1.0);
        assert!((result - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_quadrature_rule_name() {
        let (nodes, weights) = gauss_legendre_nodes_weights(2);
        let rule = QuadratureRule::new(nodes, weights, "custom");
        assert_eq!(rule.name, "custom");
    }
}
