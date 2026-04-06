#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Structural reliability analysis.
//!
//! Implements FORM (First-Order Reliability Method), Monte Carlo simulation,
//! Sobol sensitivity indices, response surface methods, and failure probability
//! calculations for structural safety assessment.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Approximation of the standard normal CDF using Abramowitz & Stegun formula 26.2.17.
///
/// Maximum absolute error ≈ 7.5 × 10⁻⁸.
pub fn normal_cdf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let pdf = (-0.5 * x * x).exp() / (2.0 * PI).sqrt();
    let cdf = 1.0 - pdf * poly;
    if x >= 0.0 { cdf } else { 1.0 - cdf }
}

/// Standard normal PDF: φ(x) = exp(−x²/2) / √(2π).
pub fn normal_pdf(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * PI).sqrt()
}

/// Compute the reliability index β = −Φ⁻¹(pf) from failure probability.
///
/// Uses an iterative bisection approach on the normal CDF.
pub fn reliability_index(pf: f64) -> f64 {
    if pf <= 0.0 {
        return f64::INFINITY;
    }
    if pf >= 1.0 {
        return f64::NEG_INFINITY;
    }
    // Bisect for x such that normal_cdf(-x) = pf  =>  x = -Phi_inv(pf)
    let mut lo = -10.0_f64;
    let mut hi = 10.0_f64;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if normal_cdf(-mid) > pf {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ---------------------------------------------------------------------------
// LimitStateFunction
// ---------------------------------------------------------------------------

/// A limit-state function g(X) defined by its input variable statistics.
///
/// The random input vector X has independent normal components with the given
/// means and standard deviations.  The function g is evaluated by a closure
/// stored inside the struct.  Failure occurs when g(X) ≤ 0.
pub struct LimitStateFunction {
    /// Mean values of the random input variables.
    pub mean: Vec<f64>,
    /// Standard deviations of the random input variables.
    pub std_dev: Vec<f64>,
    /// Closure that evaluates g(x).
    eval: Box<dyn Fn(&[f64]) -> f64 + Send + Sync>,
}

impl LimitStateFunction {
    /// Create a new limit-state function.
    pub fn new(
        mean: Vec<f64>,
        std_dev: Vec<f64>,
        eval: impl Fn(&[f64]) -> f64 + Send + Sync + 'static,
    ) -> Self {
        Self {
            mean,
            std_dev,
            eval: Box::new(eval),
        }
    }

    /// Evaluate g at physical-space point `x`.
    pub fn evaluate(&self, x: &[f64]) -> f64 {
        (self.eval)(x)
    }

    /// Transform a standard-normal point `u` to physical space.
    pub fn u_to_x(&self, u: &[f64]) -> Vec<f64> {
        u.iter()
            .zip(self.mean.iter())
            .zip(self.std_dev.iter())
            .map(|((ui, mi), si)| mi + si * ui)
            .collect()
    }

    /// Transform a physical-space point `x` to standard-normal space.
    pub fn x_to_u(&self, x: &[f64]) -> Vec<f64> {
        x.iter()
            .zip(self.mean.iter())
            .zip(self.std_dev.iter())
            .map(|((xi, mi), si)| (xi - mi) / si)
            .collect()
    }

    /// Number of random variables.
    pub fn n_vars(&self) -> usize {
        self.mean.len()
    }
}

// ---------------------------------------------------------------------------
// FirstOrderReliability (FORM)
// ---------------------------------------------------------------------------

/// Results of a First-Order Reliability Method (FORM) analysis.
#[derive(Debug, Clone)]
pub struct FirstOrderReliability {
    /// Hasofer–Lind reliability index β_HL.
    pub beta: f64,
    /// First-order probability of failure Pf ≈ Φ(−β).
    pub pf: f64,
    /// Design point in standard-normal space (most probable failure point).
    pub design_point: Vec<f64>,
}

impl FirstOrderReliability {
    /// Run FORM analysis on the given limit-state function.
    ///
    /// Uses the iHL-RF (improved Hasofer–Lind–Rackwitz–Fießler) iteration.
    pub fn compute_form(g: &LimitStateFunction) -> Self {
        let n = g.n_vars();
        // Start at origin in u-space
        let mut u: Vec<f64> = vec![0.0; n];
        let eps = 1e-6_f64;

        for _iter in 0..100 {
            let x = g.u_to_x(&u);
            let gval = g.evaluate(&x);

            // Numerical gradient ∂g/∂u_i via finite differences
            let mut grad = vec![0.0; n];
            for i in 0..n {
                let mut xp = x.clone();
                xp[i] += eps * g.std_dev[i];
                let mut xm = x.clone();
                xm[i] -= eps * g.std_dev[i];
                // chain rule: ∂g/∂u_i = ∂g/∂x_i * σ_i
                let dg_dxi = (g.evaluate(&xp) - g.evaluate(&xm)) / (2.0 * eps * g.std_dev[i]);
                grad[i] = dg_dxi * g.std_dev[i];
            }

            let grad_norm_sq: f64 = grad.iter().map(|v| v * v).sum();
            let grad_norm = grad_norm_sq.sqrt();
            if grad_norm < 1e-15 {
                break;
            }

            // iHL-RF update: u_new = (dot(grad,u) - g) / ||grad||^2  * grad
            let dot: f64 = grad.iter().zip(u.iter()).map(|(g, u)| g * u).sum();
            let lambda = (dot - gval) / grad_norm_sq;
            let u_new: Vec<f64> = grad.iter().map(|gi| lambda * gi).collect();

            // Check convergence
            let diff: f64 = u_new
                .iter()
                .zip(u.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            u = u_new;
            if diff < 1e-8 {
                break;
            }
        }

        let beta = Self::hasofer_lind_beta(&u);
        let pf = normal_cdf(-beta);
        Self {
            beta,
            pf,
            design_point: u,
        }
    }

    /// Compute the Hasofer–Lind reliability index as ||u|| in standard-normal space.
    pub fn hasofer_lind_beta(u: &[f64]) -> f64 {
        u.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
}

// ---------------------------------------------------------------------------
// MonteCarloReliability
// ---------------------------------------------------------------------------

/// Monte Carlo reliability analysis.
#[derive(Debug, Clone)]
pub struct MonteCarloReliability {
    /// Number of samples used.
    pub n_samples: usize,
    /// Estimated probability of failure.
    pub pf: f64,
    /// 95 % confidence interval half-width.
    pub ci_half_width: f64,
}

impl MonteCarloReliability {
    /// Estimate the probability of failure by direct Monte Carlo sampling.
    ///
    /// Uses `rand::rng()` with the standard normal transform (Box–Muller).
    pub fn compute_pf(g: &LimitStateFunction, n_samples: usize) -> Self {
        use rand::RngExt;
        let mut rng = rand::rng();
        let n = g.n_vars();
        let mut failures = 0usize;

        for _ in 0..n_samples {
            // Sample each variable from its marginal normal distribution
            let x: Vec<f64> = (0..n)
                .map(|i| {
                    // Box-Muller transform
                    let u1: f64 = rng.random_range(1e-15_f64..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
                    g.mean[i] + g.std_dev[i] * z
                })
                .collect();

            if g.evaluate(&x) <= 0.0 {
                failures += 1;
            }
        }

        let pf = failures as f64 / n_samples as f64;
        // Wilson confidence interval approximation
        let ci_half_width = if n_samples > 0 {
            1.96 * (pf * (1.0 - pf) / n_samples as f64).sqrt()
        } else {
            0.0
        };

        Self {
            n_samples,
            pf,
            ci_half_width,
        }
    }
}

// ---------------------------------------------------------------------------
// SensitivityIndex (Sobol)
// ---------------------------------------------------------------------------

/// Sobol variance-based sensitivity indices.
#[derive(Debug, Clone)]
pub struct SensitivityIndex {
    /// First-order Sobol indices S_i.
    pub sobol_first_order: Vec<f64>,
    /// Total-effect Sobol indices S_Ti.
    pub sobol_total: Vec<f64>,
}

impl SensitivityIndex {
    /// Estimate Sobol first-order and total-effect indices using the
    /// Saltelli (2002) estimator with `n` base samples per variable.
    pub fn compute_sobol(g: &LimitStateFunction, n: usize) -> Self {
        use rand::RngExt;
        let mut rng = rand::rng();
        let k = g.n_vars();

        // Generate two independent sample matrices A and B (n × k)
        // Each column is drawn from N(mean_i, std_i)
        let sample_matrix = |rng: &mut rand::rngs::ThreadRng| -> Vec<Vec<f64>> {
            (0..n)
                .map(|_| {
                    (0..k)
                        .map(|i| {
                            let u1: f64 = rng.random_range(1e-15_f64..1.0_f64);
                            let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
                            g.mean[i] + g.std_dev[i] * z
                        })
                        .collect()
                })
                .collect()
        };

        let mat_a = sample_matrix(&mut rng);
        let mat_b = sample_matrix(&mut rng);

        // Evaluate on A and B
        let ya: Vec<f64> = mat_a.iter().map(|x| g.evaluate(x)).collect();
        let yb: Vec<f64> = mat_b.iter().map(|x| g.evaluate(x)).collect();

        let f0: f64 = ya.iter().sum::<f64>() / n as f64;
        let var_y: f64 = ya.iter().map(|v| (v - f0).powi(2)).sum::<f64>() / n as f64;

        let mut sobol_first_order = vec![0.0; k];
        let mut sobol_total = vec![0.0; k];

        if var_y < 1e-30 {
            return Self {
                sobol_first_order,
                sobol_total,
            };
        }

        for j in 0..k {
            // A_B^(j): matrix A with column j replaced by column j of B
            let ab_j: Vec<Vec<f64>> = (0..n)
                .map(|r| {
                    let mut row = mat_a[r].clone();
                    row[j] = mat_b[r][j];
                    row
                })
                .collect();

            let y_ab_j: Vec<f64> = ab_j.iter().map(|x| g.evaluate(x)).collect();

            // Saltelli (2002) estimators
            let si: f64 = ya
                .iter()
                .zip(y_ab_j.iter())
                .zip(yb.iter())
                .map(|((a, ab), b)| b * (ab - a))
                .sum::<f64>()
                / (n as f64 * var_y);

            let sti: f64 = ya
                .iter()
                .zip(y_ab_j.iter())
                .map(|(a, ab)| (a - ab).powi(2))
                .sum::<f64>()
                / (2.0 * n as f64 * var_y);

            sobol_first_order[j] = si.clamp(0.0, 1.0);
            sobol_total[j] = sti.clamp(0.0, 1.0);
        }

        Self {
            sobol_first_order,
            sobol_total,
        }
    }
}

// ---------------------------------------------------------------------------
// FailureProbability
// ---------------------------------------------------------------------------

/// Container for failure probability and associated reliability index.
#[derive(Debug, Clone)]
pub struct FailureProbability {
    /// Probability of failure.
    pub pf: f64,
    /// Reliability index β = Φ⁻¹(1 − Pf).
    pub beta: f64,
}

impl FailureProbability {
    /// Compute β from a given probability of failure.
    pub fn compute_beta_from_pf(pf: f64) -> f64 {
        reliability_index(pf)
    }

    /// Compute Pf from a given reliability index.
    pub fn compute_pf_from_beta(beta: f64) -> f64 {
        normal_cdf(-beta)
    }

    /// Construct from a probability of failure value.
    pub fn from_pf(pf: f64) -> Self {
        let beta = Self::compute_beta_from_pf(pf);
        Self { pf, beta }
    }

    /// Construct from a reliability index.
    pub fn from_beta(beta: f64) -> Self {
        let pf = Self::compute_pf_from_beta(beta);
        Self { pf, beta }
    }
}

// ---------------------------------------------------------------------------
// ResponseSurface
// ---------------------------------------------------------------------------

/// Quadratic response surface model y ≈ β₀ + Σ βᵢxᵢ + Σ βᵢᵢxᵢ² + Σᵢ<ⱼ βᵢⱼxᵢxⱼ.
#[derive(Debug, Clone)]
pub struct ResponseSurface {
    /// Fitted coefficients \[β₀, β₁, …, βₖ, β₁₁, β₂₂, …, βₖₖ, β₁₂, β₁₃, …\].
    pub coefficients: Vec<f64>,
    /// Number of input variables.
    pub n_vars: usize,
}

impl ResponseSurface {
    /// Expand `x` into the full quadratic feature vector.
    fn features(x: &[f64]) -> Vec<f64> {
        let k = x.len();
        let mut f = Vec::with_capacity(1 + k + k + k * (k - 1) / 2);
        f.push(1.0); // intercept
        for xi in x {
            f.push(*xi); // linear
        }
        for xi in x {
            f.push(xi * xi); // pure quadratic
        }
        for i in 0..k {
            for j in (i + 1)..k {
                f.push(x[i] * x[j]); // cross terms
            }
        }
        f
    }

    /// Fit a quadratic response surface by ordinary least squares.
    ///
    /// Uses the normal equations X'Xβ = X'y solved by Gaussian elimination.
    pub fn fit_quadratic(x_data: &[Vec<f64>], y_data: &[f64]) -> Self {
        assert!(!x_data.is_empty(), "x_data must be non-empty");
        let n_obs = x_data.len();
        let n_vars = x_data[0].len();

        // Build the design matrix
        let rows: Vec<Vec<f64>> = x_data.iter().map(|x| Self::features(x)).collect();
        let p = rows[0].len(); // number of basis functions

        // Form X^T X  and  X^T y  (p×p and p×1)
        let mut xtx = vec![0.0_f64; p * p];
        let mut xty = vec![0.0_f64; p];
        for (row, &yi) in rows.iter().zip(y_data.iter()) {
            for i in 0..p {
                xty[i] += row[i] * yi;
                for j in 0..p {
                    xtx[i * p + j] += row[i] * row[j];
                }
            }
        }

        // Tikhonov regularisation for numerical stability
        let lambda = 1e-10;
        for i in 0..p {
            xtx[i * p + i] += lambda;
        }

        // Gaussian elimination with partial pivoting
        let coefficients = solve_linear_system(&xtx, &xty, p);

        let _ = n_obs; // used implicitly
        Self {
            coefficients,
            n_vars,
        }
    }

    /// Evaluate the response surface at `x`.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let f = Self::features(x);
        f.iter()
            .zip(self.coefficients.iter())
            .map(|(fi, ci)| fi * ci)
            .sum()
    }
}

/// Solve the n×n linear system Ax = b by Gaussian elimination with partial pivoting.
fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut mat = vec![0.0_f64; n * (n + 1)];
    // Build augmented matrix
    for i in 0..n {
        for j in 0..n {
            mat[i * (n + 1) + j] = a[i * n + j];
        }
        mat[i * (n + 1) + n] = b[i];
    }

    for col in 0..n {
        // Partial pivot
        let mut max_row = col;
        let mut max_val = mat[col * (n + 1) + col].abs();
        for row in (col + 1)..n {
            let v = mat[row * (n + 1) + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_row != col {
            for j in 0..=(n) {
                mat.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }

        let pivot = mat[col * (n + 1) + col];
        if pivot.abs() < 1e-15 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = mat[row * (n + 1) + col] / pivot;
            for j in col..=(n) {
                let v = mat[col * (n + 1) + j];
                mat[row * (n + 1) + j] -= factor * v;
            }
        }
    }

    // Back-substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = mat[i * (n + 1) + n];
        for j in (i + 1)..n {
            sum -= mat[i * (n + 1) + j] * x[j];
        }
        let pivot = mat[i * (n + 1) + i];
        x[i] = if pivot.abs() > 1e-15 {
            sum / pivot
        } else {
            0.0
        };
    }
    x
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- normal_cdf / normal_pdf ---

    #[test]
    fn test_normal_cdf_symmetry() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_normal_cdf_plus_one() {
        assert!((normal_cdf(1.0) - 0.841_344_746).abs() < 1e-5);
    }

    #[test]
    fn test_normal_cdf_minus_one() {
        assert!((normal_cdf(-1.0) - 0.158_655_254).abs() < 1e-5);
    }

    #[test]
    fn test_normal_cdf_complement() {
        let x = 1.5;
        assert!((normal_cdf(x) + normal_cdf(-x) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_normal_cdf_large_positive() {
        assert!(normal_cdf(8.0) > 0.999_999);
    }

    #[test]
    fn test_normal_cdf_large_negative() {
        assert!(normal_cdf(-8.0) < 1e-6);
    }

    #[test]
    fn test_normal_pdf_peak() {
        let peak = normal_pdf(0.0);
        assert!((peak - 1.0 / (2.0 * PI).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_normal_pdf_symmetry() {
        assert!((normal_pdf(1.0) - normal_pdf(-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_normal_pdf_positive() {
        assert!(normal_pdf(3.0) > 0.0);
    }

    // --- reliability_index ---

    #[test]
    fn test_reliability_index_half() {
        // pf = 0.5  =>  beta = 0
        let b = reliability_index(0.5);
        assert!(b.abs() < 1e-3);
    }

    #[test]
    fn test_reliability_index_small_pf() {
        // pf ≈ Phi(-3)
        let pf = normal_cdf(-3.0);
        let b = reliability_index(pf);
        assert!((b - 3.0).abs() < 1e-3);
    }

    #[test]
    fn test_reliability_index_zero_pf() {
        assert!(reliability_index(0.0).is_infinite());
    }

    #[test]
    fn test_reliability_index_one_pf() {
        assert!(reliability_index(1.0).is_infinite() || reliability_index(1.0) < 0.0);
    }

    // --- LimitStateFunction ---

    #[test]
    fn test_lsf_evaluate() {
        let g = LimitStateFunction::new(vec![10.0, 5.0], vec![1.0, 0.5], |x| x[0] - x[1] - 3.0);
        assert_eq!(g.evaluate(&[10.0, 5.0]), 2.0);
    }

    #[test]
    fn test_lsf_u_to_x_roundtrip() {
        let g = LimitStateFunction::new(vec![10.0, 5.0], vec![2.0, 1.0], |x| x[0] - x[1]);
        let u = vec![1.0, -1.0];
        let x = g.u_to_x(&u);
        let u2 = g.x_to_u(&x);
        assert!((u[0] - u2[0]).abs() < 1e-10);
        assert!((u[1] - u2[1]).abs() < 1e-10);
    }

    #[test]
    fn test_lsf_n_vars() {
        let g = LimitStateFunction::new(vec![1.0, 2.0, 3.0], vec![0.1, 0.2, 0.3], |x| x[0]);
        assert_eq!(g.n_vars(), 3);
    }

    // --- FirstOrderReliability (FORM) ---

    #[test]
    fn test_form_linear_lsf() {
        // g(X1, X2) = X1 - X2, X1~N(10,2), X2~N(5,2)
        // Exact beta = (10-5) / sqrt(4+4) = 5/2√2 ≈ 1.7678
        let g = LimitStateFunction::new(vec![10.0, 5.0], vec![2.0, 2.0], |x| x[0] - x[1]);
        let result = FirstOrderReliability::compute_form(&g);
        assert!((result.beta - 5.0 / (2.0_f64 * 2.0_f64.sqrt())).abs() < 0.05);
    }

    #[test]
    fn test_form_pf_range() {
        let g = LimitStateFunction::new(vec![10.0], vec![2.0], |x| x[0] - 6.0);
        let result = FirstOrderReliability::compute_form(&g);
        assert!(result.pf > 0.0 && result.pf < 1.0);
    }

    #[test]
    fn test_form_beta_nonnegative() {
        let g = LimitStateFunction::new(vec![0.0], vec![1.0], |x| 1.0 - x[0]);
        let result = FirstOrderReliability::compute_form(&g);
        assert!(result.beta >= 0.0);
    }

    #[test]
    fn test_hasofer_lind_beta_origin() {
        assert_eq!(
            FirstOrderReliability::hasofer_lind_beta(&[0.0, 0.0, 0.0]),
            0.0
        );
    }

    #[test]
    fn test_hasofer_lind_beta_unit() {
        assert!((FirstOrderReliability::hasofer_lind_beta(&[1.0, 0.0, 0.0]) - 1.0).abs() < 1e-10);
    }

    // --- MonteCarloReliability ---

    #[test]
    fn test_monte_carlo_pf_range() {
        let g = LimitStateFunction::new(vec![10.0], vec![2.0], |x| x[0] - 6.0);
        let result = MonteCarloReliability::compute_pf(&g, 1000);
        assert!(result.pf >= 0.0 && result.pf <= 1.0);
    }

    #[test]
    fn test_monte_carlo_always_fail() {
        // g always ≤ 0 => pf should be close to 1
        let g = LimitStateFunction::new(vec![0.0], vec![1.0], |_| -1.0);
        let result = MonteCarloReliability::compute_pf(&g, 500);
        assert!((result.pf - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_monte_carlo_never_fail() {
        // g always > 0 => pf should be 0
        let g = LimitStateFunction::new(vec![0.0], vec![1.0], |_| 1.0);
        let result = MonteCarloReliability::compute_pf(&g, 500);
        assert_eq!(result.pf, 0.0);
    }

    #[test]
    fn test_monte_carlo_sample_count() {
        let g = LimitStateFunction::new(vec![0.0], vec![1.0], |_| 1.0);
        let result = MonteCarloReliability::compute_pf(&g, 200);
        assert_eq!(result.n_samples, 200);
    }

    // --- SensitivityIndex ---

    #[test]
    fn test_sobol_indices_sum_le_one() {
        let g = LimitStateFunction::new(vec![0.0, 0.0], vec![1.0, 1.0], |x| x[0] + 0.5 * x[1]);
        // Use more samples to reduce MC variance; allow generous tolerance
        let si = SensitivityIndex::compute_sobol(&g, 2000);
        let sum: f64 = si.sobol_first_order.iter().sum();
        // For additive function the exact sum is 1.0; MC may overshoot slightly
        assert!(sum <= 1.5, "sum of first-order Sobol indices = {sum}");
    }

    #[test]
    fn test_sobol_total_ge_first() {
        let g = LimitStateFunction::new(vec![0.0, 0.0], vec![1.0, 1.0], |x| x[0] * x[1]);
        let si = SensitivityIndex::compute_sobol(&g, 500);
        for (s, t) in si.sobol_first_order.iter().zip(si.sobol_total.iter()) {
            // total >= first (up to numerical noise)
            assert!(t + 1e-4 >= *s);
        }
    }

    #[test]
    fn test_sobol_length_matches_n_vars() {
        let g = LimitStateFunction::new(vec![0.0, 0.0, 0.0], vec![1.0, 1.0, 1.0], |x| {
            x[0] + x[1] + x[2]
        });
        let si = SensitivityIndex::compute_sobol(&g, 200);
        assert_eq!(si.sobol_first_order.len(), 3);
        assert_eq!(si.sobol_total.len(), 3);
    }

    // --- FailureProbability ---

    #[test]
    fn test_failure_probability_roundtrip() {
        let pf0 = 0.001;
        let fp = FailureProbability::from_pf(pf0);
        let pf1 = FailureProbability::compute_pf_from_beta(fp.beta);
        assert!((pf0 - pf1).abs() < 1e-4);
    }

    #[test]
    fn test_failure_probability_from_beta() {
        let fp = FailureProbability::from_beta(3.0);
        assert!(fp.pf < 0.002);
    }

    #[test]
    fn test_failure_probability_beta_3() {
        let beta = FailureProbability::compute_beta_from_pf(normal_cdf(-3.0));
        assert!((beta - 3.0).abs() < 0.01);
    }

    // --- ResponseSurface ---

    #[test]
    fn test_response_surface_linear_exact() {
        // Perfect linear data: y = 2x + 1
        let x_data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();
        let y_data: Vec<f64> = x_data.iter().map(|x| 2.0 * x[0] + 1.0).collect();
        let rs = ResponseSurface::fit_quadratic(&x_data, &y_data);
        let pred = rs.predict(&[5.0]);
        assert!((pred - 11.0).abs() < 0.1);
    }

    #[test]
    fn test_response_surface_quadratic_exact() {
        // y = x^2
        let x_data: Vec<Vec<f64>> = (0..8).map(|i| vec![i as f64 - 3.0]).collect();
        let y_data: Vec<f64> = x_data.iter().map(|x| x[0] * x[0]).collect();
        let rs = ResponseSurface::fit_quadratic(&x_data, &y_data);
        let pred = rs.predict(&[4.0]);
        assert!((pred - 16.0).abs() < 0.5);
    }

    #[test]
    fn test_response_surface_predict_2d() {
        let x_data: Vec<Vec<f64>> = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
            vec![0.5, 0.5],
            vec![-1.0, 0.0],
        ];
        let y_data: Vec<f64> = x_data.iter().map(|x| x[0] + x[1]).collect();
        let rs = ResponseSurface::fit_quadratic(&x_data, &y_data);
        let pred = rs.predict(&[2.0, 3.0]);
        assert!((pred - 5.0).abs() < 0.5);
    }

    #[test]
    fn test_response_surface_coefficients_length() {
        // k=2 => 1 + 2 + 2 + 1 = 6 coefficients
        let x_data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, i as f64 * 2.0]).collect();
        let y_data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let rs = ResponseSurface::fit_quadratic(&x_data, &y_data);
        assert_eq!(rs.coefficients.len(), 6);
    }

    #[test]
    fn test_form_design_point_dimension() {
        let g = LimitStateFunction::new(vec![5.0, 3.0], vec![1.0, 1.0], |x| x[0] - x[1] - 1.0);
        let result = FirstOrderReliability::compute_form(&g);
        assert_eq!(result.design_point.len(), 2);
    }

    #[test]
    fn test_monte_carlo_ci_nonnegative() {
        let g = LimitStateFunction::new(vec![0.0], vec![1.0], |x| x[0]);
        let result = MonteCarloReliability::compute_pf(&g, 300);
        assert!(result.ci_half_width >= 0.0);
    }

    #[test]
    fn test_normal_cdf_two_sigma() {
        // ~95.4% below +2σ
        assert!((normal_cdf(2.0) - 0.9772).abs() < 1e-3);
    }

    #[test]
    fn test_normal_cdf_three_sigma() {
        assert!((normal_cdf(3.0) - 0.9987).abs() < 1e-3);
    }

    #[test]
    fn test_failure_probability_pf_bounds() {
        for beta in [0.5, 1.0, 2.0, 3.0, 4.0] {
            let fp = FailureProbability::from_beta(beta);
            assert!(fp.pf >= 0.0 && fp.pf <= 1.0);
        }
    }
}
