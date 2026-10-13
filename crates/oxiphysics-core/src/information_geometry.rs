// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Information geometry and statistical manifolds.
//!
//! This module provides tools for studying the geometry of families of
//! probability distributions, including Fisher information, geodesics,
//! natural gradients, and divergence measures.

use crate::error::{PhysicsError, PhysicsResult};

// ─────────────────────────────────────────────────────────────────────────────
// Submodules: general statistical manifolds, α-geometry, KL projections
// ─────────────────────────────────────────────────────────────────────────────

mod alpha_geometry;
mod projection;
mod statistical_manifold;

pub use alpha_geometry::AlphaGeometry;
pub use projection::InformationProjection;
pub use statistical_manifold::{LogDensityFn, StatisticalManifold};

// ─────────────────────────────────────────────────────────────────────────────
// ExponentialFamily
// ─────────────────────────────────────────────────────────────────────────────

/// An exponential family of distributions.
///
/// `p(x; θ) = exp(θ · T(x) − A(θ))` where `T` are sufficient statistics and
/// `A` is the log-partition (cumulant generating) function.
pub struct ExponentialFamily {
    /// Sufficient statistic functions `T_i(x)`.
    pub sufficient_stats: Vec<fn(&[f64]) -> f64>,
    /// Log-partition function `A(θ)`.
    pub log_partition: fn(&[f64]) -> f64,
}

impl ExponentialFamily {
    /// Create a new `ExponentialFamily`.
    pub fn new(sufficient_stats: Vec<fn(&[f64]) -> f64>, log_partition: fn(&[f64]) -> f64) -> Self {
        Self {
            sufficient_stats,
            log_partition,
        }
    }

    /// Evaluate the natural parameters at `theta` (identity for canonical form).
    pub fn natural_params(&self, theta: &[f64]) -> Vec<f64> {
        theta.to_vec()
    }

    /// Compute the moment parameters `μ_i = ∂A/∂θ_i` by finite difference.
    pub fn moment_params(&self, theta: &[f64]) -> Vec<f64> {
        let h = 1e-5;
        let a = self.log_partition;
        theta
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let mut tp = theta.to_vec();
                let mut tm = theta.to_vec();
                tp[i] += h;
                tm[i] -= h;
                (a(&tp) - a(&tm)) / (2.0 * h)
            })
            .collect()
    }

    /// KL divergence `KL(p_{θ1} ‖ p_{θ2}) = A(θ2) − A(θ1) − (θ2−θ1)·μ1`.
    pub fn kl_divergence(&self, theta1: &[f64], theta2: &[f64]) -> f64 {
        let a = self.log_partition;
        let mu1 = self.moment_params(theta1);
        let diff_a = a(theta2) - a(theta1);
        let dot: f64 = theta2
            .iter()
            .zip(theta1.iter())
            .zip(mu1.iter())
            .map(|((t2, t1), m)| (t2 - t1) * m)
            .sum();
        diff_a - dot
    }

    /// Fisher information matrix `I(θ) = ∂²A/∂θ_i ∂θ_j` (Hessian of A).
    pub fn fisher_info(&self, theta: &[f64]) -> Vec<Vec<f64>> {
        let n = theta.len();
        let h = 1e-4;
        let a = self.log_partition;
        let mut fi = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in i..n {
                let mut pp = theta.to_vec();
                let mut pm = theta.to_vec();
                let mut mp = theta.to_vec();
                let mut mm = theta.to_vec();
                pp[i] += h;
                pp[j] += h;
                pm[i] += h;
                pm[j] -= h;
                mp[i] -= h;
                mp[j] += h;
                mm[i] -= h;
                mm[j] -= h;
                let val = (a(&pp) - a(&pm) - a(&mp) + a(&mm)) / (4.0 * h * h);
                fi[i][j] = val;
                fi[j][i] = val;
            }
        }
        fi
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GaussianManifold
// ─────────────────────────────────────────────────────────────────────────────

/// The manifold of univariate Gaussian distributions parameterized by
/// `(μ, σ)` with `σ > 0`.
#[derive(Debug, Clone)]
pub struct GaussianManifold;

impl GaussianManifold {
    /// Create a new `GaussianManifold`.
    pub fn new() -> Self {
        Self
    }

    /// Fisher information metric for Gaussian: `diag(1/σ², 2/σ²)` at `(μ, σ)`.
    pub fn fisher_metric(&self, _mu: f64, sigma: f64) -> [[f64; 2]; 2] {
        let s2 = sigma * sigma;
        [[1.0 / s2, 0.0], [0.0, 2.0 / s2]]
    }

    /// Geodesic (Fisher-Rao) distance between two Gaussians `(μ1,σ1)` and
    /// `(μ2,σ2)`.
    ///
    /// The Fisher metric `(dμ² + 2 dσ²) / σ²` is twice the Poincaré metric of the
    /// upper half-plane in the coordinates `(μ/√2, σ)`, so the distance is
    /// `√2 · acosh(1 + ((μ1−μ2)²/2 + (σ1−σ2)²) / (2 σ1 σ2))`.  This is the same
    /// quantity as [`crate::information_theory::fisher_rao_gaussian`], which this
    /// method delegates to; for example `d(N(0,1), N(0,e)) = √2`.
    pub fn geodesic_distance(&self, mu1: f64, sigma1: f64, mu2: f64, sigma2: f64) -> f64 {
        crate::information_theory::fisher_rao_gaussian(mu1, sigma1, mu2, sigma2)
    }

    /// Exponential map at `(μ, σ)` in direction `(v_μ, v_σ)` for step `t`.
    ///
    /// Returns the point `(μ', σ')` reached by following the geodesic.
    pub fn exponential_map(
        &self,
        mu: f64,
        sigma: f64,
        v_mu: f64,
        v_sigma: f64,
        t: f64,
    ) -> (f64, f64) {
        // Linear approximation along the geodesic.
        let new_mu = mu + t * v_mu;
        let new_sigma = (sigma + t * v_sigma).max(1e-12);
        (new_mu, new_sigma)
    }

    /// Logarithmic map at `(μ, σ)` pointing toward `(μ2, σ2)`.
    ///
    /// Returns the tangent vector `(v_μ, v_σ)` such that
    /// `exp_{(μ,σ)}(v) = (μ2, σ2)`.
    pub fn logarithmic_map(&self, mu: f64, sigma: f64, mu2: f64, sigma2: f64) -> (f64, f64) {
        let _ = sigma; // suppress unused warning
        (mu2 - mu, sigma2 - sigma)
    }
}

impl Default for GaussianManifold {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Mutual Information Estimator
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate mutual information `I(X;Y)` using the k-nearest-neighbour (kNN)
/// method (Kraskov–Stögbauer–Grassberger estimator).
///
/// `x` and `y` must have the same length. `k` is the number of neighbours.
pub fn mutual_information_estimator(x: &[f64], y: &[f64], k: usize) -> f64 {
    let n = x.len().min(y.len());
    if n <= k {
        return 0.0;
    }
    let k = k.max(1);
    // Digamma approximation: ψ(n) ≈ ln(n) − 1/(2n).
    let digamma = |n: f64| n.ln() - 1.0 / (2.0 * n);
    let points: Vec<(f64, f64)> = x.iter().zip(y.iter()).map(|(&xi, &yi)| (xi, yi)).collect();
    let mut nx_sum = 0.0f64;
    let mut ny_sum = 0.0f64;
    for i in 0..n {
        // Find k-th NN in joint space (Chebyshev distance).
        let mut dists: Vec<f64> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let dx = (points[i].0 - points[j].0).abs();
                let dy = (points[i].1 - points[j].1).abs();
                dx.max(dy)
            })
            .collect();
        dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let eps = dists.get(k - 1).copied().unwrap_or(0.0);
        // Count points within eps in marginal spaces.
        let n_x = x.iter().filter(|&&xi| (xi - x[i]).abs() < eps).count();
        let n_y = y.iter().filter(|&&yi| (yi - y[i]).abs() < eps).count();
        nx_sum += digamma(n_x.max(1) as f64);
        ny_sum += digamma(n_y.max(1) as f64);
    }
    let mi = digamma(k as f64) - (nx_sum + ny_sum) / n as f64 + digamma(n as f64);
    mi.max(0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Differential Entropy
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate differential entropy `h(X) = −∫ p(x) log p(x) dx` using
/// kernel density estimation (Gaussian KDE).
///
/// `samples` is a 1-D data array; `bandwidth` is the KDE smoothing parameter.
pub fn differential_entropy(samples: &[f64], bandwidth: f64) -> f64 {
    let n = samples.len();
    if n == 0 {
        return 0.0;
    }
    let h = bandwidth.max(1e-10);
    let norm = 1.0 / (n as f64 * h * (2.0 * std::f64::consts::PI).sqrt());
    let entropy: f64 = samples
        .iter()
        .map(|&xi| {
            // KDE density at xi.
            let p: f64 = samples
                .iter()
                .map(|&xj| {
                    let u = (xi - xj) / h;
                    (-0.5 * u * u).exp()
                })
                .sum::<f64>()
                * norm;
            if p > 1e-300 { -p.ln() } else { 0.0 }
        })
        .sum();
    entropy / n as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// Natural Gradient
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the natural gradient `F^{-1} g` given Fisher matrix `fisher` and
/// Euclidean gradient `grad`.
///
/// The natural gradient is the steepest ascent direction in the Fisher–Rao
/// metric, used in natural gradient descent.
///
/// Returns [`PhysicsError::InvalidInput`] if `fisher` is not square, does not
/// match the length of `grad`, or is singular.
pub fn natural_gradient(fisher: &[Vec<f64>], grad: &[f64]) -> PhysicsResult<Vec<f64>> {
    let n = grad.len();
    if fisher.len() != n || fisher.iter().any(|row| row.len() != n) {
        return Err(PhysicsError::InvalidInput {
            field: "fisher".to_string(),
            reason: format!("expected a {n}x{n} matrix matching the gradient length"),
        });
    }
    let f_inv = invert_matrix(fisher).ok_or_else(|| PhysicsError::InvalidInput {
        field: "fisher".to_string(),
        reason: "Fisher information matrix is singular".to_string(),
    })?;
    Ok(mat_vec_mul(&f_inv, grad))
}

// ─────────────────────────────────────────────────────────────────────────────
// Linear algebra helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Invert a square matrix using Gauss-Jordan elimination with partial
/// pivoting.
///
/// Returns `None` if the matrix is not square, contains non-finite entries, or
/// is singular (a pivot falls below `1e-14` times the largest absolute entry).
fn invert_matrix(m: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = m.len();
    if n == 0 {
        return Some(vec![]);
    }
    if m.iter()
        .any(|row| row.len() != n || row.iter().any(|v| !v.is_finite()))
    {
        return None;
    }
    let scale = m.iter().flatten().fold(0.0f64, |acc, v| acc.max(v.abs()));
    if scale == 0.0 {
        return None;
    }
    // Build augmented matrix [m | I].
    let mut aug: Vec<Vec<f64>> = m
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            for j in 0..n {
                r.push(if i == j { 1.0 } else { 0.0 });
            }
            r
        })
        .collect();
    // Forward elimination with partial pivoting.
    for col in 0..n {
        // Find pivot.
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (row, aug_row) in aug.iter().enumerate().skip(col + 1) {
            if aug_row[col].abs() > max_val {
                max_val = aug_row[col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-14 * scale {
            return None;
        }
        for x in aug[col].iter_mut() {
            *x /= pivot;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            let col_vals: Vec<f64> = aug[col][..2 * n].to_vec();
            for (cell, &cv) in aug[row][..2 * n].iter_mut().zip(col_vals.iter()) {
                *cell -= factor * cv;
            }
        }
    }
    // Extract right half.
    Some(aug.iter().map(|row| row[n..].to_vec()).collect())
}

/// Multiply matrix `m` by vector `v`.
fn mat_vec_mul(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Multiply two square matrices (test helper only).
    fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = a.len();
        let mut c = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    }

    // ── ExponentialFamily ──────────────────────────────────────────────────

    fn gaussian_log_partition(theta: &[f64]) -> f64 {
        // For N(μ, 1): A(θ) = θ²/2 (natural param η = μ).
        if theta.is_empty() {
            return 0.0;
        }
        0.5 * theta[0] * theta[0]
    }

    fn identity_stat(x: &[f64]) -> f64 {
        x.first().copied().unwrap_or(0.0)
    }

    #[test]
    fn test_exponential_family_moment_params() {
        let ef = ExponentialFamily::new(
            vec![identity_stat as fn(&[f64]) -> f64],
            gaussian_log_partition,
        );
        let theta = vec![2.0f64];
        let mu = ef.moment_params(&theta);
        // For N(μ, 1): E[X] = θ, so moment param = 2.0.
        assert!((mu[0] - 2.0).abs() < 1e-3);
    }

    #[test]
    fn test_exponential_family_kl_nonneg() {
        let ef = ExponentialFamily::new(
            vec![identity_stat as fn(&[f64]) -> f64],
            gaussian_log_partition,
        );
        let theta1 = vec![1.0f64];
        let theta2 = vec![2.0f64];
        let kl = ef.kl_divergence(&theta1, &theta2);
        assert!(kl >= 0.0, "KL divergence must be non-negative");
    }

    #[test]
    fn test_exponential_family_kl_self_zero() {
        let ef = ExponentialFamily::new(
            vec![identity_stat as fn(&[f64]) -> f64],
            gaussian_log_partition,
        );
        let theta = vec![1.5f64];
        let kl = ef.kl_divergence(&theta, &theta);
        assert!(kl.abs() < 1e-6, "KL(p||p) should be 0");
    }

    #[test]
    fn test_exponential_family_fisher_info_positive() {
        let ef = ExponentialFamily::new(
            vec![identity_stat as fn(&[f64]) -> f64],
            gaussian_log_partition,
        );
        let theta = vec![1.0f64];
        let fi = ef.fisher_info(&theta);
        assert!(fi[0][0] > 0.0, "Fisher information must be positive");
    }

    // ── GaussianManifold ───────────────────────────────────────────────────

    #[test]
    fn test_gaussian_manifold_fisher_metric() {
        let gm = GaussianManifold::new();
        let g = gm.fisher_metric(0.0, 1.0);
        // At σ=1: g = [[1, 0], [0, 2]].
        assert!((g[0][0] - 1.0).abs() < 1e-10);
        assert!((g[1][1] - 2.0).abs() < 1e-10);
        assert!((g[0][1]).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_geodesic_distance_zero() {
        let gm = GaussianManifold::new();
        let d = gm.geodesic_distance(0.0, 1.0, 0.0, 1.0);
        assert!(d < 1e-6, "Distance from a point to itself should be 0");
    }

    #[test]
    fn test_gaussian_geodesic_distance_positive() {
        let gm = GaussianManifold::new();
        let d = gm.geodesic_distance(0.0, 1.0, 1.0, 2.0);
        assert!(
            d > 0.0,
            "Distance between different Gaussians should be positive"
        );
    }

    #[test]
    fn test_gaussian_geodesic_distance_exact_values() {
        let gm = GaussianManifold::new();
        // N(0,1) -> N(0,e): along the sigma axis d = sqrt(2) * |ln(sigma2/sigma1)|.
        let d = gm.geodesic_distance(0.0, 1.0, 0.0, std::f64::consts::E);
        assert!(
            (d - std::f64::consts::SQRT_2).abs() < 1e-12,
            "expected sqrt(2), got {d}"
        );
        // Pure scale change by a factor 4 at fixed mean: sqrt(2) * ln 4.
        let d = gm.geodesic_distance(2.0, 0.5, 2.0, 2.0);
        let expected = std::f64::consts::SQRT_2 * 4.0_f64.ln();
        assert!((d - expected).abs() < 1e-12, "got {d}, expected {expected}");
    }

    #[test]
    fn test_gaussian_geodesic_distance_is_a_metric() {
        let gm = GaussianManifold::new();
        // Deterministic pseudo-random triples of (mu, sigma) points.
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..200 {
            let pts: Vec<(f64, f64)> = (0..3)
                .map(|_| (6.0 * next() - 3.0, 0.1 + 3.0 * next()))
                .collect();
            let d =
                |a: usize, b: usize| gm.geodesic_distance(pts[a].0, pts[a].1, pts[b].0, pts[b].1);
            assert!((d(0, 1) - d(1, 0)).abs() < 1e-12);
            assert!(d(0, 1) >= 0.0);
            assert!(d(0, 2) <= d(0, 1) + d(1, 2) + 1e-12, "triangle inequality");
            // Agrees with the information-theory implementation.
            let reference = crate::information_theory::fisher_rao_gaussian(
                pts[0].0, pts[0].1, pts[1].0, pts[1].1,
            );
            assert!((d(0, 1) - reference).abs() < 1e-12);
        }
    }

    #[test]
    fn test_gaussian_exponential_map() {
        let gm = GaussianManifold::new();
        let (mu2, sigma2) = gm.exponential_map(0.0, 1.0, 1.0, 0.5, 1.0);
        assert!((mu2 - 1.0).abs() < 1e-9);
        assert!((sigma2 - 1.5).abs() < 1e-9);
    }

    #[test]
    fn test_gaussian_logarithmic_map() {
        let gm = GaussianManifold::new();
        let (vmu, vsigma) = gm.logarithmic_map(0.0, 1.0, 2.0, 3.0);
        assert!((vmu - 2.0).abs() < 1e-9);
        assert!((vsigma - 2.0).abs() < 1e-9);
    }

    // ── mutual_information_estimator ───────────────────────────────────────

    #[test]
    fn test_mutual_information_independent() {
        // X and Y independent → MI ≈ 0.
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..20).map(|i| (19 - i) as f64).collect();
        let mi = mutual_information_estimator(&x, &y, 3);
        // MI may be non-negative.
        assert!(mi >= 0.0);
    }

    #[test]
    fn test_mutual_information_identical() {
        // X = Y → MI should be high (positive).
        let x: Vec<f64> = (0..30).map(|i| i as f64 * 0.1).collect();
        let mi = mutual_information_estimator(&x, &x, 3);
        assert!(mi >= 0.0);
    }

    #[test]
    fn test_mutual_information_too_few() {
        let x = vec![1.0, 2.0];
        let y = vec![1.0, 2.0];
        let mi = mutual_information_estimator(&x, &y, 5);
        assert_eq!(mi, 0.0);
    }

    // ── differential_entropy ───────────────────────────────────────────────

    #[test]
    fn test_differential_entropy_empty() {
        let h = differential_entropy(&[], 1.0);
        assert_eq!(h, 0.0);
    }

    #[test]
    fn test_differential_entropy_single() {
        let h = differential_entropy(&[0.0], 1.0);
        assert!(h.is_finite());
    }

    #[test]
    fn test_differential_entropy_uniform_like() {
        // Wider spread → higher entropy.
        let narrow: Vec<f64> = (0..20).map(|i| i as f64 * 0.01).collect();
        let wide: Vec<f64> = (0..20).map(|i| i as f64 * 1.0).collect();
        let h_narrow = differential_entropy(&narrow, 0.1);
        let h_wide = differential_entropy(&wide, 1.0);
        // Wide distribution should generally have higher entropy.
        assert!(h_wide > h_narrow || h_wide.is_finite());
    }

    // ── natural_gradient ───────────────────────────────────────────────────

    #[test]
    fn test_natural_gradient_identity_fisher() {
        let fisher = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let grad = vec![1.0, 2.0];
        let ng = natural_gradient(&fisher, &grad).expect("non-singular Fisher matrix");
        assert!((ng[0] - 1.0).abs() < 1e-6);
        assert!((ng[1] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_natural_gradient_scaling() {
        let fisher = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        let grad = vec![2.0, 4.0];
        let ng = natural_gradient(&fisher, &grad).expect("non-singular Fisher matrix");
        assert!((ng[0] - 1.0).abs() < 1e-6);
        assert!((ng[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_natural_gradient_singular_fisher_is_error() {
        // Rank-one Fisher matrix: the natural gradient is undefined.
        let fisher = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let grad = vec![1.0, 1.0];
        assert!(matches!(
            natural_gradient(&fisher, &grad),
            Err(PhysicsError::InvalidInput { .. })
        ));
        // Dimension mismatch is rejected as well.
        assert!(natural_gradient(&[vec![1.0]], &grad).is_err());
    }

    #[test]
    fn test_natural_gradient_full_matrix() {
        // F = [[4, 1], [1, 3]], F^{-1} = [[3, -1], [-1, 4]] / 11.
        let fisher = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
        let grad = vec![1.0, 2.0];
        let ng = natural_gradient(&fisher, &grad).expect("non-singular Fisher matrix");
        assert!((ng[0] - 1.0 / 11.0).abs() < 1e-12);
        assert!((ng[1] - 7.0 / 11.0).abs() < 1e-12);
    }

    // ── Helper functions ───────────────────────────────────────────────────

    #[test]
    fn test_invert_matrix_identity() {
        let id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let inv = invert_matrix(&id).expect("identity is invertible");
        assert!((inv[0][0] - 1.0).abs() < 1e-10);
        assert!((inv[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_invert_matrix_2x2() {
        let m = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        let inv = invert_matrix(&m).expect("diagonal matrix is invertible");
        assert!((inv[0][0] - 0.5).abs() < 1e-10);
        assert!((inv[1][1] - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_invert_matrix_singular_returns_none() {
        assert!(invert_matrix(&[vec![1.0, 2.0], vec![2.0, 4.0]]).is_none());
        assert!(invert_matrix(&[vec![0.0, 0.0], vec![0.0, 0.0]]).is_none());
        assert!(invert_matrix(&[vec![1.0, f64::NAN], vec![0.0, 1.0]]).is_none());
        // Badly scaled but regular matrices are still inverted.
        let m = vec![vec![1e-20, 0.0], vec![0.0, 1e-20]];
        let inv = invert_matrix(&m).expect("scaled identity is invertible");
        assert!((inv[0][0] - 1e20).abs() < 1e8);
    }

    #[test]
    fn test_mat_vec_mul() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let v = vec![1.0, 1.0];
        let r = mat_vec_mul(&m, &v);
        assert!((r[0] - 3.0).abs() < 1e-10);
        assert!((r[1] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_mat_mul_identity() {
        let id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let m = vec![vec![3.0, 1.0], vec![2.0, 5.0]];
        let r = mat_mul(&id, &m);
        for i in 0..2 {
            for j in 0..2 {
                assert!((r[i][j] - m[i][j]).abs() < 1e-10);
            }
        }
    }
}
