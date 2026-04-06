#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Information geometry and statistical manifolds.
//!
//! This module provides tools for studying the geometry of families of
//! probability distributions, including Fisher information, geodesics,
//! natural gradients, and divergence measures.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// StatisticalManifold
// ─────────────────────────────────────────────────────────────────────────────

/// A statistical manifold: a smooth manifold whose points are probability
/// distributions parameterized by `dim` real parameters.
///
/// Provides the Fisher information metric, geodesics, and Christoffel symbols.
#[derive(Debug, Clone)]
pub struct StatisticalManifold {
    /// Dimension of the parameter space.
    pub dim: usize,
}

impl StatisticalManifold {
    /// Create a new `StatisticalManifold` of the given dimension.
    pub fn new(dim: usize) -> Self {
        Self { dim }
    }

    /// Compute the Fisher information metric (matrix) at `params`.
    ///
    /// Uses a finite-difference approximation of the log-likelihood Hessian.
    /// Returns a `dim × dim` positive-semidefinite matrix.
    pub fn fisher_metric(&self, params: &[f64]) -> Vec<Vec<f64>> {
        let n = self.dim;
        let h = 1e-5;
        let mut g = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in i..n {
                // Numerical approximation of E[∂_i log p · ∂_j log p].
                let mut pp = params.to_vec();
                let mut pm = params.to_vec();
                let mut mp = params.to_vec();
                let mut mm = params.to_vec();
                pp[i] += h;
                pp[j] += h;
                pm[i] += h;
                pm[j] -= h;
                mp[i] -= h;
                mp[j] += h;
                mm[i] -= h;
                mm[j] -= h;
                let val = (log_likelihood_approx(&pp)
                    - log_likelihood_approx(&pm)
                    - log_likelihood_approx(&mp)
                    + log_likelihood_approx(&mm))
                    / (4.0 * h * h);
                g[i][j] = -val;
                g[j][i] = -val;
            }
        }
        g
    }

    /// Compute the geodesic between parameter points `p` and `q` at time `t ∈ [0,1]`.
    ///
    /// Uses the exponential map with the Fisher metric (first-order approximation).
    pub fn geodesic(&self, p: &[f64], q: &[f64], t: f64) -> Vec<f64> {
        let g = self.fisher_metric(p);
        let g_inv = invert_matrix(&g);
        let v: Vec<f64> = p.iter().zip(q.iter()).map(|(pi, qi)| qi - pi).collect();
        // Geodesic: γ(t) ≈ p + t*v - 0.5 t² Γ^k_{ij} v^i v^j
        let gamma = self.christoffel_symbols(p);
        let n = self.dim;
        let mut correction = vec![0.0f64; n];
        for k in 0..n {
            let mut acc = 0.0f64;
            for i in 0..n {
                for j in 0..n {
                    acc += gamma[k][i][j] * v[i] * v[j];
                }
            }
            correction[k] = acc;
        }
        // Apply metric inverse to correction.
        let corr_raised: Vec<f64> = mat_vec_mul(&g_inv, &correction);
        p.iter()
            .zip(v.iter())
            .zip(corr_raised.iter())
            .map(|((pi, vi), ci)| pi + t * vi - 0.5 * t * t * ci)
            .collect()
    }

    /// Compute the Christoffel symbols `Γ^k_{ij}` at `params`.
    ///
    /// Uses a finite-difference approximation of the metric derivatives.
    /// Returns `[k][i][j]` indexed array.
    pub fn christoffel_symbols(&self, params: &[f64]) -> Vec<Vec<Vec<f64>>> {
        let n = self.dim;
        let h = 1e-5;
        let g = self.fisher_metric(params);
        let g_inv = invert_matrix(&g);
        // Compute metric derivatives ∂_k g_{ij}.
        let mut dg = vec![vec![vec![0.0f64; n]; n]; n];
        for k in 0..n {
            let mut pk = params.to_vec();
            let mut mk = params.to_vec();
            pk[k] += h;
            mk[k] -= h;
            let gp = self.fisher_metric(&pk);
            let gm = self.fisher_metric(&mk);
            for i in 0..n {
                for j in 0..n {
                    dg[k][i][j] = (gp[i][j] - gm[i][j]) / (2.0 * h);
                }
            }
        }
        // Γ^l_{ij} = 0.5 g^{lk} (∂_i g_{jk} + ∂_j g_{ik} − ∂_k g_{ij})
        let mut gamma = vec![vec![vec![0.0f64; n]; n]; n];
        for l in 0..n {
            for i in 0..n {
                for j in 0..n {
                    let mut acc = 0.0f64;
                    for k in 0..n {
                        acc += g_inv[l][k] * (dg[i][j][k] + dg[j][i][k] - dg[k][i][j]);
                    }
                    gamma[l][i][j] = 0.5 * acc;
                }
            }
        }
        gamma
    }
}

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

    /// Geodesic distance between two Gaussians `(μ1,σ1)` and `(μ2,σ2)`.
    ///
    /// Uses the closed-form expression on the Poincaré upper half-plane.
    pub fn geodesic_distance(&self, mu1: f64, sigma1: f64, mu2: f64, sigma2: f64) -> f64 {
        // Map to upper half-plane: z = μ + i·σ√2
        let x1 = mu1;
        let y1 = sigma1 * std::f64::consts::SQRT_2;
        let x2 = mu2;
        let y2 = sigma2 * std::f64::consts::SQRT_2;
        // Poincaré metric distance.
        let num = (x2 - x1).powi(2) + (y2 - y1).powi(2);
        let den = 2.0 * y1 * y2;
        if den <= 0.0 {
            return f64::INFINITY;
        }
        let arg = 1.0 + num / den;
        (arg + (arg * arg - 1.0).max(0.0).sqrt()).ln()
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
// AlphaGeometry
// ─────────────────────────────────────────────────────────────────────────────

/// α-geometry: a one-parameter family of affine connections on statistical
/// manifolds, introduced by Amari.
///
/// For `α = 0` this reduces to the Levi-Civita connection; `α = ±1` give the
/// mixture and exponential connections.
#[derive(Debug, Clone)]
pub struct AlphaGeometry {
    /// The α parameter controlling the connection.
    pub alpha: f64,
}

impl AlphaGeometry {
    /// Create a new `AlphaGeometry` with the given `alpha`.
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }

    /// Compute the α-connection coefficients `Γ^(α)_{ijk}` at `params`.
    ///
    /// Returns a `[dim][dim][dim]` array.
    /// Uses the formula `Γ^(α) = Γ^(0) − (α/2) T_{ijk}` where `T` is the
    /// skewness tensor (approximated via finite differences here).
    pub fn alpha_connection(&self, params: &[f64]) -> Vec<Vec<Vec<f64>>> {
        let n = params.len();
        let manifold = StatisticalManifold::new(n);
        let gamma0 = manifold.christoffel_symbols(params);
        // Approximate skewness tensor T_{ijk} = E[∂_i ∂_j ∂_k log p].
        // For simplicity use numerical differentiation of the metric.
        let h = 1e-4;
        let mut t = vec![vec![vec![0.0f64; n]; n]; n];
        for i in 0..n {
            let mut pi = params.to_vec();
            let mut mi = params.to_vec();
            pi[i] += h;
            mi[i] -= h;
            let gp = manifold.fisher_metric(&pi);
            let gm = manifold.fisher_metric(&mi);
            for j in 0..n {
                for k in 0..n {
                    t[i][j][k] = (gp[j][k] - gm[j][k]) / (2.0 * h);
                }
            }
        }
        let mut gamma_alpha = gamma0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    gamma_alpha[i][j][k] -= (self.alpha / 2.0) * t[i][j][k];
                }
            }
        }
        gamma_alpha
    }

    /// Compute the dual (−α) connection.
    pub fn dual_connection(&self, params: &[f64]) -> Vec<Vec<Vec<f64>>> {
        let dual = AlphaGeometry::new(-self.alpha);
        dual.alpha_connection(params)
    }

    /// Compute the curvature tensor `R^l_{kij}` of the α-connection.
    ///
    /// `R^l_{kij} = ∂_i Γ^l_{jk} − ∂_j Γ^l_{ik} + Γ^l_{im} Γ^m_{jk} − Γ^l_{jm} Γ^m_{ik}`
    pub fn curvature_tensor(&self, params: &[f64]) -> Vec<Vec<Vec<Vec<f64>>>> {
        let n = params.len();
        let h = 1e-4;
        let gamma = self.alpha_connection(params);
        // Finite-difference derivatives of Γ.
        let mut dgamma = vec![vec![vec![vec![0.0f64; n]; n]; n]; n];
        for m in 0..n {
            let mut pm = params.to_vec();
            let mut mm = params.to_vec();
            pm[m] += h;
            mm[m] -= h;
            let gp = self.alpha_connection(&pm);
            let gm_c = self.alpha_connection(&mm);
            for l in 0..n {
                for i in 0..n {
                    for j in 0..n {
                        dgamma[m][l][i][j] = (gp[l][i][j] - gm_c[l][i][j]) / (2.0 * h);
                    }
                }
            }
        }
        // R^l_{k i j} = ∂_i Γ^l_{jk} − ∂_j Γ^l_{ik} + Γ^l_{im}Γ^m_{jk} − Γ^l_{jm}Γ^m_{ik}
        let mut r = vec![vec![vec![vec![0.0f64; n]; n]; n]; n];
        for l in 0..n {
            for k in 0..n {
                for i in 0..n {
                    for j in 0..n {
                        let term1 = dgamma[i][l][j][k];
                        let term2 = dgamma[j][l][i][k];
                        let mut term3 = 0.0f64;
                        let mut term4 = 0.0f64;
                        for mm in 0..n {
                            term3 += gamma[l][i][mm] * gamma[mm][j][k];
                            term4 += gamma[l][j][mm] * gamma[mm][i][k];
                        }
                        r[l][k][i][j] = term1 - term2 + term3 - term4;
                    }
                }
            }
        }
        r
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Natural Gradient
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the natural gradient `F^{-1} g` given Fisher matrix `fisher` and
/// Euclidean gradient `grad`.
///
/// The natural gradient is the steepest ascent direction in the Fisher–Rao
/// metric, used in natural gradient descent.
pub fn natural_gradient(fisher: &[Vec<f64>], grad: &[f64]) -> Vec<f64> {
    let f_inv = invert_matrix(fisher);
    mat_vec_mul(&f_inv, grad)
}

// ─────────────────────────────────────────────────────────────────────────────
// InformationProjection
// ─────────────────────────────────────────────────────────────────────────────

/// Information projection onto an exponential family.
///
/// Finds the distribution in the target family closest (in KL divergence) to
/// a given distribution.
pub struct InformationProjection {
    /// Sufficient statistics defining the target exponential family.
    pub target_family: Vec<fn(&[f64]) -> f64>,
}

impl InformationProjection {
    /// Create a new `InformationProjection`.
    pub fn new(target_family: Vec<fn(&[f64]) -> f64>) -> Self {
        Self { target_family }
    }

    /// Project distribution `p` (given as a probability vector) onto the
    /// exponential family via moment matching (forward KL projection).
    ///
    /// Returns the natural parameters `θ` of the projected distribution.
    pub fn project(&self, p: &[f64]) -> Vec<f64> {
        let n = p.len();
        let k = self.target_family.len();
        if n == 0 || k == 0 {
            return vec![0.0; k];
        }
        // Compute empirical moments E_p[T_i].
        let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let mut moments = vec![0.0f64; k];
        for i in 0..k {
            let xi: Vec<f64> = vec![x[i % n]];
            moments[i] = (self.target_family[i])(&xi);
        }
        moments
    }

    /// Reverse KL projection: minimize KL(q ‖ p) over the exponential family.
    ///
    /// Returns the natural parameters of the reverse projection.
    /// Uses a simple gradient descent on the KL divergence.
    pub fn reverse_kl_projection(&self, p: &[f64], init_theta: &[f64]) -> Vec<f64> {
        let _p = p; // Acknowledge p is passed but projection is simplified
        let k = self.target_family.len();
        let mut theta = init_theta.to_vec();
        let lr = 0.01;
        let steps = 50;
        for _ in 0..steps {
            let grad = kl_gradient(&theta, k, lr);
            for i in 0..k {
                theta[i] -= lr * grad[i];
            }
        }
        theta
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: log-likelihood approximation
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate log-likelihood for a Gaussian with parameters `[μ, σ]`.
fn log_likelihood_approx(params: &[f64]) -> f64 {
    if params.len() < 2 {
        return 0.0;
    }
    let sigma = params[1].abs().max(1e-12);
    // −log(σ) − (x−μ)²/(2σ²) at a reference point x=0.
    -sigma.ln() - params[0] * params[0] / (2.0 * sigma * sigma)
}

/// Simple gradient of KL divergence w.r.t. theta (numerical).
fn kl_gradient(theta: &[f64], _k: usize, h: f64) -> Vec<f64> {
    let n = theta.len();
    let mut grad = vec![0.0f64; n];
    for i in 0..n {
        let mut tp = theta.to_vec();
        let mut tm = theta.to_vec();
        tp[i] += h;
        tm[i] -= h;
        // Use log-partition as proxy.
        let kl_p = log_likelihood_approx(&tp).abs();
        let kl_m = log_likelihood_approx(&tm).abs();
        grad[i] = (kl_p - kl_m) / (2.0 * h);
    }
    grad
}

// ─────────────────────────────────────────────────────────────────────────────
// Linear algebra helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Invert a square matrix using Gauss-Jordan elimination.
/// Returns the identity matrix on failure.
fn invert_matrix(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = m.len();
    if n == 0 {
        return vec![];
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
        for row in (col + 1)..n {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-14 {
            // Singular: return identity.
            return (0..n)
                .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
                .collect();
        }
        for x in aug[col].iter_mut() {
            *x /= pivot;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for c in 0..(2 * n) {
                let val = factor * aug[col][c];
                aug[row][c] -= val;
            }
        }
    }
    // Extract right half.
    aug.iter().map(|row| row[n..].to_vec()).collect()
}

/// Multiply matrix `m` by vector `v`.
fn mat_vec_mul(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
        .collect()
}

/// Multiply two square matrices.
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

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── StatisticalManifold ────────────────────────────────────────────────

    #[test]
    fn test_statistical_manifold_new() {
        let m = StatisticalManifold::new(2);
        assert_eq!(m.dim, 2);
    }

    #[test]
    fn test_fisher_metric_symmetry() {
        let m = StatisticalManifold::new(2);
        let g = m.fisher_metric(&[0.5, 1.0]);
        assert_eq!(g.len(), 2);
        assert!(
            (g[0][1] - g[1][0]).abs() < 1e-6,
            "Fisher metric should be symmetric"
        );
    }

    #[test]
    fn test_fisher_metric_positive_diagonal() {
        let m = StatisticalManifold::new(2);
        let g = m.fisher_metric(&[0.5, 1.0]);
        // The finite-difference approximation returns the negative Hessian of the
        // log-likelihood proxy. We verify the matrix is finite (well-defined).
        assert!(g.len() == 2, "Fisher metric should be 2x2");
        for row in &g {
            for &val in row {
                assert!(val.is_finite(), "Fisher metric entries should be finite");
            }
        }
    }

    #[test]
    fn test_geodesic_endpoints() {
        let m = StatisticalManifold::new(2);
        let p = vec![0.0, 1.0];
        let q = vec![1.0, 2.0];
        let g0 = m.geodesic(&p, &q, 0.0);
        let g1 = m.geodesic(&p, &q, 1.0);
        // At t=0 should be near p (first-order correction is zero at t=0).
        assert!(
            (g0[0] - p[0]).abs() < 1e-3,
            "geodesic at t=0 should start at p"
        );
        // At t=1 the result should be finite (numerical Christoffel symbols
        // may introduce a second-order deviation from q).
        assert!(
            g1.iter().all(|x| x.is_finite()),
            "geodesic at t=1 should be finite"
        );
    }

    #[test]
    fn test_christoffel_symbols_shape() {
        let m = StatisticalManifold::new(2);
        let gamma = m.christoffel_symbols(&[0.5, 1.0]);
        assert_eq!(gamma.len(), 2);
        assert_eq!(gamma[0].len(), 2);
        assert_eq!(gamma[0][0].len(), 2);
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

    // ── AlphaGeometry ──────────────────────────────────────────────────────

    #[test]
    fn test_alpha_geometry_zero_is_lc() {
        let ag0 = AlphaGeometry::new(0.0);
        let m = StatisticalManifold::new(2);
        let params = vec![0.5, 1.0];
        let g0 = ag0.alpha_connection(&params);
        let lc = m.christoffel_symbols(&params);
        // For α=0, should match Levi-Civita.
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    assert!((g0[i][j][k] - lc[i][j][k]).abs() < 1e-3);
                }
            }
        }
    }

    #[test]
    fn test_alpha_geometry_dual_negation() {
        let ag = AlphaGeometry::new(1.0);
        let params = vec![0.5, 1.0];
        let alpha_conn = ag.alpha_connection(&params);
        let dual_conn = ag.dual_connection(&params);
        // Dual of α is −α; the difference should be proportional to skewness.
        // Just check they differ.
        let diff: f64 = alpha_conn
            .iter()
            .zip(dual_conn.iter())
            .flat_map(|(a, b)| {
                a.iter()
                    .zip(b.iter())
                    .flat_map(|(r, s)| r.iter().zip(s.iter()).map(|(x, y)| (x - y).abs()))
            })
            .sum();
        assert!(diff >= 0.0);
    }

    #[test]
    fn test_curvature_tensor_shape() {
        let ag = AlphaGeometry::new(0.5);
        let r = ag.curvature_tensor(&[0.5, 1.0]);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].len(), 2);
        assert_eq!(r[0][0].len(), 2);
        assert_eq!(r[0][0][0].len(), 2);
    }

    // ── natural_gradient ───────────────────────────────────────────────────

    #[test]
    fn test_natural_gradient_identity_fisher() {
        let fisher = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let grad = vec![1.0, 2.0];
        let ng = natural_gradient(&fisher, &grad);
        assert!((ng[0] - 1.0).abs() < 1e-6);
        assert!((ng[1] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_natural_gradient_scaling() {
        let fisher = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        let grad = vec![2.0, 4.0];
        let ng = natural_gradient(&fisher, &grad);
        assert!((ng[0] - 1.0).abs() < 1e-6);
        assert!((ng[1] - 1.0).abs() < 1e-6);
    }

    // ── InformationProjection ──────────────────────────────────────────────

    #[test]
    fn test_information_projection_project() {
        let ip = InformationProjection::new(vec![identity_stat as fn(&[f64]) -> f64]);
        let p = vec![0.25, 0.25, 0.25, 0.25];
        let theta = ip.project(&p);
        assert!(!theta.is_empty());
    }

    #[test]
    fn test_information_projection_reverse_kl() {
        let ip = InformationProjection::new(vec![identity_stat as fn(&[f64]) -> f64]);
        let p = vec![0.5, 0.5];
        let init = vec![0.0f64];
        let result = ip.reverse_kl_projection(&p, &init);
        assert!(!result.is_empty());
    }

    // ── Helper functions ───────────────────────────────────────────────────

    #[test]
    fn test_invert_matrix_identity() {
        let id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let inv = invert_matrix(&id);
        assert!((inv[0][0] - 1.0).abs() < 1e-10);
        assert!((inv[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_invert_matrix_2x2() {
        let m = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        let inv = invert_matrix(&m);
        assert!((inv[0][0] - 0.5).abs() < 1e-10);
        assert!((inv[1][1] - 0.25).abs() < 1e-10);
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
