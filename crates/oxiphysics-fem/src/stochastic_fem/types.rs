//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use rand::RngExt;
use std::f64::consts::PI;

/// Monte Carlo FEM solver for statistical post-processing.
///
/// Generates realisations by sampling random material parameters,
/// accumulates displacement responses, and computes statistics.
pub struct MonteCarloFemSolver {
    /// Mean stiffness value.
    pub k_mean: f64,
    /// Standard deviation of stiffness.
    pub k_std: f64,
    /// Applied load.
    pub load: f64,
    /// Collected displacement samples.
    pub samples: Vec<f64>,
}
impl MonteCarloFemSolver {
    /// Create a new [`MonteCarloFemSolver`].
    ///
    /// # Arguments
    /// * `k_mean` – mean stiffness
    /// * `k_std`  – standard deviation of stiffness
    /// * `load`   – applied force
    pub fn new(k_mean: f64, k_std: f64, load: f64) -> Self {
        Self {
            k_mean,
            k_std,
            load,
            samples: Vec::new(),
        }
    }
    /// Run `n_samples` Monte Carlo realisations and store the results.
    ///
    /// Each realisation draws k ~ N(k_mean, k_std²) (clipped to k_mean·1e-6),
    /// then computes the response u = load / k.
    pub fn run(&mut self, n_samples: usize) {
        use rand::RngExt as _;
        let mut rng = rand::rng();
        self.samples.reserve(n_samples);
        for _ in 0..n_samples {
            let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
            let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
            let z = box_muller(u1, u2);
            let k = (self.k_mean + self.k_std * z).max(self.k_mean * 1e-6);
            self.samples.push(self.load / k);
        }
    }
    /// Clear all accumulated samples.
    pub fn clear(&mut self) {
        self.samples.clear();
    }
    /// Sample mean.
    pub fn mean(&self) -> f64 {
        if self.samples.is_empty() {
            0.0
        } else {
            self.samples.iter().sum::<f64>() / self.samples.len() as f64
        }
    }
    /// Sample standard deviation.
    pub fn std_dev(&self) -> f64 {
        let n = self.samples.len();
        if n < 2 {
            return 0.0;
        }
        let mu = self.mean();
        let var = self.samples.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / n as f64;
        var.sqrt()
    }
    /// Coefficient of variation (σ / |μ|).
    pub fn cov(&self) -> f64 {
        let mu = self.mean().abs();
        if mu < f64::EPSILON {
            0.0
        } else {
            self.std_dev() / mu
        }
    }
    /// Empirical p-th percentile (0 ≤ p ≤ 100).
    pub fn percentile(&self, p: f64) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        let mut sorted = self.samples.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idx = ((p / 100.0) * (n - 1) as f64).round() as usize;
        sorted[idx.min(n - 1)]
    }
    /// Estimate failure probability P(u > threshold).
    pub fn failure_probability(&self, threshold: f64) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        let n_fail = self.samples.iter().filter(|&&s| s > threshold).count();
        n_fail as f64 / self.samples.len() as f64
    }
    /// 95% confidence interval on the mean: (lower, upper).
    pub fn confidence_interval_mean_95(&self) -> (f64, f64) {
        let n = self.samples.len();
        if n < 2 {
            return (self.mean(), self.mean());
        }
        let mu = self.mean();
        let se = self.std_dev() / (n as f64).sqrt();
        (mu - 1.96 * se, mu + 1.96 * se)
    }
}
/// Karhunen-Loève (K-L) expansion for a second-order random field.
///
/// Represents a random field as a truncated series:
///
/// ```text
/// X(x, θ) ≈ μ(x) + Σᵢ √λᵢ ξᵢ φᵢ(x)
/// ```
///
/// where `λᵢ` are eigenvalues and `φᵢ` are eigenvectors of the covariance
/// operator, and `ξᵢ` are uncorrelated standard-normal random variables.
pub struct KarhunenLoeveExpansion {
    /// Number of retained modes.
    pub n_terms: usize,
    /// Eigenvalues λᵢ (variance contributions per mode).
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors φᵢ, each of length equal to the spatial discretisation.
    pub eigenvectors: Vec<Vec<f64>>,
}
impl KarhunenLoeveExpansion {
    /// Create a new empty K-L expansion retaining `n_terms` modes.
    pub fn new(n_terms: usize) -> Self {
        Self {
            n_terms,
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
        }
    }
    /// Add a mode with given `eigenvalue` and `eigenvector`.
    ///
    /// # Arguments
    /// * `eigenvalue`  – λᵢ for this mode
    /// * `eigenvector` – φᵢ (nodal values)
    pub fn add_mode(&mut self, eigenvalue: f64, eigenvector: Vec<f64>) {
        self.eigenvalues.push(eigenvalue);
        self.eigenvectors.push(eigenvector);
    }
    /// Evaluate a realisation of the random field given standard-normal
    /// random variables `xi`.
    ///
    /// Returns the nodal field values: f(x) = Σᵢ √λᵢ ξᵢ φᵢ(x).
    ///
    /// # Arguments
    /// * `xi` – slice of standard-normal samples, one per mode
    pub fn sample(&self, xi: &[f64]) -> Vec<f64> {
        let n_modes = self.eigenvalues.len().min(xi.len());
        if n_modes == 0 {
            return Vec::new();
        }
        let n_nodes = self.eigenvectors[0].len();
        let mut result = vec![0.0f64; n_nodes];
        for ((&ev, &xi_i), eigvec) in self
            .eigenvalues
            .iter()
            .zip(xi.iter())
            .zip(self.eigenvectors.iter())
            .take(n_modes)
        {
            let scale = ev.max(0.0).sqrt() * xi_i;
            for (res_j, &e_j) in result.iter_mut().zip(eigvec.iter()) {
                *res_j += scale * e_j;
            }
        }
        result
    }
    /// Total variance captured by all retained modes (sum of eigenvalues).
    pub fn total_variance(&self) -> f64 {
        self.eigenvalues.iter().cloned().sum()
    }
    /// Energy ratio: fraction of total variance captured by mode `i`.
    ///
    /// Returns 0 if no modes have been added or total variance is zero.
    pub fn energy_ratio(&self, mode_idx: usize) -> f64 {
        let total = self.total_variance();
        if total < f64::EPSILON || mode_idx >= self.eigenvalues.len() {
            return 0.0;
        }
        self.eigenvalues[mode_idx] / total
    }
    /// Cumulative energy ratio for the first `k` modes.
    ///
    /// Returns the fraction of total variance captured by modes 0..k.
    pub fn cumulative_energy(&self, k: usize) -> f64 {
        let total = self.total_variance();
        if total < f64::EPSILON {
            return 0.0;
        }
        self.eigenvalues.iter().take(k).sum::<f64>() / total
    }
    /// Build a KL expansion from a symmetric positive semi-definite covariance
    /// matrix via the power iteration method.
    ///
    /// This extracts `n_modes` dominant eigenpairs using deflated power iteration.
    /// The matrix is given as a flat row-major array of size n×n.
    ///
    /// # Arguments
    /// * `cov`     – n×n covariance matrix (flat row-major)
    /// * `n`       – dimension (number of nodes)
    /// * `n_modes` – number of modes to retain
    /// * `max_iter`– maximum power-iteration steps per mode
    pub fn from_covariance(cov: &[f64], n: usize, n_modes: usize, max_iter: usize) -> Self {
        let mut kle = KarhunenLoeveExpansion::new(n_modes);
        if n == 0 || n_modes == 0 {
            return kle;
        }
        let mut residual: Vec<f64> = cov.to_vec();
        for _mode in 0..n_modes {
            let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
            let mut lambda = 0.0_f64;
            for _ in 0..max_iter {
                let mut w = vec![0.0f64; n];
                for i in 0..n {
                    for j in 0..n {
                        w[i] += residual[i * n + j] * v[j];
                    }
                }
                let rq: f64 = v.iter().zip(w.iter()).map(|(a, b)| a * b).sum();
                lambda = rq;
                let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm < 1e-30 {
                    break;
                }
                v = w.iter().map(|x| x / norm).collect();
            }
            if lambda > 0.0 {
                for i in 0..n {
                    for j in 0..n {
                        residual[i * n + j] -= lambda * v[i] * v[j];
                    }
                }
                kle.add_mode(lambda, v);
            }
        }
        kle
    }
}
/// Stochastic FEM descriptor with random material parameters.
///
/// Stores mean and standard deviation of the stiffness coefficient for use
/// in Monte Carlo and polynomial chaos analyses.
pub struct StochasticFEM {
    /// Number of spatial nodes in the FEM model.
    pub n_nodes: usize,
    /// Number of independent random variables.
    pub n_random: usize,
    /// Mean stiffness parameter k̄.
    pub mean_k: f64,
    /// Standard deviation of the stiffness σ_k.
    pub std_k: f64,
}
impl StochasticFEM {
    /// Create a new [`StochasticFEM`] descriptor.
    ///
    /// # Arguments
    /// * `n_nodes`  – number of nodes in the FEM mesh
    /// * `n_random` – number of independent random variables
    /// * `mean_k`   – mean stiffness
    /// * `std_k`    – standard deviation of stiffness
    pub fn new(n_nodes: usize, n_random: usize, mean_k: f64, std_k: f64) -> Self {
        Self {
            n_nodes,
            n_random,
            mean_k,
            std_k,
        }
    }
    /// Coefficient of variation (σ/μ).
    pub fn cov(&self) -> f64 {
        if self.mean_k.abs() < f64::EPSILON {
            return 0.0;
        }
        self.std_k / self.mean_k
    }
}
/// Full stochastic FEM problem with random field discretisation.
///
/// Manages the mean and variance fields over a 1-D mesh of nodes, and
/// provides utilities for constructing covariance matrices and extracting
/// statistical moments from ensemble results.
pub struct StochasticFemProblem {
    /// Node coordinates (1-D).
    pub nodes: Vec<f64>,
    /// Mean value of the field at each node.
    pub mean_field: Vec<f64>,
    /// Variance of the field at each node.
    pub variance_field: Vec<f64>,
    /// Correlation length used when building the covariance matrix.
    pub corr_length: f64,
    /// Standard deviation of the underlying random field.
    pub sigma: f64,
}
impl StochasticFemProblem {
    /// Create a [`StochasticFemProblem`] for a uniform mesh of `n_nodes` nodes
    /// spanning `[0, length]`.
    ///
    /// # Arguments
    /// * `n_nodes`     – number of nodes
    /// * `length`      – domain length
    /// * `mean_value`  – constant mean field value
    /// * `sigma`       – standard deviation of the random field
    /// * `corr_length` – correlation length for the covariance kernel
    pub fn new_uniform(
        n_nodes: usize,
        length: f64,
        mean_value: f64,
        sigma: f64,
        corr_length: f64,
    ) -> Self {
        let nodes: Vec<f64> = (0..n_nodes)
            .map(|i| i as f64 * length / (n_nodes.max(2) - 1) as f64)
            .collect();
        let mean_field = vec![mean_value; n_nodes];
        let variance_field = vec![sigma * sigma; n_nodes];
        Self {
            nodes,
            mean_field,
            variance_field,
            corr_length,
            sigma,
        }
    }
    /// Return the number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// Build the dense covariance matrix C\[i,j\] using the exponential kernel.
    ///
    /// C\[i,j\] = σ² · exp(−|x_i − x_j| / l_corr)
    pub fn covariance_matrix_exponential(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let mut c = vec![vec![0.0f64; n]; n];
        for (i, c_row) in c.iter_mut().enumerate() {
            for (j, c_ij) in c_row.iter_mut().enumerate() {
                *c_ij = covariance_exponential(
                    self.nodes[i],
                    self.nodes[j],
                    self.sigma,
                    self.corr_length,
                );
            }
        }
        c
    }
    /// Build the dense covariance matrix C\[i,j\] using the Gaussian kernel.
    ///
    /// C\[i,j\] = σ² · exp(−(x_i − x_j)² / (2 l_corr²))
    pub fn covariance_matrix_gaussian(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let mut c = vec![vec![0.0f64; n]; n];
        for (i, c_row) in c.iter_mut().enumerate() {
            for (j, c_ij) in c_row.iter_mut().enumerate() {
                *c_ij =
                    covariance_gaussian(self.nodes[i], self.nodes[j], self.sigma, self.corr_length);
            }
        }
        c
    }
    /// Update mean and variance fields from an ensemble of realisations.
    ///
    /// Each row of `ensemble` is one realisation of the nodal field.
    /// Overwrites `mean_field` and `variance_field` in place.
    pub fn update_statistics_from_ensemble(&mut self, ensemble: &[Vec<f64>]) {
        let n = self.nodes.len();
        if ensemble.is_empty() {
            return;
        }
        let m = ensemble.len() as f64;
        let mut mean = vec![0.0f64; n];
        for sample in ensemble.iter() {
            for (i, &v) in sample.iter().enumerate().take(n) {
                mean[i] += v;
            }
        }
        for v in mean.iter_mut() {
            *v /= m;
        }
        let mut var = vec![0.0f64; n];
        for sample in ensemble.iter() {
            for (i, &v) in sample.iter().enumerate().take(n) {
                var[i] += (v - mean[i]).powi(2);
            }
        }
        for v in var.iter_mut() {
            *v /= m;
        }
        self.mean_field = mean;
        self.variance_field = var;
    }
    /// Coefficient of variation at each node: CoV\[i\] = √(Var\[i\]) / |mean\[i\]|.
    ///
    /// Returns 0 wherever the mean is effectively zero.
    pub fn cov_field(&self) -> Vec<f64> {
        self.mean_field
            .iter()
            .zip(self.variance_field.iter())
            .map(|(&mu, &var)| {
                if mu.abs() < f64::EPSILON {
                    0.0
                } else {
                    var.sqrt() / mu.abs()
                }
            })
            .collect()
    }
}
/// Structural reliability analysis coupled to FEM.
///
/// Implements FORM and SORM approximations, probability of failure, and
/// reliability index β for a limit state function g(X) = capacity − demand.
pub struct ReliabilityFem {
    /// Mean values of the random variables X.
    pub means: Vec<f64>,
    /// Standard deviations of the random variables X.
    pub std_devs: Vec<f64>,
    /// Limit state function coefficients (linear model): g(X) = c₀ + Σ cᵢ Xᵢ.
    pub lsf_coeffs: Vec<f64>,
}
impl ReliabilityFem {
    /// Create a new [`ReliabilityFem`] with given statistical parameters.
    ///
    /// # Arguments
    /// * `means`      – mean values of random variables
    /// * `std_devs`   – standard deviations of random variables
    /// * `lsf_coeffs` – linear limit state function coefficients c₀, c₁, …
    pub fn new(means: Vec<f64>, std_devs: Vec<f64>, lsf_coeffs: Vec<f64>) -> Self {
        Self {
            means,
            std_devs,
            lsf_coeffs,
        }
    }
    /// Evaluate the limit state function at a point in original space.
    ///
    /// g(x) = c₀ + c₁x₁ + c₂x₂ + …
    pub fn lsf(&self, x: &[f64]) -> f64 {
        let k = self.lsf_coeffs.len().saturating_sub(1);
        let mut g = if self.lsf_coeffs.is_empty() {
            0.0
        } else {
            self.lsf_coeffs[0]
        };
        for (coeff, &xi) in self.lsf_coeffs[1..].iter().zip(x.iter()).take(k) {
            g += coeff * xi;
        }
        g
    }
    /// First-Order Reliability Method (FORM) reliability index β_HL.
    ///
    /// For a linear limit state g(U) = a₀ + Σ aᵢ Uᵢ in standard normal space,
    /// β = a₀ / √(Σ aᵢ²).
    ///
    /// Returns the Hasofer-Lind reliability index.
    pub fn form_beta(&self) -> f64 {
        let k = self.means.len().min(self.std_devs.len());
        let n_coeffs = self.lsf_coeffs.len();
        let a0 = if n_coeffs == 0 {
            0.0
        } else {
            let mut s = self.lsf_coeffs[0];
            for i in 0..k.min(n_coeffs.saturating_sub(1)) {
                s += self.lsf_coeffs[i + 1] * self.means[i];
            }
            s
        };
        let denom_sq: f64 = (0..k.min(n_coeffs.saturating_sub(1)))
            .map(|i| (self.lsf_coeffs[i + 1] * self.std_devs[i]).powi(2))
            .sum();
        if denom_sq < f64::EPSILON {
            return if a0 >= 0.0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            };
        }
        a0 / denom_sq.sqrt()
    }
    /// Probability of failure via FORM: Pf = Φ(−β).
    pub fn form_pf(&self) -> f64 {
        phi_cdf(-self.form_beta())
    }
    /// Second-Order Reliability Method (SORM) – Breitung correction.
    ///
    /// For a linear limit state the curvature corrections are zero,
    /// so SORM reduces to FORM. This method returns the Breitung correction
    /// factor: ∏ᵢ (1 + β κᵢ)^{-1/2} where κᵢ are principal curvatures.
    ///
    /// # Arguments
    /// * `curvatures` – principal curvatures κᵢ at the design point
    pub fn sorm_pf_breitung(&self, curvatures: &[f64]) -> f64 {
        let beta = self.form_beta();
        if !beta.is_finite() {
            return if beta > 0.0 { 0.0 } else { 1.0 };
        }
        let correction: f64 = curvatures
            .iter()
            .map(|&kappa| {
                let v = 1.0 + beta * kappa;
                if v > 0.0 { v.powf(-0.5) } else { 0.0 }
            })
            .product();
        phi_cdf(-beta) * correction
    }
    /// Monte Carlo reliability estimate: P(g(X) < 0).
    ///
    /// Draws `n_samples` realisations of X ~ N(μ, σ²) and counts failures.
    ///
    /// # Arguments
    /// * `n_samples` – number of Monte Carlo samples
    pub fn monte_carlo_pf(&self, n_samples: usize) -> f64 {
        let mut rng = rand::rng();
        let k = self.means.len().min(self.std_devs.len());
        if n_samples == 0 || k == 0 {
            return 0.0;
        }
        let mut n_fail = 0usize;
        let mut x = vec![0.0f64; k];
        for _ in 0..n_samples {
            for (x_i, (mean_i, std_i)) in x
                .iter_mut()
                .zip(self.means.iter().zip(self.std_devs.iter()))
            {
                let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
                let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                let z = box_muller(u1, u2);
                *x_i = mean_i + std_i * z;
            }
            if self.lsf(&x) < 0.0 {
                n_fail += 1;
            }
        }
        n_fail as f64 / n_samples as f64
    }
    /// Safety factor: μ_capacity / μ_demand (ratio of first two means).
    ///
    /// Returns 0 if fewer than 2 means are defined or demand mean is zero.
    pub fn safety_factor(&self) -> f64 {
        if self.means.len() < 2 || self.means[1].abs() < f64::EPSILON {
            return 0.0;
        }
        self.means[0] / self.means[1]
    }
    /// Importance vector α (direction cosines in standard normal space).
    ///
    /// αᵢ = cᵢ σᵢ / ‖c σ‖ for a linear limit state.
    pub fn importance_vector(&self) -> Vec<f64> {
        let k = self.means.len().min(self.std_devs.len());
        let n_coeffs = self.lsf_coeffs.len();
        let scaled: Vec<f64> = (0..k.min(n_coeffs.saturating_sub(1)))
            .map(|i| self.lsf_coeffs[i + 1] * self.std_devs[i])
            .collect();
        let norm: f64 = scaled.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < f64::EPSILON {
            return vec![0.0; k];
        }
        scaled.iter().map(|v| v / norm).collect()
    }
}
/// Variance-based global sensitivity analysis.
///
/// Implements Morris elementary effects screening and the Saltelli method
/// for first-order and total-order Sobol indices.
pub struct SensitivityAnalysis {
    /// Number of input variables.
    pub n_inputs: usize,
    /// Number of model evaluations per sensitivity estimate.
    pub n_samples: usize,
}
impl SensitivityAnalysis {
    /// Create a new [`SensitivityAnalysis`] for `n_inputs` variables.
    ///
    /// # Arguments
    /// * `n_inputs`  – number of model input variables
    /// * `n_samples` – number of sample pairs for Saltelli estimator
    pub fn new(n_inputs: usize, n_samples: usize) -> Self {
        Self {
            n_inputs,
            n_samples,
        }
    }
    /// Morris elementary effects screening (mean absolute μ* and σ).
    ///
    /// Evaluates the model function `f` on Morris OAT trajectories and
    /// returns `(mu_star, sigma)` vectors of length `n_inputs`.
    ///
    /// # Arguments
    /// * `f`     – model function mapping input slice → scalar output
    /// * `delta` – step size for OAT perturbation (typically 0.5/p for p levels)
    pub fn morris_screening<F>(&self, f: &F, delta: f64) -> (Vec<f64>, Vec<f64>)
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut rng = rand::rng();
        let k = self.n_inputs;
        let r = self.n_samples.max(1);
        let mut ee: Vec<Vec<f64>> = vec![Vec::new(); k];
        for _ in 0..r {
            let mut x: Vec<f64> = (0..k).map(|_| rng.random_range(0.0_f64..1.0_f64)).collect();
            let y0 = f(&x);
            let mut perm: Vec<usize> = (0..k).collect();
            for i in (1..k).rev() {
                let j = rng.random_range(0usize..=i);
                perm.swap(i, j);
            }
            for &dim in &perm {
                let x_old = x[dim];
                let step = if x[dim] + delta <= 1.0 { delta } else { -delta };
                x[dim] += step;
                let y1 = f(&x);
                let eff = (y1 - y0) / step;
                ee[dim].push(eff);
                x[dim] = x_old + step;
            }
        }
        let mu_star: Vec<f64> = ee
            .iter()
            .map(|v| {
                if v.is_empty() {
                    0.0
                } else {
                    v.iter().map(|e| e.abs()).sum::<f64>() / v.len() as f64
                }
            })
            .collect();
        let sigma: Vec<f64> = ee
            .iter()
            .map(|v| {
                if v.len() < 2 {
                    return 0.0;
                }
                let mu = v.iter().sum::<f64>() / v.len() as f64;
                (v.iter().map(|e| (e - mu).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
            })
            .collect();
        (mu_star, sigma)
    }
    /// Saltelli estimator for first-order (`s1`) and total-order (`st`) Sobol indices.
    ///
    /// Uses two independent sample matrices A and B of size `n_samples × n_inputs`.
    /// Returns `(s1, st)` vectors of length `n_inputs`.
    ///
    /// # Arguments
    /// * `f` – model function mapping input slice → scalar output
    pub fn saltelli_sobol<F>(&self, f: &F) -> (Vec<f64>, Vec<f64>)
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut rng = rand::rng();
        let k = self.n_inputs;
        let n = self.n_samples.max(1);
        let a_mat: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..k).map(|_| rng.random_range(0.0_f64..1.0_f64)).collect())
            .collect();
        let b_mat: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..k).map(|_| rng.random_range(0.0_f64..1.0_f64)).collect())
            .collect();
        let fa: Vec<f64> = a_mat.iter().map(|row| f(row)).collect();
        let fb: Vec<f64> = b_mat.iter().map(|row| f(row)).collect();
        let total_mean = (fa.iter().sum::<f64>() + fb.iter().sum::<f64>()) / (2 * n) as f64;
        let total_var: f64 = fa
            .iter()
            .chain(fb.iter())
            .map(|&y| (y - total_mean).powi(2))
            .sum::<f64>()
            / (2 * n) as f64;
        let mut s1 = vec![0.0f64; k];
        let mut st = vec![0.0f64; k];
        if total_var < f64::EPSILON {
            return (s1, st);
        }
        for j in 0..k {
            let ab: Vec<f64> = a_mat
                .iter()
                .zip(b_mat.iter())
                .map(|(a_row, b_row)| {
                    let mut row = a_row.clone();
                    row[j] = b_row[j];
                    f(&row)
                })
                .collect();
            let s1_num: f64 = fb
                .iter()
                .zip(ab.iter())
                .zip(fa.iter())
                .map(|((&yb, &yab), &ya)| yb * (yab - ya))
                .sum::<f64>()
                / n as f64;
            s1[j] = s1_num / total_var;
            let st_num: f64 = fa
                .iter()
                .zip(ab.iter())
                .map(|(&ya, &yab)| (ya - yab).powi(2))
                .sum::<f64>()
                / (2 * n) as f64;
            st[j] = st_num / total_var;
        }
        (s1, st)
    }
    /// Variance-based importance measure for a linear model with independent inputs.
    ///
    /// For a linear model y = Σ aᵢ xᵢ with xᵢ ~ N(0, σᵢ²),
    /// the first-order Sobol index is Sᵢ = (aᵢ σᵢ)² / Var\[y\].
    ///
    /// # Arguments
    /// * `coefficients` – linear model coefficients a₁, …, aₖ
    /// * `std_devs`     – input standard deviations σ₁, …, σₖ
    pub fn linear_sobol_indices(coefficients: &[f64], std_devs: &[f64]) -> Vec<f64> {
        let k = coefficients.len().min(std_devs.len());
        let total_var: f64 = (0..k)
            .map(|i| (coefficients[i] * std_devs[i]).powi(2))
            .sum();
        if total_var < f64::EPSILON {
            return vec![0.0; k];
        }
        (0..k)
            .map(|i| (coefficients[i] * std_devs[i]).powi(2) / total_var)
            .collect()
    }
}
/// Polynomial Chaos Expansion (PCE) for a scalar quantity of interest.
///
/// Stores PCE coefficients and provides methods for statistical moments,
/// Gauss-Hermite quadrature integration, and Sobol sensitivity indices.
pub struct PolynomialChaosExpansion {
    /// PCE coefficients c₀, c₁, c₂, …
    pub coefficients: Vec<f64>,
    /// Maximum polynomial order retained.
    pub max_order: usize,
}
impl PolynomialChaosExpansion {
    /// Create a new PCE from a vector of coefficients.
    ///
    /// # Arguments
    /// * `coefficients` – PCE coefficients c₀, c₁, …
    pub fn new(coefficients: Vec<f64>) -> Self {
        let max_order = coefficients.len().saturating_sub(1);
        Self {
            coefficients,
            max_order,
        }
    }
    /// Evaluate the PCE at a given standard-normal coordinate `xi`.
    pub fn evaluate(&self, xi: f64) -> f64 {
        polynomial_chaos_expansion(&self.coefficients, xi, self.max_order)
    }
    /// Mean of the quantity of interest: E\[f\] = c₀.
    pub fn mean(&self) -> f64 {
        pce_mean(&self.coefficients)
    }
    /// Variance: Var\[f\] = Σ_{n≥1} cₙ².
    pub fn variance(&self) -> f64 {
        pce_variance(&self.coefficients)
    }
    /// Standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
    /// Coefficient of variation (σ / |μ|).  Returns 0 if mean is zero.
    pub fn cov(&self) -> f64 {
        let mu = self.mean().abs();
        if mu < f64::EPSILON {
            0.0
        } else {
            self.std_dev() / mu
        }
    }
    /// Compute PCE coefficients via `n`-point Gauss-Hermite quadrature.
    ///
    /// Given a model function `f: ξ → y`, computes the first `max_order+1`
    /// Hermite-PCE coefficients by numerical integration.
    ///
    /// # Arguments
    /// * `f`         – scalar function of a standard-normal variable
    /// * `max_order` – maximum polynomial order
    /// * `quad_pts`  – number of Gauss-Hermite quadrature points (1..=5)
    pub fn from_quadrature<F>(f: F, max_order: usize, quad_pts: usize) -> Self
    where
        F: Fn(f64) -> f64,
    {
        let (nodes, weights) = gauss_hermite_quadrature(quad_pts);
        if nodes.is_empty() {
            return Self::new(vec![0.0; max_order + 1]);
        }
        let sqrt2 = 2.0_f64.sqrt();
        let mut coeffs = vec![0.0f64; max_order + 1];
        let factorial: Vec<f64> = {
            let mut v = vec![1.0f64; max_order + 2];
            for k in 1..=max_order + 1 {
                v[k] = v[k - 1] * k as f64;
            }
            v
        };
        for order in 0..=max_order {
            let mut sum = 0.0;
            for (xi_phys, &w) in nodes.iter().zip(weights.iter()) {
                let xi_prob = xi_phys / sqrt2;
                let he = hermite_polynomial(order, xi_prob);
                sum += w * f(xi_prob) * he;
            }
            let norm = factorial[order] * PI.sqrt() / 2.0_f64.powi(order as i32);
            coeffs[order] = sum / norm;
        }
        Self {
            coefficients: coeffs,
            max_order,
        }
    }
    /// First-order Sobol index for the single input variable.
    ///
    /// For a univariate PCE, the Sobol index is defined as:
    /// Sᵢ = cᵢ² / Var\[f\]  for each order i ≥ 1.
    ///
    /// Returns a vector of length `max_order` with entries S₁, S₂, …
    pub fn sobol_indices(&self) -> Vec<f64> {
        let var = self.variance();
        if var < f64::EPSILON {
            return vec![0.0; self.max_order];
        }
        self.coefficients[1..].iter().map(|c| c * c / var).collect()
    }
    /// Total Sobol index (sum of all first-order indices).
    pub fn total_sobol_index(&self) -> f64 {
        self.sobol_indices().iter().sum()
    }
}
/// Container for Monte Carlo response samples with statistics.
pub struct StochasticResponse {
    /// Raw Monte Carlo samples.
    pub samples: Vec<f64>,
}
impl StochasticResponse {
    /// Create a [`StochasticResponse`] from a vector of samples.
    pub fn new(samples: Vec<f64>) -> Self {
        Self { samples }
    }
    /// Sample mean.
    pub fn mean(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        self.samples.iter().sum::<f64>() / self.samples.len() as f64
    }
    /// Sample variance (unbiased, divided by N).
    pub fn variance(&self) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        let mu = self.mean();
        self.samples.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / n as f64
    }
    /// Empirical percentile (0 ≤ p ≤ 100).
    ///
    /// Returns the value below which `p` percent of samples fall.
    ///
    /// # Arguments
    /// * `p` – percentile in \[0, 100\]
    pub fn percentile(&self, p: f64) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        let mut sorted = self.samples.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idx = ((p / 100.0) * (n - 1) as f64).round() as usize;
        sorted[idx.min(n - 1)]
    }
    /// Empirical cumulative distribution function P(X ≤ x).
    ///
    /// # Arguments
    /// * `x` – evaluation point
    pub fn cdf(&self, x: f64) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        let count = self.samples.iter().filter(|&&s| s <= x).count();
        count as f64 / self.samples.len() as f64
    }
}
