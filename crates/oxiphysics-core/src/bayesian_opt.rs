// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bayesian optimization with Gaussian process surrogates.
//!
//! Provides a full Bayesian optimization loop: fit a GP surrogate, evaluate
//! an acquisition function to pick the next candidate, observe the objective,
//! and iterate.  Supports RBF, Matérn-5/2, and Periodic kernels.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::{PI, SQRT_2};

// ---------------------------------------------------------------------------
// Kernel
// ---------------------------------------------------------------------------

/// The covariance kernel used by the Gaussian process.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KernelType {
    /// Radial basis function (squared-exponential) kernel.
    Rbf,
    /// Matérn 5/2 kernel.
    Matern52,
    /// Periodic kernel.
    Periodic,
}

/// Hyper-parameters for the GP kernel and likelihood noise.
#[derive(Debug, Clone)]
pub struct KernelParams {
    /// Signal variance (amplitude squared).
    pub amplitude: f64,
    /// Length scale controlling smoothness.
    pub length_scale: f64,
    /// Observation noise variance added to the diagonal.
    pub noise_variance: f64,
    /// Period for the `Periodic` kernel (ignored by other kernels).
    pub period: f64,
}

impl Default for KernelParams {
    fn default() -> Self {
        Self {
            amplitude: 1.0,
            length_scale: 1.0,
            noise_variance: 1e-4,
            period: 1.0,
        }
    }
}

/// Evaluate the kernel between two input vectors `a` and `b`.
///
/// Returns the scalar covariance `k(a, b)`.
pub fn kernel_eval(kt: KernelType, params: &KernelParams, a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let sq_dist: f64 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).powi(2)).sum();
    let dist = sq_dist.sqrt();
    let l = params.length_scale.max(1e-12);
    let amp2 = params.amplitude * params.amplitude;
    match kt {
        KernelType::Rbf => amp2 * (-0.5 * sq_dist / (l * l)).exp(),
        KernelType::Matern52 => {
            let r = SQRT_2 * 5_f64.sqrt() * dist / l;
            amp2 * (1.0 + r + r * r / 3.0) * (-r).exp()
        }
        KernelType::Periodic => {
            let sin_arg = PI * dist / params.period;
            amp2 * (-2.0 * sin_arg.sin().powi(2) / (l * l)).exp()
        }
    }
}

// ---------------------------------------------------------------------------
// Cholesky helpers (used by GP)
// ---------------------------------------------------------------------------

/// Compute the lower-triangular Cholesky factor `L` such that `A = L Lᵀ`.
///
/// `A` is stored row-major in a flat `Vec`f64` of length `n*n`.
/// Returns `Err` if the matrix is not positive-definite.
pub fn cholesky(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut l = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let s: f64 = (0..j).map(|k| l[i * n + k] * l[j * n + k]).sum();
            if i == j {
                let val = a[i * n + i] - s;
                if val < 0.0 {
                    return Err(format!(
                        "matrix is not positive-definite at diagonal ({i},{i})"
                    ));
                }
                l[i * n + j] = val.sqrt();
            } else {
                let ljj = l[j * n + j];
                if ljj.abs() < 1e-15 {
                    return Err("near-zero diagonal in Cholesky".into());
                }
                l[i * n + j] = (a[i * n + j] - s) / ljj;
            }
        }
    }
    Ok(l)
}

/// Solve `L x = b` for `x` (forward substitution).  `L` is lower-triangular.
fn solve_lower(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut x = vec![0.0; n];
    for i in 0..n {
        let s: f64 = (0..i).map(|j| l[i * n + j] * x[j]).sum();
        x[i] = (b[i] - s) / l[i * n + i];
    }
    x
}

/// Solve `Lᵀ x = b` for `x` (back substitution).  `L` is lower-triangular.
fn solve_upper(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| l[j * n + i] * x[j]).sum();
        x[i] = (b[i] - s) / l[i * n + i];
    }
    x
}

// ---------------------------------------------------------------------------
// GaussianProcess
// ---------------------------------------------------------------------------

/// A Gaussian process regressor using a stationary kernel.
///
/// Training points are stored internally.  After [`GaussianProcess::fit`] the
/// GP can predict the posterior mean and variance at arbitrary test points.
#[derive(Debug, Clone)]
pub struct GaussianProcess {
    /// Kernel type.
    pub kernel: KernelType,
    /// Kernel hyper-parameters.
    pub params: KernelParams,
    /// Training inputs (each row is one sample).
    x_train: Vec<Vec<f64>>,
    /// Training targets.
    y_train: Vec<f64>,
    /// Cholesky factor of the covariance matrix (row-major, `n × n`).
    chol: Vec<f64>,
    /// α = K⁻¹ y (used for mean prediction).
    alpha: Vec<f64>,
    /// Whether the GP has been fitted.
    fitted: bool,
}

impl GaussianProcess {
    /// Construct a new, unfitted GP with the given kernel and parameters.
    pub fn new(kernel: KernelType, params: KernelParams) -> Self {
        Self {
            kernel,
            params,
            x_train: Vec::new(),
            y_train: Vec::new(),
            chol: Vec::new(),
            alpha: Vec::new(),
            fitted: false,
        }
    }

    /// Fit the GP to the given training data.
    ///
    /// # Panics
    /// Panics if `x` and `y` have different lengths, or if `x` is empty.
    pub fn fit(&mut self, x: Vec<Vec<f64>>, y: Vec<f64>) -> Result<(), String> {
        assert_eq!(x.len(), y.len(), "x and y must have the same length");
        assert!(!x.is_empty(), "Training set must be non-empty");

        let n = x.len();
        // Build the n×n kernel matrix + noise
        let mut k = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let kval = kernel_eval(self.kernel, &self.params, &x[i], &x[j]);
                k[i * n + j] = kval;
            }
            // Add noise to the diagonal
            k[i * n + i] += self.params.noise_variance;
        }

        let l = cholesky(&k, n)?;
        // Solve K α = y  →  L α' = y, then Lᵀ α = α'
        let alpha_tmp = solve_lower(&l, &y, n);
        let alpha = solve_upper(&l, &alpha_tmp, n);

        self.x_train = x;
        self.y_train = y;
        self.chol = l;
        self.alpha = alpha;
        self.fitted = true;
        Ok(())
    }

    /// Predict the posterior (mean, variance) at test point `x_star`.
    ///
    /// Returns `(mean, variance)`.
    ///
    /// # Panics
    /// Panics if the GP has not been fitted.
    pub fn predict(&self, x_star: &[f64]) -> (f64, f64) {
        assert!(self.fitted, "GP must be fitted before calling predict");
        let n = self.x_train.len();

        // k_star = [k(x*, x_1), …, k(x*, x_n)]
        let k_star: Vec<f64> = self
            .x_train
            .iter()
            .map(|xi| kernel_eval(self.kernel, &self.params, x_star, xi))
            .collect();

        // mean = k_starᵀ α
        let mean: f64 = k_star
            .iter()
            .zip(self.alpha.iter())
            .map(|(a, b)| a * b)
            .sum();

        // var = k(x*, x*) - k_starᵀ K⁻¹ k_star
        //     = k(x*, x*) - v·v  where v = L⁻¹ k_star
        let k_ss = kernel_eval(self.kernel, &self.params, x_star, x_star);
        let v = solve_lower(&self.chol, &k_star, n);
        let var = (k_ss - v.iter().map(|vi| vi * vi).sum::<f64>()).max(0.0);

        (mean, var)
    }

    /// Return the number of training points.
    pub fn n_train(&self) -> usize {
        self.x_train.len()
    }

    /// Return whether the GP has been fitted.
    pub fn is_fitted(&self) -> bool {
        self.fitted
    }
}

// ---------------------------------------------------------------------------
// Acquisition functions
// ---------------------------------------------------------------------------

/// Acquisition function used to select the next candidate point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AcquisitionFn {
    /// Expected Improvement.
    ///
    /// Balances exploitation (improving over the current best) and
    /// exploration (exploring uncertain regions).
    ExpectedImprovement,
    /// Upper Confidence Bound.
    ///
    /// Combines mean and standard deviation: `μ + κ σ`.
    UpperConfidenceBound,
    /// Probability of Improvement.
    ///
    /// Probability that the next point exceeds the current best.
    ProbabilityOfImprovement,
}

/// Standard normal CDF Φ(z).
fn standard_normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + libm_erf(z / SQRT_2))
}

/// Standard normal PDF φ(z).
fn standard_normal_pdf(z: f64) -> f64 {
    (-0.5 * z * z).exp() / (2.0 * PI).sqrt()
}

/// Error function approximation (Abramowitz & Stegun 7.1.26, max error < 1.5e-7).
fn libm_erf(x: f64) -> f64 {
    // Use the sign symmetry erf(-x) = -erf(x)
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// Evaluate the acquisition function at a single point.
///
/// - `mean`, `var` are the GP posterior moments at the candidate.
/// - `best_y`  is the best observed value so far (for EI and PI).
/// - `kappa`   is the exploration weight for UCB.
/// - `xi`      is the exploration–exploitation trade-off for EI/PI.
pub fn acquisition_value(
    acq: AcquisitionFn,
    mean: f64,
    var: f64,
    best_y: f64,
    kappa: f64,
    xi: f64,
) -> f64 {
    let sigma = var.sqrt().max(1e-12);
    match acq {
        AcquisitionFn::ExpectedImprovement => {
            let z = (mean - best_y - xi) / sigma;
            (mean - best_y - xi) * standard_normal_cdf(z) + sigma * standard_normal_pdf(z)
        }
        AcquisitionFn::UpperConfidenceBound => mean + kappa * sigma,
        AcquisitionFn::ProbabilityOfImprovement => {
            let z = (mean - best_y - xi) / sigma;
            standard_normal_cdf(z)
        }
    }
}

// ---------------------------------------------------------------------------
// Latin-hypercube sampling
// ---------------------------------------------------------------------------

/// Generate a Latin hypercube sample of `n` points in a `dim`-dimensional box.
///
/// Each axis is divided into `n` equal intervals; one point is sampled from
/// each interval per axis.  The result is a `Vec` of `n` points, each a
/// `Vec`f64` of length `dim`.  The `bounds` slice must have length `dim`.
pub fn latin_hypercube_sample(n: usize, dim: usize, bounds: &[(f64, f64)]) -> Vec<Vec<f64>> {
    assert_eq!(bounds.len(), dim, "bounds length must equal dim");
    if n == 0 || dim == 0 {
        return Vec::new();
    }

    let mut rng = rand::rng();
    use rand::RngExt as _;

    // For each dimension, create a permuted sequence of interval midpoints
    // then add a random jitter within each interval.
    let mut samples = vec![vec![0.0_f64; dim]; n];

    for d in 0..dim {
        let (lo, hi) = bounds[d];
        let interval = (hi - lo) / n as f64;

        // Create indices 0..n and shuffle them (Fisher-Yates)
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        for (i, &slot) in order.iter().enumerate() {
            let base = lo + slot as f64 * interval;
            let jitter = rng.random_range(0.0..interval);
            samples[i][d] = base + jitter;
        }
    }
    samples
}

// ---------------------------------------------------------------------------
// BayesianOptimizer
// ---------------------------------------------------------------------------

/// Configuration for [`BayesianOptimizer`].
#[derive(Debug, Clone)]
pub struct BayesOpts {
    /// Number of initial random points (LHS sample) before GP is used.
    pub n_initial: usize,
    /// Maximum number of optimization iterations.
    pub max_iter: usize,
    /// Number of random candidates evaluated to maximise the acquisition fn.
    pub n_candidates: usize,
    /// Acquisition function to use.
    pub acquisition: AcquisitionFn,
    /// UCB exploration weight κ.
    pub kappa: f64,
    /// EI/PI exploration offset ξ.
    pub xi: f64,
}

impl Default for BayesOpts {
    fn default() -> Self {
        Self {
            n_initial: 5,
            max_iter: 20,
            n_candidates: 512,
            acquisition: AcquisitionFn::ExpectedImprovement,
            kappa: 2.576,
            xi: 0.01,
        }
    }
}

/// Bayesian optimization over a bounded box using a GP surrogate.
///
/// The optimizer maintains a GP fitted to all observations so far, and at
/// each step proposes the point that maximises the acquisition function.
#[derive(Debug, Clone)]
pub struct BayesianOptimizer {
    /// Search space bounds: one `(lo, hi)` per dimension.
    pub bounds: Vec<(f64, f64)>,
    /// GP surrogate.
    pub gp: GaussianProcess,
    /// All observed inputs.
    pub x_observed: Vec<Vec<f64>>,
    /// All observed outputs.
    pub y_observed: Vec<f64>,
    /// Best output value observed so far.
    pub best_y: f64,
    /// Input that yielded `best_y`.
    pub best_x: Vec<f64>,
    /// Optimizer configuration.
    pub opts: BayesOpts,
}

impl BayesianOptimizer {
    /// Construct a new optimizer.
    ///
    /// - `bounds` — axis-aligned box, one `(lo, hi)` per input dimension.
    /// - `kernel` — kernel type for the GP surrogate.
    /// - `params` — kernel hyper-parameters.
    /// - `opts`   — algorithm options (number of iterations, acquisition, …).
    pub fn new(
        bounds: Vec<(f64, f64)>,
        kernel: KernelType,
        params: KernelParams,
        opts: BayesOpts,
    ) -> Self {
        let gp = GaussianProcess::new(kernel, params);
        Self {
            bounds,
            gp,
            x_observed: Vec::new(),
            y_observed: Vec::new(),
            best_y: f64::NEG_INFINITY,
            best_x: Vec::new(),
            opts,
        }
    }

    /// Incorporate a new observation `(x, y)` into the optimizer state.
    ///
    /// The GP surrogate is re-fitted after each call.
    pub fn update(&mut self, x: Vec<f64>, y: f64) -> Result<(), String> {
        if y > self.best_y {
            self.best_y = y;
            self.best_x = x.clone();
        }
        self.x_observed.push(x);
        self.y_observed.push(y);

        // Re-fit the GP
        self.gp
            .fit(self.x_observed.clone(), self.y_observed.clone())
    }

    /// Suggest the next candidate point to evaluate.
    ///
    /// Uses Latin-hypercube random candidates and picks the one with the
    /// highest acquisition value.  Falls back to a random LHS point if the
    /// GP has not been fitted yet.
    pub fn suggest_next(&self) -> Vec<f64> {
        let candidates =
            latin_hypercube_sample(self.opts.n_candidates, self.bounds.len(), &self.bounds);

        if !self.gp.is_fitted() {
            // Before any GP fit, just return the first candidate
            return candidates
                .into_iter()
                .next()
                .unwrap_or_else(|| self.bounds.iter().map(|(lo, hi)| (lo + hi) / 2.0).collect());
        }

        let best_y = self.best_y;
        let acq = self.opts.acquisition;
        let kappa = self.opts.kappa;
        let xi = self.opts.xi;

        let mut best_acq = f64::NEG_INFINITY;
        let mut best_candidate = candidates[0].clone();

        for cand in &candidates {
            let (mean, var) = self.gp.predict(cand);
            let val = acquisition_value(acq, mean, var, best_y, kappa, xi);
            if val > best_acq {
                best_acq = val;
                best_candidate = cand.clone();
            }
        }
        best_candidate
    }

    /// Run the full optimization loop, evaluating the black-box `f`.
    ///
    /// First draws `n_initial` LHS points, then iterates for `max_iter`
    /// steps.  Returns the best `(x, y)` pair found.
    pub fn optimize<F>(&mut self, f: F) -> (Vec<f64>, f64)
    where
        F: Fn(&[f64]) -> f64,
    {
        // --- Initial random exploration ---
        let init_samples =
            latin_hypercube_sample(self.opts.n_initial, self.bounds.len(), &self.bounds);
        for x in init_samples {
            let y = f(&x);
            let _ = self.update(x, y);
        }

        // --- Bayesian iterations ---
        for _ in 0..self.opts.max_iter {
            let x_next = self.suggest_next();
            let y_next = f(&x_next);
            let _ = self.update(x_next, y_next);
        }

        (self.best_x.clone(), self.best_y)
    }

    /// Number of observations collected so far.
    pub fn n_observations(&self) -> usize {
        self.x_observed.len()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a small 1-D GP fitted on a few points
    fn simple_gp() -> GaussianProcess {
        let mut gp = GaussianProcess::new(KernelType::Rbf, KernelParams::default());
        let x: Vec<Vec<f64>> = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]];
        let y: Vec<f64> = vec![0.0, 1.0, 0.0, -1.0];
        gp.fit(x, y).expect("fit should succeed");
        gp
    }

    // ---- Kernel tests ----

    #[test]
    fn test_rbf_kernel_same_point() {
        let p = KernelParams::default();
        let v = kernel_eval(KernelType::Rbf, &p, &[1.0, 2.0], &[1.0, 2.0]);
        // k(x,x) = amp^2 * exp(0) = 1.0
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rbf_kernel_decreases_with_distance() {
        let p = KernelParams::default();
        let k1 = kernel_eval(KernelType::Rbf, &p, &[0.0], &[0.5]);
        let k2 = kernel_eval(KernelType::Rbf, &p, &[0.0], &[1.5]);
        assert!(k1 > k2, "RBF should decrease with distance");
    }

    #[test]
    fn test_rbf_kernel_symmetry() {
        let p = KernelParams::default();
        let a = kernel_eval(KernelType::Rbf, &p, &[1.0, 2.0], &[3.0, 4.0]);
        let b = kernel_eval(KernelType::Rbf, &p, &[3.0, 4.0], &[1.0, 2.0]);
        assert!((a - b).abs() < 1e-14);
    }

    #[test]
    fn test_matern52_kernel_same_point() {
        let p = KernelParams::default();
        let v = kernel_eval(KernelType::Matern52, &p, &[0.0], &[0.0]);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_matern52_kernel_positive() {
        let p = KernelParams::default();
        let v = kernel_eval(KernelType::Matern52, &p, &[0.0], &[2.0]);
        assert!(v >= 0.0);
    }

    #[test]
    fn test_periodic_kernel_same_point() {
        let p = KernelParams::default();
        let v = kernel_eval(KernelType::Periodic, &p, &[0.0], &[0.0]);
        // sin(0) = 0 → exp(0) = 1
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_periodic_kernel_period_recovery() {
        // k(x, x+period) should equal k(x, x) for the periodic kernel
        let p = KernelParams {
            period: 2.0,
            ..Default::default()
        };
        let v0 = kernel_eval(KernelType::Periodic, &p, &[0.0], &[0.0]);
        let v1 = kernel_eval(KernelType::Periodic, &p, &[0.0], &[2.0]);
        assert!((v0 - v1).abs() < 1e-12);
    }

    // ---- Cholesky tests ----

    #[test]
    fn test_cholesky_2x2() {
        // A = [[4, 2],[2, 3]]  →  L = [[2, 0],[1, sqrt(2)]]
        let a = vec![4.0, 2.0, 2.0, 3.0];
        let l = cholesky(&a, 2).unwrap();
        assert!((l[0] - 2.0).abs() < 1e-10);
        assert!((l[1]).abs() < 1e-10);
        assert!((l[2] - 1.0).abs() < 1e-10);
        assert!((l[3] - 2_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_cholesky_identity() {
        let a = vec![1.0, 0.0, 0.0, 1.0];
        let l = cholesky(&a, 2).unwrap();
        // L should be the identity
        assert!((l[0] - 1.0).abs() < 1e-12);
        assert!((l[1]).abs() < 1e-12);
        assert!((l[2]).abs() < 1e-12);
        assert!((l[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_cholesky_not_pd_returns_err() {
        // Matrix [[-1, 0],[0, 1]] is not PD
        let a = vec![-1.0, 0.0, 0.0, 1.0];
        assert!(cholesky(&a, 2).is_err());
    }

    // ---- GP fit / predict tests ----

    #[test]
    fn test_gp_fit_succeeds() {
        let gp = simple_gp();
        assert!(gp.is_fitted());
        assert_eq!(gp.n_train(), 4);
    }

    #[test]
    fn test_gp_predict_at_training_points_close() {
        let gp = simple_gp();
        // With very small noise the posterior mean at training points should be close
        let (mean, var) = gp.predict(&[1.0]);
        assert!(
            (mean - 1.0).abs() < 0.1,
            "mean at x=1 should be ~1, got {mean}"
        );
        assert!(var >= 0.0);
    }

    #[test]
    fn test_gp_variance_nonnegative() {
        let gp = simple_gp();
        for x in [-1.0, 0.5, 1.5, 4.0] {
            let (_, var) = gp.predict(&[x]);
            assert!(var >= 0.0, "variance must be non-negative at x={x}");
        }
    }

    #[test]
    fn test_gp_variance_higher_far_from_data() {
        let gp = simple_gp();
        let (_, var_near) = gp.predict(&[1.5]);
        let (_, var_far) = gp.predict(&[100.0]);
        assert!(
            var_far > var_near,
            "variance should be higher far from training data"
        );
    }

    #[test]
    fn test_gp_fit_empty_panics() {
        let result = std::panic::catch_unwind(|| {
            let mut gp = GaussianProcess::new(KernelType::Rbf, KernelParams::default());
            let _ = gp.fit(vec![], vec![]);
        });
        assert!(result.is_err(), "fit with empty data should panic");
    }

    #[test]
    fn test_gp_matern_fit() {
        let mut gp = GaussianProcess::new(KernelType::Matern52, KernelParams::default());
        let x: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = x.iter().map(|v| v[0].sin()).collect();
        assert!(gp.fit(x, y).is_ok());
        assert!(gp.is_fitted());
    }

    #[test]
    fn test_gp_periodic_fit() {
        let p = KernelParams {
            period: PI,
            ..Default::default()
        };
        let mut gp = GaussianProcess::new(KernelType::Periodic, p);
        let x: Vec<Vec<f64>> = (0..6).map(|i| vec![i as f64 * 0.5]).collect();
        let y: Vec<f64> = x.iter().map(|v| v[0].sin()).collect();
        assert!(gp.fit(x, y).is_ok());
    }

    // ---- Acquisition function tests ----

    #[test]
    fn test_ei_nonnegative() {
        let val = acquisition_value(
            AcquisitionFn::ExpectedImprovement,
            1.5,
            0.25,
            1.0,
            2.0,
            0.01,
        );
        assert!(val >= 0.0, "EI must be non-negative");
    }

    #[test]
    fn test_ucb_increases_with_variance() {
        let low_var = acquisition_value(
            AcquisitionFn::UpperConfidenceBound,
            1.0,
            0.01,
            0.0,
            2.0,
            0.0,
        );
        let high_var =
            acquisition_value(AcquisitionFn::UpperConfidenceBound, 1.0, 1.0, 0.0, 2.0, 0.0);
        assert!(high_var > low_var, "UCB should increase with variance");
    }

    #[test]
    fn test_pi_in_unit_interval() {
        let val = acquisition_value(
            AcquisitionFn::ProbabilityOfImprovement,
            1.5,
            0.25,
            1.0,
            2.0,
            0.0,
        );
        assert!((0.0..=1.0).contains(&val), "PI must be in [0, 1]");
    }

    #[test]
    fn test_pi_zero_when_mean_below_best() {
        // If mean is much lower than best_y, PI should be ~0
        let val = acquisition_value(
            AcquisitionFn::ProbabilityOfImprovement,
            -100.0,
            0.01,
            1.0,
            2.0,
            0.0,
        );
        assert!(val < 0.01, "PI should be near 0 when far below best_y");
    }

    #[test]
    fn test_ei_zero_with_negative_improvement() {
        // mean < best_y + xi, and tiny variance → EI is effectively 0 (but ≥ 0)
        let val = acquisition_value(
            AcquisitionFn::ExpectedImprovement,
            0.0,
            1e-8,
            10.0,
            2.0,
            0.01,
        );
        assert!(val >= 0.0);
        assert!(val < 1e-3);
    }

    // ---- Latin-hypercube sampling tests ----

    #[test]
    fn test_lhs_shape() {
        let samples = latin_hypercube_sample(10, 3, &[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]);
        assert_eq!(samples.len(), 10);
        for s in &samples {
            assert_eq!(s.len(), 3);
        }
    }

    #[test]
    fn test_lhs_within_bounds() {
        let bounds = vec![(2.0, 5.0), (-1.0, 1.0)];
        let samples = latin_hypercube_sample(20, 2, &bounds);
        for s in &samples {
            assert!(s[0] >= 2.0 && s[0] <= 5.0);
            assert!(s[1] >= -1.0 && s[1] <= 1.0);
        }
    }

    #[test]
    fn test_lhs_zero_samples() {
        let s = latin_hypercube_sample(0, 2, &[(0.0, 1.0), (0.0, 1.0)]);
        assert!(s.is_empty());
    }

    #[test]
    fn test_lhs_coverage() {
        // With n=4 samples in 1D [0,4], each unit interval should be covered
        let samples = latin_hypercube_sample(4, 1, &[(0.0, 4.0)]);
        let mut covered = [false; 4];
        for s in &samples {
            let slot = (s[0] as usize).min(3);
            covered[slot] = true;
        }
        assert!(
            covered.iter().all(|&c| c),
            "each interval should be covered"
        );
    }

    // ---- BayesianOptimizer tests ----

    #[test]
    fn test_optimizer_update_increments_count() {
        let mut opt = BayesianOptimizer::new(
            vec![(0.0, 1.0)],
            KernelType::Rbf,
            KernelParams::default(),
            BayesOpts::default(),
        );
        opt.update(vec![0.5], 1.0).unwrap();
        opt.update(vec![0.7], 2.0).unwrap();
        assert_eq!(opt.n_observations(), 2);
    }

    #[test]
    fn test_optimizer_tracks_best() {
        let mut opt = BayesianOptimizer::new(
            vec![(0.0, 1.0)],
            KernelType::Rbf,
            KernelParams::default(),
            BayesOpts::default(),
        );
        opt.update(vec![0.1], 0.5).unwrap();
        opt.update(vec![0.9], 3.0).unwrap();
        opt.update(vec![0.5], 1.0).unwrap();
        assert!((opt.best_y - 3.0).abs() < 1e-12);
        assert!((opt.best_x[0] - 0.9).abs() < 1e-12);
    }

    #[test]
    fn test_optimizer_suggest_before_fit_returns_point() {
        let opt = BayesianOptimizer::new(
            vec![(0.0, 1.0), (0.0, 1.0)],
            KernelType::Rbf,
            KernelParams::default(),
            BayesOpts::default(),
        );
        let x = opt.suggest_next();
        assert_eq!(x.len(), 2);
        assert!(x[0] >= 0.0 && x[0] <= 1.0);
        assert!(x[1] >= 0.0 && x[1] <= 1.0);
    }

    #[test]
    fn test_optimizer_convergence_quadratic() {
        // Maximize f(x) = -(x - 0.3)^2 over [0, 1] — optimum at x=0.3
        let mut opt = BayesianOptimizer::new(
            vec![(0.0, 1.0)],
            KernelType::Rbf,
            KernelParams::default(),
            BayesOpts {
                n_initial: 5,
                max_iter: 15,
                n_candidates: 256,
                acquisition: AcquisitionFn::ExpectedImprovement,
                ..BayesOpts::default()
            },
        );
        let (best_x, best_y) = opt.optimize(|x| -(x[0] - 0.3).powi(2));
        // We expect to get reasonably close to the true optimum (y=0)
        assert!(
            best_y > -0.1,
            "optimizer should find near-optimum, got y={best_y}"
        );
        assert!(
            (best_x[0] - 0.3).abs() < 0.4,
            "optimizer should find x~0.3, got x={}",
            best_x[0]
        );
    }

    #[test]
    fn test_optimizer_convergence_sinusoidal() {
        // Maximize sin(2π x) over [0, 1] — global max at x=0.25
        let mut opt = BayesianOptimizer::new(
            vec![(0.0, 1.0)],
            KernelType::Rbf,
            KernelParams::default(),
            BayesOpts {
                n_initial: 8,
                max_iter: 20,
                n_candidates: 512,
                acquisition: AcquisitionFn::UpperConfidenceBound,
                kappa: 2.0,
                ..BayesOpts::default()
            },
        );
        let (_best_x, best_y) = opt.optimize(|x| (2.0 * PI * x[0]).sin());
        assert!(best_y > 0.9, "should reach sin peak, got {best_y}");
    }

    #[test]
    fn test_optimizer_2d_convergence() {
        // Maximize -(x^2 + y^2) over [-2,2]^2 — optimum at (0,0)
        let mut opt = BayesianOptimizer::new(
            vec![(-2.0, 2.0), (-2.0, 2.0)],
            KernelType::Rbf,
            KernelParams {
                length_scale: 1.5,
                ..KernelParams::default()
            },
            BayesOpts {
                n_initial: 6,
                max_iter: 20,
                n_candidates: 512,
                acquisition: AcquisitionFn::ExpectedImprovement,
                ..BayesOpts::default()
            },
        );
        let (_best_x, best_y) = opt.optimize(|x| -(x[0].powi(2) + x[1].powi(2)));
        assert!(best_y > -1.0, "should find near-origin, got y={best_y}");
    }

    #[test]
    fn test_standard_normal_cdf_symmetry() {
        // Φ(-z) = 1 - Φ(z)
        for z in [-2.0, -1.0, 0.0, 1.0, 2.0] {
            let sum = standard_normal_cdf(z) + standard_normal_cdf(-z);
            assert!((sum - 1.0).abs() < 1e-6, "CDF symmetry failed at z={z}");
        }
    }

    #[test]
    fn test_standard_normal_cdf_midpoint() {
        assert!((standard_normal_cdf(0.0) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_erf_known_values() {
        // erf(0) ≈ 0 (A&S approximation error < 1e-8)
        assert!(libm_erf(0.0).abs() < 1e-8);
        // erf(∞) ≈ 1
        assert!((libm_erf(5.0) - 1.0).abs() < 1e-5);
        // erf(-x) = -erf(x)
        assert!((libm_erf(-1.0) + libm_erf(1.0)).abs() < 1e-10);
    }
}
