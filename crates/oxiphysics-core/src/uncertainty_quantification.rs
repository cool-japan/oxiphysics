#![allow(clippy::if_same_then_else, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Uncertainty Quantification (UQ) methods for physics simulation.
//!
//! Provides:
//! - Monte Carlo uncertainty propagation
//! - Latin hypercube sampling (LHS)
//! - Sensitivity analysis (Sobol indices, Morris method)
//! - Polynomial chaos expansion (PCE)
//! - Stochastic collocation
//! - Gaussian process regression for UQ
//! - Reliability analysis (FORM / SORM)
//! - Bootstrap confidence intervals
//! - Bayesian inference via Metropolis-Hastings MCMC
//! - Probability of failure estimation

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::{FRAC_1_SQRT_2, PI};

// ─────────────────────────────────────────────────────────────────────────────
// Internal LCG RNG (no external rand dependency in tests)
// ─────────────────────────────────────────────────────────────────────────────

/// Lightweight linear congruential random number generator used internally.
struct UqRng {
    state: u64,
}

impl UqRng {
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Box-Muller standard normal sample N(0,1).
    fn next_normal(&mut self) -> f64 {
        loop {
            let u1 = self.next_f64();
            let u2 = self.next_f64();
            if u1 > 1e-15 {
                return (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
            }
        }
    }

    /// Sample from N(mean, std).
    fn next_gaussian(&mut self, mean: f64, std: f64) -> f64 {
        mean + std * self.next_normal()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper math functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the mean of a slice of values.
///
/// Returns 0.0 for an empty slice.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.iter().sum::<f64>() / data.len() as f64
}

/// Compute the variance (population) of a slice.
///
/// Returns 0.0 for slices with fewer than 2 elements.
pub fn variance(data: &[f64]) -> f64 {
    if data.len() < 2 {
        return 0.0;
    }
    let m = mean(data);
    data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / data.len() as f64
}

/// Compute the standard deviation (population) of a slice.
pub fn std_dev(data: &[f64]) -> f64 {
    variance(data).sqrt()
}

/// Sort a slice in-place and return the p-th percentile (linear interpolation).
///
/// `p` should be in \[0, 100\].
pub fn percentile(data: &mut Vec<f64>, p: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = (p / 100.0) * (data.len() - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = (lo + 1).min(data.len() - 1);
    let frac = idx - lo as f64;
    data[lo] * (1.0 - frac) + data[hi] * frac
}

/// Inverse CDF of the standard normal (rational approximation).
///
/// Maps p ∈ (0,1) to x such that Φ(x) = p.
pub fn probit(p: f64) -> f64 {
    let p = p.clamp(1e-15, 1.0 - 1e-15);
    let p_adj = if p < 0.5 { p } else { 1.0 - p };
    let t = (-2.0 * p_adj.ln()).sqrt();
    let c0 = 2.515_517;
    let c1 = 0.802_853;
    let c2 = 0.010_328;
    let d1 = 1.432_788;
    let d2 = 0.189_269;
    let d3 = 0.001_308;
    let num = c0 + c1 * t + c2 * t * t;
    let den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;
    let x = t - num / den;
    if p >= 0.5 { x } else { -x }
}

/// Standard normal CDF (error function approximation).
///
/// Returns Φ(x) = P(Z ≤ x) for Z ~ N(0,1).
pub fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / 2.0_f64.sqrt()))
}

/// Error function (Horner's method approximation).
pub fn erf(x: f64) -> f64 {
    // Abramowitz & Stegun 7.1.26
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    sign * (1.0 - poly * (-x * x).exp())
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Monte Carlo Uncertainty Propagation
// ─────────────────────────────────────────────────────────────────────────────

/// Description of a single uncertain input parameter.
#[derive(Debug, Clone)]
pub struct UncertainParameter {
    /// Mean value of the parameter.
    pub mean: f64,
    /// Standard deviation of the parameter.
    pub std: f64,
}

impl UncertainParameter {
    /// Create a new uncertain parameter with given mean and standard deviation.
    pub fn new(mean: f64, std: f64) -> Self {
        Self { mean, std }
    }
}

/// Results from a Monte Carlo uncertainty propagation run.
#[derive(Debug, Clone)]
pub struct McUqResult {
    /// Propagated output samples.
    pub samples: Vec<f64>,
    /// Sample mean of the output.
    pub mean: f64,
    /// Sample standard deviation of the output.
    pub std: f64,
    /// 5th percentile of the output.
    pub p05: f64,
    /// 95th percentile of the output.
    pub p95: f64,
}

/// Run Monte Carlo uncertainty propagation.
///
/// Draws `n_samples` realisations of `params` (each Gaussian), evaluates
/// `model(inputs)`, and returns statistics on the output distribution.
///
/// # Arguments
/// * `params` - Slice of uncertain input parameters.
/// * `n_samples` - Number of Monte Carlo samples.
/// * `model` - Function mapping a slice of inputs to a scalar output.
/// * `seed` - RNG seed for reproducibility.
pub fn monte_carlo_propagation(
    params: &[UncertainParameter],
    n_samples: usize,
    model: &dyn Fn(&[f64]) -> f64,
    seed: u64,
) -> McUqResult {
    let mut rng = UqRng::new(seed);
    let mut outputs = Vec::with_capacity(n_samples);
    let mut inputs = vec![0.0_f64; params.len()];
    for _ in 0..n_samples {
        for (i, p) in params.iter().enumerate() {
            inputs[i] = rng.next_gaussian(p.mean, p.std);
        }
        outputs.push(model(&inputs));
    }
    let m = mean(&outputs);
    let s = std_dev(&outputs);
    let mut sorted = outputs.clone();
    let p05 = percentile(&mut sorted, 5.0);
    let mut sorted2 = outputs.clone();
    let p95 = percentile(&mut sorted2, 95.0);
    McUqResult {
        samples: outputs,
        mean: m,
        std: s,
        p05,
        p95,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Latin Hypercube Sampling (LHS)
// ─────────────────────────────────────────────────────────────────────────────

/// A Latin Hypercube Sample: rows = samples, columns = dimensions.
#[derive(Debug, Clone)]
pub struct LatinHypercubeSample {
    /// Number of samples.
    pub n_samples: usize,
    /// Number of dimensions (input variables).
    pub n_dims: usize,
    /// Sample matrix, row-major: `data[i * n_dims + j]` is sample i, dim j.
    /// Values are in \[0, 1\].
    pub data: Vec<f64>,
}

impl LatinHypercubeSample {
    /// Access sample `i`, dimension `j`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n_dims + j]
    }
}

/// Generate a Latin Hypercube Sample with `n_samples` points in `n_dims` dimensions.
///
/// Each dimension is stratified into `n_samples` equal-width intervals;
/// one point is drawn uniformly from each interval, then the points are
/// randomly permuted per dimension (mid-point LHS variant).
///
/// Returns a [`LatinHypercubeSample`] with values in \[0, 1\].
pub fn latin_hypercube_sample(n_samples: usize, n_dims: usize, seed: u64) -> LatinHypercubeSample {
    let mut rng = UqRng::new(seed);
    let mut data = vec![0.0_f64; n_samples * n_dims];
    let step = 1.0 / n_samples as f64;

    for j in 0..n_dims {
        // Create stratified points for this dimension
        let mut col: Vec<f64> = (0..n_samples)
            .map(|i| (i as f64 + rng.next_f64()) * step)
            .collect();
        // Fisher-Yates shuffle
        for k in (1..n_samples).rev() {
            let idx = (rng.next_u64() as usize) % (k + 1);
            col.swap(k, idx);
        }
        for i in 0..n_samples {
            data[i * n_dims + j] = col[i];
        }
    }
    LatinHypercubeSample {
        n_samples,
        n_dims,
        data,
    }
}

/// Scale LHS samples from \[0,1\] to \[low, high\] for each dimension.
///
/// `bounds` must have exactly `lhs.n_dims` entries of `(low, high)`.
pub fn scale_lhs(lhs: &LatinHypercubeSample, bounds: &[(f64, f64)]) -> Vec<Vec<f64>> {
    assert_eq!(bounds.len(), lhs.n_dims);
    (0..lhs.n_samples)
        .map(|i| {
            (0..lhs.n_dims)
                .map(|j| {
                    let (lo, hi) = bounds[j];
                    lo + lhs.get(i, j) * (hi - lo)
                })
                .collect()
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Sensitivity Analysis
// ─────────────────────────────────────────────────────────────────────────────

// ── 3a. Morris Method ────────────────────────────────────────────────────────

/// Result of a Morris (elementary effects) sensitivity screening.
#[derive(Debug, Clone)]
pub struct MorrisResult {
    /// Mean of absolute elementary effects μ* per input dimension.
    pub mu_star: Vec<f64>,
    /// Standard deviation of elementary effects σ per input dimension.
    pub sigma: Vec<f64>,
}

/// Run the Morris method (elementary effects) for `n_dims` inputs.
///
/// Each trajectory perturbs one input at a time by `delta` and records the
/// elementary effect on the model output.  `r` trajectories are generated.
///
/// # Arguments
/// * `n_dims` - Number of input dimensions.
/// * `r` - Number of Morris trajectories.
/// * `delta` - Fractional step size (e.g. 0.1).
/// * `model` - Scalar model function.
/// * `bounds` - `(low, high)` for each dimension.
/// * `seed` - RNG seed.
pub fn morris_sensitivity(
    n_dims: usize,
    r: usize,
    delta: f64,
    model: &dyn Fn(&[f64]) -> f64,
    bounds: &[(f64, f64)],
    seed: u64,
) -> MorrisResult {
    let mut rng = UqRng::new(seed);
    let mut all_effects: Vec<Vec<f64>> = vec![Vec::new(); n_dims];

    for _ in 0..r {
        // Random starting point
        let mut x: Vec<f64> = bounds
            .iter()
            .map(|(lo, hi)| lo + rng.next_f64() * (hi - lo))
            .collect();
        // Random permutation of dimensions
        let mut perm: Vec<usize> = (0..n_dims).collect();
        for k in (1..n_dims).rev() {
            let idx = (rng.next_u64() as usize) % (k + 1);
            perm.swap(k, idx);
        }
        let mut y0 = model(&x);
        for &dim in &perm {
            let (lo, hi) = bounds[dim];
            let step = delta * (hi - lo);
            // Ensure we stay within bounds
            let new_val = if x[dim] + step <= hi {
                x[dim] + step
            } else {
                x[dim] - step
            };
            let old_val = x[dim];
            x[dim] = new_val;
            let y1 = model(&x);
            let ee = (y1 - y0) / (new_val - old_val);
            all_effects[dim].push(ee);
            y0 = y1;
        }
    }

    let mu_star: Vec<f64> = all_effects
        .iter()
        .map(|effects| {
            if effects.is_empty() {
                0.0
            } else {
                effects.iter().map(|e| e.abs()).sum::<f64>() / effects.len() as f64
            }
        })
        .collect();

    let sigma: Vec<f64> = all_effects
        .iter()
        .map(|effects| {
            if effects.len() < 2 {
                return 0.0;
            }
            let m = effects.iter().sum::<f64>() / effects.len() as f64;
            let v =
                effects.iter().map(|e| (e - m).powi(2)).sum::<f64>() / (effects.len() - 1) as f64;
            v.sqrt()
        })
        .collect();

    MorrisResult { mu_star, sigma }
}

// ── 3b. Sobol Indices ────────────────────────────────────────────────────────

/// First-order and total-order Sobol sensitivity indices.
#[derive(Debug, Clone)]
pub struct SobolIndices {
    /// First-order Sobol index S_i for each input dimension.
    pub first_order: Vec<f64>,
    /// Total-order Sobol index T_i for each input dimension.
    pub total_order: Vec<f64>,
}

/// Estimate Sobol sensitivity indices using the Saltelli (2010) estimator.
///
/// Uses `n_base` quasi-random base samples (LHS-generated) and the Jansen
/// estimator for total-order indices.
///
/// # Arguments
/// * `n_dims` - Number of input variables.
/// * `n_base` - Base sample size N; total model evaluations ≈ N·(n_dims+2).
/// * `model` - Scalar model function mapping inputs → output.
/// * `bounds` - `(low, high)` per dimension.
/// * `seed` - RNG seed.
pub fn sobol_indices(
    n_dims: usize,
    n_base: usize,
    model: &dyn Fn(&[f64]) -> f64,
    bounds: &[(f64, f64)],
    seed: u64,
) -> SobolIndices {
    // Generate two independent LHS sample matrices A and B
    let lhs_a = latin_hypercube_sample(n_base, n_dims, seed);
    let lhs_b = latin_hypercube_sample(n_base, n_dims, seed.wrapping_add(0xDEAD_BEEF));

    let a_mat = scale_lhs(&lhs_a, bounds);
    let b_mat = scale_lhs(&lhs_b, bounds);

    // Evaluate model on A and B
    let fa: Vec<f64> = a_mat.iter().map(|x| model(x)).collect();
    let fb: Vec<f64> = b_mat.iter().map(|x| model(x)).collect();

    let var_y = {
        let all: Vec<f64> = fa.iter().chain(fb.iter()).copied().collect();
        variance(&all)
    };

    let mut first_order = vec![0.0_f64; n_dims];
    let mut total_order = vec![0.0_f64; n_dims];

    let n = n_base as f64;

    for dim in 0..n_dims {
        // Build A_B^(i): A with column i replaced by B's column i
        let ab_i: Vec<Vec<f64>> = (0..n_base)
            .map(|k| {
                let mut row = a_mat[k].clone();
                row[dim] = b_mat[k][dim];
                row
            })
            .collect();

        let fab_i: Vec<f64> = ab_i.iter().map(|x| model(x)).collect();

        // Saltelli (2010) first-order estimator: S_i = (1/N) Σ f_B(f_AB_i - f_A) / Var
        let s1_num: f64 = (0..n_base).map(|k| fb[k] * (fab_i[k] - fa[k])).sum::<f64>() / n;

        // Jansen total-order estimator: T_i = (1/(2N)) Σ (f_A - f_AB_i)² / Var
        let st_num: f64 = (0..n_base).map(|k| (fa[k] - fab_i[k]).powi(2)).sum::<f64>() / (2.0 * n);

        first_order[dim] = if var_y > 1e-15 { s1_num / var_y } else { 0.0 };
        total_order[dim] = if var_y > 1e-15 { st_num / var_y } else { 0.0 };
    }

    SobolIndices {
        first_order,
        total_order,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Polynomial Chaos Expansion (PCE)
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the n-th Hermite polynomial He_n(x) (probabilists' convention).
///
/// He_0 = 1, He_1 = x, He_{n+1} = x·He_n − n·He_{n−1}.
pub fn hermite_poly(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut h_prev = 1.0;
            let mut h_curr = x;
            for k in 1..n {
                let h_next = x * h_curr - k as f64 * h_prev;
                h_prev = h_curr;
                h_curr = h_next;
            }
            h_curr
        }
    }
}

/// Polynomial Chaos Expansion surrogate model.
///
/// Represents an output Y as a linear combination of Hermite basis polynomials
/// evaluated at standardised Gaussian inputs.
#[derive(Debug, Clone)]
pub struct PolynomialChaosExpansion {
    /// PCE coefficients, one per multi-index.
    pub coeffs: Vec<f64>,
    /// Multi-indices: each entry is a vector of polynomial degrees per dimension.
    pub multi_indices: Vec<Vec<usize>>,
    /// Number of input dimensions.
    pub n_dims: usize,
}

impl PolynomialChaosExpansion {
    /// Evaluate the PCE at a point `x` (standardised Gaussian inputs).
    ///
    /// `x.len()` must equal `self.n_dims`.
    pub fn evaluate(&self, x: &[f64]) -> f64 {
        assert_eq!(x.len(), self.n_dims);
        self.coeffs
            .iter()
            .zip(self.multi_indices.iter())
            .map(|(c, idx)| {
                let basis: f64 = idx
                    .iter()
                    .zip(x.iter())
                    .map(|(&deg, &xi)| hermite_poly(deg, xi))
                    .product();
                c * basis
            })
            .sum()
    }

    /// Mean of the PCE (= coefficient of the zero-order term).
    pub fn pce_mean(&self) -> f64 {
        // The zero multi-index has all degrees = 0; He_0 = 1, <He_0²> = 1
        self.coeffs
            .iter()
            .zip(self.multi_indices.iter())
            .find(|(_, idx)| idx.iter().all(|&d| d == 0))
            .map(|(c, _)| *c)
            .unwrap_or(0.0)
    }

    /// Variance of the PCE (= sum of squares of non-zero-order coefficients).
    ///
    /// Uses the orthogonality of Hermite polynomials: ⟨He_α, He_α⟩ = α!.
    pub fn pce_variance(&self) -> f64 {
        self.coeffs
            .iter()
            .zip(self.multi_indices.iter())
            .filter(|(_, idx)| !idx.iter().all(|&d| d == 0))
            .map(|(c, idx)| {
                let norm_sq: f64 = idx.iter().map(|&d| factorial(d) as f64).product();
                c * c * norm_sq
            })
            .sum()
    }
}

fn factorial(n: usize) -> u64 {
    (1..=n as u64).product()
}

/// Generate all multi-indices of total degree ≤ `max_degree` for `n_dims` dimensions.
pub fn pce_multi_indices(n_dims: usize, max_degree: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    pce_multi_indices_rec(n_dims, max_degree, 0, &mut Vec::new(), &mut result);
    result
}

fn pce_multi_indices_rec(
    n_dims: usize,
    remaining: usize,
    dim: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if dim == n_dims {
        result.push(current.clone());
        return;
    }
    let max_here = if dim == n_dims - 1 {
        remaining
    } else {
        remaining
    };
    for d in 0..=max_here {
        current.push(d);
        if dim + 1 < n_dims {
            pce_multi_indices_rec(n_dims, remaining - d, dim + 1, current, result);
        } else {
            result.push(current.clone());
        }
        current.pop();
    }
}

/// Fit a PCE by non-intrusive (regression) approach using LHS-sampled training data.
///
/// Solves the least-squares problem Ψ c = Y for the coefficient vector c.
///
/// # Arguments
/// * `n_dims` - Number of Gaussian input dimensions.
/// * `max_degree` - Maximum total polynomial degree.
/// * `n_train` - Number of training samples (recommend ≥ 2 × n_basis).
/// * `model` - Model function on standardised Gaussian inputs.
/// * `seed` - RNG seed.
pub fn fit_pce(
    n_dims: usize,
    max_degree: usize,
    n_train: usize,
    model: &dyn Fn(&[f64]) -> f64,
    seed: u64,
) -> PolynomialChaosExpansion {
    let multi_indices = pce_multi_indices(n_dims, max_degree);
    let n_basis = multi_indices.len();

    // Generate training points (standard normal via probit transform of LHS)
    let lhs = latin_hypercube_sample(n_train, n_dims, seed);

    let mut x_train: Vec<Vec<f64>> = (0..n_train)
        .map(|i| {
            (0..n_dims)
                .map(|j| probit(lhs.get(i, j).clamp(1e-6, 1.0 - 1e-6)))
                .collect()
        })
        .collect();

    let y_train: Vec<f64> = x_train.iter_mut().map(|x| model(x)).collect();

    // Build Vandermonde-like matrix Ψ: n_train × n_basis
    let mut psi = vec![0.0_f64; n_train * n_basis];
    for (i, xi) in x_train.iter().enumerate() {
        for (j, idx) in multi_indices.iter().enumerate() {
            let basis: f64 = idx
                .iter()
                .zip(xi.iter())
                .map(|(&deg, &x)| hermite_poly(deg, x))
                .product();
            psi[i * n_basis + j] = basis;
        }
    }

    // Solve normal equations Ψᵀ Ψ c = Ψᵀ y (simple, not production quality)
    let coeffs = normal_equations_solve(&psi, &y_train, n_train, n_basis);

    PolynomialChaosExpansion {
        coeffs,
        multi_indices,
        n_dims,
    }
}

/// Solve the normal equations A^T A x = A^T b using Cholesky or Gaussian elimination.
fn normal_equations_solve(a: &[f64], b: &[f64], m: usize, n: usize) -> Vec<f64> {
    // Compute Aᵀ A (n×n) and Aᵀ b (n)
    let mut ata = vec![0.0_f64; n * n];
    let mut atb = vec![0.0_f64; n];
    for i in 0..m {
        for j in 0..n {
            atb[j] += a[i * n + j] * b[i];
            for k in 0..n {
                ata[j * n + k] += a[i * n + j] * a[i * n + k];
            }
        }
    }
    // Add small regularisation
    for j in 0..n {
        ata[j * n + j] += 1e-10;
    }
    // Forward Gaussian elimination with partial pivoting
    let mut x = atb.clone();
    let mut mat = ata;
    for col in 0..n {
        // Find pivot
        let mut max_val = mat[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            if mat[row * n + col].abs() > max_val {
                max_val = mat[row * n + col].abs();
                max_row = row;
            }
        }
        if max_row != col {
            for k in 0..n {
                mat.swap(col * n + k, max_row * n + k);
            }
            x.swap(col, max_row);
        }
        let pivot = mat[col * n + col];
        if pivot.abs() < 1e-18 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = mat[row * n + col] / pivot;
            for k in col..n {
                let tmp = mat[col * n + k] * factor;
                mat[row * n + k] -= tmp;
            }
            x[row] -= x[col] * factor;
        }
    }
    // Back substitution
    for col in (0..n).rev() {
        let pivot = mat[col * n + col];
        if pivot.abs() < 1e-18 {
            continue;
        }
        x[col] /= pivot;
        for row in 0..col {
            let tmp = mat[row * n + col] * x[col];
            x[row] -= tmp;
        }
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Stochastic Collocation
// ─────────────────────────────────────────────────────────────────────────────

/// A collocation node (quadrature point + weight) for one dimension.
#[derive(Debug, Clone, Copy)]
pub struct CollocationNode {
    /// Quadrature abscissa (standardised).
    pub point: f64,
    /// Quadrature weight.
    pub weight: f64,
}

/// Gauss-Hermite quadrature nodes and weights for n-point rule.
///
/// Suitable for integrals of the form ∫ f(x) exp(−x²) dx.
/// Returns nodes sorted in ascending order.
///
/// Supported orders: 1, 2, 3, 4, 5.
pub fn gauss_hermite_nodes(order: usize) -> Vec<CollocationNode> {
    // Tabulated values from Abramowitz & Stegun
    match order {
        1 => vec![CollocationNode {
            point: 0.0,
            weight: 1.772_453_851_f64, // √π ≈ 1.7724538509055159
        }],
        2 => vec![
            CollocationNode {
                point: -FRAC_1_SQRT_2,
                weight: 0.886_226_925,
            },
            CollocationNode {
                point: FRAC_1_SQRT_2,
                weight: 0.886_226_925,
            },
        ],
        3 => vec![
            CollocationNode {
                point: -1.224_744_871,
                weight: 0.295_408_975,
            },
            CollocationNode {
                point: 0.0,
                weight: 1.181_635_901,
            },
            CollocationNode {
                point: 1.224_744_871,
                weight: 0.295_408_975,
            },
        ],
        4 => vec![
            CollocationNode {
                point: -1.650_680_123,
                weight: 0.081_312_835,
            },
            CollocationNode {
                point: -0.524_647_623,
                weight: 0.804_914_090,
            },
            CollocationNode {
                point: 0.524_647_623,
                weight: 0.804_914_090,
            },
            CollocationNode {
                point: 1.650_680_123,
                weight: 0.081_312_835,
            },
        ],
        5 => vec![
            CollocationNode {
                point: -2.020_182_470,
                weight: 0.019_953_242,
            },
            CollocationNode {
                point: -0.958_572_464,
                weight: 0.394_424_314,
            },
            CollocationNode {
                point: 0.0,
                weight: 0.945_308_722,
            },
            CollocationNode {
                point: 0.958_572_464,
                weight: 0.394_424_314,
            },
            CollocationNode {
                point: 2.020_182_470,
                weight: 0.019_953_242,
            },
        ],
        _ => gauss_hermite_nodes(3), // fallback
    }
}

/// Compute mean and variance of a model via tensor-product stochastic collocation.
///
/// Integrates E\[Y\] and E\[Y²\] using Gauss-Hermite quadrature for Gaussian inputs.
///
/// # Arguments
/// * `n_dims` - Number of independent Gaussian inputs.
/// * `order` - Gauss-Hermite quadrature order per dimension (1–5).
/// * `model` - Model on standardised Gaussian inputs.
pub fn stochastic_collocation(
    n_dims: usize,
    order: usize,
    model: &dyn Fn(&[f64]) -> f64,
) -> (f64, f64) {
    let nodes_1d = gauss_hermite_nodes(order);
    let _n_pts_1d = nodes_1d.len();

    // Compute total weight (normalise by π^(n/2) for standard normal measure)
    let pi_n_half = PI.powf(n_dims as f64 / 2.0);

    // Enumerate all multi-index combinations via recursion
    let mut mean_acc = 0.0_f64;
    let mut sq_acc = 0.0_f64;
    let mut x = vec![0.0_f64; n_dims];
    collocation_recurse(
        &nodes_1d,
        n_dims,
        0,
        1.0,
        &mut x,
        model,
        &mut mean_acc,
        &mut sq_acc,
    );

    mean_acc /= pi_n_half;
    sq_acc /= pi_n_half;
    let var = (sq_acc - mean_acc * mean_acc).max(0.0);
    (mean_acc, var)
}

#[allow(clippy::too_many_arguments)]
fn collocation_recurse(
    nodes: &[CollocationNode],
    n_dims: usize,
    dim: usize,
    w_acc: f64,
    x: &mut Vec<f64>,
    model: &dyn Fn(&[f64]) -> f64,
    mean_acc: &mut f64,
    sq_acc: &mut f64,
) {
    if dim == n_dims {
        let y = model(x);
        // Gauss-Hermite weights already encode the e^{-x^2} factor; use directly.
        *mean_acc += w_acc * y;
        *sq_acc += w_acc * y * y;
        return;
    }
    for node in nodes {
        x[dim] = node.point * 2.0_f64.sqrt(); // scale to N(0,1)
        collocation_recurse(
            nodes,
            n_dims,
            dim + 1,
            w_acc * node.weight,
            x,
            model,
            mean_acc,
            sq_acc,
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Gaussian Process Regression for UQ
// ─────────────────────────────────────────────────────────────────────────────

/// Squared-exponential (RBF) kernel: k(x, x') = σ² exp(-‖x−x'‖² / (2ℓ²)).
///
/// # Arguments
/// * `x` - First input vector.
/// * `xp` - Second input vector.
/// * `sigma_f` - Signal standard deviation σ_f.
/// * `length_scale` - Length scale ℓ.
pub fn rbf_kernel(x: &[f64], xp: &[f64], sigma_f: f64, length_scale: f64) -> f64 {
    let sq_dist: f64 = x.iter().zip(xp.iter()).map(|(a, b)| (a - b).powi(2)).sum();
    sigma_f * sigma_f * (-sq_dist / (2.0 * length_scale * length_scale)).exp()
}

/// A fitted Gaussian Process regression model.
#[derive(Debug, Clone)]
pub struct GaussianProcessSurrogate {
    /// Training input points (row-major: n_train × n_dims).
    pub x_train: Vec<Vec<f64>>,
    /// Training output values.
    pub y_train: Vec<f64>,
    /// Alpha vector = (K + σ²_n I)^{-1} y  (pre-computed for prediction).
    pub alpha: Vec<f64>,
    /// Signal standard deviation σ_f.
    pub sigma_f: f64,
    /// Length scale ℓ.
    pub length_scale: f64,
    /// Noise standard deviation σ_n.
    pub noise_std: f64,
}

impl GaussianProcessSurrogate {
    /// Predict the posterior mean and variance at a new point `x_star`.
    pub fn predict(&self, x_star: &[f64]) -> (f64, f64) {
        let k_star: Vec<f64> = self
            .x_train
            .iter()
            .map(|xi| rbf_kernel(xi, x_star, self.sigma_f, self.length_scale))
            .collect();

        let mean_pred: f64 = k_star
            .iter()
            .zip(self.alpha.iter())
            .map(|(k, a)| k * a)
            .sum();

        // Predictive variance (simplified diagonal version via K_star_star - k*^T K^{-1} k*)
        let k_ss = rbf_kernel(x_star, x_star, self.sigma_f, self.length_scale);
        // We approximate K^{-1} k* via the already-computed (K + σ²I)^{-1} alpha approach
        // For a full implementation use Cholesky; here we use the diagonal approximation.
        let var_pred = (k_ss
            - k_star.iter().map(|k| k * k).sum::<f64>() / (self.sigma_f * self.sigma_f + 1e-8))
            .max(0.0);

        (mean_pred, var_pred)
    }
}

/// Fit a Gaussian Process surrogate using exact GP regression.
///
/// Solves (K + σ²_n I) α = y via simple Gaussian elimination.
///
/// # Arguments
/// * `x_train` - Training inputs (one vector per sample).
/// * `y_train` - Training outputs.
/// * `sigma_f` - Signal standard deviation.
/// * `length_scale` - RBF length scale.
/// * `noise_std` - Observation noise standard deviation.
pub fn fit_gp(
    x_train: Vec<Vec<f64>>,
    y_train: Vec<f64>,
    sigma_f: f64,
    length_scale: f64,
    noise_std: f64,
) -> GaussianProcessSurrogate {
    let n = x_train.len();
    // Build K + σ²_n I
    let mut k = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            k[i * n + j] = rbf_kernel(&x_train[i], &x_train[j], sigma_f, length_scale);
        }
        k[i * n + i] += noise_std * noise_std;
    }
    // Solve Kα = y
    let alpha = gaussian_eliminate(&k, &y_train, n);
    GaussianProcessSurrogate {
        x_train,
        y_train,
        alpha,
        sigma_f,
        length_scale,
        noise_std,
    }
}

/// Solve A x = b via Gaussian elimination with partial pivoting.
fn gaussian_eliminate(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut max_val = mat[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            if mat[row * n + col].abs() > max_val {
                max_val = mat[row * n + col].abs();
                max_row = row;
            }
        }
        if max_row != col {
            for k in 0..n {
                mat.swap(col * n + k, max_row * n + k);
            }
            rhs.swap(col, max_row);
        }
        let pivot = mat[col * n + col];
        if pivot.abs() < 1e-18 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = mat[row * n + col] / pivot;
            for k in col..n {
                let tmp = mat[col * n + k] * factor;
                mat[row * n + k] -= tmp;
            }
            rhs[row] -= rhs[col] * factor;
        }
    }
    for col in (0..n).rev() {
        let pivot = mat[col * n + col];
        if pivot.abs() < 1e-18 {
            continue;
        }
        rhs[col] /= pivot;
        for row in 0..col {
            let tmp = mat[row * n + col] * rhs[col];
            rhs[row] -= tmp;
        }
    }
    rhs
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. Reliability Analysis (FORM / SORM)
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a FORM reliability analysis.
#[derive(Debug, Clone)]
pub struct FormResult {
    /// Most probable point (MPP) in standard normal space.
    pub beta_u: Vec<f64>,
    /// Reliability index β = ‖u*‖.
    pub beta: f64,
    /// FORM probability of failure P_f ≈ Φ(−β).
    pub pf_form: f64,
}

/// First-Order Reliability Method (FORM) via Hasofer-Lind iteration.
///
/// Finds the most probable point (MPP) on the limit-state surface g(u) = 0
/// in the standard normal space by iterative gradient projection.
///
/// # Arguments
/// * `n_dims` - Number of standard normal input dimensions.
/// * `g` - Limit-state function: g(u) < 0 means failure.
/// * `max_iter` - Maximum number of HL-RF iterations.
/// * `tol` - Convergence tolerance on β.
/// * `h` - Finite-difference step for gradient.
pub fn form_analysis(
    n_dims: usize,
    g: &dyn Fn(&[f64]) -> f64,
    max_iter: usize,
    tol: f64,
    h: f64,
) -> FormResult {
    let mut u = vec![0.0_f64; n_dims]; // start at origin
    let mut beta = 0.0_f64;

    for _ in 0..max_iter {
        let g0 = g(&u);
        // Numerical gradient
        let mut grad = vec![0.0_f64; n_dims];
        for k in 0..n_dims {
            let mut u_plus = u.clone();
            u_plus[k] += h;
            grad[k] = (g(&u_plus) - g0) / h;
        }
        let grad_norm: f64 = grad.iter().map(|x| x * x).sum::<f64>().sqrt();
        if grad_norm < 1e-15 {
            break;
        }
        // HL-RF update: u_{k+1} = (∇g · u − g) / |∇g|² · ∇g
        let dot_grad_u: f64 = grad.iter().zip(u.iter()).map(|(g, uk)| g * uk).sum();
        let factor = (dot_grad_u - g0) / (grad_norm * grad_norm);
        let u_new: Vec<f64> = grad.iter().map(|gi| factor * gi).collect();
        let new_beta: f64 = u_new.iter().map(|x| x * x).sum::<f64>().sqrt();
        if (new_beta - beta).abs() < tol {
            u = u_new;
            beta = new_beta;
            break;
        }
        beta = new_beta;
        u = u_new;
    }

    let pf = normal_cdf(-beta);
    FormResult {
        beta_u: u,
        beta,
        pf_form: pf,
    }
}

/// Second-Order Reliability Method (SORM) correction factor.
///
/// Applies the Breitung correction to the FORM result.  Requires principal
/// curvatures κ_i at the MPP (positive = convex toward safe domain).
///
/// Returns the SORM-corrected probability of failure.
pub fn sorm_breitung(beta: f64, curvatures: &[f64]) -> f64 {
    let correction: f64 = curvatures
        .iter()
        .map(|kappa| {
            let denom = (1.0 + beta * kappa).max(1e-10);
            1.0 / denom.sqrt()
        })
        .product();
    normal_cdf(-beta) * correction
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. Bootstrap Confidence Intervals
// ─────────────────────────────────────────────────────────────────────────────

/// Bootstrap percentile confidence interval for a statistic.
///
/// Resamples `data` with replacement `n_boot` times, computes the statistic
/// on each resample, and returns the \[alpha/2, 1−alpha/2\] percentile interval.
///
/// # Arguments
/// * `data` - Observed data.
/// * `statistic` - Function computing the statistic from a data slice.
/// * `n_boot` - Number of bootstrap resamples.
/// * `alpha` - Significance level (e.g. 0.05 for 95% CI).
/// * `seed` - RNG seed.
///
/// Returns `(lower_bound, upper_bound)`.
pub fn bootstrap_ci(
    data: &[f64],
    statistic: &dyn Fn(&[f64]) -> f64,
    n_boot: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    let n = data.len();
    let mut rng = UqRng::new(seed);
    let mut boot_stats: Vec<f64> = Vec::with_capacity(n_boot);
    let mut sample = vec![0.0_f64; n];
    for _ in 0..n_boot {
        for s in sample.iter_mut() {
            *s = data[(rng.next_u64() as usize) % n];
        }
        boot_stats.push(statistic(&sample));
    }
    let lo = percentile(&mut boot_stats.clone(), alpha / 2.0 * 100.0);
    let hi = percentile(&mut boot_stats, (1.0 - alpha / 2.0) * 100.0);
    (lo, hi)
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. MCMC — Metropolis-Hastings Bayesian Inference
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for the Metropolis-Hastings sampler.
#[derive(Debug, Clone)]
pub struct MhConfig {
    /// Number of MCMC samples to draw.
    pub n_samples: usize,
    /// Burn-in period (samples discarded).
    pub burn_in: usize,
    /// Proposal step standard deviation.
    pub step_size: f64,
    /// RNG seed.
    pub seed: u64,
}

impl Default for MhConfig {
    fn default() -> Self {
        Self {
            n_samples: 5000,
            burn_in: 500,
            step_size: 0.1,
            seed: 42,
        }
    }
}

/// Result of Metropolis-Hastings sampling.
#[derive(Debug, Clone)]
pub struct MhResult {
    /// Posterior samples (post burn-in), shape = n_kept × n_dims.
    pub samples: Vec<Vec<f64>>,
    /// Overall acceptance rate.
    pub acceptance_rate: f64,
}

/// Metropolis-Hastings MCMC sampler.
///
/// Draws samples from a target (unnormalised) log-posterior.
///
/// # Arguments
/// * `log_posterior` - Function computing log p(θ | data) up to a constant.
/// * `initial` - Starting point in parameter space.
/// * `config` - Sampler configuration.
pub fn metropolis_hastings(
    log_posterior: &dyn Fn(&[f64]) -> f64,
    initial: &[f64],
    config: &MhConfig,
) -> MhResult {
    let _n_dims = initial.len();
    let mut rng = UqRng::new(config.seed);
    let mut current = initial.to_vec();
    let mut log_p_current = log_posterior(&current);
    let total = config.n_samples + config.burn_in;
    let mut samples = Vec::with_capacity(config.n_samples);
    let mut n_accepted = 0_usize;

    for step in 0..total {
        // Random-walk Gaussian proposal
        let proposal: Vec<f64> = current
            .iter()
            .map(|x| x + rng.next_normal() * config.step_size)
            .collect();
        let log_p_proposal = log_posterior(&proposal);
        let log_accept = log_p_proposal - log_p_current;
        let u = rng.next_f64();
        if log_accept >= 0.0 || u < log_accept.exp() {
            current = proposal;
            log_p_current = log_p_proposal;
            if step >= config.burn_in {
                n_accepted += 1;
            }
        }
        if step >= config.burn_in {
            samples.push(current.clone());
        }
    }
    let acceptance_rate = n_accepted as f64 / config.n_samples as f64;
    MhResult {
        samples,
        acceptance_rate,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. Probability of Failure Estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the probability of failure P_f = P(g(X) ≤ 0) via crude Monte Carlo.
///
/// # Arguments
/// * `params` - Input uncertain parameters.
/// * `limit_state` - Function g(x); failure when g(x) ≤ 0.
/// * `n_samples` - Number of MC samples.
/// * `seed` - RNG seed.
///
/// Returns `(pf_estimate, coefficient_of_variation)`.
pub fn probability_of_failure_mc(
    params: &[UncertainParameter],
    limit_state: &dyn Fn(&[f64]) -> f64,
    n_samples: usize,
    seed: u64,
) -> (f64, f64) {
    let mut rng = UqRng::new(seed);
    let mut n_fail = 0_usize;
    let mut inputs = vec![0.0_f64; params.len()];
    for _ in 0..n_samples {
        for (i, p) in params.iter().enumerate() {
            inputs[i] = rng.next_gaussian(p.mean, p.std);
        }
        if limit_state(&inputs) <= 0.0 {
            n_fail += 1;
        }
    }
    let pf = n_fail as f64 / n_samples as f64;
    let cov = if pf > 0.0 {
        ((1.0 - pf) / (pf * n_samples as f64)).sqrt()
    } else {
        f64::INFINITY
    };
    (pf, cov)
}

/// Importance sampling estimate of P_f with a shifted Gaussian importance density.
///
/// Shifts the sampling distribution to concentrate near the MPP `u_star` in
/// standard normal space.
///
/// # Arguments
/// * `n_dims` - Number of standard normal input dimensions.
/// * `limit_state_std` - Limit-state function in standard normal space.
/// * `u_star` - Most probable point (MPP) from FORM.
/// * `n_samples` - Number of importance samples.
/// * `seed` - RNG seed.
pub fn importance_sampling_pf(
    n_dims: usize,
    limit_state_std: &dyn Fn(&[f64]) -> f64,
    u_star: &[f64],
    n_samples: usize,
    seed: u64,
) -> f64 {
    let mut rng = UqRng::new(seed);
    let mut sum_w = 0.0_f64;
    let mut sum_iw = 0.0_f64;
    let mut u = vec![0.0_f64; n_dims];
    for _ in 0..n_samples {
        // Sample from N(u_star, I)
        for k in 0..n_dims {
            u[k] = u_star[k] + rng.next_normal();
        }
        // Likelihood ratio w(u) = φ(u) / q(u) = exp(-|u|²/2) / exp(-|u−u*|²/2)
        let norm_u_sq: f64 = u.iter().map(|x| x * x).sum();
        let norm_diff_sq: f64 = u
            .iter()
            .zip(u_star.iter())
            .map(|(ui, ui_star)| (ui - ui_star).powi(2))
            .sum();
        let log_w = -0.5 * norm_u_sq + 0.5 * norm_diff_sq;
        let w = log_w.exp();
        let indicator = if limit_state_std(&u) <= 0.0 { 1.0 } else { 0.0 };
        sum_iw += indicator * w;
        sum_w += w;
    }
    if sum_w < 1e-30 { 0.0 } else { sum_iw / sum_w }
}

/// Compute the reliability index β from a probability of failure estimate.
///
/// β = −Φ⁻¹(P_f).
pub fn reliability_index(pf: f64) -> f64 {
    -probit(pf.clamp(1e-15, 1.0 - 1e-15))
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. Moment-independent sensitivity indices (Delta indices)
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the moment-independent delta sensitivity index for dimension `dim`.
///
/// Δ_i = 0.5 · E\[|f_Y(y) − f_{Y|Xi}(y)|\] approximated via binning.
///
/// # Arguments
/// * `x_samples` - Matrix of input samples (n_samples × n_dims).
/// * `y_samples` - Corresponding output samples.
/// * `dim` - Input dimension to analyse.
/// * `n_bins` - Number of bins for density estimation.
pub fn delta_sensitivity(
    x_samples: &[Vec<f64>],
    y_samples: &[f64],
    dim: usize,
    n_bins: usize,
) -> f64 {
    let n = y_samples.len();
    if n == 0 || n_bins == 0 {
        return 0.0;
    }

    let y_min = y_samples.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = y_samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = (y_max - y_min).max(1e-15);
    let bin_w = range / n_bins as f64;

    // Marginal density (histogram)
    let mut f_y = vec![0.0_f64; n_bins];
    for &yi in y_samples {
        let bin = ((yi - y_min) / bin_w).floor() as usize;
        f_y[bin.min(n_bins - 1)] += 1.0 / (n as f64 * bin_w);
    }

    // Get unique values of x_dim and stratify
    let x_dim: Vec<f64> = x_samples.iter().map(|row| row[dim]).collect();
    let x_min = x_dim.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = x_dim.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let x_range = (x_max - x_min).max(1e-15);
    let n_x_bins = (n.isqrt()).max(2);
    let x_bin_w = x_range / n_x_bins as f64;

    let mut delta = 0.0_f64;
    for xi_bin in 0..n_x_bins {
        let x_lo = x_min + xi_bin as f64 * x_bin_w;
        let x_hi = x_lo + x_bin_w;
        let sub_y: Vec<f64> = x_dim
            .iter()
            .zip(y_samples.iter())
            .filter(|&(&xi, _)| xi >= x_lo && xi < x_hi)
            .map(|(_, &yi)| yi)
            .collect();
        if sub_y.is_empty() {
            continue;
        }
        let p_x = sub_y.len() as f64 / n as f64;
        let mut f_y_given_x = vec![0.0_f64; n_bins];
        for &yi in &sub_y {
            let bin = ((yi - y_min) / bin_w).floor() as usize;
            f_y_given_x[bin.min(n_bins - 1)] += 1.0 / (sub_y.len() as f64 * bin_w);
        }
        let tv: f64 = f_y
            .iter()
            .zip(f_y_given_x.iter())
            .map(|(a, b)| (a - b).abs() * bin_w)
            .sum();
        delta += p_x * tv;
    }
    delta * 0.5
}

// ─────────────────────────────────────────────────────────────────────────────
// 12. Coefficient of Variation and standard UQ metrics
// ─────────────────────────────────────────────────────────────────────────────

/// Coefficient of variation CoV = σ / |μ|.
///
/// Returns `f64::INFINITY` when the mean is zero.
pub fn coefficient_of_variation(data: &[f64]) -> f64 {
    let m = mean(data);
    if m.abs() < 1e-15 {
        return f64::INFINITY;
    }
    std_dev(data) / m.abs()
}

/// Skewness of a dataset = E\[(X−μ)³\] / σ³.
pub fn skewness(data: &[f64]) -> f64 {
    if data.len() < 3 {
        return 0.0;
    }
    let m = mean(data);
    let s = std_dev(data);
    if s < 1e-15 {
        return 0.0;
    }
    let n = data.len() as f64;
    data.iter().map(|x| ((x - m) / s).powi(3)).sum::<f64>() / n
}

/// Excess kurtosis = E\[(X−μ)⁴\] / σ⁴ − 3.
pub fn excess_kurtosis(data: &[f64]) -> f64 {
    if data.len() < 4 {
        return 0.0;
    }
    let m = mean(data);
    let s = std_dev(data);
    if s < 1e-15 {
        return 0.0;
    }
    let n = data.len() as f64;
    data.iter().map(|x| ((x - m) / s).powi(4)).sum::<f64>() / n - 3.0
}

/// Compute a P-box (probability box) from a set of CDFs.
///
/// Given a collection of sample sets (one per scenario), returns the lower
/// and upper CDF bounds at a vector of evaluation points.
///
/// # Arguments
/// * `scenarios` - Each element is a sample set for one scenario.
/// * `eval_points` - Points at which to evaluate the CDF bounds.
///
/// Returns `(cdf_lower, cdf_upper)` each of length `eval_points.len()`.
pub fn probability_box(scenarios: &[Vec<f64>], eval_points: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n_pts = eval_points.len();
    let mut lower = vec![f64::INFINITY; n_pts];
    let mut upper = vec![f64::NEG_INFINITY; n_pts];

    for scenario in scenarios {
        let n = scenario.len() as f64;
        for (k, &ep) in eval_points.iter().enumerate() {
            let cdf = scenario.iter().filter(|&&x| x <= ep).count() as f64 / n;
            if cdf < lower[k] {
                lower[k] = cdf;
            }
            if cdf > upper[k] {
                upper[k] = cdf;
            }
        }
    }
    (lower, upper)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── helper ───────────────────────────────────────────────────────────────

    fn linear_model(x: &[f64]) -> f64 {
        x.iter()
            .enumerate()
            .map(|(i, xi)| (i + 1) as f64 * xi)
            .sum()
    }

    // ── mean / variance / std_dev ─────────────────────────────────────────────

    #[test]
    fn test_mean_empty() {
        assert_eq!(mean(&[]), 0.0);
    }

    #[test]
    fn test_mean_single() {
        assert!((mean(&[42.0]) - 42.0).abs() < 1e-12);
    }

    #[test]
    fn test_mean_basic() {
        assert!((mean(&[1.0, 2.0, 3.0, 4.0, 5.0]) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_variance_zero() {
        assert_eq!(variance(&[5.0, 5.0, 5.0]), 0.0);
    }

    #[test]
    fn test_variance_known() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let v = variance(&data);
        assert!((v - 4.0).abs() < 0.1, "variance={v}");
    }

    #[test]
    fn test_std_dev_nonneg() {
        let data: Vec<f64> = (1..=100).map(|i| i as f64).collect();
        assert!(std_dev(&data) > 0.0);
    }

    // ── probit / normal_cdf ───────────────────────────────────────────────────

    #[test]
    fn test_probit_symmetry() {
        assert!((probit(0.5)).abs() < 1e-3);
    }

    #[test]
    fn test_normal_cdf_symmetry() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_probit_normal_cdf_inverse() {
        for p in [0.1, 0.25, 0.5, 0.75, 0.9] {
            let x = probit(p);
            let p2 = normal_cdf(x);
            assert!((p2 - p).abs() < 1e-3, "p={p} p2={p2}");
        }
    }

    // ── percentile ────────────────────────────────────────────────────────────

    #[test]
    fn test_percentile_median() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((percentile(&mut data, 50.0) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_percentile_min_max() {
        let mut data = vec![10.0, 20.0, 30.0];
        assert!((percentile(&mut data, 0.0) - 10.0).abs() < 1e-10);
        let mut data2 = vec![10.0, 20.0, 30.0];
        assert!((percentile(&mut data2, 100.0) - 30.0).abs() < 1e-10);
    }

    // ── hermite_poly ─────────────────────────────────────────────────────────

    #[test]
    fn test_hermite_poly_degree0() {
        assert_eq!(hermite_poly(0, 3.5), 1.0);
    }

    #[test]
    fn test_hermite_poly_degree1() {
        assert!((hermite_poly(1, 2.7) - 2.7).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_poly_degree2() {
        // He_2(x) = x² - 1
        let x = 2.0;
        assert!((hermite_poly(2, x) - (x * x - 1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_poly_degree3() {
        // He_3(x) = x³ - 3x
        let x = 1.5;
        assert!((hermite_poly(3, x) - (x.powi(3) - 3.0 * x)).abs() < 1e-12);
    }

    // ── LHS ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_lhs_shape() {
        let lhs = latin_hypercube_sample(10, 3, 42);
        assert_eq!(lhs.n_samples, 10);
        assert_eq!(lhs.n_dims, 3);
        assert_eq!(lhs.data.len(), 30);
    }

    #[test]
    fn test_lhs_bounds() {
        let lhs = latin_hypercube_sample(50, 4, 7);
        for v in &lhs.data {
            assert!(*v >= 0.0 && *v <= 1.0, "value out of [0,1]: {v}");
        }
    }

    #[test]
    fn test_lhs_stratification_1d() {
        let n = 20;
        let lhs = latin_hypercube_sample(n, 1, 123);
        let mut bins = vec![false; n];
        let step = 1.0 / n as f64;
        for i in 0..n {
            let v = lhs.get(i, 0);
            let b = (v / step).floor() as usize;
            let b = b.min(n - 1);
            bins[b] = true;
        }
        assert!(bins.iter().all(|&b| b), "LHS should cover all strata");
    }

    #[test]
    fn test_scale_lhs() {
        let lhs = latin_hypercube_sample(5, 2, 0);
        let bounds = [(0.0, 10.0), (-5.0, 5.0)];
        let scaled = scale_lhs(&lhs, &bounds);
        for row in &scaled {
            assert!(row[0] >= 0.0 && row[0] <= 10.0);
            assert!(row[1] >= -5.0 && row[1] <= 5.0);
        }
    }

    // ── Monte Carlo propagation ───────────────────────────────────────────────

    #[test]
    fn test_mc_propagation_mean() {
        // Y = X1 + X2, X1 ~ N(2,0.1), X2 ~ N(3,0.1) => E[Y] ≈ 5
        let params = vec![
            UncertainParameter::new(2.0, 0.1),
            UncertainParameter::new(3.0, 0.1),
        ];
        let result = monte_carlo_propagation(&params, 5000, &|x| x[0] + x[1], 1);
        assert!((result.mean - 5.0).abs() < 0.1, "mean={}", result.mean);
    }

    #[test]
    fn test_mc_propagation_std() {
        // Y = X, X ~ N(0, 1) => σ_Y ≈ 1
        let params = vec![UncertainParameter::new(0.0, 1.0)];
        let result = monte_carlo_propagation(&params, 10_000, &|x| x[0], 99);
        assert!((result.std - 1.0).abs() < 0.1, "std={}", result.std);
    }

    #[test]
    fn test_mc_propagation_percentiles() {
        let params = vec![UncertainParameter::new(0.0, 1.0)];
        let result = monte_carlo_propagation(&params, 10_000, &|x| x[0], 7);
        // p05 should be ≈ -1.645, p95 ≈ +1.645
        assert!(result.p05 < 0.0 && result.p95 > 0.0);
    }

    // ── Morris sensitivity ────────────────────────────────────────────────────

    #[test]
    fn test_morris_sensitivity_shape() {
        let bounds = vec![(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)];
        let result = morris_sensitivity(3, 20, 0.1, &linear_model, &bounds, 5);
        assert_eq!(result.mu_star.len(), 3);
        assert_eq!(result.sigma.len(), 3);
    }

    #[test]
    fn test_morris_sensitivity_rank_order() {
        // linear_model = 1*x0 + 2*x1 + 3*x2 — x2 has largest coefficient, mu_star[2] >= mu_star[1] >= mu_star[0]
        let bounds = vec![(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)];
        let result = morris_sensitivity(3, 50, 0.1, &linear_model, &bounds, 13);
        assert!(
            result.mu_star[2] >= result.mu_star[1],
            "mu*[2]={} mu*[1]={}",
            result.mu_star[2],
            result.mu_star[1]
        );
        assert!(
            result.mu_star[1] >= result.mu_star[0],
            "mu*[1]={} mu*[0]={}",
            result.mu_star[1],
            result.mu_star[0]
        );
    }

    // ── Sobol indices ─────────────────────────────────────────────────────────

    #[test]
    fn test_sobol_shape() {
        let bounds = vec![(0.0, 1.0), (0.0, 1.0)];
        let indices = sobol_indices(2, 100, &|x| x[0] + x[1], &bounds, 42);
        assert_eq!(indices.first_order.len(), 2);
        assert_eq!(indices.total_order.len(), 2);
    }

    #[test]
    fn test_sobol_additive_model() {
        // For Y = X1 + X2 (independent uniform), S1 ≈ S2 ≈ 0.5
        let bounds = vec![(0.0, 1.0), (0.0, 1.0)];
        let indices = sobol_indices(2, 500, &|x| x[0] + x[1], &bounds, 99);
        assert!(
            (indices.first_order[0] - indices.first_order[1]).abs() < 0.2,
            "S1={} S2={}",
            indices.first_order[0],
            indices.first_order[1]
        );
    }

    // ── Stochastic collocation ────────────────────────────────────────────────

    #[test]
    fn test_collocation_constant_model() {
        let (m, v) = stochastic_collocation(1, 3, &|_| 5.0);
        assert!((m - 5.0).abs() < 1e-6, "mean={m}");
        assert!(v < 1e-10, "var={v}");
    }

    #[test]
    fn test_collocation_linear_model_mean() {
        // Y = X, X ~ N(0,1) => E[Y] = 0
        let (m, _v) = stochastic_collocation(1, 5, &|x| x[0]);
        assert!(m.abs() < 0.1, "mean={m}");
    }

    // ── GP regression ─────────────────────────────────────────────────────────

    #[test]
    fn test_gp_fit_predict_interpolation() {
        let x_train: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64]).collect();
        let y_train: Vec<f64> = x_train.iter().map(|x| x[0] * x[0]).collect();
        let gp = fit_gp(x_train.clone(), y_train.clone(), 2.0, 1.5, 0.01);
        // Predict at training points — should be close to training values
        for (i, xi) in x_train.iter().enumerate() {
            let (pred_m, _pred_v) = gp.predict(xi);
            assert!(
                (pred_m - y_train[i]).abs() < 2.0,
                "i={i} pred={pred_m} true={}",
                y_train[i]
            );
        }
    }

    #[test]
    fn test_rbf_kernel_symmetry() {
        let x = vec![1.0, 2.0];
        let y = vec![3.0, 1.0];
        let k1 = rbf_kernel(&x, &y, 1.0, 1.0);
        let k2 = rbf_kernel(&y, &x, 1.0, 1.0);
        assert!((k1 - k2).abs() < 1e-12);
    }

    #[test]
    fn test_rbf_kernel_self_similarity() {
        let x = vec![1.0, 2.0, 3.0];
        let k = rbf_kernel(&x, &x, 2.0, 1.0);
        assert!((k - 4.0).abs() < 1e-12); // sigma_f^2 = 4
    }

    // ── FORM ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_form_simple_limit_state() {
        // g(u) = β_target - u1 (limit state is u1 = β_target)
        let beta_target = 2.0;
        let result = form_analysis(1, &|u| beta_target - u[0], 50, 1e-6, 1e-5);
        assert!(
            (result.beta - beta_target).abs() < 0.1,
            "beta={}",
            result.beta
        );
    }

    #[test]
    fn test_form_pf_decreases_with_beta() {
        let r1 = form_analysis(2, &|u| 3.0 - u[0] - u[1], 50, 1e-6, 1e-5);
        let r2 = form_analysis(2, &|u| 1.0 - u[0] - u[1], 50, 1e-6, 1e-5);
        assert!(
            r1.pf_form < r2.pf_form,
            "pf1={} pf2={}",
            r1.pf_form,
            r2.pf_form
        );
    }

    // ── SORM ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_sorm_zero_curvature_equals_form() {
        let beta = 2.0;
        let pf_form = normal_cdf(-beta);
        let pf_sorm = sorm_breitung(beta, &[]);
        assert!((pf_sorm - pf_form).abs() < 1e-12);
    }

    #[test]
    fn test_sorm_positive_curvature_reduces_pf() {
        let beta = 2.0;
        let pf_form = normal_cdf(-beta);
        let pf_sorm = sorm_breitung(beta, &[0.5]);
        assert!(pf_sorm < pf_form, "pf_sorm={pf_sorm} pf_form={pf_form}");
    }

    // ── Bootstrap CI ─────────────────────────────────────────────────────────

    #[test]
    fn test_bootstrap_ci_contains_true_mean() {
        let data: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let (lo, hi) = bootstrap_ci(&data, &|s| mean(s), 1000, 0.05, 42);
        let true_mean = 0.495;
        assert!(
            lo <= true_mean && true_mean <= hi,
            "CI=[{lo},{hi}] true_mean={true_mean}"
        );
    }

    #[test]
    fn test_bootstrap_ci_ordering() {
        let data: Vec<f64> = (1..=50).map(|i| i as f64).collect();
        let (lo, hi) = bootstrap_ci(&data, &|s| mean(s), 500, 0.1, 7);
        assert!(lo <= hi, "lo={lo} hi={hi}");
    }

    // ── Metropolis-Hastings ────────────────────────────────────────────────────

    #[test]
    fn test_mh_samples_count() {
        let config = MhConfig {
            n_samples: 200,
            burn_in: 50,
            step_size: 0.5,
            seed: 1,
        };
        let log_post = |theta: &[f64]| -0.5 * theta[0] * theta[0];
        let result = metropolis_hastings(&log_post, &[0.0], &config);
        assert_eq!(result.samples.len(), 200);
    }

    #[test]
    fn test_mh_normal_target_mean() {
        // Target: N(3, 1), log p ∝ -0.5(x-3)^2
        let config = MhConfig {
            n_samples: 5000,
            burn_in: 500,
            step_size: 1.0,
            seed: 42,
        };
        let log_post = |theta: &[f64]| -0.5 * (theta[0] - 3.0).powi(2);
        let result = metropolis_hastings(&log_post, &[0.0], &config);
        let vals: Vec<f64> = result.samples.iter().map(|s| s[0]).collect();
        let m = mean(&vals);
        assert!((m - 3.0).abs() < 0.3, "mean={m}");
    }

    #[test]
    fn test_mh_acceptance_rate_reasonable() {
        let config = MhConfig {
            n_samples: 1000,
            burn_in: 100,
            step_size: 0.5,
            seed: 77,
        };
        let log_post = |theta: &[f64]| -0.5 * theta[0] * theta[0];
        let result = metropolis_hastings(&log_post, &[0.0], &config);
        // Acceptance rate should be between 5% and 95%
        assert!(
            result.acceptance_rate > 0.05 && result.acceptance_rate < 0.95,
            "acceptance_rate={}",
            result.acceptance_rate
        );
    }

    // ── Probability of failure ─────────────────────────────────────────────────

    #[test]
    fn test_pof_mc_certain_failure() {
        // g(x) = -1 always — should always fail
        let params = vec![UncertainParameter::new(0.0, 1.0)];
        let (pf, _cov) = probability_of_failure_mc(&params, &|_| -1.0, 1000, 0);
        assert!((pf - 1.0).abs() < 1e-10, "pf={pf}");
    }

    #[test]
    fn test_pof_mc_certain_safe() {
        // g(x) = 1 always — should never fail
        let params = vec![UncertainParameter::new(0.0, 1.0)];
        let (pf, _cov) = probability_of_failure_mc(&params, &|_| 1.0, 1000, 0);
        assert_eq!(pf, 0.0);
    }

    #[test]
    fn test_pof_mc_reasonable_estimate() {
        // g(x) = 2 - x, X ~ N(0,1): failure when g <= 0 means x >= 2 => P_f = 1 - Φ(2) ≈ 0.0228
        let params = vec![UncertainParameter::new(0.0, 1.0)];
        let (pf, _cov) = probability_of_failure_mc(&params, &|x| 2.0 - x[0], 50_000, 5);
        assert!((pf - 0.0228).abs() < 0.01, "pf={pf}, expected≈0.0228");
    }

    // ── Statistics helpers ─────────────────────────────────────────────────────

    #[test]
    fn test_cov_finite() {
        let data = vec![1.0, 2.0, 3.0];
        let cov = coefficient_of_variation(&data);
        assert!(cov.is_finite() && cov > 0.0);
    }

    #[test]
    fn test_skewness_symmetric() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(skewness(&data).abs() < 0.1);
    }

    #[test]
    fn test_excess_kurtosis_normal() {
        // Generate many N(0,1) samples via the internal RNG
        let mut rng = UqRng::new(123);
        let data: Vec<f64> = (0..5000).map(|_| rng.next_normal()).collect();
        let k = excess_kurtosis(&data);
        assert!(k.abs() < 0.5, "excess_kurtosis={k}");
    }

    // ── P-box ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_probability_box_bounds() {
        let s1 = vec![1.0, 2.0, 3.0, 4.0];
        let s2 = vec![0.5, 1.5, 3.5, 5.0];
        let (lower, upper) = probability_box(&[s1, s2], &[2.0, 3.0]);
        assert!(lower[0] <= upper[0]);
        assert!(lower[1] <= upper[1]);
    }

    // ── Reliability index ──────────────────────────────────────────────────────

    #[test]
    fn test_reliability_index_known() {
        // P_f = Φ(-β) => β = -Φ⁻¹(P_f)
        let pf = normal_cdf(-2.0);
        let beta = reliability_index(pf);
        assert!((beta - 2.0).abs() < 0.05, "beta={beta}");
    }

    // ── PCE mean and variance ──────────────────────────────────────────────────

    #[test]
    fn test_pce_mean_constant_model() {
        let pce = fit_pce(1, 1, 10, &|_| 7.0, 42);
        assert!(
            (pce.pce_mean() - 7.0).abs() < 0.5,
            "mean={}",
            pce.pce_mean()
        );
    }

    #[test]
    fn test_pce_evaluate_consistency() {
        let pce = fit_pce(2, 2, 30, &|x| x[0] + x[1], 3);
        // E[X1 + X2] = 0 for standard normal
        let pred = pce.evaluate(&[0.0, 0.0]);
        assert!(pred.abs() < 2.0, "pred={pred}");
    }

    // ── Gauss-Hermite nodes ────────────────────────────────────────────────────

    #[test]
    fn test_gauss_hermite_order1_weight() {
        let nodes = gauss_hermite_nodes(1);
        assert_eq!(nodes.len(), 1);
        assert!((nodes[0].weight - PI.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn test_gauss_hermite_order3_symmetry() {
        let nodes = gauss_hermite_nodes(3);
        assert_eq!(nodes.len(), 3);
        // Nodes should be symmetric: nodes[0].point ≈ -nodes[2].point
        assert!((nodes[0].point + nodes[2].point).abs() < 1e-8);
    }
}
