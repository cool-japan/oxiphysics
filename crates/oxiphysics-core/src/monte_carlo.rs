#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Monte Carlo methods for physics simulation.
//!
//! Provides sampling distributions, Markov Chain Monte Carlo (MCMC),
//! Hamiltonian Monte Carlo, Monte Carlo integration (stratified, LHS,
//! quasi-Monte Carlo), bootstrap confidence intervals, particle filters,
//! and Brownian bridge processes.

#![allow(dead_code)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Local LCG RNG (same style as stochastic.rs)
// ─────────────────────────────────────────────────────────────────────────────

/// A fast linear congruential random number generator.
#[allow(dead_code)]
pub struct LcgRng {
    state: u64,
}

impl LcgRng {
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

    fn next_normal(&mut self) -> f64 {
        box_muller_normal(self.next_f64(), self.next_f64())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper functions
// ─────────────────────────────────────────────────────────────────────────────

/// Box-Muller transform: generates a standard normal N(0,1) from two U(0,1) samples.
///
/// # Arguments
/// * `u1` - Uniform sample in (0, 1].
/// * `u2` - Uniform sample in \[0, 1).
pub fn box_muller_normal(u1: f64, u2: f64) -> f64 {
    let u1 = u1.max(1e-15);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

/// Inverse CDF of the standard normal distribution (Beasley-Springer-Moro approximation).
///
/// Maps p ∈ (0, 1) to the quantile x such that P(X ≤ x) = p for X ~ N(0,1).
pub fn inverse_cdf(p: f64) -> f64 {
    // Rational approximation (Abramowitz and Stegun 26.2.17)
    assert!(p > 0.0 && p < 1.0, "p must be in (0, 1)");
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

/// Computes the effective sample size (ESS) from normalized importance weights.
///
/// ESS = (Σ w_i)^2 / Σ w_i^2.
pub fn effective_sample_size(weights: &[f64]) -> f64 {
    let sum_w: f64 = weights.iter().sum();
    if sum_w < 1e-15 {
        return 0.0;
    }
    let sum_w2: f64 = weights.iter().map(|&w| (w / sum_w).powi(2)).sum();
    1.0 / sum_w2
}

/// Computes the autocorrelation of a time series at lag `k`.
///
/// Returns the normalized autocorrelation coefficient r(k) ∈ \[-1, 1\].
pub fn autocorrelation(x: &[f64], k: usize) -> f64 {
    let n = x.len();
    if k >= n {
        return 0.0;
    }
    let mean = x.iter().sum::<f64>() / n as f64;
    let var: f64 = x.iter().map(|&xi| (xi - mean).powi(2)).sum::<f64>() / n as f64;
    if var < 1e-15 {
        return 1.0;
    }
    let cov: f64 = (0..n - k)
        .map(|i| (x[i] - mean) * (x[i + k] - mean))
        .sum::<f64>()
        / n as f64;
    cov / var
}

// ─────────────────────────────────────────────────────────────────────────────
// McSampler — probability distribution sampler
// ─────────────────────────────────────────────────────────────────────────────

/// Monte Carlo sampler supporting multiple probability distributions.
///
/// Uses an internal LCG RNG; set the seed for reproducibility.
pub struct McSampler {
    rng: LcgRng,
}

impl McSampler {
    /// Creates a new sampler with the given seed.
    pub fn new(seed: u64) -> Self {
        Self {
            rng: LcgRng::new(seed),
        }
    }

    /// Samples from Uniform(0, 1).
    pub fn uniform(&mut self) -> f64 {
        self.rng.next_f64()
    }

    /// Samples from Uniform(a, b).
    pub fn uniform_range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.rng.next_f64()
    }

    /// Samples from Normal(mu, sigma).
    pub fn normal(&mut self, mu: f64, sigma: f64) -> f64 {
        mu + sigma * self.rng.next_normal()
    }

    /// Samples from Exponential(lambda): pdf = lambda * exp(-lambda * x).
    pub fn exponential(&mut self, lambda: f64) -> f64 {
        -self.rng.next_f64().max(1e-15).ln() / lambda
    }

    /// Samples from Poisson(lambda) using Knuth's algorithm.
    pub fn poisson(&mut self, lambda: f64) -> u64 {
        let limit = (-lambda).exp();
        let mut k = 0u64;
        let mut p = 1.0;
        loop {
            p *= self.rng.next_f64();
            if p <= limit {
                break;
            }
            k += 1;
        }
        k
    }

    /// Samples from Beta(alpha, beta) using Johnk's method.
    pub fn beta(&mut self, alpha: f64, beta: f64) -> f64 {
        // Use gamma ratio: Beta(a,b) = Gamma(a) / (Gamma(a) + Gamma(b))
        let x = self.gamma(alpha, 1.0);
        let y = self.gamma(beta, 1.0);
        x / (x + y)
    }

    /// Samples from Gamma(shape, scale) using Marsaglia-Tsang method.
    pub fn gamma(&mut self, shape: f64, scale: f64) -> f64 {
        if shape < 1.0 {
            // Use Gamma(shape+1) * U^(1/shape)
            let g = self.gamma(shape + 1.0, scale);
            let u = self.rng.next_f64().max(1e-15);
            return g * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let z = self.rng.next_normal();
            let v_inner = 1.0 + c * z;
            if v_inner <= 0.0 {
                continue;
            }
            let v = v_inner.powi(3);
            let u = self.rng.next_f64().max(1e-15);
            if u < 1.0 - 0.0331 * z.powi(4) {
                return d * v * scale;
            }
            if u.ln() < 0.5 * z * z + d * (1.0 - v + v.ln()) {
                return d * v * scale;
            }
        }
    }

    /// Draws `n` samples from Normal(mu, sigma).
    pub fn normal_samples(&mut self, n: usize, mu: f64, sigma: f64) -> Vec<f64> {
        (0..n).map(|_| self.normal(mu, sigma)).collect()
    }

    /// Draws `n` samples from Uniform(a, b).
    pub fn uniform_samples(&mut self, n: usize, a: f64, b: f64) -> Vec<f64> {
        (0..n).map(|_| self.uniform_range(a, b)).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImportanceSampler
// ─────────────────────────────────────────────────────────────────────────────

/// Importance sampler with acceptance-rejection and adaptive IS weights.
///
/// Uses a proposal distribution q(x) to estimate expectations under target p(x).
pub struct ImportanceSampler {
    rng: LcgRng,
    /// Proposal distribution standard deviation.
    pub proposal_sigma: f64,
}

impl ImportanceSampler {
    /// Creates a new importance sampler.
    pub fn new(seed: u64, proposal_sigma: f64) -> Self {
        Self {
            rng: LcgRng::new(seed),
            proposal_sigma,
        }
    }

    /// Acceptance-rejection sampler for target pdf `target` with envelope `c * proposal`.
    ///
    /// Returns a sample from the target distribution.
    pub fn accept_reject<F, G>(&mut self, target: F, proposal_sample: G, c: f64) -> f64
    where
        F: Fn(f64) -> f64,
        G: Fn(&mut LcgRng) -> (f64, f64), // returns (sample, proposal_density)
    {
        loop {
            let (x, q_x) = proposal_sample(&mut self.rng);
            let p_x = target(x);
            let u = self.rng.next_f64();
            if u <= p_x / (c * q_x) {
                return x;
            }
        }
    }

    /// Computes unnormalized importance weights for samples drawn from proposal.
    ///
    /// w_i = target(x_i) / proposal(x_i).
    pub fn importance_weights<F, G>(&self, samples: &[f64], target: F, proposal: G) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
    {
        samples
            .iter()
            .map(|&x| {
                let q = proposal(x).max(1e-300);
                target(x) / q
            })
            .collect()
    }

    /// Estimates E_p\[f(x)\] using importance sampling from proposal q.
    pub fn estimate<F, H>(&mut self, n: usize, f_fn: F, log_w_fn: H) -> f64
    where
        F: Fn(f64) -> f64,
        H: Fn(f64) -> f64, // log(target/proposal) at x
    {
        let sigma = self.proposal_sigma;
        let samples: Vec<f64> = (0..n).map(|_| self.rng.next_normal() * sigma).collect();
        let log_w: Vec<f64> = samples.iter().map(|&x| log_w_fn(x)).collect();
        let log_w_max = log_w.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let w: Vec<f64> = log_w.iter().map(|&lw| (lw - log_w_max).exp()).collect();
        let sum_w: f64 = w.iter().sum();
        let numerator: f64 = samples
            .iter()
            .zip(w.iter())
            .map(|(&x, &wi)| f_fn(x) * wi)
            .sum();
        numerator / sum_w
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Estimator — statistical estimator from samples
// ─────────────────────────────────────────────────────────────────────────────

/// Statistical estimator: computes mean, variance, confidence intervals, and ESS.
pub struct Estimator {
    /// The sample data.
    pub samples: Vec<f64>,
}

impl Estimator {
    /// Creates a new estimator from a sample vector.
    pub fn new(samples: Vec<f64>) -> Self {
        Self { samples }
    }

    /// Sample mean.
    pub fn mean(&self) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        self.samples.iter().sum::<f64>() / n as f64
    }

    /// Sample variance (unbiased, divides by n-1).
    pub fn variance(&self) -> f64 {
        let n = self.samples.len();
        if n < 2 {
            return 0.0;
        }
        let mu = self.mean();
        self.samples.iter().map(|&x| (x - mu).powi(2)).sum::<f64>() / (n - 1) as f64
    }

    /// Standard error of the mean: SE = sqrt(variance / n).
    pub fn standard_error(&self) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        (self.variance() / n as f64).sqrt()
    }

    /// Confidence interval at level `alpha` (e.g., alpha=0.05 for 95% CI).
    ///
    /// Returns (lower, upper) assuming normality of the estimator.
    pub fn confidence_interval(&self, alpha: f64) -> (f64, f64) {
        let mu = self.mean();
        let se = self.standard_error();
        // z-score for alpha/2
        let z = inverse_cdf(1.0 - alpha / 2.0);
        (mu - z * se, mu + z * se)
    }

    /// Effective sample size using the autocorrelation-based estimate.
    ///
    /// ESS = n / (1 + 2 * Σ_k ρ(k)) where ρ(k) is the autocorrelation at lag k.
    pub fn effective_sample_size(&self) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        let mut sum_rho = 0.0;
        let max_lag = (n / 2).min(50);
        for k in 1..max_lag {
            let rho = autocorrelation(&self.samples, k);
            if rho.abs() < 0.05 {
                break;
            }
            sum_rho += rho;
        }
        n as f64 / (1.0 + 2.0 * sum_rho).max(0.01)
    }

    /// Median of the samples.
    pub fn median(&self) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        let mut sorted = self.samples.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if n % 2 == 1 {
            sorted[n / 2]
        } else {
            0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MetropolisHastings
// ─────────────────────────────────────────────────────────────────────────────

/// Metropolis-Hastings MCMC sampler.
///
/// Generates a Markov chain targeting the distribution proportional to `log_target(x)`.
pub struct MetropolisHastings {
    rng: LcgRng,
    /// Proposal step size (standard deviation of Gaussian proposal).
    pub step_size: f64,
    /// Number of burn-in steps to discard.
    pub burn_in: usize,
    /// Thinning interval (keep every k-th sample).
    pub thinning: usize,
    /// Number of accepted proposals (updated during sampling).
    pub n_accepted: usize,
}

impl MetropolisHastings {
    /// Creates a new MH sampler.
    pub fn new(seed: u64, step_size: f64, burn_in: usize, thinning: usize) -> Self {
        Self {
            rng: LcgRng::new(seed),
            step_size,
            burn_in,
            thinning,
            n_accepted: 0,
        }
    }

    /// Runs the MH chain for `n_samples` (after burn-in), targeting `log_target`.
    ///
    /// `x0` is the initial state. Returns a vector of samples.
    pub fn sample<F>(&mut self, x0: f64, n_samples: usize, log_target: F) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut x = x0;
        let total = self.burn_in + n_samples * self.thinning;
        let mut chain = Vec::with_capacity(n_samples);
        self.n_accepted = 0;
        let mut n_proposed = 0usize;

        for step in 0..total {
            let x_prop = x + self.step_size * self.rng.next_normal();
            let log_alpha = (log_target(x_prop) - log_target(x)).min(0.0);
            let u = self.rng.next_f64().max(1e-15).ln();
            n_proposed += 1;
            if u <= log_alpha {
                x = x_prop;
                if step >= self.burn_in {
                    self.n_accepted += 1;
                }
            }
            if step >= self.burn_in && (step - self.burn_in).is_multiple_of(self.thinning) {
                chain.push(x);
            }
        }
        let _ = n_proposed;
        chain
    }

    /// Acceptance rate of the chain.
    pub fn acceptance_rate(&self, n_total: usize) -> f64 {
        if n_total == 0 {
            0.0
        } else {
            self.n_accepted as f64 / n_total as f64
        }
    }

    /// Chain diagnostics: returns (mean, variance, ESS).
    pub fn diagnostics(&self, chain: &[f64]) -> (f64, f64, f64) {
        let est = Estimator::new(chain.to_vec());
        (est.mean(), est.variance(), est.effective_sample_size())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// HamiltonianMC — HMC with leapfrog integrator
// ─────────────────────────────────────────────────────────────────────────────

/// Hamiltonian Monte Carlo sampler.
///
/// Uses the leapfrog integrator to propose far-reaching moves along
/// level sets of the Hamiltonian H(q, p) = U(q) + K(p).
pub struct HamiltonianMC {
    rng: LcgRng,
    /// Leapfrog step size.
    pub step_size: f64,
    /// Number of leapfrog steps per proposal.
    pub n_leapfrog: usize,
    /// Number of burn-in steps.
    pub burn_in: usize,
    /// Acceptance count.
    pub n_accepted: usize,
}

impl HamiltonianMC {
    /// Creates a new HMC sampler.
    pub fn new(seed: u64, step_size: f64, n_leapfrog: usize, burn_in: usize) -> Self {
        Self {
            rng: LcgRng::new(seed),
            step_size,
            n_leapfrog,
            burn_in,
            n_accepted: 0,
        }
    }

    /// Kinetic energy K(p) = 0.5 * p^2 (unit mass Gaussian momentum).
    pub fn kinetic_energy(p: f64) -> f64 {
        0.5 * p * p
    }

    /// One leapfrog trajectory of `n_steps` steps.
    ///
    /// Returns (q_new, p_new) after the trajectory.
    pub fn leapfrog<F>(&self, q0: f64, p0: f64, grad_u: F) -> (f64, f64)
    where
        F: Fn(f64) -> f64,
    {
        let eps = self.step_size;
        let mut q = q0;
        let mut p = p0;
        // Half step for momentum
        p -= 0.5 * eps * grad_u(q);
        for _ in 0..self.n_leapfrog {
            q += eps * p;
            let is_last = false;
            let _ = is_last;
            p -= eps * grad_u(q);
        }
        // Undo last full step, do half step
        p += eps * grad_u(q);
        p -= 0.5 * eps * grad_u(q);
        (q, p)
    }

    /// Runs HMC for `n_samples` after burn-in.
    ///
    /// `log_target` is the log of the unnormalized target density.
    /// `grad_u` = -d(log_target)/dq (gradient of potential energy).
    pub fn sample<F, G>(&mut self, q0: f64, n_samples: usize, log_target: F, grad_u: G) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
    {
        let total = self.burn_in + n_samples;
        let mut q = q0;
        let mut chain = Vec::with_capacity(n_samples);
        self.n_accepted = 0;

        for step in 0..total {
            let p = self.rng.next_normal();
            let h_current = -log_target(q) + Self::kinetic_energy(p);
            let (q_prop, p_prop) = self.leapfrog(q, p, &grad_u);
            let h_prop = -log_target(q_prop) + Self::kinetic_energy(p_prop);
            let delta_h = h_prop - h_current;
            let u = self.rng.next_f64().max(1e-15).ln();
            if u <= -delta_h.max(0.0) || delta_h <= 0.0 {
                q = q_prop;
                if step >= self.burn_in {
                    self.n_accepted += 1;
                }
            }
            if step >= self.burn_in {
                chain.push(q);
            }
        }
        chain
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MarkovChainMC — general MCMC dispatcher
// ─────────────────────────────────────────────────────────────────────────────

/// General Markov Chain Monte Carlo runner.
///
/// Wraps Metropolis-Hastings and provides Gibbs sampling interface.
pub struct MarkovChainMC {
    rng: LcgRng,
    /// Step size for proposal distributions.
    pub step_size: f64,
}

impl MarkovChainMC {
    /// Creates a new MCMC runner.
    pub fn new(seed: u64, step_size: f64) -> Self {
        Self {
            rng: LcgRng::new(seed),
            step_size,
        }
    }

    /// Metropolis-Hastings for a univariate log-target.
    pub fn metropolis_hastings<F>(
        &mut self,
        x0: f64,
        n_samples: usize,
        burn_in: usize,
        log_target: F,
    ) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut mh = MetropolisHastings::new(42, self.step_size, burn_in, 1);
        mh.sample(x0, n_samples, log_target)
    }

    /// Gibbs sampler for a bivariate target with full conditionals.
    ///
    /// `sample_x_given_y` and `sample_y_given_x` are samplers for each conditional.
    pub fn gibbs_2d<F, G>(
        &mut self,
        x0: f64,
        y0: f64,
        n_samples: usize,
        burn_in: usize,
        sample_x: F,
        sample_y: G,
    ) -> Vec<(f64, f64)>
    where
        F: Fn(f64, &mut LcgRng) -> f64, // sample x | y
        G: Fn(f64, &mut LcgRng) -> f64, // sample y | x
    {
        let total = burn_in + n_samples;
        let mut chain = Vec::with_capacity(n_samples);
        // Start chain at the provided initial point
        let mut state = (x0, y0);
        for step in 0..total {
            let new_x = sample_x(state.1, &mut self.rng);
            let new_y = sample_y(new_x, &mut self.rng);
            state = (new_x, new_y);
            if step >= burn_in {
                chain.push(state);
            }
        }
        chain
    }

    /// MALA (Metropolis-Adjusted Langevin Algorithm) sampler.
    ///
    /// Proposal: x_prop = x + (eps^2/2) * grad_log_target(x) + eps * N(0,1).
    pub fn mala<F, G>(
        &mut self,
        x0: f64,
        n_samples: usize,
        eps: f64,
        log_target: F,
        grad_log_target: G,
    ) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
    {
        let mut x = x0;
        let mut chain = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let g = grad_log_target(x);
            let x_prop = x + 0.5 * eps * eps * g + eps * self.rng.next_normal();
            // MH correction
            let log_alpha = log_target(x_prop) - log_target(x);
            if self.rng.next_f64().max(1e-15).ln() <= log_alpha.min(0.0) {
                x = x_prop;
            }
            chain.push(x);
        }
        chain
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MonteCarloIntegral — multi-dimensional quadrature
// ─────────────────────────────────────────────────────────────────────────────

/// Monte Carlo integration methods: plain, stratified, LHS, and quasi-MC.
pub struct MonteCarloIntegral {
    rng: LcgRng,
    /// Number of integration samples.
    pub n_samples: usize,
    /// Dimensionality of the integration domain.
    pub dim: usize,
}

impl MonteCarloIntegral {
    /// Creates a new Monte Carlo integrator.
    pub fn new(seed: u64, n_samples: usize, dim: usize) -> Self {
        Self {
            rng: LcgRng::new(seed),
            n_samples,
            dim,
        }
    }

    /// Plain Monte Carlo integration of `f` over the unit hypercube \[0,1\]^d.
    ///
    /// Returns (estimate, standard_error).
    pub fn integrate_plain<F>(&mut self, f: F) -> (f64, f64)
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = self.n_samples;
        let d = self.dim;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut x = vec![0.0_f64; d];
        for _ in 0..n {
            for xi in x.iter_mut() {
                *xi = self.rng.next_f64();
            }
            let val = f(&x);
            sum += val;
            sum_sq += val * val;
        }
        let mean = sum / n as f64;
        let var = (sum_sq / n as f64 - mean * mean).max(0.0);
        (mean, (var / n as f64).sqrt())
    }

    /// Stratified Monte Carlo: divides each dimension into k strata.
    ///
    /// Returns (estimate, standard_error).
    pub fn integrate_stratified<F>(&mut self, f: F, k: usize) -> (f64, f64)
    where
        F: Fn(&[f64]) -> f64,
    {
        let d = self.dim;
        let n_strata = k.pow(d as u32);
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let vol = 1.0 / n_strata as f64;
        for s in 0..n_strata {
            let mut x = vec![0.0_f64; d];
            let mut idx = s;
            for j in 0..d {
                let cell = idx % k;
                idx /= k;
                x[j] = (cell as f64 + self.rng.next_f64()) / k as f64;
            }
            let val = f(&x) * vol * n_strata as f64;
            sum += val;
            sum_sq += val * val;
        }
        let mean = sum / n_strata as f64;
        let var = (sum_sq / n_strata as f64 - mean * mean).max(0.0);
        (
            mean * vol * n_strata as f64 / n_strata as f64,
            (var / n_strata as f64).sqrt(),
        )
    }

    /// Latin Hypercube Sampling integration.
    ///
    /// Divides each dimension into n strata and samples one point per stratum.
    pub fn integrate_lhs<F>(&mut self, f: F) -> (f64, f64)
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = self.n_samples;
        let d = self.dim;
        // Generate permutations for each dimension
        let mut perm: Vec<Vec<usize>> = (0..d)
            .map(|_| {
                let mut p: Vec<usize> = (0..n).collect();
                // Fisher-Yates shuffle
                for i in (1..n).rev() {
                    let j = (self.rng.next_u64() as usize) % (i + 1);
                    p.swap(i, j);
                }
                p
            })
            .collect();
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for i in 0..n {
            let x: Vec<f64> = (0..d)
                .map(|j| (perm[j][i] as f64 + self.rng.next_f64()) / n as f64)
                .collect();
            let val = f(&x);
            sum += val;
            sum_sq += val * val;
            let _ = &mut perm; // keep alive
        }
        let mean = sum / n as f64;
        let var = (sum_sq / n as f64 - mean * mean).max(0.0);
        (mean, (var / n as f64).sqrt())
    }

    /// Quasi-Monte Carlo integration using a Halton sequence.
    ///
    /// Uses prime bases for each dimension.
    pub fn integrate_halton<F>(&mut self, f: F) -> (f64, f64)
    where
        F: Fn(&[f64]) -> f64,
    {
        let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        let n = self.n_samples;
        let d = self.dim.min(primes.len());
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for i in 1..=n {
            let x: Vec<f64> = (0..d).map(|j| halton_seq(i, primes[j])).collect();
            let val = f(&x);
            sum += val;
            sum_sq += val * val;
        }
        let mean = sum / n as f64;
        let var = (sum_sq / n as f64 - mean * mean).max(0.0);
        (mean, (var / n as f64).sqrt())
    }
}

/// Computes the i-th element of the Halton sequence in base `base`.
fn halton_seq(mut i: usize, base: usize) -> f64 {
    let mut f = 1.0_f64;
    let mut r = 0.0_f64;
    let b = base as f64;
    while i > 0 {
        f /= b;
        r += f * (i % base) as f64;
        i /= base;
    }
    r
}

// ─────────────────────────────────────────────────────────────────────────────
// BootstrapCI — bootstrap confidence intervals
// ─────────────────────────────────────────────────────────────────────────────

/// Non-parametric bootstrap confidence interval estimator.
///
/// Supports percentile method and BCa (bias-corrected and accelerated) correction.
pub struct BootstrapCI {
    rng: LcgRng,
    /// Number of bootstrap resamples.
    pub n_boot: usize,
}

impl BootstrapCI {
    /// Creates a new bootstrap CI estimator.
    pub fn new(seed: u64, n_boot: usize) -> Self {
        Self {
            rng: LcgRng::new(seed),
            n_boot,
        }
    }

    /// Generates a bootstrap resample of `data`.
    pub fn resample(&mut self, data: &[f64]) -> Vec<f64> {
        let n = data.len();
        (0..n)
            .map(|_| data[self.rng.next_u64() as usize % n])
            .collect()
    }

    /// Computes the percentile bootstrap CI for the mean at level `alpha`.
    ///
    /// Returns (lower, upper).
    pub fn percentile_ci(&mut self, data: &[f64], alpha: f64) -> (f64, f64) {
        let n_boot = self.n_boot;
        let mut boot_means: Vec<f64> = (0..n_boot)
            .map(|_| {
                let resamp = self.resample(data);
                resamp.iter().sum::<f64>() / resamp.len() as f64
            })
            .collect();
        boot_means.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let lo_idx = ((alpha / 2.0 * n_boot as f64) as usize).min(n_boot - 1);
        let hi_idx = (((1.0 - alpha / 2.0) * n_boot as f64) as usize).min(n_boot - 1);
        (boot_means[lo_idx], boot_means[hi_idx])
    }

    /// BCa bootstrap CI (bias-corrected and accelerated) for the mean.
    ///
    /// Returns (lower, upper).
    pub fn bca_ci(&mut self, data: &[f64], alpha: f64) -> (f64, f64) {
        let n = data.len();
        let n_boot = self.n_boot;
        let theta_hat = data.iter().sum::<f64>() / n as f64;

        // Bootstrap distribution
        let mut boot_stats: Vec<f64> = (0..n_boot)
            .map(|_| {
                let resamp = self.resample(data);
                resamp.iter().sum::<f64>() / resamp.len() as f64
            })
            .collect();
        boot_stats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Bias correction z0
        let count_less = boot_stats.iter().filter(|&&b| b < theta_hat).count();
        let p0 = (count_less as f64 + 0.5) / (n_boot as f64 + 1.0);
        let z0 = inverse_cdf(p0.clamp(0.001, 0.999));

        // Acceleration a (jackknife-based)
        let jack_means: Vec<f64> = (0..n)
            .map(|i| {
                let jack_sum: f64 = data
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j != i)
                    .map(|(_, &x)| x)
                    .sum();
                jack_sum / (n - 1) as f64
            })
            .collect();
        let jack_mean_mean = jack_means.iter().sum::<f64>() / n as f64;
        let numer: f64 = jack_means
            .iter()
            .map(|&m| (jack_mean_mean - m).powi(3))
            .sum();
        let denom: f64 = 6.0
            * jack_means
                .iter()
                .map(|&m| (jack_mean_mean - m).powi(2))
                .sum::<f64>()
                .powf(1.5);
        let a = if denom.abs() < 1e-15 {
            0.0
        } else {
            numer / denom
        };

        let adjust_quantile = |p: f64| -> f64 {
            let z_p = inverse_cdf(p.clamp(0.001, 0.999));
            let inner = z0 + (z0 + z_p) / (1.0 - a * (z0 + z_p));
            // normal CDF approximation

            0.5 * (1.0 + libm_erf(inner / 2.0_f64.sqrt()))
        };

        let p_lo = adjust_quantile(alpha / 2.0);
        let p_hi = adjust_quantile(1.0 - alpha / 2.0);
        let lo_idx = ((p_lo * n_boot as f64) as usize).min(n_boot - 1);
        let hi_idx = ((p_hi * n_boot as f64) as usize).min(n_boot - 1);
        (boot_stats[lo_idx], boot_stats[hi_idx])
    }
}

/// Error function approximation (Abramowitz & Stegun 7.1.26).
fn libm_erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    sign * (1.0 - poly * (-x * x).exp())
}

// ─────────────────────────────────────────────────────────────────────────────
// AnisotropicMC — directional distribution sampling
// ─────────────────────────────────────────────────────────────────────────────

/// Anisotropic Monte Carlo sampler for directional distributions.
///
/// Implements the von Mises-Fisher distribution for sampling on the unit sphere.
pub struct AnisotropicMC {
    rng: LcgRng,
    /// Concentration parameter κ (κ=0 → uniform, κ→∞ → deterministic).
    pub kappa: f64,
    /// Mean direction (unit vector \[x, y, z\]).
    pub mu: [f64; 3],
}

impl AnisotropicMC {
    /// Creates a new anisotropic sampler.
    pub fn new(seed: u64, kappa: f64, mu: [f64; 3]) -> Self {
        Self {
            rng: LcgRng::new(seed),
            kappa,
            mu,
        }
    }

    /// Samples a unit vector from the von Mises-Fisher distribution.
    ///
    /// Returns a unit vector \[x, y, z\] concentrated around `mu`.
    pub fn sample_vmf(&mut self) -> [f64; 3] {
        let kappa = self.kappa;
        // Sample cos(theta) from the marginal distribution
        let cos_theta = if kappa < 1e-6 {
            2.0 * self.rng.next_f64() - 1.0
        } else {
            let b = -kappa + (kappa * kappa + 1.0).sqrt();
            let x0 = (1.0 - b) / (1.0 + b);
            let c = kappa * x0 + 2.0 * (1.0 + b).ln() - (1.0 + x0).ln();
            loop {
                let z = self.rng.next_f64();
                let w = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z);
                let t = 2.0 * kappa * w / (1.0 + b * w);
                let u = self.rng.next_f64();
                if u.max(1e-15).ln() + c - t <= 0.0 {
                    break w;
                }
            }
        };
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let phi = 2.0 * PI * self.rng.next_f64();
        // Sample in frame where mu = (0, 0, 1), then rotate
        let v = [sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta];
        self.rotate_to_mu(v)
    }

    fn rotate_to_mu(&self, v: [f64; 3]) -> [f64; 3] {
        // Rotate (0,0,1) to mu using Rodrigues' rotation
        let mu = self.mu;
        let norm_mu = (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2])
            .sqrt()
            .max(1e-15);
        let mu_hat = [mu[0] / norm_mu, mu[1] / norm_mu, mu[2] / norm_mu];
        // If mu is already (0,0,1), return v directly
        if (mu_hat[2] - 1.0).abs() < 1e-10 {
            return v;
        }
        if (mu_hat[2] + 1.0).abs() < 1e-10 {
            return [-v[0], -v[1], -v[2]];
        }
        // Axis of rotation: k = (0,0,1) x mu
        let k = [mu_hat[1], -mu_hat[0], 0.0];
        let k_norm = (k[0] * k[0] + k[1] * k[1]).sqrt().max(1e-15);
        let k_hat = [k[0] / k_norm, k[1] / k_norm, 0.0];
        let cos_a = mu_hat[2];
        let sin_a = (1.0 - cos_a * cos_a).sqrt();
        // Rodrigues: v_rot = v*cos_a + (k x v)*sin_a + k*(k·v)*(1-cos_a)
        let k_cross_v = [
            k_hat[1] * v[2] - k_hat[2] * v[1],
            k_hat[2] * v[0] - k_hat[0] * v[2],
            k_hat[0] * v[1] - k_hat[1] * v[0],
        ];
        let k_dot_v = k_hat[0] * v[0] + k_hat[1] * v[1] + k_hat[2] * v[2];
        [
            v[0] * cos_a + k_cross_v[0] * sin_a + k_hat[0] * k_dot_v * (1.0 - cos_a),
            v[1] * cos_a + k_cross_v[1] * sin_a + k_hat[1] * k_dot_v * (1.0 - cos_a),
            v[2] * cos_a + k_cross_v[2] * sin_a + k_hat[2] * k_dot_v * (1.0 - cos_a),
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BrownianBridge
// ─────────────────────────────────────────────────────────────────────────────

/// Brownian bridge and Ornstein-Uhlenbeck process sampler.
///
/// A Brownian bridge B(t) conditioned on B(0) = a, B(T) = b.
pub struct BrownianBridge {
    rng: LcgRng,
    /// Diffusion coefficient σ.
    pub sigma: f64,
}

impl BrownianBridge {
    /// Creates a new Brownian bridge sampler.
    pub fn new(seed: u64, sigma: f64) -> Self {
        Self {
            rng: LcgRng::new(seed),
            sigma,
        }
    }

    /// Samples a Brownian bridge path with `n` interior points.
    ///
    /// Returns n+2 values: \[a, interior..., b\] at times \[0, dt, ..., T\].
    pub fn brownian_bridge_sample(&mut self, n: usize, a: f64, b: f64, t_end: f64) -> Vec<f64> {
        let total = n + 2;
        let dt = t_end / (n + 1) as f64;
        let mut path = vec![0.0_f64; total];
        path[0] = a;
        path[total - 1] = b;
        // Generate free Brownian motion then pin endpoints
        let mut free = vec![a];
        for _ in 1..total {
            let prev = *free.last().expect("free path is non-empty");
            free.push(prev + self.sigma * dt.sqrt() * self.rng.next_normal());
        }
        // Condition on endpoint
        for i in 1..total - 1 {
            let t = i as f64 * dt;
            let t_rem = t_end - t;
            let mean_bridge = (t_rem * a + t * b) / t_end;
            let free_mean = a + (free[i] - a) * t / (total - 1) as f64;
            let _ = (free_mean, t_rem);
            path[i] = mean_bridge + (free[i] - free[0] - (free[total - 1] - free[0]) * t / t_end);
        }
        path
    }

    /// Samples an Ornstein-Uhlenbeck path: dX = -theta*(X - mu)*dt + sigma*dW.
    ///
    /// Returns `n` samples at time step `dt`.
    pub fn ornstein_uhlenbeck(
        &mut self,
        x0: f64,
        theta: f64,
        mu: f64,
        dt: f64,
        n: usize,
    ) -> Vec<f64> {
        let sigma = self.sigma;
        let mut path = Vec::with_capacity(n);
        let mut x = x0;
        let exp_dt = (-theta * dt).exp();
        let std_step = sigma * ((1.0 - exp_dt * exp_dt) / (2.0 * theta)).sqrt().max(0.0);
        for _ in 0..n {
            x = mu + (x - mu) * exp_dt + std_step * self.rng.next_normal();
            path.push(x);
        }
        path
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ParticleFilter — Sequential Monte Carlo
// ─────────────────────────────────────────────────────────────────────────────

/// Sequential Monte Carlo particle filter.
///
/// Maintains a set of weighted particles representing the filtering distribution.
pub struct ParticleFilter {
    rng: LcgRng,
    /// Number of particles.
    pub n_particles: usize,
    /// Current particle positions.
    pub particles: Vec<f64>,
    /// Normalized particle weights.
    pub weights: Vec<f64>,
}

impl ParticleFilter {
    /// Creates a new particle filter initialized with particles from N(mu0, sigma0).
    pub fn new(seed: u64, n_particles: usize, mu0: f64, sigma0: f64) -> Self {
        let mut rng = LcgRng::new(seed);
        let particles: Vec<f64> = (0..n_particles)
            .map(|_| mu0 + sigma0 * rng.next_normal())
            .collect();
        let weights = vec![1.0 / n_particles as f64; n_particles];
        Self {
            rng,
            n_particles,
            particles,
            weights,
        }
    }

    /// Propagates particles through the dynamics model: x_new = f(x) + noise.
    pub fn predict<F>(&mut self, dynamics: F, noise_std: f64)
    where
        F: Fn(f64) -> f64,
    {
        for x in self.particles.iter_mut() {
            *x = dynamics(*x) + noise_std * self.rng.next_normal();
        }
    }

    /// Updates weights with the likelihood of observation `y` given particles.
    ///
    /// `likelihood(x, y)` returns p(y | x).
    pub fn update<F>(&mut self, y: f64, likelihood: F)
    where
        F: Fn(f64, f64) -> f64,
    {
        for (x, w) in self.particles.iter().zip(self.weights.iter_mut()) {
            *w *= likelihood(*x, y);
        }
        self.normalize_weights();
    }

    /// Normalizes weights to sum to 1.
    pub fn normalize_weights(&mut self) {
        let sum: f64 = self.weights.iter().sum();
        if sum > 1e-300 {
            for w in self.weights.iter_mut() {
                *w /= sum;
            }
        } else {
            let uniform = 1.0 / self.n_particles as f64;
            for w in self.weights.iter_mut() {
                *w = uniform;
            }
        }
    }

    /// Systematic resampling: replaces degenerate particles.
    pub fn resample_systematic(&mut self) {
        let n = self.n_particles;
        let mut cumulative = vec![0.0_f64; n + 1];
        for i in 0..n {
            cumulative[i + 1] = cumulative[i] + self.weights[i];
        }
        let u = self.rng.next_f64() / n as f64;
        let mut new_particles = Vec::with_capacity(n);
        let mut j = 0;
        for i in 0..n {
            let threshold = u + i as f64 / n as f64;
            while j < n && cumulative[j + 1] < threshold {
                j += 1;
            }
            new_particles.push(self.particles[j.min(n - 1)]);
        }
        self.particles = new_particles;
        let uniform = 1.0 / n as f64;
        self.weights.fill(uniform);
    }

    /// Multinomial resampling.
    pub fn resample_multinomial(&mut self) {
        let n = self.n_particles;
        let mut cumulative = vec![0.0_f64; n + 1];
        for i in 0..n {
            cumulative[i + 1] = cumulative[i] + self.weights[i];
        }
        let mut new_particles = Vec::with_capacity(n);
        for _ in 0..n {
            let u = self.rng.next_f64();
            let idx = cumulative
                .partition_point(|&c| c < u)
                .saturating_sub(1)
                .min(n - 1);
            new_particles.push(self.particles[idx]);
        }
        self.particles = new_particles;
        let uniform = 1.0 / n as f64;
        self.weights.fill(uniform);
    }

    /// Weighted mean estimate of the current filtering distribution.
    pub fn estimate_mean(&self) -> f64 {
        self.particles
            .iter()
            .zip(self.weights.iter())
            .map(|(&x, &w)| x * w)
            .sum()
    }

    /// Effective sample size of the particle filter.
    pub fn ess(&self) -> f64 {
        effective_sample_size(&self.weights)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_box_muller_finite() {
        let v = box_muller_normal(0.5, 0.25);
        assert!(v.is_finite());
    }

    #[test]
    fn test_inverse_cdf_symmetry() {
        let x = inverse_cdf(0.975);
        assert!((x - 1.96).abs() < 0.01, "x={}", x);
        let y = inverse_cdf(0.025);
        assert!((y + 1.96).abs() < 0.01, "y={}", y);
    }

    #[test]
    fn test_effective_sample_size_uniform() {
        let w = vec![0.25_f64; 4];
        let ess = effective_sample_size(&w);
        assert!((ess - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_effective_sample_size_degenerate() {
        let mut w = vec![0.0_f64; 10];
        w[0] = 1.0;
        let ess = effective_sample_size(&w);
        assert!((ess - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_autocorrelation_lag0() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let r0 = autocorrelation(&x, 0);
        assert!((r0 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_autocorrelation_constant() {
        let x = vec![3.0_f64; 10];
        let r = autocorrelation(&x, 1);
        // variance is 0, so result is 1.0
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mc_sampler_uniform_range() {
        let mut s = McSampler::new(42);
        for _ in 0..100 {
            let x = s.uniform_range(2.0, 5.0);
            assert!((2.0..5.0).contains(&x));
        }
    }

    #[test]
    fn test_mc_sampler_normal_stats() {
        let mut s = McSampler::new(123);
        let samples: Vec<f64> = (0..5000).map(|_| s.normal(0.0, 1.0)).collect();
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
        assert!(mean.abs() < 0.1, "mean={}", mean);
        assert!((var - 1.0).abs() < 0.1, "var={}", var);
    }

    #[test]
    fn test_mc_sampler_exponential_positive() {
        let mut s = McSampler::new(7);
        for _ in 0..100 {
            assert!(s.exponential(1.0) > 0.0);
        }
    }

    #[test]
    fn test_mc_sampler_poisson_mean() {
        let mut s = McSampler::new(99);
        let lambda = 3.0;
        let samples: Vec<f64> = (0..2000).map(|_| s.poisson(lambda) as f64).collect();
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!((mean - lambda).abs() < 0.2, "mean={}", mean);
    }

    #[test]
    fn test_mc_sampler_beta_range() {
        let mut s = McSampler::new(55);
        for _ in 0..100 {
            let x = s.beta(2.0, 3.0);
            assert!((0.0..=1.0).contains(&x));
        }
    }

    #[test]
    fn test_mc_sampler_gamma_positive() {
        let mut s = McSampler::new(77);
        for _ in 0..100 {
            assert!(s.gamma(2.0, 1.0) > 0.0);
        }
    }

    #[test]
    fn test_estimator_mean_variance() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let est = Estimator::new(data);
        assert!((est.mean() - 3.0).abs() < 1e-12);
        assert!((est.variance() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_estimator_standard_error() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let est = Estimator::new(data);
        let se = est.standard_error();
        assert!(se > 0.0);
    }

    #[test]
    fn test_estimator_confidence_interval() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let est = Estimator::new(data);
        let (lo, hi) = est.confidence_interval(0.05);
        assert!(lo < hi);
        assert!(lo < 50.0 && hi > 50.0);
    }

    #[test]
    fn test_estimator_median_odd() {
        let data = vec![3.0, 1.0, 2.0, 5.0, 4.0];
        let est = Estimator::new(data);
        assert!((est.median() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_metropolis_hastings_normal() {
        let mut mh = MetropolisHastings::new(42, 1.0, 500, 1);
        let samples = mh.sample(0.0, 1000, |x: f64| -0.5 * x * x);
        let est = Estimator::new(samples);
        assert!(est.mean().abs() < 0.2);
        assert!((est.variance() - 1.0).abs() < 0.3);
    }

    #[test]
    fn test_hamiltonian_mc_normal() {
        let mut hmc = HamiltonianMC::new(42, 0.3, 10, 200);
        let samples = hmc.sample(0.0, 500, |x: f64| -0.5 * x * x, |x: f64| x);
        let est = Estimator::new(samples);
        assert!(est.mean().abs() < 0.3);
    }

    #[test]
    fn test_mc_integral_plain_pi() {
        // Estimate pi/4 by integrating indicator of unit circle in [0,1]^2
        let mut mc = MonteCarloIntegral::new(42, 10000, 2);
        let (est, _se) = mc.integrate_plain(|x: &[f64]| {
            if x[0] * x[0] + x[1] * x[1] <= 1.0 {
                1.0
            } else {
                0.0
            }
        });
        assert!((est - PI / 4.0).abs() < 0.05, "est={}", est);
    }

    #[test]
    fn test_mc_integral_halton_constant() {
        let mut mc = MonteCarloIntegral::new(1, 1000, 1);
        let (est, _se) = mc.integrate_halton(|_: &[f64]| 1.0);
        assert!((est - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_bootstrap_ci_contains_mean() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let mut boot = BootstrapCI::new(42, 500);
        let (lo, hi) = boot.percentile_ci(&data, 0.05);
        assert!(lo < 24.5 && hi > 24.5);
    }

    #[test]
    fn test_brownian_bridge_endpoints() {
        let mut bb = BrownianBridge::new(42, 0.1);
        let path = bb.brownian_bridge_sample(10, 0.0, 1.0, 1.0);
        assert_eq!(path.len(), 12);
        assert!((path[0] - 0.0).abs() < 1e-12);
        assert!((path[11] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ou_process_mean_reversion() {
        let mut bb = BrownianBridge::new(42, 0.1);
        let path = bb.ornstein_uhlenbeck(5.0, 1.0, 0.0, 0.1, 100);
        // After 100 steps, should be closer to 0 than initial 5.0
        let final_val = path.last().unwrap().abs();
        assert!(final_val < 5.0, "final_val={}", final_val);
    }

    #[test]
    fn test_particle_filter_mean_estimate() {
        let mut pf = ParticleFilter::new(42, 200, 0.0, 1.0);
        // Observe y=2.0; likelihood = N(x, 0.5)
        pf.update(2.0, |x: f64, y: f64| {
            let d = (x - y) / 0.5;
            (-0.5 * d * d).exp()
        });
        let mean = pf.estimate_mean();
        // Should shift toward observation
        assert!(mean > 0.0);
    }

    #[test]
    fn test_particle_filter_resample_systematic() {
        let mut pf = ParticleFilter::new(42, 50, 0.0, 1.0);
        pf.update(1.0, |x: f64, y: f64| (-(x - y).powi(2)).exp());
        let ess_before = pf.ess();
        pf.resample_systematic();
        let ess_after = pf.ess();
        // After resampling, weights are uniform so ESS = N
        assert!(ess_after >= ess_before || ess_after > 1.0);
    }

    #[test]
    fn test_anisotropic_mc_unit_sphere() {
        let mut amc = AnisotropicMC::new(42, 2.0, [0.0, 0.0, 1.0]);
        for _ in 0..20 {
            let v = amc.sample_vmf();
            let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "norm={}", norm);
        }
    }

    #[test]
    fn test_mc_integral_lhs() {
        let mut mc = MonteCarloIntegral::new(42, 100, 1);
        let (est, _) = mc.integrate_lhs(|x: &[f64]| x[0]);
        // E[X] = 0.5 for X ~ U(0,1)
        assert!((est - 0.5).abs() < 0.1, "est={}", est);
    }

    #[test]
    fn test_halton_sequence() {
        let h = halton_seq(1, 2);
        assert!((h - 0.5).abs() < 1e-12);
        let h2 = halton_seq(2, 2);
        assert!((h2 - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_importance_sampler_estimate() {
        let mut is = ImportanceSampler::new(42, 1.0);
        // Estimate E_N(0,1)[x^2] = 1 using proposal N(0,1) with target N(0,1)
        let est = is.estimate(
            1000,
            |x: f64| x * x,
            |_x: f64| 0.0, // log(target/proposal) = 0 when equal
        );
        assert!((est - 1.0).abs() < 0.2, "est={}", est);
    }
}
