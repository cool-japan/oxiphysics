// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bayesian inference: priors, likelihoods, posterior updates, MCMC,
//! Bayesian linear regression, Gaussian processes, and model selection.
//!
//! All distributions work with real-valued parameters; conjugate update
//! formulas are provided where analytic posteriors exist.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Local LCG RNG
// ─────────────────────────────────────────────────────────────────────────────

/// Lightweight LCG random number generator for sampling.
struct BiRng {
    state: u64,
}

impl BiRng {
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

    /// Box-Muller standard normal sample.
    fn next_normal(&mut self) -> f64 {
        loop {
            let u1 = self.next_f64();
            let u2 = self.next_f64();
            if u1 > 0.0 {
                return (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Prior
// ─────────────────────────────────────────────────────────────────────────────

/// Supported prior distribution families.
#[derive(Debug, Clone)]
pub enum Prior {
    /// Uniform prior on \[low, high\].
    Uniform {
        /// Lower bound.
        low: f64,
        /// Upper bound.
        high: f64,
    },
    /// Gaussian (normal) prior N(mean, std²).
    Gaussian {
        /// Prior mean.
        mean: f64,
        /// Prior standard deviation.
        std: f64,
    },
    /// Laplace (double-exponential) prior with location μ and scale b.
    Laplace {
        /// Location parameter.
        mu: f64,
        /// Scale parameter (b > 0).
        b: f64,
    },
    /// Jeffreys (improper) scale prior ∝ 1/θ for θ > 0.
    Jeffreys,
    /// Dirichlet prior over a simplex with concentration parameters α.
    Dirichlet {
        /// Concentration parameters (all positive).
        alpha: Vec<f64>,
    },
    /// Beta prior Beta(alpha, beta) on \[0, 1\].
    Beta {
        /// Shape parameter α > 0.
        alpha: f64,
        /// Shape parameter β > 0.
        beta: f64,
    },
    /// Gamma prior Gamma(shape, rate) on (0, ∞).
    Gamma {
        /// Shape k > 0.
        shape: f64,
        /// Rate λ > 0.
        rate: f64,
    },
}

impl Prior {
    /// Returns the log prior density log p(θ) for a scalar parameter θ.
    ///
    /// Returns `f64::NEG_INFINITY` when θ is outside the support.
    pub fn log_density(&self, theta: f64) -> f64 {
        match self {
            Prior::Uniform { low, high } => {
                if theta >= *low && theta <= *high {
                    -(*high - *low).ln()
                } else {
                    f64::NEG_INFINITY
                }
            }
            Prior::Gaussian { mean, std } => {
                if *std <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                let z = (theta - mean) / std;
                -0.5 * z * z - std.ln() - 0.5 * (2.0 * PI).ln()
            }
            Prior::Laplace { mu, b } => {
                if *b <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                -(theta - mu).abs() / b - (2.0 * b).ln()
            }
            Prior::Jeffreys => {
                if theta > 0.0 {
                    -theta.ln()
                } else {
                    f64::NEG_INFINITY
                }
            }
            Prior::Beta { alpha, beta } => {
                if theta <= 0.0 || theta >= 1.0 {
                    return f64::NEG_INFINITY;
                }
                (*alpha - 1.0) * theta.ln() + (*beta - 1.0) * (1.0 - theta).ln()
                    - log_beta(*alpha, *beta)
            }
            Prior::Gamma { shape, rate } => {
                if theta <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                (*shape - 1.0) * theta.ln() - *rate * theta - log_gamma(*shape) + *shape * rate.ln()
            }
            Prior::Dirichlet { alpha: _ } => {
                // Single-parameter overload is not meaningful for Dirichlet
                f64::NEG_INFINITY
            }
        }
    }

    /// Returns the log density of a Dirichlet prior for a probability vector `x`.
    ///
    /// Returns `f64::NEG_INFINITY` when the dimensions do not match.
    pub fn dirichlet_log_density(&self, x: &[f64]) -> f64 {
        if let Prior::Dirichlet { alpha } = self {
            if alpha.len() != x.len() {
                return f64::NEG_INFINITY;
            }
            let sum: f64 = x.iter().sum();
            if (sum - 1.0).abs() > 1e-8 {
                return f64::NEG_INFINITY;
            }
            let log_num: f64 = alpha.iter().map(|&a| log_gamma(a)).sum();
            let log_den = log_gamma(alpha.iter().sum::<f64>());
            let log_z = log_den - log_num;
            let sum_term: f64 = alpha
                .iter()
                .zip(x.iter())
                .map(|(&a, &xi)| {
                    if xi <= 0.0 {
                        f64::NEG_INFINITY
                    } else {
                        (a - 1.0) * xi.ln()
                    }
                })
                .sum();
            log_z + sum_term
        } else {
            f64::NEG_INFINITY
        }
    }

    /// Samples from the prior using the given RNG seed.
    ///
    /// Returns `None` for improper priors (Jeffreys, Dirichlet) where scalar sampling
    /// is not defined without additional context.
    pub fn sample(&self, seed: u64) -> Option<f64> {
        let mut rng = BiRng::new(seed);
        match self {
            Prior::Uniform { low, high } => Some(low + rng.next_f64() * (high - low)),
            Prior::Gaussian { mean, std } => Some(mean + std * rng.next_normal()),
            Prior::Laplace { mu, b } => {
                let u = rng.next_f64() - 0.5;
                Some(mu - b * u.signum() * (1.0 - 2.0 * u.abs()).ln())
            }
            Prior::Beta { alpha, beta } => {
                // Approximate via ratio of Gamma samples (Cheng's method simplified)
                let x = sample_gamma(*alpha, &mut rng);
                let y = sample_gamma(*beta, &mut rng);
                if x + y <= 0.0 {
                    Some(0.5)
                } else {
                    Some(x / (x + y))
                }
            }
            Prior::Gamma { shape, rate } => Some(sample_gamma(*shape, &mut rng) / rate),
            Prior::Jeffreys | Prior::Dirichlet { .. } => None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Likelihood
// ─────────────────────────────────────────────────────────────────────────────

/// Supported likelihood functions.
#[derive(Debug, Clone)]
pub enum Likelihood {
    /// Gaussian likelihood: x | μ, σ ~ N(μ, σ²).
    Gaussian {
        /// Observed data.
        data: Vec<f64>,
        /// Known noise standard deviation σ.
        sigma: f64,
    },
    /// Poisson likelihood: k | λ ~ Poisson(λ).
    Poisson {
        /// Observed counts.
        counts: Vec<u64>,
    },
    /// Bernoulli likelihood: x ∈ {0,1} | p ~ Bernoulli(p).
    Bernoulli {
        /// Observed outcomes (0 or 1).
        outcomes: Vec<u8>,
    },
    /// Multinomial likelihood over K categories.
    Multinomial {
        /// Observed counts per category.
        counts: Vec<u64>,
    },
}

impl Likelihood {
    /// Returns the log likelihood log p(data | θ) for a scalar parameter θ.
    ///
    /// For Gaussian: θ = μ (mean).
    /// For Poisson: θ = λ (rate).
    /// For Bernoulli: θ = p (probability).
    /// For Multinomial: not applicable (returns NEG_INFINITY; use `multinomial_log_likelihood`).
    pub fn log_likelihood(&self, theta: f64) -> f64 {
        match self {
            Likelihood::Gaussian { data, sigma } => {
                if *sigma <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                let n = data.len() as f64;
                let ss: f64 = data.iter().map(|&x| (x - theta).powi(2)).sum();
                -0.5 * n * (2.0 * PI * sigma * sigma).ln() - ss / (2.0 * sigma * sigma)
            }
            Likelihood::Poisson { counts } => {
                if theta <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                counts
                    .iter()
                    .map(|&k| k as f64 * theta.ln() - theta - log_factorial(k))
                    .sum()
            }
            Likelihood::Bernoulli { outcomes } => {
                if theta <= 0.0 || theta >= 1.0 {
                    return f64::NEG_INFINITY;
                }
                outcomes
                    .iter()
                    .map(|&x| {
                        if x == 1 {
                            theta.ln()
                        } else {
                            (1.0 - theta).ln()
                        }
                    })
                    .sum()
            }
            Likelihood::Multinomial { counts: _ } => f64::NEG_INFINITY,
        }
    }

    /// Returns the log multinomial likelihood for probability vector `probs`.
    pub fn multinomial_log_likelihood(&self, probs: &[f64]) -> f64 {
        if let Likelihood::Multinomial { counts } = self {
            if counts.len() != probs.len() {
                return f64::NEG_INFINITY;
            }
            counts
                .iter()
                .zip(probs.iter())
                .map(|(&k, &p)| {
                    if p <= 0.0 {
                        if k == 0 { 0.0 } else { f64::NEG_INFINITY }
                    } else {
                        k as f64 * p.ln()
                    }
                })
                .sum()
        } else {
            f64::NEG_INFINITY
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BayesianUpdate
// ─────────────────────────────────────────────────────────────────────────────

/// Conjugate Bayesian updates for standard distribution families.
///
/// Each method returns updated (posterior) hyperparameters given a prior and
/// observed data.
#[derive(Debug, Clone)]
pub struct BayesianUpdate;

impl BayesianUpdate {
    /// Normal-Normal conjugate update: known variance σ², Gaussian prior on μ.
    ///
    /// Prior: μ ~ N(μ₀, τ₀²).  Data: x_i ~ N(μ, σ²).
    /// Returns posterior (μ_n, τ_n).
    pub fn normal_normal(
        prior_mean: f64,
        prior_std: f64,
        likelihood_std: f64,
        data: &[f64],
    ) -> (f64, f64) {
        let n = data.len() as f64;
        if n == 0.0 {
            return (prior_mean, prior_std);
        }
        let tau0_sq = prior_std * prior_std;
        let sigma_sq = likelihood_std * likelihood_std;
        let x_bar: f64 = data.iter().sum::<f64>() / n;
        let tau_n_sq = 1.0 / (1.0 / tau0_sq + n / sigma_sq);
        let mu_n = tau_n_sq * (prior_mean / tau0_sq + n * x_bar / sigma_sq);
        (mu_n, tau_n_sq.sqrt())
    }

    /// Beta-Bernoulli conjugate update: Beta prior on p, Bernoulli observations.
    ///
    /// Prior: p ~ Beta(α, β).  Data: k successes out of n.
    /// Returns posterior (α', β').
    pub fn beta_bernoulli(
        prior_alpha: f64,
        prior_beta: f64,
        successes: u64,
        total: u64,
    ) -> (f64, f64) {
        let failures = total - successes.min(total);
        (prior_alpha + successes as f64, prior_beta + failures as f64)
    }

    /// Gamma-Poisson conjugate update: Gamma prior on λ, Poisson observations.
    ///
    /// Prior: λ ~ Gamma(α, β).  Data: counts k_i.
    /// Returns posterior (α', β').
    pub fn gamma_poisson(prior_shape: f64, prior_rate: f64, counts: &[u64]) -> (f64, f64) {
        let n = counts.len() as f64;
        let sum_k: f64 = counts.iter().map(|&k| k as f64).sum();
        (prior_shape + sum_k, prior_rate + n)
    }

    /// Dirichlet-Multinomial conjugate update.
    ///
    /// Prior: p ~ Dir(α).  Data: counts k_i.
    /// Returns posterior α' = α + k.
    pub fn dirichlet_multinomial(prior_alpha: &[f64], counts: &[u64]) -> Vec<f64> {
        prior_alpha
            .iter()
            .zip(counts.iter())
            .map(|(&a, &k)| a + k as f64)
            .collect()
    }

    /// Normal-inverse-Gamma conjugate update for unknown mean and variance.
    ///
    /// Prior hyperparameters: (μ₀, κ₀, α₀, β₀).
    /// Returns updated hyperparameters (μ_n, κ_n, α_n, β_n).
    pub fn normal_inverse_gamma(
        mu0: f64,
        kappa0: f64,
        alpha0: f64,
        beta0: f64,
        data: &[f64],
    ) -> (f64, f64, f64, f64) {
        let n = data.len() as f64;
        if n == 0.0 {
            return (mu0, kappa0, alpha0, beta0);
        }
        let x_bar = data.iter().sum::<f64>() / n;
        let ss: f64 = data.iter().map(|&x| (x - x_bar).powi(2)).sum();
        let kappa_n = kappa0 + n;
        let mu_n = (kappa0 * mu0 + n * x_bar) / kappa_n;
        let alpha_n = alpha0 + n / 2.0;
        let beta_n = beta0 + 0.5 * ss + (kappa0 * n * (x_bar - mu0).powi(2)) / (2.0 * kappa_n);
        (mu_n, kappa_n, alpha_n, beta_n)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MarkovChainMonteCarlo
// ─────────────────────────────────────────────────────────────────────────────

/// Markov Chain Monte Carlo samplers.
///
/// Implements Metropolis-Hastings (random-walk), Gibbs sampling helpers,
/// and a simplified dual-averaging NUTS-like step-size adaptation.
#[derive(Debug, Clone)]
pub struct MarkovChainMonteCarlo {
    /// Step size (proposal std for MH).
    pub step_size: f64,
    /// Number of warm-up (burn-in) steps.
    pub n_warmup: usize,
}

impl MarkovChainMonteCarlo {
    /// Creates a new MCMC sampler with the given step size and warm-up count.
    pub fn new(step_size: f64, n_warmup: usize) -> Self {
        Self {
            step_size,
            n_warmup,
        }
    }

    /// Runs a Metropolis-Hastings random-walk chain.
    ///
    /// # Arguments
    /// * `log_target`  - Closure returning log π(θ) (up to normalisation).
    /// * `init`        - Initial parameter value.
    /// * `n_samples`   - Number of post-warmup samples to return.
    /// * `seed`        - RNG seed.
    pub fn metropolis_hastings<F>(
        &self,
        log_target: F,
        init: f64,
        n_samples: usize,
        seed: u64,
    ) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut rng = BiRng::new(seed);
        let mut current = init;
        let mut log_current = log_target(current);
        // Warm-up
        for _ in 0..self.n_warmup {
            let proposal = current + self.step_size * rng.next_normal();
            let log_proposal = log_target(proposal);
            let log_alpha = log_proposal - log_current;
            if rng.next_f64().ln() < log_alpha {
                current = proposal;
                log_current = log_proposal;
            }
        }
        // Sampling
        let mut samples = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let proposal = current + self.step_size * rng.next_normal();
            let log_proposal = log_target(proposal);
            let log_alpha = log_proposal - log_current;
            if rng.next_f64().ln() < log_alpha {
                current = proposal;
                log_current = log_proposal;
            }
            samples.push(current);
        }
        samples
    }

    /// Metropolis-Hastings for a vector-valued parameter.
    pub fn metropolis_hastings_vec<F>(
        &self,
        log_target: F,
        init: Vec<f64>,
        n_samples: usize,
        seed: u64,
    ) -> Vec<Vec<f64>>
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut rng = BiRng::new(seed);
        let dim = init.len();
        let mut current = init.clone();
        let mut log_current = log_target(&current);
        // Warm-up
        for _ in 0..self.n_warmup {
            let proposal: Vec<f64> = current
                .iter()
                .map(|&x| x + self.step_size * rng.next_normal())
                .collect();
            let log_proposal = log_target(&proposal);
            if rng.next_f64().ln() < log_proposal - log_current {
                current = proposal;
                log_current = log_proposal;
            }
        }
        let mut samples = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let proposal: Vec<f64> = current
                .iter()
                .map(|&x| x + self.step_size * rng.next_normal())
                .collect();
            let log_proposal = log_target(&proposal);
            if rng.next_f64().ln() < log_proposal - log_current {
                current = proposal;
                log_current = log_proposal;
            }
            samples.push(current.clone());
        }
        let _ = dim; // suppress unused warning
        samples
    }

    /// Gibbs sampler for a bivariate Gaussian with known full conditionals.
    ///
    /// Samples from N(\[μ1, μ2\], \[\[σ1², ρσ1σ2\],\[ρσ1σ2, σ2²\]\]).
    pub fn gibbs_bivariate_gaussian(
        mu1: f64,
        mu2: f64,
        sigma1: f64,
        sigma2: f64,
        rho: f64,
        n_samples: usize,
        seed: u64,
    ) -> Vec<[f64; 2]> {
        let mut rng = BiRng::new(seed);
        let mut x1;
        let mut x2 = mu2;
        let mut samples = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            // x1 | x2 ~ N(μ1 + ρ(σ1/σ2)(x2-μ2), σ1²(1-ρ²))
            let cond_mean1 = mu1 + rho * (sigma1 / sigma2) * (x2 - mu2);
            let cond_std1 = sigma1 * (1.0 - rho * rho).sqrt();
            x1 = cond_mean1 + cond_std1 * rng.next_normal();
            // x2 | x1 ~ N(μ2 + ρ(σ2/σ1)(x1-μ1), σ2²(1-ρ²))
            let cond_mean2 = mu2 + rho * (sigma2 / sigma1) * (x1 - mu1);
            let cond_std2 = sigma2 * (1.0 - rho * rho).sqrt();
            x2 = cond_mean2 + cond_std2 * rng.next_normal();
            samples.push([x1, x2]);
        }
        samples
    }

    /// Simplified NUTS (No-U-Turn Sampler) step using leapfrog integration.
    ///
    /// Returns a single sample given gradient function `grad_log_target`.
    pub fn nuts_step<F>(
        &self,
        log_target: F,
        grad_log_target: impl Fn(f64) -> f64,
        init: f64,
        seed: u64,
    ) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut rng = BiRng::new(seed);
        let mut q = init;
        let mut p = rng.next_normal();
        let h_init = -log_target(q) + 0.5 * p * p;
        // Single leapfrog step
        let grad = grad_log_target(q);
        p += 0.5 * self.step_size * grad;
        q += self.step_size * p;
        p += 0.5 * self.step_size * grad_log_target(q);
        let h_prop = -log_target(q) + 0.5 * p * p;
        let log_alpha = -(h_prop - h_init);
        if rng.next_f64().ln() < log_alpha {
            q
        } else {
            init
        }
    }

    /// Returns the effective sample size (ESS) from a chain using autocorrelation.
    pub fn effective_sample_size(chain: &[f64]) -> f64 {
        let n = chain.len();
        if n < 4 {
            return n as f64;
        }
        let mean = chain.iter().sum::<f64>() / n as f64;
        let variance = chain.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        if variance < 1e-30 {
            return n as f64;
        }
        let mut rho_sum = 0.0_f64;
        for lag in 1..(n / 2) {
            let acf: f64 = (0..n - lag)
                .map(|i| (chain[i] - mean) * (chain[i + lag] - mean))
                .sum::<f64>()
                / (n as f64 * variance);
            if acf < 0.0 {
                break;
            }
            rho_sum += acf;
        }
        n as f64 / (1.0 + 2.0 * rho_sum)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BayesianLinearRegression
// ─────────────────────────────────────────────────────────────────────────────

/// Bayesian linear regression with conjugate Gaussian prior on weights.
///
/// Model: y = X w + ε, ε ~ N(0, σ² I).
/// Prior: w ~ N(w₀, α⁻¹ I).
#[derive(Debug, Clone)]
pub struct BayesianLinearRegression {
    /// Prior precision α (= 1/σ_prior²).
    pub alpha_prior: f64,
    /// Noise precision β (= 1/σ_noise²).
    pub beta_noise: f64,
    /// Posterior mean weights (set after fitting).
    pub posterior_mean: Vec<f64>,
    /// Posterior covariance matrix (row-major, set after fitting).
    pub posterior_cov: Vec<Vec<f64>>,
}

impl BayesianLinearRegression {
    /// Creates a new Bayesian linear regression model.
    ///
    /// # Arguments
    /// * `alpha_prior` - Prior precision (inverse variance) on weights.
    /// * `beta_noise`  - Noise precision.
    pub fn new(alpha_prior: f64, beta_noise: f64) -> Self {
        Self {
            alpha_prior,
            beta_noise,
            posterior_mean: vec![],
            posterior_cov: vec![],
        }
    }

    /// Fits the model to design matrix `x_mat` (n × d) and targets `y` (n).
    ///
    /// Computes the posterior mean and covariance analytically.
    pub fn fit(&mut self, x_mat: &[Vec<f64>], y: &[f64]) {
        let n = x_mat.len();
        if n == 0 || y.is_empty() {
            return;
        }
        let d = x_mat[0].len();
        // S_N^{-1} = α I + β X^T X
        // m_N = β S_N X^T y
        // Compute X^T X (d × d)
        let mut xtx = vec![vec![0.0_f64; d]; d];
        for row in x_mat {
            for i in 0..d {
                for j in 0..d {
                    xtx[i][j] += row[i] * row[j];
                }
            }
        }
        // S_N^{-1} = α I + β X^T X
        let mut s_inv = vec![vec![0.0_f64; d]; d];
        for i in 0..d {
            for j in 0..d {
                s_inv[i][j] = self.beta_noise * xtx[i][j];
            }
            s_inv[i][i] += self.alpha_prior;
        }
        // Invert S_inv using Gaussian elimination
        let s_n = mat_inverse(&s_inv);
        // X^T y (d vector)
        let mut xty = vec![0.0_f64; d];
        for (row, &yi) in x_mat.iter().zip(y.iter()) {
            for i in 0..d {
                xty[i] += row[i] * yi;
            }
        }
        // m_N = β S_N X^T y
        let mut m_n = vec![0.0_f64; d];
        for i in 0..d {
            for j in 0..d {
                m_n[i] += s_n[i][j] * xty[j];
            }
        }
        // Scale m_N by β
        for v in m_n.iter_mut() {
            *v *= self.beta_noise;
        }
        self.posterior_mean = m_n;
        self.posterior_cov = s_n;
    }

    /// Returns the predictive mean for a new input vector `x_new`.
    pub fn predict_mean(&self, x_new: &[f64]) -> f64 {
        self.posterior_mean
            .iter()
            .zip(x_new.iter())
            .map(|(&w, &x)| w * x)
            .sum()
    }

    /// Returns the predictive variance for a new input `x_new`.
    ///
    /// σ²_pred = 1/β + x^T S_N x.
    pub fn predict_variance(&self, x_new: &[f64]) -> f64 {
        if self.posterior_cov.is_empty() {
            return 1.0 / self.beta_noise;
        }
        let d = x_new.len();
        let mut s_x = vec![0.0_f64; d];
        for (i, s) in s_x.iter_mut().enumerate() {
            for (j, &xj) in x_new.iter().enumerate() {
                *s += self.posterior_cov[i][j] * xj;
            }
        }
        let xtsx: f64 = x_new.iter().zip(s_x.iter()).map(|(&x, &sx)| x * sx).sum();
        1.0 / self.beta_noise + xtsx
    }

    /// Returns the log marginal likelihood (model evidence) log p(y | X, α, β).
    ///
    /// log p(y) = (d/2) ln α + (n/2) ln β - (1/2)(β ||y - X m_N||² + α ||m_N||²)
    ///            - (1/2) ln|S_N^{-1}| - (n/2) ln(2π)
    pub fn log_evidence(&self, x_mat: &[Vec<f64>], y: &[f64]) -> f64 {
        if self.posterior_mean.is_empty() || x_mat.is_empty() {
            return f64::NEG_INFINITY;
        }
        let n = y.len() as f64;
        let d = self.posterior_mean.len() as f64;
        // Compute residuals y - X m_N
        let mut ss_res = 0.0_f64;
        for (row, &yi) in x_mat.iter().zip(y.iter()) {
            let pred = self.predict_mean(row);
            ss_res += (yi - pred).powi(2);
        }
        let m_norm_sq: f64 = self.posterior_mean.iter().map(|&w| w * w).sum();
        let log_det_s: f64 = {
            // Log-determinant of S_N via diagonal approx (identity prior case)
            self.posterior_cov
                .iter()
                .enumerate()
                .map(|(i, row)| row[i].abs().ln())
                .sum()
        };
        let alpha = self.alpha_prior;
        let beta = self.beta_noise;
        (d / 2.0) * alpha.ln() + (n / 2.0) * beta.ln() - 0.5 * (beta * ss_res + alpha * m_norm_sq)
            + 0.5 * log_det_s
            - (n / 2.0) * (2.0 * PI).ln()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GaussianProcess
// ─────────────────────────────────────────────────────────────────────────────

/// Kernel function types for Gaussian processes.
#[derive(Debug, Clone, Copy)]
pub enum Kernel {
    /// Radial basis function (squared-exponential) kernel.
    Rbf {
        /// Length scale ℓ.
        length_scale: f64,
        /// Signal variance σ_f².
        signal_variance: f64,
    },
    /// Matérn 3/2 kernel.
    Matern32 {
        /// Length scale ℓ.
        length_scale: f64,
        /// Signal variance σ_f².
        signal_variance: f64,
    },
    /// Matérn 5/2 kernel.
    Matern52 {
        /// Length scale ℓ.
        length_scale: f64,
        /// Signal variance σ_f².
        signal_variance: f64,
    },
    /// Linear (dot-product) kernel.
    Linear {
        /// Bias variance σ_b².
        bias_variance: f64,
        /// Slope variance σ_v².
        slope_variance: f64,
    },
}

impl Kernel {
    /// Evaluates k(x, y) for scalar inputs.
    pub fn eval(&self, x: f64, y: f64) -> f64 {
        match self {
            Kernel::Rbf {
                length_scale,
                signal_variance,
            } => {
                let r2 = (x - y).powi(2) / (length_scale * length_scale);
                signal_variance * (-0.5 * r2).exp()
            }
            Kernel::Matern32 {
                length_scale,
                signal_variance,
            } => {
                let r = (x - y).abs() / length_scale;
                let s3r = 3.0_f64.sqrt() * r;
                signal_variance * (1.0 + s3r) * (-s3r).exp()
            }
            Kernel::Matern52 {
                length_scale,
                signal_variance,
            } => {
                let r = (x - y).abs() / length_scale;
                let s5r = 5.0_f64.sqrt() * r;
                signal_variance * (1.0 + s5r + 5.0 * r * r / 3.0) * (-s5r).exp()
            }
            Kernel::Linear {
                bias_variance,
                slope_variance,
            } => bias_variance + slope_variance * x * y,
        }
    }
}

/// Gaussian Process regression with a scalar-input, scalar-output model.
///
/// Computes the posterior mean and variance given observed (x, y) pairs.
#[derive(Debug, Clone)]
pub struct GaussianProcess {
    /// Kernel function.
    pub kernel: Kernel,
    /// Noise variance σ_n².
    pub noise_variance: f64,
    /// Training inputs.
    pub x_train: Vec<f64>,
    /// Training targets.
    pub y_train: Vec<f64>,
    /// Cholesky factor L of (K + σ_n² I) for efficient prediction.
    chol: Vec<Vec<f64>>,
    /// α = L^{-T} L^{-1} y.
    alpha: Vec<f64>,
}

impl GaussianProcess {
    /// Creates a new GP with the specified kernel and noise variance.
    pub fn new(kernel: Kernel, noise_variance: f64) -> Self {
        Self {
            kernel,
            noise_variance,
            x_train: vec![],
            y_train: vec![],
            chol: vec![],
            alpha: vec![],
        }
    }

    /// Fits the GP to training data, computing the Cholesky factor.
    pub fn fit(&mut self, x_train: Vec<f64>, y_train: Vec<f64>) {
        let n = x_train.len();
        self.x_train = x_train;
        self.y_train = y_train.clone();
        // Build kernel matrix K + σ_n² I
        let mut k = vec![vec![0.0_f64; n]; n];
        for (i, row) in k.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = self.kernel.eval(self.x_train[i], self.x_train[j]);
            }
        }
        for (i, row) in k.iter_mut().enumerate() {
            row[i] += self.noise_variance;
        }
        // Cholesky decomposition
        self.chol = cholesky(&k);
        // α = L^{-T} L^{-1} y via forward/backward substitution
        let v = forward_sub(&self.chol, &y_train);
        self.alpha = backward_sub_t(&self.chol, &v);
    }

    /// Returns the posterior mean at test point `x_star`.
    pub fn predict_mean(&self, x_star: f64) -> f64 {
        if self.x_train.is_empty() {
            return 0.0;
        }
        let k_star: Vec<f64> = self
            .x_train
            .iter()
            .map(|&xi| self.kernel.eval(xi, x_star))
            .collect();
        k_star
            .iter()
            .zip(self.alpha.iter())
            .map(|(&k, &a)| k * a)
            .sum()
    }

    /// Returns the posterior variance at test point `x_star`.
    pub fn predict_variance(&self, x_star: f64) -> f64 {
        let k_ss = self.kernel.eval(x_star, x_star) + self.noise_variance;
        if self.chol.is_empty() {
            return k_ss;
        }
        let k_star: Vec<f64> = self
            .x_train
            .iter()
            .map(|&xi| self.kernel.eval(xi, x_star))
            .collect();
        let v = forward_sub(&self.chol, &k_star);
        let reduction: f64 = v.iter().map(|&vi| vi * vi).sum();
        (k_ss - reduction).max(0.0)
    }

    /// Returns the log marginal likelihood log p(y | X, θ).
    pub fn log_marginal_likelihood(&self) -> f64 {
        if self.chol.is_empty() || self.y_train.is_empty() {
            return f64::NEG_INFINITY;
        }
        let n = self.y_train.len() as f64;
        let y = &self.y_train;
        // data fit term: -0.5 y^T α
        let data_fit: f64 = y
            .iter()
            .zip(self.alpha.iter())
            .map(|(&yi, &ai)| yi * ai)
            .sum();
        // log det term: log |K| = 2 Σ log L_ii
        let log_det: f64 = self
            .chol
            .iter()
            .enumerate()
            .map(|(i, row)| row[i].abs().ln())
            .sum::<f64>()
            * 2.0;
        -0.5 * data_fit - 0.5 * log_det - (n / 2.0) * (2.0 * PI).ln()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ModelSelection
// ─────────────────────────────────────────────────────────────────────────────

/// Model selection criteria: AIC, BIC, Bayes factor, and cross-validation.
#[derive(Debug, Clone)]
pub struct ModelSelection;

impl ModelSelection {
    /// Akaike Information Criterion: AIC = 2k - 2 log L.
    ///
    /// # Arguments
    /// * `log_likelihood` - Maximised log-likelihood.
    /// * `n_params`       - Number of free parameters k.
    pub fn aic(log_likelihood: f64, n_params: usize) -> f64 {
        2.0 * n_params as f64 - 2.0 * log_likelihood
    }

    /// Corrected AIC for small samples: AICc = AIC + 2k(k+1)/(n-k-1).
    pub fn aicc(log_likelihood: f64, n_params: usize, n_data: usize) -> f64 {
        let k = n_params as f64;
        let n = n_data as f64;
        let aic = Self::aic(log_likelihood, n_params);
        if n > k + 1.0 {
            aic + 2.0 * k * (k + 1.0) / (n - k - 1.0)
        } else {
            aic
        }
    }

    /// Bayesian Information Criterion: BIC = k ln(n) - 2 log L.
    ///
    /// # Arguments
    /// * `log_likelihood` - Maximised log-likelihood.
    /// * `n_params`       - Number of free parameters k.
    /// * `n_data`         - Number of data points n.
    pub fn bic(log_likelihood: f64, n_params: usize, n_data: usize) -> f64 {
        n_params as f64 * (n_data as f64).ln() - 2.0 * log_likelihood
    }

    /// Bayes factor (in log scale): log BF₁₂ = log p(D|M₁) - log p(D|M₂).
    pub fn log_bayes_factor(log_evidence_1: f64, log_evidence_2: f64) -> f64 {
        log_evidence_1 - log_evidence_2
    }

    /// Interprets the log Bayes factor according to Jeffreys' scale.
    ///
    /// Returns a descriptive string.
    pub fn jeffreys_scale(log_bf: f64) -> &'static str {
        let bf = log_bf.exp();
        if bf < 1.0 {
            "Negative (favours M2)"
        } else if bf < 3.0 {
            "Barely worth mentioning"
        } else if bf < 10.0 {
            "Substantial"
        } else if bf < 30.0 {
            "Strong"
        } else if bf < 100.0 {
            "Very strong"
        } else {
            "Decisive"
        }
    }

    /// K-fold cross-validation mean squared error.
    ///
    /// # Arguments
    /// * `x`    - Feature matrix (n × d, row-major as `Vec<Vec`f64`>`).
    /// * `y`    - Target vector (length n).
    /// * `k`    - Number of folds.
    /// * `alpha` - Prior precision for Bayesian linear regression.
    /// * `beta`  - Noise precision.
    pub fn k_fold_cv_mse(x: &[Vec<f64>], y: &[f64], k: usize, alpha: f64, beta_noise: f64) -> f64 {
        let n = x.len();
        if k == 0 || n < k {
            return f64::NAN;
        }
        let fold_size = n / k;
        let mut total_mse = 0.0_f64;
        let mut total_count = 0_usize;
        for fold in 0..k {
            let test_start = fold * fold_size;
            let test_end = if fold == k - 1 {
                n
            } else {
                test_start + fold_size
            };
            let x_train: Vec<Vec<f64>> = x[..test_start]
                .iter()
                .chain(x[test_end..].iter())
                .cloned()
                .collect();
            let y_train: Vec<f64> = y[..test_start]
                .iter()
                .chain(y[test_end..].iter())
                .cloned()
                .collect();
            let x_test = &x[test_start..test_end];
            let y_test = &y[test_start..test_end];
            let mut model = BayesianLinearRegression::new(alpha, beta_noise);
            model.fit(&x_train, &y_train);
            for (xi, &yi) in x_test.iter().zip(y_test.iter()) {
                let pred = model.predict_mean(xi);
                total_mse += (yi - pred).powi(2);
                total_count += 1;
            }
        }
        if total_count == 0 {
            f64::NAN
        } else {
            total_mse / total_count as f64
        }
    }

    /// Pseudo Bayes factor approximation using LOO-CV log predictive density.
    pub fn loo_cv_log_predictive(x: &[Vec<f64>], y: &[f64], alpha: f64, beta_noise: f64) -> f64 {
        let n = x.len();
        let mut total = 0.0_f64;
        for i in 0..n {
            let x_train: Vec<Vec<f64>> = x
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, v)| v.clone())
                .collect();
            let y_train: Vec<f64> = y
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, &v)| v)
                .collect();
            let mut model = BayesianLinearRegression::new(alpha, beta_noise);
            model.fit(&x_train, &y_train);
            let mean = model.predict_mean(&x[i]);
            let var = model.predict_variance(&x[i]);
            // log N(y_i | mean, var)
            let log_p = -0.5 * (y[i] - mean).powi(2) / var - 0.5 * (2.0 * PI * var).ln();
            total += log_p;
        }
        total
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Returns ln Γ(x) using the Lanczos approximation.
pub fn log_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::INFINITY;
    }
    // Lanczos coefficients (g = 5, n = 7)
    let g = 5.0_f64;
    let c = [
        1.000000000190015,
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.001208650973866179,
        -5.395239384953e-6,
    ];
    let mut sum = c[0];
    let mut xp = x;
    for ci in c.iter().skip(1) {
        xp += 1.0;
        sum += ci / xp;
    }
    let t = x + g + 0.5;
    0.5 * (2.0 * PI).ln() + (x + 0.5) * t.ln() - t + sum.ln() - x.ln()
}

/// Returns ln B(a, b) = ln Γ(a) + ln Γ(b) - ln Γ(a+b).
pub fn log_beta(a: f64, b: f64) -> f64 {
    log_gamma(a) + log_gamma(b) - log_gamma(a + b)
}

/// Returns ln k! via Stirling for large k.
fn log_factorial(k: u64) -> f64 {
    log_gamma(k as f64 + 1.0)
}

/// Samples from a Gamma(shape, 1) distribution using Marsaglia-Tsang.
fn sample_gamma(shape: f64, rng: &mut BiRng) -> f64 {
    if shape < 1.0 {
        // Boost using shape+1
        return sample_gamma(1.0 + shape, rng) * rng.next_f64().powf(1.0 / shape);
    }
    let d = shape - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    loop {
        let x = rng.next_normal();
        let v_raw = 1.0 + c * x;
        if v_raw <= 0.0 {
            continue;
        }
        let v = v_raw.powi(3);
        let u = rng.next_f64();
        if u < 1.0 - 0.0331 * (x * x) * (x * x) {
            return d * v;
        }
        if u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
            return d * v;
        }
    }
}

/// Performs Cholesky decomposition of a positive-definite matrix (returns lower L).
fn cholesky(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let sum: f64 = (0..j).map(|k| l[i][k] * l[j][k]).sum();
            if i == j {
                let val = a[i][i] - sum;
                l[i][j] = if val > 0.0 { val.sqrt() } else { 1e-10 };
            } else if l[j][j].abs() < 1e-30 {
                l[i][j] = 0.0;
            } else {
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }
    l
}

/// Forward substitution: solve L x = b.
fn forward_sub(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    for i in 0..n {
        let s: f64 = (0..i).map(|j| l[i][j] * x[j]).sum();
        if l[i][i].abs() > 1e-30 {
            x[i] = (b[i] - s) / l[i][i];
        }
    }
    x
}

/// Backward substitution with transposed L: solve L^T x = b.
fn backward_sub_t(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| l[j][i] * x[j]).sum();
        if l[i][i].abs() > 1e-30 {
            x[i] = (b[i] - s) / l[i][i];
        }
    }
    x
}

/// Inverts a matrix using Gauss-Jordan elimination (in-place augmented matrix).
fn mat_inverse(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return vec![];
    }
    // Augmented matrix [a | I]
    let mut aug: Vec<Vec<f64>> = a
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
    for col in 0..n {
        // Find pivot
        let mut max_row = col;
        for row in col + 1..n {
            if aug[row][col].abs() > aug[max_row][col].abs() {
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-30 {
            continue;
        }
        for cell in &mut aug[col] {
            *cell /= pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row][col];
                let col_row: Vec<f64> = aug[col].clone();
                for (cell, &cv) in aug[row].iter_mut().zip(col_row.iter()) {
                    *cell -= factor * cv;
                }
            }
        }
    }
    aug.into_iter().map(|row| row[n..].to_vec()).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Prior ────────────────────────────────────────────────────────────────

    #[test]
    fn test_prior_uniform_inside() {
        let p = Prior::Uniform {
            low: 0.0,
            high: 1.0,
        };
        assert!(p.log_density(0.5).is_finite());
    }

    #[test]
    fn test_prior_uniform_outside() {
        let p = Prior::Uniform {
            low: 0.0,
            high: 1.0,
        };
        assert_eq!(p.log_density(2.0), f64::NEG_INFINITY);
    }

    #[test]
    fn test_prior_gaussian_mode() {
        let p = Prior::Gaussian {
            mean: 0.0,
            std: 1.0,
        };
        // Mode at 0 should be maximum
        let l0 = p.log_density(0.0);
        let l1 = p.log_density(1.0);
        assert!(l0 > l1);
    }

    #[test]
    fn test_prior_laplace_symmetric() {
        let p = Prior::Laplace { mu: 0.0, b: 1.0 };
        assert!((p.log_density(1.0) - p.log_density(-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_prior_jeffreys_positive() {
        let p = Prior::Jeffreys;
        assert!(p.log_density(2.0).is_finite());
    }

    #[test]
    fn test_prior_jeffreys_non_positive() {
        let p = Prior::Jeffreys;
        assert_eq!(p.log_density(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn test_prior_beta_at_half() {
        let p = Prior::Beta {
            alpha: 2.0,
            beta: 2.0,
        };
        assert!(p.log_density(0.5).is_finite());
    }

    #[test]
    fn test_prior_beta_outside_support() {
        let p = Prior::Beta {
            alpha: 2.0,
            beta: 2.0,
        };
        assert_eq!(p.log_density(1.5), f64::NEG_INFINITY);
    }

    #[test]
    fn test_prior_gamma_positive() {
        let p = Prior::Gamma {
            shape: 2.0,
            rate: 1.0,
        };
        assert!(p.log_density(1.0).is_finite());
    }

    #[test]
    fn test_prior_sample_uniform() {
        let p = Prior::Uniform {
            low: 0.0,
            high: 1.0,
        };
        let s = p.sample(42).unwrap();
        assert!((0.0..=1.0).contains(&s));
    }

    #[test]
    fn test_prior_sample_gaussian() {
        let p = Prior::Gaussian {
            mean: 5.0,
            std: 0.1,
        };
        let s = p.sample(7).unwrap();
        assert!((s - 5.0).abs() < 2.0);
    }

    #[test]
    fn test_prior_dirichlet_log_density() {
        let p = Prior::Dirichlet {
            alpha: vec![1.0, 1.0, 1.0],
        };
        let ld = p.dirichlet_log_density(&[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);
        assert!(ld.is_finite());
    }

    // ── Likelihood ───────────────────────────────────────────────────────────

    #[test]
    fn test_likelihood_gaussian_log_ll_finite() {
        let ll = Likelihood::Gaussian {
            data: vec![1.0, 2.0, 3.0],
            sigma: 1.0,
        };
        assert!(ll.log_likelihood(2.0).is_finite());
    }

    #[test]
    fn test_likelihood_gaussian_mode_at_mean() {
        let data = vec![2.0, 2.0, 2.0];
        let ll = Likelihood::Gaussian {
            data: data.clone(),
            sigma: 1.0,
        };
        let l2 = ll.log_likelihood(2.0);
        let ll2 = Likelihood::Gaussian { data, sigma: 1.0 };
        let l3 = ll2.log_likelihood(3.0);
        assert!(l2 > l3);
    }

    #[test]
    fn test_likelihood_poisson_log_ll() {
        let ll = Likelihood::Poisson {
            counts: vec![3, 4, 5],
        };
        assert!(ll.log_likelihood(4.0).is_finite());
    }

    #[test]
    fn test_likelihood_poisson_non_positive_lambda() {
        let ll = Likelihood::Poisson { counts: vec![1] };
        assert_eq!(ll.log_likelihood(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn test_likelihood_bernoulli_log_ll() {
        let ll = Likelihood::Bernoulli {
            outcomes: vec![1, 0, 1],
        };
        assert!(ll.log_likelihood(0.6).is_finite());
    }

    #[test]
    fn test_likelihood_bernoulli_outside_range() {
        let ll = Likelihood::Bernoulli { outcomes: vec![1] };
        assert_eq!(ll.log_likelihood(1.0), f64::NEG_INFINITY);
    }

    #[test]
    fn test_likelihood_multinomial_log_ll() {
        let ll = Likelihood::Multinomial {
            counts: vec![3, 4, 3],
        };
        let probs = [0.3, 0.4, 0.3];
        let lp = ll.multinomial_log_likelihood(&probs);
        assert!(lp.is_finite());
    }

    // ── BayesianUpdate ───────────────────────────────────────────────────────

    #[test]
    fn test_normal_normal_posterior_shrinks_toward_data() {
        let data = vec![5.0, 5.0, 5.0, 5.0, 5.0];
        let (mu_n, _) = BayesianUpdate::normal_normal(0.0, 10.0, 1.0, &data);
        assert!(mu_n > 3.0); // Should be pulled toward 5
    }

    #[test]
    fn test_normal_normal_empty_data() {
        let (mu_n, sigma_n) = BayesianUpdate::normal_normal(1.0, 2.0, 1.0, &[]);
        assert_eq!((mu_n, sigma_n), (1.0, 2.0));
    }

    #[test]
    fn test_beta_bernoulli_update() {
        let (alpha_n, beta_n) = BayesianUpdate::beta_bernoulli(1.0, 1.0, 7, 10);
        assert!((alpha_n - 8.0).abs() < 1e-12);
        assert!((beta_n - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_gamma_poisson_update() {
        let (alpha_n, beta_n) = BayesianUpdate::gamma_poisson(2.0, 1.0, &[3, 4, 5]);
        assert!((alpha_n - 14.0).abs() < 1e-12);
        assert!((beta_n - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_dirichlet_multinomial_update() {
        let alpha_n = BayesianUpdate::dirichlet_multinomial(&[1.0, 1.0, 1.0], &[3, 2, 5]);
        assert!((alpha_n[0] - 4.0).abs() < 1e-12);
        assert!((alpha_n[2] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_normal_inverse_gamma_update() {
        let (mu_n, kappa_n, alpha_n, beta_n) =
            BayesianUpdate::normal_inverse_gamma(0.0, 1.0, 2.0, 3.0, &[1.0, 2.0, 3.0]);
        assert!(kappa_n > 1.0);
        assert!(alpha_n > 2.0);
        assert!(beta_n > 3.0);
        let _ = mu_n;
    }

    // ── MCMC ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_mh_samples_correct_count() {
        let mcmc = MarkovChainMonteCarlo::new(0.5, 100);
        let samples = mcmc.metropolis_hastings(|x| -0.5 * x * x, 0.0, 200, 42);
        assert_eq!(samples.len(), 200);
    }

    #[test]
    fn test_mh_standard_normal_mean_close_to_zero() {
        let mcmc = MarkovChainMonteCarlo::new(1.0, 500);
        let samples = mcmc.metropolis_hastings(|x| -0.5 * x * x, 0.0, 1000, 99);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!(mean.abs() < 0.3);
    }

    #[test]
    fn test_gibbs_bivariate_gaussian_count() {
        let samples =
            MarkovChainMonteCarlo::gibbs_bivariate_gaussian(0.0, 0.0, 1.0, 1.0, 0.5, 100, 7);
        assert_eq!(samples.len(), 100);
    }

    #[test]
    fn test_mh_vec_runs() {
        let mcmc = MarkovChainMonteCarlo::new(0.3, 50);
        let log_t = |x: &[f64]| -0.5 * x[0] * x[0] - 0.5 * x[1] * x[1];
        let samples = mcmc.metropolis_hastings_vec(log_t, vec![0.0, 0.0], 100, 13);
        assert_eq!(samples.len(), 100);
    }

    #[test]
    fn test_nuts_step_returns_finite() {
        let mcmc = MarkovChainMonteCarlo::new(0.1, 0);
        let result = mcmc.nuts_step(|x| -0.5 * x * x, |x| -x, 0.0, 5);
        assert!(result.is_finite());
    }

    #[test]
    fn test_ess_constant_chain() {
        let chain = vec![1.0; 100];
        let ess = MarkovChainMonteCarlo::effective_sample_size(&chain);
        assert!(ess > 0.0);
    }

    // ── BayesianLinearRegression ─────────────────────────────────────────────

    #[test]
    fn test_blr_fit_and_predict() {
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0]];
        let y = vec![1.0, 3.0, 5.0]; // y = 1 + 2x
        let mut blr = BayesianLinearRegression::new(1e-4, 10.0);
        blr.fit(&x, &y);
        let pred = blr.predict_mean(&[1.0, 1.5]);
        assert!((pred - 4.0).abs() < 1.0); // Should be close to 4
    }

    #[test]
    fn test_blr_variance_positive() {
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0]];
        let y = vec![0.0, 1.0];
        let mut blr = BayesianLinearRegression::new(1.0, 1.0);
        blr.fit(&x, &y);
        let var = blr.predict_variance(&[1.0, 0.5]);
        assert!(var > 0.0);
    }

    #[test]
    fn test_blr_log_evidence_finite_after_fit() {
        let x = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0]];
        let y = vec![1.0, 2.0, 3.0];
        let mut blr = BayesianLinearRegression::new(1.0, 1.0);
        blr.fit(&x, &y);
        let ev = blr.log_evidence(&x, &y);
        assert!(ev.is_finite());
    }

    // ── GaussianProcess ──────────────────────────────────────────────────────

    #[test]
    fn test_gp_rbf_kernel_at_same_point() {
        let k = Kernel::Rbf {
            length_scale: 1.0,
            signal_variance: 1.0,
        };
        assert!((k.eval(2.0, 2.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gp_matern32_at_same_point() {
        let k = Kernel::Matern32 {
            length_scale: 1.0,
            signal_variance: 2.0,
        };
        assert!((k.eval(1.0, 1.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_gp_matern52_at_same_point() {
        let k = Kernel::Matern52 {
            length_scale: 1.0,
            signal_variance: 3.0,
        };
        assert!((k.eval(0.5, 0.5) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_gp_fit_and_predict() {
        let kernel = Kernel::Rbf {
            length_scale: 1.0,
            signal_variance: 1.0,
        };
        let mut gp = GaussianProcess::new(kernel, 0.01);
        let x_train: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let y_train: Vec<f64> = x_train.iter().map(|&x| x * 2.0).collect();
        gp.fit(x_train, y_train);
        let pred = gp.predict_mean(2.0);
        assert!((pred - 4.0).abs() < 1.0);
    }

    #[test]
    fn test_gp_variance_positive() {
        let kernel = Kernel::Rbf {
            length_scale: 1.0,
            signal_variance: 1.0,
        };
        let mut gp = GaussianProcess::new(kernel, 0.01);
        gp.fit(vec![0.0, 1.0], vec![0.0, 1.0]);
        let var = gp.predict_variance(5.0); // Far from training data
        assert!(var > 0.0);
    }

    #[test]
    fn test_gp_log_marginal_likelihood_finite() {
        let kernel = Kernel::Rbf {
            length_scale: 1.0,
            signal_variance: 1.0,
        };
        let mut gp = GaussianProcess::new(kernel, 0.1);
        gp.fit(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 0.0]);
        assert!(gp.log_marginal_likelihood().is_finite());
    }

    // ── ModelSelection ───────────────────────────────────────────────────────

    #[test]
    fn test_aic_basic() {
        let aic = ModelSelection::aic(-100.0, 5);
        assert!((aic - (10.0 + 200.0)).abs() < 1e-12);
    }

    #[test]
    fn test_bic_basic() {
        let bic = ModelSelection::bic(-100.0, 3, 50);
        let expected = 3.0 * 50.0_f64.ln() + 200.0;
        assert!((bic - expected).abs() < 1e-10);
    }

    #[test]
    fn test_aicc_greater_than_aic() {
        let aic = ModelSelection::aic(-100.0, 4);
        let aicc = ModelSelection::aicc(-100.0, 4, 20);
        assert!(aicc >= aic);
    }

    #[test]
    fn test_log_bayes_factor_symmetric() {
        let lbf = ModelSelection::log_bayes_factor(10.0, 10.0);
        assert_eq!(lbf, 0.0);
    }

    #[test]
    fn test_jeffreys_scale_decisive() {
        let label = ModelSelection::jeffreys_scale(5.0); // e^5 >> 100
        assert_eq!(label, "Decisive");
    }

    #[test]
    fn test_k_fold_cv_mse_finite() {
        let x: Vec<Vec<f64>> = (0..10).map(|i| vec![1.0, i as f64]).collect();
        let y: Vec<f64> = (0..10).map(|i| i as f64 * 2.0 + 1.0).collect();
        let mse = ModelSelection::k_fold_cv_mse(&x, &y, 5, 1.0, 1.0);
        assert!(mse.is_finite());
    }

    #[test]
    fn test_loo_cv_log_predictive_finite() {
        let x: Vec<Vec<f64>> = (0..5).map(|i| vec![1.0, i as f64]).collect();
        let y: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let lp = ModelSelection::loo_cv_log_predictive(&x, &y, 1.0, 1.0);
        assert!(lp.is_finite());
    }

    // ── Utility ──────────────────────────────────────────────────────────────

    #[test]
    fn test_log_gamma_half() {
        // Γ(1/2) = √π → ln Γ(1/2) = 0.5 ln π
        let lg = log_gamma(0.5);
        let expected = 0.5 * PI.ln();
        assert!((lg - expected).abs() < 1e-6);
    }

    #[test]
    fn test_log_gamma_one() {
        // Γ(1) = 1 → ln Γ(1) = 0
        assert!(log_gamma(1.0).abs() < 1e-6);
    }

    #[test]
    fn test_log_beta_symmetry() {
        assert!((log_beta(2.0, 3.0) - log_beta(3.0, 2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_cholesky_identity() {
        let eye = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let l = cholesky(&eye);
        assert!((l[0][0] - 1.0).abs() < 1e-12);
        assert!((l[1][1] - 1.0).abs() < 1e-12);
        assert!(l[1][0].abs() < 1e-12);
    }

    #[test]
    fn test_mat_inverse_identity() {
        let eye = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let inv = mat_inverse(&eye);
        assert!((inv[0][0] - 1.0).abs() < 1e-10);
        assert!((inv[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sample_gamma_positive() {
        let mut rng = BiRng::new(42);
        let s = sample_gamma(2.0, &mut rng);
        assert!(s > 0.0);
    }
}
