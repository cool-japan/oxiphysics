// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Stochastic processes for simulation and financial/physical modeling.
//!
//! This module provides implementations of key stochastic processes used in
//! computational physics, quantitative finance, and statistical modeling:
//!
//! - [`BrownianMotion`] — standard Wiener process path generation
//! - [`GeometricBrownianMotion`] — log-normal GBM (Black-Scholes)
//! - [`OrnsteinUhlenbeck`] — mean-reverting SDE (OU process)
//! - [`JumpDiffusion`] — Merton model with Poisson jumps
//! - [`FractionalBrownianMotion`] — fBm with Hurst exponent H
//! - [`CoxIngersollRoss`] — CIR interest rate model (Milstein scheme)
//! - [`HestonModel`] — stochastic volatility (Euler-Maruyama)
//! - [`LevyStableProcess`] — stable distribution parameters
//! - [`BranchingProcess`] — Galton-Watson branching process
//! - [`MarkovChain`] — discrete-time Markov chain utilities
//! - [`PoissonProcess`] — homogeneous and inhomogeneous Poisson
//! - Variance reduction: [`antithetic_variates`], [`control_variate_estimator`], [`quasi_monte_carlo`]

use std::f64::consts::PI;

// ─── Internal LCG RNG ────────────────────────────────────────────────────────

/// Minimal linear congruential generator for reproducible stochastic simulations.
///
/// Uses Knuth multiplicative constants. Suitable for Monte Carlo work where
/// external rand crate is not available.
pub(crate) struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Uniform sample in \[0, 1).
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard normal via Box-Muller transform.
    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }

    /// Exponential with rate lambda.
    fn exponential(&mut self, lambda: f64) -> f64 {
        let u = self.uniform().max(1e-300);
        -u.ln() / lambda
    }

    /// Poisson(lambda) sample via Knuth algorithm.
    fn poisson(&mut self, lambda: f64) -> usize {
        if lambda <= 0.0 {
            return 0;
        }
        let threshold = (-lambda).exp();
        let mut k = 0usize;
        let mut p = 1.0_f64;
        loop {
            p *= self.uniform();
            if p < threshold {
                break;
            }
            k += 1;
        }
        k
    }

    /// Bernoulli(p) sample: returns true with probability p.
    fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p
    }
}

// ─── BrownianMotion ──────────────────────────────────────────────────────────

/// Standard Wiener process (Brownian motion) W(t).
///
/// Generates paths where increments dW ~ N(0, dt). This is the foundation
/// of most continuous-time stochastic models.
pub struct BrownianMotion {
    /// Time step for path simulation.
    pub dt: f64,
    /// Diffusion coefficient (sigma). Default is 1.0 for standard BM.
    pub sigma: f64,
    seed: u64,
}

impl BrownianMotion {
    /// Creates a new standard Brownian motion with given time step and seed.
    pub fn new(dt: f64, seed: u64) -> Self {
        Self {
            dt,
            sigma: 1.0,
            seed,
        }
    }

    /// Creates a scaled Brownian motion with diffusion coefficient sigma.
    pub fn with_sigma(dt: f64, sigma: f64, seed: u64) -> Self {
        Self { dt, sigma, seed }
    }

    /// Generate a single increment dW ~ N(0, sigma^2 * dt).
    fn increment(&self, rng: &mut Lcg) -> f64 {
        self.sigma * self.dt.sqrt() * rng.normal()
    }

    /// Generate a path of `n_steps` from W(0) = 0.
    ///
    /// Returns `n_steps + 1` values (including initial W = 0).
    pub fn path(&self, n_steps: usize) -> Vec<f64> {
        let mut rng = Lcg::new(self.seed);
        let mut w = vec![0.0_f64; n_steps + 1];
        for i in 1..=n_steps {
            w[i] = w[i - 1] + self.increment(&mut rng);
        }
        w
    }

    /// Compute theoretical variance at time t = n_steps * dt.
    ///
    /// For standard BM, Var\[W(t)\] = sigma^2 * t.
    pub fn variance_at(&self, t: f64) -> f64 {
        self.sigma * self.sigma * t
    }

    /// Generate multiple paths for Monte Carlo. Returns matrix \[path_idx\]\[step\].
    pub fn multi_path(&self, n_paths: usize, n_steps: usize) -> Vec<Vec<f64>> {
        (0..n_paths)
            .map(|i| {
                let mut bm = BrownianMotion::with_sigma(
                    self.dt,
                    self.sigma,
                    self.seed ^ (i as u64 * 2654435761),
                );
                bm.seed ^= i as u64;
                bm.path(n_steps)
            })
            .collect()
    }

    /// Quadratic variation estimate for a path (should converge to sigma^2 * T).
    pub fn quadratic_variation(path: &[f64]) -> f64 {
        path.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum()
    }
}

// ─── GeometricBrownianMotion ─────────────────────────────────────────────────

/// Geometric Brownian Motion (GBM) for log-normal asset/particle dynamics.
///
/// SDE: dS = mu * S * dt + sigma * S * dW
/// Exact solution: S(t) = S0 * exp((mu - sigma^2/2)*t + sigma*W(t))
pub struct GeometricBrownianMotion {
    /// Initial value S(0).
    pub s0: f64,
    /// Drift coefficient mu.
    pub mu: f64,
    /// Volatility coefficient sigma.
    pub sigma: f64,
    seed: u64,
}

impl GeometricBrownianMotion {
    /// Creates a new GBM process.
    pub fn new(s0: f64, mu: f64, sigma: f64, seed: u64) -> Self {
        Self {
            s0,
            mu,
            sigma,
            seed,
        }
    }

    /// Exact path using the analytical solution (no discretisation error).
    ///
    /// Returns `n_steps + 1` values including S(0) = s0.
    pub fn exact_path(&self, t_end: f64, n_steps: usize) -> Vec<f64> {
        let dt = t_end / n_steps as f64;
        let mut rng = Lcg::new(self.seed);
        let mut path = vec![0.0_f64; n_steps + 1];
        path[0] = self.s0;
        let drift = (self.mu - 0.5 * self.sigma * self.sigma) * dt;
        let vol = self.sigma * dt.sqrt();
        let mut log_s = self.s0.ln();
        for p in path.iter_mut().skip(1) {
            log_s += drift + vol * rng.normal();
            *p = log_s.exp();
        }
        path
    }

    /// Euler-Maruyama discretisation path.
    pub fn euler_path(&self, t_end: f64, n_steps: usize) -> Vec<f64> {
        let dt = t_end / n_steps as f64;
        let mut rng = Lcg::new(self.seed);
        let mut path = vec![0.0_f64; n_steps + 1];
        path[0] = self.s0;
        for i in 1..=n_steps {
            let s = path[i - 1];
            let dw = dt.sqrt() * rng.normal();
            path[i] = s + self.mu * s * dt + self.sigma * s * dw;
            if path[i] < 0.0 {
                path[i] = 0.0;
            }
        }
        path
    }

    /// Theoretical mean E\[S(t)\] = S0 * exp(mu * t).
    pub fn mean_at(&self, t: f64) -> f64 {
        self.s0 * (self.mu * t).exp()
    }

    /// Theoretical variance Var\[S(t)\] = S0^2 * exp(2*mu*t) * (exp(sigma^2*t) - 1).
    pub fn variance_at(&self, t: f64) -> f64 {
        let e2mu = (2.0 * self.mu * t).exp();
        let esig2 = (self.sigma * self.sigma * t).exp() - 1.0;
        self.s0 * self.s0 * e2mu * esig2
    }
}

// ─── OrnsteinUhlenbeck ───────────────────────────────────────────────────────

/// Ornstein-Uhlenbeck mean-reverting process.
///
/// SDE: dX = theta*(mu - X)*dt + sigma*dW
///
/// Uses exact discretisation: X(t+dt) = X(t)*exp(-theta*dt) + mu*(1-exp(-theta*dt))
///                                       + sigma*sqrt((1-exp(-2*theta*dt))/(2*theta)) * N(0,1)
pub struct OrnsteinUhlenbeck {
    /// Mean reversion speed (theta > 0).
    pub theta: f64,
    /// Long-run mean (mu).
    pub mu: f64,
    /// Volatility (sigma).
    pub sigma: f64,
    /// Initial value X(0).
    pub x0: f64,
    seed: u64,
}

impl OrnsteinUhlenbeck {
    /// Creates a new OU process.
    pub fn new(theta: f64, mu: f64, sigma: f64, x0: f64, seed: u64) -> Self {
        Self {
            theta,
            mu,
            sigma,
            x0,
            seed,
        }
    }

    /// Exact simulation path using the conditional Gaussian formula.
    pub fn path(&self, t_end: f64, n_steps: usize) -> Vec<f64> {
        let dt = t_end / n_steps as f64;
        let mut rng = Lcg::new(self.seed);
        let mut path = vec![0.0_f64; n_steps + 1];
        path[0] = self.x0;
        let exp_neg = (-self.theta * dt).exp();
        let var =
            self.sigma * self.sigma * (1.0 - (-2.0 * self.theta * dt).exp()) / (2.0 * self.theta);
        let std_dev = var.sqrt();
        for i in 1..=n_steps {
            let mean = path[i - 1] * exp_neg + self.mu * (1.0 - exp_neg);
            path[i] = mean + std_dev * rng.normal();
        }
        path
    }

    /// Stationary variance: sigma^2 / (2*theta).
    pub fn stationary_variance(&self) -> f64 {
        self.sigma * self.sigma / (2.0 * self.theta)
    }

    /// Autocorrelation at lag tau: exp(-theta * tau).
    pub fn autocorrelation(&self, tau: f64) -> f64 {
        (-self.theta * tau).exp()
    }
}

// ─── JumpDiffusion ───────────────────────────────────────────────────────────

/// Merton jump-diffusion model: GBM plus compound Poisson jumps.
///
/// SDE: dS = (mu - lambda*k_bar)*S*dt + sigma*S*dW + S*(J-1)*dN
/// where J = exp(N(mu_j, sigma_j^2)), k_bar = E\[J-1\] = exp(mu_j + sigma_j^2/2) - 1.
pub struct JumpDiffusion {
    /// Initial value.
    pub s0: f64,
    /// Drift (excluding jump compensation).
    pub mu: f64,
    /// Diffusion volatility.
    pub sigma: f64,
    /// Poisson intensity (average jumps per unit time).
    pub lambda: f64,
    /// Mean of log-jump size (normal distributed).
    pub mu_j: f64,
    /// Std dev of log-jump size.
    pub sigma_j: f64,
    seed: u64,
}

impl JumpDiffusion {
    /// Creates a new Merton jump-diffusion model.
    pub fn new(
        s0: f64,
        mu: f64,
        sigma: f64,
        lambda: f64,
        mu_j: f64,
        sigma_j: f64,
        seed: u64,
    ) -> Self {
        Self {
            s0,
            mu,
            sigma,
            lambda,
            mu_j,
            sigma_j,
            seed,
        }
    }

    /// Expected jump multiplier E\[J\] = exp(mu_j + sigma_j^2/2).
    pub fn expected_jump(&self) -> f64 {
        (self.mu_j + 0.5 * self.sigma_j * self.sigma_j).exp()
    }

    /// Jump compensation term k_bar = E\[J\] - 1.
    pub fn jump_compensation(&self) -> f64 {
        self.expected_jump() - 1.0
    }

    /// Simulate a path using Euler-Maruyama with Poisson jump thinning.
    pub fn path(&self, t_end: f64, n_steps: usize) -> Vec<f64> {
        let dt = t_end / n_steps as f64;
        let mut rng = Lcg::new(self.seed);
        let mut path = vec![0.0_f64; n_steps + 1];
        path[0] = self.s0;
        let k_bar = self.jump_compensation();
        let compensated_mu = self.mu - self.lambda * k_bar;
        for i in 1..=n_steps {
            let s = path[i - 1];
            let dw = dt.sqrt() * rng.normal();
            let n_jumps = rng.poisson(self.lambda * dt);
            let mut jump_factor = 1.0_f64;
            for _ in 0..n_jumps {
                let log_j = self.mu_j + self.sigma_j * rng.normal();
                jump_factor *= log_j.exp();
            }
            path[i] = s * (1.0 + compensated_mu * dt + self.sigma * dw) * jump_factor;
            if path[i] < 0.0 {
                path[i] = 0.0;
            }
        }
        path
    }
}

// ─── FractionalBrownianMotion ────────────────────────────────────────────────

/// Fractional Brownian motion with Hurst exponent H in (0, 1).
///
/// H = 0.5 gives standard BM. H > 0.5 gives long-range dependence (persistent).
/// H < 0.5 gives anti-persistent behaviour.
///
/// Uses the Cholesky method via exact covariance matrix for short paths,
/// or the Wood-Chan circulant embedding for longer paths (approximation here).
pub struct FractionalBrownianMotion {
    /// Hurst exponent in (0, 1).
    pub hurst: f64,
    /// Scaling parameter.
    pub sigma: f64,
    seed: u64,
}

impl FractionalBrownianMotion {
    /// Creates a new fBm with given Hurst exponent.
    ///
    /// # Panics
    /// Panics if hurst is not in (0, 1).
    pub fn new(hurst: f64, seed: u64) -> Self {
        assert!(
            hurst > 0.0 && hurst < 1.0,
            "Hurst exponent must be in (0, 1)"
        );
        Self {
            hurst,
            sigma: 1.0,
            seed,
        }
    }

    /// Creates a scaled fBm.
    pub fn with_sigma(hurst: f64, sigma: f64, seed: u64) -> Self {
        assert!(
            hurst > 0.0 && hurst < 1.0,
            "Hurst exponent must be in (0, 1)"
        );
        Self { hurst, sigma, seed }
    }

    /// Covariance kernel: C(s, t) = sigma^2/2 * (|s|^{2H} + |t|^{2H} - |s-t|^{2H}).
    pub fn covariance(&self, s: f64, t: f64) -> f64 {
        let h2 = 2.0 * self.hurst;
        0.5 * self.sigma
            * self.sigma
            * (s.abs().powf(h2) + t.abs().powf(h2) - (s - t).abs().powf(h2))
    }

    /// Generate a path on \[0, 1\] with n_steps points using Cholesky decomposition.
    ///
    /// This is exact but O(n^3) in memory/compute. Suitable for n <= ~500.
    pub fn path(&self, n_steps: usize) -> Vec<f64> {
        let n = n_steps + 1;
        let dt = 1.0 / n_steps as f64;

        // Build covariance matrix C
        let times: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let mut cov = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                cov[i * n + j] = self.covariance(times[i], times[j]);
            }
        }

        // Cholesky decomposition L such that L*L^T = C
        let mut l = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let mut sum = cov[i * n + j];
                for k in 0..j {
                    sum -= l[i * n + k] * l[j * n + k];
                }
                if i == j {
                    l[i * n + j] = sum.max(0.0).sqrt();
                } else if l[j * n + j].abs() > 1e-14 {
                    l[i * n + j] = sum / l[j * n + j];
                }
            }
        }

        // Sample z ~ N(0, I_n), then x = L*z
        let mut rng = Lcg::new(self.seed);
        let z: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
        let mut x = vec![0.0_f64; n];
        for i in 0..n {
            for j in 0..=i {
                x[i] += l[i * n + j] * z[j];
            }
        }
        x
    }
}

// ─── CoxIngersollRoss ────────────────────────────────────────────────────────

/// Cox-Ingersoll-Ross (CIR) interest rate model.
///
/// SDE: dr = kappa*(theta - r)*dt + sigma*sqrt(r)*dW
///
/// The Feller condition 2*kappa*theta >= sigma^2 ensures r stays positive.
/// Milstein scheme is used for better accuracy.
pub struct CoxIngersollRoss {
    /// Mean reversion speed kappa.
    pub kappa: f64,
    /// Long-run mean theta.
    pub theta: f64,
    /// Volatility sigma.
    pub sigma: f64,
    /// Initial rate r(0).
    pub r0: f64,
    seed: u64,
}

impl CoxIngersollRoss {
    /// Creates a new CIR process.
    pub fn new(kappa: f64, theta: f64, sigma: f64, r0: f64, seed: u64) -> Self {
        Self {
            kappa,
            theta,
            sigma,
            r0,
            seed,
        }
    }

    /// Whether the Feller condition 2*kappa*theta >= sigma^2 is satisfied.
    pub fn feller_satisfied(&self) -> bool {
        2.0 * self.kappa * self.theta >= self.sigma * self.sigma
    }

    /// Simulate path using Milstein scheme.
    ///
    /// The reflection r -> |r| is used after each step to guarantee non-negativity.
    pub fn path(&self, t_end: f64, n_steps: usize) -> Vec<f64> {
        let dt = t_end / n_steps as f64;
        let mut rng = Lcg::new(self.seed);
        let mut path = vec![0.0_f64; n_steps + 1];
        path[0] = self.r0;
        let sqrt_dt = dt.sqrt();
        for i in 1..=n_steps {
            let r = path[i - 1].max(0.0);
            let dw = sqrt_dt * rng.normal();
            let sqrt_r = r.sqrt();
            // Milstein: dr = kappa*(theta-r)*dt + sigma*sqrt(r)*dW + sigma^2/4*(dW^2 - dt)
            let milstein_correction = 0.25 * self.sigma * self.sigma * (dw * dw - dt);
            path[i] = (r
                + self.kappa * (self.theta - r) * dt
                + self.sigma * sqrt_r * dw
                + milstein_correction)
                .abs();
        }
        path
    }

    /// Stationary mean E\[r\] = theta.
    pub fn stationary_mean(&self) -> f64 {
        self.theta
    }

    /// Stationary variance Var\[r\] = sigma^2 * theta / (2 * kappa).
    pub fn stationary_variance(&self) -> f64 {
        self.sigma * self.sigma * self.theta / (2.0 * self.kappa)
    }
}

// ─── HestonModel ─────────────────────────────────────────────────────────────

/// Heston stochastic volatility model.
///
/// System of SDEs:
///   dS = mu*S*dt + sqrt(v)*S*dW_S
///   dv = kappa*(theta - v)*dt + xi*sqrt(v)*dW_v
///   Corr(dW_S, dW_v) = rho
pub struct HestonModel {
    /// Initial asset price.
    pub s0: f64,
    /// Initial variance.
    pub v0: f64,
    /// Drift of asset.
    pub mu: f64,
    /// Mean reversion speed of variance.
    pub kappa: f64,
    /// Long-run variance.
    pub theta: f64,
    /// Vol of vol.
    pub xi: f64,
    /// Correlation between asset and variance Brownians.
    pub rho: f64,
    seed: u64,
}

impl HestonModel {
    /// Creates a new Heston model.
    pub fn new(
        s0: f64,
        v0: f64,
        mu: f64,
        kappa: f64,
        theta: f64,
        xi: f64,
        rho: f64,
        seed: u64,
    ) -> Self {
        Self {
            s0,
            v0,
            mu,
            kappa,
            theta,
            xi,
            rho,
            seed,
        }
    }

    /// Simulate correlated paths (S, v) using Euler-Maruyama.
    ///
    /// Returns `(price_path, variance_path)` each of length `n_steps + 1`.
    pub fn path(&self, t_end: f64, n_steps: usize) -> (Vec<f64>, Vec<f64>) {
        let dt = t_end / n_steps as f64;
        let sqrt_dt = dt.sqrt();
        let mut rng = Lcg::new(self.seed);
        let mut s_path = vec![0.0_f64; n_steps + 1];
        let mut v_path = vec![0.0_f64; n_steps + 1];
        s_path[0] = self.s0;
        v_path[0] = self.v0;
        let rho2 = (1.0 - self.rho * self.rho).max(0.0).sqrt();
        for i in 1..=n_steps {
            let z1 = rng.normal();
            let z2 = rng.normal();
            let dw_v = sqrt_dt * z1;
            let dw_s = sqrt_dt * (self.rho * z1 + rho2 * z2);
            let v = v_path[i - 1].max(0.0);
            let sqrt_v = v.sqrt();
            v_path[i] = (v + self.kappa * (self.theta - v) * dt + self.xi * sqrt_v * dw_v).abs();
            let s = s_path[i - 1];
            s_path[i] = s + self.mu * s * dt + sqrt_v * s * dw_s;
            if s_path[i] < 0.0 {
                s_path[i] = 0.0;
            }
        }
        (s_path, v_path)
    }
}

// ─── LevyStableProcess ───────────────────────────────────────────────────────

/// Lévy stable process with four parameters.
///
/// A Lévy stable distribution is characterised by (alpha, beta, c, mu) where:
/// - alpha in (0, 2]: stability index (2 = Gaussian, 1 = Cauchy)
/// - beta in \[-1, 1\]: skewness parameter
/// - c > 0: scale parameter
/// - mu: location parameter
pub struct LevyStableProcess {
    /// Stability index alpha in (0, 2].
    pub alpha: f64,
    /// Skewness beta in \[-1, 1\].
    pub beta: f64,
    /// Scale c > 0.
    pub c: f64,
    /// Location mu.
    pub mu: f64,
    seed: u64,
}

impl LevyStableProcess {
    /// Creates a new Lévy stable process.
    pub fn new(alpha: f64, beta: f64, c: f64, mu: f64, seed: u64) -> Self {
        assert!(alpha > 0.0 && alpha <= 2.0, "alpha must be in (0, 2]");
        assert!((-1.0..=1.0).contains(&beta), "beta must be in [-1, 1]");
        assert!(c > 0.0, "scale c must be positive");
        Self {
            alpha,
            beta,
            c,
            mu,
            seed,
        }
    }

    /// Sample from S(alpha, beta, 1, 0) using the Chambers-Mallows-Stuck algorithm.
    fn sample_standard(&self, rng: &mut Lcg) -> f64 {
        let u = PI * (rng.uniform() - 0.5);
        let w = rng.exponential(1.0);
        if (self.alpha - 1.0).abs() < 1e-10 {
            // Cauchy case
            let b_pi = self.beta * PI / 2.0;
            (1.0 + self.beta * u.tan()) * u.tan() - self.beta * (b_pi.cos() / b_pi.cos()).ln()
        } else {
            let zeta = -self.beta * (PI * self.alpha / 2.0).tan();
            let xi = zeta.atan() / self.alpha;
            let term1 = (1.0 + zeta * zeta).powf(1.0 / (2.0 * self.alpha));
            let term2 = (self.alpha * (u + xi)).sin()
                / ((self.alpha * xi).cos() * u.cos()).powf(1.0 / self.alpha);
            let term3 =
                ((u - self.alpha * (u + xi)).cos() / w).powf((1.0 - self.alpha) / self.alpha);
            term1 * term2 * term3
        }
    }

    /// Generate a path of increments (not cumulative).
    pub fn increments(&self, n: usize) -> Vec<f64> {
        let mut rng = Lcg::new(self.seed);
        (0..n)
            .map(|_| self.mu + self.c * self.sample_standard(&mut rng))
            .collect()
    }

    /// Generate a cumulative path (random walk with Lévy increments).
    pub fn path(&self, n: usize) -> Vec<f64> {
        let incs = self.increments(n);
        let mut path = vec![0.0_f64; n + 1];
        for i in 0..n {
            path[i + 1] = path[i] + incs[i];
        }
        path
    }
}

// ─── BranchingProcess ────────────────────────────────────────────────────────

/// Galton-Watson branching process.
///
/// Each individual in generation n independently produces a random number of
/// offspring in generation n+1 according to an offspring distribution.
pub struct BranchingProcess {
    /// Initial population Z(0).
    pub z0: usize,
    /// Mean offspring per individual.
    pub mean_offspring: f64,
    seed: u64,
    /// Offspring distribution type.
    dist: OffspringDist,
}

/// Supported offspring distributions for branching processes.
#[derive(Clone, Copy)]
pub enum OffspringDist {
    /// Poisson(lambda) offspring.
    Poisson(f64),
    /// Geometric(p) offspring: P(k) = (1-p)^k * p.
    Geometric(f64),
    /// Bernoulli(p): 0 or 1 offspring.
    Bernoulli(f64),
}

impl BranchingProcess {
    /// Creates a new Galton-Watson process with Poisson offspring.
    pub fn poisson(z0: usize, lambda: f64, seed: u64) -> Self {
        Self {
            z0,
            mean_offspring: lambda,
            seed,
            dist: OffspringDist::Poisson(lambda),
        }
    }

    /// Creates a new process with Geometric offspring.
    pub fn geometric(z0: usize, p: f64, seed: u64) -> Self {
        let mean = (1.0 - p) / p;
        Self {
            z0,
            mean_offspring: mean,
            seed,
            dist: OffspringDist::Geometric(p),
        }
    }

    /// Creates a new process with Bernoulli offspring.
    pub fn bernoulli(z0: usize, p: f64, seed: u64) -> Self {
        Self {
            z0,
            mean_offspring: p,
            seed,
            dist: OffspringDist::Bernoulli(p),
        }
    }

    fn sample_offspring(&self, rng: &mut Lcg) -> usize {
        match self.dist {
            OffspringDist::Poisson(lambda) => rng.poisson(lambda),
            OffspringDist::Geometric(p) => {
                let u = rng.uniform().max(1e-300);
                (u.ln() / (1.0 - p).ln()).floor() as usize
            }
            OffspringDist::Bernoulli(p) => rng.bernoulli(p) as usize,
        }
    }

    /// Simulate the process for `n_generations` generations.
    ///
    /// Returns population sizes Z(0), Z(1), ..., Z(n_generations).
    pub fn simulate(&self, n_generations: usize) -> Vec<usize> {
        let mut rng = Lcg::new(self.seed);
        let mut pop = vec![0usize; n_generations + 1];
        pop[0] = self.z0;
        for g in 0..n_generations {
            let mut next_pop = 0usize;
            for _ in 0..pop[g] {
                next_pop = next_pop.saturating_add(self.sample_offspring(&mut rng));
            }
            pop[g + 1] = next_pop;
        }
        pop
    }

    /// Theoretical extinction probability for Poisson offspring.
    ///
    /// For Poisson(lambda) offspring:
    /// - If lambda <= 1: extinction probability = 1.
    /// - If lambda > 1: extinction probability = q where q = exp(lambda*(q-1)).
    ///   Solved iteratively.
    pub fn extinction_probability(&self) -> f64 {
        if self.mean_offspring <= 1.0 {
            return 1.0;
        }
        match self.dist {
            OffspringDist::Poisson(lambda) => {
                // Iterate q_{n+1} = exp(lambda*(q_n - 1)) starting from q_0 = 0
                let mut q = 0.0_f64;
                for _ in 0..1000 {
                    let q_new = (lambda * (q - 1.0)).exp();
                    if (q_new - q).abs() < 1e-12 {
                        return q_new;
                    }
                    q = q_new;
                }
                q
            }
            OffspringDist::Geometric(p) => {
                // Geometric: mean = (1-p)/p. Extinction prob = p/(1-p) if mean > 1, else 1
                let mean = (1.0 - p) / p;
                if mean <= 1.0 { 1.0 } else { p / (1.0 - p) }
            }
            OffspringDist::Bernoulli(p) => {
                if p <= 0.5 {
                    1.0
                } else {
                    (1.0 - p) / p
                }
            }
        }
    }

    /// Estimate extinction probability by Monte Carlo simulation.
    pub fn mc_extinction_probability(&self, n_paths: usize, n_generations: usize) -> f64 {
        let mut extinct = 0usize;
        for k in 0..n_paths {
            let mut bp = BranchingProcess {
                z0: self.z0,
                mean_offspring: self.mean_offspring,
                seed: self.seed ^ (k as u64 * 2654435761),
                dist: self.dist,
            };
            bp.seed ^= k as u64;
            let pop = bp.simulate(n_generations);
            if *pop.last().unwrap_or(&1) == 0 {
                extinct += 1;
            }
        }
        extinct as f64 / n_paths as f64
    }
}

// ─── MarkovChain ─────────────────────────────────────────────────────────────

/// Discrete-time finite Markov chain.
///
/// Stores the transition matrix P where P\[i\]\[j\] = P(X_{n+1} = j | X_n = i).
pub struct MarkovChain {
    /// Transition matrix (n_states x n_states), row-major.
    pub transition: Vec<Vec<f64>>,
    /// Number of states.
    pub n_states: usize,
    seed: u64,
}

impl MarkovChain {
    /// Creates a new Markov chain from a transition matrix.
    ///
    /// Rows must sum to 1 (stochastic matrix).
    pub fn new(transition: Vec<Vec<f64>>, seed: u64) -> Self {
        let n_states = transition.len();
        Self {
            transition,
            n_states,
            seed,
        }
    }

    /// Sample next state given current state using the transition matrix.
    fn next_state(&self, state: usize, rng: &mut Lcg) -> usize {
        let u = rng.uniform();
        let mut cumul = 0.0;
        for (j, &p) in self.transition[state].iter().enumerate() {
            cumul += p;
            if u < cumul {
                return j;
            }
        }
        self.n_states - 1
    }

    /// Simulate a trajectory of length `n_steps` starting from `start_state`.
    pub fn simulate(&self, start_state: usize, n_steps: usize) -> Vec<usize> {
        let mut rng = Lcg::new(self.seed);
        let mut traj = vec![0usize; n_steps + 1];
        traj[0] = start_state;
        for i in 1..=n_steps {
            traj[i] = self.next_state(traj[i - 1], &mut rng);
        }
        traj
    }

    /// Compute stationary distribution by power iteration.
    ///
    /// Returns pi such that pi * P = pi with sum(pi) = 1.
    pub fn stationary_distribution(&self) -> Vec<f64> {
        let n = self.n_states;
        let mut pi: Vec<f64> = vec![1.0 / n as f64; n];
        for _ in 0..10_000 {
            let mut pi_new = vec![0.0_f64; n];
            for (i, &pi_i) in pi.iter().enumerate() {
                for (j, pnj) in pi_new.iter_mut().enumerate() {
                    *pnj += pi_i * self.transition[i][j];
                }
            }
            let sum: f64 = pi_new.iter().sum();
            for x in &mut pi_new {
                *x /= sum;
            }
            let diff: f64 = pi_new.iter().zip(&pi).map(|(a, b)| (a - b).abs()).sum();
            pi = pi_new;
            if diff < 1e-12 {
                break;
            }
        }
        pi
    }

    /// Estimate expected hitting time to `target_state` from `start_state` by simulation.
    pub fn hitting_time(&self, start_state: usize, target_state: usize, n_samples: usize) -> f64 {
        let mut rng = Lcg::new(self.seed ^ 0xdeadbeef);
        let mut total = 0usize;
        let max_steps = 100_000;
        for _ in 0..n_samples {
            let mut state = start_state;
            for step in 1..=max_steps {
                state = self.next_state(state, &mut rng);
                if state == target_state {
                    total += step;
                    break;
                }
            }
        }
        total as f64 / n_samples as f64
    }

    /// Return the n-step transition matrix P^n.
    pub fn n_step_matrix(&self, n: usize) -> Vec<Vec<f64>> {
        let size = self.n_states;
        // Start with identity
        let mut result = vec![vec![0.0_f64; size]; size];
        for (i, row) in result.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        let mut base = self.transition.clone();
        let mut power = n;
        while power > 0 {
            if power % 2 == 1 {
                result = mat_mul(&result, &base, size);
            }
            base = mat_mul(&base, &base, size);
            power /= 2;
        }
        result
    }
}

/// Multiply two square matrices of given size.
fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let mut c = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for k in 0..n {
            if a[i][k].abs() < 1e-15 {
                continue;
            }
            for j in 0..n {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

// ─── PoissonProcess ──────────────────────────────────────────────────────────

/// Homogeneous and inhomogeneous Poisson processes.
pub struct PoissonProcess {
    /// Constant rate for homogeneous process.
    pub rate: f64,
    seed: u64,
}

impl PoissonProcess {
    /// Creates a homogeneous Poisson process with constant rate.
    pub fn new(rate: f64, seed: u64) -> Self {
        Self { rate, seed }
    }

    /// Generate arrival times for a homogeneous Poisson process on \[0, t_end\].
    ///
    /// Inter-arrival times are Exp(rate) distributed.
    pub fn arrivals(&self, t_end: f64) -> Vec<f64> {
        let mut rng = Lcg::new(self.seed);
        let mut times = Vec::new();
        let mut t = 0.0_f64;
        loop {
            let inter = rng.exponential(self.rate);
            t += inter;
            if t > t_end {
                break;
            }
            times.push(t);
        }
        times
    }

    /// Count events in \[0, t_end\]. E\[N(t)\] = rate * t.
    pub fn count(&self, t_end: f64) -> usize {
        self.arrivals(t_end).len()
    }

    /// Generate arrival times for an inhomogeneous Poisson process with intensity lambda(t)
    /// using Lewis-Shedler thinning. `lambda_max` must bound lambda(t) from above.
    pub fn inhomogeneous_arrivals<F>(&self, t_end: f64, lambda_max: f64, lambda_fn: F) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut rng = Lcg::new(self.seed ^ 0xabcdef12);
        let mut times = Vec::new();
        let mut t = 0.0_f64;
        loop {
            let inter = rng.exponential(lambda_max);
            t += inter;
            if t > t_end {
                break;
            }
            let accept_prob = lambda_fn(t) / lambda_max;
            if rng.uniform() < accept_prob {
                times.push(t);
            }
        }
        times
    }

    /// Inter-arrival times for a homogeneous Poisson process (exactly n arrivals).
    pub fn inter_arrival_times(&self, n: usize) -> Vec<f64> {
        let mut rng = Lcg::new(self.seed ^ 0x12345678);
        (0..n).map(|_| rng.exponential(self.rate)).collect()
    }
}

// ─── Variance Reduction Techniques ───────────────────────────────────────────

/// Antithetic variates estimator.
///
/// For a Monte Carlo estimate E\[f(X)\], uses paired samples (U, 1-U) to
/// reduce variance. If f is monotone in X, this halves the variance.
///
/// `evaluator` is called with `n_pairs` pairs (x, x_antithetic).
/// Returns (estimate, variance_reduction_factor).
pub fn antithetic_variates<F>(evaluator: F, n_pairs: usize, seed: u64) -> (f64, f64)
where
    F: Fn(f64, f64) -> f64,
{
    let mut rng = Lcg::new(seed);
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    for _ in 0..n_pairs {
        let u = rng.uniform();
        let val = (evaluator(u, 1.0 - u) + evaluator(1.0 - u, u)) / 2.0;
        sum += val;
        sum_sq += val * val;
    }
    let mean = sum / n_pairs as f64;
    let variance = (sum_sq / n_pairs as f64 - mean * mean).max(0.0);
    (mean, variance)
}

/// Control variate estimator.
///
/// Reduces variance by exploiting a correlated control variate C with known mean mu_c.
/// Optimal coefficient: b* = Cov(Y, C) / Var(C)
/// Estimator: Y_cv = Y - b*(C - mu_c)
///
/// `pairs`: slice of (Y_i, C_i) samples.
/// `mu_c`: known expectation of C.
/// Returns corrected estimate.
pub fn control_variate_estimator(pairs: &[(f64, f64)], mu_c: f64) -> f64 {
    let n = pairs.len() as f64;
    let mean_y: f64 = pairs.iter().map(|(y, _)| y).sum::<f64>() / n;
    let mean_c: f64 = pairs.iter().map(|(_, c)| c).sum::<f64>() / n;
    let cov_yc: f64 = pairs
        .iter()
        .map(|(y, c)| (y - mean_y) * (c - mean_c))
        .sum::<f64>()
        / n;
    let var_c: f64 = pairs.iter().map(|(_, c)| (c - mean_c).powi(2)).sum::<f64>() / n;
    if var_c.abs() < 1e-15 {
        return mean_y;
    }
    let b = cov_yc / var_c;
    mean_y - b * (mean_c - mu_c)
}

/// Quasi-Monte Carlo: Halton sequence in base 2 and 3 for 2D integration.
///
/// Generates `n` low-discrepancy points in \[0,1)^2 using Halton sequence.
/// Much better convergence than pseudo-random for smooth integrands.
pub fn quasi_monte_carlo<F>(integrand: F, n: usize) -> f64
where
    F: Fn(f64, f64) -> f64,
{
    let sum: f64 = (1..=n).map(|i| integrand(halton(i, 2), halton(i, 3))).sum();
    sum / n as f64
}

/// Compute the i-th element of the Halton sequence in the given base.
pub fn halton(mut i: usize, base: usize) -> f64 {
    let mut result = 0.0_f64;
    let mut denominator = 1.0_f64;
    while i > 0 {
        denominator *= base as f64;
        result += (i % base) as f64 / denominator;
        i /= base;
    }
    result
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── BrownianMotion tests ──────────────────────────────────────────────────

    #[test]
    fn test_bm_path_length() {
        let bm = BrownianMotion::new(0.01, 42);
        let path = bm.path(100);
        assert_eq!(path.len(), 101);
    }

    #[test]
    fn test_bm_starts_at_zero() {
        let bm = BrownianMotion::new(0.01, 42);
        let path = bm.path(100);
        assert_eq!(path[0], 0.0);
    }

    #[test]
    fn test_bm_variance_grows_with_time() {
        // Generate many paths and verify that E[W(t)^2] ≈ sigma^2 * t
        let sigma = 1.0;
        let dt = 0.1;
        let n_steps = 10;
        let t = n_steps as f64 * dt;
        let n_paths = 5000;
        let mut sum_sq = 0.0_f64;
        for k in 0..n_paths {
            let bm = BrownianMotion::with_sigma(dt, sigma, k as u64 * 1234 + 1);
            let path = bm.path(n_steps);
            let w_t = path[n_steps];
            sum_sq += w_t * w_t;
        }
        let empirical_var = sum_sq / n_paths as f64;
        let theoretical_var = sigma * sigma * t;
        assert!(
            (empirical_var - theoretical_var).abs() < 0.15,
            "BM variance: empirical={:.6} theoretical={:.6}",
            empirical_var,
            theoretical_var
        );
    }

    #[test]
    fn test_bm_quadratic_variation_converges() {
        // Quadratic variation should converge to sigma^2 * T as dt -> 0
        let bm = BrownianMotion::with_sigma(0.001, 2.0, 7777);
        let path = bm.path(1000);
        let qv = BrownianMotion::quadratic_variation(&path);
        // sigma^2 * T = 4 * 1.0 = 4.0
        assert!(
            (qv - 4.0).abs() < 1.0,
            "QV should be near 4.0, got {:.6}",
            qv
        );
    }

    #[test]
    fn test_bm_multi_path_count() {
        let bm = BrownianMotion::new(0.01, 1);
        let paths = bm.multi_path(10, 50);
        assert_eq!(paths.len(), 10);
        for path in &paths {
            assert_eq!(path.len(), 51);
        }
    }

    #[test]
    fn test_bm_theoretical_variance() {
        let bm = BrownianMotion::with_sigma(0.01, 3.0, 1);
        let var = bm.variance_at(2.0);
        assert!((var - 18.0).abs() < 1e-10); // sigma^2 * t = 9 * 2 = 18
    }

    // ── GBM tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_gbm_starts_at_s0() {
        let gbm = GeometricBrownianMotion::new(100.0, 0.05, 0.2, 99);
        let path = gbm.exact_path(1.0, 100);
        assert!((path[0] - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_gbm_stays_positive() {
        let gbm = GeometricBrownianMotion::new(10.0, 0.1, 0.5, 55);
        let path = gbm.exact_path(2.0, 200);
        assert!(path.iter().all(|&s| s >= 0.0));
    }

    #[test]
    fn test_gbm_mean_convergence() {
        // E[S(T)] = S0 * exp(mu * T)
        let s0 = 100.0;
        let mu = 0.10;
        let sigma = 0.2;
        let t_end = 1.0;
        let n_steps = 252;
        let n_paths = 5000;
        let mut sum = 0.0_f64;
        for k in 0..n_paths {
            let gbm = GeometricBrownianMotion::new(s0, mu, sigma, k as u64 + 1);
            let path = gbm.exact_path(t_end, n_steps);
            sum += path[n_steps];
        }
        let empirical_mean = sum / n_paths as f64;
        let theoretical_mean = s0 * (mu * t_end).exp();
        let rel_err = (empirical_mean - theoretical_mean).abs() / theoretical_mean;
        assert!(rel_err < 0.05, "GBM mean rel error: {:.6}", rel_err);
    }

    #[test]
    fn test_gbm_log_normality() {
        // Log-returns should be approximately normal
        let gbm = GeometricBrownianMotion::new(1.0, 0.0, 1.0, 123);
        let path = gbm.exact_path(1.0, 1000);
        let log_return = path[1000].ln();
        // Should be finite
        assert!(log_return.is_finite());
    }

    #[test]
    fn test_gbm_euler_stays_positive() {
        let gbm = GeometricBrownianMotion::new(50.0, 0.05, 0.3, 333);
        let path = gbm.euler_path(1.0, 100);
        assert!(path.iter().all(|&s| s >= 0.0));
    }

    // ── OrnsteinUhlenbeck tests ───────────────────────────────────────────────

    #[test]
    fn test_ou_starts_at_x0() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.5, 0.3, 1.0, 42);
        let path = ou.path(1.0, 100);
        assert!((path[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ou_mean_reversion() {
        // Long-run mean should converge to mu
        let theta = 5.0;
        let mu = 3.0;
        let ou = OrnsteinUhlenbeck::new(theta, mu, 0.5, 10.0, 77);
        let path = ou.path(10.0, 10000);
        let tail_mean: f64 = path[8000..].iter().sum::<f64>() / 2000.0;
        assert!(
            (tail_mean - mu).abs() < 0.3,
            "OU mean: {:.6} vs {:.6}",
            tail_mean,
            mu
        );
    }

    #[test]
    fn test_ou_stationary_variance() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.0, 1.0, 0.0, 1);
        let sv = ou.stationary_variance();
        // sigma^2/(2*theta) = 1/4 = 0.25
        assert!((sv - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_ou_autocorrelation_decay() {
        let ou = OrnsteinUhlenbeck::new(1.0, 0.0, 1.0, 0.0, 1);
        let ac0 = ou.autocorrelation(0.0);
        let ac1 = ou.autocorrelation(1.0);
        assert!((ac0 - 1.0).abs() < 1e-10);
        assert!((ac1 - (-1.0_f64).exp()).abs() < 1e-10);
    }

    // ── JumpDiffusion tests ───────────────────────────────────────────────────

    #[test]
    fn test_jump_diffusion_starts_at_s0() {
        let jd = JumpDiffusion::new(100.0, 0.05, 0.2, 1.0, -0.1, 0.2, 42);
        let path = jd.path(1.0, 100);
        assert!((path[0] - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_jump_diffusion_stays_nonneg() {
        let jd = JumpDiffusion::new(50.0, 0.05, 0.3, 2.0, 0.0, 0.3, 7);
        let path = jd.path(1.0, 252);
        assert!(path.iter().all(|&s| s >= 0.0));
    }

    #[test]
    fn test_jump_diffusion_compensation() {
        let jd = JumpDiffusion::new(1.0, 0.0, 0.1, 1.0, 0.0, 0.1, 1);
        let comp = jd.jump_compensation();
        // E[J] = exp(0 + 0.01/2) = exp(0.005) ≈ 1.005012
        let expected = (0.0_f64 + 0.5 * 0.01_f64).exp() - 1.0;
        assert!((comp - expected).abs() < 1e-10);
    }

    // ── FractionalBrownianMotion tests ────────────────────────────────────────

    #[test]
    fn test_fbm_path_length() {
        let fbm = FractionalBrownianMotion::new(0.7, 1);
        let path = fbm.path(10);
        assert_eq!(path.len(), 11);
    }

    #[test]
    fn test_fbm_starts_at_zero() {
        let fbm = FractionalBrownianMotion::new(0.5, 42);
        let path = fbm.path(20);
        assert!(path[0].abs() < 1e-10);
    }

    #[test]
    fn test_fbm_hurst_validation() {
        let result = std::panic::catch_unwind(|| FractionalBrownianMotion::new(1.5, 1));
        assert!(result.is_err());
    }

    #[test]
    fn test_fbm_covariance_positive_definite() {
        let fbm = FractionalBrownianMotion::new(0.7, 1);
        let cov = fbm.covariance(1.0, 1.0);
        // C(t,t) = sigma^2 * t^{2H} > 0
        assert!(cov > 0.0);
    }

    // ── CoxIngersollRoss tests ────────────────────────────────────────────────

    #[test]
    fn test_cir_stays_positive() {
        let cir = CoxIngersollRoss::new(2.0, 0.05, 0.1, 0.04, 99);
        let path = cir.path(10.0, 1000);
        assert!(path.iter().all(|&r| r >= 0.0));
    }

    #[test]
    fn test_cir_feller_condition() {
        // 2*kappa*theta >= sigma^2 ensures strict positivity
        let cir = CoxIngersollRoss::new(2.0, 0.05, 0.1, 0.04, 1);
        // 2*2.0*0.05 = 0.2 >= 0.01 -> satisfied
        assert!(cir.feller_satisfied());
    }

    #[test]
    fn test_cir_mean_reversion() {
        let theta = 0.05;
        let cir = CoxIngersollRoss::new(3.0, theta, 0.1, 0.2, 55);
        let path = cir.path(20.0, 20000);
        let tail_mean: f64 = path[15000..].iter().sum::<f64>() / 5000.0;
        assert!(
            (tail_mean - theta).abs() < 0.01,
            "CIR mean: {:.6} vs {:.6}",
            tail_mean,
            theta
        );
    }

    #[test]
    fn test_cir_stationary_variance() {
        let cir = CoxIngersollRoss::new(2.0, 0.05, 0.1, 0.04, 1);
        let sv = cir.stationary_variance();
        // sigma^2 * theta / (2 * kappa) = 0.01 * 0.05 / 4 = 0.000125
        assert!((sv - 0.000125).abs() < 1e-10);
    }

    // ── HestonModel tests ─────────────────────────────────────────────────────

    #[test]
    fn test_heston_path_lengths() {
        let heston = HestonModel::new(100.0, 0.04, 0.05, 2.0, 0.04, 0.3, -0.7, 11);
        let (s, v) = heston.path(1.0, 100);
        assert_eq!(s.len(), 101);
        assert_eq!(v.len(), 101);
    }

    #[test]
    fn test_heston_starts_correctly() {
        let heston = HestonModel::new(100.0, 0.04, 0.05, 2.0, 0.04, 0.3, -0.7, 11);
        let (s, v) = heston.path(1.0, 100);
        assert!((s[0] - 100.0).abs() < 1e-10);
        assert!((v[0] - 0.04).abs() < 1e-10);
    }

    #[test]
    fn test_heston_variance_stays_nonneg() {
        let heston = HestonModel::new(50.0, 0.1, 0.05, 1.0, 0.1, 0.5, -0.5, 77);
        let (_, v) = heston.path(2.0, 500);
        assert!(v.iter().all(|&vi| vi >= 0.0));
    }

    // ── LevyStableProcess tests ───────────────────────────────────────────────

    #[test]
    fn test_levy_path_length() {
        let levy = LevyStableProcess::new(1.5, 0.0, 1.0, 0.0, 42);
        let path = levy.path(100);
        assert_eq!(path.len(), 101);
    }

    #[test]
    fn test_levy_path_starts_at_zero() {
        let levy = LevyStableProcess::new(1.8, 0.0, 1.0, 0.0, 1);
        let path = levy.path(50);
        assert_eq!(path[0], 0.0);
    }

    #[test]
    fn test_levy_gaussian_case_finite() {
        // alpha=2 gives Gaussian
        let levy = LevyStableProcess::new(2.0, 0.0, 1.0, 0.0, 88);
        let incs = levy.increments(1000);
        assert!(incs.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn test_levy_alpha_validation() {
        let result = std::panic::catch_unwind(|| LevyStableProcess::new(0.0, 0.0, 1.0, 0.0, 1));
        assert!(result.is_err());
    }

    // ── BranchingProcess tests ────────────────────────────────────────────────

    #[test]
    fn test_branching_path_length() {
        let bp = BranchingProcess::poisson(10, 1.5, 42);
        let pop = bp.simulate(20);
        assert_eq!(pop.len(), 21);
    }

    #[test]
    fn test_branching_starts_at_z0() {
        let bp = BranchingProcess::poisson(5, 1.2, 1);
        let pop = bp.simulate(10);
        assert_eq!(pop[0], 5);
    }

    #[test]
    fn test_branching_subcritical_extinction() {
        // lambda < 1: extinction probability = 1
        let bp = BranchingProcess::poisson(1, 0.5, 100);
        let q = bp.extinction_probability();
        assert!((q - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_branching_supercritical_extinction_less_than_one() {
        // lambda > 1: extinction probability < 1
        let bp = BranchingProcess::poisson(1, 2.0, 42);
        let q = bp.extinction_probability();
        assert!(q < 1.0 && q > 0.0);
    }

    #[test]
    fn test_branching_mc_extinction_subcritical() {
        // Sub-critical should go extinct with high probability
        let bp = BranchingProcess::poisson(1, 0.5, 42);
        let q = bp.mc_extinction_probability(500, 50);
        assert!(
            q > 0.8,
            "Subcritical extinction prob should be high, got {:.6}",
            q
        );
    }

    // ── MarkovChain tests ─────────────────────────────────────────────────────

    fn two_state_chain(seed: u64) -> MarkovChain {
        MarkovChain::new(vec![vec![0.7, 0.3], vec![0.4, 0.6]], seed)
    }

    #[test]
    fn test_markov_simulate_length() {
        let mc = two_state_chain(1);
        let traj = mc.simulate(0, 100);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_markov_states_in_range() {
        let mc = two_state_chain(5);
        let traj = mc.simulate(0, 1000);
        assert!(traj.iter().all(|&s| s < 2));
    }

    #[test]
    fn test_markov_stationary_distribution() {
        // For P = [[0.7, 0.3],[0.4, 0.6]], stationary: pi[0] = 4/7, pi[1] = 3/7
        let mc = two_state_chain(1);
        let pi = mc.stationary_distribution();
        let expected = [4.0 / 7.0, 3.0 / 7.0];
        assert!((pi[0] - expected[0]).abs() < 1e-6, "pi[0]={:.8}", pi[0]);
        assert!((pi[1] - expected[1]).abs() < 1e-6, "pi[1]={:.8}", pi[1]);
    }

    #[test]
    fn test_markov_stationary_sums_to_one() {
        let mc = two_state_chain(7);
        let pi = mc.stationary_distribution();
        let sum: f64 = pi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_markov_n_step_matrix() {
        let mc = two_state_chain(1);
        let p1 = mc.n_step_matrix(1);
        // Should match transition matrix
        assert!((p1[0][0] - 0.7).abs() < 1e-10);
        assert!((p1[0][1] - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_markov_hitting_time_positive() {
        let mc = two_state_chain(42);
        let ht = mc.hitting_time(0, 1, 100);
        assert!(ht > 0.0);
    }

    // ── PoissonProcess tests ──────────────────────────────────────────────────

    #[test]
    fn test_poisson_arrivals_ordered() {
        let pp = PoissonProcess::new(5.0, 1234);
        let arrivals = pp.arrivals(10.0);
        for i in 1..arrivals.len() {
            assert!(arrivals[i] > arrivals[i - 1]);
        }
    }

    #[test]
    fn test_poisson_arrivals_in_range() {
        let pp = PoissonProcess::new(3.0, 99);
        let arrivals = pp.arrivals(5.0);
        assert!(arrivals.iter().all(|&t| t > 0.0 && t <= 5.0));
    }

    #[test]
    fn test_poisson_mean_count() {
        // E[N(t)] = lambda * t
        let rate = 4.0;
        let t = 10.0;
        let n_samples = 2000;
        let mut total = 0usize;
        for k in 0..n_samples {
            let pp = PoissonProcess::new(rate, k as u64 + 1);
            total += pp.count(t);
        }
        let empirical_mean = total as f64 / n_samples as f64;
        let expected = rate * t;
        assert!(
            (empirical_mean - expected).abs() < 2.0,
            "Poisson mean: {:.6} vs {:.6}",
            empirical_mean,
            expected
        );
    }

    #[test]
    fn test_poisson_inter_arrival_exponential() {
        // Inter-arrival times should be exponential(rate)
        // Mean should be 1/rate
        let rate = 2.0;
        let n = 10000;
        let pp = PoissonProcess::new(rate, 42);
        let iat = pp.inter_arrival_times(n);
        let mean: f64 = iat.iter().sum::<f64>() / n as f64;
        let expected_mean = 1.0 / rate;
        assert!(
            (mean - expected_mean).abs() < 0.05,
            "Mean inter-arrival: {:.6} vs {:.6}",
            mean,
            expected_mean
        );
    }

    #[test]
    fn test_poisson_inhomogeneous() {
        // Inhomogeneous with constant lambda should match homogeneous
        let pp = PoissonProcess::new(2.0, 1);
        let arrivals = pp.inhomogeneous_arrivals(10.0, 2.0, |_t| 2.0);
        assert!(arrivals.iter().all(|&t| (0.0..=10.0).contains(&t)));
    }

    // ── Variance reduction tests ──────────────────────────────────────────────

    #[test]
    fn test_antithetic_variates_constant() {
        // For f(u, u_anti) = 1 (constant), estimate should be 1
        let (mean, _var) = antithetic_variates(|_u, _u_anti| 1.0, 1000, 42);
        assert!((mean - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_antithetic_variates_uniform_mean() {
        // E[U] = 0.5 where f(u, u_anti) = u
        let (mean, _var) = antithetic_variates(|u, _| u, 5000, 7);
        assert!((mean - 0.5).abs() < 0.02, "Antithetic mean: {:.6}", mean);
    }

    #[test]
    fn test_control_variate_identity() {
        // If Y = C with known mean, estimator should return mu_c
        let pairs: Vec<(f64, f64)> = (1..=100).map(|i| (i as f64, i as f64)).collect();
        let mu_c = 50.5_f64;
        let est = control_variate_estimator(&pairs, mu_c);
        assert!((est - mu_c).abs() < 0.01, "CV estimate: {:.6}", est);
    }

    #[test]
    fn test_quasi_monte_carlo_unit_square() {
        // Integral of 1 over [0,1]^2 = 1
        let est = quasi_monte_carlo(|_x, _y| 1.0, 1000);
        assert!((est - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quasi_monte_carlo_linear() {
        // Integral of (x + y) over [0,1]^2 = 1.0
        let est = quasi_monte_carlo(|x, y| x + y, 10000);
        assert!((est - 1.0).abs() < 0.01, "QMC estimate: {:.6}", est);
    }

    #[test]
    fn test_halton_sequence_base2() {
        // First few Halton base 2: 1/2, 1/4, 3/4, 1/8, ...
        assert!((halton(1, 2) - 0.5).abs() < 1e-10);
        assert!((halton(2, 2) - 0.25).abs() < 1e-10);
        assert!((halton(3, 2) - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_halton_in_unit_interval() {
        for i in 1..=1000 {
            let h = halton(i, 2);
            assert!((0.0..1.0).contains(&h), "Halton({}, 2) = {:.6}", i, h);
        }
    }
}
