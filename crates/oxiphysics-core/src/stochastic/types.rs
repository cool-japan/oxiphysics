//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::*;
use rand::RngExt;
use std::f64::consts::PI;

/// Estimate the mean first passage time empirically from an ensemble of paths.
pub struct EmpiricalFirstPassageTime;
impl EmpiricalFirstPassageTime {
    /// Estimate the mean first passage time for an OU process to reach `barrier`.
    ///
    /// Simulates `n_paths` realisations and records the first time each path
    /// crosses `barrier`. Returns the mean of the recorded times (ignoring
    /// paths that did not cross within `t_max`).
    ///
    /// # Arguments
    /// * `ou`      - Ornstein-Uhlenbeck parameters.
    /// * `x0`      - Starting position.
    /// * `barrier` - Barrier level to cross.
    /// * `dt`      - Time step for simulation.
    /// * `_t_min`  - (Unused) minimum time threshold for recording.
    /// * `t_max`   - Maximum simulation time.
    /// * `n_paths` - Number of Monte Carlo paths.
    /// * `seed`    - Base random seed.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_ou(
        ou: &OrnsteinUhlenbeck,
        x0: f64,
        barrier: f64,
        dt: f64,
        _t_min: f64,
        t_max: f64,
        n_paths: usize,
        seed: u64,
    ) -> f64 {
        let n_steps = (t_max / dt).ceil() as usize;
        let mut fpt_sum = 0.0_f64;
        let mut n_crossed = 0usize;
        for p in 0..n_paths {
            let mut rng =
                Rng::new(seed.wrapping_add((p as u64).wrapping_mul(6_364_136_223_846_793_005)));
            let mut x = x0;
            for step in 0..n_steps {
                let z = rng.next_normal();
                let exp_factor = (-ou.theta * dt).exp();
                let mean = ou.mu + (x - ou.mu) * exp_factor;
                let std = (ou.sigma * ou.sigma / (2.0 * ou.theta)
                    * (1.0 - exp_factor * exp_factor))
                    .sqrt();
                x = mean + std * z;
                if x >= barrier {
                    fpt_sum += (step + 1) as f64 * dt;
                    n_crossed += 1;
                    break;
                }
            }
        }
        if n_crossed == 0 {
            t_max
        } else {
            fpt_sum / n_crossed as f64
        }
    }
}
/// Hull-White one-factor interest rate model.
///
/// dr = (theta(t) - a*r)*dt + sigma*dW
///
/// where theta(t) is calibrated to fit the initial term structure.
/// Here we use a constant theta for simplicity.
#[allow(dead_code)]
pub struct HullWhiteModel {
    /// Mean-reversion speed a.
    pub a: f64,
    /// Volatility σ.
    pub sigma: f64,
    /// Long-run level theta (constant approximation).
    pub theta: f64,
}
impl HullWhiteModel {
    /// Create a new Hull-White model.
    pub fn new(a: f64, sigma: f64, theta: f64) -> Self {
        Self { a, sigma, theta }
    }
    /// Simulate a short rate path using Euler-Maruyama.
    ///
    /// Returns a vector of short rate values.
    pub fn simulate_short_rate(&self, r0: f64, t_end: f64, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut path = Vec::with_capacity(n);
        let mut r = r0;
        path.push(r);
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            let dw = rng.next_normal() * dt.sqrt();
            r += (self.theta - self.a * r) * dt + self.sigma * dw;
            t += dt;
            path.push(r);
        }
        path
    }
    /// Analytical mean of r(t) given r(0):
    /// E\[r(t)\] = r(0)*exp(-a*t) + (theta/a)*(1 - exp(-a*t))
    pub fn mean_rate(&self, r0: f64, t: f64) -> f64 {
        let e = (-self.a * t).exp();
        r0 * e + (self.theta / self.a) * (1.0 - e)
    }
    /// Analytical variance of r(t):
    /// Var\[r(t)\] = sigma^2 / (2*a) * (1 - exp(-2*a*t))
    pub fn variance_rate(&self, t: f64) -> f64 {
        self.sigma * self.sigma / (2.0 * self.a) * (1.0 - (-2.0 * self.a * t).exp())
    }
    /// Bond price P(0, T) under the Hull-White model (analytic formula).
    ///
    /// Returns the zero-coupon bond price for maturity T.
    pub fn bond_price(&self, r0: f64, t_maturity: f64) -> f64 {
        let a = self.a;
        let sigma = self.sigma;
        let theta = self.theta;
        let b = (1.0 - (-a * t_maturity).exp()) / a;
        let log_a = (theta / a - sigma * sigma / (2.0 * a * a)) * (b - t_maturity)
            - sigma * sigma * b * b / (4.0 * a);
        (log_a - b * r0).exp()
    }
}
/// Geometric Brownian motion: dS = mu*S*dt + sigma*S*dW.
///
/// Commonly used for asset price models and turbulence statistics.
pub struct GeometricBrownianMotion {
    /// Drift coefficient.
    pub mu: f64,
    /// Volatility (diffusion) coefficient.
    pub sigma: f64,
}
impl GeometricBrownianMotion {
    /// Creates a new GBM with the given drift `mu` and volatility `sigma`.
    pub fn new(mu: f64, sigma: f64) -> Self {
        Self { mu, sigma }
    }
    /// Simulates a path using the Euler-Maruyama scheme.
    ///
    /// Returns a vector of S values at each time step starting from `s0`.
    pub fn simulate(&self, s0: f64, t_end: f64, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut path = Vec::with_capacity(n);
        let mut s = s0;
        path.push(s);
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            let dw = rng.next_normal() * dt.sqrt();
            s += self.mu * s * dt + self.sigma * s * dw;
            t += dt;
            path.push(s);
        }
        path
    }
    /// Simulates a path using the exact (log-normal) solution.
    ///
    /// S(t+dt) = S(t) * exp((mu - sigma^2/2)*dt + sigma*sqrt(dt)*Z)
    pub fn simulate_exact(&self, s0: f64, t_end: f64, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut path = Vec::with_capacity(n);
        let mut s = s0;
        path.push(s);
        let drift = (self.mu - 0.5 * self.sigma * self.sigma) * dt;
        let vol = self.sigma * dt.sqrt();
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            let z = rng.next_normal();
            s *= (drift + vol * z).exp();
            t += dt;
            path.push(s);
        }
        path
    }
    /// Returns the analytical mean E\[S(t)\] = s0 * exp(mu * t).
    pub fn analytical_mean(&self, s0: f64, t: f64) -> f64 {
        s0 * (self.mu * t).exp()
    }
    /// Returns the analytical variance Var\[S(t)\] = s0^2 * exp(2*mu*t) * (exp(sigma^2*t) - 1).
    pub fn analytical_variance(&self, s0: f64, t: f64) -> f64 {
        s0 * s0 * (2.0 * self.mu * t).exp() * ((self.sigma * self.sigma * t).exp() - 1.0)
    }
}
/// Langevin dynamics integrator for a Brownian thermostat in molecular dynamics.
///
/// Implements the stochastic equation of motion:
/// `m*dv/dt = F - gamma*m*v + sqrt(2*gamma*m*k_B*T) * xi(t)`
///
/// where `xi(t)` is Gaussian white noise.
pub struct LangevinDynamics {
    /// Particle mass (kg).
    pub mass: f64,
    /// Friction coefficient gamma (1/s).
    pub friction: f64,
    /// Temperature T (K); k_B is incorporated in noise amplitude.
    pub temperature: f64,
}
impl LangevinDynamics {
    /// Creates a new Langevin integrator.
    pub fn new(mass: f64, friction: f64, temperature: f64) -> Self {
        Self {
            mass,
            friction,
            temperature,
        }
    }
    /// Returns the noise amplitude: `sqrt(2 * friction * k_B * T / mass * dt)`.
    pub fn noise_amplitude(&self, dt: f64) -> f64 {
        (2.0 * self.friction * K_B * self.temperature / self.mass * dt).sqrt()
    }
    /// Performs a single Euler-Langevin step.
    ///
    /// Returns the updated `(position, velocity)` after time `dt`.
    ///
    /// Update rule:
    /// `v_new = v + (F/m - gamma*v)*dt + noise_amp * N(0,1) / sqrt(dt)`
    /// `x_new = x + v_new * dt`
    pub fn step(&self, pos: f64, vel: f64, force: f64, dt: f64, rng: &mut Rng) -> (f64, f64) {
        let noise = self.noise_amplitude(dt) * rng.next_normal() / dt.sqrt();
        let v_new = vel + (force / self.mass - self.friction * vel) * dt + noise;
        let x_new = pos + v_new * dt;
        (x_new, v_new)
    }
    /// Performs a BAOAB Langevin integrator step (symmetric splitting).
    ///
    /// The BAOAB scheme has better sampling properties for equilibrium
    /// thermodynamics than simple Euler-Langevin.
    ///
    /// Returns updated `(position, velocity)`.
    pub fn step_baoab(&self, pos: f64, vel: f64, force: f64, dt: f64, rng: &mut Rng) -> (f64, f64) {
        let half_dt = 0.5 * dt;
        let v_half = vel + (force / self.mass) * half_dt;
        let x_half = pos + v_half * half_dt;
        let c1 = (-self.friction * dt).exp();
        let c2 = ((1.0 - c1 * c1) * K_B * self.temperature / self.mass).sqrt();
        let v_o = c1 * v_half + c2 * rng.next_normal();
        let x_new = x_half + v_o * half_dt;
        let v_new = v_o + (force / self.mass) * half_dt;
        (x_new, v_new)
    }
}
/// Underdamped Langevin (Kramers) dynamics.
///
/// Solves the pair:
/// ```text
/// dx/dt = v
/// m*dv/dt = F - gamma*v + sqrt(2*gamma*kT/m) * ξ(t)
/// ```
#[allow(non_snake_case)]
pub struct KleinmanKramers {
    /// Friction coefficient γ.
    pub gamma: f64,
    /// Particle mass m.
    pub m: f64,
    /// Thermal energy k_B * T.
    pub kT: f64,
}
impl KleinmanKramers {
    /// Create a new `KleinmanKramers` integrator.
    #[allow(non_snake_case)]
    pub fn new(gamma: f64, m: f64, kT: f64) -> Self {
        Self { gamma, m, kT }
    }
    /// Euler-Maruyama step for underdamped Langevin.
    ///
    /// Returns `(new_x, new_v)`.
    #[allow(non_snake_case)]
    pub fn step(&self, x: f64, v: f64, F: f64, dt: f64, rng: &mut impl rand::Rng) -> (f64, f64) {
        let noise_amp = (2.0 * self.gamma * self.kT / self.m * dt).sqrt();
        let u1: f64 = loop {
            let val: f64 = rng.random();
            if val > 0.0 {
                break val;
            }
        };
        let u2: f64 = rng.random();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        let dv = (F / self.m - self.gamma / self.m * v) * dt + noise_amp * z;
        let v_new = v + dv;
        let x_new = x + v * dt;
        (x_new, v_new)
    }
}
/// Variance Gamma (VG) process: a Brownian motion with drift subordinated by a
/// Gamma time change.
///
/// X(t) = theta * G(t) + sigma * W(G(t))
///
/// where G(t) is a Gamma process with mean `t` and variance `nu * t`.
#[allow(dead_code)]
pub struct VarianceGammaProcess {
    /// Drift of the Brownian motion θ.
    pub theta: f64,
    /// Volatility of the Brownian motion σ.
    pub sigma: f64,
    /// Variance rate of the Gamma subordinator ν.
    pub nu: f64,
}
impl VarianceGammaProcess {
    /// Create a new Variance Gamma process.
    pub fn new(theta: f64, sigma: f64, nu: f64) -> Self {
        Self { theta, sigma, nu }
    }
    /// Sample a Gamma(shape, scale) variate using Marsaglia-Tsang method.
    fn sample_gamma(shape: f64, scale: f64, rng: &mut Rng) -> f64 {
        if shape < 1.0 {
            let u = loop {
                let v = rng.next_f64();
                if v > 0.0 {
                    break v;
                }
            };
            return Self::sample_gamma(shape + 1.0, scale, rng) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = rng.next_normal();
            let v = 1.0 + c * x;
            if v <= 0.0 {
                continue;
            }
            let v3 = v * v * v;
            let u = loop {
                let u = rng.next_f64();
                if u > 0.0 {
                    break u;
                }
            };
            let x2 = x * x;
            if u < 1.0 - 0.0331 * x2 * x2 {
                return scale * d * v3;
            }
            if u.ln() < 0.5 * x2 + d * (1.0 - v3 + v3.ln()) {
                return scale * d * v3;
            }
        }
    }
    /// Simulate a VG path over `[0, t_end]` with `n_steps` steps.
    ///
    /// Returns a vector of cumulative VG process values.
    pub fn simulate(&self, t_end: f64, n_steps: usize, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let dt = t_end / n_steps as f64;
        let gamma_shape = dt / self.nu;
        let gamma_scale = self.nu;
        let mut path = Vec::with_capacity(n_steps + 1);
        let mut x = 0.0_f64;
        path.push(x);
        for _ in 0..n_steps {
            let dg = Self::sample_gamma(gamma_shape, gamma_scale, &mut rng);
            let dw = rng.next_normal() * dg.sqrt();
            x += self.theta * dg + self.sigma * dw;
            path.push(x);
        }
        path
    }
    /// Returns the analytical mean increment over time `dt`: E\[ΔX\] = θ * dt.
    pub fn mean_increment(&self, dt: f64) -> f64 {
        self.theta * dt
    }
    /// Returns the analytical variance increment over time `dt`:
    /// Var\[ΔX\] = (σ² + θ²*ν) * dt.
    pub fn variance_increment(&self, dt: f64) -> f64 {
        (self.sigma * self.sigma + self.theta * self.theta * self.nu) * dt
    }
}
/// SABR stochastic volatility model.
///
/// dF = alpha * F^beta * dW_1
/// dalpha = nu * alpha * dW_2
/// Corr(dW_1, dW_2) = rho
///
/// Used for interest rate derivatives and FX options.
#[allow(dead_code)]
pub struct SabrModel {
    /// Initial forward rate F_0.
    pub f0: f64,
    /// Initial volatility α_0.
    pub alpha0: f64,
    /// CEV exponent β ∈ \[0, 1\].
    pub beta: f64,
    /// Volatility of volatility ν.
    pub nu: f64,
    /// Correlation ρ between F and α processes.
    pub rho: f64,
}
impl SabrModel {
    /// Create a new SABR model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(f0: f64, alpha0: f64, beta: f64, nu: f64, rho: f64) -> Self {
        Self {
            f0,
            alpha0,
            beta,
            nu,
            rho,
        }
    }
    /// Simulate forward rate and volatility paths using Euler-Maruyama.
    ///
    /// Returns `(forward_rates, volatilities)` vectors.
    pub fn simulate(&self, t_end: f64, n_steps: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
        let mut rng = Rng::new(seed);
        let dt = t_end / n_steps as f64;
        let sqrt_dt = dt.sqrt();
        let mut f = self.f0;
        let mut alpha = self.alpha0;
        let mut forwards = Vec::with_capacity(n_steps + 1);
        let mut vols = Vec::with_capacity(n_steps + 1);
        forwards.push(f);
        vols.push(alpha);
        for _ in 0..n_steps {
            let (z1, z2_raw) = rng.next_normal_pair();
            let dw1 = z1 * sqrt_dt;
            let dw2 = (self.rho * z1 + (1.0 - self.rho * self.rho).sqrt() * z2_raw) * sqrt_dt;
            let f_power = f.abs().powf(self.beta);
            f += alpha * f_power * dw1;
            f = f.max(0.0);
            alpha += self.nu * alpha * dw2;
            alpha = alpha.max(1e-12);
            forwards.push(f);
            vols.push(alpha);
        }
        (forwards, vols)
    }
    /// Approximate implied volatility using the Hagan et al. (2002) formula.
    ///
    /// Valid for F close to K and small α.
    pub fn implied_vol_approx(&self, strike: f64, t: f64) -> f64 {
        let f = self.f0;
        let k = strike;
        let alpha = self.alpha0;
        let beta = self.beta;
        let nu = self.nu;
        let rho = self.rho;
        if (f - k).abs() < 1e-8 {
            let fmid = f.powf(1.0 - beta);
            let term1 = (1.0 - beta).powi(2) * alpha * alpha / (24.0 * fmid * fmid);
            let term2 = rho * beta * nu * alpha / (4.0 * f.powf(1.0 - beta));
            let term3 = (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
            alpha / fmid * (1.0 + (term1 + term2 + term3) * t)
        } else {
            let log_fk = (f / k).ln();
            let fk_mid = (f * k).powf((1.0 - beta) / 2.0);
            let z = (nu / alpha) * fk_mid * log_fk;
            let x_z = if z.abs() < 1e-6 {
                1.0
            } else {
                let inner = (1.0 - 2.0 * rho * z + z * z).sqrt() + z - rho;
                (inner / (1.0 - rho)).ln() / z
            };
            let denom = fk_mid
                * (1.0
                    + (1.0 - beta).powi(2) / 24.0 * log_fk.powi(2)
                    + (1.0 - beta).powi(4) / 1920.0 * log_fk.powi(4));
            let term1 = (1.0 - beta).powi(2) * alpha * alpha / (24.0 * (f * k).powf(1.0 - beta));
            let term2 = rho * beta * nu * alpha / (4.0 * fk_mid);
            let term3 = (2.0 - 3.0 * rho * rho) * nu * nu / 24.0;
            let iv = alpha / denom * (z / x_z) * (1.0 + (term1 + term2 + term3) * t);
            iv.max(1e-10)
        }
    }
}
/// Monte Carlo estimator for GBM payoffs using a control variate.
///
/// Uses the terminal asset price `S_T` as the control variate, since
/// `E[S_T] = S0 * exp(mu * T)` is known analytically.
pub struct ControlVariateGbm;
impl ControlVariateGbm {
    /// Estimate `E[payoff(S_T)]` using the control variate `S_T`.
    ///
    /// # Arguments
    /// * `s0`      - Initial asset price.
    /// * `gbm`     - GBM parameters.
    /// * `t`       - Time horizon.
    /// * `dt`      - Simulation time step.
    /// * `n_paths` - Number of Monte Carlo paths.
    /// * `payoff`  - Payoff function of terminal price.
    /// * `seed`    - Random seed.
    ///
    /// Returns `(estimate, standard_error)`.
    pub fn estimate<F>(
        s0: f64,
        gbm: &GeometricBrownianMotion,
        t: f64,
        dt: f64,
        n_paths: usize,
        payoff: F,
        seed: u64,
    ) -> (f64, f64)
    where
        F: Fn(f64) -> f64,
    {
        let n_steps = (t / dt).ceil() as usize;
        let sqrt_dt = dt.sqrt();
        let mu_analytic = gbm.analytical_mean(s0, t);
        let mut payoffs = Vec::with_capacity(n_paths);
        let mut terminals = Vec::with_capacity(n_paths);
        for p in 0..n_paths {
            let mut rng = Rng::new(seed.wrapping_add(p as u64 * 1_013_904_223));
            let mut s = s0;
            for _ in 0..n_steps {
                let z = rng.next_normal();
                s *= ((gbm.mu - 0.5 * gbm.sigma * gbm.sigma) * dt + gbm.sigma * sqrt_dt * z).exp();
            }
            payoffs.push(payoff(s));
            terminals.push(s);
        }
        let n = n_paths as f64;
        let mean_payoff: f64 = payoffs.iter().sum::<f64>() / n;
        let mean_terminal: f64 = terminals.iter().sum::<f64>() / n;
        let cov: f64 = payoffs
            .iter()
            .zip(terminals.iter())
            .map(|(p, s)| (p - mean_payoff) * (s - mean_terminal))
            .sum::<f64>()
            / (n - 1.0).max(1.0);
        let var_terminal: f64 = terminals
            .iter()
            .map(|s| (s - mean_terminal).powi(2))
            .sum::<f64>()
            / (n - 1.0).max(1.0);
        let beta = if var_terminal > 1e-15 {
            cov / var_terminal
        } else {
            0.0
        };
        let cv_estimates: Vec<f64> = payoffs
            .iter()
            .zip(terminals.iter())
            .map(|(p, s)| p - beta * (s - mu_analytic))
            .collect();
        let cv_mean: f64 = cv_estimates.iter().sum::<f64>() / n;
        let cv_var: f64 = cv_estimates
            .iter()
            .map(|v| (v - cv_mean).powi(2))
            .sum::<f64>()
            / (n - 1.0).max(1.0);
        (cv_mean, (cv_var / n).sqrt())
    }
}
/// Exact discrete-time sampler for the Ornstein-Uhlenbeck process.
///
/// Uses the exact conditional distribution:
/// `X(t+dt) | X(t) ~ N(mu + (x - mu)*exp(-theta*dt), sigma²/(2*theta) * (1 - exp(-2*theta*dt)))`
#[derive(Clone)]
pub struct OuExactSampler {
    /// OU parameters.
    pub(super) ou: OrnsteinUhlenbeck,
}
impl OuExactSampler {
    /// Create a new exact OU sampler.
    pub fn new(ou: OrnsteinUhlenbeck) -> Self {
        Self { ou }
    }
    /// Draw the next sample given current value `x` and time step `dt`.
    pub fn step(&self, x: f64, dt: f64, rng: &mut Rng) -> f64 {
        let exp_factor = (-self.ou.theta * dt).exp();
        let mean = self.ou.mu + (x - self.ou.mu) * exp_factor;
        let var =
            self.ou.sigma * self.ou.sigma / (2.0 * self.ou.theta) * (1.0 - exp_factor * exp_factor);
        let std = var.sqrt();
        mean + std * rng.next_normal()
    }
    /// Generate a trajectory of `n_steps` using the exact sampler.
    pub fn simulate(&self, x0: f64, n_steps: usize, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let mut path = Vec::with_capacity(n_steps + 1);
        path.push(x0);
        let mut x = x0;
        for _ in 0..n_steps {
            x = self.step(x, dt, &mut rng);
            path.push(x);
        }
        path
    }
}
/// Merton jump-diffusion model: dS/S = (mu - lambda*m_j)*dt + sigma*dW + J*dN.
///
/// Extends GBM with compound Poisson jumps, where each jump size is
/// log-normally distributed: ln(1+J) ~ N(mu_j, sigma_j^2).
pub struct MertonJumpDiffusion {
    /// Drift coefficient (risk-neutral rate).
    pub mu: f64,
    /// Diffusion volatility.
    pub sigma: f64,
    /// Jump intensity (average number of jumps per unit time).
    pub lambda: f64,
    /// Mean of log-jump size.
    pub mu_j: f64,
    /// Std of log-jump size.
    pub sigma_j: f64,
}
impl MertonJumpDiffusion {
    /// Creates a new Merton jump-diffusion model.
    pub fn new(mu: f64, sigma: f64, lambda: f64, mu_j: f64, sigma_j: f64) -> Self {
        Self {
            mu,
            sigma,
            lambda,
            mu_j,
            sigma_j,
        }
    }
    /// Expected jump multiplier: E\[J\] = exp(mu_j + sigma_j^2/2) - 1.
    pub fn mean_jump_size(&self) -> f64 {
        (self.mu_j + 0.5 * self.sigma_j * self.sigma_j).exp() - 1.0
    }
    /// Simulates a path with jump-diffusion dynamics.
    ///
    /// Returns a vector of S values at each time step starting from `s0`.
    pub fn simulate(&self, s0: f64, t_end: f64, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut path = Vec::with_capacity(n);
        let mut s = s0;
        path.push(s);
        let m_j = self.mean_jump_size();
        let compensated_drift = self.mu - self.lambda * m_j;
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            let dw = rng.next_normal() * dt.sqrt();
            let n_jumps = rng.next_poisson(self.lambda * dt);
            let mut jump_mult = 1.0_f64;
            for _ in 0..n_jumps {
                let log_j = self.mu_j + self.sigma_j * rng.next_normal();
                jump_mult *= log_j.exp();
            }
            s = s * (compensated_drift * dt + self.sigma * dw).exp() * jump_mult;
            t += dt;
            path.push(s);
        }
        path
    }
}
/// Lévy flight: a random walk where step sizes follow a Lévy (stable) distribution.
///
/// Uses the Chambers-Mallows-Stuck algorithm to generate Lévy-stable variates
/// with stability index `alpha ∈ (0, 2]` and scale `c`.
#[allow(dead_code)]
pub struct LevyFlight {
    /// Stability index α ∈ (0, 2]. α=2 is Gaussian, α=1 is Cauchy.
    pub alpha: f64,
    /// Scale parameter c > 0.
    pub scale: f64,
}
impl LevyFlight {
    /// Create a new Lévy flight.
    pub fn new(alpha: f64, scale: f64) -> Self {
        Self { alpha, scale }
    }
    /// Sample a single Lévy-stable variate using the CMS algorithm.
    ///
    /// For α=2 this reduces to a Gaussian; for α=1 to a Cauchy distribution.
    pub fn sample_step(&self, rng: &mut Rng) -> f64 {
        let alpha = self.alpha.clamp(0.01, 2.0);
        if (alpha - 2.0).abs() < 1e-6 {
            return self.scale * rng.next_normal() * 2.0_f64.sqrt();
        }
        let u = PI * (rng.next_f64() - 0.5);
        let w = loop {
            let v = rng.next_f64();
            if v > 0.0 {
                break -v.ln();
            }
        };
        let zeta = (alpha * u).sin() / u.cos().powf(1.0 / alpha);
        let factor = (((1.0 - alpha) * u).cos() / w).powf((1.0 - alpha) / alpha);
        self.scale * zeta * factor
    }
    /// Generate a Lévy flight path of `n_steps` steps in 1-D.
    ///
    /// Returns a vector of cumulative positions.
    pub fn path(&self, n_steps: usize, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let mut path = Vec::with_capacity(n_steps + 1);
        let mut pos = 0.0_f64;
        path.push(pos);
        for _ in 0..n_steps {
            pos += self.sample_step(&mut rng);
            path.push(pos);
        }
        path
    }
}
/// Standard Brownian motion sampler using a `rand::Rng`.
///
/// Increments over each time step `dt` are drawn from N(0, dt).
pub struct WienerSampler {
    /// Time step size.
    pub dt: f64,
}
impl WienerSampler {
    /// Create a new `WienerSampler` with the given step size.
    pub fn new(dt: f64) -> Self {
        Self { dt }
    }
    /// Sample `n` independent Brownian increments ΔW ~ N(0, dt).
    pub fn sample(&self, n: usize, rng: &mut impl rand::Rng) -> Vec<f64> {
        let std = self.dt.sqrt();
        (0..n)
            .map(|_| {
                let u1: f64 = loop {
                    let v: f64 = rng.random();
                    if v > 0.0 {
                        break v;
                    }
                };
                let u2: f64 = rng.random();
                let mag = (-2.0 * u1.ln()).sqrt();
                mag * (2.0 * PI * u2).cos() * std
            })
            .collect()
    }
}
/// Metropolis-Hastings MCMC sampler.
///
/// Uses a symmetric Gaussian proposal `x_new ~ N(x, step_size^2)` and the
/// Metropolis acceptance criterion.
pub struct MetropolisHastings {
    /// Inverse temperature β = 1 / (k_B * T).
    pub beta: f64,
}
impl MetropolisHastings {
    /// Create a new `MetropolisHastings` sampler.
    pub fn new(beta: f64) -> Self {
        Self { beta }
    }
    /// Perform one Metropolis-Hastings step.
    ///
    /// Returns `(new_x, accepted)` where `accepted` indicates whether the
    /// proposed move was accepted.
    pub fn step<E: Fn(f64) -> f64>(
        &self,
        energy_fn: E,
        x: f64,
        step_size: f64,
        rng: &mut impl rand::Rng,
    ) -> (f64, bool) {
        let u1: f64 = loop {
            let v: f64 = rng.random();
            if v > 0.0 {
                break v;
            }
        };
        let u2: f64 = rng.random();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        let x_new = x + step_size * z;
        let delta_e = energy_fn(x_new) - energy_fn(x);
        let accept_prob = (-self.beta * delta_e).exp().min(1.0);
        let u: f64 = rng.random();
        if u < accept_prob {
            (x_new, true)
        } else {
            (x, false)
        }
    }
}
/// Cox-Ingersoll-Ross (CIR) square-root diffusion process.
///
/// dX = kappa*(theta - X)*dt + sigma*sqrt(X)*dW
///
/// Used for interest rate modeling and variance processes.
#[allow(dead_code)]
pub struct CirProcess {
    /// Mean-reversion speed κ.
    pub kappa: f64,
    /// Long-run mean θ.
    pub theta: f64,
    /// Volatility σ.
    pub sigma: f64,
}
impl CirProcess {
    /// Create a new CIR process.
    pub fn new(kappa: f64, theta: f64, sigma: f64) -> Self {
        Self {
            kappa,
            theta,
            sigma,
        }
    }
    /// Feller condition: 2*kappa*theta > sigma^2.
    ///
    /// If satisfied, the process never reaches zero.
    pub fn feller_satisfied(&self) -> bool {
        2.0 * self.kappa * self.theta > self.sigma * self.sigma
    }
    /// Simulate a CIR path using the Euler-Maruyama scheme.
    ///
    /// Uses full truncation to keep X ≥ 0.
    pub fn simulate(&self, x0: f64, t_end: f64, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut path = Vec::with_capacity(n);
        let mut x = x0.max(0.0);
        path.push(x);
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            let dw = rng.next_normal() * dt.sqrt();
            x += self.kappa * (self.theta - x) * dt + self.sigma * x.max(0.0).sqrt() * dw;
            x = x.max(0.0);
            t += dt;
            path.push(x);
        }
        path
    }
    /// Analytical stationary mean: theta.
    pub fn stationary_mean(&self) -> f64 {
        self.theta
    }
    /// Analytical stationary variance: theta * sigma^2 / (2 * kappa).
    pub fn stationary_variance(&self) -> f64 {
        self.theta * self.sigma * self.sigma / (2.0 * self.kappa)
    }
    /// Conditional mean at time t given X(0) = x0:
    /// E\[X(t)|X(0)\] = x0 * exp(-kappa*t) + theta * (1 - exp(-kappa*t))
    pub fn conditional_mean(&self, x0: f64, t: f64) -> f64 {
        let e = (-self.kappa * t).exp();
        x0 * e + self.theta * (1.0 - e)
    }
}
/// Standard Brownian motion (Wiener process) W(t).
///
/// Increments over a time step `dt` are drawn from N(0, dt).
pub struct WienerProcess {
    pub(super) rng: Rng,
}
impl WienerProcess {
    /// Creates a new Wiener process with the given seed.
    pub fn new(seed: u64) -> Self {
        Self {
            rng: Rng::new(seed),
        }
    }
    /// Returns a single Brownian increment dW ~ N(0, dt).
    pub fn increment(&mut self, dt: f64) -> f64 {
        self.rng.next_normal() * dt.sqrt()
    }
    /// Generates a sample path from `t0` to `t_end` with step size `dt`.
    ///
    /// Returns a `Vec` of `(time, value)` pairs starting at W(t0) = 0.
    pub fn path(&mut self, t0: f64, t_end: f64, dt: f64) -> Vec<(f64, f64)> {
        let n = ((t_end - t0) / dt).ceil() as usize + 1;
        let mut result = Vec::with_capacity(n);
        let mut t = t0;
        let mut w = 0.0_f64;
        result.push((t, w));
        while t + dt <= t_end + 1e-12 * dt {
            w += self.increment(dt);
            t += dt;
            result.push((t, w));
        }
        result
    }
    /// Generates a Brownian bridge from `(t0, w0)` to `(t_end, w_end)`.
    ///
    /// Returns a `Vec` of `(time, value)` pairs where the path is conditioned
    /// to pass through `w_end` at `t_end`.
    pub fn bridge(&mut self, t0: f64, t_end: f64, dt: f64, w0: f64, w_end: f64) -> Vec<(f64, f64)> {
        let n = ((t_end - t0) / dt).ceil() as usize + 1;
        let mut result = Vec::with_capacity(n);
        let total_t = t_end - t0;
        let free_path = self.path(t0, t_end, dt);
        let w_free_end = free_path.last().map(|&(_, w)| w).unwrap_or(0.0);
        for (i, &(ti, wi)) in free_path.iter().enumerate() {
            let frac = if total_t > 1e-15 {
                (ti - t0) / total_t
            } else {
                if i == 0 { 0.0 } else { 1.0 }
            };
            let bridge_val = wi - frac * w_free_end + w0 + frac * (w_end - w0);
            result.push((ti, bridge_val));
        }
        result
    }
}
/// Heston stochastic volatility model.
///
/// dS = mu*S*dt + sqrt(v)*S*dW_1
/// dv = kappa*(theta - v)*dt + xi*sqrt(v)*dW_2
/// Corr(dW_1, dW_2) = rho
pub struct HestonModel {
    /// Drift of the asset price.
    pub mu: f64,
    /// Mean-reversion speed of variance.
    pub kappa: f64,
    /// Long-term variance.
    pub theta: f64,
    /// Volatility of volatility.
    pub xi: f64,
    /// Correlation between price and variance Wiener processes.
    pub rho: f64,
}
impl HestonModel {
    /// Creates a new Heston model.
    pub fn new(mu: f64, kappa: f64, theta: f64, xi: f64, rho: f64) -> Self {
        Self {
            mu,
            kappa,
            theta,
            xi,
            rho,
        }
    }
    /// Checks the Feller condition: 2*kappa*theta > xi^2.
    ///
    /// If satisfied, the variance process never reaches zero.
    pub fn feller_satisfied(&self) -> bool {
        2.0 * self.kappa * self.theta > self.xi * self.xi
    }
    /// Simulates asset price and variance paths.
    ///
    /// Returns `(prices, variances)` vectors.
    #[allow(clippy::too_many_arguments)]
    pub fn simulate(
        &self,
        s0: f64,
        v0: f64,
        t_end: f64,
        dt: f64,
        seed: u64,
    ) -> (Vec<f64>, Vec<f64>) {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut prices = Vec::with_capacity(n);
        let mut variances = Vec::with_capacity(n);
        let mut s = s0;
        let mut v = v0;
        prices.push(s);
        variances.push(v);
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            let (z1, z2_raw) = rng.next_normal_pair();
            let dw1 = z1 * dt.sqrt();
            let dw2 = (self.rho * z1 + (1.0 - self.rho * self.rho).sqrt() * z2_raw) * dt.sqrt();
            let sqrt_v = v.max(0.0).sqrt();
            s += self.mu * s * dt + sqrt_v * s * dw1;
            v += self.kappa * (self.theta - v) * dt + self.xi * sqrt_v * dw2;
            v = v.max(0.0);
            t += dt;
            prices.push(s);
            variances.push(v);
        }
        (prices, variances)
    }
}
/// Fractional Brownian motion with Hurst exponent `H ∈ (0, 1)`.
///
/// Uses the Davies-Harte (circulant embedding) approximation via the
/// Cholesky decomposition of the covariance matrix for exact sampling.
/// For simplicity this implementation uses the Hosking method (exact covariance).
#[allow(non_snake_case)]
pub struct FractionalBrownianMotion {
    /// Hurst exponent: 0 < H < 1.  H=0.5 is standard Brownian motion.
    pub H: f64,
    /// Time step size.
    pub dt: f64,
}
impl FractionalBrownianMotion {
    /// Create a new `FractionalBrownianMotion`.
    #[allow(non_snake_case)]
    pub fn new(H: f64, dt: f64) -> Self {
        Self { H, dt }
    }
    /// Autocovariance of fBm increments at lag `k`:
    /// `γ(k) = 0.5 * (|k+1|^{2H} - 2|k|^{2H} + |k-1|^{2H}) * dt^{2H}`
    fn autocov(&self, k: usize) -> f64 {
        let h2 = 2.0 * self.H;
        let k = k as f64;
        0.5 * ((k + 1.0).powf(h2) - 2.0 * k.powf(h2) + (k - 1.0).abs().powf(h2)) * self.dt.powf(h2)
    }
    /// Sample `n` fBm increments using the Hosking method.
    ///
    /// This is an exact simulation algorithm for fractional Gaussian noise.
    pub fn sample(&self, n: usize, rng: &mut impl rand::Rng) -> Vec<f64> {
        if n == 0 {
            return vec![];
        }
        let std_normals: Vec<f64> = (0..n)
            .map(|_| {
                let u1: f64 = loop {
                    let v: f64 = rng.random();
                    if v > 0.0 {
                        break v;
                    }
                };
                let u2: f64 = rng.random();
                (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
            })
            .collect();
        let gamma0 = self.autocov(0);
        if gamma0 == 0.0 {
            return vec![0.0; n];
        }
        let mut x = vec![0.0_f64; n];
        let mut v = vec![0.0_f64; n];
        x[0] = gamma0.sqrt() * std_normals[0];
        v[0] = gamma0;
        for i in 1..n {
            let phi = self.autocov(i) / v[i - 1];
            v[i] = v[i - 1] * (1.0 - phi * phi);
            let mu = phi * x[i - 1];
            x[i] = mu + v[i].sqrt().max(0.0) * std_normals[i];
        }
        x
    }
}
/// Overdamped Langevin integrator.
///
/// Solves `dx/dt = F / (gamma * m) + sqrt(2 * kT / (gamma * m)) * ξ(t)`
/// where `ξ(t)` is Gaussian white noise.
#[allow(non_snake_case)]
pub struct LangevinIntegrator {
    /// Friction coefficient γ.
    pub gamma: f64,
    /// Thermal energy k_B * T.
    pub kT: f64,
    /// Particle mass.
    pub m: f64,
}
impl LangevinIntegrator {
    /// Create a new `LangevinIntegrator`.
    #[allow(non_snake_case)]
    pub fn new(gamma: f64, kT: f64, m: f64) -> Self {
        Self { gamma, kT, m }
    }
    /// Perform a single Euler-Maruyama step.
    ///
    /// Returns the new position.
    #[allow(non_snake_case)]
    pub fn step(&self, x: f64, F: f64, dt: f64, rng: &mut impl rand::Rng) -> f64 {
        let mobility = 1.0 / (self.gamma * self.m);
        let noise_amp = (2.0 * self.kT * mobility * dt).sqrt();
        let u1: f64 = loop {
            let v: f64 = rng.random();
            if v > 0.0 {
                break v;
            }
        };
        let u2: f64 = rng.random();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        x + mobility * F * dt + noise_amp * z
    }
}
/// Ornstein-Uhlenbeck (OU) process: mean-reverting stochastic process.
///
/// Equation: `dX = theta*(mu - X)*dt + sigma*dW`
///
/// Commonly used in physics for velocity autocorrelation and in finance
/// for mean-reverting interest rate models.
#[derive(Clone, Debug)]
pub struct OrnsteinUhlenbeck {
    /// Mean-reversion rate theta (1/time).
    pub theta: f64,
    /// Long-term mean mu.
    pub mu: f64,
    /// Volatility sigma.
    pub sigma: f64,
}
impl OrnsteinUhlenbeck {
    /// Creates a new Ornstein-Uhlenbeck process.
    pub fn new(theta: f64, mu: f64, sigma: f64) -> Self {
        Self { theta, mu, sigma }
    }
    /// Simulates a path using the exact discretization.
    ///
    /// Uses: `X_{n+1} = X_n * exp(-theta*dt) + mu*(1-exp(-theta*dt)) + sigma*sqrt((1-exp(-2*theta*dt))/(2*theta))*dW`
    ///
    /// Returns a vector of X values at each time step starting from `x0`.
    pub fn simulate(&self, x0: f64, t_end: f64, dt: f64, seed: u64) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let n = ((t_end / dt).ceil() as usize) + 1;
        let mut path = Vec::with_capacity(n);
        let mut x = x0;
        path.push(x);
        let e1 = (-self.theta * dt).exp();
        let std = self.sigma * ((1.0 - (-2.0 * self.theta * dt).exp()) / (2.0 * self.theta)).sqrt();
        let mut t = 0.0_f64;
        while t + dt <= t_end + 1e-12 * dt {
            x = x * e1 + self.mu * (1.0 - e1) + std * rng.next_normal();
            t += dt;
            path.push(x);
        }
        path
    }
    /// Returns the analytical stationary variance: sigma^2 / (2*theta).
    pub fn stationary_variance(&self) -> f64 {
        self.sigma * self.sigma / (2.0 * self.theta)
    }
    /// Returns the analytical autocorrelation at lag `tau`: exp(-theta*tau).
    pub fn autocorrelation_at(&self, tau: f64) -> f64 {
        (-self.theta * tau).exp()
    }
    /// Returns the analytical mean at time `t` starting from `x0`:
    /// `x0 * exp(-theta*t) + mu*(1 - exp(-theta*t))`
    pub fn analytical_mean_at(&self, x0: f64, t: f64) -> f64 {
        let e = (-self.theta * t).exp();
        x0 * e + self.mu * (1.0 - e)
    }
}
/// A fast linear congruential generator (LCG) for random number generation.
///
/// Uses Knuth's constants for good statistical properties. Suitable for
/// physics simulations where repeatability via seeding is important.
pub struct Rng {
    pub(super) state: u64,
}
impl Rng {
    /// Creates a new RNG with the given seed.
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    /// Advances the LCG state and returns the next raw 64-bit integer.
    ///
    /// Uses the recurrence: `state = state * 6364136223846793005 + 1442695040888963407`
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }
    /// Returns a uniformly distributed `f64` in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
    /// Returns a standard normal variate N(0, 1) using the Box-Muller transform.
    ///
    /// Generates a pair and discards the second value; use `next_normal_pair`
    /// when both samples are needed for efficiency.
    pub fn next_normal(&mut self) -> f64 {
        self.next_normal_pair().0
    }
    /// Returns two independent standard normal variates N(0, 1) using Box-Muller.
    ///
    /// Both values are independent draws from N(0, 1).
    pub fn next_normal_pair(&mut self) -> (f64, f64) {
        let u1 = loop {
            let v = self.next_f64();
            if v > 0.0 {
                break v;
            }
        };
        let u2 = self.next_f64();
        let mag = (-2.0 * u1.ln()).sqrt();
        let angle = 2.0 * PI * u2;
        (mag * angle.cos(), mag * angle.sin())
    }
    /// Returns an exponentially distributed variate with rate `lambda`.
    ///
    /// Uses the inverse transform: `-ln(U) / lambda`.
    pub fn next_exponential(&mut self, lambda: f64) -> f64 {
        let u = loop {
            let v = self.next_f64();
            if v > 0.0 {
                break v;
            }
        };
        -u.ln() / lambda
    }
    /// Returns a Poisson-distributed random integer with mean `lambda`.
    ///
    /// Uses Knuth's algorithm (suitable for small `lambda`).
    pub fn next_poisson(&mut self, lambda: f64) -> u64 {
        let l = (-lambda).exp();
        let mut k = 0u64;
        let mut p = 1.0_f64;
        loop {
            k += 1;
            p *= self.next_f64();
            if p <= l {
                break;
            }
        }
        k - 1
    }
}
/// A random walk in `dimension` dimensions with a fixed step size.
pub struct RandomWalk {
    /// Number of spatial dimensions.
    pub dimension: usize,
    /// Length of each step.
    pub step_size: f64,
}
impl RandomWalk {
    /// Creates a new random walk configuration.
    pub fn new(dimension: usize, step_size: f64) -> Self {
        Self {
            dimension,
            step_size,
        }
    }
    /// Generates a random walk path of `n_steps` steps starting from the origin.
    ///
    /// At each step, one random dimension is chosen and a +/-`step_size` displacement
    /// is applied. Returns a `Vec` of position vectors (each inner `Vec` has length
    /// equal to `dimension`).
    pub fn walk(&self, n_steps: usize, seed: u64) -> Vec<Vec<f64>> {
        let mut rng = Rng::new(seed);
        let mut positions = Vec::with_capacity(n_steps + 1);
        let mut pos = vec![0.0_f64; self.dimension];
        positions.push(pos.clone());
        for _ in 0..n_steps {
            let dim_idx = (rng.next_u64() as usize) % self.dimension;
            let sign = if rng.next_u64() & 1 == 0 { 1.0 } else { -1.0 };
            pos[dim_idx] += sign * self.step_size;
            positions.push(pos.clone());
        }
        positions
    }
    /// Generates a continuous random walk with Gaussian steps in all dimensions.
    ///
    /// Each step adds N(0, step_size) to every coordinate independently.
    pub fn gaussian_walk(&self, n_steps: usize, seed: u64) -> Vec<Vec<f64>> {
        let mut rng = Rng::new(seed);
        let mut positions = Vec::with_capacity(n_steps + 1);
        let mut pos = vec![0.0_f64; self.dimension];
        positions.push(pos.clone());
        for _ in 0..n_steps {
            for d in 0..self.dimension {
                pos[d] += rng.next_normal() * self.step_size;
            }
            positions.push(pos.clone());
        }
        positions
    }
    /// Computes the mean-squared displacement from a path.
    pub fn msd(path: &[Vec<f64>]) -> Vec<f64> {
        let n = path.len();
        if n == 0 {
            return vec![];
        }
        let dim = path[0].len();
        (0..n)
            .map(|i| {
                let mut sq = 0.0_f64;
                for d in 0..dim {
                    sq += path[i][d] * path[i][d];
                }
                sq
            })
            .collect()
    }
}
/// Geometric Brownian Motion with antithetic variates variance reduction.
pub struct AntitheticGbm;
impl AntitheticGbm {
    /// Generate `n_pairs` antithetic pairs of terminal GBM values.
    ///
    /// Each pair `(S+, S-)` uses `+Z` and `-Z` for the same Wiener increment.
    /// Returns a vector of `(s_plus, s_minus)` terminal values.
    pub fn generate_pairs(
        s0: f64,
        gbm: &GeometricBrownianMotion,
        t: f64,
        dt: f64,
        n_pairs: usize,
        seed: u64,
    ) -> Vec<(f64, f64)> {
        let n_steps = (t / dt).ceil() as usize;
        let mut pairs = Vec::with_capacity(n_pairs);
        let sqrt_dt = dt.sqrt();
        for p in 0..n_pairs {
            let mut rng = Rng::new(seed.wrapping_add(p as u64 * 2_654_435_761));
            let mut s_plus = s0;
            let mut s_minus = s0;
            for _ in 0..n_steps {
                let z = rng.next_normal();
                s_plus *=
                    ((gbm.mu - 0.5 * gbm.sigma * gbm.sigma) * dt + gbm.sigma * sqrt_dt * z).exp();
                s_minus *= ((gbm.mu - 0.5 * gbm.sigma * gbm.sigma) * dt
                    + gbm.sigma * sqrt_dt * (-z))
                    .exp();
            }
            pairs.push((s_plus, s_minus));
        }
        pairs
    }
    /// Estimate `E[f(S_T)]` using antithetic variates.
    ///
    /// Returns `(estimate, standard_error)`.
    pub fn estimate_mean(
        s0: f64,
        gbm: &GeometricBrownianMotion,
        t: f64,
        dt: f64,
        n_pairs: usize,
        seed: u64,
    ) -> (f64, f64) {
        let pairs = Self::generate_pairs(s0, gbm, t, dt, n_pairs, seed);
        let estimates: Vec<f64> = pairs.iter().map(|(a, b)| 0.5 * (a + b)).collect();
        let n = estimates.len() as f64;
        let mean = estimates.iter().sum::<f64>() / n;
        let var: f64 =
            estimates.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);
        (mean, (var / n).sqrt())
    }
}
