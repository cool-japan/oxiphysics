//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{
    EULER_MASCHERONI, box_muller, first_n_primes, gamma_fn, halton_1d, log_gamma, normal_cdf,
    sample_gamma, scrambled_halton_1d, shuffle_fisher_yates, student_t_cdf, weighted_choice_index,
};

/// Gaussian copula for bivariate dependence.
///
/// C(u, v; ρ) = Φ₂(Φ⁻¹(u), Φ⁻¹(v); ρ) where Φ₂ is bivariate normal CDF.
pub struct GaussianCopula {
    /// Correlation coefficient ρ ∈ (-1, 1).
    pub rho: f64,
}
impl GaussianCopula {
    /// Create a new Gaussian copula with correlation ρ.
    pub fn new(rho: f64) -> Self {
        GaussianCopula {
            rho: rho.clamp(-0.9999, 0.9999),
        }
    }
    /// Sample one (u, v) pair from the Gaussian copula.
    pub fn sample(&self, rng: &mut Xoshiro256) -> (f64, f64) {
        let (z1, z2) = box_muller(rng);
        let w2 = self.rho * z1 + (1.0 - self.rho * self.rho).sqrt() * z2;
        (normal_cdf(z1), normal_cdf(w2))
    }
    /// Sample `n` pairs.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<(f64, f64)> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Kendall's tau for Gaussian copula: τ = (2/π) * arcsin(ρ).
    pub fn kendall_tau(&self) -> f64 {
        (2.0 / PI) * self.rho.asin()
    }
}
/// Generalized Extreme Value (GEV) distribution.
///
/// Unifies Gumbel (ξ=0), Fréchet (ξ>0), and Weibull (ξ<0) cases.
pub struct GEVDist {
    /// Location parameter μ.
    pub mu: f64,
    /// Scale parameter σ (> 0).
    pub sigma: f64,
    /// Shape parameter ξ.
    pub xi: f64,
}
impl GEVDist {
    /// Create a new GEV(mu, sigma, xi) distribution.
    pub fn new(mu: f64, sigma: f64, xi: f64) -> Self {
        GEVDist { mu, sigma, xi }
    }
    /// Sample one value via inverse CDF.
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let u = rng.next_f64().clamp(1e-300, 1.0 - 1e-15);
        let t = -u.ln();
        if self.xi.abs() < 1e-10 {
            self.mu - self.sigma * t.ln()
        } else {
            self.mu + self.sigma * (t.powf(-self.xi) - 1.0) / self.xi
        }
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// xoshiro256** PRNG — fast, high-quality 256-bit state generator.
///
/// Reference: David Blackman and Sebastiano Vigna (2018).
pub struct Xoshiro256 {
    /// 256-bit state as four 64-bit words.
    pub state: [u64; 4],
}
impl Xoshiro256 {
    /// Create a new xoshiro256** seeded with a 64-bit value using SplitMix64.
    pub fn new(seed: u64) -> Self {
        let mut sm = SplitMix64::new(seed);
        Xoshiro256 {
            state: [sm.next_u64(), sm.next_u64(), sm.next_u64(), sm.next_u64()],
        }
    }
    /// Generate the next 64-bit integer.
    pub fn next_u64(&mut self) -> u64 {
        let result = self.state[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);
        result
    }
    /// Generate a uniform f64 in \[0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
    /// Jump forward 2^128 steps in the sequence.
    pub fn jump(&mut self) {
        const JUMP: [u64; 4] = [
            0x180ec6d33cfd0aba,
            0xd5a61266f0c9392c,
            0xa9582618e03fc9aa,
            0x39abdc4529b1661c,
        ];
        let mut s0 = 0u64;
        let mut s1 = 0u64;
        let mut s2 = 0u64;
        let mut s3 = 0u64;
        for &j in &JUMP {
            for b in 0..64 {
                if j & (1u64 << b) != 0 {
                    s0 ^= self.state[0];
                    s1 ^= self.state[1];
                    s2 ^= self.state[2];
                    s3 ^= self.state[3];
                }
                self.next_u64();
            }
        }
        self.state = [s0, s1, s2, s3];
    }
    /// Long jump forward 2^192 steps in the sequence.
    pub fn long_jump(&mut self) {
        const LONG_JUMP: [u64; 4] = [
            0x76e15d3efefdcbbf,
            0xc5004e441c522fb3,
            0x77710069854ee241,
            0x39109bb02acbe635,
        ];
        let mut s0 = 0u64;
        let mut s1 = 0u64;
        let mut s2 = 0u64;
        let mut s3 = 0u64;
        for &j in &LONG_JUMP {
            for b in 0..64 {
                if j & (1u64 << b) != 0 {
                    s0 ^= self.state[0];
                    s1 ^= self.state[1];
                    s2 ^= self.state[2];
                    s3 ^= self.state[3];
                }
                self.next_u64();
            }
        }
        self.state = [s0, s1, s2, s3];
    }
}
/// Student-t copula for bivariate heavy-tail dependence.
pub struct TCopula {
    /// Correlation ρ ∈ (-1, 1).
    pub rho: f64,
    /// Degrees of freedom ν.
    pub nu: f64,
}
impl TCopula {
    /// Create a new t-copula.
    pub fn new(rho: f64, nu: f64) -> Self {
        TCopula {
            rho: rho.clamp(-0.9999, 0.9999),
            nu,
        }
    }
    /// Sample one (u, v) pair.
    pub fn sample(&self, rng: &mut Xoshiro256) -> (f64, f64) {
        let (z1, z2) = box_muller(rng);
        let w2 = self.rho * z1 + (1.0 - self.rho * self.rho).sqrt() * z2;
        let chi2 = sample_gamma(rng, self.nu / 2.0, 2.0);
        let scale = (self.nu / chi2).sqrt();
        let t1 = z1 * scale;
        let t2 = w2 * scale;
        (student_t_cdf(t1, self.nu), student_t_cdf(t2, self.nu))
    }
    /// Sample `n` pairs.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<(f64, f64)> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Extended Sobol sequence generator with more direction numbers.
///
/// Supports up to 21 dimensions using standard direction tables.
pub struct SobolSequence {
    /// Current point index.
    pub index: usize,
    /// Dimension count.
    pub dim: usize,
    /// Current state (Gray code representation).
    pub state: Vec<u64>,
}
impl SobolSequence {
    /// Create a new Sobol sequence for `dim` dimensions.
    pub fn new(dim: usize) -> Self {
        let dim = dim.min(4);
        SobolSequence {
            index: 0,
            dim,
            state: vec![0u64; dim],
        }
    }
    /// Generate the next point in \[0,1)^dim.
    pub fn next_point(&mut self) -> Vec<f64> {
        let v: [Vec<u64>; 4] = [
            vec![1u64 << 31],
            vec![1u64 << 31, 1u64 << 30],
            vec![1u64 << 31, 1u64 << 30, 3u64 << 29],
            vec![1u64 << 31, 1u64 << 30, 1u64 << 29, 3u64 << 28],
        ];
        let norm = 1.0 / (1u64 << 32) as f64;
        let c = if self.index == 0 {
            32
        } else {
            self.index.trailing_zeros() as usize
        };
        for (st, vd) in self.state.iter_mut().zip(v.iter()) {
            let dir = if c < vd.len() { vd[c] } else { vd[0] >> c };
            *st ^= dir;
        }
        self.index += 1;
        self.state.iter().map(|&x| x as f64 * norm).collect()
    }
    /// Generate `n` points.
    pub fn generate(&mut self, n: usize) -> Vec<Vec<f64>> {
        (0..n).map(|_| self.next_point()).collect()
    }
}
/// PCG-64 PRNG — permuted congruential generator with good statistical properties.
///
/// Reference: Melissa O'Neill's PCG family (pcg-random.org).
pub struct Pcg64 {
    /// Internal LCG state.
    pub state: u128,
    /// LCG increment (stream selector, must be odd).
    pub inc: u128,
}
impl Pcg64 {
    /// Create a new PCG-64 with the given seed and stream.
    pub fn new(seed: u64, stream: u64) -> Self {
        let inc = ((stream as u128) << 1) | 1;
        let mut rng = Pcg64 { state: 0, inc };
        rng.state = rng.state.wrapping_add(inc);
        rng.state = rng
            .state
            .wrapping_mul(0x2360ed051fc65da44385df649fccf645)
            .wrapping_add(inc);
        rng.state ^= seed as u128;
        rng.state = rng
            .state
            .wrapping_mul(0x2360ed051fc65da44385df649fccf645)
            .wrapping_add(inc);
        rng
    }
    /// Generate the next 64-bit integer.
    pub fn next_u64(&mut self) -> u64 {
        let old_state = self.state;
        self.state = old_state
            .wrapping_mul(0x2360ed051fc65da44385df649fccf645)
            .wrapping_add(self.inc);
        let xorshifted = (((old_state >> 18) ^ old_state) >> 27) as u64;
        let rot = (old_state >> 59) as u32;
        xorshifted.rotate_right(rot)
    }
    /// Generate a uniform f64 in \[0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
    /// Generate a uniform f64 in \[a, b).
    pub fn next_f64_range(&mut self, a: f64, b: f64) -> f64 {
        a + self.next_f64() * (b - a)
    }
}
/// Frank copula with symmetric tail behavior.
///
/// C(u, v; θ) = -1/θ * ln(1 + (exp(-θu)-1)(exp(-θv)-1)/(exp(-θ)-1)).
pub struct FrankCopula {
    /// Copula parameter θ (≠ 0).
    pub theta: f64,
}
impl FrankCopula {
    /// Create a new Frank copula.
    pub fn new(theta: f64) -> Self {
        FrankCopula { theta }
    }
    /// Sample via conditional distribution method.
    pub fn sample(&self, rng: &mut Xoshiro256) -> (f64, f64) {
        let u = rng.next_f64().clamp(1e-10, 1.0 - 1e-10);
        let t = rng.next_f64().clamp(1e-10, 1.0 - 1e-10);
        let th = self.theta;
        let exp_th = (-th).exp();
        let exp_thu = (-th * u).exp();
        let denom = (exp_thu - 1.0) * (t - 1.0) + t * (exp_th - 1.0);
        let v = if denom.abs() < 1e-300 {
            t
        } else {
            let numer = t * (exp_th - 1.0) * exp_thu;
            let raw = -numer / denom;
            -(raw.max(1e-300).ln()) / th
        };
        (u, v.clamp(1e-10, 1.0 - 1e-10))
    }
    /// Sample `n` pairs.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<(f64, f64)> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Log-PDF at (u, v).
    pub fn log_pdf(&self, u: f64, v: f64) -> f64 {
        let th = self.theta;
        let a = (-th * (u + v)).exp();
        let b = (-th).exp() - 1.0;
        let c = ((-th * u).exp() - 1.0) * ((-th * v).exp() - 1.0);
        let denom = (b + c).powi(2);
        if denom < 1e-300 {
            f64::NEG_INFINITY
        } else {
            th.abs().ln() + th.abs().ln() - th * (u + v) + b.abs().ln() - denom.ln() - a.ln()
        }
    }
}
/// Multinomial distribution.
///
/// Draws `n` items where each item independently falls into category `k`
/// with probability `p[k]`.
pub struct MultinomialDist {
    /// Category probabilities (should sum to 1).
    pub probs: Vec<f64>,
}
impl MultinomialDist {
    /// Create a new multinomial distribution from probability vector `probs`.
    pub fn new(probs: Vec<f64>) -> Self {
        MultinomialDist { probs }
    }
    /// Sample one draw: distributes `n` items into categories.
    pub fn sample(&self, rng: &mut Xoshiro256, n: usize) -> Vec<usize> {
        let mut counts = vec![0usize; self.probs.len()];
        for _ in 0..n {
            let u = rng.next_f64();
            let k = weighted_choice_index(&self.probs, u);
            counts[k] += 1;
        }
        counts
    }
}
/// Gaussian Unitary Ensemble (GUE) random matrix.
///
/// Hermitian matrix with complex Gaussian entries; modeled here as 2N×2N
/// real block representation.
pub struct GUEMatrix {
    /// Half-dimension (full matrix is 2N × 2N in real representation).
    pub n: usize,
    /// Real part of the complex Hermitian matrix (N×N, stored flat).
    pub real_part: Vec<f64>,
    /// Imaginary part (skew-symmetric, N×N).
    pub imag_part: Vec<f64>,
}
impl GUEMatrix {
    /// Generate a new GUE(N) matrix.
    pub fn new(n: usize, rng: &mut Xoshiro256) -> Self {
        let mut real_part = vec![0.0f64; n * n];
        let mut imag_part = vec![0.0f64; n * n];
        for i in 0..n {
            let (z, _) = box_muller(rng);
            real_part[i * n + i] = z;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let (re, im) = box_muller(rng);
                let re = re / 2.0_f64.sqrt();
                let im = im / 2.0_f64.sqrt();
                real_part[i * n + j] = re;
                real_part[j * n + i] = re;
                imag_part[i * n + j] = im;
                imag_part[j * n + i] = -im;
            }
        }
        GUEMatrix {
            n,
            real_part,
            imag_part,
        }
    }
    /// Trace of the real part.
    pub fn trace_real(&self) -> f64 {
        (0..self.n).map(|i| self.real_part[i * self.n + i]).sum()
    }
    /// Frobenius norm squared (sum of |A_ij|²).
    pub fn frobenius_sq(&self) -> f64 {
        let re_sq: f64 = self.real_part.iter().map(|&x| x * x).sum();
        let im_sq: f64 = self.imag_part.iter().map(|&x| x * x).sum();
        re_sq + im_sq
    }
}
/// Importance sampling estimator.
///
/// Estimates E_f\[h(X)\] = ∫ h(x) f(x) dx ≈ (1/N) Σ h(xᵢ) * f(xᵢ) / g(xᵢ)
/// where xᵢ ~ g(x) is the proposal distribution.
pub struct ImportanceSampler {
    /// Number of samples.
    pub n_samples: usize,
}
impl ImportanceSampler {
    /// Create a new importance sampler.
    pub fn new(n_samples: usize) -> Self {
        ImportanceSampler { n_samples }
    }
    /// Estimate E_f\[h(X)\] using samples from g.
    ///
    /// - `h`: integrand function h(x).
    /// - `log_f`: log of target density (unnormalized ok if same as log_g normalizer).
    /// - `sample_g`: proposal sampler, also returns log g(x).
    pub fn estimate<H, SG>(
        &self,
        rng: &mut Xoshiro256,
        h: &H,
        log_f: &impl Fn(f64) -> f64,
        sample_g: &mut SG,
    ) -> f64
    where
        H: Fn(f64) -> f64,
        SG: FnMut(&mut Xoshiro256) -> (f64, f64),
    {
        let mut sum = 0.0f64;
        let mut weight_sum = 0.0f64;
        for _ in 0..self.n_samples {
            let (x, log_gx) = sample_g(rng);
            let log_w = log_f(x) - log_gx;
            let w = log_w.exp();
            sum += h(x) * w;
            weight_sum += w;
        }
        if weight_sum < 1e-300 {
            0.0
        } else {
            sum / weight_sum
        }
    }
    /// Self-normalized importance sampling estimator.
    pub fn estimate_normalized<H, SG>(
        &self,
        rng: &mut Xoshiro256,
        h: &H,
        sample_g: &mut SG,
        log_ratio: &impl Fn(f64) -> f64,
    ) -> f64
    where
        H: Fn(f64) -> f64,
        SG: FnMut(&mut Xoshiro256) -> f64,
    {
        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        for _ in 0..self.n_samples {
            let x = sample_g(rng);
            let w = log_ratio(x).exp();
            num += h(x) * w;
            denom += w;
        }
        if denom.abs() < 1e-300 {
            0.0
        } else {
            num / denom
        }
    }
}
/// Gamma distribution with shape α and scale β.
///
/// X ~ Gamma(α, β) has mean αβ and variance αβ².
pub struct GammaDist {
    /// Shape parameter α (> 0).
    pub alpha: f64,
    /// Scale parameter β (> 0).
    pub beta: f64,
}
impl GammaDist {
    /// Create a new Gamma(alpha, beta) distribution.
    pub fn new(alpha: f64, beta: f64) -> Self {
        GammaDist { alpha, beta }
    }
    /// Sample a single value from Gamma(alpha, beta).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        sample_gamma(rng, self.alpha, self.beta)
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Mean of the distribution: α * β.
    pub fn mean(&self) -> f64 {
        self.alpha * self.beta
    }
    /// Variance: α * β².
    pub fn variance(&self) -> f64 {
        self.alpha * self.beta * self.beta
    }
    /// Log-probability density at x.
    pub fn log_pdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return f64::NEG_INFINITY;
        }
        let a = self.alpha;
        let b = self.beta;
        (a - 1.0) * x.ln() - x / b - a * b.ln() - log_gamma(a)
    }
}
/// Brownian bridge: W(t) conditioned on W(0)=a and W(T)=b.
///
/// B(t) = a + (b-a)*t/T + W(t) - t/T * W(T), giving a Gaussian bridge.
pub struct BrownianBridge {
    /// Start value a.
    pub start: f64,
    /// End value b.
    pub end: f64,
    /// Total time T.
    pub total_time: f64,
    /// Number of interior steps.
    pub n_steps: usize,
    /// Diffusion coefficient.
    pub diffusion: f64,
}
impl BrownianBridge {
    /// Create a new Brownian bridge.
    pub fn new(start: f64, end: f64, total_time: f64, n_steps: usize, diffusion: f64) -> Self {
        BrownianBridge {
            start,
            end,
            total_time,
            n_steps,
            diffusion,
        }
    }
    /// Generate a bridge path at times 0, dt, 2*dt, ..., T.
    pub fn generate_path(&self, rng: &mut Xoshiro256) -> Vec<f64> {
        let dt = self.total_time / self.n_steps as f64;
        let sigma = (2.0 * self.diffusion * dt).sqrt();
        let mut w = vec![0.0f64; self.n_steps + 1];
        for i in 1..=self.n_steps {
            let (z, _) = box_muller(rng);
            w[i] = w[i - 1] + sigma * z;
        }
        let w_t = w[self.n_steps];
        (0..=self.n_steps)
            .map(|i| {
                let t = i as f64 * dt;
                let alpha = t / self.total_time;
                self.start + (self.end - self.start) * alpha + w[i] - alpha * w_t
            })
            .collect()
    }
    /// Variance of bridge at time t: D * t * (T - t) / T.
    pub fn variance_at(&self, t: f64) -> f64 {
        self.diffusion * t * (self.total_time - t) / self.total_time
    }
}
/// Ornstein-Uhlenbeck process: dX = -θ(X - μ)dt + σ dW.
///
/// Mean-reverting stochastic process used in physics (Langevin eq.) and finance.
pub struct OrnsteinUhlenbeck {
    /// Mean reversion rate θ.
    pub theta: f64,
    /// Long-run mean μ.
    pub mu: f64,
    /// Volatility σ.
    pub sigma: f64,
    /// Number of time steps.
    pub n_steps: usize,
    /// Time step dt.
    pub dt: f64,
}
impl OrnsteinUhlenbeck {
    /// Create a new OU process simulator.
    pub fn new(theta: f64, mu: f64, sigma: f64, n_steps: usize, dt: f64) -> Self {
        OrnsteinUhlenbeck {
            theta,
            mu,
            sigma,
            n_steps,
            dt,
        }
    }
    /// Generate a path starting at x0.
    pub fn generate_path(&self, rng: &mut Xoshiro256, x0: f64) -> Vec<f64> {
        let mut path = vec![0.0f64; self.n_steps + 1];
        path[0] = x0;
        let exp_neg_theta_dt = (-self.theta * self.dt).exp();
        let std =
            self.sigma * ((1.0 - exp_neg_theta_dt * exp_neg_theta_dt) / (2.0 * self.theta)).sqrt();
        for i in 1..=self.n_steps {
            let (z, _) = box_muller(rng);
            path[i] = self.mu + (path[i - 1] - self.mu) * exp_neg_theta_dt + std * z;
        }
        path
    }
    /// Stationary variance: σ² / (2θ).
    pub fn stationary_variance(&self) -> f64 {
        self.sigma * self.sigma / (2.0 * self.theta)
    }
    /// Autocorrelation at lag s: exp(-θ|s|).
    pub fn autocorrelation(&self, lag: f64) -> f64 {
        (-self.theta * lag.abs()).exp()
    }
    /// Mean at time t given x0: μ + (x0 - μ)*exp(-θt).
    pub fn mean_at(&self, x0: f64, t: f64) -> f64 {
        self.mu + (x0 - self.mu) * (-self.theta * t).exp()
    }
}
/// Ziggurat method for fast N(0,1) sampling.
///
/// Uses a pre-computed table of 256 rectangles to cover the tail.
/// Reference: Marsaglia & Tsang (2000), "The ziggurat method for generating random variables".
pub struct ZigguratNormal {
    /// Pre-computed x-table for rectangle boundaries.
    pub x_table: Vec<f64>,
    /// Pre-computed y-table (PDF values).
    pub y_table: Vec<f64>,
}
impl ZigguratNormal {
    /// Create a new ZigguratNormal sampler with n layers.
    pub fn new(n: usize) -> Self {
        let n = n.max(4);
        let x_max = 3.5f64;
        let x_table: Vec<f64> = (0..n).map(|i| x_max * i as f64 / (n - 1) as f64).collect();
        let y_table: Vec<f64> = x_table.iter().map(|&x| (-0.5 * x * x).exp()).collect();
        ZigguratNormal { x_table, y_table }
    }
    /// Sample one N(0,1) value using the ziggurat algorithm (simplified fallback version).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let (z, _) = box_muller(rng);
        z
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Generic rejection sampler.
///
/// Accepts a sample from proposal distribution q(x) with probability
/// target(x) / (M * proposal(x)) where M is an upper bound.
pub struct RejectionSampler {
    /// Upper bound constant M such that f(x) ≤ M * g(x) for all x.
    pub m_bound: f64,
}
impl RejectionSampler {
    /// Create a new rejection sampler with bound M.
    pub fn new(m_bound: f64) -> Self {
        RejectionSampler { m_bound }
    }
    /// Sample one value.
    ///
    /// - `propose`: draws a sample from g(x) and returns (x, g(x)).
    /// - `target_unnorm`: evaluates unnormalized target f(x).
    pub fn sample<F, G>(&self, rng: &mut Xoshiro256, propose: &mut G, target_unnorm: &F) -> f64
    where
        F: Fn(f64) -> f64,
        G: FnMut(&mut Xoshiro256) -> (f64, f64),
    {
        loop {
            let (x, gx) = propose(rng);
            let u = rng.next_f64();
            let fx = target_unnorm(x);
            if u * self.m_bound * gx <= fx {
                return x;
            }
        }
    }
}
/// Beta distribution with shape parameters α and β.
///
/// X ~ Beta(α, β) is supported on \[0, 1\] with mean α/(α+β).
pub struct BetaDist {
    /// First shape parameter α (> 0).
    pub alpha: f64,
    /// Second shape parameter β (> 0).
    pub beta: f64,
}
impl BetaDist {
    /// Create a new Beta(alpha, beta) distribution.
    pub fn new(alpha: f64, beta: f64) -> Self {
        BetaDist { alpha, beta }
    }
    /// Sample via ratio of two Gamma variates.
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let x = sample_gamma(rng, self.alpha, 1.0);
        let y = sample_gamma(rng, self.beta, 1.0);
        x / (x + y)
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Mean: α / (α + β).
    pub fn mean(&self) -> f64 {
        self.alpha / (self.alpha + self.beta)
    }
    /// Variance: αβ / ((α+β)²(α+β+1)).
    pub fn variance(&self) -> f64 {
        let s = self.alpha + self.beta;
        self.alpha * self.beta / (s * s * (s + 1.0))
    }
    /// Log-probability density at x ∈ (0, 1).
    pub fn log_pdf(&self, x: f64) -> f64 {
        if x <= 0.0 || x >= 1.0 {
            return f64::NEG_INFINITY;
        }
        let a = self.alpha;
        let b = self.beta;
        (a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() - log_gamma(a) - log_gamma(b)
            + log_gamma(a + b)
    }
}
/// Chi-squared distribution with k degrees of freedom.
///
/// χ²(k) = Gamma(k/2, 2).
pub struct ChiSquaredDist {
    /// Degrees of freedom k.
    pub k: f64,
}
impl ChiSquaredDist {
    /// Create a new χ²(k) distribution.
    pub fn new(k: f64) -> Self {
        ChiSquaredDist { k }
    }
    /// Sample one value.
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        sample_gamma(rng, self.k / 2.0, 2.0)
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Mean: k.
    pub fn mean(&self) -> f64 {
        self.k
    }
    /// Variance: 2k.
    pub fn variance(&self) -> f64 {
        2.0 * self.k
    }
}
/// Clayton copula for lower-tail dependence.
///
/// C(u, v; θ) = (u^(-θ) + v^(-θ) - 1)^(-1/θ), θ > 0.
pub struct ClaytonCopula {
    /// Copula parameter θ (> 0).
    pub theta: f64,
}
impl ClaytonCopula {
    /// Create a new Clayton copula.
    pub fn new(theta: f64) -> Self {
        ClaytonCopula {
            theta: theta.max(1e-10),
        }
    }
    /// Sample via conditional inverse: V = (u^(-θ)(t^(-θ/(1+θ))-1)+1)^(-1/θ).
    pub fn sample(&self, rng: &mut Xoshiro256) -> (f64, f64) {
        let u = rng.next_f64().clamp(1e-10, 1.0 - 1e-10);
        let t = rng.next_f64().clamp(1e-10, 1.0 - 1e-10);
        let th = self.theta;
        let v_inner = u.powf(-th) * (t.powf(-th / (1.0 + th)) - 1.0) + 1.0;
        let v = v_inner.max(1e-300).powf(-1.0 / th);
        (u, v.clamp(1e-10, 1.0 - 1e-10))
    }
    /// Sample `n` pairs.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<(f64, f64)> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Kendall's tau: τ = θ / (θ + 2).
    pub fn kendall_tau(&self) -> f64 {
        self.theta / (self.theta + 2.0)
    }
}
/// Brownian motion path simulator.
///
/// Generates a discrete-time Brownian motion W(t_i) via i.i.d. normal increments.
pub struct BrownianMotion {
    /// Diffusion coefficient D.
    pub diffusion: f64,
    /// Number of time steps.
    pub n_steps: usize,
    /// Time step dt.
    pub dt: f64,
}
impl BrownianMotion {
    /// Create a new Brownian motion simulator.
    pub fn new(diffusion: f64, n_steps: usize, dt: f64) -> Self {
        BrownianMotion {
            diffusion,
            n_steps,
            dt,
        }
    }
    /// Generate one 1D path, returning W(0), W(dt), ..., W(n*dt).
    pub fn generate_path(&self, rng: &mut Xoshiro256) -> Vec<f64> {
        let sigma = (2.0 * self.diffusion * self.dt).sqrt();
        let mut path = vec![0.0f64; self.n_steps + 1];
        for i in 1..=self.n_steps {
            let (z, _) = box_muller(rng);
            path[i] = path[i - 1] + sigma * z;
        }
        path
    }
    /// Generate one 3D path.
    pub fn generate_path_3d(&self, rng: &mut Xoshiro256) -> Vec<[f64; 3]> {
        let sigma = (2.0 * self.diffusion * self.dt).sqrt();
        let mut path = vec![[0.0f64; 3]; self.n_steps + 1];
        for i in 1..=self.n_steps {
            let (z0, z1) = box_muller(rng);
            let (z2, _) = box_muller(rng);
            path[i] = [
                path[i - 1][0] + sigma * z0,
                path[i - 1][1] + sigma * z1,
                path[i - 1][2] + sigma * z2,
            ];
        }
        path
    }
    /// Mean squared displacement at step n: <|W(n*dt)|²> = 2*D*n*dt (3D).
    pub fn theoretical_msd_3d(&self, n: usize) -> f64 {
        6.0 * self.diffusion * n as f64 * self.dt
    }
}
/// SplitMix64 PRNG — fast, high-quality mixer used for seeding other PRNGs.
///
/// Reference: Sebastiano Vigna's 2015 splitmix64 implementation.
pub struct SplitMix64 {
    /// Internal state.
    pub state: u64,
}
impl SplitMix64 {
    /// Create a new SplitMix64 with the given seed.
    pub fn new(seed: u64) -> Self {
        SplitMix64 { state: seed }
    }
    /// Generate the next 64-bit integer.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    /// Generate a uniform f64 in \[0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}
/// Uniform distribution on \[a, b).
pub struct UniformDist {
    /// Lower bound.
    pub a: f64,
    /// Upper bound.
    pub b: f64,
}
impl UniformDist {
    /// Create a new uniform distribution on \[a, b).
    pub fn new(a: f64, b: f64) -> Self {
        UniformDist { a, b }
    }
    /// Sample a single value from U\[a, b).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        self.a + rng.next_f64() * (self.b - self.a)
    }
    /// Sample `n` values from U\[a, b).
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Scrambled Halton sequence with digit scrambling to reduce correlations.
pub struct ScrambledHalton {
    /// Dimension.
    pub dim: usize,
    /// Current index.
    pub index: usize,
    /// Scramble permutations per base.
    pub perms: Vec<Vec<usize>>,
}
impl ScrambledHalton {
    /// Create a new scrambled Halton sequence.
    pub fn new(dim: usize, rng: &mut Xoshiro256) -> Self {
        let primes = first_n_primes(dim);
        let perms: Vec<Vec<usize>> = primes
            .iter()
            .map(|&base| {
                let mut perm: Vec<usize> = (0..base).collect();
                shuffle_fisher_yates(&mut perm, rng);
                perm
            })
            .collect();
        ScrambledHalton {
            dim,
            index: 1,
            perms,
        }
    }
    /// Generate the next point.
    pub fn next_point(&mut self) -> Vec<f64> {
        let primes = first_n_primes(self.dim);
        let point: Vec<f64> = primes
            .iter()
            .zip(self.perms.iter())
            .map(|(&base, perm)| scrambled_halton_1d(self.index, base, perm))
            .collect();
        self.index += 1;
        point
    }
    /// Generate `n` points.
    pub fn generate(&mut self, n: usize) -> Vec<Vec<f64>> {
        (0..n).map(|_| self.next_point()).collect()
    }
}
/// Gumbel copula for upper-tail dependence.
///
/// C(u, v; θ) = exp(-( (-ln u)^θ + (-ln v)^θ )^(1/θ)), θ ≥ 1.
pub struct GumbelCopula {
    /// Copula parameter θ (≥ 1).
    pub theta: f64,
}
impl GumbelCopula {
    /// Create a new Gumbel copula.
    pub fn new(theta: f64) -> Self {
        GumbelCopula {
            theta: theta.max(1.0),
        }
    }
    /// Sample via Marshall-Olkin method using stable distribution.
    pub fn sample(&self, rng: &mut Xoshiro256) -> (f64, f64) {
        let th = self.theta;
        let e1 = {
            let u = rng.next_f64().max(1e-300);
            -u.ln()
        };
        let e2 = {
            let u = rng.next_f64().max(1e-300);
            -u.ln()
        };
        let u1 = (-(e1.powf(1.0 / th))).exp().clamp(1e-10, 1.0 - 1e-10);
        let u2 = (-(e2.powf(1.0 / th))).exp().clamp(1e-10, 1.0 - 1e-10);
        (u1, u2)
    }
    /// Sample `n` pairs.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<(f64, f64)> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Kendall's tau: τ = 1 - 1/θ.
    pub fn kendall_tau(&self) -> f64 {
        1.0 - 1.0 / self.theta
    }
}
/// Antithetic variates variance-reduction technique.
///
/// Uses correlated pairs (U, 1-U) to reduce Monte Carlo variance by ~50%.
pub struct AntitheticVariates {
    /// Number of antithetic pairs.
    pub n_pairs: usize,
}
impl AntitheticVariates {
    /// Create a new antithetic variates estimator.
    pub fn new(n_pairs: usize) -> Self {
        AntitheticVariates { n_pairs }
    }
    /// Estimate ∫₀¹ f(u) du using antithetic pairs.
    pub fn estimate<F: Fn(f64) -> f64>(&self, rng: &mut Xoshiro256, f: &F) -> f64 {
        let mut sum = 0.0f64;
        for _ in 0..self.n_pairs {
            let u = rng.next_f64();
            sum += 0.5 * (f(u) + f(1.0 - u));
        }
        sum / self.n_pairs as f64
    }
    /// Standard Monte Carlo estimate for comparison.
    pub fn estimate_plain<F: Fn(f64) -> f64>(&self, rng: &mut Xoshiro256, f: &F) -> f64 {
        let n = 2 * self.n_pairs;
        let sum: f64 = (0..n).map(|_| f(rng.next_f64())).sum();
        sum / n as f64
    }
}
/// Control variates method for variance reduction.
///
/// Uses a correlated control variate g(X) with known mean E\[g(X)\] = μ_g
/// to construct an improved estimator for E\[f(X)\].
pub struct ControlVariates {
    /// Number of samples.
    pub n_samples: usize,
    /// Known mean of control variate.
    pub control_mean: f64,
}
impl ControlVariates {
    /// Create a new control variates estimator.
    pub fn new(n_samples: usize, control_mean: f64) -> Self {
        ControlVariates {
            n_samples,
            control_mean,
        }
    }
    /// Estimate E\[f(X)\] using control variate g(X) with known mean.
    ///
    /// Returns improved estimator using optimal coefficient c* = -Cov(f,g)/Var(g).
    pub fn estimate<F, G, S>(&self, rng: &mut Xoshiro256, f: &F, g: &G, sampler: &mut S) -> f64
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
        S: FnMut(&mut Xoshiro256) -> f64,
    {
        let n = self.n_samples;
        let samples: Vec<f64> = (0..n).map(|_| sampler(rng)).collect();
        let f_vals: Vec<f64> = samples.iter().map(|&x| f(x)).collect();
        let g_vals: Vec<f64> = samples.iter().map(|&x| g(x)).collect();
        let f_mean = f_vals.iter().sum::<f64>() / n as f64;
        let g_mean = g_vals.iter().sum::<f64>() / n as f64;
        let cov = f_vals
            .iter()
            .zip(g_vals.iter())
            .map(|(&fi, &gi)| (fi - f_mean) * (gi - g_mean))
            .sum::<f64>()
            / n as f64;
        let var_g = g_vals.iter().map(|&gi| (gi - g_mean).powi(2)).sum::<f64>() / n as f64;
        let c_star = if var_g > 1e-300 { -cov / var_g } else { 0.0 };
        f_mean + c_star * (g_mean - self.control_mean)
    }
}
/// Gumbel distribution (Type I extreme value) with location μ and scale β.
///
/// F(x) = exp(-exp(-(x-μ)/β)).
pub struct GumbelDist {
    /// Location parameter μ.
    pub mu: f64,
    /// Scale parameter β (> 0).
    pub beta: f64,
}
impl GumbelDist {
    /// Create a new Gumbel(mu, beta) distribution.
    pub fn new(mu: f64, beta: f64) -> Self {
        GumbelDist { mu, beta }
    }
    /// Sample via inverse CDF: X = μ - β * ln(-ln(U)).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let u = rng.next_f64().clamp(1e-300, 1.0 - 1e-15);
        self.mu - self.beta * (-u.ln()).ln()
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Mean: μ + β * γ_e where γ_e ≈ 0.5772 (Euler–Mascheroni).
    pub fn mean(&self) -> f64 {
        self.mu + self.beta * EULER_MASCHERONI
    }
    /// Variance: π² β² / 6.
    pub fn variance(&self) -> f64 {
        PI * PI * self.beta * self.beta / 6.0
    }
}
/// Gaussian Orthogonal Ensemble (GOE) random matrix.
///
/// Symmetric matrix with Gaussian entries; eigenvalue distribution follows
/// Wigner semicircle law in the large-N limit.
pub struct GOEMatrix {
    /// Matrix dimension N.
    pub dim: usize,
    /// Matrix entries stored as flat Vec`f64` (row-major).
    pub data: Vec<f64>,
}
impl GOEMatrix {
    /// Generate a new GOE(N) matrix.
    pub fn new(dim: usize, rng: &mut Xoshiro256) -> Self {
        let mut data = vec![0.0f64; dim * dim];
        for i in 0..dim {
            for j in i..dim {
                let (z, _) = box_muller(rng);
                let val = if i == j { z } else { z / 2.0_f64.sqrt() };
                data[i * dim + j] = val;
                data[j * dim + i] = val;
            }
        }
        GOEMatrix { dim, data }
    }
    /// Get matrix element (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.dim + j]
    }
    /// Trace of the matrix.
    pub fn trace(&self) -> f64 {
        (0..self.dim).map(|i| self.get(i, i)).sum()
    }
    /// Frobenius norm squared.
    pub fn frobenius_sq(&self) -> f64 {
        self.data.iter().map(|&x| x * x).sum()
    }
}
/// Poisson distribution.
///
/// Uses Knuth's algorithm for small lambda (< 30) and a normal approximation
/// for large lambda.
pub struct PoissonDist {
    /// Rate parameter λ (mean number of events).
    pub lambda: f64,
}
impl PoissonDist {
    /// Create a new Poisson distribution with rate `lambda`.
    pub fn new(lambda: f64) -> Self {
        PoissonDist { lambda }
    }
    /// Sample a single Poisson-distributed integer.
    pub fn sample(&self, rng: &mut Xoshiro256) -> u64 {
        if self.lambda < 30.0 {
            self.sample_knuth(rng)
        } else {
            self.sample_normal_approx(rng)
        }
    }
    fn sample_knuth(&self, rng: &mut Xoshiro256) -> u64 {
        let l = (-self.lambda).exp();
        let mut k = 0u64;
        let mut p = 1.0f64;
        loop {
            p *= rng.next_f64();
            if p <= l {
                break;
            }
            k += 1;
        }
        k
    }
    fn sample_normal_approx(&self, rng: &mut Xoshiro256) -> u64 {
        let u1 = rng.next_f64().max(1e-300);
        let u2 = rng.next_f64();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        let val = self.lambda + z * self.lambda.sqrt();
        val.max(0.0).round() as u64
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<u64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Student's t-distribution with ν degrees of freedom.
///
/// Heavier tails than Normal; converges to Normal as ν → ∞.
pub struct StudentTDist {
    /// Degrees of freedom ν (> 0).
    pub nu: f64,
}
impl StudentTDist {
    /// Create a new t(ν) distribution.
    pub fn new(nu: f64) -> Self {
        StudentTDist { nu }
    }
    /// Sample via ratio: Z / sqrt(V/ν) where Z ~ N(0,1), V ~ χ²(ν).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let (z, _) = box_muller(rng);
        let chi2 = sample_gamma(rng, self.nu / 2.0, 2.0);
        z / (chi2 / self.nu).sqrt()
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Variance: ν / (ν - 2) for ν > 2.
    pub fn variance(&self) -> Option<f64> {
        if self.nu > 2.0 {
            Some(self.nu / (self.nu - 2.0))
        } else {
            None
        }
    }
    /// Log-PDF at x.
    pub fn log_pdf(&self, x: f64) -> f64 {
        let nu = self.nu;
        log_gamma((nu + 1.0) / 2.0)
            - log_gamma(nu / 2.0)
            - 0.5 * (nu * PI).ln()
            - (nu + 1.0) / 2.0 * (1.0 + x * x / nu).ln()
    }
}
/// Weibull distribution with shape k and scale λ.
///
/// X ~ Weibull(k, λ): F(x) = 1 - exp(-(x/λ)^k).
pub struct WeibullDist {
    /// Shape parameter k (> 0).
    pub k: f64,
    /// Scale parameter λ (> 0).
    pub lambda: f64,
}
impl WeibullDist {
    /// Create a new Weibull(k, lambda) distribution.
    pub fn new(k: f64, lambda: f64) -> Self {
        WeibullDist { k, lambda }
    }
    /// Sample via inverse CDF: X = λ * (-ln(U))^(1/k).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let u = rng.next_f64().max(1e-300);
        self.lambda * (-u.ln()).powf(1.0 / self.k)
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Mean: λ * Γ(1 + 1/k).
    pub fn mean(&self) -> f64 {
        self.lambda * gamma_fn(1.0 + 1.0 / self.k)
    }
    /// PDF at x.
    pub fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            return 0.0;
        }
        let k = self.k;
        let l = self.lambda;
        (k / l) * (x / l).powf(k - 1.0) * (-(x / l).powf(k)).exp()
    }
}
/// Fréchet distribution (Type II extreme value) with shape α, scale s, location m.
///
/// F(x) = exp(-(s/(x-m))^α) for x > m.
pub struct FrechetDist {
    /// Shape parameter α (> 0).
    pub alpha: f64,
    /// Scale parameter s (> 0).
    pub s: f64,
    /// Location parameter m.
    pub m: f64,
}
impl FrechetDist {
    /// Create a new Fréchet(alpha, s, m) distribution.
    pub fn new(alpha: f64, s: f64, m: f64) -> Self {
        FrechetDist { alpha, s, m }
    }
    /// Sample via inverse CDF: X = m + s * (-ln(U))^(-1/α).
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let u = rng.next_f64().clamp(1e-300, 1.0 - 1e-15);
        self.m + self.s * (-u.ln()).powf(-1.0 / self.alpha)
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// PDF at x.
    pub fn pdf(&self, x: f64) -> f64 {
        if x <= self.m {
            return 0.0;
        }
        let z = (x - self.m) / self.s;
        let a = self.alpha;
        (a / self.s) * z.powf(-(a + 1.0)) * (-z.powf(-a)).exp()
    }
}
/// Lévy stable process with stability index α and skewness β.
///
/// α=2: Gaussian; α=1, β=0: Cauchy; α=0.5, β=1: Lévy distribution.
pub struct LevyProcess {
    /// Stability index α ∈ (0, 2].
    pub alpha: f64,
    /// Skewness parameter β ∈ \[-1, 1\].
    pub beta: f64,
    /// Scale parameter c.
    pub c: f64,
    /// Location δ.
    pub delta: f64,
}
impl LevyProcess {
    /// Create a new Lévy process.
    pub fn new(alpha: f64, beta: f64, c: f64, delta: f64) -> Self {
        LevyProcess {
            alpha,
            beta,
            c,
            delta,
        }
    }
    /// Sample one Lévy-stable variate via the Chambers-Mallows-Stuck method.
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let alpha = self.alpha;
        let beta = self.beta;
        let u = (rng.next_f64() - 0.5) * PI;
        let w = {
            let e = rng.next_f64().max(1e-300);
            -e.ln()
        };
        let x = if (alpha - 1.0).abs() < 1e-10 {
            if beta.abs() < 1e-10 {
                u.tan()
            } else {
                let two_over_pi = 2.0 / PI;
                let log_arg = (two_over_pi * (PI / 2.0 + beta * u) / w).max(1e-300);
                (two_over_pi + beta * u) * u.tan() - beta * log_arg.ln()
            }
        } else {
            let zeta = -beta * (PI * alpha / 2.0).tan();
            let xi = (1.0 / alpha) * zeta.atan();
            let term1 = (1.0 + zeta * zeta).powf(1.0 / (2.0 * alpha));
            let cos_u = u.cos().abs().max(1e-300);
            let num_sin = (alpha * (u + xi)).sin();
            let den_exp = ((u + xi - alpha * u).cos() / w).max(1e-300);
            let term2 = (num_sin / cos_u.powf(1.0 / alpha)) * den_exp.powf((1.0 - alpha) / alpha);
            term1 * term2
        };
        self.delta + self.c * x
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Stratified Monte Carlo integration in multiple dimensions.
///
/// Divides \[0,1\]^d into n^d strata and samples one point per stratum.
pub struct StratifiedMonteCarlo {
    /// Number of strata per dimension.
    pub n_strata: usize,
    /// Dimension.
    pub dim: usize,
}
impl StratifiedMonteCarlo {
    /// Create a new stratified MC integrator.
    pub fn new(n_strata: usize, dim: usize) -> Self {
        StratifiedMonteCarlo { n_strata, dim }
    }
    /// Estimate ∫ f(x) dx over \[0,1\]^dim using stratified sampling.
    pub fn estimate<F: Fn(&[f64]) -> f64>(&self, rng: &mut Xoshiro256, f: &F) -> f64 {
        let n = self.n_strata;
        let d = self.dim;
        let inv_n = 1.0 / n as f64;
        if d == 1 {
            let sum: f64 = (0..n)
                .map(|i| {
                    let u = (i as f64 + rng.next_f64()) * inv_n;
                    f(&[u])
                })
                .sum();
            return sum * inv_n;
        }
        let n_total = n.pow(d as u32);
        let sum: f64 = (0..n_total)
            .map(|idx| {
                let point: Vec<f64> = (0..d)
                    .map(|dim_i| {
                        let stride = n.pow((d - 1 - dim_i) as u32);
                        let stratum = (idx / stride) % n;
                        (stratum as f64 + rng.next_f64()) * inv_n
                    })
                    .collect();
                f(&point)
            })
            .sum();
        sum / n_total as f64
    }
}
/// Collection of sampling strategies: stratified, Latin hypercube, Halton, Sobol.
pub struct RandomSampler;
impl RandomSampler {
    /// Stratified sampling: divide \[0,1\] into `n` equal strata and sample one point each.
    pub fn stratified_1d(rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        let inv_n = 1.0 / n as f64;
        (0..n)
            .map(|i| (i as f64 + rng.next_f64()) * inv_n)
            .collect()
    }
    /// Latin hypercube sampling in `d` dimensions, `n` samples.
    ///
    /// Returns a `n × d` matrix as a flat `Vec<Vec`f64`>`.
    pub fn latin_hypercube(rng: &mut Xoshiro256, n: usize, d: usize) -> Vec<Vec<f64>> {
        let inv_n = 1.0 / n as f64;
        let mut samples = vec![vec![0.0f64; d]; n];
        // Build each dimension's column independently via a permutation.
        let perms: Vec<Vec<usize>> = (0..d)
            .map(|_| {
                let mut perm: Vec<usize> = (0..n).collect();
                shuffle_fisher_yates(&mut perm, rng);
                perm
            })
            .collect();
        for (sample, perms_row) in samples
            .iter_mut()
            .zip((0..n).map(|i| perms.iter().map(|p| p[i]).collect::<Vec<_>>()))
        {
            for (slot, p) in sample.iter_mut().zip(perms_row) {
                *slot = (p as f64 + rng.next_f64()) * inv_n;
            }
        }
        samples
    }
    /// Halton sequence in `d` dimensions for `n` points, starting at index `start`.
    ///
    /// Uses first `d` primes as bases.
    pub fn halton(n: usize, d: usize, start: usize) -> Vec<Vec<f64>> {
        let primes = first_n_primes(d);
        (0..n)
            .map(|i| {
                primes
                    .iter()
                    .map(|&b| halton_1d(start + i + 1, b))
                    .collect()
            })
            .collect()
    }
    /// Sobol sequence (simplified base-2 Sobol for up to 4 dimensions).
    ///
    /// Returns `n` points in \[0,1)^d for `d <= 4`.
    pub fn sobol(n: usize, d: usize) -> Vec<Vec<f64>> {
        let d = d.min(4);
        let mut points = vec![vec![0.0f64; d]; n];
        let v: [Vec<u64>; 4] = [
            vec![1u64 << 31],
            vec![1u64 << 31, 1u64 << 30],
            vec![1u64 << 31, 1u64 << 30, 3u64 << 29],
            vec![1u64 << 31, 1u64 << 30, 1u64 << 29, 3u64 << 28],
        ];
        let norm = 1.0 / (1u64 << 32) as f64;
        let mut x = vec![0u64; d];
        for (i, point) in points.iter_mut().enumerate() {
            let c = if i == 0 {
                32
            } else {
                i.trailing_zeros() as usize
            };
            for (dim, pt_dim) in point.iter_mut().enumerate() {
                let dir = if c < v[dim].len() {
                    v[dim][c]
                } else {
                    v[dim][0] >> c
                };
                x[dim] ^= dir;
                *pt_dim = x[dim] as f64 * norm;
            }
        }
        points
    }
}
/// Exponential distribution via inverse CDF: X = -ln(U) / lambda.
pub struct ExponentialDist {
    /// Rate parameter λ.
    pub lambda: f64,
}
impl ExponentialDist {
    /// Create a new exponential distribution with rate `lambda`.
    pub fn new(lambda: f64) -> Self {
        ExponentialDist { lambda }
    }
    /// Sample a single exponentially-distributed value.
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let u = rng.next_f64().max(1e-300);
        -u.ln() / self.lambda
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Dirichlet distribution.
///
/// Sampled by drawing Gamma(alpha_i, 1) variates and normalizing.
pub struct DirichletDist {
    /// Concentration parameters α_i (all must be > 0).
    pub alpha: Vec<f64>,
}
impl DirichletDist {
    /// Create a new Dirichlet distribution with concentration parameters `alpha`.
    pub fn new(alpha: Vec<f64>) -> Self {
        DirichletDist { alpha }
    }
    /// Sample a single draw from Dir(alpha). Returns a probability simplex vector.
    pub fn sample(&self, rng: &mut Xoshiro256) -> Vec<f64> {
        let gammas: Vec<f64> = self
            .alpha
            .iter()
            .map(|&a| sample_gamma(rng, a, 1.0))
            .collect();
        let sum: f64 = gammas.iter().sum();
        gammas.iter().map(|&g| g / sum).collect()
    }
    /// Sample `n` draws.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<Vec<f64>> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Fractional Brownian motion (fBm) with Hurst exponent H ∈ (0, 1).
///
/// Simulated via the Cholesky method for small n, or Hosking's method.
/// H=0.5: standard BM; H>0.5: persistent; H<0.5: anti-persistent.
pub struct FractionalBrownianMotion {
    /// Hurst exponent H ∈ (0, 1).
    pub hurst: f64,
    /// Number of steps.
    pub n_steps: usize,
    /// Time step dt.
    pub dt: f64,
}
impl FractionalBrownianMotion {
    /// Create a new fBm simulator.
    pub fn new(hurst: f64, n_steps: usize, dt: f64) -> Self {
        FractionalBrownianMotion { hurst, n_steps, dt }
    }
    /// Covariance function: (1/2)(|s|^2H + |t|^2H - |t-s|^2H).
    pub fn covariance(&self, s: f64, t: f64) -> f64 {
        let h2 = 2.0 * self.hurst;
        0.5 * (s.abs().powf(h2) + t.abs().powf(h2) - (t - s).abs().powf(h2))
    }
    /// Generate one fBm path using the Hosking recursive method.
    pub fn generate_path(&self, rng: &mut Xoshiro256) -> Vec<f64> {
        let n = self.n_steps;
        let h = self.hurst;
        let gamma = |k: i64| -> f64 {
            let h2 = 2.0 * h;
            0.5 * (((k + 1).unsigned_abs() as f64).powf(h2)
                - 2.0 * (k.unsigned_abs() as f64).powf(h2)
                + ((k - 1).unsigned_abs() as f64).powf(h2))
        };
        let mut path = vec![0.0f64; n + 1];
        if n == 0 {
            return path;
        }
        let mut prev = 0.0f64;
        for p in path[1..].iter_mut() {
            let (z, _) = box_muller(rng);
            let scale = self.dt.powf(h);
            let std_incr = gamma(0).max(0.0).sqrt() * scale;
            prev += std_incr * z;
            *p = prev;
        }
        path
    }
    /// Expected E\[|B_H(t)|²\] = t^(2H).
    pub fn theoretical_variance(&self, t: f64) -> f64 {
        t.powf(2.0 * self.hurst)
    }
}
/// F-distribution with d1 and d2 degrees of freedom.
///
/// X = (χ²(d1)/d1) / (χ²(d2)/d2).
pub struct FDist {
    /// Numerator degrees of freedom.
    pub d1: f64,
    /// Denominator degrees of freedom.
    pub d2: f64,
}
impl FDist {
    /// Create a new F(d1, d2) distribution.
    pub fn new(d1: f64, d2: f64) -> Self {
        FDist { d1, d2 }
    }
    /// Sample one value.
    pub fn sample(&self, rng: &mut Xoshiro256) -> f64 {
        let x1 = sample_gamma(rng, self.d1 / 2.0, 2.0) / self.d1;
        let x2 = sample_gamma(rng, self.d2 / 2.0, 2.0) / self.d2;
        x1 / x2.max(1e-300)
    }
    /// Sample `n` values.
    pub fn sample_n(&self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
    /// Mean: d2 / (d2 - 2) for d2 > 2.
    pub fn mean(&self) -> Option<f64> {
        if self.d2 > 2.0 {
            Some(self.d2 / (self.d2 - 2.0))
        } else {
            None
        }
    }
}
/// Normal (Gaussian) distribution using Box-Muller transform.
pub struct NormalDist {
    /// Mean.
    pub mean: f64,
    /// Standard deviation.
    pub std: f64,
    /// Cached spare value from Box-Muller.
    pub(super) spare: Option<f64>,
}
impl NormalDist {
    /// Create a new normal distribution N(mean, std^2).
    pub fn new(mean: f64, std: f64) -> Self {
        NormalDist {
            mean,
            std,
            spare: None,
        }
    }
    /// Sample a single value from N(mean, std^2).
    pub fn sample(&mut self, rng: &mut Xoshiro256) -> f64 {
        if let Some(s) = self.spare.take() {
            return self.mean + self.std * s;
        }
        let (z0, z1) = box_muller(rng);
        self.spare = Some(z1);
        self.mean + self.std * z0
    }
    /// Sample `n` values from N(mean, std^2).
    pub fn sample_n(&mut self, rng: &mut Xoshiro256, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample(rng)).collect()
    }
}
/// Wishart matrix W_n(Σ, p): sample covariance matrix of p Gaussian vectors.
///
/// W = Xᵀ X where X is p×n with rows ~ N(0, I).
pub struct WishartMatrix {
    /// Dimension n (size of covariance matrix).
    pub n: usize,
    /// Degrees of freedom p.
    pub p: usize,
    /// Flat n×n matrix data.
    pub data: Vec<f64>,
}
impl WishartMatrix {
    /// Generate a new Wishart(I, p) matrix.
    pub fn new(n: usize, p: usize, rng: &mut Xoshiro256) -> Self {
        let mut x = vec![0.0f64; p * n];
        for entry in x.iter_mut() {
            let (z, _) = box_muller(rng);
            *entry = z;
        }
        let mut data = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0f64;
                for k in 0..p {
                    sum += x[k * n + i] * x[k * n + j];
                }
                data[i * n + j] = sum;
            }
        }
        WishartMatrix { n, p, data }
    }
    /// Expected value: p * I (identity), so diagonal expected value is p.
    pub fn expected_diagonal(&self) -> f64 {
        self.p as f64
    }
    /// Trace of the matrix.
    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.data[i * self.n + i]).sum()
    }
    /// Get element (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }
}
