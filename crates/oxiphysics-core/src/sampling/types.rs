//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Type alias for a conditional mean function used in Gibbs sampling.
pub type ConditionalMeanFn = Box<dyn Fn(usize, &[f64]) -> f64>;

/// Hammersley low-discrepancy sequence in `d` dimensions.
///
/// The first dimension is `i/n` (uniform grid); subsequent dimensions use
/// the van der Corput sequence with the first `d-1` prime bases.
///
/// Returns `n` points in `[0,1)^d`.
pub struct HammersleySequence;
impl HammersleySequence {
    /// Generate `n` Hammersley points in `d` dimensions.
    pub fn sample(n: usize, d: usize) -> Vec<Vec<f64>> {
        const PRIMES: [u32; 9] = [2, 3, 5, 7, 11, 13, 17, 19, 23];
        let n_bases = (d - 1).min(PRIMES.len());
        (0..n)
            .map(|i| {
                let mut pt = Vec::with_capacity(d);
                pt.push(i as f64 / n as f64);
                for &p in PRIMES[..n_bases].iter() {
                    pt.push(HaltonSequence::van_der_corput(i as u32 + 1, p));
                }
                while pt.len() < d {
                    pt.push(0.0);
                }
                pt
            })
            .collect()
    }
}
/// Metropolis-Hastings Markov Chain Monte Carlo sampler (1-D).
///
/// Generates samples from a target distribution specified by its log-probability.
/// Uses a Gaussian proposal distribution with standard deviation `step_size`.
pub struct MetropolisHastings {
    /// Current state of the chain.
    pub current: f64,
    /// Standard deviation of the Gaussian proposal.
    pub step_size: f64,
    /// Internal RNG.
    pub(super) rng: Lcg,
}
impl MetropolisHastings {
    /// Create a new Metropolis-Hastings sampler.
    ///
    /// # Arguments
    /// * `initial`   – initial state of the chain.
    /// * `step_size` – standard deviation of the proposal distribution.
    /// * `seed`      – RNG seed.
    pub fn new(initial: f64, step_size: f64, seed: u64) -> Self {
        Self {
            current: initial,
            step_size,
            rng: Lcg::new(seed),
        }
    }
    /// Draw `n` samples from the target log-PDF `log_target`.
    ///
    /// The chain runs for `n` steps and returns all states (including burn-in).
    /// Use `samples[burn_in..]` for post-burn-in inference.
    pub fn sample(&mut self, n: usize, log_target: impl Fn(f64) -> f64) -> Vec<f64> {
        let mut samples = Vec::with_capacity(n);
        let mut log_current = log_target(self.current);
        for _ in 0..n {
            let proposal = self.current + self.rng.next_normal() * self.step_size;
            let log_proposal = log_target(proposal);
            let log_alpha = log_proposal - log_current;
            let u = self.rng.next_f64();
            if u.ln() < log_alpha {
                self.current = proposal;
                log_current = log_proposal;
            }
            samples.push(self.current);
        }
        samples
    }
    /// Acceptance rate estimate over `n` trial steps.
    pub fn estimate_acceptance_rate(&mut self, n: usize, log_target: impl Fn(f64) -> f64) -> f64 {
        let mut accepted = 0usize;
        let mut log_current = log_target(self.current);
        for _ in 0..n {
            let proposal = self.current + self.rng.next_normal() * self.step_size;
            let log_proposal = log_target(proposal);
            let log_alpha = log_proposal - log_current;
            if self.rng.next_f64().ln() < log_alpha {
                self.current = proposal;
                log_current = log_proposal;
                accepted += 1;
            }
        }
        accepted as f64 / n as f64
    }
}
/// Linear Congruential Generator using Knuth's constants.
///
/// `state_{n+1} = a * state_n + c  (mod 2^64)`
///
/// * `a` = 6 364 136 223 846 793 005
/// * `c` = 1 442 695 040 888 963 407
pub struct Lcg {
    pub(super) state: u64,
}
impl Lcg {
    const A: u64 = 6_364_136_223_846_793_005;
    const C: u64 = 1_442_695_040_888_963_407;
    /// Create a new LCG with the given seed.
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    /// Advance the generator and return the next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(Self::A).wrapping_add(Self::C);
        self.state
    }
    /// Return a uniform float in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        let bits = self.next_u64() >> 11;
        bits as f64 / (1u64 << 53) as f64
    }
    /// Return a uniform float in `[lo, hi)`.
    pub fn next_f64_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.next_f64() * (hi - lo)
    }
    /// Return a standard-normal variate using the Box–Muller transform.
    ///
    /// Each call consumes two uniform samples and discards one of the two
    /// normals produced. Use [`next_normal_pair`](Self::next_normal_pair) when
    /// both are needed.
    pub fn next_normal(&mut self) -> f64 {
        self.next_normal_pair().0
    }
    /// Return two independent standard-normal variates via Box–Muller.
    pub fn next_normal_pair(&mut self) -> (f64, f64) {
        use std::f64::consts::PI;
        let u1 = 1.0 - self.next_f64();
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * PI * u2;
        (r * theta.cos(), r * theta.sin())
    }
}
/// Simplified 1-D Sobol quasi-random sequence.
///
/// Direction numbers are the powers of two shifted into the most-significant
/// bits of a 32-bit word: `[1<<31, 1<<30, 1<<29, …]`, which corresponds to
/// the standard `m_i = i` initialisation for dimension 1.
pub struct Sobol {
    /// Number of active dimensions (always 1 for this simplified version).
    pub(super) _dim: usize,
    /// How many points have been generated so far.
    pub(super) count: usize,
    /// Direction numbers, one `Vec`u32` per dimension.
    pub(super) direction_numbers: Vec<Vec<u32>>,
}
impl Sobol {
    const BITS: usize = 32;
    /// Create a 1-D Sobol generator.
    pub fn new_1d() -> Self {
        let v: Vec<u32> = (0..Self::BITS)
            .map(|i| 1u32 << (Self::BITS - 1 - i))
            .collect();
        Self {
            _dim: 1,
            count: 0,
            direction_numbers: vec![v],
        }
    }
    /// Return the next 1-D Sobol point in `\[0, 1)`.
    pub fn next_1d(&mut self) -> f64 {
        let v = &self.direction_numbers[0];
        let c = self.count;
        self.count += 1;
        let gray = c ^ (c >> 1);
        let x = (0..Self::BITS)
            .filter(|&i| (gray >> i) & 1 == 1)
            .fold(0u32, |acc, i| acc ^ v[i]);
        x as f64 / (1u64 << Self::BITS) as f64
    }
}
/// Multi-dimensional Sobol quasi-random sequence (dims 1-3).
///
/// Uses standard direction numbers for dimensions 1, 2, and 3.
pub struct SobolSequence {
    /// Number of dimensions (1..=3).
    pub n_dims: usize,
}
impl SobolSequence {
    const BITS: usize = 32;
    /// Create a new Sobol sequence generator for `n_dims` dimensions (max 3).
    pub fn new(n_dims: usize) -> Self {
        Self {
            n_dims: n_dims.clamp(1, 3),
        }
    }
    /// Direction numbers for dimension `dim` (0-indexed).
    pub(crate) fn direction_numbers(dim: usize) -> Vec<u32> {
        match dim {
            0 => (0..Self::BITS)
                .map(|i| 1u32 << (Self::BITS - 1 - i))
                .collect(),
            1 => {
                let mut v = vec![0u32; Self::BITS];
                v[0] = 1u32 << (Self::BITS - 1);
                for i in 1..Self::BITS {
                    v[i] = v[i - 1] ^ (v[i - 1] >> 1);
                }
                v
            }
            _ => {
                let mut v = vec![0u32; Self::BITS];
                v[0] = 1u32 << (Self::BITS - 1);
                v[1] = 1u32 << (Self::BITS - 2);
                for i in 2..Self::BITS {
                    v[i] = v[i - 2] ^ (v[i - 1] >> 1) ^ v[i - 1];
                }
                v
            }
        }
    }
    /// Generate `n` Sobol points.  Each point is a `Vec`f64` of length `n_dims`
    /// with values in `[0, 1)`.
    pub fn sample(&self, n: usize) -> Vec<Vec<f64>> {
        let dir: Vec<Vec<u32>> = (0..self.n_dims).map(Self::direction_numbers).collect();
        let mut result = Vec::with_capacity(n);
        for idx in 0..n {
            let gray = idx ^ (idx >> 1);
            let mut point = Vec::with_capacity(self.n_dims);
            for dir_d in dir.iter() {
                let x = (0..Self::BITS)
                    .filter(|&i| (gray >> i) & 1 == 1)
                    .fold(0u32, |acc, i| acc ^ dir_d[i]);
                point.push(x as f64 / (1u64 << Self::BITS) as f64);
            }
            result.push(point);
        }
        result
    }
}
/// Latin Hypercube sampler.
///
/// Divides each dimension into `n_samples` equally-probable intervals and
/// places exactly one sample in each interval per dimension.
pub struct LatinHypercube {
    /// Number of samples to generate.
    pub n_samples: usize,
    /// Number of dimensions.
    pub n_dims: usize,
}
impl LatinHypercube {
    /// Create a new Latin Hypercube sampler.
    pub fn new(n_samples: usize, n_dims: usize) -> Self {
        Self { n_samples, n_dims }
    }
    /// Generate Latin Hypercube samples.  Returns `n_samples` points, each of
    /// length `n_dims`, with values in `[0, 1)`.
    pub fn sample(&self, rng: &mut Lcg) -> Vec<Vec<f64>> {
        let n = self.n_samples;
        let inv_n = 1.0 / n as f64;
        let mut result: Vec<Vec<f64>> = (0..n).map(|_| Vec::with_capacity(self.n_dims)).collect();
        for _d in 0..self.n_dims {
            let mut perm: Vec<usize> = (0..n).collect();
            for i in (1..n).rev() {
                let j = (rng.next_u64() as usize) % (i + 1);
                perm.swap(i, j);
            }
            for i in 0..n {
                let val = (perm[i] as f64 + rng.next_f64()) * inv_n;
                result[i].push(val);
            }
        }
        result
    }
}
/// Stratified 2-D sampler.
///
/// Divides the unit square `[0,1)^2` into `n_strata_per_dim^2` equal cells and
/// places one jittered sample in each cell.
pub struct StratifiedSampler {
    /// Number of strata per dimension.
    pub n_strata_per_dim: usize,
}
impl StratifiedSampler {
    /// Create a new stratified sampler.
    pub fn new(n_strata_per_dim: usize) -> Self {
        Self { n_strata_per_dim }
    }
    /// Generate stratified 2-D samples.  Returns `n_strata_per_dim^2` points.
    pub fn sample_2d(&self, rng: &mut Lcg) -> Vec<[f64; 2]> {
        let n = self.n_strata_per_dim;
        let inv_n = 1.0 / n as f64;
        let mut samples = Vec::with_capacity(n * n);
        for iy in 0..n {
            for ix in 0..n {
                let x = (ix as f64 + rng.next_f64()) * inv_n;
                let y = (iy as f64 + rng.next_f64()) * inv_n;
                samples.push([x, y]);
            }
        }
        samples
    }
}
/// Rejection-based importance sampler.
///
/// Draws samples from a target PDF using a uniform proposal on `[0, 1)` and
/// rejection sampling.
pub struct ImportanceSampler {
    /// Target PDF (un-normalised is fine as long as `pdf(x) ≤ envelope`).
    pub pdf_fn: fn(f64) -> f64,
    /// Upper bound on `pdf_fn(x)` for `x ∈ [0, 1)`.
    pub envelope: f64,
}
impl ImportanceSampler {
    /// Create a new importance sampler.
    pub fn new(pdf_fn: fn(f64) -> f64, envelope: f64) -> Self {
        Self { pdf_fn, envelope }
    }
    /// Draw `n` samples from the target PDF via rejection sampling.
    pub fn sample(&self, n: usize, rng: &mut Lcg) -> Vec<f64> {
        let mut samples = Vec::with_capacity(n);
        while samples.len() < n {
            let x = rng.next_f64();
            let u = rng.next_f64() * self.envelope;
            if u < (self.pdf_fn)(x) {
                samples.push(x);
            }
        }
        samples
    }
}
/// Simple Gibbs sampler for multi-dimensional distributions.
///
/// Updates each dimension in turn by sampling from the conditional distribution.
pub struct GibbsSampler {
    /// Number of dimensions.
    pub n_dims: usize,
    /// Current state.
    pub state: Vec<f64>,
    /// Internal RNG.
    pub(super) rng: Lcg,
}
impl GibbsSampler {
    /// Create a new Gibbs sampler initialized at zero.
    pub fn new(n_dims: usize, seed: u64) -> Self {
        Self {
            n_dims,
            state: vec![0.0; n_dims],
            rng: Lcg::new(seed),
        }
    }
    /// Draw `n` samples.
    ///
    /// `conditional_means[i](i, state)` returns the mean of the `i`-th
    /// conditional distribution given the current state (Gaussian with std
    /// `conditional_stds[i]`).
    pub fn sample(
        &mut self,
        n: usize,
        conditional_means: &[ConditionalMeanFn],
        conditional_stds: &[f64],
    ) -> Vec<Vec<f64>> {
        let mut samples = Vec::with_capacity(n);
        for _ in 0..n {
            for d in 0..self.n_dims {
                let mean = conditional_means[d](d, &self.state);
                let std = conditional_stds[d];
                self.state[d] = mean + self.rng.next_normal() * std;
            }
            samples.push(self.state.clone());
        }
        samples
    }
}
/// Gaussian kernel density estimator (KDE).
///
/// Evaluates the estimated density at each point in `query` given training
/// `data`, using bandwidth `h`.  Uses the Gaussian kernel:
///
/// `k(x, xi, h) = 1/(h * √(2π)) * exp(-0.5 * ((x - xi)/h)²)`.
pub struct GaussianKde {
    /// Training data.
    pub data: Vec<f64>,
    /// Bandwidth (smoothing parameter).
    pub bandwidth: f64,
}
impl GaussianKde {
    /// Create a new Gaussian KDE, optionally with Silverman's rule of thumb
    /// if `bandwidth` is 0 or negative.
    pub fn new(data: Vec<f64>, bandwidth: f64) -> Self {
        let h = if bandwidth > 0.0 {
            bandwidth
        } else {
            Self::silverman_bandwidth(&data)
        };
        Self { data, bandwidth: h }
    }
    /// Silverman's rule of thumb bandwidth: `1.06 * σ * n^{-1/5}`.
    pub fn silverman_bandwidth(data: &[f64]) -> f64 {
        let n = data.len() as f64;
        if n < 2.0 {
            return 1.0;
        }
        let mean = data.iter().sum::<f64>() / n;
        let var = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
        let sigma = var.sqrt().max(1e-300);
        1.06 * sigma * n.powf(-0.2)
    }
    /// Evaluate the KDE at a single query point `x`.
    pub fn evaluate(&self, x: f64) -> f64 {
        let h = self.bandwidth;
        let n = self.data.len() as f64;
        let norm = 1.0 / (h * (2.0 * std::f64::consts::PI).sqrt());
        self.data
            .iter()
            .map(|&xi| {
                let z = (x - xi) / h;
                norm * (-0.5 * z * z).exp()
            })
            .sum::<f64>()
            / n
    }
    /// Evaluate the KDE at multiple query points.
    pub fn evaluate_batch(&self, xs: &[f64]) -> Vec<f64> {
        xs.iter().map(|&x| self.evaluate(x)).collect()
    }
}
/// Halton (van der Corput) quasi-random sequence.
pub struct HaltonSequence;
impl HaltonSequence {
    /// Generate `n` points of the van der Corput sequence in the given `base`.
    ///
    /// Each value lies in `[0, 1)`.
    pub fn sample(n: usize, base: u32) -> Vec<f64> {
        (0..n)
            .map(|i| Self::van_der_corput(i as u32 + 1, base))
            .collect()
    }
    /// Compute the van der Corput value for index `i` in the given `base`.
    pub fn van_der_corput(mut i: u32, base: u32) -> f64 {
        let mut result = 0.0_f64;
        let mut denom = 1.0_f64;
        while i > 0 {
            denom *= base as f64;
            result += (i % base) as f64 / denom;
            i /= base;
        }
        result
    }
}
