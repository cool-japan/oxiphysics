//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::collections::VecDeque;

/// Beta distribution with shape parameters `alpha` and `beta`.
pub struct BetaDistribution {
    /// First shape parameter α (must be > 0).
    pub alpha: f64,
    /// Second shape parameter β (must be > 0).
    pub beta: f64,
}
impl BetaDistribution {
    /// Creates a new [`BetaDistribution`] with shape parameters `alpha` and `beta`.
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }
    /// Probability density function evaluated at `x` ∈ (0, 1).
    ///
    /// The normalisation constant is computed via the log-gamma (Stirling)
    /// approximation for the beta function.
    pub fn pdf(&self, x: f64) -> f64 {
        if x <= 0.0 || x >= 1.0 {
            return 0.0;
        }
        let log_b = log_gamma_stirling(self.alpha) + log_gamma_stirling(self.beta)
            - log_gamma_stirling(self.alpha + self.beta);
        ((self.alpha - 1.0) * x.ln() + (self.beta - 1.0) * (1.0 - x).ln() - log_b).exp()
    }
    /// Mean of the distribution: α / (α + β).
    pub fn mean(&self) -> f64 {
        self.alpha / (self.alpha + self.beta)
    }
    /// Variance of the distribution.
    pub fn variance(&self) -> f64 {
        let s = self.alpha + self.beta;
        self.alpha * self.beta / (s * s * (s + 1.0))
    }
}
/// A simple linear-congruential pseudo-random number generator.
///
/// Uses the Knuth/MMIX parameters:
/// - multiplier = 6364136223846793005
/// - increment  = 1442695040888963407
///
/// Suitable for Monte Carlo work where a full-quality CSPRNG is not required
/// and where pulling in an external crate is undesirable.
pub struct StatRng {
    pub(super) state: u64,
}
impl StatRng {
    /// Creates a new [`StatRng`] seeded with `seed`.
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }
    /// Returns the next pseudo-random value in \[0, 1).
    pub fn next_f64(&mut self) -> f64 {
        const A: u64 = 6364136223846793005;
        const C: u64 = 1442695040888963407;
        self.state = self.state.wrapping_mul(A).wrapping_add(C);
        (self.state >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Returns the next sample from a standard normal distribution N(0, 1)
    /// using the Box–Muller transform.
    pub fn next_normal(&mut self) -> f64 {
        loop {
            let u1 = self.next_f64();
            let u2 = self.next_f64();
            if u1 > 0.0 {
                let mag = (-2.0 * u1.ln()).sqrt();
                return mag * (std::f64::consts::TAU * u2).cos();
            }
        }
    }
}
/// Poisson distribution with mean `lambda`.
pub struct PoissonDistribution {
    /// Mean (and variance) of the distribution λ (must be > 0).
    pub lambda: f64,
}
impl PoissonDistribution {
    /// Creates a new [`PoissonDistribution`] with mean `lambda`.
    pub fn new(lambda: f64) -> Self {
        Self { lambda }
    }
    /// Probability mass function: P(X = k) = λ^k e^{-λ} / k!.
    pub fn pmf(&self, k: u64) -> f64 {
        let log_p = k as f64 * self.lambda.ln() - self.lambda - log_factorial(k);
        log_p.exp()
    }
    /// Draws one sample using Knuth's algorithm.
    pub fn sample(&self, rng: &mut StatRng) -> u64 {
        let threshold = (-self.lambda).exp();
        let mut p = 1.0;
        let mut k = 0u64;
        loop {
            p *= rng.next_f64();
            if p <= threshold {
                return k;
            }
            k += 1;
        }
    }
}
/// Two-dimensional Gaussian kernel density estimate.
pub struct KernelDensityEstimate2D {
    /// Data points `(x, y)`.
    pub data: Vec<[f64; 2]>,
    /// Bandwidth along x.
    pub bw_x: f64,
    /// Bandwidth along y.
    pub bw_y: f64,
}
impl KernelDensityEstimate2D {
    /// Create a new 2-D KDE with given bandwidths.
    pub fn new(data: Vec<[f64; 2]>, bw_x: f64, bw_y: f64) -> Self {
        Self { data, bw_x, bw_y }
    }
    /// Evaluate the 2-D KDE at point `(x, y)`.
    pub fn evaluate(&self, x: f64, y: f64) -> f64 {
        if self.data.is_empty() || self.bw_x <= 0.0 || self.bw_y <= 0.0 {
            return 0.0;
        }
        let n = self.data.len() as f64;
        let norm = 1.0 / (std::f64::consts::TAU * self.bw_x * self.bw_y * n);
        self.data
            .iter()
            .map(|&[xi, yi]| {
                let zx = (x - xi) / self.bw_x;
                let zy = (y - yi) / self.bw_y;
                norm * (-0.5 * (zx * zx + zy * zy)).exp()
            })
            .sum()
    }
    /// Estimate optimal bandwidths using Scott's rule for each dimension.
    pub fn optimal_bandwidths(data: &[[f64; 2]]) -> (f64, f64) {
        let n = data.len();
        if n < 2 {
            return (1.0, 1.0);
        }
        let xs: Vec<f64> = data.iter().map(|p| p[0]).collect();
        let ys: Vec<f64> = data.iter().map(|p| p[1]).collect();
        let bw_x = KernelDensityEstimate::optimal_bandwidth(&xs);
        let bw_y = KernelDensityEstimate::optimal_bandwidth(&ys);
        (bw_x, bw_y)
    }
}
/// Maxwell–Boltzmann speed distribution parameterised by `a = sqrt(k_B T / m)`.
///
/// The PDF is f(v) = sqrt(2/π) · v² / a³ · exp(−v² / (2 a²)).
pub struct MaxwellBoltzmann {
    /// Scale parameter a = √(k_B T / m).
    pub a: f64,
}
impl MaxwellBoltzmann {
    /// Creates a new [`MaxwellBoltzmann`] distribution with scale parameter `a`.
    pub fn new(a: f64) -> Self {
        Self { a }
    }
    /// Probability density function evaluated at speed `v` ≥ 0.
    pub fn pdf(&self, v: f64) -> f64 {
        if v < 0.0 {
            return 0.0;
        }
        let a2 = self.a * self.a;
        (2.0 / std::f64::consts::PI).sqrt() * v * v / (a2 * self.a) * (-v * v / (2.0 * a2)).exp()
    }
    /// Mean speed: ⟨v⟩ = 2a √(2/π).
    pub fn mean_speed(&self) -> f64 {
        2.0 * self.a * (2.0 / std::f64::consts::PI).sqrt()
    }
    /// Root-mean-square speed: v_rms = a √3.
    pub fn rms_speed(&self) -> f64 {
        self.a * 3.0_f64.sqrt()
    }
    /// Draws one sample using the chi distribution (3 degrees of freedom).
    ///
    /// v = a · √(X₁² + X₂² + X₃²) where Xᵢ ~ N(0,1).
    pub fn sample(&self, rng: &mut StatRng) -> f64 {
        let x1 = rng.next_normal();
        let x2 = rng.next_normal();
        let x3 = rng.next_normal();
        self.a * (x1 * x1 + x2 * x2 + x3 * x3).sqrt()
    }
}
/// Normal (Gaussian) distribution N(`mean`, `std`²).
pub struct NormalDistribution {
    /// Mean of the distribution.
    pub mean: f64,
    /// Standard deviation of the distribution (must be > 0).
    pub std: f64,
}
impl NormalDistribution {
    /// Creates a new [`NormalDistribution`] with the given `mean` and `std`.
    pub fn new(mean: f64, std: f64) -> Self {
        Self { mean, std }
    }
    /// Probability density function evaluated at `x`.
    pub fn pdf(&self, x: f64) -> f64 {
        let z = (x - self.mean) / self.std;
        (-0.5 * z * z).exp() / (self.std * (std::f64::consts::TAU).sqrt())
    }
    /// Cumulative distribution function evaluated at `x`.
    ///
    /// Uses an accurate polynomial approximation of the error function.
    pub fn cdf(&self, x: f64) -> f64 {
        let z = (x - self.mean) / (self.std * std::f64::consts::SQRT_2);
        0.5 * (1.0 + erf_approx(z))
    }
    /// Draws one sample from this distribution using the supplied [`StatRng`].
    pub fn sample(&self, rng: &mut StatRng) -> f64 {
        self.mean + self.std * rng.next_normal()
    }
}
/// Gaussian kernel density estimate.
pub struct KernelDensityEstimate {
    /// The data points used to build the estimate.
    pub data: Vec<f64>,
    /// Bandwidth (smoothing parameter h).
    pub bandwidth: f64,
}
impl KernelDensityEstimate {
    /// Creates a new [`KernelDensityEstimate`] from `data` with the given `bandwidth`.
    pub fn new(data: Vec<f64>, bandwidth: f64) -> Self {
        Self { data, bandwidth }
    }
    /// Evaluates the KDE at point `x`.
    pub fn evaluate(&self, x: f64) -> f64 {
        if self.data.is_empty() || self.bandwidth <= 0.0 {
            return 0.0;
        }
        let n = self.data.len() as f64;
        let h = self.bandwidth;
        let norm = 1.0 / ((std::f64::consts::TAU).sqrt() * h * n);
        self.data
            .iter()
            .map(|&xi| {
                let z = (x - xi) / h;
                norm * (-0.5 * z * z).exp()
            })
            .sum()
    }
    /// Computes the optimal bandwidth using Scott's rule: h = 1.06 σ n^{-1/5}.
    pub fn optimal_bandwidth(data: &[f64]) -> f64 {
        let n = data.len();
        if n < 2 {
            return 1.0;
        }
        let sigma = std_dev(data);
        1.06 * sigma * (n as f64).powf(-0.2)
    }
}
/// A two-dimensional histogram.
pub struct Histogram2D {
    /// Number of bins along the x axis.
    pub nx: usize,
    /// Number of bins along the y axis.
    pub ny: usize,
    /// Range (min, max) along the x axis.
    pub x_range: (f64, f64),
    /// Range (min, max) along the y axis.
    pub y_range: (f64, f64),
    /// Flattened bin counts stored in row-major order (x-major).
    pub counts: Vec<usize>,
}
impl Histogram2D {
    /// Creates a new empty [`Histogram2D`].
    pub fn new(nx: usize, ny: usize, x_range: (f64, f64), y_range: (f64, f64)) -> Self {
        Self {
            nx,
            ny,
            x_range,
            y_range,
            counts: vec![0; nx * ny],
        }
    }
    /// Adds a data point (`x`, `y`) to the histogram.
    ///
    /// Points outside the declared ranges are silently ignored.
    pub fn add(&mut self, x: f64, y: f64) {
        let (xmin, xmax) = self.x_range;
        let (ymin, ymax) = self.y_range;
        if x < xmin || x > xmax || y < ymin || y > ymax {
            return;
        }
        let ix = (((x - xmin) / (xmax - xmin)) * self.nx as f64) as usize;
        let iy = (((y - ymin) / (ymax - ymin)) * self.ny as f64) as usize;
        let ix = ix.min(self.nx - 1);
        let iy = iy.min(self.ny - 1);
        self.counts[ix * self.ny + iy] += 1;
    }
    /// Returns the total number of data points added.
    pub fn total(&self) -> usize {
        self.counts.iter().sum()
    }
    /// Returns a density estimate: each bin count divided by (total × bin_area).
    ///
    /// Returns a vector of zeros when `total()` is zero.
    pub fn normalize(&self) -> Vec<f64> {
        let total = self.total();
        if total == 0 {
            return vec![0.0; self.nx * self.ny];
        }
        let (xmin, xmax) = self.x_range;
        let (ymin, ymax) = self.y_range;
        let dx = (xmax - xmin) / self.nx as f64;
        let dy = (ymax - ymin) / self.ny as f64;
        let bin_area = dx * dy;
        self.counts
            .iter()
            .map(|&c| c as f64 / (total as f64 * bin_area))
            .collect()
    }
}
/// Result of a PCA computation.
pub struct PcaResult {
    /// Mean vector of the data (length `d`).
    pub mean: Vec<f64>,
    /// Principal component directions (row-major, each row is a PC).
    /// Ordered from most to least variance.
    pub components: Vec<Vec<f64>>,
    /// Explained variance (eigenvalue) for each component.
    pub explained_variance: Vec<f64>,
}
/// Welford online algorithm for computing running mean and variance.
///
/// Processes one observation at a time with O(1) memory.  Numerically stable
/// for large data streams (West 1979 / Welford 1962).
pub struct WelfordOnline {
    /// Number of observations so far.
    pub n: usize,
    /// Current mean estimate.
    pub mean: f64,
    /// Current sum of squared deviations from the mean (M₂).
    pub(super) m2: f64,
}
impl WelfordOnline {
    /// Create a new online accumulator (no observations yet).
    pub fn new() -> Self {
        Self {
            n: 0,
            mean: 0.0,
            m2: 0.0,
        }
    }
    /// Incorporate a new observation `x`.
    pub fn update(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }
    /// Sample variance (divides by n − 1).  Returns `0.0` for n < 2.
    pub fn variance(&self) -> f64 {
        if self.n < 2 {
            return 0.0;
        }
        self.m2 / (self.n as f64 - 1.0)
    }
    /// Sample standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
    /// Population variance (divides by n).  Returns `0.0` for n == 0.
    pub fn population_variance(&self) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        self.m2 / self.n as f64
    }
}
/// Exponential distribution with rate parameter `lambda`.
pub struct ExponentialDistribution {
    /// Rate parameter λ (must be > 0).
    pub lambda: f64,
}
impl ExponentialDistribution {
    /// Creates a new [`ExponentialDistribution`] with rate `lambda`.
    pub fn new(lambda: f64) -> Self {
        Self { lambda }
    }
    /// Probability density function evaluated at `x`.
    pub fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            0.0
        } else {
            self.lambda * (-self.lambda * x).exp()
        }
    }
    /// Mean of the distribution: 1 / λ.
    pub fn mean(&self) -> f64 {
        1.0 / self.lambda
    }
    /// Draws one sample using the inverse-CDF (quantile) method.
    pub fn sample(&self, rng: &mut StatRng) -> f64 {
        let u = loop {
            let v = rng.next_f64();
            if v > 0.0 {
                break v;
            }
        };
        -u.ln() / self.lambda
    }
}
/// Welford accumulator for a running mean and variance over a sliding window.
///
/// Uses a `VecDeque` to drop old observations.  The window variance uses
/// a simple recalculation (O(window) cost) for correctness; for O(1) use the
/// global `WelfordOnline` on the whole stream instead.
pub struct SlidingWindowStats {
    pub(super) window: VecDeque<f64>,
    pub(super) capacity: usize,
}
impl SlidingWindowStats {
    /// Create a new sliding-window accumulator with the given `capacity`.
    pub fn new(capacity: usize) -> Self {
        Self {
            window: VecDeque::with_capacity(capacity.max(1)),
            capacity: capacity.max(1),
        }
    }
    /// Add an observation to the window, evicting the oldest if full.
    pub fn push(&mut self, x: f64) {
        if self.window.len() == self.capacity {
            self.window.pop_front();
        }
        self.window.push_back(x);
    }
    /// Current window mean.
    pub fn mean(&self) -> f64 {
        if self.window.is_empty() {
            return 0.0;
        }
        self.window.iter().sum::<f64>() / self.window.len() as f64
    }
    /// Current window sample variance (n − 1 denominator).
    pub fn variance(&self) -> f64 {
        let n = self.window.len();
        if n < 2 {
            return 0.0;
        }
        let m = self.mean();
        self.window.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n as f64 - 1.0)
    }
    /// Number of observations currently in the window.
    pub fn len(&self) -> usize {
        self.window.len()
    }
    /// Returns `true` if the window is empty.
    pub fn is_empty(&self) -> bool {
        self.window.is_empty()
    }
}
/// Uniform distribution on \[`min`, `max`\].
pub struct UniformDistribution {
    /// Lower bound of the distribution.
    pub min: f64,
    /// Upper bound of the distribution.
    pub max: f64,
}
impl UniformDistribution {
    /// Creates a new [`UniformDistribution`] over \[`min`, `max`\].
    pub fn new(min: f64, max: f64) -> Self {
        Self { min, max }
    }
    /// Probability density function evaluated at `x`.
    pub fn pdf(&self, x: f64) -> f64 {
        if x >= self.min && x <= self.max {
            1.0 / (self.max - self.min)
        } else {
            0.0
        }
    }
    /// Draws one sample from this distribution using the supplied [`StatRng`].
    pub fn sample(&self, rng: &mut StatRng) -> f64 {
        self.min + rng.next_f64() * (self.max - self.min)
    }
}
