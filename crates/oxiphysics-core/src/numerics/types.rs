//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Online (streaming) mean and variance using Welford's algorithm.
#[derive(Debug, Clone, Default)]
pub struct OnlineStats {
    pub(super) n: usize,
    pub(super) mean: f64,
    pub(super) m2: f64,
}
impl OnlineStats {
    /// Create a new empty statistics accumulator.
    pub fn new() -> Self {
        Self::default()
    }
    /// Update with a new data point.
    pub fn update(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }
    /// Current mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Population variance.
    pub fn variance(&self) -> f64 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / self.n as f64
        }
    }
    /// Sample variance.
    pub fn sample_variance(&self) -> f64 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / (self.n - 1) as f64
        }
    }
    /// Standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
    /// Number of samples.
    pub fn count(&self) -> usize {
        self.n
    }
    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: &OnlineStats) {
        if other.n == 0 {
            return;
        }
        let n = self.n + other.n;
        let delta = other.mean - self.mean;
        self.mean = (self.mean * self.n as f64 + other.mean * other.n as f64) / n as f64;
        self.m2 = self.m2 + other.m2 + delta * delta * (self.n * other.n) as f64 / n as f64;
        self.n = n;
    }
}
/// Result of a root-finding operation.
#[derive(Debug, Clone)]
pub struct RootResult {
    /// Approximate root location.
    pub root: f64,
    /// Number of iterations taken.
    pub iterations: usize,
    /// Residual |f(root)|.
    pub residual: f64,
    /// Whether the algorithm converged.
    pub converged: bool,
}
