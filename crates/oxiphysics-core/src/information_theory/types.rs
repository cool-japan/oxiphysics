//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Encapsulates a discrete probability distribution with entropy utilities.
pub struct DiscreteDistribution {
    /// Normalised probability vector.
    pub probs: Vec<f64>,
}
impl DiscreteDistribution {
    /// Construct from raw (possibly unnormalised) weights.
    ///
    /// Panics if the sum is <= 0.
    pub fn from_weights(weights: &[f64]) -> Self {
        let sum: f64 = weights.iter().sum();
        assert!(
            sum > 0.0,
            "DiscreteDistribution: weights must sum to a positive value"
        );
        let probs = weights.iter().map(|&w| w / sum).collect();
        DiscreteDistribution { probs }
    }
    /// Shannon entropy in bits.
    pub fn entropy_bits(&self) -> f64 {
        entropy_bits(&self.probs)
    }
    /// Shannon entropy in nats.
    pub fn entropy_nats(&self) -> f64 {
        entropy_nats(&self.probs)
    }
    /// Renyi entropy H_alpha in bits.
    pub fn renyi(&self, alpha: f64) -> f64 {
        renyi_entropy(&self.probs, alpha)
    }
    /// Min-entropy H_inf in bits.
    pub fn min_entropy(&self) -> f64 {
        min_entropy(&self.probs)
    }
    /// Collision entropy H_2 in bits.
    pub fn collision_entropy(&self) -> f64 {
        collision_entropy(&self.probs)
    }
    /// KL divergence to another distribution.
    pub fn kl_to(&self, other: &DiscreteDistribution) -> f64 {
        kl_divergence(&self.probs, &other.probs)
    }
    /// Huffman expected code length.
    pub fn huffman_expected_length(&self) -> f64 {
        huffman_expected_length(&self.probs)
    }
    /// Fisher-Rao distance to another distribution.
    pub fn fisher_rao_to(&self, other: &DiscreteDistribution) -> f64 {
        fisher_rao_distance(&self.probs, &other.probs)
    }
}
/// Gaussian information geometry on the (mu, sigma) half-plane.
pub struct GaussianInfoGeometry {
    /// Mean parameter.
    pub mu: f64,
    /// Standard deviation.
    pub sigma: f64,
}
impl GaussianInfoGeometry {
    /// Construct from mean and standard deviation.
    pub fn new(mu: f64, sigma: f64) -> Self {
        GaussianInfoGeometry { mu, sigma }
    }
    /// Fisher information matrix at the current point.
    ///
    /// Returns \[1/sigma^2, 0; 0, 2/sigma^2\] as row-major array.
    pub fn metric(&self) -> [f64; 4] {
        fisher_matrix_gaussian(self.sigma)
    }
    /// Fisher-Rao distance to another Gaussian.
    pub fn distance_to(&self, other: &GaussianInfoGeometry) -> f64 {
        fisher_rao_gaussian(self.mu, self.sigma, other.mu, other.sigma)
    }
    /// Geodesic interpolation between two Gaussians at parameter t in \[0,1\].
    ///
    /// Uses linear interpolation in the natural parameter space as an approximation.
    pub fn geodesic_interp(&self, other: &GaussianInfoGeometry, t: f64) -> GaussianInfoGeometry {
        let eta1_a = self.mu / (self.sigma * self.sigma);
        let eta2_a = -0.5 / (self.sigma * self.sigma);
        let eta1_b = other.mu / (other.sigma * other.sigma);
        let eta2_b = -0.5 / (other.sigma * other.sigma);
        let eta1 = (1.0 - t) * eta1_a + t * eta1_b;
        let eta2 = (1.0 - t) * eta2_a + t * eta2_b;
        let sigma = (-0.5 / eta2).sqrt();
        let mu = eta1 * sigma * sigma;
        GaussianInfoGeometry { mu, sigma }
    }
}
/// Rate-distortion calculator for common source models.
pub struct RateDistortion {
    /// Source variance (for Gaussian source).
    pub sigma_sq: f64,
}
impl RateDistortion {
    /// Construct for a Gaussian source with variance `sigma_sq`.
    pub fn new(sigma_sq: f64) -> Self {
        RateDistortion { sigma_sq }
    }
    /// Gaussian rate-distortion function R(D) = 1/2 log2(sigma^2/D).
    pub fn gaussian_rd(&self, d: f64) -> f64 {
        gaussian_rate_distortion(self.sigma_sq, d)
    }
    /// Distortion at a given rate R: D(R) = sigma^2 * 2^(-2R).
    pub fn gaussian_distortion_at_rate(&self, rate: f64) -> f64 {
        gaussian_distortion_at_rate(self.sigma_sq, rate)
    }
    /// Rate-distortion curve as (D, R) pairs for uniform sampling of D.
    pub fn rd_curve(&self, n_points: usize) -> Vec<(f64, f64)> {
        let mut curve = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let d = self.sigma_sq * (i + 1) as f64 / (n_points + 1) as f64;
            let r = self.gaussian_rd(d);
            curve.push((d, r));
        }
        curve
    }
}
/// Fisher information geometry for the categorical family.
pub struct FisherGeometry {
    /// Current probability vector (point on the simplex).
    pub probs: Vec<f64>,
}
impl FisherGeometry {
    /// Construct from a probability vector.
    pub fn new(probs: Vec<f64>) -> Self {
        FisherGeometry { probs }
    }
    /// Diagonal Fisher metric at the current point.
    pub fn metric(&self) -> Vec<f64> {
        fisher_metric_categorical(&self.probs)
    }
    /// Geodesic distance on the statistical manifold (Fisher-Rao distance).
    ///
    /// d(P, Q) = 2 * arccos(sum sqrt(p_i * q_i)).
    pub fn geodesic_distance(&self, other: &[f64]) -> f64 {
        fisher_rao_distance(&self.probs, other)
    }
    /// Natural gradient step toward gradient direction.
    pub fn natural_step(&self, grad: &[f64], lr: f64) -> Vec<f64> {
        natural_gradient_step_categorical(&self.probs, grad, lr)
    }
    /// Exponential map at the current point.
    pub fn exp_map(&self, v: &[f64], t: f64) -> Vec<f64> {
        exponential_map_simplex(&self.probs, v, t)
    }
    /// Scalar curvature of the probability simplex at the current point.
    ///
    /// For a k-simplex (k+1 outcomes), the scalar curvature under Fisher metric
    /// is constant: R = k(k-1)/4.
    pub fn scalar_curvature(&self) -> f64 {
        let k = self.probs.len() as f64 - 1.0;
        if k <= 0.0 { 0.0 } else { k * (k - 1.0) / 4.0 }
    }
    /// Volume element sqrt(det(g)) at the current point.
    ///
    /// For the Fisher metric on the simplex: sqrt(det(g)) = prod(1/sqrt(p_i)).
    pub fn volume_element(&self) -> f64 {
        self.probs
            .iter()
            .map(|&p| if p > 0.0 { 1.0 / p.sqrt() } else { 0.0 })
            .product()
    }
}
/// AIC/BIC model selector.
pub struct ModelSelector {
    /// List of (num_parameters, log_likelihood) for each candidate model.
    pub models: Vec<(usize, f64)>,
    /// Number of observations.
    pub n_samples: usize,
}
impl ModelSelector {
    /// Construct with the sample size.
    pub fn new(n_samples: usize) -> Self {
        ModelSelector {
            models: Vec::new(),
            n_samples,
        }
    }
    /// Register a model candidate.
    pub fn add(&mut self, k_params: usize, log_likelihood: f64) {
        self.models.push((k_params, log_likelihood));
    }
    /// Index of the best model by AIC (lowest).
    pub fn best_aic(&self) -> Option<usize> {
        select_by_aic(&self.models)
    }
    /// Index of the best model by BIC (lowest).
    pub fn best_bic(&self) -> Option<usize> {
        self.models
            .iter()
            .enumerate()
            .min_by(|a, b| {
                let bic_a = bic(a.1.0, self.n_samples, a.1.1);
                let bic_b = bic(b.1.0, self.n_samples, b.1.1);
                bic_a
                    .partial_cmp(&bic_b)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
    }
    /// AIC values for all registered models.
    pub fn aic_values(&self) -> Vec<f64> {
        self.models.iter().map(|&(k, ll)| aic(k, ll)).collect()
    }
    /// BIC values for all registered models.
    pub fn bic_values(&self) -> Vec<f64> {
        self.models
            .iter()
            .map(|&(k, ll)| bic(k, self.n_samples, ll))
            .collect()
    }
    /// MDL values for all registered models.
    pub fn mdl_values(&self) -> Vec<f64> {
        self.models
            .iter()
            .map(|&(k, ll)| mdl(k, self.n_samples, ll))
            .collect()
    }
}
/// Channel with a stochastic transition matrix.
pub struct DiscreteChannel {
    /// Transition matrix: `transition[i][j]` = P(output j | input i).
    pub transition: Vec<Vec<f64>>,
}
impl DiscreteChannel {
    /// Construct from a row-stochastic matrix.
    pub fn new(transition: Vec<Vec<f64>>) -> Self {
        DiscreteChannel { transition }
    }
    /// Channel capacity (Blahut-Arimoto, 200 iterations), in bits.
    pub fn capacity(&self) -> f64 {
        channel_capacity_blahut(&self.transition)
    }
    /// Symmetric capacity via Shannon-Hartley (for AWGN interpretation).
    ///
    /// Treats `bandwidth_hz` and `snr` as channel parameters.
    pub fn awgn_capacity(bandwidth_hz: f64, snr: f64) -> f64 {
        shannon_hartley_capacity(bandwidth_hz, snr)
    }
    /// Number of input symbols.
    pub fn n_inputs(&self) -> usize {
        self.transition.len()
    }
    /// Number of output symbols.
    pub fn n_outputs(&self) -> usize {
        if self.transition.is_empty() {
            0
        } else {
            self.transition[0].len()
        }
    }
}
