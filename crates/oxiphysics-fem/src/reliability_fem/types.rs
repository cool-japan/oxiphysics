//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg, clippy::type_complexity)]
#[allow(unused_imports)]
use super::functions::*;
use rand::RngExt;
use std::f64::consts::PI;

/// Karhunen–Loève expansion of a one-dimensional random field on \[0, L\].
///
/// The field is discretised at `n_pts` equally-spaced abscissae; the covariance
/// matrix is assembled, eigendecomposed by the power-iteration method, and the
/// dominant `n_terms` eigenpairs are retained.
#[derive(Debug, Clone)]
pub struct RandomField {
    /// Left endpoint of the spatial domain.
    pub x_min: f64,
    /// Right endpoint of the spatial domain.
    pub x_max: f64,
    /// Number of discretisation points.
    pub n_pts: usize,
    /// Number of KL terms retained.
    pub n_terms: usize,
    /// Correlation length scale  ℓ.
    pub length_scale: f64,
    /// Variance  σ².
    pub variance: f64,
    /// KL eigenvalues λ_k (energy of each mode).
    pub eigenvalues: Vec<f64>,
    /// KL eigenvectors φ_k(x) stored column-major: `eigenvectors[k * n_pts + i]`.
    pub eigenvectors: Vec<f64>,
}
impl RandomField {
    /// Construct a new [`RandomField`] and compute the KL decomposition.
    pub fn new(
        x_min: f64,
        x_max: f64,
        n_pts: usize,
        n_terms: usize,
        length_scale: f64,
        variance: f64,
    ) -> Self {
        assert!(n_pts >= 2, "need at least 2 points");
        let n_terms = n_terms.min(n_pts);
        let dx = (x_max - x_min) / (n_pts - 1) as f64;
        let xs: Vec<f64> = (0..n_pts).map(|i| x_min + i as f64 * dx).collect();
        let n = n_pts;
        let mut cov = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let r = (xs[i] - xs[j]).abs();
                cov[i * n + j] = cov_gaussian(r, variance, length_scale) * dx;
            }
        }
        let mut eigenvalues = vec![0.0_f64; n_terms];
        let mut eigenvectors = vec![0.0_f64; n_terms * n];
        let mut deflated = cov.clone();
        for k in 0..n_terms {
            let mut v: Vec<f64> = (0..n).map(|i| ((i + k + 1) as f64).sin()).collect();
            Self::normalise(&mut v);
            let mut lambda = 0.0_f64;
            for _ in 0..200 {
                let w = mat_vec_mul(&deflated, &v, n);
                let new_lambda: f64 = v.iter().zip(w.iter()).map(|(a, b)| a * b).sum();
                let mut w2 = w;
                Self::normalise(&mut w2);
                let diff: f64 = w2
                    .iter()
                    .zip(v.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                v = w2;
                lambda = new_lambda;
                if diff < 1e-10 {
                    break;
                }
            }
            eigenvalues[k] = lambda.max(0.0);
            for i in 0..n {
                eigenvectors[k * n + i] = v[i];
            }
            for i in 0..n {
                for j in 0..n {
                    deflated[i * n + j] -= lambda * v[i] * v[j];
                }
            }
        }
        Self {
            x_min,
            x_max,
            n_pts,
            n_terms,
            length_scale,
            variance,
            eigenvalues,
            eigenvectors,
        }
    }
    /// Sample a realisation of the random field.
    ///
    /// `xi` is a slice of `n_terms` independent N(0,1) random variables.
    pub fn sample(&self, xi: &[f64]) -> Vec<f64> {
        let n = self.n_pts;
        let mut field = vec![0.0_f64; n];
        for k in 0..self.n_terms.min(xi.len()) {
            let scale = self.eigenvalues[k].sqrt() * xi[k];
            for i in 0..n {
                field[i] += scale * self.eigenvectors[k * n + i];
            }
        }
        field
    }
    /// Energy fraction captured by the retained terms.
    pub fn energy_fraction(&self) -> f64 {
        let total = self.eigenvalues.iter().sum::<f64>();
        if total < 1e-300 {
            return 1.0;
        }
        self.eigenvalues.iter().sum::<f64>() / total
    }
    fn normalise(v: &mut Vec<f64>) {
        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-300 {
            for x in v.iter_mut() {
                *x /= norm;
            }
        }
    }
    /// Covariance at distance `r` using the field's parameters.
    pub fn covariance(&self, r: f64) -> f64 {
        cov_gaussian(r, self.variance, self.length_scale)
    }
    /// Spectral density at angular frequency `omega`.
    pub fn spectral_density(&self, omega: f64) -> f64 {
        spectral_density_gaussian(omega, self.variance, self.length_scale)
    }
    /// Mean absolute value of a sample (for sanity checking).
    pub fn mean_abs_sample(&self, xi: &[f64]) -> f64 {
        let s = self.sample(xi);
        if s.is_empty() {
            return 0.0;
        }
        s.iter().map(|v| v.abs()).sum::<f64>() / s.len() as f64
    }
}
/// First Order Reliability Method (FORM) solver.
///
/// Uses the improved Hasofer–Lind–Rackwitz–Fießler (iHL-RF) algorithm.
pub struct FormAnalysis;
impl FormAnalysis {
    /// Run FORM on limit-state `g` and return a [`FormResult`].
    ///
    /// The iHL-RF update rule is:
    ///   u_{k+1} = (∇g·u_k − g(u_k)) / ‖∇g‖² · ∇g
    pub fn solve(g: &LimitState) -> FormResult {
        let n = g.n_vars();
        let mut u = vec![0.0_f64; n];
        let h = 1e-6_f64;
        let mut iters = 0usize;
        for _iter in 0..200 {
            iters += 1;
            let x = g.u_to_x(&u);
            let gval = g.evaluate(&x);
            let mut grad = vec![0.0_f64; n];
            for i in 0..n {
                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[i] += h * g.std_dev[i];
                xm[i] -= h * g.std_dev[i];
                grad[i] =
                    (g.evaluate(&xp) - g.evaluate(&xm)) / (2.0 * h * g.std_dev[i]) * g.std_dev[i];
            }
            let grad_norm_sq: f64 = grad.iter().map(|v| v * v).sum();
            if grad_norm_sq < 1e-30 {
                break;
            }
            let dot: f64 = grad.iter().zip(u.iter()).map(|(gi, ui)| gi * ui).sum();
            let lambda = (dot - gval) / grad_norm_sq;
            let u_new: Vec<f64> = grad.iter().map(|gi| lambda * gi).collect();
            let diff: f64 = u_new
                .iter()
                .zip(u.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            u = u_new;
            if diff < 1e-9 {
                break;
            }
        }
        let beta = u.iter().map(|v| v * v).sum::<f64>().sqrt();
        let pf = phi(-beta);
        let alpha: Vec<f64> = if beta > 1e-12 {
            u.iter().map(|ui| -ui / beta).collect()
        } else {
            vec![0.0; n]
        };
        FormResult {
            beta,
            pf,
            design_point: u,
            alpha,
            iterations: iters,
        }
    }
    /// Compute the Hasofer–Lind reliability index ‖u‖ from a design point.
    pub fn hasofer_lind_beta(u: &[f64]) -> f64 {
        u.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
    /// Importance vector (direction cosines) α = −u* / β.
    pub fn importance_vector(u: &[f64]) -> Vec<f64> {
        let beta = Self::hasofer_lind_beta(u);
        if beta < 1e-12 {
            return vec![0.0; u.len()];
        }
        u.iter().map(|ui| -ui / beta).collect()
    }
    /// Sensitivity of β to the mean of variable i: ∂β/∂μ_i = −α_i / σ_i.
    pub fn beta_sensitivity_to_mean(alpha: &[f64], std_dev: &[f64]) -> Vec<f64> {
        alpha
            .iter()
            .zip(std_dev.iter())
            .map(|(&ai, &si)| -ai / si.max(1e-300))
            .collect()
    }
    /// Check whether a design point satisfies the FORM convergence criterion.
    pub fn is_design_point(u: &[f64], g: &LimitState) -> bool {
        let x = g.u_to_x(u);
        g.evaluate(&x).abs() < 1e-6
    }
}
/// System reliability bounds for series and parallel systems.
///
/// Uses the Ditlevsen (1979) narrow bounds for the failure probability.
pub struct SystemReliability;
impl SystemReliability {
    /// Probability of failure for a **series** system (weakest-link).
    ///
    /// Lower bound (first-order): `max(P_f_i)`.
    /// Upper bound (first-order): `min(1, Σ P_f_i)`.
    pub fn series_bounds(pf_components: &[f64]) -> (f64, f64) {
        if pf_components.is_empty() {
            return (0.0, 0.0);
        }
        let lower = pf_components.iter().cloned().fold(0.0_f64, f64::max);
        let upper = (pf_components.iter().sum::<f64>()).min(1.0);
        (lower, upper)
    }
    /// Probability of failure for a **parallel** system (all must fail).
    ///
    /// Upper bound (first-order): `min(P_f_i)`.
    /// Lower bound (first-order): `max(0, Σ P_f_i - (n-1))`.
    pub fn parallel_bounds(pf_components: &[f64]) -> (f64, f64) {
        if pf_components.is_empty() {
            return (0.0, 0.0);
        }
        let n = pf_components.len() as f64;
        let upper = pf_components.iter().cloned().fold(f64::INFINITY, f64::min);
        let lower = (pf_components.iter().sum::<f64>() - (n - 1.0)).max(0.0);
        (lower, upper)
    }
    /// Reliability index for a series system (conservative): β = Φ⁻¹(1 − P_f_upper).
    pub fn series_beta(pf_components: &[f64]) -> f64 {
        let (_lower, upper) = Self::series_bounds(pf_components);
        phi_inv(1.0 - upper)
    }
    /// Reliability index for a parallel system: β = Φ⁻¹(1 − P_f_upper).
    pub fn parallel_beta(pf_components: &[f64]) -> f64 {
        let (_lower, upper) = Self::parallel_bounds(pf_components);
        phi_inv(1.0 - upper)
    }
    /// Ditlevsen narrow bounds for series system given pairwise correlation.
    ///
    /// Second-order lower bound:
    ///   P_f ≥ P_f_1 + Σ_{i=2}^{n} max(0, P_f_i − Σ_{j<i} P(F_i ∩ F_j))
    ///
    /// Here we approximate P(F_i ∩ F_j) ≈ phi(-beta_i) * phi(-beta_j) * (1 + rho)
    /// where rho is given.
    pub fn ditlevsen_lower(betas: &[f64], rho: f64) -> f64 {
        if betas.is_empty() {
            return 0.0;
        }
        let pfs: Vec<f64> = betas.iter().map(|&b| phi(-b)).collect();
        let mut total = pfs[0];
        for i in 1..betas.len() {
            let mut inter_sum = 0.0_f64;
            for j in 0..i {
                inter_sum += pfs[i] * pfs[j] * (1.0 + rho).clamp(0.0, 1.0);
            }
            let term = (pfs[i] - inter_sum).max(0.0);
            total += term;
        }
        total.clamp(0.0, 1.0)
    }
    /// Mean failure rate combining independent components in series.
    pub fn series_mean_rate(failure_rates: &[f64]) -> f64 {
        failure_rates.iter().sum()
    }
    /// Mean time to failure (MTTF) for a series system.
    pub fn series_mttf(failure_rates: &[f64]) -> f64 {
        let total_rate = Self::series_mean_rate(failure_rates);
        if total_rate < 1e-300 {
            return f64::INFINITY;
        }
        1.0 / total_rate
    }
}
/// Partial safety factors for loads and resistance (Eurocode-style).
#[derive(Debug, Clone)]
pub struct PartialCoefficients {
    /// Partial factor for the resistance (γ_R).
    pub gamma_r: f64,
    /// Partial factor for permanent loads (γ_G).
    pub gamma_g: f64,
    /// Partial factor for variable loads (γ_Q).
    pub gamma_q: f64,
}
impl PartialCoefficients {
    /// Eurocode default partial coefficients.
    pub fn eurocode_default() -> Self {
        Self {
            gamma_r: 1.0 / 0.8,
            gamma_g: 1.35,
            gamma_q: 1.5,
        }
    }
    /// Custom partial coefficients.
    pub fn new(gamma_r: f64, gamma_g: f64, gamma_q: f64) -> Self {
        Self {
            gamma_r,
            gamma_g,
            gamma_q,
        }
    }
}
/// Limit-state function g(X) with independent normal random variables.
///
/// Failure occurs when g(X) ≤ 0.
pub struct LimitState {
    /// Mean values of the input random variables.
    pub mean: Vec<f64>,
    /// Standard deviations of the input random variables.
    pub std_dev: Vec<f64>,
    pub(super) eval: Box<dyn Fn(&[f64]) -> f64 + Send + Sync>,
}
impl LimitState {
    /// Create a new limit-state function.
    pub fn new(
        mean: Vec<f64>,
        std_dev: Vec<f64>,
        eval: impl Fn(&[f64]) -> f64 + Send + Sync + 'static,
    ) -> Self {
        Self {
            mean,
            std_dev,
            eval: Box::new(eval),
        }
    }
    /// Evaluate g at physical-space point `x`.
    pub fn evaluate(&self, x: &[f64]) -> f64 {
        (self.eval)(x)
    }
    /// Transform standard-normal vector `u` → physical space `x`.
    pub fn u_to_x(&self, u: &[f64]) -> Vec<f64> {
        u.iter()
            .zip(self.mean.iter())
            .zip(self.std_dev.iter())
            .map(|((ui, mi), si)| mi + si * ui)
            .collect()
    }
    /// Transform physical-space vector `x` → standard-normal space `u`.
    pub fn x_to_u(&self, x: &[f64]) -> Vec<f64> {
        x.iter()
            .zip(self.mean.iter())
            .zip(self.std_dev.iter())
            .map(|((xi, mi), si)| (xi - mi) / si)
            .collect()
    }
    /// Number of random variables.
    pub fn n_vars(&self) -> usize {
        self.mean.len()
    }
    /// Estimate failure probability by direct Monte Carlo.
    pub fn monte_carlo_pf(&self, n_samples: usize) -> f64 {
        use rand::RngExt;
        let mut rng = rand::rng();
        let mut failures = 0usize;
        for _ in 0..n_samples {
            let x: Vec<f64> = self
                .mean
                .iter()
                .zip(self.std_dev.iter())
                .map(|(&mi, &si)| {
                    let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    mi + si * box_muller(u1, u2)
                })
                .collect();
            if self.evaluate(&x) <= 0.0 {
                failures += 1;
            }
        }
        failures as f64 / n_samples as f64
    }
}
/// Results of a FORM analysis.
#[derive(Debug, Clone)]
pub struct FormResult {
    /// Hasofer–Lind reliability index β_HL = ‖u*‖.
    pub beta: f64,
    /// First-order failure probability Pf ≈ Φ(−β).
    pub pf: f64,
    /// Most probable failure point (design point) in standard-normal space.
    pub design_point: Vec<f64>,
    /// Sensitivity (direction cosines α = −u*/β).
    pub alpha: Vec<f64>,
    /// Number of iHL-RF iterations to convergence.
    pub iterations: usize,
}
/// Second Order Reliability Method (SORM) solver.
///
/// Computes the principal curvatures of g = 0 at the FORM design point and
/// applies the Breitung (1984) asymptotic correction to the FORM failure
/// probability.
pub struct SormAnalysis;
impl SormAnalysis {
    /// Run SORM on a limit-state function given an existing FORM result.
    ///
    /// The Breitung correction is:
    ///   Pf_SORM ≈ Φ(−β) ∏_i (1 + β κ_i)^{−1/2}
    pub fn solve(g: &LimitState, form: &FormResult) -> SormResult {
        let n = g.n_vars();
        let beta = form.beta;
        let u_star = &form.design_point;
        let h = 1e-5_f64;
        let x_star = g.u_to_x(u_star);
        let gval = g.evaluate(&x_star);
        let mut hess = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let (pp, pm, mp, mm) = {
                    let mut xpp = x_star.clone();
                    let mut xpm = x_star.clone();
                    let mut xmp = x_star.clone();
                    let mut xmm = x_star.clone();
                    xpp[i] += h * g.std_dev[i];
                    xpp[j] += h * g.std_dev[j];
                    xpm[i] += h * g.std_dev[i];
                    xpm[j] -= h * g.std_dev[j];
                    xmp[i] -= h * g.std_dev[i];
                    xmp[j] += h * g.std_dev[j];
                    xmm[i] -= h * g.std_dev[i];
                    xmm[j] -= h * g.std_dev[j];
                    (
                        g.evaluate(&xpp),
                        g.evaluate(&xpm),
                        g.evaluate(&xmp),
                        g.evaluate(&xmm),
                    )
                };
                hess[i * n + j] = (pp - pm - mp + mm) / (4.0 * h * h * g.std_dev[i] * g.std_dev[j]);
            }
        }
        let _ = gval;
        let mut grad = vec![0.0_f64; n];
        for i in 0..n {
            let mut xp = x_star.clone();
            let mut xm = x_star.clone();
            xp[i] += h * g.std_dev[i];
            xm[i] -= h * g.std_dev[i];
            grad[i] = (g.evaluate(&xp) - g.evaluate(&xm)) / (2.0 * h * g.std_dev[i]) * g.std_dev[i];
        }
        let grad_norm = grad.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-15);
        let curvatures: Vec<f64> = (0..n).map(|i| hess[i * n + i] / grad_norm).collect();
        let mut correction_factor = 1.0_f64;
        for &ki in &curvatures {
            let denom = 1.0 + beta * ki;
            if denom > 0.0 {
                correction_factor /= denom.sqrt();
            }
        }
        let pf_form = phi(-beta);
        let pf_breitung = pf_form * correction_factor;
        SormResult {
            beta,
            curvatures,
            pf_breitung: pf_breitung.clamp(0.0, 1.0),
            correction_factor,
        }
    }
    /// Hohenbichler–Rackwitz correction for Pf_FORM.
    pub fn hohenbichler_correction(beta: f64, curvatures: &[f64]) -> f64 {
        let mu = phi_pdf(beta) / phi(-beta).max(1e-300);
        let mut prod = 1.0_f64;
        for &ki in curvatures {
            let denom = 1.0 + mu * ki;
            if denom > 0.0 {
                prod /= denom.sqrt();
            }
        }
        prod
    }
    /// Tvedt three-term approximation (first two terms shown).
    ///
    /// Returns the failure probability using Tvedt's approximation:
    ///   Pf ≈ Φ(−β) ∏(1 + β κ_i)^{−1/2}   (same as Breitung here)
    pub fn tvedt_pf(beta: f64, curvatures: &[f64]) -> f64 {
        let pf_form = phi(-beta);
        let mut prod = 1.0_f64;
        for &ki in curvatures {
            let d = 1.0 + beta * ki;
            if d > 0.0 {
                prod /= d.sqrt();
            }
        }
        (pf_form * prod).clamp(0.0, 1.0)
    }
}
/// Monte Carlo FEM analysis driver.
///
/// Evaluates a user-supplied "FEM response function" f(k) many times where
/// k is drawn from a random field (or scalar random variable), then computes
/// response statistics and failure probabilities.
pub struct MonteCarloFEM {
    /// Number of Monte Carlo samples.
    pub n_samples: usize,
    /// Mean stiffness / material parameter.
    pub mean_k: f64,
    /// Standard deviation of the stiffness.
    pub std_k: f64,
}
impl MonteCarloFEM {
    /// Create a new [`MonteCarloFEM`] driver.
    pub fn new(n_samples: usize, mean_k: f64, std_k: f64) -> Self {
        Self {
            n_samples,
            mean_k,
            std_k,
        }
    }
    /// Run Monte Carlo sampling using a response function `f(k) -> response`.
    pub fn run<F>(&self, response_fn: F, threshold: f64) -> McFemStats
    where
        F: Fn(f64) -> f64,
    {
        let mut rng = rand::rng();
        let mut samples = Vec::with_capacity(self.n_samples);
        for _ in 0..self.n_samples {
            let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
            let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
            let z = box_muller(u1, u2);
            let k = (self.mean_k + self.std_k * z).max(self.mean_k * 1e-9);
            samples.push(response_fn(k));
        }
        Self::compute_stats(&samples, threshold, self.n_samples)
    }
    /// Run Monte Carlo on a random field model.
    pub fn run_random_field(&self, rf: &RandomField, threshold: f64) -> McFemStats {
        let mut rng = rand::rng();
        let n_terms = rf.n_terms;
        let mut samples = Vec::with_capacity(self.n_samples);
        for _ in 0..self.n_samples {
            let xi: Vec<f64> = (0..n_terms)
                .map(|_| {
                    let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    box_muller(u1, u2)
                })
                .collect();
            let field = rf.sample(&xi);
            let max_val = field.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            samples.push(max_val);
        }
        Self::compute_stats(&samples, threshold, self.n_samples)
    }
    /// Run Monte Carlo with multiple independent loads.
    ///
    /// `load_fn(k, load)` computes response given material stiffness `k` and load `load`.
    /// Both k and load are drawn independently from N(mean_k, std_k²).
    pub fn run_with_load<F>(
        &self,
        load_fn: F,
        mean_load: f64,
        std_load: f64,
        threshold: f64,
    ) -> McFemStats
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut rng = rand::rng();
        let mut samples = Vec::with_capacity(self.n_samples);
        for _ in 0..self.n_samples {
            let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
            let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
            let u3: f64 = rng.random_range(f64::EPSILON..1.0_f64);
            let u4: f64 = rng.random_range(0.0_f64..1.0_f64);
            let zk = box_muller(u1, u2);
            let zl = box_muller(u3, u4);
            let k = (self.mean_k + self.std_k * zk).max(self.mean_k * 1e-9);
            let load = mean_load + std_load * zl;
            samples.push(load_fn(k, load));
        }
        Self::compute_stats(&samples, threshold, self.n_samples)
    }
    fn compute_stats(samples: &[f64], threshold: f64, n_samples: usize) -> McFemStats {
        let n = samples.len() as f64;
        let mean = if n > 0.0 {
            samples.iter().sum::<f64>() / n
        } else {
            0.0
        };
        let var = if n > 0.0 {
            samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n
        } else {
            0.0
        };
        let std_dev = var.sqrt();
        let cov = if mean.abs() > 1e-300 {
            std_dev / mean.abs()
        } else {
            0.0
        };
        let min = samples.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let failures = samples.iter().filter(|&&v| v > threshold).count();
        let pf = failures as f64 / n_samples as f64;
        McFemStats {
            n_samples,
            mean,
            std_dev,
            cov,
            min,
            max,
            pf,
        }
    }
}
/// Results of a SORM analysis.
#[derive(Debug, Clone)]
pub struct SormResult {
    /// FORM reliability index β (same as from FORM).
    pub beta: f64,
    /// Principal curvatures κ_i of the limit-state surface at the design point.
    pub curvatures: Vec<f64>,
    /// SORM probability of failure (Breitung approximation).
    pub pf_breitung: f64,
    /// Hohenbichler–Rackwitz correction factor.
    pub correction_factor: f64,
}
/// Statistics from a Monte Carlo FEM analysis.
#[derive(Debug, Clone)]
pub struct McFemStats {
    /// Number of samples.
    pub n_samples: usize,
    /// Sample mean of the response.
    pub mean: f64,
    /// Sample standard deviation of the response.
    pub std_dev: f64,
    /// Coefficient of variation (σ/μ).
    pub cov: f64,
    /// Minimum response value observed.
    pub min: f64,
    /// Maximum response value observed.
    pub max: f64,
    /// Estimated failure probability (response > threshold).
    pub pf: f64,
}
impl McFemStats {
    /// 95 % confidence interval half-width for the failure probability estimator.
    ///
    /// `Δ = 1.96 sqrt(pf (1-pf) / n)`
    pub fn pf_confidence_halfwidth(&self) -> f64 {
        let n = self.n_samples as f64;
        if n < 1.0 {
            return 0.0;
        }
        1.96 * (self.pf * (1.0 - self.pf) / n).sqrt()
    }
    /// Coefficient of variation of the failure probability estimator.
    pub fn pf_cov(&self) -> f64 {
        if self.pf < 1e-300 {
            return 0.0;
        }
        ((1.0 - self.pf) / (self.pf * self.n_samples as f64)).sqrt()
    }
}
/// Importance Sampling reliability estimator.
///
/// Samples from a proposal distribution centred at the FORM design point u*
/// and weights samples by the likelihood ratio to estimate P_f efficiently.
pub struct ImportanceSampling {
    /// Number of IS samples.
    pub n_samples: usize,
    /// Standard deviation of the IS proposal (isotropic Gaussian).
    pub proposal_std: f64,
}
impl ImportanceSampling {
    /// Create a new importance sampling estimator.
    pub fn new(n_samples: usize, proposal_std: f64) -> Self {
        Self {
            n_samples,
            proposal_std,
        }
    }
    /// Estimate the failure probability.
    ///
    /// The IS estimator is:
    ///   P_f = (1/N) Σ I(g(x_i) ≤ 0) · p(x_i) / h(x_i)
    /// where h is the proposal (N(u*, σ_IS²) in standard-normal space) and p = φ^n.
    pub fn estimate(&self, g: &LimitState, design_point: &[f64]) -> f64 {
        let mut rng = rand::rng();
        let n = g.n_vars();
        let mut pf_sum = 0.0_f64;
        for _ in 0..self.n_samples {
            let u: Vec<f64> = (0..n)
                .map(|i| {
                    let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    let z = box_muller(u1, u2);
                    design_point[i] + self.proposal_std * z
                })
                .collect();
            let x = g.u_to_x(&u);
            if g.evaluate(&x) <= 0.0 {
                let ln_p: f64 = -0.5 * u.iter().map(|v| v * v).sum::<f64>();
                let ln_h: f64 = -0.5
                    * u.iter()
                        .zip(design_point.iter())
                        .map(|(ui, di)| (ui - di).powi(2))
                        .sum::<f64>()
                    / (self.proposal_std * self.proposal_std)
                    - n as f64 * (2.0 * PI * self.proposal_std * self.proposal_std).ln() / 2.0
                    + n as f64 * (2.0 * PI).ln() / 2.0;
                let weight = (ln_p - ln_h).exp();
                pf_sum += weight;
            }
        }
        (pf_sum / self.n_samples as f64).clamp(0.0, 1.0)
    }
    /// Coefficient of variation of the IS estimator given a preliminary pf estimate.
    pub fn estimator_cov(pf_estimate: f64, n_samples: usize) -> f64 {
        if pf_estimate < 1e-300 {
            return 0.0;
        }
        ((1.0 - pf_estimate) / (pf_estimate * n_samples as f64)).sqrt()
    }
}
/// Polynomial response surface for limit state approximation.
///
/// Fits a quadratic response surface g̃(x) = a₀ + Σ aᵢ xᵢ + Σ bᵢ xᵢ²
/// to limit-state function evaluations at a central composite design.
#[derive(Debug, Clone)]
pub struct ResponseSurfaceMethod {
    /// Number of input variables.
    pub n_vars: usize,
    /// Polynomial coefficients: \[a₀, a₁, …, a_n, b₁, …, b_n\].
    pub coefficients: Vec<f64>,
    /// R² goodness of fit.
    pub r_squared: f64,
}
impl ResponseSurfaceMethod {
    /// Fit a quadratic response surface to the limit-state function `g`.
    ///
    /// Uses a central composite design with step size `h` from the mean.
    /// Fits 2n+1 coefficients via least squares (simplified: one pass).
    pub fn fit(g: &LimitState, h: f64) -> Self {
        let n = g.n_vars();
        let mean = &g.mean;
        let n_pts = 2 * n + 1;
        let mut x_pts: Vec<Vec<f64>> = Vec::with_capacity(n_pts);
        let mut y_pts: Vec<f64> = Vec::with_capacity(n_pts);
        x_pts.push(mean.clone());
        y_pts.push(g.evaluate(mean));
        for i in 0..n {
            let mut xp = mean.clone();
            let mut xm = mean.clone();
            xp[i] += h * g.std_dev[i];
            xm[i] -= h * g.std_dev[i];
            x_pts.push(xp);
            y_pts.push(g.evaluate(&x_pts[x_pts.len() - 1]));
            x_pts.push(xm);
            y_pts.push(g.evaluate(&x_pts[x_pts.len() - 1]));
        }
        let g0 = y_pts[0];
        let mut a_coeffs = vec![0.0_f64; n];
        let mut b_coeffs = vec![0.0_f64; n];
        for i in 0..n {
            let gp = y_pts[1 + 2 * i];
            let gm = y_pts[2 + 2 * i];
            let hi = h * g.std_dev[i];
            a_coeffs[i] = (gp - gm) / (2.0 * hi);
            b_coeffs[i] = (gp + gm - 2.0 * g0) / (hi * hi);
        }
        let mut coefficients = vec![g0];
        coefficients.extend_from_slice(&a_coeffs);
        coefficients.extend_from_slice(&b_coeffs);
        let y_pred: Vec<f64> = x_pts
            .iter()
            .map(|x| {
                let mut val = g0;
                for i in 0..n {
                    let dx = x[i] - mean[i];
                    val += a_coeffs[i] * dx + b_coeffs[i] * dx * dx;
                }
                val
            })
            .collect();
        let y_mean = y_pts.iter().sum::<f64>() / y_pts.len() as f64;
        let ss_res: f64 = y_pts
            .iter()
            .zip(y_pred.iter())
            .map(|(y, yp)| (y - yp).powi(2))
            .sum();
        let ss_tot: f64 = y_pts.iter().map(|y| (y - y_mean).powi(2)).sum();
        let r_squared = if ss_tot > 1e-300 {
            1.0 - ss_res / ss_tot
        } else {
            1.0
        };
        Self {
            n_vars: n,
            coefficients,
            r_squared,
        }
    }
    /// Evaluate the fitted response surface at point `x`.
    pub fn evaluate(&self, x: &[f64], mean: &[f64]) -> f64 {
        let n = self.n_vars;
        let mut val = self.coefficients[0];
        for i in 0..n {
            let dx = x[i] - mean[i];
            val += self.coefficients[1 + i] * dx;
            val += self.coefficients[1 + n + i] * dx * dx;
        }
        val
    }
    /// Run FORM on the fitted response surface.
    ///
    /// The surface is evaluated as a `LimitState` for use with `FormAnalysis`.
    pub fn form_beta(&self, g_original: &LimitState) -> f64 {
        let mean = g_original.mean.clone();
        let coeffs = self.coefficients.clone();
        let n = self.n_vars;
        let mean_clone = mean.clone();
        let approx = LimitState::new(mean.clone(), g_original.std_dev.clone(), move |x| {
            let mut val = coeffs[0];
            for i in 0..n {
                let dx = x[i] - mean_clone[i];
                val += coeffs[1 + i] * dx;
                val += coeffs[1 + n + i] * dx * dx;
            }
            val
        });
        FormAnalysis::solve(&approx).beta
    }
}
/// Sobol variance-based global sensitivity indices.
#[derive(Debug, Clone)]
pub struct SensitivityIndex {
    /// First-order Sobol indices S_i.
    pub first_order: Vec<f64>,
    /// Total-effect Sobol indices S_Ti.
    pub total_effect: Vec<f64>,
    /// Total output variance V\[Y\].
    pub total_variance: f64,
}
impl SensitivityIndex {
    /// Estimate Sobol indices using the Saltelli (2002) estimator.
    pub fn compute(g: &LimitState, n: usize) -> Self {
        let mut rng = rand::rng();
        let k = g.n_vars();
        let mut mat_a: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut mat_b: Vec<Vec<f64>> = Vec::with_capacity(n);
        for _ in 0..n {
            let row_a: Vec<f64> = (0..k)
                .map(|i| {
                    let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    let z = box_muller(u1, u2);
                    g.mean[i] + g.std_dev[i] * z
                })
                .collect();
            let row_b: Vec<f64> = (0..k)
                .map(|i| {
                    let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
                    let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                    let z = box_muller(u1, u2);
                    g.mean[i] + g.std_dev[i] * z
                })
                .collect();
            mat_a.push(row_a);
            mat_b.push(row_b);
        }
        let ya: Vec<f64> = mat_a.iter().map(|x| g.evaluate(x)).collect();
        let yb: Vec<f64> = mat_b.iter().map(|x| g.evaluate(x)).collect();
        let f0 = ya.iter().sum::<f64>() / n as f64;
        let var_y = ya.iter().map(|v| (v - f0).powi(2)).sum::<f64>() / n as f64;
        let mut first_order = vec![0.0_f64; k];
        let mut total_effect = vec![0.0_f64; k];
        if var_y > 1e-30 {
            for j in 0..k {
                let ab_j: Vec<Vec<f64>> = (0..n)
                    .map(|r| {
                        let mut row = mat_a[r].clone();
                        row[j] = mat_b[r][j];
                        row
                    })
                    .collect();
                let y_ab_j: Vec<f64> = ab_j.iter().map(|x| g.evaluate(x)).collect();
                let si: f64 = ya
                    .iter()
                    .zip(y_ab_j.iter())
                    .zip(yb.iter())
                    .map(|((a, ab), b)| b * (ab - a))
                    .sum::<f64>()
                    / (n as f64 * var_y);
                let sti: f64 = ya
                    .iter()
                    .zip(y_ab_j.iter())
                    .map(|(a, ab)| (a - ab).powi(2))
                    .sum::<f64>()
                    / (2.0 * n as f64 * var_y);
                first_order[j] = si.clamp(0.0, 1.0);
                total_effect[j] = sti.clamp(0.0, 1.0);
            }
        }
        Self {
            first_order,
            total_effect,
            total_variance: var_y,
        }
    }
    /// Dominant variable: index of the largest first-order index.
    pub fn dominant_variable(&self) -> usize {
        self.first_order
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Sum of first-order indices (should be ≤ 1 for additive models).
    pub fn first_order_sum(&self) -> f64 {
        self.first_order.iter().sum()
    }
}
/// Structural reliability assessment container.
#[derive(Debug, Clone)]
pub struct StructuralReliability {
    /// Probability of failure P_f.
    pub pf: f64,
    /// Reliability index β = −Φ⁻¹(P_f).
    pub beta: f64,
    /// Central safety factor μ_R / μ_S.
    pub safety_factor: f64,
    /// Partial safety coefficients.
    pub partial_coefficients: PartialCoefficients,
}
impl StructuralReliability {
    /// Compute structural reliability for a resistance–load problem.
    ///
    /// Assumes R ~ N(mu_r, sigma_r²) and S ~ N(mu_s, sigma_s²) (independent).
    pub fn compute(mu_r: f64, sigma_r: f64, mu_s: f64, sigma_s: f64) -> Self {
        let beta_denom = (sigma_r * sigma_r + sigma_s * sigma_s).sqrt();
        let beta = if beta_denom > 1e-300 {
            (mu_r - mu_s) / beta_denom
        } else if mu_r > mu_s {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
        let pf = phi(-beta);
        let safety_factor = if mu_s.abs() > 1e-300 {
            mu_r / mu_s
        } else {
            f64::INFINITY
        };
        Self {
            pf,
            beta,
            safety_factor,
            partial_coefficients: PartialCoefficients::eurocode_default(),
        }
    }
    /// Compute from an existing FORM result.
    pub fn from_form(form: &FormResult, mu_r: f64, mu_s: f64) -> Self {
        let safety_factor = if mu_s.abs() > 1e-300 {
            mu_r / mu_s
        } else {
            f64::INFINITY
        };
        Self {
            pf: form.pf,
            beta: form.beta,
            safety_factor,
            partial_coefficients: PartialCoefficients::eurocode_default(),
        }
    }
    /// Design resistance R_d = μ_R / γ_R.
    pub fn design_resistance(&self, mu_r: f64) -> f64 {
        mu_r / self.partial_coefficients.gamma_r
    }
    /// Design load effect E_d = γ_G · G_k + γ_Q · Q_k.
    pub fn design_load_effect(&self, g_k: f64, q_k: f64) -> f64 {
        self.partial_coefficients.gamma_g * g_k + self.partial_coefficients.gamma_q * q_k
    }
    /// Target reliability index β_target for consequence class CC2 (β = 3.8).
    pub fn target_beta_cc2() -> f64 {
        3.8
    }
    /// Target reliability index β_target for consequence class CC3 (β = 4.3).
    pub fn target_beta_cc3() -> f64 {
        4.3
    }
    /// Check whether this structure meets the consequence-class 2 target.
    pub fn meets_cc2_target(&self) -> bool {
        self.beta >= Self::target_beta_cc2()
    }
    /// Lognormal reliability index β_LN (Ang-Tang method).
    ///
    /// For R ~ LN(μ_R, σ_R) and S ~ LN(μ_S, σ_S):
    ///   β_LN ≈ ln(μ_R / μ_S) / sqrt(V_R² + V_S²)
    /// where V = σ/μ is the CoV.
    pub fn lognormal_beta(mu_r: f64, sigma_r: f64, mu_s: f64, sigma_s: f64) -> f64 {
        let vr = (sigma_r / mu_r.abs()).max(1e-9);
        let vs = (sigma_s / mu_s.abs()).max(1e-9);
        let denom = (vr * vr + vs * vs).sqrt();
        if denom < 1e-300 {
            return f64::INFINITY;
        }
        (mu_r / mu_s.max(1e-300)).abs().ln() / denom
    }
}
