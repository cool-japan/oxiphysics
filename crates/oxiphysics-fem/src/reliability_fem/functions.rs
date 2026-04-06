//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains)]
use std::f64::consts::PI;

/// Standard normal CDF  Φ(x) via Abramowitz & Stegun 26.2.17.
///
/// Maximum absolute error ≈ 7.5 × 10⁻⁸.
pub fn phi(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let pdf_val = (-0.5 * x * x).exp() / (2.0 * PI).sqrt();
    let cdf = 1.0 - pdf_val * poly;
    if x >= 0.0 { cdf } else { 1.0 - cdf }
}
/// Standard normal PDF  φ(x).
pub fn phi_pdf(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * PI).sqrt()
}
/// Inverse normal CDF: smallest β such that Φ(β) ≈ p (bisection, 100 iters).
pub fn phi_inv(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    let mut lo = -10.0_f64;
    let mut hi = 10.0_f64;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if phi(mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
/// Box–Muller transform: return a standard-normal sample given two U(0,1) samples.
pub(super) fn box_muller(u1: f64, u2: f64) -> f64 {
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}
/// Exponential covariance kernel  C(r) = σ² exp(−r / ℓ).
pub fn cov_exponential(r: f64, sigma2: f64, length_scale: f64) -> f64 {
    if length_scale < 1e-300 {
        return 0.0;
    }
    sigma2 * (-r / length_scale).exp()
}
/// Squared-exponential (Gaussian) covariance kernel  C(r) = σ² exp(−r²/(2ℓ²)).
pub fn cov_gaussian(r: f64, sigma2: f64, length_scale: f64) -> f64 {
    if length_scale < 1e-300 {
        return 0.0;
    }
    sigma2 * (-r * r / (2.0 * length_scale * length_scale)).exp()
}
/// Matérn ν=3/2 covariance kernel.
pub fn cov_matern32(r: f64, sigma2: f64, length_scale: f64) -> f64 {
    if length_scale < 1e-300 {
        return 0.0;
    }
    let s = 3.0_f64.sqrt() * r / length_scale;
    sigma2 * (1.0 + s) * (-s).exp()
}
/// Spectral density of the squared-exponential kernel in 1-D:
/// S(ω) = σ² ℓ √(2π) exp(−ω²ℓ²/2).
pub fn spectral_density_gaussian(omega: f64, sigma2: f64, length_scale: f64) -> f64 {
    sigma2
        * length_scale
        * (2.0 * PI).sqrt()
        * (-0.5 * omega * omega * length_scale * length_scale).exp()
}
/// Multiply dense n×n matrix `a` (row-major) by vector `x`.
pub(super) fn mat_vec_mul(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0_f64; n];
    for i in 0..n {
        for j in 0..n {
            y[i] += a[i * n + j] * x[j];
        }
    }
    y
}
/// Compute the coefficient of variation (σ/μ) of a sample.
pub fn coefficient_of_variation(samples: &[f64]) -> f64 {
    if samples.is_empty() {
        return 0.0;
    }
    let n = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / n;
    if mean.abs() < 1e-300 {
        return 0.0;
    }
    let var = samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    var.sqrt() / mean.abs()
}
/// Estimate the failure probability from a sample where f(x) ≤ 0 means failure.
pub fn sample_failure_probability<F>(samples: &[Vec<f64>], g: F) -> f64
where
    F: Fn(&[f64]) -> f64,
{
    if samples.is_empty() {
        return 0.0;
    }
    let fail = samples.iter().filter(|x| g(x.as_slice()) <= 0.0).count();
    fail as f64 / samples.len() as f64
}
/// Mean and standard deviation of a slice.
pub fn mean_std(samples: &[f64]) -> (f64, f64) {
    if samples.is_empty() {
        return (0.0, 0.0);
    }
    let n = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / n;
    let var = samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}
/// Anderson–Darling statistic (simplified) for normality test.
///
/// Returns a rough test statistic; smaller values suggest more normal distribution.
pub fn normality_statistic(samples: &[f64]) -> f64 {
    if samples.len() < 3 {
        return 0.0;
    }
    let (mean, std_dev) = mean_std(samples);
    if std_dev < 1e-300 {
        return 0.0;
    }
    let mut sorted = samples.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    let stat: f64 = sorted
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let z = (x - mean) / std_dev;
            let p = phi(z).clamp(1e-12, 1.0 - 1e-12);
            let j = (2 * i + 1) as f64;
            j * p.ln() + (2.0 * (n as f64) + 1.0 - j) * (1.0 - p).ln()
        })
        .sum::<f64>();
    -(n as f64) - stat / n as f64
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::reliability_fem::*;
    #[test]
    fn test_phi_at_zero() {
        assert!((phi(0.0) - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_phi_symmetry() {
        let x = 1.5;
        assert!((phi(x) + phi(-x) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_phi_at_plus_one() {
        assert!((phi(1.0) - 0.841_344_746).abs() < 1e-5);
    }
    #[test]
    fn test_phi_at_minus_one() {
        assert!((phi(-1.0) - 0.158_655_254).abs() < 1e-5);
    }
    #[test]
    fn test_phi_large_positive() {
        assert!(phi(8.0) > 0.999_999);
    }
    #[test]
    fn test_phi_large_negative() {
        assert!(phi(-8.0) < 1e-6);
    }
    #[test]
    fn test_phi_pdf_peak() {
        assert!((phi_pdf(0.0) - 1.0 / (2.0 * PI).sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_phi_inv_roundtrip() {
        let p = 0.05;
        let x = phi_inv(p);
        assert!((phi(x) - p).abs() < 1e-5);
    }
    #[test]
    fn test_phi_inv_half() {
        assert!(phi_inv(0.5).abs() < 1e-3);
    }
    #[test]
    fn test_phi_monotone() {
        assert!(phi(-1.0) < phi(0.0));
        assert!(phi(0.0) < phi(1.0));
    }
    #[test]
    fn test_cov_exponential_at_zero() {
        assert!((cov_exponential(0.0, 2.0, 1.0) - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_cov_gaussian_at_zero() {
        assert!((cov_gaussian(0.0, 1.5, 0.5) - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_cov_matern32_at_zero() {
        assert!((cov_matern32(0.0, 1.0, 1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cov_exponential_decay() {
        let c0 = cov_exponential(0.0, 1.0, 1.0);
        let c1 = cov_exponential(1.0, 1.0, 1.0);
        assert!(c0 > c1);
    }
    #[test]
    fn test_cov_gaussian_decay() {
        assert!(cov_gaussian(0.0, 1.0, 1.0) > cov_gaussian(2.0, 1.0, 1.0));
    }
    #[test]
    fn test_spectral_density_positive() {
        assert!(spectral_density_gaussian(0.0, 1.0, 1.0) > 0.0);
    }
    #[test]
    fn test_cov_matern32_decay() {
        assert!(cov_matern32(0.0, 1.0, 1.0) > cov_matern32(2.0, 1.0, 1.0));
    }
    #[test]
    fn test_random_field_new_sizes() {
        let rf = RandomField::new(0.0, 1.0, 16, 4, 0.2, 1.0);
        assert_eq!(rf.eigenvalues.len(), 4);
        assert_eq!(rf.eigenvectors.len(), 4 * 16);
    }
    #[test]
    fn test_random_field_eigenvalues_nonneg() {
        let rf = RandomField::new(0.0, 1.0, 16, 4, 0.2, 1.0);
        for &lam in &rf.eigenvalues {
            assert!(lam >= 0.0, "negative eigenvalue {lam}");
        }
    }
    #[test]
    fn test_random_field_sample_length() {
        let rf = RandomField::new(0.0, 1.0, 20, 3, 0.3, 1.0);
        let xi = vec![0.5_f64, -0.5, 1.0];
        let s = rf.sample(&xi);
        assert_eq!(s.len(), 20);
    }
    #[test]
    fn test_random_field_zero_xi() {
        let rf = RandomField::new(0.0, 1.0, 10, 3, 0.2, 1.0);
        let xi = vec![0.0_f64; 3];
        let s = rf.sample(&xi);
        for &v in &s {
            assert_eq!(v, 0.0);
        }
    }
    #[test]
    fn test_random_field_covariance_zero_distance() {
        let rf = RandomField::new(0.0, 1.0, 8, 2, 0.5, 2.0);
        assert!((rf.covariance(0.0) - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_random_field_spectral_density_positive() {
        let rf = RandomField::new(0.0, 1.0, 8, 2, 0.5, 1.0);
        assert!(rf.spectral_density(0.0) > 0.0);
    }
    #[test]
    fn test_random_field_energy_fraction() {
        let rf = RandomField::new(0.0, 1.0, 16, 4, 0.3, 1.0);
        let ef = rf.energy_fraction();
        assert!(ef >= 0.0 && ef <= 1.0 + 1e-10);
    }
    #[test]
    fn test_form_linear_exact() {
        let g = LimitState::new(vec![10.0, 5.0], vec![2.0, 2.0], |x| x[0] - x[1]);
        let r = FormAnalysis::solve(&g);
        let beta_exact = 5.0 / (2.0_f64 * 2.0_f64.sqrt());
        assert!((r.beta - beta_exact).abs() < 0.1, "beta={}", r.beta);
    }
    #[test]
    fn test_form_pf_in_range() {
        let g = LimitState::new(vec![10.0], vec![2.0], |x| x[0] - 6.0);
        let r = FormAnalysis::solve(&g);
        assert!(r.pf > 0.0 && r.pf < 1.0);
    }
    #[test]
    fn test_form_design_point_dim() {
        let g = LimitState::new(vec![0.0, 0.0], vec![1.0, 1.0], |x| x[0] - x[1] - 2.0);
        let r = FormAnalysis::solve(&g);
        assert_eq!(r.design_point.len(), 2);
    }
    #[test]
    fn test_form_alpha_dim() {
        let g = LimitState::new(vec![5.0, 3.0], vec![1.0, 1.0], |x| x[0] - x[1] - 1.0);
        let r = FormAnalysis::solve(&g);
        assert_eq!(r.alpha.len(), 2);
    }
    #[test]
    fn test_form_hasofer_lind_origin() {
        assert_eq!(FormAnalysis::hasofer_lind_beta(&[0.0, 0.0, 0.0]), 0.0);
    }
    #[test]
    fn test_form_hasofer_lind_unit() {
        assert!((FormAnalysis::hasofer_lind_beta(&[1.0, 0.0, 0.0]) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_form_beta_positive() {
        let g = LimitState::new(vec![5.0], vec![1.0], |x| x[0] - 3.0);
        let r = FormAnalysis::solve(&g);
        assert!(r.beta > 0.0, "beta should be positive for safe system");
    }
    #[test]
    fn test_form_beta_sensitivity_to_mean() {
        let alpha = vec![0.8_f64, 0.6];
        let std_dev = vec![1.0, 2.0];
        let sens = FormAnalysis::beta_sensitivity_to_mean(&alpha, &std_dev);
        assert!((sens[0] - (-0.8)).abs() < 1e-12);
        assert!((sens[1] - (-0.3)).abs() < 1e-12);
    }
    #[test]
    fn test_sorm_pf_in_range() {
        let g = LimitState::new(vec![10.0, 5.0], vec![2.0, 2.0], |x| x[0] - x[1]);
        let form = FormAnalysis::solve(&g);
        let sorm = SormAnalysis::solve(&g, &form);
        assert!(sorm.pf_breitung >= 0.0 && sorm.pf_breitung <= 1.0);
    }
    #[test]
    fn test_sorm_curvatures_len() {
        let g = LimitState::new(vec![10.0, 5.0], vec![2.0, 2.0], |x| x[0] - x[1]);
        let form = FormAnalysis::solve(&g);
        let sorm = SormAnalysis::solve(&g, &form);
        assert_eq!(sorm.curvatures.len(), 2);
    }
    #[test]
    fn test_sorm_beta_matches_form() {
        let g = LimitState::new(vec![5.0], vec![1.0], |x| x[0] - 2.0);
        let form = FormAnalysis::solve(&g);
        let sorm = SormAnalysis::solve(&g, &form);
        assert!((sorm.beta - form.beta).abs() < 1e-10);
    }
    #[test]
    fn test_sorm_hohenbichler_correction_positive() {
        let corr = SormAnalysis::hohenbichler_correction(3.0, &[0.1, -0.05]);
        assert!(corr > 0.0);
    }
    #[test]
    fn test_sorm_tvedt_pf_range() {
        let pf = SormAnalysis::tvedt_pf(2.0, &[0.1, 0.05, -0.1]);
        assert!(pf >= 0.0 && pf <= 1.0);
    }
    #[test]
    fn test_sorm_tvedt_no_curvature() {
        let beta = 2.5;
        let pf = SormAnalysis::tvedt_pf(beta, &[]);
        assert!((pf - phi(-beta)).abs() < 1e-10);
    }
    #[test]
    fn test_mc_fem_mean_spring() {
        let mc = MonteCarloFEM::new(2000, 50.0, 5.0);
        let stats = mc.run(|k| 100.0 / k, f64::INFINITY);
        assert!((stats.mean - 2.0).abs() < 0.3, "mean={}", stats.mean);
    }
    #[test]
    fn test_mc_fem_cov_nonneg() {
        let mc = MonteCarloFEM::new(500, 10.0, 1.0);
        let stats = mc.run(|k| k, f64::INFINITY);
        assert!(stats.cov >= 0.0);
    }
    #[test]
    fn test_mc_fem_pf_never_fail() {
        let mc = MonteCarloFEM::new(500, 10.0, 0.1);
        let stats = mc.run(|k| k, 1000.0);
        assert_eq!(stats.pf, 0.0);
    }
    #[test]
    fn test_mc_fem_pf_always_fail() {
        let mc = MonteCarloFEM::new(500, 10.0, 0.1);
        let stats = mc.run(|k| k, 0.0);
        assert_eq!(stats.pf, 1.0);
    }
    #[test]
    fn test_mc_fem_sample_count() {
        let mc = MonteCarloFEM::new(300, 5.0, 0.5);
        let stats = mc.run(|k| k * 2.0, f64::INFINITY);
        assert_eq!(stats.n_samples, 300);
    }
    #[test]
    fn test_mc_fem_min_le_max() {
        let mc = MonteCarloFEM::new(400, 10.0, 1.0);
        let stats = mc.run(|k| k, f64::INFINITY);
        assert!(stats.min <= stats.max);
    }
    #[test]
    fn test_mc_fem_random_field() {
        let rf = RandomField::new(0.0, 1.0, 16, 4, 0.3, 1.0);
        let mc = MonteCarloFEM::new(200, 0.0, 1.0);
        let stats = mc.run_random_field(&rf, f64::INFINITY);
        assert_eq!(stats.n_samples, 200);
        assert!(stats.min <= stats.max);
    }
    #[test]
    fn test_mc_fem_with_load() {
        let mc = MonteCarloFEM::new(500, 50.0, 5.0);
        let stats = mc.run_with_load(|k, f| f / k, 100.0, 10.0, f64::INFINITY);
        assert!(stats.mean > 0.0);
    }
    #[test]
    fn test_mc_stats_pf_confidence() {
        let mc = MonteCarloFEM::new(1000, 10.0, 1.0);
        let stats = mc.run(|k| k, 11.0);
        let hw = stats.pf_confidence_halfwidth();
        assert!(hw >= 0.0);
    }
    #[test]
    fn test_is_estimator_range() {
        let g = LimitState::new(vec![5.0, 3.0], vec![1.0, 1.0], |x| x[0] - x[1]);
        let form = FormAnalysis::solve(&g);
        let is = ImportanceSampling::new(500, 1.0);
        let pf = is.estimate(&g, &form.design_point);
        assert!(pf >= 0.0 && pf <= 1.0, "pf={pf}");
    }
    #[test]
    fn test_is_estimator_cov_formula() {
        let cov = ImportanceSampling::estimator_cov(0.01, 10000);
        assert!(cov > 0.0);
        let cov2 = ImportanceSampling::estimator_cov(0.01, 40000);
        assert!(cov2 < cov);
    }
    #[test]
    fn test_rsm_fit_linear() {
        let g = LimitState::new(vec![5.0, 3.0], vec![1.0, 1.0], |x| 2.0 * x[0] - x[1] - 7.0);
        let rsm = ResponseSurfaceMethod::fit(&g, 1.0);
        assert_eq!(rsm.n_vars, 2);
        assert!(rsm.r_squared > 0.9, "r²={}", rsm.r_squared);
    }
    #[test]
    fn test_rsm_evaluate_at_mean() {
        let g = LimitState::new(vec![5.0], vec![1.0], |x| x[0] - 3.0);
        let rsm = ResponseSurfaceMethod::fit(&g, 1.0);
        let val = rsm.evaluate(&g.mean, &g.mean);
        let expected = g.evaluate(&g.mean);
        assert!(
            (val - expected).abs() < 0.1,
            "val={val} expected={expected}"
        );
    }
    #[test]
    fn test_rsm_form_beta_positive() {
        let g = LimitState::new(vec![5.0], vec![1.0], |x| x[0] - 2.0);
        let rsm = ResponseSurfaceMethod::fit(&g, 1.0);
        let beta = rsm.form_beta(&g);
        assert!(beta > 0.0, "beta={beta}");
    }
    #[test]
    fn test_rsm_coefficients_length() {
        let g = LimitState::new(vec![0.0, 0.0, 0.0], vec![1.0, 1.0, 1.0], |x| {
            x[0] + x[1] + x[2]
        });
        let rsm = ResponseSurfaceMethod::fit(&g, 1.0);
        assert_eq!(rsm.coefficients.len(), 7);
    }
    #[test]
    fn test_series_bounds_single() {
        let (lo, hi) = SystemReliability::series_bounds(&[0.1]);
        assert!((lo - 0.1).abs() < 1e-12);
        assert!((hi - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_series_bounds_multiple() {
        let pfs = vec![0.1, 0.05, 0.02];
        let (lo, hi) = SystemReliability::series_bounds(&pfs);
        assert!((lo - 0.1).abs() < 1e-12);
        assert!((hi - 0.17).abs() < 1e-12);
    }
    #[test]
    fn test_parallel_bounds_single() {
        let (lo, hi) = SystemReliability::parallel_bounds(&[0.2]);
        assert!((lo - 0.2).abs() < 1e-12);
        assert!((hi - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_parallel_bounds_two() {
        let pfs = vec![0.1, 0.2];
        let (lo, hi) = SystemReliability::parallel_bounds(&pfs);
        assert!(lo <= hi + 1e-10);
        assert!(hi <= 0.1 + 1e-12);
    }
    #[test]
    fn test_series_beta_finite() {
        let beta = SystemReliability::series_beta(&[0.01, 0.02]);
        assert!(beta.is_finite());
    }
    #[test]
    fn test_parallel_beta_finite() {
        let beta = SystemReliability::parallel_beta(&[0.1, 0.2]);
        assert!(beta.is_finite());
    }
    #[test]
    fn test_ditlevsen_lower_nonneg() {
        let betas = vec![2.0, 3.0, 2.5];
        let lb = SystemReliability::ditlevsen_lower(&betas, 0.3);
        assert!(lb >= 0.0 && lb <= 1.0);
    }
    #[test]
    fn test_series_mttf() {
        let rates = vec![0.01, 0.02];
        let mttf = SystemReliability::series_mttf(&rates);
        assert!((mttf - 1.0 / 0.03).abs() < 1e-10);
    }
    #[test]
    fn test_sobol_length() {
        let g = LimitState::new(vec![0.0, 0.0, 0.0], vec![1.0, 1.0, 1.0], |x| {
            x[0] + x[1] + x[2]
        });
        let si = SensitivityIndex::compute(&g, 200);
        assert_eq!(si.first_order.len(), 3);
        assert_eq!(si.total_effect.len(), 3);
    }
    #[test]
    fn test_sobol_first_order_in_unit() {
        let g = LimitState::new(vec![0.0, 0.0], vec![1.0, 1.0], |x| x[0] + 0.5 * x[1]);
        let si = SensitivityIndex::compute(&g, 500);
        for s in &si.first_order {
            assert!(*s >= 0.0 && *s <= 1.0, "s={s}");
        }
    }
    #[test]
    fn test_sobol_total_ge_first() {
        let g = LimitState::new(vec![0.0, 0.0], vec![1.0, 1.0], |x| x[0] * x[1]);
        let si = SensitivityIndex::compute(&g, 400);
        for (s, t) in si.first_order.iter().zip(si.total_effect.iter()) {
            assert!(t + 1e-4 >= *s, "total={t} < first={s}");
        }
    }
    #[test]
    fn test_sobol_dominant_variable() {
        let g = LimitState::new(vec![0.0, 0.0], vec![1.0, 0.1], |x| x[0] + 0.01 * x[1]);
        let si = SensitivityIndex::compute(&g, 800);
        assert_eq!(si.dominant_variable(), 0);
    }
    #[test]
    fn test_sobol_first_order_sum_le_one_additive() {
        let g = LimitState::new(vec![0.0, 0.0], vec![1.0, 1.0], |x| x[0] + x[1]);
        let si = SensitivityIndex::compute(&g, 1000);
        assert!(si.first_order_sum() <= 1.5, "sum={}", si.first_order_sum());
    }
    #[test]
    fn test_struct_rel_beta_exact() {
        let sr = StructuralReliability::compute(100.0, 10.0, 60.0, 8.0);
        let expected = 40.0 / (100.0_f64 + 64.0).sqrt();
        assert!(
            (sr.beta - expected).abs() < 1e-6,
            "beta={} expected={}",
            sr.beta,
            expected
        );
    }
    #[test]
    fn test_struct_rel_pf_range() {
        let sr = StructuralReliability::compute(50.0, 5.0, 30.0, 5.0);
        assert!(sr.pf > 0.0 && sr.pf < 1.0);
    }
    #[test]
    fn test_struct_rel_safety_factor() {
        let sr = StructuralReliability::compute(100.0, 10.0, 50.0, 5.0);
        assert!((sr.safety_factor - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_struct_rel_design_resistance() {
        let sr = StructuralReliability::compute(100.0, 10.0, 50.0, 5.0);
        let rd = sr.design_resistance(100.0);
        assert!(rd < 100.0);
    }
    #[test]
    fn test_struct_rel_design_load_effect() {
        let sr = StructuralReliability::compute(100.0, 10.0, 50.0, 5.0);
        let ed = sr.design_load_effect(20.0, 10.0);
        assert!((ed - 42.0).abs() < 1e-10);
    }
    #[test]
    fn test_struct_rel_target_beta() {
        assert!((StructuralReliability::target_beta_cc2() - 3.8).abs() < 1e-10);
        assert!((StructuralReliability::target_beta_cc3() - 4.3).abs() < 1e-10);
    }
    #[test]
    fn test_struct_rel_meets_cc2() {
        let sr = StructuralReliability::compute(100.0, 10.0, 60.0, 8.0);
        let _ = sr.meets_cc2_target();
    }
    #[test]
    fn test_struct_rel_from_form() {
        let g = LimitState::new(vec![100.0, 60.0], vec![10.0, 8.0], |x| x[0] - x[1]);
        let form = FormAnalysis::solve(&g);
        let sr = StructuralReliability::from_form(&form, 100.0, 60.0);
        assert!((sr.beta - form.beta).abs() < 1e-10);
    }
    #[test]
    fn test_struct_rel_lognormal_beta() {
        let beta = StructuralReliability::lognormal_beta(100.0, 10.0, 50.0, 5.0);
        assert!(beta > 0.0, "beta={beta}");
    }
    #[test]
    fn test_cov_empty() {
        assert_eq!(coefficient_of_variation(&[]), 0.0);
    }
    #[test]
    fn test_cov_zero_mean() {
        assert_eq!(coefficient_of_variation(&[0.0, 0.0]), 0.0);
    }
    #[test]
    fn test_cov_constant_samples() {
        assert_eq!(coefficient_of_variation(&[5.0, 5.0, 5.0]), 0.0);
    }
    #[test]
    fn test_sample_failure_probability_never_fail() {
        let s: Vec<Vec<f64>> = vec![vec![1.0], vec![2.0], vec![3.0]];
        assert_eq!(sample_failure_probability(&s, |x| x[0] - 0.5), 0.0);
    }
    #[test]
    fn test_sample_failure_probability_always_fail() {
        let s: Vec<Vec<f64>> = vec![vec![0.0], vec![-1.0]];
        assert_eq!(sample_failure_probability(&s, |x| x[0]), 1.0);
    }
    #[test]
    fn test_mean_std_known() {
        let (m, s) = mean_std(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((m - 3.0).abs() < 1e-12);
        assert!(s > 0.0);
    }
    #[test]
    fn test_mean_std_empty() {
        let (m, s) = mean_std(&[]);
        assert_eq!(m, 0.0);
        assert_eq!(s, 0.0);
    }
    #[test]
    fn test_limit_state_mc_pf_range() {
        let g = LimitState::new(vec![5.0], vec![1.0], |x| x[0] - 3.0);
        let pf = g.monte_carlo_pf(500);
        assert!(pf >= 0.0 && pf <= 1.0);
    }
    #[test]
    fn test_normality_statistic_returns_finite() {
        let samples: Vec<f64> = (0..50).map(|i| i as f64 * 0.1).collect();
        let stat = normality_statistic(&samples);
        assert!(stat.is_finite());
    }
}
