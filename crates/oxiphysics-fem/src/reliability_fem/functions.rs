//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

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
        assert!((0.0..=1.0 + 1e-10).contains(&ef));
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
        assert!((0.0..=1.0).contains(&sorm.pf_breitung));
    }
    #[test]
    fn test_sorm_curvatures_len() {
        // n variables give n - 1 principal curvatures (the old code returned n
        // diagonal Hessian entries).
        let g = LimitState::new(vec![10.0, 5.0], vec![2.0, 2.0], |x| x[0] - x[1]);
        let form = FormAnalysis::solve(&g);
        let sorm = SormAnalysis::solve(&g, &form);
        assert_eq!(sorm.curvatures.len(), 1);
        assert!(
            sorm.curvatures[0].abs() < 1e-6,
            "linear g: {:?}",
            sorm.curvatures
        );
        assert!((sorm.correction_factor - 1.0).abs() < 1e-6);
        let single = LimitState::new(vec![5.0], vec![1.0], |x| x[0] - 2.0);
        let single_form = FormAnalysis::solve(&single);
        assert!(
            SormAnalysis::solve(&single, &single_form)
                .curvatures
                .is_empty()
        );
        let three = LimitState::new(vec![4.0, 5.0, 6.0], vec![1.0, 1.0, 1.0], |x| {
            x[0] + x[1] + x[2] - 3.0
        });
        let three_form = FormAnalysis::solve(&three);
        assert_eq!(SormAnalysis::solve(&three, &three_form).curvatures.len(), 2);
    }
    /// A FORM result with the design point `u_star` (no iteration performed).
    fn form_at(u_star: Vec<f64>) -> FormResult {
        let beta = u_star.iter().map(|v| v * v).sum::<f64>().sqrt();
        let alpha = u_star.iter().map(|u| -u / beta).collect();
        FormResult {
            beta,
            pf: phi(-beta),
            design_point: u_star,
            alpha,
            iterations: 0,
        }
    }
    #[test]
    fn test_sorm_recovers_axis_aligned_paraboloid_curvatures_in_u_space() {
        // g = beta - u_3 + (k1 u_1^2 + k2 u_2^2) / 2 with u_i = (x_i - mu_i) / sigma_i.
        // Non-unit sigma and non-zero mu: the previous code divided the Hessian
        // by sigma_i sigma_j in x-space but multiplied the gradient by sigma,
        // mixing the two spaces.
        let beta = 2.5;
        let (k1, k2) = (0.05, 0.2);
        let mean = vec![10.0, -3.0, 2.0];
        let sigma = vec![2.0, 0.5, 1.5];
        let (m, s) = (mean.clone(), sigma.clone());
        let g = LimitState::new(mean, sigma, move |x| {
            let u: Vec<f64> = x
                .iter()
                .zip(m.iter().zip(s.iter()))
                .map(|(xi, (mi, si))| (xi - mi) / si)
                .collect();
            beta - u[2] + 0.5 * (k1 * u[0] * u[0] + k2 * u[1] * u[1])
        });
        let sorm = SormAnalysis::solve(&g, &form_at(vec![0.0, 0.0, beta]));
        assert_eq!(sorm.curvatures.len(), 2);
        assert!(
            (sorm.curvatures[0] - k1).abs() < 2e-6,
            "{:?}",
            sorm.curvatures
        );
        assert!(
            (sorm.curvatures[1] - k2).abs() < 2e-6,
            "{:?}",
            sorm.curvatures
        );
        let factor = 1.0 / ((1.0 + beta * k1) * (1.0 + beta * k2)).sqrt();
        assert!((sorm.correction_factor - factor).abs() < 1e-5);
        assert!((sorm.pf_breitung - phi(-beta) * factor).abs() < 1e-8);
    }
    #[test]
    fn test_sorm_equal_curvatures_recover_kappa_for_every_dimension() {
        // g = beta - u_n + (kappa / 2) sum_{i<n} u_i^2 has n - 1 curvatures kappa.
        let (beta, kappa) = (2.0, 0.15);
        for n in 2..=5 {
            let mean: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
            let sigma: Vec<f64> = (0..n).map(|i| 0.5 + 0.3 * i as f64).collect();
            let (m, s) = (mean.clone(), sigma.clone());
            let g = LimitState::new(mean, sigma, move |x| {
                let u: Vec<f64> = (0..x.len()).map(|i| (x[i] - m[i]) / s[i]).collect();
                let last = u.len() - 1;
                beta - u[last] + 0.5 * kappa * u[..last].iter().map(|v| v * v).sum::<f64>()
            });
            let mut u_star = vec![0.0; n];
            u_star[n - 1] = beta;
            let sorm = SormAnalysis::solve(&g, &form_at(u_star));
            assert_eq!(sorm.curvatures.len(), n - 1);
            for k in &sorm.curvatures {
                assert!((k - kappa).abs() < 2e-6, "n = {n}: {:?}", sorm.curvatures);
            }
        }
    }
    #[test]
    fn test_sorm_rotates_the_hessian_into_the_tangent_plane() {
        // The design-point direction a = (1, 2, 2) / 3 is not a coordinate axis.
        // G = beta - a.u + [k1 (t1.u)^2 + k2 (t2.u)^2 + gamma (a.u - beta)^2] / 2
        // with t1, t2 orthonormal and orthogonal to a, so the principal
        // curvatures are exactly k1 and k2 whatever gamma (a curvature along the
        // gradient direction) is.  Diagonal Hessian entries would give mixtures.
        let a = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
        let t1 = [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0];
        let t2 = [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0];
        let (beta, k1, k2, gamma) = (2.2, 0.04, 0.18, 0.3);
        let dot3 = |p: &[f64; 3], u: &[f64]| p[0] * u[0] + p[1] * u[1] + p[2] * u[2];
        for scale in [1.0, 7.5, 0.01] {
            let mean = vec![1.0, -2.0, 0.5];
            let sigma = vec![0.8, 1.7, 2.2];
            let (m, s) = (mean.clone(), sigma.clone());
            let g = LimitState::new(mean, sigma, move |x| {
                let u: Vec<f64> = (0..3).map(|i| (x[i] - m[i]) / s[i]).collect();
                let along = dot3(&a, &u);
                scale
                    * (beta - along
                        + 0.5
                            * (k1 * dot3(&t1, &u).powi(2)
                                + k2 * dot3(&t2, &u).powi(2)
                                + gamma * (along - beta).powi(2)))
            });
            let u_star = vec![beta * a[0], beta * a[1], beta * a[2]];
            let sorm = SormAnalysis::solve(&g, &form_at(u_star));
            assert_eq!(sorm.curvatures.len(), 2);
            assert!(
                (sorm.curvatures[0] - k1).abs() < 5e-6 && (sorm.curvatures[1] - k2).abs() < 5e-6,
                "scale {scale}: {:?}",
                sorm.curvatures
            );
        }
    }
    #[test]
    fn test_sorm_curvature_sign_convention_matches_failure_domain_convexity() {
        let beta = 2.0;
        let build = |kappa: f64| {
            LimitState::new(vec![0.0, 0.0], vec![1.0, 1.0], move |x| {
                beta - x[1] + 0.5 * kappa * x[0] * x[0]
            })
        };
        let convex = SormAnalysis::solve(&build(0.2), &form_at(vec![0.0, beta]));
        let concave = SormAnalysis::solve(&build(-0.2), &form_at(vec![0.0, beta]));
        assert!(convex.curvatures[0] > 0.0 && convex.correction_factor < 1.0);
        assert!(concave.curvatures[0] < 0.0 && concave.correction_factor > 1.0);
        assert!(convex.pf_breitung < phi(-beta));
        assert!(concave.pf_breitung > phi(-beta));
    }
    /// Composite Simpson rule with `n` (even) intervals.
    fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
        let h = (b - a) / n as f64;
        let mut sum = f(a) + f(b);
        for i in 1..n {
            let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
            sum += weight * f(a + i as f64 * h);
        }
        sum * h / 3.0
    }
    /// Exact failure probability of `u_n >= beta + kappa u_1^2 / 2`.
    fn exact_pf_one_curvature(beta: f64, kappa: f64) -> f64 {
        simpson(
            |u| phi_pdf(u) * phi(-(beta + 0.5 * kappa * u * u)),
            -10.0,
            10.0,
            4000,
        )
    }
    #[test]
    fn test_tvedt_is_closer_to_the_exact_probability_than_breitung() {
        // One curvature: Pf = integral of phi(u) Phi(-beta - kappa u^2 / 2) du.
        for &beta in &[2.0, 2.5] {
            for &kappa in &[0.1, 0.2, -0.1] {
                let exact = exact_pf_one_curvature(beta, kappa);
                let breitung = phi(-beta) / (1.0 + beta * kappa).sqrt();
                let tvedt = SormAnalysis::tvedt_pf(beta, &[kappa]);
                let err_b = (breitung / exact - 1.0).abs();
                let err_t = (tvedt / exact - 1.0).abs();
                assert!(
                    err_b > 5e-3,
                    "beta {beta}, kappa {kappa}: Breitung error {err_b}"
                );
                assert!(
                    err_t < 3e-3,
                    "beta {beta}, kappa {kappa}: Tvedt error {err_t}"
                );
                assert!(
                    err_t < 0.25 * err_b,
                    "beta {beta}, kappa {kappa}: {err_t} vs {err_b}"
                );
            }
        }
    }
    #[test]
    fn test_tvedt_two_equal_curvatures_matches_chi_square_integral() {
        // u_1^2 + u_2^2 = s is chi-square(2): density exp(-s / 2) / 2.
        let (beta, kappa) = (2.0, 0.1);
        let exact = simpson(
            |s| 0.5 * (-0.5 * s).exp() * phi(-(beta + 0.5 * kappa * s)),
            0.0,
            120.0,
            24000,
        );
        let breitung = phi(-beta) / (1.0 + beta * kappa);
        let tvedt = SormAnalysis::tvedt_pf(beta, &[kappa, kappa]);
        assert!((breitung / exact - 1.0).abs() > 2e-2);
        assert!(
            (tvedt / exact - 1.0).abs() < 2e-3,
            "tvedt {tvedt}, exact {exact}"
        );
    }
    /// Exact failure probability of `u_n >= beta + (k1 u_1^2 + k2 u_2^2) / 2`
    /// (nested composite Simpson rule over `[-9, 9]^2`).
    fn exact_pf_two_curvatures(beta: f64, k1: f64, k2: f64) -> f64 {
        simpson(
            |u1| {
                phi_pdf(u1)
                    * simpson(
                        |u2| phi_pdf(u2) * phi(-(beta + 0.5 * (k1 * u1 * u1 + k2 * u2 * u2))),
                        -9.0,
                        9.0,
                        180,
                    )
            },
            -9.0,
            9.0,
            180,
        )
    }
    #[test]
    fn test_tvedt_distinct_and_mixed_sign_curvatures_match_two_dimensional_integral() {
        // The factor Re prod_j (1 + (beta + i) kappa_j)^(-1/2) couples the
        // curvatures through the phase of each term; equal curvatures alone
        // cannot expose an error in how distinct (or opposite) ones combine.
        for &(beta, k1, k2) in &[
            (2.5, 0.2, 0.05),
            (2.5, 0.12, -0.06),
            (2.0, 0.1, -0.05),
            (3.0, 0.15, -0.1),
        ] {
            let exact = exact_pf_two_curvatures(beta, k1, k2);
            let breitung = phi(-beta) / ((1.0 + beta * k1) * (1.0 + beta * k2)).sqrt();
            let tvedt = SormAnalysis::tvedt_pf(beta, &[k1, k2]);
            let err_breitung = (breitung / exact - 1.0).abs();
            let err_tvedt = (tvedt / exact - 1.0).abs();
            assert!(
                err_tvedt < 2e-3,
                "beta {beta}, kappa ({k1}, {k2}): Tvedt error {err_tvedt}"
            );
            assert!(
                err_tvedt < 0.5 * err_breitung,
                "beta {beta}, kappa ({k1}, {k2}): Tvedt {err_tvedt} vs Breitung {err_breitung}"
            );
            // The result does not depend on the order of the curvatures.
            let swapped = SormAnalysis::tvedt_pf(beta, &[k2, k1]);
            assert!(
                (tvedt - swapped).abs() <= 1e-13 * tvedt,
                "{tvedt} vs {swapped}"
            );
        }
    }
    #[test]
    fn test_tvedt_reduces_to_breitung_at_zero_and_is_nan_outside_validity() {
        let beta = 2.7;
        assert!((SormAnalysis::tvedt_pf(beta, &[0.0, 0.0]) - phi(-beta)).abs() < 1e-14);
        // 1 + beta kappa <= 0: outside the asymptotic range.
        assert!(SormAnalysis::tvedt_pf(2.0, &[-0.5]).is_nan());
        assert!(
            SormAnalysis::tvedt_pf(2.0, &[-0.4]).is_nan(),
            "1 + 3 kappa < 0"
        );
        let pf = SormAnalysis::tvedt_pf(3.0, &[0.3, -0.1, 0.05]);
        assert!(pf.is_finite() && (0.0..=1.0).contains(&pf));
    }
    #[test]
    fn test_sorm_end_to_end_form_then_sorm_matches_exact_probability() {
        // x-space limit state with non-trivial mean/sigma whose standard-normal
        // image is G(u) = beta - u_2 + kappa u_1^2 / 2.
        let (beta, kappa) = (2.5, 0.1);
        let mean = vec![5.0, -1.0];
        let sigma = vec![0.5, 2.0];
        let (m, s) = (mean.clone(), sigma.clone());
        let g = LimitState::new(mean, sigma, move |x| {
            let u1 = (x[0] - m[0]) / s[0];
            let u2 = (x[1] - m[1]) / s[1];
            beta - u2 + 0.5 * kappa * u1 * u1
        });
        let form = FormAnalysis::solve(&g);
        assert!((form.beta - beta).abs() < 1e-5, "beta = {}", form.beta);
        let sorm = SormAnalysis::solve(&g, &form);
        assert_eq!(sorm.curvatures.len(), 1);
        assert!(
            (sorm.curvatures[0] - kappa).abs() < 1e-5,
            "{:?}",
            sorm.curvatures
        );
        let exact = exact_pf_one_curvature(beta, kappa);
        assert!((sorm.pf_breitung / exact - 1.0).abs() < 3e-2);
        let tvedt = SormAnalysis::tvedt_pf(sorm.beta, &sorm.curvatures);
        assert!(
            (tvedt / exact - 1.0).abs() < 3e-3,
            "tvedt {tvedt}, exact {exact}"
        );
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
        assert!((0.0..=1.0).contains(&pf));
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
        assert!((0.0..=1.0).contains(&pf), "pf={pf}");
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
        assert!((0.0..=1.0).contains(&lb));
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
            assert!((0.0..=1.0).contains(s), "s={s}");
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
        assert!((0.0..=1.0).contains(&pf));
    }
    #[test]
    fn test_normality_statistic_returns_finite() {
        let samples: Vec<f64> = (0..50).map(|i| i as f64 * 0.1).collect();
        let stat = normality_statistic(&samples);
        assert!(stat.is_finite());
    }
}
