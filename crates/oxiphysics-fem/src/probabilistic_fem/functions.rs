//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains)]
use std::f64::consts::PI;

/// Standard normal CDF Φ(x) via Abramowitz & Stegun 26.2.17.
///
/// Maximum absolute error ≈ 7.5 × 10⁻⁸.
pub fn phi_cdf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let pdf_v = (-0.5 * x * x).exp() / (2.0 * PI).sqrt();
    let cdf = 1.0 - pdf_v * poly;
    if x >= 0.0 { cdf } else { 1.0 - cdf }
}
/// Standard normal PDF φ(x).
pub fn phi_pdf(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * PI).sqrt()
}
/// Inverse normal CDF via bisection (100 iterations).
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
        if phi_cdf(mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
/// Box-Muller transform: standard normal from two U(0,1) samples.
pub(super) fn box_muller(u1: f64, u2: f64) -> f64 {
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}
/// Mean of a slice.
pub(super) fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.iter().sum::<f64>() / v.len() as f64
}
/// Sample variance of a slice.
pub(super) fn variance(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let m = mean(v);
    v.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64
}
/// Coefficient of variation (std/mean).
pub(super) fn cov_stat(v: &[f64]) -> f64 {
    let m = mean(v);
    if m.abs() < 1e-14 {
        return f64::INFINITY;
    }
    variance(v).sqrt() / m.abs()
}
/// Dot product of two equal-length slices.
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}
/// L2 norm of a slice.
pub(super) fn norm2(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}
/// Exponential covariance C(r) = σ² exp(-r / ℓ).
pub fn cov_exp(r: f64, sigma2: f64, length: f64) -> f64 {
    if length < 1e-300 {
        return 0.0;
    }
    sigma2 * (-r / length).exp()
}
/// Squared-exponential covariance C(r) = σ² exp(-r²/(2ℓ²)).
pub fn cov_sqexp(r: f64, sigma2: f64, length: f64) -> f64 {
    if length < 1e-300 {
        return 0.0;
    }
    sigma2 * (-r * r / (2.0 * length * length)).exp()
}
/// Matérn ν=5/2 covariance kernel.
pub fn cov_matern52(r: f64, sigma2: f64, length: f64) -> f64 {
    if length < 1e-300 {
        return 0.0;
    }
    let s = 5.0_f64.sqrt() * r / length;
    sigma2 * (1.0 + s + s * s / 3.0) * (-s).exp()
}
/// Power-law covariance C(r) = σ² (1 + r/ℓ)^(-α).
pub fn cov_power(r: f64, sigma2: f64, length: f64, alpha: f64) -> f64 {
    if length < 1e-300 {
        return 0.0;
    }
    sigma2 * (1.0 + r / length).powf(-alpha)
}
/// Cholesky factorization L such that A = L Lᵀ.
///
/// Returns None if matrix is not positive definite.
pub(super) fn cholesky(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let sum: f64 = (0..j).map(|k| l[i][k] * l[j][k]).sum();
            if i == j {
                let diag = a[i][i] - sum;
                if diag <= 0.0 {
                    return None;
                }
                l[i][j] = diag.sqrt();
            } else {
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }
    Some(l)
}
/// Solve L x = b (lower triangular forward substitution).
pub(super) fn forward_sub(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0f64; n];
    for i in 0..n {
        let sum: f64 = (0..i).map(|j| l[i][j] * x[j]).sum();
        x[i] = (b[i] - sum) / l[i][i].max(1e-300);
    }
    x
}
/// Solve Lᵀ x = b (upper triangular back substitution).
pub(super) fn back_sub(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let sum: f64 = (i + 1..n).map(|j| l[j][i] * x[j]).sum();
        x[i] = (b[i] - sum) / l[i][i].max(1e-300);
    }
    x
}
/// Solve A x = b using Cholesky (returns zeros if factorization fails).
pub(super) fn solve_chol(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    match cholesky(a) {
        Some(l) => {
            let y = forward_sub(&l, b);
            back_sub(&l, &y)
        }
        None => vec![0.0; b.len()],
    }
}
/// Power iteration to find dominant eigenvalue and eigenvector.
///
/// Returns (eigenvalue, eigenvector).
pub(super) fn power_iteration(a: &[Vec<f64>], max_iter: usize) -> (f64, Vec<f64>) {
    let n = a.len();
    if n == 0 {
        return (0.0, vec![]);
    }
    let mut v = vec![1.0f64; n];
    let nv = norm2(&v).max(1e-300);
    for x in &mut v {
        *x /= nv;
    }
    let mut eigenval = 0.0;
    for _ in 0..max_iter {
        let mut av = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                av[i] += a[i][j] * v[j];
            }
        }
        eigenval = dot(&v, &av);
        let nrm = norm2(&av).max(1e-300);
        for i in 0..n {
            v[i] = av[i] / nrm;
        }
    }
    (eigenval, v)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::probabilistic_fem::types::FailureMode;
    use crate::probabilistic_fem::*;
    #[test]
    fn test_phi_cdf_symmetry() {
        assert!((phi_cdf(0.0) - 0.5).abs() < 1e-6);
        assert!((phi_cdf(1.645) - 0.95).abs() < 1e-3);
    }
    #[test]
    fn test_phi_pdf_positive() {
        assert!(phi_pdf(0.0) > 0.0);
        assert!((phi_pdf(0.0) - 1.0 / (2.0 * PI).sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_phi_inv_roundtrip() {
        let p = 0.975;
        let x = phi_inv(p);
        let p_back = phi_cdf(x);
        assert!((p_back - p).abs() < 1e-5);
    }
    #[test]
    fn test_mean_variance() {
        let v = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((mean(&v) - 3.0).abs() < 1e-10);
        assert!((variance(&v) - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_cov_stat() {
        let v = vec![10.0, 12.0, 8.0, 11.0, 9.0];
        let c = cov_stat(&v);
        assert!(c > 0.0 && c < 1.0);
    }
    #[test]
    fn test_cov_exp_origin() {
        let c = cov_exp(0.0, 1.0, 1.0);
        assert!((c - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_cov_sqexp_decreasing() {
        let c0 = cov_sqexp(0.0, 1.0, 1.0);
        let c1 = cov_sqexp(1.0, 1.0, 1.0);
        assert!(c0 > c1);
    }
    #[test]
    fn test_cov_matern52_positive() {
        for r in [0.0, 0.5, 1.0, 2.0] {
            assert!(cov_matern52(r, 1.0, 1.0) > 0.0);
        }
    }
    #[test]
    fn test_cov_power_decreasing() {
        let c0 = cov_power(0.5, 1.0, 1.0, 2.0);
        let c1 = cov_power(2.0, 1.0, 1.0, 2.0);
        assert!(c0 > c1);
    }
    #[test]
    fn test_cholesky_2x2() {
        let a = vec![vec![4.0, 2.0], vec![2.0, 3.0]];
        let l = cholesky(&a).expect("should factor");
        let l00 = l[0][0];
        let l10 = l[1][0];
        let l11 = l[1][1];
        assert!((l00 * l00 - 4.0).abs() < 1e-10);
        assert!((l10 * l00 - 2.0).abs() < 1e-10);
        assert!((l10 * l10 + l11 * l11 - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_solve_chol_identity() {
        let a = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        let b = vec![4.0, 9.0];
        let x = solve_chol(&a, &b);
        assert!((x[0] - 2.0).abs() < 1e-10);
        assert!((x[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_random_field_construction() {
        let rf = RandomField::new(10, 1.0, 3, 0.1, 0.3, 0);
        assert_eq!(rf.mesh.len(), 10);
        assert!(!rf.eigenvalues.is_empty());
    }
    #[test]
    fn test_random_field_eigenvalues_positive() {
        let rf = RandomField::new(8, 1.0, 3, 0.1, 0.3, 0);
        assert!(rf.eigenvalues.iter().all(|&l| l > 0.0));
    }
    #[test]
    fn test_random_field_energy_ratio() {
        let rf = RandomField::new(8, 1.0, 4, 0.1, 0.2, 1);
        let er = rf.energy_ratio();
        assert!(er >= 0.0 && er <= 1.01);
    }
    #[test]
    fn test_random_field_sample_length() {
        let rf = RandomField::new(8, 1.0, 3, 0.1, 0.3, 0);
        let xi = vec![1.0, -0.5, 0.3];
        let s = rf.sample(&xi);
        assert_eq!(s.len(), 8);
    }
    #[test]
    fn test_random_field_sample_finite() {
        let rf = RandomField::new(8, 1.0, 3, 0.05, 0.3, 2);
        let xi = vec![0.0; 3];
        let s = rf.sample(&xi);
        assert!(s.iter().all(|&v| v.is_finite()));
    }
    #[test]
    fn test_sfem_construction() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 2e11, 0.1, 2);
        assert_eq!(sfem.n_elem, 4);
    }
    #[test]
    fn test_sfem_mean_solution_positive() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.1, 2);
        let u = sfem.solve_mean();
        let u_tip = *u.last().unwrap_or(&0.0);
        assert!(
            u_tip > 0.0,
            "tip displacement should be positive under tip load"
        );
    }
    #[test]
    fn test_sfem_analytical_match() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.0, 1);
        let u_tip = *sfem.solve_mean().last().unwrap_or(&0.0);
        let u_anal = sfem.analytical_tip_displacement();
        assert!(
            (u_tip - u_anal).abs() / u_anal < 0.01,
            "numerical tip {u_tip} vs analytical {u_anal}"
        );
    }
    #[test]
    fn test_sfem_sample_nonzero() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let xi = vec![1.0, -0.5];
        let u = sfem.solve_sample(&xi);
        assert!(!u.is_empty());
        assert!(u.iter().all(|&v| v.is_finite()));
    }
    #[test]
    fn test_mc_fem_run() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let mut mc = MonteCarloFem::new(sfem, 50);
        mc.run();
        assert_eq!(mc.samples.len(), 50);
    }
    #[test]
    fn test_mc_fem_mean_positive() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let mut mc = MonteCarloFem::new(sfem, 50);
        mc.run();
        assert!(mc.mean_tip() > 0.0);
    }
    #[test]
    fn test_mc_fem_std_nonneg() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let mut mc = MonteCarloFem::new(sfem, 50);
        mc.run();
        assert!(mc.std_tip() >= 0.0);
    }
    #[test]
    fn test_mc_fem_failure_prob() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let mut mc = MonteCarloFem::new(sfem, 50);
        mc.run();
        let pf = mc.failure_probability(1e10);
        assert!(pf >= 0.0 && pf <= 1.0);
    }
    #[test]
    fn test_mc_fem_percentile_95() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let mut mc = MonteCarloFem::new(sfem, 50);
        mc.run();
        let p95 = mc.percentile_95();
        assert!(p95 >= mc.mean_tip());
    }
    #[test]
    fn test_perturbation_fem_zeroth_order() {
        let pfem = PerturbationFem::new(4, 1.0, 1.0, 100.0, 1e6, 1e10, 0.3);
        let u0_tip = pfem.first_order_mean_tip();
        assert!(u0_tip > 0.0);
    }
    #[test]
    fn test_perturbation_fem_variance_nonneg() {
        let pfem = PerturbationFem::new(4, 1.0, 1.0, 100.0, 1e6, 1e10, 0.3);
        let var = pfem.first_order_variance_tip();
        assert!(var >= 0.0 && var.is_finite());
    }
    #[test]
    fn test_perturbation_fem_second_order_correction() {
        let pfem = PerturbationFem::new(4, 1.0, 1.0, 100.0, 1e6, 1e10, 0.3);
        let corr = pfem.second_order_mean_correction();
        assert!(corr.is_finite());
    }
    #[test]
    fn test_hermite_polynomial() {
        assert!((SpectralStochasticFem::hermite(0, 1.5) - 1.0).abs() < 1e-12);
        assert!((SpectralStochasticFem::hermite(1, 2.0) - 2.0).abs() < 1e-12);
        let h2 = SpectralStochasticFem::hermite(2, 2.0);
        assert!((h2 - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_pce_norm_sq() {
        let alpha = vec![2, 1];
        let ns = SpectralStochasticFem::pce_norm_sq(&alpha);
        assert!((ns - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_spectral_sfem_multi_indices() {
        let ssfem = SpectralStochasticFem::new(2, 2, 4, 1.0, 1.0, 100.0, 1e6);
        let indices = ssfem.multi_indices();
        assert!(!indices.is_empty());
        assert!(indices.iter().all(|a| a.iter().sum::<usize>() <= 2));
    }
    #[test]
    fn test_spectral_sfem_pce_mean() {
        let mut ssfem = SpectralStochasticFem::new(1, 1, 4, 1.0, 1.0, 100.0, 1e6);
        ssfem.compute_pce_coefficients(0.05);
        let m = ssfem.pce_mean();
        assert!(m.is_finite() && m > 0.0);
    }
    #[test]
    fn test_spectral_sfem_pce_variance_nonneg() {
        let mut ssfem = SpectralStochasticFem::new(1, 1, 4, 1.0, 1.0, 100.0, 1e6);
        ssfem.compute_pce_coefficients(0.05);
        let var = ssfem.pce_variance();
        assert!(var >= 0.0 && var.is_finite());
    }
    #[test]
    fn test_ddm_sensitivity_da_length() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let ddm = DdmSensitivity::new(sfem);
        let ds = ddm.sensitivity_da();
        assert_eq!(ds.len(), 5);
    }
    #[test]
    fn test_ddm_sensitivity_de_sign() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let ddm = DdmSensitivity::new(sfem);
        let ds = ddm.sensitivity_de();
        let tip_sens = *ds.last().unwrap_or(&0.0);
        assert!(
            tip_sens < 0.0,
            "increasing E should decrease tip displacement (negative sens)"
        );
    }
    #[test]
    fn test_ddm_adjoint_finite() {
        let sfem = StochasticFem::new(4, 1.0, 1.0, 100.0, 1e6, 0.05, 2);
        let ddm = DdmSensitivity::new(sfem);
        let adj = ddm.adjoint_sensitivity_e(4);
        assert!(adj.is_finite());
    }
    #[test]
    fn test_fuzzy_fem_construction() {
        let ffem = FuzzyFem::new(4, 1.0, 1.0, 100.0, 0.8e6, 1.0e6, 1.2e6, 5);
        assert_eq!(ffem.e_fuzzy_cuts.len(), 5);
    }
    #[test]
    fn test_fuzzy_fem_alpha_cuts_ordered() {
        let ffem = FuzzyFem::new(4, 1.0, 1.0, 100.0, 0.8e6, 1.0e6, 1.2e6, 5);
        for cut in &ffem.e_fuzzy_cuts {
            assert!(cut.lower <= cut.upper);
        }
    }
    #[test]
    fn test_fuzzy_fem_response() {
        let mut ffem = FuzzyFem::new(4, 1.0, 1.0, 100.0, 0.8e6, 1.0e6, 1.2e6, 5);
        ffem.compute_fuzzy_response();
        assert_eq!(ffem.u_tip_cuts.len(), 5);
        for cut in &ffem.u_tip_cuts {
            assert!(cut.lower <= cut.upper);
            assert!(cut.lower >= 0.0);
        }
    }
    #[test]
    fn test_interval_arithmetic() {
        let a = Interval::new(1.0, 2.0);
        let b = Interval::new(3.0, 4.0);
        let c = a.mul(b);
        assert!((c.lo - 3.0).abs() < 1e-12);
        assert!((c.hi - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_fem_vertex_method() {
        let mut ifem = IntervalFem::new(4, 1.0, 1.0, 100.0, 0.8e6, 1.2e6);
        ifem.vertex_method();
        let iv = ifem.u_tip_interval;
        assert!(iv.lo <= iv.hi);
        assert!(iv.lo > 0.0);
    }
    #[test]
    fn test_interval_fem_midpoint_in_interval() {
        let mut ifem = IntervalFem::new(4, 1.0, 1.0, 100.0, 0.8e6, 1.2e6);
        ifem.vertex_method();
        let mid = ifem.midpoint_solution();
        let iv = ifem.u_tip_interval;
        assert!(
            mid >= iv.lo - 1e-10 && mid <= iv.hi + 1e-10,
            "midpoint {mid} not in [{}, {}]",
            iv.lo,
            iv.hi
        );
    }
    #[test]
    fn test_failure_mode_construction() {
        let m = FailureMode::new("yielding", 3.0, vec![0.8, 0.6]);
        assert!((m.pf - phi_cdf(-3.0)).abs() < 1e-10);
    }
    #[test]
    fn test_failure_mode_importance() {
        let m = FailureMode::new("buckling", 2.5, vec![0.6, 0.8]);
        assert!((m.importance(0) - 0.36).abs() < 1e-10);
        assert!((m.importance(1) - 0.64).abs() < 1e-10);
    }
    #[test]
    fn test_failure_mode_analysis_series_upper() {
        let m1 = FailureMode::new("m1", 3.0, vec![1.0, 0.0]);
        let m2 = FailureMode::new("m2", 2.5, vec![0.0, 1.0]);
        let fma = FailureModeAnalysis::new(vec![m1, m2]);
        let p_sys = fma.series_pf_upper();
        assert!(p_sys >= 0.0 && p_sys <= 1.0);
    }
    #[test]
    fn test_failure_mode_analysis_parallel_pf() {
        let m1 = FailureMode::new("m1", 3.0, vec![1.0, 0.0]);
        let m2 = FailureMode::new("m2", 3.0, vec![0.0, 1.0]);
        let fma = FailureModeAnalysis::new(vec![m1, m2]);
        let p_par = fma.parallel_pf();
        let p_single = phi_cdf(-3.0);
        assert!((p_par - p_single * p_single).abs() < 1e-14);
    }
    #[test]
    fn test_failure_mode_most_likely() {
        let m1 = FailureMode::new("m1", 3.5, vec![1.0]);
        let m2 = FailureMode::new("m2", 2.0, vec![1.0]);
        let fma = FailureModeAnalysis::new(vec![m1, m2]);
        let ml = fma.most_likely_mode().unwrap();
        assert_eq!(ml.name, "m2");
    }
    #[test]
    fn test_failure_mode_global_sensitivity() {
        let m1 = FailureMode::new("m1", 3.0, vec![1.0, 0.0]);
        let m2 = FailureMode::new("m2", 3.0, vec![0.0, 1.0]);
        let fma = FailureModeAnalysis::new(vec![m1, m2]);
        let gs = fma.global_sensitivity();
        let total: f64 = gs.iter().sum();
        assert!((total - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_failure_mode_system_beta() {
        let m1 = FailureMode::new("m1", 3.0, vec![1.0]);
        let fma = FailureModeAnalysis::new(vec![m1]);
        let beta_sys = fma.system_beta_series();
        assert!((beta_sys - 3.0).abs() < 1e-3);
    }
    #[test]
    fn test_rbdo_solver_runs() {
        let sfem = StochasticFem::new(4, 0.5, 1.0, 100.0, 1e6, 0.1, 2);
        let dv = vec![DesignVar {
            name: "area".to_string(),
            value: 0.5,
            lb: 0.1,
            ub: 2.0,
        }];
        let mut rbdo = RbdoSolver::new(sfem, 2.0, 10, dv);
        rbdo.run(1e-3);
        assert!(!rbdo.beta_history.is_empty());
        assert!(rbdo.beta_history.iter().all(|&b| b.is_finite()));
    }
}
