//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains, clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::types::StochasticResponse;

/// Standard normal CDF  Φ(x) via Abramowitz & Stegun approximation.
///
/// Maximum absolute error ≈ 7.5 × 10⁻⁸.
pub fn phi_cdf(x: f64) -> f64 {
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
/// Inverse standard normal CDF (bisection, 100 iterations).
///
/// Returns the value β such that Φ(β) ≈ p.
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
/// Box-Muller transform: return one standard-normal sample from two uniform
/// samples `u1 ∈ (0,1]` and `u2 ∈ [0,1)`.
pub fn box_muller(u1: f64, u2: f64) -> f64 {
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}
/// Estimate (mean displacement, std displacement) via Monte Carlo sampling.
///
/// Models a 1-DOF spring: u = F / k, where k ~ N(k_mean, k_std²) (truncated
/// to k > 0).  Returns (mean_u, std_u) over `n_samples` samples.
///
/// # Arguments
/// * `n_samples` – number of Monte Carlo samples
/// * `k_mean`    – mean stiffness (N/m)
/// * `k_std`     – standard deviation of stiffness
/// * `load`      – applied force F (N)
/// * `n_dof`     – number of degrees of freedom (unused placeholder for extension)
pub fn monte_carlo_fem(
    n_samples: usize,
    k_mean: f64,
    k_std: f64,
    load: f64,
    _n_dof: usize,
) -> (f64, f64) {
    if n_samples == 0 || k_mean <= 0.0 {
        return (0.0, 0.0);
    }
    use rand::RngExt as _;
    let mut rng = rand::rng();
    let mut samples = Vec::with_capacity(n_samples);
    for _ in 0..n_samples {
        let u1: f64 = rng.random_range(f64::EPSILON..1.0_f64);
        let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
        let z = box_muller(u1, u2);
        let k = (k_mean + k_std * z).max(k_mean * 1e-6);
        samples.push(load / k);
    }
    let mean = samples.iter().sum::<f64>() / n_samples as f64;
    let var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n_samples as f64;
    (mean, var.sqrt())
}
/// Exponential covariance kernel.
///
/// C(x₁, x₂) = σ² · exp(−|x₁ − x₂| / l_corr)
///
/// # Arguments
/// * `x1`     – first spatial coordinate
/// * `x2`     – second spatial coordinate
/// * `sigma`  – standard deviation of the field
/// * `l_corr` – correlation length
pub fn covariance_exponential(x1: f64, x2: f64, sigma: f64, l_corr: f64) -> f64 {
    if l_corr <= 0.0 {
        return if (x1 - x2).abs() < f64::EPSILON {
            sigma * sigma
        } else {
            0.0
        };
    }
    sigma * sigma * (-(x1 - x2).abs() / l_corr).exp()
}
/// Gaussian (squared-exponential) covariance kernel.
///
/// C(x₁, x₂) = σ² · exp(−(x₁ − x₂)² / (2 l_corr²))
///
/// # Arguments
/// * `x1`     – first spatial coordinate
/// * `x2`     – second spatial coordinate
/// * `sigma`  – standard deviation of the field
/// * `l_corr` – correlation length
pub fn covariance_gaussian(x1: f64, x2: f64, sigma: f64, l_corr: f64) -> f64 {
    if l_corr <= 0.0 {
        return if (x1 - x2).abs() < f64::EPSILON {
            sigma * sigma
        } else {
            0.0
        };
    }
    let d = x1 - x2;
    sigma * sigma * (-(d * d) / (2.0 * l_corr * l_corr)).exp()
}
/// Evaluate a Hermite polynomial chaos expansion at `xi`.
///
/// f(ξ) = Σₙ coeffs\[n\] · He_n(ξ)
///
/// where He_n is the probabilist's Hermite polynomial (see [`hermite_polynomial`]).
///
/// # Arguments
/// * `coeffs`    – PCE coefficients c₀, c₁, …
/// * `xi`        – evaluation point (standard normal coordinate)
/// * `max_order` – maximum polynomial order to include
pub fn polynomial_chaos_expansion(coeffs: &[f64], xi: f64, max_order: usize) -> f64 {
    let n = coeffs.len().min(max_order + 1);
    (0..n).map(|i| coeffs[i] * hermite_polynomial(i, xi)).sum()
}
/// Probabilist's Hermite polynomial He_n(x).
///
/// Defined by the three-term recurrence:
/// - He₀(x) = 1
/// - He₁(x) = x
/// - He_{n+1}(x) = x · He_n(x) − n · He_{n-1}(x)
///
/// # Arguments
/// * `n` – polynomial order
/// * `x` – evaluation point
pub fn hermite_polynomial(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut h_prev = 1.0_f64;
            let mut h_curr = x;
            for k in 1..n {
                let h_next = x * h_curr - (k as f64) * h_prev;
                h_prev = h_curr;
                h_curr = h_next;
            }
            h_curr
        }
    }
}
/// PCE mean (zeroth-order coefficient).
///
/// E\[f\] = coeffs\[0\] since E\[He_n(ξ)\] = 0 for n ≥ 1.
///
/// # Arguments
/// * `coeffs` – PCE coefficients
pub fn pce_mean(coeffs: &[f64]) -> f64 {
    if coeffs.is_empty() { 0.0 } else { coeffs[0] }
}
/// PCE variance: Var\[f\] = Σ_{n≥1} coeffs\[n\]² · n!
///
/// Uses the orthogonality of Hermite polynomials with respect to the
/// standard Gaussian measure.
///
/// # Arguments
/// * `coeffs` – PCE coefficients
pub fn pce_variance(coeffs: &[f64]) -> f64 {
    if coeffs.len() < 2 {
        return 0.0;
    }
    coeffs[1..].iter().map(|c| c * c).sum()
}
/// Gauss-Hermite quadrature nodes and weights for `n`-point rule.
///
/// Returns `(nodes, weights)` for approximating ∫ f(x) exp(-x²) dx.
/// Supports n = 1, 2, 3, 4, 5; returns empty vecs for other values.
pub fn gauss_hermite_quadrature(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![0.0], vec![PI.sqrt()]),
        2 => (
            vec![
                -std::f64::consts::FRAC_1_SQRT_2,
                std::f64::consts::FRAC_1_SQRT_2,
            ],
            vec![0.886_226_925_452_758, 0.886_226_925_452_758],
        ),
        3 => (
            vec![-1.224_744_871_391_589, 0.0, 1.224_744_871_391_589],
            vec![
                0.295_408_975_150_920,
                1.181_635_900_603_681,
                0.295_408_975_150_920,
            ],
        ),
        4 => (
            vec![
                -1.650_680_123_885_785,
                -0.524_647_623_275_290,
                0.524_647_623_275_290,
                1.650_680_123_885_785,
            ],
            vec![
                0.081_312_835_447_245,
                0.804_914_090_005_513,
                0.804_914_090_005_513,
                0.081_312_835_447_245,
            ],
        ),
        5 => (
            vec![
                -2.020_182_470_003_108,
                -0.958_572_464_613_819,
                0.0,
                0.958_572_464_613_819,
                2.020_182_470_003_108,
            ],
            vec![
                0.019_953_242_059_046,
                0.394_424_314_739_089,
                0.945_308_720_482_942,
                0.394_424_314_739_089,
                0.019_953_242_059_046,
            ],
        ),
        _ => (vec![], vec![]),
    }
}
/// Estimate P(u > threshold) from Monte Carlo samples.
///
/// # Arguments
/// * `response`  – [`StochasticResponse`] holding the samples
/// * `threshold` – failure threshold
pub fn failure_probability_monte_carlo(response: &StochasticResponse, threshold: f64) -> f64 {
    if response.samples.is_empty() {
        return 0.0;
    }
    let n_fail = response.samples.iter().filter(|&&s| s > threshold).count();
    n_fail as f64 / response.samples.len() as f64
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::stochastic_fem::*;
    #[test]
    fn test_sfem_new() {
        let sf = StochasticFEM::new(10, 3, 1000.0, 100.0);
        assert_eq!(sf.n_nodes, 10);
        assert_eq!(sf.n_random, 3);
        assert!((sf.mean_k - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_sfem_cov() {
        let sf = StochasticFEM::new(5, 2, 200.0, 20.0);
        assert!((sf.cov() - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_sfem_cov_zero_mean() {
        let sf = StochasticFEM::new(5, 2, 0.0, 10.0);
        assert_eq!(sf.cov(), 0.0);
    }
    #[test]
    fn test_hermite_he0_is_one() {
        assert!((hermite_polynomial(0, 0.0) - 1.0).abs() < 1e-14);
        assert!((hermite_polynomial(0, 5.0) - 1.0).abs() < 1e-14);
        assert!((hermite_polynomial(0, -3.0) - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_hermite_he1_is_x() {
        assert!((hermite_polynomial(1, 0.0) - 0.0).abs() < 1e-14);
        assert!((hermite_polynomial(1, 3.0) - 3.0).abs() < 1e-14);
        assert!((hermite_polynomial(1, -2.0) - (-2.0)).abs() < 1e-14);
    }
    #[test]
    fn test_hermite_he2_is_x2_minus_1() {
        for x in [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0] {
            let expected = x * x - 1.0;
            let got = hermite_polynomial(2, x);
            assert!(
                (got - expected).abs() < 1e-12,
                "He2({x}) = {got}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_hermite_he3() {
        for x in [-2.0, 0.0, 1.0, 2.0] {
            let expected = x * x * x - 3.0 * x;
            let got = hermite_polynomial(3, x);
            assert!(
                (got - expected).abs() < 1e-10,
                "He3({x}) = {got}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_hermite_he4() {
        let x = 2.0_f64;
        let expected = x.powi(4) - 6.0 * x * x + 3.0;
        let got = hermite_polynomial(4, x);
        assert!((got - expected).abs() < 1e-10);
    }
    #[test]
    fn test_pce_mean_is_first_coeff() {
        let coeffs = [3.0, 1.0, 0.5];
        assert!((pce_mean(&coeffs) - 3.0).abs() < 1e-14);
    }
    #[test]
    fn test_pce_mean_empty() {
        assert_eq!(pce_mean(&[]), 0.0);
    }
    #[test]
    fn test_pce_variance_sum_of_squares() {
        let coeffs = [0.0, 2.0, 3.0];
        let expected = 4.0 + 9.0;
        assert!((pce_variance(&coeffs) - expected).abs() < 1e-14);
    }
    #[test]
    fn test_pce_variance_single_coeff() {
        let coeffs = [5.0];
        assert_eq!(pce_variance(&coeffs), 0.0);
    }
    #[test]
    fn test_pce_variance_empty() {
        assert_eq!(pce_variance(&[]), 0.0);
    }
    #[test]
    fn test_pce_constant_field() {
        let coeffs = [7.0, 0.0, 0.0];
        for xi in [-2.0, 0.0, 1.0, 3.0] {
            let val = polynomial_chaos_expansion(&coeffs, xi, 2);
            assert!(
                (val - 7.0).abs() < 1e-12,
                "constant PCE failed at xi={xi}: {val}"
            );
        }
    }
    #[test]
    fn test_pce_linear_field() {
        let coeffs = [0.0, 1.0];
        let xi = 2.5;
        let val = polynomial_chaos_expansion(&coeffs, xi, 1);
        assert!((val - xi).abs() < 1e-12);
    }
    #[test]
    fn test_pce_max_order_truncates() {
        let coeffs = [3.0, 1.0, 2.0];
        let val = polynomial_chaos_expansion(&coeffs, 5.0, 0);
        assert!((val - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_cov_exp_same_point_equals_variance() {
        let c = covariance_exponential(1.5, 1.5, 2.0, 0.5);
        assert!((c - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_cov_exp_positive() {
        let c = covariance_exponential(0.0, 1.0, 1.0, 1.0);
        assert!(c > 0.0);
    }
    #[test]
    fn test_cov_exp_decreases_with_distance() {
        let c1 = covariance_exponential(0.0, 0.1, 1.0, 1.0);
        let c2 = covariance_exponential(0.0, 1.0, 1.0, 1.0);
        assert!(c1 > c2);
    }
    #[test]
    fn test_cov_exp_correlation_leq_1() {
        let sigma = 2.0;
        let c_max = covariance_exponential(0.0, 0.0, sigma, 1.0);
        let c_off = covariance_exponential(0.0, 0.5, sigma, 1.0);
        assert!(c_off / c_max <= 1.0 + 1e-12);
    }
    #[test]
    fn test_cov_exp_symmetric() {
        let c1 = covariance_exponential(1.0, 3.0, 1.0, 2.0);
        let c2 = covariance_exponential(3.0, 1.0, 1.0, 2.0);
        assert!((c1 - c2).abs() < 1e-14);
    }
    #[test]
    fn test_cov_gauss_same_point_equals_variance() {
        let c = covariance_gaussian(0.0, 0.0, 3.0, 1.0);
        assert!((c - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_cov_gauss_positive() {
        let c = covariance_gaussian(0.0, 2.0, 1.0, 5.0);
        assert!(c > 0.0);
    }
    #[test]
    fn test_cov_gauss_symmetric() {
        let c1 = covariance_gaussian(1.0, 4.0, 2.0, 3.0);
        let c2 = covariance_gaussian(4.0, 1.0, 2.0, 3.0);
        assert!((c1 - c2).abs() < 1e-14);
    }
    #[test]
    fn test_kle_new_empty() {
        let kle = KarhunenLoeveExpansion::new(3);
        assert_eq!(kle.n_terms, 3);
        assert_eq!(kle.eigenvalues.len(), 0);
    }
    #[test]
    fn test_kle_add_mode_and_total_variance() {
        let mut kle = KarhunenLoeveExpansion::new(2);
        kle.add_mode(4.0, vec![1.0, 0.0]);
        kle.add_mode(1.0, vec![0.0, 1.0]);
        assert!((kle.total_variance() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_kle_sample_two_modes() {
        let mut kle = KarhunenLoeveExpansion::new(2);
        kle.add_mode(1.0, vec![1.0, 0.0, 0.0]);
        kle.add_mode(1.0, vec![0.0, 1.0, 0.0]);
        let result = kle.sample(&[1.0, 1.0]);
        assert!((result[0] - 1.0).abs() < 1e-12);
        assert!((result[1] - 1.0).abs() < 1e-12);
        assert!(result[2].abs() < 1e-12);
    }
    #[test]
    fn test_kle_sample_zero_xi_gives_zero() {
        let mut kle = KarhunenLoeveExpansion::new(1);
        kle.add_mode(5.0, vec![1.0, 2.0]);
        let result = kle.sample(&[0.0]);
        for v in result {
            assert!(v.abs() < 1e-12);
        }
    }
    #[test]
    fn test_kle_energy_ratio() {
        let mut kle = KarhunenLoeveExpansion::new(2);
        kle.add_mode(3.0, vec![1.0]);
        kle.add_mode(1.0, vec![0.0]);
        assert!((kle.energy_ratio(0) - 0.75).abs() < 1e-12);
        assert!((kle.energy_ratio(1) - 0.25).abs() < 1e-12);
    }
    #[test]
    fn test_kle_cumulative_energy() {
        let mut kle = KarhunenLoeveExpansion::new(3);
        kle.add_mode(6.0, vec![1.0]);
        kle.add_mode(3.0, vec![0.0]);
        kle.add_mode(1.0, vec![0.0]);
        assert!((kle.cumulative_energy(2) - 0.9).abs() < 1e-12);
        assert!((kle.cumulative_energy(3) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_kle_from_covariance_identity() {
        let cov = vec![1.0, 0.0, 0.0, 1.0];
        let kle = KarhunenLoeveExpansion::from_covariance(&cov, 2, 2, 100);
        let tv = kle.total_variance();
        assert!(tv > 0.5, "total variance={tv}");
    }
    #[test]
    fn test_sfem_problem_uniform_mesh() {
        let prob = StochasticFemProblem::new_uniform(5, 1.0, 100.0, 10.0, 0.5);
        assert_eq!(prob.n_nodes(), 5);
        assert!((prob.mean_field[0] - 100.0).abs() < 1e-10);
        assert!((prob.variance_field[0] - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_sfem_problem_covariance_matrix_diagonal() {
        let prob = StochasticFemProblem::new_uniform(3, 1.0, 0.0, 1.0, 1.0);
        let c = prob.covariance_matrix_exponential();
        for i in 0..3 {
            assert!((c[i][i] - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_sfem_problem_covariance_matrix_symmetric() {
        let prob = StochasticFemProblem::new_uniform(4, 2.0, 0.0, 2.0, 0.5);
        let c = prob.covariance_matrix_gaussian();
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (c[i][j] - c[j][i]).abs() < 1e-12,
                    "Unsymmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_sfem_problem_update_statistics() {
        let mut prob = StochasticFemProblem::new_uniform(3, 1.0, 0.0, 1.0, 1.0);
        let ensemble = vec![vec![1.0, 2.0, 3.0], vec![3.0, 4.0, 5.0]];
        prob.update_statistics_from_ensemble(&ensemble);
        assert!((prob.mean_field[0] - 2.0).abs() < 1e-12);
        assert!((prob.mean_field[1] - 3.0).abs() < 1e-12);
        assert!((prob.mean_field[2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_sfem_problem_cov_field() {
        let mut prob = StochasticFemProblem::new_uniform(2, 1.0, 10.0, 1.0, 1.0);
        prob.mean_field = vec![10.0, 5.0];
        prob.variance_field = vec![1.0, 1.0];
        let cov = prob.cov_field();
        assert!((cov[0] - 0.1).abs() < 1e-12);
        assert!((cov[1] - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_mc_solver_run_and_mean() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 10.0, 500.0);
        solver.run(1000);
        let mu = solver.mean();
        assert!(mu > 0.0, "mean={mu}");
        assert!((mu - 0.5).abs() < 0.1, "mean too far from 0.5: {mu}");
    }
    #[test]
    fn test_mc_solver_std_nonneg() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 50.0, 200.0);
        solver.run(500);
        assert!(solver.std_dev() >= 0.0);
    }
    #[test]
    fn test_mc_solver_cov_small_for_small_uncertainty() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 1.0, 500.0);
        solver.run(500);
        assert!(solver.cov() < 0.1, "CoV={}", solver.cov());
    }
    #[test]
    fn test_mc_solver_failure_probability() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 10.0, 500.0);
        solver.run(1000);
        let pf = solver.failure_probability(1.0);
        assert_eq!(pf, 0.0);
        let pf2 = solver.failure_probability(0.0);
        assert_eq!(pf2, 1.0);
    }
    #[test]
    fn test_mc_solver_percentile() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 10.0, 500.0);
        solver.run(1000);
        let p5 = solver.percentile(5.0);
        let p95 = solver.percentile(95.0);
        assert!(p5 <= p95, "p5={p5} p95={p95}");
    }
    #[test]
    fn test_mc_solver_confidence_interval() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 10.0, 500.0);
        solver.run(1000);
        let (lo, hi) = solver.confidence_interval_mean_95();
        let mu = solver.mean();
        assert!(lo <= mu && mu <= hi);
    }
    #[test]
    fn test_mc_solver_clear() {
        let mut solver = MonteCarloFemSolver::new(1000.0, 10.0, 500.0);
        solver.run(100);
        assert_eq!(solver.samples.len(), 100);
        solver.clear();
        assert_eq!(solver.samples.len(), 0);
    }
    #[test]
    fn test_response_mean() {
        let resp = StochasticResponse::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((resp.mean() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_response_variance() {
        let resp = StochasticResponse::new(vec![1.0, 1.0, 1.0]);
        assert!(resp.variance().abs() < 1e-12);
    }
    #[test]
    fn test_response_empty() {
        let resp = StochasticResponse::new(vec![]);
        assert_eq!(resp.mean(), 0.0);
        assert_eq!(resp.variance(), 0.0);
    }
    #[test]
    fn test_response_percentile_50() {
        let resp = StochasticResponse::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let p50 = resp.percentile(50.0);
        assert!((p50 - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_response_percentile_0_and_100() {
        let resp = StochasticResponse::new(vec![1.0, 5.0, 3.0]);
        assert!((resp.percentile(0.0) - 1.0).abs() < 1e-12);
        assert!((resp.percentile(100.0) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_response_cdf_below_all() {
        let resp = StochasticResponse::new(vec![2.0, 4.0, 6.0]);
        assert_eq!(resp.cdf(1.0), 0.0);
    }
    #[test]
    fn test_response_cdf_above_all() {
        let resp = StochasticResponse::new(vec![2.0, 4.0, 6.0]);
        assert!((resp.cdf(10.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_response_cdf_midpoint() {
        let resp = StochasticResponse::new(vec![1.0, 2.0, 3.0, 4.0]);
        assert!((resp.cdf(2.0) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_failure_prob_all_fail() {
        let resp = StochasticResponse::new(vec![5.0, 6.0, 7.0]);
        let pf = failure_probability_monte_carlo(&resp, 4.0);
        assert!((pf - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_failure_prob_none_fail() {
        let resp = StochasticResponse::new(vec![1.0, 2.0, 3.0]);
        let pf = failure_probability_monte_carlo(&resp, 10.0);
        assert_eq!(pf, 0.0);
    }
    #[test]
    fn test_failure_prob_in_unit_interval() {
        let resp = StochasticResponse::new(vec![1.0, 5.0, 3.0, 7.0, 2.0]);
        let pf = failure_probability_monte_carlo(&resp, 4.0);
        assert!(pf >= 0.0 && pf <= 1.0);
    }
    #[test]
    fn test_failure_prob_empty() {
        let resp = StochasticResponse::new(vec![]);
        assert_eq!(failure_probability_monte_carlo(&resp, 1.0), 0.0);
    }
    #[test]
    fn test_mc_fem_mean_positive() {
        let (mean_u, _std_u) = monte_carlo_fem(1000, 1000.0, 10.0, 500.0, 1);
        assert!(mean_u > 0.0, "mean displacement should be positive");
    }
    #[test]
    fn test_mc_fem_std_nonneg() {
        let (_mean_u, std_u) = monte_carlo_fem(500, 1000.0, 50.0, 200.0, 1);
        assert!(std_u >= 0.0);
    }
    #[test]
    fn test_mc_fem_zero_samples() {
        let (m, s) = monte_carlo_fem(0, 1000.0, 10.0, 100.0, 1);
        assert_eq!(m, 0.0);
        assert_eq!(s, 0.0);
    }
    #[test]
    fn test_phi_cdf_midpoint() {
        assert!((phi_cdf(0.0) - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_phi_cdf_positive_tail() {
        assert!(phi_cdf(3.0) > 0.998);
    }
    #[test]
    fn test_phi_cdf_negative_tail() {
        assert!(phi_cdf(-3.0) < 0.002);
    }
    #[test]
    fn test_phi_inv_roundtrip() {
        for p in [0.1, 0.25, 0.5, 0.75, 0.9] {
            let beta = phi_inv(p);
            let p2 = phi_cdf(beta);
            assert!((p2 - p).abs() < 1e-4, "roundtrip at p={p}: got {p2}");
        }
    }
    #[test]
    fn test_pce_struct_mean_variance() {
        let pce = PolynomialChaosExpansion::new(vec![5.0, 2.0, 1.0]);
        assert!((pce.mean() - 5.0).abs() < 1e-12);
        assert!((pce.variance() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_pce_struct_std_dev() {
        let pce = PolynomialChaosExpansion::new(vec![0.0, 3.0, 4.0]);
        assert!((pce.std_dev() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_pce_struct_sobol_indices_sum() {
        let pce = PolynomialChaosExpansion::new(vec![0.0, 3.0, 4.0]);
        let s = pce.sobol_indices();
        let total: f64 = s.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "Sobol sum={total}");
    }
    #[test]
    fn test_pce_struct_evaluate() {
        let pce = PolynomialChaosExpansion::new(vec![1.0, 0.0]);
        assert!((pce.evaluate(2.5) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_pce_from_quadrature_constant() {
        let pce = PolynomialChaosExpansion::from_quadrature(|_| 3.0, 2, 3);
        assert!((pce.mean() - 3.0).abs() < 0.1, "mean={}", pce.mean());
    }
    #[test]
    fn test_gauss_hermite_1pt() {
        let (nodes, weights) = gauss_hermite_quadrature(1);
        assert_eq!(nodes.len(), 1);
        assert!((nodes[0]).abs() < 1e-14);
        assert!((weights[0] - PI.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_gauss_hermite_2pt_sum_weights() {
        let (_, weights) = gauss_hermite_quadrature(2);
        let sum: f64 = weights.iter().sum();
        assert!((sum - PI.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_gauss_hermite_3pt_middle_node_zero() {
        let (nodes, _) = gauss_hermite_quadrature(3);
        assert!((nodes[1]).abs() < 1e-14);
    }
    #[test]
    fn test_gauss_hermite_5pt_has_5_points() {
        let (nodes, weights) = gauss_hermite_quadrature(5);
        assert_eq!(nodes.len(), 5);
        assert_eq!(weights.len(), 5);
    }
    #[test]
    fn test_morris_screening_returns_correct_length() {
        let sa = SensitivityAnalysis::new(3, 20);
        let f = |x: &[f64]| x[0] + 2.0 * x[1] + 3.0 * x[2];
        let (mu_star, sigma) = sa.morris_screening(&f, 0.5);
        assert_eq!(mu_star.len(), 3);
        assert_eq!(sigma.len(), 3);
    }
    #[test]
    fn test_morris_screening_nonneg_mu_star() {
        let sa = SensitivityAnalysis::new(2, 10);
        let f = |x: &[f64]| x[0] * x[0] + x[1];
        let (mu_star, _) = sa.morris_screening(&f, 0.5);
        for &v in &mu_star {
            assert!(v >= 0.0);
        }
    }
    #[test]
    fn test_saltelli_sobol_returns_correct_length() {
        let sa = SensitivityAnalysis::new(2, 50);
        let f = |x: &[f64]| x[0] + x[1];
        let (s1, st) = sa.saltelli_sobol(&f);
        assert_eq!(s1.len(), 2);
        assert_eq!(st.len(), 2);
    }
    #[test]
    fn test_linear_sobol_indices_sum_to_one() {
        let coeffs = [1.0, 2.0, 3.0];
        let std_devs = [1.0, 1.0, 1.0];
        let s = SensitivityAnalysis::linear_sobol_indices(&coeffs, &std_devs);
        let total: f64 = s.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "sum={total}");
    }
    #[test]
    fn test_linear_sobol_indices_proportional_to_coeff_sq() {
        let coeffs = [0.0, 1.0, 2.0];
        let std_devs = [1.0, 1.0, 1.0];
        let s = SensitivityAnalysis::linear_sobol_indices(&coeffs, &std_devs);
        assert!(s[0].abs() < 1e-12);
        assert!((s[1] - 0.2).abs() < 1e-12);
        assert!((s[2] - 0.8).abs() < 1e-12);
    }
    #[test]
    fn test_reliability_form_beta_positive() {
        let rel = ReliabilityFem::new(vec![5.0], vec![1.0], vec![10.0, -1.0]);
        let beta = rel.form_beta();
        assert!((beta - 5.0).abs() < 1e-10, "beta={beta}");
    }
    #[test]
    fn test_reliability_form_pf_small_for_large_beta() {
        let rel = ReliabilityFem::new(vec![5.0], vec![1.0], vec![10.0, -1.0]);
        let pf = rel.form_pf();
        assert!(pf < 1e-4, "pf={pf}");
    }
    #[test]
    fn test_reliability_lsf_evaluation() {
        let rel = ReliabilityFem::new(vec![1.0, 1.0], vec![0.1, 0.1], vec![1.0, 2.0, -3.0]);
        let g = rel.lsf(&[1.0, 0.0]);
        assert!((g - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_reliability_safety_factor() {
        let rel = ReliabilityFem::new(vec![200.0, 100.0], vec![10.0, 5.0], vec![0.0]);
        assert!((rel.safety_factor() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_reliability_importance_vector_unit_length() {
        let rel = ReliabilityFem::new(vec![5.0, 3.0], vec![1.0, 1.0], vec![10.0, -2.0, -3.0]);
        let alpha = rel.importance_vector();
        let len: f64 = alpha.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((len - 1.0).abs() < 1e-10, "len={len}");
    }
    #[test]
    fn test_reliability_sorm_breitung_no_curvatures() {
        let rel = ReliabilityFem::new(vec![5.0], vec![1.0], vec![10.0, -1.0]);
        let sorm = rel.sorm_pf_breitung(&[]);
        let form = rel.form_pf();
        assert!((sorm - form).abs() < 1e-12);
    }
    #[test]
    fn test_reliability_mc_pf_ballpark() {
        let rel = ReliabilityFem::new(vec![5.0], vec![1.0], vec![-3.0, 1.0]);
        let pf = rel.monte_carlo_pf(2000);
        assert!(pf < 0.15, "pf={pf}");
    }
    #[test]
    fn test_box_muller_basic() {
        let z = box_muller(0.5, 0.5);
        assert!(z.is_finite());
    }
}
