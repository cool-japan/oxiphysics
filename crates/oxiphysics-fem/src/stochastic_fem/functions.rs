//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

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
/// Gauss-Hermite quadrature nodes and weights for the `n`-point rule.
///
/// Returns `(nodes, weights)` with the nodes in ascending order, for the
/// physicists' integral
///
/// ∫ f(x) exp(−x²) dx ≈ Σᵢ wᵢ f(xᵢ),
///
/// which is exact for every polynomial `f` of degree ≤ 2n − 1.  The weights
/// sum to √π.  To integrate against the standard normal density use the
/// substitution ξ = √2·x together with the factor π^(−1/2):
/// E[g(ξ)] = π^(−1/2) Σᵢ wᵢ g(√2·xᵢ).
///
/// The rule is computed for any `n ≥ 1` with the Golub–Welsch algorithm: the
/// nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix with
/// zero diagonal and off-diagonal entries √(k/2), and the weights are
/// √π times the squared first components of the normalised eigenvectors.
/// The nodes are then polished with Newton steps on the orthonormal Hermite
/// recurrence and the weights re-evaluated from the Christoffel formula
/// wᵢ = 1 / (n · hₙ₋₁(xᵢ)²), so the tabulated low-order rules are accurate to
/// machine precision.
///
/// Returns empty vectors for `n == 0`.
pub fn gauss_hermite_quadrature(n: usize) -> (Vec<f64>, Vec<f64>) {
    if n == 0 {
        return (vec![], vec![]);
    }
    let mut diag = vec![0.0_f64; n];
    let mut off_diag: Vec<f64> = (0..n)
        .map(|k| {
            if k + 1 < n {
                ((k + 1) as f64 / 2.0).sqrt()
            } else {
                0.0
            }
        })
        .collect();
    let mut first_row = vec![0.0_f64; n];
    first_row[0] = 1.0;
    if !tridiagonal_ql_first_row(&mut diag, &mut off_diag, &mut first_row) {
        return (vec![], vec![]);
    }
    let mut pairs: Vec<(f64, f64)> = diag
        .iter()
        .zip(first_row.iter())
        .map(|(&x, &z)| (x, PI.sqrt() * z * z))
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut nodes: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let mut weights: Vec<f64> = pairs.iter().map(|p| p.1).collect();
    // Enforce the exact symmetry x_{n-1-i} = -x_i, w_{n-1-i} = w_i.
    for i in 0..n / 2 {
        let j = n - 1 - i;
        let x = 0.5 * (nodes[j] - nodes[i]);
        let w = 0.5 * (weights[i] + weights[j]);
        nodes[i] = -x;
        nodes[j] = x;
        weights[i] = w;
        weights[j] = w;
    }
    if n % 2 == 1 {
        nodes[n / 2] = 0.0;
    }
    // Newton polish on the orthonormal Hermite recurrence, then Christoffel weights.
    for (x, w) in nodes.iter_mut().zip(weights.iter_mut()) {
        let mut refined = *x;
        for _ in 0..4 {
            let (h_n, h_nm1) = orthonormal_hermite_pair(n, refined);
            let derivative = (2.0 * n as f64).sqrt() * h_nm1;
            if derivative == 0.0 || !derivative.is_finite() || !h_n.is_finite() {
                break;
            }
            refined -= h_n / derivative;
        }
        if refined.is_finite() && (refined - *x).abs() <= 1.0e-6 * (1.0 + x.abs()) {
            let (_, h_nm1) = orthonormal_hermite_pair(n, refined);
            let christoffel = 1.0 / (n as f64 * h_nm1 * h_nm1);
            if christoffel.is_finite() && christoffel > 0.0 {
                *x = refined;
                *w = christoffel;
            }
        }
    }
    (nodes, weights)
}
/// Values `(hₙ(x), hₙ₋₁(x))` of the Hermite polynomials orthonormal with
/// respect to the weight exp(−x²), by the stable three-term recurrence
/// h₀ = π^(−1/4), h_{k+1} = √(2/(k+1))·x·h_k − √(k/(k+1))·h_{k−1}.
fn orthonormal_hermite_pair(n: usize, x: f64) -> (f64, f64) {
    let mut h_prev = 0.0_f64;
    let mut h = PI.powf(-0.25);
    for k in 0..n {
        let kf = k as f64;
        let h_next = (2.0 / (kf + 1.0)).sqrt() * x * h - (kf / (kf + 1.0)).sqrt() * h_prev;
        h_prev = h;
        h = h_next;
    }
    (h, h_prev)
}
/// Eigenvalues of a symmetric tridiagonal matrix by the implicit QL algorithm
/// with Wilkinson shifts, together with the first row of the orthogonal
/// eigenvector matrix.
///
/// On entry `d` holds the diagonal, `e[i]` the off-diagonal entry coupling
/// rows `i` and `i + 1` (`e[n - 1]` is unused and must be zero) and `z` the
/// first row of the identity.  On success `d` holds the eigenvalues and
/// `z[j]` the first component of the eigenvector of `d[j]`.  Returns `false`
/// if the iteration fails to converge.
fn tridiagonal_ql_first_row(d: &mut [f64], e: &mut [f64], z: &mut [f64]) -> bool {
    const MAX_ITERATIONS: usize = 60;
    let n = d.len();
    for l in 0..n {
        let mut iterations = 0;
        loop {
            let mut m = l;
            while m + 1 < n {
                let scale = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * scale {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iterations += 1;
            if iterations > MAX_ITERATIONS {
                return false;
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + r.copysign(g));
            let (mut s, mut c, mut p) = (1.0_f64, 1.0_f64, 0.0_f64);
            let mut underflow = false;
            let mut i = m;
            while i > l {
                i -= 1;
                let f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    underflow = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                let z_next = z[i + 1];
                z[i + 1] = s * z[i] + c * z_next;
                z[i] = c * z[i] - s * z_next;
            }
            if underflow {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    true
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
        for (i, row) in c.iter().enumerate() {
            assert!((row[i] - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_sfem_problem_covariance_matrix_symmetric() {
        let prob = StochasticFemProblem::new_uniform(4, 2.0, 0.0, 2.0, 0.5);
        let c = prob.covariance_matrix_gaussian();
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - c[j][i]).abs() < 1e-12, "Unsymmetric at ({i},{j})");
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
        assert!((0.0..=1.0).contains(&pf));
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
        assert!((pce.mean() - 3.0).abs() < 1e-12, "mean={}", pce.mean());
        assert!(pce.coefficients[1].abs() < 1e-12);
        assert!(pce.coefficients[2].abs() < 1e-12);
    }
    #[test]
    fn test_pce_from_quadrature_linear_model_has_unit_first_coefficient() {
        // xi = He_1(xi): c = [0, 1, 0].  The mis-scaled projection gave c1 = 0.5.
        let pce = PolynomialChaosExpansion::from_quadrature(|xi| xi, 2, 3);
        assert!(pce.coefficients[0].abs() < 1e-13);
        assert!(
            (pce.coefficients[1] - 1.0).abs() < 1e-13,
            "c1 = {}",
            pce.coefficients[1]
        );
        assert!(pce.coefficients[2].abs() < 1e-13);
    }
    #[test]
    fn test_pce_from_quadrature_quadratic_model() {
        // xi^2 = He_2 + He_0: c = [1, 0, 1].
        let pce = PolynomialChaosExpansion::from_quadrature(|xi| xi * xi, 2, 3);
        let expected = [1.0, 0.0, 1.0];
        for (k, (got, want)) in pce.coefficients.iter().zip(expected.iter()).enumerate() {
            assert!((got - want).abs() < 1e-13, "c{k} = {got}, expected {want}");
        }
    }
    #[test]
    fn test_pce_from_quadrature_cubic_model_in_all_supported_rules() {
        // xi^3 = He_3 + 3 He_1: c = [0, 3, 0, 1, 0].  Exact once 2n - 1 >= 3 + 4.
        for quad_pts in 4..=9 {
            let pce = PolynomialChaosExpansion::from_quadrature(|xi| xi.powi(3), 4, quad_pts);
            let expected = [0.0, 3.0, 0.0, 1.0, 0.0];
            for (k, (got, want)) in pce.coefficients.iter().zip(expected.iter()).enumerate() {
                assert!(
                    (got - want).abs() < 1e-11,
                    "{quad_pts} points: c{k} = {got}, expected {want}"
                );
            }
        }
    }
    #[test]
    fn test_pce_from_quadrature_exponential_model_matches_closed_form() {
        // exp(a xi) = exp(a^2 / 2) * sum_n a^n / n! * He_n(xi).
        let a = 0.5_f64;
        let max_order = 5;
        let pce = PolynomialChaosExpansion::from_quadrature(|xi| (a * xi).exp(), max_order, 12);
        let mut factorial = 1.0;
        for n in 0..=max_order {
            if n > 0 {
                factorial *= n as f64;
            }
            let expected = (0.5 * a * a).exp() * a.powi(n as i32) / factorial;
            assert!(
                (pce.coefficients[n] - expected).abs() < 1e-12,
                "c{n} = {}, expected {expected}",
                pce.coefficients[n]
            );
        }
        // The projected mean equals E[exp(a xi)] = exp(a^2 / 2).
        assert!((pce.mean() - (0.5 * a * a).exp()).abs() < 1e-12);
    }
    #[test]
    fn test_pce_from_quadrature_zero_points_gives_zero_coefficients() {
        let pce = PolynomialChaosExpansion::from_quadrature(|xi| xi, 3, 0);
        assert_eq!(pce.coefficients, vec![0.0; 4]);
    }
    #[test]
    fn test_pce_from_quadrature_evaluates_model_once_per_node() {
        let calls = std::cell::Cell::new(0usize);
        let _ = PolynomialChaosExpansion::from_quadrature(
            |xi| {
                calls.set(calls.get() + 1);
                xi
            },
            6,
            5,
        );
        assert_eq!(calls.get(), 5);
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
    fn test_gauss_hermite_5pt_matches_reference_values() {
        let (nodes, weights) = gauss_hermite_quadrature(5);
        let ref_nodes = [
            -2.020_182_870_456_086,
            -0.958_572_464_613_818_5,
            0.0,
            0.958_572_464_613_818_5,
            2.020_182_870_456_086,
        ];
        let ref_weights = [
            0.019_953_242_059_045_91,
            0.393_619_323_152_241_2,
            0.945_308_720_482_941_9,
            0.393_619_323_152_241_2,
            0.019_953_242_059_045_91,
        ];
        for k in 0..5 {
            assert!(
                (nodes[k] - ref_nodes[k]).abs() < 1e-14,
                "node {k}: {} vs {}",
                nodes[k],
                ref_nodes[k]
            );
            assert!(
                (weights[k] - ref_weights[k]).abs() < 1e-14,
                "weight {k}: {} vs {}",
                weights[k],
                ref_weights[k]
            );
        }
    }
    #[test]
    fn test_gauss_hermite_low_order_rules_match_closed_forms() {
        let (n2, w2) = gauss_hermite_quadrature(2);
        assert!((n2[1] - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-15);
        assert!((n2[0] + std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-15);
        assert!((w2[0] - PI.sqrt() / 2.0).abs() < 1e-15);
        let (n3, w3) = gauss_hermite_quadrature(3);
        assert!((n3[2] - 1.5_f64.sqrt()).abs() < 1e-15);
        assert_eq!(n3[1], 0.0);
        assert!((w3[1] - 2.0 * PI.sqrt() / 3.0).abs() < 1e-15);
        assert!((w3[0] - PI.sqrt() / 6.0).abs() < 1e-15);
        let (n4, w4) = gauss_hermite_quadrature(4);
        assert!((n4[3] - 1.650_680_123_885_784_4).abs() < 1e-14);
        assert!((n4[2] - 0.524_647_623_275_290_4).abs() < 1e-14);
        assert!((w4[3] - 0.081_312_835_447_245_24).abs() < 1e-14);
        assert!((w4[2] - 0.804_914_090_005_512_7).abs() < 1e-14);
    }
    #[test]
    fn test_gauss_hermite_weights_sum_to_sqrt_pi_for_every_order() {
        for n in 1..=30 {
            let (nodes, weights) = gauss_hermite_quadrature(n);
            assert_eq!(nodes.len(), n);
            assert_eq!(weights.len(), n);
            let sum: f64 = weights.iter().sum();
            assert!(
                (sum - PI.sqrt()).abs() < 1e-12 * PI.sqrt(),
                "n = {n}: sum of weights = {sum}"
            );
            assert!(weights.iter().all(|&w| w > 0.0));
        }
    }
    #[test]
    fn test_gauss_hermite_nodes_ascending_and_symmetric() {
        for n in 1..=20 {
            let (nodes, weights) = gauss_hermite_quadrature(n);
            for pair in nodes.windows(2) {
                assert!(pair[0] < pair[1], "n = {n}: nodes not strictly ascending");
            }
            for i in 0..n {
                assert_eq!(nodes[i], -nodes[n - 1 - i], "n = {n}: node symmetry");
                assert_eq!(weights[i], weights[n - 1 - i], "n = {n}: weight symmetry");
            }
        }
    }
    #[test]
    fn test_gauss_hermite_integrates_even_moments_exactly() {
        // ∫ x^(2k) exp(-x²) dx = Γ(k + 1/2) = (2k - 1)!! √π / 2^k, exact for 2k <= 2n - 1.
        for n in 1..=12 {
            let (nodes, weights) = gauss_hermite_quadrature(n);
            let mut exact = PI.sqrt();
            for k in 0..n {
                if k > 0 {
                    exact *= (2 * k - 1) as f64 / 2.0;
                }
                let approx: f64 = nodes
                    .iter()
                    .zip(weights.iter())
                    .map(|(&x, &w)| w * x.powi(2 * k as i32))
                    .sum();
                assert!(
                    (approx - exact).abs() < 1e-11 * exact.max(1.0),
                    "n = {n}, moment 2k = {}: {approx} vs {exact}",
                    2 * k
                );
            }
        }
    }
    #[test]
    fn test_gauss_hermite_odd_moments_vanish() {
        let (nodes, weights) = gauss_hermite_quadrature(7);
        for power in [1, 3, 5, 7, 9, 11, 13] {
            let approx: f64 = nodes
                .iter()
                .zip(weights.iter())
                .map(|(&x, &w)| w * x.powi(power))
                .sum();
            assert!(approx.abs() < 1e-12, "odd moment {power}: {approx}");
        }
    }
    #[test]
    fn test_gauss_hermite_zero_points_is_empty() {
        let (nodes, weights) = gauss_hermite_quadrature(0);
        assert!(nodes.is_empty());
        assert!(weights.is_empty());
    }
    #[test]
    fn test_gauss_hermite_standard_normal_expectations() {
        // E[xi^2] = 1, E[xi^4] = 3, E[xi^6] = 15 through the sqrt(2) substitution.
        let (nodes, weights) = gauss_hermite_quadrature(6);
        let expect = |power: i32| -> f64 {
            nodes
                .iter()
                .zip(weights.iter())
                .map(|(&x, &w)| w * (2.0_f64.sqrt() * x).powi(power))
                .sum::<f64>()
                / PI.sqrt()
        };
        assert!((expect(0) - 1.0).abs() < 1e-13);
        assert!((expect(2) - 1.0).abs() < 1e-13);
        assert!((expect(4) - 3.0).abs() < 1e-12);
        assert!((expect(6) - 15.0).abs() < 1e-11);
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
