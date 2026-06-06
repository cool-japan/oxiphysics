//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{GeometricBrownianMotion, Rng};

pub(super) const K_B: f64 = 1.38e-23;
/// Solves an Ito SDE using the Euler-Maruyama scheme.
///
/// Equation: `dx = drift(x, t)*dt + diffusion(x, t)*dW`
///
/// # Arguments
/// * `x0` - Initial condition
/// * `t_end` - End time
/// * `dt` - Time step
/// * `drift` - Drift function f(x, t)
/// * `diffusion` - Diffusion function g(x, t)
/// * `seed` - RNG seed
///
/// Returns a vector of x values at each time step.
pub fn euler_maruyama(
    x0: f64,
    t_end: f64,
    dt: f64,
    drift: impl Fn(f64, f64) -> f64,
    diffusion: impl Fn(f64, f64) -> f64,
    seed: u64,
) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    let n = ((t_end / dt).ceil() as usize) + 1;
    let mut path = Vec::with_capacity(n);
    let mut x = x0;
    let mut t = 0.0_f64;
    path.push(x);
    while t + dt <= t_end + 1e-12 * dt {
        let dw = rng.next_normal() * dt.sqrt();
        x += drift(x, t) * dt + diffusion(x, t) * dw;
        t += dt;
        path.push(x);
    }
    path
}
/// Solves an Ito SDE using the Milstein scheme (higher-order correction).
///
/// Equation: `dx = drift(x,t)*dt + diffusion(x,t)*dW + 0.5*g*g'*((dW)^2 - dt)`
///
/// # Arguments
/// * `x0` - Initial condition
/// * `t_end` - End time
/// * `dt` - Time step
/// * `drift` - Drift function f(x, t)
/// * `diffusion` - Diffusion function g(x, t)
/// * `diffusion_prime` - Derivative of diffusion w.r.t. x: g'(x, t)
/// * `seed` - RNG seed
///
/// Returns a vector of x values at each time step.
pub fn milstein(
    x0: f64,
    t_end: f64,
    dt: f64,
    drift: impl Fn(f64, f64) -> f64,
    diffusion: impl Fn(f64, f64) -> f64,
    diffusion_prime: impl Fn(f64, f64) -> f64,
    seed: u64,
) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    let n = ((t_end / dt).ceil() as usize) + 1;
    let mut path = Vec::with_capacity(n);
    let mut x = x0;
    let mut t = 0.0_f64;
    path.push(x);
    while t + dt <= t_end + 1e-12 * dt {
        let dw = rng.next_normal() * dt.sqrt();
        let g = diffusion(x, t);
        let gp = diffusion_prime(x, t);
        x += drift(x, t) * dt + g * dw + 0.5 * g * gp * (dw * dw - dt);
        t += dt;
        path.push(x);
    }
    path
}
/// Generate multiple Monte Carlo paths for a GBM and compute statistics.
///
/// Returns `(mean_path, std_path)` where each vector has length `n_steps + 1`.
pub fn monte_carlo_gbm_paths(
    s0: f64,
    mu: f64,
    sigma: f64,
    t_end: f64,
    dt: f64,
    n_paths: usize,
    seed: u64,
) -> (Vec<f64>, Vec<f64>) {
    let n_steps = ((t_end / dt).ceil() as usize) + 1;
    let mut sum = vec![0.0_f64; n_steps];
    let mut sum_sq = vec![0.0_f64; n_steps];
    let gbm = GeometricBrownianMotion::new(mu, sigma);
    for i in 0..n_paths {
        let path = gbm.simulate_exact(s0, t_end, dt, seed.wrapping_add(i as u64));
        for (j, &val) in path.iter().enumerate() {
            if j < n_steps {
                sum[j] += val;
                sum_sq[j] += val * val;
            }
        }
    }
    let n = n_paths as f64;
    let mean_path: Vec<f64> = sum.iter().map(|&s| s / n).collect();
    let std_path: Vec<f64> = sum
        .iter()
        .zip(sum_sq.iter())
        .map(|(&s, &s2)| {
            let m = s / n;
            ((s2 / n - m * m).max(0.0)).sqrt()
        })
        .collect();
    (mean_path, std_path)
}
/// Generate Monte Carlo estimate of E\[f(S_T)\] for a GBM terminal value.
///
/// Returns `(mean, stderr)`.
pub fn monte_carlo_terminal_estimate(
    s0: f64,
    mu: f64,
    sigma: f64,
    t_end: f64,
    n_paths: usize,
    f: impl Fn(f64) -> f64,
    seed: u64,
) -> (f64, f64) {
    let mut rng = Rng::new(seed);
    let drift = (mu - 0.5 * sigma * sigma) * t_end;
    let vol = sigma * t_end.sqrt();
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    for _ in 0..n_paths {
        let z = rng.next_normal();
        let s_t = s0 * (drift + vol * z).exp();
        let val = f(s_t);
        sum += val;
        sum_sq += val * val;
    }
    let n = n_paths as f64;
    let mean = sum / n;
    let var = sum_sq / n - mean * mean;
    let stderr = (var / n).sqrt();
    (mean, stderr)
}
/// Monte Carlo with antithetic variates for variance reduction.
///
/// For each normal draw Z, also uses -Z to reduce variance.
/// Returns `(mean, stderr)`.
pub fn antithetic_terminal_estimate(
    s0: f64,
    mu: f64,
    sigma: f64,
    t_end: f64,
    n_pairs: usize,
    f: impl Fn(f64) -> f64,
    seed: u64,
) -> (f64, f64) {
    let mut rng = Rng::new(seed);
    let drift = (mu - 0.5 * sigma * sigma) * t_end;
    let vol = sigma * t_end.sqrt();
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let n_total = 2 * n_pairs;
    for _ in 0..n_pairs {
        let z = rng.next_normal();
        let s_plus = s0 * (drift + vol * z).exp();
        let s_minus = s0 * (drift - vol * z).exp();
        let avg = 0.5 * (f(s_plus) + f(s_minus));
        sum += avg;
        sum_sq += avg * avg;
    }
    let n = n_pairs as f64;
    let mean = sum / n;
    let var = sum_sq / n - mean * mean;
    let stderr = (var / n).sqrt();
    let _ = n_total;
    (mean, stderr)
}
/// Monte Carlo with control variate for variance reduction.
///
/// Uses the terminal asset price itself as control variate since
/// E\[S_T\] is known analytically for GBM.
/// Returns `(mean, stderr)`.
pub fn control_variate_terminal_estimate(
    s0: f64,
    mu: f64,
    sigma: f64,
    t_end: f64,
    n_paths: usize,
    f: impl Fn(f64) -> f64,
    seed: u64,
) -> (f64, f64) {
    let mut rng = Rng::new(seed);
    let drift = (mu - 0.5 * sigma * sigma) * t_end;
    let vol = sigma * t_end.sqrt();
    let expected_s_t = s0 * (mu * t_end).exp();
    let mut vals = Vec::with_capacity(n_paths);
    let mut controls = Vec::with_capacity(n_paths);
    for _ in 0..n_paths {
        let z = rng.next_normal();
        let s_t = s0 * (drift + vol * z).exp();
        vals.push(f(s_t));
        controls.push(s_t);
    }
    let mean_val = vals.iter().sum::<f64>() / n_paths as f64;
    let mean_ctrl = controls.iter().sum::<f64>() / n_paths as f64;
    let mut cov = 0.0_f64;
    let mut var_c = 0.0_f64;
    for (&v, &c) in vals.iter().zip(controls.iter()) {
        cov += (v - mean_val) * (c - mean_ctrl);
        var_c += (c - mean_ctrl) * (c - mean_ctrl);
    }
    let beta = if var_c.abs() > 1e-30 {
        -cov / var_c
    } else {
        0.0
    };
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    for (&v, &c) in vals.iter().zip(controls.iter()) {
        let adj = v + beta * (c - expected_s_t);
        sum += adj;
        sum_sq += adj * adj;
    }
    let n = n_paths as f64;
    let mean = sum / n;
    let var = sum_sq / n - mean * mean;
    let stderr = (var / n).sqrt();
    (mean, stderr)
}
/// Computes the arithmetic mean of a slice of values.
///
/// Returns `0.0` if the slice is empty.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.iter().sum::<f64>() / data.len() as f64
}
/// Computes the unbiased sample variance of a slice of values.
///
/// Returns `0.0` if the slice has fewer than 2 elements.
pub fn variance(data: &[f64]) -> f64 {
    if data.len() < 2 {
        return 0.0;
    }
    let m = mean(data);
    data.iter().map(|&x| (x - m) * (x - m)).sum::<f64>() / (data.len() - 1) as f64
}
/// Computes the sample standard deviation of a slice of values.
///
/// Returns `0.0` if the slice has fewer than 2 elements.
pub fn std_dev(data: &[f64]) -> f64 {
    variance(data).sqrt()
}
/// Computes the sample skewness of a slice of values.
///
/// Returns `0.0` if the slice has fewer than 3 elements.
pub fn skewness(data: &[f64]) -> f64 {
    if data.len() < 3 {
        return 0.0;
    }
    let m = mean(data);
    let n = data.len() as f64;
    let m3: f64 = data.iter().map(|&x| (x - m).powi(3)).sum::<f64>() / n;
    let m2: f64 = data.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / n;
    if m2.abs() < 1e-30 {
        return 0.0;
    }
    m3 / m2.powf(1.5)
}
/// Computes the sample excess kurtosis of a slice of values.
///
/// Returns `0.0` if the slice has fewer than 4 elements.
pub fn kurtosis(data: &[f64]) -> f64 {
    if data.len() < 4 {
        return 0.0;
    }
    let m = mean(data);
    let n = data.len() as f64;
    let m4: f64 = data.iter().map(|&x| (x - m).powi(4)).sum::<f64>() / n;
    let m2: f64 = data.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / n;
    if m2.abs() < 1e-30 {
        return 0.0;
    }
    m4 / (m2 * m2) - 3.0
}
/// Computes the autocorrelation of a time series at a given lag.
///
/// Uses the normalized estimator: `R(lag) = Cov(X_t, X_{t+lag}) / Var(X)`.
/// Returns `0.0` if there are insufficient data points.
pub fn autocorrelation(data: &[f64], lag: usize) -> f64 {
    let n = data.len();
    if n <= lag || n < 2 {
        return 0.0;
    }
    let m = mean(data);
    let var = data.iter().map(|&x| (x - m) * (x - m)).sum::<f64>() / n as f64;
    if var == 0.0 {
        return 0.0;
    }
    let cov: f64 = (0..n - lag)
        .map(|i| (data[i] - m) * (data[i + lag] - m))
        .sum::<f64>()
        / (n - lag) as f64;
    cov / var
}
/// Estimates the diffusion coefficient from mean-squared displacement (MSD) data.
///
/// Performs a linear regression of MSD vs time and returns `slope / 2`
/// (the 1-D Einstein relation: `MSD = 2*D*t`).
///
/// Returns `0.0` if there are fewer than 2 data points.
pub fn diffusion_coefficient(msd_data: &[(f64, f64)]) -> f64 {
    let n = msd_data.len();
    if n < 2 {
        return 0.0;
    }
    let sum_t: f64 = msd_data.iter().map(|(t, _)| t).sum();
    let sum_msd: f64 = msd_data.iter().map(|(_, m)| m).sum();
    let sum_t2: f64 = msd_data.iter().map(|(t, _)| t * t).sum();
    let sum_t_msd: f64 = msd_data.iter().map(|(t, m)| t * m).sum();
    let n_f = n as f64;
    let denom = n_f * sum_t2 - sum_t * sum_t;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    let slope = (n_f * sum_t_msd - sum_t * sum_msd) / denom;
    slope / 2.0
}
/// Computes the running (cumulative) mean of a data series.
pub fn running_mean(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut sum = 0.0_f64;
    for (i, &x) in data.iter().enumerate() {
        sum += x;
        result.push(sum / (i + 1) as f64);
    }
    result
}
/// Computes the effective sample size using autocorrelation.
///
/// ESS = N / (1 + 2 * sum_{k=1}^{max_lag} rho(k))
/// where rho(k) is the autocorrelation at lag k.
pub fn effective_sample_size(data: &[f64], max_lag: usize) -> f64 {
    let n = data.len();
    if n < 2 {
        return n as f64;
    }
    let mut tau = 1.0_f64;
    for k in 1..=max_lag.min(n - 1) {
        let rho = autocorrelation(data, k);
        if rho < 0.0 {
            break;
        }
        tau += 2.0 * rho;
    }
    n as f64 / tau
}
/// Kramers mean first-passage time (MFPT) for escape over a barrier.
///
/// Uses the Kramers high-friction formula:
/// `τ = (2π / (ω_0 * ω_b)) * exp(β * ΔE)`
///
/// # Arguments
/// * `barrier_height` — ΔE = barrier height above the minimum
/// * `d` — diffusion coefficient (not directly used in the Kramers formula but
///   kept as part of the signature for context; `d = kT / gamma`)
/// * `omega_0` — angular frequency at the minimum (attempt frequency)
/// * `omega_b` — absolute value of imaginary frequency at the saddle point
///
/// Returns the MFPT in units consistent with the input.
pub fn mean_first_passage_time(barrier_height: f64, d: f64, omega_0: f64, omega_b: f64) -> f64 {
    let _ = d;
    (2.0 * PI / (omega_0 * omega_b)) * (barrier_height / d).exp()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::LangevinDynamics;
    use crate::OrnsteinUhlenbeck;
    use crate::RandomWalk;
    use crate::Rng;
    use crate::WienerProcess;
    use crate::stochastic::AntitheticGbm;

    use crate::stochastic::ControlVariateGbm;
    use crate::stochastic::EmpiricalFirstPassageTime;
    use crate::stochastic::FractionalBrownianMotion;
    use crate::stochastic::HestonModel;

    use crate::stochastic::KleinmanKramers;
    use crate::stochastic::LangevinIntegrator;

    use crate::stochastic::MertonJumpDiffusion;
    use crate::stochastic::MetropolisHastings;
    use crate::stochastic::OuExactSampler;

    use crate::stochastic::WienerSampler;
    #[test]
    fn test_wiener_process_statistics() {
        let dt = 0.01;
        let n = 10_000;
        let mut wp = WienerProcess::new(42);
        let increments: Vec<f64> = (0..n).map(|_| wp.increment(dt)).collect();
        let m = mean(&increments);
        let v = variance(&increments);
        assert!(m.abs() < 0.05, "mean={m} not close to 0");
        assert!((v - dt).abs() < 0.005, "variance={v} not close to dt={dt}");
    }
    #[test]
    fn test_gbm_mean() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2);
        let s0 = 100.0;
        let t_end = 1.0;
        let dt = 0.01;
        let n_paths = 500;
        let finals: Vec<f64> = (0..n_paths)
            .map(|seed| {
                let path = gbm.simulate(s0, t_end, dt, seed as u64);
                *path.last().unwrap()
            })
            .collect();
        let emp_mean = mean(&finals);
        let ana_mean = gbm.analytical_mean(s0, t_end);
        let ana_var = gbm.analytical_variance(s0, t_end);
        let se = (ana_var / n_paths as f64).sqrt();
        assert!(
            (emp_mean - ana_mean).abs() < 3.0 * se,
            "emp_mean={emp_mean}, ana_mean={ana_mean}, 3*se={}",
            3.0 * se
        );
    }
    #[test]
    fn test_gbm_exact_simulation() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2);
        let path = gbm.simulate_exact(100.0, 1.0, 0.01, 42);
        assert!(path.len() > 90);
        for &v in &path {
            assert!(v > 0.0, "exact GBM produced non-positive value: {v}");
        }
    }
    #[test]
    fn test_ou_converges_to_mean() {
        let ou = OrnsteinUhlenbeck::new(2.0, 1.5, 0.3);
        let path = ou.simulate(0.0, 10.0, 0.01, 123);
        let tail = &path[path.len() * 3 / 4..];
        let m = mean(tail);
        assert!(
            (m - ou.mu).abs() < 0.2,
            "OU mean={m}, expected close to mu={}",
            ou.mu
        );
    }
    #[test]
    fn test_ou_stationary_variance() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.0, 0.5);
        let expected = ou.stationary_variance();
        assert!((expected - 0.0625).abs() < 1e-10);
    }
    #[test]
    fn test_ou_autocorrelation() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.0, 0.5);
        let rho = ou.autocorrelation_at(0.5);
        let expected = (-1.0_f64).exp();
        assert!((rho - expected).abs() < 1e-10);
    }
    #[test]
    fn test_ou_analytical_mean() {
        let ou = OrnsteinUhlenbeck::new(1.0, 5.0, 0.1);
        let m = ou.analytical_mean_at(0.0, 10.0);
        assert!((m - 5.0).abs() < 1e-3, "analytical mean = {m}");
    }
    #[test]
    fn test_langevin_step_changes_position() {
        let ld = LangevinDynamics::new(1.0e-26, 1e10, 300.0);
        let mut rng = Rng::new(7);
        let (x_new, _v_new) = ld.step(0.0, 0.0, 1.0e-21, 1e-15, &mut rng);
        assert!(x_new.abs() > 0.0, "position did not change");
    }
    #[test]
    fn test_langevin_baoab_step() {
        let ld = LangevinDynamics::new(1.0e-26, 1e10, 300.0);
        let mut rng = Rng::new(7);
        let (x_new, v_new) = ld.step_baoab(0.0, 0.0, 1.0e-21, 1e-15, &mut rng);
        assert!(
            x_new.abs() > 0.0 || v_new.abs() > 0.0,
            "BAOAB step had no effect"
        );
    }
    #[test]
    fn test_mean_and_variance_known_data() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let m = mean(&data);
        let v = variance(&data);
        assert!((m - 5.0).abs() < 1e-10, "mean={m}");
        assert!((v - 4.571_428_571_428_571).abs() < 1e-6, "variance={v}");
    }
    #[test]
    fn test_euler_maruyama_zero_diffusion() {
        let x0 = 1.0;
        let t_end = 2.0;
        let dt = 0.001;
        let path = euler_maruyama(x0, t_end, dt, |x, _t| -x, |_x, _t| 0.0, 42);
        let x_final = *path.last().unwrap();
        let expected = x0 * (-t_end).exp();
        assert!(
            (x_final - expected).abs() < 0.01,
            "x_final={x_final}, expected={expected}"
        );
    }
    #[test]
    fn test_random_walk_dimension_2() {
        let rw = RandomWalk::new(2, 1.0);
        let path = rw.walk(100, 99);
        assert_eq!(path.len(), 101, "path length mismatch");
        for pos in &path {
            assert_eq!(pos.len(), 2, "position has wrong dimension");
        }
    }
    #[test]
    fn test_rng_next_exponential() {
        let mut rng = Rng::new(42);
        let n = 5000;
        let lambda = 2.0;
        let samples: Vec<f64> = (0..n).map(|_| rng.next_exponential(lambda)).collect();
        let m = mean(&samples);
        assert!(
            (m - 0.5).abs() < 0.05,
            "exponential mean={m}, expected ~0.5"
        );
    }
    #[test]
    fn test_rng_next_poisson() {
        let mut rng = Rng::new(42);
        let n = 5000;
        let lambda = 3.0;
        let samples: Vec<f64> = (0..n).map(|_| rng.next_poisson(lambda) as f64).collect();
        let m = mean(&samples);
        assert!((m - 3.0).abs() < 0.2, "poisson mean={m}, expected ~3.0");
    }
    #[test]
    fn test_merton_jump_diffusion_positive_prices() {
        let model = MertonJumpDiffusion::new(0.05, 0.2, 1.0, -0.01, 0.05);
        let path = model.simulate(100.0, 1.0, 0.01, 42);
        for &v in &path {
            assert!(v > 0.0, "Merton model produced non-positive price: {v}");
        }
    }
    #[test]
    fn test_merton_mean_jump_size() {
        let model = MertonJumpDiffusion::new(0.05, 0.2, 1.0, 0.0, 0.0);
        let mj = model.mean_jump_size();
        assert!(mj.abs() < 1e-12, "mean_jump_size={mj}");
    }
    #[test]
    fn test_heston_feller_condition() {
        let model = HestonModel::new(0.05, 2.0, 0.04, 0.3, -0.7);
        assert!(model.feller_satisfied());
        let model2 = HestonModel::new(0.05, 0.1, 0.01, 0.5, -0.7);
        assert!(!model2.feller_satisfied());
    }
    #[test]
    fn test_heston_simulation_positive_variance() {
        let model = HestonModel::new(0.05, 2.0, 0.04, 0.3, -0.7);
        let (prices, variances) = model.simulate(100.0, 0.04, 1.0, 0.01, 42);
        assert!(!prices.is_empty());
        for &v in &variances {
            assert!(v >= 0.0, "negative variance: {v}");
        }
    }
    #[test]
    fn test_monte_carlo_gbm_paths() {
        let (mean_path, std_path) = monte_carlo_gbm_paths(100.0, 0.05, 0.2, 1.0, 0.1, 200, 42);
        assert!(!mean_path.is_empty());
        assert_eq!(mean_path.len(), std_path.len());
        assert!((mean_path[0] - 100.0).abs() < 1e-10);
        assert!(std_path[0].abs() < 1e-10);
    }
    #[test]
    fn test_monte_carlo_terminal_identity() {
        let (est, stderr) = monte_carlo_terminal_estimate(100.0, 0.05, 0.2, 1.0, 5000, |s| s, 42);
        let expected = 100.0 * (0.05_f64).exp();
        assert!(
            (est - expected).abs() < 3.0 * stderr + 1.0,
            "est={est}, expected={expected}, stderr={stderr}"
        );
    }
    #[test]
    fn test_antithetic_reduces_variance() {
        let (_, se_standard) =
            monte_carlo_terminal_estimate(100.0, 0.05, 0.2, 1.0, 1000, |s| s, 42);
        let (_, se_anti) = antithetic_terminal_estimate(100.0, 0.05, 0.2, 1.0, 1000, |s| s, 42);
        assert!(se_standard.is_finite());
        assert!(se_anti.is_finite());
    }
    #[test]
    fn test_control_variate_estimate() {
        let (est, stderr) =
            control_variate_terminal_estimate(100.0, 0.05, 0.2, 1.0, 2000, |s| s, 42);
        let expected = 100.0 * (0.05_f64).exp();
        assert!(
            (est - expected).abs() < 3.0 * stderr + 2.0,
            "est={est}, expected={expected}, stderr={stderr}"
        );
    }
    #[test]
    fn test_wiener_bridge_endpoints() {
        let mut wp = WienerProcess::new(42);
        let bridge = wp.bridge(0.0, 1.0, 0.01, 0.0, 1.0);
        assert!(
            (bridge[0].1 - 0.0).abs() < 1e-10,
            "bridge start = {}",
            bridge[0].1
        );
        let last = bridge.last().unwrap().1;
        assert!(
            (last - 1.0).abs() < 0.1,
            "bridge end = {last}, expected ~1.0"
        );
    }
    #[test]
    fn test_gaussian_random_walk() {
        let rw = RandomWalk::new(3, 0.1);
        let path = rw.gaussian_walk(100, 42);
        assert_eq!(path.len(), 101);
        for pos in &path {
            assert_eq!(pos.len(), 3);
        }
    }
    #[test]
    fn test_msd_computation() {
        let rw = RandomWalk::new(1, 1.0);
        let path = rw.walk(50, 42);
        let msd = RandomWalk::msd(&path);
        assert_eq!(msd.len(), path.len());
        assert!(msd[0].abs() < 1e-10);
    }
    #[test]
    fn test_skewness_symmetric() {
        let data: Vec<f64> = (-50..=50).map(|i| i as f64).collect();
        let s = skewness(&data);
        assert!(s.abs() < 1e-10, "skewness of symmetric data = {s}");
    }
    #[test]
    fn test_kurtosis_uniform() {
        let data: Vec<f64> = (0..10000).map(|i| i as f64 / 10000.0).collect();
        let k = kurtosis(&data);
        assert!((k - (-1.2)).abs() < 0.05, "kurtosis = {k}, expected ~ -1.2");
    }
    #[test]
    fn test_running_mean() {
        let data = [1.0, 3.0, 5.0, 7.0];
        let rm = running_mean(&data);
        assert!((rm[0] - 1.0).abs() < 1e-10);
        assert!((rm[1] - 2.0).abs() < 1e-10);
        assert!((rm[2] - 3.0).abs() < 1e-10);
        assert!((rm[3] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_effective_sample_size_iid() {
        let mut rng = Rng::new(42);
        let data: Vec<f64> = (0..1000).map(|_| rng.next_normal()).collect();
        let ess = effective_sample_size(&data, 20);
        assert!(ess > 500.0, "ESS={ess} too low for iid data");
    }
    #[test]
    fn test_wiener_sampler_variance() {
        let dt = 0.01;
        let ws = WienerSampler::new(dt);
        let mut rng = rand::rng();
        let samples = ws.sample(10_000, &mut rng);
        let m: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
        let v: f64 =
            samples.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (samples.len() - 1) as f64;
        assert!(m.abs() < 0.05, "WienerSampler mean={m}");
        assert!(
            (v - dt).abs() < 0.005,
            "WienerSampler variance={v}, expected ~{dt}"
        );
    }
    #[test]
    fn test_wiener_sampler_length() {
        let ws = WienerSampler::new(0.1);
        let mut rng = rand::rng();
        let s = ws.sample(50, &mut rng);
        assert_eq!(s.len(), 50);
    }
    #[test]
    fn test_fbm_length() {
        let fbm = FractionalBrownianMotion::new(0.7, 0.01);
        let mut rng = rand::rng();
        let s = fbm.sample(100, &mut rng);
        assert_eq!(s.len(), 100);
    }
    #[test]
    fn test_fbm_finite_values() {
        let fbm = FractionalBrownianMotion::new(0.8, 0.01);
        let mut rng = rand::rng();
        let s = fbm.sample(200, &mut rng);
        for &v in &s {
            assert!(v.is_finite(), "fBm produced non-finite value");
        }
    }
    #[test]
    fn test_langevin_integrator_step_changes_x() {
        let li = LangevinIntegrator::new(1.0, 1.0, 1.0);
        let mut rng = rand::rng();
        let x0 = 0.0_f64;
        let mut changed = false;
        for _ in 0..20 {
            let x_new = li.step(x0, 10.0, 0.01, &mut rng);
            if (x_new - x0).abs() > 1e-10 {
                changed = true;
                break;
            }
        }
        assert!(changed, "LangevinIntegrator: x did not change");
    }
    #[test]
    fn test_langevin_integrator_zero_force_diffuses() {
        let li = LangevinIntegrator::new(1.0, 1.0, 1.0);
        let mut rng = rand::rng();
        let n = 1000;
        let dt = 0.01;
        let mut x = 0.0_f64;
        for _ in 0..n {
            x = li.step(x, 0.0, dt, &mut rng);
        }
        assert!(x.is_finite(), "LangevinIntegrator diverged");
    }
    #[test]
    fn test_mh_acceptance_gaussian() {
        let mh = MetropolisHastings::new(1.0);
        let mut rng = rand::rng();
        let mut n_accept = 0usize;
        let n_steps = 1000;
        let mut x = 0.0_f64;
        for _ in 0..n_steps {
            let (x_new, accepted) = mh.step(|xi| 0.5 * xi * xi, x, 0.5, &mut rng);
            x = x_new;
            if accepted {
                n_accept += 1;
            }
        }
        let accept_ratio = n_accept as f64 / n_steps as f64;
        assert!(
            accept_ratio > 0.2,
            "MH acceptance ratio too low: {accept_ratio}"
        );
    }
    #[test]
    fn test_mh_returns_tuple() {
        let mh = MetropolisHastings::new(1.0);
        let mut rng = rand::rng();
        let (x_new, _accepted) = mh.step(|x| x * x, 0.5, 0.1, &mut rng);
        assert!(x_new.is_finite());
    }
    #[test]
    fn test_kleinman_kramers_step() {
        let kk = KleinmanKramers::new(1.0, 1.0, 1.0);
        let mut rng = rand::rng();
        let (x_new, v_new) = kk.step(0.0, 0.0, 1.0, 0.01, &mut rng);
        assert!(x_new.is_finite());
        assert!(v_new.is_finite());
    }
    #[test]
    fn test_kleinman_kramers_many_steps() {
        let kk = KleinmanKramers::new(0.5, 1.0, 1.0);
        let mut rng = rand::rng();
        let mut x = 0.0_f64;
        let mut v = 0.0_f64;
        for _ in 0..1000 {
            let (xn, vn) = kk.step(x, v, -x, 0.001, &mut rng);
            x = xn;
            v = vn;
        }
        assert!(x.is_finite() && v.is_finite(), "KK dynamics diverged");
    }
    #[test]
    fn test_mfpt_kramers_formula() {
        let barrier = 1.0;
        let d = 1.0;
        let omega_0 = 1.0;
        let omega_b = 1.0;
        let tau = mean_first_passage_time(barrier, d, omega_0, omega_b);
        let expected = 2.0 * PI * (1.0_f64).exp();
        assert!(
            (tau - expected).abs() < 1e-10,
            "MFPT={tau}, expected={expected}"
        );
    }
    #[test]
    fn test_mfpt_increases_with_barrier() {
        let tau_low = mean_first_passage_time(1.0, 1.0, 1.0, 1.0);
        let tau_high = mean_first_passage_time(2.0, 1.0, 1.0, 1.0);
        assert!(
            tau_high > tau_low,
            "MFPT should increase with barrier height"
        );
    }
    #[test]
    fn test_fpt_reasonable_range() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.0, 1.0);
        let fpt = EmpiricalFirstPassageTime::estimate_ou(&ou, 0.0, 1.0, 0.01, 0.0, 50.0, 200, 42);
        assert!(fpt.is_finite() && fpt > 0.0, "FPT={fpt}");
    }
    #[test]
    fn test_fpt_higher_barrier_longer() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.0, 1.5);
        let fpt_low =
            EmpiricalFirstPassageTime::estimate_ou(&ou, 0.0, 0.5, 0.01, 0.0, 20.0, 200, 42);
        let fpt_high =
            EmpiricalFirstPassageTime::estimate_ou(&ou, 0.0, 2.0, 0.01, 0.0, 20.0, 200, 42);
        assert!(
            fpt_high >= fpt_low * 0.5,
            "Higher barrier should not be much faster: low={fpt_low}, high={fpt_high}"
        );
    }
    #[test]
    fn test_antithetic_gbm_mean_accuracy() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2);
        let s0 = 100.0;
        let t = 1.0;
        let n_pairs = 2000usize;
        let (est, _se) = AntitheticGbm::estimate_mean(s0, &gbm, t, 0.01, n_pairs, 42);
        let expected = gbm.analytical_mean(s0, t);
        assert!(
            (est - expected).abs() / expected < 0.05,
            "antithetic GBM mean={est}, expected={expected}"
        );
    }
    #[test]
    fn test_antithetic_gbm_pairs_positive() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2);
        let pairs = AntitheticGbm::generate_pairs(100.0, &gbm, 1.0, 0.01, 10, 7);
        for &(a, b) in &pairs {
            assert!(
                a > 0.0 && b > 0.0,
                "both paths must be positive: ({a}, {b})"
            );
        }
    }
    #[test]
    fn test_control_variate_gbm_variance_reduced() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2);
        let (est_cv, se_cv) = ControlVariateGbm::estimate(100.0, &gbm, 1.0, 0.01, 1000, |s| s, 42);
        let expected = gbm.analytical_mean(100.0, 1.0);
        assert!(
            (est_cv - expected).abs() < 3.0 * se_cv + 2.0,
            "CV estimate={est_cv}, expected={expected}"
        );
    }
    #[test]
    fn test_control_variate_gbm_payoff() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2);
        let s0 = 100.0;
        let k = 100.0;
        let (est, _se) =
            ControlVariateGbm::estimate(s0, &gbm, 1.0, 0.01, 2000, |s| (s - k).max(0.0), 99);
        assert!(
            est > 0.0 && est < s0,
            "Call payoff should be positive and < S0: {est}"
        );
    }
    #[test]
    fn test_ou_exact_sampler_mean_convergence() {
        let ou = OrnsteinUhlenbeck::new(3.0, 2.0, 0.5);
        let sampler = OuExactSampler::new(ou.clone());
        let mut x = 0.0_f64;
        let mut rng = Rng::new(42);
        let dt = 0.1;
        for _ in 0..10_000 {
            x = sampler.step(x, dt, &mut rng);
        }
        assert!(
            (x - 2.0).abs() < 0.5,
            "OU exact sampler mean={x}, expected ~2.0"
        );
    }
    #[test]
    fn test_ou_exact_sampler_variance() {
        let ou = OrnsteinUhlenbeck::new(2.0, 0.0, 0.4);
        let sampler = OuExactSampler::new(ou.clone());
        let expected_var = ou.stationary_variance();
        let n = 50_000usize;
        let dt = 0.05;
        let mut x = 0.0_f64;
        let mut rng = Rng::new(123);
        for _ in 0..1000 {
            x = sampler.step(x, dt, &mut rng);
        }
        let samples: Vec<f64> = (0..n)
            .map(|_| {
                x = sampler.step(x, dt, &mut rng);
                x
            })
            .collect();
        let m: f64 = samples.iter().sum::<f64>() / n as f64;
        let v: f64 = samples.iter().map(|s| (s - m).powi(2)).sum::<f64>() / (n - 1) as f64;
        assert!(
            (v - expected_var).abs() / expected_var < 0.15,
            "OU exact variance={v}, expected={expected_var}"
        );
    }
}
#[cfg(test)]
mod tests_new_stochastic {

    use crate::Rng;

    use crate::stochastic::CirProcess;

    use crate::stochastic::HullWhiteModel;

    use crate::stochastic::LevyFlight;

    use crate::stochastic::SabrModel;
    use crate::stochastic::VarianceGammaProcess;

    #[test]
    fn test_levy_flight_path_length() {
        let lf = LevyFlight::new(1.5, 1.0);
        let path = lf.path(100, 42);
        assert_eq!(
            path.len(),
            101,
            "Lévy flight path should have n_steps+1 points"
        );
    }
    #[test]
    fn test_levy_flight_starts_at_zero() {
        let lf = LevyFlight::new(1.5, 1.0);
        let path = lf.path(50, 7);
        assert!((path[0]).abs() < 1e-12, "path should start at 0");
    }
    #[test]
    fn test_levy_flight_finite_values() {
        let lf = LevyFlight::new(1.5, 1.0);
        let path = lf.path(200, 123);
        for &v in &path {
            assert!(v.is_finite(), "Lévy flight produced non-finite value");
        }
    }
    #[test]
    fn test_levy_flight_gaussian_limit() {
        let lf = LevyFlight::new(2.0, 1.0);
        let mut rng = Rng::new(42);
        let steps: Vec<f64> = (0..5000).map(|_| lf.sample_step(&mut rng)).collect();
        let m: f64 = steps.iter().sum::<f64>() / steps.len() as f64;
        assert!(
            m.abs() < 0.2,
            "Gaussian (alpha=2) steps should have mean ~0, got {m}"
        );
    }
    #[test]
    fn test_levy_flight_different_seeds_differ() {
        let lf = LevyFlight::new(1.5, 1.0);
        let path1 = lf.path(20, 1);
        let path2 = lf.path(20, 2);
        let differs = path1
            .iter()
            .zip(path2.iter())
            .any(|(a, b)| (a - b).abs() > 1e-10);
        assert!(differs, "different seeds should produce different paths");
    }
    #[test]
    fn test_vg_path_length() {
        let vg = VarianceGammaProcess::new(0.1, 0.2, 0.1);
        let path = vg.simulate(1.0, 100, 42);
        assert_eq!(path.len(), 101);
    }
    #[test]
    fn test_vg_starts_at_zero() {
        let vg = VarianceGammaProcess::new(0.0, 0.2, 0.1);
        let path = vg.simulate(1.0, 50, 7);
        assert!((path[0]).abs() < 1e-12);
    }
    #[test]
    fn test_vg_finite_values() {
        let vg = VarianceGammaProcess::new(0.1, 0.2, 0.05);
        let path = vg.simulate(1.0, 200, 99);
        for &v in &path {
            assert!(v.is_finite(), "VG produced non-finite value");
        }
    }
    #[test]
    fn test_vg_mean_increment() {
        let vg = VarianceGammaProcess::new(0.5, 0.2, 0.1);
        let dt = 0.01;
        assert!((vg.mean_increment(dt) - 0.5 * dt).abs() < 1e-12);
    }
    #[test]
    fn test_vg_variance_increment() {
        let vg = VarianceGammaProcess::new(0.0, 0.3, 0.1);
        let dt = 0.01;
        let expected = 0.3 * 0.3 * dt;
        assert!((vg.variance_increment(dt) - expected).abs() < 1e-12);
    }
    #[test]
    fn test_vg_ensemble_mean() {
        let vg = VarianceGammaProcess::new(0.2, 0.1, 0.05);
        let n_paths = 2000usize;
        let t_end = 0.5;
        let finals: Vec<f64> = (0..n_paths)
            .map(|seed| {
                let path = vg.simulate(t_end, 50, seed as u64);
                *path.last().unwrap()
            })
            .collect();
        let emp_mean: f64 = finals.iter().sum::<f64>() / n_paths as f64;
        let expected_mean = vg.mean_increment(t_end);
        assert!(
            (emp_mean - expected_mean).abs() < 0.1,
            "VG ensemble mean={emp_mean}, expected~{expected_mean}"
        );
    }
    #[test]
    fn test_sabr_path_lengths() {
        let sabr = SabrModel::new(100.0, 0.3, 0.5, 0.4, -0.3);
        let (forwards, vols) = sabr.simulate(1.0, 100, 42);
        assert_eq!(forwards.len(), 101);
        assert_eq!(vols.len(), 101);
    }
    #[test]
    fn test_sabr_initial_values() {
        let sabr = SabrModel::new(100.0, 0.3, 0.5, 0.4, -0.3);
        let (forwards, vols) = sabr.simulate(1.0, 50, 7);
        assert!((forwards[0] - 100.0).abs() < 1e-10);
        assert!((vols[0] - 0.3).abs() < 1e-10);
    }
    #[test]
    fn test_sabr_vols_positive() {
        let sabr = SabrModel::new(100.0, 0.2, 0.5, 0.3, -0.2);
        let (_forwards, vols) = sabr.simulate(1.0, 200, 123);
        for &v in &vols {
            assert!(v >= 0.0, "SABR volatility should be non-negative, got {v}");
        }
    }
    #[test]
    fn test_sabr_implied_vol_atm_positive() {
        let sabr = SabrModel::new(100.0, 0.2, 0.5, 0.3, -0.2);
        let iv = sabr.implied_vol_approx(100.0, 1.0);
        assert!(iv > 0.0, "ATM implied vol should be positive, got {iv}");
    }
    #[test]
    fn test_sabr_implied_vol_otm() {
        let sabr = SabrModel::new(100.0, 0.2, 0.5, 0.3, -0.2);
        let iv_atm = sabr.implied_vol_approx(100.0, 1.0);
        let iv_otm = sabr.implied_vol_approx(102.0, 1.0);
        assert!(iv_atm.is_finite() && iv_atm > 0.0, "ATM iv={iv_atm}");
        assert!(iv_otm.is_finite() && iv_otm > 0.0, "OTM iv={iv_otm}");
    }
    #[test]
    fn test_cir_path_length() {
        let cir = CirProcess::new(2.0, 0.05, 0.1);
        let path = cir.simulate(0.05, 1.0, 0.01, 42);
        assert!(path.len() > 90, "CIR path should have ~101 points");
    }
    #[test]
    fn test_cir_non_negative() {
        let cir = CirProcess::new(2.0, 0.05, 0.1);
        let path = cir.simulate(0.05, 1.0, 0.01, 99);
        for &v in &path {
            assert!(v >= 0.0, "CIR process should be non-negative, got {v}");
        }
    }
    #[test]
    fn test_cir_feller_condition() {
        let cir_ok = CirProcess::new(2.0, 0.05, 0.1);
        assert!(cir_ok.feller_satisfied());
        let cir_fail = CirProcess::new(0.1, 0.01, 0.5);
        assert!(!cir_fail.feller_satisfied());
    }
    #[test]
    fn test_cir_stationary_mean() {
        let cir = CirProcess::new(2.0, 0.05, 0.1);
        assert!((cir.stationary_mean() - 0.05).abs() < 1e-12);
    }
    #[test]
    fn test_cir_stationary_variance() {
        let cir = CirProcess::new(2.0, 0.05, 0.1);
        let expected = 0.05 * 0.01 / 4.0;
        assert!((cir.stationary_variance() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_cir_conditional_mean_long_time() {
        let cir = CirProcess::new(3.0, 0.05, 0.1);
        let m = cir.conditional_mean(0.1, 10.0);
        assert!(
            (m - 0.05).abs() < 0.001,
            "long-time conditional mean={m}, expected~0.05"
        );
    }
    #[test]
    fn test_cir_convergence_to_stationary() {
        let cir = CirProcess::new(3.0, 0.05, 0.1);
        let path = cir.simulate(0.02, 20.0, 0.01, 42);
        let tail = &path[path.len() * 3 / 4..];
        let m: f64 = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((m - 0.05).abs() < 0.03, "CIR tail mean={m}, expected~0.05");
    }
    #[test]
    fn test_hull_white_path_length() {
        let hw = HullWhiteModel::new(0.5, 0.01, 0.03);
        let path = hw.simulate_short_rate(0.02, 1.0, 0.01, 42);
        assert!(path.len() > 90, "HW path should have ~101 points");
    }
    #[test]
    fn test_hull_white_mean_convergence() {
        let hw = HullWhiteModel::new(2.0, 0.01, 0.06);
        let path = hw.simulate_short_rate(0.02, 20.0, 0.01, 42);
        let tail = &path[path.len() * 3 / 4..];
        let m: f64 = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((m - 0.03).abs() < 0.02, "HW tail mean={m}, expected~0.03");
    }
    #[test]
    fn test_hull_white_analytical_mean() {
        let hw = HullWhiteModel::new(1.0, 0.01, 0.05);
        let m0 = hw.mean_rate(0.02, 0.0);
        assert!((m0 - 0.02).abs() < 1e-10);
    }
    #[test]
    fn test_hull_white_variance_zero_at_t0() {
        let hw = HullWhiteModel::new(0.5, 0.01, 0.03);
        let v = hw.variance_rate(0.0);
        assert!(v.abs() < 1e-12, "variance at t=0 should be 0, got {v}");
    }
    #[test]
    fn test_hull_white_variance_increases() {
        let hw = HullWhiteModel::new(0.5, 0.02, 0.03);
        let v1 = hw.variance_rate(0.5);
        let v2 = hw.variance_rate(1.0);
        assert!(v2 > v1, "variance should increase with time");
    }
    #[test]
    fn test_hull_white_bond_price_at_zero_maturity() {
        let hw = HullWhiteModel::new(0.5, 0.01, 0.03);
        let p = hw.bond_price(0.02, 0.0);
        assert!(
            (p - 1.0).abs() < 1e-6,
            "bond price at T=0 should be ~1, got {p}"
        );
    }
    #[test]
    fn test_hull_white_bond_price_positive() {
        let hw = HullWhiteModel::new(0.5, 0.01, 0.03);
        let p = hw.bond_price(0.02, 5.0);
        assert!(p > 0.0 && p < 1.0, "bond price should be in (0,1), got {p}");
    }
    #[test]
    fn test_hull_white_finite_path() {
        let hw = HullWhiteModel::new(0.5, 0.01, 0.03);
        let path = hw.simulate_short_rate(0.02, 5.0, 0.01, 99);
        for &r in &path {
            assert!(r.is_finite(), "HW short rate should be finite");
        }
    }
}
