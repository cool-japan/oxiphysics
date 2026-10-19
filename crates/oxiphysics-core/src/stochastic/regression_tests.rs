//! Regression tests for the stochastic-process fixes: the Langevin thermostat
//! noise scale, the SABR smile factor and the Hosking fractional-Gaussian-noise
//! sampler.

use super::functions::K_B;
use super::types::{FractionalBrownianMotion, LangevinDynamics, Rng, SabrModel};
use rand::SeedableRng;
use rand::rngs::StdRng;

// ---------------------------------------------------------------------------
// LangevinDynamics::step
// ---------------------------------------------------------------------------

/// Stationary `m <v²> / (k_B T)` of a force-free particle after a burn-in.
fn equipartition_ratio(thermostat: &LangevinDynamics, dt: f64, seed: u64) -> f64 {
    let mut rng = Rng::new(seed);
    let mut vel = 0.0_f64;
    for _ in 0..10_000 {
        vel = thermostat.step(0.0, vel, 0.0, dt, &mut rng).1;
    }
    let steps = 1_000_000_u32;
    let mut sum_v2 = 0.0_f64;
    for _ in 0..steps {
        vel = thermostat.step(0.0, vel, 0.0, dt, &mut rng).1;
        sum_v2 += vel * vel;
    }
    thermostat.mass * (sum_v2 / f64::from(steps)) / (K_B * thermostat.temperature)
}

#[test]
fn langevin_step_thermalises_to_the_equipartition_temperature() {
    // gamma * dt = 0.01 and 0.05: the Euler-Maruyama scheme has the stationary
    // variance (k_B T / m) / (1 - gamma dt / 2), i.e. +0.5 % and +2.6 %.
    let thermostat = LangevinDynamics::new(1.0e-26, 1.0e10, 300.0);
    for dt in [1.0e-12, 5.0e-12] {
        let ratio = equipartition_ratio(&thermostat, dt, 20_240_611);
        let predicted = 1.0 / (1.0 - 0.5 * thermostat.friction * dt);
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "dt = {dt}: m<v^2>/(k_B T) = {ratio}, expected about 1"
        );
        assert!(
            (ratio - predicted).abs() < 0.05,
            "dt = {dt}: m<v^2>/(k_B T) = {ratio}, scheme predicts {predicted}"
        );
    }
}

#[test]
fn langevin_random_kick_has_the_euler_maruyama_scale() {
    // From rest and without force the new velocity is exactly the random kick,
    // whose standard deviation is sqrt(2 gamma k_B T dt / m): it must not be
    // divided by sqrt(dt) once more.
    let thermostat = LangevinDynamics::new(1.0e-26, 1.0e10, 300.0);
    for dt in [1.0e-15, 1.0e-12, 5.0e-12] {
        let mut rng = Rng::new(99);
        let samples = 40_000_u32;
        let mut sum_sq = 0.0_f64;
        for _ in 0..samples {
            let kick = thermostat.step(0.0, 0.0, 0.0, dt, &mut rng).1;
            sum_sq += kick * kick;
        }
        let rms = (sum_sq / f64::from(samples)).sqrt();
        let expected = (2.0 * thermostat.friction * K_B * thermostat.temperature / thermostat.mass
            * dt)
            .sqrt();
        assert!(
            (rms / expected - 1.0).abs() < 0.03,
            "dt = {dt}: kick rms {rms}, expected {expected}"
        );
    }
}

// ---------------------------------------------------------------------------
// SabrModel::implied_vol_approx
// ---------------------------------------------------------------------------

/// Reference implied volatilities of Hagan et al.'s formula (2.17a), evaluated
/// in 40-digit arithmetic.
fn assert_smile(sabr: &SabrModel, t: f64, reference: &[(f64, f64)]) {
    for &(strike, expected) in reference {
        let iv = sabr.implied_vol_approx(strike, t);
        assert!(
            (iv - expected).abs() <= 1e-11 * expected.abs(),
            "K = {strike}: implied vol {iv}, reference {expected}"
        );
    }
}

#[test]
fn sabr_smile_matches_the_hagan_formula() {
    let rates = SabrModel::new(0.04, 0.2, 0.5, 0.4, -0.3);
    assert_smile(
        &rates,
        1.0,
        &[
            (0.04, 1.00695),
            (0.0399, 1.00773529863871),
            (0.0401, 1.00616723008142),
            (0.039, 1.01491899463123),
            (0.041, 0.999233960637095),
            (0.03, 1.10072140328405),
            (0.06, 0.887947411358506),
        ],
    );
    let equity = SabrModel::new(100.0, 0.2, 0.5, 0.3, -0.2);
    assert_smile(
        &equity,
        1.0,
        &[
            (100.0, 0.0201380833333333),
            (102.0, 0.0197347669099528),
            (98.0, 0.0211121242317547),
            (90.0, 0.0280228245031828),
            (120.0, 0.0285592064437594),
        ],
    );
    let lognormal = SabrModel::new(1.0, 0.3, 1.0, 0.5, 0.4);
    assert_smile(
        &lognormal,
        2.0,
        &[
            (0.8, 0.301404230190628),
            (1.0, 0.3185),
            (1.25, 0.346707817780297),
        ],
    );
}

#[test]
fn sabr_smile_is_continuous_across_the_at_the_money_branch() {
    // The smile slope here is about -7.9 per unit strike: strikes closer to the
    // forward than 1e-8 take the ATM branch, all others the general branch, and
    // the two must join without a jump.
    let sabr = SabrModel::new(0.04, 0.2, 0.5, 0.4, -0.3);
    let atm = sabr.implied_vol_approx(0.04, 1.0);
    for offset in [1.0e-4, 1.0e-6, 4.0e-8, 2.0e-8] {
        for strike in [0.04 - offset, 0.04 + offset] {
            let iv = sabr.implied_vol_approx(strike, 1.0);
            assert!(
                (iv - atm).abs() < 10.0 * offset,
                "K = {strike}: implied vol {iv} jumps away from the ATM value {atm}"
            );
        }
    }
}

#[test]
fn sabr_negative_correlation_gives_a_downward_sloping_skew() {
    let sabr = SabrModel::new(0.04, 0.2, 0.5, 0.4, -0.3);
    let low = sabr.implied_vol_approx(0.039, 1.0);
    let high = sabr.implied_vol_approx(0.041, 1.0);
    assert!(
        low > high,
        "rho < 0 needs iv(0.039) = {low} > iv(0.041) = {high}"
    );
    let positive = SabrModel::new(0.04, 0.2, 0.5, 0.4, 0.3);
    let low = positive.implied_vol_approx(0.039, 1.0);
    let high = positive.implied_vol_approx(0.041, 1.0);
    assert!(
        low < high,
        "rho > 0 needs iv(0.039) = {low} < iv(0.041) = {high}"
    );
}

#[test]
fn sabr_out_of_the_money_volatility_stays_close_to_at_the_money() {
    let sabr = SabrModel::new(100.0, 0.2, 0.5, 0.3, -0.2);
    let atm = sabr.implied_vol_approx(100.0, 1.0);
    let otm = sabr.implied_vol_approx(102.0, 1.0);
    assert!(
        (otm / atm - 1.0).abs() < 0.1,
        "OTM volatility {otm} is not within 10 % of ATM {atm}"
    );
}

// ---------------------------------------------------------------------------
// FractionalBrownianMotion::sample (Hosking / Durbin-Levinson)
// ---------------------------------------------------------------------------

/// Autocovariance of fractional Gaussian noise at an integer lag.
fn fgn_autocovariance(h: f64, dt: f64, lag: usize) -> f64 {
    let k = lag as f64;
    let two_h = 2.0 * h;
    0.5 * ((k + 1.0).powf(two_h) - 2.0 * k.powf(two_h) + (k - 1.0).abs().powf(two_h))
        * dt.powf(two_h)
}

/// `L z` for the lower Cholesky factor `L` of the Toeplitz autocovariance matrix.
fn cholesky_colouring(h: f64, dt: f64, z: &[f64]) -> Vec<f64> {
    let n = z.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = fgn_autocovariance(h, dt, i - j);
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            l[i][j] = if i == j { sum.sqrt() } else { sum / l[j][j] };
        }
    }
    (0..n)
        .map(|i| (0..=i).map(|j| l[i][j] * z[j]).sum())
        .collect()
}

#[test]
fn hosking_sampler_reproduces_the_cholesky_colouring_of_the_same_normals() {
    // Durbin-Levinson is the Cholesky factorisation of the Toeplitz covariance
    // in disguise: both must map the same normals to the same sample.
    for h in [0.2, 0.3, 0.5, 0.7, 0.8, 0.95] {
        let dt = 0.37;
        let fbm = FractionalBrownianMotion::new(h, dt);
        let mut rng = Rng::new(2024);
        let z: Vec<f64> = (0..32).map(|_| rng.next_normal()).collect();
        let expected = cholesky_colouring(h, dt, &z);
        let got = fbm.correlate_normals(&z);
        assert_eq!(got.len(), expected.len());
        for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
            assert!(
                (g - e).abs() <= 1e-9 * (1.0 + e.abs()),
                "H = {h}, i = {i}: sampler {g}, Cholesky {e}"
            );
        }
    }
}

#[test]
fn fbm_increments_have_the_fractional_gaussian_noise_autocorrelation() {
    // Pooled lag-k autocorrelation of 2000 independent paths of 64 increments;
    // the sampling error of each estimate is below 0.006.
    let dt = 0.01;
    let paths = 2000;
    let n = 64;
    for (h, seed) in [(0.8, 11_u64), (0.5, 12), (0.3, 13)] {
        let fbm = FractionalBrownianMotion::new(h, dt);
        let mut rng = StdRng::seed_from_u64(seed);
        let mut lag_sums = [0.0_f64; 6];
        for _ in 0..paths {
            let x = fbm.sample(n, &mut rng);
            for (lag, sum) in lag_sums.iter_mut().enumerate() {
                *sum += x[..n - lag]
                    .iter()
                    .zip(&x[lag..])
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
            }
        }
        let gamma_hat = |lag: usize| lag_sums[lag] / (paths * (n - lag)) as f64;
        for lag in 1..=5 {
            let measured = gamma_hat(lag) / gamma_hat(0);
            let expected = fgn_autocovariance(h, dt, lag) / fgn_autocovariance(h, dt, 0);
            assert!(
                (measured - expected).abs() < 0.03,
                "H = {h}, lag {lag}: autocorrelation {measured}, expected {expected}"
            );
        }
        let variance_ratio = gamma_hat(0) / fgn_autocovariance(h, dt, 0);
        assert!(
            (variance_ratio - 1.0).abs() < 0.05,
            "H = {h}: increment variance off by a factor {variance_ratio}"
        );
    }
}

#[test]
fn fbm_partial_sums_follow_the_power_law_variance() {
    // Var[B(m dt)] = (m dt)^(2H) for fractional Brownian motion.
    let dt = 0.05;
    let paths = 4000;
    let n = 64;
    for (h, seed) in [(0.8, 21_u64), (0.3, 22)] {
        let fbm = FractionalBrownianMotion::new(h, dt);
        let mut rng = StdRng::seed_from_u64(seed);
        let checkpoints = [1_usize, 8, 32, 64];
        let mut second_moments = [0.0_f64; 4];
        for _ in 0..paths {
            let x = fbm.sample(n, &mut rng);
            for (moment, &m) in second_moments.iter_mut().zip(&checkpoints) {
                let b: f64 = x[..m].iter().sum();
                *moment += b * b;
            }
        }
        for (moment, &m) in second_moments.iter().zip(&checkpoints) {
            let measured = moment / f64::from(paths);
            let expected = (m as f64 * dt).powf(2.0 * h);
            assert!(
                (measured / expected - 1.0).abs() < 0.1,
                "H = {h}, m = {m}: Var[B] = {measured}, expected {expected}"
            );
        }
    }
}

#[test]
fn fbm_sample_handles_trivial_and_extreme_inputs() {
    let mut rng = StdRng::seed_from_u64(5);
    let fbm = FractionalBrownianMotion::new(0.7, 0.01);
    assert!(fbm.sample(0, &mut rng).is_empty());
    assert_eq!(fbm.sample(1, &mut rng).len(), 1);
    // Hurst exponents very close to 1 make the increments almost perfectly
    // correlated; the recursion must stay finite.
    let nearly_smooth = FractionalBrownianMotion::new(0.999, 0.01);
    let x = nearly_smooth.sample(512, &mut rng);
    assert_eq!(x.len(), 512);
    assert!(x.iter().all(|v| v.is_finite()));
    // H = 1/2 is white noise with variance dt: increments are independent.
    let white = FractionalBrownianMotion::new(0.5, 4.0);
    let mut state = Rng::new(31);
    let z: Vec<f64> = (0..16).map(|_| state.next_normal()).collect();
    for (value, normal) in white.correlate_normals(&z).iter().zip(&z) {
        assert!((value - 2.0 * normal).abs() < 1e-12);
    }
}
