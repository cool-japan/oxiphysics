//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::Xoshiro256;

/// Box-Muller transform: return two independent N(0,1) values.
pub fn box_muller(rng: &mut Xoshiro256) -> (f64, f64) {
    let u1 = rng.next_f64().max(1e-300);
    let u2 = rng.next_f64();
    let r = (-2.0 * u1.ln()).sqrt();
    let theta = 2.0 * PI * u2;
    (r * theta.cos(), r * theta.sin())
}
/// Sample from Gamma(shape, scale) using Marsaglia-Tsang method.
pub fn sample_gamma(rng: &mut Xoshiro256, shape: f64, scale: f64) -> f64 {
    if shape < 1.0 {
        let u = rng.next_f64().max(1e-300);
        return sample_gamma(rng, shape + 1.0, scale) * u.powf(1.0 / shape);
    }
    let d = shape - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    loop {
        let (z, _) = box_muller(rng);
        let v_raw = 1.0 + c * z;
        if v_raw <= 0.0 {
            continue;
        }
        let v = v_raw.powi(3);
        let u = rng.next_f64().max(1e-300);
        if u < 1.0 - 0.0331 * (z * z) * (z * z) {
            return d * v * scale;
        }
        if u.ln() < 0.5 * z * z + d * (1.0 - v + v.ln()) {
            return d * v * scale;
        }
    }
}
/// Choose an index from a probability vector given uniform variate `u`.
pub fn weighted_choice_index(probs: &[f64], u: f64) -> usize {
    let mut cumsum = 0.0;
    for (i, &p) in probs.iter().enumerate() {
        cumsum += p;
        if u < cumsum {
            return i;
        }
    }
    probs.len() - 1
}
/// Compute the 1D Halton sequence value at index `i` with base `base`.
pub fn halton_1d(mut i: usize, base: usize) -> f64 {
    let mut f = 1.0f64;
    let mut r = 0.0f64;
    while i > 0 {
        f /= base as f64;
        r += f * (i % base) as f64;
        i /= base;
    }
    r
}
/// Return the first `n` prime numbers.
pub fn first_n_primes(n: usize) -> Vec<usize> {
    let mut primes = Vec::with_capacity(n);
    let mut candidate = 2usize;
    while primes.len() < n {
        if primes.iter().all(|&p| !candidate.is_multiple_of(p)) {
            primes.push(candidate);
        }
        candidate += 1;
    }
    primes
}
/// Fisher-Yates shuffle of a mutable slice in-place.
pub fn shuffle_fisher_yates<T>(slice: &mut [T], rng: &mut Xoshiro256) {
    let n = slice.len();
    for i in (1..n).rev() {
        let j = (rng.next_u64() as usize) % (i + 1);
        slice.swap(i, j);
    }
}
/// Sample `k` distinct indices from `0..n` without replacement.
pub fn sample_without_replacement(rng: &mut Xoshiro256, n: usize, k: usize) -> Vec<usize> {
    let k = k.min(n);
    let mut pool: Vec<usize> = (0..n).collect();
    shuffle_fisher_yates(&mut pool, rng);
    pool.truncate(k);
    pool
}
/// Weighted random choice from a slice of `(value, weight)` pairs.
pub fn weighted_choice<T: Clone>(rng: &mut Xoshiro256, choices: &[(T, f64)]) -> Option<T> {
    if choices.is_empty() {
        return None;
    }
    let total: f64 = choices.iter().map(|(_, w)| w).sum();
    let u = rng.next_f64() * total;
    let mut cumsum = 0.0;
    for (val, w) in choices {
        cumsum += w;
        if u <= cumsum {
            return Some(val.clone());
        }
    }
    choices.last().map(|c| c.0.clone())
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::random::DirichletDist;

    use crate::random::ExponentialDist;

    use crate::random::MultinomialDist;
    use crate::random::NormalDist;

    use crate::random::Pcg64;
    use crate::random::PoissonDist;
    use crate::random::RandomSampler;

    use crate::random::SplitMix64;

    use crate::random::UniformDist;

    use crate::random::Xoshiro256;

    #[test]
    fn splitmix64_deterministic() {
        let mut a = SplitMix64::new(42);
        let mut b = SplitMix64::new(42);
        assert_eq!(a.next_u64(), b.next_u64());
    }
    #[test]
    fn splitmix64_different_seeds() {
        let mut a = SplitMix64::new(1);
        let mut b = SplitMix64::new(2);
        assert_ne!(a.next_u64(), b.next_u64());
    }
    #[test]
    fn splitmix64_f64_range() {
        let mut rng = SplitMix64::new(123);
        for _ in 0..1000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v));
        }
    }
    #[test]
    fn pcg64_deterministic() {
        let mut a = Pcg64::new(99, 1);
        let mut b = Pcg64::new(99, 1);
        for _ in 0..10 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }
    #[test]
    fn pcg64_f64_range() {
        let mut rng = Pcg64::new(42, 1);
        for _ in 0..1000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v), "v={v}");
        }
    }
    #[test]
    fn pcg64_range() {
        let mut rng = Pcg64::new(7, 3);
        for _ in 0..1000 {
            let v = rng.next_f64_range(2.0, 5.0);
            assert!((2.0..5.0).contains(&v), "v={v}");
        }
    }
    #[test]
    fn pcg64_different_streams() {
        let mut a = Pcg64::new(0, 0);
        let mut b = Pcg64::new(0, 1);
        let va: Vec<u64> = (0..10).map(|_| a.next_u64()).collect();
        let vb: Vec<u64> = (0..10).map(|_| b.next_u64()).collect();
        assert_ne!(va, vb);
    }
    #[test]
    fn xoshiro256_deterministic() {
        let mut a = Xoshiro256::new(55);
        let mut b = Xoshiro256::new(55);
        for _ in 0..20 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }
    #[test]
    fn xoshiro256_f64_range() {
        let mut rng = Xoshiro256::new(1);
        for _ in 0..1000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v));
        }
    }
    #[test]
    fn xoshiro256_jump_changes_state() {
        let mut a = Xoshiro256::new(10);
        let mut b = Xoshiro256::new(10);
        b.jump();
        assert_ne!(a.next_u64(), b.next_u64());
    }
    #[test]
    fn xoshiro256_long_jump_changes_state() {
        let mut a = Xoshiro256::new(10);
        let mut b = Xoshiro256::new(10);
        b.long_jump();
        assert_ne!(a.next_u64(), b.next_u64());
    }
    #[test]
    fn uniform_dist_range() {
        let dist = UniformDist::new(-3.0, 7.0);
        let mut rng = Xoshiro256::new(1);
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!((-3.0..7.0).contains(&v), "v={v}");
        }
    }
    #[test]
    fn uniform_dist_sample_n_count() {
        let dist = UniformDist::new(0.0, 1.0);
        let mut rng = Xoshiro256::new(2);
        let samples = dist.sample_n(&mut rng, 50);
        assert_eq!(samples.len(), 50);
    }
    #[test]
    fn normal_dist_mean_approx() {
        let mut dist = NormalDist::new(5.0, 1.0);
        let mut rng = Xoshiro256::new(3);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!((mean - 5.0).abs() < 0.1, "mean={mean}");
    }
    #[test]
    fn normal_dist_std_approx() {
        let mut dist = NormalDist::new(0.0, 2.0);
        let mut rng = Xoshiro256::new(4);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
        assert!((var.sqrt() - 2.0).abs() < 0.1, "std={}", var.sqrt());
    }
    #[test]
    fn poisson_small_lambda_mean() {
        let dist = PoissonDist::new(3.0);
        let mut rng = Xoshiro256::new(5);
        let samples: Vec<u64> = (0..10000).map(|_| dist.sample(&mut rng)).collect();
        let mean = samples.iter().sum::<u64>() as f64 / samples.len() as f64;
        assert!((mean - 3.0).abs() < 0.2, "mean={mean}");
    }
    #[test]
    fn poisson_large_lambda_mean() {
        let dist = PoissonDist::new(100.0);
        let mut rng = Xoshiro256::new(6);
        let samples: Vec<u64> = (0..5000).map(|_| dist.sample(&mut rng)).collect();
        let mean = samples.iter().sum::<u64>() as f64 / samples.len() as f64;
        assert!((mean - 100.0).abs() < 2.0, "mean={mean}");
    }
    #[test]
    fn exponential_dist_mean() {
        let dist = ExponentialDist::new(2.0);
        let mut rng = Xoshiro256::new(7);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!((mean - 0.5).abs() < 0.05, "mean={mean}");
    }
    #[test]
    fn exponential_dist_all_positive() {
        let dist = ExponentialDist::new(1.0);
        let mut rng = Xoshiro256::new(8);
        let samples = dist.sample_n(&mut rng, 1000);
        assert!(samples.iter().all(|&v| v > 0.0));
    }
    #[test]
    fn dirichlet_sums_to_one() {
        let dist = DirichletDist::new(vec![1.0, 2.0, 3.0]);
        let mut rng = Xoshiro256::new(9);
        for _ in 0..100 {
            let s = dist.sample(&mut rng);
            let sum: f64 = s.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "sum={sum}");
        }
    }
    #[test]
    fn dirichlet_all_positive() {
        let dist = DirichletDist::new(vec![0.5, 0.5, 0.5]);
        let mut rng = Xoshiro256::new(10);
        for _ in 0..100 {
            let s = dist.sample(&mut rng);
            assert!(s.iter().all(|&v| v > 0.0));
        }
    }
    #[test]
    fn dirichlet_uniform_mean() {
        let dist = DirichletDist::new(vec![1.0, 1.0, 1.0]);
        let mut rng = Xoshiro256::new(11);
        let n = 5000usize;
        let samples = dist.sample_n(&mut rng, n);
        let mean0 = samples.iter().map(|s| s[0]).sum::<f64>() / n as f64;
        assert!((mean0 - 1.0 / 3.0).abs() < 0.02, "mean0={mean0}");
    }
    #[test]
    fn multinomial_counts_sum_to_n() {
        let dist = MultinomialDist::new(vec![0.3, 0.4, 0.3]);
        let mut rng = Xoshiro256::new(12);
        let counts = dist.sample(&mut rng, 100);
        assert_eq!(counts.iter().sum::<usize>(), 100);
    }
    #[test]
    fn multinomial_distribution_approx() {
        let probs = vec![0.5, 0.3, 0.2];
        let dist = MultinomialDist::new(probs.clone());
        let mut rng = Xoshiro256::new(13);
        let n = 10000usize;
        let counts = dist.sample(&mut rng, n);
        let freq: Vec<f64> = counts.iter().map(|&c| c as f64 / n as f64).collect();
        for (f, p) in freq.iter().zip(probs.iter()) {
            assert!((f - p).abs() < 0.02, "freq={f}, prob={p}");
        }
    }
    #[test]
    fn stratified_sampling_covers_interval() {
        let mut rng = Xoshiro256::new(14);
        let samples = RandomSampler::stratified_1d(&mut rng, 10);
        assert_eq!(samples.len(), 10);
        for (i, &s) in samples.iter().enumerate() {
            let lo = i as f64 / 10.0;
            let hi = (i + 1) as f64 / 10.0;
            assert!(s >= lo && s < hi, "sample[{i}]={s} outside [{lo},{hi})");
        }
    }
    #[test]
    fn latin_hypercube_dimensions() {
        let mut rng = Xoshiro256::new(15);
        let samples = RandomSampler::latin_hypercube(&mut rng, 20, 3);
        assert_eq!(samples.len(), 20);
        assert_eq!(samples[0].len(), 3);
    }
    #[test]
    fn latin_hypercube_range() {
        let mut rng = Xoshiro256::new(16);
        let samples = RandomSampler::latin_hypercube(&mut rng, 10, 2);
        for row in &samples {
            for &v in row {
                assert!((0.0..1.0).contains(&v), "v={v}");
            }
        }
    }
    #[test]
    fn halton_sequence_dimensions() {
        let samples = RandomSampler::halton(100, 3, 0);
        assert_eq!(samples.len(), 100);
        assert_eq!(samples[0].len(), 3);
    }
    #[test]
    fn halton_1d_base2_known_values() {
        assert!((halton_1d(1, 2) - 0.5).abs() < 1e-10);
        assert!((halton_1d(2, 2) - 0.25).abs() < 1e-10);
        assert!((halton_1d(3, 2) - 0.75).abs() < 1e-10);
    }
    #[test]
    fn sobol_dimensions() {
        let samples = RandomSampler::sobol(50, 4);
        assert_eq!(samples.len(), 50);
        assert_eq!(samples[0].len(), 4);
    }
    #[test]
    fn sobol_range() {
        let samples = RandomSampler::sobol(100, 2);
        for row in &samples {
            for &v in row {
                assert!((0.0..=1.0).contains(&v), "v={v}");
            }
        }
    }
    #[test]
    fn fisher_yates_shuffle_is_permutation() {
        let mut rng = Xoshiro256::new(20);
        let mut v: Vec<usize> = (0..20).collect();
        let original = v.clone();
        shuffle_fisher_yates(&mut v, &mut rng);
        let mut sorted = v.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, original);
    }
    #[test]
    fn sample_without_replacement_count() {
        let mut rng = Xoshiro256::new(21);
        let s = sample_without_replacement(&mut rng, 100, 10);
        assert_eq!(s.len(), 10);
    }
    #[test]
    fn sample_without_replacement_unique() {
        let mut rng = Xoshiro256::new(22);
        let mut s = sample_without_replacement(&mut rng, 50, 20);
        s.sort_unstable();
        s.dedup();
        assert_eq!(s.len(), 20);
    }
    #[test]
    fn sample_without_replacement_k_gt_n() {
        let mut rng = Xoshiro256::new(23);
        let s = sample_without_replacement(&mut rng, 5, 10);
        assert_eq!(s.len(), 5);
    }
    #[test]
    fn weighted_choice_selects_valid() {
        let mut rng = Xoshiro256::new(24);
        let choices = vec![("a", 1.0), ("b", 2.0), ("c", 3.0)];
        for _ in 0..100 {
            let v = weighted_choice(&mut rng, &choices).unwrap();
            assert!(matches!(v, "a" | "b" | "c"));
        }
    }
    #[test]
    fn weighted_choice_empty_none() {
        let mut rng = Xoshiro256::new(25);
        let choices: Vec<(i32, f64)> = vec![];
        assert!(weighted_choice(&mut rng, &choices).is_none());
    }
    #[test]
    fn box_muller_produces_two_values() {
        let mut rng = Xoshiro256::new(26);
        let (z0, z1) = box_muller(&mut rng);
        assert!(z0.is_finite());
        assert!(z1.is_finite());
    }
    #[test]
    fn first_n_primes_correct() {
        let p = first_n_primes(5);
        assert_eq!(p, vec![2, 3, 5, 7, 11]);
    }
    #[test]
    fn gamma_sample_positive() {
        let mut rng = Xoshiro256::new(27);
        for _ in 0..100 {
            let v = sample_gamma(&mut rng, 2.0, 1.0);
            assert!(v > 0.0, "gamma sample should be positive");
        }
    }
    #[test]
    fn poisson_sample_n_count() {
        let dist = PoissonDist::new(5.0);
        let mut rng = Xoshiro256::new(28);
        let s = dist.sample_n(&mut rng, 30);
        assert_eq!(s.len(), 30);
    }
    #[test]
    fn weighted_choice_index_boundary() {
        let probs = vec![0.5, 0.3, 0.2];
        assert_eq!(weighted_choice_index(&probs, 0.0), 0);
        assert_eq!(weighted_choice_index(&probs, 0.6), 1);
        assert_eq!(weighted_choice_index(&probs, 0.99), 2);
    }
}
/// Euler–Mascheroni constant γ ≈ 0.5772156649.
pub const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;
/// Sample from a truncated normal distribution N(mu, sigma²) on \[a, b\]
/// using rejection sampling.
pub fn sample_truncated_normal(rng: &mut Xoshiro256, mu: f64, sigma: f64, a: f64, b: f64) -> f64 {
    loop {
        let u1 = rng.next_f64().max(1e-300);
        let u2 = rng.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * PI * u2;
        let x = mu + sigma * r * theta.cos();
        if x >= a && x <= b {
            return x;
        }
    }
}
/// Compute scrambled Halton value at index i with base and permutation.
pub fn scrambled_halton_1d(mut i: usize, base: usize, perm: &[usize]) -> f64 {
    let mut f = 1.0f64;
    let mut r = 0.0f64;
    while i > 0 {
        f /= base as f64;
        let digit = i % base;
        let scrambled_digit = if digit < perm.len() {
            perm[digit]
        } else {
            digit
        };
        r += f * scrambled_digit as f64;
        i /= base;
    }
    r
}
/// Log-gamma function using Lanczos approximation.
///
/// Accurate to ~10 decimal places for Re(z) > 0.
pub fn log_gamma(z: f64) -> f64 {
    if z <= 0.0 {
        return f64::INFINITY;
    }
    if z < 0.5 {
        return PI.ln() - (PI * z).sin().ln() - log_gamma(1.0 - z);
    }
    pub(super) const G: f64 = 7.0;
    pub(super) const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];
    let zz = z - 1.0;
    let mut x = C[0];
    for (i, &c) in C[1..].iter().enumerate() {
        x += c / (zz + (i + 1) as f64);
    }
    let t = zz + G + 0.5;
    0.5 * (2.0 * PI).ln() + (zz + 0.5) * t.ln() - t + x.ln()
}
/// Gamma function Γ(z) = exp(log_gamma(z)).
pub fn gamma_fn(z: f64) -> f64 {
    log_gamma(z).exp()
}
/// Standard normal CDF Φ(x) using rational approximation.
///
/// Accuracy: |error| < 7.5e-8 for all x (Abramowitz & Stegun 26.2.17).
pub fn normal_cdf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t
        * (0.319381530
            + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    let p = 1.0 - ((-0.5 * x * x).exp() / (2.0 * PI).sqrt()) * poly;
    if x >= 0.0 { p } else { 1.0 - p }
}
/// Inverse normal CDF (probit function) via rational approximation.
///
/// Reference: Peter Acklam's approximation (max error ~1.15e-9).
pub fn normal_ppf(p: f64) -> f64 {
    let p = p.clamp(1e-15, 1.0 - 1e-15);
    pub(super) const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383_577_518_672_69e2,
        -3.066479806614716e+01,
        2.506628277459239,
    ];
    pub(super) const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    pub(super) const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    ];
    pub(super) const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996,
        3.754408661907416,
    ];
    pub(super) const P_LOW: f64 = 0.02425;
    pub(super) const P_HIGH: f64 = 1.0 - P_LOW;
    if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (C[5] + q * (C[4] + q * (C[3] + q * (C[2] + q * (C[1] + q * C[0])))))
            / (1.0 + q * (D[3] + q * (D[2] + q * (D[1] + q * D[0]))))
    } else if p <= P_HIGH {
        let q = p - 0.5;
        let r = q * q;
        (q * (A[5] + r * (A[4] + r * (A[3] + r * (A[2] + r * (A[1] + r * A[0]))))))
            / (1.0 + r * (B[4] + r * (B[3] + r * (B[2] + r * (B[1] + r * B[0])))))
    } else {
        let q = (-(2.0 * (1.0 - p).ln())).sqrt();
        -(C[5] + q * (C[4] + q * (C[3] + q * (C[2] + q * (C[1] + q * C[0])))))
            / (1.0 + q * (D[3] + q * (D[2] + q * (D[1] + q * D[0]))))
    }
}
/// Approximate Student-t CDF using numerical integration.
///
/// Uses Simpson's rule for moderate ν; accurate to ~1e-6.
pub fn student_t_cdf(x: f64, nu: f64) -> f64 {
    if nu > 100.0 {
        return normal_cdf(x);
    }
    let t_sq = x * x;
    let z = nu / (nu + t_sq);
    let ib = regularized_incomplete_beta(nu / 2.0, 0.5, z);
    if x >= 0.0 { 1.0 - 0.5 * ib } else { 0.5 * ib }
}
/// Regularized incomplete beta function I_x(a, b).
///
/// Uses a continued fraction expansion (Lentz's method).
pub fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularized_incomplete_beta(b, a, 1.0 - x);
    }
    let ln_beta = log_gamma(a) + log_gamma(b) - log_gamma(a + b);
    let front = (a * x.ln() + b * (1.0 - x).ln() - ln_beta - a.ln()).exp();
    let mut cf = 1.0f64;
    for m in (1..=20).rev() {
        let mf = m as f64;
        let d1 = mf * (b - mf) * x / ((a + 2.0 * mf - 1.0) * (a + 2.0 * mf));
        let d2 = -(a + mf) * (a + b + mf) * x / ((a + 2.0 * mf) * (a + 2.0 * mf + 1.0));
        cf = 1.0 + d1 / (1.0 + d2 / cf);
    }
    front / (a * cf.max(1e-300))
}
#[cfg(test)]
mod extended_tests {

    use crate::bayesian_inference::log_gamma;
    use crate::normal_cdf;
    use crate::random::AntitheticVariates;
    use crate::random::BetaDist;
    use crate::random::BrownianBridge;
    use crate::random::BrownianMotion;
    use crate::random::ChiSquaredDist;
    use crate::random::ClaytonCopula;

    use crate::random::EULER_MASCHERONI;

    use crate::random::FDist;
    use crate::random::FractionalBrownianMotion;
    use crate::random::FrankCopula;
    use crate::random::FrechetDist;
    use crate::random::GEVDist;
    use crate::random::GOEMatrix;
    use crate::random::GUEMatrix;
    use crate::random::GammaDist;
    use crate::random::GaussianCopula;
    use crate::random::GumbelCopula;
    use crate::random::GumbelDist;
    use crate::random::LevyProcess;

    use crate::random::OrnsteinUhlenbeck;

    use crate::random::RejectionSampler;
    use crate::random::ScrambledHalton;
    use crate::random::SobolSequence;

    use crate::random::StratifiedMonteCarlo;
    use crate::random::StudentTDist;
    use crate::random::TCopula;

    use crate::random::WeibullDist;
    use crate::random::WishartMatrix;
    use crate::random::Xoshiro256;
    use crate::random::ZigguratNormal;
    use crate::random::functions::PI;
    use crate::random::normal_ppf;
    use crate::random::sample_truncated_normal;
    #[test]
    fn gamma_dist_positive_samples() {
        let dist = GammaDist::new(2.0, 1.0);
        let mut rng = Xoshiro256::new(100);
        for _ in 0..200 {
            assert!(dist.sample(&mut rng) > 0.0);
        }
    }
    #[test]
    fn gamma_dist_mean_approx() {
        let dist = GammaDist::new(3.0, 2.0);
        let mut rng = Xoshiro256::new(101);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!((mean - 6.0).abs() < 0.3, "mean={mean}");
    }
    #[test]
    fn gamma_dist_variance_approx() {
        let alpha = 4.0;
        let beta = 2.0;
        let dist = GammaDist::new(alpha, beta);
        let mut rng = Xoshiro256::new(102);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
        let expected_var = alpha * beta * beta;
        assert!(
            (var - expected_var).abs() / expected_var < 0.1,
            "var={var}, expected={expected_var}"
        );
    }
    #[test]
    fn beta_dist_in_range() {
        let dist = BetaDist::new(2.0, 3.0);
        let mut rng = Xoshiro256::new(110);
        for _ in 0..200 {
            let v = dist.sample(&mut rng);
            assert!(v > 0.0 && v < 1.0, "v={v}");
        }
    }
    #[test]
    fn beta_dist_mean_approx() {
        let dist = BetaDist::new(2.0, 3.0);
        let mut rng = Xoshiro256::new(111);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let expected = 2.0 / 5.0;
        assert!((mean - expected).abs() < 0.02, "mean={mean}");
    }
    #[test]
    fn beta_log_pdf_finite() {
        let dist = BetaDist::new(2.0, 2.0);
        let lp = dist.log_pdf(0.5);
        assert!(lp.is_finite(), "log_pdf={lp}");
    }
    #[test]
    fn student_t_mean_approx() {
        let dist = StudentTDist::new(10.0);
        let mut rng = Xoshiro256::new(120);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!(mean.abs() < 0.1, "mean={mean}");
    }
    #[test]
    fn student_t_variance_approx() {
        let nu = 10.0;
        let dist = StudentTDist::new(nu);
        let mut rng = Xoshiro256::new(121);
        let samples = dist.sample_n(&mut rng, 20000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
        let expected_var = nu / (nu - 2.0);
        assert!(
            (var - expected_var).abs() / expected_var < 0.15,
            "var={var}"
        );
    }
    #[test]
    fn student_t_log_pdf_finite() {
        let dist = StudentTDist::new(5.0);
        let lp = dist.log_pdf(0.0);
        assert!(lp.is_finite());
    }
    #[test]
    fn f_dist_positive_samples() {
        let dist = FDist::new(5.0, 10.0);
        let mut rng = Xoshiro256::new(130);
        for _ in 0..200 {
            assert!(dist.sample(&mut rng) > 0.0);
        }
    }
    #[test]
    fn f_dist_mean_approx() {
        let dist = FDist::new(2.0, 10.0);
        let mut rng = Xoshiro256::new(131);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let expected = 10.0 / (10.0 - 2.0);
        assert!((mean - expected).abs() / expected < 0.2, "mean={mean}");
    }
    #[test]
    fn chi2_positive_samples() {
        let dist = ChiSquaredDist::new(5.0);
        let mut rng = Xoshiro256::new(140);
        for _ in 0..200 {
            assert!(dist.sample(&mut rng) > 0.0);
        }
    }
    #[test]
    fn chi2_mean_approx() {
        let k = 4.0;
        let dist = ChiSquaredDist::new(k);
        let mut rng = Xoshiro256::new(141);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!((mean - k).abs() / k < 0.1, "mean={mean}, expected={k}");
    }
    #[test]
    fn weibull_positive_samples() {
        let dist = WeibullDist::new(2.0, 1.0);
        let mut rng = Xoshiro256::new(150);
        for _ in 0..200 {
            assert!(dist.sample(&mut rng) > 0.0);
        }
    }
    #[test]
    fn weibull_pdf_nonneg() {
        let dist = WeibullDist::new(2.0, 1.0);
        for i in 1..10 {
            assert!(dist.pdf(i as f64 * 0.1) >= 0.0);
        }
    }
    #[test]
    fn gumbel_samples_finite() {
        let dist = GumbelDist::new(0.0, 1.0);
        let mut rng = Xoshiro256::new(160);
        for _ in 0..200 {
            let v = dist.sample(&mut rng);
            assert!(v.is_finite(), "v={v}");
        }
    }
    #[test]
    fn gumbel_mean_approx() {
        let mu = 2.0;
        let beta = 1.0;
        let dist = GumbelDist::new(mu, beta);
        let mut rng = Xoshiro256::new(161);
        let samples = dist.sample_n(&mut rng, 10000);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let expected = mu + beta * EULER_MASCHERONI;
        assert!(
            (mean - expected).abs() < 0.1,
            "mean={mean}, expected={expected}"
        );
    }
    #[test]
    fn frechet_samples_finite() {
        let dist = FrechetDist::new(2.0, 1.0, 0.0);
        let mut rng = Xoshiro256::new(170);
        for _ in 0..200 {
            let v = dist.sample(&mut rng);
            assert!(v.is_finite() && v > 0.0, "v={v}");
        }
    }
    #[test]
    fn frechet_pdf_nonneg() {
        let dist = FrechetDist::new(2.0, 1.0, 0.0);
        for i in 1..10 {
            assert!(dist.pdf(i as f64 * 0.5) >= 0.0);
        }
    }
    #[test]
    fn gev_gumbel_samples_finite() {
        let dist = GEVDist::new(0.0, 1.0, 0.0);
        let mut rng = Xoshiro256::new(180);
        for _ in 0..200 {
            let v = dist.sample(&mut rng);
            assert!(v.is_finite());
        }
    }
    #[test]
    fn gev_frechet_samples_finite() {
        let dist = GEVDist::new(0.0, 1.0, 0.5);
        let mut rng = Xoshiro256::new(181);
        for _ in 0..200 {
            let v = dist.sample(&mut rng);
            assert!(v.is_finite());
        }
    }
    #[test]
    fn ziggurat_normal_samples_finite() {
        let sampler = ZigguratNormal::new(256);
        let mut rng = Xoshiro256::new(190);
        for _ in 0..200 {
            assert!(sampler.sample(&mut rng).is_finite());
        }
    }
    #[test]
    fn brownian_path_length() {
        let bm = BrownianMotion::new(1.0, 100, 0.01);
        let mut rng = Xoshiro256::new(200);
        let path = bm.generate_path(&mut rng);
        assert_eq!(path.len(), 101);
        assert_eq!(path[0], 0.0);
    }
    #[test]
    fn brownian_path_3d_finite() {
        let bm = BrownianMotion::new(1.0, 50, 0.01);
        let mut rng = Xoshiro256::new(201);
        let path = bm.generate_path_3d(&mut rng);
        for p in &path {
            assert!(p[0].is_finite() && p[1].is_finite() && p[2].is_finite());
        }
    }
    #[test]
    fn brownian_msd_formula() {
        let bm = BrownianMotion::new(1.0, 1, 0.01);
        let msd = bm.theoretical_msd_3d(1);
        assert!((msd - 0.06).abs() < 1e-10, "msd={msd}");
    }
    #[test]
    fn brownian_bridge_endpoints() {
        let bb = BrownianBridge::new(0.0, 1.0, 1.0, 100, 1.0);
        let mut rng = Xoshiro256::new(210);
        let path = bb.generate_path(&mut rng);
        assert_eq!(path.len(), 101);
        assert!((path[0] - 0.0).abs() < 1e-10);
    }
    #[test]
    fn brownian_bridge_variance_midpoint() {
        let t = 0.5;
        let bb = BrownianBridge::new(0.0, 0.0, 1.0, 10, 1.0);
        let var = bb.variance_at(t);
        assert!((var - 0.25).abs() < 1e-10, "var={var}");
    }
    #[test]
    fn ou_path_length() {
        let ou = OrnsteinUhlenbeck::new(1.0, 0.0, 1.0, 100, 0.01);
        let mut rng = Xoshiro256::new(220);
        let path = ou.generate_path(&mut rng, 5.0);
        assert_eq!(path.len(), 101);
    }
    #[test]
    fn ou_stationary_variance() {
        let theta = 2.0;
        let sigma = 1.0;
        let ou = OrnsteinUhlenbeck::new(theta, 0.0, sigma, 5000, 0.005);
        let mut rng = Xoshiro256::new(221);
        let path = ou.generate_path(&mut rng, 0.0);
        let tail = &path[2000..];
        let mean = tail.iter().sum::<f64>() / tail.len() as f64;
        let var = tail.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / tail.len() as f64;
        let expected_var = sigma * sigma / (2.0 * theta);
        assert!(
            (var - expected_var).abs() / expected_var < 0.3,
            "var={var}, expected={expected_var}"
        );
    }
    #[test]
    fn ou_mean_reversion() {
        let ou = OrnsteinUhlenbeck::new(1.0, 2.0, 0.5, 10, 0.1);
        let mean_t = ou.mean_at(10.0, 5.0);
        assert!((mean_t - 2.0).abs() < 1.0);
    }
    #[test]
    fn fbm_path_length() {
        let fbm = FractionalBrownianMotion::new(0.7, 100, 0.01);
        let mut rng = Xoshiro256::new(230);
        let path = fbm.generate_path(&mut rng);
        assert_eq!(path.len(), 101);
        assert_eq!(path[0], 0.0);
    }
    #[test]
    fn fbm_hurst_variance_scaling() {
        let fbm = FractionalBrownianMotion::new(0.8, 1, 1.0);
        let var = fbm.theoretical_variance(1.0);
        assert!((var - 1.0).abs() < 1e-10);
    }
    #[test]
    fn goe_symmetric() {
        let mut rng = Xoshiro256::new(240);
        let m = GOEMatrix::new(5, &mut rng);
        for i in 0..5 {
            for j in 0..5 {
                assert!(
                    (m.get(i, j) - m.get(j, i)).abs() < 1e-15,
                    "GOE not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn goe_frobenius_positive() {
        let mut rng = Xoshiro256::new(241);
        let m = GOEMatrix::new(4, &mut rng);
        assert!(m.frobenius_sq() > 0.0);
    }
    #[test]
    fn gue_hermitian_real_symmetric() {
        let mut rng = Xoshiro256::new(250);
        let m = GUEMatrix::new(4, &mut rng);
        for i in 0..4 {
            for j in 0..4 {
                assert!((m.real_part[i * 4 + j] - m.real_part[j * 4 + i]).abs() < 1e-15);
                assert!((m.imag_part[i * 4 + j] + m.imag_part[j * 4 + i]).abs() < 1e-15);
            }
        }
    }
    #[test]
    fn wishart_trace_approx_n_times_p() {
        let n = 4;
        let p = 100;
        let mut rng = Xoshiro256::new(260);
        let w = WishartMatrix::new(n, p, &mut rng);
        let tr = w.trace();
        assert!((tr / (n * p) as f64 - 1.0).abs() < 0.2, "trace={tr}");
    }
    #[test]
    fn wishart_symmetric() {
        let mut rng = Xoshiro256::new(261);
        let w = WishartMatrix::new(3, 5, &mut rng);
        for i in 0..3 {
            for j in 0..3 {
                assert!((w.get(i, j) - w.get(j, i)).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn sobol_sequence_range() {
        let mut sobol = SobolSequence::new(3);
        for _ in 0..100 {
            let pt = sobol.next_point();
            for &v in &pt {
                assert!((0.0..=1.0).contains(&v), "v={v}");
            }
        }
    }
    #[test]
    fn sobol_sequence_generate_count() {
        let mut sobol = SobolSequence::new(2);
        let pts = sobol.generate(50);
        assert_eq!(pts.len(), 50);
    }
    #[test]
    fn scrambled_halton_range() {
        let mut rng = Xoshiro256::new(270);
        let mut sh = ScrambledHalton::new(3, &mut rng);
        for _ in 0..100 {
            let pt = sh.next_point();
            for &v in &pt {
                assert!((0.0..=1.0).contains(&v), "v={v}");
            }
        }
    }
    #[test]
    fn gaussian_copula_in_unit_square() {
        let cop = GaussianCopula::new(0.7);
        let mut rng = Xoshiro256::new(280);
        for _ in 0..200 {
            let (u, v) = cop.sample(&mut rng);
            assert!(u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0, "u={u}, v={v}");
        }
    }
    #[test]
    fn gaussian_copula_kendall_tau() {
        let rho = 0.6;
        let cop = GaussianCopula::new(rho);
        let tau = cop.kendall_tau();
        let expected = (2.0 / PI) * rho.asin();
        assert!((tau - expected).abs() < 1e-10);
    }
    #[test]
    fn clayton_copula_in_unit_square() {
        let cop = ClaytonCopula::new(2.0);
        let mut rng = Xoshiro256::new(290);
        for _ in 0..200 {
            let (u, v) = cop.sample(&mut rng);
            assert!(u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0, "u={u}, v={v}");
        }
    }
    #[test]
    fn clayton_copula_kendall_tau() {
        let theta = 2.0;
        let cop = ClaytonCopula::new(theta);
        let tau = cop.kendall_tau();
        assert!((tau - theta / (theta + 2.0)).abs() < 1e-10);
    }
    #[test]
    fn gumbel_copula_in_unit_square() {
        let cop = GumbelCopula::new(2.0);
        let mut rng = Xoshiro256::new(300);
        for _ in 0..200 {
            let (u, v) = cop.sample(&mut rng);
            assert!(u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0, "u={u}, v={v}");
        }
    }
    #[test]
    fn gumbel_copula_kendall_tau() {
        let theta = 3.0;
        let cop = GumbelCopula::new(theta);
        assert!((cop.kendall_tau() - (1.0 - 1.0 / theta)).abs() < 1e-10);
    }
    #[test]
    fn frank_copula_in_unit_square() {
        let cop = FrankCopula::new(4.0);
        let mut rng = Xoshiro256::new(310);
        for _ in 0..200 {
            let (u, v) = cop.sample(&mut rng);
            assert!(u > 0.0 && u < 1.0 && v > 0.0 && v < 1.0, "u={u}, v={v}");
        }
    }
    #[test]
    fn antithetic_estimates_integral_linear() {
        let av = AntitheticVariates::new(5000);
        let mut rng = Xoshiro256::new(320);
        let est = av.estimate(&mut rng, &|x| x);
        assert!((est - 0.5).abs() < 0.01, "est={est}");
    }
    #[test]
    fn antithetic_estimates_integral_quadratic() {
        let av = AntitheticVariates::new(5000);
        let mut rng = Xoshiro256::new(321);
        let est = av.estimate(&mut rng, &|x| x * x);
        assert!((est - 1.0 / 3.0).abs() < 0.02, "est={est}");
    }
    #[test]
    fn stratified_mc_1d_integral() {
        let smc = StratifiedMonteCarlo::new(100, 1);
        let mut rng = Xoshiro256::new(330);
        let est = smc.estimate(&mut rng, &|x| x[0]);
        assert!((est - 0.5).abs() < 0.01, "est={est}");
    }
    #[test]
    fn log_gamma_known_values() {
        assert!(
            (log_gamma(1.0) - 0.0).abs() < 1e-6,
            "Γ(1)={}",
            log_gamma(1.0).exp()
        );
        assert!((log_gamma(2.0) - 0.0).abs() < 1e-6);
        assert!((log_gamma(3.0) - 2.0_f64.ln()).abs() < 1e-5);
        assert!((log_gamma(4.0) - 6.0_f64.ln()).abs() < 1e-4);
    }
    #[test]
    fn normal_cdf_known_values() {
        assert!(
            (normal_cdf(0.0) - 0.5).abs() < 1e-5,
            "Φ(0)={}",
            normal_cdf(0.0)
        );
        assert!((normal_cdf(1.645) - 0.95).abs() < 5e-3);
        assert!((normal_cdf(-1.645) - 0.05).abs() < 5e-3);
    }
    #[test]
    fn normal_ppf_inverse_of_cdf() {
        for p in [0.1, 0.25, 0.5, 0.75, 0.9] {
            let x = normal_ppf(p);
            let p2 = normal_cdf(x);
            assert!((p - p2).abs() < 1e-4, "p={p}, p2={p2}");
        }
    }
    #[test]
    fn truncated_normal_in_range() {
        let mut rng = Xoshiro256::new(340);
        for _ in 0..100 {
            let v = sample_truncated_normal(&mut rng, 0.0, 1.0, -2.0, 2.0);
            assert!((-2.0..=2.0).contains(&v), "v={v}");
        }
    }
    #[test]
    fn levy_process_cauchy_samples_mostly_finite() {
        let lp = LevyProcess::new(1.0, 0.0, 1.0, 0.0);
        let mut rng = Xoshiro256::new(350);
        let samples: Vec<f64> = (0..1000).map(|_| lp.sample(&mut rng)).collect();
        let count_finite = samples.iter().filter(|v| v.is_finite()).count();
        assert!(
            count_finite > 900,
            "too many non-finite samples: {count_finite}/1000"
        );
    }
    #[test]
    fn t_copula_in_unit_square() {
        let cop = TCopula::new(0.5, 5.0);
        let mut rng = Xoshiro256::new(360);
        for _ in 0..200 {
            let (u, v) = cop.sample(&mut rng);
            assert!(
                (0.0..=1.0).contains(&u) && (0.0..=1.0).contains(&v),
                "u={u}, v={v}"
            );
        }
    }
    #[test]
    fn rejection_sampler_produces_values_in_range() {
        let rs = RejectionSampler::new(2.0);
        let mut rng = Xoshiro256::new(370);
        let target = |x: f64| -> f64 {
            if (0.0..=1.0).contains(&x) {
                2.0 * (0.5 - (x - 0.5).abs())
            } else {
                0.0
            }
        };
        let mut propose = |r: &mut Xoshiro256| -> (f64, f64) {
            let x = r.next_f64();
            (x, 1.0)
        };
        for _ in 0..50 {
            let v = rs.sample(&mut rng, &mut propose, &target);
            assert!((0.0..=1.0).contains(&v));
        }
    }
}
