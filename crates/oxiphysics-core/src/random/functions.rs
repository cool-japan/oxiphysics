//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::Xoshiro256;
use crate::numerics::{inc_beta, inc_beta_xy};

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
        let samples = RandomSampler::sobol(50, 4).unwrap();
        assert_eq!(samples.len(), 50);
        assert_eq!(samples[0].len(), 4);
    }
    #[test]
    fn sobol_range() {
        let samples = RandomSampler::sobol(100, 2).unwrap();
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
/// Degrees of freedom from which [`student_t_cdf`] prefers the large-`nu` expansion.
const STUDENT_T_EXPANSION_MIN_NU: f64 = 1.0e4;
/// From this many degrees of freedom on, every `x` outside the expansion region
/// `|x| <= nu^(1/4) / 2` has a tail probability far below `1e-300`.
const STUDENT_T_EXPANSION_ONLY_NU: f64 = 1.0e8;

/// Upper tail `P(Z > a)` of the standard normal distribution.
///
/// With `z = a / sqrt(2)` the tail is `erfc(z) / 2`.  For `z < 1` it is computed
/// as `1 - erf(z)` from the Maclaurin series
/// `erf(z) = 2/sqrt(pi) exp(-z²) Σ 2^n z^(2n+1) / (2n+1)!!` (all terms positive, so
/// there is no cancellation); for `z >= 1` Lentz's continued fraction
/// `erfc(z) = exp(-z²)/sqrt(pi) / (z + (1/2)/(z + 1/(z + (3/2)/(z + ...))))`
/// keeps a relative accuracy of about `1e-13` far into the tail.  The result is
/// accurate to a few ulp for moderate `a`, unlike the `7.5e-8` approximation of
/// [`normal_cdf`].
fn standard_normal_upper_tail(a: f64) -> f64 {
    const MAX_TERMS: usize = 2000;
    const TINY: f64 = 1.0e-300;
    if a.is_nan() {
        return f64::NAN;
    }
    if a < 0.0 {
        return 1.0 - standard_normal_upper_tail(-a);
    }
    let z = a / std::f64::consts::SQRT_2;
    if z == 0.0 {
        return 0.5;
    }
    if z < 1.0 {
        let z_sq = z * z;
        let mut term = z;
        let mut sum = z;
        for n in 1..=200_u32 {
            term *= 2.0 * z_sq / f64::from(2 * n + 1);
            sum += term;
            if term < 1.0e-17 * sum {
                break;
            }
        }
        return 0.5 * (1.0 - 2.0 / PI.sqrt() * (-z_sq).exp() * sum);
    }
    let mut fraction = z;
    let mut c = z;
    let mut d = 0.0_f64;
    for k in 1..=MAX_TERMS {
        let a_k = k as f64 / 2.0;
        d = z + a_k * d;
        if d == 0.0 {
            d = TINY;
        }
        c = z + a_k / c;
        if c == 0.0 {
            c = TINY;
        }
        d = d.recip();
        let delta = c * d;
        fraction *= delta;
        if (delta - 1.0).abs() < 1.0e-16 {
            break;
        }
    }
    0.5 * (-z * z).exp() / (PI.sqrt() * fraction)
}
/// The correction `Σ_{k=1..4} g_k(x) / nu^k` in the large-`nu` expansion
/// `P(T <= x) = Φ(x) - φ(x) Σ_k g_k(x) / nu^k` of the Student-t CDF about the
/// normal law (Abramowitz & Stegun 26.7.8, carried to the fourth order), with
/// `g_1 = x (x² + 1) / 4`, `g_2 = x (x² - 3) (3 x⁴ + 2 x² + 1) / 96`, ...
fn student_t_expansion_correction(x: f64, nu: f64) -> f64 {
    let x2 = x * x;
    let inv_nu = nu.recip();
    let g1 = x * (x2 + 1.0) / 4.0;
    let g2 = x * (x2 - 3.0) * (3.0 * x2 * x2 + 2.0 * x2 + 1.0) / 96.0;
    let g3 = x * (((((x2 - 11.0) * x2 + 14.0) * x2 + 6.0) * x2 - 3.0) * x2 - 15.0) / 384.0;
    let g4 = x
        * (((((((15.0 * x2 - 375.0) * x2 + 2225.0) * x2 - 2141.0) * x2 - 939.0) * x2 - 213.0)
            * x2
            + 915.0)
            * x2
            + 945.0)
        / 92160.0;
    (((g4 * inv_nu + g3) * inv_nu + g2) * inv_nu + g1) * inv_nu
}
/// Upper tail `P(T > ax)` for `ax >= 0` from the large-`nu` expansion: the exactly
/// evaluated normal tail plus `φ(ax)` times the expansion correction.
fn student_t_upper_tail_expansion(ax: f64, nu: f64) -> f64 {
    let pdf = (-0.5 * ax * ax).exp() / (2.0 * PI).sqrt();
    (standard_normal_upper_tail(ax) + pdf * student_t_expansion_correction(ax, nu)).clamp(0.0, 1.0)
}
/// Student-t cumulative distribution function `P(T <= x)` with `nu` degrees of freedom.
///
/// * In general the CDF is evaluated exactly through the regularized incomplete
///   beta function: for `x >= 0`, `P(T <= x) = 1 - I_z(nu/2, 1/2) / 2` with
///   `z = nu / (nu + x²)`; the lower tail follows from the symmetry of the
///   density.  The complement `1 - z = x² / (nu + x²)` is formed directly, so the
///   tails keep their full relative accuracy.  The absolute error is about
///   `1e-15 * nu` (the cancellation in `ln B(nu/2, 1/2)`), below `1e-11` for
///   `nu < 1e4`.
/// * For very many degrees of freedom the incomplete beta function loses digits
///   (its argument `z` rounds to 1).  From `nu = 1e4` on and for
///   `|x| <= nu^(1/4) / 2` the CDF is therefore computed from the expansion
///   about the normal law, `Φ(x) - φ(x) Σ g_k(x) / nu^k` through the fourth
///   order (Abramowitz & Stegun 26.7.8), applied to the exactly evaluated normal
///   tail.  It is accurate to about `1e-16` absolute and `1e-12` relative in the
///   tails for every larger `nu`.
///
/// `nu = +inf` returns the normal CDF; `NaN` input or `nu <= 0` returns `NaN`.
pub fn student_t_cdf(x: f64, nu: f64) -> f64 {
    if x.is_nan() || nu.is_nan() || nu <= 0.0 {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x > 0.0 { 1.0 } else { 0.0 };
    }
    let ax = x.abs();
    if nu == f64::INFINITY {
        let tail = standard_normal_upper_tail(ax);
        return if x >= 0.0 { 1.0 - tail } else { tail };
    }
    if nu >= STUDENT_T_EXPANSION_MIN_NU {
        if ax <= 0.5 * nu.powf(0.25) {
            let tail = student_t_upper_tail_expansion(ax, nu);
            return if x >= 0.0 { 1.0 - tail } else { tail };
        }
        if nu >= STUDENT_T_EXPANSION_ONLY_NU {
            return if x >= 0.0 { 1.0 } else { 0.0 };
        }
    }
    let t_sq = ax * ax;
    let (z, one_minus_z) = if t_sq.is_finite() {
        let denom = nu + t_sq;
        (nu / denom, t_sq / denom)
    } else {
        // x² overflows: use nu / x² instead, which is finite and tiny.
        let ratio = nu / ax / ax;
        (ratio / (1.0 + ratio), 1.0 / (1.0 + ratio))
    };
    let ib = inc_beta_xy(z, one_minus_z, nu / 2.0, 0.5);
    if x >= 0.0 { 1.0 - 0.5 * ib } else { 0.5 * ib }
}
/// Regularized incomplete beta function `I_x(a, b)`.
///
/// Evaluates the Beta(`a`, `b`) CDF with the continued fraction of
/// [`crate::numerics::inc_beta`], whose argument order is `(x, a, b)`; note that
/// this function takes the shape parameters first: `(a, b, x)`.
pub fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    inc_beta(x, a, b)
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

    use crate::random::RandomSampler;
    use crate::random::RejectionSampler;
    use crate::random::ScrambledHalton;
    use crate::random::SobolSequence;

    use crate::random::StratifiedMonteCarlo;
    use crate::random::StudentTDist;
    use crate::random::TCopula;

    use crate::error::PhysicsError;
    use crate::optimization::cholesky_lower;
    use crate::random::SOBOL_MAX_DIM;
    use crate::random::WeibullDist;
    use crate::random::WishartMatrix;
    use crate::random::Xoshiro256;
    use crate::random::ZigguratNormal;
    use crate::random::box_muller;
    use crate::random::functions::PI;
    use crate::random::normal_ppf;
    use crate::random::regularized_incomplete_beta;
    use crate::random::sample_truncated_normal;
    use crate::random::student_t_cdf;
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
        let mut sobol = SobolSequence::new(3).unwrap();
        for _ in 0..100 {
            let pt = sobol.next_point();
            for &v in &pt {
                assert!((0.0..=1.0).contains(&v), "v={v}");
            }
        }
    }
    #[test]
    fn sobol_sequence_generate_count() {
        let mut sobol = SobolSequence::new(2).unwrap();
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
    /// `(x, nu, P(T <= x))` evaluated with 50-digit arithmetic (mpmath, cross-checked with SciPy).
    const STUDENT_T_REFERENCE: [(f64, f64, f64); 14] = [
        (1.0, 1.0, 0.75),
        (2.0, 10.0, 0.9633059826146299),
        (2.0, 5.0, 0.9490302605850708),
        (0.0, 3.0, 0.5),
        (-1.5, 4.0, 0.104),
        (3.0, 2.0, 0.9522670168666454),
        (-0.5, 30.0, 0.31036150244256366),
        (2.0, 100.0, 0.9758939106344332),
        (2.0, 101.0, 0.9759073943963047),
        (1.96, 1000.0, 0.9748634075221256),
        (-2.5, 7.5, 0.019410129136812778),
        (10.0, 3.0, 0.9989358004707929),
        (0.1, 1.0, 0.5317255174305535),
        (1.5, 0.5, 0.7460491046759458),
    ];
    #[test]
    fn regularized_incomplete_beta_known_values() {
        assert!((regularized_incomplete_beta(1.0, 1.0, 0.3) - 0.3).abs() < 1e-12);
        assert!((regularized_incomplete_beta(2.0, 3.0, 0.2) - 0.1808).abs() < 1e-12);
        assert!((regularized_incomplete_beta(2.0, 3.0, 0.6) - 0.8208).abs() < 1e-12);
        assert!((regularized_incomplete_beta(0.5, 0.5, 0.5) - 0.5).abs() < 1e-12);
        assert_eq!(regularized_incomplete_beta(2.0, 3.0, 0.0), 0.0);
        assert_eq!(regularized_incomplete_beta(2.0, 3.0, 1.0), 1.0);
    }
    #[test]
    fn regularized_incomplete_beta_argument_order_is_shape_first() {
        // (a, b, x) here, (x, a, b) in `numerics::inc_beta`.
        for &(a, b, x) in &[
            (2.0_f64, 5.0_f64, 0.3_f64),
            (0.5, 3.0, 0.8),
            (7.0, 0.5, 0.95),
        ] {
            let direct = regularized_incomplete_beta(a, b, x);
            let numerics = crate::numerics::inc_beta(x, a, b);
            assert!((direct - numerics).abs() < 1e-15);
        }
    }
    #[test]
    fn student_t_cdf_matches_high_precision_reference() {
        for (x, nu, expected) in STUDENT_T_REFERENCE.iter() {
            let got = student_t_cdf(*x, *nu);
            assert!(
                (got - expected).abs() < 1e-10,
                "t-cdf({x}, nu={nu}): got {got}, expected {expected}"
            );
        }
    }
    #[test]
    fn student_t_cdf_regression_values_from_audit() {
        // The old implementation returned 0.653, 0.9975 and a wrong value here.
        assert!((student_t_cdf(1.0, 1.0) - 0.75).abs() < 1e-6);
        assert!((student_t_cdf(2.0, 10.0) - 0.963306).abs() < 1e-6);
        assert!((student_t_cdf(2.0, 5.0) - 0.949030).abs() < 1e-6);
    }
    #[test]
    fn student_t_cdf_is_symmetric_and_monotone() {
        for &nu in &[0.7_f64, 1.0, 3.0, 12.0, 80.0, 400.0] {
            assert!((student_t_cdf(0.0, nu) - 0.5).abs() < 1e-14);
            let mut previous = 0.0_f64;
            for i in -60..=60 {
                let x = i as f64 / 10.0;
                let upper = student_t_cdf(x, nu);
                let lower = student_t_cdf(-x, nu);
                assert!(
                    (upper + lower - 1.0).abs() < 1e-13,
                    "symmetry at x={x}, nu={nu}"
                );
                assert!(upper >= previous, "monotone at x={x}, nu={nu}");
                previous = upper;
            }
        }
    }
    #[test]
    fn student_t_cdf_converges_to_normal_for_large_nu() {
        // Phi(1.96) = 0.9750021048517795
        let huge_nu = student_t_cdf(1.96, 1.0e8);
        assert!(
            (huge_nu - 0.9750021048517795).abs() < 1e-7,
            "nu=1e8: {huge_nu}"
        );
        let infinite = student_t_cdf(1.96, f64::INFINITY);
        assert!(
            (infinite - 0.9750021048517795).abs() < 1e-6,
            "nu=inf: {infinite}"
        );
    }
    #[test]
    fn student_t_cdf_special_inputs() {
        assert!(student_t_cdf(f64::NAN, 5.0).is_nan());
        assert!(student_t_cdf(1.0, f64::NAN).is_nan());
        assert!(student_t_cdf(1.0, 0.0).is_nan());
        assert!(student_t_cdf(1.0, -2.0).is_nan());
        assert_eq!(student_t_cdf(f64::INFINITY, 4.0), 1.0);
        assert_eq!(student_t_cdf(f64::NEG_INFINITY, 4.0), 0.0);
        let deep_tail = student_t_cdf(-50.0, 3.0);
        assert!(
            deep_tail > 0.0 && deep_tail < 1e-4,
            "deep tail = {deep_tail}"
        );
    }
    /// `(x, nu, P(T <= x))` for very many degrees of freedom, evaluated with 60-digit
    /// arithmetic (mpmath incomplete beta function).  The old implementation and a
    /// plain incomplete-beta evaluation lose 4 digits at nu = 1e8, 10 at nu = 1e12
    /// and all of them at nu = 1e15.
    const STUDENT_T_LARGE_NU_REFERENCE: [(f64, f64, f64); 45] = [
        (-3.0, 10000.0, 0.0013532240949988331),
        (-1.5, 10000.0, 0.06682298591180981),
        (0.5, 10000.0, 0.6914569603383833),
        (1.96, 10000.0, 0.9749882398840835),
        (4.0, 10000.0, 0.9999681006655884),
        (-3.0, 30000.0, 0.0013510062356349663),
        (-1.5, 30000.0, 0.06681246288629299),
        (0.5, 30000.0, 0.6914606276098914),
        (1.96, 30000.0, 0.9749974833139915),
        (4.0, 30000.0, 0.9999682528565017),
        (-3.0, 100000.0, 0.0013502304420323595),
        (-1.5, 100000.0, 0.0668087797614153),
        (0.5, 100000.0, 0.691461911172791),
        (1.96, 100000.0, 0.9750007184028492),
        (4.0, 100000.0, 0.9999683060012214),
        (-3.0, 1000000.0, 0.0013499312707108985),
        (-1.5, 1000000.0, 0.06680735911839639),
        (0.5, 1000000.0, 0.6914624062638143),
        (1.96, 1000000.0, 0.9750019662073651),
        (4.0, 1000000.0, 0.9999683264829949),
        (-3.0, 100000000.0, 0.0013498983640187472),
        (-1.5, 100000000.0, 0.06680720284735377),
        (0.5, 100000000.0, 0.691462460723911),
        (1.96, 100000000.0, 0.975002103465336),
        (4.0, 100000000.0, 0.9999683287354157),
        (-3.0, 10000000000.0, 0.0013498980349539809),
        (-1.5, 10000000000.0, 0.06680720128464303),
        (0.5, 10000000000.0, 0.6914624612685121),
        (1.96, 10000000000.0, 0.9750021048379152),
        (4.0, 10000000000.0, 0.9999683287579394),
        (-3.0, 1000000000000.0, 0.0013498980316633334),
        (-1.5, 1000000000000.0, 0.06680720126901592),
        (0.5, 1000000000000.0, 0.691462461273958),
        (1.96, 1000000000000.0, 0.975002104851641),
        (4.0, 1000000000000.0, 0.9999683287581647),
        (-3.0, 1000000000000000.0, 0.0013498980316301278),
        (-1.5, 1000000000000000.0, 0.06680720126885822),
        (0.5, 1000000000000000.0, 0.691462461274013),
        (1.96, 1000000000000000.0, 0.9750021048517794),
        (4.0, 1000000000000000.0, 0.9999683287581669),
        (-3.0, 1e+30, 0.0013498980316300946),
        (-1.5, 1e+30, 0.06680720126885807),
        (0.5, 1e+30, 0.6914624612740131),
        (1.96, 1e+30, 0.9750021048517795),
        (4.0, 1e+30, 0.9999683287581669),
    ];
    /// `(x, nu, P(T <= x))` in the far lower tail, where only a relative error test is meaningful.
    const STUDENT_T_DEEP_TAIL_REFERENCE: [(f64, f64, f64); 6] = [
        (-5.0, 10000.0, 2.9151633068313237e-07),
        (-5.0, 1000000.0, 2.866998935445371e-07),
        (-8.0, 1000000.0, 6.22753171660126e-16),
        (-12.0, 100000000.0, 1.7765754821344812e-33),
        (-20.0, 10000000000.0, 2.7536351881296547e-89),
        (-30.0, 1000000000000000.0, 4.9067139281440034e-198),
    ];
    /// `(a, P(Z > a))` for the standard normal, 40-digit reference values.
    const NORMAL_UPPER_TAIL_REFERENCE: [(f64, f64); 13] = [
        (0.0, 0.5),
        (0.1, 0.460172162722971),
        (0.5, 0.3085375387259869),
        (1.0, 0.15865525393145705),
        (1.4, 0.08075665923377107),
        (1.5, 0.06680720126885807),
        (2.0, 0.02275013194817921),
        (3.0, 0.0013498980316300946),
        (5.0, 2.866515718791939e-07),
        (8.0, 6.220960574271784e-16),
        (12.0, 1.776482112077679e-33),
        (20.0, 2.7536241186062337e-89),
        (37.0, 5.725571222524577e-300),
    ];
    #[test]
    fn standard_normal_upper_tail_matches_high_precision_reference() {
        for (a, expected) in NORMAL_UPPER_TAIL_REFERENCE.iter() {
            let got = super::standard_normal_upper_tail(*a);
            assert!(
                (got - expected).abs() <= 5e-13 * expected,
                "Q({a}): got {got:e}, expected {expected:e}"
            );
        }
        assert!(
            (super::standard_normal_upper_tail(-1.0) - (1.0 - 0.15865525393145705)).abs() < 1e-15
        );
        assert!(super::standard_normal_upper_tail(f64::NAN).is_nan());
    }
    #[test]
    fn student_t_cdf_stays_accurate_for_very_many_degrees_of_freedom() {
        for (x, nu, expected) in STUDENT_T_LARGE_NU_REFERENCE.iter() {
            let got = student_t_cdf(*x, *nu);
            assert!(
                (got - expected).abs() < 1e-14,
                "t-cdf({x}, nu={nu:e}): got {got}, expected {expected}"
            );
        }
    }
    #[test]
    fn student_t_cdf_keeps_relative_accuracy_in_the_far_tails() {
        for (x, nu, expected) in STUDENT_T_DEEP_TAIL_REFERENCE.iter() {
            let got = student_t_cdf(*x, *nu);
            assert!(
                (got - expected).abs() <= 1e-9 * expected,
                "t-cdf({x}, nu={nu:e}): got {got:e}, expected {expected:e}"
            );
        }
    }
    #[test]
    fn student_t_cdf_is_continuous_across_the_expansion_threshold() {
        // nu = 1e4 is where the evaluation switches from the incomplete beta function to
        // the large-nu expansion; both must agree to ~1e-11 (the incomplete beta route
        // loses about nu * 1e-15 to the cancellation in ln B) on either side.
        for &x in &[-4.0_f64, -1.0, 0.3, 2.5, 4.5] {
            let below = student_t_cdf(x, 9999.999999);
            let above = student_t_cdf(x, 10000.0);
            assert!((below - above).abs() < 5e-11, "x={x}: {below} vs {above}");
        }
        // Across the |x| = nu^(1/4)/2 border of the expansion region at nu = 1e4 (x = 5).
        for &nu in &[1.0e4_f64, 1.0e5] {
            let border = 0.5 * nu.powf(0.25);
            let inside = student_t_cdf(-border, nu);
            let outside = student_t_cdf(-border * (1.0 + 1e-12), nu);
            assert!(
                (inside - outside).abs() <= 1e-9 * inside,
                "nu={nu:e}: {inside:e} vs {outside:e}"
            );
        }
    }
    #[test]
    fn student_t_cdf_infinite_nu_is_the_exact_normal_cdf() {
        for (a, tail) in NORMAL_UPPER_TAIL_REFERENCE.iter().take(12) {
            let lower = student_t_cdf(-*a, f64::INFINITY);
            assert!(
                (lower - tail).abs() <= 5e-13 * tail,
                "a={a}: {lower:e} vs {tail:e}"
            );
            let upper = student_t_cdf(*a, f64::INFINITY);
            assert!((upper - (1.0 - tail)).abs() < 1e-15, "a={a}: {upper}");
        }
    }
    #[test]
    fn student_t_cdf_never_returns_nan_for_huge_arguments() {
        for &nu in &[0.5_f64, 3.0, 50.0, 2.0e4, 1.0e9, 1.0e300] {
            for &x in &[1.0e100_f64, 1.0e200, 1.0e300, f64::MAX] {
                let upper = student_t_cdf(x, nu);
                let lower = student_t_cdf(-x, nu);
                assert!(
                    (0.0..=1.0).contains(&upper) && (0.0..=1.0).contains(&lower),
                    "nu={nu:e}, x={x:e}: upper={upper}, lower={lower}"
                );
                assert!(upper > 0.999 && lower < 1e-3, "nu={nu:e}, x={x:e}");
            }
        }
    }
    #[test]
    fn t_copula_margins_are_uniform() {
        // U = F_nu(T) must be Uniform(0, 1): mean 1/2, variance 1/12.
        let cop = TCopula::new(0.6, 5.0);
        let mut rng = Xoshiro256::new(361);
        let samples = cop.sample_n(&mut rng, 6000);
        let n = samples.len() as f64;
        let mean_u = samples.iter().map(|s| s.0).sum::<f64>() / n;
        let var_u = samples.iter().map(|s| (s.0 - mean_u).powi(2)).sum::<f64>() / n;
        let low_fraction = samples.iter().filter(|s| s.1 < 0.1).count() as f64 / n;
        assert!((mean_u - 0.5).abs() < 0.02, "mean of U = {mean_u}");
        assert!(
            (var_u - 1.0 / 12.0).abs() < 0.008,
            "variance of U = {var_u}"
        );
        assert!(
            (low_fraction - 0.1).abs() < 0.02,
            "P(V < 0.1) = {low_fraction}"
        );
    }
    /// First 16 points of the 6-dimensional Sobol sequence (Joe-Kuo direction numbers),
    /// generated with SciPy's `qmc.Sobol(d=6, scramble=False, bits=32)`.
    const SOBOL_REFERENCE_6D: [[f64; 6]; 16] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
        [0.75, 0.25, 0.25, 0.25, 0.75, 0.75],
        [0.25, 0.75, 0.75, 0.75, 0.25, 0.25],
        [0.375, 0.375, 0.625, 0.875, 0.375, 0.125],
        [0.875, 0.875, 0.125, 0.375, 0.875, 0.625],
        [0.625, 0.125, 0.875, 0.625, 0.625, 0.875],
        [0.125, 0.625, 0.375, 0.125, 0.125, 0.375],
        [0.1875, 0.3125, 0.9375, 0.4375, 0.5625, 0.3125],
        [0.6875, 0.8125, 0.4375, 0.9375, 0.0625, 0.8125],
        [0.9375, 0.0625, 0.6875, 0.1875, 0.3125, 0.5625],
        [0.4375, 0.5625, 0.1875, 0.6875, 0.8125, 0.0625],
        [0.3125, 0.1875, 0.3125, 0.5625, 0.9375, 0.4375],
        [0.8125, 0.6875, 0.8125, 0.0625, 0.4375, 0.9375],
        [0.5625, 0.4375, 0.0625, 0.8125, 0.1875, 0.6875],
        [0.0625, 0.9375, 0.5625, 0.3125, 0.6875, 0.1875],
    ];
    #[test]
    fn sobol_first_two_dimensional_points_are_the_published_ones() {
        let pts = RandomSampler::sobol(4, 2).unwrap();
        let expected = [[0.0, 0.0], [0.5, 0.5], [0.75, 0.25], [0.25, 0.75]];
        for (point, want) in pts.iter().zip(expected.iter()) {
            assert_eq!(point.as_slice(), want.as_slice());
        }
    }
    #[test]
    fn sobol_six_dimensional_prefix_matches_reference_implementation() {
        let pts = RandomSampler::sobol(16, 6).unwrap();
        for (i, (point, want)) in pts.iter().zip(SOBOL_REFERENCE_6D.iter()).enumerate() {
            assert_eq!(point.as_slice(), want.as_slice(), "point {i}");
        }
    }
    #[test]
    fn sobol_columns_are_pairwise_distinct() {
        // The old tables made dimension 1 identical to dimension 0 (all points on x = y).
        let dims = 21;
        let pts = RandomSampler::sobol(64, dims).unwrap();
        for a in 0..dims {
            for b in (a + 1)..dims {
                let differing = pts.iter().filter(|p| p[a] != p[b]).count();
                assert!(differing > 32, "columns {a} and {b} nearly coincide");
            }
        }
    }
    #[test]
    fn sobol_first_two_columns_do_not_lie_on_the_diagonal() {
        let pts = RandomSampler::sobol(256, 2).unwrap();
        let off_diagonal = pts.iter().filter(|p| (p[0] - p[1]).abs() > 1e-12).count();
        assert!(
            off_diagonal > 128,
            "only {off_diagonal} of 256 points off the diagonal"
        );
    }
    #[test]
    fn sobol_one_dimensional_projections_are_permutations_of_the_dyadic_grid() {
        let n = 256;
        let pts = RandomSampler::sobol(n, SOBOL_MAX_DIM).unwrap();
        for dim in 0..SOBOL_MAX_DIM {
            let mut cells: Vec<usize> = pts.iter().map(|p| (p[dim] * n as f64) as usize).collect();
            cells.sort_unstable();
            let expected: Vec<usize> = (0..n).collect();
            assert_eq!(
                cells, expected,
                "dimension {dim} is not a (0,1)-sequence prefix"
            );
        }
    }
    #[test]
    fn sobol_first_two_dimensions_form_a_zero_net() {
        // The first 2^m points of the first two Sobol coordinates put exactly one
        // point into every dyadic rectangle of area 2^-m.
        let m = 8usize;
        let pts = RandomSampler::sobol(1 << m, 2).unwrap();
        for a in 0..=m {
            let b = m - a;
            let mut counts = vec![0usize; 1 << m];
            for p in &pts {
                let cx = (p[0] * (1u64 << a) as f64) as usize;
                let cy = (p[1] * (1u64 << b) as f64) as usize;
                counts[cx * (1 << b) + cy] += 1;
            }
            assert!(counts.iter().all(|&c| c == 1), "rectangles 2^-{a} x 2^-{b}");
        }
    }
    #[test]
    fn sobol_dimension_limit_is_an_error_not_a_truncation() {
        assert!(matches!(
            SobolSequence::new(SOBOL_MAX_DIM + 1),
            Err(PhysicsError::InvalidInput { .. })
        ));
        assert!(matches!(
            RandomSampler::sobol(10, 200),
            Err(PhysicsError::InvalidInput { .. })
        ));
        let mut full = SobolSequence::new(SOBOL_MAX_DIM).unwrap();
        assert_eq!(full.next_point().len(), SOBOL_MAX_DIM);
        let five = RandomSampler::sobol(3, 5).unwrap();
        assert!(five.iter().all(|p| p.len() == 5));
    }
    #[test]
    fn sobol_skip_to_matches_sequential_generation() {
        let mut sequential = SobolSequence::new(9).unwrap();
        let points = sequential.generate(300);
        for &start in &[0usize, 1, 2, 7, 64, 255, 299] {
            let mut jumped = SobolSequence::new(9).unwrap();
            jumped.skip_to(start);
            for (offset, want) in points[start..].iter().enumerate().take(20) {
                assert_eq!(&jumped.next_point(), want, "start {start}, offset {offset}");
            }
        }
    }
    #[test]
    fn sobol_sequence_restarts_after_two_to_the_32_points() {
        if let Ok(period) = usize::try_from(1u64 << 32) {
            let mut sobol = SobolSequence::new(4).unwrap();
            sobol.skip_to(period);
            assert!(sobol.next_point().iter().all(|&v| v == 0.0));
            assert_eq!(sobol.next_point(), vec![0.5; 4]);
        }
    }
    #[test]
    fn sobol_integrates_a_product_integrand_with_low_error() {
        // Integral of prod_i (1 + (x_i - 1/2)) over the unit cube is exactly 1.  Identical or
        // correlated columns bias the product terms (E[(x - 1/2)^2] = 1/12 per pair).
        let dims = 8;
        let n = 4096;
        let pts = RandomSampler::sobol(n, dims).unwrap();
        let estimate: f64 = pts
            .iter()
            .map(|p| p.iter().map(|&x| 1.0 + (x - 0.5)).product::<f64>())
            .sum::<f64>()
            / n as f64;
        assert!((estimate - 1.0).abs() < 2e-3, "estimate = {estimate}");
    }
    /// Lag-1 correlation of consecutive increments and mean-square terminal value
    /// of `n_paths` simulated fBm paths.
    fn fbm_statistics(
        hurst: f64,
        dt: f64,
        n_steps: usize,
        n_paths: usize,
        seed: u64,
    ) -> (f64, f64) {
        let fbm = FractionalBrownianMotion::new(hurst, n_steps, dt);
        let mut rng = Xoshiro256::new(seed);
        let (mut sum_xy, mut sum_xx, mut sum_yy, mut terminal) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for _ in 0..n_paths {
            let path = fbm.generate_path(&mut rng);
            let increments: Vec<f64> = path.windows(2).map(|w| w[1] - w[0]).collect();
            for pair in increments.windows(2) {
                sum_xy += pair[0] * pair[1];
                sum_xx += pair[0] * pair[0];
                sum_yy += pair[1] * pair[1];
            }
            terminal += path[n_steps] * path[n_steps];
        }
        (sum_xy / (sum_xx * sum_yy).sqrt(), terminal / n_paths as f64)
    }
    #[test]
    fn fbm_increment_correlation_follows_the_hurst_exponent() {
        // Corr(dB_k, dB_{k+1}) = 2^(2H-1) - 1: positive (persistent) for H > 1/2,
        // negative (anti-persistent) for H < 1/2, zero for plain Brownian motion.
        for &(hurst, seed) in &[(0.8_f64, 41_u64), (0.3, 42), (0.5, 43), (0.65, 44)] {
            let expected = 2.0_f64.powf(2.0 * hurst - 1.0) - 1.0;
            let (rho, _) = fbm_statistics(hurst, 0.1, 64, 2000, seed);
            assert!(
                (rho - expected).abs() < 0.05,
                "H={hurst}: lag-1 correlation {rho}, expected {expected}"
            );
        }
    }
    #[test]
    fn fbm_terminal_variance_matches_theory() {
        // Var[B_H(n dt)] = (n dt)^(2H); independent increments would give n dt^(2H).
        for &(hurst, seed) in &[(0.8_f64, 51_u64), (0.3, 52), (0.5, 53)] {
            let (n_steps, dt) = (64usize, 0.1_f64);
            let fbm = FractionalBrownianMotion::new(hurst, n_steps, dt);
            let expected = fbm.theoretical_variance(n_steps as f64 * dt);
            let (_, mean_square) = fbm_statistics(hurst, dt, n_steps, 4000, seed);
            assert!(
                (mean_square / expected - 1.0).abs() < 0.1,
                "H={hurst}: E[B^2] = {mean_square}, expected {expected}"
            );
        }
    }
    #[test]
    fn fbm_variance_grows_like_t_to_the_2h() {
        let hurst = 0.8;
        let (n_steps, dt) = (64usize, 1.0_f64);
        let fbm = FractionalBrownianMotion::new(hurst, n_steps, dt);
        let mut rng = Xoshiro256::new(61);
        let (mut quarter, mut full) = (0.0f64, 0.0f64);
        let n_paths = 4000;
        for _ in 0..n_paths {
            let path = fbm.generate_path(&mut rng);
            quarter += path[n_steps / 4].powi(2);
            full += path[n_steps].powi(2);
        }
        // (t_end / t_quarter)^(2H) = 4^1.6 = 9.19; independent increments would give 4.
        let ratio = full / quarter;
        let expected = 4.0_f64.powf(2.0 * hurst);
        assert!(
            (ratio / expected - 1.0).abs() < 0.15,
            "ratio {ratio}, expected {expected}"
        );
    }
    #[test]
    fn fbm_increments_have_fgn_autocovariance() {
        let (hurst, n_steps, dt) = (0.7_f64, 32usize, 0.5_f64);
        let fbm = FractionalBrownianMotion::new(hurst, n_steps, dt);
        let mut rng = Xoshiro256::new(71);
        let n_paths = 4000;
        let lags = [0usize, 1, 2, 5];
        let mut acc = [0.0f64; 4];
        let mut counts = [0usize; 4];
        for _ in 0..n_paths {
            let incr = fbm.generate_increments(&mut rng);
            for (slot, &lag) in lags.iter().enumerate() {
                for k in 0..(n_steps - lag) {
                    acc[slot] += incr[k] * incr[k + lag];
                    counts[slot] += 1;
                }
            }
        }
        for (slot, &lag) in lags.iter().enumerate() {
            let sample = acc[slot] / counts[slot] as f64;
            let expected =
                dt.powf(2.0 * hurst) * FractionalBrownianMotion::fgn_autocovariance(hurst, lag);
            assert!(
                (sample - expected).abs() < 0.06 * dt.powf(2.0 * hurst),
                "lag {lag}: sample {sample}, expected {expected}"
            );
        }
    }
    #[test]
    fn fbm_hosking_is_the_cholesky_transform_of_the_exact_covariance() {
        let m = 24usize;
        for &hurst in &[0.2_f64, 0.5, 0.65, 0.8, 0.95] {
            let mut rng = Xoshiro256::new(900 + (hurst * 100.0) as u64);
            let normals: Vec<f64> = (0..m).map(|_| box_muller(&mut rng).0).collect();
            let hosking = FractionalBrownianMotion::hosking_fgn(hurst, &normals);
            let mut cov = vec![0.0f64; m * m];
            for i in 0..m {
                for j in 0..m {
                    cov[i * m + j] =
                        FractionalBrownianMotion::fgn_autocovariance(hurst, i.abs_diff(j));
                }
            }
            let l = cholesky_lower(&cov, m);
            for i in 0..m {
                let expected: f64 = (0..=i).map(|j| l[i * m + j] * normals[j]).sum();
                assert!(
                    (hosking[i] - expected).abs() < 1e-9,
                    "H={hurst}, i={i}: Hosking {} vs Cholesky {expected}",
                    hosking[i]
                );
            }
        }
    }
    #[test]
    fn fbm_autocovariance_known_values() {
        assert_eq!(FractionalBrownianMotion::fgn_autocovariance(0.7, 0), 1.0);
        for &hurst in &[0.3_f64, 0.5, 0.8] {
            let gamma1 = FractionalBrownianMotion::fgn_autocovariance(hurst, 1);
            assert!((gamma1 - (2.0_f64.powf(2.0 * hurst - 1.0) - 1.0)).abs() < 1e-14);
        }
        // Brownian motion has uncorrelated increments.
        for lag in 1..6 {
            assert!(FractionalBrownianMotion::fgn_autocovariance(0.5, lag).abs() < 1e-14);
        }
    }
    #[test]
    fn fbm_path_is_deterministic_for_a_fixed_seed_and_scales_with_dt() {
        let fbm = FractionalBrownianMotion::new(0.7, 40, 0.25);
        let a = fbm.generate_path(&mut Xoshiro256::new(5));
        let b = fbm.generate_path(&mut Xoshiro256::new(5));
        let c = fbm.generate_path(&mut Xoshiro256::new(6));
        assert_eq!(a, b);
        assert_ne!(a, c);
        // Same normals, dt scaled by 4: the path scales by 4^H.
        let wide = FractionalBrownianMotion::new(0.7, 40, 1.0);
        let d = wide.generate_path(&mut Xoshiro256::new(5));
        for (x, y) in a.iter().zip(d.iter()) {
            assert!((y - x * 4.0_f64.powf(0.7)).abs() < 1e-9);
        }
    }
    #[test]
    fn fbm_degenerate_step_counts_and_out_of_range_hurst() {
        let zero = FractionalBrownianMotion::new(0.6, 0, 0.1);
        assert_eq!(zero.generate_path(&mut Xoshiro256::new(1)), vec![0.0]);
        let one = FractionalBrownianMotion::new(0.6, 1, 0.1);
        assert_eq!(one.generate_path(&mut Xoshiro256::new(1)).len(), 2);
        for &hurst in &[-0.5_f64, 0.0, 1.0, 2.5] {
            let path = FractionalBrownianMotion::new(hurst, 50, 0.1)
                .generate_path(&mut Xoshiro256::new(2));
            assert!(
                path.iter().all(|v| v.is_finite()),
                "H={hurst} produced non-finite values"
            );
        }
    }
}
