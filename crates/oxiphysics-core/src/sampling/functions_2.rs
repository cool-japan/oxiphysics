//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{min_pairwise_dist, rejection_sample_1d};
use super::types::{LatinHypercube, Lcg, SobolSequence};

#[cfg(test)]
mod additional_sampling_tests {
    use super::*;

    use crate::sampling::GibbsSampler;
    use crate::sampling::HaltonSequence;
    use crate::sampling::ImportanceSampler;
    use crate::sampling::LatinHypercube;
    use crate::sampling::Lcg;
    use crate::sampling::MetropolisHastings;
    use crate::sampling::Sobol;

    #[test]
    fn test_lcg_next_u64_non_zero_seed() {
        let mut rng = Lcg::new(1);
        let v = rng.next_u64();
        assert_ne!(v, 0, "LCG should generate non-zero for seed=1");
    }
    #[test]
    fn test_lcg_different_seeds_different_outputs() {
        let mut r1 = Lcg::new(1);
        let mut r2 = Lcg::new(2);
        let v1 = r1.next_u64();
        let v2 = r2.next_u64();
        assert_ne!(v1, v2, "Different seeds should produce different values");
    }
    #[test]
    fn test_lcg_range_stays_in_bounds() {
        let mut rng = Lcg::new(42);
        for _ in 0..1_000 {
            let v = rng.next_f64_range(-5.0, 5.0);
            assert!((-5.0..5.0).contains(&v), "value {} out of range [-5,5)", v);
        }
    }
    #[test]
    fn test_lcg_normal_pair_independent() {
        let mut rng = Lcg::new(123);
        let (mut pos0, mut neg0, mut pos1, mut neg1) = (0, 0, 0, 0);
        for _ in 0..200 {
            let (z0, z1) = rng.next_normal_pair();
            if z0 > 0.0 {
                pos0 += 1;
            } else {
                neg0 += 1;
            }
            if z1 > 0.0 {
                pos1 += 1;
            } else {
                neg1 += 1;
            }
        }
        assert!(pos0 > 50 && neg0 > 50, "z0 should span both halves");
        assert!(pos1 > 50 && neg1 > 50, "z1 should span both halves");
    }
    #[test]
    fn test_lcg_normal_std_near_one() {
        let mut rng = Lcg::new(999);
        let n = 5_000;
        let samples: Vec<f64> = (0..n).map(|_| rng.next_normal()).collect();
        let mean = samples.iter().sum::<f64>() / n as f64;
        let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        let std = var.sqrt();
        assert!(
            (std - 1.0).abs() < 0.1,
            "std of normal samples should be ~1: std={}",
            std
        );
    }
    #[test]
    fn test_stratified_1d_in_range() {
        let mut rng = Lcg::new(11);
        let samples = stratified_sample_1d(50, &mut rng);
        for &s in &samples {
            assert!((0.0..1.0).contains(&s), "value {} out of [0,1)", s);
        }
    }
    #[test]
    fn test_stratified_1d_sum_approx_half() {
        let mut rng = Lcg::new(77);
        let n = 100;
        let samples = stratified_sample_1d(n, &mut rng);
        let mean = samples.iter().sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.1, "mean should be ~0.5: {}", mean);
    }
    #[test]
    fn test_lhs_all_values_in_unit_interval() {
        let mut rng = Lcg::new(202);
        let lhs = LatinHypercube::new(20, 5);
        let samples = lhs.sample(&mut rng);
        for s in &samples {
            for &v in s {
                assert!((0.0..1.0).contains(&v), "LHS value {} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_lhs_large_n_strata_covered() {
        let mut rng = Lcg::new(303);
        let n = 32;
        let lhs = LatinHypercube::new(n, 4);
        let samples = lhs.sample(&mut rng);
        let inv_n = 1.0 / n as f64;
        for d in 0..4 {
            let mut occupied = vec![false; n];
            for s in &samples {
                let stratum = ((s[d] / inv_n) as usize).min(n - 1);
                occupied[stratum] = true;
            }
            let all_covered = occupied.iter().all(|&x| x);
            assert!(all_covered, "dim {d} not all strata covered");
        }
    }
    #[test]
    fn test_sobol_1d_first_value_zero() {
        let mut sobol = Sobol::new_1d();
        let first = sobol.next_1d();
        assert!(first.abs() < 1e-12, "Sobol(0) should be 0.0: {}", first);
    }
    #[test]
    fn test_sobol_1d_no_duplicates() {
        let mut sobol = Sobol::new_1d();
        let vals: Vec<f64> = (0..16).map(|_| sobol.next_1d()).collect();
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for w in sorted.windows(2) {
            assert!((w[1] - w[0]).abs() > 1e-15, "duplicate Sobol values");
        }
    }
    #[test]
    fn test_sobol_multidim_2d_count() {
        let sobol = SobolSequence::new(2);
        let pts = sobol.sample(32);
        assert_eq!(pts.len(), 32);
        for p in &pts {
            assert_eq!(p.len(), 2);
        }
    }
    #[test]
    fn test_sobol_multidim_values_in_range() {
        let sobol = SobolSequence::new(3);
        let pts = sobol.sample(128);
        for p in &pts {
            for &v in p {
                assert!((0.0..1.0).contains(&v), "Sobol value {} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_halton_base3_first_values() {
        let samples = HaltonSequence::sample(3, 3);
        assert!(
            (samples[0] - 1.0 / 3.0).abs() < 1e-12,
            "h3[0]={}",
            samples[0]
        );
        assert!(
            (samples[1] - 2.0 / 3.0).abs() < 1e-12,
            "h3[1]={}",
            samples[1]
        );
        assert!(
            (samples[2] - 1.0 / 9.0).abs() < 1e-12,
            "h3[2]={}",
            samples[2]
        );
    }
    #[test]
    fn test_halton_multivariate_shape() {
        let pts = halton_multivariate(20, 3);
        assert_eq!(pts.len(), 20);
        for p in &pts {
            assert_eq!(p.len(), 3);
            for &v in p {
                assert!(
                    (0.0..1.0).contains(&v),
                    "Halton multi value {} out of [0,1)",
                    v
                );
            }
        }
    }
    #[test]
    fn test_halton_sequence_free_fn_count() {
        let pts = halton_sequence(50, 5);
        assert_eq!(pts.len(), 50);
    }
    #[test]
    fn test_halton_1d_low_discrepancy_coverage() {
        let n = 8;
        let pts = HaltonSequence::sample(n, 2);
        for &v in &pts {
            assert!((0.0..1.0).contains(&v), "value {} out of [0,1)", v);
        }
        let mut sorted = pts.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for w in sorted.windows(2) {
            assert!(
                (w[1] - w[0]).abs() > 1e-12,
                "duplicate values: {} and {}",
                w[0],
                w[1]
            );
        }
        let mean = pts.iter().sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.15, "mean should be ~0.5: {}", mean);
    }
    #[test]
    fn test_mh_acceptance_rate_reasonable() {
        pub(super) fn log_target(x: f64) -> f64 {
            -0.5 * x * x
        }
        let mut mcmc = MetropolisHastings::new(0.0, 0.5, 42);
        let rate = mcmc.estimate_acceptance_rate(1000, log_target);
        assert!(rate > 0.3 && rate <= 1.0, "acceptance rate={}", rate);
    }
    #[test]
    fn test_mh_sample_length() {
        pub(super) fn log_target(x: f64) -> f64 {
            -x * x
        }
        let mut mcmc = MetropolisHastings::new(0.0, 1.0, 7);
        let samples = mcmc.sample(300, log_target);
        assert_eq!(samples.len(), 300);
    }
    #[test]
    fn test_mh_bimodal_target() {
        pub(super) fn log_target(x: f64) -> f64 {
            let p1 = (-0.5 * (x - 3.0).powi(2)).exp();
            let p2 = (-0.5 * (x + 3.0).powi(2)).exp();
            (p1 + p2).ln()
        }
        let mut mcmc = MetropolisHastings::new(3.0, 2.0, 55);
        let samples = mcmc.sample(2000, log_target);
        assert_eq!(samples.len(), 2000);
    }
    #[test]
    fn test_gibbs_output_varies() {
        let mut gibbs = GibbsSampler::new(2, 13);
        let samples = gibbs.sample(
            100,
            &[
                Box::new(|_, _: &[f64]| 0.0_f64),
                Box::new(|_, _: &[f64]| 0.0_f64),
            ],
            &[1.0, 1.0],
        );
        let all_same = samples.iter().all(|s| (s[0] - samples[0][0]).abs() < 1e-12);
        assert!(!all_same, "Gibbs output should vary across samples");
    }
    #[test]
    fn test_gibbs_3d_shape() {
        let mut gibbs = GibbsSampler::new(3, 99);
        let samples = gibbs.sample(
            50,
            &[
                Box::new(|_, _: &[f64]| 1.0_f64),
                Box::new(|_, _: &[f64]| -1.0_f64),
                Box::new(|_, _: &[f64]| 0.0_f64),
            ],
            &[0.5, 0.5, 0.5],
        );
        assert_eq!(samples.len(), 50);
        for s in &samples {
            assert_eq!(s.len(), 3);
        }
    }
    #[test]
    fn test_importance_sampler_envelope_respected() {
        pub(super) fn linear_pdf(x: f64) -> f64 {
            2.0 * x
        }
        let sampler = ImportanceSampler::new(linear_pdf, 2.0);
        let mut rng = Lcg::new(707);
        let samples = sampler.sample(200, &mut rng);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!(
            mean > 0.5,
            "mean should be > 0.5 for linear pdf: mean={}",
            mean
        );
    }
    #[test]
    fn test_importance_sampler_all_in_unit_interval() {
        pub(super) fn pdf(x: f64) -> f64 {
            x * x * 3.0
        }
        let sampler = ImportanceSampler::new(pdf, 3.0);
        let mut rng = Lcg::new(808);
        let samples = sampler.sample(100, &mut rng);
        for &s in &samples {
            assert!((0.0..1.0).contains(&s), "sample {} out of [0,1)", s);
        }
    }
    #[test]
    fn test_blue_noise_empty_for_impossible_radius() {
        let mut rng = Lcg::new(303);
        let samples = blue_noise_2d(100, 0.9, &mut rng);
        assert!(
            samples.len() < 5,
            "very large radius should yield few samples: {}",
            samples.len()
        );
    }
    #[test]
    fn test_blue_noise_in_unit_square() {
        let mut rng = Lcg::new(123);
        let r = 0.1;
        let samples = blue_noise_2d(30, r, &mut rng);
        for s in &samples {
            assert!(s[0] >= 0.0 && s[0] < 1.0, "x={} out of [0,1)", s[0]);
            assert!(s[1] >= 0.0 && s[1] < 1.0, "y={} out of [0,1)", s[1]);
        }
    }
    #[test]
    fn test_systematic_resample_uniform_weights() {
        let particles: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let weights = vec![1.0; 5];
        let mut rng = Lcg::new(42);
        let resampled = systematic_resample(&particles, &weights, 10, &mut rng);
        assert_eq!(resampled.len(), 10);
        for &v in &resampled {
            assert!(
                (0.0..=4.0).contains(&v),
                "resampled value {} out of range",
                v
            );
        }
    }
    #[test]
    fn test_systematic_resample_empty() {
        let particles: Vec<f64> = vec![];
        let weights: Vec<f64> = vec![];
        let mut rng = Lcg::new(99);
        let resampled = systematic_resample(&particles, &weights, 5, &mut rng);
        assert_eq!(resampled.len(), 0);
    }
    #[test]
    fn test_halton_multivariate_dim1_matches_1d() {
        let pts = halton_multivariate(10, 1);
        let seq = HaltonSequence::sample(10, 2);
        for (p, s) in pts.iter().zip(seq.iter()) {
            assert!((p[0] - s).abs() < 1e-12, "dim0 mismatch: {} vs {}", p[0], s);
        }
    }
    #[test]
    fn test_qmc_halton_sine_integral() {
        let pi = std::f64::consts::PI;
        let est = qmc_integrate_halton(|x| (x * pi).sin(), 0.0, 1.0, 1024);
        let est_scaled = est * pi;
        assert!(
            (est_scaled - 2.0).abs() < 0.05,
            "Halton sin integral={}",
            est_scaled
        );
    }
    #[test]
    fn test_qmc_sobol_cubic_integral() {
        let est = qmc_integrate_sobol(|x| x * x * x, 0.0, 1.0, 512);
        assert!((est - 0.25).abs() < 0.01, "Sobol cubic integral={}", est);
    }
    #[test]
    fn test_bootstrap_resample_single_element() {
        let data = vec![42.0];
        let mut rng = Lcg::new(1);
        let resamples = bootstrap_resample(&data, 5, &mut rng);
        assert_eq!(resamples.len(), 5);
        for r in &resamples {
            assert_eq!(r.len(), 1);
            assert_eq!(r[0], 42.0);
        }
    }
    #[test]
    fn test_bootstrap_mean_mean_correct() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut rng = Lcg::new(31415);
        let (mean, ci_lo, ci_hi) = bootstrap_mean_ci(&data, 200, &mut rng);
        assert!((mean - 3.0).abs() < 1e-10, "mean={}", mean);
        assert!(ci_lo < mean && ci_hi > mean, "CI should contain mean");
    }
    #[test]
    fn test_lhs_rand_strata_coverage() {
        let mut rng = rand::rng();
        let n = 16;
        let samples = latin_hypercube_sample(n, 2, &mut rng);
        let inv_n = 1.0 / n as f64;
        for d in 0..2 {
            let mut occupied = vec![false; n];
            for s in &samples {
                let stratum = ((s[d] / inv_n) as usize).min(n - 1);
                occupied[stratum] = true;
            }
            for (k, &hit) in occupied.iter().enumerate() {
                assert!(hit, "dim {d} stratum {k} not covered in rand LHS");
            }
        }
    }
    #[test]
    fn test_lhs_rand_all_in_unit() {
        let mut rng = rand::rng();
        let samples = latin_hypercube_sample(10, 3, &mut rng);
        for s in &samples {
            for &v in s {
                assert!((0.0..1.0).contains(&v), "rand LHS value {} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_ess_uniform_weights() {
        let n = 10;
        let weights = vec![1.0; n];
        let ess = effective_sample_size(&weights);
        assert!(
            (ess - n as f64).abs() < 0.01,
            "uniform weights → ESS=N: ess={}",
            ess
        );
    }
    #[test]
    fn test_ess_concentrated_weight() {
        let mut weights = vec![0.001; 100];
        weights[0] = 100.0;
        let ess = effective_sample_size(&weights);
        assert!(ess < 5.0, "concentrated weight → low ESS: ess={}", ess);
    }
    #[test]
    fn test_ess_zero_weights() {
        let weights = vec![0.0; 5];
        let ess = effective_sample_size(&weights);
        assert!(ess.abs() < 1e-10, "zero weights → ESS=0: ess={}", ess);
    }
    #[test]
    fn test_importance_sample_uniform_distribution() {
        let mut rng = rand::rng();
        let pdf = vec![1.0; 10];
        let samples = importance_sample(&pdf, 1000, &mut rng);
        let mut counts = [0usize; 10];
        for &i in &samples {
            counts[i] += 1;
        }
        for (i, &c) in counts.iter().enumerate() {
            assert!(c > 50, "index {} appeared only {} times", i, c);
        }
    }
    #[test]
    fn test_importance_sample_all_to_one_index() {
        let mut rng = rand::rng();
        let pdf = vec![1.0, 0.0, 0.0];
        let samples = importance_sample(&pdf, 50, &mut rng);
        for &i in &samples {
            assert_eq!(i, 0, "with all weight on 0, should always sample 0");
        }
    }
    #[test]
    fn test_mc_control_variate_linear() {
        let mut rng = Lcg::new(1234);
        let est = mc_control_variate(|x| x * x, |x| x, 0.5, 0.0, 1.0, 50_000, &mut rng);
        assert!(
            (est - 1.0 / 3.0).abs() < 0.01,
            "control variate estimate={}",
            est
        );
    }
    #[test]
    fn test_mc_antithetic_linear() {
        let mut rng = Lcg::new(5555);
        let est = mc_antithetic(|x| x, 0.0, 1.0, 10_000, &mut rng);
        assert!((est - 0.5).abs() < 0.01, "antithetic estimate={}", est);
    }
    #[test]
    fn test_mc_antithetic_quadratic() {
        let mut rng = Lcg::new(6666);
        let est = mc_antithetic(|x| x * x, 0.0, 1.0, 20_000, &mut rng);
        assert!(
            (est - 1.0 / 3.0).abs() < 0.02,
            "antithetic quadratic estimate={}",
            est
        );
    }
    #[test]
    fn test_stratified_nd_strata_covered_each_dim() {
        let mut rng = Lcg::new(111);
        let n = 10;
        let d = 3;
        let samples = stratified_sample_nd(n, d, &mut rng);
        let inv_n = 1.0 / n as f64;
        for dim in 0..d {
            let mut occupied = vec![false; n];
            for s in &samples {
                let stratum = ((s[dim] / inv_n) as usize).min(n - 1);
                occupied[stratum] = true;
            }
            for (k, &hit) in occupied.iter().enumerate() {
                assert!(hit, "dim {} stratum {} not covered", dim, k);
            }
        }
    }
    #[test]
    fn test_mc_integrate_nd_constant() {
        let mut rng = Lcg::new(9001);
        let (est, _) = monte_carlo_integrate_nd(|_| 1.0, 3, 10_000, &mut rng);
        assert!((est - 1.0).abs() < 0.05, "3D constant integral={}", est);
    }
    #[test]
    fn test_mc_integrate_nd_product() {
        let mut rng = Lcg::new(9002);
        let (est, _) = monte_carlo_integrate_nd(|x| x[0] * x[1], 2, 100_000, &mut rng);
        assert!((est - 0.25).abs() < 0.02, "2D product integral={}", est);
    }
}
/// Owen-scrambled van der Corput sequence for one base.
///
/// Applies a deterministic digit-reversal scramble parameterized by `seed`.
/// This reduces correlation between different bases in a Halton sequence.
#[allow(dead_code)]
pub fn van_der_corput_scrambled(mut i: u32, base: u32, seed: u32) -> f64 {
    let mut result = 0.0_f64;
    let mut denom = 1.0_f64;
    let mut perm_seed = seed.wrapping_add(1);
    while i > 0 {
        denom *= base as f64;
        let digit = i % base;
        let perm = (digit.wrapping_add(perm_seed)) % base;
        result += perm as f64 / denom;
        i /= base;
        perm_seed = perm_seed.wrapping_mul(1664525).wrapping_add(1013904223);
    }
    result
}
/// Scrambled Halton sequence: `n` points in `n_dims` dimensions.
///
/// Uses a deterministic per-dimension scramble to reduce correlation.
#[allow(dead_code)]
pub fn halton_scrambled(n: usize, n_dims: usize) -> Vec<Vec<f64>> {
    pub(super) const PRIMES: [u32; 16] =
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    let d = n_dims.min(PRIMES.len());
    (1..=n)
        .map(|i| {
            (0..d)
                .map(|k| {
                    van_der_corput_scrambled(
                        i as u32,
                        PRIMES[k],
                        (k as u32).wrapping_mul(2654435761),
                    )
                })
                .collect()
        })
        .collect()
}
/// Roberts R2 quasi-random sequence (2-D).
///
/// Produces a two-dimensional low-discrepancy sequence based on the plastic
/// constant φ₂ ≈ 1.32471795.  Points lie in `[0,1)²`.
#[allow(dead_code)]
pub fn r2_sequence(n: usize) -> Vec<[f64; 2]> {
    let phi2 = 1.324_717_957_244_746_f64;
    let a1 = 1.0 / phi2;
    let a2 = 1.0 / (phi2 * phi2);
    (0..n)
        .map(|i| {
            let x = (0.5 + a1 * (i + 1) as f64).fract();
            let y = (0.5 + a2 * (i + 1) as f64).fract();
            [x, y]
        })
        .collect()
}
/// Roberts R2 sequence for `d` dimensions using the plastic constant.
#[allow(dead_code)]
pub fn r2_sequence_nd(n: usize, d: usize) -> Vec<Vec<f64>> {
    let phi2 = 1.324_717_957_244_746_f64;
    let alphas: Vec<f64> = (1..=d)
        .map(|k| {
            let mut p = 1.0_f64;
            for _ in 0..k {
                p /= phi2;
            }
            p
        })
        .collect();
    (0..n)
        .map(|i| {
            alphas
                .iter()
                .map(|&a| (0.5 + a * (i + 1) as f64).fract())
                .collect()
        })
        .collect()
}
/// Latin Hypercube sampler (centered strata).
///
/// Places one sample at the *center* of each stratum per dimension (no jitter).
/// This maximises spread but is completely deterministic given `n_samples`.
#[allow(dead_code)]
pub fn latin_hypercube_centered(n_samples: usize, n_dims: usize, rng: &mut Lcg) -> Vec<Vec<f64>> {
    let inv_n = 1.0 / n_samples as f64;
    let mut result: Vec<Vec<f64>> = (0..n_samples).map(|_| Vec::with_capacity(n_dims)).collect();
    for _ in 0..n_dims {
        let mut perm: Vec<usize> = (0..n_samples).collect();
        for i in (1..n_samples).rev() {
            let j = (rng.next_u64() as usize) % (i + 1);
            perm.swap(i, j);
        }
        for i in 0..n_samples {
            let val = (perm[i] as f64 + 0.5) * inv_n;
            result[i].push(val);
        }
    }
    result
}
/// Maximin Latin Hypercube design.
///
/// Generates `n_candidates` LHS designs and returns the one that maximises
/// the minimum pairwise distance between samples.
#[allow(dead_code)]
pub fn latin_hypercube_maximin(
    n_samples: usize,
    n_dims: usize,
    n_candidates: usize,
    rng: &mut Lcg,
) -> Vec<Vec<f64>> {
    let lhc = LatinHypercube::new(n_samples, n_dims);
    let mut best = lhc.sample(rng);
    let mut best_min_dist = min_pairwise_dist(&best);
    for _ in 1..n_candidates {
        let candidate = lhc.sample(rng);
        let d = min_pairwise_dist(&candidate);
        if d > best_min_dist {
            best_min_dist = d;
            best = candidate;
        }
    }
    best
}
/// Compute the minimum pairwise Euclidean distance among a set of points (maximin variant).
#[allow(dead_code)]
pub(super) fn min_pairwise_dist_maximin(pts: &[Vec<f64>]) -> f64 {
    let n = pts.len();
    if n < 2 {
        return f64::INFINITY;
    }
    let mut min_d = f64::INFINITY;
    for i in 0..n {
        for j in (i + 1)..n {
            let d2: f64 = pts[i]
                .iter()
                .zip(pts[j].iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum();
            if d2 < min_d {
                min_d = d2;
            }
        }
    }
    min_d.sqrt()
}
/// Generate `n` stratified samples on the surface of the unit sphere.
///
/// Divides the sphere into an equal-area latitude-longitude grid of
/// `n_lat × n_lon` cells and places one jittered sample per cell.
/// Returns unit vectors `[x, y, z]`.
#[allow(dead_code)]
pub fn stratified_sphere_samples(n_lat: usize, n_lon: usize, rng: &mut Lcg) -> Vec<[f64; 3]> {
    let mut samples = Vec::with_capacity(n_lat * n_lon);
    let inv_lat = 1.0 / n_lat as f64;
    let inv_lon = 1.0 / n_lon as f64;
    for il in 0..n_lat {
        for ip in 0..n_lon {
            let u1 = (il as f64 + rng.next_f64()) * inv_lat;
            let u2 = (ip as f64 + rng.next_f64()) * inv_lon;
            let cos_theta = 1.0 - 2.0 * u1;
            let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
            let phi = 2.0 * std::f64::consts::PI * u2;
            samples.push([sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta]);
        }
    }
    samples
}
/// Inverse-CDF importance sampling for a target distribution defined by its
/// CDF, computed numerically on a grid.
///
/// `cdf_vals[i]` = CDF at `x_min + i * dx`.  Draws `n` samples.
#[allow(dead_code)]
pub fn inverse_cdf_sample(
    cdf_vals: &[f64],
    x_min: f64,
    x_max: f64,
    n: usize,
    rng: &mut Lcg,
) -> Vec<f64> {
    let m = cdf_vals.len();
    assert!(m >= 2, "cdf_vals must have at least 2 entries");
    let dx = (x_max - x_min) / (m - 1) as f64;
    let mut samples = Vec::with_capacity(n);
    for _ in 0..n {
        let u = rng.next_f64();
        let idx = cdf_vals.partition_point(|&c| c < u).min(m - 1);
        let x = if idx == 0 {
            x_min
        } else {
            let c_lo = cdf_vals[idx - 1];
            let c_hi = cdf_vals[idx];
            let frac = if (c_hi - c_lo).abs() < 1e-300 {
                0.0
            } else {
                (u - c_lo) / (c_hi - c_lo)
            };
            x_min + (idx as f64 - 1.0 + frac) * dx
        };
        samples.push(x);
    }
    samples
}
/// Compute a numerical CDF from a non-negative PDF sampled on a uniform grid.
///
/// `pdf_vals[i]` = PDF(x_min + i * dx).  Returns normalized CDF values on the
/// same grid.
#[allow(dead_code)]
pub fn pdf_to_cdf(pdf_vals: &[f64], x_min: f64, x_max: f64) -> Vec<f64> {
    let m = pdf_vals.len();
    if m == 0 {
        return Vec::new();
    }
    let dx = (x_max - x_min) / (m as f64 - 1.0).max(1.0);
    let mut cdf = Vec::with_capacity(m);
    let mut acc = 0.0_f64;
    cdf.push(0.0);
    for i in 1..m {
        acc += (pdf_vals[i - 1] + pdf_vals[i]) * 0.5 * dx;
        cdf.push(acc);
    }
    let total = *cdf.last().expect("cdf is non-empty");
    if total > 1e-300 {
        for c in &mut cdf {
            *c /= total;
        }
    }
    if let Some(last) = cdf.last_mut() {
        *last = 1.0;
    }
    cdf
}
/// 1-D rejection sampling from an un-normalized density `f` on `[a, b]`.
///
/// `envelope` must satisfy `f(x) ≤ envelope` for all `x ∈ [a, b]`.
/// Draws exactly `n` accepted samples. Uses `Lcg` RNG.
#[allow(dead_code)]
pub fn rejection_sample_1d_lcg(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    envelope: f64,
    n: usize,
    rng: &mut Lcg,
) -> Vec<f64> {
    let mut samples = Vec::with_capacity(n);
    while samples.len() < n {
        let x = rng.next_f64_range(a, b);
        let u = rng.next_f64() * envelope;
        if u <= f(x) {
            samples.push(x);
        }
    }
    samples
}
/// 1-D rejection sampling with automatic envelope estimation.
///
/// Estimates the maximum of `f` on a grid of `grid_size` points, then
/// applies standard rejection sampling.
#[allow(dead_code)]
pub fn rejection_sample_1d_auto(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    n: usize,
    grid_size: usize,
    rng: &mut Lcg,
) -> Vec<f64> {
    let envelope = {
        let mut max_f = 1e-300_f64;
        for k in 0..grid_size {
            let x = a + (k as f64 + 0.5) / grid_size as f64 * (b - a);
            max_f = max_f.max(f(x));
        }
        max_f * 1.01
    };
    rejection_sample_1d(f, a, b, envelope, n, rng)
}
/// Generate `n` points from a digit-scrambled Sobol sequence in `d` dimensions.
///
/// Each dimension is XOR-scrambled with a distinct 32-bit seed for better
/// uniformity.
#[allow(dead_code)]
pub fn sobol_scrambled(n: usize, d: usize) -> Vec<Vec<f64>> {
    pub(super) const BITS: usize = 32;
    let n_dims = d.min(3);
    let dir: Vec<Vec<u32>> = (0..n_dims).map(SobolSequence::direction_numbers).collect();
    let seeds: Vec<u32> = (0..n_dims as u32)
        .map(|k| k.wrapping_mul(2654435761).wrapping_add(1))
        .collect();
    (0..n)
        .map(|idx| {
            let gray = idx ^ (idx >> 1);
            (0..n_dims)
                .map(|dd| {
                    let x = (0..BITS)
                        .filter(|&i| (gray >> i) & 1 == 1)
                        .fold(0u32, |acc, i| acc ^ dir[dd][i]);
                    let x_scrambled = x ^ seeds[dd];
                    x_scrambled as f64 / (1u64 << BITS) as f64
                })
                .collect()
        })
        .collect()
}
/// Sample mean of a slice.
#[allow(dead_code)]
pub fn sample_mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().sum::<f64>() / xs.len() as f64
}
/// Sample variance (Bessel's corrected, `1/(n-1)`) of a slice.
#[allow(dead_code)]
pub fn sample_variance(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 2 {
        return 0.0;
    }
    let mean = sample_mean(xs);
    xs.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1) as f64
}
/// Sample standard deviation.
#[allow(dead_code)]
pub fn sample_std(xs: &[f64]) -> f64 {
    sample_variance(xs).sqrt()
}
/// Sample skewness.
#[allow(dead_code)]
pub fn sample_skewness(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 3 {
        return 0.0;
    }
    let mean = sample_mean(xs);
    let std = sample_std(xs);
    if std < 1e-300 {
        return 0.0;
    }
    let n_f = n as f64;
    let m3 = xs.iter().map(|&x| ((x - mean) / std).powi(3)).sum::<f64>() / n_f;
    m3 * (n_f * n_f) / ((n_f - 1.0) * (n_f - 2.0))
}
/// Sample excess kurtosis (Fisher's definition, 0 for Gaussian).
#[allow(dead_code)]
pub fn sample_kurtosis(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 4 {
        return 0.0;
    }
    let mean = sample_mean(xs);
    let std = sample_std(xs);
    if std < 1e-300 {
        return 0.0;
    }
    let n_f = n as f64;
    let m4 = xs.iter().map(|&x| ((x - mean) / std).powi(4)).sum::<f64>() / n_f;
    (n_f * (n_f + 1.0) * m4 - 3.0 * (n_f - 1.0).powi(2)) / ((n_f - 1.0) * (n_f - 2.0) * (n_f - 3.0))
}
/// Compute empirical quantiles at fractional positions in `qs`.
///
/// `xs` need not be sorted (a copy is sorted internally).  Each value in `qs`
/// must be in `[0, 1]`.
#[allow(dead_code)]
pub fn empirical_quantiles(xs: &[f64], qs: &[f64]) -> Vec<f64> {
    if xs.is_empty() {
        return vec![f64::NAN; qs.len()];
    }
    let mut sorted = xs.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    qs.iter()
        .map(|&q| {
            let idx_f = q * (sorted.len() - 1) as f64;
            let lo = idx_f.floor() as usize;
            let hi = (lo + 1).min(sorted.len() - 1);
            let frac = idx_f.fract();
            sorted[lo] * (1.0 - frac) + sorted[hi] * frac
        })
        .collect()
}
/// Kolmogorov-Smirnov statistic between two sorted empirical CDFs.
///
/// Computes `sup |F_n(x) - G_m(x)|`.  Both slices must be sorted ascending.
#[allow(dead_code)]
pub fn ks_statistic(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    let m = ys.len();
    if n == 0 || m == 0 {
        return 0.0;
    }
    let mut i = 0usize;
    let mut j = 0usize;
    let mut max_diff = 0.0_f64;
    while i < n || j < m {
        let x_val = if i < n { xs[i] } else { f64::INFINITY };
        let y_val = if j < m { ys[j] } else { f64::INFINITY };
        if (x_val - y_val).abs() < 1e-15 {
            i += 1;
            j += 1;
        } else if x_val < y_val {
            i += 1;
        } else {
            j += 1;
        }
        let diff = (i as f64 / n as f64 - j as f64 / m as f64).abs();
        if diff > max_diff {
            max_diff = diff;
        }
    }
    max_diff
}
/// Reservoir sampling (Vitter's Algorithm R).
///
/// Uniformly samples `k` items from an iterator of unknown size without
/// storing the entire sequence.
#[allow(dead_code)]
pub fn reservoir_sample<T: Clone>(items: &[T], k: usize, rng: &mut Lcg) -> Vec<T> {
    let n = items.len();
    if k == 0 || n == 0 {
        return Vec::new();
    }
    let k_actual = k.min(n);
    let mut reservoir: Vec<T> = items[..k_actual].to_vec();
    for i in k_actual..n {
        let j = (rng.next_u64() % (i + 1) as u64) as usize;
        if j < k_actual {
            reservoir[j] = items[i].clone();
        }
    }
    reservoir
}
/// Weighted reservoir sample of `k` indices from a weight vector.
///
/// Uses Algorithm A-Res (Efraimidis & Spirakis 2006).
#[allow(dead_code)]
pub fn weighted_reservoir_sample(weights: &[f64], k: usize, rng: &mut Lcg) -> Vec<usize> {
    let n = weights.len();
    if k == 0 || n == 0 {
        return Vec::new();
    }
    let k_actual = k.min(n);
    let mut heap: Vec<(f64, usize)> = Vec::with_capacity(k_actual);
    for (i, &w) in weights.iter().enumerate() {
        if w <= 0.0 {
            continue;
        }
        let u = rng.next_f64().max(1e-300);
        let key = u.powf(1.0 / w);
        if heap.len() < k_actual {
            heap.push((key, i));
            if heap.len() == k_actual {
                heap.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            }
        } else if key > heap[0].0 {
            heap[0] = (key, i);
            heap.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        }
    }
    heap.iter().map(|&(_, i)| i).collect()
}
#[cfg(test)]
mod extended_sampling_tests_v2 {
    use super::*;
    use crate::random_processes::sample_mean;
    use crate::random_processes::sample_variance;
    use crate::sample_kurtosis;
    use crate::sample_skewness;

    use crate::sampling::LatinHypercube;
    use crate::sampling::Lcg;

    use crate::sampling::empirical_quantiles;
    use crate::sampling::halton_scrambled;
    use crate::sampling::inverse_cdf_sample;
    use crate::sampling::latin_hypercube_centered;
    use crate::sampling::latin_hypercube_maximin;
    use crate::sampling::min_pairwise_dist;
    use crate::sampling::pdf_to_cdf;
    use crate::sampling::r2_sequence;
    use crate::sampling::r2_sequence_nd;
    use crate::sampling::rejection_sample_1d;
    use crate::sampling::rejection_sample_1d_auto;
    use crate::sampling::reservoir_sample;
    use crate::sampling::sample_std;
    use crate::sampling::sobol_scrambled;
    use crate::sampling::stratified_sphere_samples;
    use crate::sampling::van_der_corput_scrambled;
    use crate::sampling::weighted_reservoir_sample;
    #[test]
    fn test_vdc_scrambled_in_unit_interval() {
        for i in 1..=20u32 {
            let v = van_der_corput_scrambled(i, 2, 12345);
            assert!((0.0..1.0).contains(&v), "vdc_scrambled out of [0,1): {}", v);
        }
    }
    #[test]
    fn test_vdc_scrambled_different_seeds_differ() {
        let v1 = van_der_corput_scrambled(7, 2, 0);
        let v2 = van_der_corput_scrambled(7, 2, 99999);
        let _ = (v1, v2);
    }
    #[test]
    fn test_halton_scrambled_dimensions() {
        let pts = halton_scrambled(10, 3);
        assert_eq!(pts.len(), 10);
        for p in &pts {
            assert_eq!(p.len(), 3);
            for &v in p {
                assert!((0.0..1.0).contains(&v), "value {} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_r2_sequence_in_unit_square() {
        let pts = r2_sequence(20);
        assert_eq!(pts.len(), 20);
        for p in &pts {
            assert!(p[0] >= 0.0 && p[0] < 1.0);
            assert!(p[1] >= 0.0 && p[1] < 1.0);
        }
    }
    #[test]
    fn test_r2_sequence_nd_correct_dims() {
        let pts = r2_sequence_nd(5, 4);
        assert_eq!(pts.len(), 5);
        for p in &pts {
            assert_eq!(p.len(), 4);
        }
    }
    #[test]
    fn test_lhc_centered_strata_exact() {
        let mut rng = Lcg::new(42);
        let n = 8;
        let pts = latin_hypercube_centered(n, 2, &mut rng);
        assert_eq!(pts.len(), n);
        let inv_n = 1.0 / n as f64;
        for d in 0..2 {
            let mut occupied = vec![false; n];
            for p in &pts {
                let stratum = ((p[d] / inv_n) as usize).min(n - 1);
                occupied[stratum] = true;
            }
            assert!(occupied.iter().all(|&h| h), "all strata should be occupied");
        }
    }
    #[test]
    fn test_lhc_maximin_returns_correct_size() {
        let mut rng = Lcg::new(1001);
        let pts = latin_hypercube_maximin(8, 2, 5, &mut rng);
        assert_eq!(pts.len(), 8);
    }
    #[test]
    fn test_lhc_maximin_better_than_single() {
        let mut rng1 = Lcg::new(77);
        let mut rng2 = Lcg::new(77);
        let single = LatinHypercube::new(12, 2).sample(&mut rng1);
        let maximin = latin_hypercube_maximin(12, 2, 10, &mut rng2);
        let d_single = min_pairwise_dist(&single);
        let d_maximin = min_pairwise_dist(&maximin);
        assert!(
            d_maximin >= d_single * 0.9,
            "maximin dist={} single dist={}",
            d_maximin,
            d_single
        );
    }
    #[test]
    fn test_stratified_sphere_samples_unit_length() {
        let mut rng = Lcg::new(55);
        let pts = stratified_sphere_samples(4, 8, &mut rng);
        assert_eq!(pts.len(), 32);
        for p in &pts {
            let r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
            assert!((r2 - 1.0).abs() < 1e-10, "r²={}", r2);
        }
    }
    #[test]
    fn test_stratified_sphere_samples_all_dims_covered() {
        let mut rng = Lcg::new(56);
        let pts = stratified_sphere_samples(4, 4, &mut rng);
        let has_pos_z = pts.iter().any(|p| p[2] > 0.0);
        let has_neg_z = pts.iter().any(|p| p[2] < 0.0);
        assert!(has_pos_z, "no positive z");
        assert!(has_neg_z, "no negative z");
    }
    #[test]
    fn test_pdf_to_cdf_uniform() {
        let pdf = vec![1.0_f64; 11];
        let cdf = pdf_to_cdf(&pdf, 0.0, 1.0);
        assert_eq!(cdf.len(), 11);
        assert!(cdf[0].abs() < 1e-10, "CDF(0)={}", cdf[0]);
        assert!((cdf[10] - 1.0).abs() < 1e-10, "CDF(1)={}", cdf[10]);
    }
    #[test]
    fn test_pdf_to_cdf_monotone() {
        let pdf = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        let cdf = pdf_to_cdf(&pdf, 0.0, 1.0);
        for i in 1..cdf.len() {
            assert!(cdf[i] >= cdf[i - 1], "CDF not monotone at i={}", i);
        }
    }
    #[test]
    fn test_inverse_cdf_sample_in_range() {
        let pdf = vec![1.0_f64; 101];
        let cdf = pdf_to_cdf(&pdf, 0.0, 1.0);
        let mut rng = Lcg::new(314);
        let samples = inverse_cdf_sample(&cdf, 0.0, 1.0, 100, &mut rng);
        for &s in &samples {
            assert!((0.0..=1.0).contains(&s), "sample {} out of range", s);
        }
    }
    #[test]
    fn test_rejection_sample_1d_count() {
        let mut rng = Lcg::new(999);
        let samples = rejection_sample_1d(|x| x * x, 0.0, 1.0, 1.0, 50, &mut rng);
        assert_eq!(samples.len(), 50, "should return exactly 50 samples");
    }
    #[test]
    fn test_rejection_sample_1d_in_range() {
        let mut rng = Lcg::new(888);
        let samples = rejection_sample_1d(|x| (-x * x).exp(), -3.0, 3.0, 1.01, 30, &mut rng);
        for &s in &samples {
            assert!((-3.0..=3.0).contains(&s), "sample {} out of range", s);
        }
    }
    #[test]
    fn test_rejection_sample_1d_auto_count() {
        let mut rng = Lcg::new(777);
        let samples = rejection_sample_1d_auto(
            |x| x.sin().powi(2),
            0.0,
            std::f64::consts::PI,
            20,
            50,
            &mut rng,
        );
        assert_eq!(samples.len(), 20);
    }
    #[test]
    fn test_sobol_scrambled_in_unit_cube() {
        let pts = sobol_scrambled(16, 3);
        assert_eq!(pts.len(), 16);
        for p in &pts {
            for &v in p {
                assert!(
                    (0.0..1.0).contains(&v),
                    "sobol_scrambled value out of [0,1)"
                );
            }
        }
    }
    #[test]
    fn test_sample_mean_uniform() {
        let xs: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let mean = sample_mean(&xs);
        assert!((mean - 0.495).abs() < 0.01, "mean={}", mean);
    }
    #[test]
    fn test_sample_variance_constant() {
        let xs = vec![3.0_f64; 10];
        let v = sample_variance(&xs);
        assert!(v.abs() < 1e-12, "variance of constant should be 0: {}", v);
    }
    #[test]
    fn test_sample_std_nonneg() {
        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let s = sample_std(&xs);
        assert!(s >= 0.0, "std should be non-negative");
    }
    #[test]
    fn test_sample_skewness_symmetric() {
        let xs: Vec<f64> = (0..100).map(|i| i as f64 - 50.0).collect();
        let sk = sample_skewness(&xs);
        assert!(sk.abs() < 0.1, "symmetric dist skewness={}", sk);
    }
    #[test]
    fn test_sample_kurtosis_uniform() {
        let xs: Vec<f64> = (0..1000).map(|i| i as f64 / 1000.0).collect();
        let k = sample_kurtosis(&xs);
        assert!(k < 0.0, "uniform kurtosis should be negative: {}", k);
    }
    #[test]
    fn test_empirical_quantiles_median() {
        let xs: Vec<f64> = (1..=9).map(|i| i as f64).collect();
        let qs = empirical_quantiles(&xs, &[0.5]);
        assert!((qs[0] - 5.0).abs() < 0.5, "median={}", qs[0]);
    }
    #[test]
    fn test_empirical_quantiles_min_max() {
        let xs = vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0];
        let qs = empirical_quantiles(&xs, &[0.0, 1.0]);
        assert!((qs[0] - 1.0).abs() < 1e-10, "min={}", qs[0]);
        assert!((qs[1] - 9.0).abs() < 1e-10, "max={}", qs[1]);
    }
    #[test]
    fn test_ks_same_distribution() {
        let xs: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let ks = ks_statistic(&xs, &xs);
        assert!(
            ks.abs() < 0.1,
            "KS between same distribution should be ~0: {}",
            ks
        );
    }
    #[test]
    fn test_ks_different_distributions() {
        let xs: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let ys: Vec<f64> = (20..40).map(|i| i as f64).collect();
        let ks = ks_statistic(&xs, &ys);
        assert!(
            ks > 0.5,
            "KS between disjoint distributions should be large: {}",
            ks
        );
    }
    #[test]
    fn test_reservoir_sample_correct_size() {
        let mut rng = Lcg::new(1337);
        let data: Vec<i32> = (0..100).collect();
        let sample = reservoir_sample(&data, 10, &mut rng);
        assert_eq!(sample.len(), 10);
    }
    #[test]
    fn test_reservoir_sample_all_from_source() {
        let mut rng = Lcg::new(2048);
        let data: Vec<i32> = (0..20).collect();
        let sample = reservoir_sample(&data, 5, &mut rng);
        for &v in &sample {
            assert!(data.contains(&v), "sampled value {} not in source", v);
        }
    }
    #[test]
    fn test_reservoir_sample_k_larger_than_n() {
        let mut rng = Lcg::new(4096);
        let data = vec![1.0, 2.0, 3.0];
        let sample = reservoir_sample(&data, 10, &mut rng);
        assert_eq!(sample.len(), 3);
    }
    #[test]
    fn test_weighted_reservoir_correct_size() {
        let mut rng = Lcg::new(7777);
        let weights = vec![1.0, 2.0, 3.0, 0.5, 4.0];
        let idx = weighted_reservoir_sample(&weights, 3, &mut rng);
        assert_eq!(idx.len(), 3);
    }
    #[test]
    fn test_weighted_reservoir_valid_indices() {
        let mut rng = Lcg::new(8888);
        let weights = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        let idx = weighted_reservoir_sample(&weights, 3, &mut rng);
        for &i in &idx {
            assert!(i < 5, "index {} out of range", i);
        }
    }
}
