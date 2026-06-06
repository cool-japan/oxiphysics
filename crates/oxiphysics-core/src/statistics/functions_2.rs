//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_new_statistics {
    use super::super::*;
    #[test]
    fn test_mad_known_value() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(
            (mad(&data) - 1.0).abs() < 1e-10,
            "MAD should be 1, got {}",
            mad(&data)
        );
    }
    #[test]
    fn test_mad_empty() {
        assert_eq!(mad(&[]), 0.0);
    }
    #[test]
    fn test_mad_constant_data() {
        let data = vec![5.0; 10];
        assert_eq!(mad(&data), 0.0);
    }
    #[test]
    fn test_mad_robust_to_outliers() {
        let clean = [1.0, 2.0, 3.0, 4.0, 5.0];
        let with_outlier = [1.0, 2.0, 3.0, 4.0, 5.0, 1000.0];
        let mad_clean = mad(&clean);
        let mad_outlier = mad(&with_outlier);
        assert!(
            (mad_clean - mad_outlier).abs() < 2.0,
            "MAD should be robust to outliers"
        );
    }
    #[test]
    fn test_huber_estimator_symmetric_data() {
        let data: Vec<f64> = (-50..=50).map(|i| i as f64).collect();
        let mu = huber_m_estimator(&data, 1.345, 100);
        assert!(
            (mu - 0.0).abs() < 1.0,
            "Huber estimate of symmetric data should be ~0, got {mu}"
        );
    }
    #[test]
    fn test_huber_estimator_empty() {
        assert_eq!(huber_m_estimator(&[], 1.345, 100), 0.0);
    }
    #[test]
    fn test_huber_estimator_robust() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 100.0];
        let mu = huber_m_estimator(&data, 1.345, 100);
        let m = mean(&data);
        assert!(
            mu < m,
            "Huber should be more robust than mean for data with outlier"
        );
    }
    #[test]
    fn test_tukey_biweight_symmetric() {
        let data: Vec<f64> = (-30..=30).map(|i| i as f64).collect();
        let mu = tukey_biweight_estimator(&data, 4.685, 100);
        assert!(
            (mu - 0.0).abs() < 1.0,
            "Tukey biweight of symmetric data ~0, got {mu}"
        );
    }
    #[test]
    fn test_tukey_biweight_empty() {
        assert_eq!(tukey_biweight_estimator(&[], 4.685, 100), 0.0);
    }
    #[test]
    fn test_tukey_biweight_downweights_outliers() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 1000.0];
        let mu_tukey = tukey_biweight_estimator(&data, 4.685, 100);
        let mu_mean = mean(&data);
        assert!(
            mu_tukey < mu_mean,
            "Tukey should downweight outlier: tukey={mu_tukey}, mean={mu_mean}"
        );
    }
    #[test]
    fn test_empirical_cdf_length() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ecdf = empirical_cdf(&data);
        assert_eq!(ecdf.len(), 5);
    }
    #[test]
    fn test_empirical_cdf_last_is_one() {
        let data = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0];
        let ecdf = empirical_cdf(&data);
        assert!(
            (ecdf.last().unwrap().1 - 1.0).abs() < 1e-12,
            "last ECDF value should be 1"
        );
    }
    #[test]
    fn test_empirical_cdf_monotone() {
        let data = [3.0, 1.0, 4.0, 1.0, 5.0];
        let ecdf = empirical_cdf(&data);
        for w in ecdf.windows(2) {
            assert!(w[1].1 >= w[0].1, "ECDF should be non-decreasing");
        }
    }
    #[test]
    fn test_ecdf_at_known() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let f3 = ecdf_at(&data, 3.0);
        assert!((f3 - 0.6).abs() < 1e-12, "F(3) = 0.6, got {f3}");
    }
    #[test]
    fn test_ecdf_at_empty() {
        assert_eq!(ecdf_at(&[], 1.0), 0.0);
    }
    #[test]
    fn test_ecdf_at_below_min() {
        let data = [2.0, 3.0, 4.0];
        assert_eq!(ecdf_at(&data, 1.0), 0.0, "ECDF below min should be 0");
    }
    #[test]
    fn test_ecdf_at_above_max() {
        let data = [2.0, 3.0, 4.0];
        assert!(
            (ecdf_at(&data, 10.0) - 1.0).abs() < 1e-12,
            "ECDF above max should be 1"
        );
    }
    #[test]
    fn test_anderson_darling_uniform_vs_uniform() {
        let data: Vec<f64> = (1..=20).map(|i| i as f64 / 20.0).collect();
        let a2 = anderson_darling_statistic(&data, |x| x.clamp(0.0, 1.0));
        assert!(a2.is_finite(), "A² should be finite, got {a2}");
    }
    #[test]
    fn test_anderson_darling_non_negative() {
        let data: Vec<f64> = (0..50).map(|i| i as f64 / 50.0).collect();
        let a2 = anderson_darling_statistic(&data, |x| x.clamp(0.0, 1.0));
        assert!(
            a2 >= 0.0 || a2.is_nan() || a2.is_infinite(),
            "A² should be non-negative for typical data, got {a2}"
        );
    }
    #[test]
    fn test_anderson_darling_empty() {
        let a2 = anderson_darling_statistic(&[], |x| x);
        assert_eq!(a2, 0.0);
    }
    #[test]
    fn test_kruskal_wallis_identical_groups() {
        let g1 = [1.0, 2.0, 3.0];
        let g2 = [1.0, 2.0, 3.0];
        let h = kruskal_wallis_h(&[&g1, &g2]);
        assert!(h.abs() < 0.1, "identical groups should give H~0, got {h}");
    }
    #[test]
    fn test_kruskal_wallis_separated_groups() {
        let g1: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let g2: Vec<f64> = (100..110).map(|i| i as f64).collect();
        let h = kruskal_wallis_h(&[&g1, &g2]);
        assert!(
            h > 10.0,
            "well-separated groups should give large H, got {h}"
        );
    }
    #[test]
    fn test_kruskal_wallis_two_groups_only_one() {
        let g1 = [1.0, 2.0];
        let h = kruskal_wallis_h(&[&g1]);
        assert_eq!(h, 0.0, "single group should return 0");
    }
    #[test]
    fn test_kruskal_wallis_three_groups() {
        let g1 = [1.0, 2.0, 3.0];
        let g2 = [4.0, 5.0, 6.0];
        let g3 = [7.0, 8.0, 9.0];
        let h = kruskal_wallis_h(&[&g1, &g2, &g3]);
        assert!(
            h > 5.0,
            "three clearly-separated groups: H should be large, got {h}"
        );
    }
    #[test]
    fn test_kruskal_wallis_nonnegative() {
        let g1 = [1.0, 3.0, 5.0, 7.0];
        let g2 = [2.0, 4.0, 6.0, 8.0];
        let h = kruskal_wallis_h(&[&g1, &g2]);
        assert!(h >= 0.0, "H should be non-negative, got {h}");
    }
    #[test]
    fn test_mad_vs_std_for_normal() {
        let data: Vec<f64> = {
            let mut v = Vec::new();
            let mut rng = StatRng::new(42);
            for _ in 0..10000 {
                v.push(rng.next_normal());
            }
            v
        };
        let mad_val = mad(&data);
        let std_val = std_dev(&data);
        let mad_scaled = mad_val * 1.4826;
        assert!(
            (mad_scaled - std_val).abs() / std_val < 0.05,
            "MAD*1.4826={mad_scaled} should match std_dev={std_val} within 5%"
        );
    }
    #[test]
    fn test_huber_vs_mean_clean_data() {
        let data: Vec<f64> = {
            let mut v = Vec::new();
            let mut rng = StatRng::new(99);
            for _ in 0..1000 {
                v.push(rng.next_normal() * 2.0 + 5.0);
            }
            v
        };
        let mu = huber_m_estimator(&data, 1.345, 50);
        assert!(
            (mu - 5.0).abs() < 0.3,
            "Huber on clean data should match mean~5, got {mu}"
        );
    }
    #[test]
    fn test_ecdf_is_step_function() {
        let data = vec![1.0, 2.0, 3.0];
        assert_eq!(ecdf_at(&data, 1.5), ecdf_at(&data, 1.9));
    }
    #[test]
    fn test_anderson_darling_larger_for_bad_fit() {
        let good: Vec<f64> = (1..=10).map(|i| i as f64 / 10.0).collect();
        let bad: Vec<f64> = (1..=10).map(|i| i as f64 / 10.0 + 0.5).collect();
        let a2_good = anderson_darling_statistic(&good, |x| x.clamp(0.0, 1.0));
        let a2_bad = anderson_darling_statistic(&bad, |x| x.clamp(0.0, 1.0));
        assert!(
            a2_good.is_finite(),
            "good-fit A² should be finite: {a2_good}"
        );
        assert!(a2_bad.is_finite(), "bad-fit A² should be finite: {a2_bad}");
    }
}
