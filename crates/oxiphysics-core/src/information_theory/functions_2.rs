//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::information_theory::DiscreteDistribution;
    use crate::information_theory::FisherGeometry;
    use crate::information_theory::GaussianInfoGeometry;
    use crate::information_theory::ModelSelector;
    use crate::information_theory::RateDistortion;
    use std::f64::consts::PI;
    pub(super) const EPS: f64 = 1e-9;
    #[test]
    fn test_uniform_binary_entropy_is_one_bit() {
        let p = [0.5, 0.5];
        let h = entropy_bits(&p);
        assert!((h - 1.0).abs() < EPS, "H(0.5, 0.5) = {h}");
    }
    #[test]
    fn test_certain_distribution_entropy_zero() {
        let p = [1.0, 0.0, 0.0];
        assert!(entropy_bits(&p).abs() < EPS);
    }
    #[test]
    fn test_entropy_uniform_four_symbols() {
        let p = [0.25, 0.25, 0.25, 0.25];
        let h = entropy_bits(&p);
        assert!((h - 2.0).abs() < EPS, "H = {h}");
    }
    #[test]
    fn test_entropy_nats_positive() {
        let p = [0.3, 0.3, 0.4];
        assert!(entropy_nats(&p) > 0.0);
    }
    #[test]
    fn test_entropy_empty_is_zero() {
        assert_eq!(entropy_bits(&[]), 0.0);
    }
    #[test]
    fn test_joint_entropy_independent_equals_sum() {
        let joint = vec![vec![0.25, 0.25], vec![0.25, 0.25]];
        let hxy = joint_entropy(&joint);
        assert!((hxy - 2.0).abs() < EPS, "H(X,Y) = {hxy}");
    }
    #[test]
    fn test_conditional_entropy_non_negative() {
        let joint = vec![vec![0.1, 0.4], vec![0.3, 0.2]];
        assert!(conditional_entropy(&joint) >= -EPS);
    }
    #[test]
    fn test_mutual_information_non_negative() {
        let joint = vec![vec![0.1, 0.4], vec![0.3, 0.2]];
        assert!(mutual_information(&joint) >= -EPS);
    }
    #[test]
    fn test_mutual_information_zero_for_independent() {
        let joint = vec![vec![0.25, 0.25], vec![0.25, 0.25]];
        let mi = mutual_information(&joint);
        assert!(mi.abs() < 1e-10, "MI of independent = {mi}");
    }
    #[test]
    fn test_normalised_mi_bounded() {
        let joint = vec![vec![0.1, 0.4], vec![0.3, 0.2]];
        let nmi = normalised_mutual_information(&joint);
        assert!((-EPS..=1.0 + EPS).contains(&nmi));
    }
    #[test]
    fn test_kl_divergence_self_is_zero() {
        let p = vec![0.2, 0.5, 0.3];
        assert!(kl_divergence(&p, &p).abs() < EPS);
    }
    #[test]
    fn test_kl_divergence_non_negative() {
        let p = vec![0.3, 0.4, 0.3];
        let q = vec![0.1, 0.7, 0.2];
        assert!(kl_divergence(&p, &q) >= -EPS);
    }
    #[test]
    fn test_kl_divergence_infinity_on_zero_q() {
        let p = vec![0.5, 0.5];
        let q = vec![0.0, 1.0];
        assert!(kl_divergence(&p, &q).is_infinite());
    }
    #[test]
    fn test_js_divergence_symmetric() {
        let p = vec![0.3, 0.7];
        let q = vec![0.6, 0.4];
        let js_pq = js_divergence(&p, &q);
        let js_qp = js_divergence(&q, &p);
        assert!((js_pq - js_qp).abs() < EPS);
    }
    #[test]
    fn test_js_divergence_zero_for_equal() {
        let p = vec![0.4, 0.6];
        assert!(js_divergence(&p, &p).abs() < EPS);
    }
    #[test]
    fn test_js_divergence_bounded_by_ln2() {
        let p = vec![1.0, 0.0];
        let q = vec![0.0, 1.0];
        let js = js_divergence(&p, &q);
        assert!(js <= std::f64::consts::LN_2 + EPS);
    }
    #[test]
    fn test_cross_entropy_self_equals_entropy() {
        let p = vec![0.25, 0.25, 0.25, 0.25];
        let h = entropy_bits(&p);
        let ce = cross_entropy(&p, &p);
        assert!((h - ce).abs() < EPS);
    }
    #[test]
    fn test_tv_distance_bounded() {
        let p = vec![0.3, 0.7];
        let q = vec![0.6, 0.4];
        let tv = total_variation_distance(&p, &q);
        assert!((-EPS..=1.0 + EPS).contains(&tv));
    }
    #[test]
    fn test_tv_distance_zero_for_equal() {
        let p = vec![0.2, 0.8];
        assert!(total_variation_distance(&p, &p).abs() < EPS);
    }
    #[test]
    fn test_renyi_non_negative() {
        let p = [0.2, 0.3, 0.5];
        assert!(renyi_entropy(&p, 2.0) >= -EPS);
    }
    #[test]
    fn test_renyi_alpha1_equals_shannon() {
        let p = [0.2, 0.3, 0.5];
        let r = renyi_entropy(&p, 1.0);
        let h = entropy_bits(&p);
        assert!((r - h).abs() < 1e-6);
    }
    #[test]
    fn test_renyi_ordering_shannon_vs_min() {
        let p = [0.1, 0.2, 0.7];
        let h_shannon = entropy_bits(&p);
        let h_min = min_entropy(&p);
        assert!(h_min <= h_shannon + EPS);
    }
    #[test]
    fn test_min_entropy_deterministic_is_zero() {
        let p = [1.0, 0.0, 0.0];
        assert!(min_entropy(&p).abs() < EPS);
    }
    #[test]
    fn test_min_entropy_uniform_equals_log() {
        let p = [0.25, 0.25, 0.25, 0.25];
        let h_min = min_entropy(&p);
        assert!((h_min - 2.0).abs() < EPS);
    }
    #[test]
    fn test_gaussian_differential_entropy_unit_normal() {
        let expected = 0.5 * (2.0 * PI * std::f64::consts::E).ln();
        let h = gaussian_differential_entropy(1.0);
        assert!((h - expected).abs() < EPS);
    }
    #[test]
    fn test_uniform_differential_entropy_width_e() {
        let h = uniform_differential_entropy(0.0, std::f64::consts::E);
        assert!((h - 1.0).abs() < EPS);
    }
    #[test]
    fn test_exponential_differential_entropy_rate1() {
        let h = exponential_differential_entropy(1.0);
        assert!((h - 1.0).abs() < EPS);
    }
    #[test]
    fn test_transfer_entropy_non_negative() {
        let s: Vec<f64> = (0..50).map(|i| (i as f64 * 0.1).sin()).collect();
        let t: Vec<f64> = (0..50).map(|i| (i as f64 * 0.1).cos()).collect();
        let te = transfer_entropy(&s, &t, 1);
        assert!(te >= -EPS);
    }
    #[test]
    fn test_transfer_entropy_short_series_zero() {
        let s = vec![0.1, 0.2];
        let t = vec![0.3, 0.4];
        let te = transfer_entropy(&s, &t, 3);
        assert_eq!(te, 0.0);
    }
    /// Reproducible i.i.d. fair-coin series of 0.0 / 1.0 (xorshift64*).
    fn coin_flips(n: usize, seed: u64) -> Vec<f64> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                let r = state.wrapping_mul(0x2545_F491_4F6C_DD1D);
                if r >> 63 == 1 { 1.0 } else { 0.0 }
            })
            .collect()
    }
    #[test]
    fn test_conditional_transfer_entropy_detects_xor_synergy() {
        // Y_t = X_{t-1} xor Z_{t-1} with independent fair X and Z: neither
        // source alone predicts Y, but together they determine it, so
        // T(X->Y) = 0 while T(X->Y|Z) = ln 2.
        let n = 20_000;
        let x = coin_flips(n, 0x1234_5678_9ABC_DEF1);
        let z = coin_flips(n, 0x0FED_CBA9_8765_4321);
        let mut y = vec![0.0_f64; n];
        for t in 1..n {
            y[t] = if (x[t - 1] > 0.5) != (z[t - 1] > 0.5) {
                1.0
            } else {
                0.0
            };
        }
        let cte = conditional_transfer_entropy(&x, &y, &z, 1);
        let te = transfer_entropy(&x, &y, 1);
        assert!(
            (cte - std::f64::consts::LN_2).abs() < 0.05,
            "CTE(X->Y|Z) = {cte}, expected ln 2"
        );
        assert!(te < 0.01, "TE(X->Y) = {te}, expected ~0");
        // The roles of X and Z are symmetric in this construction.
        let cte_z = conditional_transfer_entropy(&z, &y, &x, 1);
        assert!((cte_z - std::f64::consts::LN_2).abs() < 0.05, "{cte_z}");
    }
    #[test]
    fn test_conditional_transfer_entropy_removes_common_driver() {
        // Z drives Y with a one-step delay and X is a same-time copy of Z, so
        // X appears to predict Y (T(X->Y) = ln 2) only through Z.
        let n = 20_000;
        let z = coin_flips(n, 0x00C0_FFEE_D00D_F00D);
        let x = z.clone();
        let mut y = vec![0.0_f64; n];
        y[1..].copy_from_slice(&z[..n - 1]);
        let te = transfer_entropy(&x, &y, 1);
        let cte = conditional_transfer_entropy(&x, &y, &z, 1);
        assert!(
            (te - std::f64::consts::LN_2).abs() < 0.05,
            "TE(X->Y) = {te}, expected ln 2"
        );
        assert!(cte < 0.01, "CTE(X->Y|Z) = {cte}, expected ~0");
    }
    #[test]
    fn test_conditional_transfer_entropy_keeps_direct_coupling() {
        // Y_t = X_{t-1}; the independent conditioning series Z must not
        // remove a genuine causal link.
        let n = 20_000;
        let x = coin_flips(n, 0xDEAD_BEEF_CAFE_F00D);
        let z = coin_flips(n, 0xABAD_1DEA_0BAD_CAFE);
        let mut y = vec![0.0_f64; n];
        y[1..].copy_from_slice(&x[..n - 1]);
        let cte = conditional_transfer_entropy(&x, &y, &z, 1);
        assert!(
            (cte - std::f64::consts::LN_2).abs() < 0.05,
            "CTE(X->Y|Z) = {cte}, expected ln 2"
        );
        // Conditioning on an unrelated series cannot create information flow
        // in the reverse direction.
        let reverse = conditional_transfer_entropy(&y, &x, &z, 1);
        assert!(reverse < 0.01, "CTE(Y->X|Z) = {reverse}");
    }
    #[test]
    fn test_conditional_transfer_entropy_respects_lag() {
        let n = 20_000;
        let x = coin_flips(n, 0x1357_9BDF_0246_8ACE);
        let z = coin_flips(n, 0x0F1E_2D3C_4B5A_6978);
        let mut y = vec![0.0_f64; n];
        for t in 2..n {
            y[t] = if (x[t - 2] > 0.5) != (z[t - 2] > 0.5) {
                1.0
            } else {
                0.0
            };
        }
        let at_lag_2 = conditional_transfer_entropy(&x, &y, &z, 2);
        let at_lag_1 = conditional_transfer_entropy(&x, &y, &z, 1);
        assert!(
            (at_lag_2 - std::f64::consts::LN_2).abs() < 0.05,
            "lag 2: {at_lag_2}"
        );
        assert!(at_lag_1 < 0.01, "lag 1: {at_lag_1}");
    }
    #[test]
    fn test_conditional_transfer_entropy_degenerate_inputs() {
        let x = coin_flips(100, 7);
        let y = coin_flips(100, 11);
        let z = coin_flips(100, 13);
        // Too short for the requested lag.
        assert_eq!(conditional_transfer_entropy(&x[..2], &y, &z, 1), 0.0);
        assert_eq!(conditional_transfer_entropy(&x, &y, &z[..5], 5), 0.0);
        // Constant data has a single bin and therefore no information flow.
        let flat = vec![1.0_f64; 100];
        assert_eq!(conditional_transfer_entropy(&flat, &flat, &flat, 1), 0.0);
        // Mismatched lengths use the shortest series.
        let truncated = conditional_transfer_entropy(&x[..60], &y[..60], &z[..60], 1);
        let mismatched = conditional_transfer_entropy(&x[..60], &y, &z, 1);
        assert_eq!(truncated, mismatched);
        // A conditional mutual information is never negative.
        assert!(conditional_transfer_entropy(&x, &y, &z, 1) >= 0.0);
    }
    #[test]
    fn test_shannon_hartley_zero_snr() {
        let c = shannon_hartley_capacity(1e6, 0.0);
        assert!(c.abs() < EPS);
    }
    #[test]
    fn test_shannon_hartley_increases_with_snr() {
        let c1 = shannon_hartley_capacity(1.0, 1.0);
        let c2 = shannon_hartley_capacity(1.0, 10.0);
        assert!(c2 > c1);
    }
    #[test]
    fn test_channel_capacity_bsc_zero_noise() {
        let t = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let c = channel_capacity_blahut(&t);
        assert!((c - 1.0).abs() < 1e-4, "BSC(0) capacity = {c}");
    }
    #[test]
    fn test_channel_capacity_bsc_half_noise() {
        let t = vec![vec![0.5, 0.5], vec![0.5, 0.5]];
        let c = channel_capacity_blahut(&t);
        assert!(c < 1e-3, "BSC(0.5) capacity = {c}");
    }
    #[test]
    fn test_huffman_lengths_single() {
        let lengths = huffman_lengths(&[1.0]);
        assert_eq!(lengths, vec![1]);
    }
    #[test]
    fn test_huffman_expected_length_geq_entropy() {
        let p = [0.5, 0.25, 0.125, 0.125];
        let el = huffman_expected_length(&p);
        let h = entropy_bits(&p);
        assert!(el >= h - EPS);
    }
    #[test]
    fn test_huffman_expected_length_near_entropy() {
        let p = [0.5, 0.25, 0.125, 0.125];
        let el = huffman_expected_length(&p);
        let h = entropy_bits(&p);
        assert!(el < h + 1.0 + EPS);
    }
    #[test]
    fn test_huffman_kraft_inequality() {
        let p = [0.4, 0.3, 0.2, 0.1];
        let lengths = huffman_lengths(&p);
        assert!(kraft_inequality(&lengths));
    }
    #[test]
    fn test_binary_rate_distortion_zero_distortion() {
        let r = binary_rate_distortion(0.0);
        assert!((r - 1.0).abs() < EPS);
    }
    #[test]
    fn test_binary_rate_distortion_half_distortion() {
        let r = binary_rate_distortion(0.5);
        assert!(r.abs() < EPS, "R(0.5) = {r}");
    }
    #[test]
    fn test_gaussian_rate_distortion_zero_when_d_geq_variance() {
        let r = gaussian_rate_distortion(1.0, 2.0);
        assert_eq!(r, 0.0);
    }
    #[test]
    fn test_gaussian_rate_distortion_positive_when_d_small() {
        let r = gaussian_rate_distortion(4.0, 1.0);
        assert!(r > 0.0, "R = {r}");
    }
    #[test]
    fn test_aic_increases_with_k() {
        let ll = -50.0;
        let a1 = aic(1, ll);
        let a2 = aic(5, ll);
        assert!(a2 > a1);
    }
    #[test]
    fn test_bic_increases_with_k() {
        let ll = -50.0;
        let b1 = bic(1, 100, ll);
        let b2 = bic(5, 100, ll);
        assert!(b2 > b1);
    }
    #[test]
    fn test_aicc_approaches_aic_for_large_n() {
        let (k, ll) = (2, -30.0);
        let a = aic(k, ll);
        let ac = aicc(k, 1000, ll);
        assert!((a - ac).abs() < 0.1);
    }
    #[test]
    fn test_model_selector_selects_simplest_among_equal_fit() {
        let mut ms = ModelSelector::new(100);
        ms.add(10, -50.0);
        ms.add(2, -50.0);
        let best = ms.best_aic().unwrap();
        assert_eq!(best, 1, "simpler model should win");
    }
    #[test]
    fn test_fisher_metric_categorical_diagonal() {
        let p = [0.5, 0.5];
        let mat = fisher_metric_categorical(&p);
        assert!((mat[0] - 2.0).abs() < EPS);
        assert!((mat[3] - 2.0).abs() < EPS);
        assert!(mat[1].abs() < EPS);
        assert!(mat[2].abs() < EPS);
    }
    #[test]
    fn test_fisher_info_gaussian_increases_with_precision() {
        let i1 = fisher_information_gaussian_mean(1.0);
        let i2 = fisher_information_gaussian_mean(0.5);
        assert!(i2 > i1);
    }
    #[test]
    fn test_fisher_info_bernoulli_half() {
        let fi = fisher_information_bernoulli(0.5);
        assert!((fi - 4.0).abs() < EPS);
    }
    #[test]
    fn test_geodesic_distance_same_distribution_is_zero() {
        let p = vec![0.3, 0.3, 0.4];
        let fg = FisherGeometry::new(p.clone());
        let d = fg.geodesic_distance(&p);
        assert!(d.abs() < EPS);
    }
    #[test]
    fn test_geodesic_distance_orthogonal_distributions() {
        let p = vec![1.0, 0.0];
        let q = vec![0.0, 1.0];
        let fg = FisherGeometry::new(p);
        let d = fg.geodesic_distance(&q);
        assert!((d - PI).abs() < EPS);
    }
    #[test]
    fn test_discrete_distribution_normalises() {
        let dd = DiscreteDistribution::from_weights(&[1.0, 1.0, 2.0]);
        let sum: f64 = dd.probs.iter().sum();
        assert!((sum - 1.0).abs() < EPS);
    }
    #[test]
    fn test_discrete_distribution_entropy_bits() {
        let dd = DiscreteDistribution::from_weights(&[1.0, 1.0]);
        assert!((dd.entropy_bits() - 1.0).abs() < EPS);
    }
    #[test]
    fn test_rate_distortion_distortion_at_zero_rate() {
        let rd = RateDistortion::new(4.0);
        let d = rd.gaussian_distortion_at_rate(0.0);
        assert!((d - 4.0).abs() < EPS);
    }
    #[test]
    fn test_rate_distortion_rd_inverse() {
        let rd = RateDistortion::new(4.0);
        let rate = 1.0;
        let d = rd.gaussian_distortion_at_rate(rate);
        let r_back = rd.gaussian_rd(d);
        assert!((r_back - rate).abs() < 1e-6);
    }
    #[test]
    fn test_code_redundancy_huffman_non_negative() {
        let p = [0.4, 0.3, 0.2, 0.1];
        let lengths = huffman_lengths(&p);
        let r = code_redundancy(&p, &lengths);
        assert!(r >= -EPS);
    }
    #[test]
    fn test_kraft_inequality_huffman() {
        let p = [0.5, 0.25, 0.125, 0.0625, 0.0625];
        let lengths = huffman_lengths(&p);
        assert!(kraft_inequality(&lengths));
    }
    #[test]
    fn test_tsallis_entropy_non_negative() {
        let p = [0.2, 0.3, 0.5];
        assert!(tsallis_entropy(&p, 2.0) >= -EPS);
    }
    #[test]
    fn test_tsallis_q1_approaches_shannon_nats() {
        let p = [0.2, 0.3, 0.5];
        let ts = tsallis_entropy(&p, 1.0);
        let sh = entropy_nats(&p);
        assert!((ts - sh).abs() < 1e-6);
    }
    #[test]
    fn test_hellinger_distance_zero_for_equal() {
        let p = vec![0.3, 0.3, 0.4];
        assert!(hellinger_distance_sq(&p, &p).abs() < EPS);
    }
    #[test]
    fn test_hellinger_distance_bounded() {
        let p = vec![1.0, 0.0];
        let q = vec![0.0, 1.0];
        let h = hellinger_distance_sq(&p, &q);
        assert!((h - 1.0).abs() < EPS);
    }
    #[test]
    fn test_bhattacharyya_distance_zero_for_equal() {
        let p = vec![0.5, 0.5];
        let d = bhattacharyya_distance(&p, &p);
        assert!(d.abs() < EPS);
    }
    #[test]
    fn test_chi_squared_divergence_zero_for_equal() {
        let p = vec![0.3, 0.7];
        assert!(chi_squared_divergence(&p, &p).abs() < EPS);
    }
    #[test]
    fn test_fisher_rao_distance_zero_for_same() {
        let p = vec![0.2, 0.3, 0.5];
        assert!(fisher_rao_distance(&p, &p).abs() < EPS);
    }
    #[test]
    fn test_fisher_rao_gaussian_same_point() {
        let d = fisher_rao_gaussian(0.0, 1.0, 0.0, 1.0);
        assert!(d.abs() < EPS);
    }
    #[test]
    fn test_fisher_rao_gaussian_different_means() {
        let d = fisher_rao_gaussian(0.0, 1.0, 3.0, 1.0);
        assert!(d > 0.0);
    }
    #[test]
    fn test_fisher_rao_gaussian_exact_scale_distance() {
        // N(0,1) -> N(0,e) lies on a vertical geodesic: d = sqrt(2) * ln(e).
        let d = fisher_rao_gaussian(0.0, 1.0, 0.0, std::f64::consts::E);
        assert!(
            (d - std::f64::consts::SQRT_2).abs() < 1e-12,
            "expected sqrt(2), got {d}"
        );
        // The same distance is reached from any common scale and mean.
        let d = fisher_rao_gaussian(3.0, 2.0, 3.0, 2.0 * std::f64::consts::E);
        assert!((d - std::f64::consts::SQRT_2).abs() < 1e-12, "got {d}");
    }
    #[test]
    fn test_fisher_rao_gaussian_exact_mean_distance() {
        // Equal sigma: cosh(d / sqrt(2)) = 1 + dmu^2 / (4 sigma^2).  Choosing
        // dmu = 2 sqrt(2) sinh(1/2) sigma gives d / sqrt(2) = 1 exactly.
        let sigma = 1.7;
        let dmu = 2.0 * std::f64::consts::SQRT_2 * 0.5_f64.sinh() * sigma;
        let d = fisher_rao_gaussian(0.4, sigma, 0.4 + dmu, sigma);
        assert!(
            (d - std::f64::consts::SQRT_2).abs() < 1e-12,
            "expected sqrt(2), got {d}"
        );
    }
    #[test]
    fn test_fisher_rao_gaussian_matches_metric_for_nearby_points() {
        // For nearby distributions d^2 = (dmu^2 + 2 dsigma^2) / sigma^2.  The
        // separation is far below sqrt(machine epsilon), where acosh(1 + x)
        // evaluated naively would return exactly zero.
        let eps = 1e-9;
        let d = fisher_rao_gaussian(0.0, 1.0, eps, 1.0);
        assert!((d - eps).abs() < 1e-17, "d = {d}");
        // `1.0 + eps` is not exactly representable, so the separation that was
        // actually passed in is the floating-point difference `sigma2 - 1`.
        let sigma2 = 1.0 + eps;
        let ds = sigma2 - 1.0;
        let d = fisher_rao_gaussian(0.0, 1.0, 0.0, sigma2);
        let expected = std::f64::consts::SQRT_2 * ds / sigma2.sqrt();
        assert!((d - expected).abs() < 1e-17, "d = {d}, expected {expected}");
    }
    #[test]
    fn test_fisher_rao_gaussian_scale_invariance_and_degenerate_inputs() {
        let a = fisher_rao_gaussian(0.3, 1.2, -1.1, 0.7);
        let k = 13.0;
        let b = fisher_rao_gaussian(k * 0.3, k * 1.2, k * -1.1, k * 0.7);
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        assert!(fisher_rao_gaussian(0.0, 0.0, 0.0, 1.0).is_infinite());
        assert!(fisher_rao_gaussian(0.0, 1.0, 0.0, 0.0).is_infinite());
        // Huge separations stay finite and monotone.
        let far = fisher_rao_gaussian(0.0, 1.0, 1.0e150, 1.0);
        assert!(far.is_finite() && far > 100.0, "far = {far}");
    }
    #[test]
    fn test_fisher_rao_gaussian_metric_axioms_on_random_triples() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..500 {
            let pts: Vec<(f64, f64)> = (0..3)
                .map(|_| (8.0 * next() - 4.0, 0.05 + 4.0 * next()))
                .collect();
            let d =
                |a: usize, b: usize| fisher_rao_gaussian(pts[a].0, pts[a].1, pts[b].0, pts[b].1);
            assert!(d(0, 0).abs() < 1e-12);
            assert!((d(0, 1) - d(1, 0)).abs() < 1e-12, "symmetry");
            assert!(
                d(0, 2) <= d(0, 1) + d(1, 2) + 1e-12,
                "triangle inequality: {} > {} + {}",
                d(0, 2),
                d(0, 1),
                d(1, 2)
            );
        }
    }
    #[test]
    fn test_gaussian_info_geometry_distance_uses_exact_fisher_rao() {
        let g1 = GaussianInfoGeometry::new(0.0, 1.0);
        let g2 = GaussianInfoGeometry::new(0.0, std::f64::consts::E);
        assert!((g1.distance_to(&g2) - std::f64::consts::SQRT_2).abs() < 1e-12);
    }
    #[test]
    fn test_collision_entropy_uniform() {
        let p = [0.25, 0.25, 0.25, 0.25];
        let h2 = collision_entropy(&p);
        assert!((h2 - 2.0).abs() < EPS);
    }
    #[test]
    fn test_hartley_entropy_four_symbols() {
        let p = [0.25, 0.25, 0.25, 0.25];
        let h0 = hartley_entropy(&p);
        assert!((h0 - 2.0).abs() < EPS);
    }
    #[test]
    fn test_bec_capacity_zero_erasure() {
        assert!((bec_capacity(0.0) - 1.0).abs() < EPS);
    }
    #[test]
    fn test_bsc_capacity_zero_error() {
        assert!((bsc_capacity(0.0) - 1.0).abs() < EPS);
    }
    #[test]
    fn test_max_entropy_uniform_sums_to_one() {
        let p = max_entropy_uniform(5);
        let s: f64 = p.iter().sum();
        assert!((s - 1.0).abs() < EPS);
    }
    #[test]
    fn test_max_entropy_mean_constraint_achieves_mean() {
        let p = max_entropy_mean_constraint(10, 3.0);
        let mean: f64 = p.iter().enumerate().map(|(i, &pi)| i as f64 * pi).sum();
        assert!((mean - 3.0).abs() < 0.1, "mean = {mean}");
    }
    #[test]
    fn test_fisher_info_poisson() {
        let fi = fisher_information_poisson(2.0);
        assert!((fi - 0.5).abs() < EPS);
    }
    #[test]
    fn test_fisher_info_exponential() {
        let fi = fisher_information_exponential(3.0);
        assert!((fi - 1.0 / 9.0).abs() < EPS);
    }
    #[test]
    fn test_gaussian_info_geometry_distance_symmetric() {
        let g1 = GaussianInfoGeometry::new(0.0, 1.0);
        let g2 = GaussianInfoGeometry::new(1.0, 2.0);
        let d1 = g1.distance_to(&g2);
        let d2 = g2.distance_to(&g1);
        assert!((d1 - d2).abs() < EPS);
    }
    #[test]
    fn test_exponential_map_at_center_stays_normalized() {
        let p = vec![0.25, 0.25, 0.25, 0.25];
        let v = vec![0.01, -0.01, 0.01, -0.01];
        let q = exponential_map_simplex(&p, &v, 0.5);
        let s: f64 = q.iter().sum();
        assert!((s - 1.0).abs() < 0.01, "sum = {s}");
    }
    #[test]
    fn test_scalar_curvature_two_simplex() {
        let fg = FisherGeometry::new(vec![0.3, 0.3, 0.4]);
        let r = fg.scalar_curvature();
        assert!((r - 0.5).abs() < EPS);
    }
    #[test]
    fn test_variation_of_information_non_negative() {
        let joint = vec![vec![0.1, 0.4], vec![0.3, 0.2]];
        assert!(variation_of_information(&joint) >= -EPS);
    }
    #[test]
    fn test_laplace_differential_entropy() {
        let h = laplace_differential_entropy(1.0);
        let expected = 1.0 + 2.0_f64.ln();
        assert!((h - expected).abs() < EPS);
    }
    #[test]
    fn test_rd_curve_length() {
        let rd = RateDistortion::new(4.0);
        let curve = rd.rd_curve(10);
        assert_eq!(curve.len(), 10);
    }
    #[test]
    fn test_kl_divergence_bits_matches_nats() {
        let p = vec![0.3, 0.7];
        let q = vec![0.5, 0.5];
        let kl_nats = kl_divergence(&p, &q);
        let kl_bits = kl_divergence_bits(&p, &q);
        assert!((kl_bits - kl_nats / 2.0_f64.ln()).abs() < 1e-10);
    }
    #[test]
    fn test_mdl_increases_with_k() {
        let m1 = mdl(1, 100, -50.0);
        let m2 = mdl(5, 100, -50.0);
        assert!(m2 > m1);
    }
    #[test]
    fn test_uniform_rate_distortion_zero_at_max_d() {
        let r = uniform_rate_distortion(4, 0.75);
        assert!(r.abs() < 0.01, "R = {r}");
    }
    #[test]
    fn test_gaussian_distortion_at_rate_zero() {
        let d = gaussian_distortion_at_rate(4.0, 0.0);
        assert!((d - 4.0).abs() < EPS);
    }
    #[test]
    fn test_shannon_fano_lengths_non_zero() {
        let p = [0.5, 0.25, 0.125, 0.125];
        let lengths = shannon_fano_lengths(&p);
        for &l in &lengths {
            assert!(l > 0);
        }
    }
    #[test]
    fn test_fisher_matrix_gaussian_positive_diagonal() {
        let m = fisher_matrix_gaussian(2.0);
        assert!(m[0] > 0.0);
        assert!(m[3] > 0.0);
        assert!(m[1].abs() < EPS);
    }
}
