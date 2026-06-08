//! Auto-generated test module (consolidated from inline `#[cfg(test)] mod` blocks)

use std::f64::consts::PI;

use super::*;

#[cfg(test)]
mod tests_2 {
    use super::*;
    #[test]
    fn test_bisect_sqrt2() {
        let result = bisect(|x| x * x - 2.0, 1.0, 2.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - 2.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_newton_sqrt2() {
        let result = newton(|x| x * x - 2.0, |x| 2.0 * x, 1.5, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 2.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_brent_cubic() {
        let result = brent(|x| x * x * x - x - 2.0, 1.0, 2.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - 1.5213797).abs() < 1e-5);
    }
    #[test]
    fn test_simpson_sin() {
        let result = simpson(f64::sin, 0.0, PI, 100);
        assert!((result - 2.0).abs() < 1e-6, "Simpson: {result}");
    }
    #[test]
    fn test_gauss_legendre5_polynomial() {
        let result = gauss_legendre5(|x| x * x * x * x, 0.0, 1.0);
        assert!((result - 0.2).abs() < 1e-12, "GL5: {result}");
    }
    #[test]
    fn test_adaptive_integrate() {
        let result = adaptive_integrate(|x| x * x, 0.0, 1.0, 1e-10, 20);
        assert!((result - 1.0 / 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_finite_diff_central_cos() {
        let d = finite_diff_central(f64::sin, 1.0, 1e-6);
        assert!((d - 1.0_f64.cos()).abs() < 1e-8);
    }
    #[test]
    fn test_erf_known_values() {
        assert!(erf(0.0).abs() < 1e-7, "erf(0)={}", erf(0.0));
        assert!((erf(1.0) - 0.842_700_792_9).abs() < 1e-5);
        assert!((erf(2.0) - 0.995_322_265_0).abs() < 1e-5);
    }
    #[test]
    fn test_normal_cdf_symmetry() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 1e-10);
        assert!((normal_cdf(-1.0) + normal_cdf(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_bessel_j0_zeros() {
        let j0 = bessel_j0(2.4048);
        assert!(j0.abs() < 0.01, "J0 near first zero: {j0}");
    }
    #[test]
    fn test_gamma_integers() {
        assert!((gamma(1.0) - 1.0).abs() < 1e-10);
        assert!((gamma(2.0) - 1.0).abs() < 1e-10);
        assert!((gamma(3.0) - 2.0).abs() < 1e-10);
        assert!((gamma(4.0) - 6.0).abs() < 1e-10);
        assert!((gamma(5.0) - 24.0).abs() < 1e-8);
    }
    #[test]
    fn test_horner() {
        assert!((horner(&[2.0, 3.0, 4.0], 2.0) - 24.0).abs() < 1e-10);
    }
    #[test]
    fn test_legendre_orthogonality() {
        assert!((legendre(0, 1.0) - 1.0).abs() < 1e-10);
        assert!((legendre(1, 1.0) - 1.0).abs() < 1e-10);
        assert!((legendre(2, 0.0) + 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_smoothstep() {
        assert!(smoothstep(-1.0).abs() < 1e-10);
        assert!((smoothstep(1.0) - 1.0).abs() < 1e-10);
        assert!((smoothstep(0.5) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_online_stats_welford() {
        let mut stats = OnlineStats::new();
        for x in [1.0, 2.0, 3.0, 4.0, 5.0] {
            stats.update(x);
        }
        assert!((stats.mean() - 3.0).abs() < 1e-10);
        assert!((stats.variance() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_cfl_number() {
        let cfl = cfl_number(10.0, 0.01, 0.1);
        assert!((cfl - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_dt_from_cfl() {
        let dt = dt_from_cfl(0.5, 100.0, 0.1);
        assert!((dt - 0.0005).abs() < 1e-10);
    }
    #[test]
    fn test_trapezoid_tabulated() {
        let x = [0.0, 1.0, 2.0, 3.0];
        let y = [0.0, 1.0, 4.0, 9.0];
        let result = trapezoid_tabulated(&x, &y);
        assert!(result > 0.0, "trapezoidal result should be positive");
    }
    #[test]
    fn test_bilinear() {
        let result = bilinear_arr([0.0, 1.0, 0.0, 1.0], 0.5, 0.5);
        assert!((result - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_romberg_sin() {
        let result = romberg(f64::sin, 0.0, PI, 10);
        assert!((result - 2.0).abs() < 1e-9, "Romberg: {result}");
    }
    #[test]
    fn test_romberg_polynomial() {
        let result = romberg(|x| x * x * x, 0.0, 1.0, 8);
        assert!((result - 0.25).abs() < 1e-10, "Romberg x^3: {result}");
    }
    #[test]
    fn test_gauss_legendre_n_sin() {
        let result = gauss_legendre_n(f64::sin, 0.0, PI, 6);
        assert!((result - 2.0).abs() < 1e-6, "GL6 sin: {result}");
    }
    #[test]
    fn test_gauss_legendre_n_degree9() {
        let result = gauss_legendre_n(|x| x.powi(8), -1.0, 1.0, 6);
        assert!((result - 2.0 / 9.0).abs() < 1e-10, "GL6 x^8: {result}");
    }
    #[test]
    fn test_gauss_legendre_nw_orders() {
        for n in [2, 3, 4, 6] {
            let result = gauss_legendre_n(|x| x * x, -1.0, 1.0, n);
            assert!((result - 2.0 / 3.0).abs() < 1e-8, "GL{n} x^2: {result}");
        }
    }
    #[test]
    fn test_richardson_derivative_cos() {
        let d = richardson_derivative(f64::sin, 1.0, 1e-4);
        assert!(
            (d - 1.0_f64.cos()).abs() < 1e-10,
            "Richardson d/dx sin: {d}"
        );
    }
    #[test]
    fn test_richardson_derivative_polynomial() {
        let d = richardson_derivative(|x: f64| x.powi(3), 2.0, 1e-4);
        assert!((d - 12.0).abs() < 1e-8, "Richardson d/dx x^3: {d}");
    }
    #[test]
    fn test_richardson_second_derivative() {
        let d2 = richardson_second_derivative(f64::sin, 1.0, 1e-4);
        assert!(
            (d2 + 1.0_f64.sin()).abs() < 1e-5,
            "Richardson d2/dx2 sin: {d2}"
        );
    }
    #[test]
    fn test_richardson_second_derivative_polynomial() {
        let d2 = richardson_second_derivative(|x: f64| x.powi(4), 2.0, 1e-3);
        assert!((d2 - 48.0).abs() < 1e-5, "Richardson d2/dx2 x^4: {d2}");
    }
    #[test]
    fn test_find_x_for_value() {
        let result = find_x_for_value(|x| x * x, 9.0, 0.0, 5.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - 3.0).abs() < 1e-9, "find_x: {}", result.root);
    }
    #[test]
    fn test_beta_symmetry() {
        let b1 = beta(2.0, 3.0);
        let b2 = beta(3.0, 2.0);
        assert!((b1 - b2).abs() < 1e-12, "Beta symmetry: {b1} vs {b2}");
    }
    #[test]
    fn test_beta_known_value() {
        assert!((beta(1.0, 1.0) - 1.0).abs() < 1e-10);
        assert!((beta(2.0, 2.0) - 1.0 / 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_beta_half_integer() {
        assert!(
            (beta(0.5, 0.5) - PI).abs() < 1e-8,
            "B(1/2,1/2)={}",
            beta(0.5, 0.5)
        );
    }
    #[test]
    fn test_digamma_integer() {
        let psi1 = digamma(1.0);
        assert!((psi1 + 0.5772156649).abs() < 1e-5, "ψ(1)={psi1}");
        let psi2 = digamma(2.0);
        assert!((psi2 - 0.4227843351).abs() < 1e-5, "ψ(2)={psi2}");
    }
    #[test]
    fn test_digamma_recurrence() {
        let x = 3.7;
        assert!((digamma(x + 1.0) - digamma(x) - 1.0 / x).abs() < 1e-8);
    }
    #[test]
    fn test_erfinv_roundtrip() {
        for &x in &[0.0_f64, 0.3, -0.3, 0.5] {
            let y = erf(x);
            let x_back = erfinv(y);
            assert!((x_back - x).abs() < 5e-3, "erfinv(erf({x}))={x_back}");
        }
    }
    #[test]
    fn test_erfinv_known() {
        assert!(erfinv(0.0).abs() < 1e-6, "erfinv(0)={}", erfinv(0.0));
    }
    #[test]
    fn test_ei_small() {
        let val = ei(1.0);
        assert!((val - 1.8951178163559366).abs() < 1e-6, "Ei(1)={val}");
    }
    #[test]
    fn test_ei_large() {
        let val = ei(10.0);
        assert!((val - 2492.2289).abs() < 0.5, "Ei(10)={val}");
    }
    #[test]
    fn test_ei_negative_returns_neg_inf() {
        assert_eq!(ei(-1.0), f64::NEG_INFINITY);
        assert_eq!(ei(0.0), f64::NEG_INFINITY);
    }
    #[test]
    fn test_secant_sqrt2() {
        let result = secant(|x| x * x - 2.0, 1.0, 2.0, 1e-12, 50);
        assert!(result.converged, "secant did not converge");
        assert!(
            (result.root - 2.0_f64.sqrt()).abs() < 1e-10,
            "root={}",
            result.root
        );
    }
    #[test]
    fn test_secant_cubic() {
        let result = secant(|x| x * x * x - x - 2.0, 1.0, 2.0, 1e-10, 50);
        assert!(result.converged);
        assert!(
            (result.root - 1.5213797).abs() < 1e-6,
            "root={}",
            result.root
        );
    }
    #[test]
    fn test_secant_sin() {
        let result = secant(f64::sin, 3.0, 4.0, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - PI).abs() < 1e-10, "root={}", result.root);
    }
    #[test]
    fn test_secant_exact_root() {
        let result = secant(|x| x - 5.0, 4.0, 6.0, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_continued_fraction_constant() {
        let val = continued_fraction(1.0, &[1.0], &[1.0]);
        assert!((val - 2.0).abs() < 1e-12, "cf={val}");
    }
    #[test]
    fn test_continued_fraction_golden_ratio() {
        let n = 30;
        let a = vec![1.0_f64; n];
        let b = vec![1.0_f64; n];
        let phi = continued_fraction(1.0, &a, &b);
        let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((phi - golden).abs() < 1e-6, "phi={phi}, golden={golden}");
    }
    #[test]
    fn test_continued_fraction_empty() {
        let val = continued_fraction(3.125, &[], &[]);
        assert!((val - 3.125).abs() < 1e-12);
    }
    #[test]
    fn test_sqrt2_cf_converges() {
        let approx = sqrt2_cf(40);
        assert!((approx - 2.0_f64.sqrt()).abs() < 1e-12, "sqrt2_cf={approx}");
    }
    #[test]
    fn test_sqrt2_cf_improves_with_depth() {
        let d5 = sqrt2_cf(5);
        let d20 = sqrt2_cf(20);
        let exact = 2.0_f64.sqrt();
        assert!(
            (d20 - exact).abs() < (d5 - exact).abs(),
            "deeper CF should be more accurate"
        );
    }
    #[test]
    fn test_bernoulli_b0_is_one() {
        let b = bernoulli_numbers(1);
        assert!((b[0] - 1.0).abs() < 1e-12, "B0={}", b[0]);
    }
    #[test]
    fn test_bernoulli_b1_is_minus_half() {
        let b = bernoulli_numbers(2);
        assert!((b[1] + 0.5).abs() < 1e-12, "B1={}", b[1]);
    }
    #[test]
    fn test_bernoulli_odd_are_zero() {
        let b = bernoulli_numbers(10);
        for k in [3, 5, 7, 9] {
            assert!(b[k].abs() < 1e-12, "B{k}={}", b[k]);
        }
    }
    #[test]
    fn test_bernoulli_b2_is_one_sixth() {
        let b = bernoulli_numbers(3);
        assert!((b[2] - 1.0 / 6.0).abs() < 1e-12, "B2={}", b[2]);
    }
    #[test]
    fn test_bernoulli_b4_is_minus_one_thirtieth() {
        let b = bernoulli_numbers(5);
        assert!((b[4] + 1.0 / 30.0).abs() < 1e-12, "B4={}", b[4]);
    }
    #[test]
    fn test_bernoulli_b6_is_one_fortysecond() {
        let b = bernoulli_numbers(7);
        assert!((b[6] - 1.0 / 42.0).abs() < 1e-12, "B6={}", b[6]);
    }
    #[test]
    fn test_bernoulli_count() {
        let b = bernoulli_numbers(8);
        assert_eq!(b.len(), 8);
    }
    #[test]
    fn test_bernoulli_empty() {
        let b = bernoulli_numbers(0);
        assert!(b.is_empty());
    }
    #[test]
    fn test_stirling_large_n_accuracy() {
        let err = stirling_relative_error(100.0);
        assert!(err < 1e-10, "Stirling relative error at n=100: {err}");
    }
    #[test]
    fn test_stirling_n10() {
        let exact_ln_10_fact = 15.104_412_573_075_518;
        let approx = stirling_ln_factorial(10.0);
        assert!(
            (approx - exact_ln_10_fact).abs() < 0.01,
            "Stirling ln(10!)={approx}"
        );
    }
    #[test]
    fn test_stirling_n1000_very_accurate() {
        let err = stirling_relative_error(1000.0);
        assert!(err < 1e-8, "Stirling relative error at n=1000: {err}");
    }
    #[test]
    fn test_stirling_factorial_50() {
        let exact_ln = lgamma(51.0);
        let approx_ln = stirling_ln_factorial(50.0);
        let rel_err = (approx_ln - exact_ln).abs() / exact_ln.abs();
        assert!(
            rel_err < 1e-10,
            "Stirling ln(50!)={approx_ln}, exact={exact_ln}"
        );
    }
    #[test]
    fn test_stirling_zero() {
        assert!(stirling_ln_factorial(0.0).abs() < 1e-12);
    }
    #[test]
    fn test_bessel_j1_at_zero() {
        assert!(bessel_j1(0.0).abs() < 1e-10, "J1(0)={}", bessel_j1(0.0));
    }
    #[test]
    fn test_bessel_j1_negative_antisymmetry() {
        let x = 2.5;
        assert!((bessel_j1(-x) + bessel_j1(x)).abs() < 1e-12);
    }
    #[test]
    fn test_bessel_j0_at_large_x() {
        assert!(bessel_j0(100.0).abs() < 1.0);
    }
    #[test]
    fn test_erfc_complement() {
        for &x in &[0.0_f64, 0.5, 1.0, 2.0] {
            let sum = erf(x) + erfc(x);
            assert!((sum - 1.0).abs() < 1e-10, "erf+erfc≠1 at x={x}: {sum}");
        }
    }
    #[test]
    fn test_gamma_half() {
        let g = gamma(0.5);
        assert!((g - PI.sqrt()).abs() < 1e-8, "Γ(0.5)={g}");
    }
    #[test]
    fn test_brent_sin_root() {
        let result = brent(f64::sin, 3.0, 4.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - PI).abs() < 1e-10);
    }
    #[test]
    fn test_chebyshev_t_known() {
        assert!((chebyshev_t(0, 0.7) - 1.0).abs() < 1e-12);
        assert!((chebyshev_t(1, 0.7) - 0.7).abs() < 1e-12);
        assert!((chebyshev_t(2, 0.7) - (2.0 * 0.7 * 0.7 - 1.0)).abs() < 1e-12);
    }
    #[test]
    fn test_finite_diff_second_sin() {
        let d2 = finite_diff_second(f64::sin, 1.0, 1e-4);
        assert!((d2 + 1.0_f64.sin()).abs() < 1e-5, "d2={d2}");
    }
    #[test]
    fn test_hessian_quadratic() {
        let f = |v: &[f64]| v[0] * v[0] + 2.0 * v[1] * v[1];
        let h = hessian(&f, &[1.0, 1.0], 1e-4);
        assert!((h[0][0] - 2.0).abs() < 1e-6, "H[0][0]={}", h[0][0]);
        assert!((h[1][1] - 4.0).abs() < 1e-6, "H[1][1]={}", h[1][1]);
        assert!(h[0][1].abs() < 1e-6, "H[0][1]={}", h[0][1]);
    }
    #[test]
    fn test_regula_falsi_sqrt3() {
        let result = regula_falsi(|x| x * x - 3.0, 1.0, 2.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - 3.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_online_stats_merge() {
        let mut s1 = OnlineStats::new();
        let mut s2 = OnlineStats::new();
        for &x in &[1.0_f64, 2.0, 3.0] {
            s1.update(x);
        }
        for &x in &[4.0_f64, 5.0] {
            s2.update(x);
        }
        s1.merge(&s2);
        assert_eq!(s1.count(), 5);
        assert!((s1.mean() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_gradient_3d_sphere() {
        let f = |v: [f64; 3]| v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let g = gradient_3d(f, [1.0, 2.0, 3.0], 1e-5);
        assert!((g[0] - 2.0).abs() < 1e-8);
        assert!((g[1] - 4.0).abs() < 1e-8);
        assert!((g[2] - 6.0).abs() < 1e-8);
    }
    #[test]
    fn test_inc_beta_symmetry() {
        let x = 0.3;
        let a = 2.0;
        let b = 3.0;
        assert!((inc_beta(x, a, b) + inc_beta(1.0 - x, b, a) - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_trilinear_arr_corners() {
        let vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        assert!((trilinear_arr(vals, 0.0, 0.0, 0.0) - 1.0).abs() < 1e-12);
        assert!((trilinear_arr(vals, 1.0, 1.0, 1.0) - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_smootherstep_boundary() {
        assert!((smootherstep(0.0)).abs() < 1e-12);
        assert!((smootherstep(1.0) - 1.0).abs() < 1e-12);
        let d = (smootherstep(0.001) - smootherstep(0.0)) / 0.001;
        assert!(d < 0.01, "derivative at 0 should be near 0: {d}");
    }
    #[test]
    fn boys_fn_n0_at_zero() {
        let f = boys_fn(0.0, 0);
        assert!((f[0] - 1.0).abs() < 1e-12, "F_0(0) = 1");
    }
    #[test]
    fn boys_fn_n1_at_zero() {
        let f = boys_fn(0.0, 1);
        assert!((f[1] - 1.0 / 3.0).abs() < 1e-12, "F_1(0) = 1/3");
    }
    #[test]
    fn boys_fn_n0_at_one() {
        let f = boys_fn(1.0, 0);
        assert!(
            (f[0] - 0.746824132812427).abs() < 1e-7,
            "F_0(1) reference value, got {}",
            f[0]
        );
    }
    #[test]
    fn boys_fn_large_x_asymptotic() {
        let x = 50.0_f64;
        let f = boys_fn(x, 2);
        let f0_expected = (std::f64::consts::PI / (4.0 * x)).sqrt();
        assert!(
            (f[0] - f0_expected).abs() < 1e-6,
            "F_0(50) asymptotic, got {}",
            f[0]
        );
        assert!(f[1] > 0.0 && f[2] > 0.0, "all boys values positive");
    }
    #[test]
    fn boys_fn_recurrence_consistency() {
        let x = 2.5;
        let f = boys_fn(x, 5);
        for n in 0..5 {
            let lhs = (2 * n + 1) as f64 * f[n] - (-x).exp();
            let rhs = 2.0 * x * f[n + 1];
            assert!((lhs - rhs).abs() < 1e-8, "recurrence at n={n}");
        }
    }
    #[test]
    fn incomplete_gamma_lower_closed_form() {
        for &x in &[0.5_f64, 1.0, 2.0, 5.0] {
            let p = incomplete_gamma_lower(1.0, x);
            let expected = 1.0 - (-x).exp();
            assert!(
                (p - expected).abs() < 1e-10,
                "P(1,{x}) = 1-exp(-{x}), got {p}"
            );
        }
    }
    #[test]
    fn incomplete_gamma_lower_half_integer() {
        let p = incomplete_gamma_lower(0.5, 1.0);
        let expected = 0.8427007929;
        assert!((p - expected).abs() < 1e-6, "P(0.5,1) ≈ erf(1), got {p}");
    }
    #[test]
    fn incomplete_gamma_lower_boundary() {
        assert_eq!(incomplete_gamma_lower(1.0, 0.0), 0.0);
        let p_large = incomplete_gamma_lower(2.0, 100.0);
        assert!((p_large - 1.0).abs() < 1e-6, "P(a,inf) → 1, got {p_large}");
    }
}
