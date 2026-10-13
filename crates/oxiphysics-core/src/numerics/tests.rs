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
    /// `F_n(x)` for `n = 0..=4`, evaluated with 50-digit arithmetic (mpmath).
    const BOYS_REFERENCE: [(f64, [f64; 5]); 11] = [
        (
            0.1,
            [
                0.9676433126355918,
                0.3140294729981613,
                0.18625500479262147,
                0.13218802963573895,
                0.10239394707106532,
            ],
        ),
        (
            1.0,
            [
                0.746824132812427,
                0.18947234582049235,
                0.10026879814501737,
                0.06673227477682225,
                0.04962324113315674,
            ],
        ),
        (
            5.0,
            [
                0.39571230961051357,
                0.03889743626114281,
                0.010995436178434296,
                0.004823923389308601,
                0.002702951672607474,
            ],
        ),
        (
            10.0,
            [
                0.2802473905066427,
                0.014010099528844014,
                0.0020992449328384776,
                0.0005225412367214952,
                0.00018061943636439907,
            ],
        ),
        (
            15.0,
            [
                0.22882279832973734,
                0.007627416414247228,
                0.0007627314446807061,
                0.00012711171070276762,
                2.9649202419962384e-5,
            ],
        ),
        (
            20.0,
            [
                0.19816636482997366,
                0.0049541590692205005,
                0.000371561878662697,
                4.644518330399656e-5,
                8.127855549358838e-6,
            ],
        ),
        (
            24.9,
            [
                0.17760094289010062,
                0.003566283993468917,
                0.00021483638484052608,
                2.1569917848475476e-5,
                3.031915855236479e-6,
            ],
        ),
        (
            25.0,
            [
                0.1772453850902791,
                0.003544907701527823,
                0.0002126944618139105,
                2.1269445903632173e-5,
                2.977722148749627e-6,
            ],
        ),
        (
            30.0,
            [
                0.16180215937964007,
                0.002696702656325775,
                0.00013483513281472914,
                1.1236261066334492e-5,
                1.3108971228460868e-6,
            ],
        ),
        (
            50.0,
            [
                0.12533141373155002,
                0.0012533141373155003,
                3.759942411946501e-5,
                1.8799712059732504e-6,
                1.3159798441812752e-7,
            ],
        ),
        (
            100.0,
            [
                0.0886226925452758,
                0.000443113462726379,
                6.646701940895685e-6,
                1.6616754852239212e-7,
                5.815864198283725e-9,
            ],
        ),
    ];
    fn assert_rel_close(actual: f64, expected: f64, rtol: f64, what: &str) {
        let scale = expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            (actual - expected).abs() <= rtol * scale,
            "{what}: got {actual:e}, expected {expected:e}"
        );
    }
    #[test]
    fn boys_fn_matches_high_precision_reference_for_all_branches() {
        for (x, expected) in BOYS_REFERENCE.iter() {
            let f = boys_fn(*x, 4);
            assert_eq!(f.len(), 5);
            for (n, want) in expected.iter().enumerate() {
                assert_rel_close(f[n], *want, 1e-12, &format!("F_{n}({x})"));
            }
        }
    }
    #[test]
    fn boys_fn_f0_alone_matches_erf_closed_form() {
        // F_0(x) = 0.5 * sqrt(pi / x) * erf(sqrt(x)); the order requested must
        // not change the lowest order (series and recursion agree).
        for (x, expected) in BOYS_REFERENCE.iter() {
            let f0_only = boys_fn(*x, 0);
            assert_eq!(f0_only.len(), 1);
            assert_rel_close(f0_only[0], expected[0], 1e-12, &format!("F_0({x})"));
        }
    }
    #[test]
    fn boys_fn_no_longer_diverges_for_moderate_arguments() {
        // The alternating Maclaurin series with 30 terms returned 4.04 at x = 15,
        // 25704 at x = 20 and 2.1e7 at x = 24.9.
        let f15 = boys_fn(15.0, 0)[0];
        let f20 = boys_fn(20.0, 0)[0];
        let f249 = boys_fn(24.9, 0)[0];
        assert!((f15 - 0.2288227983297373).abs() < 1e-12, "F_0(15) = {f15}");
        assert!((f20 - 0.19816636482997366).abs() < 1e-12, "F_0(20) = {f20}");
        assert!(
            (f249 - 0.17760094289010062).abs() < 1e-12,
            "F_0(24.9) = {f249}"
        );
    }
    #[test]
    fn boys_fn_higher_orders_match_simpson_quadrature() {
        for &x in &[0.5_f64, 3.0, 12.0, 24.9, 40.0, 80.0] {
            let f = boys_fn(x, 4);
            for (n, value) in f.iter().enumerate().skip(1) {
                let quad = simpson(
                    |t| t.powi(2 * n as i32) * (-x * t * t).exp(),
                    0.0,
                    1.0,
                    20_000,
                );
                assert_rel_close(*value, quad, 1e-9, &format!("F_{n}({x}) vs Simpson"));
            }
        }
    }
    #[test]
    fn boys_fn_large_order_uses_stable_downward_path() {
        // max_n > x makes the upward recursion unstable; the result must still be
        // accurate (reference values from 50-digit arithmetic).
        let f = boys_fn(30.0, 60);
        assert_rel_close(f[60], 1.4899170350220724e-15, 1e-12, "F_60(30)");
        assert_rel_close(f[30], 1.0430776486261033e-14, 1e-12, "F_30(30)");
        assert_rel_close(f[0], 0.16180215937964007, 1e-12, "F_0(30)");
        let f = boys_fn(50.0, 100);
        assert_rel_close(f[100], 1.8748353891302932e-24, 1e-12, "F_100(50)");
        let f = boys_fn(40.0, 45);
        assert_rel_close(f[45], 2.4032762027370097e-19, 1e-12, "F_45(40)");
    }
    #[test]
    fn boys_fn_upward_recursion_is_accurate_up_to_n_equals_x() {
        let f = boys_fn(25.0, 25);
        assert_rel_close(f[25], 1.6912722018908912e-12, 1e-12, "F_25(25)");
        let f = boys_fn(30.0, 30);
        assert_rel_close(f[30], 1.0430776486261033e-14, 1e-12, "F_30(30)");
    }
    #[test]
    fn boys_fn_is_continuous_across_the_asymptotic_boundary() {
        let delta = 1.0e-9;
        let below = boys_fn(25.0 - delta, 6);
        let above = boys_fn(25.0 + delta, 6);
        for n in 0..=5_usize {
            // dF_n/dx = -F_{n+1}, so the jump over 2*delta is 2*delta*F_{n+1}.
            let expected_jump = 2.0 * delta * boys_fn(25.0, 6)[n + 1];
            let jump = below[n] - above[n];
            assert!(
                (jump - expected_jump).abs() < 1e-13 * below[n],
                "discontinuity at n={n}: jump={jump:e}, expected={expected_jump:e}"
            );
        }
    }
    #[test]
    fn boys_fn_is_positive_and_decreasing_in_order() {
        for &x in &[0.0_f64, 0.01, 0.5, 2.0, 7.0, 13.0, 24.0, 26.0, 40.0, 300.0] {
            let f = boys_fn(x, 10);
            for (n, pair) in f.windows(2).enumerate() {
                assert!(pair[0] > 0.0 && pair[1] > 0.0, "F_n({x}) must be positive");
                assert!(
                    pair[1] < pair[0],
                    "F_{}({x}) = {:e} must be below F_{n}({x}) = {:e}",
                    n + 1,
                    pair[1],
                    pair[0]
                );
            }
        }
    }
    #[test]
    fn boys_fn_recurrence_holds_across_argument_range() {
        for &x in &[0.3_f64, 4.0, 12.0, 24.9, 25.1, 40.0] {
            let f = boys_fn(x, 8);
            for n in 0..8_usize {
                let scale = (2 * n + 1) as f64 * f[n];
                let lhs = scale - 2.0 * x * f[n + 1];
                assert!(
                    (lhs - (-x).exp()).abs() <= 1e-14 * scale,
                    "recurrence n={n} x={x}: lhs={lhs:e}, exp(-x)={:e}",
                    (-x).exp()
                );
            }
        }
    }
    #[test]
    fn boys_fn_negative_argument_matches_reference() {
        let reference: [(f64, [f64; 4]); 3] = [
            (
                -0.5,
                [
                    1.1949576619102276,
                    0.4537636087899005,
                    0.2874304443304266,
                    0.2115690490479952,
                ],
            ),
            (
                -1.0,
                [
                    1.4626517459071815,
                    0.6278150412759318,
                    0.4174183523156249,
                    0.31559503344046036,
                ],
            ),
            (
                -5.0,
                [
                    17.17215777384149,
                    13.124100132873512,
                    10.904085870395607,
                    9.389272975059857,
                ],
            ),
        ];
        for (x, expected) in reference.iter() {
            let f = boys_fn(*x, 3);
            for (n, want) in expected.iter().enumerate() {
                assert_rel_close(f[n], *want, 1e-12, &format!("F_{n}({x})"));
            }
        }
    }
    #[test]
    fn boys_fn_special_arguments() {
        assert!(boys_fn(f64::NAN, 3).iter().all(|v| v.is_nan()));
        assert_eq!(boys_fn(f64::INFINITY, 2), vec![0.0, 0.0, 0.0]);
        let at_zero = boys_fn(0.0, 4);
        for (n, v) in at_zero.iter().enumerate() {
            assert_rel_close(*v, 1.0 / (2 * n + 1) as f64, 1e-15, "F_n(0)");
        }
    }
    /// `(x, a, b, I_x(a, b))` evaluated with 50-digit arithmetic (mpmath).
    const INC_BETA_REFERENCE: [(f64, f64, f64, f64); 14] = [
        (0.3, 1.0, 1.0, 0.3),
        (0.2, 2.0, 3.0, 0.18080000000000002),
        (0.1, 0.5, 0.5, 0.20483276469913345),
        (0.5, 2.0, 3.0, 0.6875),
        (0.9, 5.0, 1.5, 0.7761721343162157),
        (0.35, 10.0, 20.0, 0.592386663663905),
        (0.999, 3.0, 0.5, 0.9407468104840537),
        (1e-5, 0.5, 40.0, 0.02249424725410791),
        (0.5, 100.0, 100.0, 0.5),
        (0.01, 0.5, 100.0, 0.843224134575559),
        (0.05, 0.5, 5.0, 0.515208786901603),
        (0.97, 0.5, 5.0, 0.999999993943697),
        (0.42, 250.0, 0.5, 3.035444351144309e-96),
        (0.9, 0.5, 0.001, 0.0036287029439785383),
    ];
    #[test]
    fn inc_beta_matches_high_precision_reference() {
        for (x, a, b, expected) in INC_BETA_REFERENCE.iter() {
            assert_rel_close(
                inc_beta(*x, *a, *b),
                *expected,
                1e-10,
                &format!("I_{x}({a}, {b})"),
            );
        }
    }
    #[test]
    fn inc_beta_uniform_case_is_identity() {
        // Beta(1, 1) is the uniform distribution: I_x(1, 1) = x.
        for &x in &[1e-8_f64, 0.01, 0.3, 0.5, 0.77, 0.999] {
            assert!(
                (inc_beta(x, 1.0, 1.0) - x).abs() < 1e-13,
                "I_{x}(1,1) = {}",
                inc_beta(x, 1.0, 1.0)
            );
        }
    }
    #[test]
    fn inc_beta_integer_parameters_match_binomial_sums() {
        // I_x(2, 3) = sum_{j=2}^{4} C(4, j) x^j (1-x)^(4-j); I_x(2, 2) = x^2 (3 - 2x).
        let x = 0.2_f64;
        let binomial =
            6.0 * x.powi(2) * (1.0 - x).powi(2) + 4.0 * x.powi(3) * (1.0 - x) + x.powi(4);
        assert!((binomial - 0.1808).abs() < 1e-15);
        assert!((inc_beta(0.2, 2.0, 3.0) - 0.1808).abs() < 1e-12);
        for &x in &[0.1_f64, 0.4, 0.5, 0.8] {
            let exact = x * x * (3.0 - 2.0 * x);
            assert!((inc_beta(x, 2.0, 2.0) - exact).abs() < 1e-13, "x = {x}");
        }
    }
    #[test]
    fn inc_beta_reflection_matches_reference_on_both_sides() {
        // Points on both sides of the (a+1)/(a+b+2) switch must agree with the
        // reflection identity I_x(a, b) = 1 - I_{1-x}(b, a) *and* the reference.
        for &(x, a, b) in &[
            (0.3_f64, 2.0_f64, 3.0_f64),
            (0.6, 2.0, 3.0),
            (0.9, 5.0, 1.5),
        ] {
            let direct = inc_beta(x, a, b);
            let reflected = 1.0 - inc_beta(1.0 - x, b, a);
            assert!((direct - reflected).abs() < 1e-12, "x={x} a={a} b={b}");
        }
        assert!((inc_beta(0.3, 2.0, 3.0) - 0.3483).abs() < 1e-12);
        assert!((inc_beta(0.6, 2.0, 3.0) - 0.8208).abs() < 1e-12);
    }
    #[test]
    fn inc_beta_boundaries_and_invalid_input() {
        assert_eq!(inc_beta(0.0, 2.0, 3.0), 0.0);
        assert_eq!(inc_beta(-0.5, 2.0, 3.0), 0.0);
        assert_eq!(inc_beta(1.0, 2.0, 3.0), 1.0);
        assert_eq!(inc_beta(1.5, 2.0, 3.0), 1.0);
        assert!(inc_beta(f64::NAN, 2.0, 3.0).is_nan());
        assert!(inc_beta(0.5, 0.0, 3.0).is_nan());
        assert!(inc_beta(0.5, 2.0, -1.0).is_nan());
    }
    #[test]
    fn inc_beta_symmetric_shape_is_one_half_at_midpoint() {
        // Beta(a, a) is symmetric about 1/2, so I_{1/2}(a, a) = 1/2 for every a.
        for &a in &[0.3_f64, 0.5, 1.0, 2.0, 3.5, 10.0, 50.0, 400.0] {
            assert!(
                (inc_beta(0.5, a, a) - 0.5).abs() < 1e-12,
                "I_0.5({a}, {a}) = {}",
                inc_beta(0.5, a, a)
            );
        }
    }
    #[test]
    fn inc_beta_unit_second_shape_is_a_power() {
        // Beta(a, 1) has CDF x^a, on both sides of the (a+1)/(a+3) switch point.
        for &a in &[0.5_f64, 1.0, 2.0, 3.7, 12.0] {
            for &x in &[0.001_f64, 0.1, 0.4, 0.7, 0.95, 0.999] {
                assert_rel_close(
                    inc_beta(x, a, 1.0),
                    x.powf(a),
                    1e-11,
                    &format!("I_{x}({a}, 1)"),
                );
            }
        }
    }
    #[test]
    fn inc_beta_unit_first_shape_is_complementary_power() {
        // Beta(1, b) has CDF 1 - (1-x)^b.
        for &b in &[0.5_f64, 1.0, 2.0, 3.7, 12.0] {
            for &x in &[0.001_f64, 0.1, 0.4, 0.7, 0.95] {
                let exact = 1.0 - (1.0 - x).powf(b);
                assert!(
                    (inc_beta(x, 1.0, b) - exact).abs() < 1e-12,
                    "I_{x}(1, {b}) = {} vs {exact}",
                    inc_beta(x, 1.0, b)
                );
            }
        }
    }
    #[test]
    fn inc_beta_never_leaves_the_unit_interval() {
        // The old continued fraction returned values such as -0.31 for I_0.5(2, 2).
        for &(a, b) in &[
            (2.0_f64, 2.0_f64),
            (1.0, 1.0),
            (2.0, 3.0),
            (0.5, 4.0),
            (7.0, 0.8),
        ] {
            for i in 0..=100 {
                let x = i as f64 / 100.0;
                let p = inc_beta(x, a, b);
                assert!((0.0..=1.0).contains(&p), "I_{x}({a}, {b}) = {p}");
            }
        }
    }
    #[test]
    fn inc_beta_is_monotone_in_x() {
        let mut previous = 0.0_f64;
        for i in 1..200 {
            let x = i as f64 / 200.0;
            let value = inc_beta(x, 3.5, 1.25);
            assert!(value >= previous, "CDF must be non-decreasing at x={x}");
            previous = value;
        }
        assert!(previous > 0.99);
    }
}
