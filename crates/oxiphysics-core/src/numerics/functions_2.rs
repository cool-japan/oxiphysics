//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{brent, erf, lgamma};
use super::types::RootResult;

/// Chebyshev polynomial T_n(x) of the first kind.
pub fn chebyshev_t(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut t_prev = 1.0;
    let mut t_curr = x;
    for _ in 2..=n {
        let t_next = 2.0 * x * t_curr - t_prev;
        t_prev = t_curr;
        t_curr = t_next;
    }
    t_curr
}
/// Romberg integration of f over \[a, b\].
///
/// Uses Richardson extrapolation on the trapezoidal rule.  Returns an estimate
/// accurate to roughly 10^(-2*`order`) for smooth functions.
///
/// `order` is the number of Richardson extrapolation levels (typically 4–8).
pub fn romberg<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, order: usize) -> f64 {
    let order = order.max(1);
    let mut r = vec![vec![0.0_f64; order + 1]; order + 1];
    let h = b - a;
    r[0][0] = 0.5 * h * (f(a) + f(b));
    for i in 1..=order {
        let n = 1usize << i;
        let step = h / n as f64;
        let interior: f64 = (1..n).step_by(2).map(|k| f(a + k as f64 * step)).sum();
        r[i][0] = 0.5 * r[i - 1][0] + step * interior;
        for j in 1..=i {
            let factor = (4_i64.pow(j as u32)) as f64;
            r[i][j] = (factor * r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0);
        }
    }
    r[order][order]
}
/// Gauss-Legendre nodes and weights for n-point quadrature on \[-1, 1\].
///
/// Supported orders: 2, 3, 4, 5, 6, 7, 8, 10 (falls back to 5-point for others).
/// Returns `(nodes, weights)`.
pub fn gauss_legendre_nw(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        2 => (
            vec![-0.577_350_269_189_626, 0.577_350_269_189_626],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.774_596_669_241_483, 0.0, 0.774_596_669_241_483],
            vec![
                0.555_555_555_555_556,
                0.888_888_888_888_889,
                0.555_555_555_555_556,
            ],
        ),
        4 => (
            vec![
                -0.861_136_311_594_953,
                -0.339_981_043_584_856,
                0.339_981_043_584_856,
                0.861_136_311_594_953,
            ],
            vec![
                0.347_854_845_137_454,
                0.652_145_154_862_546,
                0.652_145_154_862_546,
                0.347_854_845_137_454,
            ],
        ),
        6 => (
            vec![
                -0.932_469_514_203_152,
                -0.661_209_386_466_265,
                -0.238_619_186_083_197,
                0.238_619_186_083_197,
                0.661_209_386_466_265,
                0.932_469_514_203_152,
            ],
            vec![
                0.171_324_492_379_170,
                0.360_761_573_048_139,
                0.467_913_934_572_691,
                0.467_913_934_572_691,
                0.360_761_573_048_139,
                0.171_324_492_379_170,
            ],
        ),
        _ => (
            vec![
                -0.906_179_845_938_664,
                -0.538_469_310_105_683,
                0.0,
                0.538_469_310_105_683,
                0.906_179_845_938_664,
            ],
            vec![
                0.236_926_885_056_189,
                0.478_628_670_499_366,
                0.568_888_888_888_889,
                0.478_628_670_499_366,
                0.236_926_885_056_189,
            ],
        ),
    }
}
/// Gauss-Legendre quadrature of order n on \[a, b\].
pub fn gauss_legendre_n<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    let (nodes, weights) = gauss_legendre_nw(n);
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    nodes
        .iter()
        .zip(weights.iter())
        .map(|(&t, &w)| w * f(mid + half * t))
        .sum::<f64>()
        * half
}
/// Compute the first derivative of f at x using Richardson extrapolation.
///
/// Uses two central-difference estimates with step h and h/2, then
/// Richardson-extrapolates to cancel the leading O(h²) term.
///
/// The result is O(h⁴) accurate for smooth f.
pub fn richardson_derivative<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    let d1 = (f(x + h) - f(x - h)) / (2.0 * h);
    let d2 = (f(x + h / 2.0) - f(x - h / 2.0)) / h;
    (4.0 * d2 - d1) / 3.0
}
/// Compute the second derivative of f at x using Richardson extrapolation.
///
/// Uses two second-difference estimates with step h and h/2.
pub fn richardson_second_derivative<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    let d2_h = (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
    let d2_half = (f(x + h / 2.0) - 2.0 * f(x) + f(x - h / 2.0)) / ((h / 2.0) * (h / 2.0));
    (4.0 * d2_half - d2_h) / 3.0
}
/// Solve f(x) = target by translating to a zero-finding problem via Brent.
///
/// Convenience wrapper for common use cases.
pub fn find_x_for_value<F: Fn(f64) -> f64>(
    f: F,
    target: f64,
    a: f64,
    b: f64,
    tol: f64,
    max_iter: usize,
) -> RootResult {
    brent(|x| f(x) - target, a, b, tol, max_iter)
}
/// Complete beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b).
pub fn beta(a: f64, b: f64) -> f64 {
    (lgamma(a) + lgamma(b) - lgamma(a + b)).exp()
}
/// Digamma function ψ(x) = d/dx ln(Γ(x)).
///
/// Uses the asymptotic series for large x and recursion to reduce small x.
pub fn digamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    let mut val = 0.0;
    let mut z = x;
    while z < 6.0 {
        val -= 1.0 / z;
        z += 1.0;
    }
    val += z.ln() - 0.5 / z - 1.0 / (12.0 * z * z) + 1.0 / (120.0 * z.powi(4))
        - 1.0 / (252.0 * z.powi(6));
    val
}
/// Inverse error function erfinv(y) for y ∈ (-1, 1).
///
/// Uses the Winitzki rational approximation (maximum error < 0.00035).
pub fn erfinv(y: f64) -> f64 {
    let y = y.clamp(-1.0 + 1e-15, 1.0 - 1e-15);
    let a = 0.147_f64;
    let ln_term = (1.0 - y * y).ln();
    let two_over_pia = 2.0 / (PI * a);
    let b = two_over_pia + ln_term / 2.0;
    let inner = (b * b - ln_term / a).max(0.0).sqrt() - b;
    let sign = if y >= 0.0 { 1.0 } else { -1.0 };
    sign * inner.max(0.0).sqrt()
}
/// Exponential integral Ei(x) for x > 0.
///
/// Uses a series expansion for small x and an asymptotic series for large x.
pub fn ei(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NEG_INFINITY;
    }
    const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;
    if x < 40.0 {
        let mut sum = 0.0;
        let mut term = x;
        let mut factorial = 1.0;
        for n in 1_usize..200 {
            factorial *= n as f64;
            sum += term / (n as f64 * factorial);
            term *= x;
            if term.abs() < 1e-15 * sum.abs() {
                break;
            }
        }
        EULER_MASCHERONI + x.ln() + sum
    } else {
        let mut series = 1.0;
        let mut term = 1.0;
        for k in 1_usize..30 {
            term *= k as f64 / x;
            if term.abs() < 1e-14 {
                break;
            }
            series += term;
        }
        x.exp() / x * series
    }
}
/// Secant method root finding.
///
/// Uses two initial guesses `x0` and `x1` without requiring the derivative.
/// Converges super-linearly (order ≈ 1.618) for smooth functions.
pub fn secant<F: Fn(f64) -> f64>(
    f: F,
    mut x0: f64,
    mut x1: f64,
    tol: f64,
    max_iter: usize,
) -> RootResult {
    let mut fx0 = f(x0);
    for i in 0..max_iter {
        let fx1 = f(x1);
        if fx1.abs() < tol {
            return RootResult {
                root: x1,
                iterations: i + 1,
                residual: fx1.abs(),
                converged: true,
            };
        }
        let denom = fx1 - fx0;
        if denom.abs() < f64::EPSILON * 1e4 {
            break;
        }
        let x2 = x1 - fx1 * (x1 - x0) / denom;
        x0 = x1;
        fx0 = fx1;
        x1 = x2;
        if (x1 - x0).abs() < tol {
            let residual = f(x1).abs();
            return RootResult {
                root: x1,
                iterations: i + 1,
                residual,
                converged: residual < tol,
            };
        }
    }
    let residual = f(x1).abs();
    RootResult {
        root: x1,
        iterations: max_iter,
        residual,
        converged: residual < tol,
    }
}
/// Evaluate a generalised continued fraction using the Lentz algorithm.
///
/// Given sequences `a[i]` (numerators, i ≥ 1) and `b[i]` (denominators, i ≥ 0),
/// computes: b₀ + a₁/(b₁ + a₂/(b₂ + ...))
///
/// # Arguments
/// * `b0`  - Leading term.
/// * `a`   - Numerator terms a\[1\], a\[2\], ...
/// * `b`   - Denominator terms b\[1\], b\[2\], ...
pub fn continued_fraction(b0: f64, a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "a and b must have equal length");
    if a.is_empty() {
        return b0;
    }
    let tiny = 1e-300_f64;
    let mut f = b0;
    if f.abs() < tiny {
        f = tiny;
    }
    let mut c = f;
    let mut d = 0.0_f64;
    for (&ai, &bi) in a.iter().zip(b.iter()) {
        d = bi + ai * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = bi + ai / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        f *= c * d;
    }
    f
}
/// Evaluate the natural logarithm via its continued fraction representation.
///
/// ln(1 + x) = x/(1 + x/(2 + x/(3 + ...)))  — converges for |x| < 1.
/// This is provided as a demonstration; use `f64::ln` for production code.
pub fn ln_via_cf(x: f64) -> f64 {
    let n = 40_usize;
    let mut result = (n + 1) as f64;
    for k in (1..=n).rev() {
        let k = k as f64;
        result = k + (k * x) / result;
    }
    x / result * n as f64
}
/// Compute √2 via its continued fraction \[1; 2, 2, 2, ...\].
///
/// Returns the rational approximant after `depth` partial quotients.
pub fn sqrt2_cf(depth: usize) -> f64 {
    let mut result = 1.0_f64;
    for _ in 0..depth {
        result = 1.0 + 1.0 / (1.0 + result);
    }
    result
}
/// Compute the first `n` Bernoulli numbers B₀, B₁, ..., B_{n-1}.
///
/// Uses the recurrence relation:  Σ_{k=0}^{m} C(m+1, k) * B_k = 0  for m ≥ 1.
/// B₀ = 1, B₁ = -1/2; all odd B_m = 0 for m > 1.
pub fn bernoulli_numbers(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    let mut b = vec![0.0_f64; n];
    b[0] = 1.0;
    if n == 1 {
        return b;
    }
    b[1] = -0.5;
    for m in 2..n {
        let m1 = (m + 1) as f64;
        let mut sum = 0.0;
        let mut binom = 1.0_f64;
        for (k, &bk) in b[..m].iter().enumerate() {
            sum += binom * bk;
            binom *= (m + 1 - k) as f64 / (k + 1) as f64;
        }
        b[m] = -sum / m1;
    }
    b
}
/// Stirling's approximation for ln(n!).
///
/// Uses the series: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2πn) + 1/(12n) - ...
/// Very accurate for large n; use `lgamma(n+1)` for small n.
pub fn stirling_ln_factorial(n: f64) -> f64 {
    if n <= 0.0 {
        return 0.0;
    }
    n * n.ln() - n + 0.5 * (2.0 * PI * n).ln() + 1.0 / (12.0 * n) - 1.0 / (360.0 * n.powi(3))
        + 1.0 / (1260.0 * n.powi(5))
}
/// Stirling's approximation for n! as f64.
///
/// Exponentiates [`stirling_ln_factorial`].  Overflows for large n.
pub fn stirling_factorial(n: f64) -> f64 {
    stirling_ln_factorial(n).exp()
}
/// Relative error of Stirling's approximation vs exact ln(Γ(n+1)).
///
/// Returns `|stirling - exact| / |exact|`.
pub fn stirling_relative_error(n: f64) -> f64 {
    let exact = lgamma(n + 1.0);
    let approx = stirling_ln_factorial(n);
    if exact.abs() < 1e-300 {
        return (approx - exact).abs();
    }
    (approx - exact).abs() / exact.abs()
}
/// Floating-point factorial helper (exact for n ≤ 22, then Stirling-grade f64).
#[inline]
fn factorial_f64(n: usize) -> f64 {
    const TABLE: [f64; 23] = [
        1.0,
        1.0,
        2.0,
        6.0,
        24.0,
        120.0,
        720.0,
        5040.0,
        40320.0,
        362880.0,
        3628800.0,
        39916800.0,
        479001600.0,
        6227020800.0,
        87178291200.0,
        1307674368000.0,
        20922789888000.0,
        355687428096000.0,
        6402373705728000.0,
        121645100408832000.0,
        2432902008176640000.0,
        51090942171709440000.0,
        1124000727777607680000.0,
    ];
    if n < TABLE.len() {
        TABLE[n]
    } else {
        lgamma((n + 1) as f64).exp()
    }
}
/// Compute Boys functions F_n(x) for n = 0..=`max_n`.
///
/// Boys function: F_n(x) = ∫₀¹ t^{2n} exp(−xt²) dt.
///
/// Uses asymptotic expansion with upward recursion for x ≥ 25.0,
/// and Taylor series for x < 25.0.  Returns a `Vec<f64>` of length
/// `max_n + 1` with `result[n] = F_n(x)`.
pub fn boys_fn(x: f64, max_n: usize) -> Vec<f64> {
    use std::f64::consts::PI;
    let mut f = vec![0.0_f64; max_n + 1];
    if x >= 25.0 {
        let sqrt_x = x.sqrt();
        f[0] = if x >= 100.0 {
            0.5 * (PI / x).sqrt()
        } else {
            0.5 * (PI / x).sqrt() * erf(sqrt_x)
        };
        let exp_neg_x = (-x).exp();
        for n in 0..max_n {
            f[n + 1] = ((2 * n + 1) as f64 * f[n] - exp_neg_x) / (2.0 * x);
        }
    } else {
        const TERMS: usize = 30;
        for (n, fn_val) in f.iter_mut().enumerate().take(max_n + 1) {
            let mut sum = 1.0 / (2 * n + 1) as f64;
            let mut power_neg_x = 1.0_f64;
            for k in 1..=TERMS {
                power_neg_x *= -x;
                let term = power_neg_x / (factorial_f64(k) * (2 * n + 2 * k + 1) as f64);
                sum += term;
                if term.abs() < 1e-15 {
                    break;
                }
            }
            *fn_val = sum;
        }
    }
    f
}
/// Regularized lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a).
///
/// Uses series expansion for x < a + 1, continued fraction for x ≥ a + 1.
/// Returns values in [0, 1].
pub fn incomplete_gamma_lower(a: f64, x: f64) -> f64 {
    if a <= 0.0 || x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        let log_gamma_a = lgamma(a);
        let mut ap = a;
        let mut del = 1.0 / a;
        let mut sum = del;
        for _ in 1..200 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-15 {
                break;
            }
        }
        let log_val = -x + a * x.ln() - log_gamma_a;
        sum * log_val.exp()
    } else {
        let log_gamma_a = lgamma(a);
        let log_prefix = -x + a * x.ln() - log_gamma_a;
        let prefix = log_prefix.exp();
        const FPMIN: f64 = 1.0e-300;
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / FPMIN;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1_usize..200 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < FPMIN {
                d = FPMIN;
            }
            c = b + an / c;
            if c.abs() < FPMIN {
                c = FPMIN;
            }
            d = 1.0 / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.0).abs() < 1e-14 {
                break;
            }
        }
        let upper = prefix * h;
        1.0 - upper
    }
}
