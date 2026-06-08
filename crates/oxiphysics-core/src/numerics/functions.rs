//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::{FRAC_PI_4, PI};

use super::types::RootResult;

/// Bisection root finding on interval \[a, b\].
///
/// Requires `f(a)` and `f(b)` to have opposite signs.
/// Converges in at most `max_iter` iterations.
pub fn bisect<F: Fn(f64) -> f64>(
    f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> RootResult {
    let mut fa = f(a);
    let mut mid = a;
    for i in 0..max_iter {
        mid = 0.5 * (a + b);
        let fm = f(mid);
        if fm.abs() < tol || (b - a) < tol {
            return RootResult {
                root: mid,
                iterations: i + 1,
                residual: fm.abs(),
                converged: true,
            };
        }
        if fa * fm < 0.0 {
            b = mid;
        } else {
            a = mid;
            fa = fm;
        }
    }
    RootResult {
        root: mid,
        iterations: max_iter,
        residual: f(mid).abs(),
        converged: false,
    }
}
/// Regula falsi (Illinois variant) root finding.
///
/// Generally faster than bisection while maintaining bracketing.
pub fn regula_falsi<F: Fn(f64) -> f64>(
    f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> RootResult {
    let mut fa = f(a);
    let mut fb = f(b);
    let mut side = 0i32;
    for i in 0..max_iter {
        let c = (a * fb - b * fa) / (fb - fa);
        let fc = f(c);
        if fc.abs() < tol {
            return RootResult {
                root: c,
                iterations: i + 1,
                residual: fc.abs(),
                converged: true,
            };
        }
        if fa * fc < 0.0 {
            b = c;
            fb = fc;
            if side == 1 {
                fa *= 0.5;
            }
            side = 1;
        } else {
            a = c;
            fa = fc;
            if side == -1 {
                fb *= 0.5;
            }
            side = -1;
        }
        if (b - a).abs() < tol {
            return RootResult {
                root: c,
                iterations: i + 1,
                residual: fc.abs(),
                converged: true,
            };
        }
    }
    let c = (a * fb - b * fa) / (fb - fa);
    RootResult {
        root: c,
        iterations: max_iter,
        residual: f(c).abs(),
        converged: false,
    }
}
/// Newton-Raphson root finding.
///
/// Requires a good initial guess and the derivative `df`.
pub fn newton<F, DF>(f: F, df: DF, x0: f64, tol: f64, max_iter: usize) -> RootResult
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    let mut x = x0;
    for i in 0..max_iter {
        let fx = f(x);
        if fx.abs() < tol {
            return RootResult {
                root: x,
                iterations: i,
                residual: fx.abs(),
                converged: true,
            };
        }
        let dfx = df(x);
        if dfx.abs() < f64::EPSILON * 100.0 {
            break;
        }
        x -= fx / dfx;
    }
    let residual = f(x).abs();
    RootResult {
        root: x,
        iterations: max_iter,
        residual,
        converged: residual < tol,
    }
}
/// Brent's method — robust root finding combining bisection, secant, and inverse quadratic interpolation.
pub fn brent<F: Fn(f64) -> f64>(
    f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> RootResult {
    let mut fa = f(a);
    let mut fb = f(b);
    if fa * fb > 0.0 {
        return RootResult {
            root: a,
            iterations: 0,
            residual: fa.abs(),
            converged: false,
        };
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = a;
    let mut fc = fa;
    let mut s: f64;
    let mut mflag = true;
    let mut d = 0.0;
    for i in 0..max_iter {
        if fb.abs() < tol || (b - a).abs() < tol {
            return RootResult {
                root: b,
                iterations: i,
                residual: fb.abs(),
                converged: true,
            };
        }
        if (fa - fc).abs() > f64::EPSILON && (fb - fc).abs() > f64::EPSILON {
            s = a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa);
        }
        let cond1 = !(s > (3.0 * a + b) / 4.0 && s < b || s < (3.0 * a + b) / 4.0 && s > b);
        let cond2 = mflag && (s - b).abs() >= (b - c).abs() / 2.0;
        let cond3 = !mflag && (s - b).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (b - c).abs() < tol;
        let cond5 = !mflag && (c - d).abs() < tol;
        if cond1 || cond2 || cond3 || cond4 || cond5 {
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }
        let fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    RootResult {
        root: b,
        iterations: max_iter,
        residual: fb.abs(),
        converged: fb.abs() < tol,
    }
}
/// Simpson's rule integration of f over \[a, b\] with n subintervals (n must be even).
pub fn simpson<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    let n = if n % 2 == 1 { n + 1 } else { n };
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        sum += if i % 2 == 0 { 2.0 * f(x) } else { 4.0 * f(x) };
    }
    sum * h / 3.0
}
/// Gauss-Legendre 5-point quadrature on \[-1, 1\], mapped to \[a, b\].
pub fn gauss_legendre5<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> f64 {
    const NODES: [f64; 5] = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683,
        0.0,
        0.538_469_310_105_683,
        0.906_179_845_938_664,
    ];
    const WEIGHTS: [f64; 5] = [
        0.236_926_885_056_189,
        0.478_628_670_499_366,
        0.568_888_888_888_889,
        0.478_628_670_499_366,
        0.236_926_885_056_189,
    ];
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    NODES
        .iter()
        .zip(WEIGHTS.iter())
        .map(|(&t, &w)| w * f(mid + half * t))
        .sum::<f64>()
        * half
}
/// Adaptive Gaussian quadrature (recursive).
///
/// Integrates f over \[a,b\] to relative tolerance `tol`.
pub fn adaptive_integrate<F: Fn(f64) -> f64 + Clone>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
    depth: usize,
) -> f64 {
    let whole = gauss_legendre5(f.clone(), a, b);
    let mid = 0.5 * (a + b);
    let left = gauss_legendre5(f.clone(), a, mid);
    let right = gauss_legendre5(f.clone(), mid, b);
    let error = (left + right - whole).abs();
    if error < 15.0 * tol || depth == 0 {
        left + right + (whole - left - right) / 15.0
    } else {
        adaptive_integrate(f.clone(), a, mid, tol * 0.5, depth - 1)
            + adaptive_integrate(f, mid, b, tol * 0.5, depth - 1)
    }
}
/// Trapezoidal rule for tabulated data.
pub fn trapezoid_tabulated(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len());
    x.windows(2)
        .zip(y.windows(2))
        .map(|(xi, yi)| (xi[1] - xi[0]) * (yi[0] + yi[1]) * 0.5)
        .sum()
}
/// First derivative of f at x using centered differences with step h.
pub fn finite_diff_central<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    (f(x + h) - f(x - h)) / (2.0 * h)
}
/// Second derivative of f at x using centered differences.
pub fn finite_diff_second<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
}
/// Gradient of f: R^3 → R using central differences.
pub fn gradient_3d<F: Fn([f64; 3]) -> f64>(f: F, x: [f64; 3], h: f64) -> [f64; 3] {
    let fx = x;
    let gx = {
        let mut p = fx;
        p[0] += h;
        let mut m = fx;
        m[0] -= h;
        (f(p) - f(m)) / (2.0 * h)
    };
    let gy = {
        let mut p = fx;
        p[1] += h;
        let mut m = fx;
        m[1] -= h;
        (f(p) - f(m)) / (2.0 * h)
    };
    let gz = {
        let mut p = fx;
        p[2] += h;
        let mut m = fx;
        m[2] -= h;
        (f(p) - f(m)) / (2.0 * h)
    };
    [gx, gy, gz]
}
/// Hessian of f: R^n → R as a flat n×n matrix using central differences.
pub fn hessian<F: Fn(&[f64]) -> f64>(f: &F, x: &[f64], h: f64) -> Vec<Vec<f64>> {
    let n = x.len();
    let mut h_mat = vec![vec![0.0; n]; n];
    let f0 = f(x);
    for i in 0..n {
        for j in i..n {
            let mut pp = x.to_vec();
            pp[i] += h;
            pp[j] += h;
            let mut pm = x.to_vec();
            pm[i] += h;
            pm[j] -= h;
            let mut mp = x.to_vec();
            mp[i] -= h;
            mp[j] += h;
            let mut mm = x.to_vec();
            mm[i] -= h;
            mm[j] -= h;
            let val = if i == j {
                let mut xph = x.to_vec();
                xph[i] += h;
                let mut xmh = x.to_vec();
                xmh[i] -= h;
                (f(&xph) - 2.0 * f0 + f(&xmh)) / (h * h)
            } else {
                (f(&pp) - f(&pm) - f(&mp) + f(&mm)) / (4.0 * h * h)
            };
            h_mat[i][j] = val;
            h_mat[j][i] = val;
        }
    }
    h_mat
}
/// Gamma function approximation (Lanczos, g=7).
pub fn gamma(x: f64) -> f64 {
    if x < 0.5 {
        PI / ((PI * x).sin() * gamma(1.0 - x))
    } else {
        let x = x - 1.0;
        const G: f64 = 7.0;
        const C: [f64; 9] = [
            0.999_999_999_999_809_9,
            676.520_368_121_885_1,
            -1_259.139_216_722_402_8,
            771.323_428_777_653_1,
            -176.615_029_162_140_6,
            12.507_343_278_686_905,
            -0.138_571_095_265_720_12,
            9.984_369_578_019_572e-6,
            1.505_632_735_149_312e-7,
        ];
        let t = x + G + 0.5;
        let ser: f64 = C[1..]
            .iter()
            .enumerate()
            .map(|(k, &c)| c / (x + k as f64 + 1.0))
            .sum::<f64>()
            + C[0];
        (2.0 * PI).sqrt() * ser * t.powf(x + 0.5) * (-t).exp()
    }
}
/// Log gamma function.
pub fn lgamma(x: f64) -> f64 {
    if x <= 0.0 {
        return gamma(x).abs().ln();
    }
    if x < 15.0 {
        return gamma(x).abs().ln();
    }
    let z = x;
    z.ln() * (z - 0.5) - z + 0.5 * (2.0 * std::f64::consts::PI).ln() + 1.0 / (12.0 * z)
        - 1.0 / (360.0 * z.powi(3))
        + 1.0 / (1260.0 * z.powi(5))
        - 1.0 / (1680.0 * z.powi(7))
}
/// Factorial n! as f64 (using gamma).
pub fn factorial(n: u64) -> f64 {
    gamma(n as f64 + 1.0)
}
/// Regularized incomplete beta function I(x; a, b).
///
/// Uses continued fraction expansion (Lentz method).
pub fn inc_beta(x: f64, a: f64, b: f64) -> f64 {
    if !(0.0..=1.0).contains(&x) {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - inc_beta(1.0 - x, b, a);
    }
    let lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    let front = (x.powf(a) * (1.0 - x).powf(b) / a) * (-lbeta).exp();
    let mut c = 1.0;
    let mut d = 1.0 - (a + b) * x / (a + 1.0);
    if d.abs() < 1e-300 {
        d = 1e-300;
    }
    d = 1.0 / d;
    let mut cf = d;
    for m in 1_usize..200 {
        let mf = m as f64;
        for step in 0..2 {
            let num = if step == 0 {
                -(a + mf) * (a + b + mf) * x / ((a + 2.0 * mf - 1.0) * (a + 2.0 * mf))
            } else {
                mf * (b - mf) * x / ((a + 2.0 * mf - 1.0) * (a + 2.0 * mf))
            };
            d = 1.0 + num * d;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            d = 1.0 / d;
            c = 1.0 + num / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            let delta = c * d;
            cf *= delta;
            if (delta - 1.0).abs() < 1e-14 {
                break;
            }
        }
    }
    front * cf
}
/// Error function erf(x) using a rational approximation.
pub fn erf(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    sign * (1.0 - poly * (-x * x).exp())
}
/// Complementary error function erfc(x) = 1 - erf(x).
pub fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}
/// Normal CDF Φ(x).
pub fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}
/// Bessel function J0(x) (zeroth order, first kind).
pub fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let y = x * x;
        let p1 = 57568490574.0
            + y * (-13362590354.0
                + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
        let q1 = 57568490411.0
            + y * (1029532985.0
                + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
        p1 / q1
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - FRAC_PI_4;
        let p1 = 1.0
            + y * (-0.001098628627
                + y * (0.000002734511 + y * (-0.000000020761 + y * 0.000000000206)));
        let q1 = -0.01562499995
            + y * (0.000001430488
                + y * (-0.000000006911 + y * (0.000000000077 + y * -0.000000000001)));
        (2.0 / (PI * ax)).sqrt() * (xx.cos() * p1 - z * xx.sin() * q1)
    }
}
/// Bessel function J1(x) (first order, first kind).
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let y = x * x;
        let p = x
            * (72362614232.0
                + y * (-7895059235.0
                    + y * (242396853.1
                        + y * (-2972611.439 + y * (15704.48260 + y * (-30.16307633))))));
        let q = 144725228442.0
            + y * (2300535178.0
                + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
        p / q
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 2.356_194_491;
        let p1 =
            1.0 + y * (0.00183105 + y * (-0.0000766479 + y * (0.000000567 + y * -0.0000000058)));
        let q1 = 0.04687499995 + y * (-0.00000204 + y * (0.0000000869 + y * -0.000000000395));
        let ans = (2.0 / (PI * ax)).sqrt() * (xx.cos() * p1 - z * xx.sin() * q1);
        if x < 0.0 { -ans } else { ans }
    }
}
/// Linear interpolation scalar: `a + t*(b-a)`.
#[inline]
pub fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}
/// Smooth Hermite interpolation (smoothstep).
#[inline]
pub fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}
/// Smoother step (Ken Perlin's 6th degree).
#[inline]
pub fn smootherstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}
/// Bilinear interpolation on a 2×2 grid (array form).
///
/// `values` = \[f(0,0), f(1,0), f(0,1), f(1,1)\], `tx` and `ty` in \[0,1\].
pub fn bilinear_arr(values: [f64; 4], tx: f64, ty: f64) -> f64 {
    let bottom = lerp_scalar(values[0], values[1], tx);
    let top = lerp_scalar(values[2], values[3], tx);
    lerp_scalar(bottom, top, ty)
}
/// Trilinear interpolation on a 2×2×2 grid (array form).
///
/// `values` indexed as \[z0y0x0, z0y0x1, z0y1x0, z0y1x1, z1y0x0, z1y0x1, z1y1x0, z1y1x1\].
pub fn trilinear_arr(values: [f64; 8], tx: f64, ty: f64, tz: f64) -> f64 {
    let c00 = lerp_scalar(values[0], values[1], tx);
    let c10 = lerp_scalar(values[2], values[3], tx);
    let c01 = lerp_scalar(values[4], values[5], tx);
    let c11 = lerp_scalar(values[6], values[7], tx);
    let c0 = lerp_scalar(c00, c10, ty);
    let c1 = lerp_scalar(c01, c11, ty);
    lerp_scalar(c0, c1, tz)
}
/// Compute the CFL (Courant–Friedrichs–Lewy) number.
///
/// CFL = max_velocity * dt / min_cell_size
pub fn cfl_number(max_velocity: f64, dt: f64, min_cell_size: f64) -> f64 {
    if min_cell_size < f64::EPSILON {
        return f64::INFINITY;
    }
    max_velocity * dt / min_cell_size
}
/// Compute the maximum safe time step for a given CFL target.
pub fn dt_from_cfl(cfl_target: f64, max_velocity: f64, min_cell_size: f64) -> f64 {
    if max_velocity < f64::EPSILON {
        return f64::MAX;
    }
    cfl_target * min_cell_size / max_velocity
}
/// Evaluate a polynomial using Horner's method.
///
/// coeffs = \[a0, a1, a2, ...\] for a0 + a1*x + a2*x² + ...
pub fn horner(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}
/// Evaluate a Legendre polynomial P_n(x) via recurrence.
pub fn legendre(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut p_prev = 1.0;
    let mut p_curr = x;
    for k in 2..=n {
        let k = k as f64;
        let p_next = ((2.0 * k - 1.0) * x * p_curr - (k - 1.0) * p_prev) / k;
        p_prev = p_curr;
        p_curr = p_next;
    }
    p_curr
}
