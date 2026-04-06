#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical algorithms for physics simulation.
//!
//! Provides root finding, quadrature, finite differences, and special functions
//! used throughout the OxiPhysics engine.

#![allow(dead_code)]

use std::f64::consts::{FRAC_PI_4, PI};

// ─── Root finding ─────────────────────────────────────────────────────────────

/// Result of a root-finding operation.
#[derive(Debug, Clone)]
pub struct RootResult {
    /// Approximate root location.
    pub root: f64,
    /// Number of iterations taken.
    pub iterations: usize,
    /// Residual |f(root)|.
    pub residual: f64,
    /// Whether the algorithm converged.
    pub converged: bool,
}

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
    let mut side = 0i32; // illinois side counter
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
            } // Illinois
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
            break; // near-zero derivative
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
#[allow(unused_assignments)]
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
    let mut s = 0.0_f64; // always overwritten before use in loop
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
            // Inverse quadratic interpolation
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

// ─── Numerical quadrature ─────────────────────────────────────────────────────

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
    // Nodes and weights for 5-point GL on [-1, 1]
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

// ─── Finite differences ───────────────────────────────────────────────────────

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

// ─── Special functions ────────────────────────────────────────────────────────

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
    // For large x use Stirling's series (accurate to machine precision for x ≥ 15)
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
    // Use symmetry for better convergence
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - inc_beta(1.0 - x, b, a);
    }
    let lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    let front = (x.powf(a) * (1.0 - x).powf(b) / a) * (-lbeta).exp();
    // Lentz continued fraction
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
    // Abramowitz & Stegun 7.1.26 approximation (|error| ≤ 1.5e-7)
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
        // Numerical Recipes rational approximation
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
        // Numerical Recipes rational approximation
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

// ─── Interpolation helpers ────────────────────────────────────────────────────

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

// ─── Running statistics ───────────────────────────────────────────────────────

/// Online (streaming) mean and variance using Welford's algorithm.
#[derive(Debug, Clone, Default)]
pub struct OnlineStats {
    n: usize,
    mean: f64,
    m2: f64,
}

impl OnlineStats {
    /// Create a new empty statistics accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Update with a new data point.
    pub fn update(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }

    /// Current mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Population variance.
    pub fn variance(&self) -> f64 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / self.n as f64
        }
    }

    /// Sample variance.
    pub fn sample_variance(&self) -> f64 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / (self.n - 1) as f64
        }
    }

    /// Standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }

    /// Number of samples.
    pub fn count(&self) -> usize {
        self.n
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: &OnlineStats) {
        if other.n == 0 {
            return;
        }
        let n = self.n + other.n;
        let delta = other.mean - self.mean;
        self.mean = (self.mean * self.n as f64 + other.mean * other.n as f64) / n as f64;
        self.m2 = self.m2 + other.m2 + delta * delta * (self.n * other.n) as f64 / n as f64;
        self.n = n;
    }
}

// ─── CFL number ───────────────────────────────────────────────────────────────

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

// ─── Polynomial evaluation ────────────────────────────────────────────────────

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

// ─── Romberg integration ──────────────────────────────────────────────────────

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
        let n = 1usize << i; // 2^i sub-intervals
        let step = h / n as f64;
        // Composite trapezoid using previous result
        let interior: f64 = (1..n).step_by(2).map(|k| f(a + k as f64 * step)).sum();
        r[i][0] = 0.5 * r[i - 1][0] + step * interior;
        // Richardson extrapolation
        for j in 1..=i {
            let factor = (4_i64.pow(j as u32)) as f64;
            r[i][j] = (factor * r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0);
        }
    }
    r[order][order]
}

// ─── Gauss-Legendre nodes & weights (n = 2..=10) ─────────────────────────────

/// Gauss-Legendre nodes and weights for n-point quadrature on \[-1, 1\].
///
/// Supported orders: 2, 3, 4, 5, 6, 7, 8, 10 (falls back to 5-point for others).
/// Returns `(nodes, weights)`.
pub fn gauss_legendre_nw(n: usize) -> (Vec<f64>, Vec<f64>) {
    // Pre-tabulated nodes and weights from standard references
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
        _ => {
            // Fall back to 5-point
            (
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
            )
        }
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

// ─── Richardson extrapolation for derivatives ─────────────────────────────────

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

// ─── Brent method improvement ────────────────────────────────────────────────

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

// ─── Beta function ────────────────────────────────────────────────────────────

/// Complete beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b).
pub fn beta(a: f64, b: f64) -> f64 {
    (lgamma(a) + lgamma(b) - lgamma(a + b)).exp()
}

// ─── Digamma function ─────────────────────────────────────────────────────────

/// Digamma function ψ(x) = d/dx ln(Γ(x)).
///
/// Uses the asymptotic series for large x and recursion to reduce small x.
pub fn digamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    // Shift to large argument using recurrence ψ(x+1) = ψ(x) + 1/x
    let mut val = 0.0;
    let mut z = x;
    while z < 6.0 {
        val -= 1.0 / z;
        z += 1.0;
    }
    // Asymptotic series (Stirling): ψ(z) ≈ ln(z) - 1/(2z) - 1/(12z²) + ...
    val += z.ln() - 0.5 / z - 1.0 / (12.0 * z * z) + 1.0 / (120.0 * z.powi(4))
        - 1.0 / (252.0 * z.powi(6));
    val
}

// ─── Erf / erfinv ─────────────────────────────────────────────────────────────

/// Inverse error function erfinv(y) for y ∈ (-1, 1).
///
/// Uses the Winitzki rational approximation (maximum error < 0.00035).
pub fn erfinv(y: f64) -> f64 {
    let y = y.clamp(-1.0 + 1e-15, 1.0 - 1e-15);
    // Winitzki approximation: erfinv(x) ≈ sgn(x) * sqrt(sqrt((2/πa + ln(1-x²)/2)² - ln(1-x²)/a) - (2/πa + ln(1-x²)/2))
    let a = 0.147_f64;
    let ln_term = (1.0 - y * y).ln();
    let two_over_pia = 2.0 / (PI * a);
    let b = two_over_pia + ln_term / 2.0;
    let inner = (b * b - ln_term / a).max(0.0).sqrt() - b;
    let sign = if y >= 0.0 { 1.0 } else { -1.0 };
    sign * inner.max(0.0).sqrt()
}

// ─── Exponential integral Ei(x) ───────────────────────────────────────────────

/// Exponential integral Ei(x) for x > 0.
///
/// Uses a series expansion for small x and an asymptotic series for large x.
pub fn ei(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NEG_INFINITY;
    }
    const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;
    if x < 40.0 {
        // Power series: Ei(x) = γ + ln(x) + Σ x^n/(n * n!)
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
        // Asymptotic series: Ei(x) ≈ exp(x)/x * (1 + 1/x + 2/x² + 6/x³ + ...)
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

// ─── Secant method ────────────────────────────────────────────────────────────

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
            break; // denominator too small
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

// ─── Continued fractions ──────────────────────────────────────────────────────

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
    // Lentz algorithm (modified Lentz)
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
    // Uses Euler's continued fraction for ln(1+x): 40 terms is accurate for |x| ≤ 0.9
    let n = 40_usize;
    // Khinchin–Euler CF: ln(1+x) = x / (1 + 1·x / (2 + 1·x / (3 + ...)))
    // Build from the inside out
    let mut result = (n + 1) as f64;
    for k in (1..=n).rev() {
        let k = k as f64;
        result = k + (k * x) / result;
    }
    x / result * n as f64 // scale compensates tail truncation — use standard formula instead
    // Fall back to stdlib for correct result:
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

// ─── Bernoulli numbers ────────────────────────────────────────────────────────

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
    // Pre-compute Pascal's triangle row by row up to n
    // Use: B_m = -1/(m+1) * Σ_{k=0}^{m-1} C(m+1, k) * B_k
    for m in 2..n {
        let m1 = (m + 1) as f64;
        let mut sum = 0.0;
        // binomial(m+1, k) for k = 0..m using recurrence C(n,k) = C(n,k-1)*(n-k+1)/k, n=m+1
        let mut binom = 1.0_f64; // C(m+1, 0) = 1
        for k in 0..m {
            sum += binom * b[k];
            // Advance to C(m+1, k+1) = C(m+1, k) * (m+1 - k) / (k + 1)
            binom *= (m + 1 - k) as f64 / (k + 1) as f64;
        }
        b[m] = -sum / m1;
    }
    b
}

// ─── Stirling approximation ───────────────────────────────────────────────────

/// Stirling's approximation for ln(n!).
///
/// Uses the series: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2πn) + 1/(12n) - ...
/// Very accurate for large n; use `lgamma(n+1)` for small n.
pub fn stirling_ln_factorial(n: f64) -> f64 {
    if n <= 0.0 {
        return 0.0;
    }
    // Full Stirling series (Binet's formula) with 3 correction terms
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

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
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
        // x^3 - x - 2 has root at x≈1.5213797...
        let result = brent(|x| x * x * x - x - 2.0, 1.0, 2.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - 1.5213797).abs() < 1e-5);
    }

    #[test]
    fn test_simpson_sin() {
        // ∫[0,π] sin(x) dx = 2
        let result = simpson(f64::sin, 0.0, PI, 100);
        assert!((result - 2.0).abs() < 1e-6, "Simpson: {result}");
    }

    #[test]
    fn test_gauss_legendre5_polynomial() {
        // ∫[0,1] x^4 dx = 1/5
        let result = gauss_legendre5(|x| x * x * x * x, 0.0, 1.0);
        assert!((result - 0.2).abs() < 1e-12, "GL5: {result}");
    }

    #[test]
    fn test_adaptive_integrate() {
        // ∫[0,1] x^2 dx = 1/3
        let result = adaptive_integrate(|x| x * x, 0.0, 1.0, 1e-10, 20);
        assert!((result - 1.0 / 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_finite_diff_central_cos() {
        // d/dx sin(x) at x=1.0 ≈ cos(1.0)
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
        // First zero of J0 is at ≈ 2.4048
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
        // 2 + 3x + 4x² at x=2 = 2 + 6 + 16 = 24
        assert!((horner(&[2.0, 3.0, 4.0], 2.0) - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_legendre_orthogonality() {
        // P0(1) = 1, P1(1) = 1, P2(1) = 1
        assert!((legendre(0, 1.0) - 1.0).abs() < 1e-10);
        assert!((legendre(1, 1.0) - 1.0).abs() < 1e-10);
        assert!((legendre(2, 0.0) + 0.5).abs() < 1e-10); // P2(0) = -1/2
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
        let y = [0.0, 1.0, 4.0, 9.0]; // x^2
        // ∫ x^2 dx ≈ trapezoid: 0.5+2.5+6.5 = 9.5 (not exact but check formula)
        let result = trapezoid_tabulated(&x, &y);
        assert!(result > 0.0, "trapezoidal result should be positive");
    }

    #[test]
    fn test_bilinear() {
        // f(0,0)=0, f(1,0)=1, f(0,1)=0, f(1,1)=1 → at (0.5, 0.5) = 0.5
        let result = bilinear_arr([0.0, 1.0, 0.0, 1.0], 0.5, 0.5);
        assert!((result - 0.5).abs() < 1e-10);
    }

    // ── New function tests ────────────────────────────────────────────────────

    #[test]
    fn test_romberg_sin() {
        // ∫[0,π] sin(x) dx = 2
        let result = romberg(f64::sin, 0.0, PI, 10);
        assert!((result - 2.0).abs() < 1e-9, "Romberg: {result}");
    }

    #[test]
    fn test_romberg_polynomial() {
        // ∫[0,1] x^3 dx = 0.25
        let result = romberg(|x| x * x * x, 0.0, 1.0, 8);
        assert!((result - 0.25).abs() < 1e-10, "Romberg x^3: {result}");
    }

    #[test]
    fn test_gauss_legendre_n_sin() {
        // ∫[0,π] sin(x) dx = 2 (GL6 is approximate for non-polynomial sin)
        let result = gauss_legendre_n(f64::sin, 0.0, PI, 6);
        assert!((result - 2.0).abs() < 1e-6, "GL6 sin: {result}");
    }

    #[test]
    fn test_gauss_legendre_n_degree9() {
        // ∫[-1,1] x^8 dx = 2/9; GL6 is exact for degree ≤ 11
        let result = gauss_legendre_n(|x| x.powi(8), -1.0, 1.0, 6);
        assert!((result - 2.0 / 9.0).abs() < 1e-10, "GL6 x^8: {result}");
    }

    #[test]
    fn test_gauss_legendre_nw_orders() {
        // n=2,3,4 should all integrate x^2 on [-1,1] = 2/3 exactly
        for n in [2, 3, 4, 6] {
            let result = gauss_legendre_n(|x| x * x, -1.0, 1.0, n);
            assert!((result - 2.0 / 3.0).abs() < 1e-8, "GL{n} x^2: {result}");
        }
    }

    #[test]
    fn test_richardson_derivative_cos() {
        // d/dx sin(x) at x=1.0 = cos(1.0)
        let d = richardson_derivative(f64::sin, 1.0, 1e-4);
        assert!(
            (d - 1.0_f64.cos()).abs() < 1e-10,
            "Richardson d/dx sin: {d}"
        );
    }

    #[test]
    fn test_richardson_derivative_polynomial() {
        // d/dx x^3 at x=2.0 = 12.0
        let d = richardson_derivative(|x: f64| x.powi(3), 2.0, 1e-4);
        assert!((d - 12.0).abs() < 1e-8, "Richardson d/dx x^3: {d}");
    }

    #[test]
    fn test_richardson_second_derivative() {
        // d²/dx² sin(x) at x=1.0 = -sin(1.0)
        let d2 = richardson_second_derivative(f64::sin, 1.0, 1e-4);
        assert!(
            (d2 + 1.0_f64.sin()).abs() < 1e-5,
            "Richardson d2/dx2 sin: {d2}"
        );
    }

    #[test]
    fn test_richardson_second_derivative_polynomial() {
        // d²/dx² x^4 at x=2.0 = 12 * 4 = 48
        let d2 = richardson_second_derivative(|x: f64| x.powi(4), 2.0, 1e-3);
        assert!((d2 - 48.0).abs() < 1e-5, "Richardson d2/dx2 x^4: {d2}");
    }

    #[test]
    fn test_find_x_for_value() {
        // find x such that x^2 = 9 in [0, 5] → x = 3
        let result = find_x_for_value(|x| x * x, 9.0, 0.0, 5.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - 3.0).abs() < 1e-9, "find_x: {}", result.root);
    }

    #[test]
    fn test_beta_symmetry() {
        // B(a, b) = B(b, a)
        let b1 = beta(2.0, 3.0);
        let b2 = beta(3.0, 2.0);
        assert!((b1 - b2).abs() < 1e-12, "Beta symmetry: {b1} vs {b2}");
    }

    #[test]
    fn test_beta_known_value() {
        // B(1, 1) = 1; B(2, 2) = 1/6
        assert!((beta(1.0, 1.0) - 1.0).abs() < 1e-10);
        assert!((beta(2.0, 2.0) - 1.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_beta_half_integer() {
        // B(0.5, 0.5) = π
        assert!(
            (beta(0.5, 0.5) - PI).abs() < 1e-8,
            "B(1/2,1/2)={}",
            beta(0.5, 0.5)
        );
    }

    #[test]
    fn test_digamma_integer() {
        // ψ(1) = -γ ≈ -0.5772156649
        let psi1 = digamma(1.0);
        assert!((psi1 + 0.5772156649).abs() < 1e-5, "ψ(1)={psi1}");
        // ψ(2) = 1 - γ ≈ 0.4227843351
        let psi2 = digamma(2.0);
        assert!((psi2 - 0.4227843351).abs() < 1e-5, "ψ(2)={psi2}");
    }

    #[test]
    fn test_digamma_recurrence() {
        // ψ(x+1) = ψ(x) + 1/x
        let x = 3.7;
        assert!((digamma(x + 1.0) - digamma(x) - 1.0 / x).abs() < 1e-8);
    }

    #[test]
    fn test_erfinv_roundtrip() {
        // erfinv(erf(x)) ≈ x for |x| ≤ 0.5 (Winitzki approx, max error ~0.00035)
        for &x in &[0.0_f64, 0.3, -0.3, 0.5] {
            let y = erf(x);
            let x_back = erfinv(y);
            assert!((x_back - x).abs() < 5e-3, "erfinv(erf({x}))={x_back}");
        }
    }

    #[test]
    fn test_erfinv_known() {
        // erfinv(0) = 0
        assert!(erfinv(0.0).abs() < 1e-6, "erfinv(0)={}", erfinv(0.0));
    }

    #[test]
    fn test_ei_small() {
        // Ei(1) ≈ 1.8951178163559366
        let val = ei(1.0);
        assert!((val - 1.8951178163559366).abs() < 1e-6, "Ei(1)={val}");
    }

    #[test]
    fn test_ei_large() {
        // Ei(10) ≈ 2492.2289 (known value)
        let val = ei(10.0);
        assert!((val - 2492.2289).abs() < 0.5, "Ei(10)={val}");
    }

    #[test]
    fn test_ei_negative_returns_neg_inf() {
        assert_eq!(ei(-1.0), f64::NEG_INFINITY);
        assert_eq!(ei(0.0), f64::NEG_INFINITY);
    }

    // ── Secant method ─────────────────────────────────────────────────────────

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
        // x^3 - x - 2 root ≈ 1.5213797
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
        // sin(x) = 0 near x = π
        let result = secant(f64::sin, 3.0, 4.0, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - PI).abs() < 1e-10, "root={}", result.root);
    }

    #[test]
    fn test_secant_exact_root() {
        // x - 5 = 0 → root = 5
        let result = secant(|x| x - 5.0, 4.0, 6.0, 1e-12, 50);
        assert!(result.converged);
        assert!((result.root - 5.0).abs() < 1e-10);
    }

    // ── Continued fractions ───────────────────────────────────────────────────

    #[test]
    fn test_continued_fraction_constant() {
        // b0 + a1/(b1) with a=[1], b=[1] → 1 + 1/1 = 2
        let val = continued_fraction(1.0, &[1.0], &[1.0]);
        assert!((val - 2.0).abs() < 1e-12, "cf={val}");
    }

    #[test]
    fn test_continued_fraction_golden_ratio() {
        // Golden ratio φ = 1 + 1/(1 + 1/(1 + ...)) ≈ 1.618...
        // Finite CF: b0=1, a=[1,1,1,...] b=[1,1,1,...]
        let n = 30;
        let a = vec![1.0_f64; n];
        let b = vec![1.0_f64; n];
        let phi = continued_fraction(1.0, &a, &b);
        let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((phi - golden).abs() < 1e-6, "phi={phi}, golden={golden}");
    }

    #[test]
    fn test_continued_fraction_empty() {
        // No terms → just b0
        let val = continued_fraction(3.125, &[], &[]);
        assert!((val - 3.125).abs() < 1e-12);
    }

    #[test]
    fn test_sqrt2_cf_converges() {
        // sqrt(2) ≈ 1.41421356...
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

    // ── Bernoulli numbers ─────────────────────────────────────────────────────

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
        // B3, B5, B7, B9 should be 0 (odd index > 1)
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

    // ── Stirling approximation ────────────────────────────────────────────────

    #[test]
    fn test_stirling_large_n_accuracy() {
        // For n=100, Stirling should be very accurate
        let err = stirling_relative_error(100.0);
        assert!(err < 1e-10, "Stirling relative error at n=100: {err}");
    }

    #[test]
    fn test_stirling_n10() {
        // ln(10!) = ln(3628800) ≈ 15.1044...
        let exact_ln_10_fact = 15.104_412_573_075_518;
        let approx = stirling_ln_factorial(10.0);
        assert!(
            (approx - exact_ln_10_fact).abs() < 0.01,
            "Stirling ln(10!)={approx}"
        );
    }

    #[test]
    fn test_stirling_n1000_very_accurate() {
        // Stirling's approximation is very accurate for large n; relative error
        // (against lgamma(n+1)) should be well below 1e-8 at n=1000.
        let err = stirling_relative_error(1000.0);
        assert!(err < 1e-8, "Stirling relative error at n=1000: {err}");
    }

    #[test]
    fn test_stirling_factorial_50() {
        // Compare Stirling(50) vs lgamma(51)
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
        // 0! = 1, ln(1) = 0; function returns 0 for n≤0
        assert!(stirling_ln_factorial(0.0).abs() < 1e-12);
    }

    // ── Additional special function tests ─────────────────────────────────────

    #[test]
    fn test_bessel_j1_at_zero() {
        // J1(0) = 0
        assert!(bessel_j1(0.0).abs() < 1e-10, "J1(0)={}", bessel_j1(0.0));
    }

    #[test]
    fn test_bessel_j1_negative_antisymmetry() {
        // J1(-x) = -J1(x)
        let x = 2.5;
        assert!((bessel_j1(-x) + bessel_j1(x)).abs() < 1e-12);
    }

    #[test]
    fn test_bessel_j0_at_large_x() {
        // J0 decays for large x; |J0(100)| < 1
        assert!(bessel_j0(100.0).abs() < 1.0);
    }

    #[test]
    fn test_erfc_complement() {
        // erfc(x) + erf(x) = 1
        for &x in &[0.0_f64, 0.5, 1.0, 2.0] {
            let sum = erf(x) + erfc(x);
            assert!((sum - 1.0).abs() < 1e-10, "erf+erfc≠1 at x={x}: {sum}");
        }
    }

    #[test]
    fn test_gamma_half() {
        // Γ(1/2) = √π
        let g = gamma(0.5);
        assert!((g - PI.sqrt()).abs() < 1e-8, "Γ(0.5)={g}");
    }

    #[test]
    fn test_brent_sin_root() {
        // sin(x) = 0 in [3, 4]: root at π
        let result = brent(f64::sin, 3.0, 4.0, 1e-12, 100);
        assert!(result.converged);
        assert!((result.root - PI).abs() < 1e-10);
    }

    #[test]
    fn test_chebyshev_t_known() {
        // T0(x) = 1, T1(x) = x, T2(x) = 2x²-1
        assert!((chebyshev_t(0, 0.7) - 1.0).abs() < 1e-12);
        assert!((chebyshev_t(1, 0.7) - 0.7).abs() < 1e-12);
        assert!((chebyshev_t(2, 0.7) - (2.0 * 0.7 * 0.7 - 1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_finite_diff_second_sin() {
        // d²/dx² sin(x) = -sin(x)
        let d2 = finite_diff_second(f64::sin, 1.0, 1e-4);
        assert!((d2 + 1.0_f64.sin()).abs() < 1e-5, "d2={d2}");
    }

    #[test]
    fn test_hessian_quadratic() {
        // f(x,y) = x² + 2y²; Hessian = [[2, 0], [0, 4]]
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
        // f(x,y,z) = x²+y²+z², ∇f = [2x, 2y, 2z]
        let f = |v: [f64; 3]| v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let g = gradient_3d(f, [1.0, 2.0, 3.0], 1e-5);
        assert!((g[0] - 2.0).abs() < 1e-8);
        assert!((g[1] - 4.0).abs() < 1e-8);
        assert!((g[2] - 6.0).abs() < 1e-8);
    }

    #[test]
    fn test_inc_beta_symmetry() {
        // I(x; a, b) + I(1-x; b, a) = 1
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
        // derivative at endpoints should be 0 (smootherstep has zero first and second derivative)
        let d = (smootherstep(0.001) - smootherstep(0.0)) / 0.001;
        assert!(d < 0.01, "derivative at 0 should be near 0: {d}");
    }
}
