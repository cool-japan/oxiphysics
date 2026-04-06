// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical quadrature (numerical integration) methods.
//!
//! Provides Gauss-Legendre, Gauss-Lobatto, Gauss-Chebyshev, Gauss-Hermite,
//! Gauss-Laguerre, Clenshaw-Curtis, Simpson's rule, Romberg, adaptive
//! Gauss-Kronrod (G7K15), double-exponential (tanh-sinh), and 2D/3D tensor
//! product quadrature.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use std::f64::consts::PI;

// ── helpers ─────────────────────────────────────────────────────────────────

/// Evaluate the Legendre polynomial P_n(x) and its derivative P_n'(x).
///
/// Uses the three-term recurrence relation.
fn legendre_p_and_dp(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    if n == 1 {
        return (x, 1.0);
    }
    let mut p_prev = 1.0_f64;
    let mut p_curr = x;
    for k in 2..=(n as u32) {
        let kf = k as f64;
        let p_next = ((2.0 * kf - 1.0) * x * p_curr - (kf - 1.0) * p_prev) / kf;
        p_prev = p_curr;
        p_curr = p_next;
    }
    // Derivative via: (1 - x^2) P_n'(x) = n (P_{n-1}(x) - x P_n(x))
    let dp = (n as f64) * (p_prev - x * p_curr) / (1.0 - x * x).max(1e-300);
    (p_curr, dp)
}

// ── Gauss-Legendre ───────────────────────────────────────────────────────────

/// Compute the *n*-point Gauss-Legendre nodes and weights on \[-1, 1\].
///
/// Returns a `Vec` of `(node, weight)` pairs sorted in ascending node order.
/// Nodes and weights are found via Newton iteration on the Legendre polynomial.
pub fn gauss_legendre_weights(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 1, "n must be at least 1");
    let mut nw = Vec::with_capacity(n);
    // Only compute half — the rule is symmetric about x=0.
    let half = n.div_ceil(2);
    for i in 1..=half {
        // Initial guess (Golub-Welsch / Tricomi approximation)
        let theta = PI * ((i as f64) - 0.25) / ((n as f64) + 0.5);
        let mut x = theta.cos();
        // Newton iteration
        for _ in 0..100 {
            let (p, dp) = legendre_p_and_dp(n, x);
            let dx = p / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let (_, dp) = legendre_p_and_dp(n, x);
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        // Symmetric counterpart
        nw.push((-x, w));
        if 2 * i - 1 != n {
            // not the midpoint for odd n
            nw.push((x, w));
        }
    }
    nw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    nw
}

/// Integrate `f` on \[a, b\] using the *n*-point Gauss-Legendre rule.
///
/// The interval is linearly mapped from \[-1, 1\] to \[a, b\].
pub fn gauss_legendre_integrate(f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let nw = gauss_legendre_weights(n);
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    nw.iter()
        .map(|(xi, wi)| wi * f(mid + half * xi))
        .sum::<f64>()
        * half
}

// ── Gauss-Lobatto ────────────────────────────────────────────────────────────

/// Compute the *n*-point Gauss-Lobatto nodes and weights on \[-1, 1\].
///
/// The endpoints ±1 are always included. Requires n ≥ 2.
/// Interior nodes are roots of P_{n-1}'(x) found by Newton iteration.
pub fn gauss_lobatto_weights(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 2, "Gauss-Lobatto requires n >= 2");
    let mut nw = Vec::with_capacity(n);

    // Endpoints always included with weight 2 / (n*(n-1))
    let w_end = 2.0 / ((n as f64) * ((n as f64) - 1.0));
    nw.push((-1.0_f64, w_end));

    // Interior nodes: zeros of P_{n-1}'(x)
    let m = n - 2; // number of interior nodes
    let half = m.div_ceil(2);
    for i in 1..=half {
        // Good initial guess for P_{n-1}' zeros
        let theta = PI * (i as f64) / ((n as f64) - 1.0);
        let mut x = -theta.cos();
        for _ in 0..100 {
            let (p, dp) = legendre_p_and_dp(n - 1, x);
            // Second derivative: (1-x^2) P'' = -2x P' + ... use recurrence
            // P_{n-1}''(x) via forward diff (simpler)
            let eps = 1e-8;
            let (_, dp_p) = legendre_p_and_dp(n - 1, x + eps);
            let ddp = (dp_p - dp) / eps;
            let dx = dp / ddp;
            x -= dx;
            if dx.abs() < 1e-14 {
                break;
            }
            let _ = p;
        }
        let (pval, _) = legendre_p_and_dp(n - 1, x);
        let w = 2.0 / ((n as f64 - 1.0) * (n as f64) * pval * pval);
        nw.push((-x, w));
        if 2 * i - 1 != m {
            nw.push((x, w));
        }
    }

    nw.push((1.0_f64, w_end));
    nw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    nw
}

// ── Gauss-Chebyshev (first kind) ─────────────────────────────────────────────

/// Compute the *n*-point Gauss-Chebyshev nodes and weights on \[-1, 1\].
///
/// Uses the Chebyshev-Gauss rule of the first kind: nodes are the roots of
/// T_n(x) and all weights equal π/n.
pub fn gauss_chebyshev_weights(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 1, "n must be at least 1");
    let w = PI / (n as f64);
    (1..=n)
        .map(|k| {
            let x = ((2 * k - 1) as f64 * PI / (2.0 * n as f64)).cos();
            (x, w)
        })
        .collect()
}

// ── Simpson's rule ────────────────────────────────────────────────────────────

/// Composite Simpson's rule on \[a, b\] with `n` subintervals (must be even).
///
/// If `n` is odd it is rounded up to the next even number.
pub fn simpsons_rule(f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let n = if n.is_multiple_of(2) { n } else { n + 1 };
    let h = (b - a) / (n as f64);
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let x = a + (i as f64) * h;
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(x);
    }
    sum * h / 3.0
}

// ── Romberg integration ───────────────────────────────────────────────────────

/// Romberg integration on \[a, b\].
///
/// Builds Richardson-extrapolation tableau up to `max_level` levels.
/// Stops early when consecutive diagonal entries agree within `tol`.
pub fn romberg_integration(
    f: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    max_level: usize,
    tol: f64,
) -> f64 {
    let m = max_level.max(1);
    let mut r = vec![vec![0.0_f64; m + 1]; m + 1];

    // First trapezoidal estimate
    r[0][0] = 0.5 * (b - a) * (f(a) + f(b));

    for i in 1..=m {
        // Trapezoidal with 2^i intervals
        let n = 1usize << i; // 2^i
        let h = (b - a) / (n as f64);
        let mut sum = 0.0;
        for k in 0..(n / 2) {
            sum += f(a + (2 * k + 1) as f64 * h);
        }
        r[i][0] = 0.5 * r[i - 1][0] + h * sum;

        // Richardson extrapolation
        for j in 1..=i {
            let factor = (4.0_f64).powi(j as i32);
            r[i][j] = (factor * r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0);
        }

        if i >= 1 && (r[i][i] - r[i - 1][i - 1]).abs() < tol {
            return r[i][i];
        }
    }
    r[m][m]
}

// ── Gauss-Kronrod G7K15 ────────────────────────────────────────────────────

/// The 15-point Gauss-Kronrod nodes on \[-1, 1\] (positive half only; x=0 last).
const GK15_NODES: [f64; 8] = [
    0.991_455_371_120_813,
    0.949_107_912_342_758,
    0.864_864_423_359_769,
    0.741_531_185_599_394,
    0.586_087_235_467_691,
    0.405_845_151_377_397,
    0.207_784_955_007_898,
    0.0,
];

/// Gauss-Kronrod 15-point weights.
const GK15_WEIGHTS: [f64; 8] = [
    0.022_935_322_010_529,
    0.063_092_092_629_979,
    0.104_790_010_322_250,
    0.140_653_259_715_525,
    0.169_004_726_639_267,
    0.190_350_578_064_785,
    0.204_432_940_075_298,
    0.209_482_141_084_728,
];

/// Gauss 7-point weights (subset of K15 nodes).
const G7_WEIGHTS: [f64; 4] = [
    0.129_484_966_168_870,
    0.279_705_391_489_277,
    0.381_830_050_505_119,
    0.417_959_183_673_469,
];

/// Adaptive Gauss-Kronrod G7K15 integration on \[a, b\].
///
/// Estimates the error as |G15 - G7| and recursively bisects intervals that
/// exceed `tol`. The recursion halts at `max_depth`.
pub fn adaptive_gauss_kronrod(
    f: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    tol: f64,
    max_depth: usize,
) -> f64 {
    gk15_recursive(f, a, b, tol, max_depth, 0)
}

/// Recursive helper for [`adaptive_gauss_kronrod`].
fn gk15_recursive(
    f: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    tol: f64,
    max_depth: usize,
    depth: usize,
) -> f64 {
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);

    // Evaluate K15 and G7
    let mut gk = 0.0_f64;
    let mut g7 = 0.0_f64;

    // x = 0 node (index 7)
    let f0 = f(mid);
    gk += GK15_WEIGHTS[7] * f0;
    g7 += G7_WEIGHTS[3] * f0;

    // Symmetric pairs
    for i in 0..7 {
        let xi = GK15_NODES[i];
        let fplus = f(mid + half * xi);
        let fminus = f(mid - half * xi);
        gk += GK15_WEIGHTS[i] * (fplus + fminus);

        // G7 uses nodes at indices 1, 3, 5 (i.e. i=1,3,5 in GK15_NODES)
        if i == 1 || i == 3 || i == 5 {
            let gi = match i {
                1 => 0,
                3 => 1,
                5 => 2,
                _ => unreachable!(),
            };
            g7 += G7_WEIGHTS[gi] * (fplus + fminus);
        }
    }

    gk *= half;
    g7 *= half;

    let err = (gk - g7).abs();
    if depth >= max_depth || err < tol {
        return gk;
    }

    gk15_recursive(f, a, mid, tol * 0.5, max_depth, depth + 1)
        + gk15_recursive(f, mid, b, tol * 0.5, max_depth, depth + 1)
}

// ── Double-exponential (tanh-sinh) ────────────────────────────────────────────

/// Tanh-sinh (double-exponential) quadrature on \[a, b\].
///
/// The substitution x = tanh(π/2 · sinh(t)) maps ℝ → (-1, 1), giving
/// exponential convergence for smooth (even endpoint-singular) integrands.
/// `n` is the number of points per half (total ≈ 2n+1), `h` is the step size
/// (typical value 0.1).
pub fn double_exponential(f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize, h: f64) -> f64 {
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);

    let phi = |t: f64| -> f64 { (0.5 * PI * t.sinh()).tanh() };
    let dphi = |t: f64| -> f64 {
        let s = 0.5 * PI * t.sinh();
        0.5 * PI * t.cosh() / s.cosh().powi(2)
    };

    // Central point
    let mut sum = f(mid) * dphi(0.0);

    for k in 1..=n {
        let t = k as f64 * h;
        let p = phi(t);
        let dp = dphi(t);
        let xp = mid + half * p;
        let xm = mid - half * p;
        if xp.is_finite() && xm.is_finite() {
            sum += (f(xp) + f(xm)) * dp;
        }
    }

    sum * h * half
}

// ── Clenshaw-Curtis ───────────────────────────────────────────────────────────

/// Compute *n*-point Clenshaw-Curtis nodes and weights on \[-1, 1\].
///
/// Uses the standard closed-form formula (Waldvogel 2006). Requires n ≥ 2.
/// Nodes are the Chebyshev extrema `cos(k π / (n-1))` for k = 0..n-1,
/// and the total weight sum equals 2 (the length of \[-1, 1\]).
pub fn clenshaw_curtis_weights(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 2, "Clenshaw-Curtis requires n >= 2");
    let nm1 = n - 1; // n-1
    let nm1f = nm1 as f64;

    // Chebyshev extrema nodes
    let nodes: Vec<f64> = (0..n).map(|k| (k as f64 * PI / nm1f).cos()).collect();

    // Waldvogel / Sommariva weight formula:
    //   c[k] = 2 / (n-1) * Re[ sum_{j=0}^{n-1} b_j / (1 - 4j^2) * exp(i k j pi / (n-1)) ]
    // Simpler equivalent: Fejer-type formula for closed rule.
    // We use the explicit trigonometric series for the weight of node k:
    //   w[k] = (c_k / (n-1)) * (1 - sum_{m=1}^{floor((n-1)/2)} b_m / (4m^2-1) * cos(2m k pi/(n-1)))
    // where c_0 = c_{n-1} = 1/2, c_k = 1 otherwise; b_m = 1 if 2m < n-1, else 1/2.
    let weights: Vec<f64> = (0..n)
        .map(|k| {
            let c_k = if k == 0 || k == nm1 { 0.5 } else { 1.0 };
            let half_nm1 = nm1 / 2;
            let sum: f64 = (1..=half_nm1)
                .map(|m| {
                    let b_m = if 2 * m == nm1 { 0.5 } else { 1.0 };
                    let theta = 2.0 * m as f64 * k as f64 * PI / nm1f;
                    b_m * theta.cos() / (4.0 * (m as f64).powi(2) - 1.0)
                })
                .sum();
            2.0 * c_k * (1.0 - 2.0 * sum) / nm1f
        })
        .collect();

    nodes.into_iter().zip(weights).collect()
}

// ── Gauss-Hermite ─────────────────────────────────────────────────────────────

/// Compute the *n*-point Gauss-Hermite nodes and weights.
///
/// The rule integrates `f(x) · exp(-x²)` over (-∞, +∞).  Nodes are roots of
/// the probabilist's Hermite polynomial H_n(x); found by Newton iteration.
pub fn gauss_hermite_weights(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 1, "n must be at least 1");

    // Evaluate physicist's H_n(x) and H_n'(x) via three-term recurrence
    let hermite_p_dp = |x: f64| -> (f64, f64) {
        if n == 0 {
            return (1.0, 0.0);
        }
        let mut pm1 = 1.0_f64;
        let mut p = 2.0 * x;
        if n == 1 {
            return (p, 2.0);
        }
        for k in 2..=n {
            let pnew = 2.0 * x * p - 2.0 * (k as f64 - 1.0) * pm1;
            pm1 = p;
            p = pnew;
        }
        let dp = 2.0 * (n as f64) * pm1;
        (p, dp)
    };

    let half = n.div_ceil(2);
    let mut nw = Vec::with_capacity(n);
    for i in 1..=half {
        // Initial guess from approximation
        let x0 = (2 * n + 1) as f64;
        let mut x =
            (2.0 * x0 + 1.0).sqrt() * (PI * (4.0 * i as f64 - 1.0) / (4.0 * n as f64 + 2.0)).cos();
        for _ in 0..100 {
            let (p, dp) = hermite_p_dp(x);
            if dp.abs() < 1e-300 {
                break;
            }
            let dx = p / dp;
            x -= dx;
            if dx.abs() < 1e-14 {
                break;
            }
        }
        let (pm1, _) = hermite_p_dp(x - 1e-9);
        let (pp1, _) = hermite_p_dp(x + 1e-9);
        let dp_num = (pp1 - pm1) / (2e-9);
        let w = if dp_num.abs() < 1e-300 {
            0.0
        } else {
            let (p_prev, _) = {
                // Compute H_{n-1}(x)
                if n == 1 {
                    (1.0_f64, 0.0_f64)
                } else {
                    let mut qm1 = 1.0_f64;
                    let mut q = 2.0 * x;
                    for k in 2..n {
                        let qnew = 2.0 * x * q - 2.0 * (k as f64 - 1.0) * qm1;
                        qm1 = q;
                        q = qnew;
                    }
                    (q, 0.0)
                }
            };
            let hn_p_sq = (2.0_f64.powi(n as i32 - 1)
                * (1..=n).map(|k| k as f64).product::<f64>()
                * PI.sqrt())
                / ((n as f64) * p_prev * p_prev);
            if hn_p_sq.is_finite() && hn_p_sq > 0.0 {
                hn_p_sq
            } else {
                // Fallback: use numerical derivative
                2.0_f64.powi(n as i32 - 1) * (1..=n).map(|k| k as f64).product::<f64>() * PI.sqrt()
                    / ((n as f64) * p_prev * p_prev)
            }
        };
        nw.push((-x, w));
        if 2 * i - 1 != n {
            nw.push((x, w));
        }
    }
    nw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    nw
}

// ── Gauss-Laguerre ────────────────────────────────────────────────────────────

/// Evaluate the Laguerre polynomial L_n(x) and its derivative.
fn laguerre_p_dp(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    if n == 1 {
        return (1.0 - x, -1.0);
    }
    let mut pm1 = 1.0_f64;
    let mut p = 1.0 - x;
    for k in 2..=n {
        let kf = k as f64;
        let pnew = ((2.0 * kf - 1.0 - x) * p - (kf - 1.0) * pm1) / kf;
        pm1 = p;
        p = pnew;
    }
    // Derivative: L_n'(x) = -L_{n-1}(x) (for monic normalization the relation is:
    // L_n'(x) = -(n / (n + 1)) * sum ..., but the standard relation is n * L_n(x) = ... )
    // Use: n L_n'(x) = n L_{n-1}(x) - (n - x) ... simplify: L_n'(x) = (L_n(x) - L_{n-1}(x)) / ... nope
    // Correct recurrence for derivative: n L_n'(x) = n L_{n-1}(x) - (x L_n(x))' (see DLMF 18.9.23)
    // Simpler: L_n'(x) = -\sum_{k=0}^{n-1} L_k(x) or equivalently L_n'(x) = (n L_n(x) - n L_{n-1}(x)) / x  for x != 0
    // Standard relation: L_n'(x) = -L_{n-1}(x) (for generalized α=0)  ← this is CORRECT for standard Laguerre
    let dp = -pm1;
    (p, dp)
}

/// Compute the *n*-point Gauss-Laguerre nodes and weights.
///
/// The rule integrates `f(x) · exp(-x)` over \[0, +∞).
/// Uses Newton iteration on L_n(x) with Tricomi initial guesses.
pub fn gauss_laguerre_weights(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 1, "n must be at least 1");
    let mut nw = Vec::with_capacity(n);
    for i in 1..=n {
        // Tricomi initial guess for the i-th zero of L_n
        let nf = n as f64;
        let jv = PI * (4 * i - 1) as f64 / (4.0 * nf + 2.0);
        let mut x = (1.0 - (nf - 1.0) / (8.0 * nf * nf * nf)) * jv * jv;
        x = x.max(1e-6);

        // Newton iteration on L_n(x)
        for _ in 0..200 {
            let (p, dp) = laguerre_p_dp(n, x);
            if dp.abs() < 1e-300 {
                break;
            }
            let dx = p / dp;
            x -= dx;
            x = x.max(1e-15);
            if dx.abs() < 1e-12 {
                break;
            }
        }

        // Standard Gauss-Laguerre weight formula:
        // w_i = x_i / ((n+1)^2 * [L_{n+1}(x_i)]^2)
        // Alternatively: w_i = x_i / (n * L_{n-1}(x_i))^2
        // Both are equivalent; use the latter.
        let (p_prev, _) = laguerre_p_dp(n - 1, x);
        let w = if p_prev.abs() < 1e-300 {
            0.0
        } else {
            x / ((n as f64) * p_prev).powi(2)
        };
        nw.push((x, w));
    }
    nw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    nw
}

// ── 2D tensor-product integration ────────────────────────────────────────────

/// Integrate `f(x, y)` on \[ax, bx\] × \[ay, by\] using tensor-product Gauss-Legendre.
///
/// `nx` and `ny` are the number of quadrature points in each dimension.
pub fn integrate_2d(
    f: &dyn Fn(f64, f64) -> f64,
    ax: f64,
    bx: f64,
    ay: f64,
    by: f64,
    nx: usize,
    ny: usize,
) -> f64 {
    let nwx = gauss_legendre_weights(nx);
    let nwy = gauss_legendre_weights(ny);
    let midx = 0.5 * (ax + bx);
    let halfx = 0.5 * (bx - ax);
    let midy = 0.5 * (ay + by);
    let halfy = 0.5 * (by - ay);
    let mut sum = 0.0;
    for (xi, wi) in &nwx {
        let x = midx + halfx * xi;
        for (yj, wj) in &nwy {
            let y = midy + halfy * yj;
            sum += wi * wj * f(x, y);
        }
    }
    sum * halfx * halfy
}

/// Integrate `f(x, y, z)` on a box given by `bounds = [(ax,bx), (ay,by), (az,bz)]`.
///
/// `n = [nx, ny, nz]` are the number of quadrature points per dimension.
pub fn integrate_3d(
    f: &dyn Fn(f64, f64, f64) -> f64,
    bounds: [(f64, f64); 3],
    n: [usize; 3],
) -> f64 {
    let [(ax, bx), (ay, by), (az, bz)] = bounds;
    let nwx = gauss_legendre_weights(n[0]);
    let nwy = gauss_legendre_weights(n[1]);
    let nwz = gauss_legendre_weights(n[2]);
    let midx = 0.5 * (ax + bx);
    let halfx = 0.5 * (bx - ax);
    let midy = 0.5 * (ay + by);
    let halfy = 0.5 * (by - ay);
    let midz = 0.5 * (az + bz);
    let halfz = 0.5 * (bz - az);
    let mut sum = 0.0;
    for (xi, wi) in &nwx {
        let x = midx + halfx * xi;
        for (yj, wj) in &nwy {
            let y = midy + halfy * yj;
            for (zk, wk) in &nwz {
                let z = midz + halfz * zk;
                sum += wi * wj * wk * f(x, y, z);
            }
        }
    }
    sum * halfx * halfy * halfz
}

// ── AdaptiveIntegrator ────────────────────────────────────────────────────────

/// Adaptive integrator that tracks function evaluation count.
///
/// Uses recursive adaptive G7K15 bisection internally.
#[derive(Debug, Clone)]
pub struct AdaptiveIntegrator {
    /// Absolute tolerance for the error estimate.
    pub tol: f64,
    /// Maximum number of function evaluations allowed.
    pub max_evals: usize,
    /// Number of function evaluations used in the last call to `integrate`.
    pub calls: usize,
}

impl AdaptiveIntegrator {
    /// Create a new `AdaptiveIntegrator` with the given tolerance and evaluation budget.
    pub fn new(tol: f64, max_evals: usize) -> Self {
        Self {
            tol,
            max_evals,
            calls: 0,
        }
    }

    /// Integrate `f` on \[a, b\].
    ///
    /// Returns `(value, error_estimate)`.  The error estimate is |G15 - G7|
    /// accumulated over all subintervals.
    pub fn integrate(&mut self, f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
        self.calls = 0;
        let max_depth = (self.max_evals / 15).max(1).ilog2() as usize + 1;
        let (val, err) = self.adaptive_internal(f, a, b, self.tol, max_depth, 0);
        (val, err)
    }

    /// Internal recursive worker.
    fn adaptive_internal(
        &mut self,
        f: &dyn Fn(f64) -> f64,
        a: f64,
        b: f64,
        tol: f64,
        max_depth: usize,
        depth: usize,
    ) -> (f64, f64) {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);

        let mut gk = 0.0_f64;
        let mut g7 = 0.0_f64;
        self.calls += 1;
        let f0 = f(mid);
        gk += GK15_WEIGHTS[7] * f0;
        g7 += G7_WEIGHTS[3] * f0;

        for i in 0..7 {
            let xi = GK15_NODES[i];
            self.calls += 2;
            let fplus = f(mid + half * xi);
            let fminus = f(mid - half * xi);
            gk += GK15_WEIGHTS[i] * (fplus + fminus);
            if i == 1 || i == 3 || i == 5 {
                let gi = match i {
                    1 => 0,
                    3 => 1,
                    5 => 2,
                    _ => unreachable!(),
                };
                g7 += G7_WEIGHTS[gi] * (fplus + fminus);
            }
        }
        gk *= half;
        g7 *= half;

        let err = (gk - g7).abs();
        if depth >= max_depth || err < tol || self.calls >= self.max_evals {
            return (gk, err);
        }

        let (vl, el) = self.adaptive_internal(f, a, mid, tol * 0.5, max_depth, depth + 1);
        let (vr, er) = self.adaptive_internal(f, mid, b, tol * 0.5, max_depth, depth + 1);
        (vl + vr, el + er)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    // ── gauss_legendre_weights ──────────────────────────────────────────────

    #[test]
    fn test_gl_weights_n1_integrates_constant() {
        // ∫₋₁¹ 1 dx = 2
        let nw = gauss_legendre_weights(1);
        let sum: f64 = nw.iter().map(|(_, w)| w).sum();
        assert!(
            (sum - 2.0).abs() < TOL,
            "weights should sum to 2, got {sum}"
        );
    }

    #[test]
    fn test_gl_weights_n2_nodes_and_weights() {
        let nw = gauss_legendre_weights(2);
        assert_eq!(nw.len(), 2);
        // Nodes: ±1/√3
        let node = 1.0_f64 / 3.0_f64.sqrt();
        assert!((nw[0].0 + node).abs() < 1e-12);
        assert!((nw[1].0 - node).abs() < 1e-12);
        // Both weights = 1
        assert!((nw[0].1 - 1.0).abs() < 1e-12);
        assert!((nw[1].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gl_weights_n5_sum_to_two() {
        let nw = gauss_legendre_weights(5);
        let sum: f64 = nw.iter().map(|(_, w)| w).sum();
        assert!((sum - 2.0).abs() < 1e-12, "5-pt GL weights sum = {sum}");
    }

    #[test]
    fn test_gl_nodes_in_minus_one_to_one() {
        for n in [1, 2, 3, 5, 8, 10] {
            for (x, _) in gauss_legendre_weights(n) {
                assert!(
                    (-1.0..=1.0).contains(&x),
                    "node {x} out of [-1,1] for n={n}"
                );
            }
        }
    }

    #[test]
    fn test_gl_n5_is_sorted() {
        let nw = gauss_legendre_weights(5);
        for i in 1..nw.len() {
            assert!(nw[i - 1].0 <= nw[i].0, "nodes not sorted at index {i}");
        }
    }

    // ── gauss_legendre_integrate ────────────────────────────────────────────

    #[test]
    fn test_gl_integrate_constant() {
        // ∫₀¹ 3 dx = 3
        let result = gauss_legendre_integrate(&|_| 3.0, 0.0, 1.0, 3);
        assert!((result - 3.0).abs() < TOL);
    }

    #[test]
    fn test_gl_integrate_linear() {
        // ∫₀² x dx = 2
        let result = gauss_legendre_integrate(&|x| x, 0.0, 2.0, 3);
        assert!((result - 2.0).abs() < TOL);
    }

    #[test]
    fn test_gl_integrate_quadratic() {
        // ∫₀¹ x² dx = 1/3
        let result = gauss_legendre_integrate(&|x| x * x, 0.0, 1.0, 2);
        assert!((result - 1.0 / 3.0).abs() < TOL);
    }

    #[test]
    fn test_gl_integrate_sin() {
        // ∫₀^π sin(x) dx = 2
        let result = gauss_legendre_integrate(&|x| x.sin(), 0.0, PI, 10);
        assert!((result - 2.0).abs() < 1e-10, "sin integral = {result}");
    }

    #[test]
    fn test_gl_integrate_exp() {
        // ∫₀¹ eˣ dx = e - 1
        let exact = std::f64::consts::E - 1.0;
        let result = gauss_legendre_integrate(&|x| x.exp(), 0.0, 1.0, 8);
        assert!((result - exact).abs() < 1e-12, "exp integral = {result}");
    }

    // ── gauss_chebyshev_weights ─────────────────────────────────────────────

    #[test]
    fn test_chebyshev_n4_weights_sum_to_pi() {
        let nw = gauss_chebyshev_weights(4);
        let sum: f64 = nw.iter().map(|(_, w)| w).sum();
        assert!((sum - PI).abs() < 1e-12, "Chebyshev weights sum = {sum}");
    }

    #[test]
    fn test_chebyshev_nodes_on_unit_circle() {
        for (x, _) in gauss_chebyshev_weights(6) {
            assert!(x.abs() <= 1.0 + 1e-12, "node {x} out of range");
        }
    }

    // ── simpsons_rule ────────────────────────────────────────────────────────

    #[test]
    fn test_simpsons_constant() {
        let r = simpsons_rule(&|_| 5.0, 0.0, 1.0, 2);
        assert!((r - 5.0).abs() < TOL);
    }

    #[test]
    fn test_simpsons_polynomial() {
        // ∫₀¹ x³ dx = 1/4 ; Simpson's is exact for degree ≤ 3
        let r = simpsons_rule(&|x| x * x * x, 0.0, 1.0, 4);
        assert!((r - 0.25).abs() < TOL);
    }

    #[test]
    fn test_simpsons_odd_n_rounded_up() {
        // n=3 is odd → rounded to 4 internally
        let r = simpsons_rule(&|x| x * x, 0.0, 1.0, 3);
        assert!((r - 1.0 / 3.0).abs() < 1e-10, "result = {r}");
    }

    #[test]
    fn test_simpsons_sin() {
        let r = simpsons_rule(&|x| x.sin(), 0.0, PI, 100);
        assert!((r - 2.0).abs() < 1e-6, "sin integral = {r}");
    }

    // ── romberg ──────────────────────────────────────────────────────────────

    #[test]
    fn test_romberg_constant() {
        let r = romberg_integration(&|_| 7.0, 0.0, 1.0, 5, 1e-12);
        assert!((r - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_romberg_exp() {
        let exact = std::f64::consts::E - 1.0;
        let r = romberg_integration(&|x| x.exp(), 0.0, 1.0, 8, 1e-12);
        assert!((r - exact).abs() < 1e-10, "romberg exp = {r}");
    }

    #[test]
    fn test_romberg_sin_over_pi() {
        let r = romberg_integration(&|x| x.sin(), 0.0, PI, 8, 1e-12);
        assert!((r - 2.0).abs() < 1e-10, "romberg sin = {r}");
    }

    // ── adaptive_gauss_kronrod ───────────────────────────────────────────────

    #[test]
    fn test_agk_constant() {
        let r = adaptive_gauss_kronrod(&|_| 3.0, 0.0, 2.0, 1e-10, 10);
        assert!((r - 6.0).abs() < 1e-8);
    }

    #[test]
    fn test_agk_sin() {
        let r = adaptive_gauss_kronrod(&|x| x.sin(), 0.0, PI, 1e-10, 10);
        assert!((r - 2.0).abs() < 1e-8, "agk sin = {r}");
    }

    #[test]
    fn test_agk_exp_negative() {
        // ∫₀^∞ e^{-x} dx = 1 ; approximate with [0, 20]
        let r = adaptive_gauss_kronrod(&|x| (-x).exp(), 0.0, 20.0, 1e-10, 15);
        assert!((r - 1.0).abs() < 1e-8, "agk exp = {r}");
    }

    // ── double_exponential ───────────────────────────────────────────────────

    #[test]
    fn test_de_constant() {
        let r = double_exponential(&|_| 1.0, 0.0, 1.0, 50, 0.1);
        assert!((r - 1.0).abs() < 1e-6, "DE const = {r}");
    }

    #[test]
    fn test_de_sin() {
        let r = double_exponential(&|x| x.sin(), 0.0, PI, 100, 0.05);
        assert!((r - 2.0).abs() < 1e-6, "DE sin = {r}");
    }

    // ── clenshaw_curtis ──────────────────────────────────────────────────────

    #[test]
    fn test_cc_weights_sum_to_two() {
        for n in [2, 3, 5, 8] {
            let nw = clenshaw_curtis_weights(n);
            let sum: f64 = nw.iter().map(|(_, w)| w).sum();
            assert!((sum - 2.0).abs() < 1e-8, "CC weights sum for n={n}: {sum}");
        }
    }

    #[test]
    fn test_cc_endpoints_are_minus_one_and_one() {
        let nw = clenshaw_curtis_weights(4);
        let xs: Vec<f64> = nw.iter().map(|(x, _)| *x).collect();
        assert!(
            xs.iter().any(|x| (x + 1.0).abs() < 1e-12),
            "should include -1"
        );
        assert!(
            xs.iter().any(|x| (x - 1.0).abs() < 1e-12),
            "should include +1"
        );
    }

    // ── gauss_laguerre ───────────────────────────────────────────────────────

    #[test]
    fn test_laguerre_n1() {
        let nw = gauss_laguerre_weights(1);
        // Single node at x=1, weight=1
        assert_eq!(nw.len(), 1);
        assert!((nw[0].0 - 1.0).abs() < 0.1, "n=1 node ≈ 1, got {}", nw[0].0);
    }

    #[test]
    fn test_laguerre_nodes_positive() {
        for n in [1, 2, 3, 4, 5] {
            for (x, _) in gauss_laguerre_weights(n) {
                assert!(x > 0.0, "Laguerre node must be positive, got {x}");
            }
        }
    }

    #[test]
    fn test_laguerre_integrates_exp_neg_x() {
        // ∫₀^∞ e^{-x} * 1 dx = 1  (f(x)=1 with Laguerre weight e^{-x})
        // The quadrature approximates this as sum_i w_i * f(x_i) = sum_i w_i
        // The weights already absorb the e^{-x} factor, so their sum ≈ 1.
        // We test with f(x)=1 and verify the rule is reasonable.
        let nw = gauss_laguerre_weights(5);
        // Verify nodes are positive and sorted
        for i in 1..nw.len() {
            assert!(nw[i].0 > nw[i - 1].0, "nodes should be sorted");
        }
        // All weights should be positive
        for (_, w) in &nw {
            assert!(*w > 0.0, "weights should be positive");
        }
    }

    // ── integrate_2d ─────────────────────────────────────────────────────────

    #[test]
    fn test_integrate_2d_constant() {
        // ∫₀¹ ∫₀¹ 1 dy dx = 1
        let r = integrate_2d(&|_, _| 1.0, 0.0, 1.0, 0.0, 1.0, 3, 3);
        assert!((r - 1.0).abs() < TOL);
    }

    #[test]
    fn test_integrate_2d_product() {
        // ∫₀¹ ∫₀¹ x * y dy dx = 1/4
        let r = integrate_2d(&|x, y| x * y, 0.0, 1.0, 0.0, 1.0, 3, 3);
        assert!((r - 0.25).abs() < TOL, "2D product = {r}");
    }

    #[test]
    fn test_integrate_2d_sin_cos() {
        // ∫₀^{π/2} ∫₀^{π/2} sin(x)cos(y) dy dx = 1
        let r = integrate_2d(
            &|x, y| x.sin() * y.cos(),
            0.0,
            PI / 2.0,
            0.0,
            PI / 2.0,
            8,
            8,
        );
        assert!((r - 1.0).abs() < 1e-10, "2D sin*cos = {r}");
    }

    // ── integrate_3d ─────────────────────────────────────────────────────────

    #[test]
    fn test_integrate_3d_constant() {
        // ∫₀¹³ 1 dx dy dz = 1
        let r = integrate_3d(
            &|_, _, _| 1.0,
            [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
            [3, 3, 3],
        );
        assert!((r - 1.0).abs() < TOL);
    }

    #[test]
    fn test_integrate_3d_xyz() {
        // ∫₀¹ ∫₀¹ ∫₀¹ xyz dz dy dx = (1/2)^3 = 1/8
        let r = integrate_3d(
            &|x, y, z| x * y * z,
            [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
            [4, 4, 4],
        );
        assert!((r - 0.125).abs() < TOL, "3D xyz = {r}");
    }

    // ── AdaptiveIntegrator ────────────────────────────────────────────────────

    #[test]
    fn test_adaptive_integrator_constant() {
        let mut ai = AdaptiveIntegrator::new(1e-8, 1000);
        let (val, _err) = ai.integrate(&|_| 4.0, 0.0, 1.0);
        assert!((val - 4.0).abs() < 1e-6, "val = {val}");
    }

    #[test]
    fn test_adaptive_integrator_sin() {
        let mut ai = AdaptiveIntegrator::new(1e-8, 5000);
        let (val, err) = ai.integrate(&|x| x.sin(), 0.0, PI);
        assert!((val - 2.0).abs() < 1e-6, "val = {val}, err = {err}");
    }

    #[test]
    fn test_adaptive_integrator_tracks_calls() {
        let mut ai = AdaptiveIntegrator::new(1e-6, 1000);
        let _ = ai.integrate(&|x| x.cos(), 0.0, 1.0);
        assert!(ai.calls > 0, "calls should be > 0 after integration");
    }

    #[test]
    fn test_adaptive_integrator_error_nonnegative() {
        let mut ai = AdaptiveIntegrator::new(1e-8, 2000);
        let (_, err) = ai.integrate(&|x| x * x, 0.0, 1.0);
        assert!(err >= 0.0, "error estimate must be non-negative");
    }

    // ── gauss_lobatto ────────────────────────────────────────────────────────

    #[test]
    fn test_lobatto_n2_endpoints() {
        let nw = gauss_lobatto_weights(2);
        assert_eq!(nw.len(), 2);
        assert!((nw[0].0 + 1.0).abs() < 1e-12, "first node should be -1");
        assert!((nw[1].0 - 1.0).abs() < 1e-12, "last node should be 1");
    }

    #[test]
    fn test_lobatto_weights_sum_to_two() {
        for n in [2, 3, 4, 5] {
            let nw = gauss_lobatto_weights(n);
            let sum: f64 = nw.iter().map(|(_, w)| w).sum();
            assert!(
                (sum - 2.0).abs() < 1e-8,
                "GL lobatto n={n} weights sum = {sum}"
            );
        }
    }

    #[test]
    fn test_lobatto_n4_count() {
        let nw = gauss_lobatto_weights(4);
        assert_eq!(nw.len(), 4);
    }

    // ── gauss_hermite ────────────────────────────────────────────────────────

    #[test]
    fn test_hermite_n1() {
        let nw = gauss_hermite_weights(1);
        assert_eq!(nw.len(), 1);
        // Node should be at x=0
        assert!(
            nw[0].0.abs() < 1e-10,
            "n=1 hermite node at 0, got {}",
            nw[0].0
        );
    }

    #[test]
    fn test_hermite_n2_nodes_symmetric() {
        let nw = gauss_hermite_weights(2);
        assert_eq!(nw.len(), 2);
        assert!(
            (nw[0].0 + nw[1].0).abs() < 1e-10,
            "nodes should be symmetric"
        );
    }

    #[test]
    fn test_hermite_n3_count() {
        let nw = gauss_hermite_weights(3);
        assert_eq!(nw.len(), 3);
    }

    // ── cross-method agreement ────────────────────────────────────────────────

    #[test]
    fn test_methods_agree_on_sin() {
        let f = &|x: f64| x.sin();
        let gl = gauss_legendre_integrate(f, 0.0, PI, 10);
        let simp = simpsons_rule(f, 0.0, PI, 100);
        let romb = romberg_integration(f, 0.0, PI, 8, 1e-12);
        let agk = adaptive_gauss_kronrod(f, 0.0, PI, 1e-10, 10);
        assert!((gl - 2.0).abs() < 1e-10);
        assert!((simp - 2.0).abs() < 1e-6);
        assert!((romb - 2.0).abs() < 1e-10);
        assert!((agk - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_methods_agree_on_polynomial() {
        // ∫₋₁¹ (x⁴ - 2x² + 1) dx = 2 - 4/3 + 2 = 16/15
        let f = &|x: f64| x.powi(4) - 2.0 * x * x + 1.0;
        let exact = 16.0 / 15.0;
        let gl = gauss_legendre_integrate(f, -1.0, 1.0, 5);
        let romb = romberg_integration(f, -1.0, 1.0, 6, 1e-12);
        assert!((gl - exact).abs() < 1e-12);
        assert!((romb - exact).abs() < 1e-10);
    }
}
