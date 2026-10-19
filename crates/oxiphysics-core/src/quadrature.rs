// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical quadrature (numerical integration) methods.
//!
//! Provides Gauss-Legendre, Gauss-Lobatto, Gauss-Chebyshev, Gauss-Hermite,
//! Gauss-Laguerre, Clenshaw-Curtis, Simpson's rule, Romberg, adaptive
//! Gauss-Kronrod (G7K15), double-exponential (tanh-sinh), and 2D/3D tensor
//! product quadrature.

use nalgebra::DMatrix;
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

/// Magnitude above which the Hermite recurrence is rescaled to avoid overflow.
const HERMITE_RESCALE_LIMIT: f64 = 1.0e150;
/// Factor applied to the Hermite recurrence when it exceeds [`HERMITE_RESCALE_LIMIT`].
const HERMITE_RESCALE_FACTOR: f64 = 1.0e-150;

/// Three-term recurrence of the orthonormal Hermite polynomials (weight `exp(-x²)`),
/// scaled by `π^(1/4)`.
///
/// `q_k = π^(1/4) p_k` obeys `q_0 = 1` and
/// `q_(k+1) = sqrt(2/(k+1)) x q_k - sqrt(k/(k+1)) q_(k-1)`, where `p_k` is the
/// polynomial of degree `k` with unit norm for the weight `exp(-x²)`; `p_n`
/// vanishes exactly at the Gauss-Hermite nodes.
///
/// Returns `(q_n, q_(n-1), r)` for `n >= 1`.  Whenever the running pair would
/// exceed [`HERMITE_RESCALE_LIMIT`] it is multiplied by [`HERMITE_RESCALE_FACTOR`]
/// and `r` is incremented, so that the true values are the returned ones times
/// `1e150^r`.  This keeps the recurrence finite for any `n` and `x`.
fn hermite_recurrence(n: usize, x: f64) -> (f64, f64, i32) {
    let mut q_prev = 0.0_f64;
    let mut q = 1.0_f64;
    let mut rescales = 0_i32;
    for k in 0..n {
        let kf = k as f64;
        let q_next = (2.0 / (kf + 1.0)).sqrt() * x * q - (kf / (kf + 1.0)).sqrt() * q_prev;
        q_prev = q;
        q = q_next;
        if q.abs() > HERMITE_RESCALE_LIMIT {
            q *= HERMITE_RESCALE_FACTOR;
            q_prev *= HERMITE_RESCALE_FACTOR;
            rescales += 1;
        }
    }
    (q, q_prev, rescales)
}

/// Refine an approximate root `x_start` of the physicist's Hermite polynomial
/// `H_n` with Newton's method on the orthonormal recurrence.
///
/// The derivative follows from `p_n' = sqrt(2n) p_(n-1)`.  The starting value
/// comes from an eigenvalue solver and is already accurate to a few ulp, so the
/// iteration only removes the last rounding errors; if it ever moves the root
/// by more than `1e-8 (1 + |x_start|)` the starting value is kept.
fn polish_hermite_node(n: usize, x_start: f64) -> f64 {
    let sqrt_two_n = (2.0 * n as f64).sqrt();
    let mut x = x_start;
    for _ in 0..8 {
        let (q, q_prev, _) = hermite_recurrence(n, x);
        let step = q / (sqrt_two_n * q_prev);
        if !step.is_finite() {
            break;
        }
        x -= step;
        if step.abs() <= 4.0 * f64::EPSILON * x.abs() {
            break;
        }
    }
    if (x - x_start).abs() <= 1.0e-8 * (1.0 + x_start.abs()) {
        x
    } else {
        x_start
    }
}

/// Gauss-Hermite weight of the node `x`, a root of `H_n`:
/// `w = 2^(n-1) n! sqrt(π) / (n² H_(n-1)(x)²) = sqrt(π) / (n q_(n-1)(x)²)`.
///
/// Evaluated with the scaled orthonormal recurrence, so it keeps full relative
/// accuracy even for the tiny weights of the outermost nodes.
fn hermite_weight(n: usize, x: f64) -> f64 {
    let (_, q_prev, rescales) = hermite_recurrence(n, x);
    let inverse = q_prev.recip();
    PI.sqrt() / n as f64 * inverse * inverse * HERMITE_RESCALE_FACTOR.powi(2 * rescales)
}

/// Compute the *n*-point Gauss-Hermite nodes and weights.
///
/// The rule integrates `f(x) · exp(-x²)` over (-∞, +∞) and is exact for
/// polynomials `f` of degree up to `2n - 1`.  The nodes are the roots of the
/// physicist's Hermite polynomial `H_n(x)` and the weights are
/// `w_i = 2^(n-1) n! √π / (n² H_(n-1)(x_i)²)`, so that `Σ w_i = √π`.
///
/// # Algorithm
///
/// The nodes start as the eigenvalues of the symmetric tridiagonal Jacobi
/// matrix of the Hermite polynomials (zero diagonal, off-diagonal entries
/// `sqrt(k/2)` for `k = 1..n`; Golub-Welsch), computed with a symmetric
/// eigen-solver.  Unlike Newton iteration from asymptotic guesses, this cannot
/// collapse several starting points onto the same root.  Each node is then
/// polished with Newton's method and the weight follows from the orthonormal
/// recurrence, which gives machine-precision nodes and weights with full
/// relative accuracy in the tails.  The rule is symmetric: the computation is
/// done for the non-negative nodes and mirrored, so `x_(n-1-i) = -x_i` and
/// `w_(n-1-i) = w_i` hold exactly, and the middle node of an odd rule is
/// exactly `0`.
///
/// Returns `(node, weight)` pairs in ascending node order, or an empty vector
/// for `n = 0`.  The cost is dominated by the dense eigen-solve (`O(n³)`).
pub fn gauss_hermite_weights(n: usize) -> Vec<(f64, f64)> {
    if n == 0 {
        return Vec::new();
    }
    let jacobi = DMatrix::<f64>::from_fn(n, n, |i, j| {
        if i + 1 == j {
            (j as f64 / 2.0).sqrt()
        } else if j + 1 == i {
            (i as f64 / 2.0).sqrt()
        } else {
            0.0
        }
    });
    let mut eigenvalues: Vec<f64> = jacobi.symmetric_eigenvalues().iter().copied().collect();
    eigenvalues.sort_by(|a, b| a.total_cmp(b));
    let mut nw = Vec::with_capacity(n);
    for i in 0..n / 2 {
        // The spectrum is symmetric about 0: pair the i-th smallest eigenvalue
        // with the i-th largest one and keep the non-negative representative.
        let magnitude = 0.5 * (eigenvalues[n - 1 - i] - eigenvalues[i]);
        let x = polish_hermite_node(n, magnitude);
        let w = hermite_weight(n, x);
        nw.push((-x, w));
        nw.push((x, w));
    }
    if n % 2 == 1 {
        nw.push((0.0, hermite_weight(n, 0.0)));
    }
    nw.sort_by(|a, b| a.0.total_cmp(&b.0));
    nw
}

// ── Gauss-Laguerre ────────────────────────────────────────────────────────────

/// Magnitude above which the Laguerre recurrence is rescaled to avoid overflow.
const LAGUERRE_RESCALE_LIMIT: f64 = 1.0e150;
/// Factor applied to the Laguerre recurrence when it exceeds [`LAGUERRE_RESCALE_LIMIT`].
const LAGUERRE_RESCALE_FACTOR: f64 = 1.0e-150;

/// Difference form of the Laguerre three-term recurrence.
///
/// With `d_k = L_k - L_(k-1)` the recurrence
/// `(k + 1) L_(k+1) = (2k + 1 - x) L_k - k L_(k-1)` is equivalent to
/// `d_(k+1) = (k d_k - x L_k) / (k + 1)` and `L_(k+1) = L_k + d_(k+1)`, started
/// from `L_1 = 1 - x`, `d_1 = -x`.  Carrying the differences avoids the
/// cancellation of the plain recurrence near `x = 0`, so the zeros of `L_n` are
/// resolved to a few ulp, and it yields the derivative without a further
/// subtraction: `L_n'(x) = n (L_n - L_(n-1)) / x = n d_n / x`.
///
/// Returns `(L_n, d_n, r)` for `n >= 1`.  Whenever the running pair would exceed
/// [`LAGUERRE_RESCALE_LIMIT`] it is multiplied by [`LAGUERRE_RESCALE_FACTOR`] and `r`
/// is incremented, so that the true values are the returned ones times `1e150^r`.
/// This keeps the recurrence finite for any `n` and `x` (`|L_n(x)|` grows like
/// `exp(x / 2)`).
fn laguerre_recurrence(n: usize, x: f64) -> (f64, f64, i32) {
    let mut value = 1.0 - x;
    let mut difference = -x;
    let mut rescales = 0_i32;
    for k in 1..n {
        let kf = k as f64;
        difference = (kf * difference - x * value) / (kf + 1.0);
        value += difference;
        if value.abs() > LAGUERRE_RESCALE_LIMIT || difference.abs() > LAGUERRE_RESCALE_LIMIT {
            value *= LAGUERRE_RESCALE_FACTOR;
            difference *= LAGUERRE_RESCALE_FACTOR;
            rescales += 1;
        }
    }
    (value, difference, rescales)
}

/// Refine an approximate zero `x_start` of `L_n` with Newton's method on the
/// difference recurrence.
///
/// Since `L_n' = n d_n / x` the Newton step is `x L_n / (n d_n)`.  The starting
/// value comes from an eigenvalue solver and is already accurate to about
/// `1e-13 (1 + x)`, so the iteration only removes the last rounding errors; if it
/// ever moves the root by more than `1e-8 (1 + |x_start|)` the starting value is
/// kept.
fn polish_laguerre_node(n: usize, x_start: f64) -> f64 {
    let nf = n as f64;
    let mut x = x_start;
    for _ in 0..8 {
        let (value, difference, _) = laguerre_recurrence(n, x);
        let step = x * value / (nf * difference);
        if !step.is_finite() {
            break;
        }
        x -= step;
        if step.abs() <= 4.0 * f64::EPSILON * x.abs() {
            break;
        }
    }
    if (x - x_start).abs() <= 1.0e-8 * (1.0 + x_start.abs()) {
        x
    } else {
        x_start
    }
}

/// Gauss-Laguerre weight of the node `x`, a zero of `L_n`:
/// `w = x / (n L_(n-1)(x))²`.
///
/// `L_(n-1)(x) = L_n(x) - d_n(x)` comes out of the rescaled difference
/// recurrence, so the weight keeps full relative accuracy for the tiny weights of
/// the outermost nodes and never overflows.  A weight below the smallest positive
/// `f64` underflows to `0`.
fn laguerre_weight(n: usize, x: f64) -> f64 {
    let (value, difference, rescales) = laguerre_recurrence(n, x);
    let inverse = (value - difference).recip();
    let nf = n as f64;
    x / (nf * nf) * inverse * inverse * LAGUERRE_RESCALE_FACTOR.powi(2 * rescales)
}

/// Compute the *n*-point Gauss-Laguerre nodes and weights.
///
/// The rule integrates `f(x) · exp(-x)` over \[0, +∞) and is exact for
/// polynomials `f` of degree up to `2n - 1`.  The nodes are the zeros of the
/// Laguerre polynomial `L_n(x)` and the weights are
/// `w_i = x_i / (n L_(n-1)(x_i))²`, so that `Σ w_i = 1` and, more generally,
/// `Σ w_i x_i^k = k!` for `k <= 2n - 1`.
///
/// # Algorithm
///
/// The nodes start as the eigenvalues of the symmetric tridiagonal Jacobi matrix
/// of the Laguerre polynomials (diagonal `2k + 1`, off-diagonal entries `k` for
/// `k = 1..n`; Golub-Welsch), computed with a symmetric eigen-solver.  Unlike
/// Newton iteration from asymptotic starting values, this cannot collapse several
/// starting points onto the same zero.  Each node is then polished with Newton's
/// method on the difference form of the Laguerre recurrence (where
/// `L_n' = n (L_n - L_(n-1)) / x`) and the weight follows from the same
/// recurrence, which gives full relative accuracy for the nodes and the tiny
/// outer weights.
///
/// Returns `(node, weight)` pairs in ascending node order, or an empty vector for
/// `n = 0`.  The cost is dominated by the dense eigen-solve (`O(n³)`).  The
/// weight of the outermost node is roughly `exp(-4n)`; once it drops below the
/// smallest positive `f64` (for `n` above about 180) it is returned as `0`.
pub fn gauss_laguerre_weights(n: usize) -> Vec<(f64, f64)> {
    if n == 0 {
        return Vec::new();
    }
    let jacobi = DMatrix::<f64>::from_fn(n, n, |i, j| {
        if i == j {
            (2 * i + 1) as f64
        } else if i + 1 == j {
            j as f64
        } else if j + 1 == i {
            i as f64
        } else {
            0.0
        }
    });
    let mut eigenvalues: Vec<f64> = jacobi.symmetric_eigenvalues().iter().copied().collect();
    eigenvalues.sort_by(|a, b| a.total_cmp(b));
    eigenvalues
        .into_iter()
        .map(|start| {
            let x = polish_laguerre_node(n, start);
            (x, laguerre_weight(n, x))
        })
        .collect()
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

    /// `∫₀^∞ x^k exp(-x) dx = k!`.
    fn laguerre_moment(k: usize) -> f64 {
        (1..=k).fold(1.0_f64, |acc, j| acc * j as f64)
    }

    /// `L_n(x)` and `L_n'(x)` from the difference recurrence for arguments that
    /// need no rescaling.
    fn laguerre_value_and_derivative(n: usize, x: f64) -> (f64, f64) {
        let (value, difference, rescales) = laguerre_recurrence(n, x);
        assert_eq!(rescales, 0, "n={n}, x={x} must not trigger rescaling");
        (value, n as f64 * difference / x)
    }

    fn assert_laguerre_rule(nodes: &[f64], weights: &[f64], node_rtol: f64, weight_rtol: f64) {
        let n = nodes.len();
        let nw = gauss_laguerre_weights(n);
        assert_eq!(nw.len(), n);
        for i in 0..n {
            assert_rel_close(nw[i].0, nodes[i], node_rtol, &format!("n={n} node {i}"));
            assert_rel_close(
                nw[i].1,
                weights[i],
                weight_rtol,
                &format!("n={n} weight {i}"),
            );
        }
    }

    #[test]
    fn test_laguerre_recurrence_matches_closed_forms() {
        // L_2 = (x² - 4x + 2)/2 with L_2' = x - 2 (the old code used -L_1 = x - 1);
        // L_3 = (-x³ + 9x² - 18x + 6)/6 with L_3' = (-x² + 6x - 6)/2.
        for &x in &[0.05_f64, 0.5, 1.0, 2.5, 7.0, 13.0] {
            let scale = 1.0 + x * x * x;
            let (p2, dp2) = laguerre_value_and_derivative(2, x);
            assert!(
                (p2 - (x * x - 4.0 * x + 2.0) / 2.0).abs() <= 1e-14 * scale,
                "L_2({x})"
            );
            assert!((dp2 - (x - 2.0)).abs() <= 1e-14 * scale, "L_2'({x})");
            let (p3, dp3) = laguerre_value_and_derivative(3, x);
            let l3 = (-x * x * x + 9.0 * x * x - 18.0 * x + 6.0) / 6.0;
            assert!((p3 - l3).abs() <= 1e-14 * scale, "L_3({x})");
            let dl3 = (-x * x + 6.0 * x - 6.0) / 2.0;
            assert!((dp3 - dl3).abs() <= 1e-14 * scale, "L_3'({x})");
        }
    }

    #[test]
    fn test_laguerre_derivative_matches_central_differences() {
        for n in [1_usize, 2, 5, 12, 30] {
            for &x in &[0.3_f64, 1.7, 6.0, 20.0] {
                let h = 1.0e-5 * x;
                let upper = laguerre_value_and_derivative(n, x + h).0;
                let lower = laguerre_value_and_derivative(n, x - h).0;
                let numeric = (upper - lower) / (2.0 * h);
                let (_, analytic) = laguerre_value_and_derivative(n, x);
                assert!(
                    (numeric - analytic).abs() <= 1e-6 * (1.0 + analytic.abs()),
                    "n={n}, x={x}: finite difference {numeric}, analytic {analytic}"
                );
            }
        }
    }

    #[test]
    fn test_laguerre_n0_is_empty() {
        assert!(gauss_laguerre_weights(0).is_empty());
    }

    #[test]
    fn test_laguerre_n1_is_one_at_one() {
        let nw = gauss_laguerre_weights(1);
        assert_eq!(nw.len(), 1);
        assert_rel_close(nw[0].0, 1.0, 1e-15, "x(n=1)");
        assert_rel_close(nw[0].1, 1.0, 1e-15, "w(n=1)");
    }

    #[test]
    fn test_laguerre_n2_matches_closed_form() {
        // Nodes 2 ∓ √2 with weights (2 ± √2) / 4.
        let nw = gauss_laguerre_weights(2);
        let root_two = 2.0_f64.sqrt();
        assert_eq!(nw.len(), 2);
        assert_rel_close(nw[0].0, 2.0 - root_two, 1e-15, "x0");
        assert_rel_close(nw[1].0, 2.0 + root_two, 1e-15, "x1");
        assert_rel_close(nw[0].1, (2.0 + root_two) / 4.0, 1e-14, "w0");
        assert_rel_close(nw[1].1, (2.0 - root_two) / 4.0, 1e-14, "w1");
    }

    #[test]
    fn test_laguerre_matches_high_precision_references() {
        // 70-digit values (mpmath), rounded to f64; they agree with
        // scipy.special.roots_laguerre to its own accuracy.
        assert_laguerre_rule(
            &[0.4157745567834791, 2.294280360279042, 6.2899450829374794],
            &[0.711093009929173, 0.27851773356924087, 0.010389256501586135],
            1e-14,
            1e-13,
        );
        assert_laguerre_rule(
            &[
                0.2635603197181409,
                1.4134030591065168,
                3.596425771040722,
                7.085810005858837,
                12.640800844275782,
            ],
            &[
                0.5217556105828086,
                0.3986668110831759,
                0.0759424496817076,
                0.0036117586799220484,
                2.3369972385776228e-5,
            ],
            1e-14,
            1e-13,
        );
        assert_laguerre_rule(
            &[
                0.13779347054049243,
                0.7294545495031705,
                1.808342901740316,
                3.4014336978548996,
                5.552496140063804,
                8.330152746764497,
                11.843785837900066,
                16.279257831378104,
                21.99658581198076,
                29.92069701227389,
            ],
            &[
                0.30844111576502015,
                0.40111992915527356,
                0.2180682876118094,
                0.062087456098677746,
                0.0095015169751811,
                0.0007530083885875388,
                2.8259233495995656e-5,
                4.2493139849626863e-7,
                1.8395648239796308e-9,
                9.911827219609008e-13,
            ],
            1e-14,
            1e-12,
        );
    }

    #[test]
    fn test_laguerre_nodes_are_distinct_sorted_and_weights_positive() {
        // Newton iteration from the old Tricomi guesses with L_n' = -L_(n-1)
        // returned [3.60, 7.09, 12.64, 12.64, 12.64] at n = 5 (Σw = 0.08).
        for n in 1..=60 {
            let nw = gauss_laguerre_weights(n);
            assert_eq!(nw.len(), n);
            for pair in nw.windows(2) {
                assert!(
                    pair[1].0 > pair[0].0,
                    "n={n}: nodes not strictly increasing ({} then {})",
                    pair[0].0,
                    pair[1].0
                );
            }
            for (x, w) in &nw {
                assert!(*x > 0.0 && *x < 4.0 * n as f64, "n={n}: node {x}");
                assert!(*w > 0.0 && w.is_finite(), "n={n}: weight {w} at x={x}");
            }
        }
    }

    #[test]
    fn test_laguerre_weights_sum_to_one_and_reproduce_the_cubic_moment() {
        for n in [2_usize, 5, 10, 20] {
            let nw = gauss_laguerre_weights(n);
            let sum: f64 = nw.iter().map(|(_, w)| w).sum();
            assert_rel_close(sum, 1.0, 1e-12, &format!("n={n} Σw"));
            let cubic: f64 = nw.iter().map(|(x, w)| w * x.powi(3)).sum();
            assert_rel_close(cubic, 6.0, 1e-12, &format!("n={n} Σw·x³"));
        }
    }

    #[test]
    fn test_laguerre_is_exact_up_to_degree_two_n_minus_one() {
        for n in 1..=30 {
            let nw = gauss_laguerre_weights(n);
            for k in 0..2 * n {
                let moment: f64 = nw.iter().map(|(x, w)| w * x.powi(k as i32)).sum();
                assert_rel_close(
                    moment,
                    laguerre_moment(k),
                    1e-11,
                    &format!("n={n} degree {k}"),
                );
            }
        }
    }

    #[test]
    fn test_laguerre_outer_weights_keep_relative_accuracy() {
        // Outermost nodes and weights of the 20-point rule (70-digit reference values).
        let nw = gauss_laguerre_weights(20);
        let nodes = [
            (0, 0.07053988969198875),
            (1, 0.37212681800161146),
            (17, 47.6199940473465),
            (18, 55.810795750063896),
            (19, 66.52441652561575),
        ];
        for (i, node) in nodes {
            assert_rel_close(nw[i].0, node, 1e-14, &format!("x{i}"));
        }
        let weights = [
            (0, 0.16874680185111388),
            (1, 0.2912543620060683),
            (17, 1.5395221405823435e-20),
            (18, 5.286442725569158e-24),
            (19, 1.6564566124990233e-28),
        ];
        for (i, weight) in weights {
            assert_rel_close(nw[i].1, weight, 1e-12, &format!("w{i}"));
        }
    }

    #[test]
    fn test_laguerre_large_rule_is_finite_and_normalised() {
        let n = 150;
        let nw = gauss_laguerre_weights(n);
        assert_eq!(nw.len(), n);
        for pair in nw.windows(2) {
            assert!(pair[1].0 > pair[0].0);
        }
        for (x, w) in &nw {
            assert!(x.is_finite() && w.is_finite() && *w > 0.0, "x={x}, w={w}");
        }
        for (k, expected) in [(0_i32, 1.0), (1, 1.0), (2, 2.0)] {
            let moment: f64 = nw.iter().map(|(x, w)| w * x.powi(k)).sum();
            assert_rel_close(moment, expected, 1e-11, &format!("Σw·x^{k}"));
        }
    }

    #[test]
    fn test_laguerre_beyond_the_f64_range_stays_finite() {
        // L_399(x) ~ exp(x / 2) exceeds f64::MAX for the outermost nodes of the
        // 400-point rule; the unscaled recurrence gives NaN there, the rescaled
        // one underflows the (true) weights below 1e-308 to zero.
        let nw = gauss_laguerre_weights(400);
        assert_eq!(nw.len(), 400);
        for pair in nw.windows(2) {
            assert!(pair[1].0 > pair[0].0);
        }
        for (x, w) in &nw {
            assert!(x.is_finite() && w.is_finite() && *w >= 0.0, "x={x}, w={w}");
        }
        assert!(nw[0].1 > 0.0);
        assert_eq!(nw[399].1, 0.0);
        let sum: f64 = nw.iter().map(|(_, w)| w).sum();
        assert_rel_close(sum, 1.0, 1e-10, "Σw");
    }

    #[test]
    fn test_laguerre_recurrence_stays_finite_for_huge_arguments() {
        let (value, difference, rescales) = laguerre_recurrence(2000, 1500.0);
        assert!(value.is_finite() && difference.is_finite());
        assert!(rescales > 0, "rescales={rescales}");
        let w = laguerre_weight(2000, 1500.0);
        assert!(w.is_finite() && w >= 0.0, "w={w}");
    }

    #[test]
    fn test_laguerre_integrates_smooth_functions() {
        let nw = gauss_laguerre_weights(20);
        // ∫₀^∞ exp(-x) exp(-x/2) dx = 2/3
        let decay: f64 = nw.iter().map(|(x, w)| w * (-0.5 * x).exp()).sum();
        assert_rel_close(decay, 2.0 / 3.0, 1e-13, "∫exp(-3x/2)");
        // ∫₀^∞ exp(-x) cos(x) dx = 1/2
        let oscillation: f64 = nw.iter().map(|(x, w)| w * x.cos()).sum();
        assert_rel_close(oscillation, 0.5, 1e-11, "∫exp(-x)cos(x)");
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

    /// `∫ x^(2k) exp(-x²) dx = Γ(k + 1/2) = √π (2k-1)!! / 2^k`.
    fn gaussian_even_moment(k: usize) -> f64 {
        let mut value = PI.sqrt();
        for j in 1..=k {
            value *= (2 * j - 1) as f64 / 2.0;
        }
        value
    }

    fn assert_rel_close(actual: f64, expected: f64, rtol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= rtol * expected.abs(),
            "{what}: got {actual:e}, expected {expected:e}"
        );
    }

    #[test]
    fn test_hermite_n0_is_empty() {
        assert!(gauss_hermite_weights(0).is_empty());
    }

    #[test]
    fn test_hermite_n1_is_sqrt_pi_at_the_origin() {
        let nw = gauss_hermite_weights(1);
        assert_eq!(nw.len(), 1);
        assert_eq!(nw[0].0, 0.0);
        assert_rel_close(nw[0].1, PI.sqrt(), 1e-15, "w(n=1)");
    }

    #[test]
    fn test_hermite_n2_matches_closed_form() {
        // Nodes ±1/√2, weights √π/2 each; the old weights summed to 2√π.
        let nw = gauss_hermite_weights(2);
        let node = 1.0 / 2.0_f64.sqrt();
        assert_rel_close(nw[0].0, -node, 1e-15, "x0");
        assert_rel_close(nw[1].0, node, 1e-15, "x1");
        assert_rel_close(nw[0].1, PI.sqrt() / 2.0, 1e-14, "w0");
        assert_rel_close(nw[1].1, PI.sqrt() / 2.0, 1e-14, "w1");
    }

    #[test]
    fn test_hermite_n3_matches_closed_form() {
        // Nodes 0, ±√(3/2); weights 2√π/3 (centre) and √π/6.
        let nw = gauss_hermite_weights(3);
        let outer = 1.5_f64.sqrt();
        assert_rel_close(nw[0].0, -outer, 1e-15, "x0");
        assert_eq!(nw[1].0, 0.0);
        assert_rel_close(nw[2].0, outer, 1e-15, "x2");
        assert_rel_close(nw[0].1, PI.sqrt() / 6.0, 1e-14, "w0");
        assert_rel_close(nw[1].1, 2.0 * PI.sqrt() / 3.0, 1e-14, "w1");
        assert_rel_close(nw[2].1, PI.sqrt() / 6.0, 1e-14, "w2");
    }

    #[test]
    fn test_hermite_n5_matches_high_precision_reference() {
        // 60-digit Golub-Welsch values (mpmath), rounded to f64.
        let nodes = [
            -2.0201828704560856,
            -0.9585724646138185,
            0.0,
            0.9585724646138185,
            2.0201828704560856,
        ];
        let weights = [
            0.019953242059045913,
            0.3936193231522412,
            0.9453087204829419,
            0.3936193231522412,
            0.019953242059045913,
        ];
        let nw = gauss_hermite_weights(5);
        assert_eq!(nw.len(), 5);
        for i in 0..5 {
            assert!((nw[i].0 - nodes[i]).abs() < 2e-15, "node {i}: {}", nw[i].0);
            assert_rel_close(nw[i].1, weights[i], 2e-15, &format!("weight {i}"));
        }
    }

    #[test]
    fn test_hermite_n10_matches_high_precision_reference() {
        let nodes = [
            3.4361591188377374,
            2.5327316742327897,
            1.7566836492998819,
            1.0366108297895138,
            0.3429013272237046,
        ];
        let weights = [
            7.64043285523262e-6,
            0.0013436457467812326,
            0.033874394455481065,
            0.24013861108231468,
            0.6108626337353258,
        ];
        let nw = gauss_hermite_weights(10);
        assert_eq!(nw.len(), 10);
        for i in 0..5 {
            let upper = &nw[9 - i];
            let lower = &nw[i];
            assert!((upper.0 - nodes[i]).abs() < 3e-15, "node {i}: {}", upper.0);
            assert!((lower.0 + nodes[i]).abs() < 3e-15, "node -{i}: {}", lower.0);
            assert_rel_close(upper.1, weights[i], 2e-14, &format!("weight {i}"));
            assert_rel_close(lower.1, weights[i], 2e-14, &format!("weight -{i}"));
        }
    }

    #[test]
    fn test_hermite_nodes_are_distinct_sorted_and_weights_positive() {
        // Newton iteration from the old asymptotic guesses produced only 3 distinct
        // nodes at n = 5, 6 at n = 10 and 10 at n = 20.
        for n in 1..=60 {
            let nw = gauss_hermite_weights(n);
            assert_eq!(nw.len(), n);
            for pair in nw.windows(2) {
                assert!(
                    pair[1].0 > pair[0].0,
                    "n={n}: nodes not strictly increasing"
                );
            }
            for (x, w) in &nw {
                assert!(*w > 0.0 && w.is_finite(), "n={n}: weight {w} at x={x}");
                assert!(x.abs() < (2.0 * n as f64 + 1.0).sqrt(), "n={n}: node {x}");
            }
        }
    }

    #[test]
    fn test_hermite_weights_sum_to_sqrt_pi_and_reproduce_low_moments() {
        // An n-point rule is exact up to degree 2n-1, so Σ w x^(2k) = Γ(k + 1/2) only
        // for 2k <= 2n-1 (x^4 first from n = 3, x^6 from n = 4).
        for n in 2..=40 {
            let nw = gauss_hermite_weights(n);
            let sum: f64 = nw.iter().map(|(_, w)| w).sum();
            assert_rel_close(sum, PI.sqrt(), 1e-12, &format!("n={n} Σw"));
            for k in 1..=3_usize {
                if 2 * k < 2 * n {
                    let moment: f64 = nw.iter().map(|(x, w)| w * x.powi(2 * k as i32)).sum();
                    assert_rel_close(
                        moment,
                        gaussian_even_moment(k),
                        1e-12,
                        &format!("n={n} Σw·x^{}", 2 * k),
                    );
                }
            }
            let odd_1: f64 = nw.iter().map(|(x, w)| w * x).sum();
            let odd_3: f64 = nw.iter().map(|(x, w)| w * x.powi(3)).sum();
            assert!(odd_1.abs() < 1e-13, "n={n}: Σw·x = {odd_1}");
            assert!(odd_3.abs() < 1e-12, "n={n}: Σw·x³ = {odd_3}");
        }
    }

    #[test]
    fn test_hermite_is_exact_up_to_degree_two_n_minus_one() {
        for n in 1..=25 {
            let nw = gauss_hermite_weights(n);
            for k in 0..n {
                let moment: f64 = nw.iter().map(|(x, w)| w * x.powi(2 * k as i32)).sum();
                assert_rel_close(
                    moment,
                    gaussian_even_moment(k),
                    1e-11,
                    &format!("n={n} degree {}", 2 * k),
                );
            }
        }
    }

    #[test]
    fn test_hermite_rule_is_exactly_symmetric() {
        for n in [2, 3, 8, 9, 20, 31] {
            let nw = gauss_hermite_weights(n);
            for i in 0..n {
                assert_eq!(nw[i].0, -nw[n - 1 - i].0, "n={n}: node {i}");
                assert_eq!(nw[i].1, nw[n - 1 - i].1, "n={n}: weight {i}");
            }
            if n % 2 == 1 {
                assert_eq!(nw[n / 2].0, 0.0);
            }
        }
    }

    #[test]
    fn test_hermite_outer_weights_keep_relative_accuracy() {
        // Smallest weights of the 20- and 40-point rules (50-digit reference values).
        let nw20 = gauss_hermite_weights(20);
        assert_rel_close(nw20[0].1, 2.2293936455341513e-13, 1e-11, "w(n=20)");
        assert!((nw20[19].0 - 5.387480890011233).abs() < 1e-14);
        let nw40 = gauss_hermite_weights(40);
        assert_rel_close(nw40[0].1, 2.591043713847082e-29, 1e-10, "w(n=40)");
        assert!((nw40[39].0 - 8.09876113925085).abs() < 1e-13);
    }

    #[test]
    fn test_hermite_large_rule_is_finite_and_normalised() {
        let n = 150;
        let nw = gauss_hermite_weights(n);
        assert_eq!(nw.len(), n);
        let sum: f64 = nw.iter().map(|(_, w)| w).sum();
        assert_rel_close(sum, PI.sqrt(), 1e-12, "Σw");
        for pair in nw.windows(2) {
            assert!(pair[1].0 > pair[0].0);
        }
        for (x, w) in &nw {
            assert!(x.is_finite() && w.is_finite() && *w >= 0.0, "x={x}, w={w}");
        }
        let second: f64 = nw.iter().map(|(x, w)| w * x * x).sum();
        assert_rel_close(second, PI.sqrt() / 2.0, 1e-12, "Σw·x²");
    }

    #[test]
    fn test_hermite_integrates_smooth_functions() {
        // ∫ cos(x) exp(-x²) dx = √π exp(-1/4)
        let nw = gauss_hermite_weights(12);
        let value: f64 = nw.iter().map(|(x, w)| w * x.cos()).sum();
        assert_rel_close(value, PI.sqrt() * (-0.25_f64).exp(), 1e-13, "∫cos·gauss");
        // Gaussian expectation: E[X²] = 1 for X ~ N(0, 1), substituting x = √2 t.
        let nw = gauss_hermite_weights(5);
        let variance: f64 = nw.iter().map(|(t, w)| w * 2.0 * t * t).sum::<f64>() / PI.sqrt();
        assert_rel_close(variance, 1.0, 1e-13, "E[X²]");
    }

    #[test]
    fn test_hermite_recurrence_stays_finite_for_huge_arguments() {
        // exp(x²/2) overflows f64 here; the rescaled recurrence must stay finite.
        let (q, q_prev, rescales) = hermite_recurrence(2000, 60.0);
        assert!(q.is_finite() && q_prev.is_finite());
        assert!(rescales > 0, "rescales={rescales}");
        let w = hermite_weight(2000, 60.0);
        assert!(w.is_finite() && w >= 0.0, "w={w}");
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
