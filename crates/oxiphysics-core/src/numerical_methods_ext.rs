// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Extended numerical methods: integral equations, BEM, radial basis functions.
//!
//! Provides high-order quadrature, spline interpolation, RBF interpolation,
//! barycentric interpolation, and Richardson-extrapolated differentiation.

use std::f64::consts::FRAC_1_SQRT_2;

/// Gauss-Legendre quadrature rule with `n` nodes on \[-1, 1\].
///
/// Exact for polynomials of degree up to 2n-1.
pub struct GaussLegendreQuad {
    /// Number of quadrature nodes.
    pub n: usize,
    /// Quadrature nodes on \[-1, 1\].
    pub nodes: Vec<f64>,
    /// Quadrature weights summing to 2.
    pub weights: Vec<f64>,
}

impl GaussLegendreQuad {
    /// Construct a Gauss-Legendre rule with `n` points (2 ≤ n ≤ 20).
    pub fn new(n: usize) -> Self {
        let (nodes, weights) = gauss_legendre_nodes_weights(n);
        Self { n, nodes, weights }
    }

    /// Integrate `f` over \[a, b\] using the n-point rule.
    pub fn integrate<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F) -> f64 {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);
        self.nodes
            .iter()
            .zip(self.weights.iter())
            .map(|(&x, &w)| w * f(mid + half * x))
            .sum::<f64>()
            * half
    }

    /// Composite Gauss-Legendre integration with `m` subintervals.
    pub fn integrate_composite<F: Fn(f64) -> f64>(&self, a: f64, b: f64, m: usize, f: &F) -> f64 {
        let h = (b - a) / m as f64;
        (0..m)
            .map(|i| {
                let ai = a + i as f64 * h;
                let bi = ai + h;
                self.integrate(ai, bi, f)
            })
            .sum()
    }
}

/// Compute Gauss-Legendre nodes and weights using hardcoded tables for n ≤ 9,
/// and Newton iteration on Legendre polynomials for larger n.
fn gauss_legendre_nodes_weights(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![0.0], vec![2.0]),
        2 => (
            vec![-0.5773502691896257, 0.5773502691896257],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.7745966692414834, 0.0, 0.7745966692414834],
            vec![0.5555555555555556, 0.8888888888888888, 0.5555555555555556],
        ),
        4 => (
            vec![
                -0.8611363115940526,
                -0.3399810435848563,
                0.3399810435848563,
                0.8611363115940526,
            ],
            vec![
                0.3478548451374538,
                0.6521451548625461,
                0.6521451548625461,
                0.3478548451374538,
            ],
        ),
        5 => (
            vec![
                -0.906_179_845_938_664,
                -0.5384693101056831,
                0.0,
                0.5384693101056831,
                0.906_179_845_938_664,
            ],
            vec![
                0.2369268850561891,
                0.4786286704993665,
                0.5688888888888889,
                0.4786286704993665,
                0.2369268850561891,
            ],
        ),
        6 => (
            vec![
                -0.932_469_514_203_152,
                -0.6612093864662645,
                -0.2386191860831969,
                0.2386191860831969,
                0.6612093864662645,
                0.932_469_514_203_152,
            ],
            vec![
                0.1713244923791704,
                0.3607615730481386,
                0.467_913_934_572_691,
                0.467_913_934_572_691,
                0.3607615730481386,
                0.1713244923791704,
            ],
        ),
        7 => (
            vec![
                -0.9491079123427585,
                -0.7415311855993945,
                -0.4058451513773972,
                0.0,
                0.4058451513773972,
                0.7415311855993945,
                0.9491079123427585,
            ],
            vec![
                0.1294849661688697,
                0.2797053914892767,
                0.3818300505051189,
                0.4179591836734694,
                0.3818300505051189,
                0.2797053914892767,
                0.1294849661688697,
            ],
        ),
        8 => (
            vec![
                -0.9602898564975363,
                -0.7966664774136267,
                -0.525_532_409_916_329,
                -0.1834346424956498,
                0.1834346424956498,
                0.525_532_409_916_329,
                0.7966664774136267,
                0.9602898564975363,
            ],
            vec![
                0.1012285362903763,
                0.2223810344533745,
                0.3137066458778873,
                0.362_683_783_378_362,
                0.362_683_783_378_362,
                0.3137066458778873,
                0.2223810344533745,
                0.1012285362903763,
            ],
        ),
        9 => (
            vec![
                -0.9681602395076261,
                -0.8360311073266358,
                -0.6133714327005904,
                -0.3242534234038089,
                0.0,
                0.3242534234038089,
                0.6133714327005904,
                0.8360311073266358,
                0.9681602395076261,
            ],
            vec![
                0.0812743883615744,
                0.1806481606948574,
                0.2606106964029354,
                0.3123470770400029,
                0.3302393550012598,
                0.3123470770400029,
                0.2606106964029354,
                0.1806481606948574,
                0.0812743883615744,
            ],
        ),
        _ => gauss_legendre_gw(n),
    }
}

/// Compute Gauss-Legendre nodes/weights for n > 9 using Newton iteration
/// on Legendre polynomials (Golub-Welsch style).
fn gauss_legendre_gw(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = Vec::with_capacity(n);
    let mut weights = Vec::with_capacity(n);
    let m = n.div_ceil(2);
    for i in 0..m {
        // Initial guess via asymptotic formula
        let theta = std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5);
        let mut x = theta.cos();
        // Newton iteration
        for _ in 0..100 {
            let (p, dp) = legendre_pn_dpn(n, x);
            let dx = -p / dp;
            x += dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let (_, dp) = legendre_pn_dpn(n, x);
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        nodes.push(-x);
        weights.push(w);
    }
    // Mirror for symmetric pairs
    let mut nodes_full: Vec<f64> = nodes.clone();
    let mut weights_full: Vec<f64> = weights.clone();
    let start = if n % 2 == 1 { m - 1 } else { m };
    for i in (0..start).rev() {
        nodes_full.push(-nodes[i]);
        weights_full.push(weights[i]);
    }
    nodes_full.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // Recompute weights in sorted order to align
    let mut result_nodes = Vec::with_capacity(n);
    let mut result_weights = Vec::with_capacity(n);
    for &x in &nodes_full {
        let (_, dp) = legendre_pn_dpn(n, x);
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        result_nodes.push(x);
        result_weights.push(w);
    }
    (result_nodes, result_weights)
}

/// Evaluate Legendre polynomial P_n(x) and its derivative P_n'(x).
fn legendre_pn_dpn(n: usize, x: f64) -> (f64, f64) {
    let mut p0 = 1.0_f64;
    let mut p1 = x;
    for k in 1..n {
        let kf = k as f64;
        let p2 = ((2.0 * kf + 1.0) * x * p1 - kf * p0) / (kf + 1.0);
        p0 = p1;
        p1 = p2;
    }
    if n == 0 {
        (1.0, 0.0)
    } else {
        let pn = p1;
        let nf = n as f64;
        let dpn = nf * (p0 - x * p1) / (1.0 - x * x);
        (pn, dpn)
    }
}

/// Gauss-Hermite quadrature for ∫f(x)e^{-x²}dx on (-∞, +∞).
pub struct GaussHermiteQuad {
    /// Number of quadrature nodes.
    pub n: usize,
    /// Quadrature nodes.
    pub nodes: Vec<f64>,
    /// Quadrature weights.
    pub weights: Vec<f64>,
}

impl GaussHermiteQuad {
    /// Construct Gauss-Hermite rule with `n` points (up to 20).
    pub fn new(n: usize) -> Self {
        let (nodes, weights) = gauss_hermite_nw(n);
        Self { n, nodes, weights }
    }

    /// Evaluate ∫f(x)e^{-x²}dx ≈ Σ w_i f(x_i).
    pub fn integrate<F: Fn(f64) -> f64>(&self, f: &F) -> f64 {
        self.nodes
            .iter()
            .zip(self.weights.iter())
            .map(|(&x, &w)| w * f(x))
            .sum()
    }
}

/// Compute Gauss-Hermite nodes and weights using hardcoded tables for small n,
/// and Newton iteration on Hermite polynomials for larger n.
fn gauss_hermite_nw(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        2 => (
            vec![-FRAC_1_SQRT_2, FRAC_1_SQRT_2],
            vec![0.886_226_925_452_758, 0.886_226_925_452_758],
        ),
        3 => (
            vec![-1.224_744_871_391_589, 0.0, 1.224_744_871_391_589],
            vec![0.2954089751509193, 1.1816359006036772, 0.2954089751509193],
        ),
        4 => (
            vec![
                -1.6506801238857845,
                -0.5246476232752903,
                0.5246476232752903,
                1.6506801238857845,
            ],
            vec![
                0.08131283544724518,
                0.8049140900055128,
                0.8049140900055128,
                0.08131283544724518,
            ],
        ),
        5 => (
            vec![
                -2.0201828704560856,
                -0.9585724646138185,
                0.0,
                0.9585724646138185,
                2.0201828704560856,
            ],
            vec![
                0.01995324205904592,
                0.3936193231522412,
                0.9453087204829419,
                0.3936193231522412,
                0.01995324205904592,
            ],
        ),
        _ => gauss_hermite_iter(n),
    }
}

/// Compute Gauss-Hermite nodes/weights for n > 5 via Newton iteration on
/// physicist's Hermite polynomials H_n(x).
fn gauss_hermite_iter(n: usize) -> (Vec<f64>, Vec<f64>) {
    let nf = n as f64;
    let m = n.div_ceil(2);
    let mut nodes = Vec::with_capacity(n);
    let mut weights = Vec::with_capacity(n);
    for i in 0..m {
        // Initial guess
        let mut x = (2.0 * nf + 1.0).sqrt()
            - 1.85575
                * (2.0 * nf + 1.0_f64).powf(-1.0 / 6.0)
                * ((i as f64 + 1.0) * std::f64::consts::PI / (nf + 0.5)).cos();
        for _ in 0..100 {
            let (hn, hn1) = hermite_hn_hn1(n, x);
            let dx = -hn / (2.0 * nf * hn1);
            x += dx;
            if dx.abs() < 1e-14 {
                break;
            }
        }
        let (_, hn1) = hermite_hn_hn1(n, x);
        let w = (2.0_f64.powi(n as i32 - 1) * factorial_approx(n) * std::f64::consts::PI.sqrt())
            / (nf * nf * hn1 * hn1);
        nodes.push(-x);
        weights.push(w);
    }
    let start = if n % 2 == 1 { m - 1 } else { m };
    let mut full_nodes = nodes.clone();
    let mut full_weights = weights.clone();
    for i in (0..start).rev() {
        full_nodes.push(-nodes[i]);
        full_weights.push(weights[i]);
    }
    let mut pairs: Vec<(f64, f64)> = full_nodes.into_iter().zip(full_weights).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let (sorted_nodes, sorted_weights) = pairs.into_iter().unzip();
    (sorted_nodes, sorted_weights)
}

/// Evaluate physicist's Hermite polynomial H_n(x) and H_{n-1}(x).
fn hermite_hn_hn1(n: usize, x: f64) -> (f64, f64) {
    let mut h0 = 1.0_f64;
    let mut h1 = 2.0 * x;
    if n == 0 {
        return (h0, 0.0);
    }
    if n == 1 {
        return (h1, h0);
    }
    for k in 1..n {
        let kf = k as f64;
        let h2 = 2.0 * x * h1 - 2.0 * kf * h0;
        h0 = h1;
        h1 = h2;
    }
    (h1, h0)
}

/// Approximate n! via Stirling for weight computation.
fn factorial_approx(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

/// Adaptive trapezoidal rule with Richardson extrapolation.
pub struct TrapezoidalRule {
    /// Absolute tolerance for adaptive refinement.
    pub tol: f64,
    /// Maximum number of refinement levels.
    pub max_levels: usize,
}

impl TrapezoidalRule {
    /// Create with given tolerance and max refinement levels.
    pub fn new(tol: f64, max_levels: usize) -> Self {
        Self { tol, max_levels }
    }

    /// Integrate `f` over \[a, b\] adaptively.
    pub fn integrate<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F) -> f64 {
        let mut h = b - a;
        let mut t = 0.5 * h * (f(a) + f(b));
        let mut prev = t;
        for k in 1..=self.max_levels {
            let n = 1usize << k;
            let new_h = h / 2.0;
            // Add midpoints
            let sum_mid: f64 = (0..n / 2).map(|i| f(a + (2 * i + 1) as f64 * new_h)).sum();
            t = 0.5 * t + new_h * sum_mid;
            h = new_h;
            // Richardson extrapolation: T_rich = (4T - T_prev) / 3
            let t_rich = (4.0 * t - prev) / 3.0;
            if (t_rich - prev).abs() < self.tol * (1.0 + t_rich.abs()) {
                return t_rich;
            }
            prev = t;
        }
        t
    }
}

/// Adaptive Simpson's rule (1/3 rule).
pub struct SimpsonRule {
    /// Absolute tolerance.
    pub tol: f64,
    /// Maximum recursion depth.
    pub max_depth: usize,
}

/// Internal interval state for the adaptive Simpson recursion.
struct SimpsonInterval {
    /// Left endpoint.
    a: f64,
    /// Right endpoint.
    b: f64,
    /// Function value at `a`.
    fa: f64,
    /// Function value at midpoint `(a+b)/2`.
    fm: f64,
    /// Function value at `b`.
    fb: f64,
    /// Current whole-interval Simpson estimate.
    s: f64,
    /// Remaining tolerance for this sub-interval.
    tol: f64,
    /// Remaining recursion depth.
    depth: usize,
}

impl SimpsonRule {
    /// Create with given tolerance.
    pub fn new(tol: f64, max_depth: usize) -> Self {
        Self { tol, max_depth }
    }

    /// Integrate `f` over \[a, b\] using adaptive Simpson.
    pub fn integrate<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F) -> f64 {
        let fa = f(a);
        let fb = f(b);
        let fm = f(0.5 * (a + b));
        let s = simpson13(a, b, fa, fm, fb);
        let iv = SimpsonInterval {
            a,
            b,
            fa,
            fm,
            fb,
            s,
            tol: self.tol,
            depth: self.max_depth,
        };
        Self::recursive(&iv, f)
    }

    fn recursive<F: Fn(f64) -> f64>(iv: &SimpsonInterval, f: &F) -> f64 {
        let mid = 0.5 * (iv.a + iv.b);
        let fml = f(0.5 * (iv.a + mid));
        let fmr = f(0.5 * (mid + iv.b));
        let sl = simpson13(iv.a, mid, iv.fa, fml, iv.fm);
        let sr = simpson13(mid, iv.b, iv.fm, fmr, iv.fb);
        let err = ((sl + sr) - iv.s).abs() / 15.0;
        if iv.depth == 0 || err < iv.tol {
            sl + sr + err
        } else {
            let left = SimpsonInterval {
                a: iv.a,
                b: mid,
                fa: iv.fa,
                fm: fml,
                fb: iv.fm,
                s: sl,
                tol: iv.tol / 2.0,
                depth: iv.depth - 1,
            };
            let right = SimpsonInterval {
                a: mid,
                b: iv.b,
                fa: iv.fm,
                fm: fmr,
                fb: iv.fb,
                s: sr,
                tol: iv.tol / 2.0,
                depth: iv.depth - 1,
            };
            Self::recursive(&left, f) + Self::recursive(&right, f)
        }
    }

    /// Simpson 3/8 rule on \[a, b\] with four points.
    pub fn integrate38<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F) -> f64 {
        let h = (b - a) / 3.0;
        (3.0 * h / 8.0) * (f(a) + 3.0 * f(a + h) + 3.0 * f(a + 2.0 * h) + f(b))
    }
}

#[inline]
fn simpson13(a: f64, b: f64, fa: f64, fm: f64, fb: f64) -> f64 {
    (b - a) / 6.0 * (fa + 4.0 * fm + fb)
}

/// Romberg integration using Richardson extrapolation table.
pub struct RombergIntegration {
    /// Tolerance for convergence.
    pub tol: f64,
    /// Maximum table size (rows).
    pub max_rows: usize,
}

impl RombergIntegration {
    /// Create with given tolerance and max rows.
    pub fn new(tol: f64, max_rows: usize) -> Self {
        Self { tol, max_rows }
    }

    /// Integrate `f` over \[a, b\] using Romberg method. Returns (result, table).
    pub fn integrate<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F) -> (f64, Vec<Vec<f64>>) {
        let mut table: Vec<Vec<f64>> = Vec::new();
        let mut h = b - a;
        // T(0,0)
        table.push(vec![0.5 * h * (f(a) + f(b))]);

        for k in 1..self.max_rows {
            let n = 1usize << k;
            let new_h = h / 2.0;
            let sum: f64 = (0..n / 2).map(|i| f(a + (2 * i + 1) as f64 * new_h)).sum();
            let t = 0.5 * table[k - 1][0] + new_h * sum;
            h = new_h;

            let mut row = vec![t];
            for j in 1..=k {
                let prev = row[j - 1];
                let older = table[k - 1][j - 1];
                let fac = (4u64.pow(j as u32)) as f64;
                row.push((fac * prev - older) / (fac - 1.0));
            }

            let last = *row.last().expect("row is non-empty");
            let prev_last = *table[k - 1].last().expect("table row is non-empty");
            table.push(row);

            if (last - prev_last).abs() < self.tol * (1.0 + last.abs()) {
                return (last, table);
            }
        }
        let last = *table
            .last()
            .expect("table is non-empty")
            .last()
            .expect("row is non-empty");
        (last, table)
    }
}

/// Gauss-Kronrod G7K15 adaptive quadrature.
pub struct GaussKronrod {
    /// Absolute tolerance.
    pub tol: f64,
    /// Maximum recursion depth.
    pub max_depth: usize,
}

impl GaussKronrod {
    /// G7K15 nodes on \[-1,1\] (Kronrod extension of 7-point Gauss).
    const K15_NODES: [f64; 15] = [
        -0.991_455_371_120_813,
        -0.949_107_912_342_758,
        -0.864_864_423_359_769,
        -0.741_531_185_599_394,
        -0.586_087_235_467_691,
        -0.405_845_151_377_397,
        -0.207_784_955_007_898,
        0.0,
        0.207_784_955_007_898,
        0.405_845_151_377_397,
        0.586_087_235_467_691,
        0.741_531_185_599_394,
        0.864_864_423_359_769,
        0.949_107_912_342_758,
        0.991_455_371_120_813,
    ];

    /// K15 weights.
    const K15_WEIGHTS: [f64; 15] = [
        0.022_935_322_010_529,
        0.063_092_092_629_979,
        0.104_790_010_322_250,
        0.140_653_259_715_525,
        0.169_004_726_639_267,
        0.190_350_578_064_785,
        0.204_432_940_075_298,
        0.209_482_141_084_728,
        0.204_432_940_075_298,
        0.190_350_578_064_785,
        0.169_004_726_639_267,
        0.140_653_259_715_525,
        0.104_790_010_322_250,
        0.063_092_092_629_979,
        0.022_935_322_010_529,
    ];

    /// G7 weights (embedded in K15, using even indices 1,3,5,7,9,11,13).
    const G7_WEIGHTS: [f64; 7] = [
        0.129_484_966_168_870,
        0.279_705_391_489_277,
        0.381_830_050_505_119,
        0.417_959_183_673_469,
        0.381_830_050_505_119,
        0.279_705_391_489_277,
        0.129_484_966_168_870,
    ];

    /// G7 node indices in K15.
    const G7_IDX: [usize; 7] = [1, 3, 5, 7, 9, 11, 13];

    /// Create with given tolerance.
    pub fn new(tol: f64, max_depth: usize) -> Self {
        Self { tol, max_depth }
    }

    /// Integrate `f` over \[a, b\] adaptively.
    pub fn integrate<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F) -> f64 {
        self.adaptive(a, b, f, self.max_depth)
    }

    fn adaptive<F: Fn(f64) -> f64>(&self, a: f64, b: f64, f: &F, depth: usize) -> f64 {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);
        let fvals: Vec<f64> = Self::K15_NODES.iter().map(|&x| f(mid + half * x)).collect();

        let k15: f64 = Self::K15_WEIGHTS
            .iter()
            .zip(fvals.iter())
            .map(|(&w, &fv)| w * fv)
            .sum::<f64>()
            * half;

        let g7: f64 = Self::G7_IDX
            .iter()
            .zip(Self::G7_WEIGHTS.iter())
            .map(|(&idx, &w)| w * fvals[idx])
            .sum::<f64>()
            * half;

        let err = (k15 - g7).abs();
        if depth == 0 || err < self.tol * (1.0 + k15.abs()) {
            k15
        } else {
            self.adaptive(a, mid, f, depth - 1) + self.adaptive(mid, b, f, depth - 1)
        }
    }
}

/// Radial basis function types.
#[derive(Debug, Clone, Copy)]
pub enum RbfKind {
    /// Multiquadric: √(r² + c²).
    Multiquadric,
    /// Inverse multiquadric: 1/√(r² + c²).
    InverseMultiquadric,
    /// Gaussian: exp(-r²/(2σ²)).
    Gaussian,
    /// Thin-plate spline: r²·ln(r).
    ThinPlateSpline,
}

/// Radial basis function interpolation in arbitrary dimension.
pub struct RadialBasisFunction {
    /// RBF kernel type.
    pub kind: RbfKind,
    /// Shape parameter c or σ.
    pub shape: f64,
    /// Data points (each row is one point).
    pub centers: Vec<Vec<f64>>,
    /// Fitted coefficients.
    pub coeffs: Vec<f64>,
}

impl RadialBasisFunction {
    /// Create an RBF interpolant fitted to data.
    ///
    /// `points`: N×D matrix of data locations.
    /// `values`: N values at those locations.
    pub fn fit(kind: RbfKind, shape: f64, points: &[Vec<f64>], values: &[f64]) -> Self {
        let n = points.len();
        // Build kernel matrix K[i,j] = phi(||xi - xj||)
        let mut k_mat = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let r = euclidean_dist(&points[i], &points[j]);
                k_mat[i * n + j] = rbf_kernel(kind, r, shape);
            }
        }
        // Solve K * alpha = values via Gaussian elimination
        let coeffs = solve_linear_system(&k_mat, values, n);
        Self {
            kind,
            shape,
            centers: points.to_vec(),
            coeffs,
        }
    }

    /// Evaluate interpolant at a new point `x`.
    pub fn eval(&self, x: &[f64]) -> f64 {
        self.centers
            .iter()
            .zip(self.coeffs.iter())
            .map(|(c, &alpha)| {
                let r = euclidean_dist(x, c);
                alpha * rbf_kernel(self.kind, r, self.shape)
            })
            .sum()
    }
}

fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&ai, &bi)| (ai - bi).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn rbf_kernel(kind: RbfKind, r: f64, c: f64) -> f64 {
    match kind {
        RbfKind::Multiquadric => (r * r + c * c).sqrt(),
        RbfKind::InverseMultiquadric => 1.0 / (r * r + c * c).sqrt(),
        RbfKind::Gaussian => (-(r * r) / (2.0 * c * c)).exp(),
        RbfKind::ThinPlateSpline => {
            if r < 1e-14 {
                0.0
            } else {
                r * r * r.ln()
            }
        }
    }
}

/// Solve A*x = b via Gaussian elimination with partial pivoting.
pub fn solve_linear_system(a_flat: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut mat = a_flat.to_vec();
    let mut rhs = b.to_vec();

    for col in 0..n {
        // Pivot
        let mut max_row = col;
        let mut max_val = mat[col * n + col].abs();
        for row in (col + 1)..n {
            let v = mat[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_row != col {
            for k in 0..n {
                mat.swap(col * n + k, max_row * n + k);
            }
            rhs.swap(col, max_row);
        }
        let pivot = mat[col * n + col];
        if pivot.abs() < 1e-14 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = mat[row * n + col] / pivot;
            for k in col..n {
                let v = mat[col * n + k];
                mat[row * n + k] -= factor * v;
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    // Back substitution
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n {
            s -= mat[i * n + j] * x[j];
        }
        x[i] = s / mat[i * n + i];
    }
    x
}

/// Natural cubic spline interpolation.
pub struct SplineInterpolation {
    /// Breakpoints (sorted ascending).
    pub xs: Vec<f64>,
    /// Values at breakpoints.
    pub ys: Vec<f64>,
    /// Second derivatives at breakpoints (from tridiagonal solve).
    pub m: Vec<f64>,
    /// Spline boundary condition type.
    pub kind: SplineKind,
}

/// Spline boundary conditions.
#[derive(Debug, Clone, Copy)]
pub enum SplineKind {
    /// Natural: second derivative = 0 at endpoints.
    Natural,
    /// Clamped: prescribed first derivative at endpoints.
    Clamped(f64, f64),
    /// Not-a-knot: third derivative continuous at second and second-to-last knots.
    NotAKnot,
}

impl SplineInterpolation {
    /// Fit cubic spline to data `(xs, ys)` with given boundary condition.
    pub fn fit(xs: &[f64], ys: &[f64], kind: SplineKind) -> Self {
        let n = xs.len();
        assert!(n >= 3, "Need at least 3 points for cubic spline");
        let m = compute_spline_second_derivs(xs, ys, kind);
        Self {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
            m,
            kind,
        }
    }

    /// Evaluate spline at x.
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.xs.len();
        // Find interval
        let mut idx = 0;
        for i in 1..n - 1 {
            if x >= self.xs[i] {
                idx = i;
            }
        }
        let h = self.xs[idx + 1] - self.xs[idx];
        let t = (x - self.xs[idx]) / h;
        let a = self.ys[idx];
        let b_val = self.ys[idx + 1];
        let ma = self.m[idx];
        let mb = self.m[idx + 1];
        // Cubic Hermite form using second derivatives
        (1.0 - t) * a
            + t * b_val
            + t * (1.0 - t)
                * ((1.0 - t) * h * h * ma / 6.0 - t * h * h * mb / 6.0 - h * h * (ma - mb) / 6.0
                    + h * h * ma / 6.0 * (2.0 * t - 1.0))
    }

    /// Evaluate spline using proper cubic formula.
    pub fn eval_cubic(&self, x: f64) -> f64 {
        let n = self.xs.len();
        let mut idx = 0;
        for i in 1..n - 1 {
            if x >= self.xs[i] {
                idx = i;
            }
        }
        let h = self.xs[idx + 1] - self.xs[idx];
        let dx = x - self.xs[idx];
        let dx2 = self.xs[idx + 1] - x;
        let a = self.ys[idx];
        let b_val = self.ys[idx + 1];
        let ma = self.m[idx];
        let mb = self.m[idx + 1];
        (dx2 / h) * a
            + (dx / h) * b_val
            + (dx2 / h) * ((dx2 * dx2 / (h * h) - 1.0) * h * h * ma / 6.0)
            + (dx / h) * ((dx * dx / (h * h) - 1.0) * h * h * mb / 6.0)
    }
}

/// Solve the tridiagonal system for cubic spline second derivatives.
fn compute_spline_second_derivs(xs: &[f64], ys: &[f64], kind: SplineKind) -> Vec<f64> {
    let n = xs.len();
    let mut h = vec![0.0f64; n - 1];
    for i in 0..n - 1 {
        h[i] = xs[i + 1] - xs[i];
    }

    // Build tridiagonal system
    let mut diag = vec![0.0f64; n];
    let mut upper = vec![0.0f64; n];
    let mut lower = vec![0.0f64; n];
    let mut rhs = vec![0.0f64; n];

    for i in 1..n - 1 {
        lower[i] = h[i - 1];
        diag[i] = 2.0 * (h[i - 1] + h[i]);
        upper[i] = h[i];
        rhs[i] = 6.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]);
    }

    match kind {
        SplineKind::Natural => {
            diag[0] = 1.0;
            upper[0] = 0.0;
            rhs[0] = 0.0;
            diag[n - 1] = 1.0;
            lower[n - 1] = 0.0;
            rhs[n - 1] = 0.0;
        }
        SplineKind::Clamped(d0, dn) => {
            diag[0] = 2.0 * h[0];
            upper[0] = h[0];
            rhs[0] = 6.0 * ((ys[1] - ys[0]) / h[0] - d0);
            diag[n - 1] = 2.0 * h[n - 2];
            lower[n - 1] = h[n - 2];
            rhs[n - 1] = 6.0 * (dn - (ys[n - 1] - ys[n - 2]) / h[n - 2]);
        }
        SplineKind::NotAKnot => {
            // Third derivative continuous at xs[1]: h[1]*m[0] = h[0]*m[2] - (h[0]+h[1])*m[1]
            diag[0] = h[1];
            upper[0] = -(h[0] + h[1]);
            rhs[0] = 0.0;
            // ... also at xs[n-2]
            diag[n - 1] = h[n - 3];
            lower[n - 1] = -(h[n - 3] + h[n - 2]);
            rhs[n - 1] = 0.0;
        }
    }

    // Solve tridiagonal via Thomas algorithm
    thomas_algorithm(&diag, &lower, &upper, &rhs, n)
}

/// Thomas algorithm (tridiagonal matrix algorithm).
fn thomas_algorithm(diag: &[f64], lower: &[f64], upper: &[f64], rhs: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; n];
    let mut d = rhs.to_vec();
    let mut w = vec![0.0f64; n];

    w[0] = upper[0] / diag[0];
    d[0] /= diag[0];

    for i in 1..n {
        let denom = diag[i] - lower[i] * w[i - 1];
        if denom.abs() < 1e-14 {
            w[i] = 0.0;
            d[i] = 0.0;
            continue;
        }
        w[i] = upper[i] / denom;
        d[i] = (d[i] - lower[i] * d[i - 1]) / denom;
    }

    c[n - 1] = d[n - 1];
    for i in (0..n - 1).rev() {
        c[i] = d[i] - w[i] * c[i + 1];
    }
    c
}

/// Barycentric interpolation (Berrut second form).
pub struct BarycentricInterpolation {
    /// Interpolation nodes.
    pub nodes: Vec<f64>,
    /// Function values at nodes.
    pub values: Vec<f64>,
    /// Barycentric weights.
    pub weights: Vec<f64>,
}

impl BarycentricInterpolation {
    /// Fit barycentric interpolant to `(nodes, values)`.
    /// Uses Chebyshev-based weights for numerical stability.
    pub fn fit(nodes: &[f64], values: &[f64]) -> Self {
        let n = nodes.len();
        let weights = compute_barycentric_weights(nodes, n);
        Self {
            nodes: nodes.to_vec(),
            values: values.to_vec(),
            weights,
        }
    }

    /// Evaluate interpolant at x using second barycentric formula.
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.nodes.len();
        // Check if x is a node
        for i in 0..n {
            if (x - self.nodes[i]).abs() < 1e-14 {
                return self.values[i];
            }
        }
        let numer: f64 = (0..n)
            .map(|i| self.weights[i] / (x - self.nodes[i]) * self.values[i])
            .sum();
        let denom: f64 = (0..n).map(|i| self.weights[i] / (x - self.nodes[i])).sum();
        numer / denom
    }
}

fn compute_barycentric_weights(nodes: &[f64], n: usize) -> Vec<f64> {
    let mut w = vec![1.0f64; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                w[i] /= nodes[i] - nodes[j];
            }
        }
    }
    w
}

/// Richardson-extrapolated numerical differentiation.
pub struct NumericalDifferentiation {
    /// Base step size.
    pub h: f64,
    /// Number of Richardson extrapolation levels.
    pub levels: usize,
}

impl NumericalDifferentiation {
    /// Create with base step and extrapolation levels.
    pub fn new(h: f64, levels: usize) -> Self {
        Self { h, levels }
    }

    /// First derivative via central differences + Richardson extrapolation.
    pub fn first_deriv<F: Fn(f64) -> f64>(&self, x: f64, f: &F) -> f64 {
        richardson_deriv(x, self.h, self.levels, f, 1)
    }

    /// Second derivative via central differences + Richardson extrapolation.
    pub fn second_deriv<F: Fn(f64) -> f64>(&self, x: f64, f: &F) -> f64 {
        richardson_deriv2(x, self.h, self.levels, f)
    }

    /// Mixed partial ∂²f/∂x∂y at (x, y).
    pub fn mixed_partial<F: Fn(f64, f64) -> f64>(&self, x: f64, y: f64, f: &F) -> f64 {
        let h = self.h;
        let fpp = f(x + h, y + h);
        let fpm = f(x + h, y - h);
        let fmp = f(x - h, y + h);
        let fmm = f(x - h, y - h);
        (fpp - fpm - fmp + fmm) / (4.0 * h * h)
    }

    /// Arbitrary-order derivative (order 1..4) via finite differences.
    pub fn nth_deriv<F: Fn(f64) -> f64>(&self, x: f64, order: usize, f: &F) -> f64 {
        match order {
            1 => self.first_deriv(x, f),
            2 => self.second_deriv(x, f),
            3 => {
                let h = self.h;
                (f(x + 2.0 * h) - 2.0 * f(x + h) + 2.0 * f(x - h) - f(x - 2.0 * h))
                    / (2.0 * h * h * h)
            }
            4 => {
                let h = self.h;
                (f(x + 2.0 * h) - 4.0 * f(x + h) + 6.0 * f(x) - 4.0 * f(x - h) + f(x - 2.0 * h))
                    / (h * h * h * h)
            }
            _ => panic!("Only orders 1-4 supported"),
        }
    }
}

fn richardson_deriv<F: Fn(f64) -> f64>(
    x: f64,
    h0: f64,
    levels: usize,
    f: &F,
    _order: usize,
) -> f64 {
    let mut table: Vec<Vec<f64>> = Vec::new();
    let mut h = h0;
    for k in 0..levels {
        let d = (f(x + h) - f(x - h)) / (2.0 * h);
        if k == 0 {
            table.push(vec![d]);
        } else {
            let mut row = vec![d];
            for j in 1..=k {
                let fac = (4u64.pow(j as u32)) as f64;
                let val = (fac * row[j - 1] - table[k - 1][j - 1]) / (fac - 1.0);
                row.push(val);
            }
            table.push(row);
        }
        h /= 2.0;
    }
    *table
        .last()
        .expect("table is non-empty")
        .last()
        .expect("row is non-empty")
}

fn richardson_deriv2<F: Fn(f64) -> f64>(x: f64, h0: f64, levels: usize, f: &F) -> f64 {
    let mut table: Vec<Vec<f64>> = Vec::new();
    let mut h = h0;
    for k in 0..levels {
        let d = (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
        if k == 0 {
            table.push(vec![d]);
        } else {
            let mut row = vec![d];
            for j in 1..=k {
                let fac = (4u64.pow(j as u32)) as f64;
                let val = (fac * row[j - 1] - table[k - 1][j - 1]) / (fac - 1.0);
                row.push(val);
            }
            table.push(row);
        }
        h /= 2.0;
    }
    *table
        .last()
        .expect("table is non-empty")
        .last()
        .expect("row is non-empty")
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper trait for tests only
    trait DiffExact {
        fn diff_exact(self) -> f64;
    }

    impl DiffExact for f64 {
        fn diff_exact(self) -> f64 {
            3.0 * self * self
        }
    }

    // ── Gauss-Legendre tests ──────────────────────────────────────────────

    #[test]
    fn test_gl_exact_degree3() {
        // 2-point GL is exact for degree 2*2-1=3 polynomials
        let gl = GaussLegendreQuad::new(2);
        let result = gl.integrate(0.0, 1.0, &|x| x * x * x);
        assert!((result - 0.25).abs() < 1e-12, "GL2 x^3: {result:.6}");
    }

    #[test]
    fn test_gl_exact_degree5() {
        // 3-point GL exact for degree 5
        let gl = GaussLegendreQuad::new(3);
        let result = gl.integrate(-1.0, 1.0, &|x| x.powi(5));
        assert!(result.abs() < 1e-12, "GL3 x^5 on [-1,1]: {result:.6}");
    }

    #[test]
    fn test_gl_sine() {
        let gl = GaussLegendreQuad::new(5);
        let result = gl.integrate(0.0, std::f64::consts::PI, &|x| x.sin());
        assert!((result - 2.0).abs() < 1e-4, "GL5 sin: {result:.12}");
    }

    #[test]
    fn test_gl_composite_sine() {
        let gl = GaussLegendreQuad::new(3);
        let result = gl.integrate_composite(0.0, std::f64::consts::PI, 10, &|x| x.sin());
        assert!((result - 2.0).abs() < 1e-8, "GL composite sin: {result:.6}");
    }

    #[test]
    fn test_gl_nodes_count() {
        for n in [2, 3, 5, 10] {
            let gl = GaussLegendreQuad::new(n);
            assert_eq!(gl.nodes.len(), n);
            assert_eq!(gl.weights.len(), n);
        }
    }

    #[test]
    fn test_gl_weights_sum_to_2() {
        for n in [2, 3, 4, 5] {
            let gl = GaussLegendreQuad::new(n);
            let s: f64 = gl.weights.iter().sum();
            assert!((s - 2.0).abs() < 1e-10, "GL{n} weights sum {s:.6}");
        }
    }

    // ── Gauss-Hermite tests ───────────────────────────────────────────────

    #[test]
    fn test_hermite_constant() {
        // ∫ e^{-x²} dx = sqrt(pi)
        let gh = GaussHermiteQuad::new(5);
        let result = gh.integrate(&|_x| 1.0);
        assert!(
            (result - std::f64::consts::PI.sqrt()).abs() < 1e-8,
            "GH constant: {result:.6}"
        );
    }

    #[test]
    fn test_hermite_x_squared() {
        // ∫ x² e^{-x²} dx = sqrt(pi)/2
        let gh = GaussHermiteQuad::new(5);
        let result = gh.integrate(&|x| x * x);
        let exact = std::f64::consts::PI.sqrt() / 2.0;
        assert!(
            (result - exact).abs() < 1e-8,
            "GH x^2: {result:.6} vs {exact:.6}"
        );
    }

    // ── Trapezoidal tests ─────────────────────────────────────────────────

    #[test]
    fn test_trapezoidal_exp() {
        let trap = TrapezoidalRule::new(1e-10, 20);
        let result = trap.integrate(0.0, 1.0, &|x| x.exp());
        let exact = std::f64::consts::E - 1.0;
        assert!((result - exact).abs() < 1e-9, "Trap exp: {result:.6}");
    }

    #[test]
    fn test_trapezoidal_polynomial() {
        let trap = TrapezoidalRule::new(1e-12, 20);
        let result = trap.integrate(0.0, 1.0, &|x| x * x);
        assert!((result - 1.0 / 3.0).abs() < 1e-10, "Trap x^2: {result:.6}");
    }

    // ── Simpson tests ─────────────────────────────────────────────────────

    #[test]
    fn test_simpson_exp() {
        let simp = SimpsonRule::new(1e-10, 20);
        let result = simp.integrate(0.0, 1.0, &|x| x.exp());
        let exact = std::f64::consts::E - 1.0;
        assert!((result - exact).abs() < 1e-8, "Simpson exp: {result:.6}");
    }

    #[test]
    fn test_simpson38() {
        let simp = SimpsonRule::new(1e-10, 20);
        let result = simp.integrate38(0.0, 1.0, &|x| x * x * x);
        assert!((result - 0.25).abs() < 1e-10, "Simpson38 x^3: {result:.6}");
    }

    // ── Romberg tests ─────────────────────────────────────────────────────

    #[test]
    fn test_romberg_converges_faster() {
        let rom = RombergIntegration::new(1e-12, 15);
        let (result, table) = rom.integrate(0.0, 1.0, &|x| x.exp());
        let exact = std::f64::consts::E - 1.0;
        assert!((result - exact).abs() < 1e-11, "Romberg exp: {result:.12}");
        // Romberg table should have multiple rows
        assert!(table.len() >= 2);
    }

    #[test]
    fn test_romberg_sine() {
        let rom = RombergIntegration::new(1e-12, 15);
        let (result, _) = rom.integrate(0.0, std::f64::consts::PI, &|x| x.sin());
        assert!((result - 2.0).abs() < 1e-11, "Romberg sin: {result:.12}");
    }

    #[test]
    fn test_romberg_table_diagonal_improves() {
        let rom = RombergIntegration::new(1e-14, 8);
        let (_, table) = rom.integrate(0.0, 1.0, &|x| x * x);
        // Last diagonal entry should be very close to 1/3
        let last = *table.last().unwrap().last().unwrap();
        assert!(
            (last - 1.0 / 3.0).abs() < 1e-12,
            "Romberg table: {last:.14}"
        );
    }

    // ── Gauss-Kronrod tests ───────────────────────────────────────────────

    #[test]
    fn test_gauss_kronrod_exp() {
        let gk = GaussKronrod::new(1e-10, 10);
        let result = gk.integrate(0.0, 1.0, &|x| x.exp());
        let exact = std::f64::consts::E - 1.0;
        assert!((result - exact).abs() < 1e-9, "GK exp: {result:.6}");
    }

    #[test]
    fn test_gauss_kronrod_oscillatory() {
        let gk = GaussKronrod::new(1e-8, 15);
        let result = gk.integrate(0.0, 10.0 * std::f64::consts::PI, &|x| x.sin());
        // Integral over full periods = 0
        assert!(result.abs() < 1e-6, "GK oscillatory: {result:.6}");
    }

    // ── RBF tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_rbf_exact_at_data_points_gaussian() {
        let pts: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64]).collect();
        let vals: Vec<f64> = pts.iter().map(|p| p[0] * p[0]).collect();
        let rbf = RadialBasisFunction::fit(RbfKind::Gaussian, 1.0, &pts, &vals);
        for (i, pt) in pts.iter().enumerate() {
            let v = rbf.eval(pt);
            assert!(
                (v - vals[i]).abs() < 1e-8,
                "RBF Gaussian at data[{i}]: {v:.6} vs {:.6}",
                vals[i]
            );
        }
    }

    #[test]
    fn test_rbf_exact_at_data_points_mq() {
        let pts: Vec<Vec<f64>> = vec![vec![0.0], vec![0.5], vec![1.0]];
        let vals = vec![0.0, 0.25, 1.0];
        let rbf = RadialBasisFunction::fit(RbfKind::Multiquadric, 0.5, &pts, &vals);
        for (i, pt) in pts.iter().enumerate() {
            let v = rbf.eval(pt);
            assert!((v - vals[i]).abs() < 1e-8, "RBF MQ at data[{i}]: {v:.6}");
        }
    }

    #[test]
    fn test_rbf_inv_mq() {
        let pts: Vec<Vec<f64>> = vec![vec![0.0], vec![1.0], vec![2.0]];
        let vals = vec![1.0, 2.0, 4.0];
        let rbf = RadialBasisFunction::fit(RbfKind::InverseMultiquadric, 0.5, &pts, &vals);
        for (i, pt) in pts.iter().enumerate() {
            let v = rbf.eval(pt);
            assert!((v - vals[i]).abs() < 1e-8, "RBF InvMQ at data[{i}]: {v:.6}");
        }
    }

    // ── Spline tests ──────────────────────────────────────────────────────

    #[test]
    fn test_spline_natural_at_nodes() {
        let xs: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| x.sin()).collect();
        let sp = SplineInterpolation::fit(&xs, &ys, SplineKind::Natural);
        for i in 0..xs.len() {
            let v = sp.eval_cubic(xs[i]);
            assert!(
                (v - ys[i]).abs() < 1e-10,
                "Natural spline at node {i}: {v:.6} vs {:.6}",
                ys[i]
            );
        }
    }

    #[test]
    fn test_spline_natural_zero_second_deriv_endpoints() {
        let xs: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| x * x).collect();
        let sp = SplineInterpolation::fit(&xs, &ys, SplineKind::Natural);
        assert!(sp.m[0].abs() < 1e-10, "Natural spline M[0]={:.6}", sp.m[0]);
        let n = sp.m.len();
        assert!(
            sp.m[n - 1].abs() < 1e-10,
            "Natural spline M[n-1]={:.6}",
            sp.m[n - 1]
        );
    }

    #[test]
    fn test_spline_clamped_at_nodes() {
        let xs: Vec<f64> = (0..5).map(|i| i as f64 * 0.5).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| x * x * x).collect();
        let sp = SplineInterpolation::fit(&xs, &ys, SplineKind::Clamped(0.0, 3.0 * 2.0 * 2.0));
        for i in 0..xs.len() {
            let v = sp.eval_cubic(xs[i]);
            assert!(
                (v - ys[i]).abs() < 1e-8,
                "Clamped spline at node {i}: {v:.6}"
            );
        }
    }

    #[test]
    fn test_spline_not_a_knot_at_nodes() {
        let xs: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| x.exp()).collect();
        let sp = SplineInterpolation::fit(&xs, &ys, SplineKind::NotAKnot);
        for i in 0..xs.len() {
            let v = sp.eval_cubic(xs[i]);
            assert!(
                (v - ys[i]).abs() < 1e-8,
                "NotAKnot spline at node {i}: {v:.6}"
            );
        }
    }

    // ── Barycentric tests ─────────────────────────────────────────────────

    #[test]
    fn test_barycentric_exact_at_nodes() {
        let nodes: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let values: Vec<f64> = nodes.iter().map(|&x| x * x).collect();
        let bary = BarycentricInterpolation::fit(&nodes, &values);
        for i in 0..nodes.len() {
            let v = bary.eval(nodes[i]);
            assert!(
                (v - values[i]).abs() < 1e-10,
                "Barycentric at node {i}: {v:.6}"
            );
        }
    }

    #[test]
    fn test_barycentric_polynomial_reconstruction() {
        // Polynomial p(x) = x^3 - 2x + 1, interpolated at 4 nodes
        let nodes = vec![0.0, 1.0, 2.0, 3.0];
        let values: Vec<f64> = nodes.iter().map(|&x| x * x * x - 2.0 * x + 1.0).collect();
        let bary = BarycentricInterpolation::fit(&nodes, &values);
        let test_pts = [0.5, 1.5, 2.5];
        for &x in &test_pts {
            let v = bary.eval(x);
            let exact = x * x * x - 2.0 * x + 1.0;
            assert!(
                (v - exact).abs() < 1e-10,
                "Barycentric poly at {x}: {v:.6} vs {exact:.6}"
            );
        }
    }

    // ── Richardson differentiation tests ─────────────────────────────────

    #[test]
    fn test_first_deriv_sine() {
        let nd = NumericalDifferentiation::new(1e-4, 4);
        let x = 1.0;
        let d = nd.first_deriv(x, &|t| t.sin());
        let exact = x.cos();
        assert!(
            (d - exact).abs() < 1e-9,
            "d/dx sin(1): {d:.9} vs {exact:.9}"
        );
    }

    #[test]
    fn test_second_deriv_exp() {
        let nd = NumericalDifferentiation::new(1e-4, 4);
        let x = 0.5;
        let d2 = nd.second_deriv(x, &|t| t.exp());
        assert!((d2 - x.exp()).abs() < 1e-5, "d^2/dx^2 exp(0.5): {d2:.6}");
    }

    #[test]
    fn test_mixed_partial_xy() {
        let nd = NumericalDifferentiation::new(1e-4, 4);
        // f(x,y) = x*y^2, df/dxdy = 2y
        let y = 2.0;
        let d = nd.mixed_partial(1.0, y, &|x, yy| x * yy * yy);
        let exact = 2.0 * y;
        assert!(
            (d - exact).abs() < 1e-7,
            "mixed partial: {d:.6} vs {exact:.6}"
        );
    }

    #[test]
    fn test_first_deriv_accuracy_order4() {
        // Richardson should give O(h^4) accuracy
        let nd_coarse = NumericalDifferentiation::new(1e-2, 1);
        let nd_fine = NumericalDifferentiation::new(1e-2, 4);
        let x = 1.2;
        let exact = x.diff_exact(); // 3x^2
        let d_coarse = nd_coarse.first_deriv(x, &|t| t * t * t);
        let d_fine = nd_fine.first_deriv(x, &|t| t * t * t);
        // Fine should be much more accurate
        assert!(
            (d_fine - exact).abs() < (d_coarse - exact).abs(),
            "Richardson should improve: coarse err={:.2e} fine err={:.2e}",
            (d_coarse - exact).abs(),
            (d_fine - exact).abs()
        );
    }

    #[test]
    fn test_nth_deriv_third_order() {
        let nd = NumericalDifferentiation::new(1e-3, 4);
        // f(x) = x^4, f'''(x) = 24x
        let x = 1.0;
        let d3 = nd.nth_deriv(x, 3, &|t| t.powi(4));
        let exact = 24.0 * x;
        assert!(
            (d3 - exact).abs() < 1e-5,
            "d^3/dx^3 x^4 at 1: {d3:.6} vs {exact:.6}"
        );
    }

    #[test]
    fn test_nth_deriv_fourth_order() {
        let nd = NumericalDifferentiation::new(1e-3, 4);
        // f(x) = x^4, f''''(x) = 24
        let x = 1.0;
        let d4 = nd.nth_deriv(x, 4, &|t| t.powi(4));
        assert!((d4 - 24.0).abs() < 1e-2, "d^4/dx^4 x^4 at 1: {d4:.6}");
    }
}
