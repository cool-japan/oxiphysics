//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{dist3, hermite, hermite_deriv, lerp};

/// A piecewise cubic Hermite spline defined by control points and tangents.
///
/// The spline is parameterized by `t` ∈ \[0, n−1\] where `n` is the number of
/// control points.  At each integer `t = i` the curve passes exactly through
/// `points[i]` with derivative `tangents[i]`.
pub struct HermiteSpline {
    /// Control point values.
    pub points: Vec<f64>,
    /// Tangent (derivative) at each control point.
    pub tangents: Vec<f64>,
}
impl HermiteSpline {
    /// Construct a new `HermiteSpline`.
    ///
    /// Returns an error if `points` and `tangents` have different lengths or
    /// if fewer than two points are provided.
    pub fn new(points: Vec<f64>, tangents: Vec<f64>) -> Result<Self, String> {
        if points.len() != tangents.len() {
            return Err(format!(
                "points ({}) and tangents ({}) must have equal length",
                points.len(),
                tangents.len()
            ));
        }
        if points.len() < 2 {
            return Err("HermiteSpline requires at least 2 points".into());
        }
        Ok(Self { points, tangents })
    }
    /// Evaluate the spline at parameter `t` ∈ \[0, n−1\].
    pub fn evaluate(&self, t: f64) -> f64 {
        let n = self.points.len();
        let t = t.clamp(0.0, (n - 1) as f64);
        let seg = (t.floor() as usize).min(n - 2);
        let lt = t - seg as f64;
        hermite(
            self.points[seg],
            self.tangents[seg],
            self.points[seg + 1],
            self.tangents[seg + 1],
            lt,
        )
    }
    /// Evaluate the derivative of the spline at parameter `t` ∈ \[0, n−1\].
    pub fn derivative(&self, t: f64) -> f64 {
        let n = self.points.len();
        let t = t.clamp(0.0, (n - 1) as f64);
        let seg = (t.floor() as usize).min(n - 2);
        let lt = t - seg as f64;
        hermite_deriv(
            self.points[seg],
            self.tangents[seg],
            self.points[seg + 1],
            self.tangents[seg + 1],
            lt,
        )
    }
}
/// A Catmull-Rom spline through an ordered sequence of 3-D control points.
///
/// The `alpha` parameter controls the parameterization:
/// - `0.0` — uniform (classic)
/// - `0.5` — centripetal (avoids self-intersections)
/// - `1.0` — chordal
pub struct CatmullRomSpline {
    /// Control points in 3-D space.
    pub points: Vec<[f64; 3]>,
    /// Parameterization exponent (0 = uniform, 0.5 = centripetal, 1 = chordal).
    pub alpha: f64,
}
impl CatmullRomSpline {
    /// Create a new spline.  At least two points are required for a meaningful
    /// result; fewer produce degenerate (constant) output.
    pub fn new(points: Vec<[f64; 3]>, alpha: f64) -> Self {
        Self { points, alpha }
    }
    /// Evaluate the spline at global parameter `t` ∈ \[0, n−1\].
    ///
    /// Each integer value of `t` corresponds to a control point.
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let n = self.points.len();
        if n == 0 {
            return [0.0, 0.0, 0.0];
        }
        if n == 1 {
            return self.points[0];
        }
        let t = t.clamp(0.0, (n - 1) as f64);
        let seg = (t.floor() as usize).min(n - 2);
        let local_t = t - seg as f64;
        let i0 = seg.saturating_sub(1);
        let i1 = seg;
        let i2 = (seg + 1).min(n - 1);
        let i3 = (seg + 2).min(n - 1);
        let t0 = 0.0_f64;
        let t1 = t0 + dist3(self.points[i0], self.points[i1]).powf(self.alpha);
        let t2 = t1 + dist3(self.points[i1], self.points[i2]).powf(self.alpha);
        let t3 = t2 + dist3(self.points[i2], self.points[i3]).powf(self.alpha);
        let u = t1 + local_t * (t2 - t1);
        self.barry_goldman(
            self.points[i0],
            self.points[i1],
            self.points[i2],
            self.points[i3],
            t0,
            t1,
            t2,
            t3,
            u,
        )
    }
    /// Barry-Goldman algorithm for parameterized Catmull-Rom evaluation.
    #[allow(clippy::too_many_arguments)]
    fn barry_goldman(
        &self,
        p0: [f64; 3],
        p1: [f64; 3],
        p2: [f64; 3],
        p3: [f64; 3],
        t0: f64,
        t1: f64,
        t2: f64,
        t3: f64,
        t: f64,
    ) -> [f64; 3] {
        fn blend(a: [f64; 3], b: [f64; 3], ta: f64, tb: f64, t: f64) -> [f64; 3] {
            if (tb - ta).abs() < f64::EPSILON {
                return a;
            }
            let f = (t - ta) / (tb - ta);
            [
                lerp(a[0], b[0], f),
                lerp(a[1], b[1], f),
                lerp(a[2], b[2], f),
            ]
        }
        let a1 = blend(p0, p1, t0, t1, t);
        let a2 = blend(p1, p2, t1, t2, t);
        let a3 = blend(p2, p3, t2, t3, t);
        let b1 = blend(a1, a2, t0, t2, t);
        let b2 = blend(a2, a3, t1, t3, t);
        blend(b1, b2, t1, t2, t)
    }
    /// Approximate arc length by summing chord lengths over `n_samples` samples.
    pub fn arc_length(&self, n_samples: usize) -> f64 {
        let n = self.points.len();
        if n < 2 || n_samples < 2 {
            return 0.0;
        }
        let t_max = (n - 1) as f64;
        let mut total = 0.0;
        let mut prev = self.evaluate(0.0);
        for i in 1..n_samples {
            let t = t_max * (i as f64) / ((n_samples - 1) as f64);
            let cur = self.evaluate(t);
            total += dist3(prev, cur);
            prev = cur;
        }
        total
    }
    /// Evaluate the spline at arc-length parameter `s` (approximate).
    ///
    /// Uses `n_samples` samples to build a look-up table, then linearly
    /// interpolates between bracketing samples.
    pub fn evaluate_at_arc_length(&self, s: f64, n_samples: usize) -> [f64; 3] {
        let n = self.points.len();
        if n < 2 || n_samples < 2 {
            return self.points.first().copied().unwrap_or([0.0; 3]);
        }
        let t_max = (n - 1) as f64;
        let mut ts = Vec::with_capacity(n_samples);
        let mut lengths = Vec::with_capacity(n_samples);
        let mut prev = self.evaluate(0.0);
        let mut cum = 0.0;
        ts.push(0.0_f64);
        lengths.push(0.0_f64);
        for i in 1..n_samples {
            let t = t_max * (i as f64) / ((n_samples - 1) as f64);
            let cur = self.evaluate(t);
            cum += dist3(prev, cur);
            ts.push(t);
            lengths.push(cum);
            prev = cur;
        }
        let total = cum;
        let s = s.clamp(0.0, total);
        let idx = lengths.partition_point(|&l| l <= s).saturating_sub(1);
        let idx = idx.min(n_samples - 2);
        let l0 = lengths[idx];
        let l1 = lengths[idx + 1];
        let frac = if (l1 - l0).abs() < f64::EPSILON {
            0.0
        } else {
            (s - l0) / (l1 - l0)
        };
        let t = lerp(ts[idx], ts[idx + 1], frac);
        self.evaluate(t)
    }
    /// Backward-compatible arc-length approximation (same as `arc_length`).
    pub fn arc_length_approx(&self, n_samples: usize) -> f64 {
        self.arc_length(n_samples)
    }
}
/// A natural (free-end) cubic spline that passes through a set of (x, y) knots.
///
/// Knot x-positions must be strictly increasing.  The second-derivative
/// boundary conditions are zero at both ends (natural spline).
pub struct NaturalCubicSpline {
    /// Sorted knot x-positions.
    pub xs: Vec<f64>,
    /// Values at the knots.
    pub ys: Vec<f64>,
    /// Per-segment cubic coefficients `[a, b, c, d]` such that on the i-th
    /// segment the polynomial is `a + b·h + c·h² + d·h³` with `h = x − xs[i]`.
    pub coeffs: Vec<[f64; 4]>,
}
impl NaturalCubicSpline {
    /// Fit a natural cubic spline to the given knots.
    ///
    /// Returns an error if fewer than two knots are provided or if the x-values
    /// are not strictly increasing.
    pub fn fit(xs: Vec<f64>, ys: Vec<f64>) -> Result<Self, String> {
        let n = xs.len();
        if n < 2 {
            return Err("NaturalCubicSpline requires at least 2 knots".into());
        }
        if xs.len() != ys.len() {
            return Err("xs and ys must have the same length".into());
        }
        for i in 1..n {
            if xs[i] <= xs[i - 1] {
                return Err(format!(
                    "xs must be strictly increasing: xs[{}]={} <= xs[{}]={}",
                    i,
                    xs[i],
                    i - 1,
                    xs[i - 1]
                ));
            }
        }
        let m = n - 1;
        let mut h: Vec<f64> = (0..m).map(|i| xs[i + 1] - xs[i]).collect();
        if n == 2 {
            let slope = (ys[1] - ys[0]) / h[0];
            let coeffs = vec![[ys[0], slope, 0.0, 0.0]];
            return Ok(Self { xs, ys, coeffs });
        }
        let inner = n - 2;
        let mut diag: Vec<f64> = Vec::with_capacity(inner);
        let mut upper: Vec<f64> = Vec::with_capacity(inner - 1);
        let mut lower: Vec<f64> = Vec::with_capacity(inner - 1);
        let mut rhs: Vec<f64> = Vec::with_capacity(inner);
        for i in 0..inner {
            let j = i + 1;
            diag.push(2.0 * (h[i] + h[j]));
            if i < inner - 1 {
                upper.push(h[j]);
                lower.push(h[j]);
            }
            rhs.push(3.0 * ((ys[j + 1] - ys[j]) / h[j] - (ys[j] - ys[j - 1]) / h[i]));
        }
        let mut c_prime: Vec<f64> = vec![0.0; inner];
        let mut d_prime: Vec<f64> = vec![0.0; inner];
        c_prime[0] = if inner > 1 { upper[0] / diag[0] } else { 0.0 };
        d_prime[0] = rhs[0] / diag[0];
        for i in 1..inner {
            let denom = diag[i] - lower[i - 1] * c_prime[i - 1];
            c_prime[i] = if i < inner - 1 { upper[i] / denom } else { 0.0 };
            d_prime[i] = (rhs[i] - lower[i - 1] * d_prime[i - 1]) / denom;
        }
        let mut sigma: Vec<f64> = vec![0.0; n];
        sigma[inner] = d_prime[inner - 1];
        for i in (0..inner - 1).rev() {
            sigma[i + 1] = d_prime[i] - c_prime[i] * sigma[i + 2];
        }
        let mut coeffs: Vec<[f64; 4]> = Vec::with_capacity(m);
        for i in 0..m {
            let a = ys[i];
            let b = (ys[i + 1] - ys[i]) / h[i] - h[i] * (2.0 * sigma[i] + sigma[i + 1]) / 3.0;
            let c = sigma[i];
            let d = (sigma[i + 1] - sigma[i]) / (3.0 * h[i]);
            coeffs.push([a, b, c, d]);
        }
        let _ = h.pop();
        Ok(Self { xs, ys, coeffs })
    }
    /// Find the segment index containing `x` (binary search).
    fn find_segment(&self, x: f64) -> usize {
        let m = self.coeffs.len();
        if x <= self.xs[0] {
            return 0;
        }
        if x >= self.xs[self.xs.len() - 1] {
            return m - 1;
        }
        let idx = self.xs.partition_point(|&xi| xi <= x).saturating_sub(1);
        idx.min(m - 1)
    }
    /// Evaluate the spline at `x`.
    ///
    /// Outside the knot range the value is extrapolated from the nearest segment.
    pub fn evaluate(&self, x: f64) -> f64 {
        let i = self.find_segment(x);
        let h = x - self.xs[i];
        let [a, b, c, d] = self.coeffs[i];
        a + h * (b + h * (c + h * d))
    }
    /// Evaluate the first derivative of the spline at `x`.
    pub fn derivative(&self, x: f64) -> f64 {
        let i = self.find_segment(x);
        let h = x - self.xs[i];
        let [_a, b, c, d] = self.coeffs[i];
        b + h * (2.0 * c + h * 3.0 * d)
    }
    /// Definite integral of the spline from `x0` to `x1`.
    ///
    /// The integral is exact (analytic) over each segment.
    pub fn integral(&self, x0: f64, x1: f64) -> f64 {
        if x0 > x1 {
            return -self.integral(x1, x0);
        }
        let seg_integral = |xi: f64, xj: f64, i: usize| -> f64 {
            let h0 = xi - self.xs[i];
            let h1 = xj - self.xs[i];
            let [a, b, c, d] = self.coeffs[i];
            let antideriv =
                |h: f64| a * h + b * h * h / 2.0 + c * h * h * h / 3.0 + d * h * h * h * h / 4.0;
            antideriv(h1) - antideriv(h0)
        };
        let i0 = self.find_segment(x0);
        let i1 = self.find_segment(x1);
        if i0 == i1 {
            return seg_integral(x0, x1, i0);
        }
        let mut total = seg_integral(x0, self.xs[i0 + 1], i0);
        for i in (i0 + 1)..i1 {
            total += seg_integral(self.xs[i], self.xs[i + 1], i);
        }
        total += seg_integral(self.xs[i1], x1, i1);
        total
    }
}
/// A Non-Uniform Rational B-spline (NURBS) curve in 3-D space.
pub struct NurbsCurve {
    /// Control points.
    pub control_points: Vec<[f64; 3]>,
    /// Weight associated with each control point.
    pub weights: Vec<f64>,
    /// Polynomial degree.
    pub degree: usize,
    /// Knot vector.
    pub knots: Vec<f64>,
}
impl NurbsCurve {
    /// Construct a NURBS curve with a clamped uniform knot vector.
    pub fn new(points: Vec<[f64; 3]>, weights: Vec<f64>, degree: usize) -> Self {
        let n = points.len();
        let knots = BSplineCurve::uniform_knots(n, degree);
        Self {
            control_points: points,
            weights,
            degree,
            knots,
        }
    }
    /// Evaluate the NURBS curve at parameter `t` ∈ \[0, 1\].
    ///
    /// Computes the rational B-spline:
    /// `P(t) = Σ(w_i · N_i(t) · P_i) / Σ(w_i · N_i(t))`.
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let mut num = [0.0_f64; 3];
        let mut denom = 0.0_f64;
        for i in 0..n {
            let b = BSplineBasis::basis(i, self.degree, t, &self.knots);
            let w = self.weights[i];
            let wn = w * b;
            num[0] += wn * self.control_points[i][0];
            num[1] += wn * self.control_points[i][1];
            num[2] += wn * self.control_points[i][2];
            denom += wn;
        }
        if denom.abs() < f64::EPSILON {
            return [0.0; 3];
        }
        [num[0] / denom, num[1] / denom, num[2] / denom]
    }
}
/// Cox-de Boor B-spline basis functions.
pub struct BSplineBasis;
impl BSplineBasis {
    /// Evaluate the i-th B-spline basis function of degree `k` at parameter `t`.
    ///
    /// Uses the Cox-de Boor recursion.
    pub fn basis(i: usize, k: usize, t: f64, knots: &[f64]) -> f64 {
        if k == 0 {
            if i + 1 < knots.len() && knots[i] <= t && t < knots[i + 1] {
                return 1.0;
            }
            if i + 1 < knots.len()
                && (t - knots[i + 1]).abs() < f64::EPSILON
                && (t - knots[knots.len() - 1]).abs() < f64::EPSILON
            {
                return 1.0;
            }
            return 0.0;
        }
        let mut result = 0.0;
        if i + k < knots.len() {
            let denom = knots[i + k] - knots[i];
            if denom.abs() > f64::EPSILON {
                result += (t - knots[i]) / denom * Self::basis(i, k - 1, t, knots);
            }
        }
        if i + k + 1 < knots.len() {
            let denom = knots[i + k + 1] - knots[i + 1];
            if denom.abs() > f64::EPSILON {
                result += (knots[i + k + 1] - t) / denom * Self::basis(i + 1, k - 1, t, knots);
            }
        }
        result
    }
    /// Evaluate the derivative of the i-th B-spline basis of degree `k` at `t`.
    pub fn basis_derivative(i: usize, k: usize, t: f64, knots: &[f64]) -> f64 {
        if k == 0 {
            return 0.0;
        }
        let mut result = 0.0;
        if i + k < knots.len() {
            let denom = knots[i + k] - knots[i];
            if denom.abs() > f64::EPSILON {
                result += (k as f64) / denom * Self::basis(i, k - 1, t, knots);
            }
        }
        if i + k + 1 < knots.len() {
            let denom = knots[i + k + 1] - knots[i + 1];
            if denom.abs() > f64::EPSILON {
                result -= (k as f64) / denom * Self::basis(i + 1, k - 1, t, knots);
            }
        }
        result
    }
}
/// Fritsch-Carlson monotone cubic spline.
///
/// Fits a piecewise cubic that is monotone within each interval:
/// no artificial oscillations are introduced.
/// `xs` must be strictly increasing; `xs` and `ys` must have equal length ≥ 2.
pub struct MonotoneCubicSpline {
    /// Knot x-positions (strictly increasing).
    pub xs: Vec<f64>,
    /// Values at the knots.
    pub ys: Vec<f64>,
    /// Slopes at each knot (Fritsch-Carlson modified tangents).
    pub(super) ms: Vec<f64>,
}
impl MonotoneCubicSpline {
    /// Fit a monotone cubic spline using the Fritsch-Carlson algorithm.
    pub fn fit(xs: Vec<f64>, ys: Vec<f64>) -> Result<Self, String> {
        let n = xs.len();
        if n < 2 {
            return Err("MonotoneCubicSpline requires at least 2 knots".into());
        }
        if xs.len() != ys.len() {
            return Err("xs and ys must have the same length".into());
        }
        for i in 1..n {
            if xs[i] <= xs[i - 1] {
                return Err(format!("xs must be strictly increasing at index {i}"));
            }
        }
        let m = n - 1;
        let mut delta: Vec<f64> = (0..m)
            .map(|i| (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]))
            .collect();
        let mut tangent: Vec<f64> = vec![0.0; n];
        tangent[0] = delta[0];
        tangent[m] = delta[m - 1];
        for i in 1..m {
            tangent[i] = (delta[i - 1] + delta[i]) * 0.5;
        }
        for i in 0..m {
            if delta[i].abs() < f64::EPSILON {
                tangent[i] = 0.0;
                tangent[i + 1] = 0.0;
                continue;
            }
            let alpha = tangent[i] / delta[i];
            let beta = tangent[i + 1] / delta[i];
            let h = alpha.hypot(beta);
            if h > 3.0 {
                let tau = 3.0 / h;
                tangent[i] = tau * alpha * delta[i];
                tangent[i + 1] = tau * beta * delta[i];
            }
        }
        let _ = delta.pop();
        Ok(Self {
            xs,
            ys,
            ms: tangent,
        })
    }
    /// Evaluate the spline at `x` (clamps outside the knot range).
    pub fn evaluate(&self, x: f64) -> f64 {
        let n = self.xs.len();
        if x <= self.xs[0] {
            return self.ys[0];
        }
        if x >= self.xs[n - 1] {
            return self.ys[n - 1];
        }
        let i = self
            .xs
            .partition_point(|&xi| xi <= x)
            .saturating_sub(1)
            .min(n - 2);
        let h = self.xs[i + 1] - self.xs[i];
        let t = (x - self.xs[i]) / h;
        hermite(
            self.ys[i],
            self.ms[i] * h,
            self.ys[i + 1],
            self.ms[i + 1] * h,
            t,
        )
    }
}
/// Radial basis function (RBF) interpolation in arbitrary dimensions.
///
/// Supports the following kernels via [`RbfKernel`]:
/// * `Gaussian(eps)` — `exp(-eps² * r²)`
/// * `Multiquadric(eps)` — `sqrt(1 + eps² * r²)`
/// * `InverseMultiquadric(eps)` — `1 / sqrt(1 + eps² * r²)`
/// * `ThinPlateSpline` — `r² * ln(r)` (2D classical TPS)
#[derive(Debug, Clone, Copy)]
pub enum RbfKernel {
    /// Gaussian RBF: `exp(-eps^2 * r^2)`.
    Gaussian(f64),
    /// Multiquadric: `sqrt(1 + eps^2 * r^2)`.
    Multiquadric(f64),
    /// Inverse multiquadric: `1 / sqrt(1 + eps^2 * r^2)`.
    InverseMultiquadric(f64),
    /// Thin-plate spline: `r^2 * ln(r)` (returns 0 at r=0).
    ThinPlateSpline,
}
impl RbfKernel {
    /// Evaluate the kernel at distance `r`.
    pub fn eval(self, r: f64) -> f64 {
        match self {
            RbfKernel::Gaussian(eps) => (-(eps * eps * r * r)).exp(),
            RbfKernel::Multiquadric(eps) => (1.0 + eps * eps * r * r).sqrt(),
            RbfKernel::InverseMultiquadric(eps) => 1.0 / (1.0 + eps * eps * r * r).sqrt(),
            RbfKernel::ThinPlateSpline => {
                if r < f64::EPSILON {
                    0.0
                } else {
                    r * r * r.ln()
                }
            }
        }
    }
}
/// Radial Basis Function (RBF) interpolant in 3-D space.
///
/// Fits scattered data `(centers[i], values[i])` by solving for weights such
/// that `Σ w_j · φ(‖p − centers_j‖)` reproduces the data exactly.
pub struct RBFInterpolation {
    /// Data center points.
    pub centers: Vec<[f64; 3]>,
    /// Fitted weights.
    pub weights: Vec<f64>,
    /// Shape parameter ε.
    pub epsilon: f64,
}
impl RBFInterpolation {
    /// Gaussian RBF: `exp(−(ε r)²)`.
    pub fn gaussian(r: f64, eps: f64) -> f64 {
        let er = eps * r;
        (-er * er).exp()
    }
    /// Multiquadric RBF: `sqrt(1 + (ε r)²)`.
    pub fn multiquadric(r: f64, eps: f64) -> f64 {
        let er = eps * r;
        (1.0 + er * er).sqrt()
    }
    /// Fit the RBF interpolant to the given data points.
    ///
    /// Uses Gaussian RBF and solves the linear system with partial-pivot
    /// Gaussian elimination.
    ///
    /// Returns an error if `points` and `values` have different lengths or
    /// if the interpolation matrix is (near-)singular.
    pub fn fit(points: &[[f64; 3]], values: &[f64], epsilon: f64) -> Result<Self, String> {
        let n = points.len();
        if n == 0 {
            return Err("RBFInterpolation requires at least one point".into());
        }
        if values.len() != n {
            return Err(format!(
                "points ({}) and values ({}) must have equal length",
                n,
                values.len()
            ));
        }
        let mut mat: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| Self::gaussian(dist3(points[i], points[j]), epsilon))
                    .collect()
            })
            .collect();
        let mut rhs: Vec<f64> = values.to_vec();
        for col in 0..n {
            let pivot = (col..n).max_by(|&a, &b| {
                mat[a][col]
                    .abs()
                    .partial_cmp(&mat[b][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let pivot = pivot.ok_or("RBF matrix is singular")?;
            if mat[pivot][col].abs() < 1e-14 {
                return Err("RBF matrix is singular or near-singular".into());
            }
            mat.swap(col, pivot);
            rhs.swap(col, pivot);
            let diag = mat[col][col];
            for j in col..n {
                mat[col][j] /= diag;
            }
            rhs[col] /= diag;
            for row in 0..n {
                if row == col {
                    continue;
                }
                let factor = mat[row][col];
                for j in col..n {
                    let v = mat[col][j] * factor;
                    mat[row][j] -= v;
                }
                let rv = rhs[col] * factor;
                rhs[row] -= rv;
            }
        }
        Ok(Self {
            centers: points.to_vec(),
            weights: rhs,
            epsilon,
        })
    }
    /// Evaluate the RBF interpolant at point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        self.centers
            .iter()
            .zip(self.weights.iter())
            .map(|(c, w)| w * Self::gaussian(dist3(*c, p), self.epsilon))
            .sum()
    }
}
/// Akima spline — a piecewise cubic with locally computed slopes that avoids
/// the Runge oscillation effect by using a weighted average of secant slopes.
///
/// Knots must be given in strictly increasing x-order with at least 2 points.
/// Fewer than 5 knots use simpler fallback slope estimation.
pub struct AkimaSpline {
    /// Knot x-positions.
    pub xs: Vec<f64>,
    /// Values at knots.
    pub ys: Vec<f64>,
    /// Slopes at each knot.
    pub(super) ms: Vec<f64>,
}
impl AkimaSpline {
    /// Fit an Akima spline to the given knots.
    pub fn fit(xs: Vec<f64>, ys: Vec<f64>) -> Result<Self, String> {
        let n = xs.len();
        if n < 2 {
            return Err("AkimaSpline requires at least 2 knots".into());
        }
        if xs.len() != ys.len() {
            return Err("xs and ys must have the same length".into());
        }
        for i in 1..n {
            if xs[i] <= xs[i - 1] {
                return Err(format!("xs must be strictly increasing at index {i}"));
            }
        }
        let m_count = n - 1;
        let mut delta: Vec<f64> = (0..m_count)
            .map(|i| (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]))
            .collect();
        let d0 = 2.0 * delta[0] - delta[1];
        let d1 = 2.0 * d0 - delta[0];
        let dn1 = 2.0 * delta[m_count - 1] - delta[m_count - 2];
        let dn2 = 2.0 * dn1 - delta[m_count - 1];
        let mut ext: Vec<f64> = Vec::with_capacity(m_count + 4);
        ext.push(d1);
        ext.push(d0);
        ext.extend_from_slice(&delta);
        ext.push(dn1);
        ext.push(dn2);
        let mut slopes: Vec<f64> = Vec::with_capacity(n);
        for i in 0..n {
            let m_im2 = ext[i];
            let m_im1 = ext[i + 1];
            let m_i = ext[i + 2];
            let m_ip1 = ext[i + 3];
            let w1 = (m_ip1 - m_i).abs();
            let w2 = (m_im1 - m_im2).abs();
            let denom = w1 + w2;
            if denom < f64::EPSILON {
                slopes.push(0.5 * (m_im1 + m_i));
            } else {
                slopes.push((w1 * m_im1 + w2 * m_i) / denom);
            }
        }
        let _ = delta.pop();
        Ok(Self { xs, ys, ms: slopes })
    }
    /// Evaluate the Akima spline at `x` (clamps outside the knot range).
    pub fn evaluate(&self, x: f64) -> f64 {
        let n = self.xs.len();
        if x <= self.xs[0] {
            return self.ys[0];
        }
        if x >= self.xs[n - 1] {
            return self.ys[n - 1];
        }
        let i = self
            .xs
            .partition_point(|&xi| xi <= x)
            .saturating_sub(1)
            .min(n - 2);
        let h = self.xs[i + 1] - self.xs[i];
        let t = (x - self.xs[i]) / h;
        hermite(
            self.ys[i],
            self.ms[i] * h,
            self.ys[i + 1],
            self.ms[i + 1] * h,
            t,
        )
    }
}
/// A B-spline curve in 3-D space.
pub struct BSplineCurve {
    /// Control points.
    pub control_points: Vec<[f64; 3]>,
    /// Polynomial degree.
    pub degree: usize,
    /// Knot vector.
    pub knots: Vec<f64>,
}
impl BSplineCurve {
    /// Construct a B-spline curve with a clamped uniform knot vector.
    pub fn new(control_points: Vec<[f64; 3]>, degree: usize) -> Self {
        let n = control_points.len();
        let knots = Self::uniform_knots(n, degree);
        Self {
            control_points,
            degree,
            knots,
        }
    }
    /// Generate a clamped uniform knot vector for `n` control points of given `degree`.
    ///
    /// The vector has `n + degree + 1` entries: `degree + 1` zeros, then
    /// `n - degree - 1` evenly spaced interior knots, then `degree + 1` ones.
    pub fn uniform_knots(n: usize, degree: usize) -> Vec<f64> {
        let m = n + degree + 1;
        let mut knots = vec![0.0_f64; m];
        let interior = if n > degree + 1 { n - degree - 1 } else { 0 };
        let total_interior = interior + 2;
        for i in 0..m {
            if i <= degree {
                knots[i] = 0.0;
            } else if i >= m - degree - 1 {
                knots[i] = 1.0;
            } else {
                knots[i] = (i - degree) as f64 / (total_interior - 1 + 1) as f64;
            }
        }
        knots
    }
    /// Evaluate the B-spline curve at parameter `t` ∈ \[0, 1\].
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let mut result = [0.0_f64; 3];
        for i in 0..n {
            let b = BSplineBasis::basis(i, self.degree, t, &self.knots);
            result[0] += b * self.control_points[i][0];
            result[1] += b * self.control_points[i][1];
            result[2] += b * self.control_points[i][2];
        }
        result
    }
    /// Insert a knot at parameter `t` using Boehm's algorithm.
    ///
    /// The curve shape is preserved; a new control point is added.
    pub fn insert_knot(&mut self, t: f64) {
        let n = self.control_points.len();
        let p = self.degree;
        let knots = &self.knots;
        let k = knots
            .windows(2)
            .enumerate()
            .find(|(_, w)| w[0] <= t && t < w[1])
            .map(|(i, _)| i)
            .unwrap_or(knots.len() - 2);
        let mut new_pts: Vec<[f64; 3]> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            if i <= k - p {
                new_pts.push(self.control_points[i]);
            } else if i <= k {
                let alpha = if (knots[i + p] - knots[i]).abs() < f64::EPSILON {
                    0.0
                } else {
                    (t - knots[i]) / (knots[i + p] - knots[i])
                };
                let prev = if i == 0 {
                    [0.0; 3]
                } else {
                    self.control_points[i - 1]
                };
                let curr = self.control_points[i];
                new_pts.push([
                    lerp(prev[0], curr[0], alpha),
                    lerp(prev[1], curr[1], alpha),
                    lerp(prev[2], curr[2], alpha),
                ]);
            } else {
                new_pts.push(self.control_points[i - 1]);
            }
        }
        let mut new_knots = self.knots.clone();
        new_knots.insert(k + 1, t);
        self.control_points = new_pts;
        self.knots = new_knots;
    }
}
