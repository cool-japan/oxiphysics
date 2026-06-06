//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::f64::consts::PI;

/// Generate a swept surface by moving a profile curve along a spine curve.
///
/// The profile is defined in a local frame; this function sweeps it along
/// sample points of the spine, generating a grid of 3D points.
pub struct SweptSurface {
    /// Spine curve points.
    pub spine: Vec<[f64; 3]>,
    /// Profile curve (2D in local frame, stored as \[u, v\]).
    pub profile: Vec<[f64; 2]>,
}
impl SweptSurface {
    /// Create a new swept surface.
    pub fn new(spine: Vec<[f64; 3]>, profile: Vec<[f64; 2]>) -> Self {
        Self { spine, profile }
    }
    /// Compute the swept surface grid.
    ///
    /// At each spine point, constructs a Frenet-Serret-like frame and places
    /// the profile in the local normal plane.
    pub fn compute(&self) -> Vec<Vec<[f64; 3]>> {
        let ns = self.spine.len();
        let np = self.profile.len();
        let mut grid = vec![vec![[0.0f64; 3]; np]; ns];
        for (i, &sp) in self.spine.iter().enumerate() {
            let tangent = if i == 0 && ns > 1 {
                vec3_normalize(vec3_sub(self.spine[1], self.spine[0]))
            } else if i == ns - 1 && ns > 1 {
                vec3_normalize(vec3_sub(self.spine[ns - 1], self.spine[ns - 2]))
            } else if ns > 2 {
                vec3_normalize(vec3_sub(self.spine[i + 1], self.spine[i - 1]))
            } else {
                [0.0, 0.0, 1.0]
            };
            let up = if tangent[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let normal = vec3_normalize(vec3_cross(tangent, up));
            let binormal = vec3_cross(tangent, normal);
            for (j, &[u, v]) in self.profile.iter().enumerate() {
                grid[i][j] = vec3_add(vec3_add(sp, vec3_scale(normal, u)), vec3_scale(binormal, v));
            }
        }
        grid
    }
}
/// A NURBS surface in 3D space.
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    /// Control net in homogeneous coordinates, row-major (u direction first).
    pub control_net: Vec<Vec<[f64; 4]>>,
    /// Degree in u direction.
    pub degree_u: usize,
    /// Degree in v direction.
    pub degree_v: usize,
    /// Knot vector in u direction.
    pub knots_u: Vec<f64>,
    /// Knot vector in v direction.
    pub knots_v: Vec<f64>,
}
impl NurbsSurface {
    /// Create a NURBS surface.
    pub fn new(
        control_net: Vec<Vec<[f64; 4]>>,
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
    ) -> Self {
        Self {
            control_net,
            degree_u,
            degree_v,
            knots_u,
            knots_v,
        }
    }
    /// Evaluate the surface at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let nu = self.control_net.len();
        let nv = self.control_net[0].len();
        let span_u = find_knot_span(nu, self.degree_u, u, &self.knots_u);
        let span_v = find_knot_span(nv, self.degree_v, v, &self.knots_v);
        let bu = bspline_basis(span_u, self.degree_u, u, &self.knots_u);
        let bv = bspline_basis(span_v, self.degree_v, v, &self.knots_v);
        let mut hw = [0.0f64; 4];
        for (i, &bui) in bu.iter().enumerate().take(self.degree_u + 1) {
            for (j, &bvj) in bv.iter().enumerate().take(self.degree_v + 1) {
                let cpw = self.control_net[span_u - self.degree_u + i][span_v - self.degree_v + j];
                let b = bui * bvj;
                for k in 0..4 {
                    hw[k] += b * cpw[k];
                }
            }
        }
        if hw[3].abs() < 1e-300 {
            [hw[0], hw[1], hw[2]]
        } else {
            [hw[0] / hw[3], hw[1] / hw[3], hw[2] / hw[3]]
        }
    }
    /// Sample the surface on a `nu x nv` grid.
    pub fn sample_grid(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A NURBS (Non-Uniform Rational B-Spline) curve in 3D space.
#[derive(Debug, Clone)]
pub struct NurbsCurve {
    /// Control points in homogeneous coordinates `[wx, wy, wz, w]`.
    pub control_points_w: Vec<[f64; 4]>,
    /// Polynomial degree.
    pub degree: usize,
    /// Knot vector.
    pub knots: Vec<f64>,
}
impl NurbsCurve {
    /// Create a NURBS curve from weighted control points and a knot vector.
    pub fn new(control_points_w: Vec<[f64; 4]>, degree: usize, knots: Vec<f64>) -> Self {
        Self {
            control_points_w,
            degree,
            knots,
        }
    }
    /// Create a NURBS curve from 3D control points, weights, and a clamped knot vector.
    pub fn from_points_and_weights(
        points: Vec<[f64; 3]>,
        weights: Vec<f64>,
        degree: usize,
    ) -> Self {
        assert_eq!(points.len(), weights.len());
        let knots = uniform_clamped_knots(points.len(), degree);
        let control_points_w = points
            .iter()
            .zip(weights.iter())
            .map(|(p, &w)| [p[0] * w, p[1] * w, p[2] * w, w])
            .collect();
        Self::new(control_points_w, degree, knots)
    }
    /// Evaluate the rational curve at parameter `t`.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let n = self.control_points_w.len();
        let p = self.degree;
        let span = find_knot_span(n, p, t, &self.knots);
        let basis = bspline_basis(span, p, t, &self.knots);
        let mut hw = [0.0f64; 4];
        for (i, &bi) in basis.iter().enumerate().take(p + 1) {
            let cpw = self.control_points_w[span - p + i];
            for k in 0..4 {
                hw[k] += bi * cpw[k];
            }
        }
        if hw[3].abs() < 1e-300 {
            [hw[0], hw[1], hw[2]]
        } else {
            [hw[0] / hw[3], hw[1] / hw[3], hw[2] / hw[3]]
        }
    }
    /// Sample the NURBS curve at `num_samples` parameter values.
    pub fn sample(&self, num_samples: usize) -> Vec<[f64; 3]> {
        (0..num_samples)
            .map(|i| {
                let t = i as f64 / (num_samples - 1).max(1) as f64;
                self.eval(t)
            })
            .collect()
    }
    /// Construct a NURBS circle in the XY plane with given radius.
    ///
    /// Uses 9 control points (degree 2) with quarter-arc weights `cos(π/4)`.
    pub fn circle(radius: f64) -> Self {
        let w = (PI / 4.0).cos();
        let r = radius;
        let pts_w: Vec<[f64; 4]> = vec![
            [r, 0.0, 0.0, 1.0],
            [r * w, r * w, 0.0, w],
            [0.0, r, 0.0, 1.0],
            [-r * w, r * w, 0.0, w],
            [-r, 0.0, 0.0, 1.0],
            [-r * w, -r * w, 0.0, w],
            [0.0, -r, 0.0, 1.0],
            [r * w, -r * w, 0.0, w],
            [r, 0.0, 0.0, 1.0],
        ];
        let knots = vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ];
        Self::new(pts_w, 2, knots)
    }
}
/// Result of cubic spline interpolation: piece-wise cubic coefficients.
///
/// Each segment `[x[i\], x[i+1]]` has coefficients `(a, b, c, d)` such that
/// `S(x) = a + b*(x - x[i]) + c*(x - x[i])^2 + d*(x - x[i])^3`.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    /// x-knots (must be strictly increasing).
    pub x: Vec<f64>,
    /// Coefficients a (= y values at left knot).
    pub a: Vec<f64>,
    /// Coefficients b.
    pub b: Vec<f64>,
    /// Coefficients c.
    pub c: Vec<f64>,
    /// Coefficients d.
    pub d: Vec<f64>,
}
impl CubicSpline {
    /// Build a **natural** cubic spline through the given `(x, y)` data.
    ///
    /// Natural boundary: second derivative = 0 at both endpoints.
    pub fn natural(x: Vec<f64>, y: Vec<f64>) -> Self {
        let n = x.len();
        assert!(n >= 2, "need at least 2 points");
        assert_eq!(x.len(), y.len());
        let nm1 = n - 1;
        let mut h = vec![0.0f64; nm1];
        for i in 0..nm1 {
            h[i] = x[i + 1] - x[i];
            assert!(h[i] > 0.0, "x must be strictly increasing");
        }
        let mut alpha = vec![0.0f64; nm1];
        for i in 1..nm1 {
            alpha[i] = 3.0 / h[i] * (y[i + 1] - y[i]) - 3.0 / h[i - 1] * (y[i] - y[i - 1]);
        }
        let mut l = vec![1.0f64; n];
        let mut mu = vec![0.0f64; n];
        let mut z = vec![0.0f64; n];
        for i in 1..nm1 {
            l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        let mut c = vec![0.0f64; n];
        let mut b = vec![0.0f64; nm1];
        let mut d = vec![0.0f64; nm1];
        for j in (0..nm1).rev() {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        let a = y[..nm1].to_vec();
        CubicSpline {
            x,
            a,
            b,
            c: c[..nm1].to_vec(),
            d,
        }
    }
    /// Build a **clamped** cubic spline with given endpoint slopes `d0` and `dn`.
    pub fn clamped(x: Vec<f64>, y: Vec<f64>, d0: f64, dn: f64) -> Self {
        let n = x.len();
        assert!(n >= 2);
        assert_eq!(x.len(), y.len());
        let nm1 = n - 1;
        let mut h = vec![0.0f64; nm1];
        for i in 0..nm1 {
            h[i] = x[i + 1] - x[i];
        }
        let mut alpha = vec![0.0f64; n];
        alpha[0] = 3.0 * (y[1] - y[0]) / h[0] - 3.0 * d0;
        alpha[nm1] = 3.0 * dn - 3.0 * (y[nm1] - y[nm1 - 1]) / h[nm1 - 1];
        for i in 1..nm1 {
            alpha[i] = 3.0 / h[i] * (y[i + 1] - y[i]) - 3.0 / h[i - 1] * (y[i] - y[i - 1]);
        }
        let mut l = vec![0.0f64; n];
        let mut mu = vec![0.0f64; n];
        let mut z = vec![0.0f64; n];
        l[0] = 2.0 * h[0];
        mu[0] = 0.5;
        z[0] = alpha[0] / l[0];
        for i in 1..nm1 {
            l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[nm1] = h[nm1 - 1] * (2.0 - mu[nm1 - 1]);
        z[nm1] = (alpha[nm1] - h[nm1 - 1] * z[nm1 - 1]) / l[nm1];
        let mut c = vec![0.0f64; n];
        c[nm1] = z[nm1];
        let mut b = vec![0.0f64; nm1];
        let mut d = vec![0.0f64; nm1];
        for j in (0..nm1).rev() {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        let a = y[..nm1].to_vec();
        CubicSpline {
            x,
            a,
            b,
            c: c[..nm1].to_vec(),
            d,
        }
    }
    /// Evaluate the spline at `t`.
    pub fn eval(&self, t: f64) -> f64 {
        let n = self.x.len() - 1;
        let i = if t <= self.x[0] {
            0
        } else if t >= self.x[n] {
            n - 1
        } else {
            let mut lo = 0usize;
            let mut hi = n;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if self.x[mid] <= t {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        let dx = t - self.x[i];
        self.a[i] + self.b[i] * dx + self.c[i] * dx * dx + self.d[i] * dx * dx * dx
    }
    /// Evaluate the first derivative at `t`.
    pub fn eval_derivative(&self, t: f64) -> f64 {
        let n = self.x.len() - 1;
        let i = if t <= self.x[0] {
            0
        } else if t >= self.x[n] {
            n - 1
        } else {
            let mut lo = 0usize;
            let mut hi = n;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if self.x[mid] <= t {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        let dx = t - self.x[i];
        self.b[i] + 2.0 * self.c[i] * dx + 3.0 * self.d[i] * dx * dx
    }
    /// Evaluate the second derivative at `t`.
    pub fn eval_second_derivative(&self, t: f64) -> f64 {
        let n = self.x.len() - 1;
        let i = if t <= self.x[0] {
            0
        } else if t >= self.x[n] {
            n - 1
        } else {
            let mut lo = 0usize;
            let mut hi = n;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if self.x[mid] <= t {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        let dx = t - self.x[i];
        2.0 * self.c[i] + 6.0 * self.d[i] * dx
    }
}
/// A blending surface connecting two boundary curves.
///
/// The surface is parameterized by `(u, v)` where `u ∈ [0,1]` interpolates
/// between the two boundary curves, and `v ∈ [0,1]` walks along them.
#[derive(Debug, Clone)]
pub struct BlendingSurface {
    /// First boundary curve control points.
    pub curve0: BezierCurve,
    /// Second boundary curve control points.
    pub curve1: BezierCurve,
    /// Tangent scale factor at curve0 (for G1/G2).
    pub tangent_scale0: f64,
    /// Tangent scale factor at curve1 (for G1/G2).
    pub tangent_scale1: f64,
    /// Continuity order.
    pub continuity: ContinuityOrder,
}
impl BlendingSurface {
    /// Create a blending surface with specified continuity.
    pub fn new(
        curve0: BezierCurve,
        curve1: BezierCurve,
        tangent_scale0: f64,
        tangent_scale1: f64,
        continuity: ContinuityOrder,
    ) -> Self {
        Self {
            curve0,
            curve1,
            tangent_scale0,
            tangent_scale1,
            continuity,
        }
    }
    /// Evaluate the blending surface at `(u, v)`.
    ///
    /// Uses Hermite blending functions to achieve the requested continuity.
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let p0 = self.curve0.eval(v);
        let p1 = self.curve1.eval(v);
        match self.continuity {
            ContinuityOrder::G0 => vec3_lerp(p0, p1, u),
            ContinuityOrder::G1 => {
                let h00 = 2.0 * u * u * u - 3.0 * u * u + 1.0;
                let h10 = u * u * u - 2.0 * u * u + u;
                let h01 = -2.0 * u * u * u + 3.0 * u * u;
                let h11 = u * u * u - u * u;
                let t0 = self.curve0.tangent(v);
                let t1 = self.curve1.tangent(v);
                let mut result = [0.0f64; 3];
                for k in 0..3 {
                    result[k] = h00 * p0[k]
                        + h10 * self.tangent_scale0 * t0[k]
                        + h01 * p1[k]
                        + h11 * self.tangent_scale1 * t1[k];
                }
                result
            }
            ContinuityOrder::G2 => {
                let h00 = 1.0 - 10.0 * u * u * u + 15.0 * u * u * u * u - 6.0 * u * u * u * u * u;
                let h10 = u - 6.0 * u * u * u + 8.0 * u * u * u * u - 3.0 * u * u * u * u * u;
                let h20 =
                    0.5 * u * u - 1.5 * u * u * u + 1.5 * u * u * u * u - 0.5 * u * u * u * u * u;
                let h01 = 10.0 * u * u * u - 15.0 * u * u * u * u + 6.0 * u * u * u * u * u;
                let h11 = -4.0 * u * u * u + 7.0 * u * u * u * u - 3.0 * u * u * u * u * u;
                let h21 = 0.5 * u * u * u - u * u * u * u + 0.5 * u * u * u * u * u;
                let t0 = self.curve0.tangent(v);
                let t1 = self.curve1.tangent(v);
                let a0 = self.curve0.second_derivative(v);
                let a1 = self.curve1.second_derivative(v);
                let mut result = [0.0f64; 3];
                for k in 0..3 {
                    result[k] = h00 * p0[k]
                        + h10 * self.tangent_scale0 * t0[k]
                        + h20 * self.tangent_scale0 * self.tangent_scale0 * a0[k]
                        + h01 * p1[k]
                        + h11 * self.tangent_scale1 * t1[k]
                        + h21 * self.tangent_scale1 * self.tangent_scale1 * a1[k];
                }
                result
            }
        }
    }
    /// Sample the blending surface on a `nu x nv` grid.
    pub fn sample_grid(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A non-rational B-spline surface.
#[derive(Debug, Clone)]
pub struct BSplineSurface {
    /// Control net (nu x nv), stored as rows of v-direction points.
    pub control_net: Vec<Vec<[f64; 3]>>,
    /// Degree in u direction.
    pub degree_u: usize,
    /// Degree in v direction.
    pub degree_v: usize,
    /// Knot vector in u direction.
    pub knots_u: Vec<f64>,
    /// Knot vector in v direction.
    pub knots_v: Vec<f64>,
}
impl BSplineSurface {
    /// Create a B-spline surface with explicit knot vectors.
    pub fn new(
        control_net: Vec<Vec<[f64; 3]>>,
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
    ) -> Self {
        Self {
            control_net,
            degree_u,
            degree_v,
            knots_u,
            knots_v,
        }
    }
    /// Create a B-spline surface with uniform clamped knot vectors.
    pub fn with_clamped_knots(
        control_net: Vec<Vec<[f64; 3]>>,
        degree_u: usize,
        degree_v: usize,
    ) -> Self {
        let nu = control_net.len();
        let nv = control_net[0].len();
        let knots_u = uniform_clamped_knots(nu, degree_u);
        let knots_v = uniform_clamped_knots(nv, degree_v);
        Self::new(control_net, degree_u, degree_v, knots_u, knots_v)
    }
    /// Evaluate the surface at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let nu = self.control_net.len();
        let nv = self.control_net[0].len();
        let span_u = find_knot_span(nu, self.degree_u, u, &self.knots_u);
        let span_v = find_knot_span(nv, self.degree_v, v, &self.knots_v);
        let bu = bspline_basis(span_u, self.degree_u, u, &self.knots_u);
        let bv = bspline_basis(span_v, self.degree_v, v, &self.knots_v);
        let mut point = [0.0f64; 3];
        for (i, &bui) in bu.iter().enumerate().take(self.degree_u + 1) {
            for (j, &bvj) in bv.iter().enumerate().take(self.degree_v + 1) {
                let cp = self.control_net[span_u - self.degree_u + i][span_v - self.degree_v + j];
                let b = bui * bvj;
                point = vec3_add(point, vec3_scale(cp, b));
            }
        }
        point
    }
    /// Compute the surface normal at `(u, v)` using partial derivatives.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-6;
        let pu = vec3_sub(
            self.eval((u + eps).min(1.0), v),
            self.eval((u - eps).max(0.0), v),
        );
        let pv = vec3_sub(
            self.eval(u, (v + eps).min(1.0)),
            self.eval(u, (v - eps).max(0.0)),
        );
        vec3_normalize(vec3_cross(pu, pv))
    }
    /// Sample the surface on a `nu x nv` grid.
    pub fn sample_grid(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A tensor-product Bézier surface patch.
///
/// Defined by a `(m+1) x (n+1)` grid of control points (degree `m` in u, degree `n` in v).
#[derive(Debug, Clone)]
pub struct BezierPatch {
    /// Control point grid, row-major (u direction first), each row has `n+1` points.
    pub control_grid: Vec<Vec<[f64; 3]>>,
}
impl BezierPatch {
    /// Create a new Bézier patch.
    pub fn new(control_grid: Vec<Vec<[f64; 3]>>) -> Self {
        Self { control_grid }
    }
    /// Evaluate the patch at `(u, v)` using de Casteljau in both directions.
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let row_pts: Vec<[f64; 3]> = self
            .control_grid
            .iter()
            .map(|row| BezierCurve::new(row.clone()).eval(v))
            .collect();
        BezierCurve::new(row_pts).eval(u)
    }
    /// Sample the patch on an `nu x nv` grid.
    pub fn sample_grid(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A Bézier curve in 3D space defined by its control polygon.
#[derive(Debug, Clone)]
pub struct BezierCurve {
    /// Control points.
    pub control_points: Vec<[f64; 3]>,
}
impl BezierCurve {
    /// Create a new Bézier curve.
    pub fn new(control_points: Vec<[f64; 3]>) -> Self {
        Self { control_points }
    }
    /// Evaluate the curve at `t ∈ [0,1]` using de Casteljau's algorithm.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let mut pts = self.control_points.clone();
        let n = pts.len();
        for r in 1..n {
            for i in 0..n - r {
                pts[i] = vec3_lerp(pts[i], pts[i + 1], t);
            }
        }
        pts[0]
    }
    /// Compute the tangent (first derivative) at `t`.
    pub fn tangent(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        if n < 2 {
            return [0.0, 0.0, 0.0];
        }
        let d: Vec<[f64; 3]> = (0..n - 1)
            .map(|i| {
                vec3_scale(
                    vec3_sub(self.control_points[i + 1], self.control_points[i]),
                    (n - 1) as f64,
                )
            })
            .collect();
        BezierCurve::new(d).eval(t)
    }
    /// Evaluate the second derivative at `t`.
    pub fn second_derivative(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        if n < 3 {
            return [0.0, 0.0, 0.0];
        }
        let d1: Vec<[f64; 3]> = (0..n - 1)
            .map(|i| {
                vec3_scale(
                    vec3_sub(self.control_points[i + 1], self.control_points[i]),
                    (n - 1) as f64,
                )
            })
            .collect();
        let d2: Vec<[f64; 3]> = (0..n - 2)
            .map(|i| vec3_scale(vec3_sub(d1[i + 1], d1[i]), (n - 2) as f64))
            .collect();
        BezierCurve::new(d2).eval(t)
    }
    /// Split the Bézier curve at `t` returning two sub-curves (de Casteljau).
    pub fn split(&self, t: f64) -> (BezierCurve, BezierCurve) {
        let n = self.control_points.len();
        let mut triangle = vec![self.control_points.clone()];
        for r in 1..n {
            let prev = &triangle[r - 1];
            let row: Vec<[f64; 3]> = (0..n - r)
                .map(|i| vec3_lerp(prev[i], prev[i + 1], t))
                .collect();
            triangle.push(row);
        }
        let left: Vec<[f64; 3]> = (0..n).map(|r| triangle[r][0]).collect();
        let right: Vec<[f64; 3]> = (0..n).map(|r| triangle[n - 1 - r][r]).collect();
        (BezierCurve::new(left), BezierCurve::new(right))
    }
    /// Sample the curve at `num_samples` parameter values.
    pub fn sample(&self, num_samples: usize) -> Vec<[f64; 3]> {
        (0..num_samples)
            .map(|i| {
                let t = i as f64 / (num_samples - 1).max(1) as f64;
                self.eval(t)
            })
            .collect()
    }
    /// Compute approximate arc length using `num_segments` sub-intervals.
    pub fn arc_length(&self, num_segments: usize) -> f64 {
        let pts = self.sample(num_segments + 1);
        pts.windows(2)
            .map(|w| vec3_norm(vec3_sub(w[1], w[0])))
            .sum()
    }
    /// Elevate the degree of the Bézier curve by 1.
    pub fn degree_elevate(&self) -> BezierCurve {
        let n = self.control_points.len();
        let mut new_pts = Vec::with_capacity(n + 1);
        new_pts.push(self.control_points[0]);
        for i in 1..n {
            let alpha = i as f64 / n as f64;
            new_pts.push(vec3_add(
                vec3_scale(self.control_points[i - 1], alpha),
                vec3_scale(self.control_points[i], 1.0 - alpha),
            ));
        }
        new_pts.push(self.control_points[n - 1]);
        BezierCurve::new(new_pts)
    }
}
/// A periodic (closed) B-spline curve.
#[derive(Debug, Clone)]
pub struct PeriodicBSpline {
    /// Control points (the curve wraps back to the first point).
    pub control_points: Vec<[f64; 3]>,
    /// Polynomial degree.
    pub degree: usize,
}
impl PeriodicBSpline {
    /// Create a periodic B-spline.
    pub fn new(control_points: Vec<[f64; 3]>, degree: usize) -> Self {
        Self {
            control_points,
            degree,
        }
    }
    /// Evaluate the periodic curve at `t ∈ [0,1]`.
    ///
    /// Wraps the control points to create a uniform periodic knot vector.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let p = self.degree;
        let mut ext_cp = self.control_points.clone();
        for i in 0..p {
            ext_cp.push(self.control_points[i % n]);
        }
        let m = ext_cp.len() + p + 1;
        let knots: Vec<f64> = (0..m).map(|i| i as f64 / (m - 1) as f64).collect();
        let num_ctrl = ext_cp.len();
        let t_min = knots[p];
        let t_max = knots[num_ctrl];
        let t_inner = t_min + t * (t_max - t_min);
        let span = find_knot_span(num_ctrl, p, t_inner.min(t_max - 1e-12), &knots);
        let basis = bspline_basis(span, p, t_inner.min(t_max - 1e-12), &knots);
        let mut point = [0.0f64; 3];
        for i in 0..=p {
            let cp = ext_cp[span - p + i];
            point = vec3_add(point, vec3_scale(cp, basis[i]));
        }
        point
    }
}
/// A B-spline curve in 3D space defined by control points, degree, and knot vector.
#[derive(Debug, Clone)]
pub struct BSplineCurve {
    /// Control points.
    pub control_points: Vec<[f64; 3]>,
    /// Polynomial degree.
    pub degree: usize,
    /// Knot vector.
    pub knots: Vec<f64>,
}
impl BSplineCurve {
    /// Create a new B-spline curve.
    ///
    /// # Panics
    /// Panics if the knot vector length is not `control_points.len() + degree + 1`.
    pub fn new(control_points: Vec<[f64; 3]>, degree: usize, knots: Vec<f64>) -> Self {
        assert_eq!(
            knots.len(),
            control_points.len() + degree + 1,
            "knot vector length must equal n + p + 2"
        );
        Self {
            control_points,
            degree,
            knots,
        }
    }
    /// Create a B-spline curve with a uniform clamped knot vector.
    pub fn with_clamped_knots(control_points: Vec<[f64; 3]>, degree: usize) -> Self {
        let knots = uniform_clamped_knots(control_points.len(), degree);
        Self::new(control_points, degree, knots)
    }
    /// Evaluate the curve at parameter `t` in `[0, 1]`.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let p = self.degree;
        let span = find_knot_span(n, p, t, &self.knots);
        let basis = bspline_basis(span, p, t, &self.knots);
        let mut point = [0.0f64; 3];
        for (i, &bi) in basis.iter().enumerate().take(p + 1) {
            let cp = self.control_points[span - p + i];
            point = vec3_add(point, vec3_scale(cp, bi));
        }
        point
    }
    /// Compute the first derivative (tangent) at parameter `t`.
    pub fn derivative(&self, t: f64) -> [f64; 3] {
        let p = self.degree;
        if p == 0 {
            return [0.0, 0.0, 0.0];
        }
        let n = self.control_points.len();
        let mut d_cp = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let denom = self.knots[i + p + 1] - self.knots[i + 1];
            if denom.abs() < 1e-300 {
                d_cp.push([0.0, 0.0, 0.0]);
            } else {
                let diff = vec3_sub(self.control_points[i + 1], self.control_points[i]);
                d_cp.push(vec3_scale(diff, p as f64 / denom));
            }
        }
        let d_knots = self.knots[1..self.knots.len() - 1].to_vec();
        let d_curve = BSplineCurve::new(d_cp, p - 1, d_knots);
        d_curve.eval(t)
    }
    /// Sample the curve at `num_samples` equally-spaced parameter values.
    pub fn sample(&self, num_samples: usize) -> Vec<[f64; 3]> {
        (0..num_samples)
            .map(|i| {
                let t = i as f64 / (num_samples - 1).max(1) as f64;
                self.eval(t)
            })
            .collect()
    }
    /// Compute approximate arc length using `num_segments` samples.
    pub fn arc_length(&self, num_segments: usize) -> f64 {
        let pts = self.sample(num_segments + 1);
        pts.windows(2)
            .map(|w| vec3_norm(vec3_sub(w[1], w[0])))
            .sum()
    }
}
/// Continuity order for blending surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuityOrder {
    /// G0: positional continuity.
    G0,
    /// G1: tangent continuity.
    G1,
    /// G2: curvature continuity.
    G2,
}
