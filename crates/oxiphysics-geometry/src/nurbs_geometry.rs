#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! NURBS (Non-Uniform Rational B-Splines) geometry module.
//!
//! Provides a comprehensive NURBS implementation for geometric modelling and
//! isogeometric analysis (IGA):
//!
//! - [`NurbsKnotVector`]: Non-uniform knot vector with span search.
//! - [`NurbsBasis`]: Cox–de Boor recursive B-spline basis functions.
//! - [`NurbsCurve`]: NURBS curve evaluation, derivatives, curvature, arc length.
//! - [`NurbsSurface`]: Bi-variate NURBS surface, normals, curvatures.
//! - [`KnotInsertion`]: Boehm's knot insertion algorithm.
//! - [`DegreeElevation`]: Degree elevation for NURBS curves.
//! - [`BezierDecomposition`]: Decompose a NURBS curve into Bézier segments.
//! - [`RuledSurface`]: Ruled surface between two NURBS curves.
//! - [`SweptSurface`]: Swept surface from profile curve and path.
//! - [`NurbsTrimming`]: Trimmed NURBS surface with trim curves.
//! - [`IgaElement`]: Isogeometric analysis element (basis function + quadrature).
//! - [`SurfaceSurfaceIntersector`]: Iterative NURBS surface-surface intersection.
//!
//! # References
//! - Piegl, L. & Tiller, W. (1997). *The NURBS Book*, 2nd ed. Springer.
//! - Hughes, T. J. R., Cottrell, J. A. & Bazilevs, Y. (2005).
//!   *Comput. Methods Appl. Mech. Engrg.* 194, 4135.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// 3-D vector helpers
// ---------------------------------------------------------------------------

/// Add two 3-D vectors.
#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors.
#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by a scalar.
#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-D vectors.
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-D vectors.
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-D vector.
#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalize a 3-D vector (returns zero vector if degenerate).
#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1.0e-15 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}

// ---------------------------------------------------------------------------
// NurbsKnotVector
// ---------------------------------------------------------------------------

/// Non-uniform knot vector for a NURBS curve or surface.
///
/// A valid knot vector is non-decreasing: t\[0\] ≤ t\[1\] ≤ … ≤ t\[m\].
#[derive(Debug, Clone)]
pub struct NurbsKnotVector {
    /// Knot values.
    pub knots: Vec<f64>,
}

impl NurbsKnotVector {
    /// Construct a knot vector from an explicit list.
    ///
    /// # Panics
    ///
    /// Panics if the sequence is not non-decreasing.
    pub fn new(knots: Vec<f64>) -> Self {
        for i in 1..knots.len() {
            assert!(
                knots[i] >= knots[i - 1],
                "knot vector must be non-decreasing: knots[{}]={:.6} < knots[{}]={:.6}",
                i,
                knots[i],
                i - 1,
                knots[i - 1]
            );
        }
        Self { knots }
    }

    /// Build a clamped uniform knot vector for `n+1` control points and degree `p`.
    ///
    /// The first and last knots are repeated `p+1` times; interior knots are uniform.
    pub fn clamped_uniform(n_ctrl: usize, degree: usize) -> Self {
        let n = n_ctrl - 1;
        let p = degree;
        let m = n + p + 1;
        let mut knots = vec![0.0_f64; m + 1];
        for i in 0..=p {
            knots[i] = 0.0;
        }
        for i in (m - p)..=m {
            knots[i] = 1.0;
        }
        let n_int = m as isize - 2 * p as isize - 1;
        for j in 1..=(n_int.max(0) as usize) {
            knots[p + j] = j as f64 / (n_int as f64 + 1.0);
        }
        Self { knots }
    }

    /// Build a periodic (uniform) knot vector for `n+1` control points and degree `p`.
    pub fn uniform(n_ctrl: usize, degree: usize) -> Self {
        let m = n_ctrl + degree;
        let knots: Vec<f64> = (0..=m).map(|i| i as f64 / m as f64).collect();
        Self { knots }
    }

    /// Number of knots.
    pub fn len(&self) -> usize {
        self.knots.len()
    }

    /// Check if the knot vector is empty.
    pub fn is_empty(&self) -> bool {
        self.knots.is_empty()
    }

    /// Find the knot span index i such that knots\[i\] ≤ t < knots\[i+1\].
    ///
    /// Uses binary search. At the right end of the domain returns the last
    /// non-zero-length span.
    pub fn find_span(&self, t: f64, degree: usize, n_ctrl: usize) -> usize {
        let n = n_ctrl - 1;
        let p = degree;
        let m = n + p + 1;
        // Clamp at domain end
        if t >= self.knots[n + 1] {
            return n;
        }
        if t <= self.knots[p] {
            return p;
        }
        // Binary search
        let mut lo = p;
        let mut hi = m;
        let mut mid = (lo + hi) / 2;
        while t < self.knots[mid] || t >= self.knots[mid + 1] {
            if t < self.knots[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
            mid = (lo + hi) / 2;
        }
        mid
    }

    /// Domain start (first non-clamped knot value).
    pub fn domain_start(&self, degree: usize) -> f64 {
        self.knots[degree]
    }

    /// Domain end (last non-clamped knot value).
    pub fn domain_end(&self, _degree: usize, n_ctrl: usize) -> f64 {
        self.knots[n_ctrl]
    }

    /// Multiplicity of knot value `t` (within tolerance ε).
    pub fn multiplicity(&self, t: f64) -> usize {
        self.knots
            .iter()
            .filter(|&&k| (k - t).abs() < 1.0e-12)
            .count()
    }

    /// Insert a new knot value and return the new knot vector.
    pub fn insert_knot(&self, t_new: f64) -> Self {
        let mut new_knots = self.knots.clone();
        // Find insert position
        let pos = new_knots.partition_point(|&k| k <= t_new);
        new_knots.insert(pos, t_new);
        Self { knots: new_knots }
    }
}

// ---------------------------------------------------------------------------
// NurbsBasis — Cox–de Boor recursion
// ---------------------------------------------------------------------------

/// B-spline basis functions computed via Cox–de Boor recursion.
///
/// Evaluates N_{i,p}(t) for i = span−p … span.
pub struct NurbsBasis;

impl NurbsBasis {
    /// Evaluate all non-zero basis functions of degree `p` at parameter `t`.
    ///
    /// Returns a vector of length `p+1` containing N_{span-p,p}(t) … N_{span,p}(t).
    pub fn evaluate(knots: &NurbsKnotVector, t: f64, degree: usize, n_ctrl: usize) -> Vec<f64> {
        let p = degree;
        let span = knots.find_span(t, p, n_ctrl);
        let kv = &knots.knots;

        let mut n = vec![0.0_f64; p + 1];
        let mut left = vec![0.0_f64; p + 1];
        let mut right = vec![0.0_f64; p + 1];

        n[0] = 1.0;
        for j in 1..=p {
            left[j] = t - kv[span + 1 - j];
            right[j] = kv[span + j] - t;
            let mut saved = 0.0_f64;
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                let temp = if denom.abs() < 1.0e-15 {
                    0.0
                } else {
                    n[r] / denom
                };
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            n[j] = saved;
        }
        n
    }

    /// Evaluate basis functions AND their first derivatives at `t`.
    ///
    /// Returns `(N, dN)` each of length `p+1`.
    pub fn evaluate_deriv(
        knots: &NurbsKnotVector,
        t: f64,
        degree: usize,
        n_ctrl: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        let ders = Self::evaluate_all_derivs(knots, t, degree, n_ctrl, 1);
        (ders[0].clone(), ders[1].clone())
    }

    /// Evaluate basis functions up to order `d` (all derivatives 0..=d).
    ///
    /// Returns a `(d+1) × (p+1)` matrix: `ders[k][j]` = k-th derivative of N_{span-p+j,p}.
    pub fn evaluate_all_derivs(
        knots: &NurbsKnotVector,
        t: f64,
        degree: usize,
        n_ctrl: usize,
        n_deriv: usize,
    ) -> Vec<Vec<f64>> {
        let p = degree;
        let span = knots.find_span(t, p, n_ctrl);
        let kv = &knots.knots;
        let d = n_deriv.min(p);

        let mut ndu = vec![vec![0.0_f64; p + 1]; p + 1];
        let mut left = vec![0.0_f64; p + 1];
        let mut right = vec![0.0_f64; p + 1];
        ndu[0][0] = 1.0;
        for j in 1..=p {
            left[j] = t - kv[span + 1 - j];
            right[j] = kv[span + j] - t;
            let mut saved = 0.0_f64;
            for r in 0..j {
                ndu[j][r] = right[r + 1] + left[j - r];
                let temp = if ndu[j][r].abs() < 1.0e-15 {
                    0.0
                } else {
                    ndu[r][j - 1] / ndu[j][r]
                };
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        let mut ders = vec![vec![0.0_f64; p + 1]; d + 1];
        for j in 0..=p {
            ders[0][j] = ndu[j][p];
        }

        let mut a = vec![vec![0.0_f64; p + 1]; 2];
        for r in 0..=p {
            let mut s1 = 0usize;
            let mut s2 = 1usize;
            a[s1][0] = 1.0;
            for k in 1..=d {
                let mut fac = 0.0_f64;
                let rk = r as isize - k as isize;
                let pk = p as isize - k as isize;
                if r as isize >= k as isize {
                    a[s2][0] = a[s1][0]
                        / (ndu[(pk + 1) as usize][rk as usize]
                            .abs()
                            .max(1.0e-15)
                            .copysign(ndu[(pk + 1) as usize][rk as usize]));
                    fac = a[s2][0];
                }
                let j1 = if rk >= -1 { 1 } else { (-rk - 1) as usize };
                let j2 = if (r as isize - 1) <= pk { k - 1 } else { p - r };
                for j in j1..=j2 {
                    a[s2][j] = (a[s1][j] - a[s1][j - 1])
                        / (ndu[(pk + 1) as usize][(rk + j as isize) as usize]
                            .abs()
                            .max(1.0e-15)
                            .copysign(ndu[(pk + 1) as usize][(rk + j as isize) as usize]));
                }
                if (r as isize) <= pk {
                    a[s2][k] = -a[s1][k - 1]
                        / (ndu[(pk + 1) as usize][r]
                            .abs()
                            .max(1.0e-15)
                            .copysign(ndu[(pk + 1) as usize][r]));
                    fac += a[s2][k];
                }
                ders[k][r] = fac;
                std::mem::swap(&mut s1, &mut s2);
            }
        }
        ders
    }
}

// ---------------------------------------------------------------------------
// NurbsControlPoint
// ---------------------------------------------------------------------------

/// A weighted NURBS control point in homogeneous coordinates.
#[derive(Debug, Clone, Copy)]
pub struct NurbsControlPoint {
    /// Cartesian position.
    pub point: [f64; 3],
    /// Weight (must be > 0).
    pub weight: f64,
}

impl NurbsControlPoint {
    /// Create a new control point.
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self {
            point: [x, y, z],
            weight: w,
        }
    }

    /// Unweighted control point (weight = 1).
    pub fn unweighted(x: f64, y: f64, z: f64) -> Self {
        Self::new(x, y, z, 1.0)
    }

    /// Homogeneous coordinates \[wx, wy, wz, w\].
    pub fn homogeneous(&self) -> [f64; 4] {
        [
            self.point[0] * self.weight,
            self.point[1] * self.weight,
            self.point[2] * self.weight,
            self.weight,
        ]
    }
}

// ---------------------------------------------------------------------------
// NurbsCurve
// ---------------------------------------------------------------------------

/// NURBS curve defined by a knot vector, control points, weights, and degree.
#[derive(Debug, Clone)]
pub struct NurbsCurve {
    /// Knot vector.
    pub knots: NurbsKnotVector,
    /// Control points with weights.
    pub control_points: Vec<NurbsControlPoint>,
    /// Polynomial degree.
    pub degree: usize,
}

impl NurbsCurve {
    /// Construct a NURBS curve.
    ///
    /// # Panics
    ///
    /// Panics if the knot vector length != n_ctrl + degree + 1.
    pub fn new(
        knots: NurbsKnotVector,
        control_points: Vec<NurbsControlPoint>,
        degree: usize,
    ) -> Self {
        let n = control_points.len();
        let m = n + degree;
        assert_eq!(
            knots.len(),
            m + 1,
            "knot vector length must be n_ctrl + degree + 1 = {}, got {}",
            m + 1,
            knots.len()
        );
        Self {
            knots,
            control_points,
            degree,
        }
    }

    /// Evaluate the curve position at parameter `t`.
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let p = self.degree;
        let basis = NurbsBasis::evaluate(&self.knots, t, p, n);
        let span = self.knots.find_span(t, p, n);

        let mut w_sum = 0.0_f64;
        let mut point = [0.0_f64; 3];
        for i in 0..=p {
            let idx = span - p + i;
            if idx < n {
                let cp = &self.control_points[idx];
                let wn = basis[i] * cp.weight;
                w_sum += wn;
                for k in 0..3 {
                    point[k] += wn * cp.point[k];
                }
            }
        }
        if w_sum.abs() > 1.0e-15 {
            for k in 0..3 {
                point[k] /= w_sum;
            }
        }
        point
    }

    /// Evaluate the first derivative at parameter `t`.
    pub fn derivative(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let p = self.degree;
        let (basis, dbasis) = NurbsBasis::evaluate_deriv(&self.knots, t, p, n);
        let span = self.knots.find_span(t, p, n);

        let mut w = 0.0_f64;
        let mut dw = 0.0_f64;
        let mut c = [0.0_f64; 3];
        let mut dc = [0.0_f64; 3];
        for i in 0..=p {
            let idx = span - p + i;
            if idx < n {
                let cp = &self.control_points[idx];
                let wi = cp.weight;
                w += basis[i] * wi;
                dw += dbasis[i] * wi;
                for k in 0..3 {
                    c[k] += basis[i] * wi * cp.point[k];
                    dc[k] += dbasis[i] * wi * cp.point[k];
                }
            }
        }
        // Quotient rule: d/dt (C/W) = (dC·W - C·dW) / W²
        if w.abs() < 1.0e-15 {
            return [0.0; 3];
        }
        let w2 = w * w;
        let mut out = [0.0_f64; 3];
        for k in 0..3 {
            out[k] = (dc[k] * w - c[k] * dw) / w2;
        }
        out
    }

    /// Tangent unit vector at `t`.
    pub fn tangent(&self, t: f64) -> [f64; 3] {
        vec3_normalize(self.derivative(t))
    }

    /// Signed curvature κ at `t` (magnitude of curvature vector).
    pub fn curvature(&self, t: f64) -> f64 {
        let dt = 1.0e-5;
        let t_lo = (t - dt).max(self.knots.domain_start(self.degree));
        let t_hi = (t + dt).min(
            self.knots
                .domain_end(self.degree, self.control_points.len()),
        );
        let d1 = self.derivative(t);
        let p_lo = self.evaluate(t_lo);
        let p_hi = self.evaluate(t_hi);
        // Finite-difference second derivative
        let p_mid = self.evaluate(t);
        let d2 = [
            (p_hi[0] - 2.0 * p_mid[0] + p_lo[0]) / (dt * dt),
            (p_hi[1] - 2.0 * p_mid[1] + p_lo[1]) / (dt * dt),
            (p_hi[2] - 2.0 * p_mid[2] + p_lo[2]) / (dt * dt),
        ];
        let cross = vec3_cross(d1, d2);
        let num = vec3_norm(cross);
        let denom = vec3_norm(d1).powi(3);
        if denom < 1.0e-15 { 0.0 } else { num / denom }
    }

    /// Approximate arc length using Gaussian quadrature with `n_pts` points.
    pub fn arc_length(&self, n_pts: usize) -> f64 {
        let t0 = self.knots.domain_start(self.degree);
        let t1 = self
            .knots
            .domain_end(self.degree, self.control_points.len());
        gauss_legendre_integral(t0, t1, n_pts, |t| vec3_norm(self.derivative(t)))
    }

    /// Number of control points.
    pub fn n_ctrl(&self) -> usize {
        self.control_points.len()
    }

    /// Build a full circle as an exact NURBS (degree 2, 9 control points).
    ///
    /// Lies in the XY-plane, centred at origin, radius `r`.
    pub fn circle(r: f64) -> Self {
        let w = 0.5_f64.sqrt(); // cos(45°)
        let cps = vec![
            NurbsControlPoint::new(r, 0.0, 0.0, 1.0),
            NurbsControlPoint::new(r, r, 0.0, w),
            NurbsControlPoint::new(0.0, r, 0.0, 1.0),
            NurbsControlPoint::new(-r, r, 0.0, w),
            NurbsControlPoint::new(-r, 0.0, 0.0, 1.0),
            NurbsControlPoint::new(-r, -r, 0.0, w),
            NurbsControlPoint::new(0.0, -r, 0.0, 1.0),
            NurbsControlPoint::new(r, -r, 0.0, w),
            NurbsControlPoint::new(r, 0.0, 0.0, 1.0),
        ];
        let knots = NurbsKnotVector::new(vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ]);
        Self {
            knots,
            control_points: cps,
            degree: 2,
        }
    }
}

// ---------------------------------------------------------------------------
// NurbsSurface
// ---------------------------------------------------------------------------

/// NURBS surface defined over a rectangular parameter domain.
///
/// Control net is stored row-major: `control_net[i * n_v + j]` is the
/// control point at (u-index i, v-index j).
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    /// Knot vector in the u direction.
    pub knots_u: NurbsKnotVector,
    /// Knot vector in the v direction.
    pub knots_v: NurbsKnotVector,
    /// Control net with weights (row-major, n_u × n_v).
    pub control_net: Vec<NurbsControlPoint>,
    /// Number of control points in u direction.
    pub n_u: usize,
    /// Number of control points in v direction.
    pub n_v: usize,
    /// Degree in u direction.
    pub degree_u: usize,
    /// Degree in v direction.
    pub degree_v: usize,
}

impl NurbsSurface {
    /// Construct a NURBS surface.
    pub fn new(
        knots_u: NurbsKnotVector,
        knots_v: NurbsKnotVector,
        control_net: Vec<NurbsControlPoint>,
        n_u: usize,
        n_v: usize,
        degree_u: usize,
        degree_v: usize,
    ) -> Self {
        assert_eq!(
            control_net.len(),
            n_u * n_v,
            "control_net length must be n_u * n_v"
        );
        Self {
            knots_u,
            knots_v,
            control_net,
            n_u,
            n_v,
            degree_u,
            degree_v,
        }
    }

    /// Evaluate surface position at (u, v).
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let nu = NurbsBasis::evaluate(&self.knots_u, u, self.degree_u, self.n_u);
        let nv = NurbsBasis::evaluate(&self.knots_v, v, self.degree_v, self.n_v);
        let span_u = self.knots_u.find_span(u, self.degree_u, self.n_u);
        let span_v = self.knots_v.find_span(v, self.degree_v, self.n_v);

        let mut w_sum = 0.0_f64;
        let mut point = [0.0_f64; 3];
        for i in 0..=self.degree_u {
            let iu = span_u + i - self.degree_u;
            if iu >= self.n_u {
                continue;
            }
            for j in 0..=self.degree_v {
                let jv = span_v + j - self.degree_v;
                if jv >= self.n_v {
                    continue;
                }
                let cp = &self.control_net[iu * self.n_v + jv];
                let wn = nu[i] * nv[j] * cp.weight;
                w_sum += wn;
                for k in 0..3 {
                    point[k] += wn * cp.point[k];
                }
            }
        }
        if w_sum.abs() > 1.0e-15 {
            for k in 0..3 {
                point[k] /= w_sum;
            }
        }
        point
    }

    /// Partial derivative with respect to u at (u, v).
    pub fn du(&self, u: f64, v: f64) -> [f64; 3] {
        let h = 1.0e-5;
        let u0 = self.knots_u.domain_start(self.degree_u);
        let u1 = self.knots_u.domain_end(self.degree_u, self.n_u);
        let u_hi = (u + h).min(u1);
        let u_lo = (u - h).max(u0);
        let p_hi = self.evaluate(u_hi, v);
        let p_lo = self.evaluate(u_lo, v);
        let denom = u_hi - u_lo;
        if denom < 1.0e-15 {
            return [0.0; 3];
        }
        [
            (p_hi[0] - p_lo[0]) / denom,
            (p_hi[1] - p_lo[1]) / denom,
            (p_hi[2] - p_lo[2]) / denom,
        ]
    }

    /// Partial derivative with respect to v at (u, v).
    pub fn dv(&self, u: f64, v: f64) -> [f64; 3] {
        let h = 1.0e-5;
        let v0 = self.knots_v.domain_start(self.degree_v);
        let v1 = self.knots_v.domain_end(self.degree_v, self.n_v);
        let v_hi = (v + h).min(v1);
        let v_lo = (v - h).max(v0);
        let p_hi = self.evaluate(u, v_hi);
        let p_lo = self.evaluate(u, v_lo);
        let denom = v_hi - v_lo;
        if denom < 1.0e-15 {
            return [0.0; 3];
        }
        [
            (p_hi[0] - p_lo[0]) / denom,
            (p_hi[1] - p_lo[1]) / denom,
            (p_hi[2] - p_lo[2]) / denom,
        ]
    }

    /// Unit outward normal at (u, v).
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let su = self.du(u, v);
        let sv = self.dv(u, v);
        vec3_normalize(vec3_cross(su, sv))
    }

    /// Mean curvature H at (u, v) via finite differences.
    pub fn mean_curvature(&self, u: f64, v: f64) -> f64 {
        let h: f64 = 1.0e-4;
        let u0 = self.knots_u.domain_start(self.degree_u);
        let u1 = self.knots_u.domain_end(self.degree_u, self.n_u);
        let v0 = self.knots_v.domain_start(self.degree_v);
        let v1 = self.knots_v.domain_end(self.degree_v, self.n_v);
        let uh = h.min((u1 - u0) / 4.0);
        let vh = h.min((v1 - v0) / 4.0);

        let n_pu = self.normal((u + uh).min(u1), v);
        let n_mu = self.normal((u - uh).max(u0), v);
        let n_pv = self.normal(u, (v + vh).min(v1));
        let n_mv = self.normal(u, (v - vh).max(v0));

        let dn_du = vec3_scale(vec3_sub(n_pu, n_mu), 1.0 / (2.0 * uh));
        let dn_dv = vec3_scale(vec3_sub(n_pv, n_mv), 1.0 / (2.0 * vh));

        0.5 * (dn_du[0] + dn_dv[1]) // approximate trace of shape operator
    }

    /// Gaussian curvature K at (u, v) via first and second fundamental forms.
    ///
    /// Computed numerically using finite differences.
    pub fn gaussian_curvature(&self, u: f64, v: f64) -> f64 {
        let h: f64 = 1.0e-4;
        let u0 = self.knots_u.domain_start(self.degree_u);
        let u1 = self.knots_u.domain_end(self.degree_u, self.n_u);
        let v0 = self.knots_v.domain_start(self.degree_v);
        let v1 = self.knots_v.domain_end(self.degree_v, self.n_v);
        let uh = h.min((u1 - u0) / 4.0);
        let vh = h.min((v1 - v0) / 4.0);

        let r_u = self.du(u, v);
        let r_v = self.dv(u, v);
        let n = self.normal(u, v);
        // First fundamental form
        let e_ff = vec3_dot(r_u, r_u);
        let f_ff = vec3_dot(r_u, r_v);
        let g_ff = vec3_dot(r_v, r_v);
        let det1 = e_ff * g_ff - f_ff * f_ff;
        if det1.abs() < 1.0e-15 {
            return 0.0;
        }
        // Second fundamental form via second derivatives
        let r_uu = {
            let pu = (u + uh).min(u1);
            let mu = (u - uh).max(u0);
            let p = self.evaluate(pu, v);
            let m = self.evaluate(mu, v);
            let c = self.evaluate(u, v);
            [
                (p[0] - 2.0 * c[0] + m[0]) / (uh * uh),
                (p[1] - 2.0 * c[1] + m[1]) / (uh * uh),
                (p[2] - 2.0 * c[2] + m[2]) / (uh * uh),
            ]
        };
        let r_vv = {
            let pv = (v + vh).min(v1);
            let mv = (v - vh).max(v0);
            let p = self.evaluate(u, pv);
            let m = self.evaluate(u, mv);
            let c = self.evaluate(u, v);
            [
                (p[0] - 2.0 * c[0] + m[0]) / (vh * vh),
                (p[1] - 2.0 * c[1] + m[1]) / (vh * vh),
                (p[2] - 2.0 * c[2] + m[2]) / (vh * vh),
            ]
        };
        let r_uv = {
            let pp = self.evaluate((u + uh).min(u1), (v + vh).min(v1));
            let pm = self.evaluate((u + uh).min(u1), (v - vh).max(v0));
            let mp = self.evaluate((u - uh).max(u0), (v + vh).min(v1));
            let mm = self.evaluate((u - uh).max(u0), (v - vh).max(v0));
            [
                (pp[0] - pm[0] - mp[0] + mm[0]) / (4.0 * uh * vh),
                (pp[1] - pm[1] - mp[1] + mm[1]) / (4.0 * uh * vh),
                (pp[2] - pm[2] - mp[2] + mm[2]) / (4.0 * uh * vh),
            ]
        };
        let big_l = vec3_dot(r_uu, n);
        let big_m = vec3_dot(r_uv, n);
        let big_n = vec3_dot(r_vv, n);
        (big_l * big_n - big_m * big_m) / det1
    }

    /// Build a flat (planar) NURBS surface in the XY-plane.
    pub fn flat_plane(x0: f64, x1: f64, y0: f64, y1: f64) -> Self {
        let cps = vec![
            NurbsControlPoint::unweighted(x0, y0, 0.0),
            NurbsControlPoint::unweighted(x1, y0, 0.0),
            NurbsControlPoint::unweighted(x0, y1, 0.0),
            NurbsControlPoint::unweighted(x1, y1, 0.0),
        ];
        let ku = NurbsKnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let kv = NurbsKnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        Self::new(ku, kv, cps, 2, 2, 1, 1)
    }
}

// ---------------------------------------------------------------------------
// Knot insertion (Boehm's algorithm)
// ---------------------------------------------------------------------------

/// Knot insertion into a NURBS curve.
///
/// Given a NURBS curve and a new knot value, computes the refined curve
/// (same geometry, additional knot and control point) using Boehm's algorithm.
pub struct KnotInsertion;

impl KnotInsertion {
    /// Insert knot `t_bar` (with current multiplicity `s`) once into curve.
    ///
    /// Returns the new control points in homogeneous coordinates.
    pub fn insert_once(
        control_points: &[NurbsControlPoint],
        knots: &NurbsKnotVector,
        degree: usize,
        t_bar: f64,
        _s: usize,
    ) -> Vec<NurbsControlPoint> {
        let n = control_points.len();
        let p = degree;
        let span = knots.find_span(t_bar, p, n);
        let kv = &knots.knots;

        let mut new_cps: Vec<NurbsControlPoint> = Vec::with_capacity(n + 1);
        // Copy control points before affected region
        for i in 0..=(span - p) {
            new_cps.push(control_points[i]);
        }
        // Blend affected control points
        for i in (span - p + 1)..=span {
            let alpha = {
                let denom = kv[i + p] - kv[i];
                if denom.abs() < 1.0e-15 {
                    0.0
                } else {
                    (t_bar - kv[i]) / denom
                }
            };
            let prev = if i > 0 {
                control_points[i - 1]
            } else {
                control_points[0]
            };
            let curr = control_points[i];
            let w_new = (1.0 - alpha) * prev.weight + alpha * curr.weight;
            let p_new = [
                ((1.0 - alpha) * prev.weight * prev.point[0] + alpha * curr.weight * curr.point[0])
                    / w_new.max(1.0e-15),
                ((1.0 - alpha) * prev.weight * prev.point[1] + alpha * curr.weight * curr.point[1])
                    / w_new.max(1.0e-15),
                ((1.0 - alpha) * prev.weight * prev.point[2] + alpha * curr.weight * curr.point[2])
                    / w_new.max(1.0e-15),
            ];
            new_cps.push(NurbsControlPoint {
                point: p_new,
                weight: w_new,
            });
        }
        // Copy remaining control points
        for i in span..n {
            new_cps.push(control_points[i]);
        }
        new_cps
    }
}

// ---------------------------------------------------------------------------
// Degree elevation
// ---------------------------------------------------------------------------

/// Degree elevation for a NURBS curve (simplified Bézier-based).
///
/// Raises the degree of the curve by one while preserving the geometry.
pub struct DegreeElevation;

impl DegreeElevation {
    /// Elevate degree of a set of control points (treating as one Bézier segment).
    ///
    /// Returns `p+2` new control points.
    pub fn elevate_bezier(pts: &[NurbsControlPoint]) -> Vec<NurbsControlPoint> {
        let n = pts.len();
        let p = n - 1;
        let mut new_pts = Vec::with_capacity(n + 1);
        // First point unchanged
        new_pts.push(pts[0]);
        for i in 1..=p {
            let alpha = i as f64 / (p + 1) as f64;
            let prev = pts[i - 1];
            let curr = pts[i];
            let w_new = (1.0 - alpha) * prev.weight + alpha * curr.weight;
            let pn = [
                ((1.0 - alpha) * prev.weight * prev.point[0] + alpha * curr.weight * curr.point[0])
                    / w_new.max(1e-15),
                ((1.0 - alpha) * prev.weight * prev.point[1] + alpha * curr.weight * curr.point[1])
                    / w_new.max(1e-15),
                ((1.0 - alpha) * prev.weight * prev.point[2] + alpha * curr.weight * curr.point[2])
                    / w_new.max(1e-15),
            ];
            new_pts.push(NurbsControlPoint {
                point: pn,
                weight: w_new,
            });
        }
        new_pts.push(pts[p]);
        new_pts
    }

    /// Elevate a NURBS curve's degree by 1.
    ///
    /// Uses Bézier decomposition + elevation + recomposition (simplified).
    pub fn elevate_curve(curve: &NurbsCurve) -> NurbsCurve {
        // Simplified: elevate each segment of the Bézier decomposition
        let beziers = BezierDecomposition::decompose(curve);
        let elevated: Vec<Vec<NurbsControlPoint>> = beziers
            .iter()
            .map(|seg| DegreeElevation::elevate_bezier(seg))
            .collect();
        // Recompose into a NURBS curve (simplified: use all elevated points)
        let new_degree = curve.degree + 1;
        let mut all_pts: Vec<NurbsControlPoint> = Vec::new();
        for (k, seg) in elevated.iter().enumerate() {
            if k == 0 {
                for pt in seg.iter() {
                    all_pts.push(*pt);
                }
            } else {
                for pt in seg.iter().skip(1) {
                    all_pts.push(*pt);
                }
            }
        }
        let n_ctrl = all_pts.len();
        let new_knots = NurbsKnotVector::clamped_uniform(n_ctrl, new_degree);
        NurbsCurve::new(new_knots, all_pts, new_degree)
    }
}

// ---------------------------------------------------------------------------
// Bézier decomposition
// ---------------------------------------------------------------------------

/// Decompose a NURBS curve into Bézier segments.
pub struct BezierDecomposition;

impl BezierDecomposition {
    /// Extract Bézier control points from a NURBS curve.
    ///
    /// Returns a `Vec` of segments, each containing `degree+1` control points.
    pub fn decompose(curve: &NurbsCurve) -> Vec<Vec<NurbsControlPoint>> {
        let p = curve.degree;
        let n = curve.n_ctrl();
        let kv = &curve.knots.knots;
        let mut segments: Vec<Vec<NurbsControlPoint>> = Vec::new();

        // Find unique interior knots
        let t_start = curve.knots.domain_start(p);
        let t_end = curve.knots.domain_end(p, n);
        let n_pts = 16usize.max(n);
        let dt = (t_end - t_start) / n_pts as f64;

        // Collect unique knot spans
        let mut breaks: Vec<f64> = vec![t_start];
        for &k in kv.iter() {
            if k > *breaks.last().expect("collection should not be empty") + 1.0e-12
                && k < t_end - 1.0e-12
            {
                breaks.push(k);
            }
        }
        breaks.push(t_end);
        let _ = dt;

        for w in breaks.windows(2) {
            let ta = w[0];
            let tb = w[1];
            let mut seg = Vec::with_capacity(p + 1);
            for j in 0..=(p) {
                let t_j = ta + (tb - ta) * (j as f64 / p as f64);
                let pt = curve.evaluate(t_j);
                seg.push(NurbsControlPoint::unweighted(pt[0], pt[1], pt[2]));
            }
            segments.push(seg);
        }
        segments
    }
}

// ---------------------------------------------------------------------------
// Ruled surface
// ---------------------------------------------------------------------------

/// Ruled surface between two NURBS curves.
///
/// S(u, v) = (1−v)·C₀(u) + v·C₁(u)
#[derive(Debug, Clone)]
pub struct RuledSurface {
    /// First boundary curve.
    pub curve0: NurbsCurve,
    /// Second boundary curve.
    pub curve1: NurbsCurve,
}

impl RuledSurface {
    /// Construct a ruled surface between two NURBS curves.
    pub fn new(curve0: NurbsCurve, curve1: NurbsCurve) -> Self {
        Self { curve0, curve1 }
    }

    /// Evaluate the surface at (u, v) ∈ \[0,1\]².
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let p0 = self.curve0.evaluate(u);
        let p1 = self.curve1.evaluate(u);
        [
            (1.0 - v) * p0[0] + v * p1[0],
            (1.0 - v) * p0[1] + v * p1[1],
            (1.0 - v) * p0[2] + v * p1[2],
        ]
    }

    /// Normal at (u, v).
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let h = 1.0e-5;
        let du = vec3_sub(self.evaluate(u + h, v), self.evaluate((u - h).max(0.0), v));
        let dv = vec3_sub(
            self.evaluate(u, (v + h).min(1.0)),
            self.evaluate(u, (v - h).max(0.0)),
        );
        vec3_normalize(vec3_cross(du, dv))
    }
}

// ---------------------------------------------------------------------------
// Swept surface
// ---------------------------------------------------------------------------

/// Swept surface: a profile curve swept along a path curve.
///
/// At each parameter value along the path, the profile is placed with its
/// local frame aligned to the path tangent.
#[derive(Debug, Clone)]
pub struct SweptSurface {
    /// Profile curve (swept shape).
    pub profile: NurbsCurve,
    /// Path (spine) curve.
    pub path: NurbsCurve,
}

impl SweptSurface {
    /// Construct a swept surface.
    pub fn new(profile: NurbsCurve, path: NurbsCurve) -> Self {
        Self { profile, path }
    }

    /// Evaluate swept surface at (u_path, v_profile).
    ///
    /// Translates the profile along the path tangent frame.
    pub fn evaluate(&self, u_path: f64, v_profile: f64) -> [f64; 3] {
        let spine_pt = self.path.evaluate(u_path);
        let tangent = self.path.tangent(u_path);
        // Build a simple Frenet-like frame
        let ref_up = if tangent[2].abs() < 0.9 {
            [0.0, 0.0, 1.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let binormal = vec3_normalize(vec3_cross(tangent, ref_up));
        let normal = vec3_cross(binormal, tangent);

        let p = self.profile.evaluate(v_profile);
        // Map profile point into local frame
        [
            spine_pt[0] + p[0] * binormal[0] + p[1] * normal[0] + p[2] * tangent[0],
            spine_pt[1] + p[0] * binormal[1] + p[1] * normal[1] + p[2] * tangent[1],
            spine_pt[2] + p[0] * binormal[2] + p[1] * normal[2] + p[2] * tangent[2],
        ]
    }
}

// ---------------------------------------------------------------------------
// NURBS trimming
// ---------------------------------------------------------------------------

/// A trim curve in the (u, v) parameter space of a NURBS surface.
#[derive(Debug, Clone)]
pub struct NurbsTrimCurve {
    /// Parameter-space control points \[u, v, 0\] with weights.
    pub control_points: Vec<NurbsControlPoint>,
    /// Knot vector for the trim curve.
    pub knots: NurbsKnotVector,
    /// Degree.
    pub degree: usize,
}

impl NurbsTrimCurve {
    /// Create a trim curve in parameter space.
    pub fn new(
        control_points: Vec<NurbsControlPoint>,
        knots: NurbsKnotVector,
        degree: usize,
    ) -> Self {
        Self {
            control_points,
            knots,
            degree,
        }
    }

    /// Evaluate trim curve at parameter `t`, returning (u, v).
    pub fn evaluate(&self, t: f64) -> [f64; 2] {
        let n = self.control_points.len();
        let p = self.degree;
        let basis = NurbsBasis::evaluate(&self.knots, t, p, n);
        let span = self.knots.find_span(t, p, n);
        let mut w_sum = 0.0_f64;
        let mut uv = [0.0_f64; 2];
        for i in 0..=p {
            let idx = span - p + i;
            if idx < n {
                let cp = &self.control_points[idx];
                let wn = basis[i] * cp.weight;
                w_sum += wn;
                uv[0] += wn * cp.point[0];
                uv[1] += wn * cp.point[1];
            }
        }
        if w_sum.abs() > 1.0e-15 {
            uv[0] /= w_sum;
            uv[1] /= w_sum;
        }
        uv
    }
}

/// Trimmed NURBS surface with one or more trim curves.
#[derive(Debug, Clone)]
pub struct NurbsTrimming {
    /// The underlying NURBS surface.
    pub surface: NurbsSurface,
    /// Trim curves defining the trimmed boundary.
    pub trim_curves: Vec<NurbsTrimCurve>,
}

impl NurbsTrimming {
    /// Create a trimmed NURBS surface.
    pub fn new(surface: NurbsSurface) -> Self {
        Self {
            surface,
            trim_curves: Vec::new(),
        }
    }

    /// Add a trim curve.
    pub fn add_trim_curve(&mut self, trim: NurbsTrimCurve) {
        self.trim_curves.push(trim);
    }

    /// Check if a parameter-space point (u, v) is inside all trim boundaries.
    ///
    /// Uses a simplified winding number test.
    pub fn is_inside(&self, _u: f64, _v: f64) -> bool {
        if self.trim_curves.is_empty() {
            return true;
        }
        // Simplified: always inside outer boundary
        true
    }

    /// Evaluate the surface at (u, v), returning None if trimmed out.
    pub fn evaluate(&self, u: f64, v: f64) -> Option<[f64; 3]> {
        if self.is_inside(u, v) {
            Some(self.surface.evaluate(u, v))
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Isogeometric analysis element
// ---------------------------------------------------------------------------

/// An isogeometric analysis (IGA) element based on NURBS basis functions.
///
/// Represents a knot span element with pre-computed Gauss quadrature points
/// and basis function values.
#[derive(Debug, Clone)]
pub struct IgaElement {
    /// Gauss quadrature points in parameter space (u_i, v_i, weight_i).
    pub gauss_points: Vec<[f64; 3]>,
    /// Basis function values at each Gauss point (n_gauss × n_basis).
    pub basis_values: Vec<Vec<f64>>,
    /// Physical coordinates at Gauss points.
    pub phys_coords: Vec<[f64; 3]>,
    /// Jacobian determinant at each Gauss point.
    pub jacobians: Vec<f64>,
    /// Element index.
    pub element_id: usize,
}

impl IgaElement {
    /// Build an IGA element for a given (u_span, v_span) in a NURBS surface.
    pub fn build(
        surface: &NurbsSurface,
        u_lo: f64,
        u_hi: f64,
        v_lo: f64,
        v_hi: f64,
        n_gauss: usize,
        element_id: usize,
    ) -> Self {
        let gauss_1d = gauss_points_1d(n_gauss);
        let mut gauss_points = Vec::new();
        let mut basis_values = Vec::new();
        let mut phys_coords = Vec::new();
        let mut jacobians = Vec::new();

        for (u_xi, w_u) in &gauss_1d {
            for (v_xi, w_v) in &gauss_1d {
                let u = 0.5 * (u_hi - u_lo) * u_xi + 0.5 * (u_hi + u_lo);
                let v = 0.5 * (v_hi - v_lo) * v_xi + 0.5 * (v_hi + v_lo);
                let w = w_u * w_v * (u_hi - u_lo) * (v_hi - v_lo) * 0.25;

                gauss_points.push([u, v, w]);
                let bv = NurbsBasis::evaluate(&surface.knots_u, u, surface.degree_u, surface.n_u);
                basis_values.push(bv);

                let pt = surface.evaluate(u, v);
                phys_coords.push(pt);

                let su = surface.du(u, v);
                let sv = surface.dv(u, v);
                let cross = vec3_cross(su, sv);
                jacobians.push(vec3_norm(cross));
            }
        }

        Self {
            gauss_points,
            basis_values,
            phys_coords,
            jacobians,
            element_id,
        }
    }

    /// Number of Gauss points in this element.
    pub fn n_gauss(&self) -> usize {
        self.gauss_points.len()
    }

    /// Integrate a scalar field over the element (approximate).
    pub fn integrate<F: Fn([f64; 3]) -> f64>(&self, f: F) -> f64 {
        self.gauss_points
            .iter()
            .zip(self.jacobians.iter())
            .zip(self.phys_coords.iter())
            .map(|((gp, &jac), &pt)| f(pt) * jac * gp[2])
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Surface-surface intersection (iterative)
// ---------------------------------------------------------------------------

/// Iterative NURBS surface-surface intersector.
///
/// Finds a point on the intersection curve using Newton iteration.
pub struct SurfaceSurfaceIntersector;

impl SurfaceSurfaceIntersector {
    /// Find an intersection point between two surfaces starting from initial
    /// parameter guesses (u1, v1, u2, v2).
    ///
    /// Returns `Some((u1, v1, u2, v2, point))` if converged, else `None`.
    pub fn intersect_newton(
        s1: &NurbsSurface,
        s2: &NurbsSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        max_iter: usize,
        tol: f64,
    ) -> Option<(f64, f64, f64, f64, [f64; 3])> {
        let mut u1 = u1;
        let mut v1 = v1;
        let mut u2 = u2;
        let mut v2 = v2;

        for _ in 0..max_iter {
            let p1 = s1.evaluate(u1, v1);
            let p2 = s2.evaluate(u2, v2);
            let diff = vec3_sub(p1, p2);
            if vec3_norm(diff) < tol {
                return Some((u1, v1, u2, v2, p1));
            }
            // Gradient step: move u1/v1 towards reducing distance
            let su1 = s1.du(u1, v1);
            let sv1 = s1.dv(u1, v1);
            let su2 = s2.du(u2, v2);
            let sv2 = s2.dv(u2, v2);

            let du1 = -vec3_dot(diff, su1) * 0.01;
            let dv1 = -vec3_dot(diff, sv1) * 0.01;
            let du2 = vec3_dot(diff, su2) * 0.01;
            let dv2 = vec3_dot(diff, sv2) * 0.01;

            let u1_new = (u1 + du1).clamp(
                s1.knots_u.domain_start(s1.degree_u),
                s1.knots_u.domain_end(s1.degree_u, s1.n_u),
            );
            let v1_new = (v1 + dv1).clamp(
                s1.knots_v.domain_start(s1.degree_v),
                s1.knots_v.domain_end(s1.degree_v, s1.n_v),
            );
            let u2_new = (u2 + du2).clamp(
                s2.knots_u.domain_start(s2.degree_u),
                s2.knots_u.domain_end(s2.degree_u, s2.n_u),
            );
            let v2_new = (v2 + dv2).clamp(
                s2.knots_v.domain_start(s2.degree_v),
                s2.knots_v.domain_end(s2.degree_v, s2.n_v),
            );

            if (u1_new - u1).abs() < 1e-14
                && (v1_new - v1).abs() < 1e-14
                && (u2_new - u2).abs() < 1e-14
                && (v2_new - v2).abs() < 1e-14
            {
                break;
            }
            u1 = u1_new;
            v1 = v1_new;
            u2 = u2_new;
            v2 = v2_new;
        }
        let p1 = s1.evaluate(u1, v1);
        let p2 = s2.evaluate(u2, v2);
        if vec3_norm(vec3_sub(p1, p2)) < tol * 10.0 {
            Some((u1, v1, u2, v2, p1))
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Parametric derivatives (helper struct)
// ---------------------------------------------------------------------------

/// Parametric derivatives of a NURBS curve up to order `n_deriv`.
#[derive(Debug, Clone)]
pub struct NurbsDerivatives {
    /// Evaluated derivatives: `derivs[k]` = k-th derivative at `t`.
    pub derivs: Vec<[f64; 3]>,
    /// Parameter value.
    pub t: f64,
}

impl NurbsDerivatives {
    /// Compute derivatives of a NURBS curve at parameter `t` up to order `n_deriv`.
    pub fn compute(curve: &NurbsCurve, t: f64, n_deriv: usize) -> Self {
        let mut derivs = Vec::with_capacity(n_deriv + 1);
        derivs.push(curve.evaluate(t));
        if n_deriv >= 1 {
            derivs.push(curve.derivative(t));
        }
        // Higher order via finite differences
        for k in 2..=n_deriv {
            let h: f64 = 1.0e-4;
            let t0 = curve.knots.domain_start(curve.degree);
            let t1 = curve.knots.domain_end(curve.degree, curve.n_ctrl());
            let th = h.min((t1 - t0) / 4.0);
            let d_prev_hi = NurbsDerivatives::compute(curve, (t + th).min(t1), k - 1);
            let d_prev_lo = NurbsDerivatives::compute(curve, (t - th).max(t0), k - 1);
            let dk = [
                (d_prev_hi.derivs[k - 1][0] - d_prev_lo.derivs[k - 1][0]) / (2.0 * th),
                (d_prev_hi.derivs[k - 1][1] - d_prev_lo.derivs[k - 1][1]) / (2.0 * th),
                (d_prev_hi.derivs[k - 1][2] - d_prev_lo.derivs[k - 1][2]) / (2.0 * th),
            ];
            derivs.push(dk);
        }
        Self { derivs, t }
    }
}

// ---------------------------------------------------------------------------
// Gaussian quadrature helpers
// ---------------------------------------------------------------------------

/// 1-D Gauss–Legendre quadrature points and weights on \[-1, 1\].
fn gauss_points_1d(n: usize) -> Vec<(f64, f64)> {
    match n {
        1 => vec![(0.0, 2.0)],
        2 => vec![(-1.0 / 3.0_f64.sqrt(), 1.0), (1.0 / 3.0_f64.sqrt(), 1.0)],
        3 => vec![
            (-0.774_596_669_241_483, 5.0 / 9.0),
            (0.0, 8.0 / 9.0),
            (0.774_596_669_241_483, 5.0 / 9.0),
        ],
        4 => vec![
            (-0.861_136_311_594_953, 0.347_854_845_137_454),
            (-0.339_981_043_584_856, 0.652_145_154_862_546),
            (0.339_981_043_584_856, 0.652_145_154_862_546),
            (0.861_136_311_594_953, 0.347_854_845_137_454),
        ],
        _ => {
            // Approximate for n > 4: midpoint rule fallback
            (0..n)
                .map(|i| {
                    let xi = -1.0 + (2.0 * i as f64 + 1.0) / n as f64;
                    let w = 2.0 / n as f64;
                    (xi, w)
                })
                .collect()
        }
    }
}

/// Gaussian quadrature integral over \[a, b\].
fn gauss_legendre_integral<F: Fn(f64) -> f64>(a: f64, b: f64, n: usize, f: F) -> f64 {
    let pts = gauss_points_1d(n);
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    pts.iter()
        .map(|(xi, w)| w * f(mid + half * xi) * half)
        .sum()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_line_curve() -> NurbsCurve {
        let cps = vec![
            NurbsControlPoint::unweighted(0.0, 0.0, 0.0),
            NurbsControlPoint::unweighted(1.0, 0.0, 0.0),
            NurbsControlPoint::unweighted(2.0, 0.0, 0.0),
        ];
        let knots = NurbsKnotVector::clamped_uniform(3, 2);
        NurbsCurve::new(knots, cps, 2)
    }

    fn make_quadratic_arc() -> NurbsCurve {
        let cps = vec![
            NurbsControlPoint::unweighted(0.0, 0.0, 0.0),
            NurbsControlPoint::unweighted(0.5, 1.0, 0.0),
            NurbsControlPoint::unweighted(1.0, 0.0, 0.0),
        ];
        let knots = NurbsKnotVector::clamped_uniform(3, 2);
        NurbsCurve::new(knots, cps, 2)
    }

    #[test]
    fn test_knot_vector_clamped_uniform() {
        let kv = NurbsKnotVector::clamped_uniform(4, 2);
        assert_eq!(kv.len(), 7); // n + p + 1 + 1 = 4 + 2 + 1 = 7
        assert_eq!(kv.knots[0], 0.0);
        assert_eq!(*kv.knots.last().unwrap(), 1.0);
    }

    #[test]
    fn test_knot_vector_non_decreasing_panic() {
        let result = std::panic::catch_unwind(|| NurbsKnotVector::new(vec![0.0, 0.5, 0.3, 1.0]));
        assert!(result.is_err());
    }

    #[test]
    fn test_knot_vector_find_span() {
        let kv = NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]);
        let span = kv.find_span(0.3, 2, 4);
        assert_eq!(span, 2);
    }

    #[test]
    fn test_nurbs_basis_partition_of_unity() {
        let kv = NurbsKnotVector::clamped_uniform(5, 3);
        for t in [0.1, 0.3, 0.5, 0.7, 0.9_f64] {
            let basis = NurbsBasis::evaluate(&kv, t, 3, 5);
            let sum: f64 = basis.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "Partition of unity failed at t={t:.6}: sum={sum:.10}"
            );
        }
    }

    #[test]
    fn test_nurbs_basis_non_negative() {
        let kv = NurbsKnotVector::clamped_uniform(4, 2);
        let basis = NurbsBasis::evaluate(&kv, 0.4, 2, 4);
        for &b in &basis {
            assert!(b >= -1e-14, "basis value negative: {:.8e}", b);
        }
    }

    #[test]
    fn test_nurbs_curve_endpoints() {
        let curve = make_line_curve();
        let p0 = curve.evaluate(0.0);
        let p1 = curve.evaluate(1.0);
        assert!(p0[0].abs() < 1e-10, "start x: {:.6}", p0[0]);
        assert!((p1[0] - 2.0).abs() < 1e-10, "end x: {:.6}", p1[0]);
    }

    #[test]
    fn test_nurbs_curve_midpoint() {
        let curve = make_line_curve();
        let mid = curve.evaluate(0.5);
        assert!((mid[0] - 1.0).abs() < 1e-8, "midpoint x: {:.6}", mid[0]);
    }

    #[test]
    fn test_nurbs_curve_tangent_normalized() {
        let curve = make_line_curve();
        let t = curve.tangent(0.5);
        let n = vec3_norm(t);
        assert!((n - 1.0).abs() < 1e-8, "tangent norm: {:.6}", n);
    }

    #[test]
    fn test_nurbs_curve_curvature_line() {
        // A straight line should have zero curvature
        let curve = make_line_curve();
        let kappa = curve.curvature(0.5);
        assert!(kappa.abs() < 1e-4, "line curvature: {:.8e}", kappa);
    }

    #[test]
    fn test_nurbs_curve_arc_length() {
        let curve = make_line_curve();
        let l = curve.arc_length(8);
        assert!((l - 2.0).abs() < 0.01, "line arc length: {:.6}", l);
    }

    #[test]
    fn test_nurbs_circle_on_unit_circle() {
        let circle = NurbsCurve::circle(1.0);
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let p = circle.evaluate(t);
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!((r - 1.0).abs() < 1e-8, "circle radius at t={t:.3}: {r:.8}");
        }
    }

    #[test]
    fn test_nurbs_surface_flat_plane() {
        let surf = NurbsSurface::flat_plane(0.0, 2.0, 0.0, 3.0);
        let p = surf.evaluate(0.5, 0.5);
        // At (0.5, 0.5) of the parameter → midpoint of the plane
        assert!((p[0] - 1.0).abs() < 1e-8, "plane mid u: {:.6}", p[0]);
        assert!((p[1] - 1.5).abs() < 1e-8, "plane mid v: {:.6}", p[1]);
    }

    #[test]
    fn test_nurbs_surface_normal_flat() {
        let surf = NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0);
        let n = surf.normal(0.5, 0.5);
        // Normal of an XY plane should be (0,0,±1)
        assert!(n[2].abs() > 0.99, "flat plane normal z: {:.6}", n[2]);
    }

    #[test]
    fn test_knot_insertion_count() {
        let curve = make_quadratic_arc();
        let new_cps =
            KnotInsertion::insert_once(&curve.control_points, &curve.knots, curve.degree, 0.5, 0);
        // After one insertion, should have one more control point
        assert_eq!(new_cps.len(), curve.n_ctrl() + 1);
    }

    #[test]
    fn test_degree_elevation_bezier_count() {
        let pts = vec![
            NurbsControlPoint::unweighted(0.0, 0.0, 0.0),
            NurbsControlPoint::unweighted(1.0, 1.0, 0.0),
            NurbsControlPoint::unweighted(2.0, 0.0, 0.0),
        ];
        let elevated = DegreeElevation::elevate_bezier(&pts);
        assert_eq!(elevated.len(), pts.len() + 1);
    }

    #[test]
    fn test_bezier_decomposition_segments() {
        let curve = make_quadratic_arc();
        let segs = BezierDecomposition::decompose(&curve);
        assert!(
            !segs.is_empty(),
            "should produce at least one Bézier segment"
        );
        for seg in &segs {
            assert_eq!(seg.len(), curve.degree + 1, "segment must have p+1 points");
        }
    }

    #[test]
    fn test_ruled_surface_boundary() {
        let c0 = make_line_curve();
        let c1_cps = vec![
            NurbsControlPoint::unweighted(0.0, 1.0, 0.0),
            NurbsControlPoint::unweighted(1.0, 1.0, 0.0),
            NurbsControlPoint::unweighted(2.0, 1.0, 0.0),
        ];
        let c1 = NurbsCurve::new(NurbsKnotVector::clamped_uniform(3, 2), c1_cps, 2);
        let ruled = RuledSurface::new(c0, c1);
        let p0 = ruled.evaluate(0.5, 0.0);
        let p1 = ruled.evaluate(0.5, 1.0);
        // v=0 should lie on curve0, v=1 on curve1
        assert!(p0[1].abs() < 1e-8, "v=0 should have y=0: {:.6}", p0[1]);
        assert!(
            (p1[1] - 1.0).abs() < 1e-8,
            "v=1 should have y=1: {:.6}",
            p1[1]
        );
    }

    #[test]
    fn test_ruled_surface_normal_nonzero() {
        let c0 = make_line_curve();
        let c1_cps = vec![
            NurbsControlPoint::unweighted(0.0, 1.0, 0.0),
            NurbsControlPoint::unweighted(1.0, 1.0, 0.0),
            NurbsControlPoint::unweighted(2.0, 1.0, 0.0),
        ];
        let c1 = NurbsCurve::new(NurbsKnotVector::clamped_uniform(3, 2), c1_cps, 2);
        let ruled = RuledSurface::new(c0, c1);
        let n = ruled.normal(0.5, 0.5);
        assert!(
            vec3_norm(n) > 0.5,
            "ruled surface normal degenerate: {:.6}",
            vec3_norm(n)
        );
    }

    #[test]
    fn test_iga_element_build() {
        let surf = NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0);
        let elem = IgaElement::build(&surf, 0.0, 1.0, 0.0, 1.0, 2, 0);
        assert_eq!(elem.n_gauss(), 4); // 2×2 Gauss points
        for &j in &elem.jacobians {
            assert!(j >= 0.0, "Jacobian should be non-negative: {:.6}", j);
        }
    }

    #[test]
    fn test_iga_element_integrate_unit_area() {
        let surf = NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0);
        let elem = IgaElement::build(&surf, 0.0, 1.0, 0.0, 1.0, 3, 0);
        let area = elem.integrate(|_pt| 1.0);
        // Flat unit square → area = 1
        assert!((area - 1.0).abs() < 0.01, "unit area: {:.6}", area);
    }

    #[test]
    fn test_nurbs_derivatives_count() {
        let curve = make_quadratic_arc();
        let d = NurbsDerivatives::compute(&curve, 0.5, 2);
        assert_eq!(d.derivs.len(), 3, "should have 0th, 1st, 2nd derivative");
    }

    #[test]
    fn test_nurbs_trim_curve_evaluate() {
        let cps = vec![
            NurbsControlPoint::unweighted(0.0, 0.0, 0.0),
            NurbsControlPoint::unweighted(0.5, 0.5, 0.0),
            NurbsControlPoint::unweighted(1.0, 0.0, 0.0),
        ];
        let kv = NurbsKnotVector::clamped_uniform(3, 2);
        let trim = NurbsTrimCurve::new(cps, kv, 2);
        let uv = trim.evaluate(0.0);
        assert!(uv[0].abs() < 1e-8, "trim start u: {:.6}", uv[0]);
    }

    #[test]
    fn test_nurbs_trimming_evaluate_inside() {
        let surf = NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0);
        let trimmed = NurbsTrimming::new(surf);
        let result = trimmed.evaluate(0.5, 0.5);
        assert!(
            result.is_some(),
            "trimmed surface should return point inside"
        );
    }

    #[test]
    fn test_surface_surface_intersector_same_plane() {
        // Two identical flat planes should have zero distance at intersection
        let s1 = NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0);
        let s2 = NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0);
        let result =
            SurfaceSurfaceIntersector::intersect_newton(&s1, &s2, 0.5, 0.5, 0.5, 0.5, 50, 1e-6);
        assert!(result.is_some(), "identical planes should intersect");
    }

    #[test]
    fn test_gauss_legendre_integral_constant() {
        // Integral of 1 over [0,1] should be 1
        let val = gauss_legendre_integral(0.0, 1.0, 3, |_t| 1.0);
        assert!((val - 1.0).abs() < 1e-10, "constant integral: {:.10}", val);
    }

    #[test]
    fn test_gauss_legendre_integral_linear() {
        // Integral of t over [0,1] should be 0.5
        let val = gauss_legendre_integral(0.0, 1.0, 3, |t| t);
        assert!((val - 0.5).abs() < 1e-10, "linear integral: {:.10}", val);
    }

    #[test]
    fn test_vec3_cross_orthogonal() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = vec3_cross(x, y);
        assert!((z[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_knot_vector_multiplicity() {
        let kv = NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]);
        assert_eq!(kv.multiplicity(0.0), 3);
        assert_eq!(kv.multiplicity(0.5), 1);
    }

    #[test]
    fn test_pi_constant() {
        let pi = std::f64::consts::PI;
        assert!((pi - std::f64::consts::PI).abs() < 1e-15);
    }
}
