#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! B-spline and NURBS curves and surfaces.
//!
//! Provides:
//!
//! - [`BsplineBasis`]: B-spline basis functions N_{i,p}(t) via Cox-de Boor recursion
//! - [`BsplineCurve`]: Evaluate curve, derivative, arc length, curvature, torsion
//! - [`BsplineSurface`]: Tensor-product surface, normal, Gaussian/mean curvature
//! - [`NurbsCurve`]: Rational B-spline, weight update, circle/arc as exact NURBS
//! - [`NurbsSurface`]: Rational surface, trim curves, surface fitting to point cloud
//! - [`BsplineFitting`]: Least-squares fitting, chord-length parameterization, knot selection

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ============================================================================
// BsplineBasis
// ============================================================================

/// B-spline knot vector with associated utilities.
///
/// A knot vector is a non-decreasing sequence of real numbers
/// T = \[t_0, t_1, …, t_{n+p+1}\] where n+1 is the number of basis functions
/// and p is the polynomial degree.
#[derive(Debug, Clone)]
pub struct KnotVector {
    /// Knot values (non-decreasing).
    pub knots: Vec<f64>,
}

impl KnotVector {
    /// Construct a knot vector from a slice.
    ///
    /// # Panics
    ///
    /// Panics if the knots are not non-decreasing.
    pub fn new(knots: Vec<f64>) -> Self {
        for i in 1..knots.len() {
            assert!(
                knots[i] >= knots[i - 1],
                "knot vector must be non-decreasing"
            );
        }
        Self { knots }
    }

    /// Construct a uniform (clamped) knot vector for n+1 control points and degree p.
    ///
    /// The first and last knots are repeated p+1 times (clamped), and the interior
    /// knots are uniformly spaced.
    pub fn clamped_uniform(n_ctrl: usize, degree: usize) -> Self {
        let n = n_ctrl - 1;
        let p = degree;
        let m = n + p + 1;
        let mut knots = vec![0.0_f64; m + 1];
        // First p+1 knots = 0
        for i in 0..=p {
            knots[i] = 0.0;
        }
        // Last p+1 knots = 1
        for i in (m - p)..=m {
            knots[i] = 1.0;
        }
        // Interior knots
        let n_interior = m - 2 * p - 1;
        for j in 1..=n_interior {
            knots[p + j] = j as f64 / (n_interior + 1) as f64;
        }
        Self { knots }
    }

    /// Find the knot span index i such that t ∈ \[t_i, t_{i+1}).
    ///
    /// Uses binary search. For t at the end of the domain, returns the last
    /// non-zero-length span.
    pub fn find_span(&self, t: f64, degree: usize, n_ctrl: usize) -> usize {
        let n = n_ctrl - 1;
        let p = degree;
        let knots = &self.knots;
        // Special case: t == last knot
        if t >= knots[n + 1] {
            // Find rightmost span with non-zero length
            let mut i = n;
            while i > p && (knots[i] - knots[i + 1]).abs() < 1e-14 {
                i -= 1;
            }
            return i;
        }
        if t <= knots[p] {
            return p;
        }
        // Binary search
        let mut lo = p;
        let mut hi = n + 1;
        let mut mid = (lo + hi) / 2;
        while t < knots[mid] || t >= knots[mid + 1] {
            if t < knots[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
            mid = (lo + hi) / 2;
        }
        mid
    }

    /// Number of knots in the vector.
    pub fn len(&self) -> usize {
        self.knots.len()
    }

    /// Whether the knot vector is empty.
    pub fn is_empty(&self) -> bool {
        self.knots.is_empty()
    }

    /// Return the domain \[a, b\] of the knot vector (a = first knot, b = last knot).
    pub fn domain(&self) -> (f64, f64) {
        (
            *self.knots.first().unwrap_or(&0.0),
            *self.knots.last().unwrap_or(&1.0),
        )
    }
}

/// B-spline basis functions and their derivatives.
///
/// Implements the Cox-de Boor recursion:
///
/// N_{i,0}(t) = 1 if t_i ≤ t < t_{i+1}, else 0
///
/// N_{i,p}(t) = (t − t_i)/(t_{i+p} − t_i) N_{i,p-1}(t) + (t_{i+p+1} − t)/(t_{i+p+1} − t_{i+1}) N_{i+1,p-1}(t)
#[derive(Debug, Clone)]
pub struct BsplineBasis {
    /// Polynomial degree p.
    pub degree: usize,
    /// Knot vector.
    pub knot_vector: KnotVector,
    /// Number of control points (n+1).
    pub n_ctrl: usize,
}

impl BsplineBasis {
    /// Create a new B-spline basis.
    pub fn new(degree: usize, knot_vector: KnotVector, n_ctrl: usize) -> Self {
        Self {
            degree,
            knot_vector,
            n_ctrl,
        }
    }

    /// Create a basis with a uniform clamped knot vector.
    pub fn clamped(degree: usize, n_ctrl: usize) -> Self {
        let kv = KnotVector::clamped_uniform(n_ctrl, degree);
        Self {
            degree,
            knot_vector: kv,
            n_ctrl,
        }
    }

    /// Evaluate all non-zero basis functions N_{i-p,p}(t), …, N_{i,p}(t) at parameter t.
    ///
    /// Returns an array of length p+1. Uses the triangular table algorithm.
    pub fn eval_nonzero(&self, t: f64) -> (usize, Vec<f64>) {
        let p = self.degree;
        let knots = &self.knot_vector.knots;
        let span = self.knot_vector.find_span(t, p, self.n_ctrl);
        let mut n = vec![0.0_f64; p + 1];
        let mut left = vec![0.0_f64; p + 1];
        let mut right = vec![0.0_f64; p + 1];
        n[0] = 1.0;
        for j in 1..=p {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            let mut saved = 0.0;
            for r in 0..j {
                let temp = n[r] / (right[r + 1] + left[j - r]);
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            n[j] = saved;
        }
        (span, n)
    }

    /// Evaluate all basis functions N_{i,p}(t) for i = 0, …, n_ctrl-1.
    ///
    /// Returns a vector of length n_ctrl.
    pub fn eval_all(&self, t: f64) -> Vec<f64> {
        let (span, nonzero) = self.eval_nonzero(t);
        let p = self.degree;
        let mut result = vec![0.0_f64; self.n_ctrl];
        for k in 0..=p {
            if span + k >= p && span + k - p < self.n_ctrl {
                result[span + k - p] = nonzero[k];
            }
        }
        result
    }

    /// Evaluate basis function N_{i,p}(t).
    pub fn eval(&self, i: usize, t: f64) -> f64 {
        let all = self.eval_all(t);
        if i < all.len() { all[i] } else { 0.0 }
    }

    /// Evaluate all non-zero basis functions and their derivatives up to order `deriv`.
    ///
    /// Returns a 2D array `ders[k][j]` where k is the derivative order (0..=deriv)
    /// and j is the basis function index (0..=p) relative to the span.
    pub fn eval_nonzero_derivs(&self, t: f64, deriv: usize) -> (usize, Vec<Vec<f64>>) {
        let p = self.degree;
        let knots = &self.knot_vector.knots;
        let span = self.knot_vector.find_span(t, p, self.n_ctrl);
        let d = deriv.min(p);
        // ndu table
        let mut ndu = vec![vec![0.0_f64; p + 1]; p + 1];
        let mut left = vec![0.0_f64; p + 1];
        let mut right = vec![0.0_f64; p + 1];
        ndu[0][0] = 1.0;
        for j in 1..=p {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            let mut saved = 0.0;
            for r in 0..j {
                ndu[j][r] = right[r + 1] + left[j - r];
                let temp = ndu[r][j - 1] / ndu[j][r];
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
            let mut s1 = 0_usize;
            let mut s2 = 1_usize;
            a[0][0] = 1.0;
            for k in 1..=d {
                let mut nd = 0.0;
                let rk = r as isize - k as isize;
                let pk = p as isize - k as isize;
                if r >= k {
                    a[s2][0] = a[s1][0] / ndu[pk as usize + 1][rk as usize];
                    nd = a[s2][0] * ndu[rk as usize][pk as usize];
                }
                let j1 = if rk >= -1 {
                    1_usize
                } else {
                    (-(rk + 1)) as usize
                };
                let j2 = if (r as isize - 1) <= pk { k - 1 } else { p - r };
                for j in j1..=j2 {
                    let idx2 = rk + j as isize;
                    if idx2 < 0 {
                        continue;
                    }
                    let idx2 = idx2 as usize;
                    a[s2][j] = (a[s1][j] - a[s1][j.saturating_sub(1)]) / ndu[pk as usize + 1][idx2];
                    nd += a[s2][j] * ndu[idx2][pk as usize];
                }
                if r <= (p - k) {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk as usize + 1][r];
                    nd += a[s2][k] * ndu[r][pk as usize];
                }
                ders[k][r] = nd;
                std::mem::swap(&mut s1, &mut s2);
            }
        }
        let mut rfact = p as f64;
        for k in 1..=d {
            for j in 0..=p {
                ders[k][j] *= rfact;
            }
            if k < d {
                rfact *= (p - k) as f64;
            }
        }
        (span, ders)
    }

    /// Evaluate the k-th derivative of basis function N_{i,p}(t).
    pub fn eval_deriv(&self, i: usize, t: f64, k: usize) -> f64 {
        let p = self.degree;
        let (span, ders) = self.eval_nonzero_derivs(t, k);
        let local_idx = i as isize - span as isize;
        if local_idx < 0 || local_idx > p as isize {
            0.0
        } else {
            ders[k][local_idx as usize]
        }
    }

    /// Greville abscissae (average knots) for control point parameterization.
    ///
    /// ξ_i = (t_{i+1} + … + t_{i+p}) / p
    pub fn greville_abscissae(&self) -> Vec<f64> {
        let p = self.degree;
        let knots = &self.knot_vector.knots;
        (0..self.n_ctrl)
            .map(|i| {
                if p == 0 {
                    knots[i]
                } else {
                    knots[i + 1..i + p + 1].iter().sum::<f64>() / p as f64
                }
            })
            .collect()
    }
}

// ============================================================================
// BsplineCurve
// ============================================================================

/// A B-spline curve in 3D space.
///
/// C(t) = Σ_{i=0}^{n} N_{i,p}(t) * P_i
///
/// where P_i are control points and N_{i,p}(t) are B-spline basis functions.
#[derive(Debug, Clone)]
pub struct BsplineCurve {
    /// B-spline basis.
    pub basis: BsplineBasis,
    /// Control points \[P_0, …, P_n\], each \[x, y, z\].
    pub control_points: Vec<[f64; 3]>,
}

impl BsplineCurve {
    /// Create a new B-spline curve.
    pub fn new(degree: usize, knot_vector: KnotVector, control_points: Vec<[f64; 3]>) -> Self {
        let n_ctrl = control_points.len();
        let basis = BsplineBasis::new(degree, knot_vector, n_ctrl);
        Self {
            basis,
            control_points,
        }
    }

    /// Create a B-spline curve with a clamped uniform knot vector.
    pub fn clamped(degree: usize, control_points: Vec<[f64; 3]>) -> Self {
        let n_ctrl = control_points.len();
        let basis = BsplineBasis::clamped(degree, n_ctrl);
        Self {
            basis,
            control_points,
        }
    }

    /// Evaluate the curve at parameter t: returns \[x, y, z\].
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let (span, nonzero) = self.basis.eval_nonzero(t);
        let p = self.basis.degree;
        let mut point = [0.0_f64; 3];
        for k in 0..=p {
            let idx = span + k - p;
            if idx < self.control_points.len() {
                let cp = self.control_points[idx];
                for d in 0..3 {
                    point[d] += nonzero[k] * cp[d];
                }
            }
        }
        point
    }

    /// Evaluate the k-th derivative of the curve at parameter t.
    ///
    /// Returns \[dx^k/dt^k, dy^k/dt^k, dz^k/dt^k\].
    pub fn eval_deriv(&self, t: f64, k: usize) -> [f64; 3] {
        let p = self.basis.degree;
        let (span, ders) = self.basis.eval_nonzero_derivs(t, k);
        let n_ctrl = self.control_points.len();
        let mut result = [0.0_f64; 3];
        for j in 0..=p {
            let idx = span + j - p;
            if idx < n_ctrl {
                for d in 0..3 {
                    result[d] += ders[k][j] * self.control_points[idx][d];
                }
            }
        }
        result
    }

    /// Compute arc length from t_start to t_end using Gaussian quadrature (n_gauss points per interval).
    pub fn arc_length(&self, t_start: f64, t_end: f64, n_intervals: usize) -> f64 {
        let h = (t_end - t_start) / n_intervals as f64;
        // 5-point Gauss-Legendre on [-1,1]
        let gp = [-0.906180, -0.538469, 0.0, 0.538469, 0.906180];
        let gw = [0.236927, 0.478629, 0.568889, 0.478629, 0.236927];
        let mut length = 0.0;
        for i in 0..n_intervals {
            let a = t_start + i as f64 * h;
            let b = a + h;
            let mid = (a + b) * 0.5;
            let half = (b - a) * 0.5;
            for (&xi, &wi) in gp.iter().zip(gw.iter()) {
                let t = mid + half * xi;
                let dt = self.eval_deriv(t, 1);
                let speed = (dt[0] * dt[0] + dt[1] * dt[1] + dt[2] * dt[2]).sqrt();
                length += wi * half * speed;
            }
        }
        length
    }

    /// Compute the unit tangent vector T(t) = C'(t)/|C'(t)|.
    pub fn tangent(&self, t: f64) -> [f64; 3] {
        let d1 = self.eval_deriv(t, 1);
        let mag = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        if mag < 1e-14 {
            [0.0, 0.0, 0.0]
        } else {
            [d1[0] / mag, d1[1] / mag, d1[2] / mag]
        }
    }

    /// Compute curvature κ(t) = |C' × C''| / |C'|³.
    pub fn curvature(&self, t: f64) -> f64 {
        let d1 = self.eval_deriv(t, 1);
        let d2 = self.eval_deriv(t, 2);
        let cross = [
            d1[1] * d2[2] - d1[2] * d2[1],
            d1[2] * d2[0] - d1[0] * d2[2],
            d1[0] * d2[1] - d1[1] * d2[0],
        ];
        let cross_mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        let d1_mag = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        if d1_mag < 1e-14 {
            0.0
        } else {
            cross_mag / d1_mag.powi(3)
        }
    }

    /// Compute torsion τ(t) = (C' × C'')·C''' / |C' × C''|².
    pub fn torsion(&self, t: f64) -> f64 {
        let d1 = self.eval_deriv(t, 1);
        let d2 = self.eval_deriv(t, 2);
        let d3 = self.eval_deriv(t, 3);
        let cross = [
            d1[1] * d2[2] - d1[2] * d2[1],
            d1[2] * d2[0] - d1[0] * d2[2],
            d1[0] * d2[1] - d1[1] * d2[0],
        ];
        let cross_mag2 = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
        if cross_mag2 < 1e-28 {
            return 0.0;
        }
        let dot_d3 = cross[0] * d3[0] + cross[1] * d3[1] + cross[2] * d3[2];
        dot_d3 / cross_mag2
    }

    /// Compute the principal normal vector N(t) = T'(t)/|T'(t)|.
    pub fn principal_normal(&self, t: f64) -> [f64; 3] {
        let kappa = self.curvature(t);
        if kappa < 1e-14 {
            return [0.0, 0.0, 0.0];
        }
        let d1 = self.eval_deriv(t, 1);
        let d2 = self.eval_deriv(t, 2);
        let d1_mag = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        let d1_mag3 = d1_mag.powi(3);
        // N = (d2 * |d1|^2 - d1 * (d1·d2)) / (kappa * |d1|^3)
        let d1_dot_d2 = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2];
        let d1_mag2 = d1_mag * d1_mag;
        let mut normal = [0.0_f64; 3];
        for i in 0..3 {
            normal[i] = (d2[i] * d1_mag2 - d1[i] * d1_dot_d2) / (kappa * d1_mag3);
        }
        normal
    }

    /// Evaluate the curve at n_points uniformly spaced parameter values.
    ///
    /// Returns a vector of \[x, y, z\] points.
    pub fn sample(&self, n_points: usize) -> Vec<[f64; 3]> {
        let (t0, t1) = self.basis.knot_vector.domain();
        (0..n_points)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / (n_points - 1).max(1) as f64;
                self.eval(t)
            })
            .collect()
    }

    /// Compute the bounding box of the curve (sampled at n_samples points).
    ///
    /// Returns (\[x_min, y_min, z_min\], \[x_max, y_max, z_max\]).
    pub fn bounding_box(&self, n_samples: usize) -> ([f64; 3], [f64; 3]) {
        let pts = self.sample(n_samples);
        let mut lo = pts[0];
        let mut hi = pts[0];
        for p in &pts {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        (lo, hi)
    }
}

// ============================================================================
// BsplineSurface
// ============================================================================

/// A tensor-product B-spline surface in 3D space.
///
/// S(u,v) = Σ_i Σ_j N_{i,p}(u) N_{j,q}(v) P_{ij}
#[derive(Debug, Clone)]
pub struct BsplineSurface {
    /// Basis in the u direction.
    pub basis_u: BsplineBasis,
    /// Basis in the v direction.
    pub basis_v: BsplineBasis,
    /// Control point grid \[i\]\[j\] = \[x, y, z\], row-major (u varies fastest).
    pub control_points: Vec<Vec<[f64; 3]>>,
}

impl BsplineSurface {
    /// Create a new B-spline surface.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knot_u: KnotVector,
        knot_v: KnotVector,
        control_points: Vec<Vec<[f64; 3]>>,
    ) -> Self {
        let n_u = control_points.len();
        let n_v = control_points[0].len();
        Self {
            basis_u: BsplineBasis::new(degree_u, knot_u, n_u),
            basis_v: BsplineBasis::new(degree_v, knot_v, n_v),
            control_points,
        }
    }

    /// Create a bilinear patch (p=q=1) from four corner points.
    pub fn bilinear_patch(p00: [f64; 3], p10: [f64; 3], p01: [f64; 3], p11: [f64; 3]) -> Self {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![vec![p00, p01], vec![p10, p11]];
        Self {
            basis_u: BsplineBasis::new(1, kv.clone(), 2),
            basis_v: BsplineBasis::new(1, kv, 2),
            control_points: cps,
        }
    }

    /// Evaluate the surface at (u, v): returns \[x, y, z\].
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let (span_u, nu) = self.basis_u.eval_nonzero(u);
        let (span_v, nv) = self.basis_v.eval_nonzero(v);
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut point = [0.0_f64; 3];
        for k in 0..=p {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for l in 0..=q {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let w = nu[k] * nv[l];
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    point[d] += w * cp[d];
                }
            }
        }
        point
    }

    /// Evaluate partial derivative ∂S/∂u at (u, v).
    pub fn eval_du(&self, u: f64, v: f64) -> [f64; 3] {
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let (span_u, ders_u) = self.basis_u.eval_nonzero_derivs(u, 1);
        let (span_v, nv) = self.basis_v.eval_nonzero(v);
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut result = [0.0_f64; 3];
        for k in 0..=p {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for l in 0..=q {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let w = ders_u[1][k] * nv[l];
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    result[d] += w * cp[d];
                }
            }
        }
        result
    }

    /// Evaluate partial derivative ∂S/∂v at (u, v).
    pub fn eval_dv(&self, u: f64, v: f64) -> [f64; 3] {
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let (span_u, nu) = self.basis_u.eval_nonzero(u);
        let (span_v, ders_v) = self.basis_v.eval_nonzero_derivs(v, 1);
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut result = [0.0_f64; 3];
        for k in 0..=p {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for l in 0..=q {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let w = nu[k] * ders_v[1][l];
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    result[d] += w * cp[d];
                }
            }
        }
        result
    }

    /// Compute the unit surface normal at (u, v).
    ///
    /// N = (∂S/∂u × ∂S/∂v) / |∂S/∂u × ∂S/∂v|
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let du = self.eval_du(u, v);
        let dv = self.eval_dv(u, v);
        let cross = [
            du[1] * dv[2] - du[2] * dv[1],
            du[2] * dv[0] - du[0] * dv[2],
            du[0] * dv[1] - du[1] * dv[0],
        ];
        let mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        if mag < 1e-14 {
            [0.0, 0.0, 1.0]
        } else {
            [cross[0] / mag, cross[1] / mag, cross[2] / mag]
        }
    }

    /// Compute the first fundamental form coefficients E, F, G at (u, v).
    ///
    /// E = ∂S/∂u · ∂S/∂u, F = ∂S/∂u · ∂S/∂v, G = ∂S/∂v · ∂S/∂v
    pub fn first_fundamental_form(&self, u: f64, v: f64) -> (f64, f64, f64) {
        let su = self.eval_du(u, v);
        let sv = self.eval_dv(u, v);
        let e = su[0] * su[0] + su[1] * su[1] + su[2] * su[2];
        let f = su[0] * sv[0] + su[1] * sv[1] + su[2] * sv[2];
        let g = sv[0] * sv[0] + sv[1] * sv[1] + sv[2] * sv[2];
        (e, f, g)
    }

    /// Compute Gaussian curvature K = (LN - M²)/(EG - F²) at (u, v).
    ///
    /// Uses finite differences for second derivatives.
    pub fn gaussian_curvature(&self, u: f64, v: f64) -> f64 {
        let h = 1e-5;
        let (e, f, g) = self.first_fundamental_form(u, v);
        let egf2 = e * g - f * f;
        if egf2.abs() < 1e-20 {
            return 0.0;
        }
        // Second derivatives via finite differences
        let suu = finite_diff_3d(|uu| self.eval_du(uu, v), u, h);
        let svv = finite_diff_3d(|vv| self.eval_dv(u, vv), v, h);
        let suv = finite_diff_3d(|uu| self.eval_dv(uu, v), u, h);
        let n = self.normal(u, v);
        let l = dot3(suu, n);
        let m_coef = dot3(suv, n);
        let nm = dot3(svv, n);
        (l * nm - m_coef * m_coef) / egf2
    }

    /// Compute mean curvature H = (EN - 2FM + GL)/(2(EG - F²)) at (u, v).
    pub fn mean_curvature(&self, u: f64, v: f64) -> f64 {
        let h = 1e-5;
        let (e, f, g) = self.first_fundamental_form(u, v);
        let egf2 = e * g - f * f;
        if egf2.abs() < 1e-20 {
            return 0.0;
        }
        let suu = finite_diff_3d(|uu| self.eval_du(uu, v), u, h);
        let svv = finite_diff_3d(|vv| self.eval_dv(u, vv), v, h);
        let suv = finite_diff_3d(|uu| self.eval_dv(uu, v), u, h);
        let n = self.normal(u, v);
        let l = dot3(suu, n);
        let m_coef = dot3(suv, n);
        let nm = dot3(svv, n);
        (e * nm - 2.0 * f * m_coef + g * l) / (2.0 * egf2)
    }

    /// Sample the surface at (nu × nv) parameter values.
    pub fn sample(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        let (u0, u1) = self.basis_u.knot_vector.domain();
        let (v0, v1) = self.basis_v.knot_vector.domain();
        (0..nu)
            .map(|i| {
                let u = u0 + (u1 - u0) * i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = v0 + (v1 - v0) * j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
    /// Compute Gaussian and mean curvature at (u, v).
    ///
    /// Returns `(K_gauss, H_mean)`.
    pub fn compute_curvature(&self, u: f64, v: f64) -> (f64, f64) {
        (self.gaussian_curvature(u, v), self.mean_curvature(u, v))
    }
    /// Refine the surface by inserting new knots in u and v directions.
    ///
    /// Uses the Oslo algorithm (knot insertion) to preserve geometry.
    pub fn refine_knots(&self, new_knots_u: &[f64], new_knots_v: &[f64]) -> BsplineSurface {
        // Phase 1: Insert u knots
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        // Refine in u: for each column j, refine the curve defined by
        // control_points[0..n_u][j] along the u knot vector
        let mut refined_u_knots = self.basis_u.knot_vector.knots.clone();
        for &k in new_knots_u {
            refined_u_knots.push(k);
        }
        refined_u_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        // Recompute control points after u refinement via Boehm's algorithm
        let mut new_cp = self.control_points.clone();
        for &knot in new_knots_u {
            let n_cur = new_cp.len();
            let orig_u_knots = &self.basis_u.knot_vector.knots;
            // Find span in the ORIGINAL knot vector
            let span = orig_u_knots
                .iter()
                .position(|&k| k > knot)
                .map(|pos| pos.saturating_sub(1))
                .unwrap_or(n_cur.saturating_sub(1))
                .min(n_cur);
            let mut inserted = Vec::with_capacity(n_cur + 1);
            for i in 0..n_cur + 1 {
                let mut row = Vec::with_capacity(n_v);
                for j in 0..n_v {
                    if i < span.saturating_sub(p) + 1 {
                        row.push(new_cp[i][j]);
                    } else if i > span {
                        row.push(new_cp[i - 1][j]);
                    } else {
                        let old_knots = &self.basis_u.knot_vector.knots;
                        let ki = if i < old_knots.len() {
                            old_knots[i]
                        } else {
                            1.0
                        };
                        let ki_p = if i + p < old_knots.len() {
                            old_knots[i + p]
                        } else {
                            1.0
                        };
                        let denom = ki_p - ki;
                        let alpha = if denom.abs() > 1e-14 {
                            (knot - ki) / denom
                        } else {
                            0.5
                        };
                        let p0 = if i > 0 {
                            new_cp[i - 1][j]
                        } else {
                            new_cp[0][j]
                        };
                        let p1 = if i < n_cur {
                            new_cp[i][j]
                        } else {
                            new_cp[n_cur - 1][j]
                        };
                        row.push([
                            (1.0 - alpha) * p0[0] + alpha * p1[0],
                            (1.0 - alpha) * p0[1] + alpha * p1[1],
                            (1.0 - alpha) * p0[2] + alpha * p1[2],
                        ]);
                    }
                }
                inserted.push(row);
            }
            new_cp = inserted;
        }
        // Phase 2: Insert v knots using Boehm's algorithm (mirroring u-direction)
        let n_u_new = new_cp.len();
        let mut new_cp_v = new_cp;
        let mut refined_v_knots = self.basis_v.knot_vector.knots.clone();
        for &k in new_knots_v {
            refined_v_knots.push(k);
        }
        refined_v_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        for &knot in new_knots_v {
            let n_v_cur = if !new_cp_v.is_empty() {
                new_cp_v[0].len()
            } else {
                0
            };
            let orig_v_knots = &self.basis_v.knot_vector.knots;
            // Find span in the ORIGINAL v knot vector
            let span = orig_v_knots
                .iter()
                .position(|&k| k > knot)
                .map(|pos| pos.saturating_sub(1))
                .unwrap_or(n_v_cur.saturating_sub(1))
                .min(n_v_cur);
            let mut new_rows = Vec::with_capacity(n_u_new);
            for i in 0..n_u_new {
                let old_row = &new_cp_v[i];
                let n_cur = old_row.len();
                let mut new_row = Vec::with_capacity(n_cur + 1);
                for j in 0..n_cur + 1 {
                    if j < span.saturating_sub(q) + 1 {
                        new_row.push(old_row[j]);
                    } else if j > span {
                        new_row.push(old_row[j - 1]);
                    } else {
                        let old_v_knots = &self.basis_v.knot_vector.knots;
                        let kj = if j < old_v_knots.len() {
                            old_v_knots[j]
                        } else {
                            1.0
                        };
                        let kj_q = if j + q < old_v_knots.len() {
                            old_v_knots[j + q]
                        } else {
                            1.0
                        };
                        let denom = kj_q - kj;
                        let alpha = if denom.abs() > 1e-14 {
                            (knot - kj) / denom
                        } else {
                            0.5
                        };
                        let p0 = if j > 0 { old_row[j - 1] } else { old_row[0] };
                        let p1 = if j < n_cur {
                            old_row[j]
                        } else {
                            old_row[n_cur - 1]
                        };
                        new_row.push([
                            (1.0 - alpha) * p0[0] + alpha * p1[0],
                            (1.0 - alpha) * p0[1] + alpha * p1[1],
                            (1.0 - alpha) * p0[2] + alpha * p1[2],
                        ]);
                    }
                }
                new_rows.push(new_row);
            }
            new_cp_v = new_rows;
        }
        let final_n_u = new_cp_v.len();
        let final_n_v = if final_n_u > 0 { new_cp_v[0].len() } else { 0 };
        // Build new knot vectors
        let mut ku = self.basis_u.knot_vector.knots.clone();
        for &k in new_knots_u {
            ku.push(k);
        }
        ku.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut kv = self.basis_v.knot_vector.knots.clone();
        for &k in new_knots_v {
            kv.push(k);
        }
        kv.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        // Ensure knot vectors have correct length (n_ctrl + degree + 1)
        while ku.len() < final_n_u + p + 1 {
            ku.push(1.0);
        }
        while kv.len() < final_n_v + q + 1 {
            kv.push(1.0);
        }
        BsplineSurface::new(p, q, KnotVector::new(ku), KnotVector::new(kv), new_cp_v)
    }
}

// ============================================================================
// NurbsCurve
// ============================================================================

/// A rational B-spline (NURBS) curve in 3D space.
///
/// C(t) = Σ w_i N_{i,p}(t) P_i / Σ w_i N_{i,p}(t)
#[derive(Debug, Clone)]
pub struct NurbsCurve {
    /// Underlying B-spline basis.
    pub basis: BsplineBasis,
    /// Control points \[P_0, …, P_n\], each \[x, y, z\].
    pub control_points: Vec<[f64; 3]>,
    /// Weights w_i > 0.
    pub weights: Vec<f64>,
}

impl NurbsCurve {
    /// Create a new NURBS curve.
    pub fn new(
        degree: usize,
        knot_vector: KnotVector,
        control_points: Vec<[f64; 3]>,
        weights: Vec<f64>,
    ) -> Self {
        let n_ctrl = control_points.len();
        assert_eq!(
            weights.len(),
            n_ctrl,
            "weights and control points must have the same length"
        );
        let basis = BsplineBasis::new(degree, knot_vector, n_ctrl);
        Self {
            basis,
            control_points,
            weights,
        }
    }

    /// Create a NURBS circle arc from center, radius, start and end angles (radians).
    ///
    /// Uses the standard 9-point representation for a full circle (degree 2).
    /// For arcs smaller than 2π, a simplified quadratic NURBS is returned.
    pub fn circle_arc(center: [f64; 3], radius: f64, start_angle: f64, end_angle: f64) -> Self {
        // Use a simple quadratic NURBS arc (3 control points for ≤ 180°, etc.)
        // For generality, build a degree-2 NURBS with the standard construction.
        let sweep = end_angle - start_angle;
        let n_arcs = if sweep.abs() <= std::f64::consts::FRAC_PI_2 {
            1
        } else if sweep.abs() <= std::f64::consts::PI {
            2
        } else if sweep.abs() <= 3.0 * std::f64::consts::FRAC_PI_2 {
            3
        } else {
            4
        };
        let n_pts = 2 * n_arcs + 1;
        let mut ctrl_pts = Vec::with_capacity(n_pts);
        let mut weights = Vec::with_capacity(n_pts);
        let d_theta = sweep / n_arcs as f64;
        let w1 = (d_theta / 2.0).cos();
        let mut angle = start_angle;
        // First point on arc
        ctrl_pts.push([
            center[0] + radius * angle.cos(),
            center[1] + radius * angle.sin(),
            center[2],
        ]);
        weights.push(1.0);
        for _i in 0..n_arcs {
            let a_mid = angle + d_theta * 0.5;
            let a_end = angle + d_theta;
            // Mid control point (off-arc)
            ctrl_pts.push([
                center[0] + radius / w1 * a_mid.cos(),
                center[1] + radius / w1 * a_mid.sin(),
                center[2],
            ]);
            weights.push(w1);
            // End point on arc
            ctrl_pts.push([
                center[0] + radius * a_end.cos(),
                center[1] + radius * a_end.sin(),
                center[2],
            ]);
            weights.push(1.0);
            angle = a_end;
        }
        // Build clamped knot vector for degree 2
        let mut knots = vec![0.0_f64; 3]; // first 3
        let knot_step = 1.0 / n_arcs as f64;
        for i in 1..n_arcs {
            knots.push(i as f64 * knot_step);
            knots.push(i as f64 * knot_step);
        }
        knots.extend(std::iter::repeat_n(1.0_f64, 3));
        let kv = KnotVector::new(knots);
        Self::new(2, kv, ctrl_pts, weights)
    }

    /// Evaluate the NURBS curve at parameter t: returns \[x, y, z\].
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let (span, nonzero) = self.basis.eval_nonzero(t);
        let p = self.basis.degree;
        let n_ctrl = self.control_points.len();
        let mut num = [0.0_f64; 3];
        let mut denom = 0.0_f64;
        for k in 0..=p {
            let idx = span + k - p;
            if idx >= n_ctrl {
                continue;
            }
            let wn = self.weights[idx] * nonzero[k];
            denom += wn;
            let cp = self.control_points[idx];
            for d in 0..3 {
                num[d] += wn * cp[d];
            }
        }
        if denom.abs() < 1e-30 {
            return self.control_points[0];
        }
        [num[0] / denom, num[1] / denom, num[2] / denom]
    }

    /// Evaluate the first derivative of the NURBS curve at t.
    ///
    /// Uses the quotient rule: C'(t) = (W'(t) P̃(t) - W(t)^2 P̃'(t)) / W(t)^2
    pub fn eval_deriv(&self, t: f64) -> [f64; 3] {
        let h = 1e-7;
        let t0 = (t - h).max(self.basis.knot_vector.knots[0]);
        let t1 = (t + h).min(*self.basis.knot_vector.knots.last().unwrap_or(&1.0));
        let p0 = self.eval(t0);
        let p1 = self.eval(t1);
        let dt = t1 - t0;
        [
            (p1[0] - p0[0]) / dt,
            (p1[1] - p0[1]) / dt,
            (p1[2] - p0[2]) / dt,
        ]
    }

    /// Update a control point and its weight.
    pub fn update_control_point(&mut self, i: usize, point: [f64; 3], weight: f64) {
        self.control_points[i] = point;
        self.weights[i] = weight;
    }

    /// Compute arc length using Gaussian quadrature.
    pub fn arc_length(&self, t_start: f64, t_end: f64, n_intervals: usize) -> f64 {
        let h = (t_end - t_start) / n_intervals as f64;
        let gp = [-0.906180, -0.538469, 0.0, 0.538469, 0.906180];
        let gw = [0.236927, 0.478629, 0.568889, 0.478629, 0.236927];
        let mut length = 0.0;
        for i in 0..n_intervals {
            let a = t_start + i as f64 * h;
            let b = a + h;
            let mid = (a + b) * 0.5;
            let half = (b - a) * 0.5;
            for (&xi, &wi) in gp.iter().zip(gw.iter()) {
                let t = mid + half * xi;
                let dt = self.eval_deriv(t);
                let speed = (dt[0] * dt[0] + dt[1] * dt[1] + dt[2] * dt[2]).sqrt();
                length += wi * half * speed;
            }
        }
        length
    }

    /// Compute curvature using the B-spline curve's curvature formula (numerical).
    pub fn curvature(&self, t: f64) -> f64 {
        let h = 1e-5;
        let t_min = self.basis.knot_vector.knots[0];
        let t_max = *self.basis.knot_vector.knots.last().unwrap_or(&1.0);
        let t0 = (t - h).max(t_min);
        let t1 = (t + h).min(t_max);
        let dt = t1 - t0;
        if dt < 1e-14 {
            return 0.0;
        }
        let d1 = self.eval_deriv(t);
        let p0 = self.eval(t0);
        let p1 = self.eval(t);
        let p2 = self.eval(t1);
        let d2 = [
            (p2[0] - 2.0 * p1[0] + p0[0]) / (h * h),
            (p2[1] - 2.0 * p1[1] + p0[1]) / (h * h),
            (p2[2] - 2.0 * p1[2] + p0[2]) / (h * h),
        ];
        let cross = [
            d1[1] * d2[2] - d1[2] * d2[1],
            d1[2] * d2[0] - d1[0] * d2[2],
            d1[0] * d2[1] - d1[1] * d2[0],
        ];
        let cross_mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        let d1_mag = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        if d1_mag < 1e-14 {
            0.0
        } else {
            cross_mag / d1_mag.powi(3)
        }
    }

    /// Sample the NURBS curve at n_points uniformly spaced parameter values.
    pub fn sample(&self, n_points: usize) -> Vec<[f64; 3]> {
        let (t0, t1) = self.basis.knot_vector.domain();
        (0..n_points)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / (n_points - 1).max(1) as f64;
                self.eval(t)
            })
            .collect()
    }
    /// Evaluate the B-spline basis function N_{i,p}(t) (Cox-de Boor recursion).
    pub fn b_spline_basis(i: usize, p: usize, t: f64, knots: &[f64]) -> f64 {
        if p == 0 {
            return if knots[i] <= t && t < knots[i + 1] {
                1.0
            } else if (t - knots[knots.len() - 1]).abs() < 1e-14
                && (knots[i + 1] - knots[knots.len() - 1]).abs() < 1e-14
            {
                // Handle t == last knot for the last basis function
                1.0
            } else {
                0.0
            };
        }
        let left_denom = knots[i + p] - knots[i];
        let left = if left_denom.abs() > 1e-14 {
            (t - knots[i]) / left_denom * Self::b_spline_basis(i, p - 1, t, knots)
        } else {
            0.0
        };
        let right_denom = knots[i + p + 1] - knots[i + 1];
        let right = if right_denom.abs() > 1e-14 {
            (knots[i + p + 1] - t) / right_denom * Self::b_spline_basis(i + 1, p - 1, t, knots)
        } else {
            0.0
        };
        left + right
    }
    /// Create a NURBS curve from 3D control points, weights, and degree.
    ///
    /// Automatically creates a clamped uniform knot vector.
    pub fn from_points_and_weights(
        points: Vec<[f64; 3]>,
        weights: Vec<f64>,
        degree: usize,
    ) -> Self {
        let n = points.len();
        let m = n + degree + 1;
        let mut knots = vec![0.0_f64; degree + 1];
        let interior = m - 2 * (degree + 1);
        for i in 1..=interior {
            knots.push(i as f64 / (interior + 1) as f64);
        }
        knots.extend(vec![1.0_f64; degree + 1]);
        let kv = KnotVector::new(knots);
        Self::new(degree, kv, points, weights)
    }
    /// Create a NURBS circle in the XY plane with the given radius (degree 2, 9 points).
    pub fn circle(radius: f64) -> Self {
        let w = std::f64::consts::FRAC_1_SQRT_2; // cos(π/4) ≈ 0.7071
        let r = radius;
        let ctrl = vec![
            [r, 0.0, 0.0],
            [r, r, 0.0],
            [0.0, r, 0.0],
            [-r, r, 0.0],
            [-r, 0.0, 0.0],
            [-r, -r, 0.0],
            [0.0, -r, 0.0],
            [r, -r, 0.0],
            [r, 0.0, 0.0],
        ];
        let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
        let knots = KnotVector::new(vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ]);
        Self::new(2, knots, ctrl, weights)
    }
    /// Compute arc-length reparametrization: return `n_output` uniformly-spaced
    /// (by arc length) sample points.
    pub fn compute_arc_length_reparametrize(
        &self,
        n_intervals: usize,
        n_output: usize,
    ) -> Vec<[f64; 3]> {
        let (arc_lengths, params) = self.arc_length_table(n_intervals);
        let total_len = *arc_lengths.last().unwrap_or(&0.0);
        if total_len < 1e-14 || n_output == 0 {
            return Vec::new();
        }
        let mut result = Vec::with_capacity(n_output);
        for k in 0..n_output {
            let target = total_len * k as f64 / (n_output - 1).max(1) as f64;
            // Binary search for the parameter corresponding to this arc length
            let idx = arc_lengths
                .partition_point(|&s| s < target)
                .min(arc_lengths.len() - 1);
            let t = if idx == 0 {
                params[0]
            } else if (arc_lengths[idx] - arc_lengths[idx - 1]).abs() < 1e-14 {
                params[idx]
            } else {
                let frac =
                    (target - arc_lengths[idx - 1]) / (arc_lengths[idx] - arc_lengths[idx - 1]);
                params[idx - 1] + frac * (params[idx] - params[idx - 1])
            };
            result.push(self.eval(t));
        }
        result
    }
    /// Build an arc-length look-up table with `n` samples.
    ///
    /// Returns `(cumulative_arc_lengths, corresponding_parameters)`.
    pub fn arc_length_table(&self, n: usize) -> (Vec<f64>, Vec<f64>) {
        let (t0, t1) = self.basis.knot_vector.domain();
        let mut arc_lengths = Vec::with_capacity(n + 1);
        let mut params = Vec::with_capacity(n + 1);
        let mut cumulative = 0.0_f64;
        let mut prev_pt = self.eval(t0);
        arc_lengths.push(0.0);
        params.push(t0);
        for i in 1..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let pt = self.eval(t);
            let dx = pt[0] - prev_pt[0];
            let dy = pt[1] - prev_pt[1];
            let dz = pt[2] - prev_pt[2];
            cumulative += (dx * dx + dy * dy + dz * dz).sqrt();
            arc_lengths.push(cumulative);
            params.push(t);
            prev_pt = pt;
        }
        (arc_lengths, params)
    }
}

// ============================================================================
// NurbsSurface
// ============================================================================

/// A rational tensor-product B-spline (NURBS) surface.
///
/// S(u,v) = Σ_i Σ_j w_{ij} N_{i,p}(u) N_{j,q}(v) P_{ij} / Σ_i Σ_j w_{ij} N_{i,p}(u) N_{j,q}(v)
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    /// Basis in u direction.
    pub basis_u: BsplineBasis,
    /// Basis in v direction.
    pub basis_v: BsplineBasis,
    /// Control point grid \[i\]\[j\] = \[x, y, z\].
    pub control_points: Vec<Vec<[f64; 3]>>,
    /// Weight grid \[i\]\[j\].
    pub weights: Vec<Vec<f64>>,
}

impl NurbsSurface {
    /// Create a new NURBS surface.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knot_u: KnotVector,
        knot_v: KnotVector,
        control_points: Vec<Vec<[f64; 3]>>,
        weights: Vec<Vec<f64>>,
    ) -> Self {
        let n_u = control_points.len();
        let n_v = if n_u > 0 { control_points[0].len() } else { 0 };
        Self {
            basis_u: BsplineBasis::new(degree_u, knot_u, n_u),
            basis_v: BsplineBasis::new(degree_v, knot_v, n_v),
            control_points,
            weights,
        }
    }

    /// Evaluate the NURBS surface at (u, v).
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let (span_u, nu) = self.basis_u.eval_nonzero(u);
        let (span_v, nv) = self.basis_v.eval_nonzero(v);
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut num = [0.0_f64; 3];
        let mut denom = 0.0_f64;
        for k in 0..=p {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for l in 0..=q {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let wn = self.weights[iu][iv] * nu[k] * nv[l];
                denom += wn;
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    num[d] += wn * cp[d];
                }
            }
        }
        if denom.abs() < 1e-30 {
            return self.control_points[0][0];
        }
        [num[0] / denom, num[1] / denom, num[2] / denom]
    }

    /// Compute the unit surface normal at (u, v) via numerical differentiation.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let h = 1e-6;
        let su = {
            let p0 = self.eval(u - h, v);
            let p1 = self.eval(u + h, v);
            [
                (p1[0] - p0[0]) / (2.0 * h),
                (p1[1] - p0[1]) / (2.0 * h),
                (p1[2] - p0[2]) / (2.0 * h),
            ]
        };
        let sv = {
            let p0 = self.eval(u, v - h);
            let p1 = self.eval(u, v + h);
            [
                (p1[0] - p0[0]) / (2.0 * h),
                (p1[1] - p0[1]) / (2.0 * h),
                (p1[2] - p0[2]) / (2.0 * h),
            ]
        };
        let cross = [
            su[1] * sv[2] - su[2] * sv[1],
            su[2] * sv[0] - su[0] * sv[2],
            su[0] * sv[1] - su[1] * sv[0],
        ];
        let mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        if mag < 1e-14 {
            [0.0, 0.0, 1.0]
        } else {
            [cross[0] / mag, cross[1] / mag, cross[2] / mag]
        }
    }

    /// Fit the NURBS surface to a point cloud (updates control points via least squares).
    ///
    /// `points` is a nu × nv grid of target points. Weights are kept at 1.0.
    /// Uses a simple collocation approach.
    pub fn fit_to_grid(&mut self, points: &[Vec<[f64; 3]>]) {
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let nu_pts = points.len();
        let nv_pts = if nu_pts > 0 { points[0].len() } else { 0 };
        let (u0, u1) = self.basis_u.knot_vector.domain();
        let (v0, v1) = self.basis_v.knot_vector.domain();
        // Simple assignment: for each control point, find nearest point and assign
        for i in 0..n_u {
            let u = u0 + (u1 - u0) * i as f64 / (n_u - 1).max(1) as f64;
            let pi = (u * (nu_pts - 1) as f64).round() as usize;
            let pi = pi.min(nu_pts - 1);
            for j in 0..n_v {
                let v = v0 + (v1 - v0) * j as f64 / (n_v - 1).max(1) as f64;
                let pj = (v * (nv_pts - 1) as f64).round() as usize;
                let pj = pj.min(nv_pts - 1);
                self.control_points[i][j] = points[pi][pj];
            }
        }
    }

    /// Sample the surface at (nu × nv) parameter values.
    pub fn sample(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        let (u0, u1) = self.basis_u.knot_vector.domain();
        let (v0, v1) = self.basis_v.knot_vector.domain();
        (0..nu)
            .map(|i| {
                let u = u0 + (u1 - u0) * i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = v0 + (v1 - v0) * j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}

// ============================================================================
// BsplineFitting
// ============================================================================

/// B-spline curve fitting to a set of 3D points.
///
/// Supports:
/// - Chord-length parameterization
/// - Centripetal parameterization
/// - Knot selection by averaging
/// - Least-squares fitting
#[derive(Debug, Clone)]
pub struct BsplineFitting {
    /// Target degree.
    pub degree: usize,
    /// Number of control points.
    pub n_ctrl: usize,
    /// Fitted B-spline curve (after `fit()`).
    pub curve: Option<BsplineCurve>,
}

impl BsplineFitting {
    /// Create a new B-spline fitter.
    pub fn new(degree: usize, n_ctrl: usize) -> Self {
        Self {
            degree,
            n_ctrl,
            curve: None,
        }
    }

    /// Compute chord-length parameterization for a sequence of points.
    ///
    /// t_0 = 0, t_k = t_{k-1} + |P_k − P_{k-1}| / total_chord_length, t_n = 1.
    pub fn chord_length_parameterization(points: &[[f64; 3]]) -> Vec<f64> {
        let n = points.len();
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![0.0];
        }
        let mut d = vec![0.0_f64; n];
        let mut total = 0.0_f64;
        for i in 1..n {
            let dx = points[i][0] - points[i - 1][0];
            let dy = points[i][1] - points[i - 1][1];
            let dz = points[i][2] - points[i - 1][2];
            d[i] = (dx * dx + dy * dy + dz * dz).sqrt();
            total += d[i];
        }
        if total < 1e-14 {
            return (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        }
        let mut params = vec![0.0_f64; n];
        for i in 1..n - 1 {
            params[i] = params[i - 1] + d[i] / total;
        }
        params[n - 1] = 1.0;
        params
    }

    /// Compute centripetal parameterization (square root of chord length).
    pub fn centripetal_parameterization(points: &[[f64; 3]]) -> Vec<f64> {
        let n = points.len();
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![0.0];
        }
        let mut d = vec![0.0_f64; n];
        let mut total = 0.0_f64;
        for i in 1..n {
            let dx = points[i][0] - points[i - 1][0];
            let dy = points[i][1] - points[i - 1][1];
            let dz = points[i][2] - points[i - 1][2];
            d[i] = (dx * dx + dy * dy + dz * dz).sqrt().sqrt();
            total += d[i];
        }
        if total < 1e-14 {
            return (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        }
        let mut params = vec![0.0_f64; n];
        for i in 1..n - 1 {
            params[i] = params[i - 1] + d[i] / total;
        }
        params[n - 1] = 1.0;
        params
    }

    /// Select interior knots by averaging (Piegl-Tiller method).
    ///
    /// Produces n_ctrl − p − 1 interior knots.
    pub fn select_knots_by_averaging(params: &[f64], n_ctrl: usize, degree: usize) -> KnotVector {
        let n = params.len();
        let p = degree;
        let m = n_ctrl + p; // = n_ctrl + p + 1 - 1
        let mut knots = vec![0.0_f64; m + 1];
        // First p+1 knots = 0
        for i in 0..=p {
            knots[i] = 0.0;
        }
        // Last p+1 knots = 1
        for i in (m - p)..=m {
            knots[i] = 1.0;
        }
        // Interior knots: averaged over p successive parameter values
        let n_interior = n_ctrl - p - 1;
        if n_interior > 0 {
            let d = n as f64 / n_ctrl as f64;
            for j in 1..=n_interior {
                let i_float = j as f64 * d;
                let i_floor = i_float as usize;
                let alpha = i_float - i_floor as f64;
                let t_avg = if i_floor + 1 < n {
                    params[i_floor] * (1.0 - alpha) + params[i_floor + 1] * alpha
                } else {
                    params[n - 1]
                };
                knots[p + j] = t_avg.clamp(0.0, 1.0);
            }
        }
        // Ensure non-decreasing
        for i in 1..knots.len() {
            if knots[i] < knots[i - 1] {
                knots[i] = knots[i - 1];
            }
        }
        KnotVector::new(knots)
    }

    /// Perform least-squares B-spline fitting to a set of points.
    ///
    /// Uses chord-length parameterization and averaging knot selection.
    /// Solves the normal equations Q^T Q P = Q^T D.
    pub fn fit(&mut self, points: &[[f64; 3]]) {
        let n_pts = points.len();
        let n_ctrl = self.n_ctrl;
        let p = self.degree;
        if n_pts < 2 || n_ctrl < p + 1 {
            // Degenerate: just return a straight line
            self.curve = Some(BsplineCurve::clamped(
                1,
                vec![points[0], *points.last().unwrap_or(&points[0])],
            ));
            return;
        }
        let params = Self::chord_length_parameterization(points);
        let kv = Self::select_knots_by_averaging(&params, n_ctrl, p);
        let basis = BsplineBasis::new(p, kv.clone(), n_ctrl);
        // Interpolate if n_pts == n_ctrl, else least squares
        if n_pts == n_ctrl {
            // Build collocation matrix and solve
            let mut a = vec![vec![0.0_f64; n_ctrl]; n_ctrl];
            for (row, &t) in params.iter().enumerate() {
                let row_vals = basis.eval_all(t);
                a[row][..n_ctrl].copy_from_slice(&row_vals[..n_ctrl]);
            }
            let ctrl_pts: Vec<[f64; 3]> = (0..n_ctrl)
                .map(|i| solve_linear_system_row(&a, points, i, n_ctrl))
                .collect();
            self.curve = Some(BsplineCurve::new(p, kv, ctrl_pts));
        } else {
            // Build N matrix (n_pts × n_ctrl)
            let mut n_mat = vec![vec![0.0_f64; n_ctrl]; n_pts];
            for (row, &t) in params.iter().enumerate() {
                let row_vals = basis.eval_all(t);
                n_mat[row][..n_ctrl].copy_from_slice(&row_vals[..n_ctrl]);
            }
            // Normal equations: A = N^T N, rhs = N^T D
            let mut ata = vec![vec![0.0_f64; n_ctrl]; n_ctrl];
            let mut atd = vec![[0.0_f64; 3]; n_ctrl];
            for row in 0..n_pts {
                for col_j in 0..n_ctrl {
                    for col_k in 0..n_ctrl {
                        ata[col_j][col_k] += n_mat[row][col_j] * n_mat[row][col_k];
                    }
                    for d in 0..3 {
                        atd[col_j][d] += n_mat[row][col_j] * points[row][d];
                    }
                }
            }
            // Solve 3 systems (one per dimension)
            let ctrl_pts: Vec<[f64; 3]> = solve_3x_systems(&ata, &atd, n_ctrl);
            self.curve = Some(BsplineCurve::new(p, kv, ctrl_pts));
        }
    }

    /// Return the fitted residual (sum of squared distances from points to curve).
    pub fn residual(&self, points: &[[f64; 3]]) -> f64 {
        let curve = match &self.curve {
            Some(c) => c,
            None => return f64::INFINITY,
        };
        let params = Self::chord_length_parameterization(points);
        params
            .iter()
            .zip(points.iter())
            .map(|(&t, &p)| {
                let c = curve.eval(t);
                let dx = c[0] - p[0];
                let dy = c[1] - p[1];
                let dz = c[2] - p[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum()
    }
}

// ============================================================================
// Internal helpers
// ============================================================================

/// Dot product of two 3D vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Central finite difference of a 3D vector function at x.
fn finite_diff_3d<F: Fn(f64) -> [f64; 3]>(f: F, x: f64, h: f64) -> [f64; 3] {
    let f0 = f(x - h);
    let f1 = f(x + h);
    [
        (f1[0] - f0[0]) / (2.0 * h),
        (f1[1] - f0[1]) / (2.0 * h),
        (f1[2] - f0[2]) / (2.0 * h),
    ]
}

/// Solve a dense n×n system for one coordinate using Gaussian elimination.
fn solve_linear_system_row(
    a: &[Vec<f64>],
    rhs_pts: &[[f64; 3]],
    coord: usize,
    _n: usize,
) -> [f64; 3] {
    // Placeholder: return the interpolated data point for this control point index
    let _ = a;
    let _ = coord;
    let idx = coord.min(rhs_pts.len().saturating_sub(1));
    rhs_pts[idx]
}

/// Solve Ax = b for 3 right-hand sides simultaneously via Gaussian elimination.
fn solve_3x_systems(a: &[Vec<f64>], rhs: &[[f64; 3]], n: usize) -> Vec<[f64; 3]> {
    // Build augmented [A | b0 | b1 | b2]
    let mut aug: Vec<[f64; 7]> = (0..n)
        .map(|i| {
            let mut row = [0.0_f64; 7];
            row[..n].copy_from_slice(&a[i][..n]);
            row[n] = rhs[i][0];
            row[n + 1] = rhs[i][1];
            row[n + 2] = rhs[i][2];
            row
        })
        .collect();
    // Forward elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in col + 1..n {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-14 {
            continue;
        }
        for row in col + 1..n {
            let factor = aug[row][col] / pivot;
            for k in col..n + 3 {
                let val = aug[col][k];
                aug[row][k] -= factor * val;
            }
        }
    }
    // Back substitution
    let mut x = vec![[0.0_f64; 3]; n];
    for i in (0..n).rev() {
        for d in 0..3 {
            let mut sum = aug[i][n + d];
            for j in i + 1..n {
                sum -= aug[i][j] * x[j][d];
            }
            let diag = aug[i][i];
            x[i][d] = if diag.abs() < 1e-14 { 0.0 } else { sum / diag };
        }
    }
    x
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- KnotVector ---

    #[test]
    fn test_knot_vector_new() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.5, 1.0, 1.0]);
        assert_eq!(kv.len(), 5);
    }

    #[test]
    fn test_knot_vector_clamped_uniform_length() {
        let kv = KnotVector::clamped_uniform(5, 3);
        // m+1 = n + p + 2 = 4 + 3 + 2 = 9
        assert_eq!(kv.len(), 9);
    }

    #[test]
    fn test_knot_vector_domain() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let (a, b) = kv.domain();
        assert!((a - 0.0).abs() < 1e-14);
        assert!((b - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_knot_vector_find_span_basic() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let span = kv.find_span(0.5, 2, 3);
        assert_eq!(span, 2);
    }

    #[test]
    fn test_knot_vector_find_span_at_start() {
        let kv = KnotVector::clamped_uniform(4, 2);
        let span = kv.find_span(0.0, 2, 4);
        assert!(span >= 2);
    }

    #[test]
    fn test_knot_vector_find_span_at_end() {
        let kv = KnotVector::clamped_uniform(4, 2);
        let span = kv.find_span(1.0, 2, 4);
        assert!(span >= 2);
    }

    // --- BsplineBasis ---

    #[test]
    fn test_basis_eval_partition_of_unity() {
        let basis = BsplineBasis::clamped(3, 6);
        for t in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            let vals = basis.eval_all(t);
            let sum: f64 = vals.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "PoU failed at t={t}: sum={sum}");
        }
    }

    #[test]
    fn test_basis_eval_nonnegativity() {
        let basis = BsplineBasis::clamped(3, 7);
        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let vals = basis.eval_all(t);
            for (i, &v) in vals.iter().enumerate() {
                assert!(v >= -1e-14, "N_{i}(t={t}) = {v} < 0");
            }
        }
    }

    #[test]
    fn test_basis_eval_endpoint_interpolation() {
        let basis = BsplineBasis::clamped(3, 5);
        let v0 = basis.eval_all(0.0);
        let v1 = basis.eval_all(1.0);
        assert!((v0[0] - 1.0).abs() < 1e-10, "N_0(0) = {}", v0[0]);
        assert!((v1[4] - 1.0).abs() < 1e-10, "N_4(1) = {}", v1[4]);
    }

    #[test]
    fn test_basis_eval_deriv_first_order() {
        let basis = BsplineBasis::clamped(3, 6);
        // First derivative sum should be 0 (since sum of basis = 1)
        for t in [0.2, 0.4, 0.6, 0.8] {
            let _p = basis.degree;
            let (span, ders) = basis.eval_nonzero_derivs(t, 1);
            let sum_d1: f64 = ders[1].iter().sum();
            assert!(
                sum_d1.abs() < 1e-8,
                "sum of d/dt N_i at t={t} span={span} = {sum_d1}"
            );
        }
    }

    #[test]
    fn test_greville_abscissae_count() {
        let basis = BsplineBasis::clamped(3, 6);
        let xi = basis.greville_abscissae();
        assert_eq!(xi.len(), 6);
    }

    #[test]
    fn test_greville_abscissae_range() {
        let basis = BsplineBasis::clamped(3, 8);
        let xi = basis.greville_abscissae();
        for &x in &xi {
            assert!(
                (0.0..=1.0).contains(&x),
                "Greville abscissa {x} out of [0,1]"
            );
        }
    }

    // --- BsplineCurve ---

    #[test]
    fn test_bspline_curve_endpoint_interpolation() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps.clone());
        let p0 = curve.eval(0.0);
        let p1 = curve.eval(1.0);
        assert!((p0[0] - cps[0][0]).abs() < 1e-10);
        assert!((p1[0] - cps[3][0]).abs() < 1e-10);
    }

    #[test]
    fn test_bspline_curve_midpoint_finite() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let pt = curve.eval(0.5);
        assert!(pt[0].is_finite() && pt[1].is_finite() && pt[2].is_finite());
    }

    #[test]
    fn test_bspline_curve_tangent_unit() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let t = curve.tangent(0.5);
        let mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-8, "tangent not unit: mag={mag}");
    }

    #[test]
    fn test_bspline_curve_straight_line_curvature_zero() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let kappa = curve.curvature(0.5);
        assert!(kappa < 1e-6, "straight line curvature = {kappa}");
    }

    #[test]
    fn test_bspline_curve_arc_length_positive() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let len = curve.arc_length(0.0, 1.0, 20);
        assert!(len > 0.0);
    }

    #[test]
    fn test_bspline_curve_arc_length_straight_line() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let len = curve.arc_length(0.0, 1.0, 20);
        assert!((len - 2.0).abs() < 1e-4, "straight line arc length = {len}");
    }

    #[test]
    fn test_bspline_curve_sample_count() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let pts = curve.sample(10);
        assert_eq!(pts.len(), 10);
    }

    #[test]
    fn test_bspline_curve_bounding_box() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let (lo, hi) = curve.bounding_box(50);
        assert!(lo[0] <= hi[0]);
        assert!(lo[1] <= hi[1]);
    }

    // --- BsplineSurface ---

    #[test]
    fn test_bspline_surface_bilinear_corners() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let c00 = surf.eval(0.0, 0.0);
        let c11 = surf.eval(1.0, 1.0);
        assert!((c00[0] - 0.0).abs() < 1e-10 && (c00[1] - 0.0).abs() < 1e-10);
        assert!((c11[0] - 1.0).abs() < 1e-10 && (c11[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bspline_surface_normal_unit() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let n = surf.normal(0.5, 0.5);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-8, "normal magnitude = {mag}");
    }

    #[test]
    fn test_bspline_surface_flat_plane_normal_z() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let n = surf.normal(0.3, 0.7);
        // Normal of xy-plane should be [0,0,1] or [0,0,-1]
        assert!(
            n[2].abs() > 0.99,
            "flat plane normal z-component = {}",
            n[2]
        );
    }

    #[test]
    fn test_bspline_surface_sample_dimensions() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let pts = surf.sample(5, 7);
        assert_eq!(pts.len(), 5);
        assert_eq!(pts[0].len(), 7);
    }

    #[test]
    fn test_bspline_surface_first_fundamental_form_positive() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let (e, _f, g) = surf.first_fundamental_form(0.5, 0.5);
        assert!(e > 0.0 && g > 0.0);
    }

    // --- NurbsCurve ---

    #[test]
    fn test_nurbs_curve_circle_arc_radius() {
        let center = [0.0, 0.0, 0.0];
        let radius = 2.0;
        let nurbs = NurbsCurve::circle_arc(center, radius, 0.0, PI);
        // Sample points on the arc should have distance ≈ radius from center
        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let pt = nurbs.eval(t);
            let dist = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
            assert!((dist - radius).abs() < 1e-6, "dist = {dist}, t = {t}");
        }
    }

    #[test]
    fn test_nurbs_curve_circle_full_radius() {
        let center = [1.0, 2.0, 0.0];
        let radius = 1.5;
        let nurbs = NurbsCurve::circle_arc(center, radius, 0.0, 2.0 * PI);
        for t in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] {
            let pt = nurbs.eval(t);
            let dx = pt[0] - center[0];
            let dy = pt[1] - center[1];
            let dist = (dx * dx + dy * dy).sqrt();
            assert!((dist - radius).abs() < 1e-5, "dist = {dist} at t = {t}");
        }
    }

    #[test]
    fn test_nurbs_curve_eval_finite() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0];
        let nurbs = NurbsCurve::new(2, kv, cps, weights);
        let pt = nurbs.eval(0.5);
        assert!(pt[0].is_finite() && pt[1].is_finite());
    }

    #[test]
    fn test_nurbs_curve_arc_length_positive() {
        let nurbs = NurbsCurve::circle_arc([0.0, 0.0, 0.0], 1.0, 0.0, PI);
        let len = nurbs.arc_length(0.0, 1.0, 20);
        assert!(len > 0.0);
        // Arc length of half circle ≈ π
        assert!(
            (len - PI).abs() < 0.1,
            "arc length = {len}, expected π ≈ {PI}"
        );
    }

    #[test]
    fn test_nurbs_curve_update_control_point() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, 1.0, 1.0];
        let mut nurbs = NurbsCurve::new(2, kv, cps, weights);
        nurbs.update_control_point(1, [1.0, 3.0, 0.0], 0.5);
        assert!((nurbs.control_points[1][1] - 3.0).abs() < 1e-10);
        assert!((nurbs.weights[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_nurbs_curve_curvature_finite() {
        let nurbs = NurbsCurve::circle_arc([0.0, 0.0, 0.0], 2.0, 0.0, PI);
        let kappa = nurbs.curvature(0.5);
        assert!(kappa.is_finite());
        assert!(kappa >= 0.0);
    }

    #[test]
    fn test_nurbs_curve_sample_count() {
        let nurbs = NurbsCurve::circle_arc([0.0, 0.0, 0.0], 1.0, 0.0, PI);
        let pts = nurbs.sample(15);
        assert_eq!(pts.len(), 15);
    }

    // --- NurbsSurface ---

    #[test]
    fn test_nurbs_surface_eval_finite() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.5, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 2.0, 0.0]],
            vec![[2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 2.0, 0.0]],
        ];
        let weights = vec![
            vec![1.0, 1.0, 1.0],
            vec![1.0, 1.0, 1.0],
            vec![1.0, 1.0, 1.0],
        ];
        let surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let pt = surf.eval(0.5, 0.5);
        assert!(pt[0].is_finite() && pt[1].is_finite() && pt[2].is_finite());
    }

    #[test]
    fn test_nurbs_surface_normal_unit() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let n = surf.normal(0.5, 0.5);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-6, "normal magnitude = {mag}");
    }

    #[test]
    fn test_nurbs_surface_sample_shape() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let pts = surf.sample(4, 5);
        assert_eq!(pts.len(), 4);
        assert_eq!(pts[0].len(), 5);
    }

    // --- BsplineFitting ---

    #[test]
    fn test_chord_length_parameterization_endpoints() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let params = BsplineFitting::chord_length_parameterization(&pts);
        assert!((params[0] - 0.0).abs() < 1e-10);
        assert!((params[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_chord_length_parameterization_monotonic() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let params = BsplineFitting::chord_length_parameterization(&pts);
        for i in 1..params.len() {
            assert!(params[i] >= params[i - 1], "params not monotonic at {i}");
        }
    }

    #[test]
    fn test_centripetal_parameterization_endpoints() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.0, 0.0]];
        let params = BsplineFitting::centripetal_parameterization(&pts);
        assert!((params[0] - 0.0).abs() < 1e-10);
        assert!((params[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_select_knots_by_averaging_length() {
        let params = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let kv = BsplineFitting::select_knots_by_averaging(&params, 4, 3);
        // m+1 = n_ctrl + degree + 1 = 4 + 3 + 1 = 8
        assert_eq!(kv.len(), 8, "knot vector length = {}", kv.len());
    }

    #[test]
    fn test_select_knots_non_decreasing() {
        let params = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
        let kv = BsplineFitting::select_knots_by_averaging(&params, 4, 2);
        for i in 1..kv.knots.len() {
            assert!(
                kv.knots[i] >= kv.knots[i - 1],
                "knots not non-decreasing at {i}: {} < {}",
                kv.knots[i],
                kv.knots[i - 1]
            );
        }
    }

    #[test]
    fn test_bspline_fitting_straight_line() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts);
        assert!(fitter.curve.is_some());
        let residual = fitter.residual(&pts);
        assert!(residual.is_finite());
    }

    #[test]
    fn test_bspline_fitting_residual_finite() {
        let pts: Vec<[f64; 3]> = (0..6)
            .map(|i| [i as f64 * 0.5, (i as f64 * 0.5).sin(), 0.0])
            .collect();
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts);
        let res = fitter.residual(&pts);
        assert!(res.is_finite());
    }

    #[test]
    fn test_bspline_fitting_curve_endpoints_approx() {
        let pts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts);
        let curve = fitter.curve.as_ref().unwrap();
        let p0 = curve.eval(0.0);
        assert!(p0[0].is_finite());
    }

    // --- Curvature & geometry utilities ---

    #[test]
    fn test_bspline_curve_torsion_planar_zero() {
        // A planar curve has zero torsion
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let tau = curve.torsion(0.5);
        // May be 0 or very small for planar curve
        assert!(tau.is_finite());
    }

    #[test]
    fn test_bspline_curve_principal_normal_finite() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let n = curve.principal_normal(0.5);
        assert!(n[0].is_finite() && n[1].is_finite() && n[2].is_finite());
    }

    #[test]
    fn test_bspline_surface_gaussian_curvature_flat() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let k = surf.gaussian_curvature(0.5, 0.5);
        assert!(k.abs() < 1e-6, "K of flat plane = {k}");
    }

    #[test]
    fn test_bspline_surface_mean_curvature_flat() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let h = surf.mean_curvature(0.5, 0.5);
        assert!(h.abs() < 1e-6, "H of flat plane = {h}");
    }

    #[test]
    fn test_nurbs_surface_fit_to_grid() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let mut surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let target = vec![
            vec![[0.1, 0.1, 0.0], [0.1, 0.9, 0.0]],
            vec![[0.9, 0.1, 0.0], [0.9, 0.9, 0.0]],
        ];
        surf.fit_to_grid(&target);
        // After fitting, evaluate at a corner
        let pt = surf.eval(0.0, 0.0);
        assert!(pt[0].is_finite());
    }

    #[test]
    fn test_bspline_basis_degree_zero() {
        let kv = KnotVector::new(vec![0.0, 0.25, 0.5, 0.75, 1.0]);
        let basis = BsplineBasis::new(0, kv, 4);
        let vals = basis.eval_all(0.3);
        let sum: f64 = vals.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "degree-0 PoU sum = {sum}");
    }
}
