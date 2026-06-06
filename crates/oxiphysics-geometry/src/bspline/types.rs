//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{solve_3x_systems, solve_linear_system_row};

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
        ctrl_pts.push([
            center[0] + radius * angle.cos(),
            center[1] + radius * angle.sin(),
            center[2],
        ]);
        weights.push(1.0);
        for _i in 0..n_arcs {
            let a_mid = angle + d_theta * 0.5;
            let a_end = angle + d_theta;
            ctrl_pts.push([
                center[0] + radius / w1 * a_mid.cos(),
                center[1] + radius / w1 * a_mid.sin(),
                center[2],
            ]);
            weights.push(w1);
            ctrl_pts.push([
                center[0] + radius * a_end.cos(),
                center[1] + radius * a_end.sin(),
                center[2],
            ]);
            weights.push(1.0);
            angle = a_end;
        }
        let mut knots = vec![0.0_f64; 3];
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
        for (k, &n_val) in nonzero[..=p].iter().enumerate() {
            let idx = span + k - p;
            if idx >= n_ctrl {
                continue;
            }
            let wn = self.weights[idx] * n_val;
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
            let at_end = (t - knots[knots.len() - 1]).abs() < 1e-14
                && (knots[i + 1] - knots[knots.len() - 1]).abs() < 1e-14;
            return if (knots[i] <= t && t < knots[i + 1]) || at_end {
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
        let w = std::f64::consts::FRAC_1_SQRT_2;
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
        for (k, ders_row) in ders[1..=d].iter_mut().enumerate().map(|(i, v)| (i + 1, v)) {
            for val in ders_row[..=p].iter_mut() {
                *val *= rfact;
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
        let m = n_ctrl + p;
        let mut knots = vec![0.0_f64; m + 1];
        for val in knots[..=p].iter_mut() {
            *val = 0.0;
        }
        for val in knots[(m - p)..=m].iter_mut() {
            *val = 1.0;
        }
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
            self.curve = Some(BsplineCurve::clamped(
                1,
                vec![points[0], *points.last().unwrap_or(&points[0])],
            ));
            return;
        }
        let params = Self::chord_length_parameterization(points);
        let kv = Self::select_knots_by_averaging(&params, n_ctrl, p);
        let basis = BsplineBasis::new(p, kv.clone(), n_ctrl);
        if n_pts == n_ctrl {
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
            let mut n_mat = vec![vec![0.0_f64; n_ctrl]; n_pts];
            for (row, &t) in params.iter().enumerate() {
                let row_vals = basis.eval_all(t);
                n_mat[row][..n_ctrl].copy_from_slice(&row_vals[..n_ctrl]);
            }
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
        for (k, &n_val) in nonzero[..=p].iter().enumerate() {
            let idx = span + k - p;
            if idx < self.control_points.len() {
                let cp = self.control_points[idx];
                for (pt, &cp_d) in point.iter_mut().zip(cp.iter()) {
                    *pt += n_val * cp_d;
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
        for (j, &d_val) in ders[k][..=p].iter().enumerate() {
            let idx = span + j - p;
            if idx < n_ctrl {
                let cp = self.control_points[idx];
                for (res_d, &cp_d) in result.iter_mut().zip(cp.iter()) {
                    *res_d += d_val * cp_d;
                }
            }
        }
        result
    }
    /// Compute arc length from t_start to t_end using Gaussian quadrature (n_gauss points per interval).
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
        for val in knots[..=p].iter_mut() {
            *val = 0.0;
        }
        for val in knots[(m - p)..=m].iter_mut() {
            *val = 1.0;
        }
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
        if t >= knots[n + 1] {
            let mut i = n;
            while i > p && (knots[i] - knots[i + 1]).abs() < 1e-14 {
                i -= 1;
            }
            return i;
        }
        if t <= knots[p] {
            return p;
        }
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
