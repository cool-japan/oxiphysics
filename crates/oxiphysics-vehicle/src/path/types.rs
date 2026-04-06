//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// An ordered sequence of waypoints describing a vehicle path.
#[derive(Debug, Clone)]
pub struct Path {
    /// Ordered waypoints that define the path.
    pub waypoints: Vec<Waypoint>,
    /// When `true` the path wraps from the last waypoint back to the first.
    pub closed: bool,
}
impl Path {
    /// Create an empty path.
    pub fn new(closed: bool) -> Self {
        Self {
            waypoints: Vec::new(),
            closed,
        }
    }
    /// Append a waypoint.
    pub fn add_waypoint(&mut self, pos: [f64; 3], speed: f64) {
        self.waypoints.push(Waypoint {
            position: pos,
            speed_limit: speed,
            heading: 0.0,
        });
    }
    /// Euclidean distance between two 3-D points.
    fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let dz = b[2] - a[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Sum of all inter-waypoint segment lengths.
    pub fn total_length(&self) -> f64 {
        let n = self.waypoints.len();
        if n < 2 {
            return 0.0;
        }
        let mut len = 0.0;
        for i in 0..n - 1 {
            len += Self::dist3(self.waypoints[i].position, self.waypoints[i + 1].position);
        }
        if self.closed {
            len += Self::dist3(self.waypoints[n - 1].position, self.waypoints[0].position);
        }
        len
    }
    /// Return the index of the waypoint nearest to `pos` and its distance.
    pub fn nearest_waypoint(&self, pos: [f64; 3]) -> (usize, f64) {
        let mut best_idx = 0;
        let mut best_dist = f64::MAX;
        for (i, wp) in self.waypoints.iter().enumerate() {
            let d = Self::dist3(pos, wp.position);
            if d < best_dist {
                best_dist = d;
                best_idx = i;
            }
        }
        if self.waypoints.is_empty() {
            (0, 0.0)
        } else {
            (best_idx, best_dist)
        }
    }
    /// Return the world position at arc-length `d` measured from the first waypoint.
    pub fn interpolate_at_distance(&self, d: f64) -> Option<[f64; 3]> {
        let n = self.waypoints.len();
        if n == 0 {
            return None;
        }
        if n == 1 {
            return Some(self.waypoints[0].position);
        }
        if d <= 0.0 {
            return Some(self.waypoints[0].position);
        }
        let mut accumulated = 0.0;
        let segments = if self.closed { n } else { n - 1 };
        for i in 0..segments {
            let a = self.waypoints[i].position;
            let b = self.waypoints[(i + 1) % n].position;
            let seg_len = Self::dist3(a, b);
            if accumulated + seg_len >= d {
                let t = (d - accumulated) / seg_len;
                return Some([
                    a[0] + t * (b[0] - a[0]),
                    a[1] + t * (b[1] - a[1]),
                    a[2] + t * (b[2] - a[2]),
                ]);
            }
            accumulated += seg_len;
        }
        Some(self.waypoints[n - 1].position)
    }
    /// Return the lookahead point on the path for the pure-pursuit algorithm.
    pub fn lookahead_point(&self, vehicle_pos: [f64; 3], lookahead_dist: f64) -> Option<[f64; 3]> {
        let n = self.waypoints.len();
        if n < 2 {
            return None;
        }
        let (start_idx, _) = self.nearest_waypoint(vehicle_pos);
        let segments = if self.closed { n } else { n - 1 };
        for step in 0..segments {
            let i = (start_idx + step) % n;
            let a = self.waypoints[i].position;
            let b = self.waypoints[(i + 1) % n].position;
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];
            let dz = b[2] - a[2];
            let fx = a[0] - vehicle_pos[0];
            let fy = a[1] - vehicle_pos[1];
            let fz = a[2] - vehicle_pos[2];
            let aa = dx * dx + dy * dy + dz * dz;
            let bb = 2.0 * (fx * dx + fy * dy + fz * dz);
            let cc = fx * fx + fy * fy + fz * fz - lookahead_dist * lookahead_dist;
            let discriminant = bb * bb - 4.0 * aa * cc;
            if discriminant < 0.0 {
                continue;
            }
            let sqrt_disc = discriminant.sqrt();
            let t2 = (-bb + sqrt_disc) / (2.0 * aa);
            let t1 = (-bb - sqrt_disc) / (2.0 * aa);
            let t = if (0.0..=1.0).contains(&t2) {
                t2
            } else if (0.0..=1.0).contains(&t1) {
                t1
            } else {
                continue;
            };
            return Some([a[0] + t * dx, a[1] + t * dy, a[2] + t * dz]);
        }
        Some(self.waypoints[(start_idx + segments) % n].position)
    }
    /// Compute cumulative arc lengths from the first waypoint.
    ///
    /// Returns a vector of the same length as `waypoints` where entry `i` is the
    /// arc-length distance from waypoint 0 to waypoint `i`.
    pub fn cumulative_arc_lengths(&self) -> Vec<f64> {
        let n = self.waypoints.len();
        if n == 0 {
            return vec![];
        }
        let mut arcs = vec![0.0; n];
        for i in 1..n {
            arcs[i] = arcs[i - 1]
                + Self::dist3(self.waypoints[i - 1].position, self.waypoints[i].position);
        }
        arcs
    }
    /// Compute the heading at each waypoint from finite differences.
    ///
    /// The heading is the angle (radians from +X) of the vector from this
    /// waypoint to the next.
    pub fn compute_headings(&mut self) {
        let n = self.waypoints.len();
        if n < 2 {
            return;
        }
        for i in 0..n {
            let next = if i + 1 < n {
                i + 1
            } else if self.closed {
                0
            } else {
                i
            };
            if next == i {
                if i > 0 {
                    self.waypoints[i].heading = self.waypoints[i - 1].heading;
                }
            } else {
                let dx = self.waypoints[next].position[0] - self.waypoints[i].position[0];
                let dy = self.waypoints[next].position[1] - self.waypoints[i].position[1];
                self.waypoints[i].heading = dy.atan2(dx);
            }
        }
    }
}
/// Pure pursuit lateral controller.
#[derive(Debug, Clone)]
pub struct PurePursuitController {
    /// Distance ahead of the vehicle at which the target point is chosen (m).
    pub lookahead_distance: f64,
    /// Distance between front and rear axles (m).
    pub wheelbase: f64,
}
impl PurePursuitController {
    /// Create a new controller.
    pub fn new(lookahead_dist: f64, wheelbase: f64) -> Self {
        Self {
            lookahead_distance: lookahead_dist,
            wheelbase,
        }
    }
    /// Compute the required front-wheel steering angle (radians).
    pub fn steering_angle(
        &self,
        vehicle_pos: [f64; 3],
        vehicle_heading: f64,
        target: [f64; 3],
    ) -> f64 {
        let alpha = self.alpha(vehicle_pos, vehicle_heading, target);
        f64::atan2(2.0 * self.wheelbase * alpha.sin(), self.lookahead_distance)
    }
    /// Compute path curvature kappa = 2*sin(alpha) / l_d.
    pub fn curvature(&self, vehicle_pos: [f64; 3], vehicle_heading: f64, target: [f64; 3]) -> f64 {
        let alpha = self.alpha(vehicle_pos, vehicle_heading, target);
        2.0 * alpha.sin() / self.lookahead_distance
    }
    /// Signed angle alpha from vehicle heading to the target (radians).
    fn alpha(&self, vehicle_pos: [f64; 3], vehicle_heading: f64, target: [f64; 3]) -> f64 {
        let dx = target[0] - vehicle_pos[0];
        let dy = target[1] - vehicle_pos[1];
        let angle_to_target = dy.atan2(dx);
        let mut alpha = angle_to_target - vehicle_heading;
        while alpha > std::f64::consts::PI {
            alpha -= 2.0 * std::f64::consts::PI;
        }
        while alpha < -std::f64::consts::PI {
            alpha += 2.0 * std::f64::consts::PI;
        }
        alpha
    }
}
/// Stanley lateral controller.
#[derive(Debug, Clone)]
pub struct StanleyController {
    /// Gain on the cross-track error term.
    pub k: f64,
    /// Softening constant that prevents division-by-zero at low speeds (m/s).
    pub k_soft: f64,
}
impl StanleyController {
    /// Create a new Stanley controller with the given gain.
    pub fn new(k: f64) -> Self {
        Self { k, k_soft: 1.0 }
    }
    /// Compute the steering angle (radians).
    pub fn steering_angle(&self, cross_track_error: f64, heading_error: f64, speed: f64) -> f64 {
        heading_error + f64::atan2(self.k * cross_track_error, self.k_soft + speed)
    }
}
/// Natural cubic spline through a sequence of 2-D (X, Y) points.
///
/// Computes spline coefficients once, then allows fast evaluation at any
/// parameter `t in [0, n-1]` where integer values correspond to waypoints.
pub struct CubicSpline2D {
    /// X coordinates of control points.
    pub(super) xs: Vec<f64>,
    /// Y coordinates of control points.
    pub(super) ys: Vec<f64>,
    /// Second derivatives for X spline.
    pub(super) mx: Vec<f64>,
    /// Second derivatives for Y spline.
    pub(super) my: Vec<f64>,
}
impl CubicSpline2D {
    /// Build a natural cubic spline from the given (X, Y) waypoints.
    ///
    /// Requires at least 2 points; returns `None` otherwise.
    pub fn new(xs: &[f64], ys: &[f64]) -> Option<Self> {
        let n = xs.len();
        if n < 2 || ys.len() != n {
            return None;
        }
        let mx = Self::compute_second_derivatives(xs);
        let my = Self::compute_second_derivatives(ys);
        Some(Self {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
            mx,
            my,
        })
    }
    /// Compute second derivatives for a natural cubic spline.
    ///
    /// Natural boundary conditions: M_0 = M_{n-1} = 0.
    fn compute_second_derivatives(vals: &[f64]) -> Vec<f64> {
        let n = vals.len();
        if n < 3 {
            return vec![0.0; n];
        }
        let mut m = vec![0.0; n];
        let mut d = vec![0.0; n];
        let mut z = vec![0.0; n];
        for i in 1..n - 1 {
            let rhs = 6.0 * (vals[i + 1] - 2.0 * vals[i] + vals[i - 1]);
            let denom = 4.0 - d[i - 1];
            d[i] = 1.0 / denom;
            z[i] = (rhs - z[i - 1]) / denom;
        }
        for i in (1..n - 1).rev() {
            m[i] = z[i] - d[i] * m[i + 1];
        }
        m
    }
    /// Evaluate the spline at parameter `t`.
    ///
    /// `t = 0` corresponds to the first point; `t = n-1` to the last.
    /// Values outside `[0, n-1]` are clamped.
    pub fn evaluate(&self, t: f64) -> [f64; 2] {
        let n = self.xs.len();
        let t_clamped = t.clamp(0.0, (n - 1) as f64);
        let i = (t_clamped.floor() as usize).min(n - 2);
        let s = t_clamped - i as f64;
        let x = Self::spline_segment(&self.xs, &self.mx, i, s);
        let y = Self::spline_segment(&self.ys, &self.my, i, s);
        [x, y]
    }
    /// Evaluate one cubic segment: val\[i\]*(1-s) + val\[i+1\]*s
    ///   + (1/6)*((m\[i\]*(1-s)^3 - m\[i\]*(1-s)) + (m\[i+1\]*s^3 - m\[i+1\]*s))
    fn spline_segment(vals: &[f64], m: &[f64], i: usize, s: f64) -> f64 {
        let a = 1.0 - s;
        vals[i] * a + vals[i + 1] * s + (m[i] * (a * a * a - a) + m[i + 1] * (s * s * s - s)) / 6.0
    }
    /// Number of control points.
    pub fn len(&self) -> usize {
        self.xs.len()
    }
    /// Whether the spline has no points.
    pub fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }
    /// Compute the first derivative (tangent) at parameter `t`.
    pub fn tangent(&self, t: f64) -> [f64; 2] {
        let n = self.xs.len();
        let t_clamped = t.clamp(0.0, (n - 1) as f64);
        let i = (t_clamped.floor() as usize).min(n - 2);
        let s = t_clamped - i as f64;
        let tx = Self::spline_derivative(&self.xs, &self.mx, i, s);
        let ty = Self::spline_derivative(&self.ys, &self.my, i, s);
        [tx, ty]
    }
    /// Derivative of one cubic segment w.r.t. s.
    fn spline_derivative(vals: &[f64], m: &[f64], i: usize, s: f64) -> f64 {
        let a = 1.0 - s;
        (vals[i + 1] - vals[i])
            + (m[i] * (-3.0 * a * a + 1.0) + m[i + 1] * (3.0 * s * s - 1.0)) / 6.0
    }
}
/// A simplified racing line built from track edge points.
#[allow(dead_code)]
pub struct RacingLine {
    /// Racing line points `[x, y, z]`.
    pub points: Vec<[f64; 3]>,
    /// Indices of apex points (local minima of distance from the inside edge).
    pub apex_indices: Vec<usize>,
}
impl RacingLine {
    /// Build a simple geometric racing line from left and right track edges.
    ///
    /// The initial line is the midpoint between each pair of edge points.
    /// For the apices, points where the midpoint is locally closest to the
    /// inside of a curve are detected by looking for curvature sign changes.
    ///
    /// `n_pts` is the desired number of output points (resampled uniformly).
    #[allow(dead_code)]
    pub fn from_track_edges(left: &[[f64; 3]], right: &[[f64; 3]], n_pts: usize) -> Self {
        let m = left.len().min(right.len());
        assert!(m >= 2, "need at least 2 edge points");
        let n_pts = n_pts.max(2);
        let midpoints: Vec<[f64; 3]> = (0..m)
            .map(|i| {
                [
                    (left[i][0] + right[i][0]) * 0.5,
                    (left[i][1] + right[i][1]) * 0.5,
                    (left[i][2] + right[i][2]) * 0.5,
                ]
            })
            .collect();
        let mut cum = vec![0.0f64; m];
        for i in 1..m {
            let dx = midpoints[i][0] - midpoints[i - 1][0];
            let dy = midpoints[i][1] - midpoints[i - 1][1];
            let dz = midpoints[i][2] - midpoints[i - 1][2];
            cum[i] = cum[i - 1] + (dx * dx + dy * dy + dz * dz).sqrt();
        }
        let total = cum[m - 1];
        let mut points = Vec::with_capacity(n_pts);
        for k in 0..n_pts {
            let s = total * (k as f64) / (n_pts - 1) as f64;
            let idx = cum.partition_point(|&c| c < s).saturating_sub(1).min(m - 2);
            let seg_start = cum[idx];
            let seg_end = cum[idx + 1];
            let seg_len = seg_end - seg_start;
            let t = if seg_len > 1e-15 {
                (s - seg_start) / seg_len
            } else {
                0.0
            };
            let a = midpoints[idx];
            let b = midpoints[idx + 1];
            points.push([
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ]);
        }
        let mut apex_indices = Vec::new();
        for i in 1..n_pts.saturating_sub(1) {
            let prev = points[i - 1];
            let cur = points[i];
            let next = points[i + 1];
            let dx = 0.5 * (next[0] - prev[0]);
            let dy = 0.5 * (next[1] - prev[1]);
            let ddx = next[0] - 2.0 * cur[0] + prev[0];
            let ddy = next[1] - 2.0 * cur[1] + prev[1];
            let cross = dx * ddy - dy * ddx;
            let denom = (dx * dx + dy * dy).sqrt();
            if denom > 1e-10 {
                let kappa = cross.abs() / denom.powi(3);
                if kappa > 0.1 {
                    apex_indices.push(i);
                }
            }
        }
        Self {
            points,
            apex_indices,
        }
    }
}
/// 3-D cubic spline path built from Catmull-Rom tangents.
///
/// Supports uniform arc-length sampling, tangent queries, and curvature
/// queries via `t ∈ [0, 1]`.
pub struct CubicSplinePath {
    /// Control waypoints `[x, y, z]`.
    pub waypoints: Vec<[f64; 3]>,
    /// Catmull-Rom tangent at each waypoint.
    pub tangents: Vec<[f64; 3]>,
    /// Pre-computed total arc length (m).
    pub total_length: f64,
    /// Cumulative arc lengths at each waypoint boundary.
    pub(super) arc_lengths: Vec<f64>,
}
impl CubicSplinePath {
    fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let dz = b[2] - a[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    fn normalize3(v: [f64; 3]) -> [f64; 3] {
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-15);
        [v[0] / len, v[1] / len, v[2] / len]
    }
    /// Evaluate the Hermite cubic at `t` for one segment.
    ///
    /// `t ∈ [0, 1]`.  `p0`, `p1` are endpoints; `m0`, `m1` are tangents
    /// scaled to segment length.
    fn hermite_eval(p0: [f64; 3], p1: [f64; 3], m0: [f64; 3], m1: [f64; 3], t: f64) -> [f64; 3] {
        let t2 = t * t;
        let t3 = t2 * t;
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;
        [
            h00 * p0[0] + h10 * m0[0] + h01 * p1[0] + h11 * m1[0],
            h00 * p0[1] + h10 * m0[1] + h01 * p1[1] + h11 * m1[1],
            h00 * p0[2] + h10 * m0[2] + h01 * p1[2] + h11 * m1[2],
        ]
    }
    /// Derivative of the Hermite cubic w.r.t. `t`.
    fn hermite_deriv(p0: [f64; 3], p1: [f64; 3], m0: [f64; 3], m1: [f64; 3], t: f64) -> [f64; 3] {
        let t2 = t * t;
        let h00 = 6.0 * t2 - 6.0 * t;
        let h10 = 3.0 * t2 - 4.0 * t + 1.0;
        let h01 = -6.0 * t2 + 6.0 * t;
        let h11 = 3.0 * t2 - 2.0 * t;
        [
            h00 * p0[0] + h10 * m0[0] + h01 * p1[0] + h11 * m1[0],
            h00 * p0[1] + h10 * m0[1] + h01 * p1[1] + h11 * m1[1],
            h00 * p0[2] + h10 * m0[2] + h01 * p1[2] + h11 * m1[2],
        ]
    }
    /// Second derivative of the Hermite cubic w.r.t. `t`.
    fn hermite_second_deriv(
        p0: [f64; 3],
        p1: [f64; 3],
        m0: [f64; 3],
        m1: [f64; 3],
        t: f64,
    ) -> [f64; 3] {
        let h00 = 12.0 * t - 6.0;
        let h10 = 6.0 * t - 4.0;
        let h01 = -12.0 * t + 6.0;
        let h11 = 6.0 * t - 2.0;
        [
            h00 * p0[0] + h10 * m0[0] + h01 * p1[0] + h11 * m1[0],
            h00 * p0[1] + h10 * m0[1] + h01 * p1[1] + h11 * m1[1],
            h00 * p0[2] + h10 * m0[2] + h01 * p1[2] + h11 * m1[2],
        ]
    }
    /// Map a normalised arc-length parameter `t ∈ [0,1]` to a segment index
    /// and local parameter `s ∈ [0,1]`.
    fn arc_param(&self, t: f64) -> (usize, f64) {
        let n = self.waypoints.len();
        if n < 2 {
            return (0, 0.0);
        }
        let target = (t.clamp(0.0, 1.0) * self.total_length).min(self.total_length);
        let idx = self
            .arc_lengths
            .partition_point(|&a| a < target)
            .saturating_sub(1);
        let i = idx.min(n - 2);
        let seg_start = self.arc_lengths[i];
        let seg_end = self.arc_lengths[i + 1];
        let seg_len = seg_end - seg_start;
        let s = if seg_len > 1e-15 {
            (target - seg_start) / seg_len
        } else {
            0.0
        };
        (i, s.clamp(0.0, 1.0))
    }
    /// Build a `CubicSplinePath` from an ordered list of 3-D waypoints.
    ///
    /// Uses Catmull-Rom tangents: at interior waypoints the tangent is the
    /// centred finite difference of the neighbouring positions.  At the
    /// endpoints the one-sided difference is used.
    #[allow(dead_code)]
    pub fn from_waypoints(pts: &[[f64; 3]]) -> Self {
        let n = pts.len();
        assert!(n >= 2, "CubicSplinePath requires at least 2 waypoints");
        let mut tangents = vec![[0.0f64; 3]; n];
        for i in 0..n {
            let prev = if i == 0 { pts[0] } else { pts[i - 1] };
            let next = if i == n - 1 { pts[n - 1] } else { pts[i + 1] };
            tangents[i] = [
                (next[0] - prev[0]) * 0.5,
                (next[1] - prev[1]) * 0.5,
                (next[2] - prev[2]) * 0.5,
            ];
        }
        let steps = 16usize;
        let mut arc_lengths = vec![0.0f64; n];
        for i in 0..n - 1 {
            let p0 = pts[i];
            let p1 = pts[i + 1];
            let m0 = tangents[i];
            let m1 = tangents[i + 1];
            let mut seg_len = 0.0f64;
            let mut prev_pt = p0;
            for k in 1..=steps {
                let s = k as f64 / steps as f64;
                let cur = Self::hermite_eval(p0, p1, m0, m1, s);
                seg_len += Self::dist3(prev_pt, cur);
                prev_pt = cur;
            }
            arc_lengths[i + 1] = arc_lengths[i] + seg_len;
        }
        let total_length = *arc_lengths.last().unwrap_or(&0.0);
        Self {
            waypoints: pts.to_vec(),
            tangents,
            total_length,
            arc_lengths,
        }
    }
    /// Sample the path at normalised arc-length parameter `t ∈ [0, 1]`.
    #[allow(dead_code)]
    pub fn sample(&self, t: f64) -> [f64; 3] {
        let n = self.waypoints.len();
        if n == 1 {
            return self.waypoints[0];
        }
        let (i, s) = self.arc_param(t);
        Self::hermite_eval(
            self.waypoints[i],
            self.waypoints[i + 1],
            self.tangents[i],
            self.tangents[i + 1],
            s,
        )
    }
    /// Normalised tangent direction at `t ∈ [0, 1]`.
    #[allow(dead_code)]
    pub fn tangent_at(&self, t: f64) -> [f64; 3] {
        let n = self.waypoints.len();
        if n == 1 {
            return [1.0, 0.0, 0.0];
        }
        let (i, s) = self.arc_param(t);
        let d = Self::hermite_deriv(
            self.waypoints[i],
            self.waypoints[i + 1],
            self.tangents[i],
            self.tangents[i + 1],
            s,
        );
        Self::normalize3(d)
    }
    /// Unsigned curvature (1/m) at `t ∈ [0, 1]`.
    ///
    /// Uses the standard formula κ = |r' × r''| / |r'|³.
    #[allow(dead_code)]
    pub fn curvature_at(&self, t: f64) -> f64 {
        let n = self.waypoints.len();
        if n == 1 {
            return 0.0;
        }
        let (i, s) = self.arc_param(t);
        let d1 = Self::hermite_deriv(
            self.waypoints[i],
            self.waypoints[i + 1],
            self.tangents[i],
            self.tangents[i + 1],
            s,
        );
        let d2 = Self::hermite_second_deriv(
            self.waypoints[i],
            self.waypoints[i + 1],
            self.tangents[i],
            self.tangents[i + 1],
            s,
        );
        let cx = d1[1] * d2[2] - d1[2] * d2[1];
        let cy = d1[2] * d2[0] - d1[0] * d2[2];
        let cz = d1[0] * d2[1] - d1[1] * d2[0];
        let cross_mag = (cx * cx + cy * cy + cz * cz).sqrt();
        let d1_mag = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        let denom = d1_mag.powi(3);
        if denom > 1e-15 {
            cross_mag / denom
        } else {
            0.0
        }
    }
}
/// A single target point along a planned path.
#[derive(Debug, Clone)]
pub struct Waypoint {
    /// 3-D position \[x, y, z\] in world space.
    pub position: [f64; 3],
    /// Advisory speed limit at this waypoint (m/s).
    pub speed_limit: f64,
    /// Desired heading at this waypoint (radians, measured from +X axis).
    pub heading: f64,
}
/// A complete track description made from an ordered list of segments.
#[allow(dead_code)]
pub struct TrackLayout {
    /// Track name.
    pub name: String,
    /// Ordered segments that make up the track.
    pub segments: Vec<TrackSegment>,
}
impl TrackLayout {
    /// Create an empty track with the given name.
    #[allow(dead_code)]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            segments: Vec::new(),
        }
    }
    /// Append a segment to the track.
    #[allow(dead_code)]
    pub fn add_segment(&mut self, seg: TrackSegment) {
        self.segments.push(seg);
    }
    /// Total track length in metres.
    #[allow(dead_code)]
    pub fn total_length(&self) -> f64 {
        self.segments.iter().map(|s| s.length()).sum()
    }
    /// Number of segments.
    #[allow(dead_code)]
    pub fn num_segments(&self) -> usize {
        self.segments.len()
    }
    /// Number of corner segments (arcs and chicanes).
    #[allow(dead_code)]
    pub fn num_corners(&self) -> usize {
        self.segments
            .iter()
            .filter(|s| !matches!(s, TrackSegment::Straight { .. }))
            .count()
    }
    /// Total corner length (sum of arc and chicane lengths) in metres.
    #[allow(dead_code)]
    pub fn corner_length(&self) -> f64 {
        self.segments
            .iter()
            .filter(|s| !matches!(s, TrackSegment::Straight { .. }))
            .map(|s| s.length())
            .sum()
    }
    /// Fraction of the track made up of corners (0–1).
    #[allow(dead_code)]
    pub fn corner_fraction(&self) -> f64 {
        let total = self.total_length();
        if total < 1e-12 {
            return 0.0;
        }
        self.corner_length() / total
    }
    /// Sample the track at arc-length distance `s` (measured from the start
    /// of the first segment).
    ///
    /// Returns the 2-D position `[x, y]`, or `None` if the track has no segments.
    #[allow(dead_code)]
    pub fn sample_at_distance(&self, s: f64) -> Option<[f64; 2]> {
        if self.segments.is_empty() {
            return None;
        }
        let mut remaining = s.max(0.0);
        for seg in &self.segments {
            let len = seg.length();
            if remaining <= len || len < 1e-15 {
                let t = if len > 1e-15 { remaining / len } else { 0.0 };
                return Some(seg.sample(t.min(1.0)));
            }
            remaining -= len;
        }
        Some(self.segments.last()?.sample(1.0))
    }
}
/// Result of the minimum-curvature racing line optimisation.
///
/// Each entry in `points` and `alphas` corresponds to one input row
/// (left/right edge pair).  `alpha = 0` is the left edge; `alpha = 1` is
/// the right edge.
#[allow(dead_code)]
pub struct MinCurvatureRacingLine {
    /// Optimised 2-D racing line points.
    pub points: Vec<[f64; 2]>,
    /// Lateral position parameter α ∈ \[0, 1\] for each row.
    pub alphas: Vec<f64>,
    /// Curvature at each racing line point (1/m).
    pub curvatures: Vec<f64>,
}
impl MinCurvatureRacingLine {
    /// Compute the minimum-curvature racing line by iterative gradient descent
    /// on the curvature energy.
    ///
    /// # Arguments
    /// * `left`  – ordered left-edge points `[x, y]`
    /// * `right` – ordered right-edge points `[x, y]` (same count as `left`)
    /// * `iters` – number of gradient-descent iterations (e.g. 50)
    ///
    /// The optimisation minimises `Σ κᵢ²` subject to each point lying on the
    /// segment between the corresponding left and right edge points.
    #[allow(dead_code)]
    pub fn optimize(left: &[[f64; 2]], right: &[[f64; 2]], iters: usize) -> Self {
        let n = left.len().min(right.len());
        assert!(n >= 3, "need at least 3 edge pairs");
        let mut alphas = vec![0.5f64; n];
        let step_size = 0.02_f64;
        for _ in 0..iters {
            let pts = Self::alphas_to_points(left, right, &alphas);
            let kappas = compute_path_curvature(&pts);
            let mut grad = vec![0.0f64; n];
            let delta = 1e-4;
            for i in 1..n - 1 {
                let mut ap = alphas.clone();
                ap[i] = (alphas[i] + delta).clamp(0.0, 1.0);
                let pp = Self::alphas_to_points(left, right, &ap);
                let kp = compute_path_curvature(&pp);
                let mut am = alphas.clone();
                am[i] = (alphas[i] - delta).clamp(0.0, 1.0);
                let pm = Self::alphas_to_points(left, right, &am);
                let km = compute_path_curvature(&pm);
                let e_plus: f64 = kp.iter().map(|k| k * k).sum();
                let e_minus: f64 = km.iter().map(|k| k * k).sum();
                grad[i] = (e_plus - e_minus) / (2.0 * delta);
                let _ = &kappas;
            }
            for i in 1..n - 1 {
                alphas[i] = (alphas[i] - step_size * grad[i]).clamp(0.0, 1.0);
            }
        }
        let points = Self::alphas_to_points(left, right, &alphas);
        let curvatures = compute_path_curvature(&points);
        MinCurvatureRacingLine {
            points,
            alphas,
            curvatures,
        }
    }
    /// Convert α values to 2-D points.
    fn alphas_to_points(left: &[[f64; 2]], right: &[[f64; 2]], alphas: &[f64]) -> Vec<[f64; 2]> {
        let n = left.len().min(right.len()).min(alphas.len());
        (0..n)
            .map(|i| {
                let a = alphas[i].clamp(0.0, 1.0);
                [
                    left[i][0] + a * (right[i][0] - left[i][0]),
                    left[i][1] + a * (right[i][1] - left[i][1]),
                ]
            })
            .collect()
    }
    /// Mean curvature of the optimised line.
    #[allow(dead_code)]
    pub fn mean_curvature(&self) -> f64 {
        if self.curvatures.is_empty() {
            return 0.0;
        }
        self.curvatures.iter().sum::<f64>() / self.curvatures.len() as f64
    }
    /// Maximum curvature on the optimised line.
    #[allow(dead_code)]
    pub fn max_curvature(&self) -> f64 {
        self.curvatures.iter().cloned().fold(0.0_f64, f64::max)
    }
}
/// A feedforward+feedback controller that tracks a `CubicSplinePath`
/// using an arc-length parametrisation.
///
/// The controller computes:
/// 1. The current arc-length fraction `s/L`.
/// 2. A speed set-point limited by curvature and grip.
/// 3. A steering command using pure-pursuit with adaptive look-ahead.
#[allow(dead_code)]
pub struct ArcLengthController {
    /// The reference path.
    pub path: CubicSplinePath,
    /// Wheelbase (m).
    pub wheelbase: f64,
    /// Minimum look-ahead distance (m).
    pub lookahead_min: f64,
    /// Maximum speed allowed (m/s).
    pub v_max: f64,
    /// Internal pure-pursuit controller.
    pub(super) pursuit: PurePursuitController,
}
impl ArcLengthController {
    /// Create a new controller.
    ///
    /// * `path`         – reference path
    /// * `wheelbase`    – vehicle wheelbase (m)
    /// * `lookahead_min` – minimum look-ahead distance (m)
    /// * `v_max`        – maximum speed set-point (m/s)
    #[allow(dead_code)]
    pub fn new(path: CubicSplinePath, wheelbase: f64, lookahead_min: f64, v_max: f64) -> Self {
        let pursuit = PurePursuitController::new(lookahead_min, wheelbase);
        Self {
            path,
            wheelbase,
            lookahead_min,
            v_max,
            pursuit,
        }
    }
    /// Convert an arc-length distance `s` to a normalised fraction `t ∈ [0, 1]`.
    #[allow(dead_code)]
    pub fn arc_fraction(&self, s: f64) -> f64 {
        if self.path.total_length < 1e-15 {
            return 0.0;
        }
        (s / self.path.total_length).clamp(0.0, 1.0)
    }
    /// Compute the speed set-point at arc distance `s` (m/s).
    ///
    /// Limited by curvature and grip: `v = sqrt(mu * g / kappa)`.
    #[allow(dead_code)]
    pub fn target_speed_at_arc(&self, s: f64, friction: f64, gravity: f64) -> f64 {
        let t = self.arc_fraction(s);
        let kappa = self.path.curvature_at(t).max(1e-9);
        let v_lat = (friction * gravity / kappa).sqrt();
        v_lat.min(self.v_max).max(1.0)
    }
    /// Compute the steering command (rad) to follow the path from the given
    /// vehicle state.
    ///
    /// * `vehicle_pos`     – 3-D world position `[x, y, z]`
    /// * `vehicle_heading` – yaw angle (rad from +X)
    /// * `current_arc`     – current arc-length position along the path (m)
    ///
    /// The look-ahead distance is adapted based on curvature:
    /// shorter in tight corners, longer on straights.
    #[allow(dead_code)]
    pub fn steering_command(
        &mut self,
        vehicle_pos: [f64; 3],
        vehicle_heading: f64,
        current_arc: f64,
    ) -> f64 {
        let t = self.arc_fraction(current_arc);
        let kappa = self.path.curvature_at(t).max(1e-12);
        let k_curv = 5.0;
        let lookahead = curvature_dependent_lookahead(kappa, 20.0, k_curv, self.lookahead_min);
        self.pursuit.lookahead_distance = lookahead;
        let lookahead_arc = (current_arc + lookahead).min(self.path.total_length);
        let t_ahead = self.arc_fraction(lookahead_arc);
        let target_3d = {
            let p = self.path.sample(t_ahead);
            [p[0], p[1], p[2]]
        };
        self.pursuit
            .steering_angle(vehicle_pos, vehicle_heading, target_3d)
    }
    /// Nearest arc-length on the path to the given world position.
    ///
    /// Searches over `n_samples` uniformly spaced samples and returns the arc
    /// distance `s ∈ [0, total_length]`.
    #[allow(dead_code)]
    pub fn nearest_arc(&self, world_pos: [f64; 3], n_samples: usize) -> f64 {
        let n = n_samples.max(2);
        let mut best_s = 0.0f64;
        let mut best_d2 = f64::INFINITY;
        for k in 0..n {
            let t = k as f64 / (n - 1) as f64;
            let p = self.path.sample(t);
            let dx = p[0] - world_pos[0];
            let dy = p[1] - world_pos[1];
            let dz = p[2] - world_pos[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_d2 {
                best_d2 = d2;
                best_s = t * self.path.total_length;
            }
        }
        best_s
    }
}
/// Frenet–Serret frame at a point on a `CubicSplinePath`.
///
/// The three orthonormal vectors are:
/// - `tangent`  – unit tangent T (direction of travel)
/// - `normal`   – unit principal normal N (pointing toward the centre of curvature)
/// - `binormal` – unit binormal B = T × N
#[allow(dead_code)]
pub struct FrenetFrame {
    /// Unit tangent vector.
    pub tangent: [f64; 3],
    /// Unit principal normal vector (zero if curvature ≈ 0).
    pub normal: [f64; 3],
    /// Unit binormal vector B = T × N.
    pub binormal: [f64; 3],
    /// Curvature κ (1/m) at this point.
    pub curvature: f64,
}
impl FrenetFrame {
    /// Compute the Frenet frame at the given arc-length parameter `t ∈ [0, 1]`.
    ///
    /// Uses the first and second derivatives of the spline at `t` to construct
    /// the Frenet–Serret basis.
    #[allow(dead_code)]
    pub fn compute(path: &CubicSplinePath, t: f64) -> Self {
        let n = path.waypoints.len();
        if n < 2 {
            return FrenetFrame {
                tangent: [1.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                curvature: 0.0,
            };
        }
        let tangent = vec3_norm(path.tangent_at(t));
        let eps = 1e-5_f64;
        let t1 = (t + eps).clamp(0.0, 1.0);
        let t0 = (t - eps).clamp(0.0, 1.0);
        let tan1 = path.tangent_at(t1);
        let tan0 = path.tangent_at(t0);
        let dt_vec = [
            (tan1[0] - tan0[0]) / (t1 - t0 + 1e-30),
            (tan1[1] - tan0[1]) / (t1 - t0 + 1e-30),
            (tan1[2] - tan0[2]) / (t1 - t0 + 1e-30),
        ];
        let curvature_mag = vec3_len(dt_vec);
        let normal = if curvature_mag > 1e-15 {
            vec3_norm(dt_vec)
        } else {
            let perp = if tangent[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            vec3_norm(vec3_cross(vec3_cross(tangent, perp), tangent))
        };
        let binormal = vec3_norm(vec3_cross(tangent, normal));
        FrenetFrame {
            tangent,
            normal,
            binormal,
            curvature: curvature_mag,
        }
    }
    /// Returns `true` if the Frenet basis is right-handed (det ≈ +1).
    #[allow(dead_code)]
    pub fn is_right_handed(&self) -> bool {
        let det_approx = vec3_cross(self.tangent, self.normal)[0] * self.binormal[0]
            + vec3_cross(self.tangent, self.normal)[1] * self.binormal[1]
            + vec3_cross(self.tangent, self.normal)[2] * self.binormal[2];
        det_approx > 0.0
    }
}
/// A segment of a racing circuit.
///
/// Tracks are composed of straights, arcs (constant-radius corners), and
/// chicanes (S-bends or other combined-corner sections).
#[derive(Debug, Clone)]
pub enum TrackSegment {
    /// A straight section between two 2-D points.
    Straight {
        /// Start point `[x, y]`.
        start: [f64; 2],
        /// End point `[x, y]`.
        end: [f64; 2],
    },
    /// A constant-radius arc.
    Arc {
        /// Centre of the arc `[cx, cy]`.
        center: [f64; 2],
        /// Radius in metres.
        radius: f64,
        /// Start angle in radians (measured from +X axis).
        start_angle: f64,
        /// Signed sweep angle in radians (positive = counter-clockwise).
        sweep: f64,
    },
    /// A chicane defined by a polyline of 2-D waypoints.
    Chicane {
        /// Ordered waypoints forming the chicane.
        waypoints: Vec<[f64; 2]>,
    },
}
impl TrackSegment {
    /// Arc-length of this segment in metres.
    #[allow(dead_code)]
    pub fn length(&self) -> f64 {
        match self {
            TrackSegment::Straight { start, end } => {
                let dx = end[0] - start[0];
                let dy = end[1] - start[1];
                (dx * dx + dy * dy).sqrt()
            }
            TrackSegment::Arc { radius, sweep, .. } => radius * sweep.abs(),
            TrackSegment::Chicane { waypoints } => {
                let n = waypoints.len();
                if n < 2 {
                    return 0.0;
                }
                let mut len = 0.0;
                for i in 0..n - 1 {
                    let dx = waypoints[i + 1][0] - waypoints[i][0];
                    let dy = waypoints[i + 1][1] - waypoints[i][1];
                    len += (dx * dx + dy * dy).sqrt();
                }
                len
            }
        }
    }
    /// Sample a point on the segment at fraction `t ∈ [0, 1]`.
    #[allow(dead_code)]
    pub fn sample(&self, t: f64) -> [f64; 2] {
        let t = t.clamp(0.0, 1.0);
        match self {
            TrackSegment::Straight { start, end } => [
                start[0] + t * (end[0] - start[0]),
                start[1] + t * (end[1] - start[1]),
            ],
            TrackSegment::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => {
                let angle = start_angle + t * sweep;
                [
                    center[0] + radius * angle.cos(),
                    center[1] + radius * angle.sin(),
                ]
            }
            TrackSegment::Chicane { waypoints } => {
                let n = waypoints.len();
                if n == 0 {
                    return [0.0; 2];
                }
                if n == 1 {
                    return waypoints[0];
                }
                let total = self.length();
                let target = t * total;
                let mut acc = 0.0;
                for i in 0..n - 1 {
                    let dx = waypoints[i + 1][0] - waypoints[i][0];
                    let dy = waypoints[i + 1][1] - waypoints[i][1];
                    let seg = (dx * dx + dy * dy).sqrt();
                    if acc + seg >= target {
                        let s = if seg > 1e-15 {
                            (target - acc) / seg
                        } else {
                            0.0
                        };
                        return [waypoints[i][0] + s * dx, waypoints[i][1] + s * dy];
                    }
                    acc += seg;
                }
                *waypoints.last().expect("collection should not be empty")
            }
        }
    }
    /// Mean curvature over this segment (1/m).
    ///
    /// Straights have zero curvature; arcs have constant curvature = 1/R;
    /// chicanes report the average of three-point curvatures.
    #[allow(dead_code)]
    pub fn mean_curvature(&self) -> f64 {
        match self {
            TrackSegment::Straight { .. } => 0.0,
            TrackSegment::Arc { radius, .. } => {
                if radius.abs() < 1e-12 {
                    0.0
                } else {
                    1.0 / radius.abs()
                }
            }
            TrackSegment::Chicane { waypoints } => {
                let n = waypoints.len();
                if n < 3 {
                    return 0.0;
                }
                let k: Vec<f64> = compute_path_curvature(waypoints);
                k.iter().sum::<f64>() / k.len() as f64
            }
        }
    }
}
