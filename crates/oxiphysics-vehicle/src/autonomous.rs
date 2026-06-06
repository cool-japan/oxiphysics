// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Autonomous driving algorithms.
//!
//! Provides path representation, lateral control laws, obstacle avoidance
//! via potential fields, and longitudinal speed planning.
//!
//! # Overview
//!
//! - [`Waypoint`] — a path waypoint with position, heading, and speed limit
//! - [`Path`] — ordered list of waypoints with geometry helpers
//! - [`PurePursuitController`] — classic look-ahead lateral controller
//! - [`StanleyController`] — front-axle lateral controller (Stanley method)
//! - [`LaneKeeping`] — simple lane-keeping correction from lateral/angular error
//! - [`ObstacleAvoidance`] — potential-field repulsive force computation
//! - [`SpeedPlanner`] — curvature-aware longitudinal speed planning

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 2-D vector helpers (plain [f64;2] arrays)
// ---------------------------------------------------------------------------

#[inline]
fn v2_add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn v2_sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn v2_scale(a: [f64; 2], s: f64) -> [f64; 2] {
    [a[0] * s, a[1] * s]
}

#[inline]
fn v2_dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[inline]
fn v2_norm(a: [f64; 2]) -> f64 {
    v2_dot(a, a).sqrt()
}

#[inline]
fn v2_dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    v2_norm(v2_sub(a, b))
}

#[inline]
fn v2_normalize(a: [f64; 2]) -> [f64; 2] {
    let n = v2_norm(a);
    if n > 1e-12 {
        v2_scale(a, 1.0 / n)
    } else {
        [0.0; 2]
    }
}

/// Wrap an angle to \[-π, π\].
#[inline]
fn wrap_angle(a: f64) -> f64 {
    let mut r = a;
    while r > PI {
        r -= 2.0 * PI;
    }
    while r < -PI {
        r += 2.0 * PI;
    }
    r
}

// ---------------------------------------------------------------------------
// Waypoint
// ---------------------------------------------------------------------------

/// A single waypoint in a planned path.
#[derive(Debug, Clone)]
pub struct Waypoint {
    /// 2-D position `[x, y]` in metres.
    pub position: [f64; 2],
    /// Maximum speed at this waypoint in m/s.
    pub speed_limit: f64,
    /// Desired heading at this waypoint in radians.
    pub heading: f64,
}

impl Waypoint {
    /// Create a new `Waypoint`.
    pub fn new(position: [f64; 2], speed_limit: f64, heading: f64) -> Self {
        Self {
            position,
            speed_limit,
            heading,
        }
    }
}

// ---------------------------------------------------------------------------
// Path
// ---------------------------------------------------------------------------

/// An ordered sequence of [`Waypoint`]s forming a planned path.
#[derive(Debug, Clone, Default)]
pub struct Path {
    /// Ordered list of waypoints.
    pub waypoints: Vec<Waypoint>,
}

impl Path {
    /// Create an empty path.
    pub fn new() -> Self {
        Self {
            waypoints: Vec::new(),
        }
    }

    /// Add a waypoint.
    pub fn push(&mut self, wp: Waypoint) {
        self.waypoints.push(wp);
    }

    /// Total arc length of the path (sum of Euclidean segment lengths).
    pub fn total_length(&self) -> f64 {
        if self.waypoints.len() < 2 {
            return 0.0;
        }
        self.waypoints
            .windows(2)
            .map(|w| v2_dist(w[0].position, w[1].position))
            .sum()
    }

    /// Estimate signed curvature at each interior waypoint using a
    /// three-point formula.  Returns a `Vec` with length `waypoints.len()`;
    /// the first and last entries are copies of their nearest interior value
    /// (or zero if fewer than 3 waypoints exist).
    pub fn compute_curvature(&self) -> Vec<f64> {
        let n = self.waypoints.len();
        if n < 3 {
            return vec![0.0; n];
        }
        let mut curvatures = vec![0.0; n];
        for (idx, w) in self.waypoints.windows(3).enumerate() {
            let p0 = w[0].position;
            let p1 = w[1].position;
            let p2 = w[2].position;
            curvatures[idx + 1] = three_point_curvature(p0, p1, p2);
        }
        curvatures[0] = curvatures[1];
        curvatures[n - 1] = curvatures[n - 2];
        curvatures
    }

    /// Find the index of the waypoint closest to `pos`.
    pub fn nearest_index(&self, pos: [f64; 2]) -> usize {
        self.waypoints
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                v2_dist(a.position, pos)
                    .partial_cmp(&v2_dist(b.position, pos))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// Find the first waypoint ahead of `pos` at or beyond `lookahead` distance.
    ///
    /// Returns `None` if no such waypoint exists.
    pub fn lookahead_point(&self, pos: [f64; 2], lookahead: f64) -> Option<[f64; 2]> {
        for wp in &self.waypoints {
            if v2_dist(pos, wp.position) >= lookahead {
                return Some(wp.position);
            }
        }
        self.waypoints.last().map(|w| w.position)
    }

    /// Compute the signed lateral error of `pos` from the nearest path segment.
    ///
    /// Positive when `pos` is to the left of the path direction.
    pub fn lateral_error(&self, pos: [f64; 2]) -> f64 {
        if self.waypoints.len() < 2 {
            return 0.0;
        }
        let idx = self.nearest_index(pos);
        let (p0, p1) = if idx + 1 < self.waypoints.len() {
            (
                self.waypoints[idx].position,
                self.waypoints[idx + 1].position,
            )
        } else {
            (
                self.waypoints[idx - 1].position,
                self.waypoints[idx].position,
            )
        };
        let seg = v2_sub(p1, p0);
        let to_pos = v2_sub(pos, p0);
        // cross product (z-component) gives signed area / |seg|
        let cross = seg[0] * to_pos[1] - seg[1] * to_pos[0];
        let seg_len = v2_norm(seg);
        if seg_len < 1e-12 {
            0.0
        } else {
            cross / seg_len
        }
    }
}

/// Menger curvature from three consecutive points.
///
/// Formula: κ = 4·Area / (a·b·c) where Area is the triangle area.
fn three_point_curvature(p0: [f64; 2], p1: [f64; 2], p2: [f64; 2]) -> f64 {
    let a = v2_dist(p0, p1);
    let b = v2_dist(p1, p2);
    let c = v2_dist(p2, p0);
    // cross product magnitude = 2 * triangle area
    let cross = {
        let d01 = v2_sub(p1, p0);
        let d02 = v2_sub(p2, p0);
        (d01[0] * d02[1] - d01[1] * d02[0]).abs()
    };
    // area = cross / 2 → κ = 4 * area / (a*b*c) = 2 * cross / (a*b*c)
    let denom = a * b * c;
    if denom < 1e-20 {
        0.0
    } else {
        2.0 * cross / denom
    }
}

// ---------------------------------------------------------------------------
// PurePursuitController
// ---------------------------------------------------------------------------

/// Pure pursuit lateral controller.
///
/// Computes a steering angle to steer a vehicle with a given wheelbase toward
/// a look-ahead point on the path.
#[derive(Debug, Clone)]
pub struct PurePursuitController {
    /// Look-ahead distance in metres.
    pub lookahead_distance: f64,
    /// Axle-to-axle wheelbase in metres.
    pub wheelbase: f64,
}

impl PurePursuitController {
    /// Create a new `PurePursuitController`.
    pub fn new(lookahead_distance: f64, wheelbase: f64) -> Self {
        Self {
            lookahead_distance,
            wheelbase,
        }
    }

    /// Compute the required front-wheel steering angle (radians) to track
    /// `path` from the current vehicle position `pos` and heading `heading`.
    ///
    /// Returns a steering angle in `[-π/2, π/2]`.
    pub fn compute_steering(&self, pos: [f64; 2], heading: f64, path: &Path) -> f64 {
        let target = match path.lookahead_point(pos, self.lookahead_distance) {
            Some(t) => t,
            None => return 0.0,
        };
        // Transform target into vehicle frame
        let dx = target[0] - pos[0];
        let dy = target[1] - pos[1];
        let cos_h = heading.cos();
        let sin_h = heading.sin();
        // Lateral offset in vehicle frame
        let local_y = -dx * sin_h + dy * cos_h;
        let dist = v2_dist(pos, target).max(1e-6);
        // Pure pursuit: delta = atan(2 * L * sin(alpha) / ld)
        // where sin(alpha) ≈ local_y / dist
        let alpha = (local_y / dist).clamp(-1.0, 1.0).asin();
        (2.0 * self.wheelbase * alpha.sin() / dist).atan()
    }
}

// ---------------------------------------------------------------------------
// StanleyController
// ---------------------------------------------------------------------------

/// Stanley lateral controller (front-axle method).
///
/// Computes steering as the sum of heading error and cross-track error term:
/// `δ = ψ_e + atan(k * e / v)`
#[derive(Debug, Clone)]
pub struct StanleyController {
    /// Cross-track error gain.
    pub k_gain: f64,
    /// Softening parameter to avoid division by zero at low speeds (m/s).
    pub speed_softening: f64,
    /// Maximum steering angle magnitude (rad).
    pub max_steer: f64,
}

impl StanleyController {
    /// Create a new `StanleyController` with given gain.
    pub fn new(k_gain: f64) -> Self {
        Self {
            k_gain,
            speed_softening: 1.0,
            max_steer: PI / 4.0,
        }
    }

    /// Set a custom softening constant and return `self`.
    pub fn with_softening(mut self, s: f64) -> Self {
        self.speed_softening = s;
        self
    }

    /// Compute front-wheel steering angle (radians).
    ///
    /// * `pos` — current vehicle position `[x, y]`
    /// * `heading` — current heading (yaw) in radians
    /// * `speed` — current speed in m/s (used for normalization)
    /// * `path` — reference path
    pub fn compute_steering(&self, pos: [f64; 2], heading: f64, speed: f64, path: &Path) -> f64 {
        if path.waypoints.len() < 2 {
            return 0.0;
        }
        // Heading error: difference between path heading at nearest point and vehicle heading
        let idx = path.nearest_index(pos);
        let path_heading = if idx + 1 < path.waypoints.len() {
            let p0 = path.waypoints[idx].position;
            let p1 = path.waypoints[idx + 1].position;
            let d = v2_sub(p1, p0);
            d[1].atan2(d[0])
        } else {
            path.waypoints[idx].heading
        };
        let heading_error = wrap_angle(path_heading - heading);
        // Cross-track error (signed lateral deviation)
        let cte = path.lateral_error(pos);
        // Stanley formula: negate cte because our lateral_error is positive-left,
        // but the standard convention for front-axle error is positive-right.
        let steer =
            heading_error + (-self.k_gain * cte / (speed.abs() + self.speed_softening)).atan();
        steer.clamp(-self.max_steer, self.max_steer)
    }
}

// ---------------------------------------------------------------------------
// LaneKeeping
// ---------------------------------------------------------------------------

/// Simple lane-keeping controller using proportional lateral and angular error.
///
/// Outputs a steering correction proportional to:
/// `correction = k_lateral * lateral_error + k_angular * angular_error`
#[derive(Debug, Clone)]
pub struct LaneKeeping {
    /// Nominal lane width in metres (informational).
    pub lane_width: f64,
    /// Signed lateral distance from lane centre (m); positive = left of centre.
    pub lateral_error: f64,
    /// Signed angular error relative to lane direction (rad).
    pub angular_error: f64,
    /// Proportional gain on lateral error.
    pub k_lateral: f64,
    /// Proportional gain on angular error.
    pub k_angular: f64,
}

impl LaneKeeping {
    /// Create a new `LaneKeeping` controller.
    pub fn new(lane_width: f64, k_lateral: f64, k_angular: f64) -> Self {
        Self {
            lane_width,
            lateral_error: 0.0,
            angular_error: 0.0,
            k_lateral,
            k_angular,
        }
    }

    /// Update current errors.
    pub fn update(&mut self, lateral_error: f64, angular_error: f64) {
        self.lateral_error = lateral_error;
        self.angular_error = angular_error;
    }

    /// Compute steering correction (rad).
    pub fn compute_correction(&self) -> f64 {
        -(self.k_lateral * self.lateral_error + self.k_angular * self.angular_error)
    }
}

// ---------------------------------------------------------------------------
// ObstacleAvoidance
// ---------------------------------------------------------------------------

/// Potential-field obstacle avoidance.
///
/// Obstacles generate a repulsive force that decays with distance, steering
/// the vehicle away when it enters the influence radius.
#[derive(Debug, Clone)]
pub struct ObstacleAvoidance {
    /// List of obstacle positions `[x, y]`.
    pub obstacles: Vec<[f64; 2]>,
    /// Repulsive gain.
    pub repulsive_gain: f64,
    /// Influence radius beyond which obstacles are ignored (m).
    pub influence_radius: f64,
}

impl ObstacleAvoidance {
    /// Create an empty `ObstacleAvoidance`.
    pub fn new(repulsive_gain: f64, influence_radius: f64) -> Self {
        Self {
            obstacles: Vec::new(),
            repulsive_gain,
            influence_radius,
        }
    }

    /// Add an obstacle at `pos`.
    pub fn add_obstacle(&mut self, pos: [f64; 2]) {
        self.obstacles.push(pos);
    }

    /// Compute the total repulsive force vector `[fx, fy]` at `pos`.
    ///
    /// Each obstacle within `influence_radius` contributes a force proportional
    /// to `(1/d - 1/d_inf)² / d²` pointing away from the obstacle.
    pub fn repulsive_force(&self, pos: [f64; 2]) -> [f64; 2] {
        let mut total = [0.0_f64; 2];
        for &obs in &self.obstacles {
            let diff = v2_sub(pos, obs);
            let dist = v2_norm(diff).max(1e-6);
            if dist >= self.influence_radius {
                continue;
            }
            // Potential: 0.5 * k * (1/d - 1/d_inf)^2
            // Gradient (force = -dU): k * (1/d - 1/d_inf) * (1/d^2) * (pos-obs)/d
            let factor = self.repulsive_gain
                * (1.0 / dist - 1.0 / self.influence_radius)
                * (1.0 / (dist * dist));
            let dir = v2_normalize(diff);
            total = v2_add(total, v2_scale(dir, factor));
        }
        total
    }

    /// Compute the repulsive steering direction (heading change) at `pos`
    /// given the current vehicle heading.  Projects force onto lateral axis.
    pub fn steering_correction(&self, pos: [f64; 2], heading: f64) -> f64 {
        let force = self.repulsive_force(pos);
        // Lateral component in vehicle frame

        -force[0] * heading.sin() + force[1] * heading.cos()
    }
}

// ---------------------------------------------------------------------------
// SpeedPlanner
// ---------------------------------------------------------------------------

/// Longitudinal speed planner based on path curvature.
///
/// Slows the vehicle in curves to maintain lateral acceleration below a
/// comfort limit.
#[derive(Debug, Clone)]
pub struct SpeedPlanner {
    /// Maximum longitudinal deceleration (m/s²; positive value).
    pub comfort_decel: f64,
    /// Maximum allowed speed (m/s).
    pub max_speed: f64,
    /// Maximum allowed lateral acceleration (m/s²).
    pub max_lat_accel: f64,
}

impl SpeedPlanner {
    /// Create a new `SpeedPlanner`.
    pub fn new(comfort_decel: f64, max_speed: f64, max_lat_accel: f64) -> Self {
        Self {
            comfort_decel,
            max_speed,
            max_lat_accel,
        }
    }

    /// Plan a comfortable speed for the given path curvature (1/m).
    ///
    /// Returns the maximum speed consistent with `max_lat_accel`:
    /// `v = sqrt(a_lat_max / |κ|)` clamped to `[0, max_speed]`.
    pub fn plan_speed(&self, curvature: f64) -> f64 {
        let kappa = curvature.abs();
        if kappa < 1e-9 {
            return self.max_speed;
        }
        let v_curve = (self.max_lat_accel / kappa).sqrt();
        v_curve.min(self.max_speed).max(0.0)
    }

    /// Plan speeds for every waypoint of a path.
    ///
    /// Returns a `Vec`f64` with one speed per waypoint.
    pub fn plan_path_speeds(&self, path: &Path) -> Vec<f64> {
        let curvatures = path.compute_curvature();
        curvatures.iter().map(|&k| self.plan_speed(k)).collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// Build a straight path along X from 0 to `len` with `n` waypoints.
    fn straight_path(len: f64, n: usize) -> Path {
        let mut p = Path::new();
        for i in 0..n {
            let x = len * i as f64 / (n - 1).max(1) as f64;
            p.push(Waypoint::new([x, 0.0], 20.0, 0.0));
        }
        p
    }

    /// Build a circular arc in the XY plane with `n` waypoints, radius `r`.
    fn circle_path(r: f64, n: usize) -> Path {
        let mut p = Path::new();
        for i in 0..n {
            let theta = 2.0 * PI * i as f64 / n as f64;
            let x = r * theta.cos();
            let y = r * theta.sin();
            let heading = theta + PI / 2.0;
            p.push(Waypoint::new([x, y], 10.0, heading));
        }
        p
    }

    // ── v2 helpers ──────────────────────────────────────────────────────────

    #[test]
    fn v2_add_basic() {
        assert_eq!(v2_add([1.0, 2.0], [3.0, 4.0]), [4.0, 6.0]);
    }

    #[test]
    fn v2_norm_pythagoras() {
        assert!(approx_eq(v2_norm([3.0, 4.0]), 5.0, 1e-10));
    }

    #[test]
    fn v2_normalize_unit() {
        let n = v2_normalize([3.0, 4.0]);
        assert!(approx_eq(v2_norm(n), 1.0, 1e-10));
    }

    #[test]
    fn v2_normalize_zero_is_zero() {
        assert_eq!(v2_normalize([0.0; 2]), [0.0; 2]);
    }

    #[test]
    fn wrap_angle_plus_pi() {
        assert!(approx_eq(wrap_angle(3.0 * PI), PI, 1e-10));
    }

    #[test]
    fn wrap_angle_minus_pi() {
        assert!(approx_eq(wrap_angle(-3.0 * PI), -PI, 1e-10));
    }

    #[test]
    fn wrap_angle_identity() {
        assert!(approx_eq(wrap_angle(0.5), 0.5, 1e-12));
    }

    // ── Waypoint ─────────────────────────────────────────────────────────────

    #[test]
    fn waypoint_construction() {
        let wp = Waypoint::new([1.0, 2.0], 30.0, PI / 4.0);
        assert_eq!(wp.position, [1.0, 2.0]);
        assert!(approx_eq(wp.speed_limit, 30.0, 1e-10));
    }

    // ── Path ─────────────────────────────────────────────────────────────────

    #[test]
    fn path_total_length_straight() {
        let p = straight_path(10.0, 11);
        assert!(approx_eq(p.total_length(), 10.0, 1e-10));
    }

    #[test]
    fn path_total_length_empty() {
        let p = Path::new();
        assert_eq!(p.total_length(), 0.0);
    }

    #[test]
    fn path_total_length_single_point() {
        let mut p = Path::new();
        p.push(Waypoint::new([0.0, 0.0], 10.0, 0.0));
        assert_eq!(p.total_length(), 0.0);
    }

    #[test]
    fn path_curvature_straight_is_near_zero() {
        let p = straight_path(10.0, 5);
        let curvs = p.compute_curvature();
        for k in curvs {
            assert!(k.abs() < 1e-10, "expected zero curvature, got {k}");
        }
    }

    #[test]
    fn path_curvature_circle_radius_10() {
        let r = 10.0;
        let p = circle_path(r, 100);
        let curvs = p.compute_curvature();
        // Interior curvatures should be close to 1/r
        let interior: Vec<f64> = curvs[1..curvs.len() - 1].to_vec();
        for k in &interior {
            assert!(
                approx_eq(*k, 1.0 / r, 0.02),
                "expected curvature ≈{}, got {k}",
                1.0 / r
            );
        }
    }

    #[test]
    fn path_nearest_index() {
        let p = straight_path(10.0, 11);
        let idx = p.nearest_index([3.0, 0.5]);
        assert_eq!(idx, 3); // nearest to x=3.0
    }

    #[test]
    fn path_lookahead_point() {
        let p = straight_path(10.0, 11);
        let pt = p.lookahead_point([0.0, 0.0], 5.0);
        assert!(pt.is_some());
        let x = pt.unwrap()[0];
        assert!(x >= 5.0, "lookahead x={x} should be >= 5");
    }

    #[test]
    fn path_lateral_error_on_path_is_zero() {
        let p = straight_path(10.0, 5);
        let err = p.lateral_error([5.0, 0.0]);
        assert!(approx_eq(err, 0.0, 1e-10));
    }

    #[test]
    fn path_lateral_error_left_positive() {
        let p = straight_path(10.0, 5);
        // Path along +X; point above path (positive Y in right-hand system = left)
        let err = p.lateral_error([5.0, 1.0]);
        assert!(
            err > 0.0,
            "lateral error should be positive (left), got {err}"
        );
    }

    // ── PurePursuitController ─────────────────────────────────────────────────

    #[test]
    fn pure_pursuit_straight_path_no_steer() {
        let pp = PurePursuitController::new(2.0, 2.5);
        let p = straight_path(20.0, 21);
        let steer = pp.compute_steering([0.0, 0.0], 0.0, &p);
        assert!(
            steer.abs() < 0.01,
            "straight path → near-zero steering, got {steer}"
        );
    }

    #[test]
    fn pure_pursuit_target_left_gives_positive_steer() {
        let pp = PurePursuitController::new(2.0, 2.5);
        // Path that is offset in +Y (to the left when heading = 0)
        let mut p = Path::new();
        p.push(Waypoint::new([5.0, 2.0], 20.0, 0.0)); // only one point far left
        let steer = pp.compute_steering([0.0, 0.0], 0.0, &p);
        assert!(steer > 0.0, "target to left → positive steer, got {steer}");
    }

    #[test]
    fn pure_pursuit_target_right_gives_negative_steer() {
        let pp = PurePursuitController::new(2.0, 2.5);
        let mut p = Path::new();
        p.push(Waypoint::new([5.0, -2.0], 20.0, 0.0));
        let steer = pp.compute_steering([0.0, 0.0], 0.0, &p);
        assert!(steer < 0.0, "target to right → negative steer, got {steer}");
    }

    #[test]
    fn pure_pursuit_empty_path_returns_zero() {
        let pp = PurePursuitController::new(2.0, 2.5);
        let p = Path::new();
        let steer = pp.compute_steering([0.0, 0.0], 0.0, &p);
        assert_eq!(steer, 0.0);
    }

    #[test]
    fn pure_pursuit_larger_wheelbase_more_steer() {
        // In pure pursuit δ = atan(2L·sin α / Ld): larger wheelbase → larger steering angle.
        let mut p = Path::new();
        p.push(Waypoint::new([4.0, 2.0], 20.0, 0.0));
        let pp_small = PurePursuitController::new(2.0, 1.5);
        let pp_large = PurePursuitController::new(2.0, 3.0);
        let s1 = pp_small.compute_steering([0.0, 0.0], 0.0, &p).abs();
        let s2 = pp_large.compute_steering([0.0, 0.0], 0.0, &p).abs();
        assert!(
            s2 > s1,
            "larger wheelbase produces larger steering angle: s1={s1}, s2={s2}"
        );
    }

    // ── StanleyController ─────────────────────────────────────────────────────

    #[test]
    fn stanley_straight_path_near_zero() {
        let sc = StanleyController::new(0.5);
        let p = straight_path(20.0, 21);
        let steer = sc.compute_steering([0.0, 0.0], 0.0, 5.0, &p);
        assert!(
            steer.abs() < 0.01,
            "straight path → near-zero steer, got {steer}"
        );
    }

    #[test]
    fn stanley_lateral_error_produces_steer() {
        let sc = StanleyController::new(1.0);
        let p = straight_path(20.0, 21);
        // Vehicle displaced 1 m to the right of straight path
        let steer = sc.compute_steering([5.0, -1.0], 0.0, 5.0, &p);
        // Should steer left (positive) to correct
        assert!(steer > 0.0, "positive correction expected, got {steer}");
    }

    #[test]
    fn stanley_heading_error_dominates_at_high_speed() {
        let sc = StanleyController::new(0.1).with_softening(0.1);
        let p = straight_path(20.0, 21);
        // Large heading error, small CTE
        let steer_fast = sc.compute_steering([5.0, 0.0], 0.3, 100.0, &p);
        // At high speed CTE term → 0, only heading error remains
        assert!(steer_fast.abs() < PI / 4.0 + 0.01);
    }

    #[test]
    fn stanley_empty_path_returns_zero() {
        let sc = StanleyController::new(1.0);
        let p = Path::new();
        assert_eq!(sc.compute_steering([0.0, 0.0], 0.0, 5.0, &p), 0.0);
    }

    #[test]
    fn stanley_clamps_to_max_steer() {
        let sc = StanleyController::new(100.0); // huge gain
        let p = straight_path(20.0, 21);
        let steer = sc.compute_steering([5.0, -3.0], 0.5, 0.01, &p);
        assert!(steer.abs() <= PI / 4.0 + 1e-10);
    }

    // ── LaneKeeping ──────────────────────────────────────────────────────────

    #[test]
    fn lane_keeping_zero_error_zero_correction() {
        let lk = LaneKeeping::new(3.5, 1.0, 1.0);
        assert!(approx_eq(lk.compute_correction(), 0.0, 1e-12));
    }

    #[test]
    fn lane_keeping_lateral_error_correction() {
        let mut lk = LaneKeeping::new(3.5, 2.0, 0.0);
        lk.update(1.0, 0.0); // 1 m left of centre
        // correction = -(2 * 1 + 0) = -2 (steer right)
        assert!(approx_eq(lk.compute_correction(), -2.0, 1e-10));
    }

    #[test]
    fn lane_keeping_angular_error_correction() {
        let mut lk = LaneKeeping::new(3.5, 0.0, 3.0);
        lk.update(0.0, 0.5);
        assert!(approx_eq(lk.compute_correction(), -1.5, 1e-10));
    }

    #[test]
    fn lane_keeping_combined_errors() {
        let mut lk = LaneKeeping::new(3.5, 2.0, 3.0);
        lk.update(1.0, 0.5);
        // -(2*1 + 3*0.5) = -(2 + 1.5) = -3.5
        assert!(approx_eq(lk.compute_correction(), -3.5, 1e-10));
    }

    // ── ObstacleAvoidance ─────────────────────────────────────────────────────

    #[test]
    fn obstacle_avoidance_no_obstacles_zero_force() {
        let oa = ObstacleAvoidance::new(1.0, 5.0);
        let f = oa.repulsive_force([0.0, 0.0]);
        assert!(approx_eq(f[0], 0.0, 1e-12));
        assert!(approx_eq(f[1], 0.0, 1e-12));
    }

    #[test]
    fn obstacle_outside_influence_radius_no_force() {
        let mut oa = ObstacleAvoidance::new(1.0, 2.0);
        oa.add_obstacle([10.0, 0.0]); // far away
        let f = oa.repulsive_force([0.0, 0.0]);
        assert!(approx_eq(f[0], 0.0, 1e-12));
    }

    #[test]
    fn obstacle_inside_radius_produces_force() {
        let mut oa = ObstacleAvoidance::new(1.0, 5.0);
        oa.add_obstacle([1.0, 0.0]); // 1 m to the right
        let f = oa.repulsive_force([0.0, 0.0]);
        // force should point away from obstacle (negative x)
        assert!(f[0] < 0.0, "should repel in -x direction, got f={f:?}");
    }

    #[test]
    fn obstacle_closer_gives_stronger_force() {
        let mut oa_far = ObstacleAvoidance::new(1.0, 10.0);
        oa_far.add_obstacle([3.0, 0.0]);
        let mut oa_close = ObstacleAvoidance::new(1.0, 10.0);
        oa_close.add_obstacle([1.0, 0.0]);
        let f_far = oa_far.repulsive_force([0.0, 0.0]);
        let f_close = oa_close.repulsive_force([0.0, 0.0]);
        assert!(
            f_close[0].abs() > f_far[0].abs(),
            "closer obstacle should give stronger force"
        );
    }

    #[test]
    fn obstacle_steering_correction_is_finite() {
        let mut oa = ObstacleAvoidance::new(1.0, 5.0);
        oa.add_obstacle([1.0, 0.0]);
        let sc = oa.steering_correction([0.0, 0.0], 0.0);
        assert!(sc.is_finite());
    }

    // ── SpeedPlanner ──────────────────────────────────────────────────────────

    #[test]
    fn speed_planner_zero_curvature_returns_max_speed() {
        let sp = SpeedPlanner::new(3.0, 30.0, 2.0);
        assert!(approx_eq(sp.plan_speed(0.0), 30.0, 1e-10));
    }

    #[test]
    fn speed_planner_high_curvature_reduces_speed() {
        let sp = SpeedPlanner::new(3.0, 30.0, 2.0);
        let v = sp.plan_speed(0.1); // 10 m radius
        assert!(v < 30.0, "speed should be limited in curve: {v}");
        // expected: sqrt(2.0 / 0.1) ≈ 4.47 m/s
        assert!(approx_eq(v, (2.0_f64 / 0.1).sqrt(), 1e-6));
    }

    #[test]
    fn speed_planner_tighter_curve_slower() {
        let sp = SpeedPlanner::new(3.0, 30.0, 2.0);
        let v_gentle = sp.plan_speed(0.05);
        let v_tight = sp.plan_speed(0.2);
        assert!(v_tight < v_gentle);
    }

    #[test]
    fn speed_planner_never_exceeds_max() {
        let sp = SpeedPlanner::new(3.0, 20.0, 100.0);
        assert!(approx_eq(sp.plan_speed(0.0), 20.0, 1e-10));
        assert!(approx_eq(sp.plan_speed(0.001), 20.0, 1e-3));
    }

    #[test]
    fn speed_planner_path_speeds_straight() {
        let sp = SpeedPlanner::new(3.0, 20.0, 2.0);
        let p = straight_path(50.0, 6);
        let speeds = sp.plan_path_speeds(&p);
        assert_eq!(speeds.len(), 6);
        for v in speeds {
            assert!(
                approx_eq(v, 20.0, 1e-8),
                "straight path should give max speed {v}"
            );
        }
    }

    #[test]
    fn speed_planner_circle_path_speeds_limited() {
        let sp = SpeedPlanner::new(3.0, 30.0, 2.0);
        let p = circle_path(10.0, 20);
        let speeds = sp.plan_path_speeds(&p);
        for v in &speeds {
            assert!(*v < 30.0, "speed should be limited on circle: {v}");
        }
    }

    // ── three_point_curvature ─────────────────────────────────────────────────

    #[test]
    fn three_point_curvature_collinear_is_zero() {
        let k = three_point_curvature([0.0, 0.0], [1.0, 0.0], [2.0, 0.0]);
        assert!(approx_eq(k, 0.0, 1e-12));
    }

    #[test]
    fn three_point_curvature_unit_circle() {
        // Three equidistant points on unit circle → curvature ≈ 1
        let theta = 2.0 * PI / 3.0;
        let p0 = [theta.cos(), theta.sin()];
        let p1 = [1.0, 0.0];
        let p2 = [(2.0 * theta).cos(), (2.0 * theta).sin()];
        let k = three_point_curvature(p0, p1, p2);
        assert!(
            approx_eq(k, 1.0, 0.02),
            "curvature ≈ 1 for unit circle, got {k}"
        );
    }
}
