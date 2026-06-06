//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CubicSplinePath, FrenetFrame, Path, Waypoint};

/// Compute the curvature (1/m) at each waypoint using three-point finite differences.
///
/// Uses the formula: kappa = |dx' * dy'' - dy' * dx''| / (dx'^2 + dy'^2)^(3/2)
/// where primes denote finite differences with respect to the waypoint index.
///
/// Returns a vector of curvatures (same length as the waypoints list).
pub fn compute_path_curvature(waypoints: &[[f64; 2]]) -> Vec<f64> {
    let n = waypoints.len();
    if n < 3 {
        return vec![0.0; n];
    }
    let mut kappa = vec![0.0; n];
    for i in 1..n - 1 {
        let xm = waypoints[i - 1][0];
        let ym = waypoints[i - 1][1];
        let x0 = waypoints[i][0];
        let y0 = waypoints[i][1];
        let xp = waypoints[i + 1][0];
        let yp = waypoints[i + 1][1];
        let dx = 0.5 * (xp - xm);
        let dy = 0.5 * (yp - ym);
        let ddx = xp - 2.0 * x0 + xm;
        let ddy = yp - 2.0 * y0 + ym;
        let denom = (dx * dx + dy * dy).powf(1.5);
        if denom > 1e-15 {
            kappa[i] = (dx * ddy - dy * ddx).abs() / denom;
        }
    }
    kappa[0] = kappa[1];
    kappa[n - 1] = kappa[n - 2];
    kappa
}
/// Compute signed curvature (positive = left turn in standard orientation).
pub fn compute_signed_curvature(waypoints: &[[f64; 2]]) -> Vec<f64> {
    let n = waypoints.len();
    if n < 3 {
        return vec![0.0; n];
    }
    let mut kappa = vec![0.0; n];
    for i in 1..n - 1 {
        let xm = waypoints[i - 1][0];
        let ym = waypoints[i - 1][1];
        let x0 = waypoints[i][0];
        let y0 = waypoints[i][1];
        let xp = waypoints[i + 1][0];
        let yp = waypoints[i + 1][1];
        let dx = 0.5 * (xp - xm);
        let dy = 0.5 * (yp - ym);
        let ddx = xp - 2.0 * x0 + xm;
        let ddy = yp - 2.0 * y0 + ym;
        let denom = (dx * dx + dy * dy).powf(1.5);
        if denom > 1e-15 {
            kappa[i] = (dx * ddy - dy * ddx) / denom;
        }
    }
    kappa[0] = kappa[1];
    kappa[n - 1] = kappa[n - 2];
    kappa
}
/// Compute a maximum cornering speed for a given curvature and friction.
///
/// v_max = sqrt(mu * g / kappa), clamped to `speed_limit`.
pub fn max_cornering_speed(curvature: f64, friction: f64, gravity: f64, speed_limit: f64) -> f64 {
    if curvature.abs() < 1e-12 {
        return speed_limit;
    }
    let v = (friction * gravity / curvature.abs()).sqrt();
    v.min(speed_limit)
}
/// Simple speed profile along a path: for each waypoint, compute the maximum
/// speed limited by the local curvature and friction.
///
/// Returns a vector of speeds (m/s) for each waypoint.
pub fn curvature_limited_speed_profile(
    waypoints: &[[f64; 2]],
    friction: f64,
    gravity: f64,
    speed_limit: f64,
) -> Vec<f64> {
    let curvatures = compute_path_curvature(waypoints);
    curvatures
        .iter()
        .map(|&k| max_cornering_speed(k, friction, gravity, speed_limit))
        .collect()
}
/// Forward-backward pass to smooth the speed profile accounting for
/// acceleration and braking limits.
///
/// This is a simplified version of the optimal racing line speed integration.
pub fn smooth_speed_profile(
    speeds: &[f64],
    segment_lengths: &[f64],
    max_accel: f64,
    max_decel: f64,
) -> Vec<f64> {
    let n = speeds.len();
    if n == 0 {
        return vec![];
    }
    let mut result = speeds.to_vec();
    for i in 1..n {
        let ds = if i - 1 < segment_lengths.len() {
            segment_lengths[i - 1]
        } else {
            1.0
        };
        let v_max_sq = result[i - 1] * result[i - 1] + 2.0 * max_accel * ds;
        let v_max = v_max_sq.max(0.0).sqrt();
        result[i] = result[i].min(v_max);
    }
    for i in (0..n - 1).rev() {
        let ds = if i < segment_lengths.len() {
            segment_lengths[i]
        } else {
            1.0
        };
        let v_max_sq = result[i + 1] * result[i + 1] + 2.0 * max_decel * ds;
        let v_max = v_max_sq.max(0.0).sqrt();
        result[i] = result[i].min(v_max);
    }
    result
}
/// Smooth a 2-D path using iterative averaging (Chaikin-like).
///
/// Each iteration replaces interior points with a weighted average of
/// their neighbours. The `weight` parameter (0..1) controls how much
/// smoothing is applied per iteration.
pub fn smooth_path_iterative(points: &[[f64; 2]], iterations: usize, weight: f64) -> Vec<[f64; 2]> {
    let n = points.len();
    if n < 3 {
        return points.to_vec();
    }
    let w = weight.clamp(0.0, 1.0);
    let mut current = points.to_vec();
    for _ in 0..iterations {
        let mut next = current.clone();
        for i in 1..n - 1 {
            next[i][0] =
                current[i][0] + w * (0.5 * (current[i - 1][0] + current[i + 1][0]) - current[i][0]);
            next[i][1] =
                current[i][1] + w * (0.5 * (current[i - 1][1] + current[i + 1][1]) - current[i][1]);
        }
        current = next;
    }
    current
}
/// Smooth with a penalty on deviation from the original path.
///
/// Balances smoothness against fidelity to the original trajectory using
/// two weights: `alpha` (data fidelity) and `beta` (smoothness).
pub fn smooth_path_with_penalty(
    original: &[[f64; 2]],
    iterations: usize,
    alpha: f64,
    beta: f64,
) -> Vec<[f64; 2]> {
    let n = original.len();
    if n < 3 {
        return original.to_vec();
    }
    let mut current = original.to_vec();
    for _ in 0..iterations {
        let prev = current.clone();
        for i in 1..n - 1 {
            for dim in 0..2 {
                let data_term = alpha * (original[i][dim] - prev[i][dim]);
                let smooth_term = beta * (prev[i - 1][dim] + prev[i + 1][dim] - 2.0 * prev[i][dim]);
                current[i][dim] = prev[i][dim] + data_term + smooth_term;
            }
        }
    }
    current
}
/// Linearly interpolate between two waypoints by parameter `t` in \[0, 1\].
pub fn interpolate_waypoints(a: &Waypoint, b: &Waypoint, t: f64) -> Waypoint {
    let t = t.clamp(0.0, 1.0);
    Waypoint {
        position: [
            a.position[0] + t * (b.position[0] - a.position[0]),
            a.position[1] + t * (b.position[1] - a.position[1]),
            a.position[2] + t * (b.position[2] - a.position[2]),
        ],
        speed_limit: a.speed_limit + t * (b.speed_limit - a.speed_limit),
        heading: a.heading + t * (b.heading - a.heading),
    }
}
/// Resample a path to have uniformly spaced waypoints.
///
/// Returns a new `Path` with `n_points` waypoints evenly distributed
/// along the arc length of the original path.
pub fn resample_path(path: &Path, n_points: usize) -> Path {
    let mut resampled = Path::new(path.closed);
    if path.waypoints.is_empty() || n_points == 0 {
        return resampled;
    }
    if n_points == 1 || path.waypoints.len() == 1 {
        resampled.add_waypoint(path.waypoints[0].position, path.waypoints[0].speed_limit);
        return resampled;
    }
    let total = path.total_length();
    for i in 0..n_points {
        let d = total * (i as f64) / (n_points - 1) as f64;
        if let Some(pos) = path.interpolate_at_distance(d) {
            let (idx, _) = path.nearest_waypoint(pos);
            let speed = path.waypoints[idx].speed_limit;
            resampled.add_waypoint(pos, speed);
        }
    }
    resampled
}
/// Compute speed-dependent look-ahead distance for pure pursuit.
///
/// `l_d = l_min + k * v`
///
/// Clamped to `[l_min, l_max]`.
pub fn speed_dependent_lookahead(speed: f64, k: f64, l_min: f64, l_max: f64) -> f64 {
    (l_min + k * speed.abs()).clamp(l_min, l_max)
}
/// Compute curvature-dependent look-ahead distance.
///
/// At high curvature, use shorter look-ahead to track tight corners.
/// At low curvature, extend look-ahead for stability.
///
/// `l_d = l_base / (1 + k_curv * |kappa|)`
pub fn curvature_dependent_lookahead(curvature: f64, l_base: f64, k_curv: f64, l_min: f64) -> f64 {
    let l = l_base / (1.0 + k_curv * curvature.abs());
    l.max(l_min)
}
/// Signed lateral (cross-track) distance from `vehicle_pos` to a path point.
pub fn cross_track_error(vehicle_pos: [f64; 3], path_pos: [f64; 3], path_heading: f64) -> f64 {
    let dx = vehicle_pos[0] - path_pos[0];
    let dy = vehicle_pos[1] - path_pos[1];
    -dx * path_heading.sin() + dy * path_heading.cos()
}
/// Estimate minimum lap time along a `CubicSplinePath`.
///
/// Uses a simplified velocity-profile integration:
/// 1. Compute the lateral-g limited speed at each sample: `v = √(a_lat * R)`
/// 2. Limit by braking: backward pass using `v² = v_next² + 2 * a_brake * ds`
/// 3. Limit by power: `v_power = P / (F = m * a_lat)`, simplified as max power
///    constraint `v ≤ P / (m * a_lat_g)`.
/// 4. Integrate `dt = ds / v` over the path.
pub fn minimum_lap_time_estimate(
    path: &CubicSplinePath,
    max_lateral_g: f64,
    max_brake_g: f64,
    max_power: f64,
    mass: f64,
) -> f64 {
    let g = 9.81;
    let n = 200usize;
    if path.total_length < 1e-9 || n < 2 {
        return 0.0;
    }
    let ds = path.total_length / (n - 1) as f64;
    let mut v: Vec<f64> = (0..n)
        .map(|k| {
            let t = k as f64 / (n - 1) as f64;
            let kappa = path.curvature_at(t).max(1e-9);
            let r = 1.0 / kappa;
            let v_lat = (max_lateral_g * g * r).sqrt();
            let _ = max_power;
            let _ = mass;
            v_lat
        })
        .collect();
    let a_brake = max_brake_g * g;
    for k in (0..n - 1).rev() {
        let v_max = (v[k + 1] * v[k + 1] + 2.0 * a_brake * ds).sqrt();
        if v[k] > v_max {
            v[k] = v_max;
        }
    }
    let a_drive = if mass > 0.0 {
        max_power / (mass * v[0].max(0.1))
    } else {
        0.0
    };
    let _ = a_drive;
    for k in 0..n - 1 {
        let v_cur = v[k];
        let f_drive = if v_cur > 0.01 {
            max_power / v_cur
        } else {
            max_power / 0.01
        };
        let a_fwd = if mass > 0.0 { f_drive / mass } else { 0.0 };
        let v_max_fwd = (v_cur * v_cur + 2.0 * a_fwd * ds).sqrt();
        if v[k + 1] > v_max_fwd {
            v[k + 1] = v_max_fwd;
        }
    }
    let mut time = 0.0;
    for k in 0..n - 1 {
        let v_avg = ((v[k] + v[k + 1]) * 0.5).max(0.01);
        time += ds / v_avg;
    }
    time
}
/// Compute cross-track error statistics from a recorded trajectory vs a
/// reference path.
///
/// Returns `(mean_error, max_error, rms_error)` in metres.
pub fn cross_track_statistics(trajectory: &[[f64; 2]], reference: &[[f64; 2]]) -> (f64, f64, f64) {
    if trajectory.is_empty() || reference.len() < 2 {
        return (0.0, 0.0, 0.0);
    }
    let mut sum = 0.0f64;
    let mut sum_sq = 0.0f64;
    let mut max = 0.0f64;
    for tp in trajectory {
        let d = min_dist_to_polyline(tp, reference);
        sum += d;
        sum_sq += d * d;
        if d > max {
            max = d;
        }
    }
    let n = trajectory.len() as f64;
    let mean = sum / n;
    let rms = (sum_sq / n).sqrt();
    (mean, max, rms)
}
/// Minimum Euclidean distance from a 2-D point to a polyline.
pub(super) fn min_dist_to_polyline(pt: &[f64; 2], poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    if n == 0 {
        return f64::INFINITY;
    }
    if n == 1 {
        let dx = pt[0] - poly[0][0];
        let dy = pt[1] - poly[0][1];
        return (dx * dx + dy * dy).sqrt();
    }
    let mut min_d2 = f64::INFINITY;
    for i in 0..n - 1 {
        let ax = poly[i][0];
        let ay = poly[i][1];
        let bx = poly[i + 1][0];
        let by = poly[i + 1][1];
        let abx = bx - ax;
        let aby = by - ay;
        let len2 = abx * abx + aby * aby;
        let t = if len2 > 1e-15 {
            ((pt[0] - ax) * abx + (pt[1] - ay) * aby) / len2
        } else {
            0.0
        }
        .clamp(0.0, 1.0);
        let cx = ax + t * abx - pt[0];
        let cy = ay + t * aby - pt[1];
        let d2 = cx * cx + cy * cy;
        if d2 < min_d2 {
            min_d2 = d2;
        }
    }
    min_d2.sqrt()
}
/// Inline 3-D cross product.
pub(super) fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Inline 3-D vector length.
pub fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Normalize a 3-D vector; returns `[1,0,0]` if near-zero.
pub(super) fn vec3_norm(v: [f64; 3]) -> [f64; 3] {
    let len = vec3_len(v);
    if len < 1e-15 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}
/// Reparameterise a set of 3-D polyline points by arc-length.
///
/// Returns a vector of cumulative arc-length values corresponding to each
/// input point (same length as `pts`).
pub fn arc_length_parameter(pts: &[[f64; 3]]) -> Vec<f64> {
    let n = pts.len();
    if n == 0 {
        return vec![];
    }
    let mut s = vec![0.0_f64; n];
    for i in 1..n {
        let dx = pts[i][0] - pts[i - 1][0];
        let dy = pts[i][1] - pts[i - 1][1];
        let dz = pts[i][2] - pts[i - 1][2];
        s[i] = s[i - 1] + (dx * dx + dy * dy + dz * dz).sqrt();
    }
    s
}
/// Sample a polyline at a specific arc-length `target_s`.
///
/// Linearly interpolates between the surrounding points.
/// Clamps to the first/last point if `target_s` is out of range.
pub fn sample_polyline_at_arc_length(pts: &[[f64; 3]], target_s: f64) -> [f64; 3] {
    let n = pts.len();
    if n == 0 {
        return [0.0; 3];
    }
    if n == 1 {
        return pts[0];
    }
    let s = arc_length_parameter(pts);
    if target_s <= 0.0 {
        return pts[0];
    }
    if target_s >= s[n - 1] {
        return pts[n - 1];
    }
    let idx = s
        .partition_point(|&si| si < target_s)
        .saturating_sub(1)
        .min(n - 2);
    let seg_len = s[idx + 1] - s[idx];
    let t = if seg_len > 1e-15 {
        (target_s - s[idx]) / seg_len
    } else {
        0.0
    };
    let a = pts[idx];
    let b = pts[idx + 1];
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}
/// Generate a clothoid (Euler spiral) as a sequence of 2-D points.
///
/// A clothoid satisfies κ(s) = s / A², where A is the clothoid parameter
/// (units: √m).  The curvature increases linearly from zero, making it ideal
/// for smooth transitions from straight to curved sections.
///
/// # Arguments
/// * `kappa_start` – curvature at the start (1/m; 0 for entry clothoid)
/// * `a_sq`        – clothoid parameter A² (m²); larger → gentler transition
/// * `heading0`    – initial heading (radians from +X)
/// * `ds`          – arc-length step between sample points (m)
/// * `n_pts`       – number of output points
///
/// Returns a `Vec` of 2-D points `[x, y]` starting at the origin.
pub fn clothoid_points(
    kappa_start: f64,
    a_sq: f64,
    heading0: f64,
    ds: f64,
    n_pts: usize,
) -> Vec<[f64; 2]> {
    if n_pts == 0 || ds <= 0.0 || a_sq <= 0.0 {
        return vec![];
    }
    let mut pts = Vec::with_capacity(n_pts);
    let mut x = 0.0_f64;
    let mut y = 0.0_f64;
    let mut theta = heading0;
    let mut kappa = kappa_start;
    pts.push([x, y]);
    for _ in 1..n_pts {
        let dkappa = 1.0 / a_sq;
        theta += kappa * ds;
        kappa += dkappa * ds;
        x += ds * theta.cos();
        y += ds * theta.sin();
        pts.push([x, y]);
    }
    pts
}
/// Radius of curvature at arc-length `s` along a clothoid with parameter A².
///
/// R(s) = A² / s  (tends to ∞ at s=0, decreases linearly).
pub fn clothoid_radius_at(a_sq: f64, s: f64) -> f64 {
    if s.abs() < 1e-15 {
        f64::INFINITY
    } else {
        a_sq / s.abs()
    }
}
/// Approximate the clothoid length required to reach radius `target_radius` (m)
/// starting from a given parameter A².
///
/// `s = A² / R`
pub fn clothoid_length_for_radius(a_sq: f64, target_radius: f64) -> f64 {
    if target_radius <= 0.0 {
        return 0.0;
    }
    a_sq / target_radius
}
/// Compute the signed heading error between two angles, wrapped to `(-π, π]`.
///
/// Returns `desired_heading - actual_heading` normalised to `(-π, π]`.
pub fn heading_error_wrapped(actual: f64, desired: f64) -> f64 {
    let mut err = desired - actual;
    while err > std::f64::consts::PI {
        err -= 2.0 * std::f64::consts::PI;
    }
    while err <= -std::f64::consts::PI {
        err += 2.0 * std::f64::consts::PI;
    }
    err
}
/// Project a world position onto a `CubicSplinePath` and return both the
/// signed cross-track error and heading error.
///
/// # Arguments
/// * `path`         – reference path
/// * `vehicle_pos`  – vehicle 3-D world position
/// * `vehicle_yaw`  – vehicle heading (radians from +X)
/// * `n_samples`    – number of uniform samples to use in the projection search
///
/// # Returns
/// `(cross_track_error_m, heading_error_rad)`
/// CTE is positive when the vehicle is to the left of the path (in the
/// path-normal frame).
pub fn frenet_errors(
    path: &CubicSplinePath,
    vehicle_pos: [f64; 3],
    vehicle_yaw: f64,
    n_samples: usize,
) -> (f64, f64) {
    let n = n_samples.max(10);
    let mut best_t = 0.0_f64;
    let mut best_d2 = f64::INFINITY;
    for k in 0..n {
        let t = k as f64 / (n - 1) as f64;
        let p = path.sample(t);
        let dx = p[0] - vehicle_pos[0];
        let dy = p[1] - vehicle_pos[1];
        let dz = p[2] - vehicle_pos[2];
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < best_d2 {
            best_d2 = d2;
            best_t = t;
        }
    }
    let frame = FrenetFrame::compute(path, best_t);
    let p_on_path = path.sample(best_t);
    let diff = [
        vehicle_pos[0] - p_on_path[0],
        vehicle_pos[1] - p_on_path[1],
        vehicle_pos[2] - p_on_path[2],
    ];
    let cte = diff[0] * frame.normal[0] + diff[1] * frame.normal[1] + diff[2] * frame.normal[2];
    let path_yaw = frame.tangent[1].atan2(frame.tangent[0]);
    let he = heading_error_wrapped(vehicle_yaw, path_yaw);
    (cte, he)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::path::ArcLengthController;
    use crate::path::CubicSpline2D;
    use crate::path::MinCurvatureRacingLine;
    use crate::path::PurePursuitController;
    use crate::path::RacingLine;
    use crate::path::StanleyController;
    use crate::path::TrackLayout;
    use crate::path::TrackSegment;
    use std::f64::consts::PI;
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    #[test]
    fn test_path_total_length_two_waypoints() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([3.0, 4.0, 0.0], 10.0);
        assert!(approx_eq(path.total_length(), 5.0, 1e-10));
    }
    #[test]
    fn test_nearest_waypoint() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([10.0, 0.0, 0.0], 10.0);
        path.add_waypoint([20.0, 0.0, 0.0], 10.0);
        let (idx, dist) = path.nearest_waypoint([9.0, 1.0, 0.0]);
        assert_eq!(idx, 1);
        assert!(approx_eq(dist, (1.0_f64 + 1.0_f64).sqrt(), 1e-10));
    }
    #[test]
    fn test_pure_pursuit_target_directly_ahead() {
        let ctrl = PurePursuitController::new(5.0, 2.5);
        let vehicle_pos = [0.0, 0.0, 0.0];
        let vehicle_heading = 0.0;
        let target = [5.0, 0.0, 0.0];
        let angle = ctrl.steering_angle(vehicle_pos, vehicle_heading, target);
        assert!(approx_eq(angle, 0.0, 1e-10));
    }
    #[test]
    fn test_stanley_zero_cte_zero_heading_error() {
        let ctrl = StanleyController::new(1.0);
        let angle = ctrl.steering_angle(0.0, 0.0, 5.0);
        assert!(approx_eq(angle, 0.0, 1e-10));
    }
    #[test]
    fn test_interpolate_at_distance_endpoints() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([1.0, 0.0, 0.0], 10.0);
        path.add_waypoint([2.0, 0.0, 0.0], 10.0);
        let start = path.interpolate_at_distance(0.0).unwrap();
        assert!(approx_eq(start[0], 0.0, 1e-10));
        let total = path.total_length();
        let end = path.interpolate_at_distance(total).unwrap();
        assert!(approx_eq(end[0], 2.0, 1e-10));
    }
    #[test]
    fn test_cross_track_error_on_path() {
        let cte = cross_track_error([5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.0);
        assert!(approx_eq(cte, 0.0, 1e-10));
    }
    #[test]
    fn test_cross_track_error_offset() {
        let cte = cross_track_error([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], 0.0);
        assert!(approx_eq(cte, 1.0, 1e-10));
    }
    #[test]
    fn test_pure_pursuit_curvature_zero_alpha() {
        let ctrl = PurePursuitController::new(5.0, 2.5);
        let kappa = ctrl.curvature([0.0, 0.0, 0.0], 0.0, [5.0, 0.0, 0.0]);
        assert!(approx_eq(kappa, 0.0, 1e-10));
    }
    #[test]
    fn test_stanley_nonzero_heading_error() {
        let ctrl = StanleyController::new(1.0);
        let he = PI / 6.0;
        let angle = ctrl.steering_angle(0.0, he, 10.0);
        assert!(approx_eq(angle, he, 1e-10));
    }
    #[test]
    fn test_cubic_spline_passes_through_points() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 0.0, 1.0];
        let spline = CubicSpline2D::new(&xs, &ys).unwrap();
        for i in 0..xs.len() {
            let p = spline.evaluate(i as f64);
            assert!(
                approx_eq(p[0], xs[i], 1e-6),
                "x mismatch at t={i}: got {}, expected {}",
                p[0],
                xs[i]
            );
            assert!(
                approx_eq(p[1], ys[i], 1e-6),
                "y mismatch at t={i}: got {}, expected {}",
                p[1],
                ys[i]
            );
        }
    }
    #[test]
    fn test_cubic_spline_midpoint_reasonable() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 0.0];
        let spline = CubicSpline2D::new(&xs, &ys).unwrap();
        let mid = spline.evaluate(0.5);
        assert!(mid[0] > 0.0 && mid[0] < 1.0);
        assert!(mid[1] > 0.0 && mid[1] < 1.0);
    }
    #[test]
    fn test_cubic_spline_tangent() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 0.0, 0.0, 0.0];
        let spline = CubicSpline2D::new(&xs, &ys).unwrap();
        let t = spline.tangent(1.0);
        assert!(approx_eq(t[0], 1.0, 1e-6));
        assert!(approx_eq(t[1], 0.0, 1e-6));
    }
    #[test]
    fn test_cubic_spline_clamped_endpoints() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 4.0];
        let spline = CubicSpline2D::new(&xs, &ys).unwrap();
        let left = spline.evaluate(-1.0);
        assert!(approx_eq(left[0], 0.0, 1e-10));
        let right = spline.evaluate(10.0);
        assert!(approx_eq(right[0], 2.0, 1e-10));
    }
    #[test]
    fn test_cubic_spline_len() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 4.0];
        let spline = CubicSpline2D::new(&xs, &ys).unwrap();
        assert_eq!(spline.len(), 3);
        assert!(!spline.is_empty());
    }
    #[test]
    fn test_curvature_straight_line() {
        let pts: Vec<[f64; 2]> = (0..10).map(|i| [i as f64, 0.0]).collect();
        let k = compute_path_curvature(&pts);
        for &ki in &k {
            assert!(ki.abs() < 1e-10, "straight line should have zero curvature");
        }
    }
    #[test]
    fn test_curvature_circle() {
        let r = 10.0;
        let n = 100;
        let pts: Vec<[f64; 2]> = (0..n)
            .map(|i| {
                let theta = 2.0 * PI * (i as f64) / (n as f64);
                [r * theta.cos(), r * theta.sin()]
            })
            .collect();
        let k = compute_path_curvature(&pts);
        for (i, &ki) in k.iter().enumerate().take(n - 5).skip(5) {
            assert!(
                approx_eq(ki, 1.0 / r, 0.05 / r),
                "curvature at {i}: got {}, expected {}",
                ki,
                1.0 / r
            );
        }
    }
    #[test]
    fn test_signed_curvature_direction() {
        let r = 10.0;
        let n = 50;
        let pts: Vec<[f64; 2]> = (0..n)
            .map(|i| {
                let theta = 2.0 * PI * (i as f64) / (n as f64);
                [r * theta.cos(), r * theta.sin()]
            })
            .collect();
        let k = compute_signed_curvature(&pts);
        for &ki in k.iter().take(n - 5).skip(5) {
            assert!(ki > 0.0, "CCW circle should have positive signed curvature");
        }
    }
    #[test]
    fn test_max_cornering_speed_straight() {
        let v = max_cornering_speed(0.0, 1.0, 9.81, 100.0);
        assert!((v - 100.0).abs() < 1e-10, "straight → speed limit");
    }
    #[test]
    fn test_max_cornering_speed_curve() {
        let kappa = 1.0 / 50.0;
        let v = max_cornering_speed(kappa, 1.0, 9.81, 200.0);
        let expected = (1.0_f64 * 9.81 * 50.0).sqrt();
        assert!(approx_eq(v, expected, 1e-6));
    }
    #[test]
    fn test_curvature_limited_speed_profile() {
        let pts: Vec<[f64; 2]> = (0..10).map(|i| [i as f64, 0.0]).collect();
        let speeds = curvature_limited_speed_profile(&pts, 1.0, 9.81, 50.0);
        assert_eq!(speeds.len(), 10);
        for &s in &speeds {
            assert!(approx_eq(s, 50.0, 1e-6));
        }
    }
    #[test]
    fn test_smooth_speed_profile() {
        let speeds = vec![50.0, 50.0, 50.0, 20.0, 50.0, 50.0];
        let seg_lens = vec![10.0, 10.0, 10.0, 10.0, 10.0];
        let smoothed = smooth_speed_profile(&speeds, &seg_lens, 5.0, 5.0);
        assert!(smoothed[2] < 50.0, "should slow down before the dip");
        assert!(smoothed[4] < 50.0, "should accelerate out of the dip");
    }
    #[test]
    fn test_smooth_path_preserves_endpoints() {
        let pts = vec![[0.0, 0.0], [1.0, 2.0], [3.0, 1.0], [4.0, 0.0]];
        let smoothed = smooth_path_iterative(&pts, 10, 0.5);
        assert!(approx_eq(smoothed[0][0], 0.0, 1e-10));
        assert!(approx_eq(smoothed[0][1], 0.0, 1e-10));
        assert!(approx_eq(smoothed[3][0], 4.0, 1e-10));
        assert!(approx_eq(smoothed[3][1], 0.0, 1e-10));
    }
    #[test]
    fn test_smooth_path_reduces_curvature() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0], [2.0, -1.0], [3.0, 1.0], [4.0, 0.0]];
        let original_k = compute_path_curvature(&pts);
        let smoothed = smooth_path_iterative(&pts, 20, 0.5);
        let smoothed_k = compute_path_curvature(&smoothed);
        let orig_max: f64 = original_k.iter().cloned().fold(0.0, f64::max);
        let smooth_max: f64 = smoothed_k.iter().cloned().fold(0.0, f64::max);
        assert!(
            smooth_max < orig_max,
            "smoothing should reduce max curvature"
        );
    }
    #[test]
    fn test_smooth_path_with_penalty() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0], [2.0, -1.0], [3.0, 1.0], [4.0, 0.0]];
        let smoothed = smooth_path_with_penalty(&pts, 50, 0.1, 0.3);
        assert!(approx_eq(smoothed[0][0], 0.0, 1e-10));
        assert!(approx_eq(smoothed[4][0], 4.0, 1e-10));
        let orig_k = compute_path_curvature(&pts);
        let smooth_k = compute_path_curvature(&smoothed);
        let orig_var: f64 = orig_k.iter().map(|k| k * k).sum();
        let smooth_var: f64 = smooth_k.iter().map(|k| k * k).sum();
        assert!(
            smooth_var < orig_var,
            "penalty smoothing should reduce curvature energy"
        );
    }
    #[test]
    fn test_interpolate_waypoints_endpoints() {
        let a = Waypoint {
            position: [0.0, 0.0, 0.0],
            speed_limit: 10.0,
            heading: 0.0,
        };
        let b = Waypoint {
            position: [10.0, 0.0, 0.0],
            speed_limit: 20.0,
            heading: PI,
        };
        let start = interpolate_waypoints(&a, &b, 0.0);
        assert!(approx_eq(start.position[0], 0.0, 1e-10));
        let end = interpolate_waypoints(&a, &b, 1.0);
        assert!(approx_eq(end.position[0], 10.0, 1e-10));
    }
    #[test]
    fn test_interpolate_waypoints_midpoint() {
        let a = Waypoint {
            position: [0.0, 0.0, 0.0],
            speed_limit: 10.0,
            heading: 0.0,
        };
        let b = Waypoint {
            position: [10.0, 0.0, 0.0],
            speed_limit: 20.0,
            heading: 1.0,
        };
        let mid = interpolate_waypoints(&a, &b, 0.5);
        assert!(approx_eq(mid.position[0], 5.0, 1e-10));
        assert!(approx_eq(mid.speed_limit, 15.0, 1e-10));
        assert!(approx_eq(mid.heading, 0.5, 1e-10));
    }
    #[test]
    fn test_resample_path() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([10.0, 0.0, 0.0], 10.0);
        let resampled = resample_path(&path, 5);
        assert_eq!(resampled.waypoints.len(), 5);
        assert!(approx_eq(resampled.waypoints[0].position[0], 0.0, 1e-10));
        assert!(approx_eq(resampled.waypoints[4].position[0], 10.0, 1e-10));
    }
    #[test]
    fn test_speed_dependent_lookahead() {
        let l = speed_dependent_lookahead(10.0, 0.5, 3.0, 20.0);
        assert!(approx_eq(l, 8.0, 1e-10));
    }
    #[test]
    fn test_speed_dependent_lookahead_clamped() {
        let l = speed_dependent_lookahead(100.0, 0.5, 3.0, 20.0);
        assert!(approx_eq(l, 20.0, 1e-10), "should clamp to l_max");
    }
    #[test]
    fn test_curvature_dependent_lookahead() {
        let l = curvature_dependent_lookahead(0.0, 10.0, 5.0, 2.0);
        assert!(approx_eq(l, 10.0, 1e-10), "zero curvature → l_base");
    }
    #[test]
    fn test_curvature_dependent_lookahead_tight_curve() {
        let l = curvature_dependent_lookahead(1.0, 10.0, 5.0, 2.0);
        assert!(approx_eq(l, 2.0, 1e-10));
    }
    #[test]
    fn test_cumulative_arc_lengths() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([3.0, 4.0, 0.0], 10.0);
        path.add_waypoint([6.0, 8.0, 0.0], 10.0);
        let arcs = path.cumulative_arc_lengths();
        assert_eq!(arcs.len(), 3);
        assert!(approx_eq(arcs[0], 0.0, 1e-10));
        assert!(approx_eq(arcs[1], 5.0, 1e-10));
        assert!(approx_eq(arcs[2], 10.0, 1e-10));
    }
    #[test]
    fn test_compute_headings() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([1.0, 0.0, 0.0], 10.0);
        path.add_waypoint([2.0, 0.0, 0.0], 10.0);
        path.compute_headings();
        for wp in &path.waypoints {
            assert!(approx_eq(wp.heading, 0.0, 1e-10));
        }
    }
    #[test]
    fn test_compute_headings_diagonal() {
        let mut path = Path::new(false);
        path.add_waypoint([0.0, 0.0, 0.0], 10.0);
        path.add_waypoint([1.0, 1.0, 0.0], 10.0);
        path.compute_headings();
        assert!(approx_eq(path.waypoints[0].heading, PI / 4.0, 1e-10));
    }
    #[test]
    fn test_cubic_spline_path_from_waypoints_count() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        assert_eq!(sp.waypoints.len(), 5, "should have 5 waypoints");
        assert_eq!(sp.tangents.len(), 5, "should have 5 tangents");
    }
    #[test]
    fn test_cubic_spline_path_sample_t0_first_point() {
        let pts: &[[f64; 3]] = &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let sp = CubicSplinePath::from_waypoints(pts);
        let p = sp.sample(0.0);
        assert!(
            approx_eq(p[0], 0.0, 1e-6),
            "t=0 should give first point x=0, got {}",
            p[0]
        );
        assert!(approx_eq(p[1], 0.0, 1e-6));
    }
    #[test]
    fn test_cubic_spline_path_sample_t1_last_point() {
        let pts: &[[f64; 3]] = &[
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 3.0, 0.0],
        ];
        let sp = CubicSplinePath::from_waypoints(pts);
        let p = sp.sample(1.0);
        let last = pts[pts.len() - 1];
        assert!(
            approx_eq(p[0], last[0], 1e-5),
            "t=1 should give last point x={}, got {}",
            last[0],
            p[0]
        );
        assert!(
            approx_eq(p[1], last[1], 1e-5),
            "t=1 should give last point y={}, got {}",
            last[1],
            p[1]
        );
    }
    #[test]
    fn test_cubic_spline_path_curvature_straight_near_zero() {
        let pts: Vec<[f64; 3]> = (0..6).map(|i| [i as f64 * 2.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        for k in 0..=10 {
            let t = k as f64 / 10.0;
            let kappa = sp.curvature_at(t);
            assert!(
                kappa < 1e-4,
                "straight line curvature should be ~0, got {kappa} at t={t}"
            );
        }
    }
    #[test]
    fn test_cubic_spline_path_total_length_straight() {
        let pts: &[[f64; 3]] = &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let sp = CubicSplinePath::from_waypoints(pts);
        assert!(
            approx_eq(sp.total_length, 3.0, 0.05),
            "total_length={}",
            sp.total_length
        );
    }
    #[test]
    fn test_racing_line_midpoints() {
        let left: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 1.0, 0.0]).collect();
        let right: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, -1.0, 0.0]).collect();
        let rl = RacingLine::from_track_edges(&left, &right, 5);
        assert_eq!(rl.points.len(), 5);
        for pt in &rl.points {
            assert!(
                approx_eq(pt[1], 0.0, 1e-9),
                "midpoint y should be 0, got {}",
                pt[1]
            );
        }
    }
    #[test]
    fn test_minimum_lap_time_positive() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        let t = minimum_lap_time_estimate(&sp, 1.5, 1.5, 150_000.0, 800.0);
        assert!(t > 0.0, "lap time must be positive, got {t}");
    }
    #[test]
    fn test_straight_segment_length() {
        let s = TrackSegment::Straight {
            start: [0.0, 0.0],
            end: [100.0, 0.0],
        };
        assert!((s.length() - 100.0).abs() < 1e-9);
    }
    #[test]
    fn test_arc_segment_length() {
        let r = 20.0;
        let s = TrackSegment::Arc {
            center: [0.0, r],
            radius: r,
            start_angle: -PI / 2.0,
            sweep: PI / 2.0,
        };
        let expected = r * PI / 2.0;
        assert!(
            (s.length() - expected).abs() < 1e-6,
            "arc length={} expected={}",
            s.length(),
            expected
        );
    }
    #[test]
    fn test_chicane_segment_length_positive() {
        let s = TrackSegment::Chicane {
            waypoints: vec![
                [0.0, 0.0],
                [10.0, 5.0],
                [20.0, 0.0],
                [30.0, -5.0],
                [40.0, 0.0],
            ],
        };
        assert!(s.length() > 0.0, "chicane length must be positive");
    }
    #[test]
    fn test_track_layout_total_length() {
        let mut layout = TrackLayout::new("test_track");
        layout.add_segment(TrackSegment::Straight {
            start: [0.0, 0.0],
            end: [50.0, 0.0],
        });
        layout.add_segment(TrackSegment::Straight {
            start: [50.0, 0.0],
            end: [100.0, 0.0],
        });
        assert!((layout.total_length() - 100.0).abs() < 1e-9);
    }
    #[test]
    fn test_track_layout_num_segments() {
        let mut layout = TrackLayout::new("circuit");
        assert_eq!(layout.num_segments(), 0);
        layout.add_segment(TrackSegment::Straight {
            start: [0.0, 0.0],
            end: [10.0, 0.0],
        });
        assert_eq!(layout.num_segments(), 1);
    }
    #[test]
    fn test_min_curvature_racing_line_output_count() {
        let left: Vec<[f64; 2]> = (0..10).map(|i| [i as f64 * 5.0, 2.0]).collect();
        let right: Vec<[f64; 2]> = (0..10).map(|i| [i as f64 * 5.0, -2.0]).collect();
        let line = MinCurvatureRacingLine::optimize(&left, &right, 50);
        assert_eq!(
            line.points.len(),
            10,
            "should return one point per input row"
        );
    }
    #[test]
    fn test_min_curvature_alphas_in_range() {
        let left: Vec<[f64; 2]> = (0..8).map(|i| [i as f64 * 10.0, 3.0]).collect();
        let right: Vec<[f64; 2]> = (0..8).map(|i| [i as f64 * 10.0, -3.0]).collect();
        let line = MinCurvatureRacingLine::optimize(&left, &right, 30);
        for &a in &line.alphas {
            assert!((0.0..=1.0).contains(&a), "alpha={a} out of [0,1]");
        }
    }
    #[test]
    fn test_arc_length_controller_full_lap_fraction_one() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        let ctrl = ArcLengthController::new(sp, 30.0, 5.0, 50.0);
        let frac = ctrl.arc_fraction(ctrl.path.total_length);
        assert!((frac - 1.0).abs() < 1e-9, "end of path → fraction=1");
    }
    #[test]
    fn test_arc_length_controller_zero_start() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        let ctrl = ArcLengthController::new(sp, 30.0, 5.0, 50.0);
        let frac = ctrl.arc_fraction(0.0);
        assert!((frac).abs() < 1e-9, "start → fraction=0");
    }
    #[test]
    fn test_arc_length_controller_target_speed_at_zero() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        let ctrl = ArcLengthController::new(sp, 30.0, 5.0, 50.0);
        let v = ctrl.target_speed_at_arc(0.0, 1.0, 9.81);
        assert!(v > 0.0, "target speed must be positive");
        assert!(v <= 50.0, "speed must not exceed v_max");
    }
    #[test]
    fn test_arc_length_controller_steering_ahead_on_line() {
        let pts: Vec<[f64; 3]> = (0..6).map(|i| [i as f64 * 5.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        let mut ctrl = ArcLengthController::new(sp, 5.0, 2.5, 30.0);
        let steer = ctrl.steering_command([0.0, 0.0, 0.0], 0.0, 0.0);
        assert!(
            steer.abs() < 0.5,
            "on-path steering should be near zero, got {steer}"
        );
    }
    #[test]
    fn test_frenet_frame_straight_line() {
        let pts: Vec<[f64; 3]> = (0..6).map(|i| [i as f64 * 5.0, 0.0, 0.0]).collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        let frenet = FrenetFrame::compute(&sp, 0.5);
        assert!(frenet.tangent[0] > 0.9, "tangent should be mostly +X");
        let n_len = vec3_len(frenet.normal);
        let b_len = vec3_len(frenet.binormal);
        assert!(
            (n_len - 1.0).abs() < 1e-9 || n_len < 1e-9,
            "normal should be unit"
        );
        assert!(
            (b_len - 1.0).abs() < 1e-9 || b_len < 1e-9,
            "binormal should be unit"
        );
    }
    #[test]
    fn test_frenet_tangent_is_unit_vector() {
        let pts: Vec<[f64; 3]> = (0..5)
            .map(|i| {
                let t = i as f64 * 0.4;
                [t.cos() * 10.0, t.sin() * 10.0, 0.0]
            })
            .collect();
        let sp = CubicSplinePath::from_waypoints(&pts);
        for k in 0..=5 {
            let t = k as f64 / 5.0;
            let frenet = FrenetFrame::compute(&sp, t);
            let len = vec3_len(frenet.tangent);
            assert!(
                (len - 1.0).abs() < 1e-9,
                "tangent at t={t} not unit: len={len}"
            );
        }
    }
    #[test]
    fn test_clothoid_start_is_origin() {
        let pts = clothoid_points(0.0, 20.0, 0.0, 0.1, 10);
        assert!((pts[0][0]).abs() < 1e-9);
        assert!((pts[0][1]).abs() < 1e-9);
    }
    #[test]
    fn test_clothoid_length_positive() {
        let pts = clothoid_points(0.0, 20.0, 0.0, 0.1, 20);
        assert!(pts.len() == 20, "should produce exactly n_pts points");
    }
    #[test]
    fn test_clothoid_curvature_at_end() {
        let a_sq = 100.0_f64;
        let s = 10.0_f64;
        let kappa = s / a_sq;
        assert!((kappa - 0.1).abs() < 1e-9);
    }
    #[test]
    fn test_cross_track_error_function_zero() {
        let err = cross_track_error([5.0, 3.0, 0.0], [5.0, 3.0, 0.0], PI / 4.0);
        assert!(err.abs() < 1e-9);
    }
    #[test]
    fn test_heading_error_wrap() {
        let err = heading_error_wrapped(0.1, 0.1 + 2.0 * PI);
        assert!(err.abs() < 1e-9, "wrapping should give ~0, got {err}");
    }
    #[test]
    fn test_heading_error_small_positive() {
        let err = heading_error_wrapped(0.3, 0.1);
        assert!(
            (err - (-0.2)).abs() < 1e-9,
            "heading error should be -0.2, got {err}"
        );
    }
    #[test]
    fn test_adaptive_lookahead_increases_with_speed() {
        let la1 = speed_dependent_lookahead(5.0, 0.5, 3.0, 30.0);
        let la2 = speed_dependent_lookahead(15.0, 0.5, 3.0, 30.0);
        assert!(la2 > la1, "lookahead should increase with speed");
    }
    #[test]
    fn test_adaptive_lookahead_min_clamp() {
        let la = speed_dependent_lookahead(0.0, 0.5, 5.0, 30.0);
        assert!((la - 5.0).abs() < 1e-9, "at zero speed lookahead = l_min");
    }
    #[test]
    fn test_stanley_front_axle_cte() {
        let ctrl = StanleyController {
            k: 2.0,
            k_soft: 1.0,
        };
        let angle = ctrl.steering_angle(1.0, 0.0, 10.0);
        assert!(angle > 0.0, "positive CTE → positive steering, got {angle}");
    }
    #[test]
    fn test_stanley_front_axle_high_speed_reduces_cte_term() {
        let ctrl = StanleyController {
            k: 2.0,
            k_soft: 1.0,
        };
        let angle_slow = ctrl.steering_angle(1.0, 0.0, 1.0);
        let angle_fast = ctrl.steering_angle(1.0, 0.0, 50.0);
        assert!(
            angle_fast < angle_slow,
            "higher speed → smaller CTE contribution"
        );
    }
    #[test]
    fn test_cross_track_statistics_zero_error() {
        let ref_path: Vec<[f64; 2]> = (0..5).map(|i| [i as f64 * 10.0, 0.0]).collect();
        let traj: Vec<[f64; 2]> = (0..5).map(|i| [i as f64 * 10.0, 0.0]).collect();
        let (mean, max, rms) = cross_track_statistics(&traj, &ref_path);
        assert!(mean < 1e-9);
        assert!(max < 1e-9);
        assert!(rms < 1e-9);
    }
    #[test]
    fn test_cross_track_statistics_constant_offset() {
        let ref_path: Vec<[f64; 2]> = (0..5).map(|i| [i as f64 * 10.0, 0.0]).collect();
        let traj: Vec<[f64; 2]> = (0..5).map(|i| [i as f64 * 10.0, 2.0]).collect();
        let (mean, max, _rms) = cross_track_statistics(&traj, &ref_path);
        assert!(
            (mean - 2.0).abs() < 0.01,
            "mean error should be ~2.0, got {mean}"
        );
        assert!(
            (max - 2.0).abs() < 0.01,
            "max error should be ~2.0, got {max}"
        );
    }
    #[test]
    fn test_min_curvature_lower_mean_than_midline_on_bend() {
        let n = 12;
        let left: Vec<[f64; 2]> = (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64 * PI;
                [5.0 * t.cos(), 5.0 * t.sin() + 4.0]
            })
            .collect();
        let right: Vec<[f64; 2]> = (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64 * PI;
                [3.0 * t.cos(), 3.0 * t.sin() + 4.0]
            })
            .collect();
        let line = MinCurvatureRacingLine::optimize(&left, &right, 20);
        for &k in &line.curvatures {
            assert!(k.is_finite());
            assert!(k >= 0.0);
        }
    }
}
