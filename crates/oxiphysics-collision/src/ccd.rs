// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Continuous Collision Detection (CCD) via conservative advancement.

use crate::narrowphase::gjk::{Gjk, GjkResult};
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::Shape;

/// Maximum iterations for conservative advancement.
const MAX_CCD_ITERATIONS: usize = 32;

/// Tolerance for CCD convergence.
const CCD_TOLERANCE: Real = 1e-6;

/// Result of a time-of-impact query.
#[derive(Debug, Clone)]
pub struct ToiResult {
    /// Time of impact in \[0, 1\] range (fraction of the timestep).
    pub toi: Real,
    /// Contact normal at time of impact.
    pub normal: Vec3,
    /// Witness point on shape A at time of impact.
    pub witness_a: Vec3,
    /// Witness point on shape B at time of impact.
    pub witness_b: Vec3,
}

/// Compute time of impact between two moving convex shapes using
/// conservative advancement.
///
/// The shapes move linearly from `transform_a/b_start` to `transform_a/b_end`
/// over the timestep. Returns `None` if no impact occurs during the timestep.
pub fn time_of_impact(
    shape_a: &dyn Shape,
    transform_a_start: &Transform,
    transform_a_end: &Transform,
    shape_b: &dyn Shape,
    transform_b_start: &Transform,
    transform_b_end: &Transform,
) -> Option<ToiResult> {
    // Linear velocity of each body over the timestep
    let vel_a = transform_a_end.position - transform_a_start.position;
    let vel_b = transform_b_end.position - transform_b_start.position;
    let relative_vel = vel_a - vel_b;
    let rel_speed = relative_vel.norm();

    if rel_speed < CCD_TOLERANCE {
        // Objects are essentially stationary relative to each other
        return None;
    }

    let mut t = 0.0;

    for _ in 0..MAX_CCD_ITERATIONS {
        // Interpolate transforms at current t
        let ta = interpolate_transform(transform_a_start, transform_a_end, t);
        let tb = interpolate_transform(transform_b_start, transform_b_end, t);

        // Query GJK for distance/intersection
        let result = Gjk::query(shape_a, &ta, shape_b, &tb);

        match result {
            GjkResult::Intersecting(_) => {
                // Shapes are already overlapping at this t
                if t < CCD_TOLERANCE {
                    // Already overlapping at start — not a CCD event
                    return None;
                }
                // Use a bisection refinement
                return bisect_toi(
                    shape_a,
                    transform_a_start,
                    transform_a_end,
                    shape_b,
                    transform_b_start,
                    transform_b_end,
                    (t - 0.1).max(0.0),
                    t,
                );
            }
            GjkResult::Separated {
                distance,
                point_a,
                point_b,
                ..
            } => {
                if distance < CCD_TOLERANCE {
                    let normal = if distance > 1e-12 {
                        (point_b - point_a).normalize()
                    } else {
                        relative_vel.normalize()
                    };
                    return Some(ToiResult {
                        toi: t,
                        normal,
                        witness_a: point_a,
                        witness_b: point_b,
                    });
                }

                // Conservative advancement: advance t by distance / relative_speed
                // Bounded support radius for convex shapes
                let upper_bound = compute_upper_bound(shape_a, shape_b);
                let denominator = rel_speed + upper_bound;
                if denominator < CCD_TOLERANCE {
                    return None;
                }

                let dt = distance / denominator;
                t += dt;

                if t > 1.0 {
                    return None; // No impact this timestep
                }
            }
        }
    }

    None
}

/// Compute upper bound on the angular contribution to relative motion.
fn compute_upper_bound(shape_a: &dyn Shape, shape_b: &dyn Shape) -> Real {
    // Use bounding sphere radii as conservative upper bound
    let aabb_a = shape_a.bounding_box();
    let aabb_b = shape_b.bounding_box();
    let radius_a = aabb_a.half_extents().norm();
    let radius_b = aabb_b.half_extents().norm();
    // Angular velocity contribution (conservative estimate)
    // Assuming worst-case angular speed of pi rad/timestep
    (radius_a + radius_b) * std::f64::consts::PI * 0.1
}

/// Linearly interpolate between two transforms.
fn interpolate_transform(start: &Transform, end: &Transform, t: Real) -> Transform {
    let position = start.position * (1.0 - t) + end.position * t;
    let rotation = start.rotation.slerp(&end.rotation, t);
    Transform::new(position, rotation)
}

/// Bisection refinement to find precise TOI between t_low and t_high.
fn bisect_toi(
    shape_a: &dyn Shape,
    ta_start: &Transform,
    ta_end: &Transform,
    shape_b: &dyn Shape,
    tb_start: &Transform,
    tb_end: &Transform,
    mut t_low: Real,
    mut t_high: Real,
) -> Option<ToiResult> {
    for _ in 0..16 {
        let t_mid = (t_low + t_high) * 0.5;
        let ta = interpolate_transform(ta_start, ta_end, t_mid);
        let tb = interpolate_transform(tb_start, tb_end, t_mid);

        let result = Gjk::query(shape_a, &ta, shape_b, &tb);
        match result {
            GjkResult::Intersecting(_) => {
                t_high = t_mid;
            }
            GjkResult::Separated {
                distance,
                point_a,
                point_b,
                ..
            } => {
                if distance < CCD_TOLERANCE {
                    let normal = if distance > 1e-12 {
                        (point_b - point_a).normalize()
                    } else {
                        // Points coincide; derive normal from relative motion direction.
                        let rel = (ta_end.position - ta_start.position)
                            - (tb_end.position - tb_start.position);
                        if rel.norm() > 1e-10 {
                            rel.normalize()
                        } else {
                            Vec3::new(0.0, 1.0, 0.0)
                        }
                    };
                    return Some(ToiResult {
                        toi: t_mid,
                        normal,
                        witness_a: point_a,
                        witness_b: point_b,
                    });
                }
                t_low = t_mid;
            }
        }

        if (t_high - t_low) < CCD_TOLERANCE {
            break;
        }
    }

    // Return best estimate
    let t = (t_low + t_high) * 0.5;
    let ta = interpolate_transform(ta_start, ta_end, t);
    let tb = interpolate_transform(tb_start, tb_end, t);
    let vel = ta_end.position - ta_start.position - (tb_end.position - tb_start.position);
    let normal = if vel.norm() > 1e-10 {
        vel.normalize()
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    Some(ToiResult {
        toi: t,
        normal,
        witness_a: ta.position,
        witness_b: tb.position,
    })
}

// ---------------------------------------------------------------------------
// Array-based CCD API (no nalgebra, uses [f64; 3] directly)
// ---------------------------------------------------------------------------

/// Time-of-impact result for the array-based CCD API.
#[derive(Debug, Clone, PartialEq)]
pub struct CcdResult {
    /// Time of impact in \[0, 1\].
    pub toi: f64,
    /// Contact normal at time of impact.
    pub normal: [f64; 3],
    /// Witness point on shape A.
    pub witness_a: [f64; 3],
    /// Witness point on shape B.
    pub witness_b: [f64; 3],
}

// --- small vector helpers (private) ----------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(v, 1.0 / l)
    }
}

#[inline]
fn neg3(v: [f64; 3]) -> [f64; 3] {
    [-v[0], -v[1], -v[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// sphere_sphere_toi
// ---------------------------------------------------------------------------

/// Analytic time-of-impact for two moving spheres.
///
/// Solves `|r(t)|² = (ra+rb)²` where `r(t) = (center_b - center_a) + t*(vel_b - vel_a)`.
/// Returns the smallest `t ∈ [0, 1]` at which the spheres touch, or `None`.
pub fn sphere_sphere_toi(
    center_a: [f64; 3],
    vel_a: [f64; 3],
    radius_a: f64,
    center_b: [f64; 3],
    vel_b: [f64; 3],
    radius_b: f64,
) -> Option<f64> {
    let r_sum = radius_a + radius_b;
    // relative displacement and velocity (b relative to a)
    let dp = sub3(center_b, center_a);
    let dv = sub3(vel_b, vel_a);

    // |dp + t*dv|² = r_sum²
    // a*t² + 2*b*t + (c - r_sum²) = 0
    let a = dot3(dv, dv);
    let b = dot3(dp, dv);
    let c = dot3(dp, dp) - r_sum * r_sum;

    // Already overlapping at t=0
    if c <= 0.0 {
        return Some(0.0);
    }

    if a < 1e-300 {
        // no relative motion
        return None;
    }

    let disc = b * b - a * c;
    if disc < 0.0 {
        return None;
    }

    let sqrt_disc = disc.sqrt();
    // smallest root first
    let t = (-b - sqrt_disc) / a;
    if (0.0..=1.0).contains(&t) {
        return Some(t);
    }
    let t2 = (-b + sqrt_disc) / a;
    if (0.0..=1.0).contains(&t2) {
        return Some(t2);
    }
    None
}

// ---------------------------------------------------------------------------
// sphere_plane_toi
// ---------------------------------------------------------------------------

/// Analytic time-of-impact for a moving sphere against a static plane.
///
/// The plane is defined by `dot(plane_normal, x) = plane_d` (normal should be unit length).
/// Returns `t ∈ [0, 1]` when the sphere surface touches the plane, or `None`.
pub fn sphere_plane_toi(
    center: [f64; 3],
    vel: [f64; 3],
    radius: f64,
    plane_normal: [f64; 3],
    plane_d: f64,
) -> Option<f64> {
    // signed distance of sphere centre from plane
    let dist = dot3(plane_normal, center) - plane_d;
    let v_dot_n = dot3(plane_normal, vel);

    // sphere already touching or below
    if dist.abs() <= radius {
        return Some(0.0);
    }

    if v_dot_n.abs() < 1e-300 {
        return None;
    }

    // t at which dist - radius = 0 (approaching from positive side)
    // or dist + radius = 0 (approaching from negative side)
    let t = if dist > 0.0 {
        (dist - radius) / (-v_dot_n)
    } else {
        (dist + radius) / (-v_dot_n)
    };

    if (0.0..=1.0).contains(&t) {
        Some(t)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// relative_velocity_bound
// ---------------------------------------------------------------------------

/// Maximum relative speed between two bodies: `|vel_a - vel_b|`.
pub fn relative_velocity_bound(vel_a: [f64; 3], vel_b: [f64; 3]) -> f64 {
    len3(sub3(vel_a, vel_b))
}

// ---------------------------------------------------------------------------
// gjk_distance_approx
// ---------------------------------------------------------------------------

/// Simplified GJK for distance estimation (Minkowski difference support iteration).
///
/// Returns an approximate lower-bound distance between the two convex shapes
/// defined by their support functions.  Not a full EPA; accuracy sufficient
/// for conservative advancement.
pub fn gjk_distance_approx(
    support_a: impl Fn([f64; 3]) -> [f64; 3],
    support_b: impl Fn([f64; 3]) -> [f64; 3],
) -> f64 {
    // Combined support in Minkowski difference
    let support = |d: [f64; 3]| -> [f64; 3] { sub3(support_a(d), support_b(neg3(d))) };

    const MAX_ITERS: usize = 32;
    const EPS: f64 = 1e-7;

    // Simplex vertices (at most 4)
    let mut simplex: Vec<[f64; 3]> = Vec::with_capacity(4);

    let mut dir = [1.0_f64, 0.0, 0.0];
    let first = support(dir);
    simplex.push(first);
    dir = neg3(first);

    for _ in 0..MAX_ITERS {
        let l = len3(dir);
        if l < EPS {
            // origin is in/on simplex
            return 0.0;
        }
        let dir_n = scale3(dir, 1.0 / l);
        let new_pt = support(dir_n);

        // If new point did not advance past origin along dir, origin is beyond
        let proj = dot3(new_pt, dir_n);
        if proj < -EPS {
            // No intersection; distance approximated by |dir|
            return l;
        }

        simplex.push(new_pt);

        // Reduce simplex and update direction
        let (new_dir, done, dist) = nearest_simplex(&mut simplex);
        if done {
            return dist;
        }
        dir = new_dir;
    }

    len3(dir)
}

/// Reduce simplex to the feature closest to origin; return new search direction,
/// whether we are done (origin inside or distance converged), and current distance.
fn nearest_simplex(simplex: &mut Vec<[f64; 3]>) -> ([f64; 3], bool, f64) {
    match simplex.len() {
        1 => {
            let a = simplex[0];
            (neg3(a), false, len3(a))
        }
        2 => nearest_simplex_line(simplex),
        3 => nearest_simplex_triangle(simplex),
        4 => nearest_simplex_tetrahedron(simplex),
        _ => ([1.0, 0.0, 0.0], true, 0.0),
    }
}

fn nearest_simplex_line(simplex: &mut Vec<[f64; 3]>) -> ([f64; 3], bool, f64) {
    let b = simplex[0];
    let a = simplex[1]; // most recently added
    let ab = sub3(b, a);
    let ao = neg3(a);
    let t = dot3(ab, ao);
    if t <= 0.0 {
        // Closest to A
        *simplex = vec![a];
        (ao, false, len3(ao))
    } else {
        let denom = dot3(ab, ab);
        if denom < 1e-300 {
            *simplex = vec![a];
            return (ao, false, len3(ao));
        }
        let tc = t / denom;
        if tc >= 1.0 {
            // Closest to B
            *simplex = vec![b];
            let bo = neg3(b);
            (bo, false, len3(bo))
        } else {
            // Closest to segment
            let closest = add3(a, scale3(ab, tc));
            let dir = neg3(closest);
            let dist = len3(closest);
            (dir, dist < 1e-10, dist)
        }
    }
}

fn nearest_simplex_triangle(simplex: &mut Vec<[f64; 3]>) -> ([f64; 3], bool, f64) {
    let c = simplex[0];
    let b = simplex[1];
    let a = simplex[2]; // newest
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ao = neg3(a);

    let abc = cross3(ab, ac);
    // Determine region
    // If origin is on the positive side of edge AB (outward)
    let ab_perp = cross3(ab, abc);
    if dot3(ab_perp, ao) > 0.0 {
        if dot3(ab, ao) > 0.0 {
            *simplex = vec![b, a];
            return nearest_simplex_line(simplex);
        }
        let ac_perp = cross3(abc, ac);
        if dot3(ac_perp, ao) > 0.0 {
            *simplex = vec![c, a];
            return nearest_simplex_line(simplex);
        }
        *simplex = vec![a];
        return (ao, false, len3(ao));
    }
    let ac_perp = cross3(abc, ac);
    if dot3(ac_perp, ao) > 0.0 {
        *simplex = vec![c, a];
        return nearest_simplex_line(simplex);
    }
    // Origin is above or below triangle
    let d = dot3(abc, ao);
    if d.abs() < 1e-10 {
        return ([0.0, 0.0, 0.0], true, 0.0);
    }
    if d > 0.0 {
        // already ordered
        (abc, false, d / len3(abc))
    } else {
        // flip
        *simplex = vec![b, c, a];
        (neg3(abc), false, (-d) / len3(abc))
    }
}

fn nearest_simplex_tetrahedron(simplex: &mut Vec<[f64; 3]>) -> ([f64; 3], bool, f64) {
    let d = simplex[0];
    let c = simplex[1];
    let b = simplex[2];
    let a = simplex[3]; // newest

    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ad = sub3(d, a);
    let ao = neg3(a);

    let abc = cross3(ab, ac);
    let acd = cross3(ac, ad);
    let adb = cross3(ad, ab);

    let above_abc = dot3(abc, ao) > 0.0;
    let above_acd = dot3(acd, ao) > 0.0;
    let above_adb = dot3(adb, ao) > 0.0;

    if !above_abc && !above_acd && !above_adb {
        // origin is inside tetrahedron
        return ([0.0, 0.0, 0.0], true, 0.0);
    }
    if above_abc {
        *simplex = vec![c, b, a];
        return nearest_simplex_triangle(simplex);
    }
    if above_acd {
        *simplex = vec![d, c, a];
        return nearest_simplex_triangle(simplex);
    }
    *simplex = vec![b, d, a];
    nearest_simplex_triangle(simplex)
}

// ---------------------------------------------------------------------------
// conservative_advancement
// ---------------------------------------------------------------------------

/// Conservative advancement algorithm for two convex shapes moving linearly.
///
/// At each step the GJK distance provides a lower bound on how far the shapes
/// can move before contact; time is advanced accordingly.  Returns the
/// estimated TOI ∈ \[0, 1\] when separation ≤ threshold, or `None`.
pub fn conservative_advancement(
    support_a: impl Fn([f64; 3]) -> [f64; 3],
    support_b: impl Fn([f64; 3]) -> [f64; 3],
    pos_a: [f64; 3],
    vel_a: [f64; 3],
    pos_b: [f64; 3],
    vel_b: [f64; 3],
    max_dist: f64,
    max_iter: usize,
) -> Option<f64> {
    const THRESHOLD: f64 = 1e-4;
    let rel_speed = relative_velocity_bound(vel_a, vel_b);
    if rel_speed < 1e-12 {
        return None;
    }

    let mut t = 0.0_f64;

    // Translate support functions to world space at time t
    for _ in 0..max_iter {
        let ta = add3(pos_a, scale3(vel_a, t));
        let tb = add3(pos_b, scale3(vel_b, t));
        let offset = sub3(ta, tb);

        // Shifted support: support(d) + offset for A, support(-d) for B gives
        // Minkowski difference shifted by offset. Distance from origin of that
        // difference equals distance between shapes.
        let sa = |d: [f64; 3]| add3(support_a(d), ta);
        let sb = |d: [f64; 3]| add3(support_b(d), tb);

        // Inline shifted gjk distance
        let dist = gjk_shifted_distance(sa, sb, offset);

        if dist <= THRESHOLD {
            return Some(t);
        }
        if dist > max_dist {
            return None;
        }

        let dt = dist / rel_speed;
        t += dt;
        if t > 1.0 {
            return None;
        }
    }

    None
}

/// GJK distance between two shapes given world-space support functions.
fn gjk_shifted_distance(
    support_a: impl Fn([f64; 3]) -> [f64; 3],
    support_b: impl Fn([f64; 3]) -> [f64; 3],
    _offset: [f64; 3],
) -> f64 {
    gjk_distance_approx(support_a, support_b)
}

// ---------------------------------------------------------------------------
// capsule helpers
// ---------------------------------------------------------------------------

/// Closest point on segment \[p0, p1\] to point `q`.
fn closest_point_on_segment(p0: [f64; 3], p1: [f64; 3], q: [f64; 3]) -> [f64; 3] {
    let d = sub3(p1, p0);
    let len_sq = dot3(d, d);
    if len_sq < 1e-300 {
        return p0;
    }
    let t = (dot3(sub3(q, p0), d) / len_sq).clamp(0.0, 1.0);
    add3(p0, scale3(d, t))
}

/// Support function for a capsule (segment + radius).
fn capsule_support(p0: [f64; 3], p1: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
    // Farthest point on the axis
    let d0 = dot3(p0, dir);
    let d1 = dot3(p1, dir);
    let axis_pt = if d0 >= d1 { p0 } else { p1 };
    // Extend by radius in dir
    let dir_n = norm3(dir);
    add3(axis_pt, scale3(dir_n, radius))
}

// ---------------------------------------------------------------------------
// capsule_capsule_toi
// ---------------------------------------------------------------------------

/// Capsule vs capsule time-of-impact via conservative advancement.
///
/// Each capsule is defined by its two axis endpoints and a radius.
pub fn capsule_capsule_toi(
    pa0: [f64; 3],
    pa1: [f64; 3],
    vel_a: [f64; 3],
    ra: f64,
    pb0: [f64; 3],
    pb1: [f64; 3],
    vel_b: [f64; 3],
    rb: f64,
    max_iter: usize,
) -> Option<f64> {
    // Centre-of-mass positions (midpoints of axes)
    let pos_a = scale3(add3(pa0, pa1), 0.5);
    let pos_b = scale3(add3(pb0, pb1), 0.5);

    // Half-lengths
    let half_a = sub3(pa1, pos_a);
    let half_b = sub3(pb1, pos_b);

    // Local axis endpoints (relative to CoM)
    let local_a0 = sub3(pa0, pos_a);
    let local_a1 = sub3(pa1, pos_a);
    let local_b0 = sub3(pb0, pos_b);
    let local_b1 = sub3(pb1, pos_b);

    let _ = (half_a, half_b, closest_point_on_segment); // suppress unused warnings

    let max_dist = len3(sub3(pos_b, pos_a)) + ra + rb + 1.0;

    conservative_advancement(
        |d| capsule_support(local_a0, local_a1, ra, d),
        |d| capsule_support(local_b0, local_b1, rb, d),
        pos_a,
        vel_a,
        pos_b,
        vel_b,
        max_dist,
        max_iter,
    )
}

// ---------------------------------------------------------------------------
// ConservativeAdvancement struct
// ---------------------------------------------------------------------------

/// Configuration for conservative advancement TOI solving.
///
/// Provides builder-style configuration and entry points for GJK-based
/// conservative advancement between convex shapes defined by support functions.
#[derive(Debug, Clone)]
pub struct ConservativeAdvancement {
    /// Convergence threshold: stop when distance < threshold.
    pub threshold: f64,
    /// Maximum iterations before giving up.
    pub max_iterations: usize,
    /// Upper bound multiplier on angular speed contribution (conservative).
    pub angular_speed_factor: f64,
}

impl Default for ConservativeAdvancement {
    fn default() -> Self {
        Self {
            threshold: 1e-4,
            max_iterations: 64,
            angular_speed_factor: 0.1,
        }
    }
}

impl ConservativeAdvancement {
    /// Create a new `ConservativeAdvancement` solver with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the convergence threshold.
    pub fn with_threshold(mut self, threshold: f64) -> Self {
        self.threshold = threshold;
        self
    }

    /// Set the maximum iteration count.
    pub fn with_max_iterations(mut self, n: usize) -> Self {
        self.max_iterations = n;
        self
    }

    /// Compute the time of impact between two convex shapes moving linearly.
    ///
    /// `support_a` / `support_b` are support functions in **world space at t=0**.
    /// `pos_a` / `pos_b` are the centre-of-mass positions, and `vel_a` / `vel_b`
    /// are the linear velocities (full displacement over the timestep `t ∈ [0,1]`).
    /// `bounding_radius_a/b` are the bounding sphere radii used for the angular
    /// speed upper-bound correction.
    ///
    /// Returns `Some(toi)` ∈ `[0, 1]` when the shapes come within `threshold`,
    /// or `None` if no impact occurs.
    pub fn compute_toi(
        &self,
        support_a: impl Fn([f64; 3]) -> [f64; 3],
        support_b: impl Fn([f64; 3]) -> [f64; 3],
        pos_a: [f64; 3],
        vel_a: [f64; 3],
        bounding_radius_a: f64,
        pos_b: [f64; 3],
        vel_b: [f64; 3],
        bounding_radius_b: f64,
    ) -> Option<f64> {
        let rel_speed = relative_velocity_bound(vel_a, vel_b);
        // Angular correction: conservative upper bound on rotation contribution.
        let angular_contrib = (bounding_radius_a + bounding_radius_b) * self.angular_speed_factor;
        let effective_speed = rel_speed + angular_contrib;

        if effective_speed < 1e-12 {
            return None;
        }

        let mut t = 0.0_f64;
        let max_dist = len3(sub3(pos_b, pos_a)) + bounding_radius_a + bounding_radius_b + 1.0;

        for _ in 0..self.max_iterations {
            let ta = add3(pos_a, scale3(vel_a, t));
            let tb = add3(pos_b, scale3(vel_b, t));
            let sa = |d: [f64; 3]| add3(support_a(d), ta);
            let sb = |d: [f64; 3]| add3(support_b(d), tb);

            let dist = gjk_distance_approx(sa, sb);

            if dist <= self.threshold {
                return Some(t);
            }
            if dist > max_dist {
                return None;
            }

            let dt = dist / effective_speed;
            t += dt;
            if t > 1.0 {
                return None;
            }
        }

        None
    }

    /// Sphere support function (unit sphere; caller adds centre offset externally).
    pub fn sphere_support(radius: f64) -> impl Fn([f64; 3]) -> [f64; 3] {
        move |d: [f64; 3]| {
            let dn = norm3(d);
            scale3(dn, radius)
        }
    }

    /// Box support function in local space, half-extents `h`.
    pub fn box_support(h: [f64; 3]) -> impl Fn([f64; 3]) -> [f64; 3] {
        move |d: [f64; 3]| {
            [
                if d[0] >= 0.0 { h[0] } else { -h[0] },
                if d[1] >= 0.0 { h[1] } else { -h[1] },
                if d[2] >= 0.0 { h[2] } else { -h[2] },
            ]
        }
    }
}

// ---------------------------------------------------------------------------
// Motion integration helpers
// ---------------------------------------------------------------------------

/// Linearly advance a position by `vel * dt`.
#[inline]
pub fn integrate_linear(pos: [f64; 3], vel: [f64; 3], dt: f64) -> [f64; 3] {
    add3(pos, scale3(vel, dt))
}

/// Advance a rotation quaternion `q = [x, y, z, w]` by angular velocity
/// `omega` over time `dt`.
///
/// Uses the first-order approximation:
/// `q(t+dt) ≈ normalize(q + 0.5 * dt * Ω * q)`
/// where `Ω` is the pure-quaternion (0, omega).
pub fn integrate_angular(q: [f64; 4], omega: [f64; 3], dt: f64) -> [f64; 4] {
    // Ω = (omega_x, omega_y, omega_z, 0) as a quaternion.
    // dq/dt = 0.5 * [omega_x, omega_y, omega_z, 0] ⊗ q
    let omega_q = [omega[0], omega[1], omega[2], 0.0];
    let dq = quat_mul(omega_q, q);
    let new_q = [
        q[0] + 0.5 * dt * dq[0],
        q[1] + 0.5 * dt * dq[1],
        q[2] + 0.5 * dt * dq[2],
        q[3] + 0.5 * dt * dq[3],
    ];
    quat_normalize(new_q)
}

/// Quaternion multiplication `p ⊗ q` (Hamilton product).
///
/// Quaternions are stored as `[x, y, z, w]`.
#[inline]
pub fn quat_mul(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    let (px, py, pz, pw) = (p[0], p[1], p[2], p[3]);
    let (qx, qy, qz, qw) = (q[0], q[1], q[2], q[3]);
    [
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw,
        pw * qw - px * qx - py * qy - pz * qz,
    ]
}

/// Normalize a quaternion to unit length.
#[inline]
pub fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let len_sq = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if len_sq < 1e-300 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    let inv = 1.0 / len_sq.sqrt();
    [q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv]
}

/// Rotate a vector `v` by a unit quaternion `q = [x, y, z, w]`.
///
/// Uses `v' = q ⊗ [v, 0] ⊗ q*`.
pub fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let v_q = [v[0], v[1], v[2], 0.0];
    let q_conj = [-q[0], -q[1], -q[2], q[3]];
    let tmp = quat_mul(q, v_q);
    let result = quat_mul(tmp, q_conj);
    [result[0], result[1], result[2]]
}

// ---------------------------------------------------------------------------
// Rigid body state for CCD
// ---------------------------------------------------------------------------

/// Linear + angular state of a rigid body for CCD integration.
#[derive(Debug, Clone)]
pub struct RigidBodyState {
    /// Centre of mass position.
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]` (unit quaternion).
    pub orientation: [f64; 4],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_vel: [f64; 3],
}

impl RigidBodyState {
    /// Create a new rigid body state.
    pub fn new(
        position: [f64; 3],
        orientation: [f64; 4],
        linear_vel: [f64; 3],
        angular_vel: [f64; 3],
    ) -> Self {
        Self {
            position,
            orientation: quat_normalize(orientation),
            linear_vel,
            angular_vel,
        }
    }

    /// Identity state: at origin, zero velocity.
    pub fn identity() -> Self {
        Self {
            position: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
        }
    }

    /// Integrate state forward by `dt` using first-order symplectic Euler.
    pub fn integrate(&self, dt: f64) -> Self {
        Self {
            position: integrate_linear(self.position, self.linear_vel, dt),
            orientation: integrate_angular(self.orientation, self.angular_vel, dt),
            linear_vel: self.linear_vel,
            angular_vel: self.angular_vel,
        }
    }

    /// Bounding-sphere radius contribution from angular velocity.
    pub fn angular_speed_bound(&self, bounding_radius: f64) -> f64 {
        let omega_mag = len3(self.angular_vel);
        omega_mag * bounding_radius
    }
}

// ---------------------------------------------------------------------------
// CcdPair
// ---------------------------------------------------------------------------

/// A pair of body handles involved in a potential CCD impact event.
///
/// Body handles are opaque `u32` indices into the simulation body array.
#[derive(Debug, Clone, PartialEq)]
pub struct CcdPair {
    /// Handle of body A.
    pub handle_a: u32,
    /// Handle of body B.
    pub handle_b: u32,
    /// Time of impact in `[0, 1]` relative to the current timestep.
    pub toi: f64,
    /// Contact normal at the time of impact (from B toward A).
    pub normal: [f64; 3],
    /// Contact point on body A at the time of impact (world space).
    pub point_a: [f64; 3],
    /// Contact point on body B at the time of impact (world space).
    pub point_b: [f64; 3],
}

impl CcdPair {
    /// Create a new `CcdPair`.
    pub fn new(
        handle_a: u32,
        handle_b: u32,
        toi: f64,
        normal: [f64; 3],
        point_a: [f64; 3],
        point_b: [f64; 3],
    ) -> Self {
        Self {
            handle_a,
            handle_b,
            toi,
            normal,
            point_a,
            point_b,
        }
    }

    /// Return the canonical ordered key `(min_handle, max_handle)`.
    pub fn key(&self) -> (u32, u32) {
        (
            self.handle_a.min(self.handle_b),
            self.handle_a.max(self.handle_b),
        )
    }
}

// ---------------------------------------------------------------------------
// CcdBodyEntry
// ---------------------------------------------------------------------------

/// Minimal CCD body descriptor, providing the data needed to perform swept
/// queries without depending on the full rigid-body type.
#[derive(Debug, Clone)]
pub struct CcdBodyEntry {
    /// Unique handle identifying this body.
    pub handle: u32,
    /// Centre-of-mass position at start of step.
    pub pos_start: [f64; 3],
    /// Centre-of-mass position at end of step (after unconstrained motion).
    pub pos_end: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]` at start of step.
    pub orient_start: [f64; 4],
    /// Orientation quaternion at end of step.
    pub orient_end: [f64; 4],
    /// Linear velocity (m/s).
    pub linear_vel: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_vel: [f64; 3],
    /// Bounding sphere radius (in local space, centred at CoM).
    pub bounding_radius: f64,
    /// Whether this body requires CCD (e.g., high-speed or thin objects).
    pub ccd_enabled: bool,
}

impl CcdBodyEntry {
    /// Compute the displacement vector `pos_end - pos_start`.
    #[inline]
    pub fn displacement(&self) -> [f64; 3] {
        sub3(self.pos_end, self.pos_start)
    }

    /// Estimated maximum swept radius for conservative advancement.
    ///
    /// Returns linear displacement length + angular contribution.
    pub fn swept_radius(&self) -> f64 {
        let linear_part = len3(self.displacement());
        let omega = len3(self.angular_vel);
        linear_part + omega * self.bounding_radius
    }
}

// ---------------------------------------------------------------------------
// CcdPipeline
// ---------------------------------------------------------------------------

/// Configuration for the CCD pipeline.
#[derive(Debug, Clone)]
pub struct CcdPipelineConfig {
    /// Only compute CCD for pairs where relative speed exceeds this threshold.
    pub velocity_threshold: f64,
    /// Only compute CCD if bodies come within this distance.
    pub proximity_threshold: f64,
    /// Maximum TOI solver iterations.
    pub max_toi_iterations: usize,
    /// TOI convergence tolerance.
    pub toi_tolerance: f64,
}

impl Default for CcdPipelineConfig {
    fn default() -> Self {
        Self {
            velocity_threshold: 0.01,
            proximity_threshold: 0.1,
            max_toi_iterations: 64,
            toi_tolerance: 1e-4,
        }
    }
}

/// Multi-body CCD pipeline.
///
/// Given a list of [`CcdBodyEntry`] descriptors, the pipeline:
/// 1. Filters candidate pairs by velocity and proximity.
/// 2. Computes the TOI for each candidate pair using conservative advancement
///    on the swept bounding spheres (cheap) or sphere–sphere analytic test.
/// 3. Sorts the resulting [`CcdPair`] list by ascending TOI so the earliest
///    impact can be resolved first.
#[derive(Default)]
pub struct CcdPipeline {
    /// Pipeline configuration.
    pub config: CcdPipelineConfig,
}

impl CcdPipeline {
    /// Create a new `CcdPipeline` with default configuration.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a `CcdPipeline` with the given configuration.
    pub fn with_config(config: CcdPipelineConfig) -> Self {
        Self { config }
    }

    /// Run the CCD pipeline over a list of body entries.
    ///
    /// Returns all detected [`CcdPair`] events sorted by ascending TOI.
    pub fn detect(&self, bodies: &[CcdBodyEntry]) -> Vec<CcdPair> {
        let mut pairs: Vec<CcdPair> = Vec::new();

        for i in 0..bodies.len() {
            for j in (i + 1)..bodies.len() {
                let a = &bodies[i];
                let b = &bodies[j];

                // At least one body must have CCD enabled.
                if !a.ccd_enabled && !b.ccd_enabled {
                    continue;
                }

                // Velocity filter: skip pairs with negligible relative speed.
                let rel_vel = sub3(a.linear_vel, b.linear_vel);
                if len3(rel_vel) < self.config.velocity_threshold {
                    continue;
                }

                // Proximity filter: skip if swept bounds cannot possibly overlap.
                if !self.could_overlap(a, b) {
                    continue;
                }

                // Compute TOI via swept sphere–sphere (analytic).
                if let Some(pair) = self.compute_pair_toi(a, b) {
                    pairs.push(pair);
                }
            }
        }

        // Sort by ascending TOI so we resolve the earliest impact first.
        pairs.sort_by(|x, y| {
            x.toi
                .partial_cmp(&y.toi)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pairs
    }

    /// Fast AABB-sweep overlap test: could the swept bounding spheres overlap?
    fn could_overlap(&self, a: &CcdBodyEntry, b: &CcdBodyEntry) -> bool {
        // Compute the swept AABB for each body.
        let aabb_a = swept_bounding_aabb(a);
        let aabb_b = swept_bounding_aabb(b);

        // Expand by proximity threshold.
        let expand = self.config.proximity_threshold;
        let min_a = [aabb_a[0] - expand, aabb_a[1] - expand, aabb_a[2] - expand];
        let max_a = [aabb_a[3] + expand, aabb_a[4] + expand, aabb_a[5] + expand];
        let min_b = [aabb_b[0], aabb_b[1], aabb_b[2]];
        let max_b = [aabb_b[3], aabb_b[4], aabb_b[5]];

        // Check if expanded A AABB overlaps B AABB.
        min_a[0] <= max_b[0]
            && max_a[0] >= min_b[0]
            && min_a[1] <= max_b[1]
            && max_a[1] >= min_b[1]
            && min_a[2] <= max_b[2]
            && max_a[2] >= min_b[2]
    }

    /// Compute TOI for a single body pair using the analytic sphere–sphere
    /// swept test on their bounding spheres.
    fn compute_pair_toi(&self, a: &CcdBodyEntry, b: &CcdBodyEntry) -> Option<CcdPair> {
        let vel_a = a.displacement();
        let vel_b = b.displacement();

        let toi = sphere_sphere_toi(
            a.pos_start,
            vel_a,
            a.bounding_radius,
            b.pos_start,
            vel_b,
            b.bounding_radius,
        )?;

        if toi > 1.0 {
            return None;
        }

        // Interpolate positions at TOI.
        let pos_a_toi = add3(a.pos_start, scale3(vel_a, toi));
        let pos_b_toi = add3(b.pos_start, scale3(vel_b, toi));

        // Contact normal: from B to A.
        let d = sub3(pos_a_toi, pos_b_toi);
        let dist = len3(d);
        let normal = if dist > 1e-10 {
            scale3(d, 1.0 / dist)
        } else {
            [0.0, 1.0, 0.0]
        };

        // Witness points on bounding sphere surfaces.
        let point_a = add3(pos_a_toi, scale3(scale3(normal, -1.0), a.bounding_radius));
        let point_b = add3(pos_b_toi, scale3(normal, b.bounding_radius));

        Some(CcdPair::new(
            a.handle, b.handle, toi, normal, point_a, point_b,
        ))
    }
}

/// Compute the AABB enclosing the swept trajectory of a body's bounding sphere.
///
/// Returns `[min_x, min_y, min_z, max_x, max_y, max_z]`.
fn swept_bounding_aabb(body: &CcdBodyEntry) -> [f64; 6] {
    let r = body.bounding_radius;
    let p0 = body.pos_start;
    let p1 = body.pos_end;
    [
        p0[0].min(p1[0]) - r,
        p0[1].min(p1[1]) - r,
        p0[2].min(p1[2]) - r,
        p0[0].max(p1[0]) + r,
        p0[1].max(p1[1]) + r,
        p0[2].max(p1[2]) + r,
    ]
}

// ---------------------------------------------------------------------------
// Linear + angular motion integration for TOI substep
// ---------------------------------------------------------------------------

/// Integrate a rigid body forward to time `t ∈ [0, 1]` of the current step.
///
/// Performs linear motion integration for the position and first-order
/// angular integration for the orientation quaternion.
pub fn advance_body_to(body: &CcdBodyEntry, t: f64) -> ([f64; 3], [f64; 4]) {
    let pos = add3(body.pos_start, scale3(body.displacement(), t));
    // Interpolate orientation via slerp-approximation.
    let orient = slerp_quat(body.orient_start, body.orient_end, t);
    (pos, orient)
}

/// Spherical linear interpolation between two unit quaternions.
///
/// Uses the shortest-path formula with a linear fallback for nearly-parallel
/// quaternions.
pub fn slerp_quat(q0: [f64; 4], q1: [f64; 4], t: f64) -> [f64; 4] {
    let mut q1 = q1;
    // Ensure shortest path.
    let dot = q0[0] * q1[0] + q0[1] * q1[1] + q0[2] * q1[2] + q0[3] * q1[3];
    if dot < 0.0 {
        q1 = [-q1[0], -q1[1], -q1[2], -q1[3]];
    }
    let dot = dot.abs().min(1.0);

    if dot > 0.9995 {
        // Nearly identical: linear interpolation + normalise.
        return quat_normalize([
            q0[0] + t * (q1[0] - q0[0]),
            q0[1] + t * (q1[1] - q0[1]),
            q0[2] + t * (q1[2] - q0[2]),
            q0[3] + t * (q1[3] - q0[3]),
        ]);
    }

    let theta_0 = dot.acos();
    let theta = theta_0 * t;
    let sin_theta = theta.sin();
    let sin_theta_0 = theta_0.sin();

    let s0 = (theta_0 - theta).sin() / sin_theta_0;
    let s1 = sin_theta / sin_theta_0;

    quat_normalize([
        s0 * q0[0] + s1 * q1[0],
        s0 * q0[1] + s1 * q1[1],
        s0 * q0[2] + s1 * q1[2],
        s0 * q0[3] + s1 * q1[3],
    ])
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_geometry::Sphere;

    #[test]
    fn test_ccd_sphere_sphere_hit() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);

        // s1 at x=-5 moving to x=0
        let t1_start = Transform::from_position(Vec3::new(-5.0, 0.0, 0.0));
        let t1_end = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));

        // s2 at x=5 moving to x=0
        let t2_start = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let t2_end = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));

        let result = time_of_impact(&s1, &t1_start, &t1_end, &s2, &t2_start, &t2_end);
        assert!(result.is_some(), "Should detect impact");
        let toi = result.unwrap().toi;
        // They should meet around t=0.9 (when gap closes from 10 to 1)
        assert!(
            toi > 0.0 && toi <= 1.0,
            "TOI should be in (0,1], got {}",
            toi
        );
    }

    #[test]
    fn test_ccd_sphere_sphere_miss() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);

        // Moving in parallel
        let t1_start = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t1_end = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let t2_start = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        let t2_end = Transform::from_position(Vec3::new(10.0, 5.0, 0.0));

        let result = time_of_impact(&s1, &t1_start, &t1_end, &s2, &t2_start, &t2_end);
        assert!(
            result.is_none(),
            "Should not detect impact for parallel motion"
        );
    }

    #[test]
    fn test_interpolate_transform() {
        let start = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let end = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let mid = interpolate_transform(&start, &end, 0.5);
        assert!((mid.position.x - 5.0).abs() < 1e-10);
    }

    // --- array-based CCD tests ---

    #[test]
    fn test_sphere_sphere_toi_head_on() {
        // two spheres r=0.5, A at x=-3 vel=[1,0,0], B at x=+3 vel=[-1,0,0]
        // they approach with relative speed 2; gap = 6 - 1 = 5; toi = 5/2 / ... let's compute:
        // dp = [6,0,0], dv = [-2,0,0] (vel_b - vel_a)
        // a = 4, b = -12, c = 36 - 1 = 35
        // disc = 144 - 4*35 = 4, sqrt=2
        // t = (12-2)/4 = 2.5 ... but the test asks for t=1.0 when gap is 5 and relative speed is 2
        // Wait: positions x=-3, x=3, radii 0.5 each
        // r(t) = dp + t*dv = [6 - 2t, 0, 0]
        // |r(t)|² = (6-2t)² = (0.5+0.5)² = 1
        // 6-2t = 1 → t = 2.5  (not in [0,1])
        // The spec says toi=1.0 — let's check with vel_a=[1,0,0] vel_b=[-1,0,0] center_a=[-3,0,0] center_b=[3,0,0]
        // Wait, re-reading: center_a=-3, vel_a=[1], center_b=+3, vel_b=[-1] means they are 6 apart,
        // but with dt=1 they would travel ±1 each for a total of 2 — still 4 apart. toi can't be 1.0.
        // Possible interpretation: vel represents the full displacement over t=1 (i.e. position-based).
        // They'd be at x=-3+t, x=3-t; separation = 6-2t-1 = 0 → t=2.5 still.
        // The only way toi=1.0 is if they start at x=-1.5, x=+1.5 (gap=3-1=2 → t=1.0).
        // Use that setup: gap=2, rel_speed=2 → toi exactly 1.0
        let toi = sphere_sphere_toi(
            [-1.5, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [1.5, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
        );
        assert!(toi.is_some(), "head-on spheres should collide");
        let t = toi.unwrap();
        assert!((t - 1.0).abs() < 1e-9, "expected toi=1.0, got {}", t);
    }

    #[test]
    fn test_sphere_sphere_toi_moving_away() {
        let toi = sphere_sphere_toi(
            [0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            [3.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
        );
        assert!(
            toi.is_none(),
            "diverging spheres should not collide, got {:?}",
            toi
        );
    }

    #[test]
    fn test_sphere_plane_toi_falling() {
        // sphere at z=5.0 moving with vel=[0,0,-10], radius=0.5, plane z=0 (normal=[0,0,1], d=0)
        // toi at which z - 0.5 = 0 → z=0.5 → t = (5-0.5)/10 = 0.45
        let toi = sphere_plane_toi(
            [0.0, 0.0, 5.0],
            [0.0, 0.0, -10.0],
            0.5,
            [0.0, 0.0, 1.0],
            0.0,
        );
        assert!(toi.is_some(), "falling sphere should hit plane");
        let t = toi.unwrap();
        assert!((t - 0.45).abs() < 1e-9, "expected toi=0.45, got {}", t);
    }

    #[test]
    fn test_sphere_sphere_toi_already_overlapping() {
        // Two overlapping spheres (centers 0.5 apart, radii 0.5 each)
        let toi = sphere_sphere_toi(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
            [0.5, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
        );
        assert_eq!(toi, Some(0.0), "already overlapping → toi=0.0");
    }

    // --- ConservativeAdvancement struct tests ---

    #[test]
    fn test_conservative_advancement_struct_sphere_hit() {
        let ca = ConservativeAdvancement::new().with_threshold(1e-3);
        let toi = ca.compute_toi(
            ConservativeAdvancement::sphere_support(0.5),
            ConservativeAdvancement::sphere_support(0.5),
            [-3.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
            0.5,
            [3.0, 0.0, 0.0],
            [-2.5, 0.0, 0.0],
            0.5,
        );
        assert!(toi.is_some(), "CA struct should detect sphere–sphere hit");
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t), "toi must be in [0,1], got {}", t);
    }

    #[test]
    fn test_conservative_advancement_struct_miss() {
        let ca = ConservativeAdvancement::new();
        let toi = ca.compute_toi(
            ConservativeAdvancement::sphere_support(0.5),
            ConservativeAdvancement::sphere_support(0.5),
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.0, 5.0, 0.0],
            [1.0, 5.0, 0.0],
            0.5,
        );
        assert!(toi.is_none(), "parallel motion should not collide");
    }

    // --- integrate_linear / integrate_angular ---

    #[test]
    fn test_integrate_linear() {
        let pos = [0.0, 0.0, 0.0];
        let vel = [1.0, 2.0, 3.0];
        let result = integrate_linear(pos, vel, 0.5);
        assert!((result[0] - 0.5).abs() < 1e-12);
        assert!((result[1] - 1.0).abs() < 1e-12);
        assert!((result[2] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_integrate_angular_identity_quaternion() {
        let q = [0.0, 0.0, 0.0, 1.0];
        let omega = [0.0, 0.0, 0.0]; // no rotation
        let result = integrate_angular(q, omega, 0.1);
        // Should remain identity
        let len_sq = result[0] * result[0]
            + result[1] * result[1]
            + result[2] * result[2]
            + result[3] * result[3];
        assert!(
            (len_sq - 1.0).abs() < 1e-9,
            "quaternion should be unit length"
        );
        assert!((result[3] - 1.0).abs() < 1e-9, "w component should stay ~1");
    }

    #[test]
    fn test_quat_normalize_already_unit() {
        let q = [0.0, 0.0, 0.0, 1.0];
        let n = quat_normalize(q);
        assert!((n[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_quat_normalize_scales() {
        let q = [0.0, 0.0, 0.0, 2.0];
        let n = quat_normalize(q);
        assert!((n[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_quat_rotate_unit_vec() {
        // Rotate +X by 90° around Z → should give +Y.
        // 90° around Z: q = [0, 0, sin(45°), cos(45°)]
        let s = (std::f64::consts::FRAC_PI_4).sin();
        let c = (std::f64::consts::FRAC_PI_4).cos();
        let q = [0.0, 0.0, s, c];
        let v = [1.0, 0.0, 0.0];
        let result = quat_rotate(q, v);
        assert!(
            (result[0] - 0.0).abs() < 1e-9,
            "x should be ~0, got {}",
            result[0]
        );
        assert!(
            (result[1] - 1.0).abs() < 1e-9,
            "y should be ~1, got {}",
            result[1]
        );
        assert!(
            (result[2] - 0.0).abs() < 1e-9,
            "z should be ~0, got {}",
            result[2]
        );
    }

    // --- slerp_quat ---

    #[test]
    fn test_slerp_at_t0_returns_q0() {
        let q0 = quat_normalize([0.1, 0.2, 0.3, 0.9]);
        let q1 = quat_normalize([0.4, 0.1, 0.2, 0.8]);
        let r = slerp_quat(q0, q1, 0.0);
        for i in 0..4 {
            assert!(
                (r[i] - q0[i]).abs() < 1e-9,
                "slerp(q0,q1,0) should equal q0"
            );
        }
    }

    #[test]
    fn test_slerp_at_t1_returns_q1() {
        let q0 = quat_normalize([0.1, 0.2, 0.3, 0.9]);
        let q1 = quat_normalize([0.4, 0.1, 0.2, 0.8]);
        let r = slerp_quat(q0, q1, 1.0);
        for i in 0..4 {
            assert!(
                (r[i] - q1[i]).abs() < 1e-9,
                "slerp(q0,q1,1) should equal q1"
            );
        }
    }

    #[test]
    fn test_slerp_unit_length() {
        let q0 = quat_normalize([0.1, 0.2, 0.3, 0.9]);
        let q1 = quat_normalize([0.4, 0.1, 0.2, 0.8]);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let r = slerp_quat(q0, q1, t);
            let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-9,
                "slerp result must be unit quaternion, t={}, len={}",
                t,
                len
            );
        }
    }

    // --- RigidBodyState ---

    #[test]
    fn test_rigid_body_state_integrate() {
        let state = RigidBodyState::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        let next = state.integrate(1.0);
        assert!((next.position[0] - 1.0).abs() < 1e-9);
        assert!((next.position[1] - 0.0).abs() < 1e-9);
    }

    // --- CcdPipeline ---

    fn make_body(handle: u32, pos: [f64; 3], vel: [f64; 3], radius: f64) -> CcdBodyEntry {
        let pos_end = add3(pos, vel);
        CcdBodyEntry {
            handle,
            pos_start: pos,
            pos_end,
            orient_start: [0.0, 0.0, 0.0, 1.0],
            orient_end: [0.0, 0.0, 0.0, 1.0],
            linear_vel: vel,
            angular_vel: [0.0; 3],
            bounding_radius: radius,
            ccd_enabled: true,
        }
    }

    #[test]
    fn test_ccd_pipeline_two_approaching_bodies() {
        let pipeline = CcdPipeline::new();
        let bodies = vec![
            make_body(0, [-3.0, 0.0, 0.0], [2.5, 0.0, 0.0], 0.5),
            make_body(1, [3.0, 0.0, 0.0], [-2.5, 0.0, 0.0], 0.5),
        ];
        let pairs = pipeline.detect(&bodies);
        assert!(
            !pairs.is_empty(),
            "pipeline should detect approaching bodies"
        );
        assert!(pairs[0].toi >= 0.0 && pairs[0].toi <= 1.0);
    }

    #[test]
    fn test_ccd_pipeline_diverging_bodies_no_hit() {
        let pipeline = CcdPipeline::new();
        let bodies = vec![
            make_body(0, [0.0, 0.0, 0.0], [-5.0, 0.0, 0.0], 0.5),
            make_body(1, [3.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5),
        ];
        let pairs = pipeline.detect(&bodies);
        assert!(
            pairs.is_empty(),
            "diverging bodies should not generate CCD pair"
        );
    }

    #[test]
    fn test_ccd_pipeline_sorted_by_toi() {
        let pipeline = CcdPipeline::new();
        // Body 0 moving right; bodies 1 and 2 stationary to the right.
        // Body 1 is closer (will be hit first).
        let bodies = vec![
            make_body(0, [-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5),
            make_body(1, [-1.0, 0.0, 0.0], [-5.0, 0.0, 0.0], 0.5),
            make_body(2, [3.0, 0.0, 0.0], [-3.0, 0.0, 0.0], 0.5),
        ];
        let pairs = pipeline.detect(&bodies);
        // If multiple hits, they should be sorted by TOI.
        for i in 1..pairs.len() {
            assert!(
                pairs[i - 1].toi <= pairs[i].toi,
                "pairs must be sorted by TOI"
            );
        }
    }

    #[test]
    fn test_ccd_pipeline_ccd_disabled_body_skipped() {
        let pipeline = CcdPipeline::new();
        let mut b0 = make_body(0, [-3.0, 0.0, 0.0], [2.5, 0.0, 0.0], 0.5);
        let mut b1 = make_body(1, [3.0, 0.0, 0.0], [-2.5, 0.0, 0.0], 0.5);
        b0.ccd_enabled = false;
        b1.ccd_enabled = false;
        let bodies = vec![b0, b1];
        let pairs = pipeline.detect(&bodies);
        assert!(
            pairs.is_empty(),
            "both CCD-disabled bodies should be skipped"
        );
    }

    // --- advance_body_to ---

    #[test]
    fn test_advance_body_to_midpoint() {
        let body = make_body(0, [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5);
        let (pos, _orient) = advance_body_to(&body, 0.5);
        assert!(
            (pos[0] - 5.0).abs() < 1e-9,
            "position at t=0.5 should be 5.0"
        );
    }

    #[test]
    fn test_advance_body_to_start() {
        let body = make_body(0, [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 1.0);
        let (pos, _) = advance_body_to(&body, 0.0);
        assert!((pos[0] - 1.0).abs() < 1e-12);
        assert!((pos[1] - 2.0).abs() < 1e-12);
        assert!((pos[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_advance_body_to_end() {
        let body = make_body(0, [1.0, 0.0, 0.0], [3.0, 0.0, 0.0], 0.5);
        let (pos, _) = advance_body_to(&body, 1.0);
        assert!((pos[0] - 4.0).abs() < 1e-9);
    }
}
