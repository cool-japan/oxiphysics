// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shape casting (continuous collision detection / swept queries).
//!
//! Provides time-of-impact (TOI) computation for moving shapes: spheres,
//! capsules, boxes, and general convex shapes. Uses `[f64; 3]` for 3D vectors.
//!
//! The casts translate one shape along a straight line (`position(t) = pos +
//! t * vel`) against a static one and find the first `t` at which their gap
//! reaches [`ShapeCastConfig::target_distance`] by conservative advancement:
//! the exact distance between the shapes is measured, and time advances by the
//! largest step that cannot skip the contact. The step is the gap divided by
//! the speed at which the gap closes along the direction of the closest
//! features (a Newton step on the convex gap function), see the
//! [`advancement`] module, so oblique approaches converge as fast as head-on
//! ones. The distances are exact: segment-box and box-box queries are closed
//! form, and convex hulls given as point sets use the
//! Gilbert-Johnson-Keerthi algorithm ([`hull_distance`]).

mod advancement;
mod box_queries;
mod hull_distance;

use advancement::{Probe, advance, separation_direction, unit_vector};
use box_queries::{
    aabb_overlap_normal, aabb_separation, box_corners, box_exit_normal, clamp_to_box,
    segment_box_distance,
};
use hull_distance::{HullDistance, hull_distance};

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a).max(1e-15);
    scale3(a, 1.0 / l)
}

fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    add3(scale3(a, 1.0 - t), scale3(b, t))
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeCastConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for shape casting queries.
#[derive(Debug, Clone, Copy)]
pub struct ShapeCastConfig {
    /// Maximum time of contact to consider (normalized to \[0, 1\]).
    pub max_toc: f64,
    /// Target separation distance (shapes are in contact when dist <= target_distance).
    pub target_distance: f64,
    /// Distance tolerance for time-of-impact convergence: a cast reports a hit
    /// as soon as the gap between the shapes is at most
    /// `target_distance + toi_tolerance`.
    pub toi_tolerance: f64,
}

impl Default for ShapeCastConfig {
    fn default() -> Self {
        Self {
            max_toc: 1.0,
            target_distance: 0.0,
            toi_tolerance: 1e-6,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeCastResult
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a shape cast query.
#[derive(Debug, Clone)]
pub struct ShapeCastResult {
    /// Normalized time of impact ∈ \[0, max_toc\], or `max_toc` if no hit.
    pub toi: f64,
    /// Contact normal at the time of impact (pointing from B to A).
    pub normal_at_toi: [f64; 3],
    /// Witness point on shape A at TOI: the closest point of A to B (for a
    /// capsule cast, the closest point of the capsule's axis segment).
    pub witness_a: [f64; 3],
    /// Witness point on shape B at TOI: the closest point of B to A (for a
    /// capsule-capsule cast, the closest point of the other axis segment).
    pub witness_b: [f64; 3],
    /// True if the shapes were already penetrating at t=0.
    pub is_penetrating: bool,
    /// True if a hit was found within max_toc.
    pub hit: bool,
}

impl ShapeCastResult {
    /// Creates a hit result.
    pub fn hit(
        toi: f64,
        normal_at_toi: [f64; 3],
        witness_a: [f64; 3],
        witness_b: [f64; 3],
    ) -> Self {
        Self {
            toi,
            normal_at_toi,
            witness_a,
            witness_b,
            is_penetrating: false,
            hit: true,
        }
    }

    /// Creates a penetrating (initial overlap) result.
    pub fn penetrating(witness_a: [f64; 3], witness_b: [f64; 3], normal: [f64; 3]) -> Self {
        Self {
            toi: 0.0,
            normal_at_toi: normal,
            witness_a,
            witness_b,
            is_penetrating: true,
            hit: true,
        }
    }

    /// Creates a no-hit result.
    pub fn no_hit(max_toc: f64) -> Self {
        Self {
            toi: max_toc,
            normal_at_toi: [0.0, 1.0, 0.0],
            witness_a: [0.0; 3],
            witness_b: [0.0; 3],
            is_penetrating: false,
            hit: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LinearMotion
// ─────────────────────────────────────────────────────────────────────────────

/// Linear motion of a rigid body: position(t) = pos + t * vel.
#[derive(Debug, Clone, Copy)]
pub struct LinearMotion {
    /// Initial position.
    pub pos: [f64; 3],
    /// Velocity vector.
    pub vel: [f64; 3],
}

impl LinearMotion {
    /// Creates a new linear motion.
    pub fn new(pos: [f64; 3], vel: [f64; 3]) -> Self {
        Self { pos, vel }
    }

    /// Position at time t.
    pub fn position_at(&self, t: f64) -> [f64; 3] {
        add3(self.pos, scale3(self.vel, t))
    }

    /// Relative velocity between two motions.
    pub fn relative_velocity(a: &Self, b: &Self) -> [f64; 3] {
        sub3(a.vel, b.vel)
    }

    /// Conservative TOI upper bound for two linearly moving spheres.
    ///
    /// Returns None if they will not collide within \[0, max_t\].
    pub fn compute_toi_linear_sphere(
        motion_a: &LinearMotion,
        radius_a: f64,
        motion_b: &LinearMotion,
        radius_b: f64,
        max_t: f64,
    ) -> Option<f64> {
        let rel_pos = sub3(motion_a.pos, motion_b.pos);
        let rel_vel = sub3(motion_a.vel, motion_b.vel);
        let r = radius_a + radius_b;

        let a = dot3(rel_vel, rel_vel);
        let b = 2.0 * dot3(rel_pos, rel_vel);
        let c = dot3(rel_pos, rel_pos) - r * r;

        if c <= 0.0 {
            // Already overlapping
            return Some(0.0);
        }
        if a < 1e-15 {
            return None;
        }
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let t = (-b - disc.sqrt()) / (2.0 * a);
        if t >= 0.0 && t <= max_t {
            Some(t)
        } else {
            None
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AngularMotion
// ─────────────────────────────────────────────────────────────────────────────

/// Combined translation + rotation motion.
#[derive(Debug, Clone, Copy)]
pub struct AngularMotion {
    /// Initial position.
    pub pos: [f64; 3],
    /// Linear velocity.
    pub vel: [f64; 3],
    /// Angular velocity (rotation axis * angular speed).
    pub angular_vel: [f64; 3],
    /// Initial orientation (unit quaternion \[x, y, z, w\]).
    pub orientation: [f64; 4],
}

impl AngularMotion {
    /// Creates a new angular motion (no rotation by default).
    pub fn new(pos: [f64; 3], vel: [f64; 3]) -> Self {
        Self {
            pos,
            vel,
            angular_vel: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Creates a motion with both translation and rotation.
    pub fn with_rotation(pos: [f64; 3], vel: [f64; 3], angular_vel: [f64; 3]) -> Self {
        Self {
            pos,
            vel,
            angular_vel,
            orientation: [0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Position at time t (linear component only).
    pub fn position_at(&self, t: f64) -> [f64; 3] {
        add3(self.pos, scale3(self.vel, t))
    }

    /// Computes an upper bound on the displacement magnitude over \[0, t\].
    pub fn motion_bound(&self, t: f64, point_offset: f64) -> f64 {
        let linear = len3(self.vel) * t;
        let angular = len3(self.angular_vel) * t * point_offset;
        linear + angular
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Probe helpers shared by the casts
// ─────────────────────────────────────────────────────────────────────────────

/// Probe of a sphere (centre `pos`, `radius`) against the point `nearest` of a
/// static convex shape that is closest to the centre.
///
/// The witness on the sphere is the point of its surface facing `nearest`.
fn sphere_probe(
    pos: [f64; 3],
    radius: f64,
    nearest: [f64; 3],
    fallback: Option<[f64; 3]>,
) -> Probe {
    let direction = separation_direction(pos, nearest);
    let facing = direction.or(fallback);
    Probe {
        gap: dist3(pos, nearest) - radius,
        direction,
        fallback,
        witness_a: facing.map_or(pos, |n| sub3(pos, scale3(n, radius))),
        witness_b: nearest,
    }
}

/// Probe of two convex hulls given by the distance query `nearest`.
fn hull_probe(nearest: &HullDistance, fallback: Option<[f64; 3]>) -> Probe {
    Probe {
        gap: nearest.distance,
        direction: if nearest.intersecting {
            None
        } else {
            separation_direction(nearest.point_a, nearest.point_b)
        },
        fallback,
        witness_a: nearest.point_a,
        witness_b: nearest.point_b,
    }
}

/// Unit normal of the triangle `(v0, v1, v2)`, oriented against the motion
/// `vel` (so that it points towards the side a mover comes from), or `None` for
/// a degenerate triangle.
fn triangle_normal_against(
    vel: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> Option<[f64; 3]> {
    let normal = unit_vector(cross3(sub3(v1, v0), sub3(v2, v0)))?;
    Some(if dot3(normal, vel) > 0.0 {
        scale3(normal, -1.0)
    } else {
        normal
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// SphereCast — moving sphere cast
// ─────────────────────────────────────────────────────────────────────────────

/// Cast a moving sphere against static primitives.
pub struct SphereCast;

impl SphereCast {
    /// Sphere vs sphere TOI computation.
    pub fn vs_sphere(
        motion_a: &LinearMotion,
        radius_a: f64,
        center_b: [f64; 3],
        radius_b: f64,
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let static_motion = LinearMotion::new(center_b, [0.0; 3]);
        match LinearMotion::compute_toi_linear_sphere(
            motion_a,
            radius_a,
            &static_motion,
            radius_b,
            cfg.max_toc,
        ) {
            Some(toi) if toi <= 0.0 => {
                let n = normalize3(sub3(motion_a.pos, center_b));
                ShapeCastResult::penetrating(motion_a.pos, center_b, n)
            }
            Some(toi) => {
                let pos_a = motion_a.position_at(toi);
                let n = normalize3(sub3(pos_a, center_b));
                let wa = add3(pos_a, scale3(n, -radius_a));
                let wb = add3(center_b, scale3(n, radius_b));
                ShapeCastResult::hit(toi, n, wa, wb)
            }
            None => ShapeCastResult::no_hit(cfg.max_toc),
        }
    }

    /// Sphere vs axis-aligned box TOI via conservative advancement.
    ///
    /// The witness on the sphere is the point of its surface facing the box.
    pub fn vs_box(
        motion_a: &LinearMotion,
        radius_a: f64,
        box_center: [f64; 3],
        half_extents: [f64; 3],
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let lo = sub3(box_center, half_extents);
        let hi = add3(box_center, half_extents);
        advance(motion_a.vel, cfg, |t| {
            let pos = motion_a.position_at(t);
            let nearest = clamp_to_box(pos, lo, hi);
            Some(sphere_probe(
                pos,
                radius_a,
                nearest,
                Some(box_exit_normal(pos, lo, hi)),
            ))
        })
    }

    /// Sphere vs triangle TOI via conservative advancement.
    ///
    /// The witness on the sphere is the point of its surface facing the triangle.
    pub fn vs_triangle(
        motion_a: &LinearMotion,
        radius_a: f64,
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let fallback = triangle_normal_against(motion_a.vel, v0, v1, v2);
        advance(motion_a.vel, cfg, |t| {
            let pos = motion_a.position_at(t);
            let (nearest, _) = crate_point_tri(pos, v0, v1, v2);
            Some(sphere_probe(pos, radius_a, nearest, fallback))
        })
    }
}

/// Local copy of the point-triangle closest point logic (avoids cross-crate dep).
fn crate_point_tri(p: [f64; 3], v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let ab = sub3(v1, v0);
    let ac = sub3(v2, v0);
    let ap = sub3(p, v0);

    let d1 = dot3(ab, ap);
    let d2 = dot3(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (v0, [1.0, 0.0, 0.0]);
    }
    let bp = sub3(p, v1);
    let d3 = dot3(ab, bp);
    let d4 = dot3(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (v1, [0.0, 1.0, 0.0]);
    }
    let cp = sub3(p, v2);
    let d5 = dot3(ab, cp);
    let d6 = dot3(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (v2, [0.0, 0.0, 1.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let q = add3(v0, scale3(ab, v));
        return (q, [1.0 - v, v, 0.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let q = add3(v0, scale3(ac, w));
        return (q, [1.0 - w, 0.0, w]);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let q = lerp3(v1, v2, w);
        return (q, [0.0, 1.0 - w, w]);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let u = 1.0 - v - w;
    let q = add3(add3(scale3(v0, u), scale3(v1, v)), scale3(v2, w));
    (q, [u, v, w])
}

// ─────────────────────────────────────────────────────────────────────────────
// CapsuleCast — moving capsule cast
// ─────────────────────────────────────────────────────────────────────────────

/// Cast a moving capsule against static primitives.
pub struct CapsuleCast;

impl CapsuleCast {
    /// Capsule vs capsule TOI via conservative advancement.
    ///
    /// Each capsule is given by its axis segment and radius. The witnesses are
    /// the closest points of the two axes at the time of impact; moving
    /// `radius` along the normal from them reaches the capsule surfaces.
    pub fn vs_capsule(
        motion: &LinearMotion,
        local_a: [f64; 3],
        local_b: [f64; 3],
        radius_a: f64,
        static_a: [f64; 3],
        static_b: [f64; 3],
        radius_b: f64,
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let r = radius_a + radius_b;
        advance(motion.vel, cfg, |t| {
            let offset = motion.position_at(t);
            let (seg_dist, wa, wb) = seg_seg_dist(
                add3(local_a, offset),
                add3(local_b, offset),
                static_a,
                static_b,
            );
            Some(Probe {
                gap: seg_dist - r,
                direction: separation_direction(wa, wb),
                fallback: None,
                witness_a: wa,
                witness_b: wb,
            })
        })
    }

    /// Capsule vs axis-aligned box TOI via conservative advancement.
    ///
    /// The capsule is given by its axis segment (`local_a`..`local_b`, relative
    /// to the motion's position) and `radius`. The distance used for the
    /// advancement is the exact distance between the whole axis segment and the
    /// box, so a capsule whose middle hits a face or an edge is detected just
    /// like one whose end caps do. `witness_a` is the closest point of the axis
    /// segment and `witness_b` the closest point of the box; moving `radius`
    /// from `witness_a` along `-normal_at_toi` reaches the capsule surface.
    pub fn vs_box(
        motion: &LinearMotion,
        local_a: [f64; 3],
        local_b: [f64; 3],
        radius: f64,
        box_center: [f64; 3],
        half_extents: [f64; 3],
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let lo = sub3(box_center, half_extents);
        let hi = add3(box_center, half_extents);
        advance(motion.vel, cfg, |t| {
            let offset = motion.position_at(t);
            let nearest =
                segment_box_distance(add3(local_a, offset), add3(local_b, offset), lo, hi);
            Some(Probe {
                gap: nearest.distance - radius,
                direction: separation_direction(nearest.segment_point, nearest.box_point),
                fallback: Some(box_exit_normal(nearest.segment_point, lo, hi)),
                witness_a: nearest.segment_point,
                witness_b: nearest.box_point,
            })
        })
    }
}

/// Local segment-segment closest distance.
///
/// Returns `(distance, point_on_p1p2, point_on_p3p4)`. This is the clamped
/// closest-point computation from Ericson, *Real-Time Collision Detection*
/// (5.1.9): when the parameter on the second segment has to be clamped to one of
/// its end points, the parameter on the first segment is recomputed against that
/// end point, so the result is the true minimum-distance pair.
fn seg_seg_dist(
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    p4: [f64; 3],
) -> (f64, [f64; 3], [f64; 3]) {
    let d1 = sub3(p2, p1);
    let d2 = sub3(p4, p3);
    let r = sub3(p1, p3);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    if a < 1e-15 && e < 1e-15 {
        return (dist3(p1, p3), p1, p3);
    }
    // `s` parameterises p1..p2 and `t` parameterises p3..p4.
    let (s, t) = if a < 1e-15 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = dot3(d1, r);
        if e < 1e-15 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = dot3(d1, d2);
            let denom = a * e - b * b;
            // Point on segment 1 closest to the infinite line through segment 2
            // (arbitrary when the segments are parallel).
            let s_line = if denom > 1e-12 * a * e {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t_nom = b * s_line + f;
            if t_nom < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t_nom > e {
                (((b - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s_line, t_nom / e)
            }
        }
    };
    let wa = lerp3(p1, p2, s);
    let wb = lerp3(p3, p4, t);
    (dist3(wa, wb), wa, wb)
}

// ─────────────────────────────────────────────────────────────────────────────
// BoxCast — oriented box cast
// ─────────────────────────────────────────────────────────────────────────────

/// Cast a moving axis-aligned box against static primitives.
pub struct BoxCast;

impl BoxCast {
    /// Box vs box TOI via conservative advancement.
    ///
    /// The witnesses are the closest points of the two boxes; where their
    /// extents overlap on an axis both sit in the middle of the overlap.
    pub fn vs_box(
        motion: &LinearMotion,
        half_a: [f64; 3],
        center_b: [f64; 3],
        half_b: [f64; 3],
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let lo_b = sub3(center_b, half_b);
        let hi_b = add3(center_b, half_b);
        advance(motion.vel, cfg, |t| {
            let pos_a = motion.position_at(t);
            let nearest = aabb_separation(sub3(pos_a, half_a), add3(pos_a, half_a), lo_b, hi_b);
            Some(Probe {
                gap: nearest.distance,
                direction: nearest.direction,
                fallback: Some(aabb_overlap_normal(pos_a, half_a, center_b, half_b)),
                witness_a: nearest.point_a,
                witness_b: nearest.point_b,
            })
        })
    }

    /// Box vs triangle TOI via conservative advancement.
    ///
    /// The advancement uses the exact distance between the box and the whole
    /// triangle (face, edges and vertices), computed with GJK on the box
    /// corners and the triangle vertices, so a box that lands on the interior
    /// of a large triangle is detected at the correct time. The witnesses are
    /// the closest points of the box and of the triangle.
    pub fn vs_triangle(
        motion: &LinearMotion,
        half_extents: [f64; 3],
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let corners = box_corners(half_extents);
        let triangle = [v0, v1, v2];
        let fallback = triangle_normal_against(motion.vel, v0, v1, v2);
        advance(motion.vel, cfg, |t| {
            let nearest = hull_distance(&corners, motion.position_at(t), &triangle)?;
            Some(hull_probe(&nearest, fallback))
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ConvexCast — GJK-based convex shape cast
// ─────────────────────────────────────────────────────────────────────────────

/// GJK-based convex shape cast using conservative advancement.
///
/// Shapes are represented as point clouds: the convex hull of the points is the
/// shape. The distance between the hulls is exact (Gilbert-Johnson-Keerthi on
/// the Minkowski difference), so the cast does not depend on how densely the
/// hulls are sampled.
pub struct ConvexCast;

impl ConvexCast {
    /// Convex cast via conservative advancement.
    ///
    /// `shape_a` is moved by `motion.position_at(t)` (its points are relative to
    /// the motion's position); `shape_b` is static. The witnesses are the
    /// closest points of the two hulls. If either point set is empty there is
    /// nothing to hit and the result is a miss.
    pub fn cast(
        motion: &LinearMotion,
        shape_a: &[[f64; 3]],
        shape_b: &[[f64; 3]],
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        advance(motion.vel, cfg, |t| {
            let nearest = hull_distance(shape_a, motion.position_at(t), shape_b)?;
            Some(hull_probe(&nearest, None))
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SubstepCast — tunneling prevention via substep subdivision
// ─────────────────────────────────────────────────────────────────────────────

/// Prevents tunneling by subdividing the motion into smaller substeps.
pub struct SubstepCast {
    /// Minimum substep size (dt) to prevent infinite subdivision.
    pub min_substep: f64,
    /// Maximum number of substeps.
    pub max_substeps: usize,
}

impl SubstepCast {
    /// Creates a new substep caster.
    pub fn new(min_substep: f64, max_substeps: usize) -> Self {
        Self {
            min_substep,
            max_substeps,
        }
    }

    /// Runs a shape cast query with tunneling prevention.
    ///
    /// `cast_fn` performs the actual cast for a sub-interval.
    pub fn cast<F>(&self, t_start: f64, t_end: f64, cast_fn: F) -> Option<ShapeCastResult>
    where
        F: Fn(f64, f64) -> ShapeCastResult,
    {
        let dt = t_end - t_start;
        let n_steps = ((dt / self.min_substep).ceil() as usize)
            .min(self.max_substeps)
            .max(1);
        let step = dt / n_steps as f64;
        for i in 0..n_steps {
            let ta = t_start + i as f64 * step;
            let tb = ta + step;
            let result = cast_fn(ta, tb);
            if result.hit {
                return Some(result);
            }
        }
        None
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NonLinearCast — general non-linear motion cast
// ─────────────────────────────────────────────────────────────────────────────

/// Non-linear motion cast via time integration and bisection.
///
/// Handles arbitrary smooth motions (e.g., ballistic, pendulum) by
/// stepping forward and bisecting when contact is detected.
pub struct NonLinearCast {
    /// Initial time step for integration.
    pub dt: f64,
    /// Tolerance for contact detection.
    pub tol: f64,
    /// Maximum number of bisection steps.
    pub max_bisect: usize,
}

impl NonLinearCast {
    /// Creates a new non-linear cast solver.
    pub fn new(dt: f64, tol: f64, max_bisect: usize) -> Self {
        Self {
            dt,
            tol,
            max_bisect,
        }
    }

    /// Casts using a provided distance function `dist_fn(t) -> distance`.
    ///
    /// Returns the approximate TOI or None if no contact within \[0, t_max\].
    pub fn cast<F>(&self, t_max: f64, dist_fn: F) -> Option<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut t = 0.0;
        let d0 = dist_fn(0.0);
        if d0 <= self.tol {
            return Some(0.0);
        }

        while t < t_max {
            t += self.dt;
            let t_cur = t.min(t_max);
            let d = dist_fn(t_cur);
            if d <= self.tol {
                // Bisect in [t - dt, t_cur]
                let mut lo = (t - self.dt).max(0.0);
                let mut hi = t_cur;
                for _ in 0..self.max_bisect {
                    let mid = 0.5 * (lo + hi);
                    let d_mid = dist_fn(mid);
                    if d_mid <= self.tol {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                    if hi - lo < self.tol * 0.1 {
                        break;
                    }
                }
                return Some(hi);
            }
        }
        None
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MeshBvhCast — BVH-traversal cast against mesh
// ─────────────────────────────────────────────────────────────────────────────

/// Shape cast against a complex mesh using a simple BVH traversal.
pub struct MeshBvhCast {
    /// Mesh vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices.
    pub indices: Vec<usize>,
}

impl MeshBvhCast {
    /// Creates a new mesh BVH cast structure.
    pub fn new(vertices: Vec<[f64; 3]>, indices: Vec<usize>) -> Self {
        assert!(indices.len().is_multiple_of(3));
        Self { vertices, indices }
    }

    /// Casts a moving sphere against the mesh.
    ///
    /// Returns the earliest TOI ShapeCastResult.
    pub fn sphere_cast(
        &self,
        motion: &LinearMotion,
        radius: f64,
        cfg: &ShapeCastConfig,
    ) -> ShapeCastResult {
        let n_tris = self.indices.len() / 3;
        let mut best = ShapeCastResult::no_hit(cfg.max_toc);

        for tri in 0..n_tris {
            let v0 = self.vertices[self.indices[3 * tri]];
            let v1 = self.vertices[self.indices[3 * tri + 1]];
            let v2 = self.vertices[self.indices[3 * tri + 2]];
            let result = SphereCast::vs_triangle(motion, radius, v0, v1, v2, cfg);
            if result.hit && (!best.hit || result.toi < best.toi) {
                best = result;
            }
        }
        best
    }

    /// Number of triangles in the mesh.
    pub fn n_triangles(&self) -> usize {
        self.indices.len() / 3
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod regression_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn default_cfg() -> ShapeCastConfig {
        ShapeCastConfig::default()
    }

    #[test]
    fn test_shape_cast_config_default() {
        let cfg = ShapeCastConfig::default();
        assert!((cfg.max_toc - 1.0).abs() < 1e-12);
        assert_eq!(cfg.target_distance, 0.0);
    }

    #[test]
    fn test_linear_motion_position_at() {
        let m = LinearMotion::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = m.position_at(2.0);
        assert!((p[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_linear_motion_toi_sphere_hit() {
        let a = LinearMotion::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let b = LinearMotion::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let toi = LinearMotion::compute_toi_linear_sphere(&a, 0.5, &b, 0.5, 10.0);
        assert!(toi.is_some());
        assert!((toi.unwrap() - 4.0).abs() < 0.1);
    }

    #[test]
    fn test_linear_motion_toi_sphere_miss() {
        let a = LinearMotion::new([10.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let b = LinearMotion::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let toi = LinearMotion::compute_toi_linear_sphere(&a, 0.5, &b, 0.5, 10.0);
        assert!(toi.is_none());
    }

    #[test]
    fn test_linear_motion_toi_sphere_initial_overlap() {
        let a = LinearMotion::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let b = LinearMotion::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let toi = LinearMotion::compute_toi_linear_sphere(&a, 1.0, &b, 1.0, 10.0);
        assert!(toi.is_some());
        assert!(toi.unwrap() <= 0.0);
    }

    #[test]
    fn test_sphere_cast_vs_sphere_hit() {
        // Sphere at [0,0,-1.5] moving +z at speed 1. Both radii 0.5.
        // Contact at t = 1.5 - 1.0 = 0.5 (within max_toc=1.0).
        let motion = LinearMotion::new([0.0, 0.0, -1.5], [0.0, 0.0, 1.0]);
        let cfg = default_cfg();
        let result = SphereCast::vs_sphere(&motion, 0.5, [0.0, 0.0, 0.0], 0.5, &cfg);
        assert!(result.hit);
        assert!(result.toi > 0.0 && result.toi <= 1.0);
    }

    #[test]
    fn test_sphere_cast_vs_sphere_miss() {
        let motion = LinearMotion::new([10.0, 0.0, -5.0], [0.0, 0.0, -1.0]);
        let cfg = default_cfg();
        let result = SphereCast::vs_sphere(&motion, 0.5, [0.0, 0.0, 0.0], 0.5, &cfg);
        assert!(!result.hit);
    }

    #[test]
    fn test_sphere_cast_vs_sphere_penetrating() {
        let motion = LinearMotion::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let cfg = default_cfg();
        let result = SphereCast::vs_sphere(&motion, 1.0, [0.0, 0.0, 0.0], 1.0, &cfg);
        assert!(result.hit);
        assert!(result.is_penetrating);
    }

    #[test]
    fn test_sphere_cast_vs_box_hit() {
        let motion = LinearMotion::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cfg = ShapeCastConfig {
            max_toc: 10.0,
            ..default_cfg()
        };
        let result = SphereCast::vs_box(&motion, 0.5, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], &cfg);
        assert!(result.hit);
    }

    #[test]
    fn test_capsule_cast_vs_capsule() {
        let motion = LinearMotion::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cfg = ShapeCastConfig {
            max_toc: 10.0,
            ..default_cfg()
        };
        let result = CapsuleCast::vs_capsule(
            &motion,
            [-0.5, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.25,
            [-0.5, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.25,
            &cfg,
        );
        assert!(result.hit);
    }

    #[test]
    fn test_box_cast_vs_box_hit() {
        let motion = LinearMotion::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cfg = ShapeCastConfig {
            max_toc: 10.0,
            ..default_cfg()
        };
        let result = BoxCast::vs_box(
            &motion,
            [0.5, 0.5, 0.5],
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            &cfg,
        );
        assert!(result.hit);
    }

    #[test]
    fn test_box_cast_vs_box_miss() {
        let motion = LinearMotion::new([10.0, 10.0, -5.0], [0.0, 0.0, -1.0]);
        let cfg = default_cfg();
        let result = BoxCast::vs_box(
            &motion,
            [0.1, 0.1, 0.1],
            [0.0, 0.0, 0.0],
            [0.1, 0.1, 0.1],
            &cfg,
        );
        assert!(!result.hit);
    }

    #[test]
    fn test_convex_cast_hit() {
        let motion = LinearMotion::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cfg = ShapeCastConfig {
            max_toc: 10.0,
            ..default_cfg()
        };
        let shape_a = vec![[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let shape_b = vec![[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let result = ConvexCast::cast(&motion, &shape_a, &shape_b, &cfg);
        assert!(result.hit);
    }

    #[test]
    fn test_angular_motion_position_at() {
        let m = AngularMotion::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = m.position_at(1.0);
        assert!((p[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_angular_motion_bound() {
        let m = AngularMotion::with_rotation([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let bound = m.motion_bound(1.0, 0.5);
        assert!(bound > 1.0); // at least linear part
    }

    #[test]
    fn test_substep_cast_finds_hit() {
        let caster = SubstepCast::new(0.1, 100);
        let result = caster.cast(0.0, 1.0, |ta, _tb| {
            if ta > 0.5 {
                ShapeCastResult::hit(ta, [0.0, 1.0, 0.0], [0.0; 3], [0.0; 3])
            } else {
                ShapeCastResult::no_hit(1.0)
            }
        });
        assert!(result.is_some());
    }

    #[test]
    fn test_substep_cast_no_hit() {
        let caster = SubstepCast::new(0.1, 10);
        let result = caster.cast(0.0, 1.0, |_ta, _tb| ShapeCastResult::no_hit(1.0));
        assert!(result.is_none());
    }

    #[test]
    fn test_nonlinear_cast_hit() {
        let caster = NonLinearCast::new(0.01, 0.05, 20);
        // Distance decreases linearly: d(t) = 1.0 - t
        let toi = caster.cast(2.0, |t| 1.0 - t);
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!((t - 1.0).abs() < 0.1, "toi={}", t);
    }

    #[test]
    fn test_nonlinear_cast_no_hit() {
        let caster = NonLinearCast::new(0.1, 0.05, 10);
        // Distance always > 1
        let toi = caster.cast(1.0, |_t| 2.0);
        assert!(toi.is_none());
    }

    #[test]
    fn test_mesh_bvh_cast_sphere() {
        let verts = vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let idx = vec![0usize, 1, 2];
        let mesh = MeshBvhCast::new(verts, idx);
        let motion = LinearMotion::new([0.0, 0.5, -3.0], [0.0, 0.0, 1.0]);
        let cfg = ShapeCastConfig {
            max_toc: 10.0,
            ..default_cfg()
        };
        let result = mesh.sphere_cast(&motion, 0.2, &cfg);
        assert!(result.hit);
    }

    #[test]
    fn test_shape_cast_result_no_hit() {
        let r = ShapeCastResult::no_hit(1.0);
        assert!(!r.hit);
        assert!(!r.is_penetrating);
    }

    #[test]
    fn test_shape_cast_result_penetrating() {
        let r = ShapeCastResult::penetrating([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]);
        assert!(r.hit);
        assert!(r.is_penetrating);
        assert_eq!(r.toi, 0.0);
    }

    #[test]
    fn test_relative_velocity() {
        let a = LinearMotion::new([0.0; 3], [1.0, 0.0, 0.0]);
        let b = LinearMotion::new([0.0; 3], [0.5, 0.0, 0.0]);
        let rel = LinearMotion::relative_velocity(&a, &b);
        assert!((rel[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_sphere_vs_triangle_no_hit() {
        let motion = LinearMotion::new([10.0, 10.0, -5.0], [0.0, 0.0, -1.0]);
        let cfg = default_cfg();
        let result = SphereCast::vs_triangle(
            &motion,
            0.1,
            [-0.1, -0.1, 0.0],
            [0.1, -0.1, 0.0],
            [0.0, 0.1, 0.0],
            &cfg,
        );
        assert!(!result.hit);
    }

    #[test]
    fn test_seg_seg_dist_recomputes_first_parameter_after_clamp() {
        // Segment 1: (0,0,0)-(10,0,0). Segment 2: (5,1,0)-(6,2,0).
        // The infinite-line closest pair clamps segment 2 to its start point;
        // segment 1 must then be re-projected onto (5,1,0), giving (5,0,0) and a
        // distance of exactly 1 (not sqrt(2) at (4,0,0)).
        let (dist, wa, wb) = seg_seg_dist(
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
            [6.0, 2.0, 0.0],
        );
        assert!((dist - 1.0).abs() < 1e-12, "dist={dist}");
        assert!(len3(sub3(wa, [5.0, 0.0, 0.0])) < 1e-12, "wa={wa:?}");
        assert!(len3(sub3(wb, [5.0, 1.0, 0.0])) < 1e-12, "wb={wb:?}");
    }

    #[test]
    fn test_seg_seg_dist_degenerate_and_parallel_inputs() {
        // Both segments are points.
        let (d, _, _) = seg_seg_dist([0.0; 3], [0.0; 3], [3.0, 4.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-12);
        // First segment is a point; second is clamped at its end point.
        let (d, _, wb) = seg_seg_dist([9.0, 1.0, 0.0], [9.0, 1.0, 0.0], [0.0; 3], [2.0, 0.0, 0.0]);
        assert!((d - (49.0_f64 + 1.0).sqrt()).abs() < 1e-12, "d={d}");
        assert!(len3(sub3(wb, [2.0, 0.0, 0.0])) < 1e-12);
        // Parallel overlapping segments one unit apart.
        let (d, _, _) = seg_seg_dist([0.0; 3], [4.0, 0.0, 0.0], [1.0, 1.0, 0.0], [3.0, 1.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-12, "d={d}");
    }

    #[test]
    fn test_capsule_cast_vs_capsule_end_point_clamped_pair() {
        // The long capsule A (axis (0,0,0)-(10,0,0), radius 0.25) sweeps up +y
        // towards the short static capsule B (axis (5,1,0)-(6,2,0), radius 0.25).
        // B's lower end point (5,1,0) is the nearest feature, so the surfaces
        // touch when the axes are 0.5 apart, i.e. at t = 0.5. The buggy
        // segment-segment distance overestimated the gap by sqrt(2) and reported
        // a late hit at t ~ 0.91.
        let motion = LinearMotion::new([0.0; 3], [0.0, 1.0, 0.0]);
        let cfg = ShapeCastConfig {
            max_toc: 10.0,
            ..default_cfg()
        };
        let result = CapsuleCast::vs_capsule(
            &motion,
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            0.25,
            [5.0, 1.0, 0.0],
            [6.0, 2.0, 0.0],
            0.25,
            &cfg,
        );
        assert!(result.hit);
        assert!((result.toi - 0.5).abs() < 1e-9, "toi={}", result.toi);
        assert!(len3(sub3(result.witness_a, [5.0, 0.5, 0.0])) < 1e-9);
        assert!(len3(sub3(result.witness_b, [5.0, 1.0, 0.0])) < 1e-9);
        assert!(
            result.normal_at_toi[1] < -0.99,
            "{:?}",
            result.normal_at_toi
        );
    }
}
