// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ray casting against placed collider shapes.
//!
//! Every solid shape is intersected analytically.  The ray is expressed as
//! `origin + t * dir` with a **unit** direction, so `t` is the distance along
//! the ray.  For a solid the ray crosses the shape in one interval
//! `[enter, exit]`:
//!
//! * if the ray starts outside the shape the hit is the entry point and the
//!   reported normal is the outward surface normal there;
//! * if it starts inside the shape the hit is the exit point (again with the
//!   outward normal), so a ray fired from within a collider reports the
//!   surface it leaves through instead of missing;
//! * crossings at `t <= MIN_HIT_DISTANCE` are ignored.
//!
//! Planes are one-sided: they are hit only when the ray travels against the
//! plane normal (from the free half-space toward the solid one).

use super::geom::{Vec3, add, dot, madd, neg, normalized, scale, sub};
use super::shapes::{Shape, WorldShape};

/// Crossings closer to the origin than this are not reported.
pub(super) const MIN_HIT_DISTANCE: f64 = 1e-6;

/// Direction components smaller than this count as "parallel" to a slab.
const PARALLEL_EPS: f64 = 1e-12;

/// The nearest hit of a ray on one shape.
#[derive(Debug, Clone, Copy)]
pub(super) struct RayHit {
    /// Distance from the ray origin to the hit point.
    pub distance: f64,
    /// Outward unit surface normal at the hit point (world space).
    pub normal: Vec3,
}

/// The stretch of a ray inside a convex solid.
#[derive(Debug, Clone, Copy)]
struct Span {
    enter: f64,
    enter_normal: Vec3,
    exit: f64,
    exit_normal: Vec3,
}

/// Select the reported hit of a span (entry from outside, exit from inside).
fn pick_hit(span: Span) -> Option<(f64, Vec3)> {
    if span.enter > MIN_HIT_DISTANCE {
        Some((span.enter, span.enter_normal))
    } else if span.exit > MIN_HIT_DISTANCE {
        Some((span.exit, span.exit_normal))
    } else {
        None
    }
}

/// Intersect a ray with a placed shape.
///
/// `origin` and `dir` are in world space and `dir` must be a unit vector.
pub(super) fn ray_shape(origin: Vec3, dir: Vec3, shape: &WorldShape) -> Option<RayHit> {
    match shape.shape {
        Shape::Plane { normal, offset } => {
            let denom = dot(normal, dir);
            if denom >= -PARALLEL_EPS {
                return None;
            }
            let t = (offset - dot(normal, origin)) / denom;
            (t > MIN_HIT_DISTANCE).then_some(RayHit {
                distance: t,
                normal,
            })
        }
        _ => {
            // All remaining shapes are intersected in their local frame.
            let o = shape.to_local_point(origin);
            let d = shape.to_local_dir(dir);
            let span = match shape.shape {
                Shape::Sphere { radius } => sphere_span(o, d, radius),
                Shape::Box { half_extents } => box_span(o, d, half_extents),
                Shape::Capsule {
                    radius,
                    half_height,
                } => capsule_span(o, d, radius, half_height),
                Shape::Cylinder {
                    radius,
                    half_height,
                } => cylinder_span(o, d, radius, half_height),
                Shape::Plane { .. } => None,
            }?;
            let (distance, local_normal) = pick_hit(span)?;
            Some(RayHit {
                distance,
                normal: shape.to_world_dir(local_normal),
            })
        }
    }
}

/// Span of a ray inside a sphere centred on the local origin.
fn sphere_span(o: Vec3, d: Vec3, radius: f64) -> Option<Span> {
    let (enter, exit) = sphere_interval(o, d, [0.0; 3], radius)?;
    Some(Span {
        enter,
        enter_normal: radial_normal(o, d, enter, radius),
        exit,
        exit_normal: radial_normal(o, d, exit, radius),
    })
}

/// Outward normal of a sphere (centred on the origin) at ray parameter `t`.
fn radial_normal(o: Vec3, d: Vec3, t: f64, radius: f64) -> Vec3 {
    scale(madd(o, d, t), 1.0 / radius)
}

/// Ray parameters at which a unit-direction ray crosses a sphere.
fn sphere_interval(o: Vec3, d: Vec3, center: Vec3, radius: f64) -> Option<(f64, f64)> {
    let oc = sub(o, center);
    let b = dot(oc, d);
    let c = dot(oc, oc) - radius * radius;
    let disc = b * b - c;
    if disc < 0.0 {
        return None;
    }
    let root = disc.sqrt();
    Some((-b - root, -b + root))
}

/// Span of a ray inside an axis-aligned box centred on the local origin.
fn box_span(o: Vec3, d: Vec3, half: Vec3) -> Option<Span> {
    let mut enter = f64::NEG_INFINITY;
    let mut exit = f64::INFINITY;
    let mut enter_normal = [0.0; 3];
    let mut exit_normal = [0.0; 3];
    for axis in 0..3 {
        let h = half[axis];
        if d[axis].abs() < PARALLEL_EPS {
            if o[axis].abs() > h {
                return None;
            }
            continue;
        }
        let inv = 1.0 / d[axis];
        let mut unit = [0.0; 3];
        unit[axis] = 1.0;
        // Crossing of the slab's lower (-h) and upper (+h) planes.
        let mut t_low = (-h - o[axis]) * inv;
        let mut t_high = (h - o[axis]) * inv;
        let mut n_low = neg(unit);
        let mut n_high = unit;
        if t_low > t_high {
            std::mem::swap(&mut t_low, &mut t_high);
            std::mem::swap(&mut n_low, &mut n_high);
        }
        if t_low > enter {
            enter = t_low;
            enter_normal = n_low;
        }
        if t_high < exit {
            exit = t_high;
            exit_normal = n_high;
        }
        if enter > exit {
            return None;
        }
    }
    Some(Span {
        enter,
        enter_normal,
        exit,
        exit_normal,
    })
}

/// Ray parameters at which a unit-direction ray crosses the infinite
/// cylinder `x^2 + z^2 = r^2`.
///
/// `Some((-inf, inf))` for a ray parallel to the axis and inside the radius,
/// `None` when the ray misses the cylinder.
fn infinite_cylinder_interval(o: Vec3, d: Vec3, radius: f64) -> Option<(f64, f64)> {
    let a = d[0] * d[0] + d[2] * d[2];
    let b = o[0] * d[0] + o[2] * d[2];
    let c = o[0] * o[0] + o[2] * o[2] - radius * radius;
    if a < 1e-18 {
        return (c <= 0.0).then_some((f64::NEG_INFINITY, f64::INFINITY));
    }
    let disc = b * b - a * c;
    if disc < 0.0 {
        return None;
    }
    let root = disc.sqrt();
    Some(((-b - root) / a, (-b + root) / a))
}

/// Ray parameters at which a ray crosses the slab `|y| <= half_height`.
fn y_slab_interval(o: Vec3, d: Vec3, half_height: f64) -> Option<(f64, f64)> {
    if d[1].abs() < PARALLEL_EPS {
        return (o[1].abs() <= half_height).then_some((f64::NEG_INFINITY, f64::INFINITY));
    }
    let t0 = (-half_height - o[1]) / d[1];
    let t1 = (half_height - o[1]) / d[1];
    Some((t0.min(t1), t0.max(t1)))
}

/// Span of a ray inside a capsule around the local Y axis.
///
/// The capsule is the union of the cylindrical middle section and the two
/// end-cap spheres; because it is convex, the ray's crossing interval is the
/// hull of the intervals obtained for the three pieces.
fn capsule_span(o: Vec3, d: Vec3, radius: f64, half_height: f64) -> Option<Span> {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let mut include = |interval: Option<(f64, f64)>| {
        if let Some((a, b)) = interval
            && a <= b
        {
            lo = lo.min(a);
            hi = hi.max(b);
        }
    };
    let side = infinite_cylinder_interval(o, d, radius).and_then(|(a, b)| {
        let (s0, s1) = y_slab_interval(o, d, half_height)?;
        Some((a.max(s0), b.min(s1)))
    });
    include(side);
    include(sphere_interval(o, d, [0.0, half_height, 0.0], radius));
    include(sphere_interval(o, d, [0.0, -half_height, 0.0], radius));
    if lo > hi {
        return None;
    }
    let normal_at = |t: f64| -> Vec3 {
        let p = madd(o, d, t);
        let axis_point = [0.0, p[1].clamp(-half_height, half_height), 0.0];
        normalized(sub(p, axis_point)).unwrap_or([1.0, 0.0, 0.0])
    };
    Some(Span {
        enter: lo,
        enter_normal: normal_at(lo),
        exit: hi,
        exit_normal: normal_at(hi),
    })
}

/// Span of a ray inside a solid cylinder around the local Y axis.
fn cylinder_span(o: Vec3, d: Vec3, radius: f64, half_height: f64) -> Option<Span> {
    let (r0, r1) = infinite_cylinder_interval(o, d, radius)?;
    let (s0, s1) = y_slab_interval(o, d, half_height)?;
    let enter = r0.max(s0);
    let exit = r1.min(s1);
    if enter > exit {
        return None;
    }
    let radial_at = |t: f64| -> Vec3 {
        let p = madd(o, d, t);
        normalized([p[0], 0.0, p[2]]).unwrap_or([1.0, 0.0, 0.0])
    };
    // A ray travelling toward +y enters through the bottom cap and leaves
    // through the top cap (and the other way round when it travels toward -y).
    let (bottom, top) = ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0]);
    let (cap_enter, cap_exit) = if d[1] > 0.0 {
        (bottom, top)
    } else {
        (top, bottom)
    };
    let enter_normal = if r0 >= s0 {
        radial_at(enter)
    } else {
        cap_enter
    };
    let exit_normal = if r1 <= s1 { radial_at(exit) } else { cap_exit };
    Some(Span {
        enter,
        enter_normal,
        exit,
        exit_normal,
    })
}

/// Point on a ray at parameter `t` (helper shared with the engine).
pub(super) fn point_at(origin: Vec3, dir: Vec3, t: f64) -> Vec3 {
    add(origin, scale(dir, t))
}
