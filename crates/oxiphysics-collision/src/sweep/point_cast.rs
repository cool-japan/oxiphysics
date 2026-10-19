// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time of impact of a moving point against the pieces of a swept sphere's
//! offset surface.
//!
//! A sphere of radius `r` whose centre follows `c(t) = origin + t * dir` for
//! `t` in `[0, 1]` touches a convex body exactly when the centre reaches the
//! body's offset surface at distance `r`. That surface is assembled from
//! spherical caps around the vertices, finite cylinders around the edges and,
//! for a triangle, two offset planes. The functions in this module find the
//! first time at which the moving point enters one such piece from outside;
//! the caller takes the earliest entry over all pieces of the body.
//!
//! All quadratic roots use the cancellation-free form `c / (-b + sqrt(disc))`
//! for the smaller root, which stays accurate for grazing and near-parallel
//! motion and never divides by the (possibly zero) leading coefficient.

use super::{add, dot, len_sq, scale, sub};

/// Slack on the axis parameter of a cylinder hit. A hit just past an end point
/// of the segment is accepted; the neighbouring end-point sphere reports the
/// same time, so no gap opens between the pieces of a capsule.
const AXIS_PARAM_SLACK: f64 = 1e-9;

/// Relative slack of the "already touching" test (see [`within_contact`]).
const CONTACT_SLACK: f64 = 1e-9;

/// True when a point at squared distance `dist_sq` from a body counts as
/// touching a sphere of `radius` around it.
///
/// A small relative slack makes a sphere that rests on a surface (and whose
/// computed distance is `radius` up to rounding) report contact at the start of
/// the step instead of a rounding-dependent time slightly before it.
pub(super) fn within_contact(dist_sq: f64, radius: f64) -> bool {
    let limit = radius + CONTACT_SLACK * (1.0 + radius);
    dist_sq <= limit * limit
}

/// The earlier of two optional entry times.
pub(super) fn earliest(current: Option<f64>, candidate: Option<f64>) -> Option<f64> {
    match (current, candidate) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (None, other) => other,
        (kept, None) => kept,
    }
}

/// Closest point on the segment `a`-`b` to `p`, with its parameter in `[0, 1]`.
pub(super) fn closest_point_on_segment(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> ([f64; 3], f64) {
    let ab = sub(b, a);
    let ab_len_sq = len_sq(ab);
    if ab_len_sq <= 0.0 {
        return (a, 0.0);
    }
    let u = (dot(sub(p, a), ab) / ab_len_sq).clamp(0.0, 1.0);
    (add(a, scale(ab, u)), u)
}

/// First time in `[0, 1]` at which `origin + t * dir` is within `radius` of
/// `centre`.
///
/// Returns `Some(0.0)` for a point that starts inside (or on) the sphere and
/// `None` when the point never reaches it during the step.
pub(super) fn sphere_entry(
    origin: [f64; 3],
    dir: [f64; 3],
    centre: [f64; 3],
    radius: f64,
) -> Option<f64> {
    let offset = sub(origin, centre);
    let c = len_sq(offset) - radius * radius;
    if c <= 0.0 {
        return Some(0.0);
    }
    // |offset + t dir|^2 = radius^2  <=>  a t^2 + 2 b t + c = 0.
    let b = dot(offset, dir);
    if b >= 0.0 {
        // Moving away from the centre (or not at all): the distance only grows.
        return None;
    }
    let a = len_sq(dir);
    let disc = b * b - a * c;
    if disc < 0.0 {
        return None;
    }
    let t = c / (-b + disc.sqrt());
    (t <= 1.0).then_some(t)
}

/// First time in `[0, 1]` at which `origin + t * dir` reaches the lateral
/// surface of the cylinder of `radius` around the segment `a`-`b`, coming from
/// outside the cylinder, at a position along the axis that lies within the
/// segment.
///
/// A point that starts inside the infinite cylinder yields `None`: whether it
/// touches the segment is decided by the end-point spheres (see
/// [`capsule_entry`]), because it can only reach the segment's lateral surface
/// from inside after it has crossed an end cap. A zero-length segment has no
/// lateral surface and also yields `None`.
pub(super) fn cylinder_entry(
    origin: [f64; 3],
    dir: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    radius: f64,
) -> Option<f64> {
    let axis = sub(b, a);
    let axis_len_sq = len_sq(axis);
    if axis_len_sq <= 0.0 {
        return None;
    }
    let rel = sub(origin, a);
    let rel_along = dot(rel, axis);
    let dir_along = dot(dir, axis);
    // Components perpendicular to the axis, formed explicitly so that motion
    // nearly parallel to the axis does not cancel catastrophically.
    let rel_perp = sub(rel, scale(axis, rel_along / axis_len_sq));
    let dir_perp = sub(dir, scale(axis, dir_along / axis_len_sq));
    let c = len_sq(rel_perp) - radius * radius;
    if c <= 0.0 {
        return None;
    }
    let b_coef = dot(rel_perp, dir_perp);
    if b_coef >= 0.0 {
        // The distance to the axis line does not shrink.
        return None;
    }
    let a_coef = len_sq(dir_perp);
    let disc = b_coef * b_coef - a_coef * c;
    if disc < 0.0 {
        return None;
    }
    let t = c / (-b_coef + disc.sqrt());
    if t > 1.0 {
        return None;
    }
    let along = (rel_along + t * dir_along) / axis_len_sq;
    (-AXIS_PARAM_SLACK..=1.0 + AXIS_PARAM_SLACK)
        .contains(&along)
        .then_some(t)
}

/// First time in `[0, 1]` at which `origin + t * dir` is within `radius` of the
/// segment `a`-`b`: the first contact of a sphere moving along that path with a
/// capsule of that radius around the segment.
///
/// Returns `Some(0.0)` when the point starts within `radius` of the segment.
pub(super) fn capsule_entry(
    origin: [f64; 3],
    dir: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    radius: f64,
) -> Option<f64> {
    let (closest, _) = closest_point_on_segment(a, b, origin);
    if within_contact(len_sq(sub(origin, closest)), radius) {
        return Some(0.0);
    }
    let lateral = cylinder_entry(origin, dir, a, b, radius);
    let cap_a = sphere_entry(origin, dir, a, radius);
    let cap_b = sphere_entry(origin, dir, b, radius);
    earliest(earliest(lateral, cap_a), cap_b)
}
