// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact generation between a [`ShapeKind`] and an infinite plane.
//!
//! A [`ShapeKind::Plane`] `n . x = offset` is the boundary of the solid
//! half-space `n . x <= offset` (see `point_query`: points below the plane are
//! inside). A convex shape penetrates the plane when its point of least signed
//! distance `s(p) = n . p - offset` is negative, and that point is the contact
//! witness. This is evaluated in closed form from the shape's extreme points
//! (box corners, capsule end caps, convex vertices) instead of going through a
//! general support-function based algorithm, which cannot represent an
//! unbounded plane.
//!
//! Every contact follows the [`NarrowPhaseContact`] convention: the normal
//! points from shape A towards shape B, `point_a` lies on A and `point_b` on B.
//! The functions in this module take the shape as A and the plane as B, so the
//! normal is `-n` (into the plane); callers that have the plane first use
//! [`NarrowPhaseContact::flipped`].

use super::functions::{dot3, len3, scale3, sub3};
use super::types::{NarrowPhaseContact, ShapeKind};

/// Relative tolerance used to decide which extreme points are tied for the
/// deepest penetration (e.g. the four bottom corners of a box lying flat).
const TIE_TOLERANCE: f64 = 1e-9;

/// Contact between `shape` (body A) and the plane `normal . x = offset` (body B).
///
/// The plane normal need not be unit length (it is normalised together with the
/// offset); a zero or non-finite normal yields `None`, as does a plane paired
/// with another plane (two unbounded half-spaces have no meaningful contact
/// point). A shape that merely touches the plane is in contact with depth `0`,
/// matching the sphere-plane routine.
///
/// For a box lying flat on the plane, or a capsule parallel to it, several extreme
/// points tie for the deepest penetration; the witness is then their centroid
/// so a single contact does not bias the shape towards one corner. Otherwise the
/// witness is the unique deepest point.
pub(super) fn shape_vs_plane(
    shape: &ShapeKind,
    normal: [f64; 3],
    offset: f64,
) -> Option<NarrowPhaseContact> {
    let n_len = len3(normal);
    if n_len <= 1e-12 || !n_len.is_finite() {
        return None;
    }
    let n = scale3(normal, 1.0 / n_len);
    let plane_offset = offset / n_len;
    let (min_distance, deepest_point) = match shape {
        ShapeKind::Plane { .. } => return None,
        ShapeKind::Sphere { center, radius } => {
            let distance = dot3(*center, n) - plane_offset - radius;
            (distance, sub3(*center, scale3(n, *radius)))
        }
        ShapeKind::Box {
            center,
            half_extents,
        } => deepest_extreme_point(&box_corners(*center, *half_extents), n, plane_offset, 0.0)?,
        ShapeKind::Capsule { p0, p1, radius } => {
            deepest_extreme_point(&[*p0, *p1], n, plane_offset, *radius)?
        }
        ShapeKind::Convex { vertices } => deepest_extreme_point(vertices, n, plane_offset, 0.0)?,
    };
    if min_distance > 0.0 || min_distance.is_nan() {
        return None;
    }
    let depth = -min_distance;
    // `deepest_point` lies `depth` below the plane along `n`; moving by
    // `n * depth` projects it onto the plane.
    let plane_point = [
        deepest_point[0] + n[0] * depth,
        deepest_point[1] + n[1] * depth,
        deepest_point[2] + n[2] * depth,
    ];
    Some(NarrowPhaseContact {
        normal: scale3(n, -1.0),
        depth,
        point_a: deepest_point,
        point_b: plane_point,
    })
}

/// The eight corners of an axis-aligned box.
fn box_corners(center: [f64; 3], half: [f64; 3]) -> [[f64; 3]; 8] {
    let mut corners = [[0.0_f64; 3]; 8];
    for (i, corner) in corners.iter_mut().enumerate() {
        *corner = [
            center[0] + if i & 1 == 0 { -half[0] } else { half[0] },
            center[1] + if i & 2 == 0 { -half[1] } else { half[1] },
            center[2] + if i & 4 == 0 { -half[2] } else { half[2] },
        ];
    }
    corners
}

/// Deepest point of a swept point set against the plane `n . x = plane_offset`.
///
/// `points` are the extreme points of a convex shape; `inflate` is a radius
/// swept around them (capsule radius), so the shape's deepest point lies
/// `inflate` below the deepest extreme point along `n`. Returns the minimum
/// signed distance of the shape to the plane together with the deepest surface
/// point, or `None` for an empty set. When several points are within a small
/// tolerance of the minimum their centroid is returned.
fn deepest_extreme_point(
    points: &[[f64; 3]],
    n: [f64; 3],
    plane_offset: f64,
    inflate: f64,
) -> Option<(f64, [f64; 3])> {
    let distances = |p: &[f64; 3]| dot3(*p, n) - plane_offset;
    let min_s = points.iter().map(distances).fold(f64::INFINITY, f64::min);
    let max_s = points
        .iter()
        .map(distances)
        .fold(f64::NEG_INFINITY, f64::max);
    if !min_s.is_finite() || !max_s.is_finite() {
        return None;
    }
    let tolerance = TIE_TOLERANCE * (1.0 + (max_s - min_s));
    let mut sum = [0.0_f64; 3];
    let mut count = 0.0_f64;
    for p in points.iter().filter(|p| distances(p) <= min_s + tolerance) {
        sum = [sum[0] + p[0], sum[1] + p[1], sum[2] + p[2]];
        count += 1.0;
    }
    if count < 1.0 {
        return None;
    }
    let centroid = scale3(sum, 1.0 / count);
    // The centroid of the tied points has signed distance `min_s` (up to the
    // tie tolerance); sweep it by the radius to reach the shape's surface.
    let deepest = sub3(centroid, scale3(n, inflate));
    Some((min_s - inflate, deepest))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_box(center_y: f64) -> ShapeKind {
        ShapeKind::Box {
            center: [0.0, center_y, 0.0],
            half_extents: [0.5, 0.5, 0.5],
        }
    }

    #[test]
    fn box_straddling_the_plane_penetrates_by_the_lowest_corner() {
        let c = shape_vs_plane(&flat_box(0.4), [0.0, 1.0, 0.0], 0.0).expect("box straddles");
        assert!((c.depth - 0.1).abs() < 1e-12, "depth={}", c.depth);
        assert_eq!(c.normal, [0.0, -1.0, 0.0]);
        // The four bottom corners tie; the witness is their centroid.
        assert!(
            (c.point_a[1] + 0.1).abs() < 1e-12,
            "point_a={:?}",
            c.point_a
        );
        assert!(c.point_a[0].abs() < 1e-12 && c.point_a[2].abs() < 1e-12);
        assert!(c.point_b[1].abs() < 1e-12, "point_b={:?}", c.point_b);
    }

    #[test]
    fn box_above_the_plane_has_no_contact_and_touching_counts() {
        assert!(shape_vs_plane(&flat_box(0.6), [0.0, 1.0, 0.0], 0.0).is_none());
        let touching = shape_vs_plane(&flat_box(0.5), [0.0, 1.0, 0.0], 0.0)
            .expect("touching box is a contact");
        assert!(touching.depth.abs() < 1e-12);
    }

    #[test]
    fn non_unit_normal_and_offset_are_normalised_together() {
        // Plane 2y = 2 -> y = 1; a unit box centred at y = 1.3 reaches y = 0.8.
        let c = shape_vs_plane(&flat_box(1.3), [0.0, 2.0, 0.0], 2.0).expect("box straddles");
        assert!((c.depth - 0.2).abs() < 1e-12, "depth={}", c.depth);
        assert!((c.point_b[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn degenerate_inputs_yield_no_contact() {
        let plane = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        assert!(shape_vs_plane(&plane, [0.0, 1.0, 0.0], 0.0).is_none());
        assert!(shape_vs_plane(&flat_box(0.0), [0.0, 0.0, 0.0], 0.0).is_none());
        assert!(shape_vs_plane(&flat_box(0.0), [f64::NAN, 1.0, 0.0], 0.0).is_none());
        let empty = ShapeKind::Convex {
            vertices: Vec::new(),
        };
        assert!(shape_vs_plane(&empty, [0.0, 1.0, 0.0], 0.0).is_none());
    }
}
