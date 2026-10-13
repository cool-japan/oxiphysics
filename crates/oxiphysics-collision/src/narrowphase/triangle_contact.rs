// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact between a convex [`ShapeKind`] and a single triangle.
//!
//! This is the per-triangle kernel behind [`super::functions::convex_vs_mesh`].
//! The triangle is treated as the thin plate it is, never as a bounding sphere:
//!
//! * a **sphere** is resolved from the closest point of the triangle to its
//!   centre,
//! * a **capsule** from the closest pair of points between its axis segment and
//!   the triangle,
//! * every other convex shape (box, convex hull, plane) goes through the general
//!   pair router with the triangle expressed as a three-vertex convex hull, which
//!   ends in GJK + EPA.
//!
//! Contacts follow the [`NarrowPhaseContact`] convention: `normal` is the unit
//! direction from the convex shape (A) towards the triangle (B), `point_a` lies
//! on the convex shape, `point_b` on the triangle and `point_a - point_b ==
//! normal * depth`.

use super::functions::{add3, cross3, dot3, len3, scale3, shape_shape_contact, sub3};
use super::simple::{closest_point_on_triangle, segment_segment_closest};
use super::types::{NarrowPhaseContact, ShapeKind};

/// Distances below this (relative to the shape size) count as coincident, where
/// the direction from one closest point to the other is undefined.
const COINCIDENT: f64 = 1e-12;

/// Contact between `convex` and the triangle `tri`, or `None` if they do not
/// touch (or the triangle has no area).
pub(super) fn convex_vs_triangle(
    convex: &ShapeKind,
    tri: [[f64; 3]; 3],
) -> Option<NarrowPhaseContact> {
    let [v0, v1, v2] = tri;
    let face = cross3(sub3(v1, v0), sub3(v2, v0));
    let face_len = len3(face);
    let longest = len3(sub3(v1, v0))
        .max(len3(sub3(v2, v1)))
        .max(len3(sub3(v0, v2)));
    // Zero-area (or non-finite) triangles have no surface to collide with.
    if !face_len.is_finite() || face_len <= 1e-12 * longest * longest {
        return None;
    }
    let face_normal = scale3(face, 1.0 / face_len);

    let (convex_min, convex_max) = convex.aabb();
    for axis in 0..3 {
        let tri_min = v0[axis].min(v1[axis]).min(v2[axis]);
        let tri_max = v0[axis].max(v1[axis]).max(v2[axis]);
        if convex_max[axis] < tri_min || tri_max < convex_min[axis] {
            return None;
        }
    }

    match convex {
        ShapeKind::Sphere { center, radius } => {
            sphere_vs_triangle(*center, *radius, tri, face_normal)
        }
        ShapeKind::Capsule { p0, p1, radius } => {
            capsule_vs_triangle(convex, (*p0, *p1), *radius, tri, face_normal)
        }
        _ => {
            let plate = ShapeKind::Convex {
                vertices: tri.to_vec(),
            };
            shape_shape_contact(convex, &plate)
        }
    }
}

/// Sphere against the triangle `tri` (unit face normal `face_normal`).
fn sphere_vs_triangle(
    center: [f64; 3],
    radius: f64,
    tri: [[f64; 3]; 3],
    face_normal: [f64; 3],
) -> Option<NarrowPhaseContact> {
    let [v0, v1, v2] = tri;
    let closest = closest_point_on_triangle(center, v0, v1, v2);
    let offset = sub3(closest, center);
    let distance = len3(offset);
    if !distance.is_finite() || distance > radius {
        return None;
    }
    let normal = if distance > COINCIDENT * (1.0 + radius) {
        scale3(offset, 1.0 / distance)
    } else {
        // The centre lies on the triangle: leave through the front face.
        scale3(face_normal, -1.0)
    };
    Some(NarrowPhaseContact {
        normal,
        depth: radius - distance,
        point_a: add3(center, scale3(normal, radius)),
        point_b: closest,
    })
}

/// Capsule (`axis` endpoints and `radius`) against the triangle `tri`.
fn capsule_vs_triangle(
    capsule: &ShapeKind,
    axis: ([f64; 3], [f64; 3]),
    radius: f64,
    tri: [[f64; 3]; 3],
    face_normal: [f64; 3],
) -> Option<NarrowPhaseContact> {
    let (on_axis, on_triangle) = segment_triangle_closest(axis.0, axis.1, tri, face_normal);
    let offset = sub3(on_triangle, on_axis);
    let distance = len3(offset);
    if !distance.is_finite() || distance > radius {
        return None;
    }
    if distance <= COINCIDENT * (1.0 + radius) {
        // The axis pierces (or lies in) the triangle, so there is no closest-point
        // direction; the minimum translation needs the general solver.
        let plate = ShapeKind::Convex {
            vertices: tri.to_vec(),
        };
        return shape_shape_contact(capsule, &plate);
    }
    let normal = scale3(offset, 1.0 / distance);
    Some(NarrowPhaseContact {
        normal,
        depth: radius - distance,
        point_a: add3(on_axis, scale3(normal, radius)),
        point_b: on_triangle,
    })
}

/// Closest pair `(point on segment, point on triangle)` between the segment
/// `p0 p1` and the triangle `tri`; both points coincide when they intersect.
fn segment_triangle_closest(
    p0: [f64; 3],
    p1: [f64; 3],
    tri: [[f64; 3]; 3],
    face_normal: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let [v0, v1, v2] = tri;
    // A segment that crosses the triangle's plane inside the triangle intersects it.
    let d0 = dot3(face_normal, sub3(p0, v0));
    let d1 = dot3(face_normal, sub3(p1, v0));
    if d0 * d1 <= 0.0 && d0 != d1 {
        let t = d0 / (d0 - d1);
        let hit = add3(p0, scale3(sub3(p1, p0), t));
        if point_in_triangle(hit, tri, face_normal) {
            return (hit, hit);
        }
    }
    // Otherwise the closest points involve an endpoint or a triangle edge.
    let mut best = (p0, closest_point_on_triangle(p0, v0, v1, v2));
    let mut best_sq = dist_sq(best.0, best.1);
    let end = closest_point_on_triangle(p1, v0, v1, v2);
    if dist_sq(p1, end) < best_sq {
        best = (p1, end);
        best_sq = dist_sq(p1, end);
    }
    for (a, b) in [(v0, v1), (v1, v2), (v2, v0)] {
        let (_, _, on_segment, on_edge) = segment_segment_closest(p0, p1, a, b);
        let d = dist_sq(on_segment, on_edge);
        if d < best_sq {
            best = (on_segment, on_edge);
            best_sq = d;
        }
    }
    best
}

/// Whether `point` (assumed to lie in the triangle's plane) is inside `tri`.
fn point_in_triangle(point: [f64; 3], tri: [[f64; 3]; 3], face_normal: [f64; 3]) -> bool {
    let [v0, v1, v2] = tri;
    let tolerance = -1e-12 * (len3(sub3(v1, v0)) + len3(sub3(v2, v0)));
    [(v0, v1), (v1, v2), (v2, v0)]
        .iter()
        .all(|&(a, b)| dot3(face_normal, cross3(sub3(b, a), sub3(point, a))) >= tolerance)
}

fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    dot3(d, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GROUND: [[f64; 3]; 3] = [[-10.0, 0.0, -10.0], [10.0, 0.0, -10.0], [0.0, 0.0, 10.0]];

    fn sphere(center: [f64; 3], radius: f64) -> ShapeKind {
        ShapeKind::Sphere { center, radius }
    }

    fn assert_close(actual: [f64; 3], expected: [f64; 3]) {
        for (a, e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-12, "{actual:?} != {expected:?}");
        }
    }

    #[test]
    fn sphere_above_the_surface_is_not_a_contact_even_inside_the_bounding_sphere() {
        // The bounding sphere of GROUND reaches far above the plane.
        assert!(convex_vs_triangle(&sphere([0.0, 0.6, 0.0], 0.5), GROUND).is_none());
    }

    #[test]
    fn sphere_penetrating_the_face_reports_depth_and_normal_towards_the_triangle() {
        let c = convex_vs_triangle(&sphere([1.0, 0.4, 0.5], 0.5), GROUND).expect("overlap");
        assert!((c.depth - 0.1).abs() < 1e-12, "depth {}", c.depth);
        assert_close(c.normal, [0.0, -1.0, 0.0]);
        assert!((c.point_a[1] - -0.1).abs() < 1e-12, "{:?}", c.point_a);
        assert!(c.point_b[1].abs() < 1e-12 && (c.point_b[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sphere_beside_an_edge_uses_the_edge_direction() {
        let tri = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        // Beyond the x = 0 edge, 0.3 away, in the plane of the triangle.
        let c = convex_vs_triangle(&sphere([-0.3, 0.0, 1.0], 0.5), tri).expect("overlap");
        assert!((c.depth - 0.2).abs() < 1e-12);
        assert!((c.normal[0] - 1.0).abs() < 1e-12 && c.normal[1].abs() < 1e-12);
    }

    #[test]
    fn sphere_centred_on_the_triangle_leaves_through_the_front_face() {
        // GROUND winds so that its face normal is -Y, i.e. the front side is
        // below. The A -> B normal therefore points up and A moves down to leave.
        let c = convex_vs_triangle(&sphere([0.0, 0.0, 0.0], 0.5), GROUND).expect("overlap");
        assert!((c.depth - 0.5).abs() < 1e-12);
        assert_close(c.normal, [0.0, 1.0, 0.0]);
        assert!((c.point_a[1] - 0.5).abs() < 1e-12 && c.point_b[1].abs() < 1e-12);
    }

    #[test]
    fn capsule_lying_over_the_surface_touches_along_its_axis() {
        let capsule = ShapeKind::Capsule {
            p0: [-1.0, 0.4, 0.0],
            p1: [1.0, 0.4, 0.0],
            radius: 0.5,
        };
        let c = convex_vs_triangle(&capsule, GROUND).expect("overlap");
        assert!((c.depth - 0.1).abs() < 1e-12, "depth {}", c.depth);
        assert_close(c.normal, [0.0, -1.0, 0.0]);
    }

    #[test]
    fn capsule_piercing_the_triangle_resolves_through_the_general_solver() {
        let capsule = ShapeKind::Capsule {
            p0: [0.0, -1.0, 0.0],
            p1: [0.0, 1.0, 0.0],
            radius: 0.25,
        };
        let c = convex_vs_triangle(&capsule, GROUND).expect("overlap");
        // Clearing the plate along Y costs half the axis plus the radius; leaving
        // sideways would cost more than four units.
        assert!((c.depth - 1.25).abs() < 1e-6, "depth {}", c.depth);
        assert!(
            (c.normal[1].abs() - 1.0).abs() < 1e-6,
            "normal {:?}",
            c.normal
        );
    }

    #[test]
    fn capsule_beyond_the_triangle_misses() {
        let capsule = ShapeKind::Capsule {
            p0: [-1.0, 0.8, 0.0],
            p1: [1.0, 0.8, 0.0],
            radius: 0.5,
        };
        assert!(convex_vs_triangle(&capsule, GROUND).is_none());
    }

    #[test]
    fn box_resting_on_the_surface_penetrates_by_its_lowest_face() {
        let cube = ShapeKind::Box {
            center: [0.5, 0.4, 0.0],
            half_extents: [0.5, 0.5, 0.5],
        };
        let c = convex_vs_triangle(&cube, GROUND).expect("overlap");
        assert!((c.depth - 0.1).abs() < 1e-9, "depth {}", c.depth);
        assert!(c.normal[1] < -0.999, "normal {:?}", c.normal);
        let lifted = ShapeKind::Box {
            center: [0.5, 0.6, 0.0],
            half_extents: [0.5, 0.5, 0.5],
        };
        assert!(convex_vs_triangle(&lifted, GROUND).is_none());
    }

    #[test]
    fn degenerate_triangles_never_collide() {
        let sliver = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(convex_vs_triangle(&sphere([1.0, 0.0, 0.0], 1.0), sliver).is_none());
    }
}
