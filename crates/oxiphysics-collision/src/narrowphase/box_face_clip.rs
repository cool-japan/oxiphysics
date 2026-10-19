// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact manifolds for two oriented boxes by reference/incident face clipping.
//!
//! [`build_box_manifold`] turns the result of a separating-axis test (a contact
//! normal and a penetration depth) into the contact points of the overlap
//! region:
//!
//! * **Face contacts** (the normal is, or is close to, a face axis of one of the
//!   boxes). The box with the face best aligned to the normal supplies the
//!   *reference face*; the face of the other box that is most anti-parallel to
//!   it is the *incident face*. The incident quad is clipped, in 3-D, against
//!   the four side planes of the reference face (Sutherland-Hodgman) and every
//!   clipped vertex that lies on or below the reference plane becomes a contact
//!   whose depth is its own distance below that plane. A box lying flat gives
//!   the four corners of the contact patch, a box standing on an edge two
//!   points and a box standing on a vertex one point, each with the depth the
//!   geometry dictates.
//! * **Edge-edge contacts** (the normal is perpendicular to an edge of each box
//!   and parallel to no face): a single contact at the midpoint of the closest
//!   points of the two supporting edges, with the depth reported by the
//!   separating-axis test.
//!
//! The normal follows the convention of [`ContactPoint`]: it points from body B
//! towards body A, so A is pushed out along `+normal` and B along `-normal`.

use super::functions::{add3, closest_segment_segment, dot3, len3, scale3, sub3};
use super::manifold::{ContactManifold, ContactPoint, clip_polygon_by_halfspace};

type Vec3 = [f64; 3];

/// Tolerance, relative to the largest box half-extent, for accepting a clipped
/// vertex that sits on the reference plane (touching) and for merging
/// coincident contact points.
const SIZE_TOLERANCE: f64 = 1e-9;
/// Absolute cosine below which a box axis counts as perpendicular to the normal.
const PERPENDICULAR: f64 = 1e-6;
/// Absolute cosine above which a box axis counts as parallel to the normal.
const PARALLEL: f64 = 1.0 - 1e-6;

/// `v / |v|`, or `None` for a zero or non-finite vector.
fn unit(v: Vec3) -> Option<Vec3> {
    let length = len3(v);
    (length.is_finite() && length > 1e-12).then(|| scale3(v, 1.0 / length))
}

/// A box in world space: centre, orthonormal axes and half-extents.
struct OrientedBox {
    centre: Vec3,
    axes: [Vec3; 3],
    half: [f64; 3],
}

impl OrientedBox {
    /// Read a box from its half-extents and a row-major homogeneous transform
    /// (rotation in the upper-left 3x3 block, translation in the last column).
    /// Returns `None` for negative or non-finite extents and for a transform
    /// with a degenerate rotation.
    fn new(half: [f64; 3], transform: [[f64; 4]; 4]) -> Option<Self> {
        if !half.iter().all(|h| h.is_finite() && *h >= 0.0) {
            return None;
        }
        let centre = [transform[0][3], transform[1][3], transform[2][3]];
        if !centre.iter().all(|c| c.is_finite()) {
            return None;
        }
        let column = |c: usize| [transform[0][c], transform[1][c], transform[2][c]];
        Some(Self {
            centre,
            axes: [unit(column(0))?, unit(column(1))?, unit(column(2))?],
            half,
        })
    }

    /// The axis index whose direction is most aligned with `direction`, and the
    /// signed cosine between that axis and `direction`.
    fn best_axis(&self, direction: Vec3) -> (usize, f64) {
        let mut best = (0, dot3(self.axes[0], direction));
        for (index, axis) in self.axes.iter().enumerate().skip(1) {
            let cosine = dot3(*axis, direction);
            if cosine.abs() > best.1.abs() {
                best = (index, cosine);
            }
        }
        best
    }

    /// The four corners, in loop order, of the face whose outward normal is
    /// `sign * axes[axis]`.
    fn face_corners(&self, axis: usize, sign: f64) -> [Vec3; 4] {
        let centre = add3(self.centre, scale3(self.axes[axis], sign * self.half[axis]));
        let u = (axis + 1) % 3;
        let w = (axis + 2) % 3;
        let along_u = scale3(self.axes[u], self.half[u]);
        let along_w = scale3(self.axes[w], self.half[w]);
        [
            add3(add3(centre, along_u), along_w),
            add3(sub3(centre, along_u), along_w),
            sub3(sub3(centre, along_u), along_w),
            add3(sub3(centre, along_w), along_u),
        ]
    }

    /// The edge parallel to `axes[free]` that is furthest along `direction`.
    fn support_edge(&self, free: usize, direction: Vec3) -> (Vec3, Vec3) {
        let mut base = self.centre;
        for (index, axis) in self.axes.iter().enumerate() {
            if index == free {
                continue;
            }
            let sign = if dot3(*axis, direction) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            base = add3(base, scale3(*axis, sign * self.half[index]));
        }
        let along = scale3(self.axes[free], self.half[free]);
        (sub3(base, along), add3(base, along))
    }
}

/// Clip `polygon` to the prism over the face of `reference` that is
/// perpendicular to `reference.axes[axis]`: the four planes through the sides
/// of that face, parallel to `axis`.
fn clip_to_face_prism(mut polygon: Vec<Vec3>, reference: &OrientedBox, axis: usize) -> Vec<Vec3> {
    for tangent in [(axis + 1) % 3, (axis + 2) % 3] {
        for sign in [1.0, -1.0] {
            if polygon.is_empty() {
                return polygon;
            }
            let side = scale3(reference.axes[tangent], sign);
            // Keep `side . p <= side . centre + half`, i.e. `-side . p >= -limit`.
            let limit = dot3(side, reference.centre) + reference.half[tangent];
            polygon = clip_polygon_by_halfspace(&polygon, scale3(side, -1.0), -limit);
        }
    }
    polygon
}

/// Contact of two crossing edges: the midpoint of their closest points.
fn edge_edge_contact(
    box_a: &OrientedBox,
    box_b: &OrientedBox,
    normal: Vec3,
    edge_axes: (usize, usize),
    depth: f64,
) -> ContactPoint {
    // A reaches towards B along `-normal`, B towards A along `+normal`.
    let (a0, a1) = box_a.support_edge(edge_axes.0, scale3(normal, -1.0));
    let (b0, b1) = box_b.support_edge(edge_axes.1, normal);
    let (on_a, on_b) = closest_segment_segment(a0, a1, b0, b1);
    ContactPoint {
        point: scale3(add3(on_a, on_b), 0.5),
        normal,
        depth,
    }
}

/// Build the contact manifold of two overlapping oriented boxes.
///
/// * `half_extents_a` / `half_extents_b`: half-extents along each box's local
///   axes.
/// * `transform_a` / `transform_b`: **row-major** 4x4 homogeneous transforms
///   (`transform[row][col]`) with the rotation in the upper-left 3x3 block and
///   the translation in the last column. The columns of the rotation are the
///   box's local axes in world space.
/// * `contact_normal`: the contact normal from the separating-axis test, pointing
///   from box B towards box A (it is normalised internally).
/// * `penetration_depth`: the overlap of the boxes along `contact_normal`. It is
///   the depth of an edge-edge contact, where no face plane measures it; for
///   face contacts every point carries its own geometric depth instead. A
///   negative or non-finite value (the boxes are separated) yields an empty
///   manifold.
///
/// The manifold holds up to eight points (the clipped incident face is a convex
/// polygon of at most eight vertices); reduce it with
/// [`reduce_manifold`](super::manifold::reduce_manifold) or
/// [`build_4point_manifold`](super::manifold::build_4point_manifold) if the
/// solver wants at most four. Each `ContactPoint::point` is the midpoint between
/// the penetrating vertex and its projection on the reference face, so the
/// contacts straddle the interface symmetrically.
///
/// The manifold is empty for invalid input (a zero or non-finite normal, a
/// degenerate transform, negative extents) and when no vertex of the incident
/// face reaches the reference face.
pub fn build_box_manifold(
    half_extents_a: [f64; 3],
    transform_a: [[f64; 4]; 4],
    half_extents_b: [f64; 3],
    transform_b: [[f64; 4]; 4],
    contact_normal: [f64; 3],
    penetration_depth: f64,
) -> ContactManifold {
    let Some(normal) = unit(contact_normal) else {
        return ContactManifold::new(contact_normal);
    };
    let mut manifold = ContactManifold::new(normal);
    if !penetration_depth.is_finite() || penetration_depth < 0.0 {
        return manifold;
    }
    let (Some(box_a), Some(box_b)) = (
        OrientedBox::new(half_extents_a, transform_a),
        OrientedBox::new(half_extents_b, transform_b),
    ) else {
        return manifold;
    };

    let (axis_a, cosine_a) = box_a.best_axis(normal);
    let (axis_b, cosine_b) = box_b.best_axis(normal);
    let face_aligned = cosine_a.abs() >= PARALLEL || cosine_b.abs() >= PARALLEL;
    if !face_aligned {
        let perpendicular = |b: &OrientedBox| {
            b.axes
                .iter()
                .position(|axis| dot3(*axis, normal).abs() <= PERPENDICULAR)
        };
        if let (Some(edge_a), Some(edge_b)) = (perpendicular(&box_a), perpendicular(&box_b)) {
            manifold.contacts.push(edge_edge_contact(
                &box_a,
                &box_b,
                normal,
                (edge_a, edge_b),
                penetration_depth,
            ));
            return manifold;
        }
    }

    // The box whose face axis is closer to the normal owns the reference face;
    // A wins a tie so that the result does not depend on rounding noise.
    let reference_is_a = cosine_a.abs() >= cosine_b.abs();
    let (reference, incident, reference_axis, towards_incident) = if reference_is_a {
        (&box_a, &box_b, axis_a, scale3(normal, -1.0))
    } else {
        (&box_b, &box_a, axis_b, normal)
    };
    // The reference face is the one that faces the incident box.
    let facing = if dot3(reference.axes[reference_axis], towards_incident) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let reference_normal = scale3(reference.axes[reference_axis], facing);
    let reference_offset =
        dot3(reference_normal, reference.centre) + reference.half[reference_axis];

    // The incident face is the one most anti-parallel to the reference face.
    let (incident_axis, incident_cosine) = incident.best_axis(reference_normal);
    let incident_sign = if incident_cosine >= 0.0 { -1.0 } else { 1.0 };
    let incident_face = incident.face_corners(incident_axis, incident_sign).to_vec();
    let clipped = clip_to_face_prism(incident_face, reference, reference_axis);

    let size = box_a
        .half
        .iter()
        .chain(box_b.half.iter())
        .fold(0.0_f64, |largest, half| largest.max(*half));
    let tolerance = SIZE_TOLERANCE * size;
    for vertex in clipped {
        // Signed distance of the vertex from the reference plane; negative
        // below it, i.e. inside the reference box.
        let separation = dot3(reference_normal, vertex) - reference_offset;
        if separation > tolerance {
            continue;
        }
        let on_reference = sub3(vertex, scale3(reference_normal, separation));
        let point = scale3(add3(vertex, on_reference), 0.5);
        let duplicate = manifold
            .contacts
            .iter()
            .any(|contact| len3(sub3(contact.point, point)) <= tolerance);
        if !duplicate {
            manifold.contacts.push(ContactPoint {
                point,
                normal,
                depth: (-separation).max(0.0),
            });
        }
    }
    manifold
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::narrowphase::functions::cross3;

    type Matrix3 = [[f64; 3]; 3];

    const IDENTITY: Matrix3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    const UP: Vec3 = [0.0, 1.0, 0.0];
    const DOWN: Vec3 = [0.0, -1.0, 0.0];
    const UNIT_CUBE: Vec3 = [0.5, 0.5, 0.5];

    fn transform(rotation: Matrix3, centre: Vec3) -> [[f64; 4]; 4] {
        [
            [rotation[0][0], rotation[0][1], rotation[0][2], centre[0]],
            [rotation[1][0], rotation[1][1], rotation[1][2], centre[1]],
            [rotation[2][0], rotation[2][1], rotation[2][2], centre[2]],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }

    fn rotation_x(angle: f64) -> Matrix3 {
        let (s, c) = angle.sin_cos();
        [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]
    }

    fn rotation_y(angle: f64) -> Matrix3 {
        let (s, c) = angle.sin_cos();
        [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]
    }

    fn rotation_z(angle: f64) -> Matrix3 {
        let (s, c) = angle.sin_cos();
        [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
    }

    /// Whether every expected point is produced by exactly one contact.
    fn assert_points(manifold: &ContactManifold, expected: &[Vec3]) {
        assert_eq!(
            manifold.contacts.len(),
            expected.len(),
            "contacts: {:?}",
            manifold.contacts
        );
        for wanted in expected {
            let matches = manifold
                .contacts
                .iter()
                .filter(|c| len3(sub3(c.point, *wanted)) < 1e-9)
                .count();
            assert_eq!(matches, 1, "{wanted:?} in {:?}", manifold.contacts);
        }
    }

    #[test]
    fn stacked_cubes_overlapping_by_a_tenth_give_four_points_of_that_depth() {
        // A rests on B and sinks 0.1 into it; the normal points from B to A.
        let manifold = build_box_manifold(
            UNIT_CUBE,
            transform(IDENTITY, [0.0, 0.9, 0.0]),
            UNIT_CUBE,
            transform(IDENTITY, [0.0; 3]),
            UP,
            0.1,
        );
        assert_eq!(manifold.normal, UP);
        assert_points(
            &manifold,
            &[
                [0.5, 0.45, 0.5],
                [-0.5, 0.45, 0.5],
                [-0.5, 0.45, -0.5],
                [0.5, 0.45, -0.5],
            ],
        );
        for contact in &manifold.contacts {
            assert!(
                (contact.depth - 0.1).abs() < 1e-12,
                "depth {}",
                contact.depth
            );
            assert_eq!(contact.normal, UP);
        }
    }

    #[test]
    fn swapping_the_boxes_gives_the_same_points_with_the_opposite_normal() {
        let swapped = build_box_manifold(
            UNIT_CUBE,
            transform(IDENTITY, [0.0; 3]),
            UNIT_CUBE,
            transform(IDENTITY, [0.0, 0.9, 0.0]),
            DOWN,
            0.1,
        );
        assert_eq!(swapped.normal, DOWN);
        assert_points(
            &swapped,
            &[
                [0.5, 0.45, 0.5],
                [-0.5, 0.45, 0.5],
                [-0.5, 0.45, -0.5],
                [0.5, 0.45, -0.5],
            ],
        );
        for contact in &swapped.contacts {
            assert!((contact.depth - 0.1).abs() < 1e-12);
            assert_eq!(contact.normal, DOWN);
        }
    }

    #[test]
    fn shifted_cubes_clip_to_the_common_square() {
        let manifold = build_box_manifold(
            UNIT_CUBE,
            transform(IDENTITY, [0.25, 0.9, 0.25]),
            UNIT_CUBE,
            transform(IDENTITY, [0.0; 3]),
            UP,
            0.1,
        );
        // The footprints share x and z in [-0.25, 0.5].
        assert_points(
            &manifold,
            &[
                [0.5, 0.45, 0.5],
                [-0.25, 0.45, 0.5],
                [-0.25, 0.45, -0.25],
                [0.5, 0.45, -0.25],
            ],
        );
        for contact in &manifold.contacts {
            assert!((contact.depth - 0.1).abs() < 1e-12);
        }
    }

    #[test]
    fn a_cube_turned_by_45_degrees_about_the_vertical_gives_an_octagon() {
        let manifold = build_box_manifold(
            UNIT_CUBE,
            transform(rotation_y(std::f64::consts::FRAC_PI_4), [0.0, 0.9, 0.0]),
            UNIT_CUBE,
            transform(IDENTITY, [0.0; 3]),
            UP,
            0.1,
        );
        assert_eq!(manifold.contacts.len(), 8, "{:?}", manifold.contacts);
        let reach = 0.5 * std::f64::consts::SQRT_2;
        for contact in &manifold.contacts {
            let [x, y, z] = contact.point;
            assert!((y - 0.45).abs() < 1e-9, "height {y}");
            assert!((contact.depth - 0.1).abs() < 1e-9);
            // Inside the square footprint of B and the diamond footprint of A,
            // and on the boundary of at least one of them.
            assert!(x.abs() <= 0.5 + 1e-9 && z.abs() <= 0.5 + 1e-9, "{x} {z}");
            assert!(x.abs() + z.abs() <= reach + 1e-9, "{x} {z}");
            let on_square = (x.abs() - 0.5).abs() < 1e-9 || (z.abs() - 0.5).abs() < 1e-9;
            let on_diamond = (x.abs() + z.abs() - reach).abs() < 1e-9;
            assert!(on_square || on_diamond, "{x} {z}");
        }
    }

    #[test]
    fn a_tilted_cube_has_a_depth_for_each_corner() {
        // A narrow cube tilted by 5 degrees about z sinks into a wide slab; its
        // two lower corner pairs penetrate by different amounts.
        let tilt = 5.0_f64.to_radians();
        let half_a = [0.4, 0.5, 0.4];
        let rotation = rotation_z(tilt);
        let centre = [0.0, 0.9, 0.0];
        let manifold = build_box_manifold(
            half_a,
            transform(rotation, centre),
            [1.0, 0.5, 1.0],
            transform(IDENTITY, [0.0; 3]),
            UP,
            0.15,
        );
        assert_eq!(manifold.contacts.len(), 4, "{:?}", manifold.contacts);
        let mut expected = Vec::new();
        for sx in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                let local = [sx * half_a[0], -half_a[1], sz * half_a[2]];
                let height = centre[1] + rotation[1][0] * local[0] + rotation[1][1] * local[1];
                expected.push(0.5 - height);
            }
        }
        let mut depths: Vec<f64> = manifold.contacts.iter().map(|c| c.depth).collect();
        depths.sort_by(f64::total_cmp);
        expected.sort_by(f64::total_cmp);
        for (found, wanted) in depths.iter().zip(&expected) {
            assert!((found - wanted).abs() < 1e-12, "{depths:?} vs {expected:?}");
        }
        assert!(
            depths[3] - depths[0] > 0.05,
            "depths should differ: {depths:?}"
        );
        for contact in &manifold.contacts {
            // The contact lies halfway between the corner and the slab's top.
            let height = 0.5 - contact.depth * 0.5;
            assert!((contact.point[1] - height).abs() < 1e-12);
        }
    }

    #[test]
    fn a_cube_standing_on_a_vertex_touches_in_one_point() {
        // Orient the cube so that a body diagonal is vertical.
        let third = 1.0 / 3.0_f64.sqrt();
        let spread = (2.0_f64 / 3.0).sqrt();
        let axis = |phi: f64| [spread * phi.cos(), third, spread * phi.sin()];
        let a1 = axis(0.0);
        let a2 = axis(2.0 * std::f64::consts::FRAC_PI_3);
        let a3 = axis(4.0 * std::f64::consts::FRAC_PI_3);
        let rotation = [
            [a1[0], a2[0], a3[0]],
            [a1[1], a2[1], a3[1]],
            [a1[2], a2[2], a3[2]],
        ];
        // The lowest vertex is half a body diagonal below the centre; sink it by
        // 0.05 into the top face of B.
        let drop = 3.0_f64.sqrt() * 0.5;
        let manifold = build_box_manifold(
            UNIT_CUBE,
            transform(rotation, [0.0, 0.45 + drop, 0.0]),
            UNIT_CUBE,
            transform(IDENTITY, [0.0; 3]),
            UP,
            0.05,
        );
        assert_points(&manifold, &[[0.0, 0.475, 0.0]]);
        assert!((manifold.contacts[0].depth - 0.05).abs() < 1e-12);
    }

    #[test]
    fn a_cube_resting_on_an_edge_touches_in_two_points() {
        // Turned by 45 degrees about x, the cube stands on an edge along x.
        let half_a = [0.4, 0.5, 0.5];
        let drop = 0.5 * std::f64::consts::SQRT_2;
        let manifold = build_box_manifold(
            half_a,
            transform(
                rotation_x(std::f64::consts::FRAC_PI_4),
                [0.0, 0.45 + drop, 0.0],
            ),
            UNIT_CUBE,
            transform(IDENTITY, [0.0; 3]),
            UP,
            0.05,
        );
        assert_points(&manifold, &[[0.4, 0.475, 0.0], [-0.4, 0.475, 0.0]]);
        for contact in &manifold.contacts {
            assert!((contact.depth - 0.05).abs() < 1e-12);
        }
    }

    #[test]
    fn crossing_edges_give_one_contact_with_the_reported_depth() {
        // A is turned 45 degrees about x, B 45 degrees about z; A's lowest edge
        // (along x) crosses B's highest edge (along z) in a plane normal to y.
        let diagonal = 0.5 * std::f64::consts::SQRT_2;
        let overlap = 2.0 * diagonal - 1.4;
        let manifold = build_box_manifold(
            UNIT_CUBE,
            transform(rotation_x(std::f64::consts::FRAC_PI_4), [0.0, 1.4, 0.0]),
            UNIT_CUBE,
            transform(rotation_z(std::f64::consts::FRAC_PI_4), [0.0; 3]),
            UP,
            overlap,
        );
        assert_points(&manifold, &[[0.0, 0.7, 0.0]]);
        assert_eq!(manifold.contacts[0].depth, overlap);
        assert_eq!(manifold.contacts[0].normal, UP);
    }

    #[test]
    fn separated_touching_and_invalid_inputs_are_handled() {
        let b = transform(IDENTITY, [0.0; 3]);
        // A hovering 0.1 above B: nothing reaches the reference face.
        let hovering = build_box_manifold(
            UNIT_CUBE,
            transform(IDENTITY, [0.0, 1.1, 0.0]),
            UNIT_CUBE,
            b,
            UP,
            0.0,
        );
        assert!(hovering.contacts.is_empty(), "{:?}", hovering.contacts);
        // Exactly touching: four contacts of depth zero.
        let touching = build_box_manifold(
            UNIT_CUBE,
            transform(IDENTITY, [0.0, 1.0, 0.0]),
            UNIT_CUBE,
            b,
            UP,
            0.0,
        );
        assert_eq!(touching.contacts.len(), 4);
        assert!(touching.contacts.iter().all(|c| c.depth == 0.0));
        let a = transform(IDENTITY, [0.0, 0.9, 0.0]);
        // A negative or non-finite depth means the SAT found a gap.
        for depth in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(
                build_box_manifold(UNIT_CUBE, a, UNIT_CUBE, b, UP, depth)
                    .contacts
                    .is_empty()
            );
        }
        // A zero or non-finite normal, negative extents and a degenerate
        // rotation cannot be clipped.
        for bad_normal in [[0.0; 3], [f64::NAN, 1.0, 0.0]] {
            assert!(
                build_box_manifold(UNIT_CUBE, a, UNIT_CUBE, b, bad_normal, 0.1)
                    .contacts
                    .is_empty()
            );
        }
        assert!(
            build_box_manifold([-0.5, 0.5, 0.5], a, UNIT_CUBE, b, UP, 0.1)
                .contacts
                .is_empty()
        );
        let flat = transform([[0.0; 3]; 3], [0.0, 0.9, 0.0]);
        assert!(
            build_box_manifold(UNIT_CUBE, flat, UNIT_CUBE, b, UP, 0.1)
                .contacts
                .is_empty()
        );
    }

    /// Deterministic generator for the randomised footprint check.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }

        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.next_f64()
        }
    }

    /// Whether `point` lies inside the footprint of a box of half-extents `half`
    /// turned by `yaw` about the vertical axis at `centre`.
    fn in_footprint(point: Vec3, centre: Vec3, yaw: f64, half: Vec3) -> bool {
        let (s, c) = yaw.sin_cos();
        let dx = point[0] - centre[0];
        let dz = point[2] - centre[2];
        let local_x = c * dx - s * dz;
        let local_z = s * dx + c * dz;
        local_x.abs() <= half[0] + 1e-9 && local_z.abs() <= half[2] + 1e-9
    }

    #[test]
    fn random_yawed_stacks_lie_in_both_footprints_at_the_overlap_depth() {
        let mut rng = Lcg(0xB0C5_F0C5);
        let mut with_patch = 0;
        for _ in 0..200 {
            let half_a = [
                rng.range(0.3, 0.8),
                rng.range(0.3, 0.6),
                rng.range(0.3, 0.8),
            ];
            let half_b = [
                rng.range(0.3, 0.8),
                rng.range(0.3, 0.6),
                rng.range(0.3, 0.8),
            ];
            let yaw_a = rng.range(0.0, 3.0);
            let yaw_b = rng.range(0.0, 3.0);
            let offset = [rng.range(-0.3, 0.3), rng.range(-0.3, 0.3)];
            let overlap = rng.range(0.01, 0.25);
            // B sits at the origin; A is above it and sinks `overlap` into B.
            let centre_b = [0.0, 0.0, 0.0];
            let centre_a = [offset[0], half_a[1] + half_b[1] - overlap, offset[1]];
            let manifold = build_box_manifold(
                half_a,
                transform(rotation_y(yaw_a), centre_a),
                half_b,
                transform(rotation_y(yaw_b), centre_b),
                UP,
                overlap,
            );
            if manifold.contacts.len() >= 3 {
                with_patch += 1;
            }
            for contact in &manifold.contacts {
                assert!(
                    (contact.depth - overlap).abs() < 1e-9,
                    "depth {} vs overlap {overlap}",
                    contact.depth
                );
                assert!(in_footprint(contact.point, centre_a, yaw_a, half_a));
                assert!(in_footprint(contact.point, centre_b, yaw_b, half_b));
                assert_eq!(contact.normal, UP);
            }
        }
        assert!(with_patch > 150, "only {with_patch} patches");
    }

    #[test]
    fn clipping_keeps_the_polygon_inside_the_reference_prism() {
        // A tall incident quad is cut down to the 1 x 1 reference face.
        let reference = OrientedBox::new(UNIT_CUBE, transform(IDENTITY, [0.0; 3]))
            .expect("valid reference box");
        let quad = vec![
            [-2.0, 0.0, -2.0],
            [2.0, 0.0, -2.0],
            [2.0, 0.0, 2.0],
            [-2.0, 0.0, 2.0],
        ];
        let clipped = clip_to_face_prism(quad, &reference, 1);
        assert_eq!(clipped.len(), 4);
        for vertex in &clipped {
            assert!(vertex[0].abs() <= 0.5 + 1e-12 && vertex[2].abs() <= 0.5 + 1e-12);
        }
        // A quad that misses the prism entirely clips to nothing.
        let far = vec![
            [3.0, 0.0, 3.0],
            [4.0, 0.0, 3.0],
            [4.0, 0.0, 4.0],
            [3.0, 0.0, 4.0],
        ];
        assert!(clip_to_face_prism(far, &reference, 1).is_empty());
    }

    /// Rotation matrix of the quaternion `(w, x, y, z)` (normalised first).
    fn rotation_from_quaternion(q: [f64; 4]) -> Matrix3 {
        let length = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
        let [w, x, y, z] = q.map(|component| component / length);
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    }

    /// Whether `point` lies inside `b`, up to a rounding tolerance.
    fn inside_box(b: &OrientedBox, point: Vec3) -> bool {
        let offset = sub3(point, b.centre);
        (0..3).all(|i| dot3(b.axes[i], offset).abs() <= b.half[i] + 1e-9)
    }

    /// The kind of axis a separating-axis test ended on.
    #[derive(Clone, Copy, Debug, PartialEq)]
    enum AxisKind {
        /// A face normal of box A, with the local axis index.
        FaceOfA(usize),
        /// A face normal of box B, with the local axis index.
        FaceOfB(usize),
        /// The cross product of an edge direction of each box.
        EdgePair,
    }

    /// Axis of least overlap of two oriented boxes, oriented from B towards A,
    /// with the overlap along it; `None` if the boxes are separated.
    fn least_overlap(box_a: &OrientedBox, box_b: &OrientedBox) -> Option<(Vec3, f64, AxisKind)> {
        let mut candidates: Vec<(Vec3, AxisKind)> = Vec::new();
        for (index, axis) in box_a.axes.iter().enumerate() {
            candidates.push((*axis, AxisKind::FaceOfA(index)));
        }
        for (index, axis) in box_b.axes.iter().enumerate() {
            candidates.push((*axis, AxisKind::FaceOfB(index)));
        }
        for axis_a in &box_a.axes {
            for axis_b in &box_b.axes {
                let crossed = cross3(*axis_a, *axis_b);
                let length = len3(crossed);
                // Nearly parallel edges have no well-defined crossing axis.
                if length > 1e-3 {
                    candidates.push((scale3(crossed, 1.0 / length), AxisKind::EdgePair));
                }
            }
        }
        let offset = sub3(box_a.centre, box_b.centre);
        let mut best: Option<(Vec3, f64, AxisKind)> = None;
        for (axis, kind) in candidates {
            let radius = |b: &OrientedBox| -> f64 {
                (0..3)
                    .map(|i| b.half[i] * dot3(b.axes[i], axis).abs())
                    .sum()
            };
            let along = dot3(offset, axis);
            let overlap = radius(box_a) + radius(box_b) - along.abs();
            if overlap < 0.0 {
                return None;
            }
            // Face axes come first, so an edge axis has to beat them clearly.
            if best.is_none_or(|(_, least, _)| overlap < least - 1e-9) {
                let towards_a = if along >= 0.0 {
                    axis
                } else {
                    scale3(axis, -1.0)
                };
                best = Some((towards_a, overlap, kind));
            }
        }
        best
    }

    #[test]
    fn random_oriented_boxes_give_contacts_on_the_reference_face_within_the_overlap() {
        let mut rng = Lcg(0x0B0C_5F0C_0183);
        let (mut face_cases, mut with_contacts, mut edge_cases) = (0, 0, 0);
        for _ in 0..900 {
            let mut half = || {
                [
                    rng.range(0.3, 1.0),
                    rng.range(0.3, 1.0),
                    rng.range(0.3, 1.0),
                ]
            };
            let half_a = half();
            let half_b = half();
            let mut orientation = || {
                rotation_from_quaternion([
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                ])
            };
            let rotation_a = orientation();
            let rotation_b = orientation();
            let centre_a = [
                rng.range(-1.2, 1.2),
                rng.range(-1.2, 1.2),
                rng.range(-1.2, 1.2),
            ];
            let transform_a = transform(rotation_a, centre_a);
            let transform_b = transform(rotation_b, [0.0; 3]);
            let (Some(box_a), Some(box_b)) = (
                OrientedBox::new(half_a, transform_a),
                OrientedBox::new(half_b, transform_b),
            ) else {
                panic!("rotations from quaternions are valid");
            };
            let Some((normal, overlap, kind)) = least_overlap(&box_a, &box_b) else {
                continue;
            };
            let manifold =
                build_box_manifold(half_a, transform_a, half_b, transform_b, normal, overlap);
            assert!(manifold.contacts.len() <= 8, "{:?}", manifold.contacts);
            for contact in &manifold.contacts {
                assert!(len3(sub3(contact.normal, normal)) < 1e-12);
                assert!(contact.depth >= 0.0 && contact.depth <= overlap + 1e-9);
            }
            for (i, first) in manifold.contacts.iter().enumerate() {
                for second in &manifold.contacts[i + 1..] {
                    assert!(len3(sub3(first.point, second.point)) > 1e-12, "duplicate");
                }
            }
            let (reference, axis, outward) = match kind {
                AxisKind::EdgePair => {
                    edge_cases += 1;
                    assert_eq!(manifold.contacts.len(), 1, "{:?}", manifold.contacts);
                    assert_eq!(manifold.contacts[0].depth, overlap);
                    // Midway between the two crossing edges, i.e. inside both boxes.
                    let point = manifold.contacts[0].point;
                    assert!(inside_box(&box_a, point), "{point:?} is outside box A");
                    assert!(inside_box(&box_b, point), "{point:?} is outside box B");
                    continue;
                }
                // The normal points from B to A, so A's reference face looks
                // along its negative and B's along its positive direction.
                AxisKind::FaceOfA(axis) => (&box_a, axis, scale3(normal, -1.0)),
                AxisKind::FaceOfB(axis) => (&box_b, axis, normal),
            };
            face_cases += 1;
            if !manifold.contacts.is_empty() {
                with_contacts += 1;
            }
            let facing = if dot3(reference.axes[axis], outward) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let reference_normal = scale3(reference.axes[axis], facing);
            let plane = dot3(reference_normal, reference.centre) + reference.half[axis];
            for contact in &manifold.contacts {
                // Halfway between the penetrating vertex and the reference face.
                let height = dot3(reference_normal, contact.point) - plane;
                assert!(
                    (height + 0.5 * contact.depth).abs() < 1e-9,
                    "height {height} for depth {}",
                    contact.depth
                );
                for tangent in [(axis + 1) % 3, (axis + 2) % 3] {
                    let lateral = dot3(
                        reference.axes[tangent],
                        sub3(contact.point, reference.centre),
                    );
                    assert!(lateral.abs() <= reference.half[tangent] + 1e-9, "{lateral}");
                }
            }
        }
        assert!(face_cases > 100, "only {face_cases} face cases");
        assert!(edge_cases > 20, "only {edge_cases} edge cases");
        assert!(
            4 * with_contacts >= 3 * face_cases,
            "{with_contacts} of {face_cases} face cases produced contacts"
        );
    }
}
