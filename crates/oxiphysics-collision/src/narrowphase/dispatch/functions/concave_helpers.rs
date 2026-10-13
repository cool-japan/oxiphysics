// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Helpers of [`super::dispatch_concave_vs_convex`]: the mesh-frame bounds used
//! to cull convex pieces and the spatial reduction of the collected contacts.

use crate::types::Contact;
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::Shape;

/// Axis-aligned bounds of `convex`, expressed in the local frame of the mesh.
///
/// The shape's local bounding box is carried through `convex_transform` into the
/// world and then through the inverse of `mesh_transform`; the bounds of the
/// eight transformed corners enclose the shape for any rotation. If a corner
/// does not come out finite (an unbounded shape such as a plane), the bounds
/// are infinite so that no piece is culled.
pub(super) fn convex_bounds_in_mesh_frame(
    convex: &dyn Shape,
    convex_transform: &Transform,
    mesh_transform: &Transform,
) -> ([f64; 3], [f64; 3]) {
    let local = convex.bounding_box();
    let to_mesh = mesh_transform.inverse().compose(convex_transform);
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for corner in 0..8_u32 {
        let pick = |low: f64, high: f64, bit: u32| if corner & bit == 0 { low } else { high };
        let point = to_mesh.transform_point(&Vec3::new(
            pick(local.min.x, local.max.x, 1),
            pick(local.min.y, local.max.y, 2),
            pick(local.min.z, local.max.z, 4),
        ));
        let coordinates = [point.x, point.y, point.z];
        if coordinates.iter().any(|c| !c.is_finite()) {
            return ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]);
        }
        for (axis, value) in coordinates.iter().enumerate() {
            min[axis] = min[axis].min(*value);
            max[axis] = max[axis].max(*value);
        }
    }
    (min, max)
}

/// `true` when the boxes `[a_min, a_max]` and `[b_min, b_max]` overlap or touch.
pub(super) fn bounds_overlap(
    a_min: [f64; 3],
    a_max: [f64; 3],
    b_min: [f64; 3],
    b_max: [f64; 3],
) -> bool {
    (0..3).all(|axis| a_min[axis] <= b_max[axis] && a_max[axis] >= b_min[axis])
}

/// Squared distance from `points[candidate]` to the nearest chosen point.
fn distance_to_set_sq(points: &[Vec3], chosen: &[usize], candidate: usize) -> f64 {
    chosen
        .iter()
        .map(|&c| (points[candidate] - points[c]).norm_squared())
        .fold(f64::INFINITY, f64::min)
}

/// Unchosen point furthest from the chosen set (lowest index on ties).
fn furthest_from_set(points: &[Vec3], chosen: &[usize]) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for candidate in 0..points.len() {
        if chosen.contains(&candidate) {
            continue;
        }
        let score = distance_to_set_sq(points, chosen, candidate);
        if best.is_none_or(|(_, best_score)| score > best_score) {
            best = Some((candidate, score));
        }
    }
    best.map(|(index, _)| index)
}

/// Unchosen point forming the largest triangle with the first two chosen points.
fn largest_triangle_with_first_two(points: &[Vec3], chosen: &[usize]) -> Option<usize> {
    let origin = points[chosen[0]];
    let edge = points[chosen[1]] - origin;
    let mut best: Option<(usize, f64)> = None;
    for (candidate, point) in points.iter().enumerate() {
        if chosen.contains(&candidate) {
            continue;
        }
        let score = edge.cross(&(point - origin)).norm_squared();
        if best.is_none_or(|(_, best_score)| score > best_score) {
            best = Some((candidate, score));
        }
    }
    best.map(|(index, _)| index)
}

/// Fourth point of the reduction: the unchosen point on the opposite side of the
/// edge `chosen[0]`-`chosen[1]` from `chosen[2]` that spans the largest triangle
/// with that edge, which turns the triangle into a quadrilateral. When every
/// remaining point is on the third point's side, the point furthest from the
/// chosen set is used instead.
fn quadrilateral_extension(points: &[Vec3], chosen: &[usize]) -> Option<usize> {
    let origin = points[chosen[0]];
    let edge = points[chosen[1]] - origin;
    let normal = edge.cross(&(points[chosen[2]] - origin));
    let mut best: Option<(usize, f64)> = None;
    for (candidate, point) in points.iter().enumerate() {
        if chosen.contains(&candidate) {
            continue;
        }
        let spanned = edge.cross(&(point - origin));
        if spanned.dot(&normal) >= 0.0 {
            continue;
        }
        let score = spanned.norm_squared();
        if best.is_none_or(|(_, best_score)| score > best_score) {
            best = Some((candidate, score));
        }
    }
    match best {
        Some((index, _)) => Some(index),
        None => furthest_from_set(points, chosen),
    }
}

/// Reduce `contacts` to at most `max` points that cover the contact area.
///
/// Fewer than `max` contacts are returned unchanged. Otherwise the deepest
/// contact is kept, then the contact furthest from it, then the one that
/// maximises the triangle area, then the one that completes the largest
/// quadrilateral (further points, if `max > 4`, are the ones furthest from the
/// kept set). Positions are the midpoints of the two witness points. The kept
/// contacts retain their original relative order, so the result is
/// deterministic.
pub(super) fn reduce_contacts(contacts: Vec<Contact>, max: usize) -> Vec<Contact> {
    if contacts.len() <= max {
        return contacts;
    }
    if max == 0 {
        return Vec::new();
    }
    let points: Vec<Vec3> = contacts.iter().map(Contact::point).collect();
    let mut deepest = 0;
    for (index, contact) in contacts.iter().enumerate() {
        if contact.depth > contacts[deepest].depth {
            deepest = index;
        }
    }
    let mut chosen = vec![deepest];
    while chosen.len() < max {
        let next = match chosen.len() {
            2 => largest_triangle_with_first_two(&points, &chosen),
            3 => quadrilateral_extension(&points, &chosen),
            _ => furthest_from_set(&points, &chosen),
        };
        match next {
            Some(index) => chosen.push(index),
            None => break,
        }
    }
    chosen.sort_unstable();
    chosen
        .into_iter()
        .map(|index| contacts[index].clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_geometry::{BoxShape, Sphere};

    fn contact_at(x: f64, z: f64, depth: f64) -> Contact {
        let p = Vec3::new(x, 0.0, z);
        Contact::new(p, p, Vec3::new(0.0, -1.0, 0.0), depth)
    }

    #[test]
    fn few_contacts_pass_through_unchanged() {
        let contacts = vec![contact_at(0.0, 0.0, 0.1), contact_at(1.0, 0.0, 0.2)];
        let kept = reduce_contacts(contacts.clone(), 4);
        assert_eq!(kept.len(), 2);
        assert_eq!(kept[0].point_a, contacts[0].point_a);
        assert_eq!(kept[1].point_a, contacts[1].point_a);
        assert!(reduce_contacts(contacts, 0).is_empty());
    }

    #[test]
    fn duplicated_footprint_corners_reduce_to_the_four_corners() {
        // Two pieces each report the corners of their half of a 1 x 1 footprint,
        // so the seam corners appear twice (8 contacts, equal depth).
        let mut contacts = Vec::new();
        for (x0, x1) in [(-0.5, 0.0), (0.0, 0.5)] {
            for x in [x0, x1] {
                for z in [-0.5, 0.5] {
                    contacts.push(contact_at(x, z, 0.1));
                }
            }
        }
        let kept = reduce_contacts(contacts, 4);
        assert_eq!(kept.len(), 4);
        let mut corners: Vec<(i32, i32)> = kept
            .iter()
            .map(|c| {
                (
                    (c.point_a.x * 2.0).round() as i32,
                    (c.point_a.z * 2.0).round() as i32,
                )
            })
            .collect();
        corners.sort_unstable();
        assert_eq!(corners, [(-1, -1), (-1, 1), (1, -1), (1, 1)]);
    }

    #[test]
    fn the_deepest_contact_is_always_kept() {
        let mut contacts: Vec<Contact> = (0..6)
            .map(|k| contact_at(k as f64 * 0.01, 0.0, 0.05))
            .collect();
        contacts[3].depth = 0.5;
        let kept = reduce_contacts(contacts, 3);
        assert_eq!(kept.len(), 3);
        assert!(kept.iter().any(|c| (c.depth - 0.5).abs() < 1e-15));
    }

    #[test]
    fn collinear_contacts_still_reduce_to_distinct_points() {
        let contacts: Vec<Contact> = (0..7).map(|k| contact_at(k as f64, 0.0, 0.1)).collect();
        let kept = reduce_contacts(contacts, 4);
        assert_eq!(kept.len(), 4);
        let mut xs: Vec<i64> = kept.iter().map(|c| c.point_a.x.round() as i64).collect();
        xs.sort_unstable();
        xs.dedup();
        assert_eq!(xs.len(), 4);
        // The extremes of the line are kept.
        assert_eq!(xs.first(), Some(&0));
        assert_eq!(xs.last(), Some(&6));
    }

    #[test]
    fn reduction_beyond_four_points_keeps_spreading() {
        let contacts: Vec<Contact> = (0..9).map(|k| contact_at(k as f64, 0.0, 0.1)).collect();
        let kept = reduce_contacts(contacts, 6);
        assert_eq!(kept.len(), 6);
    }

    #[test]
    fn bounds_of_a_translated_sphere_in_the_mesh_frame() {
        let sphere = Sphere::new(0.5);
        let convex_tf = Transform::from_position(Vec3::new(3.0, 1.0, -2.0));
        let mesh_tf = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let (min, max) = convex_bounds_in_mesh_frame(&sphere, &convex_tf, &mesh_tf);
        let expected_min = [1.5, 0.5, -2.5];
        let expected_max = [2.5, 1.5, -1.5];
        for (got, want) in min.iter().zip(expected_min.iter()) {
            assert!((got - want).abs() < 1e-12, "min {min:?}");
        }
        for (got, want) in max.iter().zip(expected_max.iter()) {
            assert!((got - want).abs() < 1e-12, "max {max:?}");
        }
    }

    #[test]
    fn bounds_of_a_rotated_box_enclose_its_rotated_corners() {
        // Box (half extents 1, 0.5, 0.25) turned 90 degrees about z: the x and y
        // extents swap.
        let shape = BoxShape::new(Vec3::new(1.0, 0.5, 0.25));
        let convex_tf = Transform::from_axis_angle(
            Vec3::zeros(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::FRAC_PI_2,
        );
        let mesh_tf = Transform::from_position(Vec3::zeros());
        let (min, max) = convex_bounds_in_mesh_frame(&shape, &convex_tf, &mesh_tf);
        assert!(
            (max[0] - 0.5).abs() < 1e-12 && (min[0] + 0.5).abs() < 1e-12,
            "{min:?} {max:?}"
        );
        assert!(
            (max[1] - 1.0).abs() < 1e-12 && (min[1] + 1.0).abs() < 1e-12,
            "{min:?} {max:?}"
        );
        assert!((max[2] - 0.25).abs() < 1e-12 && (min[2] + 0.25).abs() < 1e-12);
    }

    #[test]
    fn bounds_overlap_is_inclusive_and_per_axis() {
        let (lo, hi) = ([0.0; 3], [1.0; 3]);
        assert!(
            bounds_overlap(lo, hi, [1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
            "touching"
        );
        assert!(!bounds_overlap(lo, hi, [1.1, 0.0, 0.0], [2.0, 1.0, 1.0]));
        assert!(!bounds_overlap(lo, hi, [0.0, 0.0, -2.0], [1.0, 1.0, -0.5]));
        assert!(bounds_overlap(
            lo,
            hi,
            [f64::NEG_INFINITY; 3],
            [f64::INFINITY; 3]
        ));
        // An empty (inverted) piece box never overlaps.
        assert!(!bounds_overlap(lo, hi, [f64::MAX; 3], [f64::MIN; 3]));
    }
}
