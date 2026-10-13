// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Narrow-phase contact generation for every supported shape pair.
//!
//! ## Contact convention
//!
//! [`collide`] returns a [`Contact`] whose `normal` is the unit vector pointing
//! **from shape B toward shape A**: translating A along `+normal` (or B along
//! `-normal`) by `depth` separates the pair.  This is the convention
//! documented on [`ContactResult`](crate::types::ContactResult) and consumed by
//! the impulse resolver.
//!
//! For contacts that touch in a single point (sphere, capsule, edge-edge and
//! convex-solver contacts) `point_a` is the point of A's surface that
//! penetrates deepest into B and `point_b` the matching point on B's surface,
//! so `(point_a - point_b) . normal == -depth`.  Contacts that are spread over
//! a face or an edge (box-box face contacts and every shape against a plane)
//! report the penetration-weighted centre of the touching region instead; the
//! two points are then separated along the normal by the penetration at that
//! centre, which does not exceed `depth` (the deepest penetration).
//!
//! The dispatcher reorders every pair so that the lower-ranked shape (see
//! [`Shape::rank`]) is A; the (A, B) routines below therefore only exist for
//! one ordering of each pair and the opposite ordering is obtained by
//! flipping the result.
//!
//! | pair                                   | algorithm                         |
//! |----------------------------------------|-----------------------------------|
//! | sphere-sphere, sphere-capsule          | closest point on the axis         |
//! | sphere-box                             | closest point in the box frame    |
//! | capsule-capsule                        | closest points of the two axes    |
//! | box-box                                | SAT (15 axes) + face clipping     |
//! | sphere/capsule/box/cylinder - plane    | support points below the plane    |
//! | every other convex pair (cylinders,    | GJK + EPA on rounded cores        |
//! | capsule-box)                           | ([`gjk`](super::gjk))             |

use super::geom::{
    Vec3, add, any_perpendicular, clamp_to_box, closest_point_on_segment, closest_points_segments,
    cross, dist_sq, dot, length_sq, madd, neg, normalized, scale, sub, unit_perpendicular_part,
};
use super::gjk::convex_contact;
use super::shapes::{Shape, WorldShape};

/// Normal used when two shapes are exactly concentric and no direction is
/// defined by the geometry.
const FALLBACK_NORMAL: Vec3 = [0.0, 1.0, 0.0];

/// A single contact point between two shapes.
#[derive(Debug, Clone, Copy)]
pub(super) struct Contact {
    /// Unit normal pointing from shape B toward shape A.
    pub normal: Vec3,
    /// Penetration depth (always positive).
    pub depth: f64,
    /// Deepest point of A inside B, on A's surface.
    pub point_a: Vec3,
    /// Matching point on B's surface.
    pub point_b: Vec3,
}

impl Contact {
    /// The same contact seen with the roles of A and B exchanged.
    pub(super) fn flipped(self) -> Self {
        Self {
            normal: neg(self.normal),
            depth: self.depth,
            point_a: self.point_b,
            point_b: self.point_a,
        }
    }
}

/// Conservative bounding-sphere test: `false` guarantees that the two shapes
/// cannot touch (planes have an infinite bound and always pass).
pub(super) fn bounds_overlap(a: &WorldShape, b: &WorldShape) -> bool {
    if a.bound_radius.is_finite() && b.bound_radius.is_finite() {
        let reach = a.bound_radius + b.bound_radius;
        dist_sq(a.position, b.position) <= reach * reach
    } else {
        true
    }
}

/// Generate the contact between two placed shapes, if they penetrate.
///
/// The returned contact refers to `a` as body A and `b` as body B.
pub(super) fn collide(a: &WorldShape, b: &WorldShape) -> Option<Contact> {
    if !bounds_overlap(a, b) {
        return None;
    }
    if a.shape.rank() > b.shape.rank() {
        collide_ordered(b, a).map(Contact::flipped)
    } else {
        collide_ordered(a, b)
    }
}

/// Dispatch for pairs already ordered by [`Shape::rank`].
fn collide_ordered(a: &WorldShape, b: &WorldShape) -> Option<Contact> {
    match (a.shape, b.shape) {
        (Shape::Sphere { radius: ra }, Shape::Sphere { radius: rb }) => {
            sphere_sphere(a.position, ra, b.position, rb, FALLBACK_NORMAL)
        }
        (Shape::Sphere { radius }, Shape::Capsule { radius: rb, .. }) => {
            sphere_capsule(a.position, radius, b, rb)
        }
        (Shape::Sphere { radius }, Shape::Box { half_extents }) => {
            sphere_box(a.position, radius, b, half_extents)
        }
        (Shape::Sphere { radius }, Shape::Plane { normal, offset }) => {
            sphere_plane(a.position, radius, normal, offset)
        }
        (Shape::Capsule { radius: ra, .. }, Shape::Capsule { radius: rb, .. }) => {
            capsule_capsule(a, ra, b, rb)
        }
        (Shape::Capsule { radius, .. }, Shape::Plane { normal, offset }) => {
            capsule_plane(a, radius, normal, offset)
        }
        (
            Shape::Box { half_extents },
            Shape::Box {
                half_extents: other,
            },
        ) => box_box(a, half_extents, b, other),
        (Shape::Box { half_extents }, Shape::Plane { normal, offset }) => {
            box_plane(a, half_extents, normal, offset)
        }
        (
            Shape::Cylinder {
                radius,
                half_height,
            },
            Shape::Plane { normal, offset },
        ) => cylinder_plane(a, radius, half_height, normal, offset),
        // Every pair with a plane as shape B is matched above, and a plane is
        // the highest-ranked shape, so only plane-plane pairs reach this arm
        // (an infinite half-space has no contact with another half-space).
        (Shape::Plane { .. }, _) => None,
        // Remaining pairs are all convex-convex and involve a capsule, box or
        // cylinder: sphere-cylinder, capsule-box, capsule-cylinder,
        // box-cylinder and cylinder-cylinder.
        _ => {
            let (core_a, margin_a) = a.core()?;
            let (core_b, margin_b) = b.core()?;
            convex_contact(&core_a, margin_a, &core_b, margin_b)
        }
    }
}

/// Sphere A against sphere B.
///
/// `fallback_normal` is used when the centres coincide.
fn sphere_sphere(
    pos_a: Vec3,
    radius_a: f64,
    pos_b: Vec3,
    radius_b: f64,
    fallback_normal: Vec3,
) -> Option<Contact> {
    let d = sub(pos_a, pos_b);
    let dist_sq = length_sq(d);
    let reach = radius_a + radius_b;
    if dist_sq >= reach * reach {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = if dist > 1e-12 {
        scale(d, 1.0 / dist)
    } else {
        fallback_normal
    };
    Some(Contact {
        normal,
        depth: reach - dist,
        point_a: madd(pos_a, normal, -radius_a),
        point_b: madd(pos_b, normal, radius_b),
    })
}

/// Sphere A against a half-space plane B with unit `normal`.
fn sphere_plane(pos: Vec3, radius: f64, normal: Vec3, offset: f64) -> Option<Contact> {
    let signed = dot(normal, pos) - offset;
    let depth = radius - signed;
    if depth <= 0.0 {
        return None;
    }
    Some(Contact {
        normal,
        depth,
        point_a: madd(pos, normal, -radius),
        point_b: madd(pos, normal, -signed),
    })
}

/// Sphere A against capsule B.
fn sphere_capsule(
    sphere_pos: Vec3,
    sphere_radius: f64,
    capsule: &WorldShape,
    capsule_radius: f64,
) -> Option<Contact> {
    let (e0, e1) = capsule.axis_endpoints()?;
    let closest = closest_point_on_segment(sphere_pos, e0, e1);
    sphere_sphere(
        sphere_pos,
        sphere_radius,
        closest,
        capsule_radius,
        capsule.any_perpendicular_to_axis(),
    )
}

/// Sphere A against an oriented box B.
fn sphere_box(
    sphere_pos: Vec3,
    radius: f64,
    box_shape: &WorldShape,
    half: Vec3,
) -> Option<Contact> {
    let local = box_shape.to_local_point(sphere_pos);
    let clamped = clamp_to_box(local, half);
    let diff = sub(local, clamped);
    let dist_sq = length_sq(diff);
    let (normal_local, depth, surface_local) = if dist_sq > 0.0 {
        // Centre outside the box: the closest point of the box is the contact.
        if dist_sq >= radius * radius {
            return None;
        }
        let dist = dist_sq.sqrt();
        (scale(diff, 1.0 / dist), radius - dist, clamped)
    } else {
        // Centre inside (or on the surface of) the box: leave through the
        // nearest face.
        let gaps = [
            half[0] - local[0].abs(),
            half[1] - local[1].abs(),
            half[2] - local[2].abs(),
        ];
        let mut axis = 0;
        if gaps[1] < gaps[axis] {
            axis = 1;
        }
        if gaps[2] < gaps[axis] {
            axis = 2;
        }
        let sign = if local[axis] >= 0.0 { 1.0 } else { -1.0 };
        let mut normal = [0.0; 3];
        normal[axis] = sign;
        let mut surface = local;
        surface[axis] = sign * half[axis];
        (normal, radius + gaps[axis], surface)
    };
    let normal = box_shape.to_world_dir(normal_local);
    Some(Contact {
        normal,
        depth,
        point_a: madd(sphere_pos, normal, -radius),
        point_b: box_shape.to_world_point(surface_local),
    })
}

/// Capsule A against capsule B.
fn capsule_capsule(
    a: &WorldShape,
    radius_a: f64,
    b: &WorldShape,
    radius_b: f64,
) -> Option<Contact> {
    let (a0, a1) = a.axis_endpoints()?;
    let (b0, b1) = b.axis_endpoints()?;
    let (ca, cb) = closest_points_segments(a0, a1, b0, b1);
    let d = sub(ca, cb);
    let dist_sq = length_sq(d);
    let reach = radius_a + radius_b;
    if dist_sq >= reach * reach {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = if dist > 1e-12 {
        scale(d, 1.0 / dist)
    } else {
        // The axes intersect: the cheapest way out is perpendicular to both.
        // Parallel axes have no such direction (their cross product is only
        // rounding noise), any direction perpendicular to the axis will do.
        let axes_cross = cross(a.axis(1), b.axis(1));
        let mut n = if length_sq(axes_cross) > 1e-24 {
            normalized(axes_cross)
        } else {
            None
        }
        .unwrap_or_else(|| any_perpendicular(a.axis(1)));
        if dot(n, sub(a.position, b.position)) < 0.0 {
            n = neg(n);
        }
        n
    };
    Some(Contact {
        normal,
        depth: reach - dist,
        point_a: madd(ca, normal, -radius_a),
        point_b: madd(cb, normal, radius_b),
    })
}

/// Contact of a set of surface points (all belonging to shape A) with a plane.
///
/// Every point whose signed distance to the plane is negative penetrates; the
/// reported contact point is the penetration-weighted centroid of those points
/// (the centre of pressure of a linear spring model) and the depth is the
/// deepest penetration.
fn plane_contact(points: &[Vec3], normal: Vec3, offset: f64) -> Option<Contact> {
    let mut weight_sum = 0.0;
    let mut weighted = [0.0; 3];
    let mut depth: f64 = 0.0;
    for p in points {
        let pen = offset - dot(normal, *p);
        if pen > 0.0 {
            weight_sum += pen;
            weighted = madd(weighted, *p, pen);
            depth = depth.max(pen);
        }
    }
    if weight_sum <= 0.0 {
        return None;
    }
    let point_a = scale(weighted, 1.0 / weight_sum);
    let pen_at_point = offset - dot(normal, point_a);
    Some(Contact {
        normal,
        depth,
        point_a,
        point_b: madd(point_a, normal, pen_at_point),
    })
}

/// Capsule A against a half-space plane B.
fn capsule_plane(a: &WorldShape, radius: f64, normal: Vec3, offset: f64) -> Option<Contact> {
    let (e0, e1) = a.axis_endpoints()?;
    // Lowest point of each end cap along the plane normal.
    let lowest = [madd(e0, normal, -radius), madd(e1, normal, -radius)];
    plane_contact(&lowest, normal, offset)
}

/// Box A against a half-space plane B (all eight corners are tested).
fn box_plane(a: &WorldShape, half: Vec3, normal: Vec3, offset: f64) -> Option<Contact> {
    let axes = [a.axis(0), a.axis(1), a.axis(2)];
    let mut corners = [[0.0; 3]; 8];
    for (index, corner) in corners.iter_mut().enumerate() {
        let mut p = a.position;
        for (bit, (axis, h)) in axes.iter().zip(half.iter()).enumerate() {
            let sign = if ((index >> bit) & 1) == 1 { 1.0 } else { -1.0 };
            p = madd(p, *axis, sign * h);
        }
        *corner = p;
    }
    plane_contact(&corners, normal, offset)
}

/// Cylinder A against a half-space plane B.
///
/// The lowest point of the cylinder is used: a rim point for a tilted
/// cylinder, the cap centre for one standing on its cap and the middle of the
/// side line for one lying on its side.
fn cylinder_plane(
    a: &WorldShape,
    radius: f64,
    half_height: f64,
    normal: Vec3,
    offset: f64,
) -> Option<Contact> {
    const TIE: f64 = 1e-9;
    let axis = a.axis(1);
    let along = dot(normal, axis);
    let end = if along > TIE {
        -half_height
    } else if along < -TIE {
        half_height
    } else {
        0.0
    };
    let mut lowest = madd(a.position, axis, end);
    // For a cylinder standing on the plane (axis parallel to the normal) the
    // whole cap is lowest and its centre represents it.
    if let Some(radial) = unit_perpendicular_part(normal, axis) {
        lowest = madd(lowest, radial, -radius);
    }
    plane_contact(&[lowest], normal, offset)
}

/// Which separating axis produced the minimum overlap of two boxes.
#[derive(Debug, Clone, Copy)]
enum BoxAxis {
    /// Face normal `i` of box A.
    FaceA(usize),
    /// Face normal `j` of box B.
    FaceB(usize),
    /// Cross product of edge `i` of A and edge `j` of B.
    Edge(usize, usize),
}

/// An oriented box described by its centre, axes and half-extents.
#[derive(Debug, Clone, Copy)]
struct BoxFrame {
    /// Box centre.
    position: Vec3,
    /// Unit direction of each box axis in world space.
    axes: [Vec3; 3],
    /// Half-extents along `axes`.
    half: Vec3,
}

impl BoxFrame {
    fn new(shape: &WorldShape, half: Vec3) -> Self {
        Self {
            position: shape.position,
            axes: [shape.axis(0), shape.axis(1), shape.axis(2)],
            half,
        }
    }

    /// Radius of the box projected on the unit direction `dir`.
    fn projected_radius(&self, dir: Vec3) -> f64 {
        self.axes
            .iter()
            .zip(self.half.iter())
            .map(|(axis, h)| h * dot(*axis, dir).abs())
            .sum()
    }

    /// Centre of the edge parallel to axis `skip` that lies furthest along
    /// `dir`.
    fn support_edge_center(&self, skip: usize, dir: Vec3) -> Vec3 {
        let mut center = self.position;
        for (k, (axis, h)) in self.axes.iter().zip(self.half.iter()).enumerate() {
            if k != skip {
                let s = if dot(dir, *axis) >= 0.0 { 1.0 } else { -1.0 };
                center = madd(center, *axis, s * h);
            }
        }
        center
    }
}

/// Oriented box A against oriented box B using the separating axis theorem.
///
/// The axis of least overlap decides the contact type: a face axis yields a
/// face contact (the incident face is clipped against the side planes of the
/// reference face and the contact is the penetration-weighted centroid of the
/// clipped polygon), an edge-edge axis yields the closest points of the two
/// edges.  Face axes are preferred over edge axes of nearly equal overlap so
/// that resting contacts do not flicker between the two.
fn box_box(a: &WorldShape, half_a: Vec3, b: &WorldShape, half_b: Vec3) -> Option<Contact> {
    let box_a = BoxFrame::new(a, half_a);
    let box_b = BoxFrame::new(b, half_b);
    let t = sub(box_b.position, box_a.position);

    let mut best_overlap = f64::INFINITY;
    let mut best_axis = [0.0; 3];
    let mut best_kind = BoxAxis::FaceA(0);

    for (i, axis) in box_a.axes.iter().enumerate() {
        let overlap = box_a.half[i] + box_b.projected_radius(*axis) - dot(t, *axis).abs();
        if overlap < 0.0 {
            return None;
        }
        if overlap < best_overlap {
            best_overlap = overlap;
            best_axis = *axis;
            best_kind = BoxAxis::FaceA(i);
        }
    }
    for (j, axis) in box_b.axes.iter().enumerate() {
        let overlap = box_b.half[j] + box_a.projected_radius(*axis) - dot(t, *axis).abs();
        if overlap < 0.0 {
            return None;
        }
        if overlap < best_overlap {
            best_overlap = overlap;
            best_axis = *axis;
            best_kind = BoxAxis::FaceB(j);
        }
    }
    for (i, axis_a) in box_a.axes.iter().enumerate() {
        for (j, axis_b) in box_b.axes.iter().enumerate() {
            let c = cross(*axis_a, *axis_b);
            let len_sq = length_sq(c);
            if len_sq < 1e-10 {
                // Parallel edges: the axis is covered by the face axes.
                continue;
            }
            let axis = scale(c, 1.0 / len_sq.sqrt());
            let overlap =
                box_a.projected_radius(axis) + box_b.projected_radius(axis) - dot(t, axis).abs();
            if overlap < 0.0 {
                return None;
            }
            if overlap < 0.95 * best_overlap - 1e-12 {
                best_overlap = overlap;
                best_axis = axis;
                best_kind = BoxAxis::Edge(i, j);
            }
        }
    }

    // Orient the axis from B toward A (`t` points from A to B).
    let normal = if dot(best_axis, t) > 0.0 {
        neg(best_axis)
    } else {
        best_axis
    };

    match best_kind {
        BoxAxis::FaceA(i) => Some(box_face_contact(
            &box_a,
            i,
            &box_b,
            neg(normal),
            best_overlap,
            true,
        )),
        BoxAxis::FaceB(j) => Some(box_face_contact(
            &box_b,
            j,
            &box_a,
            normal,
            best_overlap,
            false,
        )),
        BoxAxis::Edge(i, j) => {
            // Support edge of A toward B (direction -normal) and of B toward A.
            let edge_a = box_a.support_edge_center(i, neg(normal));
            let edge_b = box_b.support_edge_center(j, normal);
            let (ca, cb) = closest_points_segments(
                madd(edge_a, box_a.axes[i], -box_a.half[i]),
                madd(edge_a, box_a.axes[i], box_a.half[i]),
                madd(edge_b, box_b.axes[j], -box_b.half[j]),
                madd(edge_b, box_b.axes[j], box_b.half[j]),
            );
            Some(Contact {
                normal,
                depth: best_overlap,
                point_a: ca,
                point_b: cb,
            })
        }
    }
}

/// Clip a convex polygon against the half-space `normal . x <= offset`
/// (Sutherland-Hodgman).
fn clip_polygon(polygon: &[Vec3], normal: Vec3, offset: f64) -> Vec<Vec3> {
    let mut out = Vec::with_capacity(polygon.len() + 2);
    for (idx, current) in polygon.iter().enumerate() {
        let next = polygon[(idx + 1) % polygon.len()];
        let dc = dot(normal, *current) - offset;
        let dn = dot(normal, next) - offset;
        let current_inside = dc <= 0.0;
        let next_inside = dn <= 0.0;
        if current_inside {
            out.push(*current);
        }
        if current_inside != next_inside {
            let t = dc / (dc - dn);
            out.push(madd(*current, sub(next, *current), t));
        }
    }
    out
}

/// Face contact of two boxes.
///
/// `reference` owns the face with index `face` (its outward normal is
/// `ref_normal`, pointing toward `incident`); the incident face of the other
/// box that is most anti-parallel to it is clipped against the reference
/// face's side planes.  `depth` is the SAT overlap along the face normal and
/// `reference_is_a` tells which of the two boxes plays shape A.
fn box_face_contact(
    reference: &BoxFrame,
    face: usize,
    incident: &BoxFrame,
    ref_normal: Vec3,
    depth: f64,
    reference_is_a: bool,
) -> Contact {
    // Incident face: the one whose outward normal is most anti-parallel to the
    // reference normal.
    let mut inc_axis = 0;
    let mut best_dot = 0.0_f64;
    for (j, axis) in incident.axes.iter().enumerate() {
        let d = dot(ref_normal, *axis).abs();
        if d > best_dot {
            best_dot = d;
            inc_axis = j;
        }
    }
    let inc_sign = if dot(ref_normal, incident.axes[inc_axis]) > 0.0 {
        -1.0
    } else {
        1.0
    };
    let face_center = madd(
        incident.position,
        incident.axes[inc_axis],
        inc_sign * incident.half[inc_axis],
    );
    let k = (inc_axis + 1) % 3;
    let l = (inc_axis + 2) % 3;
    let u = scale(incident.axes[k], incident.half[k]);
    let v = scale(incident.axes[l], incident.half[l]);
    let incident_face = [
        add(add(face_center, u), v),
        add(sub(face_center, u), v),
        sub(sub(face_center, u), v),
        sub(add(face_center, u), v),
    ];

    // Clip against the four side planes of the reference face.
    let mut polygon: Vec<Vec3> = incident_face.to_vec();
    for side in 0..3 {
        if side == face {
            continue;
        }
        let axis = reference.axes[side];
        let center = dot(axis, reference.position);
        polygon = clip_polygon(&polygon, axis, center + reference.half[side]);
        if polygon.is_empty() {
            break;
        }
        polygon = clip_polygon(&polygon, neg(axis), -center + reference.half[side]);
        if polygon.is_empty() {
            break;
        }
    }
    if polygon.is_empty() {
        polygon = incident_face.to_vec();
    }

    // Signed distance of each point above the reference face plane.
    let plane_offset = dot(ref_normal, reference.position) + reference.half[face];
    let mut weight_sum = 0.0;
    let mut weighted = [0.0; 3];
    let mut deepest = polygon[0];
    let mut deepest_sep = f64::INFINITY;
    for p in &polygon {
        let sep = dot(ref_normal, *p) - plane_offset;
        if sep < deepest_sep {
            deepest_sep = sep;
            deepest = *p;
        }
        if sep <= 1e-9 {
            let w = (-sep).max(0.0) + 1e-12;
            weight_sum += w;
            weighted = madd(weighted, *p, w);
        }
    }
    let on_incident = if weight_sum > 0.0 {
        scale(weighted, 1.0 / weight_sum)
    } else {
        deepest
    };
    let lift = plane_offset - dot(ref_normal, on_incident);
    let on_reference = madd(on_incident, ref_normal, lift);
    if reference_is_a {
        Contact {
            normal: neg(ref_normal),
            depth,
            point_a: on_reference,
            point_b: on_incident,
        }
    } else {
        Contact {
            normal: ref_normal,
            depth,
            point_a: on_incident,
            point_b: on_reference,
        }
    }
}
