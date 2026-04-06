// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact manifold generation and management.
//!
//! Covers manifold point reduction (4-point), contact normal computation,
//! contact depth (penetration distance), tangent-plane basis, manifold
//! merging/splitting, persistent contact caching, Jacobian computation for
//! contacts, contact manifolds for spheres/boxes/capsules, and manifold
//! quality metrics.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Vec3 helpers (no nalgebra – plain [f64; 3])
// ─────────────────────────────────────────────────────────────────────────────

/// Adds two 3-vectors component-wise.
#[inline]
pub fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtracts `b` from `a` component-wise.
#[inline]
pub fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scales a 3-vector by a scalar.
#[inline]
pub fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Computes the dot product of two 3-vectors.
#[inline]
pub fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Computes the cross product of two 3-vectors.
#[inline]
pub fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Returns the squared length of a 3-vector.
#[inline]
pub fn v3_len_sq(a: [f64; 3]) -> f64 {
    v3_dot(a, a)
}

/// Returns the length of a 3-vector.
#[inline]
pub fn v3_len(a: [f64; 3]) -> f64 {
    v3_len_sq(a).sqrt()
}

/// Normalises a 3-vector. Returns `[1,0,0]` for near-zero vectors.
#[inline]
pub fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let l = v3_len(a);
    if l < f64::EPSILON {
        [1.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

/// Negates a 3-vector.
#[inline]
pub fn v3_neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

/// Linear interpolation between two vectors: `a + t*(b - a)`.
#[inline]
pub fn v3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    v3_add(a, v3_scale(v3_sub(b, a), t))
}

/// Computes the distance between two points.
#[inline]
pub fn v3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    v3_len(v3_sub(a, b))
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact point
// ─────────────────────────────────────────────────────────────────────────────

/// A single contact point within a manifold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactPoint {
    /// World-space position of the contact.
    pub position: [f64; 3],
    /// Local position in body A's frame.
    pub local_a: [f64; 3],
    /// Local position in body B's frame.
    pub local_b: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Accumulated normal impulse for warm starting.
    pub normal_impulse: f64,
    /// Accumulated tangent impulse (friction), component 0.
    pub tangent_impulse0: f64,
    /// Accumulated tangent impulse (friction), component 1.
    pub tangent_impulse1: f64,
    /// Frame counter when this point was last updated.
    pub age: u32,
}

impl ContactPoint {
    /// Creates a new contact point with all impulses zeroed.
    pub fn new(position: [f64; 3], local_a: [f64; 3], local_b: [f64; 3], depth: f64) -> Self {
        Self {
            position,
            local_a,
            local_b,
            depth,
            normal_impulse: 0.0,
            tangent_impulse0: 0.0,
            tangent_impulse1: 0.0,
            age: 0,
        }
    }

    /// Returns `true` if this contact is still penetrating (depth > tolerance).
    pub fn is_penetrating(&self, tolerance: f64) -> bool {
        self.depth > tolerance
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tangent plane basis
// ─────────────────────────────────────────────────────────────────────────────

/// An orthonormal frame attached to a contact.
///
/// `normal` points from body B into body A.
/// `tangent0` and `tangent1` span the contact plane.
#[derive(Debug, Clone, Copy)]
pub struct ContactBasis {
    /// Contact normal (unit vector).
    pub normal: [f64; 3],
    /// First tangent vector (unit vector, perpendicular to `normal`).
    pub tangent0: [f64; 3],
    /// Second tangent vector (unit vector, perpendicular to both).
    pub tangent1: [f64; 3],
}

impl ContactBasis {
    /// Constructs the contact basis from a (possibly un-normalised) normal.
    ///
    /// Builds an arbitrary but stable orthonormal frame using the
    /// "smallest component" method.
    pub fn from_normal(n: [f64; 3]) -> Self {
        let normal = v3_normalize(n);
        let tangent0 = compute_tangent(normal);
        let tangent1 = v3_cross(normal, tangent0);
        Self {
            normal,
            tangent0,
            tangent1,
        }
    }

    /// Returns `true` if the basis vectors are mutually orthogonal (within `tol`).
    pub fn is_orthonormal(&self, tol: f64) -> bool {
        let d01 = v3_dot(self.normal, self.tangent0).abs();
        let d02 = v3_dot(self.normal, self.tangent1).abs();
        let d12 = v3_dot(self.tangent0, self.tangent1).abs();
        let l0 = (v3_len(self.normal) - 1.0).abs();
        let l1 = (v3_len(self.tangent0) - 1.0).abs();
        let l2 = (v3_len(self.tangent1) - 1.0).abs();
        d01 < tol && d02 < tol && d12 < tol && l0 < tol && l1 < tol && l2 < tol
    }
}

/// Builds a unit tangent vector perpendicular to `n` using the smallest component.
pub fn compute_tangent(n: [f64; 3]) -> [f64; 3] {
    let abs_x = n[0].abs();
    let abs_y = n[1].abs();
    let abs_z = n[2].abs();
    let perp = if abs_x <= abs_y && abs_x <= abs_z {
        [1.0_f64, 0.0, 0.0]
    } else if abs_y <= abs_z {
        [0.0_f64, 1.0, 0.0]
    } else {
        [0.0_f64, 0.0, 1.0]
    };
    let t = v3_cross(n, perp);
    v3_normalize(t)
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact manifold
// ─────────────────────────────────────────────────────────────────────────────

/// Maximum number of contact points per manifold.
pub const MAX_MANIFOLD_POINTS: usize = 4;

/// A contact manifold between two bodies.
///
/// Maintains up to `MAX_MANIFOLD_POINTS` contact points and the shared
/// contact normal / tangent basis.
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Body A index.
    pub body_a: usize,
    /// Body B index.
    pub body_b: usize,
    /// Shared contact basis (normal + tangents).
    pub basis: ContactBasis,
    /// Active contact points (up to `MAX_MANIFOLD_POINTS`).
    pub points: [Option<ContactPoint>; MAX_MANIFOLD_POINTS],
    /// Number of active contact points.
    pub point_count: usize,
    /// Frame counter for the last update.
    pub last_update: u32,
}

impl ContactManifold {
    /// Creates an empty manifold between `body_a` and `body_b` with the given normal.
    pub fn new(body_a: usize, body_b: usize, normal: [f64; 3]) -> Self {
        Self {
            body_a,
            body_b,
            basis: ContactBasis::from_normal(normal),
            points: [None; MAX_MANIFOLD_POINTS],
            point_count: 0,
            last_update: 0,
        }
    }

    /// Returns an iterator over the active contact points.
    pub fn active_points(&self) -> impl Iterator<Item = &ContactPoint> {
        self.points[..self.point_count]
            .iter()
            .filter_map(|p| p.as_ref())
    }

    /// Returns the number of active contact points.
    pub fn len(&self) -> usize {
        self.point_count
    }

    /// Returns `true` if the manifold has no contact points.
    pub fn is_empty(&self) -> bool {
        self.point_count == 0
    }

    /// Adds a contact point, reducing to 4 if needed.
    ///
    /// If there are already `MAX_MANIFOLD_POINTS` points the new point
    /// is reduced into the set using [`reduce_to_4`].
    pub fn add_point(&mut self, cp: ContactPoint) {
        if self.point_count < MAX_MANIFOLD_POINTS {
            self.points[self.point_count] = Some(cp);
            self.point_count += 1;
        } else {
            // collect existing active points + new
            let mut pts: Vec<ContactPoint> = self.points.iter().filter_map(|p| *p).collect();
            pts.push(cp);
            let reduced = reduce_to_4(&pts);
            for (i, rp) in reduced.iter().enumerate() {
                self.points[i] = Some(*rp);
            }
            self.point_count = reduced.len();
        }
    }

    /// Clears all contact points.
    pub fn clear(&mut self) {
        self.points = [None; MAX_MANIFOLD_POINTS];
        self.point_count = 0;
    }

    /// Updates the contact normal and rebuilds the basis.
    pub fn update_normal(&mut self, normal: [f64; 3]) {
        self.basis = ContactBasis::from_normal(normal);
    }

    /// Computes the average contact position across all active points.
    pub fn average_position(&self) -> [f64; 3] {
        let count = self.point_count;
        if count == 0 {
            return [0.0; 3];
        }
        let sum = self
            .active_points()
            .fold([0.0_f64; 3], |acc, p| v3_add(acc, p.position));
        v3_scale(sum, 1.0 / count as f64)
    }

    /// Returns the deepest penetrating contact point, if any.
    pub fn deepest_point(&self) -> Option<&ContactPoint> {
        self.active_points().max_by(|a, b| {
            a.depth
                .partial_cmp(&b.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Increments the `age` of all active contact points by 1.
    pub fn age_points(&mut self) {
        for p in self.points[..self.point_count].iter_mut().flatten() {
            p.age += 1;
        }
    }

    /// Removes contact points older than `max_age` frames.
    pub fn prune_old_points(&mut self, max_age: u32) {
        let mut new_count = 0;
        let mut new_points: [Option<ContactPoint>; MAX_MANIFOLD_POINTS] =
            [None; MAX_MANIFOLD_POINTS];
        for slot in &self.points[..self.point_count] {
            if let Some(p) = slot
                && p.age <= max_age
            {
                new_points[new_count] = Some(*p);
                new_count += 1;
            }
        }
        self.points = new_points;
        self.point_count = new_count;
    }

    /// Computes the signed separation along `basis.normal` for a given point.
    ///
    /// Positive values indicate penetration.
    pub fn signed_depth_at(&self, pos: [f64; 3]) -> f64 {
        v3_dot(v3_sub(pos, self.average_position()), self.basis.normal)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold point reduction to 4 points
// ─────────────────────────────────────────────────────────────────────────────

/// Reduces a list of contact points to at most 4 by maximising the area of the
/// convex hull of the contact positions projected onto the contact plane.
///
/// Algorithm:
/// 1. Keep the deepest point.
/// 2. Keep the point farthest from the deepest.
/// 3. Keep the point maximising area with the first two.
/// 4. Keep the point maximising area with the current triangle.
pub fn reduce_to_4(points: &[ContactPoint]) -> Vec<ContactPoint> {
    if points.len() <= MAX_MANIFOLD_POINTS {
        return points.to_vec();
    }

    let mut result: Vec<ContactPoint> = Vec::with_capacity(MAX_MANIFOLD_POINTS);

    // 1. deepest point
    let Some(&p0) = points.iter().max_by(|a, b| {
        a.depth
            .partial_cmp(&b.depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    }) else {
        return result;
    };
    result.push(p0);

    // 2. farthest from p0
    let Some(&p1) = points.iter().max_by(|a, b| {
        v3_dist(a.position, p0.position)
            .partial_cmp(&v3_dist(b.position, p0.position))
            .unwrap_or(std::cmp::Ordering::Equal)
    }) else {
        return result;
    };
    result.push(p1);

    // 3. point maximising area of triangle with p0, p1
    if let Some(p2) = points.iter().max_by(|a, b| {
        triangle_area_sq(p0.position, p1.position, a.position)
            .partial_cmp(&triangle_area_sq(p0.position, p1.position, b.position))
            .unwrap_or(std::cmp::Ordering::Equal)
    }) {
        result.push(*p2);

        // 4. point maximising area with existing triangle
        if let Some(p3) = points.iter().max_by(|a, b| {
            let area_a = triangle_area_sq(p0.position, p1.position, a.position)
                + triangle_area_sq(p1.position, p2.position, a.position);
            let area_b = triangle_area_sq(p0.position, p1.position, b.position)
                + triangle_area_sq(p1.position, p2.position, b.position);
            area_a
                .partial_cmp(&area_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
            result.push(*p3);
        }
    }

    result
}

/// Returns the squared area of the triangle with vertices `a`, `b`, `c`.
pub fn triangle_area_sq(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = v3_sub(b, a);
    let ac = v3_sub(c, a);
    let cross = v3_cross(ab, ac);
    v3_len_sq(cross) * 0.25
}

/// Returns the area of the triangle with vertices `a`, `b`, `c`.
pub fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    triangle_area_sq(a, b, c).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold generation: sphere vs sphere
// ─────────────────────────────────────────────────────────────────────────────

/// Generates a contact manifold between two spheres.
///
/// Returns `None` if the spheres are separated.
/// `centre_a`, `centre_b` are world-space positions; `ra`, `rb` are radii.
pub fn manifold_sphere_sphere(
    body_a: usize,
    body_b: usize,
    centre_a: [f64; 3],
    ra: f64,
    centre_b: [f64; 3],
    rb: f64,
) -> Option<ContactManifold> {
    let diff = v3_sub(centre_b, centre_a);
    let dist_sq = v3_len_sq(diff);
    let sum_r = ra + rb;
    if dist_sq >= sum_r * sum_r {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = if dist < f64::EPSILON {
        [0.0, 1.0, 0.0]
    } else {
        v3_normalize(diff)
    };
    let depth = sum_r - dist;
    let contact_pos = v3_add(centre_a, v3_scale(normal, ra - depth * 0.5));
    let local_a = v3_scale(normal, ra);
    let local_b = v3_scale(v3_neg(normal), rb);

    let mut manifold = ContactManifold::new(body_a, body_b, normal);
    manifold.add_point(ContactPoint::new(contact_pos, local_a, local_b, depth));
    Some(manifold)
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold generation: sphere vs box (AABB)
// ─────────────────────────────────────────────────────────────────────────────

/// Generates a contact manifold between a sphere and an axis-aligned box (AABB).
///
/// `box_centre` is the centre of the AABB; `half_extents` is `[hx, hy, hz]`.
pub fn manifold_sphere_aabb(
    body_sphere: usize,
    body_box: usize,
    sphere_centre: [f64; 3],
    radius: f64,
    box_centre: [f64; 3],
    half_extents: [f64; 3],
) -> Option<ContactManifold> {
    // Closest point on the AABB to sphere centre
    let closest = [
        sphere_centre[0].clamp(
            box_centre[0] - half_extents[0],
            box_centre[0] + half_extents[0],
        ),
        sphere_centre[1].clamp(
            box_centre[1] - half_extents[1],
            box_centre[1] + half_extents[1],
        ),
        sphere_centre[2].clamp(
            box_centre[2] - half_extents[2],
            box_centre[2] + half_extents[2],
        ),
    ];

    let diff = v3_sub(sphere_centre, closest);
    let dist_sq = v3_len_sq(diff);
    if dist_sq >= radius * radius {
        return None;
    }

    let dist = dist_sq.sqrt();
    let normal = if dist < f64::EPSILON {
        // centre inside box – find minimum overlap axis
        penetration_axis_aabb(sphere_centre, box_centre, half_extents)
    } else {
        v3_normalize(diff)
    };

    let depth = radius - dist;
    let contact_pos = closest;
    let local_a = v3_scale(v3_neg(normal), radius);
    let local_b = v3_sub(closest, box_centre);

    let mut manifold = ContactManifold::new(body_sphere, body_box, normal);
    manifold.add_point(ContactPoint::new(contact_pos, local_a, local_b, depth));
    Some(manifold)
}

/// Returns the axis of minimum penetration for a point inside an AABB.
fn penetration_axis_aabb(
    point: [f64; 3],
    box_centre: [f64; 3],
    half_extents: [f64; 3],
) -> [f64; 3] {
    let local = v3_sub(point, box_centre);
    let overlaps = [
        half_extents[0] - local[0].abs(),
        half_extents[1] - local[1].abs(),
        half_extents[2] - local[2].abs(),
    ];
    let min_axis = overlaps
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let mut axis = [0.0_f64; 3];
    axis[min_axis] = if local[min_axis] >= 0.0 { 1.0 } else { -1.0 };
    axis
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold generation: capsule vs sphere
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the closest point on the line segment `[p, q]` to point `t`.
pub fn closest_point_on_segment(p: [f64; 3], q: [f64; 3], t: [f64; 3]) -> [f64; 3] {
    let pq = v3_sub(q, p);
    let len_sq = v3_len_sq(pq);
    if len_sq < f64::EPSILON {
        return p;
    }
    let s = v3_dot(v3_sub(t, p), pq) / len_sq;
    let s_clamped = s.clamp(0.0, 1.0);
    v3_add(p, v3_scale(pq, s_clamped))
}

/// Generates a contact manifold between a capsule and a sphere.
///
/// `cap_a` and `cap_b` are the endpoint centres of the capsule with radius `cap_r`.
/// `sphere_centre` with `sphere_r` is the sphere.
pub fn manifold_capsule_sphere(
    body_cap: usize,
    body_sphere: usize,
    cap_a: [f64; 3],
    cap_b: [f64; 3],
    cap_r: f64,
    sphere_centre: [f64; 3],
    sphere_r: f64,
) -> Option<ContactManifold> {
    let closest = closest_point_on_segment(cap_a, cap_b, sphere_centre);
    let diff = v3_sub(sphere_centre, closest);
    let dist_sq = v3_len_sq(diff);
    let sum_r = cap_r + sphere_r;
    if dist_sq >= sum_r * sum_r {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = if dist < f64::EPSILON {
        [0.0, 1.0, 0.0]
    } else {
        v3_normalize(diff)
    };
    let depth = sum_r - dist;
    let contact_pos = v3_add(closest, v3_scale(normal, cap_r - depth * 0.5));
    let local_a = v3_scale(normal, cap_r);
    let local_b = v3_scale(v3_neg(normal), sphere_r);

    let mut manifold = ContactManifold::new(body_cap, body_sphere, normal);
    manifold.add_point(ContactPoint::new(contact_pos, local_a, local_b, depth));
    Some(manifold)
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold generation: capsule vs capsule
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the closest points on two line segments `[p1,p2]` and `[p3,p4]`.
///
/// Returns `(s, t, c1, c2)` where `c1` and `c2` are the closest points.
pub fn closest_points_on_segments(
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    p4: [f64; 3],
) -> (f64, f64, [f64; 3], [f64; 3]) {
    let d1 = v3_sub(p2, p1);
    let d2 = v3_sub(p4, p3);
    let r = v3_sub(p1, p3);

    let a = v3_dot(d1, d1);
    let e = v3_dot(d2, d2);
    let f = v3_dot(d2, r);

    let (s, t) = if a < f64::EPSILON && e < f64::EPSILON {
        (0.0_f64, 0.0_f64)
    } else if a < f64::EPSILON {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = v3_dot(d1, r);
        if e < f64::EPSILON {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = v3_dot(d1, d2);
            let denom = a * e - b * b;
            let s = if denom.abs() > f64::EPSILON {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + f) / e;
            let (s, t) = if t < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t > 1.0 {
                (((b - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s, t)
            };
            (s, t)
        }
    };

    let c1 = v3_add(p1, v3_scale(d1, s));
    let c2 = v3_add(p3, v3_scale(d2, t));
    (s, t, c1, c2)
}

/// Generates a contact manifold between two capsules.
#[allow(clippy::too_many_arguments)]
pub fn manifold_capsule_capsule(
    body_a: usize,
    body_b: usize,
    a1: [f64; 3],
    a2: [f64; 3],
    ra: f64,
    b1: [f64; 3],
    b2: [f64; 3],
    rb: f64,
) -> Option<ContactManifold> {
    let (_s, _t, ca, cb) = closest_points_on_segments(a1, a2, b1, b2);
    let diff = v3_sub(cb, ca);
    let dist_sq = v3_len_sq(diff);
    let sum_r = ra + rb;
    if dist_sq >= sum_r * sum_r {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = if dist < f64::EPSILON {
        [0.0, 1.0, 0.0]
    } else {
        v3_normalize(diff)
    };
    let depth = sum_r - dist;
    let contact_pos = v3_add(ca, v3_scale(normal, ra - depth * 0.5));
    let local_a = v3_scale(normal, ra);
    let local_b = v3_scale(v3_neg(normal), rb);

    let mut manifold = ContactManifold::new(body_a, body_b, normal);
    manifold.add_point(ContactPoint::new(contact_pos, local_a, local_b, depth));
    Some(manifold)
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold merging
// ─────────────────────────────────────────────────────────────────────────────

/// Merges two manifolds between the same body pair.
///
/// Returns a new manifold containing the union of both point sets, reduced
/// to `MAX_MANIFOLD_POINTS` if necessary.
pub fn merge_manifolds(a: &ContactManifold, b: &ContactManifold) -> ContactManifold {
    debug_assert_eq!(a.body_a, b.body_a);
    debug_assert_eq!(a.body_b, b.body_b);

    let mut all_points: Vec<ContactPoint> = Vec::new();
    for p in a.active_points() {
        all_points.push(*p);
    }
    for p in b.active_points() {
        all_points.push(*p);
    }

    // Average normal
    let n = v3_normalize(v3_add(a.basis.normal, b.basis.normal));
    let mut merged = ContactManifold::new(a.body_a, a.body_b, n);

    for cp in reduce_to_4(&all_points) {
        merged.add_point(cp);
    }
    merged
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold cache (persistent contacts)
// ─────────────────────────────────────────────────────────────────────────────

/// Key for uniquely identifying a body pair in the cache.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ManifoldKey {
    /// Lower body index.
    pub a: usize,
    /// Higher body index.
    pub b: usize,
}

impl ManifoldKey {
    /// Creates a canonical key with `a <= b`.
    pub fn new(i: usize, j: usize) -> Self {
        if i <= j {
            Self { a: i, b: j }
        } else {
            Self { a: j, b: i }
        }
    }
}

/// Cache mapping body pairs to their current contact manifold.
///
/// Used for warm-starting and persistent contact tracking.
#[derive(Debug, Default)]
pub struct ManifoldCache {
    entries: std::collections::HashMap<ManifoldKey, ContactManifold>,
}

impl ManifoldCache {
    /// Creates an empty cache.
    pub fn new() -> Self {
        Self {
            entries: std::collections::HashMap::new(),
        }
    }

    /// Inserts or replaces the manifold for a body pair.
    pub fn insert(&mut self, manifold: ContactManifold) {
        let key = ManifoldKey::new(manifold.body_a, manifold.body_b);
        self.entries.insert(key, manifold);
    }

    /// Returns a reference to the manifold for the given body pair, if any.
    pub fn get(&self, a: usize, b: usize) -> Option<&ContactManifold> {
        self.entries.get(&ManifoldKey::new(a, b))
    }

    /// Returns a mutable reference to the manifold for the given body pair, if any.
    pub fn get_mut(&mut self, a: usize, b: usize) -> Option<&mut ContactManifold> {
        self.entries.get_mut(&ManifoldKey::new(a, b))
    }

    /// Removes the manifold for the given body pair.
    pub fn remove(&mut self, a: usize, b: usize) {
        self.entries.remove(&ManifoldKey::new(a, b));
    }

    /// Returns the number of cached manifolds.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Removes manifolds that have no active contact points.
    pub fn purge_empty(&mut self) {
        self.entries.retain(|_, v| !v.is_empty());
    }

    /// Increments `age` on all contact points in all manifolds.
    pub fn age_all(&mut self) {
        for m in self.entries.values_mut() {
            m.age_points();
        }
    }

    /// Removes contact points older than `max_age` from all manifolds.
    pub fn prune_all(&mut self, max_age: u32) {
        for m in self.entries.values_mut() {
            m.prune_old_points(max_age);
        }
        self.purge_empty();
    }

    /// Iterates over all cached manifolds.
    pub fn iter(&self) -> impl Iterator<Item = &ContactManifold> {
        self.entries.values()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact Jacobian
// ─────────────────────────────────────────────────────────────────────────────

/// The Jacobian rows for a single contact constraint.
///
/// For rigid bodies, the constraint velocity is:
/// `J_n * [v_a; ω_a; v_b; ω_b]ᵀ ≥ 0`
#[derive(Debug, Clone, Copy)]
pub struct ContactJacobian {
    /// Linear part for body A.
    pub lin_a: [f64; 3],
    /// Angular part for body A.
    pub ang_a: [f64; 3],
    /// Linear part for body B (negated normal).
    pub lin_b: [f64; 3],
    /// Angular part for body B.
    pub ang_b: [f64; 3],
}

impl ContactJacobian {
    /// Computes the normal contact Jacobian.
    ///
    /// `r_a` is the vector from body A's centre to the contact point.
    /// `r_b` is the vector from body B's centre to the contact point.
    /// `normal` points from B to A.
    pub fn normal(r_a: [f64; 3], r_b: [f64; 3], normal: [f64; 3]) -> Self {
        Self {
            lin_a: normal,
            ang_a: v3_cross(r_a, normal),
            lin_b: v3_neg(normal),
            ang_b: v3_neg(v3_cross(r_b, normal)),
        }
    }

    /// Computes a tangent contact Jacobian (for friction).
    pub fn tangent(r_a: [f64; 3], r_b: [f64; 3], tangent: [f64; 3]) -> Self {
        Self {
            lin_a: tangent,
            ang_a: v3_cross(r_a, tangent),
            lin_b: v3_neg(tangent),
            ang_b: v3_neg(v3_cross(r_b, tangent)),
        }
    }

    /// Evaluates the Jacobian for the given velocities.
    ///
    /// Returns the constraint velocity (scalar).
    pub fn evaluate(
        &self,
        v_a: [f64; 3],
        omega_a: [f64; 3],
        v_b: [f64; 3],
        omega_b: [f64; 3],
    ) -> f64 {
        v3_dot(self.lin_a, v_a)
            + v3_dot(self.ang_a, omega_a)
            + v3_dot(self.lin_b, v_b)
            + v3_dot(self.ang_b, omega_b)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold quality metrics
// ─────────────────────────────────────────────────────────────────────────────

/// Quality metrics for a contact manifold.
#[derive(Debug, Clone, Copy)]
pub struct ManifoldQuality {
    /// Average penetration depth across all contact points.
    pub avg_depth: f64,
    /// Maximum penetration depth.
    pub max_depth: f64,
    /// Spread of contact points (largest pairwise distance).
    pub spread: f64,
    /// Area of the contact polygon (0 for fewer than 3 points).
    pub area: f64,
    /// Number of active contact points.
    pub point_count: usize,
}

/// Computes quality metrics for a contact manifold.
pub fn manifold_quality(m: &ContactManifold) -> ManifoldQuality {
    let pts: Vec<&ContactPoint> = m.active_points().collect();
    let n = pts.len();
    if n == 0 {
        return ManifoldQuality {
            avg_depth: 0.0,
            max_depth: 0.0,
            spread: 0.0,
            area: 0.0,
            point_count: 0,
        };
    }

    let avg_depth = pts.iter().map(|p| p.depth).sum::<f64>() / n as f64;
    let max_depth = pts
        .iter()
        .map(|p| p.depth)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut spread = 0.0_f64;
    for i in 0..n {
        for j in i + 1..n {
            let d = v3_dist(pts[i].position, pts[j].position);
            if d > spread {
                spread = d;
            }
        }
    }

    let area = if n >= 3 {
        triangle_area(pts[0].position, pts[1].position, pts[2].position)
    } else {
        0.0
    };

    ManifoldQuality {
        avg_depth,
        max_depth,
        spread,
        area,
        point_count: n,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Normal computation helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the best-fit normal for a set of contact points using PCA.
///
/// Falls back to the input `hint_normal` if the point set is degenerate (< 3 points).
pub fn best_fit_normal(points: &[[f64; 3]], hint_normal: [f64; 3]) -> [f64; 3] {
    if points.len() < 3 {
        return v3_normalize(hint_normal);
    }
    // Centroid
    let n = points.len() as f64;
    let centroid = points.iter().fold([0.0_f64; 3], |acc, p| v3_add(acc, *p));
    let centroid = v3_scale(centroid, 1.0 / n);

    // Covariance matrix (upper triangle)
    let mut cov = [[0.0_f64; 3]; 3];
    for p in points {
        let d = v3_sub(*p, centroid);
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }

    // Use the smallest eigenvector as normal (power iteration on (I - cov/max))
    // Simple heuristic: cross product of first two independent difference vectors
    let d0 = v3_sub(points[1], points[0]);
    let d1 = v3_sub(points[2], points[0]);
    let cross = v3_cross(d0, d1);
    let candidate = v3_normalize(cross);

    // orient with hint
    if v3_dot(candidate, hint_normal) < 0.0 {
        v3_neg(candidate)
    } else {
        candidate
    }
}

/// Recomputes the manifold normal from the current contact points.
///
/// Useful after incrementally adding points.
pub fn recompute_manifold_normal(manifold: &mut ContactManifold, hint: [f64; 3]) {
    let positions: Vec<[f64; 3]> = manifold.active_points().map(|p| p.position).collect();
    if positions.len() >= 3 {
        let n = best_fit_normal(&positions, hint);
        manifold.update_normal(n);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Manifold splitting
// ─────────────────────────────────────────────────────────────────────────────

/// Splits a manifold into two groups based on a separation plane.
///
/// Contact points with `dot(position - plane_point, plane_normal) >= 0`
/// go into the "positive" manifold; the rest into the "negative" manifold.
pub fn split_manifold(
    m: &ContactManifold,
    plane_point: [f64; 3],
    plane_normal: [f64; 3],
) -> (ContactManifold, ContactManifold) {
    let mut pos = ContactManifold::new(m.body_a, m.body_b, m.basis.normal);
    let mut neg = ContactManifold::new(m.body_a, m.body_b, m.basis.normal);

    for cp in m.active_points() {
        let side = v3_dot(v3_sub(cp.position, plane_point), plane_normal);
        if side >= 0.0 {
            pos.add_point(*cp);
        } else {
            neg.add_point(*cp);
        }
    }
    (pos, neg)
}

// ─────────────────────────────────────────────────────────────────────────────
// Warm-starting helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Attempts to match the new manifold's contact points to cached warm-start data.
///
/// For each new point, finds the cached point whose `local_a` is within `radius`
/// and copies the accumulated impulses.
pub fn warm_start_manifold(
    new_manifold: &mut ContactManifold,
    cached: &ContactManifold,
    radius: f64,
) {
    for new_pt in new_manifold.points[..new_manifold.point_count]
        .iter_mut()
        .flatten()
    {
        // find best match in cached manifold
        if let Some(cached_pt) = cached.active_points().min_by(|a, b| {
            v3_dist(a.local_a, new_pt.local_a)
                .partial_cmp(&v3_dist(b.local_a, new_pt.local_a))
                .unwrap_or(std::cmp::Ordering::Equal)
        }) && v3_dist(cached_pt.local_a, new_pt.local_a) < radius
        {
            new_pt.normal_impulse = cached_pt.normal_impulse;
            new_pt.tangent_impulse0 = cached_pt.tangent_impulse0;
            new_pt.tangent_impulse1 = cached_pt.tangent_impulse1;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AABB-AABB contact manifold
// ─────────────────────────────────────────────────────────────────────────────

/// Generates a contact manifold between two axis-aligned boxes.
///
/// Uses the SAT (separating axis theorem) along the 3 coordinate axes.
/// Returns `None` if separated.
pub fn manifold_aabb_aabb(
    body_a: usize,
    body_b: usize,
    centre_a: [f64; 3],
    half_a: [f64; 3],
    centre_b: [f64; 3],
    half_b: [f64; 3],
) -> Option<ContactManifold> {
    let diff = v3_sub(centre_b, centre_a);

    let overlap = [
        (half_a[0] + half_b[0]) - diff[0].abs(),
        (half_a[1] + half_b[1]) - diff[1].abs(),
        (half_a[2] + half_b[2]) - diff[2].abs(),
    ];

    if overlap[0] <= 0.0 || overlap[1] <= 0.0 || overlap[2] <= 0.0 {
        return None;
    }

    // minimum overlap axis
    let axis = if overlap[0] <= overlap[1] && overlap[0] <= overlap[2] {
        0
    } else if overlap[1] <= overlap[2] {
        1
    } else {
        2
    };

    let depth = overlap[axis];
    let mut normal = [0.0_f64; 3];
    normal[axis] = if diff[axis] >= 0.0 { 1.0 } else { -1.0 };

    // contact point: face centre on B towards A
    let mut contact_pos = centre_b;
    contact_pos[axis] -= normal[axis] * half_b[axis];

    let local_a = v3_sub(contact_pos, centre_a);
    let local_b = v3_sub(contact_pos, centre_b);

    let mut manifold = ContactManifold::new(body_a, body_b, normal);
    manifold.add_point(ContactPoint::new(contact_pos, local_a, local_b, depth));
    Some(manifold)
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Vec3 helpers ─────────────────────────────────────────────────────────

    #[test]
    fn test_v3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let r = v3_add(a, b);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_v3_dot() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!((v3_dot(a, b)).abs() < 1e-12);
        assert!((v3_dot(a, a) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_v3_cross_perpendicular() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = v3_cross(a, b);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_v3_normalize_unit() {
        let a = [3.0, 4.0, 0.0];
        let n = v3_normalize(a);
        let len = v3_len(n);
        assert!((len - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_v3_normalize_zero_returns_fallback() {
        let n = v3_normalize([0.0, 0.0, 0.0]);
        assert_eq!(n, [1.0, 0.0, 0.0]);
    }

    // ── ContactBasis ─────────────────────────────────────────────────────────

    #[test]
    fn test_contact_basis_orthonormal_x() {
        let basis = ContactBasis::from_normal([1.0, 0.0, 0.0]);
        assert!(basis.is_orthonormal(1e-10), "basis not orthonormal");
    }

    #[test]
    fn test_contact_basis_orthonormal_y() {
        let basis = ContactBasis::from_normal([0.0, 1.0, 0.0]);
        assert!(basis.is_orthonormal(1e-10));
    }

    #[test]
    fn test_contact_basis_orthonormal_diagonal() {
        let n = v3_normalize([1.0, 1.0, 1.0]);
        let basis = ContactBasis::from_normal(n);
        assert!(basis.is_orthonormal(1e-10));
    }

    // ── ContactPoint ─────────────────────────────────────────────────────────

    #[test]
    fn test_contact_point_penetrating() {
        let cp = ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.05);
        assert!(cp.is_penetrating(0.01));
        assert!(!cp.is_penetrating(0.1));
    }

    #[test]
    fn test_contact_point_impulse_default_zero() {
        let cp = ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.1);
        assert_eq!(cp.normal_impulse, 0.0);
        assert_eq!(cp.tangent_impulse0, 0.0);
        assert_eq!(cp.tangent_impulse1, 0.0);
    }

    // ── ContactManifold ───────────────────────────────────────────────────────

    #[test]
    fn test_manifold_add_up_to_four_points() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        for k in 0..4 {
            let pos = [k as f64, 0.0, 0.0];
            m.add_point(ContactPoint::new(pos, pos, pos, 0.1));
        }
        assert_eq!(m.len(), 4);
    }

    #[test]
    fn test_manifold_clear() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.1));
        m.clear();
        assert!(m.is_empty());
    }

    #[test]
    fn test_manifold_deepest_point() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.1));
        m.add_point(ContactPoint::new([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.5));
        m.add_point(ContactPoint::new([2.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.2));
        let dp = m.deepest_point().unwrap();
        assert!((dp.depth - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_average_position() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        m.add_point(ContactPoint::new([2.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        let avg = m.average_position();
        assert!((avg[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_age_and_prune() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.1));
        m.age_points();
        m.age_points();
        m.prune_old_points(1); // age > 1 → pruned
        assert!(m.is_empty());
    }

    #[test]
    fn test_manifold_prune_keeps_young() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.1));
        m.prune_old_points(5);
        assert_eq!(m.len(), 1);
    }

    // ── reduce_to_4 ──────────────────────────────────────────────────────────

    #[test]
    fn test_reduce_to_4_fewer_unchanged() {
        let pts: Vec<ContactPoint> = (0..3)
            .map(|i| ContactPoint::new([i as f64, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1))
            .collect();
        let reduced = reduce_to_4(&pts);
        assert_eq!(reduced.len(), 3);
    }

    #[test]
    fn test_reduce_to_4_reduces_large_set() {
        let pts: Vec<ContactPoint> = (0..8)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::PI / 4.0;
                ContactPoint::new(
                    [angle.cos(), angle.sin(), 0.0],
                    [0.0; 3],
                    [0.0; 3],
                    0.1 + i as f64 * 0.01,
                )
            })
            .collect();
        let reduced = reduce_to_4(&pts);
        assert!(reduced.len() <= MAX_MANIFOLD_POINTS);
        assert!(!reduced.is_empty());
    }

    // ── triangle area ─────────────────────────────────────────────────────────

    #[test]
    fn test_triangle_area_unit() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let area = triangle_area(a, b, c);
        assert!((area - 0.5).abs() < 1e-10, "area={area}");
    }

    #[test]
    fn test_triangle_area_degenerate() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [2.0, 0.0, 0.0]; // collinear
        let area = triangle_area(a, b, c);
        assert!(area.abs() < 1e-12, "area={area}");
    }

    // ── sphere-sphere manifold ────────────────────────────────────────────────

    #[test]
    fn test_sphere_sphere_separated_returns_none() {
        let result = manifold_sphere_sphere(0, 1, [0.0; 3], 1.0, [5.0, 0.0, 0.0], 1.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_sphere_sphere_overlapping_returns_some() {
        let result = manifold_sphere_sphere(0, 1, [0.0; 3], 1.0, [1.5, 0.0, 0.0], 1.0);
        assert!(result.is_some());
        let m = result.unwrap();
        assert_eq!(m.len(), 1);
        let dp = m.deepest_point().unwrap();
        assert!(dp.depth > 0.0);
    }

    #[test]
    fn test_sphere_sphere_normal_unit_length() {
        let result = manifold_sphere_sphere(0, 1, [0.0; 3], 1.0, [1.5, 0.0, 0.0], 1.0).unwrap();
        let nl = v3_len(result.basis.normal);
        assert!((nl - 1.0).abs() < 1e-10, "normal length={nl}");
    }

    #[test]
    fn test_sphere_sphere_depth_correct() {
        let ra = 1.0_f64;
        let rb = 1.0_f64;
        let dist = 1.5_f64;
        let expected_depth = ra + rb - dist;
        let result = manifold_sphere_sphere(0, 1, [0.0; 3], ra, [dist, 0.0, 0.0], rb).unwrap();
        let dp = result.deepest_point().unwrap();
        assert!((dp.depth - expected_depth).abs() < 1e-10);
    }

    // ── sphere-AABB manifold ──────────────────────────────────────────────────

    #[test]
    fn test_sphere_aabb_outside_returns_none() {
        let result = manifold_sphere_aabb(0, 1, [5.0, 0.0, 0.0], 0.5, [0.0; 3], [1.0; 3]);
        assert!(result.is_none());
    }

    #[test]
    fn test_sphere_aabb_touching_returns_some() {
        // sphere at [2, 0, 0] with r=1, box at origin with half=1 → touching at x=1
        let result = manifold_sphere_aabb(0, 1, [1.9, 0.0, 0.0], 1.0, [0.0; 3], [1.0; 3]);
        assert!(result.is_some());
    }

    // ── capsule-sphere manifold ───────────────────────────────────────────────

    #[test]
    fn test_capsule_sphere_separated() {
        let result = manifold_capsule_sphere(
            0,
            1,
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.2,
            [5.0, 0.0, 0.0],
            0.2,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_capsule_sphere_overlapping() {
        let result = manifold_capsule_sphere(
            0,
            1,
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.0, 0.0, 0.0],
            0.5,
        );
        assert!(result.is_some());
        let m = result.unwrap();
        assert!(m.deepest_point().unwrap().depth > 0.0);
    }

    // ── closest point on segment ─────────────────────────────────────────────

    #[test]
    fn test_closest_point_midpoint() {
        let p = [0.0, 0.0, 0.0];
        let q = [2.0, 0.0, 0.0];
        let t = [1.0, 1.0, 0.0];
        let c = closest_point_on_segment(p, q, t);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!(c[1].abs() < 1e-12);
    }

    #[test]
    fn test_closest_point_clamp_start() {
        let p = [0.0, 0.0, 0.0];
        let q = [1.0, 0.0, 0.0];
        let t = [-1.0, 0.0, 0.0];
        let c = closest_point_on_segment(p, q, t);
        assert!((c[0]).abs() < 1e-12); // clamped to p
    }

    // ── capsule-capsule manifold ──────────────────────────────────────────────

    #[test]
    fn test_capsule_capsule_separated() {
        let result = manifold_capsule_capsule(
            0,
            1,
            [-5.0, 0.0, 0.0],
            [-3.0, 0.0, 0.0],
            0.2,
            [3.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            0.2,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_capsule_capsule_overlapping() {
        let result = manifold_capsule_capsule(
            0,
            1,
            [0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            0.5,
            [0.3, 0.0, 0.0],
            [0.3, 2.0, 0.0],
            0.5,
        );
        assert!(result.is_some());
    }

    // ── AABB-AABB manifold ────────────────────────────────────────────────────

    #[test]
    fn test_aabb_aabb_separated() {
        let result = manifold_aabb_aabb(0, 1, [0.0; 3], [1.0; 3], [5.0, 0.0, 0.0], [1.0; 3]);
        assert!(result.is_none());
    }

    #[test]
    fn test_aabb_aabb_overlapping() {
        let result = manifold_aabb_aabb(
            0,
            1,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.5, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_some());
        let m = result.unwrap();
        assert!(m.deepest_point().unwrap().depth > 0.0);
    }

    #[test]
    fn test_aabb_aabb_normal_axis_aligned() {
        let result = manifold_aabb_aabb(
            0,
            1,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.5, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        )
        .unwrap();
        // minimum overlap is along x
        let n = result.basis.normal;
        assert!(n[0].abs() > 0.9, "normal should be along x, got {:?}", n);
    }

    // ── Contact Jacobian ──────────────────────────────────────────────────────

    #[test]
    fn test_jacobian_normal_constraint_velocity() {
        let r_a = [0.0, -1.0, 0.0];
        let r_b = [0.0, 1.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let jac = ContactJacobian::normal(r_a, r_b, normal);
        // v_a = [0,1,0], v_b = [0,-1,0], omega = [0,0,0]
        let cv = jac.evaluate([0.0, 1.0, 0.0], [0.0; 3], [0.0, -1.0, 0.0], [0.0; 3]);
        assert!((cv - 2.0).abs() < 1e-10, "cv={cv}");
    }

    #[test]
    fn test_jacobian_tangent_non_zero() {
        let r_a = [0.0, 0.0, 0.0];
        let r_b = [0.0, 0.0, 0.0];
        let tangent = [1.0, 0.0, 0.0];
        let jac = ContactJacobian::tangent(r_a, r_b, tangent);
        let cv = jac.evaluate([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3], [0.0; 3]);
        assert!((cv - 1.0).abs() < 1e-10, "cv={cv}");
    }

    // ── Manifold quality ──────────────────────────────────────────────────────

    #[test]
    fn test_manifold_quality_empty() {
        let m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        let q = manifold_quality(&m);
        assert_eq!(q.point_count, 0);
        assert_eq!(q.avg_depth, 0.0);
    }

    #[test]
    fn test_manifold_quality_one_point() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.3));
        let q = manifold_quality(&m);
        assert!((q.avg_depth - 0.3).abs() < 1e-12);
        assert!((q.max_depth - 0.3).abs() < 1e-12);
        assert_eq!(q.spread, 0.0);
    }

    #[test]
    fn test_manifold_quality_spread() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        m.add_point(ContactPoint::new([3.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        let q = manifold_quality(&m);
        assert!((q.spread - 3.0).abs() < 1e-10);
    }

    // ── Manifold cache ────────────────────────────────────────────────────────

    #[test]
    fn test_manifold_cache_insert_and_get() {
        let mut cache = ManifoldCache::new();
        let m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        cache.insert(m);
        assert!(cache.get(0, 1).is_some());
        assert!(cache.get(1, 0).is_some()); // canonical key
    }

    #[test]
    fn test_manifold_cache_remove() {
        let mut cache = ManifoldCache::new();
        cache.insert(ContactManifold::new(0, 1, [0.0, 1.0, 0.0]));
        cache.remove(0, 1);
        assert!(cache.get(0, 1).is_none());
    }

    #[test]
    fn test_manifold_cache_purge_empty() {
        let mut cache = ManifoldCache::new();
        cache.insert(ContactManifold::new(0, 1, [0.0, 1.0, 0.0])); // no points → empty
        cache.insert(ContactManifold::new(2, 3, [0.0, 1.0, 0.0]));
        // add a point to 2-3
        if let Some(m) = cache.get_mut(2, 3) {
            m.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0; 3], 0.1));
        }
        cache.purge_empty();
        assert!(cache.get(0, 1).is_none());
        assert!(cache.get(2, 3).is_some());
    }

    // ── Merging manifolds ─────────────────────────────────────────────────────

    #[test]
    fn test_merge_manifolds_combines_points() {
        let mut m1 = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m1.add_point(ContactPoint::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        let mut m2 = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m2.add_point(ContactPoint::new([2.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.2));
        let merged = merge_manifolds(&m1, &m2);
        assert!(!merged.is_empty() && merged.len() <= MAX_MANIFOLD_POINTS);
    }

    // ── Manifold split ─────────────────────────────────────────────────────────

    #[test]
    fn test_split_manifold_separates_correctly() {
        let mut m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        m.add_point(ContactPoint::new([-1.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1));
        let (pos, neg) = split_manifold(&m, [0.0; 3], [1.0, 0.0, 0.0]);
        assert_eq!(pos.len(), 1);
        assert_eq!(neg.len(), 1);
    }

    // ── Warm starting ─────────────────────────────────────────────────────────

    #[test]
    fn test_warm_start_copies_impulse() {
        let mut cached = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        let mut cp = ContactPoint::new([0.0; 3], [0.1, 0.0, 0.0], [0.0; 3], 0.1);
        cp.normal_impulse = 5.0;
        cached.add_point(cp);

        let mut new_m = ContactManifold::new(0, 1, [0.0, 1.0, 0.0]);
        new_m.add_point(ContactPoint::new([0.0; 3], [0.1, 0.0, 0.0], [0.0; 3], 0.09));
        warm_start_manifold(&mut new_m, &cached, 0.05);
        let pt = new_m.active_points().next().unwrap();
        assert!((pt.normal_impulse - 5.0).abs() < 1e-12);
    }

    // ── best_fit_normal ──────────────────────────────────────────────────────

    #[test]
    fn test_best_fit_normal_xy_plane() {
        let pts = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let hint = [0.0, 0.0, 1.0];
        let n = best_fit_normal(&pts, hint);
        assert!(v3_dot(n, [0.0, 0.0, 1.0]).abs() > 0.99, "got {n:?}");
    }

    // ── ManifoldKey ──────────────────────────────────────────────────────────

    #[test]
    fn test_manifold_key_canonical() {
        let k1 = ManifoldKey::new(3, 5);
        let k2 = ManifoldKey::new(5, 3);
        assert_eq!(k1, k2);
    }

    // ── v3_lerp ──────────────────────────────────────────────────────────────

    #[test]
    fn test_v3_lerp_midpoint() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 4.0, 6.0];
        let m = v3_lerp(a, b, 0.5);
        assert!((m[0] - 1.0).abs() < 1e-12);
        assert!((m[1] - 2.0).abs() < 1e-12);
        assert!((m[2] - 3.0).abs() < 1e-12);
    }
}
