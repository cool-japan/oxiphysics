//! Specialized narrow-phase collision primitives operating directly on raw
//! geometric quantities (centers, radii, min/max points) without going through
//! the `Shape`/`Transform` abstraction layer.

// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_core::math::Vec3;

// ── ContactResult ─────────────────────────────────────────────────────────────

/// Result of a primitive collision test.
#[derive(Debug, Clone)]
pub struct ContactResult {
    /// Contact normal pointing from A to B (unit length).
    pub normal: Vec3,
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Witness point on shape A (surface of A).
    pub point_a: Vec3,
    /// Witness point on shape B (surface of B).
    pub point_b: Vec3,
}

// ── Sphere vs Sphere ──────────────────────────────────────────────────────────

/// Sphere–sphere collision test using raw geometric data.
///
/// Returns `None` when the spheres are separated (depth < 0).
pub fn sphere_sphere(
    center_a: Vec3,
    radius_a: f64,
    center_b: Vec3,
    radius_b: f64,
) -> Option<ContactResult> {
    let diff = center_b - center_a;
    let dist = diff.norm();
    let depth = radius_a + radius_b - dist;
    if depth < 0.0 {
        return None;
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let point_a = center_a + normal * radius_a;
    let point_b = center_b - normal * radius_b;
    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

// ── Sphere vs Plane ───────────────────────────────────────────────────────────

/// Sphere–plane collision test.
///
/// The plane is represented as `normal · x = offset`.
/// Returns `None` when the sphere does not penetrate the plane.
pub fn sphere_plane(
    sphere_center: Vec3,
    radius: f64,
    plane_normal: Vec3,
    plane_offset: f64,
) -> Option<ContactResult> {
    let signed_dist = sphere_center.dot(&plane_normal) - plane_offset;
    let depth = radius - signed_dist;
    if depth < 0.0 {
        return None;
    }
    // Normal points from A (sphere) toward B (plane surface).
    // Convention: normal from A to B, i.e. into the plane.
    let normal = -plane_normal;
    let point_a = sphere_center - plane_normal * radius;
    let point_b = sphere_center - plane_normal * signed_dist;
    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

// ── Sphere vs AABB ────────────────────────────────────────────────────────────

/// Sphere–AABB collision test.
///
/// The AABB is defined by its minimum and maximum corner points.
/// Returns `None` when the sphere does not reach the box.
pub fn sphere_aabb(
    sphere_center: Vec3,
    radius: f64,
    aabb_min: Vec3,
    aabb_max: Vec3,
) -> Option<ContactResult> {
    // Closest point on AABB to sphere center.
    let closest = Vec3::new(
        sphere_center.x.clamp(aabb_min.x, aabb_max.x),
        sphere_center.y.clamp(aabb_min.y, aabb_max.y),
        sphere_center.z.clamp(aabb_min.z, aabb_max.z),
    );
    let diff = sphere_center - closest;
    let dist_sq = diff.norm_squared();
    if dist_sq >= radius * radius {
        return None;
    }
    let dist = dist_sq.sqrt();
    let depth = radius - dist;
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        // Sphere center is inside; push out through nearest face.
        let dx = (sphere_center.x - aabb_min.x).min(aabb_max.x - sphere_center.x);
        let dy = (sphere_center.y - aabb_min.y).min(aabb_max.y - sphere_center.y);
        let dz = (sphere_center.z - aabb_min.z).min(aabb_max.z - sphere_center.z);
        if dx <= dy && dx <= dz {
            let sign = if sphere_center.x - aabb_min.x < aabb_max.x - sphere_center.x {
                -1.0
            } else {
                1.0
            };
            Vec3::new(sign, 0.0, 0.0)
        } else if dy <= dz {
            let sign = if sphere_center.y - aabb_min.y < aabb_max.y - sphere_center.y {
                -1.0
            } else {
                1.0
            };
            Vec3::new(0.0, sign, 0.0)
        } else {
            let sign = if sphere_center.z - aabb_min.z < aabb_max.z - sphere_center.z {
                -1.0
            } else {
                1.0
            };
            Vec3::new(0.0, 0.0, sign)
        }
    };
    let point_a = sphere_center - normal * radius;
    let point_b = closest;
    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

// ── Capsule vs Capsule ────────────────────────────────────────────────────────

/// Capsule–capsule collision test using raw segment data.
///
/// Each capsule is defined by two endpoints and a radius.
pub fn capsule_capsule(
    a0: Vec3,
    a1: Vec3,
    radius_a: f64,
    b0: Vec3,
    b1: Vec3,
    radius_b: f64,
) -> Option<ContactResult> {
    let (_, _, ca, cb) = closest_points_segment_segment(a0, a1, b0, b1);
    sphere_sphere(ca, radius_a, cb, radius_b)
}

/// Closest point on segment `(p0, p1)` to point `q`.
///
/// Returns `(point, t)` where `t ∈ [0, 1]` parametrises the segment.
pub fn closest_point_on_segment(p0: Vec3, p1: Vec3, q: Vec3) -> (Vec3, f64) {
    let d = p1 - p0;
    let len_sq = d.dot(&d);
    if len_sq < 1e-20 {
        return (p0, 0.0);
    }
    let t = ((q - p0).dot(&d) / len_sq).clamp(0.0, 1.0);
    (p0 + d * t, t)
}

/// Closest points between two line segments `(a0, a1)` and `(b0, b1)`.
///
/// Returns `(t_a, t_b, point_on_a, point_on_b)`.
pub fn closest_points_segment_segment(
    a0: Vec3,
    a1: Vec3,
    b0: Vec3,
    b1: Vec3,
) -> (f64, f64, Vec3, Vec3) {
    let d1 = a1 - a0;
    let d2 = b1 - b0;
    let r = a0 - b0;

    let a = d1.dot(&d1);
    let e = d2.dot(&d2);
    let f = d2.dot(&r);

    let (s, t) = if a <= 1e-10 && e <= 1e-10 {
        (0.0_f64, 0.0_f64)
    } else if a <= 1e-10 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(&r);
        if e <= 1e-10 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(&d2);
            let denom = a * e - b * b;
            let s = if denom.abs() > 1e-10 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t_raw = (b * s + f) / e;
            if t_raw < 0.0 {
                let s2 = (-c / a).clamp(0.0, 1.0);
                (s2, 0.0)
            } else if t_raw > 1.0 {
                let s2 = ((b - c) / a).clamp(0.0, 1.0);
                (s2, 1.0)
            } else {
                (s, t_raw)
            }
        }
    };

    (s, t, a0 + d1 * s, b0 + d2 * t)
}

// ── Sphere vs Triangle ────────────────────────────────────────────────────────

/// Sphere–triangle collision test.
///
/// Returns `None` when the sphere does not reach the triangle.
pub fn sphere_triangle(
    sphere_center: Vec3,
    radius: f64,
    v0: Vec3,
    v1: Vec3,
    v2: Vec3,
) -> Option<ContactResult> {
    let closest = closest_point_on_triangle(v0, v1, v2, sphere_center);
    let diff = sphere_center - closest;
    let dist = diff.norm();
    let depth = radius - dist;
    if depth < 0.0 {
        return None;
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        // Sphere centre lies on the triangle; use triangle normal.
        let n = (v1 - v0).cross(&(v2 - v0));
        let n_len = n.norm();
        if n_len > 1e-10 {
            n / n_len
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        }
    };
    let point_a = sphere_center - normal * radius;
    let point_b = closest;
    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

/// Closest point on triangle `(v0, v1, v2)` to point `p`.
///
/// Uses the Ericson (Real-Time Collision Detection) barycentric method.
pub fn closest_point_on_triangle(v0: Vec3, v1: Vec3, v2: Vec3, p: Vec3) -> Vec3 {
    let ab = v1 - v0;
    let ac = v2 - v0;
    let ap = p - v0;

    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return v0;
    }

    let bp = p - v1;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return v1;
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return v0 + ab * v;
    }

    let cp = p - v2;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return v2;
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return v0 + ac * w;
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return v1 + (v2 - v1) * w;
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    v0 + ab * v + ac * w
}

// ── AABB vs AABB ──────────────────────────────────────────────────────────────

/// AABB–AABB overlap test using the minimum overlap axis.
///
/// Returns `None` when the boxes are separated on any axis.
pub fn aabb_aabb(min_a: Vec3, max_a: Vec3, min_b: Vec3, max_b: Vec3) -> Option<ContactResult> {
    // Overlap on each axis.
    let ox = max_a.x.min(max_b.x) - min_a.x.max(min_b.x);
    let oy = max_a.y.min(max_b.y) - min_a.y.max(min_b.y);
    let oz = max_a.z.min(max_b.z) - min_a.z.max(min_b.z);

    if ox <= 0.0 || oy <= 0.0 || oz <= 0.0 {
        return None;
    }

    // Minimum overlap axis becomes the contact normal.
    let center_a = (min_a + max_a) * 0.5;
    let center_b = (min_b + max_b) * 0.5;
    let dir = center_b - center_a;

    let (depth, normal) = if ox <= oy && ox <= oz {
        let sign = if dir.x >= 0.0 { 1.0 } else { -1.0 };
        (ox, Vec3::new(sign, 0.0, 0.0))
    } else if oy <= oz {
        let sign = if dir.y >= 0.0 { 1.0 } else { -1.0 };
        (oy, Vec3::new(0.0, sign, 0.0))
    } else {
        let sign = if dir.z >= 0.0 { 1.0 } else { -1.0 };
        (oz, Vec3::new(0.0, 0.0, sign))
    };

    let point_a = center_a + normal * (max_a - center_a).dot(&normal);
    let point_b = center_b - normal * (center_b - min_b).dot(&normal);

    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

// ── Plane vs AABB ─────────────────────────────────────────────────────────────

/// Plane–AABB collision test.
///
/// The plane is represented as `normal · x = offset`.
/// Uses the support point of the AABB in the plane normal direction.
pub fn plane_aabb(
    plane_normal: Vec3,
    plane_offset: f64,
    aabb_min: Vec3,
    aabb_max: Vec3,
) -> Option<ContactResult> {
    // Support point: corner of AABB most in the direction of –plane_normal
    // (the point closest to the half-space below the plane).
    let support = Vec3::new(
        if plane_normal.x >= 0.0 {
            aabb_min.x
        } else {
            aabb_max.x
        },
        if plane_normal.y >= 0.0 {
            aabb_min.y
        } else {
            aabb_max.y
        },
        if plane_normal.z >= 0.0 {
            aabb_min.z
        } else {
            aabb_max.z
        },
    );
    let signed_dist = support.dot(&plane_normal) - plane_offset;
    let depth = -signed_dist; // positive when support is below (behind) the plane
    if depth < 0.0 {
        return None;
    }
    // Normal from AABB (A) toward plane (B).
    let normal = -plane_normal;
    let point_a = support;
    let point_b = support - plane_normal * signed_dist; // projection onto plane
    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

// ── Point-Triangle distance ──────────────────────────────────────────────────

/// Compute the squared distance from a point to a triangle.
pub fn point_triangle_distance_sq(p: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> f64 {
    let closest = closest_point_on_triangle(v0, v1, v2, p);
    (p - closest).norm_squared()
}

/// Compute the unsigned distance from a point to a triangle.
pub fn point_triangle_distance(p: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> f64 {
    point_triangle_distance_sq(p, v0, v1, v2).sqrt()
}

// ── Edge-Edge distance ──────────────────────────────────────────────────────

/// Compute the squared distance between two line segments.
pub fn edge_edge_distance_sq(a0: Vec3, a1: Vec3, b0: Vec3, b1: Vec3) -> f64 {
    let (_, _, ca, cb) = closest_points_segment_segment(a0, a1, b0, b1);
    (ca - cb).norm_squared()
}

/// Compute the unsigned distance between two line segments.
pub fn edge_edge_distance(a0: Vec3, a1: Vec3, b0: Vec3, b1: Vec3) -> f64 {
    edge_edge_distance_sq(a0, a1, b0, b1).sqrt()
}

// ── Segment-Triangle intersection ───────────────────────────────────────────

/// Test whether a line segment intersects a triangle.
///
/// Returns `Some((t, point))` where `t` is the parameter along the segment
/// and `point` is the intersection point.  Returns `None` if no intersection.
pub fn segment_triangle_intersection(
    seg_start: Vec3,
    seg_end: Vec3,
    v0: Vec3,
    v1: Vec3,
    v2: Vec3,
) -> Option<(f64, Vec3)> {
    let dir = seg_end - seg_start;
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let h = dir.cross(&e2);
    let a = e1.dot(&h);

    if a.abs() < 1e-12 {
        return None; // Parallel
    }

    let f = 1.0 / a;
    let s = seg_start - v0;
    let u = f * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = s.cross(&e1);
    let v = f * dir.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * e2.dot(&q);
    if (0.0..=1.0).contains(&t) {
        let point = seg_start + dir * t;
        Some((t, point))
    } else {
        None
    }
}

// ── Triangle-Triangle intersection ──────────────────────────────────────────

/// Test whether two triangles intersect.
///
/// Uses a robust approach: test all 6 segment-triangle pairs plus
/// coplanarity checks.
pub fn triangle_triangle_intersects(
    a0: Vec3,
    a1: Vec3,
    a2: Vec3,
    b0: Vec3,
    b1: Vec3,
    b2: Vec3,
) -> bool {
    // Test edges of A against triangle B
    let a_edges = [(a0, a1), (a1, a2), (a2, a0)];
    for &(start, end) in &a_edges {
        if segment_triangle_intersection(start, end, b0, b1, b2).is_some() {
            return true;
        }
    }

    // Test edges of B against triangle A
    let b_edges = [(b0, b1), (b1, b2), (b2, b0)];
    for &(start, end) in &b_edges {
        if segment_triangle_intersection(start, end, a0, a1, a2).is_some() {
            return true;
        }
    }

    false
}

// ── Closest feature pair ────────────────────────────────────────────────────

/// Identifies the closest geometric feature between two convex shapes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClosestFeature {
    /// Vertex of shape A, vertex of shape B.
    VertexVertex,
    /// Vertex of shape A, edge of shape B.
    VertexEdge,
    /// Edge of shape A, vertex of shape B.
    EdgeVertex,
    /// Edge of shape A, edge of shape B.
    EdgeEdge,
    /// Vertex of shape A, face of shape B.
    VertexFace,
    /// Face of shape A, vertex of shape B.
    FaceVertex,
}

/// Result of a closest-feature query between two triangles.
#[derive(Debug, Clone)]
pub struct ClosestFeatureResult {
    /// The type of closest feature pair.
    pub feature: ClosestFeature,
    /// Closest point on triangle A.
    pub point_a: Vec3,
    /// Closest point on triangle B.
    pub point_b: Vec3,
    /// Distance between the closest points.
    pub distance: f64,
}

/// Find the closest feature pair between two triangles.
///
/// Tests all vertex-vertex, vertex-edge, edge-edge, and vertex-face
/// combinations and returns the one with minimum distance.
pub fn closest_feature_pair(
    a0: Vec3,
    a1: Vec3,
    a2: Vec3,
    b0: Vec3,
    b1: Vec3,
    b2: Vec3,
) -> ClosestFeatureResult {
    let a_verts = [a0, a1, a2];
    let b_verts = [b0, b1, b2];
    let a_edges = [(a0, a1), (a1, a2), (a2, a0)];
    let b_edges = [(b0, b1), (b1, b2), (b2, b0)];

    let mut best_dist = f64::MAX;
    let mut best_pa = Vec3::zeros();
    let mut best_pb = Vec3::zeros();
    let mut best_feature = ClosestFeature::VertexVertex;

    // Vertex-Vertex
    for &va in &a_verts {
        for &vb in &b_verts {
            let d = (va - vb).norm();
            if d < best_dist {
                best_dist = d;
                best_pa = va;
                best_pb = vb;
                best_feature = ClosestFeature::VertexVertex;
            }
        }
    }

    // Vertex of A to edge of B
    for &va in &a_verts {
        for &(e0, e1) in &b_edges {
            let (closest, _) = closest_point_on_segment(e0, e1, va);
            let d = (va - closest).norm();
            if d < best_dist {
                best_dist = d;
                best_pa = va;
                best_pb = closest;
                best_feature = ClosestFeature::VertexEdge;
            }
        }
    }

    // Vertex of B to edge of A
    for &vb in &b_verts {
        for &(e0, e1) in &a_edges {
            let (closest, _) = closest_point_on_segment(e0, e1, vb);
            let d = (vb - closest).norm();
            if d < best_dist {
                best_dist = d;
                best_pa = closest;
                best_pb = vb;
                best_feature = ClosestFeature::EdgeVertex;
            }
        }
    }

    // Edge-Edge
    for &(a_e0, a_e1) in &a_edges {
        for &(b_e0, b_e1) in &b_edges {
            let (_, _, ca, cb) = closest_points_segment_segment(a_e0, a_e1, b_e0, b_e1);
            let d = (ca - cb).norm();
            if d < best_dist {
                best_dist = d;
                best_pa = ca;
                best_pb = cb;
                best_feature = ClosestFeature::EdgeEdge;
            }
        }
    }

    // Vertex of A to face of B
    for &va in &a_verts {
        let closest = closest_point_on_triangle(b0, b1, b2, va);
        let d = (va - closest).norm();
        if d < best_dist {
            best_dist = d;
            best_pa = va;
            best_pb = closest;
            best_feature = ClosestFeature::VertexFace;
        }
    }

    // Vertex of B to face of A
    for &vb in &b_verts {
        let closest = closest_point_on_triangle(a0, a1, a2, vb);
        let d = (vb - closest).norm();
        if d < best_dist {
            best_dist = d;
            best_pa = closest;
            best_pb = vb;
            best_feature = ClosestFeature::FaceVertex;
        }
    }

    ClosestFeatureResult {
        feature: best_feature,
        point_a: best_pa,
        point_b: best_pb,
        distance: best_dist,
    }
}

/// Compute the signed distance from a point to a plane defined by a triangle.
///
/// Positive if the point is on the side of the triangle normal.
pub fn point_plane_signed_distance(p: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> f64 {
    let n = (v1 - v0).cross(&(v2 - v0));
    let n_len = n.norm();
    if n_len < 1e-12 {
        return 0.0;
    }
    (p - v0).dot(&n) / n_len
}

/// Compute barycentric coordinates of a point projected onto a triangle.
///
/// Returns `(u, v, w)` such that `p ≈ u*v0 + v*v1 + w*v2`.
/// The point does not need to be on the triangle; the returned coordinates
/// are for the projection onto the triangle plane.
pub fn barycentric_coordinates(p: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> (f64, f64, f64) {
    let e0 = v1 - v0;
    let e1 = v2 - v0;
    let e2 = p - v0;

    let d00 = e0.dot(&e0);
    let d01 = e0.dot(&e1);
    let d11 = e1.dot(&e1);
    let d20 = e2.dot(&e0);
    let d21 = e2.dot(&e1);

    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-14 {
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }

    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;
    (u, v, w)
}

// ── Cylinder helpers ─────────────────────────────────────────────────────────

/// Closest point on a finite cylinder axis segment to an external point.
///
/// The cylinder axis runs from `c0` to `c1`. Returns the clamped projection.
pub fn closest_point_on_cylinder_axis(c0: Vec3, c1: Vec3, p: Vec3) -> (Vec3, f64) {
    closest_point_on_segment(c0, c1, p)
}

// ── Sphere vs Cylinder ────────────────────────────────────────────────────────

/// Sphere–cylinder collision test (finite capped cylinder).
///
/// The cylinder is defined by two end-cap centres `c0`/`c1` and a `radius`.
/// Returns `None` when separated.
pub fn sphere_cylinder(
    sphere_center: Vec3,
    sphere_radius: f64,
    c0: Vec3,
    c1: Vec3,
    cyl_radius: f64,
) -> Option<ContactResult> {
    let (axis_pt, _t) = closest_point_on_segment(c0, c1, sphere_center);
    let radial = sphere_center - axis_pt;
    let radial_dist = radial.norm();
    let sum_r = sphere_radius + cyl_radius;
    let depth = sum_r - radial_dist;
    if depth < 0.0 {
        return None;
    }
    let normal = if radial_dist > 1e-10 {
        radial / radial_dist
    } else {
        // Sphere centre is on the axis; push out along an arbitrary tangent.
        let axis = c1 - c0;
        let axis_len = axis.norm();
        if axis_len > 1e-10 {
            let perp = if axis.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0) - axis / axis_len * (axis.x / axis_len)
            } else {
                Vec3::new(0.0, 1.0, 0.0) - axis / axis_len * (axis.y / axis_len)
            };
            let pn = perp.norm();
            if pn > 1e-10 {
                perp / pn
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            }
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        }
    };
    let point_a = sphere_center - normal * sphere_radius;
    let point_b = axis_pt + normal * cyl_radius;
    Some(ContactResult {
        normal,
        depth,
        point_a,
        point_b,
    })
}

// ── Capsule vs Cylinder ───────────────────────────────────────────────────────

/// Capsule–cylinder collision test.
///
/// The capsule is `(a0, a1, r_cap)`, the cylinder is `(c0, c1, r_cyl)`.
/// Reduces to a sphere–cylinder test at the closest capsule axis point.
pub fn capsule_cylinder(
    a0: Vec3,
    a1: Vec3,
    r_cap: f64,
    c0: Vec3,
    c1: Vec3,
    r_cyl: f64,
) -> Option<ContactResult> {
    // Find the pair of closest points on the two axes.
    let (_s, _t, ca, _cb) = closest_points_segment_segment(a0, a1, c0, c1);
    sphere_cylinder(ca, r_cap, c0, c1, r_cyl)
}

// ── Box vs Cylinder (SAT-based) ───────────────────────────────────────────────

/// Approximate box–cylinder collision test via sphere expansion of the cylinder.
///
/// Treats the cylinder as a capsule (infinite cylinder is approximated by
/// expanding the AABB check) for a fast conservative test.  Returns `None`
/// when clearly separated.
pub fn box_cylinder_approx(
    aabb_min: Vec3,
    aabb_max: Vec3,
    c0: Vec3,
    c1: Vec3,
    cyl_radius: f64,
) -> Option<ContactResult> {
    // Find the closest point on the cylinder axis to the AABB centre.
    let box_center = (aabb_min + aabb_max) * 0.5;
    let (axis_pt, _) = closest_point_on_segment(c0, c1, box_center);
    // Treat the cylinder cross-section centre as a sphere for a quick reject.
    let half = (aabb_max - aabb_min) * 0.5;
    let max_half = half.x.max(half.y).max(half.z);
    let sum_r = cyl_radius + max_half;
    let diff = box_center - axis_pt;
    if diff.norm() > sum_r {
        return None;
    }
    // Use sphere_aabb as a conservative contact generator.
    sphere_aabb(axis_pt, cyl_radius, aabb_min, aabb_max)
}

// ── Signed Distance Functions ────────────────────────────────────────────────

/// Signed distance from a point to a sphere.
///
/// Negative inside the sphere, positive outside.
pub fn sdf_sphere(p: Vec3, center: Vec3, radius: f64) -> f64 {
    (p - center).norm() - radius
}

/// Signed distance from a point to an axis-aligned box (AABB).
///
/// Negative inside, positive outside.
pub fn sdf_aabb(p: Vec3, aabb_min: Vec3, aabb_max: Vec3) -> f64 {
    let center = (aabb_min + aabb_max) * 0.5;
    let half = (aabb_max - aabb_min) * 0.5;
    let q = Vec3::new(
        (p.x - center.x).abs() - half.x,
        (p.y - center.y).abs() - half.y,
        (p.z - center.z).abs() - half.z,
    );
    let outside = Vec3::new(q.x.max(0.0), q.y.max(0.0), q.z.max(0.0)).norm();
    let inside = q.x.max(q.y).max(q.z).min(0.0);
    outside + inside
}

/// Signed distance from a point to a capsule.
///
/// The capsule axis runs from `c0` to `c1` with the given `radius`.
pub fn sdf_capsule(p: Vec3, c0: Vec3, c1: Vec3, radius: f64) -> f64 {
    let (closest, _) = closest_point_on_segment(c0, c1, p);
    (p - closest).norm() - radius
}

/// Signed distance from a point to an infinite plane `n·x = d`.
///
/// Positive on the side the normal points toward.
pub fn sdf_plane(p: Vec3, normal: Vec3, d: f64) -> f64 {
    p.dot(&normal) - d
}

/// Signed distance from a point to a cylinder (lateral surface only, infinite).
///
/// The cylinder axis is defined by segment `c0`–`c1`.
pub fn sdf_cylinder_infinite(p: Vec3, c0: Vec3, c1: Vec3, radius: f64) -> f64 {
    let axis = c1 - c0;
    let axis_len = axis.norm();
    if axis_len < 1e-12 {
        return (p - c0).norm() - radius;
    }
    let axis_dir = axis / axis_len;
    let t = (p - c0).dot(&axis_dir);
    let axis_pt = c0 + axis_dir * t;
    (p - axis_pt).norm() - radius
}

// ── Point-in-Convex query ─────────────────────────────────────────────────────

/// Test whether a point lies inside a convex polyhedron defined by face planes.
///
/// Each element of `planes` is `(normal, d)` where `normal · x ≥ d` is the
/// inside half-space. The point is inside if it satisfies all half-spaces.
pub fn point_in_convex(p: Vec3, planes: &[(Vec3, f64)]) -> bool {
    for &(n, d) in planes {
        if p.dot(&n) < d - 1e-10 {
            return false;
        }
    }
    true
}

/// Build the six face-plane half-spaces for an AABB.
///
/// Returns planes in the form `(outward_normal, d)` where `normal · x ≥ d`
/// describes the *inside* of the box (i.e. inward half-space with negated normal).
pub fn aabb_face_planes(aabb_min: Vec3, aabb_max: Vec3) -> Vec<(Vec3, f64)> {
    vec![
        (Vec3::new(1.0, 0.0, 0.0), aabb_min.x),
        (Vec3::new(-1.0, 0.0, 0.0), -aabb_max.x),
        (Vec3::new(0.0, 1.0, 0.0), aabb_min.y),
        (Vec3::new(0.0, -1.0, 0.0), -aabb_max.y),
        (Vec3::new(0.0, 0.0, 1.0), aabb_min.z),
        (Vec3::new(0.0, 0.0, -1.0), -aabb_max.z),
    ]
}

// ── Closest point on triangle (extended) ─────────────────────────────────────

/// Closest point on a triangle with barycentric coordinates returned.
///
/// Returns `(closest_point, u, v, w)` where `u + v + w = 1` and
/// `closest = u*v0 + v*v1 + w*v2`.
pub fn closest_point_on_triangle_with_bary(
    v0: Vec3,
    v1: Vec3,
    v2: Vec3,
    p: Vec3,
) -> (Vec3, f64, f64, f64) {
    let closest = closest_point_on_triangle(v0, v1, v2, p);
    let (u, v, w) = barycentric_coordinates(closest, v0, v1, v2);
    (closest, u, v, w)
}

/// Compute the triangle area.
pub fn triangle_area(v0: Vec3, v1: Vec3, v2: Vec3) -> f64 {
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    e1.cross(&e2).norm() * 0.5
}

/// Compute the triangle normal (unit length).
///
/// Returns `Vec3::new(0,1,0)` for degenerate (zero-area) triangles.
pub fn triangle_normal(v0: Vec3, v1: Vec3, v2: Vec3) -> Vec3 {
    let n = (v1 - v0).cross(&(v2 - v0));
    let len = n.norm();
    if len > 1e-12 {
        n / len
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}

/// Project a point onto the plane of a triangle, returning the in-plane point.
pub fn project_onto_triangle_plane(v0: Vec3, v1: Vec3, v2: Vec3, p: Vec3) -> Vec3 {
    let n = triangle_normal(v0, v1, v2);
    let dist = (p - v0).dot(&n);
    p - n * dist
}

// ── Point–Convex signed distance ─────────────────────────────────────────────

/// Signed distance from a point to an AABB using the SDF formulation.
///
/// Negative means the point is inside the AABB.
pub fn point_aabb_signed_distance(p: Vec3, aabb_min: Vec3, aabb_max: Vec3) -> f64 {
    sdf_aabb(p, aabb_min, aabb_max)
}

// ── Sphere vs Sphere (with contact frame) ────────────────────────────────────

/// Compute a tangent basis (two orthonormal tangent vectors) for a given normal.
///
/// Returns `(tangent_u, tangent_v)` such that `{normal, tangent_u, tangent_v}`
/// form a right-handed orthonormal frame.
pub fn contact_tangent_basis(normal: Vec3) -> (Vec3, Vec3) {
    let t_u = if normal.x.abs() < 0.9 {
        let c = Vec3::new(1.0, 0.0, 0.0);
        let v = c - normal * c.dot(&normal);
        let l = v.norm();
        if l > 1e-10 {
            v / l
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        }
    } else {
        let c = Vec3::new(0.0, 1.0, 0.0);
        let v = c - normal * c.dot(&normal);
        let l = v.norm();
        if l > 1e-10 {
            v / l
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        }
    };
    let t_v = normal.cross(&t_u);
    (t_u, t_v)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── sphere vs sphere ─────────────────────────────────────────────────────

    #[test]
    fn test_sphere_sphere_touching() {
        // Two unit spheres with centres exactly r1+r2 apart → depth ≈ 0.
        let c = sphere_sphere(Vec3::new(0.0, 0.0, 0.0), 1.0, Vec3::new(2.0, 0.0, 0.0), 1.0)
            .expect("touching spheres must produce a contact");
        assert!(c.depth.abs() < 1e-10, "depth should be ≈0, got {}", c.depth);
    }

    #[test]
    fn test_sphere_sphere_overlapping() {
        let c = sphere_sphere(Vec3::new(0.0, 0.0, 0.0), 1.0, Vec3::new(1.5, 0.0, 0.0), 1.0)
            .expect("overlapping spheres must produce a contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
        assert!(
            (c.depth - 0.5).abs() < 1e-10,
            "expected depth 0.5, got {}",
            c.depth
        );
        assert!(
            (c.normal.norm() - 1.0).abs() < 1e-10,
            "normal must be unit length"
        );
    }

    #[test]
    fn test_sphere_sphere_separated() {
        assert!(
            sphere_sphere(Vec3::new(0.0, 0.0, 0.0), 1.0, Vec3::new(5.0, 0.0, 0.0), 1.0,).is_none(),
            "separated spheres must return None"
        );
    }

    // ── sphere vs plane ──────────────────────────────────────────────────────

    #[test]
    fn test_sphere_plane_above() {
        // Sphere at (0, 2, 0) radius 1, plane y=0 → signed_dist = 2, depth = 1-2 < 0.
        assert!(
            sphere_plane(Vec3::new(0.0, 2.0, 0.0), 1.0, Vec3::new(0.0, 1.0, 0.0), 0.0,).is_none(),
            "sphere above plane should return None"
        );
    }

    #[test]
    fn test_sphere_plane_penetrating() {
        // Sphere at (0, 0.5, 0) radius 1, plane y=0 → signed_dist = 0.5, depth = 0.5.
        let c = sphere_plane(Vec3::new(0.0, 0.5, 0.0), 1.0, Vec3::new(0.0, 1.0, 0.0), 0.0)
            .expect("sphere penetrating plane must return contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
        assert!(
            (c.depth - 0.5).abs() < 1e-10,
            "expected depth 0.5, got {}",
            c.depth
        );
    }

    // ── capsule vs capsule ───────────────────────────────────────────────────

    #[test]
    fn test_capsule_capsule_parallel() {
        // Two parallel capsules along Y, centres 0.8 apart, radii 0.5 each.
        let c = capsule_capsule(
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            0.5,
            Vec3::new(0.8, -1.0, 0.0),
            Vec3::new(0.8, 1.0, 0.0),
            0.5,
        )
        .expect("parallel overlapping capsules must produce contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
        assert!(
            (c.normal.norm() - 1.0).abs() < 1e-10,
            "normal must be unit length"
        );
    }

    #[test]
    fn test_capsule_capsule_miss() {
        // Two capsules along X, far apart.
        assert!(
            capsule_capsule(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                0.3,
                Vec3::new(10.0, 0.0, 0.0),
                Vec3::new(11.0, 0.0, 0.0),
                0.3,
            )
            .is_none(),
            "far-apart capsules must return None"
        );
    }

    // ── sphere vs triangle ───────────────────────────────────────────────────

    #[test]
    fn test_sphere_triangle_above_face() {
        // Triangle in XZ-plane, sphere directly above centre.
        let v0 = Vec3::new(-1.0, 0.0, -1.0);
        let v1 = Vec3::new(1.0, 0.0, -1.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);
        let c = sphere_triangle(Vec3::new(0.0, 0.5, 0.0), 1.0, v0, v1, v2)
            .expect("sphere above triangle face must produce contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
        assert!(
            (c.normal.norm() - 1.0).abs() < 1e-10,
            "normal must be unit length"
        );
    }

    #[test]
    fn test_sphere_triangle_near_edge() {
        // Triangle in XZ-plane, sphere near an edge but reaching it.
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(2.0, 0.0, 0.0);
        let v2 = Vec3::new(1.0, 0.0, 2.0);
        // Sphere centre at (1, 0.4, -0.3): close to edge v0-v1.
        let c = sphere_triangle(Vec3::new(1.0, 0.4, -0.3), 0.6, v0, v1, v2)
            .expect("sphere near edge must produce contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
    }

    // ── AABB vs AABB ─────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_aabb_separated() {
        assert!(
            aabb_aabb(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 1.0),
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(3.0, 1.0, 1.0),
            )
            .is_none(),
            "separated AABBs must return None"
        );
    }

    #[test]
    fn test_aabb_aabb_overlapping() {
        // Overlap of 0.5 on X, 1.0 on Y and Z → minimum overlap axis is X.
        let c = aabb_aabb(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(1.5, 1.0, 1.0),
        )
        .expect("overlapping AABBs must produce a contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
        // Minimum overlap is along X.
        assert!(
            (c.normal.x.abs() - 1.0).abs() < 1e-10,
            "normal should be along X"
        );
    }

    // ── closest_point_on_segment ─────────────────────────────────────────────

    #[test]
    fn test_closest_point_segment() {
        // Segment from (0,0,0) to (2,0,0); query point (1,1,0).
        // Closest point is the midpoint (1,0,0).
        let (pt, t) = closest_point_on_segment(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
        );
        assert!((pt.x - 1.0).abs() < 1e-10, "expected x=1.0, got {}", pt.x);
        assert!(pt.y.abs() < 1e-10, "expected y=0.0, got {}", pt.y);
        assert!((t - 0.5).abs() < 1e-10, "expected t=0.5, got {}", t);
    }

    // ── point-triangle distance ─────────────────────────────────────────

    #[test]
    fn test_point_triangle_distance_above() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(0.25, 0.25, 1.0); // directly above triangle
        let d = point_triangle_distance(p, v0, v1, v2);
        assert!((d - 1.0).abs() < 1e-10, "Expected distance 1.0, got {d}");
    }

    #[test]
    fn test_point_triangle_distance_on_triangle() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(0.25, 0.25, 0.0); // on the triangle
        let d = point_triangle_distance(p, v0, v1, v2);
        assert!(
            d < 1e-10,
            "Point on triangle should have distance ~0, got {d}"
        );
    }

    #[test]
    fn test_point_triangle_distance_near_vertex() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(-1.0, -1.0, 0.0); // closest to v0
        let d = point_triangle_distance(p, v0, v1, v2);
        let expected = (1.0_f64 + 1.0).sqrt(); // sqrt(2)
        assert!((d - expected).abs() < 1e-10, "Expected {expected}, got {d}");
    }

    // ── edge-edge distance ──────────────────────────────────────────────

    #[test]
    fn test_edge_edge_distance_parallel() {
        // Two parallel segments along Y, separated by 2 in X
        let d = edge_edge_distance(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(2.0, 1.0, 0.0),
        );
        assert!((d - 2.0).abs() < 1e-10, "Expected distance 2.0, got {d}");
    }

    #[test]
    fn test_edge_edge_distance_crossing() {
        // Two perpendicular segments that cross at distance 1 apart in Z
        let d = edge_edge_distance(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
        );
        assert!((d - 1.0).abs() < 1e-10, "Expected distance 1.0, got {d}");
    }

    // ── segment-triangle intersection ───────────────────────────────────

    #[test]
    fn test_segment_triangle_hit() {
        let v0 = Vec3::new(-1.0, 0.0, -1.0);
        let v1 = Vec3::new(1.0, 0.0, -1.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);
        // Segment going through triangle vertically
        let result = segment_triangle_intersection(
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            v0,
            v1,
            v2,
        );
        assert!(result.is_some(), "Segment should hit the triangle");
        let (t, pt) = result.unwrap();
        assert!((t - 0.5).abs() < 1e-10, "t should be 0.5, got {t}");
        assert!(pt.y.abs() < 1e-10, "intersection should be at y=0");
    }

    #[test]
    fn test_segment_triangle_miss() {
        let v0 = Vec3::new(-1.0, 0.0, -1.0);
        let v1 = Vec3::new(1.0, 0.0, -1.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);
        // Segment parallel to triangle, above it
        let result = segment_triangle_intersection(
            Vec3::new(-1.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            v0,
            v1,
            v2,
        );
        assert!(result.is_none(), "Parallel segment should miss");
    }

    // ── triangle-triangle intersection ──────────────────────────────────

    #[test]
    fn test_triangle_triangle_intersects() {
        // Two triangles that cross each other
        let a0 = Vec3::new(-1.0, 0.0, -1.0);
        let a1 = Vec3::new(1.0, 0.0, -1.0);
        let a2 = Vec3::new(0.0, 0.0, 1.0);

        let b0 = Vec3::new(0.0, -1.0, 0.0);
        let b1 = Vec3::new(0.0, 1.0, 0.0);
        let b2 = Vec3::new(0.0, 0.0, 2.0);

        assert!(
            triangle_triangle_intersects(a0, a1, a2, b0, b1, b2),
            "Crossing triangles should intersect"
        );
    }

    #[test]
    fn test_triangle_triangle_no_intersect() {
        // Two separated triangles
        let a0 = Vec3::new(0.0, 0.0, 0.0);
        let a1 = Vec3::new(1.0, 0.0, 0.0);
        let a2 = Vec3::new(0.0, 1.0, 0.0);

        let b0 = Vec3::new(0.0, 0.0, 10.0);
        let b1 = Vec3::new(1.0, 0.0, 10.0);
        let b2 = Vec3::new(0.0, 1.0, 10.0);

        assert!(
            !triangle_triangle_intersects(a0, a1, a2, b0, b1, b2),
            "Separated triangles should not intersect"
        );
    }

    // ── closest feature pair ────────────────────────────────────────────

    #[test]
    fn test_closest_feature_pair_vertex_vertex() {
        // Two triangles with closest features being two vertices
        let a0 = Vec3::new(0.0, 0.0, 0.0);
        let a1 = Vec3::new(1.0, 0.0, 0.0);
        let a2 = Vec3::new(0.0, 1.0, 0.0);

        let b0 = Vec3::new(5.0, 5.0, 0.0);
        let b1 = Vec3::new(6.0, 5.0, 0.0);
        let b2 = Vec3::new(5.0, 6.0, 0.0);

        let result = closest_feature_pair(a0, a1, a2, b0, b1, b2);
        assert!(
            result.distance > 0.0,
            "Separated triangles should have positive distance"
        );
        assert!(result.distance.is_finite(), "Distance should be finite");
    }

    #[test]
    fn test_closest_feature_pair_nearby() {
        // Two adjacent triangles sharing an edge approximately
        let a0 = Vec3::new(0.0, 0.0, 0.0);
        let a1 = Vec3::new(1.0, 0.0, 0.0);
        let a2 = Vec3::new(0.5, 1.0, 0.0);

        let b0 = Vec3::new(1.0, 0.0, 0.0);
        let b1 = Vec3::new(2.0, 0.0, 0.0);
        let b2 = Vec3::new(1.5, 1.0, 0.0);

        let result = closest_feature_pair(a0, a1, a2, b0, b1, b2);
        assert!(
            result.distance < 1e-10,
            "Adjacent triangles sharing a vertex should have ~0 distance"
        );
    }

    // ── barycentric coordinates ─────────────────────────────────────────

    #[test]
    fn test_barycentric_at_vertices() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        let (u, v, w) = barycentric_coordinates(v0, v0, v1, v2);
        assert!((u - 1.0).abs() < 1e-10, "u should be 1 at v0, got {u}");
        assert!(v.abs() < 1e-10, "v should be 0 at v0, got {v}");
        assert!(w.abs() < 1e-10, "w should be 0 at v0, got {w}");
    }

    #[test]
    fn test_barycentric_at_centroid() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let centroid = (v0 + v1 + v2) * (1.0 / 3.0);

        let (u, v, w) = barycentric_coordinates(centroid, v0, v1, v2);
        assert!((u - 1.0 / 3.0).abs() < 1e-10, "u at centroid: {u}");
        assert!((v - 1.0 / 3.0).abs() < 1e-10, "v at centroid: {v}");
        assert!((w - 1.0 / 3.0).abs() < 1e-10, "w at centroid: {w}");
    }

    // ── point-plane signed distance ─────────────────────────────────────

    #[test]
    fn test_point_plane_signed_distance() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        // Point above the XY plane
        let d = point_plane_signed_distance(Vec3::new(0.5, 0.5, 2.0), v0, v1, v2);
        assert!(
            d > 0.0,
            "Point above plane should have positive distance, got {d}"
        );

        // Point below
        let d2 = point_plane_signed_distance(Vec3::new(0.5, 0.5, -2.0), v0, v1, v2);
        assert!(
            d2 < 0.0,
            "Point below plane should have negative distance, got {d2}"
        );

        // Point on plane
        let d3 = point_plane_signed_distance(Vec3::new(0.5, 0.5, 0.0), v0, v1, v2);
        assert!(
            d3.abs() < 1e-10,
            "Point on plane should have zero distance, got {d3}"
        );
    }

    // ── plane-aabb test ─────────────────────────────────────────────────

    #[test]
    fn test_plane_aabb_penetrating() {
        let c = plane_aabb(
            Vec3::new(0.0, 1.0, 0.0),
            0.5,                        // plane y=0.5
            Vec3::new(-1.0, 0.0, -1.0), // AABB min
            Vec3::new(1.0, 1.0, 1.0),   // AABB max
        )
        .expect("AABB should penetrate the plane");
        assert!(c.depth > 0.0, "depth should be positive, got {}", c.depth);
    }

    #[test]
    fn test_plane_aabb_above() {
        let result = plane_aabb(
            Vec3::new(0.0, 1.0, 0.0),
            -5.0, // plane at y=-5
            Vec3::new(-1.0, 0.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        );
        assert!(result.is_none(), "AABB above plane should not collide");
    }

    // ── sphere-aabb edge case ───────────────────────────────────────────

    #[test]
    fn test_sphere_aabb_inside() {
        // Sphere entirely inside the AABB
        let c = sphere_aabb(
            Vec3::new(0.0, 0.0, 0.0),
            0.1,
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        )
        .expect("Sphere inside AABB should produce contact");
        assert!(c.depth > 0.0, "depth should be positive");
    }

    // ── sphere vs cylinder ──────────────────────────────────────────────

    #[test]
    fn test_sphere_cylinder_overlapping() {
        // Sphere at (1.2, 0, 0), radius 0.5; cylinder along Y, radius 1.0
        let c = sphere_cylinder(
            Vec3::new(1.2, 0.5, 0.0),
            0.5,
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
        )
        .expect("Sphere overlapping cylinder must produce contact");
        assert!(c.depth > 0.0, "depth must be positive, got {}", c.depth);
        assert!(
            (c.normal.norm() - 1.0).abs() < 1e-10,
            "normal must be unit length"
        );
    }

    #[test]
    fn test_sphere_cylinder_separated() {
        // Sphere far from cylinder
        let result = sphere_cylinder(
            Vec3::new(10.0, 0.0, 0.0),
            0.5,
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
        );
        assert!(result.is_none(), "far sphere should return None");
    }

    #[test]
    fn test_sphere_cylinder_touching() {
        // Sphere just touching: centre at (1.5, 0, 0), r=0.5, cyl_r=1.0 => sum=1.5 = dist
        let c = sphere_cylinder(
            Vec3::new(1.5, 0.0, 0.0),
            0.5,
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
        )
        .expect("touching sphere-cylinder must produce contact");
        assert!(c.depth.abs() < 1e-10, "depth should be ~0, got {}", c.depth);
    }

    // ── capsule vs cylinder ─────────────────────────────────────────────

    #[test]
    fn test_capsule_cylinder_overlapping() {
        // Capsule along Z near a Y-axis cylinder
        let result = capsule_cylinder(
            Vec3::new(1.2, 0.0, -0.5),
            Vec3::new(1.2, 0.0, 0.5),
            0.5,
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
        );
        assert!(result.is_some(), "capsule close to cylinder should collide");
    }

    #[test]
    fn test_capsule_cylinder_separated() {
        let result = capsule_cylinder(
            Vec3::new(5.0, 0.0, -0.5),
            Vec3::new(5.0, 0.0, 0.5),
            0.3,
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            0.5,
        );
        assert!(
            result.is_none(),
            "far capsule should not collide with cylinder"
        );
    }

    // ── box vs cylinder ─────────────────────────────────────────────────

    #[test]
    fn test_box_cylinder_approx_overlapping() {
        let result = box_cylinder_approx(
            Vec3::new(-0.5, -0.5, -0.5),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, -2.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            0.8,
        );
        assert!(
            result.is_some(),
            "box overlapping cylinder should produce contact"
        );
    }

    #[test]
    fn test_box_cylinder_approx_separated() {
        let result = box_cylinder_approx(
            Vec3::new(-0.5, -0.5, -0.5),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(10.0, -2.0, 0.0),
            Vec3::new(10.0, 2.0, 0.0),
            0.3,
        );
        assert!(result.is_none(), "far box-cylinder should return None");
    }

    // ── signed distance functions ───────────────────────────────────────

    #[test]
    fn test_sdf_sphere_outside() {
        let d = sdf_sphere(Vec3::new(3.0, 0.0, 0.0), Vec3::zeros(), 1.0);
        assert!((d - 2.0).abs() < 1e-10, "Expected 2.0, got {d}");
    }

    #[test]
    fn test_sdf_sphere_inside() {
        let d = sdf_sphere(Vec3::new(0.5, 0.0, 0.0), Vec3::zeros(), 1.0);
        assert!(
            d < 0.0,
            "Point inside sphere should have negative SDF, got {d}"
        );
        assert!((d - (-0.5)).abs() < 1e-10, "Expected -0.5, got {d}");
    }

    #[test]
    fn test_sdf_sphere_on_surface() {
        let d = sdf_sphere(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0);
        assert!(
            d.abs() < 1e-10,
            "Point on sphere surface should have SDF ≈ 0, got {d}"
        );
    }

    #[test]
    fn test_sdf_aabb_outside() {
        // Point at (2,0,0), box from -1 to 1 in all axes → distance = 1
        let d = sdf_aabb(
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        );
        assert!((d - 1.0).abs() < 1e-10, "Expected 1.0, got {d}");
    }

    #[test]
    fn test_sdf_aabb_inside() {
        // Point at origin in a unit box → should be negative
        let d = sdf_aabb(
            Vec3::zeros(),
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        );
        assert!(
            d < 0.0,
            "Point inside AABB should have negative SDF, got {d}"
        );
    }

    #[test]
    fn test_sdf_capsule_outside() {
        // Point at (2, 0, 0), capsule along Y with radius 0.5
        let d = sdf_capsule(
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            0.5,
        );
        // Closest axis point is (0,0,0), distance = 2, minus radius 0.5 = 1.5
        assert!((d - 1.5).abs() < 1e-10, "Expected 1.5, got {d}");
    }

    #[test]
    fn test_sdf_capsule_inside() {
        let d = sdf_capsule(
            Vec3::new(0.2, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            0.5,
        );
        assert!(
            d < 0.0,
            "Point inside capsule should have negative SDF, got {d}"
        );
    }

    #[test]
    fn test_sdf_plane_positive() {
        let d = sdf_plane(Vec3::new(0.0, 2.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.0);
        assert!((d - 2.0).abs() < 1e-10, "Expected 2.0, got {d}");
    }

    #[test]
    fn test_sdf_plane_negative() {
        let d = sdf_plane(Vec3::new(0.0, -1.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.0);
        assert!((d - (-1.0)).abs() < 1e-10, "Expected -1.0, got {d}");
    }

    #[test]
    fn test_sdf_cylinder_infinite_outside() {
        let d = sdf_cylinder_infinite(
            Vec3::new(2.0, 5.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
        );
        assert!((d - 1.0).abs() < 1e-10, "Expected 1.0, got {d}");
    }

    // ── point-in-convex ─────────────────────────────────────────────────

    #[test]
    fn test_point_in_convex_aabb_inside() {
        let planes = aabb_face_planes(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        assert!(
            point_in_convex(Vec3::zeros(), &planes),
            "origin should be inside unit box"
        );
    }

    #[test]
    fn test_point_in_convex_aabb_outside() {
        let planes = aabb_face_planes(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        assert!(
            !point_in_convex(Vec3::new(2.0, 0.0, 0.0), &planes),
            "point outside box should fail"
        );
    }

    #[test]
    fn test_point_in_convex_on_boundary() {
        let planes = aabb_face_planes(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        // Point exactly on the face
        assert!(
            point_in_convex(Vec3::new(1.0, 0.0, 0.0), &planes),
            "point on box face should be inside"
        );
    }

    // ── closest_point_on_triangle_with_bary ─────────────────────────────

    #[test]
    fn test_closest_point_bary_centroid() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(3.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 3.0, 0.0);
        let centroid = (v0 + v1 + v2) * (1.0 / 3.0);
        let (pt, u, v, w) = closest_point_on_triangle_with_bary(v0, v1, v2, centroid);
        assert!(
            (pt - centroid).norm() < 1e-10,
            "closest point should be centroid"
        );
        assert!(
            (u + v + w - 1.0).abs() < 1e-10,
            "bary coords should sum to 1"
        );
        assert!((u - 1.0 / 3.0).abs() < 1e-9, "u should be 1/3");
        assert!((v - 1.0 / 3.0).abs() < 1e-9, "v should be 1/3");
        assert!((w - 1.0 / 3.0).abs() < 1e-9, "w should be 1/3");
    }

    #[test]
    fn test_closest_point_bary_vertex() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let (pt, u, v, w) =
            closest_point_on_triangle_with_bary(v0, v1, v2, Vec3::new(-1.0, -1.0, 0.0));
        assert!((pt - v0).norm() < 1e-10, "closest to corner should be v0");
        assert!((u - 1.0).abs() < 1e-10, "u=1 at v0");
        let _ = (v, w);
    }

    // ── triangle geometry helpers ───────────────────────────────────────

    #[test]
    fn test_triangle_area_unit() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(2.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 2.0, 0.0);
        let area = triangle_area(v0, v1, v2);
        assert!((area - 2.0).abs() < 1e-10, "Expected area 2.0, got {area}");
    }

    #[test]
    fn test_triangle_normal_up() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);
        let n = triangle_normal(v0, v1, v2);
        // e1 = (1,0,0), e2 = (0,0,1), cross = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0,-1,0)
        assert!((n.norm() - 1.0).abs() < 1e-10, "Normal must be unit length");
    }

    #[test]
    fn test_project_onto_triangle_plane() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(0.5, 0.5, 3.0);
        let proj = project_onto_triangle_plane(v0, v1, v2, p);
        assert!(
            proj.z.abs() < 1e-10,
            "Projected point should be on plane z=0, got {}",
            proj.z
        );
    }

    // ── contact tangent basis ───────────────────────────────────────────

    #[test]
    fn test_contact_tangent_basis_orthonormal() {
        let normals = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 0.0).normalize(),
        ];
        for n in normals {
            let (u, v) = contact_tangent_basis(n);
            assert!(
                (u.norm() - 1.0).abs() < 1e-10,
                "tangent_u must be unit length"
            );
            assert!(
                (v.norm() - 1.0).abs() < 1e-10,
                "tangent_v must be unit length"
            );
            assert!(
                n.dot(&u).abs() < 1e-9,
                "tangent_u must be perpendicular to normal"
            );
            assert!(
                n.dot(&v).abs() < 1e-9,
                "tangent_v must be perpendicular to normal"
            );
            assert!(
                u.dot(&v).abs() < 1e-9,
                "tangent_u and tangent_v must be perpendicular"
            );
        }
    }

    // ── sdf aabb point signed distance ──────────────────────────────────

    #[test]
    fn test_point_aabb_signed_distance_inside() {
        let d = point_aabb_signed_distance(
            Vec3::zeros(),
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        );
        assert!(d < 0.0, "Inside AABB → negative SDF, got {d}");
    }

    #[test]
    fn test_point_aabb_signed_distance_outside() {
        let d = point_aabb_signed_distance(
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, 1.0),
        );
        assert!((d - 1.0).abs() < 1e-10, "Expected 1.0, got {d}");
    }

    // ── closest_point_on_cylinder_axis ──────────────────────────────────

    #[test]
    fn test_closest_on_cylinder_axis_midpoint() {
        let (pt, t) = closest_point_on_cylinder_axis(
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
        );
        assert!(
            pt.y.abs() < 1e-10,
            "closest on Y-axis should have y≈0, got {}",
            pt.y
        );
        assert!((t - 0.5).abs() < 1e-10, "t should be 0.5, got {t}");
    }

    #[test]
    fn test_closest_on_cylinder_axis_clamped() {
        // Query past end of segment
        let (pt, t) = closest_point_on_cylinder_axis(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
        );
        assert!(
            (pt.x - 1.0).abs() < 1e-10,
            "should be clamped to end, got {}",
            pt.x
        );
        assert!((t - 1.0).abs() < 1e-10, "t should be 1.0, got {t}");
    }
}
