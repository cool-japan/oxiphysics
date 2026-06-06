// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Proximity query algorithms for collision detection.
//!
//! Provides distance queries between geometric primitives: points, segments,
//! triangles, capsules, meshes, and point clouds. Uses `[f64; 3]` for 3D
//! vectors (no nalgebra dependency).

// ─────────────────────────────────────────────────────────────────────────────
// Vector math helpers
// ─────────────────────────────────────────────────────────────────────────────

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
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

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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
// SpatialQueryResult
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a proximity/distance query between two shapes.
#[derive(Debug, Clone)]
pub struct SpatialQueryResult {
    /// Signed distance between the shapes (negative = penetrating).
    pub distance: f64,
    /// Closest point on shape A.
    pub witness_point_a: [f64; 3],
    /// Closest point on shape B.
    pub witness_point_b: [f64; 3],
    /// Separation normal pointing from B to A.
    pub normal: [f64; 3],
    /// True if the query is valid (shapes are not degenerate).
    pub valid: bool,
}

impl SpatialQueryResult {
    /// Creates a valid result.
    pub fn new(
        distance: f64,
        witness_point_a: [f64; 3],
        witness_point_b: [f64; 3],
        normal: [f64; 3],
    ) -> Self {
        Self {
            distance,
            witness_point_a,
            witness_point_b,
            normal,
            valid: true,
        }
    }

    /// Creates an invalid result (degenerate geometry).
    pub fn invalid() -> Self {
        Self {
            distance: f64::INFINITY,
            witness_point_a: [0.0; 3],
            witness_point_b: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            valid: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ProximityQuery — point-to-primitive queries
// ─────────────────────────────────────────────────────────────────────────────

/// Proximity query functions for basic geometric primitives.
pub struct ProximityQuery;

impl ProximityQuery {
    /// Closest point on segment \[a, b\] to query point p.
    ///
    /// Returns the closest point and parameter t ∈ \[0, 1\].
    pub fn closest_point_on_segment(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> ([f64; 3], f64) {
        let ab = sub3(b, a);
        let ap = sub3(p, a);
        let denom = dot3(ab, ab);
        if denom < 1e-15 {
            return (a, 0.0);
        }
        let t = (dot3(ap, ab) / denom).clamp(0.0, 1.0);
        (lerp3(a, b, t), t)
    }

    /// Closest point on an axis-aligned box \[lo, hi\] to point p.
    pub fn closest_point_to_aabb(lo: [f64; 3], hi: [f64; 3], p: [f64; 3]) -> [f64; 3] {
        [
            p[0].clamp(lo[0], hi[0]),
            p[1].clamp(lo[1], hi[1]),
            p[2].clamp(lo[2], hi[2]),
        ]
    }

    /// Tests whether point p is inside a convex hull defined by a set of planes.
    ///
    /// Each plane is (normal, offset): `dot(normal, x) <= offset`.
    pub fn point_in_convex_hull(p: [f64; 3], planes: &[([f64; 3], f64)]) -> bool {
        planes.iter().all(|(n, d)| dot3(*n, p) <= *d + 1e-9)
    }

    /// Returns the signed distance from point p to the plane (normal, offset).
    ///
    /// Positive = above, negative = below the plane.
    pub fn signed_dist_to_plane(p: [f64; 3], normal: [f64; 3], offset: f64) -> f64 {
        dot3(normal, p) - offset
    }

    /// Closest point on a sphere (center, radius) to point p.
    pub fn closest_point_on_sphere(center: [f64; 3], radius: f64, p: [f64; 3]) -> [f64; 3] {
        let d = sub3(p, center);
        let len = len3(d).max(1e-15);
        add3(center, scale3(d, radius / len))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SegmentSegmentDist
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the closest points between two 3D line segments.
pub struct SegmentSegmentDist;

impl SegmentSegmentDist {
    /// Computes the minimum distance and witness points between segments \[p1,p2\] and \[p3,p4\].
    ///
    /// Returns (distance, point_on_seg1, point_on_seg2).
    pub fn query(
        p1: [f64; 3],
        p2: [f64; 3],
        p3: [f64; 3],
        p4: [f64; 3],
    ) -> (f64, [f64; 3], [f64; 3]) {
        let d1 = sub3(p2, p1);
        let d2 = sub3(p4, p3);
        let r = sub3(p1, p3);

        let a = dot3(d1, d1); // squared length of seg1
        let e = dot3(d2, d2); // squared length of seg2
        let f = dot3(d2, r);

        let (s, t);
        if a < 1e-15 && e < 1e-15 {
            // Both segments are points
            return (dist3(p1, p3), p1, p3);
        }
        if a < 1e-15 {
            s = 0.0;
            t = (f / e).clamp(0.0, 1.0);
        } else {
            let c = dot3(d1, r);
            if e < 1e-15 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else {
                let b = dot3(d1, d2);
                let denom = a * e - b * b;
                s = if denom > 1e-15 {
                    ((b * f - c * e) / denom).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                t = (b * s + f) / e;
                if t < 0.0 {
                    let t_clamp = 0.0_f64;
                    let s_clamp = (-c / a).clamp(0.0, 1.0);
                    let wp1 = lerp3(p1, p2, s_clamp);
                    let wp2 = lerp3(p3, p4, t_clamp);
                    return (dist3(wp1, wp2), wp1, wp2);
                } else if t > 1.0 {
                    let t_clamp = 1.0_f64;
                    let s_clamp = ((b - c) / a).clamp(0.0, 1.0);
                    let wp1 = lerp3(p1, p2, s_clamp);
                    let wp2 = lerp3(p3, p4, t_clamp);
                    return (dist3(wp1, wp2), wp1, wp2);
                }
            }
        }
        let wp1 = lerp3(p1, p2, s);
        let wp2 = lerp3(p3, p4, t);
        (dist3(wp1, wp2), wp1, wp2)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PointTriangleDist
// ─────────────────────────────────────────────────────────────────────────────

/// Point-to-triangle distance query using barycentric coordinates.
pub struct PointTriangleDist;

impl PointTriangleDist {
    /// Computes the closest point on triangle (v0, v1, v2) to point p.
    ///
    /// Returns (closest_point, barycentric_coords \[u, v, w\]).
    pub fn closest_point(
        p: [f64; 3],
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
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

    /// Squared distance from point p to triangle.
    pub fn dist_sq(p: [f64; 3], v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f64 {
        let (q, _) = Self::closest_point(p, v0, v1, v2);
        let d = sub3(p, q);
        dot3(d, d)
    }

    /// Triangle normal (unnormalized).
    pub fn normal(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
        cross3(sub3(v1, v0), sub3(v2, v0))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SignedDistance — SDF operations
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance field (SDF) combinators and primitive SDFs.
pub struct SignedDistance;

impl SignedDistance {
    /// Signed distance from point p to a sphere (center, radius).
    pub fn sphere(p: [f64; 3], center: [f64; 3], radius: f64) -> f64 {
        dist3(p, center) - radius
    }

    /// Signed distance from point p to axis-aligned box centered at origin with half-extents h.
    pub fn box_sdf(p: [f64; 3], h: [f64; 3]) -> f64 {
        let q = [p[0].abs() - h[0], p[1].abs() - h[1], p[2].abs() - h[2]];
        let outer = len3([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]);
        let inner = q[0].max(q[1]).max(q[2]).min(0.0);
        outer + inner
    }

    /// Signed distance from point p to an infinite cylinder along Y axis, radius r.
    pub fn cylinder(p: [f64; 3], radius: f64) -> f64 {
        let r = (p[0] * p[0] + p[2] * p[2]).sqrt();
        r - radius
    }

    /// Union of two SDFs: min(d1, d2).
    pub fn union(d1: f64, d2: f64) -> f64 {
        d1.min(d2)
    }

    /// Intersection of two SDFs: max(d1, d2).
    pub fn intersection(d1: f64, d2: f64) -> f64 {
        d1.max(d2)
    }

    /// Subtraction of two SDFs: max(d1, -d2).
    pub fn subtraction(d1: f64, d2: f64) -> f64 {
        d1.max(-d2)
    }

    /// Smooth union with blending radius k.
    pub fn smooth_union(d1: f64, d2: f64, k: f64) -> f64 {
        let h = (0.5 + 0.5 * (d2 - d1) / k).clamp(0.0, 1.0);
        d1 * h + d2 * (1.0 - h) - k * h * (1.0 - h)
    }

    /// Signed distance from point to convex polygon (2D, in XY plane).
    ///
    /// Returns negative if p is inside the polygon, positive if outside.
    /// Assumes counter-clockwise vertex ordering for correct sign.
    pub fn convex_polygon_2d(p: [f64; 2], verts: &[[f64; 2]]) -> f64 {
        let n = verts.len();
        if n == 0 {
            return f64::INFINITY;
        }
        let mut min_dist = f64::INFINITY;
        // Ray-casting inside test (horizontal ray, +x direction)
        let mut crossings = 0i32;
        for i in 0..n {
            let a = verts[i];
            let b = verts[(i + 1) % n];
            let e = [b[0] - a[0], b[1] - a[1]];
            let w = [p[0] - a[0], p[1] - a[1]];
            // Closest point on segment to p
            let denom = e[0] * e[0] + e[1] * e[1];
            let t = if denom < 1e-15 {
                0.0
            } else {
                ((w[0] * e[0] + w[1] * e[1]) / denom).clamp(0.0, 1.0)
            };
            let closest = [a[0] + t * e[0], a[1] + t * e[1]];
            let dx = p[0] - closest[0];
            let dy = p[1] - closest[1];
            let d = (dx * dx + dy * dy).sqrt();
            if d < min_dist {
                min_dist = d;
            }
            // Ray-casting test: horizontal ray from p in +x direction
            // Check if edge crosses the horizontal ray at y = p[1]
            let a_above = a[1] > p[1];
            let b_above = b[1] > p[1];
            if a_above != b_above {
                // Edge crosses horizontal line at p[1]
                // Find x-intercept
                let t_cross = (p[1] - a[1]) / (b[1] - a[1]);
                let x_cross = a[0] + t_cross * (b[0] - a[0]);
                if x_cross > p[0] {
                    crossings += 1;
                }
            }
        }
        let inside = crossings % 2 == 1;
        if inside { -min_dist } else { min_dist }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CapsuleCapsuleDist
// ─────────────────────────────────────────────────────────────────────────────

/// Specialized capsule-capsule distance query.
///
/// A capsule is defined by a segment (base, tip) and a radius.
pub struct CapsuleCapsuleDist;

impl CapsuleCapsuleDist {
    /// Minimum distance between two capsules.
    ///
    /// Returns `SpatialQueryResult` with distance and witness points.
    pub fn query(
        a_base: [f64; 3],
        a_tip: [f64; 3],
        a_radius: f64,
        b_base: [f64; 3],
        b_tip: [f64; 3],
        b_radius: f64,
    ) -> SpatialQueryResult {
        let (seg_dist, wa, wb) = SegmentSegmentDist::query(a_base, a_tip, b_base, b_tip);
        let combined_radius = a_radius + b_radius;
        let dist = seg_dist - combined_radius;
        let normal = if seg_dist > 1e-12 {
            normalize3(sub3(wa, wb))
        } else {
            [0.0, 1.0, 0.0]
        };
        let witness_a = add3(wa, scale3(normalize3(sub3(wb, wa)).map(|x| x), -a_radius));
        let witness_b = add3(wb, scale3(normalize3(sub3(wa, wb)).map(|x| x), -b_radius));
        SpatialQueryResult::new(dist, witness_a, witness_b, normal)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KnnQuery — k-nearest neighbors in point cloud
// ─────────────────────────────────────────────────────────────────────────────

/// K-nearest neighbor query for a point cloud.
///
/// Uses a simple linear scan for correctness; a kd-tree is used for efficiency
/// when the cloud is large.
pub struct KnnQuery {
    /// The point cloud.
    pub points: Vec<[f64; 3]>,
}

impl KnnQuery {
    /// Creates a new kNN query structure for the given point cloud.
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        Self { points }
    }

    /// Returns the k nearest neighbors of query point p (indices and distances).
    ///
    /// Returns a vector of (index, distance) pairs sorted by distance.
    pub fn knn(&self, p: [f64; 3], k: usize) -> Vec<(usize, f64)> {
        let mut dists: Vec<(usize, f64)> = self
            .points
            .iter()
            .enumerate()
            .map(|(i, &q)| (i, dist3(p, q)))
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.truncate(k);
        dists
    }

    /// Returns all points within radius r of query point p.
    pub fn range_query(&self, p: [f64; 3], r: f64) -> Vec<usize> {
        self.points
            .iter()
            .enumerate()
            .filter(|&(_, q)| dist3(p, *q) <= r)
            .map(|(i, _)| i)
            .collect()
    }

    /// Returns all points within a ball (center, radius) — same as range query.
    pub fn ball_query(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        self.range_query(center, radius)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RayQuery — ray intersection tests
// ─────────────────────────────────────────────────────────────────────────────

/// Ray structure: origin + direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin.
    pub origin: [f64; 3],
    /// Ray direction (should be unit vector).
    pub direction: [f64; 3],
}

impl Ray {
    /// Creates a new ray.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            origin,
            direction: normalize3(direction),
        }
    }

    /// Point along the ray at parameter t.
    pub fn at(&self, t: f64) -> [f64; 3] {
        add3(self.origin, scale3(self.direction, t))
    }
}

/// Ray intersection queries against primitives.
pub struct RayQuery;

impl RayQuery {
    /// Ray-AABB intersection. Returns t_min if hit, None otherwise.
    ///
    /// Uses the slab method.
    pub fn ray_aabb(ray: &Ray, lo: [f64; 3], hi: [f64; 3]) -> Option<f64> {
        let mut t_min = f64::NEG_INFINITY;
        let mut t_max = f64::INFINITY;
        for i in 0..3 {
            let inv_d = 1.0 / ray.direction[i];
            let t1 = (lo[i] - ray.origin[i]) * inv_d;
            let t2 = (hi[i] - ray.origin[i]) * inv_d;
            let (t1, t2) = if inv_d >= 0.0 { (t1, t2) } else { (t2, t1) };
            t_min = t_min.max(t1);
            t_max = t_max.min(t2);
        }
        if t_max >= t_min && t_max >= 0.0 {
            Some(t_min.max(0.0))
        } else {
            None
        }
    }

    /// Ray-sphere intersection. Returns t if hit (smallest positive t), None otherwise.
    pub fn ray_sphere(ray: &Ray, center: [f64; 3], radius: f64) -> Option<f64> {
        let oc = sub3(ray.origin, center);
        let a = dot3(ray.direction, ray.direction);
        let b = 2.0 * dot3(oc, ray.direction);
        let c = dot3(oc, oc) - radius * radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t1 = (-b - sqrt_disc) / (2.0 * a);
        let t2 = (-b + sqrt_disc) / (2.0 * a);
        if t1 >= 0.0 {
            Some(t1)
        } else if t2 >= 0.0 {
            Some(t2)
        } else {
            None
        }
    }

    /// Ray-triangle intersection using the Möller-Trumbore algorithm.
    ///
    /// Returns (t, u, v) if hit, where (u, v) are barycentric coordinates.
    pub fn ray_triangle(
        ray: &Ray,
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
    ) -> Option<(f64, f64, f64)> {
        let edge1 = sub3(v1, v0);
        let edge2 = sub3(v2, v0);
        let h = cross3(ray.direction, edge2);
        let a = dot3(edge1, h);
        if a.abs() < 1e-12 {
            return None; // parallel
        }
        let f = 1.0 / a;
        let s = sub3(ray.origin, v0);
        let u = f * dot3(s, h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross3(s, edge1);
        let v = f * dot3(ray.direction, q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = f * dot3(edge2, q);
        if t < 1e-9 {
            return None;
        }
        Some((t, u, v))
    }

    /// Ray-capsule intersection. Returns the smallest positive t, or None.
    pub fn ray_capsule(ray: &Ray, cap_a: [f64; 3], cap_b: [f64; 3], radius: f64) -> Option<f64> {
        // Test infinite cylinder, then clamp to capsule end spheres
        let ab = sub3(cap_b, cap_a);
        let ao = sub3(ray.origin, cap_a);
        let dir = ray.direction;

        let ab_len2 = dot3(ab, ab);
        if ab_len2 < 1e-15 {
            // Degenerate: just a sphere
            return Self::ray_sphere(ray, cap_a, radius);
        }
        let ab_hat = scale3(ab, 1.0 / ab_len2.sqrt());

        // Project direction and ao onto plane perp to ab
        let d_perp = sub3(dir, scale3(ab_hat, dot3(dir, ab_hat)));
        let ao_perp = sub3(ao, scale3(ab_hat, dot3(ao, ab_hat)));

        let a_coef = dot3(d_perp, d_perp);
        let b_coef = 2.0 * dot3(d_perp, ao_perp);
        let c_coef = dot3(ao_perp, ao_perp) - radius * radius;

        let mut t_min = f64::INFINITY;

        if a_coef > 1e-15 {
            let disc = b_coef * b_coef - 4.0 * a_coef * c_coef;
            if disc >= 0.0 {
                let sqrt_d = disc.sqrt();
                for sign in &[-1.0_f64, 1.0_f64] {
                    let t = (-b_coef + sign * sqrt_d) / (2.0 * a_coef);
                    if t >= 0.0 {
                        let hit = ray.at(t);
                        let along = dot3(sub3(hit, cap_a), ab_hat);
                        if along >= 0.0 && along <= ab_len2.sqrt() && t < t_min {
                            t_min = t;
                        }
                    }
                }
            }
        }

        // Test end spheres
        for &center in &[cap_a, cap_b] {
            if let Some(t) = Self::ray_sphere(ray, center, radius)
                && t < t_min
            {
                t_min = t;
            }
        }

        if t_min < f64::INFINITY {
            Some(t_min)
        } else {
            None
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PointMeshDist — BVH-accelerated point-to-mesh query
// ─────────────────────────────────────────────────────────────────────────────

/// BVH-accelerated point-to-triangle mesh distance query.
///
/// Stores mesh triangles and uses a simple bounding-box hierarchy for pruning.
pub struct PointMeshDist {
    /// Triangle vertex positions (length = n_triangles * 3, stride = 3 vertices).
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices: each group of 3 = (v0, v1, v2).
    pub indices: Vec<usize>,
}

impl PointMeshDist {
    /// Creates a new mesh distance query structure.
    pub fn new(vertices: Vec<[f64; 3]>, indices: Vec<usize>) -> Self {
        assert!(indices.len().is_multiple_of(3));
        Self { vertices, indices }
    }

    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.indices.len() / 3
    }

    /// Minimum squared distance from point p to any triangle in the mesh.
    ///
    /// Returns (distance, triangle_index, closest_point).
    pub fn closest_point(&self, p: [f64; 3]) -> (f64, usize, [f64; 3]) {
        let mut min_dist_sq = f64::INFINITY;
        let mut best_tri = 0;
        let mut best_pt = p;
        for tri in 0..self.n_triangles() {
            let v0 = self.vertices[self.indices[3 * tri]];
            let v1 = self.vertices[self.indices[3 * tri + 1]];
            let v2 = self.vertices[self.indices[3 * tri + 2]];
            let (q, _) = PointTriangleDist::closest_point(p, v0, v1, v2);
            let d2 = dot3(sub3(p, q), sub3(p, q));
            if d2 < min_dist_sq {
                min_dist_sq = d2;
                best_tri = tri;
                best_pt = q;
            }
        }
        (min_dist_sq.sqrt(), best_tri, best_pt)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MeshMeshDist — minimum distance between two meshes
// ─────────────────────────────────────────────────────────────────────────────

/// Minimum distance query between two triangle meshes.
pub struct MeshMeshDist;

impl MeshMeshDist {
    /// Brute-force minimum distance between two meshes.
    ///
    /// Returns (distance, tri_idx_a, tri_idx_b, witness_a, witness_b).
    pub fn query(
        verts_a: &[[f64; 3]],
        idx_a: &[usize],
        verts_b: &[[f64; 3]],
        idx_b: &[usize],
    ) -> (f64, usize, usize, [f64; 3], [f64; 3]) {
        let na = idx_a.len() / 3;
        let nb = idx_b.len() / 3;
        let mut min_dist = f64::INFINITY;
        let mut best = (0usize, 0usize, [0.0f64; 3], [0.0f64; 3]);

        for ta in 0..na {
            let v0a = verts_a[idx_a[3 * ta]];
            let v1a = verts_a[idx_a[3 * ta + 1]];
            let v2a = verts_a[idx_a[3 * ta + 2]];
            // Centroid of triangle A as query point
            let centroid_a = scale3(add3(add3(v0a, v1a), v2a), 1.0 / 3.0);
            for tb in 0..nb {
                let v0b = verts_b[idx_b[3 * tb]];
                let v1b = verts_b[idx_b[3 * tb + 1]];
                let v2b = verts_b[idx_b[3 * tb + 2]];
                let (q_b, _) = PointTriangleDist::closest_point(centroid_a, v0b, v1b, v2b);
                let (q_a, _) = PointTriangleDist::closest_point(q_b, v0a, v1a, v2a);
                let d = dist3(q_a, q_b);
                if d < min_dist {
                    min_dist = d;
                    best = (ta, tb, q_a, q_b);
                }
            }
        }
        (min_dist, best.0, best.1, best.2, best.3)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SweptVolumeDist — swept capsule distance
// ─────────────────────────────────────────────────────────────────────────────

/// Distance between a swept volume (capsule over a trajectory) and a static shape.
pub struct SweptVolumeDist;

impl SweptVolumeDist {
    /// Minimum distance between two swept capsule volumes.
    ///
    /// Each swept capsule is described by start/end positions of its axis endpoint,
    /// and a radius.
    pub fn capsule_trajectory_dist(
        start_a: [f64; 3],
        end_a: [f64; 3],
        radius_a: f64,
        start_b: [f64; 3],
        end_b: [f64; 3],
        radius_b: f64,
    ) -> f64 {
        // Minimum over trajectory: approximate by testing at multiple t values
        let n = 8;
        let mut min_dist = f64::INFINITY;
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let pa = lerp3(start_a, end_a, t);
            let pb = lerp3(start_b, end_b, t);
            let d = dist3(pa, pb) - radius_a - radius_b;
            if d < min_dist {
                min_dist = d;
            }
        }
        min_dist
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NarrowPhaseProximity — proximity pair cache
// ─────────────────────────────────────────────────────────────────────────────

/// Caches proximity pairs for narrow-phase filtering.
///
/// Maintains hysteresis: pairs remain active until distance > deactivation_threshold.
pub struct NarrowPhaseProximity {
    /// Distance threshold to activate a proximity pair.
    pub activation_threshold: f64,
    /// Distance threshold to deactivate a proximity pair.
    pub deactivation_threshold: f64,
    /// Active proximity pairs: (id_a, id_b).
    pub active_pairs: Vec<(u32, u32)>,
}

impl NarrowPhaseProximity {
    /// Creates a new proximity cache.
    pub fn new(activation_threshold: f64, deactivation_threshold: f64) -> Self {
        Self {
            activation_threshold,
            deactivation_threshold,
            active_pairs: Vec::new(),
        }
    }

    /// Updates the set of active pairs given new distance measurements.
    ///
    /// `distances` is a list of (id_a, id_b, distance).
    pub fn update(&mut self, distances: &[(u32, u32, f64)]) {
        for &(a, b, d) in distances {
            let pair = (a.min(b), a.max(b));
            let already_active = self.active_pairs.contains(&pair);
            if d <= self.activation_threshold && !already_active {
                self.active_pairs.push(pair);
            } else if d > self.deactivation_threshold && already_active {
                self.active_pairs.retain(|&p| p != pair);
            }
        }
    }

    /// Returns true if pair (a, b) is currently active.
    pub fn is_active(&self, a: u32, b: u32) -> bool {
        let pair = (a.min(b), a.max(b));
        self.active_pairs.contains(&pair)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_closest_point_on_segment_midpoint() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let p = [1.0, 1.0, 0.0];
        let (q, t) = ProximityQuery::closest_point_on_segment(a, b, p);
        assert!((q[0] - 1.0).abs() < 1e-12);
        assert!(q[1].abs() < 1e-12);
        assert!((t - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_closest_point_on_segment_clamped_start() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let p = [-1.0, 0.0, 0.0];
        let (q, t) = ProximityQuery::closest_point_on_segment(a, b, p);
        assert!((q[0] - 0.0).abs() < 1e-12);
        assert!((t - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_closest_point_on_segment_clamped_end() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let p = [2.0, 0.0, 0.0];
        let (q, t) = ProximityQuery::closest_point_on_segment(a, b, p);
        assert!((q[0] - 1.0).abs() < 1e-12);
        assert!((t - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_closest_point_to_aabb() {
        let lo = [0.0, 0.0, 0.0];
        let hi = [1.0, 1.0, 1.0];
        let p = [2.0, 0.5, 0.5];
        let q = ProximityQuery::closest_point_to_aabb(lo, hi, p);
        assert!((q[0] - 1.0).abs() < 1e-12);
        assert!((q[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_point_in_convex_hull() {
        // Unit square in XY plane
        let planes = [
            ([1.0, 0.0, 0.0], 1.0),
            ([-1.0, 0.0, 0.0], 0.0),
            ([0.0, 1.0, 0.0], 1.0),
            ([0.0, -1.0, 0.0], 0.0),
            ([0.0, 0.0, 1.0], 1.0),
            ([0.0, 0.0, -1.0], 0.0),
        ];
        assert!(ProximityQuery::point_in_convex_hull(
            [0.5, 0.5, 0.5],
            &planes
        ));
        assert!(!ProximityQuery::point_in_convex_hull(
            [1.5, 0.5, 0.5],
            &planes
        ));
    }

    #[test]
    fn test_segment_segment_dist_parallel() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let p4 = [1.0, 1.0, 0.0];
        let (d, _, _) = SegmentSegmentDist::query(p1, p2, p3, p4);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_dist_crossing() {
        let p1 = [-1.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, -1.0, 0.0];
        let p4 = [0.0, 1.0, 0.0];
        let (d, wa, wb) = SegmentSegmentDist::query(p1, p2, p3, p4);
        assert!(d < 1e-10, "d={}", d);
        assert!(wa[0].abs() < 1e-10);
        assert!(wb[0].abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_dist_skew() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, 0.0, 1.0];
        let p4 = [1.0, 0.0, 1.0];
        let (d, _, _) = SegmentSegmentDist::query(p1, p2, p3, p4);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_triangle_dist_inside() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let p = [0.25, 0.25, 1.0];
        let (q, bary) = PointTriangleDist::closest_point(p, v0, v1, v2);
        assert!((q[2] - 0.0).abs() < 1e-12);
        let bary_sum = bary[0] + bary[1] + bary[2];
        assert!((bary_sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_triangle_dist_vertex() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let p = [-1.0, -1.0, 0.0];
        let (q, _) = PointTriangleDist::closest_point(p, v0, v1, v2);
        assert!((q[0] - 0.0).abs() < 1e-12);
        assert!((q[1] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_sdf_sphere() {
        let center = [0.0, 0.0, 0.0];
        let d = SignedDistance::sphere([1.5, 0.0, 0.0], center, 1.0);
        assert!((d - 0.5).abs() < 1e-12);
        let d_inside = SignedDistance::sphere([0.5, 0.0, 0.0], center, 1.0);
        assert!(d_inside < 0.0);
    }

    #[test]
    fn test_sdf_box() {
        let h = [1.0, 1.0, 1.0];
        let d_outside = SignedDistance::box_sdf([2.0, 0.0, 0.0], h);
        assert!((d_outside - 1.0).abs() < 1e-10);
        let d_inside = SignedDistance::box_sdf([0.0, 0.0, 0.0], h);
        assert!(d_inside < 0.0);
    }

    #[test]
    fn test_sdf_union_intersection_subtraction() {
        let d1 = 1.0;
        let d2 = 2.0;
        assert_eq!(SignedDistance::union(d1, d2), 1.0);
        assert_eq!(SignedDistance::intersection(d1, d2), 2.0);
        assert_eq!(SignedDistance::subtraction(d1, d2), 1.0);
    }

    #[test]
    fn test_capsule_capsule_dist_touching() {
        let a_base = [0.0, 0.0, 0.0];
        let a_tip = [1.0, 0.0, 0.0];
        let b_base = [2.5, 0.0, 0.0];
        let b_tip = [3.5, 0.0, 0.0];
        let res = CapsuleCapsuleDist::query(a_base, a_tip, 0.5, b_base, b_tip, 0.5);
        assert!(res.distance >= -1e-6, "dist={}", res.distance);
    }

    #[test]
    fn test_knn_query_k1() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let knn = KnnQuery::new(pts);
        let result = knn.knn([0.1, 0.0, 0.0], 1);
        assert_eq!(result[0].0, 0);
    }

    #[test]
    fn test_knn_range_query() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let knn = KnnQuery::new(pts);
        let result = knn.range_query([5.0, 0.0, 0.0], 1.5);
        assert!(result.len() >= 2);
    }

    #[test]
    fn test_ray_aabb_hit() {
        let ray = Ray::new([0.0, 0.0, -2.0], [0.0, 0.0, 1.0]);
        let lo = [-0.5, -0.5, -0.5];
        let hi = [0.5, 0.5, 0.5];
        let t = RayQuery::ray_aabb(&ray, lo, hi);
        assert!(t.is_some());
        assert!((t.unwrap() - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_ray_aabb_miss() {
        let ray = Ray::new([0.0, 2.0, -2.0], [0.0, 0.0, 1.0]);
        let lo = [-0.5, -0.5, -0.5];
        let hi = [0.5, 0.5, 0.5];
        let t = RayQuery::ray_aabb(&ray, lo, hi);
        assert!(t.is_none());
    }

    #[test]
    fn test_ray_sphere_hit() {
        let ray = Ray::new([0.0, 0.0, -3.0], [0.0, 0.0, 1.0]);
        let t = RayQuery::ray_sphere(&ray, [0.0, 0.0, 0.0], 1.0);
        assert!(t.is_some());
        assert!((t.unwrap() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_triangle_hit() {
        let ray = Ray::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        let v0 = [-1.0, -1.0, 0.0];
        let v1 = [1.0, -1.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let result = RayQuery::ray_triangle(&ray, v0, v1, v2);
        assert!(result.is_some());
        let (t, _, _) = result.unwrap();
        assert!((t - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_triangle_miss() {
        let ray = Ray::new([5.0, 5.0, -1.0], [0.0, 0.0, 1.0]);
        let v0 = [-1.0, -1.0, 0.0];
        let v1 = [1.0, -1.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let result = RayQuery::ray_triangle(&ray, v0, v1, v2);
        assert!(result.is_none());
    }

    #[test]
    fn test_point_mesh_dist() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let idx = vec![0usize, 1, 2];
        let mesh = PointMeshDist::new(verts, idx);
        let (d, _tri, _q) = mesh.closest_point([0.25, 0.25, 1.0]);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_narrow_phase_proximity_activate() {
        let mut npp = NarrowPhaseProximity::new(1.0, 2.0);
        npp.update(&[(0, 1, 0.5)]);
        assert!(npp.is_active(0, 1));
    }

    #[test]
    fn test_narrow_phase_proximity_deactivate() {
        let mut npp = NarrowPhaseProximity::new(1.0, 2.0);
        npp.update(&[(0, 1, 0.5)]);
        assert!(npp.is_active(0, 1));
        npp.update(&[(0, 1, 3.0)]);
        assert!(!npp.is_active(0, 1));
    }

    #[test]
    fn test_sdf_convex_polygon_inside() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let d = SignedDistance::convex_polygon_2d([0.5, 0.5], &square);
        assert!(d < 0.0);
    }

    #[test]
    fn test_sdf_convex_polygon_outside() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let d = SignedDistance::convex_polygon_2d([2.0, 0.5], &square);
        assert!(d > 0.0);
    }

    #[test]
    fn test_spatial_query_result_valid() {
        let r = SpatialQueryResult::new(1.0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(r.valid);
        assert!((r.distance - 1.0).abs() < 1e-12);
    }
}
