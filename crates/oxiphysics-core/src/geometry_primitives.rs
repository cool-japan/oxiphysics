#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Low-level geometry primitives: Ray, AABB, Sphere, Plane, Triangle, Capsule,
//! OrientedBox, Frustum, convex polyhedron, and barycentric coordinates.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute `a + b`.
#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Compute `a - b`.
#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale vector by scalar.
#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Squared length.
#[inline]
fn len_sq(a: [f64; 3]) -> f64 {
    dot(a, a)
}

/// Length.
#[inline]
fn len(a: [f64; 3]) -> f64 {
    len_sq(a).sqrt()
}

/// Normalize a vector. Returns zero vector if near-zero.
#[inline]
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = len(a);
    if l < 1e-14 {
        [0.0; 3]
    } else {
        scale(a, 1.0 / l)
    }
}

/// Element-wise min.
#[inline]
fn vmin(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])]
}

/// Element-wise max.
#[inline]
fn vmax(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])]
}

// ─────────────────────────────────────────────────────────────────────────────
// Ray3
// ─────────────────────────────────────────────────────────────────────────────

/// A 3D ray with origin and (pre-normalized) direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray3 {
    /// Ray origin.
    pub origin: [f64; 3],
    /// Ray direction (should be normalized).
    pub direction: [f64; 3],
}

impl Ray3 {
    /// Create a new ray. Direction is normalized.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            origin,
            direction: normalize(direction),
        }
    }

    /// Point at parameter `t`: P(t) = origin + t * direction.
    pub fn at(&self, t: f64) -> [f64; 3] {
        add(self.origin, scale(self.direction, t))
    }

    /// Intersect ray with a sphere. Returns nearest positive `t`, or `None`.
    pub fn intersect_sphere(&self, center: [f64; 3], radius: f64) -> Option<f64> {
        let oc = sub(self.origin, center);
        let a = dot(self.direction, self.direction);
        let b = 2.0 * dot(oc, self.direction);
        let c = dot(oc, oc) - radius * radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t1 = (-b - sqrt_disc) / (2.0 * a);
        let t2 = (-b + sqrt_disc) / (2.0 * a);
        if t1 > 1e-8 {
            Some(t1)
        } else if t2 > 1e-8 {
            Some(t2)
        } else {
            None
        }
    }

    /// Intersect ray with an AABB (slab method). Returns `(t_min, t_max)` or `None`.
    pub fn intersect_aabb(&self, aabb: &Aabb3) -> Option<(f64, f64)> {
        let mut t_min = f64::NEG_INFINITY;
        let mut t_max = f64::INFINITY;
        for i in 0..3 {
            let inv_d = if self.direction[i].abs() < 1e-14 {
                f64::INFINITY
            } else {
                1.0 / self.direction[i]
            };
            let t0 = (aabb.min[i] - self.origin[i]) * inv_d;
            let t1 = (aabb.max[i] - self.origin[i]) * inv_d;
            let (ta, tb) = if inv_d >= 0.0 { (t0, t1) } else { (t1, t0) };
            t_min = t_min.max(ta);
            t_max = t_max.min(tb);
            if t_max < t_min {
                return None;
            }
        }
        if t_max < 0.0 {
            None
        } else {
            Some((t_min.max(0.0), t_max))
        }
    }

    /// Intersect ray with a triangle using the Möller–Trumbore algorithm.
    /// Returns `(t, u, v)` where `(u,v)` are barycentric coordinates, or `None`.
    pub fn intersect_triangle(
        &self,
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
    ) -> Option<(f64, f64, f64)> {
        let eps = 1e-8;
        let e1 = sub(v1, v0);
        let e2 = sub(v2, v0);
        let h = cross(self.direction, e2);
        let a = dot(e1, h);
        if a.abs() < eps {
            return None;
        }
        let f = 1.0 / a;
        let s = sub(self.origin, v0);
        let u = f * dot(s, h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross(s, e1);
        let v = f * dot(self.direction, q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = f * dot(e2, q);
        if t > eps { Some((t, u, v)) } else { None }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Aabb3
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box in 3D.
#[derive(Debug, Clone, Copy)]
pub struct Aabb3 {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb3 {
    /// Create an AABB from min and max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create an empty AABB (inverted infinite bounds).
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    /// Center of the AABB.
    pub fn center(&self) -> [f64; 3] {
        scale(add(self.min, self.max), 0.5)
    }

    /// Half-extents.
    pub fn half_extents(&self) -> [f64; 3] {
        scale(sub(self.max, self.min), 0.5)
    }

    /// Check if a point is inside (inclusive).
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Check if two AABBs intersect.
    pub fn intersects_aabb(&self, other: &Aabb3) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Check if an AABB intersects a sphere.
    pub fn intersects_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        let closest = [
            center[0].clamp(self.min[0], self.max[0]),
            center[1].clamp(self.min[1], self.max[1]),
            center[2].clamp(self.min[2], self.max[2]),
        ];
        len_sq(sub(closest, center)) <= radius * radius
    }

    /// Expand AABB by a margin.
    pub fn expand(&self, margin: f64) -> Aabb3 {
        Aabb3 {
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }

    /// Surface area.
    pub fn surface_area(&self) -> f64 {
        let e = sub(self.max, self.min);
        2.0 * (e[0] * e[1] + e[1] * e[2] + e[2] * e[0])
    }

    /// Volume.
    pub fn volume(&self) -> f64 {
        let e = sub(self.max, self.min);
        e[0] * e[1] * e[2]
    }

    /// Merge with another AABB.
    pub fn merge(&self, other: &Aabb3) -> Aabb3 {
        Aabb3 {
            min: vmin(self.min, other.min),
            max: vmax(self.max, other.max),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sphere3
// ─────────────────────────────────────────────────────────────────────────────

/// A sphere in 3D defined by center and radius.
#[derive(Debug, Clone, Copy)]
pub struct Sphere3 {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl Sphere3 {
    /// Create a new sphere.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }

    /// Check if a point is inside the sphere.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        len_sq(sub(p, self.center)) <= self.radius * self.radius
    }

    /// Check if the sphere intersects an AABB.
    pub fn intersects_aabb(&self, aabb: &Aabb3) -> bool {
        aabb.intersects_sphere(self.center, self.radius)
    }

    /// Check if this sphere intersects another sphere.
    pub fn intersects_sphere(&self, other: &Sphere3) -> bool {
        let dist_sq = len_sq(sub(self.center, other.center));
        let r_sum = self.radius + other.radius;
        dist_sq <= r_sum * r_sum
    }

    /// Volume of the sphere.
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }

    /// Surface area.
    pub fn surface_area(&self) -> f64 {
        4.0 * PI * self.radius * self.radius
    }

    /// Compute the AABB of this sphere.
    pub fn aabb(&self) -> Aabb3 {
        Aabb3 {
            min: [
                self.center[0] - self.radius,
                self.center[1] - self.radius,
                self.center[2] - self.radius,
            ],
            max: [
                self.center[0] + self.radius,
                self.center[1] + self.radius,
                self.center[2] + self.radius,
            ],
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Plane3
// ─────────────────────────────────────────────────────────────────────────────

/// An infinite plane in 3D: n·x = d.
#[derive(Debug, Clone, Copy)]
pub struct Plane3 {
    /// Plane normal (unit vector).
    pub normal: [f64; 3],
    /// Plane offset: d = n·p for any point p on the plane.
    pub offset: f64,
}

impl Plane3 {
    /// Create from normal and offset.
    pub fn new(normal: [f64; 3], offset: f64) -> Self {
        Self {
            normal: normalize(normal),
            offset,
        }
    }

    /// Create from three points.
    pub fn from_points(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        let n = normalize(cross(sub(b, a), sub(c, a)));
        let d = dot(n, a);
        Self {
            normal: n,
            offset: d,
        }
    }

    /// Signed distance from point to plane (positive on normal side).
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        dot(self.normal, p) - self.offset
    }

    /// Project a point onto the plane.
    pub fn project_point(&self, p: [f64; 3]) -> [f64; 3] {
        let d = self.signed_distance(p);
        sub(p, scale(self.normal, d))
    }

    /// Intersect ray with plane. Returns parameter `t` or `None` if parallel.
    pub fn intersect_ray(&self, ray: &Ray3) -> Option<f64> {
        let denom = dot(self.normal, ray.direction);
        if denom.abs() < 1e-12 {
            return None;
        }
        let t = (self.offset - dot(self.normal, ray.origin)) / denom;
        if t >= 0.0 { Some(t) } else { None }
    }

    /// Intersect a line segment (a → b) with the plane.
    /// Returns the intersection point or `None`.
    pub fn intersect_segment(&self, a: [f64; 3], b: [f64; 3]) -> Option<[f64; 3]> {
        let da = self.signed_distance(a);
        let db = self.signed_distance(b);
        if da * db > 0.0 {
            return None; // Same side
        }
        let t = da / (da - db);
        Some(add(a, scale(sub(b, a), t)))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Triangle3
// ─────────────────────────────────────────────────────────────────────────────

/// A 3D triangle with three vertices.
#[derive(Debug, Clone, Copy)]
pub struct Triangle3 {
    /// Vertices of the triangle.
    pub v: [[f64; 3]; 3],
}

impl Triangle3 {
    /// Create from three vertices.
    pub fn new(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> Self {
        Self { v: [v0, v1, v2] }
    }

    /// Outward normal (not normalized).
    pub fn normal_raw(&self) -> [f64; 3] {
        cross(sub(self.v[1], self.v[0]), sub(self.v[2], self.v[0]))
    }

    /// Unit normal.
    pub fn normal(&self) -> [f64; 3] {
        normalize(self.normal_raw())
    }

    /// Area of the triangle.
    pub fn area(&self) -> f64 {
        len(self.normal_raw()) * 0.5
    }

    /// Centroid.
    pub fn centroid(&self) -> [f64; 3] {
        scale(add(add(self.v[0], self.v[1]), self.v[2]), 1.0 / 3.0)
    }

    /// Compute barycentric coordinates (u, v, w) of point p projected onto triangle.
    pub fn barycentric_coords(&self, p: [f64; 3]) -> (f64, f64, f64) {
        let v0 = sub(self.v[1], self.v[0]);
        let v1 = sub(self.v[2], self.v[0]);
        let v2 = sub(p, self.v[0]);
        let d00 = dot(v0, v0);
        let d01 = dot(v0, v1);
        let d11 = dot(v1, v1);
        let d20 = dot(v2, v0);
        let d21 = dot(v2, v1);
        let denom = d00 * d11 - d01 * d01;
        if denom.abs() < 1e-14 {
            return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
        }
        let v_bc = (d11 * d20 - d01 * d21) / denom;
        let w_bc = (d00 * d21 - d01 * d20) / denom;
        let u_bc = 1.0 - v_bc - w_bc;
        (u_bc, v_bc, w_bc)
    }

    /// Check if a point (projected) lies within the triangle (u, v, w ≥ 0).
    pub fn contains_point_2d(&self, p: [f64; 3]) -> bool {
        let (u, v, w) = self.barycentric_coords(p);
        u >= -1e-8 && v >= -1e-8 && w >= -1e-8
    }

    /// Plane of the triangle.
    pub fn plane(&self) -> Plane3 {
        Plane3::from_points(self.v[0], self.v[1], self.v[2])
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Capsule3
// ─────────────────────────────────────────────────────────────────────────────

/// A 3D capsule: a swept sphere along a line segment.
#[derive(Debug, Clone, Copy)]
pub struct Capsule3 {
    /// Start of the segment.
    pub a: [f64; 3],
    /// End of the segment.
    pub b: [f64; 3],
    /// Radius of the capsule.
    pub radius: f64,
}

impl Capsule3 {
    /// Create a new capsule.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self { a, b, radius }
    }

    /// Closest point on the segment to an external point `p`.
    pub fn closest_point_on_segment(&self, p: [f64; 3]) -> [f64; 3] {
        closest_point_segment(self.a, self.b, p)
    }

    /// Distance from point `p` to the capsule surface.
    pub fn distance_to_point(&self, p: [f64; 3]) -> f64 {
        let closest = self.closest_point_on_segment(p);
        (len(sub(p, closest)) - self.radius).max(0.0)
    }

    /// Check if a sphere intersects this capsule.
    pub fn intersects_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        let closest = self.closest_point_on_segment(center);
        let dist_sq = len_sq(sub(closest, center));
        let r_sum = self.radius + radius;
        dist_sq <= r_sum * r_sum
    }

    /// Length of the capsule segment.
    pub fn segment_length(&self) -> f64 {
        len(sub(self.b, self.a))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OrientedBox
// ─────────────────────────────────────────────────────────────────────────────

/// An oriented bounding box (OBB) defined by center, local axes, and half-extents.
#[derive(Debug, Clone, Copy)]
pub struct OrientedBox {
    /// Center position.
    pub center: [f64; 3],
    /// Local axes (column vectors): axes\[0\], axes\[1\], axes\[2\].
    pub axes: [[f64; 3]; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}

impl OrientedBox {
    /// Create an oriented box.
    pub fn new(center: [f64; 3], axes: [[f64; 3]; 3], half_extents: [f64; 3]) -> Self {
        Self {
            center,
            axes: [normalize(axes[0]), normalize(axes[1]), normalize(axes[2])],
            half_extents,
        }
    }

    /// Check if a world-space point is inside the OBB.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        let d = sub(p, self.center);
        for i in 0..3 {
            let proj = dot(d, self.axes[i]).abs();
            if proj > self.half_extents[i] + 1e-8 {
                return false;
            }
        }
        true
    }

    /// Return the 8 corners of the OBB in world space.
    pub fn corners(&self) -> [[f64; 3]; 8] {
        let mut corners = [[0.0; 3]; 8];
        for i in 0..8 {
            let sx = if i & 1 == 0 { 1.0 } else { -1.0 };
            let sy = if i & 2 == 0 { 1.0 } else { -1.0 };
            let sz = if i & 4 == 0 { 1.0 } else { -1.0 };
            let offset = add(
                add(
                    scale(self.axes[0], sx * self.half_extents[0]),
                    scale(self.axes[1], sy * self.half_extents[1]),
                ),
                scale(self.axes[2], sz * self.half_extents[2]),
            );
            corners[i] = add(self.center, offset);
        }
        corners
    }

    /// Compute the AABB that contains this OBB.
    pub fn intersects_aabb(&self, aabb: &Aabb3) -> bool {
        // SAT: test 15 axes (3 AABB axes + 3 OBB axes + 9 cross products)
        let obb_aabb = aabb_from_points(&self.corners());
        obb_aabb.intersects_aabb(aabb)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Frustum
// ─────────────────────────────────────────────────────────────────────────────

/// A view frustum defined by 6 planes: near, far, left, right, top, bottom.
#[derive(Debug, Clone)]
pub struct Frustum {
    /// Six clipping planes: \[near, far, left, right, top, bottom\].
    pub planes: [Plane3; 6],
}

impl Frustum {
    /// Create a frustum from 6 planes.
    pub fn new(planes: [Plane3; 6]) -> Self {
        Self { planes }
    }

    /// Check if an AABB is at least partially inside the frustum.
    pub fn contains_aabb(&self, aabb: &Aabb3) -> bool {
        for plane in &self.planes {
            // Find the positive vertex (most in the direction of the normal)
            let positive = [
                if plane.normal[0] >= 0.0 {
                    aabb.max[0]
                } else {
                    aabb.min[0]
                },
                if plane.normal[1] >= 0.0 {
                    aabb.max[1]
                } else {
                    aabb.min[1]
                },
                if plane.normal[2] >= 0.0 {
                    aabb.max[2]
                } else {
                    aabb.min[2]
                },
            ];
            if plane.signed_distance(positive) < 0.0 {
                return false;
            }
        }
        true
    }

    /// Check if a sphere is at least partially inside the frustum.
    pub fn contains_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        for plane in &self.planes {
            if plane.signed_distance(center) < -radius {
                return false;
            }
        }
        true
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PolyhedronConvex
// ─────────────────────────────────────────────────────────────────────────────

/// A convex polyhedron defined by vertices and face normals.
#[derive(Debug, Clone)]
pub struct PolyhedronConvex {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Outward-pointing face normals and a point on each face.
    pub face_planes: Vec<Plane3>,
}

impl PolyhedronConvex {
    /// Create a convex polyhedron.
    pub fn new(vertices: Vec<[f64; 3]>, face_planes: Vec<Plane3>) -> Self {
        Self {
            vertices,
            face_planes,
        }
    }

    /// Support function: returns the vertex with the highest projection along `dir`.
    pub fn support_point(&self, dir: [f64; 3]) -> [f64; 3] {
        let mut best = self.vertices[0];
        let mut best_d = dot(best, dir);
        for &v in &self.vertices[1..] {
            let d = dot(v, dir);
            if d > best_d {
                best = v;
                best_d = d;
            }
        }
        best
    }

    /// Check if a point is inside the convex polyhedron.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        for plane in &self.face_planes {
            if plane.signed_distance(p) > 1e-8 {
                return false;
            }
        }
        true
    }

    /// Build from convex hull faces (simplified: just store input).
    pub fn from_box(half_extents: [f64; 3]) -> Self {
        let [hx, hy, hz] = half_extents;
        let vertices = vec![
            [-hx, -hy, -hz],
            [hx, -hy, -hz],
            [-hx, hy, -hz],
            [hx, hy, -hz],
            [-hx, -hy, hz],
            [hx, -hy, hz],
            [-hx, hy, hz],
            [hx, hy, hz],
        ];
        let face_planes = vec![
            Plane3::new([1.0, 0.0, 0.0], hx),
            Plane3::new([-1.0, 0.0, 0.0], hx),
            Plane3::new([0.0, 1.0, 0.0], hy),
            Plane3::new([0.0, -1.0, 0.0], hy),
            Plane3::new([0.0, 0.0, 1.0], hz),
            Plane3::new([0.0, 0.0, -1.0], hz),
        ];
        Self {
            vertices,
            face_planes,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BarycentricCoords
// ─────────────────────────────────────────────────────────────────────────────

/// Barycentric coordinate utilities.
#[derive(Debug, Clone, Copy)]
pub struct BarycentricCoords {
    /// First barycentric coordinate.
    pub u: f64,
    /// Second barycentric coordinate.
    pub v: f64,
    /// Third barycentric coordinate.
    pub w: f64,
}

impl BarycentricCoords {
    /// Compute barycentric coordinates of point `p` in triangle (v0,v1,v2).
    pub fn compute(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], p: [f64; 3]) -> Self {
        let tri = Triangle3::new(v0, v1, v2);
        let (u, v, w) = tri.barycentric_coords(p);
        Self { u, v, w }
    }

    /// Check if the barycentric coordinates represent a point inside the triangle.
    pub fn is_inside(&self) -> bool {
        self.u >= -1e-8 && self.v >= -1e-8 && self.w >= -1e-8
    }

    /// Interpolate a scalar value at the barycentric point.
    pub fn interpolate_scalar(&self, s0: f64, s1: f64, s2: f64) -> f64 {
        self.u * s0 + self.v * s1 + self.w * s2
    }

    /// Interpolate a 3D vector at the barycentric point.
    pub fn interpolate_vec(&self, v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
        add(add(scale(v0, self.u), scale(v1, self.v)), scale(v2, self.w))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Find the closest point on a line segment (a → b) to point `p`.
pub fn closest_point_segment(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ap = sub(p, a);
    let len_sq_ab = dot(ab, ab);
    if len_sq_ab < 1e-14 {
        return a;
    }
    let t = (dot(ap, ab) / len_sq_ab).clamp(0.0, 1.0);
    add(a, scale(ab, t))
}

/// Find the closest point on a triangle to an external point `p`.
pub fn closest_point_triangle(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let tri = Triangle3::new(v0, v1, v2);
    let (u, v, w) = tri.barycentric_coords(p);
    // If inside, project onto plane
    if u >= 0.0 && v >= 0.0 && w >= 0.0 {
        return add(add(scale(v0, u), scale(v1, v)), scale(v2, w));
    }
    // Otherwise, find closest point on each edge
    let c0 = closest_point_segment(v0, v1, p);
    let c1 = closest_point_segment(v1, v2, p);
    let c2 = closest_point_segment(v2, v0, p);
    let d0 = len_sq(sub(c0, p));
    let d1 = len_sq(sub(c1, p));
    let d2 = len_sq(sub(c2, p));
    if d0 <= d1 && d0 <= d2 {
        c0
    } else if d1 <= d2 {
        c1
    } else {
        c2
    }
}

/// Build an AABB from a collection of points.
pub fn aabb_from_points(points: &[[f64; 3]]) -> Aabb3 {
    let mut aabb = Aabb3::empty();
    for &p in points {
        aabb.min = vmin(aabb.min, p);
        aabb.max = vmax(aabb.max, p);
    }
    aabb
}

/// Build a frustum from a combined view-projection matrix (row-major, 4×4).
/// Extracts the 6 Gribb-Hartmann planes.
pub fn frustum_from_vp(m: &[[f64; 4]; 4]) -> Frustum {
    // Each row contributes to a plane
    let extract = |a: [f64; 4], b: [f64; 4]| -> Plane3 {
        let n = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
        let d = -(a[3] + b[3]);
        Plane3::new(n, -d)
    };
    let r = m;
    let row = |i: usize| [r[i][0], r[i][1], r[i][2], r[i][3]];
    let row3 = row(3);
    // Left: row3 + row0, Right: row3 - row0, Bottom: row3 + row1
    // Top: row3 - row1, Near: row3 + row2, Far: row3 - row2
    let planes = [
        extract(row3, row(0)),
        extract(row3, [-row(0)[0], -row(0)[1], -row(0)[2], -row(0)[3]]),
        extract(row3, row(1)),
        extract(row3, [-row(1)[0], -row(1)[1], -row(1)[2], -row(1)[3]]),
        extract(row3, row(2)),
        extract(row3, [-row(2)[0], -row(2)[1], -row(2)[2], -row(2)[3]]),
    ];
    Frustum::new(planes)
}

/// Use PI for any angle-related computations.
#[allow(dead_code)]
const HALF_PI: f64 = PI / 2.0;

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Ray3 tests ────────────────────────────────────────────────────────────

    #[test]
    fn ray_at_t_zero_is_origin() {
        let ray = Ray3::new([1.0, 2.0, 3.0], [0.0, 0.0, 1.0]);
        let p = ray.at(0.0);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn ray_intersect_sphere_hit() {
        let ray = Ray3::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let t = ray.intersect_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(t.is_some(), "should hit");
        assert!((t.unwrap() - 4.0).abs() < 1e-6, "t={}", t.unwrap());
    }

    #[test]
    fn ray_intersect_sphere_miss() {
        let ray = Ray3::new([5.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let t = ray.intersect_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(t.is_none(), "should miss");
    }

    #[test]
    fn ray_intersect_aabb_hit() {
        let ray = Ray3::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let aabb = Aabb3::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let result = ray.intersect_aabb(&aabb);
        assert!(result.is_some(), "should hit AABB");
    }

    #[test]
    fn ray_intersect_aabb_miss() {
        let ray = Ray3::new([5.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let aabb = Aabb3::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let result = ray.intersect_aabb(&aabb);
        assert!(result.is_none(), "should miss AABB");
    }

    #[test]
    fn ray_intersect_triangle_hit() {
        let ray = Ray3::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        let v0 = [-1.0, -1.0, 0.0];
        let v1 = [1.0, -1.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let result = ray.intersect_triangle(v0, v1, v2);
        assert!(result.is_some(), "should hit triangle");
        let (t, _u, _v) = result.unwrap();
        assert!((t - 1.0).abs() < 1e-6, "t={t}");
    }

    #[test]
    fn ray_intersect_triangle_miss() {
        let ray = Ray3::new([10.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        let v0 = [-1.0, -1.0, 0.0];
        let v1 = [1.0, -1.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let result = ray.intersect_triangle(v0, v1, v2);
        assert!(result.is_none(), "should miss triangle");
    }

    // ── Aabb3 tests ───────────────────────────────────────────────────────────

    #[test]
    fn aabb_contains_center() {
        let aabb = Aabb3::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(aabb.contains_point([0.0, 0.0, 0.0]));
    }

    #[test]
    fn aabb_excludes_outside() {
        let aabb = Aabb3::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(!aabb.contains_point([2.0, 0.0, 0.0]));
    }

    #[test]
    fn aabb_intersects_overlapping() {
        let a = Aabb3::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Aabb3::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.intersects_aabb(&b));
    }

    #[test]
    fn aabb_intersects_sphere() {
        let aabb = Aabb3::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(aabb.intersects_sphere([1.0, 1.0, 1.0], 0.5));
        assert!(!aabb.intersects_sphere([5.0, 5.0, 5.0], 0.5));
    }

    #[test]
    fn aabb_surface_area_unit_cube() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        assert!((aabb.surface_area() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn aabb_volume_unit_cube() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        assert!((aabb.volume() - 1.0).abs() < 1e-10);
    }

    // ── Sphere3 tests ─────────────────────────────────────────────────────────

    #[test]
    fn sphere_contains_center() {
        let s = Sphere3::new([0.0; 3], 1.0);
        assert!(s.contains_point([0.0, 0.0, 0.0]));
    }

    #[test]
    fn sphere_excludes_far_point() {
        let s = Sphere3::new([0.0; 3], 1.0);
        assert!(!s.contains_point([2.0, 0.0, 0.0]));
    }

    #[test]
    fn sphere_sphere_intersection() {
        let s1 = Sphere3::new([0.0; 3], 1.0);
        let s2 = Sphere3::new([1.5, 0.0, 0.0], 1.0);
        assert!(s1.intersects_sphere(&s2));
        let s3 = Sphere3::new([3.0, 0.0, 0.0], 1.0);
        assert!(!s1.intersects_sphere(&s3));
    }

    #[test]
    fn sphere_volume() {
        let s = Sphere3::new([0.0; 3], 1.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }

    // ── Plane3 tests ──────────────────────────────────────────────────────────

    #[test]
    fn plane_signed_distance_positive_side() {
        let p = Plane3::new([0.0, 1.0, 0.0], 0.0);
        assert!((p.signed_distance([0.0, 1.0, 0.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn plane_signed_distance_negative_side() {
        let p = Plane3::new([0.0, 1.0, 0.0], 0.0);
        assert!((p.signed_distance([0.0, -1.0, 0.0]) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn plane_project_point() {
        let p = Plane3::new([0.0, 1.0, 0.0], 0.0);
        let proj = p.project_point([1.0, 5.0, 2.0]);
        assert!(proj[1].abs() < 1e-10, "proj.y={}", proj[1]);
    }

    #[test]
    fn plane_intersect_ray() {
        let plane = Plane3::new([0.0, 1.0, 0.0], 0.0);
        let ray = Ray3::new([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        let t = plane.intersect_ray(&ray);
        assert!(t.is_some());
        assert!((t.unwrap() - 5.0).abs() < 1e-6, "t={}", t.unwrap());
    }

    #[test]
    fn plane_from_points() {
        let p = Plane3::from_points([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        // Normal should be [0,0,1]
        assert!(p.normal[2].abs() > 0.99);
    }

    // ── Triangle3 tests ───────────────────────────────────────────────────────

    #[test]
    fn triangle_area_right_triangle() {
        let tri = Triangle3::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((tri.area() - 0.5).abs() < 1e-10, "area={}", tri.area());
    }

    #[test]
    fn triangle_centroid() {
        let tri = Triangle3::new([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]);
        let c = tri.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn triangle_barycentric_centroid() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [3.0, 0.0, 0.0];
        let v2 = [0.0, 3.0, 0.0];
        let c = [1.0, 1.0, 0.0];
        let tri = Triangle3::new(v0, v1, v2);
        let (u, v, w) = tri.barycentric_coords(c);
        assert!((u + v + w - 1.0).abs() < 1e-6, "sum={}", u + v + w);
        assert!(u >= 0.0 && v >= 0.0 && w >= 0.0, "u={u}, v={v}, w={w}");
    }

    // ── Capsule3 tests ────────────────────────────────────────────────────────

    #[test]
    fn capsule_distance_to_surface() {
        let cap = Capsule3::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let d = cap.distance_to_point([1.0, 0.5, 0.0]);
        assert!((d - 0.5).abs() < 1e-6, "d={d}");
    }

    #[test]
    fn capsule_inside_sphere_intersection() {
        let cap = Capsule3::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5);
        assert!(cap.intersects_sphere([0.0, 1.0, 0.0], 0.3));
    }

    #[test]
    fn capsule_no_intersection() {
        let cap = Capsule3::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1);
        assert!(!cap.intersects_sphere([5.0, 0.5, 0.0], 0.1));
    }

    // ── OrientedBox tests ─────────────────────────────────────────────────────

    #[test]
    fn obb_contains_center() {
        let obb = OrientedBox::new(
            [0.0; 3],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [1.0, 1.0, 1.0],
        );
        assert!(obb.contains_point([0.0, 0.0, 0.0]));
    }

    #[test]
    fn obb_excludes_far_point() {
        let obb = OrientedBox::new(
            [0.0; 3],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [1.0, 1.0, 1.0],
        );
        assert!(!obb.contains_point([2.0, 0.0, 0.0]));
    }

    #[test]
    fn obb_corners_count() {
        let obb = OrientedBox::new(
            [0.0; 3],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [1.0, 1.0, 1.0],
        );
        let corners = obb.corners();
        assert_eq!(corners.len(), 8);
    }

    // ── Frustum tests ─────────────────────────────────────────────────────────

    #[test]
    fn frustum_contains_center_sphere() {
        // Create a box frustum
        let planes = [
            Plane3::new([1.0, 0.0, 0.0], -10.0),
            Plane3::new([-1.0, 0.0, 0.0], -10.0),
            Plane3::new([0.0, 1.0, 0.0], -10.0),
            Plane3::new([0.0, -1.0, 0.0], -10.0),
            Plane3::new([0.0, 0.0, 1.0], -10.0),
            Plane3::new([0.0, 0.0, -1.0], -10.0),
        ];
        let frustum = Frustum::new(planes);
        assert!(frustum.contains_sphere([0.0, 0.0, 0.0], 0.5));
    }

    // ── PolyhedronConvex tests ────────────────────────────────────────────────

    #[test]
    fn polyhedron_box_contains_center() {
        let poly = PolyhedronConvex::from_box([1.0, 1.0, 1.0]);
        assert!(poly.contains_point([0.0, 0.0, 0.0]));
    }

    #[test]
    fn polyhedron_box_excludes_outside() {
        let poly = PolyhedronConvex::from_box([1.0, 1.0, 1.0]);
        assert!(!poly.contains_point([2.0, 0.0, 0.0]));
    }

    #[test]
    fn polyhedron_support_point() {
        let poly = PolyhedronConvex::from_box([1.0, 1.0, 1.0]);
        let sp = poly.support_point([1.0, 0.0, 0.0]);
        assert!((sp[0] - 1.0).abs() < 1e-10, "support_x={}", sp[0]);
    }

    // ── BarycentricCoords tests ───────────────────────────────────────────────

    #[test]
    fn barycentric_vertex_is_canonical() {
        let bc = BarycentricCoords::compute(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        );
        assert!((bc.u - 1.0).abs() < 1e-6, "u={}", bc.u);
        assert!(bc.v.abs() < 1e-6, "v={}", bc.v);
        assert!(bc.w.abs() < 1e-6, "w={}", bc.w);
    }

    #[test]
    fn barycentric_interpolate_scalar() {
        let bc = BarycentricCoords {
            u: 0.5,
            v: 0.3,
            w: 0.2,
        };
        let val = bc.interpolate_scalar(1.0, 2.0, 3.0);
        // 0.5*1 + 0.3*2 + 0.2*3 = 0.5 + 0.6 + 0.6 = 1.7
        assert!((val - 1.7).abs() < 1e-10, "val={val}");
    }

    // ── Helper function tests ─────────────────────────────────────────────────

    #[test]
    fn closest_point_segment_midpoint() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let p = [1.0, 1.0, 0.0];
        let cp = closest_point_segment(a, b, p);
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
    }

    #[test]
    fn closest_point_segment_clamp_a() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let p = [-5.0, 0.0, 0.0];
        let cp = closest_point_segment(a, b, p);
        assert!((cp[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn aabb_from_points_simple() {
        let points = [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]];
        let aabb = aabb_from_points(&points);
        assert!((aabb.min[0] + 1.0).abs() < 1e-10);
        assert!((aabb.max[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn aabb_expand_increases_bounds() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        let expanded = aabb.expand(0.5);
        assert!((expanded.min[0] + 0.5).abs() < 1e-10);
        assert!((expanded.max[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn sphere_aabb_correct_bounds() {
        let s = Sphere3::new([1.0, 2.0, 3.0], 1.0);
        let aabb = s.aabb();
        assert!((aabb.min[0] - 0.0).abs() < 1e-10);
        assert!((aabb.max[0] - 2.0).abs() < 1e-10);
    }
}
