// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compound and complex shape handling for collision detection.
//!
//! Provides:
//! - `CompoundShape` — collection of child shapes with local transforms
//! - `ConvexDecomposition` — iterative approximate convex decomposition
//! - `TriangleMeshShape` — triangle mesh with BVH and normals
//! - `MeshShapeQuery` — point-in-mesh, closest point, normal
//! - `GjkPolytope` — convex polytope with support function
//! - `MinkowskiSum` — Minkowski sum for GJK distance
//! - `ShapeProxy` — dispatch wrapper for shape types
//! - `CompoundBvh` — BVH over compound sub-shapes
//! - `ShapeOffset` — Minkowski sum with sphere (rounded shapes)
//! - `ImportedCollisionMesh` — mesh from vertex data

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Squared norm of a 3-vector.
#[inline]
pub fn norm_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

/// Norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    norm_sq3(a).sqrt()
}

/// Normalize a 3-vector. Returns zero vector if degenerate.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-15 {
        return [0.0; 3];
    }
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ---------------------------------------------------------------------------
// AABB (Axis-Aligned Bounding Box)
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb3 {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb3 {
    /// Create a new AABB.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create a degenerate (empty) AABB.
    pub fn empty() -> Self {
        Self {
            min: [f64::MAX; 3],
            max: [f64::MIN; 3],
        }
    }

    /// Expand AABB to include point p.
    pub fn include_point(&mut self, p: [f64; 3]) {
        self.min
            .iter_mut()
            .zip(p.iter())
            .for_each(|(m, &pi)| *m = m.min(pi));
        self.max
            .iter_mut()
            .zip(p.iter())
            .for_each(|(m, &pi)| *m = m.max(pi));
    }

    /// Expand AABB to include another AABB.
    pub fn include_aabb(&mut self, other: &Aabb3) {
        self.min
            .iter_mut()
            .zip(other.min.iter())
            .for_each(|(m, &o)| *m = m.min(o));
        self.max
            .iter_mut()
            .zip(other.max.iter())
            .for_each(|(m, &o)| *m = m.max(o));
    }

    /// Center of this AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Half-extents of this AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            0.5 * (self.max[0] - self.min[0]),
            0.5 * (self.max[1] - self.min[1]),
            0.5 * (self.max[2] - self.min[2]),
        ]
    }

    /// Volume of this AABB.
    pub fn volume(&self) -> f64 {
        let e = [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ];
        e[0].max(0.0) * e[1].max(0.0) * e[2].max(0.0)
    }

    /// Check if two AABBs overlap.
    pub fn overlaps(&self, other: &Aabb3) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Check if point is inside AABB.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Inflate AABB by margin on all sides.
    pub fn inflate(&self, margin: f64) -> Aabb3 {
        Aabb3::new(
            [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        )
    }
}

// ---------------------------------------------------------------------------
// Local Transform (position + orientation as 3x3 rotation)
// ---------------------------------------------------------------------------

/// Simple local transform: position + 3x3 rotation matrix.
#[derive(Debug, Clone, Copy)]
pub struct LocalTransform {
    /// Translation.
    pub position: [f64; 3],
    /// Rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
}

impl LocalTransform {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            position: [0.0; 3],
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Translation-only transform.
    pub fn from_translation(pos: [f64; 3]) -> Self {
        let mut t = Self::identity();
        t.position = pos;
        t
    }

    /// Apply transform to point p: p' = R * p + t.
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        let rp = [
            r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
            r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
            r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
        ];
        add3(rp, self.position)
    }

    /// Apply inverse transform: p_local = R^T * (p - t).
    pub fn inverse_apply(&self, p: [f64; 3]) -> [f64; 3] {
        let dp = sub3(p, self.position);
        let r = &self.rotation;
        [
            r[0][0] * dp[0] + r[1][0] * dp[1] + r[2][0] * dp[2],
            r[0][1] * dp[0] + r[1][1] * dp[1] + r[2][1] * dp[2],
            r[0][2] * dp[0] + r[1][2] * dp[1] + r[2][2] * dp[2],
        ]
    }
}

// ---------------------------------------------------------------------------
// Primitive shape types (for compound)
// ---------------------------------------------------------------------------

/// Primitive shape variant.
#[derive(Debug, Clone)]
pub enum PrimitiveShape {
    /// Sphere with radius.
    Sphere(f64),
    /// Box with half-extents \[hx, hy, hz\].
    Box([f64; 3]),
    /// Capsule with radius and half-height.
    Capsule(f64, f64),
    /// Cylinder with radius and half-height.
    Cylinder(f64, f64),
    /// Convex hull with vertices.
    ConvexHull(Vec<[f64; 3]>),
}

impl PrimitiveShape {
    /// Compute AABB of this primitive shape.
    pub fn aabb(&self) -> Aabb3 {
        match self {
            PrimitiveShape::Sphere(r) => Aabb3::new([-r, -r, -r], [*r, *r, *r]),
            PrimitiveShape::Box(he) => Aabb3::new([-he[0], -he[1], -he[2]], [he[0], he[1], he[2]]),
            PrimitiveShape::Capsule(r, hh) => Aabb3::new([-r, -hh - r, -r], [*r, hh + r, *r]),
            PrimitiveShape::Cylinder(r, hh) => Aabb3::new([-r, -hh, -r], [*r, *hh, *r]),
            PrimitiveShape::ConvexHull(verts) => {
                let mut aabb = Aabb3::empty();
                for &v in verts {
                    aabb.include_point(v);
                }
                aabb
            }
        }
    }

    /// Compute support point in direction d.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        match self {
            PrimitiveShape::Sphere(r) => {
                let n = norm3(d);
                if n < 1e-15 {
                    return [0.0; 3];
                }
                scale3(d, r / n)
            }
            PrimitiveShape::Box(he) => [
                he[0] * d[0].signum(),
                he[1] * d[1].signum(),
                he[2] * d[2].signum(),
            ],
            PrimitiveShape::Capsule(r, hh) => {
                let n = norm3(d);
                let dn = if n < 1e-15 {
                    [0.0, 1.0, 0.0]
                } else {
                    scale3(d, 1.0 / n)
                };
                let cap_center = [0.0, dn[1].signum() * hh, 0.0];
                add3(cap_center, scale3(dn, *r))
            }
            PrimitiveShape::Cylinder(r, hh) => {
                let xy_len = (d[0] * d[0] + d[2] * d[2]).sqrt();
                let x = if xy_len > 1e-15 {
                    r * d[0] / xy_len
                } else {
                    0.0
                };
                let z = if xy_len > 1e-15 {
                    r * d[2] / xy_len
                } else {
                    0.0
                };
                [x, hh * d[1].signum(), z]
            }
            PrimitiveShape::ConvexHull(verts) => {
                let mut best = [0.0f64; 3];
                let mut best_d = f64::NEG_INFINITY;
                for &v in verts {
                    let proj = dot3(v, d);
                    if proj > best_d {
                        best_d = proj;
                        best = v;
                    }
                }
                best
            }
        }
    }

    /// Volume of this primitive.
    ///
    /// For [`PrimitiveShape::ConvexHull`] this is the volume of the convex hull
    /// of the vertices; it is zero when the vertices do not span three
    /// dimensions (fewer than four points, or all coplanar).
    pub fn volume(&self) -> f64 {
        let pi = std::f64::consts::PI;
        match self {
            PrimitiveShape::Sphere(r) => 4.0 / 3.0 * pi * r * r * r,
            PrimitiveShape::Box(he) => 8.0 * he[0] * he[1] * he[2],
            PrimitiveShape::Capsule(r, hh) => pi * r * r * 2.0 * hh + 4.0 / 3.0 * pi * r * r * r,
            PrimitiveShape::Cylinder(r, hh) => pi * r * r * 2.0 * hh,
            PrimitiveShape::ConvexHull(verts) => hull_volume(verts),
        }
    }
}

// ---------------------------------------------------------------------------
// CompoundShape
// ---------------------------------------------------------------------------

/// Compound shape: a collection of primitive shapes with local transforms.
///
/// Provides AABB, volume, and center-of-mass computation over all children.
#[derive(Debug, Clone)]
pub struct CompoundShape {
    /// Child shapes with their local transforms.
    pub children: Vec<(PrimitiveShape, LocalTransform)>,
}

impl CompoundShape {
    /// Create an empty compound shape.
    pub fn new() -> Self {
        Self {
            children: Vec::new(),
        }
    }

    /// Add a child shape with a local transform.
    pub fn add_child(&mut self, shape: PrimitiveShape, transform: LocalTransform) {
        self.children.push((shape, transform));
    }

    /// Compute the world AABB of the compound shape.
    pub fn aabb(&self) -> Aabb3 {
        let mut aabb = Aabb3::empty();
        for (shape, xform) in &self.children {
            let local_aabb = shape.aabb();
            // Transform all 8 corners
            let corners = [
                [local_aabb.min[0], local_aabb.min[1], local_aabb.min[2]],
                [local_aabb.max[0], local_aabb.min[1], local_aabb.min[2]],
                [local_aabb.min[0], local_aabb.max[1], local_aabb.min[2]],
                [local_aabb.max[0], local_aabb.max[1], local_aabb.min[2]],
                [local_aabb.min[0], local_aabb.min[1], local_aabb.max[2]],
                [local_aabb.max[0], local_aabb.min[1], local_aabb.max[2]],
                [local_aabb.min[0], local_aabb.max[1], local_aabb.max[2]],
                [local_aabb.max[0], local_aabb.max[1], local_aabb.max[2]],
            ];
            for c in &corners {
                aabb.include_point(xform.apply(*c));
            }
        }
        aabb
    }

    /// Compute total volume (sum of child volumes).
    pub fn total_volume(&self) -> f64 {
        self.children.iter().map(|(s, _)| s.volume()).sum()
    }

    /// Compute center of mass (mass-weighted centroid, assuming uniform density).
    pub fn center_of_mass(&self) -> [f64; 3] {
        let mut total_vol = 0.0f64;
        let mut com = [0.0f64; 3];
        for (shape, xform) in &self.children {
            let v = shape.volume();
            let c = xform.position; // center of child at transform origin
            com[0] += v * c[0];
            com[1] += v * c[1];
            com[2] += v * c[2];
            total_vol += v;
        }
        if total_vol < 1e-20 {
            return [0.0; 3];
        }
        scale3(com, 1.0 / total_vol)
    }

    /// Support point for the compound shape in direction d.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        let mut best = [0.0f64; 3];
        let mut best_proj = f64::NEG_INFINITY;
        for (shape, xform) in &self.children {
            // Transform direction to local frame
            let d_local = xform.inverse_apply(add3(d, xform.position));
            let local_sup = shape.support(d_local);
            let world_sup = xform.apply(local_sup);
            let proj = dot3(world_sup, d);
            if proj > best_proj {
                best_proj = proj;
                best = world_sup;
            }
        }
        best
    }

    /// Number of child shapes.
    pub fn num_children(&self) -> usize {
        self.children.len()
    }
}

impl Default for CompoundShape {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute the AABB of a compound shape (free function).
pub fn compute_compound_aabb(compound: &CompoundShape) -> Aabb3 {
    compound.aabb()
}

// ---------------------------------------------------------------------------
// ConvexDecomposition
// ---------------------------------------------------------------------------

mod convex_decomposition;
pub use convex_decomposition::ConvexDecomposition;
use convex_decomposition::{hull_vertices, hull_volume};
#[cfg(test)]
mod test_support;

/// Merge multiple convex hulls into one (free function).
pub fn merge_hulls(hulls: &[Vec<[f64; 3]>]) -> Vec<[f64; 3]> {
    let mut all_verts = Vec::new();
    for h in hulls {
        all_verts.extend_from_slice(h);
    }
    all_verts
}

// ---------------------------------------------------------------------------
// TriangleMeshShape
// ---------------------------------------------------------------------------

/// Triangle mesh as a collision shape.
///
/// Stores triangle data with precomputed normals and a BVH for fast queries.
#[derive(Debug, Clone)]
pub struct TriangleMeshShape {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices (triplets).
    pub indices: Vec<[usize; 3]>,
    /// Precomputed per-face normals.
    pub normals: Vec<[f64; 3]>,
    /// Per-face AABB.
    pub face_aabbs: Vec<Aabb3>,
    /// Overall AABB.
    pub aabb: Aabb3,
}

impl TriangleMeshShape {
    /// Create a triangle mesh shape from vertices and indices.
    pub fn new(vertices: Vec<[f64; 3]>, indices: Vec<[usize; 3]>) -> Self {
        let normals = Self::compute_normals(&vertices, &indices);
        let face_aabbs = Self::compute_face_aabbs(&vertices, &indices);
        let mut aabb = Aabb3::empty();
        for &v in &vertices {
            aabb.include_point(v);
        }
        Self {
            vertices,
            indices,
            normals,
            face_aabbs,
            aabb,
        }
    }

    /// Compute per-face normals.
    fn compute_normals(verts: &[[f64; 3]], indices: &[[usize; 3]]) -> Vec<[f64; 3]> {
        indices
            .iter()
            .map(|tri| {
                let a = verts[tri[0]];
                let b = verts[tri[1]];
                let c = verts[tri[2]];
                let ab = sub3(b, a);
                let ac = sub3(c, a);
                normalize3(cross3(ab, ac))
            })
            .collect()
    }

    /// Compute per-face AABBs.
    fn compute_face_aabbs(verts: &[[f64; 3]], indices: &[[usize; 3]]) -> Vec<Aabb3> {
        indices
            .iter()
            .map(|tri| {
                let mut aabb = Aabb3::empty();
                for &vi in tri {
                    aabb.include_point(verts[vi]);
                }
                aabb
            })
            .collect()
    }

    /// Number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.indices.len()
    }

    /// Get a triangle by index.
    pub fn triangle(&self, i: usize) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let tri = &self.indices[i];
        (
            self.vertices[tri[0]],
            self.vertices[tri[1]],
            self.vertices[tri[2]],
        )
    }

    /// Compute total surface area.
    pub fn surface_area(&self) -> f64 {
        (0..self.num_triangles())
            .map(|i| {
                let (a, b, c) = self.triangle(i);
                let ab = sub3(b, a);
                let ac = sub3(c, a);
                0.5 * norm3(cross3(ab, ac))
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// MeshShapeQuery
// ---------------------------------------------------------------------------

/// Query operations on triangle mesh shapes.
///
/// Provides point-in-mesh test, closest point, and normal at hit.
#[derive(Debug, Clone)]
pub struct MeshShapeQuery<'a> {
    /// Reference to the mesh.
    pub mesh: &'a TriangleMeshShape,
}

impl<'a> MeshShapeQuery<'a> {
    /// Create a new mesh query.
    pub fn new(mesh: &'a TriangleMeshShape) -> Self {
        Self { mesh }
    }

    /// Point-in-mesh test using ray-parity method.
    ///
    /// Casts a ray in +x direction and counts crossings.
    pub fn point_inside(&self, point: [f64; 3]) -> bool {
        let dir = [1.0, 0.0, 0.0];
        let count = (0..self.mesh.num_triangles())
            .filter(|&i| {
                let (a, b, c) = self.mesh.triangle(i);
                Self::ray_triangle_intersect(point, dir, a, b, c).is_some()
            })
            .count();
        count % 2 == 1
    }

    /// Möller-Trumbore ray-triangle intersection.
    ///
    /// Returns t parameter if hit, None otherwise.
    pub fn ray_triangle_intersect(
        origin: [f64; 3],
        direction: [f64; 3],
        a: [f64; 3],
        b: [f64; 3],
        c: [f64; 3],
    ) -> Option<f64> {
        let edge1 = sub3(b, a);
        let edge2 = sub3(c, a);
        let h = cross3(direction, edge2);
        let det = dot3(edge1, h);
        if det.abs() < 1e-14 {
            return None;
        }
        let inv_det = 1.0 / det;
        let s = sub3(origin, a);
        let u = inv_det * dot3(s, h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross3(s, edge1);
        let v = inv_det * dot3(direction, q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = inv_det * dot3(edge2, q);
        if t > 1e-10 { Some(t) } else { None }
    }

    /// Find closest point on mesh surface to query point.
    ///
    /// Returns (closest_point, face_index, distance).
    pub fn closest_point(&self, query: [f64; 3]) -> ([f64; 3], usize, f64) {
        let mut best_pt = [0.0f64; 3];
        let mut best_face = 0;
        let mut best_dist = f64::MAX;

        for (i, _) in (0..self.mesh.num_triangles()).enumerate() {
            let (a, b, c) = self.mesh.triangle(i);
            let pt = Self::closest_point_triangle(query, a, b, c);
            let d = norm3(sub3(query, pt));
            if d < best_dist {
                best_dist = d;
                best_pt = pt;
                best_face = i;
            }
        }
        (best_pt, best_face, best_dist)
    }

    /// Closest point on triangle (a, b, c) to point p.
    pub fn closest_point_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let ap = sub3(p, a);

        let d1 = dot3(ab, ap);
        let d2 = dot3(ac, ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            return a;
        }

        let bp = sub3(p, b);
        let d3 = dot3(ab, bp);
        let d4 = dot3(ac, bp);
        if d3 >= 0.0 && d4 <= d3 {
            return b;
        }

        let cp = sub3(p, c);
        let d5 = dot3(ab, cp);
        let d6 = dot3(ac, cp);
        if d6 >= 0.0 && d5 <= d6 {
            return c;
        }

        // Edge AB
        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return add3(a, scale3(ab, v));
        }

        // Edge AC
        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            return add3(a, scale3(ac, w));
        }

        // Edge BC
        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return add3(b, scale3(sub3(c, b), w));
        }

        // Inside triangle
        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        add3(a, add3(scale3(ab, v), scale3(ac, w)))
    }

    /// Get surface normal at query point (closest face normal).
    pub fn normal_at(&self, query: [f64; 3]) -> [f64; 3] {
        let (_, face_idx, _) = self.closest_point(query);
        *self.mesh.normals.get(face_idx).unwrap_or(&[0.0, 1.0, 0.0])
    }
}

// ---------------------------------------------------------------------------
// GjkPolytope
// ---------------------------------------------------------------------------

/// Generic convex polytope for GJK distance computation.
///
/// Defined by a finite set of vertices and a support function.
#[derive(Debug, Clone)]
pub struct GjkPolytope {
    /// Vertex set.
    pub vertices: Vec<[f64; 3]>,
    /// Center (for numerical stability).
    pub center: [f64; 3],
}

impl GjkPolytope {
    /// Create a new GJK polytope from vertices.
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        let n = vertices.len();
        let center = if n > 0 {
            let s: [f64; 3] = vertices.iter().fold([0.0; 3], |a, &v| add3(a, v));
            scale3(s, 1.0 / n as f64)
        } else {
            [0.0; 3]
        };
        Self { vertices, center }
    }

    /// Compute support point in direction d.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        let mut best = self.center;
        let mut best_proj = f64::NEG_INFINITY;
        for &v in &self.vertices {
            let proj = dot3(v, d);
            if proj > best_proj {
                best_proj = proj;
                best = v;
            }
        }
        best
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// AABB of this polytope.
    pub fn aabb(&self) -> Aabb3 {
        let mut aabb = Aabb3::empty();
        for &v in &self.vertices {
            aabb.include_point(v);
        }
        aabb
    }
}

// ---------------------------------------------------------------------------
// MinkowskiSum
// ---------------------------------------------------------------------------

/// Minkowski sum of two convex shapes (for GJK distance computation).
///
/// Support function: h_A⊕B(d) = h_A(d) + h_B(d)
#[derive(Debug, Clone)]
pub struct MinkowskiSum {
    /// First convex shape.
    pub shape_a: GjkPolytope,
    /// Second convex shape (translated by offset_b).
    pub shape_b: GjkPolytope,
    /// Translation of shape B.
    pub offset_b: [f64; 3],
}

impl MinkowskiSum {
    /// Create a new Minkowski sum.
    pub fn new(shape_a: GjkPolytope, shape_b: GjkPolytope, offset_b: [f64; 3]) -> Self {
        Self {
            shape_a,
            shape_b,
            offset_b,
        }
    }

    /// Compute support of A ⊕ (-B) in direction d (for GJK).
    ///
    /// Support of Minkowski difference: s_A(d) - s_B(-d)
    pub fn support_difference(&self, d: [f64; 3]) -> [f64; 3] {
        let s_a = self.shape_a.support(d);
        let neg_d = [-d[0], -d[1], -d[2]];
        let s_b = self.shape_b.support(neg_d);
        // Translate B by offset_b
        let s_b_world = add3(s_b, self.offset_b);
        sub3(s_a, s_b_world)
    }

    /// Check if origin is contained in Minkowski difference (shapes overlap).
    pub fn origin_in_md(&self) -> bool {
        // Simple check: sample some directions
        let dirs = [
            [1.0, 0.0, 0.0_f64],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        // If origin is inside, all support points project positively
        for d in &dirs {
            let s = self.support_difference(*d);
            if dot3(s, *d) < 0.0 {
                return false;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// ShapeProxy
// ---------------------------------------------------------------------------

/// Thin wrapper for dispatching different shape types.
#[derive(Debug, Clone)]
pub enum ShapeProxy {
    /// Primitive shape.
    Primitive(PrimitiveShape),
    /// Compound shape.
    Compound(CompoundShape),
    /// Triangle mesh.
    TriangleMesh(TriangleMeshShape),
    /// GJK polytope.
    Polytope(GjkPolytope),
}

impl ShapeProxy {
    /// Compute AABB for this shape.
    pub fn aabb(&self) -> Aabb3 {
        match self {
            ShapeProxy::Primitive(p) => p.aabb(),
            ShapeProxy::Compound(c) => c.aabb(),
            ShapeProxy::TriangleMesh(m) => m.aabb,
            ShapeProxy::Polytope(p) => p.aabb(),
        }
    }

    /// Compute support point in direction d.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        match self {
            ShapeProxy::Primitive(p) => p.support(d),
            ShapeProxy::Compound(c) => c.support(d),
            ShapeProxy::TriangleMesh(m) => {
                let mut best = [0.0f64; 3];
                let mut best_proj = f64::NEG_INFINITY;
                for &v in &m.vertices {
                    let p = dot3(v, d);
                    if p > best_proj {
                        best_proj = p;
                        best = v;
                    }
                }
                best
            }
            ShapeProxy::Polytope(p) => p.support(d),
        }
    }

    /// Shape name string.
    pub fn name(&self) -> &str {
        match self {
            ShapeProxy::Primitive(PrimitiveShape::Sphere(_)) => "Sphere",
            ShapeProxy::Primitive(PrimitiveShape::Box(_)) => "Box",
            ShapeProxy::Primitive(PrimitiveShape::Capsule(_, _)) => "Capsule",
            ShapeProxy::Primitive(PrimitiveShape::Cylinder(_, _)) => "Cylinder",
            ShapeProxy::Primitive(PrimitiveShape::ConvexHull(_)) => "ConvexHull",
            ShapeProxy::Compound(_) => "Compound",
            ShapeProxy::TriangleMesh(_) => "TriangleMesh",
            ShapeProxy::Polytope(_) => "Polytope",
        }
    }
}

// ---------------------------------------------------------------------------
// CompoundBvh
// ---------------------------------------------------------------------------

/// BVH over compound sub-shapes for fast overlap tests.
#[derive(Debug, Clone)]
pub struct CompoundBvhNode {
    /// AABB for this node.
    pub aabb: Aabb3,
    /// Child indices (up to 2) or leaf shape index.
    pub children: [Option<usize>; 2],
    /// Leaf shape index (None if internal node).
    pub leaf_shape_idx: Option<usize>,
}

/// BVH over compound shapes.
#[derive(Debug, Clone)]
pub struct CompoundBvh {
    /// BVH nodes.
    pub nodes: Vec<CompoundBvhNode>,
    /// Root node index.
    pub root: usize,
    /// Shape proxies at leaves.
    pub shapes: Vec<(ShapeProxy, LocalTransform)>,
}

impl CompoundBvh {
    /// Build a BVH over a list of shapes.
    pub fn build(shapes: Vec<(ShapeProxy, LocalTransform)>) -> Self {
        let n = shapes.len();
        let mut nodes = Vec::new();

        if n == 0 {
            let empty_node = CompoundBvhNode {
                aabb: Aabb3::empty(),
                children: [None, None],
                leaf_shape_idx: None,
            };
            nodes.push(empty_node);
            return Self {
                nodes,
                root: 0,
                shapes,
            };
        }

        // Build a simple binary tree by AABB subdivision
        let shape_aabbs: Vec<Aabb3> = shapes
            .iter()
            .map(|(s, xform)| {
                let local_aabb = s.aabb();
                let mut world_aabb = Aabb3::empty();
                for corner in Self::aabb_corners(&local_aabb) {
                    world_aabb.include_point(xform.apply(corner));
                }
                world_aabb
            })
            .collect();

        let root = Self::build_recursive(&mut nodes, &(0..n).collect::<Vec<_>>(), &shape_aabbs);

        Self {
            nodes,
            root,
            shapes,
        }
    }

    fn aabb_corners(aabb: &Aabb3) -> [[f64; 3]; 8] {
        [
            [aabb.min[0], aabb.min[1], aabb.min[2]],
            [aabb.max[0], aabb.min[1], aabb.min[2]],
            [aabb.min[0], aabb.max[1], aabb.min[2]],
            [aabb.max[0], aabb.max[1], aabb.min[2]],
            [aabb.min[0], aabb.min[1], aabb.max[2]],
            [aabb.max[0], aabb.min[1], aabb.max[2]],
            [aabb.min[0], aabb.max[1], aabb.max[2]],
            [aabb.max[0], aabb.max[1], aabb.max[2]],
        ]
    }

    fn build_recursive(
        nodes: &mut Vec<CompoundBvhNode>,
        indices: &[usize],
        aabbs: &[Aabb3],
    ) -> usize {
        let mut total_aabb = Aabb3::empty();
        for &i in indices {
            total_aabb.include_aabb(&aabbs[i]);
        }

        if indices.len() == 1 {
            let node_idx = nodes.len();
            nodes.push(CompoundBvhNode {
                aabb: total_aabb,
                children: [None, None],
                leaf_shape_idx: Some(indices[0]),
            });
            return node_idx;
        }

        // Split by longest axis
        let extents = total_aabb.half_extents();
        let axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
            0
        } else if extents[1] >= extents[2] {
            1
        } else {
            2
        };
        let mid = total_aabb.center()[axis];
        let (left, right): (Vec<_>, Vec<_>) = indices
            .iter()
            .partition(|&&i| 0.5 * (aabbs[i].min[axis] + aabbs[i].max[axis]) <= mid);
        let left = if left.is_empty() {
            indices[..indices.len() / 2].to_vec()
        } else {
            left
        };
        let right = if right.is_empty() {
            indices[indices.len() / 2..].to_vec()
        } else {
            right
        };

        let left_idx = Self::build_recursive(nodes, &left, aabbs);
        let right_idx = Self::build_recursive(nodes, &right, aabbs);

        let node_idx = nodes.len();
        nodes.push(CompoundBvhNode {
            aabb: total_aabb,
            children: [Some(left_idx), Some(right_idx)],
            leaf_shape_idx: None,
        });
        node_idx
    }

    /// Query shapes that overlap a given AABB.
    pub fn query_aabb(&self, query: &Aabb3) -> Vec<usize> {
        let mut result = Vec::new();
        if self.nodes.is_empty() {
            return result;
        }
        self.query_recursive(self.root, query, &mut result);
        result
    }

    fn query_recursive(&self, node_idx: usize, query: &Aabb3, result: &mut Vec<usize>) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        if !node.aabb.overlaps(query) {
            return;
        }
        if let Some(leaf_idx) = node.leaf_shape_idx {
            result.push(leaf_idx);
            return;
        }
        for &child in &node.children {
            if let Some(c) = child {
                self.query_recursive(c, query, result);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ShapeOffset
// ---------------------------------------------------------------------------

/// Outward offset (inflation) of a shape by a radius.
///
/// Equivalent to Minkowski sum of shape with sphere of radius r.
#[derive(Debug, Clone)]
pub struct ShapeOffset {
    /// Inner shape.
    pub inner: ShapeProxy,
    /// Inflation radius.
    pub offset_radius: f64,
}

impl ShapeOffset {
    /// Create a new offset shape.
    pub fn new(inner: ShapeProxy, offset_radius: f64) -> Self {
        Self {
            inner,
            offset_radius,
        }
    }

    /// Compute support point: s_inner(d) + r * d̂.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        let n = norm3(d);
        let s = self.inner.support(d);
        if n < 1e-15 {
            return s;
        }
        add3(s, scale3(d, self.offset_radius / n))
    }

    /// Compute AABB (inflated by offset_radius).
    pub fn aabb(&self) -> Aabb3 {
        self.inner.aabb().inflate(self.offset_radius)
    }
}

// ---------------------------------------------------------------------------
// ImportedCollisionMesh
// ---------------------------------------------------------------------------

/// Collision mesh mode.
#[derive(Debug, Clone, PartialEq)]
pub enum CollisionMeshMode {
    /// Use as triangle mesh directly.
    TriangleMesh,
    /// Use convex hull.
    ConvexHull,
    /// Use approximate convex decomposition.
    ConvexDecomposition,
}

/// Maximum number of hulls for imported meshes in decomposition mode.
const DECOMPOSITION_MAX_HULLS: usize = 16;

/// Concavity threshold (mesh units) for imported meshes in decomposition mode.
const DECOMPOSITION_MAX_CONCAVITY: f64 = 0.01;

/// Surface sampling subdivisions for imported meshes in decomposition mode.
const DECOMPOSITION_SAMPLING: usize = 4;

/// Collision mesh imported from vertex/index data.
///
/// Can operate in triangle mesh, convex hull, or convex decomposition mode.
#[derive(Debug, Clone)]
pub struct ImportedCollisionMesh {
    /// Mesh mode.
    pub mode: CollisionMeshMode,
    /// Original vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices.
    pub indices: Vec<[usize; 3]>,
    /// Triangle mesh (for TriangleMesh mode).
    pub mesh: Option<TriangleMeshShape>,
    /// Convex hull vertices (for ConvexHull mode).
    pub convex_hull: Vec<[f64; 3]>,
    /// Decomposed hulls (for ConvexDecomposition mode).
    pub decomposed_hulls: Vec<Vec<[f64; 3]>>,
}

impl ImportedCollisionMesh {
    /// Create a new imported collision mesh.
    pub fn new(vertices: Vec<[f64; 3]>, indices: Vec<[usize; 3]>, mode: CollisionMeshMode) -> Self {
        let mut mesh = None;
        let mut convex_hull = Vec::new();
        let mut decomposed_hulls = Vec::new();

        match &mode {
            CollisionMeshMode::TriangleMesh => {
                mesh = Some(TriangleMeshShape::new(vertices.clone(), indices.clone()));
            }
            CollisionMeshMode::ConvexHull => {
                convex_hull = Self::compute_convex_hull(&vertices);
            }
            CollisionMeshMode::ConvexDecomposition => {
                let mut decomp = ConvexDecomposition::new(
                    DECOMPOSITION_MAX_HULLS,
                    DECOMPOSITION_MAX_CONCAVITY,
                    DECOMPOSITION_SAMPLING,
                );
                decomp.decompose(&vertices, &indices);
                decomposed_hulls = decomp.hulls;
            }
        }

        Self {
            mode,
            vertices,
            indices,
            mesh,
            convex_hull,
            decomposed_hulls,
        }
    }

    /// Compute the convex hull of the mesh vertices.
    ///
    /// Returns the corner points of the exact convex hull (QuickHull); vertices
    /// inside the hull or in the interior of one of its faces or edges are
    /// dropped. When the vertices do not span three dimensions (fewer than four
    /// points, or all coplanar) there is no volumetric hull and all distinct
    /// vertices are kept, which keeps [`ImportedCollisionMesh::support`] exact.
    /// Non-finite vertices are ignored.
    fn compute_convex_hull(verts: &[[f64; 3]]) -> Vec<[f64; 3]> {
        hull_vertices(verts)
    }

    /// Get AABB.
    pub fn aabb(&self) -> Aabb3 {
        match &self.mesh {
            Some(m) => m.aabb,
            None => {
                let mut aabb = Aabb3::empty();
                for &v in &self.vertices {
                    aabb.include_point(v);
                }
                aabb
            }
        }
    }

    /// Compute support point.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        let verts = match &self.mode {
            CollisionMeshMode::ConvexHull => &self.convex_hull,
            _ => &self.vertices,
        };
        if verts.is_empty() {
            return [0.0; 3];
        }
        verts
            .iter()
            .max_by(|a, b| {
                dot3(**a, d)
                    .partial_cmp(&dot3(**b, d))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or([0.0; 3])
    }
}

/// Compute voxelized collision mask (simplified).
///
/// Returns a flat boolean grid of occupied voxels.
pub fn voxelized_collision(mesh: &TriangleMeshShape, resolution: usize) -> Vec<bool> {
    let aabb = mesh.aabb;
    let n = resolution;
    let mut grid = vec![false; n * n * n];
    let query = MeshShapeQuery::new(mesh);

    for xi in 0..n {
        for yi in 0..n {
            for zi in 0..n {
                let p = [
                    aabb.min[0] + (xi as f64 + 0.5) / n as f64 * (aabb.max[0] - aabb.min[0]),
                    aabb.min[1] + (yi as f64 + 0.5) / n as f64 * (aabb.max[1] - aabb.min[1]),
                    aabb.min[2] + (zi as f64 + 0.5) / n as f64 * (aabb.max[2] - aabb.min[2]),
                ];
                grid[xi * n * n + yi * n + zi] = query.point_inside(p);
            }
        }
    }
    grid
}

/// Compute approximate shape volume using sampling.
pub fn shape_volume(shape: &ShapeProxy, resolution: usize) -> f64 {
    let aabb = shape.aabb();
    let vol = aabb.volume();
    if vol < 1e-20 || resolution == 0 {
        return vol;
    }
    let n = resolution;
    let mut inside_count = 0usize;
    let total = n * n * n;

    for xi in 0..n {
        for yi in 0..n {
            for zi in 0..n {
                let p = [
                    aabb.min[0] + (xi as f64 + 0.5) / n as f64 * (aabb.max[0] - aabb.min[0]),
                    aabb.min[1] + (yi as f64 + 0.5) / n as f64 * (aabb.max[1] - aabb.min[1]),
                    aabb.min[2] + (zi as f64 + 0.5) / n as f64 * (aabb.max[2] - aabb.min[2]),
                ];
                // Simple test: point inside support
                let s = shape.support(p);
                if dot3(s, p) > dot3(p, p) - 1e-6 {
                    inside_count += 1;
                }
            }
        }
    }
    vol * inside_count as f64 / total as f64
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::test_support::{cube_mesh, extrude_polygon, icosahedron, icosphere, signed_volume};
    use super::*;

    fn make_unit_sphere() -> PrimitiveShape {
        PrimitiveShape::Sphere(1.0)
    }

    fn make_box_shape() -> PrimitiveShape {
        PrimitiveShape::Box([1.0, 1.0, 1.0])
    }

    fn make_tetra() -> TriangleMeshShape {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let idxs = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]];
        TriangleMeshShape::new(verts, idxs)
    }

    #[test]
    fn test_aabb_volume() {
        let aabb = Aabb3::new([0.0; 3], [2.0, 3.0, 4.0]);
        assert!((aabb.volume() - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_overlaps() {
        let a = Aabb3::new([0.0; 3], [1.0; 3]);
        let b = Aabb3::new([0.5; 3], [1.5; 3]);
        assert!(a.overlaps(&b));
        let c = Aabb3::new([2.0; 3], [3.0; 3]);
        assert!(!a.overlaps(&c));
    }

    #[test]
    fn test_aabb_contains_point() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        assert!(aabb.contains_point([0.5, 0.5, 0.5]));
        assert!(!aabb.contains_point([1.5, 0.5, 0.5]));
    }

    #[test]
    fn test_sphere_support() {
        let sphere = make_unit_sphere();
        let s = sphere.support([1.0, 0.0, 0.0]);
        assert!(
            (s[0] - 1.0).abs() < 1e-10,
            "sphere support in +x should be [1,0,0]"
        );
    }

    #[test]
    fn test_box_support() {
        let b = make_box_shape();
        let s = b.support([1.0, 1.0, 1.0]);
        assert!((s[0] - 1.0).abs() < 1e-10);
        assert!((s[1] - 1.0).abs() < 1e-10);
        assert!((s[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_volume() {
        let s = make_unit_sphere();
        let v = s.volume();
        let expected = 4.0 / 3.0 * std::f64::consts::PI;
        assert!((v - expected).abs() < 1e-6);
    }

    #[test]
    fn test_compound_shape_add_child() {
        let mut compound = CompoundShape::new();
        compound.add_child(make_unit_sphere(), LocalTransform::identity());
        compound.add_child(
            make_box_shape(),
            LocalTransform::from_translation([3.0, 0.0, 0.0]),
        );
        assert_eq!(compound.num_children(), 2);
    }

    #[test]
    fn test_compound_aabb_contains_children() {
        let mut compound = CompoundShape::new();
        compound.add_child(
            PrimitiveShape::Sphere(1.0),
            LocalTransform::from_translation([0.0, 0.0, 0.0]),
        );
        compound.add_child(
            PrimitiveShape::Sphere(1.0),
            LocalTransform::from_translation([5.0, 0.0, 0.0]),
        );
        let aabb = compound.aabb();
        assert!(
            aabb.max[0] >= 6.0,
            "AABB should extend to 5+1=6 in x: {}",
            aabb.max[0]
        );
        assert!(
            aabb.min[0] <= -1.0,
            "AABB should extend to -1 in x: {}",
            aabb.min[0]
        );
    }

    #[test]
    fn test_compound_center_of_mass() {
        let mut compound = CompoundShape::new();
        compound.add_child(
            PrimitiveShape::Sphere(1.0),
            LocalTransform::from_translation([0.0, 0.0, 0.0]),
        );
        compound.add_child(
            PrimitiveShape::Sphere(1.0),
            LocalTransform::from_translation([2.0, 0.0, 0.0]),
        );
        let com = compound.center_of_mass();
        assert!(
            (com[0] - 1.0).abs() < 1e-6,
            "COM should be at x=1 for equal spheres at 0 and 2: {}",
            com[0]
        );
    }

    #[test]
    fn test_triangle_mesh_normals() {
        let mesh = make_tetra();
        assert_eq!(mesh.normals.len(), 4, "tetrahedron should have 4 normals");
        for n in &mesh.normals {
            let len = norm3(*n);
            assert!(
                (len - 1.0).abs() < 1e-6,
                "normal should be unit length: {}",
                len
            );
        }
    }

    #[test]
    fn test_triangle_mesh_surface_area() {
        // Flat unit triangle
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let idxs = vec![[0, 1, 2]];
        let mesh = TriangleMeshShape::new(verts, idxs);
        let area = mesh.surface_area();
        assert!(
            (area - 0.5).abs() < 1e-10,
            "unit right triangle area should be 0.5: {}",
            area
        );
    }

    #[test]
    fn test_closest_point_triangle() {
        let a = [0.0f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [2.0, 0.0, 0.0]; // outside triangle, closest to b
        let cp = MeshShapeQuery::closest_point_triangle(p, a, b, c);
        assert!(
            (cp[0] - 1.0).abs() < 1e-6,
            "closest point should be at b=[1,0,0]: {:?}",
            cp
        );
    }

    #[test]
    fn test_ray_triangle_intersect_hit() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let origin = [0.25, 0.25, -1.0];
        let dir = [0.0, 0.0, 1.0];
        let result = MeshShapeQuery::ray_triangle_intersect(origin, dir, a, b, c);
        assert!(result.is_some(), "ray should hit triangle");
    }

    #[test]
    fn test_ray_triangle_intersect_miss() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let origin = [2.0, 2.0, -1.0];
        let dir = [0.0, 0.0, 1.0];
        let result = MeshShapeQuery::ray_triangle_intersect(origin, dir, a, b, c);
        assert!(result.is_none(), "ray should miss triangle");
    }

    #[test]
    fn test_gjk_polytope_support() {
        let verts = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let poly = GjkPolytope::new(verts);
        let s = poly.support([1.0, 0.0, 0.0]);
        assert!((s[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_minkowski_sum_support() {
        let a = GjkPolytope::new(vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]);
        let b = GjkPolytope::new(vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]);
        let mk = MinkowskiSum::new(a, b, [3.0, 0.0, 0.0]);
        // s_A([1,0,0]) = [1,0,0], s_B([-1,0,0]) = [-0.5,0,0] (local) + [3,0,0] = [2.5,0,0]
        // support_diff = [1,0,0] - [2.5,0,0] = [-1.5,0,0]
        let s = mk.support_difference([1.0, 0.0, 0.0]);
        let _ = s; // Just check no panic
    }

    #[test]
    fn test_shape_proxy_name() {
        let proxy = ShapeProxy::Primitive(PrimitiveShape::Sphere(1.0));
        assert_eq!(proxy.name(), "Sphere");
        let box_proxy = ShapeProxy::Primitive(PrimitiveShape::Box([1.0; 3]));
        assert_eq!(box_proxy.name(), "Box");
    }

    #[test]
    fn test_compound_bvh_build() {
        let shapes = vec![
            (
                ShapeProxy::Primitive(PrimitiveShape::Sphere(1.0)),
                LocalTransform::from_translation([0.0, 0.0, 0.0]),
            ),
            (
                ShapeProxy::Primitive(PrimitiveShape::Sphere(1.0)),
                LocalTransform::from_translation([5.0, 0.0, 0.0]),
            ),
        ];
        let bvh = CompoundBvh::build(shapes);
        assert!(!bvh.nodes.is_empty());
    }

    #[test]
    fn test_compound_bvh_query() {
        let shapes = vec![
            (
                ShapeProxy::Primitive(PrimitiveShape::Sphere(1.0)),
                LocalTransform::from_translation([0.0, 0.0, 0.0]),
            ),
            (
                ShapeProxy::Primitive(PrimitiveShape::Sphere(1.0)),
                LocalTransform::from_translation([10.0, 0.0, 0.0]),
            ),
        ];
        let bvh = CompoundBvh::build(shapes);
        let query = Aabb3::new([-2.0; 3], [2.0; 3]);
        let result = bvh.query_aabb(&query);
        assert!(
            !result.is_empty(),
            "query near origin should find first sphere"
        );
    }

    #[test]
    fn test_shape_offset_inflates_aabb() {
        let inner = ShapeProxy::Primitive(PrimitiveShape::Sphere(1.0));
        let offset = ShapeOffset::new(inner, 0.5);
        let aabb = offset.aabb();
        assert!(
            aabb.max[0] >= 1.5,
            "offset should inflate AABB to at least 1.5"
        );
    }

    #[test]
    fn test_imported_mesh_convex_hull() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let idxs = vec![[0, 1, 2], [0, 1, 3]];
        let mesh = ImportedCollisionMesh::new(verts, idxs, CollisionMeshMode::ConvexHull);
        assert_eq!(
            mesh.convex_hull.len(),
            4,
            "all four tetrahedron corners are hull vertices"
        );
    }

    #[test]
    fn test_imported_mesh_triangle_mode() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let idxs = vec![[0, 1, 2]];
        let mesh = ImportedCollisionMesh::new(verts, idxs, CollisionMeshMode::TriangleMesh);
        assert!(mesh.mesh.is_some(), "TriangleMesh mode should build mesh");
    }

    #[test]
    fn test_convex_decomposition_basic() {
        let (verts, idxs) = cube_mesh(0.5, [0.5; 3]);
        let mut decomp = ConvexDecomposition::new(4, 0.1, 32);
        decomp.decompose(&verts, &idxs);
        assert_eq!(decomp.num_hulls(), 1, "a cube is already convex");
        assert_eq!(decomp.hulls[0].len(), 8);
    }

    #[test]
    fn test_imported_mesh_convex_hull_keeps_every_extreme_vertex() {
        // All 12 vertices of an icosahedron (and all 42 of the subdivided
        // icosphere) are extreme points; the old six-direction filter kept at
        // most six of them, so a sphere mesh collided as an octahedron.
        for (mesh, expected) in [(icosahedron(1.0), 12), (icosphere(1.0), 42)] {
            let imported =
                ImportedCollisionMesh::new(mesh.0.clone(), mesh.1, CollisionMeshMode::ConvexHull);
            assert_eq!(imported.convex_hull.len(), expected);
            for v in &mesh.0 {
                assert!(
                    imported
                        .convex_hull
                        .iter()
                        .any(|h| norm3(sub3(*h, *v)) < 1e-12),
                    "vertex {v:?} missing from the hull"
                );
            }
        }
    }

    #[test]
    fn test_convex_hull_does_not_depend_on_the_vertex_order() {
        // A sorted icosphere is the classic input on which an order-sensitive
        // hull builder produced an open, self-overlapping surface.
        let mesh = icosphere(2.0);
        let expected = signed_volume(&mesh);
        let mut sorted = mesh.0.clone();
        sorted.sort_by(|a, b| {
            a[0].total_cmp(&b[0])
                .then(a[1].total_cmp(&b[1]))
                .then(a[2].total_cmp(&b[2]))
        });
        let mut reversed = sorted.clone();
        reversed.reverse();
        for vertices in [mesh.0.clone(), sorted, reversed] {
            let shape = PrimitiveShape::ConvexHull(vertices.clone());
            assert!(
                (shape.volume() - expected).abs() < 1e-9 * expected,
                "volume {} vs {expected}",
                shape.volume()
            );
            let imported =
                ImportedCollisionMesh::new(vertices, Vec::new(), CollisionMeshMode::ConvexHull);
            assert_eq!(imported.convex_hull.len(), 42);
        }
    }

    #[test]
    fn test_imported_mesh_hull_drops_interior_and_face_points() {
        let (mut verts, idxs) = cube_mesh(1.0, [0.0; 3]);
        verts.push([0.0, 0.0, 0.0]); // interior
        verts.push([1.0, 0.0, 0.0]); // centre of the +x face
        verts.push([0.3, -0.2, 1.0]); // inside the +z face
        let mesh = ImportedCollisionMesh::new(verts, idxs, CollisionMeshMode::ConvexHull);
        assert_eq!(mesh.convex_hull.len(), 8);
    }

    #[test]
    fn test_imported_mesh_hull_ignores_duplicated_vertices() {
        // Meshes exported with one vertex per face corner repeat every cube
        // corner three times; the duplicates must not change the hull.
        let (corners, _) = cube_mesh(1.0, [0.0; 3]);
        let mut verts = Vec::new();
        for corner in &corners {
            verts.extend([*corner; 3]);
        }
        let mesh =
            ImportedCollisionMesh::new(verts.clone(), Vec::new(), CollisionMeshMode::ConvexHull);
        assert_eq!(mesh.convex_hull.len(), 8);
        let volume = PrimitiveShape::ConvexHull(verts).volume();
        assert!((volume - 8.0).abs() < 1e-9, "hull volume {volume}");
    }

    #[test]
    fn test_imported_mesh_support_matches_brute_force() {
        let (verts, idxs) = icosphere(2.0);
        let mesh = ImportedCollisionMesh::new(verts.clone(), idxs, CollisionMeshMode::ConvexHull);
        let mut seed = 12345_u64;
        let mut next = || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) as f64 / (1u64 << 31) as f64 * 2.0 - 1.0
        };
        let mut directions = vec![[1.0, 1.0, 1.0]];
        for _ in 0..50 {
            directions.push([next(), next(), next()]);
        }
        for d in directions {
            let best = verts
                .iter()
                .map(|v| dot3(*v, d))
                .fold(f64::NEG_INFINITY, f64::max);
            let support = mesh.support(d);
            assert!(
                (dot3(support, d) - best).abs() < 1e-12,
                "support is not extremal along {d:?}"
            );
        }
    }

    #[test]
    fn test_imported_mesh_flat_hull_keeps_all_points() {
        let quad = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let mesh = ImportedCollisionMesh::new(
            quad.clone(),
            vec![[0, 1, 2], [0, 2, 3]],
            CollisionMeshMode::ConvexHull,
        );
        assert_eq!(mesh.convex_hull.len(), 4);
        let support = mesh.support([1.0, 1.0, 0.0]);
        assert!((support[0] - 1.0).abs() < 1e-12 && (support[1] - 1.0).abs() < 1e-12);

        // Fewer than four points, duplicates and non-finite points.
        let degenerate = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [f64::NAN, 0.0, 0.0],
        ];
        let mesh =
            ImportedCollisionMesh::new(degenerate, Vec::new(), CollisionMeshMode::ConvexHull);
        assert_eq!(mesh.convex_hull, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let empty =
            ImportedCollisionMesh::new(Vec::new(), Vec::new(), CollisionMeshMode::ConvexHull);
        assert!(empty.convex_hull.is_empty());
    }

    #[test]
    fn test_imported_mesh_decomposition_mode_splits_concave_mesh() {
        let l_shape = [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let (verts, idxs) = extrude_polygon(&l_shape, 1.0);
        let mesh = ImportedCollisionMesh::new(verts, idxs, CollisionMeshMode::ConvexDecomposition);
        assert_eq!(mesh.decomposed_hulls.len(), 2);
        let volume: f64 = mesh
            .decomposed_hulls
            .iter()
            .map(|h| PrimitiveShape::ConvexHull(h.clone()).volume())
            .sum();
        assert!((volume - 3.0).abs() < 1e-6, "volume {volume}");
    }

    #[test]
    fn test_convex_hull_primitive_volume_is_the_hull_volume() {
        let cube = PrimitiveShape::ConvexHull(cube_mesh(1.0, [0.0; 3]).0);
        assert!((cube.volume() - 8.0).abs() < 1e-9);

        let tetra = PrimitiveShape::ConvexHull(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        // The bounding box would give 1.
        assert!((tetra.volume() - 1.0 / 6.0).abs() < 1e-12);

        let octahedron = PrimitiveShape::ConvexHull(vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ]);
        // The bounding box would give 8.
        assert!((octahedron.volume() - 4.0 / 3.0).abs() < 1e-12);

        let flat = PrimitiveShape::ConvexHull(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]);
        assert_eq!(flat.volume(), 0.0);
        assert_eq!(PrimitiveShape::ConvexHull(Vec::new()).volume(), 0.0);
    }

    #[test]
    fn test_compound_total_volume_uses_hull_volume() {
        let mut compound = CompoundShape::new();
        compound.add_child(
            PrimitiveShape::ConvexHull(cube_mesh(1.0, [0.0; 3]).0),
            LocalTransform::identity(),
        );
        compound.add_child(PrimitiveShape::Box([1.0; 3]), LocalTransform::identity());
        assert!((compound.total_volume() - 16.0).abs() < 1e-9);
    }

    #[test]
    fn test_merge_hulls() {
        let hulls = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        ];
        let merged = merge_hulls(&hulls);
        assert_eq!(merged.len(), 4, "merged should have all 4 vertices");
    }

    #[test]
    fn test_local_transform_apply() {
        let t = LocalTransform::from_translation([1.0, 2.0, 3.0]);
        let p = t.apply([0.0; 3]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_local_transform_round_trip() {
        let t = LocalTransform::from_translation([5.0, -3.0, 2.0]);
        let p = [1.0, 2.0, 3.0];
        let p2 = t.apply(p);
        let p3 = t.inverse_apply(p2);
        for i in 0..3 {
            assert!(
                (p3[i] - p[i]).abs() < 1e-10,
                "round trip failed at axis {}",
                i
            );
        }
    }
}
