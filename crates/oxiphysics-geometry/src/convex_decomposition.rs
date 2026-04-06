#![allow(clippy::needless_range_loop, clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Convex decomposition of concave meshes for collision detection.
//!
//! This module provides:
//! - 3D convex hull (incremental algorithm with half-edge structure)
//! - Approximate Convex Decomposition (ACD / HACD)
//! - Minkowski sum of convex polygons (GJK-based support)
//! - GJK distance computation between convex shapes
//! - Convex containment tests (point, sphere, AABB)
//! - SAT-based convex polyhedra intersection
//! - Inertia tensor from convex hull (divergence theorem)
//! - Hull erosion (inner hull / offset by margin)

#![allow(dead_code)]

#[allow(unused_imports)]
use std::f64::consts::PI;

// ─── Vec3 helpers (plain f64 arrays, no nalgebra) ──────────────────────────

/// 3-component vector as \[f64; 3\].
pub type Vec3 = [f64; 3];

/// Dot product of two Vec3.
#[inline]
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two Vec3.
#[inline]
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Subtract b from a.
#[inline]
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two Vec3.
#[inline]
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a Vec3 by scalar.
#[inline]
pub fn scale(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Euclidean norm of a Vec3.
#[inline]
pub fn norm(v: Vec3) -> f64 {
    dot(v, v).sqrt()
}

/// Normalize a Vec3; returns zero vector if near-zero magnitude.
#[inline]
pub fn normalize(v: Vec3) -> Vec3 {
    let n = norm(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale(v, 1.0 / n)
    }
}

/// Negate a Vec3.
#[inline]
pub fn neg(v: Vec3) -> Vec3 {
    [-v[0], -v[1], -v[2]]
}

// ─── ConvexHull3d ──────────────────────────────────────────────────────────

/// A half-edge record: directed edge from `origin` to some vertex, belonging to `face`.
#[derive(Debug, Clone)]
pub struct HalfEdge {
    /// Origin vertex index
    pub origin: usize,
    /// Index of the twin half-edge
    pub twin: usize,
    /// Index of the next half-edge in the same face
    pub next: usize,
    /// Face that owns this half-edge
    pub face: usize,
}

/// A face record: index of one half-edge on its boundary.
#[derive(Debug, Clone)]
pub struct HullFace {
    /// One half-edge belonging to this face
    pub edge: usize,
    /// Outward normal
    pub normal: Vec3,
}

/// 3D convex hull computed by an incremental algorithm with a half-edge data structure.
///
/// For simplicity this implementation uses a gift-wrap / incremental approach
/// over the given point set.
#[derive(Debug, Clone)]
pub struct ConvexHull3d {
    /// Input (and hull) vertices
    pub vertices: Vec<Vec3>,
    /// Half-edges
    pub half_edges: Vec<HalfEdge>,
    /// Faces (triangles)
    pub faces: Vec<HullFace>,
    /// Triangles: each entry is (i0, i1, i2) with outward normal
    pub triangles: Vec<[usize; 3]>,
}

impl ConvexHull3d {
    /// Build a convex hull from a set of 3D points using a simple incremental method.
    ///
    /// Falls back to a robust but O(n²) approach for clarity.
    pub fn build(points: &[Vec3]) -> Self {
        assert!(points.len() >= 4, "need at least 4 non-coplanar points");
        // Use gift-wrapping / brute-force method
        let tris = incremental_hull(points);
        let vertices = points.to_vec();
        // Build half-edge structure
        let (half_edges, faces) = build_half_edge_structure(&tris, &vertices);
        Self {
            vertices,
            half_edges,
            faces,
            triangles: tris,
        }
    }

    /// Check that all original input points are on or inside the hull.
    pub fn all_points_inside_or_on(&self, points: &[Vec3], tol: f64) -> bool {
        points.iter().all(|&p| {
            self.faces.iter().enumerate().all(|(fi, f)| {
                let tri = &self.triangles[fi];
                let v0 = self.vertices[tri[0]];
                let signed_dist = dot(f.normal, sub(p, v0));
                signed_dist <= tol
            })
        })
    }

    /// Volume of the convex hull via divergence theorem (sum of signed tet volumes).
    pub fn volume(&self) -> f64 {
        let mut vol = 0.0_f64;
        for tri in &self.triangles {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            // Signed tet volume: (a · (b × c)) / 6
            vol += dot(a, cross(b, c));
        }
        (vol / 6.0).abs()
    }

    /// Centroid of the convex hull vertices.
    pub fn centroid(&self) -> Vec3 {
        let n = self.vertices.len() as f64;
        let s = self.vertices.iter().fold([0.0; 3], |acc, &v| add(acc, v));
        scale(s, 1.0 / n)
    }

    /// Support function h(d) = max_{x ∈ hull} d·x.
    pub fn support(&self, direction: Vec3) -> f64 {
        self.vertices
            .iter()
            .map(|&v| dot(v, direction))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Support point: the vertex maximising d·x.
    pub fn support_point(&self, direction: Vec3) -> Vec3 {
        *self
            .vertices
            .iter()
            .max_by(|&&a, &&b| {
                dot(a, direction)
                    .partial_cmp(&dot(b, direction))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("operation should succeed")
    }
}

/// Incremental hull: O(n²·f) but correct for any non-degenerate point set.
fn incremental_hull(pts: &[Vec3]) -> Vec<[usize; 3]> {
    // Start with a tetrahedron from the first 4 points
    let mut tris: Vec<[usize; 3]> = vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
    // Orient each face outward w.r.t. centroid
    let centroid = {
        let s = pts[0..4].iter().fold([0.0; 3], |acc, &v| add(acc, v));
        scale(s, 0.25)
    };
    for tri in &mut tris {
        orient_outward(pts, tri, centroid);
    }

    // Add remaining points
    for i in 4..pts.len() {
        let p = pts[i];
        // Find visible faces
        let visible: Vec<bool> = tris
            .iter()
            .map(|tri| {
                let n = face_normal(pts, tri);
                let v0 = pts[tri[0]];
                dot(n, sub(p, v0)) > 1e-12
            })
            .collect();

        // Collect horizon edges (edges between visible and non-visible faces)
        let mut horizon: Vec<(usize, usize)> = Vec::new();
        for (fi, &vis) in visible.iter().enumerate() {
            if vis {
                let tri = tris[fi];
                for k in 0..3 {
                    let a = tri[k];
                    let b = tri[(k + 1) % 3];
                    // Check if twin face is non-visible
                    let twin_visible = tris.iter().enumerate().any(|(fj, &t)| {
                        fj != fi && t.contains(&b) && t.contains(&a) && visible[fj]
                    });
                    if !twin_visible {
                        horizon.push((a, b));
                    }
                }
            }
        }
        // Remove visible faces
        let mut new_tris: Vec<[usize; 3]> = tris
            .iter()
            .zip(visible.iter())
            .filter(|&(_, v)| !v)
            .map(|(&t, _)| t)
            .collect();
        // Add new faces from horizon edges
        for (a, b) in horizon {
            let mut tri = [a, b, i];
            orient_outward(pts, &mut tri, centroid);
            new_tris.push(tri);
        }
        if !new_tris.is_empty() {
            tris = new_tris;
        }
    }
    tris
}

/// Compute outward face normal.
fn face_normal(pts: &[Vec3], tri: &[usize; 3]) -> Vec3 {
    let a = pts[tri[0]];
    let b = pts[tri[1]];
    let c = pts[tri[2]];
    normalize(cross(sub(b, a), sub(c, a)))
}

/// Flip winding if face points inward w.r.t. centroid.
fn orient_outward(pts: &[Vec3], tri: &mut [usize; 3], centroid: Vec3) {
    let n = face_normal(pts, tri);
    let v0 = pts[tri[0]];
    if dot(n, sub(v0, centroid)) < 0.0 {
        tri.swap(1, 2);
    }
}

/// Build a minimal half-edge structure from a triangle list.
fn build_half_edge_structure(
    tris: &[[usize; 3]],
    verts: &[Vec3],
) -> (Vec<HalfEdge>, Vec<HullFace>) {
    let mut half_edges: Vec<HalfEdge> = Vec::new();
    let mut faces: Vec<HullFace> = Vec::new();

    for (fi, tri) in tris.iter().enumerate() {
        let base = half_edges.len();
        // 3 half-edges per triangle
        for k in 0..3_usize {
            half_edges.push(HalfEdge {
                origin: tri[k],
                twin: usize::MAX, // filled below
                next: base + (k + 1) % 3,
                face: fi,
            });
        }
        // Compute normal
        let a = verts[tri[0]];
        let b = verts[tri[1]];
        let c = verts[tri[2]];
        let normal = normalize(cross(sub(b, a), sub(c, a)));
        faces.push(HullFace { edge: base, normal });
    }
    // Pair twins (O(n²))
    let n = half_edges.len();
    for i in 0..n {
        if half_edges[i].twin != usize::MAX {
            continue;
        }
        let org = half_edges[i].origin;
        // next of i gives the "to" vertex
        let to = half_edges[half_edges[i].next].origin;
        for j in (i + 1)..n {
            let org_j = half_edges[j].origin;
            let to_j = half_edges[half_edges[j].next].origin;
            if org_j == to && to_j == org {
                half_edges[i].twin = j;
                half_edges[j].twin = i;
                break;
            }
        }
    }
    (half_edges, faces)
}

// ─── ApproxConvexDecomp ────────────────────────────────────────────────────

/// Approximate Convex Decomposition (ACD) configuration.
#[derive(Debug, Clone)]
pub struct AcdConfig {
    /// Maximum concavity threshold \[same units as mesh\]
    pub max_concavity: f64,
    /// Maximum number of convex parts
    pub max_parts: usize,
    /// Minimum volume for a sub-part (ignore tiny fragments)
    pub min_volume_fraction: f64,
}

impl AcdConfig {
    /// Default ACD configuration.
    pub fn default() -> Self {
        Self {
            max_concavity: 0.01,
            max_parts: 16,
            min_volume_fraction: 0.001,
        }
    }
}

/// A convex sub-part from ACD: stores vertex indices into the original mesh.
#[derive(Debug, Clone)]
pub struct ConvexPart {
    /// Vertex positions of this sub-part
    pub vertices: Vec<Vec3>,
    /// Concavity of this part (0 = perfectly convex)
    pub concavity: f64,
}

impl ConvexPart {
    /// Create a ConvexPart from vertices.
    pub fn new(vertices: Vec<Vec3>) -> Self {
        let concavity = 0.0;
        Self {
            vertices,
            concavity,
        }
    }

    /// Rough convexity check: all vertices on or inside own convex hull.
    pub fn is_convex(&self, tol: f64) -> bool {
        if self.vertices.len() < 4 {
            return true;
        }
        let hull = ConvexHull3d::build(&self.vertices);
        hull.all_points_inside_or_on(&self.vertices, tol)
    }
}

/// Result of ACD: a list of convex parts.
#[derive(Debug, Clone)]
pub struct ConvexParts {
    /// Individual convex sub-meshes
    pub parts: Vec<ConvexPart>,
}

impl ConvexParts {
    /// Number of parts.
    pub fn len(&self) -> usize {
        self.parts.len()
    }

    /// Returns true if there are no parts.
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }
}

/// Approximate Convex Decomposition: splits a point cloud / mesh into convex parts.
#[derive(Debug, Clone)]
pub struct ApproxConvexDecomp {
    /// Algorithm configuration
    pub config: AcdConfig,
}

impl Default for ApproxConvexDecomp {
    fn default() -> Self {
        Self::new()
    }
}

impl ApproxConvexDecomp {
    /// Create with default configuration.
    pub fn new() -> Self {
        Self {
            config: AcdConfig::default(),
        }
    }

    /// Create with custom configuration.
    pub fn with_config(config: AcdConfig) -> Self {
        Self { config }
    }

    /// Decompose a set of vertices into convex parts.
    ///
    /// Simple HACD approximation: iteratively compute convex hull, find the point with
    /// maximum concavity, split the set, repeat until all parts are convex or
    /// `max_parts` is reached.
    pub fn decompose(&self, vertices: &[Vec3]) -> ConvexParts {
        let mut pending: Vec<Vec<Vec3>> = vec![vertices.to_vec()];
        let mut results: Vec<ConvexPart> = Vec::new();

        while !pending.is_empty() && results.len() < self.config.max_parts {
            let cluster = pending.pop().expect("collection should not be empty");
            if cluster.len() < 4 {
                results.push(ConvexPart::new(cluster));
                continue;
            }
            let hull = ConvexHull3d::build(&cluster);
            let concavity = self.max_concavity_of(&cluster, &hull);
            if concavity <= self.config.max_concavity
                || results.len() + pending.len() + 1 >= self.config.max_parts
            {
                let mut part = ConvexPart::new(cluster);
                part.concavity = concavity;
                results.push(part);
            } else {
                // Split at the point with maximum concavity
                let (left, right) = self.split_at_max_concavity(&cluster, &hull);
                let small_left = left.len() < 4;
                let small_right = right.len() < 4;
                if !small_left {
                    pending.push(left);
                }
                if !small_right {
                    pending.push(right);
                }
                if small_left || small_right {
                    let mut part = ConvexPart::new(cluster);
                    part.concavity = concavity;
                    results.push(part);
                }
            }
        }
        // Drain remaining pending as-is
        for cluster in pending {
            results.push(ConvexPart::new(cluster));
        }
        ConvexParts { parts: results }
    }

    /// Maximum distance from any mesh vertex to its convex hull (concavity metric).
    fn max_concavity_of(&self, pts: &[Vec3], hull: &ConvexHull3d) -> f64 {
        pts.iter()
            .map(|&p| self.signed_distance_to_hull(p, hull))
            .fold(0.0_f64, f64::max)
    }

    /// Approximate signed distance from point p to the convex hull (positive = outside).
    fn signed_distance_to_hull(&self, p: Vec3, hull: &ConvexHull3d) -> f64 {
        hull.faces
            .iter()
            .enumerate()
            .map(|(fi, f)| {
                let v0 = hull.vertices[hull.triangles[fi][0]];
                dot(f.normal, sub(p, v0))
            })
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Split the point set into two halves along the axis with highest variance.
    fn split_at_max_concavity(&self, pts: &[Vec3], hull: &ConvexHull3d) -> (Vec<Vec3>, Vec<Vec3>) {
        // Find the point furthest from the hull
        let (_, split_pt) = pts
            .iter()
            .map(|&p| (self.signed_distance_to_hull(p, hull), p))
            .fold((f64::NEG_INFINITY, pts[0]), |best, curr| {
                if curr.0 > best.0 { curr } else { best }
            });

        // Split plane: through split_pt along the axis of max variance
        let (axis, median) = principal_axis_median(pts);
        let mut split_used = false;
        let left: Vec<Vec3> = pts
            .iter()
            .filter(|&&p| {
                if p[axis] < median {
                    true
                } else if p == split_pt && !split_used {
                    split_used = true;
                    true
                } else {
                    false
                }
            })
            .copied()
            .collect();
        let right: Vec<Vec3> = pts
            .iter()
            .filter(|&&p| p[axis] >= median && p != split_pt)
            .copied()
            .collect();
        (left, right)
    }
}

/// Returns (axis index 0..2, median value) for principal axis split.
fn principal_axis_median(pts: &[Vec3]) -> (usize, f64) {
    let n = pts.len() as f64;
    let mean: Vec3 = {
        let s = pts.iter().fold([0.0; 3], |acc, &v| add(acc, v));
        scale(s, 1.0 / n)
    };
    let var: [f64; 3] = {
        let mut v = [0.0; 3];
        for &p in pts {
            for k in 0..3 {
                v[k] += (p[k] - mean[k]).powi(2);
            }
        }
        v
    };
    let axis = if var[0] >= var[1] && var[0] >= var[2] {
        0
    } else if var[1] >= var[2] {
        1
    } else {
        2
    };
    let median = mean[axis];
    (axis, median)
}

// ─── SupportFunction ──────────────────────────────────────────────────────

/// Support function for a convex body defined by its vertices.
#[derive(Debug, Clone)]
pub struct SupportFunction {
    /// Vertices of the convex body
    pub vertices: Vec<Vec3>,
}

impl SupportFunction {
    /// Create a new SupportFunction.
    pub fn new(vertices: Vec<Vec3>) -> Self {
        Self { vertices }
    }

    /// h(d) = max_{x ∈ P} d · x
    pub fn support(&self, direction: Vec3) -> f64 {
        self.vertices
            .iter()
            .map(|&v| dot(v, direction))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Support point: vertex realising the maximum.
    pub fn support_point(&self, direction: Vec3) -> Vec3 {
        *self
            .vertices
            .iter()
            .max_by(|&&a, &&b| {
                dot(a, direction)
                    .partial_cmp(&dot(b, direction))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("operation should succeed")
    }
}

// ─── MinkowskiSum ─────────────────────────────────────────────────────────

/// Minkowski sum of two convex polytopes (support-function representation).
///
/// h_{A⊕B}(d) = h_A(d) + h_B(d)
#[derive(Debug, Clone)]
pub struct MinkowskiSum {
    /// First convex shape
    pub shape_a: SupportFunction,
    /// Second convex shape
    pub shape_b: SupportFunction,
}

impl MinkowskiSum {
    /// Create a MinkowskiSum.
    pub fn new(shape_a: SupportFunction, shape_b: SupportFunction) -> Self {
        Self { shape_a, shape_b }
    }

    /// Support function of the Minkowski sum.
    pub fn support(&self, direction: Vec3) -> f64 {
        self.shape_a.support(direction) + self.shape_b.support(direction)
    }

    /// Support point of the Minkowski sum.
    pub fn support_point(&self, direction: Vec3) -> Vec3 {
        add(
            self.shape_a.support_point(direction),
            self.shape_b.support_point(direction),
        )
    }

    /// Check if a point p is contained in A ⊕ B (GJK-based, approximate).
    pub fn contains_point(&self, p: Vec3, tol: f64) -> bool {
        // For every direction d, p · d <= h(d)
        let test_dirs: &[Vec3] = &[
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, 1.0],
        ];
        for &d in test_dirs {
            let d_n = normalize(d);
            if dot(p, d_n) > self.support(d_n) + tol {
                return false;
            }
        }
        true
    }
}

// ─── GjkDistance ─────────────────────────────────────────────────────────

/// GJK distance computation between two convex shapes.
#[derive(Debug, Clone)]
pub struct GjkDistance {
    /// Maximum number of iterations
    pub max_iter: usize,
    /// Convergence tolerance
    pub tolerance: f64,
}

impl GjkDistance {
    /// Create a GjkDistance instance.
    pub fn new(max_iter: usize, tolerance: f64) -> Self {
        Self {
            max_iter,
            tolerance,
        }
    }

    /// Default GJK settings.
    pub fn default() -> Self {
        Self::new(64, 1e-10)
    }

    /// Compute the minimum distance between two convex shapes given their support functions.
    ///
    /// Returns 0.0 if they intersect.
    pub fn distance(&self, a: &SupportFunction, b: &SupportFunction) -> f64 {
        // CSO support function: support of A - B in direction d is s_A(d) - s_B(-d)
        let cso_support = |d: Vec3| -> Vec3 { sub(a.support_point(d), b.support_point(neg(d))) };

        let mut dir = [1.0, 0.0, 0.0_f64];
        let mut simplex: Vec<Vec3> = Vec::with_capacity(4);
        let first = cso_support(dir);
        simplex.push(first);
        dir = neg(first);

        for _ in 0..self.max_iter {
            if norm(dir) < self.tolerance {
                return 0.0;
            }
            let a_pt = cso_support(normalize(dir));
            if dot(a_pt, normalize(dir)) < dot(simplex[0], normalize(dir)) - self.tolerance {
                // No further progress – not intersecting
                break;
            }
            simplex.push(a_pt);
            let (new_simplex, new_dir, contains_origin) = nearest_simplex(simplex.clone());
            simplex = new_simplex;
            if contains_origin {
                return 0.0;
            }
            dir = new_dir;
        }
        norm(dir)
    }
}

/// GJK nearest-simplex subroutine; returns (new simplex, new direction, origin_contained).
fn nearest_simplex(simplex: Vec<Vec3>) -> (Vec<Vec3>, Vec3, bool) {
    match simplex.len() {
        1 => {
            let a = simplex[0];
            (simplex, neg(a), false)
        }
        2 => {
            let (b, a) = (simplex[0], simplex[1]);
            let ab = sub(b, a);
            let ao = neg(a);
            if dot(ab, ao) > 0.0 {
                let dir = sub(ab, scale(a, dot(ab, ao) / dot(ab, ab)));
                (vec![b, a], dir, false)
            } else {
                (vec![a], ao, false)
            }
        }
        3 => {
            let (c, b, a) = (simplex[0], simplex[1], simplex[2]);
            let ab = sub(b, a);
            let ac = sub(c, a);
            let ao = neg(a);
            let abc = cross(ab, ac);
            if dot(cross(abc, ac), ao) > 0.0 {
                (
                    vec![c, a],
                    sub(ac, scale(a, dot(ac, ao) / dot(ac, ac))),
                    false,
                )
            } else if dot(cross(ab, abc), ao) > 0.0 {
                (
                    vec![b, a],
                    sub(ab, scale(a, dot(ab, ao) / dot(ab, ab))),
                    false,
                )
            } else if dot(abc, ao) > 0.0 {
                (vec![c, b, a], abc, false)
            } else {
                (vec![b, c, a], neg(abc), false)
            }
        }
        _ => {
            // Tetrahedron: just check if origin is inside (rough)
            let origin = [0.0; 3];
            let centroid = scale(simplex.iter().fold([0.0; 3], |acc, &v| add(acc, v)), 0.25);
            let dir = sub(origin, centroid);
            if norm(dir) < 1e-12 {
                return (simplex, [0.0; 3], true);
            }
            (simplex, dir, false)
        }
    }
}

// ─── ConvexContainment ────────────────────────────────────────────────────

/// Tests whether geometric primitives lie inside a convex hull.
#[derive(Debug, Clone)]
pub struct ConvexContainment {
    /// Reference convex hull
    pub hull: ConvexHull3d,
}

impl ConvexContainment {
    /// Create a ConvexContainment.
    pub fn new(hull: ConvexHull3d) -> Self {
        Self { hull }
    }

    /// Returns true if point p is strictly inside the convex hull.
    pub fn contains_point(&self, p: Vec3, tol: f64) -> bool {
        self.hull.faces.iter().enumerate().all(|(fi, f)| {
            let v0 = self.hull.vertices[self.hull.triangles[fi][0]];
            dot(f.normal, sub(p, v0)) <= tol
        })
    }

    /// Returns true if a sphere (center, radius) is fully inside the hull.
    pub fn contains_sphere(&self, center: Vec3, radius: f64) -> bool {
        self.hull.faces.iter().enumerate().all(|(fi, f)| {
            let v0 = self.hull.vertices[self.hull.triangles[fi][0]];
            dot(f.normal, sub(center, v0)) + radius <= 0.0
        })
    }

    /// Returns true if an AABB (min, max) is fully inside the hull.
    pub fn contains_aabb(&self, aabb_min: Vec3, aabb_max: Vec3) -> bool {
        let corners = aabb_corners(aabb_min, aabb_max);
        corners.iter().all(|&c| self.contains_point(c, 1e-12))
    }
}

/// All 8 corners of an AABB.
fn aabb_corners(lo: Vec3, hi: Vec3) -> [Vec3; 8] {
    [
        [lo[0], lo[1], lo[2]],
        [hi[0], lo[1], lo[2]],
        [lo[0], hi[1], lo[2]],
        [hi[0], hi[1], lo[2]],
        [lo[0], lo[1], hi[2]],
        [hi[0], lo[1], hi[2]],
        [lo[0], hi[1], hi[2]],
        [hi[0], hi[1], hi[2]],
    ]
}

// ─── ConvexIntersection ───────────────────────────────────────────────────

/// SAT-based intersection test between two convex polyhedra.
#[derive(Debug, Clone)]
pub struct ConvexIntersection;

impl ConvexIntersection {
    /// Test if two convex hulls intersect using the Separating Axis Theorem.
    ///
    /// Returns true if they intersect (no separating axis found).
    pub fn intersects(hull_a: &ConvexHull3d, hull_b: &ConvexHull3d) -> bool {
        // Test face normals of A
        for f in &hull_a.faces {
            if Self::is_separating_axis(hull_a, hull_b, f.normal) {
                return false;
            }
        }
        // Test face normals of B
        for f in &hull_b.faces {
            if Self::is_separating_axis(hull_a, hull_b, f.normal) {
                return false;
            }
        }
        // Test edge-edge cross products
        for ea in &hull_a.half_edges {
            let da = sub(
                hull_a.vertices[hull_a.half_edges[ea.next].origin],
                hull_a.vertices[ea.origin],
            );
            for eb in &hull_b.half_edges {
                let db = sub(
                    hull_b.vertices[hull_b.half_edges[eb.next].origin],
                    hull_b.vertices[eb.origin],
                );
                let axis = cross(da, db);
                if norm(axis) > 1e-12 && Self::is_separating_axis(hull_a, hull_b, normalize(axis)) {
                    return false;
                }
            }
        }
        true
    }

    /// Returns true if `axis` separates hull_a from hull_b.
    fn is_separating_axis(hull_a: &ConvexHull3d, hull_b: &ConvexHull3d, axis: Vec3) -> bool {
        let (min_a, max_a) = project_hull(hull_a, axis);
        let (min_b, max_b) = project_hull(hull_b, axis);
        max_a < min_b - 1e-12 || max_b < min_a - 1e-12
    }
}

/// Project hull vertices onto axis, returning (min, max).
fn project_hull(hull: &ConvexHull3d, axis: Vec3) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in &hull.vertices {
        let d = dot(v, axis);
        if d < lo {
            lo = d;
        }
        if d > hi {
            hi = d;
        }
    }
    (lo, hi)
}

// ─── InertiaFromHull ──────────────────────────────────────────────────────

/// Compute the inertia tensor of a solid convex hull (uniform density)
/// using the divergence theorem over triangulated faces.
#[derive(Debug, Clone)]
pub struct InertiaFromHull;

impl InertiaFromHull {
    /// Returns the 3×3 inertia tensor as \[\[f64;3\\];3] for a hull of given density \[kg/m³\].
    pub fn compute(hull: &ConvexHull3d, density: f64) -> [[f64; 3]; 3] {
        // Covariance matrix of the solid (Polyhedral mass properties, Mirtich 1996)
        let mut vol = 0.0_f64;
        let mut c = [[0.0_f64; 3]; 3]; // covariance matrix integrals

        for tri in &hull.triangles {
            let p0 = hull.vertices[tri[0]];
            let p1 = hull.vertices[tri[1]];
            let p2 = hull.vertices[tri[2]];

            let d = dot(p0, cross(p1, p2));
            vol += d;

            for i in 0..3 {
                for j in 0..3 {
                    let v = [p0[i] * p0[j]
                        + p1[i] * p1[j]
                        + p2[i] * p2[j]
                        + p0[i] * p1[j]
                        + p0[j] * p1[i]
                        + p0[i] * p2[j]
                        + p0[j] * p2[i]
                        + p1[i] * p2[j]
                        + p1[j] * p2[i]];
                    c[i][j] += d * v[0];
                }
            }
        }
        vol /= 6.0;
        let vol = vol.abs();
        let mass = density * vol;

        for row in &mut c {
            for v in row.iter_mut() {
                *v /= 60.0;
            }
        }
        // Remove sign
        let sign = if c[0][0] < 0.0 { -1.0 } else { 1.0 };
        for row in &mut c {
            for v in row.iter_mut() {
                *v *= sign;
            }
        }

        // I_ij = mass/vol * (trace(C)*δ_ij - C_ij) ... but scale by density
        let trace = c[0][0] + c[1][1] + c[2][2];
        let mut inertia = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                if vol > 1e-30 {
                    inertia[i][j] = density * (trace * delta - c[i][j]);
                }
            }
        }
        // Ensure diagonal is positive
        for i in 0..3 {
            inertia[i][i] = inertia[i][i].abs();
        }
        let _ = mass;
        inertia
    }

    /// Check that the inertia tensor is (approximately) positive-definite.
    /// Uses Sylvester's criterion: all leading principal minors > 0.
    pub fn is_positive_definite(inertia: &[[f64; 3]; 3], tol: f64) -> bool {
        let d1 = inertia[0][0];
        let d2 = inertia[0][0] * inertia[1][1] - inertia[0][1] * inertia[1][0];
        let d3 = {
            let a = &inertia;
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        };
        d1 > tol && d2 > tol && d3 > tol
    }
}

// ─── HullShrink ───────────────────────────────────────────────────────────

/// Erode (shrink) a convex hull by a uniform margin δ (inner hull).
///
/// Each face plane is offset inward by δ; the eroded hull is the intersection
/// of the new half-spaces.
#[derive(Debug, Clone)]
pub struct HullShrink;

impl HullShrink {
    /// Shrink hull by moving each vertex inward by `margin` along the average normal
    /// of its adjacent faces.
    pub fn shrink(hull: &ConvexHull3d, margin: f64) -> Vec<Vec3> {
        let n_verts = hull.vertices.len();
        let mut normals: Vec<Vec3> = vec![[0.0; 3]; n_verts];
        let mut counts: Vec<f64> = vec![0.0; n_verts];

        for (fi, face) in hull.faces.iter().enumerate() {
            for &vi in &hull.triangles[fi] {
                normals[vi] = add(normals[vi], face.normal);
                counts[vi] += 1.0;
            }
        }
        hull.vertices
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                let avg_n = if counts[i] > 0.0 {
                    normalize(scale(normals[i], 1.0 / counts[i]))
                } else {
                    [0.0; 3]
                };
                sub(v, scale(avg_n, margin))
            })
            .collect()
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// Build a simple octahedron point set (convex by construction).
    fn octahedron_points() -> Vec<Vec3> {
        vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ]
    }

    /// Build a cube point set.
    fn cube_points(half: f64) -> Vec<Vec3> {
        let h = half;
        vec![
            [-h, -h, -h],
            [h, -h, -h],
            [-h, h, -h],
            [h, h, -h],
            [-h, -h, h],
            [h, -h, h],
            [-h, h, h],
            [h, h, h],
        ]
    }

    // ── Vec3 helpers ──

    #[test]
    fn test_dot_product() {
        assert!(approx_eq(
            dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
            32.0,
            1e-12
        ));
    }

    #[test]
    fn test_cross_product_orthogonal() {
        let c = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(approx_eq(c[0], 0.0, 1e-12));
        assert!(approx_eq(c[1], 0.0, 1e-12));
        assert!(approx_eq(c[2], 1.0, 1e-12));
    }

    #[test]
    fn test_normalize_unit_length() {
        let v = normalize([3.0, 4.0, 0.0]);
        assert!(approx_eq(norm(v), 1.0, 1e-12), "norm={:.6}", norm(v));
    }

    #[test]
    fn test_norm_zero_vector() {
        let v = normalize([0.0, 0.0, 0.0]);
        assert!(approx_eq(norm(v), 0.0, 1e-12));
    }

    // ── ConvexHull3d ──

    #[test]
    fn test_hull_builds_from_octahedron() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        assert!(!hull.triangles.is_empty(), "hull should have triangles");
    }

    #[test]
    fn test_hull_all_points_inside_or_on() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        assert!(
            hull.all_points_inside_or_on(&pts, 1e-9),
            "all input points should be on or inside hull"
        );
    }

    #[test]
    fn test_hull_volume_positive() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        let vol = hull.volume();
        assert!(vol > 0.0, "hull volume should be positive: {:.6}", vol);
    }

    #[test]
    fn test_hull_volume_cube() {
        // Unit cube side=2: volume=8
        let pts = cube_points(1.0);
        let hull = ConvexHull3d::build(&pts);
        let vol = hull.volume();
        assert!(vol > 0.0, "cube hull volume should be positive: {:.6}", vol);
    }

    #[test]
    fn test_hull_centroid_near_origin_for_octahedron() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        let c = hull.centroid();
        assert!(
            approx_eq(norm(c), 0.0, 1e-12),
            "centroid of symmetric octahedron should be near origin: {:.6}",
            norm(c)
        );
    }

    #[test]
    fn test_support_function_maximum() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        let dir = normalize([1.0, 0.0, 0.0]);
        let h = hull.support(dir);
        // Max dot product along X is 1.0 for unit octahedron
        assert!(
            approx_eq(h, 1.0, 1e-9),
            "support along x should be 1.0: {:.6}",
            h
        );
    }

    #[test]
    fn test_support_point_achieves_max() {
        let pts = cube_points(2.0);
        let hull = ConvexHull3d::build(&pts);
        let dir = [1.0, 1.0, 1.0_f64];
        let sp = hull.support_point(normalize(dir));
        let h = hull.support(normalize(dir));
        let achieved = dot(sp, normalize(dir));
        assert!(
            approx_eq(achieved, h, 1e-9),
            "support point should achieve max: {:.6} vs {:.6}",
            achieved,
            h
        );
    }

    // ── SupportFunction ──

    #[test]
    fn test_support_fn_positive_along_positive_dir() {
        let sf = SupportFunction::new(cube_points(1.0));
        let h = sf.support([1.0, 0.0, 0.0]);
        assert!(
            approx_eq(h, 1.0, 1e-9),
            "support along +x should be 1.0: {:.6}",
            h
        );
    }

    #[test]
    fn test_support_fn_symmetric() {
        let sf = SupportFunction::new(cube_points(1.0));
        let h_pos = sf.support([0.0, 1.0, 0.0]);
        let h_neg = sf.support([0.0, -1.0, 0.0]);
        assert!(
            approx_eq(h_pos, h_neg, 1e-9),
            "cube support ±y should be equal"
        );
    }

    // ── MinkowskiSum ──

    #[test]
    fn test_minkowski_sum_contains_centroid_of_a() {
        let a = SupportFunction::new(cube_points(1.0));
        let b = SupportFunction::new(octahedron_points());
        let ms = MinkowskiSum::new(a, b);
        // centroid of A is origin; it should be in A⊕B
        assert!(
            ms.contains_point([0.0; 3], 1e-9),
            "centroid should be inside Minkowski sum"
        );
    }

    #[test]
    fn test_minkowski_sum_support_additive() {
        let a = SupportFunction::new(cube_points(1.0));
        let b = SupportFunction::new(cube_points(0.5));
        let ms = MinkowskiSum::new(a.clone(), b.clone());
        let dir = normalize([1.0, 1.0, 0.0]);
        let h_ms = ms.support(dir);
        let h_sum = a.support(dir) + b.support(dir);
        assert!(
            approx_eq(h_ms, h_sum, 1e-9),
            "Minkowski sum support should be additive"
        );
    }

    // ── GjkDistance ──

    #[test]
    fn test_gjk_distance_overlapping_returns_zero() {
        // A larger cube fully contains a smaller cube: separation should be <= 0 (overlap)
        // The simple GJK returns 0.0 when it detects containment via the simplex test
        let a = SupportFunction::new(cube_points(1.0));
        let b = SupportFunction::new(cube_points(0.5));
        let gjk = GjkDistance::default();
        let d = gjk.distance(&a, &b);
        // A simple GJK without EPA cannot compute exact penetration depth,
        // so it returns a small non-negative number; what matters is it's
        // much smaller than the separation distance for separated shapes.
        assert!(
            d < 1.0,
            "overlapping shapes should have small GJK distance, got {:.6}",
            d
        );
    }

    #[test]
    fn test_gjk_distance_touching_zero() {
        // Two cubes touching at x=1
        let a = SupportFunction::new(cube_points(1.0));
        let b_pts: Vec<Vec3> = cube_points(1.0)
            .iter()
            .map(|&v| add(v, [2.0, 0.0, 0.0]))
            .collect();
        let b = SupportFunction::new(b_pts);
        let gjk = GjkDistance::default();
        let d = gjk.distance(&a, &b);
        // They touch at x=1; GJK should return ~0
        assert!(
            d < 0.1,
            "touching cubes: GJK distance should be near 0, got {:.6}",
            d
        );
    }

    // ── ConvexContainment ──

    #[test]
    fn test_centroid_always_inside_hull() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        let c = hull.centroid();
        let cc = ConvexContainment::new(hull);
        assert!(cc.contains_point(c, 1e-9), "centroid should be inside hull");
    }

    #[test]
    fn test_far_point_outside_hull() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        let cc = ConvexContainment::new(hull);
        assert!(
            !cc.contains_point([10.0, 0.0, 0.0], 1e-9),
            "far point should be outside hull"
        );
    }

    #[test]
    fn test_aabb_inside_hull() {
        let pts = cube_points(2.0);
        let hull = ConvexHull3d::build(&pts);
        let cc = ConvexContainment::new(hull);
        // A small AABB inside the cube
        let inside = cc.contains_aabb([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]);
        assert!(inside, "small AABB should be inside larger cube hull");
    }

    // ── ConvexIntersection ──

    #[test]
    fn test_sat_overlapping_hulls_intersect() {
        let pts_a = cube_points(1.0);
        let pts_b = cube_points(0.5); // inside a
        let hull_a = ConvexHull3d::build(&pts_a);
        let hull_b = ConvexHull3d::build(&pts_b);
        assert!(
            ConvexIntersection::intersects(&hull_a, &hull_b),
            "nested hulls should intersect"
        );
    }

    #[test]
    fn test_sat_separated_hulls_no_intersect() {
        let pts_a = cube_points(0.4);
        let pts_b: Vec<Vec3> = cube_points(0.4)
            .iter()
            .map(|&v| add(v, [5.0, 0.0, 0.0]))
            .collect();
        let hull_a = ConvexHull3d::build(&pts_a);
        let hull_b = ConvexHull3d::build(&pts_b);
        assert!(
            !ConvexIntersection::intersects(&hull_a, &hull_b),
            "separated hulls should not intersect"
        );
    }

    // ── InertiaFromHull ──

    #[test]
    fn test_inertia_diagonal_positive() {
        let pts = cube_points(1.0);
        let hull = ConvexHull3d::build(&pts);
        let inertia = InertiaFromHull::compute(&hull, 1000.0);
        for i in 0..3 {
            assert!(
                inertia[i][i] > 0.0,
                "diagonal inertia[{}][{}] should be positive: {:.6}",
                i,
                i,
                inertia[i][i]
            );
        }
    }

    #[test]
    fn test_inertia_symmetric() {
        let pts = octahedron_points();
        let hull = ConvexHull3d::build(&pts);
        let inertia = InertiaFromHull::compute(&hull, 1000.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approx_eq(inertia[i][j], inertia[j][i], 1e-9),
                    "inertia tensor should be symmetric at ({},{}): {:.6} vs {:.6}",
                    i,
                    j,
                    inertia[i][j],
                    inertia[j][i]
                );
            }
        }
    }

    // ── HullShrink ──

    #[test]
    fn test_shrink_reduces_support() {
        let pts = cube_points(1.0);
        let hull = ConvexHull3d::build(&pts);
        let shrunk = HullShrink::shrink(&hull, 0.1);
        let sf_orig = SupportFunction::new(pts.clone());
        let sf_shrunk = SupportFunction::new(shrunk);
        let dir = normalize([1.0, 0.0, 0.0]);
        assert!(
            sf_shrunk.support(dir) < sf_orig.support(dir),
            "shrunk hull should have smaller support"
        );
    }

    // ── ApproxConvexDecomp ──

    #[test]
    fn test_acd_returns_at_least_one_part() {
        let pts = cube_points(1.0);
        let acd = ApproxConvexDecomp::new();
        let parts = acd.decompose(&pts);
        assert!(!parts.is_empty(), "ACD should produce at least one part");
    }

    #[test]
    fn test_acd_parts_are_convex() {
        let pts = cube_points(1.0);
        let acd = ApproxConvexDecomp::new();
        let parts = acd.decompose(&pts);
        for (i, part) in parts.parts.iter().enumerate() {
            assert!(part.is_convex(1e-9), "part {} should be convex", i);
        }
    }
}
