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

use crate::Result;

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

mod acd_mesh;
#[cfg(test)]
mod acd_tests;

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
    /// The hull starts from a non-degenerate tetrahedron chosen from the
    /// whole point set (so leading coplanar points are harmless) and grows
    /// one point at a time, which costs O(n²) overall.  All input points are
    /// kept in [`ConvexHull3d::vertices`]; the surface is described by
    /// [`ConvexHull3d::triangles`] (outward, counter-clockwise).  The
    /// tolerances scale with the extent of the point set.
    ///
    /// Points that span no volume (fewer than four points, or all points
    /// coincident, collinear or coplanar, or a non-finite coordinate) have no
    /// hull surface: the result then has no triangles, faces or half-edges
    /// and zero volume.
    pub fn build(points: &[Vec3]) -> Self {
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

    /// Check that every point lies on the hull boundary: within `tol` of at
    /// least one facet plane, so that none is buried strictly inside the hull.
    ///
    /// Together with [`ConvexHull3d::all_points_inside_or_on`] this says that
    /// the points are in convex position.  A hull without facets (degenerate
    /// input) has no interior, so every point passes.
    pub fn all_points_on_boundary(&self, points: &[Vec3], tol: f64) -> bool {
        if self.faces.is_empty() {
            return true;
        }
        points.iter().all(|&p| {
            self.faces.iter().enumerate().any(|(fi, f)| {
                let tri = &self.triangles[fi];
                let v0 = self.vertices[tri[0]];
                dot(f.normal, sub(p, v0)).abs() <= tol
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

    /// Support point: the vertex maximising d·x (the origin for an empty hull).
    pub fn support_point(&self, direction: Vec3) -> Vec3 {
        self.vertices
            .iter()
            .max_by(|&&a, &&b| {
                dot(a, direction)
                    .partial_cmp(&dot(b, direction))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or([0.0; 3])
    }
}

/// Relative tolerance of the incremental hull.
///
/// A point only sees a facet when it lies further than this fraction of the
/// largest extent of the point set in front of the facet plane, and a seed
/// vertex must be at least that far from the line (plane) through the
/// previous seed vertices.
const HULL_RELATIVE_TOLERANCE: f64 = 1e-10;

/// Component-wise minimum and maximum of a non-empty point set, `None` if it
/// is empty or has a non-finite coordinate.
fn point_cloud_bounds(pts: &[Vec3]) -> Option<(Vec3, Vec3)> {
    let first = *pts.first()?;
    let mut lo = first;
    let mut hi = first;
    for p in pts {
        for ((low, high), &c) in lo.iter_mut().zip(hi.iter_mut()).zip(p.iter()) {
            if !c.is_finite() {
                return None;
            }
            *low = low.min(c);
            *high = high.max(c);
        }
    }
    Some((lo, hi))
}

/// Index and value of the point that maximises `distance`.
fn farthest(pts: &[Vec3], distance: impl Fn(Vec3) -> f64) -> Option<(usize, f64)> {
    pts.iter()
        .enumerate()
        .map(|(i, &p)| (i, distance(p)))
        .max_by(|x, y| x.1.total_cmp(&y.1))
}

/// Four affinely independent points to start the hull from, or `None` if the
/// points are coincident, collinear or coplanar within `tolerance`.
///
/// The first two are the extremes along `axis` (the axis of largest extent),
/// the third is the point farthest from the line through them and the fourth
/// the point farthest from the plane through the three.
fn seed_tetrahedron(pts: &[Vec3], axis: usize, tolerance: f64) -> Option<[usize; 4]> {
    let mut lowest = 0_usize;
    let mut highest = 0_usize;
    for (i, p) in pts.iter().enumerate() {
        if p[axis] < pts[lowest][axis] {
            lowest = i;
        }
        if p[axis] > pts[highest][axis] {
            highest = i;
        }
    }
    let origin = pts[lowest];
    let along = sub(pts[highest], origin);
    let along_length = norm(along);
    if lowest == highest || along_length <= tolerance {
        return None;
    }
    let (third, line_distance) =
        farthest(pts, |p| norm(cross(sub(p, origin), along)) / along_length)?;
    if line_distance <= tolerance {
        return None;
    }
    let normal = cross(along, sub(pts[third], origin));
    let normal_length = norm(normal);
    let (fourth, plane_distance) =
        farthest(pts, |p| dot(normal, sub(p, origin)).abs() / normal_length)?;
    if plane_distance <= tolerance {
        return None;
    }
    Some([lowest, highest, third, fourth])
}

/// Incremental hull: O(n²·f).  The seed tetrahedron is chosen from the whole
/// point set, so the order of the points (for example four leading coplanar
/// points) does not matter.  Returns outward counter-clockwise triangles, or
/// an empty list if the points span no volume.
fn incremental_hull(pts: &[Vec3]) -> Vec<[usize; 3]> {
    let Some((lo, hi)) = point_cloud_bounds(pts) else {
        return Vec::new();
    };
    let extents = sub(hi, lo);
    let axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
        0
    } else if extents[1] >= extents[2] {
        1
    } else {
        2
    };
    let extent = extents[axis];
    if extent <= 0.0 || !extent.is_finite() {
        return Vec::new();
    }
    let tolerance = HULL_RELATIVE_TOLERANCE * extent;
    let Some(seed) = seed_tetrahedron(pts, axis, tolerance) else {
        return Vec::new();
    };
    // The centroid of the seed tetrahedron lies strictly inside every hull
    // grown from it, so it fixes the outward side of every facet.
    let centroid = scale(seed.iter().fold([0.0; 3], |acc, &i| add(acc, pts[i])), 0.25);
    let [s0, s1, s2, s3] = seed;
    let mut tris: Vec<[usize; 3]> = vec![[s0, s1, s2], [s0, s3, s1], [s0, s2, s3], [s1, s3, s2]];
    for tri in &mut tris {
        orient_outward(pts, tri, centroid);
    }

    for (i, &p) in pts.iter().enumerate() {
        if seed.contains(&i) {
            continue;
        }
        // Facets that see the point.
        let visible: Vec<bool> = tris
            .iter()
            .map(|tri| dot(face_normal(pts, tri), sub(p, pts[tri[0]])) > tolerance)
            .collect();
        if !visible.contains(&true) {
            continue;
        }

        // Horizon: directed edges of visible facets whose neighbour across
        // the edge does not see the point.
        let mut horizon: Vec<(usize, usize)> = Vec::new();
        for (fi, tri) in tris.iter().enumerate() {
            if !visible[fi] {
                continue;
            }
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let twin_visible = tris
                    .iter()
                    .enumerate()
                    .any(|(fj, t)| fj != fi && visible[fj] && t.contains(&a) && t.contains(&b));
                if !twin_visible {
                    horizon.push((a, b));
                }
            }
        }

        // Replace the visible facets by a cone from the point to the horizon.
        let mut next: Vec<[usize; 3]> = tris
            .iter()
            .zip(visible.iter())
            .filter(|&(_, &sees)| !sees)
            .map(|(&t, _)| t)
            .collect();
        for (a, b) in horizon {
            // A point on the line through the edge cannot form a facet with it.
            let edge = sub(pts[b], pts[a]);
            if norm(cross(edge, sub(p, pts[a]))) <= tolerance * norm(edge) {
                continue;
            }
            let mut tri = [a, b, i];
            orient_outward(pts, &mut tri, centroid);
            next.push(tri);
        }
        tris = next;
    }
    tris
}

/// Unit normal of the triangle, following its winding (the zero vector for a
/// zero-area triangle).
///
/// Unlike [`normalize`] this has no absolute length threshold, so the result
/// does not depend on the scale of the point set.
fn face_normal(pts: &[Vec3], tri: &[usize; 3]) -> Vec3 {
    let a = pts[tri[0]];
    let b = pts[tri[1]];
    let c = pts[tri[2]];
    let normal = cross(sub(b, a), sub(c, a));
    let length = norm(normal);
    if length > 0.0 && length.is_finite() {
        scale(normal, 1.0 / length)
    } else {
        [0.0; 3]
    }
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
        for (k, &orig) in tri.iter().enumerate().take(3_usize) {
            half_edges.push(HalfEdge {
                origin: orig,
                twin: usize::MAX, // filled below
                next: base + (k + 1) % 3,
                face: fi,
            });
        }
        let normal = face_normal(verts, tri);
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
    ///
    /// A part is accepted once the deepest gap between its convex hull and
    /// the real surface is at most this large.
    pub max_concavity: f64,
    /// Maximum number of convex parts (at least 1)
    pub max_parts: usize,
    /// Minimum volume for a sub-part (ignore tiny fragments)
    ///
    /// Parts whose hull volume is below this fraction (in `[0, 1]`) of the
    /// volume of the convex hull of the whole mesh are discarded.
    pub min_volume_fraction: f64,
}

impl Default for AcdConfig {
    fn default() -> Self {
        Self {
            max_concavity: 0.01,
            max_parts: 16,
            min_volume_fraction: 0.001,
        }
    }
}

impl AcdConfig {
    /// Construct an AcdConfig with default values.
    pub fn new() -> Self {
        Self::default()
    }
}

/// A convex sub-part of a decomposition: the convex hull of one cluster.
#[derive(Debug, Clone)]
pub struct ConvexPart {
    /// Vertex positions of this sub-part (hull vertices for ACD)
    pub vertices: Vec<Vec3>,
    /// Concavity of this part (0 = perfectly convex)
    pub concavity: f64,
    /// Outward counter-clockwise hull triangles as indices into `vertices`
    /// (empty for parts built with [`ConvexPart::new`]).
    pub indices: Vec<[u32; 3]>,
    /// Volume of the hull of this part (0.0 for parts built with
    /// [`ConvexPart::new`]).
    pub volume: f64,
    /// Centroid of the hull solid of this part ([0, 0, 0] for parts built
    /// with [`ConvexPart::new`]).
    pub centroid: Vec3,
}

impl ConvexPart {
    /// Create a ConvexPart from bare vertices: no hull triangles, zero concavity,
    /// volume 0 and centroid at the origin.
    pub fn new(vertices: Vec<Vec3>) -> Self {
        let concavity = 0.0;
        Self {
            vertices,
            concavity,
            indices: Vec::new(),
            volume: 0.0,
            centroid: [0.0; 3],
        }
    }

    /// Check that the vertices form a convex polytope: every vertex lies on
    /// the boundary of the convex hull of all the vertices (within `tol`), so
    /// none of them is buried inside it.
    ///
    /// Fewer than four vertices, and vertex sets that span no volume, count
    /// as convex.
    pub fn is_convex(&self, tol: f64) -> bool {
        if self.vertices.len() < 4 {
            return true;
        }
        let hull = ConvexHull3d::build(&self.vertices);
        hull.all_points_inside_or_on(&self.vertices, tol)
            && hull.all_points_on_boundary(&self.vertices, tol)
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

/// Approximate Convex Decomposition: splits a closed triangle mesh into
/// parts that are each close to their own convex hull.
///
/// Concavity is measured against the mesh surface (how deep the hull of a
/// part bulges over the real surface), and the part with the largest
/// concavity is split by the best of several candidate planes until every
/// part is below [`AcdConfig::max_concavity`] or [`AcdConfig::max_parts`] is
/// reached; see [`ApproxConvexDecomp::decompose`].
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

    /// Decompose the triangle mesh `(vertices, triangles)` into convex parts.
    ///
    /// The mesh should be closed (watertight); either triangle winding is
    /// accepted.  A part of the decomposition is the piece of the solid cut
    /// out by the split planes chosen so far, represented by the mesh
    /// triangles clipped to it; the result replaces each part by its convex
    /// hull.
    ///
    /// The *concavity* of a part measures how far its convex hull sticks out
    /// beyond the part: it is the largest distance from a probe point on a
    /// hull facet to the part, where the probes form a barycentric grid on
    /// every facet (so the value is a sampled lower bound of the exact
    /// one-sided Hausdorff distance).  Probes that lie on the mesh surface or
    /// inside the solid, which includes the planar faces created by earlier
    /// cuts, are at distance zero.
    ///
    /// Parts are split, worst first, until every part is below
    /// [`AcdConfig::max_concavity`] or [`AcdConfig::max_parts`] is reached.
    /// A part is split by the candidate plane for which the concavities of
    /// the two children add up to the least (ties go to the smaller total
    /// hull volume).  The candidates have the coordinate axes or the
    /// principal axes of the part as normal and pass through the most concave
    /// point (the surface point nearest to the worst probe) or through the
    /// hull centroid; the planes of the mesh face nearest to the worst probe
    /// and of the faces that share a corner with it (the extensions of the
    /// walls of a notch) are tried as well.
    ///
    /// Every returned part is the convex hull of its cluster, with the hull
    /// triangles in [`ConvexPart::indices`] and the hull volume and centroid
    /// filled in.  Parts below [`AcdConfig::min_volume_fraction`] of the
    /// volume of the hull of the whole mesh are dropped.  A part whose
    /// [`ConvexPart::concavity`] still exceeds [`AcdConfig::max_concavity`]
    /// signals that `max_parts` was reached (or no plane could improve it).
    ///
    /// # Performance
    ///
    /// Every candidate plane of every split measures the concavity of both
    /// children against their clipped triangles, which costs time
    /// proportional to the number of hull facets times the number of
    /// triangles.  The method is meant for offline preparation of collision
    /// meshes with up to a few thousand triangles.
    ///
    /// # Errors
    ///
    /// * [`InvalidMesh`](crate::Error::InvalidMesh) for an empty mesh or a
    ///   non-finite vertex,
    /// * [`IndexOutOfBounds`](crate::Error::IndexOutOfBounds) for a triangle
    ///   index past `vertices`,
    /// * [`DegenerateGeometry`](crate::Error::DegenerateGeometry) if the mesh
    ///   spans no volume,
    /// * [`OutOfRange`](crate::Error::OutOfRange) for an invalid configuration
    ///   (`max_parts` of zero, negative or non-finite `max_concavity`, or
    ///   `min_volume_fraction` outside `[0, 1]`).
    pub fn decompose(&self, vertices: &[Vec3], triangles: &[[usize; 3]]) -> Result<ConvexParts> {
        acd_mesh::decompose_mesh(&self.config, vertices, triangles)
    }
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

    /// Support point: vertex realising the maximum (the origin when there are
    /// no vertices).
    pub fn support_point(&self, direction: Vec3) -> Vec3 {
        self.vertices
            .iter()
            .max_by(|&&a, &&b| {
                dot(a, direction)
                    .partial_cmp(&dot(b, direction))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or([0.0; 3])
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

/// GJK distance computation between two convex polytopes.
///
/// Implements the Gilbert–Johnson–Keerthi algorithm on the Minkowski
/// difference `A ⊖ B = { a − b }` of the two vertex sets.  The distance
/// between the hulls equals the distance from the origin to `A ⊖ B`.  Each
/// iteration keeps a simplex of at most four configuration-space points, finds
/// the point `v` of that simplex closest to the origin (dropping vertices that
/// carry zero weight), and asks for the support point of `A ⊖ B` in the
/// direction `−v`.  The loop stops when that support point cannot improve on
/// `v` (the lower bound meets the upper bound `|v|`), when the origin lies in
/// a tetrahedron (the hulls overlap) or when `max_iter` is reached.
#[derive(Debug, Clone)]
pub struct GjkDistance {
    /// Maximum number of iterations (support-point evaluations).
    ///
    /// If the limit is hit before convergence the current upper bound `|v|`
    /// on the distance is returned.
    pub max_iter: usize,
    /// Relative accuracy of the returned distance.
    ///
    /// The search stops when the GJK lower bound is within this fraction of
    /// the current distance estimate; shapes closer than `tolerance` times
    /// the size of their Minkowski difference are reported as touching
    /// (distance `0.0`).
    pub tolerance: f64,
}

impl Default for GjkDistance {
    fn default() -> Self {
        Self {
            max_iter: 64,
            tolerance: 1e-10,
        }
    }
}

impl GjkDistance {
    /// Create a GjkDistance instance.
    pub fn new(max_iter: usize, tolerance: f64) -> Self {
        Self {
            max_iter,
            tolerance,
        }
    }

    /// Compute the minimum Euclidean distance between the convex hulls of two
    /// vertex sets.
    ///
    /// * Returns `0.0` if the hulls overlap or touch (GJK cannot report the
    ///   penetration depth; use an EPA-style algorithm for that).
    /// * Returns `f64::INFINITY` if either vertex set is empty.
    /// * Returns `f64::NAN` if any vertex has a non-finite coordinate.
    pub fn distance(&self, a: &SupportFunction, b: &SupportFunction) -> f64 {
        if a.vertices.is_empty() || b.vertices.is_empty() {
            return f64::INFINITY;
        }
        let all_finite = a
            .vertices
            .iter()
            .chain(b.vertices.iter())
            .all(|v| v.iter().all(|c| c.is_finite()));
        if !all_finite {
            return f64::NAN;
        }
        // Support point of the Minkowski difference A ⊖ B in direction d:
        // s_A(d) − s_B(−d).
        let cso_support = |d: Vec3| -> Vec3 {
            sub(
                farthest_vertex(&a.vertices, d),
                farthest_vertex(&b.vertices, neg(d)),
            )
        };
        let tol = self.tolerance.max(0.0);
        let first = cso_support([1.0, 0.0, 0.0]);
        let mut simplex: [Vec3; 4] = [first; 4];
        let mut simplex_len = 1_usize;
        let mut closest = first;
        let mut extent_sq = dot(first, first);

        for _ in 0..self.max_iter {
            let closest_sq = dot(closest, closest);
            // The origin lies in (or on the boundary of) A ⊖ B.
            if closest_sq <= tol * tol * extent_sq {
                return 0.0;
            }
            let w = cso_support(neg(closest));
            extent_sq = extent_sq.max(dot(w, w));
            // v·w = min over A ⊖ B of v·x, so (|v|² − v·w) / |v| bounds the
            // remaining error of the estimate |v| from above.
            if closest_sq - dot(closest, w) <= tol * closest_sq {
                return closest_sq.sqrt();
            }
            // A support point that is already in the simplex cannot improve it.
            if simplex[..simplex_len].contains(&w) {
                return closest_sq.sqrt();
            }
            if simplex_len == simplex.len() {
                // A full tetrahedron survived the reduction: the origin is in it.
                return 0.0;
            }
            simplex[simplex_len] = w;
            simplex_len += 1;
            let (new_closest, kept, kept_len) = closest_point_on_simplex(&simplex[..simplex_len]);
            if dot(new_closest, new_closest) >= closest_sq {
                // Rounding stall: the estimate can no longer decrease.
                return closest_sq.sqrt();
            }
            simplex = kept;
            simplex_len = kept_len;
            closest = new_closest;
        }
        norm(closest)
    }
}

/// Vertex with the largest projection on `direction` (first one on ties).
fn farthest_vertex(vertices: &[Vec3], direction: Vec3) -> Vec3 {
    let mut best = vertices.first().copied().unwrap_or([0.0; 3]);
    let mut best_dot = dot(best, direction);
    for &v in vertices.iter().skip(1) {
        let d = dot(v, direction);
        if d > best_dot {
            best_dot = d;
            best = v;
        }
    }
    best
}

/// Point of the simplex `conv(points)` (at most four points) closest to the
/// origin, together with the minimal sub-simplex that carries it.
///
/// Every face of the simplex (every non-empty subset of the points) is tested:
/// the closest point of the face's affine hull is computed from the normal
/// equations and accepted when all barycentric weights are non-negative.  The
/// accepted candidate with the smallest norm is the closest point of the whole
/// simplex.  Faces are visited in order of increasing dimension, so a vertex
/// whose weight is zero is never kept.
///
/// Returns `(closest point, kept vertices, number of kept vertices)`.
fn closest_point_on_simplex(points: &[Vec3]) -> (Vec3, [Vec3; 4], usize) {
    let n = points.len().min(4);
    let mut best_point = points.first().copied().unwrap_or([0.0; 3]);
    let mut best_sq = dot(best_point, best_point);
    let mut best_mask: u8 = 1;
    for face_size in 1..=n {
        for mask in 1_u8..(1_u8 << n) {
            if mask.count_ones() as usize != face_size {
                continue;
            }
            if let Some(candidate) = closest_point_on_face(points, mask) {
                let candidate_sq = dot(candidate, candidate);
                if candidate_sq < best_sq {
                    best_sq = candidate_sq;
                    best_point = candidate;
                    best_mask = mask;
                }
            }
        }
    }
    let mut kept = [[0.0_f64; 3]; 4];
    let mut kept_len = 0_usize;
    for (i, &p) in points.iter().take(n).enumerate() {
        if best_mask & (1_u8 << i) != 0 {
            kept[kept_len] = p;
            kept_len += 1;
        }
    }
    (best_point, kept, kept_len)
}

/// Closest point to the origin of the affine hull of the face selected by
/// `mask` (bit `i` selects `points[i]`), if it lies inside the face.
///
/// Returns `None` when the selected points are affinely dependent (the normal
/// equations are singular) or when the closest point of the hull falls
/// outside the face (a negative barycentric weight).
fn closest_point_on_face(points: &[Vec3], mask: u8) -> Option<Vec3> {
    // Barycentric weights below this bound are rounding noise, not geometry.
    const WEIGHT_SLACK: f64 = 1e-12;
    let mut idx = [0_usize; 4];
    let mut k = 0_usize;
    for i in 0..points.len().min(4) {
        if mask & (1_u8 << i) != 0 {
            idx[k] = i;
            k += 1;
        }
    }
    if k == 0 {
        return None;
    }
    let base = points[idx[0]];
    if k == 1 {
        return Some(base);
    }
    // Minimise |base + Σ t_j e_j|² over the edge vectors e_j = p_j − base:
    // (EᵀE) t = −Eᵀ base.
    let m = k - 1;
    let mut edges = [[0.0_f64; 3]; 3];
    for (edge, &i) in edges.iter_mut().zip(idx[1..k].iter()) {
        *edge = sub(points[i], base);
    }
    let mut gram = [[0.0_f64; 3]; 3];
    let mut rhs = [0.0_f64; 3];
    for (i, &edge_i) in edges.iter().enumerate().take(m) {
        for (j, &edge_j) in edges.iter().enumerate().take(m) {
            gram[i][j] = dot(edge_i, edge_j);
        }
        rhs[i] = -dot(edge_i, base);
    }
    let t = solve_gram_system(&gram, &rhs, m)?;
    let sum_t: f64 = t[..m].iter().sum();
    let mut weights = [0.0_f64; 4];
    weights[0] = 1.0 - sum_t;
    weights[1..k].copy_from_slice(&t[..m]);
    if weights[..k].iter().any(|&w| w < -WEIGHT_SLACK) {
        return None;
    }
    // Clamp rounding noise and renormalise so the result is a genuine convex
    // combination of the face vertices (hence never closer than the true
    // closest point).
    let total: f64 = weights[..k].iter().map(|w| w.max(0.0)).sum();
    if total <= 0.0 {
        return None;
    }
    let mut point = [0.0_f64; 3];
    for (&w, &i) in weights[..k].iter().zip(idx[..k].iter()) {
        point = add(point, scale(points[i], w.max(0.0) / total));
    }
    Some(point)
}

/// Solve the `m × m` (m ≤ 3) symmetric positive semi-definite system
/// `gram · t = rhs` with Cramer's rule.  Returns `None` when the matrix is
/// numerically singular relative to the product of its diagonal (Hadamard
/// bound), i.e. when the edge vectors are (nearly) linearly dependent.
fn solve_gram_system(gram: &[[f64; 3]; 3], rhs: &[f64; 3], m: usize) -> Option<[f64; 3]> {
    // Smallest accepted ratio det(G) / Π G_ii; the ratio equals Π sin²(angle)
    // for orthogonal-ish edges, so 1e-10 rejects edges closer than ~1e-5 rad.
    const MIN_RELATIVE_DET: f64 = 1e-10;
    let mut t = [0.0_f64; 3];
    match m {
        1 => {
            let g00 = gram[0][0];
            if g00 <= 0.0 {
                return None;
            }
            t[0] = rhs[0] / g00;
        }
        2 => {
            let (g00, g01, g11) = (gram[0][0], gram[0][1], gram[1][1]);
            let det = g00 * g11 - g01 * g01;
            if g00 <= 0.0 || g11 <= 0.0 || det <= MIN_RELATIVE_DET * g00 * g11 {
                return None;
            }
            t[0] = (rhs[0] * g11 - g01 * rhs[1]) / det;
            t[1] = (g00 * rhs[1] - g01 * rhs[0]) / det;
        }
        3 => {
            let det = determinant_3x3(gram);
            let diag_product = gram[0][0] * gram[1][1] * gram[2][2];
            if diag_product <= 0.0 || det <= MIN_RELATIVE_DET * diag_product {
                return None;
            }
            for (col, slot) in t.iter_mut().enumerate() {
                // Cramer's rule: replace column `col` of G by the right-hand side.
                let mut replaced = *gram;
                for (row, replaced_row) in replaced.iter_mut().enumerate() {
                    replaced_row[col] = rhs[row];
                }
                *slot = determinant_3x3(&replaced) / det;
            }
        }
        _ => return None,
    }
    Some(t)
}

/// Determinant of a 3×3 matrix given as rows.
fn determinant_3x3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
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
    /// Returns true if they intersect (no separating axis found).  Both hulls
    /// should span a volume: a hull without facets (degenerate input to
    /// [`ConvexHull3d::build`]) contributes no separating axes, so the test
    /// then only uses the axes of the other hull and can report an
    /// intersection that does not exist.
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
        for (i, row) in inertia.iter_mut().enumerate() {
            row[i] = row[i].abs();
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
mod tests;
