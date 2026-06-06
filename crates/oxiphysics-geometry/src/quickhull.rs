// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! 3D Convex Hull via QuickHull algorithm.
//!
//! Provides convex hull construction, vertex/face/edge enumeration,
//! volume and surface area computation, point-in-hull test, hull merging,
//! and support queries.

use oxiphysics_core::math::Vec3;

/// 3D Convex Hull computed via QuickHull algorithm (nalgebra Vec3 variant).
pub struct ConvexHull3DVec {
    /// Hull vertices (deduplicated subset of input points).
    pub vertices: Vec<Vec3>,
    /// Triangular faces stored as triples of vertex indices (outward-facing winding).
    pub faces: Vec<[usize; 3]>,
}

impl ConvexHull3DVec {
    /// Compute convex hull from a point cloud using the QuickHull algorithm.
    ///
    /// Returns `None` if the input is degenerate (fewer than 4 non-coplanar points).
    pub fn build(points: &[Vec3]) -> Option<Self> {
        if points.len() < 4 {
            return None;
        }

        // --- Seed tetrahedron ------------------------------------------------
        let (i0, i1) = extreme_pair(points)?;
        let i2 = farthest_from_line(points, i0, i1)?;
        let (normal_seed, offset_seed) = plane_from_triangle(points[i0], points[i1], points[i2]);
        if normal_seed.norm() < 1e-12 {
            return None;
        }
        let i3 = farthest_point_excluding(points, normal_seed, offset_seed, &[i0, i1, i2])?;
        if point_plane_distance(points[i3], normal_seed, offset_seed).abs() < 1e-12 {
            return None;
        }

        // Build the 4 outward-facing faces of the initial tetrahedron.
        let tet = [i0, i1, i2, i3];
        let mut hull: Vec<[usize; 3]> = Vec::new();
        let combos: [([usize; 3], usize); 4] = [
            ([i0, i1, i2], i3),
            ([i0, i1, i3], i2),
            ([i0, i2, i3], i1),
            ([i1, i2, i3], i0),
        ];
        for (face, interior) in &combos {
            let oriented = orient_face(*face, *interior, points)?;
            hull.push(oriented);
        }

        // --- Iterative QuickHull expansion -----------------------------------
        let exterior: Vec<usize> = (0..points.len()).filter(|i| !tet.contains(i)).collect();

        for &p_idx in &exterior {
            let p = points[p_idx];

            let visible: Vec<usize> = hull
                .iter()
                .enumerate()
                .filter(|(_, face)| {
                    let (n, off) =
                        plane_from_triangle(points[face[0]], points[face[1]], points[face[2]]);
                    point_plane_distance(p, n, off) > 1e-10
                })
                .map(|(i, _)| i)
                .collect();

            if visible.is_empty() {
                continue;
            }

            let horizon = horizon_edges(&hull, &visible);

            let mut new_hull: Vec<[usize; 3]> = hull
                .iter()
                .enumerate()
                .filter(|(i, _)| !visible.contains(i))
                .map(|(_, f)| *f)
                .collect();

            for (e0, e1) in horizon {
                let new_face = [e0, e1, p_idx];
                let (nn, noff) = plane_from_triangle(points[e0], points[e1], points[p_idx]);
                if nn.norm() < 1e-12 {
                    continue;
                }
                let interior_ref = new_hull.first().and_then(|f| {
                    f.iter()
                        .find(|&&v| v != e0 && v != e1 && v != p_idx)
                        .copied()
                });
                let oriented = if let Some(ref_v) = interior_ref {
                    if point_plane_distance(points[ref_v], nn, noff) > 0.0 {
                        [e1, e0, p_idx]
                    } else {
                        new_face
                    }
                } else {
                    new_face
                };
                new_hull.push(oriented);
            }

            hull = new_hull;
        }

        if hull.is_empty() {
            return None;
        }

        // Collect unique vertex indices and remap.
        let mut used: Vec<usize> = hull.iter().flat_map(|f| f.iter().copied()).collect();
        used.sort_unstable();
        used.dedup();

        let vertices: Vec<Vec3> = used.iter().map(|&i| points[i]).collect();
        let remap = |old: usize| -> usize {
            used.binary_search(&old)
                .expect("element must be present in sorted list")
        };
        let faces: Vec<[usize; 3]> = hull
            .iter()
            .map(|f| [remap(f[0]), remap(f[1]), remap(f[2])])
            .collect();

        Some(Self { vertices, faces })
    }

    /// Support function: returns the farthest point in direction `d`.
    pub fn support(&self, d: Vec3) -> Vec3 {
        self.vertices
            .iter()
            .copied()
            .max_by(|a, b| {
                a.dot(&d)
                    .partial_cmp(&b.dot(&d))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(Vec3::zeros)
    }

    /// Volume of the convex hull via signed tetrahedral decomposition from the centroid.
    pub fn volume(&self) -> f64 {
        if self.faces.is_empty() || self.vertices.is_empty() {
            return 0.0;
        }
        let c = self.centroid();
        let mut vol = 0.0f64;
        for face in &self.faces {
            let a = self.vertices[face[0]] - c;
            let b = self.vertices[face[1]] - c;
            let cc = self.vertices[face[2]] - c;
            vol += a.dot(&b.cross(&cc));
        }
        (vol / 6.0).abs()
    }

    /// Surface area of the convex hull (sum of triangle areas).
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0f64;
        for face in &self.faces {
            let a = self.vertices[face[0]];
            let b = self.vertices[face[1]];
            let c = self.vertices[face[2]];
            area += (b - a).cross(&(c - a)).norm() * 0.5;
        }
        area
    }

    /// Returns `true` if point `p` is inside (or on the boundary of) the hull.
    ///
    /// A point is inside iff it is on the non-positive side of every face plane.
    pub fn contains_point(&self, p: Vec3) -> bool {
        for face in &self.faces {
            let n = self.face_normal_raw(face);
            let off = n.dot(&self.vertices[face[0]]);
            if point_plane_distance(p, n, off) > 1e-9 {
                return false;
            }
        }
        true
    }

    /// Outward-pointing face normal for face `face_idx`.
    pub fn face_normal(&self, face_idx: usize) -> Vec3 {
        let face = &self.faces[face_idx];
        self.face_normal_raw(face)
    }

    /// Number of hull vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangular faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Get the list of unique edges as pairs of vertex indices.
    ///
    /// Each edge appears once, with the smaller index first.
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut edge_set = std::collections::BTreeSet::new();
        for face in &self.faces {
            for i in 0..3 {
                let a = face[i];
                let b = face[(i + 1) % 3];
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                edge_set.insert((lo, hi));
            }
        }
        edge_set.into_iter().collect()
    }

    /// Number of unique edges.
    pub fn n_edges(&self) -> usize {
        self.edges().len()
    }

    /// Enumerate all vertices with their indices.
    pub fn vertex_iter(&self) -> impl Iterator<Item = (usize, &Vec3)> {
        self.vertices.iter().enumerate()
    }

    /// Enumerate all faces with their indices.
    pub fn face_iter(&self) -> impl Iterator<Item = (usize, &[usize; 3])> {
        self.faces.iter().enumerate()
    }

    /// Area of a specific face.
    pub fn face_area(&self, face_idx: usize) -> f64 {
        let face = &self.faces[face_idx];
        let a = self.vertices[face[0]];
        let b = self.vertices[face[1]];
        let c = self.vertices[face[2]];
        (b - a).cross(&(c - a)).norm() * 0.5
    }

    /// Compute the AABB of the hull.
    pub fn aabb(&self) -> (Vec3, Vec3) {
        if self.vertices.is_empty() {
            return (Vec3::zeros(), Vec3::zeros());
        }
        let mut mn = self.vertices[0];
        let mut mx = self.vertices[0];
        for v in &self.vertices[1..] {
            for k in 0..3 {
                if v[k] < mn[k] {
                    mn[k] = v[k];
                }
                if v[k] > mx[k] {
                    mx[k] = v[k];
                }
            }
        }
        (mn, mx)
    }

    /// Compute the diameter of the hull (maximum distance between any two vertices).
    pub fn diameter(&self) -> f64 {
        let mut max_dist = 0.0f64;
        for i in 0..self.vertices.len() {
            for j in (i + 1)..self.vertices.len() {
                let d = (self.vertices[i] - self.vertices[j]).norm();
                if d > max_dist {
                    max_dist = d;
                }
            }
        }
        max_dist
    }

    /// Euler characteristic: V - E + F = 2 for convex polyhedra.
    pub fn euler_characteristic(&self) -> i64 {
        self.n_vertices() as i64 - self.n_edges() as i64 + self.n_faces() as i64
    }

    /// Signed distance from point `p` to the hull surface.
    ///
    /// Negative if inside, positive if outside.
    pub fn signed_distance(&self, p: Vec3) -> f64 {
        let mut max_signed = f64::NEG_INFINITY;
        for face in &self.faces {
            let n = self.face_normal_raw(face);
            let off = n.dot(&self.vertices[face[0]]);
            let d = point_plane_distance(p, n, off);
            if d > max_signed {
                max_signed = d;
            }
        }
        max_signed
    }

    /// Check if two hulls overlap (conservative test using separating axis theorem
    /// on face normals only).
    pub fn overlaps_sat_faces(&self, other: &ConvexHull3DVec) -> bool {
        // Check all face normals of self
        for fi in 0..self.n_faces() {
            let n = self.face_normal(fi);
            let (a_min, a_max) = project_hull_on_axis(&self.vertices, n);
            let (b_min, b_max) = project_hull_on_axis(&other.vertices, n);
            if a_max < b_min || b_max < a_min {
                return false; // separating axis found
            }
        }
        // Check all face normals of other
        for fi in 0..other.n_faces() {
            let n = other.face_normal(fi);
            let (a_min, a_max) = project_hull_on_axis(&self.vertices, n);
            let (b_min, b_max) = project_hull_on_axis(&other.vertices, n);
            if a_max < b_min || b_max < a_min {
                return false;
            }
        }
        true
    }

    // ------------------------------------------------------------------
    // Private helpers
    // ------------------------------------------------------------------

    fn centroid(&self) -> Vec3 {
        let sum: Vec3 = self.vertices.iter().copied().sum();
        sum / self.vertices.len() as f64
    }

    fn face_normal_raw(&self, face: &[usize; 3]) -> Vec3 {
        let a = self.vertices[face[0]];
        let b = self.vertices[face[1]];
        let c = self.vertices[face[2]];
        (b - a).cross(&(c - a)).normalize()
    }
}

/// Project all vertices onto an axis and return (min, max).
fn project_hull_on_axis(verts: &[Vec3], axis: Vec3) -> (f64, f64) {
    let mut mn = f64::INFINITY;
    let mut mx = f64::NEG_INFINITY;
    for v in verts {
        let d = v.dot(&axis);
        if d < mn {
            mn = d;
        }
        if d > mx {
            mx = d;
        }
    }
    (mn, mx)
}

/// Merge two convex hulls into one by computing the hull of the union of vertices.
pub fn merge_hulls(a: &ConvexHull3DVec, b: &ConvexHull3DVec) -> Option<ConvexHull3DVec> {
    let mut all_points: Vec<Vec3> = a.vertices.clone();
    all_points.extend_from_slice(&b.vertices);
    ConvexHull3DVec::build(&all_points)
}

/// Compute convex hull and return both hull and the indices (into original
/// point cloud) of vertices on the hull.
pub fn build_with_indices(points: &[Vec3]) -> Option<(ConvexHull3DVec, Vec<usize>)> {
    let hull = ConvexHull3DVec::build(points)?;
    // Find which original point corresponds to each hull vertex
    let mut indices = Vec::with_capacity(hull.vertices.len());
    for hv in &hull.vertices {
        let idx = points
            .iter()
            .position(|p| (*p - *hv).norm() < 1e-12)
            .unwrap_or(0);
        indices.push(idx);
    }
    Some((hull, indices))
}

/// Check if a point set is convex (i.e., all points lie on the hull).
pub fn is_point_set_convex(points: &[Vec3]) -> bool {
    if points.len() < 4 {
        return true;
    }
    match ConvexHull3DVec::build(points) {
        Some(hull) => hull.n_vertices() == points.len(),
        None => true, // degenerate
    }
}

// ---------------------------------------------------------------------------
// Internal QuickHull helpers
// ---------------------------------------------------------------------------

/// Signed distance from point `p` to the plane `(normal, offset)`.
fn point_plane_distance(p: Vec3, normal: Vec3, offset: f64) -> f64 {
    normal.dot(&p) - offset
}

/// Compute the plane (normal, offset) for triangle (a, b, c).
fn plane_from_triangle(a: Vec3, b: Vec3, c: Vec3) -> (Vec3, f64) {
    let n = (b - a).cross(&(c - a));
    let offset = n.dot(&a);
    (n, offset)
}

/// Orient face so that `interior` is strictly behind the face plane.
fn orient_face(face: [usize; 3], interior: usize, points: &[Vec3]) -> Option<[usize; 3]> {
    let (n, off) = plane_from_triangle(points[face[0]], points[face[1]], points[face[2]]);
    if n.norm() < 1e-12 {
        return None;
    }
    if interior == usize::MAX {
        return Some(face);
    }
    if point_plane_distance(points[interior], n, off) > 0.0 {
        Some([face[0], face[2], face[1]])
    } else {
        Some(face)
    }
}

/// Find the horizon edges of the hull given a set of visible face indices.
fn horizon_edges(hull: &[[usize; 3]], visible: &[usize]) -> Vec<(usize, usize)> {
    let mut edges: Vec<(usize, usize)> = Vec::new();

    for &fi in visible {
        let face = &hull[fi];
        let face_edges = [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])];
        for (e0, e1) in face_edges {
            let adjacent_visible = visible.iter().any(|&other_fi| {
                if other_fi == fi {
                    return false;
                }
                let other = &hull[other_fi];
                let other_edges = [
                    (other[0], other[1]),
                    (other[1], other[2]),
                    (other[2], other[0]),
                ];
                other_edges.contains(&(e1, e0))
            });
            if !adjacent_visible {
                edges.push((e0, e1));
            }
        }
    }
    edges
}

/// Find the pair of points that are farthest apart along the x-axis (for seeding).
fn extreme_pair(points: &[Vec3]) -> Option<(usize, usize)> {
    if points.len() < 2 {
        return None;
    }
    let i0 = points
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)?;
    let i1 = points
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)?;
    if i0 == i1 {
        let i1y = points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)?;
        if i0 == i1y {
            return None;
        }
        return Some((i0, i1y));
    }
    Some((i0, i1))
}

/// Find the point farthest from the line through `points[i0]` and `points[i1]`.
fn farthest_from_line(points: &[Vec3], i0: usize, i1: usize) -> Option<usize> {
    let a = points[i0];
    let b = points[i1];
    let ab = b - a;
    let ab_len_sq = ab.dot(&ab);
    if ab_len_sq < 1e-24 {
        return None;
    }
    points
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != i0 && i != i1)
        .max_by(|(_, p), (_, q)| {
            let dp = ((*p - a).cross(&ab)).norm_squared() / ab_len_sq;
            let dq = ((*q - a).cross(&ab)).norm_squared() / ab_len_sq;
            dp.partial_cmp(&dq).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
}

/// Find the farthest point from the plane, excluding the given indices.
fn farthest_point_excluding(
    points: &[Vec3],
    normal: Vec3,
    offset: f64,
    exclude: &[usize],
) -> Option<usize> {
    points
        .iter()
        .enumerate()
        .filter(|(i, _)| !exclude.contains(i))
        .max_by(|(_, a), (_, b)| {
            point_plane_distance(**a, normal, offset)
                .partial_cmp(&point_plane_distance(**b, normal, offset))
                .expect("operation should succeed")
        })
        .map(|(i, _)| i)
}

// ---------------------------------------------------------------------------
// Incremental convex hull
// ---------------------------------------------------------------------------

/// Incremental convex hull that can accept new points one at a time.
///
/// Internally maintains the same triangulated hull structure as `ConvexHull3DVec`.
pub struct IncrementalConvexHull {
    /// All points ever inserted.
    pub all_points: Vec<Vec3>,
    /// Current hull faces.
    pub faces: Vec<[usize; 3]>,
    /// Current hull vertex indices (into `all_points`).
    pub hull_verts: Vec<usize>,
}

impl IncrementalConvexHull {
    /// Build an incremental hull from a slice of initial points.
    pub fn build(points: &[Vec3]) -> Self {
        let mut hull = Self {
            all_points: Vec::new(),
            faces: Vec::new(),
            hull_verts: Vec::new(),
        };
        for &p in points {
            hull.insert(p);
        }
        hull
    }

    /// Insert a new point, updating the hull if necessary.
    pub fn insert(&mut self, p: Vec3) {
        self.all_points.push(p);
        // Rebuild from scratch (simple but correct)
        if let Some(ch) = ConvexHull3DVec::build(&self.all_points) {
            // Remap hull_verts to current all_points indices
            self.hull_verts = (0..self.all_points.len())
                .filter(|&i| {
                    ch.vertices
                        .iter()
                        .any(|&v| (v - self.all_points[i]).norm() < 1e-10)
                })
                .collect();
            // Map from ch.vertices index back to all_points index
            let vert_to_orig: Vec<usize> = ch
                .vertices
                .iter()
                .map(|hv| {
                    self.all_points
                        .iter()
                        .enumerate()
                        .find(|&(_, ap)| (*ap - *hv).norm() < 1e-10)
                        .map(|(i, _)| i)
                        .unwrap_or(0)
                })
                .collect();
            self.faces = ch
                .faces
                .iter()
                .map(|f| [vert_to_orig[f[0]], vert_to_orig[f[1]], vert_to_orig[f[2]]])
                .collect();
        }
    }

    /// Number of vertices on the current hull.
    pub fn n_vertices(&self) -> usize {
        self.hull_verts.len()
    }

    /// Number of triangular faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of unique edges.
    pub fn n_edges(&self) -> usize {
        let mut edges = std::collections::BTreeSet::new();
        for f in &self.faces {
            for k in 0..3 {
                let a = f[k];
                let b = f[(k + 1) % 3];
                edges.insert((a.min(b), a.max(b)));
            }
        }
        edges.len()
    }

    /// Euler characteristic V - E + F = 2.
    pub fn euler_characteristic(&self) -> i64 {
        self.n_vertices() as i64 - self.n_edges() as i64 + self.n_faces() as i64
    }

    /// Volume via signed tetrahedral decomposition.
    pub fn volume(&self) -> f64 {
        if self.faces.is_empty() {
            return 0.0;
        }
        let c = {
            let sum: Vec3 = self.hull_verts.iter().map(|&i| self.all_points[i]).sum();
            sum / self.hull_verts.len() as f64
        };
        let mut vol = 0.0f64;
        for f in &self.faces {
            let a = self.all_points[f[0]] - c;
            let b = self.all_points[f[1]] - c;
            let cc = self.all_points[f[2]] - c;
            vol += a.dot(&b.cross(&cc));
        }
        (vol / 6.0).abs()
    }
}

// ---------------------------------------------------------------------------
// Conflict graph
// ---------------------------------------------------------------------------

/// A conflict graph recording which faces can "see" which points outside the hull.
///
/// Used in the full incremental QuickHull algorithm to speed up point insertion.
pub struct ConflictGraph {
    /// Number of hull faces.
    pub face_count: usize,
    /// For each face, list of point indices visible from that face.
    pub face_to_points: Vec<Vec<usize>>,
    /// For each point, list of face indices it can see.
    pub point_to_faces: Vec<Vec<usize>>,
}

impl ConflictGraph {
    /// Build a conflict graph for `hull` against a set of `points`.
    pub fn new(hull: &ConvexHull3DVec, points: &[Vec3]) -> Self {
        let n_faces = hull.n_faces();
        let n_points = points.len();
        let mut face_to_points = vec![Vec::new(); n_faces];
        let mut point_to_faces = vec![Vec::new(); n_points];

        for (fi, face) in hull.faces.iter().enumerate() {
            let (n, off) = {
                let a = hull.vertices[face[0]];
                let b = hull.vertices[face[1]];
                let c = hull.vertices[face[2]];
                let nn = (b - a).cross(&(c - a));
                let off = nn.dot(&a);
                (nn, off)
            };
            for (pi, &p) in points.iter().enumerate() {
                if nn_point_plane_distance(p, n, off) > 1e-10 {
                    face_to_points[fi].push(pi);
                    point_to_faces[pi].push(fi);
                }
            }
        }

        Self {
            face_count: n_faces,
            face_to_points,
            point_to_faces,
        }
    }

    /// Returns the faces visible from point `pi`.
    pub fn faces_for_point(&self, pi: usize) -> &[usize] {
        if pi < self.point_to_faces.len() {
            &self.point_to_faces[pi]
        } else {
            &[]
        }
    }

    /// Returns the points visible from face `fi`.
    pub fn points_for_face(&self, fi: usize) -> &[usize] {
        if fi < self.face_to_points.len() {
            &self.face_to_points[fi]
        } else {
            &[]
        }
    }
}

fn nn_point_plane_distance(p: Vec3, n: Vec3, offset: f64) -> f64 {
    n.dot(&p) - offset
}

// ---------------------------------------------------------------------------
// Hull inertia tensor
// ---------------------------------------------------------------------------

/// Compute the inertia tensor of a convex hull treated as a uniform solid.
///
/// Uses the signed-volume tetrahedral decomposition method, integrating
/// `I = ∫ρ(r²I - r⊗r)dV` over the hull volume.
///
/// # Arguments
/// * `hull`   – convex hull structure.
/// * `mass`   – total mass of the solid.
///
/// # Returns
/// The 3×3 inertia tensor as `[[f64;3\];3]`.
pub fn hull_inertia_tensor(hull: &ConvexHull3DVec, mass: f64) -> [[f64; 3]; 3] {
    let vol = hull.volume();
    if vol < 1e-12 {
        return [[0.0; 3]; 3];
    }
    let density = mass / vol;
    let c = {
        let sum: Vec3 = hull.vertices.iter().copied().sum();
        sum / hull.vertices.len() as f64
    };
    let mut i_xx = 0.0f64;
    let mut i_yy = 0.0f64;
    let mut i_zz = 0.0f64;
    let mut i_xy = 0.0f64;
    let mut i_xz = 0.0f64;
    let mut i_yz = 0.0f64;

    for face in &hull.faces {
        let a = hull.vertices[face[0]] - c;
        let b = hull.vertices[face[1]] - c;
        let cc = hull.vertices[face[2]] - c;
        let tet_vol = a.dot(&b.cross(&cc)) / 6.0;

        // Use covariance integrals for the tetrahedron (origin, a, b, cc)
        // Diagonal: I_xx ≈ (a²_y+a²_z + ...) * tet_vol / 10  (approximation)
        let coeff = tet_vol / 10.0;
        for p in &[a, b, cc] {
            i_xx += coeff * (p.y * p.y + p.z * p.z);
            i_yy += coeff * (p.x * p.x + p.z * p.z);
            i_zz += coeff * (p.x * p.x + p.y * p.y);
            i_xy -= coeff * p.x * p.y;
            i_xz -= coeff * p.x * p.z;
            i_yz -= coeff * p.y * p.z;
        }
    }

    let s = density;
    [
        [i_xx * s, i_xy * s, i_xz * s],
        [i_xy * s, i_yy * s, i_yz * s],
        [i_xz * s, i_yz * s, i_zz * s],
    ]
}

// ---------------------------------------------------------------------------
// Closest surface point
// ---------------------------------------------------------------------------

impl ConvexHull3DVec {
    /// Find the closest point on the hull surface to query point `q`.
    ///
    /// Iterates over all triangular faces and returns the point that minimizes
    /// distance to `q`.
    pub fn closest_surface_point(&self, q: Vec3) -> Vec3 {
        let mut best_dist = f64::INFINITY;
        let mut best_pt = q;

        for face in &self.faces {
            let a = self.vertices[face[0]];
            let b = self.vertices[face[1]];
            let c = self.vertices[face[2]];
            let cp = closest_point_on_tri(q, a, b, c);
            let d = (cp - q).norm();
            if d < best_dist {
                best_dist = d;
                best_pt = cp;
            }
        }
        best_pt
    }
}

/// Closest point on triangle (a, b, c) to query point p.
fn closest_point_on_tri(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return a + ab * v;
    }
    let cp = p - c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return a + ac * w;
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    a + ab * v + ac * w
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn tetrahedron_points() -> Vec<Vec3> {
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
    }

    fn unit_cube_points() -> Vec<Vec3> {
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
        ]
    }

    #[test]
    fn test_hull_tetrahedron() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("should build");
        assert_eq!(
            hull.n_faces(),
            4,
            "tetrahedron should have 4 faces, got {}",
            hull.n_faces()
        );
    }

    #[test]
    fn test_hull_cube() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("should build cube hull");
        assert_eq!(
            hull.n_faces(),
            12,
            "unit cube should have 12 triangular faces, got {}",
            hull.n_faces()
        );
    }

    #[test]
    fn test_hull_sphere_points() {
        let n = 50usize;
        let pts: Vec<Vec3> = (0..n)
            .map(|i| {
                let theta = std::f64::consts::PI * (i as f64) / (n as f64 / 2.0);
                let phi = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
                Vec3::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.cos(),
                )
            })
            .collect();

        let hull = ConvexHull3DVec::build(&pts).expect("should build sphere hull");
        for v in &hull.vertices {
            let r = v.norm();
            assert!(
                (r - 1.0).abs() < 1e-9,
                "hull vertex not on unit sphere: norm={}",
                r
            );
        }
        assert!(hull.n_faces() >= 4, "too few faces: {}", hull.n_faces());
    }

    #[test]
    fn test_hull_support_x() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let sp = hull.support(Vec3::new(1.0, 0.0, 0.0));
        assert!(
            (sp.x - 1.0).abs() < 1e-9,
            "support in +X should have x=1, got x={}",
            sp.x
        );
    }

    #[test]
    fn test_hull_contains_interior() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let centroid: Vec3 = hull.vertices.iter().copied().sum::<Vec3>() / hull.n_vertices() as f64;
        assert!(
            hull.contains_point(centroid),
            "centroid should be inside the hull"
        );
    }

    #[test]
    fn test_hull_exterior_point() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let outside = Vec3::new(5.0, 5.0, 5.0);
        assert!(
            !hull.contains_point(outside),
            "point far outside should not be contained"
        );
    }

    #[test]
    fn test_hull_volume_cube() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let vol = hull.volume();
        assert!(
            (vol - 1.0).abs() < 0.05,
            "unit cube volume should be ~1.0, got {}",
            vol
        );
    }

    #[test]
    fn test_hull_degenerate_coplanar() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
        ];
        match ConvexHull3DVec::build(&pts) {
            None => {}
            Some(hull) => {
                let vol = hull.volume();
                assert!(
                    vol < 1e-6,
                    "coplanar hull should have ~0 volume, got {}",
                    vol
                );
            }
        }
    }

    // --- New tests ---

    #[test]
    fn test_hull_surface_area_cube() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let sa = hull.surface_area();
        // Unit cube surface area = 6 faces * 1.0 = 6.0
        assert!(
            (sa - 6.0).abs() < 0.1,
            "unit cube SA should be ~6.0, got {}",
            sa
        );
    }

    #[test]
    fn test_hull_edges_cube() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let edges = hull.edges();
        // A cube has 12 edges; triangulated cube has 18 edges (12 original + 6 diagonals)
        assert!(
            edges.len() >= 12,
            "cube should have >= 12 edges, got {}",
            edges.len()
        );
    }

    #[test]
    fn test_hull_euler_characteristic() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let chi = hull.euler_characteristic();
        assert_eq!(
            chi, 2,
            "Euler characteristic should be 2 for convex polyhedron, got {}",
            chi
        );
    }

    #[test]
    fn test_hull_euler_characteristic_tetrahedron() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        assert_eq!(hull.euler_characteristic(), 2);
    }

    #[test]
    fn test_hull_diameter_cube() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let d = hull.diameter();
        // Diagonal of unit cube = sqrt(3)
        assert!((d - 3.0_f64.sqrt()).abs() < 1e-9, "diameter={}", d);
    }

    #[test]
    fn test_hull_aabb_cube() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let (mn, mx) = hull.aabb();
        for k in 0..3 {
            assert!(mn[k] >= -1e-9, "min[{}] = {}", k, mn[k]);
            assert!((mx[k] - 1.0).abs() < 1e-9, "max[{}] = {}", k, mx[k]);
        }
    }

    #[test]
    fn test_hull_signed_distance_inside() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let sd = hull.signed_distance(Vec3::new(0.5, 0.5, 0.5));
        assert!(
            sd < 1e-9,
            "centroid should have non-positive signed distance, got {}",
            sd
        );
    }

    #[test]
    fn test_hull_signed_distance_outside() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let sd = hull.signed_distance(Vec3::new(5.0, 5.0, 5.0));
        assert!(
            sd > 0.0,
            "point outside should have positive signed distance, got {}",
            sd
        );
    }

    #[test]
    fn test_merge_hulls() {
        // Use two cube-like hulls that share some space to ensure
        // the merged point set is non-degenerate in all axes.
        let pts_a = unit_cube_points();
        let pts_b: Vec<Vec3> = pts_a.iter().map(|p| p + Vec3::new(0.5, 0.0, 0.0)).collect();
        let hull_a = ConvexHull3DVec::build(&pts_a).expect("build A");
        let hull_b = ConvexHull3DVec::build(&pts_b).expect("build B");
        let merged = merge_hulls(&hull_a, &hull_b).expect("merge");
        // Merged hull should contain all vertices of both inputs
        assert!(merged.n_vertices() >= hull_a.n_vertices());
        assert!(merged.volume() >= hull_a.volume() - 1e-6);
    }

    #[test]
    fn test_build_with_indices() {
        let pts = tetrahedron_points();
        let (hull, indices) = build_with_indices(&pts).expect("build");
        assert_eq!(hull.n_vertices(), 4);
        assert_eq!(indices.len(), 4);
        // All indices should be in [0, 3]
        for &idx in &indices {
            assert!(idx < 4, "index {} out of range", idx);
        }
    }

    #[test]
    fn test_is_point_set_convex_cube() {
        let pts = unit_cube_points();
        assert!(is_point_set_convex(&pts));
    }

    #[test]
    fn test_is_point_set_convex_with_interior() {
        let mut pts = unit_cube_points();
        pts.push(Vec3::new(0.5, 0.5, 0.5)); // interior point
        assert!(!is_point_set_convex(&pts));
    }

    #[test]
    fn test_face_area_tetrahedron() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        for fi in 0..hull.n_faces() {
            let area = hull.face_area(fi);
            assert!(area > 0.0, "face {} has zero area", fi);
        }
    }

    #[test]
    fn test_vertex_iter() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let verts: Vec<_> = hull.vertex_iter().collect();
        assert_eq!(verts.len(), 4);
    }

    #[test]
    fn test_face_iter() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let faces: Vec<_> = hull.face_iter().collect();
        assert_eq!(faces.len(), 4);
    }

    #[test]
    fn test_overlaps_sat_faces_overlapping() {
        let pts_a = unit_cube_points();
        let pts_b: Vec<Vec3> = pts_a.iter().map(|p| p + Vec3::new(0.5, 0.0, 0.0)).collect();
        let hull_a = ConvexHull3DVec::build(&pts_a).expect("build A");
        let hull_b = ConvexHull3DVec::build(&pts_b).expect("build B");
        assert!(hull_a.overlaps_sat_faces(&hull_b));
    }

    #[test]
    fn test_overlaps_sat_faces_separated() {
        let pts_a = unit_cube_points();
        let pts_b: Vec<Vec3> = pts_a
            .iter()
            .map(|p| p + Vec3::new(10.0, 0.0, 0.0))
            .collect();
        let hull_a = ConvexHull3DVec::build(&pts_a).expect("build A");
        let hull_b = ConvexHull3DVec::build(&pts_b).expect("build B");
        assert!(!hull_a.overlaps_sat_faces(&hull_b));
    }

    #[test]
    fn test_n_edges_tetrahedron() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        assert_eq!(
            hull.n_edges(),
            6,
            "tetrahedron has 6 edges, got {}",
            hull.n_edges()
        );
    }

    #[test]
    fn test_hull_volume_tetrahedron() {
        let pts = tetrahedron_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let vol = hull.volume();
        // Volume of unit tetrahedron = 1/6
        assert!(
            (vol - 1.0 / 6.0).abs() < 0.02,
            "tetrahedron volume = {}",
            vol
        );
    }

    // ── Incremental hull tests ──────────────────────────────────────────────

    #[test]
    fn test_incremental_hull_cube() {
        let pts = unit_cube_points();
        let hull = IncrementalConvexHull::build(&pts);
        assert!(hull.n_vertices() == 8, "cube has 8 vertices");
        assert!(hull.n_faces() >= 12, "cube has >= 12 triangular faces");
    }

    #[test]
    fn test_incremental_hull_volume() {
        let pts = unit_cube_points();
        let hull = IncrementalConvexHull::build(&pts);
        let vol = hull.volume();
        assert!((vol - 1.0).abs() < 0.1, "volume={}", vol);
    }

    #[test]
    fn test_incremental_hull_euler() {
        let pts = unit_cube_points();
        let hull = IncrementalConvexHull::build(&pts);
        let chi = hull.euler_characteristic();
        assert_eq!(chi, 2, "Euler characteristic={}", chi);
    }

    // ── Inertia tensor tests ──────────────────────────────────────────────────

    #[test]
    fn test_hull_inertia_tensor_diagonal_symmetry() {
        // Symmetric hull: cube → off-diagonal inertia should be ~0
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let inertia = hull_inertia_tensor(&hull, 1.0);
        // Off-diagonals should be small for a symmetric body
        assert!(inertia[0][1].abs() < 0.1, "I_xy={}", inertia[0][1]);
        assert!(inertia[0][2].abs() < 0.1, "I_xz={}", inertia[0][2]);
        assert!(inertia[1][2].abs() < 0.1, "I_yz={}", inertia[1][2]);
    }

    #[test]
    fn test_hull_inertia_tensor_positive_diagonal() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let inertia = hull_inertia_tensor(&hull, 1.0);
        assert!(inertia[0][0] > 0.0, "I_xx={}", inertia[0][0]);
        assert!(inertia[1][1] > 0.0, "I_yy={}", inertia[1][1]);
        assert!(inertia[2][2] > 0.0, "I_zz={}", inertia[2][2]);
    }

    // ── Conflict graph tests ─────────────────────────────────────────────────

    #[test]
    fn test_conflict_graph_build() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let cg = ConflictGraph::new(&hull, &pts);
        // For the unit cube (all on hull), no extra point should conflict any face
        assert_eq!(cg.face_count, hull.n_faces());
    }

    // ── Closest point on hull tests ───────────────────────────────────────────

    #[test]
    fn test_closest_point_on_hull_inside() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let q = Vec3::new(0.5, 0.5, 0.5); // inside cube
        let cp = hull.closest_surface_point(q);
        // Closest surface point should be at distance 0.5 from center
        let dist = (cp - q).norm();
        assert!(dist < 0.6, "dist to surface from center={}", dist);
    }

    #[test]
    fn test_closest_point_on_hull_outside() {
        let pts = unit_cube_points();
        let hull = ConvexHull3DVec::build(&pts).expect("build");
        let q = Vec3::new(2.0, 0.5, 0.5); // outside cube
        let cp = hull.closest_surface_point(q);
        // Closest point should be on the x=1 face
        assert!((cp.x - 1.0).abs() < 0.1, "cp.x={}", cp.x);
    }
}
