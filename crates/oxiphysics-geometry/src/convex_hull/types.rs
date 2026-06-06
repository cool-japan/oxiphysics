//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::{Real, Vec3};

use super::functions::*;

/// A convex hull defined by a set of vertices (uses nalgebra `Vec3`).
#[derive(Debug, Clone)]
pub struct ConvexHull {
    /// Vertices of the convex hull.
    pub vertices: Vec<Vec3>,
}
impl ConvexHull {
    /// Create a new convex hull from vertices.
    pub fn new(vertices: Vec<Vec3>) -> Self {
        Self { vertices }
    }
    /// Compute volume and center of mass via tetrahedral decomposition.
    /// Returns (volume, center_of_mass).
    pub(super) fn tetra_decomposition(&self) -> (Real, Vec3) {
        if self.vertices.len() < 4 {
            return (0.0, Vec3::zeros());
        }
        let reference = self.centroid();
        let n = self.vertices.len();
        let mut total_volume = 0.0;
        let mut weighted_center = Vec3::zeros();
        for i in 1..n {
            for j in (i + 1)..n {
                let a = self.vertices[0] - reference;
                let b = self.vertices[i] - reference;
                let c = self.vertices[j] - reference;
                let vol = a.dot(&b.cross(&c)) / 6.0;
                let tet_center =
                    (reference + self.vertices[0] + self.vertices[i] + self.vertices[j]) * 0.25;
                total_volume += vol;
                weighted_center += tet_center * vol;
            }
        }
        if total_volume.abs() > 1e-12 {
            weighted_center /= total_volume;
            (total_volume.abs(), weighted_center)
        } else {
            (0.0, self.centroid())
        }
    }
    fn centroid(&self) -> Vec3 {
        if self.vertices.is_empty() {
            return Vec3::zeros();
        }
        let sum: Vec3 = self.vertices.iter().sum();
        sum / self.vertices.len() as Real
    }
    /// Find the support point (vertex with maximum dot product in the given direction).
    pub fn support_point(&self, direction: &Vec3) -> Vec3 {
        let mut best = Vec3::zeros();
        let mut best_dot = f64::NEG_INFINITY;
        for v in &self.vertices {
            let d = v.dot(direction);
            if d > best_dot {
                best_dot = d;
                best = *v;
            }
        }
        best
    }
    /// Compute the axis-aligned bounding box.
    pub fn bounding_box(&self) -> ConvexHullAabb {
        let mut min = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for v in &self.vertices {
            for k in 0..3 {
                if v[k] < min[k] {
                    min[k] = v[k];
                }
                if v[k] > max[k] {
                    max[k] = v[k];
                }
            }
        }
        ConvexHullAabb { min, max }
    }
}
/// Axis-aligned bounding box for ConvexHull.
#[derive(Debug, Clone)]
pub struct ConvexHullAabb {
    /// Minimum corner.
    pub min: Vec3,
    /// Maximum corner.
    pub max: Vec3,
}
/// 3-D convex hull computed via the Quickhull algorithm.
///
/// Uses plain `[f64; 3]` arrays — no nalgebra dependency.
#[derive(Debug, Clone)]
pub struct ConvexHull3D {
    /// Unique vertices of the hull (subset of the input point cloud).
    pub vertices: Vec<[f64; 3]>,
    /// Triangular faces as triples of indices into `vertices` (CCW outward winding).
    pub faces: Vec<[usize; 3]>,
}
impl ConvexHull3D {
    /// Build a convex hull from a point cloud using the Quickhull algorithm.
    ///
    /// Returns `None` if the input has fewer than 4 non-coplanar points.
    pub fn build(points: &[[f64; 3]]) -> Option<ConvexHull3D> {
        if points.len() < 4 {
            return None;
        }
        let (i0, i1) = ch_extreme_pair(points)?;
        let i2 = ch_farthest_from_line(points, i0, i1)?;
        let n_seed = face_normal(points[i0], points[i1], points[i2]);
        if dot(n_seed, n_seed) < 1e-24 {
            return None;
        }
        let i3 = ch_farthest_from_plane(points, n_seed, points[i0], &[i0, i1, i2])?;
        if face_plane_dist(n_seed, points[i0], points[i3]).abs() < 1e-12 {
            return None;
        }
        let tet = [i0, i1, i2, i3];
        let mut hull: Vec<[usize; 3]> = Vec::new();
        let combos: [([usize; 3], usize); 4] = [
            ([i0, i1, i2], i3),
            ([i0, i1, i3], i2),
            ([i0, i2, i3], i1),
            ([i1, i2, i3], i0),
        ];
        for (face, interior) in &combos {
            if let Some(oriented) = ch_orient_face(*face, *interior, points) {
                hull.push(oriented);
            }
        }
        let exterior: Vec<usize> = (0..points.len()).filter(|i| !tet.contains(i)).collect();
        for &p_idx in &exterior {
            let p = points[p_idx];
            let visible: Vec<usize> = hull
                .iter()
                .enumerate()
                .filter(|(_, face)| {
                    let n = face_normal(points[face[0]], points[face[1]], points[face[2]]);
                    face_plane_dist(n, points[face[0]], p) > 1e-10
                })
                .map(|(i, _)| i)
                .collect();
            if visible.is_empty() {
                continue;
            }
            let horizon = ch_horizon_edges(&hull, &visible);
            let mut new_hull: Vec<[usize; 3]> = hull
                .iter()
                .enumerate()
                .filter(|(i, _)| !visible.contains(i))
                .map(|(_, f)| *f)
                .collect();
            for (e0, e1) in horizon {
                let n_new = face_normal(points[e0], points[e1], points[p_idx]);
                if dot(n_new, n_new) < 1e-24 {
                    continue;
                }
                let interior_ref = new_hull.first().and_then(|f| {
                    f.iter()
                        .find(|&&v| v != e0 && v != e1 && v != p_idx)
                        .copied()
                });
                let oriented = match interior_ref {
                    Some(ref_v) => {
                        if face_plane_dist(n_new, points[e0], points[ref_v]) > 0.0 {
                            [e1, e0, p_idx]
                        } else {
                            [e0, e1, p_idx]
                        }
                    }
                    None => [e0, e1, p_idx],
                };
                new_hull.push(oriented);
            }
            hull = new_hull;
        }
        if hull.is_empty() {
            return None;
        }
        let mut used: Vec<usize> = hull.iter().flat_map(|f| f.iter().copied()).collect();
        used.sort_unstable();
        used.dedup();
        let vertices: Vec<[f64; 3]> = used.iter().map(|&i| points[i]).collect();
        let remap = |old: usize| -> usize {
            used.binary_search(&old)
                .expect("element must be present in sorted list")
        };
        let faces: Vec<[usize; 3]> = hull
            .iter()
            .map(|f| [remap(f[0]), remap(f[1]), remap(f[2])])
            .collect();
        Some(ConvexHull3D { vertices, faces })
    }
    /// Volume via signed tetrahedral decomposition from the centroid.
    pub fn volume(&self) -> f64 {
        if self.faces.is_empty() || self.vertices.is_empty() {
            return 0.0;
        }
        let c = self.centroid();
        let mut vol = 0.0f64;
        for face in &self.faces {
            let a = sub(self.vertices[face[0]], c);
            let b = sub(self.vertices[face[1]], c);
            let cc = sub(self.vertices[face[2]], c);
            vol += dot(a, cross(b, cc));
        }
        (vol / 6.0).abs()
    }
    /// Surface area as the sum of triangle areas.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0f64;
        for face in &self.faces {
            let n = face_normal(
                self.vertices[face[0]],
                self.vertices[face[1]],
                self.vertices[face[2]],
            );
            area += dot(n, n).sqrt() * 0.5;
        }
        area
    }
    /// Returns `true` if `p` is inside (or on the boundary of) the convex hull.
    ///
    /// A point is inside iff it is on the non-positive side of every face plane.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        for face in &self.faces {
            let n = face_normal(
                self.vertices[face[0]],
                self.vertices[face[1]],
                self.vertices[face[2]],
            );
            if face_plane_dist(n, self.vertices[face[0]], p) > 1e-9 {
                return false;
            }
        }
        true
    }
    /// Support function: returns the furthest vertex in direction `dir` (for GJK).
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        self.vertices
            .iter()
            .copied()
            .max_by(|a, b| {
                dot(*a, dir)
                    .partial_cmp(&dot(*b, dir))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or([0.0, 0.0, 0.0])
    }
    /// Number of hull vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Number of triangular faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }
    fn centroid(&self) -> [f64; 3] {
        let n = self.vertices.len() as f64;
        let mut s = [0.0f64; 3];
        for v in &self.vertices {
            s = add(s, *v);
        }
        scale(s, 1.0 / n)
    }
}
impl ConvexHull3D {
    /// Compute the volume of the convex hull via the divergence theorem.
    ///
    /// Each triangular face contributes a signed tetrahedral volume from the
    /// origin.  The absolute value gives the total enclosed volume.
    pub fn compute_volume(&self) -> f64 {
        if self.faces.is_empty() {
            return 0.0;
        }
        let mut vol = 0.0f64;
        for face in &self.faces {
            let a = self.vertices[face[0]];
            let b = self.vertices[face[1]];
            let c = self.vertices[face[2]];
            vol += dot(a, cross(b, c));
        }
        (vol / 6.0).abs()
    }
    /// Compute the center of mass (centroid) of the solid convex hull.
    ///
    /// Uses the signed-tetrahedral decomposition: each triangle + origin
    /// forms a tetrahedron whose centroid is `(a + b + c + origin) / 4`.
    /// The weighted average of tet centroids (by signed volume) is the COM.
    ///
    /// Returns `[0.0; 3]` for empty or degenerate hulls.
    pub fn compute_center_of_mass(&self) -> [f64; 3] {
        if self.faces.is_empty() {
            return [0.0; 3];
        }
        let mut total_vol = 0.0f64;
        let mut weighted = [0.0f64; 3];
        let origin = [0.0f64; 3];
        for face in &self.faces {
            let a = self.vertices[face[0]];
            let b = self.vertices[face[1]];
            let c = self.vertices[face[2]];
            let signed_vol = dot(a, cross(b, c)) / 6.0;
            let tet_com = scale(add(add(add(origin, a), b), c), 0.25);
            total_vol += signed_vol;
            weighted = add(weighted, scale(tet_com, signed_vol));
        }
        if total_vol.abs() > 1e-24 {
            scale(weighted, 1.0 / total_vol)
        } else {
            self.centroid()
        }
    }
    /// Compute the 3×3 inertia tensor of the solid convex hull, assuming
    /// uniform density with total `mass`.
    ///
    /// The tensor is returned in row-major order as `[[f64; 3\]; 3]`.
    /// Off-diagonal terms (products of inertia) are included.
    ///
    /// Uses the tetrahedral decomposition method from Polyhedral Mass Properties
    /// (David Eberly, Geometric Tools).
    pub fn compute_moment_of_inertia(&self, mass: f64) -> [[f64; 3]; 3] {
        if self.faces.is_empty() || mass <= 0.0 {
            return [[0.0; 3]; 3];
        }
        let com = self.compute_center_of_mass();
        let mut ixx = 0.0f64;
        let mut iyy = 0.0f64;
        let mut izz = 0.0f64;
        let mut ixy = 0.0f64;
        let mut ixz = 0.0f64;
        let mut iyz = 0.0f64;
        let mut total_vol = 0.0f64;
        for face in &self.faces {
            let a = sub(self.vertices[face[0]], com);
            let b = sub(self.vertices[face[1]], com);
            let c = sub(self.vertices[face[2]], com);
            let signed_vol = dot(a, cross(b, c)) / 6.0;
            total_vol += signed_vol;
            for &(p, q) in &[(a, b), (b, c), (a, c), (a, a), (b, b), (c, c)] {
                let _ = (p, q);
            }
            let ax = a[0];
            let ay = a[1];
            let az = a[2];
            let bx = b[0];
            let by = b[1];
            let bz = b[2];
            let cx = c[0];
            let cy = c[1];
            let cz = c[2];
            let f2x = ax * ax + bx * bx + cx * cx + ax * bx + ax * cx + bx * cx;
            let f2y = ay * ay + by * by + cy * cy + ay * by + ay * cy + by * cy;
            let f2z = az * az + bz * bz + cz * cz + az * bz + az * cz + bz * cz;
            let fxy = ax * ay
                + bx * by
                + cx * cy
                + 0.5 * (ax * by + ay * bx + ax * cy + ay * cx + bx * cy + by * cx);
            let fxz = ax * az
                + bx * bz
                + cx * cz
                + 0.5 * (ax * bz + az * bx + ax * cz + az * cx + bx * cz + bz * cx);
            let fyz = ay * az
                + by * bz
                + cy * cz
                + 0.5 * (ay * bz + az * by + ay * cz + az * cy + by * cz + bz * cy);
            let w = signed_vol / 10.0;
            ixx += w * (f2y + f2z);
            iyy += w * (f2x + f2z);
            izz += w * (f2x + f2y);
            ixy -= w * fxy;
            ixz -= w * fxz;
            iyz -= w * fyz;
        }
        let volume = total_vol.abs();
        if volume < 1e-24 {
            return [[0.0; 3]; 3];
        }
        let density = mass / volume;
        [
            [density * ixx, density * ixy, density * ixz],
            [density * ixy, density * iyy, density * iyz],
            [density * ixz, density * iyz, density * izz],
        ]
    }
}
impl ConvexHull3D {
    /// Compute the centroid (arithmetic mean) of hull vertices.
    pub fn vertex_centroid(&self) -> [f64; 3] {
        self.centroid()
    }
    /// Compute a signed-tetrahedral volume estimate and centroid together.
    ///
    /// More accurate than the vertex centroid for volume purposes.
    /// Returns `(volume, centroid)`.
    pub fn volume_centroid(&self) -> (f64, [f64; 3]) {
        if self.faces.is_empty() {
            return (0.0, self.centroid());
        }
        let c = self.centroid();
        let mut total_vol = 0.0f64;
        let mut weighted = [0.0f64; 3];
        for face in &self.faces {
            let a = sub(self.vertices[face[0]], c);
            let b = sub(self.vertices[face[1]], c);
            let cc = sub(self.vertices[face[2]], c);
            let signed_vol = dot(a, cross(b, cc)) / 6.0;
            let tet_c = scale(
                add(
                    add(add(c, self.vertices[face[0]]), self.vertices[face[1]]),
                    self.vertices[face[2]],
                ),
                0.25,
            );
            total_vol += signed_vol;
            weighted = add(weighted, scale(tet_c, signed_vol));
        }
        if total_vol.abs() > 1e-24 {
            (total_vol.abs(), scale(weighted, 1.0 / total_vol))
        } else {
            (0.0, c)
        }
    }
}
/// Build an approximate convex hull after optionally deduplicating and
/// thinning the input, and return it together with quality metrics.
pub struct HullWithQuality {
    /// The convex hull.
    pub hull: ConvexHull3D,
    /// Sphericity ∈ (0, 1].
    pub sphericity: f64,
    /// Aspect ratio (1 = cube-like, ∞ = degenerate).
    pub aspect_ratio: f64,
}
impl HullWithQuality {
    /// Build the hull and compute quality metrics in one step.
    pub fn build(points: &[[f64; 3]], tolerance: f64) -> Option<Self> {
        let unique = deduplicate_points(points, tolerance);
        let hull = ConvexHull3D::build(&unique)?;
        let sphericity = hull_sphericity(&hull);
        let aspect_ratio = hull_aspect_ratio(&hull);
        Some(HullWithQuality {
            hull,
            sphericity,
            aspect_ratio,
        })
    }
}
/// Incremental 3D convex hull builder.
///
/// Allows adding points one at a time and maintains a valid convex hull at all
/// times.  Based on the Quickhull increment approach: start with a seed
/// tetrahedron and add points using the existing `ConvexHull3D::build` logic
/// incrementally by maintaining a mutable point cloud and rebuilding.
///
/// For large datasets, prefer `ConvexHull3D::build` directly.  This builder is
/// useful when points arrive in a stream.
pub struct IncrementalConvexHull {
    /// All points added so far.
    pub(super) points: Vec<[f64; 3]>,
    /// Current best hull (if at least 4 non-coplanar points exist).
    pub(super) hull: Option<ConvexHull3D>,
}
impl IncrementalConvexHull {
    /// Create an empty incremental hull builder.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            hull: None,
        }
    }
    /// Add a single point and update the convex hull.
    ///
    /// If the point is already inside the current hull, the hull is unchanged.
    /// Otherwise, the hull is rebuilt from all points.
    pub fn add_point(&mut self, p: [f64; 3]) {
        if let Some(h) = &self.hull
            && h.contains_point(p)
        {
            return;
        }
        self.points.push(p);
        self.hull = ConvexHull3D::build(&self.points);
    }
    /// Add multiple points at once.
    pub fn add_points(&mut self, pts: &[[f64; 3]]) {
        for &p in pts {
            self.add_point(p);
        }
    }
    /// Get a reference to the current hull.
    pub fn hull(&self) -> Option<&ConvexHull3D> {
        self.hull.as_ref()
    }
    /// Number of points added so far.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
    /// Remove all points and reset.
    pub fn clear(&mut self) {
        self.points.clear();
        self.hull = None;
    }
    /// Compute the current volume, or 0 if not yet a valid hull.
    pub fn volume(&self) -> f64 {
        self.hull.as_ref().map(|h| h.volume()).unwrap_or(0.0)
    }
}
