//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions_2::*;
use std::collections::HashSet;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{Face, Uv, Vertex};

/// An anisotropic size field: target edge lengths along principal curvature
/// directions.
#[derive(Debug, Clone)]
pub struct AnisotropicSizeField {
    /// Per-vertex minimum principal curvature direction.
    pub k_min_dir: Vec<[f64; 3]>,
    /// Per-vertex maximum principal curvature direction.
    pub k_max_dir: Vec<[f64; 3]>,
    /// Target edge length along minimum curvature direction.
    pub len_min: Vec<f64>,
    /// Target edge length along maximum curvature direction.
    pub len_max: Vec<f64>,
}
/// Aggregate quality statistics over a whole mesh.
#[derive(Debug, Clone)]
pub struct MeshQualityStats {
    /// Mean aspect ratio across all faces.
    pub mean_aspect_ratio: f64,
    /// Maximum aspect ratio (worst triangle).
    pub max_aspect_ratio: f64,
    /// Minimum angle (degrees) in the mesh.
    pub min_angle_deg: f64,
    /// Mean area quality score (0..1; 1 = equilateral).
    pub mean_area_quality: f64,
}
/// A feature edge defined by its two endpoint vertex indices and the dihedral
/// angle between the two adjacent faces.
#[derive(Debug, Clone, Copy)]
pub struct FeatureEdge {
    /// First vertex index.
    pub v0: usize,
    /// Second vertex index.
    pub v1: usize,
    /// Dihedral angle between adjacent faces in radians.
    pub dihedral_angle: f64,
}
/// A symmetric 4×4 quadric error matrix stored in compact upper-triangle form.
/// Indices: Q\[0\]=q00, Q\[1\]=q01, Q\[2\]=q02, Q\[3\]=q03, Q\[4\]=q11, Q\[5\]=q12,
/// Q\[6\]=q13, Q\[7\]=q22, Q\[8\]=q23, Q\[9\]=q33.
#[derive(Debug, Clone, Copy, Default)]
pub struct Quadric {
    /// Upper-triangle elements of the 4×4 symmetric matrix.
    pub q: [f64; 10],
}
impl Quadric {
    /// Create a zero quadric.
    pub fn zero() -> Self {
        Self { q: [0.0; 10] }
    }
    /// Add another quadric to this one.
    pub fn add(&self, other: &Quadric) -> Quadric {
        let mut r = Quadric::zero();
        for i in 0..10 {
            r.q[i] = self.q[i] + other.q[i];
        }
        r
    }
    /// Build a fundamental error quadric for the plane (a,b,c,d) where
    /// `ax + by + cz + d = 0` and `(a,b,c)` is the unit normal.
    pub fn from_plane(a: f64, b: f64, c: f64, d: f64) -> Quadric {
        Quadric {
            q: [
                a * a,
                a * b,
                a * c,
                a * d,
                b * b,
                b * c,
                b * d,
                c * c,
                c * d,
                d * d,
            ],
        }
    }
    /// Evaluate Q(v) = vᵀQv (cost of placing a vertex at `v`).
    pub fn evaluate(&self, v: [f64; 3]) -> f64 {
        let [x, y, z] = v;
        let w = 1.0_f64;
        let q = &self.q;
        q[0] * x * x
            + 2.0 * q[1] * x * y
            + 2.0 * q[2] * x * z
            + 2.0 * q[3] * x * w
            + q[4] * y * y
            + 2.0 * q[5] * y * z
            + 2.0 * q[6] * y * w
            + q[7] * z * z
            + 2.0 * q[8] * z * w
            + q[9] * w * w
    }
}
/// Quality metrics for a single triangle.
#[derive(Debug, Clone, Copy)]
pub struct TriangleQuality {
    /// Aspect ratio (circumradius / 2·inradius); perfect equilateral = 1.
    pub aspect_ratio: f64,
    /// Minimum interior angle in radians.
    pub min_angle: f64,
    /// Maximum interior angle in radians.
    pub max_angle: f64,
    /// Normalised area (0..1 relative to equilateral with same perimeter).
    pub area_quality: f64,
}
/// A triangle mesh used as the primary data structure for all processing.
#[derive(Debug, Clone)]
pub struct ProcessMesh {
    /// Vertex positions.
    pub verts: Vec<Vertex>,
    /// Triangular faces (vertex-index triples).
    pub faces: Vec<Face>,
    /// Per-vertex normals (optional).
    pub normals: Option<Vec<Vertex>>,
    /// Per-vertex UV coordinates (optional).
    pub uvs: Option<Vec<Uv>>,
}
impl ProcessMesh {
    /// Create a new mesh from vertices and faces.
    pub fn new(verts: Vec<Vertex>, faces: Vec<Face>) -> Self {
        Self {
            verts,
            faces,
            normals: None,
            uvs: None,
        }
    }
    /// Return the number of vertices.
    pub fn num_verts(&self) -> usize {
        self.verts.len()
    }
    /// Return the number of faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }
    /// Compute per-vertex normals by averaging adjacent face normals.
    pub fn compute_normals(&mut self) {
        let nv = self.verts.len();
        let mut normals = vec![[0.0_f64; 3]; nv];
        for &[a, b, c] in &self.faces {
            let va = self.verts[a];
            let vb = self.verts[b];
            let vc = self.verts[c];
            let ab = vec3_sub(vb, va);
            let ac = vec3_sub(vc, va);
            let n = vec3_cross(ab, ac);
            for &i in &[a, b, c] {
                normals[i] = vec3_add(normals[i], n);
            }
        }
        for n in &mut normals {
            *n = vec3_normalize(*n);
        }
        self.normals = Some(normals);
    }
    /// Build a map from each vertex to its list of neighbouring vertex indices.
    pub fn build_adjacency(&self) -> Vec<Vec<usize>> {
        let nv = self.verts.len();
        let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); nv];
        for &[a, b, c] in &self.faces {
            adj[a].insert(b);
            adj[a].insert(c);
            adj[b].insert(a);
            adj[b].insert(c);
            adj[c].insert(a);
            adj[c].insert(b);
        }
        adj.into_iter().map(|s| s.into_iter().collect()).collect()
    }
}
/// UV parameterization result with per-vertex texture coordinates.
#[derive(Debug, Clone)]
pub struct UvParameterization {
    /// Mesh after parameterization (same topology).
    pub mesh: ProcessMesh,
    /// Per-vertex UV coordinates in \[0,1\]².
    pub uvs: Vec<Uv>,
}
/// A rectangular texture patch cut from the atlas.
#[derive(Debug, Clone)]
pub struct AtlasPatch {
    /// Index of the original face group.
    pub group_id: usize,
    /// UV bounding box: \[umin, vmin, umax, vmax\].
    pub bounds: [f64; 4],
    /// Vertices (UVs) of the patch in atlas space.
    pub uvs: Vec<Uv>,
}
/// Result of a boolean mesh operation.
#[derive(Debug, Clone)]
pub struct BooleanResult {
    /// The combined / intersected / subtracted mesh.
    pub mesh: ProcessMesh,
    /// Whether the operation was exact (currently always false for this stub).
    pub is_exact: bool,
}
