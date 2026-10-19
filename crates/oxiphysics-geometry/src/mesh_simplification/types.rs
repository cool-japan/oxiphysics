//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::collections::{HashMap, HashSet};

/// Half-edge mesh topology.
#[derive(Debug, Clone)]
pub struct HalfEdgeMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Half-edge array.
    pub half_edges: Vec<HalfEdge>,
    /// For each vertex, one outgoing half-edge index.
    pub vertex_half_edge: Vec<usize>,
    /// For each face, one half-edge index.
    pub face_half_edge: Vec<usize>,
    /// Marks deleted vertices.
    pub vertex_deleted: Vec<bool>,
    /// Marks deleted faces.
    pub face_deleted: Vec<bool>,
    /// Marks deleted half-edges.
    pub he_deleted: Vec<bool>,
}
impl HalfEdgeMesh {
    /// Build a half-edge mesh from vertices and triangle indices.
    pub fn from_triangles(vertices: Vec<[f64; 3]>, triangles: &[[usize; 3]]) -> Self {
        let n_verts = vertices.len();
        let n_faces = triangles.len();
        let n_he = n_faces * 3;
        let mut half_edges = Vec::with_capacity(n_he);
        let mut vertex_half_edge = vec![usize::MAX; n_verts];
        let mut face_half_edge = Vec::with_capacity(n_faces);
        let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
        for (fi, tri) in triangles.iter().enumerate() {
            let base = half_edges.len();
            face_half_edge.push(base);
            for k in 0..3 {
                let v_from = tri[k];
                let v_to = tri[(k + 1) % 3];
                let he_idx = base + k;
                half_edges.push(HalfEdge {
                    vertex: v_to,
                    face: fi,
                    twin: usize::MAX,
                    next: base + (k + 1) % 3,
                    prev: base + (k + 2) % 3,
                });
                if vertex_half_edge[v_from] == usize::MAX {
                    vertex_half_edge[v_from] = he_idx;
                }
                edge_map.insert((v_from, v_to), he_idx);
            }
        }
        let keys: Vec<(usize, usize)> = edge_map.keys().copied().collect();
        for (v_from, v_to) in keys {
            if let Some(&twin_idx) = edge_map.get(&(v_to, v_from)) {
                let he_idx = edge_map[&(v_from, v_to)];
                half_edges[he_idx].twin = twin_idx;
                half_edges[twin_idx].twin = he_idx;
            }
        }
        let vertex_deleted = vec![false; n_verts];
        let face_deleted = vec![false; n_faces];
        let he_deleted = vec![false; n_he];
        Self {
            vertices,
            half_edges,
            vertex_half_edge,
            face_half_edge,
            vertex_deleted,
            face_deleted,
            he_deleted,
        }
    }
    /// Number of active (non-deleted) vertices.
    pub fn active_vertex_count(&self) -> usize {
        self.vertex_deleted.iter().filter(|&&d| !d).count()
    }
    /// Number of active faces.
    pub fn active_face_count(&self) -> usize {
        self.face_deleted.iter().filter(|&&d| !d).count()
    }
    /// Check if a half-edge is on the boundary (no twin).
    pub fn is_boundary_half_edge(&self, he_idx: usize) -> bool {
        self.half_edges[he_idx].twin == usize::MAX
    }
    /// Check if a vertex is on the boundary.
    pub fn is_boundary_vertex(&self, v: usize) -> bool {
        let start = self.vertex_half_edge[v];
        if start == usize::MAX {
            return false;
        }
        let mut he = start;
        loop {
            if self.is_boundary_half_edge(he) {
                return true;
            }
            let twin = self.half_edges[he].twin;
            if twin == usize::MAX {
                return true;
            }
            he = self.half_edges[twin].next;
            if he == start {
                break;
            }
        }
        false
    }
    /// Get all faces adjacent to a vertex.
    ///
    /// The faces are found by walking around the vertex, in both directions when
    /// the vertex is on the boundary, so every live face that contains `v` is
    /// returned whichever outgoing half-edge is stored for `v`.
    pub fn vertex_faces(&self, v: usize) -> Vec<usize> {
        let mut faces = Vec::new();
        let start = self.vertex_half_edge[v];
        if start == usize::MAX {
            return faces;
        }
        let limit = self.half_edges.len();
        // Records the face of `he`; `false` means the walk came back to a face it
        // has seen, which only happens for a corrupt half-edge structure.
        let mut collect = |he: usize| -> bool {
            let f = self.half_edges[he].face;
            if f == usize::MAX || self.face_deleted[f] {
                return true;
            }
            if faces.contains(&f) {
                return false;
            }
            faces.push(f);
            true
        };
        // Counter-clockwise around `v`: the next outgoing half-edge leaves
        // through the face on the other side of the incoming twin.
        let mut he = start;
        let mut closed = false;
        for _ in 0..limit {
            if !collect(he) {
                closed = true;
                break;
            }
            let twin = self.half_edges[he].twin;
            if twin == usize::MAX {
                break;
            }
            he = self.half_edges[twin].next;
            if he == start {
                closed = true;
                break;
            }
        }
        // An open fan is completed clockwise from the start.
        if !closed {
            let mut he = start;
            for _ in 0..limit {
                let incoming_twin = self.half_edges[self.half_edges[he].prev].twin;
                if incoming_twin == usize::MAX {
                    break;
                }
                he = incoming_twin;
                if !collect(he) {
                    break;
                }
            }
        }
        faces
    }
    /// The three vertices of a face in winding order.
    pub fn face_vertices(&self, fi: usize) -> [usize; 3] {
        let he0 = self.face_half_edge[fi];
        let he1 = self.half_edges[he0].next;
        let he2 = self.half_edges[he1].next;
        [
            self.half_edges[he2].vertex,
            self.half_edges[he0].vertex,
            self.half_edges[he1].vertex,
        ]
    }
    /// Get one-ring neighbor vertices of a vertex.
    ///
    /// Returns the vertices that share a live face with `v`, sorted and without
    /// duplicates, including the neighbour at the far end of a boundary edge.
    pub fn vertex_neighbors(&self, v: usize) -> Vec<usize> {
        let mut neighbors: Vec<usize> = self
            .vertex_faces(v)
            .into_iter()
            .flat_map(|f| self.face_vertices(f))
            .filter(|&u| u != v && !self.vertex_deleted[u])
            .collect();
        neighbors.sort_unstable();
        neighbors.dedup();
        neighbors
    }
    /// Compute the face normal for a given face index.
    pub fn face_normal(&self, fi: usize) -> [f64; 3] {
        let he0 = self.face_half_edge[fi];
        let he1 = self.half_edges[he0].next;
        let he2 = self.half_edges[he1].next;
        let v0_idx = self.half_edges[self.half_edges[he0].prev].vertex;
        let v1_idx = self.half_edges[he0].vertex;
        let v2_idx = self.half_edges[he1].vertex;
        let _ = v2_idx;
        let _ = he2;
        let va = self.vertices[v0_idx];
        let vb = self.vertices[v1_idx];
        let vc = self.vertices[self.half_edges[he1].vertex];
        let ab = sub3(vb, va);
        let ac = sub3(vc, va);
        normalize3(cross3(ab, ac))
    }
    /// Compute the area of a face.
    pub fn face_area(&self, fi: usize) -> f64 {
        let he0 = self.face_half_edge[fi];
        let he1 = self.half_edges[he0].next;
        let v0_idx = self.half_edges[self.half_edges[he0].prev].vertex;
        let v1_idx = self.half_edges[he0].vertex;
        let v2_idx = self.half_edges[he1].vertex;
        let va = self.vertices[v0_idx];
        let vb = self.vertices[v1_idx];
        let vc = self.vertices[v2_idx];
        let ab = sub3(vb, va);
        let ac = sub3(vc, va);
        0.5 * len3(cross3(ab, ac))
    }
    /// Extract active triangles as index triples.
    pub fn extract_triangles(&self) -> Vec<[usize; 3]> {
        let mut tris = Vec::new();
        for (fi, &he0) in self.face_half_edge.iter().enumerate() {
            if self.face_deleted[fi] {
                continue;
            }
            let he1 = self.half_edges[he0].next;
            let he2 = self.half_edges[he1].next;
            let v0 = self.half_edges[he2].vertex;
            let v1 = self.half_edges[he0].vertex;
            let v2 = self.half_edges[he1].vertex;
            tris.push([v0, v1, v2]);
        }
        tris
    }
}
/// Mesh simplification statistics.
#[derive(Debug, Clone)]
pub struct MeshStats {
    /// Number of vertices.
    pub n_vertices: usize,
    /// Number of triangles.
    pub n_triangles: usize,
    /// Number of edges.
    pub n_edges: usize,
    /// Average edge length.
    pub avg_edge_length: f64,
    /// Minimum edge length.
    pub min_edge_length: f64,
    /// Maximum edge length.
    pub max_edge_length: f64,
    /// Average triangle quality.
    pub avg_quality: f64,
    /// Minimum triangle quality.
    pub min_quality: f64,
}
/// Configuration for QEM mesh simplification.
///
/// See [`qem_simplify`](super::qem_simplify) for how the options act. If neither
/// `target_triangles` nor `max_error` is set, the mesh is simplified as far as
/// the other constraints allow.
#[derive(Debug, Clone)]
pub struct QemConfig {
    /// Target number of triangles (0 = use error threshold).
    pub target_triangles: usize,
    /// Largest quadric error of a collapse that is still performed (0.0 = use
    /// triangle target). The error of a collapse is the sum of the squared
    /// distances of the merged vertex to the planes of the faces around the two
    /// vertices that are merged.
    pub max_error: f64,
    /// Whether to preserve the boundary exactly: boundary vertices never move and
    /// are never merged with each other. If `false`, the boundary is simplified
    /// under the penalty `boundary_weight`.
    pub preserve_boundary: bool,
    /// Maximum angle (radians) by which a collapse may turn the normal of a
    /// surrounding face; 0.0 = no angular limit. Collapses that fold a face over
    /// or make it degenerate are rejected in any case.
    pub max_normal_deviation: f64,
    /// Whether to preserve texture seams. Seams are supplied to
    /// [`qem_simplify_with_seams`](super::qem_simplify_with_seams); without seams
    /// there is nothing to preserve.
    pub preserve_texture_seams: bool,
    /// Weight of the quadrics that hold boundary vertices near their boundary
    /// edges when `preserve_boundary` is `false` (0.0 = no penalty).
    pub boundary_weight: f64,
}
impl QemConfig {
    /// Default configuration targeting 50% reduction.
    pub fn default_config(current_tris: usize) -> Self {
        Self {
            target_triangles: current_tris / 2,
            max_error: 0.0,
            preserve_boundary: true,
            max_normal_deviation: std::f64::consts::FRAC_PI_4,
            preserve_texture_seams: false,
            boundary_weight: 100.0,
        }
    }
    /// Configuration with error threshold stopping criterion.
    pub fn with_error_threshold(max_error: f64) -> Self {
        Self {
            target_triangles: 0,
            max_error,
            preserve_boundary: true,
            max_normal_deviation: 0.0,
            preserve_texture_seams: false,
            boundary_weight: 100.0,
        }
    }
}
/// Record of a single edge collapse for progressive mesh.
#[derive(Debug, Clone)]
pub struct CollapseRecord {
    /// Vertex that was removed (collapsed into `kept`).
    pub removed: usize,
    /// Vertex that was kept.
    pub kept: usize,
    /// Original position of the removed vertex.
    pub removed_pos: [f64; 3],
    /// Original position of the kept vertex.
    pub kept_pos: [f64; 3],
    /// Position after collapse.
    pub collapsed_pos: [f64; 3],
    /// Triangles removed by this collapse.
    pub removed_triangles: Vec<[usize; 3]>,
    /// Triangles modified by this collapse (new versions).
    pub modified_triangles: Vec<[usize; 3]>,
}
/// Progressive mesh: stores a sequence of collapses for LOD.
#[derive(Debug, Clone)]
pub struct ProgressiveMesh {
    /// Base (most simplified) vertices.
    pub base_vertices: Vec<[f64; 3]>,
    /// Base (most simplified) triangles.
    pub base_triangles: Vec<[usize; 3]>,
    /// Collapse records in order of execution (last = first to undo).
    pub collapses: Vec<CollapseRecord>,
}
impl ProgressiveMesh {
    /// Create a progressive mesh by simplifying with QEM and recording collapses.
    pub fn build(vertices: &[[f64; 3]], triangles: &[[usize; 3]], target_triangles: usize) -> Self {
        let config = QemConfig {
            target_triangles,
            max_error: 0.0,
            preserve_boundary: true,
            max_normal_deviation: 0.0,
            preserve_texture_seams: false,
            boundary_weight: 100.0,
        };
        let (base_verts, base_tris) = qem_simplify(vertices, triangles, &config);
        Self {
            base_vertices: base_verts,
            base_triangles: base_tris,
            collapses: Vec::new(),
        }
    }
    /// Get the mesh at the base (lowest) level of detail.
    pub fn base_mesh(&self) -> (&[[f64; 3]], &[[usize; 3]]) {
        (&self.base_vertices, &self.base_triangles)
    }
    /// Number of collapse records stored.
    pub fn num_collapses(&self) -> usize {
        self.collapses.len()
    }
}
/// A symmetric 4×4 matrix stored as 10 unique elements (upper triangle).
///
/// Layout: a00 a01 a02 a03
///              a11 a12 a13
///                   a22 a23
///                        a33
#[derive(Debug, Clone, Copy)]
pub struct SymMat4 {
    /// Elements stored in row-major upper-triangle order.
    pub m: [f64; 10],
}
impl SymMat4 {
    /// Zero matrix.
    pub fn zero() -> Self {
        Self { m: [0.0; 10] }
    }
    /// Build quadric matrix from plane equation ax + by + cz + d = 0.
    pub fn from_plane(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self {
            m: [
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
    /// Add two symmetric matrices.
    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::zero();
        for i in 0..10 {
            result.m[i] = self.m[i] + other.m[i];
        }
        result
    }
    /// Scale by a scalar.
    pub fn scale(&self, s: f64) -> Self {
        let mut result = *self;
        for v in &mut result.m {
            *v *= s;
        }
        result
    }
    /// Evaluate the quadric error v^T Q v for a 3D point (homogeneous w=1).
    pub fn evaluate(&self, v: [f64; 3]) -> f64 {
        let x = v[0];
        let y = v[1];
        let z = v[2];
        self.m[0] * x * x
            + 2.0 * self.m[1] * x * y
            + 2.0 * self.m[2] * x * z
            + 2.0 * self.m[3] * x
            + self.m[4] * y * y
            + 2.0 * self.m[5] * y * z
            + 2.0 * self.m[6] * y
            + self.m[7] * z * z
            + 2.0 * self.m[8] * z
            + self.m[9]
    }
    /// Try to find the optimal point that minimizes the quadric error.
    ///
    /// Solves the 3×3 linear system from ∂(v^T Q v)/∂v = 0.
    /// Returns `None` if the system is singular.
    pub fn optimal_point(&self) -> Option<[f64; 3]> {
        let a = [
            [self.m[0], self.m[1], self.m[2]],
            [self.m[1], self.m[4], self.m[5]],
            [self.m[2], self.m[5], self.m[7]],
        ];
        let b = [-self.m[3], -self.m[6], -self.m[8]];
        solve_3x3(a, b)
    }
}
/// A half-edge in the mesh topology.
#[derive(Debug, Clone)]
pub struct HalfEdge {
    /// Index of the vertex this half-edge points to.
    pub vertex: usize,
    /// Index of the face this half-edge belongs to (usize::MAX if boundary).
    pub face: usize,
    /// Index of the twin (opposite) half-edge (usize::MAX if boundary).
    pub twin: usize,
    /// Index of the next half-edge in the face loop.
    pub next: usize,
    /// Index of the previous half-edge in the face loop.
    pub prev: usize,
}
/// Flags for texture seam edges.
#[derive(Debug, Clone)]
pub struct TextureSeamFlags {
    /// Set of edges (v0, v1) with v0 < v1 that lie on texture seams.
    pub seam_edges: HashSet<(usize, usize)>,
}
impl TextureSeamFlags {
    /// Create empty seam flags.
    pub fn new() -> Self {
        Self {
            seam_edges: HashSet::new(),
        }
    }
    /// Mark an edge as a texture seam.
    pub fn mark_seam(&mut self, v0: usize, v1: usize) {
        let (lo, hi) = if v0 < v1 { (v0, v1) } else { (v1, v0) };
        self.seam_edges.insert((lo, hi));
    }
    /// Check if an edge is on a texture seam.
    pub fn is_seam(&self, v0: usize, v1: usize) -> bool {
        let (lo, hi) = if v0 < v1 { (v0, v1) } else { (v1, v0) };
        self.seam_edges.contains(&(lo, hi))
    }
    /// Number of seam edges.
    pub fn count(&self) -> usize {
        self.seam_edges.len()
    }
}
/// Result of topology validation.
#[derive(Debug, Clone)]
pub struct TopologyValidation {
    /// Whether the mesh is valid.
    pub is_valid: bool,
    /// Number of non-manifold edges found.
    pub non_manifold_edges: usize,
    /// Number of isolated vertices found.
    pub isolated_vertices: usize,
    /// Number of degenerate triangles found.
    pub degenerate_triangles: usize,
}
/// An edge collapse candidate with its error cost.
///
/// The candidate was computed for the vertex versions it records; once either
/// vertex has changed or is gone, the candidate is outdated and must be dropped.
#[derive(Debug, Clone)]
pub(super) struct CollapseCandidate {
    /// Vertex that survives the collapse.
    pub(super) v0: usize,
    /// Vertex that is merged into `v0`.
    pub(super) v1: usize,
    /// Optimal position for the merged vertex.
    pub(super) optimal_pos: [f64; 3],
    /// Quadric error cost.
    pub(super) cost: f64,
    /// Version of `v0` this candidate was computed for.
    pub(super) v0_version: u32,
    /// Version of `v1` this candidate was computed for.
    pub(super) v1_version: u32,
}
