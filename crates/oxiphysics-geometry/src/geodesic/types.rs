//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// A single geodesic isoline segment on the mesh.
#[derive(Debug, Clone)]
pub struct IsolineSegment {
    /// Start point of the segment.
    pub start: [f64; 3],
    /// End point of the segment.
    pub end: [f64; 3],
}
/// A geodesic Voronoi cell: the set of mesh vertices closest to a given source.
#[derive(Debug, Clone)]
pub struct GeoVoronoiCell {
    /// Index into the sources array.
    pub source_idx: usize,
    /// Original source vertex index.
    pub source_vertex: usize,
    /// Vertex indices assigned to this cell.
    pub vertices: Vec<usize>,
    /// Total surface area of this cell (sum of proportional triangle areas).
    pub area: f64,
}
/// Result of intrinsic Delaunay triangulation.
#[derive(Debug, Clone)]
pub struct IntrinsicDelaunay {
    /// Edge lengths for each face edge: `edge_lengths[fi][k]` is the length of
    /// edge opposite to local vertex `k` in face `fi`.
    pub edge_lengths: Vec<[f64; 3]>,
    /// Face connectivity (same as input but potentially with flipped edges).
    pub faces: Vec<[usize; 3]>,
    /// Vertex positions (unchanged from input).
    pub vertices: Vec<[f64; 3]>,
}
/// A simple triangle mesh for geodesic computation.
#[derive(Debug, Clone)]
pub struct GeoMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle face indices (each entry is \[v0, v1, v2\]).
    pub faces: Vec<[usize; 3]>,
}
impl GeoMesh {
    /// Create a new geodesic mesh.
    pub fn new(vertices: Vec<[f64; 3]>, faces: Vec<[usize; 3]>) -> Self {
        Self { vertices, faces }
    }
    /// Number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }
    /// Build an adjacency list: for each vertex, list all neighbouring vertex
    /// indices along with the edge length.
    pub fn build_adjacency(&self) -> Vec<Vec<(usize, f64)>> {
        let n = self.vertices.len();
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for face in &self.faces {
            let [a, b, c] = *face;
            let lab = dist3(self.vertices[a], self.vertices[b]);
            let lbc = dist3(self.vertices[b], self.vertices[c]);
            let lca = dist3(self.vertices[c], self.vertices[a]);
            adj[a].push((b, lab));
            adj[b].push((a, lab));
            adj[b].push((c, lbc));
            adj[c].push((b, lbc));
            adj[c].push((a, lca));
            adj[a].push((c, lca));
        }
        adj
    }
    /// Compute the area of a triangular face.
    pub fn face_area(&self, face_idx: usize) -> f64 {
        let [a, b, c] = self.faces[face_idx];
        let va = self.vertices[a];
        let vb = self.vertices[b];
        let vc = self.vertices[c];
        let ab = sub3(vb, va);
        let ac = sub3(vc, va);
        len3(cross3(ab, ac)) * 0.5
    }
    /// Compute the normal of a triangular face (not normalised).
    pub fn face_normal_raw(&self, face_idx: usize) -> [f64; 3] {
        let [a, b, c] = self.faces[face_idx];
        let ab = sub3(self.vertices[b], self.vertices[a]);
        let ac = sub3(self.vertices[c], self.vertices[a]);
        cross3(ab, ac)
    }
}
/// Result of the heat method geodesic computation.
pub struct HeatGeodesicResult {
    /// Geodesic distances from the source vertex to all vertices.
    pub distances: Vec<f64>,
    /// Normalised heat gradient (per face, not per vertex) for debugging.
    pub gradient: Vec<[f64; 3]>,
}
/// Entry for the Dijkstra priority queue.
#[derive(PartialEq)]
pub(super) struct DijkEntry {
    pub(super) dist: f64,
    pub(super) vertex: usize,
}
