//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet};

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    fill_hole_ear_clipping, find_boundary_loops, fix_non_manifold_edges, merge_duplicate_vertices,
};

/// A single directed half-edge in a half-edge mesh.
pub struct HalfEdge {
    /// Index of the vertex this half-edge points *to*.
    pub vertex: usize,
    /// Index of the face this half-edge belongs to.
    pub face: usize,
    /// Index of the next half-edge in the same face (CCW).
    pub next: usize,
    /// Index of the previous half-edge in the same face.
    pub prev: usize,
    /// Index of the opposite (twin) half-edge, if any.
    pub twin: Option<usize>,
}
/// Half-edge mesh data structure for topology queries.
pub struct HalfEdgeMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// All half-edges in the mesh.
    pub half_edges: Vec<HalfEdge>,
    /// For each face, the index of its first half-edge.
    pub face_starts: Vec<usize>,
}
impl HalfEdgeMesh {
    /// Build a half-edge mesh from a triangle soup.
    pub fn from_triangle_soup(vertices: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>) -> Self {
        let mut half_edges: Vec<HalfEdge> = Vec::new();
        let mut face_starts: Vec<usize> = Vec::new();
        let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
        for (face_idx, tri) in triangles.iter().enumerate() {
            let base = half_edges.len();
            face_starts.push(base);
            for k in 0..3 {
                let to_v = tri[(k + 1) % 3];
                half_edges.push(HalfEdge {
                    vertex: to_v,
                    face: face_idx,
                    next: base + (k + 1) % 3,
                    prev: base + (k + 2) % 3,
                    twin: None,
                });
                edge_map.insert((tri[k], to_v), base + k);
            }
        }
        let n = half_edges.len();
        for i in 0..n {
            if half_edges[i].twin.is_none() {
                let prev_idx = half_edges[i].prev;
                let from_v = half_edges[prev_idx].vertex;
                let to_v = half_edges[i].vertex;
                if let Some(&twin_idx) = edge_map.get(&(to_v, from_v)) {
                    half_edges[i].twin = Some(twin_idx);
                    half_edges[twin_idx].twin = Some(i);
                }
            }
        }
        HalfEdgeMesh {
            vertices,
            half_edges,
            face_starts,
        }
    }
    /// Returns `true` if the mesh is manifold.
    pub fn is_manifold(&self) -> bool {
        let mut seen: HashMap<(usize, usize), usize> = HashMap::new();
        for (i, he) in self.half_edges.iter().enumerate() {
            let prev = &self.half_edges[he.prev];
            let from_v = prev.vertex;
            let to_v = he.vertex;
            let count = seen.entry((from_v, to_v)).or_insert(0);
            *count += 1;
            if *count > 1 {
                let _ = i;
                return false;
            }
        }
        true
    }
    /// Returns the indices of all boundary half-edges (those with no twin).
    pub fn boundary_edges(&self) -> Vec<usize> {
        self.half_edges
            .iter()
            .enumerate()
            .filter_map(|(i, he)| if he.twin.is_none() { Some(i) } else { None })
            .collect()
    }
    /// Returns the valence (number of incident edges) of vertex `v`.
    pub fn vertex_valence(&self, v: usize) -> usize {
        self.half_edges.iter().filter(|he| he.vertex == v).count()
    }
    /// Returns the number of faces in the mesh.
    pub fn face_count(&self) -> usize {
        self.face_starts.len()
    }
    /// Returns the number of vertices in the mesh.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    /// Collect boundary loops: returns a list of loops, each loop is a list
    /// of vertex indices forming a closed boundary.
    pub fn boundary_loops(&self) -> Vec<Vec<usize>> {
        let boundary_hes = self.boundary_edges();
        if boundary_hes.is_empty() {
            return Vec::new();
        }
        let mut from_to_he: HashMap<usize, usize> = HashMap::new();
        for &he_idx in &boundary_hes {
            let prev_idx = self.half_edges[he_idx].prev;
            let from_v = self.half_edges[prev_idx].vertex;
            from_to_he.insert(from_v, he_idx);
        }
        let mut visited: HashSet<usize> = HashSet::new();
        let mut loops = Vec::new();
        for &start_he in &boundary_hes {
            if visited.contains(&start_he) {
                continue;
            }
            let mut loop_verts = Vec::new();
            let mut current = start_he;
            loop {
                visited.insert(current);
                let to_v = self.half_edges[current].vertex;
                loop_verts.push(to_v);
                match from_to_he.get(&to_v) {
                    Some(&next_he) if !visited.contains(&next_he) => {
                        current = next_he;
                    }
                    _ => break,
                }
            }
            if loop_verts.len() >= 3 {
                loops.push(loop_verts);
            }
        }
        loops
    }
}
/// Summary report from a mesh repair pass.
pub struct MeshRepairReport {
    /// Number of degenerate triangles removed.
    pub n_degenerate_removed: usize,
    /// Number of holes filled.
    pub n_holes_filled: usize,
    /// Number of duplicate vertices merged.
    pub n_duplicates_merged: usize,
    /// Number of face normals flipped for consistency.
    pub n_normals_fixed: usize,
}
/// A simple triangle mesh with an adjacency structure that supports edge flips.
///
/// Useful for converting a mesh towards a Delaunay triangulation by flipping
/// edges whose opposite angles violate the Delaunay in-circle criterion.
pub struct FlipMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle list; each triangle stores 3 vertex indices.
    pub triangles: Vec<[usize; 3]>,
    /// Number of edge-flips performed in the last call to `flip_to_delaunay`.
    pub n_flips: usize,
}
impl FlipMesh {
    /// Create a new FlipMesh from existing vertices and triangles.
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>) -> Self {
        Self {
            vertices,
            triangles,
            n_flips: 0,
        }
    }
    /// Check if the shared edge between two triangles satisfies the Delaunay
    /// in-circle criterion (i.e., neither opposite vertex lies inside the other
    /// triangle's circumcircle). Returns `true` if the edge is locally Delaunay.
    ///
    /// The two triangles share vertices `a` and `b`; the opposite vertices are
    /// `c` (triangle 0) and `d` (triangle 1).
    pub fn is_locally_delaunay(va: [f64; 3], vb: [f64; 3], vc: [f64; 3], vd: [f64; 3]) -> bool {
        let ax = va[0] - vd[0];
        let ay = va[1] - vd[1];
        let bx = vb[0] - vd[0];
        let by = vb[1] - vd[1];
        let cx = vc[0] - vd[0];
        let cy = vc[1] - vd[1];
        let det = ax * (by * (cx * cx + cy * cy) - cy * (bx * bx + by * by))
            - ay * (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by))
            + (ax * ax + ay * ay) * (bx * cy - by * cx);
        det <= 0.0
    }
    /// Attempt one pass of edge-flipping towards Delaunay. Returns the number of
    /// edges flipped. Repeat until 0 to fully Delaunay-ize a planar mesh.
    pub fn flip_pass(&mut self) -> usize {
        let n = self.triangles.len();
        let mut edge_map: std::collections::HashMap<(usize, usize), Vec<usize>> =
            std::collections::HashMap::new();
        for (ti, tri) in self.triangles.iter().enumerate() {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                edge_map.entry(key).or_default().push(ti);
            }
        }
        let mut flips = 0usize;
        let candidates: Vec<((usize, usize), usize, usize)> = edge_map
            .iter()
            .filter_map(|(&key, tris)| {
                if tris.len() == 2 {
                    Some((key, tris[0], tris[1]))
                } else {
                    None
                }
            })
            .collect();
        for (edge_key, t0, t1) in candidates {
            if t0 >= n || t1 >= n {
                continue;
            }
            let tri0 = self.triangles[t0];
            let tri1 = self.triangles[t1];
            let (ea, eb) = edge_key;
            let opp0 = tri0.iter().copied().find(|&v| v != ea && v != eb);
            let opp1 = tri1.iter().copied().find(|&v| v != ea && v != eb);
            let (c, d) = match (opp0, opp1) {
                (Some(c), Some(d)) => (c, d),
                _ => continue,
            };
            let va = self.vertices[ea];
            let vb = self.vertices[eb];
            let vc = self.vertices[c];
            let vd = self.vertices[d];
            if !Self::is_locally_delaunay(va, vb, vc, vd) {
                self.triangles[t0] = [ea, c, d];
                self.triangles[t1] = [eb, d, c];
                flips += 1;
            }
        }
        self.n_flips += flips;
        flips
    }
    /// Run up to `max_passes` flip passes until no more flips occur.
    pub fn flip_to_delaunay(&mut self, max_passes: usize) {
        for _ in 0..max_passes {
            if self.flip_pass() == 0 {
                break;
            }
        }
    }
}
/// Compute mesh statistics: vertex count, triangle count, edge count,
/// boundary edge count, non-manifold edge count, genus estimate.
pub struct MeshStats {
    /// Number of vertices.
    pub n_vertices: usize,
    /// Number of triangles.
    pub n_triangles: usize,
    /// Number of unique edges.
    pub n_edges: usize,
    /// Number of boundary edges.
    pub n_boundary_edges: usize,
    /// Number of non-manifold edges.
    pub n_non_manifold_edges: usize,
    /// Euler characteristic (V - E + F).
    pub euler_characteristic: i64,
}
/// A stateful mesh repair helper that owns the vertex and triangle lists.
///
/// Provides higher-level repair operations beyond the standalone functions:
/// - `fill_holes` – fan-triangulate every boundary loop
/// - `remove_duplicate_vertices` – merge vertices within epsilon
/// - `fix_non_manifold_edges` – remove excess triangles until manifold
/// - `compute_euler_characteristic` – V − E + F
pub struct MeshRepair {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle index triples.
    pub triangles: Vec<[usize; 3]>,
}
impl MeshRepair {
    /// Create a new `MeshRepair` from existing geometry.
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>) -> Self {
        Self {
            vertices,
            triangles,
        }
    }
    /// Fill all open boundary holes using fan triangulation.
    ///
    /// For each boundary loop a "fan" of triangles is built from the first
    /// vertex of the loop to every subsequent edge.  This is faster than
    /// ear-clipping and always produces valid triangles for convex loops.
    ///
    /// Returns the number of holes filled.
    pub fn fill_holes(&mut self) -> usize {
        let loops = find_boundary_loops(&self.triangles);
        let n_holes = loops.len();
        for lp in &loops {
            if lp.len() < 3 {
                continue;
            }
            let hub = lp[0];
            for i in 1..(lp.len() - 1) {
                self.triangles.push([hub, lp[i], lp[i + 1]]);
            }
        }
        n_holes
    }
    /// Fill all open boundary holes using ear-clipping (higher quality).
    ///
    /// Returns the number of holes filled.
    pub fn fill_holes_ear(&mut self) -> usize {
        let loops = find_boundary_loops(&self.triangles);
        let n_holes = loops.len();
        for lp in &loops {
            let new_tris = fill_hole_ear_clipping(&self.vertices, lp);
            self.triangles.extend_from_slice(&new_tris);
        }
        n_holes
    }
    /// Merge all vertex pairs closer than `epsilon` and compact the vertex
    /// array.  Triangle references are remapped accordingly.
    ///
    /// Returns the number of vertices that were merged (eliminated).
    pub fn remove_duplicate_vertices(&mut self, epsilon: f64) -> usize {
        let (new_v, new_t, n_merged) =
            merge_duplicate_vertices(&self.vertices, &self.triangles, epsilon);
        self.vertices = new_v;
        self.triangles = new_t;
        n_merged
    }
    /// Remove triangles until every edge is shared by at most 2 faces.
    ///
    /// The greedy strategy keeps the first two triangles seen for any edge;
    /// later ones are discarded.  Returns the number of triangles removed.
    pub fn fix_non_manifold_edges(&mut self) -> usize {
        let before = self.triangles.len();
        self.triangles = fix_non_manifold_edges(&self.triangles);
        before - self.triangles.len()
    }
    /// Compute the Euler characteristic χ = V − E + F.
    ///
    /// For a closed orientable surface of genus g: χ = 2 − 2g.
    /// A sphere has χ = 2; a torus has χ = 0.
    pub fn compute_euler_characteristic(&self) -> i32 {
        let v = self.vertices.len();
        let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
        for tri in &self.triangles {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                edge_set.insert(key);
            }
        }
        let e = edge_set.len();
        let f = self.triangles.len();
        v as i32 - e as i32 + f as i32
    }
}
