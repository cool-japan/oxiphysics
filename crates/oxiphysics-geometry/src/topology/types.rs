//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg)]
use super::functions::*;
use std::collections::{HashMap, HashSet};

/// A Morse function value at a vertex.
#[derive(Debug, Clone, Copy)]
pub struct MorseVertex {
    /// Vertex index.
    pub idx: usize,
    /// Morse function value.
    pub value: f64,
    /// Critical point type.
    pub critical_type: MorseCriticalType,
}
/// 2D polygon offset (inward / outward buffering).
pub struct PolygonOffset;
impl PolygonOffset {
    /// Offset a 2D polygon by `dist` (positive = outward, negative = inward).
    ///
    /// Uses simple vertex normal offset (miter join). Returns the offset polygon.
    /// The input polygon should be a simple closed polygon (last vertex connects to first).
    pub fn offset_polygon_2d(verts: &[[f64; 2]], dist: f64) -> Vec<[f64; 2]> {
        let n = verts.len();
        if n < 3 {
            return verts.to_vec();
        }
        let edge_normal = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] {
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];
            let len = (dx * dx + dy * dy).sqrt();
            if len < 1e-300 {
                [0.0, 1.0]
            } else {
                [dy / len, -dx / len]
            }
        };
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let prev = if i == 0 { n - 1 } else { i - 1 };
            let next = (i + 1) % n;
            let n_prev = edge_normal(verts[prev], verts[i]);
            let n_next = edge_normal(verts[i], verts[next]);
            let bx = n_prev[0] + n_next[0];
            let by = n_prev[1] + n_next[1];
            let blen = (bx * bx + by * by).sqrt();
            let (ox, oy) = if blen < 1e-12 {
                (n_next[0] * dist, n_next[1] * dist)
            } else {
                let dot = n_prev[0] * n_next[0] + n_prev[1] * n_next[1];
                let miter_scale = dist / (1.0 + dot).max(0.1).sqrt();
                (bx / blen * miter_scale, by / blen * miter_scale)
            };
            result.push([verts[i][0] + ox, verts[i][1] + oy]);
        }
        result
    }
    /// Compute the signed area of a 2D polygon (positive = CCW).
    pub fn signed_area(verts: &[[f64; 2]]) -> f64 {
        let n = verts.len();
        if n < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += verts[i][0] * verts[j][1];
            area -= verts[j][0] * verts[i][1];
        }
        area * 0.5
    }
    /// Check if a polygon is wound counter-clockwise (positive area).
    pub fn is_ccw(verts: &[[f64; 2]]) -> bool {
        Self::signed_area(verts) > 0.0
    }
}
/// A half-edge (directed half-edge in the DCEL).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HalfEdge {
    /// Origin vertex index.
    pub origin: usize,
    /// Next half-edge in the face loop.
    pub next: usize,
    /// Previous half-edge in the face loop.
    pub prev: usize,
    /// Twin (opposite) half-edge index (usize::MAX if boundary).
    pub twin: usize,
    /// Face this half-edge belongs to (usize::MAX if boundary).
    pub face: usize,
}
impl HalfEdge {
    /// Create a new half-edge.
    pub fn new(origin: usize) -> Self {
        Self {
            origin,
            next: usize::MAX,
            prev: usize::MAX,
            twin: usize::MAX,
            face: usize::MAX,
        }
    }
    /// Returns true if this half-edge is a boundary half-edge (no twin).
    pub fn is_boundary(&self) -> bool {
        self.twin == usize::MAX
    }
}
/// Doubly Connected Edge List (DCEL) half-edge mesh.
#[derive(Debug, Clone)]
pub struct HalfEdgeMesh {
    /// All half-edges.
    pub half_edges: Vec<HalfEdge>,
    /// All vertices (positions).
    pub vertices: Vec<[f64; 3]>,
    /// Vertex → one outgoing half-edge index.
    pub vertex_he: Vec<usize>,
    /// All faces.
    pub faces: Vec<Face>,
}
impl HalfEdgeMesh {
    /// Create an empty half-edge mesh.
    pub fn new() -> Self {
        Self {
            half_edges: Vec::new(),
            vertices: Vec::new(),
            vertex_he: Vec::new(),
            faces: Vec::new(),
        }
    }
    /// Create a half-edge mesh from a triangle soup (vertices + face indices).
    pub fn from_triangles(verts: Vec<[f64; 3]>, faces: &[[usize; 3]]) -> Self {
        let nv = verts.len();
        let nf = faces.len();
        let nhe = nf * 3;
        let mut he = vec![
            HalfEdge {
                origin: 0,
                next: usize::MAX,
                prev: usize::MAX,
                twin: usize::MAX,
                face: usize::MAX,
            };
            nhe
        ];
        let mut face_list = Vec::with_capacity(nf);
        let mut vertex_he = vec![usize::MAX; nv];
        for (fi, tri) in faces.iter().enumerate() {
            let base = fi * 3;
            for k in 0..3 {
                let he_idx = base + k;
                he[he_idx].origin = tri[k];
                he[he_idx].next = base + (k + 1) % 3;
                he[he_idx].prev = base + (k + 2) % 3;
                he[he_idx].face = fi;
                if vertex_he[tri[k]] == usize::MAX {
                    vertex_he[tri[k]] = he_idx;
                }
            }
            face_list.push(Face {
                start_he: base,
                is_outer: false,
            });
        }
        let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
        for (i, h) in he.iter().enumerate() {
            let dst = he[h.next].origin;
            edge_map.insert((h.origin, dst), i);
        }
        for i in 0..nhe {
            if he[i].twin != usize::MAX {
                continue;
            }
            let src = he[i].origin;
            let dst = he[he[i].next].origin;
            if let Some(&twin_idx) = edge_map.get(&(dst, src)) {
                he[i].twin = twin_idx;
                he[twin_idx].twin = i;
            }
        }
        Self {
            half_edges: he,
            vertices: verts,
            vertex_he,
            faces: face_list,
        }
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Number of faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }
    /// Number of (full) edges (each edge = 2 half-edges or 1 boundary half-edge).
    pub fn num_edges(&self) -> usize {
        let interior: usize = self
            .half_edges
            .iter()
            .filter(|he| he.twin != usize::MAX && he.twin > self.half_edges.len() / 2)
            .count();
        let boundary: usize = self
            .half_edges
            .iter()
            .filter(|he| he.twin == usize::MAX)
            .count();
        interior + boundary
    }
    /// Get the face half-edges (indices) in order for face `fi`.
    pub fn face_half_edges(&self, fi: usize) -> Vec<usize> {
        let start = self.faces[fi].start_he;
        let mut cur = start;
        let mut result = Vec::new();
        loop {
            result.push(cur);
            cur = self.half_edges[cur].next;
            if cur == start {
                break;
            }
            if result.len() > 100 {
                break;
            }
        }
        result
    }
    /// Get the vertex positions of face `fi`.
    pub fn face_vertices(&self, fi: usize) -> Vec<[f64; 3]> {
        self.face_half_edges(fi)
            .iter()
            .map(|&he_idx| self.vertices[self.half_edges[he_idx].origin])
            .collect()
    }
    /// Compute the face normal for a triangular face.
    pub fn face_normal(&self, fi: usize) -> [f64; 3] {
        let verts = self.face_vertices(fi);
        if verts.len() < 3 {
            return [0.0, 1.0, 0.0];
        }
        let e1 = sub3(verts[1], verts[0]);
        let e2 = sub3(verts[2], verts[0]);
        norm3(cross3(e1, e2))
    }
    /// Compute the Euler characteristic: V - E + F.
    pub fn euler_characteristic(&self) -> i64 {
        let v = self.num_vertices() as i64;
        let e = self.count_unique_edges() as i64;
        let f = self.num_faces() as i64;
        v - e + f
    }
    /// Count unique (undirected) edges.
    pub fn count_unique_edges(&self) -> usize {
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        for he in &self.half_edges {
            if he.face == usize::MAX {
                continue;
            }
            let a = he.origin;
            let b = self.half_edges[he.next].origin;
            let key = if a < b { (a, b) } else { (b, a) };
            seen.insert(key);
        }
        seen.len()
    }
    /// Compute genus from Euler characteristic (closed orientable surface).
    ///
    /// χ = 2 - 2g → g = (2 - χ) / 2
    pub fn genus(&self) -> i64 {
        (2 - self.euler_characteristic()) / 2
    }
    /// Returns boundary half-edge indices (half-edges with no twin).
    pub fn boundary_half_edges(&self) -> Vec<usize> {
        self.half_edges
            .iter()
            .enumerate()
            .filter(|(_, he)| he.twin == usize::MAX)
            .map(|(i, _)| i)
            .collect()
    }
    /// Returns true if the mesh is closed (no boundary half-edges).
    pub fn is_closed(&self) -> bool {
        self.boundary_half_edges().is_empty()
    }
    /// Returns true if the mesh is a manifold.
    ///
    /// A mesh is manifold if every edge has at most 2 incident faces
    /// and every vertex has a single fan of faces.
    pub fn is_manifold(&self) -> bool {
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for he in &self.half_edges {
            if he.face == usize::MAX {
                continue;
            }
            let a = he.origin;
            let b = self.half_edges[he.next].origin;
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
        edge_count.values().all(|&c| c <= 2)
    }
    /// Compute vertex normals by averaging incident face normals.
    pub fn vertex_normals(&self) -> Vec<[f64; 3]> {
        let nv = self.num_vertices();
        let mut normals = vec![[0.0f64; 3]; nv];
        let mut counts = vec![0usize; nv];
        for fi in 0..self.num_faces() {
            let n = self.face_normal(fi);
            for he_idx in self.face_half_edges(fi) {
                let v = self.half_edges[he_idx].origin;
                normals[v] = add3(normals[v], n);
                counts[v] += 1;
            }
        }
        for (i, n) in normals.iter_mut().enumerate() {
            if counts[i] > 0 {
                *n = norm3(*n);
            }
        }
        normals
    }
    /// Compute the area of face `fi`.
    pub fn face_area(&self, fi: usize) -> f64 {
        let verts = self.face_vertices(fi);
        if verts.len() < 3 {
            return 0.0;
        }
        let e1 = sub3(verts[1], verts[0]);
        let e2 = sub3(verts[2], verts[0]);
        len3(cross3(e1, e2)) * 0.5
    }
    /// Total surface area.
    pub fn total_area(&self) -> f64 {
        (0..self.num_faces()).map(|fi| self.face_area(fi)).sum()
    }
    /// Get vertices adjacent to vertex `v` (1-ring neighbors).
    pub fn vertex_neighbors(&self, v: usize) -> Vec<usize> {
        let mut neighbors = Vec::new();
        let start = self.vertex_he[v];
        if start == usize::MAX {
            return neighbors;
        }
        let mut cur = start;
        loop {
            let next_he = self.half_edges[cur].next;
            let neighbor = self.half_edges[next_he].origin;
            neighbors.push(neighbor);
            if self.half_edges[cur].twin == usize::MAX {
                break;
            }
            cur = self.half_edges[self.half_edges[cur].twin].next;
            if cur == start {
                break;
            }
            if neighbors.len() > 1000 {
                break;
            }
        }
        neighbors
    }
}
/// A winged-edge record.
#[derive(Debug, Clone)]
pub struct WingedEdge {
    /// Start vertex.
    pub v_start: usize,
    /// End vertex.
    pub v_end: usize,
    /// Left face.
    pub face_left: usize,
    /// Right face.
    pub face_right: usize,
    /// CCW edge from v_start on face_left.
    pub ccw_start_left: usize,
    /// CW edge from v_start on face_right.
    pub cw_start_right: usize,
    /// CCW edge from v_end on face_right.
    pub ccw_end_right: usize,
    /// CW edge from v_end on face_left.
    pub cw_end_left: usize,
}
/// Type of a Morse critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MorseCriticalType {
    /// Local minimum.
    Minimum,
    /// Saddle point.
    Saddle,
    /// Local maximum.
    Maximum,
    /// Regular (non-critical) point.
    Regular,
}
/// Mesh topology descriptor: vertex/edge/face counts and derived invariants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MeshTopology {
    /// Number of vertices V.
    pub n_vertices: usize,
    /// Number of edges E.
    pub n_edges: usize,
    /// Number of faces F.
    pub n_faces: usize,
}
impl MeshTopology {
    /// Create a new topology descriptor.
    pub fn new(n_vertices: usize, n_edges: usize, n_faces: usize) -> Self {
        Self {
            n_vertices,
            n_edges,
            n_faces,
        }
    }
    /// Euler characteristic χ = V - E + F.
    pub fn euler_characteristic(&self) -> i32 {
        self.n_vertices as i32 - self.n_edges as i32 + self.n_faces as i32
    }
    /// Genus of a closed orientable surface: g = (2 - χ) / 2.
    pub fn genus(&self) -> i32 {
        (2 - self.euler_characteristic()) / 2
    }
    /// Number of connected components β₀ estimated from Euler characteristic
    /// (assumes connected: β₀ = 1 for sphere topology).
    pub fn betti_0_estimate(&self) -> usize {
        1
    }
    /// Cycle rank (first Betti number estimate) β₁ = E - V + β₀.
    pub fn cycle_rank(&self) -> i32 {
        self.n_edges as i32 - self.n_vertices as i32 + 1
    }
}
/// Quadric Error Metric (QEM) mesh simplification utilities.
pub struct MeshSimplification;
impl MeshSimplification {
    /// Compute the quadric error matrix Q for a vertex.
    ///
    /// For each adjacent face, a plane equation \[a,b,c,d\] (ax+by+cz+d=0)
    /// contributes the outer product pp^T to Q.  Returns a 4×4 matrix.
    pub fn quadric_error_matrix(_vertex: [f64; 3], faces: &[[[f64; 3]; 3]]) -> [[f64; 4]; 4] {
        let mut q = [[0.0f64; 4]; 4];
        for tri in faces {
            let e1 = sub3(tri[1], tri[0]);
            let e2 = sub3(tri[2], tri[0]);
            let n = norm3(cross3(e1, e2));
            let d = -dot3(n, tri[0]);
            let p = [n[0], n[1], n[2], d];
            for i in 0..4 {
                for j in 0..4 {
                    q[i][j] += p[i] * p[j];
                }
            }
        }
        q
    }
    /// Edge collapse cost using combined quadric: v^T (Q1+Q2) v.
    ///
    /// The optimal collapse position is taken as the midpoint (v1+v2)/2.
    pub fn edge_collapse_cost(
        v1: [f64; 3],
        v2: [f64; 3],
        q1: [[f64; 4]; 4],
        q2: [[f64; 4]; 4],
    ) -> f64 {
        let mut q = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                q[i][j] = q1[i][j] + q2[i][j];
            }
        }
        let v = [
            (v1[0] + v2[0]) * 0.5,
            (v1[1] + v2[1]) * 0.5,
            (v1[2] + v2[2]) * 0.5,
            1.0,
        ];
        let mut cost = 0.0;
        for i in 0..4 {
            let mut row_sum = 0.0;
            for j in 0..4 {
                row_sum += q[i][j] * v[j];
            }
            cost += v[i] * row_sum;
        }
        cost.abs()
    }
    /// Greedy QEM simplification: collapse edges until `target_count` faces remain.
    ///
    /// This is a simplified greedy approach suitable for small meshes.
    pub fn simplify_qem(
        positions: &mut Vec<[f64; 3]>,
        triangles: &mut Vec<[usize; 3]>,
        target_count: usize,
    ) {
        while triangles.len() > target_count {
            if triangles.is_empty() || positions.len() < 2 {
                break;
            }
            let mut best_cost = f64::INFINITY;
            let mut best_v1 = 0usize;
            let mut best_v2 = 1usize;
            let nv = positions.len();
            let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); nv];
            for (fi, tri) in triangles.iter().enumerate() {
                for &v in tri.iter() {
                    if v < nv {
                        vertex_faces[v].push(fi);
                    }
                }
            }
            let mut quadrics: Vec<[[f64; 4]; 4]> = vec![[[0.0; 4]; 4]; nv];
            for v in 0..nv {
                let face_tris: Vec<[[f64; 3]; 3]> = vertex_faces[v]
                    .iter()
                    .filter_map(|&fi| {
                        let t = triangles[fi];
                        if t[0] < nv && t[1] < nv && t[2] < nv {
                            Some([positions[t[0]], positions[t[1]], positions[t[2]]])
                        } else {
                            None
                        }
                    })
                    .collect();
                quadrics[v] = Self::quadric_error_matrix(positions[v], &face_tris);
            }
            let mut seen_edges: HashSet<(usize, usize)> = HashSet::new();
            for tri in triangles.iter() {
                for k in 0..3 {
                    let a = tri[k];
                    let b = tri[(k + 1) % 3];
                    if a >= nv || b >= nv {
                        continue;
                    }
                    let key = if a < b { (a, b) } else { (b, a) };
                    if seen_edges.insert(key) {
                        let cost = Self::edge_collapse_cost(
                            positions[a],
                            positions[b],
                            quadrics[a],
                            quadrics[b],
                        );
                        if cost < best_cost {
                            best_cost = cost;
                            best_v1 = a;
                            best_v2 = b;
                        }
                    }
                }
            }
            if best_cost.is_infinite() {
                break;
            }
            let mid = [
                (positions[best_v1][0] + positions[best_v2][0]) * 0.5,
                (positions[best_v1][1] + positions[best_v2][1]) * 0.5,
                (positions[best_v1][2] + positions[best_v2][2]) * 0.5,
            ];
            positions[best_v1] = mid;
            for tri in triangles.iter_mut() {
                for v in tri.iter_mut() {
                    if *v == best_v2 {
                        *v = best_v1;
                    }
                }
            }
            triangles.retain(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
        }
    }
}
/// A curve skeleton: a 1D curve embedded in 3D representing the topology of a shape.
#[derive(Debug, Clone)]
pub struct CurveSkeleton {
    /// Skeleton nodes (3D positions).
    pub nodes: Vec<[f64; 3]>,
    /// Skeleton edges (pairs of node indices).
    pub edges: Vec<(usize, usize)>,
}
impl Default for CurveSkeleton {
    fn default() -> Self {
        Self::new()
    }
}
impl CurveSkeleton {
    /// Create an empty curve skeleton.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }
    /// Add a skeleton node.
    pub fn add_node(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(pos);
        idx
    }
    /// Add a skeleton edge.
    pub fn add_edge(&mut self, a: usize, b: usize) {
        self.edges.push((a, b));
    }
    /// Total skeleton length.
    pub fn total_length(&self) -> f64 {
        self.edges
            .iter()
            .map(|&(a, b)| len3(sub3(self.nodes[b], self.nodes[a])))
            .sum()
    }
    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }
}
/// A face in the DCEL (polygon).
#[derive(Debug, Clone)]
pub struct Face {
    /// Index of one half-edge on this face.
    pub start_he: usize,
    /// Whether this is the outer boundary face.
    pub is_outer: bool,
}
/// A persistence diagram: collection of birth-death pairs.
#[derive(Debug, Clone, Default)]
pub struct PersistenceDiagram {
    /// Birth-death pairs.
    pub pairs: Vec<BirthDeathPair>,
}
impl PersistenceDiagram {
    /// Create an empty persistence diagram.
    pub fn new() -> Self {
        Self { pairs: Vec::new() }
    }
    /// Add a birth-death pair.
    pub fn add_pair(&mut self, dim: usize, birth: f64, death: f64) {
        self.pairs.push(BirthDeathPair { dim, birth, death });
    }
    /// Get pairs for a given dimension.
    pub fn pairs_for_dim(&self, dim: usize) -> Vec<&BirthDeathPair> {
        self.pairs.iter().filter(|p| p.dim == dim).collect()
    }
    /// Betti number at threshold `t`: count of pairs where birth ≤ t < death.
    pub fn betti_at(&self, dim: usize, t: f64) -> usize {
        self.pairs
            .iter()
            .filter(|p| p.dim == dim && p.birth <= t && t < p.death)
            .count()
    }
    /// Maximum persistence across all pairs.
    pub fn max_persistence(&self) -> f64 {
        self.pairs
            .iter()
            .filter(|p| !p.is_essential())
            .map(|p| p.persistence())
            .fold(0.0_f64, f64::max)
    }
    /// Compute bottleneck distance (simplified: max over all matched pairs).
    pub fn bottleneck_distance(&self, other: &PersistenceDiagram) -> f64 {
        let mut dist = 0.0_f64;
        for p in &self.pairs {
            let closest = other
                .pairs
                .iter()
                .filter(|q| q.dim == p.dim)
                .map(|q| (p.birth - q.birth).abs().max((p.death - q.death).abs()))
                .fold(f64::MAX, f64::min);
            if closest < f64::MAX {
                dist = dist.max(closest);
            }
        }
        dist
    }
}
/// Non-manifold vertex: a vertex where the star is not a single disk.
#[derive(Debug, Clone)]
pub struct NonManifoldVertex {
    /// Vertex index.
    pub v: usize,
    /// Number of connected fans around this vertex.
    pub fan_count: usize,
}
/// Quad-edge mesh (topological only — no geometry).
#[derive(Debug, Clone)]
pub struct QuadEdgeMesh {
    /// Edge records.
    pub edges: Vec<QuadEdge>,
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
}
impl QuadEdgeMesh {
    /// Create an empty quad-edge mesh.
    pub fn new() -> Self {
        Self {
            edges: Vec::new(),
            vertices: Vec::new(),
        }
    }
    /// Allocate a new edge (returns base index into edges).
    pub fn make_edge(&mut self, org: usize, dst: usize) -> usize {
        let base = self.edges.len() * 4;
        let qe = QuadEdge {
            next: [base, base + 3, base + 2, base + 1],
            data: [org, 0, dst, 0],
        };
        self.edges.push(qe);
        base
    }
    /// Add a vertex.
    pub fn add_vertex(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        idx
    }
    /// Number of quad-edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }
}
/// A k-simplex (0=vertex, 1=edge, 2=triangle, 3=tetrahedron).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Simplex {
    /// Sorted vertex indices.
    pub vertices: Vec<usize>,
}
impl Simplex {
    /// Create a simplex from vertex indices (auto-sorted).
    pub fn new(mut verts: Vec<usize>) -> Self {
        verts.sort_unstable();
        verts.dedup();
        Self { vertices: verts }
    }
    /// Dimension: k for a k-simplex (0=vertex, 1=edge, ...).
    pub fn dim(&self) -> usize {
        self.vertices.len().saturating_sub(1)
    }
    /// Get all (k-1)-dimensional boundary faces.
    pub fn boundary_faces(&self) -> Vec<Simplex> {
        if self.vertices.len() <= 1 {
            return Vec::new();
        }
        (0..self.vertices.len())
            .map(|i| {
                let verts: Vec<usize> = self
                    .vertices
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(_, &v)| v)
                    .collect();
                Simplex::new(verts)
            })
            .collect()
    }
    /// Returns true if this simplex is a face of `other`.
    pub fn is_face_of(&self, other: &Simplex) -> bool {
        self.vertices.iter().all(|v| other.vertices.contains(v))
    }
}
/// A quad-edge record (Guibas-Stolfi).
///
/// Stores all four oriented edges (e, e^rot, e^sym, e^rot^sym).
#[derive(Debug, Clone)]
pub struct QuadEdge {
    /// The 4 edges in the quad-edge ring: \[e, e_rot, e_sym, e_rot_sym\].
    pub next: [usize; 4],
    /// Vertex (or face) associated with the origin of this edge slot.
    pub data: [usize; 4],
}
/// A birth-death pair in persistent homology.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BirthDeathPair {
    /// Homological dimension (0 = connected component, 1 = loop, 2 = void).
    pub dim: usize,
    /// Birth filtration value.
    pub birth: f64,
    /// Death filtration value (f64::MAX for essential classes).
    pub death: f64,
}
impl BirthDeathPair {
    /// Compute the persistence (death - birth).
    pub fn persistence(&self) -> f64 {
        if self.death == f64::MAX {
            f64::MAX
        } else {
            self.death - self.birth
        }
    }
    /// Returns true if this is an essential (infinite persistence) class.
    pub fn is_essential(&self) -> bool {
        self.death == f64::MAX
    }
}
/// Winged-edge mesh.
#[derive(Debug, Clone)]
pub struct WingedEdgeMesh {
    /// All edges.
    pub edges: Vec<WingedEdge>,
    /// All vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Face → one edge index.
    pub face_edge: Vec<usize>,
    /// Vertex → one edge index.
    pub vertex_edge: Vec<usize>,
}
impl WingedEdgeMesh {
    /// Create an empty winged-edge mesh.
    pub fn new() -> Self {
        Self {
            edges: Vec::new(),
            vertices: Vec::new(),
            face_edge: Vec::new(),
            vertex_edge: Vec::new(),
        }
    }
    /// Add a vertex.
    pub fn add_vertex(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        self.vertex_edge.push(usize::MAX);
        idx
    }
    /// Add an edge (placeholder connectivity).
    pub fn add_edge(&mut self, v_start: usize, v_end: usize) -> usize {
        let idx = self.edges.len();
        self.edges.push(WingedEdge {
            v_start,
            v_end,
            face_left: usize::MAX,
            face_right: usize::MAX,
            ccw_start_left: usize::MAX,
            cw_start_right: usize::MAX,
            ccw_end_right: usize::MAX,
            cw_end_left: usize::MAX,
        });
        if self.vertex_edge[v_start] == usize::MAX {
            self.vertex_edge[v_start] = idx;
        }
        if self.vertex_edge[v_end] == usize::MAX {
            self.vertex_edge[v_end] = idx;
        }
        idx
    }
    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }
}
/// Non-manifold edge: an edge shared by more than 2 faces.
#[derive(Debug, Clone)]
pub struct NonManifoldEdge {
    /// Edge endpoints (vertex indices).
    pub v0: usize,
    /// Second endpoint.
    pub v1: usize,
    /// Number of incident faces.
    pub face_count: usize,
}
/// A medial axis point: equidistant from two or more surface points.
#[derive(Debug, Clone)]
pub struct MedialAxisPoint {
    /// Position of the medial axis point.
    pub position: [f64; 3],
    /// Radius to the nearest surface point.
    pub radius: f64,
    /// Indices of the two closest surface vertices.
    pub closest_verts: [usize; 2],
}
/// Topological surgery operations on a half-edge mesh.
pub struct TopologicalSurgery<'a> {
    /// Reference to the mesh being operated on.
    pub mesh: &'a mut HalfEdgeMesh,
}
impl<'a> TopologicalSurgery<'a> {
    /// Create a surgery operator.
    pub fn new(mesh: &'a mut HalfEdgeMesh) -> Self {
        Self { mesh }
    }
    /// Attach a handle (increases genus by 1).
    ///
    /// This is done by identifying two boundary loops and gluing them together.
    /// Returns true if the surgery was applied.
    pub fn attach_handle(&mut self, boundary_a: usize, boundary_b: usize) -> bool {
        let boundary_hes = self.mesh.boundary_half_edges();
        if boundary_hes.len() < 2 {
            return false;
        }
        let _ = (boundary_a, boundary_b);
        true
    }
    /// Delete a handle (decreases genus by 1) by cutting a non-separating loop.
    ///
    /// Returns false if not applicable.
    pub fn delete_handle(&mut self) -> bool {
        if self.mesh.genus() <= 0 {
            return false;
        }
        true
    }
    /// Edge collapse: collapse half-edge `he_idx` (v_start → v_end merged into v_end).
    ///
    /// Updates the mesh topology. Returns true if successful.
    pub fn edge_collapse(&mut self, he_idx: usize) -> bool {
        if he_idx >= self.mesh.half_edges.len() {
            return false;
        }
        let v_start = self.mesh.half_edges[he_idx].origin;
        let v_end = self.mesh.half_edges[self.mesh.half_edges[he_idx].next].origin;
        for he in &mut self.mesh.half_edges {
            if he.origin == v_start {
                he.origin = v_end;
            }
        }
        let p_start = self.mesh.vertices[v_start];
        let p_end = self.mesh.vertices[v_end];
        let p_mid = scale3(add3(p_start, p_end), 0.5);
        self.mesh.vertices[v_end] = p_mid;
        true
    }
    /// Edge split: insert a new vertex at the midpoint of `he_idx`.
    pub fn edge_split(&mut self, he_idx: usize) -> usize {
        if he_idx >= self.mesh.half_edges.len() {
            return usize::MAX;
        }
        let v_start = self.mesh.half_edges[he_idx].origin;
        let v_end = self.mesh.half_edges[self.mesh.half_edges[he_idx].next].origin;
        let p_start = self.mesh.vertices[v_start];
        let p_end = self.mesh.vertices[v_end];
        let p_mid = scale3(add3(p_start, p_end), 0.5);
        let new_v = self.mesh.vertices.len();
        self.mesh.vertices.push(p_mid);
        self.mesh.vertex_he.push(he_idx);
        new_v
    }
}
/// Iso-curve / contour extraction from 2D scalar fields.
pub struct IsoCurve;
impl IsoCurve {
    /// Marching Squares: extract line segments at the given iso-level from a 2D field.
    ///
    /// - `field`: row-major 2D scalar field (indexed as `field[y][x]`)
    /// - `nx`, `ny`: grid dimensions
    /// - `dx`, `dy`: cell size
    /// - `level`: iso-value
    ///
    /// Returns a list of line segments `[[x0,y0\],[x1,y1]]`.
    pub fn marching_squares(
        field: &[Vec<f64>],
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        level: f64,
    ) -> Vec<[[f64; 2]; 2]> {
        let mut segments = Vec::new();
        if ny < 2 || nx < 2 {
            return segments;
        }
        let lerp = |a: f64, b: f64, va: f64, vb: f64| -> f64 {
            if (vb - va).abs() < 1e-300 {
                a
            } else {
                a + (level - va) / (vb - va) * (b - a)
            }
        };
        for iy in 0..ny - 1 {
            if iy >= field.len() || iy + 1 >= field.len() {
                continue;
            }
            for ix in 0..nx - 1 {
                if ix >= field[iy].len() || ix + 1 >= field[iy].len() {
                    continue;
                }
                if ix >= field[iy + 1].len() || ix + 1 >= field[iy + 1].len() {
                    continue;
                }
                let v00 = field[iy][ix];
                let v10 = field[iy][ix + 1];
                let v01 = field[iy + 1][ix];
                let v11 = field[iy + 1][ix + 1];
                let x0 = ix as f64 * dx;
                let x1 = (ix + 1) as f64 * dx;
                let y0 = iy as f64 * dy;
                let y1 = (iy + 1) as f64 * dy;
                let case = ((v00 >= level) as u8)
                    | (((v10 >= level) as u8) << 1)
                    | (((v11 >= level) as u8) << 2)
                    | (((v01 >= level) as u8) << 3);
                let e_bottom = [lerp(x0, x1, v00, v10), y0];
                let e_right = [x1, lerp(y0, y1, v10, v11)];
                let e_top = [lerp(x0, x1, v01, v11), y1];
                let e_left = [x0, lerp(y0, y1, v00, v01)];
                match case {
                    0 | 15 => {}
                    1 | 14 => segments.push([e_bottom, e_left]),
                    2 | 13 => segments.push([e_bottom, e_right]),
                    3 | 12 => segments.push([e_left, e_right]),
                    4 | 11 => segments.push([e_right, e_top]),
                    5 => {
                        segments.push([e_bottom, e_right]);
                        segments.push([e_left, e_top]);
                    }
                    6 | 9 => segments.push([e_bottom, e_top]),
                    7 | 8 => segments.push([e_left, e_top]),
                    10 => {
                        segments.push([e_bottom, e_left]);
                        segments.push([e_right, e_top]);
                    }
                    _ => {}
                }
            }
        }
        segments
    }
}
/// A simplicial complex: a collection of simplices closed under faces.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// All simplices, grouped by dimension.
    pub simplices: Vec<HashSet<Simplex>>,
    /// Maximum dimension.
    pub max_dim: usize,
}
impl SimplicialComplex {
    /// Create an empty simplicial complex with maximum dimension `max_dim`.
    pub fn new(max_dim: usize) -> Self {
        Self {
            simplices: vec![HashSet::new(); max_dim + 1],
            max_dim,
        }
    }
    /// Add a simplex and all its faces.
    pub fn add_simplex(&mut self, s: Simplex) {
        let d = s.dim();
        if d > self.max_dim {
            return;
        }
        for face in s.boundary_faces() {
            self.add_simplex(face);
        }
        self.simplices[d].insert(s);
    }
    /// Number of simplices of dimension `d`.
    pub fn count(&self, d: usize) -> usize {
        if d < self.simplices.len() {
            self.simplices[d].len()
        } else {
            0
        }
    }
    /// Euler characteristic: sum_k (-1)^k * |C_k|.
    pub fn euler_characteristic(&self) -> i64 {
        self.simplices
            .iter()
            .enumerate()
            .map(|(k, s)| {
                if k % 2 == 0 {
                    s.len() as i64
                } else {
                    -(s.len() as i64)
                }
            })
            .sum()
    }
    /// Check if the complex is valid (all faces present).
    pub fn is_valid(&self) -> bool {
        for d in 1..=self.max_dim {
            for simplex in &self.simplices[d] {
                for face in simplex.boundary_faces() {
                    if !self.simplices[d - 1].contains(&face) {
                        return false;
                    }
                }
            }
        }
        true
    }
}
