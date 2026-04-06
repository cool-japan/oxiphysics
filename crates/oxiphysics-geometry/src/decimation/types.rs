//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use std::collections::{HashMap, HashSet};

/// A record of one progressive mesh simplification step (edge collapse).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EdgeCollapseRecord {
    /// The vertex that was removed (merged into `target`).
    pub removed: usize,
    /// The vertex that `removed` was merged into.
    pub target: usize,
    /// The triangles that were deleted by this collapse.
    pub deleted_triangles: Vec<[usize; 3]>,
}
/// Progressive mesh: stores the base (simplified) mesh and a list of refinement records.
#[allow(dead_code)]
pub struct ProgressiveMeshSimple {
    /// Current simplified mesh.
    pub current: SimpleMesh,
    /// Collapse history (oldest first).
    pub history: Vec<EdgeCollapseRecord>,
}
impl ProgressiveMeshSimple {
    /// Build a progressive mesh from a base mesh.
    pub fn new(mesh: SimpleMesh) -> Self {
        Self {
            current: mesh,
            history: Vec::new(),
        }
    }
    /// Perform one vertex-to-vertex edge collapse (vertex `src` → vertex `dst`).
    ///
    /// Merges `src` into `dst`, removes all triangles that had both `src` and `dst`,
    /// and replaces `src` with `dst` in all remaining triangles.
    pub fn collapse_edge(&mut self, src: usize, dst: usize) {
        let mut deleted = Vec::new();
        let mut remaining = Vec::new();
        for tri in &self.current.triangles {
            let has_src = tri.contains(&src);
            let has_dst = tri.contains(&dst);
            if has_src && has_dst {
                deleted.push(*tri);
            } else if has_src {
                let new_tri = tri.map(|v| if v == src { dst } else { v });
                remaining.push(new_tri);
            } else {
                remaining.push(*tri);
            }
        }
        self.history.push(EdgeCollapseRecord {
            removed: src,
            target: dst,
            deleted_triangles: deleted,
        });
        self.current.triangles = remaining;
    }
    /// Number of collapses performed.
    pub fn n_collapses(&self) -> usize {
        self.history.len()
    }
    /// Current triangle count.
    pub fn n_triangles(&self) -> usize {
        self.current.triangle_count()
    }
}
/// Progressive mesh: stores the sequence of collapses so the mesh can be
/// coarsened incrementally and the collapses can be reversed.
pub struct ProgressiveMesh {
    /// Current mesh state.
    pub mesh: SimpleMesh,
    /// Collapse history (oldest first).
    pub history: Vec<CollapseRecord>,
}
impl ProgressiveMesh {
    /// Wrap an existing mesh.
    pub fn new(mesh: SimpleMesh) -> Self {
        Self {
            mesh,
            history: Vec::new(),
        }
    }
    /// Perform one QEM-based edge collapse (cheapest edge) and record it.
    ///
    /// Returns `false` when the mesh has ≤ 3 triangles and cannot be reduced
    /// further.
    pub fn collapse_one(&mut self) -> bool {
        if self.mesh.triangle_count() <= 3 {
            return false;
        }
        let quadrics = compute_vertex_quadrics(&self.mesh);
        let mut edge_set: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for &[a, b, c] in &self.mesh.triangles {
            edge_set.insert((a.min(b), a.max(b)));
            edge_set.insert((b.min(c), b.max(c)));
            edge_set.insert((a.min(c), a.max(c)));
        }
        let mut best_cost = f64::MAX;
        let mut best_v1 = 0usize;
        let mut best_v2 = 0usize;
        let mut best_pos = [0.0f64; 3];
        for (v1, v2) in &edge_set {
            let (cost, pos) = edge_collapse_cost(
                self.mesh.vertices[*v1],
                self.mesh.vertices[*v2],
                &quadrics[*v1],
                &quadrics[*v2],
            );
            if cost < best_cost {
                best_cost = cost;
                best_v1 = *v1;
                best_v2 = *v2;
                best_pos = pos;
            }
        }
        let old_pos = self.mesh.vertices[best_v1];
        let record = CollapseRecord {
            v_kept: best_v1,
            v_removed: best_v2,
            old_pos,
            new_pos: best_pos,
        };
        collapse_edge(&mut self.mesh, best_v1, best_v2, best_pos);
        self.history.push(record);
        true
    }
    /// Reduce the mesh until it has at most `target_triangles` triangles.
    pub fn decimate_to(&mut self, target_triangles: usize) {
        while self.mesh.triangle_count() > target_triangles {
            if !self.collapse_one() {
                break;
            }
        }
    }
}
/// A record of one edge collapse operation (for progressive meshes / undo).
#[derive(Debug, Clone)]
pub struct CollapseRecord {
    /// The surviving vertex index.
    pub v_kept: usize,
    /// The removed vertex index.
    pub v_removed: usize,
    /// The position before the collapse.
    pub old_pos: [f64; 3],
    /// The new position assigned to `v_kept`.
    pub new_pos: [f64; 3],
}
/// A directed edge in the mesh.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshEdge {
    /// First vertex index.
    pub v0: usize,
    /// Second vertex index.
    pub v1: usize,
    /// QEM collapse cost.
    pub cost: f64,
}
/// A min-priority queue of edges sorted by QEM cost.
#[allow(dead_code)]
pub struct EdgePriorityQueue {
    pub(super) heap: std::collections::BinaryHeap<MeshEdge>,
}
impl EdgePriorityQueue {
    /// Create an empty queue.
    pub fn new() -> Self {
        Self {
            heap: std::collections::BinaryHeap::new(),
        }
    }
    /// Insert an edge.
    pub fn push(&mut self, edge: MeshEdge) {
        self.heap.push(edge);
    }
    /// Pop the lowest-cost edge.
    pub fn pop(&mut self) -> Option<MeshEdge> {
        self.heap.pop()
    }
    /// Number of edges.
    pub fn len(&self) -> usize {
        self.heap.len()
    }
    /// Whether the queue is empty.
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }
    /// Build a priority queue from all unique edges of a mesh with QEM costs.
    pub fn from_mesh(mesh: &SimpleMesh) -> Self {
        let mut q = Self::new();
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        for tri in &mesh.triangles {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])] {
                let key = if a < b { (a, b) } else { (b, a) };
                if seen.insert(key) {
                    let va = mesh.vertices[a];
                    let vb = mesh.vertices[b];
                    let cost =
                        (va[0] - vb[0]).powi(2) + (va[1] - vb[1]).powi(2) + (va[2] - vb[2]).powi(2);
                    q.push(MeshEdge {
                        v0: key.0,
                        v1: key.1,
                        cost,
                    });
                }
            }
        }
        q
    }
}
/// A simple triangle mesh with vertices and index triples.
#[derive(Clone)]
pub struct SimpleMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle index triples (counter-clockwise winding).
    pub triangles: Vec<[usize; 3]>,
}
impl SimpleMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
        }
    }
    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, v: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(v);
        idx
    }
    /// Add a triangle by vertex indices.
    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize) {
        self.triangles.push([a, b, c]);
    }
    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }
    /// Outward-facing unit normal of triangle `t_idx`.
    pub fn triangle_normal(&self, t_idx: usize) -> [f64; 3] {
        let [a, b, c] = self.triangles[t_idx];
        let va = self.vertices[a];
        let vb = self.vertices[b];
        let vc = self.vertices[c];
        let ab = [vb[0] - va[0], vb[1] - va[1], vb[2] - va[2]];
        let ac = [vc[0] - va[0], vc[1] - va[1], vc[2] - va[2]];
        let n = cross(ab, ac);
        normalize(n)
    }
    /// Area of triangle `t_idx`.
    pub fn triangle_area(&self, t_idx: usize) -> f64 {
        let [a, b, c] = self.triangles[t_idx];
        let va = self.vertices[a];
        let vb = self.vertices[b];
        let vc = self.vertices[c];
        let ab = [vb[0] - va[0], vb[1] - va[1], vb[2] - va[2]];
        let ac = [vc[0] - va[0], vc[1] - va[1], vc[2] - va[2]];
        let cp = cross(ab, ac);
        0.5 * length(cp)
    }
    /// Total surface area.
    pub fn surface_area(&self) -> f64 {
        (0..self.triangle_count())
            .map(|i| self.triangle_area(i))
            .sum()
    }
}
/// 4×4 symmetric quadric matrix used for QEM mesh simplification.
pub struct QuadricMatrix {
    /// Row-major 4×4 matrix entries.
    pub q: [[f64; 4]; 4],
}
impl QuadricMatrix {
    /// Zero matrix.
    pub fn zero() -> Self {
        Self { q: [[0.0; 4]; 4] }
    }
    /// Fundamental error quadric Kp = pp^T for a plane ax + by + cz + d = 0.
    /// The plane coefficients must satisfy a² + b² + c² = 1 (unit normal).
    pub fn from_plane(a: f64, b: f64, c: f64, d: f64) -> Self {
        let p = [a, b, c, d];
        let mut q = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                q[i][j] = p[i] * p[j];
            }
        }
        Self { q }
    }
    /// Element-wise sum of two quadrics.
    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                result.q[i][j] = self.q[i][j] + other.q[i][j];
            }
        }
        result
    }
    /// Vertex error v^T Q v in homogeneous coordinates \[x, y, z, 1\].
    pub fn vertex_error(&self, v: [f64; 3]) -> f64 {
        let h = [v[0], v[1], v[2], 1.0];
        let mut result = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                result += h[i] * self.q[i][j] * h[j];
            }
        }
        result
    }
}
/// A normal cone: represents the range of normals in a region.
///
/// Defined by an axis (average normal) and a half-angle (max deviation).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct NormalCone {
    /// Axis of the cone (unit vector).
    pub axis: [f64; 3],
    /// Half-angle of the cone in radians.
    pub half_angle: f64,
}
impl NormalCone {
    /// Create a normal cone from a single normal vector.
    pub fn from_normal(n: [f64; 3]) -> Self {
        Self {
            axis: normalize3(n),
            half_angle: 0.0,
        }
    }
    /// Merge two normal cones into a bounding cone.
    pub fn merge(&self, other: &NormalCone) -> NormalCone {
        let axis_new = normalize3(add3(self.axis, other.axis));
        let angle1 = dot3(self.axis, axis_new).clamp(-1.0, 1.0).acos() + self.half_angle;
        let angle2 = dot3(other.axis, axis_new).clamp(-1.0, 1.0).acos() + other.half_angle;
        NormalCone {
            axis: axis_new,
            half_angle: angle1.max(angle2),
        }
    }
    /// Check whether a new normal falls within this cone (with tolerance `tol`).
    pub fn contains(&self, n: [f64; 3], tol: f64) -> bool {
        let n_hat = normalize3(n);
        let cos_angle = dot3(self.axis, n_hat).clamp(-1.0, 1.0);
        let angle = cos_angle.acos();
        angle <= self.half_angle + tol
    }
    /// Convert half-angle to degrees.
    pub fn half_angle_deg(&self) -> f64 {
        self.half_angle.to_degrees()
    }
}
/// Feature detection result for an edge.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgeFeature {
    /// Smooth interior edge.
    Smooth,
    /// Sharp crease edge (dihedral angle > threshold).
    Crease,
    /// Boundary edge (belongs to only one triangle).
    Boundary,
}
/// Statistics from a decimation run.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct DecimationMetrics {
    /// Original vertex count.
    pub original_vertices: usize,
    /// Original triangle count.
    pub original_triangles: usize,
    /// Reduced vertex count.
    pub reduced_vertices: usize,
    /// Reduced triangle count.
    pub reduced_triangles: usize,
    /// Number of edge collapses performed.
    pub n_collapses: usize,
    /// Total QEM cost accumulated.
    pub total_qem_cost: f64,
    /// Maximum QEM cost in any single collapse.
    pub max_qem_cost: f64,
}
impl DecimationMetrics {
    /// Compute reduction ratio for vertices \[0, 1\].
    pub fn vertex_reduction_ratio(&self) -> f64 {
        if self.original_vertices == 0 {
            return 0.0;
        }
        1.0 - self.reduced_vertices as f64 / self.original_vertices as f64
    }
    /// Compute reduction ratio for triangles \[0, 1\].
    pub fn triangle_reduction_ratio(&self) -> f64 {
        if self.original_triangles == 0 {
            return 0.0;
        }
        1.0 - self.reduced_triangles as f64 / self.original_triangles as f64
    }
    /// Average QEM cost per collapse.
    pub fn avg_qem_cost(&self) -> f64 {
        if self.n_collapses == 0 {
            return 0.0;
        }
        self.total_qem_cost / self.n_collapses as f64
    }
}
/// A QEM-based decimator that wraps `SimpleMesh` and provides adaptive
/// error thresholds, boundary preservation, and feature scoring.
pub struct QemDecimation {
    /// The mesh being simplified.
    pub mesh: SimpleMesh,
    /// Per-vertex quadric matrices.
    pub(super) quadrics: Vec<QuadricMatrix>,
}
impl QemDecimation {
    /// Construct a `QemDecimation` from a mesh.  Quadrics are computed once.
    pub fn new(mesh: SimpleMesh) -> Self {
        let quadrics = compute_vertex_quadrics(&mesh);
        Self { mesh, quadrics }
    }
    /// Recompute the per-vertex quadric matrices from the current mesh.
    pub fn recompute_quadrics(&mut self) {
        self.quadrics = compute_vertex_quadrics(&self.mesh);
    }
    /// Compute an *adaptive* QEM error threshold.
    ///
    /// The threshold is `scale_factor * avg_edge_length² * mean_curvature_proxy`,
    /// where the mean curvature proxy is the average edge collapse cost over
    /// a sample of edges.  This adapts the threshold to mesh scale and feature
    /// density.
    ///
    /// `scale_factor` controls aggressiveness (try 0.01 – 1.0).
    pub fn compute_error_threshold(&self, scale_factor: f64) -> f64 {
        let avg_edge = average_edge_length(&self.mesh);
        if avg_edge < 1e-12 {
            return 0.0;
        }
        let scale_sq = avg_edge * avg_edge;
        let n_verts = self.mesh.vertex_count();
        if n_verts < 2 {
            return scale_factor * scale_sq;
        }
        let mut sampled_cost = 0.0f64;
        let mut sampled_count = 0usize;
        let step = (n_verts / 32).max(1);
        for tri in self.mesh.triangles.iter().step_by(step) {
            let v0 = tri[0];
            let v1 = tri[1];
            if v0 >= self.quadrics.len() || v1 >= self.quadrics.len() {
                continue;
            }
            let (cost, _) = edge_collapse_cost(
                self.mesh.vertices[v0],
                self.mesh.vertices[v1],
                &self.quadrics[v0],
                &self.quadrics[v1],
            );
            sampled_cost += cost;
            sampled_count += 1;
        }
        let mean_cost = if sampled_count > 0 {
            sampled_cost / sampled_count as f64
        } else {
            1.0
        };
        (scale_factor * scale_sq * mean_cost).max(1e-12)
    }
    /// Collapse edges while preserving boundary geometry.
    ///
    /// Boundary edges (shared by exactly one triangle) are never collapsed.
    /// Interior edges with QEM cost below `threshold` are collapsed greedily
    /// using a priority queue.
    ///
    /// Returns the number of collapses performed.
    pub fn preserve_boundary(&mut self, threshold: f64) -> usize {
        let boundary_edges = find_boundary_edges(&self.mesh);
        let boundary_set: HashSet<(usize, usize)> = boundary_edges
            .iter()
            .flat_map(|&(a, b)| [(a.min(b), a.max(b)), (a.min(b), a.max(b))])
            .collect();
        let mut n_collapses = 0usize;
        let mut keep_going = true;
        while keep_going {
            keep_going = false;
            self.recompute_quadrics();
            let mut best: Option<(f64, [f64; 3], usize, usize)> = None;
            for tri in &self.mesh.triangles {
                for k in 0..3 {
                    let v0 = tri[k];
                    let v1 = tri[(k + 1) % 3];
                    let key = (v0.min(v1), v0.max(v1));
                    if boundary_set.contains(&key) {
                        continue;
                    }
                    if v0 >= self.quadrics.len() || v1 >= self.quadrics.len() {
                        continue;
                    }
                    let (cost, opt) = edge_collapse_cost(
                        self.mesh.vertices[v0],
                        self.mesh.vertices[v1],
                        &self.quadrics[v0],
                        &self.quadrics[v1],
                    );
                    if cost <= threshold && best.as_ref().is_none_or(|b| cost < b.0) {
                        best = Some((cost, opt, v0, v1));
                    }
                }
            }
            if let Some((_cost, opt, v0, v1)) = best {
                self.mesh.vertices[v0] = opt;
                for tri in &mut self.mesh.triangles {
                    for idx in tri.iter_mut() {
                        if *idx == v1 {
                            *idx = v0;
                        }
                    }
                }
                self.mesh
                    .triangles
                    .retain(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);
                n_collapses += 1;
                keep_going = true;
            }
        }
        n_collapses
    }
    /// Compute a feature importance score for each edge in the mesh.
    ///
    /// The score is the sum of the dihedral-angle curvature between the two
    /// incident faces.  High scores indicate sharp creases or geometric
    /// features that should be preserved.
    ///
    /// Returns a `HashMap<(usize, usize), f64>` keyed by the canonical
    /// (min, max) vertex index pair.
    pub fn compute_feature_score(&self) -> HashMap<(usize, usize), f64> {
        let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (fi, tri) in self.mesh.triangles.iter().enumerate() {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = (a.min(b), a.max(b));
                edge_faces.entry(key).or_default().push(fi);
            }
        }
        let mut scores: HashMap<(usize, usize), f64> = HashMap::new();
        for (&edge, faces) in &edge_faces {
            if faces.len() != 2 {
                scores.insert(edge, f64::INFINITY);
                continue;
            }
            let n0 = self.mesh.triangle_normal(faces[0]);
            let n1 = self.mesh.triangle_normal(faces[1]);
            let cos_a = (dot3(n0, n1)).clamp(-1.0, 1.0);
            let score = 1.0 - cos_a;
            scores.insert(edge, score);
        }
        scores
    }
}
/// Statistics about a decimation run.
#[derive(Debug, Clone)]
pub struct DecimationStats {
    /// Vertex count before decimation.
    pub original_vertices: usize,
    /// Triangle count before decimation.
    pub original_triangles: usize,
    /// Vertex count after decimation.
    pub result_vertices: usize,
    /// Triangle count after decimation.
    pub result_triangles: usize,
    /// Ratio of removed triangles to original triangle count (0..1).
    pub reduction_ratio: f64,
    /// Total QEM error accumulated over all collapses.
    pub total_qem_error: f64,
}
impl DecimationStats {
    /// Compute stats from before/after snapshots and accumulated error.
    pub fn compute(
        original_v: usize,
        original_t: usize,
        result: &SimpleMesh,
        total_qem_error: f64,
    ) -> Self {
        let reduction_ratio = if original_t > 0 {
            1.0 - result.triangle_count() as f64 / original_t as f64
        } else {
            0.0
        };
        Self {
            original_vertices: original_v,
            original_triangles: original_t,
            result_vertices: result.vertex_count(),
            result_triangles: result.triangle_count(),
            reduction_ratio,
            total_qem_error,
        }
    }
}
