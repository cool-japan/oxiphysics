//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Status of a vertex in the FMM algorithm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum FmmStatus {
    /// Vertex has not been reached yet.
    Far,
    /// Vertex is in the narrow band (tentative distance).
    Trial,
    /// Vertex distance is finalized.
    Known,
}
/// Parameters for the heat method.
#[derive(Debug, Clone)]
pub struct HeatMethodParams {
    /// Time step factor: `t = factor * mean_edge_length^2`.
    pub time_factor: f64,
}
/// A geodesic Voronoi cell on a mesh surface.
#[derive(Debug, Clone)]
pub struct GeodesicVoronoiCell {
    /// The source vertex (site) of this cell.
    pub site: usize,
    /// Indices of vertices belonging to this cell.
    pub vertices: Vec<usize>,
    /// Total area of the cell (sum of Voronoi areas of its vertices).
    pub area: f64,
}
/// Result of a discrete exponential map computation.
#[derive(Debug, Clone)]
pub struct ExpMapResult {
    /// The target vertex index on the mesh (nearest to the mapped point).
    pub vertex: usize,
    /// The 3D position of the mapped point on the surface.
    pub position: [f64; 3],
    /// The face that contains the mapped point.
    pub face: usize,
    /// Barycentric coordinates within the face.
    pub bary: [f64; 3],
}
/// Result of a discrete logarithmic map.
#[derive(Debug, Clone)]
pub struct LogMapResult {
    /// 2D coordinates in the tangent plane at the source vertex.
    pub coords: [f64; 2],
    /// Geodesic distance from source to target.
    pub distance: f64,
    /// Angle in the tangent plane (radians).
    pub angle: f64,
}
/// A node in the priority queue for distance computations.
#[derive(Debug, Clone)]
pub(super) struct DistNode {
    /// Vertex index.
    pub(super) vertex: usize,
    /// Distance from source.
    pub(super) dist: f64,
}
/// Triangle mesh for geodesic geometry computations.
///
/// Stores vertex positions and triangle connectivity. All algorithms in this
/// module operate on this structure.
#[derive(Debug, Clone)]
pub struct GeodesicMesh {
    /// Vertex positions in 3-space.
    pub vertices: Vec<[f64; 3]>,
    /// Triangular faces as index triples `[v0, v1, v2]`.
    pub faces: Vec<[usize; 3]>,
}
impl GeodesicMesh {
    /// Create a new geodesic mesh.
    pub fn new(vertices: Vec<[f64; 3]>, faces: Vec<[usize; 3]>) -> Self {
        Self { vertices, faces }
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Number of faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }
    /// Build adjacency: for each vertex, the list of neighbouring vertices.
    pub fn vertex_adjacency(&self) -> Vec<Vec<usize>> {
        let n = self.num_vertices();
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for f in &self.faces {
            for i in 0..3 {
                let a = f[i];
                let b = f[(i + 1) % 3];
                if !adj[a].contains(&b) {
                    adj[a].push(b);
                }
                if !adj[b].contains(&a) {
                    adj[b].push(a);
                }
            }
        }
        adj
    }
    /// Build face adjacency: for each vertex, which faces contain it.
    pub fn vertex_to_faces(&self) -> Vec<Vec<usize>> {
        let n = self.num_vertices();
        let mut v2f: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (fi, f) in self.faces.iter().enumerate() {
            for &vi in f {
                v2f[vi].push(fi);
            }
        }
        v2f
    }
    /// Compute the face normal of a triangle (unnormalized).
    pub fn face_normal(&self, fi: usize) -> [f64; 3] {
        let f = self.faces[fi];
        let e1 = sub3(self.vertices[f[1]], self.vertices[f[0]]);
        let e2 = sub3(self.vertices[f[2]], self.vertices[f[0]]);
        cross3(e1, e2)
    }
    /// Compute area-weighted vertex normal.
    pub fn vertex_normal(&self, vi: usize) -> [f64; 3] {
        let v2f = self.vertex_to_faces();
        let mut n = [0.0; 3];
        for &fi in &v2f[vi] {
            let fn_ = self.face_normal(fi);
            n = add3(n, fn_);
        }
        normalize3(n)
    }
    /// Compute face area.
    pub fn face_area(&self, fi: usize) -> f64 {
        0.5 * len3(self.face_normal(fi))
    }
    /// Total surface area.
    pub fn total_area(&self) -> f64 {
        (0..self.num_faces()).map(|fi| self.face_area(fi)).sum()
    }
    /// Build a local tangent frame at a vertex: returns `(tangent, bitangent, normal)`.
    pub fn tangent_frame(&self, vi: usize) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let n = self.vertex_normal(vi);
        let adj = self.vertex_adjacency();
        let t = if !adj[vi].is_empty() {
            let edge = sub3(self.vertices[adj[vi][0]], self.vertices[vi]);
            let proj = sub3(edge, scale3(n, dot3(edge, n)));
            normalize3(proj)
        } else {
            let candidate = if n[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let proj = sub3(candidate, scale3(n, dot3(candidate, n)));
            normalize3(proj)
        };
        let b = cross3(n, t);
        (t, b, n)
    }
    /// Cotangent weight for edge (i, j). Returns the sum of cot(alpha) + cot(beta)
    /// where alpha, beta are the angles opposite to edge ij in the two adjacent triangles.
    pub fn cotangent_weight(&self, i: usize, j: usize) -> f64 {
        let mut cot_sum = 0.0;
        for f in &self.faces {
            let mut found = false;
            let mut opposite = 0;
            for k in 0..3 {
                let a = f[k];
                let b = f[(k + 1) % 3];
                let c = f[(k + 2) % 3];
                if (a == i && b == j) || (a == j && b == i) {
                    found = true;
                    opposite = c;
                    break;
                }
            }
            if found {
                let pi = self.vertices[i];
                let pj = self.vertices[j];
                let po = self.vertices[opposite];
                let ei = sub3(pi, po);
                let ej = sub3(pj, po);
                let cos_a = dot3(ei, ej) / (len3(ei) * len3(ej) + 1e-30);
                let sin_a = len3(cross3(ei, ej)) / (len3(ei) * len3(ej) + 1e-30);
                if sin_a.abs() > 1e-15 {
                    cot_sum += cos_a / sin_a;
                }
            }
        }
        cot_sum
    }
    /// Voronoi area at vertex `vi` (mixed area formulation).
    pub fn voronoi_area(&self, vi: usize) -> f64 {
        let v2f = self.vertex_to_faces();
        let mut area = 0.0;
        for &fi in &v2f[vi] {
            let f = self.faces[fi];
            let local_idx = f.iter().position(|&v| v == vi).expect("element must exist");
            let a = self.vertices[f[local_idx]];
            let b = self.vertices[f[(local_idx + 1) % 3]];
            let c = self.vertices[f[(local_idx + 2) % 3]];
            let angle_a = angle_at_vertex(b, a, c);
            let angle_b = angle_at_vertex(a, b, c);
            let angle_c = angle_at_vertex(a, c, b);
            if angle_a > std::f64::consts::FRAC_PI_2 {
                area += self.face_area(fi) * 0.5;
            } else if angle_b > std::f64::consts::FRAC_PI_2 || angle_c > std::f64::consts::FRAC_PI_2
            {
                area += self.face_area(fi) * 0.25;
            } else {
                let ab2 = dot3(sub3(b, a), sub3(b, a));
                let ac2 = dot3(sub3(c, a), sub3(c, a));
                let cot_b = angle_b.cos() / (angle_b.sin() + 1e-30);
                let cot_c = angle_c.cos() / (angle_c.sin() + 1e-30);
                area += 0.125 * (cot_b * ac2 + cot_c * ab2);
            }
        }
        area.max(1e-30)
    }
}
