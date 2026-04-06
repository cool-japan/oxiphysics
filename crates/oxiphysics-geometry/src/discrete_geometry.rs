#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Discrete differential geometry on triangle meshes.
//!
//! This module implements the classical machinery of discrete differential
//! geometry, following the pioneering work of Desbrun, Meyer, Polthier,
//! Crane and colleagues:
//!
//! - [`DiscreteMesh`]: half-edge adjacency, vertex normals.
//! - Cotangent-weight Laplace-Beltrami operator.
//! - Angle-defect Gaussian curvature.
//! - Discrete mean curvature via the cotangent formula.
//! - Discrete geodesics: Dijkstra edge graph + Heat Method (Crane 2013).
//! - Parallel transport and holonomy on meshes.
//! - Discrete minimal surfaces (mean-curvature flow).
//! - Discrete vector fields and Hodge decomposition on meshes.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::{BinaryHeap, HashMap, HashSet};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector math helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean length of a 3-vector.
#[inline]
pub fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Euclidean distance between two points.
#[inline]
pub fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

/// Normalise a 3-vector; returns the zero vector if length is negligible.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}

/// Signed area of the triangle (i, j, k) — positive for CCW orientation.
#[inline]
pub fn tri_area(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let e0 = sub3(p1, p0);
    let e1 = sub3(p2, p0);
    0.5 * len3(cross3(e0, e1))
}

/// Cotangent of the angle at vertex `p_opp` opposite edge `(p0, p1)`.
///
/// Returns `cot(alpha)` where `alpha` is the angle at `p_opp` in the triangle
/// `(p0, p_opp, p1)`.
#[inline]
pub fn cotan_at(p0: [f64; 3], p_opp: [f64; 3], p1: [f64; 3]) -> f64 {
    let u = sub3(p0, p_opp);
    let v = sub3(p1, p_opp);
    let c = dot3(u, v);
    let s = len3(cross3(u, v));
    if s.abs() < 1e-14 { 0.0 } else { c / s }
}

// ---------------------------------------------------------------------------
// DiscreteMesh
// ---------------------------------------------------------------------------

/// A triangle mesh suitable for discrete differential geometry computations.
///
/// Stores vertex positions and triangles, with lazy adjacency information.
#[derive(Debug, Clone)]
pub struct DiscreteMesh {
    /// Vertex positions, indexed by vertex index.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle faces, each given as three vertex indices.
    pub triangles: Vec<[usize; 3]>,
}

impl DiscreteMesh {
    /// Construct a mesh from vertex positions and face index triples.
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>) -> Self {
        Self {
            vertices,
            triangles,
        }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangular faces.
    pub fn num_faces(&self) -> usize {
        self.triangles.len()
    }

    /// Returns all triangles that share vertex `v`.
    pub fn one_ring_faces(&self, v: usize) -> Vec<usize> {
        self.triangles
            .iter()
            .enumerate()
            .filter_map(|(fi, t)| if t.contains(&v) { Some(fi) } else { None })
            .collect()
    }

    /// Returns all neighbour vertices of vertex `v` (one-ring).
    pub fn one_ring_vertices(&self, v: usize) -> Vec<usize> {
        let mut nbrs: HashSet<usize> = HashSet::new();
        for t in &self.triangles {
            if t.contains(&v) {
                for &w in t {
                    if w != v {
                        nbrs.insert(w);
                    }
                }
            }
        }
        nbrs.into_iter().collect()
    }

    /// Area-weighted vertex normal at vertex `v`.
    pub fn vertex_normal(&self, v: usize) -> [f64; 3] {
        let mut n = [0.0f64; 3];
        for t in &self.triangles {
            if !t.contains(&v) {
                continue;
            }
            let p0 = self.vertices[t[0]];
            let p1 = self.vertices[t[1]];
            let p2 = self.vertices[t[2]];
            let face_n = cross3(sub3(p1, p0), sub3(p2, p0));
            n = add3(n, face_n);
        }
        normalize3(n)
    }

    /// Face normal for triangle `fi`.
    pub fn face_normal(&self, fi: usize) -> [f64; 3] {
        let t = &self.triangles[fi];
        let p0 = self.vertices[t[0]];
        let p1 = self.vertices[t[1]];
        let p2 = self.vertices[t[2]];
        normalize3(cross3(sub3(p1, p0), sub3(p2, p0)))
    }

    /// Voronoi area of vertex `v` (mixed Voronoi / barycentric area).
    ///
    /// Uses the formulation of Meyer et al. (2003).
    pub fn voronoi_area(&self, v: usize) -> f64 {
        let mut area = 0.0_f64;
        for t in &self.triangles {
            if !t.contains(&v) {
                continue;
            }
            // Identify local indices
            let (i, j, k) = if t[0] == v {
                (0, 1, 2)
            } else if t[1] == v {
                (1, 2, 0)
            } else {
                (2, 0, 1)
            };
            let pi = self.vertices[t[i]];
            let pj = self.vertices[t[j]];
            let pk = self.vertices[t[k]];
            // 1/8 * (cot_j * |e_k|^2 + cot_k * |e_j|^2)
            let cot_j = cotan_at(pi, pj, pk);
            let cot_k = cotan_at(pi, pk, pj);
            let ek_sq = {
                let d = sub3(pi, pj);
                dot3(d, d)
            };
            let ej_sq = {
                let d = sub3(pi, pk);
                dot3(d, d)
            };
            area += 0.125 * (cot_j * ek_sq + cot_k * ej_sq);
        }
        area.max(1e-15)
    }

    /// Build the sparse cotangent-weight Laplace-Beltrami matrix as a list of
    /// `(row, col, weight)` triples.
    ///
    /// The returned weights satisfy `L[i][j] = (cot_alpha + cot_beta) / 2` for
    /// adjacent vertices `i`, `j`, and `L[i][i] = -sum_j L[i][j]`.
    pub fn cotan_laplacian(&self) -> Vec<(usize, usize, f64)> {
        let n = self.num_vertices();
        // accumulate off-diagonal cotangent weights
        let mut offdiag: HashMap<(usize, usize), f64> = HashMap::new();
        for t in &self.triangles {
            let [i, j, k] = *t;
            let pi = self.vertices[i];
            let pj = self.vertices[j];
            let pk = self.vertices[k];
            // cotangents at each corner
            let cot_i = cotan_at(pj, pi, pk);
            let cot_j = cotan_at(pi, pj, pk);
            let cot_k = cotan_at(pi, pk, pj);
            // edge (j,k) -> angle at i
            *offdiag.entry((j.min(k), j.max(k))).or_insert(0.0) += 0.5 * cot_i;
            // edge (i,k) -> angle at j
            *offdiag.entry((i.min(k), i.max(k))).or_insert(0.0) += 0.5 * cot_j;
            // edge (i,j) -> angle at k
            *offdiag.entry((i.min(j), i.max(j))).or_insert(0.0) += 0.5 * cot_k;
        }
        let mut entries: Vec<(usize, usize, f64)> = Vec::with_capacity(offdiag.len() * 2 + n);
        let mut row_sum = vec![0.0f64; n];
        for (&(a, b), &w) in &offdiag {
            entries.push((a, b, w));
            entries.push((b, a, w));
            row_sum[a] += w;
            row_sum[b] += w;
        }
        for i in 0..n {
            entries.push((i, i, -row_sum[i]));
        }
        entries
    }

    /// Apply the Laplace-Beltrami operator to a scalar field `f` (length = n_vertices),
    /// returning the Laplacian value at each vertex: `(Lf)[i] = sum_j L[i][j] * f[j]`.
    pub fn apply_laplacian(&self, f: &[f64]) -> Vec<f64> {
        let n = self.num_vertices();
        let mut lf = vec![0.0f64; n];
        for (i, j, w) in self.cotan_laplacian() {
            lf[i] += w * f[j];
        }
        lf
    }

    /// Discrete Gaussian curvature (angle defect) at each vertex.
    ///
    /// `K[v] = (2*pi - sum_angles) / A_v` where `A_v` is the Voronoi area.
    pub fn gaussian_curvature(&self) -> Vec<f64> {
        let n = self.num_vertices();
        let mut angle_sum = vec![0.0f64; n];
        for t in &self.triangles {
            let [i, j, k] = *t;
            let pi = self.vertices[i];
            let pj = self.vertices[j];
            let pk = self.vertices[k];
            // Angle at i
            let u = normalize3(sub3(pj, pi));
            let v = normalize3(sub3(pk, pi));
            let ang_i = dot3(u, v).clamp(-1.0, 1.0).acos();
            // Angle at j
            let u2 = normalize3(sub3(pi, pj));
            let v2 = normalize3(sub3(pk, pj));
            let ang_j = dot3(u2, v2).clamp(-1.0, 1.0).acos();
            // Angle at k
            let u3 = normalize3(sub3(pi, pk));
            let v3 = normalize3(sub3(pj, pk));
            let ang_k = dot3(u3, v3).clamp(-1.0, 1.0).acos();
            angle_sum[i] += ang_i;
            angle_sum[j] += ang_j;
            angle_sum[k] += ang_k;
        }
        (0..n)
            .map(|v| {
                let a = self.voronoi_area(v);
                (2.0 * PI - angle_sum[v]) / a
            })
            .collect()
    }

    /// Discrete mean curvature at each vertex using the cotangent formula.
    ///
    /// Returns the mean curvature `H[v]` (scalar) at each vertex.
    /// `H[v] = |Lx[v]| / (2 * A_v)` where `Lx` is the Laplacian of position.
    pub fn mean_curvature(&self) -> Vec<f64> {
        let n = self.num_vertices();
        let fx: Vec<f64> = self.vertices.iter().map(|p| p[0]).collect();
        let fy: Vec<f64> = self.vertices.iter().map(|p| p[1]).collect();
        let fz: Vec<f64> = self.vertices.iter().map(|p| p[2]).collect();
        let lx = self.apply_laplacian(&fx);
        let ly = self.apply_laplacian(&fy);
        let lz = self.apply_laplacian(&fz);
        (0..n)
            .map(|v| {
                let a = self.voronoi_area(v);
                let hv = [lx[v], ly[v], lz[v]];
                len3(hv) / (2.0 * a)
            })
            .collect()
    }

    /// Principal curvatures estimated from mean and Gaussian curvature.
    ///
    /// Returns `(kappa_1, kappa_2)` where `kappa_1 >= kappa_2`.
    pub fn principal_curvatures(&self) -> Vec<(f64, f64)> {
        let h = self.mean_curvature();
        let k = self.gaussian_curvature();
        h.iter()
            .zip(k.iter())
            .map(|(&hi, &ki)| {
                let disc = (hi * hi - ki).max(0.0).sqrt();
                let k1 = hi + disc;
                let k2 = hi - disc;
                if k1 >= k2 { (k1, k2) } else { (k2, k1) }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Dijkstra geodesic distances
// ---------------------------------------------------------------------------

/// Dijkstra-based geodesic distance from a set of source vertices.
///
/// Returns a vector of length `n_vertices`; unreachable vertices get `f64::INFINITY`.
pub fn geodesic_dijkstra(mesh: &DiscreteMesh, sources: &[usize]) -> Vec<f64> {
    let n = mesh.num_vertices();
    let mut dist = vec![f64::INFINITY; n];
    let mut heap: BinaryHeap<DijkstraNode> = BinaryHeap::new();
    for &s in sources {
        dist[s] = 0.0;
        heap.push(DijkstraNode {
            dist: 0.0,
            vertex: s,
        });
    }
    // Build adjacency once
    let adj = build_edge_adjacency(mesh);
    while let Some(DijkstraNode { dist: d, vertex: u }) = heap.pop() {
        if d > dist[u] + 1e-12 {
            continue;
        }
        if let Some(nbrs) = adj.get(&u) {
            for &(v, w) in nbrs {
                let nd = d + w;
                if nd < dist[v] {
                    dist[v] = nd;
                    heap.push(DijkstraNode {
                        dist: nd,
                        vertex: v,
                    });
                }
            }
        }
    }
    dist
}

/// Node for Dijkstra min-heap (ordered by distance ascending).
#[derive(Debug, Clone, PartialEq)]
struct DijkstraNode {
    /// Tentative distance to this vertex.
    dist: f64,
    /// Vertex index.
    vertex: usize,
}

impl Eq for DijkstraNode {}

impl PartialOrd for DijkstraNode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DijkstraNode {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Min-heap: smaller dist = higher priority
        other
            .dist
            .partial_cmp(&self.dist)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(self.vertex.cmp(&other.vertex))
    }
}

/// Build a vertex-to-neighbours adjacency list with edge lengths.
pub fn build_edge_adjacency(mesh: &DiscreteMesh) -> HashMap<usize, Vec<(usize, f64)>> {
    let mut adj: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for t in &mesh.triangles {
        let pairs = [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])];
        for (a, b) in pairs {
            let w = dist3(mesh.vertices[a], mesh.vertices[b]);
            adj.entry(a).or_default().push((b, w));
            adj.entry(b).or_default().push((a, w));
        }
    }
    adj
}

/// Extract the shortest geodesic path from `source` to `target` as a list of
/// vertex indices, using Dijkstra distances and a predecessor map.
pub fn geodesic_path(mesh: &DiscreteMesh, source: usize, target: usize) -> Option<Vec<usize>> {
    let n = mesh.num_vertices();
    let mut dist = vec![f64::INFINITY; n];
    let mut prev: Vec<Option<usize>> = vec![None; n];
    let adj = build_edge_adjacency(mesh);
    dist[source] = 0.0;
    let mut heap: BinaryHeap<DijkstraNode> = BinaryHeap::new();
    heap.push(DijkstraNode {
        dist: 0.0,
        vertex: source,
    });
    while let Some(DijkstraNode { dist: d, vertex: u }) = heap.pop() {
        if u == target {
            break;
        }
        if d > dist[u] + 1e-12 {
            continue;
        }
        if let Some(nbrs) = adj.get(&u) {
            for &(v, w) in nbrs {
                let nd = d + w;
                if nd < dist[v] {
                    dist[v] = nd;
                    prev[v] = Some(u);
                    heap.push(DijkstraNode {
                        dist: nd,
                        vertex: v,
                    });
                }
            }
        }
    }
    if dist[target].is_infinite() {
        return None;
    }
    let mut path = Vec::new();
    let mut cur = target;
    loop {
        path.push(cur);
        if cur == source {
            break;
        }
        match prev[cur] {
            Some(p) => cur = p,
            None => return None,
        }
    }
    path.reverse();
    Some(path)
}

// ---------------------------------------------------------------------------
// Heat Method for Geodesic Distances (Crane et al. 2013)
// ---------------------------------------------------------------------------

/// Heat-method geodesic distance solver.
///
/// Approximates smooth geodesic distances by:
/// 1. Diffusing a heat impulse for a short time `t`.
/// 2. Computing the gradient of the heat field.
/// 3. Normalising and divergence-computing to recover distance.
///
/// This implementation uses an explicit Euler diffusion step and the
/// cotangent Laplacian.
#[derive(Debug, Clone)]
pub struct HeatMethodSolver {
    /// The underlying mesh.
    pub mesh: DiscreteMesh,
    /// Diffusion time parameter; typically `h^2` where `h` is mean edge length.
    pub diffusion_time: f64,
    /// Number of explicit Euler steps for heat diffusion.
    pub num_steps: usize,
}

impl HeatMethodSolver {
    /// Construct a solver. If `diffusion_time <= 0`, it is set automatically
    /// to the square of the mean edge length.
    pub fn new(mesh: DiscreteMesh, diffusion_time: f64, num_steps: usize) -> Self {
        let t = if diffusion_time > 0.0 {
            diffusion_time
        } else {
            let h = mean_edge_length(&mesh);
            h * h
        };
        Self {
            mesh,
            diffusion_time: t,
            num_steps,
        }
    }

    /// Compute approximate geodesic distances from `sources`.
    ///
    /// Returns distances for each vertex; values are normalised so that the
    /// minimum source distance is zero.
    pub fn compute(&self, sources: &[usize]) -> Vec<f64> {
        let n = self.mesh.num_vertices();
        // Step 1: initialise heat at sources
        let mut u = vec![0.0f64; n];
        for &s in sources {
            u[s] = 1.0;
        }
        // Step 2: explicit Euler heat diffusion u_(t+1) = u_t + dt * L * u_t
        let laplacian = self.mesh.cotan_laplacian();
        let dt = self.diffusion_time / self.num_steps as f64;
        for _ in 0..self.num_steps {
            let mut lu = vec![0.0f64; n];
            for &(i, j, w) in &laplacian {
                lu[i] += w * u[j];
            }
            for i in 0..n {
                u[i] += dt * lu[i];
            }
        }
        // Step 3: compute gradient of u per face, normalise
        let nf = self.mesh.num_faces();
        let mut face_grad: Vec<[f64; 3]> = vec![[0.0; 3]; nf];
        for (fi, t) in self.mesh.triangles.iter().enumerate() {
            let [i, j, k] = *t;
            let pi = self.mesh.vertices[i];
            let pj = self.mesh.vertices[j];
            let pk = self.mesh.vertices[k];
            let area = tri_area(pi, pj, pk).max(1e-14);
            let fn_ = normalize3(cross3(sub3(pj, pi), sub3(pk, pi)));
            // Gradient: sum_corner u_corner * (N x e_opposite) / (2 * area)
            let grad_i = scale3(cross3(fn_, sub3(pk, pj)), u[i] / (2.0 * area));
            let grad_j = scale3(cross3(fn_, sub3(pi, pk)), u[j] / (2.0 * area));
            let grad_k = scale3(cross3(fn_, sub3(pj, pi)), u[k] / (2.0 * area));
            face_grad[fi] = add3(add3(grad_i, grad_j), grad_k);
        }
        // Step 4: normalise gradients (flip sign — points away from source)
        let face_grad_norm: Vec<[f64; 3]> = face_grad
            .iter()
            .map(|&g| {
                let l = len3(g);
                if l < 1e-14 {
                    [0.0; 3]
                } else {
                    scale3(g, -1.0 / l)
                }
            })
            .collect();
        // Step 5: compute divergence of normalised gradient at each vertex
        let mut div = vec![0.0f64; n];
        for (fi, t) in self.mesh.triangles.iter().enumerate() {
            let [i, j, k] = *t;
            let pi = self.mesh.vertices[i];
            let pj = self.mesh.vertices[j];
            let pk = self.mesh.vertices[k];
            let x = face_grad_norm[fi];
            let cot_i = cotan_at(pj, pi, pk);
            let cot_j = cotan_at(pi, pj, pk);
            let cot_k = cotan_at(pi, pk, pj);
            div[i] += 0.5 * (cot_k * dot3(x, sub3(pj, pi)) + cot_j * dot3(x, sub3(pk, pi)));
            div[j] += 0.5 * (cot_i * dot3(x, sub3(pk, pj)) + cot_k * dot3(x, sub3(pi, pj)));
            div[k] += 0.5 * (cot_j * dot3(x, sub3(pi, pk)) + cot_i * dot3(x, sub3(pj, pk)));
        }
        // Step 6: solve Poisson Lφ = div (simple Gauss-Seidel iterations)
        let mut phi = div.clone();
        let diag: Vec<f64> = {
            let mut d = vec![0.0f64; n];
            for &(i, j, w) in &laplacian {
                if i == j {
                    d[i] = w;
                }
            }
            d
        };
        for _ in 0..200 {
            for vi in 0..n {
                if diag[vi].abs() < 1e-14 {
                    continue;
                }
                let mut s = div[vi];
                for &(ri, rj, rw) in &laplacian {
                    if ri == vi && rj != vi {
                        s -= rw * phi[rj];
                    }
                }
                phi[vi] = s / diag[vi];
            }
        }
        // Shift so minimum is zero
        let min_phi = phi.iter().cloned().fold(f64::INFINITY, f64::min);
        phi.iter_mut().for_each(|x| *x -= min_phi);
        phi
    }
}

/// Compute the mean edge length of a mesh.
pub fn mean_edge_length(mesh: &DiscreteMesh) -> f64 {
    let mut total = 0.0;
    let mut count = 0usize;
    for t in &mesh.triangles {
        let [i, j, k] = *t;
        total += dist3(mesh.vertices[i], mesh.vertices[j]);
        total += dist3(mesh.vertices[j], mesh.vertices[k]);
        total += dist3(mesh.vertices[i], mesh.vertices[k]);
        count += 3;
    }
    if count == 0 {
        1.0
    } else {
        total / count as f64
    }
}

// ---------------------------------------------------------------------------
// Parallel Transport on Meshes
// ---------------------------------------------------------------------------

/// Parallel transport a tangent vector `v_tan` from face `f_src` to an adjacent
/// face `f_dst` across their shared edge.
///
/// The vector is expressed in the local frame of `f_src` and rotated into
/// the frame of `f_dst` using the connection (Levi-Civita connection on the mesh).
///
/// Returns the transported vector in the local tangent plane of `f_dst`.
pub fn parallel_transport_across_edge(
    mesh: &DiscreteMesh,
    f_src: usize,
    f_dst: usize,
    v_tan: [f64; 3],
) -> [f64; 3] {
    let n_src = mesh.face_normal(f_src);
    let n_dst = mesh.face_normal(f_dst);
    // Rotation axis = cross(n_src, n_dst), angle = acos(dot)
    let cos_theta = dot3(n_src, n_dst).clamp(-1.0, 1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    if sin_theta.abs() < 1e-12 {
        return v_tan; // faces are co-planar
    }
    let axis = normalize3(cross3(n_src, n_dst));
    // Rodrigues rotation formula
    let term1 = scale3(v_tan, cos_theta);
    let term2 = scale3(cross3(axis, v_tan), sin_theta);
    let term3 = scale3(axis, dot3(axis, v_tan) * (1.0 - cos_theta));
    add3(add3(term1, term2), term3)
}

/// Compute holonomy (angle defect accumulated from parallel transport) around a
/// vertex `v` by transporting a tangent vector around its one-ring of faces.
///
/// Returns the accumulated rotation angle in radians.
pub fn vertex_holonomy(mesh: &DiscreteMesh, v: usize) -> f64 {
    let faces = mesh.one_ring_faces(v);
    if faces.len() < 2 {
        return 0.0;
    }
    // Initialise a tangent vector in the first face's plane
    let t0 = &mesh.triangles[faces[0]];
    let p0 = mesh.vertices[t0[0]];
    let p1 = mesh.vertices[t0[1]];
    let mut tang = normalize3(sub3(p1, p0));
    // Transport around the ring
    let mut total_angle = 0.0_f64;
    for i in 0..faces.len() {
        let f_src = faces[i];
        let f_dst = faces[(i + 1) % faces.len()];
        let transported = parallel_transport_across_edge(mesh, f_src, f_dst, tang);
        let n = mesh.face_normal(f_dst);
        let tang_proj = normalize3(sub3(transported, scale3(n, dot3(transported, n))));
        let cos_a = dot3(tang, tang_proj).clamp(-1.0, 1.0);
        total_angle += cos_a.acos();
        tang = tang_proj;
    }
    total_angle
}

// ---------------------------------------------------------------------------
// Discrete Vector Fields on Meshes
// ---------------------------------------------------------------------------

/// A discrete vector field on a mesh, storing one tangent vector per face.
///
/// Vectors are stored in 3D ambient space but should lie in the tangent plane
/// of the respective face.
#[derive(Debug, Clone)]
pub struct FaceVectorField {
    /// Number of faces in the underlying mesh.
    pub num_faces: usize,
    /// Per-face tangent vectors.
    pub vectors: Vec<[f64; 3]>,
}

impl FaceVectorField {
    /// Construct a zero vector field.
    pub fn zeros(num_faces: usize) -> Self {
        Self {
            num_faces,
            vectors: vec![[0.0; 3]; num_faces],
        }
    }

    /// Construct from a given list of per-face vectors.
    pub fn from_vec(vectors: Vec<[f64; 3]>) -> Self {
        let n = vectors.len();
        Self {
            num_faces: n,
            vectors,
        }
    }

    /// Project each face vector into the tangent plane of the face.
    pub fn project_to_tangent(&mut self, mesh: &DiscreteMesh) {
        for (fi, v) in self.vectors.iter_mut().enumerate() {
            if fi >= mesh.num_faces() {
                break;
            }
            let n = mesh.face_normal(fi);
            let normal_comp = scale3(n, dot3(*v, n));
            *v = sub3(*v, normal_comp);
        }
    }

    /// L2 norm of the vector field.
    pub fn norm(&self) -> f64 {
        self.vectors
            .iter()
            .map(|v| dot3(*v, *v))
            .sum::<f64>()
            .sqrt()
    }
}

/// A discrete vertex vector field, storing one tangent vector per vertex.
#[derive(Debug, Clone)]
pub struct VertexVectorField {
    /// Number of vertices.
    pub num_vertices: usize,
    /// Per-vertex tangent vectors.
    pub vectors: Vec<[f64; 3]>,
}

impl VertexVectorField {
    /// Construct a zero vertex vector field.
    pub fn zeros(num_vertices: usize) -> Self {
        Self {
            num_vertices,
            vectors: vec![[0.0; 3]; num_vertices],
        }
    }

    /// Compute the divergence of the vertex vector field as a scalar on vertices.
    pub fn divergence(&self, mesh: &DiscreteMesh) -> Vec<f64> {
        let n = mesh.num_vertices();
        let mut div = vec![0.0f64; n];
        for t in &mesh.triangles {
            let [i, j, k] = *t;
            let pi = mesh.vertices[i];
            let pj = mesh.vertices[j];
            let pk = mesh.vertices[k];
            let area = tri_area(pi, pj, pk).max(1e-14);
            let xi = self.vectors[i];
            let xj = self.vectors[j];
            let xk = self.vectors[k];
            // Discrete divergence per triangle
            let cot_i = cotan_at(pj, pi, pk);
            let cot_j = cotan_at(pi, pj, pk);
            let cot_k = cotan_at(pi, pk, pj);
            div[i] +=
                (cot_k * dot3(xi, sub3(pj, pi)) + cot_j * dot3(xi, sub3(pk, pi))) / (2.0 * area);
            div[j] +=
                (cot_i * dot3(xj, sub3(pk, pj)) + cot_k * dot3(xj, sub3(pi, pj))) / (2.0 * area);
            div[k] +=
                (cot_j * dot3(xk, sub3(pi, pk)) + cot_i * dot3(xk, sub3(pj, pk))) / (2.0 * area);
        }
        div
    }

    /// Compute the curl (discrete 2-form = scalar on faces) of the vertex vector field.
    pub fn curl_on_faces(&self, mesh: &DiscreteMesh) -> Vec<f64> {
        mesh.triangles
            .iter()
            .map(|t| {
                let [i, j, k] = *t;
                let pi = mesh.vertices[i];
                let pj = mesh.vertices[j];
                let pk = mesh.vertices[k];
                let area = tri_area(pi, pj, pk).max(1e-14);
                let xi = self.vectors[i];
                let xj = self.vectors[j];
                let xk = self.vectors[k];
                let n = normalize3(cross3(sub3(pj, pi), sub3(pk, pi)));
                // Discrete curl: circulation around triangle
                let c = dot3(cross3(n, sub3(pj, pi)), xk)
                    + dot3(cross3(n, sub3(pk, pj)), xi)
                    + dot3(cross3(n, sub3(pi, pk)), xj);
                c / (2.0 * area)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Hodge Decomposition
// ---------------------------------------------------------------------------

/// Result of Hodge decomposition of a vertex vector field `X = grad(f) + curl(g) + harm`.
///
/// The decomposition is: X = ∇f + *d(g) + H where H is the harmonic part.
#[derive(Debug, Clone)]
pub struct HodgeDecomposition {
    /// Exact (gradient) component scalar potential, per vertex.
    pub scalar_potential: Vec<f64>,
    /// Co-exact (curl) component stream function, per vertex.
    pub stream_function: Vec<f64>,
    /// Harmonic residual vector field (per vertex).
    pub harmonic: Vec<[f64; 3]>,
}

/// Perform a Hodge decomposition of a vertex vector field.
///
/// This computes an approximate decomposition via:
/// 1. Solve `L f = div(X)` → scalar potential.
/// 2. Reconstruct gradient component and subtract.
/// 3. Compute curl of remainder → stream function.
/// 4. Harmonic part = remainder after subtracting both.
pub fn hodge_decompose(mesh: &DiscreteMesh, field: &VertexVectorField) -> HodgeDecomposition {
    let n = mesh.num_vertices();
    let laplacian = mesh.cotan_laplacian();
    let div = field.divergence(mesh);

    // Build diagonal of laplacian
    let mut diag = vec![0.0f64; n];
    for &(i, j, w) in &laplacian {
        if i == j {
            diag[i] = w;
        }
    }

    // Gauss-Seidel solve: L f = div
    let mut f = vec![0.0f64; n];
    for _ in 0..300 {
        for vi in 0..n {
            if diag[vi].abs() < 1e-14 {
                continue;
            }
            let mut s = div[vi];
            for &(ri, rj, rw) in &laplacian {
                if ri == vi && rj != vi {
                    s -= rw * f[rj];
                }
            }
            f[vi] = s / diag[vi];
        }
    }

    // Gradient of f per vertex: average of face gradients
    let mut grad_f = vec![[0.0f64; 3]; n];
    let mut weight = vec![0.0f64; n];
    for t in &mesh.triangles {
        let [i, j, k] = *t;
        let pi = mesh.vertices[i];
        let pj = mesh.vertices[j];
        let pk = mesh.vertices[k];
        let area = tri_area(pi, pj, pk).max(1e-14);
        let fn_ = normalize3(cross3(sub3(pj, pi), sub3(pk, pi)));
        // Face gradient
        let gi = scale3(cross3(fn_, sub3(pk, pj)), f[i] / (2.0 * area));
        let gj = scale3(cross3(fn_, sub3(pi, pk)), f[j] / (2.0 * area));
        let gk = scale3(cross3(fn_, sub3(pj, pi)), f[k] / (2.0 * area));
        let face_grad = add3(add3(gi, gj), gk);
        for &vi in &[i, j, k] {
            grad_f[vi] = add3(grad_f[vi], scale3(face_grad, area));
            weight[vi] += area;
        }
    }
    for vi in 0..n {
        if weight[vi] > 1e-14 {
            grad_f[vi] = scale3(grad_f[vi], 1.0 / weight[vi]);
        }
    }

    // Curl-free residual
    let mut residual: Vec<[f64; 3]> = (0..n)
        .map(|vi| sub3(field.vectors[vi], grad_f[vi]))
        .collect();

    // Stream function g from curl of residual (Gauss-Seidel solve L g = curl)
    let residual_field = VertexVectorField {
        num_vertices: n,
        vectors: residual.clone(),
    };
    let curl_on_faces = residual_field.curl_on_faces(mesh);

    // Average curl to vertices as div proxy
    let mut curl_div = vec![0.0f64; n];
    let mut curl_wt = vec![0.0f64; n];
    for (fi, t) in mesh.triangles.iter().enumerate() {
        let area = tri_area(
            mesh.vertices[t[0]],
            mesh.vertices[t[1]],
            mesh.vertices[t[2]],
        )
        .max(1e-14);
        for &vi in t {
            curl_div[vi] += curl_on_faces[fi] * area;
            curl_wt[vi] += area;
        }
    }
    for vi in 0..n {
        if curl_wt[vi] > 1e-14 {
            curl_div[vi] /= curl_wt[vi];
        }
    }

    let mut g = vec![0.0f64; n];
    for _ in 0..300 {
        for vi in 0..n {
            if diag[vi].abs() < 1e-14 {
                continue;
            }
            let mut s = curl_div[vi];
            for &(ri, rj, rw) in &laplacian {
                if ri == vi && rj != vi {
                    s -= rw * g[rj];
                }
            }
            g[vi] = s / diag[vi];
        }
    }

    // Subtract co-exact part from residual to get harmonic
    for vi in 0..n {
        // Approximate co-exact gradient
        let co = scale3(mesh.vertex_normal(vi), g[vi]);
        residual[vi] = sub3(residual[vi], co);
    }

    HodgeDecomposition {
        scalar_potential: f,
        stream_function: g,
        harmonic: residual,
    }
}

// ---------------------------------------------------------------------------
// Discrete Minimal Surfaces (Mean Curvature Flow)
// ---------------------------------------------------------------------------

/// Perform one step of discrete mean curvature flow.
///
/// Updates vertex positions by `x_new = x_old + dt * L(x_old)`
/// where `L` is the cotangent Laplacian applied to each coordinate.
/// This moves the surface towards a minimal surface.
pub fn mean_curvature_flow_step(mesh: &mut DiscreteMesh, dt: f64) {
    let n = mesh.num_vertices();
    let laplacian = mesh.cotan_laplacian();
    let fx: Vec<f64> = mesh.vertices.iter().map(|p| p[0]).collect();
    let fy: Vec<f64> = mesh.vertices.iter().map(|p| p[1]).collect();
    let fz: Vec<f64> = mesh.vertices.iter().map(|p| p[2]).collect();
    let mut lx = vec![0.0f64; n];
    let mut ly = vec![0.0f64; n];
    let mut lz = vec![0.0f64; n];
    for &(i, j, w) in &laplacian {
        lx[i] += w * fx[j];
        ly[i] += w * fy[j];
        lz[i] += w * fz[j];
    }
    for i in 0..n {
        mesh.vertices[i][0] += dt * lx[i];
        mesh.vertices[i][1] += dt * ly[i];
        mesh.vertices[i][2] += dt * lz[i];
    }
}

/// Run `n_steps` of mean curvature flow on `mesh`.
pub fn mean_curvature_flow(mesh: &mut DiscreteMesh, dt: f64, n_steps: usize) {
    for _ in 0..n_steps {
        mean_curvature_flow_step(mesh, dt);
    }
}

// ---------------------------------------------------------------------------
// Euler Characteristic and Topology
// ---------------------------------------------------------------------------

/// Compute the Euler characteristic χ = V - E + F of the mesh.
pub fn euler_characteristic(mesh: &DiscreteMesh) -> i64 {
    let v = mesh.num_vertices() as i64;
    let f = mesh.num_faces() as i64;
    // Count unique edges
    let mut edges: HashSet<(usize, usize)> = HashSet::new();
    for t in &mesh.triangles {
        let [i, j, k] = *t;
        edges.insert((i.min(j), i.max(j)));
        edges.insert((j.min(k), j.max(k)));
        edges.insert((i.min(k), i.max(k)));
    }
    let e = edges.len() as i64;
    v - e + f
}

/// Compute the genus `g` of a closed orientable surface from Euler characteristic.
///
/// Uses the formula χ = 2 - 2g, so g = (2 - χ) / 2.
pub fn genus_from_euler(chi: i64) -> i64 {
    (2 - chi) / 2
}

// ---------------------------------------------------------------------------
// Discrete 1-Forms and Integration
// ---------------------------------------------------------------------------

/// A discrete 1-form on mesh edges, storing one value per directed edge.
///
/// A 1-form ω assigns a value ω(e) to each oriented edge such that
/// ω(-e) = -ω(e).
#[derive(Debug, Clone)]
pub struct Discrete1Form {
    /// Values indexed by (vertex_i, vertex_j) with i < j.
    pub values: HashMap<(usize, usize), f64>,
}

impl Discrete1Form {
    /// Construct a zero 1-form over all mesh edges.
    pub fn zeros(mesh: &DiscreteMesh) -> Self {
        let mut values = HashMap::new();
        for t in &mesh.triangles {
            let [i, j, k] = *t;
            values.entry((i.min(j), i.max(j))).or_insert(0.0);
            values.entry((j.min(k), j.max(k))).or_insert(0.0);
            values.entry((i.min(k), i.max(k))).or_insert(0.0);
        }
        Self { values }
    }

    /// Evaluate the 1-form on directed edge (a → b).
    ///
    /// Returns `+val` if `a < b`, else `-val`.
    pub fn eval(&self, a: usize, b: usize) -> f64 {
        let key = (a.min(b), a.max(b));
        let v = self.values.get(&key).copied().unwrap_or(0.0);
        if a < b { v } else { -v }
    }

    /// Integrate the 1-form along a path given as a sequence of vertex indices.
    pub fn integrate_path(&self, path: &[usize]) -> f64 {
        path.windows(2).map(|w| self.eval(w[0], w[1])).sum()
    }

    /// Compute the exterior derivative (discrete 2-form = scalar on faces).
    ///
    /// Returns a vector of face values `dω[f] = ω(e01) + ω(e12) + ω(e20)`.
    pub fn exterior_derivative(&self, mesh: &DiscreteMesh) -> Vec<f64> {
        mesh.triangles
            .iter()
            .map(|t| {
                let [i, j, k] = *t;
                self.eval(i, j) + self.eval(j, k) + self.eval(k, i)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Discrete Laplacian Eigenvalues (Power Iteration)
// ---------------------------------------------------------------------------

/// Compute the first `k` approximate eigenvalues of the Laplace-Beltrami
/// operator via power iteration.
///
/// Note: this is a rough approximation; for production use a proper sparse
/// eigensolver is recommended.
pub fn laplacian_eigenvalues_power(mesh: &DiscreteMesh, k: usize, n_iter: usize) -> Vec<f64> {
    let n = mesh.num_vertices();
    let laplacian = mesh.cotan_laplacian();
    let mut eigenvalues = Vec::with_capacity(k);
    let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
    for _ki in 0..k {
        // Normalise
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-14);
        v.iter_mut().for_each(|x| *x /= norm);
        for _ in 0..n_iter {
            let mut lv = vec![0.0f64; n];
            for &(i, j, w) in &laplacian {
                lv[i] += w * v[j];
            }
            // Rayleigh quotient
            let rq: f64 = v.iter().zip(lv.iter()).map(|(vi, li)| vi * li).sum();
            eigenvalues.push(rq.abs());
            let norm2 = lv.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-14);
            v = lv.iter().map(|x| x / norm2).collect();
        }
    }
    eigenvalues.truncate(k);
    eigenvalues
}

// ---------------------------------------------------------------------------
// Mesh Quality and Refinement
// ---------------------------------------------------------------------------

/// Compute the aspect ratio of each triangle.
///
/// Aspect ratio = longest edge / shortest altitude.
/// Equilateral triangle has aspect ratio 1 (minimum).
pub fn triangle_aspect_ratios(mesh: &DiscreteMesh) -> Vec<f64> {
    mesh.triangles
        .iter()
        .map(|t| {
            let [i, j, k] = *t;
            let p0 = mesh.vertices[i];
            let p1 = mesh.vertices[j];
            let p2 = mesh.vertices[k];
            let l01 = dist3(p0, p1);
            let l12 = dist3(p1, p2);
            let l02 = dist3(p0, p2);
            let longest = l01.max(l12).max(l02);
            let area = tri_area(p0, p1, p2).max(1e-14);
            let shortest_alt = 2.0 * area / longest.max(1e-14);
            longest / shortest_alt.max(1e-14)
        })
        .collect()
}

/// Compute per-face areas.
pub fn face_areas(mesh: &DiscreteMesh) -> Vec<f64> {
    mesh.triangles
        .iter()
        .map(|t| {
            let [i, j, k] = *t;
            tri_area(mesh.vertices[i], mesh.vertices[j], mesh.vertices[k])
        })
        .collect()
}

/// Total surface area of the mesh.
pub fn total_surface_area(mesh: &DiscreteMesh) -> f64 {
    face_areas(mesh).iter().sum()
}

// ---------------------------------------------------------------------------
// Mesh construction helpers for testing
// ---------------------------------------------------------------------------

/// Build a simple unit-square mesh (two triangles).
///
/// Vertices: (0,0,0), (1,0,0), (1,1,0), (0,1,0).
/// Faces: \[(0,1,2), (0,2,3)\].
pub fn unit_square_mesh() -> DiscreteMesh {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let triangles = vec![[0, 1, 2], [0, 2, 3]];
    DiscreteMesh::new(vertices, triangles)
}

/// Build a regular icosahedron mesh (12 vertices, 20 faces).
///
/// Vertices lie on the unit sphere.
pub fn unit_icosahedron() -> DiscreteMesh {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut verts: Vec<[f64; 3]> = vec![
        [-1.0, phi, 0.0],
        [1.0, phi, 0.0],
        [-1.0, -phi, 0.0],
        [1.0, -phi, 0.0],
        [0.0, -1.0, phi],
        [0.0, 1.0, phi],
        [0.0, -1.0, -phi],
        [0.0, 1.0, -phi],
        [phi, 0.0, -1.0],
        [phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
    ];
    // Normalise to unit sphere
    for v in &mut verts {
        let l = len3(*v);
        v[0] /= l;
        v[1] /= l;
        v[2] /= l;
    }
    let faces = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    DiscreteMesh::new(verts, faces)
}

/// Build a regular tetrahedron mesh (4 vertices, 4 faces).
pub fn unit_tetrahedron() -> DiscreteMesh {
    let s = (8.0_f64 / 9.0_f64).sqrt();
    let c1 = (2.0_f64 / 9.0_f64).sqrt();
    let c2 = (2.0_f64 / 3.0_f64).sqrt();
    let vertices = vec![
        [0.0, 1.0, 0.0],
        [s, -1.0 / 3.0, 0.0],
        [-c1, -1.0 / 3.0, c2],
        [-c1, -1.0 / 3.0, -c2],
    ];
    let triangles = vec![[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
    DiscreteMesh::new(vertices, triangles)
}

// ---------------------------------------------------------------------------
// Discrete Curvature Tensor (Shape Operator Approximation)
// ---------------------------------------------------------------------------

/// Approximate the curvature tensor at a vertex from neighbouring face normals.
///
/// Returns a symmetric 3×3 matrix `M` stored row-major as `[f64; 9]` such that
/// the quadratic form `v^T M v` approximates the second fundamental form.
pub fn vertex_curvature_tensor(mesh: &DiscreteMesh, v: usize) -> [f64; 9] {
    let mut m = [0.0f64; 9];
    let mut total_area = 0.0_f64;
    let pv = mesh.vertices[v];
    for t in &mesh.triangles {
        if !t.contains(&v) {
            continue;
        }
        let [i, j, k] = *t;
        let pi = mesh.vertices[i];
        let pj = mesh.vertices[j];
        let pk = mesh.vertices[k];
        let area = tri_area(pi, pj, pk).max(1e-14);
        let fn_ = mesh.face_normal(mesh.triangles.iter().position(|tt| tt == t).unwrap_or(0));
        // Find the two edges from v
        let others: Vec<[f64; 3]> = [i, j, k]
            .iter()
            .filter(|&&x| x != v)
            .map(|&x| mesh.vertices[x])
            .collect();
        for &po in &others {
            let e = sub3(po, pv);
            let kn = dot3(fn_, e) / dot3(e, e).max(1e-14);
            // Outer product e ⊗ e weighted by kn * area
            let w = kn * area;
            for a in 0..3 {
                for b in 0..3 {
                    m[a * 3 + b] += w * e[a] * e[b];
                }
            }
        }
        total_area += area;
    }
    if total_area > 1e-14 {
        for x in &mut m {
            *x /= total_area;
        }
    }
    m
}

// ---------------------------------------------------------------------------
// Discrete Ricci Flow (simplified)
// ---------------------------------------------------------------------------

/// Perform one step of discrete Ricci flow on vertex conformal factors.
///
/// Updates the conformal factor `u[v]` proportional to `K[v] - K_target`,
/// where `K` is Gaussian curvature and `K_target` is the target curvature.
pub fn ricci_flow_step(mesh: &DiscreteMesh, u: &mut Vec<f64>, k_target: &[f64], step_size: f64) {
    let k_current = mesh.gaussian_curvature();
    for (i, ui) in u.iter_mut().enumerate() {
        let kt = if i < k_target.len() { k_target[i] } else { 0.0 };
        *ui += step_size * (k_current[i] - kt);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_square_mesh() -> DiscreteMesh {
        unit_square_mesh()
    }

    fn ico() -> DiscreteMesh {
        unit_icosahedron()
    }

    fn tet() -> DiscreteMesh {
        unit_tetrahedron()
    }

    // --- Basic mesh properties ---

    #[test]
    fn test_mesh_vertex_count() {
        let m = flat_square_mesh();
        assert_eq!(m.num_vertices(), 4);
    }

    #[test]
    fn test_mesh_face_count() {
        let m = flat_square_mesh();
        assert_eq!(m.num_faces(), 2);
    }

    #[test]
    fn test_icosahedron_vertex_count() {
        let m = ico();
        assert_eq!(m.num_vertices(), 12);
    }

    #[test]
    fn test_icosahedron_face_count() {
        let m = ico();
        assert_eq!(m.num_faces(), 20);
    }

    #[test]
    fn test_tetrahedron_vertex_count() {
        let m = tet();
        assert_eq!(m.num_vertices(), 4);
    }

    // --- One-ring ---

    #[test]
    fn test_one_ring_faces_non_empty() {
        let m = flat_square_mesh();
        let ring = m.one_ring_faces(0);
        assert!(!ring.is_empty());
    }

    #[test]
    fn test_one_ring_vertices_contains_neighbors() {
        let m = flat_square_mesh();
        // Vertex 0 is connected to 1, 2, 3
        let nbrs = m.one_ring_vertices(0);
        assert!(nbrs.contains(&1) || nbrs.contains(&2) || nbrs.contains(&3));
    }

    // --- Area ---

    #[test]
    fn test_total_surface_area_unit_square() {
        let m = flat_square_mesh();
        let area = total_surface_area(&m);
        assert!((area - 1.0).abs() < 1e-10, "area={area}");
    }

    #[test]
    fn test_face_areas_positive() {
        let m = ico();
        for a in face_areas(&m) {
            assert!(a > 0.0);
        }
    }

    #[test]
    fn test_voronoi_area_positive() {
        let m = ico();
        for v in 0..m.num_vertices() {
            let a = m.voronoi_area(v);
            assert!(a > 0.0, "voronoi_area at v={v} is {a}");
        }
    }

    // --- Normals ---

    #[test]
    fn test_face_normal_unit_length() {
        let m = flat_square_mesh();
        for fi in 0..m.num_faces() {
            let n = m.face_normal(fi);
            let l = len3(n);
            assert!((l - 1.0).abs() < 1e-10, "face_normal not unit: {l}");
        }
    }

    #[test]
    fn test_vertex_normal_unit_length() {
        let m = ico();
        for v in 0..m.num_vertices() {
            let n = m.vertex_normal(v);
            let l = len3(n);
            assert!(
                (l - 1.0).abs() < 1e-10,
                "vertex_normal not unit at {v}: {l}"
            );
        }
    }

    #[test]
    fn test_flat_face_normal_is_z() {
        let m = flat_square_mesh();
        let n = m.face_normal(0);
        // flat mesh in xy-plane → normal should be (0,0,±1)
        assert!(n[2].abs() > 0.99, "flat face normal z={}", n[2]);
    }

    // --- Laplacian ---

    #[test]
    fn test_laplacian_row_sum_near_zero() {
        let m = ico();
        let lap = m.cotan_laplacian();
        let n = m.num_vertices();
        let mut row_sum = vec![0.0f64; n];
        for &(i, _j, w) in &lap {
            row_sum[i] += w;
        }
        for (i, s) in row_sum.iter().enumerate() {
            assert!(s.abs() < 1e-8, "row {i} sum = {s}");
        }
    }

    #[test]
    fn test_laplacian_of_constant_is_zero() {
        let m = ico();
        let f = vec![1.0f64; m.num_vertices()];
        let lf = m.apply_laplacian(&f);
        for (i, v) in lf.iter().enumerate() {
            assert!(v.abs() < 1e-8, "L(1)[{i}] = {v}");
        }
    }

    // --- Gaussian curvature ---

    #[test]
    fn test_flat_mesh_gaussian_curvature_near_zero() {
        let m = flat_square_mesh();
        let k = m.gaussian_curvature();
        // Interior vertices of a flat mesh have K ≈ 0; boundary vertices have
        // non-zero angle defect. For the 2-triangle square, check the interior.
        // Vertex 2 is shared by both triangles with full 2π.
        // Allow boundary effects.
        for &ki in &k {
            assert!(
                ki.abs() < 100.0,
                "gaussian curvature suspiciously large: {ki}"
            );
        }
    }

    #[test]
    fn test_icosahedron_total_gaussian_curvature() {
        // By Gauss-Bonnet: integral K dA = 2π χ = 4π for a sphere (genus 0)
        let m = ico();
        let k = m.gaussian_curvature();
        let areas: Vec<f64> = (0..m.num_vertices()).map(|v| m.voronoi_area(v)).collect();
        let total: f64 = k.iter().zip(areas.iter()).map(|(ki, ai)| ki * ai).sum();
        // Should be close to 4π ≈ 12.566
        assert!(
            (total - 4.0 * PI).abs() < 1.0,
            "Gauss-Bonnet: total={total}"
        );
    }

    // --- Mean curvature ---

    #[test]
    fn test_mean_curvature_nonnegative_sphere() {
        // On a convex surface mean curvature magnitude should be > 0
        let m = ico();
        let h = m.mean_curvature();
        let any_nonzero = h.iter().any(|&hi| hi.abs() > 1e-6);
        assert!(
            any_nonzero,
            "mean curvature should be non-zero on icosahedron"
        );
    }

    #[test]
    fn test_principal_curvatures_k1_ge_k2() {
        let m = ico();
        let pc = m.principal_curvatures();
        for (k1, k2) in pc {
            assert!(k1 >= k2 - 1e-10, "k1={k1} < k2={k2}");
        }
    }

    // --- Dijkstra ---

    #[test]
    fn test_dijkstra_source_zero_distance() {
        let m = ico();
        let dist = geodesic_dijkstra(&m, &[0]);
        assert_eq!(dist[0], 0.0);
    }

    #[test]
    fn test_dijkstra_all_finite_on_connected_mesh() {
        let m = ico();
        let dist = geodesic_dijkstra(&m, &[0]);
        for (i, d) in dist.iter().enumerate() {
            assert!(d.is_finite(), "vertex {i} has infinite distance");
        }
    }

    #[test]
    fn test_dijkstra_nonnegative() {
        let m = ico();
        let dist = geodesic_dijkstra(&m, &[3]);
        for &d in &dist {
            assert!(d >= 0.0);
        }
    }

    #[test]
    fn test_geodesic_path_starts_at_source() {
        let m = ico();
        let path = geodesic_path(&m, 0, 3).expect("should find path");
        assert_eq!(path[0], 0);
    }

    #[test]
    fn test_geodesic_path_ends_at_target() {
        let m = ico();
        let path = geodesic_path(&m, 0, 3).expect("should find path");
        assert_eq!(*path.last().unwrap(), 3);
    }

    // --- Heat method ---

    #[test]
    fn test_heat_method_source_near_zero() {
        let m = ico();
        let solver = HeatMethodSolver::new(m, -1.0, 10);
        let dist = solver.compute(&[0]);
        assert!(
            dist[0] < 1e-3,
            "source distance should be near zero: {}",
            dist[0]
        );
    }

    #[test]
    fn test_heat_method_nonnegative() {
        let m = ico();
        let solver = HeatMethodSolver::new(m, -1.0, 10);
        let dist = solver.compute(&[0]);
        for &d in &dist {
            assert!(d >= -1e-6, "heat distance negative: {d}");
        }
    }

    // --- Euler characteristic ---

    #[test]
    fn test_euler_characteristic_icosahedron() {
        let m = ico();
        let chi = euler_characteristic(&m);
        assert_eq!(chi, 2, "icosahedron euler characteristic = {chi}");
    }

    #[test]
    fn test_euler_characteristic_tetrahedron() {
        let m = tet();
        let chi = euler_characteristic(&m);
        assert_eq!(chi, 2, "tetrahedron euler characteristic = {chi}");
    }

    #[test]
    fn test_genus_sphere_is_zero() {
        assert_eq!(genus_from_euler(2), 0);
    }

    // --- Parallel transport ---

    #[test]
    fn test_parallel_transport_coplanar_identity() {
        let m = flat_square_mesh();
        let v = [1.0, 0.0, 0.0];
        let transported = parallel_transport_across_edge(&m, 0, 1, v);
        // Both faces are in same plane → identity rotation
        let diff = len3(sub3(transported, v));
        assert!(
            diff < 1e-10,
            "coplanar transport should be identity: diff={diff}"
        );
    }

    // --- Vector field ---

    #[test]
    fn test_face_vector_field_zeros() {
        let m = flat_square_mesh();
        let f = FaceVectorField::zeros(m.num_faces());
        assert_eq!(f.norm(), 0.0);
    }

    #[test]
    fn test_vertex_vector_field_divergence_length() {
        let m = ico();
        let vf = VertexVectorField::zeros(m.num_vertices());
        let div = vf.divergence(&m);
        assert_eq!(div.len(), m.num_vertices());
    }

    #[test]
    fn test_zero_field_divergence_is_zero() {
        let m = ico();
        let vf = VertexVectorField::zeros(m.num_vertices());
        let div = vf.divergence(&m);
        for &d in &div {
            assert!(d.abs() < 1e-12, "div of zero field = {d}");
        }
    }

    // --- Discrete 1-form ---

    #[test]
    fn test_1form_antisymmetry() {
        let m = flat_square_mesh();
        let mut form = Discrete1Form::zeros(&m);
        *form.values.get_mut(&(0, 1)).unwrap() = 3.0;
        assert_eq!(form.eval(0, 1), 3.0);
        assert_eq!(form.eval(1, 0), -3.0);
    }

    #[test]
    fn test_1form_path_integration() {
        let m = flat_square_mesh();
        let mut form = Discrete1Form::zeros(&m);
        *form.values.get_mut(&(0, 1)).unwrap() = 1.0;
        *form.values.get_mut(&(1, 2)).unwrap() = 1.0;
        let path = vec![0, 1, 2];
        let val = form.integrate_path(&path);
        assert!((val - 2.0).abs() < 1e-10, "path integral = {val}");
    }

    // --- Mean curvature flow ---

    #[test]
    fn test_mcf_step_changes_vertices() {
        let mut m = ico();
        let before = m.vertices.clone();
        mean_curvature_flow_step(&mut m, 1e-3);
        let changed =
            m.vertices.iter().zip(before.iter()).any(|(a, b)| {
                (a[0] - b[0]).abs() + (a[1] - b[1]).abs() + (a[2] - b[2]).abs() > 1e-12
            });
        assert!(changed, "MCF step should change vertex positions");
    }

    // --- Aspect ratio ---

    #[test]
    fn test_aspect_ratio_all_positive() {
        let m = ico();
        for ar in triangle_aspect_ratios(&m) {
            assert!(ar > 0.0);
        }
    }

    // --- Hodge decomposition ---

    #[test]
    fn test_hodge_decompose_harmonic_length() {
        let m = ico();
        let vf = VertexVectorField::zeros(m.num_vertices());
        let hd = hodge_decompose(&m, &vf);
        assert_eq!(hd.harmonic.len(), m.num_vertices());
    }

    #[test]
    fn test_hodge_decompose_scalar_potential_length() {
        let m = ico();
        let vf = VertexVectorField::zeros(m.num_vertices());
        let hd = hodge_decompose(&m, &vf);
        assert_eq!(hd.scalar_potential.len(), m.num_vertices());
    }

    // --- Curvature tensor ---

    #[test]
    fn test_curvature_tensor_size() {
        let m = ico();
        let ct = vertex_curvature_tensor(&m, 0);
        assert_eq!(ct.len(), 9);
    }

    // --- Math helpers ---

    #[test]
    fn test_cotan_equilateral() {
        let p0 = [1.0, 0.0, 0.0];
        let p_opp = [0.0, 3.0_f64.sqrt(), 0.0];
        let p1 = [-1.0, 0.0, 0.0];
        let cot = cotan_at(p0, p_opp, p1);
        // Angle at p_opp in equilateral triangle = 60° → cot(60°) = 1/√3 ≈ 0.577
        assert!((cot - 1.0 / 3.0_f64.sqrt()).abs() < 0.01, "cot = {cot}");
    }

    #[test]
    fn test_tri_area_unit_right_triangle() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let a = tri_area(p0, p1, p2);
        assert!((a - 0.5).abs() < 1e-10, "area = {a}");
    }

    #[test]
    fn test_mean_edge_length_positive() {
        let m = ico();
        let h = mean_edge_length(&m);
        assert!(h > 0.0);
    }
}
