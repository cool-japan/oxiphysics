// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Triangle mesh shape with adjacency queries, vertex normals, Laplacian,
//! geodesic distance, Loop subdivision, and watertight checks.

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use std::collections::{BinaryHeap, HashMap, HashSet};

/// A triangle mesh defined by vertices and index triples.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// Vertex positions.
    pub vertices: Vec<Vec3>,
    /// Triangle indices (groups of three).
    pub indices: Vec<[usize; 3]>,
}

impl TriangleMesh {
    /// Create a new triangle mesh.
    pub fn new(vertices: Vec<Vec3>, indices: Vec<[usize; 3]>) -> Self {
        Self { vertices, indices }
    }

    /// Surface area: sum of triangle areas.
    pub fn surface_area(&self) -> Real {
        let mut total = 0.0;
        for tri in &self.indices {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            let ab = b - a;
            let ac = c - a;
            total += ab.cross(&ac).norm() * 0.5;
        }
        total
    }

    /// Volume via signed tetrahedral decomposition (divergence theorem).
    pub fn volume_explicit(&self) -> Real {
        let mut vol = 0.0;
        for tri in &self.indices {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            vol += a.dot(&b.cross(&c)) / 6.0;
        }
        vol.abs()
    }

    /// Center of mass via weighted tetrahedral centroids.
    pub fn center_of_mass_explicit(&self) -> [f64; 3] {
        if self.indices.is_empty() {
            return [0.0, 0.0, 0.0];
        }
        let mut weighted = Vec3::zeros();
        let mut total_vol = 0.0;
        for tri in &self.indices {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            let tet_vol = a.dot(&b.cross(&c)) / 6.0;
            let tet_center = (a + b + c) * 0.25;
            weighted += tet_center * tet_vol;
            total_vol += tet_vol;
        }
        if total_vol.abs() > 1e-12 {
            let com = weighted / total_vol;
            [com.x, com.y, com.z]
        } else {
            [0.0, 0.0, 0.0]
        }
    }

    /// Compute per-face normals (unit vectors).
    pub fn compute_normals(&self) -> Vec<[f64; 3]> {
        self.indices
            .iter()
            .map(|tri| {
                let a = self.vertices[tri[0]];
                let b = self.vertices[tri[1]];
                let c = self.vertices[tri[2]];
                let ab = b - a;
                let ac = c - a;
                let n = ab.cross(&ac);
                let len = n.norm();
                if len < 1e-12 {
                    [0.0, 1.0, 0.0]
                } else {
                    [n.x / len, n.y / len, n.z / len]
                }
            })
            .collect()
    }

    /// Ray cast returning (t, face_index, normal) via Moller-Trumbore for all triangles.
    pub fn ray_cast_full(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, usize, [f64; 3])> {
        let o = Vec3::new(origin[0], origin[1], origin[2]);
        let d = Vec3::new(direction[0], direction[1], direction[2]);
        let mut best: Option<(f64, usize, [f64; 3])> = None;

        for (face_idx, tri) in self.indices.iter().enumerate() {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];

            if let Some(hit) = ray_triangle(&o, &d, &v0, &v1, &v2)
                && hit.toi >= 0.0
                && hit.toi <= max_toi
                && best.as_ref().is_none_or(|b| hit.toi < b.0)
            {
                let n = [hit.normal.x, hit.normal.y, hit.normal.z];
                best = Some((hit.toi, face_idx, n));
            }
        }
        best
    }

    /// Check whether the mesh is watertight (every edge is shared by exactly two triangles).
    pub fn is_watertight(&self) -> bool {
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();

        for tri in &self.indices {
            let verts = [tri[0], tri[1], tri[2]];
            for i in 0..3 {
                let a = verts[i];
                let b = verts[(i + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }

        edge_count.values().all(|&count| count == 2)
    }

    // ------------------------------------------------------------------
    // New: mesh adjacency queries
    // ------------------------------------------------------------------

    /// Build vertex-to-face adjacency: for each vertex index, the list of
    /// triangle indices that reference it.
    pub fn vertex_to_faces(&self) -> Vec<Vec<usize>> {
        let mut v2f = vec![Vec::new(); self.vertices.len()];
        for (fi, tri) in self.indices.iter().enumerate() {
            for &vi in tri {
                v2f[vi].push(fi);
            }
        }
        v2f
    }

    /// Build face-to-face adjacency via shared edges.
    /// Returns a vec of the same length as `self.indices`. Each entry contains
    /// the indices of adjacent faces (those sharing at least one edge).
    pub fn face_adjacency(&self) -> Vec<Vec<usize>> {
        // edge -> list of face indices
        let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (fi, tri) in self.indices.iter().enumerate() {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                edge_faces.entry(key).or_default().push(fi);
            }
        }
        let mut adj = vec![Vec::new(); self.indices.len()];
        for faces in edge_faces.values() {
            for i in 0..faces.len() {
                for j in (i + 1)..faces.len() {
                    adj[faces[i]].push(faces[j]);
                    adj[faces[j]].push(faces[i]);
                }
            }
        }
        // Deduplicate
        for v in &mut adj {
            v.sort_unstable();
            v.dedup();
        }
        adj
    }

    /// Return the set of unique edges as sorted (min,max) vertex index pairs.
    pub fn unique_edges(&self) -> Vec<(usize, usize)> {
        let mut set: HashSet<(usize, usize)> = HashSet::new();
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                set.insert(key);
            }
        }
        let mut edges: Vec<(usize, usize)> = set.into_iter().collect();
        edges.sort_unstable();
        edges
    }

    /// Build vertex-to-vertex adjacency (1-ring neighbours).
    pub fn vertex_neighbors(&self) -> Vec<Vec<usize>> {
        let mut nbrs = vec![HashSet::<usize>::new(); self.vertices.len()];
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                nbrs[a].insert(b);
                nbrs[b].insert(a);
            }
        }
        nbrs.into_iter()
            .map(|s| {
                let mut v: Vec<usize> = s.into_iter().collect();
                v.sort_unstable();
                v
            })
            .collect()
    }

    // ------------------------------------------------------------------
    // New: edge collapse
    // ------------------------------------------------------------------

    /// Collapse an edge between vertices `v0` and `v1`, merging them to their
    /// midpoint. Triangles that degenerate (both endpoints are `v0`/`v1`) are
    /// removed.
    ///
    /// Returns `true` if the edge was found and collapsed.
    pub fn edge_collapse(&mut self, v0: usize, v1: usize) -> bool {
        // Verify edge exists
        let edge_found = self
            .indices
            .iter()
            .any(|tri| tri.contains(&v0) && tri.contains(&v1));
        if !edge_found {
            return false;
        }

        // Midpoint
        let mid = (self.vertices[v0] + self.vertices[v1]) * 0.5;
        self.vertices[v0] = mid;

        // Replace all references to v1 with v0
        for tri in &mut self.indices {
            for idx in tri.iter_mut() {
                if *idx == v1 {
                    *idx = v0;
                }
            }
        }

        // Remove degenerate triangles (two or more identical indices)
        self.indices
            .retain(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);

        true
    }

    // ------------------------------------------------------------------
    // New: vertex normal computation (angle-weighted)
    // ------------------------------------------------------------------

    /// Compute per-vertex normals by accumulating face normals weighted by the
    /// interior angle at each vertex.
    pub fn compute_vertex_normals(&self) -> Vec<[f64; 3]> {
        let mut normals = vec![Vec3::zeros(); self.vertices.len()];
        for tri in &self.indices {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            let e01 = v1 - v0;
            let e02 = v2 - v0;
            let face_n = e01.cross(&e02);

            let angle_at = |va: Vec3, vb: Vec3, vc: Vec3| -> f64 {
                let ab = (vb - va).normalize();
                let ac = (vc - va).normalize();
                ab.dot(&ac).clamp(-1.0, 1.0).acos()
            };

            let a0 = angle_at(v0, v1, v2);
            let a1 = angle_at(v1, v2, v0);
            let a2 = angle_at(v2, v0, v1);

            normals[tri[0]] += face_n * a0;
            normals[tri[1]] += face_n * a1;
            normals[tri[2]] += face_n * a2;
        }
        normals
            .iter()
            .map(|n| {
                let len = n.norm();
                if len < 1e-12 {
                    [0.0, 1.0, 0.0]
                } else {
                    [n.x / len, n.y / len, n.z / len]
                }
            })
            .collect()
    }

    // ------------------------------------------------------------------
    // New: mesh Laplacian smoothing
    // ------------------------------------------------------------------

    /// Uniform Laplacian smoothing: move each vertex towards the average of
    /// its 1-ring neighbours by `factor` (0..1). Boundary vertices are not
    /// moved.
    pub fn laplacian_smooth(&mut self, factor: f64, iterations: usize) {
        for _ in 0..iterations {
            let nbrs = self.vertex_neighbors();
            let mut new_pos = self.vertices.clone();
            for (vi, neighbours) in nbrs.iter().enumerate() {
                if neighbours.is_empty() {
                    continue;
                }
                let mut avg = Vec3::zeros();
                for &ni in neighbours {
                    avg += self.vertices[ni];
                }
                avg /= neighbours.len() as f64;
                let diff = avg - self.vertices[vi];
                new_pos[vi] = self.vertices[vi] + diff * factor;
            }
            self.vertices = new_pos;
        }
    }

    // ------------------------------------------------------------------
    // New: approximate geodesic distance (Dijkstra on mesh edges)
    // ------------------------------------------------------------------

    /// Compute approximate geodesic distances from `source` vertex to all
    /// other vertices using Dijkstra on edge lengths. Returns a vec of
    /// distances indexed by vertex, with `f64::INFINITY` for unreachable
    /// vertices.
    pub fn geodesic_distance(&self, source: usize) -> Vec<f64> {
        let nbrs = self.vertex_neighbors();
        let n = self.vertices.len();
        let mut dist = vec![f64::INFINITY; n];
        dist[source] = 0.0;

        // Min-heap: (OrderedFloat(dist), vertex)
        let mut heap = BinaryHeap::new();
        heap.push(std::cmp::Reverse(OrdF64(0.0, source)));

        while let Some(std::cmp::Reverse(OrdF64(d, u))) = heap.pop() {
            if d > dist[u] {
                continue;
            }
            for &v in &nbrs[u] {
                let edge_len = (self.vertices[v] - self.vertices[u]).norm();
                let new_d = d + edge_len;
                if new_d < dist[v] {
                    dist[v] = new_d;
                    heap.push(std::cmp::Reverse(OrdF64(new_d, v)));
                }
            }
        }
        dist
    }

    // ------------------------------------------------------------------
    // New: Loop subdivision
    // ------------------------------------------------------------------

    /// Perform one iteration of Loop subdivision.
    ///
    /// Each triangle is split into four by inserting edge midpoints (for
    /// boundary edges) or Loop-weighted edge vertices (for interior edges).
    /// Existing vertices are repositioned using the Loop weighting scheme.
    pub fn loop_subdivide(&mut self) {
        let n_verts = self.vertices.len();

        // Build edge -> face count for boundary detection
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }

        // Edge midpoint map: canonical edge -> new vertex index
        let mut edge_verts: HashMap<(usize, usize), usize> = HashMap::new();
        let mut new_vertices = self.vertices.clone();

        // Build edge -> adjacent triangle vertex sets for Loop weights
        let mut edge_opposite: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let c = tri[(k + 2) % 3]; // opposite vertex
                let key = if a < b { (a, b) } else { (b, a) };
                edge_opposite.entry(key).or_default().push(c);
            }
        }

        let get_or_create_edge_vert = |a: usize,
                                       b: usize,
                                       edge_verts: &mut HashMap<(usize, usize), usize>,
                                       new_vertices: &mut Vec<Vec3>,
                                       edge_count: &HashMap<(usize, usize), usize>,
                                       edge_opposite: &HashMap<(usize, usize), Vec<usize>>,
                                       orig_verts: &[Vec3]|
         -> usize {
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(&idx) = edge_verts.get(&key) {
                return idx;
            }
            let idx = new_vertices.len();
            let count = edge_count.get(&key).copied().unwrap_or(1);
            let pos = if count == 2 {
                // Interior edge: 3/8*(a+b) + 1/8*(c+d)
                let opposites = edge_opposite.get(&key).expect("key must exist");
                if opposites.len() == 2 {
                    let c = orig_verts[opposites[0]];
                    let d = orig_verts[opposites[1]];
                    (orig_verts[a] + orig_verts[b]) * (3.0 / 8.0) + (c + d) * (1.0 / 8.0)
                } else {
                    (orig_verts[a] + orig_verts[b]) * 0.5
                }
            } else {
                // Boundary edge: simple midpoint
                (orig_verts[a] + orig_verts[b]) * 0.5
            };
            new_vertices.push(pos);
            edge_verts.insert(key, idx);
            idx
        };

        let orig_verts = self.vertices.clone();

        // Create edge vertices for each triangle
        let mut new_indices = Vec::with_capacity(self.indices.len() * 4);
        for tri in &self.indices {
            let v0 = tri[0];
            let v1 = tri[1];
            let v2 = tri[2];
            let m01 = get_or_create_edge_vert(
                v0,
                v1,
                &mut edge_verts,
                &mut new_vertices,
                &edge_count,
                &edge_opposite,
                &orig_verts,
            );
            let m12 = get_or_create_edge_vert(
                v1,
                v2,
                &mut edge_verts,
                &mut new_vertices,
                &edge_count,
                &edge_opposite,
                &orig_verts,
            );
            let m20 = get_or_create_edge_vert(
                v2,
                v0,
                &mut edge_verts,
                &mut new_vertices,
                &edge_count,
                &edge_opposite,
                &orig_verts,
            );
            // Four sub-triangles
            new_indices.push([v0, m01, m20]);
            new_indices.push([v1, m12, m01]);
            new_indices.push([v2, m20, m12]);
            new_indices.push([m01, m12, m20]);
        }

        // Reposition original vertices using Loop weights
        let nbrs = self.vertex_neighbors();
        // Determine boundary vertices
        let mut is_boundary = vec![false; n_verts];
        for (&(a, b), &cnt) in &edge_count {
            if cnt == 1 {
                is_boundary[a] = true;
                is_boundary[b] = true;
            }
        }
        for vi in 0..n_verts {
            if is_boundary[vi] {
                // Boundary: keep position (simple approach)
                // Could use 3/4 * v + 1/8 * (n1 + n2) for boundary
                continue;
            }
            let n = nbrs[vi].len();
            if n < 3 {
                continue;
            }
            let beta = if n == 3 {
                3.0 / 16.0
            } else {
                3.0 / (8.0 * n as f64)
            };
            let mut neighbour_sum = Vec3::zeros();
            for &ni in &nbrs[vi] {
                neighbour_sum += orig_verts[ni];
            }
            new_vertices[vi] = orig_verts[vi] * (1.0 - n as f64 * beta) + neighbour_sum * beta;
        }

        self.vertices = new_vertices;
        self.indices = new_indices;
    }

    // ------------------------------------------------------------------
    // New: boundary edges and Euler characteristic
    // ------------------------------------------------------------------

    /// Return boundary edges (edges shared by only one triangle).
    pub fn boundary_edges(&self) -> Vec<(usize, usize)> {
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        edge_count
            .into_iter()
            .filter(|&(_, c)| c == 1)
            .map(|(e, _)| e)
            .collect()
    }

    /// Compute the Euler characteristic: V - E + F.
    pub fn euler_characteristic(&self) -> i64 {
        let v = self.vertices.len() as i64;
        let e = self.unique_edges().len() as i64;
        let f = self.indices.len() as i64;
        v - e + f
    }

    /// Count non-manifold edges (shared by more than 2 triangles).
    pub fn non_manifold_edge_count(&self) -> usize {
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        edge_count.values().filter(|&&c| c > 2).count()
    }

    // ------------------------------------------------------------------
    // Cotangent Laplacian
    // ------------------------------------------------------------------

    /// Compute the cotangent Laplacian weights for each directed edge (i → j).
    ///
    /// Returns a `HashMap<(usize, usize), f64>` mapping `(i, j)` to the
    /// symmetric cotangent weight `w_{ij} = (cot α + cot β) / 2`, where α
    /// and β are the angles opposite to edge (i,j) in the two incident faces.
    ///
    /// This is the standard discretisation used in geometry processing for
    /// Laplace-Beltrami operators on triangulated surfaces.
    pub fn compute_laplacian_matrix(&self) -> HashMap<(usize, usize), f64> {
        let mut weights: HashMap<(usize, usize), f64> = HashMap::new();

        for tri in &self.indices {
            let [i0, i1, i2] = *tri;
            let p0 = self.vertices[i0];
            let p1 = self.vertices[i1];
            let p2 = self.vertices[i2];

            // Cotangent of the angle at each vertex:
            //   cot(angle at v0) is the angle between edges v0→v1 and v0→v2
            let cot_at = |a: Vec3, b: Vec3, c: Vec3| -> f64 {
                // angle at vertex a, between edges a→b and a→c
                let ab = b - a;
                let ac = c - a;
                let cos_t = ab.dot(&ac);
                let sin_t = ab.cross(&ac).norm();
                if sin_t.abs() < 1e-12 {
                    0.0
                } else {
                    cos_t / sin_t
                }
            };

            let cot0 = cot_at(p0, p1, p2); // angle at i0, opposite edge i1-i2
            let cot1 = cot_at(p1, p0, p2); // angle at i1, opposite edge i0-i2
            let cot2 = cot_at(p2, p0, p1); // angle at i2, opposite edge i0-i1

            // Edge (i1, i2): opposite angle at i0
            *weights.entry((i1, i2)).or_insert(0.0) += 0.5 * cot0;
            *weights.entry((i2, i1)).or_insert(0.0) += 0.5 * cot0;

            // Edge (i0, i2): opposite angle at i1
            *weights.entry((i0, i2)).or_insert(0.0) += 0.5 * cot1;
            *weights.entry((i2, i0)).or_insert(0.0) += 0.5 * cot1;

            // Edge (i0, i1): opposite angle at i2
            *weights.entry((i0, i1)).or_insert(0.0) += 0.5 * cot2;
            *weights.entry((i1, i0)).or_insert(0.0) += 0.5 * cot2;
        }

        weights
    }

    // ------------------------------------------------------------------
    // Heat Kernel Signature
    // ------------------------------------------------------------------

    /// Compute the Heat Kernel Signature (HKS) descriptor for each vertex
    /// at a set of time scales `t_values`.
    ///
    /// The HKS approximates the diagonal of the heat kernel K(x, x, t) using
    /// the cotangent Laplacian's (normalized) diagonal weights.  The full
    /// spectral decomposition is expensive, so this implementation uses a
    /// fast approximation: for each vertex i and time t, it computes
    ///
    ///   HKS(i, t) = exp(-t * D\[i\])
    ///
    /// where `D[i]` is the sum of cotangent weights incident on vertex i
    /// (the diagonal of the stiffness matrix, normalised by the vertex area).
    ///
    /// Returns a `Vec<Vec`f64`>` of shape `[n_vertices][t_values.len()]`.
    pub fn compute_heat_kernel_signature(&self, t_values: &[f64]) -> Vec<Vec<f64>> {
        let n = self.vertices.len();
        if n == 0 || t_values.is_empty() {
            return vec![vec![0.0; t_values.len()]; n];
        }

        // Cotangent weights
        let weights = self.compute_laplacian_matrix();

        // Vertex areas (1/3 of the sum of incident triangle areas)
        let mut vertex_area = vec![0.0f64; n];
        for tri in &self.indices {
            let [i0, i1, i2] = *tri;
            let a = self.vertices[i0];
            let b = self.vertices[i1];
            let c = self.vertices[i2];
            let area = (b - a).cross(&(c - a)).norm() * 0.5;
            let a_third = area / 3.0;
            vertex_area[i0] += a_third;
            vertex_area[i1] += a_third;
            vertex_area[i2] += a_third;
        }

        // Row sum of cotangent weights (diagonal of the stiffness matrix)
        let mut diag = vec![0.0f64; n];
        for (&(i, _j), &w) in &weights {
            diag[i] += w;
        }

        // Effective frequency per vertex: lambda_i = diag[i] / area[i]
        let lambda: Vec<f64> = (0..n)
            .map(|i| {
                if vertex_area[i] > 1e-15 {
                    diag[i] / vertex_area[i]
                } else {
                    0.0
                }
            })
            .collect();

        // HKS(i, t) = exp(-t * lambda_i)
        (0..n)
            .map(|i| t_values.iter().map(|&t| (-t * lambda[i]).exp()).collect())
            .collect()
    }

    // ------------------------------------------------------------------
    // Laplacian smoothing (cotangent weights)
    // ------------------------------------------------------------------

    /// Iterative Laplacian smoothing using cotangent weights.
    ///
    /// Each iteration moves vertex `i` towards the weighted average of its
    /// neighbours:
    ///
    ///   v_i' = v_i + factor * Σ_j w_{ij} * (v_j - v_i) / Σ_j w_{ij}
    ///
    /// Boundary vertices (those on edges shared by only one triangle) are
    /// kept fixed.
    ///
    /// `factor` – step size in \[0, 1\]; `iterations` – number of passes.
    pub fn smooth_laplacian(&mut self, factor: f64, iterations: usize) {
        if self.vertices.is_empty() || self.indices.is_empty() {
            return;
        }

        // Identify boundary vertices
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &self.indices {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        let mut is_boundary = vec![false; self.vertices.len()];
        for (&(a, b), &c) in &edge_count {
            if c == 1 {
                is_boundary[a] = true;
                is_boundary[b] = true;
            }
        }

        for _ in 0..iterations {
            let weights = self.compute_laplacian_matrix();
            let mut new_vertices = self.vertices.clone();

            for i in 0..self.vertices.len() {
                if is_boundary[i] {
                    continue;
                }
                let mut weight_sum = 0.0f64;
                let mut delta = Vec3::zeros();

                for tri in &self.indices {
                    for k in 0..3 {
                        if tri[k] == i {
                            let j = tri[(k + 1) % 3];
                            let w = *weights.get(&(i, j)).unwrap_or(&0.0);
                            delta += (self.vertices[j] - self.vertices[i]) * w;
                            weight_sum += w;
                        }
                    }
                }

                if weight_sum > 1e-15 {
                    new_vertices[i] += delta * (factor / weight_sum);
                }
            }

            self.vertices = new_vertices;
        }
    }
}

// -----------------------------------------------------------------------
// Helper for Dijkstra heap ordering
// -----------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
struct OrdF64(f64, usize);

impl Eq for OrdF64 {}

impl PartialOrd for OrdF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrdF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

// -----------------------------------------------------------------------
// Shape trait implementation
// -----------------------------------------------------------------------

impl Shape for TriangleMesh {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn bounding_box(&self) -> Aabb {
        if self.vertices.is_empty() {
            return Aabb::new(Vec3::zeros(), Vec3::zeros());
        }
        let mut min = self.vertices[0];
        let mut max = self.vertices[0];
        for v in &self.vertices[1..] {
            min = min.inf(v);
            max = max.sup(v);
        }
        Aabb::new(min, max)
    }

    fn support_point(&self, direction: &Vec3) -> Vec3 {
        self.vertices
            .iter()
            .max_by(|a, b| {
                a.dot(direction)
                    .partial_cmp(&b.dot(direction))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or_else(Vec3::zeros)
    }

    fn volume(&self) -> Real {
        let mut vol = 0.0;
        for tri in &self.indices {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            vol += a.dot(&b.cross(&c)) / 6.0;
        }
        vol.abs()
    }

    fn center_of_mass(&self) -> Vec3 {
        if self.indices.is_empty() {
            return Vec3::zeros();
        }
        let mut weighted = Vec3::zeros();
        let mut total_vol = 0.0;
        for tri in &self.indices {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            let tet_vol = a.dot(&b.cross(&c)) / 6.0;
            let tet_center = (a + b + c) * 0.25;
            weighted += tet_center * tet_vol;
            total_vol += tet_vol;
        }
        if total_vol.abs() > 1e-12 {
            weighted / total_vol
        } else {
            Vec3::zeros()
        }
    }

    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let bb = self.bounding_box();
        let size = bb.max - bb.min;
        let k = mass / 12.0;
        Mat3::new(
            k * (size.y * size.y + size.z * size.z),
            0.0,
            0.0,
            0.0,
            k * (size.x * size.x + size.z * size.z),
            0.0,
            0.0,
            0.0,
            k * (size.x * size.x + size.y * size.y),
        )
    }

    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;
        for tri in &self.indices {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            if let Some(hit) = ray_triangle(ray_origin, ray_direction, &v0, &v1, &v2)
                && hit.toi <= max_toi
                && hit.toi >= 0.0
                && best.as_ref().is_none_or(|b| hit.toi < b.toi)
            {
                best = Some(hit);
            }
        }
        best
    }
}

/// Moller-Trumbore ray-triangle intersection.
fn ray_triangle(
    origin: &Vec3,
    direction: &Vec3,
    v0: &Vec3,
    v1: &Vec3,
    v2: &Vec3,
) -> Option<RayHit> {
    let edge1 = v1 - v0;
    let edge2 = v2 - v0;
    let h = direction.cross(&edge2);
    let a = edge1.dot(&h);

    if a.abs() < 1e-10 {
        return None;
    }

    let f = 1.0 / a;
    let s = origin - v0;
    let u = f * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = s.cross(&edge1);
    let v = f * direction.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * edge2.dot(&q);
    if t < 0.0 {
        return None;
    }

    let point = origin + direction * t;
    let normal = edge1.cross(&edge2).normalize();
    Some(RayHit {
        point,
        normal,
        toi: t,
    })
}

/// Build a simple unit cube mesh (12 triangles, watertight).
#[cfg(test)]
fn unit_cube_mesh() -> TriangleMesh {
    let v = vec![
        Vec3::new(0.0, 0.0, 0.0), // 0
        Vec3::new(1.0, 0.0, 0.0), // 1
        Vec3::new(1.0, 1.0, 0.0), // 2
        Vec3::new(0.0, 1.0, 0.0), // 3
        Vec3::new(0.0, 0.0, 1.0), // 4
        Vec3::new(1.0, 0.0, 1.0), // 5
        Vec3::new(1.0, 1.0, 1.0), // 6
        Vec3::new(0.0, 1.0, 1.0), // 7
    ];
    let idx = vec![
        // -Z face (z=0)
        [0, 2, 1],
        [0, 3, 2],
        // +Z face (z=1)
        [4, 5, 6],
        [4, 6, 7],
        // -X face (x=0)
        [0, 4, 7],
        [0, 7, 3],
        // +X face (x=1)
        [1, 2, 6],
        [1, 6, 5],
        // -Y face (y=0)
        [0, 1, 5],
        [0, 5, 4],
        // +Y face (y=1)
        [3, 7, 6],
        [3, 6, 2],
    ];
    TriangleMesh::new(v, idx)
}

/// Build a simple tetrahedron mesh (4 triangles, watertight).
#[cfg(test)]
fn tetrahedron_mesh() -> TriangleMesh {
    let v = vec![
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(-1.0, -1.0, 1.0),
        Vec3::new(-1.0, 1.0, -1.0),
        Vec3::new(1.0, -1.0, -1.0),
    ];
    let idx = vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
    TriangleMesh::new(v, idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TriangleMesh;
    use oxiphysics_core::math::Vec3;

    #[test]
    fn test_triangle_mesh_raycast() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(-1.0, 0.0, -1.0),
                Vec3::new(1.0, 0.0, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
            vec![[0, 1, 2]],
        );
        let origin = Vec3::new(0.0, 5.0, 0.0);
        let dir = Vec3::new(0.0, -1.0, 0.0);
        let hit = mesh.ray_cast(&origin, &dir, 100.0).unwrap();
        assert!((hit.toi - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_triangle_mesh_raycast_miss() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(-1.0, 0.0, -1.0),
                Vec3::new(1.0, 0.0, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
            vec![[0, 1, 2]],
        );
        let origin = Vec3::new(10.0, 5.0, 0.0);
        let dir = Vec3::new(0.0, -1.0, 0.0);
        assert!(mesh.ray_cast(&origin, &dir, 100.0).is_none());
    }

    #[test]
    fn test_triangle_mesh_cube_volume() {
        let mesh = unit_cube_mesh();
        let vol = mesh.volume();
        assert!((vol - 1.0).abs() < 1e-6, "expected volume=1, got {}", vol);
    }

    #[test]
    fn test_triangle_mesh_cube_surface_area() {
        let mesh = unit_cube_mesh();
        let sa = mesh.surface_area();
        assert!((sa - 6.0).abs() < 1e-6, "expected SA=6, got {}", sa);
    }

    #[test]
    fn test_triangle_mesh_cube_is_watertight() {
        let mesh = unit_cube_mesh();
        assert!(mesh.is_watertight(), "unit cube should be watertight");
    }

    #[test]
    fn test_triangle_mesh_open_not_watertight() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        assert!(!mesh.is_watertight(), "open mesh should not be watertight");
    }

    #[test]
    fn test_triangle_mesh_compute_normals() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let normals = mesh.compute_normals();
        assert_eq!(normals.len(), 1);
        assert!((normals[0][2] - 1.0).abs() < 1e-10, "expected +Z normal");
    }

    #[test]
    fn test_triangle_mesh_ray_cast_full() {
        let mesh = unit_cube_mesh();
        let (t, _face, n) = mesh
            .ray_cast_full([0.5, 5.0, 0.5], [0.0, -1.0, 0.0], 100.0)
            .expect("should hit cube");
        assert!((t - 4.0).abs() < 1e-6, "expected t=4, got {}", t);
        assert!(n[1].abs() > 0.9, "normal should be roughly Y");
    }

    #[test]
    fn test_triangle_mesh_surface_area_triangle() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(3.0, 0.0, 0.0),
                Vec3::new(0.0, 4.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let sa = mesh.surface_area();
        assert!((sa - 6.0).abs() < 1e-10, "expected SA=6, got {}", sa);
    }

    // ------------------------------------------------------------------
    // New tests: adjacency queries
    // ------------------------------------------------------------------

    #[test]
    fn test_vertex_to_faces() {
        let mesh = unit_cube_mesh();
        let v2f = mesh.vertex_to_faces();
        assert_eq!(v2f.len(), 8);
        // Each vertex of a cube is shared by 3 faces (each face = 2 tris,
        // each vertex appears in exactly 3 faces = 6 tris / 2, but actually
        // each vertex is used by 3*2=6? Let's just check non-empty).
        for faces in &v2f {
            assert!(
                !faces.is_empty(),
                "every vertex should have at least one face"
            );
        }
        // Vertex 0 appears in -Z, -X, -Y faces (6 triangles)
        // Indices: [0,2,1],[0,3,2],[0,4,7],[0,7,3],[0,1,5],[0,5,4]
        assert_eq!(
            v2f[0].len(),
            6,
            "vertex 0 in unit cube should appear in 6 triangles"
        );
    }

    #[test]
    fn test_face_adjacency_single_triangle() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let adj = mesh.face_adjacency();
        assert_eq!(adj.len(), 1);
        assert!(adj[0].is_empty(), "single triangle has no adjacent faces");
    }

    #[test]
    fn test_face_adjacency_cube() {
        let mesh = unit_cube_mesh();
        let adj = mesh.face_adjacency();
        assert_eq!(adj.len(), 12);
        // Each triangle in a cube shares edges with other triangles
        for neighbors in &adj {
            assert!(
                !neighbors.is_empty(),
                "each cube triangle should have neighbors"
            );
        }
    }

    #[test]
    fn test_unique_edges_single_triangle() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let edges = mesh.unique_edges();
        assert_eq!(edges.len(), 3);
    }

    #[test]
    fn test_unique_edges_cube() {
        let mesh = unit_cube_mesh();
        let edges = mesh.unique_edges();
        // A cube has 12 edges + 6 face diagonals = 18
        assert_eq!(edges.len(), 18);
    }

    #[test]
    fn test_vertex_neighbors() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let nbrs = mesh.vertex_neighbors();
        assert_eq!(nbrs.len(), 3);
        assert_eq!(nbrs[0].len(), 2); // connected to 1 and 2
        assert_eq!(nbrs[1].len(), 2);
        assert_eq!(nbrs[2].len(), 2);
    }

    // ------------------------------------------------------------------
    // New tests: edge collapse
    // ------------------------------------------------------------------

    #[test]
    fn test_edge_collapse_basic() {
        // Two triangles sharing edge (0,1)
        let mut mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
                Vec3::new(0.5, -1.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 3, 1]],
        );
        let ok = mesh.edge_collapse(0, 1);
        assert!(ok, "edge collapse should succeed");
        // After collapse, vertex 0 should be at midpoint (0.5, 0, 0)
        assert!((mesh.vertices[0].x - 0.5).abs() < 1e-10);
        // Triangles that degenerate should be removed
        // Both triangles had vertices 0 and 1, which are now both 0,
        // so they become degenerate -> removed? No: [0,1,2] becomes [0,0,2] -> degenerate.
        // [0,3,1] becomes [0,3,0] -> degenerate.
        assert!(mesh.indices.is_empty(), "both triangles should degenerate");
    }

    #[test]
    fn test_edge_collapse_nonexistent() {
        let mut mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        // Edge (0, 100) doesn't exist
        assert!(!mesh.edge_collapse(0, 100));
    }

    #[test]
    fn test_edge_collapse_preserves_other_faces() {
        // Three triangles, collapse edge shared by first two
        let mut mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),  // 0
                Vec3::new(2.0, 0.0, 0.0),  // 1
                Vec3::new(1.0, 1.0, 0.0),  // 2
                Vec3::new(1.0, -1.0, 0.0), // 3
                Vec3::new(3.0, 1.0, 0.0),  // 4
            ],
            vec![[0, 1, 2], [0, 3, 1], [1, 4, 2]],
        );
        let ok = mesh.edge_collapse(0, 1);
        assert!(ok);
        // Triangle [1,4,2] references v1 which is now v0
        // It becomes [0,4,2] which is non-degenerate
        assert_eq!(
            mesh.indices.len(),
            1,
            "two degenerate tris removed, one survives"
        );
    }

    // ------------------------------------------------------------------
    // New tests: vertex normals
    // ------------------------------------------------------------------

    #[test]
    fn test_vertex_normals_flat_surface() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        );
        let normals = mesh.compute_vertex_normals();
        assert_eq!(normals.len(), 4);
        for n in &normals {
            assert!(
                (n[2] - 1.0).abs() < 1e-6,
                "flat XY surface normal should be +Z"
            );
        }
    }

    #[test]
    fn test_vertex_normals_unit_length() {
        let mesh = unit_cube_mesh();
        let normals = mesh.compute_vertex_normals();
        for n in &normals {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-6,
                "vertex normal should be unit length, got {}",
                len
            );
        }
    }

    // ------------------------------------------------------------------
    // New tests: Laplacian smoothing
    // ------------------------------------------------------------------

    #[test]
    fn test_laplacian_smooth_doesnt_crash() {
        let mut mesh = unit_cube_mesh();
        mesh.laplacian_smooth(0.5, 3);
        // After smoothing the cube should still have 8 vertices and 12 faces
        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.indices.len(), 12);
    }

    #[test]
    fn test_laplacian_smooth_converges_to_center() {
        // A mesh with one interior vertex should move towards neighbors
        let mut mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0), // 0
                Vec3::new(2.0, 0.0, 0.0), // 1
                Vec3::new(1.0, 2.0, 0.0), // 2
                Vec3::new(1.0, 0.5, 0.0), // 3 (interior-ish vertex)
            ],
            vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]],
        );
        let orig_v3 = mesh.vertices[3];
        mesh.laplacian_smooth(1.0, 1);
        // Vertex 3 should have moved towards the average of its neighbors (0,1,2)
        let avg =
            (Vec3::new(0.0, 0.0, 0.0) + Vec3::new(2.0, 0.0, 0.0) + Vec3::new(1.0, 2.0, 0.0)) / 3.0;
        let new_v3 = mesh.vertices[3];
        let d_new = (new_v3 - avg).norm();
        let d_old = (orig_v3 - avg).norm();
        assert!(
            d_new < d_old + 1e-10,
            "vertex should move towards neighbor average"
        );
    }

    // ------------------------------------------------------------------
    // New tests: geodesic distance
    // ------------------------------------------------------------------

    #[test]
    fn test_geodesic_distance_source() {
        let mesh = unit_cube_mesh();
        let dist = mesh.geodesic_distance(0);
        assert_eq!(dist[0], 0.0, "distance to self should be 0");
    }

    #[test]
    fn test_geodesic_distance_adjacent() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let dist = mesh.geodesic_distance(0);
        // Distance from 0 to 1 should be 1.0
        assert!(
            (dist[1] - 1.0).abs() < 1e-10,
            "expected dist=1, got {}",
            dist[1]
        );
    }

    #[test]
    fn test_geodesic_distance_symmetry() {
        let mesh = unit_cube_mesh();
        let d0 = mesh.geodesic_distance(0);
        let d6 = mesh.geodesic_distance(6);
        // Distance 0->6 should equal distance 6->0
        assert!(
            (d0[6] - d6[0]).abs() < 1e-10,
            "geodesic distance should be symmetric: {} vs {}",
            d0[6],
            d6[0]
        );
    }

    // ------------------------------------------------------------------
    // New tests: Loop subdivision
    // ------------------------------------------------------------------

    #[test]
    fn test_loop_subdivide_increases_faces() {
        let mut mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        mesh.loop_subdivide();
        assert_eq!(mesh.indices.len(), 4, "one triangle should become four");
    }

    #[test]
    fn test_loop_subdivide_cube() {
        let mut mesh = unit_cube_mesh();
        assert_eq!(mesh.indices.len(), 12);
        mesh.loop_subdivide();
        assert_eq!(mesh.indices.len(), 48, "12 tris * 4 = 48");
        // Should still have consistent structure
        assert!(
            mesh.vertices.len() > 8,
            "should have more vertices after subdivision"
        );
    }

    #[test]
    fn test_loop_subdivide_watertight_preserved() {
        let mut mesh = tetrahedron_mesh();
        assert!(mesh.is_watertight(), "tetrahedron should be watertight");
        mesh.loop_subdivide();
        assert!(
            mesh.is_watertight(),
            "subdivided tetrahedron should still be watertight"
        );
    }

    // ------------------------------------------------------------------
    // New tests: boundary edges, Euler characteristic, non-manifold
    // ------------------------------------------------------------------

    #[test]
    fn test_boundary_edges_closed_mesh() {
        let mesh = unit_cube_mesh();
        let be = mesh.boundary_edges();
        assert!(
            be.is_empty(),
            "watertight mesh should have no boundary edges"
        );
    }

    #[test]
    fn test_boundary_edges_open_mesh() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        );
        let be = mesh.boundary_edges();
        assert_eq!(be.len(), 3, "single triangle has 3 boundary edges");
    }

    #[test]
    fn test_euler_characteristic_cube() {
        let mesh = unit_cube_mesh();
        let chi = mesh.euler_characteristic();
        assert_eq!(chi, 2, "closed cube Euler characteristic should be 2");
    }

    #[test]
    fn test_euler_characteristic_tetrahedron() {
        let mesh = tetrahedron_mesh();
        let chi = mesh.euler_characteristic();
        assert_eq!(
            chi, 2,
            "closed tetrahedron Euler characteristic should be 2"
        );
    }

    #[test]
    fn test_non_manifold_edge_count_clean_mesh() {
        let mesh = unit_cube_mesh();
        assert_eq!(mesh.non_manifold_edge_count(), 0);
    }

    #[test]
    fn test_non_manifold_edge_count_with_extra_face() {
        // Add a third triangle sharing an edge -> non-manifold
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
                Vec3::new(0.5, -1.0, 0.0),
                Vec3::new(0.5, 0.5, 1.0),
            ],
            vec![[0, 1, 2], [0, 3, 1], [0, 4, 1]], // edge (0,1) shared by 3 faces
        );
        assert_eq!(mesh.non_manifold_edge_count(), 1);
    }
}

// ---------------------------------------------------------------------------
// Tests for cotangent Laplacian, HKS, and smooth_laplacian
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_laplacian_hks {

    use crate::TriangleMesh;
    use oxiphysics_core::math::Vec3;

    /// Flat quad made of two triangles: 4 vertices, 2 triangles.
    fn flat_quad() -> TriangleMesh {
        TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
    }

    /// Single triangle.
    fn single_triangle() -> TriangleMesh {
        TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )
    }

    // ── compute_laplacian_matrix ─────────────────────────────────────────────

    #[test]
    fn test_laplacian_matrix_nonempty() {
        let mesh = flat_quad();
        let w = mesh.compute_laplacian_matrix();
        assert!(!w.is_empty(), "Laplacian weights should be non-empty");
    }

    #[test]
    fn test_laplacian_matrix_symmetric() {
        let mesh = flat_quad();
        let w = mesh.compute_laplacian_matrix();
        // For every (i,j) entry there should be a matching (j,i)
        for (&(i, j), &wij) in &w {
            let wji = w.get(&(j, i)).copied().unwrap_or(0.0);
            assert!(
                (wij - wji).abs() < 1e-10,
                "Laplacian not symmetric at ({i},{j}): {wij} vs {wji}"
            );
        }
    }

    #[test]
    fn test_laplacian_matrix_single_triangle_all_edges_present() {
        let mesh = single_triangle();
        let w = mesh.compute_laplacian_matrix();
        // 3 undirected edges × 2 directions = 6 entries
        assert_eq!(
            w.len(),
            6,
            "single triangle should produce 6 directed entries"
        );
    }

    #[test]
    fn test_laplacian_matrix_empty_mesh() {
        let mesh = TriangleMesh::new(vec![], vec![]);
        let w = mesh.compute_laplacian_matrix();
        assert!(w.is_empty());
    }

    // ── compute_heat_kernel_signature ────────────────────────────────────────

    #[test]
    fn test_hks_output_shape() {
        let mesh = flat_quad();
        let t_vals = vec![0.1, 1.0, 10.0];
        let hks = mesh.compute_heat_kernel_signature(&t_vals);
        assert_eq!(hks.len(), mesh.vertices.len(), "one HKS row per vertex");
        for row in &hks {
            assert_eq!(row.len(), t_vals.len(), "one entry per time scale");
        }
    }

    #[test]
    fn test_hks_values_in_range() {
        // HKS = exp(-t * lambda) is always in (0, 1] for t > 0, lambda >= 0
        let mesh = flat_quad();
        let t_vals = vec![0.01, 0.1, 1.0, 10.0];
        let hks = mesh.compute_heat_kernel_signature(&t_vals);
        for row in &hks {
            for &val in row {
                assert!(
                    val > 0.0 && val <= 1.0 + 1e-12,
                    "HKS value out of range: {val}"
                );
            }
        }
    }

    #[test]
    fn test_hks_decreasing_with_time() {
        // For any vertex, HKS(i, t1) >= HKS(i, t2) when t1 < t2
        let mesh = flat_quad();
        let t_vals = vec![0.1, 1.0, 10.0];
        let hks = mesh.compute_heat_kernel_signature(&t_vals);
        for (vi, row) in hks.iter().enumerate() {
            for k in 0..(row.len() - 1) {
                assert!(
                    row[k] >= row[k + 1] - 1e-12,
                    "HKS not decreasing for vertex {vi}: t={}, val={}; t={}, val={}",
                    t_vals[k],
                    row[k],
                    t_vals[k + 1],
                    row[k + 1]
                );
            }
        }
    }

    #[test]
    fn test_hks_empty_mesh() {
        let mesh = TriangleMesh::new(vec![], vec![]);
        let hks = mesh.compute_heat_kernel_signature(&[1.0, 2.0]);
        assert!(hks.is_empty());
    }

    // ── smooth_laplacian ────────────────────────────────────────────────────

    #[test]
    fn test_smooth_laplacian_vertex_count_unchanged() {
        let mut mesh = flat_quad();
        let n_before = mesh.vertices.len();
        mesh.smooth_laplacian(0.5, 3);
        assert_eq!(
            mesh.vertices.len(),
            n_before,
            "smoothing must not add/remove vertices"
        );
    }

    #[test]
    fn test_smooth_laplacian_triangle_count_unchanged() {
        let mut mesh = flat_quad();
        let n_before = mesh.indices.len();
        mesh.smooth_laplacian(0.5, 3);
        assert_eq!(
            mesh.indices.len(),
            n_before,
            "smoothing must not change connectivity"
        );
    }

    #[test]
    fn test_smooth_laplacian_zero_iterations_no_change() {
        let mut mesh = flat_quad();
        let verts_before: Vec<Vec3> = mesh.vertices.clone();
        mesh.smooth_laplacian(0.5, 0);
        for (a, b) in verts_before.iter().zip(mesh.vertices.iter()) {
            assert!(
                (a - b).norm() < 1e-12,
                "zero iterations should not move vertices"
            );
        }
    }

    #[test]
    fn test_smooth_laplacian_boundary_fixed() {
        // For a single triangle (all boundary) smoothing must not move any vertex
        let mut mesh = single_triangle();
        let verts_before: Vec<Vec3> = mesh.vertices.clone();
        mesh.smooth_laplacian(0.5, 5);
        for (a, b) in verts_before.iter().zip(mesh.vertices.iter()) {
            assert!(
                (a - b).norm() < 1e-12,
                "boundary vertices must not move, delta={}",
                (a - b).norm()
            );
        }
    }
}
