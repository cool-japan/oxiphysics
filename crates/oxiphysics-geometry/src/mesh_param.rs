// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Triangle mesh parameterization (UV unwrapping) and mesh processing utilities.
//!
//! Implements LSCM and Tutte parameterization, UV quality metrics,
//! mesh subdivision, and mesh repair operations.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn length3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = length3(a);
    if l < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

// ---------------------------------------------------------------------------
// TriMesh (parameterization mesh)
// ---------------------------------------------------------------------------

/// A triangle mesh with UV coordinates and normals, used for parameterization.
///
/// Note: This is distinct from `mesh_ops::TriMesh` which uses nalgebra Vec3.
#[derive(Debug, Clone)]
pub struct ParamTriMesh {
    /// 3D positions of vertices.
    pub positions: Vec<[f64; 3]>,
    /// UV texture coordinates per vertex.
    pub uvs: Vec<[f64; 2]>,
    /// Per-vertex normals.
    pub normals: Vec<[f64; 3]>,
    /// Triangle indices.
    pub triangles: Vec<[usize; 3]>,
}

impl ParamTriMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            uvs: Vec::new(),
            normals: Vec::new(),
            triangles: Vec::new(),
        }
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// Add a vertex with position and UV; returns the vertex index.
    pub fn add_vertex(&mut self, pos: [f64; 3], uv: [f64; 2]) -> usize {
        let idx = self.positions.len();
        self.positions.push(pos);
        self.uvs.push(uv);
        self.normals.push([0.0, 0.0, 0.0]);
        idx
    }

    /// Add a triangle by vertex indices.
    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize) {
        self.triangles.push([a, b, c]);
    }

    /// Compute per-vertex area-weighted average normals.
    pub fn compute_normals(&mut self) {
        let n = self.positions.len();
        let mut acc = vec![[0.0f64; 3]; n];

        for tri in &self.triangles {
            let p0 = self.positions[tri[0]];
            let p1 = self.positions[tri[1]];
            let p2 = self.positions[tri[2]];
            let e1 = sub3(p1, p0);
            let e2 = sub3(p2, p0);
            let face_n = cross3(e1, e2); // magnitude = 2 * area
            for &v in tri {
                acc[v][0] += face_n[0];
                acc[v][1] += face_n[1];
                acc[v][2] += face_n[2];
            }
        }

        for (i, n_acc) in acc.into_iter().enumerate() {
            self.normals[i] = normalize3(n_acc);
        }
    }
}

impl Default for ParamTriMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Topology helpers
// ---------------------------------------------------------------------------

/// Build a half-edge table for a triangle mesh.
///
/// Returns one entry per directed edge (3 per triangle).
/// Each entry is `[next_he, prev_he, twin_he]` where indices refer to
/// half-edges numbered as `tri_idx * 3 + local_edge`.
/// Twin is `-1` if the edge is on the boundary.
pub fn build_halfedge(triangles: &[[usize; 3]]) -> Vec<[i64; 3]> {
    let n_tri = triangles.len();
    let n_he = n_tri * 3;
    let mut table = vec![[0i64; 3]; n_he];

    // Fill next / prev within each triangle
    for (t, _tri) in triangles.iter().enumerate() {
        let base = t * 3;
        table[base][0] = (base + 1) as i64; // next
        table[base][1] = (base + 2) as i64; // prev
        table[base + 1][0] = (base + 2) as i64;
        table[base + 1][1] = base as i64;
        table[base + 2][0] = base as i64;
        table[base + 2][1] = (base + 1) as i64;
    }

    // Build edge -> half-edge map for twin lookup
    // Edge key: (min_v, max_v, direction: from < to)
    let mut edge_map: HashMap<(usize, usize), i64> = HashMap::new();

    for (t, tri) in triangles.iter().enumerate() {
        for local in 0..3usize {
            let v_from = tri[local];
            let v_to = tri[(local + 1) % 3];
            let he_idx = (t * 3 + local) as i64;
            edge_map.insert((v_from, v_to), he_idx);
        }
    }

    // Set twins
    for (t, tri) in triangles.iter().enumerate() {
        for local in 0..3usize {
            let he_idx = t * 3 + local;
            let v_from = tri[local];
            let v_to = tri[(local + 1) % 3];
            // Twin has opposite direction
            if let Some(&twin_idx) = edge_map.get(&(v_to, v_from)) {
                table[he_idx][2] = twin_idx;
            } else {
                table[he_idx][2] = -1; // boundary
            }
        }
    }

    table
}

/// Return all vertices adjacent to `vertex_idx`.
pub fn vertex_neighbors(vertex_idx: usize, triangles: &[[usize; 3]]) -> Vec<usize> {
    let mut neighbors = Vec::new();
    for tri in triangles {
        for (local, &v) in tri.iter().enumerate() {
            if v == vertex_idx {
                let next = tri[(local + 1) % 3];
                let prev = tri[(local + 2) % 3];
                if !neighbors.contains(&next) {
                    neighbors.push(next);
                }
                if !neighbors.contains(&prev) {
                    neighbors.push(prev);
                }
            }
        }
    }
    neighbors
}

/// Return the set of boundary vertices (those on at least one boundary edge).
pub fn boundary_vertices(n_verts: usize, triangles: &[[usize; 3]]) -> Vec<usize> {
    // An edge is a boundary edge if it has no twin.
    // Build directed-edge count; boundary edges appear exactly once.
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles {
        for local in 0..3usize {
            let v_from = tri[local];
            let v_to = tri[(local + 1) % 3];
            *edge_count.entry((v_from, v_to)).or_insert(0) += 1;
        }
    }

    let mut on_boundary = vec![false; n_verts];
    for (&(v_from, v_to), &cnt) in &edge_count {
        if cnt == 1 {
            // Check if reverse edge exists
            if !edge_count.contains_key(&(v_to, v_from)) {
                on_boundary[v_from] = true;
                on_boundary[v_to] = true;
            }
        }
    }

    on_boundary
        .into_iter()
        .enumerate()
        .filter(|&(_, b)| b)
        .map(|(i, _)| i)
        .collect()
}

/// Return true if vertex `v` lies on the boundary of the mesh.
pub fn is_boundary_vertex(v: usize, triangles: &[[usize; 3]]) -> bool {
    let n_verts = triangles
        .iter()
        .flat_map(|t| t.iter())
        .copied()
        .max()
        .map(|m| m + 1)
        .unwrap_or(0);
    boundary_vertices(n_verts, triangles).contains(&v)
}

// ---------------------------------------------------------------------------
// LSCM Parameterization
// ---------------------------------------------------------------------------

/// Least Squares Conformal Maps parameterization.
pub struct LscmParameterization;

impl LscmParameterization {
    /// Compute UV coordinates for the mesh using a simplified LSCM approach.
    ///
    /// Pins the first two boundary vertices and solves iteratively.
    pub fn compute(mesh: &ParamTriMesh) -> Vec<[f64; 2]> {
        let n = mesh.vertex_count();
        if n == 0 {
            return Vec::new();
        }

        let bv = boundary_vertices(n, &mesh.triangles);

        // Pin two boundary vertices
        let pin0 = if bv.is_empty() { 0 } else { bv[0] };
        let pin1 = if bv.len() < 2 {
            (pin0 + 1).min(n - 1)
        } else {
            bv[1]
        };

        let mut uvs = vec![[0.0f64; 2]; n];
        // Set pinned vertices
        uvs[pin0] = [0.0, 0.0];
        uvs[pin1] = [1.0, 0.0];

        // Iterative Gauss-Seidel relaxation using conformal weights
        for _iter in 0..200 {
            for vi in 0..n {
                if vi == pin0 || vi == pin1 {
                    continue;
                }
                // Collect triangles containing vi
                let mut sum_u = 0.0f64;
                let mut sum_v = 0.0f64;
                let mut weight_total = 0.0f64;

                for tri in &mesh.triangles {
                    let local_opt = tri.iter().position(|&x| x == vi);
                    let local = match local_opt {
                        Some(l) => l,
                        None => continue,
                    };

                    let vj = tri[(local + 1) % 3];
                    let vk = tri[(local + 2) % 3];

                    let pi = mesh.positions[vi];
                    let pj = mesh.positions[vj];
                    let pk = mesh.positions[vk];

                    // Cotangent weights
                    let ej = sub3(pj, pi);
                    let ek = sub3(pk, pi);
                    let area2 = length3(cross3(ej, ek));
                    let w = if area2 > 1e-15 { 1.0 / area2 } else { 0.0 };

                    sum_u += w * (uvs[vj][0] + uvs[vk][0]);
                    sum_v += w * (uvs[vj][1] + uvs[vk][1]);
                    weight_total += 2.0 * w;
                }

                if weight_total > 1e-15 {
                    uvs[vi] = [sum_u / weight_total, sum_v / weight_total];
                }
            }
        }

        uvs
    }

    /// Compute conformal energy: sum of |J - R|^2 per triangle.
    ///
    /// J is the Jacobian of the UV mapping, R is the closest rotation.
    pub fn conformal_energy(
        positions: &[[f64; 3]],
        uvs: &[[f64; 2]],
        triangles: &[[usize; 3]],
    ) -> f64 {
        let mut energy = 0.0f64;

        for tri in triangles {
            let p0 = positions[tri[0]];
            let p1 = positions[tri[1]];
            let p2 = positions[tri[2]];

            let u0 = uvs[tri[0]];
            let u1 = uvs[tri[1]];
            let u2 = uvs[tri[2]];

            // 3D edge vectors
            let e1 = sub3(p1, p0);
            let e2 = sub3(p2, p0);

            // UV edge vectors
            let f1 = sub2(u1, u0);
            let f2 = sub2(u2, u0);

            // Area in 3D
            let cross = cross3(e1, e2);
            let area_3d = length3(cross) * 0.5;
            if area_3d < 1e-15 {
                continue;
            }

            // Jacobian entries (2x2 matrix mapping UV -> 3D)
            // J = [e1, e2] * inv([f1, f2])
            let det_f = cross2(f1, f2);
            if det_f.abs() < 1e-15 {
                energy += area_3d; // degenerate UV triangle counts as full error
                continue;
            }
            let inv_det = 1.0 / det_f;

            // J (3x2): columns are J*e_u, J*e_v
            let ju = [
                (e1[0] * f2[1] - e2[0] * f1[1]) * inv_det,
                (e1[1] * f2[1] - e2[1] * f1[1]) * inv_det,
                (e1[2] * f2[1] - e2[2] * f1[1]) * inv_det,
            ];
            let jv = [
                (e2[0] * f1[0] - e1[0] * f2[0]) * inv_det,
                (e2[1] * f1[0] - e1[1] * f2[0]) * inv_det,
                (e2[2] * f1[0] - e1[2] * f2[0]) * inv_det,
            ];

            // Conformal energy: |Ju|^2 + |Jv|^2 - 2 (should be 0 for conformal)
            let norm_ju = dot3(ju, ju);
            let norm_jv = dot3(jv, jv);
            // LSCM energy = |Ju - perp(Jv)|^2; simplified as |Ju|^2 + |Jv|^2 - 2*area_3d/area_uv
            let area_uv = det_f.abs() * 0.5;
            let stretch = (norm_ju + norm_jv) * area_uv - 2.0 * area_3d;
            energy += stretch.abs();
        }

        energy
    }
}

// ---------------------------------------------------------------------------
// Tutte Parameterization
// ---------------------------------------------------------------------------

/// Tutte (barycentric) parameterization.
pub struct TutteParameterization;

impl TutteParameterization {
    /// Compute UV coordinates using Tutte's method.
    ///
    /// Maps boundary vertices to a circle and solves the Laplacian
    /// system for interior vertices using Gauss-Seidel iteration.
    pub fn compute(mesh: &ParamTriMesh) -> Vec<[f64; 2]> {
        let n = mesh.vertex_count();
        if n == 0 {
            return Vec::new();
        }

        let bv = boundary_vertices(n, &mesh.triangles);
        if bv.is_empty() {
            // No boundary: return zeros
            return vec![[0.0; 2]; n];
        }

        let boundary_uvs = Self::map_boundary_to_circle(&bv, n);

        let mut uvs = vec![[0.0f64; 2]; n];
        // Assign boundary UVs
        for &bvi in bv.iter() {
            uvs[bvi] = boundary_uvs[bvi];
        }

        let is_boundary: Vec<bool> = (0..n).map(|v| bv.contains(&v)).collect();

        // Gauss-Seidel for interior vertices
        for _iter in 0..500 {
            for vi in 0..n {
                if is_boundary[vi] {
                    continue;
                }
                let neighbors = vertex_neighbors(vi, &mesh.triangles);
                if neighbors.is_empty() {
                    continue;
                }
                let k = neighbors.len() as f64;
                let sum_u: f64 = neighbors.iter().map(|&nb| uvs[nb][0]).sum();
                let sum_v: f64 = neighbors.iter().map(|&nb| uvs[nb][1]).sum();
                uvs[vi] = [sum_u / k, sum_v / k];
            }
        }

        uvs
    }

    /// Map boundary vertices equally spaced on the unit circle.
    ///
    /// Returns a Vec of length `n_verts` where only boundary indices are set;
    /// interior entries are `[0.0, 0.0]`.
    pub fn map_boundary_to_circle(boundary: &[usize], n_verts: usize) -> Vec<[f64; 2]> {
        use std::f64::consts::TAU;
        let mut uvs = vec![[0.0f64; 2]; n_verts];
        let nb = boundary.len();
        for (i, &bv) in boundary.iter().enumerate() {
            let angle = TAU * i as f64 / nb as f64;
            uvs[bv] = [angle.cos(), angle.sin()];
        }
        uvs
    }
}

// ---------------------------------------------------------------------------
// UV Quality Metrics
// ---------------------------------------------------------------------------

/// Compute per-triangle texture distortion: ratio of UV area to 3D area.
///
/// A value of 1.0 means no distortion.
pub fn texture_distortion(
    positions: &[[f64; 3]],
    uvs: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> Vec<f64> {
    triangles
        .iter()
        .map(|tri| {
            let p0 = positions[tri[0]];
            let p1 = positions[tri[1]];
            let p2 = positions[tri[2]];

            let u0 = uvs[tri[0]];
            let u1 = uvs[tri[1]];
            let u2 = uvs[tri[2]];

            let e1_3d = sub3(p1, p0);
            let e2_3d = sub3(p2, p0);
            let area_3d = length3(cross3(e1_3d, e2_3d)) * 0.5;

            let e1_uv = sub2(u1, u0);
            let e2_uv = sub2(u2, u0);
            let area_uv = cross2(e1_uv, e2_uv).abs() * 0.5;

            if area_3d < 1e-15 || area_uv < 1e-15 {
                0.0
            } else {
                area_uv / area_3d
            }
        })
        .collect()
}

/// Mean UV distortion across all triangles.
pub fn uv_stretch(positions: &[[f64; 3]], uvs: &[[f64; 2]], triangles: &[[usize; 3]]) -> f64 {
    let distortions = texture_distortion(positions, uvs, triangles);
    if distortions.is_empty() {
        return 0.0;
    }
    distortions.iter().sum::<f64>() / distortions.len() as f64
}

/// Count triangles with inverted UV mapping (negative signed area).
pub fn uv_overlap_check(uvs: &[[f64; 2]], triangles: &[[usize; 3]]) -> usize {
    triangles
        .iter()
        .filter(|tri| {
            let u0 = uvs[tri[0]];
            let u1 = uvs[tri[1]];
            let u2 = uvs[tri[2]];
            let e1 = sub2(u1, u0);
            let e2 = sub2(u2, u0);
            cross2(e1, e2) < 0.0
        })
        .count()
}

// ---------------------------------------------------------------------------
// Mesh Subdivision
// ---------------------------------------------------------------------------

/// Loop subdivision: split each triangle into 4 sub-triangles.
///
/// Interior vertices use 3/8 weight from neighbors; boundary uses 1/8 from ends.
pub fn loop_subdivision(mesh: &ParamTriMesh) -> ParamTriMesh {
    let positions = &mesh.positions;
    let triangles = &mesh.triangles;
    let uvs = &mesh.uvs;
    let n_orig = positions.len();

    // Map edge -> midpoint vertex index
    let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_positions = positions.clone();
    let mut new_uvs = uvs.clone();

    // Ensure uvs has enough entries
    while new_uvs.len() < new_positions.len() {
        new_uvs.push([0.0, 0.0]);
    }

    // Create edge midpoints
    for tri in triangles {
        for local in 0..3usize {
            let v0 = tri[local];
            let v1 = tri[(local + 1) % 3];
            let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
            edge_midpoints.entry(key).or_insert_with(|| {
                let mid_pos = [
                    (positions[v0][0] + positions[v1][0]) * 0.5,
                    (positions[v0][1] + positions[v1][1]) * 0.5,
                    (positions[v0][2] + positions[v1][2]) * 0.5,
                ];
                let mid_uv = if !uvs.is_empty() {
                    [
                        (uvs[v0][0] + uvs[v1][0]) * 0.5,
                        (uvs[v0][1] + uvs[v1][1]) * 0.5,
                    ]
                } else {
                    [0.0, 0.0]
                };
                let idx = new_positions.len();
                new_positions.push(mid_pos);
                new_uvs.push(mid_uv);
                idx
            });
        }
    }

    // Build new triangle list (4 per original)
    let mut new_triangles: Vec<[usize; 3]> = Vec::with_capacity(triangles.len() * 4);
    for tri in triangles {
        let v0 = tri[0];
        let v1 = tri[1];
        let v2 = tri[2];

        let key01 = if v0 < v1 { (v0, v1) } else { (v1, v0) };
        let key12 = if v1 < v2 { (v1, v2) } else { (v2, v1) };
        let key20 = if v2 < v0 { (v2, v0) } else { (v0, v2) };

        let m01 = edge_midpoints[&key01];
        let m12 = edge_midpoints[&key12];
        let m20 = edge_midpoints[&key20];

        new_triangles.push([v0, m01, m20]);
        new_triangles.push([m01, v1, m12]);
        new_triangles.push([m20, m12, v2]);
        new_triangles.push([m01, m12, m20]);
    }

    // Update original vertex positions using Loop weights
    let n_new = new_positions.len();
    let mut smoothed = new_positions.clone();

    let bv_set: std::collections::HashSet<usize> =
        boundary_vertices(n_orig, triangles).into_iter().collect();

    for vi in 0..n_orig {
        let neighbors = vertex_neighbors(vi, triangles);
        let k = neighbors.len();
        if k == 0 {
            continue;
        }

        if bv_set.contains(&vi) {
            // Boundary rule: 3/4 own + 1/8 * (two boundary neighbors)
            let boundary_neighbors: Vec<usize> = neighbors
                .iter()
                .copied()
                .filter(|&nb| bv_set.contains(&nb))
                .collect();
            if boundary_neighbors.len() >= 2 {
                let nb0 = boundary_neighbors[0];
                let nb1 = boundary_neighbors[1];
                smoothed[vi] = [
                    0.75 * new_positions[vi][0]
                        + 0.125 * (new_positions[nb0][0] + new_positions[nb1][0]),
                    0.75 * new_positions[vi][1]
                        + 0.125 * (new_positions[nb0][1] + new_positions[nb1][1]),
                    0.75 * new_positions[vi][2]
                        + 0.125 * (new_positions[nb0][2] + new_positions[nb1][2]),
                ];
            }
        } else {
            // Interior rule: beta weight
            let beta = if k == 3 {
                3.0 / 16.0
            } else {
                3.0 / (8.0 * k as f64)
            };
            let neighbor_sum: [f64; 3] = neighbors.iter().fold([0.0; 3], |acc, &nb| {
                [
                    acc[0] + new_positions[nb][0],
                    acc[1] + new_positions[nb][1],
                    acc[2] + new_positions[nb][2],
                ]
            });
            smoothed[vi] = [
                (1.0 - beta * k as f64) * new_positions[vi][0] + beta * neighbor_sum[0],
                (1.0 - beta * k as f64) * new_positions[vi][1] + beta * neighbor_sum[1],
                (1.0 - beta * k as f64) * new_positions[vi][2] + beta * neighbor_sum[2],
            ];
        }
    }

    // Pad smoothed to full length
    for pos in new_positions.iter().take(n_new).skip(n_orig) {
        smoothed.push(*pos);
    }

    // Build result mesh
    let mut result = ParamTriMesh {
        positions: smoothed,
        uvs: new_uvs,
        normals: vec![[0.0; 3]; n_new],
        triangles: new_triangles,
    };
    result.compute_normals();
    result
}

/// Simple midpoint subdivision: each edge gets a midpoint, 1 triangle -> 4.
pub fn midpoint_subdivision(
    positions: &[[f64; 3]],
    triangles: &[[usize; 3]],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut new_positions = positions.to_vec();
    let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_triangles: Vec<[usize; 3]> = Vec::with_capacity(triangles.len() * 4);

    for tri in triangles {
        let v0 = tri[0];
        let v1 = tri[1];
        let v2 = tri[2];

        let get_or_create =
            |v_a: usize,
             v_b: usize,
             positions: &[[f64; 3]],
             new_positions: &mut Vec<[f64; 3]>,
             edge_midpoints: &mut HashMap<(usize, usize), usize>| {
                let key = if v_a < v_b { (v_a, v_b) } else { (v_b, v_a) };
                if let Some(&idx) = edge_midpoints.get(&key) {
                    idx
                } else {
                    let mid = [
                        (positions[v_a][0] + positions[v_b][0]) * 0.5,
                        (positions[v_a][1] + positions[v_b][1]) * 0.5,
                        (positions[v_a][2] + positions[v_b][2]) * 0.5,
                    ];
                    let idx = new_positions.len();
                    new_positions.push(mid);
                    edge_midpoints.insert(key, idx);
                    idx
                }
            };

        let m01 = get_or_create(v0, v1, positions, &mut new_positions, &mut edge_midpoints);
        let m12 = get_or_create(v1, v2, positions, &mut new_positions, &mut edge_midpoints);
        let m20 = get_or_create(v2, v0, positions, &mut new_positions, &mut edge_midpoints);

        new_triangles.push([v0, m01, m20]);
        new_triangles.push([m01, v1, m12]);
        new_triangles.push([m20, m12, v2]);
        new_triangles.push([m01, m12, m20]);
    }

    (new_positions, new_triangles)
}

// ---------------------------------------------------------------------------
// Mesh Repair
// ---------------------------------------------------------------------------

/// Merge vertices within `tol` distance and return the cleaned mesh.
pub fn remove_duplicate_vertices(
    positions: &[[f64; 3]],
    triangles: &[[usize; 3]],
    tol: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let n = positions.len();
    let mut remap = vec![0usize; n];
    let mut unique: Vec<[f64; 3]> = Vec::new();

    for (i, &p) in positions.iter().enumerate() {
        let found = unique.iter().position(|&q| {
            let dx = p[0] - q[0];
            let dy = p[1] - q[1];
            let dz = p[2] - q[2];
            (dx * dx + dy * dy + dz * dz).sqrt() <= tol
        });
        if let Some(j) = found {
            remap[i] = j;
        } else {
            remap[i] = unique.len();
            unique.push(p);
        }
    }

    let new_triangles: Vec<[usize; 3]> = triangles
        .iter()
        .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
        .filter(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2])
        .collect();

    (unique, new_triangles)
}

/// Fill holes by fan-triangulating boundary loops.
///
/// For each detected boundary loop, creates a fan of triangles from the
/// first boundary vertex to all others in the loop.
pub fn fill_holes(triangles: &mut Vec<[usize; 3]>, n_verts: usize) {
    // Find boundary edges (appear exactly once)
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles.iter() {
        for local in 0..3usize {
            let v0 = tri[local];
            let v1 = tri[(local + 1) % 3];
            *edge_count.entry((v0, v1)).or_insert(0) += 1;
        }
    }

    // Collect boundary edges (no reverse)
    let mut boundary_edges: Vec<(usize, usize)> = Vec::new();
    for (&(v0, v1), &cnt) in &edge_count {
        if cnt == 1 && !edge_count.contains_key(&(v1, v0)) {
            boundary_edges.push((v0, v1));
        }
    }

    if boundary_edges.is_empty() {
        return;
    }

    // Build adjacency from boundary edges to find loops
    let mut next_map: HashMap<usize, usize> = HashMap::new();
    for &(v0, v1) in &boundary_edges {
        next_map.insert(v0, v1);
    }

    let mut visited = vec![false; n_verts];
    for start in 0..n_verts {
        if !next_map.contains_key(&start) || visited[start] {
            continue;
        }

        // Trace loop
        let mut loop_verts = vec![start];
        let mut cur = start;
        loop {
            match next_map.get(&cur) {
                Some(&nxt) if nxt != start && !visited[nxt] => {
                    loop_verts.push(nxt);
                    cur = nxt;
                }
                _ => break,
            }
        }

        for &v in &loop_verts {
            visited[v] = true;
        }

        // Fan-fill from first vertex
        if loop_verts.len() >= 3 {
            let apex = loop_verts[0];
            for i in 1..loop_verts.len() - 1 {
                triangles.push([apex, loop_verts[i], loop_verts[i + 1]]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Harmonic map parameterization
// ---------------------------------------------------------------------------

/// Harmonic map UV parameterization (also known as harmonic embedding).
///
/// Maps boundary vertices to a fixed convex boundary (circle) and solves
/// the discrete harmonic (Laplacian) system with cotangent weights for
/// interior vertices.  This provides a more accurate conformal-ish mapping
/// than uniform Tutte.
pub struct HarmonicMapParameterization;

impl HarmonicMapParameterization {
    /// Compute cotangent-weighted UV coordinates for the mesh.
    ///
    /// Boundary vertices are fixed on a unit circle.
    /// Interior vertices are solved by Gauss-Seidel iteration with
    /// cotangent weights.
    pub fn compute(mesh: &ParamTriMesh) -> Vec<[f64; 2]> {
        let n = mesh.vertex_count();
        if n == 0 {
            return Vec::new();
        }

        let bv = boundary_vertices(n, &mesh.triangles);
        if bv.is_empty() {
            return vec![[0.0; 2]; n];
        }

        // Map boundary to circle.
        let boundary_uvs = TutteParameterization::map_boundary_to_circle(&bv, n);
        let mut uvs = vec![[0.0f64; 2]; n];
        let is_boundary: Vec<bool> = (0..n).map(|v| bv.contains(&v)).collect();

        for &bvi in &bv {
            uvs[bvi] = boundary_uvs[bvi];
        }

        // Precompute cotangent weights.
        // w_ij = cot(alpha_ij) + cot(beta_ij) where alpha, beta are opposite angles.
        // For simplicity use area-based weights (simpler to compute but still SPD).
        let cot_weight = |vi: usize, vj: usize| -> f64 {
            let mut w = 0.0f64;
            for tri in &mesh.triangles {
                // Find if vi and vj are both in this triangle.
                let pos_i = tri.iter().position(|&v| v == vi);
                let pos_j = tri.iter().position(|&v| v == vj);
                if let (Some(pi), Some(pj)) = (pos_i, pos_j) {
                    // The opposite vertex.
                    let pk = 3 - pi - pj; // since 0+1+2=3
                    let vk = tri[pk];
                    let pk_pos = mesh.positions[vk];
                    let pi_pos = mesh.positions[vi];
                    let pj_pos = mesh.positions[vj];
                    let ei = sub3(pi_pos, pk_pos);
                    let ej = sub3(pj_pos, pk_pos);
                    let dot = dot3(ei, ej);
                    let cross_n = length3(cross3(ei, ej));
                    if cross_n > 1e-15 {
                        w += dot / cross_n; // = cot(angle at vk)
                    }
                }
            }
            w.max(0.0) // clamp to avoid negative weights
        };

        // Gauss-Seidel with cotangent weights.
        for _iter in 0..500 {
            for vi in 0..n {
                if is_boundary[vi] {
                    continue;
                }
                let neighbors = vertex_neighbors(vi, &mesh.triangles);
                if neighbors.is_empty() {
                    continue;
                }
                let mut sum_u = 0.0f64;
                let mut sum_v = 0.0f64;
                let mut sum_w = 0.0f64;
                for &vj in &neighbors {
                    let w = cot_weight(vi, vj).max(1e-10);
                    sum_u += w * uvs[vj][0];
                    sum_v += w * uvs[vj][1];
                    sum_w += w;
                }
                if sum_w > 1e-15 {
                    uvs[vi] = [sum_u / sum_w, sum_v / sum_w];
                }
            }
        }

        uvs
    }
}

// ---------------------------------------------------------------------------
// UV seam handling
// ---------------------------------------------------------------------------

/// A UV seam: a set of edges where the UV mapping has a discontinuity.
#[derive(Debug, Clone)]
pub struct UvSeam {
    /// Edge pairs `(v_from, v_to)` that are seam edges.
    pub edges: Vec<(usize, usize)>,
}

impl UvSeam {
    /// Create an empty seam.
    pub fn new() -> Self {
        Self { edges: Vec::new() }
    }

    /// Add a seam edge (undirected).
    pub fn add_edge(&mut self, v0: usize, v1: usize) {
        let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
        if !self.edges.contains(&key) {
            self.edges.push(key);
        }
    }

    /// Check whether a directed edge is a seam edge.
    pub fn is_seam_edge(&self, v0: usize, v1: usize) -> bool {
        let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
        self.edges.contains(&key)
    }

    /// Detect seam edges automatically: boundary edges are always seams;
    /// interior edges with a large UV discontinuity are also seams.
    pub fn detect_seams(
        mesh: &ParamTriMesh,
        uvs: &[[f64; 2]],
        discontinuity_threshold: f64,
    ) -> Self {
        let mut seam = Self::new();
        let bverts = boundary_vertices(mesh.vertex_count(), &mesh.triangles);
        // Mark boundary edges as seams.
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &mesh.triangles {
            for local in 0..3usize {
                let v0 = tri[local];
                let v1 = tri[(local + 1) % 3];
                *edge_count.entry((v0, v1)).or_insert(0) += 1;
            }
        }
        for (&(v0, v1), &cnt) in &edge_count {
            if cnt == 1 && !edge_count.contains_key(&(v1, v0)) {
                seam.add_edge(v0, v1);
            }
        }
        // Also mark edges with large UV discontinuity.
        if uvs.len() >= mesh.vertex_count() {
            for &(v0, v1) in edge_count.keys() {
                let du = uvs[v0][0] - uvs[v1][0];
                let dv = uvs[v0][1] - uvs[v1][1];
                if (du * du + dv * dv).sqrt() > discontinuity_threshold {
                    seam.add_edge(v0, v1);
                }
            }
        }
        let _ = bverts;
        seam
    }
}

impl Default for UvSeam {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Texture atlas chart
// ---------------------------------------------------------------------------

/// A single chart (UV island) for texture atlas packing.
#[derive(Debug, Clone)]
pub struct AtlasChart {
    /// Vertex indices belonging to this chart.
    pub vertices: Vec<usize>,
    /// UV coordinates for these vertices (local to \[0,1\]^2).
    pub uvs: Vec<[f64; 2]>,
    /// Axis-aligned bounding box in UV space: `(min_u, min_v, max_u, max_v)`.
    pub bbox: (f64, f64, f64, f64),
}

impl AtlasChart {
    /// Create a chart from a set of vertex indices and their UV coordinates.
    pub fn new(vertices: Vec<usize>, uvs: Vec<[f64; 2]>) -> Self {
        let bbox = if uvs.is_empty() {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let min_u = uvs.iter().map(|uv| uv[0]).fold(f64::INFINITY, f64::min);
            let min_v = uvs.iter().map(|uv| uv[1]).fold(f64::INFINITY, f64::min);
            let max_u = uvs.iter().map(|uv| uv[0]).fold(f64::NEG_INFINITY, f64::max);
            let max_v = uvs.iter().map(|uv| uv[1]).fold(f64::NEG_INFINITY, f64::max);
            (min_u, min_v, max_u, max_v)
        };
        Self {
            vertices,
            uvs,
            bbox,
        }
    }

    /// Width of the chart in UV space.
    pub fn width(&self) -> f64 {
        self.bbox.2 - self.bbox.0
    }

    /// Height of the chart in UV space.
    pub fn height(&self) -> f64 {
        self.bbox.3 - self.bbox.1
    }

    /// Normalize UVs to the \[0, 1\] range.
    pub fn normalize(&mut self) {
        let (min_u, min_v, max_u, max_v) = self.bbox;
        let du = (max_u - min_u).max(1e-15);
        let dv = (max_v - min_v).max(1e-15);
        for uv in &mut self.uvs {
            uv[0] = (uv[0] - min_u) / du;
            uv[1] = (uv[1] - min_v) / dv;
        }
        self.bbox = (0.0, 0.0, 1.0, 1.0);
    }
}

// ---------------------------------------------------------------------------
// Chart packing (texture atlas generation)
// ---------------------------------------------------------------------------

/// Result of packing multiple charts into a texture atlas.
#[derive(Debug, Clone)]
pub struct AtlasPackResult {
    /// Packed UV coordinates per vertex (indexed by original vertex index).
    pub packed_uvs: Vec<[f64; 2]>,
    /// Atlas width (normalised, always 1.0 in this implementation).
    pub atlas_width: f64,
    /// Atlas height (normalised).
    pub atlas_height: f64,
    /// Per-chart placement: `(offset_u, offset_v, scale)`.
    pub placements: Vec<(f64, f64, f64)>,
}

/// Pack a list of `AtlasChart`s into a texture atlas using a simple row-based packing.
///
/// Charts are normalized, sorted by height, and placed in rows.
pub fn pack_atlas_charts(
    charts: &mut [AtlasChart],
    n_verts: usize,
    padding: f64,
) -> AtlasPackResult {
    for chart in charts.iter_mut() {
        chart.normalize();
    }
    // Sort charts by height descending.
    charts.sort_by(|a, b| {
        b.height()
            .partial_cmp(&a.height())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut packed_uvs = vec![[0.0f64; 2]; n_verts];
    let mut placements = Vec::new();

    let mut cursor_u = 0.0f64;
    let mut cursor_v = 0.0f64;
    let mut row_height = 0.0f64;
    let scale = 1.0f64 / (charts.len() as f64).sqrt().max(1.0);

    for chart in charts.iter() {
        let w = (chart.width() * scale + padding).min(1.0);
        let h = (chart.height() * scale + padding).min(1.0);

        if cursor_u + w > 1.0 {
            cursor_v += row_height + padding;
            cursor_u = 0.0;
            row_height = 0.0;
        }

        // Place the chart.
        let offset_u = cursor_u;
        let offset_v = cursor_v;
        placements.push((offset_u, offset_v, scale));

        for (&vi, &uv) in chart.vertices.iter().zip(chart.uvs.iter()) {
            if vi < packed_uvs.len() {
                packed_uvs[vi] = [
                    offset_u + uv[0] * chart.width() * scale,
                    offset_v + uv[1] * chart.height() * scale,
                ];
            }
        }

        row_height = row_height.max(h);
        cursor_u += w;
    }

    let atlas_height = cursor_v + row_height;

    AtlasPackResult {
        packed_uvs,
        atlas_width: 1.0,
        atlas_height: atlas_height.max(1.0),
        placements,
    }
}

// ---------------------------------------------------------------------------
// Parameterization quality metrics
// ---------------------------------------------------------------------------

/// Compute the angle distortion per triangle.
///
/// Angle distortion is the maximum relative change in angles from 3D to UV space.
/// A value of 0 means perfect conformality.
pub fn angle_distortion(
    positions: &[[f64; 3]],
    uvs: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> Vec<f64> {
    triangles
        .iter()
        .map(|tri| {
            let p = [positions[tri[0]], positions[tri[1]], positions[tri[2]]];
            let u = [uvs[tri[0]], uvs[tri[1]], uvs[tri[2]]];

            // Compute 3D angles.
            let angles_3d = triangle_angles_3d(&p);
            // Compute UV angles.
            let angles_uv = triangle_angles_2d(&u);

            // Max angle difference.
            let mut max_diff = 0.0f64;
            for i in 0..3 {
                let diff = (angles_3d[i] - angles_uv[i]).abs();
                if diff > max_diff {
                    max_diff = diff;
                }
            }
            max_diff
        })
        .collect()
}

fn triangle_angles_3d(p: &[[f64; 3]; 3]) -> [f64; 3] {
    let mut angles = [0.0f64; 3];
    for i in 0..3 {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let a = sub3(p[j], p[i]);
        let b = sub3(p[k], p[i]);
        let la = length3(a);
        let lb = length3(b);
        if la > 1e-15 && lb > 1e-15 {
            let cos_a = (dot3(a, b) / (la * lb)).clamp(-1.0, 1.0);
            angles[i] = cos_a.acos();
        }
    }
    angles
}

fn triangle_angles_2d(u: &[[f64; 2]; 3]) -> [f64; 3] {
    let mut angles = [0.0f64; 3];
    for i in 0..3 {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let a = [u[j][0] - u[i][0], u[j][1] - u[i][1]];
        let b = [u[k][0] - u[i][0], u[k][1] - u[i][1]];
        let la = (a[0] * a[0] + a[1] * a[1]).sqrt();
        let lb = (b[0] * b[0] + b[1] * b[1]).sqrt();
        if la > 1e-15 && lb > 1e-15 {
            let cos_a = ((a[0] * b[0] + a[1] * b[1]) / (la * lb)).clamp(-1.0, 1.0);
            angles[i] = cos_a.acos();
        }
    }
    angles
}

/// Mean angle distortion across all triangles.
pub fn mean_angle_distortion(
    positions: &[[f64; 3]],
    uvs: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> f64 {
    let dists = angle_distortion(positions, uvs, triangles);
    if dists.is_empty() {
        return 0.0;
    }
    dists.iter().sum::<f64>() / dists.len() as f64
}

/// Compute the isometric distortion: ratio of singular values of the
/// Jacobian for each triangle.
///
/// Returns 0 for triangles with zero area.  A value of 1 means isometric (ideal).
pub fn isometric_distortion(
    positions: &[[f64; 3]],
    uvs: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> Vec<f64> {
    triangles
        .iter()
        .map(|tri| {
            let p0 = positions[tri[0]];
            let p1 = positions[tri[1]];
            let p2 = positions[tri[2]];
            let u0 = uvs[tri[0]];
            let u1 = uvs[tri[1]];
            let u2 = uvs[tri[2]];

            let e1_3d = sub3(p1, p0);
            let e2_3d = sub3(p2, p0);
            let e1_uv = [u1[0] - u0[0], u1[1] - u0[1]];
            let e2_uv = [u2[0] - u0[0], u2[1] - u0[1]];

            let area_3d = length3(cross3(e1_3d, e2_3d)) * 0.5;
            let area_uv = (cross2(e1_uv, e2_uv)).abs() * 0.5;

            if area_3d < 1e-15 || area_uv < 1e-15 {
                return 0.0;
            }

            // Stretch metric: sqrt(area_3d / area_uv) (or its inverse)
            let ratio = (area_3d / area_uv).sqrt();
            // Normalize so ideal = 1: return distance from 1.
            (ratio - 1.0).abs()
        })
        .collect()
}

/// Summary quality report for a UV parameterization.
#[derive(Debug, Clone)]
pub struct ParamQualityReport {
    /// Mean UV stretch (area ratio).
    pub mean_stretch: f64,
    /// Number of inverted triangles.
    pub n_inverted: usize,
    /// Mean angle distortion (radians).
    pub mean_angle_distortion: f64,
    /// Mean isometric distortion.
    pub mean_isometric: f64,
    /// Total UV area (should be non-zero for a valid parameterization).
    pub total_uv_area: f64,
}

/// Compute a comprehensive quality report for a UV parameterization.
pub fn parameterization_quality_report(
    positions: &[[f64; 3]],
    uvs: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> ParamQualityReport {
    let mean_stretch = uv_stretch(positions, uvs, triangles);
    let n_inverted = uv_overlap_check(uvs, triangles);
    let mean_ang = mean_angle_distortion(positions, uvs, triangles);
    let iso = isometric_distortion(positions, uvs, triangles);
    let mean_iso = if iso.is_empty() {
        0.0
    } else {
        iso.iter().sum::<f64>() / iso.len() as f64
    };
    let total_uv_area: f64 = triangles
        .iter()
        .map(|tri| {
            let u0 = uvs[tri[0]];
            let u1 = uvs[tri[1]];
            let u2 = uvs[tri[2]];
            let e1 = [u1[0] - u0[0], u1[1] - u0[1]];
            let e2 = [u2[0] - u0[0], u2[1] - u0[1]];
            cross2(e1, e2).abs() * 0.5
        })
        .sum();

    ParamQualityReport {
        mean_stretch,
        n_inverted,
        mean_angle_distortion: mean_ang,
        mean_isometric: mean_iso,
        total_uv_area,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_triangle_mesh() -> ParamTriMesh {
        let mut mesh = ParamTriMesh::new();
        mesh.add_vertex([0.0, 0.0, 0.0], [0.0, 0.0]);
        mesh.add_vertex([1.0, 0.0, 0.0], [1.0, 0.0]);
        mesh.add_vertex([0.0, 1.0, 0.0], [0.0, 1.0]);
        mesh.add_triangle(0, 1, 2);
        mesh
    }

    #[test]
    fn test_trimesh_vertices_and_triangles() {
        let mesh = flat_triangle_mesh();
        assert_eq!(mesh.vertex_count(), 3);
        assert_eq!(mesh.triangle_count(), 1);
        assert_eq!(mesh.positions[0], [0.0, 0.0, 0.0]);
        assert_eq!(mesh.triangles[0], [0, 1, 2]);
    }

    #[test]
    fn test_compute_normals_flat_triangle_points_up() {
        let mut mesh = flat_triangle_mesh();
        mesh.compute_normals();
        // Normal of a triangle in the XY plane should be +Z = [0, 0, 1]
        for i in 0..3 {
            let n = mesh.normals[i];
            assert!(n[2] > 0.99, "normal.z should be ~1.0, got {}", n[2]);
            assert!(n[0].abs() < 1e-10, "normal.x should be ~0, got {}", n[0]);
            assert!(n[1].abs() < 1e-10, "normal.y should be ~0, got {}", n[1]);
        }
    }

    #[test]
    fn test_boundary_vertices_simple_mesh() {
        // A simple strip: two triangles sharing an interior edge.
        // Vertices: 0(0,0,0), 1(1,0,0), 2(0,1,0), 3(1,1,0)
        // Triangles: [0,1,2], [1,3,2]
        let triangles: Vec<[usize; 3]> = vec![[0, 1, 2], [1, 3, 2]];
        let bv = boundary_vertices(4, &triangles);
        // All 4 vertices are on the boundary of this open patch
        assert_eq!(bv.len(), 4, "expected 4 boundary vertices, got {:?}", bv);
        for v in 0..4usize {
            assert!(bv.contains(&v), "vertex {} should be on boundary", v);
        }
    }

    #[test]
    fn test_midpoint_subdivision_one_to_four() {
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triangles = vec![[0usize, 1, 2]];
        let (new_pos, new_tri) = midpoint_subdivision(&positions, &triangles);
        assert_eq!(new_tri.len(), 4, "1 triangle should become 4");
        // 3 original + 3 edge midpoints = 6 vertices
        assert_eq!(new_pos.len(), 6, "should have 6 vertices after subdivision");
    }

    #[test]
    fn test_uv_overlap_check_zero_for_valid_uv() {
        // CCW UV triangle: no overlap
        let uvs = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = vec![[0usize, 1, 2]];
        let overlaps = uv_overlap_check(&uvs, &triangles);
        assert_eq!(overlaps, 0, "no inverted triangles expected");
    }

    #[test]
    fn test_uv_overlap_check_detects_inverted() {
        // CW UV triangle: inverted
        let uvs = vec![[0.0f64, 0.0], [0.0, 1.0], [1.0, 0.0]];
        let triangles = vec![[0usize, 1, 2]];
        let overlaps = uv_overlap_check(&uvs, &triangles);
        assert_eq!(overlaps, 1, "one inverted triangle expected");
    }

    #[test]
    fn test_map_boundary_to_circle() {
        let boundary = vec![0usize, 1, 2, 3];
        let uvs = TutteParameterization::map_boundary_to_circle(&boundary, 4);
        for &vi in &boundary {
            let u = uvs[vi][0];
            let v = uvs[vi][1];
            let r = (u * u + v * v).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-10,
                "vertex {} not on unit circle: r={}",
                vi,
                r
            );
        }
    }

    #[test]
    fn test_texture_distortion_uniform() {
        // A triangle in 3D that matches its UV triangle exactly should have distortion = 1.
        // 3D triangle: right triangle with legs = 1 in XY plane -> area = 0.5
        // UV triangle: same shape -> area = 0.5
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let uvs = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = vec![[0usize, 1, 2]];
        let d = texture_distortion(&positions, &uvs, &triangles);
        assert_eq!(d.len(), 1);
        assert!(
            (d[0] - 1.0).abs() < 1e-10,
            "expected distortion 1.0, got {}",
            d[0]
        );
    }

    // ── Harmonic map parameterization ─────────────────────────────────────────

    #[test]
    fn test_harmonic_map_single_triangle() {
        // Single triangle: 3 boundary vertices -> harmonic should place them on circle.
        let mesh = flat_triangle_mesh();
        let uvs = HarmonicMapParameterization::compute(&mesh);
        assert_eq!(uvs.len(), 3);
    }

    #[test]
    fn test_harmonic_map_grid_mesh() {
        // Build a small 2×2 grid.
        let mut mesh = ParamTriMesh::new();
        for j in 0..3 {
            for i in 0..3 {
                mesh.add_vertex([i as f64, j as f64, 0.0], [0.0, 0.0]);
            }
        }
        for j in 0..2 {
            for i in 0..2 {
                let v00 = j * 3 + i;
                let v10 = j * 3 + i + 1;
                let v01 = (j + 1) * 3 + i;
                let v11 = (j + 1) * 3 + i + 1;
                mesh.add_triangle(v00, v10, v11);
                mesh.add_triangle(v00, v11, v01);
            }
        }
        let uvs = HarmonicMapParameterization::compute(&mesh);
        assert_eq!(uvs.len(), 9);
        // No UV should be NaN.
        for uv in &uvs {
            assert!(
                uv[0].is_finite() && uv[1].is_finite(),
                "UV should be finite"
            );
        }
    }

    // ── UV seam handling ──────────────────────────────────────────────────────

    #[test]
    fn test_uv_seam_add_and_query() {
        let mut seam = UvSeam::new();
        seam.add_edge(0, 1);
        assert!(seam.is_seam_edge(0, 1));
        assert!(seam.is_seam_edge(1, 0)); // undirected
        assert!(!seam.is_seam_edge(0, 2));
    }

    #[test]
    fn test_uv_seam_no_duplicates() {
        let mut seam = UvSeam::new();
        seam.add_edge(0, 1);
        seam.add_edge(1, 0); // same edge, reversed
        assert_eq!(seam.edges.len(), 1, "Should not store duplicate edges");
    }

    #[test]
    fn test_uv_seam_detect_boundary() {
        let mesh = flat_triangle_mesh();
        let uvs = vec![[0.0f64, 0.0]; 3];
        let seam = UvSeam::detect_seams(&mesh, &uvs, 0.5);
        // Boundary edges of a single triangle should be seams.
        assert!(
            seam.edges.len() >= 3,
            "All 3 edges of single tri are boundary seams"
        );
    }

    // ── AtlasChart ────────────────────────────────────────────────────────────

    #[test]
    fn test_atlas_chart_bbox() {
        let chart = AtlasChart::new(vec![0, 1, 2], vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        assert!((chart.width() - 1.0).abs() < 1e-10);
        assert!((chart.height() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_atlas_chart_normalize() {
        let mut chart = AtlasChart::new(vec![0, 1, 2], vec![[2.0, 3.0], [4.0, 3.0], [2.0, 5.0]]);
        chart.normalize();
        // After normalization, UVs should be in [0, 1].
        for uv in &chart.uvs {
            assert!(uv[0] >= -1e-10 && uv[0] <= 1.0 + 1e-10);
            assert!(uv[1] >= -1e-10 && uv[1] <= 1.0 + 1e-10);
        }
    }

    // ── Chart packing ─────────────────────────────────────────────────────────

    #[test]
    fn test_pack_atlas_charts_basic() {
        let mut charts = vec![
            AtlasChart::new(vec![0, 1, 2], vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            AtlasChart::new(vec![3, 4, 5], vec![[0.0, 0.0], [0.5, 0.0], [0.0, 0.5]]),
        ];
        let result = pack_atlas_charts(&mut charts, 6, 0.01);
        assert_eq!(result.packed_uvs.len(), 6);
        assert_eq!(result.placements.len(), 2);
    }

    #[test]
    fn test_pack_atlas_charts_empty() {
        let mut charts: Vec<AtlasChart> = Vec::new();
        let result = pack_atlas_charts(&mut charts, 0, 0.0);
        assert!(result.packed_uvs.is_empty());
    }

    // ── Parameterization quality metrics ─────────────────────────────────────

    #[test]
    fn test_angle_distortion_zero_for_identical() {
        // When the UV and 3D triangle have the same shape, angle distortion = 0.
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let uvs = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = vec![[0usize, 1, 2]];
        let d = angle_distortion(&positions, &uvs, &triangles);
        assert_eq!(d.len(), 1);
        assert!(
            d[0] < 1e-8,
            "angle distortion should be 0 for same shape, got {}",
            d[0]
        );
    }

    #[test]
    fn test_mean_angle_distortion_empty() {
        let d = mean_angle_distortion(&[], &[], &[]);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_isometric_distortion_zero_for_same_shape() {
        // Same 3D and UV triangle: stretch is 1 -> isometric distortion = 0.
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let uvs = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = vec![[0usize, 1, 2]];
        let iso = isometric_distortion(&positions, &uvs, &triangles);
        assert_eq!(iso.len(), 1);
        assert!(
            iso[0] < 1e-8,
            "isometric distortion should be 0, got {}",
            iso[0]
        );
    }

    #[test]
    fn test_parameterization_quality_report() {
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let uvs = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = vec![[0usize, 1, 2]];
        let report = parameterization_quality_report(&positions, &uvs, &triangles);
        assert_eq!(report.n_inverted, 0);
        assert!(report.total_uv_area > 0.0);
        assert!(report.mean_stretch > 0.0);
    }

    #[test]
    fn test_uv_stretch_single_triangle_identity() {
        let positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let uvs = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = vec![[0usize, 1, 2]];
        let s = uv_stretch(&positions, &uvs, &triangles);
        assert!(
            (s - 1.0).abs() < 1e-10,
            "UV stretch for identity mapping should be 1, got {s}"
        );
    }
}
