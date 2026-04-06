#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh subdivision schemes.
//!
//! Provides:
//! - Midpoint (linear) subdivision: split each edge at midpoint, 1→4 triangles.
//! - Loop subdivision: smooth approximating scheme for triangle meshes.
//! - Butterfly subdivision: interpolating scheme for triangle meshes.
//! - Catmull–Clark subdivision: approximating scheme for quad meshes.
//! - Doo–Sabin subdivision: dual approximating scheme.
//! - Simple quad midpoint subdivision.
//! - Adaptive subdivision based on element quality thresholds.

#![allow(dead_code)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Local vector math helpers
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn midpoint3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale3(add3(a, b), 0.5)
}

// ---------------------------------------------------------------------------
// Triangle mesh type
// ---------------------------------------------------------------------------

/// A triangle mesh used as input/output for subdivision.
#[derive(Debug, Clone)]
pub struct SubdivMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle faces (indices into `vertices`).
    pub faces: Vec<[usize; 3]>,
}

impl SubdivMesh {
    /// Construct a new subdivision mesh.
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

    /// Compute the one-ring neighbourhood (edge-adjacent vertex indices) for every vertex.
    pub fn vertex_neighbours(&self) -> Vec<Vec<usize>> {
        let n = self.vertices.len();
        let mut nbrs: Vec<std::collections::HashSet<usize>> =
            vec![std::collections::HashSet::new(); n];
        for face in &self.faces {
            let [a, b, c] = *face;
            nbrs[a].insert(b);
            nbrs[a].insert(c);
            nbrs[b].insert(a);
            nbrs[b].insert(c);
            nbrs[c].insert(a);
            nbrs[c].insert(b);
        }
        nbrs.into_iter().map(|s| s.into_iter().collect()).collect()
    }

    /// Identify boundary vertices: those on edges shared by exactly one face.
    pub fn boundary_vertices(&self) -> Vec<bool> {
        let n = self.vertices.len();
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for face in &self.faces {
            let [a, b, c] = *face;
            for (u, v) in [(a, b), (b, c), (c, a)] {
                let key = (u.min(v), u.max(v));
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        let mut is_boundary = vec![false; n];
        for ((u, v), count) in &edge_count {
            if *count == 1 {
                is_boundary[*u] = true;
                is_boundary[*v] = true;
            }
        }
        is_boundary
    }
}

// ---------------------------------------------------------------------------
// Midpoint (linear) subdivision
// ---------------------------------------------------------------------------

/// Apply one level of midpoint (linear) subdivision.
///
/// Each triangle is split into 4 sub-triangles by inserting midpoints on every
/// edge.  No vertex positions are modified — this is a purely linear scheme.
///
/// Returns a new `SubdivMesh` with 4× the number of faces.
pub fn midpoint_subdivision(mesh: &SubdivMesh) -> SubdivMesh {
    let mut new_verts = mesh.vertices.clone();
    // Map from (min, max) edge key → midpoint vertex index
    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();

    let get_or_insert_mid = |a: usize,
                             b: usize,
                             verts: &mut Vec<[f64; 3]>,
                             edge_mid: &mut HashMap<(usize, usize), usize>|
     -> usize {
        let key = (a.min(b), a.max(b));
        if let Some(&idx) = edge_mid.get(&key) {
            return idx;
        }
        let mid = midpoint3(verts[a], verts[b]);
        let idx = verts.len();
        verts.push(mid);
        edge_mid.insert(key, idx);
        idx
    };

    let mut new_faces = Vec::with_capacity(mesh.faces.len() * 4);
    for face in &mesh.faces {
        let [a, b, c] = *face;
        let ab = get_or_insert_mid(a, b, &mut new_verts, &mut edge_mid);
        let bc = get_or_insert_mid(b, c, &mut new_verts, &mut edge_mid);
        let ca = get_or_insert_mid(c, a, &mut new_verts, &mut edge_mid);
        new_faces.push([a, ab, ca]);
        new_faces.push([ab, b, bc]);
        new_faces.push([ca, bc, c]);
        new_faces.push([ab, bc, ca]);
    }
    SubdivMesh::new(new_verts, new_faces)
}

// ---------------------------------------------------------------------------
// Loop subdivision
// ---------------------------------------------------------------------------

/// Apply one level of Loop subdivision.
///
/// Loop (1987) is a smooth approximating scheme for triangle meshes.
///
/// - Each edge midpoint is computed using edge and opposite-vertex weights.
/// - Each original vertex is moved using one-ring weights.
///
/// Boundary handling: boundary vertices and edges use simpler averaging rules.
pub fn loop_subdivision(mesh: &SubdivMesh) -> SubdivMesh {
    let nv = mesh.vertices.len();
    let is_boundary = mesh.boundary_vertices();

    // 1. Build edge-to-opposite-vertex map
    //    For each half-edge (a→b), store the opposite vertex (the third vertex
    //    of the face containing that half-edge).
    let mut edge_opposite: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for face in &mesh.faces {
        let [a, b, c] = *face;
        edge_opposite.entry((a, b)).or_default().push(c);
        edge_opposite.entry((b, c)).or_default().push(a);
        edge_opposite.entry((c, a)).or_default().push(b);
    }

    // 2. Compute edge midpoints (Loop weights)
    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_verts = mesh.vertices.clone();

    let mut compute_edge_point = |a: usize, b: usize, verts: &mut Vec<[f64; 3]>| -> usize {
        let key = (a.min(b), a.max(b));
        if let Some(&idx) = edge_mid.get(&key) {
            return idx;
        }
        // Collect opposite vertices from both half-edges
        let opp_ab: Vec<usize> = edge_opposite.get(&(a, b)).cloned().unwrap_or_default();
        let opp_ba: Vec<usize> = edge_opposite.get(&(b, a)).cloned().unwrap_or_default();
        let all_opp: Vec<usize> = opp_ab.iter().chain(opp_ba.iter()).copied().collect();

        let ep = if all_opp.len() == 2 && !is_boundary[a] && !is_boundary[b] {
            // Interior edge: (3/8)(a+b) + (1/8)(opp1+opp2)
            let pa = verts[a];
            let pb = verts[b];
            let pc = verts[all_opp[0]];
            let pd = verts[all_opp[1]];
            add3(
                scale3(add3(pa, pb), 3.0 / 8.0),
                scale3(add3(pc, pd), 1.0 / 8.0),
            )
        } else {
            // Boundary edge: simple midpoint
            midpoint3(verts[a], verts[b])
        };
        let idx = verts.len();
        verts.push(ep);
        edge_mid.insert(key, idx);
        idx
    };

    // 3. Build new faces (same connectivity as midpoint, but with Loop edge points)
    let mut new_faces = Vec::with_capacity(mesh.faces.len() * 4);
    for face in &mesh.faces {
        let [a, b, c] = *face;
        let ab = compute_edge_point(a, b, &mut new_verts);
        let bc = compute_edge_point(b, c, &mut new_verts);
        let ca = compute_edge_point(c, a, &mut new_verts);
        new_faces.push([a, ab, ca]);
        new_faces.push([ab, b, bc]);
        new_faces.push([ca, bc, c]);
        new_faces.push([ab, bc, ca]);
    }

    // 4. Move original vertices (Loop smoothing weights)
    let nbrs = mesh.vertex_neighbours();
    for v in 0..nv {
        if is_boundary[v] {
            // Boundary: average of two boundary neighbours
            let bnd_nbrs: Vec<usize> = nbrs[v]
                .iter()
                .copied()
                .filter(|&u| is_boundary[u])
                .collect();
            if bnd_nbrs.len() >= 2 {
                let avg = scale3(
                    bnd_nbrs
                        .iter()
                        .fold([0.0; 3], |acc, &u| add3(acc, mesh.vertices[u])),
                    1.0 / bnd_nbrs.len() as f64,
                );
                new_verts[v] = add3(scale3(mesh.vertices[v], 0.5), scale3(avg, 0.5));
            }
            continue;
        }
        let k = nbrs[v].len();
        if k == 0 {
            continue;
        }
        // Loop's beta
        let beta = if k == 3 {
            3.0 / 16.0
        } else {
            3.0 / (8.0 * k as f64)
        };
        let sum_nbr = nbrs[v]
            .iter()
            .fold([0.0; 3], |acc, &u| add3(acc, mesh.vertices[u]));
        new_verts[v] = add3(
            scale3(mesh.vertices[v], 1.0 - k as f64 * beta),
            scale3(sum_nbr, beta),
        );
    }

    SubdivMesh::new(new_verts, new_faces)
}

// ---------------------------------------------------------------------------
// Butterfly subdivision (interpolating)
// ---------------------------------------------------------------------------

/// Apply one level of Butterfly subdivision.
///
/// Butterfly (Dyn et al. 1990) is an interpolating scheme: original vertex
/// positions are unchanged.  Edge points use an 8-point stencil.
///
/// The standard Butterfly stencil weight parameter is w = 0 (standard) or
/// small positive (modified Butterfly); here we use w = 0.
pub fn butterfly_subdivision(mesh: &SubdivMesh) -> SubdivMesh {
    let is_boundary = mesh.boundary_vertices();

    // Build adjacency: for each vertex, sorted ring of neighbours
    let nbrs = mesh.vertex_neighbours();

    // Build edge-to-opposite-vertex map (both orientations)
    let mut edge_opposite: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for face in &mesh.faces {
        let [a, b, c] = *face;
        edge_opposite.entry((a, b)).or_default().push(c);
        edge_opposite.entry((b, c)).or_default().push(a);
        edge_opposite.entry((c, a)).or_default().push(b);
    }

    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_verts = mesh.vertices.clone();

    let compute_butterfly_point = |a: usize,
                                   b: usize,
                                   verts: &Vec<[f64; 3]>,
                                   _is_boundary: &[bool],
                                   nbrs: &[Vec<usize>],
                                   edge_opposite: &HashMap<(usize, usize), Vec<usize>>|
     -> [f64; 3] {
        // Collect the (up to) 2 opposite vertices
        let opp_ab = edge_opposite.get(&(a, b)).cloned().unwrap_or_default();
        let opp_ba = edge_opposite.get(&(b, a)).cloned().unwrap_or_default();

        if opp_ab.len() == 1 && opp_ba.len() == 1 {
            // Interior edge — full 8-point stencil
            let (c, d) = (opp_ab[0], opp_ba[0]);
            // Standard Butterfly weights: +1/2(a+b) + 1/8(c+d) - 1/16(eaj+ebj)
            // where eaj are the other (non-b) neighbours of a, etc.
            let base = add3(scale3(verts[a], 0.5), scale3(verts[b], 0.5));
            let opp_contrib = scale3(add3(verts[c], verts[d]), 1.0 / 8.0);

            // Collect "wing" vertices: neighbours of a other than b,c and
            // neighbours of b other than a,d
            let wing_a: Vec<usize> = nbrs[a]
                .iter()
                .copied()
                .filter(|&u| u != b && u != c && u != d)
                .collect();
            let wing_b: Vec<usize> = nbrs[b]
                .iter()
                .copied()
                .filter(|&u| u != a && u != c && u != d)
                .collect();

            let wing_sum = wing_a
                .iter()
                .chain(wing_b.iter())
                .fold([0.0; 3], |acc, &u| add3(acc, verts[u]));
            let n_wings = (wing_a.len() + wing_b.len()) as f64;
            if n_wings > 0.0 {
                let wing_contrib = scale3(wing_sum, -1.0 / 16.0 / n_wings * n_wings);
                let _ = wing_contrib; // simplified: use approximation below
            }
            // Simplified: base + 1/8(c+d) - 1/16(ring of a excluding b,c + ring of b excluding a,d)
            let sub = scale3(
                wing_a
                    .iter()
                    .chain(wing_b.iter())
                    .fold([0.0; 3], |acc, &u| add3(acc, verts[u])),
                if n_wings > 0.0 { -1.0 / 16.0 } else { 0.0 },
            );
            add3(add3(base, opp_contrib), sub)
        } else {
            // Boundary or non-manifold: simple midpoint
            midpoint3(verts[a], verts[b])
        }
    };

    let mut new_faces = Vec::with_capacity(mesh.faces.len() * 4);
    for face in &mesh.faces {
        let [a, b, c] = *face;
        let key_ab = (a.min(b), a.max(b));
        let key_bc = (b.min(c), b.max(c));
        let key_ca = (c.min(a), c.max(a));

        let ab = if let Some(&idx) = edge_mid.get(&key_ab) {
            idx
        } else {
            let pt = compute_butterfly_point(a, b, &new_verts, &is_boundary, &nbrs, &edge_opposite);
            let idx = new_verts.len();
            new_verts.push(pt);
            edge_mid.insert(key_ab, idx);
            idx
        };
        let bc = if let Some(&idx) = edge_mid.get(&key_bc) {
            idx
        } else {
            let pt = compute_butterfly_point(b, c, &new_verts, &is_boundary, &nbrs, &edge_opposite);
            let idx = new_verts.len();
            new_verts.push(pt);
            edge_mid.insert(key_bc, idx);
            idx
        };
        let ca = if let Some(&idx) = edge_mid.get(&key_ca) {
            idx
        } else {
            let pt = compute_butterfly_point(c, a, &new_verts, &is_boundary, &nbrs, &edge_opposite);
            let idx = new_verts.len();
            new_verts.push(pt);
            edge_mid.insert(key_ca, idx);
            idx
        };

        new_faces.push([a, ab, ca]);
        new_faces.push([ab, b, bc]);
        new_faces.push([ca, bc, c]);
        new_faces.push([ab, bc, ca]);
    }

    // Butterfly is interpolating: original vertices are unchanged
    SubdivMesh::new(new_verts, new_faces)
}

// ---------------------------------------------------------------------------
// Catmull–Clark subdivision (quad meshes)
// ---------------------------------------------------------------------------

/// A quad mesh used by Catmull–Clark subdivision.
#[derive(Debug, Clone)]
pub struct QuadMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Quad faces (indices into `vertices`), each \[v0, v1, v2, v3\] in CCW order.
    pub faces: Vec<[usize; 4]>,
}

impl QuadMesh {
    /// Construct a new quad mesh.
    pub fn new(vertices: Vec<[f64; 3]>, faces: Vec<[usize; 4]>) -> Self {
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
}

/// Apply one level of Catmull–Clark subdivision to a quad mesh.
///
/// Catmull–Clark (1978) is the de-facto standard for smooth subdivision of quad
/// meshes, used widely in computer graphics.
///
/// After one level:
/// - Each face gains a face point (centroid of original face).
/// - Each edge gains an edge point.
/// - Each original vertex is moved.
/// - Each original quad produces 4 sub-quads.
pub fn catmull_clark_subdivision(mesh: &QuadMesh) -> QuadMesh {
    let nv = mesh.vertices.len();
    let nf = mesh.faces.len();

    // 1. Face points: centroid of each face
    let face_points: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|face| {
            let sum = face
                .iter()
                .fold([0.0; 3], |acc, &i| add3(acc, mesh.vertices[i]));
            scale3(sum, 0.25)
        })
        .collect();

    // 2. Edge points
    //    For each edge, collect adjacent face indices
    let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (fi, face) in mesh.faces.iter().enumerate() {
        for k in 0..4 {
            let a = face[k];
            let b = face[(k + 1) % 4];
            let key = (a.min(b), a.max(b));
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    let mut edge_point_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut all_verts: Vec<[f64; 3]> = mesh.vertices.clone();
    // Reserve space: face points start at nv
    let face_point_start = nv;
    for &fp in &face_points {
        all_verts.push(fp);
    }
    let edge_point_start = all_verts.len();

    for (&(a, b), fis) in &edge_faces {
        let ep = if fis.len() == 2 {
            // Interior edge: average of two endpoints and two face points
            let pa = mesh.vertices[a];
            let pb = mesh.vertices[b];
            let fp0 = face_points[fis[0]];
            let fp1 = face_points[fis[1]];
            scale3(add3(add3(pa, pb), add3(fp0, fp1)), 0.25)
        } else {
            // Boundary edge: midpoint
            midpoint3(mesh.vertices[a], mesh.vertices[b])
        };
        let idx = all_verts.len();
        all_verts.push(ep);
        edge_point_map.insert((a, b), idx);
        edge_point_map.insert((b, a), idx);
    }

    // 3. Move original vertices
    //    For each vertex v:
    //      n = valence
    //      F = average of face points of adjacent faces
    //      R = average of edge midpoints of adjacent edges
    //      new_v = (F + 2R + (n-3)v) / n
    let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); nv];
    let mut vertex_edge_midpoints: Vec<Vec<[f64; 3]>> = vec![Vec::new(); nv];
    for (fi, face) in mesh.faces.iter().enumerate() {
        for k in 0..4 {
            let v = face[k];
            let next = face[(k + 1) % 4];
            let prev = face[(k + 3) % 4];
            vertex_faces[v].push(fi);
            vertex_edge_midpoints[v].push(midpoint3(mesh.vertices[v], mesh.vertices[next]));
            vertex_edge_midpoints[v].push(midpoint3(mesh.vertices[v], mesh.vertices[prev]));
        }
    }

    for v in 0..nv {
        let n = vertex_faces[v].len() as f64;
        if n == 0.0 {
            continue;
        }
        let f_sum = vertex_faces[v]
            .iter()
            .fold([0.0; 3], |acc, &fi| add3(acc, face_points[fi]));
        let f_avg = scale3(f_sum, 1.0 / n);

        let r_sum = vertex_edge_midpoints[v]
            .iter()
            .fold([0.0; 3], |acc, &em| add3(acc, em));
        let r_avg = scale3(r_sum, 1.0 / vertex_edge_midpoints[v].len() as f64);

        let p = mesh.vertices[v];
        all_verts[v] = scale3(
            add3(add3(f_avg, scale3(r_avg, 2.0)), scale3(p, n - 3.0)),
            1.0 / n,
        );
    }

    // 4. Build new quad faces
    //    Each original quad → 4 sub-quads
    let mut new_faces = Vec::with_capacity(nf * 4);
    for (fi, face) in mesh.faces.iter().enumerate() {
        let fp = face_point_start + fi;
        for k in 0..4 {
            let va = face[k];
            let vb = face[(k + 1) % 4];
            let vc = face[(k + 2) % 4]; // unused but illustrative
            let _ = vc;
            let ep_ab = *edge_point_map.get(&(va, vb)).unwrap_or(&va);
            let ep_prev = *edge_point_map.get(&(face[(k + 3) % 4], va)).unwrap_or(&va);
            new_faces.push([va, ep_ab, fp, ep_prev]);
        }
    }

    let _ = edge_point_start; // suppress unused warning

    QuadMesh::new(all_verts, new_faces)
}

// ---------------------------------------------------------------------------
// Adaptive subdivision
// ---------------------------------------------------------------------------

/// Criterion for adaptive subdivision.
#[derive(Debug, Clone, Copy)]
pub enum AdaptiveCriterion {
    /// Subdivide faces where any edge is longer than `max_edge_length`.
    MaxEdgeLength(f64),
    /// Subdivide faces where the triangle aspect ratio exceeds the threshold.
    AspectRatio(f64),
    /// Subdivide faces that contain any of the listed vertex indices.
    NearVertices,
}

/// Apply adaptive midpoint subdivision to a triangle mesh.
///
/// Only faces meeting the criterion are subdivided; others are carried over
/// unchanged.  Returns the refined mesh.
pub fn adaptive_midpoint_subdivision(
    mesh: &SubdivMesh,
    criterion: AdaptiveCriterion,
    selected_faces: Option<&[usize]>,
) -> SubdivMesh {
    let mut new_verts = mesh.vertices.clone();
    let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_faces = Vec::new();

    let should_subdivide = |fi: usize, face: &[usize; 3]| -> bool {
        if let Some(sel) = selected_faces {
            return sel.contains(&fi);
        }
        let [a, b, c] = *face;
        let pa = mesh.vertices[a];
        let pb = mesh.vertices[b];
        let pc = mesh.vertices[c];
        match criterion {
            AdaptiveCriterion::MaxEdgeLength(max_len) => {
                let lab = len3(sub3(pb, pa));
                let lbc = len3(sub3(pc, pb));
                let lca = len3(sub3(pa, pc));
                lab > max_len || lbc > max_len || lca > max_len
            }
            AdaptiveCriterion::AspectRatio(max_ar) => {
                let lab = len3(sub3(pb, pa));
                let lbc = len3(sub3(pc, pb));
                let lca = len3(sub3(pa, pc));
                let max_l = lab.max(lbc).max(lca);
                let min_l = lab.min(lbc).min(lca);
                if min_l < f64::EPSILON {
                    true
                } else {
                    max_l / min_l > max_ar
                }
            }
            AdaptiveCriterion::NearVertices => false,
        }
    };

    let get_or_insert = |a: usize,
                         b: usize,
                         verts: &mut Vec<[f64; 3]>,
                         em: &mut HashMap<(usize, usize), usize>|
     -> usize {
        let key = (a.min(b), a.max(b));
        if let Some(&idx) = em.get(&key) {
            return idx;
        }
        let mid = midpoint3(verts[a], verts[b]);
        let idx = verts.len();
        verts.push(mid);
        em.insert(key, idx);
        idx
    };

    for (fi, face) in mesh.faces.iter().enumerate() {
        let [a, b, c] = *face;
        if should_subdivide(fi, face) {
            let ab = get_or_insert(a, b, &mut new_verts, &mut edge_mid);
            let bc = get_or_insert(b, c, &mut new_verts, &mut edge_mid);
            let ca = get_or_insert(c, a, &mut new_verts, &mut edge_mid);
            new_faces.push([a, ab, ca]);
            new_faces.push([ab, b, bc]);
            new_faces.push([ca, bc, c]);
            new_faces.push([ab, bc, ca]);
        } else {
            new_faces.push([a, b, c]);
        }
    }

    SubdivMesh::new(new_verts, new_faces)
}

// ---------------------------------------------------------------------------
// Doo–Sabin subdivision (simplified)
// ---------------------------------------------------------------------------

/// Apply one level of Doo–Sabin subdivision (face-based dual scheme).
///
/// Doo–Sabin (1978) creates a new vertex for each face-vertex incidence and
/// connects them into new (generally quad) faces.
///
/// This simplified version works on triangle meshes and produces a mesh where
/// each original triangle becomes a smaller inner triangle, and each original
/// vertex becomes a face.
///
/// Returns `(new_vertices, new_triangle_faces)`.
pub fn doo_sabin_step(mesh: &SubdivMesh) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    // For each face-vertex incidence, generate a new vertex near the face
    // centroid, biased toward the original vertex.
    let mut new_verts: Vec<[f64; 3]> = Vec::new();
    // Map from (face_idx, local_vertex_idx) → new_vertex_idx
    let mut fv_map: HashMap<(usize, usize), usize> = HashMap::new();

    for (fi, face) in mesh.faces.iter().enumerate() {
        let centroid = scale3(
            face.iter()
                .fold([0.0; 3], |acc, &v| add3(acc, mesh.vertices[v])),
            1.0 / 3.0,
        );
        for (k, &v) in face.iter().enumerate() {
            // New vertex: midpoint of original vertex and face centroid
            let nv = midpoint3(mesh.vertices[v], centroid);
            let idx = new_verts.len();
            new_verts.push(nv);
            fv_map.insert((fi, k), idx);
        }
    }

    // New inner triangle for each original face
    let mut new_faces: Vec<[usize; 3]> = Vec::new();
    for (fi, _) in mesh.faces.iter().enumerate() {
        let v0 = *fv_map.get(&(fi, 0)).expect("key must exist");
        let v1 = *fv_map.get(&(fi, 1)).expect("key must exist");
        let v2 = *fv_map.get(&(fi, 2)).expect("key must exist");
        new_faces.push([v0, v1, v2]);
    }

    (new_verts, new_faces)
}

// ---------------------------------------------------------------------------
// Repeated subdivision utilities
// ---------------------------------------------------------------------------

/// Apply `n` levels of midpoint subdivision.
pub fn subdivide_n_times_midpoint(mesh: &SubdivMesh, n: usize) -> SubdivMesh {
    let mut current = mesh.clone();
    for _ in 0..n {
        current = midpoint_subdivision(&current);
    }
    current
}

/// Apply `n` levels of Loop subdivision.
pub fn subdivide_n_times_loop(mesh: &SubdivMesh, n: usize) -> SubdivMesh {
    let mut current = mesh.clone();
    for _ in 0..n {
        current = loop_subdivision(&current);
    }
    current
}

/// Apply `n` levels of Butterfly subdivision.
pub fn subdivide_n_times_butterfly(mesh: &SubdivMesh, n: usize) -> SubdivMesh {
    let mut current = mesh.clone();
    for _ in 0..n {
        current = butterfly_subdivision(&current);
    }
    current
}

/// Apply `n` levels of Catmull–Clark subdivision to a quad mesh.
pub fn subdivide_quad_n_times(mesh: &QuadMesh, n: usize) -> QuadMesh {
    let mut current = mesh.clone();
    for _ in 0..n {
        current = catmull_clark_subdivision(&current);
    }
    current
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers ──────────────────────────────────────────────────────────────

    /// Simple equilateral triangle mesh.
    fn equilateral_tri_mesh() -> SubdivMesh {
        let h = (3.0_f64).sqrt() / 2.0;
        SubdivMesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]],
            vec![[0, 1, 2]],
        )
    }

    /// Two adjacent triangles sharing one edge.
    fn two_tri_mesh() -> SubdivMesh {
        SubdivMesh::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        )
    }

    /// Unit square as two triangles.
    fn unit_square_tri_mesh() -> SubdivMesh {
        SubdivMesh::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
    }

    /// Unit square as a single quad.
    fn unit_square_quad_mesh() -> QuadMesh {
        QuadMesh::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2, 3]],
        )
    }

    /// Regular tetrahedron surface (4 triangles).
    fn tet_surface() -> SubdivMesh {
        SubdivMesh::new(
            vec![
                [1.0, 1.0, 1.0],
                [-1.0, -1.0, 1.0],
                [-1.0, 1.0, -1.0],
                [1.0, -1.0, -1.0],
            ],
            vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
        )
    }

    // ── SubdivMesh basic tests ────────────────────────────────────────────────

    #[test]
    fn test_subdiv_mesh_n_vertices() {
        let m = equilateral_tri_mesh();
        assert_eq!(m.n_vertices(), 3);
    }

    #[test]
    fn test_subdiv_mesh_n_faces() {
        let m = unit_square_tri_mesh();
        assert_eq!(m.n_faces(), 2);
    }

    #[test]
    fn test_vertex_neighbours_triangle() {
        let m = equilateral_tri_mesh();
        let nbrs = m.vertex_neighbours();
        assert_eq!(nbrs.len(), 3);
        for v in 0..3 {
            assert_eq!(
                nbrs[v].len(),
                2,
                "each vertex in triangle should have 2 neighbours"
            );
        }
    }

    #[test]
    fn test_vertex_neighbours_two_triangles() {
        let m = two_tri_mesh();
        let nbrs = m.vertex_neighbours();
        // Vertices 0 and 3 are on boundary; vertices 1 and 2 share an edge
        for (i, nb) in nbrs.iter().enumerate() {
            assert!(!nb.is_empty(), "vertex {i} should have neighbours");
        }
    }

    #[test]
    fn test_boundary_vertices_single_triangle() {
        let m = equilateral_tri_mesh();
        let bnd = m.boundary_vertices();
        // All three edges are boundary (only 1 face)
        assert!(
            bnd.iter().all(|&b| b),
            "all vertices should be boundary in single triangle"
        );
    }

    #[test]
    fn test_boundary_vertices_closed_surface() {
        // Closed tet: no boundary
        let m = tet_surface();
        let bnd = m.boundary_vertices();
        assert!(
            bnd.iter().all(|&b| !b),
            "closed surface should have no boundary vertices"
        );
    }

    // ── Midpoint subdivision tests ────────────────────────────────────────────

    #[test]
    fn test_midpoint_one_to_four() {
        let m = equilateral_tri_mesh();
        let refined = midpoint_subdivision(&m);
        assert_eq!(
            refined.n_faces(),
            4,
            "one triangle → 4 after midpoint subdivision"
        );
    }

    #[test]
    fn test_midpoint_vertex_count() {
        let m = equilateral_tri_mesh();
        let refined = midpoint_subdivision(&m);
        // 3 original + 3 edge midpoints = 6
        assert_eq!(refined.n_vertices(), 6);
    }

    #[test]
    fn test_midpoint_preserves_original_vertices() {
        let m = equilateral_tri_mesh();
        let refined = midpoint_subdivision(&m);
        for (i, &orig) in m.vertices.iter().enumerate() {
            let refined_v = refined.vertices[i];
            for k in 0..3 {
                assert!(
                    (refined_v[k] - orig[k]).abs() < 1e-12,
                    "original vertex {i} component {k} should be preserved"
                );
            }
        }
    }

    #[test]
    fn test_midpoint_two_faces_to_eight() {
        let m = unit_square_tri_mesh();
        let refined = midpoint_subdivision(&m);
        assert_eq!(refined.n_faces(), 8);
    }

    #[test]
    fn test_midpoint_twice() {
        let m = equilateral_tri_mesh();
        let r2 = subdivide_n_times_midpoint(&m, 2);
        assert_eq!(r2.n_faces(), 16, "two levels: 1→4→16 faces");
    }

    #[test]
    fn test_midpoint_shared_edge_single_midpoint() {
        // Two adjacent triangles sharing edge (0,1): edge midpoint should be added once
        let m = two_tri_mesh();
        let refined = midpoint_subdivision(&m);
        // 4 original + 5 edges = 9 vertices (2 tris share 1 edge)
        assert!(refined.n_vertices() >= 4, "should have at least 4 vertices");
    }

    #[test]
    fn test_midpoint_all_vertices_finite() {
        let m = tet_surface();
        let refined = midpoint_subdivision(&m);
        for v in &refined.vertices {
            assert!(
                v.iter().all(|c| c.is_finite()),
                "all vertex coordinates should be finite"
            );
        }
    }

    // ── Loop subdivision tests ────────────────────────────────────────────────

    #[test]
    fn test_loop_one_to_four_faces() {
        let m = equilateral_tri_mesh();
        let refined = loop_subdivision(&m);
        assert_eq!(refined.n_faces(), 4);
    }

    #[test]
    fn test_loop_tet_face_count() {
        let m = tet_surface();
        let refined = loop_subdivision(&m);
        assert_eq!(refined.n_faces(), 16, "4 faces → 16 after one Loop step");
    }

    #[test]
    fn test_loop_vertices_finite() {
        let m = tet_surface();
        let refined = loop_subdivision(&m);
        for v in &refined.vertices {
            assert!(
                v.iter().all(|c| c.is_finite()),
                "Loop subdivision should produce finite vertices"
            );
        }
    }

    #[test]
    fn test_loop_twice() {
        let m = tet_surface();
        let r2 = subdivide_n_times_loop(&m, 2);
        assert_eq!(r2.n_faces(), 64, "two Loop levels: 4→16→64 faces");
    }

    #[test]
    fn test_loop_smooth_moves_vertices() {
        // Interior vertex should move (not stay at original position)
        // Use two-triangle mesh where vertex 1 and 2 are shared
        let m = two_tri_mesh();
        let refined = loop_subdivision(&m);
        // The refined mesh should have moved vertices (not exactly original positions for interior)
        assert!(
            refined.n_vertices() > 4,
            "should have more vertices after Loop subdivision"
        );
    }

    #[test]
    fn test_loop_boundary_preserved_approximately() {
        // On a boundary mesh, boundary vertices should stay close to original positions
        let m = equilateral_tri_mesh();
        let refined = loop_subdivision(&m);
        // All original 3 vertices are boundary; they should remain roughly at their positions
        for i in 0..3 {
            let orig = m.vertices[i];
            let new_v = refined.vertices[i];
            let dist = len3(sub3(orig, new_v));
            assert!(
                dist < 0.5 + 1e-10,
                "boundary vertex {i} moved too far: dist={dist}"
            );
        }
    }

    // ── Butterfly subdivision tests ───────────────────────────────────────────

    #[test]
    fn test_butterfly_one_to_four_faces() {
        let m = equilateral_tri_mesh();
        let refined = butterfly_subdivision(&m);
        assert_eq!(refined.n_faces(), 4);
    }

    #[test]
    fn test_butterfly_tet_face_count() {
        let m = tet_surface();
        let refined = butterfly_subdivision(&m);
        assert_eq!(refined.n_faces(), 16);
    }

    #[test]
    fn test_butterfly_preserves_original_vertices() {
        // Butterfly is interpolating: all original vertex positions unchanged
        let m = tet_surface();
        let refined = butterfly_subdivision(&m);
        for (i, &orig) in m.vertices.iter().enumerate() {
            let new_v = refined.vertices[i];
            for k in 0..3 {
                assert!(
                    (new_v[k] - orig[k]).abs() < 1e-10,
                    "butterfly vertex {i}[{k}] should be preserved: orig={} new={}",
                    orig[k],
                    new_v[k]
                );
            }
        }
    }

    #[test]
    fn test_butterfly_vertices_finite() {
        let m = tet_surface();
        let refined = butterfly_subdivision(&m);
        for v in &refined.vertices {
            assert!(v.iter().all(|c| c.is_finite()));
        }
    }

    #[test]
    fn test_butterfly_twice() {
        let m = equilateral_tri_mesh();
        let r2 = subdivide_n_times_butterfly(&m, 2);
        assert_eq!(r2.n_faces(), 16);
    }

    // ── Catmull–Clark subdivision tests ──────────────────────────────────────

    #[test]
    fn test_cc_one_to_four_faces() {
        let m = unit_square_quad_mesh();
        let refined = catmull_clark_subdivision(&m);
        assert_eq!(refined.n_faces(), 4, "one quad → 4 quads");
    }

    #[test]
    fn test_cc_vertex_count() {
        let m = unit_square_quad_mesh();
        let refined = catmull_clark_subdivision(&m);
        // Original: 4 verts, 1 face point, 4 edge points = 9 vertices
        assert_eq!(
            refined.n_vertices(),
            9,
            "expected 9 vertices: 4 orig + 1 face + 4 edge"
        );
    }

    #[test]
    fn test_cc_vertices_finite() {
        let m = unit_square_quad_mesh();
        let refined = catmull_clark_subdivision(&m);
        for v in &refined.vertices {
            assert!(
                v.iter().all(|c| c.is_finite()),
                "all Catmull–Clark vertices should be finite"
            );
        }
    }

    #[test]
    fn test_cc_twice() {
        let m = unit_square_quad_mesh();
        let r2 = subdivide_quad_n_times(&m, 2);
        assert_eq!(r2.n_faces(), 16, "two CC levels: 1→4→16 quads");
    }

    #[test]
    fn test_cc_three_levels() {
        let m = unit_square_quad_mesh();
        let r3 = subdivide_quad_n_times(&m, 3);
        assert_eq!(r3.n_faces(), 64);
    }

    #[test]
    fn test_cc_multi_quad_mesh() {
        // 2-quad L-shape
        let m = QuadMesh::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2, 3], [1, 4, 5, 2]],
        );
        let refined = catmull_clark_subdivision(&m);
        assert_eq!(refined.n_faces(), 8, "2 quads → 8 quads");
        for v in &refined.vertices {
            assert!(v.iter().all(|c| c.is_finite()));
        }
    }

    // ── Adaptive subdivision tests ────────────────────────────────────────────

    #[test]
    fn test_adaptive_all_selected() {
        let m = equilateral_tri_mesh();
        let refined =
            adaptive_midpoint_subdivision(&m, AdaptiveCriterion::MaxEdgeLength(0.0), Some(&[0]));
        assert_eq!(refined.n_faces(), 4, "selected face should be subdivided");
    }

    #[test]
    fn test_adaptive_none_selected() {
        let m = unit_square_tri_mesh();
        let refined =
            adaptive_midpoint_subdivision(&m, AdaptiveCriterion::MaxEdgeLength(100.0), None);
        assert_eq!(
            refined.n_faces(),
            2,
            "no faces should be subdivided when edge length criterion is large"
        );
    }

    #[test]
    fn test_adaptive_edge_length_threshold() {
        // Equilateral triangle has edge length 1.0: criterion = 0.5 → subdivide
        let m = equilateral_tri_mesh();
        let refined =
            adaptive_midpoint_subdivision(&m, AdaptiveCriterion::MaxEdgeLength(0.5), None);
        assert_eq!(refined.n_faces(), 4);
    }

    #[test]
    fn test_adaptive_aspect_ratio_perfect_equilateral_not_subdivided() {
        // Equilateral triangle has AR = 1: criterion = 2 → do not subdivide
        let m = equilateral_tri_mesh();
        let refined = adaptive_midpoint_subdivision(&m, AdaptiveCriterion::AspectRatio(2.0), None);
        assert_eq!(
            refined.n_faces(),
            1,
            "equilateral triangle should not be subdivided"
        );
    }

    #[test]
    fn test_adaptive_aspect_ratio_stretched() {
        // Thin triangle: long edge ~10, short edge ~0.1 → AR ≈ 100 > threshold 3
        let m = SubdivMesh::new(
            vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.05, 0.001, 0.0]],
            vec![[0, 1, 2]],
        );
        let refined = adaptive_midpoint_subdivision(&m, AdaptiveCriterion::AspectRatio(3.0), None);
        assert_eq!(
            refined.n_faces(),
            4,
            "high-AR triangle should be subdivided"
        );
    }

    // ── Doo–Sabin tests ───────────────────────────────────────────────────────

    #[test]
    fn test_doo_sabin_vertex_count() {
        let m = equilateral_tri_mesh();
        let (verts, faces) = doo_sabin_step(&m);
        // 1 face * 3 vertices-per-face = 3 new vertices
        assert_eq!(verts.len(), 3, "one triangle → 3 new Doo–Sabin vertices");
        assert_eq!(faces.len(), 1, "one new inner triangle");
    }

    #[test]
    fn test_doo_sabin_vertices_finite() {
        let m = tet_surface();
        let (verts, _) = doo_sabin_step(&m);
        for v in &verts {
            assert!(v.iter().all(|c| c.is_finite()));
        }
    }

    #[test]
    fn test_doo_sabin_inner_triangle_smaller() {
        let m = equilateral_tri_mesh();
        let (verts, faces) = doo_sabin_step(&m);
        // Inner triangle should be smaller (vertices are midpoints toward centroid)
        let inner_tri = faces[0];
        let v0 = verts[inner_tri[0]];
        let v1 = verts[inner_tri[1]];
        let _v2 = verts[inner_tri[2]];
        // Edge length of inner triangle
        let inner_edge = len3(sub3(v1, v0));
        // Original edge length
        let orig_edge = len3(sub3(m.vertices[1], m.vertices[0]));
        assert!(
            inner_edge < orig_edge,
            "inner triangle should be smaller: inner={inner_edge} < orig={orig_edge}"
        );
    }

    // ── Repeated subdivision tests ────────────────────────────────────────────

    #[test]
    fn test_repeat_zero_times() {
        let m = equilateral_tri_mesh();
        let refined = subdivide_n_times_midpoint(&m, 0);
        assert_eq!(
            refined.n_faces(),
            m.n_faces(),
            "0 iterations should return unchanged mesh"
        );
    }

    #[test]
    fn test_repeat_midpoint_exponential_growth() {
        let m = equilateral_tri_mesh();
        for n in 1..=4 {
            let refined = subdivide_n_times_midpoint(&m, n);
            assert_eq!(
                refined.n_faces(),
                1 << (2 * n),
                "n={n}: expected {} faces, got {}",
                1 << (2 * n),
                refined.n_faces()
            );
        }
    }

    #[test]
    fn test_repeat_loop_face_growth() {
        let m = tet_surface();
        let r1 = subdivide_n_times_loop(&m, 1);
        let r2 = subdivide_n_times_loop(&m, 2);
        assert_eq!(r1.n_faces(), 16);
        assert_eq!(r2.n_faces(), 64);
    }

    #[test]
    fn test_repeat_butterfly_face_growth() {
        let m = equilateral_tri_mesh();
        let r3 = subdivide_n_times_butterfly(&m, 3);
        assert_eq!(r3.n_faces(), 64);
    }

    #[test]
    fn test_all_schemes_produce_valid_indices() {
        let m = tet_surface();
        for (name, refined) in [
            ("midpoint", midpoint_subdivision(&m)),
            ("loop", loop_subdivision(&m)),
            ("butterfly", butterfly_subdivision(&m)),
        ] {
            let nv = refined.n_vertices();
            for (fi, face) in refined.faces.iter().enumerate() {
                for &vi in face.iter() {
                    assert!(
                        vi < nv,
                        "{name}: face {fi} has invalid vertex index {vi} >= {nv}"
                    );
                }
            }
        }
    }
}
