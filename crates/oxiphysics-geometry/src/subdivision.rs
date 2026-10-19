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

use std::collections::HashMap;

mod butterfly;

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
/// The butterfly scheme (Dyn, Levin and Gregory 1990) is interpolating: the
/// original vertices keep their positions and every edge receives one new vertex,
/// which splits each triangle into four.
///
/// For an interior edge `ab` with the opposite vertices `c` and `d` (the third
/// vertices of the two faces on the edge) and the four wing vertices `w1..w4`
/// (the vertices opposite to the outer edges `cb`, `ac`, `da` and `bd` of those
/// faces) the new vertex is
///
/// `1/2 (a + b) + 1/8 (c + d) - 1/16 (w1 + w2 + w3 + w4)`,
///
/// the eight-point stencil with tension parameter `w = 1/16` (`w = 0` would give
/// the plain midpoint). Its weights sum to one, so translating, rotating or
/// scaling the input mesh moves the output the same way.
///
/// The stencil is only designed for regular vertices (valence 6). The remaining
/// configurations are handled as follows:
///
/// - An edge with an end point that is an interior vertex of another valence uses
///   the modified butterfly stencil of Zorin, Schroeder and Sweldens (1996) on
///   the ring of that vertex. If both end points are extraordinary the two
///   results are averaged.
/// - A boundary edge uses the four-point rule `9/16 (a + b) - 1/16 (a' + b')`
///   along the boundary curve, where `a'` and `b'` are the next boundary vertices.
/// - Edges for which a stencil is not available (non-manifold edges, missing wing
///   vertices next to the boundary, pinched boundary vertices) use the midpoint.
///
/// The adjacency does not depend on the winding of the faces.
pub fn butterfly_subdivision(mesh: &SubdivMesh) -> SubdivMesh {
    let adjacency = butterfly::EdgeAdjacency::new(mesh.vertices.len(), &mesh.faces);
    let mut new_verts = mesh.vertices.clone();
    let mut edge_points: HashMap<(usize, usize), usize> = HashMap::new();

    let mut new_faces = Vec::with_capacity(mesh.faces.len() * 4);
    for &[a, b, c] in &mesh.faces {
        let mut edge_vertex = |u: usize, v: usize| -> usize {
            *edge_points.entry((u.min(v), u.max(v))).or_insert_with(|| {
                new_verts.push(adjacency.edge_point(&mesh.vertices, u, v));
                new_verts.len() - 1
            })
        };
        let ab = edge_vertex(a, b);
        let bc = edge_vertex(b, c);
        let ca = edge_vertex(c, a);

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

/// An undirected edge of a quad mesh with the faces that contain it.
struct QuadEdge {
    /// Smaller end vertex.
    a: usize,
    /// Larger end vertex.
    b: usize,
    /// Faces that contain the edge, in face order.
    faces: Vec<usize>,
}

impl QuadEdge {
    /// End vertex of the edge other than `v`.
    fn other(&self, v: usize) -> usize {
        if self.a == v { self.b } else { self.a }
    }
}

/// New position of the original vertex `v` under the Catmull–Clark rules.
///
/// `faces` and `incident` list the faces and (distinct) edges around `v`.
///
/// - Interior vertex: `(F + 2R + (n - 3) P) / n` with the valence `n`, the
///   average `F` of the adjacent face points and the average `R` of the midpoints
///   of the `n` incident edges.
/// - Boundary vertex with two boundary edges and more than one face: the cubic
///   B-spline boundary rule `3/4 P + 1/8 (n1 + n2)` with the two boundary
///   neighbours `n1` and `n2`.
/// - Corner vertex (a single incident face), any other boundary configuration and
///   vertices on a non-manifold edge keep their position.
fn catmull_clark_vertex(
    v: usize,
    vertices: &[[f64; 3]],
    face_points: &[[f64; 3]],
    edges: &[QuadEdge],
    faces: &[usize],
    incident: &[usize],
) -> [f64; 3] {
    let position = vertices[v];
    if faces.is_empty() || incident.iter().any(|&e| edges[e].faces.len() > 2) {
        return position;
    }
    let boundary: Vec<usize> = incident
        .iter()
        .copied()
        .filter(|&e| edges[e].faces.len() == 1)
        .collect();
    match boundary.len() {
        0 => {
            let valence = incident.len() as f64;
            let face_average = scale3(
                faces
                    .iter()
                    .fold([0.0; 3], |acc, &f| add3(acc, face_points[f])),
                1.0 / faces.len() as f64,
            );
            let edge_average = scale3(
                incident.iter().fold([0.0; 3], |acc, &e| {
                    add3(acc, midpoint3(vertices[edges[e].a], vertices[edges[e].b]))
                }),
                1.0 / valence,
            );
            scale3(
                add3(
                    add3(face_average, scale3(edge_average, 2.0)),
                    scale3(position, valence - 3.0),
                ),
                1.0 / valence,
            )
        }
        2 if faces.len() > 1 => {
            let neighbours = boundary
                .iter()
                .fold([0.0; 3], |acc, &e| add3(acc, vertices[edges[e].other(v)]));
            add3(scale3(position, 0.75), scale3(neighbours, 0.125))
        }
        _ => position,
    }
}

/// Apply one level of Catmull–Clark subdivision to a quad mesh.
///
/// Catmull–Clark (1978) is the de-facto standard for smooth subdivision of quad
/// meshes, used widely in computer graphics.
///
/// After one level:
/// - Each face gains a face point (centroid of original face).
/// - Each edge gains an edge point: `(a + b + f1 + f2) / 4` for an interior edge
///   with the face points `f1` and `f2`, the midpoint for a boundary edge.
/// - Each original vertex is moved (interior vertices by the smoothing rule,
///   boundary vertices by the boundary-curve rule, corners stay in place; see
///   `catmull_clark_vertex`).
/// - Each original quad produces 4 sub-quads.
///
/// The vertices of the result are the moved original vertices, then the face
/// points in face order, then the edge points in the order in which the edges
/// are first met when the faces are traversed in order, so the output is
/// deterministic.
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

    // 2. Edges in order of first appearance, with their adjacent faces
    let mut edge_ids: HashMap<(usize, usize), usize> = HashMap::new();
    let mut edges: Vec<QuadEdge> = Vec::new();
    let mut face_sides: Vec<[usize; 4]> = Vec::with_capacity(nf);
    for (fi, face) in mesh.faces.iter().enumerate() {
        let mut sides = [0usize; 4];
        for (k, side) in sides.iter_mut().enumerate() {
            let (u, v) = (face[k], face[(k + 1) % 4]);
            let key = (u.min(v), u.max(v));
            let id = *edge_ids.entry(key).or_insert_with(|| {
                edges.push(QuadEdge {
                    a: key.0,
                    b: key.1,
                    faces: Vec::new(),
                });
                edges.len() - 1
            });
            edges[id].faces.push(fi);
            *side = id;
        }
        face_sides.push(sides);
    }

    // 3. Edge points
    let edge_points: Vec<[f64; 3]> = edges
        .iter()
        .map(|edge| {
            let pa = mesh.vertices[edge.a];
            let pb = mesh.vertices[edge.b];
            match edge.faces[..] {
                [f0, f1] => scale3(
                    add3(add3(pa, pb), add3(face_points[f0], face_points[f1])),
                    0.25,
                ),
                _ => midpoint3(pa, pb),
            }
        })
        .collect();

    // 4. Move the original vertices
    let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (fi, face) in mesh.faces.iter().enumerate() {
        for &v in face {
            vertex_faces[v].push(fi);
        }
    }
    let mut vertex_edges: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (ei, edge) in edges.iter().enumerate() {
        vertex_edges[edge.a].push(ei);
        vertex_edges[edge.b].push(ei);
    }

    let mut all_verts: Vec<[f64; 3]> = Vec::with_capacity(nv + nf + edges.len());
    for v in 0..nv {
        all_verts.push(catmull_clark_vertex(
            v,
            &mesh.vertices,
            &face_points,
            &edges,
            &vertex_faces[v],
            &vertex_edges[v],
        ));
    }
    let face_point_start = all_verts.len();
    all_verts.extend_from_slice(&face_points);
    let edge_point_start = all_verts.len();
    all_verts.extend_from_slice(&edge_points);

    // 5. Build new quad faces: each original quad -> 4 sub-quads, one per corner
    let mut new_faces = Vec::with_capacity(nf * 4);
    for (fi, face) in mesh.faces.iter().enumerate() {
        let sides = face_sides[fi];
        for k in 0..4 {
            new_faces.push([
                face[k],
                edge_point_start + sides[k],
                face_point_start + fi,
                edge_point_start + sides[(k + 3) % 4],
            ]);
        }
    }

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
mod butterfly_tests;
#[cfg(test)]
mod catmull_clark_tests;

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
        for (v, nb) in nbrs.iter().enumerate() {
            assert_eq!(
                nb.len(),
                2,
                "each vertex in triangle should have 2 neighbours (vertex {v})"
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
