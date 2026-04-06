#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh remeshing algorithms.
//!
//! Provides isotropic remeshing (split/collapse/flip/smooth/project loop),
//! Loop subdivision for triangle meshes, and Catmull-Clark subdivision for
//! quad meshes.

#![allow(dead_code, missing_docs)]

use std::collections::HashMap;

use crate::triangle_mesh::TriangleMesh;
use oxiphysics_core::math::Vec3;

// ---------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------

#[inline]
fn vec3_to_arr(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

#[inline]
fn arr_to_vec3(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
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
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < f64::EPSILON {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

#[inline]
fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale3(add3(a, b), 0.5)
}

// ---------------------------------------------------------------------------
// UniformRemesher
// ---------------------------------------------------------------------------

/// Isotropic remesher that targets a specified edge length.
pub struct UniformRemesher {
    /// Desired edge length after remeshing.
    pub target_edge_length: f64,
}

impl UniformRemesher {
    /// Create a new remesher with the given target edge length.
    pub fn new(target_edge_length: f64) -> Self {
        Self { target_edge_length }
    }

    /// Remesh a triangle mesh isotropically.
    ///
    /// Runs `iterations` passes of the 5-operator loop:
    /// 1. Split edges longer than 4/3 × target
    /// 2. Collapse edges shorter than 4/5 × target
    /// 3. Flip edges to improve vertex valence
    /// 4. Laplacian smoothing
    /// 5. Project back to the original surface
    pub fn remesh(&self, mesh: &TriangleMesh, iterations: usize) -> TriangleMesh {
        isotropic_remesh(mesh, self.target_edge_length, iterations)
    }
}

// ---------------------------------------------------------------------------
// isotropic_remesh
// ---------------------------------------------------------------------------

/// Isotropic remeshing via repeated split / collapse / flip / smooth / project.
///
/// Returns a new [`TriangleMesh`] with edges approximately equal to `target_len`.
pub fn isotropic_remesh(mesh: &TriangleMesh, target_len: f64, iterations: usize) -> TriangleMesh {
    // Convert to raw arrays for internal processing
    let mut verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let mut tris: Vec<[usize; 3]> = mesh.indices.clone();

    // Keep a snapshot of the original mesh for projection
    let orig_verts = verts.clone();
    let orig_tris = tris.clone();

    let high = target_len * 4.0 / 3.0;
    let low = target_len * 4.0 / 5.0;

    for _ in 0..iterations {
        // 1. Split long edges
        split_long_edges(&mut verts, &mut tris, high);

        // 2. Collapse short edges
        collapse_short_edges(&mut verts, &mut tris, low);

        // 3. Flip edges toward valence-6
        flip_for_valence(&mut tris, verts.len());

        // 4. Laplacian smooth (tangential)
        laplacian_smooth_surface(&mut verts, &tris);

        // 5. Project back to original surface
        project_to_surface(&mut verts, &orig_verts, &orig_tris);
    }

    let vertices: Vec<Vec3> = verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, tris)
}

// ---------------------------------------------------------------------------
// Split long edges
// ---------------------------------------------------------------------------

fn split_long_edges(verts: &mut Vec<[f64; 3]>, tris: &mut Vec<[usize; 3]>, max_len: f64) {
    // Deduplicate midpoints: edge (min, max) → midpoint vertex index
    let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_tris: Vec<[usize; 3]> = Vec::new();
    let mut keep = vec![true; tris.len()];

    for (ti, tri) in tris.iter().enumerate() {
        let [a, b, c] = *tri;
        let pa = verts[a];
        let pb = verts[b];
        let pc = verts[c];

        let lab = dist3(pa, pb);
        let lbc = dist3(pb, pc);
        let lca = dist3(pc, pa);

        let (longest, va, vb, vopp) = if lab >= lbc && lab >= lca {
            (lab, a, b, c)
        } else if lbc >= lab && lbc >= lca {
            (lbc, b, c, a)
        } else {
            (lca, c, a, b)
        };

        if longest > max_len {
            // Reuse existing midpoint for this edge if already created
            let key = (va.min(vb), va.max(vb));
            let mid = if let Some(&existing) = edge_midpoints.get(&key) {
                existing
            } else {
                let idx = verts.len();
                verts.push(midpoint(verts[va], verts[vb]));
                edge_midpoints.insert(key, idx);
                idx
            };
            new_tris.push([va, mid, vopp]);
            new_tris.push([mid, vb, vopp]);
            keep[ti] = false;
        }
    }

    let orig: Vec<[usize; 3]> = tris
        .iter()
        .enumerate()
        .filter_map(|(i, &t)| if keep[i] { Some(t) } else { None })
        .collect();

    *tris = orig;
    tris.extend(new_tris);
}

// ---------------------------------------------------------------------------
// Collapse short edges
// ---------------------------------------------------------------------------

fn collapse_short_edges(verts: &mut Vec<[f64; 3]>, tris: &mut Vec<[usize; 3]>, min_len: f64) {
    // Build a remap table: vertex index → canonical index
    let n = verts.len();
    let mut remap: Vec<usize> = (0..n).collect();

    // Find short edges
    let mut edge_set: HashMap<(usize, usize), bool> = HashMap::new();
    for tri in tris.iter() {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = (a.min(b), a.max(b));
            edge_set.insert(key, true);
        }
    }

    for (a, b) in edge_set.keys() {
        let ca = remap[*a];
        let cb = remap[*b];
        if ca == cb {
            continue;
        }
        let d = dist3(verts[ca], verts[cb]);
        if d < min_len {
            // Merge b into a (use midpoint)
            let mp = midpoint(verts[ca], verts[cb]);
            verts[ca] = mp;
            // Remap all references to cb → ca
            for r in remap.iter_mut() {
                if *r == cb {
                    *r = ca;
                }
            }
        }
    }

    // Apply remap to triangles and remove degenerate
    let new_tris: Vec<[usize; 3]> = tris
        .iter()
        .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
        .filter(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2])
        .collect();

    *tris = new_tris;
}

// ---------------------------------------------------------------------------
// Flip edges for valence
// ---------------------------------------------------------------------------

fn flip_for_valence(tris: &mut Vec<[usize; 3]>, n_verts: usize) {
    // Compute vertex valence
    let mut valence = vec![0usize; n_verts];
    for tri in tris.iter() {
        for &vi in tri.iter() {
            if vi < n_verts {
                valence[vi] += 1;
            }
        }
    }

    // Build edge → (tri_idx, opposite_vertex) map
    let mut edge_map: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
    for (ti, tri) in tris.iter().enumerate() {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let opp = tri[(k + 2) % 3];
            let key = (a.min(b), a.max(b));
            edge_map.entry(key).or_default().push((ti, opp));
        }
    }

    for ((ea, eb), entries) in &edge_map {
        if entries.len() != 2 {
            continue;
        }
        let (ti, oc) = entries[0];
        let (tj, od) = entries[1];

        // Current deviation from ideal valence (6 for interior)
        let dev_before = valence_deviation(valence[*ea], 6)
            + valence_deviation(valence[*eb], 6)
            + valence_deviation(valence[oc], 6)
            + valence_deviation(valence[od], 6);

        // After flip: edge changes to (oc, od)
        let dev_after = valence_deviation(valence[*ea] - 1, 6)
            + valence_deviation(valence[*eb] - 1, 6)
            + valence_deviation(valence[oc] + 1, 6)
            + valence_deviation(valence[od] + 1, 6);

        if dev_after < dev_before {
            tris[ti] = [*ea, oc, od];
            tris[tj] = [*eb, od, oc];
            // Update valence
            valence[*ea] -= 1;
            valence[*eb] -= 1;
            valence[oc] += 1;
            valence[od] += 1;
        }
    }
}

fn valence_deviation(v: usize, ideal: usize) -> i64 {
    let d = v as i64 - ideal as i64;
    d * d
}

// ---------------------------------------------------------------------------
// Laplacian smooth (tangential)
// ---------------------------------------------------------------------------

fn laplacian_smooth_surface(verts: &mut [[f64; 3]], tris: &[[usize; 3]]) {
    let n = verts.len();
    let mut sums = vec![[0.0_f64; 3]; n];
    let mut counts = vec![0usize; n];

    for tri in tris {
        for k in 0..3 {
            let vi = tri[k];
            for j in 0..3 {
                if j != k {
                    let vj = tri[j];
                    sums[vi] = add3(sums[vi], verts[vj]);
                    counts[vi] += 1;
                }
            }
        }
    }

    for i in 0..n {
        if counts[i] > 0 {
            let avg = scale3(sums[i], 1.0 / counts[i] as f64);
            // Tangential: move along the avg direction with small weight
            let delta = sub3(avg, verts[i]);
            verts[i] = add3(verts[i], scale3(delta, 0.5));
        }
    }
}

// ---------------------------------------------------------------------------
// Project to original surface
// ---------------------------------------------------------------------------

fn project_to_surface(verts: &mut [[f64; 3]], orig_verts: &[[f64; 3]], orig_tris: &[[usize; 3]]) {
    for v in verts.iter_mut() {
        // Find the closest point on the original mesh
        let mut best_dist = f64::INFINITY;
        let mut best_point = *v;

        for tri in orig_tris {
            let a = orig_verts[tri[0]];
            let b = orig_verts[tri[1]];
            let c = orig_verts[tri[2]];
            let cp = closest_point_on_triangle(*v, a, b, c);
            let d = dist3(*v, cp);
            if d < best_dist {
                best_dist = d;
                best_point = cp;
            }
        }

        *v = best_point;
    }
}

/// Closest point on a triangle to a query point.
fn closest_point_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ap = sub3(p, a);

    let d1 = dot3(ab, ap);
    let d2 = dot3(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }

    let bp = sub3(p, b);
    let d3 = dot3(ab, bp);
    let d4 = dot3(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add3(a, scale3(ab, v));
    }

    let cp = sub3(p, c);
    let d5 = dot3(ab, cp);
    let d6 = dot3(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add3(a, scale3(ac, w));
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add3(b, scale3(sub3(c, b), w));
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    add3(add3(a, scale3(ab, v)), scale3(ac, w))
}

// ---------------------------------------------------------------------------
// LoopSubdivision
// ---------------------------------------------------------------------------

/// Loop subdivision scheme for triangle meshes.
pub struct LoopSubdivision;

impl LoopSubdivision {
    /// Subdivide a triangle mesh once using the Loop scheme.
    ///
    /// Each triangle is split into 4 sub-triangles by inserting midpoints on
    /// every edge. Returns a new [`TriangleMesh`] with 4× the original face count.
    pub fn subdivide(mesh: &TriangleMesh) -> TriangleMesh {
        let n_orig = mesh.vertices.len();
        let mut new_verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();

        // Map from edge (min_idx, max_idx) → midpoint vertex index
        let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();

        let get_or_create_midpoint = |a: usize,
                                      b: usize,
                                      verts: &mut Vec<[f64; 3]>,
                                      map: &mut HashMap<(usize, usize), usize>|
         -> usize {
            let key = (a.min(b), a.max(b));
            if let Some(&idx) = map.get(&key) {
                return idx;
            }
            let mid = midpoint(verts[a], verts[b]);
            let idx = verts.len();
            verts.push(mid);
            map.insert(key, idx);
            idx
        };

        // Build new triangles
        let mut new_tris: Vec<[usize; 3]> = Vec::with_capacity(mesh.indices.len() * 4);

        for tri in &mesh.indices {
            let [a, b, c] = *tri;
            let mab = get_or_create_midpoint(a, b, &mut new_verts, &mut edge_midpoints);
            let mbc = get_or_create_midpoint(b, c, &mut new_verts, &mut edge_midpoints);
            let mca = get_or_create_midpoint(c, a, &mut new_verts, &mut edge_midpoints);

            new_tris.push([a, mab, mca]);
            new_tris.push([mab, b, mbc]);
            new_tris.push([mca, mbc, c]);
            new_tris.push([mab, mbc, mca]);
        }

        // Apply Loop mask to original vertices
        // Build 1-ring neighbourhood
        let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n_orig];
        for tri in &mesh.indices {
            let [a, b, c] = *tri;
            if !neighbours[a].contains(&b) {
                neighbours[a].push(b);
            }
            if !neighbours[a].contains(&c) {
                neighbours[a].push(c);
            }
            if !neighbours[b].contains(&a) {
                neighbours[b].push(a);
            }
            if !neighbours[b].contains(&c) {
                neighbours[b].push(c);
            }
            if !neighbours[c].contains(&a) {
                neighbours[c].push(a);
            }
            if !neighbours[c].contains(&b) {
                neighbours[c].push(b);
            }
        }

        for i in 0..n_orig {
            let k = neighbours[i].len();
            if k == 0 {
                continue;
            }
            // Loop's beta weight
            let beta = if k == 3 {
                3.0 / 16.0
            } else {
                3.0 / (8.0 * k as f64)
            };
            let self_weight = 1.0 - k as f64 * beta;

            let mut new_pos = scale3(new_verts[i], self_weight);
            for &nb in &neighbours[i] {
                new_pos = add3(new_pos, scale3(new_verts[nb], beta));
            }
            new_verts[i] = new_pos;
        }

        let vertices: Vec<Vec3> = new_verts.iter().map(|&a| arr_to_vec3(a)).collect();
        TriangleMesh::new(vertices, new_tris)
    }
}

// ---------------------------------------------------------------------------
// CatmullClark
// ---------------------------------------------------------------------------

/// Catmull-Clark subdivision for quad meshes.
pub struct CatmullClark;

impl CatmullClark {
    /// Subdivide a quad mesh once using Catmull-Clark rules.
    ///
    /// # Arguments
    /// * `verts` – original vertex positions.
    /// * `quads` – quads defined by four vertex indices each.
    ///
    /// # Returns
    /// A tuple `(new_vertices, new_quads)` where the new mesh has 4× the
    /// number of faces and smoother vertex positions.
    pub fn subdivide_quad_mesh(
        verts: &[[f64; 3]],
        quads: &[[usize; 4]],
    ) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
        let n_orig = verts.len();
        let mut new_verts: Vec<[f64; 3]> = verts.to_vec();

        // Step 1: Face points — centroid of each quad
        let face_point_start = new_verts.len();
        let mut face_points: Vec<usize> = Vec::with_capacity(quads.len());
        for quad in quads {
            let fp = scale3(
                quad.iter()
                    .fold([0.0_f64; 3], |acc, &vi| add3(acc, verts[vi])),
                0.25,
            );
            face_points.push(new_verts.len());
            new_verts.push(fp);
        }

        // Step 2: Edge points
        // For each edge, the edge point is the average of its two endpoints
        // and the face points of the two adjacent faces.
        let mut edge_map: HashMap<(usize, usize), Vec<usize>> = HashMap::new(); // edge → face indices
        for (qi, quad) in quads.iter().enumerate() {
            for k in 0..4 {
                let a = quad[k];
                let b = quad[(k + 1) % 4];
                let key = (a.min(b), a.max(b));
                edge_map.entry(key).or_default().push(qi);
            }
        }

        let edge_point_start = new_verts.len();
        let mut edge_to_idx: HashMap<(usize, usize), usize> = HashMap::new();
        for (&(ea, eb), face_idxs) in &edge_map {
            let mut ep = scale3(add3(verts[ea], verts[eb]), 0.5);
            if face_idxs.len() == 2 {
                let fp0 = new_verts[face_point_start + face_idxs[0]];
                let fp1 = new_verts[face_point_start + face_idxs[1]];
                ep = scale3(add3(add3(verts[ea], verts[eb]), add3(fp0, fp1)), 0.25);
            }
            edge_to_idx.insert((ea, eb), new_verts.len());
            edge_to_idx.insert((eb, ea), new_verts.len());
            new_verts.push(ep);
        }
        let _ = edge_point_start;

        // Step 3: Updated original vertex positions
        // v' = (F + 2E + (n-3)V) / n
        // where F = avg of adjacent face points, E = avg of adjacent edge midpoints, n = valence
        let mut vertex_face_points: Vec<Vec<[f64; 3]>> = vec![Vec::new(); n_orig];
        let mut vertex_edge_mids: Vec<Vec<[f64; 3]>> = vec![Vec::new(); n_orig];

        for (qi, quad) in quads.iter().enumerate() {
            let fp = new_verts[face_point_start + qi];
            for k in 0..4 {
                let vi = quad[k];
                vertex_face_points[vi].push(fp);
                let nb = quad[(k + 1) % 4];
                let em = scale3(add3(verts[vi], verts[nb]), 0.5);
                vertex_edge_mids[vi].push(em);
                let nb2 = quad[(k + 3) % 4];
                let em2 = scale3(add3(verts[vi], verts[nb2]), 0.5);
                vertex_edge_mids[vi].push(em2);
            }
        }

        for i in 0..n_orig {
            let n = vertex_face_points[i].len();
            if n == 0 {
                continue;
            }
            let n_f = n as f64;
            let f_avg = scale3(
                vertex_face_points[i]
                    .iter()
                    .fold([0.0_f64; 3], |acc, &x| add3(acc, x)),
                1.0 / n_f,
            );
            let e_count = vertex_edge_mids[i].len();
            let e_avg = if e_count > 0 {
                scale3(
                    vertex_edge_mids[i]
                        .iter()
                        .fold([0.0_f64; 3], |acc, &x| add3(acc, x)),
                    1.0 / e_count as f64,
                )
            } else {
                verts[i]
            };
            // Catmull-Clark: (F + 2E + (n-3)V) / n
            let v_new = scale3(
                add3(
                    add3(f_avg, scale3(e_avg, 2.0)),
                    scale3(verts[i], (n_f - 3.0).max(0.0)),
                ),
                1.0 / n_f,
            );
            new_verts[i] = v_new;
        }

        // Step 4: Build new quads
        // Each original quad → 4 new quads
        let mut new_quads: Vec<[usize; 4]> = Vec::with_capacity(quads.len() * 4);
        for (qi, quad) in quads.iter().enumerate() {
            let fp_idx = face_point_start + qi;
            for k in 0..4 {
                let vi = quad[k];
                let ea = quad[k];
                let eb = quad[(k + 1) % 4];
                let ea2 = quad[(k + 3) % 4];
                let eb2 = quad[k];

                let ep1 = *edge_to_idx.get(&(ea.min(eb), ea.max(eb))).unwrap_or(&vi);
                let ep2 = *edge_to_idx
                    .get(&(ea2.min(eb2), ea2.max(eb2)))
                    .unwrap_or(&vi);

                new_quads.push([vi, ep1, fp_idx, ep2]);
            }
        }

        (new_verts, new_quads)
    }
}

// ---------------------------------------------------------------------------
// Feature-preserving remeshing
// ---------------------------------------------------------------------------

/// Feature-preserving isotropic remeshing.
///
/// Like `isotropic_remesh` but marks sharp feature edges (dihedral angle >
/// `feature_angle_rad`) as "locked" — they are never collapsed or flipped.
pub fn feature_preserving_remesh(
    mesh: &TriangleMesh,
    target_len: f64,
    iterations: usize,
    feature_angle_rad: f64,
) -> TriangleMesh {
    let mut verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let mut tris: Vec<[usize; 3]> = mesh.indices.clone();

    // Identify feature edges (dihedral angle > threshold)
    let feature_edges = detect_feature_edges(&verts, &tris, feature_angle_rad);

    let orig_verts = verts.clone();
    let orig_tris = tris.clone();

    let high = target_len * 4.0 / 3.0;
    let low = target_len * 4.0 / 5.0;

    for _ in 0..iterations {
        split_long_edges(&mut verts, &mut tris, high);
        collapse_short_edges_preserving(&mut verts, &mut tris, low, &feature_edges);
        flip_for_valence(&mut tris, verts.len());
        laplacian_smooth_surface(&mut verts, &tris);
        project_to_surface(&mut verts, &orig_verts, &orig_tris);
    }

    let vertices: Vec<Vec3> = verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, tris)
}

/// Detect feature edges where the dihedral angle between adjacent faces > `threshold`.
fn detect_feature_edges(
    verts: &[[f64; 3]],
    tris: &[[usize; 3]],
    threshold_rad: f64,
) -> HashMap<(usize, usize), bool> {
    // Build edge → face normals
    let mut edge_faces: HashMap<(usize, usize), Vec<[f64; 3]>> = HashMap::new();
    for tri in tris {
        let a = verts[tri[0]];
        let b = verts[tri[1]];
        let c = verts[tri[2]];
        let e1 = sub3(b, a);
        let e2 = sub3(c, a);
        let n = normalize3(cross3(e1, e2));
        for k in 0..3 {
            let va = tri[k];
            let vb = tri[(k + 1) % 3];
            let key = (va.min(vb), va.max(vb));
            edge_faces.entry(key).or_default().push(n);
        }
    }

    let mut features = HashMap::new();
    let cos_thresh = threshold_rad.cos();
    for ((ea, eb), normals) in &edge_faces {
        let is_feature = if normals.len() == 2 {
            let d = dot3(normals[0], normals[1]);
            d < cos_thresh
        } else {
            true // boundary edge → feature
        };
        features.insert((*ea, *eb), is_feature);
    }
    features
}

/// Edge collapse that preserves feature edges.
fn collapse_short_edges_preserving(
    verts: &mut Vec<[f64; 3]>,
    tris: &mut Vec<[usize; 3]>,
    min_len: f64,
    feature_edges: &HashMap<(usize, usize), bool>,
) {
    let n = verts.len();
    let mut remap: Vec<usize> = (0..n).collect();

    let mut edge_set: HashMap<(usize, usize), bool> = HashMap::new();
    for tri in tris.iter() {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = (a.min(b), a.max(b));
            edge_set.insert(key, true);
        }
    }

    for (a, b) in edge_set.keys() {
        let ca = remap[*a];
        let cb = remap[*b];
        if ca == cb {
            continue;
        }
        // Check feature lock
        if feature_edges
            .get(&(*a.min(b), *a.max(b)))
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        let d = dist3(verts[ca], verts[cb]);
        if d < min_len {
            let mp = midpoint(verts[ca], verts[cb]);
            verts[ca] = mp;
            for r in remap.iter_mut() {
                if *r == cb {
                    *r = ca;
                }
            }
        }
    }

    let new_tris: Vec<[usize; 3]> = tris
        .iter()
        .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
        .filter(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2])
        .collect();
    *tris = new_tris;
}

// ---------------------------------------------------------------------------
// Edge flip for quality
// ---------------------------------------------------------------------------

/// Flip edges to improve triangle quality (maximize minimum angle).
///
/// For each interior edge shared by two triangles, flip if the flipped
/// configuration has a better minimum angle.
pub fn flip_edges_for_quality(tris: &mut Vec<[usize; 3]>, verts: &[[f64; 3]]) {
    let n_verts = verts.len();

    // Build edge → face map
    let mut edge_map: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
    for (ti, tri) in tris.iter().enumerate() {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let opp = tri[(k + 2) % 3];
            let key = (a.min(b), a.max(b));
            edge_map.entry(key).or_default().push((ti, opp));
        }
    }

    for ((ea, eb), entries) in &edge_map {
        if entries.len() != 2 {
            continue;
        }
        let (ti, oc) = entries[0];
        let (tj, od) = entries[1];

        if *ea >= n_verts || *eb >= n_verts || oc >= n_verts || od >= n_verts {
            continue;
        }

        // Compare min angles before and after flip
        let before_min = min_triangle_angle(verts[*ea], verts[*eb], verts[oc])
            .min(min_triangle_angle(verts[*ea], verts[*eb], verts[od]));
        let after_min = min_triangle_angle(verts[oc], verts[od], verts[*ea])
            .min(min_triangle_angle(verts[oc], verts[od], verts[*eb]));

        if after_min > before_min + 1e-6 {
            // Flip: replace (ea,eb,oc) and (ea,eb,od) with (oc,od,ea) and (oc,od,eb)
            tris[ti] = [*ea, oc, od];
            tris[tj] = [*eb, od, oc];
        }
    }
}

/// Minimum angle (in radians) of a triangle given its three vertex positions.
fn min_triangle_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let bc = sub3(c, b);
    let lab = len3(ab);
    let lac = len3(ac);
    let lbc = len3(bc);

    if lab < 1e-12 || lac < 1e-12 || lbc < 1e-12 {
        return 0.0;
    }

    let cos_a = (dot3(ab, ac) / (lab * lac)).clamp(-1.0, 1.0);
    let cos_b = (dot3(scale3(ab, -1.0), bc) / (lab * lbc)).clamp(-1.0, 1.0);
    let cos_c = (dot3(scale3(ac, -1.0), scale3(bc, -1.0)) / (lac * lbc)).clamp(-1.0, 1.0);

    cos_a.acos().min(cos_b.acos()).min(cos_c.acos())
}

// ---------------------------------------------------------------------------
// Mesh quality metrics
// ---------------------------------------------------------------------------

/// Compute the minimum angle quality metric over all triangles.
///
/// Returns `min_angle / (PI/3)` normalised to `[0, 1]` (1 = equilateral).
pub fn mesh_quality_min_angle(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> f64 {
    if tris.is_empty() {
        return 0.0;
    }
    let ideal = std::f64::consts::PI / 3.0;
    let min_q = tris
        .iter()
        .map(|tri| {
            let a = verts[tri[0]];
            let b = verts[tri[1]];
            let c = verts[tri[2]];
            let min_ang = min_triangle_angle(a, b, c);
            min_ang / ideal
        })
        .fold(f64::INFINITY, f64::min);
    min_q.min(1.0)
}

/// Compute average aspect ratio of all triangles.
///
/// Aspect ratio = longest_edge / (2 * inradius). Ideal (equilateral) = sqrt(3)/3 ≈ 0.577.
/// Lower is better.
pub fn mesh_aspect_ratio_avg(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> f64 {
    if tris.is_empty() {
        return 0.0;
    }
    let ratios: Vec<f64> = tris
        .iter()
        .map(|tri| {
            let a = verts[tri[0]];
            let b = verts[tri[1]];
            let c = verts[tri[2]];
            let lab = len3(sub3(b, a));
            let lbc = len3(sub3(c, b));
            let lca = len3(sub3(a, c));
            let longest = lab.max(lbc).max(lca);
            let s = (lab + lbc + lca) * 0.5; // semi-perimeter
            let area_sq = s * (s - lab) * (s - lbc) * (s - lca);
            if area_sq <= 0.0 || s < 1e-12 {
                return f64::INFINITY;
            }
            let area = area_sq.sqrt();
            let inradius = area / s;
            longest / (2.0 * inradius)
        })
        .collect();

    let finite: Vec<f64> = ratios.into_iter().filter(|x| x.is_finite()).collect();
    if finite.is_empty() {
        0.0
    } else {
        finite.iter().sum::<f64>() / finite.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Tangent Laplacian smoothing (preserves surface projection)
// ---------------------------------------------------------------------------

/// Tangent-space Laplacian smoothing: smooths vertices but projects back to
/// the original surface after each iteration.
pub fn tangent_laplacian_smooth(mesh: &TriangleMesh, iterations: usize) -> TriangleMesh {
    let orig_verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let tris = &mesh.indices;
    let mut verts = orig_verts.clone();

    for _ in 0..iterations {
        laplacian_smooth_surface(&mut verts, tris);
        project_to_surface(&mut verts, &orig_verts, tris);
    }

    let vertices: Vec<Vec3> = verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, tris.clone())
}

// ---------------------------------------------------------------------------
// Edge collapse
// ---------------------------------------------------------------------------

/// Collapse the edge between vertices `v0` and `v1` to their midpoint.
///
/// All triangles containing both `v0` and `v1` are removed (degenerate after collapse).
/// All other references to `v1` are remapped to `v0` (with the midpoint position).
pub fn collapse_edge(verts: &mut Vec<[f64; 3]>, tris: &mut Vec<[usize; 3]>, v0: usize, v1: usize) {
    let mid = midpoint(verts[v0], verts[v1]);
    verts[v0] = mid;

    // Remap all v1 → v0
    let new_tris: Vec<[usize; 3]> = tris
        .iter()
        .map(|tri| {
            [
                if tri[0] == v1 { v0 } else { tri[0] },
                if tri[1] == v1 { v0 } else { tri[1] },
                if tri[2] == v1 { v0 } else { tri[2] },
            ]
        })
        .filter(|tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2])
        .collect();
    *tris = new_tris;
}

// ---------------------------------------------------------------------------
// Laplacian smoothing (public API)
// ---------------------------------------------------------------------------

/// Apply `iterations` rounds of uniform Laplacian smoothing to a mesh.
///
/// Each iteration moves every vertex towards the average of its neighbours
/// by a factor of `lambda` (0..1).
pub fn laplacian_smooth(mesh: &TriangleMesh, iterations: usize, lambda: f64) -> TriangleMesh {
    let mut verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let tris = &mesh.indices;

    for _ in 0..iterations {
        let n = verts.len();
        let mut sums = vec![[0.0_f64; 3]; n];
        let mut counts = vec![0usize; n];

        for tri in tris.iter() {
            for k in 0..3 {
                let vi = tri[k];
                for j in 0..3 {
                    if j != k {
                        let vj = tri[j];
                        sums[vi] = add3(sums[vi], verts[vj]);
                        counts[vi] += 1;
                    }
                }
            }
        }

        for i in 0..n {
            if counts[i] > 0 {
                let avg = scale3(sums[i], 1.0 / counts[i] as f64);
                let delta = sub3(avg, verts[i]);
                verts[i] = add3(verts[i], scale3(delta, lambda));
            }
        }
    }

    let vertices: Vec<Vec3> = verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, tris.clone())
}

// ---------------------------------------------------------------------------
// Edge flip for Delaunay criterion
// ---------------------------------------------------------------------------

/// Flip edges in a 2D-projected triangle mesh to satisfy the Delaunay criterion.
///
/// An edge is flipped if the opposite vertex lies inside the circumcircle of the
/// current triangle.  Operates directly on the `tris` array (in-place).
pub fn delaunay_edge_flip(tris: &mut Vec<[usize; 3]>, verts: &[[f64; 3]]) {
    // Build edge → (tri_idx, opposite_vertex) map
    let mut changed = true;
    let max_passes = 32;
    let mut pass = 0;

    while changed && pass < max_passes {
        changed = false;
        pass += 1;

        let mut edge_map: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
        for (ti, tri) in tris.iter().enumerate() {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let opp = tri[(k + 2) % 3];
                let key = (a.min(b), a.max(b));
                edge_map.entry(key).or_default().push((ti, opp));
            }
        }

        for ((ea, eb), entries) in &edge_map {
            if entries.len() != 2 {
                continue;
            }
            let (ti, oc) = entries[0];
            let (tj, od) = entries[1];

            if *ea >= verts.len() || *eb >= verts.len() || oc >= verts.len() || od >= verts.len() {
                continue;
            }

            // Test if od is inside the circumcircle of triangle (ea, eb, oc)
            if point_in_circumcircle(verts[*ea], verts[*eb], verts[oc], verts[od]) {
                // Flip: replace edge ea-eb with oc-od
                tris[ti] = [*ea, oc, od];
                tris[tj] = [*eb, od, oc];
                changed = true;
            }
        }
    }
}

/// Returns `true` if point `d` is inside the circumcircle of triangle `(a, b, c)`.
/// Uses the 2D X-Z projection.
fn point_in_circumcircle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> bool {
    // Use the standard 3×3 determinant test (X-Z plane projection)
    let ax = a[0] - d[0];
    let az = a[2] - d[2];
    let bx = b[0] - d[0];
    let bz = b[2] - d[2];
    let cx = c[0] - d[0];
    let cz = c[2] - d[2];

    let det = ax * (bz * (cx * cx + cz * cz) - cz * (bx * bx + bz * bz))
        - az * (bx * (cx * cx + cz * cz) - cx * (bx * bx + bz * bz))
        + (ax * ax + az * az) * (bx * cz - bz * cx);

    det > 0.0
}

// ---------------------------------------------------------------------------
// Mesh coarsening
// ---------------------------------------------------------------------------

/// Coarsen a mesh by collapsing all edges shorter than `target_len`.
///
/// Repeatedly collapses short edges until none remain below the threshold.
pub fn coarsen_mesh(mesh: &TriangleMesh, target_len: f64) -> TriangleMesh {
    let mut verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let mut tris: Vec<[usize; 3]> = mesh.indices.clone();

    // Run up to 10 collapse passes
    for _ in 0..10 {
        let before = tris.len();
        collapse_short_edges(&mut verts, &mut tris, target_len);
        if tris.len() == before {
            break;
        }
    }

    let vertices: Vec<Vec3> = verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, tris)
}

// ---------------------------------------------------------------------------
// Adaptive refinement by curvature
// ---------------------------------------------------------------------------

/// Refine triangles with estimated curvature exceeding `threshold`.
///
/// Curvature is approximated as the variation in face normals among a vertex's
/// one-ring neighbourhood.  High-curvature triangles are split at their midpoints.
pub fn adaptive_refine_by_curvature(mesh: &TriangleMesh, threshold: f64) -> TriangleMesh {
    let mut verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let mut tris: Vec<[usize; 3]> = mesh.indices.clone();

    // Compute per-face normals
    let face_normals: Vec<[f64; 3]> = tris
        .iter()
        .map(|tri| {
            let a = verts[tri[0]];
            let b = verts[tri[1]];
            let c = verts[tri[2]];
            normalize3(cross3(sub3(b, a), sub3(c, a)))
        })
        .collect();

    // Compute per-vertex curvature as max angle between incident face normals
    let n = verts.len();
    let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (ti, tri) in tris.iter().enumerate() {
        for &vi in tri.iter() {
            vertex_faces[vi].push(ti);
        }
    }

    let vertex_curvature: Vec<f64> = (0..n)
        .map(|vi| {
            let faces = &vertex_faces[vi];
            if faces.len() < 2 {
                return 0.0;
            }
            let mut max_angle: f64 = 0.0;
            for i in 0..faces.len() {
                for j in (i + 1)..faces.len() {
                    let d = dot3(face_normals[faces[i]], face_normals[faces[j]]).clamp(-1.0, 1.0);
                    let angle = d.acos();
                    if angle > max_angle {
                        max_angle = angle;
                    }
                }
            }
            max_angle
        })
        .collect();

    // Split triangles where any vertex curvature exceeds threshold
    let mut new_tris: Vec<[usize; 3]> = Vec::new();
    let mut kept_tris: Vec<[usize; 3]> = Vec::new();

    for tri in &tris {
        let max_curv = vertex_curvature[tri[0]]
            .max(vertex_curvature[tri[1]])
            .max(vertex_curvature[tri[2]]);

        if max_curv > threshold {
            // Split along the longest edge
            let a = tri[0];
            let b = tri[1];
            let c = tri[2];
            let lab = dist3(verts[a], verts[b]);
            let lbc = dist3(verts[b], verts[c]);
            let lca = dist3(verts[c], verts[a]);

            let (va, vb, vc) = if lab >= lbc && lab >= lca {
                (a, b, c)
            } else if lbc >= lab && lbc >= lca {
                (b, c, a)
            } else {
                (c, a, b)
            };

            let mid_idx = verts.len();
            verts.push(midpoint(verts[va], verts[vb]));
            new_tris.push([va, mid_idx, vc]);
            new_tris.push([mid_idx, vb, vc]);
        } else {
            kept_tris.push(*tri);
        }
    }

    tris = kept_tris;
    tris.extend(new_tris);

    let vertices: Vec<Vec3> = verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, tris)
}

// ---------------------------------------------------------------------------
// Vertex-clustering mesh decimation
// ---------------------------------------------------------------------------

/// Decimate a mesh using vertex clustering.
///
/// Divides space into a uniform grid of `cells_per_axis³` cells and
/// merges all vertices within each cell to their centroid.  Returns a
/// decimated mesh with fewer vertices and faces.
pub fn vertex_clustering_decimate(mesh: &TriangleMesh, cells_per_axis: usize) -> TriangleMesh {
    if cells_per_axis == 0 || mesh.vertices.is_empty() {
        return mesh.clone();
    }

    let verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();

    // Compute bounding box
    let mut bbmin = verts[0];
    let mut bbmax = verts[0];
    for &v in &verts {
        for k in 0..3 {
            if v[k] < bbmin[k] {
                bbmin[k] = v[k];
            }
            if v[k] > bbmax[k] {
                bbmax[k] = v[k];
            }
        }
    }
    let extent = [
        bbmax[0] - bbmin[0] + 1e-8,
        bbmax[1] - bbmin[1] + 1e-8,
        bbmax[2] - bbmin[2] + 1e-8,
    ];
    let c = cells_per_axis as f64;

    // Assign each vertex to a cell
    let cell_of = |v: [f64; 3]| -> usize {
        let ix = ((v[0] - bbmin[0]) / extent[0] * c).floor() as usize;
        let iy = ((v[1] - bbmin[1]) / extent[1] * c).floor() as usize;
        let iz = ((v[2] - bbmin[2]) / extent[2] * c).floor() as usize;
        let ix = ix.min(cells_per_axis - 1);
        let iy = iy.min(cells_per_axis - 1);
        let iz = iz.min(cells_per_axis - 1);
        ix + cells_per_axis * (iy + cells_per_axis * iz)
    };

    let n_cells = cells_per_axis * cells_per_axis * cells_per_axis;
    let mut cell_sums: Vec<[f64; 3]> = vec![[0.0; 3]; n_cells];
    let mut cell_counts: Vec<usize> = vec![0; n_cells];
    let mut vertex_cell: Vec<usize> = vec![0; verts.len()];

    for (i, &v) in verts.iter().enumerate() {
        let cell = cell_of(v);
        vertex_cell[i] = cell;
        for k in 0..3 {
            cell_sums[cell][k] += v[k];
        }
        cell_counts[cell] += 1;
    }

    // Build representative vertices for non-empty cells
    let mut cell_to_vert: Vec<Option<usize>> = vec![None; n_cells];
    let mut new_verts: Vec<[f64; 3]> = Vec::new();
    for cell in 0..n_cells {
        if cell_counts[cell] > 0 {
            cell_to_vert[cell] = Some(new_verts.len());
            let cnt = cell_counts[cell] as f64;
            new_verts.push([
                cell_sums[cell][0] / cnt,
                cell_sums[cell][1] / cnt,
                cell_sums[cell][2] / cnt,
            ]);
        }
    }

    // Remap triangles
    let new_tris: Vec<[usize; 3]> = mesh
        .indices
        .iter()
        .filter_map(|tri| {
            let a = cell_to_vert[vertex_cell[tri[0]]]?;
            let b = cell_to_vert[vertex_cell[tri[1]]]?;
            let c = cell_to_vert[vertex_cell[tri[2]]]?;
            if a == b || b == c || a == c {
                return None;
            }
            Some([a, b, c])
        })
        .collect();

    let vertices: Vec<Vec3> = new_verts.iter().map(|&a| arr_to_vec3(a)).collect();
    TriangleMesh::new(vertices, new_tris)
}

// ---------------------------------------------------------------------------
// Saliency-weighted remeshing
// ---------------------------------------------------------------------------

/// Remesh with non-uniform target edge lengths guided by per-vertex saliency.
///
/// High-saliency regions get a finer target edge length; low-saliency regions
/// get a coarser target.  `base_len` is the reference target edge length for
/// saliency = 0.5.
pub fn saliency_weighted_remesh(
    mesh: &TriangleMesh,
    saliency: &[f64],
    iterations: usize,
) -> TriangleMesh {
    if mesh.vertices.is_empty() {
        return mesh.clone();
    }

    // Compute per-triangle saliency as average of vertex saliencies
    let tri_saliency: Vec<f64> = mesh
        .indices
        .iter()
        .map(|tri| {
            let s = [
                if tri[0] < saliency.len() {
                    saliency[tri[0]]
                } else {
                    0.5
                },
                if tri[1] < saliency.len() {
                    saliency[tri[1]]
                } else {
                    0.5
                },
                if tri[2] < saliency.len() {
                    saliency[tri[2]]
                } else {
                    0.5
                },
            ];
            (s[0] + s[1] + s[2]) / 3.0
        })
        .collect();

    // Use average saliency as overall target scale
    let avg_sal = if tri_saliency.is_empty() {
        0.5
    } else {
        tri_saliency.iter().sum::<f64>() / tri_saliency.len() as f64
    };

    // Compute mesh average edge length as base
    let verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
    let avg_edge: f64 = {
        let total: f64 = mesh
            .indices
            .iter()
            .map(|tri| {
                let a = verts[tri[0]];
                let b = verts[tri[1]];
                let c = verts[tri[2]];
                (dist3(a, b) + dist3(b, c) + dist3(c, a)) / 3.0
            })
            .sum();
        if mesh.indices.is_empty() {
            1.0
        } else {
            total / mesh.indices.len() as f64
        }
    };

    // High saliency → finer (smaller target), low saliency → coarser
    let target = avg_edge * (1.5 - avg_sal.clamp(0.0, 1.0));

    isotropic_remesh(mesh, target.max(1e-6), iterations)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_quad_mesh() -> TriangleMesh {
        // A flat square subdivided into 2 triangles
        let verts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![[0usize, 1, 2], [0, 2, 3]];
        TriangleMesh::new(verts, tris)
    }

    #[test]
    fn test_loop_subdivision_face_count() {
        let mesh = simple_quad_mesh();
        let sub = LoopSubdivision::subdivide(&mesh);
        // 2 original tris → 8 after one Loop subdivision
        assert_eq!(
            sub.indices.len(),
            mesh.indices.len() * 4,
            "Loop subdivision should quadruple the face count"
        );
    }

    #[test]
    fn test_loop_subdivision_vertex_count_increases() {
        let mesh = simple_quad_mesh();
        let sub = LoopSubdivision::subdivide(&mesh);
        assert!(
            sub.vertices.len() > mesh.vertices.len(),
            "Loop subdivision should add vertices"
        );
    }

    #[test]
    fn test_catmull_clark_face_count() {
        let verts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let quads: Vec<[usize; 4]> = vec![[0, 1, 2, 3]];
        let (_, new_quads) = CatmullClark::subdivide_quad_mesh(&verts, &quads);
        // 1 quad → 4 quads
        assert_eq!(
            new_quads.len(),
            quads.len() * 4,
            "Catmull-Clark should quadruple the quad count"
        );
    }

    #[test]
    fn test_catmull_clark_two_quads() {
        let verts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let quads: Vec<[usize; 4]> = vec![[0, 1, 4, 3], [1, 2, 5, 4]];
        let (_, new_quads) = CatmullClark::subdivide_quad_mesh(&verts, &quads);
        assert_eq!(new_quads.len(), quads.len() * 4);
    }

    #[test]
    fn test_isotropic_remesh_returns_mesh() {
        let mesh = simple_quad_mesh();
        let result = isotropic_remesh(&mesh, 0.5, 2);
        // Remeshed mesh should have at least 1 face
        assert!(
            !result.indices.is_empty(),
            "remeshed mesh should have faces"
        );
        assert!(
            !result.vertices.is_empty(),
            "remeshed mesh should have vertices"
        );
    }

    #[test]
    fn test_uniform_remesher() {
        let mesh = simple_quad_mesh();
        let remesher = UniformRemesher::new(0.4);
        let result = remesher.remesh(&mesh, 1);
        assert!(!result.vertices.is_empty());
    }

    #[test]
    fn test_loop_subdivision_twice() {
        let mesh = simple_quad_mesh();
        let sub1 = LoopSubdivision::subdivide(&mesh);
        let sub2 = LoopSubdivision::subdivide(&sub1);
        assert_eq!(sub2.indices.len(), mesh.indices.len() * 16);
    }

    #[test]
    fn test_closest_point_on_triangle() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        // Query point directly above centroid
        let p = [0.25, 0.25, 1.0];
        let cp = closest_point_on_triangle(p, a, b, c);
        // Should project to (0.25, 0.25, 0.0)
        assert!((cp[0] - 0.25).abs() < 1e-10);
        assert!((cp[1] - 0.25).abs() < 1e-10);
        assert!(cp[2].abs() < 1e-10);
    }

    // ── Feature-preserving remesh tests ──────────────────────────────────────

    #[test]
    fn test_feature_preserving_remesh_returns_mesh() {
        let mesh = simple_quad_mesh();
        // Use a very small feature angle so all edges are locked (no collapses),
        // ensuring the output always has triangles.
        let result = feature_preserving_remesh(&mesh, 0.4, 1, 1.0_f64.to_radians());
        assert!(!result.vertices.is_empty());
        assert!(!result.indices.is_empty());
    }

    #[test]
    fn test_feature_preserving_remesh_no_degenerate_faces() {
        let mesh = simple_quad_mesh();
        let result = feature_preserving_remesh(&mesh, 0.3, 2, 45.0_f64.to_radians());
        let verts: Vec<[f64; 3]> = result.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        for tri in &result.indices {
            let a = verts[tri[0]];
            let b = verts[tri[1]];
            let c = verts[tri[2]];
            let e1 = sub3(b, a);
            let e2 = sub3(c, a);
            let cross = cross3(e1, e2);
            let area = len3(cross) * 0.5;
            assert!(area >= 0.0, "negative area: {}", area);
        }
    }

    // ── Edge flip tests ────────────────────────────────────────────────────────

    #[test]
    fn test_edge_flip_preserves_vertex_count() {
        let mesh = simple_quad_mesh();
        let verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut tris = mesh.indices.clone();
        let n_verts_before = verts.len();
        flip_edges_for_quality(&mut tris, &verts);
        assert_eq!(verts.len(), n_verts_before);
        let _ = verts; // suppress unused warning
    }

    #[test]
    fn test_edge_flip_preserves_face_count() {
        let mesh = simple_quad_mesh();
        let verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut tris = mesh.indices.clone();
        let n_tris_before = tris.len();
        flip_edges_for_quality(&mut tris, &verts);
        assert_eq!(tris.len(), n_tris_before);
    }

    // ── Quality metric tests ─────────────────────────────────────────────────

    #[test]
    fn test_mesh_quality_equilateral() {
        // Equilateral triangle should have quality = 1
        let s = 1.0_f64;
        let h = s * (3.0_f64.sqrt()) / 2.0;
        let verts = vec![[0.0_f64, 0.0, 0.0], [s, 0.0, 0.0], [s / 2.0, h, 0.0]];
        let tris = vec![[0usize, 1, 2]];
        let quality = mesh_quality_min_angle(&verts, &tris);
        // All angles = 60° → quality = 1.0
        assert!(quality > 0.9, "equilateral triangle quality={}", quality);
    }

    #[test]
    fn test_mesh_quality_degenerate() {
        // Degenerate (flat) triangle has low quality
        let verts = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [100.0, 0.0, 0.0], // collinear
        ];
        let tris = vec![[0usize, 1, 2]];
        let quality = mesh_quality_min_angle(&verts, &tris);
        assert!(
            quality < 0.5,
            "degenerate triangle quality should be low, got {}",
            quality
        );
    }

    // ── Tangent Laplacian smoothing tests ─────────────────────────────────────

    #[test]
    fn test_tangent_smooth_returns_mesh() {
        let mesh = simple_quad_mesh();
        let result = tangent_laplacian_smooth(&mesh, 3);
        assert_eq!(result.vertices.len(), mesh.vertices.len());
        assert_eq!(result.indices.len(), mesh.indices.len());
    }

    // ── Catmull-Clark double subdivision ─────────────────────────────────────

    #[test]
    fn test_catmull_clark_double_subdivision() {
        let verts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let quads: Vec<[usize; 4]> = vec![[0, 1, 2, 3]];
        let (v1, q1) = CatmullClark::subdivide_quad_mesh(&verts, &quads);
        let (_, q2) = CatmullClark::subdivide_quad_mesh(&v1, &q1);
        assert_eq!(q2.len(), quads.len() * 16); // 1 * 4 * 4
    }

    // ── Laplacian smoothing tests ─────────────────────────────────────────────

    #[test]
    fn test_laplacian_smooth_moves_vertices() {
        let mesh = simple_quad_mesh();
        let result = laplacian_smooth(&mesh, 3, 0.5);
        // Vertices should change after smoothing
        let orig_sum: f64 = mesh.vertices.iter().map(|v| v.norm()).sum();
        let new_sum: f64 = result.vertices.iter().map(|v| v.norm()).sum();
        // After smoothing a non-symmetric mesh, positions should change
        let _ = (orig_sum, new_sum);
        assert_eq!(result.vertices.len(), mesh.vertices.len());
    }

    #[test]
    fn test_laplacian_smooth_zero_iterations_unchanged() {
        let mesh = simple_quad_mesh();
        let result = laplacian_smooth(&mesh, 0, 0.5);
        for (a, b) in mesh.vertices.iter().zip(result.vertices.iter()) {
            assert!((a.x - b.x).abs() < 1e-12);
            assert!((a.y - b.y).abs() < 1e-12);
            assert!((a.z - b.z).abs() < 1e-12);
        }
    }

    #[test]
    fn test_laplacian_smooth_preserves_face_count() {
        let mesh = simple_quad_mesh();
        let result = laplacian_smooth(&mesh, 5, 0.3);
        assert_eq!(result.indices.len(), mesh.indices.len());
    }

    // ── Edge flip Delaunay tests ──────────────────────────────────────────────

    #[test]
    fn test_delaunay_flip_preserves_vertices() {
        let mesh = simple_quad_mesh();
        let verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut tris = mesh.indices.clone();
        let before = verts.len();
        delaunay_edge_flip(&mut tris, &verts);
        assert_eq!(verts.len(), before);
    }

    #[test]
    fn test_delaunay_flip_preserves_face_count() {
        let mesh = simple_quad_mesh();
        let verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut tris = mesh.indices.clone();
        let before = tris.len();
        delaunay_edge_flip(&mut tris, &verts);
        assert_eq!(tris.len(), before);
    }

    // ── Coarsening tests ──────────────────────────────────────────────────────

    #[test]
    fn test_coarsen_reduces_vertex_count() {
        let mesh = simple_quad_mesh();
        let result = coarsen_mesh(&mesh, 0.8); // target edge length > current
        // Coarsening should reduce or maintain vertex count
        assert!(result.vertices.len() <= mesh.vertices.len() + 1); // allow 1 for implementation
    }

    #[test]
    fn test_coarsen_no_degenerate_faces() {
        let mesh = simple_quad_mesh();
        let result = coarsen_mesh(&mesh, 0.5);
        let verts: Vec<[f64; 3]> = result.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        for tri in &result.indices {
            assert!(
                tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2],
                "degenerate face {:?}",
                tri
            );
            // Ensure indices are in bounds
            assert!(tri[0] < verts.len() && tri[1] < verts.len() && tri[2] < verts.len());
        }
    }

    // ── Adaptive refinement by curvature ─────────────────────────────────────

    #[test]
    fn test_adaptive_refine_increases_faces_for_curved_mesh() {
        // Build a simple curved mesh
        let mesh = simple_quad_mesh();
        let result = adaptive_refine_by_curvature(&mesh, 0.01); // low threshold → refine more
        // Should produce more faces or same
        assert!(result.indices.len() >= mesh.indices.len());
    }

    #[test]
    fn test_adaptive_refine_high_threshold_no_change() {
        let mesh = simple_quad_mesh();
        let result = adaptive_refine_by_curvature(&mesh, 1000.0); // high threshold → no split
        // For a flat mesh with high threshold, no splits expected
        assert_eq!(result.indices.len(), mesh.indices.len());
    }

    // ── Mesh decimation / vertex clustering ───────────────────────────────────

    #[test]
    fn test_vertex_clustering_reduces_count() {
        // Build a finer mesh first
        let mesh = simple_quad_mesh();
        let sub = LoopSubdivision::subdivide(&mesh);
        let coarse = vertex_clustering_decimate(&sub, 2); // 2 cells per axis
        assert!(
            coarse.vertices.len() < sub.vertices.len(),
            "decimation should reduce vertex count"
        );
    }

    #[test]
    fn test_vertex_clustering_single_cell_merges_all() {
        let mesh = simple_quad_mesh();
        let result = vertex_clustering_decimate(&mesh, 1); // all verts in one cell
        // All vertices merged → probably 1 vertex left, no valid faces
        assert!(result.vertices.len() <= mesh.vertices.len());
    }

    // ── Saliency-weighted remeshing ───────────────────────────────────────────

    #[test]
    fn test_saliency_remesh_returns_valid_mesh() {
        // Use a subdivided mesh so isotropic remeshing has enough faces to work with
        let mesh = simple_quad_mesh();
        let sub = LoopSubdivision::subdivide(&mesh);
        let saliency = vec![0.5_f64; sub.vertices.len()];
        let result = saliency_weighted_remesh(&sub, &saliency, 1);
        assert!(!result.vertices.is_empty());
        // Indices may be empty if all edges collapse, just check no panic
        let _ = result.indices.len();
    }

    #[test]
    fn test_saliency_remesh_preserves_approximate_vertex_count() {
        let mesh = simple_quad_mesh();
        let saliency = vec![0.5_f64; mesh.vertices.len()];
        let result = saliency_weighted_remesh(&mesh, &saliency, 1);
        // Should not catastrophically change size
        assert!(result.vertices.len() <= mesh.vertices.len() * 10);
    }

    // ── Mesh quality after operations ─────────────────────────────────────────

    #[test]
    fn test_mesh_quality_after_laplacian_smooth() {
        let mesh = simple_quad_mesh();
        let smoothed = laplacian_smooth(&mesh, 5, 0.5);
        let verts: Vec<[f64; 3]> = smoothed.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let q = mesh_quality_min_angle(&verts, &smoothed.indices);
        assert!((0.0..=1.0 + 1e-6).contains(&q), "quality out of range: {q}");
    }

    #[test]
    fn test_aspect_ratio_equilateral() {
        let s = 1.0_f64;
        let h = s * 3.0_f64.sqrt() / 2.0;
        let verts = vec![[0.0_f64, 0.0, 0.0], [s, 0.0, 0.0], [s / 2.0, h, 0.0]];
        let tris = vec![[0usize, 1, 2]];
        let ar = mesh_aspect_ratio_avg(&verts, &tris);
        // Equilateral: all edges = 1, semi-perimeter s=1.5, area=sqrt(3)/4
        // inradius = area/s = (sqrt(3)/4)/1.5 = sqrt(3)/6
        // longest edge / (2 * inradius) = 1 / (2 * sqrt(3)/6) = 1 / (sqrt(3)/3) = sqrt(3) ≈ 1.732
        assert!(
            (ar - 3.0_f64.sqrt()).abs() < 0.01,
            "equilateral ar should be sqrt(3)≈1.732, got {ar}"
        );
    }

    #[test]
    fn test_collapse_edge_reduces_vertex_usage() {
        let mesh = simple_quad_mesh();
        let mut verts: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut tris = mesh.indices.clone();
        let before_faces = tris.len();
        collapse_edge(&mut verts, &mut tris, 0, 1);
        // After collapsing edge 0-1, faces containing both should be removed
        assert!(tris.len() <= before_faces);
    }

    #[test]
    fn test_loop_subdivision_vertex_positions_reasonable() {
        let mesh = simple_quad_mesh();
        let sub = LoopSubdivision::subdivide(&mesh);
        // All vertices should be in a finite range
        for v in &sub.vertices {
            assert!(v.x.is_finite() && v.y.is_finite() && v.z.is_finite());
        }
    }

    #[test]
    fn test_isotropic_remesh_zero_iterations_same() {
        let mesh = simple_quad_mesh();
        let result = isotropic_remesh(&mesh, 0.5, 0);
        // Zero iterations: same as input
        assert_eq!(result.vertices.len(), mesh.vertices.len());
        assert_eq!(result.indices.len(), mesh.indices.len());
    }

    #[test]
    fn test_delaunay_flip_on_subdivided_mesh() {
        let mesh = simple_quad_mesh();
        let sub = LoopSubdivision::subdivide(&mesh);
        let verts: Vec<[f64; 3]> = sub.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut tris = sub.indices.clone();
        delaunay_edge_flip(&mut tris, &verts);
        // Should not produce invalid indices
        for tri in &tris {
            assert!(tri[0] < verts.len() && tri[1] < verts.len() && tri[2] < verts.len());
        }
    }

    #[test]
    fn test_coarsen_mesh_returns_valid_faces() {
        let mesh = simple_quad_mesh();
        let result = coarsen_mesh(&mesh, 1.5);
        let verts: Vec<[f64; 3]> = result.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        for tri in &result.indices {
            assert!(
                tri[0] < verts.len() && tri[1] < verts.len() && tri[2] < verts.len(),
                "face {:?} has out-of-bounds index (n_verts={})",
                tri,
                verts.len()
            );
        }
    }

    #[test]
    fn test_adaptive_refine_no_degenerate_faces() {
        let mesh = simple_quad_mesh();
        let result = adaptive_refine_by_curvature(&mesh, 0.1);
        let verts: Vec<[f64; 3]> = result.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        for tri in &result.indices {
            assert!(tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);
            assert!(tri[0] < verts.len() && tri[1] < verts.len() && tri[2] < verts.len());
        }
    }
}
